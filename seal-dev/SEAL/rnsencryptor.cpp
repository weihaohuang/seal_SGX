#include <algorithm>
#include <stdexcept>
#include "rnsencryptor.h"
#include "util/common.h"
#include "util/uintarith.h"
#include "util/polyarith.h"
#include "util/smallpolyarith.h"
#include "util/polyarithsmallmod.h"
#include "util/polyarithmod.h"
#include "util/polyfftmultmod.h"
#include "util/polyfftmultsmallmod.h"
#include "util/clipnormal.h"
#include "util/randomtostd.h"
#include "util/smallntt.h"
#include "smallmodulus.h"
#include "bigpoly.h"
#include "seal.h"
#include "util/uintarithsmallmod.h"

using namespace std;
using namespace seal::util;

namespace seal
{
    
    RNSEncryptor::RNSEncryptor(const RNSContext &context, const PublicKey &public_key, const MemoryPoolHandle &pool) :
        pool_(pool), parms_(context.get_parms()), qualifiers_(context.get_qualifiers()),
        coeff_mod_array_(context.coeff_mod_array())
    {
        // Verify parameters
        if (!qualifiers_.parameters_set)
        {
            throw invalid_argument("encryption parameters are not valid");
        }
        
        // Verify parameters.
        if (public_key.rns_hash_block_ != parms_.get_hash_block())
        {
            throw invalid_argument("public key is not valid for encryption parameters");
        }
        int coeff_count = parms_.poly_modulus().coeff_count();
        int poly_coeff_uint64_count = parms_.poly_modulus().coeff_uint64_count();
        int coeff_uint64_count = context.coeff_mod_array().size();
        int coeff_mod_count = coeff_mod_array_.size();

        // Set SmallNTTTables
        small_ntt_tables_.resize(coeff_mod_count, pool_);
        small_ntt_tables_ = context.small_ntt_tables_;
        
        // Allocate space and copy over key
        public_key_ = allocate_poly(2 * coeff_count, coeff_uint64_count, pool_);
        set_poly_poly(public_key.get_array().pointer(0), 2 * coeff_count, coeff_uint64_count, public_key_.get());

        // Calculate coeff_modulus / plain_modulus and upper_half_increment.
        coeff_div_plain_modulus_ = allocate_uint(coeff_uint64_count, pool_);
        upper_half_increment_ = allocate_uint(coeff_uint64_count, pool_);
        ConstPointer wide_plain_modulus(duplicate_uint_if_needed(parms_.plain_modulus().pointer(), parms_.plain_modulus().uint64_count(), coeff_uint64_count, false, pool_));
        divide_uint_uint(context.coeff_modulus().pointer(), wide_plain_modulus.get(), coeff_uint64_count, coeff_div_plain_modulus_.get(), upper_half_increment_.get(), pool_);

        // Calculate (plain_modulus + 1) / 2 * coeff_div_plain_modulus.
        Pointer temp(allocate_uint(coeff_uint64_count, pool_));
        upper_half_threshold_ = allocate_uint(coeff_uint64_count, pool_);
        half_round_up_uint(wide_plain_modulus.get(), coeff_uint64_count, temp.get());
        multiply_truncate_uint_uint(temp.get(), coeff_div_plain_modulus_.get(), coeff_uint64_count, upper_half_threshold_.get());

        // Initialize moduli.
        polymod_ = PolyModulus(parms_.poly_modulus().pointer(), coeff_count, poly_coeff_uint64_count);
    }

    void RNSEncryptor::rns_decompose(uint64_t *value)
    {
#ifdef _DEBUG
        if (value == nullptr)
        {
            throw invalid_argument("input cannot be null");
        }
#endif
        int coeff_count = parms_.poly_modulus().coeff_count();
        int coeff_mod_count = coeff_mod_array_.size();
        int total_uint64_count = coeff_mod_count * coeff_count;
        
        Pointer coefficients(allocate_uint(total_uint64_count, pool_));
        //uint64_t *coefficients_ptr = coefficients.get();
        uint64_t mod_ptr;

        // Copy value pointer
        uint64_t *value_ptr = value;

        set_zero_uint(total_uint64_count, coefficients.get());

        for (int i = 0; i < coeff_count; i++)
        {
            for (int j = 0; j < coeff_mod_count; j++)
            {
                mod_ptr = 0;
                small_modulo_uint(value_ptr, coeff_mod_count, coeff_mod_array_[j], &mod_ptr, pool_);
                //coefficients_ptr[(j * coeff_count) + i] = mod_ptr;
                *(coefficients.get() + i + (j * coeff_count)) = mod_ptr;
            }
            value_ptr += coeff_mod_count;
        }
        set_uint_uint(coefficients.get(), total_uint64_count, value);
    }

    void RNSEncryptor::rns_encrypt(const Plaintext &plain, Ciphertext &destination)
    {
        const BigPoly &plain_poly = plain.get_poly();
        BigPolyArray &destination_array = destination.get_mutable_array();

        int coeff_count = parms_.poly_modulus().coeff_count();
        int poly_coeff_uint64_count = parms_.poly_modulus().coeff_uint64_count();
        int coeff_mod_count = coeff_mod_array_.size();
#ifdef _DEBUG
        if (plain_poly.significant_coeff_count() >= coeff_count || !are_poly_coefficients_less_than(plain_poly.pointer(), plain_poly.coeff_count(), plain_poly.coeff_uint64_count(), parms_.plain_modulus().pointer(), parms_.plain_modulus().uint64_count()))
        {
            throw invalid_argument("plain is too large to be represented by encryption parameters");
        }
#endif
        // Make destination have right size
        destination_array.resize(2, coeff_count, coeff_mod_count * bits_per_uint64);

        /*
        Ciphertext (c_0,c_1) should be a BigPolyArray
        c_0 = Delta * m + public_key_[0] * u + e_1 where u sampled from R_2 and e_1 sampled from chi.
        c_1 = public_key_[1] * u + e_2 where e_2 sampled from chi.
        */

        // Multiply plain by scalar coeff_div_plaintext and reposition if in upper-half.
        // Result gets added into the c_0 term of ciphertext (c_0,c_1).
        rns_preencrypt(plain_poly.pointer(), plain_poly.coeff_count(), plain_poly.coeff_uint64_count(), destination_array.pointer(0));

        // Generate u 
        Pointer u(allocate_poly(coeff_count, coeff_mod_count, pool_));
        unique_ptr<UniformRandomGenerator> random(parms_.random_generator()->create());
        
        set_poly_coeffs_zero_one_negone_rns(u.get(), random.get());
        //set_poly_coeffs_zero_one(u.get(), random.get());
        
        // Calculate public_key_[0] * u.
        // Since we (may) need both public_key_[0] and u later, need a temp variable to store the solution.
        // Add temp into destination_array[0].
        Pointer temp(allocate_poly(coeff_count, coeff_mod_count, pool_));

        // Multiply both u*public_key_[0] and u*public_key_[1] using the same FFT
        set_zero_uint(coeff_mod_count, get_poly_coeff(temp.get(), coeff_count - 1, 1));
        set_zero_uint(coeff_mod_count, get_poly_coeff(destination_array.pointer(1), coeff_count - 1, 1));
        
        for (int i = 0; i < coeff_mod_count; i++)
        {
            smallntt_double_multiply_poly_nttpoly(u.get() + (i * coeff_count), public_key_.get() + (i * coeff_count), public_key_.get() + (coeff_count * coeff_mod_count) + (i * coeff_count), small_ntt_tables_[i], temp.get() + (i * coeff_count), destination_array.pointer(1) + (i * coeff_count), pool_);
            add_poly_poly_coeff_smallmod(temp.get() + (i * coeff_count), destination_array.pointer(0) + (i * coeff_count), coeff_count, coeff_mod_array_[i], destination_array.pointer(0) + (i * coeff_count));
        }
        // Generate e_0, add this value into destination_array[0].
        set_poly_coeffs_normal_rns(temp.get(), random.get());
        for (int i = 0; i < coeff_mod_count; i++)
        {
            add_poly_poly_coeff_smallmod(temp.get() + (i * coeff_count), destination_array.pointer(0) + (i * coeff_count), coeff_count, coeff_mod_array_[i], destination_array.pointer(0) + (i * coeff_count));
        }
        // Generate e_1, add this value into destination_array[1].
        set_poly_coeffs_normal_rns(temp.get(), random.get());
        for (int i = 0; i < coeff_mod_count; i++)
        {
            add_poly_poly_coeff_smallmod(temp.get() + (i * coeff_count), destination_array.pointer(1) + (i * coeff_count), coeff_count, coeff_mod_array_[i], destination_array.pointer(1) + (i * coeff_count));
        }

        // Set the hash block
        destination.rns_hash_block_ = parms_.get_hash_block();
    }

    void RNSEncryptor::rns_preencrypt(const uint64_t *plain, int plain_coeff_count, int plain_coeff_uint64_count, uint64_t *destination)
    {
        int coeff_count = parms_.poly_modulus().coeff_count();
        int coeff_uint64_count = coeff_mod_array_.size();

        if (plain_coeff_count > coeff_count)
        {
            plain_coeff_count = coeff_count;
        }
        uint64_t *destination_cpy = destination;
        // Multiply plain by scalar coeff_div_plain_modulus_ and reposition if in upper-half.
        if (plain == destination)
        {
            // If plain and destination are same poly, then need another storage for multiply output.
            Pointer temp(allocate_uint(coeff_uint64_count, pool_));
            for (int i = 0; i < plain_coeff_count; i++)
            {
                multiply_uint_uint(plain, plain_coeff_uint64_count, coeff_div_plain_modulus_.get(), coeff_uint64_count, coeff_uint64_count, temp.get());
                bool is_upper_half = is_greater_than_or_equal_uint_uint(temp.get(), upper_half_threshold_.get(), coeff_uint64_count);
                if (is_upper_half)
                {
                    add_uint_uint(temp.get(), upper_half_increment_.get(), coeff_uint64_count, destination);
                }
                else
                {
                    set_uint_uint(temp.get(), coeff_uint64_count, destination);
                }
                plain += plain_coeff_uint64_count;
                destination += coeff_uint64_count;
            }
        }
        else
        {
            for (int i = 0; i < plain_coeff_count; i++)
            {
                //Multiply plain by coeff_div_plain_modulus_ and put the result in destination
                multiply_uint_uint(plain, plain_coeff_uint64_count, coeff_div_plain_modulus_.get(), coeff_uint64_count, coeff_uint64_count, destination);

                // check if destination >= upper half threshold
                bool is_upper_half = is_greater_than_or_equal_uint_uint(destination, upper_half_threshold_.get(), coeff_uint64_count);
                if (is_upper_half)
                {
                    add_uint_uint(destination, upper_half_increment_.get(), coeff_uint64_count, destination);
                }
                plain += plain_coeff_uint64_count;
                destination += coeff_uint64_count;
            }
        }
        // Zero any remaining coefficients.
        for (int i = plain_coeff_count; i < coeff_count; i++)
        {
            set_zero_uint(coeff_uint64_count, destination);
            destination += coeff_uint64_count;
        }

        // Compute the RNS decomposition form of the destination result
        rns_decompose(destination_cpy);
    }

    void RNSEncryptor::set_poly_coeffs_zero_one_negone_rns(uint64_t *poly, UniformRandomGenerator *random) const
    {
        int coeff_count = parms_.poly_modulus().coeff_count();
        int coeff_mod_count = coeff_mod_array_.size();

        Pointer coeff_modulus_minus_one(allocate_uint(coeff_mod_count, pool_));
        for (int i = 0; i < coeff_mod_count; i++)
        {
            decrement_uint(coeff_mod_array_[i].pointer(), coeff_mod_array_[i].uint64_count(), coeff_modulus_minus_one.get() + i);
        }

        RandomToStandardAdapter engine(random);
        uniform_int_distribution<int> dist(-1, 1);
        for (int i = 0; i < coeff_count - 1; i++)
        {
            int rand_index = dist(engine);
            if (rand_index == 1)
            {
                for (int j = 0; j < coeff_mod_count; j++)
                {
                    poly[i + (j * coeff_count)] = 1;
                }
            }
            else if (rand_index == -1)
            {
                for (int j = 0; j < coeff_mod_count; j++)
                {
                    poly[i + (j * coeff_count)] = *(coeff_modulus_minus_one.get() + j);
                }
            }
            else
            {
                for (int j = 0; j < coeff_mod_count; j++)
                {
                    poly[i + (j * coeff_count)] = 0;
                }
            }
        }

        // Set the last coefficient equal to zero in RNS representation
        for (int i = 0; i < coeff_mod_count; i++)
        {
            poly[(coeff_count - 1) + (i * coeff_count)] = 0;
        }
    }

    void RNSEncryptor::set_poly_coeffs_zero_one_rns(uint64_t *poly, UniformRandomGenerator *random) const
    {
        int coeff_count = parms_.poly_modulus().coeff_count();
        int coeff_mod_count = coeff_mod_array_.size();

        RandomToStandardAdapter engine(random);
        uniform_int_distribution<int> dist(0, 1);

        set_zero_poly(coeff_count, coeff_mod_count, poly);
        for (int i = 0; i < coeff_count; i++)
        {
            int rand_index = dist(engine);
            for (int j = 0; j < coeff_mod_count; j++)
            {
                poly[i + (j * coeff_count)] = rand_index;
            }
        }
    }

    void RNSEncryptor::set_poly_coeffs_normal_rns(uint64_t *poly, UniformRandomGenerator *random) const
    {
        int coeff_count = parms_.poly_modulus().coeff_count();
        int coeff_mod_count = coeff_mod_array_.size();
        if (parms_.noise_standard_deviation() == 0 || parms_.noise_max_deviation() == 0)
        {
            set_zero_poly(coeff_count, coeff_mod_count, poly);
            return;
        }
        RandomToStandardAdapter engine(random);
        ClippedNormalDistribution dist(0, parms_.noise_standard_deviation(), parms_.noise_max_deviation());
        for (int i = 0; i < coeff_count - 1; i++)
        {
            int64_t noise = static_cast<int64_t>(dist(engine));
            if (noise > 0)
            {
                for (int j = 0; j < coeff_mod_count; j++)
                {
                    poly[i + (j * coeff_count)] = static_cast<uint64_t>(noise);
                }
            }
            else if (noise < 0)
            {
                noise = -noise;
                for (int j = 0; j < coeff_mod_count; j++)
                {
                    poly[i + (j * coeff_count)] = coeff_mod_array_[j].value() - static_cast<uint64_t>(noise);
                }
            }
            else
            {
                for (int j = 0; j < coeff_mod_count; j++)
                {
                    poly[i + (j * coeff_count)] = 0;
                }
            }
        }
        
        // Set the last coefficient equal to zero in RNS representation
        for (int i = 0; i < coeff_mod_count; i++)
        {
            poly[(coeff_count - 1) + (i * coeff_count)] = 0;
        }
    }

    RNSEncryptor::RNSEncryptor(const RNSEncryptor &copy) :
        pool_(copy.pool_), parms_(copy.parms_), qualifiers_(copy.qualifiers_),
        small_ntt_tables_(copy.small_ntt_tables_),
        coeff_mod_array_(copy.coeff_mod_array_)
    {
        int coeff_count = parms_.poly_modulus().coeff_count();
        int poly_coeff_uint64_count = parms_.poly_modulus().coeff_uint64_count();
        int coeff_uint64_count = coeff_mod_array_.size();

        // Allocate and copy over data
        upper_half_threshold_ = allocate_uint(coeff_uint64_count, pool_);
        set_uint_uint(copy.upper_half_threshold_.get(), coeff_uint64_count, upper_half_threshold_.get());

        upper_half_increment_ = allocate_uint(coeff_uint64_count, pool_);
        set_uint_uint(copy.upper_half_increment_.get(), coeff_uint64_count, upper_half_increment_.get());

        coeff_div_plain_modulus_ = allocate_uint(coeff_uint64_count, pool_);
        set_uint_uint(copy.coeff_div_plain_modulus_.get(), coeff_uint64_count, coeff_div_plain_modulus_.get());

        public_key_ = allocate_poly(2 * coeff_count, coeff_uint64_count, pool_);
        set_poly_poly(copy.public_key_.get(), 2 * coeff_count, coeff_uint64_count, public_key_.get());

        // Initialize moduli.
        polymod_ = PolyModulus(parms_.poly_modulus().pointer(), coeff_count, poly_coeff_uint64_count);
    }
}