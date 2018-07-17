#include <algorithm>
#include <stdexcept>
#include "rnsdecryptor.h"
#include "util/common.h"
#include "util/uintcore.h"
#include "util/uintarith.h"
#include "util/polycore.h"
#include "util/smallpolyarith.h"
#include "util/polyarithsmallmod.h"
#include "bigpoly.h"
#include "util/uintarithmod.h"
#include "util/polyextras.h"
#include "util/polyfftmultsmallmod.h"

using namespace std;
using namespace seal::util;

namespace seal
{
    RNSDecryptor::RNSDecryptor(const RNSContext &context, const SecretKey &secret_key, const MemoryPoolHandle &pool) :
        pool_(pool), parms_(context.get_parms()), qualifiers_(context.get_qualifiers()),
        base_convertor_(context.base_convertor()), coeff_mod_array_(context.coeff_mod_array())
    {
        // Verify parameters
        if (!qualifiers_.parameters_set)
        {
            throw invalid_argument("encryption parameters are not set correctly");
        }

        // Verify parameters.
        if (secret_key.rns_hash_block_ != parms_.get_hash_block())
        {
            throw invalid_argument("secret key is not valid for encryption parameters");
        }
        
        int coeff_count = parms_.poly_modulus().coeff_count();
        int poly_coeff_uint64_count = parms_.poly_modulus().coeff_uint64_count();
        int coeff_mod_count = base_convertor_.coeff_base_mod_count();

        // Set SmallNTTTables
        small_ntt_tables_.resize(coeff_mod_count, pool_);
        small_ntt_tables_ = context.small_ntt_tables_;

        // Populate coeff products array for compose functions (used in noise budget)
        coeff_products_array_ = allocate_uint(coeff_mod_count * coeff_mod_count, pool_);
        Pointer tmp_coeff(allocate_uint(coeff_mod_count, pool_));
        set_zero_uint(coeff_mod_count * coeff_mod_count, coeff_products_array_.get());

        for (int i = 0; i < coeff_mod_count; i++)
        {
            *(coeff_products_array_.get() + (i * coeff_mod_count)) = 1;
            for (int j = 0; j < coeff_mod_count; j++)
            {
                if (i != j)
                {
                    multiply_uint_uint64(coeff_products_array_.get() + (i * coeff_mod_count), coeff_mod_count, coeff_mod_array_[j].value(), coeff_mod_count, tmp_coeff.get());
                    set_uint_uint(tmp_coeff.get(), coeff_mod_count, coeff_products_array_.get() + (i * coeff_mod_count));
                }
            }
        }

        // Set the negative inverse of coeff moduli products mod plain modulus U gamma array
        neg_inv_coeff_mod_ = base_convertor_.neg_inv_coeff();

        // Set the plain modulus gamma product mod coeff moduli array
        plain_gamma_product_ = base_convertor_.plain_gamma_product();

        // Set plain_modulus, gamma array 
        plain_gamma_array_ = base_convertor_.get_plain_gamma_array();

        // Set gamma inverse mod plain_modulus 
        inv_gamma_mod_plain_modulus_ = base_convertor_.get_inv_gamma();

        // Set inverse coeff mod coeff array
        inv_coeff_mod_coeff_array_ = base_convertor_.get_inv_coeff_mod_coeff_array();

        // Allocate secret_key_ and copy over value
        secret_key_ = allocate_poly(coeff_count, coeff_mod_count, pool_);
        set_poly_poly(secret_key.get_poly().pointer(), coeff_count, coeff_mod_count, secret_key_.get());

        // Set the secret_key_array to have size 1 (first power of secret) 
        secret_key_array_.resize(1, coeff_count, coeff_mod_count * bits_per_uint64);
        set_poly_poly(secret_key_.get(), coeff_count, coeff_mod_count, secret_key_array_.pointer(0));

        // Set the big coeff modulus for noise computation
        product_modulus_ = allocate_uint(coeff_mod_count, pool_);
        set_uint_uint(context.coeff_modulus_.pointer(), coeff_mod_count, product_modulus_.get());

        // Initialize moduli.
        mod_ = Modulus(product_modulus_.get(), coeff_mod_count);
        polymod_ = PolyModulus(parms_.poly_modulus().pointer(), coeff_count, poly_coeff_uint64_count);
    }

    RNSDecryptor::RNSDecryptor(const RNSDecryptor &copy) :
        pool_(copy.pool_), parms_(copy.parms_), qualifiers_(copy.qualifiers_),
        base_convertor_(copy.base_convertor_), 
        small_ntt_tables_(copy.small_ntt_tables_),
        inv_gamma_mod_plain_modulus_(copy.inv_gamma_mod_plain_modulus_), 
        coeff_mod_array_(copy.coeff_mod_array_),
        plain_gamma_array_(copy.plain_gamma_array_),
        plain_gamma_product_(copy.plain_gamma_product_),
        inv_coeff_mod_coeff_array_(copy.inv_coeff_mod_coeff_array_),
        neg_inv_coeff_mod_(copy.neg_inv_coeff_mod_),
        secret_key_array_(copy.secret_key_array_)
    {
        int coeff_count = parms_.poly_modulus().coeff_count();
        int poly_coeff_uint64_count = parms_.poly_modulus().coeff_uint64_count();
        int coeff_mod_count = base_convertor_.coeff_base_mod_count();

        // Populate coeff products array for compose functions (used in noise budget)
        coeff_products_array_ = allocate_uint(coeff_mod_count * coeff_mod_count, pool_);
        set_uint_uint(copy.coeff_products_array_.get(), coeff_mod_count * coeff_mod_count, coeff_products_array_.get());

        // Allocate secret_key_ and copy over value
        secret_key_ = allocate_poly(coeff_count, coeff_mod_count, pool_);
        set_poly_poly(copy.secret_key_.get(), coeff_count, coeff_mod_count, secret_key_.get());

        // Set the big coeff modulus for noise computation
        product_modulus_ = allocate_uint(coeff_mod_count, pool_);
        set_uint_uint(copy.product_modulus_.get(), coeff_mod_count, product_modulus_.get());

        // Initialize moduli.
        mod_ = Modulus(product_modulus_.get(), coeff_mod_count);
        polymod_ = PolyModulus(parms_.poly_modulus().pointer(), coeff_count, poly_coeff_uint64_count);
    }

    void RNSDecryptor::rns_decrypt(const Ciphertext &encrypted, Plaintext &destination)
    {
        int coeff_count = parms_.poly_modulus().coeff_count();
        int coeff_mod_count = base_convertor_.coeff_base_mod_count();
        int coeff_bit_count = coeff_mod_count * bits_per_uint64;
        int array_poly_uint64_count = coeff_count * coeff_mod_count;
        
        // The number of uint64 count for plain_modulus and gamma together
        int plain_gamma_uint64_count = 2;

        const BigPolyArray &encrypted_array = encrypted.get_array();
        BigPoly &destination_poly = destination.get_poly();

        // Verify parameters.
        if (encrypted.rns_hash_block_ != parms_.get_hash_block())
        {
            throw invalid_argument("encrypted is not valid for encryption parameters");
        }

        // Make sure destination is of right size to perform all computations. At the end 
        // we will resize the coefficients to be the size of plain_modulus.
        if (destination_poly.coeff_count() != coeff_count || destination_poly.coeff_bit_count() != coeff_bit_count)
        {
            destination_poly.resize(coeff_count, coeff_bit_count);
        }

        // Make sure we have enough secret keys computed
        compute_secret_key_array(encrypted_array.size() - 1);

        /*
        Firstly find c_0 + c_1 *s + ... + c_{count-1} * s^{count-1} mod q
        This is equal to Delta m + v where ||v|| < Delta/2.
        So, add Delta / 2 and now we have something which is Delta * (m + epsilon) where epsilon < 1
        Therefore, we can (integer) divide by Delta and the answer will round down to m.
        */

        // Make a temp destination for all the arithmetic mod qi before calling FastBConverse
        Pointer tmp_dest_modq(allocate_poly(coeff_count, coeff_mod_count, pool_));
        set_zero_poly(coeff_count, coeff_mod_count, tmp_dest_modq.get());

        // put < (c_1 , c_2, ... , c_{count-1}) , (s,s^2,...,s^{count-1}) > mod q in destination
        if (qualifiers_.enable_ntt)
        {
            // Make a copy of the encrypted BigPolyArray for NTT (except the first polynomial is not needed)
            Pointer encrypted_copy(allocate_poly((encrypted_array.size() - 1) * coeff_count, coeff_mod_count, pool_));
            set_poly_poly(encrypted_array.pointer(1), (encrypted_array.size() - 1) * coeff_count, coeff_mod_count, encrypted_copy.get());


            // Now do the dot product of encrypted_copy and the secret key array using NTT. The secret key powers are already NTT transformed.
            for (int i = 0; i < coeff_mod_count; i++)
            {
                smallntt_dot_product_bigpolyarray_nttbigpolyarray(encrypted_copy.get() + (i * coeff_count), secret_key_array_.pointer(0) + (i * coeff_count), encrypted_array.size() - 1, array_poly_uint64_count, small_ntt_tables_[i], tmp_dest_modq.get() + (i * coeff_count), pool_);
            }
        }
        else
        {
            // This branch should never be reached
            throw logic_error("invalid encryption parameters");
        }

        // add c_0 into destination
        for (int i = 0; i < coeff_mod_count; i++)
        {
            add_poly_poly_coeff_smallmod(tmp_dest_modq.get() + (i * coeff_count), encrypted_array.pointer(0) + (i * coeff_count), coeff_count, coeff_mod_array_[i], tmp_dest_modq.get() + (i * coeff_count));

            // Compute |gamma * plain|qi * ct(s)
            multiply_poly_scalar_coeff_smallmod(tmp_dest_modq.get() + (i * coeff_count), coeff_count, &plain_gamma_product_[i], coeff_mod_array_[i], tmp_dest_modq.get() + (i * coeff_count));
        }
        
        // Make another temp destination to get the poly in mod {gamma U plain_modulus}
        Pointer tmp_dest_plain_gamma(allocate_poly(coeff_count, plain_gamma_uint64_count, pool_));

        // Compute FastBConvert from q to {gamma, plain_modulus}
        base_convertor_.fastbconv_plain_gamma(tmp_dest_modq.get(), tmp_dest_plain_gamma.get());
        
        // Compute result multiply by coeff_modulus inverse in mod {gamma U plain_modulus}
        for (int i = 0; i < plain_gamma_uint64_count; i++)
        {
            multiply_poly_scalar_coeff_smallmod(tmp_dest_plain_gamma.get() + (i * coeff_count), coeff_count, &neg_inv_coeff_mod_[i], plain_gamma_array_[i], tmp_dest_plain_gamma.get() + (i * coeff_count));
        }

        // Resize the coefficient to the original plain_modulus size
        destination_poly.resize(coeff_count, parms_.plain_modulus().bit_count());

        // First correct the values which are larger than floor(gamma/2)
        uint64_t gamma_div_2 = plain_gamma_array_[1].value() >> 1;

        // Now compute the subtraction to remove error and perform final multiplication by gamma inverse mod plain_modulus
        for (int i = 0; i < coeff_count; i++)
        {
            // Need correction beacuse of center mod
            if (*(tmp_dest_plain_gamma.get() + (i + coeff_count)) > gamma_div_2)
            {
                // Compute -(gamma - a) instead of (a - gamma)
                *(tmp_dest_plain_gamma.get() + (i + coeff_count)) = plain_gamma_array_[1].value() - *(tmp_dest_plain_gamma.get() + (i + coeff_count));
                *(tmp_dest_plain_gamma.get() + (i + coeff_count)) %= plain_gamma_array_[0].value();
                add_uint_uint_smallmod(tmp_dest_plain_gamma.get() + i, tmp_dest_plain_gamma.get() + (i + coeff_count), plain_gamma_array_[0], destination_poly.pointer() + i);
            }
            // No correction needed
            else
            {
                *(tmp_dest_plain_gamma.get() + (i + coeff_count)) %= plain_gamma_array_[0].value();
                sub_uint_uint_smallmod(tmp_dest_plain_gamma.get() + i, tmp_dest_plain_gamma.get() + (i + coeff_count), plain_gamma_array_[0], destination_poly.pointer() + i);
            }
        }

        // Perform final multiplication by gamma inverse mod plain_modulus
        multiply_poly_scalar_coeff_smallmod(destination_poly.pointer(), coeff_count, &inv_gamma_mod_plain_modulus_, plain_gamma_array_[0], destination_poly.pointer());
    }

    void RNSDecryptor::compute_secret_key_array(int max_power)
    {
        int old_count = secret_key_array_.size();
        int new_count = max(max_power, secret_key_array_.size());

        if (old_count == new_count)
        {
            return;
        }

        int coeff_count = parms_.poly_modulus().coeff_count();
        int coeff_mod_count = base_convertor_.coeff_base_mod_count();
        int coeff_bit_count = coeff_mod_count * bits_per_uint64;


        // Compute powers of secret key until max_power
        secret_key_array_.resize(new_count, coeff_count, coeff_bit_count);

        int poly_ptr_increment = coeff_count * coeff_mod_count;
        uint64_t *prev_poly_ptr = secret_key_array_.pointer(old_count - 1);
        uint64_t *next_poly_ptr = prev_poly_ptr + poly_ptr_increment;

        if (qualifiers_.enable_ntt)
        {
            // Since all of the key powers in secret_key_array_ are already NTT transformed, to get the next one 
            // we simply need to compute a dyadic product of the last one with the first one [which is equal to NTT(secret_key_)].
            for (int i = old_count; i < new_count; i++)
            {
                for (int j = 0; j < coeff_mod_count; j++)
                {
                    dyadic_product_coeff_smallmod(prev_poly_ptr + (j * coeff_count), secret_key_array_.pointer(0) + (j * coeff_count), coeff_count, coeff_mod_array_[j], next_poly_ptr + (j * coeff_count), pool_);
                }
                prev_poly_ptr = next_poly_ptr;
                next_poly_ptr += poly_ptr_increment;
            }
        }
        else
        {
            // This branch should never be reached
            throw logic_error("invalid encryption parameters");
        }
    }

    void RNSDecryptor::rns_compose(uint64_t *value)
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

        // Set temporary coefficients_ptr pointer to point to either an existing allocation given as parameter,
        // or else to a new allocation from the memory pool.
        Pointer coefficients(allocate_uint(total_uint64_count, pool_));
        uint64_t *coefficients_ptr = coefficients.get();

        set_zero_uint(total_uint64_count, coefficients_ptr);

        // Re-merge the coefficients first
        for (int i = 0; i < coeff_count; i++)
        {
            for (int j = 0; j < coeff_mod_count; j++)
            {
                coefficients_ptr[(i * coeff_mod_count) + j] = value[(j * coeff_count) + i];
            }
        }

        Pointer big_alloc(allocate_zero_uint(coeff_mod_count, pool_));
        set_zero_uint(total_uint64_count, value);

        uint64_t* value_ptr = value;
        for (int i = 0; i < coeff_count; i++)
        {
            for (int j = 0; j < coeff_mod_count; j++)
            {
                uint64_t tmp;
                multiply_uint64_smallmod(coefficients_ptr + j, &inv_coeff_mod_coeff_array_[j], coeff_mod_array_[j], &tmp);
                multiply_uint_uint64(coeff_products_array_.get() + (j * coeff_mod_count), coeff_mod_count, tmp, coeff_mod_count, big_alloc.get());
                add_uint_uint_mod(big_alloc.get(), value + (i * coeff_mod_count), mod_.get(), coeff_mod_count, value + (i * coeff_mod_count));
            }
            set_zero_uint(coeff_mod_count, big_alloc.get());
            coefficients_ptr += coeff_mod_count;
        }
    }

    int RNSDecryptor::invariant_noise_budget(const Ciphertext &encrypted)
    {
        int coeff_count = parms_.poly_modulus().coeff_count();
        int coeff_mod_count = coeff_mod_array_.size();
        int array_poly_uint64_count = coeff_count * coeff_mod_count;

        const BigPolyArray &encrypted_array = encrypted.get_array();

        // Verify parameters.
        if (encrypted.rns_hash_block_ != parms_.get_hash_block())
        {
            throw invalid_argument("encrypted is not valid for encryption parameters");
        }

        // Make sure destination is of right size.
        Pointer destination(allocate_uint(coeff_mod_count, pool_));

        // Now need to compute c(s) - Delta*m (mod q)

        // Make sure we have enough secret keys computed
        compute_secret_key_array(encrypted_array.size() - 1);

        Pointer noise_poly(allocate_poly(coeff_count, coeff_mod_count, pool_));
        Pointer plain_poly(allocate_poly(coeff_count, coeff_mod_count, pool_));
        Pointer plain_poly_copy(allocate_poly(coeff_count, coeff_mod_count, pool_));
        Pointer result_poly(allocate_poly(coeff_count, coeff_mod_count, pool_));

        /*
        Firstly find c_0 + c_1 *s + ... + c_{count-1} * s^{count-1} mod q
        This is equal to Delta m + v where ||v|| < Delta/2.
        */
        // put < (c_1 , c_2, ... , c_{count-1}) , (s,s^2,...,s^{count-1}) > mod q in destination_poly
        if (qualifiers_.enable_ntt)
        {
            // Make a copy of the encrypted BigPolyArray for NTT (except the first polynomial is not needed)
            Pointer encrypted_copy(allocate_poly((encrypted_array.size() - 1) * coeff_count, coeff_mod_count, pool_));
            set_poly_poly(encrypted_array.pointer(1), (encrypted_array.size() - 1) * coeff_count, coeff_mod_count, encrypted_copy.get());

            // Now do the dot product of encrypted_copy and the secret key array using NTT. The secret key powers are already NTT transformed.
            for (int i = 0; i < coeff_mod_count; i++)
            {
                smallntt_dot_product_bigpolyarray_nttbigpolyarray(encrypted_copy.get() + (i * coeff_count), secret_key_array_.pointer(0) + (i * coeff_count), encrypted_array.size() - 1, array_poly_uint64_count, small_ntt_tables_[i], noise_poly.get() + (i * coeff_count), pool_);
                
                // add c_0 into noise_poly
                add_poly_poly_coeff_smallmod(noise_poly.get() + (i * coeff_count), encrypted_array.pointer(0) + (i * coeff_count), coeff_count, coeff_mod_array_[i], noise_poly.get() + (i * coeff_count));

                // Multiply by parms_.plain_modulus() and reduce mod parms_.coeff_modulus() to get parms_.coeff_modulus()*noise
                ConstPointer wide_plain_modulus(duplicate_uint_if_needed(parms_.plain_modulus().pointer(), parms_.plain_modulus().uint64_count(), coeff_mod_array_[i].uint64_count(), false, pool_));
                multiply_poly_scalar_coeff_smallmod(noise_poly.get() + (i * coeff_count), coeff_count, wide_plain_modulus.get(), coeff_mod_array_[i], noise_poly.get() + (i * coeff_count));
            }
        }
        else
        {
            // This branch should never be reached
            throw logic_error("invalid encryption parameters");
        }

        // Compose the noise
        rns_compose(noise_poly.get());
        
        // Next we compute the infinity norm mod parms_.coeff_modulus()
        poly_infty_norm_coeffmod(noise_poly.get(), coeff_count, coeff_mod_count, mod_, destination.get(), pool_);

        // The -1 accounts for scaling the invariant noise by 2 
		return max(0, mod_.significant_bit_count() - get_significant_bit_count_uint(destination.get(), coeff_mod_count) - 1);
    }
}
