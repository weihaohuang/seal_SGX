#include <algorithm>
#include <stdexcept>
#include <chrono>
#include "rnsevaluator.h"
#include "util/common.h"
#include "util/uintcore.h"
#include "util/uintarith.h"
#include "util/uintarithmod.h"
#include "util/uintarithsmallmod.h"
#include "util/polycore.h"
#include "util/smallpolyarith.h"
#include "util/polyarithsmallmod.h"
#include "util/polyfftmult.h"
#include "util/polyfftmultsmallmod.h"
#include "bigpoly.h"
#include "util/smallntt.h"
#include "util/modulus.h"

using namespace std;
using namespace seal::util;

namespace seal
{
    namespace
    {
        ConstPointer duplicate_bigpolyarray_if_needed(const BigPolyArray &operand, bool force, MemoryPool &pool_)
        {
            return duplicate_if_needed(operand.pointer(0), operand.size() * operand.coeff_count() * operand.coeff_uint64_count(), force, pool_);
        }
    }

    RNSEvaluator::RNSEvaluator(const RNSContext &context, const RNSEvaluationKeys &evaluation_keys, const MemoryPoolHandle &pool) :
        pool_(pool), parms_(context.get_parms()), qualifiers_(context.get_qualifiers()), 
        base_convertor_(context.base_convertor()), 
        evaluation_keys_(evaluation_keys), 
        coeff_mod_array_(context.coeff_mod_array()) 
    {
        // Verify parameters
        if (!qualifiers_.parameters_set)
        {
            throw invalid_argument("encryption parameters are not set correctly");
        }
        if (evaluation_keys.size() > 0)
        {
            if (evaluation_keys.get_hash_block() != parms_.get_hash_block())
            {
                throw invalid_argument("evaluation keys are not valid for encryption parameters");
            }
        }

        int coeff_count = parms_.poly_modulus().coeff_count();
        int poly_coeff_uint64_count = parms_.poly_modulus().coeff_uint64_count();
        int coeff_mod_count = coeff_mod_array_.size();
        int coeff_bit_count = coeff_mod_count * bits_per_uint64;
        bsk_base_mod_count_ = base_convertor_.bsk_base_mod_count();
        
        // Set SmallNTTTables
        bsk_small_ntt_tables_.resize(bsk_base_mod_count_, pool_);
        bsk_small_ntt_tables_ = base_convertor_.get_bsk_small_ntt_table();

        coeff_small_ntt_tables_.resize(coeff_mod_count, pool_);
        coeff_small_ntt_tables_ = context.small_ntt_tables_;

        // Copy over bsk moduli array
        bsk_mod_array_ = base_convertor_.get_bsk_mod_array();

        // Copy over inverse of coeff moduli products mod each coeff moduli
        inv_coeff_products_mod_coeff_array_ = base_convertor_.get_inv_coeff_mod_coeff_array();

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

        // Calculate coeff_modulus / plain_modulus.
        coeff_div_plain_modulus_ = allocate_uint(coeff_mod_count, pool_);
        ConstPointer wide_plain_modulus(duplicate_uint_if_needed(parms_.plain_modulus().pointer(), parms_.plain_modulus().uint64_count(), coeff_mod_count, false, pool_));
        Pointer temp(allocate_uint(coeff_mod_count, pool_));
        divide_uint_uint(context.coeff_modulus().pointer(), wide_plain_modulus.get(), coeff_mod_count, coeff_div_plain_modulus_.get(), temp.get(), pool_);

        // Calculate (plain_modulus + 1) / 2.
        plain_upper_half_threshold_ = allocate_uint(coeff_mod_count, pool_);
        half_round_up_uint(wide_plain_modulus.get(), coeff_mod_count, plain_upper_half_threshold_.get());

        // Calculate coeff_modulus - plain_modulus.
        plain_upper_half_increment_ = allocate_uint(coeff_mod_count, pool_);
        sub_uint_uint(context.coeff_modulus().pointer(), wide_plain_modulus.get(), coeff_mod_count, plain_upper_half_increment_.get());

		// Calculate coeff_modulus[i] - plain_modulus if not enable_plain_lift
		if (qualifiers_.enable_fast_plain_lift)
		{
			plain_upper_half_increment_array_.resize(coeff_mod_count);
			for (int i = 0; i < coeff_mod_count; i++)
			{
				plain_upper_half_increment_array_[i] = coeff_mod_array_[i].value() - parms_.plain_modulus().value();
			}
		}

        // Calculate (plain_modulus + 1) / 2 * coeff_div_plain_modulus.
        upper_half_threshold_ = allocate_uint(coeff_mod_count, pool_);
        multiply_truncate_uint_uint(plain_upper_half_threshold_.get(), coeff_div_plain_modulus_.get(), coeff_mod_count, upper_half_threshold_.get());

        // Calculate upper_half_increment.
        upper_half_increment_ = allocate_uint(coeff_mod_count, pool_);
        multiply_truncate_uint_uint(wide_plain_modulus.get(), coeff_div_plain_modulus_.get(), coeff_mod_count, upper_half_increment_.get());
        sub_uint_uint(context.coeff_modulus().pointer(), upper_half_increment_.get(), coeff_mod_count, upper_half_increment_.get());

        // Calculate coeff_modulus_ / 2.
        coeff_modulus_div_two_ = allocate_uint(coeff_mod_count, pool_);
        right_shift_uint(context.coeff_modulus().pointer(), 1, coeff_mod_count, coeff_modulus_div_two_.get());

        // Set the big coeff modulus for noise computation
        product_modulus_ = allocate_uint(coeff_mod_count, pool_);
        set_uint_uint(context.coeff_modulus_.pointer(), coeff_mod_count, product_modulus_.get());

        // Initialize moduli.
        mod_ = Modulus(product_modulus_.get(), coeff_mod_count);
        polymod_ = PolyModulus(parms_.poly_modulus().pointer(), coeff_count, poly_coeff_uint64_count);

		// Calculate map from Zmstar to generate representation
		int n = coeff_count - 1;
		int m = n << 1;
		int galois_elt;
		for (int i = 0; i < n / 2; i++)
		{
			galois_elt = (power(3, i)) % m;
			pair<int, int> temp_pair1 = { i, 0 };
			map_from_Zmstar_to_generator.emplace(galois_elt, temp_pair1);
			int galois_elt = (power(3, i) * (m - 1)) % m;
			pair<int, int> temp_pair2 = { i, 1 };
			map_from_Zmstar_to_generator.emplace(galois_elt, temp_pair2);
		}
    }

    RNSEvaluator::RNSEvaluator(const RNSEvaluator &copy) :
        pool_(copy.pool_), parms_(copy.parms_), qualifiers_(copy.qualifiers_), 
        base_convertor_(copy.base_convertor_), 
        coeff_small_ntt_tables_(copy.coeff_small_ntt_tables_),
        bsk_small_ntt_tables_(copy.bsk_small_ntt_tables_), 
        evaluation_keys_(copy.evaluation_keys_),
        coeff_mod_array_(copy.coeff_mod_array_), 
        bsk_mod_array_(copy.bsk_mod_array_),
        inv_coeff_products_mod_coeff_array_(copy.inv_coeff_products_mod_coeff_array_),
        bsk_base_mod_count_(copy.bsk_base_mod_count_)
    {
        int coeff_count = parms_.poly_modulus().coeff_count();
        int poly_coeff_uint64_count = parms_.poly_modulus().coeff_uint64_count();
        int coeff_mod_count = coeff_mod_array_.size();

        // Allocate memory and copy over values
        // Calculate (plain_modulus + 1) / 2 * coeff_div_plain_modulus.
        upper_half_threshold_ = allocate_uint(coeff_mod_count, pool_);
        set_uint_uint(copy.upper_half_threshold_.get(), coeff_mod_count, upper_half_threshold_.get());

        // Calculate upper_half_increment.
        upper_half_increment_ = allocate_uint(coeff_mod_count, pool_);
        set_uint_uint(copy.upper_half_increment_.get(), coeff_mod_count, upper_half_increment_.get());

        // Calculate coeff_modulus / plain_modulus.
        coeff_div_plain_modulus_ = allocate_uint(coeff_mod_count, pool_);
        set_uint_uint(copy.coeff_div_plain_modulus_.get(), coeff_mod_count, coeff_div_plain_modulus_.get());

        // Calculate (plain_modulus + 1) / 2.
        plain_upper_half_threshold_ = allocate_uint(coeff_mod_count, pool_);
        set_uint_uint(copy.plain_upper_half_threshold_.get(), coeff_mod_count, plain_upper_half_threshold_.get());

        // Calculate coeff_modulus - plain_modulus.
        plain_upper_half_increment_ = allocate_uint(coeff_mod_count, pool_);
        set_uint_uint(copy.plain_upper_half_increment_.get(), coeff_mod_count, plain_upper_half_increment_.get());

        // Calculate coeff_modulus_ / 2.
        coeff_modulus_div_two_ = allocate_uint(coeff_mod_count, pool_);
        set_uint_uint(copy.coeff_div_plain_modulus_.get(), coeff_mod_count, coeff_div_plain_modulus_.get());

        // Populate coeff products array for compose functions (used in noise budget)
        coeff_products_array_ = allocate_uint(coeff_mod_count * coeff_mod_count, pool_);
        set_uint_uint(copy.coeff_products_array_.get(), coeff_mod_count * coeff_mod_count, coeff_products_array_.get());

        // Set the big coeff modulus for noise computation
        product_modulus_ = allocate_uint(coeff_mod_count, pool_);
        set_uint_uint(copy.product_modulus_.get(), coeff_mod_count, product_modulus_.get());

        // Initialize moduli.
        mod_ = Modulus(product_modulus_.get(), coeff_mod_count);
        polymod_ = PolyModulus(parms_.poly_modulus().pointer(), coeff_count, poly_coeff_uint64_count);
    }

    void RNSEvaluator::negate(const Ciphertext &encrypted, Ciphertext &destination)
    {
        const BigPolyArray &encrypted_array = encrypted.get_array();
        BigPolyArray &destination_array = destination.get_mutable_array();

        // Extract encryption parameters.
        int coeff_count = parms_.poly_modulus().coeff_count();
        int coeff_mod_count = coeff_mod_array_.size();
        int coeff_bit_count = coeff_mod_count * bits_per_uint64;
        int count = encrypted_array.size();

        // Verify parameters.
        if (encrypted.rns_hash_block_ != parms_.get_hash_block())
        {
            throw invalid_argument("encrypted is not valid for encryption parameters");
        }

        ConstPointer encrypted_safe = duplicate_bigpolyarray_if_needed(encrypted_array, &encrypted == &destination, pool_);

        // Prepare destination
        destination_array.resize(count, coeff_count, coeff_bit_count);
        destination.rns_hash_block_ = parms_.get_hash_block();

        int encrypted_ptr_increment = coeff_count * coeff_mod_count;

        // Negate each bigpoly in the array
        for (int j = 0; j < count; j++)
        {
            for (int i = 0; i < coeff_mod_count; i++)
            {
                negate_poly_coeff_smallmod(encrypted_safe.get() + (i * coeff_count) + (j * encrypted_ptr_increment), coeff_count, coeff_mod_array_[i], coeff_mod_array_[i].uint64_count(), destination_array.pointer(0) + (i * coeff_count) + (j * encrypted_ptr_increment));
            }
        }
    }

    void RNSEvaluator::add(const Ciphertext &encrypted1, const Ciphertext &encrypted2, Ciphertext &destination)
    {
        const BigPolyArray &encrypted1_array = encrypted1.get_array();
        const BigPolyArray &encrypted2_array = encrypted2.get_array();
        BigPolyArray &destination_array = destination.get_mutable_array();

        // Extract encryption parameters.
        int coeff_count = parms_.poly_modulus().coeff_count();
        int coeff_mod_count = coeff_mod_array_.size();
        int coeff_bit_count = coeff_mod_count * bits_per_uint64;
        int encrypted1_count = encrypted1_array.size();
        int encrypted2_count = encrypted2_array.size();
        int max_count = max(encrypted1_count, encrypted2_count);
        int min_count = min(encrypted1_count, encrypted2_count);

        // Verify parameters.
        if (encrypted1.rns_hash_block_ != parms_.get_hash_block())
        {
            throw invalid_argument("encrypted1 is not valid for encryption parameters");
        }
        if (encrypted2.rns_hash_block_ != parms_.get_hash_block())
        {
            throw invalid_argument("encrypted2 is not valid for encryption parameters");
        }

        ConstPointer encrypted1_safe = duplicate_bigpolyarray_if_needed(encrypted1_array, &encrypted1 == &destination, pool_);
        ConstPointer encrypted2_safe = duplicate_bigpolyarray_if_needed(encrypted2_array, &encrypted2 == &destination, pool_);

        // Prepare destination
        destination_array.resize(max_count, coeff_count, coeff_bit_count);
        destination.rns_hash_block_ = parms_.get_hash_block();

        int encrypted_ptr_increment = coeff_count * coeff_mod_count;
        // Add BigPolyArrays
        for (int j = 0; j < min_count; j++)
        {
            for (int i = 0; i < coeff_mod_count; i++)
            {
                add_poly_poly_coeff_smallmod(encrypted1_safe.get() + (i * coeff_count) + (j * encrypted_ptr_increment), encrypted2_safe.get() + (i * coeff_count) + (j * encrypted_ptr_increment), coeff_count, coeff_mod_array_[i], destination_array.pointer(0) + (i * coeff_count) + (j * encrypted_ptr_increment));
            }
        }

        // If the sizes are the same, we are done
        if (max_count == min_count)
        {
            return;
        }

        // Copy the remainding BigPolys of the array with larger count into destination
        set_poly_poly((encrypted1_count == max_count ? encrypted1_safe.get() : encrypted2_safe.get()) + min_count * coeff_count * coeff_mod_count,
            coeff_count * (max_count - min_count), coeff_mod_count, destination_array.pointer(min_count));
    }

    void RNSEvaluator::add_many(const vector<Ciphertext> &encrypteds, Ciphertext &destination)
    {
        if (encrypteds.empty())
        {
            throw invalid_argument("encrypteds cannot be empty");
        }

        destination = encrypteds[0];
        for (size_t i = 1; i < encrypteds.size(); i++)
        {
            add(destination, encrypteds[i], destination);
        }
    }

    void RNSEvaluator::sub(const Ciphertext &encrypted1, const Ciphertext &encrypted2, Ciphertext &destination)
    {
        const BigPolyArray &encrypted1_array = encrypted1.get_array();
        const BigPolyArray &encrypted2_array = encrypted2.get_array();
        BigPolyArray &destination_array = destination.get_mutable_array();

        // Extract encryption parameters.
        int coeff_count = parms_.poly_modulus().coeff_count();
        int coeff_mod_count = coeff_mod_array_.size();
        int coeff_bit_count = coeff_mod_count * bits_per_uint64;
        int encrypted1_count = encrypted1_array.size();
        int encrypted2_count = encrypted2_array.size();
        int max_count = max(encrypted1_count, encrypted2_count);
        int min_count = min(encrypted1_count, encrypted2_count);

        // Verify parameters.
        if (encrypted1.rns_hash_block_ != parms_.get_hash_block())
        {
            throw invalid_argument("encrypted1 is not valid for encryption parameters");
        }
        if (encrypted2.rns_hash_block_ != parms_.get_hash_block())
        {
            throw invalid_argument("encrypted2 is not valid for encryption parameters");
        }

        ConstPointer encrypted1_safe = duplicate_bigpolyarray_if_needed(encrypted1_array, &encrypted1 == &destination, pool_);
        ConstPointer encrypted2_safe = duplicate_bigpolyarray_if_needed(encrypted2_array, &encrypted2 == &destination, pool_);

        // Prepare destination
        destination_array.resize(max_count, coeff_count, coeff_bit_count);
        destination.rns_hash_block_ = parms_.get_hash_block();

        // Subtract polynomials.
        for (int j = 0; j < min_count; j++)
        {
            for (int i = 0; i < coeff_mod_count; i++)
            {
                sub_poly_poly_coeff_smallmod(encrypted1_safe.get() + (i * coeff_count) + (j * coeff_count * coeff_mod_count), encrypted2_safe.get() + (i * coeff_count) + (j * coeff_count * coeff_mod_count), coeff_count, coeff_mod_array_[i], destination_array.pointer(0) + (i * coeff_count) + (j * coeff_count * coeff_mod_count));
            }
        }

        // If the sizes are the same, we are done
        if (max_count == min_count)
        {
            return;
        }

        // If encrypted1 has larger count, copy remaining entries into destination_array.
        if (encrypted1_count == max_count)
        {
            // Copy the remainding BigPolys of the array with larger count into destination
            set_poly_poly(encrypted1_safe.get() + min_count * coeff_count * coeff_mod_count, coeff_count * (max_count - min_count), coeff_mod_count, destination_array.pointer(min_count));
        }
        // Otherwise, encrypted2 has larger count, negate remaining entries and copy into destination_array.
        else
        {
            for (int i = 0; i < coeff_mod_count; i++)
            {
                negate_poly_coeff_smallmod(encrypted2_safe.get() + (i * coeff_count) + min_count * coeff_count * coeff_mod_count, coeff_count * (max_count - min_count), coeff_mod_array_[i], coeff_mod_array_[i].uint64_count(), destination_array.pointer(min_count) + (i * coeff_count));
            }
        }
    }

    void RNSEvaluator::multiply(const Ciphertext &encrypted1, const Ciphertext &encrypted2, Ciphertext &destination)
    {
        const BigPolyArray &encrypted1_array = encrypted1.get_array();
        const BigPolyArray &encrypted2_array = encrypted2.get_array();
        BigPolyArray &destination_array = destination.get_mutable_array();

        // Extract encryption parameters.
        int coeff_count = parms_.poly_modulus().coeff_count();
        int coeff_mod_count = coeff_mod_array_.size();
        int coeff_bit_count = coeff_mod_count * bits_per_uint64;
        int bsk_mtilde_count = bsk_base_mod_count_ + 1;

        int encrypted1_count = encrypted1_array.size();
        int encrypted2_count = encrypted2_array.size();

        // Determine destination_array.size()
        int dest_count = encrypted1_count + encrypted2_count - 1; // default is 3 (c_0, c_1, c_2)

                                                                  // Verify parameters.
        if (encrypted1.rns_hash_block_ != parms_.get_hash_block())
        {
            throw invalid_argument("encrypted1 is not valid for encryption parameters");
        }
        if (encrypted2.rns_hash_block_ != parms_.get_hash_block())
        {
            throw invalid_argument("encrypted2 is not valid for encryption parameters");
        }

        ConstPointer encrypted1_safe = duplicate_bigpolyarray_if_needed(encrypted1_array, &encrypted1 == &destination, pool_);
        ConstPointer encrypted2_safe = duplicate_bigpolyarray_if_needed(encrypted2_array, &encrypted2 == &destination, pool_);

        // Prepare destination
        destination_array.resize(dest_count, coeff_count, coeff_bit_count);
        destination.rns_hash_block_ = parms_.get_hash_block();
        destination_array.set_zero();

        int encrypted_ptr_increment = coeff_count * coeff_mod_count;
        int encrypted_bsk_mtilde_ptr_increment = coeff_count * bsk_mtilde_count;
        int encrypted_bsk_ptr_increment = coeff_count * bsk_base_mod_count_;
        // Make temp polys for FastBConvertor result from q ---> Bsk U {m_tilde}
        Pointer tmp_encrypted1_bsk_mtilde(allocate_poly(coeff_count * encrypted1_count, bsk_mtilde_count, pool_));
        Pointer tmp_encrypted2_bsk_mtilde(allocate_poly(coeff_count * encrypted2_count, bsk_mtilde_count, pool_));

        // Make temp polys for FastBConvertor result from Bsk U {m_tilde} -----> Bsk
        Pointer tmp_encrypted1_bsk(allocate_poly(coeff_count * encrypted1_count, bsk_base_mod_count_, pool_));
        Pointer tmp_encrypted2_bsk(allocate_poly(coeff_count * encrypted2_count, bsk_base_mod_count_, pool_));

        // Step 0: fast base convert from q to Bsk U {m_tilde}
        // Step 1: reduce q-overflows in Bsk
        // Iterate over all the ciphertexts inside encrypted1
        for (int i = 0; i < encrypted1_count; i++)
        {
            base_convertor_.fastbconv_mtilde(encrypted1_safe.get() + (i * encrypted_ptr_increment), tmp_encrypted1_bsk_mtilde.get() + (i * encrypted_bsk_mtilde_ptr_increment));
            base_convertor_.mont_rq(tmp_encrypted1_bsk_mtilde.get() + (i * encrypted_bsk_mtilde_ptr_increment), tmp_encrypted1_bsk.get() + (i * encrypted_bsk_ptr_increment));
        }
        
        // Iterate over all the ciphertexts inside encrypted2
        for (int i = 0; i < encrypted2_count; i++)
        {
            base_convertor_.fastbconv_mtilde(encrypted2_safe.get() + (i * encrypted_ptr_increment), tmp_encrypted2_bsk_mtilde.get() + (i * encrypted_bsk_mtilde_ptr_increment));
            base_convertor_.mont_rq(tmp_encrypted2_bsk_mtilde.get() + (i * encrypted_bsk_mtilde_ptr_increment), tmp_encrypted2_bsk.get() + (i * encrypted_bsk_ptr_increment));
        }
        
        // Step 2: compute product and multiply plain modulus to the result
        // We need to multiply both in q and Bsk. Values in encrypted_safe are in base q and values in tmp_encrypted_bsk are in base Bsk
        // We iterates over destination BigPolyArray and generates each BigPoly based on the indeces of inputs (arbitrary sizes for ciphertexts)
        // First allocate two temp polys: one for results in base q and the other for the result in base Bsk
        Pointer tmp_des_coeff_base(allocate_poly(coeff_count * dest_count, coeff_mod_count, pool_));
        Pointer tmp_des_bsk_base(allocate_poly(coeff_count * dest_count, bsk_base_mod_count_, pool_));

        // Reset the dest polys 
        set_zero_poly(coeff_count * dest_count, coeff_mod_count, tmp_des_coeff_base.get());
        set_zero_poly(coeff_count * dest_count, bsk_base_mod_count_, tmp_des_bsk_base.get());

        // Allocate two tmp polys: one for NTT multiplication results in base q and one for result in base Bsk
        Pointer tmp1_poly_coeff_base(allocate_poly(coeff_count, coeff_mod_count, pool_));
        Pointer tmp1_poly_bsk_base(allocate_poly(coeff_count, bsk_base_mod_count_, pool_));
        Pointer tmp2_poly_coeff_base(allocate_poly(coeff_count, coeff_mod_count, pool_));
        Pointer tmp2_poly_bsk_base(allocate_poly(coeff_count, bsk_base_mod_count_, pool_));

        int current_encrypted1_limit = 0;

        // First convert all the inputs into NTT form
        Pointer copy_encrypted1_ntt_coeff_mod(allocate_poly(coeff_count * encrypted1_count, coeff_mod_count, pool_));
        Pointer copy_encrypted1_ntt_bsk_base_mod(allocate_poly(coeff_count * encrypted1_count, bsk_base_mod_count_, pool_));
        Pointer copy_encrypted2_ntt_coeff_mod(allocate_poly(coeff_count * encrypted2_count, coeff_mod_count, pool_));
        Pointer copy_encrypted2_ntt_bsk_base_mod(allocate_poly(coeff_count * encrypted2_count, bsk_base_mod_count_, pool_));

        set_poly_poly(encrypted1_safe.get(), coeff_count * encrypted1_count, coeff_mod_count, copy_encrypted1_ntt_coeff_mod.get());
        set_poly_poly(tmp_encrypted1_bsk.get(), coeff_count * encrypted1_count, bsk_base_mod_count_, copy_encrypted1_ntt_bsk_base_mod.get());
        set_poly_poly(encrypted2_safe.get(), coeff_count * encrypted2_count, coeff_mod_count, copy_encrypted2_ntt_coeff_mod.get());
        set_poly_poly(tmp_encrypted2_bsk.get(), coeff_count * encrypted2_count, bsk_base_mod_count_, copy_encrypted2_ntt_bsk_base_mod.get());

        for (int i = 0; i < encrypted1_count; i++)
        {
            for (int j = 0; j < coeff_mod_count; j++)
            {
                smallntt_negacyclic_harvey_lazy(copy_encrypted1_ntt_coeff_mod.get() + (j * coeff_count) + (i * encrypted_ptr_increment), coeff_small_ntt_tables_[j], pool_);
            }
            for (int j = 0; j < bsk_base_mod_count_; j++)
            {
                smallntt_negacyclic_harvey_lazy(copy_encrypted1_ntt_bsk_base_mod.get() + (j * coeff_count) + (i * encrypted_bsk_ptr_increment), bsk_small_ntt_tables_[j], pool_);
            }
        }

        for (int i = 0; i < encrypted2_count; i++)
        {
            for (int j = 0; j < coeff_mod_count; j++)
            {
                smallntt_negacyclic_harvey_lazy(copy_encrypted2_ntt_coeff_mod.get() + (j * coeff_count) + (i * encrypted_ptr_increment), coeff_small_ntt_tables_[j], pool_);
            }
            for (int j = 0; j < bsk_base_mod_count_; j++)
            {
                smallntt_negacyclic_harvey_lazy(copy_encrypted2_ntt_bsk_base_mod.get() + (j * coeff_count) + (i * encrypted_bsk_ptr_increment), bsk_small_ntt_tables_[j], pool_);
            }
        }

        // Perform Karatsuba multiplication on size 2 ciphertexts
        if (encrypted1_count == 2 && encrypted2_count == 2)
        {
            Pointer tmp_first_mul_coeff_base(allocate_poly(coeff_count, coeff_mod_count, pool_));
            // Compute c0 + c1 and c0*d0 in base q
            for (int i = 0; i < coeff_mod_count; i++)
            {
                add_poly_poly_coeff_smallmod(copy_encrypted1_ntt_coeff_mod.get() + (i * coeff_count), copy_encrypted1_ntt_coeff_mod.get() + (i * coeff_count) + encrypted_ptr_increment, coeff_count, coeff_mod_array_[i], tmp1_poly_coeff_base.get() + (i * coeff_count));
                dyadic_product_coeff_smallmod(copy_encrypted1_ntt_coeff_mod.get() + (i * coeff_count), copy_encrypted2_ntt_coeff_mod.get() + (i * coeff_count), coeff_count, coeff_mod_array_[i], tmp_first_mul_coeff_base.get() + (i * coeff_count), pool_);
            }

            Pointer tmp_first_mul_bsk_base(allocate_poly(coeff_count, bsk_base_mod_count_, pool_));
            // Compute c0 + c1 and c0*d0 in base bsk
            for (int i = 0; i < bsk_base_mod_count_; i++)
            {
                add_poly_poly_coeff_smallmod(copy_encrypted1_ntt_bsk_base_mod.get() + (i * coeff_count), copy_encrypted1_ntt_bsk_base_mod.get() + (i * coeff_count) + encrypted_bsk_ptr_increment, coeff_count, bsk_mod_array_[i], tmp1_poly_bsk_base.get() + (i * coeff_count));
                dyadic_product_coeff_smallmod(copy_encrypted1_ntt_bsk_base_mod.get() + (i * coeff_count), copy_encrypted2_ntt_bsk_base_mod.get() + (i * coeff_count), coeff_count, bsk_mod_array_[i], tmp_first_mul_bsk_base.get() + (i * coeff_count), pool_);
            }

            Pointer tmp_second_mul_coeff_base(allocate_poly(coeff_count, coeff_mod_count, pool_));
            // Compute d0 + d1 and c1*d1 in base q
            for (int i = 0; i < coeff_mod_count; i++)
            {
                add_poly_poly_coeff_smallmod(copy_encrypted2_ntt_coeff_mod.get() + (i * coeff_count), copy_encrypted2_ntt_coeff_mod.get() + (i * coeff_count) + encrypted_ptr_increment, coeff_count, coeff_mod_array_[i], tmp2_poly_coeff_base.get() + (i * coeff_count));
                dyadic_product_coeff_smallmod(copy_encrypted1_ntt_coeff_mod.get() + (i * coeff_count) + encrypted_ptr_increment, copy_encrypted2_ntt_coeff_mod.get() + (i * coeff_count) + encrypted_ptr_increment, coeff_count, coeff_mod_array_[i], tmp_second_mul_coeff_base.get() + (i * coeff_count), pool_);
            }

            Pointer tmp_second_mul_bsk_base(allocate_poly(coeff_count, bsk_base_mod_count_, pool_));
            // Compute d0 + d1 and c1*d1 in base bsk
            for (int i = 0; i < bsk_base_mod_count_; i++)
            {
                add_poly_poly_coeff_smallmod(copy_encrypted2_ntt_bsk_base_mod.get() + (i * coeff_count), copy_encrypted2_ntt_bsk_base_mod.get() + (i * coeff_count) + encrypted_bsk_ptr_increment, coeff_count, bsk_mod_array_[i], tmp2_poly_bsk_base.get() + (i * coeff_count));
                dyadic_product_coeff_smallmod(copy_encrypted1_ntt_bsk_base_mod.get() + (i * coeff_count) + encrypted_bsk_ptr_increment, copy_encrypted2_ntt_bsk_base_mod.get() + (i * coeff_count) + encrypted_bsk_ptr_increment, coeff_count, bsk_mod_array_[i], tmp_second_mul_bsk_base.get() + (i * coeff_count), pool_);
            }

            Pointer tmp_mul_poly_coeff_base(allocate_poly(coeff_count, coeff_mod_count, pool_));
            Pointer tmp_mul_poly_bsk_base(allocate_poly(coeff_count, bsk_base_mod_count_, pool_));
            set_zero_poly(coeff_count, coeff_mod_count, tmp_mul_poly_coeff_base.get());
            set_zero_poly(coeff_count, bsk_base_mod_count_, tmp_mul_poly_bsk_base.get());

            // Set destination first and third polys in base q
            set_poly_poly(tmp_first_mul_coeff_base.get(), coeff_count, coeff_mod_count, tmp_des_coeff_base.get()); // Des[0] in base q
            set_poly_poly(tmp_second_mul_coeff_base.get(), coeff_count, coeff_mod_count, tmp_des_coeff_base.get() + 2 * encrypted_ptr_increment); // Des[2] in base q
            // Compute (c0 + c1)*(d0 + d1) - c0*d0 - c1*d1 in base q
            for (int i = 0; i < coeff_mod_count; i++)
            {
                dyadic_product_coeff_smallmod(tmp1_poly_coeff_base.get() + (i * coeff_count), tmp2_poly_coeff_base.get() + (i * coeff_count), coeff_count, coeff_mod_array_[i], tmp_mul_poly_coeff_base.get() + (i * coeff_count), pool_);
                sub_poly_poly_coeff_smallmod(tmp_mul_poly_coeff_base.get() + (i * coeff_count), tmp_first_mul_coeff_base.get() + (i * coeff_count), coeff_count, coeff_mod_array_[i], tmp_mul_poly_coeff_base.get() + (i * coeff_count));
                sub_poly_poly_coeff_smallmod(tmp_mul_poly_coeff_base.get() + (i * coeff_count), tmp_second_mul_coeff_base.get() + (i * coeff_count), coeff_count, coeff_mod_array_[i], tmp_des_coeff_base.get() + (i * coeff_count) + encrypted_ptr_increment); // Des[1] in base q
            }

            // Set destination first and third polys in base bsk
            set_poly_poly(tmp_first_mul_bsk_base.get(), coeff_count, bsk_base_mod_count_, tmp_des_bsk_base.get()); // Des[0] in base q
            set_poly_poly(tmp_second_mul_bsk_base.get(), coeff_count, bsk_base_mod_count_, tmp_des_bsk_base.get() + 2 * encrypted_bsk_ptr_increment); // Des[2] in base q
            // Compute (c0 + c1)*(d0 + d1)  - c0d0 - c1d1 in base bsk
            for (int i = 0; i < bsk_base_mod_count_; i++)
            {
                dyadic_product_coeff_smallmod(tmp1_poly_bsk_base.get() + (i * coeff_count), tmp2_poly_bsk_base.get() + (i * coeff_count), coeff_count, bsk_mod_array_[i], tmp_mul_poly_bsk_base.get() + (i * coeff_count), pool_);
                sub_poly_poly_coeff_smallmod(tmp_mul_poly_bsk_base.get() + (i * coeff_count), tmp_first_mul_bsk_base.get() + (i * coeff_count), coeff_count, bsk_mod_array_[i], tmp_mul_poly_bsk_base.get() + (i * coeff_count));
                sub_poly_poly_coeff_smallmod(tmp_mul_poly_bsk_base.get() + (i * coeff_count), tmp_second_mul_bsk_base.get() + (i * coeff_count), coeff_count, bsk_mod_array_[i], tmp_des_bsk_base.get() + (i * coeff_count) + encrypted_bsk_ptr_increment); // Des[1] in bsk
            }
        }
        else
        {
            // Perform multiplication on arbitrary size ciphertexts
            for (int secret_power_index = 0; secret_power_index < dest_count; ++secret_power_index)
            {
                // loop over encrypted1 components [i], seeing if a match exists with an encrypted2 component [j] such that [i+j]=[secret_power_index]
                // only need to check encrypted1 components up to and including [secret_power_index], and strictly less than [encrypted_array.size()]
                current_encrypted1_limit = min(encrypted1_count, secret_power_index + 1);

                for (int encrypted1_index = 0; encrypted1_index < current_encrypted1_limit; ++encrypted1_index)
                {
                    // check if a corresponding component in encrypted2 exists
                    if (encrypted2_count > secret_power_index - encrypted1_index)
                    {
                        int encrypted2_index = secret_power_index - encrypted1_index;

                        // reset temp_poly
                        set_zero_poly(coeff_count, coeff_mod_count, tmp1_poly_coeff_base.get());
                        set_zero_poly(coeff_count, bsk_base_mod_count_, tmp1_poly_bsk_base.get());

                        // NTT Multiplication and addition for results in q
                        for (int i = 0; i < coeff_mod_count; i++)
                        {
                            dyadic_product_coeff_smallmod(copy_encrypted1_ntt_coeff_mod.get() + (i * coeff_count) + (encrypted_ptr_increment * encrypted1_index), copy_encrypted2_ntt_coeff_mod.get() + (i * coeff_count) + (encrypted_ptr_increment * encrypted2_index), coeff_small_ntt_tables_[i].coeff_count(), coeff_small_ntt_tables_[i].smallmodulus(), tmp1_poly_coeff_base.get() + (i * coeff_count), pool_);
                            add_poly_poly_coeff_smallmod(tmp1_poly_coeff_base.get() + (i * coeff_count), tmp_des_coeff_base.get() + (i * coeff_count) + (secret_power_index * coeff_count * coeff_mod_count), coeff_count, coeff_mod_array_[i], tmp_des_coeff_base.get() + (i * coeff_count) + (secret_power_index * coeff_count * coeff_mod_count));
                        }

                        // NTT Multiplication and addition for results in Bsk
                        for (int i = 0; i < bsk_base_mod_count_; i++)
                        {
                            dyadic_product_coeff_smallmod(copy_encrypted1_ntt_bsk_base_mod.get() + (i * coeff_count) + (encrypted_bsk_ptr_increment * encrypted1_index), copy_encrypted2_ntt_bsk_base_mod.get() + (i * coeff_count) + (encrypted_bsk_ptr_increment * encrypted2_index), bsk_small_ntt_tables_[i].coeff_count(), bsk_mod_array_[i], tmp1_poly_bsk_base.get() + (i * coeff_count), pool_);
                            add_poly_poly_coeff_smallmod(tmp1_poly_bsk_base.get() + (i * coeff_count), tmp_des_bsk_base.get() + (i * coeff_count) + (secret_power_index * coeff_count * bsk_base_mod_count_), coeff_count, bsk_mod_array_[i], tmp_des_bsk_base.get() + (i * coeff_count) + (secret_power_index * coeff_count * bsk_base_mod_count_));
                        }
                    }
                }
            }
        }
        // Convert back outputs from NTT form
        for (int i = 0; i < dest_count; i++)
        {
            for (int j = 0; j < coeff_mod_count; j++)
            {
                inverse_smallntt_negacyclic_harvey(tmp_des_coeff_base.get() + (i * (encrypted_ptr_increment)) + (j * coeff_count), coeff_small_ntt_tables_[j], pool_);
            }
            for (int j = 0; j < bsk_base_mod_count_; j++)
            {
                inverse_smallntt_negacyclic_harvey(tmp_des_bsk_base.get() + (i * (encrypted_bsk_ptr_increment)) + (j * coeff_count), bsk_small_ntt_tables_[j], pool_);
            }
        }

        // Now we multiply plain modulus to both results in base q and Bsk and allocate them together in one container as (te0)q(te'0)Bsk | ... |te count)q (te' count)Bsk to make it ready for fast_floor 
        Pointer tmp_coeff_bsk_together(allocate_poly(coeff_count, dest_count * (coeff_mod_count + bsk_base_mod_count_), pool_));
        uint64_t *tmp_coeff_bsk_together_ptr = tmp_coeff_bsk_together.get();

        // Base q 
        for (int i = 0; i < dest_count; i++)
        {
            for (int j = 0; j < coeff_mod_count; j++)
            {
                multiply_poly_scalar_coeff_smallmod(tmp_des_coeff_base.get() + (j * coeff_count) + (i * encrypted_ptr_increment), coeff_count, parms_.plain_modulus().pointer(), coeff_mod_array_[j], tmp_coeff_bsk_together_ptr + (j * coeff_count));
            }
            tmp_coeff_bsk_together_ptr += encrypted_ptr_increment;
            
            for (int k = 0; k < bsk_base_mod_count_; k++)
            {
                multiply_poly_scalar_coeff_smallmod(tmp_des_bsk_base.get() + (k * coeff_count) + (i * encrypted_bsk_ptr_increment), coeff_count, parms_.plain_modulus().pointer(), bsk_mod_array_[k], tmp_coeff_bsk_together_ptr + (k * coeff_count));
            }
            tmp_coeff_bsk_together_ptr += encrypted_bsk_ptr_increment;
        }

        // Allocate a new poly for fast floor result in Bsk
        Pointer tmp_result_bsk(allocate_poly(coeff_count, dest_count * bsk_base_mod_count_, pool_));
        for (int i = 0; i < dest_count; i++)
        {
            // Step 3: fast floor from q U {Bsk} to Bsk 
            base_convertor_.fast_floor(tmp_coeff_bsk_together.get() + (i * (encrypted_ptr_increment + encrypted_bsk_ptr_increment)), tmp_result_bsk.get() + (i * encrypted_bsk_ptr_increment));

            // Step 4: fast base convert from Bsk to q
            base_convertor_.fastbconv_sk(tmp_result_bsk.get() + (i * encrypted_bsk_ptr_increment), destination_array.pointer(i));
        }
    }

    void RNSEvaluator::relinearize(const Ciphertext &encrypted, Ciphertext &destination, int destination_size)
    {
        const BigPolyArray &encrypted_array = encrypted.get_array();
        BigPolyArray &destination_array = destination.get_mutable_array();

        // Extract encryption parameters.
        int coeff_count = parms_.poly_modulus().coeff_count();
        int coeff_mod_count = coeff_mod_array_.size();
        int coeff_bit_count = coeff_mod_count * bits_per_uint64;
        int encrypted_count = encrypted_array.size();

        // Verify parameters.
        if (encrypted.rns_hash_block_ != parms_.get_hash_block())
        {
            throw invalid_argument("encrypted is not valid for encryption parameters");
        }
        if (destination_size < 2 || destination_size > encrypted_count)
        {
            throw invalid_argument("destination_size must be greater than or equal to 2 and less than or equal to current count");
        }
        if (evaluation_keys_.size() < encrypted_count - 2)
        {
            throw invalid_argument("not enough evaluation keys");
        }

        // If encrypted is already at the desired level, return
        if (destination_size == encrypted_count)
        {
            destination = encrypted;
            return;
        }

        ConstPointer encrypted_safe = duplicate_bigpolyarray_if_needed(encrypted_array, &encrypted == &destination, pool_);

        // Prepare destination
        destination_array.resize(destination_size, coeff_count, coeff_bit_count);
        destination.rns_hash_block_ = parms_.get_hash_block();

        // Relinearize one step at a time
        Pointer temp1(allocate_poly(encrypted_count * coeff_count, coeff_mod_count, pool_));
        Pointer temp2(allocate_poly((encrypted_count - 1) * coeff_count, coeff_mod_count, pool_));
        set_uint_uint(encrypted_safe.get(), encrypted_count * coeff_count * coeff_mod_count, temp1.get());

        // Calculate number of relinearize_one_step calls needed
        int relins_needed = encrypted_count - destination_size;

        // Update temp to store the current result after relinearization
        for (int i = 0; i < relins_needed; i++)
        {
            relinearize_one_step(temp1.get(), encrypted_count, temp2.get(), pool_);
            set_uint_uint(temp2.get(), (encrypted_count - 1) * coeff_count * coeff_mod_count, temp1.get());
            encrypted_count--;
        }

        // Put the output of final relinearization into destination_array.
        // We assume here that destination has been resized by the calling function.
        set_poly_poly(temp2.get(), destination_size * coeff_count, coeff_mod_count, destination_array.pointer(0));
    }

    void RNSEvaluator::relinearize_one_step(const uint64_t *encrypted, int encrypted_size, uint64_t *destination, MemoryPool &pool_)
    {
        // extract encryption parameters.
        int coeff_count = parms_.poly_modulus().coeff_count();
        int coeff_mod_count = coeff_mod_array_.size();
        int array_poly_uint64_count = coeff_count * coeff_mod_count;

        // encrypted_array[2], encrypted_array[3], ... encrypted_array[count - 2] all stay the same
        if (encrypted_size > 3)
        {
            set_poly_poly(encrypted + 2 * array_poly_uint64_count, (encrypted_size - 3) * coeff_count, coeff_mod_count, destination + 2 * array_poly_uint64_count);
        }

        const uint64_t *encrypted_coeff = encrypted + (encrypted_size - 1) * array_poly_uint64_count;
        Pointer encrypted_coeff_prod_inv_coeff(allocate_zero_uint(coeff_count, pool_));

        // decompose encrypted_array[count-1] into base w
        // want to create an array of polys, each of whose components i is (encrypted_array[count-1])^(i) - in the notation of FV paper
        // This allocation stores one of the decomposed factors modulo one of the primes
        Pointer decomp_encrypted_last(allocate_uint(coeff_count, pool_));

        // Lazy reduction   
        Pointer wide_innerresult0(allocate_zero_poly(coeff_count, 2 * coeff_mod_count, pool_));
        Pointer wide_innerresult1(allocate_zero_poly(coeff_count, 2 * coeff_mod_count, pool_));
        Pointer innerresult(allocate_zero_poly(coeff_count, coeff_mod_count, pool_));
        Pointer temp_decomp_coeff(allocate_uint(coeff_count, pool_));

        /*
        For lazy reduction to work here, we need to ensure that the 128-bit accumulators (wide_innerresult0 and wide_innerresult1)
        do not overflow. Since the modulus primes are at most 60 bits, if the total number of summands is K, then the size of the
        total sum of products (without reduction) is at most 62 + 60 + bit_length(K). We need this to be at most 128, thus we need
        bit_length(K) <= 6. Thus, we need K <= 63. In this case, this means sum_i evaluation_keys.get_data()[0][i].first.size() <= 63.
        */
        for (int i = 0; i < coeff_mod_count; i++)
        {
            multiply_poly_scalar_coeff_smallmod(encrypted_coeff + (i * coeff_count), coeff_count, &inv_coeff_products_mod_coeff_array_[i], coeff_mod_array_[i], encrypted_coeff_prod_inv_coeff.get());

            // Populate one poly at a time.
            // Create a polynomial to store the current decomposition value which will be copied into the array to populate it at the current index.
            int shift = 0;

            for (int k = 0; k < evaluation_keys_.get_data()[0][i].first.size(); k++)
            {
                // Decompose here
                for (int coeff_index = 0; coeff_index < coeff_count; coeff_index++)
                {
                    *(decomp_encrypted_last.get() + coeff_index) = *(encrypted_coeff_prod_inv_coeff.get() + coeff_index) >> shift;
                    *(decomp_encrypted_last.get() + coeff_index) &= (1ULL << parms_.decomposition_bit_count()) - 1;
                }

                for (int j = 0; j < coeff_mod_count; j++)
                {
                    set_uint_uint(decomp_encrypted_last.get(), coeff_count, temp_decomp_coeff.get());

                    // We don't reduce here, so might get up to two extra bits. Thus 62 bits at most.
                    smallntt_negacyclic_harvey_lazy(temp_decomp_coeff.get(), coeff_small_ntt_tables_[j], pool_);

                    // Lazy reduction
                    uint64_t wide_innerproduct[2];
                    for (int m = 0; m < coeff_count; m++)
                    {
                        multiply_uint64(*(temp_decomp_coeff.get() + m), *(evaluation_keys_.get_data()[encrypted_size - 3][i].first.pointer(k) + (m + j * coeff_count)), wide_innerproduct);
                        unsigned char carry = add_uint64(*(wide_innerresult0.get() + 2 * (m + j * coeff_count)), wide_innerproduct[0], 0, wide_innerresult0.get() + 2 * (m + j * coeff_count));
                        *(wide_innerresult0.get() + 2 * (m + j * coeff_count) + 1) += wide_innerproduct[1] + carry;
                    }
                    for (int m = 0; m < coeff_count; m++)
                    {
                        multiply_uint64(*(temp_decomp_coeff.get() + m), *(evaluation_keys_.get_data()[encrypted_size - 3][i].second.pointer(k) + (m + j * coeff_count)), wide_innerproduct);
                        unsigned char carry = add_uint64(*(wide_innerresult1.get() + 2 * (m + j * coeff_count)), wide_innerproduct[0], 0, wide_innerresult1.get() + 2 * (m + j * coeff_count));
                        *(wide_innerresult1.get() + 2 * (m + j * coeff_count) + 1) += wide_innerproduct[1] + carry;
                    }

                    //dyadic_product_coeff_smallmod(temp_decomp_coeff.get(), evaluation_keys.get_key(encrypted_size - 1)[i].first.pointer(k) + (j * coeff_count), coeff_count, coeff_mod_array_[j], innerproduct0.get() + (j * coeff_count), pool_);
                    //dyadic_product_coeff_smallmod(temp_decomp_coeff.get(), evaluation_keys.get_key(encrypted_size - 1)[i].second.pointer(k) + (j * coeff_count), coeff_count, coeff_mod_array_[j], innerproduct1.get() + (j * coeff_count), pool_);
                    //add_poly_poly_coeff_smallmod(innerresult0.get() + (j * coeff_count), innerproduct0.get() + (j * coeff_count), coeff_count, coeff_mod_array_[j], innerresult0.get() + (j * coeff_count));
                    //add_poly_poly_coeff_smallmod(innerresult1.get() + (j * coeff_count), innerproduct1.get() + (j * coeff_count), coeff_count, coeff_mod_array_[j], innerresult1.get() + (j * coeff_count));
                }

                shift += parms_.decomposition_bit_count();
            }
        }

        for (int i = 0; i < coeff_mod_count; ++i)
        {
            for (int m = 0; m < coeff_count; m++)
            {
                barrett_reduce_128(wide_innerresult0.get() + 2 * (m + i * coeff_count), coeff_mod_array_[i], innerresult.get() + m + (i * coeff_count));
            }
            inverse_smallntt_negacyclic_harvey(innerresult.get() + (i * coeff_count), coeff_small_ntt_tables_[i], pool_);
            add_poly_poly_coeff_smallmod(encrypted + (i * coeff_count), innerresult.get() + (i * coeff_count), coeff_count, coeff_mod_array_[i], destination + (i * coeff_count));

            for (int m = 0; m < coeff_count; m++)
            {
                barrett_reduce_128(wide_innerresult1.get() + 2 * (m + i * coeff_count), coeff_mod_array_[i], innerresult.get() + m + (i * coeff_count));
            }
            inverse_smallntt_negacyclic_harvey(innerresult.get() + (i * coeff_count), coeff_small_ntt_tables_[i], pool_);
            add_poly_poly_coeff_smallmod(encrypted + (i * coeff_count) + array_poly_uint64_count, innerresult.get() + (i * coeff_count), coeff_count, coeff_mod_array_[i], destination + (i * coeff_count) + array_poly_uint64_count);
        }

        //for (int i = 0; i < coeff_mod_count; i++)
        //{
        //    inverse_smallntt_negacyclic_harvey(innerresult0.get() + (i * coeff_count), coeff_small_ntt_tables_[i], pool_);
        //    add_poly_poly_coeff_smallmod(encrypted + (i * coeff_count), innerresult0.get() + (i * coeff_count), coeff_count, coeff_mod_array_[i], destination + (i * coeff_count));
        //}
        //for (int i = 0; i < coeff_mod_count; i++)
        //{
        //    inverse_smallntt_negacyclic_harvey(innerresult1.get() + (i * coeff_count), coeff_small_ntt_tables_[i], pool_);
        //    add_poly_poly_coeff_smallmod(encrypted + (i * coeff_count) + array_poly_uint64_count, innerresult1.get() + (i * coeff_count), coeff_count, coeff_mod_array_[i], destination + (i * coeff_count) + array_poly_uint64_count);
        //}
    }



    Ciphertext RNSEvaluator::multiply_many(vector<Ciphertext> encrypteds)
    {
        size_t original_size = encrypteds.size();

        // Verify parameters.
        if (original_size == 0)
        {
            throw invalid_argument("encrypteds vector must not be empty");
        }

        // If there is only one ciphertext, return it after checking validity.
        if (original_size == 1)
        {
            // Verify parameters.
            if (encrypteds[0].rns_hash_block_ != parms_.get_hash_block())
            {
                throw invalid_argument("encrypteds is not valid for encryption parameters");
            }
            return encrypteds[0];
        }

        // Repeatedly multiply and add to the back of the vector until the end is reached
        for (size_t i = 0; i < encrypteds.size() - 1; i += 2)
        {
            encrypteds.emplace_back(relinearize(multiply(encrypteds[i], encrypteds[i + 1])));
        }

        return encrypteds[encrypteds.size() - 1];
    }

    void RNSEvaluator::multiply_plain_ntt(const Ciphertext &encrypted_ntt, const Plaintext &plain_ntt, Ciphertext &destination_ntt)
    {
        const BigPolyArray &encrypted_ntt_array = encrypted_ntt.get_array();
        const BigPoly &plain_ntt_poly = plain_ntt.get_poly();
        BigPolyArray &destination_ntt_array = destination_ntt.get_mutable_array();

        // Do the encryption parameters support NTT?
        if (!qualifiers_.enable_ntt)
        {
            throw logic_error("encryption parameters do not support NTT");
        }

        // Extract encryption parameters.
        int coeff_count = parms_.poly_modulus().coeff_count();
        int coeff_mod_count = coeff_mod_array_.size();
        int coeff_bit_count = coeff_mod_count * bits_per_uint64;
        int encrypted_count = encrypted_ntt.size();

        // Verify parameters.
        if (encrypted_ntt.rns_hash_block_ != parms_.get_hash_block())
        {
            throw invalid_argument("encrypted_ntt is not valid for encryption parameters");
        }
        if (plain_ntt_poly.is_zero())
        {
            throw invalid_argument("plain cannot be zero");
        }

        // Prepare destination
        destination_ntt_array.resize(encrypted_count, coeff_count, coeff_bit_count);
        destination_ntt.rns_hash_block_ = parms_.get_hash_block();

        Pointer temp(allocate_poly(coeff_count, coeff_mod_count, pool_));
        for (int i = 0; i < encrypted_count; i++)
        {
            set_poly_poly(encrypted_ntt_array.pointer(i), coeff_count, coeff_mod_count, temp.get());
            for (int j = 0; j < coeff_mod_count; j++)
            {
                dyadic_product_coeff_smallmod(temp.get() + (j * coeff_count), plain_ntt_poly.pointer() + (j * coeff_count) , coeff_count - 1, coeff_mod_array_[j], destination_ntt_array.pointer(i) + (j * coeff_count), pool_);
            }
        }
    }

    void RNSEvaluator::exponentiate(const Ciphertext &encrypted, uint64_t exponent, Ciphertext &destination)
    {
        // Verify parameters.
        if (encrypted.rns_hash_block_ != parms_.get_hash_block())
        {
            throw invalid_argument("encrypted is not valid for encryption parameters");
        }
        if (exponent == 0)
        {
            throw invalid_argument("exponent cannot be 0");
        }

        if (exponent == 1)
        {
            destination = encrypted;
            return;
        }

        vector<Ciphertext> exp_vector(exponent, encrypted);
        destination = multiply_many(exp_vector);
    }

    void RNSEvaluator::rns_decompose(uint64_t *value)
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

        // Copy value pointer
        uint64_t *value_ptr = value;

        for (int i = 0; i < coeff_count; i++)
        {
            for (int j = 0; j < coeff_mod_count; j++)
            {
                uint64_t reduced_value;
                small_modulo_uint(value_ptr, coeff_mod_count, coeff_mod_array_[j], &reduced_value, pool_);
                *(coefficients.get() + i + (j * coeff_count)) = reduced_value;
            }
            value_ptr += coeff_mod_count;
        }
        set_uint_uint(coefficients.get(), total_uint64_count, value);
    }

    void RNSEvaluator::rns_preencrypt(const uint64_t *plain, int plain_coeff_count, int plain_coeff_uint64_count, uint64_t *destination)
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

    void RNSEvaluator::add_plain(const Ciphertext &encrypted, const Plaintext &plain, Ciphertext &destination)
    {
        const BigPolyArray &encrypted_array = encrypted.get_array();
        const BigPoly &plain_poly = plain.get_poly();
        BigPolyArray &destination_array = destination.get_mutable_array();

        // Extract encryption parameters.
        int coeff_count = parms_.poly_modulus().coeff_count();
        int coeff_mod_count = coeff_mod_array_.size();
        int coeff_bit_count = coeff_mod_count * bits_per_uint64;
        int encrypted_count = encrypted_array.size();
        int plain_coeff_uint64_count = divide_round_up(plain_poly.coeff_bit_count(), bits_per_uint64);

        // Verify parameters.
        if (encrypted.rns_hash_block_ != parms_.get_hash_block())
        {
            throw invalid_argument("encrypted is not valid for encryption parameters");
        }
#ifdef _DEBUG
        if (plain_poly.significant_coeff_count() >= coeff_count || !are_poly_coefficients_less_than(plain_poly.pointer(), plain_poly.coeff_count(), plain_poly.coeff_uint64_count(), parms_.plain_modulus().pointer(), parms_.plain_modulus().uint64_count()))
        {
            throw invalid_argument("plain is too large to be represented by encryption parameters");
        }
#endif
        ConstPointer encrypted_safe = duplicate_bigpolyarray_if_needed(encrypted_array, &encrypted == &destination, pool_);

        // Prepare destination
        destination_array.resize(encrypted_count, coeff_count, coeff_bit_count);
        destination.rns_hash_block_ = parms_.get_hash_block();

        // Multiply plain by scalar coeff_div_plaintext and reposition if in upper-half, store in destination_array[0]
        rns_preencrypt(plain_poly.pointer(), plain_poly.coeff_count(), plain_coeff_uint64_count, destination_array.pointer(0));

        // Add preencrypted-version of plain with encrypted_array[0], store in destination_array[0].
        for (int i = 0; i < coeff_mod_count; i++)
        {
            add_poly_poly_coeff_smallmod(encrypted_safe.get() + (i * coeff_count), destination_array.pointer(0) + (i * coeff_count), coeff_count, coeff_mod_array_[i], destination_array.pointer(0) + (i * coeff_count));
        }
        
        // Set the remainder of destination to be as in encrypted_array.
        set_poly_poly(encrypted_safe.get() + coeff_count * coeff_mod_count, coeff_count * (encrypted_count - 1), coeff_mod_count, destination_array.pointer(1));
    }

    void RNSEvaluator::sub_plain(const Ciphertext &encrypted, const Plaintext &plain, Ciphertext &destination)
    {
        const BigPolyArray &encrypted_array = encrypted.get_array();
        const BigPoly &plain_poly = plain.get_poly();
        BigPolyArray &destination_array = destination.get_mutable_array();

        // Extract encryption parameters.
        int coeff_count = parms_.poly_modulus().coeff_count();
        int coeff_mod_count = coeff_mod_array_.size();
        int coeff_bit_count = coeff_mod_count * bits_per_uint64;
        int encrypted_count = encrypted_array.size();
        int plain_coeff_uint64_count = divide_round_up(plain_poly.coeff_bit_count(), bits_per_uint64);

        // Verify parameters.
        if (encrypted.rns_hash_block_ != parms_.get_hash_block())
        {
            throw invalid_argument("encrypted is not valid for encryption parameters");
        }
#ifdef _DEBUG
        if (plain_poly.significant_coeff_count() >= coeff_count || !are_poly_coefficients_less_than(plain_poly.pointer(), plain_poly.coeff_count(), plain_poly.coeff_uint64_count(), parms_.plain_modulus().pointer(), parms_.plain_modulus().uint64_count()))
        {
            throw invalid_argument("plain is too large to be represented by encryption parameters");
        }
#endif
        ConstPointer encrypted_safe = duplicate_bigpolyarray_if_needed(encrypted_array, &encrypted == &destination, pool_);

        // Prepare destination
        destination_array.resize(encrypted_count, coeff_count, coeff_bit_count);
        destination.rns_hash_block_ = parms_.get_hash_block();

        // Multiply plain by scalar coeff_div_plaintext and reposition if in upper-half, store in destination_array[0]
        rns_preencrypt(plain_poly.pointer(), plain_poly.coeff_count(), plain_coeff_uint64_count, destination_array.pointer(0));

        // Subtract preencrypted-version of plain2 with encrypted_array[0], store in destination_array[0].
        for (int i = 0; i < coeff_mod_count; i++)
        {
            sub_poly_poly_coeff_smallmod(encrypted_safe.get() + (i * coeff_count), destination_array.pointer(0) + (i * coeff_count), coeff_count, coeff_mod_array_[i], destination_array.pointer(0) + (i * coeff_count));
        }

        // Set the remainder of destination to be as in encrypted_array.
        set_poly_poly(encrypted_safe.get() + coeff_count * coeff_mod_count, coeff_count * (encrypted_count - 1), coeff_mod_count, destination_array.pointer(1));
    }

    
	void RNSEvaluator::multiply_plain(const Ciphertext &encrypted, const Plaintext &plain, Ciphertext &destination)
	{
		const BigPolyArray &encrypted_array = encrypted.get_array();
		const BigPoly &plain_poly = plain.get_poly();
		BigPolyArray &destination_array = destination.get_mutable_array();

		// Extract encryption parameters.
		int coeff_count = parms_.poly_modulus().coeff_count();
		int coeff_mod_count = coeff_mod_array_.size();
		int coeff_bit_count = coeff_mod_count * bits_per_uint64;
		int encrypted_count = encrypted_array.size();

		// Verify parameters.
		if (encrypted.rns_hash_block_ != parms_.get_hash_block())
		{
			throw invalid_argument("encrypted is not valid for encryption parameters");
		}
		if (plain_poly.is_zero())
		{
			throw invalid_argument("plain cannot be zero");
		}
#ifdef _DEBUG
		if (plain_poly.significant_coeff_count() >= coeff_count || !are_poly_coefficients_less_than(plain_poly.pointer(), plain_poly.coeff_count(), plain_poly.coeff_uint64_count(), parms_.plain_modulus().pointer(), parms_.plain_modulus().uint64_count()))
		{
			throw invalid_argument("plain is too large to be represented by encryption parameters");
		}
#endif
		// Prepare destination
		destination_array.resize(encrypted_count, coeff_count, coeff_bit_count);
		destination.rns_hash_block_ = parms_.get_hash_block();

		// Multiplying just by a constant?
		int plain_coeff_count = min(plain_poly.coeff_count(), coeff_count);
		int plain_coeff_uint64_count = plain_poly.coeff_uint64_count();
		const uint64_t *plain_coeff = plain_poly.pointer();
		if (plain_coeff_count == 1)
		{
			if (!qualifiers_.enable_fast_plain_lift)
			{
				Pointer moved2ptr(allocate_uint(coeff_mod_count, pool_));
				uint64_t *moved2_coeff = moved2ptr.get();
				bool is_upper_half = is_greater_than_or_equal_uint_uint(plain_coeff, plain_coeff_uint64_count, plain_upper_half_threshold_.get(), coeff_mod_count);
				if (is_upper_half)
				{
					add_uint_uint(plain_coeff, plain_coeff_uint64_count, plain_upper_half_increment_.get(), coeff_mod_count, 0, coeff_mod_count, moved2_coeff);
				}
				else
				{
					set_uint_uint(plain_coeff, plain_coeff_uint64_count, coeff_mod_count, moved2_coeff);
				}

				rns_decompose(moved2_coeff);

				for (int i = 0; i < encrypted_count; i++)
				{
					for (int j = 0; j < coeff_mod_count; j++)
					{
						multiply_poly_scalar_coeff_smallmod(encrypted_array.pointer(i) + (j * coeff_count), coeff_count, moved2_coeff + (j * coeff_count), coeff_mod_array_[j], destination_array.pointer(i) + (j * coeff_count));
					}
				}
				return;
			}
			else
			{
				Pointer moved2ptr(allocate_uint(coeff_count * coeff_mod_count, pool_));
				set_zero_uint(coeff_count * coeff_mod_count, moved2ptr.get());

				// Need for lift plain coefficient in RNS form regarding to each qi
				if (*(plain_coeff) >= *(plain_upper_half_threshold_.get()))
				{
					for (int j = 0; j < coeff_mod_count; j++)
					{
						add_uint_uint_smallmod(plain_coeff, &plain_upper_half_increment_array_[j], coeff_mod_array_[j], moved2ptr.get() + (j * coeff_count));
					}
				}
				// No need for lifting
				else
				{
					for (int j = 0; j < coeff_mod_count; j++)
					{
						*(moved2ptr.get() + (j * coeff_count)) = *(plain_coeff);
					}
				}

				for (int i = 0; i < encrypted_count; i++)
				{
					for (int j = 0; j < coeff_mod_count; j++)
					{
						multiply_poly_scalar_coeff_smallmod(encrypted_array.pointer(i) + (j * coeff_count), coeff_count, moved2ptr.get() + (j * coeff_count), coeff_mod_array_[j], destination_array.pointer(i) + (j * coeff_count));
					}
				}
				return;
			}
		}

		if (!qualifiers_.enable_fast_plain_lift)
		{
			// Reposition coefficients.
			Pointer moved2ptr(allocate_poly(coeff_count, coeff_mod_count, pool_));
			uint64_t *moved2_coeff = moved2ptr.get();
			for (int i = 0; i < plain_coeff_count; i++)
			{
				bool is_upper_half = is_greater_than_or_equal_uint_uint(plain_coeff, plain_coeff_uint64_count, plain_upper_half_threshold_.get(), coeff_mod_count);
				if (is_upper_half)
				{
					add_uint_uint(plain_coeff, plain_coeff_uint64_count, plain_upper_half_increment_.get(), coeff_mod_count, 0, coeff_mod_count, moved2_coeff);
				}
				else
				{
					set_uint_uint(plain_coeff, plain_coeff_uint64_count, coeff_mod_count, moved2_coeff);
				}
				moved2_coeff += coeff_mod_count;
				plain_coeff += plain_coeff_uint64_count;
			}
			set_zero_poly(coeff_count - plain_coeff_count, coeff_mod_count, moved2_coeff);

			rns_decompose(moved2ptr.get());
			// Need to multiply each component in encrypted with moved2ptr (plain poly)
			if (qualifiers_.enable_ntt)
			{
				// Transform plain poly only once
				for (int i = 0; i < coeff_mod_count; i++)
				{
					smallntt_negacyclic_harvey(moved2ptr.get() + (i * coeff_count), coeff_small_ntt_tables_[i], pool_);
				}

				for (int i = 0; i < encrypted_count; i++)
				{
					for (int j = 0; j < coeff_mod_count; j++)
					{
						smallntt_multiply_poly_nttpoly(encrypted_array.pointer(i) + (j * coeff_count), moved2ptr.get() + (j * coeff_count), coeff_small_ntt_tables_[j], destination_array.pointer(i) + (j * coeff_count), pool_);
					}
				}
			}
			else
			{
				// This branch should never be reached
				throw logic_error("invalid encryption parameters");
			}
		}
		else
		{
			Pointer moved2ptr(allocate_poly(coeff_count, coeff_mod_count, pool_));
			set_zero_uint(coeff_count * coeff_mod_count, moved2ptr.get());

			for (int i = 0; i < plain_coeff_count; i++)
			{
				// Need to lift the coefficient in each qi
				if (*(plain_coeff + i) >= *(plain_upper_half_threshold_.get()))
				{
					for (int j = 0; j < coeff_mod_count; j++)
					{
						add_uint_uint_smallmod(plain_coeff + i, &plain_upper_half_increment_array_[j], coeff_mod_array_[j], moved2ptr.get() + i + (j * coeff_count));
					}
				}
				// No need for lifting
				else
				{
					for (int j = 0; j < coeff_mod_count; j++)
					{
						*(moved2ptr.get() + i + (j * coeff_count)) = *(plain_coeff + i);
					}
				}
			}
			// Need to multiply each component in encrypted with moved2ptr (plain poly)
			if (qualifiers_.enable_ntt)
			{
				// Transform plain poly only once
				for (int i = 0; i < coeff_mod_count; i++)
				{
					smallntt_negacyclic_harvey(moved2ptr.get() + (i * coeff_count), coeff_small_ntt_tables_[i], pool_);
				}

				for (int i = 0; i < encrypted_count; i++)
				{
					for (int j = 0; j < coeff_mod_count; j++)
					{
						smallntt_multiply_poly_nttpoly(encrypted_array.pointer(i) + (j * coeff_count), moved2ptr.get() + (j * coeff_count), coeff_small_ntt_tables_[j], destination_array.pointer(i) + (j * coeff_count), pool_);
					}
				}
			}
			else
			{
				// This branch should never be reached
				throw logic_error("invalid encryption parameters");
			}
		}
	}

	void RNSEvaluator::transform_to_ntt(Plaintext &plain)
	{
		BigPoly &plain_poly = plain.get_poly();

		// Do the encryption parameters support NTT?
		if (!qualifiers_.enable_ntt)
		{
			throw logic_error("encryption parameters do not support NTT");
		}

		// Extract encryption parameters.
		int coeff_count = parms_.poly_modulus().coeff_count();
		int coeff_mod_count = coeff_mod_array_.size();
		int coeff_bit_count = coeff_mod_count * bits_per_uint64;
		int plain_coeff_count = min(plain_poly.coeff_count(), coeff_count);
		int plain_coeff_uint64_count = plain_poly.coeff_uint64_count();


		// Verify parameters.
		if (plain_poly.coeff_count() > coeff_count || plain_poly.coeff_bit_count() > coeff_bit_count)
		{
			throw invalid_argument("plain is not valid for encryption parameters");
		}
#ifdef _DEBUG
		if (plain_poly.significant_coeff_count() >= coeff_count || !are_poly_coefficients_less_than(plain_poly.pointer(), plain_poly.coeff_count(), plain_poly.coeff_uint64_count(), parms_.plain_modulus().pointer(), parms_.plain_modulus().uint64_count()))
		{
			throw invalid_argument("plain is too large to be represented by encryption parameters");
		}
#endif
		// Verify if plain lift is needed
		if (!qualifiers_.enable_fast_plain_lift)
		{
			const uint64_t *plain_coeff = plain_poly.pointer();
			Pointer moved2ptr(allocate_poly(coeff_count, coeff_mod_count, pool_));
			uint64_t *moved2_coeff = moved2ptr.get();
			for (int i = 0; i < plain_coeff_count; i++)
			{
				bool is_upper_half = is_greater_than_or_equal_uint_uint(plain_coeff, plain_coeff_uint64_count, plain_upper_half_threshold_.get(), coeff_mod_count);
				if (is_upper_half)
				{
					add_uint_uint(plain_coeff, plain_coeff_uint64_count, plain_upper_half_increment_.get(), coeff_mod_count, 0, coeff_mod_count, moved2_coeff);
				}
				else
				{
					set_uint_uint(plain_coeff, plain_coeff_uint64_count, coeff_mod_count, moved2_coeff);
				}
				moved2_coeff += coeff_mod_count;
				plain_coeff += plain_coeff_uint64_count;
			}
			set_zero_poly(coeff_count - plain_coeff_count, coeff_mod_count, moved2_coeff);

			rns_decompose(moved2ptr.get());

			// Transform to NTT domain
			for (int i = 0; i < coeff_mod_count; i++)
			{
				smallntt_negacyclic_harvey(moved2ptr.get() + (i * coeff_count), coeff_small_ntt_tables_[i], pool_);
			}
			plain_poly.resize(coeff_count, coeff_bit_count);
			set_poly_poly(moved2ptr.get(), coeff_count, coeff_mod_count, plain_poly.pointer());
		}
		// No need for composed plain lift and rns_decomposition
		else
		{
			uint64_t *plain_coeff = plain_poly.pointer();
			Pointer moved2ptr(allocate_poly(coeff_count, coeff_mod_count, pool_));
			set_zero_uint(coeff_count * coeff_mod_count, moved2ptr.get());

			for (int i = 0; i < plain_coeff_count; i++)
			{
				// Need to lift the coefficient in each qi
				if (*(plain_coeff + i) >= *(plain_upper_half_threshold_.get()))
				{
					for (int j = 0; j < coeff_mod_count; j++)
					{
						add_uint_uint_smallmod(plain_coeff + i, &plain_upper_half_increment_array_[j], coeff_mod_array_[j], moved2ptr.get() + i + (j * coeff_count));
					}
				}
				// No need for lifting
				else
				{
					for (int j = 0; j < coeff_mod_count; j++)
					{
						*(moved2ptr.get() + i + (j * coeff_count)) = *(plain_coeff + i);
					}
				}
			}

			// Transform to NTT domain
			for (int i = 0; i < coeff_mod_count; i++)
			{
				smallntt_negacyclic_harvey(moved2ptr.get() + (i * coeff_count), coeff_small_ntt_tables_[i], pool_);
			}
			plain_poly.resize(coeff_count, coeff_bit_count);
			set_poly_poly(moved2ptr.get(), coeff_count, coeff_mod_count, plain_poly.pointer());
		}
	}

    void RNSEvaluator::transform_to_ntt(Ciphertext &encrypted)
    {
        BigPolyArray &encrypted_array = encrypted.get_mutable_array();

        // Do the encryption parameters support NTT?
        if (!qualifiers_.enable_ntt)
        {
            throw logic_error("encryption parameters do not support NTT");
        }

        // Extract encryption parameters.
        int coeff_count = parms_.poly_modulus().coeff_count();
        int coeff_mod_count = coeff_mod_array_.size();

        // Verify parameters.
        if (encrypted.rns_hash_block_ != parms_.get_hash_block())
        {
            throw invalid_argument("encrypted is not valid for encryption parameters");
        }

        // Transform each polynomial to NTT domain
        for (int i = 0; i < encrypted_array.size(); i++)
        {
            for (int j = 0; j < coeff_mod_count; j++)
            {
                smallntt_negacyclic_harvey(encrypted_array.pointer(i) + (j * coeff_count), coeff_small_ntt_tables_[j], pool_);
            }
        }
    }

    void RNSEvaluator::transform_from_ntt(Plaintext &plain_ntt)
    {
        BigPoly &plain_ntt_poly = plain_ntt.get_poly();

        // Do the encryption parameters support NTT?
        if (!qualifiers_.enable_ntt)
        {
            throw logic_error("encryption parameters do not support NTT");
        }

        // Extract encryption parameters.
        int coeff_count = parms_.poly_modulus().coeff_count();
        int coeff_mod_count = coeff_mod_array_.size();
        int coeff_bit_count = coeff_mod_count * bits_per_uint64;

        // Verify parameters.
        if (plain_ntt_poly.coeff_count() != coeff_count || plain_ntt_poly.coeff_bit_count() != coeff_bit_count)
        {
            throw invalid_argument("plain_ntt is not valid for encryption parameters");
        }

        // Transform back to polynomial domain
        for (int i = 0; i < coeff_mod_count; i++)
        {
            inverse_smallntt_negacyclic_harvey(plain_ntt_poly.pointer() + (i * coeff_count), coeff_small_ntt_tables_[i], pool_);
        }

        uint64_t *plain_coeff = plain_ntt_poly.pointer();

        rns_compose(plain_coeff);
        
        for (int i = 0; i < coeff_count; i++)
        {
            bool is_upper_half = is_greater_than_or_equal_uint_uint(plain_coeff, upper_half_threshold_.get(), coeff_mod_count);
            if (is_upper_half)
            {
                sub_uint_uint(plain_coeff, plain_upper_half_increment_.get(), coeff_mod_count, plain_coeff);
            }
            plain_coeff += coeff_mod_count;
        }
        
        // Verify output.
#ifdef _DEBUG
        if (plain_ntt_poly.significant_coeff_count() >= coeff_count || !are_poly_coefficients_less_than(plain_ntt_poly.pointer(), plain_ntt_poly.coeff_count(), plain_ntt_poly.coeff_uint64_count(), parms_.plain_modulus().pointer(), parms_.plain_modulus().uint64_count()))
        {
            throw invalid_argument("plain_ntt is not valid for encryption parameters");
        }
#endif

    }

    void RNSEvaluator::transform_from_ntt(Ciphertext &encrypted_ntt)
    {
        BigPolyArray &encrypted_ntt_array = encrypted_ntt.get_mutable_array();

        // Do the encryption parameters support NTT?
        if (!qualifiers_.enable_ntt)
        {
            throw logic_error("encryption parameters do not support NTT");
        }

        // Extract encryption parameters.
        int coeff_count = parms_.poly_modulus().coeff_count();
        int coeff_mod_count = coeff_mod_array_.size();

        // Verify parameters.
        if (encrypted_ntt.rns_hash_block_ != parms_.get_hash_block())
        {
            throw invalid_argument("encrypted_ntt is not valid for encryption parameters");
        }

        // Transform each polynomial from NTT domain
        for (int i = 0; i < encrypted_ntt.size(); i++)
        {
            for (int j = 0; j < coeff_mod_count; j++)
            {
                inverse_smallntt_negacyclic_harvey(encrypted_ntt_array.pointer(i) + (j * coeff_count), coeff_small_ntt_tables_[j], pool_);
            }
        }
    }

    void RNSEvaluator::square(const Ciphertext &encrypted, Ciphertext &destination)
    {
        const BigPolyArray &encrypted_array = encrypted.get_array();
        BigPolyArray &destination_array = destination.get_mutable_array();

        int encrypted_count = encrypted_array.size();

        // Optimization implemented currently only for size 2 ciphertexts
        if (encrypted_count != 2)
        {
            multiply(encrypted, encrypted, destination);
            return;
        }

        int coeff_count = parms_.poly_modulus().coeff_count();
        int coeff_mod_count = coeff_mod_array_.size();
        int bsk_mtilde_count = bsk_base_mod_count_ + 1;
        int coeff_bit_count = coeff_mod_count * bits_per_uint64;
        int encrypted_ptr_increment = coeff_count * coeff_mod_count;
        int encrypted_bsk_mtilde_ptr_increment = coeff_count * bsk_mtilde_count;
        int encrypted_bsk_ptr_increment = coeff_count * bsk_base_mod_count_;

        // Determine destination_array.size()
        int dest_count = (encrypted_count << 1) - 1;

        // Verify parameters.
        if (encrypted.rns_hash_block_ != parms_.get_hash_block())
        {
            throw invalid_argument("encrypted is not valid for encryption parameters");
        }

        ConstPointer encrypted_safe = duplicate_bigpolyarray_if_needed(encrypted_array, &encrypted == &destination, pool_);

        // Prepare destination
        destination_array.resize(dest_count, coeff_count, coeff_bit_count);
        destination.rns_hash_block_ = parms_.get_hash_block();
        destination_array.set_zero();

        // Make temp poly for FastBConvertor result from q ---> Bsk U {m_tilde}
        Pointer tmp_encrypted_bsk_mtilde(allocate_poly(coeff_count * encrypted_count, bsk_mtilde_count, pool_));

        // Make temp poly for FastBConvertor result from Bsk U {m_tilde} -----> Bsk
        Pointer tmp_encrypted_bsk(allocate_poly(coeff_count * encrypted_count, bsk_base_mod_count_, pool_));

        // Step 0: fast base convert from q to Bsk U {m_tilde}
        // Step 1: reduce q-overflows in Bsk
        // Iterate over all the ciphertexts inside encrypted1
        for (int i = 0; i < encrypted_count; i++)
        {
            base_convertor_.fastbconv_mtilde(encrypted_safe.get() + (i * encrypted_ptr_increment), tmp_encrypted_bsk_mtilde.get() + (i * encrypted_bsk_mtilde_ptr_increment));
            base_convertor_.mont_rq(tmp_encrypted_bsk_mtilde.get() + (i * encrypted_bsk_mtilde_ptr_increment), tmp_encrypted_bsk.get() + (i * encrypted_bsk_ptr_increment));
        }

        // Step 2: compute product and multiply plain modulus to the result
        // We need to multiply both in q and Bsk. Values in encrypted_safe are in base q and values in tmp_encrypted_bsk are in base Bsk
        // We iterates over destination BigPolyArray and generates each BigPoly based on the indeces of inputs (arbitrary sizes for ciphertexts)
        // First allocate two temp polys: one for results in base q and the other for the result in base Bsk
        Pointer tmp_des_coeff_base(allocate_poly(coeff_count * dest_count, coeff_mod_count, pool_));
        Pointer tmp_des_bsk_base(allocate_poly(coeff_count * dest_count, bsk_base_mod_count_, pool_));

        // Reset the dest polys 
        set_zero_poly(coeff_count * dest_count, coeff_mod_count, tmp_des_coeff_base.get());
        set_zero_poly(coeff_count * dest_count, bsk_base_mod_count_, tmp_des_bsk_base.get());
        
        // First convert all the inputs into NTT form
        Pointer copy_encrypted_ntt_coeff_mod(allocate_poly(coeff_count * encrypted_count, coeff_mod_count, pool_));
        Pointer copy_encrypted_ntt_bsk_base_mod(allocate_poly(coeff_count * encrypted_count, bsk_base_mod_count_, pool_));

        set_poly_poly(encrypted_safe.get(), coeff_count * encrypted_count, coeff_mod_count, copy_encrypted_ntt_coeff_mod.get());
        set_poly_poly(tmp_encrypted_bsk.get(), coeff_count * encrypted_count, bsk_base_mod_count_, copy_encrypted_ntt_bsk_base_mod.get());

        for (int i = 0; i < encrypted_count; i++)
        {
            for (int j = 0; j < coeff_mod_count; j++)
            {
                smallntt_negacyclic_harvey_lazy(copy_encrypted_ntt_coeff_mod.get() + (j * coeff_count) + (i * encrypted_ptr_increment), coeff_small_ntt_tables_[j], pool_);
            }
            for (int j = 0; j < bsk_base_mod_count_; j++)
            {
                smallntt_negacyclic_harvey_lazy(copy_encrypted_ntt_bsk_base_mod.get() + (j * coeff_count) + (i * encrypted_bsk_ptr_increment), bsk_small_ntt_tables_[j], pool_);
            }
        }

        // Perform fast squaring
        // Compute c0^2 in base q
        for (int i = 0; i < coeff_mod_count; i++)
        {
            dyadic_product_coeff_smallmod(copy_encrypted_ntt_coeff_mod.get() + (i * coeff_count), copy_encrypted_ntt_coeff_mod.get() + (i * coeff_count), coeff_count, coeff_mod_array_[i], tmp_des_coeff_base.get() + (i * coeff_count), pool_); // Des[0] in q
            dyadic_product_coeff_smallmod(copy_encrypted_ntt_coeff_mod.get() + (i * coeff_count) + encrypted_ptr_increment, copy_encrypted_ntt_coeff_mod.get() + (i * coeff_count) + encrypted_ptr_increment, coeff_count, coeff_mod_array_[i], tmp_des_coeff_base.get() + (i * coeff_count) + (2 * encrypted_ptr_increment), pool_); // Des[2] in q
        }

        // Compute c0^2 in base bsk
        for (int i = 0; i < bsk_base_mod_count_; i++)
        {
            dyadic_product_coeff_smallmod(copy_encrypted_ntt_bsk_base_mod.get() + (i * coeff_count), copy_encrypted_ntt_bsk_base_mod.get() + (i * coeff_count), coeff_count, bsk_mod_array_[i], tmp_des_bsk_base.get() + (i * coeff_count), pool_); // Des[0] in bsk
            dyadic_product_coeff_smallmod(copy_encrypted_ntt_bsk_base_mod.get() + (i * coeff_count) + encrypted_bsk_ptr_increment, copy_encrypted_ntt_bsk_base_mod.get() + (i * coeff_count) + encrypted_bsk_ptr_increment, coeff_count, bsk_mod_array_[i], tmp_des_bsk_base.get() + (i * coeff_count) + (2 * encrypted_bsk_ptr_increment), pool_); // Des[2] in bsk
        }

        Pointer tmp_second_mul_coeff_base(allocate_poly(coeff_count, coeff_mod_count, pool_));
        // Compute 2*c0*c1 in base q
        for (int i = 0; i < coeff_mod_count; i++)
        {
            dyadic_product_coeff_smallmod(copy_encrypted_ntt_coeff_mod.get() + (i * coeff_count), copy_encrypted_ntt_coeff_mod.get() + (i * coeff_count) + encrypted_ptr_increment, coeff_count, coeff_mod_array_[i], tmp_second_mul_coeff_base.get() + (i * coeff_count), pool_);
            add_poly_poly_coeff_smallmod(tmp_second_mul_coeff_base.get() + (i * coeff_count), tmp_second_mul_coeff_base.get() + (i * coeff_count), coeff_count, coeff_mod_array_[i], tmp_des_coeff_base.get() + (i * coeff_count) + encrypted_ptr_increment); // D[1] in q
        }

        Pointer tmp_second_mul_bsk_base(allocate_poly(coeff_count, bsk_base_mod_count_, pool_));
        // Compute 2*c0*c1 in base bsk
        for (int i = 0; i < bsk_base_mod_count_; i++)
        {
            dyadic_product_coeff_smallmod(copy_encrypted_ntt_bsk_base_mod.get() + (i * coeff_count), copy_encrypted_ntt_bsk_base_mod.get() + (i * coeff_count) + encrypted_bsk_ptr_increment, coeff_count, bsk_mod_array_[i], tmp_second_mul_bsk_base.get() + (i * coeff_count), pool_);
            add_poly_poly_coeff_smallmod(tmp_second_mul_bsk_base.get() + (i * coeff_count), tmp_second_mul_bsk_base.get() + (i * coeff_count), coeff_count, bsk_mod_array_[i], tmp_des_bsk_base.get() + (i * coeff_count) + encrypted_bsk_ptr_increment); // D[1] in bsk
        }
        // Convert back outputs from NTT form
        for (int i = 0; i < dest_count; i++)
        {
            for (int j = 0; j < coeff_mod_count; j++)
            {
                inverse_smallntt_negacyclic_harvey(tmp_des_coeff_base.get() + (i * (encrypted_ptr_increment)) + (j * coeff_count), coeff_small_ntt_tables_[j], pool_);
            }
            for (int j = 0; j < bsk_base_mod_count_; j++)
            {
                inverse_smallntt_negacyclic_harvey(tmp_des_bsk_base.get() + (i * (encrypted_bsk_ptr_increment)) + (j * coeff_count), bsk_small_ntt_tables_[j], pool_);
            }
        }

        // Now we multiply plain modulus to both results in base q and Bsk and allocate them together in one container as (te0)q(te'0)Bsk | ... |te count)q (te' count)Bsk to make it ready for fast_floor 
        Pointer tmp_coeff_bsk_together(allocate_poly(coeff_count, dest_count * (coeff_mod_count + bsk_base_mod_count_), pool_));
        uint64_t *tmp_coeff_bsk_together_ptr = tmp_coeff_bsk_together.get();

        // Base q 
        for (int i = 0; i < dest_count; i++)
        {
            for (int j = 0; j < coeff_mod_count; j++)
            {
                multiply_poly_scalar_coeff_smallmod(tmp_des_coeff_base.get() + (j * coeff_count) + (i * encrypted_ptr_increment), coeff_count, parms_.plain_modulus().pointer(), coeff_mod_array_[j], tmp_coeff_bsk_together_ptr + (j * coeff_count));
            }
            tmp_coeff_bsk_together_ptr += encrypted_ptr_increment;

            for (int k = 0; k < bsk_base_mod_count_; k++)
            {
                multiply_poly_scalar_coeff_smallmod(tmp_des_bsk_base.get() + (k * coeff_count) + (i * encrypted_bsk_ptr_increment), coeff_count, parms_.plain_modulus().pointer(), bsk_mod_array_[k], tmp_coeff_bsk_together_ptr + (k * coeff_count));
            }
            tmp_coeff_bsk_together_ptr += encrypted_bsk_ptr_increment;
        }

        // Allocate a new poly for fast floor result in Bsk
        Pointer tmp_result_bsk(allocate_poly(coeff_count, dest_count * bsk_base_mod_count_, pool_));
        for (int i = 0; i < dest_count; i++)
        {
            // Step 3: fast floor from q U {Bsk} to Bsk 
            base_convertor_.fast_floor(tmp_coeff_bsk_together.get() + (i * (encrypted_ptr_increment + encrypted_bsk_ptr_increment)), tmp_result_bsk.get() + (i * encrypted_bsk_ptr_increment));

            // Step 4: fast base convert from Bsk to q
            base_convertor_.fastbconv_sk(tmp_result_bsk.get() + (i * encrypted_bsk_ptr_increment), destination_array.pointer(i));
        }
    }

    void RNSEvaluator::rns_compose(uint64_t *value)
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

        Pointer coefficients(allocate_zero_uint(total_uint64_count, pool_));
        uint64_t *coefficients_ptr = coefficients.get();

        // Re-merge the coefficients first
        for (int i = 0; i < coeff_count; i++)
        {
            for (int j = 0; j < coeff_mod_count; j++)
            {
                coefficients_ptr[(i * coeff_mod_count) + j] = value[(j * coeff_count) + i];
            }
        }

        Pointer temp(allocate_uint(coeff_mod_count, pool_));
        set_zero_uint(total_uint64_count, value);

        for (int i = 0; i < coeff_count; i++)
        {
            for (int j = 0; j < coeff_mod_count; j++)
            {
                uint64_t tmp;
                multiply_uint64_smallmod(coefficients_ptr + j, &inv_coeff_products_mod_coeff_array_[j], coeff_mod_array_[j], &tmp);
                multiply_uint_uint64(coeff_products_array_.get() + (j * coeff_mod_count), coeff_mod_count, tmp, coeff_mod_count, temp.get());
                add_uint_uint_mod(temp.get(), value + (i * coeff_mod_count), mod_.get(), coeff_mod_count, value + (i * coeff_mod_count));
            }
            set_zero_uint(coeff_mod_count, temp.get());
            coefficients_ptr += coeff_mod_count;
        }
    }

	//////////////////////////////////// New Functions ////////////////////////////////////////////////

	void RNSEvaluator::Polynomial_Evaluation(const Ciphertext &encrypted, vector<uint64_t> &coeff, Ciphertext &destination, const SmallModulus &eval_plain_modulus_)
	{
		// Save original paramters for recover at the end of this algorithm
		RNSEncryptionParameters origin_parms_(parms_);

		// Temperary change the plain_modulus to eval_plain_modulus
		Ciphertext encrypted_copy = encrypted;
		parms_.set_plain_modulus(eval_plain_modulus_);
		encrypted_copy.rns_hash_block_ = parms_.get_hash_block();

		int degree = coeff.size() - 1;

		// Compute parameter k = sqrt(degree/2), m = smallest interger with degree_prime = k * (2^m - 1) > degree
		//int k = (int)(sqrt(((double)degree) / 2.0) + 0.5);
        int k = optimal_parameter_paterson(degree);
		int degree_prime;
		int m = 1;
		while (1)
		{
			degree_prime = k * (power(2, m) - 1);
			if (degree_prime > degree) break;
			else m++;
		}
		//cout << "degree (old)  = " << degree << ", degree (new) = " << degree_prime << endl;
		//cout << "size of baby step = " << k + 1 << " size of giant step = " << m - 1 << endl;
		//cout << "number of non-scalar multiplications  = " << k - 1  + m - 1 + (1 << (m-1)) - 1 << endl;
		// If degree is low just compute all power and dot product
		if (degree < (k - 1) + 2 * (m - 1) + (1 << (m - 1)) - 1 || degree < 3) {
			//cout << "base case reached" << endl;
			vector<Ciphertext> power_of_all(degree + 1);
			compute_all_powers(encrypted_copy, degree, power_of_all);
			Polynomial_Evaluation(power_of_all, coeff, destination);
			//cout << "coeffs = "; 
			//for (int i = 0; i < coeff.size(); i++) {
			//	cout << coeff[i] << ", "; 
			//}
			//cout << endl;
			// Go back to original parms_
			parms_ = origin_parms_;
			destination.rns_hash_block_ = parms_.get_hash_block();
			return;
		}

		// Compute encryption of baby step (1, msg, msg^2, ... , msg^k)
		vector<Ciphertext> baby_step(k + 1);
		compute_all_powers(encrypted_copy, k, baby_step);

		// Compute encryption of giant step (msg^(2k), ... , msg^(2^{m-1}k))
		vector<Ciphertext> giant_step(m - 1);
		giant_step[0] = square(baby_step[k]);
		relinearize(giant_step[0], giant_step[0]);
		for (int i = 0; i < m - 2; i++)
		{
			square(giant_step[i], giant_step[i + 1]);
			relinearize(giant_step[i + 1], giant_step[i + 1]);
		}

		// Set f'(x) = X^degree_prime + f(x)
		vector<uint64_t> f_prime(degree_prime + 1);
		for (int i = 0; i < degree_prime + 1; i++)
		{
			if (i < coeff.size())
			{
				f_prime[i] = coeff[i];
			}
			else if (i < degree_prime)
			{
				f_prime[i] = 0;
			}
		}
		f_prime[degree_prime] = 1;

		//Evaluate f_prime
		Ciphertext encrypted_f_prime;
		PatersonStockmeyer(baby_step, giant_step, encrypted_copy, f_prime, encrypted_f_prime, k, m);

		//Compute encryption of X^degree_prime = X^k * X^2k ... X^{2^{m-1}k}
		Ciphertext power_of_degree_prime;
		multiply(baby_step[k], giant_step[0], power_of_degree_prime);
		relinearize(power_of_degree_prime, power_of_degree_prime, 2);
		for (int i = 0; i < m - 2; i++)
		{
			multiply(power_of_degree_prime, giant_step[i + 1], power_of_degree_prime);
			relinearize(power_of_degree_prime, power_of_degree_prime, 2);
		}

		//Substract X^degree_prime from temp
		sub(encrypted_f_prime, power_of_degree_prime, destination);

		// Go back to original parms_
		parms_ = origin_parms_;
		destination.rns_hash_block_ = parms_.get_hash_block();
	}

	void RNSEvaluator::DigitExtraction(Ciphertext &encrypted, vector<Ciphertext> &destination, long p, long e, long r, bool centered)
	{
		if (e < 2)
		{
			throw invalid_argument("paramter e should be larger than 1");
		}

		if (parms_.plain_modulus_base() != p || parms_.plain_modulus_exponent() != e)
		{
			throw invalid_argument("plain modulus is invalid for digit extraction");
		}

		//Extract encryption parameters.
		int coeff_count = parms_.poly_modulus().coeff_count();
		int coeff_mod_count = parms_.coeff_mod_array().size();
		int coeff_bit_count = coeff_mod_count * bits_per_uint64;
		int encrypted_count = encrypted.size();

		//Find k such that p^k > (e-1)(p-1)+1
		long k = 1;
		long pow_of_p = 1;
		while (1)
		{
			pow_of_p *= p;
			if (pow_of_p < (e - 1)*(p - 1) + 1) k++;
			else break;
		}


		//Initialize destination
		destination.resize(r);
		for (int i = 0; i < r; i++)
		{
			destination[i].get_mutable_array().resize(encrypted_count, coeff_count, coeff_bit_count);
			destination[i].hash_block_ = parms_.get_hash_block();
		}

		//lift[i] = poly_lift(3, i+1) which makes x + O(3^{i+1}) to x + O(3^{i+2})
		vector<vector<uint64_t> > lift(r - 1);
		for (int i = 0; i < r - 1; i++) {
			lift[i].resize(p + 1);
			lift[i] = compute_lift_uint64(p, i + 1);
		}

		//remainlsd[i](x) = [x]_3 in mod 3^{i+2}
		vector<vector<uint64_t> > remainlsd;
		remainlsd.resize(e - 1);
		for (int i = 0; i < e - 1; i++) {
			int degree = (i + 1)*(p - 1) + 1;
			remainlsd[i].resize(degree + 1);
			if (i >= (e - r - 1)) remainlsd[i] = compute_remainlsd_uint64(p, i + 2);
		}

		//Temp ciphertexts for computation
		vector<vector<Ciphertext> > lift_encrypted;
		lift_encrypted.resize(r);
		for (int i = 0; i < r; i++)
		{
			lift_encrypted[i].resize(k - 1);
		}

		Ciphertext cipher;
		cipher.get_mutable_array().resize(encrypted_count, coeff_count, coeff_bit_count);
		cipher = encrypted;

		//BigUInt plain_modulus_ = parms_.plain_modulus();
		uint64_t plain_modulus_ = parms_.plain_modulus().value();
		uint64_t halfp = p / 2;
		for (int i = 0; i < r; i++)
		{
			// p**i * p/2;
			uint64_t shift = power(p, i) *halfp;
			Plaintext const_shift(BigPoly(1, parms_.plain_modulus().bit_count(), &shift));

			if (centered) add_plain(cipher, const_shift, cipher);
			SmallModulus eval_plain_modulus(plain_modulus_);
			if (i < r - 1)
			{
				Polynomial_Evaluation(cipher, lift[0], lift_encrypted[i][0], eval_plain_modulus);
				if (centered) sub_plain(lift_encrypted[i][0], const_shift, lift_encrypted[i][0]);

			}
			Polynomial_Evaluation(cipher, remainlsd[e - i - 2], destination[i], eval_plain_modulus);
			if (centered) sub_plain(destination[i], const_shift, destination[i]);
			for (int j = 0; j < k - 2; j++)
			{
				if (i + j < r - 1)
				{
					if (centered) add_plain(lift_encrypted[i][j], const_shift, lift_encrypted[i][j]);
					Polynomial_Evaluation(lift_encrypted[i][j], lift[j + 1], lift_encrypted[i][j + 1], eval_plain_modulus);
					if (centered) sub_plain(lift_encrypted[i][j+1], const_shift, lift_encrypted[i][j+1]);
				}
			}
			cipher = encrypted;
			if (i < r - 1)
			{
				for (int j = 0; j < i + 1; j++)
				{
					if (i - j > k - 2)
					{
						sub(cipher, destination[j], cipher);
					}
					else
					{
						sub(cipher, lift_encrypted[j][i - j], cipher);
					}
				}
			}
			plain_modulus_ /= p;
		}
	}

	void RNSEvaluator::apply_galois(const Ciphertext &encrypted, uint64_t galois_elt, Ciphertext &destination, bool ntt_output)
	{
		const BigPolyArray &encrypted_array = encrypted.get_array();
		BigPolyArray &destination_array = destination.get_mutable_array();

		// Extract paramters
		int coeff_count = parms_.poly_modulus().coeff_count();
		int coeff_mod_count = parms_.coeff_mod_array().size();
		int coeff_bit_count = coeff_mod_count * bits_per_uint64;
		int encrypted_count = encrypted.size();
		int encrypted_ptr_increment = coeff_count * coeff_mod_count;

		// Verify coprime conditions. 
		if (galois_elt % 2 == 0 || galois_elt >= 2 * coeff_count)
		{
			throw invalid_argument("galois element is not valid");
		}

		// Verify parameters.
		if (encrypted.rns_hash_block_ != parms_.get_hash_block())
		{
			throw invalid_argument("encrypted is not valid for encryption parameters");
		}

		if (encrypted_count > 2)
		{
			throw invalid_argument("ciphertext size bigger than 2 is not yet implemented.");
		}

		// Check whethere galois key is generated or not.
		if (galois_keys_.has_key(galois_elt) == false)
		{
			int n = coeff_count - 1;
			int m = n << 1;

			// galois_etl = 3^order1 * (-1)^order2
			int order1 = map_from_Zmstar_to_generator.at(galois_elt).first;
			int order2 = map_from_Zmstar_to_generator.at(galois_elt).second;

			// Split order1 to bit
			vector<int> bit_order1;
			while (1)
			{
				if (order1 == 0)
				{
					break;
				}
				else
				{
					bit_order1.push_back(order1 % 2);
					order1 >>= 1;
				}
			}

			destination = encrypted;
			uint64_t two_power_of_three = 3;
			for (int i = 0; i < bit_order1.size(); i++)
			{
				if (bit_order1[i] != 0)
				{
					if (galois_keys_.has_key(two_power_of_three) == false)
					{
						throw invalid_argument("Turn on the auto_galois_key_gen option in generate");
					}
					Ciphertext temp;
					apply_galois(destination, two_power_of_three, temp);
					destination = temp;
				}
				two_power_of_three *= two_power_of_three;
				two_power_of_three %= m;
			}
			if (order2 == 1)
			{
				if (galois_keys_.has_key(m - 1) == false)
				{
					throw invalid_argument("Turn on the auto_galois_key_gen option in generate");
				}
				Ciphertext temp;
				apply_galois(destination, m - 1, temp);
				destination = temp;
			}
			destination.rns_hash_block_ = parms_.get_hash_block();
			return;
		}

		ConstPointer encrypted_safe = duplicate_bigpolyarray_if_needed(encrypted_array, &encrypted == &destination, pool_);

		// Prepare destination
		destination_array.resize(encrypted_count, coeff_count, coeff_bit_count);
		destination.rns_hash_block_ = parms_.get_hash_block();

		// Apply Galois for each ciphertext
		Pointer temp0(allocate_zero_uint(coeff_count * coeff_mod_count, pool_));
		Pointer temp1(allocate_zero_uint(coeff_count * coeff_mod_count, pool_));

		for (int i = 0; i < coeff_mod_count; i++)
		{
			apply_galois_poly_smallmod(encrypted_safe.get() + (i * coeff_count), coeff_count - 1, galois_elt, coeff_mod_array_[i], temp0.get() + i * coeff_count);
			apply_galois_poly_smallmod(encrypted_safe.get() + (i * coeff_count) + encrypted_ptr_increment, coeff_count - 1, galois_elt, coeff_mod_array_[i], temp1.get() + i * coeff_count);
		}

		// Calculate (temp1 * galois_key.first, temp1 * galois_key.second) + (temp0, temp1)
		const uint64_t *encrypted_coeff = temp1.get();
		Pointer encrypted_coeff_prod_inv_coeff(allocate_uint(coeff_count * coeff_mod_count, pool_));

		// decompose encrypted_array[count-1] into base w
		// want to create an array of polys, each of whose components i is (encrypted_array[count-1])^(i) - in the notation of fv paper
		Pointer decomp_encrypted_last(allocate_zero_poly(coeff_count, evaluation_keys_.get_data()[0][0].first.size() * coeff_mod_count, pool_));

		Pointer innerproduct0(allocate_poly(coeff_count, coeff_mod_count, pool_));
		Pointer innerproduct1(allocate_poly(coeff_count, coeff_mod_count, pool_));
		Pointer innerresult0(allocate_zero_poly(coeff_count, coeff_mod_count, pool_));
		Pointer innerresult1(allocate_zero_poly(coeff_count, coeff_mod_count, pool_));
		Pointer temp_decomp_coeff(allocate_poly(coeff_count, 1, pool_));

		int array_poly_uint64_count = coeff_count * coeff_mod_count;

		if (qualifiers_.enable_ntt)
		{
			for (int i = 0; i < coeff_mod_count; i++)
			{
				multiply_poly_scalar_coeff_smallmod(encrypted_coeff + (i * coeff_count), coeff_count, &inv_coeff_products_mod_coeff_array_[i], coeff_mod_array_[i], encrypted_coeff_prod_inv_coeff.get() + (i * coeff_count));

				// populate one poly at a time.
				// create a polynomial to store the current decomposition value which will be copied into the array to populate it at the current index.
				int shift = 0;
				for (int k = 0; k < evaluation_keys_.get_data()[0][0].first.size(); k++)
				{
					// Decompose here
					uint64_t *decomp_coeff = decomp_encrypted_last.get() + k * array_poly_uint64_count;
					for (int coeff_index = 0; coeff_index < coeff_count; ++coeff_index)
					{
						right_shift_uint(encrypted_coeff_prod_inv_coeff.get() + coeff_index + (i * coeff_count), shift, 1, decomp_coeff + coeff_index + (i * coeff_count));
						filter_highbits_uint(decomp_coeff + coeff_index + (i * coeff_count), 1, parms_.decomposition_bit_count());
					}

					for (int j = 0; j < coeff_mod_count; j++)
					{
						set_poly_poly(decomp_coeff + (i * coeff_count), coeff_count, 1, temp_decomp_coeff.get());
						smallntt_negacyclic_harvey(temp_decomp_coeff.get(), coeff_small_ntt_tables_[j], pool_);
						dyadic_product_coeff_smallmod(temp_decomp_coeff.get(), galois_keys_.get_mutable_data()[galois_elt][i].first.pointer(k) + (j * coeff_count), coeff_count, coeff_mod_array_[j], innerproduct0.get() + (j * coeff_count), pool_);
						dyadic_product_coeff_smallmod(temp_decomp_coeff.get(), galois_keys_.get_mutable_data()[galois_elt][i].second.pointer(k) + (j * coeff_count), coeff_count, coeff_mod_array_[j], innerproduct1.get() + (j * coeff_count), pool_);
						add_poly_poly_coeff_smallmod(innerresult0.get() + (j * coeff_count), innerproduct0.get() + (j * coeff_count), coeff_count, coeff_mod_array_[j], innerresult0.get() + (j * coeff_count));
						add_poly_poly_coeff_smallmod(innerresult1.get() + (j * coeff_count), innerproduct1.get() + (j * coeff_count), coeff_count, coeff_mod_array_[j], innerresult1.get() + (j * coeff_count));
					}
					shift += parms_.decomposition_bit_count();
				}
			}

			for (int i = 0; i < coeff_mod_count; i++)
			{
				if (ntt_output)
				{
					destination.in_ntt_form = true;
					smallntt_negacyclic_harvey(temp0.get() + (i * coeff_count), coeff_small_ntt_tables_[i], pool_);
					add_poly_poly_coeff_smallmod(temp0.get() + (i * coeff_count), innerresult0.get() + (i * coeff_count), coeff_count, coeff_mod_array_[i], destination.get_mutable_array().pointer(0) + (i * coeff_count));
					set_poly_poly(innerresult1.get() + (i * coeff_count), coeff_count, 1, destination.get_mutable_array().pointer(1) + (i * coeff_count));
				}
				else
				{
					inverse_smallntt_negacyclic_harvey(innerresult0.get() + (i * coeff_count), coeff_small_ntt_tables_[i], pool_);
					inverse_smallntt_negacyclic_harvey(innerresult1.get() + (i * coeff_count), coeff_small_ntt_tables_[i], pool_);
					add_poly_poly_coeff_smallmod(temp0.get() + (i * coeff_count), innerresult0.get() + (i * coeff_count), coeff_count, coeff_mod_array_[i], destination.get_mutable_array().pointer(0) + (i * coeff_count));
					set_poly_poly(innerresult1.get() + (i * coeff_count), coeff_count, 1, destination.get_mutable_array().pointer(1) + (i * coeff_count));
				}
			}
		}
		else
		{
			// this branch should never be reached
			throw logic_error("invalid encryption parameters");
		}
	}

	void RNSEvaluator::apply_galois_many(const Ciphertext &encrypted, vector<uint64_t> galois_elts, vector<Ciphertext> &destination, bool ntt_output)
	{
		const BigPolyArray &encrypted_array = encrypted.get_array();
		
		// Extract paramters
		int coeff_count = parms_.poly_modulus().coeff_count();
		int coeff_mod_count = parms_.coeff_mod_array().size();
		int coeff_bit_count = coeff_mod_count * bits_per_uint64;
		int encrypted_count = encrypted.size();
		int encrypted_ptr_increment = coeff_count * coeff_mod_count;
		int decompose_number = galois_keys_.get_mutable_data()[galois_elts[0]][0].first.size();
		int galois_elts_size = galois_elts.size();
	
		for (int i = 0; i < galois_elts_size; i++)
		{
			// Verify coprime conditions.
			if (galois_elts[i] % 2 == 0 || galois_elts[i] >= 2 * coeff_count)
			{
				throw invalid_argument("galois element is not valid");
			}
			// Check whethere galois key is generated or not.
			if (galois_keys_.has_key(galois_elts[i]) == false)
			{
				throw invalid_argument("galois key is not generated");
			}
		}

		// Verify parameters.
		if (encrypted.rns_hash_block_ != parms_.get_hash_block())
		{
			throw invalid_argument("encrypted is not valid for encryption parameters");
		}

		if (encrypted_count > 2)
		{
			throw invalid_argument("ciphertext size bigger than 2 is not yet implemented.");
		}

		Pointer encrypted_safe(allocate_poly(encrypted_count * coeff_count, coeff_mod_count, pool_));
		set_poly_poly(encrypted.get_array().pointer(0), encrypted_count * coeff_count, coeff_mod_count, encrypted_safe.get());

		// Compute permute index for galois elts
		map<int, vector<int> > permute_index;
		int logn = get_power_of_two(coeff_count - 1);
		for (int l = 0; l < galois_elts.size(); l++)
		{
			vector<int> index(coeff_count - 1);
			for (int i = 0; i < coeff_count - 1; i++)
			{
				int reversed = bit_reversal_general(i, logn);
				int index_raw = (2 * reversed + 1) * galois_elts[l];
				index_raw %= (2 * (coeff_count - 1));
				index[i] = bit_reversal_general((index_raw - 1) >> 1, logn);
			}
			permute_index.emplace(galois_elts[l], index);
		}

		// Prepare destination
		destination.resize(galois_elts_size);
		for (int i = 0; i < galois_elts_size; i++)
		{
			destination[i].get_mutable_array().resize(encrypted_count, coeff_count, coeff_bit_count);
			destination[i].rns_hash_block_ = parms_.get_hash_block();
		}

		const uint64_t *encrypted_coeff = encrypted_safe.get() + encrypted_ptr_increment;
		Pointer encrypted_coeff_prod_inv_coeff(allocate_uint(coeff_count * coeff_mod_count, pool_));

		// decompose encrypted_array[count-1] into base w
		// want to create an array of polys, each of whose components i is (encrypted_array[count-1])^(i) - in the notation of fv paper
		Pointer decomp_encrypted_last(allocate_zero_poly(coeff_count, evaluation_keys_.get_data()[0][0].first.size() * coeff_mod_count, pool_));

		Pointer innerproduct0(allocate_poly(galois_elts_size * coeff_count, coeff_mod_count, pool_));
		Pointer innerproduct1(allocate_poly(galois_elts_size * coeff_count, coeff_mod_count, pool_));
		Pointer innerresult0(allocate_zero_poly(galois_elts_size * coeff_count, coeff_mod_count, pool_));
		Pointer innerresult1(allocate_zero_poly(galois_elts_size * coeff_count, coeff_mod_count, pool_));
		Pointer temp_decomp_coeff(allocate_uint(coeff_count, pool_));
		Pointer temp_permuted_coeff(allocate_zero_poly(coeff_count, 1, pool_));
		
		int array_poly_uint64_count = coeff_count * coeff_mod_count;
		int galois_elts_size_time_coeff_count = galois_elts_size * coeff_count;

		if (qualifiers_.enable_ntt)
		{
			for (int i = 0; i < coeff_mod_count; i++)
			{
				multiply_poly_scalar_coeff_smallmod(encrypted_coeff + (i * coeff_count), coeff_count, &inv_coeff_products_mod_coeff_array_[i], coeff_mod_array_[i], encrypted_coeff_prod_inv_coeff.get() + (i * coeff_count));

				// populate one poly at a time.
				// create a polynomial to store the current decomposition value which will be copied into the array to populate it at the current index.
				int shift = 0;
				for (int k = 0; k < decompose_number; k++)
				{
					// Decompose here
					uint64_t *decomp_coeff = decomp_encrypted_last.get() + k * array_poly_uint64_count;
					for (int coeff_index = 0; coeff_index < coeff_count; ++coeff_index)
					{
						right_shift_uint(encrypted_coeff_prod_inv_coeff.get() + coeff_index + (i * coeff_count), shift, 1, decomp_coeff + coeff_index + (i * coeff_count));
						filter_highbits_uint(decomp_coeff + coeff_index + (i * coeff_count), 1, parms_.decomposition_bit_count());
					}

					for (int j = 0; j < coeff_mod_count; j++)
					{
						set_poly_poly(decomp_coeff + (i * coeff_count), coeff_count, 1, temp_decomp_coeff.get());
						smallntt_negacyclic_harvey(temp_decomp_coeff.get(), coeff_small_ntt_tables_[j], pool_);
						for (int l = 0; l < galois_elts_size; l++)
						{
							permute_ntt_poly_smallmod_with_index(temp_decomp_coeff.get(), coeff_count - 1, permute_index.at(galois_elts[l]), temp_permuted_coeff.get());
							dyadic_product_coeff_smallmod(temp_permuted_coeff.get(), galois_keys_.get_mutable_data()[galois_elts[l]][i].first.pointer(k) + (j * coeff_count), coeff_count, coeff_mod_array_[j], innerproduct0.get() + (j * coeff_count) + (l * array_poly_uint64_count), pool_);
							dyadic_product_coeff_smallmod(temp_permuted_coeff.get(), galois_keys_.get_mutable_data()[galois_elts[l]][i].second.pointer(k) + (j * coeff_count), coeff_count, coeff_mod_array_[j], innerproduct1.get() + (j * coeff_count) + (l * array_poly_uint64_count), pool_);
							add_poly_poly_coeff_smallmod(innerresult0.get() + (j * coeff_count) + (l * array_poly_uint64_count), innerproduct0.get() + (j * coeff_count) + (l * array_poly_uint64_count), coeff_count, coeff_mod_array_[j], innerresult0.get() + (j * coeff_count) + (l * array_poly_uint64_count));
							add_poly_poly_coeff_smallmod(innerresult1.get() + (j * coeff_count) + (l * array_poly_uint64_count), innerproduct1.get() + (j * coeff_count) + (l * array_poly_uint64_count), coeff_count, coeff_mod_array_[j], innerresult1.get() + (j * coeff_count) + (l * array_poly_uint64_count));
						}
					}
					shift += parms_.decomposition_bit_count();
				}
			}
			
			// rotate encrypted[0]
			Pointer rotated_encrypted0(allocate_zero_uint(coeff_count, pool_));
			for (int l = 0; l < galois_elts_size; l++)
			{
				for (int i = 0; i < coeff_mod_count; i++)
				{
					if (ntt_output)
					{
						destination[l].in_ntt_form = true;
						apply_galois_poly_smallmod(encrypted_safe.get() + (i * coeff_count), coeff_count - 1, galois_elts[l], coeff_mod_array_[i], rotated_encrypted0.get());
						smallntt_negacyclic_harvey(rotated_encrypted0.get(), coeff_small_ntt_tables_[i], pool_);
						add_poly_poly_coeff_smallmod(rotated_encrypted0.get(), innerresult0.get() + (i * coeff_count) + (l * array_poly_uint64_count), coeff_count, coeff_mod_array_[i], destination[l].get_mutable_array().pointer(0) + (i * coeff_count));
						set_poly_poly(innerresult1.get() + (i * coeff_count) + (l * array_poly_uint64_count), coeff_count, 1, destination[l].get_mutable_array().pointer(1) + (i * coeff_count));
					}
					else
					{
						apply_galois_poly_smallmod(encrypted_safe.get() + (i * coeff_count), coeff_count - 1, galois_elts[l], coeff_mod_array_[i], rotated_encrypted0.get());
						inverse_smallntt_negacyclic_harvey(innerresult0.get() + (i * coeff_count) + (l * array_poly_uint64_count), coeff_small_ntt_tables_[i], pool_);
						inverse_smallntt_negacyclic_harvey(innerresult1.get() + (i * coeff_count) + (l * array_poly_uint64_count), coeff_small_ntt_tables_[i], pool_);
						add_poly_poly_coeff_smallmod(rotated_encrypted0.get(), innerresult0.get() + (i * coeff_count) + (l * array_poly_uint64_count), coeff_count, coeff_mod_array_[i], destination[l].get_mutable_array().pointer(0) + (i * coeff_count));
						set_poly_poly(innerresult1.get() + (i * coeff_count) + (l * array_poly_uint64_count), coeff_count, 1, destination[l].get_mutable_array().pointer(1) + (i * coeff_count));
					}
				}

			}
		}
		else
		{
			// this branch should never be reached
			throw logic_error("invalid encryption parameters");
		}
	}

	void RNSEvaluator::matrix_vector_product_plain(std::vector<Ciphertext>& encrypted, const std::vector<std::vector<Plaintext>>&coeff_matrix, std::vector<Ciphertext>& destination, bool ntt_output)
	{
		BigPolyArray encarray = encrypted[0].get_array();
		int count = encarray.size();
		// int coeff_bit_count = encarray.coeff_bit_count();
		int coeff_count = parms_.poly_modulus().coeff_count();
		int coeff_mod_count = parms_.coeff_mod_array().size();
		int coeff_bit_count = coeff_mod_count * bits_per_uint64;
		int encrypted_count = encrypted.size();
		int encrypted_ptr_increment = coeff_count * coeff_mod_count;
		// assert encrypted.size() == coeff_matrix[0].size();
		destination.resize(coeff_matrix.size());
		//for (int i = 0; i < destination.size(); i++) {
		//	BigPolyArray &destination_array = destination[i].get_mutable_array();
		//	destination_array.resize(encrypted_count, coeff_count, coeff_bit_count);
		//	destination[i].rns_hash_block_ = parms_.get_hash_block();
		//}
		int din = encrypted.size();
		int dout = destination.size();

		// cout << "performing a matrix vector product with input dimension =  "
		//    << din << " and output dimension = " << dout << endl;

		bool use_ntt = qualifiers_.enable_ntt;
		//bool use_ntt = false;
		if (use_ntt) {
			// cout << "   going through the ntt branch..." << endl;


			// First we transform all the plain coefficients to ntt. This includes transforming from t to q size. 
			// This needs to be .... independently timed. 
			//vector<vector<Plaintext>> coeff_matrix_ntt(dout);
			//for (int j = 0; j < dout; j++) {
			//	coeff_matrix_ntt[j].resize(din);
			//	for (int i = 0; i < din; i++) {
			//		coeff_matrix_ntt[j][i] = transform_to_ntt(coeff_matrix[j][i]);
			//	}
			//}


			auto time_Mat_start = chrono::high_resolution_clock::now();

			// Next transform the input to NTT domain if necessary. 
			// vector<Ciphertext> encrypted_ntt(encrypted.size());
			for (int i = 0; i < din; i++) {
				//encrypted_ntt[i] = encrypted[i];
				// transform to NTT only when necessary
				if (!encrypted[i].is_ntt_form()) {
					transform_to_ntt(encrypted[i]);
				}
				//else {
				//	if (i == 0) cout << " skipped an NTT in matrix vector product ..." << endl;
				//}
			}

			// Next, we perform dot_product and sum over ntt domain. 
			for (int j = 0; j < dout; j++) {
				destination[j].get_mutable_array().resize(count, coeff_count, coeff_bit_count);
				destination[j].get_mutable_array().set_zero();
				destination[j].rns_hash_block_ = parms_.get_hash_block();
				Ciphertext temp;
				for (int i = 0; i < din; i++) {
					if (coeff_matrix[j][i].is_zero() == false)
					{
						Plaintext temp_plain;
						transform_to_ntt(coeff_matrix[j][i], temp_plain);
						multiply_plain_ntt(encrypted[i], temp_plain, temp);
						add(destination[j], temp, destination[j]);
					}
				}
				if (!ntt_output) {
					transform_from_ntt(destination[j]);
				}
				destination[j].set_ntt_form(ntt_output);
			}
			auto time_Mat_end = chrono::high_resolution_clock::now();
			chrono::microseconds time_pre = chrono::duration_cast<chrono::microseconds>(time_Mat_end - time_Mat_start);
			cout << "+matrix vector product online time : " << (double)time_pre.count() / 1000000 << " seconds" << endl;

		}
		else {
			cout << "matrix vector product going through non-NTT branch ..." << endl;
			for (int j = 0; j < dout; j++) {
				//destination[j].get_mutable_array().resize(count, coeff_count, coeff_bit_count);
				//destination[j].get_mutable_array().set_zero();
				//destination[j].hash_block_ = parms_.get_hash_block();
				Ciphertext temp;
				for (int i = 0; i < din; i++) {
					multiply_plain(encrypted[i], coeff_matrix[j][i], temp);
					add(destination[j], temp, destination[j]);
				}
			}
		}
		return;
	}

	void RNSEvaluator::matrix_vector_product_plain_ntt(std::vector<Ciphertext>& encrypted, const std::vector<std::vector<Plaintext>>&coeff_matrix_ntt, std::vector<Ciphertext>& destination, bool ntt_output)
	{
		BigPolyArray encarray = encrypted[0].get_array();
		int count = encarray.size();
		// int coeff_bit_count = encarray.coeff_bit_count();
		int coeff_count = parms_.poly_modulus().coeff_count();
		int coeff_mod_count = parms_.coeff_mod_array().size();
		int coeff_bit_count = coeff_mod_count * bits_per_uint64;
		int encrypted_count = encrypted.size();
		int encrypted_ptr_increment = coeff_count * coeff_mod_count;
		// assert encrypted.size() == coeff_matrix[0].size();
		destination.resize(coeff_matrix_ntt.size());

		int din = encrypted.size();
		int dout = destination.size();

		if (qualifiers_.enable_ntt) {
			auto time_Mat_start = chrono::high_resolution_clock::now();
			for (int i = 0; i < din; i++)
			{
				if (!encrypted[i].is_ntt_form())
				{
					transform_to_ntt(encrypted[i]);
				}
			}
			for (int j = 0; j < dout; j++)
			{
				destination[j].get_mutable_array().resize(count, coeff_count, coeff_bit_count);
				destination[j].get_mutable_array().set_zero();
				destination[j].rns_hash_block_ = parms_.get_hash_block();
				Ciphertext temp;
				for (int i = 0; i < din; i++) {
					if (coeff_matrix_ntt[j][i].is_zero() == false)
					{
						multiply_plain_ntt(encrypted[i], coeff_matrix_ntt[j][i], temp);
						add(destination[j], temp, destination[j]);
					}
				}
				if (!ntt_output) {
					transform_from_ntt(destination[j]);
				}
				destination[j].set_ntt_form(ntt_output);
			}
			auto time_Mat_end = chrono::high_resolution_clock::now();
			chrono::microseconds time_pre = chrono::duration_cast<chrono::microseconds>(time_Mat_end - time_Mat_start);
			cout << "+matrix vector product online time : " << (double)time_pre.count() / 1000000 << " seconds" << endl;
		}
		else
		{
			throw;
		}
		return;
	}

	void RNSEvaluator::apply_galois_plain(const Plaintext &poly, uint64_t galois_elt, Plaintext &destination) {
		// TODO
		int coeff_count = parms_.poly_modulus().coeff_count();
		int coeff_mod_count = parms_.coeff_mod_array().size();
		int coeff_bit_count = coeff_mod_count * bits_per_uint64;
		//int coeff_uint64_count = parms_.coeff_modulus().uint64_count();
		int plain_coeff_bit_count = parms_.plain_modulus().bit_count();
		int plain_coeff_uint64_count = parms_.plain_modulus().uint64_count();
		if (galois_elt % 2 == 0 || galois_elt >= 2 * coeff_count)
		{
			throw invalid_argument("galois element is not valid");
		}

		// Verify parameters.
		if (poly.get_poly().coeff_count() != coeff_count || poly.get_poly().coeff_bit_count() != plain_coeff_bit_count)
		{
			throw invalid_argument("parameters");
		}
		// First preencrypt, then .... 
		destination.get_poly().resize(coeff_count, plain_coeff_bit_count);
		//Pointer widepoly(allocate_poly(coeff_count, coeff_uint64_count, pool_));
		// Fixme: transfer the poly to wide poly
		//preencrypt(poly.pointer(), coeff_count, plain_coeff_uint64_count, widepoly.get());
		apply_galois_poly_smallmod(poly.get_poly().pointer(), coeff_count - 1, galois_elt, parms_.plain_modulus(), destination.get_poly().pointer());
		set_zero_uint(plain_coeff_uint64_count, get_poly_coeff(destination.get_poly().pointer(), coeff_count - 1, plain_coeff_uint64_count));
		return;
	}

	// Inplace apply Galois.
	void RNSEvaluator::apply_galois_plain(Plaintext & poly, uint64_t galois_elt)
	{	
		int coeff_count = parms_.poly_modulus().coeff_count();
		int plain_coeff_bit_count = parms_.plain_modulus().bit_count();
		int plain_coeff_uint64_count = parms_.plain_modulus().uint64_count();
		if (galois_elt % 2 == 0 || galois_elt >= 2 * coeff_count)
		{
			throw invalid_argument("galois element is not valid");
		}

		// Verify parameters.
		if (poly.get_poly().coeff_count() != coeff_count || poly.get_poly().coeff_bit_count() != plain_coeff_bit_count)
		{
			throw invalid_argument("parameters");
		}
		// First preencrypt, then .... 
		//.get_poly().resize(coeff_count, plain_coeff_bit_count);
		Pointer temp(allocate_poly(coeff_count, plain_coeff_uint64_count, pool_));
		// Fixme: transfer the poly to wide poly
		//preencrypt(poly.pointer(), coeff_count, plain_coeff_uint64_count, widepoly.get());
		apply_galois_poly_smallmod(poly.get_poly().pointer(), coeff_count - 1, galois_elt, parms_.plain_modulus(), temp.get());
		set_zero_uint(plain_coeff_uint64_count, get_poly_coeff(temp.get(), coeff_count - 1, plain_coeff_uint64_count));
		set_poly_poly(temp.get(), coeff_count, plain_coeff_uint64_count, poly.get_poly().pointer());
	}

	void RNSEvaluator::rotations_weighted_sum(const Ciphertext & encrypted, std::vector<std::vector<Plaintext>>& coeffs, std::vector<std::uint64_t> babystep, std::vector<std::uint64_t> giantstep, Ciphertext & destination, bool ntt_output)
	{
		int coeff_count = parms_.poly_modulus().coeff_count();
		int coeff_mod_count = parms_.coeff_mod_array().size();
		int coeff_bit_count = coeff_mod_count * bits_per_uint64;
		int encrypted_count = encrypted.size();

		int n = coeff_count - 1;
		int m = 2 * n;
		BigPolyArray temp, encrypted_bs, temp2;

		destination.get_mutable_array().resize(encrypted_count, coeff_count, coeff_bit_count);
		destination.get_mutable_array().set_zero();
		destination.rns_hash_block_ = parms_.get_hash_block();

		cout << " evaluating rotations_weighted_sum using babystep = " << babystep.size() << " and giant step = " << giantstep.size() << endl;

		// TODO: rewrite this as a matrix-vector product followed by rotation-sum.



		// chrono::microseconds rotations_time(0);
		int gs, bs;

		// TODO: we can reduce some memory usage here? 
		// vector<vector<Plaintext>> plain_coeffs(giantstep.size());
		for (int j = 0; j < giantstep.size(); j++) {
			gs = giantstep[j];
			//  plain_coeffs[j].resize(babystep.size());
			for (int i = 0; i < babystep.size(); i++) {
				bs = babystep[i];
				//apply_galois(encrypted, bs, encrypted_bs);
				// int coeff_index = (bs*gs % m - 1) / 2;
				apply_galois_plain(coeffs[j][i], mod_inverse(gs, m));
			}
		}

		// Compute all rotations in babystep.
		vector<Ciphertext> babystep_rotations(babystep.size());
		cout << "computing babystep rotations using hoisting ... " << endl;
		auto time_start = chrono::high_resolution_clock::now();
		apply_galois_many(encrypted, babystep, babystep_rotations, true);
		auto time_end = chrono::high_resolution_clock::now();
		chrono::microseconds online_time = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
		cout << "+rotations weighted sum hoisting time: " << (double)online_time.count() / 1000000 << " seconds" << endl;

		//for (int i = 0; i < babystep.size();  i++) {
		//    babystep_rotations[i].get_array().resize(encrypted_count, coeff_count, coeff_bit_count);
		//    babystep_rotations[i].get_array().set_zero();
		//    babystep_rotations[i].hash_block_ = parms_.get_hash_block();
		//    bs = babystep[i]; 
		//    apply_galois(encrypted, bs, babystep_rotations[i]); 
		//}

		vector<Ciphertext> intermediate_vector(giantstep.size());;

		// Populate plain coefficients; in this case, the output has dimension giantstep
		// and the input has dimension babystep. So this is correct. 

		matrix_vector_product_plain(babystep_rotations, coeffs, intermediate_vector, false);


		auto time_gs_start = chrono::high_resolution_clock::now();
		for (int j = 0; j < giantstep.size(); j++) {
			Ciphertext temp;
			apply_galois(intermediate_vector[j], giantstep[j], temp, true);
			//add(destination, intermediate_vector[j], destination);
			add(destination, temp, destination);
		}
		if (!ntt_output) {
			transform_from_ntt(destination);
		}
		auto time_gs_end = chrono::high_resolution_clock::now();
		chrono::microseconds time_pre = chrono::duration_cast<chrono::microseconds>(time_gs_end - time_gs_start);
		cout << "+giant step time : " << (double)time_pre.count() / 1000000 << " seconds" << endl;
		return;
	}

	void RNSEvaluator::rotations_weighted_sum(const Ciphertext & encrypted, const vector<Plaintext>& coeffs, vector<uint64_t> babystep, vector<uint64_t> giantstep, Ciphertext & destination, bool ntt_output)
	{

		int coeff_count = parms_.poly_modulus().coeff_count();
		int coeff_mod_count = parms_.coeff_mod_array().size();
		int coeff_bit_count = coeff_mod_count * bits_per_uint64;
		int encrypted_count = encrypted.size();

		int n = coeff_count - 1;
		int m = 2 * n;
		BigPolyArray temp, encrypted_bs, temp2;

		destination.get_mutable_array().resize(encrypted_count, coeff_count, coeff_bit_count);
		destination.get_mutable_array().set_zero();
		destination.rns_hash_block_ = parms_.get_hash_block();

		cout << " evaluating rotations_weighted_sum using babystep = " << babystep.size() << " and giant step = " << giantstep.size() << endl;

		int gs, bs;

		// TODO: we can reduce some memory usage here? 
		vector<vector<Plaintext>> plain_coeffs(giantstep.size());
		for (int j = 0; j < giantstep.size(); j++) {
			gs = giantstep[j];
			plain_coeffs[j].resize(babystep.size());
			for (int i = 0; i < babystep.size(); i++) {
				bs = babystep[i];
				//apply_galois(encrypted, bs, encrypted_bs);
				int coeff_index = (bs*gs % m - 1) / 2;
				apply_galois_plain(coeffs[coeff_index], mod_inverse(gs, m), plain_coeffs[j][i]);
			}
		}

		// Compute all rotations in babystep.
		vector<Ciphertext> babystep_rotations(babystep.size());
		cout << "computing babystep rotations using hoisting ... " << endl;
		auto time_start = chrono::high_resolution_clock::now();
		apply_galois_many(encrypted, babystep, babystep_rotations, true);
		auto time_end = chrono::high_resolution_clock::now();
		chrono::microseconds online_time = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
		cout << "+rotations weighted sum hoisting time: " << (double)online_time.count() / 1000000 << " seconds" << endl;

		//for (int i = 0; i < babystep.size();  i++) {
		//    babystep_rotations[i].get_array().resize(encrypted_count, coeff_count, coeff_bit_count);
		//    babystep_rotations[i].get_array().set_zero();
		//    babystep_rotations[i].hash_block_ = parms_.get_hash_block();
		//    bs = babystep[i]; 
		//    apply_galois(encrypted, bs, babystep_rotations[i]); 
		//}

		vector<Ciphertext> intermediate_vector(giantstep.size());;

		// Populate plain coefficients; in this case, the output has dimension giantstep
		// and the input has dimension babystep. So this is correct. 

		matrix_vector_product_plain(babystep_rotations, plain_coeffs, intermediate_vector, false);

		auto time_gs_start = chrono::high_resolution_clock::now();
		for (int j = 0; j < giantstep.size(); j++) {
			Ciphertext temp;
			apply_galois(intermediate_vector[j], giantstep[j], temp, true);
			//add(destination, intermediate_vector[j], destination);
			add(destination, temp, destination);
		}
		if (!ntt_output) {
			transform_from_ntt(destination);
		}
		auto time_gs_end = chrono::high_resolution_clock::now();
		chrono::microseconds time_pre = chrono::duration_cast<chrono::microseconds>(time_gs_end - time_gs_start);
		cout << "+giant step time : " << (double)time_pre.count() / 1000000 << " seconds" << endl;
		return;
	}

	void RNSEvaluator::rotations_weighted_sum_coeffs_ntt(const Ciphertext & encrypted, const vector<vector<Plaintext> > &coeffs_ntt, vector<uint64_t> babystep, vector<uint64_t> giantstep, Ciphertext & destination, bool ntt_output)
	{
		int coeff_count = parms_.poly_modulus().coeff_count();
		int coeff_mod_count = parms_.coeff_mod_array().size();
		int coeff_bit_count = coeff_mod_count * bits_per_uint64;
		int encrypted_count = encrypted.size();

		int n = coeff_count - 1;
		int m = 2 * n;
		BigPolyArray temp, encrypted_bs, temp2;

		destination.get_mutable_array().resize(encrypted_count, coeff_count, coeff_bit_count);
		destination.get_mutable_array().set_zero();
		destination.rns_hash_block_ = parms_.get_hash_block();

		// Compute all rotations in babystep.
		vector<Ciphertext> babystep_rotations(babystep.size());
		cout << "computing babystep rotations using hoisting ... " << endl;
		auto time_start = chrono::high_resolution_clock::now();
		apply_galois_many(encrypted, babystep, babystep_rotations, true);
		auto time_end = chrono::high_resolution_clock::now();
		chrono::microseconds online_time = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
		cout << "+rotations weighted sum hoisting time: " << (double)online_time.count() / 1000000 << " seconds" << endl;

		vector<Ciphertext> intermediate_vector(giantstep.size());;
		matrix_vector_product_plain_ntt(babystep_rotations, coeffs_ntt, intermediate_vector, false);

		auto time_gs_start = chrono::high_resolution_clock::now();
		for (int j = 0; j < giantstep.size(); j++) {
			Ciphertext temp;
			apply_galois(intermediate_vector[j], giantstep[j], temp, true);
			//add(destination, intermediate_vector[j], destination);
			add(destination, temp, destination);
		}
		if (!ntt_output) {
			transform_from_ntt(destination);
		}
		auto time_gs_end = chrono::high_resolution_clock::now();
		chrono::microseconds time_pre = chrono::duration_cast<chrono::microseconds>(time_gs_end - time_gs_start);
		cout << "+giant step time : " << (double)time_pre.count() / 1000000 << " seconds" << endl;
		return;
	}

	/**
	Recursive polynomial evaution methdod,
	First, f(x) = ((X^{k2^{step-1}}) + a(x)) * q(x) + X^(k2^{step-1}-1) + b(x)
	Second, Compute q9x) and X^(k2^{step-1}-1) + b(x) using recursive function all
	Third, This the step = 1, jus-t do dot product and return
	*/
	void RNSEvaluator::PatersonStockmeyer(vector<Ciphertext> &baby_step, vector<Ciphertext> &giant_step, Ciphertext &encrypted, vector<uint64_t> &coeff, Ciphertext &destination, int k, int step)
	{
		//Check input coeff degree = k* (2^step - 1)?
		int current_degree = k * (power(2, step) - 1);
		if (coeff.size() != current_degree + 1)
		{
			invalid_argument("degree is incorrect for PatersonStockmeyer method");
		}

		//If the degree is k, do just dot product
		if (step == 1) {
			Polynomial_Evaluation(baby_step, coeff, destination);
			return;
		}

		uint64_t plain_modulus_ = parms_.plain_modulus().value();

		// Set q(x) and r(x) to be a f(x) = q(x) * X^{k2^{step-1}} + r(x)
		int temp_degree = (k * power(2, step - 1)) - 1;
		int next_degree = k * (power(2, step - 1) - 1);
		vector<uint64_t> qx(next_degree + 1);
		vector<uint64_t> rx(temp_degree + 1);
		for (int i = 0; i < coeff.size(); i++)
		{
			if (i < temp_degree + 1)
			{
				rx[i] = coeff[i];
			}
			else
			{
				qx[i - temp_degree - 1] = coeff[i];
			}
		}

		// r(x) = r(x) - X^(k2^{step-1}-1)
		if (rx[next_degree] > 0)
		{
			rx[next_degree] -= 1;
		}
		else
		{
			rx[next_degree] += plain_modulus_ - 1;
		}

		// convert from vector to BigPoly
		BigPoly poly_qx(rx.size(), 64);
		for (int i = 0; i < qx.size(); i++) {
			set_uint_uint(&qx[i], 1, poly_qx[i].pointer());
		}
		for (int i = qx.size(); i < rx.size(); i++)
		{
			set_zero_uint(1, poly_qx[i].pointer());
		}

		// convert from vector to BigPoly
		BigPoly poly_rx(rx.size(), 64);
		for (int i = 0; i < rx.size(); i++) {
			set_uint_uint(&(rx[i]), 1, poly_rx[i].pointer());
		}

		// r(x) = a(x) * q(x) + b(x)
		BigPoly poly_ax(poly_rx.coeff_count(), parms_.plain_modulus().bit_count());
		BigPoly poly_bx(poly_rx.coeff_count(), parms_.plain_modulus().bit_count());
		divide_poly_poly_coeff_smallmod(poly_rx.pointer(), poly_qx.pointer(), poly_rx.coeff_count(), parms_.plain_modulus(), poly_ax.pointer(), poly_bx.pointer(), pool_);

		//Convert bigpoly to vector<uint64_t>
		vector<uint64_t> ax(poly_ax.significant_coeff_count());
		for (int i = 0; i < poly_ax.significant_coeff_count(); i++)
		{
			ax[i] = *(poly_ax[i].pointer());
		}

		//Evaluate ax using baby_step (dot product)
		Ciphertext encrypted_ax;
		Polynomial_Evaluation(baby_step, ax, encrypted_ax);

		//Compute (X^{k2^{step-1}}) + a(x) using addition
		Ciphertext encrypted_ax_add_power_of_two;
		add(encrypted_ax, giant_step[step - 2], encrypted_ax_add_power_of_two);

		//Compute bx polynomial = X^(k2^{step-1}-1) + b(x)
		vector<uint64_t> bx(next_degree + 1);
		for (int i = 0; i < poly_bx.significant_coeff_count(); i++)
		{
			bx[i] = *(poly_bx[i].pointer());
		}
		for (int i = poly_bx.significant_coeff_count(); i < next_degree; i++)
		{
			bx[i] = 0;
		}
		bx[next_degree] = 1;

		//Evaluate q(x) recursively
		Ciphertext encrypted_qx;
		PatersonStockmeyer(baby_step, giant_step, encrypted, qx, encrypted_qx, k, step - 1);

		//Evaluate X^(k2^{step-1}-1) + b(x) recursively
		Ciphertext encrypted_bx;
		PatersonStockmeyer(baby_step, giant_step, encrypted, bx, encrypted_bx, k, step - 1);

		//Compute ((X^{k2^{step-1}}) + a(x)) * q(x)
		Ciphertext encrypted_part;
		multiply(encrypted_ax_add_power_of_two, encrypted_qx, encrypted_part);
		relinearize(encrypted_part, encrypted_part, 2);

		//Compute f(x) = [((X^{k2^{step-1}}) + a(x)) * q(x)] + X^(k2^{step-1}-1) + b(x)
		add(encrypted_part, encrypted_bx, destination);
	}

	void RNSEvaluator::Polynomial_Evaluation(const vector<Ciphertext> &all_powers_encrypted, vector<uint64_t> &coeff, Ciphertext &destination)
	{
		// Check whether we have enough all powers
		if (all_powers_encrypted.size() < coeff.size())
		{
			throw invalid_argument("There are not enough powers_encrypted!!!");
		}

		// Extract encryption parameters.
		int coeff_count = parms_.poly_modulus().coeff_count();
		int coeff_mod_count = parms_.coeff_mod_array().size();
		int coeff_bit_count = coeff_mod_count * bits_per_uint64;
		int encrypted_count = all_powers_encrypted[0].size();

		// Initialize destination
		destination.get_mutable_array().resize(encrypted_count, coeff_count, coeff_bit_count);
		destination.get_mutable_array().set_zero();
		destination.rns_hash_block_ = parms_.get_hash_block();

		// Dot product
		int degree = coeff.size() - 1;
		Pointer temp(allocate_uint(coeff_count, pool_));

		for (int i = 0; i < degree + 1; i++)
		{
			if (coeff[i] != 0)
			{
				for (int j = 0; j < encrypted_count; j++)
				{
					const uint64_t *all_powers_encrypted_ptr = all_powers_encrypted[i].get_array().pointer(j);
					uint64_t *destination_ptr = destination.get_mutable_array().pointer(j);
					for (int k = 0; k < coeff_mod_count; k++)
					{
						multiply_poly_scalar_coeff_smallmod(all_powers_encrypted_ptr + k * coeff_count, coeff_count, &coeff[i], coeff_mod_array_[k], temp.get());
						add_poly_poly_coeff_smallmod(destination_ptr + k * coeff_count, temp.get(), coeff_count, coeff_mod_array_[k], destination_ptr + k * coeff_count);
					}
				}
			}
		}
	}

	void RNSEvaluator::combine(vector<Ciphertext> &encrypted, Ciphertext & destination, vector<Plaintext> &unit_vecs) {

		int coeff_count = parms_.poly_modulus().coeff_count();
		int coeff_mod_count = parms_.coeff_mod_array().size();
		int coeff_bit_count = coeff_mod_count * bits_per_uint64;
		int coeff_uint64_count = coeff_mod_count;
		int encrypted_ptr_increment = coeff_count * coeff_mod_count;
		int encrypted_count = encrypted.size();


		// Multiply old unit vector
		destination = multiply_plain(encrypted[0], unit_vecs[0]);
		for (int i = 0; i < encrypted.size(); i++)
		{
			Ciphertext temp = multiply_plain(encrypted[i], unit_vecs[i]);
			add(destination, temp, destination);
		}
		return;
	}


	void RNSEvaluator::expand(const Ciphertext &encrypted, vector<Ciphertext> & destination, vector<Plaintext> &unit_vecs, vector<uint64_t> rotations_index)
	{
		int coeff_count = parms_.poly_modulus().coeff_count();
		int coeff_mod_count = parms_.coeff_mod_array().size();
		int coeff_bit_count = coeff_mod_count * bits_per_uint64;
		int coeff_uint64_count = coeff_mod_count;
		int encrypted_ptr_increment = coeff_count * coeff_mod_count;
		int encrypted_count = encrypted.size();

		for (int i = 0; i < destination.size(); i++) {
			// Multiply old unit vector to get (m1, 0, ..., 0).
			Ciphertext temp = multiply_plain(encrypted, unit_vecs[i]);

			// log(numSlots_) rotations to send it to (m1, m1, ...)
			int lognumSlot = rotations_index.size();
			for (int j = 0; j < lognumSlot; j++)
			{
				Ciphertext rotate;
				rotate = apply_galois(temp, rotations_index[j]);
				temp = add(temp, rotate);
			}
			destination[i] = temp;
		}
		return;
	}

	void RNSEvaluator::compute_all_powers(const Ciphertext &encrypted, int degree, vector<Ciphertext> &destination)
	{
		// Verify parameters.
		if (encrypted.rns_hash_block_ != parms_.get_hash_block())
		{
			throw invalid_argument("encrypted is not valid for encryption parameters");
		}

		// Extract encryption parameters.
		int coeff_count = parms_.poly_modulus().coeff_count();
		int coeff_mod_count = parms_.coeff_mod_array().size();
		int coeff_bit_count = coeff_mod_count * bits_per_uint64;
		int encrypted_count = encrypted.size();

		// Resize destination and clear
		destination.resize(degree + 1);
		destination[0].get_mutable_array().resize(encrypted_count, coeff_count, coeff_bit_count);
		destination[0].get_mutable_array().set_zero();
		destination[0].rns_hash_block_ = parms_.get_hash_block();

		//Set destination[0] = (Delta, 0, ..., 0)
		BigUInt coeff_div_current_plain_modulus_(coeff_bit_count);
		Pointer temp(allocate_uint(coeff_mod_count, pool_));
		Pointer wide_plain_modulus(allocate_uint(coeff_mod_count, pool_));
		set_uint_uint(parms_.plain_modulus().pointer(), parms_.plain_modulus().uint64_count(), coeff_mod_count, wide_plain_modulus.get());
		divide_uint_uint(product_modulus_.get(), wide_plain_modulus.get(), coeff_mod_count, coeff_div_current_plain_modulus_.pointer(), temp.get(), pool_);
		set_uint_uint(coeff_div_current_plain_modulus_.pointer(), coeff_mod_count, destination[0].get_mutable_array().pointer(0));
		rns_decompose(destination[0].get_mutable_array().pointer(0));

		//Compute all X^{2^i} using square function
		int log_2_degree = (int)(log2((double)degree));
		vector<Ciphertext> all_two_powers_encrypted(log_2_degree + 1);
		all_two_powers_encrypted[0] = encrypted;
		for (int i = 1; i <= log_2_degree; i++)
		{
			square(all_two_powers_encrypted[i - 1], all_two_powers_encrypted[i]);
			relinearize(all_two_powers_encrypted[i], all_two_powers_encrypted[i]);
		}

		// Compute all X^{i} using multiplication tree
		for (int i = 1; i <= degree; i++)
		{
			int i1 = optimal_slpit(i);
			int i2 = i - i1;
			if (i1 == 0 || i2 == 0)
			{
				destination[i].get_mutable_array().resize(encrypted_count, coeff_count, coeff_bit_count);
				destination[i] = all_two_powers_encrypted[(int)log2((double)i)];
			}
			else
			{
				destination[i].get_mutable_array().resize(encrypted_count, coeff_count, coeff_bit_count);
				multiply(destination[i1], destination[i2], destination[i]);
				relinearize(destination[i], destination[i]);
			}
		}
		return;
	}

}
