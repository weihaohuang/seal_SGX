#include <algorithm>
#include <random>
#include "rnskeygenerator.h"
#include "util/uintcore.h"
#include "util/uintarith.h"
#include "util/uintarithmod.h"
#include "util/smallpolyarith.h"
#include "util/polyarithsmallmod.h"
#include "util/polyfftmultmod.h"
#include "util/randomtostd.h"
#include "util/clipnormal.h"
#include "util/polyextras.h"
#include "util/polycore.h"
#include "util/smallntt.h"

using namespace std;
using namespace seal::util;

namespace seal
{
    void RNSKeyGenerator::generate(int evaluation_keys_count, int hammingweight, bool auto_galois_key_gen)
    {
        // if decomposition bit count is zero, evaluation keys must be empty
        if (!qualifiers_.enable_relinearization && evaluation_keys_count != 0)
        {
            throw invalid_argument("cannot generate evaluation keys for specified encryption parameters");
        }

        // If already generated, reset everything.
        if (generated_)
        {
            evaluation_keys_.get_mutable_data().clear();
            secret_key_.get_mutable_poly().set_zero();
            public_key_.get_mutable_array().set_zero();
            generated_ = false;
        }

        // Extract encryption parameters.
        int coeff_count = parms_.poly_modulus().coeff_count();
        int coeff_uint64_count = coeff_mod_count_;
        int coeff_bit_count = bits_per_uint64 * coeff_mod_count_;

        unique_ptr<UniformRandomGenerator> random(random_generator_->create());

        // generate secret key
        uint64_t *secret_key = secret_key_.get_mutable_poly().pointer();
        set_poly_coeffs_zero_one_negone_rns(secret_key, random.get(), hammingweight);

        //generate public key: (pk[0],pk[1]) = ([-(as+e)]_q, a)

        // sample a uniformly at random
        // set pk[1] = a
        uint64_t *public_key_1 = public_key_.get_mutable_array().pointer(1);
        set_poly_coeffs_uniform_rns(public_key_1, random.get());

        // calculate a*s + e (mod q) and store in pk[0]
        if (qualifiers_.enable_ntt)
        {
            for (int i = 0; i < coeff_mod_count_; i++)
            {
                // transform the secret s into NTT representation. 
                smallntt_negacyclic_harvey(secret_key + (i * coeff_count), small_ntt_tables_[i], pool_);

                // transform the uniform random polynomial a into NTT representation. 
                smallntt_negacyclic_harvey(public_key_1 + (i * coeff_count), small_ntt_tables_[i], pool_);
            }
                Pointer noise(allocate_poly(coeff_count, coeff_uint64_count, pool_));
                set_poly_coeffs_normal_rns(noise.get(), random.get());
            for (int i = 0; i < coeff_mod_count_; i++)
            {
                // transform the noise e into NTT representation.
                smallntt_negacyclic_harvey(noise.get() + (i * coeff_count), small_ntt_tables_[i], pool_);

                dyadic_product_coeff_smallmod(secret_key + (i * coeff_count), public_key_1 + (i * coeff_count), coeff_count, parms_.coeff_mod_array()[i], public_key_.get_mutable_array().pointer(0) + (i * coeff_count), pool_);
                add_poly_poly_coeff_smallmod(noise.get() + (i * coeff_count), public_key_.get_mutable_array().pointer(0) + (i * coeff_count), coeff_count, parms_.coeff_mod_array()[i], public_key_.get_mutable_array().pointer(0) + (i * coeff_count));
            }
        }
        else
        {
            // This branch should never be reached
            throw logic_error("invalid encryption parameters");
        }

        // negate and set this value to pk[0]
        // pk[0] is now -(as+e) mod q
        for (int i = 0; i < coeff_mod_count_; i++)
        {
            negate_poly_coeff_smallmod(public_key_.get_mutable_array().pointer(0) + (i * coeff_count), coeff_count, parms_.coeff_mod_array()[i], parms_.coeff_mod_array()[i].uint64_count(), public_key_.get_mutable_array().pointer(0) + (i * coeff_count));
        }

        // Set the secret_key_array to have size 1 (first power of secret) 
        secret_key_array_.resize(1, coeff_count, coeff_bit_count);
        set_poly_poly(secret_key_.get_mutable_poly().pointer(), coeff_count, coeff_uint64_count, secret_key_array_.pointer(0));

        // Set the parameter hashes for public and secret key
        public_key_.get_rns_mutable_hash_block() = parms_.get_hash_block();
        secret_key_.get_rns_mutable_hash_block() = parms_.get_hash_block();
        
        // Secret and public keys have been generated
        generated_ = true;

        // generate the requested number of evaluation keys
        generate_rns_evaluation_keys(evaluation_keys_count);

		if (auto_galois_key_gen)
		{
			// generate the galois keys which is needed for apply_galois
			generate_rns_galois_keys_logn(galois_keys_);
		}
    }


    void RNSKeyGenerator::generate_rns_evaluation_keys(int count)
    {
        // if decomposition bit count is zero, evaluation keys must be empty
        if (!qualifiers_.enable_relinearization && count != 0)
        {
            throw invalid_argument("cannot generate evaluation keys for specified encryption parameters");
        }

        // Extract encryption parameters.
        int coeff_count = parms_.poly_modulus().coeff_count();
        int coeff_uint64_count = coeff_mod_count_;
        int coeff_bit_count = bits_per_uint64 * coeff_mod_count_;

        // Check to see if secret key and public key have been generated
        if (!generated_)
        {
            throw logic_error("cannot generate evaluation keys for unspecified secret key");
        }

        // Check to see if count is not negative
        if (count < 0)
        {
            throw invalid_argument("count must be non-negative");
        }

        // Check if the specified number of evaluation keys have already been generated.
        // This would be the case if count is less than the current evaluation_keys_.size().
        // In this case, do nothing.
        if (count <= evaluation_keys_.size())
        {
            return;
        }

        // In the constructor, evaluation_keys_ is initialized to have evaluation_keys_.size() = 0.
        // In a previous call to generate_evaluation_keys, evaluation_keys_ was only initialized to contain evaluation_keys_.size() entries.
        // Therefore need to initialize further if count > evaluation_keys_.size().
        int initial_evaluation_key_size = evaluation_keys_.size();
        int rns_evaluation_factors_count = rns_evaluation_factors_[0].size();
        for (int i = initial_evaluation_key_size; i < count; i++)
        {
            vector<pair<BigPolyArray, BigPolyArray> > temp_vec_pair;
            for (int j = 0; j < coeff_mod_count_; j++)
            {
                temp_vec_pair.emplace_back(BigPolyArray(rns_evaluation_factors_count, coeff_count, coeff_bit_count), BigPolyArray(rns_evaluation_factors_count, coeff_count, coeff_bit_count));
            }
            evaluation_keys_.get_mutable_data().emplace_back(temp_vec_pair);
        }

        unique_ptr<UniformRandomGenerator> random(random_generator_->create());

        // Create evaluation keys.
        Pointer noise(allocate_poly(coeff_count, coeff_uint64_count, pool_));
        Pointer power(allocate_uint(coeff_uint64_count, pool_));
        Pointer secret_key_power(allocate_poly(coeff_count, coeff_uint64_count, pool_));
        Pointer temp(allocate_poly(coeff_count, coeff_uint64_count, pool_));

        int poly_ptr_increment = coeff_count * coeff_uint64_count;

        // Make sure we have enough secret keys computed
        compute_rns_secret_key_array(count + 1);

        uint64_t *secret_key = secret_key_.get_mutable_poly().pointer();
        if (qualifiers_.enable_ntt)
        {
            // assume the secret key is already transformed into NTT form. 
            for (int k = initial_evaluation_key_size; k < count; k++)
            {
                set_poly_poly(secret_key_array_.pointer(0) + (k + 1)*poly_ptr_increment, coeff_count, coeff_uint64_count, secret_key_power.get());
                for (int l = 0; l < coeff_mod_count_; l++)
                {
                    //populate evaluate_keys_[k]
                    for (int i = 0; i < rns_evaluation_factors_[0].size(); i++)
                    {
                        //generate NTT(a_i) and store in evaluation_keys_[k][l].second[i]
                        uint64_t *eval_keys_second = evaluation_keys_.get_mutable_data()[k][l].second.pointer(i);
                        uint64_t *eval_keys_first = evaluation_keys_.get_mutable_data()[k][l].first.pointer(i);

                        set_poly_coeffs_uniform_rns(eval_keys_second, random.get());
                        for (int j = 0; j < coeff_mod_count_; j++)
                        {
                            smallntt_negacyclic_harvey(eval_keys_second + (j * coeff_count), small_ntt_tables_[j], pool_);

                            // calculate a_i*s and store in evaluation_keys_[k].first[i]
                            dyadic_product_coeff_smallmod(eval_keys_second + (j * coeff_count), secret_key + (j * coeff_count), coeff_count, parms_.coeff_mod_array()[j], eval_keys_first + (j * coeff_count), pool_);
                        }

                        //generate NTT(e_i) 
                        set_poly_coeffs_normal_rns(noise.get(), random.get());
                        for (int j = 0; j < coeff_mod_count_; j++)
                        {
                            smallntt_negacyclic_harvey(noise.get() + (j * coeff_count), small_ntt_tables_[j], pool_);

                            //add e_i into evaluation_keys_[k].first[i]
                            add_poly_poly_coeff_smallmod(noise.get() + (j * coeff_count), eval_keys_first + (j * coeff_count), coeff_count, parms_.coeff_mod_array()[j], eval_keys_first + (j * coeff_count));

                            // negate value in evaluation_keys_[k].first[i]
                            negate_poly_coeff_smallmod(eval_keys_first + (j * coeff_count), coeff_count, parms_.coeff_mod_array()[j], parms_.coeff_mod_array()[j].uint64_count(), eval_keys_first + (j * coeff_count));

                            //multiply w^i * s^(k+2)
                            uint64_t rns_evaluation_factors_mod_ = (l == j) * rns_evaluation_factors_[l][i];
                            multiply_poly_scalar_coeff_smallmod(secret_key_power.get() + (j * coeff_count), coeff_count, &rns_evaluation_factors_mod_, parms_.coeff_mod_array()[j], temp.get() + (j * coeff_count));

                            //add w^i . s^(k+2) into evaluation_keys_[k].first[i]
                            add_poly_poly_coeff_smallmod(eval_keys_first + (j * coeff_count), temp.get() + (j * coeff_count), coeff_count, parms_.coeff_mod_array()[j], eval_keys_first + (j * coeff_count));
                        }
                    }
                }
            }
        }
        else
        {
            // This branch should never be reached
            throw logic_error("invalid encryption parameters");
        }
        
        // Set the parameter hash
        evaluation_keys_.get_mutable_hash_block() = parms_.get_hash_block();
    }

    void RNSKeyGenerator::set_poly_coeffs_zero_one_negone_rns(uint64_t *poly, UniformRandomGenerator *random, int hammingweight) const
    {
        int coeff_count = parms_.poly_modulus().coeff_count();
        int coeff_mod_count = parms_.coeff_mod_array().size();

        RandomToStandardAdapter engine(random);
        uniform_int_distribution<int> dist(-1, 1);
		
		int bound;
		if (hammingweight == 0)
		{
			bound = coeff_count;
		}
		else
		{
			bound = hammingweight;
		}

		int count = 0;
        for (int i = 0; i < coeff_count - 1; i++)
        {
            int rand_index = dist(engine);
            if (rand_index == 1 && count < bound)
            {
                for (int j = 0; j < coeff_mod_count; j++)
                {
                    poly[i + (j * coeff_count)] = 1;
                }
				count++;
            }
            else if (rand_index == -1 && count < bound)
            {
                for (int j = 0; j < coeff_mod_count; j++)
                {
                    poly[i + (j * coeff_count)] = parms_.coeff_mod_array()[j].value() - 1;
                }
				count++;
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

    void RNSKeyGenerator::set_poly_coeffs_normal_rns(uint64_t *poly, UniformRandomGenerator *random) const
    {
        int coeff_count = parms_.poly_modulus().coeff_count();
        int coeff_mod_count = parms_.coeff_mod_array().size();
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
                    poly[i + (j * coeff_count)] = parms_.coeff_mod_array()[j].value() - static_cast<uint64_t>(noise);
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

    /*Set the coeffs of a BigPoly to be uniform modulo coeff_mod*/
    void RNSKeyGenerator::set_poly_coeffs_uniform_rns(uint64_t *poly, UniformRandomGenerator *random)
    {
        //get parameters
        int coeff_count = parms_.poly_modulus().coeff_count();
        int coeff_mod_count = parms_.coeff_mod_array().size();

        // set up source of randomness which produces random things of size 32 bit
        RandomToStandardAdapter engine(random);

        // pointer to increment to fill in all coeffs
        uint32_t * value_ptr = reinterpret_cast<uint32_t *> (poly);

        // number of 32 bit blocks to cover all coeffs
        int significant_value_uint32_count = (coeff_count - 1) * 2;
        int value_uint32_count = significant_value_uint32_count + 2;

        // Sample randomness to all but last coefficient
        for (int j = 0; j < coeff_mod_count; j++)
        {
            for (int i = 0; i < value_uint32_count; i++)
            {
                *value_ptr++ = i < significant_value_uint32_count ? engine() : 0;
            }
        }
        
        // when poly is fully populated, reduce all coefficient modulo coeff_modulus
        for (int i = 0; i < coeff_mod_count; i++)
        {
            smallmod_poly_coeffs(poly + (i * coeff_count), coeff_count, parms_.coeff_mod_array()[i], pool_);
        }
        
    }

    RNSKeyGenerator::RNSKeyGenerator(const RNSContext &context, const MemoryPoolHandle &pool) :
        pool_(pool), parms_(context.get_parms()),
        random_generator_(parms_.random_generator()),
        qualifiers_(context.get_qualifiers())
    {
        // Verify parameters
        if (!qualifiers_.parameters_set)
        {
            throw invalid_argument("encryption parameters are not set correctly");
        }

        // Resize encryption parameters to consistent size.
        int coeff_count = parms_.poly_modulus().coeff_count();
        int poly_coeff_uint64_count = parms_.poly_modulus().coeff_uint64_count();
        int coeff_bit_count = context.coeff_modulus().bit_count();
        coeff_mod_count_ = context.coeff_mod_array().size();
        int coeff_products_uint64_count = coeff_mod_count_;

        for (int i = 0; i < coeff_mod_count_; i++)
        {
            small_ntt_tables_.emplace_back(util::SmallNTTTables(pool_));
            
        }
        small_ntt_tables_ = context.small_ntt_tables_;

        // Initialize public and secret key.
        public_key_.get_mutable_array().resize(2, coeff_count, coeff_mod_count_ * bits_per_uint64);
        secret_key_.get_mutable_poly().resize(coeff_count, coeff_mod_count_ * bits_per_uint64);

        // Initialize evaluation_factors_, if required
        if (qualifiers_.enable_relinearization)
        {
            populate_rns_evaluation_factors();
        }

        // Initialize evaluation keys to empty
        evaluation_keys_.get_mutable_data().clear();

        // Initialize moduli.
        polymod_ = PolyModulus(parms_.poly_modulus().pointer(), coeff_count, poly_coeff_uint64_count);

        // Secret key and public key have not been generated
        generated_ = false;
    }

    RNSKeyGenerator::RNSKeyGenerator(const RNSContext &context, const SecretKey &secret_key, const PublicKey &public_key, const RNSEvaluationKeys &evaluation_keys, const MemoryPoolHandle &pool) :
        pool_(pool), parms_(context.get_parms()), qualifiers_(context.get_qualifiers()),
        evaluation_keys_(evaluation_keys), public_key_(public_key), secret_key_(secret_key),
        random_generator_(parms_.random_generator()), generated_(false)
    {
        // Verify parameters
        if (!qualifiers_.parameters_set)
        {
            throw invalid_argument("encryption parameters are not set correctly");
        }

        // Decomposition bit count should only be zero if evaluation keys are empty
        if (!qualifiers_.enable_relinearization && evaluation_keys.size() != 0)
        {
            throw invalid_argument("evaluation keys are not valid for encryption parameters");
        }

        // Verify parameters.
        if (secret_key.get_rns_hash_block() != parms_.get_hash_block())
        {
            throw invalid_argument("secret key is not valid for encryption parameters");
        }
        if (public_key.get_rns_hash_block() != parms_.get_hash_block())
        {
            throw invalid_argument("public key is not valid for encryption parameters");
        }
        if (evaluation_keys.size() > 0)
        {
            if (evaluation_keys.get_hash_block() != parms_.get_hash_block())
            {
                throw invalid_argument("evaluation keys are not valid for encryption parameters");
            }
        }

        // Resize encryption parameters to consistent size.
        int coeff_count = parms_.poly_modulus().coeff_count();
        int poly_coeff_uint64_count = parms_.poly_modulus().coeff_uint64_count();
        int coeff_mod_count_ = context.coeff_mod_array().size();

        // Set SmallNTTTables
        small_ntt_tables_.resize(coeff_mod_count_, pool_);
        small_ntt_tables_ = context.small_ntt_tables_;

        // Initialize public and secret key.
        public_key_.get_mutable_array().resize(2, coeff_count, coeff_mod_count_);
        secret_key_.get_mutable_poly().resize(coeff_count, coeff_mod_count_);

        // Initialize evaluation_factors_, if required
        if (qualifiers_.enable_relinearization)
        {
            populate_rns_evaluation_factors();
        }

        // Initialize moduli.
        polymod_ = PolyModulus(parms_.poly_modulus().pointer(), coeff_count, poly_coeff_uint64_count);

        // Secret key and public key are generated
        generated_ = true;
    }

    // rns_evaluation_factors[i][j] = 2^(w*j) * hat-q_i mod q_i
    void RNSKeyGenerator::populate_rns_evaluation_factors()
    {
        rns_evaluation_factors_.clear();

        int coeff_mod_count = parms_.coeff_mod_array().size();
        int rns_evaluation_factors_size = 64 / parms_.decomposition_bit_count();
        // Initialize evaluation_factors_
        rns_evaluation_factors_.resize(coeff_mod_count_);
        uint64_t power_of_w = 1ULL << parms_.decomposition_bit_count();

        // Compute hat-q_i mod q_i
        vector<uint64_t> coeff_prod_mod(coeff_mod_count_);
        for (int i = 0; i < coeff_mod_count_; i++)
        {
            coeff_prod_mod[i] = 1;
            for (int j = 0; j < coeff_mod_count_; j++)
            {
                if (i != j)
                {
                    multiply_uint64_smallmod(&coeff_prod_mod[i], parms_.coeff_mod_array()[j].pointer(), parms_.coeff_mod_array()[i], &coeff_prod_mod[i]);
                }
            }
        }

        for (int i = 0; i < coeff_mod_count_; i++)
        {
            uint64_t current_rns_evaluation_factor = coeff_prod_mod[i];
            for (int j = 0; j < rns_evaluation_factors_size; j++)
            {
                rns_evaluation_factors_[i].emplace_back(current_rns_evaluation_factor);
                //multiply 2^w mod q_i
                if (j < rns_evaluation_factors_size - 1)
                {
                    multiply_uint64_smallmod(&current_rns_evaluation_factor, &power_of_w, parms_.coeff_mod_array()[i], &current_rns_evaluation_factor);
                }
            }
        }
    }

    const SecretKey &RNSKeyGenerator::secret_key() const
    {
        if (!generated_)
        {
            throw logic_error("encryption keys have not been generated");
        }
        return secret_key_;
    }

    const PublicKey &RNSKeyGenerator::public_key() const
    {
        if (!generated_)
        {
            throw logic_error("encryption keys have not been generated");
        }
        return public_key_;
    }

    const RNSEvaluationKeys &RNSKeyGenerator::evaluation_keys() const
    {
        if (!generated_)
        {
            throw logic_error("encryption keys have not been generated");
        }
        if (evaluation_keys_.size() == 0)
        {
            throw logic_error("no evaluation keys have been generated");
        }
        return evaluation_keys_;
    }

    void RNSKeyGenerator::compute_rns_secret_key_array(int max_power)
    {
        int old_count = secret_key_array_.size();
        int new_count = max(max_power, secret_key_array_.size());

        if (old_count == new_count)
        {
            return;
        }

        int coeff_count = parms_.poly_modulus().coeff_count();
        int coeff_bit_count = bits_per_uint64 * coeff_mod_count_;
        
        // Compute powers of secret key until max_power
        secret_key_array_.resize(new_count, coeff_count, coeff_bit_count);

        int poly_ptr_increment = coeff_count * coeff_mod_count_;
        uint64_t *prev_poly_ptr = secret_key_array_.pointer(old_count - 1);
        uint64_t *next_poly_ptr = prev_poly_ptr + poly_ptr_increment;

        if (qualifiers_.enable_ntt)
        {
            // Since all of the key powers in secret_key_array_ are already NTT transformed, to get the next one 
            // we simply need to compute a dyadic product of the last one with the first one [which is equal to NTT(secret_key_)].
            for (int i = old_count; i < new_count; i++)
            {
                for (int j = 0; j < coeff_mod_count_; j++)
                {
                    dyadic_product_coeff_smallmod(prev_poly_ptr + (j * coeff_count), secret_key_array_.pointer(0) + (j * coeff_count), coeff_count, parms_.coeff_mod_array()[j], next_poly_ptr + (j * coeff_count), pool_);
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

	//////////////////////////////////////////////////////////////

	void RNSKeyGenerator::generate(const SecretKey &secret_key, int evaluation_keys_count)
	{
		// if decomposition bit count is zero, evaluation keys must be empty
		if (!qualifiers_.enable_relinearization && evaluation_keys_count != 0)
		{
			throw invalid_argument("cannot generate evaluation keys for specified encryption parameters");
		}

		// Verity parameter of secretkey
		if (secret_key.rns_hash_block_ != parms_.get_hash_block())
		{
			throw invalid_argument("SecretKey does not have correct paramter");
		}

		// If already generated, reset everything.
		if (generated_)
		{
			evaluation_keys_.get_mutable_data().clear();
			secret_key_.get_mutable_poly().set_zero();
			public_key_.get_mutable_array().set_zero();
			generated_ = false;
		}

		// Extract encryption parameters.
		int coeff_count = parms_.poly_modulus().coeff_count();
		int coeff_uint64_count = coeff_mod_count_;
		int coeff_bit_count = bits_per_uint64 * coeff_mod_count_;

		unique_ptr<UniformRandomGenerator> random(random_generator_->create());

		// Copy secret_key (input secret_key is in NTT form)
		secret_key_ = secret_key;
		uint64_t *secret_key_ptr = secret_key_.get_mutable_poly().pointer();

		//generate public key: (pk[0],pk[1]) = ([-(as+e)]_q, a)

		// sample a uniformly at random
		// set pk[1] = a
		uint64_t *public_key_1 = public_key_.get_mutable_array().pointer(1);
		set_poly_coeffs_uniform_rns(public_key_1, random.get());

		// calculate a*s + e (mod q) and store in pk[0]
		if (qualifiers_.enable_ntt)
		{
			for (int i = 0; i < coeff_mod_count_; i++)
			{
				// transform the uniform random polynomial a into NTT representation. 
				smallntt_negacyclic_harvey(public_key_1 + (i * coeff_count), small_ntt_tables_[i], pool_);
			}
			Pointer noise(allocate_poly(coeff_count, coeff_uint64_count, pool_));
			set_poly_coeffs_normal_rns(noise.get(), random.get());
			for (int i = 0; i < coeff_mod_count_; i++)
			{
				// transform the noise e into NTT representation.
				smallntt_negacyclic_harvey(noise.get() + (i * coeff_count), small_ntt_tables_[i], pool_);

				dyadic_product_coeff_smallmod(secret_key_ptr + (i * coeff_count), public_key_1 + (i * coeff_count), coeff_count, parms_.coeff_mod_array()[i], public_key_.get_mutable_array().pointer(0) + (i * coeff_count), pool_);
				add_poly_poly_coeff_smallmod(noise.get() + (i * coeff_count), public_key_.get_mutable_array().pointer(0) + (i * coeff_count), coeff_count, parms_.coeff_mod_array()[i], public_key_.get_mutable_array().pointer(0) + (i * coeff_count));
			}
		}
		else
		{
			// This branch should never be reached
			throw logic_error("invalid encryption parameters");
		}

		// negate and set this value to pk[0]
		// pk[0] is now -(as+e) mod q
		for (int i = 0; i < coeff_mod_count_; i++)
		{
			negate_poly_coeff_smallmod(public_key_.get_mutable_array().pointer(0) + (i * coeff_count), coeff_count, parms_.coeff_mod_array()[i], parms_.coeff_mod_array()[i].uint64_count(), public_key_.get_mutable_array().pointer(0) + (i * coeff_count));
		}

		// Set the secret_key_array to have size 1 (first power of secret) 
		secret_key_array_.resize(1, coeff_count, coeff_bit_count);
		set_poly_poly(secret_key_.get_mutable_poly().pointer(), coeff_count, coeff_uint64_count, secret_key_array_.pointer(0));

		// Set the parameter hashes for public and secret key
		public_key_.get_rns_mutable_hash_block() = parms_.get_hash_block();
		secret_key_.get_rns_mutable_hash_block() = parms_.get_hash_block();

		// Secret and public keys have been generated
		generated_ = true;

		// generate the requested number of evaluation keys
		generate_rns_evaluation_keys(evaluation_keys_count);
	}

	void RNSKeyGenerator::generate_rns_galois_keys(uint64_t galois_elt)
	{
		// Extract encryption parameters.
		int coeff_count = parms_.poly_modulus().coeff_count();
		int coeff_mod_count = parms_.coeff_mod_array().size();
		int coeff_bit_count = bits_per_uint64 * coeff_mod_count_;

		// Verify coprime conditions. 
		if (galois_elt % 2 == 0 || galois_elt >= 2 * coeff_count)
		{
			throw invalid_argument("galois element is not valid");
		}

		// Check to see if secret key and public key have been generated
		if (!generated_)
		{
			throw logic_error("cannot generate evaluation keys for unspecified secret key");
		}

		// Check to see if galois key for this galois_elt is generated or not
		if (galois_keys_.has_key(galois_elt))
		{
			return;
		}

		// Rotate secret key for each coeff_mod_array
		Pointer rotated_secret_key(allocate_poly(coeff_count, coeff_mod_count_, pool_));
		for (int i = 0; i < coeff_mod_count_; i++)
		{
			permute_ntt_poly_smallmod(secret_key_.get_poly().pointer() + i * coeff_count, coeff_count - 1, galois_elt, rotated_secret_key.get() + i * coeff_count);
		}

		// Initialize galois key
		int rns_evaluation_factors_count = rns_evaluation_factors_[0].size();
		vector<pair<BigPolyArray, BigPolyArray> > temp_vec_pair(coeff_mod_count_);
		for (int i = 0; i < coeff_mod_count_; i++)
		{
			temp_vec_pair[i].first.resize(rns_evaluation_factors_count, coeff_count, coeff_bit_count);
			temp_vec_pair[i].second.resize(rns_evaluation_factors_count, coeff_count, coeff_bit_count);
		}
		galois_keys_.get_mutable_data().emplace(galois_elt, temp_vec_pair);

		unique_ptr<UniformRandomGenerator> random(random_generator_->create());

		// Create evaluation keys.
		Pointer noise(allocate_poly(coeff_count, coeff_mod_count, pool_));
		Pointer power(allocate_uint(coeff_mod_count, pool_));
		Pointer temp(allocate_poly(coeff_count, coeff_mod_count, pool_));

		int poly_ptr_increment = coeff_count * coeff_mod_count;

		const uint64_t *secret_key = secret_key_.get_poly().pointer();
		if (qualifiers_.enable_ntt)
		{
			for (int l = 0; l < coeff_mod_count_; l++)
			{
				//populate evaluate_keys_[k]
				for (int i = 0; i < rns_evaluation_factors_[0].size(); i++)
				{
					//generate NTT(a_i) and store in evaluation_keys_[k][l].second[i]
					uint64_t *eval_keys_second = galois_keys_.get_mutable_data()[galois_elt][l].second.pointer(i);
					uint64_t *eval_keys_first = galois_keys_.get_mutable_data()[galois_elt][l].first.pointer(i);

					set_poly_coeffs_uniform_rns(eval_keys_second, random.get());
					for (int j = 0; j < coeff_mod_count_; j++)
					{
						// a_i in NTT form
						smallntt_negacyclic_harvey(eval_keys_second + (j * coeff_count), small_ntt_tables_[j], pool_);
						// calculate a_i*s and store in evaluation_keys_[k].first[i]
						dyadic_product_coeff_smallmod(eval_keys_second + (j * coeff_count), secret_key + (j * coeff_count), coeff_count, parms_.coeff_mod_array()[j], eval_keys_first + (j * coeff_count), pool_);
					}

					//generate NTT(e_i) 
					set_poly_coeffs_normal_rns(noise.get(), random.get());
					for (int j = 0; j < coeff_mod_count_; j++)
					{
						smallntt_negacyclic_harvey(noise.get() + (j * coeff_count), small_ntt_tables_[j], pool_);

						//add NTT(e_i) into evaluation_keys_[k].first[i]
						add_poly_poly_coeff_smallmod(noise.get() + (j * coeff_count), eval_keys_first + (j * coeff_count), coeff_count, parms_.coeff_mod_array()[j], eval_keys_first + (j * coeff_count));

						// negate value in evaluation_keys_[k].first[i]
						negate_poly_coeff_smallmod(eval_keys_first + (j * coeff_count), coeff_count, parms_.coeff_mod_array()[j], parms_.coeff_mod_array()[j].uint64_count(), eval_keys_first + (j * coeff_count));

						//multiply w^i * rotated_secret_key
						uint64_t rns_evaluation_factors_mod_ = (l == j) * rns_evaluation_factors_[l][i];
						multiply_poly_scalar_coeff_smallmod(rotated_secret_key.get() + (j * coeff_count), coeff_count, &rns_evaluation_factors_mod_, parms_.coeff_mod_array()[j], temp.get() + (j * coeff_count));

						//add w^i * rotated_secret_key into evaluation_keys_[k].first[i]
						add_poly_poly_coeff_smallmod(eval_keys_first + (j * coeff_count), temp.get() + (j * coeff_count), coeff_count, parms_.coeff_mod_array()[j], eval_keys_first + (j * coeff_count));
					}
				}
			}
		}
		else
		{
			// This branch should never be reached
			throw logic_error("invalid encryption parameters");
		}

		// Set the parameter hash
		galois_keys_.rns_hash_block_ = parms_.get_hash_block();
	}

	void RNSKeyGenerator::generate_rns_galois_keys(uint64_t galois_elt, RNSGaloisKeys &galois_keys)
	{
		// Extract encryption parameters.
		int coeff_count = parms_.poly_modulus().coeff_count();
		int coeff_mod_count = parms_.coeff_mod_array().size();
		int coeff_bit_count = bits_per_uint64 * coeff_mod_count_;

		// Verify coprime conditions. 
		if (galois_elt % 2 == 0 || galois_elt >= 2 * coeff_count)
		{
			throw invalid_argument("galois element is not valid");
		}

		// Check to see if secret key and public key have been generated
		if (!generated_)
		{
			throw logic_error("cannot generate evaluation keys for unspecified secret key");
		}

		// Check to see if galois key for this galois_elt is generated or not
		if (galois_keys.has_key(galois_elt))
		{
			return;
		}

		// Rotate secret key for each coeff_mod_array
		Pointer rotated_secret_key(allocate_poly(coeff_count, coeff_mod_count_, pool_));
		for (int i = 0; i < coeff_mod_count_; i++)
		{
			permute_ntt_poly_smallmod(secret_key_.get_poly().pointer() + i * coeff_count, coeff_count - 1, galois_elt, rotated_secret_key.get() + i * coeff_count);
		}

		// Initialize galois key
		int rns_evaluation_factors_count = rns_evaluation_factors_[0].size();
		vector<pair<BigPolyArray, BigPolyArray> > temp_vec_pair(coeff_mod_count_);
		for (int i = 0; i < coeff_mod_count_; i++)
		{
			temp_vec_pair[i].first.resize(rns_evaluation_factors_count, coeff_count, coeff_bit_count);
			temp_vec_pair[i].second.resize(rns_evaluation_factors_count, coeff_count, coeff_bit_count);
		}
		galois_keys.get_mutable_data().emplace(galois_elt, temp_vec_pair);

		unique_ptr<UniformRandomGenerator> random(random_generator_->create());

		// Create evaluation keys.
		Pointer noise(allocate_poly(coeff_count, coeff_mod_count, pool_));
		Pointer power(allocate_uint(coeff_mod_count, pool_));
		Pointer temp(allocate_poly(coeff_count, coeff_mod_count, pool_));

		int poly_ptr_increment = coeff_count * coeff_mod_count;

		const uint64_t *secret_key = secret_key_.get_poly().pointer();
		if (qualifiers_.enable_ntt)
		{
			for (int l = 0; l < coeff_mod_count_; l++)
			{
				//populate evaluate_keys_[k]
				for (int i = 0; i < rns_evaluation_factors_[0].size(); i++)
				{
					//generate NTT(a_i) and store in evaluation_keys_[k][l].second[i]
					uint64_t *eval_keys_second = galois_keys.get_mutable_data()[galois_elt][l].second.pointer(i);
					uint64_t *eval_keys_first = galois_keys.get_mutable_data()[galois_elt][l].first.pointer(i);

					set_poly_coeffs_uniform_rns(eval_keys_second, random.get());
					for (int j = 0; j < coeff_mod_count_; j++)
					{
						// a_i in NTT form
						smallntt_negacyclic_harvey(eval_keys_second + (j * coeff_count), small_ntt_tables_[j], pool_);
						// calculate a_i*s and store in evaluation_keys_[k].first[i]
						dyadic_product_coeff_smallmod(eval_keys_second + (j * coeff_count), secret_key + (j * coeff_count), coeff_count, parms_.coeff_mod_array()[j], eval_keys_first + (j * coeff_count), pool_);
					}

					//generate NTT(e_i) 
					set_poly_coeffs_normal_rns(noise.get(), random.get());
					for (int j = 0; j < coeff_mod_count_; j++)
					{
						smallntt_negacyclic_harvey(noise.get() + (j * coeff_count), small_ntt_tables_[j], pool_);

						//add NTT(e_i) into evaluation_keys_[k].first[i]
						add_poly_poly_coeff_smallmod(noise.get() + (j * coeff_count), eval_keys_first + (j * coeff_count), coeff_count, parms_.coeff_mod_array()[j], eval_keys_first + (j * coeff_count));

						// negate value in evaluation_keys_[k].first[i]
						negate_poly_coeff_smallmod(eval_keys_first + (j * coeff_count), coeff_count, parms_.coeff_mod_array()[j], parms_.coeff_mod_array()[j].uint64_count(), eval_keys_first + (j * coeff_count));

						//multiply w^i * rotated_secret_key
						uint64_t rns_evaluation_factors_mod_ = (l == j) * rns_evaluation_factors_[l][i];
						multiply_poly_scalar_coeff_smallmod(rotated_secret_key.get() + (j * coeff_count), coeff_count, &rns_evaluation_factors_mod_, parms_.coeff_mod_array()[j], temp.get() + (j * coeff_count));

						//add w^i * rotated_secret_key into evaluation_keys_[k].first[i]
						add_poly_poly_coeff_smallmod(eval_keys_first + (j * coeff_count), temp.get() + (j * coeff_count), coeff_count, parms_.coeff_mod_array()[j], eval_keys_first + (j * coeff_count));
					}
				}
			}
		}
		else
		{
			// This branch should never be reached
			throw logic_error("invalid encryption parameters");
		}

		// Set the parameter hash
		galois_keys.rns_hash_block_ = parms_.get_hash_block();
	}

	void RNSKeyGenerator::generate_rns_galois_keys_logn(RNSGaloisKeys &galois_key)
	{
		int coeff_count = parms_.poly_modulus().coeff_count();
		int n = coeff_count - 1;
		int m = n << 1;
		int logn = get_power_of_two(n);

		// Generate galois key for m - 1 (X -> X^{m-1})
		generate_rns_galois_keys(m - 1, galois_key);

		// Generate galois key for power of 3 mod m (X -> X^{3^{2^k}})
		uint64_t two_power_of_three = 3;
		for (int i = 0; i < logn - 1; i++)
		{
			generate_rns_galois_keys(two_power_of_three, galois_key);
			two_power_of_three *= two_power_of_three;
			two_power_of_three %= m;
		}
	}

	const RNSGaloisKeys &RNSKeyGenerator::galois_keys() const
	{
		if (!generated_)
		{
			throw logic_error("encryption keys have not been generated");
		}
		if (galois_keys_.size() == 0)
		{
			throw logic_error("no galois keys have been generated");
		}
		return galois_keys_;
	}

}
