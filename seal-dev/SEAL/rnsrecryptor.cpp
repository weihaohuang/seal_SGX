#include <algorithm>
#include <stdexcept>
#include <chrono>
#include "evaluator.h"
#include "primes.h"
#include "util/common.h"
#include "util/uintcore.h"
#include "util/uintarith.h"
#include "util/uintarithmod.h"
#include "util/polycore.h"
#include "util/polyarith.h"
#include "util/smallpolyarith.h"
#include "util/polyarithmod.h"
#include "util/polyfftmult.h"
#include "util/polyfftmultmod.h"
#include "bigpoly.h"
#include "rnsrecryptor.h"
#include "ciphertext.h"
#include "util/ntt.h"
#include "bigpolyarith.h"

#include "rnskeygenerator.h"
#include "util/expolycrt.h"
#include <fstream>

//#define _DEBUG
#define NUM_THREAD 3

using namespace seal;
using namespace seal::util;

namespace seal {
	RNSRecryptor::RNSRecryptor(const RNSContext &oldcontext, SecretKey &secret_key, const RNSEvaluationKeys &evaluation_key, const RNSGaloisKeys &galois_key, const std::string &hex_poly, int e, bool nonBatch, const MemoryPoolHandle &pool) :
		pool_(pool), parms_(oldcontext.get_parms()), nonBatch_(nonBatch)
	{
		// Set e_
		e_ = e;

		// Set oldevaluator
		oldcontext_.reset(new RNSContext(parms_, pool_));
		oldevaluator_.reset(new RNSEvaluator(*(oldcontext_), evaluation_key, pool_));
		oldevaluator_->set_galois_keys(galois_key);

		// Compute new parameters
		cout << "\nCompute new paramters...";
		compute_new_parameters();
		cout << "done" << endl;

		// Set SEALContext
		newcontext_.reset(new RNSContext(newparms_, pool_));

		// Save secret_key
		secret_key_ = secret_key;
		secret_key_.rns_hash_block_ = newparms_.get_hash_block();
		decryptor_.reset(new RNSDecryptor(*(newcontext_), secret_key_));

		if (!nonBatch)
		{
			// Compute cosets
			cout << "Compute cosets...";
			compute_cosets_classic();
			cout << "done" << endl;

			babystep_ = cosets_;
			int elt = 1;
			giantstep_.clear();
			for (int i = 0; i < extensionDegree_; i++) {
				giantstep_.push_back(elt);
				elt *= p_;
				elt %= m_;
			}
		}

		// Set Exring and initialize
		cout << "Set Exring and ExPolyCRTBuilder...";
		BigUInt prime(64, p_);
		ExRing_ = ExRing::acquire_ring(prime, e_, hex_poly);
		ExRing_->init_frobe_table(m_);
		ExPolyCRTBuilder_.reset(new ExPolyCRTBuilder(ExRing_, get_power_of_two(n_)));
		oldExRing_ = ExRing::acquire_ring(prime, r_, hex_poly); 
		oldExRing_->init_frobe_table(m_);
		oldExPolyCRTBuilder_.reset(new ExPolyCRTBuilder(oldExRing_, get_power_of_two(n_)));
		cout << "done" << endl;

		///////////////////////////////////////////////////////////////////

		vector<int64_t> result = xgcd(3, m_);
		if (result[0] != 1)
		{
			throw "Inverse does not exist!";
		}
		if (result[1] < 0)
		{
			result[1] += m_;
		}
		uint64_t three_inv_mod_m = result[1];

		int log_slot_num = get_power_of_two(ExPolyCRTBuilder_->slot_count());
		power_of_three_coset_index.resize(log_slot_num);
		for (int i = 0; i < log_slot_num; i++)
		{
			power_of_three_coset_index[log_slot_num - i - 1] = three_inv_mod_m;
			three_inv_mod_m *= three_inv_mod_m;
			three_inv_mod_m %= m_;
		}

		///////////////////////////////////////////////////////////////////

		int logn = get_power_of_two(n_);
		for (int galois_elt = 1; galois_elt < m_; galois_elt += 2)
		{
			vector<int> index(n_);
			for (int i = 0; i < n_; i++)
			{
				int reversed = bit_reversal_general(i, logn);
				int index_raw = (2 * reversed + 1) * galois_elt;
				index_raw %= m_;
				index[i] = bit_reversal_general((index_raw - 1) >> 1, logn);
			}
			permute_index_.emplace(galois_elt, index);
		}

		unit_plain_.resize(numSlots_);
		for (int i = 0; i < numSlots_; i++)
		{
			vector<ExRingElement> unit_vector(numSlots_);
			for (int j = 0; j < numSlots_; j++)
			{
				if (i == j)
				{
					unit_vector[j] = ExRingElement(ExRing_, "1");
				}
				else
				{
					unit_vector[j] = ExRingElement(ExRing_, "0");
				}
			}
			unit_plain_[i] = ExPolyCRTBuilder_->compose(unit_vector);
		}

		old_unit_plain_.resize(numSlots_);
		for (int i = 0; i < numSlots_; i++)
		{
			vector<ExRingElement> unit_vector(numSlots_);
			for (int j = 0; j < numSlots_; j++)
			{
				if (i == j)
				{
					unit_vector[j] = ExRingElement(oldExRing_, "1");
				}
				else
				{
					unit_vector[j] = ExRingElement(oldExRing_, "0");
				}
			}
			old_unit_plain_[i] = oldExPolyCRTBuilder_->compose(unit_vector);
		}

		///////////////////////////////////////////////////////////////////

		// Generate keys
		cout << "Start to generate galois keys..." << endl;

		generate_recrypt_key_and_new_evaluator(secret_key_);
		
		cout << "... Galois Keys and Recrypt Key are computed!" << endl;

		// Compute lift and remainlsd poly
		compute_lift_removelsd_poly();

		// Compute transformation polys
		if(!nonBatch)
			compute_transformation_polys(hex_poly);
		cout << "... Transformation polys are computed!" << endl;

		if (nonBatch_ == false)
		{
			string filename = "linpoly-materials-n" + std::to_string(n_) + "p" + std::to_string(p_) + "e" + std::to_string(e_) + ".txt";
			compute_unpack_polys(filename);
			cout << "unpack polys computed!" << endl;
		}
	}

	void RNSRecryptor::compute_new_parameters() {
		// Notice!!! parms_.plain_modulus_base() = 2 case is not supported
		if (parms_.plain_modulus_base() == 2)
		{
			invalid_argument("Recrypt can not be done with this parameter");
		}

		// Extract paramters from parms_
		p_ = parms_.plain_modulus_base();
		n_ = parms_.poly_modulus().coeff_count() - 1;
		r_ = parms_.plain_modulus_exponent();
		m_ = n_ << 1;

		// Compute numSlots and extensionDegree
		int k = 1;
		int q = p_;
		while (q % m_ != 1) {
			q = q*q % m_;
			k <<= 1;
		}

		numSlots_ = (int)(n_ / k);
		extensionDegree_ = k;

		int logn = get_power_of_two(n_);

		uint64_t newt = power(p_, e_);
		vector<int64_t> result = xgcd(n_, newt);
		if (result[0] != 1)
		{
			throw "Inverse does not exist!";
		}
		if (result[1] < 0)
		{
			result[1] += newt;
		}
		nInverse_ = result[1];

		uint64_t new_plain_modulus_uint64_ = power(p_, e_);
		SmallModulus new_plain_modulus_;
		new_plain_modulus_ = new_plain_modulus_uint64_;

		newparms_.set_coeff_modulus(parms_.coeff_mod_array());
		newparms_.set_plain_modulus(new_plain_modulus_);
		newparms_.set_poly_modulus(parms_.poly_modulus());
		newparms_.set_decomposition_bit_count(62);
		newparms_.set_plain_modulus_base(p_);
		newparms_.set_plain_modulus_exponent(e_);
	}

	void RNSRecryptor::generate_recrypt_key_and_new_evaluator(SecretKey &secret_key)
	{
		//Extract new parameters
		int coeff_count = newparms_.poly_modulus().coeff_count();
		int coeff_bit_count = newparms_.coeff_mod_array().size()*64;
		int coeff_mod_count = newparms_.coeff_mod_array().size();

		//Generate Key
		RNSKeyGenerator newgenerator(*(newcontext_));
		newgenerator.generate(secret_key, 1);

		newevaluator_.reset(new RNSEvaluator(*(newcontext_), newgenerator.evaluation_keys(), pool_));

		// Generate Galoise key for additon between slots
		if (nonBatch_)
		{
			for (int i = 0; i < power_of_three_coset_index.size(); i++)
			{
				newgenerator.generate_rns_galois_keys(power_of_three_coset_index[i], newevaluator_->galois_keys_);
			}
			vector<int> galois_elts;
			int n = newparms_.poly_modulus().coeff_count() - 1;
			int logn = get_power_of_two(n);
			for (int i = 0; i < logn; i++)
			{
				galois_elts.push_back((n + power(2, i)) / power(2, i));
			}
			for (int i = 0; i < galois_elts.size(); i++)
			{
				newgenerator.generate_rns_galois_keys(galois_elts[i], newevaluator_->galois_keys_);
			}
			cout << "... Galois Keys are generated (" << power_of_three_coset_index.size() + galois_elts.size() << ")" << endl;
		}
		else
		{
			// Generate the keys corresponding to baby step giant step
			// TODO: Uncommenting the next two lines triggers some operator delete fault. 
			for (int i = 0; i < babystep_.size(); i++)
			{
				newgenerator.generate_rns_galois_keys(babystep_[i], newevaluator_->galois_keys_);
			}
			for (int i = 0; i < giantstep_.size(); i++)
			{
				newgenerator.generate_rns_galois_keys(giantstep_[i], newevaluator_->galois_keys_);
			}
			cout << "... Galois Keys are generated (" << babystep_.size() + giantstep_.size() + power_of_three_coset_index.size() << ")" << endl;
		}
		

		//Extract plain_modulus parameter
		int plain_bit_count = newparms_.plain_modulus().bit_count();
		
		//Calculate - 1 (mod new_plain_modulus)
		uint64_t new_plain_modulus_minus_one_ = newparms_.plain_modulus().value() - 1;

		//Transform secret key to nonNTT
		if (newcontext_->get_qualifiers().enable_ntt)
		{
			for (int i = 0; i < coeff_mod_count; i++)
			{
				inverse_smallntt_negacyclic_harvey(secret_key.get_mutable_poly().pointer() + i * coeff_count, newcontext_->small_ntt_tables_[i], pool_);
			}
		}

		//Change the secret_key to new plain modulus (to make recrypt_key)
		Pointer secret_key_new_plain_mod_pointer(allocate_poly(coeff_count, 1, pool_));
		const uint64_t *secret_key_ptr = secret_key.get_poly().pointer();
		for (int j = 0; j < coeff_count - 1; j++)
		{
			if (*(secret_key_ptr + j) == 1)
			{
				set_uint(1, 1, secret_key_new_plain_mod_pointer.get() + j);
			}
			else if (*(secret_key_ptr + j) == 0)
			{
				set_uint(0, 1, secret_key_new_plain_mod_pointer.get() + j);
			}
			else
			{
				set_uint(new_plain_modulus_minus_one_, 1, secret_key_new_plain_mod_pointer.get() + j);
			}
		}
		set_zero_uint(1, secret_key_new_plain_mod_pointer.get() + (coeff_count - 1));
		
		// Generate and save recrypt key
		RNSEncryptor newencryptor(*(newcontext_), newgenerator.public_key());
		Plaintext secret_key_new_plain_mod_plain_;
		secret_key_new_plain_mod_plain_.get_poly().resize(coeff_count, plain_bit_count);
		set_poly_poly(secret_key_new_plain_mod_pointer.get(), coeff_count, 1, secret_key_new_plain_mod_plain_.get_poly().pointer());
		recrypt_key_ = newencryptor.rns_encrypt(secret_key_new_plain_mod_plain_);

		// Save newevaluator
		// newevaluator_->set_galois_keys(newgenerator.galois_keys());

		//Go back to NTT form if needed
		if (oldcontext_->get_qualifiers().enable_ntt == true)
		{
			for (int i = 0; i < coeff_mod_count; i++)
			{
				smallntt_negacyclic_harvey(secret_key.get_mutable_poly().pointer(), newcontext_->small_ntt_tables_[i], pool_);
			}
		}
	}

	void RNSRecryptor::compute_cosets_classic() {
		cosets_.clear();
		vector<uint64_t> classes = conjugate_classes(m_, p_);
		vector<uint64_t> orders = multiplicative_orders(classes, m_);
		for (int j = 1; j < classes.size(); j++) {
			if (classes[j] == j) cosets_.push_back(j);
			if (cosets_.size() >= numSlots_) break;
		}
		classes_ = classes;

#ifdef _DEBUG
		cout << "cosets_ = ";
		for (int i = 0; i < numSlots_; i++) {
			cout << cosets_[i] << " ";
		}
		cout << endl;
#endif // _DEBUG
		return;
	}

	void RNSRecryptor::compute_cosets() {
		// TODO: check that these powers of g does form the entire coset. 
		int g = 5;
		cosets_.clear();
		uint64_t elt = 1;
		for (int i = 0; i < numSlots_; i++) {
			cosets_.push_back(elt);
			elt *= g;
			elt %= m_;
		}

#ifdef _DEBUG
		cout << "cosets_ = ";
		for (int i = 0; i < numSlots_; i++) {
			cout << cosets_[i] << " ";
		}
		cout << endl;
#endif // _DEBUG

	}

	int RNSRecryptor::coset_permutation_action(int j, int i) {
		uint64_t gj = cosets_[j];
		uint64_t gi = cosets_[i];
		uint64_t element = gj*gi;
		element %= m_;
		int x = std::distance(cosets_.begin(), std::find(cosets_.begin(), cosets_.end(), classes_[element]));
		if (cosets_[x] != classes_[element]) throw;
		return x;
	}

	//void Recryptor::compute_exring(const std::string &hex_poly)
	//{
	//	if (exring_ != nullptr)
	//	{
	//		return;
	//	}
	//	BigUInt p(64, p_);
	//	std::shared_ptr<ExRing> E = ExRing::acquire_ring(p, e_, hex_poly);
	//	E->init_frobe_table();
	//	exring_ = E;
	//}

	void RNSRecryptor::compute_transformation_polys(const std::string &hex_poly)
	{
		BigUInt p(64, p_);
		//std::shared_ptr<ExRing> E = ExRing::acquire_ring(p, e_, hex_poly);
		//E->init_frobe_table(m_); // TODO: this takes too much time (80% of construct)
		ExRingElement zeta(ExRing_, "1x^1");
		ExRingElement ninv(ExRing_);
		int exring_coeff_count = ninv.ex_ring().get()->coeff_count();
		int exring_coeff_uint64_count = ninv.ex_ring().get()->coeff_uint64_count();
		int exring_coeff_bit_count = ninv.ex_ring().get()->coeff_modulus().significant_bit_count();

		set_uint_uint(&nInverse_, exring_coeff_uint64_count, ninv.pointer(0));

		vector<ExRingElement> power_trace(m_);
		for (int i = 0; i < m_; i++) {
			ExRingElement power = zeta^i;
			power_trace[i] = power.trace();
		}

		vector<BigUInt> g1coeffs;
		for (int k = 0; k < n_; k++) {
			ExRingElement coeff(ExRing_);
			for (int i = 0; i < numSlots_; i++) {
				int exponent = (m_ - k - i) % m_;
				int modified_exponent = exponent * cosets_[i] % m_;
				coeff += power_trace[modified_exponent];
			}
			coeff *= ninv;
			BigUInt temp(exring_coeff_bit_count, coeff.pointer(0));
			g1coeffs.push_back(temp);
		}
		BigPoly g1_ =  biguint_vector_to_bigpoly(g1coeffs);

		hpolys_.clear();
		for (int j = 0; j < numSlots_; j++)
		{
			vector<BigUInt> newhcoeffs;
			for (int k = 0; k < n_; k++)
			{
				ExRingElement coeff(ExRing_);
				for (int i = 0; i < numSlots_; i++)
				{
					int exponent = (m_ - k + coset_permutation_action(j, i)) % m_;
					int modified_exponent = exponent * cosets_[i] % m_;
					coeff += power_trace[modified_exponent];
				}
				coeff *= ninv;
				BigUInt temp(exring_coeff_bit_count, coeff.pointer(0));
				newhcoeffs.push_back(temp);
			}
			hpolys_.push_back(biguint_vector_to_bigpoly(newhcoeffs));
		}

		if (nonBatch_ == false)
		{
			vector<BigPoly> pre_linear_coeffs_;
			pre_linear_coeffs_.clear();
			SmallModulus plainmod = newparms_.plain_modulus();
			for (int i = 0; i < extensionDegree_; i++)
			{
				vector<BigUInt> coeffs(n_);
				for (int kk = 0; kk< n_; kk++)
				{
					coeffs[kk].resize(64);
					coeffs[kk].set_zero();
				}
				for (int j = 0; j < extensionDegree_; j++)
				{
					int exponent_raw = (1 - (2 * i + 1)*numSlots_)*j; //FIXME;
					int exponent = exponent_raw % n_;
					if (exponent < 0) exponent += n_;
					int multiples = (exponent_raw - exponent) / n_;
					if (multiples % 2 == 0)
					{
						coeffs[exponent] = 1;
					}
					else
					{
						coeffs[exponent] = plainmod.value() - 1; // the new plain modulus minus one.
					}
				}
				pre_linear_coeffs_.push_back(biguint_vector_to_bigpoly(coeffs));
			}
			BigUInt plain_modulus_(newparms_.plain_modulus().bit_count(), newparms_.plain_modulus().value());
			BigPoly g1squared;
			BigPolyArith arith;
			arith.multiply(g1_, g1_, newparms_.poly_modulus(), plain_modulus_, g1squared);
			prelinear_polys_.resize(giantstep_.size());
			for (int i = 0; i < prelinear_polys_.size(); i++)
			{
				prelinear_polys_[i].resize(babystep_.size());
			}
			BigPoly tempg = g1_;
			BigPoly temp;
			Modulus plainMod(newparms_.plain_modulus().pointer(), 1, pool_);
			for (int i = 0; i < n_; i++)
			{
				auto index = decompose_bs_gs(m_, (2 * i + 1) % m_, babystep_, giantstep_);
				// Verify destination size.
				prelinear_polys_[index.first][index.second].get_poly().resize(n_ + 1, plainMod.significant_bit_count());
				nussbaumer_multiply_poly_poly_coeffmod(tempg.pointer(), pre_linear_coeffs_[i % extensionDegree_].pointer(), get_power_of_two(n_), plainMod, prelinear_polys_[index.first][index.second].get_poly().pointer(), pool_);
				if (i < n_ - 1)
				{
					arith.multiply(g1squared, tempg, newparms_.poly_modulus(), plain_modulus_, temp);
					tempg = temp;
				}
			}
		}
		else
		{
			BigPolyArith arith;
			BigUInt plain_modulus_(newparms_.plain_modulus().bit_count(), newparms_.plain_modulus().value());

			BigPoly g1 = g1_;
			BigPoly g1squared;
			BigPoly const_ninverse(newparms_.poly_modulus().coeff_count(), newparms_.plain_modulus().bit_count());
			const_ninverse.set_zero();
			const_ninverse[0] = nInverse_;
			arith.multiply(g1, g1, newparms_.poly_modulus(), plain_modulus_, g1squared);

			vector<BigPoly> g_polys;
			g_polys.push_back(g1);
			for (int i = 1; i < n_; i++)
			{
				BigPoly temp;
				arith.multiply(g1squared, g_polys.back(), newparms_.poly_modulus(), plain_modulus_, temp);
				g_polys.push_back(temp);
			}
			for (int i = 0; i < g_polys.size(); i++) {
				BigPoly temp;
				arith.multiply(g_polys[i], const_ninverse, newparms_.poly_modulus(), plain_modulus_, temp);
				g_polys[i].duplicate_from(temp);
			}

			g_polys_ptxt_.resize(n_);
			for (int i = 0; i < n_; i++)
			{
				g_polys_ptxt_[i] = g_polys[i];
			}

			h_polys_ptxt_.resize(n_);
			for (int i = 0; i < n_; i++)
			{
				h_polys_ptxt_[i].get_poly().resize(newparms_.poly_modulus().coeff_count(), newparms_.plain_modulus().bit_count());
				h_polys_ptxt_[i].get_poly().set_zero();
			}

			for (int i = 0; i < numSlots_; i++)
			{
				int gal = cosets_[i];
				h_polys_ptxt_[(gal - 1) / 2].get_poly() = hpolys_[i];
			}

			g1_ntt_ = g_polys[0];
			g1_squared_ntt_ = g1squared;
			newevaluator_->transform_to_ntt(g1_ntt_);
			newevaluator_->transform_to_ntt(g1_squared_ntt_);

			/*int coeff_count = newparms_.poly_modulus().coeff_count();
			g_poly_ptxt_plain_coeffs_ntt_.resize(giantstep_.size());
			for (int j = 0; j < giantstep_.size(); j++) {
				int gs = giantstep_[j];
				g_poly_ptxt_plain_coeffs_ntt_[j].resize(babystep_.size());
				for (int i = 0; i < babystep_.size(); i++) {
					int bs = babystep_[i];
					int coeff_index = (bs*gs % m_ - 1) / 2;
					g_poly_ptxt_plain_coeffs_ntt_[j][i].get_poly().resize(coeff_count, newparms_.plain_modulus().bit_count());
					apply_galois_poly_smallmod(g_polys_ptxt_[coeff_index].get_poly().pointer(), coeff_count - 1, mod_inverse(gs, m_), newparms_.plain_modulus(), g_poly_ptxt_plain_coeffs_ntt_[j][i].get_poly().pointer());
					set_zero_uint(1, get_poly_coeff(g_poly_ptxt_plain_coeffs_ntt_[j][i].get_poly().pointer(), coeff_count - 1, 1));
					newevaluator_->transform_to_ntt(g_poly_ptxt_plain_coeffs_ntt_[j][i]);
				}
			}*/

		}
		return;
	}

	void RNSRecryptor::compute_lift_removelsd_poly()
	{
		int num = e_ - r_;
		lift_poly_.resize(num);
		for (int i = 0; i < num; i++) {
			lift_poly_[i].resize(p_ + 1);
			lift_poly_[i] = compute_lift_uint64(p_, i + 1); //TODO: for p_ is large case compute_lift poly takes much time (we need to use some recursive interpolyation)
		}

		remainlsd_poly_.resize(e_ - 1);
		for (int i = 0; i < e_ - 1; i++) {
			int degree = (i + 1)*(p_ - 1) + 1;
			remainlsd_poly_[i].resize(degree + 1);
			if (i >= r_ - 2) remainlsd_poly_[i] = compute_remainlsd_uint64(p_, i + 2);
		}

#ifdef _DEBUG
		for (int i = 0; i < num; i++)
		{
			cout << i << "-th lift poly : ";
			for (int j = 0; j < lift_poly_[i].size(); j++)
			{
				cout << lift_poly_[i][j] << " ";
			}
			cout << endl;
		}
		for (int i = 0; i < e_ - 1; i++)
		{
			cout << i << "-th remainlsd poly : ";
			for (int j = 0; j < remainlsd_poly_[i].size(); j++)
			{
				cout << remainlsd_poly_[i][j] << " ";
			}
			cout << endl;
		}
#endif // _DEBUG

	}

	void RNSRecryptor::InverseLinearTransformation(const vector<Ciphertext> &encrypted, Ciphertext &destination)
	{
		int coeff_count = parms_.poly_modulus().coeff_count();
		int coeff_mod_count = parms_.coeff_mod_array().size();
		int coeff_bit_count = coeff_mod_count * bits_per_uint64;
		int coeff_uint64_count = coeff_mod_count;
		int encrypted_count = encrypted.size();
		vector<BigPolyArray> encrypted_array(encrypted.size());
		for (int i = 0; i < encrypted_array.size(); i++)
		{
			encrypted_array[i] = encrypted[i].get_array();
		}
		destination.get_mutable_array().resize(encrypted_count, coeff_count, coeff_bit_count);
		destination.get_mutable_array().set_zero();
		destination.rns_hash_block_ = newparms_.get_hash_block();

		int shift;
		BigPoly hcoeff;
		Ciphertext temp2;
		temp2.get_mutable_array().resize(2, encrypted_array[0].coeff_count(), encrypted_array[0].coeff_bit_count());
		temp2.hash_block_ = newparms_.get_hash_block();

		// TODO: use Pointer for rotated_poly and temp
		BigPolyArray rotated_poly;
		rotated_poly.resize(2, encrypted_array[0].coeff_count(), encrypted_array[0].coeff_bit_count());
		BigPolyArray temp;
		temp.resize(2, encrypted_array[0].coeff_count(), encrypted_array[0].coeff_bit_count());

		int n = n_;
		int ninv = nInverse_;
		int m = m_;
		int d = extensionDegree_;
		int k = numSlots_;

		vector<Plaintext> hpolys_ntt(k);
		for (int i = 0; i < k; i++)
		{
			newevaluator_->transform_to_ntt(Plaintext(hpolys_[i]), hpolys_ntt[i]);
		}
		
		// Online part
		auto time_start = chrono::high_resolution_clock::now();
		vector<Ciphertext> partialsums(k);
		for (int i = 0; i < k; i++)
		{
			partialsums[i].get_mutable_array().resize(2, encrypted_array[0].coeff_count(), encrypted_array[0].coeff_bit_count());
			partialsums[i].hash_block_ = newparms_.get_hash_block();
			int gal = cosets_[i];
			vector<int64_t> result = xgcd(gal, m);
			if (result[0] != 1)
			{
				throw "Inverse does not exist!";
			}
			if (result[1] < 0) {
				result[1] += m;
			}
			int galinv = result[1];
			int k_adjusted = k * galinv;
			k_adjusted %= m;
			for (int j = 0; j < d; j++)
			{
				shift = j * k_adjusted;
				for (int l = 0; l < coeff_mod_count; l++){
					negacyclic_shift_poly_smallmod(encrypted_array[j].pointer(0) + l * coeff_count, coeff_count - 1, shift, temp.pointer(0) + l * coeff_count, newparms_.coeff_mod_array()[l]);
					negacyclic_shift_poly_smallmod(encrypted_array[j].pointer(1) + l * coeff_count, coeff_count - 1, shift, temp.pointer(1) + l * coeff_count, newparms_.coeff_mod_array()[l]);
					add_poly_poly_coeff_smallmod(partialsums[i].get_array().pointer(0) + l * coeff_count, temp.pointer(0) + l * coeff_count, coeff_count, newparms_.coeff_mod_array()[l], partialsums[i].get_mutable_array().pointer(0) + l * coeff_count);
					add_poly_poly_coeff_smallmod(partialsums[i].get_array().pointer(1) + l * coeff_count, temp.pointer(1) + l * coeff_count, coeff_count, newparms_.coeff_mod_array()[l], partialsums[i].get_mutable_array().pointer(1) + l * coeff_count);
				}
			}
			Ciphertext partialsum_rotated;
			// partial sum = gal(c1) + gal(c2) *x^k + ... + gal(cd)*x^((d-1)k)
			partialsums[i].rns_hash_block_ = newparms_.get_hash_block();
			newevaluator_->apply_galois(partialsums[i], gal, partialsum_rotated, true);
			temp2 = newevaluator_->multiply_plain_ntt(partialsum_rotated, hpolys_ntt[i]);
			newevaluator_->add(destination, temp2, destination);
		}
		newevaluator_->transform_from_ntt(destination);
		auto time_end = chrono::high_resolution_clock::now();
		chrono::microseconds online_time = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
		cout << "+ inverse linear transform online time : " << (double)online_time.count() / 1000000 << " seconds" << endl;
		return;
	}

	//////////////////////////// Non-Batch version ///////////////////////////////////
	void RNSRecryptor::LinearTransformation_nonBatch(const Ciphertext &encrypted, Ciphertext &destination, int num_Batch, int num_thread)
	{
		auto time_start = chrono::high_resolution_clock::now();

		int coeff_count = newparms_.poly_modulus().coeff_count();
		int coeff_mod_count = newparms_.coeff_mod_array().size();
		int coeff_bit_count = coeff_mod_count * bits_per_uint64;
		int encrypted_count = encrypted.size();

		if (num_Batch > numSlots_)
		{
			throw invalid_argument("number of batch should be smaller than numSlots");
		}

		int bound;
		if (num_Batch == 0)
		{
			bound = numSlots_;
		}
		else
		{
			bound = num_Batch;
		}

		// Suppose input is Enc(m_0 + m_1*X + ... m_{r-1}*X^{r-1})
		// Make a vector of ciphertext (each plaintext of ciphertexts's constant term is m_i)
		vector<Ciphertext> rotated_encrypted(bound);
		rotated_encrypted[0] = encrypted;
		int n = coeff_count - 1;

		for (int i = 1; i < bound; i++)
		{
			// Multiply X^{n-i}
			string temp_string = "1x^" + to_string(n - i);
			Plaintext power_of_X = temp_string;
			rotated_encrypted[i] = newevaluator_->multiply_plain(encrypted, power_of_X);
			rotated_encrypted[i] = newevaluator_->negate(rotated_encrypted[i]);
		}
		vector<int> galois_elts;
		int logn = get_power_of_two(n);
		for (int i = 0; i < logn; i++)
		{
			galois_elts.push_back((n + power(2, i)) / power(2, i));
		}

		// Make a compose of unit
		//vector<Plaintext> unit_plain(bound);
		//for (int i = 0; i < bound; i++)
		//{
		//	vector<ExRingElement> unit(numSlots_);
		//	for (int j = 0; j < numSlots_; j++)
		//	{
		//		if (i == j)
		//		{
		//			unit[j] = ExRingElement(ExRing_, "1");
		//		}
		//		else
		//		{
		//			unit[j] = ExRingElement(ExRing_, "0");
		//		}
		//	}
		//	unit_plain[i] = ExPolyCRTBuilder_->compose(unit);
		//}

		Plaintext nInv_ptxt(BigPoly(1, newparms_.plain_modulus().bit_count(), &nInverse_));
		atomic<int> index = ATOMIC_VAR_INIT(0);
		auto mthread = [&]()
		{
			while (true)
			{
				int my_index = index++;
				if (my_index >= bound)
					break;
				// Function for each index
				for (int i = 0; i < logn; i++)
				{
					Ciphertext temp;
					temp = newevaluator_->apply_galois(rotated_encrypted[my_index], galois_elts[i]);
					rotated_encrypted[my_index] = newevaluator_->add(temp, rotated_encrypted[my_index]);
				}
				rotated_encrypted[my_index] = newevaluator_->multiply_plain(rotated_encrypted[my_index], nInv_ptxt);
				rotated_encrypted[my_index] = newevaluator_->multiply_plain(rotated_encrypted[my_index], unit_plain_[my_index]);
			}
		};

		vector<thread> pool;
		for (int i = 0; i < num_thread; i++)
		{
			pool.emplace_back(thread(mthread));
		}

		for (int i = 0; i < num_thread; i++)
		{
			pool[i].join();
		}
		pool.clear();

		//Add all rotated_encrypted
		destination.get_mutable_array().resize(encrypted_count, coeff_count, coeff_bit_count);
		destination.get_mutable_array().set_zero();
		destination.rns_hash_block_ = newparms_.get_hash_block();

		for (int i = 0; i < bound; i++)
		{
			newevaluator_->add(destination, rotated_encrypted[i], destination);
		}

		auto time_end = chrono::high_resolution_clock::now();
		cout << "... LT_nonBatch time : " << (double)chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() / 1000000. << " seconds" << endl;
	}


	void RNSRecryptor::InverseLinearTransformation_nonBatch(const Ciphertext &encrypted, Ciphertext &destination)
	{
		auto time_start = chrono::high_resolution_clock::now();
		newevaluator_->rotations_weighted_sum(encrypted, h_polys_ptxt_, babystep_, giantstep_, destination);
		destination.hash_block_ = newparms_.get_hash_block();
		auto time_end = chrono::high_resolution_clock::now();
		cout << "... invLT_nonBatch time : " << (double)chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() / 1000000. << " seconds" << endl;
		return;
	}

	void RNSRecryptor::DigitExtraction(Ciphertext &encrypted, Ciphertext &destination)
	{
		long p = p_;
		long e = e_;
		long r = e_ - r_;

		if (e < 2)
		{
			throw invalid_argument("paramter e should be larger than 1");
		}

		// Extract encryption parameters.
		int coeff_count = newparms_.poly_modulus().coeff_count(); // Fixed to newparms_, Kay
		int coeff_mod_count = newparms_.coeff_mod_array().size();
		int encrypted_count = encrypted.size();

		// Find k such that p^k > (e-1)(p-1)+1
		long k = 1;
		long power_of_p = 1;
		while (1)
		{
			power_of_p *= p;
			if (power_of_p < (e - 1)*(p - 1) + 1)
			{
				k++;
			}
			else break;
		}

		//Initialize destination
		destination.hash_block_ = newparms_.get_hash_block();

		// Temp ciphertexts for computation
		vector<vector<Ciphertext> > lift_encrypted;
		lift_encrypted.resize(r);
		for (int i = 0; i < r; i++)
		{
			lift_encrypted[i].resize(k - 1);
		}
		vector<Ciphertext> remainlsd_encrypted;
		remainlsd_encrypted.resize(r);

		Ciphertext cipher;
		cipher = encrypted;
		destination = encrypted;

		uint64_t plain_modulus_ = newparms_.plain_modulus().value();
		SmallModulus eval_plain_modulus_(plain_modulus_);
		for (int i = 0; i < r; i++)
		{
			if (i < r - 1)
			{
				newevaluator_->Polynomial_Evaluation(cipher, lift_poly_[0], lift_encrypted[i][0], eval_plain_modulus_);

#ifdef _DEBUG
				cout << "Noise bubget after lift : " << decryptor_->invariant_noise_budget(lift_encrypted[i][0]) << endl;
				cout << "Msg after lift: ";
				for (int j = 0; j < ExPolyCRTBuilder_->slot_count(); j++)
				{
					cout << ExPolyCRTBuilder_->decompose(decryptor_->rns_decrypt(lift_encrypted[i][0]))[j].to_string() << ", ";
				}
				cout << endl;
#endif

			}
			newevaluator_->Polynomial_Evaluation(cipher, remainlsd_poly_[e - i - 2], remainlsd_encrypted[i], eval_plain_modulus_);

#ifdef _DEBUG
			cout << "Noise bubget after remainlsd : " << decryptor_->invariant_noise_budget(remainlsd_encrypted[i]) << endl;
			cout << "Msg after remainlsd: ";
			for (int j = 0; j < ExPolyCRTBuilder_->slot_count(); j++)
			{
				cout << ExPolyCRTBuilder_->decompose(decryptor_->rns_decrypt(remainlsd_encrypted[i]))[j].to_string() << ", ";
			}
			cout << endl;
#endif

			for (int j = 0; j < k - 2; j++)
			{
				if (i + j < r - 1)
				{
					newevaluator_->Polynomial_Evaluation(lift_encrypted[i][j], lift_poly_[j + 1], lift_encrypted[i][j + 1], eval_plain_modulus_);

#ifdef _DEBUG
					cout << "Noise bubget after lift : " << decryptor_->invariant_noise_budget(lift_encrypted[i][j + 1]) << endl;
#endif

				}
			}
			if (i < r - 1)
			{
				cipher = encrypted;
				for (int j = 0; j < i + 1; j++)
				{
					if (i - j > k - 2)
					{
						newevaluator_->sub(cipher, remainlsd_encrypted[j], cipher);
					}
					else
					{
						newevaluator_->sub(cipher, lift_encrypted[j][i - j], cipher);
					}
				}
			}
			newevaluator_->sub(destination, remainlsd_encrypted[i], destination);
			plain_modulus_ /= p;
			eval_plain_modulus_ = plain_modulus_;
		}
	}

	void RNSRecryptor::recrypt(Ciphertext &encrypted, Ciphertext &desitination)
	{
		// Extract parameters
		int coeff_count = newparms_.poly_modulus().coeff_count();
		int coeff_mod_count = newparms_.coeff_mod_array().size();
		int encrypted_size = encrypted.size();
		int encrypted_ptr_increment = coeff_count * coeff_mod_count;

		if (encrypted_size > 2)
		{
			throw invalid_argument("Input of recrypt should be ciphertext with size 2");
		}

		// Copy encrypted
		Pointer encrypted_copy(allocate_poly(encrypted_size * coeff_count, coeff_mod_count, pool_));
		for (int i = 0; i < encrypted_size; i++)
		{
			set_poly_poly(encrypted.get_array().pointer(i), coeff_count, coeff_mod_count, encrypted_copy.get() + i * encrypted_ptr_increment);
		}

		// Compose the encrypted_copy
		for (int i = 0; i < encrypted_size; i++)
		{
			oldevaluator_->rns_compose(encrypted_copy.get() + i * encrypted_ptr_increment);
		}

		// Multiply newparms_.plain_modulus() and divide in the coeff_products_array
		uint64_t new_plain_modulus = newparms_.plain_modulus().value();
		int big_coeff_mod_count = coeff_mod_count + 1;
		Pointer temp(allocate_zero_uint(big_coeff_mod_count, pool_));
		Pointer quotient(allocate_zero_uint(big_coeff_mod_count, pool_));
		Pointer wide_coeff_products_array(allocate_uint(big_coeff_mod_count, pool_));
		Pointer encrypted_after_MS_ptr(allocate_zero_uint(encrypted_size * coeff_count, pool_));
		set_uint_uint(oldcontext_->coeff_modulus().pointer(), coeff_mod_count, big_coeff_mod_count, wide_coeff_products_array.get());

		for (int i = 0; i < encrypted_size; i++)
		{
			const uint64_t *encrypted_coeff_ptr = encrypted_copy.get() + i * encrypted_ptr_increment;
			for (int j = 0; j < coeff_count - 1; j++)
			{
				multiply_uint_uint64(encrypted_coeff_ptr, coeff_mod_count, new_plain_modulus, big_coeff_mod_count, temp.get());
				divide_uint_uint_inplace(temp.get(), wide_coeff_products_array.get(), big_coeff_mod_count, quotient.get(), pool_);
				set_uint_uint(quotient.get(), big_coeff_mod_count, 1, encrypted_after_MS_ptr.get() + i * coeff_count + j);
				encrypted_coeff_ptr += coeff_mod_count;
			}
		}

		//Dot product
		Ciphertext encrypted_after_DP;
		BigPoly ptxt_c0_poly(coeff_count - 1, 1, encrypted_after_MS_ptr.get());
		BigPoly ptxt_c1_poly(coeff_count - 1, 1, encrypted_after_MS_ptr.get() + coeff_count);
		Plaintext ptxt_c0(ptxt_c0_poly);
		Plaintext ptxt_c1(ptxt_c1_poly);
		newevaluator_->multiply_plain(recrypt_key_, ptxt_c1, encrypted_after_DP);
		newevaluator_->add_plain(encrypted_after_DP, ptxt_c0, encrypted_after_DP);

#ifdef _DEBUG
		Plaintext plain = decryptor_->rns_decrypt(encrypted_after_DP);
		cout << "Message after DP : " << plain.to_string() << endl;
		cout << "Noise budget after DP : " << decryptor_->invariant_noise_budget(encrypted_after_DP) << endl;
		uint64_t Delta = power(p_, e_ - r_);
		uint64_t t = power(p_, r_);
		cout << " Delta =  " << Delta << endl;
		cout << "underlying message = ... ";
		for (int i = 0; i < plain.get_poly().coeff_count() - 1; i++)
		{
			BigUInt coeff = plain[i];
			BigUInt plaincoeff = ((coeff + (Delta / 2)) / Delta) % t;
			if (plaincoeff != 0) {
				cout << i << "-th coefficient = " << plaincoeff.to_dec_string() << endl;
			}
		}
		cout << "*****" << endl;
#endif

		//LinearTransformation
		vector<Ciphertext> encrypted_after_LT(extensionDegree_);
		//LinearTransformation(encrypted_after_DP, encrypted_after_LT, newevaluator);
		recrypt_linear_transform(encrypted_after_DP, encrypted_after_LT);

#ifdef _DEBUG
		for (int i = 0; i < extensionDegree_; i++) {
			Plaintext plain_after_LT = decryptor_->rns_decrypt(encrypted_after_LT[i]);
			vector<ExRingElement> slots_after_LT = ExPolyCRTBuilder_->decompose(plain_after_LT);
			cout << "Message after LT : ";
			for (int j = 0; j < ExPolyCRTBuilder_->slot_count(); j++)
			{
				cout << slots_after_LT[j].to_string() << ", ";
			}
			cout << endl;
			cout << "Noise budget after LT : " << decryptor_->invariant_noise_budget(encrypted_after_LT[i]) << endl;
		}
#endif

		//Shift before DigitExtraction
		uint64_t shift = power(p_, e_ - r_) / 2;
		Plaintext const_p2kover2(BigPoly(1, newparms_.plain_modulus().bit_count(), &shift));
		
		//DigitExtraction
		vector<Ciphertext> encrypted_after_DE(extensionDegree_);

		auto time_DE_start = chrono::high_resolution_clock::now();
		for (int i = 0; i < extensionDegree_; i++) {
			newevaluator_->add_plain(encrypted_after_LT[i], const_p2kover2, encrypted_after_LT[i]);
			DigitExtraction(encrypted_after_LT[i], encrypted_after_DE[i]);
			cout << i + 1 << "/" << extensionDegree_ << " digit extraction done ..." << endl;
		}
		auto time_DE_end = chrono::high_resolution_clock::now();

		chrono::microseconds time_DE = chrono::duration_cast<chrono::microseconds>(time_DE_end - time_DE_start);
		cout << "Digit Extraction time (each) : " << (double)time_DE.count() / (1000000 * extensionDegree_) << " seconds " << endl;
		cout << "Digit Extraction time (total) : " << (double)time_DE.count() / 1000000 << " seconds " << endl;

#ifdef _DEBUG
		for (int i = 0; i < extensionDegree_; i++) {
			Plaintext plain_after_DE = decryptor_->rns_decrypt(encrypted_after_DE[i]);
			vector<ExRingElement> slots_after_DE = ExPolyCRTBuilder_->decompose(plain_after_DE);
			cout << "Message after DE : ";
			for (int j = 0; j < ExPolyCRTBuilder_->slot_count(); j++)
			{
				cout << slots_after_DE[j].to_string() << ", ";
			}
			cout << endl;
			cout << "Noise budget after DE : " << decryptor_->invariant_noise_budget(encrypted_after_DE[i]) << endl;
		}
#endif

		//InverseLinearTransformaiton
		Ciphertext encrypted_after_invLT;
		auto time_inv_start = chrono::high_resolution_clock::now();
		InverseLinearTransformation(encrypted_after_DE, encrypted_after_invLT);
		auto time_inv_end = chrono::high_resolution_clock::now();
		chrono::microseconds time_inv = chrono::duration_cast<chrono::microseconds>(time_inv_end - time_inv_start);
		cout << "Inverse linear transform time : " << (double)time_inv.count() / 1000000 << " seconds " << endl;

#ifdef _DEBUG
		cout << "Message after Inverse linear transform : " << decryptor_->rns_decrypt(encrypted_after_invLT).to_string() << endl;
		cout << "Noise budget after inverse linear transform : " << decryptor_->invariant_noise_budget(encrypted_after_invLT) << endl;
#endif

		//Modulus Swith to original coefficient
		desitination = encrypted_after_invLT;
		desitination.rns_hash_block_ = parms_.get_hash_block();
		//newevaluator_->modulus_switching(encrypted_after_invLT, parms_.coeff_modulus(), desitination);
	}

	Ciphertext RNSRecryptor::recrypt(Ciphertext &encrypted)
	{
		Ciphertext result;
		recrypt(encrypted, result);
		return result;
	}

	void RNSRecryptor::recrypt_nonBatch(Ciphertext &encrypted, Ciphertext &desitination)
	{
		// Extract parameters
		int coeff_count = newparms_.poly_modulus().coeff_count();
		int coeff_mod_count = newparms_.coeff_mod_array().size();
		int encrypted_size = encrypted.size();
		int encrypted_ptr_increment = coeff_count * coeff_mod_count;

		if (encrypted_size > 2)
		{
			throw invalid_argument("Input of recrypt should be ciphertext with size 2");
		}

		// Copy encrypted
		Pointer encrypted_copy(allocate_poly(encrypted_size * coeff_count, coeff_mod_count, pool_));
		for (int i = 0; i < encrypted_size; i++)
		{
			set_poly_poly(encrypted.get_array().pointer(i), coeff_count, coeff_mod_count, encrypted_copy.get() + i * encrypted_ptr_increment);
		}

		// Compose the encrypted_copy
		for (int i = 0; i < encrypted_size; i++)
		{
			oldevaluator_->rns_compose(encrypted_copy.get() + i * encrypted_ptr_increment);
		}

		// Multiply newparms_.plain_modulus() and divide in the coeff_products_array
		uint64_t new_plain_modulus = newparms_.plain_modulus().value();
		int big_coeff_mod_count = coeff_mod_count + 1;
		Pointer temp(allocate_zero_uint(big_coeff_mod_count, pool_));
		Pointer quotient(allocate_zero_uint(big_coeff_mod_count, pool_));
		Pointer wide_coeff_products_array(allocate_uint(big_coeff_mod_count, pool_));
		Pointer encrypted_after_MS_ptr(allocate_zero_uint(encrypted_size * coeff_count, pool_));
		set_uint_uint(oldcontext_->coeff_modulus().pointer(), coeff_mod_count, big_coeff_mod_count, wide_coeff_products_array.get());

		for (int i = 0; i < encrypted_size; i++)
		{
			const uint64_t *encrypted_coeff_ptr = encrypted_copy.get() + i * encrypted_ptr_increment;
			for (int j = 0; j < coeff_count - 1; j++)
			{
				multiply_uint_uint64(encrypted_coeff_ptr, coeff_mod_count, new_plain_modulus, big_coeff_mod_count, temp.get());
				divide_uint_uint_inplace(temp.get(), wide_coeff_products_array.get(), big_coeff_mod_count, quotient.get(), pool_);
				set_uint_uint(quotient.get(), big_coeff_mod_count, 1, encrypted_after_MS_ptr.get() + i * coeff_count + j);
				encrypted_coeff_ptr += coeff_mod_count;
			}
		}

		//Dot product
		Ciphertext encrypted_after_DP;
		BigPoly ptxt_c0_poly(coeff_count - 1, 1, encrypted_after_MS_ptr.get());
		BigPoly ptxt_c1_poly(coeff_count - 1, 1, encrypted_after_MS_ptr.get() + coeff_count);
		Plaintext ptxt_c0(ptxt_c0_poly);
		Plaintext ptxt_c1(ptxt_c1_poly);
		newevaluator_->multiply_plain(recrypt_key_, ptxt_c1, encrypted_after_DP);
		newevaluator_->add_plain(encrypted_after_DP, ptxt_c0, encrypted_after_DP);

#ifdef _DEBUG
		cout << "Message after DP : " << decryptor_->rns_decrypt(encrypted_after_DP).to_string() << endl;
		cout << "Noise bubget after DP : " << decryptor_->invariant_noise_budget(encrypted_after_DP) << endl;
#endif

		//LinearTransformation
		Ciphertext encrypted_after_LT;
		LinearTransformation_nonBatch(encrypted_after_DP, encrypted_after_LT);

#ifdef _DEBUG
		Plaintext plain_after_LT = decryptor_->rns_decrypt(encrypted_after_LT);
		vector<ExRingElement> slots_after_LT = ExPolyCRTBuilder_->decompose(plain_after_LT);
		cout << "Message after LT : ";
		int Slot_num = ExPolyCRTBuilder_->slot_count();
		for (int i = 0; i < Slot_num; i++)
		{
			cout << slots_after_LT[i].to_string() << ", ";
		}
		cout << endl;
		cout << "Noise bubget after LT : " << decryptor_->invariant_noise_budget(encrypted_after_LT) << endl;
#endif

		//Shift before DigitExtraction
		uint64_t shift = power(p_, e_ - r_) / 2;
		Plaintext const_p2kover2(BigPoly(1, newparms_.plain_modulus().bit_count(), &shift));
		newevaluator_->add_plain(encrypted_after_LT, const_p2kover2, encrypted_after_LT);

		//DigitExtraction
		Ciphertext encrypted_after_DE;
		DigitExtraction(encrypted_after_LT, encrypted_after_DE);

#ifdef _DEBUG
		Plaintext plain_after_DE = decryptor_->rns_decrypt(encrypted_after_DE);
		vector<ExRingElement> slots_after_DE = ExPolyCRTBuilder_->decompose(plain_after_DE);
		cout << "Message after DE : ";
		for (int i = 0; i < ExPolyCRTBuilder_->slot_count(); i++)
		{
			cout << slots_after_DE[i].to_string() << ", ";
		}
		cout << endl;
		cout << "Noise bubget after DE : " << decryptor_->invariant_noise_budget(encrypted_after_DE) << endl;
#endif

		//InverseLinearTransformaiton
		Ciphertext encrypted_after_invLT;
		InverseLinearTransformation_nonBatch(encrypted_after_DE, encrypted_after_invLT);

#ifdef _DEBUG
		Plaintext plain_after_invLT = decryptor_->rns_decrypt(encrypted_after_invLT);
		vector<ExRingElement> slots_after_invLT = ExPolyCRTBuilder_->decompose(plain_after_invLT);
		cout << "Message after invLT : ";
		for (int i = 0; i < ExPolyCRTBuilder_->slot_count(); i++)
		{
			cout << slots_after_invLT[i].to_string() << ", ";
		}
		cout << endl;
		cout << "Noise bubget after invLT : " << decryptor_->invariant_noise_budget(encrypted_after_invLT) << endl;
#endif

		//Modulus Swith to original coefficient
		desitination = encrypted_after_invLT;
		desitination.rns_hash_block_ = parms_.get_hash_block();
		//newevaluator_->modulus_switching(encrypted_after_invLT, parms_.coeff_modulus(), desitination);
	}

	Ciphertext RNSRecryptor::recrypt_nonBatch(Ciphertext &encrypted)
	{
		Ciphertext result;
		recrypt_nonBatch(encrypted, result);
		return result;
	}

	void RNSRecryptor::floor(Ciphertext &encrypted, Ciphertext &destination, bool invLT)
	{
		int coeff_count = parms_.poly_modulus().coeff_count();
		int coeff_mod_count = parms_.coeff_mod_array().size();
		int coeff_bit_count = coeff_mod_count * bits_per_uint64;
		int coeff_uint64_count = coeff_mod_count;
		int encrypted_ptr_increment = coeff_count * coeff_mod_count;
		int encrypted_count = encrypted.size();

		// Copy encrypted
		Pointer encrypted_copy(allocate_poly(encrypted_count * coeff_count, coeff_mod_count, pool_));
		for (int i = 0; i < encrypted_count; i++)
		{
			set_poly_poly(encrypted.get_array().pointer(i), coeff_count, coeff_mod_count, encrypted_copy.get() + i * encrypted_ptr_increment);
		}

		// Compose the encrypted_copy
		for (int i = 0; i < encrypted_count; i++)
		{
			oldevaluator_->rns_compose(encrypted_copy.get() + i * encrypted_ptr_increment);
		}

		// Coefficient wise division to p^{r-1}
		Pointer encrypted_floor_ptr(allocate_uint(encrypted_count * coeff_count * coeff_mod_count, pool_));

		// Move power_r_minus_one to wide pointer
		BigUInt power_r_minus_one;
		Pointer wide_power_of_r_minus_one(allocate_zero_uint(coeff_mod_count, pool_));
		power_r_minus_one = power(p_, r_ - 1);
		set_uint_uint(power_r_minus_one.pointer(), power_r_minus_one.uint64_count(), coeff_mod_count, wide_power_of_r_minus_one.get());

		// Pointer for quotien and remainder
		Pointer quotient(allocate_zero_uint(coeff_mod_count, pool_));
		Pointer remainder(allocate_zero_uint(coeff_mod_count, pool_));
		Pointer temp(allocate_zero_uint(coeff_mod_count, pool_));

		// Calculate coeff_div_plain_modulus_
		Pointer coeff_div_plain_modulus_(allocate_uint(coeff_mod_count, pool_));
		Pointer wide_plain_modulus(allocate_uint(coeff_mod_count, pool_));
		set_uint_uint(parms_.plain_modulus().pointer(), parms_.plain_modulus().uint64_count(), coeff_mod_count, wide_plain_modulus.get());
		divide_uint_uint(oldcontext_->coeff_modulus().pointer(), wide_plain_modulus.get(), coeff_mod_count, coeff_div_plain_modulus_.get(), remainder.get(), pool_);

		// Calculate coeff_div_plain_modulus / 2
		uint64_t two = 2;
		Pointer delta_div_two(allocate_zero_uint(coeff_mod_count, pool_));
		Pointer wide_two(allocate_zero_uint(coeff_mod_count, pool_));
		set_uint_uint(&two, 1, coeff_mod_count, wide_two.get());
		divide_uint_uint(coeff_div_plain_modulus_.get(), wide_two.get(), coeff_mod_count, delta_div_two.get(), remainder.get(), pool_);

		// Calculate power_r_minus_one * coeff_div_plain_modulus
		Pointer power_r_minus_one_time_delta(allocate_zero_uint(coeff_mod_count, pool_));
		multiply_uint_uint(coeff_div_plain_modulus_.get(), wide_power_of_r_minus_one.get(), coeff_mod_count, power_r_minus_one_time_delta.get());

		for (int i = 0; i < encrypted_count; i++)
		{
			uint64_t *encrypted_coeff_ptr = encrypted_copy.get() + i * encrypted_ptr_increment;
			uint64_t *encrypted_floor_coeff_ptr = encrypted_floor_ptr.get() + i * encrypted_ptr_increment;
			for (int j = 0; j < coeff_count - 1; j++)
			{
				divide_uint_uint(encrypted_coeff_ptr, wide_power_of_r_minus_one.get(), coeff_mod_count, quotient.get(), remainder.get(), pool_);
				if (i == 0)
				{
					add_uint_uint_mod(quotient.get(), power_r_minus_one_time_delta.get(), oldcontext_->coeff_modulus().pointer(), coeff_mod_count, temp.get());
					sub_uint_uint_mod(temp.get(), delta_div_two.get(), oldcontext_->coeff_modulus().pointer(), coeff_mod_count, encrypted_floor_coeff_ptr);
				}
				else
				{
					set_uint_uint(quotient.get(), coeff_uint64_count, coeff_mod_count, encrypted_floor_coeff_ptr);
				}
				encrypted_coeff_ptr += coeff_mod_count;
				encrypted_floor_coeff_ptr += coeff_mod_count;
			}
			// Set zero for last term
			set_zero_uint(coeff_mod_count, encrypted_floor_coeff_ptr);
		}

		// Multiply newparms_.plain_modulus() and divide in the coeff_products_array
		uint64_t new_plain_modulus = newparms_.plain_modulus().value();
		int big_coeff_mod_count = coeff_mod_count + 1;
		Pointer temp2(allocate_zero_uint(big_coeff_mod_count, pool_));
		Pointer quotient2(allocate_zero_uint(big_coeff_mod_count, pool_));
		Pointer wide_coeff_products_array(allocate_uint(big_coeff_mod_count, pool_));
		Pointer encrypted_after_MS_ptr(allocate_zero_uint(encrypted_count * coeff_count, pool_));
		set_uint_uint(oldcontext_->coeff_modulus().pointer(), coeff_mod_count, big_coeff_mod_count, wide_coeff_products_array.get());

		for (int i = 0; i < encrypted_count; i++)
		{
			const uint64_t *encrypted_coeff_ptr = encrypted_floor_ptr.get() + i * encrypted_ptr_increment;
			for (int j = 0; j < coeff_count - 1; j++)
			{
				multiply_uint_uint64(encrypted_coeff_ptr, coeff_mod_count, new_plain_modulus, big_coeff_mod_count, temp2.get());
				divide_uint_uint_inplace(temp2.get(), wide_coeff_products_array.get(), big_coeff_mod_count, quotient2.get(), pool_);
				set_uint_uint(quotient2.get(), big_coeff_mod_count, 1, encrypted_after_MS_ptr.get() + i * coeff_count + j);
				encrypted_coeff_ptr += coeff_mod_count;
			}
		}

		//Dot product
		Ciphertext encrypted_after_DP;
		BigPoly ptxt_c0_poly(coeff_count - 1, 1, encrypted_after_MS_ptr.get());
		BigPoly ptxt_c1_poly(coeff_count - 1, 1, encrypted_after_MS_ptr.get() + coeff_count);
		Plaintext ptxt_c0(ptxt_c0_poly);
		Plaintext ptxt_c1(ptxt_c1_poly);
		newevaluator_->multiply_plain(recrypt_key_, ptxt_c1, encrypted_after_DP);
		newevaluator_->add_plain(encrypted_after_DP, ptxt_c0, encrypted_after_DP);

#ifdef _DEBUG
		cout << "Noise bubget after DP : " << decryptor_->invariant_noise_budget(encrypted_after_DP) << endl;
#endif

		// LinearTransformation
		Ciphertext encrypted_after_LT;
		LinearTransformation_nonBatch(encrypted_after_DP, encrypted_after_LT);

#ifdef _DEBUG
		Plaintext ptxt_after_LT = decryptor_->rns_decrypt(encrypted_after_LT);
		vector<ExRingElement> slots_after_LT = ExPolyCRTBuilder_->decompose(ptxt_after_LT);
		cout << "Message after LT : ";
		for (int i = 0; i < ExPolyCRTBuilder_->slot_count(); i++)
		{
			cout << slots_after_LT[i].to_string() << ", ";
		}
		cout << endl;
		cout << "Noise bubget after LT : " << decryptor_->invariant_noise_budget(encrypted_after_LT) << endl;
#endif

		//Shift before DigitExtraction
		uint64_t shift = power(p_, e_ - r_) / 2;
		Plaintext const_p2kover2(BigPoly(1, newparms_.plain_modulus().bit_count(), &shift));
		newevaluator_->add_plain(encrypted_after_LT, const_p2kover2, encrypted_after_LT);

		// Digit Extraction
		vector<Ciphertext> encrypted_after_DE(e_ - r_ + 1);
		// Special Case: use only 1 lift and 1 removelsd
		if (e_ - r_ == 1)
		{
			// Set eval_plain_modulus for Polynomial Evaluation
			uint64_t plain_modulus_special = power(p_, e_);
			SmallModulus eval_plain_modulus_special = power(p_, e_);

			// Calculate lift_poly
			vector<uint64_t> lift_poly_special;
			lift_poly_special = compute_lift_uint64(p_, 1);

			// Calculate remainlsd_poly
			vector<uint64_t> remainlsd_poly_special;
			remainlsd_poly_special = compute_remainlsd_uint64(p_, r_);

			// Lift in encrypted state
			Ciphertext encrypted_after_lift = newevaluator_->Polynomial_Evaluation(encrypted_after_LT, lift_poly_special, eval_plain_modulus_special);

#ifdef _DEBUG
			cout << "Noise bubget after lift : " << decryptor_->invariant_noise_budget(encrypted_after_lift) << endl;
#endif

			// Substract from encrypted_after_LT (least significant digit is removed!!!)
			Ciphertext encrypted_origin_minus_lift = newevaluator_->sub(encrypted_after_LT, encrypted_after_lift);

			// Change the eval_plain_modulus
			plain_modulus_special /= p_;
			eval_plain_modulus_special = plain_modulus_special;

			// Evaluate remainlsd_poly for remove errors above the message
			encrypted_after_DE[e_ - r_] = newevaluator_->Polynomial_Evaluation(encrypted_origin_minus_lift, remainlsd_poly_special, eval_plain_modulus_special);
		}
		else
		{
			// In general just use DigitiExtraction Algorithm in newevaluator
			encrypted_after_DE = newevaluator_->DigitExtraction(encrypted_after_LT, p_, e_, e_ - r_ + 1);
		}

#ifdef _DEBUG
		Plaintext ptxt_after_DE = decryptor_->rns_decrypt(encrypted_after_DE[e_ - r_]);
		vector<ExRingElement> slots_after_DE = ExPolyCRTBuilder_->decompose(ptxt_after_DE);
		cout << "Message after DE : ";
		for (int i = 0; i < ExPolyCRTBuilder_->slot_count(); i++)
		{
			cout << slots_after_DE[i].to_string() << ", ";
		}
		cout << endl;
		cout << "Noise bubget after DE : " << decryptor_->invariant_noise_budget(encrypted_after_DE[e_ - r_]) << endl;
#endif

		// If user does not want to do invLT at last step go to this branch
		if (invLT == false)
		{
			destination = encrypted_after_DE[e_ - r_];
			destination.rns_hash_block_ = parms_.get_hash_block();
			return;
		}

		//InverseLinearTransformaiton
		Ciphertext encrypted_after_invLT;
		InverseLinearTransformation_nonBatch(encrypted_after_DE[e_ - r_], encrypted_after_invLT);

#ifdef _DEBUG
		Plaintext ptxt_after_invLT = decryptor_->rns_decrypt(encrypted_after_invLT);
		vector<ExRingElement> slots_after_invLT = ExPolyCRTBuilder_->decompose(ptxt_after_invLT);
		cout << "Message after invLT : ";
		cout << ptxt_after_invLT.to_string() << endl;
		cout << "Slots after invLT : ";
		for (int i = 0; i < ExPolyCRTBuilder_->slot_count(); i++)
		{
			cout << slots_after_invLT[i].to_string() << ", ";
		}
		cout << endl;
		cout << "Noise bubget after invLT : " << decryptor_->invariant_noise_budget(encrypted_after_invLT) << endl;
#endif

		//Modulus Switch to original coefficient
		if (newparms_.coeff_mod_array() == parms_.coeff_mod_array())
		{
			destination = encrypted_after_invLT;
		}
		else
		{
			//newevaluator_->modulus_switching(encrypted_after_invLT, parms_.coeff_modulus(), destination);
		}
		destination.rns_hash_block_ = parms_.get_hash_block();
	}

	Ciphertext RNSRecryptor::floor(Ciphertext &encrypted, bool invLT)
	{
		Ciphertext result;
		floor(encrypted, result, invLT);
		return result;
	}

	void RNSRecryptor::floor_batch(Ciphertext &encrypted, vector<Ciphertext> &destination, int num_Batch, int num_thread)
	{
		int coeff_count = parms_.poly_modulus().coeff_count();
		int coeff_mod_count = parms_.coeff_mod_array().size();
		int coeff_bit_count = coeff_mod_count * bits_per_uint64;
		int coeff_uint64_count = coeff_mod_count;
		int encrypted_ptr_increment = coeff_count * coeff_mod_count;
		int encrypted_count = encrypted.size();

		// Copy encrypted
		Pointer encrypted_copy(allocate_poly(encrypted_count * coeff_count, coeff_mod_count, pool_));
		for (int i = 0; i < encrypted_count; i++)
		{
			set_poly_poly(encrypted.get_array().pointer(i), coeff_count, coeff_mod_count, encrypted_copy.get() + i * encrypted_ptr_increment);
		}

		// Compose the encrypted_copy
		for (int i = 0; i < encrypted_count; i++)
		{
			oldevaluator_->rns_compose(encrypted_copy.get() + i * encrypted_ptr_increment);
		}

		// Coefficient wise division to p^{r-1}
		Pointer encrypted_floor_ptr(allocate_uint(encrypted_count * coeff_count * coeff_mod_count, pool_));

		// Move power_r_minus_one to wide pointer
		BigUInt power_r_minus_one;
		Pointer wide_power_of_r_minus_one(allocate_zero_uint(coeff_mod_count, pool_));
		power_r_minus_one = power(p_, r_ - 1);
		set_uint_uint(power_r_minus_one.pointer(), power_r_minus_one.uint64_count(), coeff_mod_count, wide_power_of_r_minus_one.get());

		// Pointer for quotien and remainder
		Pointer quotient(allocate_zero_uint(coeff_mod_count, pool_));
		Pointer remainder(allocate_zero_uint(coeff_mod_count, pool_));
		Pointer temp(allocate_zero_uint(coeff_mod_count, pool_));

		// Calculate coeff_div_plain_modulus_
		Pointer coeff_div_plain_modulus_(allocate_uint(coeff_mod_count, pool_));
		Pointer wide_plain_modulus(allocate_uint(coeff_mod_count, pool_));
		set_uint_uint(parms_.plain_modulus().pointer(), parms_.plain_modulus().uint64_count(), coeff_mod_count, wide_plain_modulus.get());
		divide_uint_uint(oldcontext_->coeff_modulus().pointer(), wide_plain_modulus.get(), coeff_mod_count, coeff_div_plain_modulus_.get(), remainder.get(), pool_);

		// Calculate coeff_div_plain_modulus / 2
		uint64_t two = 2;
		Pointer delta_div_two(allocate_zero_uint(coeff_mod_count, pool_));
		Pointer wide_two(allocate_zero_uint(coeff_mod_count, pool_));
		set_uint_uint(&two, 1, coeff_mod_count, wide_two.get());
		divide_uint_uint(coeff_div_plain_modulus_.get(), wide_two.get(), coeff_mod_count, delta_div_two.get(), remainder.get(), pool_);

		// Calculate power_r_minus_one * coeff_div_plain_modulus
		Pointer power_r_minus_one_time_delta(allocate_zero_uint(coeff_mod_count, pool_));
		multiply_uint_uint(coeff_div_plain_modulus_.get(), wide_power_of_r_minus_one.get(), coeff_mod_count, power_r_minus_one_time_delta.get());

		for (int i = 0; i < encrypted_count; i++)
		{
			uint64_t *encrypted_coeff_ptr = encrypted_copy.get() + i * encrypted_ptr_increment;
			uint64_t *encrypted_floor_coeff_ptr = encrypted_floor_ptr.get() + i * encrypted_ptr_increment;
			for (int j = 0; j < coeff_count - 1; j++)
			{
				divide_uint_uint(encrypted_coeff_ptr, wide_power_of_r_minus_one.get(), coeff_mod_count, quotient.get(), remainder.get(), pool_);
				if (i == 0)
				{
					sub_uint_uint_mod(quotient.get(), delta_div_two.get(), oldcontext_->coeff_modulus().pointer(), coeff_mod_count, encrypted_floor_coeff_ptr);
				}
				else
				{
					set_uint_uint(quotient.get(), coeff_mod_count, coeff_mod_count, encrypted_floor_coeff_ptr);
				}
				encrypted_coeff_ptr += coeff_mod_count;
				encrypted_floor_coeff_ptr += coeff_mod_count;
			}
			// Set zero for last term
			set_zero_uint(coeff_mod_count, encrypted_floor_coeff_ptr);
		}

		// Multiply newparms_.plain_modulus() and divide in the coeff_products_array
		uint64_t new_plain_modulus = newparms_.plain_modulus().value();
		int big_coeff_mod_count = coeff_mod_count + 1;
		Pointer temp2(allocate_zero_uint(big_coeff_mod_count, pool_));
		Pointer quotient2(allocate_zero_uint(big_coeff_mod_count, pool_));
		Pointer wide_coeff_products_array(allocate_uint(big_coeff_mod_count, pool_));
		Pointer encrypted_after_MS_ptr(allocate_zero_uint(encrypted_count * coeff_count, pool_));
		set_uint_uint(oldcontext_->coeff_modulus().pointer(), coeff_mod_count, big_coeff_mod_count, wide_coeff_products_array.get());

		for (int i = 0; i < encrypted_count; i++)
		{
			const uint64_t *encrypted_coeff_ptr = encrypted_floor_ptr.get() + i * encrypted_ptr_increment;
			for (int j = 0; j < coeff_count - 1; j++)
			{
				multiply_uint_uint64(encrypted_coeff_ptr, coeff_mod_count, new_plain_modulus, big_coeff_mod_count, temp2.get());
				divide_uint_uint_inplace(temp2.get(), wide_coeff_products_array.get(), big_coeff_mod_count, quotient2.get(), pool_);
				set_uint_uint(quotient2.get(), big_coeff_mod_count, 1, encrypted_after_MS_ptr.get() + i * coeff_count + j);
				encrypted_coeff_ptr += coeff_mod_count;
			}
		}

		// Dot product
		Ciphertext encrypted_after_DP;
		BigPoly ptxt_c0_poly(coeff_count, 1, encrypted_after_MS_ptr.get());
		BigPoly ptxt_c1_poly(coeff_count, 1, encrypted_after_MS_ptr.get() + coeff_count);
		Plaintext ptxt_c0(ptxt_c0_poly);
		Plaintext ptxt_c1(ptxt_c1_poly);
		newevaluator_->multiply_plain(recrypt_key_, ptxt_c1, encrypted_after_DP);
		newevaluator_->add_plain(encrypted_after_DP, ptxt_c0, encrypted_after_DP);

#ifdef _DEBUG
		cout << "Noise bubget after DP : " << decryptor_->invariant_noise_budget(encrypted_after_DP) << endl;
#endif

		// LinearTransformation
		Ciphertext encrypted_after_LT;
		LinearTransformation_nonBatch(encrypted_after_DP, encrypted_after_LT, num_Batch, num_thread);

#ifdef _DEBUG
		Plaintext ptxt_after_LT = decryptor_->rns_decrypt(encrypted_after_LT);
		vector<ExRingElement> slots_after_LT = ExPolyCRTBuilder_->decompose(ptxt_after_LT);
		cout << "Message after LT : ";
		for (int i = 0; i < ExPolyCRTBuilder_->slot_count(); i++)
		{
			cout << slots_after_LT[i].to_string() << ", ";
		}
		cout << endl;
		cout << "Noise bubget after LT : " << decryptor_->invariant_noise_budget(encrypted_after_LT) << endl;
#endif

		auto time_DE_start = chrono::high_resolution_clock::now();
		// Shift before DigitExtraction
		uint64_t shift = power(p_, e_ - r_) / 2;
		Plaintext const_p2kover2(BigPoly(1, newparms_.plain_modulus().bit_count(), &shift));
		newevaluator_->add_plain(encrypted_after_LT, const_p2kover2, encrypted_after_LT);

		// Digit Extraction
		vector<Ciphertext> encrypted_after_DE(e_ - r_ + 1);
		// Special Case: use only 1 lift and 1 removelsd
		if (e_ - r_ == 1)
		{
			// Set eval_plain_modulus for Polynomial Evaluation
			uint64_t plain_modulus_special = power(p_, e_);
			SmallModulus eval_plain_modulus_special = power(p_, e_);

			// Calculate lift_poly
			vector<uint64_t> lift_poly_special;
			lift_poly_special = compute_lift_uint64(p_, 1);

			// Calculate remainlsd_poly
			vector<uint64_t> remainlsd_poly_special;
			remainlsd_poly_special = compute_remainlsd_uint64(p_, r_);

			// Lift in encrypted state
			Ciphertext encrypted_after_lift = newevaluator_->Polynomial_Evaluation(encrypted_after_LT, lift_poly_special, eval_plain_modulus_special);

#ifdef _DEBUG
			cout << "Noise bubget after lift : " << decryptor_->invariant_noise_budget(encrypted_after_lift) << endl;
#endif

			// Substract from encrypted_after_LT (least significant digit is removed!!!)
			Ciphertext encrypted_origin_minus_lift = newevaluator_->sub(encrypted_after_LT, encrypted_after_lift);

			// Change the eval_plain_modulus
			plain_modulus_special /= p_;
			eval_plain_modulus_special = plain_modulus_special;

			// Evaluate remainlsd_poly for remove errors above the message
			encrypted_after_DE[e_ - r_] = newevaluator_->Polynomial_Evaluation(encrypted_origin_minus_lift, remainlsd_poly_special, eval_plain_modulus_special);
		}
		else
		{
			// In general just use DigitiExtraction Algorithm in newevaluator
			encrypted_after_DE = newevaluator_->DigitExtraction(encrypted_after_LT, p_, e_, e_ - r_ + 1);
		}
		auto time_DE_end = chrono::high_resolution_clock::now();
		cout << "... DE time : " << (double)chrono::duration_cast<chrono::microseconds>(time_DE_end - time_DE_start).count() / 1000000. << " seconds" << endl;

#ifdef _DEBUG
		Plaintext ptxt_after_DE = decryptor_->rns_decrypt(encrypted_after_DE[e_ - r_]);
		vector<ExRingElement> slots_after_DE = ExPolyCRTBuilder_->decompose(ptxt_after_DE);
		cout << "Message after DE : ";
		for (int i = 0; i < ExPolyCRTBuilder_->slot_count(); i++)
		{
			cout << slots_after_DE[i].to_string() << ", ";
		}
		cout << endl;
		cout << "Noise bubget after DE : " << decryptor_->invariant_noise_budget(encrypted_after_DE[e_ - r_]) << endl;
#endif

		// Make a unit plaintext which is generalizaed batching of unit vector
		int num_Slot = ExPolyCRTBuilder_->slot_count();
		int log_num_slot = get_power_of_two(num_Slot);

		auto time_sp_start = chrono::high_resolution_clock::now();
		atomic<int> index = ATOMIC_VAR_INIT(0);
		auto mthread = [&]()
		{
			while (true)
			{
				int my_index = index++;
				if (my_index >= num_Batch)
					break;
				// Functions for each my_index;
				destination[my_index] = newevaluator_->multiply_plain(encrypted_after_DE[e_ - r_], unit_plain_[my_index]);
				for (int i = 0; i < log_num_slot; i++)
				{
					Ciphertext rotate;
					rotate = newevaluator_->apply_galois(destination[my_index], power_of_three_coset_index[i]);
					destination[my_index] = newevaluator_->add(destination[my_index], rotate);
				}
				destination[my_index].rns_hash_block_ = parms_.get_hash_block();
			}
		};

		destination.resize(num_Batch);
		vector<thread> pool_split;
		for (int i = 0; i < num_thread; i++)
		{
			pool_split.emplace_back(mthread);
		}
		for (int i = 0; i < num_thread; i++)
		{
			pool_split[i].join();
		}

		auto time_sp_end = chrono::high_resolution_clock::now();
		cout << "... split time : " << (double)chrono::duration_cast<chrono::microseconds>(time_sp_end - time_sp_start).count() / 1000000. << " seconds" << endl;
	}

	vector<Ciphertext> RNSRecryptor::floor_batch(Ciphertext &encrypted, int num_Batch, int num_thread)
	{
		vector<Ciphertext> result;
		floor_batch(encrypted, result, num_Batch, num_thread);
		return result;
	}


	// Send encryption of (Delta1, .., Delta1), ..., (Deltak, ..., Deltak) to 
	// (Delta1, .., Deltak ,0, ..., 0). 
	void RNSRecryptor::combine(vector<Ciphertext> &encrypted, Ciphertext & destination) {
		
		int coeff_count = parms_.poly_modulus().coeff_count();
		int coeff_mod_count = parms_.coeff_mod_array().size();
		int coeff_bit_count = coeff_mod_count * bits_per_uint64;
		int coeff_uint64_count = coeff_mod_count;
		int encrypted_ptr_increment = coeff_count * coeff_mod_count;
		int encrypted_count = encrypted.size();


		// Multiply old unit vector
		destination = oldevaluator_->multiply_plain(encrypted[0], old_unit_plain_[0]);
		for (int i = 0; i < encrypted.size(); i++)
		{
			Ciphertext temp = oldevaluator_->multiply_plain(encrypted[i], old_unit_plain_[i]);
			oldevaluator_->add(destination, temp, destination);
		}
		return;
	}

	// This is the inverse map of RNSRecryptor::combine
	// Send encryption of (Delta1, ..., Deltak, 0,..,0) to 
	// (Delta1, ..., Delta1), ..., (Deltak ,..., Deltak).
	void RNSRecryptor::expand(const Ciphertext &encrypted, vector<Ciphertext> & destination) 
	{
		int coeff_count = parms_.poly_modulus().coeff_count();
		int coeff_mod_count = parms_.coeff_mod_array().size();
		int coeff_bit_count = coeff_mod_count * bits_per_uint64;
		int coeff_uint64_count = coeff_mod_count;
		int encrypted_ptr_increment = coeff_count * coeff_mod_count;
		int encrypted_count = encrypted.size();

		for (int i = 0; i < destination.size(); i++) {
			// Multiply old unit vector to get (m1, 0, ..., 0).
			Ciphertext temp = oldevaluator_->multiply_plain(encrypted, old_unit_plain_[i]);

			// log(numSlots_) rotations to send it to (m1, m1, ...)
			int lognumSlot = get_power_of_two(oldExPolyCRTBuilder_->slot_count());
			for (int j = 0; j < lognumSlot; j++)
			{
				Ciphertext rotate;
				rotate = oldevaluator_->apply_galois(temp, power_of_three_coset_index[j]);
				temp = oldevaluator_->add(temp, rotate);
			}
			destination[i] = temp;
		}
		return; 
	}


	void RNSRecryptor::floor_simd(Ciphertext &encrypted, Ciphertext &destination, int k, int num_thread)
	{
		if (k >= parms_.plain_modulus_exponent())
		{
			throw invalid_argument("k should be smaller than plain modulus exponent");
		}

		vector<thread> pool;

		int coeff_count = parms_.poly_modulus().coeff_count();
		int coeff_mod_count = parms_.coeff_mod_array().size();
		int coeff_bit_count = coeff_mod_count * bits_per_uint64;
		int coeff_uint64_count = coeff_mod_count;
		int encrypted_ptr_increment = coeff_count * coeff_mod_count;
		int encrypted_count = encrypted.size();

		// Step1. Send Encryption of (m1, m2, m3, ...) to (m1 + m2 * X + m3 * X^2 + ...)

		// Multiply old unit vector
		Ciphertext encrypted_after_invLT;
		encrypted_after_invLT = oldevaluator_->multiply_plain(encrypted, old_unit_plain_[0]);

		// log(numSlots_) rotations to send it to (m1, m1, ...)
		int lognumSlot = get_power_of_two(oldExPolyCRTBuilder_->slot_count());
		for (int i = 0; i < lognumSlot; i++)
		{
			Ciphertext rotate;
			rotate = oldevaluator_->apply_galois(encrypted_after_invLT, power_of_three_coset_index[i]);
			encrypted_after_invLT = oldevaluator_->add(encrypted_after_invLT, rotate);
		}

		vector<Ciphertext> encrypted_ith_slot(numSlots_);
		atomic<int> index_invLT = ATOMIC_VAR_INIT(0);
		auto mthread_invLT = [&]()
		{
			while (true)
			{
				int myindex = index_invLT++;
				if (myindex >= numSlots_)
				{
					break;
				}
				encrypted_ith_slot[myindex] = oldevaluator_->multiply_plain(encrypted, old_unit_plain_[myindex]);
				for (int i = 0; i < lognumSlot; i++)
				{
					Ciphertext rotate;
					rotate = oldevaluator_->apply_galois(encrypted_ith_slot[myindex], power_of_three_coset_index[i]);
					encrypted_ith_slot[myindex] = oldevaluator_->add(encrypted_ith_slot[myindex], rotate);
				}
				string temp_string = "1x^" + to_string(myindex);
				Plaintext power_of_X = temp_string;
				encrypted_ith_slot[myindex] = oldevaluator_->multiply_plain(encrypted_ith_slot[myindex], power_of_X);
			}
		};

		for (int i = 0; i < num_thread; i++)
		{
			pool.emplace_back(mthread_invLT);
		}

		for (int i = 0; i < num_thread; i++)
		{
			pool[i].join();
		}
		pool.clear();

		// Add all encrypted_ith_slot
		encrypted_after_invLT = encrypted_ith_slot[0];
		for (int i = 1; i < numSlots_; i++)
		{
			encrypted_after_invLT = oldevaluator_->add(encrypted_after_invLT, encrypted_ith_slot[i]);
		}

		// Copy encrypted
		Pointer encrypted_copy(allocate_poly(encrypted_count * coeff_count, coeff_mod_count, pool_));
		for (int i = 0; i < encrypted_count; i++)
		{
			set_poly_poly(encrypted_after_invLT.get_array().pointer(i), coeff_count, coeff_mod_count, encrypted_copy.get() + i * encrypted_ptr_increment);
		}

		// Compose the encrypted_copy
		for (int i = 0; i < encrypted_count; i++)
		{
			oldevaluator_->rns_compose(encrypted_copy.get() + i * encrypted_ptr_increment);
		}

		// Step2. Divide each coeffcient to p^k

		Pointer encrypted_floor_ptr(allocate_zero_uint(encrypted_count * coeff_count * coeff_mod_count, pool_));

		// Move power_r_minus_one to wide pointer
		BigUInt power_k;
		Pointer wide_power_k(allocate_zero_uint(coeff_mod_count, pool_));
		power_k = power(p_, k);
		set_uint_uint(power_k.pointer(), power_k.uint64_count(), coeff_mod_count, wide_power_k.get());

		// Pointer for quotien and remainder
		Pointer quotient(allocate_zero_uint(coeff_mod_count, pool_));
		Pointer remainder(allocate_zero_uint(coeff_mod_count, pool_));
		Pointer temp(allocate_zero_uint(coeff_mod_count, pool_));

		// Calculate coeff_div_plain_modulus_
		Pointer coeff_div_plain_modulus_(allocate_uint(coeff_mod_count, pool_));
		Pointer wide_plain_modulus(allocate_uint(coeff_mod_count, pool_));
		set_uint_uint(parms_.plain_modulus().pointer(), parms_.plain_modulus().uint64_count(), coeff_mod_count, wide_plain_modulus.get());
		divide_uint_uint(oldcontext_->coeff_modulus().pointer(), wide_plain_modulus.get(), coeff_mod_count, coeff_div_plain_modulus_.get(), remainder.get(), pool_);

		// Calculate coeff_div_plain_modulus / 2
		uint64_t two = 2;
		Pointer delta_div_two(allocate_zero_uint(coeff_mod_count, pool_));
		Pointer wide_two(allocate_zero_uint(coeff_mod_count, pool_));
		set_uint_uint(&two, 1, coeff_mod_count, wide_two.get());
		divide_uint_uint(coeff_div_plain_modulus_.get(), wide_two.get(), coeff_mod_count, delta_div_two.get(), remainder.get(), pool_);

		// Calculate power_k * coeff_div_plain_modulus
		Pointer power_r_minus_one_time_delta(allocate_zero_uint(coeff_mod_count, pool_));
		multiply_uint_uint(coeff_div_plain_modulus_.get(), wide_power_k.get(), coeff_mod_count, power_r_minus_one_time_delta.get());

		for (int i = 0; i < encrypted_count; i++)
		{
			uint64_t *encrypted_coeff_ptr = encrypted_copy.get() + i * encrypted_ptr_increment;
			uint64_t *encrypted_floor_coeff_ptr = encrypted_floor_ptr.get() + i * encrypted_ptr_increment;
			for (int j = 0; j < coeff_count - 1; j++)
			{
				divide_uint_uint(encrypted_coeff_ptr, wide_power_k.get(), coeff_mod_count, quotient.get(), remainder.get(), pool_);
				if (i == 0)
				{
					sub_uint_uint_mod(quotient.get(), delta_div_two.get(), oldcontext_->coeff_modulus().pointer(), coeff_mod_count, encrypted_floor_coeff_ptr);
				}
				else
				{
					set_uint_uint(quotient.get(), coeff_mod_count, coeff_mod_count, encrypted_floor_coeff_ptr);
				}
				encrypted_coeff_ptr += coeff_mod_count;
				encrypted_floor_coeff_ptr += coeff_mod_count;
			}
			// Set zero for last term
			set_zero_uint(coeff_mod_count, encrypted_floor_coeff_ptr);
		}

		// Multiply newparms_.plain_modulus() and divide in the coeff_products_array
		uint64_t new_plain_modulus = newparms_.plain_modulus().value();
		int big_coeff_mod_count = coeff_mod_count + 1;
		Pointer temp2(allocate_zero_uint(big_coeff_mod_count, pool_));
		Pointer quotient2(allocate_zero_uint(big_coeff_mod_count, pool_));
		Pointer wide_coeff_products_array(allocate_uint(big_coeff_mod_count, pool_));
		Pointer encrypted_after_MS_ptr(allocate_zero_uint(encrypted_count * coeff_count, pool_));
		set_uint_uint(oldcontext_->coeff_modulus().pointer(), coeff_mod_count, big_coeff_mod_count, wide_coeff_products_array.get());

		for (int i = 0; i < encrypted_count; i++)
		{
			const uint64_t *encrypted_coeff_ptr = encrypted_floor_ptr.get() + i * encrypted_ptr_increment;
			for (int j = 0; j < coeff_count - 1; j++)
			{
				multiply_uint_uint64(encrypted_coeff_ptr, coeff_mod_count, new_plain_modulus, big_coeff_mod_count, temp2.get());
				divide_uint_uint_inplace(temp2.get(), wide_coeff_products_array.get(), big_coeff_mod_count, quotient2.get(), pool_);
				set_uint_uint(quotient2.get(), big_coeff_mod_count, 1, encrypted_after_MS_ptr.get() + i * coeff_count + j);
				encrypted_coeff_ptr += coeff_mod_count;
			}
		}

		// Dot product
		Ciphertext encrypted_after_DP;
		BigPoly ptxt_c0_poly(coeff_count, 1, encrypted_after_MS_ptr.get());
		BigPoly ptxt_c1_poly(coeff_count, 1, encrypted_after_MS_ptr.get() + coeff_count);
		Plaintext ptxt_c0(ptxt_c0_poly);
		Plaintext ptxt_c1(ptxt_c1_poly);
		newevaluator_->multiply_plain(recrypt_key_, ptxt_c1, encrypted_after_DP);
		newevaluator_->add_plain(encrypted_after_DP, ptxt_c0, encrypted_after_DP);

		//#ifdef _DEBUG
		//cout << "Message after DP: "; 
		//auto plainDP = decryptor_->rns_decrypt(encrypted_after_DP);
		//uint64_t ptoe = power(p_, e_); 
		//for (int i = 0; i < n_; i++) {
		//	int64_t value = static_cast<int64_t> (plainDP[i].pointer()[0]); 
		//	// if (value > ptoe / 2) value -= ptoe;
		//	cout << i << "-th value = " << value; 
		//	vector<int> digits = centered_digits( value, p_, e_);
		//	cout << ", digits = ["; 
		//	for (int j = 0; j <  digits.size(); j++) {
		//		cout << digits[j] << ", ";
		//	}
		//	cout << "] " << endl;;
		//}
		//cout << endl;
		cout << "Noise bubget after DP : " << decryptor_->invariant_noise_budget(encrypted_after_DP) << endl;
		//#endif

		// LinearTransformation
		Ciphertext encrypted_after_LT;
		LinearTransformation_nonBatch(encrypted_after_DP, encrypted_after_LT, num_thread);

		//#ifdef _DEBUG
		Plaintext ptxt_after_LT = decryptor_->rns_decrypt(encrypted_after_LT);
		vector<ExRingElement> slots_after_LT = ExPolyCRTBuilder_->decompose(ptxt_after_LT);
		//cout << "Message after LT : ";
		//for (int i = 0; i < ExPolyCRTBuilder_->slot_count(); i++)
		//{
		//	int64_t value = static_cast<int64_t> (slots_after_LT[i].pointer()[0]);
		//	// if (value > ptoe / 2) value -= ptoe;
		//	cout << i << "-th value = " << value;
		//	vector<int> digits = centered_digits(value, p_, e_);
		//	cout << ", digits = [";
		//	for (int j = 0; j < digits.size(); j++) {
		//		cout << digits[j] << ", ";
		//	}
		//	cout << "] " << endl;;
		//}
		//cout << endl;
		cout << "Noise bubget after LT : " << decryptor_->invariant_noise_budget(encrypted_after_LT) << endl;
		//#endif

		auto time_DE_start = chrono::high_resolution_clock::now();

		// Shift before DigitExtraction
		//uint64_t shift = power(p_, e_ - r_) / 2;
		//Plaintext const_p2kover2(BigPoly(1, newparms_.plain_modulus().bit_count(), &shift));
		//newevaluator_->add_plain(encrypted_after_LT, const_p2kover2, encrypted_after_LT);

		// Digit Extraction
		vector<Ciphertext> encrypted_after_DE(e_ - k);
		encrypted_after_DE = newevaluator_->DigitExtraction(encrypted_after_LT, p_, e_, e_ - k, true);

		// Debug:

		//if (true)
		//{
		//	for (int i = 0; i < encrypted_after_DE.size(); i++) {
		//		Plaintext ptxt = decryptor_->rns_decrypt(encrypted_after_DE[i]);
		//		vector<ExRingElement> slots = ExPolyCRTBuilder_->decompose(ptxt);
		//		cout << i << "-th digit after DE : ";
		//		for (int rr = 0; rr < ExPolyCRTBuilder_->slot_count(); rr++)
		//		{
		//			int64_t value = static_cast<int64_t> (slots[rr].pointer()[0]);
		//			if (value > ptoe / 2) value -= ptoe;
		//			cout << value << ", ";
		//		}
		//		cout << endl;
		//	}
		//}



		auto time_DE_end = chrono::high_resolution_clock::now();
		cout << "... DE time : " << (double)chrono::duration_cast<chrono::microseconds>(time_DE_end - time_DE_start).count() / 1000000. << " seconds" << endl;

		// Add encrypted_after_DE[e_ - r_] + ... encrypted_after_DE[e_ - k_- 1]
		destination = encrypted_after_DE[e_ - r_];
		for (int i = 1; i < r_ - k; i++)
		{
			cout << "noise budget =  " << decryptor_->invariant_noise_budget(encrypted_after_DE[e_ - r_ + i]) << endl;
			destination = newevaluator_->add(destination, encrypted_after_DE[e_ - r_ + i]);
		}

		// Final debug
		//cout << "intermediate result = ";
		//Plaintext ptxt = decryptor_->rns_decrypt(destination);
		//vector<ExRingElement> slots = ExPolyCRTBuilder_->decompose(ptxt);
		//for (int rr = 0; rr < ExPolyCRTBuilder_->slot_count(); rr++)
		//{
		//	int64_t value = static_cast<int64_t> (slots[rr].pointer()[0]);
		//	if (value > ptoe / 2) value -= ptoe;
		//	cout << value << ", ";
		//}
		//cout << endl;
		// 


		destination.rns_hash_block_ = parms_.get_hash_block();
	}

	Ciphertext RNSRecryptor::floor_simd(Ciphertext &encrypted, int k, int num_thread)
	{
		Ciphertext result;
		floor_simd(encrypted, result, k, num_thread);
		return result;
	}

	void RNSRecryptor::recrypt_pre_linear_transform(const Ciphertext & encrypted, Ciphertext & destination)
	{
		BigUInt plain_modulus_(newparms_.plain_modulus().bit_count(), newparms_.plain_modulus().value());
		//BigPoly g1 = g1_;
		//BigPoly g1squared;
		//BigPolyArith arith;
		BigPoly const_ninverse;
		const_ninverse.resize(newparms_.poly_modulus().coeff_count(), newparms_.plain_modulus().bit_count());
		const_ninverse.set_zero();
		const_ninverse[0] = nInverse_;

		//arith.multiply(g1, g1, newparms_.poly_modulus(), plain_modulus_, g1squared);
		////vector<Plaintext> result_polys;
		//vector<vector<Plaintext>> result_polys(giantstep_.size());
		//for (int i = 0; i < result_polys.size(); i++) {
		//	result_polys[i].resize(babystep_.size());
		//}
		//BigPoly tempg = g1_;
		//BigPoly temp;
		//Modulus plainMod(newparms_.plain_modulus().pointer(), 1, pool_);
		//for (int i = 0; i < n_; i++)
		//{
		//	auto index = decompose_bs_gs(m_, (2 * i + 1) % m_, babystep_, giantstep_);
		//	// first is gs,second is bs.
		//	// result_polys[index.first][index.second] = Plaintext();
		//	//ConstPointer poly1ptr = duplicate_poly_if_needed(poly1, coeff_count, coeff_uint64_count, poly1.pointer() == result.pointer(), pool_);
		//	//ConstPointer poly2ptr = duplicate_poly_if_needed(poly2, coeff_count, coeff_uint64_count, poly2.pointer() == result.pointer(), pool_);
		//	//ConstPointer polymodptr = duplicate_poly_if_needed(poly_mod, coeff_count, coeff_uint64_count, poly_mod.pointer() == result.pointer(), pool_);

		//	// Verify destination size.
		//	result_polys[index.first][index.second].get_poly().resize(n_+1, plainMod.significant_bit_count());
		//	// Multiply polynomials.
		//	//PolyModulus polymod(polymodptr.get(), coeff_count, coeff_uint64_count);
		//	nussbaumer_multiply_poly_poly_coeffmod(tempg.pointer(), pre_linear_coeffs_[i % extensionDegree_].pointer(), get_power_of_two(n_) , plainMod, result_polys[index.first][index.second].get_poly().pointer(), pool_);
		//	// result_polys[index.first][index.second].get_poly().duplicate_from(arith.multiply(tempg, pre_linear_coeffs_[i % extensionDegree_], newparms_.poly_modulus(), plain_modulus_));
		//	//result_polys.push_back(Plaintext(arith.multiply(tempg, pre_linear_coeffs_[i % extensionDegree_], newparms_.poly_modulus(), plain_modulus_)));
		//	if (i < n_ - 1) {
		//		arith.multiply(g1squared, tempg, newparms_.poly_modulus(), plain_modulus_, temp);
		//		// g_polys.push_back(temp);
		//		tempg = temp;
		//	}
		//	// if (i== 0) cout << 2 * i + 1 << " , " << result_polys[index.first][index.second].get_poly().to_string() << endl;
		//}
		// final step
#ifdef _DEBUG
		//for (int i = 0; i < result_polys.size(); i++) {
		//    cout << i << "-th result poly = " << result_polys[i].to_string() << endl;
		//}
#endif
		auto time_pre_start = chrono::high_resolution_clock::now();
		Ciphertext tempdest;
		//newevaluator_->rotations_weighted_sum(encrypted, result_polys, babystep_, giantstep_, tempdest);
		newevaluator_->rotations_weighted_sum(encrypted, prelinear_polys_, babystep_, giantstep_, tempdest, false);
		// Final scaling by ninverse.  
		// Fixme: change this to multiply_plain_uint
		newevaluator_->multiply_plain(tempdest, const_ninverse, destination);
		auto time_pre_end = chrono::high_resolution_clock::now();
		chrono::microseconds time_pre = chrono::duration_cast<chrono::microseconds>(time_pre_end - time_pre_start);
		cout << "Rotations weighted sum total time: " << (double)time_pre.count() / 1000000 << " seconds" << endl;
		return;
	}

	void RNSRecryptor::recrypt_unpack(const Ciphertext & encrypted, vector<Ciphertext>& destination)
	{
		BigPolyArith arith;
		int d = extensionDegree_;
		vector<Ciphertext> encrypted_rotations(d);
		BigPolyArray temp;
		int p = p_;
		int gal = 1;
		int m = m_;
		destination.resize(d);
		// keep track of all rotations. 
		vector<uint64_t> rotations(d);

		int coeff_count = parms_.poly_modulus().coeff_count();
		int coeff_mod_count = parms_.coeff_mod_array().size();
		int coeff_bit_count = coeff_mod_count * bits_per_uint64;
		int coeff_uint64_count = divide_round_up(coeff_bit_count, bits_per_uint64);
		int encrypted_count = encrypted.size();

		rotations.clear();

		for (int i = 0; i < d; i++) {
			//destination[i].get_mutable_array().resize(2, coeff_count, coeff_bit_count);
			//destination[i].get_mutable_array().set_zero();
			//destination[i].hash_block_ = newparms_.get_hash_block();
			//encrypted_rotations[i].get_mutable_array().resize(2, coeff_count, coeff_bit_count);
			//encrypted_rotations[i].get_mutable_array().set_zero();
			//encrypted_rotations[i].hash_block_ = newparms_.get_hash_block();
			rotations.push_back(gal);
			gal *= p;
			gal %= m;
		}

		// Modify below to use matrix vector product....
		// Populate the matrix coefficients. 
		vector<vector<Plaintext>> coeffs(d);
		for (int i = 0; i < d; i++) {
			coeffs[i].resize(d);
			for (int j = 0; j < d; j++) {
				newevaluator_->apply_galois_plain(unpack_polys_[i], rotations[j], coeffs[i][j]);
			}
		}

		cout << "unpack: computing rotations using hoisting ..." << endl;
		auto time_start = chrono::high_resolution_clock::now();
		newevaluator_->apply_galois_many(encrypted, rotations, encrypted_rotations, true);
		auto time_end = chrono::high_resolution_clock::now();
		chrono::microseconds online_time = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
		cout << "+unpack hoisting time: " << (double)online_time.count() / 1000000 << " seconds" << endl;

		//for (int i = 0; i < d; i++) {
		//	newevaluator.apply_galois(encrypted, rotations[i], cipher_rotations[i]);
		//}
#ifdef _DEBUG
		//cout << i << "-th rotation in recrypt_unpack " << decryptor.decrypt(cipher_rotations[i]).to_string() << ", noise budget: " <<
		//    decryptor.invariant_noise_budget(cipher_rotations[i]) << endl;
#endif

		// populate all rotations. 
		//vector<Ciphertext> encrypted_rotations(d);
		//for (int i = 0; i < d; i++) {
		//    newevaluator.apply_galois(encrypted, rotations[i], encrypted_rotations[i]);
		//}

		// apply linear transform. Here we do not want the output to be in ntt form. 
		newevaluator_->matrix_vector_product_plain(encrypted_rotations, coeffs, destination, false);
		return;

		//for (int i = 0; i < d; i++) {
		//    // compute Tr(c[i] x)
		//    destination[i].resize(encrypted.size(), encrypted.coeff_count(), encrypted.coeff_bit_count()); // FIXME
		//    destination[i].set_zero();
		//    for (int j = 0; j < d; j++) {
		//        // sigma_
		//        BigPoly unpack_poly_rotated;    
		//        // Compute Frob^j(c[i]). (rotations[j] = p**j) 
		//        apply_galois_plain(data.unpack_polys()[i], rotations[j], unpack_poly_rotated);
		//        // Compute Frob^j(c[i]) * Frob^j(x)
		//        multiply_plain(cipher_rotations[j], unpack_poly_rotated, temp); 
		//        add(destination[i], temp, destination[i]);
		//    }
		//}
		//return; 
	}

	// read the linearized polys from file. 
	// Assume the file already contains 
	void RNSRecryptor::set_unpack_polys(string filename)
	{
		unpack_polys_.resize(extensionDegree_);

		ifstream ifs(filename.c_str());
		for (int i = 0; i < extensionDegree_; i++) {
			unpack_polys_[i].resize(n_ + 1, newparms_.plain_modulus().bit_count());
			string poly_str;
			getline(ifs, poly_str);
			unpack_polys_[i] = BigPoly(poly_str);
			cout << "read unpack poly " << unpack_polys_[i].to_string() << endl;
		}
		ifs.close();
	}

	void RNSRecryptor::evaluate_twotimes_sigmoid(const Ciphertext & encrypted, Ciphertext & destination, RNSDecryptor &decryptor)
	{
		// This function evaluates the degree 3 approximation to the sigmoid function, 
		// scaled up by 2. 

		// The polynomial is f(x) = -0.008x^3 + 0.394 * 2 x + 1
		// We assume p = 127. Other cases can be implemented... 
		if (parms_.plain_modulus_base() != 127 || parms_.plain_modulus_exponent() != 3) throw;

		uint64_t A0, A1, A3;

		A0 = 127 * 127; // A0 = 1, but  scaled twice. 
		A1 = 50;   // A1 ~= 0.394*127
		A3 = 1; // A3 ~= 0.008*127 

		int n = parms_.poly_modulus().coeff_count() - 1;
		BigPoly constA0;
		constA0.resize(n + 1, parms_.plain_modulus().bit_count());
		constA0[0] = A0;
		BigPoly constA1;
		constA1.resize(n + 1, parms_.plain_modulus().bit_count());
		constA1[0] = A1;

		Ciphertext x = encrypted;

		Ciphertext A1x = oldevaluator_->multiply_plain(encrypted, constA1);

		// We compute the non-scaled version of x^2. 
		Ciphertext xsquare_nonscaled = oldevaluator_->relinearize(oldevaluator_->square(encrypted));




		// The scaled version of x^2
		Ciphertext xsquare = floor_simd(xsquare_nonscaled, 1);

		Ciphertext xfloor = floor_simd(x, 1); // scaled version of (x*1)

		Ciphertext xcube_nonscaled = oldevaluator_->relinearize(oldevaluator_->multiply(xfloor, xsquare));

		auto x_ptxt = oldExPolyCRTBuilder_->decompose(decryptor.rns_decrypt(encrypted));
		auto xsquare_ptxt = oldExPolyCRTBuilder_->decompose(decryptor.rns_decrypt(xsquare));
		auto xfloor_ptxt = oldExPolyCRTBuilder_->decompose(decryptor.rns_decrypt(xfloor));
		auto xcube_nonscaled_ptxt = oldExPolyCRTBuilder_->decompose(decryptor.rns_decrypt(xcube_nonscaled));
		auto a1x_nonscaled_ptxt = oldExPolyCRTBuilder_->decompose(decryptor.rns_decrypt(A1x));


		Ciphertext temp, tempdest;
		oldevaluator_->sub(A1x, xcube_nonscaled, temp); // non-scaled version of a1x - a3x^3
		oldevaluator_->add_plain(temp, constA0, tempdest); // non-scaled version of f(x) = a0 + a1x - a3x^3.

		floor_simd(tempdest, destination, 1);

		auto temp_ptxt = oldExPolyCRTBuilder_->decompose(decryptor.rns_decrypt(temp));
		auto tempdest_ptxt = oldExPolyCRTBuilder_->decompose(decryptor.rns_decrypt(tempdest));
		auto dest_ptxt = oldExPolyCRTBuilder_->decompose(decryptor.rns_decrypt(destination));

		//for (int i = 0; i < numSlots_; i++)
		//{
		//	cout << "x =  " << x_ptxt[i].pointer()[0] << ", ";
		//	cout << "x^2 = " << xsquare_ptxt[i].pointer()[0] << ", ";
		//	cout << "x floor = " << xfloor_ptxt[i].pointer()[0] << ", ";
		//	cout << "x^3 (not scaled) =" << xcube_nonscaled_ptxt[i].pointer()[0] << ", ";
		//	cout << "a1x (not scaled) = " << a1x_nonscaled_ptxt[i].pointer()[0] << ", ";
		//	cout << "a1x  - x^3 (not scaled) = " << temp_ptxt[i].pointer()[0] << ", ";
		//	cout << "a0 + a1x - x^3 (nonscaled) = " << tempdest_ptxt[i].pointer()[0] << ", ";
		//	cout << "a0 + a1x - x^3 = " << dest_ptxt[i].pointer()[0] << ", ";

		//	cout << endl;
		//}
		return;
	}


	void RNSRecryptor::multiply_by_learningrate(const Ciphertext &encrypted, Ciphertext &destination) {
		// In iDASH task the learning rate is 0.0032. 
		// Because we scaled up the sigmoid function by 2, we need to scale down the learning rate
		// so alpha' = 	0.0016.
		// When p = 127, we can achieve this by multiplying with sqrt(alpha')  = 0. 04  and floor 2 times. 

		if (parms_.plain_modulus_base() != 127 || parms_.plain_modulus_exponent() != 3) throw;

		BigPoly constsqrtAlpha;
		constsqrtAlpha.resize(parms_.poly_modulus().coeff_count() - 1, parms_.plain_modulus().bit_count());
		constsqrtAlpha[0] = 25; // 5 ~= 0.04*127, so 25 = 0.0016*127**2.

		Ciphertext first = oldevaluator_->multiply_plain(encrypted, constsqrtAlpha);
		floor_simd(first, destination, 2); //floor down by p^2. This is possible because the input is bounded by p^2/2. 
		// So after multiplying by 25 < p, it is bounded by p^3/2. 
	}

	void RNSRecryptor::compute_unpack_polys(string file_storing_exringelts)
	{
		ifstream ifs(file_storing_exringelts.c_str());
		if (!ifs)
		{
			throw runtime_error("file not found");
		}
		int logn = get_power_of_two(n_);
		unpack_polys_.resize(extensionDegree_);
		vector<ExRingElement> vec(numSlots_);
		for (int i = 0; i < extensionDegree_; i++)
		{
			unpack_polys_[i].resize(n_ + 1, newparms_.plain_modulus().bit_count());
			for (int j = 0; j < numSlots_; j++)
			{
				string poly_str;
				getline(ifs, poly_str);
				vec[j] = ExRingElement(ExRing_, poly_str);
			}
			unpack_polys_[i] = ExPolyCRTBuilder_->compose(vec).get_poly();
			//cout << "computed unpack poly " << unpack_polys_[i].to_string() << endl;
		}
		ifs.close();
	}

	BigPoly biguint_vector_to_bigpoly(vector<BigUInt> coeffs) {
		BigPoly output;
		int coeff_count = coeffs.size() + 1;
		int coeff_bit_count = coeffs[0].bit_count();
		output.resize(coeff_count, coeff_bit_count);
		for (int i = 0; i < coeffs.size(); i++) {
			set_uint_uint(coeffs[i].pointer(), coeffs[i].uint64_count(), output[i].pointer());
		}
		output[coeff_count - 1].set_zero();
		return output;
	}
}