#pragma once

#include <vector>
#include <utility>
#include <map>
#include <chrono>
#include <thread>
#include "util/mempool.h"
#include "util/modulus.h"
#include "util/polymodulus.h"
#include "util/ntt.h"
#include "util/exring.h"
#include "util/expolycrt.h"
#include "utilities.h"
#include "encryptionparams.h"
#include "evaluationkeys.h"
#include "secretkey.h"
#include "publickey.h"
#include "memorypoolhandle.h"
#include "ciphertext.h"
#include "bigpoly.h"
#include "keygenerator.h"
#include "rnsencryptor.h"
#include "rnsdecryptor.h"
#include "rnsevaluator.h"

using namespace std;
using namespace seal::util;

namespace seal
{
	/**
	Contstructing an Recryptor requires the encryption parameters (set through an SEALContext object)
	and the secret key and evaluation key.
	*/
	class RNSRecryptor
	{
	public:
		/**
		Creates an Recryptor instance initialized with the specified SEALContext, secret_key and
		evaluation_keys.

		n = poly_modulus degree
		p = plain_modulus_base
		r = plain_modulus_exponent
		You should give f(x) and e to contruct Recryptor
		(e.g Recryptor recryptor(context, secret_key, evaluation_key, f(x), e))
		(f(x) is 1st factor of 1X^n + 1 in modulo p^e, and input should be hex_poly form)
		Following is example paramter tables for recrypt.
		/-----------------------------------------------------------------------\
		| n     | p                    | r | e  | f(x)                          |
		| ------|----------------------|---|----|-------------------------------|
		| 16    | 3					   | 1 | 5  | 1x^8 + DDx^4 + F2             |
		| 32    | 3					   | 1 | 5  | 1x^16 + DDx^8 + F2            |
		| 64    | 3					   | 1 | 5  | 1x^32 + DDx^16 + F2           |
		| 128   | 3					   | 1 | 5  | 1x^64 + DDx^32 + F2           |
		| 256   | 3					   | 1 | 5  | 1x^128 + DDx^64 + F2          |
		| 512   | 3					   | 1 | 5  | 1x^256 + DDx^128 + F2         |
		| 1024  | 3					   | 1 | 5  | 1x^512 + DDx^256 + F2         |
		| 2048  | 3					   | 1 | 5  | 1x^1024 + DDx^512 + F2        |

		| 16    | 3					   | 2 | 6  | 1x^8 + 1FCx^4 + 2D8           |
		| 32    | 3					   | 2 | 6  | 1x^16 + 1FCx^8 + 2D8          |
		| 64    | 3					   | 2 | 6  | 1x^32 + 1FCx^16 + 2D8         |
		| 128   | 3					   | 2 | 6  | 1x^64 + 1FCx^32 + 2D8         |
		| 256   | 3					   | 2 | 6  | 1x^128 + 1FCx^64 + 2D8        |
		| 512   | 3					   | 2 | 6  | 1x^256 + 1FCx^128 + 2D8       |
		| 1024  | 3					   | 2 | 6  | 1x^512 + 1FCx^256 + 2D8       |
		| 2048  | 3					   | 2 | 6  | 1x^1024 + 1FCx^512 + 2D8      |

		| 32    | 31                   | 1 | 3  | 1x^2 + 169Bx^1 + 745E         |
		| 64    | 31                   | 1 | 3  | 1x^4 + 169Bx^2 + 745E         |
		| 128   | 31                   | 1 | 3  | 1x^8 + 169Bx^4 + 745E         |
		| 256   | 31                   | 1 | 3  | 1x^16 + 169Bx^8 + 745E        |
		| 512   | 31                   | 1 | 3  | 1x^32 + 169Bx^16 + 745E       |
		| 1024  | 31                   | 1 | 3  | 1x^64 + 169Bx^32 + 745E       |
		| 2048  | 31                   | 1 | 3  | 1x^128 + 169Bx^64 + 745E      |

		| 128   | 127                   | 1 | 2  | 1x^2 + 2C3x^1 + 3F00         |
		| 256   | 127                   | 1 | 2  | 1x^4 + 2C3x^2 + 3F00         |
		| 512   | 127                   | 1 | 2  | 1x^8 + 2C3x^4 + 3F00         |
		| 1024  | 127                   | 1 | 2  | 1x^16 + 2C3x^8 + 3F00        |
		| 2048  | 127                   | 1 | 2  | 1x^32 + 2C3x^16 + 3F00       |


		| 2048  | 257					| 1 | 2 |  1x^16 + 869A					|
		| 4096  | 257					| 1 | 2 |  1x^32 + 869A					|
		| 8192  | 257					| 1 | 2 |  1x^64 + 869A					|
		| 16384 | 257					| 1 | 2 |  1x^128 + 869A				|
		\----------------------------------------------------------------------/
		*/
		RNSRecryptor(const RNSContext &context, SecretKey &secret_key, const RNSEvaluationKeys &evaluation_key, const RNSGaloisKeys &galois_key,
			const std::string &hex_poly, int e, bool nonBatch = false, const MemoryPoolHandle &pool = MemoryPoolHandle::acquire_global());

		/**
		Creates a new Recryptor by moving an old one.

		@param[in] source The Recryptor to move from
		*/
		RNSRecryptor(RNSRecryptor &&source) = default;

		void recrypt(Ciphertext &encrypted, Ciphertext &destination);
		Ciphertext recrypt(Ciphertext &encrypted);

		void recrypt_nonBatch(Ciphertext &encrypted, Ciphertext &destination);
		Ciphertext recrypt_nonBatch(Ciphertext &encrypted);

		void floor(Ciphertext &ecnrypted, Ciphertext &destination, bool invLT = true);
		Ciphertext floor(Ciphertext &encrypted, bool invLT = true);


		void floor_batch(Ciphertext &ecnrypted, vector<Ciphertext> &destination, int num_Batch, int num_thread = 1);

		vector<Ciphertext> floor_batch(Ciphertext &encrypted, int num_Batch, int num_thread = 1);

		void combine(vector<Ciphertext>& encrypted, Ciphertext & destination);

		void expand(const Ciphertext & encrypted, vector<Ciphertext> & destination);
		
		void floor_simd(Ciphertext &ecnrypted, Ciphertext &destination, int k=1, int num_thread = 1);
		Ciphertext floor_simd(Ciphertext &encrypted, int k=1, int num_thread = 1);

		void multiply_by_learningrate(const Ciphertext &encrypted, Ciphertext &destination);


		void recrypt_pre_linear_transform(const Ciphertext & encrypted, Ciphertext & destination);

		void recrypt_unpack(const Ciphertext & encrypted, vector<Ciphertext>& destination);


		void set_unpack_polys(string filename);


		void evaluate_twotimes_sigmoid(const Ciphertext & encrypted, Ciphertext & destination, RNSDecryptor & decryptor);

		// Computing the unpack polynomials using generalized batching. 
		// Input is a file of numSlots_*d rows, and each consecutive numSlots_ rows 
		// is batched in a plaintext polynomial f_i. 
		// Output is plaintext f_1, ..., f_d. 
		void compute_unpack_polys(string file_storing_exringelts);

		void recrypt_linear_transform(const Ciphertext & encrypted, vector<Ciphertext>& destination) {
			Ciphertext intermediate_dest;
			auto time_pre_start = chrono::high_resolution_clock::now();
			recrypt_pre_linear_transform(encrypted, intermediate_dest);
			auto time_pre_end = chrono::high_resolution_clock::now();
			chrono::microseconds time_pre = chrono::duration_cast<chrono::microseconds>(time_pre_end - time_pre_start);
			cout << "Pre linear Transform time : " << (double)time_pre.count() / 1000000 << " seconds" << endl;
#ifdef DEBUG_
			BigPoly plain = decryptor.decrypt(intermediate_dest);
			cout << "intermediate poly = " << plain.to_string() << endl;
#endif 
			recrypt_unpack(intermediate_dest, destination);
			auto time_unpack_end = chrono::high_resolution_clock::now();
			chrono::microseconds time_unpack = chrono::duration_cast<chrono::microseconds>(time_unpack_end - time_pre_end);
			cout << "Unpack  time : " << (double)time_unpack.count() / 1000000 << " seconds" << endl;
		}

		vector<uint64_t> cosets()
		{
			return cosets_;
		}

		vector<uint64_t> classes()
		{
			return classes_;
		}

	private:
		void compute_new_parameters();
		void compute_cosets_classic();
		void compute_cosets();
		// left multiplication of cosets[j] on cosets[i]. 
		int coset_permutation_action(int j, int i);
		void compute_transformation_polys(const std::string &hex_poly);
		void generate_recrypt_key_and_new_evaluator(SecretKey &secret_key);
		void compute_lift_removelsd_poly();
		//void LinearTransformation(const Ciphertext &encrypted, vector<Ciphertext> &destination);
		void LinearTransformation_nonBatch(const Ciphertext &encrypted, Ciphertext &destination, int num_Batch = 0, int num_thread = 1);
		void InverseLinearTransformation(const vector<Ciphertext> &encrypted, Ciphertext &destination);
		void InverseLinearTransformation_nonBatch(const Ciphertext &encrypted, Ciphertext &destination);
		void DigitExtraction(Ciphertext &encrypted, Ciphertext &destination);

		bool nonBatch_;
		int n_;
		int m_;
		int p_;
		int e_;
		int r_;
		int numSlots_;
		int extensionDegree_;
		uint64_t nInverse_;
		vector<uint64_t> cosets_;
		vector<uint64_t> classes_;
		map<int, vector<int> > permute_index_;

		//Polynomial for Digit Extraction
		vector<vector<uint64_t> > lift_poly_;
		vector<vector<uint64_t> > remainlsd_poly_;

		//Polynomial for Linear and invLinear Transformation
		vector<BigPoly> hpolys_;
		vector<BigPoly> unpack_polys_;
		vector<Plaintext> g_polys_ptxt_;
		vector<Plaintext> h_polys_ptxt_;
		vector<vector<Plaintext>> prelinear_polys_;
		vector<vector<Plaintext>> g_poly_ptxt_plain_coeffs_ntt_;
		Plaintext g1_ntt_;
		Plaintext g1_squared_ntt_;

		//For Smart rotations (BSGS technique)
		vector<uint64_t> babystep_;
		vector<uint64_t> giantstep_;

		SecretKey secret_key_;
		unique_ptr<RNSDecryptor> decryptor_;

		Ciphertext recrypt_key_;

		shared_ptr<ExRing> oldExRing_;
		unique_ptr<ExPolyCRTBuilder> oldExPolyCRTBuilder_;
		vector<Plaintext> old_unit_plain_;
		shared_ptr<ExRing> ExRing_;
		unique_ptr<ExPolyCRTBuilder> ExPolyCRTBuilder_;
		vector<Plaintext> unit_plain_;

		unique_ptr<RNSEvaluator> oldevaluator_;
		unique_ptr<RNSEvaluator> newevaluator_;

		RNSEncryptionParameters parms_;
		RNSEncryptionParameters newparms_;

		unique_ptr<RNSContext> oldcontext_;
		unique_ptr<RNSContext> newcontext_;

		MemoryPoolHandle pool_;

		//////////////////////////////////////////
		vector<uint64_t> power_of_three_coset_index;

	};
	BigPoly biguint_vector_to_bigpoly(vector<BigUInt> coeffs);
}