#pragma once

#include <vector>
#include <utility>
#include "rnsencryptionparams.h"
#include "rnscontext.h"
#include "rnsevaluationkeys.h"
#include "util/mempool.h"
#include "smallmodulus.h"
#include "util/modulus.h"
#include "util/polymodulus.h"
#include "util/smallntt.h"
#include "memorypoolhandle.h"
#include "ciphertext.h"
#include "plaintext.h"
#include "util/baseconvertor.h"

///////// New header file ////////
#include "rnsgaloiskeys.h"
#include "utilities.h"

using namespace std;

namespace seal
{
    /**
    Provides arithmetic functions for operating on ciphertexts. The add, subtract, and multiply function
    variants allow both operands to be encrypted. The "_plain" variants allow one of the inputs to be
    encrypted and the other unencrypted.

    Every valid ciphertext consists of at least two polynomials. Homomorphic multiplication increases
    the size of the ciphertext in such a way that if the input ciphertexts have size M and N, then the
    output ciphertext will have size M+N-1. The multiplication operation will require M*N polynomial
    multiplications to be performed. To read the current size of a ciphertext the user can use
    BigPolyArray::size().

    A relinearization operation can be used to reduce the size of a ciphertext to any smaller size
    (but at least 2), potentially improving the performance of a subsequent multiplication using it.
    However, relinearization consumes the invariant noise budget in a ciphertext by an additive factor
    proportional to 2^DBC, and relinearizing from size K to L will require 2*(K-L)*[floor(log_2(coeff_modulus)/DBC)+1]
    polynomial multiplications, where DBC denotes the decomposition bit count set in the encryption parameters.
    Note that the larger the decomposition bit count is, the faster relinearization will be, but also the
    more invariant noise budget will be consumed.

    Relinearization requires the key generator to generate evaluation keys. More specifically, to relinearize
    a ciphertext of size K down to any size smaller than K (but at least 2), at least K-2 evaluation keys will
    be needed. These have to be given as an input parameter to the constructor of RNSEvaluator.

    @par Invariant Noise Budget
    The invariant noise polynomial of a ciphertext is a rational coefficient polynomial, such that
    a ciphertext decrypts correctly as long as the coefficients of the invariant noise polynomial are
    of absolute value less than 1/2. Thus, we call the infinity-norm of the invariant noise polynomial
    the invariant noise, and for correct decryption require it to be less than 1/2. If v denotes the
    invariant noise, we define the invariant noise budget as -log2(2v). Thus, the invariant noise budget
    starts from some initial value, which depends on the encryption parameters, and decreases to 0 when
    computations are performed. When the budget reaches 0, the ciphertext becomes too noisy to decrypt
    correctly.
    */
    class RNSEvaluator
    {
    public:
        /**
        Creates an RNSEvaluator instance initialized with the specified RNSContext and evaluation
        keys. If no evaluation keys will be needed, one can simply pass a newly created empty
        instance of EvaluationKeys to the function. Optionally, the user can give a reference
        to a MemoryPoolHandle object to use a custom memory pool instead of the global memory
        pool (default).

        @param[in] context The RNSContext
        @param[in] evaluation_keys The evaluation keys
        @param[in] pool The memory pool handle
        @throws std::invalid_argument if encryption parameters or evaluation keys are not valid
        @see RNSContext for more details on valid encryption parameters.
        @see MemoryPoolHandle for more details on memory pool handles.
        */
        RNSEvaluator(const RNSContext &context, const RNSEvaluationKeys &evaluation_keys,
            const MemoryPoolHandle &pool = MemoryPoolHandle::acquire_global());

        /**
        Creates an RNSEvaluator instance initialized with the specified RNSContext. Optionally,
        the user can give a reference to a MemoryPoolHandle object to use a custom memory pool
        instead of the global memory pool (default).

        @param[in] context The RNSContext
        @param[in] pool The memory pool handle
        @throws std::invalid_argument if encryption parameters are not valid
        @see RNSContext for more details on valid encryption parameters.
        @see MemoryPoolHandle for more details on memory pool handles.
        */
        RNSEvaluator(const RNSContext &context, const MemoryPoolHandle &pool = MemoryPoolHandle::acquire_global())
            : RNSEvaluator(context, RNSEvaluationKeys(), pool)
        {
        }

        /**
        Creates a copy of a RNSEvaluator.

        @param[in] copy The RNSEvaluator to copy from
        */
        RNSEvaluator(const RNSEvaluator &copy);

        /**
        Creates a new RNSEvaluator by moving an old one.

        @param[in] source The RNSEvaluator to move from
        */
        RNSEvaluator(RNSEvaluator &&source) = default;

        /**
        Negates a ciphertext and stores the result in the destination parameter.

        @param[in] encrypted The ciphertext to negate
        @param[out] destination The ciphertext to overwrite with the negated result
        @throws std::invalid_argument if the ciphertext is not valid for the encryption parameters
        */
        void negate(const Ciphertext &encrypted, Ciphertext &destination);

        /**
        Negates a ciphertext and returns the result.

        @param[in] encrypted The ciphertext to negate
        @throws std::invalid_argument if the ciphertext is not valid for the encryption parameters
        */
        Ciphertext negate(const Ciphertext &encrypted)
        {
            Ciphertext result;
            negate(encrypted, result);
            return result;
        }

        /**
        Adds two ciphertexts and stores the result in the destination parameter.

        @param[in] encrypted1 The first ciphertext to add
        @param[in] encrypted2 The second ciphertext to add
        @param[out] destination The ciphertext to overwrite with the addition result
        @throws std::invalid_argument if the ciphertexts are not valid for the encryption parameters
        */
        void add(const Ciphertext &encrypted1, const Ciphertext &encrypted2, Ciphertext &destination);

        /**
        Adds two ciphertexts and returns the result.

        @param[in] encrypted1 The first ciphertext to add
        @param[in] encrypted2 The second ciphertext to add
        @throws std::invalid_argument if the ciphertexts are not valid for the encryption parameters
        */
        Ciphertext add(const Ciphertext &encrypted1, const Ciphertext &encrypted2)
        {
            Ciphertext result;
            add(encrypted1, encrypted2, result);
            return result;
        }

        /**
        Adds together an number of ciphertexts stored as elements of std::vector<Ciphertext>
        and stores the result in the destination parameter.

        @param[in] encrypteds The ciphertexts to add
        @param[out] destination The ciphertext to overwrite with the addition result
        @throws std::invalid_argument if encrypteds is empty
        @throws std::invalid_argument if the ciphertexts are not valid for the encryption parameters
        */
        void add_many(const std::vector<Ciphertext> &encrypteds, Ciphertext &destination);

        /**
        Adds together an number of ciphertexts stored as elements of std::vector<Ciphertext>
        and returns the result.

        @param[in] encrypteds The ciphertexts to add
        @throws std::invalid_argument if encrypteds is empty
        @throws std::invalid_argument if the ciphertexts are not valid for the encryption parameters
        */
        Ciphertext add_many(const std::vector<Ciphertext> &encrypteds)
        {
            Ciphertext result;
            add_many(encrypteds, result);
            return result;
        }

        /**
        Subtracts two ciphertexts and stores the result in the destination parameter.

        @param[in] encrypted1 The ciphertext to subtract from
        @param[in] encrypted2 The ciphertext to subtract
        @param[out] destination The ciphertext to overwrite with the subtraction result
        @throws std::invalid_argument if the ciphertexts are not valid for the encryption parameters
        */
        void sub(const Ciphertext &encrypted1, const Ciphertext &encrypted2, Ciphertext &destination);

        /**
        Subtracts two ciphertexts and returns the result.

        @param[in] encrypted1 The ciphertext to subtract from
        @param[in] encrypted2 The ciphertext to subtract
        @throws std::invalid_argument if the ciphertexts are not valid for the encryption parameters
        */
        Ciphertext sub(const Ciphertext &encrypted1, const Ciphertext &encrypted2)
        {
            Ciphertext result;
            sub(encrypted1, encrypted2, result);
            return result;
        }

        /**
        Multiplies two ciphertexts and stores the result in the destination parameter.

        @param[in] encrypted1 The first ciphertext to multiply
        @param[in] encrypted2 The second ciphertext to multiply
        @param[out] destination The ciphertext to overwrite with the multiplication result
        @throws std::invalid_argument if the ciphertexts are not valid for the encryption parameters
        */
        void multiply(const Ciphertext &encrypted1, const Ciphertext &encrypted2, Ciphertext &destination);

        /**
        Multiplies two ciphertexts without performing relinearization, and returns the result.

        @param[in] encrypted1 The first ciphertext to multiply
        @param[in] encrypted2 The second ciphertext to multiply
        @throws std::invalid_argument if the ciphertexts are not valid for the encryption parameters
        */
        Ciphertext multiply(const Ciphertext &encrypted1, const Ciphertext &encrypted2)
        {
            Ciphertext result;
            multiply(encrypted1, encrypted2, result);
            return result;
        }

        /**
        Squares a ciphertext and stores the result in the destination parameter.

        @param[in] encrypted The ciphertext to square
        @param[out] destination The ciphertext to overwrite with the result
        @throws std::invalid_argument if the ciphertexts are not valid for the encryption parameters
        */
        void square(const Ciphertext &encrypted, Ciphertext &destination);

        /**
        Squares a ciphertext and returns the result.

        @param[in] encrypted The ciphertext to square
        @throws std::invalid_argument if the ciphertexts are not valid for the encryption parameters
        */
        Ciphertext square(const Ciphertext &encrypted)
        {
            Ciphertext result;
            square(encrypted, result);
            return result;
        }

        /**
        Relinearizes a ciphertext and stores the result in the destination parameter.

        @param[in] encrypted The ciphertext to relinearize
        @param[in] destination_size The size of the output ciphertext (defaults to 2)
        @param[out] destination The ciphertext to overwrite with the relinearized result
        @throws std::invalid_argument if the ciphertext is not valid for the encryption parameters
        @throws std::invalid_argument if destination_size is less than 2 or greater than number of elements currently in ciphertext
        @throws std::invalid_argument if not enough evaluation keys have been generated
        */
        void relinearize(const Ciphertext &encrypted, Ciphertext &destination, int destination_size = 2);

        /**
        Relinearizes a ciphertext and returns the result.

        @param[in] encrypted The ciphertext to relinearize
        @param[in] destination_size The size of the output ciphertext (defaults to 2)
        @throws std::invalid_argument if the ciphertext is not valid for the encryption parameters
        @throws std::invalid_argument if destination_size is less than 2 or greater than number of elements currently in ciphertext
        @throws std::invalid_argument if not enough evaluation keys have been generated
        */
        Ciphertext relinearize(const Ciphertext &encrypted, int destination_size = 2)
        {
            Ciphertext result;
            relinearize(encrypted, result, destination_size);
            return result;
        }

        /**
        Multiplies a std::vector of ciphertexts together and returns the result. Relinearization is performed after
        every multiplication, so enough encryption keys must have been given to the constructor of the RNSEvaluator.

        @param[in] encrypteds The std::vector of ciphertexts to multiply
        @throws std::invalid_argument if the encrypteds std::vector is empty
        @throws std::invalid_argument if the ciphertexts are not valid ciphertexts for the encryption parameters
        */
        Ciphertext multiply_many(std::vector<Ciphertext> encrypteds);

        /**
        Multiplies a std::vector of ciphertexts together and stores the result in the destination parameter.
        Relinearization is performed after every multiplication, so enough encryption keys must have been given
        to the constructor of the RNSEvaluator.

        @param[in] encrypteds The std::vector of ciphertexts to multiply
        @param[out] destination The ciphertext to overwrite with the multiplication result
        @throws std::invalid_argument if the encrypteds std::vector is empty
        @throws std::invalid_argument if the ciphertexts are not valid ciphertexts for the encryption parameters
        */
        void multiply_many(std::vector<Ciphertext> &encrypteds, Ciphertext &destination)
        {
            destination = multiply_many(encrypteds);
        }

        /**
        Raises a ciphertext to the specified power and stores the result in the destination parameter.

        Exponentiation to power 0 is not allowed and will result in the library throwing an invalid argument
        exception. The reason behind this design choice is that the result should be a fresh encryption
        of 1, but creating fresh encryptions should not be something this class does. Instead the user
        should separately handle the cases where the exponent is 0. Relinearization is performed after
        every multiplication, so enough encryption keys must have been given to the constructor of the RNSEvaluator.

        @param[in] encrypted The ciphertext to raise to a power
        @param[in] exponent The power to raise the ciphertext to
        @param[out] destination The ciphertext to overwrite with the exponentiation result
        @throws std::invalid_argument if the ciphertext is not valid for the encryption parameters
        @throws std::invalid_argument if the exponent is zero
        */
        void exponentiate(const Ciphertext &encrypted, std::uint64_t exponent, Ciphertext &destination);

        /**
        Raises a ciphertext to the specified power and returns the result.

        Exponentiation to power 0 is not allowed and will result in the library throwing
        an invalid argument exception. The reason behind this design choice is that the
        result should be a fresh encryption of 1, but creating fresh encryptions should
        not be something this class does. Instead the user has to separately handle
        the cases where the exponent is 0. Relinearization is performed after every multiplication,
        so enough encryption keys must have been given to the constructor of the RNSEvaluator.

        @param[in] encrypted The ciphertext to raise to a power
        @param[in] exponent The power to raise the ciphertext to
        @throws std::invalid_argument if the ciphertext is not valid for the encryption parameters
        @throws std::invalid_argument if the exponent is zero
        */
        Ciphertext exponentiate(const Ciphertext &encrypted, std::uint64_t exponent)
        {
            Ciphertext result;
            exponentiate(encrypted, exponent, result);
            return result;
        }

        /**
        Adds a ciphertext with a plaintext, and stores the result in the destination
        parameter. The plaintext must have a significant coefficient count smaller than
        the coefficient count specified by the encryption parameters, and with
        coefficient values less than the plain modulus (RNSContext::plain_modulus()).

        @param[in] encrypted The ciphertext to add
        @param[in] plain The plaintext to add
        @param[out] destination The ciphertext to overwrite with the addition result
        @throws std::invalid_argument if the ciphertexts are not valid for the encryption parameters
        @throws std::invalid_argument if the plaintext's significant coefficient count or coefficient
        values are too large to represent with the encryption parameters
        */
        void add_plain(const Ciphertext &encrypted, const Plaintext &plain, Ciphertext &destination);

        /**
        Adds a ciphertext with a plaintext, and returns the result. The plaintext
        must have a significant coefficient count smaller than the coefficient count specified by the
        encryption parameters, and with coefficient values less than the plain modulus
        (RNSContext::plain_modulus()).

        @param[in] encrypted The ciphertext to add
        @param[in] plain The plaintext to add
        @throws std::invalid_argument if the ciphertexts are not valid for the encryption parameters
        @throws std::invalid_argument if the plaintext's significant coefficient count or coefficient
        values are too large to represent with the encryption parameters
        */
        Ciphertext add_plain(const Ciphertext &encrypted, const Plaintext &plain)
        {
            Ciphertext result;
            add_plain(encrypted, plain, result);
            return result;
        }

        /**
        Subtracts a ciphertext with a plaintext, and stores the result in the destination
        parameter. The plaintext must have a significant coefficient count smaller than
        the coefficient count specified by the encryption parameters, and with
        coefficient values less than the plain modulus (RNSContext::plain_modulus()).

        @param[in] encrypted The ciphertext to subtract from
        @param[in] plain The plaintext to subtract
        @param[out] destination The ciphertext to overwrite with the subtraction result
        @throws std::invalid_argument if the ciphertexts are not valid for the encryption parameters
        @throws std::invalid_argument if the plaintext's significant coefficient count or coefficient
        values are too large to represent with the encryption parameters
        */
        void sub_plain(const Ciphertext &encrypted, const Plaintext &plain, Ciphertext &destination);

        /**
        Subtracts a ciphertext with a plaintext, and returns the result. The plaintext
        must have a significant coefficient count smaller than the coefficient count specified by the
        encryption parameters, and with coefficient values less than the plain modulus
        (RNSContext::plain_modulus()).

        @param[in] encrypted The ciphertext to subtract from
        @param[in] plain The plaintext to subtract
        @throws std::invalid_argument if the ciphertexts are not valid for the encryption parameters
        @throws std::invalid_argument if the plaintext's significant coefficient count or coefficient
        values are too large to represent with the encryption parameters
        */
        Ciphertext sub_plain(const Ciphertext &encrypted, const Plaintext &plain)
        {
            Ciphertext result;
            sub_plain(encrypted, plain, result);
            return result;
        }

        /**
        Multiplies a ciphertext with a plaintext, and stores the result in the destination
        parameter. The plaintext must have a significant coefficient count smaller than the
        coefficient count specified by the encryption parameters, and with coefficient values
        less than the plain modulus (RNSContext::plain_modulus()).

        Multiplying by a plaintext 0 is not allowed and will result in the library throwing an invalid
        argument exception. The reason behind this design choice is that the result should
        be a fresh encryption of 0, but creating fresh encryptions should not be something
        this class does. Instead the user should separately handle the cases where the
        plain multiplier is 0.

        @param[in] encrypted The ciphertext to multiply
        @param[in] plain The plaintext to multiply
        @param[out] destination The ciphertext to overwrite with the multiplication result
        @throws std::invalid_argument if the encrypted is not valid for the encryption parameters
        @throws std::invalid_argument if the plaintext's significant coefficient count or coefficient
        values are too large to represent with the encryption parameters
        @throws std::invalid_argument if the plaintext multiplier is zero
        */
        void multiply_plain(const Ciphertext &encrypted, const Plaintext &plain, Ciphertext &destination);

        /**
        Multiplies a ciphertext with a plaintext, and returns the result. The plaintext
        must have a significant coefficient count smaller than the coefficient count specified by the
        encryption parameters, and with coefficient values less than the plain modulus
        (RNSContext::plain_modulus()).

        Multiplying by a plaintext 0 is not allowed and will result in the library throwing an invalid
        argument exception. The reason behind this design choice is that the result should
        be a fresh encryption of 0, but creating fresh encryptions should not be something
        this class does. Instead the user should separately handle the cases where the
        plain multiplier is 0.

        @param[in] encrypted The ciphertext to multiply
        @param[in] plain The plaintext to multiply
        @throws std::invalid_argument if the encrypted is not valid for the encryption parameters
        @throws std::invalid_argument if the plaintext's significant coefficient count or coefficient
        values are too large to represent with the encryption parameters
        @throws std::invalid_argument if the plaintext multiplier is zero
        */
        Ciphertext multiply_plain(const Ciphertext &encrypted, const Plaintext &plain)
        {
            Ciphertext result;
            multiply_plain(encrypted, plain, result);
            return result;
        }

        /**
        Returns the evaluation keys used by the RNSEvaluator.
        */
        const RNSEvaluationKeys &evaluation_keys() const
        {
            return evaluation_keys_;
        }

        /*
        Transform a plaintext from coefficient domain to NTT domain, with respect to the coefficient
        modulus. This function first embeds integers modulo the plaintext modulus to integers modulo
        the coefficient modulus, and then performs David Harvey's NTT on the resulting polynomial.

        Ciphertexts in the NTT domain can be added as usual, and multiplied by plaintext polynomials
        (also in the NTT domain) using multiply_plain_ntt, but cannot be homomorphically multiplied
        with other ciphertexts without first transforming both inputs to coefficient domain with
        transform_from_ntt.

        @param[in] plain The plaintext to transform
        @throws std::logic_error if the encryption parameters do not support NTT
        @throws std::invalid_argument if plain is not valid for the encryption parameters
        */
        void transform_to_ntt(Plaintext &plain);

		void transform_to_ntt(const Plaintext &plain, Plaintext & destination) {
			destination = plain;
			transform_to_ntt(destination);
			return;
		}

        /*
        Transform a plaintext from NTT domain to coefficient domain, with respect to the coefficient
        modulus. This function first performs David Harvey's inverse NTT, and follows it by an inverse
        of the coefficient embedding performed by transform_to_ntt(BigPoly &plain).

        Ciphertexts in the NTT domain can be added as usual, and multiplied by plaintext polynomials
        (also in the NTT domain) using multiply_plain_ntt, but cannot be homomorphically multiplied
        with other ciphertexts without first transforming both inputs to coefficient domain with
        transform_from_ntt.

        @param[in] plain_ntt The plaintext to transform
        @throws std::logic_error if the encryption parameters do not support NTT
        @throws std::invalid_argument if plain_ntt is not valid for the encryption parameters
        */
        void transform_from_ntt(Plaintext &plain_ntt);

        /*
        Transform a ciphertext from coefficient domain to NTT domain, with respect to the coefficient
        modulus. This function performs David Harvey's NTT separately on each of the polynomials
        in the given BigPolyArray.

        Ciphertexts in the NTT domain can be added as usual, and multiplied
        by plaintext polynomials (also in the NTT domain) using multiply_plain_ntt, but cannot be
        homomorphically multiplied with other ciphertexts without first transforming both inputs
        to coefficient domain with transform_from_ntt.

        @param[in] encrypted The ciphertext to transform
        @throws std::logic_error if the encryption parameters do not support NTT
        @throws std::invalid_argument if encrypted is not valid for the encryption parameters
        */
        void transform_to_ntt(Ciphertext &encrypted);

        /*
        Transform a ciphertext from NTT domain to coefficient domain, with respect to the coefficient
        modulus. This function performs David Harvey's inverse NTT separately on each of the polynomials
        in the given BigPolyArray.

        Ciphertexts in the NTT domain can be added as usual, and multiplied by plaintext polynomials
        (also in the NTT domain) using multiply_plain_ntt, but cannot be homomorphically multiplied
        with other ciphertexts without first transforming both inputs to coefficient domain with
        transform_from_ntt.

        @param[in] encrypted_ntt The ciphertext to transform
        @throws std::logic_error if the encryption parameters do not support NTT
        @throws std::invalid_argument if encrypted_ntt is not valid for the encryption parameters
        */
        void transform_from_ntt(Ciphertext &encrypted_ntt);

        /*
        Multiplies a ciphertext with a plaintext, assuming both are already transformed to NTT domain.
        The result ciphertext remains in the NTT domain, and can be transformed back to coefficient
        domain with transform_from_ntt.

        Ciphertexts in the NTT domain can be added as usual, and multiplied by plaintext polynomials
        (also in the NTT domain) using multiply_plain_ntt, but cannot be homomorphically multiplied
        with other ciphertexts without first transforming both inputs to coefficient domain with
        transform_from_ntt.

        @param[in] encrypted_ntt The ciphertext to multiply (in NTT domain)
        @param[in] plain_ntt The plaintext to multiply (in NTT domain)
        @param[out] destination_ntt The ciphertext to overwrite with the multiplication result (in NTT domain)
        @throws std::logic_error if the encryption parameters do not support NTT
        @throws std::invalid_argument if encrypted_ntt is not valid for the encryption parameters
        @throws std::invalid_argument if the plaintext's significant coefficient count or coefficient
        values are too large to represent with the encryption parameters
        @throws std::invalid_argument if the plaintext multiplier is zero
        */
        void multiply_plain_ntt(const Ciphertext &encrypted_ntt, const Plaintext &plain_ntt, Ciphertext &destination_ntt);

        /*
        Multiplies a ciphertext with a plaintext, assuming both are already transformed to NTT domain,
        and returns the result. The result ciphertext remains in the NTT domain, and can be transformed
        back to coefficient domain with transform_from_ntt.

        Ciphertexts in the NTT domain can be added as usual, and multiplied by plaintext polynomials
        (also in the NTT domain) using multiply_plain_ntt, but cannot be homomorphically multiplied
        with other ciphertexts without first transforming both inputs to coefficient domain with
        transform_from_ntt.

        @param[in] encrypted_ntt The ciphertext to multiply (in NTT domain)
        @param[in] plain_ntt The plaintext to multiply (in NTT domain)
        @throws std::logic_error if the encryption parameters do not support NTT
        @throws std::invalid_argument if encrypted_ntt is not valid for the encryption parameters
        @throws std::invalid_argument if the plaintext's significant coefficient count or coefficient
        values are too large to represent with the encryption parameters
        @throws std::invalid_argument if the plaintext multiplier is zero
        */
        Ciphertext multiply_plain_ntt(const Ciphertext &encrypted_ntt, const Plaintext &plain_ntt)
        {
            Ciphertext result_ntt;
            multiply_plain_ntt(encrypted_ntt, plain_ntt, result_ntt);
            return result_ntt;
        }

		////////////// New Functions //////////////

		void set_galois_keys(const RNSGaloisKeys &galois_keys)
		{
			galois_keys_ = galois_keys;
		}

		void Polynomial_Evaluation(const Ciphertext &encrypted, vector<uint64_t> &coeff, Ciphertext &destination, const SmallModulus &eval_plain_modulus_);

		Ciphertext Polynomial_Evaluation(const Ciphertext &encrypted, vector<uint64_t> &coeff, const SmallModulus &eval_plain_modulus_)
		{
			Ciphertext result;
			Polynomial_Evaluation(encrypted, coeff, result, eval_plain_modulus_);
			return result;
		}

		void DigitExtraction(Ciphertext &encrypted, vector<Ciphertext> &destination, long p, long e, long r, bool centered = false);

		vector<Ciphertext> DigitExtraction(Ciphertext &encrypted, long p, long e, long r, bool centered = false)
		{
			vector<Ciphertext> result;
			DigitExtraction(encrypted, result, p, e, r, centered);
			return result;
		}

		void apply_galois(const Ciphertext &encrypted, uint64_t galois_elt, Ciphertext &destination, bool ntt_output = false);

		Ciphertext apply_galois(const Ciphertext &encrypted, uint64_t galois_elt, bool ntt_output = false)
		{
			Ciphertext result;
			apply_galois(encrypted, galois_elt, result, ntt_output);
			return result;
		}

		void apply_galois_many(const Ciphertext &encrypted, vector<uint64_t> galois_elt, vector<Ciphertext> &destination, bool ntt_output = false);

		vector<Ciphertext> apply_galois_many(const Ciphertext &encrypted, vector<uint64_t> galois_elt, bool ntt_output = false)
		{
			vector<Ciphertext> result;
			apply_galois_many(encrypted, galois_elt, result, ntt_output);
			return result;
		}

        void matrix_vector_product_plain(std::vector<Ciphertext> &encrypted, const std::vector<std::vector<Plaintext>> &coeff_matrix, std::vector<Ciphertext> &destination, bool ntt_output);

		void matrix_vector_product_plain_ntt(std::vector<Ciphertext> &encrypted, const std::vector<std::vector<Plaintext>> &coeff_matrix_ntt, std::vector<Ciphertext> &destination, bool ntt_output);

		void apply_galois_plain(const Plaintext & poly, std::uint64_t galois_elt, Plaintext & destination);

		void apply_galois_plain(Plaintext & poly, std::uint64_t galois_elt);

        void rotations_weighted_sum(const Ciphertext &encrypted, std::vector<std::vector<Plaintext>> &coeffs, std::vector<std::uint64_t> babystep, std::vector<std::uint64_t> giantstep, Ciphertext &destination, bool ntt_output);

		void rotations_weighted_sum(const Ciphertext & encrypted, const vector<Plaintext> &coeffs, vector<uint64_t> babystep, vector<uint64_t> giantstep, Ciphertext &destination, bool ntt_output = false);

		void rotations_weighted_sum_coeffs_ntt(const Ciphertext & encrypted, const vector<vector<Plaintext> > &coeffs_ntt, vector<uint64_t> babystep, vector<uint64_t> giantstep, Ciphertext &destination, bool ntt_output = false);

		void expand(const Ciphertext & encrypted, vector<Ciphertext>& destination, vector<Plaintext>& unit_vecs, vector<uint64_t> rotations_index);


		void combine(vector<Ciphertext>& encrypted, Ciphertext & destination, vector<Plaintext> &unit_vecs);
    private:
        RNSEvaluator &operator =(const RNSEvaluator &assign) = delete;

        RNSEvaluator &operator =(RNSEvaluator &&assign) = delete;

        void rns_preencrypt(const std::uint64_t *plain, int plain_coeff_count, int plain_coeff_uint64_count, std::uint64_t *destination);

        void rns_decompose(uint64_t *value);

        void rns_compose(uint64_t *value);

        void relinearize_one_step(const std::uint64_t *encrypted, int encrypted_size, std::uint64_t *destination, util::MemoryPool &pool);

        MemoryPoolHandle pool_;

        RNSEncryptionParameters parms_;

        RNSEncryptionParameterQualifiers qualifiers_;
        
        util::BaseConvertor base_convertor_;
        
        std::vector<util::SmallNTTTables> coeff_small_ntt_tables_;

        std::vector<util::SmallNTTTables> bsk_small_ntt_tables_;

        RNSEvaluationKeys evaluation_keys_;

        util::Pointer upper_half_threshold_;

        util::Pointer upper_half_increment_;

        util::Pointer coeff_div_plain_modulus_;

        util::Pointer plain_upper_half_threshold_;

        util::Pointer plain_upper_half_increment_;

		std::vector<uint64_t> plain_upper_half_increment_array_;

        util::Pointer coeff_modulus_div_two_;

        util::Pointer coeff_products_array_;

        util::Pointer product_modulus_;

        util::Modulus mod_;

        util::PolyModulus polymod_;

        std::vector<SmallModulus> coeff_mod_array_;

        std::vector<SmallModulus> bsk_mod_array_;

        std::vector<std::uint64_t> inv_coeff_products_mod_coeff_array_;

        int bsk_base_mod_count_;

        ////////////// New Functions and variables //////////////

        RNSGaloisKeys galois_keys_;

        map<int, pair<int, int> > map_from_Zmstar_to_generator;

        void PatersonStockmeyer(vector<Ciphertext> &baby_step, vector<Ciphertext> &giant_step, Ciphertext &encrypted, vector<uint64_t> &coeff, Ciphertext &destination, int k, int step);


		void compute_all_powers(const Ciphertext &encrypted, int degree, vector<Ciphertext> &destination);

        void Polynomial_Evaluation(const vector<Ciphertext> &all_powers_encrypted, vector<uint64_t> &coeff, Ciphertext &destination);


        friend class RNSRecryptor;
    };
}
