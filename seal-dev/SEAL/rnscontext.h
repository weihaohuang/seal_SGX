#pragma once
#pragma once

#include <iostream>
#include <utility>
#include <string>
#include <array>
#include "rnsencryptionparams.h"
#include "biguint.h"
#include "bigpoly.h"
#include "util/smallntt.h"
#include "randomgen.h"
#include "memorypoolhandle.h"
#include "util/baseconvertor.h"

namespace seal
{
	/**
	Stores a set of attributes (qualifiers) of a set of encryption parameters.
	These parameters are used in various parts of the library, e.g. to determine which
	algorithms can be used. The qualifiers are silently passed on to classes such as
	Encryptor, Evaluator, and Decryptor, and the only way to change them is by changing
	the encryption parameters accordingly. In other words, a user will never have to
	create their own instance of EncryptionParameterQualifiers.

	@see RNSEncryptionParameters::GetQualifiers for obtaining the EncryptionParameterQualifiers
	corresponding to a certain parameter set.
	*/
	struct RNSEncryptionParameterQualifiers
	{
		/**
		If the encryption parameters are set in a way that is considered valid by SEAL, the
		variable parameters_set will be set to true.
		*/
		bool parameters_set;

		/**
		If RNSEncryptionParameters::decomposition_bit_count is set to a positive value, the variable
		enable_relinearization will be set to true.
		*/
		bool enable_relinearization;

		/**
		If the polynomial modulus is of the form X^N+1, where N is a power of two, then
		Nussbaumer convolution can be used for fast multiplication of polynomials modulo
		the polynomial modulus. In this case the variable enable_nussbaumer will be set to true.
		However, currently SEAL requires the polynomial modulus to be of this form to even
		consider the parameters to be valid. Therefore, parameters_set can only be true if
		enable_nussbaumer is true.
		*/
		bool enable_nussbaumer;

		/**
		If the coefficient modulus is congruent to 1 modulo 2N, where X^N+1 is the polynomial
		modulus and N is a power of two, then the number-theoretic transform (NTT) can be used
		for fast multiplications of polynomials modulo the polynomial modulus and coefficient
		modulus. In this case the variable enable_ntt will be set to true.
		*/
		bool enable_ntt;

		/**
		If the plaintext modulus is congruent to 1 modulo 2N, where X^N+1 is the polynomial
		modulus and N is a power of two, then it is possible to use PolyCRTBuilder to do batching,
		which is a fundamental technique in homomorphic encryption to enable powerful SIMD
		functionality, often called "batching" in homomorphic encryption literature. In this
		case the variable enable_batching will be set to true.
		*/
		bool enable_batching;

		bool enable_fast_plain_lift;

	private:
		RNSEncryptionParameterQualifiers() :
			parameters_set(false),
			enable_relinearization(false),
			enable_nussbaumer(false),
			enable_ntt(false),
			enable_batching(false)
		{
		}

		friend class RNSContext;
	};

	/**
	Performs sanity checks (validation) and pre-computations for a given set of encryption
	parameters.

	After the user has set at least the poly_modulus, the coeff_modulus, and the plain_modulus
	parameters in a given RNSEncryptionParameters instance, the parameters need to be validated for
	correctness and functionality, and appropriate pre-computations need to be performed,
	before the parameters can be used for encryption. The constructor of RNSContext
	does this automatically, and concludes by constructing and storing an instance of the
	EncryptionParameterQualifiers class. If the returned instance of EncryptionParameterQualifiers
	has the EncryptionParameterQualifiers::parameters_set flag set to true, the parameter set
	is valid and ready to be used. If the parameters were for some reason not appropriately set,
	the returned EncryptionParameterQualifiers instance will have
	EncryptionParameterQualifiers::parameters_set set to false.

	@see RNSEncryptionParameters for more details on the parameters.
	*/
	class RNSContext
	{
	public:
		/**
		Creates an instance of RNSContext and sets the parameters according to the
		given encryption parameters. Optionally, the user can give a reference to a MemoryPoolHandle
		object to use a custom memory pool instead of the global memory pool (default).

		@param[in] parms The encryption parameters
		@param[in] pool The memory pool handle
		@see RNSEncryptionParameters for more details on the parameters.
		@see MemoryPoolHandle for more details on memory pool handles.
		*/
		RNSContext(const RNSEncryptionParameters &parms,
			const MemoryPoolHandle &pool = MemoryPoolHandle::acquire_global());

		/**
		Creates a copy of a RNSContext.

		@param[in] copy The RNSContext to copy from
		*/
		RNSContext(const RNSContext &copy) = default;

		/**
		Overwrites the RNSContext instance with the specified instance.

		@param[in] assign The RNSContext instance that should be assigned to the current instance
		*/
		RNSContext &operator =(const RNSContext &assign) = default;

		/**
		Creates a new RNSContext by moving an old one.

		@param[in] source The RNSContext to move from
		*/
		RNSContext(RNSContext &&source) = default;

		/**
		Moves an old RNSContext instance to the current one.

		@param[in] assign The RNSContext to move from
		*/
		RNSContext &operator =(RNSContext &&assign) = default;

		/**
		Returns the underlying encryption parameters in a RNSEncryptionParameters instance.

		@see RNSEncryptionParameters for more details on the parameters.
		*/
		const RNSEncryptionParameters &get_parms() const
		{
			return parms_;
		}

		/**
		Returns the set of qualifiers (as an instance of EncryptionParameterQualifiers) for the current
		encryption parameters. Note that to get an updated set of qualifiers it is necessary to call
		validate() after any change to the encryption parameters.

		@see EncryptionParameterQualifiers for more details on the qualifiers.
		*/
		inline RNSEncryptionParameterQualifiers get_qualifiers() const
		{
			return qualifiers_;
		}

		/**
		Returns a constant reference to the polynomial modulus parameter (represented by BigPoly).
		Note that the polynomial modulus directly determines the number of coefficients of encrypted
		polynomials, and the maximum number of coefficients for plaintext polynomials that are
		representable by the library.
		*/
		inline const BigPoly &poly_modulus() const
		{
			return parms_.poly_modulus();
		}

		/**
		Returns a constant reference to the coefficient modulus parameter (represented by a BigUInt).
		This value is computed as the products of small moduli and represents the actual value of the 
		coefficient modulus. Note that the coefficient modulus directly determines the number of
		bits-per-coefficient of encrypted polynomials and the maximum value allowed for plain_modulus()
		(which should be significantly smaller than coeff_modulus()).
		*/
		inline const BigUInt &coeff_modulus() const
		{
			return coeff_modulus_;
		}

		/**
		Returns a constant reference to the plaintext modulus parameter (represented by a BigUInt).
		Note that the plaintext modulus is one greater than the maximum value allowed for any
		plaintext coefficient that the library can encrypt or represent.
		*/
		inline const SmallModulus &plain_modulus() const
		{
			return parms_.plain_modulus();
		}

		/**
		Returns the standard deviation of normalized noise used during key generation and encryption.
		*/
		inline double noise_standard_deviation() const
		{
			return parms_.noise_standard_deviation();
		}

		/**
		Returns the maximum deviation of normalized noise used during key generation and encryption.
		*/
		inline double noise_max_deviation() const
		{
			return parms_.noise_max_deviation();
		}

		/**
		Returns the decomposition bit count which directly determines the number of evaluation keys
		required by the scheme. Smaller decomposition bit count reduces the accumulation of noise
		during multiplication operations, but can also significantly increase the time required
		to perform multiplication.
		*/
		inline int decomposition_bit_count() const
		{
			return parms_.decomposition_bit_count();
		}

		/**
		Returns the random number generator factory to use for encryption. By default, the random
		generator is set to UniformRandomGeneratorFactory::default_factory(). Setting this value
		allows a user to specify a custom random number generator source.
		*/
		inline UniformRandomGeneratorFactory *random_generator() const
		{
			return parms_.random_generator();
		}

		/**
		Returns the vector of coefficient moduli. This array is set by user as the input of coefficient modulus
		*/
		inline const std::vector<SmallModulus> coeff_mod_array() const
		{
			return parms_.coeff_mod_array();
		}

		/**
		Returns the base convertor object. The base convertor object includes all the required pre-computed tables
		for RNS computations. Both RNSDecryptor and RNSEvaluator need to access these tables during the computations.
		*/
		inline const util::BaseConvertor &base_convertor() const
		{
			return base_convertor_;
		}

	private:
		RNSEncryptionParameterQualifiers validate();

		MemoryPoolHandle pool_;

		RNSEncryptionParameters parms_;

		RNSEncryptionParameterQualifiers qualifiers_;

		util::BaseConvertor base_convertor_;

		std::vector<util::SmallNTTTables> small_ntt_tables_;

		util::SmallNTTTables plain_ntt_tables_;

		BigUInt coeff_modulus_;

		friend class RNSDecryptor;

		friend class RNSEncryptor;

		friend class RNSEvaluator;

		friend class RNSRecryptor;

		friend class PolyCRTBuilder;

		friend class RNSKeyGenerator;

        friend class RNSPolyCRTBuilder;
	};
}