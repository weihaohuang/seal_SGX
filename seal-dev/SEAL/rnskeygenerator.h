#pragma once

#include <memory>
#include <utility>
#include "rnsencryptionparams.h"
#include "rnscontext.h"
#include "util/modulus.h"
#include "util/polymodulus.h"
#include "rnsevaluationkeys.h"
#include "util/smallntt.h"
#include "memorypoolhandle.h"
#include "publickey.h"
#include "secretkey.h"

#include "rnsgaloiskeys.h"

namespace seal
{
    /**
    Generates matching secret key, public key, and evaluation keys for encryption, decryption,
    and evaluation functions.

    Constructing an RNSKeyGenerator requires the encryption parameters (set through a RNSContext object).
    Invoking the generate() function will generate a new secret key (which can be read from secret_key()),
    public key (which can be read from public_key()), and evaluation keys (which can be read from
    evaluation_keys()).

    @warning RNSKeyGenerator is not thread-safe and a separate instance is needed for each potentially
    concurrent usage.
    */
    class RNSKeyGenerator
    {
    public:
        /**
        Creates an RNSKeyGenerator instance initialized with the specified RNSContext. Optionally, the
        user can give a reference to a MemoryPoolHandle object to use a custom memory pool instead of
        the global memory pool (default).

        @param[in] context The RNSContext
        @param[in] pool The memory pool handle
        @throws std::invalid_argument if encryption parameters are not valid
        @see RNSContext for more details on valid encryption parameters.
        @see MemoryPoolHandle for more details on memory pool handles.
        */
        RNSKeyGenerator(const RNSContext &context, const MemoryPoolHandle &pool = MemoryPoolHandle::acquire_global());

        /**
        Creates an RNSKeyGenerator instance initialized with the specified RNSContext and specified
        previously generated keys. This can be used to increase the number of evaluation keys
        (using GenerateEvaluationKeys()) from what had earlier been generated. If no evaluation keys
        had been generated earlier, one can simply pass a newly created empty instance of EvaluationKeys
        to the function. Optionally, the user can give a reference to a MemoryPoolHandle object to use
        a custom memory pool instead of the global memory pool (default).

        @param[in] context The RNSContext
        @param[in] secret_key A previously generated secret key
        @param[in] public_key A previously generated public key
        @param[in] evaluation_keys A previously generated set of evaluation keys
        @param[in] pool The memory pool handle
        @throws std::invalid_argument if encryption parameters are not valid
        @throws std::invalid_argument if secret, public or evaluation keys are not valid
        @see RNSContext for more details on valid encryption parameters.
        @see MemoryPoolHandle for more details on memory pool handles.
        */
        RNSKeyGenerator(const RNSContext &context, const SecretKey &secret_key, const PublicKey &public_key, 
            const RNSEvaluationKeys &evaluation_keys, const MemoryPoolHandle &pool = MemoryPoolHandle::acquire_global());

        /**
        Generates new matching set of secret key, public key, and any number of evaluation keys.
        The number of evaluation keys that will be generated can be specified by the input
        parameter evaluation_keys_count, which defaults to 0.

        @param[in] count The number of evaluation keys to generate
        @see secret_key() to read generated secret key
        @see public_key() to read generated public key
        @see evaluation_keys() to read generated evaluation keys
        @throws std::invalid_argument if evaluation keys cannot be generated for specified encryption parameters
        */
        void generate(int evaluation_keys_count = 0, int hammingweight = 0, bool auto_galois_key_gen = false);

        /**
        Generates evaluation keys so that there are count many in total. Each key is added as a new entry to the std::vector of evaluation keys.
        This function is automatically called by generate() to generate evaluation keys, but can be later called by the user to increase
        the number of evaluation keys on top of what has already been generated. An error is thrown if the user tries to generate evaluation
        keys before a secret key and public key have been generated.

        @param[in] count The total number of evaluation keys to have been generated.
        @throws std::invalid_argument if count is less than 0
        @throws std::invalid_argument if evaluation keys cannot be generated for specified encryption parameters
        @throws std::logic_error if called before the secret key and public key have been generated
        @see evaluation_keys() to read generated evaluation keys.
        @see generated() to see if the secret and public keys have already been generated.
        */
        void generate_rns_evaluation_keys(int count);

        /**
        Returns true or false depending on whether secret key and public key have been generated.
        */
        bool is_generated() const
        {
            return generated_;
        }

        /**
        Returns the generated secret key after a generate() invocation.
        @throws std::logic_error if encryption keys have not been generated
        */
        const SecretKey &secret_key() const;

        /**
        Returns the generated public key after a generate() invocation.
        @throws std::logic_error if encryption keys have not been generated
        */
        const PublicKey &public_key() const;

        /**
        Returns evaluation keys after a generate_evaluation_keys() or generate() invocation.
        @throws std::logic_error if encryption keys have not been generated
        @throws std::logic_error if evaluation keys have not been generated
        */
        const RNSEvaluationKeys &evaluation_keys() const;

		//////////////////////////////////////////////////////////////

		void generate(const SecretKey &secret_key, int evaluation_keys_count = 0);
		
		void generate_rns_galois_keys(std::uint64_t galois_elt);

		void generate_rns_galois_keys(std::uint64_t galois_elt, RNSGaloisKeys &galois_keys);

		void generate_rns_galois_keys_logn(RNSGaloisKeys &galois_keys);

		const RNSGaloisKeys &galois_keys() const;

    private:
        RNSKeyGenerator(const RNSKeyGenerator &copy) = delete;

        RNSKeyGenerator &operator =(const RNSKeyGenerator &assign) = delete;

        RNSKeyGenerator(RNSKeyGenerator &&source) = delete;

        RNSKeyGenerator &operator =(RNSKeyGenerator &&assign) = delete;

        void set_poly_coeffs_zero_one_negone_rns(std::uint64_t *poly, UniformRandomGenerator *random, int hammingweight) const;

        void set_poly_coeffs_normal_rns(std::uint64_t *poly, UniformRandomGenerator *random) const;

        void set_poly_coeffs_uniform_rns(std::uint64_t *poly, UniformRandomGenerator *random);

        void compute_rns_secret_key_array(int max_power);

        void populate_rns_evaluation_factors();

        MemoryPoolHandle pool_;

        RNSEncryptionParameters parms_;

        RNSEncryptionParameterQualifiers qualifiers_;

        std::vector<util::SmallNTTTables> small_ntt_tables_;

        RNSEvaluationKeys evaluation_keys_;

        PublicKey public_key_;

        SecretKey secret_key_;

        UniformRandomGeneratorFactory *random_generator_;

        util::PolyModulus polymod_;

        std::vector<std::vector<std::uint64_t> > rns_evaluation_factors_;

        BigPolyArray secret_key_array_;

        bool generated_;

        int coeff_mod_count_;

		//////////////////////////////////////////////////////////////

		RNSGaloisKeys galois_keys_;
    };
}
