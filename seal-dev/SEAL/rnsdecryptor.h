#pragma once

#include <utility>
#include "rnsencryptionparams.h"
#include "rnscontext.h"
#include "util/modulus.h"
#include "util/polymodulus.h"
#include "util/smallntt.h"
#include "memorypoolhandle.h"
#include "ciphertext.h"
#include "plaintext.h"
#include "secretkey.h"
#include "util/baseconvertor.h"
#include "smallmodulus.h"

namespace seal
{
    /**
    Decrypts Ciphertext objects into Plaintext objects. Constructing an RNSDecryptor requires the encryption
    parameters (set through a RNSContext object) and the secret key. The public and evaluation keys
    are not needed for decryption.
    */
    class RNSDecryptor
    {
    public:
        /**
        Creates an RNSDecryptor instance initialized with the specified RNSContext and secret key.
        Optionally, the user can give a reference to a MemoryPoolHandle object to use a custom memory
        pool instead of the global memory pool (default).

        @param[in] context The RNSContext
        @param[in] secret_key The secret key
        @param[in] pool The memory pool handle
        @throws std::invalid_argument if encryption parameters or secret key are not valid
        @see RNSContext for more details on valid encryption parameters.
        @see MemoryPoolHandle for more details on memory pool handles.
        */
        RNSDecryptor(const RNSContext &context, const SecretKey &secret_key,
            const MemoryPoolHandle &pool = MemoryPoolHandle::acquire_global());

        /**
        Creates a copy of a RNSDecryptor.

        @param[in] copy The RNSDecryptor to copy from
        */
        RNSDecryptor(const RNSDecryptor &copy);

        /**
        Creates a new RNSDecryptor by moving an old one.

        @param[in] source The RNSDecryptor to move from
        */
        RNSDecryptor(RNSDecryptor &&source) = default;

        /*
        Decrypts a Ciphertext and stores the result in the destination parameter.

        @param[in] encrypted The ciphertext to decrypt
        @param[out] destination The plaintext to overwrite with the decrypted ciphertext
        @throws std::invalid_argument if the ciphertext is not a valid ciphertext for the encryption parameters
        */
        void rns_decrypt(const Ciphertext &encrypted, Plaintext &destination);

        /**
        Decrypts a Ciphertext and returns the result.

        @param[in] encrypted The ciphertext to decrypt
        @throws std::invalid_argument if the ciphertext is not a valid ciphertext for the encryption parameters
        */
        Plaintext rns_decrypt(const Ciphertext &encrypted)
        {
            Plaintext result;
            rns_decrypt(encrypted, result);
            return result;
        }

        /*
        Computes the invariant noise budget (in bits) of a ciphertext. The invariant noise budget
        measures the amount of room there is for the noise to grow while ensuring correct decryptions.

        @par Invariant Noise Budget
        The invariant noise polynomial of a ciphertext is a rational coefficient polynomial, such that
        a ciphertext decrypts correctly as long as the coefficients of the invariant noise polynomial are
        of absolute value less than 1/2. Thus, we call the infinity-norm of the invariant noise polynomial
        the invariant noise, and for correct decryption require it to be less than 1/2. If v denotes the
        invariant noise, we define the invariant noise budget as -log2(2v). Thus, the invariant noise budget
        starts from some initial value, which depends on the encryption parameters, and decreases to 0 when
        computations are performed. When the budget reaches 0, the ciphertext becomes too noisy to decrypt
        correctly.

        @param[in] encrypted The ciphertext
        @throws std::invalid_argument if encrypted is not a valid ciphertext for the encryption parameters
        */
        int invariant_noise_budget(const Ciphertext &encrypted);

    private:
        RNSDecryptor &operator =(const RNSDecryptor &assign) = delete;

        RNSDecryptor &operator =(RNSDecryptor &&assign) = delete;

        void compute_secret_key_array(int max_power);

        void rns_compose(std::uint64_t *value);

        MemoryPoolHandle pool_;

        RNSEncryptionParameters parms_;

        RNSEncryptionParameterQualifiers qualifiers_;

        util::BaseConvertor base_convertor_;

        std::vector<util::SmallNTTTables> small_ntt_tables_;

        std::uint64_t inv_gamma_mod_plain_modulus_;

        util::Pointer coeff_products_array_;

        util::Pointer secret_key_;

        util::Pointer product_modulus_;

        util::Modulus mod_;

        util::PolyModulus polymod_;

        std::vector<SmallModulus> coeff_mod_array_;

        std::vector<SmallModulus> plain_gamma_array_;

        std::vector<std::uint64_t> plain_gamma_product_;

        std::vector<std::uint64_t> inv_coeff_mod_coeff_array_;

        std::vector<std::uint64_t> neg_inv_coeff_mod_;

        BigPolyArray secret_key_array_;
    };
}