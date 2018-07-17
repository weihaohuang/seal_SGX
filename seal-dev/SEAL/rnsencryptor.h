#pragma once

#include <memory>
#include <utility>
#include "rnsencryptionparams.h"
#include "util/modulus.h"
#include "util/polymodulus.h"
#include "bigpolyarray.h"
#include "plaintext.h"
#include "ciphertext.h"
#include "util/ntt.h"
#include "memorypoolhandle.h"
#include "rnscontext.h"
#include "util/smallntt.h"
#include "publickey.h"

namespace seal
{
    /**
    Encrypts Plaintext objects into Ciphertext objects in RNS form. Constructing an Encryptor requires
    the encryption parameters (set through an SEALContext object) and the public
    key. The secret and evaluation keys are not needed for encryption.
    */
    class RNSEncryptor
    {
    public:
        /**
        Creates an Encryptor instances initialized with the specified SEALContext
        and public key. Optionally, the user can give a reference to a MemoryPoolHandle
        object to use a custom memory pool instead of the global memory pool (default).

        @param[in] context The SEALContext
        @param[in] public_key The public key
        @param[in] pool The memory pool handle
        @throws std::invalid_argument if encryption parameters or public key are not valid
        @see SEALContext for more details on valid encryption parameters.
        @see MemoryPoolHandle for more details on memory pool handles.
        */
        RNSEncryptor(const RNSContext &context, const PublicKey &public_key,
            const MemoryPoolHandle &pool = MemoryPoolHandle::acquire_global());

        /**
        Creates a copy of a Encryptor.

        @param[in] copy The Encryptor to copy from
        */
        RNSEncryptor(const RNSEncryptor &copy);

        /**
        Creates a new Encryptor by moving an old one.

        @param[in] source The Encryptor to move from
        */
        RNSEncryptor(RNSEncryptor &&source) = default;

        /**
        Encrypts a plaintext and stores the result in the destination parameter. The
        destination parameter is resized if and only if its coefficient count or
        coefficient bit count does not match the encryption parameters. The plaintext
        polynomial must have a significant coefficient count smaller than the coefficient
        count specified by the encryption parameters, and with coefficient values less-than
        the plain modulus (SEALContext::plain_modulus()).

        @param[in] plain The plaintext to encrypt
        @param[out] destination The ciphertext to overwrite with the encrypted plaintext
        @throws std::invalid_argument if the plaintext polynomial's significant coefficient
        count or coefficient values are too large to represent with the encryption parameters
        */
        void rns_encrypt(const Plaintext &plain, Ciphertext &destination);

        /**
        Encrypts a plaintext and returns the result. The plaintext polynomial must have
        a significant coefficient count smaller than the coefficient count specified by
        the encryption parameters, and with coefficient values less-than the plain modulus
        (SEALContext::plain_modulus()).

        @param[in] plain The plaintext to encrypt
        @throws std::invalid_argument if the plaintext polynomial's significant coefficient
        count or coefficient values are too large to represent with the encryption parameters
        */
        Ciphertext rns_encrypt(const Plaintext &plain)
        {
            Ciphertext result;
            rns_encrypt(plain, result);
            return result;
        }

    private:
        RNSEncryptor &operator =(const RNSEncryptor &assign) = delete;

        RNSEncryptor &operator =(RNSEncryptor &&assign) = delete;

        void rns_preencrypt(const std::uint64_t *plain, int plain_coeff_count, int plain_coeff_uint64_count, std::uint64_t *destination);

        void set_poly_coeffs_normal_rns(std::uint64_t *poly, UniformRandomGenerator *random) const;

        void set_poly_coeffs_zero_one_negone_rns(uint64_t *poly, UniformRandomGenerator *random) const;

        void set_poly_coeffs_zero_one_rns(uint64_t *poly, UniformRandomGenerator *random) const;

        void rns_decompose(uint64_t *value);

        MemoryPoolHandle pool_;

        RNSEncryptionParameters parms_;

        RNSEncryptionParameterQualifiers qualifiers_;

        std::vector<util::SmallNTTTables> small_ntt_tables_;

        util::Pointer upper_half_threshold_;

        util::Pointer upper_half_increment_;

        util::Pointer coeff_div_plain_modulus_;

        util::Pointer public_key_;

        util::PolyModulus polymod_;

        std::vector<SmallModulus> coeff_mod_array_;
    };
}