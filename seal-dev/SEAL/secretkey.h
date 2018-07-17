#pragma once

#include "bigpoly.h"
#include "encryptionparams.h"
#include "rnsencryptionparams.h"
#include <iostream>

namespace seal
{
    /**
    Class to store a secret key. Internally, the secret key is represented by a BigPoly object,
    and is created by KeyGenerator. The secret key can be saved and loaded from a stream with the
    save() and load() functions.

    @par Thread Safety
    In general, reading from SecretKey is thread-safe as long as no other thread is concurrently
    mutating it. This is due to the underlying data structure storing the secret key not being
    thread-safe.

    @see PublicKey for the class that stores the public key.
    @see EvaluationKeys for the class that stores the evaluation keys.
    @see KeyGenerator for the class that generates the secret key.
    */
    class SecretKey
    {
    public:
        /**
        Creates an empty secret key. No memory is allocated by this constructor.
        */
        SecretKey() = default;

        /**
        Creates a new SecretKey by copying an old one.

        @param[in] copy The SecretKey to copy from
        */
        SecretKey(const SecretKey &copy) = default;

        /**
        Creates a new SecretKey by moving an old one.

        @param[in] source The SecretKey to move from
        */
        SecretKey(SecretKey &&source) = default;

        /**
        Copies an old SecretKey to the current one.

        @param[in] assign The SecretKey to copy from
        */
        SecretKey &operator =(const SecretKey &assign)
        {
            sk_poly_.duplicate_from(assign.sk_poly_);
            hash_block_ = assign.hash_block_;
            rns_hash_block_ = assign.rns_hash_block_;
            return *this;
        }

        /**
        Moves an old SecretKey to the current one.

        @param[in] assign The SecretKey to move from
        */
        SecretKey &operator =(SecretKey &&assign) = default;

        /**
        Returns a constant reference to the underlying BigPolyArray.
        */
        inline const BigPoly &get_poly() const
        {
            return sk_poly_;
        }

        /**
        Saves the SecretKey to an output stream. The output is in binary format and not human-readable.
        The output stream must have the "binary" flag set.

        @param[in] stream The stream to save the SecretKey to
        @see load() to load a saved SecretKey.
        */
        void save(std::ostream &stream) const
        {
            stream.write(reinterpret_cast<const char*>(&hash_block_), sizeof(EncryptionParameters::hash_block_type));
            stream.write(reinterpret_cast<const char*>(&rns_hash_block_), sizeof(EncryptionParameters::hash_block_type));
            sk_poly_.save(stream);
        }

        /**
        Loads a SecretKey from an input stream overwriting the current SecretKey.

        @param[in] stream The stream to load the SecretKey from
        @see save() to save a SecretKey.
        */
        void load(std::istream &stream)
        {
            stream.read(reinterpret_cast<char*>(&hash_block_), sizeof(EncryptionParameters::hash_block_type));
            stream.read(reinterpret_cast<char*>(&rns_hash_block_), sizeof(EncryptionParameters::hash_block_type));
            sk_poly_.load(stream);
        }

        /**
        Returns a constant reference to the hash block.

        @see EncryptionParameters for more information about the hash block.
        */
        inline const EncryptionParameters::hash_block_type &get_hash_block() const
        {
            return hash_block_;
        }

        inline const RNSEncryptionParameters::hash_block_type &get_rns_hash_block() const
        {
            return rns_hash_block_;
        }

        /**
        Enables access to private members of seal::SecretKey for .NET wrapper.
        */
        struct SecretKeyPrivateHelper;

    private:
        /**
        Returns a reference to the underlying BigPoly. The user should never have
        a reason to modify the secret key by hand.
        */
        inline BigPoly &get_mutable_poly()
        {
            return sk_poly_;
        }

        /**
        Returns a reference to the hash block. The user should never have a reason to
        modify the hash block by hand.

        @see EncryptionParameters for more information about the hash block.
        */
        inline EncryptionParameters::hash_block_type &get_mutable_hash_block()
        {
            return hash_block_;
        }

        inline RNSEncryptionParameters::hash_block_type &get_rns_mutable_hash_block()
        {
            return rns_hash_block_;
        }

        BigPoly sk_poly_;

        EncryptionParameters::hash_block_type hash_block_;

        RNSEncryptionParameters::hash_block_type rns_hash_block_;

        friend class KeyGenerator;

        friend class Decryptor;

        friend class RNSKeyGenerator;

        friend class RNSDecryptor;


		friend class RNSRecryptor;
    };
}