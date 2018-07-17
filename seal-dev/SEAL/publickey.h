#pragma once

#include "bigpolyarray.h"
#include "encryptionparams.h"
#include "rnsencryptionparams.h"
#include <iostream>

namespace seal
{
    /**
    Class to store a public key. Internally, the public key is represented by a BigPolyArray object,
    and is created by KeyGenerator. The public key can be saved and loaded from a stream with the 
    save() and load() functions.

    @par Thread Safety
    In general, reading from PublicKey is thread-safe as long as no other thread is concurrently
    mutating it. This is due to the underlying data structure storing the public key not being
    thread-safe.

    @see SecretKey for the class that stores the secret key.
    @see EvaluationKeys for the class that stores the evaluation keys.
    @see KeyGenerator for the class that generates the public key.
    */
    class PublicKey
    {
    public:
        /**
        Creates an empty public key. No memory is allocated by this constructor.
        */
        PublicKey() = default;

        /**
        Creates a new PublicKey by copying an old one.

        @param[in] copy The PublicKey to copy from
        */
        PublicKey(const PublicKey &copy) = default;

        /**
        Creates a new PublicKey by moving an old one.

        @param[in] source The PublicKey to move from
        */
        PublicKey(PublicKey &&source) = default;

        /**
        Copies an old PublicKey to the current one.

        @param[in] assign The PublicKey to copy from
        */
        PublicKey &operator =(const PublicKey &assign) = default;

        /**
        Moves an old PublicKey to the current one.

        @param[in] assign The PublicKey to move from
        */
        PublicKey &operator =(PublicKey &&assign) = default;

        /**
        Returns a constant reference to the underlying BigPolyArray.
        */
        inline const BigPolyArray &get_array() const
        {
            return pk_array_;
        }

        /**
        Saves the PublicKey to an output stream. The output is in binary format and not 
        human-readable. The output stream must have the "binary" flag set.

        @param[in] stream The stream to save the PublicKey to
        @see load() to load a saved PublicKey.
        */
        void save(std::ostream &stream) const
        {
            stream.write(reinterpret_cast<const char*>(&hash_block_), sizeof(EncryptionParameters::hash_block_type));
            stream.write(reinterpret_cast<const char*>(&rns_hash_block_), sizeof(EncryptionParameters::hash_block_type));
            pk_array_.save(stream);
        }

        /**
        Loads a PublicKey from an input stream overwriting the current PublicKey.

        @param[in] stream The stream to load the PublicKey from
        @see save() to save a PublicKey.
        */
        void load(std::istream &stream)
        {
            stream.read(reinterpret_cast<char*>(&hash_block_), sizeof(EncryptionParameters::hash_block_type));
            stream.read(reinterpret_cast<char*>(&rns_hash_block_), sizeof(EncryptionParameters::hash_block_type));
            pk_array_.load(stream);
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
        Enables access to private members of seal::PublicKey for .NET wrapper.
        */
        struct PublicKeyPrivateHelper;

    private:
        /**
        Returns a reference to the underlying BigPolyArray. The user should never have
        a reason to modify the public key by hand.
        */
        inline BigPolyArray &get_mutable_array()
        {
            return pk_array_;
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

        BigPolyArray pk_array_;

        EncryptionParameters::hash_block_type hash_block_;

        RNSEncryptionParameters::hash_block_type rns_hash_block_;

        friend class KeyGenerator;

        friend class Encryptor;

        friend class RNSKeyGenerator;

        friend class RNSEncryptor;
    };
}