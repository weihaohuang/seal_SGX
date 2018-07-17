#pragma once

#include <string>
#include <iostream>
#include "bigpolyarray.h"
#include "encryptionparams.h"
#include "rnsencryptionparams.h"

namespace seal
{
    /**
    Class to store a ciphertext element. The Ciphertext class wraps an instance of the BigPolyArray 
    class, and stores a hash of the encryption parameters that it was created with. 

    @par Thread Safety
    In general, reading from Ciphertext is thread-safe as long as no other thread is concurrently
    mutating it. This is due to the underlying data structure storing the ciphertext not being
    thread-safe.

    @see Plaintext for the class that stores plaintexts.
    */
    class Ciphertext
    {
    public:
        /**
        Creates an empty ciphertext. No memory is allocated by this constructor.
        */
        Ciphertext() = default;

        /**
        Creates a new Ciphertext by copying an old one.

        @param[in] copy The Ciphertext to copy from
        */
        Ciphertext(const Ciphertext &copy) = default;

        /**
        Creates a new Ciphertext by moving an old one.

        @param[in] source The Ciphertext to move from
        */
        Ciphertext(Ciphertext &&source) = default;

        /**
        Copies an old Ciphertext to the current one.

        @param[in] assign The Ciphertext to copy from
        */
        Ciphertext &operator =(const Ciphertext &assign) = default;

        /**
        Moves an old Ciphertext to the current one.

        @param[in] assign The Ciphertext to move from
        */
        Ciphertext &operator =(Ciphertext &&assign) = default;

        /**
        Returns a constant reference to the underlying BigPolyArray.
        */
        inline const BigPolyArray &get_array() const
        {
            return ciphertext_array_;
        }

        /**
        Returns the size of the ciphertext.
        */
        inline int size() const
        {
            return ciphertext_array_.size();
        }

        /**
        Saves the Ciphertext to an output stream. The output is in binary format and not 
        human-readable. The output stream must have the "binary" flag set.

        @param[in] stream The stream to save the Ciphertext to
        @see load() to load a saved Ciphertext.
        */
        void save(std::ostream &stream) const
        {
            stream.write(reinterpret_cast<const char*>(&hash_block_), sizeof(EncryptionParameters::hash_block_type));
            stream.write(reinterpret_cast<const char*>(&rns_hash_block_), sizeof(EncryptionParameters::hash_block_type));
            ciphertext_array_.save(stream);
        }

        /**
        Loads a Ciphertext from an input stream overwriting the current Ciphertext.

        @param[in] stream The stream to load the Ciphertext from
        @see save() to save a Ciphertext.
        */
        void load(std::istream &stream)
        {
            stream.read(reinterpret_cast<char*>(&hash_block_), sizeof(EncryptionParameters::hash_block_type));
            stream.read(reinterpret_cast<char*>(&rns_hash_block_), sizeof(EncryptionParameters::hash_block_type));
            ciphertext_array_.load(stream);
        }

        /**
        Returns a constant reference to the hash block.

        @see EncryptionParameters for more information about the hash block.
        */
        inline const EncryptionParameters::hash_block_type &get_hash_block() const
        {
            return hash_block_;
        }

        /**
        Enables access to private members of seal::Ciphertext for .NET wrapper.
        */
        struct CiphertextPrivateHelper;

        
        inline const RNSEncryptionParameters::hash_block_type &get_rns_hash_block() const
        {
            return rns_hash_block_;
        }


		inline const bool is_ntt_form() const {
			return in_ntt_form;
		}

		inline void set_ntt_form(bool form) {
			in_ntt_form = form;
		}

    private:


		bool in_ntt_form; 

        /**
        Returns a reference to the underlying BigPolyArray. The user should never have
        a reason to modify the ciphertext by hand.
        */
        inline BigPolyArray &get_mutable_array()
        {
            return ciphertext_array_;
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

        BigPolyArray ciphertext_array_;

        EncryptionParameters::hash_block_type hash_block_;

        RNSEncryptionParameters::hash_block_type rns_hash_block_;

        friend class Decryptor;

        friend class Encryptor;

        friend class Evaluator;

        friend class RNSEncryptor;

        friend class RNSDecryptor;

        friend class RNSEvaluator;

		friend class RNSRecryptor;
    };
}