#pragma once

#include "seal/bigpoly.h"
#include "seal/encryptionparams.h"
//#include <iostream>

namespace seal
{
    /**
    Class to store a secret key. Internally, the secret key is represented 
    by a BigPoly object, and is created by KeyGenerator.

    @par Thread Safety
    In general, reading from SecretKey is thread-safe as long as no other 
    thread is concurrently mutating it. This is due to the underlying data 
    structure storing the secret key not being thread-safe.

    @see KeyGenerator for the class that generates the secret key.
    @see PublicKey for the class that stores the public key.
    @see EvaluationKeys for the class that stores the evaluation keys.
    @see GaloisKeys for the class that stores the Galois keys.
    */
    class SecretKey
    {
    public:
        /**
        Creates an empty secret key.
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
        inline const BigPoly &data() const
        {
            return sk_poly_;
        }

        /**
        Saves the SecretKey to an output stream. The output is in binary format and 
        not human-readable. The output stream must have the "binary" flag set.

        @param[in] stream The stream to save the SecretKey to
        @see load() to load a saved SecretKey.
        */
//        void save(std::ostream &stream) const
//        {
//            stream.write(reinterpret_cast<const char*>(&hash_block_), 
//                sizeof(EncryptionParameters::hash_block_type));
 //           sk_poly_.save(stream);
 //       }
        char* save(int &buffer_size)
        {
            int poly_size;
            char *sk_buffer = sk_poly_.save(poly_size);
            
            buffer_size = sizeof(EncryptionParameters::hash_block_type) + poly_size; // not sure whether + 1 or not
            char *buffer = (char*)malloc(buffer_size);
            
            int offset = 0;
        		memcpy(buffer + offset, (char*)&hash_block_, sizeof(EncryptionParameters::hash_block_type));
        		offset += sizeof(EncryptionParameters::hash_block_type);
        		
        		memcpy(buffer + offset, sk_buffer, poly_size);
           
            delete [] sk_buffer;
        		return buffer;
        }

        /**
        Loads a SecretKey from an input stream overwriting the current SecretKey.

        @param[in] stream The stream to load the SecretKey from
        @see save() to save a SecretKey.
        */
 //       void load(std::istream &stream)
 //       {
  //          stream.read(reinterpret_cast<char*>(&hash_block_), 
  //              sizeof(EncryptionParameters::hash_block_type));
 //           sk_poly_.load(stream);
  //      }
        void load(char* buffer)
        {
            int offset = 0;
        		int32_t read_coeff_count = 0;
            memcpy(&hash_block_, buffer + offset, sizeof(EncryptionParameters::hash_block_type));
            offset += sizeof(EncryptionParameters::hash_block_type);
            
            sk_poly_.load(buffer + offset);
            
            //memcpy(plaintext_poly_.get(), buffer + offset, read_coeff_count * bytes_per_uint64);
        }
        
        /**
        Returns a constant reference to the hash block.

        @see EncryptionParameters for more information about the hash block.
        */
        inline const EncryptionParameters::hash_block_type &hash_block() const
        {
            return hash_block_;
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
        inline BigPoly &mutable_data()
        {
            return sk_poly_;
        }
#ifdef SEAL_EXPOSE_MUTABLE_HASH_BLOCK
    public:
#endif
        /**
        Returns a reference to the hash block. The user should normally never have
        a reason to modify the hash block by hand.

        @see EncryptionParameters for more information about the hash block.
        */
        inline EncryptionParameters::hash_block_type &mutable_hash_block()
        {
            return hash_block_;
        }
#ifdef SEAL_EXPOSE_MUTABLE_HASH_BLOCK
    private:
#endif
        // C++11 compatibility
        EncryptionParameters::hash_block_type hash_block_{ { 0 } };

        BigPoly sk_poly_;

        friend class KeyGenerator;

        friend class Decryptor;

        friend class KeyGenerator;

        friend class Decryptor;
    };
}
