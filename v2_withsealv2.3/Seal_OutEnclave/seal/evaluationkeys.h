#pragma once

#include <iostream>
#include <vector>
#include "seal/ciphertext.h"
#include "seal/encryptionparams.h"

namespace seal
{
    /**
    Class to store evaluation keys. An evaluation key has type std::vector<Ciphertext>.
    An instance of the EvaluationKeys class stores internally an std::vector of evaluation
    keys.
    
    @par Relinearization
    Concretely, an evaluation key corresponding to a power K of the secret key can be used
    in the relinearization operation to change a ciphertext of size K+1 to size K. Recall
    that the smallest possible size for a ciphertext is 2, so the first evaluation key is
    corresponds to the square of the secret key. The second evaluation key corresponds to
    the cube of the secret key, and so on. For example, to relinearize a ciphertext of size
    7 back to size 2, one would need 5 evaluation keys, although it is hard to imagine 
    a situation where it makes sense to have size 7 ciphertexts, as operating on such objects
    would be very slow. Most commonly only one evaluation key is needed, and relinearization
    is performed after every multiplication.

    @par Decomposition Bit Count
    Decomposition bit count (dbc) is a parameter that describes a performance trade-off in
    the relinearization process. Namely, in the relinearization process the polynomials in 
    the ciphertexts (with large coefficients) get decomposed into a smaller base 2^dbc,
    coefficient-wise. Each of the decomposition factors corresponds to a piece of data in
    the evaluation key, so the smaller the dbc is, the larger the evaluation keys are.
    Moreover, a smaller dbc results in less invariant noise budget being consumed in the
    relinearization process. However, using a large dbc is much faster, and often one 
    would want to optimize the dbc to be as large as possible for performance. The dbc is 
    upper-bounded by the value of 60, and lower-bounded by the value of 1.

    @par Thread Safety
    In general, reading from EvaluationKeys is thread-safe as long as no other thread is
    concurrently mutating it. This is due to the underlying data structure storing the 
    evaluation keys not being thread-safe.

    @see SecretKey for the class that stores the secret key.
    @see PublicKey for the class that stores the public key.
    @see GaloisKeys for the class that stores the Galois keys.
    @see KeyGenerator for the class that generates the evaluation keys.
    */
    class EvaluationKeys
    {
    public:
        /**
        Creates an empty set of evaluation keys.
        */
        EvaluationKeys() = default;

        /**
        Creates a new EvaluationKeys instance by copying a given instance.

        @param[in] copy The EvaluationKeys to copy from
        */
        EvaluationKeys(const EvaluationKeys &copy) = default;

        /**
        Creates a new EvaluationKeys instance by moving a given instance.

        @param[in] source The EvaluationKeys to move from
        */
        EvaluationKeys(EvaluationKeys &&source) = default;

        /**
        Copies a given EvaluationKeys instance to the current one.

        @param[in] assign The EvaluationKeys to copy from
        */
        EvaluationKeys &operator =(const EvaluationKeys &assign) = default;

        /**
        Moves a given EvaluationKeys instance to the current one.

        @param[in] assign The EvaluationKeys to move from
        */
        EvaluationKeys &operator =(EvaluationKeys &&assign) = default;

        /**
        Returns the current number of evaluation keys.
        */
        inline int size() const
        {
            return keys_.size();
        }

        /**
        Returns the decomposition bit count.
        */
        inline int decomposition_bit_count() const
        {
            return decomposition_bit_count_;
        }

        /**
        Returns a constant reference to the evaluation keys data.
        */
        inline const std::vector<std::vector<Ciphertext> > &data() const
        {
            return keys_;
        }

        /**
        Returns a constant reference to an evaluation key. The returned evaluation key 
        corresponds to the given power of the secret key.

        @param[in] key_power The power of the secret key
        @throw std::invalid_argument if the key corresponding to key_power does not exist
        */
        inline const std::vector<Ciphertext> &key(int key_power) const
        {
            if (!has_key(key_power))
            {
                throw std::invalid_argument("requested key does not exist");
            }
            return keys_[key_power - 2];
        }

        /**
        Returns whether an evaluation key corresponding to a given power of the secret key 
        exists.

        @param[in] key_power The power of the secret key
        */
        inline bool has_key(int key_power) const
        {
            return (key_power >= 2) && (keys_.size() >= (key_power - 1));
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
        Saves the EvaluationKeys instance to an output stream. The output is in binary format 
        and not human-readable. The output stream must have the "binary" flag set.

        @param[in] stream The stream to save the EvaluationKeys to
        @see load() to load a saved EvaluationKeys instance.
        */
        void save(std::ostream &stream) const
        {
            // Save the hash block
            stream.write(reinterpret_cast<const char*>(&hash_block_), sizeof(EncryptionParameters::hash_block_type));
    
            // Save the decomposition bit count
            int32_t decomposition_bit_count32 = static_cast<int32_t>(decomposition_bit_count_);
            stream.write(reinterpret_cast<const char*>(&decomposition_bit_count32), sizeof(int32_t));
    
            // Save the size of keys_
            int32_t keys_dim1 = static_cast<int32_t>(keys_.size());
            stream.write(reinterpret_cast<const char*>(&keys_dim1), sizeof(int32_t));
    
            // Now loop again over keys_dim1
            for (int32_t index = 0; index < keys_dim1; index++)
            {
                // Save second dimension of keys_
                int32_t keys_dim2 = keys_[index].size();
                stream.write(reinterpret_cast<const char*>(&keys_dim2), sizeof(int32_t));
    
                // Loop over keys_dim2 and save all (or none)
                for (int32_t j = 0; j < keys_dim2; j++)
                {
                    // Save the key
                    keys_[index][j].save(stream);
                }
            }
        }

        /**
        Loads an EvaluationKeys instance from an input stream overwriting the current
        EvaluationKeys instance.

        @param[in] stream The stream to load the EvaluationKeys instance from
        @see save() to save an EvaluationKeys instance.
        */
        void load(std::istream &stream)
        {
            // Clear current keys
            keys_.clear();
    
            // Read the hash block
            stream.read(reinterpret_cast<char*>(&hash_block_), sizeof(EncryptionParameters::hash_block_type));
    
            // Read the decomposition_bit_count
            int32_t decomposition_bit_count32 = 0;
            stream.read(reinterpret_cast<char*>(&decomposition_bit_count32), sizeof(int32_t));
            decomposition_bit_count_ = decomposition_bit_count32;
    
            // Read in the size of keys_
            int32_t keys_dim1 = 0;
            stream.read(reinterpret_cast<char*>(&keys_dim1), sizeof(int32_t));
    
            // Resize first dimension of keys_
            keys_.resize(keys_dim1);
    
            // Loop over the first dimension of keys_
            for (int32_t index = 0; index < keys_dim1; index++)
            {
                // Read the size of the second dimension
                int32_t keys_dim2 = 0;
                stream.read(reinterpret_cast<char*>(&keys_dim2), sizeof(int32_t));
    
                // Resize
                keys_[index].resize(keys_dim2);
                for (int32_t j = 0; j < keys_dim2; j++)
                {
                    keys_[index][j].load(stream);
                }
            }
        }

        /**
        Enables access to private members of seal::EvaluationKeys for .NET wrapper.
        */
        struct EvaluationKeysPrivateHelper;

    private:
        /**
        Returns a reference to the vector of evaluation keys. The user should never have 
        a reason to modify the evaluation keys by hand.
        */
        inline std::vector<std::vector<Ciphertext> > &mutable_data()
        {
            return keys_;
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

        /**
        The vector of evaluation keys.
        */
        std::vector<std::vector<Ciphertext> > keys_;

        int decomposition_bit_count_ = 0;

        friend class KeyGenerator;

        friend class Evaluator;

        friend class KeyGenerator;
    };
}
