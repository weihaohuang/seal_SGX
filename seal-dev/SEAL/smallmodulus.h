#pragma once

#include <cstdint>
#include "util/mempool.h"
#include <math.h>
#include <iostream>
#include <array>

namespace seal
{
    /**
    This class is used to represent small integers moduli. In particular, the encryption
    parameter plain_modulus is represented by an instance of SmallModulus, and the encryption 
    parameter coeff_modulus is represented by one or more instances of SmallModulus.

    @par Barrett reduction
    The purpose of this class is to store the value of a pre-computation required by Barrett 
    reduction. The class has a constructor that sets the value of the modulus to the given 
    std::uint64_t, and performs the pre-computation. The pre-computation is the constant ratio 
    given by floor(4^k / m), where k is the smallest integer where 2^k > m and m is the modulus. 
    For Barrett reduction and the pre-computation to work correctly, the value of the modulus 
    is restricted to be at most 63 bits, and cannot be exactly 2^62. However, for convenience,
    SmallModulus also has a constructor that takes no parameters, and sets its value, the 
    Barrett ratio, and the bit count to zero.

    @par Thread Safety
    In general, reading from SmallModulus is thread-safe as long as no other thread is 
    concurrently mutating it.

    @see EncryptionParameters for a description of the encryption parameters.
    */
    class SmallModulus
    {
    public:
        /**
        Creates a SmallModulus instance with given value represented by std::uint64_t. 

        @param[in] value The integer modulus
        @throws std::invalid_argument if value is exactly 2^62, or 64 bits long
        */
        SmallModulus(std::uint64_t value = 0);

        /**
        Creates a new SmallModulus by copying an old one.

        @param[in] copy The SmallModulus to copy from
        */
        SmallModulus(const SmallModulus &copy) = default;

        /**
        Creates a new SmallModulus by copying an old one.

        @param[in] source The SmallModulus to move from
        */
        SmallModulus(SmallModulus &&source) = default;

        /**
        Copies an old SmallModulus to the current one.

        @param[in] assign The SmallModulus to copy from
        */
        SmallModulus &operator =(const SmallModulus &assign) = default;

        /**
        Sets the value of the SmallModulus to a given one represented by std::uint64_t.

        @param[in] value The integer modulus
        @throws std::invalid_argument if value is exactly 2^62, or 64 bits long
        */
        inline SmallModulus &operator =(std::uint64_t value)
        {
            set_value(value);
            return *this;
        }

        /**
        Moves an old SmallModulus to the current one.

        @param[in] assign The SmallModulus to move from
        */
        SmallModulus &operator =(SmallModulus &&assign) = default;

        /**
        Returns the significant bit count of the value of the current SmallModulus.
        */
        inline int bit_count() const
        {
            return bit_count_;
        }

        /**
        Returns the length in std::uint64_t of the value of the current SmallModulus.
        */
        inline int uint64_count() const
        {
            return uint64_count_;
        }

        /**
        Returns a constant pointer to the value of the current SmallModulus.
        */
        inline const uint64_t *pointer() const
        {
            return &value_;
        }

        /**
        Returns the value of the current SmallModulus.
        */
        inline std::uint64_t value() const
        {
            return value_;
        }

        /**
        Returns the Barrett ratio computed for the value of the current SmallModulus.
        */
        inline const std::array<std::uint64_t, 3> &const_ratio() const
        {
            return const_ratio_;
        }

        /**
        Returns whether or not the value of the current SmallModulus is zero.
        */
        inline bool is_zero() const
        {
            return value_ == 0;
        }

        /**
        Returns whether or not a SmallModulus is equal to a second SmallModulus. The input operands 
        are not modified.

        @param[in] compare The SmallModulus to compare against
        */
        inline bool operator ==(const SmallModulus &compare) const
        {
            return value_ == compare.value_;
        }

        /**
        Returns whether or not the value of a SmallModulus is equal to a given value represented
        by std::uint64_t. The input operands are not modified.

        @param[in] compare The value to compare against
        */
        inline bool operator ==(std::uint64_t compare) const
        {
            return value_ == compare;
        }

        /**
        Returns whether or not a SmallModulus is different from a second SmallModulus. The input operands
        are not modified.

        @param[in] compare The SmallModulus to compare against
        */
        inline bool operator !=(const SmallModulus &compare) const
        {
            return !(value_ == compare.value_);
        }

        /**
        Returns whether or not the value of a SmallModulus is different from a given value represented
        by std::uint64_t. The input operands are not modified.

        @param[in] compare The value to compare against
        */
        inline bool operator !=(std::uint64_t compare) const
        {
            return !(value_ == compare);
        }

        /**
        Returns the value of the current SmallModulus as a hexadecimal string.
        */
        std::string to_string() const;

        /**
        Returns the value of the current SmallModulus as a decimal string.
        */
        std::string to_dec_string() const;

        /**
        Saves the SmallModulus to an output stream. The full state of the modulus is serialized. 
        The output is in binary format and not human-readable. The output stream must have the 
        "binary" flag set.

        @param[in] stream The stream to save the SmallModulus to
        @see load() to load a saved SmallModulus.
        */
        void save(std::ostream &stream) const;

        /**
        Loads a SmallModulus from an input stream overwriting the current SmallModulus.

        @param[in] stream The stream to load the SmallModulus from
        @see save() to save a SmallModulus.
        */
        void load(std::istream &stream);

        /**
        Enables access to private members of seal::Ciphertext for .NET wrapper.
        */
        struct SmallModulusPrivateHelper;

    private:
        SmallModulus(std::uint64_t value, std::array<std::uint64_t, 3> const_ratio, int bit_count, int uint64_count) :
            value_(value), const_ratio_(const_ratio), bit_count_(bit_count), uint64_count_(uint64_count)
        {
        }

        void set_value(std::uint64_t value);

        std::uint64_t value_;

        std::array<std::uint64_t, 3> const_ratio_;

        int bit_count_;

        int uint64_count_;
    };
}
