#include "smallmodulus.h"
#include "util/common.h"
#include "util/uintarith.h"
#include "memorypoolhandle.h"
#include <stdexcept>

using namespace seal::util;
using namespace std;

namespace seal
{
    SmallModulus::SmallModulus(std::uint64_t value)
    {
        set_value(value);
    }

    string SmallModulus::to_string() const
    {
        return uint_to_hex_string(&value_, uint64_count_);
    }

    string SmallModulus::to_dec_string() const
    {
        return uint_to_dec_string(&value_, uint64_count_, MemoryPoolHandle::acquire_global());
    }

    void SmallModulus::save(std::ostream &stream) const
    {
        stream.write(reinterpret_cast<const char*>(&bit_count_), sizeof(int32_t));
        stream.write(reinterpret_cast<const char*>(&uint64_count_), sizeof(int32_t));
        stream.write(reinterpret_cast<const char*>(&value_), bytes_per_uint64);
        stream.write(reinterpret_cast<const char*>(&const_ratio_), 3 * bytes_per_uint64);
    }

    void SmallModulus::load(std::istream &stream)
    {
        stream.read(reinterpret_cast<char*>(&bit_count_), sizeof(int32_t));
        stream.read(reinterpret_cast<char*>(&uint64_count_), sizeof(int32_t));
        stream.read(reinterpret_cast<char*>(&value_), bytes_per_uint64);
        stream.read(reinterpret_cast<char *>(&const_ratio_), 3 * bytes_per_uint64);
    }

    void SmallModulus::set_value(uint64_t value)
    {
        if (value == 0)
        {
            // Zero settings

            value_ = 0;
            const_ratio_ = { 0 };
            bit_count_ = 0;
            uint64_count_ = 1;
        }
        else if (value >> 63 != 0 || value == 0x4000000000000000)
        {
            throw std::invalid_argument("value can be at most 63 bits, and cannot be 2^62");
        }
        else
        {
            // All normal, compute const_ratio and set everything
            value_ = value;
            bit_count_ = get_significant_bit_count(value_);

            // Compute Barrett ratios for 64-bit words (barrett_reduce_128)
            uint64_t numerator[3]{ 0, 0, 1 };
            uint64_t tmp_coeff_mod[3]{ value_, 0, 0 };
            uint64_t remainder[3]{ 0 };
            uint64_t quotient[3]{ 0 };
            divide_uint_uint(numerator, tmp_coeff_mod, 3, quotient, remainder, MemoryPoolHandle::acquire_global());
            const_ratio_[0] = quotient[0];
            const_ratio_[1] = quotient[1];

            // We need the remainder for barrett_reduce_256
            const_ratio_[2] = remainder[0];

            uint64_count_ = 1;
        }
    }
}