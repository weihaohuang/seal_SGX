#pragma once

#include <cstdint>
#include "util/mempool.h"
#include "util/modulus.h"
#include "smallmodulus.h"

namespace seal
{
    namespace util
    {
        void small_modulo_uint_inplace(std::uint64_t *value, int value_uint64_count, const SmallModulus &smallmodulus, MemoryPool &pool);

        void small_modulo_uint(const std::uint64_t *value, int value_uint64_count, const SmallModulus &smallmodulus, std::uint64_t *result, MemoryPool &pool);

        void increment_uint_smallmod(const std::uint64_t *operand, const SmallModulus &smallmodulus, std::uint64_t *result);

        void decrement_uint_smallmod(const std::uint64_t *operand, const SmallModulus &smallmodulus, std::uint64_t *result);

        void negate_uint_smallmod(const std::uint64_t *operand, const SmallModulus &smallmodulus, std::uint64_t *result);

        void div2_uint_smallmod(const std::uint64_t *operand, const SmallModulus &smallmodulus, std::uint64_t *result);

        void add_uint_uint_smallmod(const std::uint64_t *operand1, const std::uint64_t *operand2, const SmallModulus &smallmodulus, std::uint64_t *result);

        void sub_uint_uint_smallmod(const std::uint64_t *operand1, const std::uint64_t *operand2, const SmallModulus &smallmodulus, std::uint64_t *result);

        void barrett_reduce_128(const std::uint64_t *input, const SmallModulus &smallmodulus, std::uint64_t *result);

        void multiply_uint64_smallmod(const std::uint64_t *operand1, const std::uint64_t *operand2, const SmallModulus &smallmodulus, std::uint64_t *result);

        bool is_primitive_root_smallmod(const std::uint64_t *root, std::uint64_t degree, const SmallModulus &prime_modulus);

        // Try to find a primitive degree-th root of unity modulo small prime modulus, where degree must be a power of two.
        bool try_primitive_root_smallmod(std::uint64_t degree, const SmallModulus &prime_modulus, MemoryPool &pool, std::uint64_t *destination);

        // Try to find the smallest (as integer) primitive degree-th root of unity modulo small prime modulus, where degree must be a power of two.
        bool try_minimal_primitive_root_smallmod(std::uint64_t degree, const SmallModulus &prime_modulus, MemoryPool &pool, std::uint64_t *destination);
    }
}
