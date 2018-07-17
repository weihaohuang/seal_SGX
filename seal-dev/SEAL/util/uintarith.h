#pragma once

#include <stdexcept>
#include <cstdint>
#include <functional>
#include "util/uintcore.h"
#include "util/mempool.h"
#include "util/modulus.h"
#include "util/defines.h"

namespace seal
{
    namespace util
    {
        unsigned char increment_uint(const std::uint64_t *operand, int uint64_count, std::uint64_t *result);

        unsigned char decrement_uint(const std::uint64_t *operand, int uint64_count, std::uint64_t *result);

        void negate_uint(const std::uint64_t *operand, int uint64_count, std::uint64_t *result);

        void left_shift_uint(const std::uint64_t *operand, int shift_amount, int uint64_count, std::uint64_t *result);

        void mul_by_shift(const uint64_t *operand, int shift_amount, int uint64_count, uint64_t *result);

        void right_shift_uint(const std::uint64_t *operand, int shift_amount, int uint64_count, std::uint64_t *result);

        void right_shift_sign_extend_uint(const std::uint64_t *operand, int shift_amount, int uint64_count, std::uint64_t *result);

        void half_round_up_uint(const std::uint64_t *operand, int uint64_count, std::uint64_t *result);

        void not_uint(const std::uint64_t *operand, int uint64_count, std::uint64_t *result);

        void and_uint_uint(const std::uint64_t *operand1, const std::uint64_t *operand2, int uint64_count, std::uint64_t *result);

        void or_uint_uint(const std::uint64_t *operand1, const std::uint64_t *operand2, int uint64_count, std::uint64_t *result);

        void xor_uint_uint(const std::uint64_t *operand1, const std::uint64_t *operand2, int uint64_count, std::uint64_t *result);

        inline unsigned char add_uint64_generic(std::uint64_t operand1, std::uint64_t operand2, unsigned char carry, std::uint64_t *result)
        {
#ifdef _DEBUG
            if (result == nullptr)
            {
                throw std::invalid_argument("result cannot be null");
            }
#endif
            operand1 += operand2;
            *result = operand1 + carry;
            return (operand1 < operand2) || (~operand1 < carry);

            //// Same as
            //uint64_t sum = operand1 + operand2;
            //*result = sum + carry;
            //return (~sum < carry) | (sum < operand1);
        }

        inline unsigned char add_uint64(std::uint64_t operand1, std::uint64_t operand2, unsigned char carry, std::uint64_t *result)
        {
            return ADD_CARRY_UINT64(operand1, operand2, carry, result);
        }

        unsigned char add_uint_uint(const std::uint64_t *operand1, int operand1_uint64_count, const std::uint64_t *operand2, int operand2_uint64_count, unsigned char carry, int result_uint64_count, std::uint64_t *result);

        unsigned char add_uint_uint(const std::uint64_t *operand1, const std::uint64_t *operand2, int uint64_count, std::uint64_t *result);

        inline unsigned char sub_uint64_generic(std::uint64_t operand1, std::uint64_t operand2, unsigned char borrow, std::uint64_t *result)
        {
#ifdef _DEBUG
            if (result == nullptr)
            {
                throw std::invalid_argument("result cannot be null");
            }
#endif
            uint64_t diff = operand1 - operand2;
            *result = diff - (borrow != 0);
            return (diff > operand1) || (diff < borrow);
        }

        inline unsigned char sub_uint64(std::uint64_t operand1, std::uint64_t operand2, unsigned char borrow, std::uint64_t *result)
        {
            return SUB_BORROW_UINT64(operand1, operand2, borrow, result);
        }

        unsigned char sub_uint_uint(const std::uint64_t *operand1, int operand1_uint64_count, const std::uint64_t *operand2, int operand2_uint64_count, unsigned char borrow, int result_uint64_count, std::uint64_t *result);

        unsigned char sub_uint_uint(const std::uint64_t *operand1, const std::uint64_t *operand2, int uint64_count, std::uint64_t *result);

        inline void multiply_uint64_generic(std::uint64_t operand1, std::uint64_t operand2, std::uint64_t *result128)
        {
#ifdef _DEBUG
            if (result128 == nullptr)
            {
                throw std::invalid_argument("result128 cannot be null");
            }
#endif
            uint64_t operand1_coeff_right = operand1 & 0x00000000FFFFFFFF;
            uint64_t operand2_coeff_right = operand2 & 0x00000000FFFFFFFF;
            operand1 >>= 32;
            operand2 >>= 32;

            uint64_t middle1 = operand1 * operand2_coeff_right;
            uint64_t middle;
            uint64_t left = operand1 * operand2
                + (static_cast<uint64_t>(add_uint64(middle1, operand2 * operand1_coeff_right, 0, &middle)) << 32);
            uint64_t right = operand1_coeff_right * operand2_coeff_right;
            uint64_t temp_sum = (right >> 32) + (middle & 0x00000000FFFFFFFF);
            
            result128[1] = left + (middle >> 32) + (temp_sum >> 32);
            result128[0] = (temp_sum << 32) | (right & 0x00000000FFFFFFFF);
        }

        inline void multiply_uint64(std::uint64_t operand1, std::uint64_t operand2, std::uint64_t *result128)
        {
            MULTIPLY_UINT64(operand1, operand2, result128);
        }

        inline void multiply_uint64_hw64_generic(std::uint64_t operand1, std::uint64_t operand2, std::uint64_t *hw64)
        {
#ifdef _DEBUG
            if (hw64 == nullptr)
            {
                throw std::invalid_argument("hw64 cannot be null");
            }
#endif
            uint64_t operand1_coeff_right = operand1 & 0x00000000FFFFFFFF;
            uint64_t operand2_coeff_right = operand2 & 0x00000000FFFFFFFF;
            operand1 >>= 32;
            operand2 >>= 32;

            uint64_t middle1 = operand1 * operand2_coeff_right;
            uint64_t middle;
            uint64_t left = operand1 * operand2
                + (static_cast<uint64_t>(add_uint64(middle1, operand2 * operand1_coeff_right, 0, &middle)) << 32);
            uint64_t right = operand1_coeff_right * operand2_coeff_right;
            uint64_t temp_sum = (right >> 32) + (middle & 0x00000000FFFFFFFF);

            *hw64 = left + (middle >> 32) + (temp_sum >> 32);
        }

        inline void multiply_uint64_hw64(std::uint64_t operand1, std::uint64_t operand2, std::uint64_t *hw64)
        {
            MULTIPLY_UINT64_HW64(operand1, operand2, hw64);
        }

        void multiply_uint_uint(const std::uint64_t *operand1, int operand1_uint64_count, const std::uint64_t *operand2, int operand2_uint64_count, int result_uint64_count, std::uint64_t *result);

        inline void multiply_uint_uint(const std::uint64_t *operand1, const std::uint64_t *operand2, int uint64_count, std::uint64_t *result)
        {
            multiply_uint_uint(operand1, uint64_count, operand2, uint64_count, uint64_count * 2, result);
        }

        void multiply_uint_uint64(const std::uint64_t *operand1, int operand1_uint64_count, std::uint64_t operand2, int result_uint64_count, std::uint64_t *result);

        inline void multiply_truncate_uint_uint(const std::uint64_t *operand1, const std::uint64_t *operand2, int uint64_count, std::uint64_t *result)
        {
            multiply_uint_uint(operand1, uint64_count, operand2, uint64_count, uint64_count, result);
        }

        void divide_uint_uint_inplace(std::uint64_t *numerator, const std::uint64_t *denominator, int uint64_count, std::uint64_t *quotient, MemoryPool &pool);

        void divide_uint_uint_inplace(std::uint64_t *numerator, const Modulus &denominator, int uint64_count, std::uint64_t *quotient, MemoryPool &pool);

        inline void divide_uint_uint(const std::uint64_t *numerator, const std::uint64_t *denominator, int uint64_count, std::uint64_t *quotient, std::uint64_t *remainder, MemoryPool &pool)
        {
            set_uint_uint(numerator, uint64_count, remainder);
            divide_uint_uint_inplace(remainder, denominator, uint64_count, quotient, pool);
        }
    }
}
