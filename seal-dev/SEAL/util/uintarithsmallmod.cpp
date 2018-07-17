#include "util/uintcore.h"
#include "util/uintarith.h"
#include "util/uintarithmod.h"
#include "util/common.h"
#include "smallmodulus.h"
#include "util/uintextras.h"
#include <random>
#include "uintarithsmallmod.h"

using namespace std;

namespace seal
{
    namespace util
    {
        void small_modulo_uint_inplace(uint64_t *value, int value_uint64_count, const SmallModulus &smallmodulus, MemoryPool &pool)
        {
#ifdef _DEBUG
            if (value == nullptr && value_uint64_count > 0)
            {
                throw invalid_argument("value");
            }
            if (value_uint64_count < 0)
            {
                throw invalid_argument("value_uint64_count");
            }
#endif
            // Determine significant bits in value and modulus.
            int value_bits = get_significant_bit_count_uint(value, value_uint64_count);
            int modulus_bits = smallmodulus.bit_count();

            // If value has fewer bits than modulus, then done.
            if (value_bits < modulus_bits)
            {
                return;
            }

            // Only perform computation on non-zero uint64s.
            int uint64_count = divide_round_up(value_bits, bits_per_uint64);

            // If value is smaller, then done.
            if (value_bits == modulus_bits && *value < smallmodulus.value())
            {
                return;
            }

            // Handle fast case.
            if (uint64_count == 1)
            {
                *value %= smallmodulus.value();
                return;
            }

            Pointer shifted(allocate_uint(uint64_count, pool));

            // Store mutable copy of modulus.
            set_uint(smallmodulus.value(), uint64_count, shifted.get());

            // Create temporary space to store difference calculation.
            Pointer difference(allocate_uint(uint64_count, pool));

            // Shift modulus to bring MSB in alignment with MSB of value.
            int modulus_shift = value_bits - modulus_bits;
            left_shift_uint(shifted.get(), modulus_shift, uint64_count, shifted.get());
            modulus_bits += modulus_shift;

            // Perform bit-wise division algorithm.
            int remaining_shifts = modulus_shift;
            while (value_bits == modulus_bits)
            {
                // NOTE: MSBs of value and shifted modulus are aligned.

                // Even though MSB of value and modulus are aligned, still possible value < shifted_modulus.
                if (sub_uint_uint(value, shifted.get(), uint64_count, difference.get()))
                {
                    // value < shifted_modulus, so current quotient bit is zero and next one is definitely one.
                    if (remaining_shifts == 0)
                    {
                        // No shifts remain and value < modulus so done.
                        break;
                    }

                    // Effectively shift value left by 1 by instead adding value to difference (to prevent overflow in value).
                    add_uint_uint(difference.get(), value, uint64_count, difference.get());

                    // Adjust remaining shifts as a result of shifting value.
                    remaining_shifts--;
                }
                // Difference is the new value with modulus subtracted.

                // Determine amount to shift value to bring MSB in alignment with modulus.
                value_bits = get_significant_bit_count_uint(difference.get(), uint64_count);
                int value_shift = modulus_bits - value_bits;
                if (value_shift > remaining_shifts)
                {
                    // Clip the maximum shift to determine only the integer (as opposed to fractional) bits.
                    value_shift = remaining_shifts;
                }

                // Shift and update value.
                if (value_bits > 0)
                {
                    left_shift_uint(difference.get(), value_shift, uint64_count, value);
                    value_bits += value_shift;
                }
                else
                {
                    // Value is zero so no need to shift, just set to zero.
                    set_zero_uint(uint64_count, value);
                }

                // Adjust remaining shifts as a result of shifting value.
                remaining_shifts -= value_shift;
            }

            // Correct value (which is also the remainder) for shifting of modulus.
            right_shift_uint(value, modulus_shift, uint64_count, value);
        }

        void small_modulo_uint(const uint64_t *value, int value_uint64_count, const SmallModulus &smallmodulus, uint64_t *result, MemoryPool &pool)
        {
#ifdef _DEBUG
            if (value == nullptr && value_uint64_count > 0)
            {
                throw invalid_argument("value");
            }
            if (value_uint64_count <= 0)
            {
                throw invalid_argument("value_uint64_count");
            }
            if (result == nullptr && value_uint64_count > 0)
            {
                throw invalid_argument("result");
            }
#endif
            if (value_uint64_count == 1)
            {
                // If value < smallmodulus no operation is needed
                if (*value < smallmodulus.value())	
                {
                    *result = *value;
                }
                else
                {
                    *result = *value % smallmodulus.value();
                }
                return;
            }

            Pointer value_copy_anchor;
            uint64_t *value_copy;
            value_copy_anchor = allocate_uint(value_uint64_count, pool);
            value_copy = value_copy_anchor.get();
            
            set_uint_uint(value, value_uint64_count, value_copy);
            small_modulo_uint_inplace(value_copy, value_uint64_count, smallmodulus, pool);
            set_uint_uint(value_copy, value_uint64_count, 1, result);
        }
        
        void increment_uint_smallmod(const uint64_t *operand, const SmallModulus &smallmodulus, uint64_t *result)
        {
#ifdef _DEBUG
            if (operand == nullptr)
            {
                throw invalid_argument("operand");
            }
            if (smallmodulus.value() == 0)
            {
                throw invalid_argument("modulus");
            }
            if (result == nullptr)
            {
                throw invalid_argument("result");
            }
            if (*operand >= smallmodulus.value())
            {
                throw out_of_range("operand");
            }
#endif
            //unsigned char carry = increment_uint(operand, 1, result);
            // Faster:
            *result = *operand + 1;
            *result -= smallmodulus.value() & static_cast<uint64_t>(-static_cast<int64_t>(*result >= smallmodulus.value()));
        }

        void decrement_uint_smallmod(const uint64_t *operand, const SmallModulus &smallmodulus, uint64_t *result)
        {
#ifdef _DEBUG
            if (operand == nullptr)
            {
                throw invalid_argument("operand");
            }
            if (smallmodulus.value() == 0)
            {
                throw invalid_argument("modulus");
            }
            if (result == nullptr)
            {
                throw invalid_argument("result");
            }
            if (*operand >= smallmodulus.value())
            {
                throw out_of_range("operand");
            }
#endif
            int64_t carry = (*operand == 0);
            *result = *operand - 1 + (smallmodulus.value() & static_cast<uint64_t>(-carry));
        }

        void negate_uint_smallmod(const uint64_t *operand, const SmallModulus &smallmodulus, uint64_t *result)
        {
#ifdef _DEBUG
            if (operand == nullptr)
            {
                throw invalid_argument("operand");
            }
            if (smallmodulus.value() == 0)
            {
                throw invalid_argument("modulus");
            }
            if (result == nullptr)
            {
                throw invalid_argument("result");
            }
            if (*operand >= smallmodulus.value())
            {
                throw out_of_range("operand");
            }
#endif
            *result = (smallmodulus.value() - *operand) & static_cast<uint64_t>(-static_cast<int64_t>(*operand != 0));
        }

        void div2_uint_smallmod(const uint64_t *operand, const SmallModulus &smallmodulus, uint64_t *result)
        {
#ifdef _DEBUG
            if (operand == nullptr)
            {
                throw invalid_argument("operand");
            }
            if (smallmodulus.value() == 0)
            {
                throw invalid_argument("modulus");
            }
            if (result == nullptr)
            {
                throw invalid_argument("result");
            }
            if (smallmodulus.value() == 0)
            {
                throw invalid_argument("modulus");
            }
            //if (is_greater_than_or_equal_uint_uint(operand, modulus, uint64_count))
            //{
            //    throw out_of_range("operand");
            //}
#endif
            if (*operand & 1)
            {
                int64_t carry = add_uint64(*operand, smallmodulus.value(), 0, result);
                *result >>= 1;
                if (carry)
                {
                    *result |= 1ULL << (bits_per_uint64 - 1);
                }
            }
            else
            {
                *result = *operand >> 1;
            }
        }

        void add_uint_uint_smallmod(const uint64_t *operand1, const uint64_t *operand2, const SmallModulus &smallmodulus, uint64_t *result)
        {
#ifdef _DEBUG
            if (operand1 == nullptr)
            {
                throw invalid_argument("operand1");
            }
            if (operand2 == nullptr)
            {
                throw invalid_argument("operand2");
            }
            if (smallmodulus.value() == 0)
            {
                throw invalid_argument("modulus");
            }
            if (result == nullptr)
            {
                throw invalid_argument("result");
            }
            if (*operand1 >= smallmodulus.value())
            {
                throw out_of_range("operand1");
            }
            if (*operand2 >= smallmodulus.value())
            {
                throw out_of_range("operand2");
            }
#endif
            bool carry = add_uint64(*operand1, *operand2, 0, result);
            *result -= smallmodulus.value() & static_cast<uint64_t>(-static_cast<int64_t>(carry || (*result >= smallmodulus.value())));
        }

        void sub_uint_uint_smallmod(const uint64_t *operand1, const uint64_t *operand2, const SmallModulus &smallmodulus, uint64_t *result)
        {
#ifdef _DEBUG
            if (operand1 == nullptr)
            {
                throw invalid_argument("operand1");
            }
            if (operand2 == nullptr)
            {
                throw invalid_argument("operand2");
            }
            if (smallmodulus.value() == 0)
            {
                throw invalid_argument("modulus");
            }
            if (result == nullptr)
            {
                throw invalid_argument("result");
            }
            if (*operand1 >= smallmodulus.value())
            {
                throw out_of_range("operand1");
            }
            if (*operand2 >= smallmodulus.value())
            {
                throw out_of_range("operand2");
            }
#endif
            int64_t borrow = sub_uint64(*operand1, *operand2, 0, result);
            *result += smallmodulus.value() & static_cast<uint64_t>(-borrow);
        }

        void barrett_reduce_128(const uint64_t *input, const SmallModulus &smallmodulus, uint64_t *result)
        {
            // Reduces input using base 2^64 Barrett reduction
            // input allocation size must be 128 bits

            uint64_t tmp1, tmp2[2], tmp3, carry;
            const uint64_t *const_ratio = smallmodulus.const_ratio().data();

            // Multiply input and const_ratio
            // Round 1
            multiply_uint64_hw64(input[0], const_ratio[0], &carry);

            multiply_uint64(input[0], const_ratio[1], tmp2);
            tmp3 = tmp2[1] + add_uint64(tmp2[0], carry, 0, &tmp1);

            // Round 2
            multiply_uint64(input[1], const_ratio[0], tmp2);
            carry = tmp2[1] + add_uint64(tmp1, tmp2[0], 0, &tmp1);

            // This is all we care about
            tmp1 = input[1] * const_ratio[1] + tmp3 + carry;

            // Barrett subtraction
            *result = input[0] - tmp1 * smallmodulus.value();

            // Claim: One more subtraction is enough
            *result -= smallmodulus.value() & static_cast<uint64_t>(-static_cast<int64_t>(*result >= smallmodulus.value()));
        }

        void multiply_uint64_smallmod(const uint64_t *operand1, const uint64_t *operand2, const SmallModulus &smallmodulus, uint64_t *result)
        {
#ifdef _DEBUG
            if (&operand1 == nullptr)
            {
                throw invalid_argument("operand1");
            }
            if (&operand2 == nullptr)
            {
                throw invalid_argument("operand2");
            }
            if (smallmodulus.value() == 0)
            {
                throw invalid_argument("modulus");
            }
            if (result == nullptr)
            {
                throw invalid_argument("result");
            }
#endif
            uint64_t z[2];
            multiply_uint64(*operand1, *operand2, z);
            barrett_reduce_128(z, smallmodulus, result);
        }

        bool is_primitive_root_smallmod(const uint64_t *root, uint64_t degree, const SmallModulus &smallmodulus)
        {
#ifdef _DEBUG
            if (root == nullptr)
            {
                throw invalid_argument("root");
            }
            if (smallmodulus.bit_count() < 2)
            {
                throw invalid_argument("modulus");
            }
            if (*root >= smallmodulus.value())
            {
                throw out_of_range("operand");
            }
            if (get_power_of_two(degree) < 1)
            {
                throw invalid_argument("degree must be a power of two and at least two");
            }
#endif
            if (*root == 0)
            {
                return false;
            }

            // We check if root is a degree-th root of unity in integers modulo modulus, where degree is a power of two.
            // It suffices to check that root^(degree/2) is -1 modulo modulus.
            uint64_t power;
            degree >>= 1;
            exponentiate_uint_smallmod(root, degree, smallmodulus, &power);
            increment_uint_smallmod(&power, smallmodulus, &power);

            return (power == 0);
        }
        
        bool try_primitive_root_smallmod(uint64_t degree, const SmallModulus &smallmodulus, MemoryPool &pool, uint64_t *destination)
        {
#ifdef _DEBUG
            if (destination == nullptr)
            {
                throw invalid_argument("destination");
            }
            if (smallmodulus.bit_count() < 2)
            {
                throw invalid_argument("modulus");
            }
            if (get_power_of_two(degree) < 1)
            {
                throw invalid_argument("degree must be a power of two and at least two");
            }
#endif
            // We need to divide modulus-1 by degree to get the size of the quotient group
            uint64_t size_entire_group = smallmodulus.value() - 1;

            // Compute size of quotient group
            uint64_t size_quotient_group = (smallmodulus.value() - 1) / degree;

            // size_entire_group must be divisible by degree, or otherwise the primitive root does not exist in integers modulo modulus
            if (size_entire_group - size_quotient_group * degree != 0)
            {
                return false;
            }

            // For randomness
            random_device rd;

            int attempt_counter = 0;
            int attempt_counter_max = 100;
            do
            {
                attempt_counter++;

                // Set destination to be a random number modulo modulus
                *reinterpret_cast<uint32_t*>(destination) = rd();
                *(reinterpret_cast<uint32_t*>(destination) + 1) = rd();

                small_modulo_uint_inplace(destination, 1, smallmodulus, pool);

                // Raise the random number to power the size of the quotient to get rid of irrelevant part
                exponentiate_uint_smallmod(destination, size_quotient_group, smallmodulus, destination);
            } while (!is_primitive_root_smallmod(destination, degree, smallmodulus) && attempt_counter < attempt_counter_max);

            return is_primitive_root_smallmod(destination, degree, smallmodulus);
        }
        
        bool try_minimal_primitive_root_smallmod(uint64_t degree, const SmallModulus &smallmodulus, MemoryPool &pool, uint64_t *destination)
        {
            if (!try_primitive_root_smallmod(degree, smallmodulus, pool, destination))
            {
                return false;
            }
            uint64_t generator_sq;
            multiply_uint64_smallmod(destination, destination, smallmodulus, &generator_sq);

            uint64_t current_generator = *destination;

            // destination is going to always contain the smallest generator found
            for (size_t i = 0; i < degree; i++)
            {
                // If our current generator is strictly smaller than destination, update
                if (current_generator < *destination)
                {
                    *destination = current_generator;
                }

                // Then move on to the next generator
                multiply_uint64_smallmod(&current_generator, &generator_sq, smallmodulus, &current_generator);
            }

            return true;
        }

    }
}
