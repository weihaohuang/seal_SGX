#include "util/common.h"
#include "util/uintcore.h"
#include <algorithm>

using namespace std;

namespace seal
{
    namespace util
    {
        void set_uint_uint(const uint64_t *value, int value_uint64_count, int result_uint64_count, uint64_t* result)
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
            if (result_uint64_count < 0)
            {
                throw invalid_argument("result_uint64_count");
            }
            if (result == nullptr && result_uint64_count > 0)
            {
                throw invalid_argument("result");
            }
#endif
            if (value == result || value_uint64_count == 0)
            {
                // Fast path to handle self assignment.
                result += value_uint64_count;
                for (int i = value_uint64_count; i < result_uint64_count; i++)
                {
                    *result++ = 0;
                }
                return;
            }
            for (int i = 0; i < result_uint64_count; i++)
            {
                *result++ = (i < value_uint64_count) ? *value++ : 0;
            }
        }

        int get_power_of_two_uint(const uint64_t *operand, int uint64_count)
        {
#ifdef _DEBUG
            if (operand == nullptr && uint64_count > 0)
            {
                throw invalid_argument("operand");
            }
            if (uint64_count < 0)
            {
                throw invalid_argument("uint64_count");
            }
#endif
            operand += uint64_count;
            int long_index = uint64_count, local_result = -1;
            for (; (long_index >= 1) && (local_result == -1); long_index--)
            {
                operand--;
                local_result = get_power_of_two(*operand);
            }

            // If local_result != -1, we've found a power-of-two highest order block,
            // in which case need to check that rest are zero.
            // If local_result == -1, operand is not power of two.
            if (local_result == -1)
            {
                return -1;
            }

            int zeros = 1;
            for (int j = long_index; j >= 1; j--)
            {
                zeros &= (*--operand == 0);
            }

            return zeros * (local_result + long_index * bits_per_uint64) + zeros - 1;
        }

        int get_power_of_two_minus_one_uint(const uint64_t *operand, int uint64_count)
        {
#ifdef _DEBUG
            if (operand == nullptr && uint64_count > 0)
            {
                throw invalid_argument("operand");
            }
            if (uint64_count < 0)
            {
                throw invalid_argument("uint64_count");
            }
#endif
            operand += uint64_count;
            int long_index = uint64_count, local_result = 0;
            for (; (long_index >= 1) && (local_result == 0); long_index--)
            {
                operand--;
                local_result = get_power_of_two_minus_one(*operand);
            }

            // If local_result != -1, we've found a power-of-two-minus-one highest order block,
            // in which case need to check that rest are ~0.
            // If local_result == -1, operand is not power of two minus one.
            if (local_result == -1)
            {
                return -1;
            }

            int ones = 1;
            for (int j = long_index; j >= 1; j--)
            {
                ones &= (~*--operand == static_cast<uint64_t>(0));
            }

            return ones * (local_result + long_index * bits_per_uint64) + ones - 1;
        }

        void filter_highbits_uint(uint64_t *operand, int uint64_count, int bit_count)
        {
#ifdef _DEBUG
            if (operand == nullptr && uint64_count > 0)
            {
                throw invalid_argument("operand");
            }
            if (bit_count < 0 || bit_count > uint64_count * bits_per_uint64)
            {
                throw invalid_argument("bit_count");
            }
            if (uint64_count < 0)
            {
                throw invalid_argument("uint64_count");
            }
#endif
            if (bit_count == uint64_count * bits_per_uint64)
            {
                return;
            }
            int uint64_index = bit_count / bits_per_uint64;
            int subbit_index = bit_count - uint64_index * bits_per_uint64;
            operand += uint64_index;
            *operand++ &= (1ULL << subbit_index) - 1;
            for (int long_index = uint64_index + 1; long_index < uint64_count; ++long_index)
            {
                *operand++ = 0;
            }
        }

        ConstPointer duplicate_uint_if_needed(const uint64_t *uint, int uint64_count, int new_uint64_count, bool force, MemoryPool &pool)
        {
#ifdef _DEBUG
            if (uint == nullptr && uint64_count > 0)
            {
                throw invalid_argument("uint");
            }
            if (uint64_count < 0)
            {
                throw invalid_argument("uint64_count");
            }
            if (new_uint64_count < 0)
            {
                throw invalid_argument("new_uint64_count");
            }
#endif
            if (!force && uint64_count >= new_uint64_count)
            {
                return ConstPointer::Aliasing(uint);
            }
            Pointer allocation(pool.get_for_uint64_count(new_uint64_count));
            set_uint_uint(uint, uint64_count, new_uint64_count, allocation.get());
            ConstPointer const_allocation;
            const_allocation.acquire(allocation);
            return const_allocation;
        }

        int compare_uint_uint(const uint64_t *operand1, const uint64_t *operand2, int uint64_count)
        {
#ifdef _DEBUG
            if (operand1 == nullptr && uint64_count > 0)
            {
                throw invalid_argument("operand1");
            }
            if (operand2 == nullptr && uint64_count > 0)
            {
                throw invalid_argument("operand2");
            }
            if (uint64_count < 0)
            {
                throw invalid_argument("uint64_count");
            }
#endif
            int result = 0;
            operand1 += uint64_count - 1;
            operand2 += uint64_count - 1;

            for (int i = uint64_count - 1; (result == 0) && (i >= 0); i--)
            {
                result = (*operand1 > *operand2) - (*operand1 < *operand2);
                operand1--;
                operand2--;
            }
            return result;
        }

        int compare_uint_uint(const uint64_t *operand1, int operand1_uint64_count, const uint64_t *operand2, int operand2_uint64_count)
        {
#ifdef _DEBUG
            if (operand1 == nullptr && operand1_uint64_count > 0)
            {
                throw invalid_argument("operand1");
            }
            if (operand1_uint64_count < 0)
            {
                throw invalid_argument("operand1_uint64_count");
            }
            if (operand2 == nullptr && operand2_uint64_count > 0)
            {
                throw invalid_argument("operand2");
            }
            if (operand2_uint64_count < 0)
            {
                throw invalid_argument("operand2_uint64_count");
            }
#endif
            int result = 0;
            operand1 += operand1_uint64_count - 1;
            operand2 += operand2_uint64_count - 1;

            int min_uint64_count = min(operand1_uint64_count, operand2_uint64_count);
            for (int i = operand1_uint64_count - 1; (result == 0) && (i >= min_uint64_count); i--)
            {
                result = (*operand1-- > 0);
            }
            for (int i = operand2_uint64_count - 1; (result == 0) && (i >= min_uint64_count); i--)
            {
                result = -(*operand2-- > 0);
            }
            for (int i = min_uint64_count - 1; (result == 0) && (i >= 0); i--)
            {
                result = (*operand1 > *operand2) - (*operand1 < *operand2);
                operand1--;
                operand2--;
            }
            return result;
        }
    }
}