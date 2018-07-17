#include <algorithm>
#include <string>
#include "util/common.h"
#include "util/uintarith.h"
#include <sstream>

using namespace std;

namespace seal
{
    namespace util
    {
        int get_power_of_two(uint64_t value)
        {
            if (value == 0 || (value & (value - 1)) != 0)
            {
                return -1;
            }
            int power = 0;
            while (value != 1)
            {
                power++;
                value >>= 1;
            }
            return power;
        }

		int get_close_power_of_two(uint64_t value)
		{
			if (value == 0)
			{
				return -1;
			}
			int power = 0;
			while (value != 1)
			{
				power++;
				value >>= 1;
			}
			return power;
		}


        int get_power_of_two_minus_one(uint64_t value)
        {
            if (value == 0xFFFFFFFFFFFFFFFF)
            {
                return bits_per_uint64;
            }
            return get_power_of_two(value + 1);
        }

        string uint_to_hex_string(const uint64_t *value, int uint64_count)
        {
#ifdef _DEBUG
            if (uint64_count < 0)
            {
                throw invalid_argument("uint64_count");
            }
            if (uint64_count > 0 && value == nullptr)
            {
                throw invalid_argument("value");
            }
#endif
            // Start with a string with a zero for each nibble in the array.
            int num_nibbles = uint64_count * nibbles_per_uint64;
            string output(num_nibbles, '0');

            // Iterate through each uint64 in array and set string with correct nibbles in hex.
            int nibble_index = num_nibbles;
            int leftmost_non_zero_pos = num_nibbles;
            for (int i = 0; i < uint64_count; ++i)
            {
                uint64_t part = *value++;

                // Iterate through each nibble in the current uint64.
                for (int j = 0; j < nibbles_per_uint64; ++j)
                {
                    int nibble = static_cast<int>(part & 0x0F);
                    int pos = --nibble_index;
                    if (nibble != 0)
                    {
                        // If nibble is not zero, then update string and save this pos to determine
                        // number of leading zeros.
                        output[pos] = nibble_to_upper_hex(nibble);
                        leftmost_non_zero_pos = pos;
                    }
                    part >>= 4;
                }
            }

            // Trim string to remove leading zeros.
            output.erase(0, leftmost_non_zero_pos);

            // Return 0 if nothing remains.
            if (output.empty())
            {
                return "0";
            }

            return output;
        }

        string poly_to_hex_string(const uint64_t *value, int coeff_count, int coeff_uint64_count)
        {
#ifdef _DEBUG
            if (coeff_count < 0)
            {
                throw invalid_argument("coeff_count");
            }
            if (coeff_uint64_count < 0)
            {
                throw invalid_argument("coeff_uint64_count");
            }
            if (coeff_uint64_count > 0 && coeff_count > 0 && value == nullptr)
            {
                throw invalid_argument("value");
            }
#endif
            ostringstream result;
            bool empty = true;
            value += (coeff_count - 1) * coeff_uint64_count;
            for (int i = coeff_count - 1; i >= 0; i--)
            {
                if (is_zero_uint(value, coeff_uint64_count))
                {
                    value -= coeff_uint64_count;
                    continue;
                }
                if (!empty)
                {
                    result << " + ";
                }
                result << uint_to_hex_string(value, coeff_uint64_count);
                if (i > 0)
                {
                    result << "x^" << i;
                }
                empty = false;
                value -= coeff_uint64_count;
            }
            if (empty)
            {
                result << "0";
            }
            return result.str();
        }

        string uint_to_dec_string(const uint64_t *value, int uint64_count, MemoryPool &pool)
        {
#ifdef _DEBUG
            if (uint64_count < 0)
            {
                throw invalid_argument("uint64_count");
            }
            if (uint64_count > 0 && value == nullptr)
            {
                throw invalid_argument("value");
            }
#endif
            if (uint64_count == 0)
            {
                return "0";
            }
            Pointer remainder(allocate_uint(uint64_count, pool));
            Pointer quotient(allocate_uint(uint64_count, pool));
            Pointer base(allocate_uint(uint64_count, pool));
            uint64_t *remainderptr = remainder.get();
            uint64_t *quotientptr = quotient.get();
            uint64_t *baseptr = base.get();
            set_uint(10, uint64_count, baseptr);
            set_uint_uint(value, uint64_count, remainderptr);
            string output;
            while (!is_zero_uint(remainderptr, uint64_count))
            {
                divide_uint_uint_inplace(remainderptr, baseptr, uint64_count, quotientptr, pool);
                char digit = static_cast<char>(remainderptr[0]) + '0';
                output += digit;
                swap(remainderptr, quotientptr);
            }
            reverse(output.begin(), output.end());

            // Return 0 if nothing remains.
            if (output.empty())
            {
                return "0";
            }

            return output;
        }

        string poly_to_dec_string(const uint64_t *value, int coeff_count, int coeff_uint64_count, MemoryPool &pool)
        {
#ifdef _DEBUG
            if (coeff_count < 0)
            {
                throw invalid_argument("coeff_count");
            }
            if (coeff_uint64_count < 0)
            {
                throw invalid_argument("coeff_uint64_count");
            }
            if (coeff_uint64_count > 0 && coeff_count > 0 && value == nullptr)
            {
                throw invalid_argument("value");
            }
#endif
            ostringstream result;
            bool empty = true;
            value += coeff_count - 1;
            for (int i = coeff_count - 1; i >= 0; i--)
            {
                if (is_zero_uint(value, coeff_uint64_count))
                {
                    value -= coeff_uint64_count;
                    continue;
                }
                if (!empty)
                {
                    result << " + ";
                }
                result << uint_to_dec_string(value, coeff_uint64_count, pool);
                if (i > 0)
                {
                    result << "x^" << i;
                }
                empty = false;
                value -= coeff_uint64_count;
            }
            if (empty)
            {
                result << "0";
            }
            return result.str();
        }

        void hex_string_to_uint(const char *hex_string, int char_count, int uint64_count, uint64_t *result)
        {
#ifdef _DEBUG
            if (hex_string == nullptr && char_count > 0)
            {
                throw invalid_argument("hex_string");
            }
            if (char_count < 0)
            {
                throw invalid_argument("char_count");
            }
            if (uint64_count < 0)
            {
                throw invalid_argument("uint64_count");
            }
            if (uint64_count > 0 && result == nullptr)
            {
                throw invalid_argument("result");
            }
            if (get_hex_string_bit_count(hex_string, char_count) > uint64_count * bits_per_uint64)
            {
                throw out_of_range("hex_string");
            }
#endif
            const char *hex_string_ptr = hex_string + char_count;
            for (int uint64_index = 0; uint64_index < uint64_count; ++uint64_index)
            {
                uint64_t value = 0;
                for (int bit_index = 0; bit_index < bits_per_uint64; bit_index += bits_per_nibble)
                {
                    if (hex_string_ptr == hex_string)
                    {
                        break;
                    }
                    char hex = *--hex_string_ptr;
                    int nibble = hex_to_nibble(hex);
                    if (nibble == -1)
                    {
                        throw invalid_argument("hex_value");
                    }
                    value |= static_cast<uint64_t>(nibble) << bit_index;
                }
                result[uint64_index] = value;
            }
        }

        int get_hex_string_bit_count(const char *hex_string, int char_count)
        {
#ifdef _DEBUG
            if (hex_string == nullptr && char_count > 0)
            {
                throw invalid_argument("hex_string");
            }
            if (char_count < 0)
            {
                throw invalid_argument("char_count");
            }
#endif
            for (int i = 0; i < char_count; ++i)
            {
                char hex = *hex_string++;
                int nibble = hex_to_nibble(hex);
                if (nibble != 0)
                {
                    int nibble_bits = get_significant_bit_count(nibble);
                    int remaining_nibbles = (char_count - i - 1) * bits_per_nibble;
                    return nibble_bits + remaining_nibbles;
                }
            }
            return 0;
        }

		/////////////////////////// New functions /////////////////////////////////

		int bit_reversal_general(int x, int k)
		{
			// Assume k <= 32. Use the 32-bit bit-reversal function 
			uint32_t res = reverse_bits(static_cast<uint32_t>(x));
			return static_cast<int> (res >> (32 - k));
		}

		void babystep_giantstep(int n, std::vector<uint64_t> &baby_step, std::vector<uint64_t> &giant_step)
		{
			int exponent = get_power_of_two(n);
			if (exponent < 0) {
				throw invalid_argument("n must be a power of 2");
			}
			int k, l;
			k = 1 << (exponent / 2);
			l = n / k;
			// k stores the baby steps. 
			baby_step.clear();
			giant_step.clear();

			int m = n << 1;
			int g = 5; // the generator; 
			int kprime = k >> 1;
			int value = 1;
			for (int i = 0; i < kprime; i++) {
				baby_step.push_back(value);
				baby_step.push_back(m - value);
				value *= g;
				value %= m;
			}

			// now value1 should equal to g**kprime 
			int value2 = value;
			for (int j = 0; j < l; j++) {
				giant_step.push_back(value2);
				value2 *= value;
				value2 %= m;
			}
			return;
		}

		std::pair<int, int>  decompose_bs_gs(uint64_t m, uint64_t x, std::vector<uint64_t>& baby_step, std::vector<uint64_t>& giant_step)
		{
			for (int i = 0; i < giant_step.size(); i++) {
				uint64_t gs = giant_step[i];
				for (int j = 0; j < baby_step.size(); j++) {
					uint64_t bs = baby_step[j];
					if ((gs*bs) % m == x) return { i,j };
				}
			}
			throw;
		}

		std::vector<uint64_t> conjugate_classes(uint64_t m, uint64_t p)
		{
			std::vector<uint64_t> classes;
			int j;
			for (int i = 0; i < m; i++) {
				if (gcd(i, m) > 1) {
					classes.push_back(0);
				}
				else {
					classes.push_back(i);
				}
			}
			for (int i = 0; i < m; i++) {
				if (classes[i] == 0) {
					continue;
				}
				if (classes[i] < i) {  // i is not a pivot, updated its pivot
					classes[i] = classes[classes[i]];
					continue;
				}
				//       # If i is a pivot, update other pivots to point to it
				j = i*p % m;
				while (classes[j] != i) {
					classes[classes[j]] = i; // # Merge the equivalence classes of j and i
											 // # Note: if classes[j] != j then classes[j] will be updated later, #     
											 //when we get to i = j and use the code for "i not pivot".
					j = j*p % m;
				}
			}
			return classes;
		}

		std::vector<uint64_t> multiplicative_orders(std::vector<uint64_t> classes, uint64_t m)
		{
			std::vector<uint64_t> orders;
			orders.push_back(0);
			orders.push_back(1);
			int j, order;
			for (int i = 2; i < m; i++) {
				if (classes[i] <= 1) {
					orders.push_back(classes[i]);
					continue;
				}
				if (classes[i] < i) {
					orders.push_back(orders[classes[i]]);
					continue;
				}
				j = i * i % m;
				order = 2;
				while (classes[j] != 1) {
					j = j*i % m;
					order++;
				}
				orders.push_back(order);

			}
			return orders;
		}

		int gcd(int x, int y)
		{
			if (x < y)
				return gcd(y, x);
			if (y == 0) {
				return x;
			}
			int f = x % y;
			if (f == 0)
				return y;
			else
				return gcd(y, f);
		}

		std::vector<int64_t> xgcd(int64_t a, int64_t b)
		{
			/* Extended GCD:
			Returns(gcd, x, y) where gcd is the greatest common divisor of a and b
			with the sign of b if b is nonzero, and with the sign of a if b is 0.
			The numbers x, y are such that gcd = ax + by."""
			*/
			int64_t prevx = 1;
			int64_t x = 0;
			int64_t prevy = 0;
			int64_t y = 1;
			int64_t q;
			int64_t temp1, temp2, temp3;
			while (b != 0) {
				q = a / b;
				temp1 = a % b;
				a = b;
				b = temp1;

				temp2 = x;
				x = prevx - q*x;
				prevx = temp2;

				temp3 = y;
				y = prevy - q*y;
				prevy = temp3;
			}
			std::vector<int64_t> result = { a, prevx, prevy };
			return result;
		}

		uint64_t mod_inverse(int64_t a, int64_t n)
		{
			std::vector<int64_t> result = xgcd(a, n);
			if (result[0] != 1) {
				throw "Inverse does not exist!";
			}
			if (result[1] < 0) {
				result[1] += n;
			}
			return result[1];
		}


		uint64_t power(uint64_t base, uint64_t exponent)
		{
			uint64_t result = 1;
			for (size_t i = 0; i < exponent; i++) {
				result *= base;
			}
			return result;
		}

		int HW(uint64_t x)
		{
			int res = 0;
			while (x)
			{
				res++;
				x &= x - 1;
			}
			return res;
		}

		uint64_t optimal_slpit(uint64_t x)
		{
			int hwx = HW(x);
			int target = hwx / 2;
			int now = 0;
			int xbit;
			uint64_t result = 0;
			for (int i = 0; i < 64; i++)
			{
				xbit = x & 1;
				x = x >> 1;
				now += xbit;
				result += (xbit << i);
				if (now >= target) { break; }
			}
			return result;
		}
    }
}
