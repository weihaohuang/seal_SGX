#include "util/uintcore.h"
#include "util/uintarith.h"
#include "uintarithsmallmod.h"
#include "util/polycore.h"
#include "smallpolyarith.h"

using namespace std;

namespace seal
{
    namespace util
    {

        void smallmod_poly_coeffs(uint64_t *poly, int coeff_count, const SmallModulus &smallmodulus, MemoryPool &pool)
        {
#ifdef _DEBUG
            if (poly == nullptr && coeff_count > 0)
            {
                throw invalid_argument("poly");
            }
            if (coeff_count < 0)
            {
                throw invalid_argument("coeff_count");
            }
#endif
            for (int i = 0; i < coeff_count; i++, poly++)
            {
                small_modulo_uint_inplace(poly, 1, smallmodulus, pool);
            }
        }

        void negate_poly_coeff_smallmod(const uint64_t *poly, int coeff_count, const SmallModulus &smallmodulus, int coeff_uint64_count, uint64_t *result)
        {
#ifdef _DEBUG
            if (poly == nullptr && coeff_count > 0)
            {
                throw invalid_argument("poly");
            }
            if (coeff_count < 0)
            {
                throw invalid_argument("coeff_count");
            }
            if (smallmodulus.value() == 0)
            {
                throw invalid_argument("smallmodulus");
            }
            if (coeff_uint64_count <= 0)
            {
                throw invalid_argument("coeff_uint64_count");
            }
            if (result == nullptr && coeff_count > 0)
            {
                throw invalid_argument("result");
            }
#endif
            for (int i = 0; i < coeff_count; i++)
            {
                negate_uint_smallmod(poly, smallmodulus, result);
                poly += coeff_uint64_count;
                result += coeff_uint64_count;
            }
        }

        void add_poly_poly_coeff_smallmod(const uint64_t *operand1, const uint64_t *operand2, int coeff_count, const SmallModulus &smallmodulus, uint64_t *result)
        {
#ifdef _DEBUG
            if (operand1 == nullptr && coeff_count > 0)
            {
                throw invalid_argument("operand1");
            }
            if (operand2 == nullptr && coeff_count > 0)
            {
                throw invalid_argument("operand2");
            }
            if (coeff_count < 0)
            {
                throw invalid_argument("coeff_count");
            }
            if (smallmodulus.value() == 0)
            {
                throw invalid_argument("smallmodulus");
            }
            if (result == nullptr && coeff_count > 0)
            {
                throw invalid_argument("result");
            }
#endif
            for (int i = 0; i < coeff_count; i++)
            {
                add_uint_uint_smallmod(operand1++, operand2++, smallmodulus, result++);
            }
        }

        void sub_poly_poly_coeff_smallmod(const uint64_t *operand1, const uint64_t *operand2, int coeff_count, const SmallModulus &smallmodulus, uint64_t *result)
        {
#ifdef _DEBUG
            if (operand1 == nullptr && coeff_count > 0)
            {
                throw invalid_argument("operand1");
            }
            if (operand2 == nullptr && coeff_count > 0)
            {
                throw invalid_argument("operand2");
            }
            if (coeff_count < 0)
            {
                throw invalid_argument("coeff_count");
            }
            if (smallmodulus.value() == 0)
            {
                throw invalid_argument("smallmodulus");
            }
            if (result == nullptr && coeff_count > 0)
            {
                throw invalid_argument("result");
            }
#endif
            for (int i = 0; i < coeff_count; i++)
            {
                sub_uint_uint_smallmod(operand1, operand2, smallmodulus, result);
                operand1++;
                operand2++;
                result++;
            }
        }

        void multiply_poly_scalar_coeff_smallmod(const uint64_t *poly, int coeff_count, const uint64_t *scalar, const SmallModulus &smallmodulus, uint64_t *result)
        {
#ifdef _DEBUG
            if (poly == nullptr && coeff_count > 0)
            {
                throw invalid_argument("poly");
            }
            if (coeff_count < 0)
            {
                throw invalid_argument("coeff_count");
            }
            if (scalar == nullptr)
            {
                throw invalid_argument("scalar");
            }
            if (result == nullptr && coeff_count > 0)
            {
                throw invalid_argument("result");
            }
#endif
            for (int i = 0; i < coeff_count; i++)
            {
                multiply_uint64_smallmod(poly, scalar, smallmodulus, result);
                poly++;
                result++;
            }
        }

        void multiply_poly_poly_coeff_smallmod(const uint64_t *operand1, int operand1_coeff_count, const uint64_t *operand2, int operand2_coeff_count,
            const SmallModulus &smallmodulus, int result_coeff_count, uint64_t *result)
        {
#ifdef _DEBUG
            if (operand1 == nullptr && operand1_coeff_count > 0)
            {
                throw invalid_argument("operand1");
            }
            if (operand1_coeff_count < 0)
            {
                throw invalid_argument("operand1_coeff_count");
            }
            if (operand2 == nullptr && operand2_coeff_count > 0)
            {
                throw invalid_argument("operand2");
            }
            if (operand2_coeff_count < 0)
            {
                throw invalid_argument("operand2_coeff_count");
            }
            if (result_coeff_count < 0)
            {
                throw invalid_argument("result_coeff_count");
            }
            if (result == nullptr && result_coeff_count > 0)
            {
                throw invalid_argument("result");
            }
            if (result != nullptr && (operand1 == result || operand2 == result))
            {
                throw invalid_argument("result cannot point to the same value as operand1, operand2, or modulus");
            }
#endif
            // Clear product.
            int result_coeff_uint64_count = smallmodulus.uint64_count();
            set_zero_poly(result_coeff_count, result_coeff_uint64_count, result);

            const uint64_t *modulusptr = smallmodulus.pointer();
            operand1_coeff_count = get_significant_coeff_count_poly(operand1, operand1_coeff_count, 1);
            operand2_coeff_count = get_significant_coeff_count_poly(operand2, operand2_coeff_count, 1);
            for (int operand1_index = 0; operand1_index < operand1_coeff_count; operand1_index++)
            {
                const uint64_t *operand1_coeff = operand1 + operand1_index;
                if (*operand1_coeff == 0)
                {
                    // If coefficient is 0, then move on to next coefficient.
                    continue;
                }
                // Lastly, do more expensive add if other cases don't handle it.
                for (int operand2_index = 0; operand2_index < operand2_coeff_count; operand2_index++)
                {
                    int product_coeff_index = operand1_index + operand2_index;
                    if (product_coeff_index >= result_coeff_count)
                    {
                        break;
                    }

                    const uint64_t *operand2_coeff = operand2 + operand2_index;
                    if (*operand2_coeff == 0)
                    {
                        // If coefficient is 0, then move on to next coefficient.
                        continue;
                    }

                    uint64_t temp;
                    multiply_uint64_smallmod(operand1_coeff, operand2_coeff, smallmodulus, &temp);

                    uint64_t *result_coeff = get_poly_coeff(result, product_coeff_index, result_coeff_uint64_count);
                    add_uint_uint_smallmod(result_coeff, &temp, smallmodulus, result_coeff);
                }
            }
        }

        void multiply_poly_poly_coeff_smallmod(const uint64_t *operand1, const uint64_t *operand2, int coeff_count, const SmallModulus &smallmodulus, uint64_t *result)
        {
            int result_coeff_count = coeff_count + coeff_count - 1;

            // Clear product.
            set_zero_poly(result_coeff_count, 1, result);

            for (int operand1_index = 0; operand1_index < coeff_count; operand1_index++)
            {
                //const uint64_t *operand1_coeff = get_poly_coeff(operand1, operand1_index, smallmodulus.uint64_count());
                uint64_t operand1_coeff = operand1[operand1_index];
                if (operand1_coeff == 0)
                {
                    // If coefficient is 0, then move on to next coefficient.
                    continue;
                }
                // Lastly, do more expensive add if other cases don't handle it.
                for (int operand2_index = 0; operand2_index < coeff_count; operand2_index++)
                {
                    int product_coeff_index = operand1_index + operand2_index;

                    //const uint64_t *operand2_coeff = get_poly_coeff(operand2, operand2_index, smallmodulus.uint64_count());
                    uint64_t operand2_coeff = operand2[operand2_index];
                    if (operand2_coeff == 0)
                    {
                        // If coefficient is 0, then move on to next coefficient.
                        continue;
                    }
                    uint64_t temp;
                    multiply_uint64_smallmod(&operand1_coeff, &operand2_coeff, smallmodulus, &temp);
                    uint64_t *result_coeff = result + product_coeff_index;
                    add_uint_uint_smallmod(result_coeff, &temp, smallmodulus, result_coeff);
                }
            }
        }

        void divide_poly_poly_coeff_smallmod_inplace(uint64_t *numerator, const uint64_t *denominator, int coeff_count, const SmallModulus &smallmodulus, uint64_t *quotient, MemoryPool &pool)
        {
#ifdef _DEBUG
            if (numerator == nullptr)
            {
                throw invalid_argument("numerator");
            }
            if (denominator == nullptr)
            {
                throw invalid_argument("denominator");
            }
            if (is_zero_poly(denominator, coeff_count, smallmodulus.uint64_count()))
            {
                throw invalid_argument("denominator");
            }
            if (coeff_count < 0)
            {
                throw invalid_argument("coeff_count");
            }
            if (quotient == nullptr)
            {
                throw invalid_argument("quotient");
            }
            if (numerator == quotient || denominator == quotient)
            {
                throw invalid_argument("quotient cannot point to same value as numerator or denominator");
            }
            if (numerator == denominator)
            {
                throw invalid_argument("numerator cannot point to same value as denominator");
            }
#endif
            // Clear quotient.
            int coeff_uint64_count = smallmodulus.uint64_count();
            set_zero_poly(coeff_count, coeff_uint64_count, quotient);

            // Determine most significant coefficients of numerator and denominator.
            int numerator_coeffs = get_significant_coeff_count_poly(numerator, coeff_count, coeff_uint64_count);
            int denominator_coeffs = get_significant_coeff_count_poly(denominator, coeff_count, coeff_uint64_count);

            // If numerator has lesser degree than denominator, then done.
            if (numerator_coeffs < denominator_coeffs)
            {
                return;
            }

            int intermediate_uint64_count = coeff_uint64_count * 2;
            Pointer big_alloc(allocate_uint(coeff_uint64_count + intermediate_uint64_count + intermediate_uint64_count + 7 * coeff_uint64_count, pool));

            // Create scalar to store value that makes denominator monic.
            uint64_t *monic_denominator_scalar = big_alloc.get();

            // Create temporary scalars used during calculation of quotient.
            // Both are purposely twice as wide to store intermediate product prior to modulo operation.
            uint64_t *temp_quotient = monic_denominator_scalar + coeff_uint64_count;
            uint64_t *subtrahend = temp_quotient + intermediate_uint64_count;

            // We still have 7 x coeff_uint64_count of memory left in the big allocation
            uint64_t *alloc_ptr = subtrahend + intermediate_uint64_count;

            // Determine scalar necessary to make denominator monic.
            const uint64_t *modulusptr = smallmodulus.pointer();
            const uint64_t *leading_denominator_coeff = get_poly_coeff(denominator, denominator_coeffs - 1, coeff_uint64_count);
            if (!try_invert_uint_mod(leading_denominator_coeff, modulusptr, coeff_uint64_count, monic_denominator_scalar, pool))
            {
                throw invalid_argument("smallmodulus is not coprime with leading denominator coefficient");
            }

            // Perform coefficient-wise division algorithm.
            while (numerator_coeffs >= denominator_coeffs)
            {
                // Determine leading numerator coefficient.
                const uint64_t *leading_numerator_coeff = get_poly_coeff(numerator, numerator_coeffs - 1, coeff_uint64_count);

                // If leading numerator coefficient is not zero, then need to make zero by subtraction.
                if (!is_zero_uint(leading_numerator_coeff, coeff_uint64_count))
                {
                    // Determine shift necesarry to bring significant coefficients in alignment.
                    int denominator_shift = numerator_coeffs - denominator_coeffs;

                    // Determine quotient's coefficient, which is scalar that makes denominator's leading coefficient one
                    // multiplied by leading coefficient of denominator (which when subtracted will zero out the topmost
                    // denominator coefficient).
                    uint64_t *quotient_coeff = get_poly_coeff(quotient, denominator_shift, coeff_uint64_count);
                    multiply_uint64_smallmod(monic_denominator_scalar, leading_numerator_coeff, smallmodulus, temp_quotient);
                    set_uint_uint(temp_quotient, coeff_uint64_count, quotient_coeff);

                    // Subtract numerator and quotient*denominator (shifted by denominator_shift).
                    for (int denominator_coeff_index = 0; denominator_coeff_index < denominator_coeffs; denominator_coeff_index++)
                    {
                        // Multiply denominator's coefficient by quotient.
                        const uint64_t *denominator_coeff = get_poly_coeff(denominator, denominator_coeff_index, coeff_uint64_count);
                        multiply_uint64_smallmod(temp_quotient, denominator_coeff, smallmodulus, subtrahend);

                        // Subtract numerator with resulting product, appropriately shifted by denominator shift.
                        uint64_t *numerator_coeff = get_poly_coeff(numerator, denominator_coeff_index + denominator_shift, coeff_uint64_count);
                        sub_uint_uint_smallmod(numerator_coeff, subtrahend, smallmodulus, numerator_coeff);
                    }
                }

                // Top numerator coefficient must now be zero, so adjust coefficient count.
                numerator_coeffs--;
            }
        }

        void add_bigpolyarray_coeff_smallmod(const uint64_t *array1, const uint64_t *array2, int count, int coeff_count, const SmallModulus &smallmodulus, uint64_t *result)
        {
            // Check validity of inputs
#ifdef _DEBUG
            if (array1 == nullptr)
            {
                throw invalid_argument("array1");
            }
            if (array2 == nullptr)
            {
                throw invalid_argument("array2");
            }
            if (result == nullptr)
            {
                throw invalid_argument("result");
            }
            if (count < 1)
            {
                throw invalid_argument("count");
            }
            if (coeff_count < 1)
            {
                throw invalid_argument("coeff_count");
            }
#endif
            // initialise pointers for addition
            const uint64_t *current_array1 = array1;
            const uint64_t *current_array2 = array2;
            uint64_t *current_result = result;

            for (int i = 0; i < count; i++)
            {
                add_poly_poly_coeff_smallmod(current_array1, current_array2, coeff_count, smallmodulus, current_result);
                current_array1 += coeff_count;
                current_array2 += coeff_count;
                current_result += coeff_count;
            }
        }

		void apply_galois_poly_smallmod(const uint64_t *input, int coeff_count, int galois_elt, const SmallModulus &smallmodulus, uint64_t *result)
		{
			const uint64_t *value_coeff = input;
			int index_raw, index, multiples;
			for (int i = 0; i < coeff_count; i++)
			{
				index_raw = i * galois_elt;
				index = index_raw % coeff_count;
				multiples = index_raw / coeff_count;
				if (multiples % 2 == 0)
				{
					set_uint_uint(value_coeff, 1, result + index);
				}
				else
				{
					negate_uint_smallmod(value_coeff, smallmodulus, result + index);
				}
				value_coeff++;
			}
		}

		void negacyclic_shift_poly_smallmod(const std::uint64_t *value, int coeff_count, int shift, uint64_t * result, const SmallModulus &smallmodulus)
		{
			if (shift < 0)
			{
				throw invalid_argument("shift must be non-negative");
			}
			const uint64_t *value_coeff = value;
			int index_raw, index, multiples;
			for (int i = 0; i < coeff_count; i++)
			{
				index_raw = i + shift;
				index = index_raw % coeff_count;
				multiples = index_raw / coeff_count;
				if (multiples % 2 == 0) {
					set_uint_uint(value_coeff, 1, result + index);
				}
				else {
					negate_uint_smallmod(value_coeff, smallmodulus, result + index);
				}
				value_coeff += 1;
			}
		}

		void permute_ntt_poly_smallmod(const std::uint64_t * value, int coeff_count, int galois_elt, uint64_t * result)
		{
			const uint64_t *value_coeff = value;
			int index_raw, index, reversed;
			int n = coeff_count;
			int m = 2 * n;
			//int index_raw, index, multiples;
			int logn = get_power_of_two(coeff_count);
			for (int i = 0; i < coeff_count; i++)
			{
				reversed = bit_reversal_general(i, logn);
				index_raw = (2 * reversed + 1)*galois_elt;
				index_raw %= m;
				index = bit_reversal_general((index_raw - 1) >> 1, logn);
				set_uint_uint(value_coeff + index, 1, result + i);
			}
		}

		void permute_ntt_poly_smallmod_with_index(const std::uint64_t * value, int coeff_count, vector<int> &index, uint64_t * result)
		{
			const uint64_t *value_coeff = value;
			for (int i = 0; i < coeff_count; i++)
			{
				set_uint_uint(value_coeff + index[i], 1, result + i);
			}
		}

    }
}
