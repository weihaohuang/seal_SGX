#include "util/uintcore.h"
#include "uintarithsmallmod.h"
#include "util/polycore.h"
#include "smallpolyarith.h"
#include "util/polyarith.h"
#include "polyarithsmallmod.h"
#include "util/polyfftmultmod.h"

using namespace std;

namespace seal
{
	namespace util
	{
		void dyadic_product_coeff_smallmod(const uint64_t *operand1, const uint64_t *operand2, int coeff_count, const SmallModulus &smallmodulus, uint64_t *result, MemoryPool &pool)
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
			if (result == nullptr)
			{
				throw invalid_argument("result");
			}
			if (coeff_count <= 0)
			{
				throw invalid_argument("coeff_count");
			}
			if (smallmodulus.uint64_count() <= 0 || smallmodulus.value() == 0)
			{
				throw invalid_argument("modulus");
			}
#endif
			// Use the same allocation for every instance of multiply_uint_uint_mod.
			for (int i = 0; i < coeff_count; i++)
			{
				multiply_uint64_smallmod(operand1, operand2, smallmodulus, result);
				operand1++;
				operand2++;
				result++;
			}
		}

		void smallmodulo_poly_inplace(uint64_t *value, int value_coeff_count, const PolyModulus &poly_modulus, const SmallModulus &smallmodulus, MemoryPool &pool)
		{
#ifdef _DEBUG
			if (value == nullptr)
			{
				throw invalid_argument("value");
			}
			if (value_coeff_count <= 0)
			{
				throw invalid_argument("value_coeff_count");
			}
			if (value == poly_modulus.get())
			{
				throw invalid_argument("value cannot point to same value as poly_modulus");
			}
#endif
			// Determine most significant coefficients of value and poly_modulus.
			int coeff_uint64_count = smallmodulus.uint64_count();
			//int value_coeffs = get_significant_coeff_count_poly(value, value_coeff_count, coeff_uint64_count);
			int value_coeffs = value_coeff_count;
			while (value_coeffs > 0 && value[value_coeffs - 1] == 0)
				value_coeffs--;
			int poly_modulus_coeff_count = poly_modulus.coeff_count();

			// If value has lesser degree than poly_modulus, then done.
			if (value_coeffs < poly_modulus_coeff_count)
			{
				return;
			}

			// Handle 1x^n + 1 polynomials more efficiently.
			const uint64_t *coeff_modulus = smallmodulus.pointer();
			if (poly_modulus.is_one_zero_one())
			{
				// Perform coefficient-wise division algorithm.
				while (value_coeffs >= poly_modulus_coeff_count)
				{
					// Determine leading value coefficient.
					uint64_t *leading_value_coeff = get_poly_coeff(value, value_coeffs - 1, coeff_uint64_count);

					// If leading value coefficient is not zero, then need to make zero by subtraction.
					if (!is_zero_uint(leading_value_coeff, coeff_uint64_count))
					{
						// Determine shift necesarry to bring significant coefficients in alignment.
						int poly_modulus_shift = value_coeffs - poly_modulus_coeff_count;

						// Subtract top coefficient from bottom-shifted coefficient.
						uint64_t *value_coeff = get_poly_coeff(value, poly_modulus_shift, coeff_uint64_count);
						sub_uint_uint_smallmod(value_coeff, leading_value_coeff, smallmodulus, value_coeff);

						// Zero-out leading coefficient.
						set_zero_uint(coeff_uint64_count, leading_value_coeff);
					}

					// Top value coefficient must now be zero, so adjust coefficient count.
					value_coeffs--;
				}
				return;
			}

			uint64_t monic_poly_modulus_scalar = 0, temp_quotient = 0, subtrahend = 0;

			// Determine scalar necessary to make poly_modulus monic.
			const uint64_t *polymodptr = poly_modulus.get();
			uint64_t leading_poly_modulus_coeff = polymodptr[poly_modulus_coeff_count - 1];
			if (!try_invert_uint_mod(&leading_poly_modulus_coeff, coeff_modulus, coeff_uint64_count, &monic_poly_modulus_scalar, pool))
			{
				throw invalid_argument("coeff_modulus is not coprime with leading poly_modulus coefficient");
			}

			// Perform coefficient-wise division algorithm.
			while (value_coeffs >= poly_modulus_coeff_count)
			{
				// Determine leading value coefficient.
				uint64_t leading_value_coeff = value[value_coeffs - 1];

				// If leading value coefficient is not zero, then need to make zero by subtraction.
				if (leading_value_coeff != 0)
				{
					// Determine shift necessary to bring significant coefficients in alignment.
					int poly_modulus_shift = value_coeffs - poly_modulus_coeff_count;

					// Determine quotient's coefficient, which is scalar that makes poly_modulus's leading coefficient one
					// multiplied by leading coefficient of poly_modulus (which when subtracted will zero out the topmost
					// poly_modulus coefficient).
					multiply_uint64_smallmod(&monic_poly_modulus_scalar, &leading_value_coeff, smallmodulus, &temp_quotient);

					// Subtract value and quotient*poly_modulus (shifted by poly_modulus_shift).
					for (int poly_modulus_coeff_index = 0; poly_modulus_coeff_index < poly_modulus_coeff_count; poly_modulus_coeff_index++)
					{
						// Multiply poly_modulus's coefficient by quotient.
						uint64_t poly_modulus_coeff = polymodptr[poly_modulus_coeff_index];
						if (poly_modulus_coeff != 0)
						{
							multiply_uint64_smallmod(&temp_quotient, &poly_modulus_coeff, smallmodulus, &subtrahend);

							// Subtract value with resulting product, appropriately shifted by poly_modulus shift.
							uint64_t *value_coeff = value + (poly_modulus_coeff_index + poly_modulus_shift);
							sub_uint_uint_smallmod(value_coeff, &subtrahend, smallmodulus, value_coeff);
						}
					}
				}

				// Top value coefficient must now be zero, so adjust coefficient count.
				value_coeffs--;
			}
		}

		void smallmodulo_poly(const uint64_t *value, int value_coeff_count, const PolyModulus &poly_modulus, const SmallModulus &smallmodulus, uint64_t *result, MemoryPool &pool)
		{
#ifdef _DEBUG
			if (value == nullptr)
			{
				throw invalid_argument("value");
			}
			if (value_coeff_count <= 0)
			{
				throw invalid_argument("value_coeff_count");
			}
			if (result == nullptr)
			{
				throw invalid_argument("result");
			}
#endif
			int coeff_uint64_count = smallmodulus.uint64_count();
			Pointer value_copy(allocate_poly(value_coeff_count, coeff_uint64_count, pool));
			set_poly_poly(value, value_coeff_count, coeff_uint64_count, value_copy.get());
			smallmodulo_poly_inplace(value_copy.get(), value_coeff_count, poly_modulus, smallmodulus, pool);
			set_poly_poly(value_copy.get(), poly_modulus.coeff_count(), coeff_uint64_count, result);
		}

		void nonfft_multiply_poly_poly_polymod_coeff_smallmod(const uint64_t *operand1, const uint64_t *operand2, const PolyModulus &poly_modulus, const SmallModulus &smallmodulus, uint64_t *result, MemoryPool &pool)
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
			if (result == nullptr)
			{
				throw invalid_argument("result");
			}
			if (get_significant_coeff_count_poly(operand1, poly_modulus.coeff_count(), poly_modulus.coeff_uint64_count()) >= poly_modulus.coeff_count())
			{
				throw out_of_range("operand1");
			}
			if (get_significant_coeff_count_poly(operand2, poly_modulus.coeff_count(), poly_modulus.coeff_uint64_count()) >= poly_modulus.coeff_count())
			{
				throw out_of_range("operand2");
			}
#endif
			// Calculate normal product.
			int coeff_count = poly_modulus.coeff_count();
			int coeff_uint64_count = poly_modulus.coeff_uint64_count();
			int intermediate_coeff_count = coeff_count * 2 - 1;
			Pointer intermediate(allocate_poly(intermediate_coeff_count, coeff_uint64_count, pool));
			multiply_poly_poly_coeff_smallmod(operand1, operand2, coeff_count, smallmodulus, intermediate.get());

			// Perform modulo operation.
			smallmodulo_poly_inplace(intermediate.get(), intermediate_coeff_count, poly_modulus, smallmodulus, pool);

			// Copy to result.
			set_poly_poly(intermediate.get(), coeff_count, coeff_uint64_count, result);
		}

		void nonfft_multiply_poly_poly_polymod_coeff_smallmod_inplace(const uint64_t *operand1, const uint64_t *operand2, const PolyModulus &poly_modulus, const SmallModulus &smallmodulus, uint64_t *result, MemoryPool &pool)
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
			if (result == nullptr)
			{
				throw invalid_argument("result");
			}
			if (get_significant_coeff_count_poly(operand1, poly_modulus.coeff_count(), poly_modulus.coeff_uint64_count()) >= poly_modulus.coeff_count())
			{
				throw out_of_range("operand1");
			}
			if (get_significant_coeff_count_poly(operand2, poly_modulus.coeff_count(), poly_modulus.coeff_uint64_count()) >= poly_modulus.coeff_count())
			{
				throw out_of_range("operand2");
			}
#endif
			// Calculate normal product.
			int coeff_count = poly_modulus.coeff_count();
			int result_coeff_count = coeff_count * 2 - 1;
			multiply_poly_poly_coeff_smallmod(operand1, operand2, coeff_count, smallmodulus, result);

			// Perform modulo operation.
			smallmodulo_poly_inplace(result, result_coeff_count, poly_modulus, smallmodulus, pool);
		}

		bool try_invert_poly_coeff_smallmod(const uint64_t *operand, const uint64_t *poly_modulus, int coeff_count, const SmallModulus &smallmodulus, uint64_t *result, MemoryPool &pool)
		{
#ifdef _DEBUG
			if (operand == nullptr)
			{
				throw invalid_argument("operand");
			}
			if (poly_modulus == nullptr)
			{
				throw invalid_argument("poly_modulus");
			}
			if (coeff_count <= 0)
			{
				throw invalid_argument("coeff_count");
			}
			if (result == nullptr)
			{
				throw invalid_argument("result");
			}
			if (get_significant_coeff_count_poly(operand, coeff_count, smallmodulus.uint64_count()) >= get_significant_coeff_count_poly(poly_modulus, coeff_count, smallmodulus.uint64_count()))
			{
				throw out_of_range("operand");
			}
#endif
			// Cannot invert 0 poly.
			int coeff_uint64_count = smallmodulus.uint64_count();
			if (is_zero_poly(operand, coeff_count, coeff_uint64_count))
			{
				return false;
			}

			// Construct a mutable copy of operand and modulus, with numerator being modulus
			// and operand being denominator. Notice that degree(numerator) >= degree(denominator).
			Pointer numerator_anchor(allocate_poly(coeff_count, coeff_uint64_count, pool));
			uint64_t *numerator = numerator_anchor.get();
			set_poly_poly(poly_modulus, coeff_count, coeff_uint64_count, numerator);
			Pointer denominator_anchor(allocate_poly(coeff_count, coeff_uint64_count, pool));
			uint64_t *denominator = denominator_anchor.get();
			set_poly_poly(operand, coeff_count, coeff_uint64_count, denominator);

			// Determine most significant coefficients of each.
			int numerator_coeffs = get_significant_coeff_count_poly(numerator, coeff_count, coeff_uint64_count);
			int denominator_coeffs = get_significant_coeff_count_poly(denominator, coeff_count, coeff_uint64_count);

			// Create poly to store quotient.
			Pointer quotient(allocate_poly(coeff_count, coeff_uint64_count, pool));

			// Create scalar to store value that makes denominator monic.
			Pointer monic_denominator_scalar(allocate_uint(coeff_uint64_count, pool));

			// Create temporary scalars used during calculation of quotient.
			// Both are purposely twice as wide to store intermediate product prior to modulo operation.
			int intermediate_uint64_count = coeff_uint64_count * 2;
			Pointer temp_quotient(allocate_uint(intermediate_uint64_count, pool));
			Pointer subtrahend(allocate_uint(intermediate_uint64_count, pool));

			// Create three polynomials to store inverse.
			// Initialize invert_prior to 0 and invert_curr to 1.
			Pointer invert_prior_anchor(allocate_poly(coeff_count, coeff_uint64_count, pool));
			uint64_t *invert_prior = invert_prior_anchor.get();
			set_zero_poly(coeff_count, coeff_uint64_count, invert_prior);
			Pointer invert_curr_anchor(allocate_poly(coeff_count, coeff_uint64_count, pool));
			uint64_t *invert_curr = invert_curr_anchor.get();
			set_zero_poly(coeff_count, coeff_uint64_count, invert_curr);
			uint64_t *invert_curr_first_coeff = get_poly_coeff(invert_curr, 0, coeff_uint64_count);
			set_uint(1, coeff_uint64_count, invert_curr_first_coeff);
			Pointer invert_next_anchor(allocate_poly(coeff_count, coeff_uint64_count, pool));
			uint64_t *invert_next = invert_next_anchor.get();

			//Pointer big_alloc(allocate_uint(7 * coeff_uint64_count, pool));

			// Perform extended Euclidean algorithm.
			const uint64_t *modulusptr = smallmodulus.pointer();
			while (true)
			{
				// NOTE: degree(numerator) >= degree(denominator).

				// Determine scalar necessary to make denominator monic.
				const uint64_t *leading_denominator_coeff = get_poly_coeff(denominator, denominator_coeffs - 1, coeff_uint64_count);
				if (!try_invert_uint_mod(leading_denominator_coeff, modulusptr, coeff_uint64_count, monic_denominator_scalar.get(), pool))
				{
					throw invalid_argument("coeff_modulus is not coprime with leading denominator coefficient");
				}

				// Clear quotient.
				set_zero_poly(coeff_count, coeff_uint64_count, quotient.get());

				// Perform coefficient-wise division algorithm.
				while (numerator_coeffs >= denominator_coeffs)
				{
					// Determine leading numerator coefficient.
					const uint64_t *leading_numerator_coeff = get_poly_coeff(numerator, numerator_coeffs - 1, coeff_uint64_count);

					// If leading numerator coefficient is not zero, then need to make zero by subtraction.
					if (!is_zero_uint(leading_numerator_coeff, coeff_uint64_count))
					{
						// Determine shift necessary to bring significant coefficients in alignment.
						int denominator_shift = numerator_coeffs - denominator_coeffs;

						// Determine quotient's coefficient, which is scalar that makes denominator's leading coefficient one
						// multiplied by leading coefficient of denominator (which when subtracted will zero out the topmost
						// denominator coefficient).
						uint64_t *quotient_coeff = get_poly_coeff(quotient.get(), denominator_shift, coeff_uint64_count);
						multiply_uint64_smallmod(monic_denominator_scalar.get(), leading_numerator_coeff, smallmodulus, temp_quotient.get());
						set_uint_uint(temp_quotient.get(), coeff_uint64_count, quotient_coeff);

						// Subtract numerator and quotient*denominator (shifted by denominator_shift).
						for (int denominator_coeff_index = 0; denominator_coeff_index < denominator_coeffs; denominator_coeff_index++)
						{
							// Multiply denominator's coefficient by quotient.
							const uint64_t *denominator_coeff = get_poly_coeff(denominator, denominator_coeff_index, coeff_uint64_count);
							multiply_uint64_smallmod(temp_quotient.get(), denominator_coeff, smallmodulus, subtrahend.get());

							// Subtract numerator with resulting product, appropriately shifted by denominator shift.
							uint64_t *numerator_coeff = get_poly_coeff(numerator, denominator_coeff_index + denominator_shift, coeff_uint64_count);
							sub_uint_uint_smallmod(numerator_coeff, subtrahend.get(), smallmodulus, numerator_coeff);
						}
					}

					// Top numerator coefficient must now be zero, so adjust coefficient count.
					numerator_coeffs--;
				}

				// Double check that numerator coefficients is correct because possible other coefficients are zero.
				numerator_coeffs = get_significant_coeff_count_poly(numerator, coeff_count, coeff_uint64_count);

				// We are done if numerator is zero.
				if (numerator_coeffs == 0)
				{
					break;
				}

				// Integrate quotient with invert coefficients.
				// Calculate: invert_next = invert_prior + -quotient * invert_curr
				multiply_truncate_poly_poly_coeff_smallmod(quotient.get(), invert_curr, coeff_count, smallmodulus, invert_next);
				sub_poly_poly_coeff_smallmod(invert_prior, invert_next, coeff_count, smallmodulus, invert_next);

				// Swap prior and curr, and then curr and next.
				swap(invert_prior, invert_curr);
				swap(invert_curr, invert_next);

				// Swap numerator and denominator.
				swap(numerator, denominator);
				swap(numerator_coeffs, denominator_coeffs);
			}

			// Polynomial is invertible only if denominator is just a scalar.
			if (denominator_coeffs != 1)
			{
				return false;
			}

			// Determine scalar necessary to make denominator monic.
			const uint64_t *leading_denominator_coeff = get_poly_coeff(denominator, 0, coeff_uint64_count);
			if (!try_invert_uint_mod(leading_denominator_coeff, modulusptr, coeff_uint64_count, monic_denominator_scalar.get(), pool))
			{
				throw invalid_argument("coeff_modulus is not coprime with leading denominator coefficient");
			}

			// Multiply inverse by scalar and done.
			multiply_poly_scalar_coeff_smallmod(invert_curr, coeff_count, monic_denominator_scalar.get(), smallmodulus, result);
			return true;
		}
	}
}
