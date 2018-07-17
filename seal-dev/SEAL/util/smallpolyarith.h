#pragma once

#include <cstdint>
#include "util/polycore.h"
#include "util/mempool.h"
#include "util/uintarithmod.h"
#include "uintarithsmallmod.h"

namespace seal
{
	namespace util
	{
		void smallmod_poly_coeffs(std::uint64_t *poly, int coeff_count, const SmallModulus &modulus, MemoryPool &pool);

		void negate_poly_coeff_smallmod(const std::uint64_t *poly, int coeff_count, const SmallModulus &modulus, int coeff_uint64_count, std::uint64_t *result);

		void add_poly_poly_coeff_smallmod(const std::uint64_t *operand1, const std::uint64_t *operand2, int coeff_count, const SmallModulus &coeff_modulus, std::uint64_t *result);

		void sub_poly_poly_coeff_smallmod(const std::uint64_t *operand1, const std::uint64_t *operand2, int coeff_count, const SmallModulus &coeff_modulus, std::uint64_t *result);

		void multiply_poly_scalar_coeff_smallmod(const std::uint64_t *poly, int coeff_count, const std::uint64_t *scalar, const SmallModulus &modulus, std::uint64_t *result);

		void multiply_poly_poly_coeff_smallmod(const std::uint64_t *operand1, int operand1_coeff_count, const std::uint64_t *operand2, int operand2_coeff_count,
			const SmallModulus &modulus, int result_coeff_count, std::uint64_t *result);

		void multiply_poly_poly_coeff_smallmod(const std::uint64_t *operand1, const std::uint64_t *operand2, int coeff_count, const SmallModulus &modulus, std::uint64_t *result);

		inline void multiply_truncate_poly_poly_coeff_smallmod(const std::uint64_t *operand1, const std::uint64_t *operand2, int coeff_count, const SmallModulus &modulus, std::uint64_t *result)
		{
			multiply_poly_poly_coeff_smallmod(operand1, coeff_count, operand2, coeff_count, modulus, coeff_count, result);
		}

		void divide_poly_poly_coeff_smallmod_inplace(std::uint64_t *numerator, const std::uint64_t *denominator, int coeff_count, const SmallModulus &modulus, std::uint64_t *quotient, MemoryPool &pool);

		inline void divide_poly_poly_coeff_smallmod(const std::uint64_t *numerator, const std::uint64_t *denominator, int coeff_count, const SmallModulus &modulus, std::uint64_t *quotient, std::uint64_t *remainder, MemoryPool &pool)
		{
			int coeff_uint64_count = modulus.uint64_count();
			set_poly_poly(numerator, coeff_count, coeff_uint64_count, remainder);
			divide_poly_poly_coeff_smallmod_inplace(remainder, denominator, coeff_count, modulus, quotient, pool);
		}

		void add_bigpolyarray_coeff_smallmod(const std::uint64_t *array1, const std::uint64_t *array2, int count, int coeff_count, const SmallModulus &modulus, std::uint64_t *result);

		////////////////////////////////////////////////////

		void apply_galois_poly_smallmod(const uint64_t *input, int coeff_count, int galois_elt, const SmallModulus &smallmodulus, uint64_t *result);

		void permute_ntt_poly_smallmod(const std::uint64_t *value, int coeff_count, int galois_elt, uint64_t *result);

		void permute_ntt_poly_smallmod_with_index(const std::uint64_t *value, int coeff_count, std::vector<int> &index, uint64_t *result);

		void negacyclic_shift_poly_smallmod(const std::uint64_t *value, int coeff_count, int shift, uint64_t *result, const SmallModulus &smallmodulus);

	}
}