#pragma once

#include <cstdint>
#include "smallmodulus.h"
#include "util/mempool.h"
#include "smallntt.h"
#include "util/polymodulus.h"

namespace seal
{
	namespace util
	{
		void smallntt_multiply_poly_poly(const std::uint64_t *operand1, const std::uint64_t *operand2, const SmallNTTTables &tables, std::uint64_t *result, MemoryPool &pool);

		void smallntt_multiply_poly_nttpoly(const std::uint64_t *operand1, const std::uint64_t *operand2, const SmallNTTTables &tables, std::uint64_t *result, MemoryPool &pool);

		void smallntt_double_multiply_poly_nttpoly(const std::uint64_t *operand1, const std::uint64_t *operand2, const std::uint64_t *operand3, const SmallNTTTables &tables, std::uint64_t *result1, std::uint64_t *result2, MemoryPool &pool);

		void nussbaumer_multiply_poly_poly_coeff_smallmod(const std::uint64_t *operand1, const std::uint64_t *operand2, int coeff_count_power, const SmallModulus &modulus, std::uint64_t *result, MemoryPool &pool);

		void smallntt_dot_product_bigpolyarray_nttbigpolyarray(const std::uint64_t *array1, const std::uint64_t *array2, int count, int array_poly_uint64_count, const SmallNTTTables &tables, std::uint64_t *result, MemoryPool &pool);

		void smallntt_double_dot_product_bigpolyarray_nttbigpolyarrays(const std::uint64_t *array1, const std::uint64_t *array2, const std::uint64_t *array3, int count, int array_poly_uint64_count, const SmallNTTTables &tables, std::uint64_t *result1, std::uint64_t *result2, MemoryPool &pool);

		void nussbaumer_dot_product_bigpolyarray_coeff_smallmod(const std::uint64_t *array1, const std::uint64_t *array2, int count, const PolyModulus &poly_modulus, const SmallModulus &modulus, std::uint64_t *result, MemoryPool &pool);
	}
}