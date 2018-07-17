#pragma once

#include <cstdint>
#include "smallmodulus.h"
#include "util/mempool.h"
#include "util/polymodulus.h"

namespace seal
{
	namespace util
	{
		void dyadic_product_coeff_smallmod(const std::uint64_t *operand1, const std::uint64_t *operand2, int coeff_count, const SmallModulus &modulus, std::uint64_t *result, MemoryPool &pool);

		void smallmodulo_poly_inplace(std::uint64_t *value, int value_coeff_count, const PolyModulus &poly_modulus, const SmallModulus &modulus, MemoryPool &pool);

		void smallmodulo_poly(const std::uint64_t *value, int value_coeff_count, const PolyModulus &poly_modulus, const SmallModulus &modulus, std::uint64_t *result, MemoryPool &pool);

		void nonfft_multiply_poly_poly_polymod_coeff_smallmod(const std::uint64_t *operand1, const std::uint64_t *operand2, const PolyModulus &poly_modulus, const SmallModulus &modulus, std::uint64_t *result, MemoryPool &pool);

		void nonfft_multiply_poly_poly_polymod_coeff_smallmod_inplace(const std::uint64_t *operand1, const std::uint64_t *operand2, const PolyModulus &poly_modulus, const SmallModulus &modulus, std::uint64_t *result, MemoryPool &pool);

		// Not currently used
		bool try_invert_poly_coeff_smallmod(const std::uint64_t *operand, const std::uint64_t *poly_modulus, int coeff_count, const SmallModulus &modulus, std::uint64_t *result, MemoryPool &pool);
	}
}