#include "util/exring.h"
#include "util/uintextras.h"
#include "util/uintcore.h"
#include "util/polyarith.h"
#include "util/polyfftmult.h"
#include "util/polyfftmultmod.h"
#include "util/polyarithmod.h"
#include "util/polycore.h"
#include "util/uintarith.h"
#include <random>
#include "util/defines.h"
#include "util/polyarithsmallmod.h"
#include "util/uintarithsmallmod.h"
#include "util/smallpolyarith.h"

using namespace std;

namespace seal 
{
    namespace util 
    {
        namespace 
        {
            /**
            This function is almost a direct copy of "exponentiate_poly_polymod_coeffmod" in polyextras.cpp.
            The only difference is that I change the line
            "if (poly_modulus.is_fft_modulus())"
            to
            "if (poly_modulus.is_fft_modulus() && poly_modulus.coeff_count_power_of_two() > 1)"
            in order to handle the case when the polynomial modulus is "x^2 + 1".
            */
            void exring_exponentiate(const uint64_t *poly, const uint64_t *exponent, int exponent_uint64_count,
                const PolyModulus &poly_modulus, const Modulus &coeff_modulus, uint64_t *result, MemoryPool &pool)
            {
                int coeff_modulus_uint64_count = coeff_modulus.uint64_count();
                int poly_modulus_coeff_count = poly_modulus.coeff_count();
#ifdef _DEBUG
                int poly_modulus_coeff_uint64_count = poly_modulus.coeff_uint64_count();
                if (poly == nullptr)
                {
                    throw invalid_argument("poly");
                }
                if (exponent == nullptr)
                {
                    throw invalid_argument("exponent");
                }
                if (exponent_uint64_count <= 0)
                {
                    throw invalid_argument("exponent_uint64_count");
                }
                if (is_zero_uint(coeff_modulus.get(), coeff_modulus_uint64_count))
                {
                    throw invalid_argument("coeff_modulus");
                }
                if (is_zero_poly(poly_modulus.get(), poly_modulus_coeff_count, poly_modulus_coeff_uint64_count))
                {
                    throw invalid_argument("poly_modulus");
                }
                if (result == nullptr)
                {
                    throw invalid_argument("result");
                }
#endif
                // Fast cases
                if (is_zero_uint(exponent, exponent_uint64_count))
                {
                    set_zero_poly(poly_modulus_coeff_count, coeff_modulus_uint64_count, result);
                    *result = 1;
                    return;
                }

                modulo_poly(poly, poly_modulus_coeff_count, poly_modulus, coeff_modulus, result, pool);

                if (is_equal_uint(exponent, exponent_uint64_count, 1))
                {
                    return;
                }

                // Need to make a copy of exponent
                Pointer exponent_copy(allocate_uint(exponent_uint64_count, pool));
                set_uint_uint(exponent, exponent_uint64_count, exponent_copy.get());

                // Perform binary exponentiation.
                Pointer alloc_anchor(allocate_poly(poly_modulus_coeff_count + poly_modulus_coeff_count + poly_modulus_coeff_count, coeff_modulus_uint64_count, pool));

                uint64_t *powerptr = alloc_anchor.get();
                uint64_t *productptr = get_poly_coeff(powerptr, poly_modulus_coeff_count, coeff_modulus_uint64_count);
                uint64_t *intermediateptr = get_poly_coeff(productptr, poly_modulus_coeff_count, coeff_modulus_uint64_count);

                set_poly_poly(result, poly_modulus_coeff_count, coeff_modulus_uint64_count, powerptr);
                set_zero_poly(poly_modulus_coeff_count, coeff_modulus_uint64_count, intermediateptr);
                *intermediateptr = 1;

                // Initially: power = operand and intermediate = 1, product is not initialized.
                if (poly_modulus.is_fft_modulus() && poly_modulus.coeff_count_power_of_two() > 1)
                {
                    while (true)
                    {
                        if ((*exponent_copy.get() % 2) == 1)
                        {
                            nussbaumer_multiply_poly_poly_coeffmod(powerptr, intermediateptr, poly_modulus.coeff_count_power_of_two(), coeff_modulus, productptr, pool);
                            swap(productptr, intermediateptr);
                        }
                        right_shift_uint(exponent_copy.get(), 1, exponent_uint64_count, exponent_copy.get());
                        if (is_zero_uint(exponent_copy.get(), exponent_uint64_count))
                        {
                            break;
                        }
                        nussbaumer_multiply_poly_poly_coeffmod(powerptr, powerptr, poly_modulus.coeff_count_power_of_two(), coeff_modulus, productptr, pool);
                        swap(productptr, powerptr);
                    }
                }
                else
                {
                    while (true)
                    {
                        if ((*exponent_copy.get() % 2) == 1)
                        {
                            nonfft_multiply_poly_poly_polymod_coeffmod(powerptr, intermediateptr, poly_modulus, coeff_modulus, productptr, pool);
                            swap(productptr, intermediateptr);
                        }
                        right_shift_uint(exponent_copy.get(), 1, exponent_uint64_count, exponent_copy.get());
                        if (is_zero_uint(exponent_copy.get(), exponent_uint64_count))
                        {
                            break;
                        }
                        nonfft_multiply_poly_poly_polymod_coeffmod(powerptr, powerptr, poly_modulus, coeff_modulus, productptr, pool);
                        swap(productptr, powerptr);
                    }
                }
                set_poly_poly(intermediateptr, poly_modulus_coeff_count, coeff_modulus_uint64_count, result);
            }

            /* Returns (operand1 + operand2) mod modulus. Assuming the two operands are already reduced. */
            uint64_t fast_add_mod(uint64_t operand1, uint64_t operand2, uint64_t modulus)
            {
                if (operand1 >= modulus - operand2)
                    operand2 -= modulus;
                return operand1 + operand2;
            }

            /* Return (operand1 - operand2) mod modulus. Assuming the two operands are already reduced. */
            uint64_t fast_sub_mod(uint64_t operand1, uint64_t operand2, uint64_t modulus)
            {
                if (operand1 < operand2)
                    return operand1 + (modulus - operand2);
                return operand1 - operand2;
            }

        }

        /***** ExRing *****/

        ExRing::ExRing(
            const BigUInt &characteristic,
            uint64_t exponent,
            const BigPoly &poly_mod,
            const MemoryPoolHandle &pool)
            : characteristic_(characteristic), exponent_(exponent), pool_(pool), frobe_initialized_(false)
        {
            // 1. Do we need to verify whether "prime" is indeed a prime number? And does "poly_mod" need to be irreducible?
            // 2. There is no way to construct Modulus and PolyModulus without defining another corresponding variable to hold the data;
            // namely, for a Modulus (PolyModulus) member variable, there must be another member variable holding the data. If it is just
            // a local variable, its memory is released after this function, and Modulus will point to dangling address.

            // Calculate coefficient modulus p^r.
            int modulus_uint64_count = divide_round_up(characteristic_.significant_bit_count() * exponent_, bits_per_uint64);
            modulus_pointer_ = allocate_uint(modulus_uint64_count, pool_);
            exponentiate_uint(characteristic_.pointer(), characteristic_.uint64_count(), &exponent_, 1, modulus_uint64_count, modulus_pointer_.get(), pool_);
            coeff_mod_ = Modulus(modulus_pointer_.get(), modulus_uint64_count, pool_);

            // Each coefficient of PolyModulus should have the same size as Modulus.
            int coeff_count = poly_mod.significant_coeff_count();
            poly_modulus_pointer_ = allocate_zero_poly(coeff_count, modulus_uint64_count, pool_);
            set_poly_poly(poly_mod.pointer(), coeff_count, poly_mod.coeff_uint64_count(), coeff_count, modulus_uint64_count, poly_modulus_pointer_.get());
            poly_mod_ = PolyModulus(poly_modulus_pointer_.get(), coeff_count, modulus_uint64_count);

            /* If the modulus is less than 64 bits, we can have faster arithmetic. (64 bits or 2^62 will NOT work) */
            if (coeff_mod_.significant_bit_count() < 64 && coeff_mod_.get()[0] != (1ULL << 62))
            {
                small_mod_ = SmallModulus(coeff_mod_.get()[0]);
            }

            compute_hash();
        }

        void ExRing::init_frobe_table()
        {
            // Construct the table for frobenius computation
            ExRingElement x(get_shared_ptr(), "1x^1"), temp(get_shared_ptr());
            frobe_table_ = allocate_elements(poly_mod_.coeff_count(), poly_mod_.coeff_count(), frobe_table_backing_);
            for (uint64_t k = 1; k < poly_mod_.coeff_count(); k++) // The first vector will never be used.
            {
                /* Compute p^k. (BigUInt should provide an exponentiation operation) */
                BigUInt exp(characteristic_.significant_bit_count() * k);
                exponentiate_uint(characteristic_.pointer(), characteristic_.uint64_count(), &k, 1, exp.uint64_count(), exp.pointer(), pool_);
                for (int i = 0; i < poly_mod_.coeff_count(); i++)
                {
                    // frobe_table_[k][i] = x^(exp * i); // x^{i*p^k} (mod f(x))
                    exponentiate(x, exp * i, frobe_table_[k][i]);
                }
            }
            frobe_initialized_ = true;
        }

		void ExRing::init_frobe_table(int m_)
		{
			// Construct the table for frobenius computation
			ExRingElement x(get_shared_ptr(), "1x^1"), temp(get_shared_ptr());
			frobe_table_ = allocate_elements(poly_mod_.coeff_count(), poly_mod_.coeff_count(), frobe_table_backing_);
			uint64_t base = *(characteristic_.pointer());
			uint64_t exp = base;
			for (uint64_t k = 1; k < poly_mod_.coeff_count(); k++) // The first vector will never be used.
			{
				for (int i = 0; i < poly_mod_.coeff_count(); i++)
				{
					exponentiate(x, (exp * i) % m_, frobe_table_[k][i]);
				}
				exp *= base;
				exp %= m_;
			}
			frobe_initialized_ = true;
		}

        void ExRing::set_frobe_table(const vector<vector<ExRingElement> > &table)
        {
            frobe_table_ = allocate_elements(poly_mod_.coeff_count(), poly_mod_.coeff_count(), frobe_table_backing_);
            for (int k = 1; k < poly_mod_.coeff_count(); k++)
            {
                for (int i = 0; i < poly_mod_.coeff_count(); i++)
                {
                    frobe_table_[k][i] = table[k][i];
                }
            }
            frobe_initialized_ = true;
        }

        void ExRing::add(const ExRingElement &operand1, const ExRingElement &operand2, ExRingElement &result) const
        {
            if (coeff_mod_.uint64_count() == 1)
            {
                add_poly_poly_coeff_smallmod(operand1.pointer(), operand2.pointer(), poly_mod_.coeff_count(), small_mod_, result.pointer());
            }
            else
            {
                add_poly_poly_coeffmod(operand1.pointer(), operand2.pointer(), poly_mod_.coeff_count(), coeff_mod_.get(), coeff_mod_.uint64_count(), result.pointer());
            }
        }

        void ExRing::sub(const ExRingElement &operand1, const ExRingElement &operand2, ExRingElement &result) const
        {
            if (coeff_mod_.uint64_count() == 1)
            {
                sub_poly_poly_coeff_smallmod(operand1.pointer(), operand2.pointer(), poly_mod_.coeff_count(), small_mod_, result.pointer());
            }
            else
            {
                sub_poly_poly_coeffmod(operand1.pointer(), operand2.pointer(), poly_mod_.coeff_count(), coeff_mod_.get(), coeff_mod_.uint64_count(), result.pointer());
            }
        }

        void ExRing::negate(const ExRingElement &operand, ExRingElement &result) const
        {
            if (coeff_mod_.uint64_count() == 1)
            {
                negate_poly_coeff_smallmod(operand.pointer(), poly_mod_.coeff_count(), small_mod_, 1, result.pointer());
            }
            else
            {
                negate_poly_coeffmod(operand.pointer(), poly_mod_.coeff_count(), coeff_mod_.get(), coeff_mod_.uint64_count(), result.pointer());
            }
        }

        void ExRing::multiply(const ExRingElement &operand1, const ExRingElement &operand2, ExRingElement &result) const
        {
            if (poly_mod_.coeff_count() == 2 && poly_mod_.coeff_uint64_count() == 1)  // The fast case: elements are just uint64_t integers.
            {
                multiply_uint64_smallmod(operand1.pointer(0), operand2.pointer(0), small_mod_, result.pointer(0));
            }
            else if (poly_mod_.is_fft_modulus() && poly_mod_.coeff_count_power_of_two() > 1)
            {
                nussbaumer_multiply_poly_poly_coeffmod(operand1.pointer(), operand2.pointer(), poly_mod_.coeff_count_power_of_two(), coeff_mod_, result.pointer(), pool_);
            }
            else
            {
                nonfft_multiply_poly_poly_polymod_coeff_smallmod(operand1.pointer(), operand2.pointer(), poly_mod_, small_mod_, result.pointer(), pool_);
            }
        }

        void ExRing::multiply_constant(const ExRingElement &operand, const uint64_t *constant, ExRingElement &result) const
        {
            int coeff_count = poly_mod_.coeff_count();
            int coeff_uint64_count = poly_mod_.coeff_uint64_count();

            if (coeff_uint64_count == 1) // The fast case
            {
                uint64_t *result_coeff_ptr = result.pointer();
                const uint64_t *operand_coeff_ptr = operand.pointer();
                for (int i = 0; i < coeff_count; i++)
                {
                    if (*operand_coeff_ptr == 0)
                    {
                        *result_coeff_ptr = 0;
                    }
                    else
                    {
                        multiply_uint64_smallmod(operand_coeff_ptr, constant, small_mod_, result_coeff_ptr);
                    }
                    result_coeff_ptr++;
                    operand_coeff_ptr++;
                }
            }
            else
            {
                // Reduce constant modulo coeff_mod_, because "multiply_poly_scalar_coeffmod" assumes the scalar is reduced.
                Pointer reduced_constant(allocate_uint(coeff_uint64_count, pool_));
                modulo_uint(constant, coeff_uint64_count, coeff_mod_, reduced_constant.get(), pool_);

                multiply_poly_scalar_coeffmod(operand.pointer(), coeff_count, reduced_constant.get(), coeff_mod_, result.pointer(), pool_);
            }
        }

        void ExRing::exponentiate(const ExRingElement &base, uint64_t exponent, ExRingElement &result) const
        {
            exring_exponentiate(base.pointer(), &exponent, 1, poly_mod_, coeff_mod_, result.pointer(), pool_);
        }

        void ExRing::exponentiate(const ExRingElement &base, const BigUInt &exponent, ExRingElement &result) const
        {
            exring_exponentiate(base.pointer(), exponent.pointer(), exponent.uint64_count(), poly_mod_, coeff_mod_, result.pointer(), pool_);
        }

        void ExRing::frobenius(const ExRingElement &operand, int repetition, ExRingElement &result)
        {
            if (repetition < 0)
            {
                throw invalid_argument("repetition cannot be negative");
            }
            if (repetition == 0)
            {
                result = operand;
                return;
            }
            else if (repetition >= poly_mod_.coeff_count() || repetition < 0)
            {
                throw invalid_argument("Unsupported repetition of frobenius operations");
            }

            if (!frobe_initialized_)
            {
                throw logic_error("Frobenius table is not intialized");
            }

            ExRingElement temp_result(get_shared_ptr());
            multiply_constant(frobe_table_[repetition][0], operand.pointer(0), temp_result);
            ExRingElement temp(get_shared_ptr());
            for (int i = 1; i < poly_mod_.coeff_count(); i++)
            {
                if (is_zero_uint(operand.pointer(i), coeff_mod_.uint64_count()))
                    continue;
                multiply_constant(frobe_table_[repetition][i], operand.pointer(i), temp); // x^{i*p^k} * a_i, where a_i is the i-th coefficient of the polynomial (operand), k = repetition
                temp_result += temp;
            }

            result = temp_result;
        }

        void ExRing::trace(const ExRingElement &operand, ExRingElement &result)
        {
            int degree = poly_mod_.coeff_count() - 1;
            ExRingElement intermediate = operand;
            result = operand;
            for (int i = 1; i < degree; i++)
            {
                frobenius(operand, i, intermediate);
                result += intermediate;
            }
        }

        ExRingElement ExRing::random_element()
        {
            ExRingElement random(get_shared_ptr());
            uint64_t* ptr = random.pointer();
            random_device rd;
            for (int i = 0; i < coeff_count() - 1; i++)
            {
                uint32_t * random_uint32_ptr = reinterpret_cast<uint32_t*>(ptr);
                for (int j = 0; j < 2 * coeff_uint64_count(); j++)
                {
                    *random_uint32_ptr++ = rd();
                }
                modulo_uint_inplace(ptr, coeff_uint64_count(), coeff_mod_, pool_);
                ptr += coeff_uint64_count();
            }
            return random;
        }

        vector<ExRingElement> ExRing::dyadic_multiply(const vector<ExRingElement> &input1, const vector<ExRingElement> &input2)
        {
            if (input1.size() != input2.size())
            {
                throw invalid_argument("Two input vectors do not have the same length");
            }

            vector<ExRingElement> result(input1.size());
            dyadic_multiply(input1, input2, result);
            return result;
        }

        void ExRing::dyadic_multiply(const vector<ExRingElement> &input1, const vector<ExRingElement> &input2, vector<ExRingElement> &result)
        {
            if (input1.size() != result.size() || input2.size() != result.size())
            {
                throw invalid_argument("Vectors do not have the same length");
            }

            for (int index = 0; index < input1.size(); index++)
            {
                multiply(input1[index], input2[index], result[index]);
            }
        }

        void ExRing::dyadic_square_inplace(vector<ExRingElement> &input)
        {
            for (int index = 0; index < input.size(); index++)
            {
                multiply(input[index], input[index], input[index]);
            }
        }

        bool ExRing::operator ==(const ExRing& other)
        {
            /*if (characteristic_ == other.characteristic_
                && exponent_ == other.exponent_
                && poly_mod_.coeff_count() == other.poly_mod_.coeff_count()
                && is_equal_uint_uint(poly_mod_.get(), other.poly_mod_.get(), poly_mod_.coeff_count() * poly_mod_.coeff_uint64_count()))*/
            return (hash_block_ == other.hash_block_);
        }

        bool ExRing::is_integer_ring()
        {
            int coeff_uint64_count = poly_mod_.coeff_uint64_count();
            return (poly_mod_.coeff_count() == 2)
                && is_equal_uint(poly_mod_.get(), coeff_uint64_count, 0)
                && is_equal_uint(poly_mod_.get() + coeff_uint64_count, coeff_uint64_count, 1);
        }

        vector<ExRingElement> ExRing::allocate_elements(int dim, Pointer &backing)
        {
            vector<ExRingElement> result(dim);
            int uint64_count = element_uint64_count();
            backing = allocate_zero_uint(uint64_count * dim, pool_);
            uint64_t *backing_ptr = backing.get();
            shared_ptr<ExRing> this_ring = get_shared_ptr();
            for (int i = 0; i < dim; i++)
            {
                result[i] = ExRingElement(this_ring, backing_ptr);
                backing_ptr += uint64_count;
            }
            return result;
        }

        vector<vector<ExRingElement> > ExRing::allocate_elements(int dim1, int dim2, Pointer &backing)
        {
            vector<vector<ExRingElement>> result(dim1, vector<ExRingElement>(dim2));
            int uint64_count = element_uint64_count();
            backing = allocate_zero_uint(uint64_count * dim1 * dim2, pool_);
            uint64_t *backing_ptr = backing.get();
            shared_ptr<ExRing> this_ring = get_shared_ptr();
            for (int i = 0; i < dim1; i++)
            {
                for (int j = 0; j < dim2; j++)
                {
                    result[i][j] = ExRingElement(this_ring, backing_ptr);
                    backing_ptr += uint64_count;
                }
            }
            return result;
        }

        vector<vector<vector<ExRingElement> > > ExRing::allocate_elements(int dim1, int dim2, int dim3, Pointer &backing)
        {
            vector<vector<vector<ExRingElement>>> result(dim1, vector<vector<ExRingElement>>(dim2, vector<ExRingElement>(dim3)));
            int uint64_count = element_uint64_count();
            backing = allocate_zero_uint(uint64_count * dim1 * dim2 * dim3, pool_);
            uint64_t *backing_ptr = backing.get();
            shared_ptr<ExRing> this_ring = get_shared_ptr();
            for (int i = 0; i < dim1; i++)
            {
                for (int j = 0; j < dim2; j++)
                {
                    for (int k = 0; k < dim3; k++)
                    {
                        result[i][j][k] = ExRingElement(this_ring, backing_ptr);
                        backing_ptr += uint64_count;
                    }
                }

            }
            return result;
        }

        void ExRing::compute_hash()
        {
            int total_uint64_count = poly_mod_.coeff_uint64_count() * poly_mod_.coeff_count()
                + coeff_mod_.uint64_count();

            Pointer param_data(allocate_uint(total_uint64_count, pool_));
            uint64_t *param_data_ptr = param_data.get();

            set_poly_poly(poly_mod_.get(), poly_mod_.coeff_count(), poly_mod_.coeff_uint64_count(), param_data_ptr);
            param_data_ptr += poly_mod_.coeff_uint64_count() * poly_mod_.coeff_count();
            set_uint_uint(coeff_mod_.get(), coeff_mod_.uint64_count(), param_data_ptr);

            HashFunction::sha3_hash(param_data.get(), total_uint64_count, hash_block_);
        }

        /***** ExRingElement *****/

        ExRingElement::ExRingElement(shared_ptr<ExRing> ex_ring) : 
            ex_ring_(move(ex_ring))
        {
            if (!ex_ring_)
            {
                throw invalid_argument("ExRing is null. Use default constructor instead");
            }
            poly_values_ = allocate_zero_uint(ex_ring_->element_uint64_count(), ex_ring_->pool_);
        }

        ExRingElement::ExRingElement(shared_ptr<ExRing> ex_ring, const string &hex_poly) : 
            ex_ring_(move(ex_ring))
        {
            if (!ex_ring_)
            {
                throw invalid_argument("ExRing is null");
            }
            poly_values_ = allocate_uint(ex_ring_->element_uint64_count(), ex_ring_->pool_);
            BigPoly(ex_ring_->coeff_count(), ex_ring_->coeff_uint64_count() * bits_per_uint64, poly_values_.get()) = hex_poly;
        }

        ExRingElement::ExRingElement(shared_ptr<ExRing> ex_ring, uint64_t *pointer) : 
            ex_ring_(move(ex_ring)), poly_values_(Pointer::Aliasing(pointer))
        {
            if (!ex_ring_)
            {
                throw invalid_argument("ExRing is null");
            }
        }

        ExRingElement::ExRingElement(shared_ptr<ExRing> ex_ring, const string &hex_poly, uint64_t *pointer) : 
            ex_ring_(move(ex_ring)), poly_values_(Pointer::Aliasing(pointer))
        {
            if (!ex_ring_)
            {
                throw invalid_argument("ExRing is null");
            }
            BigPoly(ex_ring_->coeff_count(), ex_ring_->coeff_uint64_count() * bits_per_uint64, poly_values_.get()) = hex_poly;
        }

        ExRingElement::ExRingElement(const ExRingElement &copy)
        {
            operator =(copy);
        }

        ExRingElement& ExRingElement::operator =(const ExRingElement &assign)
        {
            if (assign.empty())
            {
                poly_values_.release();
                ex_ring_.reset();
                return *this;
            }

            if(empty() || *ex_ring_ != *assign.ex_ring_) // Only allocates new memory when the ExRings are different.
            {
                ex_ring_ = assign.ex_ring_;
                poly_values_ = allocate_uint(ex_ring_->element_uint64_count(), ex_ring_->pool_);
            }
            set_uint_uint(assign.poly_values_.get(), ex_ring_->element_uint64_count(), poly_values_.get());

            return *this;
        }

        ExRingElement::ExRingElement(ExRingElement &&source)
        {
            operator =(move(source));
        }

        ExRingElement& ExRingElement::operator =(ExRingElement &&assign)
        {
            if (assign.empty())
            {
                poly_values_.release();
                ex_ring_.reset();
                return *this;
            }

            if(empty() || *ex_ring_ != *assign.ex_ring_) 
            {
                /* A naive move of everything, but it is quite inefficient because of Pointer's thread-safe acquire
                operation in the move, and because of BigPoly's alias operation.
                */

                ex_ring_ = move(assign.ex_ring_);
                poly_values_ = move(assign.poly_values_);
            }
            else
            {
                /* If the two exrings are the same, then it means the size of allocation is also the same. In this case,
                copying the memory is more efficient than moving. */
                set_uint_uint(assign.poly_values_.get(), ex_ring_->element_uint64_count(), poly_values_.get());
            }
            return *this;
        }

        const uint64_t *ExRingElement::pointer(int coeff_index) const
        {
            if (coeff_index < 0 || coeff_index >= ex_ring_->coeff_count())
            {
                throw out_of_range("coeff_index must be within [0, coeff_count)");
            }
            return poly_values_.get() + coeff_index * ex_ring_->coeff_uint64_count();
        }

        uint64_t *ExRingElement::pointer(int coeff_index)
        {
            if (coeff_index < 0 || coeff_index >= ex_ring_->coeff_count())
            {
                throw out_of_range("coeff_index must be within [0, coeff_count)");
            }
            return poly_values_.get() + coeff_index * ex_ring_->coeff_uint64_count();
        }

        bool ExRingElement::operator ==(const ExRingElement &operand2) const
        {
            return (*ex_ring_ == *operand2.ex_ring_) && 
                is_equal_poly_poly(poly_values_.get(), operand2.poly_values_.get(), ex_ring_->coeff_count(), ex_ring_->coeff_uint64_count());
        }

        ExRingElement ExRingElement::operator +(const ExRingElement &operand2) const
        {
            ExRingElement result(ex_ring_);
            ex_ring_->add(*this, operand2, result);
            return result;
        }

        ExRingElement& ExRingElement::operator +=(const ExRingElement &operand2)
        {
            ex_ring_->add(*this, operand2, *this);
            return *this;
        }

        ExRingElement ExRingElement::operator -(const ExRingElement &operand2) const
        {
            ExRingElement result(ex_ring_);
            ex_ring_->sub(*this, operand2, result);
            return result;
        }

        ExRingElement& ExRingElement::operator -=(const ExRingElement &operand2)
        {
            ex_ring_->sub(*this, operand2, *this);
            return *this;
        }

        ExRingElement ExRingElement::operator -() const
        {
            ExRingElement result(ex_ring_);
            ex_ring_->negate(*this, result);
            return result;
        }

        ExRingElement ExRingElement::operator *(const ExRingElement &operand2) const
        {
            ExRingElement result(ex_ring_);
            ex_ring_->multiply(*this, operand2, result);
            return result;
        }

        ExRingElement ExRingElement::operator *(const BigUInt &operand2) const
        {
            ExRingElement result(ex_ring_);
            ex_ring_->multiply_constant(*this, operand2, result);
            return result;
        }

        ExRingElement& ExRingElement::operator *=(const ExRingElement &operand2)
        {
            ex_ring_->multiply(*this, operand2, *this);
            return *this;
        }

        ExRingElement& ExRingElement::operator *=(const BigUInt &operand2)
        {
            ex_ring_->multiply_constant(*this, operand2, *this);
            return *this;
        }

        ExRingElement ExRingElement::operator ^(const uint64_t exponent) const
        {
            ExRingElement result(ex_ring_);
            ex_ring_->exponentiate(*this, exponent, result);
            return result;
        }

        ExRingElement ExRingElement::operator ^(const BigUInt &exponent) const
        {
            ExRingElement result(ex_ring_);
            ex_ring_->exponentiate(*this, exponent, result);
            return result;
        }

        ExRingElement ExRingElement::frobenius() const
        {
            ExRingElement result(ex_ring_);
            ex_ring_->frobenius(*this, result);
            return result;
        }

        ExRingElement ExRingElement::frobenius_k(int k) const
        {
            ExRingElement result(ex_ring_);
            ex_ring_->frobenius(*this, k, result);
            return result;
        }

        ExRingElement ExRingElement::trace() const
        {
            ExRingElement result(ex_ring_);
            ex_ring_->trace(*this, result);
            return result;
        }
    }
}