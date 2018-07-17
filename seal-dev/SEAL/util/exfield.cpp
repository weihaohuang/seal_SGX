#include "util/exfield.h"
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
#include "util/polyfftmultsmallmod.h"

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
            void exfield_exponentiate_smallmod(const uint64_t *poly, const uint64_t *exponent, int exponent_uint64_count,
                const PolyModulus &poly_modulus, const SmallModulus &coeff_modulus, uint64_t *result, MemoryPool &pool)
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
                if (is_zero_uint(coeff_modulus.pointer(), coeff_modulus_uint64_count))
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

                smallmodulo_poly(poly, poly_modulus_coeff_count, poly_modulus, coeff_modulus, result, pool);

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
                            nussbaumer_multiply_poly_poly_coeff_smallmod(powerptr, intermediateptr, poly_modulus.coeff_count_power_of_two(), coeff_modulus, productptr, pool);
                            swap(productptr, intermediateptr);
                        }
                        right_shift_uint(exponent_copy.get(), 1, exponent_uint64_count, exponent_copy.get());
                        if (is_zero_uint(exponent_copy.get(), exponent_uint64_count))
                        {
                            break;
                        }
                        nussbaumer_multiply_poly_poly_coeff_smallmod(powerptr, powerptr, poly_modulus.coeff_count_power_of_two(), coeff_modulus, productptr, pool);
                        swap(productptr, powerptr);
                    }
                }
                else
                {
                    while (true)
                    {
                        if ((*exponent_copy.get() % 2) == 1)
                        {
                            nonfft_multiply_poly_poly_polymod_coeff_smallmod(powerptr, intermediateptr, poly_modulus, coeff_modulus, productptr, pool);
                            swap(productptr, intermediateptr);
                        }
                        right_shift_uint(exponent_copy.get(), 1, exponent_uint64_count, exponent_copy.get());
                        if (is_zero_uint(exponent_copy.get(), exponent_uint64_count))
                        {
                            break;
                        }
                        nonfft_multiply_poly_poly_polymod_coeff_smallmod(powerptr, powerptr, poly_modulus, coeff_modulus, productptr, pool);
                        swap(productptr, powerptr);
                    }
                }
                set_poly_poly(intermediateptr, poly_modulus_coeff_count, coeff_modulus_uint64_count, result);
            }

        }

        /***** ExField *****/

        ExField::ExField(
            uint64_t characteristic,
            const BigPoly &poly_mod,
            const MemoryPoolHandle &pool)
            : pool_(pool), frobe_initialized_(false)
        {
            small_mod_ = SmallModulus(characteristic);

            // Each coefficient of PolyModulus should have the same size as Modulus.
            int coeff_count = poly_mod.significant_coeff_count();
            poly_modulus_pointer_ = allocate_zero_poly(coeff_count, 1, pool_);
            set_poly_poly(poly_mod.pointer(), coeff_count, poly_mod.coeff_uint64_count(), coeff_count, 1, poly_modulus_pointer_.get());
            poly_mod_ = PolyModulus(poly_modulus_pointer_.get(), coeff_count, 1);

            compute_hash();
        }

        void ExField::init_frobe_table()
        {
            // Construct the table for frobenius computation
            ExFieldElement x(get_shared_ptr(), "1x^1"), temp(get_shared_ptr());
            frobe_table_ = allocate_elements(poly_mod_.coeff_count(), poly_mod_.coeff_count(), frobe_table_backing_);
            for (uint64_t k = 1; k < poly_mod_.coeff_count(); k++) // The first vector will never be used.
            {
                /* Compute p^k. (BigUInt should provide an exponentiation operation) */
                BigUInt exp(small_mod_.bit_count() * k);
                exponentiate_uint(small_mod_.pointer(), 1, &k, 1, exp.uint64_count(), exp.pointer(), pool_);
                for (int i = 0; i < poly_mod_.coeff_count(); i++)
                {
                    // frobe_table_[k][i] = x^(exp * i); // x^{i*p^k} (mod f(x))
                    exponentiate(x, exp * i, frobe_table_[k][i]);
                }
            }
            frobe_initialized_ = true;
        }

        void ExField::set_frobe_table(const vector<vector<ExFieldElement> > &table)
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

        void ExField::add(const ExFieldElement &operand1, const ExFieldElement &operand2, ExFieldElement &result) const
        {

            add_poly_poly_coeff_smallmod(operand1.pointer(), operand2.pointer(), poly_mod_.coeff_count(), small_mod_, result.pointer());

        }

        void ExField::sub(const ExFieldElement &operand1, const ExFieldElement &operand2, ExFieldElement &result) const
        {

            sub_poly_poly_coeff_smallmod(operand1.pointer(), operand2.pointer(), poly_mod_.coeff_count(), small_mod_, result.pointer());

        }

        void ExField::negate(const ExFieldElement &operand, ExFieldElement &result) const
        {
            negate_poly_coeff_smallmod(operand.pointer(), poly_mod_.coeff_count(), small_mod_, 1, result.pointer());
        }

        void ExField::multiply(const ExFieldElement &operand1, const ExFieldElement &operand2, ExFieldElement &result) const
        {
            if (poly_mod_.coeff_count() == 2 && poly_mod_.coeff_uint64_count() == 1)  // The fast case: elements are just uint64_t integers.
            {
                multiply_uint64_smallmod(operand1.pointer(), operand2.pointer(), small_mod_, result.pointer());
            }
            else if (poly_mod_.is_fft_modulus() && poly_mod_.coeff_count_power_of_two() > 1)
            {
                nussbaumer_multiply_poly_poly_coeff_smallmod(operand1.pointer(), operand2.pointer(), poly_mod_.coeff_count_power_of_two(), small_mod_, result.pointer(), pool_);
            }
            else
            {
                nonfft_multiply_poly_poly_polymod_coeff_smallmod(operand1.pointer(), operand2.pointer(), poly_mod_, small_mod_, result.pointer(), pool_);
            }
        }

        void ExField::multiply_constant(const ExFieldElement &operand, const uint64_t *constant, ExFieldElement &result) const
        {
            int coeff_count = poly_mod_.coeff_count();

            uint64_t *result_coeff_ptr = result.pointer();
            const uint64_t *operand_coeff_ptr = operand.pointer();
            for (int i = 0; i < coeff_count; i++)
            {
                if (*operand_coeff_ptr == 0)
                    *result_coeff_ptr++ = 0;
                else
                    multiply_uint64_smallmod(operand_coeff_ptr, constant, small_mod_, result_coeff_ptr++);
                operand_coeff_ptr++;
            }
        }

        void ExField::exponentiate(const ExFieldElement &base, uint64_t exponent, ExFieldElement &result) const
        {
            exfield_exponentiate_smallmod(base.pointer(), &exponent, 1, poly_mod_, small_mod_, result.pointer(), pool_);
        }

        void ExField::exponentiate(const ExFieldElement &base, const BigUInt &exponent, ExFieldElement &result) const
        {
            exfield_exponentiate_smallmod(base.pointer(), exponent.pointer(), exponent.uint64_count(), poly_mod_, small_mod_, result.pointer(), pool_);
        }

        void ExField::frobenius(const ExFieldElement &operand, int repetition, ExFieldElement &result)
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

            ExFieldElement temp_result(get_shared_ptr());
            multiply_constant(frobe_table_[repetition][0], operand.pointer(), temp_result);
            ExFieldElement temp(get_shared_ptr());
            for (int i = 1; i < poly_mod_.coeff_count(); i++)
            {
                if (operand[i] == 0)
                    continue;
                multiply_constant(frobe_table_[repetition][i], operand.pointer(i), temp); // x^{i*p^k} * a_i, where a_i is the i-th coefficient of the polynomial (operand), k = repetition
                temp_result += temp;
            }

            result = temp_result;
        }

        void ExField::trace(const ExFieldElement &operand, ExFieldElement &result)
        {
            int degree = poly_mod_.coeff_count() - 1;
            ExFieldElement intermediate = operand;
            result = operand;
            for (int i = 1; i < degree; i++)
            {
                frobenius(operand, i, intermediate);
                result += intermediate;
            }
        }

        ExFieldElement ExField::random_element()
        {
            ExFieldElement random(get_shared_ptr());
            random_device rd;
            for (int i = 0; i < coeff_count() - 1; i++)
            {
                random[i] = (uint64_t)rd();
                random[i] <<= 32;
                random[i] = random[i] | (uint64_t)rd();

                small_modulo_uint_inplace(&random[i], 1, small_mod_, pool_);
            }
            return random;
        }

        vector<ExFieldElement> ExField::dyadic_multiply(const vector<ExFieldElement> &input1, const vector<ExFieldElement> &input2)
        {
            if (input1.size() != input2.size())
            {
                throw invalid_argument("Two input vectors do not have the same length");
            }

            vector<ExFieldElement> result(input1.size());
            dyadic_multiply(input1, input2, result);
            return result;
        }

        void ExField::dyadic_multiply(const vector<ExFieldElement> &input1, const vector<ExFieldElement> &input2, vector<ExFieldElement> &result)
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

        void ExField::dyadic_square_inplace(vector<ExFieldElement> &input)
        {
            for (int index = 0; index < input.size(); index++)
            {
                multiply(input[index], input[index], input[index]);
            }
        }

        bool ExField::operator ==(const ExField& other)
        {
            /*if (characteristic_ == other.characteristic_
            && exponent_ == other.exponent_
            && poly_mod_.coeff_count() == other.poly_mod_.coeff_count()
            && is_equal_uint_uint(poly_mod_.get(), other.poly_mod_.get(), poly_mod_.coeff_count() * poly_mod_.coeff_uint64_count()))*/
            return (hash_block_ == other.hash_block_);
        }

        vector<ExFieldElement> ExField::allocate_elements(int dim, Pointer &backing)
        {
            vector<ExFieldElement> result(dim);
            int uint64_count = element_uint64_count();
            backing = allocate_zero_uint(uint64_count * dim, pool_);
            uint64_t *backing_ptr = backing.get();
            shared_ptr<ExField> this_field = get_shared_ptr();
            for (int i = 0; i < dim; i++)
            {
                result[i] = ExFieldElement(this_field, backing_ptr);
                backing_ptr += uint64_count;
            }
            return result;
        }

        vector<vector<ExFieldElement> > ExField::allocate_elements(int dim1, int dim2, Pointer &backing)
        {
            vector<vector<ExFieldElement>> result(dim1, vector<ExFieldElement>(dim2));
            int uint64_count = element_uint64_count();
            backing = allocate_zero_uint(uint64_count * dim1 * dim2, pool_);
            uint64_t *backing_ptr = backing.get();
            shared_ptr<ExField> this_field = get_shared_ptr();
            for (int i = 0; i < dim1; i++)
            {
                for (int j = 0; j < dim2; j++)
                {
                    result[i][j] = ExFieldElement(this_field, backing_ptr);
                    backing_ptr += uint64_count;
                }
            }
            return result;
        }

        vector<vector<vector<ExFieldElement> > > ExField::allocate_elements(int dim1, int dim2, int dim3, Pointer &backing)
        {
            vector<vector<vector<ExFieldElement>>> result(dim1, vector<vector<ExFieldElement>>(dim2, vector<ExFieldElement>(dim3)));
            int uint64_count = element_uint64_count();
            backing = allocate_zero_uint(uint64_count * dim1 * dim2 * dim3, pool_);
            uint64_t *backing_ptr = backing.get();
            shared_ptr<ExField> this_field = get_shared_ptr();
            for (int i = 0; i < dim1; i++)
            {
                for (int j = 0; j < dim2; j++)
                {
                    for (int k = 0; k < dim3; k++)
                    {
                        result[i][j][k] = ExFieldElement(this_field, backing_ptr);
                        backing_ptr += uint64_count;
                    }
                }

            }
            return result;
        }

        void ExField::compute_hash()
        {
            int total_uint64_count = poly_mod_.coeff_count() + 1;

            Pointer param_data(allocate_uint(total_uint64_count, pool_));
            uint64_t *param_data_ptr = param_data.get();

            set_poly_poly(poly_mod_.get(), poly_mod_.coeff_count(), poly_mod_.coeff_uint64_count(), param_data_ptr);
            param_data_ptr += poly_mod_.coeff_count();
            set_uint_uint(small_mod_.pointer(), 1, param_data_ptr);

            HashFunction::sha3_hash(param_data.get(), total_uint64_count, hash_block_);
        }

        /***** ExFieldElement *****/

        ExFieldElement::ExFieldElement(shared_ptr<ExField> ex_field) :
            ex_field_(move(ex_field))
        {
            if (!ex_field_)
            {
                throw invalid_argument("ExField is null. Use default constructor instead");
            }
            poly_values_ = allocate_zero_uint(ex_field_->element_uint64_count(), ex_field_->pool_);
        }

        ExFieldElement::ExFieldElement(shared_ptr<ExField> ex_field, const string &hex_poly) :
            ex_field_(move(ex_field))
        {
            if (!ex_field_)
            {
                throw invalid_argument("ExField is null");
            }
            poly_values_ = allocate_uint(ex_field_->element_uint64_count(), ex_field_->pool_);
            BigPoly(ex_field_->coeff_count(), ex_field_->coeff_uint64_count() * bits_per_uint64, poly_values_.get()) = hex_poly;
        }

        ExFieldElement::ExFieldElement(shared_ptr<ExField> ex_field, uint64_t *pointer) :
            ex_field_(move(ex_field)), poly_values_(Pointer::Aliasing(pointer))
        {
            if (!ex_field_)
            {
                throw invalid_argument("ExField is null");
            }
        }

        ExFieldElement::ExFieldElement(shared_ptr<ExField> ex_field, const string &hex_poly, uint64_t *pointer) :
            ex_field_(move(ex_field)), poly_values_(Pointer::Aliasing(pointer))
        {
            if (!ex_field_)
            {
                throw invalid_argument("ExField is null");
            }
            BigPoly(ex_field_->coeff_count(), ex_field_->coeff_uint64_count() * bits_per_uint64, poly_values_.get()) = hex_poly;
        }

        ExFieldElement::ExFieldElement(const ExFieldElement &copy)
        {
            operator =(copy);
        }

        ExFieldElement& ExFieldElement::operator =(const ExFieldElement &assign)
        {
            if (assign.empty())
            {
                poly_values_.release();
                ex_field_.reset();
                return *this;
            }

            if (empty() || *ex_field_ != *assign.ex_field_) // Only allocates new memory when the ExFields are different.
            {
                ex_field_ = assign.ex_field_;
                poly_values_ = allocate_uint(ex_field_->element_uint64_count(), ex_field_->pool_);
            }
            set_uint_uint(assign.poly_values_.get(), ex_field_->element_uint64_count(), poly_values_.get());

            return *this;
        }

        ExFieldElement::ExFieldElement(ExFieldElement &&source)
        {
            operator =(move(source));
        }

        ExFieldElement& ExFieldElement::operator =(ExFieldElement &&assign)
        {
            if (assign.empty())
            {
                poly_values_.release();
                ex_field_.reset();
                return *this;
            }

            if (empty() || *ex_field_ != *assign.ex_field_)
            {
                /* A naive move of everything, but it is quite inefficient because of Pointer's thread-safe acquire
                operation in the move, and because of BigPoly's alias operation.
                */

                ex_field_ = move(assign.ex_field_);
                poly_values_ = move(assign.poly_values_);
            }
            else
            {
                /* If the two exfields are the same, then it means the size of allocation is also the same. In this case,
                copying the memory is more efficient than moving. */
                set_uint_uint(assign.poly_values_.get(), ex_field_->element_uint64_count(), poly_values_.get());
            }
            return *this;
        }

        uint64_t ExFieldElement::operator[](int coeff_index) const
        {
            if (coeff_index < 0 || coeff_index >= ex_field_->coeff_count())
            {
                throw out_of_range("coeff_index must be within [0, coeff_count)");
            }
            return *(poly_values_.get() + coeff_index);
        }

        uint64_t& ExFieldElement::operator[](int coeff_index)
        {
            if (coeff_index < 0 || coeff_index >= ex_field_->coeff_count())
            {
                throw out_of_range("coeff_index must be within [0, coeff_count)");
            }
            return *(poly_values_.get() + coeff_index);
        }

        const uint64_t *ExFieldElement::pointer(int coeff_index) const
        {
            if (coeff_index < 0 || coeff_index >= ex_field_->coeff_count())
            {
                throw out_of_range("coeff_index must be within [0, coeff_count)");
            }
            return poly_values_.get() + coeff_index;
        }

        uint64_t *ExFieldElement::pointer(int coeff_index)
        {
            if (coeff_index < 0 || coeff_index >= ex_field_->coeff_count())
            {
                throw out_of_range("coeff_index must be within [0, coeff_count)");
            }
            return poly_values_.get() + coeff_index;
        }

        bool ExFieldElement::operator ==(const ExFieldElement &operand2) const
        {
            return (*ex_field_ == *operand2.ex_field_) &&
                is_equal_poly_poly(poly_values_.get(), operand2.poly_values_.get(), ex_field_->coeff_count(), ex_field_->coeff_uint64_count());
        }

        ExFieldElement ExFieldElement::operator +(const ExFieldElement &operand2) const
        {
            ExFieldElement result(ex_field_);
            ex_field_->add(*this, operand2, result);
            return result;
        }

        ExFieldElement& ExFieldElement::operator +=(const ExFieldElement &operand2)
        {
            ex_field_->add(*this, operand2, *this);
            return *this;
        }

        ExFieldElement ExFieldElement::operator -(const ExFieldElement &operand2) const
        {
            ExFieldElement result(ex_field_);
            ex_field_->sub(*this, operand2, result);
            return result;
        }

        ExFieldElement& ExFieldElement::operator -=(const ExFieldElement &operand2)
        {
            ex_field_->sub(*this, operand2, *this);
            return *this;
        }

        ExFieldElement ExFieldElement::operator -() const
        {
            ExFieldElement result(ex_field_);
            ex_field_->negate(*this, result);
            return result;
        }

        ExFieldElement ExFieldElement::operator *(const ExFieldElement &operand2) const
        {
            ExFieldElement result(ex_field_);
            ex_field_->multiply(*this, operand2, result);
            return result;
        }

        ExFieldElement ExFieldElement::operator *(const BigUInt &operand2) const
        {
            ExFieldElement result(ex_field_);
            ex_field_->multiply_constant(*this, operand2, result);
            return result;
        }

        ExFieldElement& ExFieldElement::operator *=(const ExFieldElement &operand2)
        {
            ex_field_->multiply(*this, operand2, *this);
            return *this;
        }

        ExFieldElement& ExFieldElement::operator *=(const BigUInt &operand2)
        {
            ex_field_->multiply_constant(*this, operand2, *this);
            return *this;
        }

        ExFieldElement ExFieldElement::operator ^(const uint64_t exponent) const
        {
            ExFieldElement result(ex_field_);
            ex_field_->exponentiate(*this, exponent, result);
            return result;
        }

        ExFieldElement ExFieldElement::operator ^(const BigUInt &exponent) const
        {
            ExFieldElement result(ex_field_);
            ex_field_->exponentiate(*this, exponent, result);
            return result;
        }

        ExFieldElement ExFieldElement::frobenius() const
        {
            ExFieldElement result(ex_field_);
            ex_field_->frobenius(*this, result);
            return result;
        }

        ExFieldElement ExFieldElement::frobenius_k(int k) const
        {
            /*ExFieldElement result = *this;
            for (int i = 0; i < k; i++)
            ex_field_->frobenius(result, result);*/
            ExFieldElement result(ex_field_);
            ex_field_->frobenius(*this, k, result);
            return result;
        }

        ExFieldElement ExFieldElement::trace() const
        {
            ExFieldElement result(ex_field_);
            ex_field_->trace(*this, result);
            return result;
        }
    }
}