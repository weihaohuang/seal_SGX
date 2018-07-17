#pragma once

#include "util/modulus.h"
#include "util/polymodulus.h"
#include "memorypoolhandle.h"
#include "bigpoly.h"
#include "biguint.h"
#include <memory>
#include <vector>
#include "smallmodulus.h"
#include <array>
#include "hash.h"

namespace seal
{
    namespace util
    {

        class ExFieldElement;

        /**
        Represents an extension field Z_{p}/(f(x)), where f(x) is the polynomial modulus, and each element of this field is a polynomial
        with coefficients modulo p. p is the charateristic (prime) that satisfies the constraint of SmallModulus.
        */
        class ExField : public std::enable_shared_from_this<ExField>
        {
        public:
            typedef HashFunction::sha3_block_type hash_block_type;

            /**
            Creates an extension field with the specified characteristic, exponent and polynomial modulus. Users are advised to wrap ExField into
            a shared_ptr (in order to avoid unexpected behavior of calling get_shared_ptr). Uses acquire_field instead to get an extension field.

            @param[in] characteristic Characteristic (p) of the field Z_{p}.
            @param[in] poly_mod Polynomial modulus (f(x)) of the extension field.
            */
            ExField(
                uint64_t characteristic,
                const BigPoly &poly_mod,
                const MemoryPoolHandle &pool = MemoryPoolHandle::acquire_global()
            );

            /**
            Static factory function to acquire a shared_ptr of the specified extension field. This is essentially equivalent to the constructor but
            it returns a shared_ptr which can be passed to construct extension field elements (ExFieldElement).

            @param[in] characteristic Characteristic (p) of the field Z_{p}.
            @param[in] poly_mod Polynomial modulus (f(x)) of the extension field.
            */
            inline static std::shared_ptr<ExField> acquire_field(
                uint64_t characteristic,
                const BigPoly &poly_mod,
                const MemoryPoolHandle &pool = MemoryPoolHandle::acquire_global())
            {
                return std::make_shared<ExField>(characteristic, poly_mod, pool);
            }

            void add(const ExFieldElement &operand1, const ExFieldElement &operand2, ExFieldElement &result) const;

            void sub(const ExFieldElement &operand1, const ExFieldElement &operand2, ExFieldElement &result) const;

            void negate(const ExFieldElement &operand, ExFieldElement &result) const;

            void multiply(const ExFieldElement &operand1, const ExFieldElement &operand2, ExFieldElement &result) const;

            inline void multiply_constant(const ExFieldElement &operand, const BigUInt &constant, ExFieldElement &result) const
            {
                multiply_constant(operand, constant.pointer(), result);
            }

            void exponentiate(const ExFieldElement &base, std::uint64_t exponent, ExFieldElement &result) const;

            void exponentiate(const ExFieldElement &base, const BigUInt &exponent, ExFieldElement &result) const;

            void frobenius(const ExFieldElement &operand, int repetition, ExFieldElement &result);

            /**
            Computes the frobenius of the specified ExFieldElement (g, as a polynomial) and returns the result, i.e., frobe(g) = g(x^p), where p is the
            characteristic of the extension field. The result is reduced modulo the polynomial modulus of the extension field.

            @param[in] operand The specified ExFieldElement (g).
            @param[out] result The element where we store the result (frobe(g)).
            */
            void frobenius(const ExFieldElement &operand, ExFieldElement &result)
            {
                frobenius(operand, 1, result);
            }

            /**
            Computes the trace of the specified ExFieldElement (g, as a polynomial) and returns the result, i.e.,
            trace(g) = g + frobe(g) + frobe(frobe(g)) + ... = g + g(x^p) + g(x^{p^2}) + ... + g(x^{p^{d-1}}),
            where p is the characteristic of the extension field, and d is the highest degree of the polynomial modulus.
            The result is reduced modulo the polynomial modulus of the extension field.

            @param[in] operand The specified ExFieldElement (g).
            @param[out] result The element where we store the result (trace(g)).
            */
            void trace(const ExFieldElement &operand, ExFieldElement &result);

            std::vector<ExFieldElement> dyadic_multiply(const std::vector<ExFieldElement> &input1, const std::vector<ExFieldElement> &input2);

            void dyadic_multiply(const std::vector<ExFieldElement> &input1, const std::vector<ExFieldElement> &input2, std::vector<ExFieldElement> &result);

            void dyadic_square_inplace(std::vector<ExFieldElement> &input);

            /**
            Gets a shared_ptr of this object. This should only be called by a previously shared object, i.e., an object managed by shared_ptr.
            Otherwise, the behavior is undefined.
            */
            std::shared_ptr<ExField> get_shared_ptr()
            {
                return shared_from_this();
            }

            /**
            The total uint64 count to store an element (polynomial) in this extension field.
            */
            inline int element_uint64_count() const
            {
                return poly_mod_.coeff_count();
            }

            inline int coeff_count() const
            {
                return poly_mod_.coeff_count();
            }

            inline int coeff_uint64_count() const
            {
                return poly_mod_.coeff_uint64_count();
            }

            inline const SmallModulus &coeff_modulus() const
            {
                return small_mod_;
            }

            inline const PolyModulus &poly_modulus() const
            {
                return poly_mod_;
            }

            inline uint64_t characteristic() const
            {
                return small_mod_.value();
            }

            inline MemoryPoolHandle &pool()
            {
                return pool_;
            }

            /**
            Returns a random element in this field.
            */
            ExFieldElement random_element();

            void init_frobe_table();

            void set_frobe_table(const std::vector<std::vector<ExFieldElement> > &table);

            inline const std::vector<std::vector<ExFieldElement> > &frobe_table() const
            {
                return frobe_table_;
            }

            bool operator ==(const ExField &other);

            bool operator !=(const ExField &other)
            {
                return !operator==(other);
            }


            std::vector<ExFieldElement> allocate_elements(int dim, Pointer &backing);

            std::vector<std::vector<ExFieldElement> > allocate_elements(int dim1, int dim2, Pointer &backing);

            std::vector<std::vector<std::vector<ExFieldElement> > > allocate_elements(int dim1, int dim2, int dim3, Pointer &backing);

        private:
            void multiply_constant(const ExFieldElement &operand, const std::uint64_t *constant, ExFieldElement &result) const;

            void compute_hash();

            SmallModulus small_mod_;

            PolyModulus poly_mod_;

            Pointer poly_modulus_pointer_;

            /* Precomputed frobenius table. frobe_table_[k][i] represents 'x^{i*p^k}'. */
            std::vector<std::vector<ExFieldElement> > frobe_table_;

            Pointer frobe_table_backing_;

            bool frobe_initialized_;

            MemoryPoolHandle pool_;

            hash_block_type hash_block_;

            friend class ExFieldElement;
        };

        class ExFieldElement
        {
        public:
            /**
            Creates an empty element with no information.
            */
            ExFieldElement() : ex_field_(nullptr)
            {
            }

            /**
            Creates a zero element in the extension field, i.e., a polynomial with coeffcients all equal to 0.

            @param[in] ex_field The extension field.
            */
            ExFieldElement(std::shared_ptr<ExField> ex_field);

            /**
            Creates an extension field element represented by the hexadecimal polynomial.

            @param[in] ex_field The extension field.
            @param[in] hex_poly The hexadecimal polynomial.
            */
            ExFieldElement(std::shared_ptr<ExField> ex_field, const std::string &hex_poly);

            /**
            Creates a zero element in the extension field, i.e., a polynomial with coeffcients all equal to 0. Uses the
            specified pointer for the backing array of this element. Users are responsible to make sure that the pointer
            points to an array of required size for an element (i.e., the size should be "ex_field->element_uint64_count()").

            @param[in] ex_field The extension field.
            @param[in] pointer Pointer to the backing array.
            */
            ExFieldElement(std::shared_ptr<ExField> ex_field, std::uint64_t *pointer);

            /**
            Creates an extension field element represented by the hexadecimal polynomial. Uses the specified pointer for
            the backing array of this element. Users are responsible to make sure that the pointer points to an array
            of required size for an element (i.e., the size should be "ex_field->element_uint64_count()").

            @param[in] ex_field The extension field
            @param[in] hex_poly The hexadecimal polynomial
            @param[in] pointer Pointer to the backing array.
            */
            ExFieldElement(std::shared_ptr<ExField> ex_field, const std::string &hex_poly, std::uint64_t *pointer);

            /**
            Makes a deep copy of the extension field element.

            @param[in] copy The element to be copied.
            */
            ExFieldElement(const ExFieldElement &copy);

            /**
            Overwrites the current ExFieldElement with the specified one, essentially making a deep copy of the specified one.

            @param[in] assign The element whose values should be assigned to the current ExFieldElement.
            */
            ExFieldElement& operator =(const ExFieldElement &assign);

            ExFieldElement(ExFieldElement &&source);

            ExFieldElement &operator =(ExFieldElement &&assign);

            std::uint64_t operator[](int coeff_index) const;

            std::uint64_t& operator[](int coeff_index);

            const std::uint64_t *pointer(int coeff_index = 0) const;

            std::uint64_t *pointer(int coeff_index = 0);

            bool operator ==(const ExFieldElement &operand2) const;

            /*TODO: Check whether two operands are in the same field. */

            /**
            Adds two ExFieldElement and returns the sum. The input operands are not modified. The result is guaranteed to be reduced modulo
            the polynomial modulus of the extension field.

            @param[in] operand2 The second operand to add
            */
            ExFieldElement operator +(const ExFieldElement &operand2) const;

            /**
            Adds two ExFieldElement and stores the sum in the current ExFieldElement. The result is guaranteed to be reduced modulo
            the polynomial modulus of the extension field.

            @param[in] operand2 The second operand to add
            */
            ExFieldElement &operator +=(const ExFieldElement &operand2);

            ExFieldElement operator -(const ExFieldElement &operand2) const;

            ExFieldElement &operator -=(const ExFieldElement &operand2);

            /**
            Returns a negated copy of the element.
            */
            ExFieldElement operator -() const;

            /**
            Multiplies two ExFieldElement and returns the product. The input operands are not modified. The result is guaranteed to be reduced modulo
            the polynomial modulus of the extension field.

            @param[in] operand2 The second operand to multiply
            */
            ExFieldElement operator *(const ExFieldElement &operand2) const;

            /**
            Multiplies a constant.
            */
            ExFieldElement operator *(const BigUInt &operand2) const;

            /**
            Multiplies two ExFieldElement and stores the product in the current ExFieldElement. The result is guaranteed to be reduced modulo
            the polynomial modulus of the extension field.

            @param[in] operand2 The second operand to multiply
            */
            ExFieldElement &operator *=(const ExFieldElement &operand2);

            /**
            Multiplies a constant and stores the product in the current ExFieldElement. The result is guaranteed to be reduced modulo
            the polynomial modulus of the extension field.

            @param[in] operand2 The second operand to multiply
            */
            ExFieldElement &operator *=(const BigUInt &operand2);

            /**
            Exponentiation.

            @param[in] exponent The exponent.
            */
            ExFieldElement operator ^(const std::uint64_t exponent) const;

            /**
            Exponentiation.

            @param[in] exponent The exponent.
            */
            ExFieldElement operator ^(const BigUInt &exponent) const;

            /**
            Computes the frobenius of the current ExFieldElement (g, as a polynomial) and returns the result, i.e., frobe(g) = g(x^p), where p is the
            characteristic of the extension field. The current ExFieldElement is not modified. The result is guaranteed to be reduced modulo the
            polynomial modulus of the extension field.
            */
            ExFieldElement frobenius() const;

            ExFieldElement frobenius_k(int k) const;

            /**
            Computes the trace of the current ExFieldElement (g, as a polynomial) and returns the result, i.e.,
            trace(g) = g + frobe(g) + frobe(frobe(g)) + ... = g + g(x^p) + g(x^{p^2}) + ... + g(x^{p^{d-1}}),
            where p is the characteristic of the extension field, and d is the highest degree of the polynomial modulus. The current ExFieldElement is
            not modified. The result is guaranteed to be reduced modulo the polynomial modulus of the extension field.
            */
            ExFieldElement trace() const;

            std::string to_string() const
            {
                return poly_to_hex_string(poly_values_.get(), ex_field_->coeff_count(), ex_field_->coeff_uint64_count());
            }

            std::shared_ptr<ExField> ex_field()
            {
                return ex_field_;
            }

            /**
            Returns whether this element is an empty element with no information.
            */
            inline bool empty() const
            {
                return (ex_field_.get() == nullptr);
            }

        private:
            Pointer poly_values_;

            std::shared_ptr<ExField> ex_field_;
        };

    }
}
