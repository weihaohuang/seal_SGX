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

        class ExRingElement;

        /**
        Represents an extension ring Z_{p^r}/(f(x)), where f(x) is the polynomial modulus, and each element of this ring is a polynomial 
        with coefficients modulo p^r. p is the charateristic (prime) and r is the exponent.
        
        Question: Is there any requirement about f(x) that a user of this API needs to know?
        */
        class ExRing : public std::enable_shared_from_this<ExRing>
        {
        public:
            typedef HashFunction::sha3_block_type hash_block_type;

            /**
            Creates an extension ring with the specified characteristic, exponent and polynomial modulus. Users are advised to wrap ExRing into 
            a shared_ptr (in order to avoid unexpected behavior of calling get_shared_ptr). Uses acquire_ring instead to get an extension ring.

            @param[in] characteristic Characteristic (p) of the ring Z_{p^r}.
            @param[in] exponent Exponent (r) of the ring Z_{p^r}.
            @param[in] poly_mod Polynomial modulus (f(x)) of the extension ring.
            */
            ExRing(
                const BigUInt &characteristic,
                std::uint64_t exponent,
                const BigPoly &poly_mod,
                const MemoryPoolHandle &pool = MemoryPoolHandle::acquire_global()
            );

            /**
            Static factory function to acquire a shared_ptr of the specified extension ring. This is essentially equivalent to the constructor but
            it returns a shared_ptr which can be passed to construct extension ring elements (ExRingElement).

            @param[in] characteristic Characteristic (p) of the ring Z_{p^r}.
            @param[in] exponent Exponent (r) of the ring Z_{p^r}.
            @param[in] poly_mod Polynomial modulus (f(x)) of the extension ring.
            */
            inline static std::shared_ptr<ExRing> acquire_ring(
                const BigUInt &characteristic,
                std::uint64_t exponent,
                const BigPoly &poly_mod,
                const MemoryPoolHandle &pool = MemoryPoolHandle::acquire_global())
            {
                return std::make_shared<ExRing>(characteristic, exponent, poly_mod, pool);
            }

            void add(const ExRingElement &operand1, const ExRingElement &operand2, ExRingElement &result) const;

            void sub(const ExRingElement &operand1, const ExRingElement &operand2, ExRingElement &result) const;

            void negate(const ExRingElement &operand, ExRingElement &result) const;
            
            void multiply(const ExRingElement &operand1, const ExRingElement &operand2, ExRingElement &result) const;

            inline void multiply_constant(const ExRingElement &operand, const BigUInt &constant, ExRingElement &result) const
            {
                multiply_constant(operand, constant.pointer(), result);
            }

            void exponentiate(const ExRingElement &base, std::uint64_t exponent, ExRingElement &result) const;

            void exponentiate(const ExRingElement &base, const BigUInt &exponent, ExRingElement &result) const;

            void frobenius(const ExRingElement &operand, int repetition, ExRingElement &result);

            /**
            Computes the frobenius of the specified ExRingElement (g, as a polynomial) and returns the result, i.e., frobe(g) = g(x^p), where p is the
            characteristic of the extension ring. The result is reduced modulo the polynomial modulus of the extension ring.
            
            @param[in] operand The specified ExRingElement (g).
            @param[out] result The element where we store the result (frobe(g)).
            */
            void frobenius(const ExRingElement &operand, ExRingElement &result)
            {
                frobenius(operand, 1, result);
            }

            /**
            Computes the trace of the specified ExRingElement (g, as a polynomial) and returns the result, i.e., 
            trace(g) = g + frobe(g) + frobe(frobe(g)) + ... = g + g(x^p) + g(x^{p^2}) + ... + g(x^{p^{d-1}}), 
            where p is the characteristic of the extension ring, and d is the highest degree of the polynomial modulus. 
            The result is reduced modulo the polynomial modulus of the extension ring.

            @param[in] operand The specified ExRingElement (g).
            @param[out] result The element where we store the result (trace(g)).
            */
            void trace(const ExRingElement &operand, ExRingElement &result);

            std::vector<ExRingElement> dyadic_multiply(const std::vector<ExRingElement> &input1, const std::vector<ExRingElement> &input2);

            void dyadic_multiply(const std::vector<ExRingElement> &input1, const std::vector<ExRingElement> &input2, std::vector<ExRingElement> &result);

            void dyadic_square_inplace(std::vector<ExRingElement> &input);

            /**
            Gets a shared_ptr of this object. This should only be called by a previously shared object, i.e., an object managed by shared_ptr.
            Otherwise, the behavior is undefined.
            */
            std::shared_ptr<ExRing> get_shared_ptr()
            {
                return shared_from_this();
            }

            /**
            The total uint64 count to store an element (polynomial) in this extension ring.
            */
            inline int element_uint64_count() const
            {
                return poly_mod_.coeff_count() * poly_mod_.coeff_uint64_count();
            }

            inline int coeff_count() const
            {
                return poly_mod_.coeff_count();
            }

            inline int coeff_uint64_count() const
            {
                return poly_mod_.coeff_uint64_count();
            }

            inline const Modulus &coeff_modulus() const
            {
                return coeff_mod_;
            }

            inline const PolyModulus &poly_modulus() const
            {
                return poly_mod_;
            }

            inline const BigUInt &characteristic() const
            {
                return characteristic_;
            }

            inline std::uint64_t exponent() const
            {
                return exponent_;
            }

            inline MemoryPoolHandle &pool()
            {
                return pool_;
            }

            /**
            Returns a random element in this ring.
            */
            ExRingElement random_element();

            void init_frobe_table();

			void init_frobe_table(int m);

            void set_frobe_table(const std::vector<std::vector<ExRingElement> > &table);

            inline const std::vector<std::vector<ExRingElement> > &frobe_table() const
            {
                return frobe_table_;
            }

            bool operator ==(const ExRing &other);

            bool operator !=(const ExRing &other)
            {
                return !operator==(other);
            }

            /* If f(x)=x, then this is simply an integer ring. */
            bool is_integer_ring();

            /* Returns whether this extension ring is actually a prime integer field. */
            inline bool is_integer_field()
            {
                return is_integer_ring() && exponent_ == 1;
            }

            std::vector<ExRingElement> allocate_elements(int dim, Pointer &backing);

            std::vector<std::vector<ExRingElement> > allocate_elements(int dim1, int dim2, Pointer &backing);

            std::vector<std::vector<std::vector<ExRingElement> > > allocate_elements(int dim1, int dim2, int dim3, Pointer &backing);

        private:
            void multiply_constant(const ExRingElement &operand, const std::uint64_t *constant, ExRingElement &result) const;

            void compute_hash();

            BigUInt characteristic_;

            std::uint64_t exponent_;

            Modulus coeff_mod_;

            Pointer modulus_pointer_;

            SmallModulus small_mod_;

            PolyModulus poly_mod_;

            Pointer poly_modulus_pointer_;

            /* Precomputed frobenius table. frobe_table_[k][i] represents 'x^{i*p^k}'. */
            std::vector<std::vector<ExRingElement> > frobe_table_;

            Pointer frobe_table_backing_;

            bool frobe_initialized_;

            MemoryPoolHandle pool_;

            hash_block_type hash_block_;

            friend class ExRingElement;
        };

        class ExRingElement
        {
        public:
            /**
            Creates an empty element with no information.
            */
            ExRingElement() : ex_ring_(nullptr)
            {
            }

            /**
            Creates a zero element in the extension ring, i.e., a polynomial with coeffcients all equal to 0.

            @param[in] ex_ring The extension ring.
            */
            ExRingElement(std::shared_ptr<ExRing> ex_ring);

            /**
            Creates an extension ring element represented by the hexadecimal polynomial.

            @param[in] ex_ring The extension ring.
            @param[in] hex_poly The hexadecimal polynomial.
            */
            ExRingElement(std::shared_ptr<ExRing> ex_ring, const std::string &hex_poly);

            /**
            Creates a zero element in the extension ring, i.e., a polynomial with coeffcients all equal to 0. Uses the 
            specified pointer for the backing array of this element. Users are responsible to make sure that the pointer 
            points to an array of required size for an element (i.e., the size should be "ex_ring->element_uint64_count()"). 

            @param[in] ex_ring The extension ring.
            @param[in] pointer Pointer to the backing array.
            */
            ExRingElement(std::shared_ptr<ExRing> ex_ring, std::uint64_t *pointer);

            /**
            Creates an extension ring element represented by the hexadecimal polynomial. Uses the specified pointer for
            the backing array of this element. Users are responsible to make sure that the pointer points to an array
            of required size for an element (i.e., the size should be "ex_ring->element_uint64_count()").

            @param[in] ex_ring The extension ring
            @param[in] hex_poly The hexadecimal polynomial
            @param[in] pointer Pointer to the backing array.
            */
            ExRingElement(std::shared_ptr<ExRing> ex_ring, const std::string &hex_poly, std::uint64_t *pointer);

            /**
            Makes a deep copy of the extension ring element.

            @param[in] copy The element to be copied.
            */
            ExRingElement(const ExRingElement &copy);

            /**
            Overwrites the current ExRingElement with the specified one, essentially making a deep copy of the specified one.

            @param[in] assign The element whose values should be assigned to the current ExRingElement.
            */
            ExRingElement& operator =(const ExRingElement &assign);

            ExRingElement(ExRingElement &&source);

            ExRingElement &operator =(ExRingElement &&assign);

            const std::uint64_t *pointer(int coeff_index = 0) const;

            std::uint64_t *pointer(int coeff_index = 0);

            bool operator ==(const ExRingElement &operand2) const;

            /*TODO: Check whether two operands are in the same ring. */

            /**
            Adds two ExRingElement and returns the sum. The input operands are not modified. The result is guaranteed to be reduced modulo
            the polynomial modulus of the extension ring.

            @param[in] operand2 The second operand to add
            */
            ExRingElement operator +(const ExRingElement &operand2) const;

            /**
            Adds two ExRingElement and stores the sum in the current ExRingElement. The result is guaranteed to be reduced modulo
            the polynomial modulus of the extension ring.

            @param[in] operand2 The second operand to add
            */
            ExRingElement &operator +=(const ExRingElement &operand2);

            ExRingElement operator -(const ExRingElement &operand2) const;

            ExRingElement &operator -=(const ExRingElement &operand2);

            /**
            Returns a negated copy of the element.
            */
            ExRingElement operator -() const;

            /**
            Multiplies two ExRingElement and returns the product. The input operands are not modified. The result is guaranteed to be reduced modulo
            the polynomial modulus of the extension ring.

            @param[in] operand2 The second operand to multiply
            */
            ExRingElement operator *(const ExRingElement &operand2) const;

            /**
            Multiplies a constant.
            */
            ExRingElement operator *(const BigUInt &operand2) const;

            /**
            Multiplies two ExRingElement and stores the product in the current ExRingElement. The result is guaranteed to be reduced modulo
            the polynomial modulus of the extension ring.

            @param[in] operand2 The second operand to multiply
            */
            ExRingElement &operator *=(const ExRingElement &operand2);

            /**
            Multiplies a constant and stores the product in the current ExRingElement. The result is guaranteed to be reduced modulo
            the polynomial modulus of the extension ring.

            @param[in] operand2 The second operand to multiply
            */
            ExRingElement &operator *=(const BigUInt &operand2);

            /**
            Exponentiation.

            @param[in] exponent The exponent.
            */
            ExRingElement operator ^(const std::uint64_t exponent) const;

            /**
            Exponentiation.

            @param[in] exponent The exponent.
            */
            ExRingElement operator ^(const BigUInt &exponent) const;

            /**
            Computes the frobenius of the current ExRingElement (g, as a polynomial) and returns the result, i.e., frobe(g) = g(x^p), where p is the 
            characteristic of the extension ring. The current ExRingElement is not modified. The result is guaranteed to be reduced modulo the 
            polynomial modulus of the extension ring.
            */
            ExRingElement frobenius() const;

            ExRingElement frobenius_k(int k) const;

            /**
            Computes the trace of the current ExRingElement (g, as a polynomial) and returns the result, i.e., 
            trace(g) = g + frobe(g) + frobe(frobe(g)) + ... = g + g(x^p) + g(x^{p^2}) + ... + g(x^{p^{d-1}}), 
            where p is the characteristic of the extension ring, and d is the highest degree of the polynomial modulus. The current ExRingElement is 
            not modified. The result is guaranteed to be reduced modulo the polynomial modulus of the extension ring.
            */
            ExRingElement trace() const;

            std::string to_string() const
            {
                return poly_to_hex_string(poly_values_.get(), ex_ring_->coeff_count(), ex_ring_->coeff_uint64_count());
            }

            std::shared_ptr<ExRing> ex_ring()
            {
                return ex_ring_;
            }

            /**
            Returns whether this element is an empty element with no information.
            */
            inline bool empty() const
            {
                return (ex_ring_.get() == nullptr);
            }

        private:
            Pointer poly_values_;

            std::shared_ptr<ExRing> ex_ring_;
        };

    }
}
