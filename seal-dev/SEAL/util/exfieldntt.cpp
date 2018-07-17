#include "util/exfieldntt.h"
#include "util/polyarithmod.h"
#include "util/uintcore.h"
#include "util/uintarithsmallmod.h"
#include "util/uintarithmod.h"

using namespace std;

namespace seal
{
    namespace util
    {
        ExFieldNTT::ExFieldNTT(shared_ptr<ExField> ex_field, int coeff_count_power) :
            ex_field_(move(ex_field)),
            coeff_count_power_(coeff_count_power),
            coeff_count_(1ULL << coeff_count_power)
        {
            /* For efficiency consideration, we allocate big chunks of memory for the roots of powers array,
            instead of making independent allocation for each ExFieldElement.
            */
            root_powers_ = ex_field_->allocate_elements(coeff_count_, root_powers_backing_);
            inv_root_powers_ = ex_field_->allocate_elements(coeff_count_, inv_root_powers_backing_);

            generate();
        }

        bool ExFieldNTT::generate()
        {
            ExFieldElement root = compute_primitive_root();
            if (root.empty())
            {
                return false;
            }
            ntt_powers_of_primitive_root(root, root_powers_);

            ExFieldElement inv_root = inv_primitive_root(root);
            ntt_powers_of_primitive_root(inv_root, inv_root_powers_);

            const SmallModulus& coeff_modulus = ex_field_->coeff_modulus();
            inv_n_.resize(coeff_modulus.uint64_count() * bits_per_uint64);
            inv_n_ = coeff_count_;
            small_modulo_uint_inplace(inv_n_.pointer(), coeff_modulus.uint64_count(), coeff_modulus, ex_field_->pool());
            try_invert_uint_mod(inv_n_.pointer(), coeff_modulus.pointer(), coeff_modulus.uint64_count(), inv_n_.pointer(), ex_field_->pool());

            return true;
        }

        ExFieldElement ExFieldNTT::compute_primitive_root()
        {
            if ((ex_field_->characteristic() - 1) % (2 * coeff_count_) == 0) // p=1(mod 2n)
            {
                Pointer root = allocate_uint(ex_field_->coeff_uint64_count(), ex_field_->pool());
                if (!try_minimal_primitive_root_smallmod(2 * coeff_count_, ex_field_->coeff_modulus(), ex_field_->pool(), root.get()))
                {
                    return ExFieldElement();
                }
                else
                {
                    ExFieldElement e(ex_field_);
                    set_uint_uint(root.get(), ex_field_->coeff_uint64_count(), e.pointer(0));
                    return e;
                }
            }
            else // Else, we assume the user has chosen the parameters specified in the header file and thus the 2n-th primitive root of unity is 'x'.
            {
                return ExFieldElement(ex_field_, "1x^1");
            }
        }

        ExFieldElement ExFieldNTT::inv_primitive_root(ExFieldElement root)
        {
            return (root ^ (2 * coeff_count_ - 1));
        }

        void ExFieldNTT::ntt_powers_of_primitive_root(const ExFieldElement &root, vector<ExFieldElement> &powers)
        {
            powers[0] = ExFieldElement(ex_field_, string("1"));
            for (int i = 1; i < coeff_count_; i++)
            {
                uint32_t reversed_i = (reverse_bits(static_cast<uint32_t>(i)) >> (32 - coeff_count_power_));
                ex_field_->exponentiate(root, reversed_i, powers[i]);
            }
        }

        void ExFieldNTT::ntt_negacyclic(vector<ExFieldElement> &sequence)
        {
            if (sequence.size() != coeff_count_)
            {
                throw invalid_argument("sequence length does not match");
            }

            int n = coeff_count_;
            int t = n;
            ExFieldElement temp(ex_field_); /* Don't move this inside the loop. It is inefficient to allocate new memory for new element. */
            for (int m = 1; m < n; m <<= 1)
            {
                t >>= 1;
                for (int i = 0; i < m; i++)
                {
                    int j1 = 2 * i * t;
                    int j2 = j1 + t - 1;
                    const ExFieldElement &S = root_powers_[m + i];
                    for (int j = j1; j <= j2; j++)
                    {
                        /*U = sequence[j];
                        V = sequence[j + t] * S;
                        sequence[j] = U + V;
                        sequence[j + t] = U - V;*/
                        ex_field_->multiply(sequence[j + t], S, temp);
                        ex_field_->sub(sequence[j], temp, sequence[j + t]);
                        ex_field_->add(sequence[j], temp, sequence[j]);
                    }
                }
            }
        }

        void ExFieldNTT::inverse_ntt_negacyclic(vector<ExFieldElement> &sequence)
        {
            if (sequence.size() != coeff_count_)
            {
                throw invalid_argument("sequence length does not match");
            }

            int n = coeff_count_;
            int t = 1;
            ExFieldElement temp(ex_field_); /* Don't move this inside the loop. It is inefficient to allocate new memory for new element. */
            for (int m = n; m > 1; m >>= 1)
            {
                int j1 = 0;
                int h = m >> 1;
                for (int i = 0; i < h; i++)
                {
                    int j2 = j1 + t - 1;
                    const ExFieldElement& S = inv_root_powers_[h + i];
                    for (int j = j1; j <= j2; j++)
                    {
                        /*
                        U = sequence[j];
                        V = sequence[j + t];
                        sequence[j] = U + V;
                        sequence[j + t] = (U - V) * S;
                        */
                        ex_field_->sub(sequence[j], sequence[j + t], temp);
                        ex_field_->add(sequence[j], sequence[j + t], sequence[j]);
                        ex_field_->multiply(temp, S, sequence[j + t]);
                    }
                    j1 += (t << 1);
                }
                t <<= 1;
            }

            for (int j = 0; j < n; j++)
            {
                sequence[j] *= inv_n_;
            }
        }
    }
}