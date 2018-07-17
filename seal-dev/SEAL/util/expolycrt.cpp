#include "util/expolycrt.h"
#include "util/uintcore.h"
#include "util/polycore.h"

using namespace std;

namespace seal
{
    namespace util
    {
        ExPolyCRTBuilder::ExPolyCRTBuilder(shared_ptr<ExRing> ex_ring, int coeff_count_power)
            : ex_ring_(move(ex_ring)), coeff_count_power_(coeff_count_power), 
            coeff_count_(1 << coeff_count_power), cosets_(1 << coeff_count_power),
            exntt_(ex_ring_, coeff_count_power_)
        {
            if (coeff_count_ % (ex_ring_->coeff_count() - 1) != 0)
            {
                throw invalid_argument("coeff_count_ is not divisible by degree of ex_ring poly modulus.");
            }
            data_slots_ = coeff_count_ / (ex_ring_->coeff_count() - 1);

            compute_coset();
        }

        void ExPolyCRTBuilder::compute_coset()
        {
            for (int i = 0; i < coeff_count_; i++)
            {
                cosets_[i].odd = 2 * i + 1;
                cosets_[i].rep = 2 * i + 1;
                cosets_[i].hop = 0;
            }
            int idx = 0;
            for (int i = 0; i < coeff_count_; i++)
            {
                if (cosets_[i].rep < (2 * i + 1))
                {
                    continue;
                }
                data_idx_[2 * i + 1] = idx++;
                int j = (((ex_ring_->characteristic() * (2 * i + 1)) % (2 * coeff_count_)).pointer()[0] - 1) / 2; // 2*j+1 = (p*(2*i+1)) % m, where m = 2*n
                int hop = 1;
                while (cosets_[j].rep != (2 * i + 1))
                {
                    cosets_[j].rep = 2 * i + 1;
                    cosets_[j].hop = hop;
                    j = (((ex_ring_->characteristic() * (2 * j + 1)) % (2 * coeff_count_)).pointer()[0] - 1) / 2;
                    hop += 1;
                }
            }
        }

        void ExPolyCRTBuilder::bit_reversal(vector<ExRingElement> &values)
        {
            for (int i = 0; i < coeff_count_; i++)
            {
                int reversed_i = reverse_bits(i) >> (32 - coeff_count_power_);
                if (i < reversed_i)
                {
                    swap(values[i], values[reversed_i]);
                }
            }
        }

        void ExPolyCRTBuilder::compose(const vector<ExRingElement> &values, Plaintext &destination)
        {
            Pointer expanded_backing;
            vector<ExRingElement> expanded(ex_ring_->allocate_elements(coeff_count_, expanded_backing));
            expand(values, expanded);
            bit_reversal(expanded);
            exntt_.inverse_ntt_negacyclic(expanded);
#ifdef _DEBUG
            for (int i = 0; i < coeff_count_; i++)
            {
                if (get_significant_coeff_count_poly(expanded[i].pointer(), ex_ring_->coeff_count(), ex_ring_->coeff_uint64_count()) > 1)
                {
                    throw logic_error("composed element not in integer ring");
                }
            }
#endif
            BigPoly &destination_poly = destination.get_poly();
            if (destination_poly.coeff_count() != (coeff_count_ + 1) || 
                destination_poly.coeff_bit_count() != ex_ring_->coeff_uint64_count() * bits_per_uint64)
            {
                destination_poly.resize(coeff_count_ + 1, ex_ring_->coeff_uint64_count() * bits_per_uint64);
            }
            destination_poly.set_zero();
            
            for (int i = 0; i < coeff_count_; i++)
            {
                destination[i] = BigUInt(expanded[i].ex_ring()->coeff_uint64_count() * bits_per_uint64, expanded[i].pointer(0));
            }
        }

        void ExPolyCRTBuilder::decompose(const Plaintext &plain, vector<ExRingElement> &destination)
        {
            Pointer expanded_backing;
            vector<ExRingElement> expanded = ex_ring_->allocate_elements(coeff_count_, expanded_backing);
            
            for (int i = 0; i < coeff_count_; i++)
            {
                set_uint_uint(plain.get_poly().pointer(i), plain.get_poly().coeff_uint64_count(), ex_ring_->coeff_uint64_count(), expanded[i].pointer(0));
            }

            exntt_.ntt_negacyclic(expanded);
            bit_reversal(expanded);
            contract(expanded, destination);
        }

        void ExPolyCRTBuilder::expand(const vector<ExRingElement> &data, vector<ExRingElement> &expanded)
        {
            if (data.size() != data_slots_)
            {
                throw invalid_argument("invalid number of data values.");
            }

            expanded.resize(coeff_count_);

            for (int i = 0; i < coeff_count_; i++)
            {
                int power = 2 * i + 1,
                    reduced_power = cosets_[i].rep,
                    frob_index = cosets_[i].hop;
                ex_ring_->frobenius(data[data_idx_[reduced_power]], frob_index, expanded[i]);
            }
        }

        void ExPolyCRTBuilder::contract(const vector<ExRingElement> &expanded, vector<ExRingElement> &data)
        {
            if (expanded.size() != coeff_count_)
            {
                throw invalid_argument("invalid number of expanded values.");
            }

            data.resize(data_slots_);

            for (int i = 0; i < coeff_count_; i++)
            {
                if (cosets_[i].rep == (2 * i + 1))
                {
                    data[data_idx_[2 * i + 1]] = expanded[i];
                }
            }
        }
    }
}