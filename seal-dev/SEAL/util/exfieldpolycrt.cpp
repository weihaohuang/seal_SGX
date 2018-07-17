#include "util/exfieldpolycrt.h"
#include "util/uintcore.h"
#include "util/polycore.h"

using namespace std;

namespace seal
{
    namespace util
    {
        ExFieldPolyCRTBuilder::ExFieldPolyCRTBuilder(shared_ptr<ExField> ex_field, int coeff_count_power)
            : ex_field_(move(ex_field)), coeff_count_power_(coeff_count_power),
            coeff_count_(1ULL << coeff_count_power), cosets_(1ULL << coeff_count_power),
            exfieldntt_(ex_field_, coeff_count_power_)
        {
            if (coeff_count_ % (ex_field_->coeff_count() - 1) != 0)
            {
                throw invalid_argument("coeff_count_ is not divisible by degree of ex_field poly modulus.");
            }
            data_slots_ = coeff_count_ / (ex_field_->coeff_count() - 1);

            compute_coset();
        }

        void ExFieldPolyCRTBuilder::compute_coset()
        {
            for (int i = 0; i < coeff_count_; i++)
            {
                cosets_[i].odd = 2 * i + 1;
                cosets_[i].rep = 2 * i + 1;
                cosets_[i].hop = 0;
            }
            int idx = 0;
            uint64_t reduced_characteristic = ex_field_->characteristic() % (2 * coeff_count_);
            for (int i = 0; i < coeff_count_; i++)
            {
                if (cosets_[i].rep < (2 * i + 1))
                {
                    continue;
                }
                data_idx_[2 * i + 1] = idx++;
                int j = (((reduced_characteristic * (2 * i + 1)) % (2 * coeff_count_)) - 1) / 2; // 2*j+1 = (p*(2*i+1)) % m, where m = 2*n
                int hop = 1;
                while (cosets_[j].rep != (2 * i + 1))
                {
                    cosets_[j].rep = 2 * i + 1;
                    cosets_[j].hop = hop;
                    j = (((reduced_characteristic * (2 * j + 1)) % (2 * coeff_count_)) - 1) / 2;
                    hop += 1;
                }
            }
        }

        void ExFieldPolyCRTBuilder::bit_reversal(vector<ExFieldElement> &values)
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

        void ExFieldPolyCRTBuilder::compose(const vector<ExFieldElement> &values, Plaintext &destination)
        {
            Pointer expanded_backing;
            vector<ExFieldElement> expanded(ex_field_->allocate_elements(coeff_count_, expanded_backing));
            expand(values, expanded);
            bit_reversal(expanded);
            exfieldntt_.inverse_ntt_negacyclic(expanded);
#ifdef _DEBUG
            for (int i = 0; i < coeff_count_; i++)
            {
                if (get_significant_coeff_count_poly(expanded[i].pointer(), ex_field_->coeff_count(), ex_field_->coeff_uint64_count()) > 1)
                {
                    throw logic_error("composed element not in integer field");
                }
            }
#endif
            BigPoly &destination_poly = destination.get_poly();
            if (destination_poly.coeff_count() != (coeff_count_ + 1) ||
                destination_poly.coeff_bit_count() != ex_field_->coeff_uint64_count() * bits_per_uint64)
            {
                destination_poly.resize(coeff_count_ + 1, ex_field_->coeff_uint64_count() * bits_per_uint64);
            }
            destination_poly.set_zero();

            for (int i = 0; i < coeff_count_; i++)
            {
                destination[i] = BigUInt(expanded[i].ex_field()->coeff_uint64_count() * bits_per_uint64, expanded[i].pointer(0));
            }
        }

        void ExFieldPolyCRTBuilder::decompose(const Plaintext &plain, vector<ExFieldElement> &destination)
        {
            Pointer expanded_backing;
            vector<ExFieldElement> expanded = ex_field_->allocate_elements(coeff_count_, expanded_backing);

            for (int i = 0; i < coeff_count_; i++)
            {
                set_uint_uint(plain.get_poly().pointer(i), plain.get_poly().coeff_uint64_count(), ex_field_->coeff_uint64_count(), expanded[i].pointer(0));
            }

            exfieldntt_.ntt_negacyclic(expanded);
            bit_reversal(expanded);
            contract(expanded, destination);
        }

        void ExFieldPolyCRTBuilder::expand(const vector<ExFieldElement> &data, vector<ExFieldElement> &expanded)
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
                ex_field_->frobenius(data[data_idx_[reduced_power]], frob_index, expanded[i]);
            }
        }

        void ExFieldPolyCRTBuilder::contract(const vector<ExFieldElement> &expanded, vector<ExFieldElement> &data)
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