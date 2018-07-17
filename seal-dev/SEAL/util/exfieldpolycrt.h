#pragma once

#include "encryptionparams.h"
#include "plaintext.h"
#include "util/exfield.h"
#include "util/exfieldntt.h"
#include <map>

namespace seal
{
    namespace util
    {
        
        /**
        Generalized batching.
        */
        class ExFieldPolyCRTBuilder
        {
        public:
            ExFieldPolyCRTBuilder(std::shared_ptr<ExField> ex_field,
                int coeff_count_power);

            void compose(const std::vector<ExFieldElement> &values, Plaintext &destination);

            Plaintext compose(const std::vector<ExFieldElement> &values)
            {
                Plaintext destination(BigPoly(coeff_count_ + 1, ex_field_->coeff_uint64_count() * bits_per_uint64));
                compose(values, destination);
                return destination;
            }

            void decompose(const Plaintext &plain, std::vector<ExFieldElement> &destination);

            std::vector<ExFieldElement> decompose(const Plaintext &plain)
            {
                std::vector<ExFieldElement> destination;
                decompose(plain, destination);
                return destination;
            }

            int slot_count()
            {
                return data_slots_;
            }

        private:
            void expand(const std::vector<ExFieldElement> &data, std::vector<ExFieldElement> &expanded);

            void contract(const std::vector<ExFieldElement> &expanded, std::vector<ExFieldElement> &data);

            void compute_coset();

            void bit_reversal(std::vector<ExFieldElement> &values);

            std::shared_ptr<ExField> ex_field_;

            int coeff_count_power_;

            int coeff_count_;

            int data_slots_;

            /*
            Equivalence classes. The i-th element corresponds to the odd number (2*i + 1).
            */
            struct coset_element
            {
                int odd; // the odd number 'a'
                int rep; // representative 'b' of the equivalence class
                int hop; // the 'j' in "a = b*(t^j) (mod m)".
            };
            std::vector<coset_element> cosets_;

            /*
            Index map from equivalence class representative to data index.
            */
            std::map<int, int> data_idx_;

            /*
            ExField NTT.
            */
            ExFieldNTT exfieldntt_;

        };


    }
}