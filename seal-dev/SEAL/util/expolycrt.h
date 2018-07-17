#pragma once

#include "encryptionparams.h"
#include "plaintext.h"
#include "util/exring.h"
#include "util/exntt.h"
#include <map>

namespace seal
{
    namespace util
    {
        struct coset_element
        {
            int odd; // the odd number 'a'
            int rep; // representative 'b' of the equivalence class
            int hop; // the 'j' in "a = b*(t^j) (mod m)".
        };

        /**
        Generalized batching.
        */
        class ExPolyCRTBuilder
        {
        public:
            ExPolyCRTBuilder(std::shared_ptr<ExRing> ex_ring,
                int coeff_count_power);
        
            void compose(const std::vector<ExRingElement> &values, Plaintext &destination);

            Plaintext compose(const std::vector<ExRingElement> &values)
            {
                Plaintext destination(BigPoly(coeff_count_ + 1, ex_ring_->coeff_uint64_count() * bits_per_uint64));
                compose(values, destination);
                return destination;
            }

            void decompose(const Plaintext &plain, std::vector<ExRingElement> &destination);

            std::vector<ExRingElement> decompose(const Plaintext &plain)
            {
                std::vector<ExRingElement> destination;
                decompose(plain, destination);
                return destination;
            }

            int slot_count()
            {
                return data_slots_;
            }

        private:
            void expand(const std::vector<ExRingElement> &data, std::vector<ExRingElement> &expanded);

            void contract(const std::vector<ExRingElement> &expanded, std::vector<ExRingElement> &data);

            void compute_coset();

            void bit_reversal(std::vector<ExRingElement> &values);

            std::shared_ptr<ExRing> ex_ring_;

            int coeff_count_power_;

            int coeff_count_;

            int data_slots_;

            /*
            Equivalence classes. The i-th element corresponds to the odd number (2*i + 1).
            */
            std::vector<coset_element> cosets_;

            /*
            Index map from equivalence class representative to data index.
            */
            std::map<int, int> data_idx_;

            /*
            ExRing NTT.
            */
            ExNTT exntt_;

        };

        
    }
}