#pragma once

#include <stdexcept>
#include "memorypoolhandle.h"
#include "smallmodulus.h"
#include "util/mempool.h"

namespace seal
{
    namespace util
    {
        class SmallNTTTables
        {
        public:
            SmallNTTTables(const MemoryPoolHandle &pool);

            SmallNTTTables(int coeff_count_power, const SmallModulus &modulus, const MemoryPoolHandle &pool);

            SmallNTTTables(const SmallNTTTables &copy);

            SmallNTTTables &operator =(const SmallNTTTables &assign);

            SmallNTTTables(SmallNTTTables &&source) = default;

            SmallNTTTables &operator =(SmallNTTTables &&assign) = default;

            inline bool is_generated() const
            {
                return generated_;
            }

            bool generate(int coeff_count_power, const SmallModulus &modulus);

            void reset();

            inline const std::uint64_t *get_root() const
            {
#ifdef _DEBUG
                if (!generated_)
                {
                    throw std::logic_error("tables are not generated");
                }
#endif
                return &root_;
            }

            inline const std::uint64_t *get_from_root_powers(int index) const
            {
#ifdef _DEBUG
                if (index >= coeff_count_)
                {
                    throw std::out_of_range("index");
                }
                if (!generated_)
                {
                    throw std::logic_error("tables are not generated");
                }
#endif
                return root_powers_.get() + index;
            }

            inline const std::uint64_t *get_from_scaled_root_powers(int index) const
            {
#ifdef _DEBUG
                if (index >= coeff_count_)
                {
                    throw std::out_of_range("index");
                }
                if (!generated_)
                {
                    throw std::logic_error("tables are not generated");
                }
#endif
                return scaled_root_powers_.get() + index;
            }

            inline const std::uint64_t *get_from_inv_root_powers(int index) const
            {
#ifdef _DEBUG
                if (index >= coeff_count_)
                {
                    throw std::out_of_range("index");
                }
                if (!generated_)
                {
                    throw std::logic_error("tables are not generated");
                }
#endif
                return inv_root_powers_.get() + index;
            }

            inline const std::uint64_t *get_from_scaled_inv_root_powers(int index) const
            {
#ifdef _DEBUG
                if (index >= coeff_count_)
                {
                    throw std::out_of_range("index");
                }
                if (!generated_)
                {
                    throw std::logic_error("tables are not generated");
                }
#endif
                return scaled_inv_root_powers_.get() + index;
            }

            inline const std::uint64_t *get_from_inv_root_powers_div_two(int index) const
            {
#ifdef _DEBUG
                if (index >= coeff_count_)
                {
                    throw std::out_of_range("index");
                }
                if (!generated_)
                {
                    throw std::logic_error("tables are not generated");
                }
#endif
                return inv_root_powers_div_two_.get() + index;
            }

            inline const std::uint64_t *get_from_scaled_inv_root_powers_div_two(int index) const {
#ifdef _DEBUG
                if (index >= coeff_count_)
                {
                    throw std::out_of_range("index");
                }
                if (!generated_)
                {
                    throw std::logic_error("tables are not generated");
                }
#endif
                return scaled_inv_root_powers_div_two_.get() + index;
            }

            inline const std::uint64_t *get_inv_degree_modulo() const
            {
#ifdef _DEBUG
                if (!generated_)
                {
                    throw std::logic_error("tables are not generated");
                }
#endif
                return &inv_degree_modulo_;
            }

            inline const SmallModulus &smallmodulus() const
            {
                return modulus_;
            }

            inline int coeff_count_power() const
            {
                return coeff_count_power_;
            }

            inline int coeff_count() const
            {
                return coeff_count_;
            }

        private:
            // Computed bit-scrambled vector of first 1 << coeff_count_power powers of a primitive root.
            void ntt_powers_of_primitive_root(const std::uint64_t *root, std::uint64_t *destination) const;

            // Scales the elements of a vector returned by powers_of_primitive_root(...) by word_size/modulus and rounds down.
            void ntt_scale_powers_of_primitive_root(const std::uint64_t *input, std::uint64_t *destination) const;

            MemoryPoolHandle pool_;

            // Size coeff_count_
            Pointer root_powers_;

            // Size coeff_count_
            Pointer scaled_root_powers_;

            // Size coeff_count_
            Pointer inv_root_powers_div_two_;

            // Size coeff_count_
            Pointer scaled_inv_root_powers_div_two_;

            bool generated_;

            int coeff_count_power_;

            int coeff_count_;

            SmallModulus modulus_;

            std::uint64_t root_;

            // Size coeff_count_
            Pointer inv_root_powers_;

            // Size coeff_count_
            Pointer scaled_inv_root_powers_;

            std::uint64_t inv_degree_modulo_;

        };

        void smallntt_negacyclic_harvey_lazy(std::uint64_t *operand, const SmallNTTTables &tables, MemoryPool &pool);

        void smallntt_negacyclic_harvey(std::uint64_t *operand, const SmallNTTTables &tables, MemoryPool &pool);

        void inverse_smallntt_negacyclic_harvey_lazy(std::uint64_t *operand, const SmallNTTTables &tables, MemoryPool &pool);

        void inverse_smallntt_negacyclic_harvey(std::uint64_t *operand, const SmallNTTTables &tables, MemoryPool &pool);
    }
}
