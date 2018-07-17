#pragma once
 
#include <cstdint>
#include <array>

namespace seal
{
    namespace util
    {
        class HashFunction
        {
        public:
            HashFunction() = delete;

            static const int sha3_block_uint64_count = 4;

            typedef std::array<std::uint64_t, sha3_block_uint64_count> sha3_block_type;

            static const sha3_block_type sha3_zero_block;

            static void sha3_hash(const std::uint64_t *input, int uint64_count, sha3_block_type &destination);

            inline static void sha3_hash(std::uint64_t input, sha3_block_type &destination)
            {
                sha3_hash(&input, 1, destination);
            }

        private:
            static const std::uint8_t sha3_round_count = 24;

            static const std::uint8_t sha3_rate_uint64_count = 17; // Rate 1088 = 17 * 64 bits

            static const std::uint8_t sha3_capacity_uint64_count = 8; // Capacity 512 = 8 * 64 bits

            static const std::uint8_t sha3_state_uint64_count = 25; // State size = 1600 = 25 * 64 bits

            typedef std::uint64_t sha3_state_type[5][5];

            static const std::uint8_t rho[24];

            static const std::uint64_t sha3_round_consts[sha3_round_count];

            inline static std::uint64_t rot(std::uint64_t input, std::uint8_t s)
            {
                return (input << s) | (input >> (64 - s));
            }

            static void keccak_1600(sha3_state_type &state);

            static void sponge_absorb(const std::uint64_t sha3_block[sha3_rate_uint64_count], sha3_state_type &state);

            static void sponge_squeeze(const sha3_state_type &sha3_state, sha3_block_type &sha3_block);
        };
    }
}