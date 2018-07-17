#pragma once

#include "memorypoolhandle.h"
#include "util/exring.h"

namespace seal
{
	namespace util
	{
		/**
		NTT for sequences of elements in extension rings Z_{p^r}[x]/(f(x)).
		Negacyclic NTT requires finding a 2n-th primitive root of unity. This version does not support automatic searching for root of unity. 
		Instead of requiring users to manually supply the root of unity, we provide several parameter settings here such that the 2n-th primitive
		root of unity for these settings are simply 'x'.
		
		n: The length of the sequence to transform
		p: p in the extension ring
		r: r in the extension ring
		k: Degree of the irreducible polynomial f(x)

		/-----------------------------------------------------------------------------------------------------\
		| n     | p                    | r | k  | f(x)                         | Slots  | Element Bit Length  |
		| ------|----------------------|---|----|------------------------------|--------|---------------------|
		| 128   | 281 (10 bits)        | 1 | 2  | 1x^2 + 15                     | 64     | 20                 |
		| 128   | C1 (8 bits)          | 1 | 4  | 1x^4 + B                      | 32     | 32                 |
		| 128   | C1 (8 bits)          | 2 | 4  | 1x^4 + 563C                   | 32     | 32                 |
		
		| 2048  | DA19AE01 (32 bits)   | 1 | 8  | 1x^8 + 37595A                 | 256    | 256                |
		| 2048  | B201 (16 bits)       | 2 | 8  | 1x^8 + 36227D7D               | 256    | 128                |
		| 2048  | C101 (16 bits)       | 1 | 16 | 1x^16 + 392                   | 128    | 256                |
		| 2048  | 101 (9 bits)         | 1 | 16 | 1x^16 + 3                     | 128    | 144                |
		| 2048  | 1E01 (13 bits)       | 1 | 8  | 1x^8 + 3E                     | 256    | 104                |
		| 2048  | DFF (12 bits)        | 1 | 8  | 1x^8 + 1Ax^4 + DFE            | 256    | 96                 |
		| 2048  | 3401 (14 bits)       | 1 | 4  | 1x^4 + 7                      | 512    | 56                 |
		| 2048  | 13FF (13 bits)       | 1 | 4  | 1x^4 + 7x^2 + 13FE            | 512    | 52                 |

		| 4096  | 1E01 (13 bits)	   | 1 | 16 | 1x^16 + 3E                    | 256    | 208                |
		| 4096  | 3401 (14 bits)       | 1 | 8  | 1x^8 + 7                      | 512    | 112                |
		| 4096  | 17FF (13 bits)       | 1 | 4  | 1x^4 + 4x^2 + 17FE            | 1024   | 52                 |
		| 4096  | 4801 (15 bits)       | 1 | 4  | 1x^4 + 13                     | 1024   | 60                 |

		| 8192  | C401 (16 bits)       | 1 | 16 | 1x^16 + 23                    | 512    | 256                |
		| 8192  | E801 (16 bits)       | 1 | 8  | 1x^8 + 3                      | 1024   | 128                |
		| 8192  | 101 (9 bits)         | 1 | 64 | 1x^64 + 3                     | 128    | 576                |

		| 16384 | 83FF (16 bits)       | 1 | 32 | 1x^32 + 1Ax^16 + 83FE         | 512    | 512                |
		| 16384 | F7FF (16 bits)       | 1 | 16 | 1x^16 + 25x^8 + F7FE          | 1024   | 256                |

		| 32768 | F7FF (16 bits)       | 1 | 32 | 1x^32 + 25x^16 + F7FE         | 1024   | 512                |
		\-----------------------------------------------------------------------------------------------------/

		Special case: when f(x) = x, r = 1, and p=1(mod 2n), it falls back to the integer case, which is also handled in this class.
		/------------------\
		| n     | p        |
		|------------------|
		| 4096  | A001     |
		| 8192  | 820001   |
		\------------------/

		*/
		class ExNTT
		{
		public:
			ExNTT(std::shared_ptr<ExRing> ex_ring, int coeff_count_power);

			void ntt_negacyclic(std::vector<ExRingElement> &sequence);

			void inverse_ntt_negacyclic(std::vector<ExRingElement> &sequence);

		private:
			bool generate();

			ExRingElement compute_primitive_root();

			ExRingElement inv_primitive_root(ExRingElement root);

			void ntt_powers_of_primitive_root(const ExRingElement &root, std::vector<ExRingElement> &powers);

			std::shared_ptr<ExRing> ex_ring_;

			int coeff_count_power_;

			int coeff_count_;

			Pointer root_powers_backing_;

			Pointer inv_root_powers_backing_;
			
			std::vector<ExRingElement> root_powers_;

			std::vector<ExRingElement> inv_root_powers_;

			BigUInt inv_n_;

		};
	}
}