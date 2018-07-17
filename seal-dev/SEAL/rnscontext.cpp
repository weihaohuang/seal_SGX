#include "rnscontext.h"
#include "rnsencryptionparams.h"
#include "util/polycore.h"
#include "util/uintarith.h"
#include "util/modulus.h"
#include "util/polymodulus.h"
#include <stdexcept>
#include "primes.h"

using namespace std;
using namespace seal::util;

namespace seal
{
    RNSEncryptionParameterQualifiers RNSContext::validate()
    {
        qualifiers_ = RNSEncryptionParameterQualifiers();
        int coeff_mod_count = parms_.coeff_mod_count_;
        int coeff_count = parms_.poly_modulus().coeff_count();
        int coeff_modulus_bit_count = 0;
        for (int i = 0; i < coeff_mod_count; i++)
        {
            coeff_modulus_bit_count += parms_.coeff_mod_array()[i].bit_count();
        }

        // Compute the actual value of coeff modulus
        coeff_modulus_.resize(coeff_mod_count * bits_per_uint64);
        Pointer tmp_products_all(allocate_uint(coeff_mod_count, pool_));
        set_zero_uint(coeff_mod_count, coeff_modulus_.pointer());
        *(coeff_modulus_.pointer()) = 1;

        // Compute the product of all coeff moduli
        for (int i = 0; i < coeff_mod_count; i++)
        {
            multiply_uint_uint64(coeff_modulus_.pointer(), coeff_mod_count, parms_.coeff_mod_array()[i].value(), coeff_mod_count, tmp_products_all.get());
            set_uint_uint(tmp_products_all.get(), coeff_mod_count, coeff_modulus_.pointer());
        }

        // Verify required parameters
        if (!parms_.poly_modulus().is_zero() && !is_zero_uint(coeff_modulus_.pointer(), coeff_mod_count)
            && !parms_.plain_modulus().is_zero() && parms_.decomposition_bit_count() >= 0
            && parms_.noise_standard_deviation() >= 0 && parms_.noise_max_deviation() >= 0
            && is_less_than_uint_uint(parms_.plain_modulus().pointer(), parms_.plain_modulus().uint64_count(), coeff_modulus_.pointer(), coeff_mod_count)
            && are_poly_coefficients_less_than(parms_.poly_modulus().pointer(), parms_.poly_modulus().coeff_count(), parms_.poly_modulus().coeff_uint64_count(), coeff_modulus_.pointer(), coeff_mod_count))
        {
            // The parameters look good so far
            qualifiers_.parameters_set = true;
        }
        else
        {
            // Parameters are not valid
            qualifiers_.parameters_set = false;
            return qualifiers_;
        }

        PolyModulus poly_mod(parms_.poly_modulus().pointer(), parms_.poly_modulus().coeff_count(),
            parms_.poly_modulus().coeff_uint64_count());

        // We will additionally require that poly_modulus is of the form x^N+1, where N is a power of two
        if (poly_mod.is_fft_modulus())
        {
            qualifiers_.enable_nussbaumer = true;
        }
        else
        {
            // Parameters are not valid
            qualifiers_.parameters_set = false;
            return qualifiers_;
        }

        int coeff_count_power = poly_mod.coeff_count_power_of_two();
        //Can relinearization be done? Note that also evaluation keys will have to be generated
        if (parms_.decomposition_bit_count() > 0)
        {
        	qualifiers_.enable_relinearization = true;
        	for (int i = 0; i < coeff_mod_count; i++)
        	{
        		if (parms_.decomposition_bit_count() >= 63)
        		{
        			// Parameters are not valid
        			qualifiers_.parameters_set = false;
        			return qualifiers_;
        		}
        	}
        }

        // Can we use NTT with coeff_modulus?
        qualifiers_.enable_ntt = true;
        for (int i = 0; i < coeff_mod_count; i++)
        {
            if (!small_ntt_tables_[i].generate(coeff_count_power, parms_.coeff_mod_array()[i]))
            {
                // Parameters are not valid
                qualifiers_.enable_ntt = false;
                qualifiers_.parameters_set = false;
                return qualifiers_;
            }
        }

        // Can we use batching? (NTT with plain_modulus)
        if (plain_ntt_tables_.generate(coeff_count_power, parms_.plain_modulus()))
        {
            qualifiers_.enable_batching = true;
        }

        base_convertor_ = BaseConvertor(parms_.coeff_mod_array(), coeff_count, coeff_count_power, parms_.plain_modulus(), pool_);
        if (!base_convertor_.is_generated())
        {
            // Parameters are not valid
            qualifiers_.parameters_set = false;
            return qualifiers_;
        }

		qualifiers_.enable_fast_plain_lift = true;
		// Check for plain_lift 
		// If all the small coeff moduli are larger than plain modulus we can lift plain coefficients in RNS form regarding each coefficient modulus
		for (int i = 0; i < coeff_mod_count; i++)
		{
			if (parms_.coeff_mod_array()[i].value() <= parms_.plain_modulus().value())
			{
				qualifiers_.enable_fast_plain_lift = false;
			}
		}

        // Done with validation and pre-computations
        return qualifiers_;
    }

    RNSContext::RNSContext(const RNSEncryptionParameters &parms, const MemoryPoolHandle &pool) :
        pool_(pool), parms_(parms), plain_ntt_tables_(pool_), base_convertor_(pool_)
    {
        int coeff_mod_count = parms_.coeff_mod_array().size();
        small_ntt_tables_.resize(coeff_mod_count, pool_);

        // Set random generator
        if (parms_.random_generator() == nullptr)
        {
            parms_.set_random_generator(UniformRandomGeneratorFactory::default_factory());
        }

        qualifiers_ = validate();
    }
}