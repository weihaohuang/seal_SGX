#pragma once

#include <iostream>
#include <vector>
#include <map>
#include "bigpolyarray.h"
#include "rnsencryptionparams.h"

namespace seal
{
	/**
	Class to store GaloisKeys keys.
	*/
	class RNSGaloisKeys
	{
	public:
		/**
		Creates an empty set of RNSGaloisKeys. No memory is allocated by this
		constructor.
		*/
		RNSGaloisKeys() = default;

		/**
		Creates a new RNSGaloisKeys by copying an old one.

		@param[in] copy The RNSGaloisKeys to copy from
		*/
		RNSGaloisKeys(const RNSGaloisKeys &copy) = default;

		/**
		Creates a new RNSGaloisKeys by movin an old one.

		@param[in] source The RNSGaloisKeys to move from
		*/
		RNSGaloisKeys(RNSGaloisKeys &&source) = default;

		/**
		Copies an old RNSGaloisKeys to the current one.

		@param[in] assign The RNSGaloisKeys to copy from
		*/
		RNSGaloisKeys &operator =(const RNSGaloisKeys &assign) = default;

		/**
		Moves an old RNSGaloisKeys to the current one.

		@param[in] assign The RNSGaloisKeys to move from
		*/
		RNSGaloisKeys &operator =(RNSGaloisKeys &&assign) = default;

		/**
		Returns the current number of Galois keys.
		*/
		int size() const
		{
			return keys_.size();
		}

		const std::map<int, std::vector<std::pair<BigPolyArray, BigPolyArray> > > &get_data() const
		{
			return keys_;
		}

		const std::vector<std::pair<BigPolyArray, BigPolyArray> > &get_key(int galois_elt) {
			return keys_.at(galois_elt);
		}

		bool has_key(int galois_elt) {
			return keys_.find(galois_elt) != keys_.end();
		}

		inline const RNSEncryptionParameters::hash_block_type &get_hash_block() const
		{
			return rns_hash_block_;
		}

	private:
		/**
		Returns a reference to the vector of galois keys. The user should never
		have a reason to modify the evaluation keys by hand.
		*/
		inline std::map<int, std::vector<std::pair<BigPolyArray, BigPolyArray> > > &get_mutable_data()
		{
			return keys_;
		}

		/**
		Returns a reference to the hash block. The user should never have a reason to
		modify the hash block by hand.

		@see EncryptionParameters for more information about the hash block.
		*/
		inline RNSEncryptionParameters::hash_block_type &get_mutable_hash_block()
		{
			return rns_hash_block_;
		}

		/**
		The map of galois keys (each of type std::vector<std::pair<BigPolyArray, BigPolyArray>>).
		Each key corresponds to one particular rotation of the secret key whose influence in
		a ciphertext is removed using the apply_galois() function of Evaluator.
		*/
		std::map<int, std::vector<std::pair<BigPolyArray, BigPolyArray> > >  keys_;

		RNSEncryptionParameters::hash_block_type rns_hash_block_;

		friend class RNSKeyGenerator;

		friend class RNSEvaluator;

	};
}