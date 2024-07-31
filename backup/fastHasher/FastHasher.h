#ifndef FastHasher_h
#define FastHasher_h

#include <cstdint>
#include <cstring>
#include <algorithm>
#include <vector>
#include <mutex>
#include "MBGCommon.h"

namespace MBG
{

class FastHasher
{
public:
	FastHasher(size_t kmerSize, uint64_t fwHash, uint64_t bwHash);
	FastHasher(size_t kmerSize);
	__attribute__((always_inline))
	inline void addChar(CharType c)
	{
		fwHash = rotlone(fwHash) ^ fwAddHash(c);
		bwHash = rotrone(bwHash) ^ bwAddHash(c);
	}
	__attribute__((always_inline))
	inline void removeChar(CharType c)
	{
		fwHash ^= fwRemoveHash(c);
		bwHash ^= bwRemoveHash(c);
	}
	__attribute__((always_inline))
	inline uint64_t hash() const
	{
		return std::min(fwHash, bwHash);
	}
	__attribute__((always_inline))
	inline uint64_t getFwHash() const
	{
		return fwHash;
	}
	__attribute__((always_inline))
	inline uint64_t getBwHash() const
	{
		return bwHash;
	}
private:
	__attribute__((always_inline))
	inline uint64_t rotlone(uint64_t val) const
	{
		return (val << 1) | (val >> (64-1));
	};
	__attribute__((always_inline))
	inline uint64_t rotrone(uint64_t val) const
	{
		return (val >> 1) | (val << (64-1));
	};
	__attribute__((always_inline))
	inline uint64_t rotlk(uint64_t val) const
	{
		return (val << kmerSize) | (val >> (64-kmerSize));
	};
	__attribute__((always_inline))
	inline uint64_t rotlkmin1(uint64_t val) const
	{
		return (val << (kmerSize-1)) | (val >> (64-(kmerSize-1)));
	};
	__attribute__((always_inline))
	inline uint64_t fwAddHash(CharType c)
	{
		return fwAdd[(int)c];
	}
	__attribute__((always_inline))
	inline uint64_t fwRemoveHash(CharType c)
	{
		return fwRemove[(int)c];
	}
	__attribute__((always_inline))
	inline uint64_t bwAddHash(CharType c)
	{
		return bwAdd[(int)c];
	}
	__attribute__((always_inline))
	inline uint64_t bwRemoveHash(CharType c)
	{
		return bwRemove[(int)c];
	}
	std::vector<uint64_t>& fwAdd;
	std::vector<uint64_t>& fwRemove;
	std::vector<uint64_t>& bwAdd;
	std::vector<uint64_t>& bwRemove;
	uint64_t fwHash;
	uint64_t bwHash;
	size_t kmerSize;
};

}

#endif
