#include <vector>
#include <cassert>
#include "FastHasher.h"
#include "MBGCommon.h"
#include "ErrorMaskHelper.h"

using namespace MBG;

// pointers so modifying map doesn't invalidate references
phmap::flat_hash_map<size_t, std::vector<uint64_t>*> fwAddsPerKmerSize;
phmap::flat_hash_map<size_t, std::vector<uint64_t>*> bwAddsPerKmerSize;
phmap::flat_hash_map<size_t, std::vector<uint64_t>*> fwRemovesPerKmerSize;
phmap::flat_hash_map<size_t, std::vector<uint64_t>*> bwRemovesPerKmerSize;
std::mutex precalcMutex;
std::vector<uint64_t> charHashes;

__attribute__((always_inline))
inline uint64_t rotlone(uint64_t val, size_t kmerSize)
{
	return (val << 1) | (val >> (64-1));
};

__attribute__((always_inline))
inline uint64_t rotrone(uint64_t val, size_t kmerSize)
{
	return (val >> 1) | (val << (64-1));
};

__attribute__((always_inline))
inline uint64_t rotlk(uint64_t val, size_t kmerSize)
{
	return (val << kmerSize) | (val >> (64-kmerSize));
};

__attribute__((always_inline))
inline uint64_t rotlkmin1(uint64_t val, size_t kmerSize)
{
	return (val << (kmerSize-1)) | (val >> (64-(kmerSize-1)));
};

// https://naml.us/post/inverse-of-a-hash-function/
uint64_t getHash(uint64_t key) {
	key = (~key) + (key << 21); // key = (key << 21) - key - 1;
	key = key ^ (key >> 24);
	key = (key + (key << 3)) + (key << 8); // key * 265
	key = key ^ (key >> 14);
	key = (key + (key << 2)) + (key << 4); // key * 21
	key = key ^ (key >> 28);
	key = key + (key << 31);
	return key;
}

void makePrecalcs(size_t kmerSize)
{
	if (charHashes.size() == 0)
	{
		charHashes.resize(maxCode());
		for (size_t i = 0; i < maxCode(); i++)
		{
			charHashes[i] = getHash(i);
		}
	}
	std::vector<uint64_t>* fwAddPtr = new std::vector<uint64_t>;
	std::vector<uint64_t>* bwAddPtr = new std::vector<uint64_t>;
	std::vector<uint64_t>* fwRemovePtr = new std::vector<uint64_t>;
	std::vector<uint64_t>* bwRemovePtr = new std::vector<uint64_t>;
	fwAddsPerKmerSize[kmerSize] = fwAddPtr;
	bwAddsPerKmerSize[kmerSize] = bwAddPtr;
	fwRemovesPerKmerSize[kmerSize] = fwRemovePtr;
	bwRemovesPerKmerSize[kmerSize] = bwRemovePtr;
	std::vector<uint64_t>& fwAdd = *fwAddPtr;
	std::vector<uint64_t>& bwAdd = *bwAddPtr;
	std::vector<uint64_t>& fwRemove = *fwRemovePtr;
	std::vector<uint64_t>& bwRemove = *bwRemovePtr;
	fwAdd.resize(maxCode());
	fwRemove.resize(maxCode());
	bwAdd.resize(maxCode());
	bwRemove.resize(maxCode());
	for (size_t i = 0; i < maxCode(); i++)
	{
		fwAdd[i] = charHashes[i];
	}
	for (size_t i = 0; i < maxCode(); i++)
	{
		fwRemove[i] = rotlk(charHashes[i], kmerSize);
	}
	for (size_t i = 0; i < maxCode(); i++)
	{
		bwAdd[i] = rotlkmin1(charHashes[(int)complement(i)], kmerSize);
	}
	for (size_t i = 0; i < maxCode(); i++)
	{
		bwRemove[i] = rotrone(charHashes[(int)complement(i)], kmerSize);
	}
}

// weird roundabout way of deallocating the pointers when program quits
class PrecalcDeallocator
{
public:
	PrecalcDeallocator() = default;
	~PrecalcDeallocator()
	{
		std::lock_guard<std::mutex> lock { precalcMutex };
		for (auto pair : fwAddsPerKmerSize) delete pair.second;
		for (auto pair : bwAddsPerKmerSize) delete pair.second;
		for (auto pair : fwRemovesPerKmerSize) delete pair.second;
		for (auto pair : bwRemovesPerKmerSize) delete pair.second;
		fwAddsPerKmerSize.clear();
		bwAddsPerKmerSize.clear();
		fwRemovesPerKmerSize.clear();
		bwRemovesPerKmerSize.clear();
	}
};
PrecalcDeallocator precalcDeallocator;

std::vector<uint64_t>& getFwAdd(size_t kmerSize)
{
	std::lock_guard<std::mutex> lock { precalcMutex };
	if (fwAddsPerKmerSize.count(kmerSize) == 1) return *fwAddsPerKmerSize.at(kmerSize);
	assert((size_t)maxCode() < (size_t)std::numeric_limits<CharType>::max());
	makePrecalcs(kmerSize);
	return *fwAddsPerKmerSize.at(kmerSize);
}

std::vector<uint64_t>& getBwAdd(size_t kmerSize)
{
	std::lock_guard<std::mutex> lock { precalcMutex };
	if (bwAddsPerKmerSize.count(kmerSize) == 1) return *bwAddsPerKmerSize.at(kmerSize);
	assert((size_t)maxCode() < (size_t)std::numeric_limits<CharType>::max());
	makePrecalcs(kmerSize);
	return *bwAddsPerKmerSize.at(kmerSize);
}

std::vector<uint64_t>& getFwRemove(size_t kmerSize)
{
	std::lock_guard<std::mutex> lock { precalcMutex };
	if (fwRemovesPerKmerSize.count(kmerSize) == 1) return *fwRemovesPerKmerSize.at(kmerSize);
	assert((size_t)maxCode() < (size_t)std::numeric_limits<CharType>::max());
	makePrecalcs(kmerSize);
	return *fwRemovesPerKmerSize.at(kmerSize);
}

std::vector<uint64_t>& getBwRemove(size_t kmerSize)
{
	std::lock_guard<std::mutex> lock { precalcMutex };
	if (bwRemovesPerKmerSize.count(kmerSize) == 1) return *bwRemovesPerKmerSize.at(kmerSize);
	assert((size_t)maxCode() < (size_t)std::numeric_limits<CharType>::max());
	makePrecalcs(kmerSize);
	return *bwRemovesPerKmerSize.at(kmerSize);
}

FastHasher::FastHasher(size_t kmerSize, uint64_t fwHash, uint64_t bwHash) :
fwAdd(getFwAdd(kmerSize)),
fwRemove(getFwRemove(kmerSize)),
bwAdd(getBwAdd(kmerSize)),
bwRemove(getBwRemove(kmerSize)),
fwHash(fwHash),
bwHash(bwHash),
kmerSize(kmerSize % 64)
{
}

FastHasher::FastHasher(size_t kmerSize) :
fwAdd(getFwAdd(kmerSize)),
fwRemove(getFwRemove(kmerSize)),
bwAdd(getBwAdd(kmerSize)),
bwRemove(getBwRemove(kmerSize)),
fwHash(0),
bwHash(0),
kmerSize(kmerSize % 64)
{
}
