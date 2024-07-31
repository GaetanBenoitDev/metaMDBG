#ifndef CompressedSequence_h
#define CompressedSequence_h

#include <vector>
#include <string>
#include <cstdint>
#include "../parallel_hashmap/phmap.h"
#include "StringIndex.h"
#include "TwobitLittleBigVector.h"

namespace MBG
{

class CompressedSequence
{
public:
	CompressedSequence() = default;
	CompressedSequence(TwobitLittleBigVector<uint16_t>&& compressed, std::vector<uint8_t>&& simpleExpanded, std::vector<std::pair<uint32_t, uint32_t>>&& complexExpanded);
	void setCompressed(size_t i, uint16_t c);
	void setCompressedAndClearInputVectorAndResizeExpanded(TwobitLittleBigVector<uint16_t>& compressed); // try not to use this if possible
	uint16_t getCompressed(size_t i) const;
	uint32_t getExpanded(size_t i) const;
	std::string getExpandedStr(size_t i, const StringIndex& index) const;
	void setExpanded(size_t i, uint32_t c);
	size_t compressedSize() const;
	std::vector<size_t> getExpandedPositions(const StringIndex& index) const;
	std::string getExpandedSequence(const StringIndex& index) const;
	CompressedSequence substr(size_t start, size_t len) const;
	CompressedSequence revComp(const StringIndex& index) const;
	std::vector<uint16_t> compressedSubstr(size_t start, size_t len) const;
	void insertEnd(const CompressedSequence& seq);
	void resize(size_t size);
private:
	TwobitLittleBigVector<uint16_t> compressed;
	std::vector<uint8_t> simpleExpanded;
	phmap::flat_hash_map<uint32_t, uint32_t> complexExpanded;
};

}

#endif
