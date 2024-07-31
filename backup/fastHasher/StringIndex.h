#ifndef StringIndex_h
#define StringIndex_h

#include <variant>
#include <vector>
#include <string>
#include "../parallel_hashmap/phmap.h"

namespace MBG
{

class StringIndex
{
public:
	void init(size_t maxCode);
	uint32_t getIndex(uint16_t compressed, std::variant<size_t, std::string> expanded);
	std::string getString(uint16_t compressed, uint32_t index) const;
	void buildReverseIndex();
	uint32_t getReverseIndex(uint16_t compressed, uint32_t index) const;
private:
	std::vector<phmap::flat_hash_map<std::string, uint32_t>> index;
	std::vector<std::vector<std::string>> reverseIndex;
};

}

#endif
