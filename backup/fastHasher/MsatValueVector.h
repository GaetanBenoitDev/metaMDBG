#ifndef MsatValueVector_h
#define MsatValueVector_h

#include <cstdint>
#include <vector>

namespace MBG
{

class MsatValueVector
{
	class MsatValueChunk
	{
	public:
		MsatValueChunk();
		MsatValueChunk(const MsatValueChunk& other);
		MsatValueChunk(MsatValueChunk&& other);
		MsatValueChunk& operator=(const MsatValueChunk& other);
		MsatValueChunk& operator=(MsatValueChunk&& other);
		~MsatValueChunk();
		uint16_t get(uint8_t index) const;
		size_t size() const;
		void set(uint8_t index, uint16_t value);
		void erase(uint8_t index);
	private:
		size_t capacity() const;
		uint64_t filledIndices;
		uint64_t* values;
	};
public:
	MsatValueVector() = default;
	uint16_t get(size_t index) const;
	void set(size_t index, uint16_t val);
	void resize(size_t newSize);
	void erase(size_t index);
	void printSizeInfo() const;
private:
	std::vector<MsatValueChunk> chunks;
};

}

#endif
