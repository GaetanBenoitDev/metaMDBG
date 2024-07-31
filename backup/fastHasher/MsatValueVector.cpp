#include <cassert>
#include <cstring>
#include "MsatValueVector.h"

using namespace MBG;

int popcount(uint64_t x);

MsatValueVector::MsatValueChunk::MsatValueChunk() :
	filledIndices(0),
	values(0)
{
}

MsatValueVector::MsatValueChunk::MsatValueChunk(MsatValueChunk&& other) :
	filledIndices(other.filledIndices),
	values(other.values)
{
	other.filledIndices = 0;
	other.values = 0;
}

MsatValueVector::MsatValueChunk::MsatValueChunk(const MsatValueChunk& other) :
	filledIndices(0),
	values(0)
{
	*this = other;
}

MsatValueVector::MsatValueChunk::~MsatValueChunk()
{
	if (size() > 4) delete [] values;
}

MsatValueVector::MsatValueChunk& MsatValueVector::MsatValueChunk::operator=(MsatValueChunk&& other)
{
	if (size() > 4) delete [] values;
	values = other.values;
	filledIndices = other.filledIndices;
	other.values = 0;
	other.filledIndices = 0;
	return *this;
}

MsatValueVector::MsatValueChunk& MsatValueVector::MsatValueChunk::operator=(const MsatValueChunk& other)
{
	if (size() > 4) delete [] values;
	if (other.size() <= 3)
	{
		values = other.values;
		filledIndices = other.filledIndices;
		return *this;
	}
	values = new uint64_t[other.capacity()/4];
	filledIndices = other.filledIndices;
	memcpy(values, other.values, capacity()/4);
	return *this;
}

void MsatValueVector::MsatValueChunk::set(uint8_t index, uint16_t val)
{
	uint16_t got = get(index);
	if (got == val) return;
	if (got != 65535) erase(index);
	uint64_t fillIndex = 1ull << (uint64_t)index;
	assert((fillIndex & filledIndices) == 0);
	uint64_t mask = fillIndex-1;
	size_t position = popcount(mask & filledIndices);
	if (size() < 4)
	{
		filledIndices |= fillIndex;
		switch(position)
		{
		case 0:
			values = (uint64_t*)((uint64_t)values << 16ull);
			values = (uint64_t*)((uint64_t)values + (uint64_t)val);
			return;
		case 1:
			values = (uint64_t*)((((uint64_t)values << 16ull) & 0xFFFFFFFF00000000ull) + ((uint64_t)values & 0x000000000000FFFFull));
			values = (uint64_t*)((uint64_t)values + ((uint64_t)val << 16ull));
			return;
		case 2:
			values = (uint64_t*)((((uint64_t)values << 16ull) & 0xFFFF000000000000ull) + ((uint64_t)values & 0x00000000FFFFFFFFull));
			values = (uint64_t*)((uint64_t)values + ((uint64_t)val << 32ull));
			return;
		case 3:
			values = (uint64_t*)((uint64_t)values + ((uint64_t)val << 48ull));
			return;
		default:
			assert(false);
		}
	}
	if (size() == 4)
	{
		uint64_t* newValues = new uint64_t[2];
		newValues[0] = (uint64_t)values;
		newValues[1] = 0;
		filledIndices |= fillIndex;
		values = newValues;
	}
	else if (size() == capacity())
	{
		uint64_t* newValues = new uint64_t[(size_t)(capacity()/4)*2];
		for (size_t i = 0; i < capacity()/4; i++)
		{
			newValues[i] = values[i];
		}
		for (size_t i = capacity()/4; i < (size_t)(capacity()/4)*2; i++)
		{
			newValues[i] = 0;
		}
		filledIndices |= fillIndex;
		delete [] values;
		values = newValues;
	}
	else
	{
		filledIndices |= fillIndex;
	}
	size_t wordIndex = position/4;
	size_t wordOffset = position%4;
	for (size_t i = (size()+3)/4-1; i > wordIndex; i--)
	{
		values[i] <<= 16ull;
		values[i] += values[i-1] >> 48ull;
	}
	switch(wordOffset)
	{
	case 0:
		values[wordIndex] = ((values[wordIndex] << 16ull) & 0xFFFFFFFFFFFF0000ull);
		values[wordIndex] += (uint64_t)val;
		return;
	case 1:
		values[wordIndex] = ((values[wordIndex] << 16ull) & 0xFFFFFFFF00000000ull) + (values[wordIndex] & 0x000000000000FFFFull);
		values[wordIndex] += (uint64_t)val << 16ull;
		return;
	case 2:
		values[wordIndex] = ((values[wordIndex] << 16ull) & 0xFFFF000000000000ull) + (values[wordIndex] & 0x00000000FFFFFFFFull);
		values[wordIndex] += (uint64_t)val << 32ull;
		return;
	case 3:
		values[wordIndex] = (values[wordIndex] & 0x0000FFFFFFFFFFFFull);
		values[wordIndex] += (uint64_t)val << 48ull;
		return;
	default:
		assert(false);
	}
}

void MsatValueVector::MsatValueChunk::erase(uint8_t index)
{
	//shouldn't ever be called??
	assert(false);
}

uint16_t MsatValueVector::MsatValueChunk::get(uint8_t index) const
{
	uint64_t checkIndex = (1ull) << (uint64_t)index;
	if ((filledIndices & checkIndex) == 0) return 65535;
	size_t position = popcount((checkIndex-1) & filledIndices);
	if (size() <= 4)
	{
		return (((uint64_t)values) >> (16ull * position)) & 0xFFFFull;
	}
	return ((values[position/4]) >> (16ull * (position % 4))) & 0xFFFFull;
}

size_t MsatValueVector::MsatValueChunk::capacity() const
{
	size_t items = size();
	if (items <= 4) return 4;
	if (items <= 8) return 8;
	if (items <= 16) return 16;
	if (items <= 32) return 32;
	if (items <= 64) return 64;
	assert(false);
}

size_t MsatValueVector::MsatValueChunk::size() const
{
	return popcount(filledIndices);
}

uint16_t MsatValueVector::get(size_t index) const
{
	size_t vecIndex = index / 64;
	size_t vecOffset = index % 64;
	auto result = chunks[vecIndex].get(vecOffset);
	assert(result >= 4);
	return result;
}

void MsatValueVector::set(size_t index, uint16_t val)
{
	assert(val >= 4);
	size_t vecIndex = index / 64;
	size_t vecOffset = index % 64;
	chunks[vecIndex].set(vecOffset, val);
}

void MsatValueVector::resize(size_t size)
{
	chunks.resize((size+63)/64);
}

void MsatValueVector::erase(size_t index)
{
	size_t vecIndex = index / 64;
	size_t vecOffset = index % 64;
	chunks[vecIndex].erase(vecOffset);
}
