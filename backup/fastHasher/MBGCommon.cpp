#include <cassert>
#include <string_view>
#include <algorithm>
#include <vector>
#include "MBGCommon.h"

namespace MBG
{

uint16_t complement(const uint16_t original);

HashType hash(VectorView<uint16_t> sequence)
{
	assert(sequence.size() % 2 == 1);
	size_t half = (sequence.size()+1) / 2;
	VectorView<uint16_t> firstHalf { sequence.data, sequence.startpos, sequence.startpos + half };
	std::vector<uint16_t> secondHalf;
	secondHalf.resize(half);
	for (size_t i = 0; i < half; i++)
	{
		secondHalf[i] = complement(sequence[sequence.size()-1-i]);
	}
	size_t low = std::hash<std::string_view>{}( firstHalf.getView());
	size_t high = std::hash<std::string_view>{}(VectorView<uint16_t> { secondHalf, 0, secondHalf.size() }.getView());
	return (HashType)low + (((HashType)high) << 64);
}

HashType hash(VectorView<uint16_t> sequence, VectorView<uint16_t> reverseSequence)
{
	assert(sequence.size() % 2 == 1);
	assert(sequence.size() == reverseSequence.size());
	size_t half = (sequence.size()+1) / 2;
	VectorView<uint16_t> firstHalf { sequence.data, sequence.startpos, sequence.startpos + half };
	VectorView<uint16_t> secondHalf { reverseSequence.data, reverseSequence.startpos, reverseSequence.startpos + half };
	size_t low = std::hash<std::string_view>{}(firstHalf.getView());
	size_t high = std::hash<std::string_view>{}(secondHalf.getView());
	return (HashType)low + (((HashType)high) << 64);
}

HashType hash(std::vector<uint16_t> sequence)
{
	return hash(VectorView<uint16_t> { sequence, 0, sequence.size() });
}

std::ostream& operator<<(std::ostream& os, HashType t)
{
	if (t == 0)
	{
		os << "0";
		return os;
	}
	std::string decimal;
	while (t != 0)
	{
		decimal += "0123456789"[t % 10];
		t /= 10;
	}
	std::reverse(decimal.begin(), decimal.end());
	os << decimal;
	return os;
}

std::istream& operator>>(std::istream& is, HashType& t)
{
	std::string decimal;
	is >> decimal;
	t = 0;
	for (size_t i = 0; i < decimal.size(); i++)
	{
		t *= 10;
		if (decimal[i] >= '0' && decimal[i] <= '9') t += decimal[i]-'0';
	}
	return is;
}

std::pair<size_t, bool> reverse(std::pair<size_t, bool> pos)
{
	return std::make_pair(pos.first, !pos.second);
}

std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>> canon(std::pair<size_t, bool> from, std::pair<size_t, bool> to)
{
	if (to.first < from.first)
	{
		return std::make_pair(reverse(to), reverse(from));
	}
	if (to.first == from.first && !to.second && !from.second)
	{
		return std::make_pair(reverse(to), reverse(from));
	}
	return std::make_pair(from, to);
}

std::string revCompRaw(const std::string& raw)
{
	std::string result { raw.rbegin(), raw.rend() };
	for (size_t i = 0; i < result.size(); i++)
	{
		switch (result[i])
		{
			case 'a':
			case 'A':
				result[i] = 'T';
				break;
			case 'c':
			case 'C':
				result[i] = 'G';
				break;
			case 'g':
			case 'G':
				result[i] = 'C';
				break;
			case 't':
			case 'T':
				result[i] = 'A';
				break;
			default:
				assert(false);
				break;
		}
	}
	return result;
}

std::vector<std::pair<size_t, bool>> revCompPath(const std::vector<std::pair<size_t, bool>>& original)
{
	std::vector<std::pair<size_t, bool>> result = original;
	std::reverse(result.begin(), result.end());
	for (size_t i = 0; i < result.size(); i++)
	{
		result[i] = reverse(result[i]);
	}
	return result;
}

}
