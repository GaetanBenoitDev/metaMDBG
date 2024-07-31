#include <algorithm>
#include <limits>
#include <cassert>
#include <cmath>
#include "ErrorMaskHelper.h"

using namespace MBG;

namespace MBG
{

constexpr size_t MaxMotifLength = 6;

std::vector<uint16_t> calculateNumsBefore()
{
	std::vector<uint16_t> lengths;
	for (size_t motifLength = 0; motifLength <= MaxMotifLength+1; motifLength++)
	{
		size_t result = 0;
		for (uint16_t i = 1; i < motifLength; i++)
		{
			result += pow(4, i) * i;
		}
		assert(result < std::numeric_limits<uint16_t>::max());
		lengths.push_back(result);
	}
	return lengths;
}

std::vector<uint16_t> numsBefore = calculateNumsBefore();

uint16_t getNumBefore(uint16_t motifLength)
{
	return numsBefore[motifLength];
}

size_t codeMotifLength(const uint16_t code)
{
	uint16_t motifLength = 1;
	while (getNumBefore(motifLength+1) <= code)
	{
		motifLength += 1;
		assert(motifLength <= MaxMotifLength);
	}
	return motifLength;
}

CharType getReverseComplement(CharType code)
{
	uint16_t motifLength = 1;
	while (getNumBefore(motifLength+1) <= code)
	{
		motifLength += 1;
		assert(motifLength <= MaxMotifLength);
	}
	code -= getNumBefore(motifLength);
	assert(code < pow(4, motifLength) * motifLength);
	uint16_t motif = code / motifLength;
	assert(motif < pow(4, motifLength));
	uint16_t overhang = code % motifLength;
	assert(overhang < motifLength);
	if (overhang > 0)
	{
		// rotate to start with the overhang
		motif = (motif >> ((motifLength - overhang) * 2)) + ((motif & ((1 << ((motifLength - overhang) * 2)) - 1)) << ((overhang) * 2));
	}
	assert(motif < pow(4, motifLength));
	uint16_t newMotif = 0;
	for (uint16_t i = 0; i < motifLength; i++)
	{
		newMotif <<= 2;
		newMotif |= (~motif) & 3;
		motif >>= 2;
	}
	assert(newMotif < pow(4, motifLength));
	CharType newCode = getNumBefore(motifLength) + newMotif * motifLength + overhang;
	return newCode;
}

std::vector<CharType> getReverseComplements()
{
	SequenceCharType result;
	for (uint16_t motifLength = 1; motifLength <= MaxMotifLength; motifLength++)
	{
		for (uint16_t i = 0; i < getNumBefore(motifLength+1) - getNumBefore(motifLength); i++)
		{
			result.push_back(getReverseComplement(getNumBefore(motifLength) + i));
		}
	}
	assert(result.size() == getNumBefore(MaxMotifLength+1));
	assert(result.size() < std::numeric_limits<CharType>::max());
	for (size_t i = 0; i < result.size(); i++)
	{
		assert(result[result[i]] == i);
	}
	return result;
}

std::vector<CharType> multiRLEReverseComplements = getReverseComplements();

CharType complement(CharType code)
{
	assert(code < multiRLEReverseComplements.size());
	return multiRLEReverseComplements[code];
}

size_t maxCode()
{
	return multiRLEReverseComplements.size();
}

std::pair<CharType, LengthType> getCodeAndRunlength(const SequenceCharType& str, size_t start, size_t end, uint16_t motifLength)
{
	assert(end > start);
	assert(end >= start + motifLength);
	uint16_t overhang = (end - start) % motifLength;
	uint16_t motif = 0;
	for (size_t i = start; i < start+motifLength; i++)
	{
		uint16_t mask = 0;
		assert(str[i] >= 0 && str[i] <= 3);
		mask = str[i];
		motif <<= 2;
		motif |= mask;
	}
	CharType code = getNumBefore(motifLength) + motif * motifLength + overhang;
	assert((end - start) / motifLength < std::numeric_limits<LengthType>::max()-1);
	LengthType runLength = ((end - start) - overhang) / motifLength;
	assert(runLength > 0);
	assert((((end - start) - overhang) % motifLength) == 0);
	return std::make_pair(code, runLength);
}

std::string getSequence(CharType code, size_t runLength)
{
	uint16_t motifLength = 1;
	while (getNumBefore(motifLength+1) <= code)
	{
		motifLength += 1;
		assert(motifLength <= MaxMotifLength);
	}
	code -= getNumBefore(motifLength);
	uint16_t motif = code / motifLength;
	uint16_t overhang = code % motifLength;
	std::string motifStr;
	for (size_t i = 0; i < motifLength; i++)
	{
		uint16_t charCode = (motif >> ((motifLength - i - 1) * 2)) & 3;
		switch(charCode)
		{
			case 0:
				motifStr += 'A';
				break;
			case 1:
				motifStr += 'C';
				break;
			case 2:
				motifStr += 'G';
				break;
			case 3:
				motifStr += 'T';
				break;
		}
	}
	std::string result;
	for (size_t i = 0; i < runLength; i++)
	{
		result += motifStr;
	}
	if (overhang > 0) result += motifStr.substr(0, overhang);
	return result;
}

template <typename F>
void iterateRuns(const SequenceCharType& str, const SequenceLengthType& poses, const size_t maxMaskLength, F callback)
{
	size_t lastRunEnd = 0;
	uint64_t runChecker = 0;
	assert(str.size() >= 32);
	for (size_t i = 0; i < 31; i++)
	{
		runChecker >>= 2;
		assert(str[i] >= 0 && str[i] <= 3);
		runChecker += ((uint64_t)str[i]) << 62LL;
	}
	for (size_t i = 0; i < str.size(); i++)
	{
		runChecker >>= 2;
		if (i+31 < str.size())
		{
			assert(str[i+31] >= 0 && str[i+31] <= 3);
			runChecker += ((uint64_t)str[i+31]) << 62LL;
		}
		if (i + 7 < lastRunEnd) continue;
		std::tuple<size_t, size_t, uint8_t> currentBestRun = std::make_tuple(i, i+1, 1);
		for (size_t motifLength = 2; motifLength <= maxMaskLength; motifLength++)
		{
			if (((runChecker ^ (runChecker >> (motifLength*2LL))) & ((1LL << (motifLength*2LL)) - 1LL)) != 0LL) continue;
			if (i + motifLength * 2 > str.size()) break;
			size_t runLength = 2;
			size_t overhang = 0;
			while (i+motifLength*runLength+motifLength <= str.size())
			{
				bool match = true;
				size_t j = 0;
				for (; j < motifLength; j++)
				{
					if (str[i+j] != str[i + motifLength * runLength + j])
					{
						match = false;
						break;
					}
				}
				if (match)
				{
					runLength += 1;
				}
				else
				{
					overhang = j;
					break;
				}
			}
			assert(runLength >= 2);
			if (i+motifLength*runLength+motifLength > str.size())
			{
				overhang = 0;
				for (size_t j = 0; i+motifLength*runLength+j < str.size(); j++)
				{
					if (str[i+j] != str[i+motifLength*runLength+j]) break;
					overhang += 1;
				}
			}
			else if (overhang == motifLength) overhang = 0;
			assert(runLength >= 2);
			assert(overhang < motifLength);
			size_t lengthHere = motifLength*runLength+overhang;
			if (i + lengthHere <= lastRunEnd)
			{
				continue;
			}
			if (i + lengthHere > std::get<1>(currentBestRun))
			{
				currentBestRun = std::make_tuple(i, i + lengthHere, motifLength);
			}
		}
		if (std::get<1>(currentBestRun) > lastRunEnd)
		{
			callback(currentBestRun);
			lastRunEnd = std::get<1>(currentBestRun);
		}
	}
}

template <typename F>
void iterateNonOverlappingRuns(const SequenceCharType& str, const SequenceLengthType& poses, const size_t maxMaskLength, F callback)
{
	std::tuple<size_t, size_t, uint8_t> lastRun { 0, 0, 0 };
	size_t lastOneChar = 0;
	bool first = true;
	iterateRuns(str, poses, maxMaskLength, [&callback, &lastRun, &first, &lastOneChar](const std::tuple<size_t, size_t, uint8_t> currentRun)
	{
		if (first)
		{
			lastRun = currentRun;
			first = false;
			return;
		}
		assert(std::get<0>(currentRun) > std::get<0>(lastRun));
		assert(std::get<1>(currentRun) > std::get<1>(lastRun));
		if (lastOneChar > std::get<0>(currentRun))
		{
			for (size_t j = lastOneChar; j < std::get<1>(lastRun); j++)
			{
				// maybe todo: hpc in overlapping parts
				callback(std::make_tuple(j, j+1, 1));
			}
			lastOneChar = std::get<1>(lastRun);
			lastRun = currentRun;
			return;
		}
		assert(lastOneChar <= std::get<0>(currentRun));
		if (std::get<0>(currentRun) >= std::get<1>(lastRun))
		{
			assert(std::get<0>(currentRun) == std::get<1>(lastRun));
			assert(std::max(lastOneChar, std::get<0>(lastRun)) < std::get<1>(lastRun));
			callback(std::make_tuple(std::max(lastOneChar, std::get<0>(lastRun)), std::get<1>(lastRun), std::get<2>(lastRun)));
			lastRun = currentRun;
			return;
		}
		assert(std::get<0>(currentRun) < std::get<1>(lastRun));
		assert(std::max(lastOneChar, std::get<0>(lastRun)) <= std::get<0>(currentRun));
		if (std::max(lastOneChar, std::get<0>(lastRun)) < std::get<0>(currentRun)) callback(std::make_tuple(std::max(lastOneChar, std::get<0>(lastRun)), std::get<0>(currentRun), std::get<2>(lastRun)));
		lastOneChar = std::get<1>(lastRun);
		for (size_t j = std::get<0>(currentRun); j < std::get<1>(lastRun); j++)
		{
			// maybe todo: hpc in overlapping parts
			callback(std::make_tuple(j, j+1, 1));
		}
		lastRun = currentRun;
	});
	if (lastOneChar > std::get<0>(lastRun))
	{
		if (lastOneChar != std::get<1>(lastRun)) callback(std::make_tuple(lastOneChar, std::get<1>(lastRun), std::get<2>(lastRun)));
	}
	else
	{
		callback(lastRun);
	}
}

std::pair<SequenceCharType, SequenceLengthType> multiRLECompressOne(const SequenceCharType& str, const SequenceLengthType& poses, const size_t maxMaskLength)
{
	assert(maxMaskLength <= MaxMotifLength);
	std::pair<SequenceCharType, SequenceLengthType> result;
	result.first.reserve(str.size());
	result.second.reserve(str.size());
	std::tuple<size_t, size_t, uint8_t> lastRun { 0, 0, 0 };
	iterateNonOverlappingRuns(str, poses, maxMaskLength, [&result, &lastRun, &str, &poses](const std::tuple<size_t, size_t, uint8_t> run)
	{
		assert(std::get<1>(lastRun) == std::get<0>(run));
		assert(std::get<1>(run) > std::get<0>(run));
		CharType code;
		LengthType runLength;
		uint8_t motifLength = std::get<2>(run);
		if (motifLength > std::get<1>(run) - std::get<0>(run)) motifLength = std::get<1>(run) - std::get<0>(run);
		std::tie(code, runLength) = getCodeAndRunlength(str, std::get<0>(run), std::get<1>(run), motifLength);
		result.first.emplace_back(code);
		result.second.emplace_back(poses[std::get<0>(run)]);
		lastRun = run;
	});
	result.second.emplace_back(poses[std::get<1>(lastRun)]);
	return result;
}

std::pair<SequenceCharType, SequenceLengthType> multiRLECompress(const SequenceCharType& str, const SequenceLengthType& poses, const size_t maxMaskLength)
{
	assert(maxMaskLength <= MaxMotifLength);
	assert(str.size() >= 32);
	std::vector<std::pair<SequenceCharType, SequenceLengthType>> result;
	size_t lastBreak = 0;
	for (size_t i = 0; i < str.size(); i++)
	{
		assert(str[i] <= 3);
	}
	return multiRLECompressOne(SequenceCharType { str.begin() + lastBreak, str.end() }, poses, maxMaskLength);
}

SequenceCharType revCompRLE(const SequenceCharType& codes)
{
	SequenceCharType result;
	result.resize(codes.size(), 0);
	for (size_t i = 0; i < result.size(); i++)
	{
		result[i] = complement(codes[codes.size()-1-i]);
	}
	return result;
}

}
