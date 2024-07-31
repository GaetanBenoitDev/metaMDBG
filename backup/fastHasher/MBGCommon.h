#ifndef MBGCommon_h
#define MBGCommon_h

#include <fstream>
#include <tuple>
#include <vector>
#include "VectorView.h"
#include "CompressedSequence.h"

namespace MBG
{

using HashType = unsigned __int128;
using NodeType = size_t;
using CharType = uint16_t;
using LengthType = size_t;
using SequenceCharType = std::vector<CharType>;
using SequenceLengthType = std::vector<LengthType>;
using CompressedSequenceType = CompressedSequence;
using ReadName = std::pair<std::string, size_t>;

HashType hash(VectorView<uint16_t> sequence);
HashType hash(VectorView<uint16_t> sequence, VectorView<uint16_t> reverseSequence);
HashType hash(std::vector<uint16_t> sequence);
std::ostream& operator<<(std::ostream& os, HashType t);
std::istream& operator>>(std::istream& is, HashType& t);
std::pair<size_t, bool> reverse(std::pair<size_t, bool> pos);
std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>> canon(std::pair<size_t, bool> from, std::pair<size_t, bool> to);
std::string revCompRaw(const std::string& seq);
std::vector<std::pair<size_t, bool>> revCompPath(const std::vector<std::pair<size_t, bool>>& original);

class PalindromicKmer : std::exception {};

}

namespace std
{
	template <> struct hash<const std::vector<size_t>&>
	{
		size_t operator()(const std::vector<size_t>& x) const
		{
			size_t h = 0;
			for (auto i : x)
			{
				h = h ^ std::hash<size_t>{}(i);
				h = (h << 3) + (h >> 61);
			}
			return h;
		}
	};
/*
#if !defined(__clang__)
	template <> struct hash<MBG::HashType>
	{
		size_t operator()(MBG::HashType x) const
		{
			return (size_t)x ^ (size_t)(x >> 64);
		}
	};
#endif*/
	template <> struct hash<std::pair<std::string, size_t>>
	{
		size_t operator()(const std::pair<std::string, size_t>& x) const
		{
			return hash<std::string>{}(x.first) ^ hash<size_t>{}(x.second);
		}
	};
	template <> struct hash<std::pair<MBG::HashType, bool>>
	{
		size_t operator()(std::pair<MBG::HashType, bool> x) const
		{
			return hash<MBG::HashType>{}(x.first);
		}
	};
	template <> struct hash<std::pair<MBG::HashType, MBG::HashType>>
	{
		size_t operator()(std::pair<MBG::HashType, MBG::HashType> x) const
		{
			return (size_t)x.first ^ (size_t)x.second;
		}
	};
	template <> struct hash<std::pair<size_t, bool>>
	{
		size_t operator()(std::pair<size_t, bool> x) const
		{
			return (size_t)x.first;
		}
	};
	template <> struct hash<std::pair<size_t, size_t>>
	{
		size_t operator()(std::pair<size_t, size_t> x) const
		{
			return (size_t)x.first ^ (size_t)x.second;
		}
	};
	template <> struct hash<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
	{
		size_t operator()(std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>> x) const
		{
			return (hash<std::pair<size_t, bool>>{}(x.first) << 16) + hash<std::pair<size_t, bool>>{}(x.second);
		}
	};
}
#endif
