

#ifndef MDBG_METAG_COMMONS
#define MDBG_METAG_COMMONS

//#include <gatb/gatb_core.hpp>
#include "utils/MurmurHash3.h"
#include "zlib.h"
#include <string>
#include <sstream>
#include <unordered_set>
#include <unordered_map>
#include <regex>
#include <algorithm>
#include <libgen.h>
#include <set>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <cmath>
#include "./utils/ArgParse.hpp"
#include <random>
#include <chrono>
namespace fs = std::filesystem;
//namespace fs = std::filesystem;

//#include <sys/types.h>
//#include <sys/stat.h>

#include <functional>
#include <memory>
#include "./utils/ntHashIterator.hpp"
#include "./utils/kseq.h"
KSEQ_INIT(gzFile, gzread)

using namespace std;

/*
#define STR_OUTPUT "-o"
#define STR_INPUT "-i"
#define STR_MINIM_SIZE "-l"
#define STR_KMINMER_SIZE "-k"
#define STR_DENSITY "-d"
#define STR_INPUT_DIR "-idir"
#define STR_INPUT_TRUTH "-itruth"
#define STR_HIFIASM_DEBUG "--debug"
*/

//typedef uint64_t u_int64_t;
//typedef uint32_t u_int64_t;
//typedef uint16_t u_int64_t;
//typedef uint8_t u_int8_t;

//KMER_SPAN(1)
//typedef Kmer<>::ModelDirect    ModelDirect;
//typedef Kmer<>::ModelCanonical ModelCanonical;
//typedef Kmer<>::Type  kmer_type;
//typedef Kmer<>::Count kmer_count;
//typedef typename Kmer<>::Type  Type;
//typedef typename Kmer<>::Count Count;
//typedef Kmer<>::ModelMinimizer<ModelCanonical> ModelMinimizer;

typedef u_int32_t ReadIndexType;
//complement of one NT
const unsigned char comp_NT[4] = {  2,3,0,1  };
const char bin2NT[4] = {'A','C','T','G'};
const char binrev[4] = {2,3,0,1};

//reverse complement of 4NT,  ie one byte
const unsigned char revcomp_4NT[256] = {
 0xaa,
 0xea,
 0x2a,
 0x6a,
 0xba,
 0xfa,
 0x3a,
 0x7a,
 0x8a,
 0xca,
 0xa,
 0x4a,
 0x9a,
 0xda,
 0x1a,
 0x5a,
 0xae,
 0xee,
 0x2e,
 0x6e,
 0xbe,
 0xfe,
 0x3e,
 0x7e,
 0x8e,
 0xce,
 0xe,
 0x4e,
 0x9e,
 0xde,
 0x1e,
 0x5e,
 0xa2,
 0xe2,
 0x22,
 0x62,
 0xb2,
 0xf2,
 0x32,
 0x72,
 0x82,
 0xc2,
 0x2,
 0x42,
 0x92,
 0xd2,
 0x12,
 0x52,
 0xa6,
 0xe6,
 0x26,
 0x66,
 0xb6,
 0xf6,
 0x36,
 0x76,
 0x86,
 0xc6,
 0x6,
 0x46,
 0x96,
 0xd6,
 0x16,
 0x56,
 0xab,
 0xeb,
 0x2b,
 0x6b,
 0xbb,
 0xfb,
 0x3b,
 0x7b,
 0x8b,
 0xcb,
 0xb,
 0x4b,
 0x9b,
 0xdb,
 0x1b,
 0x5b,
 0xaf,
 0xef,
 0x2f,
 0x6f,
 0xbf,
 0xff,
 0x3f,
 0x7f,
 0x8f,
 0xcf,
 0xf,
 0x4f,
 0x9f,
 0xdf,
 0x1f,
 0x5f,
 0xa3,
 0xe3,
 0x23,
 0x63,
 0xb3,
 0xf3,
 0x33,
 0x73,
 0x83,
 0xc3,
 0x3,
 0x43,
 0x93,
 0xd3,
 0x13,
 0x53,
 0xa7,
 0xe7,
 0x27,
 0x67,
 0xb7,
 0xf7,
 0x37,
 0x77,
 0x87,
 0xc7,
 0x7,
 0x47,
 0x97,
 0xd7,
 0x17,
 0x57,
 0xa8,
 0xe8,
 0x28,
 0x68,
 0xb8,
 0xf8,
 0x38,
 0x78,
 0x88,
 0xc8,
 0x8,
 0x48,
 0x98,
 0xd8,
 0x18,
 0x58,
 0xac,
 0xec,
 0x2c,
 0x6c,
 0xbc,
 0xfc,
 0x3c,
 0x7c,
 0x8c,
 0xcc,
 0xc,
 0x4c,
 0x9c,
 0xdc,
 0x1c,
 0x5c,
 0xa0,
 0xe0,
 0x20,
 0x60,
 0xb0,
 0xf0,
 0x30,
 0x70,
 0x80,
 0xc0,
 0x0,
 0x40,
 0x90,
 0xd0,
 0x10,
 0x50,
 0xa4,
 0xe4,
 0x24,
 0x64,
 0xb4,
 0xf4,
 0x34,
 0x74,
 0x84,
 0xc4,
 0x4,
 0x44,
 0x94,
 0xd4,
 0x14,
 0x54,
 0xa9,
 0xe9,
 0x29,
 0x69,
 0xb9,
 0xf9,
 0x39,
 0x79,
 0x89,
 0xc9,
 0x9,
 0x49,
 0x99,
 0xd9,
 0x19,
 0x59,
 0xad,
 0xed,
 0x2d,
 0x6d,
 0xbd,
 0xfd,
 0x3d,
 0x7d,
 0x8d,
 0xcd,
 0xd,
 0x4d,
 0x9d,
 0xdd,
 0x1d,
 0x5d,
 0xa1,
 0xe1,
 0x21,
 0x61,
 0xb1,
 0xf1,
 0x31,
 0x71,
 0x81,
 0xc1,
 0x1,
 0x41,
 0x91,
 0xd1,
 0x11,
 0x51,
 0xa5,
 0xe5,
 0x25,
 0x65,
 0xb5,
 0xf5,
 0x35,
 0x75,
 0x85,
 0xc5,
 0x5,
 0x45,
 0x95,
 0xd5,
 0x15,
 0x55
};

const string ARG_INPUT_FILENAME = "i";
const string ARG_INPUT_FILENAME_TRUTH = "itruth";
const string ARG_OUTPUT_DIR = "o";
const string ARG_MINIMIZER_LENGTH = "l";
const string ARG_KMINMER_LENGTH = "k";
const string ARG_MINIMIZER_DENSITY = "d";
const string ARG_DEBUG = "debug";
const string ARG_INPUT_FILENAME_CONTIG = "c";
const string ARG_INPUT_FILENAME_UNITIG_NT = "unitigNt";
const string ARG_INPUT_FILENAME_UNITIG_CLUSTER = "cluster";
const string ARG_INPUT_FILENAME_ABUNDANCE = "a";

struct UnitigData{
	u_int32_t _index;
	vector<u_int64_t> _readIndexes;
	//float _meanAbundance;
	//u_int16_t _nbKminmers;
    //vector<float> _compositionMean;
    //u_int32_t _compositionNb;
	//unordered_set<ReadIndexType> _readIndexes_exists;
};

struct ReadKminmer{
	u_int32_t _read_pos_start;
	u_int32_t _read_pos_end;
	u_int16_t _length;
	bool _isReversed;
	u_int16_t _position_of_second_minimizer;
	u_int16_t _position_of_second_to_last_minimizer;
	u_int16_t _position_of_second_minimizer_seq;
	u_int16_t _position_of_second_to_last_minimizer_seq;
	u_int16_t _seq_length_start;
	u_int16_t _seq_length_end;
};


struct DbgEdge{
	u_int32_t _from;
	u_int32_t _to;

	bool operator==(const DbgEdge &other) const{
		return ((_from == other._from) && (_to == other._to));
	}

    size_t hash() const{
		std::size_t seed = 2;
		seed ^= _from + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		seed ^= _to + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        //auto hash1 = hash<T1>{}(p._from);
        //auto hash2 = hash<T2>{}(p._to);
        return seed;
	}

    DbgEdge normalize(){
        DbgEdge edge1 = {_from, _to};
        DbgEdge edge2 = {_to, _from};
        if(edge1.hash() < edge2.hash()){
            return edge1;
        }
        else{
            return edge2;
        }
    }

};


struct hash_pair {
    size_t operator()(const DbgEdge& p) const
    {
        return p.hash();
		//std::size_t seed = 2;
		//seed ^= p._from + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		//seed ^= p._to + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        //auto hash1 = hash<T1>{}(p._from);
        //auto hash2 = hash<T2>{}(p._to);
        //return seed;
    }
};


struct UnitigEdgeScore{
	u_int32_t _from;
	u_int32_t _to;
	u_int16_t _score;
};

struct MinimizerPair{
	u_int64_t _first;
	u_int64_t _second;

	bool operator==(const MinimizerPair &other) const{
		return _first == other._first && _second == other._second;
	}
};


struct DbgNode{
	u_int32_t _index;
	u_int32_t _abundance;
	u_int16_t _length;
	u_int16_t _overlapLength_start;
	u_int16_t _overlapLength_end;
	bool _isReversed;
};



struct ReadData{
	u_int32_t _length;
	vector<float> _composition;
};

struct KmerVec{
	vector<u_int64_t> _kmers;

	/*
	KmerVec clone(){
		KmerVec vec;
		vec._kmers = _kmers;
		return vec;
	}*/

	bool operator==(const KmerVec &other) const{
		return _kmers == other._kmers;
	}

	KmerVec prefix() const{
		KmerVec vec = (*this);
		vec._kmers.pop_back();
		//vec._kmers.erase(vec._kmers.begin());
		return vec;
	}

	KmerVec suffix() const{
		KmerVec vec = (*this);
		vec._kmers.erase(vec._kmers.begin());
		return vec;
	}

	KmerVec reverse() const{
		KmerVec vec_reverse = (*this);
		std::reverse(vec_reverse._kmers.begin(), vec_reverse._kmers.end());
		return vec_reverse;
	}

	KmerVec normalize(){

		KmerVec vec_reverse = reverse();

		if(h() < vec_reverse.h()){
			return *this;
		}
		else{
			return vec_reverse;
		}
		
	}

	KmerVec normalize(bool& isReversed){

		KmerVec vec_reverse = reverse();

		if(h() < vec_reverse.h()){
			isReversed = false;
			return *this;
		}
		else{
			isReversed = true;
			return vec_reverse;
		}
		
	}

	bool isPalindrome() const{
		//vector<u_int64_t> kmers = _kmers;
		return  equal(_kmers.begin(), _kmers.begin() + _kmers.size()/2, _kmers.rbegin()); //suffix() == prefix().reverse();
	}

	size_t h() const{
		std::size_t seed = _kmers.size();
		for(auto& i : _kmers) {
			seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		}
		return seed;
	}
	
	string toString(){
		string s = "";
		for(u_int64_t m : _kmers){
			s += to_string(m) + " ";
		}
		return s;
	}
};



  
struct ReadKminmerComplete{
	KmerVec _vec;
	bool _isReversed;
	u_int32_t _read_pos_start;
	u_int32_t _read_pos_end;
	u_int32_t _seq_length_start;
	u_int32_t _seq_length_end;
	u_int32_t _length;
};

struct ContigNode{
	u_int32_t _nodeIndex;
	u_int64_t _supportingReadIndex;

	bool operator==(const ContigNode &other) const{
		return _nodeIndex == other._nodeIndex && _supportingReadIndex == other._supportingReadIndex;
	}

};

struct MinimizerPair_Edge{
	MinimizerPair _from;
	MinimizerPair _to;
};

/*
namespace std {
	template <>
	struct hash<DbgEdge>{
		std::size_t operator()(const DbgEdge& edge) const{
			std::size_t seed = 2;
			//for(auto& i : _kmers) {
			seed ^= edge._from + 0x9e3779b9 + (seed << 6) + (seed >> 2);
			seed ^= edge._to + 0x9e3779b9 + (seed << 6) + (seed >> 2);
			//}
			return seed;
		}
	};
}*/

namespace std {
	template <>
	struct hash<MinimizerPair>{
		std::size_t operator()(const MinimizerPair& k) const{
			using std::size_t;
			using std::hash;
			using std::string;

			return ((hash<u_int64_t>()(k._first) ^ (hash<u_int64_t>()(k._second) << 1)) >> 1);
		}
	};


	template <>
	struct hash<ContigNode>{
		std::size_t operator()(const ContigNode& k) const{
			//using std::size_t;
			//using std::hash;
			//using std::string;

			std::size_t seed = 2;
			seed ^= k._nodeIndex + 0x9e3779b9 + (seed << 6) + (seed >> 2);
			seed ^= k._supportingReadIndex + 0x9e3779b9 + (seed << 6) + (seed >> 2);
			//auto hash1 = hash<T1>{}(p._from);
			//auto hash2 = hash<T2>{}(p._to);
			return seed;

			//return ((hash<u_int64_t>()(k._first) ^ (hash<u_int64_t>()(k._second) << 1)) >> 1);
		}
	};
}

namespace std {
	template <>
	struct hash<KmerVec>{
		std::size_t operator()(const KmerVec& vec) const{
			return vec.h();
			/*
			std::size_t seed = vec._kmers.size();
			for(auto& i : vec._kmers) {
				seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
			}
			return seed;
			*/
		}
	};
}

/*
struct ReadData{
	u_int64_t _length;
	vector<float> _composition;
	vector<u_int64_t> _minimizers;
	//vertex_t _graphVertex;
	//bool _vertexCreated;
	//ReadIndexType _readBestMatch;
	//int _readBestMatch_dist;
	//vector<ReadIndexType> _overlaps;
};
*/
struct Overlap {
    ReadIndexType _r1;
    ReadIndexType _r2;
    u_int16_t _nbMinimizers;
};

static const unsigned char basemap[256] = {
	0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
	16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
	32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
	48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
	64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
	'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',  91,  92,  93,  94,  95,
	96, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
	'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127,
	128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
	144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
	160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
	176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
	192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
	208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223,
	224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
	240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255
};

class Utils{




public:

	static double compute_first_quartile(vector<u_int64_t> scores){
		size_t size = scores.size();

		if (size == 0){
			return 0;  // Undefined, really.
		}
		else{
			sort(scores.begin(), scores.end());
			return scores[size * 0.1];
		}
	}

	static double compute_median(vector<u_int32_t> scores){
		size_t size = scores.size();

		if (size == 0){
			return 0;  // Undefined, really.
		}
		else{
			sort(scores.begin(), scores.end());
			if (size % 2 == 0){
				return (scores[size / 2 - 1] + scores[size / 2]) / 2; //scores[size / 2 - 1];
			}
			else {
				return scores[size / 2];
			}
		}
	}

	static double compute_median_float(vector<float> scores){
		size_t size = scores.size();

		if (size == 0){
			return 0;  // Undefined, really.
		}
		else{
			sort(scores.begin(), scores.end());
			if (size % 2 == 0){
				return (scores[size / 2 - 1] + scores[size / 2]) / 2; //scores[size / 2 - 1];
			}
			else {
				return scores[size / 2];
			}
		}
	}

	static void revcomp(string& sequence){

		std::reverse(sequence.begin(), sequence.end());
		
		for(size_t i=0; i<sequence.size(); i++){
			sequence[i] = basemap[sequence[i]];
		}
	}

	static float computeJaccardDistance(const vector<u_int32_t>& v1, const vector<u_int32_t>& v2){

		size_t i=0;
		size_t j=0;
		u_int64_t nbShared = 0;
		u_int64_t nbElements = 0;


		while(i < v1.size() && j < v2.size()){


			if(v1[i] == v2[j]){
				nbShared += 1;
				i += 1;
				j += 1;
				nbElements += 1;
			}
			else if(v1[i] < v2[j]){
				i += 1;
				nbElements += 1;
			}
			else{
				j += 1;
				nbElements += 1;
			}

		}

		return 1 - ((float)nbShared) / nbElements;
	}

	static u_int64_t computeSharedElements(const vector<u_int32_t>& reads1, const vector<u_int32_t>& reads2){

		size_t i=0;
		size_t j=0;
		u_int64_t nbShared = 0;

		while(i < reads1.size() && j < reads2.size()){

			if(reads1[i] == reads2[j]){
				nbShared += 1;
				i += 1;
				j += 1;
			}
			else if(reads1[i] < reads2[j]){
				i += 1;
			}
			else{
				j += 1;
			}

		}

		return nbShared;
	}

	static u_int64_t collectSharedElements(const vector<u_int32_t>& reads1, const vector<u_int32_t>& reads2, unordered_set<u_int32_t>& sharedReads){

		sharedReads.clear();

		size_t i=0;
		size_t j=0;
		u_int64_t nbShared = 0;

		while(i < reads1.size() && j < reads2.size()){
			if(reads1[i] == reads2[j]){
				sharedReads.insert(reads1[i]);
				nbShared += 1;
				i += 1;
				j += 1;
			}
			else if(reads1[i] < reads2[j]){
				i += 1;
			}
			else{
				j += 1;
			}

		}

		return nbShared;
	}

	static u_int64_t computeSharedReads(const vector<u_int64_t>& reads1, const vector<u_int64_t>& reads2){

		size_t i=0;
		size_t j=0;
		u_int64_t nbShared = 0;

		while(i < reads1.size() && j < reads2.size()){

			if(reads1[i] == reads2[j]){
				nbShared += 1;
				i += 1;
				j += 1;
			}
			else if(reads1[i] < reads2[j]){
				i += 1;
			}
			else{
				j += 1;
			}

		}

		return nbShared;
	}

	static u_int64_t computeSharedReads(const UnitigData& utg1, const UnitigData& utg2){

		//if(utg1._readIndexes.size() == 0 || utg2._readIndexes.size() == 0) return 1;
		//cout << "------------------- " << utg1._index << endl;
		//for(size_t i=0; i<utg1._readIndexes.size(); i++){
		//	cout << "| " << utg1._readIndexes[i] << endl;
		//}
		//cout << "- " << utg2._index << endl;
		//for(size_t i=0; i<utg2._readIndexes.size(); i++){
		//	cout << "| " << utg2._readIndexes[i] << endl;
		//}

		size_t i=0;
		size_t j=0;
		u_int64_t nbShared = 0;

		while(i < utg1._readIndexes.size() && j < utg2._readIndexes.size()){
			//if(removedReadIndex.find(utg1._readIndexes[i]) != removedReadIndex.end()){
			//	i += 1;
			//	continue;
			//}
			//if(removedReadIndex.find(utg2._readIndexes[j]) != removedReadIndex.end()){
			//	j += 1;
			//	continue;
			//}

			if(utg1._readIndexes[i] == utg2._readIndexes[j]){
				nbShared += 1;
				i += 1;
				j += 1;
			}
			else if(utg1._readIndexes[i] < utg2._readIndexes[j]){
				i += 1;
			}
			else{
				j += 1;
			}

		}

		return nbShared;
	}

	static u_int64_t collectSharedReads(const vector<u_int64_t>& reads1, const vector<u_int64_t>& reads2, vector<u_int64_t>& sharedReads){

		sharedReads.clear();

		size_t i=0;
		size_t j=0;
		u_int64_t nbShared = 0;

		while(i < reads1.size() && j < reads2.size()){
			if(reads1[i] == reads2[j]){
				sharedReads.push_back(reads1[i]);
				nbShared += 1;
				i += 1;
				j += 1;
			}
			else if(reads1[i] < reads2[j]){
				i += 1;
			}
			else{
				j += 1;
			}

		}

		return nbShared;
	}

	static u_int64_t collectSharedReads(const UnitigData& utg1, const UnitigData& utg2, vector<u_int64_t>& sharedReads){

		sharedReads.clear();

		//if(utg1._readIndexes.size() == 0 || utg2._readIndexes.size() == 0) return;
		//cout << "------------------- " << utg1._index << endl;
		//for(size_t i=0; i<utg1._readIndexes.size(); i++){
		//	cout << "| " << utg1._readIndexes[i] << endl;
		//}
		//cout << "- " << utg2._index << endl;
		//for(size_t i=0; i<utg2._readIndexes.size(); i++){
		//	cout << "| " << utg2._readIndexes[i] << endl;
		//}


		size_t i=0;
		size_t j=0;
		u_int64_t nbShared = 0;

		while(i < utg1._readIndexes.size() && j < utg2._readIndexes.size()){
			if(utg1._readIndexes[i] == utg2._readIndexes[j]){
				sharedReads.push_back(utg1._readIndexes[i]);
				nbShared += 1;
				i += 1;
				j += 1;
			}
			else if(utg1._readIndexes[i] < utg2._readIndexes[j]){
				i += 1;
			}
			else{
				j += 1;
			}

		}

		return nbShared;
	}
	

	static bool shareAnyRead(const UnitigData& utg1, const UnitigData& utg2){

		if(utg1._readIndexes.size() == 0 || utg2._readIndexes.size() == 0) return true;
		//cout << "------------------- " << utg1._index << endl;
		//for(size_t i=0; i<utg1._readIndexes.size(); i++){
		//	cout << "| " << utg1._readIndexes[i] << endl;
		//}
		//cout << "- " << utg2._index << endl;
		//for(size_t i=0; i<utg2._readIndexes.size(); i++){
		//	cout << "| " << utg2._readIndexes[i] << endl;
		//}

		size_t i=0;
		size_t j=0;

		while(i < utg1._readIndexes.size() && j < utg2._readIndexes.size()){

			//cout << i << " " << j << endl;
			if(utg1._readIndexes[i] == utg2._readIndexes[j]){
				return true;
			}
			else if(utg1._readIndexes[i] < utg2._readIndexes[j]){
				i += 1;
			}
			else{
				j += 1;
			}

		}

		return false;
	}

};


class Encoder{

public: 

	static void encode_rle(const char* sequence, size_t length, string& rleSequence, vector<u_int64_t>& rlePositions) {
		/*
		string sequence_str;

				char lastChar = '0';
				for(size_t i=0; i<sequence.getDataSize(); i++){
					if(readseq[i] == lastChar) continue;
					sequence_str += readseq[i];
					lastChar = readseq[i];
				}
		*/
		rleSequence = "";
		char lastChar = '#';
		u_int64_t lastPos = 0;

		for(size_t i=0; i<length; i++){
			//cout << i << " " << length << endl;
			char c = sequence[i];
			if(c == lastChar) continue;
			if(lastChar != '#'){
				rleSequence += lastChar;
				rlePositions.push_back(lastPos);
				lastPos = i;
				//cout << lastChar << endl;
			}
			lastChar = c;
		}
		//cout << lastChar << endl;
		rleSequence += lastChar;
		rlePositions.push_back(lastPos);
		rlePositions.push_back(length);
	}

};

class MDBG{

public:

	size_t _k;
	u_int32_t _node_id;
	unordered_map<KmerVec, DbgNode> _dbg_nodes;
	unordered_map<KmerVec, vector<KmerVec>> _dbg_edges;

	MDBG(size_t k){
		_k = k;
		_node_id = 0;
	}

	void addNode(const KmerVec& vec, u_int16_t length, u_int16_t overlapLength_start, u_int16_t overlapLength_end, bool isReversed){




		if(_dbg_nodes.find(vec) != _dbg_nodes.end()){
			//if(_dbg_nodes[vec]._abundance > 1000){
			//	cout << _dbg_nodes[vec]._index << " " << _dbg_nodes[vec]._abundance << endl;
			//}
			_dbg_nodes[vec]._abundance += 1;
			return;
		}


		DbgNode node = {_node_id, 1, length, overlapLength_start, overlapLength_end, isReversed};
		_dbg_nodes[vec] = node;

		_dbg_edges[vec.prefix().normalize()].push_back(vec);
		_dbg_edges[vec.suffix().normalize()].push_back(vec);
		
		//cout << _dbg_edges.size() << endl;

		//if(_node_id == 2285){
		//	cout << vec.isPalindrome() << endl;
		//	cout << vec._kmers[0] << endl;
		//	cout << vec._kmers[1] << endl;
		//	cout << vec._kmers[2] << endl;
		//}

		_node_id += 1;

	}
	
	void dump(const string& filename){
		gzFile file = gzopen(filename.c_str(),"wb");

		//bool lala = true;
		for(auto it : _dbg_nodes){
			const KmerVec& vec = it.first;
			const DbgNode& node = it.second;

			gzwrite(file, (const char*)&vec._kmers[0], _k * sizeof(u_int64_t));
			//cout << sizeof(DbgNode) << endl;
			gzwrite(file, (const char*)&node, sizeof(DbgNode));

			/*
			if(lala){
				lala = false;
				for(size_t i=0; i<vec._kmers.size(); i++){
					cout << vec._kmers[i] << endl;
				}
				cout << node._index << endl;
				cout << node._abundance << endl;
				cout << node._length << endl;
			}*/

		}

		gzclose(file);
	}

	void load(const string& filename){

		//cout << _k << endl;
		gzFile file = gzopen(filename.c_str(),"rb");

		while(true){

			KmerVec vec;
			vec._kmers.resize(_k);
			DbgNode node;

			//cout << vec._kmers.size() << endl;
			gzread(file, (char*)&vec._kmers[0], _k * sizeof(u_int64_t));




			if(gzeof(file)) break;
			

			gzread(file, (char*)&node, sizeof(node));

			//cout << node._index << endl;
			_dbg_nodes[vec] = node;
			/*
			if(_node_id == 0){
				for(size_t i=0; i<vec._kmers.size(); i++){
					cout << vec._kmers[i] << endl;
				}
				cout << node._index << endl;
				cout << node._abundance << endl;
				cout << node._length << endl;
			}*/

			_node_id += 1;
		}
		/*
		for(auto it : _dbg_nodes){
			const KmerVec& vec = it.first;
			const DbgNode& node = it.second;

			gzwrite(file, (const char*)&vec._kmers[0], _k * sizeof(u_int64_t));
			//cout << sizeof(DbgNode) << endl;
			gzwrite(file, (const char*)&node, sizeof(DbgNode));
		}
		*/

		gzclose(file);
	}



	static void getKminmers(const size_t l, const size_t k, const vector<u_int64_t>& minimizers, const vector<u_int64_t>& minimizersPos, vector<KmerVec>& kminmers, vector<ReadKminmer>& kminmersLength, const vector<u_int64_t>& rlePositions, int readIndex, bool allowPalindrome){

		/*
		bool isLala = false;
		if(std::find(minimizers.begin(), minimizers.end(), 56801217747741349) !=  minimizers.end()){
			isLala = true;
		}*/
		
        kminmersLength.clear();
        bool doesComputeLength = minimizersPos.size() > 0;
		/*
		vector<u_int64_t> minimizers_filtered;

		for(u_int64_t minimizer : minimizers){
			if(filteredMinimizers.find(minimizer) != filteredMinimizers.end()) continue;

			minimizers_filtered.push_back(minimizer);
		}
		*/

		/*
		if(minimizers.size() < _kminmerSize) return;

		int i_max = ((int)minimizers.size()) - (int)_kminmerSize + 1;
		for(int i=0; i<i_max; i++){

			KmerVec vec;

			bool valid = true;
			for(int j=i; j<i+_kminmerSize; j++){
				u_int64_t minimizer = minimizers[j];
				if(filteredMinimizers.find(minimizer) != filteredMinimizers.end()){
					valid = false;
					break;
				}

				vec._kmers.push_back(minimizer);
			}

			if(valid) kminmers.push_back(vec.normalize());
			//mdbg_repeatFree->addNode(vec.normalize());
			//kminmers.push_back(vec.normalize());
		}*/
		
		/*
		if(readIndex == 39679){
			for(size_t i=0; i<minimizers.size(); i++){
				cout << i << ": " << minimizers[i] << endl;
			}
		}
		*/
		
		if(minimizers.size() < k) return;

		//unordered_set<u_int64_t> bannedMinimizers;
		vector<bool> bannedPositions(minimizers.size(), false);






		/*
		size_t banned_k = 3;

		while(true){

			//if(readIndex == 39679){
			//	cout << "------------------------" << endl;
			//}

			bool hasPalindrome = false;
			KmerVec prevVec;

			int i_max = ((int)minimizers.size()) - (int)banned_k + 1;
			for(int i=0; i<i_max; i++){

				//if(readIndex == 94931){
				//	cout << i << ": " << minimizers[i] << endl;
				//}

				if(bannedPositions[i]) continue;

				KmerVec vec;
				vector<u_int32_t> currentMinimizerIndex;

				int j=i;
				while(true){
					
					if(j >= minimizers.size()) break;
					if(bannedPositions[j]){
						j += 1;
						continue;
					}

					u_int64_t minimizer = minimizers[j];


					vec._kmers.push_back(minimizer);
					currentMinimizerIndex.push_back(j);

					if(vec._kmers.size() == banned_k){

						
						if((vec.isPalindrome() || (i > 0 && vec.normalize() == prevVec.normalize()))){ //Palindrome: 121 (créé un cycle), Large palindrome = 122 221 (créé une tip)

							for(size_t m=0; m<banned_k; m++){
								bannedPositions[currentMinimizerIndex[m]] = true;
								//cout << "Banned: " << currentMinimizerIndex[m] << endl;
							}
							
							
							for(size_t i=0; i<bannedPositions.size()-banned_k+1; i++){
								if(bannedPositions[i]) continue;
								bool isBanned = false;
								for(size_t j=i; j<i+banned_k; j++){
									if(bannedPositions[j]){
										isBanned = true;
									}
								}

								if(isBanned){
									for(size_t j=i; j<i+banned_k; j++){
										cout << "Banned: " << j << endl;
										bannedPositions[j] = true;
									}
								}
							}
							
							

							hasPalindrome = true;
							break;
						}
					}

					j += 1;
				}

				if(hasPalindrome) break;
			}

			
			if(!hasPalindrome) break;
			//kminmers.clear();
        	//kminmersLength.clear();
		}
		*/
		/*
		if(isLala){
			for(size_t i=0; i<minimizers.size(); i++){
				cout << i << " " << minimizers[i] << "     " << bannedPositions[i] << endl;
				if(minimizers[i] == 56801217747741349){
					cout << "\tb" << endl; 
				}
				//if(bannedPositions[i]){
				//	cout << "Banned: " << i << endl;
				//}
			}
			getchar();
		}*/


		while(true){

			//if(readIndex == 39679){
			//	cout << "------------------------" << endl;
			//}

			bool hasPalindrome = false;
			KmerVec prevVec;

			int i_max = ((int)minimizers.size()) - (int)k + 1;
			for(int i=0; i<i_max; i++){

				//if(readIndex == 94931){
				//	cout << i << ": " << minimizers[i] << endl;
				//}

				if(bannedPositions[i]){
					//cout << "ban i" << endl;
					continue;
				}
				KmerVec vec;
				vector<u_int32_t> currentMinimizerIndex;

				int j=i;
				while(true){
					
					if(j >= minimizers.size()) break;
					if(bannedPositions[j]){
						//cout << "ban j" << " " << j << endl;
						j += 1;
						continue;
					}

					u_int64_t minimizer = minimizers[j];
					//cout << "MU: " << j << " " << minimizer << endl;

					//if(bannedMinimizers.find(minimizer) != bannedMinimizers.end()) continue;

					vec._kmers.push_back(minimizer);
					currentMinimizerIndex.push_back(j);

					if(vec._kmers.size() == k){

						/*
						if(readIndex == 39679){
							cout << "-----" << endl;
							cout << vec._kmers[0] << endl;
							cout << vec._kmers[1] << endl;
							cout << vec._kmers[2] << endl;
						}
						*/

						
						if(vec.isPalindrome() || (i > 0 && vec.normalize() == prevVec.normalize())){ //Palindrome: 121 (créé un cycle), Large palindrome = 122 221 (créé une tip)

							//if(readIndex == 96573){
							//	cout << "\tPalouf" << endl;
							//}

							//cout << "Palindrome!" << endl;
							//for(size_t p=j-_kminmerSize+1; p<=j-1; p++){
							//	bannedPositions[p] = true;
							//}
							for(size_t m=0; m<k; m++){
								bannedPositions[currentMinimizerIndex[m]] = true;

								//if(readIndex == 39679){
								//	cout << "Banned: " << currentMinimizerIndex[m] << endl;
								//}
								//cout << "Banned: " << currentMinimizerIndex[m] << endl;
							}
							/*
							for(size_t i=0; i<bannedPositions.size()-k+1; i++){
								if(bannedPositions[i]) continue;
								bool isBanned = false;
								for(size_t j=i; j<i+k; j++){
									if(bannedPositions[j]){
										isBanned = true;
									}
								}

								if(isBanned){
									for(size_t j=i; j<i+k; j++){
										bannedPositions[j] = true;
									}
								}
							}*/
							

							hasPalindrome = true;
							break;
						}
						else{

							bool isReversed;
							vec = vec.normalize(isReversed);

							//cout << "Is reversed: " << isReversed << endl;
                            if(doesComputeLength){

								u_int32_t indexFirstMinimizer = currentMinimizerIndex[0];
								u_int32_t indexSecondMinimizer = currentMinimizerIndex[1];
								u_int32_t indexSecondLastMinimizer = currentMinimizerIndex[currentMinimizerIndex.size()-2];
								u_int32_t indexLastMinimizer = currentMinimizerIndex[currentMinimizerIndex.size()-1];


								u_int32_t read_pos_start = minimizersPos[indexFirstMinimizer];
								u_int32_t read_pos_end = minimizersPos[indexLastMinimizer];
								if(rlePositions.size() > 0){
									read_pos_start = rlePositions[minimizersPos[indexFirstMinimizer]];
									read_pos_end = rlePositions[minimizersPos[indexLastMinimizer]];
									read_pos_end +=  (rlePositions[minimizersPos[indexLastMinimizer] + l] - rlePositions[minimizersPos[indexLastMinimizer]]); //l-1 a check
								}
								else{
									read_pos_end +=  l;//(minimizersPos[indexLastMinimizer] + l - minimizersPos[indexLastMinimizer]); //l-1 a check

								}
								//u_int32_t read_pos_start = rlePositions[minimizersPos[indexFirstMinimizer]];
								//u_int32_t read_pos_end = rlePositions[minimizersPos[indexLastMinimizer]]; // + (rlePositions[minimizersPos[i+k-1]] - rlePositions[minimizersPos[i+k-1] + l]); //+ l;
								//read_pos_end +=  (rlePositions[minimizersPos[indexLastMinimizer] + l] - rlePositions[minimizersPos[indexLastMinimizer]]); //l-1 a check

								//cout << "HI: " << rlePositions[minimizersPos[i+k-1] + l - 1] << endl;
								//cout << "HI: " << rlePositions[minimizersPos[i+k-1] + l] << endl;
								u_int16_t length = read_pos_end - read_pos_start;
								//cout << read_pos_start << " " << read_pos_end << " " << length << endl;
								/*
								if(readIndex == 159){
												cout << "----------------" << endl;
												cout << vec._kmers[0] << endl;
												cout << vec._kmers[1] << endl;
												cout << vec._kmers[2] << endl;
												//cout << read_pos_start << " " << read_pos_end << endl;
												//cout << currentMinimizerIndex[0] << endl;
												//cout << currentMinimizerIndex[1] << endl;
												//cout << currentMinimizerIndex[2] << endl;

												cout << read_pos_start << " " << read_pos_end << endl;
												//cout << "huuu" << endl;
								}*/
								
								// seqSize = read_pos_end - read_pos_start;

								//if(isReversed){
								//	read_pos_start = rlePositions[minimizersPos[i+k-1]];
								//	read_pos_end = rlePositions[minimizersPos[i]];
								//}

								/*
								cout << "------------" << endl;
								cout << rlePositions[minimizersPos[i]] << endl;
								cout << rlePositions[minimizersPos[i+1]] << endl;
								cout << rlePositions[minimizersPos[i+2]] << endl;
								cout << "Length: " << length << endl;
								*/

								u_int16_t seq_length_start = 0;
								u_int16_t seq_length_end = 0;

								u_int16_t position_of_second_minimizer = 0;
								u_int16_t position_of_second_minimizer_seq = 0;
								if(isReversed){
									//position_of_second_minimizer = rlePositions[minimizersPos[indexSecondLastMinimizer]]; //rlePositions[minimizersPos[i+k-1]] - rlePositions[minimizersPos[i+k-2]];
									//position_of_second_minimizer_seq = position_of_second_minimizer;
									//position_of_second_minimizer_seq += (rlePositions[minimizersPos[i+k-2] + l - 1] - rlePositions[minimizersPos[i+k-2]]);
									//exit(1);


									u_int16_t pos_last_minimizer = minimizersPos[indexSecondLastMinimizer] + l;
									if(rlePositions.size() > 0) pos_last_minimizer = rlePositions[pos_last_minimizer];
									//u_int16_t pos_last_minimizer = rlePositions[minimizersPos[indexSecondLastMinimizer] + l];
									seq_length_start = read_pos_end - pos_last_minimizer; //rlePositions[minimizersPos[i+k-1]] - rlePositions[minimizersPos[i+k-2]];
								}
								else{

									//u_int16_t pos_last_minimizer = rlePositions[minimizersPos[indexSecondMinimizer]];
									u_int16_t pos_last_minimizer = minimizersPos[indexSecondMinimizer];
									if(rlePositions.size() > 0) pos_last_minimizer = rlePositions[pos_last_minimizer];
									seq_length_start = pos_last_minimizer - read_pos_start;

									//cout << "todo" << endl;
									//seq_length_start = 
									//position_of_second_minimizer = rlePositions[minimizersPos[i+1]];// - rlePositions[minimizersPos[i]];
								}

								u_int16_t position_of_second_to_last_minimizer = 0;
								u_int16_t position_of_second_to_last_minimizer_seq = 0;
								if(isReversed){
									
									u_int16_t pos_last_minimizer = minimizersPos[indexSecondMinimizer];
									if(rlePositions.size() > 0) pos_last_minimizer = rlePositions[pos_last_minimizer];
									seq_length_end = pos_last_minimizer - read_pos_start;

									//position_of_second_to_last_minimizer = rlePositions[minimizersPos[i+1]]; // - rlePositions[minimizersPos[i]];
								}
								else{
									//position_of_second_to_last_minimizer = rlePositions[minimizersPos[indexSecondLastMinimizer]]; //rlePositions[minimizersPos[i+k-1]] - rlePositions[minimizersPos[i+k-2]];
									//position_of_second_to_last_minimizer_seq = position_of_second_to_last_minimizer;
									//position_of_second_to_last_minimizer_seq += (rlePositions[minimizersPos[i+k-2] + l - 1] - rlePositions[minimizersPos[i+k-2]]);

									//u_int16_t pos_last_minimizer = rlePositions[minimizersPos[indexSecondLastMinimizer] + l];
									u_int16_t pos_last_minimizer = minimizersPos[indexSecondLastMinimizer] + l;
									if(rlePositions.size() > 0) pos_last_minimizer = rlePositions[pos_last_minimizer];
									//cout << "lala: " << rlePositions.size() << " " << (minimizersPos[i+k-1] + l) << endl;
									//cout << read_pos_end << " " << pos_last_minimizer << endl;
									//position_of_second_to_last_minimizer_seq = rlePositions[minimizersPos[i+k-2] + l - 1];
									seq_length_end = read_pos_end - pos_last_minimizer; //rlePositions[minimizersPos[i+k-1]] - rlePositions[minimizersPos[i+k-2]];
									//position_of_second_to_last_minimizer_seq = length - seq_length_end;
									//seq_length_end = rlePositions[minimizersPos[i+k-1]] - rlePositions[minimizersPos[i+k-2]];
									//position_of_second_to_last_minimizer_seq =  length - seq_length_end + 1; //(rlePositions[minimizersPos[i+k-2]] - read_pos_start);//rlePositions[minimizersPos[i+k-2]] - read_pos_start;
								}

								//cout << read_pos_start << " " <<  read_pos_end << " " << length << "      " << seq_length_start << " " << seq_length_end << endl;
								//u_int16_t lala = seq_length_end;
								//lala -= (rlePositions[minimizersPos[i+k-2] + l])
								//cout << "lala: " << rlePositions[minimizersPos[i+k-2] + l] << endl;

								/*
								cout << "\t" << read_pos_start << endl;
								cout << "\t" << seq_length_start << " " << seq_length_end << endl;
								cout << "\t" << position_of_second_to_last_minimizer_seq << endl;
								cout << "\t" << read_pos_end << endl;
								*/
							
								//cout << minimizersPos[i] << endl;
								//cout << rlePositions.size() << endl;

                                


								//cout << "HI: " << read_pos_start << " " << read_pos_end << " " << position_of_second_to_last_minimizer << endl;
								//cout << "HIIIII: " << rlePositions[minimizersPos[i+k-1] + l] << " " << rlePositions[minimizersPos[i+k-1]] << endl;


								//cout << "HO: " << read_pos_start << " " << read_pos_end << " " << position_of_second_to_last_minimizer << endl;

								//cout << read_pos_start << endl;
								//cout << read_pos_end << endl;

								//u_int16_t length = (rlePositions[minimizersPos[i+k-1]] + 1 - rlePositions[minimizersPos[i]] + 1);
								//cout << "HAAAA: " << length << endl;
                                kminmersLength.push_back({read_pos_start, read_pos_end, length, isReversed, position_of_second_minimizer, position_of_second_to_last_minimizer, position_of_second_minimizer_seq, position_of_second_to_last_minimizer_seq, seq_length_start, seq_length_end});
                            }
							else{
                                kminmersLength.push_back({0, 0, 0, isReversed, 0, 0, 0, 0, 0, 0});
							}

							//if(currentMinimizerIndex[0] + 1 != currentMinimizerIndex[1] || currentMinimizerIndex[1] + 1 != currentMinimizerIndex[2]){
							//	cout << vec._kmers[0] << " " << vec._kmers[1] << " " << vec._kmers[2] << endl;
							//}
							//if(readIndex == 30539){
							//	cout << vec._kmers[0] << " " << vec._kmers[1] << " " << vec._kmers[2] << endl;
							//}

							prevVec = vec;
							kminmers.push_back(vec);
							break;
						}
					}

					j += 1;
				}

				if(hasPalindrome) break;
			}

			
			if(!hasPalindrome) break;
			kminmers.clear();
        	kminmersLength.clear();
		}


	}

	static void getKminmers_complete(const size_t k, const vector<u_int64_t>& minimizers, const vector<u_int64_t>& minimizersPos, vector<ReadKminmerComplete>& kminmers, int readIndex){

		bool hasPalindromeTotal = false;

        kminmers.clear();

		if(minimizers.size() < k) return;

		vector<bool> bannedPositions(minimizers.size(), false);



		while(true){


			bool hasPalindrome = false;
			KmerVec prevVec;

			int i_max = ((int)minimizers.size()) - (int)k + 1;
			for(int i=0; i<i_max; i++){


				if(bannedPositions[i]){
					continue;
				}
				KmerVec vec;
				vector<u_int32_t> currentMinimizerIndex;

				int j=i;
				while(true){
					
					if(j >= minimizers.size()) break;
					if(bannedPositions[j]){
						j += 1;
						continue;
					}

					u_int64_t minimizer = minimizers[j];

					vec._kmers.push_back(minimizer);
					currentMinimizerIndex.push_back(j);

					if(vec._kmers.size() == k){

	

						
						if(vec.isPalindrome() || (i > 0 && vec.normalize() == prevVec.normalize())){ //Palindrome: 121 (créé un cycle), Large palindrome = 122 221 (créé une tip)

							for(size_t m=0; m<k; m++){
								bannedPositions[currentMinimizerIndex[m]] = true;
								hasPalindromeTotal = true;
								//cout << "banned " << currentMinimizerIndex[m] << endl;
							}

							hasPalindrome = true;
							break;
						}
						else{

							bool isReversed;
							vec = vec.normalize(isReversed);

							u_int32_t indexFirstMinimizer = currentMinimizerIndex[0];
							u_int32_t indexSecondMinimizer = currentMinimizerIndex[1];
							u_int32_t indexSecondLastMinimizer = currentMinimizerIndex[currentMinimizerIndex.size()-2];
							u_int32_t indexLastMinimizer = currentMinimizerIndex[currentMinimizerIndex.size()-1];


							u_int32_t read_pos_start = indexFirstMinimizer;
							u_int32_t read_pos_end = indexLastMinimizer;
							u_int32_t length = read_pos_end - read_pos_start + 1;

							
							u_int32_t seq_length_start = 0;
							u_int32_t seq_length_end = 0;

							if(isReversed){

								u_int32_t pos_last_minimizer = indexSecondLastMinimizer;
								seq_length_start = read_pos_end - pos_last_minimizer;
							}
							else{

								u_int32_t pos_last_minimizer = indexSecondMinimizer;
								seq_length_start = pos_last_minimizer - read_pos_start;
							}

							if(isReversed){
								
								u_int32_t pos_last_minimizer = indexSecondMinimizer;
								seq_length_end = pos_last_minimizer - read_pos_start;

							}
							else{
								u_int32_t pos_last_minimizer = indexSecondLastMinimizer;
								seq_length_end = read_pos_end - pos_last_minimizer; 
							}

							//cout << read_pos_start << " " << read_pos_end << " " << isReversed << "     " << seq_length_start << " " << seq_length_end << endl;
							//vector<u_int64_t> kminmerMinimizers;

							//for(size_t i=indexFirstMinimizer; i<=indexLastMinimizer; i++){
							//	kminmerMinimizers.push_back(minimizers[i]);
							//}


							prevVec = vec;
							kminmers.push_back({vec, isReversed, read_pos_start, read_pos_end, seq_length_start, seq_length_end, length});
							break;
						}
					}

					j += 1;
				}

				if(hasPalindrome) break;
			}

			
			if(!hasPalindrome) break;
        	kminmers.clear();
		}

		//if(hasPalindromeTotal){
		//	cout << "mioum" << endl;
		//	getchar();
		//}
	}

};





class KmerDirect
{
public:
	/** Returns the value of the kmer.
	 * \return the kmer value as a Type object. */
	const u_int64_t& value  () const { return _value;   }

	/** This is a dummy function that always returns the value of the kmer in the forward direction, even when "which" is 1.
	 * it's provided for API compatibility with KmerCanonical
	 * We could compute revcomp(_value) when which=1, but KmerDirect doesn't know about its k-mer size.
	 * \param[in] which: dummy parameter
	 * \return the kmer value as a Type object. */
	const u_int64_t& value  (int which) const { if (which==1){ std::cout << "unsupported call to value(which) for KmerDirect" << std::endl; exit(1); } 
											return _value;   }

	/*_ Comparison operator between two instances.
		* \param[in] t : object to be compared to
		* \return true if the values are the same, false otherwise. */
	bool operator< (const KmerDirect& t) const  { return this->_value < t._value; };

	/** Set the value of the kmer
	 * \param[in] val : value to be set. */
	void set (const u_int64_t& val) { _value=val; }

	/** Tells whether the kmer is valid or not. It may be invalid if some unwanted
	 * nucleotides characters (like N) have been used to build it.
	 * \return true if valid, false otherwise. */
	bool isValid () const { return _isValid; }

	/* compatibility with KmerCanonical API */
	bool which () const { return true; }

	/* compatibility with KmerCanonical API */
	//Strand strand() const { return STRAND_FORWARD;  }

	/* compatibility with KmerCanonical API */
	const u_int64_t& forward() const { return value(); }

	/** Returns the reverse complement value of this canonical kmer.
	 * \return the reverse complement value */
	//const Type& revcomp() const { return  }


	u_int64_t _value;
	bool _isValid;
	//friend class ModelDirect;

	/** Extract a mmer from a kmer. This is done by using a mask on the kmer.
	 * \param[in] mask : mask to be applied to the current kmer
	 * \param[in] size : shift size (needed for some kmer classes but not all)
	 * \param[in] mmer_lut : lookup table of minimizers
	 * \return the extracted kmer.
	 */

	KmerDirect extract      (const u_int64_t& mask, size_t size, u_int64_t * mmer_lut)  {  KmerDirect output;  output.set (mmer_lut[(this->value() & mask)]);  return output;  }
	KmerDirect extractShift (const u_int64_t& mask, size_t size, u_int64_t * mmer_lut)  {  KmerDirect output = extract(mask,size,mmer_lut);  _value = _value >> 2;  return output;  }
};



class KmerCanonical
{
public:

	/** Returns the value of the kmer.
	 * \return the kmer value as a Type object. */
	const u_int64_t& value  () const { return table[(int)choice];   }

	/** Returns the value of the kmer.
	 * \param[in] which: forward or reverse strand
	 * \return the kmer value as a Type object. */
	const u_int64_t& value  (int which) const { return table[which];   }

	/** Comparison operator between two instances.
	 * \param[in] t : object to be compared to
	 * \return true if the values are the same, false otherwise. */
	bool operator< (const KmerDirect& t) const  { return this->value() < t.value(); };

	/** Set the value of the kmer. IMPORTANT: Not really a forward/revcomp couple,
	 * but may be useful for the minimizer default value.
	 * \param[in] val : value to be set (set to both forward and reverse complement). */
	void set (const u_int64_t& val)
	{
		table[0]=val;
		table[1]=val;
		choice = 0;
	}

	/** Set the forward/revcomp attributes. The canonical form is computed here.
	 * \param[in] forward : forward value
	 * \param[in] revcomp : reverse complement value.
	 */
	void set (const u_int64_t& forward, const u_int64_t& revcomp)
	{
		table[0]=forward;
		table[1]=revcomp;
		updateChoice ();
	}
	
	/** Tells whether the kmer is valid or not. It may be invalid if some unwanted
	 * nucleotides characters (like N) have been used to build it.
	 * \return true if valid, false otherwise. */
	bool isValid () const { return _isValid; }

	/** Returns the forward value of this canonical kmer.
	 * \return the forward value */
	const u_int64_t& forward() const { return table[0]; }

	/** Returns the reverse complement value of this canonical kmer.
	 * \return the reverse complement value */
	const u_int64_t& revcomp() const { return table[1]; }

	/** Tells which strand is used for the kmer.
	 * \return true if the kmer value is the forward value, false if it is the reverse complement value
	 */
	bool which () const { return choice==0 ? true : false; }

	/** Tells which strand is used.
	 * \return the used strand. */
	//Strand strand() const { return which() ? STRAND_FORWARD : STRAND_REVCOMP; }

	/* tells whether a kmer and its revcomp are identical */
	bool isPalindrome () const { return table[0] == table[1]; }

	u_int64_t table[2];  char choice;
	
	bool _isValid;
	void updateChoice () { choice = (table[0] < table[1]) ? 0 : 1; }
	//friend class ModelCanonical;

	/** Extract a mmer from a kmer. This is done by using a mask on the kmer.
	 * \param[in] mask : mask to be applied to the current kmer
	 * \param[in] size : shift size (needed for some kmer classes but not all)
	 * \param[in] mmer_lut : lookup table of minimizers
	 * \return the extracted kmer.
	 */
	KmerCanonical extract (const u_int64_t& mask, size_t size, u_int64_t* mmer_lut)
	{

		KmerCanonical output;
		
		output.set(mmer_lut[(this->table[0] & mask)]); //no need to recomp updateChoice with this
		//mmer_lut takes care of revcomp and forbidden mmers
		//output.set (this->table[0] & mask, (this->table[1] >> size) & mask);
		//output.updateChoice();
		return output;
	}


	KmerCanonical extractShift (const u_int64_t& mask, size_t size, u_int64_t * mmer_lut)
	{
		KmerCanonical output = extract (mask, size,mmer_lut);
		table[0] = table[0] >> 2;   table[1] = table[1] << 2;  updateChoice();
		return output;
	}
	
};

class KmerModel{
public:

    typedef std::pair<char,char> ConvertChar;
    struct ConvertASCII    { static ConvertChar get (const char* buffer, size_t idx)  { return ConvertChar((buffer[idx]>>1) & 3, (buffer[idx]>>3) & 1); }};

	
	size_t  _kmerSize;
	u_int64_t  _kmerMask;
	u_int64_t _revcompTable[4];

	typedef KmerCanonical Kmer;
	//typedef Kmer<span>::KmerCanonical Kmer;

	KmerModel(size_t kmerSize){

		_kmerSize = kmerSize;

		u_int64_t un;
		un = 1;
		_kmerMask = (un << (_kmerSize*2)) - un;

		size_t shift = 2*(_kmerSize-1);

		/** The _revcompTable is a shortcut used while computing revcomp recursively. */
		/** Important: don't forget the Type cast, otherwise the result in only on 32 bits. */
		for (size_t i=0; i<4; i++)   {  u_int64_t tmp; tmp = comp_NT[i];  _revcompTable[i] = tmp << shift;  }
	}

	//template<class Convert>
	int polynom (const char* seq, u_int64_t& kmer, size_t startIndex)  const
	{
		ConvertChar c;
		int badIndex = -1;

		kmer = 0;
		for (size_t i=0; i<_kmerSize; ++i)
		{
			//cout << seq[i] << endl;
			c = ConvertASCII::get(seq,i+startIndex);

			kmer = (kmer<<2) + c.first;

			if (c.second)  { badIndex = i; }
		}

		return badIndex;
	}

	inline static u_int64_t revcomp(const u_int64_t& x, size_t sizeKmer)
    {
        u_int64_t res = x;

        unsigned char* kmerrev  = (unsigned char *) (&(res));
        unsigned char* kmer     = (unsigned char *) (&(x));

        for (size_t i=0; i<8; ++i)  {  kmerrev[8-1-i] = revcomp_4NT [kmer[i]];  }

        return (res >> (2*( 32 - sizeKmer))) ;
    }

	u_int64_t reverse (const u_int64_t& kmer)  const  { return revcomp (kmer, this->_kmerSize); }
	/*
	std::string toString (u_int64_t val, size_t sizeKmer) const
    {
        char seq[sizeKmer+1];
        char bin2NT[4] = {'A','C','T','G'};

        for (size_t i=0; i<sizeKmer; i++)  {  seq[sizeKmer-i-1] = bin2NT [(*this)[i]];  }
        seq[sizeKmer]='\0';
        return seq;
    }*/

	bool iterate (const char* seq, size_t length, vector<u_int64_t>& tmpKmers) const{

		int32_t nbKmers = length - _kmerSize + 1;
		if (nbKmers <= 0)  { return false; }

		tmpKmers.resize(nbKmers);

		Kmer result;
		int indexBadChar = first(seq, result, 0);
		//int indexBadChar = static_cast<const ModelCanonical*>(this)->template first<Convert> (seq, result, 0);

		size_t idxComputed = 0;

		/*
		cout << result.value() << endl;

		u_int64_t kmerVal = result.value();
		Type un;
		un.setVal(kmerVal);

		cout << un.toString(_kmerSize) << endl;

		Type deux;
		deux.setVal(revcomp(kmerVal, _kmerSize));

		cout << deux.toString(_kmerSize) << endl;

		exit(1);
		//this->notification<Callback> (result, idxComputed, callback);
		*/
		tmpKmers[idxComputed] = result.value();
		if(!result.isValid()) tmpKmers[idxComputed] = -1; //Kmer with max value will be skipped as minimizers
		idxComputed += 1;
		
		for (size_t idx=_kmerSize; idx<length; idx++)
		{
			ConvertChar c = ConvertASCII::get (seq, idx);

			if (c.second)  { indexBadChar = _kmerSize-1; }
			else           { indexBadChar--;     }

			next(c.first, result, indexBadChar<0);
			tmpKmers[idxComputed] = result.value();
			if(!result.isValid()) tmpKmers[idxComputed] = -1; //Kmer with max value will be skipped as minimizers

			//if(!result.isValid()) cout << "img" << endl;
			//tmpKmers.push_back(result.value());
			//cout << result.value() << endl;
			//static_cast<const ModelCanonical*>(this)->template next<Convert> (c.first, result, indexBadChar<0);

			//this->notification<Callback> (result, ++idxComputed, callback);
			idxComputed += 1;
		}

		return true;
	}

	int first (const char* seq, Kmer& value, size_t startIndex)   const
	{

		int result = polynom(seq, value.table[0], startIndex);
		value._isValid = result < 0;
		value.table[1] = this->reverse (value.table[0]);
		value.updateChoice();
		return result;
	}

	void  next (char c, Kmer& value, bool isValid)   const
	{
		value.table[0] = ( (value.table[0] << 2) +  c                          ) & this->_kmerMask;
		value.table[1] = ( (value.table[1] >> 2) +  this->_revcompTable[(int)c]) & this->_kmerMask;
		value._isValid = isValid;

		value.updateChoice();
	}

	/*
	uint64_t getHash(const u_int64_t &k) const
	{
		return hash1(k, 0);
	}

	void getHash2(const u_int64_t &k) const
	{
		hash2(k, 1LL);
	}*/

};


class MinimizerParser{

public:

	u_int16_t _minimizerSize;
	KmerModel* _kmerModel;
	u_int32_t _seed;
	u_int64_t* _hash_otpt;
	double _minimizerBound;

	MinimizerParser(u_int16_t minimizerSize, double minimizerDensity){
		_minimizerSize = minimizerSize;
		_kmerModel = new KmerModel(_minimizerSize);
		_seed = 42;
		_hash_otpt = new u_int64_t[2];
		u_int64_t maxHashValue = -1;
		_minimizerBound = minimizerDensity * maxHashValue;
	}

	~MinimizerParser(){
		delete[] _hash_otpt;
	}

	void parse(const string& seq, vector<u_int64_t>& minimizers, vector<u_int64_t>& minimizersPos){

		minimizers.clear();
		minimizersPos.clear();

		//0 0
		//100000 3278692
		//200000 6550207
		//300000 9827077

		vector<u_int64_t> kmers;
		_kmerModel->iterate(seq.c_str(), seq.size(), kmers);

		if(kmers.size() == 0) return;

		for(u_int64_t pos=1; pos<kmers.size()-1; pos++){

			//cout << itKmer->value().getVal() << endl;
			//cout << itKmer->value() << endl;
			//cout << pos << " " << (rleSequence.size()-_minimizerSize) << endl;

			/*
			if(pos == 0){
				//cout << "lala1" << endl;
				//pos += 1;
				continue;
			}
			else if(pos == kmers.size()-1; //seq.size()-_minimizerSize){
				//cout << "lala2" << endl;
				continue;
			}*/

			//if(!itKmer->value().isValid()) continue;
			//kmer_type kmerMin = min(itKmer->value(), revcomp(itKmer->value(), _kmerSize));
			//if(lala < 100 ) cout << model.toString(itKmer->value()) << endl;
			//lala += 1;
			u_int64_t kmerValue = kmers[pos];
			MurmurHash3_x64_128 ((const char*)&kmerValue, sizeof(kmerValue), _seed, _hash_otpt);
			u_int64_t minimizer = _hash_otpt[0];



			//if(minimizerCounts[minimizer] > 1000) cout << minimizer << endl;
			//double kmerHashed_norm = ((double) minimizer) / maxHashValue;
			if(minimizer < _minimizerBound){//_minimizerDensity){


				minimizers.push_back(minimizer);
				minimizersPos.push_back(pos);

				
				//minimizerCounts[minimizer] += 1;
				//nbMinimizers += 1;
				//"reprise: systeme de functor pour les reads et kmers?"
				

			}
			
			//pos += 1;
		}

	}

};


class ReadParser{
public:

	string _inputFilename;
	bool _isFile;	
	size_t _l;
	size_t _k;
	float _density;
	vector<string> _filenames;
	u_int64_t _nbDatasets;

	ReadParser(const string& inputFilename, bool isFile){
		_inputFilename = inputFilename;
		_isFile = isFile;
		
		if(_isFile){
			_nbDatasets = 1;
			_filenames.push_back(inputFilename);
		}
		else{
			parseFilenames();
		}
	}


	ReadParser(const string& inputFilename, bool isFile, size_t l, size_t k, float density){
		_inputFilename = inputFilename;
		_isFile = isFile;
		_l = l;
		_k = k;
		_density = density;

		if(_isFile){
			_nbDatasets = 1;
			_filenames.push_back(inputFilename);
		}
		else{
			parseFilenames();
		}

	}

	void parseFilenames(){

		_nbDatasets = 0;

		std::ifstream infile(_inputFilename.c_str());
		std::string line;

		while (std::getline(infile, line)){

    		line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());

			fs::path path(line);

			if(fs::exists (path)){
				_nbDatasets += 1;
				_filenames.push_back(line);
			}
			else{
				cout << "File not found: " << line << endl;
			}
		}


	}



	void parse(const std::function<void(kseq_t*, u_int64_t)>& fun){

		u_int64_t readIndex = 0;

		for(const string& filename : _filenames){

			cout << filename << endl;

			gzFile fp;
			kseq_t *seq;
			int slen = 0, qlen = 0;
			fp = gzopen(filename.c_str(), "r");
			seq = kseq_init(fp);

			while (kseq_read(seq) >= 0){
				fun(seq, readIndex);
				readIndex += 1;
			}
				
			kseq_destroy(seq);
			gzclose(fp);
		}

		/*
		if(_isFile){
			gzFile fp;
			kseq_t *seq;
			int slen = 0, qlen = 0;
			fp = gzopen(_inputFilename.c_str(), "r");
			seq = kseq_init(fp);

			while (kseq_read(seq) >= 0){
				fun(seq, readIndex);
				readIndex += 1;
			}
				
			gzclose(fp);
		}
		else{
			std::ifstream infile(_inputFilename.c_str());
			std::string line;

			while (std::getline(infile, line))
			{
				cout << line << endl;

				gzFile fp;
				kseq_t *seq;
				int slen = 0, qlen = 0;
				fp = gzopen(line.c_str(), "r");
				seq = kseq_init(fp);

				while (kseq_read(seq) >= 0){
					fun(seq, readIndex);
					readIndex += 1;
				}
					
				gzclose(fp);
			}
		}
		*/

	}


	void parseKminmers(const std::function<void(vector<KmerVec>, vector<ReadKminmer>, u_int64_t, u_int64_t, string, string)>& fun){

		MinimizerParser* _minimizerParser = new MinimizerParser(_l, _density);

		u_int64_t readIndex = 0;
		u_int64_t datasetIndex = 0;

		for(const string& filename : _filenames){
			cout << filename << endl;

			readIndex = 0;

			gzFile fp;
			kseq_t *read;
			int slen = 0, qlen = 0;
			fp = gzopen(filename.c_str(), "r");
			read = kseq_init(fp);

			while (kseq_read(read) >= 0){

				string rleSequence;
				vector<u_int64_t> rlePositions;
				Encoder::encode_rle(read->seq.s, strlen(read->seq.s), rleSequence, rlePositions);

				vector<u_int64_t> minimizers;
				vector<u_int64_t> minimizers_pos;
				_minimizerParser->parse(rleSequence, minimizers, minimizers_pos);

				vector<KmerVec> kminmers; 
				vector<ReadKminmer> kminmersInfo;
				MDBG::getKminmers(_l, _k, minimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, readIndex, false);


				fun(kminmers, kminmersInfo, readIndex, datasetIndex, string(read->name.s, strlen(read->name.s)), string(read->seq.s, strlen(read->seq.s)));

				//cout << readIndex << endl;

				readIndex += 1;

				//if(readIndex > 50000) break;
			}
				
			kseq_destroy(read);
			gzclose(fp);

			datasetIndex += 1;

		}


		/*
		if(_isFile){
			gzFile fp;
			kseq_t *read;
			int slen = 0, qlen = 0;
			fp = gzopen(_inputFilename.c_str(), "r");
			read = kseq_init(fp);

			while (kseq_read(read) >= 0){

				//cout << readIndex << " " << read->name.s << endl;
				string rleSequence;
				vector<u_int64_t> rlePositions;
				Encoder::encode_rle(read->seq.s, strlen(read->seq.s), rleSequence, rlePositions);

				vector<u_int64_t> minimizers;
				vector<u_int64_t> minimizers_pos;
				_minimizerParser->parse(rleSequence, minimizers, minimizers_pos);

				vector<KmerVec> kminmers; 
				vector<ReadKminmer> kminmersInfo;
				MDBG::getKminmers(_l, _k, minimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, readIndex, false);

				fun(kminmers, kminmersInfo, readIndex, datasetIndex, string(read->name.s, strlen(read->name.s)));

				//fun(seq, readIndex);
				readIndex += 1;
				
			}
				
			gzclose(fp);
		}
		else{
			
			std::ifstream infile(_inputFilename.c_str());
			std::string line;

			while (std::getline(infile, line))
			{
				cout << line << endl;

				readIndex = 0;

				gzFile fp;
				kseq_t *read;
				int slen = 0, qlen = 0;
				fp = gzopen(line.c_str(), "r");
				read = kseq_init(fp);

				while (kseq_read(read) >= 0){

					string rleSequence;
					vector<u_int64_t> rlePositions;
					Encoder::encode_rle(read->seq.s, strlen(read->seq.s), rleSequence, rlePositions);

					vector<u_int64_t> minimizers;
					vector<u_int64_t> minimizers_pos;
					_minimizerParser->parse(rleSequence, minimizers, minimizers_pos);

					vector<KmerVec> kminmers; 
					vector<ReadKminmer> kminmersInfo;
					MDBG::getKminmers(_l, _k, minimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, readIndex, false);


					//fun(kminmers, kminmersInfo, readIndex, datasetIndex, string(read->name.s, strlen(read->name.s)));

					cout << readIndex << endl;

					readIndex += 1;

					//if(readIndex > 50000) break;
				}
					
				gzclose(fp);

				datasetIndex += 1;
			}
		}
		*/

		delete _minimizerParser;
	}

	void extractSubsample(const string& outputFilename, unordered_set<u_int64_t>& selectedReads){

		_extractSubsample_outputFile = gzopen(outputFilename.c_str(),"wb");
		_selectedReads = selectedReads;

		//ReadParser parser(readFilename, false);
		auto fp = std::bind(&ReadParser::extractSubsample_read, this, std::placeholders::_1, std::placeholders::_2);
		parse(fp);

		gzclose(_extractSubsample_outputFile);

	}

	gzFile _extractSubsample_outputFile;
	unordered_set<u_int64_t> _selectedReads;

	void extractSubsample_read(kseq_t* read, u_int64_t readIndex){

		if(_selectedReads.find(readIndex) != _selectedReads.end()){
			string header = ">" + string(read->name.s, strlen(read->name.s)) + "\n";
			string seq = string(read->seq.s, strlen(read->seq.s)) + "\n";
			gzwrite(_extractSubsample_outputFile, (const char*)&header[0], header.size());
			gzwrite(_extractSubsample_outputFile, (const char*)&seq[0], seq.size());
		}

	}

};






class KminmerParser{

public:

	string _inputFilename;
	size_t _l;
	size_t _k;

	KminmerParser(const string& inputFilename, size_t l, size_t k){
		_inputFilename = inputFilename;
		_l = l;
		_k = k;
	}

	void parse(const std::function<void(vector<KmerVec>, vector<ReadKminmer>, u_int64_t)>& fun){

		gzFile file_readData = gzopen(_inputFilename.c_str(),"rb");

		u_int64_t readIndex = 0;

		while(true){
			
			u_int16_t size;
			vector<u_int64_t> minimizers;
			vector<u_int16_t> minimizersPosOffsets; 
			gzread(file_readData, (char*)&size, sizeof(size));

			if(gzeof(file_readData)) break;
			
			minimizers.resize(size);
			minimizersPosOffsets.resize(size);
			gzread(file_readData, (char*)&minimizers[0], size * sizeof(u_int64_t));
			gzread(file_readData, (char*)&minimizersPosOffsets[0], size * sizeof(u_int16_t));

			
			//cout << "----" << endl;
			vector<u_int64_t> minimizersPos; 
			if(size > 0){
				u_int64_t pos = minimizersPosOffsets[0];
				minimizersPos.push_back(pos);
				for(size_t i=1; i<minimizersPosOffsets.size(); i++){
					pos += minimizersPosOffsets[i];
					minimizersPos.push_back(pos);
					//cout << minimizersPosOffsets[i] << " " << pos << endl;
				}
			}


			vector<KmerVec> kminmers; 
			vector<ReadKminmer> kminmersInfo;
			vector<u_int64_t> rlePositions;
			MDBG::getKminmers(_l, _k, minimizers, minimizersPos, kminmers, kminmersInfo, rlePositions, 0, false);

			fun(kminmers, kminmersInfo, readIndex);

			readIndex += 1;
		}

		gzclose(file_readData);

	}

	void parseMinspace(const std::function<void(vector<u_int64_t>, vector<ReadKminmerComplete>, u_int64_t)>& fun){

		gzFile file_readData = gzopen(_inputFilename.c_str(),"rb");

		u_int64_t readIndex = 0;

		while(true){
			
			u_int16_t size;
			vector<u_int64_t> minimizers;
			vector<u_int16_t> minimizersPosOffsets; 
			gzread(file_readData, (char*)&size, sizeof(size));

			if(gzeof(file_readData)) break;
			
			minimizers.resize(size);
			minimizersPosOffsets.resize(size);
			gzread(file_readData, (char*)&minimizers[0], size * sizeof(u_int64_t));
			gzread(file_readData, (char*)&minimizersPosOffsets[0], size * sizeof(u_int16_t));

			
			//cout << "----" << endl;
			vector<u_int64_t> minimizersPos; 
			if(size > 0){
				u_int64_t pos = minimizersPosOffsets[0];
				minimizersPos.push_back(pos);
				for(size_t i=1; i<minimizersPosOffsets.size(); i++){
					pos += minimizersPosOffsets[i];
					minimizersPos.push_back(pos);
					//cout << minimizersPosOffsets[i] << " " << pos << endl;
				}
			}


			vector<ReadKminmerComplete> kminmersInfo;
			MDBG::getKminmers_complete(_k, minimizers, minimizersPos, kminmersInfo, readIndex);

			fun(minimizers, kminmersInfo, readIndex);

			readIndex += 1;
		}

		gzclose(file_readData);

	}

	void parseDuo(const std::function<void(vector<KmerVec>, vector<ReadKminmer>, u_int64_t, vector<KmerVec>, vector<ReadKminmer>)>& fun){

		gzFile file_readData = gzopen(_inputFilename.c_str(),"rb");

		u_int64_t readIndex = 0;

		while(true){
			
			u_int16_t size;
			vector<u_int64_t> minimizers;
			vector<u_int16_t> minimizersPosOffsets; 
			gzread(file_readData, (char*)&size, sizeof(size));

			if(gzeof(file_readData)) break;
			
			minimizers.resize(size);
			minimizersPosOffsets.resize(size);
			gzread(file_readData, (char*)&minimizers[0], size * sizeof(u_int64_t));
			gzread(file_readData, (char*)&minimizersPosOffsets[0], size * sizeof(u_int16_t));

			
			//cout << "----" << endl;
			vector<u_int64_t> minimizersPos; 
			if(size > 0){
				u_int64_t pos = minimizersPosOffsets[0];
				minimizersPos.push_back(pos);
				for(size_t i=1; i<minimizersPosOffsets.size(); i++){
					pos += minimizersPosOffsets[i];
					minimizersPos.push_back(pos);
					//cout << minimizersPosOffsets[i] << " " << pos << endl;
				}
			}



			vector<KmerVec> kminmers; 
			vector<ReadKminmer> kminmersInfo;
			vector<u_int64_t> rlePositions;
			MDBG::getKminmers(_l, _k, minimizers, minimizersPos, kminmers, kminmersInfo, rlePositions, 0, false);


			vector<KmerVec> kminmers2; 
			vector<ReadKminmer> kminmersInfo2;
			vector<u_int64_t> rlePositions2;
			vector<u_int64_t> minimizersPos2; 
			MDBG::getKminmers(_l, 3, minimizers, minimizersPos2, kminmers2, kminmersInfo2, rlePositions2, 0, false);

			fun(kminmers, kminmersInfo, readIndex, kminmers2, kminmersInfo2);

			readIndex += 1;
		}

		gzclose(file_readData);

	}

	void parse_mContigs(const std::function<void(vector<KmerVec>, vector<ReadKminmer>, u_int64_t)>& fun){

		gzFile file_mContigs = gzopen(_inputFilename.c_str(), "rb");

		u_int64_t readIndex = 2000000000;

		while(true){
			
			//cout << readIndex << endl;
			u_int64_t size;
			gzread(file_mContigs, (char*)&size, sizeof(size));
			
			if(gzeof(file_mContigs)) break;

			vector<u_int64_t> minimizers;
			minimizers.resize(size);
			gzread(file_mContigs, (char*)&minimizers[0], size * sizeof(u_int64_t));

			//if(minimizers.size() < 50) continue;

			vector<u_int64_t> minimizersPos(minimizers.size(), 0);
			for(size_t i=0; i<minimizersPos.size(); i++){
				minimizersPos[i] = 200 +i*200;
			}


			vector<KmerVec> kminmers; 
			vector<ReadKminmer> kminmersInfos;
			vector<u_int64_t> rlePositions;
			MDBG::getKminmers(_l, _k, minimizers, minimizersPos, kminmers, kminmersInfos, rlePositions, 0, false);

			fun(kminmers, kminmersInfos, readIndex);

			readIndex += 1;
		}

		gzclose(file_mContigs);
	
	}

	void parse_mContigs_minSpace(const std::function<void(vector<u_int64_t>, vector<ReadKminmerComplete>, u_int64_t)>& fun){

		gzFile file_mContigs = gzopen(_inputFilename.c_str(), "rb");

		u_int64_t readIndex = 2000000000;

		while(true){
			
			//cout << readIndex << endl;
			u_int64_t size;
			gzread(file_mContigs, (char*)&size, sizeof(size));
			
			if(gzeof(file_mContigs)) break;

			vector<u_int64_t> minimizers;
			minimizers.resize(size);
			gzread(file_mContigs, (char*)&minimizers[0], size * sizeof(u_int64_t));

			//if(minimizers.size() < 50) continue;

			vector<u_int64_t> minimizersPos(minimizers.size(), 0);
			for(size_t i=0; i<minimizersPos.size(); i++){
				minimizersPos[i] = 200 +i*200;
			}



			vector<ReadKminmerComplete> kminmersInfo;
			MDBG::getKminmers_complete(_k, minimizers, minimizersPos, kminmersInfo, readIndex);

			fun(minimizers, kminmersInfo, readIndex);
		}

		gzclose(file_mContigs);
	
	}

};

class Tool{

public:

	void run(int argc, char* argv[]){
		parseArgs(argc, argv);
		execute();
	}

    virtual void parseArgs (int argc, char* argv[]) = 0;
    virtual void execute () = 0;
};

#endif