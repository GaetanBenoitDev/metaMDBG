/*

- add "polish" dans les commandes du tool si ça fonctionne bien
- tester les différents paramètres (surtout minimizer size)

TODO:
- toBasespace: il faudrait utiliser le processus de mapping pour tous les kminmer en fait, pas seulement ceux repeter (ne aps etre repeter dasn les contigs ne veux pas dire que la sequence du kminmer n'est pas ambigue dans les reads)
- Argument a changer dans COntigPolisheer (liste de read en argument, fichier de sortie en argument)

Paralelisation:
	- ToBasespace: read indexing phmap
	- ToMinsapce: phmap

GenerateContigs:
	- a test: pas besoin de générer le graph initial (t=1), de le compacter et générer les contigs (normalement on veut pas générer les side des bubbles enlevé etc)

- ajuster valeur par defaut de l (minimizer size, 13 ou 21 a test)

- Polisher:  si on a plus besoin de low mmeory correction, enlever l'exe "mapper" eventuellement, et utiliser gzip pour compresser els sortie de minimap2 (on en aura besoin dans derep je pense) 
- enelever tous els code de abpoa et spoa, et sortir tous ce qui est contigpolisher
- Add tool version in program help
- enelever cxxopts (old arg system)

ContigPolisher:
	- determiner le coverage d'un contig, utiliser cette valeur pour estimater la memory total du contig _windowByteSize = (contigLength*windowLength*NbWindows) (puis _windowByteSize a enelever)
	- flag "c" rempalcer par "l" parfois a check

*/

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
#include "./utils/args.hxx"
#include <random>
#include <chrono>
#include "utils/parallel_hashmap/phmap.h"
#include "utils/kmer/Kmer.hpp"
#include <omp.h>
namespace fs = std::filesystem;
using namespace std::chrono;
//using std::chrono::high_resolution_clock;
//using std::chrono::duration_cast;
//using std::chrono::duration;
//using std::chrono::milliseconds;
//namespace fs = std::filesystem;

typedef unsigned __int128 u_int128_t;

//#define FIX_PALINDROME

//#include <sys/types.h>
//#include <sys/stat.h>

#include <functional>
#include <memory>
//#include "./utils/ntHashIterator.hpp"
#include "./utils/kseq.h"
#include "./utils/DnaBitset.hpp"
KSEQ_INIT(gzFile, gzread)

using namespace std;

struct Read{
	u_int64_t _index;
	string _header;
	string _seq;
	string _qual;
};



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

const string ARG_INPUT_FILENAME = "i";
const string ARG_INPUT_FILENAME_TRUTH = "itruth";
const string ARG_OUTPUT_DIR = "o";
const string ARG_OUTPUT_FILENAME = "f";
const string ARG_MINIMIZER_LENGTH = "l";
const string ARG_KMINMER_LENGTH = "k";
const string ARG_MINIMIZER_DENSITY = "d";
const string ARG_DEBUG = "debug";
const string ARG_FINAL = "final";
const string ARG_INPUT_FILENAME_CONTIG = "c";
const string ARG_INPUT_FILENAME_CONTIG_FASTA = "cf";
const string ARG_INPUT_FILENAME_UNITIG_NT = "unitigNt";
const string ARG_INPUT_FILENAME_UNITIG_CLUSTER = "cluster";
const string ARG_INPUT_FILENAME_ABUNDANCE = "a";
const string ARG_INPUT_FILENAME_BINNING = "bi";
const string ARG_OUTPUT_FILENAME_BINNING = "bo";
const string ARG_FIRST_PASS = "firstpass";
const string ARG_FASTA = "fasta";
const string ARG_NB_CORES = "t";
const string ARG_EVAL = "eval";
const string ARG_BLOOM_FILTER = "nofilter";

const string NB_CORES_DEFAULT = "3";
const int NB_CORES_DEFAULT_INT = 3;
//const string FILENAME_NO_KMINMER_READS = "reads_noKminmers.bin";




const char ARG_INPUT_FILENAME2 = 'i';
const char ARG_OUTPUT_DIR2 = 'o';
const char ARG_OUTPUT_FILENAME2 = 'f';
const char ARG_MINIMIZER_LENGTH2 = 'l';
const char ARG_KMINMER_LENGTH2 = 'k';
const char ARG_MINIMIZER_DENSITY2 = 'd';
const char ARG_INPUT_FILENAME_CONTIG2 = 'c';
const char ARG_INPUT_FILENAME_ABUNDANCE2 = 'a';
const char ARG_NB_CORES2 = 't';


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
	u_int32_t _length;
	bool _isReversed;
	u_int32_t _position_of_second_minimizer;
	u_int32_t _position_of_second_to_last_minimizer;
	u_int32_t _position_of_second_minimizer_seq;
	u_int32_t _position_of_second_to_last_minimizer_seq;
	u_int32_t _seq_length_start;
	u_int32_t _seq_length_end;
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

		if(_from < _to){
			return edge1;
		}
		else{
			return edge2;
		}
		/*
        DbgEdge edge1 = {_from, _to};
        DbgEdge edge2 = {_to, _from};
		size_t e1 = edge1.hash();
		size_t e2 = edge2.hash();
		if(e1 == e2){
			if(_from < _to){
				return edge1;
			}
			else{
				return edge2;
			}
		}
        else if(e1 < e2){
            return edge1;
        }
        else{
            return edge2;
        }
		*/
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
	u_int32_t _quality;
	//u_int16_t _length;
	//u_int16_t _overlapLength_start;
	//u_int16_t _overlapLength_end;
	//bool _isReversed;
	//u_int32_t _unitigNbNodes;
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
		/*
		KmerVec vec_reverse = reverse();

		if(h() < vec_reverse.h()){
			return *this;
		}
		else{
			return vec_reverse;
		}
		*/
		bool dummy;
		return normalize(dummy);
		
	}

	KmerVec normalize(bool& isReversed){

		KmerVec vec_reverse = reverse();
		
		for(size_t i=0; i<_kmers.size(); i++){
			if(_kmers[i] == vec_reverse._kmers[i]){
				continue;
			}
			else if(_kmers[i] < vec_reverse._kmers[i]){
				isReversed = false;
				return *this;
			}
			else{
				isReversed = true;
				return vec_reverse;
			}
		}

		isReversed = true;
		return vec_reverse;

		//if(h() < vec_reverse.h()){
		//	isReversed = false;
		//	return *this;
		//}
		//else{
		//	isReversed = true;
		//	return vec_reverse;
		//}
		
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
	
	/*
    inline static u_int64_t simplehash16_64(u_int64_t key, int  shift)
    {
        u_int64_t input = key >> shift;
        u_int64_t res = random_values[input & 255]   ;
        
        input = input  >> 8;
        res  ^= random_values[input & 255] ;

        
        res  ^= random_values[key & 255] ;

        return res;
    }
	*/

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
	u_int8_t _quality;
};

struct ReadNodeName{
	u_int32_t _nodeIndex;
	u_int64_t _supportingReadIndex;

	bool operator==(const ReadNodeName &other) const{
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
	struct hash<ReadNodeName>{
		std::size_t operator()(const ReadNodeName& k) const{
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

struct KminmerEdge2{
	u_int32_t _nodeName;
	u_int64_t _minimizer;
	bool _isReversed;
	bool _isPrefix;
};

struct KminmerEdge{
	u_int32_t _nodeName;
	KmerVec _vec;
};

//typedef phmap::parallel_flat_hash_set<KmerVec, phmap::priv::hash_default_hash<KmerVec>, phmap::priv::hash_default_eq<KmerVec>, std::allocator<KmerVec>, 4, std::mutex> KmerVecSet;
typedef phmap::parallel_flat_hash_map<KmerVec, DbgNode, phmap::priv::hash_default_hash<KmerVec>, phmap::priv::hash_default_eq<KmerVec>, std::allocator<std::pair<KmerVec, DbgNode>>, 4, std::mutex> MdbgNodeMap;
typedef phmap::parallel_flat_hash_map<KmerVec, vector<KminmerEdge>, phmap::priv::hash_default_hash<KmerVec>, phmap::priv::hash_default_eq<KmerVec>, std::allocator<std::pair<KmerVec, vector<KminmerEdge>>>, 4, std::mutex> MdbgEdgeMap;
typedef phmap::parallel_flat_hash_map<KmerVec, vector<KminmerEdge2>, phmap::priv::hash_default_hash<KmerVec>, phmap::priv::hash_default_eq<KmerVec>, std::allocator<std::pair<KmerVec, vector<KminmerEdge2>>>, 4, std::mutex> MdbgEdgeMap2;


//unordered_map<KmerVec, DbgNode> _dbg_nodes;
//unordered_map<KmerVec, vector<KmerVec>> _dbg_edges;

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

struct KminmerList{
	u_int64_t _readIndex;
	vector<u_int64_t> _readMinimizers;
	vector<ReadKminmerComplete> _kminmersInfo;
	Read _read;
	//vector<KmerVec> _kminmers;
	//vector<ReadKminmer> _kminmersInfo;
};

class Commons{

public:

	static void createInputFile(auto& paths, const string& filename){

		//string inputFilenames = _inputFilename;
		ofstream inputFile(filename);

		//cout << inputFilenames << endl;
		//vector<string>* fields = new vector<string>();
		//GfaParser::tokenize(inputFilenames, fields, ',');

		for(auto &&path : paths){
			//string filenameAbs = fs::absolute(filename);
			//cout << path << endl; //" " << fs::canonical(filename) << " " << fs::weakly_canonical(filename) << endl;
			inputFile << path << endl;
		}

		//delete fields;
		inputFile.close();


	}

};

class Utils{




public:

	static string shortenHeader(const string& header){

		string shortenName;

		auto find = header.find(' ');
		if(find == std::string::npos){
			shortenName = header;
		}
		else{
			shortenName = header.substr(0, find);
		}

		return shortenName;
	}

	static u_int64_t contigName_to_contigIndex(const string& header){
		string name = header;
		size_t pos = name.find("ctg");
		name.erase(pos, 3);

		//remove circular or linear indicator
		if(name[name.size()-1] == 'l' || name[name.size()-1] == 'c'){
			name.pop_back(); 
		}

		u_int64_t contigIndex = stoull(name);
		return contigIndex;
	}

	static void toReverseComplement(string& seq) {

		//string seq_rc;

		//seq_rc.clear();
		//seq_rc.reserve(seq.size());

		size_t size = seq.size();

		std::reverse(seq.begin(), seq.end());
		for (std::size_t i = 0; i < size; ++i){
			
			switch (seq[i]){
			case 'A':
				seq[i] = 'T';
				break;    
			case 'C':
				seq[i] = 'G';
				break;
			case 'G':
				seq[i] = 'C';
				break;
			case 'T':
				seq[i] = 'A';
				break;
			//default:
			//	seq[i] = seq[i];
			//	break;
			}
		}

		/*
		for (int32_t i = seq.size() - 1; i >= 0; --i) {
			switch (seq[i]) {
				case 'A':
					seq_rc += 'T';
					break;
				case 'T':
					seq_rc += 'A';
					break;
				case 'C':
					seq_rc += 'G';
					break;
				case 'G':
					seq_rc += 'C';
					break;
				default:
					seq_rc += seq[i];
					break;
			}
		}

		//reverse_quality_.clear();
		//reverse_quality_.reserve(quality_.size());

		//for (int32_t i = quality_.size() - 1; i >= 0; --i) {
		//	reverse_quality_ += quality_[i];
		//}
		*/
	}


	static void executeCommand(const string& command, const string& outputDir){

		static double _maxMemoryUsage = 0;
		static string s = "Maximumresidentsetsize(kbytes):";

		cout << command << endl;

		string command2 = "{ time -v " + command + "; } 2> " + outputDir + "/time.txt";

		int ret = system(command2.c_str());
		if(ret != 0){
			cerr << "Command failed: " << ret << endl;
			exit(ret);
		}

		ifstream infile(outputDir + "/time.txt");
		string line;
		while(std::getline(infile, line)){

    		//line.erase(std::remove_if(line.begin(), line.end(), ::istab), line.end());
    		line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
			if(line.empty()) continue;

			size_t pos = line.find(s);
			if (pos != std::string::npos){
				line.erase(pos, s.length());
			}
			else{
				continue;
			}

			double maxMem = stod(line);
			maxMem /= 1024;
			maxMem /= 1024;
			_maxMemoryUsage = max(_maxMemoryUsage, maxMem);
		}
		infile.close();
		//cout << "Max memory: " << _maxMemoryUsage << " GB" << endl;
		//getchar();

		ofstream outfile(outputDir + "/perf.txt");
		outfile << "Peak memory (GB): " << _maxMemoryUsage << endl;
		outfile.close();

	}

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

	
	template<typename T>
	static u_int64_t collectSharedElements(const vector<T>& reads1, const vector<T>& reads2, unordered_set<T>& sharedReads){

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

	template<typename T>
	static u_int64_t collectSharedElements(const vector<T>& reads1, const vector<T>& reads2, vector<T>& sharedReads){

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

	static u_int64_t computeSharedReads(const vector<u_int32_t>& nodeNames, const vector<UnitigData>& unitigDatas){

		vector<u_int64_t> sharedReads;
		for(u_int64_t readIndex : unitigDatas[nodeNames[0]]._readIndexes){
			sharedReads.push_back(readIndex);
		}

		for(size_t i=1; i<nodeNames.size(); i++){

			vector<u_int64_t> sharedReadsTmp;
			collectSharedReads(sharedReads, unitigDatas[nodeNames[i]]._readIndexes, sharedReadsTmp);
			sharedReads = sharedReadsTmp;
		}

		return sharedReads.size();
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

	static bool shareAny(const vector<u_int64_t>& utg1, const vector<u_int64_t>& utg2){

		if(utg1.size() == 0 || utg2.size() == 0) return true;
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

		while(i < utg1.size() && j < utg2.size()){

			//cout << i << " " << j << endl;
			if(utg1[i] == utg2[j]){
				return true;
			}
			else if(utg1[i] < utg2[j]){
				i += 1;
			}
			else{
				j += 1;
			}

		}

		return false;
	}

};


class EncoderRLE{

public: 

	void execute(const char* sequence, size_t length, string& rleSequence, vector<u_int64_t>& rlePositions) {
		/*
		string sequence_str;

				char lastChar = '0';
				for(size_t i=0; i<sequence.getDataSize(); i++){
					if(readseq[i] == lastChar) continue;
					sequence_str += readseq[i];
					lastChar = readseq[i];
				}
		*/
		rlePositions.size();
		rleSequence.size();
		
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
	//unordered_map<KmerVec, DbgNode> _dbg_nodes;
	MdbgNodeMap _dbg_nodes;
	MdbgEdgeMap _dbg_edges;
	//unordered_map<KmerVec, vector<KmerVec>> _dbg_edges;

	MDBG(size_t k){
		_k = k;
		_node_id = 0;
	}

	void addNode(const KmerVec& vec, u_int16_t length, u_int16_t overlapLength_start, u_int16_t overlapLength_end, bool isReversed){

		/*
		if(_dbg_nodes.find(vec) != _dbg_nodes.end()){
			//if(_dbg_nodes[vec]._abundance > 1000){
			//	cout << _dbg_nodes[vec]._index << " " << _dbg_nodes[vec]._abundance << endl;
			//}
			_dbg_nodes[vec]._abundance += 1;
			return;
		}

		//Abundance = 2 since we first check that the kminmer is not unique
		DbgNode node = {_node_id, 2, length, overlapLength_start, overlapLength_end, isReversed};
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
		*/

	}
	
	/*
	void dump(const string& filename){
		gzFile file = gzopen(filename.c_str(),"wb");

		//bool lala = true;
		for(auto it : _dbg_nodes){
			const KmerVec& vec = it.first;
			const DbgNode& node = it.second;

			gzwrite(file, (const char*)&vec._kmers[0], _k * sizeof(u_int64_t));
			//cout << sizeof(DbgNode) << endl;
			gzwrite(file, (const char*)&node, sizeof(DbgNode));



		}

		gzclose(file);
	}
	*/

	void dump(const string& filename){

		ofstream kminmerFile = ofstream(filename);

		for(auto& it : _dbg_nodes){
			KmerVec vec = it.first;

			//it.second._index = _nodeName_to_deterministicNodeName[it.second._index];
			
			//KmerVecSorterData& d = kmerVecs[i];
			u_int32_t nodeName = it.second._index;
			u_int32_t abundance = it.second._abundance;
			u_int32_t quality = it.second._quality;

			//if(quality==0){
			//	cout << "omg " << nodeName << endl;
			//	getchar();
			//}
			//bool isReversed;
			//d._kmerVec.normalize(isReversed);

			//vector<u_int64_t> minimizerSeq = d._kmerVec.normalize()._kmers;

			vector<u_int64_t> minimizerSeq = vec._kmers;
			//if(kminmerInfo._isReversed){
			//	std::reverse(minimizerSeq.begin(), minimizerSeq.end());
			//}

			u_int16_t size = minimizerSeq.size();
			//_kminmerFile.write((const char*)&size, sizeof(size));
			kminmerFile.write((const char*)&minimizerSeq[0], size*sizeof(uint64_t));

			kminmerFile.write((const char*)&nodeName, sizeof(nodeName));
			kminmerFile.write((const char*)&abundance, sizeof(abundance));
			kminmerFile.write((const char*)&quality, sizeof(quality));
		}

		kminmerFile.close();

	}

	void load(const string& filename, bool removeUniqueKminmer){

		ifstream kminmerFile(filename);

		while (true) {

			u_int16_t size = _k;
			//kminmerFile.read((char*)&size, sizeof(size));


			vector<u_int64_t> minimizerSeq;
			minimizerSeq.resize(size);
			kminmerFile.read((char*)&minimizerSeq[0], size*sizeof(u_int64_t));

			if(kminmerFile.eof())break;

			u_int32_t nodeName;
			u_int32_t abundance;
			u_int32_t quality;
			//bool isReversed = false;

			kminmerFile.read((char*)&nodeName, sizeof(nodeName));
			kminmerFile.read((char*)&abundance, sizeof(abundance));
			kminmerFile.read((char*)&quality, sizeof(quality));


			//if(removeUniqueKminmer && abundance == 1) continue;
			KmerVec vec;
			vec._kmers = minimizerSeq;

			_dbg_nodes[vec] = {nodeName, abundance, quality};
		}

		kminmerFile.close();
	}
	
	/*
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


			_node_id += 1;
		}

		gzclose(file);
	}
	*/



	static void getKminmers(const size_t l, const size_t k, const vector<u_int64_t>& minimizers, const vector<u_int64_t>& minimizersPos, vector<KmerVec>& kminmers, vector<ReadKminmer>& kminmersLength, const vector<u_int64_t>& rlePositions, int readIndex, bool allowPalindrome){

		kminmers.clear();
        kminmersLength.clear();
        bool doesComputeLength = minimizersPos.size() > 0;
		if(minimizers.size() < k) return;

		#ifdef FIX_PALINDROME
		/*
		bool isLala = false;
		if(std::find(minimizers.begin(), minimizers.end(), 56801217747741349) !=  minimizers.end()){
			isLala = true;
		}*/
		
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
		#else



		int i_max = ((int)minimizers.size()) - (int)k + 1;
		for(int i=0; i<i_max; i++){

			KmerVec vec;
			vector<u_int32_t> currentMinimizerIndex;

			int j=i;
			while(true){
				
				if(j >= minimizers.size()) break;

				u_int64_t minimizer = minimizers[j];

				vec._kmers.push_back(minimizer);
				currentMinimizerIndex.push_back(j);

				if(vec._kmers.size() == k){

					bool isReversed;
					vec = vec.normalize(isReversed);

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

						u_int32_t length = read_pos_end - read_pos_start;

						//cout << read_pos_start << " " << read_pos_end << " " << 
						u_int32_t seq_length_start = 0;
						u_int32_t seq_length_end = 0;

						u_int32_t position_of_second_minimizer = 0;
						u_int32_t position_of_second_minimizer_seq = 0;
						if(isReversed){

							u_int32_t pos_last_minimizer = minimizersPos[indexSecondLastMinimizer] + l;
							if(rlePositions.size() > 0) pos_last_minimizer = rlePositions[pos_last_minimizer];
							//u_int16_t pos_last_minimizer = rlePositions[minimizersPos[indexSecondLastMinimizer] + l];
							seq_length_start = read_pos_end - pos_last_minimizer; //rlePositions[minimizersPos[i+k-1]] - rlePositions[minimizersPos[i+k-2]];
						}
						else{

							//u_int16_t pos_last_minimizer = rlePositions[minimizersPos[indexSecondMinimizer]];
							u_int32_t pos_last_minimizer = minimizersPos[indexSecondMinimizer];
							if(rlePositions.size() > 0) pos_last_minimizer = rlePositions[pos_last_minimizer];
							seq_length_start = pos_last_minimizer - read_pos_start;

							//cout << "todo" << endl;
							//seq_length_start = 
							//position_of_second_minimizer = rlePositions[minimizersPos[i+1]];// - rlePositions[minimizersPos[i]];
						}

						u_int32_t position_of_second_to_last_minimizer = 0;
						u_int32_t position_of_second_to_last_minimizer_seq = 0;
						if(isReversed){
							
							u_int32_t pos_last_minimizer = minimizersPos[indexSecondMinimizer];
							if(rlePositions.size() > 0) pos_last_minimizer = rlePositions[pos_last_minimizer];
							seq_length_end = pos_last_minimizer - read_pos_start;

							//position_of_second_to_last_minimizer = rlePositions[minimizersPos[i+1]]; // - rlePositions[minimizersPos[i]];
						}
						else{
							//position_of_second_to_last_minimizer = rlePositions[minimizersPos[indexSecondLastMinimizer]]; //rlePositions[minimizersPos[i+k-1]] - rlePositions[minimizersPos[i+k-2]];
							//position_of_second_to_last_minimizer_seq = position_of_second_to_last_minimizer;
							//position_of_second_to_last_minimizer_seq += (rlePositions[minimizersPos[i+k-2] + l - 1] - rlePositions[minimizersPos[i+k-2]]);

							//u_int16_t pos_last_minimizer = rlePositions[minimizersPos[indexSecondLastMinimizer] + l];
							u_int32_t pos_last_minimizer = minimizersPos[indexSecondLastMinimizer] + l;
							if(rlePositions.size() > 0) pos_last_minimizer = rlePositions[pos_last_minimizer];
							//cout << "lala: " << rlePositions.size() << " " << (minimizersPos[i+k-1] + l) << endl;
							//cout << read_pos_end << " " << pos_last_minimizer << endl;
							//position_of_second_to_last_minimizer_seq = rlePositions[minimizersPos[i+k-2] + l - 1];
							seq_length_end = read_pos_end - pos_last_minimizer; //rlePositions[minimizersPos[i+k-1]] - rlePositions[minimizersPos[i+k-2]];
							//position_of_second_to_last_minimizer_seq = length - seq_length_end;
							//seq_length_end = rlePositions[minimizersPos[i+k-1]] - rlePositions[minimizersPos[i+k-2]];
							//position_of_second_to_last_minimizer_seq =  length - seq_length_end + 1; //(rlePositions[minimizersPos[i+k-2]] - read_pos_start);//rlePositions[minimizersPos[i+k-2]] - read_pos_start;
						}
						//cout << i << ": " << seq_length_start << " " << seq_length_end << endl;
						kminmersLength.push_back({read_pos_start, read_pos_end, length, isReversed, position_of_second_minimizer, position_of_second_to_last_minimizer, position_of_second_minimizer_seq, position_of_second_to_last_minimizer_seq, seq_length_start, seq_length_end});
					}
					else{
						kminmersLength.push_back({0, 0, 0, isReversed, 0, 0, 0, 0, 0, 0});
					}

					kminmers.push_back(vec);
					break;
					
				}

				j += 1;
			}

		}

			
		
		#endif


	}

	static void getKminmers_complete(const size_t k, const vector<u_int64_t>& minimizers, const vector<u_int64_t>& minimizersPos, vector<ReadKminmerComplete>& kminmers, int readIndex, const vector<u_int8_t>& minimizerQualities){

        kminmers.clear();
		if(minimizers.size() < k) return;


		#ifdef FIX_PALINDROME
		bool hasPalindromeTotal = false;



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
		#else

		int i_max = ((int)minimizers.size()) - (int)k + 1;
		for(int i=0; i<i_max; i++){


			KmerVec vec;
			vector<u_int32_t> currentMinimizerIndex;

			int j=i;
			while(true){
				
				if(j >= minimizers.size()) break;

				u_int64_t minimizer = minimizers[j];

				vec._kmers.push_back(minimizer);
				currentMinimizerIndex.push_back(j);

				if(vec._kmers.size() == k){

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

					u_int8_t minQuality = -1;
					for(size_t i=indexFirstMinimizer; i<=indexLastMinimizer; i++){
						if(minimizerQualities[i] < minQuality){
							minQuality = minimizerQualities[i];
						}
					}

					kminmers.push_back({vec, isReversed, read_pos_start, read_pos_end, seq_length_start, seq_length_end, length, minQuality});
					break;
				}

				j += 1;
			}

		}

			
		#endif

		//if(hasPalindromeTotal){
		//	cout << "mioum" << endl;
		//	getchar();
		//}
	}

};




class ReadParserParallel{

public:

	string _inputFilename;
	bool _isFile;
	bool _isBitset;
	size_t _l;
	size_t _k;
	float _density;
	vector<string> _filenames;
	u_int64_t _nbDatasets;
	int _nbCores;

	ReadParserParallel(const string& inputFilename, bool isFile, bool isBitset, int nbCores){

		if(!fs::exists(inputFilename)){
			cout << "File not found: " << inputFilename << endl;
			exit(1);
		}

		_inputFilename = inputFilename;
		_isFile = isFile;
		_isBitset = isBitset;
		_nbCores = nbCores;
		
		if(_isFile){
			_nbDatasets = 1;
			_filenames.push_back(inputFilename);
		}
		else{
			parseFilenames();
		}
	}


	ReadParserParallel(const string& inputFilename, bool isFile, size_t l, size_t k, float density){
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
			if(line.empty()) continue;

			fs::path path(line);

			if(fs::exists (path)){
				_nbDatasets += 1;
				_filenames.push_back(line);
			}
			else{
				cout << "File not found: " << line << endl;
				exit(1);
			}
		}


	}



	template<typename Functor>
	void parse(const Functor& functor){

		//#pragma omp for
		//cout << _nbCores << endl;

		u_int64_t readIndex = -1;

		for(const string& filename : _filenames){

			cout << filename << endl;

			gzFile fp = gzopen(filename.c_str(), "r");
			kseq_t *seq;
			seq = kseq_init(fp);

			#pragma omp parallel num_threads(_nbCores)
			{

				bool isEOF = false;
				Functor functorSub(functor);
				/*
				cout << "t" << endl;
				kseq_t *seq;

					#pragma omp critical
					{
				seq = kseq_init(fp);
					}


				int result = 1;
				*/
				Read read;

				while(true){


					#pragma omp critical
					{
						int result = kseq_read(seq);
						isEOF = result < 0;

						if(!isEOF){
							readIndex += 1;

							if(seq->qual.l == 0){
								read = {readIndex, string(seq->name.s), string(seq->seq.s)};
							}
							else{
								read = {readIndex, string(seq->name.s), string(seq->seq.s), string(seq->qual.s)};
							}
						}

						//cout << "allo" << endl;
						//cout << seq->name.s << endl;
						//cout << read._seq << endl;
						//cout << read._qual << endl;
						//cout << seq->name.s << endl;
						//cout << "1" << endl;

						//getchar();
						
						//cout << result << endl;
					}

					if(isEOF) break;
					functorSub(read);

				}
				

			}
			
			kseq_destroy(seq);	
			gzclose(fp);

		}

		cout << readIndex << endl;
	}
};

class ReadParser{
public:

	string _inputFilename;
	bool _isFile;
	bool _isBitset;
	size_t _l;
	size_t _k;
	float _density;
	vector<string> _filenames;
	u_int64_t _nbDatasets;
	u_int64_t _maxReads;

	ReadParser(const string& inputFilename, bool isFile, bool isBitset){
		_maxReads = 0;
		_inputFilename = inputFilename;
		_isFile = isFile;
		_isBitset = isBitset;
		
		if(_isFile){
			_nbDatasets = 1;
			_filenames.push_back(inputFilename);
		}
		else{
			parseFilenames();
		}
	}


	ReadParser(const string& inputFilename, bool isFile, size_t l, size_t k, float density){
		_maxReads = 0;
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
			if(line.empty()) continue;

			fs::path path(line);

			if(fs::exists (path)){
				_nbDatasets += 1;
				_filenames.push_back(line);
			}
			else{
				cout << "File not found: " << line << endl;
				exit(1);
			}
		}


	}



	void parse(const std::function<void(Read)>& fun){

		u_int64_t readIndex = 0;

		vector<Read> sequenceBuffer;

		for(const string& filename : _filenames){

			cout << filename << endl;

			/*
			if(_isBitset){

				ifstream fp(filename);

				kseq_t *read;
				read = new kseq_t();
				//read = kseq_init(fp);

				while (true) {

					u_int32_t sizeData = -1;
					fp.read((char*)&sizeData, sizeof(sizeData));

    				if(fp.eof())break;

					u_int32_t sizeSequence = -1;
					fp.read((char*)&sizeSequence, sizeof(sizeSequence));

					uint8_t* m_data = new uint8_t[sizeData];
					fp.read((char*)&m_data[0], sizeData*sizeof(uint8_t));

					DnaBitset* dnaBitset = new DnaBitset(m_data, sizeData, sizeSequence);

					char* seq = dnaBitset->to_string();
					read->seq.s = seq;
					
					fun(read, readIndex);

					free(seq);
					delete dnaBitset;

					readIndex += 1;
				}

				delete read;
				//kseq_destroy(read);
				fp.close();
			}
			else{
				*/
				gzFile fp;
				kseq_t *seq;
				int slen = 0, qlen = 0;
				fp = gzopen(filename.c_str(), "r");
				seq = kseq_init(fp);

				while (kseq_read(seq) >= 0){

					if(seq->qual.l == 0){
						fun({readIndex, string(seq->name.s), string(seq->seq.s)});
					}
					else{
						fun({readIndex, string(seq->name.s), string(seq->seq.s), string(seq->qual.s)});
					}
					readIndex += 1;

					if(_maxReads > 0 && readIndex > _maxReads) break;
					/*
					sequenceBuffer.push_back({readIndex, string(seq->name.s), string(seq->seq.s)});
					readIndex += 1;

					if(sequenceBuffer.size() > 100){
						for(const Read& read : sequenceBuffer){
							fun(read);
						}
						sequenceBuffer.clear();
					}

					//fun(seq, readIndex);
					*/
				}
					
				kseq_destroy(seq);
				gzclose(fp);
			//}

		}

	}
	


	void parseKminmers(const std::function<void(vector<KmerVec>, vector<ReadKminmer>, u_int64_t, u_int64_t, string, string)>& fun){

		MinimizerParser* _minimizerParser = new MinimizerParser(_l, _density);
		EncoderRLE encoderRLE;

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
				encoderRLE.execute(read->seq.s, strlen(read->seq.s), rleSequence, rlePositions);

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
		auto fp = std::bind(&ReadParser::extractSubsample_read, this, std::placeholders::_1);
		parse(fp);

		gzclose(_extractSubsample_outputFile);

	}

	gzFile _extractSubsample_outputFile;
	unordered_set<u_int64_t> _selectedReads;

	void extractSubsample_read(const Read& read){

		u_int64_t readIndex = read._index;

		if(_selectedReads.find(readIndex) != _selectedReads.end()){
			string header = ">" + read._header + "\n";
			string seq = read._seq + "\n";
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
	bool _usePos;
	bool _hasQuality;

	unordered_set<u_int64_t> _isReadProcessed;

	KminmerParser(){
	}

	KminmerParser(const string& inputFilename, size_t l, size_t k, bool usePos, bool hasQuality){
		_inputFilename = inputFilename;
		_l = l;
		_k = k;
		_usePos = usePos;
		_hasQuality = hasQuality;
	}

	void parse(const std::function<void(vector<u_int64_t>, vector<KmerVec>, vector<ReadKminmer>, bool, u_int64_t)>& fun){

		ifstream file_readData(_inputFilename, std::ios::binary);

		u_int64_t readIndex = 0;

		while(true){
			
			
			u_int32_t size;
			vector<u_int64_t> minimizers;
			vector<u_int16_t> minimizersPosOffsets; 
			vector<u_int8_t> minimizerQualities;
			
			file_readData.read((char*)&size, sizeof(size));

			if(file_readData.eof())break;

			minimizers.resize(size);
			minimizersPosOffsets.resize(size);
			minimizerQualities.resize(size, 0);


			bool isCircular;
			file_readData.read((char*)&isCircular, sizeof(isCircular));

			//cout << size << " " << isCircular << endl;
			file_readData.read((char*)&minimizers[0], size*sizeof(u_int64_t));
			if(_usePos) file_readData.read((char*)&minimizersPosOffsets[0], size*sizeof(u_int16_t));
			if(_hasQuality) file_readData.read((char*)&minimizerQualities[0], size*sizeof(u_int8_t));

			//if(_isReadProcessed.size() > 0 && _isReadProcessed.find(readIndex) != _isReadProcessed.end()){
			//	readIndex += 1;
			//	continue;
			//}

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

			fun(minimizers, kminmers, kminmersInfo, isCircular, readIndex);

			readIndex += 1;

		}

		file_readData.close();

	}

	
	//template<typename SizeType>
	void parseMinspace(const std::function<void(vector<u_int64_t>, vector<ReadKminmerComplete>, bool, u_int64_t)>& fun){

		ifstream file_readData(_inputFilename, std::ios::binary);

		u_int64_t readIndex = 0;

		while(true){
			
			u_int32_t size;
			vector<u_int64_t> minimizers;
			vector<u_int16_t> minimizersPosOffsets; 
			vector<u_int8_t> minimizerQualities;
			
			file_readData.read((char*)&size, sizeof(size));

			if(file_readData.eof())break;

			minimizers.resize(size);
			minimizersPosOffsets.resize(size);
			minimizerQualities.resize(size, 0);


			bool isCircular;
			file_readData.read((char*)&isCircular, sizeof(isCircular));

			//cout << size << " " << isCircular << endl;

			file_readData.read((char*)&minimizers[0], size*sizeof(u_int64_t));
			if(_usePos){
				file_readData.read((char*)&minimizersPosOffsets[0], size*sizeof(u_int16_t));
			}
			if(_hasQuality) file_readData.read((char*)&minimizerQualities[0], size*sizeof(u_int8_t));
			
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
			MDBG::getKminmers_complete(_k, minimizers, minimizersPos, kminmersInfo, readIndex, minimizerQualities);

			fun(minimizers, kminmersInfo, isCircular, readIndex);

			readIndex += 1;
		}

		file_readData.close();

	}
	
	void parseSequences(const std::function<void(vector<u_int64_t>, bool, u_int64_t)>& fun){

		ifstream file_readData(_inputFilename, std::ios::binary);

		u_int64_t readIndex = 0;

		while(true){
			
			u_int32_t size;
			vector<u_int64_t> minimizers;
			vector<u_int16_t> minimizersPosOffsets; 
			vector<u_int8_t> minimizerQualities;
			
			file_readData.read((char*)&size, sizeof(size));

			if(file_readData.eof())break;

			minimizers.resize(size);
			minimizersPosOffsets.resize(size);
			minimizerQualities.resize(size, 0);


			bool isCircular;
			file_readData.read((char*)&isCircular, sizeof(isCircular));


			file_readData.read((char*)&minimizers[0], size*sizeof(u_int64_t));
			if(_usePos){
				file_readData.read((char*)&minimizersPosOffsets[0], size*sizeof(u_int16_t));
			}
			if(_hasQuality) file_readData.read((char*)&minimizerQualities[0], size*sizeof(u_int8_t));
			
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

			fun(minimizers, isCircular, readIndex);

			readIndex += 1;
		}

		file_readData.close();

	}

	/*
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
	*/
	/*
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
	*/

};




class KminmerParserParallel{

public:

	string _inputFilename;
	size_t _l;
	size_t _k;
	bool _usePos;
	int _nbCores;
	bool _hasQuality;

	unordered_set<u_int64_t> _isReadProcessed;

	KminmerParserParallel(){
	}

	KminmerParserParallel(const string& inputFilename, size_t l, size_t k, bool usePos,bool hasQuality, int nbCores){

		if(!fs::exists(inputFilename)){
			cout << "File not found: " << inputFilename << endl;
			exit(1);
		}

		_inputFilename = inputFilename;
		_l = l;
		_k = k;
		_usePos = usePos;
		_hasQuality = hasQuality;
		_nbCores = nbCores;
	}

	template<typename Functor>
	void parse(const Functor& functor){
	//void parse(const std::function<void(vector<u_int64_t>, vector<KmerVec>, vector<ReadKminmer>, u_int64_t)>& fun){

		ifstream file_readData(_inputFilename);

		u_int64_t readIndex = -1;

		#pragma omp parallel num_threads(_nbCores)
		{

			bool isEOF = false;
			Functor functorSub(functor);
			vector<u_int64_t> minimizers;
			vector<u_int16_t> minimizersPosOffsets; 
			vector<u_int8_t> minimizerQualities; 
			u_int32_t size;
			KminmerList kminmerList;
			//KminmerList kminmer;

			while(true){
				

				#pragma omp critical
				{

					
					file_readData.read((char*)&size, sizeof(size));

					if(file_readData.eof()) isEOF = true;

					if(!isEOF){

						readIndex += 1;

						kminmerList = {readIndex};

						minimizers.resize(size);
						minimizersPosOffsets.resize(size);
						minimizerQualities.resize(size, -1);

						
						bool isCircular;
						file_readData.read((char*)&isCircular, sizeof(isCircular));

						file_readData.read((char*)&minimizers[0], size*sizeof(u_int64_t));
						//if(_usePos) file_readData.read((char*)&minimizersPosOffsets[0], size*sizeof(u_int16_t));
						if(_hasQuality) file_readData.read((char*)&minimizerQualities[0], size*sizeof(u_int8_t));
					}

				}

				

				//if(_isReadProcessed.size() > 0 && _isReadProcessed.find(readIndex) != _isReadProcessed.end()){
				//	readIndex += 1;
				//	continue;
				//}

				//cout << "----" << endl;

				if(isEOF) break;

				
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
				

				//vector<KmerVec> kminmers; 
				//vector<ReadKminmer> kminmersInfo;
				vector<u_int64_t> rlePositions;
				vector<ReadKminmerComplete> kminmersInfo;
				//MDBG::getKminmers(_l, _k, minimizers, minimizersPos, kminmers, kminmersInfo, rlePositions, 0, false);
				MDBG::getKminmers_complete(_k, minimizers, minimizersPos, kminmersInfo, readIndex, minimizerQualities);
				
				//fun(minimizers, kminmers, kminmersInfo, readIndex);
				kminmerList._readMinimizers = minimizers;
				//kminmerList._kminmers = kminmers;
				kminmerList._kminmersInfo = kminmersInfo;
				functorSub(kminmerList);
			}
		}

		file_readData.close();

	}

	/*
	template<typename Functor>
	void parse(const Functor& functor){

		//#pragma omp for
		cout << _nbCores << endl;

		u_int64_t readIndex = -1;

		for(const string& filename : _filenames){

			cout << filename << endl;

			gzFile fp = gzopen(filename.c_str(), "r");
			kseq_t *seq;
			seq = kseq_init(fp);

			#pragma omp parallel num_threads(_nbCores)
			{

				bool isEOF = false;
				Functor functorSub(functor);

				Read read;

				while(true){


					#pragma omp critical
					{
						int result = kseq_read(seq);
						readIndex += 1;
						isEOF = result < 0;

						read = {readIndex, string(seq->name.s), string(seq->seq.s)};
						//cout << seq->name.s << endl;
						//cout << "1" << endl;

						//getchar();
						
						//cout << result << endl;
					}

					if(isEOF) break;
					functorSub(read);

				}
				

			}
			
			kseq_destroy(seq);	
			gzclose(fp);

		}

		cout << readIndex << endl;
	}
	*/

	/*
	//template<typename SizeType>
	void parseMinspace(const std::function<void(vector<u_int64_t>, vector<ReadKminmerComplete>, u_int64_t)>& fun){

		ifstream file_readData(_inputFilename, std::ios::binary);

		u_int64_t readIndex = 0;

		while(true){
			
			u_int32_t size;
			vector<u_int64_t> minimizers;
			vector<u_int16_t> minimizersPosOffsets; 
			
			file_readData.read((char*)&size, sizeof(size));

			if(file_readData.eof())break;

			minimizers.resize(size);
			minimizersPosOffsets.resize(size);

			file_readData.read((char*)&minimizers[0], size*sizeof(u_int64_t));
			if(_usePos){
				file_readData.read((char*)&minimizersPosOffsets[0], size*sizeof(u_int16_t));
			}

			
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

		file_readData.close();

	}
	*/

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