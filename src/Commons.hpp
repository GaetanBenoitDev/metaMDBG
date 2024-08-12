

/*

- segfault dans final n50
- COntigPolisher: autoriser plus de mémoire, en fonction du maximum utilisé jusqu'a là
- essayer d'enlever les contigs avec une faible abondance calculer par semibon2, et refaire le binning
- dans la table finale, enlever directement les single contig complete des complete MAGs

- humanGut: subsampler par pas de 10G x 5

- build: penser a enlever le first graph (dans pipeline et cleaning)

- dorado correct: besoin de bgzip (on aurait du tout comrpesser en bgzip directement)
- meilleur progress dans correction (genere trop de ligne la)

- option: --high-quality-correction ?

- maxMappedRead: trouver un truc malin
	- une fois qu'on a atteint une certaine quantité de mapping, on doit pouvoir jauger si un nouvel alignment va apporter qqch, en gros check s'il peut theoriquement augmenter le score d'alignment ou si c'est juste impossible en fonction du nb de matche etc

- "high density chaining: faire un petit index de minimizer pour eviter de tout query? pour les read chargé"

- chaining:
	- trouver la bonne valeur pour la banding

- truc fait a tester:
	- estimate_alignmentCOverage (juste decommenté)


*/


/*

- papier: refaire une simulation plus grosse avec une dizaine de ecoli.


- isRepeatSide: tester > _kminmerSize*2 au lieu de 50k length
- last superbubble in circular component not removed (cycle including source)

- unitig graph:
	- memory a optimiser (le graphSUccessor et graphSImplify peuvent etre detruit une fois que le unitg graph est construit)
	

- truc desactivé:
	- GraphSimplify: auto currentAbundance 


- contig polisher: plus de window pour les espece plus abondante ? (fraction du contig coverage)
- OverlapRemover : ligne 595: tester d'enlever le +1 de if(contig._minimizers.size() <= _kminmerSize+1){
- determinstic node a un imapct sur la qualité des resultats sur donées simulé, voir si c'est le cas sur AD

- Dans les première phase du multi, ne pas ollapse les bubbles si leur source et sink sont très abondante (on ne mergera donc pas deux organismes, et il ne devrait pas break vu que c'est abondante)
- polisher circ detection incomplete (il faudrait ajouter la partie manquante au backbone)
- polisher circ detection: update pour gerer les petit contig sur lequel des read entier pourrait mapper sur toute leur longeur

- check whitin contig contamination

- tout ce qui est écrit dans le cerr doit aussi etre ecrit dans les logs sinon c'est illisible

Paralelisation:
	- OverlapRemover: " parallelisation ligne 475 for(long i=0; i<_contigs.size(); i++){ en lisant sequenciellement avec une result queue comme dans readSelection"
		- prototype dans backup mais fonctionne pas :/

ContigPolisher:
	- read mal mapper sur les bord des contig circulaire (double mapping)

- ajouter une progress bar global

---------------------------------- FIXED

- [Fixed] CreateGraph: determinstic sutff always enabled, remove for speed gain

- [Fixed temporaire] pb avec le abundance filter: les repeat se font output en premier et peuvent créé des faux unitig circulaires
	- on pourrait empecher les petit unitig d'etre circulaire dans les premiere passes, sinon il faut reussir a savoir si se sont des repeats

- multi-k: en fait le +1 est nécessaire que au debut probablement, on pourrait faire un plus grand pas dans les grande valeur de k
	- update:on trouve moins de component circulaire 

*/


#ifndef MDBG_METAG_COMMONS
#define MDBG_METAG_COMMONS

//#include <gatb/gatb_core.hpp>
#include "utils/MurmurHash3.h"
#include <zlib.h>
#include <string>
#include <sstream>
#include <unordered_set>
#include <unordered_map>
//#include <regex>
#include <algorithm>
#include <libgen.h>
#include <set>
#include <list>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <cmath>
//#include "./utils/ArgParse.hpp"
#include "./utils/args.hxx"
#include <random>
#include <chrono>
#include "utils/parallel_hashmap/phmap.h"
#include "utils/kmer/Kmer.hpp"
#include <omp.h>
#include <queue>
#include "utils/PafParser.hpp"
#include "readSelection/MinimizerAligner.hpp"
#include "utils/Logger.h"
#include <exception>

namespace fs = std::filesystem;
using namespace std::chrono;
//using std::chrono::high_resolution_clock;
//using std::chrono::duration_cast;
//using std::chrono::duration;
//using std::chrono::milliseconds;
//namespace fs = std::filesystem;

typedef unsigned __int128 u_int128_t;
typedef u_int32_t ReadType;


//#include <sys/types.h>
//#include <sys/stat.h>

#include <functional>
#include <memory>
//#include "./utils/ntHashIterator.hpp"
#include "./utils/kseq.h"
#include "./utils/DnaBitset.hpp"
#include "./utils/BloomFilter.hpp"

KSEQ_INIT(gzFile, gzread)

using namespace std;

struct Read{
	u_int64_t _index;
	string _header;
	string _seq;
	string _qual;
	u_int64_t _datasetIndex;
};

struct ReadMatch{
	u_int32_t _readIndex;
	bool _isReversed;
	u_int32_t _referenceStart;
	u_int32_t _referenceEnd;
	u_int32_t _queryStart;
	u_int32_t _queryEnd;
	u_int32_t _hangLeft;
	u_int32_t _hangRight;
	u_int32_t _alignNbMatches;
	u_int32_t _alignLength;
	float _divergence;

};


struct ReadMatchBound{
	ReadType _queryReadIndex;
	u_int32_t _nbMatches;
	u_int16_t _minIndex;
	u_int16_t _maxIndex;

	ReadMatchBound(){
		_queryReadIndex = -1;
		_nbMatches = 0;
		_minIndex = -1;
		_maxIndex = 0;
	}
};


struct KminmerSequenceVariant{
	u_int16_t _editDistance;
	DnaBitset* _sequence;
};

struct KminmerSequenceVariant_Comparator {
	bool operator()(KminmerSequenceVariant const& p1, KminmerSequenceVariant const& p2)
	{
		return p1._editDistance < p2._editDistance;
	}
};

typedef priority_queue<KminmerSequenceVariant, vector<KminmerSequenceVariant>, KminmerSequenceVariant_Comparator> VariantQueue;

struct ReadSequence{
	u_int64_t _readIndex;
	DnaBitset* _sequence;
	VariantQueue _variants;
};

/*
struct AlignmentResult{

	public:

	u_int64_t _readIndex;
	int64_t _nbMatches;
	int64_t _nbMissmatches;
	int64_t _nbInsertions;
	int64_t _nbDeletions;
	int64_t _alignLengthBps;

	int64_t score() const{
		int64_t nbErrors = _nbMissmatches + _nbInsertions + _nbDeletions;
		return _nbMatches - nbErrors;
	}

	int64_t nbErrors() const{
		return _nbMissmatches + _nbInsertions + _nbDeletions;
	}

	double divergence() const{

		double nbSeeds = _nbMatches + _nbMissmatches + _nbInsertions;// + _nbDeletions;
		
		if(_nbMatches == nbSeeds) return 0;
		if(_nbMatches == 0) return 1;
		//(95/100)^(1/13)
		//if(_nbMatches == 0) return 1;
		return 1.0 - pow((_nbMatches / nbSeeds), 1.0/13.0);
	}
	
};
*/


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


const string ARG_HOMOPOLYMER_COMPRESSION = "homopolymer-compression";
//const string ARG_INPUT_FILENAME = "i";
//const string ARG_INPUT_FILENAME_TRUTH = "itruth";
//const string ARG_OUTPUT_DIR = "o";
const string ARG_OUTPUT_FILENAME = "f";
//const string ARG_MINIMIZER_LENGTH = "k";
//const string ARG_KMINMER_LENGTH = "k";
//const string ARG_MINIMIZER_DENSITY = "d";
const string ARG_DEBUG = "debug";
//const string ARG_FINAL = "final";
const string ARG_INPUT_FILENAME_CONTIG = "c";
const string ARG_INPUT_FILENAME_CONTIG_FASTA = "cf";
//const string ARG_INPUT_FILENAME_UNITIG_NT = "unitigNt";
//const string ARG_INPUT_FILENAME_UNITIG_CLUSTER = "cluster";
//const string ARG_INPUT_FILENAME_ABUNDANCE = "a";
//const string ARG_INPUT_FILENAME_BINNING = "bi";
//const string ARG_OUTPUT_FILENAME_BINNING = "bo";
const string ARG_FIRST_PASS = "firstpass";
//const string ARG_FASTA = "fasta";
const string ARG_NB_CORES = "t";
//const string ARG_EVAL = "eval";
const string ARG_BLOOM_FILTER = "nofilter";		
const char ARG_MIN_IDENTITY = 'i';
//const char ARG_CIRCULAR_LENGTH = 'c';
//const char ARG_LINEAR_LENGTH = 'l';
const char ARG_NB_WINDOWS = 'n';
const string ARG_MIN_READ_QUALITY = "min-read-quality";
//const string ARG_HOMOPOLYMER_COMPRESSION = "hpc";		
//const string ARG_CORRECTION = "correction";	

const string NB_CORES_DEFAULT = "3";
const int NB_CORES_DEFAULT_INT = 1;
//const string FILENAME_NO_KMINMER_READS = "reads_noKminmers.bin";


const u_int8_t CONTIG_LINEAR = 0; 
const u_int8_t CONTIG_CIRCULAR = 1;
const u_int8_t CONTIG_CIRCULAR_RESCUED = 2;

//const char ARG_INPUT_FILENAME2 = 'i';
const string ARG_INPUT_HIFI = "in-hifi";
const string ARG_INPUT_NANOPORE = "in-ont";
const string ARG_OUTPUT_DIR2 = "out-dir";
const char ARG_OUTPUT_FILENAME2 = 'f';
const string ARG_MINIMIZER_LENGTH2 = "kmer-size";
//const char ARG_KMINMER_LENGTH2 = 'k';
const string ARG_MINIMIZER_DENSITY_ASSEMBLY = "density-assembly";
const string ARG_MINIMIZER_DENSITY_CORRECTION = "density-correction";
const string ARG_MAXK = "max-k";
//const char ARG_INPUT_FILENAME_CONTIG2 = 'c';
//const char ARG_INPUT_FILENAME_ABUNDANCE2 = 'a';
const string ARG_NB_CORES2 = "threads";
//const string ARG_CIRCULARIZE = "circ";
//const int MIN_NB_MINIMIZER_REFERENCE = 4000;
//const int MIN_NB_MINIMIZER_QUERY = 4000;

const vector<string> possibleInputArguments = {"--" + ARG_INPUT_HIFI, "--" + ARG_INPUT_NANOPORE};

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

struct MinimizerRead{

	public:

	ReadType _readIndex;
	vector<MinimizerType> _minimizers;
	vector<u_int32_t> _minimizersPos;
	vector<u_int8_t> _qualities;
	vector<u_int8_t> _readMinimizerDirections;
	//float _meanReadQuality;
	u_int32_t _readLength;
	//float _debugNbMatches;

	void print() const{
		for(size_t i=0; i<_minimizers.size(); i++){
			cout << i << "\t" << _minimizers[i] << "\t" << _minimizersPos[i] << "\t" << (int) _qualities[i] << "\t" << (int) _readMinimizerDirections[i] << endl;
		}
	}


	/*
	void compare(const MinimizerRead& otherRead){
		if(_readIndex != otherRead._readIndex){
			cout << "pb 1" << endl;
			exit(1);
		}
		if(_minimizers.size() != otherRead._minimizers.size()){
			cout << "pb 2" << endl;
			exit(1);
		}
		if(_minimizersPos.size() != otherRead._minimizersPos.size()){
			cout << "pb 3" << endl;
			exit(1);
		}
		if(_qualities.size() != otherRead._qualities.size()){
			cout << "pb 4" << endl;
			exit(1);
		}
		if(_readMinimizerDirections.size() != otherRead._readMinimizerDirections.size()){
			cout << "pb 5" << endl;
			exit(1);
		}
		for(size_t i=0; i<_minimizers.size(); i++){
			if(_minimizers[i] != otherRead._minimizers[i]){
				cout << "pb 6" << endl;
				exit(1);
			}
			if(_minimizersPos[i] != otherRead._minimizersPos[i]){
				cout << otherRead._readIndex << endl;
				cout << i << "\t" << _minimizersPos[i] << "\t" << otherRead._minimizersPos[i] << endl;
				cout << "pb 7" << endl;
				exit(1);
			}
			if(_qualities[i] != otherRead._qualities[i]){
				cout << "pb 8" << endl;
				exit(1);
			}
			if(_readMinimizerDirections[i] != otherRead._readMinimizerDirections[i]){
				cout << "pb 9" << endl;
				exit(1);
			}
		}
	}
	*/
};

/*
struct MinimizerReadBinary{
	ReadType _readIndex;
	vector<u_int64_t> _minimizers;
	//float _debugNbMatches;
};

struct MinimizerRead{
	u_int64_t _readIndex;
	vector<MinimizerType> _minimizers;
	vector<u_int32_t> _minimizersPos;
	vector<u_int8_t> _qualities;
	vector<u_int8_t> _readMinimizerDirections;
	float _debugNbMatches;
	u_int32_t _hangLeft;
	u_int32_t _hangRight;
};
*/
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


//typedef phmap::parallel_flat_hash_map<MinimizerType, vector<u_int32_t>, phmap::priv::hash_default_hash<MinimizerType>, phmap::priv::hash_default_eq<MinimizerType>, std::allocator<std::pair<MinimizerType, vector<u_int32_t>>>, 4, std::mutex> MinimizerReadMap;





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

/*
struct UnitigEdgeScore{
	u_int32_t _from;
	u_int32_t _to;
	u_int16_t _score;
};

struct MinimizerPair{
	MinimizerType _first;
	MinimizerType _second;

	bool operator==(const MinimizerPair &other) const{
		return _first == other._first && _second == other._second;
	}
};

*/
struct DbgNode{
	u_int32_t _index;
	u_int32_t _abundance;
	//u_int32_t _quality;
	//u_int16_t _length;
	//u_int16_t _overlapLength_start;
	//u_int16_t _overlapLength_end;
	//bool _isReversed;
	//u_int32_t _unitigNbNodes;
};



//struct ReadData{
//	u_int32_t _length;
//	vector<float> _composition;
//};

struct KmerVec{
	vector<MinimizerType> _kmers;

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
		for(MinimizerType m : _kmers){
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

//struct MinimizerPair_Edge{
//	MinimizerPair _from;
//	MinimizerPair _to;
//};


struct ContigPosition{
	u_int32_t _contigIndex;
	u_int32_t _contigPosition;

	bool operator==(const ContigPosition &other) const{
		return _contigIndex == other._contigIndex && _contigPosition == other._contigPosition;
	}

};

struct ReadPosition{
	u_int32_t _readIndex;
	u_int32_t _readPosition;

	bool operator==(const ReadPosition &other) const{
		return _readIndex == other._readIndex && _readPosition == other._readPosition;
	}

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
	//template <>
	//struct hash<MinimizerPair>{
	//	std::size_t operator()(const MinimizerPair& k) const{
	//		using std::size_t;
	//		using std::hash;
	//		using std::string;

	//		return ((hash<MinimizerType>()(k._first) ^ (hash<MinimizerType>()(k._second) << 1)) >> 1);
	//	}
	//};


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

	template <>
	struct hash<ContigPosition>{
		std::size_t operator()(const ContigPosition& k) const{
			//using std::size_t;
			//using std::hash;
			//using std::string;

			std::size_t seed = 2;
			seed ^= k._contigIndex + 0x9e3779b9 + (seed << 6) + (seed >> 2);
			seed ^= k._contigPosition + 0x9e3779b9 + (seed << 6) + (seed >> 2);
			//auto hash1 = hash<T1>{}(p._from);
			//auto hash2 = hash<T2>{}(p._to);
			return seed;

			//return ((hash<u_int64_t>()(k._first) ^ (hash<u_int64_t>()(k._second) << 1)) >> 1);
		}
	};

	template <>
	struct hash<ReadPosition>{
		std::size_t operator()(const ReadPosition& k) const{
			//using std::size_t;
			//using std::hash;
			//using std::string;

			std::size_t seed = 2;
			seed ^= k._readIndex + 0x9e3779b9 + (seed << 6) + (seed >> 2);
			seed ^= k._readPosition + 0x9e3779b9 + (seed << 6) + (seed >> 2);
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
struct Overlap {
    ReadIndexType _r1;
    ReadIndexType _r2;
    u_int16_t _nbMinimizers;
};

*/
struct KminmerEdge2{
	u_int32_t _nodeName;
	MinimizerType _minimizer;
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
//typedef phmap::parallel_flat_hash_set<KmerVec, phmap::priv::hash_default_hash<KmerVec>, phmap::priv::hash_default_eq<KmerVec>, std::allocator<KmerVec>, 4, std::mutex> KminmerSet;


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

/*
struct KminmerListPaired{
	u_int64_t _readIndex;
	vector<MinimizerType> _readMinimizers1;
	vector<MinimizerType> _readMinimizers2;
	vector<u_int32_t> _readMinimizers1pos;
};
*/


struct NodePath{
	u_int64_t _readIndex;
	vector<u_int32_t> _nodePath;
	bool _isCircular;
	//vector<KmerVec> _kminmers;
	//vector<ReadKminmer> _kminmersInfo;
};

class Commons{

public:

	static string inputFileToFilenames(const string& inputFilename){

		string filenames = " ";
		ifstream inputFile(inputFilename);
		std::string line;

		while (std::getline(inputFile, line)){

			filenames += line + " ";
		}

		return filenames;

	}

	static void createInputFile(auto& paths, const string& filename){


		//ofstream inputFileHiFi(filename + ".hifi");
		//ofstream inputFileNanopore(filename + ".ont");

		//string inputFilenames = _inputFilename;
		ofstream inputFile(filename);

		int  i = 0;
		//cout << inputFilenames << endl;
		//vector<string>* fields = new vector<string>();
		//GfaParser::tokenize(inputFilenames, fields, ',');

		for(auto &&path : paths){

			string filename = path;


    		//filename.erase(std::remove_if(filename.begin(), filename.end(), ::isspace), filename.end());
			//if(filename.empty()) continue;
			
			if(!fs::exists (filename)){
				Logger::get().error() << "File not found: " << filename;
				exit(1);
			}

			if(fs::is_directory(filename)){
				Logger::get().error() << "Input file (" << filename << ") is a directory (need fasta/fastq filenames)";
				exit(1);
			}


			filename = fs::canonical(fs::absolute(filename));
			//else{
			//	cerr << "File not found: " << line << endl;
			//	exit(1);
			//}

			//string filenameAbs = fs::absolute(filename);
			//cout << path << endl; //" " << fs::canonical(filename) << " " << fs::weakly_canonical(filename) << endl;
			inputFile << filename << endl;

			//if(i == 0){
			//	inputFileHiFi << filename << endl;
			//}
			//else if(i == 1){
			//	inputFileNanopore << filename << endl;
			//}

			i += 1;
		}

		//delete fields;c
		inputFile.close();
		//inputFileHiFi.close();
		//inputFileNanopore.close();


	}

	static size_t getMultikStep(size_t k){
		return 1;

		if(k < 20){
			return 1;
		}
		else if(k < 40){
			return 2;
		}
		else{
			return 5;
		}
	}

};


class Utils{




public:

	//static bool isReadTooShort(const MinimizerRead& minimizerReadLowDensity){
	//	if(minimizerReadLowDensity._minimizers.size() < 10) return true;
	//	if(minimizerReadLowDensity._readLength < 500) return true;

	//	return false;
	//}

	static string getProgress(const ReadType& nbReadProcessed, const ReadType& totalNbReads){
		int progess = ((((long double) nbReadProcessed) / ((long double) totalNbReads)) * 100);

   		//progess = (int)( 100 * progess )  /  100.0;
		//progess *= 100;
		return to_string(progess) + "%";
	}
	static string createContigHeader(const u_int32_t& contigIndex, const u_int32_t& length, const u_int32_t& coverage, bool isCircular){

		string circularStr = "no";
		if(isCircular) circularStr = "yes";

		return "ctg" + to_string(contigIndex) + " length=" + to_string(length) + " coverage=" + to_string(coverage) + " circular=" + circularStr;
	}

	struct ContigHeader{
		u_int32_t _contigIndex;
		u_int32_t _length;
		u_int32_t _coverage;
		bool _isCircular;
	};

	static ContigHeader extractContigHeader(const string& header){
		vector<string> fields = Utils::split(header, ' ');

		ContigHeader contigHeader;

		string contigIndexField = fields[0];
		contigIndexField.erase(0,3); //"remove "ctg" letters

		contigHeader._contigIndex = stoull(contigIndexField);
		contigHeader._length = stoull(Utils::split(fields[1], '=')[1]);
		contigHeader._coverage = stoull(Utils::split(fields[2], '=')[1]);
		contigHeader._isCircular = Utils::split(fields[3], '=')[1] == "yes";

		return contigHeader;
	}

	static bool isContigCircular(const string& header){

		const ContigHeader& contigHeader = Utils::extractContigHeader(header);
		return contigHeader._isCircular;
		//string circularStr = Utils::split(contigHeader, ' ')[3];
		
		//cout << circularStr << endl;
		//cout << Utils::split(circularStr, '=')[1] << endl;
		//getchar();
		//return Utils::split(circularStr, '=')[1] == "yes";
	}

	static vector<string> split(const string& str, const char& delimiter){
		vector<string> fields;
		
		istringstream iss(str);
		string s;
		while ( getline( iss, s, delimiter) ) {
			fields.push_back(s);
		}

		return fields;
	}
	static u_int64_t computeN50(vector<u_int32_t> allReadLengths){

		if(allReadLengths.size() == 0) return 0;

		std::sort(allReadLengths.begin(),  allReadLengths.end(), std::greater<u_int32_t>());

		//u_int32_t n=_allReadSizes.size();
		//u_int32_t max=_allReadSizes[0];                 	
		//u_int32_t  sum = accumulate(_allReadSizes.begin(), _allReadSizes.end(), 0.0);
		vector<u_int64_t> readLengthCumuls;
		u_int64_t cumul = 0;

		for(size_t i=0; i<allReadLengths.size(); i++){
			cumul += allReadLengths[i];
			readLengthCumuls.push_back(cumul);
		}

		std::reverse(allReadLengths.begin(), allReadLengths.end());
		std::reverse(readLengthCumuls.begin(), readLengthCumuls.end());

		u_int32_t n50 = allReadLengths[allReadLengths.size()-1];
		u_int64_t halfsize = readLengthCumuls[0]/2;

		for(size_t i=0; i<allReadLengths.size(); i++){
			if(readLengthCumuls[i] < halfsize){
				n50 = allReadLengths[i];
				break;
			}
		}

		return n50;
	}

	static float transformQuality(u_int8_t quality){
		float q = quality;
		return pow(10.0f, -q/10.0f);
	}

	/*
	enum CigarType{
		Undefined,
		Match,
		Insertion,
		Deletion,
	};

	struct CigarElement{
		u_int32_t _nbOccurences;
		CigarType _cigarType;
	};



	static vector<CigarElement> extractCigarSequence(const string& cigar){

		//cout << cigar << endl;
		vector<CigarElement> cigarSequence;

		string numberStr = "";

		for(size_t i=0; i<cigar.size(); i++){

			char c = cigar[i];

			if(isdigit(c)){
				numberStr += c;
			}
			else{

				u_int32_t nbOccurences = stoull(numberStr);
				numberStr = "";

				CigarType cigarType = CigarType::Undefined;

				if(c == 'M'){
					cigarType = CigarType::Match;
				}
				else if(c == 'I'){
					cigarType = CigarType::Insertion;
				}
				else if(c == 'D'){
					cigarType = CigarType::Deletion;
				}

				CigarElement cigarElement = {nbOccurences, cigarType};
				cigarSequence.push_back(cigarElement);
			}
		}

		return cigarSequence;
	}

	const static u_int16_t maxInt16Value = -1;
	
	static bool isValidPositions(const vector<u_int32_t>& minimizerPositions){

		u_int64_t prevMinimizerPos = 0;

		for(size_t i=0; i<minimizerPositions.size(); i++){

			if(minimizerPositions[i] - prevMinimizerPos >= Utils::maxInt16Value){
				return false;
			}

			prevMinimizerPos = minimizerPositions[i];
		}

		return true;
	}

	static MinimizerReadBinary compressMinimizerRead(const MinimizerRead& read){

		//bool isValid = true;

		vector<u_int64_t> minimizerDatas(read._minimizers.size(), 0);
		u_int16_t prevMinimizerPos = 0;

		for(size_t i=0; i<read._minimizers.size(); i++){

			MinimizerType minimizer = read._minimizers[i];

			//if(read._minimizersPos[i] - prevMinimizerPos >= maxInt16Value){
			//	isValid = false;
			//	break;
			//}

			u_int16_t minimizerPosDiff = read._minimizersPos[i] - prevMinimizerPos;
			u_int8_t quality = read._qualities[i];
			u_int8_t direction = read._readMinimizerDirections[i];

			//cout << "loulou: " << minimizer << " " << minimizerPosDiff << endl;
			u_int64_t minimizerData = 0;
			//minimizerData |= minimizerData & 0xffffffff00000000UL;
			minimizerData += (((u_int64_t) minimizer) << 32);// | (((u_int16_t) minimizerPosDiff) << 16) | (((u_int8_t) minimizer) << 8) | (((u_int8_t) direction));
			minimizerData += (((u_int64_t) minimizerPosDiff) << 16);
			minimizerData += (((u_int64_t) quality) << 8);
			//cout << (int) quality << " " << ((int) quality << 8) << " " << minimizerData << endl;
			minimizerData += (((u_int64_t) direction));

			minimizerDatas[i] = minimizerData;
			//cout << minimizerData << endl;

			prevMinimizerPos = read._minimizersPos[i];
		}

		//if(!isValid){
		//	return {read._readIndex, {}};
		//}

		return {read._readIndex, minimizerDatas};
	}


	static vector<MinimizerType> decompressMinimizerRead_minimizerOnly(const MinimizerReadBinary& read){

		vector<MinimizerType> minimizers(read._minimizers.size(), 0);

		for(size_t i=0; i<read._minimizers.size(); i++){
			u_int64_t minimizerData = read._minimizers[i];
			minimizers[i] = (u_int64_t) (minimizerData >> 32);
		}

		return minimizers;
	}


	static MinimizerRead decompressMinimizerRead(const MinimizerReadBinary& read){

		vector<MinimizerType> minimizers(read._minimizers.size(), 0);
		vector<u_int32_t> minimizersPos(read._minimizers.size()+1, 0);
		vector<u_int8_t> qualities(read._minimizers.size(), 0);
		vector<u_int8_t> readMinimizerDirections(read._minimizers.size(), 0);

		//vector<u_int64_t> minimizerDatas(read._minimizers.size(), 0);
		u_int32_t prevMinimizerPos = 0;


		for(size_t i=0; i<read._minimizers.size(); i++){

			u_int64_t minimizerData = read._minimizers[i];
			//cout << minimizerData << endl;

			MinimizerType minimizer = (u_int64_t) (minimizerData >> 32);//& 0xffffffff00000000UL);
			u_int16_t minimizerPosDiff = (u_int64_t) ((minimizerData & 0xffffffffUL) >> 16);
			u_int8_t quality = (u_int64_t) ((minimizerData & 0xffffUL) >> 8);
			u_int8_t direction = (u_int64_t) ((minimizerData & 0xffUL));
			//cout << (int) quality << " " << (minimizerData & 0xffffUL) << " " << ((int) (minimizerData & 0xffffUL) >> 8) << endl;
			//cout << "lala " << minimizer << " " << minimizerPosDiff << endl;
			minimizers[i] = minimizer;
			minimizersPos[i] = minimizerPosDiff+prevMinimizerPos;
			qualities[i] = quality;
			readMinimizerDirections[i] = direction;

			prevMinimizerPos = minimizerPosDiff+prevMinimizerPos;
		}

		minimizersPos[minimizersPos.size()-1] = 0;

		return {read._readIndex, minimizers, minimizersPos, qualities, readMinimizerDirections};
	}
	*/

	static void applyDensityThreshold(const float densityThreshold, const vector<MinimizerType>& minimizers, const vector<u_int32_t>& minimizerPos, const vector<u_int8_t>& minimizerDirections, const vector<u_int8_t>& minimizerQualities, vector<MinimizerType>& minimizersFiltered, vector<u_int32_t>& minimizerPosFiltered, vector<u_int8_t>& minimizerDirectionsFiltered, vector<u_int8_t>& minimizerQualitiesFiltered){

		minimizersFiltered.clear();
		minimizerPosFiltered.clear();
		minimizerDirectionsFiltered.clear();
		minimizerQualitiesFiltered.clear();

		//if(_densityThreshold != -1){

		//u_int32_t maxHashValue_int32 = -1;
		//double minimizerBound_int32 = densityThreshold * maxHashValue_int32;
		//_minimizerBound_int32 = minimizerDensity * maxHashValue_int32;

		MinimizerType maxHashValue = -1;
		double minimizerBound = densityThreshold * maxHashValue;

		for(size_t i=0; i<minimizers.size(); i++){
			MinimizerType m = minimizers[i];

			u_int32_t minimizer_int32 = m;
			if(m < minimizerBound){
			//if(minimizer_int32 < minimizerBound_int32){
				minimizersFiltered.push_back(minimizers[i]);
				minimizerPosFiltered.push_back(minimizerPos[i]);
				minimizerDirectionsFiltered.push_back(minimizerDirections[i]);
				minimizerQualitiesFiltered.push_back(minimizerQualities[i]);
			}


		}

		//u_int32_t readLength = minimizerPos[minimizerPos.size()-1];
		//minimizerPosFiltered.push_back(readLength);
		//minimizers = minimizersFiltered;
		//minimizerPos = minimizerPosFiltered;
		//minimizerPos.push_back(readLength);
		//minimizerDirections = minimizerDirectionsFiltered; 
		//minimizerQualities = minimizerQualitiesFiltered; 

		//}
	}

	static MinimizerRead getLowDensityMinimizerRead(const MinimizerRead& highDensityRead, float minimizerDensity){

		//return highDensityRead;
		//cout << highDensityRead._minimizers.size() << endl;
		//return highDensityRead;

		vector<MinimizerType> minimizersFiltered;
		vector<u_int32_t> minimizerPosFiltered; 
		vector<u_int8_t> minimizerDirectionsFiltered; 
		vector<u_int8_t> minimizerQualitiesFiltered; 

		Utils::applyDensityThreshold(minimizerDensity, highDensityRead._minimizers, highDensityRead._minimizersPos, highDensityRead._readMinimizerDirections, highDensityRead._qualities, minimizersFiltered, minimizerPosFiltered, minimizerDirectionsFiltered, minimizerQualitiesFiltered);
		

		MinimizerRead lowDensityRead = {highDensityRead._readIndex, minimizersFiltered, minimizerPosFiltered, minimizerQualitiesFiltered, minimizerDirectionsFiltered, highDensityRead._readLength};
		
		//if(highDensityRead._minimizers.size() != lowDensityRead._minimizers.size()) exit(1);
		//cout << lowDensityRead._minimizers.size() << endl;
		//getchar();

		
		return lowDensityRead;
	}

	static vector<u_int64_t> vec32_to_vec64(const vector<u_int32_t>& vec32) {

		//return vec32;
		
		vector<u_int64_t> vec64(vec32.size(), 0);

		for(size_t i=0; i<vec32.size(); i++){
			vec64[i] = vec32[i];
		}

		return vec64;
		
		
	}

	static string formatTimeToHumanReadable(std::chrono::seconds secs)
	{
		using namespace std;
		using namespace std::chrono;
		//bool neg = secs < 0s;
		//if (neg)
		//	secs = -secs;
		auto h = duration_cast<hours>(secs);
		secs -= h;
		auto m = duration_cast<minutes>(secs);
		secs -= m;
		std::string result;
		//if (neg)
		//	result.push_back('-');
		//if (h < 10h)
		//	result.push_back('0');
		if(h/1h != 0){
			result += to_string(h/1h);
			result += "h ";
		}
		//if (m < 10min)
		//	result.push_back('0');
		if(m/1min != 0){
			result += to_string(m/1min);
			result += "min ";
		}
		//if (secs < 10s)
		//	result.push_back('0');
		result += to_string(secs/1s);
		result += "sec";
		return result;
	}

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

	/*
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
	*/

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


		//const string& failedFilename = outputDir + "/failed.txt";
		static double _maxMemoryUsage = 0;
		static string s = "Maximumresidentsetsize(kbytes):";
		static string stringCpu = "PercentofCPUthisjobgot:";
		static string stringTime = "Elapsed(wallclock)time(h:mm:ssorm:ss):";


		string command2 = "{ \\time -v " + command + "; } 2> " + outputDir + "/time.txt";

		Logger::get().debug() << "";
		Logger::get().debug() << "";
		Logger::get().debug() << command;
		//logFile << endl;
		//logFile << command2 << endl;
		//cout << command2 << endl;
		
		//auto timeStart = high_resolution_clock::now();

		int ret = system(command2.c_str());
		
		//auto timeEnd = high_resolution_clock::now();

		//cout << ret  << " " << (ret != 0) << endl;
		if(ret != 0){
			//logFile.close();

			Logger::get().error() << "";
			Logger::get().error() << "ERROR (see logs: " << Logger::get()._logFilename << ")";
			//cerr << endl;
			//cerr << "ERROR (logs: " << outputDir + "/logs.txt)" << endl;
			//cerr << "Command failed: " << ret << endl;
			//cerr << "Logs: " << outputDir + "/logs.txt" << endl;

			//ofstream outfile(failedFilename);
			//outfile.close();

			exit(ret);
		}

		//if(fs::exists(failedFilename)){
		
		//	logFile << endl;
		//	logFile << command << endl;
		//	logFile.close();
		//	cerr << endl;
		//	cerr << "ERROR (logs: " << outputDir + "/logs.txt)" << endl;
			//cerr << "Logs: " << outputDir + "/logs.txt" << endl;
		//	exit(1);
		//}
		
		ifstream infile(outputDir + "/time.txt");
		string line;
		double maxMem = 0;
		string cpu;
		string time;

		while(std::getline(infile, line)){
			//cout << line << endl;
    		//line.erase(std::remove_if(line.begin(), line.end(), ::istab), line.end());
    		line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
			if(line.empty()) continue;

			size_t pos = line.find(s);
			if (pos != std::string::npos){
				line.erase(pos, s.length());
				maxMem = stod(line);
				maxMem /= 1000;
				maxMem /= 1000;
				_maxMemoryUsage = max(_maxMemoryUsage, maxMem);
			}

			pos = line.find(stringCpu);
			if (pos != std::string::npos){
				line.erase(pos, stringCpu.length());
				cpu = line;
			}

			pos = line.find(stringTime);
			if (pos != std::string::npos){
				line.erase(pos, stringTime.length());
				time = line;
			}
		}
		infile.close();
		//cout << "Max memory: " << _maxMemoryUsage << " GB" << endl;
		//getchar();
		

		//float runtime = duration_cast<seconds>(timeEnd - timeStart).count();


		//cout << outputDir + "/memoryTrack.txt" << endl;
		ofstream outfileMem(outputDir + "/memoryTrack.txt", std::ios_base::app);
		outfileMem << command << endl;
		outfileMem << "\tTime: " << time << endl;
		outfileMem << "\tMemory (GB): " << maxMem << endl;
		outfileMem << "\tCPU: " << cpu << endl;
		outfileMem.close();

		ofstream outfile(outputDir + "/perf.txt");
		outfile << _maxMemoryUsage;
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

	template<typename T>
	static double compute_median(vector<T> scores){
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

	/*
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
	*/

	template<typename T>
	static bool sharedAllElements(const vector<T>& reads1, const vector<T>& reads2){

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
				if(nbShared == reads2.size()) return true;
				return false;
			}

		}

		if(nbShared == reads2.size()) return true;
		return false;
	}

	template<typename T>
	static u_int64_t computeSharedElements(const vector<T>& reads1, const vector<T>& reads2){

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
	/*
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
	*/



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

	/*
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
	*/

	template<typename T>
	static bool shareAny(const vector<T>& utg1, const vector<T>& utg2){

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

/*
struct AlignmentResult2{

	public:

	ReadType _referenceReadIndex;
	ReadType _queryReadIndex;
	//u_int32_t _cigarReferenceStart;
	//u_int32_t _cigarQueryStart;
	//string _cigar;
	//bool _isQueryReversed;
	float _chainingScore;
	int32_t _nbMatches;
	int32_t _nbMissmatches;
	int32_t _nbDeletions;
	int32_t _nbInsertions;
	float _identity;
	u_int32_t _overHangStart;
	u_int32_t _overHangEnd;
	u_int32_t _alignLength;
	vector<pair<int32_t, int32_t>> _alignments;

	AlignmentResult2(){

		_referenceReadIndex = 0;
		_queryReadIndex = 0;
		_chainingScore = 0;
		_nbMatches = 0;
		_nbMissmatches = 0;
		_nbDeletions = 0;
		_nbInsertions = 0;
		_identity = 0;
		_overHangStart = 0;
		_overHangEnd = 0;
		_alignLength = 0;
	}
	
	int64_t getScore() const{
		return _nbMatches - _nbMissmatches - _nbDeletions - _nbInsertions;
	}

	float getSimilarity() const{
		return _nbMatches / (_nbMatches+_nbMissmatches+_nbDeletions+_nbInsertions);// - _nbMissmatches - _nbDeletions - _nbInsertions;
	}

	AlignmentResult2 toSymmetrical(u_int32_t referenceSize, u_int32_t querySize){

		AlignmentResult2 alignment;
		alignment._cigarReferenceStart = 0;
		alignment._cigarQueryStart = 0;
		alignment._cigar = getSymmetricalCigar();

		if(_isQueryReversed){

			if(_cigarQueryStart == 0 && _cigarReferenceStart == 0){
				int64_t r1  = querySize - _cigarQueryStart - getAlignLengthQuery(); //cigarReferenceStart
				int64_t q1 = referenceSize  - getAlignLengthReference(); //alignment._cigarQueryStart
				alignment._cigarReferenceStart = r1;
				alignment._cigarQueryStart = q1;
			}
			else if(_cigarQueryStart > 0){
				int64_t r1  = querySize - _cigarQueryStart - getAlignLengthQuery(); //cigarReferenceStart
				int64_t q1 = referenceSize  - getAlignLengthReference(); //alignment._cigarQueryStart
				alignment._cigarReferenceStart = r1;
				alignment._cigarQueryStart = q1;
			}
			else if(_cigarReferenceStart > 0){
				int64_t r2 = referenceSize - _cigarReferenceStart - getAlignLengthReference(); //alignment._cigarQueryStart
				int64_t q2 = querySize - getAlignLengthQuery(); //alignment._cigarReferenceStart
				alignment._cigarQueryStart = r2;
				alignment._cigarReferenceStart = q2;
			}
			//cout << r1 << " " << q1 << " " << r2 << " " << q2 << endl;
			//int64_t alignLength = getAlignLengthReference();
		}
		else{
			alignment._cigarReferenceStart = _cigarQueryStart;
			alignment._cigarQueryStart = _cigarReferenceStart;
		}


		alignment._referenceReadIndex = _queryReadIndex;
		alignment._queryReadIndex = _referenceReadIndex;
		alignment._isQueryReversed = _isQueryReversed;

		return alignment;
	}

	string getSymmetricalCigar(){
		
		string cigar = "";
		
		vector<Utils::CigarElement> cigarSequence = Utils::extractCigarSequence(_cigar);

		if(_isQueryReversed){
			std::reverse(cigarSequence.begin(), cigarSequence.end());
		}

		for(const Utils::CigarElement& cigarElement : cigarSequence){

			cigar += to_string(cigarElement._nbOccurences);

			if(cigarElement._cigarType == Utils::CigarType::Match){
				cigar += 'M';
			}
			else if(cigarElement._cigarType == Utils::CigarType::Insertion){
				cigar += 'D';
			}
			else if(cigarElement._cigarType == Utils::CigarType::Deletion){
				cigar += 'I';
			}
			
		}

		return cigar;

	}

	int64_t getAlignLengthReference(){

		int64_t alignLength = 0;
		vector<Utils::CigarElement> cigarSequence = Utils::extractCigarSequence(_cigar);

		for(const Utils::CigarElement& cigarElement : cigarSequence){
			if(cigarElement._cigarType == Utils::CigarType::Insertion) continue;
			alignLength += cigarElement._nbOccurences;
		}

		return alignLength;
	}

	int64_t getAlignLengthQuery(){

		int64_t alignLength = 0;
		vector<Utils::CigarElement> cigarSequence = Utils::extractCigarSequence(_cigar);

		for(const Utils::CigarElement& cigarElement : cigarSequence){
			if(cigarElement._cigarType == Utils::CigarType::Deletion) continue;
			alignLength += cigarElement._nbOccurences;
		}

		return alignLength;
	}

	void print(){
		cout << "---" << endl;
		cout << "\tReference index: " << _referenceReadIndex << endl;
		cout << "\tQuery index    : " << _queryReadIndex << endl;
		cout << "\tIs reversed    : " << _isQueryReversed << endl;
		cout << "\tReference start: " << _cigarReferenceStart << endl;
		cout << "\tQuery start    : " << _cigarQueryStart << endl;
		cout << "\tAlign length reference: " << getAlignLengthReference() << endl;
		cout << "\tAlign length query    : " << getAlignLengthQuery() << endl;
		cout << "\tCigar: " << _cigar << endl;

	}

	void getStats(const vector<MinimizerType>& referenceMinimizers, const vector<MinimizerType>& queryMinimizers, int64_t& nbMatches, int64_t& nbMissmatches, int64_t& nbInsertions, int64_t& nbDeletions) const{

		nbMatches = 0;
		nbMissmatches = 0;
		nbInsertions = 0;
		nbDeletions = 0;

		//vector<MinimizerType> queryMinimizers = queryMinimizersOriginal;

		//if(_isQueryReversed){
		//	std::reverse(queryMinimizers.begin(), queryMinimizers.end());
		//}
		
		vector<Utils::CigarElement> cigarSequence = Utils::extractCigarSequence(_cigar);

		size_t referencePosition = _cigarReferenceStart;
		size_t queryPosition = _cigarQueryStart;

		for(const Utils::CigarElement& cigarElement : cigarSequence){

			//cout << endl;
			//cout << referencePosition << " " << referenceMinimizers.size() << endl;
			//cout << queryPosition << " " << queryMinimizers.size() << endl;
			//cout << endl;
			if(cigarElement._cigarType == Utils::CigarType::Match){

				for(size_t i=0; i<cigarElement._nbOccurences; i++){

					if(referenceMinimizers[referencePosition] == queryMinimizers[queryPosition]){ //Match
						nbMatches += 1;
					}
					else{ //missmatch
						nbMissmatches += 1;
					}

					referencePosition += 1;
					if(_isQueryReversed){
						queryPosition -= 1;
					}	
					else{
						queryPosition += 1;
					}					
				}
			}
			else if(cigarElement._cigarType == Utils::CigarType::Insertion){

				for(size_t i=0; i<cigarElement._nbOccurences; i++){

					nbInsertions += 1;

					if(_isQueryReversed){
						queryPosition -= 1;
					}	
					else{
						queryPosition += 1;
					}	
				}
				
			}
			else if(cigarElement._cigarType == Utils::CigarType::Deletion){

				for(size_t i=0; i<cigarElement._nbOccurences; i++){
					nbDeletions += 1;
					referencePosition += 1;
				}


			}

			//getchar();
			//cout << cigarElement._nbOccurences << endl;
			//cout << cigarElement._cigarType << endl;
		}

		

	}
	
	vector<std::pair<int64_t, int64_t>> getSpoaAlignment(const vector<MinimizerType>& referenceMinimizers, const vector<MinimizerType>& queryMinimizersOriginal) const{

		vector<MinimizerType> queryMinimizers = queryMinimizersOriginal;

		if(_isQueryReversed){
			std::reverse(queryMinimizers.begin(), queryMinimizers.end());
		}

		vector<std::pair<int64_t, int64_t>> alignments;

		vector<Utils::CigarElement> cigarSequence = Utils::extractCigarSequence(_cigar);

		size_t referencePosition = _cigarReferenceStart;
		size_t queryPosition = _cigarQueryStart;

		for(const Utils::CigarElement& cigarElement : cigarSequence){

			//cout << endl;
			//cout << referencePosition << " " << referenceMinimizers.size() << endl;
			//cout << queryPosition << " " << queryMinimizers.size() << endl;
			//cout << endl;
			if(cigarElement._cigarType == Utils::CigarType::Match){

				for(size_t i=0; i<cigarElement._nbOccurences; i++){

					if(referenceMinimizers[referencePosition] == queryMinimizers[queryPosition]){ //Match
						alignments.push_back(std::pair(referencePosition, queryPosition));
					}
					else{ //missmatch
						alignments.push_back(std::pair(referencePosition, queryPosition));
					}

					referencePosition += 1;
					//if(_isQueryReversed){
					//	queryPosition -= 1;
					//}	
					//else{
						queryPosition += 1;
					//}					
				}
			}
			else if(cigarElement._cigarType == Utils::CigarType::Insertion){

				for(size_t i=0; i<cigarElement._nbOccurences; i++){

					alignments.push_back(std::pair(-1, queryPosition));

					//if(_isQueryReversed){
					//	queryPosition -= 1;
					//}	
					//else{
						queryPosition += 1;
					//}	
				}
				
			}
			else if(cigarElement._cigarType == Utils::CigarType::Deletion){

				for(size_t i=0; i<cigarElement._nbOccurences; i++){
					alignments.push_back(std::pair(referencePosition, -1));
					referencePosition += 1;
				}


			}

			//getchar();
			//cout << cigarElement._nbOccurences << endl;
			//cout << cigarElement._cigarType << endl;
		}

		return alignments;

	}
	
};
*/

struct KminmerList{
	ReadType _readIndex;
	vector<MinimizerType> _readMinimizers;
	vector<u_int32_t> _minimizerPos;
	vector<u_int8_t> _readMinimizerDirections;
	vector<u_int8_t> _readQualities;
	vector<ReadKminmerComplete> _kminmersInfo;
	Read _read;
	u_int8_t _isCircular;
	//vector<AlignmentResult2> _alignments;
	float _meanReadQuality;
	u_int32_t _readLength;
	//vector<KmerVec> _kminmers;
	//vector<ReadKminmer> _kminmersInfo;
};
/*
struct AlignmentResult3{
	ReadType _queryReadIndex;
	u_int32_t _cigarReferenceStart;
	u_int32_t _cigarQueryStart;
	string _cigar;
	bool _isQueryReversed;
	float _debug_mapScore;
	float _debug_originalOrder;

	void getStats(const vector<MinimizerType>& referenceMinimizers, const vector<MinimizerType>& queryMinimizers, int64_t& nbMatches, int64_t& nbMissmatches, int64_t& nbInsertions, int64_t& nbDeletions) const{

		nbMatches = 0;
		nbMissmatches = 0;
		nbInsertions = 0;
		nbDeletions = 0;

		//vector<MinimizerType> queryMinimizers = queryMinimizersOriginal;

		//if(_isQueryReversed){
		//	std::reverse(queryMinimizers.begin(), queryMinimizers.end());
		//}
		
		vector<Utils::CigarElement> cigarSequence = Utils::extractCigarSequence(_cigar);

		size_t referencePosition = _cigarReferenceStart;
		size_t queryPosition = _cigarQueryStart;

		for(const Utils::CigarElement& cigarElement : cigarSequence){

			//cout << endl;
			//cout << referencePosition << " " << referenceMinimizers.size() << endl;
			//cout << queryPosition << " " << queryMinimizers.size() << endl;
			//cout << endl;
			if(cigarElement._cigarType == Utils::CigarType::Match){

				for(size_t i=0; i<cigarElement._nbOccurences; i++){

					if(referenceMinimizers[referencePosition] == queryMinimizers[queryPosition]){ //Match
						nbMatches += 1;
					}
					else{ //missmatch
						nbMissmatches += 1;
					}

					referencePosition += 1;
					if(_isQueryReversed){
						queryPosition -= 1;
					}	
					else{
						queryPosition += 1;
					}					
				}
			}
			else if(cigarElement._cigarType == Utils::CigarType::Insertion){

				for(size_t i=0; i<cigarElement._nbOccurences; i++){

					nbInsertions += 1;

					if(_isQueryReversed){
						queryPosition -= 1;
					}	
					else{
						queryPosition += 1;
					}	
				}
				
			}
			else if(cigarElement._cigarType == Utils::CigarType::Deletion){

				for(size_t i=0; i<cigarElement._nbOccurences; i++){
					nbDeletions += 1;
					referencePosition += 1;
				}


			}

			//getchar();
			//cout << cigarElement._nbOccurences << endl;
			//cout << cigarElement._cigarType << endl;
		}

		

	}

	void getStatsOld(const vector<MinimizerType>& referenceMinimizers, const vector<MinimizerType>& queryMinimizersOriginal, int64_t& nbMatches, int64_t& nbMissmatches, int64_t& nbInsertions, int64_t& nbDeletions) const{

		nbMatches = 0;
		nbMissmatches = 0;
		nbInsertions = 0;
		nbDeletions = 0;

		vector<MinimizerType> queryMinimizers = queryMinimizersOriginal;

		if(_isQueryReversed){
			std::reverse(queryMinimizers.begin(), queryMinimizers.end());
		}
		
		vector<Utils::CigarElement> cigarSequence = Utils::extractCigarSequence(_cigar);

		
		size_t referencePosition = _cigarReferenceStart;
		size_t queryPosition = _cigarQueryStart;


		for(const Utils::CigarElement& cigarElement : cigarSequence){

			if(cigarElement._cigarType == Utils::CigarType::Match){

				for(size_t i=0; i<cigarElement._nbOccurences; i++){

					if(referenceMinimizers[referencePosition] == queryMinimizers[queryPosition]){ //Match
						nbMatches += 1;
					}
					else{ //missmatch
						nbMissmatches += 1;
					}

					referencePosition += 1;
					queryPosition += 1;
				}
			}
			else if(cigarElement._cigarType == Utils::CigarType::Insertion){

				for(size_t i=0; i<cigarElement._nbOccurences; i++){
					queryPosition += 1;
					nbInsertions += 1;
				}
				
			}
			else if(cigarElement._cigarType == Utils::CigarType::Deletion){

				for(size_t i=0; i<cigarElement._nbOccurences; i++){
					referencePosition += 1;
					nbDeletions += 1;
				}


			}

		}
		


	}
	*/
	/*
	//u_int64_t _readIndex;
	//int64_t _nbMatches;
	//int64_t _nbMissmatches;
	//int64_t _nbInsertions;
	//int64_t _nbDeletions;
	//int64_t _alignLengthBps;

	int64_t score() const{
		int64_t nbErrors = _nbMissmatches + _nbInsertions + _nbDeletions;
		return _nbMatches - nbErrors;
	}

	int64_t nbErrors() const{
		return _nbMissmatches + _nbInsertions + _nbDeletions;
	}

	double divergence() const{

		double nbSeeds = _nbMatches + _nbMissmatches + _nbInsertions;// + _nbDeletions;
		
		if(_nbMatches == nbSeeds) return 0;
		if(_nbMatches == 0) return 1;
		//(95/100)^(1/13)
		//if(_nbMatches == 0) return 1;
		return 1.0 - pow((_nbMatches / nbSeeds), 1.0/13.0);
	}
	*/

//};
/*
//typedef phmap::parallel_flat_hash_map<ReadType, vector<AlignmentResult3>, phmap::priv::hash_default_hash<ReadType>, phmap::priv::hash_default_eq<ReadType>, std::allocator<std::pair<ReadType, vector<AlignmentResult3>>>, 4, std::mutex> AlignmentMap;

class AlignmentFile{

	public:

	//enum Mode{
	//	Read,
	//	Write,
	//};

	string _outputDir;
	//Mode _mode;
	//size_t _nbPartitions;
	ReadType _nbReads;
	ReadType _maxReadPerPartition;
	//int _nbCores;
	//unordered_map<size_t, omp_lock_t> _partitionIndex_to_lock;
	vector<omp_lock_t> _locks;
	vector<gzFile> _files;
	vector<bool> _isFileOpened;
	size_t _nbPartitions;
	//vector<gzFile> _partitionIndex_to_file;
	//unordered_map<size_t, gzFile> _partitionIndex_to_file;
	//size_t _nbOpenedPartitions;
	//gzFile _readFile;
	//vector<omp_lock_t> _locks;
	//vector<gzFile> _files;

	AlignmentFile(const string& outputDir, ReadType nbReads, ReadType maxReadPerPartition, char mode){

		//_readFile = nullptr;
		_outputDir = outputDir;
		_nbReads = nbReads;
		_maxReadPerPartition = maxReadPerPartition;
		//_nbOpenedPartitions = 0;
		_higherReferenceReadIndex = -1;
		_readingCurrentPartitionIndex = 0;
		//_mode = mode;
		//_nbCores = nbCores;
		//_nbPartitions = ceil((long double)totalNbReads / (long double)readPerPartition);
		//_locks.resize(nbCores);

		//for(size_t i=0; i<205; i++){
		//	createPartition(i);
		//}
		_nbPartitions = ceil((long double) nbReads / (long double) maxReadPerPartition);
		cout << "Nb partitions: " << _nbPartitions << endl;

		if(mode == 'w'){ //Write mode


			_locks.resize(_nbPartitions);
			_files.resize(_nbPartitions);
			_isFileOpened.resize(_nbPartitions);

			//cout << _locks.size() << endl;
			for(size_t i=0; i<_locks.size(); i++){
				omp_init_lock(&_locks[i]);
				_isFileOpened[i] = false;
			}
		}
	} 

	~AlignmentFile(){
		for(size_t i=0; i<_locks.size(); i++){
			omp_destroy_lock(&_locks[i]);
		}
	}

	void writeAlignment(const AlignmentResult2& alignment){


		size_t partitionIndex = readIndex_to_partitionIndex(alignment._referenceReadIndex);


		omp_set_lock(&_locks[partitionIndex]);
		
		if(!_isFileOpened[partitionIndex]){
			_isFileOpened[partitionIndex] = true;
			string filename = _outputDir + "/alignmentPart_" + to_string(partitionIndex) + ".gz";
			_files[partitionIndex] = gzopen(filename.c_str(),"wb");
		}
		//if(_partitionIndex_to_file.find(partitionIndex) == _partitionIndex_to_file.end()){
		//	string filename = _outputDir + "/alignmentPart_" + to_string(partitionIndex) + ".gz";
		//	_partitionIndex_to_file[partitionIndex] = gzopen(filename.c_str(),"wb");
		//}
		//cout << "Write a: " << alignment._referenceReadIndex << " " << alignment._queryReadIndex << endl;
		//cout << partitionIndex << endl;
		//gzFile& file = _partitionIndex_to_file[partitionIndex];

		writeAlignment(alignment, _files[partitionIndex]);

		//cout << "Write b: " << alignment._referenceReadIndex << " " << alignment._queryReadIndex << endl;
		omp_unset_lock(&_locks[partitionIndex]);

	}


	void writeAlignment(const AlignmentResult2& alignment, gzFile& file){

		gzwrite(file, (const char*)&alignment._referenceReadIndex, sizeof(alignment._referenceReadIndex));
		gzwrite(file, (const char*)&alignment._queryReadIndex, sizeof(alignment._queryReadIndex));
		gzwrite(file, (const char*)&alignment._cigarReferenceStart, sizeof(alignment._cigarReferenceStart));
		gzwrite(file, (const char*)&alignment._cigarQueryStart, sizeof(alignment._cigarQueryStart));
		u_int32_t cigarSize = alignment._cigar.size();
		gzwrite(file, (const char*)&cigarSize, sizeof(cigarSize));
		gzwrite(file, (const char*)&alignment._cigar[0], cigarSize);
		gzwrite(file, (const char*)&alignment._isQueryReversed, sizeof(alignment._isQueryReversed));
		gzwrite(file, (const char*)&alignment._chainingScore, sizeof(alignment._chainingScore));

	}

	size_t readIndex_to_partitionIndex(const ReadType readIndex){
		return readIndex / _maxReadPerPartition;
	}


	void close(){

		for(size_t i=0; i<_files.size(); i++){
			if(_isFileOpened[i]){
				gzclose(_files[i]);
			}
			_isFileOpened[i] = false;
		}


	}


	size_t _readingCurrentReadIndex;
	size_t _readingCurrentPartitionIndex;
	ReadType _higherReferenceReadIndex;
	vector<AlignmentResult2> _readingPartition;

	void readNextPartition(){

		cout << "Reading partition: " << _readingCurrentPartitionIndex  << " / " << _nbPartitions << endl;

		_readingCurrentReadIndex = 0;
		_readingPartition.clear();

		if(_readingCurrentPartitionIndex >= _nbPartitions) return;
		string filename = _outputDir + "/alignmentPart_" + to_string(_readingCurrentPartitionIndex) + ".gz";
		//if(!fs::exists(filename)) return;

		//cout << "Reading partition: " << _readingCurrentPartitionIndex  << " / " << _nbPartitions << endl;

		gzFile file = gzopen(filename.c_str(),"rb");
		
		while(true){

			ReadType referenceReadIndex;
			ReadType queryReadIndex;
			u_int32_t cigarReferenceStart;
			u_int32_t cigarQueryStart;
			u_int32_t cigarSize;
			string cigar;
			bool isQueryReversed;
			float chainingScore;


			gzread(file, (char*)&referenceReadIndex, sizeof(referenceReadIndex));

			if(gzeof(file)) break;


			gzread(file, (char*)&queryReadIndex, sizeof(queryReadIndex));
			gzread(file, (char*)&cigarReferenceStart, sizeof(cigarReferenceStart));
			gzread(file, (char*)&cigarQueryStart, sizeof(cigarQueryStart));
			//u_int32_t cigarSize = alignment._cigar.size();
			gzread(file, (char*)&cigarSize, sizeof(cigarSize));
			cigar.resize(cigarSize);
			gzread(file, (char*)&cigar[0], cigarSize);
			gzread(file, (char*)&isQueryReversed, sizeof(isQueryReversed));
			gzread(file, (char*)&chainingScore, sizeof(chainingScore));

			//cout << referenceReadIndex << " " << queryReadIndex << " " << cigar << endl;
			AlignmentResult2 alignmentResult = {referenceReadIndex, queryReadIndex, cigarReferenceStart, cigarQueryStart, cigar, isQueryReversed, chainingScore};
			_readingPartition.push_back(alignmentResult);
			
		}

		gzclose(file);

		std::sort(_readingPartition.begin(), _readingPartition.end(), [](const AlignmentResult2 & a, const AlignmentResult2 & b){
			if(a._referenceReadIndex == b._referenceReadIndex){
				return a._queryReadIndex < b._queryReadIndex;
			}
			return a._referenceReadIndex < b._referenceReadIndex;
		});

		if(_readingPartition.size() > 0){
			_higherReferenceReadIndex = _readingPartition[_readingPartition.size()-1]._referenceReadIndex;
		}

		_readingCurrentPartitionIndex += 1;
	}



	//ReadType _readCurrentReferenceIndex;
	//AlignmentResult2 _readAlignmentTmp;
	//vector<AlignmentResult2> _readBuffer;

	void readNext(vector<AlignmentResult2>& alignments, size_t expectedReferenceReadIndex){

		//cout << "Try read: " << expectedReferenceReadIndex << " " << _higherReferenceReadIndex << " " << _higherReferenceReadIndex << endl;
		//cout << _readingCurrentPartitionIndex << " " << _nbOpenedPartitions << endl;
		alignments.clear();

		//if(_readingCurrentPartitionIndex >= _nbOpenedPartitions) return;

		while(_higherReferenceReadIndex == -1 || _higherReferenceReadIndex < expectedReferenceReadIndex){
			readNextPartition();
		}
		
		if(_readingPartition.size() == 0) return; //End of reading

		for(size_t i=_readingCurrentReadIndex; _readingCurrentReadIndex<_readingPartition.size(); _readingCurrentReadIndex++){
			AlignmentResult2& alignmentResult = _readingPartition[_readingCurrentReadIndex];

			if(alignmentResult._referenceReadIndex == expectedReferenceReadIndex){
				alignments.push_back(alignmentResult);
			}
			else{
				break;
			}
		}

		//if(alignments.size() > 0) _readingCurrentReadIndex -= 1;

		//cout << expectedReferenceReadIndex << " " << alignments.size() << "    " << _readingCurrentReadIndex << " " << _readingPartition.size() << endl;
		//getchar();
	

	}


	
};
*/

class EncoderRLE{

public: 

	void execute(const char* sequence, size_t length, string& rleSequence, vector<u_int64_t>& rlePositions, bool useHomopolymerCompression) {
		

		rlePositions.clear();
		rleSequence.clear();

		if(useHomopolymerCompression){

		
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
		else{

			rleSequence = string(sequence);
			for(size_t i=0; i<rleSequence.size(); i++){
				rlePositions.push_back(i);
			}

		}


	}
	/*
	void executeMSR(const char* sequence, size_t length, string& rleSequence, vector<u_int64_t>& rlePositions){


		
		
		rlePositions.clear();
		rleSequence.clear();
		
		char lastChar = '#';
		u_int64_t lastPos = 0;

		for(size_t i=0; i<length; i++){
			//cout << i << " " << length << endl;
			char c = sequence[i];
			//if(c == lastChar) continue;



			if(lastChar != '#'){

				char msr = getMSR(lastChar, c);

				if(msr == ' '){

				}
				else{
					rleSequence += msr;
					rlePositions.push_back(lastPos);
					lastPos = i;
				}
				//rleSequence += lastChar;
				//rlePositions.push_back(lastPos);
				//lastPos = i;
				//cout << lastChar << endl;
			}
			lastChar = c;
		}
		//cout << lastChar << endl;
		//rleSequence += lastChar;
		//rlePositions.push_back(lastPos);
		//rlePositions.push_back(length);

		//cout << rleSequence << endl;

		//getchar();
	}

	char getMSR(char prevChar, char currentChar){

		if(prevChar == 'A' && currentChar == 'A'){
			return 'A';
		}
		else if(prevChar == 'A' && currentChar == 'C'){
			return 'T';
		}
		else if(prevChar == 'A' && currentChar == 'G'){
			return 'T';
		}
		else if(prevChar == 'A' && currentChar == 'T'){
			return ' ';
		}
		else if(prevChar == 'C' && currentChar == 'A'){
			return 'G';
		}
		else if(prevChar == 'C' && currentChar == 'C'){
			return ' ';
		}
		else if(prevChar == 'C' && currentChar == 'G'){
			return ' ';
		}
		else if(prevChar == 'C' && currentChar == 'T'){
			return 'A';
		}
		else if(prevChar == 'G' && currentChar == 'A'){
			return 'A';
		}
		else if(prevChar == 'G' && currentChar == 'C'){
			return ' ';
		}
		else if(prevChar == 'G' && currentChar == 'G'){
			return ' ';
		}
		else if(prevChar == 'G' && currentChar == 'T'){
			return 'A';
		}
		else if(prevChar == 'T' && currentChar == 'A'){
			return ' ';
		}
		else if(prevChar == 'T' && currentChar == 'C'){
			return 'T';
		}
		else if(prevChar == 'T' && currentChar == 'G'){
			return 'C';
		}
		else if(prevChar == 'T' && currentChar == 'T'){
			return 'T';
		}

		cout << "Wrong MSR" << endl;
		return ' ';
	}
	*/
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

			//u_int32_t quality = it.second._quality;

			//if(quality==0){
			//	cout << "omg " << nodeName << endl;
			//	getchar();
			//}
			//bool isReversed;
			//d._kmerVec.normalize(isReversed);

			//vector<u_int64_t> minimizerSeq = d._kmerVec.normalize()._kmers;

			vector<MinimizerType> minimizerSeq = vec._kmers;
			//if(kminmerInfo._isReversed){
			//	std::reverse(minimizerSeq.begin(), minimizerSeq.end());
			//}
			


			u_int16_t size = minimizerSeq.size();
			//_kminmerFile.write((const char*)&size, sizeof(size));
			kminmerFile.write((const char*)&minimizerSeq[0], size*sizeof(MinimizerType));

			kminmerFile.write((const char*)&nodeName, sizeof(nodeName));
			kminmerFile.write((const char*)&abundance, sizeof(abundance));
			//kminmerFile.write((const char*)&quality, sizeof(quality));
		}

		kminmerFile.close();

	}

	void load(const string& filename, bool removeUniqueKminmer){

		ifstream kminmerFile(filename);

		while (true) {

			u_int16_t size = _k;
			//kminmerFile.read((char*)&size, sizeof(size));


			vector<MinimizerType> minimizerSeq;
			minimizerSeq.resize(size);
			kminmerFile.read((char*)&minimizerSeq[0], size*sizeof(MinimizerType));

			if(kminmerFile.eof())break;

			u_int32_t nodeName;
			u_int32_t abundance;
			//u_int32_t quality;
			//bool isReversed = false;

			kminmerFile.read((char*)&nodeName, sizeof(nodeName));
			kminmerFile.read((char*)&abundance, sizeof(abundance));
			//kminmerFile.read((char*)&quality, sizeof(quality));


			//if(removeUniqueKminmer && abundance == 1) continue;
			KmerVec vec;
			vec._kmers = minimizerSeq;

			_dbg_nodes[vec] = {nodeName, abundance};
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



	static void getKminmers(const size_t l, const size_t k, const vector<MinimizerType>& minimizers, const vector<u_int32_t>& minimizersPos, vector<KmerVec>& kminmers, vector<ReadKminmer>& kminmersLength, const vector<u_int64_t>& rlePositions, int readIndex, bool allowPalindrome){

		kminmers.clear();
        kminmersLength.clear();
        bool doesComputeLength = minimizersPos.size() > 0;
		if(minimizers.size() < k) return;




		int i_max = ((int)minimizers.size()) - (int)k + 1;
		for(int i=0; i<i_max; i++){

			KmerVec vec;
			vector<u_int32_t> currentMinimizerIndex;

			int j=i;
			while(true){
				
				if(j >= minimizers.size()) break;

				MinimizerType minimizer = minimizers[j];

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

			
		


	}


	static void getKminmers_complete(const size_t k, const vector<MinimizerType>& minimizers, const vector<u_int32_t>& minimizersPos, vector<ReadKminmerComplete>& kminmers, int readIndex, const vector<u_int8_t>& minimizerQualities){

        kminmers.clear();
		if(minimizers.size() < k) return;

		int i_max = ((int)minimizers.size()) - (int)k + 1;
		for(int i=0; i<i_max; i++){


			KmerVec vec;
			vector<u_int32_t> currentMinimizerIndex;

			int j=i;
			while(true){
				
				if(j >= minimizers.size()) break;

				MinimizerType minimizer = minimizers[j];

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

			

	}

};


/*
class ReadParserParallelPaired{

public:

	string _inputFilename1;
	string _inputFilename2;
	int _nbCores;
	//ofstream& _logFile;
	u_int64_t _maxReads;

	ReadParserParallelPaired(const string& inputFilename1, const string& inputFilename2, int nbCores){

		if(!fs::exists(inputFilename1)){
			Logger::get().error() << "File not found: " << inputFilename1;
			exit(1);
		}

		if(!fs::exists(inputFilename2)){
			Logger::get().error() << "File not found: " << inputFilename2;
			exit(1);
		}

		_inputFilename1 = inputFilename1;
		_inputFilename2 = inputFilename2;
		_nbCores = nbCores;
		_maxReads = 0;
	}

	template<typename Functor>
	void parse(const Functor& functor){

		u_int64_t readIndex = -1;

		Logger::get().debug() << "Parsing file: " << _inputFilename1 << " " << _inputFilename2 << endl;

		gzFile fp1 = gzopen(_inputFilename1.c_str(), "r");
		kseq_t *seq1;
		seq1 = kseq_init(fp1);

		gzFile fp2 = gzopen(_inputFilename2.c_str(), "r");
		kseq_t *seq2;
		seq2 = kseq_init(fp2);

		#pragma omp parallel num_threads(_nbCores)
		{

			bool isEOF = false;
			Functor functorSub(functor);


			Read read1;
			Read read2;

			while(true){


				#pragma omp critical(ReadParserParallelPaired_parse)
				{
					int result1 = kseq_read(seq1);
					int result2 = kseq_read(seq2);
					isEOF = result1 < 0 || result2 < 0;

					if(string(seq1->name.s) != string(seq2->name.s)){
						cerr << "paired read have different headers: " << seq1->name.s << " " << seq2->name.s << endl;
						isEOF = true;
					}

					if(!isEOF){
						readIndex += 1;


						if(seq1->qual.l == 0){
							read1 = {readIndex, string(seq1->name.s), string(seq1->seq.s), "", 0};
							read2 = {readIndex, string(seq2->name.s), string(seq2->seq.s), "", 0};
						}
						else{
							read1 = {readIndex, string(seq1->name.s), string(seq1->seq.s), string(seq1->qual.s), 0};
							read2 = {readIndex, string(seq2->name.s), string(seq2->seq.s), string(seq2->qual.s), 0};
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
				if(_maxReads > 0 && readIndex >= _maxReads) break;

				functorSub(read1, read2);

			}
			

		}
		
		kseq_destroy(seq1);	
		gzclose(fp1);
		kseq_destroy(seq2);	
		gzclose(fp2);
		
		_logFile << "Parsing file done (nb reads: " << (readIndex+1) << ")" << endl;
	}

};
*/


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

	ReadParserParallel(const string& inputFilename, bool isFile, bool isBitset, int nbCores) {

		if(!fs::exists(inputFilename)){
			Logger::get().error() << "File not found: " << inputFilename;
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


			_nbDatasets += 1;
			_filenames.push_back(line);

			/*
    		line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
			if(line.empty()) continue;

			fs::path path(line);

			if(fs::exists (path)){
				
				if(fs::is_directory(path)){
					cerr << line << " is a directory (need fasta/fastq filename)" << line << endl;
					exit(1);
				}

			}
			else{
				cerr << "File not found: " << line << endl;
				exit(1);
			}
			*/
		}


	}



	template<typename Functor>
	void parse(const Functor& functor){

		//#pragma omp for
		//cout << _nbCores << endl;

		u_int64_t datasetIndex = 0;
		u_int64_t readIndex = -1;

		for(const string& filename : _filenames){

			Logger::get().debug() << "Parsing file: " << filename;

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


					#pragma omp critical(ReadParserParallel_parse)
					{
						int result = kseq_read(seq);
						isEOF = result < 0;

						if(!isEOF){
							readIndex += 1;

							string header = "";
							if(seq->comment.l == 0){
								header = string(seq->name.s);
							}
							else{
								header = string(seq->name.s) + " " + string(seq->comment.s);
							}

							if(seq->qual.l == 0){
								read = {readIndex, header, string(seq->seq.s), "", datasetIndex};
							}
							else{
								read = {readIndex, header, string(seq->seq.s), string(seq->qual.s), datasetIndex};
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
			datasetIndex += 1;
		}

		Logger::get().debug() << "Parsing file done (nb reads: " << (readIndex+1) << ")";
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

			//fs::path path(line);

			//if(fs::exists (path)){
				_nbDatasets += 1;
				_filenames.push_back(line);
			//}
			//else{
			//	cerr << "File not found: " << line << endl;
			//	exit(1);
			//}
		}


	}



	void parse(const std::function<void(Read)>& fun){

		u_int64_t readIndex = 0;

		vector<Read> sequenceBuffer;

		for(const string& filename : _filenames){

			Logger::get().debug() << "Parsing file: " << filename;
			//cout << filename << endl;

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
				//int slen = 0, qlen = 0;
				fp = gzopen(filename.c_str(), "r");
				seq = kseq_init(fp);

				while (kseq_read(seq) >= 0){

					string header = "";
					if(seq->comment.l == 0){
						header = string(seq->name.s);
					}
					else{
						header = string(seq->name.s) + " " + string(seq->comment.s);
					}
					
					if(seq->qual.l == 0){
						fun({readIndex, header, string(seq->seq.s)});
					}
					else{
						fun({readIndex, header, string(seq->seq.s), string(seq->qual.s)});
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
	

	/*
	void parseKminmers(const std::function<void(vector<KmerVec>, vector<ReadKminmer>, u_int64_t, u_int64_t, string, string)>& fun){

		MinimizerParser* _minimizerParser = new MinimizerParser(_l, _density);
		EncoderRLE encoderRLE;

		u_int64_t readIndex = 0;
		u_int64_t datasetIndex = 0;

		for(const string& filename : _filenames){
			_logFile << "Parsing file: " << filename << endl;
			//cout << filename << endl;

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

				vector<MinimizerType> minimizers;
				vector<u_int32_t> minimizerPos;
				vector<u_int8_t> minimizers_direction;
				_minimizerParser->parse(rleSequence, minimizers, minimizerPos, minimizers_direction);

				vector<KmerVec> kminmers; 
				vector<ReadKminmer> kminmersInfo;
				MDBG::getKminmers(_l, _k, minimizers, minimizerPos, kminmers, kminmersInfo, rlePositions, readIndex, false);


				fun(kminmers, kminmersInfo, readIndex, datasetIndex, string(read->name.s, strlen(read->name.s)), string(read->seq.s, strlen(read->seq.s)));

				//cout << readIndex << endl;

				readIndex += 1;

				//if(readIndex > 50000) break;
			}
				
			kseq_destroy(read);
			gzclose(fp);

			datasetIndex += 1;

		}

		delete _minimizerParser;
	}
	*/

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

	void parse(const std::function<void(vector<MinimizerType>, vector<KmerVec>, vector<ReadKminmer>, u_int8_t, u_int64_t)>& fun){

		ifstream file_readData(_inputFilename, std::ios::binary);

		u_int64_t readIndex = 0;

		while(true){
			
			
			u_int32_t size;
			vector<MinimizerType> minimizers;
			vector<u_int32_t> minimizerPos; 
			vector<u_int8_t> minimizerQualities;
			
			file_readData.read((char*)&size, sizeof(size));

			if(file_readData.eof())break;

			minimizers.resize(size);
			minimizerPos.resize(size);
			minimizerQualities.resize(size, 0);


			u_int8_t isCircular;
			file_readData.read((char*)&isCircular, sizeof(isCircular));

			//cout << size << " " << isCircular << endl;
			file_readData.read((char*)&minimizers[0], size*sizeof(MinimizerType));
			if(_usePos) file_readData.read((char*)&minimizerPos[0], size*sizeof(u_int32_t));
			//if(_hasQuality) file_readData.read((char*)&minimizerQualities[0], size*sizeof(u_int8_t));

			//if(_isReadProcessed.size() > 0 && _isReadProcessed.find(readIndex) != _isReadProcessed.end()){
			//	readIndex += 1;
			//	continue;
			//}
			/*
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
			*/

			vector<KmerVec> kminmers; 
			vector<ReadKminmer> kminmersInfo;
			vector<u_int64_t> rlePositions;
			MDBG::getKminmers(_l, _k, minimizers, minimizerPos, kminmers, kminmersInfo, rlePositions, 0, false);

			fun(minimizers, kminmers, kminmersInfo, isCircular, readIndex);

			readIndex += 1;

		}

		file_readData.close();

	}

	
	//template<typename SizeType>
	void parseMinspace(const std::function<void(vector<MinimizerType>, vector<ReadKminmerComplete>, u_int8_t, u_int64_t)>& fun){

		ifstream file_readData(_inputFilename, std::ios::binary);

		u_int64_t readIndex = 0;

		while(true){
			
			u_int32_t size;
			vector<MinimizerType> minimizers;
			vector<u_int32_t> minimizerPos; 
			vector<u_int8_t> minimizerQualities;
			
			file_readData.read((char*)&size, sizeof(size));

			if(file_readData.eof())break;

			minimizers.resize(size);
			minimizerPos.resize(size);
			minimizerQualities.resize(size, 0);


			u_int8_t isCircular;
			file_readData.read((char*)&isCircular, sizeof(isCircular));

			//cout << size << " " << isCircular << endl;

			file_readData.read((char*)&minimizers[0], size*sizeof(MinimizerType));
			if(_usePos){
				file_readData.read((char*)&minimizerPos[0], size*sizeof(u_int32_t));
			}
			//if(_hasQuality) file_readData.read((char*)&minimizerQualities[0], size*sizeof(u_int8_t));
			
			/*
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
			*/

			vector<ReadKminmerComplete> kminmersInfo;
			MDBG::getKminmers_complete(_k, minimizers, minimizerPos, kminmersInfo, readIndex, minimizerQualities);

			fun(minimizers, kminmersInfo, isCircular, readIndex);

			readIndex += 1;
		}

		file_readData.close();

	}
	
	void parseSequences(const std::function<void(vector<MinimizerType>, u_int8_t, u_int64_t)>& fun){

		ifstream file_readData(_inputFilename, std::ios::binary);

		u_int64_t readIndex = 0;

		while(true){
			
			u_int32_t size;
			vector<MinimizerType> minimizers;
			vector<u_int16_t> minimizersPosOffsets; 
			vector<u_int8_t> minimizerQualities;
			
			file_readData.read((char*)&size, sizeof(size));

			if(file_readData.eof())break;

			minimizers.resize(size);
			minimizersPosOffsets.resize(size);
			minimizerQualities.resize(size, 0);


			u_int8_t isCircular;
			file_readData.read((char*)&isCircular, sizeof(isCircular));


			file_readData.read((char*)&minimizers[0], size*sizeof(MinimizerType));
			if(_usePos){
				file_readData.read((char*)&minimizersPosOffsets[0], size*sizeof(u_int16_t));
			}
			//if(_hasQuality) file_readData.read((char*)&minimizerQualities[0], size*sizeof(u_int8_t));
			
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


class NodePathParserParallel{

public:

	string _inputFilename;
	int _nbCores;

	NodePathParserParallel(){
	}

	NodePathParserParallel(const string& inputFilename, int nbCores){
		_inputFilename = inputFilename;
		_nbCores = nbCores;
	}

	template<typename Functor>
	void parse(const Functor& functor){

		ifstream file_readData(_inputFilename);

		u_int64_t readIndex = -1;

		#pragma omp parallel num_threads(_nbCores)
		{

			bool isEOF = false;
			Functor functorSub(functor);
			
			u_int32_t size;
			vector<u_int32_t> nodePath;
			u_int8_t isCircular;
			NodePath nodePathObject;

			while(true){
				

				#pragma omp critical(NodePathParserParallel_parse)
				{
					/*
					//vector<u_int64_t> supportingReads;
					u_int64_t size;
					contigFile.read((char*)&size, sizeof(size));
					

					if(contigFile.eof()) break;

					u_int8_t isCircular;
					contigFile.read((char*)&isCircular, sizeof(isCircular));

					nodePath.resize(size);
					//supportingReads.resize(size);
					contigFile.read((char*)&nodePath[0], size * sizeof(u_int32_t));
					*/


					file_readData.read((char*)&size, sizeof(size));

					if(file_readData.eof()) isEOF = true;

					if(!isEOF){

						readIndex += 1;

						nodePathObject = {readIndex};

						nodePath.resize(size);

						
						file_readData.read((char*)&isCircular, sizeof(isCircular));

						file_readData.read((char*)&nodePath[0], size*sizeof(u_int32_t));
						//if(_usePos) file_readData.read((char*)&minimizersPosOffsets[0], size*sizeof(u_int16_t));
						//if(_hasQuality) file_readData.read((char*)&minimizerQualities[0], size*sizeof(u_int8_t));
					}

				}

				

				//if(_isReadProcessed.size() > 0 && _isReadProcessed.find(readIndex) != _isReadProcessed.end()){
				//	readIndex += 1;
				//	continue;
				//}

				//cout << "----" << endl;

				if(isEOF) break;

				nodePathObject._nodePath = nodePath;
				nodePathObject._isCircular = isCircular;

				functorSub(nodePathObject);
			}
		}

		file_readData.close();

	}
};

/*
class KminmerParserParallelCorrection{

public:

	string _inputFilename;
	size_t _l;
	size_t _k;
	bool _usePos;
	int _nbCores;
	bool _hasQuality;

	vector<vector<u_int64_t>>& _mReads;
	MinimizerReadMap& _minimizer_to_readIndex;
	//const MinimizerReadMap& _minimizer_to_readIndex;


	KminmerParserParallelCorrection(const string& inputFilename, size_t l, size_t k, bool usePos,bool hasQuality, int nbCores, vector<vector<u_int64_t>>& mReads, MinimizerReadMap& minimizer_to_readIndex) : _mReads(mReads), _minimizer_to_readIndex(minimizer_to_readIndex){

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
			u_int8_t isCircular;
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

						
						file_readData.read((char*)&isCircular, sizeof(isCircular));

						file_readData.read((char*)&minimizers[0], size*sizeof(u_int64_t));
						//if(_usePos) file_readData.read((char*)&minimizersPosOffsets[0], size*sizeof(u_int16_t));
						//if(_hasQuality) file_readData.read((char*)&minimizerQualities[0], size*sizeof(u_int8_t));
					}

				}

				

				//if(_isReadProcessed.size() > 0 && _isReadProcessed.find(readIndex) != _isReadProcessed.end()){
				//	readIndex += 1;
				//	continue;
				//}

				//cout << "----" << endl;

				if(isEOF) break;

				vector<u_int64_t> minimizersPos; 


				
				unordered_map<u_int32_t, u_int32_t> readIndex_to_matchCount;

				for(u_int64_t minimizer : minimizers){

					if(_minimizer_to_readIndex.find(minimizer) == _minimizer_to_readIndex.end()) continue;

					for(u_int32_t readIndex : _minimizer_to_readIndex[minimizer]){
						readIndex_to_matchCount[readIndex] += 1;
					}
				}

				//cout << "Nb minimizers: " << minimizers.size() << endl;
				//cout << "Total read matches: " << readIndex_to_matchCount.size() << endl;

				vector<u_int32_t> matchingReadIndexes;

				for(const auto& it : readIndex_to_matchCount){

					u_int32_t nbMatches = it.second;
					if(nbMatches < 6) continue;
					//float sim = nbMatches / minimizers.size();
					
					//if(sim < 0.2) continue;

					//cout << "Read index: " << it.first << " " << it.second << endl; 

					matchingReadIndexes.push_back(it.first);
				}

				//cout << "Nb matching reads: " << matchingReadIndexes.size() << endl;


				unordered_map<u_int64_t, u_int32_t> minimizer_to_abundance;

				for(u_int32_t readIndex : matchingReadIndexes){

					const vector<u_int64_t>& readMinimizers = _mReads[readIndex];
					for(u_int64_t minimizer : readMinimizers){
						minimizer_to_abundance[minimizer] += 1;
					}
				}


				vector<float> readAbundances;
				for(u_int64_t minimizer : minimizers){
					readAbundances.push_back(minimizer_to_abundance[minimizer]);
				}

				float readAbundance = Utils::compute_median_float(readAbundances);
				//cout << "Rad abundance: " << Utils::compute_median_float(readAbundances) << endl;

				//getchar();


				unordered_set<u_int64_t> solidMinimizers;

				//cout << "Nb minimizers: " << minimizer_to_abundance.size() << endl;
				float minAbundance =  readAbundance * 0.1;

				for(const auto& it : minimizer_to_abundance){

					u_int64_t minimizer = it.first;
					u_int64_t abundance = it.second;

					if(abundance < minAbundance) continue;

					//cout << minimizer << ": " << abundance << endl;

					solidMinimizers.insert(minimizer);
				}

				vector<u_int64_t> minimizersCorrected;
				for(u_int64_t minimizer : minimizers){
					if(solidMinimizers.find(minimizer) == solidMinimizers.end()) continue;
					minimizersCorrected.push_back(minimizer);
				}


				//for(size_t i=0; i<minimizers.size() && i<minimizersCorrected.size() ; i++){
				//	cout << minimizers[i] << " " << minimizersCorrected[i] << " " << (minimizers[i] == minimizersCorrected[i]) << endl;
				//}

				//getchar();
				//cout << "Nb minimizer (original): " << minimizers.size() << endl;
				//cout << "Nb minimizer (corrected): " << minimizersCorrected.size() << endl;
				//getchar();
				//vector<KmerVec> kminmers; 
				//vector<ReadKminmer> kminmersInfo;
				//vector<u_int64_t> rlePositions;
				vector<ReadKminmerComplete> kminmersInfo;
				//MDBG::getKminmers(_l, _k, minimizers, minimizersPos, kminmers, kminmersInfo, rlePositions, 0, false);
				MDBG::getKminmers_complete(_k, minimizersCorrected, minimizersPos, kminmersInfo, readIndex, minimizerQualities);
				
				//fun(minimizers, kminmers, kminmersInfo, readIndex);
				kminmerList._readMinimizers = minimizersCorrected;
				//kminmerList._kminmers = kminmers;
				kminmerList._kminmersInfo = kminmersInfo;
				kminmerList._isCircular = isCircular;
				functorSub(kminmerList);
			}
		}

		file_readData.close();

	}
};
*/
/*
class MinimizerReadParserParallel{

public:

	string _inputFilename;
	size_t _k;
	bool _usePos;
	int _nbCores;
	bool _hasQuality;
	u_int64_t _chunkSize;
	float _densityThreshold;

	MinimizerReadParserParallel(){
	}

	MinimizerReadParserParallel(const string& inputFilename, size_t k, bool usePos, bool hasQuality, u_int64_t chunkSize, int nbCores){

		if(!fs::exists(inputFilename)){
			cout << "File not found: " << inputFilename << endl;
			exit(1);
		}

		_inputFilename = inputFilename;
		_k = k;
		_usePos = usePos;
		_hasQuality = hasQuality;
		_nbCores = nbCores;
		_chunkSize = chunkSize;
		_densityThreshold = -1;
	}


	template<typename ChuckFunctor>
	void execute(const ChuckFunctor& functor){

		ifstream file_readData(_inputFilename);

		u_int64_t readIndex = 0;


		//bool isEOF = false;
		//Functor functorSub(functor);
		vector<MinimizerType> minimizers;
		vector<u_int32_t> minimizerPos; 
		vector<u_int8_t> minimizerDirections; 
		vector<u_int8_t> minimizerQualities; 
		u_int32_t size;
		//KminmerList kminmerList;
		u_int8_t isCircular;
		//KminmerList kminmer;
		float meanReadQuality;
		u_int32_t readLength;

		vector<MinimizerRead> reads;

		while(true){
			
			file_readData.read((char*)&size, sizeof(size));

			if(file_readData.eof()) break;
			//if(file_readData.eof()) isEOF = true;

			//kminmerList = {readIndex};

			minimizers.resize(size);
			minimizerPos.resize(size);
			minimizerQualities.resize(size);
			minimizerDirections.resize(size);
			
			file_readData.read((char*)&isCircular, sizeof(isCircular));

			file_readData.read((char*)&minimizers[0], size*sizeof(MinimizerType));
			if(_hasQuality) file_readData.read((char*)&minimizerPos[0], size*sizeof(u_int32_t));
			if(_hasQuality) file_readData.read((char*)&minimizerDirections[0], size*sizeof(u_int8_t));
			if(_hasQuality) file_readData.read((char*)&minimizerQualities[0], size*sizeof(u_int8_t));
			if(_hasQuality) file_readData.read((char*)&meanReadQuality, sizeof(meanReadQuality));
			if(_hasQuality) file_readData.read((char*)&readLength, sizeof(readLength));

			if(_densityThreshold != -1){

				vector<MinimizerType> minimizersFiltered;
				vector<u_int32_t> minimizerPosFiltered; 
				vector<u_int8_t> minimizerDirectionsFiltered; 
				vector<u_int8_t> minimizerQualitiesFiltered; 

				Utils::applyDensityThreshold(_densityThreshold, minimizers, minimizerPos, minimizerDirections, minimizerQualities, minimizersFiltered, minimizerPosFiltered, minimizerDirectionsFiltered, minimizerQualitiesFiltered);
				
				minimizers = minimizersFiltered;
				minimizerPos = minimizerPosFiltered;
				minimizerDirections = minimizerDirectionsFiltered; 
				minimizerQualities = minimizerQualitiesFiltered; 

	

			}

			//vector<KmerVec> kminmers; 
			//vector<ReadKminmer> kminmersInfo;
			//vector<u_int64_t> rlePositions;
			//vector<ReadKminmerComplete> kminmersInfo;
			//MDBG::getKminmers(_l, _k, minimizers, minimizersPos, kminmers, kminmersInfo, rlePositions, 0, false);
			//MDBG::getKminmers_complete(_k, minimizers, minimizerPos, kminmersInfo, readIndex, minimizerQualities);
			
			//fun(minimizers, kminmers, kminmersInfo, readIndex);
			//kminmerList._readMinimizers = minimizers;
			//kminmerList._minimizerPos = minimizerPos;
			//kminmerList._readMinimizerDirections = minimizerDirections;
			//kminmerList._readQualities = minimizerQualities;
			//kminmerList._kminmers = kminmers;
			//kminmerList._kminmersInfo = kminmersInfo;
			//kminmerList._isCircular = isCircular;
			//functorSub(kminmerList);

			//if(readIndex > 1000000)
			reads.push_back({readIndex, minimizers, minimizerPos, minimizerQualities, minimizerDirections, meanReadQuality, readLength});

			if(reads.size() >= _chunkSize){
				processChunk(functor, reads);
				reads.clear();
			}

			readIndex += 1;
		}

		file_readData.close();


		if(reads.size() > 0){
			processChunk(functor, reads);
			reads.clear();
		}


	}


	template<typename Functor>
	void processChunk(const Functor& functor, vector<MinimizerRead>& reads){
		functor(reads);
	}



};
*/

class KminmerParserParallel{

public:

	string _inputFilename;
	size_t _l;
	size_t _k;
	bool _usePos;
	int _nbCores;
	bool _hasQuality;
	float _densityThreshold;
	//AlignmentFile* _alignmentFile;

	unordered_set<u_int64_t> _isReadProcessed;

	KminmerParserParallel(){
	}

	KminmerParserParallel(const string& inputFilename, size_t l, size_t k, bool usePos, bool hasQuality, int nbCores){

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
		_densityThreshold = -1;
		//_alignmentFile = nullptr;
	}

	class FunctorChuckDummy
	{
		public:

		void operator () () const {

		}
	};



	template<typename Functor, typename FunctorChunk=FunctorChuckDummy>
	void parse(const Functor& functor, const FunctorChunk& functorChunk=FunctorChuckDummy()){
	//void parse(const std::function<void(vector<u_int64_t>, vector<KmerVec>, vector<ReadKminmer>, u_int64_t)>& fun){

		ifstream file_readData(_inputFilename);

		u_int64_t readIndex = -1;
		u_int64_t nbReadsChunk = 0;

		#pragma omp parallel num_threads(_nbCores)
		{

			bool isEOF = false;
			Functor functorSub(functor);
			vector<MinimizerType> minimizers;
			vector<u_int32_t> minimizerPos; 
			vector<u_int8_t> minimizerDirections; 
			vector<u_int8_t> minimizerQualities; 
			u_int32_t size;
			KminmerList kminmerList;
			u_int8_t isCircular;
			float meanReadQuality;
			u_int32_t readLength;
			//KminmerList kminmer;

			while(true){
				

				#pragma omp critical(KminmerParserParallel_parse)
				{

					
					file_readData.read((char*)&size, sizeof(size));

					if(file_readData.eof()) isEOF = true;

					if(!isEOF){

						if(nbReadsChunk >= 1000){
							functorChunk();
							nbReadsChunk = 0;
						}

						readIndex += 1;
						nbReadsChunk += 1;

						kminmerList = {readIndex};

						minimizers.resize(size);
						minimizerPos.resize(size);
						minimizerQualities.resize(size);
						minimizerDirections.resize(size);
						
						file_readData.read((char*)&isCircular, sizeof(isCircular));

						file_readData.read((char*)&minimizers[0], size*sizeof(MinimizerType));
						if(_hasQuality) file_readData.read((char*)&minimizerPos[0], size*sizeof(u_int32_t));
						if(_hasQuality) file_readData.read((char*)&minimizerDirections[0], size*sizeof(u_int8_t));
						if(_hasQuality) file_readData.read((char*)&minimizerQualities[0], size*sizeof(u_int8_t));
						if(_hasQuality) file_readData.read((char*)&meanReadQuality, sizeof(meanReadQuality));
						if(_hasQuality) file_readData.read((char*)&readLength, sizeof(readLength));
					}

				}

				

				//if(_isReadProcessed.size() > 0 && _isReadProcessed.find(readIndex) != _isReadProcessed.end()){
				//	readIndex += 1;
				//	continue;
				//}

				//cout << "----" << endl;

				if(isEOF) break;

				/*
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
				*/

				if(_densityThreshold != -1){

					vector<MinimizerType> minimizersFiltered;
					vector<u_int32_t> minimizerPosFiltered; 
					vector<u_int8_t> minimizerDirectionsFiltered; 
					vector<u_int8_t> minimizerQualitiesFiltered; 

					Utils::applyDensityThreshold(_densityThreshold, minimizers, minimizerPos, minimizerDirections, minimizerQualities, minimizersFiltered, minimizerPosFiltered, minimizerDirectionsFiltered, minimizerQualitiesFiltered);
					
					minimizers = minimizersFiltered;
					minimizerPos = minimizerPosFiltered;
					minimizerDirections = minimizerDirectionsFiltered; 
					minimizerQualities = minimizerQualitiesFiltered; 

					
					/*
					MinimizerType maxHashValue = -1;
					MinimizerType minimizerBound = _densityThreshold * maxHashValue;


					for(size_t i=0; i<kminmerList._readMinimizers.size(); i++){
						u_int64_t m = kminmerList._readMinimizers[i];
						if(m > minimizerBound) continue;

						minimizersFiltered.push_back(minimizers[i]);
						minimizerPosFiltered.push_back(minimizerPos[i]);
						minimizerDirectionsFiltered.push_back(minimizerDirections[i]);
						minimizerQualitiesFiltered.push_back(minimizerQualities[i]);

					}

					u_int32_t readLength = minimizerPos[minimizerPos.size()-1];
					minimizers = minimizersFiltered;
					minimizerPos = minimizerPosFiltered;
					minimizerPos.push_back(readLength);
					minimizerDirections = minimizerDirectionsFiltered; 
					minimizerQualities = minimizerQualitiesFiltered; 
					*/

				}


				//vector<KmerVec> kminmers; 
				//vector<ReadKminmer> kminmersInfo;
				vector<u_int64_t> rlePositions;
				vector<ReadKminmerComplete> kminmersInfo;
				//MDBG::getKminmers(_l, _k, minimizers, minimizersPos, kminmers, kminmersInfo, rlePositions, 0, false);
				MDBG::getKminmers_complete(_k, minimizers, minimizerPos, kminmersInfo, readIndex, minimizerQualities);
				
				
				//fun(minimizers, kminmers, kminmersInfo, readIndex);
				kminmerList._readMinimizers = minimizers;
				kminmerList._minimizerPos = minimizerPos;
				kminmerList._readMinimizerDirections = minimizerDirections;
				kminmerList._readQualities = minimizerQualities;
				//kminmerList._kminmers = kminmers;
				kminmerList._kminmersInfo = kminmersInfo;
				kminmerList._isCircular = isCircular;
				kminmerList._meanReadQuality = meanReadQuality;
				kminmerList._readLength = readLength;
				functorSub(kminmerList);
			}
		}

		file_readData.close();


		if(nbReadsChunk > 0){
			functorChunk();
			nbReadsChunk = 0;
		}


	}

	template<typename Functor>
	void parseSequences(const Functor& functor){
	//void parse(const std::function<void(vector<u_int64_t>, vector<KmerVec>, vector<ReadKminmer>, u_int64_t)>& fun){

		ifstream file_readData(_inputFilename);

		u_int64_t readIndex = -1;

		#pragma omp parallel num_threads(_nbCores)
		{

			bool isEOF = false;
			Functor functorSub(functor);
			vector<MinimizerType> minimizers;
			vector<u_int8_t> minimizerDirections; 
			vector<u_int32_t> minimizerPos; 
			vector<u_int8_t> minimizerQualities; 
			u_int32_t size;
			KminmerList kminmerList;
			u_int8_t isCircular;
			//vector<AlignmentResult2> alignments;
			float meanReadQuality;
			u_int32_t readLength;
			//KminmerList kminmer;

			while(true){
				

				#pragma omp critical(KminmerParserParallel_parseSequences)
				{

					
					file_readData.read((char*)&size, sizeof(size));

					if(file_readData.eof()) isEOF = true;

					if(!isEOF){

						readIndex += 1;

						kminmerList = {readIndex};

						minimizers.resize(size);
						minimizerPos.resize(size);
						minimizerQualities.resize(size);
						minimizerDirections.resize(size);

						
						file_readData.read((char*)&isCircular, sizeof(isCircular));

						file_readData.read((char*)&minimizers[0], size*sizeof(MinimizerType));
						if(_hasQuality) file_readData.read((char*)&minimizerPos[0], size*sizeof(u_int32_t));
						if(_hasQuality) file_readData.read((char*)&minimizerDirections[0], size*sizeof(u_int8_t));
						if(_hasQuality) file_readData.read((char*)&minimizerQualities[0], size*sizeof(u_int8_t));
						if(_hasQuality) file_readData.read((char*)&meanReadQuality, sizeof(meanReadQuality));
						if(_hasQuality) file_readData.read((char*)&readLength, sizeof(readLength));


						//if(_alignmentFile != nullptr){
						//	_alignmentFile->readNext(alignments, readIndex);
							//cout << "Read al: " << readIndex << " " << alignments.size() << endl;
							//getchar();
						//}

						//cout << readIndex << " " << minimizers.size() << endl;
					}

				}

				

				//if(_isReadProcessed.size() > 0 && _isReadProcessed.find(readIndex) != _isReadProcessed.end()){
				//	readIndex += 1;
				//	continue;
				//}

				//cout << "----" << endl;

				if(isEOF) break;

				/*
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
				kminmerList._isCircular = isCircular;
				*/
				if(_densityThreshold != -1){

					vector<MinimizerType> minimizersFiltered;
					vector<u_int32_t> minimizerPosFiltered; 
					vector<u_int8_t> minimizerDirectionsFiltered; 
					vector<u_int8_t> minimizerQualitiesFiltered; 

					Utils::applyDensityThreshold(_densityThreshold, minimizers, minimizerPos, minimizerDirections, minimizerQualities, minimizersFiltered, minimizerPosFiltered, minimizerDirectionsFiltered, minimizerQualitiesFiltered);
					
					minimizers = minimizersFiltered;
					minimizerPos = minimizerPosFiltered;
					minimizerDirections = minimizerDirectionsFiltered; 
					minimizerQualities = minimizerQualitiesFiltered; 
				}

				kminmerList._readMinimizers = minimizers;
				kminmerList._minimizerPos = minimizerPos;
				kminmerList._readMinimizerDirections = minimizerDirections;
				kminmerList._readQualities = minimizerQualities;
				kminmerList._isCircular = isCircular;
				//kminmerList._alignments = alignments;
				kminmerList._meanReadQuality = meanReadQuality;
				kminmerList._readLength = readLength;

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
/*
class KminmerParserParallelPaired{

public:

	string _inputFilename1;
	string _inputFilename2;
	bool _hasQuality1;
	bool _hasQuality2;
	int _nbCores;


	KminmerParserParallelPaired(){
	}

	KminmerParserParallelPaired(const string& inputFilename1, const string& inputFilename2, bool hasQuality1, bool hasQuality2, int nbCores){

		if(!fs::exists(inputFilename1)){
			cout << "File not found: " << inputFilename1 << endl;
			exit(1);
		}
		
		if(!fs::exists(inputFilename2)){
			cout << "File not found: " << inputFilename2 << endl;
			exit(1);
		}

		_inputFilename1 = inputFilename1;
		_inputFilename2 = inputFilename2;
		_hasQuality1 = hasQuality1;
		_hasQuality2 = hasQuality2;
		_nbCores = nbCores;
	}


	template<typename Functor>
	void parseSequences(const Functor& functor){
		
		ifstream file_readData1(_inputFilename1);
		ifstream file_readData2(_inputFilename2);

		u_int64_t readIndex = -1;

		#pragma omp parallel num_threads(_nbCores)
		{

			bool isEOF = false;
			Functor functorSub(functor);

			u_int32_t size1;
			vector<MinimizerType> minimizers1;
			vector<u_int32_t> minimizerPos1; 
			vector<u_int8_t> minimizerDirections1; 
			vector<u_int8_t> minimizerQualities1; 
			u_int8_t isCircular1;

			u_int32_t size2;
			vector<MinimizerType> minimizers2;
			vector<u_int32_t> minimizerPos2; 
			vector<u_int8_t> minimizerDirections2; 
			vector<u_int8_t> minimizerQualities2; 
			u_int8_t isCircular2;

			KminmerListPaired kminmerList;
			//KminmerList kminmer;

			while(true){
				

				#pragma omp critical(KminmerParserParallelPaired_parseSequences)
				{

					
					file_readData1.read((char*)&size1, sizeof(size1));
					file_readData2.read((char*)&size2, sizeof(size2));

					if(file_readData1.eof()) isEOF = true;
					if(file_readData2.eof()) isEOF = true;

					if(!isEOF){

						readIndex += 1;

						kminmerList = {readIndex};

						minimizers1.resize(size1);
						minimizerPos1.resize(size1+1);
						minimizerQualities1.resize(size1);
						minimizerDirections1.resize(size1);


						minimizers2.resize(size2);
						minimizerPos2.resize(size2+1);
						minimizerQualities2.resize(size2);
						minimizerDirections2.resize(size2);

						//cout << size1 << " " << size2 << endl;
						//getchar();
						file_readData1.read((char*)&isCircular1, sizeof(isCircular1));
						file_readData1.read((char*)&minimizers1[0], size1*sizeof(MinimizerType));
						file_readData1.read((char*)&minimizerPos1[0], (size1+1)*sizeof(u_int32_t));
						if(_hasQuality1) file_readData1.read((char*)&minimizerDirections1[0], size1*sizeof(u_int8_t));
						if(_hasQuality1) file_readData1.read((char*)&minimizerQualities1[0], size1*sizeof(u_int8_t));

						file_readData2.read((char*)&isCircular2, sizeof(isCircular2));
						file_readData2.read((char*)&minimizers2[0], (size2)*sizeof(MinimizerType));
						file_readData2.read((char*)&minimizerPos2[0], (size2+1)*sizeof(u_int32_t));
						if(_hasQuality2) file_readData2.read((char*)&minimizerDirections2[0], size2*sizeof(u_int8_t));
						if(_hasQuality2) file_readData2.read((char*)&minimizerQualities2[0], size2*sizeof(u_int8_t));
					}

				}


				if(isEOF) break;

				kminmerList._readMinimizers1 = minimizers1;
				kminmerList._readMinimizers1pos = minimizerPos1;
				//kminmerList._readMinimizerDirections1 = minimizerDirections1;
				//kminmerList._readQualities1 = minimizerQualities1;
				//kminmerList._isCircular1 = isCircular1;

				kminmerList._readMinimizers2 = minimizers2;
				//kminmerList._readMinimizerDirections2 = minimizerDirections2;
				//kminmerList._readQualities2 = minimizerQualities2;
				//kminmerList._isCircular2 = isCircular2;

				functorSub(kminmerList);
			}
		}

		file_readData1.close();
		file_readData2.close();

	}



};
*/

class Tool{

public:

	//ofstream _logFile;
	//string _tmpDir;
	string _logFilename;

	void run(int argc, char* argv[]){

		try
		{
			parseArgs(argc, argv);
			execute();
			end();
		}
		catch (std::exception& e)
		{
			Logger::get().error() << "ERROR: " << e.what();
			exit(1);
		}

	}

    virtual void parseArgs (int argc, char* argv[]) = 0;
    virtual void execute () = 0;

	void openLogFile(const string& dir){

		string dirNorm = dir;
		while(dirNorm[dirNorm.size()-1] == '/'){
			dirNorm.pop_back();
		}

        fs::path pathNorm = dirNorm;
        std::string parentDir = pathNorm.parent_path().string();
		
		_logFilename = parentDir + "/metaMDBG.log";
		Logger::get().setOutputFile(_logFilename);
	}

	void end(){

		//ofstream file(_tmpDir + "/peakMemory.txt");
		//file << getPeakRSS() << endl;
		//file.close();
		//cout << getPeakRSS()  << " Gb" << endl;
		//getchar();
		//_logFile.close();
	}
	
	

	//void closeLogFile(){
	//	cout << getPeakRSS()  << " Gb" << endl;
	//	getchar();
	//	_logFile.close();
	//}

};




#endif
