

#ifndef MDBG_METAG_COMMONS
#define MDBG_METAG_COMMONS

#include <gatb/gatb_core.hpp>
#include "utils/MurmurHash3.h"

#include <string>
#include <sstream>
#include <unordered_set>
#include <unordered_map>
#include <regex>
#include <algorithm>
#include <libgen.h>
#include <set>


using namespace std;


#define STR_OUTPUT "-o"
#define STR_INPUT "-i"
#define STR_MINIM_SIZE "-l"
#define STR_KMINMER_SIZE "-k"
#define STR_DENSITY "-d"
#define STR_INPUT_DIR "-idir"
#define STR_INPUT_TRUTH "-itruth"
#define STR_HIFIASM_DEBUG "--debug"


//KMER_SPAN(1)
typedef Kmer<>::ModelDirect    ModelDirect;
typedef Kmer<>::ModelCanonical ModelCanonical;
typedef Kmer<>::Type  kmer_type;
typedef Kmer<>::Count kmer_count;
typedef typename Kmer<>::Type  Type;
typedef typename Kmer<>::Count Count;
typedef Kmer<>::ModelMinimizer<ModelCanonical> ModelMinimizer;
//typedef gatb::core::tools::collections::impl OaHash;

typedef u_int32_t ReadIndexType;

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
		return suffix() == prefix().reverse();
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

	static void revcomp(string& sequence){

		std::reverse(sequence.begin(), sequence.end());
		
		for(size_t i=0; i<sequence.size(); i++){
			sequence[i] = basemap[sequence[i]];
		}
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



	static void getKminmers(const size_t l, const size_t k, const vector<u_int64_t>& minimizers, const vector<u_int64_t>& minimizersPos, vector<KmerVec>& kminmers, vector<ReadKminmer>& kminmersLength, const vector<u_int64_t>& rlePositions, int readIndex){


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
							}
							*/

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




								u_int32_t read_pos_start = rlePositions[minimizersPos[indexFirstMinimizer]];
								u_int32_t read_pos_end = rlePositions[minimizersPos[indexLastMinimizer]]; // + (rlePositions[minimizersPos[i+k-1]] - rlePositions[minimizersPos[i+k-1] + l]); //+ l;
								read_pos_end +=  (rlePositions[minimizersPos[indexLastMinimizer] + l] - rlePositions[minimizersPos[indexLastMinimizer]]); //l-1 a check

								//cout << "HI: " << rlePositions[minimizersPos[i+k-1] + l - 1] << endl;
								//cout << "HI: " << rlePositions[minimizersPos[i+k-1] + l] << endl;
								u_int16_t length = read_pos_end - read_pos_start;

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
									position_of_second_minimizer = rlePositions[minimizersPos[indexSecondLastMinimizer]]; //rlePositions[minimizersPos[i+k-1]] - rlePositions[minimizersPos[i+k-2]];
									position_of_second_minimizer_seq = position_of_second_minimizer;
									//position_of_second_minimizer_seq += (rlePositions[minimizersPos[i+k-2] + l - 1] - rlePositions[minimizersPos[i+k-2]]);
									//exit(1);

									u_int16_t pos_last_minimizer = rlePositions[minimizersPos[indexSecondLastMinimizer] + l];
									seq_length_start = read_pos_end - pos_last_minimizer; //rlePositions[minimizersPos[i+k-1]] - rlePositions[minimizersPos[i+k-2]];
								}
								else{

									u_int16_t pos_last_minimizer = rlePositions[minimizersPos[indexSecondMinimizer]];
									seq_length_start = pos_last_minimizer - read_pos_start;

									//cout << "todo" << endl;
									//seq_length_start = 
									//position_of_second_minimizer = rlePositions[minimizersPos[i+1]];// - rlePositions[minimizersPos[i]];
								}

								u_int16_t position_of_second_to_last_minimizer = 0;
								u_int16_t position_of_second_to_last_minimizer_seq = 0;
								if(isReversed){
									
									u_int16_t pos_last_minimizer = rlePositions[minimizersPos[indexSecondMinimizer]];
									seq_length_end = pos_last_minimizer - read_pos_start;

									//position_of_second_to_last_minimizer = rlePositions[minimizersPos[i+1]]; // - rlePositions[minimizersPos[i]];
								}
								else{
									position_of_second_to_last_minimizer = rlePositions[minimizersPos[indexSecondLastMinimizer]]; //rlePositions[minimizersPos[i+k-1]] - rlePositions[minimizersPos[i+k-2]];
									position_of_second_to_last_minimizer_seq = position_of_second_to_last_minimizer;
									//position_of_second_to_last_minimizer_seq += (rlePositions[minimizersPos[i+k-2] + l - 1] - rlePositions[minimizersPos[i+k-2]]);

									u_int16_t pos_last_minimizer = rlePositions[minimizersPos[indexSecondLastMinimizer] + l];
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

};


#endif