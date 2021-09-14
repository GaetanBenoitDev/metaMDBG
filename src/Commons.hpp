

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
	u_int16_t _abundance;
	u_int32_t _length;
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

	void addNode(const KmerVec& vec, u_int32_t length){




		if(_dbg_nodes.find(vec) != _dbg_nodes.end()){
			//if(_dbg_nodes[vec]._abundance > 1000){
			//	cout << _dbg_nodes[vec]._index << " " << _dbg_nodes[vec]._abundance << endl;
			//}
			_dbg_nodes[vec]._abundance += 1;
			return;
		}


		DbgNode node = {_node_id, 1, length};
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

	static void getKminmers(size_t kminmerSize, const vector<u_int64_t>& minimizers, vector<KmerVec>& kminmers){

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
		

		
		if(minimizers.size() < kminmerSize) return;

		//unordered_set<u_int64_t> bannedMinimizers;
		vector<bool> bannedPositions(minimizers.size(), false);

		while(true){

			bool hasPalindrome = false;

			int i_max = ((int)minimizers.size()) - (int)kminmerSize + 1;
			for(int i=0; i<i_max; i++){
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

					if(vec._kmers.size() == kminmerSize){
						if(vec.isPalindrome()){
							//for(size_t p=j-_kminmerSize+1; p<=j-1; p++){
							//	bannedPositions[p] = true;
							//}
							for(size_t m=0; m<kminmerSize-1; m++){
								bannedPositions[currentMinimizerIndex[m]] = true;
								//cout << "Banned: " << currentMinimizerIndex[m] << endl;
							}

							hasPalindrome = true;
							break;
						}
						else{
							kminmers.push_back(vec.normalize());
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


	}

};


#endif