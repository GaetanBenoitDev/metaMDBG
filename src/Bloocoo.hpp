
//Real data
//./bin/Bloocoo  -i ~/workspace/data/overlap_test/input.txt -l 21 -k 3 -d 0.005 -o ~/workspace/run/overlap_test/ -ihifiasm ~/workspace/run/hifiasm_meta/AD_components/big/component_3.fasta -idir ~/workspace/run/overlap_test/

//Simulation
///bin/Bloocoo  -i ~/workspace/data/overlap_test/input.txt -l 16 -k 3 -d 0.005 -o ~/workspace/run/overlap_test_3/ -ihifiasm ~/workspace/data/overlap_test/genome_2371_20x/truth_input.txt -idir ~/workspace/run/overlap_test_3/

#ifndef _BLOOCOO_HPP_
#define _BLOOCOO_HPP_

#include "Commons.hpp"

#include <string>
#include <sstream>
#include <unordered_set>
#include <unordered_map>
#include <regex>
#include <algorithm>
#include <libgen.h>
#include <set>
//#include "graph/Graph.hpp"
#include "graph/GfaParser.hpp"
#include "graph/GraphSimplify.hpp"

//#include <boost/graph/adjacency_list.hpp>
//#include <boost/graph/graph_utility.hpp>
//#include <boost/property_map/property_map.hpp>
//#include <boost/graph/connected_components.hpp>
//#include <boost/graph/dijkstra_shortest_paths.hpp>

// force BOOST ublas optimizations
#define BOOST_UBLAS_INLINE inline
#define BOOST_UBLAS_CHECK_ENABLE 0
#define BOOST_UBLAS_USE_FAST_SAME
#define BOOST_UBLAS_TYPE_CHECK 0
#include <boost/math/distributions.hpp>

#include "utils/MurmurHash3.h"

#define BLOOCOO_VERSION_MAJOR 1
#define BLOOCOO_VERSION_MINOR 0
#define BLOOCOO_VERSION_PATCH 6


//#include <vector>
//#include <bitset>

using namespace std;


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

typedef double Distance;
typedef double Similarity;
typedef boost::math::normal_distribution<Distance> Normal;
typedef boost::math::poisson_distribution<Distance> Poisson;
#define LOG log
#define LOG10 log10
#define SQRT sqrt
#define EXP exp
#define POW pow
#define FABS fabs
static size_t minSamples = 10; //minimum number of sample sizes for considering correlation based recruiting
static int B = 0;
static size_t nABD = 0;
static const size_t nTNF = 136;
static Distance minCV = 1;
static Distance minCVSum = 2;
static double LOG101 = log(101);
static bool sumLowCV = false;

#define STR_MINIM_SIZE "-l"
#define STR_KMINMER_SIZE "-k"
#define STR_OUTPUT "-o"
#define STR_INPUT "-i"
#define STR_DENSITY "-d"
#define STR_INPUT_DIR "-idir"
#define STR_INPUT_EXTRACT_KMINMERS "-ihifiasm"
//typedef boost::adjacency_list<>::vertex_descriptor vertex_t;
//struct VertexData{
//	ReadIndexType _readIndex;
//};

//typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, VertexData> OverlapGraph;

//static boost::numeric::ublas::matrix<float> ABD;
//static boost::numeric::ublas::matrix<float> ABD_VAR;
//static boost::numeric::ublas::matrix<float> TNF;
//typedef gatb::core::tools::collections::impl OaHash;

using namespace std;

/*
std::size_t operator()(std::vector<uint32_t> const& vec) const {

}
*/

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





class MinimizerPairMap{

public:

    MinimizerPairMap(){
        _current_node_id = 0;
    }

    u_int64_t _current_node_id;
    unordered_map<MinimizerPair, u_int32_t> _readIndex_to_id;
    vector<MinimizerPair> _id_to_readIndex;


    //vector<u_int64_t> _unitigs_length;

    void addNode(MinimizerPair& readIndex){
        if (_readIndex_to_id.find(readIndex) != _readIndex_to_id.end()) return;

        //cout << "Add node: " << readIndex << " " << _current_node_id << endl;y
        _readIndex_to_id[readIndex] = _current_node_id;
        _id_to_readIndex.push_back(readIndex);
        _current_node_id += 1;
    }

    u_int64_t pair_to_id(MinimizerPair& readIndex){
        return _readIndex_to_id[readIndex];
    }

    MinimizerPair& id_to_pair(u_int64_t id){
        return _id_to_readIndex[id];
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

};


class CompositionManager{

public:

    ModelCanonical _model;
    ModelCanonical::Iterator _itKmer;
	vector<size_t> _kmer_to_compositionIndex;
	u_int64_t _compositionVectorSize;
    
    CompositionManager(int kmerSize) : _model(kmerSize), _itKmer(_model) {
		
        
        unordered_map<u_int64_t, u_int64_t> setlala;


        int compositionIndex = 0;

        for(int i=0 ;i<pow(4, kmerSize); i++){
            kmer_type kmer;
            //kmer += (u_int64_t) i;
            kmer.setVal((u_int64_t) i);
            kmer_type kmer_min = min(revcomp(kmer, kmerSize), kmer);

            if(setlala.find(kmer_min.getVal()) == setlala.end()){
                setlala[kmer_min.getVal()] = compositionIndex;
                _kmer_to_compositionIndex.push_back(compositionIndex);
                compositionIndex += 1;
            }
            else{
                _kmer_to_compositionIndex.push_back(setlala[kmer_min.getVal()]);
            }

        }


        _compositionVectorSize = setlala.size();
        cout << "Kmer composition size: " << _compositionVectorSize << endl;

    }
	/*
	void readToComposition(Sequence& sequence, vector<float>& composition){
		
		composition.resize(_compositionVectorSize);
		u_int64_t unitigLength = sequence.getDataSize();

		size_t i=0;
		_itKmer.setData (sequence.getData());

		//size_t overlap_size = 0;

		for (_itKmer.first(); !_itKmer.isDone(); _itKmer.next()){
			composition[_kmer_to_compositionIndex[_itKmer->value().getVal()]] += 1;
		}

		Distance rsum = 0;
		for(size_t i = 0; i < composition.size(); ++i) {
			rsum += composition[i] * composition[i];
		}
		rsum = SQRT(rsum);
		for(size_t i = 0; i < composition.size(); ++i) {
			composition[i] /= rsum;
		}


	}
	*/
    
	void readToComposition(Data& sequenceData, u_int32_t sequenceLength, vector<float>& composition){
		
		composition.resize(_compositionVectorSize);
		//u_int64_t unitigLength = sequence.getDataSize();

		size_t i=0;
		_itKmer.setData (sequenceData);

		//size_t overlap_size = 0;

		for (_itKmer.first(); !_itKmer.isDone(); _itKmer.next()){
			composition[_kmer_to_compositionIndex[_itKmer->value().getVal()]] += 1;
		}

		//for(size_t i = 0; i < composition.size(); ++i) {
		//	composition[i] /= sequenceLength;
		//}

		
		Distance rsum = 0;
		for(size_t i = 0; i < composition.size(); ++i) {
			rsum += composition[i] * composition[i];
		}
		rsum = SQRT(rsum);
		for(size_t i = 0; i < composition.size(); ++i) {
			composition[i] /= rsum;
		}
		

	}

};

struct SuccessorData{
	u_int32_t _nodeIndex;
	u_int32_t _abundance;
	u_int32_t _sourceAbundance;
	u_int32_t _prevRank;
	bool _prevRankFinished;
	vector<u_int32_t> _processedNodeIndex;
};






class PathExplorer{

public: 

	//u_int32_t _index;
	//unordered_set<DbgEdge, hash_pair> isEdgeVisited;
	vector<u_int32_t> _prevNodes;
	u_int32_t _source_abundance;
	u_int32_t _source_nodeIndex;
	u_int32_t _start_nodeIndex;
	float _abundanceCutoff_min;
	unordered_set<u_int32_t>& _visitedNodes;

	vector<u_int32_t> _exploredNodes;

	PathExplorer(const vector<u_int32_t>& prevNodes, u_int32_t source_abundance, u_int32_t source_nodeIndex, u_int32_t start_nodeIndex, float abundanceCutoff_min, unordered_set<u_int32_t>& visitedNodes) : _visitedNodes(visitedNodes){
		_prevNodes = prevNodes;
		_source_abundance = source_abundance;
		_source_nodeIndex = source_nodeIndex;
		_start_nodeIndex = start_nodeIndex;
		_abundanceCutoff_min = abundanceCutoff_min;
	}

	u_int32_t getNextNode(u_int32_t current_nodeIndex, GraphSimplify* graph, vector<UnitigData>& _unitigDatas, bool forward, u_int32_t currentDepth){

		if(currentDepth > 5) return -1;

		//u_int64_t iter = 0;
		bool orient_dummy = false;
		vector<u_int32_t> successors;

		//u_int32_t current_nodeIndex = _start_nodeIndex;



		u_int32_t current_nodeName = graph->_graphSuccessors->nodeIndex_to_nodeName(current_nodeIndex, orient_dummy);
		u_int32_t current_abundance = graph->_nodeAbundances[current_nodeName]; //_unitigDatas[current_unitigIndex]._meanAbundance;

		//if(_iter > 10000) return;

		//cout << "----------- " << iter << endl;
		//adjNode* node = graph->_nodes[utg_nodeIndex];
		vector<SuccessorData> data_successors;

		//if(!canExplorePath){
		//	cout << "HAAAAA " << current_nodeName << " " << graph->_graphSuccessors->nodeIndex_to_nodeName(_source_nodeIndex, orient_dummy) << endl;
		//}

		//cout << _prevNodes.size() << endl;
		//if(_prevNodes.size() > 1 && current_nodeIndex == _source_nodeIndex){
		//	cout << "Path complete! " << endl;
		//	//_pathDatas.push_back(pathData);
		//	return -2;
		//}
		//else if(iter > maxIter){
		//	return -1;
		//}


		if(forward){
			graph->getSuccessors(current_nodeIndex, _abundanceCutoff_min, successors);
		}
		else{
			graph->getPredecessors(current_nodeIndex, _abundanceCutoff_min, successors);
		}

		//bool isBranchingNode = successors.size() > 1;
		


		if(currentDepth == 0){
			if(successors.size() > 1){
				for(size_t i=0; i<currentDepth; i++) cout << "  ";
				cout << "----------- " << endl;
				for(u_int32_t utg_n : successors){
				for(size_t i=0; i<currentDepth; i++) cout << "  ";
					cout << graph->_graphSuccessors->nodeToString(current_nodeIndex) << " " << " -> " <<  graph->_graphSuccessors->nodeToString(utg_n) << " " << computeSharedReads(_unitigDatas[current_nodeName], _unitigDatas[graph->_graphSuccessors->nodeIndex_to_nodeName(utg_n, orient_dummy)]) << endl;
				
				}
			}
		}

		for(u_int32_t utg_n : successors){

			//u_int64_t utg_n = node->val;


			u_int32_t successor_nodeName = graph->_graphSuccessors->nodeIndex_to_nodeName(utg_n, orient_dummy);
			u_int32_t successor_abundance = graph->_nodeAbundances[successor_nodeName]; //_unitigDatas[unitigIndex]._meanAbundance;

			if(successors.size() > 1){
				//if(currentDepth == 0){
				if(isPathAlreadyExplored(utg_n, current_nodeIndex, graph, _unitigDatas, forward, currentDepth+1, 20)){
					//cout << "Already explored: " << current_nodeName << " " << successor_nodeName << endl;
					continue;
				}
				//}
			}


			SuccessorData successor = {utg_n, successor_abundance, 0, false};
			successor._sourceAbundance = abs((int)successor._abundance - (int)_source_abundance);
			data_successors.push_back(successor);

		}
		
		if(data_successors.size() == 0){
			if(currentDepth == 0){
				for(size_t i=0; i<currentDepth; i++) cout << "  ";
				cout << "No successors" << endl;
			}
			return -1;
		}
		else if(data_successors.size() == 1){
			current_nodeIndex = data_successors[0]._nodeIndex;
			return current_nodeIndex;
		}
		else{


			/*
			//-------------------------------------------------------------------------------
			//Solve multiple simple cycle
			vector<SuccessorData> successors_nonVisited;
			for(SuccessorData& successor : data_successors){
				if(isPathAlreadyExplored(successor._nodeIndex, current_nodeIndex, graph, _unitigDatas, forward, currentDepth+1, 0)) continue;	
				successors_nonVisited.push_back(successor);
			}

			if(successors_nonVisited.size() == 1){
				for(SuccessorData& successor : successors_nonVisited){
					if(isSmallCycle(successor._nodeIndex, current_nodeIndex, graph, _unitigDatas, forward, currentDepth+1)){
						return successor._nodeIndex;
					}
				}	
			}
	
			//-------------------------------------------------------------------------------
			*/



			u_int32_t currentUnitigIndex = graph->_nodeToUnitig[_prevNodes[_prevNodes.size()-1]];

			u_int32_t prevRank = 0;
			u_int32_t prevRank_unitig = 0;

			while(true){
				

				bool isFinished = true;
				for(size_t i=0; i<data_successors.size(); i++){
					if(!data_successors[i]._prevRankFinished){
						isFinished = false;
					}
				}
				//cout << "        " << isFinished << endl;
				if(isFinished) break;

				int prevIndex = _prevNodes.size() - prevRank - 1;
				if(prevIndex < 0 ) break;

				u_int32_t prev_nodeIndex = _prevNodes[prevIndex];
				u_int32_t prev_nodeName = graph->_graphSuccessors->nodeIndex_to_nodeName(prev_nodeIndex, orient_dummy);

				//cout << prevIndex << " " << prev_nodeIndex << endl;
				//u_int32_t current_unitigIndex = graph->_graphSuccessors->nodeIndex_to_nodeName(current_nodeIndex, orient_dummy);
				if(currentUnitigIndex != graph->_nodeToUnitig[prev_nodeIndex]){
					prevRank_unitig += 1;
					currentUnitigIndex = graph->_nodeToUnitig[prev_nodeIndex];
				}

				//cout << current_nodeName << " " << _node_to_unitig[current_nodeName] << endl;
				string str_debug = "    " + to_string(prevRank_unitig) + ": " +  graph->_graphSuccessors->nodeToString(prev_nodeIndex) + " utg" + to_string(graph->_nodeToUnitig[prev_nodeIndex]);

				for(SuccessorData& successor : data_successors){
					if(successor._prevRankFinished){
						str_debug += "    " + graph->_graphSuccessors->nodeToString(successor._nodeIndex) + " " + "-";
						continue;
					}

					//if(std::find(successor._processedNodeIndex.begin(), successor._processedNodeIndex.end(), prev_nodeIndex) != successor._processedNodeIndex.end()){
					//	successor._prevRankFinished = true;
					//	str_debug += "    " + graph->_graphSuccessors->nodeToString(successor._nodeIndex) + " " + "-";
					//	continue;
					//} 

					u_int32_t successor_nodeName = graph->_graphSuccessors->nodeIndex_to_nodeName(successor._nodeIndex, orient_dummy);
					
					u_int32_t nbSharedReads = computeSharedReads(_unitigDatas[prev_nodeName], _unitigDatas[successor_nodeName]);
					//if(nbSharedReads > _abundanceCutoff_min/2){
					//if(nbSharedReads > successor._abundance/5){
					//if(nbSharedReads > 0){
					if(nbSharedReads > _abundanceCutoff_min/2){
						successor._prevRank = prevRank; //prevRank_unitig; //prevRank; //TODO better comparison in number of nucletoides
						//successor._processedNodeIndex.push_back(prev_nodeIndex);
					}
					else{
						successor._prevRankFinished = true;
					}

					//if(nbSharedReads == 0) continue;
					//if(nbSharedReads > 0 && nbSharedReads >=  data_successors[j]/5){
					//	successor_foundPath[j] = true;
					//	foundPath = true;
					//}

					str_debug += "    " + graph->_graphSuccessors->nodeToString(successor._nodeIndex) + " " + to_string(nbSharedReads);
					//cout << "    " << prevRank << ": " <<  graph->nodeToString(current_nodeIndex) << "     " << graph->nodeToString(successor._nodeIndex)  << ": " << nbSharedReads;
					
					//if(foundPath) break;
				}
				//cout << canExplorePath << endl;
				//cout << current_nodeIndex << " "  <<  _source_nodeIndex << endl;

				if(currentDepth == 0){
					for(size_t i=0; i<currentDepth; i++) cout << "  ";
					cout << str_debug << endl;
				}

				prevRank += 1;

				//cout << currentUnitigIndex << "  " << _node_to_unitig[current_nodeName] << endl;


			}

			std::sort(data_successors.begin(), data_successors.end(), SuccessorComparator_byPrevRank);
			u_int32_t maxPrevRank = data_successors[0]._prevRank;
			vector<SuccessorData> successors_bestPrevRank;
			for(SuccessorData& successor : data_successors){
				if(successor._prevRank > maxPrevRank-5){
					successors_bestPrevRank.push_back(successor);
				}
			}

			if(successors_bestPrevRank.size() == 1){


				//DbgEdge edge = {current_nodeIndex, successors_bestPrevRank[0]._nodeIndex};
				//edge = edge.normalize();
				//isEdgeVisited.insert(edge);


				current_nodeIndex = successors_bestPrevRank[0]._nodeIndex;
				//nodeExplored(current_nodeIndex, graph);

				if(currentDepth == 0){
					for(size_t i=0; i<currentDepth; i++) cout << "  ";
					cout << "Node chosen: " << graph->_graphSuccessors->nodeToString(current_nodeIndex) << " " << graph->_nodeAbundances[graph->_graphSuccessors->nodeIndex_to_nodeName(current_nodeIndex, orient_dummy)]  << endl;
				}

				return current_nodeIndex;
			}
			else{
				



				
				if(currentDepth == 0){

					if(currentDepth == 0){
						for(size_t i=0; i<currentDepth; i++) cout << "  ";
						cout << "Check simple cycle" << endl;
					}

					for(SuccessorData& successor : successors_bestPrevRank){

						if(isSmallCycle(successor._nodeIndex, current_nodeIndex, graph, _unitigDatas, forward, currentDepth+1)){
							return successor._nodeIndex;
						}
						
					}



					for(size_t i=0; i<successors_bestPrevRank.size(); i++){
						u_int32_t to_nodeIndex = successors_bestPrevRank[i]._nodeIndex;
						for(size_t j=0; j<successors_bestPrevRank.size(); j++){
							if(i == j) continue;
							
							u_int32_t from_nodeIndex = successors_bestPrevRank[j]._nodeIndex;

							if(currentDepth == 0){
								for(size_t i=0; i<currentDepth; i++) cout << "  ";
								cout << "Check simple cycle from: " << graph->_graphSuccessors->nodeToString(from_nodeIndex) << "    to:    " << graph->_graphSuccessors->nodeToString(to_nodeIndex) << endl;
							}

							if(isSmallCycle(from_nodeIndex, to_nodeIndex, graph, _unitigDatas, forward, currentDepth+1)){
								//exit(1);
								return from_nodeIndex;
							}

						}	
					}
				}
				

				return -1;
				
			}

		}


	}

	bool isPathAlreadyExplored(u_int32_t current_nodeIndex, u_int32_t source_nodeIndex, GraphSimplify* graph, vector<UnitigData>& _unitigDatas, bool forward, u_int32_t currentDepth, u_int64_t maxIter){
		
		//bool lala = true;

		if(currentDepth == 1){
			for(size_t i=0; i<currentDepth; i++) cout << "  ";
			cout << "Is path explored: " << graph->_graphSuccessors->nodeToString(current_nodeIndex) << " ?    ";
		}

		if(_visitedNodes.find(current_nodeIndex) == _visitedNodes.end()){
			if(currentDepth == 1){
				cout << " No" << endl;
			}
			//lala = false;
			return false;
		}

		if(maxIter == 0){
			if(_visitedNodes.find(current_nodeIndex) == _visitedNodes.end()){
				if(currentDepth == 1){
					cout << " No" << endl;
				}
				return false;
			}
			else{
				if(currentDepth == 1){
					cout << " Yes" << endl;
				}
				return true;
			}
		}

		//u_int32_t source_nodeIndex = current_nodeIndex;
		PathExplorer pathExplorer(_prevNodes, _source_abundance, current_nodeIndex, current_nodeIndex, _abundanceCutoff_min, _visitedNodes);

		//u_int64_t maxIter = 10;


		u_int64_t iter = 0;
		//u_int32_t current_nodeIndex = _start_nodeIndex;
		pathExplorer.nodeExplored(current_nodeIndex, graph);

		//cout << "Start extension: " << graph->_graphSuccessors->nodeToString(current_nodeIndex) << endl;

		while(true){
		
			current_nodeIndex = pathExplorer.getNextNode(current_nodeIndex, graph, _unitigDatas, forward, currentDepth);
			//cout <<  " " << graph->_graphSuccessors->nodeToString(current_nodeIndex);
			if(current_nodeIndex == source_nodeIndex) break;
			if(current_nodeIndex == -1){ //dead end or multiple braching path
				if(currentDepth == 1){
					for(size_t i=0; i<currentDepth; i++) cout << "  ";
					cout << " No" << endl;
				}
				//lala = false;
				return false;
			}

			if(_visitedNodes.find(current_nodeIndex) == _visitedNodes.end()){
				if(currentDepth == 1){
					for(size_t i=0; i<currentDepth; i++) cout << "  ";
					cout << " No" << endl;
				}
				//lala = false;
				return false;
			}

			pathExplorer.nodeExplored(current_nodeIndex, graph);

			if(iter >= maxIter) break;

			iter += 1;
		}
		
		if(currentDepth == 1){
			for(size_t i=0; i<currentDepth; i++) cout << "  ";
			cout << " Yes" << endl;
		}

		//if(!lala) return false;
		return true;
	}

	/*
	bool isSimpleCycle(u_int32_t current_nodeIndex, u_int32_t source_nodeIndex, GraphSimplify* graph, vector<UnitigData>& _unitigDatas, bool forward){
		
		//bool lala = true;

		cout << "\tIs simple cycle: " << graph->_graphSuccessors->nodeToString(current_nodeIndex) << " ?    ";

		PathExplorer pathExplorer(_prevNodes, _source_abundance, current_nodeIndex, current_nodeIndex, _abundanceCutoff_min, _visitedNodes);

		u_int64_t maxIter = 100;
		u_int64_t iter = 0;
		pathExplorer.nodeExplored(current_nodeIndex, graph);

		//cout << "Start extension: " << graph->_graphSuccessors->nodeToString(current_nodeIndex) << endl;

		while(true){
		
			current_nodeIndex = pathExplorer.getNextNode(current_nodeIndex, graph, _unitigDatas, forward, false);
			//cout <<  " " << graph->_graphSuccessors->nodeToString(current_nodeIndex);
			if(current_nodeIndex == source_nodeIndex){
				cout << " Yes" << endl;
				return true;
			}
			if(current_nodeIndex == -1){ //dead end or multiple braching path
				cout << " No" << endl;
				//lala = false;
				return false;
			}



			pathExplorer.nodeExplored(current_nodeIndex, graph);

			if(iter > maxIter) break;

			iter += 1;
		}
		
		cout << " No" << endl;

		//if(!lala) return false;
		return false;
	}
	*/

	bool isSmallCycle(u_int32_t from_nodeIndex, u_int32_t to_nodeIndex, GraphSimplify* graph, vector<UnitigData>& _unitigDatas, bool forward, u_int32_t currentDepth){
		
		unordered_set<u_int32_t> visitedNodes = _visitedNodes;
		//bool lala = true;


		if(currentDepth == 1){
			for(size_t i=0; i<currentDepth; i++) cout << "  ";
			cout << "Is simple cycle: " << graph->_graphSuccessors->nodeToString(from_nodeIndex) << " ?    ";
		}

		PathExplorer pathExplorer(_prevNodes, _source_abundance, from_nodeIndex, from_nodeIndex, _abundanceCutoff_min, visitedNodes);


		u_int64_t maxIter = 100;
		u_int64_t iter = 0;
		pathExplorer.nodeExplored(from_nodeIndex, graph);
		pathExplorer._visitedNodes.insert(from_nodeIndex);

		//cout << "Start extension: " << graph->_graphSuccessors->nodeToString(from_nodeIndex) << endl;

		while(true){
		
			from_nodeIndex = pathExplorer.getNextNode(from_nodeIndex, graph, _unitigDatas, forward, currentDepth);
			//cout << graph->_graphSuccessors->nodeToString(from_nodeIndex) << endl;
			pathExplorer._visitedNodes.insert(from_nodeIndex);

			//cout <<  " " << graph->_graphSuccessors->nodeToString(current_nodeIndex);
			if(from_nodeIndex == to_nodeIndex){
				if(currentDepth == 1){
					for(size_t i=0; i<currentDepth; i++) cout << "  ";
					cout << " Yes" << endl;
				}
				return true;
			}
			if(from_nodeIndex == -1){ //dead end or multiple braching path
				if(currentDepth == 1){
					for(size_t i=0; i<currentDepth; i++) cout << "  ";
					cout << " No" << endl;
				}
				//lala = false;
				return false;
			}



			pathExplorer.nodeExplored(from_nodeIndex, graph);

			if(iter > maxIter) break;

			iter += 1;
		}
		
		if(currentDepth == 1){
			for(size_t i=0; i<currentDepth; i++) cout << "  ";
			cout << " No" << endl;
		}

		//if(!lala) return false;
		return false;
	}

	/*
	bool isPathAlreadyExplored(u_int32_t start_nodeIndex, GraphSimplify* graph, vector<UnitigData>& _unitigDatas, u_int64_t maxIter, bool forward){

		//bool dummy = false;
		//if(graph->_graphSuccessors->nodeIndex_to_nodeName(start_nodeIndex, dummy) == 6307){
		//	cout << "hey" << endl;
		//	exit(1);
		//}

		PathExplorer pathExplorer(_prevNodes, _source_abundance, start_nodeIndex, start_nodeIndex, _abundanceCutoff_min, _visitedNodes);
		pathExplorer.extend(graph, _unitigDatas, maxIter, forward, false);
		
		//cout << "explored" << endl;
		for(u_int32_t nodeIndex : pathExplorer._exploredNodes){
			//cout << "\t\t" << (_visitedNodes.find(nodeIndex) == _visitedNodes.end()) << endl;

			if(_visitedNodes.find(nodeIndex) == _visitedNodes.end()){
				cout << "No" << endl;
				return false;
			}
			//binNode(nodeIndex, pathData.prevNodes, graph, pathData._index);
			//visitedNodes.insert(nodeIndex);
		}

		cout << "Yes" << endl;
		//cout << "Is explored" << endl;
		return true;
	}
	*/

	void nodeExplored(u_int32_t nodeIndex, GraphSimplify* graph){
		
		bool orient_dummy;
		_prevNodes.push_back(nodeIndex);
		_exploredNodes.push_back(nodeIndex);

		u_int32_t nodeName = graph->_graphSuccessors->nodeIndex_to_nodeName(nodeIndex, orient_dummy);


		//_binnedNodes.insert(current_unitigIndex);
		//cout << "Node explored: " << graph->_graphSuccessors->nodeToString(nodeIndex) << " " << graph->_nodeAbundances[nodeName]  << endl;

		//_nbVisitedTimes[current_unitigIndex] += 1;
		//cout << _nbVisitedTimes[current_unitigIndex] << endl;
		//return utg_nodeIndex;
	}

	u_int64_t computeSharedReads(const UnitigData& utg1, const UnitigData& utg2){

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


	bool shareAnyRead(const UnitigData& utg1, const UnitigData& utg2){

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

	static bool SuccessorComparator_byPrevRank(const SuccessorData &a, const SuccessorData &b){
		return a._prevRank > b._prevRank;
	}

};




	/*
	bool extend(GraphSimplify* graph, vector<UnitigData>& _unitigDatas, u_int64_t maxIter, bool forward, bool canExplorePath){

		//if(_start_nodeIndex == 6307){
		//	cout << "HAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" << endl;
		//}
		u_int64_t iter = 0;
		bool orient_dummy = false;
		vector<u_int32_t> successors;

		u_int32_t current_nodeIndex = _start_nodeIndex;

		nodeExplored(current_nodeIndex, graph);
		//unordered_set<DbgEdge, hash_pair>& isEdgeVisited = pathData.isEdgeVisited;
		//vector<u_int32_t>& prevNodes = pathData.prevNodes;
		//u_int32_t source_nodeIndex = pathData.source_nodeIndex_path;
		//u_int32_t source_abundance = pathData.source_abundance;
		//prevNodes.push_back(utg_nodeIndex);

		//binNode(utg_nodeIndex, prevNodes, graph, pathData._index);


		while(true){

			u_int32_t current_nodeName = graph->_graphSuccessors->nodeIndex_to_nodeName(current_nodeIndex, orient_dummy);
			u_int32_t current_abundance = graph->_nodeAbundances[current_nodeName]; //_unitigDatas[current_unitigIndex]._meanAbundance;

			//if(_iter > 10000) return;

			//cout << "----------- " << iter << endl;
			//adjNode* node = graph->_nodes[utg_nodeIndex];
			vector<SuccessorData> data_successors;

			//if(!canExplorePath){
			//	cout << "HAAAAA " << current_nodeName << " " << graph->_graphSuccessors->nodeIndex_to_nodeName(_source_nodeIndex, orient_dummy) << endl;
			//}

			cout << _prevNodes.size() << endl;
			if(_prevNodes.size() > 1 && current_nodeIndex == _source_nodeIndex){
				cout << "Path complete! " << endl;
				//_pathDatas.push_back(pathData);
				return false;
			}
			else if(iter > maxIter){
				return false;
			}


			if(forward){
				graph->getSuccessors(current_nodeIndex, _abundanceCutoff_min, successors);
			}
			else{
				graph->getPredecessors(current_nodeIndex, _abundanceCutoff_min, successors);
			}

			//bool isBranchingNode = successors.size() > 1;
			


			if(successors.size() > 1){
				cout << "----------- " << iter << endl;
				for(u_int32_t utg_n : successors){
					cout << graph->_graphSuccessors->nodeToString(current_nodeIndex) << " " << " -> " <<  graph->_graphSuccessors->nodeToString(utg_n) << " " << computeSharedReads(_unitigDatas[current_nodeName], _unitigDatas[graph->_graphSuccessors->nodeIndex_to_nodeName(utg_n, orient_dummy)]) << endl;
				
				}

			}

			for(u_int32_t utg_n : successors){

				//u_int64_t utg_n = node->val;


				u_int32_t successor_nodeName = graph->_graphSuccessors->nodeIndex_to_nodeName(utg_n, orient_dummy);
				u_int32_t successor_abundance = graph->_nodeAbundances[successor_nodeName]; //_unitigDatas[unitigIndex]._meanAbundance;

				if(successors.size() > 1){
					if(canExplorePath){
						if(isPathAlreadyExplored(utg_n, graph, _unitigDatas, 500, forward)){
							//cout << "Already explored: " << current_nodeName << " " << successor_nodeName << endl;
							continue;
						}
					}
				}

				//graph->_unitigs[graph->_nodeToUnitig[source_nodeIndex]]._abundance

				//u_int32_t nbSharedReads = computeSharedReads(_unitigDatas[current_nodeName], _unitigDatas[successor_nodeName]);
				//cout << current_nodeName << " " << successor_nodeName << " " << nbSharedReads << endl;
				//if(successor_nodeName % 2 == 1){}
				//cout << graph->_graphSuccessors->nodeToString(current_nodeIndex) << " " << current_abundance << " " << " -> " <<  graph->_graphSuccessors->nodeToString(utg_n) << " " << successor_abundance << " " << nbSharedReads << endl;
				//cout << _node_to_unitig[successor_nodeName] << endl;
				//cout << "Neigh: " << node_neighbor << endl;

				
				//if(nbSharedReads == 0) continue;
				
				//if(_node_to_unitig[successor_nodeName] == -1){ //Cleaned
				//if(graph->_isNodeRemoved[utg_n]){ //Cleaned
				//	cout << "\t\tCleaned" << endl;
					//node = node->next;
				//	continue;
				//}

				//if(nbSharedReads <=  current_abundance / (nbNeighbors*4)){
				//	cout << "\t\tAbundance cutoff" << endl;
				//	node = node->next;
				//	continue;
				//}

				//if(nbSharedReads <=  current_abundance / (nbNeighbors*4)){
				//	cout << "\t\tAbundance cutoff" << endl;
				//	node = node->next;
				//	continue;
				//}

				//if(_nbVisitedTimes[current_nodeName] > 100){
				//	cout << "\t\tExiting Infinite cycle" << endl;
				//	continue;
				//}




				SuccessorData successor = {utg_n, successor_abundance, 0, false};
				successor._sourceAbundance = abs((int)successor._abundance - (int)_source_abundance);
				data_successors.push_back(successor);
				//successors.push_back(utg_n);
				//successors_abundance.push_back(successor_abundance);
				//predecessor_rank.push_back(0);
				//successor_finished.push_back(false);

				//node = node->next;
			}
			
			if(data_successors.size() == 0){
				cout << "No successors" << endl;
				return false;
			}
			else if(data_successors.size() == 1){

				//if(isBranchingNode){
				//	DbgEdge edge = {current_nodeIndex, data_successors[0]._nodeIndex};
				//	edge = edge.normalize();
				//	isEdgeVisited.insert(edge);
				//}

				current_nodeIndex = data_successors[0]._nodeIndex;
				nodeExplored(current_nodeIndex, graph);
				cout << graph->_graphSuccessors->nodeToString(current_nodeIndex) << " " << graph->_graphSuccessors->nodeToString(_source_nodeIndex) << endl;
				return true;
			}
			else{


				

				//u_int32_t currentUnitigIndex = _node_to_unitig[graph->_graphSuccessors->nodeIndex_to_nodeName(prevNodes[prevNodes.size()-1], orient_dummy)];
				u_int32_t currentUnitigIndex = graph->_nodeToUnitig[_prevNodes[_prevNodes.size()-1]];

				u_int32_t prevRank = 0;
				u_int32_t prevRank_unitig = 0;

				while(true){
					

					bool isFinished = true;
					for(size_t i=0; i<data_successors.size(); i++){
						if(!data_successors[i]._prevRankFinished){
							isFinished = false;
						}
					}
					//cout << "        " << isFinished << endl;
					if(isFinished) break;

					int prevIndex = _prevNodes.size() - prevRank - 1;
					if(prevIndex < 0 ) break;

					u_int32_t prev_nodeIndex = _prevNodes[prevIndex];
					u_int32_t prev_nodeName = graph->_graphSuccessors->nodeIndex_to_nodeName(prev_nodeIndex, orient_dummy);

					//cout << prevIndex << " " << prev_nodeIndex << endl;
					//u_int32_t current_unitigIndex = graph->_graphSuccessors->nodeIndex_to_nodeName(current_nodeIndex, orient_dummy);
					if(currentUnitigIndex != graph->_nodeToUnitig[prev_nodeIndex]){
						prevRank_unitig += 1;
						currentUnitigIndex = graph->_nodeToUnitig[prev_nodeIndex];
					}

					//cout << current_nodeName << " " << _node_to_unitig[current_nodeName] << endl;
					string str_debug = "    " + to_string(prevRank_unitig) + ": " +  graph->_graphSuccessors->nodeToString(prev_nodeIndex) + " utg" + to_string(graph->_nodeToUnitig[prev_nodeIndex]);

					for(SuccessorData& successor : data_successors){
						if(successor._prevRankFinished){
							str_debug += "    " + graph->_graphSuccessors->nodeToString(successor._nodeIndex) + " " + "-";
							continue;
						}

						//if(std::find(successor._processedNodeIndex.begin(), successor._processedNodeIndex.end(), prev_nodeIndex) != successor._processedNodeIndex.end()){
						//	successor._prevRankFinished = true;
						//	str_debug += "    " + graph->_graphSuccessors->nodeToString(successor._nodeIndex) + " " + "-";
						//	continue;
						//} 

						u_int32_t successor_nodeName = graph->_graphSuccessors->nodeIndex_to_nodeName(successor._nodeIndex, orient_dummy);
						
						u_int32_t nbSharedReads = computeSharedReads(_unitigDatas[prev_nodeName], _unitigDatas[successor_nodeName]);
						//if(nbSharedReads > _abundanceCutoff_min/2){
						//if(nbSharedReads > successor._abundance/5){
						//if(nbSharedReads > 0){
						if(nbSharedReads > 0){
							successor._prevRank = prevRank; //prevRank_unitig; //prevRank; //TODO better comparison in number of nucletoides
							//successor._processedNodeIndex.push_back(prev_nodeIndex);
						}
						else{
							successor._prevRankFinished = true;
						}

						//if(nbSharedReads == 0) continue;
						//if(nbSharedReads > 0 && nbSharedReads >=  data_successors[j]/5){
						//	successor_foundPath[j] = true;
						//	foundPath = true;
						//}

						str_debug += "    " + graph->_graphSuccessors->nodeToString(successor._nodeIndex) + " " + to_string(nbSharedReads);
						//cout << "    " << prevRank << ": " <<  graph->nodeToString(current_nodeIndex) << "     " << graph->nodeToString(successor._nodeIndex)  << ": " << nbSharedReads;
						
						//if(foundPath) break;
					}
					//cout << canExplorePath << endl;
					//cout << current_nodeIndex << " "  <<  _source_nodeIndex << endl;
					cout << str_debug << endl;


					prevRank += 1;

					//cout << currentUnitigIndex << "  " << _node_to_unitig[current_nodeName] << endl;


				}

				std::sort(data_successors.begin(), data_successors.end(), SuccessorComparator_byPrevRank);
				u_int32_t maxPrevRank = data_successors[0]._prevRank;
				vector<SuccessorData> successors_bestPrevRank;
				for(SuccessorData& successor : data_successors){
					if(successor._prevRank > maxPrevRank-5){
						successors_bestPrevRank.push_back(successor);
					}
				}

				if(successors_bestPrevRank.size() == 1){


					//DbgEdge edge = {current_nodeIndex, successors_bestPrevRank[0]._nodeIndex};
					//edge = edge.normalize();
					//isEdgeVisited.insert(edge);


					current_nodeIndex = successors_bestPrevRank[0]._nodeIndex;
					nodeExplored(current_nodeIndex, graph);
					cout << "Node chosen: " << graph->_graphSuccessors->nodeToString(current_nodeIndex) << " " << graph->_nodeAbundances[graph->_graphSuccessors->nodeIndex_to_nodeName(current_nodeIndex, orient_dummy)]  << endl;
					
					return true;
				}
				else{
					
					return false;
					
				}

			}

			iter += 1;

		}

	}
	*/




class Bloocoo : public Tool{
    
public:

	string _inputFilename;
	string _outputDir;
	float _minimizerDensity;
    size_t _minimizerSize;
    size_t _kminmerSize;
    IBank* _inputBank;
	//vector<ReadData> _readData;
	CompositionManager* _compositionManager;
	AdjGraph* _overlapGraph;
	string _inputDir;
	string _input_extractKminmers;
	
	string _filename_readMinimizers;
	string _filename_filteredMinimizers;
	string _filename_hifiasmGroundtruth;
	string _filename_readCompositions;

	vector<u_int32_t> _evaluation_readToDataset;
	unordered_map<KmerVec, u_int16_t> _evaluation_hifiasmGroundTruth;
	unordered_map<KmerVec, u_int32_t> _evaluation_hifiasmGroundTruth_position;
	MinimizerPairMap* _minimizerPairMap;

    Bloocoo ();
    void execute ();
	void createSimilarityGraph(GraphInfo* graphInfo);
	void createGroundTruth();
	void execute_binning();
	void execute_binning_cleanGraph();
	void extract_kminmers();

	static bool Ralalalala(pair<u_int64_t, u_int64_t> a, pair<u_int64_t,u_int64_t> b) {
		return a.second > b.second;
	}



    inline float euclidianDistance(const vector<float>& v1, const vector<float>& v2){
        float sum = 0;
        for(int i=0; i<v1.size(); i++) {
            sum += pow(v2[i] - v1[i], 2);
        }
        return sqrt(sum);
    }

	float getCompositionDistance(ReadData& ni, ReadData& nj){
		
		return euclidianDistance(ni._composition, nj._composition);
        /*
		if(_unitigCompositions.find(ni) == _unitigCompositions.end()) return 1;
		if(_unitigCompositions.find(nj) == _unitigCompositions.end()) return 1;


		float distance = -1;
		u_int64_t edge_1 = (u_int64_t) ni << 32 | nj;
		
		if (_cache_linkDistances.find(edge_1) != _cache_linkDistances.end()){
			distance = _cache_linkDistances[edge_1];
		}
		else{
			u_int64_t edge_2 = (u_int64_t) nj << 32 | ni;
			distance = euclidianDistance(_unitigCompositions[ni], _unitigCompositions[nj]);
			_cache_linkDistances[edge_1] = distance;
			_cache_linkDistances[edge_2] = distance;
		}

		return distance;
        */
	} 

    Distance computeDistanceTNF(ReadData& ni, ReadData& nj) {

		/*
		Distance d = 0;

		for (size_t i = 0; i < nTNF; ++i) {
			d += (TNF(r1,i) - TNF(r2,i)) * (TNF(r1,i) - TNF(r2,i)); //euclidean distance
		}

		d = SQRT(d);
		*/

		Distance d = getCompositionDistance(ni, nj);
		//cout << d << endl;
		Distance b,c; //parameters

		u_int32_t ctg1 = std::min(ni._length, (u_int32_t)500000);
		u_int32_t ctg2 = std::min(nj._length, (u_int32_t)500000);
		//cout << ctg1 << " " << ctg2 << endl;
		Distance lw11 = LOG10(std::min(ctg1, ctg2));
		Distance lw21 = LOG10(std::max(ctg1, ctg2));
		Distance lw12 = lw11 * lw11;
		Distance lw13 = lw12 * lw11;
		Distance lw14 = lw13 * lw11;
		Distance lw15 = lw14 * lw11;
		Distance lw16 = lw15 * lw11;
		Distance lw17 = lw16 * lw11;
		Distance lw22 = lw21 * lw21;
		Distance lw23 = lw22 * lw21;
		Distance lw24 = lw23 * lw21;
		Distance lw25 = lw24 * lw21;
		Distance lw26 = lw25 * lw21;

		Distance prob;

		b = 46349.1624324381 + -76092.3748553155*lw11 + -639.918334183*lw21 + 53873.3933743949*lw12 + -156.6547554844*lw22 + -21263.6010657275*lw13 + 64.7719132839*lw23 +
				5003.2646455284*lw14 + -8.5014386744*lw24 + -700.5825500292*lw15 + 0.3968284526*lw25 + 54.037542743*lw16 + -1.7713972342*lw17 + 474.0850141891*lw11*lw21 +
				-23.966597785*lw12*lw22 + 0.7800219061*lw13*lw23 + -0.0138723693*lw14*lw24 + 0.0001027543*lw15*lw25;
		c = -443565.465710869 + 718862.10804858*lw11 + 5114.1630934534*lw21 + -501588.206183097*lw12 + 784.4442123743*lw22 + 194712.394138513*lw13 + -377.9645994741*lw23 +
				-45088.7863182741*lw14 + 50.5960513287*lw24 + 6220.3310639927*lw15 + -2.3670776453*lw25 + -473.269785487*lw16 + 15.3213264134*lw17 + -3282.8510348085*lw11*lw21 +
				164.0438603974*lw12*lw22 + -5.2778800755*lw13*lw23 + 0.0929379305*lw14*lw24 + -0.0006826817*lw15*lw25;

		//logistic model
		prob = 1.0 / ( 1 + EXP(-(b + c * d)) );

		if(prob >= .1) { //second logistic model
			b = 6770.9351457442 + -5933.7589419767*lw11 + -2976.2879986855*lw21 + 3279.7524685865*lw12 + 1602.7544794819*lw22 + -967.2906583423*lw13 + -462.0149190219*lw23 +
					159.8317289682*lw14 + 74.4884405822*lw24 + -14.0267151808*lw15 + -6.3644917671*lw25 + 0.5108811613*lw16 + 0.2252455343*lw26 + 0.965040193*lw12*lw22 +
					-0.0546309127*lw13*lw23 + 0.0012917084*lw14*lw24 + -1.14383e-05*lw15*lw25;
			c = 39406.5712626297 + -77863.1741143294*lw11 + 9586.8761567725*lw21 + 55360.1701572325*lw12 + -5825.2491611377*lw22 + -21887.8400068324*lw13 + 1751.6803621934*lw23 +
					5158.3764225203*lw14 + -290.1765894829*lw24 + -724.0348081819*lw15 + 25.364646181*lw25 + 56.0522105105*lw16 + -0.9172073892*lw26 + -1.8470088417*lw17 +
					449.4660736502*lw11*lw21 + -24.4141920625*lw12*lw22 + 0.8465834103*lw13*lw23 + -0.0158943762*lw14*lw24 + 0.0001235384*lw15*lw25;
			prob = 1.0 / ( 1 + EXP(-(b + c * d)) );
			prob = prob < .1 ? .1 : prob;
		}

		//cout << prob << " " << d << " " << b << " " << c << endl;
		return prob;
	}


	bool checkRemoveNode(u_int32_t n, AdjGraph* graph, MinimizerPairMap* minimizerMap, unordered_map<MinimizerPair, u_int16_t>& minimizerCount){
		
		u_int16_t abundance_n = minimizerCount[minimizerMap->id_to_pair(n)];
		if(abundance_n <= 1) return true;

		float meanAbundance = 0;
		u_int16_t nbNeighbors = 0;

		adjNode* node = graph->_nodes[n];
        while (node != nullptr) {

			u_int32_t nn = node->val;

			u_int16_t abundance_nn = minimizerCount[minimizerMap->id_to_pair(nn)];

			meanAbundance += abundance_nn;
			nbNeighbors += 1;
			//if(abundance_nn > maxAbundance) maxAbundance = abundance_nn;
            node = node->next;
        }

		if(nbNeighbors == 0){
			if(abundance_n <= 1){
				return true;
			}
			return false;
		}

		meanAbundance /= nbNeighbors;

		if(abundance_n <= 2 && meanAbundance > 10) return true;
		//if(abundance_n < (meanAbundance-15)) return true;
		//int diff = abs(abundance_n - meanAbundance);
		//cout << abundance_n << " " << meanAbundance << " " << diff << endl;
		//if(diff > 10) return true;
		return false;
	}

	void print_stats(vector<float>& elems){

		//cout << "allo " << elems.size() << endl;
		float mean = 0;
		
		for(float e : elems){
			mean += e;
		}
		mean /= elems.size();

		float var = 0;
		for(float e : elems){
			var += pow((e - mean), 2);
		}
		var /= (elems.size() - 1);
		float sd = SQRT(var);
		
		cout << "Mean: " << mean << endl;
		cout << "Sd: " << sd << endl;
		cout << "Var: " << var << endl;
	}

	static bool OverlapComparator(const Overlap &a, const Overlap &b)
	{
		return a._nbMinimizers > b._nbMinimizers;
	}

	u_int64_t computeSharedReads(const UnitigData& utg1, const UnitigData& utg2){

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


	bool shareAnyRead(const UnitigData& utg1, const UnitigData& utg2){

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

	static bool UnitigEdgeScoreComparator(const UnitigEdgeScore& a, const UnitigEdgeScore& b)
	{
		return a._score > b._score;
	}



	void getKminmers_filterRepeatedEdge(const vector<u_int64_t>& minimizers, const unordered_set<u_int64_t>& filteredMinimizers, vector<KmerVec>& kminmers, unordered_map<KmerVec, u_int16_t>& kminmers_counts){

		
		
		if(minimizers.size() < _kminmerSize) return;

		


		
		kminmers.clear();

		//unordered_set<u_int64_t> bannedMinimizers;
		vector<bool> bannedPositions(minimizers.size(), false);

		while(true){



			//cout << "------------------------" << endl;
			/*
			getKminmers(minimizers, filteredMinimizers, kminmers);
			int lala = 0;
			for(KmerVec& vec : kminmers){
				cout << lala << " " << vec.toString() << endl;
				u_int16_t kminmer_abundance = kminmers_counts[vec.normalize()];
				u_int16_t prefix_abundance = kminmersPreSuf_counts[vec.normalize().prefix()];
				u_int16_t suffix_abundance = kminmersPreSuf_counts[vec.normalize().suffix()];
				cout << kminmer_abundance << " " << prefix_abundance << " " << suffix_abundance << endl;
				lala += 1;
			}*/


			bool hasPalindrome = false;

			int i_max = ((int)minimizers.size()) - (int)_kminmerSize + 1;
			for(int i=0; i<i_max; i++){
				if(bannedPositions[i]) continue;

				KmerVec vec;

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

					//cout << vec.toString() << endl;
					if(vec._kmers.size() == _kminmerSize){
						if(vec.isPalindrome()){
							for(size_t p=j-_kminmerSize+1; p<=j; p++){
								bannedPositions[p] = true;
							}
							hasPalindrome = true;
							break;
						}
						else{

							/*
							bool canCheckRepeat = true;
							u_int16_t kminmer_abundance = 0;
							u_int16_t prefix_abundance = 0;
							u_int16_t suffix_abundance = 0;
							if(kminmers_counts.find(vec.normalize()) != kminmers_counts.end()){
								kminmer_abundance = kminmers_counts[vec.normalize()];
							}
							else{
								//cout << "cant 1" << endl;
								canCheckRepeat = false;
							}
							if(kminmersPreSuf_counts.find(vec.normalize().prefix()) != kminmersPreSuf_counts.end()){
								prefix_abundance = kminmersPreSuf_counts[vec.normalize().prefix()];
							}
							else{
								//cout << "cant 2" << endl;
								canCheckRepeat = false;
							}
							if(kminmersPreSuf_counts.find(vec.normalize().suffix()) != kminmersPreSuf_counts.end()){
								suffix_abundance = kminmersPreSuf_counts[vec.normalize().suffix()];
							}
							else{
								//cout << "cant 3" << endl;
								canCheckRepeat = false;
							}

							//cout << "test" << endl;
							//cout << kminmer_abundance << " " << prefix_abundance << " " << suffix_abundance << endl;
							//cout << (prefix_abundance > kminmer_abundance*10) << " " << (suffix_abundance > kminmer_abundance*10) << endl;
							//u_int16_t kminmer_abundance = kminmers_counts[vec.normalize()];
							//u_int16_t prefix_abundance = kminmersPreSuf_counts[vec.normalize().prefix().normalize()];
							//u_int16_t suffix_abundance = kminmersPreSuf_counts[vec.normalize().suffix().normalize()];

							if(canCheckRepeat){
								if(prefix_abundance > kminmer_abundance*20){
									
									if(vec.h() < vec.reverse().h()){
										//cout << "banned 1: " << j << endl;
										bannedPositions[j] = true;
										bannedPositions[j-1] = true;
									}
									else{
										//cout << "banned 2: " << (j-_kminmerSize+1) << endl;
										bannedPositions[j] = true;
										bannedPositions[j-1] = true;
									}

									//for(size_t p=j-_kminmerSize+1; p<=j; p++){
									//	bannedPositions[p] = true;
									//	cout << "banned: " << p << endl;
									//}
									hasPalindrome = true;
									break;

								} 
								else if(suffix_abundance > kminmer_abundance*20){

									if(vec.h() < vec.reverse().h()){
										//cout << "banned 3: " << j << endl;
										bannedPositions[j] = true;
										bannedPositions[j-1] = true;
									}
									else{
										//cout << "banned 4: " << (j-_kminmerSize+1) << endl;
										bannedPositions[j] = true;
										bannedPositions[j-1] = true;
									}

									//for(size_t p=j-_kminmerSize+1; p<=j; p++){
									//	bannedPositions[p] = true;
									//	cout << "banned: " << p << endl;
									//}
									hasPalindrome = true;
									break;
								}
							}
							*/
							if(kminmers_counts.find(vec.normalize()) != kminmers_counts.end()){
								u_int16_t kminmer_abundance = kminmers_counts[vec.normalize()];
								if(kminmer_abundance > 100){
									cout << kminmer_abundance << endl;
									for(size_t p=j-_kminmerSize+1; p<=j; p++){
										bannedPositions[p] = true;
									}
									hasPalindrome = true;
									break;
								}
							}



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

	void getKminmers(const vector<u_int64_t>& minimizers, const unordered_set<u_int64_t>& filteredMinimizers, vector<KmerVec>& kminmers){

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
		

		
		if(minimizers.size() < _kminmerSize) return;

		//unordered_set<u_int64_t> bannedMinimizers;
		vector<bool> bannedPositions(minimizers.size(), false);

		while(true){

			bool hasPalindrome = false;

			int i_max = ((int)minimizers.size()) - (int)_kminmerSize + 1;
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

					if(vec._kmers.size() == _kminmerSize){
						if(vec.isPalindrome()){
							//for(size_t p=j-_kminmerSize+1; p<=j-1; p++){
							//	bannedPositions[p] = true;
							//}
							for(size_t m=0; m<_kminmerSize-1; m++){
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

	/*
	void findpaths(vector<vector<int> >&g, int src, int dst, int v){
		// create a queue which stores
		// the paths
		queue<vector<int> > q;
	
		// path vector to store the current path
		vector<int> path;
		path.push_back(src);
		q.push(path);
		while (!q.empty()) {
			path = q.front();
			q.pop();
			int last = path[path.size() - 1];
	
			// if last vertex is the desired destination
			// then print the path
			if (last == dst)
				printpath(path);       
	
			// traverse to all the nodes connected to
			// current vertex and push new path to queue
			for (int i = 0; i < g[last].size(); i++) {
				if (isNotVisited(g[last][i], path)) {
					vector<int> newpath(path);
					newpath.push_back(g[last][i]);
					q.push(newpath);
				}
			}
		}
	}*/

	/*
	vector<bool> isVisited;
    vector<int> distance;
    vector<int> prev;

    int shortest_path(const OverlapGraph& graph, u_int64_t src, u_int64_t dest, vector<u_int64_t>& path){

		path.clear();
        for(size_t i=0; i<boost::num_vertices(graph); i++){
            isVisited[i] = false;
            distance[i] = 0;
            prev[i] = -1;
        }
        // Keep track of visited nodes

        // Initialize initial distances as 0 for all nodes
		queue <int> queue1;

        distance[src] = 0;
        queue1.push(src);
        isVisited[src] = true;
        bool found = false;

        while (!queue1.empty() && !found){

            u_int64_t node_current = queue1.front();
            //cout << id_to_name(node_current) << endl;

            queue1.pop();

			auto neighbours = boost::adjacent_vertices(node_current, graph);
			for (auto nn : make_iterator_range(neighbours)){


                //u_int64_t node_neighbor = node->val;
                //cout << id_to_name(node_neighbor) << endl;

                if (isVisited[nn]){
                    continue;
                }

                distance[nn] = distance[node_current] + 1;
                queue1.push(nn);
                isVisited[nn] = true;
                prev[nn] = node_current;

                if(nn == dest){
                    found = true;
                    break;
                }

            }



        }

        if(found){

            path.clear();
            u_int64_t n = dest;
            while(n != src){
                path.push_back(n);
                //cout << id_to_name(n) << endl;
                n = prev[n];

            }
            path.push_back(src);

            return distance[dest];
        }

        return -1;
    }

    void collectNeighbors(const OverlapGraph& graph, u_int64_t src, u_int64_t maxDistance, vector<u_int64_t>& neighbors){

		neighbors.clear();

        for(size_t i=0; i<boost::num_vertices(graph); i++){
            isVisited[i] = false;
            distance[i] = -1;
            //prev[i] = -1;
        }
        // Keep track of visited nodes

        // Initialize initial distances as 0 for all nodes
		queue <int> queue1;

        distance[src] = 0;
        queue1.push(src);
        isVisited[src] = true;

		int iter = 0;
        while (!queue1.empty()){

            u_int64_t node_current = queue1.front();
            //cout << id_to_name(node_current) << endl;

            queue1.pop();

			auto neighbours = boost::adjacent_vertices(node_current, graph);
			for (auto nn : make_iterator_range(neighbours)){

				iter += 1;

                //u_int64_t node_neighbor = node->val;
                //cout << id_to_name(node_neighbor) << endl;

                if (isVisited[nn]){
                    continue;
                }

                distance[nn] = distance[node_current] + 1;

				//cout << distance[nn] << " " << maxDistance << endl;

                isVisited[nn] = true;
                //prev[nn] = node_current;

                neighbors.push_back(nn);

				//cout << neighbors.size() << endl;
				//cout << queue1.size() << endl;

				
				if(distance[nn] >= maxDistance) continue;
                queue1.push(nn);


            }



        }

		//cout << "iter: " <<  iter << << endl;

    }*/




	/*
	unordered_set<u_int32_t> isVisited;
	//vector<bool> isVisited;
	vector<u_int16_t> distance;

	u_int32_t shortest_path(u_int32_t src, u_int32_t dest, u_int16_t maxDistance){

        //for(size_t i=0; i<_dbg_nodes.size(); i++){
        //    isVisited[i] = false;
        //    distance[i] = 0;
        //}
        // Keep track of visited nodes

        // Initialize initial distances as 0 for all nodes

        queue <int> queue1;
        distance[src] = 0;
        queue1.push(src);
        isVisited[src] = true;
        bool found = false;

        while (!queue1.empty() && !found){

            u_int64_t node_current = queue1.front();
            //cout << id_to_name(node_current) << endl;

            queue1.pop();

            adjNode* node = _nodes[node_current];
            while (node != nullptr) {

                u_int64_t node_neighbor = node->val;
                //cout << id_to_name(node_neighbor) << endl;

                if (isVisited[node_neighbor]){
                    node = node->next;
                    continue;
                }

                distance[node_neighbor] = distance[node_current] + 1;
                queue1.push(node_neighbor);
                isVisited[node_neighbor] = true;

                if(node_neighbor == dest){
                    found = true;
                    break;
                }

                node = node->next;
            }



        }

        return -1;
    }
	*/

	double compute_median(vector<u_int16_t> scores){
		size_t size = scores.size();

		if (size == 0){
			return 0;  // Undefined, really.
		}
		else{
			sort(scores.begin(), scores.end());
			if (size % 2 == 0){
				return (scores[size / 2 - 1] + scores[size / 2]) / 2;
			}
			else {
				return scores[size / 2];
			}
		}
	}

	double compute_median_float(vector<float> scores){
		size_t size = scores.size();

		if (size == 0){
			return 0;  // Undefined, really.
		}
		else{
			sort(scores.begin(), scores.end());
			if (size % 2 == 0){
				return (scores[size / 2 - 1] + scores[size / 2]) / 2;
			}
			else {
				return scores[size / 2];
			}
		}
	}

	vector<ReadData> _readDatas;
	vector<UnitigData> _unitigDatas;
	//vector<u_int32_t> _nodeAbundances;

	void load_read_compositions(){
		
		//CompositionManager* compositionManager = new CompositionManager(4);
		gzFile file_readComposition = gzopen(_filename_readCompositions.c_str(),"rb");

		while(true){
		
			//cout << readIndex << endl;

			u_int32_t sequenceLength;
			vector<float> composition;
			gzread(file_readComposition, (char*)&sequenceLength, sizeof(sequenceLength));

			if(gzeof(file_readComposition)) break;
			
			composition.resize(_compositionManager->_compositionVectorSize);
			gzread(file_readComposition, (char*)&composition[0], composition.size() * sizeof(float));
			
			_readDatas.push_back({sequenceLength, composition});
			//cout << sequenceLength << endl;
		}
	}
		
	unordered_set<u_int64_t> _filteredMinimizers;

	void load_filtered_minimizers(){
		gzFile file_filteredMinimiers = gzopen(_filename_filteredMinimizers.c_str(),"rb");

		while(true){
			u_int64_t minimizer;
			gzread(file_filteredMinimiers, (char*)&minimizer, sizeof(minimizer));
			if(gzeof(file_filteredMinimiers)) break;
			_filteredMinimizers.insert(minimizer);
			//cout << minimizer << endl;
		}

		gzclose(file_filteredMinimiers);
		cout << "Nb filtered minimizers: " << _filteredMinimizers.size() << endl;
	}

	/*
    void collectNeighbors_sharingReads(AdjGraph* graph, u_int32_t utg, const vector<UnitigData>& unitigDatas, vector<u_int32_t>& neighbors, unordered_set<u_int32_t>& binNodes){

		//results.clear();
		unordered_set<u_int32_t> isVisited;
        neighbors.clear();
        //if(_nodes[s] == nullptr) return;

        //neighbors.push_back(s);

        //for(size_t i=0; i<_nbNodes; i++){
        //    isVisited[i] = false;
        //    distance[i] = 0;
            //prev[i] = -1;
        //}

    
        queue<u_int32_t> queue;
    
        //isVisited[s] = true;
        //distance[s] = 0;
        queue.push(utg);
		isVisited.insert(utg);
		//results.insert(utg_n);
    
    
        while(!queue.empty()){
            
            u_int32_t n = queue.front();
            queue.pop();

            adjNode* node = graph->_nodes[n];

            while (node != nullptr) {

                u_int32_t utg_n = node->val;

                if (isVisited.find(utg_n) != isVisited.end()){
                    node = node->next;
                    continue;
                }

				isVisited.insert(utg_n);

				if(!shareAnyRead(unitigDatas[utg], unitigDatas[utg_n])){
                    node = node->next;
                    continue;
				}
                //isVisited[nn] = true;
                //distance[nn] = distance[n] + 1;

				if(binNodes.find(utg_n) == binNodes.end()){
                	neighbors.push_back(utg_n);
				}


                //if(distance[nn] >= maxDistance) continue;

                //cout << "Push: " << nn << " " << distance[nn] << endl;
                queue.push(utg_n);
				//results.insert(utg_n);

                node = node->next;
            }

        }

    }*/

	/*
	void solveBin(u_int32_t utg, AdjGraph* graph, const vector<UnitigData>& unitigDatas, ofstream& file, int color){

		//cout << computeSharedReads(unitigDatas[924320], unitigDatas[924321]) << endl;
		//cout << computeSharedReads(unitigDatas[924328], unitigDatas[924327]) << endl;
		//cout << computeSharedReads(unitigDatas[924328], unitigDatas[4451136]) << endl;

		unordered_set<u_int32_t> binNodes;

		stack<ReadIndexType> stack;
		stack.push(utg);
		binNodes.insert(utg);
		//unordered_set<u_int32_t> isNodeVisisted;
		unordered_map<DbgEdge, float, hash_pair> cache_distanceTnf;
		unordered_set<DbgEdge, hash_pair> isEdgeVisited;





		while (!stack.empty()){
			
			utg = stack.top();
			stack.pop();


		

			file << utg << "," << color << endl;





			vector<u_int32_t> neighbors;
			//unordered_set<u_int32_t> neighborsShar

			//collectNeighbors_sharingReads(graph, utg, unitigDatas, neighbors, binNodes);
			graph->collectNeighbors(utg, 1, neighbors, 0);
			//if(neighbors.size() == 0) continue;

			//cout << "---------------" << endl;
			//cout << utg << endl;

			for(u_int32_t utg_n : neighbors){
			
			
				DbgEdge edge = {utg, utg_n};
				edge = edge.normalize();

				if(isEdgeVisited.find(edge) != isEdgeVisited.end()){
					continue;
				}

				isEdgeVisited.insert(edge);


				
				vector<float> tnf_distances;
				
				for(u_int32_t r1 : unitigDatas[utg]._readIndexes){
					for(u_int32_t r2 : unitigDatas[utg_n]._readIndexes){
						//if(r1 == r2) continue;
						
						DbgEdge readPair = {r1, r2};
						readPair = readPair.normalize();

						float distance_tnf = -1;
						if(cache_distanceTnf.find(readPair) != cache_distanceTnf.end()){
							//cout << "oyo" << endl;
							distance_tnf = cache_distanceTnf[readPair];
						}
						else{
							distance_tnf = computeDistanceTNF(_readDatas[r1], _readDatas[r2]);
							//distance_tnf = euclidianDistance(_readDatas[r1]._composition, _readDatas[r2]._composition);
							cache_distanceTnf[readPair] = distance_tnf;
						}
						
					
						tnf_distances.push_back(distance_tnf);

					}
					
				}

				float median = compute_median_float(tnf_distances);

				if(median > 0.1){
					//float dist = euclidianDistance(unitigDatas[utg]._compositionMean, unitigDatas[utg_n]._compositionMean);
					cout << utg << " " << utg_n << "      " << median << endl;
				}

				
				u_int64_t nbSharedReads = computeSharedReads(unitigDatas[utg], unitigDatas[utg_n]);
				if(nbSharedReads == 0) continue;
				if(nbSharedReads <= 30){
					//if(median > 0.001){
						continue;
					//}
				}

				stack.push(utg_n);
				binNodes.insert(utg_n);
				cout << "Bin nodes: " << binNodes.size() << " " << utg_n << " " << median << endl;
				


			}





		}

		cout << "Bin nodes: " << binNodes.size() << endl;
	}
	*/


	//void getUnitigLengths(const string& gfa_filename){
	//	u_int64_t nbUnitigs = GfaParser::getUnitigLengths(gfa_filename, _unitigLengths);
	//}

	/*
	void computeNodeAbundance(MDBG* mdbg, const string& gfa_filename){

		string gfa_filename_unitig = gfa_filename + "_tmp.gfa";

		string command = "gfatools asm -u " + gfa_filename + " > " + gfa_filename_unitig;
		
		cout << command << endl;
		int ret = system(command.c_str());
		if(ret != 0){
			cout << "ERROR IN GFA TOOLS" << endl;
			exit(ret);
		}

		vector<int32_t> node_to_unitig(mdbg->_dbg_nodes.size(), -1);
		u_int64_t nbUnitigs = GfaParser::getNodeToUnitig(gfa_filename_unitig, node_to_unitig);

		cout << nbUnitigs << endl;
		
		_nodeAbundances.resize(mdbg->_dbg_nodes.size(), 0);
		vector<u_int32_t> unitigAbundances(nbUnitigs, 0);
		vector<u_int32_t> unitigAbundancesSize(nbUnitigs, 0);

		for(auto it : mdbg->_dbg_nodes){

			const KmerVec& vec = it.first;

			u_int32_t kminmer_index = mdbg->_dbg_nodes[vec]._index;
			//u_int32_t unitigIndex = node_to_unitig[kminmer_index];
			//if(unitigIndex == -1) continue;

			//unitigDatas[unitigIndex]._nbKminmers += 1;
			//unitigDatas[unitigIndex]._meanAbundance += mdbg->_dbg_nodes[vec]._abundance;
			u_int32_t unitigIndex = node_to_unitig[kminmer_index];
			//if(unitigIndex == -1) continue;

			if(mdbg->_dbg_nodes[vec]._abundance > unitigAbundances[unitigIndex]){
				unitigAbundances[unitigIndex] = mdbg->_dbg_nodes[vec]._abundance;
			}
			//unitigAbundances[unitigIndex] = max(unitigAbundances[unitigIndex], mdbg->_dbg_nodes[vec]._abundance);
			//unitigAbundancesSize[unitigIndex] += 1;

			//_unitigDatas[kminmer_index]._meanAbundance += mdbg->_dbg_nodes[vec]._abundance;
			//_unitigDatas[kminmer_index]._nbKminmers += 1;
		}
		
		for(auto it : mdbg->_dbg_nodes){
			const KmerVec& vec = it.first;
			u_int32_t kminmer_index = mdbg->_dbg_nodes[vec]._index;

			u_int32_t unitigIndex = node_to_unitig[kminmer_index];
			//_unitigDatas[unitigIndex]._meanAbundance = unitigAbundances[unitigIndex];
			_nodeAbundances[kminmer_index] = unitigAbundances[unitigIndex];
			//if(unitigIndex == -1) continue;

			//_unitigDatas[kminmer_index]._meanAbundance = unitigAbundances[unitigIndex]; // ((float) unitigAbundances[unitigIndex]) / unitigAbundancesSize[unitigIndex];
		}

		std::remove(gfa_filename_unitig.c_str());

		
		cout << mdbg->_dbg_nodes.size() << endl;
		float abundanceCutoff_min = 4.2;
		unordered_set<u_int32_t> nodes;
		for(auto it : mdbg->_dbg_nodes){
			const KmerVec& vec = it.first;
			u_int32_t kminmer_index = mdbg->_dbg_nodes[vec]._index;

			u_int32_t unitigIndex = node_to_unitig[kminmer_index];
			//if(unitigIndex == -1) continue;

			//_unitigDatas[kminmer_index]._meanAbundance = unitigAbundances[unitigIndex]; // ((float) unitigAbundances[unitigIndex]) / unitigAbundancesSize[unitigIndex];
			//if(_unitigDatas[unitigIndex]._meanAbundance > abundanceCutoff_min){
			if(_nodeAbundances[kminmer_index] > abundanceCutoff_min){
				nodes.insert(kminmer_index);
			}
		}
		GfaParser::rewriteGfa_withNodes(gfa_filename, gfa_filename + "_debug_filterByAbundance.gfa", nodes);
		cout << nodes.size() << endl;
		
		//size_t i=0;
		//for(UnitigData& data : _unitigDatas){
		//	cout << i << " " << data._meanAbundance << endl;
		//	i += 1;
		//}
	}*/




	struct PathData{
		u_int32_t _index;
		unordered_set<DbgEdge, hash_pair> isEdgeVisited;
		vector<u_int32_t> prevNodes;
		u_int32_t source_abundance;
		u_int32_t source_nodeIndex;
		u_int32_t source_nodeIndex_path;
		float _abundanceCutoff_min;
	};

	//vector<PathData> _pathDatas;
	//vector<int32_t> _node_to_unitig;

	//unordered_map<u_int32_t, u_int32_t> _nodeLabel;
	//unordered_set<u_int32_t> _globalVisitedNodes;

	static bool SuccessorComparator_byPrevRank(const SuccessorData &a, const SuccessorData &b){
		return a._prevRank > b._prevRank;
	}

	static bool SuccessorComparator_byAbundance(const SuccessorData &a, const SuccessorData &b){
		return a._sourceAbundance < b._sourceAbundance;
	}

	static bool UnitigComparator_ByLength(const UnitigLength &a, const UnitigLength &b){
		return a._length > b._length;
	}

	size_t _iter;
	ofstream file_groundTruth;

	/*
	void getSuccessors(u_int32_t nodeIndex, const PathData& pathData, BiGraph* graph, vector<u_int32_t>& successors){

		bool orient_dummy = false;
		successors.clear();

		adjNode* node = graph->_nodes[nodeIndex];

		while(node != nullptr){

			u_int64_t utg_n = node->val;
			u_int32_t nodeName = graph->nodeIndex_to_nodeName(utg_n, orient_dummy);
			//cout << unitigIndex << " " << _node_to_unitig[unitigIndex] << endl;
			
			u_int32_t unitigIndex = _node_to_unitig[nodeName];

			//cout << _nodeAbundances[nodeName] << endl;
			if(unitigIndex == -1){ //Cleaned
				node = node->next;
				continue;
			}

			if(_nodeAbundances[nodeName] < pathData._abundanceCutoff_min){ //Abundance min cutoff
				node = node->next;
				continue;
			}

			//cout << _nodeAbundances[nodeName] << " " <<  pathData._abundanceCutoff_min << endl;
			successors.push_back(node->val);

			node = node->next;
		}

	}
	*/

	unordered_set<u_int32_t> _binnedNodes;
	unordered_map<u_int32_t, u_int32_t> _nbVisitedTimes;

	void binNode(u_int32_t nodeIndex, vector<u_int32_t>& prevNodes, GraphSimplify* graph, u_int32_t pathIndex){

		bool orient_dummy;
		//u_int32_t utg_nodeIndex = nodeIndex; //successors[0]._nodeIndex;
		prevNodes.push_back(nodeIndex);

		//_nodeLabel[nodeIndex] = pathIndex;
		u_int32_t current_unitigIndex = graph->_graphSuccessors->nodeIndex_to_nodeName(nodeIndex, orient_dummy);
		file_groundTruth << current_unitigIndex << "," << pathIndex << endl;

		_binnedNodes.insert(current_unitigIndex);
		//cout << "Add node: " << graph->_graphSuccessors->nodeToString(nodeIndex) << " " << graph->_nodeAbundances[current_unitigIndex]  << endl;

		_nbVisitedTimes[current_unitigIndex] += 1;
		//cout << _nbVisitedTimes[current_unitigIndex] << endl;
		//return utg_nodeIndex;
	}

	float computeAbundanceCutoff_min(u_int32_t abundance){
		return abundance / 4.0;
	}

	void solveBin(u_int32_t source_nodeIndex, u_int32_t abundance, GraphSimplify* graph, int pathIndex){


		bool orient_dummy = false;
		u_int32_t nodeName = graph->_graphSuccessors->nodeIndex_to_nodeName(source_nodeIndex, orient_dummy);

		cout << "----- Start solve bin ------" << endl;
		cout << "Source: " << nodeName << " " << graph->_graphSuccessors->nodeToString(source_nodeIndex) << endl;

		//cout << graph->_nodeToUnitig[source_nodeIndex] << endl;
		//cout << graph->_unitigs[graph->_nodeToUnitig[source_nodeIndex]]._abundance << endl;
		//bool orient_dummy;
		//u_int32_t utg_nodeIndex = source_nodeIndex; //graph->nodeName_to_nodeIndex(utg, orient_dummy);
		
		//u_int32_t utg_nodeIndex = graph->nodeName_to_nodeIndex(utg, orient_dummy);
		

		//float abundanceCutoff_min = graph->_nodeAbundances[nodeName] / 5.0;
		float abundanceCutoff_min = computeAbundanceCutoff_min(abundance);
		cout << "Abundance cutoff min: " << abundanceCutoff_min << endl;
		//if(abundanceCutoff_min < 30) return;

		cout << "Simplifying graph local" << endl;
		graph->clear();
		graph->debug_writeGfaErrorfree(abundanceCutoff_min);
		//graph->clear();
		//graph->execute(abundanceCutoff_min);
		//graph = graph->clone();
		//cout << "clone done" << endl;
		//graph->execute(abundanceCutoff_min);


		//u_int32_t source_unitigIndex = graph->nodeIndex_to_nodeName(utg, orient_dummy);
		if(graph->_nodeToUnitig[source_nodeIndex] == -1){
			cout << "Source node removed :(" << endl;
			return; //????
		}
		u_int32_t source_abundance = graph->_unitigs[graph->_nodeToUnitig[source_nodeIndex]]._abundance;
		cout << graph->_unitigs[graph->_nodeToUnitig[source_nodeIndex]]._abundance << endl;
		//cout << "PAS BON CA: utiliser successeur puis predeesseur" << endl;
		cout << "----- Forward ------" << endl;
		_iter = 0;
		//u_int32_t source_nodeIndex = graph->_graphPredecessors->nodeName_to_nodeIndex(utg, false);
		PathData pathData = {pathIndex, {}, {}, source_abundance, source_nodeIndex, source_nodeIndex, abundanceCutoff_min};
		bool pathSolved = solveBin_path(pathData, graph, true);
		
		if(pathSolved) return;

		cout << "----- Backward ------" << endl;
		_iter = 0;
		pathData = {pathIndex, {}, {}, source_abundance, source_nodeIndex, source_nodeIndex, abundanceCutoff_min};
		solveBin_path(pathData, graph, false);

		/*
		cout << "----- Backward ------" << endl;
		_iter = 0;
		//u_int32_t source_nodeIndex = graph->_graphPredecessors->nodeName_to_nodeIndex(utg, false);
		PathData pathData = {pathIndex, {}, {}, source_abundance, source_nodeIndex, source_nodeIndex, abundanceCutoff_min};
		solveBin_path(pathData, graph);
		*/

		/*
		cout << "----- Start extending left ------" << endl;
		_iter = 0;
		source_nodeIndex = graph->nodeName_to_nodeIndex(utg, true);
		pathData = {pathIndex, {}, {}, source_abundance, source_nodeIndex, source_nodeIndex, abundanceCutoff_min};
		solveBin_path(pathData, graph);
		//solveBin_path(pathData, graph, false);
		*/
		//cout << _pathDatas.size() << endl;
	}

	bool solveBin_path(PathData& pathData, GraphSimplify* graph, bool forward){

		unordered_set<u_int32_t> visitedNodes;

		u_int32_t current_nodeIndex = pathData.source_nodeIndex;
		binNode(current_nodeIndex, pathData.prevNodes, graph, pathData._index);

		while(true){
			PathExplorer pathExplorer(pathData.prevNodes, pathData.source_abundance, pathData.source_nodeIndex, current_nodeIndex, pathData._abundanceCutoff_min, visitedNodes);
			//u_int32_t nextNodeIndex = pathExplorer.getNextNode( graph, _unitigDatas, 100000, forward, true);
			current_nodeIndex = pathExplorer.getNextNode(current_nodeIndex, graph, _unitigDatas, forward, 0);
			//cout << graph->_graphSuccessors->nodeToString(current_nodeIndex) << endl;

			if(current_nodeIndex == -1) return false; //No more successors, or no branching solution
			
			if(current_nodeIndex == pathData.source_nodeIndex){ //Path complete
				cout << "Path complete!" << endl;
				return true; 
			}
			//if(current_nodeIndex == -2) return true; //Path complete

			binNode(current_nodeIndex, pathData.prevNodes, graph, pathData._index);
			visitedNodes.insert(current_nodeIndex);

			/*
			for(size_t i=1; i<pathExplorer._exploredNodes.size(); i++){
				u_int32_t nodeIndex = pathExplorer._exploredNodes[i];
				current_nodeIndex = nodeIndex;
			}

			if(!foundPath) break;
			*/
		}
	}



	

};	

/********************************************************************************/

#endif /* _BLOOCOO_HPP_ */



/*
void solveBin(u_int32_t utg, AdjGraph* graph, const vector<UnitigData>& unitigDatas, ofstream& file, int color){



		unordered_set<u_int32_t> binNodes;

		stack<ReadIndexType> stack;
		stack.push(utg);
		binNodes.insert(utg);
		unordered_set<u_int32_t> isNodeVisisted;
		unordered_map<DbgEdge, float, hash_pair> cache_distanceTnf;




		while (!stack.empty()){
			
			utg = stack.top();
			stack.pop();


			file << GfaParser::unitigIndex_to_unitigName(utg) << "," << color << endl;





			vector<u_int32_t> neighbors;
			//unordered_set<u_int32_t> neighborsShar

			//collectNeighbors_sharingReads(graph, utg, unitigDatas, neighbors, binNodes);
			graph->collectNeighbors(utg, 1, neighbors);
			//if(neighbors.size() == 0) continue;

			cout << "---------------" << endl;
			cout << GfaParser::unitigIndex_to_unitigName(utg) << endl;

			for(u_int32_t utg_n : neighbors){
			
				if(isNodeVisisted.find(utg_n) != isNodeVisisted.end()){
					continue;
				}

				isNodeVisisted.insert(utg_n);


				float dist = euclidianDistance(unitigDatas[utg]._compositionMean, unitigDatas[utg_n]._compositionMean);
				cout << GfaParser::unitigIndex_to_unitigName(utg) << " " << GfaParser::unitigIndex_to_unitigName(utg_n) << " " << dist << endl;
				
				//if(binNodes.find(utg_n) != binNodes.end()){
					//node = node->next;
				//	continue;
				//}
				
				vector<float> tnf_distances;
				
				for(u_int32_t r1 : unitigDatas[utg]._readIndexes){
					for(u_int32_t r2 : unitigDatas[utg_n]._readIndexes){
						//if(r1 == r2) continue;
						
						DbgEdge readPair = {r1, r2};
						readPair = readPair.normalize();

						float distance_tnf = -1;
						if(cache_distanceTnf.find(readPair) != cache_distanceTnf.end()){
							//cout << "oyo" << endl;
							distance_tnf = cache_distanceTnf[readPair];
						}
						else{
							distance_tnf = computeDistanceTNF(_readDatas[r1], _readDatas[r2]);
							cache_distanceTnf[readPair] = distance_tnf;
						}
						
					
						tnf_distances.push_back(distance_tnf);

					}
					
				}

				float median = compute_median_float(tnf_distances);
				cout << median << endl;

				u_int64_t nbSharedReads = computeSharedReads(unitigDatas[utg], unitigDatas[utg_n]);
				if(nbSharedReads <= 1) continue;

				stack.push(utg_n);
				binNodes.insert(utg_n);
				cout << "Bin nodes: " << binNodes.size() << " " << GfaParser::unitigIndex_to_unitigName(utg_n) << endl;
				


			}

		}

		cout << "Bin nodes: " << binNodes.size() << endl;
	}
*/
