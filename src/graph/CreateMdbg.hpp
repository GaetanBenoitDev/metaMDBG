
#ifndef MDBG_METAG_CREATEMDBG
#define MDBG_METAG_CREATEMDBG

#include "Commons.hpp"
#include "graph/ProgressiveAbundanceFilter.hpp"

//#include "graph/Graph.hpp"
//#include "GraphSimplify.hpp"
//#include "../utils/BloomFilter.hpp"
//#include "./utils/ntHashIterator.hpp"

//#include <boost/graph/adjacency_list.hpp>
//#include <boost/graph/graph_utility.hpp>
//#include <boost/property_map/property_map.hpp>
//#include <boost/graph/connected_components.hpp>
//#include <boost/graph/dijkstra_shortest_paths.hpp>




//struct KminmerAbundanceMapNode{
//	u_int32_t _abundance;
//	u_int32_t _contigIndex;
//};

struct DbgNodeLight{
	bool _isWritten;
	u_int32_t _index;
	u_int32_t _abundance;
};

//typedef phmap::parallel_flat_hash_map<u_int128_t, u_int32_t> KminmerAbundanceMap;
typedef phmap::parallel_flat_hash_map<u_int128_t, DbgNodeLight, phmap::priv::hash_default_hash<u_int128_t>, phmap::priv::hash_default_eq<u_int128_t>, std::allocator<std::pair<u_int128_t, DbgNodeLight>>, 10, std::mutex> MdbgNodeMapLight;
typedef phmap::parallel_flat_hash_map<u_int128_t, u_int32_t, phmap::priv::hash_default_hash<u_int128_t>, phmap::priv::hash_default_eq<u_int128_t>, std::allocator<std::pair<u_int128_t, u_int32_t>>, 10, std::mutex> KminmerAbundanceMap;

//typedef phmap::parallel_flat_hash_map<KmerVec, u_int32_t, phmap::priv::hash_default_hash<KmerVec>, phmap::priv::hash_default_eq<KmerVec>, std::allocator<std::pair<KmerVec, u_int32_t>>, 10, std::mutex> UnitigEdgeMap;

struct UnitigEdge{
	UnitigType _unitigIndexFrom;
	UnitigType _unitigIndexTo;
};

class CreateMdbg : public Tool{
    
	//struct KminmerData{
	//	u_int32_t _count;
	//	u_int16_t _length;
	//	u_int16_t _overlapLength_start;
	//	u_int16_t _overlapLength_end;
	//	bool _isReversed;
	//};

	//struct KminmerAbundance{
	//	KmerVec _vec;
	//	u_int32_t _count;
	//};


	/*
	class Parallel_MdbgEdgeMap3{

		public: 
		
		int _nbPartitions;
		vector<MdbgEdgeMap3> _maps;
		vector<omp_lock_t> _mutex;

		void setup(){

			_nbPartitions = 32;

			_maps.resize(_nbPartitions);
			_mutex.resize(_nbPartitions);

			for(size_t i=0; i<_nbPartitions; i++){
				omp_init_lock(&_mutex[i]);
			}
		}

		void clear(){

			for(size_t i=0; i<_nbPartitions; i++){
				omp_destroy_lock(&_mutex[i]);
			}
		}

		void insert(const KmerVec& edgeVec, const u_int128_t& key, const KminmerEdge3& value){

			int partition = key % _nbPartitions;

			omp_set_lock(&_mutex[partition]);

			#pragma omp critical
			{
				cout << partition << endl;
			}

			cout << partition << endl;
			MdbgEdgeMap3& map = _maps[partition];

			if(!successorExists(map, edgeVec, key, value._isReversed, value._isPrefix)){
				map[key].push_back(value);
				if(map[key].size() > 2) cout << map[key].size() << " " << map.size() << endl;
			}

			omp_unset_lock(&_mutex[partition]);
		}

		bool get(const u_int128_t& key, KminmerEdge3& value){

			int partition = key % _nbPartitions;

			const auto& map = _maps[partition];

			if(map.find(key) == map.end()) return false;
			
			//value = map[key];
			return true;

		}


		bool successorExists(MdbgEdgeMap3& map, const KmerVec& edgeVec, const u_int128_t& key, const bool isReversed, const bool isPrefix){

			//return false;
			if(map.find(key) == map.end()) return false;

			bool isPalindrome = edgeVec.isPalindrome();

			for(KminmerEdge3& edge : map[key]){
				
				if(isPalindrome){
					edge._hasMultipleSuccessors = true;
					return true;
				}

				if(edge._isReversed == isReversed && edge._isPrefix == isPrefix){
					edge._hasMultipleSuccessors = true;
					return true;
				}

				if(edge._isReversed == !isReversed && edge._isPrefix == !isPrefix){
					edge._hasMultipleSuccessors = true;
					return true;
				}

			}

			return false;


		}

	};*/


public:

	/*
	//string _inputFilename;
	string _outputDir;
	float _minimizerDensity;
    size_t _minimizerSize;
    size_t _kminmerSize;
    size_t _kminmerSizeFirst;
    size_t _kminmerSizePrev;
	size_t _kminmerSizeLast;
	string _filename_noKminmerReads;
	//ofstream _file_noKminmerReads;
	bool _useBloomFilter;
	bool _useCorrectedRead;
	u_int32_t _unitigIndex;
	u_int64_t _nbKminmersTotal;
	//unordered_set<u_int128_t> _isKminmerIndexed;
	//UnitigEdgeMap _startNode_to_unitigIndex;
	//UnitigEdgeMap _endNode_to_unitigIndex;
	//unordered_map<KmerVec, u_int32_t> _startNode_to_unitigIndex;
	//unordered_map<KmerVec, u_int32_t> _endNode_to_unitigIndex;
	
	//BloomCacheCoherent<u_int64_t>* _bloomFilter;
	//unordered_set<KmerVec> _bloomFilterExact;
    //IBank* _inputBank;

	//vector<ReadData> _readData;
	//CompositionManager* _compositionManager;
	//AdjGraph* _overlapGraph;
	//string _inputDir;
	//string _input_extractKminmers;
	//KmerVec _lala;

	//u_int64_t _nbReads;
	//string _filename_readMinimizers;
	//string _filename_contigMinimizers;
	//string _filename_inputContigs;
	//string _filename_solidKminmers;
	*/
	string _outputDir;
	bool _isFirstPass;
	Parameters _params;
	ReadStats _readStats;

    size_t _kminmerSize;
    size_t _kminmerSizeFirst;
    size_t _kminmerSizePrev;
	int _minAbundance;
	int _nbCores;
	UnitigType _unitigIndex;
	u_int64_t _nbKminmersTotal;
	//size_t _kminmerSizeLast;

	string _filename_smallContigs;
	//string _filename_smallContigsNew;
	//string _filename_smallContigsPrev;
	//string _filename_smallContigs;

	//unordered_set<u_int32_t> writtenNodeNames;
	ofstream _kminmerFile;
	ofstream _kminmerAbundanceFile;
	//ofstream _kminmerFile2;
	//ofstream _readFile;
	ofstream _fileSmallContigs;
	//string _filename_filteredMinimizers;
	//string _filename_readCompositions;

	//vector<u_int32_t> _evaluation_readToDataset;
	//bool _parseReads;
	//ankerl::unordered_dense::map<u_int128_t, u_int32_t> _mdbgNodesLight;
	KminmerAbundanceMap _mdbgNodesLight;
	//MdbgNodeMapLight _mdbgNodesLightUnique;
	//unordered_map<KmerVec, u_int32_t> _kmervec_to_nodeName;
	//unordered_set<u_int128_t> _chimericKminmers;
	//MDBG* _mdbg;
	//MDBG* _mdbgFirst;
	//MDBG* _mdbgSaved;
	//MDBG* _mdbgInit;
	KminmerAbundanceMap _kminmerAbundances;
	//MdbgEdgeMap _mdbgEdges;
	//MdbgEdgeMap2 _mdbgEdges2;
	//Parallel_MdbgEdgeMap3 _mdbgEdges5;
	MdbgEdgeMap3 _mdbgEdges3;
	MdbgEdgeMap4 _mdbgEdges4;
	//MdbgEdgeMapTest _mdbgEdgesTest;
	//MdbgEdgeMap6 _mdbgEdges5;
	//MDBG* _mdbgNoFilter;
	//bool _parsingContigs;
	//bool _savingSmallContigs;
	//u_int32_t _node_id;
	//MinimizerPairMap* _minimizerPairMap;

	float _minimizerSpacingMean;
	float _kminmerLengthMean;
	float _kminmerOverlapMean;

	double _minimizerSpacing_sum;
	double _minimizerSpacing_n;
	double _kminmerLength_sum;
	double _kminmerLength_n;
	u_int64_t _nbEdges;
	int _nbPartitions;

	//int _lalaLol;
	//vector<MinimizerType> _lalaLol2;
	//unordered_map<KmerVec, u_int32_t> _kmerVec_to_nodeName;
	vector<vector<MinimizerType>> _unitigName_to_minimizers;
	MdbgNodeMapLight _solidTriplets;

	u_int64_t _checksum_unitigNodes;
	u_int64_t _checksum_unitigEdges;
	u_int64_t _checksum_unitigAbundances;

	typedef boomphf::SingleHashFunctor<u_int128_t>  hasher_t;
	typedef boomphf::mphf<u_int128_t, hasher_t> boophf_t;

	u_int64_t _nbUnitigNodes;
	u_int64_t _nbUnitigEdges;
	vector<omp_lock_t> _mutexes;
	/*
	// iterator from disk file of uint64_t with buffered read,   todo template
	template <typename basetype>
	class bfile_iterator : public std::iterator<std::forward_iterator_tag, basetype>{
	public:
		
		bfile_iterator()
		: _is(nullptr)
		, _pos(0) ,_inbuff (0), _cptread(0)
		{
			_buffsize = 10000;
			_buffer = (basetype *) malloc(_buffsize*sizeof(basetype));
		}
		
		bfile_iterator(const bfile_iterator& cr)
		{
			_buffsize = cr._buffsize;
			_pos = cr._pos;
			_is = cr._is;
			_buffer = (basetype *) malloc(_buffsize*sizeof(basetype));
			 memcpy(_buffer,cr._buffer,_buffsize*sizeof(basetype) );
			_inbuff = cr._inbuff;
			_cptread = cr._cptread;
			_elem = cr._elem;
		}
		
		bfile_iterator(FILE* is): _is(is) , _pos(0) ,_inbuff (0), _cptread(0)
		{
			//printf("bf it %p\n",_is);
			_buffsize = 10000;
			_buffer = (basetype *) malloc(_buffsize*sizeof(basetype));
			int reso = fseek(_is,0,SEEK_SET);
			advance();
		}
		
		~bfile_iterator()
		{
			if(_buffer!=NULL)
				free(_buffer);
		}
		
		
		basetype const& operator*()  {  return _elem;  }
		
		bfile_iterator& operator++()
		{
			advance();
			return *this;
		}
		
		friend bool operator==(bfile_iterator const& lhs, bfile_iterator const& rhs)
		{
			if (!lhs._is || !rhs._is)  {  if (!lhs._is && !rhs._is) {  return true; } else {  return false;  } }
			assert(lhs._is == rhs._is);
			return rhs._pos == lhs._pos;
		}
		
		friend bool operator!=(bfile_iterator const& lhs, bfile_iterator const& rhs)  {  return !(lhs == rhs);  }
	private:
		void advance()
		{
			
			//printf("_cptread %i _inbuff %i \n",_cptread,_inbuff);
			
			_pos++;
			
			if(_cptread >= _inbuff)
			{

				int res = fread(_buffer,sizeof(basetype),_buffsize,_is);

				//printf("read %i new elem last %llu  %p\n",res,_buffer[res-1],_is);
				_inbuff = res; _cptread = 0;
				
				if(res == 0)
				{
					_is = nullptr;
					_pos = 0;
					return;
				}
			}
			
			_elem = _buffer[_cptread];
			_cptread ++;
		}
		basetype _elem;
		FILE * _is;
		unsigned long _pos;
		
		basetype * _buffer; // for buffered read
		int _inbuff, _cptread;
		int _buffsize;
	};
	
	
	template <typename type_elem>
	class file_binary{
	public:
		
		file_binary(const char* filename)
		{
			_is = fopen(filename, "rb");

			if (!_is) {
				throw std::invalid_argument("Error opening " + std::string(filename));
			}
		}
		
		~file_binary()
		{
			fclose(_is);
		}
		
		bfile_iterator<type_elem> begin() const
		{
			return bfile_iterator<type_elem>(_is);
		}
		
		bfile_iterator<type_elem> end() const {return bfile_iterator<type_elem>(); }
		
		size_t        size () const  {  return 0;  }//todo ?
		
	private:
		FILE * _is;
	};
	*/




	//ofstream _debugFile;
    CreateMdbg ();
    void execute ();
	void createMDBG();
	//void createMDBG_read(kseq_t* read, u_int64_t readIndex);
	void createMDBG_collectKminmers(const vector<MinimizerType>& minimizers, const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex);
	void createMDBG_collectKminmers_contig(const vector<MinimizerType>& minimizers, const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex);
	void removeErroneousKminmers(const vector<MinimizerType>& minimizers, const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex);
	void loadSolidKminmers();
	void dumpNodes();

	void createMDBG_collectKminmers_read(kseq_t* read, u_int64_t readIndex);
	void extractKminmerSequence(const char* sequenceOriginal, const ReadKminmer& kminmerInfo, string& sequence);
	//void computeDeterministicNodeNames();
	void indexEdgeUnitig(const KmerVec& vec, UnitigType unitigIndex, bool isStartNode);
	void computeNextUnitigGraph();
	void removeChimericUnitigs(UnitigGraph2* unitigGraph);
	bool isUnsupportedUnitig(const vector<MinimizerType>& minimizers);
	void solveEdges(UnitigGraph2* unitigGraph);
	void solveEdge(UnitigGraph2* unitigGraph, const UnitigType unitigName, const UnitigType unitigIndex, unordered_set<UnitigType>& processedUnitigNames);
	void createDoubletNode(UnitigGraph2* unitigGraph, const UnitigType unitigIndex, const UnitigType unitigIndexSuccessor, const vector<MinimizerType>& minimizers, unordered_set<UnitigType>& processedUnitigNames);
	void solveSmallUnitigs(UnitigGraph2* unitigGraph);
	void solveSmallUnitigsSub(UnitigGraph2* unitigGraph, UnitigGraph2::UnitigNode* node);
	void solveSmallUnitigsSub2(UnitigGraph2* unitigGraph, UnitigGraph2::UnitigNode* node);
	void removeUnsupportedUnitigs(UnitigGraph2* unitigGraph);
	void createEdgeNodePredecessor(UnitigGraph2* unitigGraph, const UnitigType unitigIndex, const UnitigType unitigIndexPredecessor, phmap::parallel_flat_hash_map<pair<UnitigType, UnitigType>, UnitigGraph2::UnitigNode*>& edgeNodes);
	void createEdgeNodeSuccessor(UnitigGraph2* unitigGraph, const UnitigType unitigIndex, const UnitigType unitigIndexSuccessor, phmap::parallel_flat_hash_map<pair<UnitigType, UnitigType>, UnitigGraph2::UnitigNode*>& edgeNodes);	
	UnitigGraph2::UnitigNode* createEdgeNode(UnitigGraph2* unitigGraph, const vector<MinimizerType>& minimizers);
	void writeUnitigs(UnitigGraph2* unitigGraph);
	void extendNode(UnitigGraph2* unitigGraph, const UnitigType unitigIndex, bool checkOnly);
	void extendNodeBranching(UnitigGraph2* unitigGraph, const UnitigType unitigIndex, bool checkOnly);
	void computeUnitigSequences(UnitigGraph2* unitigGraph, const string outputFilename);

	void createMDBG_index(const vector<MinimizerType>& minimizers, const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex);
	void createGfa();
	void parseArgs(int argc, char* argv[]);
	void parseContigs();
	float computeKmerVecAbundance(const vector<MinimizerType>& minimizers, bool isContig);
	void computeContigAbundance();
	void indexEdge(const KmerVec& vec, u_int32_t nodeName);
	void indexEdges();
	void indexUnitigEdges();
	void indexUnitigEdge(const UnitigType unitigIndex, const vector<MinimizerType>& minimizers);
	void computeUnitigEdges();
	void computeUnitigEdge(const UnitigType& unitigIndex, const vector<MinimizerType>& minimizers);
	void computeUnitigNodes();
	//void computeUnitigNode(const KmerVec& vec);
	void getSuccessors(const KmerVec& vec, vector<KmerVec>& successors);
	void getPredecessors(const KmerVec& vec, vector<KmerVec>& predecessors);
	bool hasSingleSuccessor(const KmerVec& vec);
	bool hasSinglePredecessor(const KmerVec& vec);
	void getSuccessors_unitig(const KmerVec& vec, const UnitigType& unitigIndexFrom, vector<UnitigType>& successors);
	void getPredecessors_unitig(const KmerVec& vec, const UnitigType& unitigIndexFrom, vector<UnitigType>& predecessors);
	void dumpUnitigNode(const UnitigType& unitigIndex, const vector<MinimizerType>& unitig, ofstream& outputFile);
	void dumpUnitigEdge(const UnitigType& unitigIndexFrom, const vector<UnitigType>& successors, const vector<UnitigType>& predecessors);

	//bool isBranchingNode(const KmerVec& vec);
	//void setNodeUnitigged(const KmerVec& vec);
	//bool isNodeUnitigged(const KmerVec& vec);
	void dumpUnitigAbundances();
	void loadRefinedAbundances();
	void computeDeterministicUnitigs();
	bool successorExists(const KmerVec& edgeVec, const u_int128_t& key, const bool isReversed, const bool isPrefix);
	//int getNbSuccessors(const KmerVec& vec);
	//int getNbPredecessors(const KmerVec& vec);
	bool getNbSuccessors(const KmerVec& vec, KmerVec& singleSuccessor);
	bool getNbPredecessors(const KmerVec& vec, KmerVec& singlePredecessor);
	void writeUnitigGraphStat(const u_int64_t nbUnitigNodes, const u_int64_t nbUnitigEdges);

	//u_int64_t _edgeIndexElements;
	ofstream _unitigGraphFile_nodes;
	ofstream _unitigGraphFile_nodes_abundances;
	ofstream _unitigGraphFile_edges_successors;
	//ofstream _kminmerFile;
	/*
	vector<MinimizerType> normalizeUnitig(const vector<MinimizerType>& minimizers) const{

		vector<MinimizerType> minimizers_reversed = minimizers;
		std::reverse(minimizers_reversed.begin(), minimizers_reversed.end());
		
		for(size_t i=0; i<minimizers.size(); i++){
			if(minimizers[i] == minimizers_reversed[i]){
				continue;
			}
			else if(minimizers[i] < minimizers_reversed[i]){
				return minimizers;
				//isReversed = false;
				//return *this;
			}
			else{
				return minimizers_reversed;
				//isReversed = true;
				//return vec_reverse;
			}
		}

		//isReversed = true;
		//return vec_reverse;
		return minimizers_reversed;

		//if(h() < vec_reverse.h()){
		//	isReversed = false;
		//	return *this;
		//}
		//else{
		//	isReversed = true;
		//	return vec_reverse;
		//}
		
	}
	*/
	/*
	struct Contig{
		u_int64_t _readIndex;
		vector<u_int32_t> _nodepath;
		vector<u_int64_t> _minimizers;
		//vector<vector<u_int64_t>> _kmerVecs;
		//vector<u_int32_t> _nodepath_sorted;
	};

	vector<Contig> _contigs_current;
	vector<Contig> _contigs_prev;

	static bool ContigComparator_ByLength(const Contig &a, const Contig &b){

		if(a._nodepath.size() == b._nodepath.size()){
			for(size_t i=0; i<a._nodepath.size() && i<b._nodepath.size(); i++){
				if(BiGraph::nodeIndex_to_nodeName(a._nodepath[i]) == BiGraph::nodeIndex_to_nodeName(b._nodepath[i])){
					continue;
				}
				else{
					return BiGraph::nodeIndex_to_nodeName(a._nodepath[i]) > BiGraph::nodeIndex_to_nodeName(b._nodepath[i]);
				}
			}
		}


		return a._nodepath.size() > b._nodepath.size();
	}
	*/
	/*
	//phmap::parallel_flat_hash_map<u_int32_t, vector<u_int32_t>> _nodeName_to_contigs;
	//phmap::parallel_flat_hash_set<u_int32_t> _invalidContigIndex;
	unordered_set<u_int32_t> isUnitigDuplicated;

	void checkDuplication(){


		ofstream bannedKminmers(_outputDir + "/bannedKminmers/bannedKminmers_" + to_string(_kminmerSizePrev) + ".txt");
		cout << "Removing duplication" << endl;
		
		ofstream colorFile(_outputDir + "/color.csv");
		colorFile << "Name,Color" << endl;

		cout << "wirte: " << _kminmerSizePrev << endl;
		bannedKminmers.write((const char*)&_kminmerSizePrev, sizeof(_kminmerSizePrev));

		//KminmerParserParallel parser(_filename_output, _minimizerSize, _kminmerSizeFirst, false, false);
		//auto fp = std::bind(&ToMinspace::loadContigs_min_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
		//parser.parseMinspace(fp);

		KminmerParserParallel parser(_outputDir + "/unitig_data.txt", _minimizerSize, _kminmerSizeFirst, false, false, _nbCores);
		parser.parse(LoadContigKminmerFunctor(*this, _contigs_current));

		KminmerParserParallel parser2(_outputDir + "/unitig_data_prev.txt", _minimizerSize, _kminmerSizeFirst, false, false, _nbCores);
		parser2.parse(LoadContigKminmerFunctor(*this, _contigs_prev));

		for(size_t i=0; i<_contigs_current.size(); i++){

			const Contig& contigCurrent = _contigs_current[i];

			for(size_t j=0; j<_contigs_prev.size(); j++){

				const Contig& contigPrev = _contigs_prev[j];

				if(contigCurrent._nodepath == contigPrev._nodepath){

					isUnitigDuplicated.insert(contigCurrent._readIndex);
					//for(const vector<u_int64_t>& vec : contigCurrent._kmerVecs){
					//	bannedKminmers.write((const char*)&vec[0], vec.size()*sizeof(uint64_t));
					//}

					cout << "duplicate " << contigCurrent._nodepath.size() << endl;

					vector<u_int64_t> rlePositions;
					vector<u_int64_t> minimizers_pos;//(minimizers.size());
					vector<KmerVec> kminmers; 
					vector<ReadKminmer> kminmersInfo;
					MDBG::getKminmers(_minimizerSize, _kminmerSizePrev, contigCurrent._minimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, 0, false);

					for(const KmerVec& vec : kminmers){
						//bannedKminmers.write((const char*)&vec._kmers[0], vec._kmers.size()*sizeof(uint64_t));
						u_int32_t nodeName = _kmerVec_to_nodeName[vec];
						colorFile << nodeName << ",red" << endl;
					}

					//if(contigCurrent._nodepath.size() > 500){
					//	cout << "duplicate " << contigCurrent._nodepath.size() << endl;
					//}
				}
			}
		}

		bannedKminmers.close();
		colorFile.close();
		//getchar();
	}

	class LoadContigKminmerFunctor {

		public:

		CreateMdbg& _graph;
		vector<Contig>& _contigs;

		LoadContigKminmerFunctor(CreateMdbg& graph, vector<Contig>& contigs) : _graph(graph), _contigs(contigs){

		}

		LoadContigKminmerFunctor(const LoadContigKminmerFunctor& copy) : _graph(copy._graph), _contigs(copy._contigs){

		}

		~LoadContigKminmerFunctor(){
			//delete _minimizerParser;
		}



		void operator () (const KminmerList& kminmerList) {


			u_int64_t readIndex = kminmerList._readIndex;
			const vector<u_int64_t>& readMinimizers = kminmerList._readMinimizers;
			//const vector<KmerVec>& kminmers = kminmerList._kminmers;
			const vector<ReadKminmerComplete>& kminmersInfos = kminmerList._kminmersInfo;

			vector<u_int32_t> nodepath;
			//vector<vector<u_int64_t>> kmerVecs;

			for(size_t i=0; i<kminmersInfos.size(); i++){
				
				const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];

				KmerVec vec = kminmerInfo._vec;
				
				if(_graph._mdbgFirst->_dbg_nodes.find(vec) == _graph._mdbgFirst->_dbg_nodes.end()){
					_graph._logFile << "Not found kminmer" << endl;
					//getchar();
					continue;
				}



				u_int32_t nodeName = _graph._mdbgFirst->_dbg_nodes[vec]._index;
				nodepath.push_back(nodeName);
				//kmerVecs.push_back(vec._kmers);

				//#pragma omp critical
				//{
					
				//	vector<u_int32_t>& contigIndexes = _graph._nodeName_to_contigs[nodeName]; 
				//	if(std::find(contigIndexes.begin(), contigIndexes.end(), readIndex) == contigIndexes.end()){
				//		contigIndexes.push_back(readIndex);
				//	}
				//}
				

			}

			//_nbContigs += 1;

			#pragma omp critical
			{
				//vector<u_int32_t> nodepath_sorted = nodepath;
				//std::sort(nodepath_sorted.begin(), nodepath_sorted.end());
				_contigs.push_back({readIndex, nodepath, readMinimizers});
				//cout << "load contig: " << nodepath.size() << endl;
			}
			
		}
	};


	vector<phmap::parallel_flat_hash_set<KmerVec>> _bannedKminmers;
	unordered_map<KmerVec, u_int32_t> _kmerVec_to_nodeName;

	void loadBannedKminmers(){
		
		//phmap::parallel_flat_hash_set<KmerVec> lol;
		cout << "loading banned kminmers" << endl;

		_bannedKminmers.resize(_kminmerSize+1);

		string bannedDir = _outputDir + "/bannedKminmers/";
		string ext = ".txt";
		

		for (auto &p : fs::recursive_directory_iterator(bannedDir)){
			if (p.path().extension() == ext){
				string filename =  p.path();
				//cout << filename << endl;
				ifstream bannedKminmers(filename);

				ifstream file(filename);

    			size_t kminmerSize;

				file.read((char*)&kminmerSize, sizeof(kminmerSize));
				//cout << "load kminmer size: " << kminmerSize << endl;

				while(true){

					vector<u_int64_t> minimizerSeq;

					minimizerSeq.resize(kminmerSize);
					file.read((char*)&minimizerSeq[0], kminmerSize*sizeof(u_int64_t));

					if(file.eof()) break;

					KmerVec vec;
					vec._kmers = minimizerSeq;
					//cout << vec._kmers[0] << endl;
					//cout << vec._kmers.size() << endl;
					//lol.insert(vec.normalize());
					_bannedKminmers[kminmerSize].insert(vec.normalize());
				}
			}
		}

		//cout << "done" << endl;

	}
	*/

	//void createMDBG_collectKminmers_minspace_read(const vector<u_int64_t>& readMinimizers, const vector<ReadKminmerComplete>& kminmersInfos, u_int64_t readIndex);

	//void createSimilarityGraph(GraphInfo* graphInfo);
	//void execute_binning();
	//void execute_binning_cleanGraph();
	//void extract_kminmers();

	u_int64_t _debug_nbMinimizers;
	gzFile _file_readData;
	MinimizerParser* _minimizerParser;
	//unordered_map<KmerVec, u_int32_t> _solidKminmerAbundances;
	u_int64_t _nbSolidKminmers;

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
		float sd = sqrt(var);
		
		cout << "Mean: " << mean << endl;
		cout << "Sd: " << sd << endl;
		cout << "Var: " << var << endl;
	}

	double compute_median(vector<u_int32_t> scores){
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


	/*
	class FillBloomFilter {

		public:

		CreateMdbg& _graph;
		BloomCacheCoherent<u_int64_t>* _bloomFilter;
		u_int64_t _nbSolidKminmers;

		FillBloomFilter(CreateMdbg& graph) : _graph(graph){
			_bloomFilter = graph._bloomFilter;
		}

		FillBloomFilter(const FillBloomFilter& copy) : _graph(copy._graph){
			_bloomFilter = copy._bloomFilter;
			_nbSolidKminmers = 0;
		}

		~FillBloomFilter(){
			cout << _nbSolidKminmers << endl;
			//delete _minimizerParser;
		}

		void operator () (const KminmerList& kminmerList) {

			u_int64_t readIndex = kminmerList._readIndex;
			const vector<u_int64_t>& readMinimizers = kminmerList._readMinimizers;
			//const vector<KmerVec>& kminmers = kminmerList._kminmers;
			const vector<ReadKminmerComplete>& kminmersInfos = kminmerList._kminmersInfo;

			for(size_t i=0; i<kminmersInfos.size(); i++){

				const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];
				const KmerVec& vec = kminmerInfo._vec;

				u_int64_t val = 0;

				if(_bloomFilter->contains(val)){
					_nbSolidKminmers += 1;
					__sync_fetch_and_add(&_graph._nbSolidKminmers, 1);
				}
				else{
					_bloomFilter->insert(val);
				}
			}

		}
	};
	*/

	vector<u_int32_t> _contigsNbNodes;
	u_int64_t _nbSmallContigs;
	/*
	class IndexContigFunctor { 

		public:

		CreateMdbg& _graph;

		IndexContigFunctor(CreateMdbg& graph) : _graph(graph){
		}

		IndexContigFunctor(const IndexContigFunctor& copy) : _graph(copy._graph){
		}

		~IndexContigFunctor(){
		}
		
		void operator () (const KminmerList& kminmerList) {
			
			//#pragma omp critical
			//{
			_graph._contigsNbNodes.push_back(kminmerList._kminmersInfo.size());
			
			for(const ReadKminmerComplete& rkc : kminmerList._kminmersInfo){

				if(_graph._kminmerAbundances.find(rkc._vec) != _graph._kminmerAbundances.end()){
					_graph._kminmerAbundances[rkc._vec]._contigIndex = kminmerList._readIndex;
					//cout << _graph._kminmerAbundances[rkc._vec]._contigIndex << " " << _graph._kminmerAbundances[rkc._vec]._abundance << endl;
					//getchar();
				}
			}
			//}
			//_mdbg->_dbg_nodes[vec] = {nodeName, abundance, quality};
		}

	};
	*/
	
    struct ReadWriter{
        u_int64_t _readIndex;
		vector<MinimizerType> _contigSequence;
		bool _isCircular;
    };

    struct ReadWriter_Comparator {
        bool operator()(ReadWriter const& p1, ReadWriter const& p2){
            return p1._readIndex > p2._readIndex;
        }
    };

	priority_queue<ReadWriter, vector<ReadWriter> , ReadWriter_Comparator> _readWriterQueue;
	u_int64_t _nextReadIndexWriter;

	void writeSmallContig(u_int64_t readIndex, const vector<MinimizerType>& contigSequence, bool isCircular){

		#pragma omp critical(writeSmallContig)
		{
			_readWriterQueue.push({readIndex, contigSequence, isCircular});

			while(!_readWriterQueue.empty()){

				const ReadWriter& readWriter = _readWriterQueue.top();

				if(readWriter._readIndex == _nextReadIndexWriter){

					if(readWriter._contigSequence.size() > 0){
						u_int32_t contigSize = readWriter._contigSequence.size();
						_fileSmallContigs.write((const char*)&contigSize, sizeof(contigSize));
						_fileSmallContigs.write((const char*)&readWriter._isCircular, sizeof(readWriter._isCircular));
						_fileSmallContigs.write((const char*)&readWriter._contigSequence[0], contigSize*sizeof(MinimizerType));

						_nbSmallContigs += 1;
					}

					_readWriterQueue.pop();
					_nextReadIndexWriter += 1;
				}
				else{
					break;
				}
			}
			
		}

	}


	class IndexKminmerFunctor {

		public:

		CreateMdbg& _parent;
		bool _extractingContigs;

		IndexKminmerFunctor(CreateMdbg& parent, bool extractingContigs) : _parent(parent){
			_extractingContigs = extractingContigs;
		}

		IndexKminmerFunctor(const IndexKminmerFunctor& copy) : _parent(copy._parent){
			_extractingContigs = copy._extractingContigs;
		}

		~IndexKminmerFunctor(){
			//delete _minimizerParser;
		}

		/*
		float getAbundance(const vector<MinimizerType>& readMinimizers, const ReadKminmerComplete& kminmerInfo){

			const KmerVec& vec = kminmerInfo._vec;

			vector<MinimizerType> minimizerSeq;
			for(size_t i=kminmerInfo._read_pos_start; i<=kminmerInfo._read_pos_end; i++){
				minimizerSeq.push_back(readMinimizers[i]);
			}
			if(kminmerInfo._isReversed){
				std::reverse(minimizerSeq.begin(), minimizerSeq.end());
			}

			return getAbundance(minimizerSeq);
		}
		*/

		//float getAbundance(const vector<MinimizerType>& readMinimizers){
		float getAbundance(const size_t pos, const vector<u_int32_t>& prevAbundances){

			if(prevAbundances.size() <= 1) return prevAbundances[0];

			size_t nbSubKminmers = 2; //_parent._kminmerSize - (_parent._kminmerSizeFirst+1) + 1;
			
			//cout << pos << " " << _parent._kminmerSize << " " << nbSubKminmers << endl;

			u_int32_t minAbundance = -1;


			for(size_t i=pos; i < pos+nbSubKminmers; i++){

				const u_int32_t abundance = prevAbundances[i];

				//cout << "\t" << i << endl;
				if(abundance < minAbundance){
					minAbundance = abundance;
				}

			}

			return minAbundance;

			/*
			vector<u_int64_t> rlePositions;
			vector<u_int32_t> minimizers_pos;//(minimizers.size());
			vector<KmerVec> kminmers; 
			vector<ReadKminmer> kminmersInfo;
			MDBG::getKminmers(0, _parent._kminmerSizePrev, minimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, 0, false);
			//MDBG::getKminmers(0, _parent._kminmerSizeFirst+1, minimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, 0, false);

			//string s = "";

			u_int32_t minAbundance = -1;
			//u_int32_t maxAbundance = 0;
			//float sum = 0;
			//float n = 0;
			//vector<float> abundances;
			//cout << kminmers.size() << " " << _kminmerAbundances.size() << endl;
			for(const KmerVec& vec : kminmers){

				//cout << vec._kmers.size() << endl;
				//cout << vec.toString() << endl;

				const u_int128_t hash = vec.hash128();

				//cout << (_mdbgInit->_dbg_nodes.find(vec) != _mdbgInit->_dbg_nodes.end()) << endl;
				if(_parent._kminmerAbundances.find(hash) != _parent._kminmerAbundances.end()){

					//u_int32_t abundance = _kminmerAbundances[vec];
					//const KminmerAbundanceMapNode& node = _kminmerAbundances[vec];
					const u_int32_t abundance = _parent._kminmerAbundances[hash];
					//if(nodeName == 17430) cout << abundance << endl;

					if(abundance == 0){
						minAbundance = 1;
						break;
					}

					if(abundance < minAbundance){
						minAbundance = abundance;
					}

					//if(abundance > maxAbundance){
					//	maxAbundance = abundance;
					//}

					//s += to_string(abundance) + " ";
					//cout << abundance << " ";
					//abundances.push_back(abundance);

					//sum += abundance;
					//n += 1;
				}
				else{
					minAbundance = 1;
					break;
				}
			}


			return minAbundance;
			*/
		}

		/*
		float getAbundance2(const size_t pos, const vector<u_int32_t>& prevAbundances){

			size_t nbSubKminmers = 2; //_parent._kminmerSize - (_parent._kminmerSizeFirst+1) + 1;

			//cout << pos << " " << _parent._kminmerSize << " " << nbSubKminmers << endl;

			u_int32_t minAbundance = -1;


			for(size_t i=pos; i < pos+nbSubKminmers; i++){

				const u_int32_t abundance = prevAbundances[i];

				//cout << "\t" << i << endl;
				if(abundance < minAbundance){
					minAbundance = abundance;
				}

			}

			return minAbundance;

			
		}
		*/
		/*
		void indexSolidTriplets(const vector<MinimizerType>& readMinimizers){

			vector<u_int64_t> rlePositions;
			vector<u_int32_t> minimizers_pos;//(minimizers.size());
			vector<KmerVec> kminmers; 
			vector<ReadKminmer> kminmersInfo;
			MDBG::getKminmers(_minimizerSize, _graph._kminmerSize+1, readMinimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, 0, false);


			for(const KmerVec& vec : kminmers){

				//cout << vec._kmer

				u_int128_t vecHash = vec.hash128();



				_graph._solidTriplets.lazy_emplace_l(vecHash, 
				[this](MdbgNodeMapLight::value_type& v) { // key exist
				},           
				[&vecHash, this](const MdbgNodeMapLight::constructor& ctor) { // key inserted
					
					
					DbgNodeLight node = {};

					ctor(vecHash, node); 

				}); // construct value_type in place when key not present
						
			}

		}
		*/
		/*
		struct SumN{
			double _sum;
			double _n;
			double _nbNodes_sum;
			double _nbNodes_n;
		};

		SumN getAbundanceMean(const vector<u_int64_t>& readMinimizers, const ReadKminmerComplete& kminmerInfo){

			SumN sumN = {0, 0, 0, 0};

			const KmerVec& vec = kminmerInfo._vec;

			vector<u_int64_t> minimizerSeq;
			for(size_t i=kminmerInfo._read_pos_start; i<=kminmerInfo._read_pos_end; i++){
				minimizerSeq.push_back(readMinimizers[i]);
			}
			if(kminmerInfo._isReversed){
				std::reverse(minimizerSeq.begin(), minimizerSeq.end());
			}


			vector<u_int64_t> rlePositions;
			vector<u_int64_t> minimizers_pos;//(minimizers.size());
			vector<KmerVec> kminmers; 
			vector<ReadKminmer> kminmersInfo;
			MDBG::getKminmers(_minimizerSize, _kminmerSizeFirst, minimizerSeq, minimizers_pos, kminmers, kminmersInfo, rlePositions, 0, false);


			//u_int32_t minAbundance = -1;
			//float sum = 0;
			//float n = 0;
			//vector<float> abundances;
			//cout << kminmers.size() << endl;
			
			float lastAbundance = -1;
			float lastN = -1;

			//cout << kminmers.size() << endl;

			for(const KmerVec& vec : kminmers){
				//cout << (_mdbgInit->_dbg_nodes.find(vec) != _mdbgInit->_dbg_nodes.end()) << endl;


				if(_mdbgInit->_dbg_nodes.find(vec) != _mdbgInit->_dbg_nodes.end()){

					u_int32_t abundance = _mdbgInit->_dbg_nodes[vec]._abundance;
					size_t n = _mdbgInit->_dbg_nodes[vec]._unitigNbNodes;

					//if(abundance == lastAbundance && lastN == n) continue;

					for(size_t i=0; i<n; i++){
						sumN._n += 1;
						sumN._sum += abundance;
					}

					sumN._nbNodes_sum += n;
					sumN._nbNodes_n += 1;

					cout << abundance << " " << n << endl;

					lastAbundance = abundance;
					lastN = n;
					//if(nodeName == 17430) cout << abundance << endl;

					//if(abundance < minAbundance){
					//	minAbundance = abundance;
					//}
					//cout << abundance << " ";
					//abundances.push_back(abundance);

					//sum += abundance;
					//n += 1;
				}
				//else{
					//minAbundance = 1;
					//break;
				//}
			}

			return sumN;
		}
		*/
		/*
		bool isInGraph(const vector<ReadKminmerComplete>& kminmersInfos){

			if(kminmersInfos.size() == 0) return false;

			for(size_t i=0; i<kminmersInfos.size(); i++){
					

				//const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];

				const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];
				const KmerVec& vec = kminmerInfo._vec;

				if(_parent._mdbgNodesLight.find(vec.hash128()) == _parent._mdbgNodesLight.end()){
					return false;
				}
			}

			return true;
		}
		*/
		
		vector<u_int32_t> getPrevAbundances(const vector<MinimizerType>& readMinimizers){

			vector<u_int32_t> prevAbundances;

			vector<u_int64_t> rlePositions;
			vector<u_int32_t> minimizers_pos;//(minimizers.size());
			vector<KmerVec> kminmers; 
			vector<ReadKminmer> kminmersInfo;
			MDBG::getKminmers(0, _parent._kminmerSizePrev, readMinimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, 0, false);
			
			for(const KmerVec& vec : kminmers){

				const u_int128_t vecHash = vec.hash128();

				if(_parent._kminmerAbundances.find(vecHash) != _parent._kminmerAbundances.end()){

					prevAbundances.push_back(_parent._kminmerAbundances[vecHash]);
					
				}
				else{
					prevAbundances.push_back(1);
				}
			}

			return prevAbundances;
		}
		

		void operator () (const KminmerList& kminmerList) {


			u_int64_t readIndex = kminmerList._readIndex;
			const vector<MinimizerType>& readMinimizers = kminmerList._readMinimizers;
			//const vector<KmerVec>& kminmers = kminmerList._kminmers;
			const vector<ReadKminmerComplete>& kminmersInfos = kminmerList._kminmersInfo;

			vector<u_int32_t> prevAbundances = getPrevAbundances(readMinimizers);
			/*
			cout << "---" << endl;
			for(size_t i=0; i<prevAbundances.size(); i++){
				cout << i << "\t" << prevAbundances[i] << endl;
			}
			cout << endl;
			for(size_t i=0; i<kminmerList._abundances.size(); i++){
				cout << i << "\t" << kminmerList._abundances[i] << endl;
			}
			cout << prevAbundances.size() << " " << kminmerList._abundances.size() << endl;
			*/

			//for(size_t i=0; i<prevAbundances.size(); i++){
			//	if(prevAbundances[i] != kminmerList._abundances[i]){
			//		cout << "derp" << endl;
			//		getchar();
			//	}
			//}
			
			//if(_graph._savingSmallContigs){

			/*
			if(_extractingContigs && _kminmerSize > 8 && kminmersInfos.size() <= 0 && getAbundance(readMinimizers) > 1){

				if(!isInGraph(kminmersInfos)){
					_graph.writeSmallContig(readIndex, readMinimizers, false);
				}
				else{
					_graph.writeSmallContig(readIndex, {}, false);
				}
				
				return;
			}

			if(_extractingContigs){
				_graph.writeSmallContig(readIndex, {}, false);
			}
			*/
			//if(_graph._savingSmallContigs){
			/*
			#pragma omp critical(writeSmallContig)
			{
				if(_extractingContigs && _parent._kminmerSize > 8 && kminmersInfos.size() <= 0){
					cout << "---" << endl;
					cout << getAbundance(0, kminmerList._abundances) << endl;
					cout << kminmerList._abundances.size() << endl;
					for(size_t i=0; i<kminmerList._abundances.size(); i++){
						cout << "\t" << kminmerList._abundances[i] << endl;
					}
				}
			}
			*/
		
			if(_extractingContigs && _parent._kminmerSize > 8 && kminmersInfos.size() <= 0 && getAbundance(0, prevAbundances) > 1){

				//if(!isInGraph(kminmersInfos)){
				#pragma omp critical(writeSmallContig)
				{
					u_int32_t contigSize = readMinimizers.size();
					_parent._fileSmallContigs.write((const char*)&contigSize, sizeof(contigSize));

					u_int8_t isCircular = kminmerList._isCircular;
					_parent._fileSmallContigs.write((const char*)&isCircular, sizeof(isCircular));
					_parent._fileSmallContigs.write((const char*)&readMinimizers[0], contigSize*sizeof(MinimizerType));
					//cout << "small contig" << endl;
					//getchar();

					//if(kminmersInfos.size() > 0){
					_parent._nbSmallContigs += 1;
					//}
				}

				//}

				return;
			}




			//cout << "-----" << endl;

			for(size_t i=0; i<kminmersInfos.size(); i++){
				



				const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];
				const KmerVec& vec = kminmerInfo._vec;
				u_int128_t vecHash = vec.hash128();
				//if(vec.isPalindrome()) continue;

				/*
				int n=0;
				for(MinimizerType m : vec._kmers){
					if(m == 68874849684788112){
						n += 1;
					}
				}

				if(n > 3){
					cout << endl;
					for(auto m : readMinimizers){
						cout << m << endl;
					}
				}
				*/
				/*
				int n=0;
				for(MinimizerType m : vec._kmers){
					if(m == 92102290589771269){
						n += 1;
					}
					if(m == 54785415165247199){
						n += 1;
					}
					if(m == 1904695471279270){
						n += 1;
					}
					if(m == 92170037183223943){
						n += 1;
					}
					if(m == 69724447621269186){
						n += 1;
					}
				}


				if(n == 5){
					cout << "huhuhu" << endl;
					getchar();
				}
				*/





				u_int32_t kminmerAbundance = getAbundance(i, prevAbundances);

				/*
				if(_extractingContigs){

					//if(kminmerAbundance != getAbundance2(i, kminmerList._abundances)){
					//	cout << "derp: " << kminmerAbundance << " " << getAbundance2(i, kminmerList._abundances) << endl;
					//	getchar();
					//}

					kminmerAbundance = getAbundance2(i, kminmerList._abundances);
				}
				else{

					//if(kminmerAbundance != getAbundance(i, kminmerList._abundances)){
					//	cout << "derp: " << kminmerAbundance << " " << getAbundance2(i, kminmerList._abundances) << endl;
					//	getchar();
					//}
					
					kminmerAbundance = getAbundance(i, kminmerList._abundances);
				}
				*/

				if(kminmerAbundance <= 1){
					continue;
				}

				//if(_parent._mdbgNodesLight.find(vecHash) != _parent._mdbgNodesLight.end()) continue;

				//#pragma omp critical(IndexKminmerFunctor)
				//{
				//	_parent._mdbgNodesLight[vecHash] = kminmerAbundance;
				//}

				
				_parent._mdbgNodesLight.lazy_emplace_l(vecHash, 
				[](KminmerAbundanceMap::value_type& v) { // key exist
					//v.second += 1;
				},           
				[&vecHash, &kminmerAbundance](const KminmerAbundanceMap::constructor& ctor) { // key inserted
					
					//ctor(vecHash, 1); 
					ctor(vecHash, kminmerAbundance); 

				}); // construct value_type in place when key not present
				
				
			}

		}
	};
	



	/*
	class FilterKminmerFunctor2 {

		public:

		CreateMdbg& _graph;
		bool _isFirstPass;
		bool _parsingContigs;
		ofstream& _readFile;
		unordered_set<KmerVec>& _kminmerExist;
		MDBG* _mdbg;
		KminmerAbundanceMap& _kminmerAbundances;
		ofstream& _kminmerFile;
		float _minimizerSpacingMean;
		float _kminmerLengthMean;
		float _kminmerOverlapMean;
		size_t _minimizerSize;
		size_t _kminmerSize;
		size_t _kminmerSizeFirst;
		//BloomCacheCoherent<u_int64_t>* _bloomFilter;

		FilterKminmerFunctor2(CreateMdbg& graph) : _graph(graph), _readFile(graph._readFile), _kminmerExist(graph._kminmerExist), _kminmerAbundances(graph._kminmerAbundances), _kminmerFile(graph._kminmerFile){
			_isFirstPass = graph._isFirstPass;
			_parsingContigs = graph._parsingContigs;
			_mdbg = graph._mdbg;
			//_mdbgInit = graph._mdbgInit;
			_minimizerSpacingMean = graph._minimizerSpacingMean;
			_kminmerLengthMean = graph._kminmerLengthMean;
			_kminmerOverlapMean = graph._kminmerOverlapMean;
			_minimizerSize = graph._minimizerSize;
			_kminmerSize = graph._kminmerSize;
			_kminmerSizeFirst = graph._kminmerSizeFirst;
			//_bloomFilter = graph._bloomFilter;
		}

		FilterKminmerFunctor2(const FilterKminmerFunctor2& copy) : _graph(copy._graph), _readFile(copy._readFile), _kminmerExist(copy._kminmerExist), _kminmerAbundances(copy._kminmerAbundances), _kminmerFile(copy._kminmerFile){
			_isFirstPass = copy._isFirstPass;
			_parsingContigs = copy._parsingContigs;
			_mdbg = copy._mdbg;
			//_mdbgInit = copy._mdbgInit;
			_minimizerSpacingMean = copy._minimizerSpacingMean;
			_kminmerLengthMean = copy._kminmerLengthMean;
			_kminmerOverlapMean = copy._kminmerOverlapMean;
			_minimizerSize = copy._minimizerSize;
			_kminmerSize = copy._kminmerSize;
			_kminmerSizeFirst = copy._kminmerSizeFirst;
			//_bloomFilter = copy._bloomFilter;
		}

		~FilterKminmerFunctor2(){
			//delete _minimizerParser;
		}

	

		void operator () (const KminmerList& kminmerList) {



			u_int64_t readIndex = kminmerList._readIndex;
			const vector<MinimizerType>& readMinimizers = kminmerList._readMinimizers;
			//const vector<KmerVec>& kminmers = kminmerList._kminmers;
			const vector<ReadKminmerComplete>& kminmersInfos = kminmerList._kminmersInfo;

			
		
			//float readAbundance = getReadAbundance(readMinimizers);

			int nbKminmers = 0;
			float readAbundance = getReadAbundance(readMinimizers, nbKminmers);
			double cutoff = readAbundance * 0.1f;

			


			
			//if(cutoff > 1) return;

			for(size_t i=0; i<kminmersInfos.size(); i++){
					

				//const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];

				const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];
				const KmerVec& vec = kminmerInfo._vec;
				//if(vec.isPalindrome()) continue;

				float kminmerAbundance = getKminmerAbundance(readMinimizers, kminmerInfo);
				
				//cout << kminmerAbundance << " " << readAbundance << endl;
				
				//if(kminmerAbundance == 1) getchar();
				//if(kminmerAbundance == 1 && cutoff > 1) continue;
				//if(kminmerAbundance == 1) continue;
				//if(kminmerAbundance == 1) continue; // && kminmerAbundance < cutoff) continue;
				if(kminmerAbundance == 0) continue;
				if(kminmerAbundance < cutoff) continue;
				//if(kminmerAbundance <= 1 && kminmerAbundance < cutoff) continue;


				bool isNewKey = false;

				_mdbg->_dbg_nodes.lazy_emplace_l(vec, 
				[this, &kminmerInfo](MdbgNodeMap::value_type& v) { // key exist
					//v.second._quality += kminmerInfo._quality;
					//if(_isFirstPass){
					//	v.second._abundance += 1;
					//}
				},           
				[&vec, this, &kminmerInfo, &readMinimizers, &readIndex, &isNewKey, &kminmerAbundance](const MdbgNodeMap::constructor& ctor) { // key inserted
					
					
					isNewKey = true;


					u_int32_t nodeName = -1;
					DbgNode node = {nodeName, 2};

					node._abundance = kminmerAbundance; 

					//node._quality = 0;
					//node._quality += kminmerInfo._quality;

					ctor(vec, node); 

					//cout << "jinsere" << endl;
				}); // construct value_type in place when key not present

				
				if(isNewKey){

					u_int32_t nodeName;

					#pragma omp critical
					{

						
						nodeName = _graph._node_id;
						_graph._node_id += 1;
						

					}

					auto set_value = [&nodeName](MdbgNodeMap::value_type& v) { v.second._index = nodeName; };
					_mdbg->_dbg_nodes.modify_if(vec, set_value);
				}


			}


		}
		
		float getReadAbundance(const vector<MinimizerType>& readMinimizers, int& nbKminmers){

			vector<u_int64_t> rlePositions;
			vector<u_int32_t> minimizers_pos;//(minimizers.size());
			vector<KmerVec> kminmers; 
			vector<ReadKminmer> kminmersInfo;
			MDBG::getKminmers(_minimizerSize, _graph._kminmerSizePrev, readMinimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, 0, false);



			//u_int32_t minAbundance = -1;
			//float sum = 0;
			//float n = 0;
			vector<float> abundances;
			unordered_set<u_int32_t> processedContigIndex;
			//cout << "----" << endl;
			
			for(const KmerVec& vec : kminmers){

				//cout << (_mdbgInit->_dbg_nodes.find(vec) != _mdbgInit->_dbg_nodes.end()) << endl;
				if(_kminmerAbundances.find(vec) != _kminmerAbundances.end()){

					KminmerAbundanceMapNode node = _kminmerAbundances[vec];
					if(processedContigIndex.find(node._contigIndex) != processedContigIndex.end()) continue;

					processedContigIndex.insert(node._contigIndex);

					//cout << "ctg: " << node._contigIndex << " " << _graph._contigsNbNodes[node._contigIndex] << " " << node._abundance << endl;
					
					for(size_t i=0; i<_graph._contigsNbNodes[node._contigIndex]; i++){
						abundances.push_back(node._abundance);
					}
				}
				else{
					abundances.push_back(1);
				}
			}

			nbKminmers = abundances.size();
			float abundance = Utils::compute_median_float(abundances);

			return abundance;
		}
		

		float getKminmerAbundance(const vector<MinimizerType>& readMinimizers, const ReadKminmerComplete& kminmerInfo){

			const KmerVec& vec = kminmerInfo._vec;

			vector<MinimizerType> minimizerSeq;
			for(size_t i=kminmerInfo._read_pos_start; i<=kminmerInfo._read_pos_end; i++){
				minimizerSeq.push_back(readMinimizers[i]);
			}
			if(kminmerInfo._isReversed){
				std::reverse(minimizerSeq.begin(), minimizerSeq.end());
			}

			return getAbundance(minimizerSeq);
		}

		float getAbundance(const vector<MinimizerType>& readMinimizers){

			vector<u_int64_t> rlePositions;
			vector<u_int32_t> minimizers_pos;//(minimizers.size());
			vector<KmerVec> kminmers; 
			vector<ReadKminmer> kminmersInfo;
			MDBG::getKminmers(_minimizerSize, _graph._kminmerSizePrev, readMinimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, 0, false);


			//cout << "-----" << endl;
			//cout << kminmers.size() << endl;

			u_int32_t minAbundance = -1;
			float sum = 0;
			float n = 0;
			//vector<float> abundances;
			//cout << kminmers.size() << endl;
			for(const KmerVec& vec : kminmers){
				//cout << (_mdbgInit->_dbg_nodes.find(vec) != _mdbgInit->_dbg_nodes.end()) << endl;
				if(_kminmerAbundances.find(vec) != _kminmerAbundances.end()){

					//u_int32_t abundance = _kminmerAbundances[vec];
					//const KminmerAbundanceMapNode& node = _kminmerAbundances[vec];
					u_int32_t abundance = _kminmerAbundances[vec]._abundance;

					//cout << abundance << endl;
					//if(nodeName == 17430) cout << abundance << endl;

					if(abundance < minAbundance){
						minAbundance = abundance;
					}
					//cout << abundance << " ";
					//abundances.push_back(abundance);

					sum += abundance;
					n += 1;
				}
				else{
					minAbundance = 0;
					break;
				}
			}

			return minAbundance;
		}

	};
	*/


	/*
	class DumpKminmerFunctor {

		public:

		CreateMdbg& _graph;
		ofstream& _kminmerFile;

		DumpKminmerFunctor(CreateMdbg& graph, ofstream& kminmerFile) : _graph(graph), _kminmerFile(kminmerFile){
		}

		DumpKminmerFunctor(const DumpKminmerFunctor& copy) : _graph(copy._graph), _kminmerFile(copy._kminmerFile){
		}

		~DumpKminmerFunctor(){
		}


		void operator () (const KminmerList& kminmerList) {


			u_int64_t readIndex = kminmerList._readIndex;
			const vector<MinimizerType>& readMinimizers = kminmerList._readMinimizers;
			//const vector<KmerVec>& kminmers = kminmerList._kminmers;
			const vector<ReadKminmerComplete>& kminmersInfos = kminmerList._kminmersInfo;

			#pragma omp critical(DumpKminmerFunctor)
			{

				for(size_t i=0; i<kminmersInfos.size(); i++){
					
					const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];
					const KmerVec& vec = kminmerInfo._vec;

					if(_graph._mdbgNodesLight.find(vec.hash128()) == _graph._mdbgNodesLight.end()) continue;

					DbgNodeLight& node = _graph._mdbgNodesLight[vec.hash128()];
					if(node._isWritten) continue;

					node._isWritten = true;
					MDBG::writeKminmer(node._index, node._abundance, vec._kmers, _kminmerFile);

					
				}

			}
		}
	};
	*/

	
	class ComputeUnitigFunctor {

		public:

		CreateMdbg& _graph;

		ComputeUnitigFunctor(CreateMdbg& graph) : _graph(graph){
		}

		ComputeUnitigFunctor(const ComputeUnitigFunctor& copy) : _graph(copy._graph){
		}

		~ComputeUnitigFunctor(){
		}

		/*
		void operator () (const KminmerList& kminmerList) {


			u_int64_t readIndex = kminmerList._readIndex;
			const vector<MinimizerType>& readMinimizers = kminmerList._readMinimizers;
			//const vector<KmerVec>& kminmers = kminmerList._kminmers;
			const vector<ReadKminmerComplete>& kminmersInfos = kminmerList._kminmersInfo;

			//#pragma omp critical(DumpKminmerFunctor)
			//{

			for(size_t i=0; i<kminmersInfos.size(); i++){
				
				const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];
				const KmerVec& vec = kminmerInfo._vec;

				computeUnitigNode(vec);
				
			}

			//}
		}
		*/

		void setNodeUnitigged(const KmerVec& vec){

			//#pragma omp critical
			//{
			//KmerVec lol;
			//lol._kmers = {77646226911515249, 90273684072333150, 83465560366503996};

			bool isReversedPrefix;
			bool isReversedSuffix;
			KmerVec prefix = vec.prefix().normalize(isReversedPrefix);
			KmerVec suffix = vec.suffix().normalize(isReversedSuffix);
			MinimizerType prefix_minimizer = vec._kmers[vec._kmers.size()-1];
			MinimizerType suffix_minimizer = vec._kmers[0];
			
			/*
			u_int64_t m1 = -1;
			u_int64_t m2 = -1;


			for(KminmerEdge3& edge : _graph._mdbgEdges3[suffix.hash128()]){
				if(edge._isReversed != isReversedSuffix) continue;
				if(edge._isPrefix) continue;

				if(edge._minimizer == suffix_minimizer){
					//if(suffix == lol){
					//	cout << "omg suffix " << suffix.toString() << " " << isReversedSuffix << " " << suffix_minimizer << endl;
					//	getchar();
					//}
					edge._isUnitigged = true;
					m1 = edge._minimizer;

					break;
				}
			}
			
			

			
			for(KminmerEdge3& edge : _graph._mdbgEdges3[prefix.hash128()]){
				if(edge._isReversed != isReversedPrefix) continue;
				if(!edge._isPrefix) continue;

				if(edge._minimizer == prefix_minimizer){
					//if(prefix == lol){
						//cout << "omg prefix " << prefix.toString() << " " << isReversedPrefix << " " << prefix_minimizer << endl;
						//getchar();
					//}
					edge._isUnitigged = true;
					m2 = edge._minimizer;

					break;
				}
			}
			*/

			u_int64_t m3 = -1;
			u_int64_t m4 = -1;


			//KminmerEdge33& edge1 = _graph._mdbgEdges5[suffix.hash128()];
			u_int128_t key1 = _graph._mdbgEdges10._keys->lookup(suffix.hash128());
			KminmerEdge33& edge1 = _graph._mdbgEdges10._values[key1];

			if(edge1._minimizer1 != -1){
				if(edge1._isReversed1 == isReversedSuffix && !edge1._isPrefix1 && edge1._minimizer1 == suffix_minimizer){
					m3 = edge1._minimizer1;
					edge1._isUnitigged1 = true;
				}

				//if(suffix.isPalindrome()){
				//	edge1._isUnitigged1 = true;
				//}
			}

			if(edge1._minimizer2 != -1){
				if(edge1._isReversed2 == isReversedSuffix && !edge1._isPrefix2 && edge1._minimizer2 == suffix_minimizer){
					if(m3 == -1) m3 = edge1._minimizer2;
					//m4 = edge1._minimizer2;
					edge1._isUnitigged2 = true;
				}
				
			}

			u_int128_t key2 = _graph._mdbgEdges10._keys->lookup(prefix.hash128());
			KminmerEdge33& edge2 = _graph._mdbgEdges10._values[key2];
			//KminmerEdge33& edge2 = _graph._mdbgEdges5[prefix.hash128()];

			if(edge2._minimizer1 != -1){
				if(edge2._isReversed1 == isReversedPrefix && edge2._isPrefix1 && edge2._minimizer1 == prefix_minimizer){
					m4 = edge2._minimizer1;
					edge2._isUnitigged1 = true;
				}
				
				//if(prefix.isPalindrome()){
				//	edge2._isUnitigged1 = true;
				//}
			}

			if(edge2._minimizer2 != -1){
				if(edge2._isReversed2 == isReversedPrefix && edge2._isPrefix2 && edge2._minimizer2 == prefix_minimizer){
					if(m4 == -1) m4 = edge2._minimizer2;
					edge2._isUnitigged2 = true;
				}
			}
			
			/*
			cout << "Available: " << edge1._minimizer1 << " " << edge1._minimizer2 << " " << edge2._minimizer1 << " " << edge2._minimizer2 << endl;

			for(KminmerEdge3& edge : _graph._mdbgEdges3[suffix.hash128()]){
				cout << "\tAvailable a: " << edge._minimizer << endl;
			}
			for(KminmerEdge3& edge : _graph._mdbgEdges3[prefix.hash128()]){
				cout << "\tAvailable b: " << edge._minimizer << endl;
			}

			cout << m1 << " " << m2 << "    " << m3 << " " << m4 << endl;

			if(m1 != -1 && m3 == -1){
				cout << "derp a: " << vec.toString() << endl;
				getchar();
			}

			if(m2 != -1 && m4 == -1){
				cout << "derp b: " << vec.toString() << endl;
				getchar();
			}

			//}
			
			*/
		}

		bool isNodeUnitigged(const KmerVec& vec){

			bool isReversedPrefix;
			bool isReversedSuffix;
			KmerVec prefix = vec.prefix().normalize(isReversedPrefix);
			KmerVec suffix = vec.suffix().normalize(isReversedSuffix);
			MinimizerType prefix_minimizer = vec._kmers[vec._kmers.size()-1];
			MinimizerType suffix_minimizer = vec._kmers[0];
			
			
			u_int128_t key1 = _graph._mdbgEdges10._keys->lookup(suffix.hash128());
			KminmerEdge33& edge1 = _graph._mdbgEdges10._values[key1];
			//KminmerEdge33& edge1 = _graph._mdbgEdges5[suffix.hash128()];

			if(edge1._minimizer1 != -1){
				if(edge1._isReversed1 == isReversedSuffix && !edge1._isPrefix1 && edge1._minimizer1 == suffix_minimizer){
					return edge1._isUnitigged1;
				}

				//if(suffix.isPalindrome()){
				//	return edge1._isUnitigged1;
				//}
			}

			if(edge1._minimizer2 != -1){
				if(edge1._isReversed2 == isReversedSuffix && !edge1._isPrefix2 && edge1._minimizer2 == suffix_minimizer){
					return edge1._isUnitigged2;
				}
				
			}





			u_int128_t key2 = _graph._mdbgEdges10._keys->lookup(prefix.hash128());
			KminmerEdge33& edge2 = _graph._mdbgEdges10._values[key2];
			//KminmerEdge33& edge2 = _graph._mdbgEdges5[prefix.hash128()];

			if(edge2._minimizer1 != -1){
				if(edge2._isReversed1 == isReversedPrefix && edge2._isPrefix1 && edge2._minimizer1 == prefix_minimizer){
					return edge2._isUnitigged1;
				}
				
				//if(prefix.isPalindrome()){
				//	return edge2._isUnitigged1;
				//}
			}

			if(edge2._minimizer2 != -1){
				if(edge2._isReversed2 == isReversedPrefix && edge2._isPrefix2 && edge2._minimizer2 == prefix_minimizer){
					return edge2._isUnitigged2;
				}
			}
			

			
			/*
			for(KminmerEdge3& edge : _graph._mdbgEdges3[suffix.hash128()]){
				if(edge._isReversed != isReversedSuffix) continue;
				if(edge._isPrefix) continue;

				if(edge._minimizer == suffix_minimizer){
					//cout << "Suffix: " << suffix.toString() << " " << isReversedSuffix << " " << suffix_minimizer << endl;
					return edge._isUnitigged;
				}
			}

			for(KminmerEdge3& edge : _graph._mdbgEdges3[prefix.hash128()]){
				if(edge._isReversed != isReversedPrefix) continue;
				if(!edge._isPrefix) continue;

				if(edge._minimizer == prefix_minimizer){
					//cout << "Prefix: " << prefix.toString() << " " << isReversedPrefix << " " << prefix_minimizer << endl;
					return edge._isUnitigged;
				}
			}
			*/
			return false;
		}

		
		bool isNodeExists(const KmerVec& vec){

			
			//cout << "Node exists?: " << vec.toString() << " " << _graph._kminmerAbundances[vec.normalize().hash128()] << endl;
			bool isReversedPrefix;
			bool isReversedSuffix;
			KmerVec prefix = vec.prefix().normalize(isReversedPrefix);
			KmerVec suffix = vec.suffix().normalize(isReversedSuffix);
			MinimizerType prefix_minimizer = vec._kmers[vec._kmers.size()-1];
			MinimizerType suffix_minimizer = vec._kmers[0];
		
			/*
 			if(_graph._mdbgEdges5.find(suffix.hash128()) == _graph._mdbgEdges5.end()){
				//cout << "derp 1" << endl;
				return false;
			}
			if(_graph._mdbgEdges5.find(prefix.hash128()) == _graph._mdbgEdges5.end()){
				//cout << "derp 2" << endl;
				return false;
			}

			KminmerEdge33& edge1 = _graph._mdbgEdges5[suffix.hash128()];

			if(edge1._minimizer1 != -1){
				if(edge1._isReversed1 == isReversedSuffix && !edge1._isPrefix1 && edge1._minimizer1 == suffix_minimizer){
					return true;
				}
			}

			if(edge1._minimizer2 != -1){
				if(edge1._isReversed2 == isReversedSuffix && !edge1._isPrefix2 && edge1._minimizer2 == suffix_minimizer){
					return true;
				}
			}




			KminmerEdge33& edge2 = _graph._mdbgEdges5[prefix.hash128()];

			if(edge2._minimizer1 != -1){
				if(edge2._isReversed1 == isReversedSuffix && edge2._isPrefix1 && edge2._minimizer1 == prefix_minimizer){
					return true;
				}
			}

			if(edge2._minimizer2 != -1){
				if(edge2._isReversed2 == isReversedPrefix && edge2._isPrefix2 && edge2._minimizer2 == prefix_minimizer){
					return true;
				}
			}
			*/
			
 			if(_graph._mdbgEdges3.find(suffix.hash128()) == _graph._mdbgEdges3.end()){
				//cout << "derp 1" << endl;
				return false;
			}
			if(_graph._mdbgEdges3.find(prefix.hash128()) == _graph._mdbgEdges3.end()){
				//cout << "derp 2" << endl;
				return false;
			}

			//cout << "\tSuffix: " << suffix.toString() << " " << suffix_minimizer << endl;

			for(KminmerEdge3& edge : _graph._mdbgEdges3[suffix.hash128()]){
				if(edge._isReversed != isReversedSuffix) continue;
				if(edge._isPrefix) continue;


				if(edge._minimizer == suffix_minimizer){
					//cout << "derp 3" << endl;
					return true;
				}
			}


			//cout << "\tPrefix: " << prefix.toString() << " " << prefix_minimizer << endl;

			for(KminmerEdge3& edge : _graph._mdbgEdges3[prefix.hash128()]){

				if(edge._isReversed != isReversedPrefix) continue;
				if(!edge._isPrefix) continue;
				
				//cout << edge._minimizer << " " << edge._isReversed << " " << isReversedSuffix << endl;
				//if(edge._isReversed != isReversedSuffix) continue;
				//if(!edge._isReversed && isReversedSuffix) continue;

				if(edge._minimizer == prefix_minimizer){
					//cout << "derp 4" << endl;
					return true;
				}
			}

			//cout << "derp 5" << endl;
			
			

			return false;
		}
		


		/*
		bool isBranchingNode(const KmerVec& vec){

			vector<KmerVec> successors;
			vector<KmerVec> predecessors;


			//Check successors
			_graph.getSuccessors(vec, successors);
			if(successors.size() != 1) true;

			KmerVec successor = successors[0];
			
			_graph.getPredecessors(successor, predecessors);
			if(predecessors.size() != 1) true;


			//Check predecessors
			_graph.getPredecessors(vec, predecessors);
			if(predecessors.size() != 1) return true;

			KmerVec predecessor = predecessors[0];

			_graph.getSuccessors(predecessor, successors);
			
			if(successors.size() != 1) return true;


			return false;
		}
		*/

		/*
		void computeUnitigNode(const KmerVec& sourceNode){


			const KmerVec& sourceNodeRC = sourceNode.reverse(); //.hash128();
			//bool lala1 = isNodeUnitigged(sourceNode);
			//bool lala2 = isNodeUnitigged(sourceNodeRC);
			//cout << "a" << endl;
			//if(!isBranchingNode(sourceNode)) return;

			//u_int128_t sourceNodeHash = sourceNode.hash128();

			//cout << "\t" << sourceNode.toString() << " " <<  lala1 << " " << lala2 << endl;

			if(!isNodeExists(sourceNode)) return;
			//cout << "b" << endl;
			//if(!isNodeExists(sourceNodeRC)) return;
			//cout << "c" << endl;
			

			//cout << "--------------" << endl;
			//cout << sourceNode.toString() << endl;
			bool exist = false;
			#pragma omp critical(unitig)
			{

				if(isNodeUnitigged(sourceNode)){
					//exist = true;
				}
				else if(isNodeUnitigged(sourceNodeRC)){
					//exist = true;
				}

				//if(_isKminmerIndexed.find(sourceNodeHash) != _isKminmerIndexed.end()){
				//	exist = true;
				//}
				//else if(_isKminmerIndexed.find(sourceNodeRCHash) != _isKminmerIndexed.end()){
				//	exist = true;
				//}

				//if(_isKminmerIndexed.find(sourceNodeHash) != _isKminmerIndexed.end()) exist = true;
				//if(_nodeToUnitig.find(nodeIndex) != _nodeToUnitig.end()) exist = true;
			}



			if(exist) return;

			//bool dummy = false;

			KmerVec startNode = sourceNode;
			KmerVec endNode = sourceNode;

			u_int8_t isCircular = CONTIG_LINEAR;
			//u_int32_t startNode = nodeIndex;
			//u_int32_t endNode = nodeIndex;

			//u_int8_t isCircular = CONTIG_LINEAR;



			KmerVec smallestNode;
			u_int128_t smallestNodeHash = -1;

			//forward
			while(true){

				//cout << "a: " << endNode.toString() << endl;
				//_logFile << endNode << " " << nodeIndex << endl;

				if(endNode.hash128() < smallestNodeHash){
					smallestNode = endNode;
					smallestNodeHash = endNode.hash128();
				}

				vector<KmerVec> successors;
				_graph.getSuccessors(endNode, successors);

				//cout << "a: Successors: " << successors.size() << endl;

				//cout << "---" << endl;
				//cout << successors.size() << endl;
				//cout << endNode.toString() << endl;

				//if(successors.size() > 0){
				//	cout << "\t" << successors[0].toString() << endl;
				//}

				//_logFile << "lala1: " << neighbors.size() << endl;
				if(successors.size() != 1) break;

				//_logFile << endNode << " " << neighbors[0] << endl;

				KmerVec successor = successors[0];

				
				//cout << "\t" << successor.toString() << endl;
				
				vector<KmerVec> predecessors;
				_graph.getPredecessors(successor, predecessors);
				
				//cout << "a: Successors: " << predecessors.size() << endl;
				//for(KmerVec vec : predecessors){
				//	cout << "\t\tPredecessors: " << vec.toString() << endl;
				//}

				//cout << predecessors.size() << endl;
				//if(predecessors.size() > 0){
				//	cout << "\t" << predecessors[0].toString() << endl;
				//}

				if(predecessors.size() != 1) break;

				if(successor == sourceNode){
					endNode = successor;
					isCircular = CONTIG_CIRCULAR;
					//_logFile << "circular" << endl;
					break; //Circular
				}


				endNode = successor;

			}

			//getchar();

			if(!isCircular){

				//backward
				while(true){
					
					//cout << "b: " << startNode.toString() << endl;

					vector<KmerVec> predecessors;
					_graph.getPredecessors(startNode, predecessors);
					

					//_logFile << "loulou1: " << neighbors.size() << endl;
					if(predecessors.size() != 1) break;


					KmerVec predecessor = predecessors[0];

					
					vector<KmerVec> successors;
					_graph.getSuccessors(predecessor, successors);



					//_logFile << "loulou2: " << neighbors.size() << endl;
					if(successors.size() != 1) break;

					startNode = predecessor;


				}

			}
			else{
				startNode = smallestNode;
				endNode = smallestNode;

			}
			//u_int32_t length = 0;
			KmerVec node = startNode;


			size_t i=0;
			u_int32_t nbNodes = 0;
			vector<MinimizerType> unitig;// = node._kmers;
			//u_int32_t lastNode = -1;

			//vector<u_int32_t> nodes;
			//vector<u_int32_t> nodesRC;
			//u_int32_t qualitySum = 0;
			//u_int32_t minQuality = -1;
			//u_int32_t maxQuality = 0;

			//vector<float> abundances;

			while(true){

				//u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(node);


				//if(i == 0){
				//    length += _graph->_kminmerLengthMean;
				//}
				//else{
				//    length += _graph->_kminmerLengthMean - _graph->_kminmerOverlapMean;
					
				//}

				//abundances.push_back(_graph->_nodeName_to_abundance[nodeName]);

				vector<KmerVec> successors;
				_graph.getSuccessors(node, successors);

				if(i == 0){
					unitig = node._kmers;
				}
				else{
					unitig.push_back(node._kmers[node._kmers.size()-1]);
				}
				
				//nodes.push_back(node);
				//nodesRC.push_back(UnitigGraph::nodeIndex_toReverseDirection(node));
				//nbNodes += 1;

				if(isCircular){
					if(i == 0){ 
						if(successors[0] == endNode) break; //Unitig with single node circular, or single node without edges
					}
					else{
						if(node == endNode) break;
					}
				}
				else{
					if(node == endNode) break;
				}

				//lastNode = node;
				node = successors[0];

				i += 1;


			}

			KmerVec startNodeRC = endNode.reverse();
			KmerVec endNodeRC = startNode.reverse();
			//if(isCircular){
			//    abundances.pop_back();
			//}


			//std::reverse(nodesRC.begin(), nodesRC.end());
			//float median = Utils::compute_median_float(abundances);// _graph->compute_median_float(abundances);

			
			#pragma omp critical(unitig)
			{


				bool isValid = true; 
				
				
				if(isNodeUnitigged(startNode)){
					//isValid = false;;
				}
				else if(isNodeUnitigged(endNode)){
					//isValid = false;;
				}
				else if(isNodeUnitigged(startNodeRC)){
					//isValid = false;;
				}
				else if(isNodeUnitigged(endNodeRC)){
					//isValid = false;;
				}


				if(isValid){

					
					if(isCircular && unitig.size() > _graph._kminmerSize){
						unitig.pop_back();
						//nodes.pop_back();
						//nodesRC.erase(nodesRC.begin());
					}

					//for(MinimizerType minimizer : unitig){
					//	cout << "\t" << minimizer << endl;
					//}

					for(size_t i=0; i<unitig.size()-_graph._kminmerSize+1; i++){
						KmerVec vec;
						for(size_t j=0; j<_graph._kminmerSize; j++){
							vec._kmers.push_back(unitig[i+j]);
						}
						
						setNodeUnitigged(vec);
						setNodeUnitigged(vec.reverse());
						//_isKminmerIndexed.insert(vec.hash128());
						//cout << "\t" << i << " " << vec.toString() << endl;
					}
					
					//_isKminmerIndexed.insert(startNode.hash128());
					//_isKminmerIndexed.insert(endNode.hash128());
					//_isKminmerIndexed.insert(startNodeRC.hash128());
					//_isKminmerIndexed.insert(endNodeRC.hash128());
					//u_int32_t unitigIndexRC = _graph->_nextUnitigIndex + 1;

					//for(u_int32_t nodeIndex : nodes){

					//	_graph->_isNodeIndexIndexed[nodeIndex] = true;
					//	_graph->_isNodeIndexIndexed[UnitigGraph::nodeIndex_toReverseDirection(nodeIndex)] = true;


					//}




					_graph.dumpUnitigNode(_graph._unitigIndex, nodeName, unitig);


					//u_int32_t startNode = nodes[0];
					//u_int32_t endNode = nodes[nodes.size()-1];
					//u_int32_t startNodeRC = nodesRC[0];
					//u_int32_t endNodeRC = nodesRC[nodesRC.size()-1];


					//_graph._startNode_to_unitigIndex[startNode] = _graph._unitigIndex;
					//_graph._endNode_to_unitigIndex[endNode] = _graph._unitigIndex;
					//_graph._startNode_to_unitigIndex[startNodeRC] = _graph._unitigIndex+1;
					//_graph._endNode_to_unitigIndex[endNodeRC] = _graph._unitigIndex+1;

					//cout << "utg" << _graph._unitigIndex << " " << unitig.size() << endl;
					//cout << _graph._unitigIndex << " " << unitig.size() << endl;
					//<< " " << isBranchingNode(startNode) << " " << isBranchingNode(endNode) << " " << isBranchingNode(startNodeRC) << " " << isBranchingNode(endNodeRC)  << endl;


					//getchar();
					_graph._unitigIndex += 2;
					//getchar();
					//cout << "\t" << (int) isCircular << endl;
					//cout << "\t" << startNode.toString() << endl;
					//cout << "\t" << endNode.toString() << endl;
					//cout << "\t" << startNodeRC.toString() << endl;
					//cout << "\t" << endNodeRC.toString() << endl;

					//getchar();
				}

			}
			
			
		}
		*/

		void computeUnitigNode2(const KmerVec& sourceNode){
			
			//cout << endl << "Source node: " << sourceNode.toString() << endl;
			//vector<MinimizerType> ms = {95488693, 303676104, 331395039, 303676104};
			//bool print = (sourceNode._kmers == ms);

			const KmerVec& sourceNodeRC = sourceNode.reverse();

			//if(!isNodeExists(sourceNode)){
			//	cout << "No exists: " << sourceNode.toString() << endl;
				//_graph._debugFile << "\tNo exists" << endl;
			//	return;
			//}

		
		
			bool exist = false;
			#pragma omp critical(unitig)
			{

				if(isNodeUnitigged(sourceNode)){
					exist = true;
				}
				else if(isNodeUnitigged(sourceNodeRC)){
					exist = true;
				}

			
			}



			if(exist){
				//_graph._debugFile << "\tIs unitiged" << endl;
				return;
			}

			//if(sourceNode.prefix().isPalindrome()) return;
			//if(sourceNode.suffix().isPalindrome()) return;
			
			

			vector<MinimizerType> unitig;
			KmerVec startNode = sourceNode;
			KmerVec endNode = sourceNode;

			u_int8_t isCircular = CONTIG_LINEAR;
			
			unitig = sourceNode._kmers;

			//cout << "a" << endl;
			//forward
			while(true){

				//cout << endNode.toString() << endl;
				
				KmerVec successor;
				bool hasSingleSuccessor = _graph.getNbSuccessors(endNode, successor);

				if(!hasSingleSuccessor){
					if(_graph._mdbgEdges3.size() > 0 && _graph.hasSingleSuccessor(endNode)){
						cout << "derp 11" << endl;
						getchar();
					}
					//cout << "\t" << endNode.toString() << endl;
					//cout << "\tcc" << endl;
					break;
				}

				//cout << endNode.toString() << " " << hasSingleSuccessor << " " << successor.toString() << endl;
				
				if(_graph._mdbgEdges3.size() > 0 && !_graph.hasSingleSuccessor(endNode)){
					cout << "derp 22" << endl;
					getchar();
				}

				/*
				vector<KmerVec> successors;
				_graph.getSuccessors(endNode, successors);

				const KmerVec& successor = successors[0];
				*/

				KmerVec predecessor;
				bool hasSinglePredecessor = _graph.getNbPredecessors(successor, predecessor);

				if(!hasSinglePredecessor){
					if(_graph._mdbgEdges3.size() > 0 && _graph.hasSinglePredecessor(successor)){
						cout << "derp 33" << endl;
						getchar();
					}
					
					break;
				}

				
				if(_graph._mdbgEdges3.size() > 0 && !_graph.hasSinglePredecessor(successor)){
					cout << "derp 44" << endl;
					getchar();
				}
				

				if(successor == sourceNode){ //Circular
					//cout << "\tcirc" << endl;
					//for(size_t i=0; i<_graph._kminmerSize-1; i++){
					//	unitig.push_back(sourceNode._kmers[i]);
					//}

					endNode = successor;
					startNode = successor;
					isCircular = CONTIG_CIRCULAR;
					break; 
				}

				endNode = successor;
				/*
				u_int128_t hash = endNode.hash128();

				if(hash < smallestVecHash){
					smallestVecHash = hash;
					smallestNode = endNode;
				}
				*/
				unitig.push_back(endNode._kmers[endNode._kmers.size()-1]);

				//cout << endNode.toString() << endl;
			}

			//if(print){
			//cout << "b" << endl;
			//getchar();
			//}
			//cout << unitig.size() << endl;
			//getchar();



			if(!isCircular){

				//cout << "c" << endl;
				vector<MinimizerType> unitigBackward;

				startNode = sourceNode;

				//backward
				while(true){
					

					KmerVec predecessor;
					bool hasSinglePredecessor = _graph.getNbPredecessors(startNode, predecessor);

					if(!hasSinglePredecessor){

						if(_graph._mdbgEdges3.size() > 0 && _graph.hasSinglePredecessor(startNode)){
							cout << "derp 55" << endl;
							getchar();
						}

						break;
					}
					
					if(_graph._mdbgEdges3.size() > 0 && !_graph.hasSinglePredecessor(startNode)){
						cout << "derp 66" << endl;
						getchar();
					}
					
					//vector<KmerVec> predecessors;
					//_graph.getPredecessors(startNode, predecessors);
					
					//KmerVec predecessor = predecessors[0];

					KmerVec successor;
					bool hasSingleSuccessor = _graph.getNbSuccessors(predecessor, successor);

					if(!hasSingleSuccessor){
						if(_graph._mdbgEdges3.size() > 0 &&  _graph.hasSingleSuccessor(predecessor)){
							cout << "derp 77" << endl;
							getchar();
						}
						break;
					}

					
					if(_graph._mdbgEdges3.size() > 0 && !_graph.hasSingleSuccessor(predecessor)){
							cout << "derp 88" << endl;
							getchar();
					}
					
					
					//vector<KmerVec> successors2;
					//_graph.getSuccessors(predecessor, successors2);



					startNode = predecessor;

					unitigBackward.push_back(startNode._kmers[0]);

				}

				std::reverse(unitigBackward.begin(), unitigBackward.end());

				unitig.insert(unitig.begin(), unitigBackward.begin(), unitigBackward.end());

				//if(print) getchar();
			}
			else{

				//cout << "d" << endl;

				//#pragma omp critical(unitig)
				//{

				
				//for(size_t i=0; i<unitig.size(); i++){
				//	cout << unitig[i] << endl;
				//}
				
			

				u_int128_t smallestVecHash = -1; 
				size_t smallestI = -1;
				bool smallestIsReversed = false;

				vector<KmerVec> kminmers = MDBG::minimizersToKminmers(unitig, _graph._kminmerSize);
				
				for(size_t i=0; i<kminmers.size(); i++){

					const KmerVec& kminmer = kminmers[i];

					bool isReversed;
					const KmerVec& kminmerNorm = kminmer.normalize(isReversed);

					const u_int128_t vecHash = kminmerNorm.hash128();

					if(vecHash < smallestVecHash){
						smallestIsReversed = isReversed;
						smallestVecHash = vecHash;
						smallestI = i;
					}

				}

				if(smallestIsReversed){ //For a circular unitig, the overlaping minimizers will changed depending on the unitig orientation, so we normalize the unitig sequence here
					

					KmerVec vecOriginal = kminmers[smallestI];

					std::reverse(unitig.begin(), unitig.end());
					kminmers = MDBG::minimizersToKminmers(unitig, _graph._kminmerSize);

					for(size_t i=0; i<kminmers.size(); i++){

						const KmerVec& kminmer = kminmers[i];
						const u_int128_t vecHash = kminmer.hash128();

						if(vecHash == smallestVecHash){
							smallestI = i;
							break;
						}

					}

					KmerVec vecNorm = kminmers[smallestI];

					//cout << vecOriginal.toString() << endl;
					//cout << vecNorm.toString() << endl;
					
				}

				unitig.clear();
				unitig = kminmers[smallestI]._kmers;

				for(size_t i=smallestI+1; i<kminmers.size(); i++){
					unitig.push_back(kminmers[i]._kmers[kminmers[i]._kmers.size()-1]);
				}

				for(size_t i=0; i<smallestI; i++){
					unitig.push_back(kminmers[i]._kmers[kminmers[i]._kmers.size()-1]);
				}

				startNode = kminmers[smallestI].normalize();
				endNode = kminmers[smallestI].normalize();

				//cout << endl;
				//for(size_t i=0; i<unitig.size(); i++){
				//	cout << unitig[i] << endl;
				//}
				
				//cout << "Found circular unitig: " << unitig.size() << endl;

				//}

				//return;
				/*
				#pragma omp critical(unitig)
				{
					cout << "Found circular unitig: " << unitig.size() << endl;

					for(size_t i=0; i<unitig.size(); i++){
						cout << unitig[i] << endl;
					}
					
				}

				//deterministic circular forward
				//sourceNode = smallestNode;

				startNode = smallestNode;
				endNode = smallestNode;
				
				unitig = smallestNode._kmers;


				while(true){

					vector<KmerVec> successors;
					_graph.getSuccessors(endNode, successors);

					const KmerVec& successor = successors[0];

					if(successor == smallestNode){ //Circular
						break; 
					}

					endNode = successor;
					
					unitig.push_back(endNode._kmers[endNode._kmers.size()-1]);
				}


				#pragma omp critical(unitig)
				{
					
					cout << "Found circular unitig 2: " << unitig.size() << endl;
					for(size_t i=0; i<unitig.size(); i++){
						cout << unitig[i] << endl;
					}
				}
				*/
			}
			
			//cout << "e" << endl;

			KmerVec startNodeRC = endNode.reverse();
			KmerVec endNodeRC = startNode.reverse();

			//cout << unitig.size() << endl;
			//getchar();

			//if(isCircular){
			//	cout << "skip circular" << endl;
			//	return;
			//}

			#pragma omp critical(unitig)
			{


				bool isValid = true; 
				
				
				if(isNodeUnitigged(startNode)){
					isValid = false;
				}
				else if(isNodeUnitigged(endNode)){
					isValid = false;
				}
				else if(isNodeUnitigged(startNodeRC)){
					isValid = false;
				}
				else if(isNodeUnitigged(endNodeRC)){
					isValid = false;
				}


				if(isValid){

					for(size_t i=0; i<unitig.size()-_graph._kminmerSize+1; i++){
						KmerVec vec;
						for(size_t j=0; j<_graph._kminmerSize; j++){
							vec._kmers.push_back(unitig[i+j]);
						}
						
						setNodeUnitigged(vec);
						setNodeUnitigged(vec.reverse());
					}


					//if(isCircular){
						//cout << "Write circular unitig: " << unitig.size() << endl;
						//for(size_t i=0; i<unitig.size(); i++){
						//	cout << unitig[i] << endl;
						//}
					//}

					_graph.dumpUnitigNode(_graph._unitigIndex, unitig, _graph._unitigGraphFile_nodes);

					_graph._unitigIndex += 2;
					
				}

			}
			
			
		}

	};

	bool isPalindrome(const vector<MinimizerType>& minimizers) const{
		//vector<u_int64_t> kmers = _kmers;
		return  equal(minimizers.begin(), minimizers.begin() + minimizers.size()/2, minimizers.rbegin());
	}
	/*
	class RemapUnitigFunctor {

		public:

		CreateMdbg& _graph;
		ofstream& _outputFile;
		vector<UnitigType>& _unitigNameMap;
		bool _write;

		RemapUnitigFunctor(CreateMdbg& graph, ofstream& outputFile, vector<UnitigType>& unitigNameMap, bool write) : _graph(graph), _outputFile(outputFile), _unitigNameMap(unitigNameMap), _write(write){
		}

		RemapUnitigFunctor(const RemapUnitigFunctor& copy) : _graph(copy._graph), _outputFile(copy._outputFile), _unitigNameMap(copy._unitigNameMap), _write(copy._write){
		}

		~RemapUnitigFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {


			const ReadType readIndex = kminmerList._readIndex;
			
			UnitigType unitigName = _unitigNameMap[readIndex];

			//if(unitigGraph->_unitigs[unitigName] == nullptr) continue;

			const vector<MinimizerType>& minimizers = kminmerList._readMinimizers;
			u_int32_t size = minimizers.size();

			_graph._unitigName_to_minimizers[unitigName] = minimizers;
			//cout << unitigName << " " << minimizers.size() << endl;

			//if(unitigName == 195){
			//	cout << minimizers.size() << endl;
			//}

			unitigName *= 2;

			if(_write){
				_outputFile.write((const char*)&size, sizeof(size));
				_outputFile.write((const char*)&minimizers[0], size * sizeof(MinimizerType));
				_outputFile.write((const char*)&unitigName, sizeof(unitigName));
			}
			

			//if(minimizers.size() == 0){
			//	cout << "derp 2: " << unitigName << " " << minimizers.size() << endl;
			//}

		}
	};
	*/

	int longestOverlap2(const vector<MinimizerType>& list1, const bool isEdgeNode1, const vector<MinimizerType>& list2, const bool isEdgeNode2) {

		if(list1.size() == _kminmerSizePrev && list2.size() == _kminmerSizePrev) return _kminmerSizePrev-1;
		if(isEdgeNode1 || isEdgeNode2) return _kminmerSize-1;
		return _kminmerSizePrev-1;
	}

	int longestOverlap(const UnitigType unitigIndex1, const vector<MinimizerType>& list1, const bool isEdgeNode1, const UnitigType unitigIndex2, const vector<MinimizerType>& list2, const bool isEdgeNode2) {

		return longestOverlap2(list1, isEdgeNode1, list2, isEdgeNode2);

		int x = 0;

		//if(list1.size() == _kminmerSizePrev && list2.size() == _kminmerSizePrev) return _kminmerSizePrev-1;
		//if(isEdgeNode1 || isEdgeNode2) return _kminmerSize-1;
		//return _kminmerSizePrev-1;
		//cout << "over1: " << endl;
		//for(auto m : list1){
		//	cout << m << endl;
		//}
		//cout << "over2: " << endl;
		//for(auto m : list2){
		//	cout << m << endl;
		//}

		size_t n = list1.size();
		size_t m = list2.size();
		size_t maxk = min(n, m);

		for (size_t k = maxk; k >= 1; --k) {
			bool ok = true;
			for (size_t i = 0; i < k; ++i) {
				if (list1[n - k + i] != list2[i]) {
					ok = false;
					break;
				}
			}
			if (ok){
				/*
				cout << _kminmerSize << " " << k << "    " << list1.size() << " " << list2.size() << endl;
				if(k == _kminmerSize-1 || k == _kminmerSize-2){
					if(k == _kminmerSize-1){
					}
				}
				else{
					cout << "derp" << endl;
					getchar();
				}

				return k;
				*/
				x = k;
				break;
			}
		}


		int x2 = longestOverlap2(list1, isEdgeNode1, list2, isEdgeNode2);

		if(unitigIndex1 == unitigIndex2 && list1.size() > 1000){

			cout << endl;
			cout << endl;
			cout << endl;
			cout << endl;
			cout << endl;
			cout << endl;
			cout << "lala: " << list1.size() << endl;

			for(size_t i=0; i<100; i++){
				cout << list1[i] << endl;
			}
			cout << endl;
			for(size_t i=list1.size()-100; i<list1.size(); i++){
				cout << list1[i] << endl;
			}
			cout << x2 << endl;
		}


		/*
		if(x != x2){

			cout << "over1: " << endl;
			for(auto m : list1){
				cout << m << endl;
			}
			cout << "over2: " << endl;
			for(auto m : list2){
				cout << m << endl;
			}


			cout << UnitigGraph2::unitigIndexToString(unitigIndex1) << " -> " << UnitigGraph2::unitigIndexToString(unitigIndex2) << endl;
			cout << x << " " << x2 << endl;
			//getchar();
		}
		*/

		//cout << "none" << endl;
		return x2;
	}

	/*
	void loadReference(size_t kminmerSize){

		_debug_vecHash_to_referencePosition.clear();
		_debug_reference_vecHashs.clear();

		_debug_kminmerSize = kminmerSize;

		cout << "Debug: load reference" << endl;
		const string filename = "/pasteur/appa/scratch/gbenoit/data/genomes/genomes/201/GCF_000816365.1_ASM81636v1_genomic.fna";


		auto fp = std::bind(&CreateMdbg::extract_truth_kminmers_bin_read, this, std::placeholders::_1);
		ReadParser readParser(filename, true, false);
		readParser.parse(fp);
	}


	vector<u_int128_t> _debug_reference_vecHashs;
	unordered_map<u_int128_t, u_int32_t> _debug_vecHash_to_referencePosition;
	size_t _debug_kminmerSize;

	void extract_truth_kminmers_bin_read(const Read& read){

		//ottalSize += strlen(read->seq.s);

		unordered_set<MinimizerType> isRepetitiveMinimizers = Commons::loadRepetitiveMinimizers(_outputDir + "/repetitiveMinimizers.bin", 0);
		
		EncoderRLE _encoderRLE;
		MinimizerParser* _minimizerParser = new MinimizerParser(_minimizerSize, _minimizerDensity, isRepetitiveMinimizers);

		u_int64_t readIndex = read._index;
		string seq = read._seq; //.substr(0, 100);

		string rleSequence;
		vector<u_int64_t> rlePositions;
		_encoderRLE.execute(seq.c_str(), seq.size(), rleSequence, rlePositions, false);

		vector<MinimizerType> minimizers;
		vector<u_int32_t> minimizerPos;
		vector<u_int8_t> minimizerDirections;
		_minimizerParser->parse(rleSequence, minimizers, minimizerPos, minimizerDirections);

		const vector<KmerVec>& kminmers = MDBG::minimizersToKminmers(minimizers, _debug_kminmerSize);


		for(size_t i=0; i<kminmers.size(); i++){

			const KmerVec& kminmer = kminmers[i];
			const u_int128_t vecHash = kminmer.normalize().hash128();
			_debug_vecHash_to_referencePosition[vecHash] = i;
			_debug_reference_vecHashs.push_back(vecHash);
		}

		


	}

	void writeReferencePositions(UnitigGraph2* unitigGraph, const string outputFilename){

		ofstream outputFile(outputFilename);
		outputFile << "Name,Pos" << endl;


		unordered_map<u_int32_t, vector<u_int32_t>> vecHash_to_unitigName;
		
		for(UnitigGraph2::UnitigNode* node : unitigGraph->_unitigs){

			if(node == nullptr) continue;

			vector<MinimizerType> minimizers = _unitigName_to_minimizers[node->_unitigName];

			//cout << node->_unitigName << " " << minimizers.size() << endl;

			const vector<KmerVec>& kminmers = MDBG::minimizersToKminmers(minimizers, _debug_kminmerSize);

			bool isError = false;

			for(size_t i=0; i<kminmers.size(); i++){

				const KmerVec& kminmer = kminmers[i];
				const u_int128_t vecHash = kminmer.normalize().hash128();
				vecHash_to_unitigName[vecHash].push_back(node->_unitigName);

				if(_debug_vecHash_to_referencePosition.find(vecHash) == _debug_vecHash_to_referencePosition.end()){
					isError = true;
					//cout << "utg" << node->_unitigName << " is erroneous " << i << " " << kminmer.toString() << endl;
				}
			}

			if(isError){

				//cout << "utg" << node->_unitigName << endl;
				//for(MinimizerType m : minimizers){
				//	cout << "\t" << m << endl;
				//}

				//outputFile << "utg" << node->_unitigName << "," << "red" << endl;
			}
			else{
				//outputFile << "utg" << node->_unitigName << "," << "green" << endl;
			}


		}

		
		u_int32_t referencePos = 0;
		u_int32_t prevUnitigName = -1;

		for(size_t i=0; i<_debug_reference_vecHashs.size(); i++){

			const u_int128_t vecHash = _debug_reference_vecHashs[i];
			if(vecHash_to_unitigName.find(vecHash) == vecHash_to_unitigName.end()) continue;

			for(const u_int32_t unitigName : vecHash_to_unitigName[vecHash]){

				//u_int32_t unitigName = ;
				
				bool isNew = false;

				if(unitigName != prevUnitigName){
					
					isNew = true;
					outputFile << "utg" << unitigName << "," << referencePos << endl;

				}

				if(isNew){
					prevUnitigName = unitigName;
					referencePos += 1;
				}
			}

		}
		


		
		outputFile.close();
	}
	*/

	//bool isEdgeSupported(UnitigGraph2* unitigGraph, const UnitigType unitigIndexPred, const vector<MinimizerType>& predMinimizers, const UnitigType unitigIndexSucc, const vector<MinimizerType>& succMinimizers, const size_t kminmerSize, auto& kminmers){
	bool isEdgeSupported(const u_int128_t& vecHash, auto& kminmers){
		/*
		//cout << unitigIndexPred << " " << unitigIndexSucc << " " << predMinimizers.size() << endl;
		bool print = false;
		//if(UnitigGraph2::unitigIndex_to_unitigName(unitigIndexPred) == 88 && UnitigGraph2::unitigIndex_to_unitigName(unitigIndexSucc) == 62) print = true;
		//if(UnitigGraph2::unitigIndex_to_unitigName(unitigIndexPred) == 88 && UnitigGraph2::unitigIndex_to_unitigName(unitigIndexSucc) == 85) print = true;


		UnitigGraph2::UnitigNode* predNode = unitigGraph->_unitigs[UnitigGraph2::unitigIndex_to_unitigName(unitigIndexPred)];
		UnitigGraph2::UnitigNode* succNode = unitigGraph->_unitigs[UnitigGraph2::unitigIndex_to_unitigName(unitigIndexSucc)];

		int overlapSize = longestOverlap(unitigIndexPred, predMinimizers, predNode->_isEdgeNode, unitigIndexSucc, succMinimizers, succNode->_isEdgeNode);


		if(print) cout << "\t\tOverlap size: " << overlapSize << " " << kminmerSize << endl; 
		//cout << overlapSize << endl;

		//if(overlapSize <= 2){
		//	cout << "Invalid edge (self loop)" << endl;
		//	return false;
		//}

		//if(unitigIndexPred == unitigIndexSucc && overlapSize == predMinimizers.size()){
		//	cout << "Invalid edge (self loop)" << endl;
		//	return false;
		//}

		if(overlapSize >= kminmerSize) {
			cout << "Issue in isEdgeSupported: overlapSize == kminmerSize" << endl;
			return false;
		}

		vector<MinimizerType> doublet;

		for(size_t i=0; i<overlapSize; i++){
			doublet.push_back(succMinimizers[i]);
		}

		int left =0;
		int right = 0;

		while(true){

			if(doublet.size() == kminmerSize) break;

			if(overlapSize+right >= succMinimizers.size()){
				cout << overlapSize << endl;
				cout << "getDoublet derp 1" << endl;
				return false;
				//getchar();
			}

			doublet.push_back(succMinimizers[overlapSize+right]);
			right += 1;

			if(doublet.size() == kminmerSize) break;

			if(predMinimizers.size()-overlapSize-1-left >= predMinimizers.size()){
				cout << overlapSize << endl;
				cout << "getDoublet derp 2" << endl;
				return false;
				//getchar();
			}

			doublet.insert(doublet.begin(), predMinimizers[predMinimizers.size()-overlapSize-1-left]);
			left += 1;
			
			if(doublet.size() == kminmerSize) break;
		}

		//cout << "\tDoublet:" << endl;
		//for(MinimizerType m : doublet) {
		//	cout << "\t\t" << m << endl;
		//}


		KmerVec kminmer;
		kminmer._kmers = doublet;

		if(print){
			cout << kminmer.toString() << endl;
			cout << (kminmers.find(kminmer.hash128()) != kminmers.end()) << endl;
			cout << kminmers[kminmer.hash128()]._abundance << endl;
		}

		//cout << kminmer.toString() << endl;
		
		kminmer = kminmer.normalize();


		//
		//getchar();
		*/

		if(kminmers.find(vecHash) == kminmers.end()) return false;

		return kminmers[vecHash] >= 2;

	}

	void getDoublet2(UnitigGraph2* unitigGraph, const UnitigType unitigIndexPred, const UnitigType unitigIndexSucc, vector<MinimizerType>& minimizers){

		minimizers.clear();

		vector<MinimizerType> succMinimizers;
		UnitigGraph2::getStartNode(unitigIndexSucc, _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(unitigIndexSucc)], _kminmerSizePrev, succMinimizers);

		vector<MinimizerType> predMinimizers;
		UnitigGraph2::getEndNode(unitigIndexPred, _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(unitigIndexPred)], _kminmerSizePrev, predMinimizers);

		//minimizers = succMinimizers;
		minimizers.push_back(predMinimizers[0]);
		for(const MinimizerType m : succMinimizers){
			minimizers.push_back(m);
		}

		//cout << minimizers.size() << endl;

	}
	
	const vector<MinimizerType> getDoublet(UnitigGraph2* unitigGraph, const UnitigType unitigIndexPred, const vector<MinimizerType>& predMinimizers, const UnitigType unitigIndexSucc, const vector<MinimizerType>& succMinimizers, const size_t kminmerSize){

		UnitigGraph2::UnitigNode* predNode = unitigGraph->_unitigs[UnitigGraph2::unitigIndex_to_unitigName(unitigIndexPred)];
		UnitigGraph2::UnitigNode* succNode = unitigGraph->_unitigs[UnitigGraph2::unitigIndex_to_unitigName(unitigIndexSucc)];
		//int overlapSize = longestOverlap(predMinimizers, succMinimizers);
		int overlapSize = longestOverlap(unitigIndexPred, predMinimizers, predNode->_isEdgeNode, unitigIndexSucc, succMinimizers, succNode->_isEdgeNode);

		vector<MinimizerType> doublet;

		for(size_t i=0; i<overlapSize; i++){
			doublet.push_back(succMinimizers[i]);
		}

		int left =0;
		int right = 0;

		while(true){

			if(doublet.size() == kminmerSize) break;

			if(overlapSize+right >= succMinimizers.size()){
				cout << "getDoublet derp 1" << endl;
				//getchar();
			}

			doublet.push_back(succMinimizers[overlapSize+right]);
			right += 1;

			if(doublet.size() == kminmerSize) break;

			if(predMinimizers.size()-overlapSize-1-left >= predMinimizers.size()){
				cout << "getDoublet derp 2" << endl;
				//getchar();
			}

			doublet.insert(doublet.begin(), predMinimizers[predMinimizers.size()-overlapSize-1-left]);
			left += 1;
			
			if(doublet.size() == kminmerSize) break;
		}

		return doublet;
	}

	vector<MinimizerType> kminmersToMinimizers(const vector<KmerVec> vecs){

		vector<MinimizerType> minimizers;

		for(size_t i=0; i<vecs.size(); i++){

			const KmerVec& vec = vecs[i];

			if(i==0){

				for(const MinimizerType minimizer : vec._kmers){
					minimizers.push_back(minimizer);
				}
			}
			else{
				minimizers.push_back(vec._kmers[vec._kmers.size()-1]);
			}
		}

		return minimizers;
	}

	void writeUnitigSequences(UnitigGraph2* unitigGraph, const vector<UnitigType>& currentUnitigName_to_newUnitigName){


		ofstream unitigGraphFile_nodes(_outputDir + "/unitigGraph.nodes.bin");

		//for(UnitigGraph2::UnitigNode* node : unitigGraph->_unitigs){


		#pragma omp parallel for num_threads(_nbCores)
		for(size_t i=0; i<unitigGraph->_unitigs.size(); i++){

			UnitigGraph2::UnitigNode* node = unitigGraph->_unitigs[i];
			if(node == nullptr) continue;

			vector<UnitigType> unitigs;
			if(node->_unitigMerge.size() == 0){
				unitigs.push_back(UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, false));    
			}
			else{
				unitigs = node->_unitigMerge;
			}

			if(node->_isReversed){
				vector<UnitigType> unitigsRC;
				UnitigGraph2::reverseComplementUnitigs(unitigs, unitigsRC);

				unitigs = unitigsRC;
			}

			bool isCircular = false;
			//u_int32_t size = unitigs.size();

			//outputContigFile.write((const char*)&size, sizeof(size));
			//outputContigFile.write((const char*)&isCircular, sizeof(isCircular));
			//outputContigFile.write((const char*)&unitigs[0], size * sizeof(UnitigType));

			//nbUnitigsWritten += 1;

			vector<MinimizerType> minimizers;
			unitigsToMinimizers(unitigs, minimizers);

			writeUnitigSequence(unitigGraphFile_nodes, currentUnitigName_to_newUnitigName[node->_unitigName], minimizers);

		}

		unitigGraphFile_nodes.close();
	}



	void unitigsToMinimizers(const vector<UnitigType>& unitigs, vector<MinimizerType>& outputMinimizers){

		outputMinimizers.clear();
		vector<MinimizerType> prevMinimizers;

		for(size_t i=0; i<unitigs.size(); i++){

			UnitigType unitigIndex = unitigs[i];


			bool isReversed;
			const UnitigType& unitigName = UnitigGraph2::unitigIndex_to_unitigName(unitigIndex, isReversed);

			//cout << "\t" << unitigName << " " << isReversed << endl;
			//cout << "\t" << unitigName << " " << _parent._unitigName_to_minimizers.size() << " " << endl;
			//cout << "\t" << _parent._unitigName_to_minimizers[unitigName].size() << " " << endl;

			vector<MinimizerType> minimizers = _unitigName_to_minimizers[unitigName];

			if(isReversed){
				std::reverse(minimizers.begin(), minimizers.end());
			}

			if(i==0){

				outputMinimizers = minimizers;
			}
			else{

				int overlapSize = getUnitigOverlap(prevMinimizers, minimizers);
				
				outputMinimizers.insert(outputMinimizers.end(), minimizers.begin()+overlapSize, minimizers.end());

			}

			prevMinimizers = minimizers;
		}

	}
	
	int getUnitigOverlap(const vector<MinimizerType>& list1, const vector<MinimizerType>& list2) {

		if(list1.size() == _kminmerSize && list2.size() == _kminmerSize){
			if(list1 == list2) return _kminmerSize;
		}

		return _kminmerSize - 1;
	
	}


	void writeUnitigSequence(ofstream& outputFile, const UnitigType unitigName, const vector<MinimizerType>& minimizers){
	
		#pragma omp critical(writeUnitigSequence)
		{

			if(minimizers.size() == _kminmerSizePrev){
				cout << "Remaining small unitigs" << endl;
				getchar();
			}
			u_int32_t size = minimizers.size();
			UnitigType unitigIndex = unitigName*2;

			outputFile.write((const char*)&size, sizeof(size));
			outputFile.write((const char*)&minimizers[0], size * sizeof(MinimizerType));
			outputFile.write((const char*)&unitigIndex, sizeof(unitigIndex));
			_nbUnitigNodes += 1;
			
			//cout << "write unitig: " << unitigIndex << " " << minimizers.size() << endl;
		}

	}

	/*
	void countSolidKminmerInGraph(UnitigGraph2* unitigGraph){

		for(UnitigGraph2::UnitigNode* node : unitigGraph->_unitigs){

			if(node == nullptr) continue;

			const vector<MinimizerType>& minimizers = _unitigName_to_minimizers[node->_unitigName];

			for(const KmerVec& kminmer : MDBG::minimizersToKminmers(minimizers, _kminmerSize)){

				const u_int128_t hash = kminmer.normalize().hash128();
				
				if(_mdbgNodesLight.find(hash) != _mdbgNodesLight.end()){
					if(_kminmerSize == _kminmerSizeFirst+1){
						if (_mdbgNodesLight[hash]._abundance > 10 && _mdbgNodesLight[hash]._index <= 1) {
							cout << "Found chimeric kminmer" << endl;
							return true;
						}
					}
				}
				else{
					cout << "Unitig has kminmers not in set" << endl;
					return true;
				}
			}

		}
	}
	*/

	class PartitionFile{

		public:

		//gzFile _file;
		ofstream _file;
		omp_lock_t _mutex;

		PartitionFile(u_int32_t partition, const string& partitionFilename){
			_file = ofstream(partitionFilename.c_str());

			omp_init_lock(&_mutex);
		}

		~PartitionFile(){
			_file.close();
			omp_destroy_lock(&_mutex);
		}
	};



	class KminmerCounter{

		public:

		CreateMdbg& _parent;
		vector<ofstream> _partitionFiles;
		vector<omp_lock_t> _mutex;
		int _nbPartitions;
		string _partitionDir;
		ofstream _outputFile;
		ofstream _outputFileAbundance;
		u_int64_t _nbKminmers;
		u_int64_t _nbSolidKminmers;

		vector<KmerVec> _kminmersCache;

		KminmerCounter(CreateMdbg& parent) : _parent(parent){


			_nbPartitions = parent._nbPartitions;

			_partitionDir = _parent._outputDir + "/kminmerCounter_tmp/";
			if(!fs::exists(_partitionDir)){
				fs::create_directories(_partitionDir); 
			}

				
			_mutex.resize(_nbPartitions);

			for(size_t i=0; i<_mutex.size(); i++){
				omp_init_lock(&_mutex[i]);
			}

			_nbKminmers = 0;
			_nbSolidKminmers = 0;
		}

		~KminmerCounter(){
			for(size_t i=0; i<_mutex.size(); i++){
				omp_destroy_lock(&_mutex[i]);
			}
		}

		void execute(){

			
			partitionKminmers();
			dereplicatePartitions();
			
			fs::remove_all(_partitionDir);
			
		}
		
		//string getOutputFilename(){
		//	return _parent._outputDir + "/kminmerData_min.txt";
		//}

		string getPartitionFilename(int partition){
			return _partitionDir + "/part_" + to_string(partition) + ".bin";
		}

		void partitionKminmers(){

			_partitionFiles.resize(_nbPartitions);

			for(size_t i=0; i<_partitionFiles.size(); i++){
				_partitionFiles[i] = ofstream(getPartitionFilename(i));
			}

			KminmerParserParallel parser1(_parent._outputDir + "/read_data_corrected.txt", 0, _parent._kminmerSize, false, false, _parent._nbCores);
			parser1.parse(KminmerParserFunctor(*this));

			if(!_parent._isFirstPass){
				KminmerParserParallel parser2(_parent._outputDir + "/unitig_data.txt", 0, _parent._kminmerSize, false, false, _parent._nbCores);
				parser2.parse(KminmerParserFunctor(*this));
			}


			for(ofstream& outputFile : _partitionFiles){
				outputFile.close();
			}
			
		}


		class KminmerParserFunctor {

			public:

			KminmerCounter& _parent;

			KminmerParserFunctor(KminmerCounter& parent) : _parent(parent){
			}

			KminmerParserFunctor(const KminmerParserFunctor& copy) : _parent(copy._parent){
				
			}

			~KminmerParserFunctor(){
			}

			void operator () (const KminmerList& kminmerList) {


				u_int64_t readIndex = kminmerList._readIndex;
				const vector<MinimizerType>& readMinimizers = kminmerList._readMinimizers;
				//const vector<KmerVec>& kminmers = kminmerList._kminmers;
				const vector<ReadKminmerComplete>& kminmersInfos = kminmerList._kminmersInfo;

				for(size_t i=0; i<kminmersInfos.size(); i++){
					


					const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];
					const KmerVec& vec = kminmerInfo._vec;
					
					_parent.partitionKminmer(vec);
				}


			}
		};

		void partitionKminmer(const KmerVec& vec){

			u_int128_t vecHash = vec.hash128();
			int partition = vecHash % _nbPartitions;

			omp_set_lock(&_mutex[partition]);

			MDBG::writeKminmer(vec._kmers, _partitionFiles[partition]);

			omp_unset_lock(&_mutex[partition]);
		}

		void dereplicatePartitions(){
			
			//_outputFile = ofstream(getOutputFilename());
			//_outputFileAbundance = ofstream(_parent._outputDir + "/kminmerData_abundance.txt");

			for(size_t i=0; i<_nbPartitions; i++){
				dereplicatePartition(i);
			}

			//_outputFile.close();
			//_outputFileAbundance.close();
		}

		struct KminmerAbundance{
			KmerVec _vec;
			u_int32_t _abundance;
		};

		void dereplicatePartition(int partition){

			auto start2 = high_resolution_clock::now();

			Logger::get().debug() << "\tDereplicate partition: " << partition;
			auto start = high_resolution_clock::now();
			//unordered_set<KmerVec> truth;
			//int nbKminmers = 0;

			//vector<vector<MinimizerType>> kminmers2;
			//vector<KmerVec> kminmers;

			ifstream partitionFile(getPartitionFilename(partition));
			vector<MinimizerType> minimizers(_parent._kminmerSize, 0);

			size_t i=0;

			while (true) {


				bool isEOF = MDBG::readKminmer2(minimizers, partitionFile);
				
				if(isEOF) break;

				if(i < _kminmersCache.size()){
					for(size_t j=0; j<minimizers.size(); j++){
						_kminmersCache[i]._kmers[j] = minimizers[j];
					}
				}
				else{
					
					KmerVec vec;
					vec._kmers = minimizers;

					_kminmersCache.push_back(vec);
				}

				i += 1;

				//if(vec._kmers.size() != 4) cout << "allo??" << endl;
				//kminmers2.push_back(minimizers);

				//if(kminmers.size() > 100000) break;
			}

			u_int64_t size = i;
			//cout << _kminmersCache.size() << " " << size << endl;

			partitionFile.close();

			Logger::get().debug() << "\tLoading: " << _kminmersCache.size() << " " << duration_cast<seconds>(high_resolution_clock::now() - start).count(); 

			start = high_resolution_clock::now();
			Commons::sortParallel(_kminmersCache, size, _parent._nbCores);
			Logger::get().debug() << "\tSorting: " << duration_cast<seconds>(high_resolution_clock::now() - start).count();
			
			
			start = high_resolution_clock::now();

			vector<KminmerAbundance> buffer;

			KmerVec prevKminmer;
			u_int32_t abundance;

			for (size_t i=0; i<size; i++) {
				
				const KmerVec& kminmer = _kminmersCache[i];

				//truth.insert(kminmer);

				//cout << nbKminmers << " " << truth.size() << endl;
				if(i == 0){

					prevKminmer = kminmer;
					abundance = 1;
				}
				else if(kminmer != prevKminmer){
					
					buffer.push_back({prevKminmer, abundance});

					if(buffer.size() > 10000){
						dumpBuffer(buffer);
						buffer.clear();
					}
					//dumpKminmer(prevKminmer, abundance);

					prevKminmer = kminmer;
					abundance = 1;
				}
				else{
					abundance += 1;
				}

				
			}
			
			//Last kminmer
			buffer.push_back({prevKminmer, abundance});
			dumpBuffer(buffer);
			buffer.clear();
			//dumpKminmer(prevKminmer, abundance);


			Logger::get().debug() << "\tDumping: " << duration_cast<seconds>(high_resolution_clock::now() - start).count();

			Logger::get().debug() << "\tDone: " << _nbKminmers << " " << duration_cast<seconds>(high_resolution_clock::now() - start2).count() << "s " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";
			
		}

		void dumpBuffer(const vector<KminmerAbundance>& kminmerAbundances){


			#pragma omp parallel for num_threads(_parent._nbCores)
			for(size_t i=0; i<kminmerAbundances.size(); i++){
				dumpKminmer(kminmerAbundances[i]._vec, kminmerAbundances[i]._abundance);
			}
		}

		void dumpKminmer(const KmerVec& vec, u_int32_t abundance){

			if(!_parent._isFirstPass){
				abundance = getRefinedAbundance(vec);
			}

			if(abundance <= 1) return;
			if(_parent._isFirstPass && abundance < _parent._minAbundance) return;

			#pragma omp critical(dumpKminmer)
			{
				_nbKminmers += 1;
				MDBG::writeKminmer(vec._kmers, _parent._kminmerFile);

				const u_int128_t vecHash = vec.hash128();
				MDBG::writeKminmerAbundance(vecHash, abundance, _parent._kminmerAbundanceFile);

				_nbSolidKminmers += 1;
			}
			

		}

		/*

		static void sortParallel(vector<KmerVec>& vec, const size_t size, const int nbCores) {

			const u_int64_t data_count = size;
			auto get_block_edge = [data_count](int i, int n) {
				return data_count * i / n;
			};

			//cout << "a" << endl;

			u_int64_t blocks = nbCores;
			#pragma omp parallel num_threads(nbCores)
			{
				//blocks = omp_get_num_threads();
				u_int64_t block = omp_get_thread_num();
				u_int64_t start = get_block_edge(block, blocks);
				u_int64_t finish = get_block_edge(block + 1, blocks);
			
				
				std::sort(std::begin(vec) + start, std::begin(vec) + finish, [](const auto& a, const auto& b){
					return a < b;
				});
				
			}

			//cout << "b" << endl;
			//cout << "merging" << endl;
			for (int merge_step = 1; merge_step < blocks; merge_step *= 2) {
				#pragma omp parallel for num_threads(nbCores)
				for (u_int64_t i = 0; i < blocks; i += 2 * merge_step) {
					u_int64_t start = get_block_edge(i, blocks);
					u_int64_t mid = std::min(get_block_edge(i + merge_step, blocks), data_count);
					u_int64_t finish = std::min(get_block_edge(i + 2 * merge_step, blocks), data_count);
					if (mid < finish)
						std::inplace_merge(std::begin(vec) + start, std::begin(vec) + mid, std::begin(vec) + finish, [](const auto& a, const auto& b){
							return a < b;
						});
				}
			}
			
			//cout << "c" << endl;
			//cout << "done" << endl;
			
		}
		*/


		float getRefinedAbundance(const KmerVec& vec){

			vector<u_int64_t> rlePositions;
			vector<u_int32_t> minimizers_pos;//(minimizers.size());
			vector<KmerVec> kminmers; 
			vector<ReadKminmer> kminmersInfo;
			MDBG::getKminmers(0, _parent._kminmerSizePrev, vec._kmers, minimizers_pos, kminmers, kminmersInfo, rlePositions, 0, false);

			//string s = "";

			u_int32_t minAbundance = -1;
			//u_int32_t maxAbundance = 0;
			//float sum = 0;
			//float n = 0;
			//vector<float> abundances;
			//cout << kminmers.size() << " " << _kminmerAbundances.size() << endl;
			for(const KmerVec& vec : kminmers){

				//cout << vec._kmers.size() << endl;
				//cout << vec.toString() << endl;

				const u_int128_t vecHash = vec.hash128();

				//cout << (_mdbgInit->_dbg_nodes.find(vec) != _mdbgInit->_dbg_nodes.end()) << endl;
				if(_parent._kminmerAbundances.find(vecHash) != _parent._kminmerAbundances.end()){

					//u_int32_t abundance = _kminmerAbundances[vec];
					//const KminmerAbundanceMapNode& node = _kminmerAbundances[vec];
					u_int32_t abundance = _parent._kminmerAbundances[vecHash];
					//if(nodeName == 17430) cout << abundance << endl;

					if(abundance == 0){
						minAbundance = 1;
						break;
					}

					if(abundance < minAbundance){
						minAbundance = abundance;
					}

					//if(abundance > maxAbundance){
					//	maxAbundance = abundance;
					//}

					//s += to_string(abundance) + " ";
					//cout << abundance << " ";
					//abundances.push_back(abundance);

					//sum += abundance;
					//n += 1;
				}
				else{
					minAbundance = 1;
					break;
				}
			}

			/*
			#pragma omp critical
			{

				if(_graph._kminmerSize >= 18){
					cout << endl;
					cout << _graph._kminmerSize << " " << _graph._kminmerSizeFirst << endl;
					cout << s << endl;
					cout << minAbundance << endl;
					if(maxAbundance > 65) getchar();
				}
			}
			*/

			return minAbundance;
		}

	};


	class EdgeIndexer{

		public:

		CreateMdbg& _parent;
		vector<ofstream> _partitionFiles;
		vector<omp_lock_t> _mutex;
		int _nbPartitions;
		string _partitionDir;
		ofstream _outputFile;
		u_int64_t _nbEdges;
		u_int64_t _checksum;
		vector<u_int128_t> _edgeCache;

		EdgeIndexer(CreateMdbg& parent) : _parent(parent){

			_nbPartitions = parent._nbPartitions;


			_partitionDir = _parent._outputDir + "/edgeIndexer_tmp/";
			if(!fs::exists(_partitionDir)){
				fs::create_directories(_partitionDir); 
			}

				
			_mutex.resize(_nbPartitions);

			for(size_t i=0; i<_mutex.size(); i++){
				omp_init_lock(&_mutex[i]);
			}



			_nbEdges = 0;
			_checksum = 0;
			
		}

		~EdgeIndexer(){
			for(size_t i=0; i<_mutex.size(); i++){
				omp_destroy_lock(&_mutex[i]);
			}
		}

		void execute(){

			partitionEdges();
			dereplicatePartitions();
			
			fs::remove_all(_partitionDir);
		}
		
		string getOutputFilename(){
			return _parent._outputDir + "/edges.bin";
		}

		string getPartitionFilename(int partition){
			return _partitionDir + "/part_" + to_string(partition) + ".bin";
		}

		void partitionEdges(){

			_partitionFiles.resize(_nbPartitions);

			for(size_t i=0; i<_partitionFiles.size(); i++){
				_partitionFiles[i] = ofstream(getPartitionFilename(i));
			}

			size_t i=0;
			ifstream kminmerFile(_parent._outputDir + "/kminmerData_min.txt");

			while (true) {


				vector<MinimizerType> minimizers;
				bool isEOF = MDBG::readKminmer(minimizers, kminmerFile, _parent._kminmerSize);
				
				if(isEOF) break;

				KmerVec vec;
				vec._kmers = minimizers;

				partitionNode(vec);

				i += 1;
				//if(i % 1000000 == 0) cout << "\tPartitionning edges: " << i << endl;

			}

			kminmerFile.close();


			for(ofstream& outputFile : _partitionFiles){
				outputFile.close();
			}
			
		}

		void partitionNode(const KmerVec& vec){

			bool needprint = false;

			bool isReversedPrefix;
			bool isReversedSuffix;
			KmerVec prefix = vec.prefix().normalize(isReversedPrefix);
			KmerVec suffix = vec.suffix().normalize(isReversedSuffix);

			const u_int128_t suffixHash = suffix.hash128();
			const u_int128_t prefixHash = prefix.hash128();

			partitionEdge(suffixHash);
			partitionEdge(prefixHash);
		}

		void partitionEdge(const u_int128_t& vecHash){


			int partition = vecHash % _nbPartitions;

			omp_set_lock(&_mutex[partition]);

			_partitionFiles[partition].write((const char*)&vecHash, sizeof(vecHash));

			omp_unset_lock(&_mutex[partition]);
		}

		void dereplicatePartitions(){
			
			_outputFile = ofstream(getOutputFilename());

			for(size_t i=0; i<_nbPartitions; i++){
				dereplicatePartition(i);
			}

			_outputFile.close();
		}

		void dereplicatePartition(int partition){

			Logger::get().debug() << "\tDereplicate partition: " << partition;
			//unordered_set<u_int128_t> truth;

			//vector<u_int128_t> edges;
			ifstream partitionFile(getPartitionFilename(partition));

			size_t i=0;

			while (true) {

				u_int128_t edgeHash;
				partitionFile.read((char*)&edgeHash, sizeof(edgeHash));

				if(partitionFile.eof()) break;

				//edges.push_back(edgeHash);



				if(i < _edgeCache.size()){
					_edgeCache[i] = edgeHash;
				}
				else{
					_edgeCache.push_back(edgeHash);
				}

				i += 1;
			}

			partitionFile.close();

			u_int64_t size = i;
			
			Commons::sortParallel(_edgeCache, size, _parent._nbCores);

			size_t nbEdges = 0;
			u_int128_t prevEdge = -1;

			for (size_t i=0; i<size; i++) {
				
				const u_int128_t& edge = _edgeCache[i];
				//truth.insert(edge);

				if(i == 0){
					
					_checksum += edge;
					_nbEdges += 1;
					nbEdges += 1;
					_outputFile.write((const char*)&edge, sizeof(edge));
					prevEdge = edge;
				}
				else if(edge != prevEdge){
					
					_checksum += edge;
					_nbEdges += 1;
					nbEdges += 1;
					_outputFile.write((const char*)&edge, sizeof(edge));
					prevEdge = edge;
				}
				
			}

			/*
			u_int64_t nbPositions = 0;

			for(const auto& it : _minimizerLookupTable){
				for(size_t i=it.second.first; i <= it.second.second; i++){
					nbPositions += 1;
				}
			}

			cout << "Lookup table check: " << _allMinimizerPositions.size() << " " << nbPositions << endl;
			if(nbPositions != _allMinimizerPositions.size()){
				cout << "Lookup table issue" << endl;
				getchar();
			}
			*/
			
			Logger::get().debug() << "\tDone: " << nbEdges << " " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";
		}

	};



	class UnitigEdgeIndexer{

		public:

		CreateMdbg& _parent;
		vector<ofstream> _partitionFiles;
		vector<omp_lock_t> _mutex;
		int _nbPartitions;
		string _partitionDir;
		ofstream _outputFile;
		u_int64_t _nbEdges;
		vector<u_int128_t> _edgeCache;

		UnitigEdgeIndexer(CreateMdbg& parent) : _parent(parent){

			_nbPartitions = parent._nbPartitions;


			_partitionDir = _parent._outputDir + "/unitigEdgeIndexer_tmp/";
			if(!fs::exists(_partitionDir)){
				fs::create_directories(_partitionDir); 
			}

				
			_mutex.resize(_nbPartitions);

			for(size_t i=0; i<_mutex.size(); i++){
				omp_init_lock(&_mutex[i]);
			}



			_nbEdges = 0;
			
		}

		~UnitigEdgeIndexer(){
			for(size_t i=0; i<_mutex.size(); i++){
				omp_destroy_lock(&_mutex[i]);
			}
		}

		void execute(){

			partitionEdges();
			dereplicatePartitions();
			
			fs::remove_all(_partitionDir);
		}
		
		string getOutputFilename(){
			return _parent._outputDir + "/unitig_edges.bin";
		}

		string getPartitionFilename(int partition){
			return _partitionDir + "/part_" + to_string(partition) + ".bin";
		}

		void partitionEdges(){

			_partitionFiles.resize(_nbPartitions);

			for(size_t i=0; i<_partitionFiles.size(); i++){
				_partitionFiles[i] = ofstream(getPartitionFilename(i));
			}

			size_t i=0;
			/*
			ifstream kminmerFile(_parent._outputDir + "/kminmerData_min.txt");

			while (true) {


				vector<MinimizerType> minimizers;
				bool isEOF = MDBG::readKminmer(minimizers, kminmerFile, _parent._kminmerSize);
				
				if(isEOF) break;

				KmerVec vec;
				vec._kmers = minimizers;

				partitionNode(vec);

				i += 1;
				if(i % 1000000 == 0) cout << "\tPartitionning edges: " << i << endl;

			}

			kminmerFile.close();
			*/
			ifstream nodeFile(_parent._outputDir + "/unitigGraph.nodes.bin");

			#pragma omp parallel num_threads(_parent._nbCores)
			{

				bool isEOF = false;
				vector<MinimizerType> unitig;
				UnitigType unitigIndex;
				//u_int32_t nodeName;

				while (true) {

					#pragma omp critical(indexEdge)
					{

						u_int32_t size;
						nodeFile.read((char*)&size, sizeof(size));

						isEOF = nodeFile.eof();

						

						if(!isEOF){
							unitig.resize(size);
							nodeFile.read((char*)&unitig[0], size * sizeof(MinimizerType));


							nodeFile.read((char*)&unitigIndex, sizeof(unitigIndex));
							//nodeFile.read((char*)&nodeName, sizeof(nodeName));
						}

						//cout << nodeName << endl;
					}

					if(isEOF) break;
					
					partitionUnitig(unitigIndex, unitig);

				}
			}

			nodeFile.close();


			for(ofstream& outputFile : _partitionFiles){
				outputFile.close();
			}
			
		}

		void partitionUnitig(const u_int32_t unitigIndex, const vector<MinimizerType>& minimizers){

			u_int32_t unitigIndexRC = unitigIndex + 1;
			const vector<KmerVec>& nodes = MDBG::minimizersToKminmers(minimizers, _parent._kminmerSize);

			const KmerVec& startNode = nodes[0];
			const KmerVec& endNode = nodes[nodes.size()-1];
			const KmerVec& startNodeRC = endNode.reverse();
			const KmerVec& endNodeRC = startNode.reverse();


			bool isReversed;
			const KmerVec& startNodeNorm = startNode.normalize(isReversed);
			
			
			if(isReversed){
				partitionUnitigEdge(startNodeNorm, unitigIndexRC);
			}
			else{
				partitionUnitigEdge(startNodeNorm, unitigIndex);
			}

			if(startNode == endNode) return;
			
			const KmerVec& endNodeNorm = endNode.normalize(isReversed);

			//cout << "\tEnd node: " << endNodeNorm.toString() << " " << isReversed << endl;

			if(isReversed){
				partitionUnitigEdge(endNodeNorm, unitigIndexRC);
			}
			else{
				partitionUnitigEdge(endNodeNorm, unitigIndex);
			}

		}

		void partitionUnitigEdge(const KmerVec& vec, u_int32_t unitigIndex){


			bool isReversedPrefix;
			bool isReversedSuffix;
			KmerVec prefix = vec.prefix().normalize(isReversedPrefix);
			KmerVec suffix = vec.suffix().normalize(isReversedSuffix);

			const u_int128_t suffixHash = suffix.hash128();
			const u_int128_t prefixHash = prefix.hash128();

			partitionEdge(suffixHash);
			partitionEdge(prefixHash);
			
		}

		void partitionEdge(const u_int128_t& vecHash){

			int partition = vecHash % _nbPartitions;

			omp_set_lock(&_mutex[partition]);

			_partitionFiles[partition].write((const char*)&vecHash, sizeof(vecHash));

			omp_unset_lock(&_mutex[partition]);
		}

		void dereplicatePartitions(){
			
			_outputFile = ofstream(getOutputFilename());

			for(size_t i=0; i<_nbPartitions; i++){
				dereplicatePartition(i);
			}

			_outputFile.close();
		}

		void dereplicatePartition(int partition){

			Logger::get().debug() << "\tDereplicate partition: " << partition;
			//unordered_set<u_int128_t> truth;

			//vector<u_int128_t> edges;
			ifstream partitionFile(getPartitionFilename(partition));

			size_t i=0;

			while (true) {

				u_int128_t edgeHash;
				partitionFile.read((char*)&edgeHash, sizeof(edgeHash));

				if(partitionFile.eof()) break;

				//edges.push_back(edgeHash);

				if(i < _edgeCache.size()){
					_edgeCache[i] = edgeHash;
				}
				else{
					_edgeCache.push_back(edgeHash);
				}
				
				i += 1;
			}

			partitionFile.close();

			u_int64_t size = i;
			
			Commons::sortParallel(_edgeCache, size, _parent._nbCores);

			size_t nbEdges = 0;
			u_int128_t prevEdge = -1;

			for (size_t i=0; i<size; i++) {
				
				const u_int128_t& edge = _edgeCache[i];
				//truth.insert(edge);

				if(i == 0){
					
					_nbEdges += 1;
					nbEdges += 1;
					_outputFile.write((const char*)&edge, sizeof(edge));
					prevEdge = edge;
				}
				else if(edge != prevEdge){
					
					_nbEdges += 1;
					nbEdges += 1;
					_outputFile.write((const char*)&edge, sizeof(edge));
					prevEdge = edge;
				}
				
			}
			
			Logger::get().debug() << "\tDone: " << nbEdges << " " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";
		}

	};

	ankerl::unordered_dense::map<u_int128_t, u_int32_t> _kminmerAbundancesRescue;
	u_int64_t _nbRescuedKminmers;

	void rescueKminmers(){


		Logger::get().debug() << "\tLoading solid kminmer abundance";
		auto start = high_resolution_clock::now();
		_nbRescuedKminmers = 0;
		loadAbundanceRescuing();
		Logger::get().debug() << "\t\tDone: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";




		Logger::get().debug() << "\tRescuing kminmers";
		start = high_resolution_clock::now();

		KminmerParserParallel parser(_outputDir + "/read_data_corrected.txt", 0, _kminmerSize, false, false, _nbCores);
		parser.parse(RescueKminmerFunctor(*this));

		ankerl::unordered_dense::map<u_int128_t, u_int32_t>().swap(_kminmerAbundancesRescue);
		
		Logger::get().debug() << "\t\tDone: " << " " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";
	}

	void loadAbundanceRescuing(){

		ifstream kminmerAbundanceFile(_outputDir + "/kminmerData_abundance.txt");

		while (true) {

			u_int128_t vecHash;
			AbundanceType abundance;

			bool iseof = MDBG::readKminmerAbundance(vecHash, abundance, kminmerAbundanceFile);

			if(iseof) break;

			if(abundance == 1) continue; //kminmerData_abundance should only contain solid kminmer with abundance > 1 at this point
			_kminmerAbundancesRescue[vecHash] = abundance; //Using u_int64_t as key?


		}

	}

	
	class RescueKminmerFunctor {

		public:

		CreateMdbg& _parent;

		RescueKminmerFunctor(CreateMdbg& parent) : _parent(parent){

		}

		RescueKminmerFunctor(const RescueKminmerFunctor& copy) : _parent(copy._parent){
		}

		~RescueKminmerFunctor(){
		}

	
		void operator () (const KminmerList& kminmerList) {


			u_int64_t readIndex = kminmerList._readIndex;
			
			if(readIndex % 1000000 == 0) Logger::get().debug() << "\tRescuing kminmers: " << readIndex;

			const vector<MinimizerType>& readMinimizers = kminmerList._readMinimizers;
			const vector<ReadKminmerComplete>& kminmersInfos = kminmerList._kminmersInfo;

			bool allAbundanceOne = true;
			vector<u_int32_t> abundances;
			
			for(size_t i=0; i<kminmersInfos.size(); i++){
				
				const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];
				const KmerVec& vec = kminmerInfo._vec;
				const u_int128_t vecHash = vec.hash128();
				
				if(_parent._kminmerAbundancesRescue.find(vecHash) != _parent._kminmerAbundancesRescue.end()){
					abundances.push_back(_parent._kminmerAbundancesRescue[vecHash]);
					allAbundanceOne = false;
				}
				else{
					abundances.push_back(1);
				}

			}

			u_int32_t median = Utils::compute_median(abundances);
			double cutoff = median * 0.1f;

			
			if(cutoff > 1) return;
			if(allAbundanceOne) return;

			for(size_t i=0; i<kminmersInfos.size(); i++){
				
				const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];
				const KmerVec& vec = kminmerInfo._vec;
				const u_int128_t vecHash = vec.hash128();

				if(_parent._kminmerAbundancesRescue.find(vecHash) != _parent._kminmerAbundancesRescue.end()) continue; // ici on peut juste check le vecteur d'abondance ci dessus
 
				u_int32_t abundance = 1;

				#pragma omp critical(RescueKminmerFunctor)
				{

					MDBG::writeKminmer(vec._kmers, _parent._kminmerFile);
					MDBG::writeKminmerAbundance(vecHash, abundance, _parent._kminmerAbundanceFile);
					_parent._nbRescuedKminmers += 1;
				}
				
			}

			//getchar();

		}
		

	};

	

	struct MdbgEdgeMap10{
		
		public:

		boophf_t* _keys;
		vector<KminmerEdge33> _values;

		void clear(){
			delete _keys;
			vector<KminmerEdge33>().swap(_values);
		}
	};

	MdbgEdgeMap10 _mdbgEdges10;

	struct KminmerEdgeU{
		UnitigType _unitigIndex;
		bool _isReversed;
		bool _isPrefix;
	};

	struct UnitigEdgeMap{
		
		public:

		boophf_t* _keys;
		vector<vector<KminmerEdgeU>> _values;

		void clear(){
			delete _keys;
			vector<vector<KminmerEdgeU>>().swap(_values);
		}
	};

	UnitigEdgeMap _unitigEdgeMap;


	/*
	void rewriteMinimizerReadForMultiplex(){

		ofstream outputFile(_outputDir + "/read_data_corrected_abundance.txt");

		KminmerParserParallel parser(_outputDir + "/read_data_corrected.txt", 0, 0, false, false, _nbCores);
		parser.parseSequences(RewriteReadFunctor(*this, outputFile));

		outputFile.close();
	}


	class RewriteReadFunctor {

		public:

		CreateMdbg& _parent;
		ofstream& _outputFile;

		RewriteReadFunctor(CreateMdbg& parent, ofstream& outputFile) : _parent(parent), _outputFile(outputFile){
		}

		RewriteReadFunctor(const RewriteReadFunctor& copy) : _parent(copy._parent), _outputFile(copy._outputFile){
			
		}

		~RewriteReadFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {



			u_int64_t readIndex = kminmerList._readIndex;
			const vector<MinimizerType>& readMinimizers = kminmerList._readMinimizers;
			//const vector<KmerVec>& kminmers = kminmerList._kminmers;
			//const vector<ReadKminmerComplete>& kminmersInfos = kminmerList._kminmersInfo;

			vector<u_int32_t> prevAbundances = getPrevAbundances(readMinimizers);

			writeRead(readMinimizers, prevAbundances);
		}

		vector<u_int32_t> getPrevAbundances(const vector<MinimizerType>& readMinimizers){

			vector<u_int32_t> prevAbundances;

			vector<u_int64_t> rlePositions;
			vector<u_int32_t> minimizers_pos;//(minimizers.size());
			vector<KmerVec> kminmers; 
			vector<ReadKminmer> kminmersInfo;
			MDBG::getKminmers(0, _parent._kminmerSizeFirst+1, readMinimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, 0, false);
			
			for(const KmerVec& vec : kminmers){

				const u_int128_t vecHash = vec.hash128();

				if(_parent._kminmerAbundances.find(vecHash) != _parent._kminmerAbundances.end()){

					prevAbundances.push_back(_parent._kminmerAbundances[vecHash]);
					
				}
				else{
					prevAbundances.push_back(1);
				}
			}

			return prevAbundances;
		}
		
		void writeRead(const vector<MinimizerType>& minimizers, const vector<u_int32_t>& abundances){

			#pragma omp critical(writeRead)
			{

				u_int32_t size = minimizers.size();
				_outputFile.write((const char*)&size, sizeof(size));

				u_int8_t isCircular = CONTIG_LINEAR;
				_outputFile.write((const char*)&isCircular, sizeof(isCircular));

				_outputFile.write((const char*)&minimizers[0], size*sizeof(MinimizerType));
				
				u_int32_t abundanceSize = abundances.size();
				_outputFile.write((const char*)&abundanceSize, sizeof(abundanceSize));
				_outputFile.write((const char*)&abundances[0], abundanceSize*sizeof(u_int32_t));
			}

		}

	};
	*/

};	


#endif 

