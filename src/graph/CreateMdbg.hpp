
#ifndef MDBG_METAG_CREATEMDBG
#define MDBG_METAG_CREATEMDBG

#include "Commons.hpp"

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

typedef phmap::parallel_flat_hash_map<u_int128_t, u_int32_t> KminmerAbundanceMap;
typedef phmap::parallel_flat_hash_map<u_int128_t, DbgNodeLight, phmap::priv::hash_default_hash<u_int128_t>, phmap::priv::hash_default_eq<u_int128_t>, std::allocator<std::pair<u_int128_t, DbgNodeLight>>, 4, std::mutex> MdbgNodeMapLight;

typedef phmap::parallel_flat_hash_map<KmerVec, u_int32_t, phmap::priv::hash_default_hash<KmerVec>, phmap::priv::hash_default_eq<KmerVec>, std::allocator<std::pair<KmerVec, u_int32_t>>, 4, std::mutex> UnitigEdgeMap;

struct UnitigEdge{
	UnitigType _unitigIndexFrom;
	UnitigType _unitigIndexTo;
};

class CreateMdbg : public Tool{
    
	struct KminmerData{
		u_int32_t _count;
		u_int16_t _length;
		u_int16_t _overlapLength_start;
		u_int16_t _overlapLength_end;
		bool _isReversed;
	};

	struct KminmerAbundance{
		KmerVec _vec;
		u_int32_t _count;
	};

public:

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
	int _minAbundance;
	int _nbCores;
	bool _useBloomFilter;
	bool _useCorrectedRead;
	u_int32_t _unitigIndex;
	//unordered_set<u_int128_t> _isKminmerIndexed;
	UnitigEdgeMap _startNode_to_unitigIndex;
	UnitigEdgeMap _endNode_to_unitigIndex;
	//unordered_map<KmerVec, u_int32_t> _startNode_to_unitigIndex;
	//unordered_map<KmerVec, u_int32_t> _endNode_to_unitigIndex;

	BloomCacheCoherent<u_int64_t>* _bloomFilter;
	//unordered_set<KmerVec> _bloomFilterExact;
    //IBank* _inputBank;

	//vector<ReadData> _readData;
	//CompositionManager* _compositionManager;
	//AdjGraph* _overlapGraph;
	//string _inputDir;
	//string _input_extractKminmers;
	KmerVec _lala;

	u_int64_t _nbReads;
	bool _isFirstPass;
	string _filename_readMinimizers;
	string _filename_contigMinimizers;
	string _filename_inputContigs;
	string _filename_solidKminmers;

	string _filename_smallContigs;
	//string _filename_smallContigsNew;
	//string _filename_smallContigsPrev;
	//string _filename_smallContigs;

	//unordered_set<u_int32_t> writtenNodeNames;
	ofstream _kminmerFile;
	//ofstream _kminmerFile2;
	ofstream _readFile;
	ofstream _fileSmallContigs;
	//string _filename_filteredMinimizers;
	//string _filename_readCompositions;

	//vector<u_int32_t> _evaluation_readToDataset;
	bool _parseReads;
	MdbgNodeMapLight _mdbgNodesLight;
	MdbgNodeMapLight _mdbgNodesLightUnique;
	unordered_map<KmerVec, u_int32_t> _kmervec_to_nodeName;
	//MDBG* _mdbg;
	//MDBG* _mdbgFirst;
	//MDBG* _mdbgSaved;
	//MDBG* _mdbgInit;
	KminmerAbundanceMap _kminmerAbundances;
	//MdbgEdgeMap _mdbgEdges;
	MdbgEdgeMap2 _mdbgEdges2;
	MdbgEdgeMap3 _mdbgEdges3;
	MdbgEdgeMap4 _mdbgEdges4;
	//MDBG* _mdbgNoFilter;
	bool _parsingContigs;
	bool _savingSmallContigs;
	u_int32_t _node_id;
	//MinimizerPairMap* _minimizerPairMap;

	float _minimizerSpacingMean;
	float _kminmerLengthMean;
	float _kminmerOverlapMean;

	double _minimizerSpacing_sum;
	double _minimizerSpacing_n;
	double _kminmerLength_sum;
	double _kminmerLength_n;
	u_int64_t _nbEdges;

	int _lalaLol;
	vector<MinimizerType> _lalaLol2;

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
	void computeDeterministicNodeNames();
	void indexEdgeUnitig(const KmerVec& vec, u_int32_t unitigIndex);

	void createMDBG_index(const vector<MinimizerType>& minimizers, const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex);
	void createGfa();
	void parseArgs(int argc, char* argv[]);
	void parseContigs();
	float computeKmerVecAbundance(const vector<MinimizerType>& minimizers, bool isContig);
	void computeContigAbundance();
	void indexEdge(const KmerVec& vec, u_int32_t nodeName);
	void computeEdgeOld(const KmerVec& vec, u_int32_t id);
	void dumpEdgeOld(u_int32_t nodeNameFrom, u_int8_t nodeNameFromOri, u_int32_t nodeNameTo, u_int8_t nodeNameToOri);
	void indexEdges();
	void indexUnitigEdges();
	void indexUnitigEdge(const u_int32_t unitigIndex, const vector<MinimizerType>& minimizers);
	void computeUnitigEdges();
	void computeUnitigEdge(const u_int32_t unitigIndex, const vector<MinimizerType>& minimizers);
	void computeUnitigNodes();
	//void computeUnitigNode(const KmerVec& vec);
	void getSuccessors(const KmerVec& vec, vector<KmerVec>& successors);
	void getPredecessors(const KmerVec& vec, vector<KmerVec>& predecessors);
	bool hasSingleSuccessor(const KmerVec& vec);
	bool hasSinglePredecessor(const KmerVec& vec);
	void getSuccessors_unitig(const KmerVec& vec, UnitigType unitigIndexFrom, vector<UnitigEdge>& successors);
	void getPredecessors_unitig(const KmerVec& vec, UnitigType unitigIndexFrom, vector<UnitigEdge>& predecessors);
	void computeEdgesOld();
	void dumpUnitigNode(const UnitigType& unitigIndex, const vector<MinimizerType>& unitig);
	void dumpUnitigEdge(u_int32_t fromUnitigIndex, u_int32_t toUnitigIndex, bool isSuccessor);

	void computeUnitiStartEndNode();
	//bool isBranchingNode(const KmerVec& vec);
	//void setNodeUnitigged(const KmerVec& vec);
	//bool isNodeUnitigged(const KmerVec& vec);
	void dumpUnitigAbundances();
	void loadRefinedAbundances();

	
	u_int64_t _edgeIndexElements;
	ofstream _unitigGraphFile_nodes;
	ofstream _unitigGraphFile_nodes_abundances;
	ofstream _unitigGraphFile_edges_successors;
	ofstream _unitigGraphFile_edges_predecessors;
	//ofstream _kminmerFile;


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
	unordered_map<KmerVec, u_int32_t> _solidKminmerAbundances;
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
	
	class IndexKminmerFunctor {

		public:

		CreateMdbg& _graph;
		bool _isFirstPass;
		bool _parsingContigs;
		ofstream& _readFile;
		//MDBG* _mdbg;
		KminmerAbundanceMap& _kminmerAbundances;
		ofstream& _kminmerFile;
		float _minimizerSpacingMean;
		float _kminmerLengthMean;
		float _kminmerOverlapMean;
		size_t _minimizerSize;
		size_t _kminmerSize;
		size_t _kminmerSizeFirst;
		//BloomCacheCoherent<u_int64_t>* _bloomFilter;
		bool _extractingContigs;

		IndexKminmerFunctor(CreateMdbg& graph, bool extractingContigs) : _graph(graph), _readFile(graph._readFile), _kminmerAbundances(graph._kminmerAbundances), _kminmerFile(graph._kminmerFile){
			_isFirstPass = graph._isFirstPass;
			_parsingContigs = graph._parsingContigs;
			//_mdbg = graph._mdbg;
			//_mdbgInit = graph._mdbgInit;
			_minimizerSpacingMean = graph._minimizerSpacingMean;
			_kminmerLengthMean = graph._kminmerLengthMean;
			_kminmerOverlapMean = graph._kminmerOverlapMean;
			_minimizerSize = graph._minimizerSize;
			_kminmerSize = graph._kminmerSize;
			_kminmerSizeFirst = graph._kminmerSizeFirst;
			//_bloomFilter = graph._bloomFilter;
			_extractingContigs = extractingContigs;
		}

		IndexKminmerFunctor(const IndexKminmerFunctor& copy) : _graph(copy._graph), _readFile(copy._readFile), _kminmerAbundances(copy._kminmerAbundances), _kminmerFile(copy._kminmerFile){
			_isFirstPass = copy._isFirstPass;
			_parsingContigs = copy._parsingContigs;
			//_mdbg = copy._mdbg;
			//_mdbgInit = copy._mdbgInit;
			_minimizerSpacingMean = copy._minimizerSpacingMean;
			_kminmerLengthMean = copy._kminmerLengthMean;
			_kminmerOverlapMean = copy._kminmerOverlapMean;
			_minimizerSize = copy._minimizerSize;
			_kminmerSize = copy._kminmerSize;
			_kminmerSizeFirst = copy._kminmerSizeFirst;
			//_bloomFilter = copy._bloomFilter;
			_extractingContigs = copy._extractingContigs;
		}

		~IndexKminmerFunctor(){
			//delete _minimizerParser;
		}

		/*
		bool isValid(const vector<u_int64_t>& readMinimizers){
			
			//cout << "is valid: " << readMinimizers.size() << endl;

			for(size_t kminmerSize=_kminmerSizeFirst; kminmerSize<_kminmerSize; kminmerSize+=1){

				if(kminmerSize >= _graph._bannedKminmers.size()) return true; //pass k=4 and k =5

				const phmap::parallel_flat_hash_set<KmerVec>& bannedKminmers = _graph._bannedKminmers[kminmerSize];

				//cout << kminmerSize << " " << bannedKminmers.size() << endl;
				if(bannedKminmers.size() == 0) continue;
				//vector<phmap::parallel_flat_hash_set<KmerVec>> _bannedKminmers;
				
				vector<u_int64_t> rlePositions;
				vector<u_int64_t> minimizers_pos;//(minimizers.size());
				vector<KmerVec> kminmers; 
				vector<ReadKminmer> kminmersInfo;
				MDBG::getKminmers(_minimizerSize, kminmerSize, readMinimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, 0, false);

				for(const KmerVec& vec : kminmers){
					if(bannedKminmers.find(vec) != bannedKminmers.end()){
						//cout << "no" << endl;
						return false;
					}
				}
			}

			//cout << "yes" << endl;
			return true;
		}
		*/

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
		
		float getAbundance(const vector<MinimizerType>& readMinimizers){

			vector<u_int64_t> rlePositions;
			vector<u_int32_t> minimizers_pos;//(minimizers.size());
			vector<KmerVec> kminmers; 
			vector<ReadKminmer> kminmersInfo;
			MDBG::getKminmers(_minimizerSize, _graph._kminmerSizePrev, readMinimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, 0, false);

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

				u_int128_t hash = vec.hash128();

				//cout << (_mdbgInit->_dbg_nodes.find(vec) != _mdbgInit->_dbg_nodes.end()) << endl;
				if(_kminmerAbundances.find(hash) != _kminmerAbundances.end()){

					//u_int32_t abundance = _kminmerAbundances[vec];
					//const KminmerAbundanceMapNode& node = _kminmerAbundances[vec];
					u_int32_t abundance = _kminmerAbundances[hash];
					//if(nodeName == 17430) cout << abundance << endl;

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
		
		bool isInGraph(const vector<ReadKminmerComplete>& kminmersInfos){

			if(kminmersInfos.size() == 0) return false;

			for(size_t i=0; i<kminmersInfos.size(); i++){
					

				//const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];

				const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];
				const KmerVec& vec = kminmerInfo._vec;

				if(_graph._mdbgNodesLight.find(vec.hash128()) == _graph._mdbgNodesLight.end()){
				//if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()){
					return false;
				}
			}

			return true;
		}

		void operator () (const KminmerList& kminmerList) {


			u_int64_t readIndex = kminmerList._readIndex;
			const vector<MinimizerType>& readMinimizers = kminmerList._readMinimizers;
			//const vector<KmerVec>& kminmers = kminmerList._kminmers;
			const vector<ReadKminmerComplete>& kminmersInfos = kminmerList._kminmersInfo;

			if(_graph._savingSmallContigs){

				if(_extractingContigs && _kminmerSize > 8 && kminmersInfos.size() <= 0 && getAbundance(readMinimizers) > 1){

					if(!isInGraph(kminmersInfos)){
						#pragma omp critical
						{
							u_int32_t contigSize = readMinimizers.size();
							_graph._fileSmallContigs.write((const char*)&contigSize, sizeof(contigSize));

							u_int8_t isCircular = kminmerList._isCircular;
							_graph._fileSmallContigs.write((const char*)&isCircular, sizeof(isCircular));
							_graph._fileSmallContigs.write((const char*)&readMinimizers[0], contigSize*sizeof(MinimizerType));
							//cout << "small contig" << endl;
							//getchar();

							if(kminmersInfos.size() > 0){
								_graph._nbSmallContigs += 1;
							}
						}
					
					}
					
				}


				return;
			}


			/*
			if(_extractingContigs && _kminmerSize > 8 && kminmersInfos.size() < _kminmerSize*2 && getAbundance(readMinimizers) > 1){

				bool isInGraph = true;
				
				for(size_t i=0; i<kminmersInfos.size(); i++){
					

					//const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];

					const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];
					const KmerVec& vec = kminmerInfo._vec;

					if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()){
						isInGraph = false;
						break;
					}
				}

				//if(kminmersInfos.size() == 0){
					//cout << isInGraph << endl;
				//}

				if(!isInGraph){
					#pragma omp critical
					{
						u_int32_t contigSize = readMinimizers.size();
						_graph._fileSmallContigs.write((const char*)&contigSize, sizeof(contigSize));

						u_int8_t isCircular = kminmerList._isCircular;
						_graph._fileSmallContigs.write((const char*)&isCircular, sizeof(isCircular));
						_graph._fileSmallContigs.write((const char*)&readMinimizers[0], contigSize*sizeof(u_int64_t));
						//cout << "small contig" << endl;
						//getchar();
						_graph._nbSmallContigs += 1;
					}
				
				}

			}

			*/

			/*
			if(_extractingContigs && _kminmerSize > 8 && kminmersInfos.size() == 0 && getAbundance(readMinimizers) > 1){

				#pragma omp critical
				{
					u_int32_t contigSize = readMinimizers.size();
					_graph._fileSmallContigs.write((const char*)&contigSize, sizeof(contigSize));

					u_int8_t isCircular = kminmerList._isCircular;
					_graph._fileSmallContigs.write((const char*)&isCircular, sizeof(isCircular));
					_graph._fileSmallContigs.write((const char*)&readMinimizers[0], contigSize*sizeof(u_int64_t));
					//cout << "small contig" << endl;
					//getchar();
				}
				
				return;
			}
			*/

			
			if(_extractingContigs && kminmersInfos.size() < _kminmerSize*2 ){ //&& getAbundance(readMinimizers) > 1
				//if(isInGraph(kminmersInfos)) return;
				//return;
			}

			
			

			/*
			//if(!isValid(readMinimizers)){
			if(_extractingContigs && _graph.isUnitigDuplicated.find(readIndex) != _graph.isUnitigDuplicated.end()){	
				#pragma omp critical
				{
					//cout << "Invalid: " << readMinimizers.size() << endl;
					u_int32_t contigSize = readMinimizers.size();
					_graph._fileSmallContigs.write((const char*)&contigSize, sizeof(contigSize));

					u_int8_t isCircular = kminmerList._isCircular;
					_graph._fileSmallContigs.write((const char*)&isCircular, sizeof(isCircular));
					_graph._fileSmallContigs.write((const char*)&readMinimizers[0], contigSize*sizeof(u_int64_t));
				}

				return;
			}
			*/
			
			if(_extractingContigs){
				//u_int64_t length = _kminmerLengthMean + ((kminmersInfos.size()-1) * (_kminmerLengthMean-_kminmerOverlapMean));
				//if(length < 50000) return;
				//#pragma omp critical
				//{
				//	cout << kminmersInfos.size() << " " << kminmerList._isCircular << endl;

					//if(kminmersInfos.size() < 200) return;
				//}

			}

			//if(_extractingContigs && kminmersInfos.size() < 200) return;
			/*
			double abundanceCutoff = 0;
			if(!_isFirstPass){
				
				vector<u_int64_t> rlePositions;
				vector<u_int64_t> minimizers_pos;//(minimizers.size());
				vector<KmerVec> kminmers; 
				vector<ReadKminmer> kminmersInfo;
				MDBG::getKminmers(_minimizerSize, _kminmerSize-1, readMinimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, 0, false);


				//u_int32_t minAbundance = -1;
				double sumTotal = 0;
				double nTotal = 0;
				//vector<float> abundances;
				//cout << kminmers.size() << endl;
				
				float lastAbundance = -1;
				float lastN = -1;
				double nMax = 0;

				//cout << "-------------" << endl;

				for(const KmerVec& vec : kminmers){
					//cout << (_mdbgInit->_dbg_nodes.find(vec) != _mdbgInit->_dbg_nodes.end()) << endl;


					if(_mdbgInit->_dbg_nodes.find(vec) != _mdbgInit->_dbg_nodes.end()){

						u_int32_t abundance = _mdbgInit->_dbg_nodes[vec]._abundance;
						size_t n = _mdbgInit->_dbg_nodes[vec]._unitigNbNodes;


						//if(abundance == lastAbundance && lastN == n) continue;

						sumTotal += abundance * n;
						nTotal += n;
						
						if (n > nMax) nMax = n;
						//for(size_t i=0; i<n; i++){
							//sumN._n += 1;
							//sumN._sum += abundance;
						//}

						//sumN._nbNodes_sum += n;
						//sumN._nbNodes_n += 1;

						//cout << abundance << " " << n << endl;

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

				//cout << nTotal << endl;
				if(nTotal == 0){
					abundanceCutoff = 0;
				}
				else{
					//cout << sumTotal << " " << nTotal << endl;
					double abundanceMean = sumTotal / nTotal;
					abundanceCutoff = abundanceMean * 0.5;

					if(nMax < 300) abundanceCutoff = 0;
					//cout << abundanceMean << " " << nMax << endl;

					//if(abundanceMean > 60) getchar();
				}
				//abundanceCutoff = ;

			}
			*/
			
			/*
			bool isHere = false;

			vector<u_int64_t> rlePositions;
			vector<u_int64_t> minimizers_pos;//(minimizers.size());
			vector<KmerVec> kminmers; 
			vector<ReadKminmer> kminmersInfo;
			MDBG::getKminmers(_minimizerSize, _kminmerSizeFirst, readMinimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, readIndex, false);


			for(const KmerVec& vec : kminmers){
				if(_mdbgInit->_dbg_nodes.find(vec) != _mdbgInit->_dbg_nodes.end()){

					if(vec._kmers[0] == 66918863945726617 && vec._kmers[1] == 46622693399843280 && vec._kmers[2] == 91239340561015544){
						isHere = true;
					}
				}
			}

			if(isHere){
				cout << "yes! " << readIndex << " " << _extractingContigs << endl;
				getchar();
			}
			*/


			/*
			if(_extractingContigs){
				if(kminmersInfos.size() == 0){
					#pragma omp critical
					{
						u_int32_t contigSize = readMinimizers.size();
						//_graph._fileSmallContigs.write((const char*)&contigSize, sizeof(contigSize));
						//_graph._fileSmallContigs.write((const char*)&readMinimizers[0], contigSize*sizeof(u_int64_t));
						//cout << "small contig" << endl;
						//getchar();
					}
				}
			}
			*/
			/*
			//Ici save too short contigs
			#pragma omp critical
			{
				if(_kminmerSize > 30){
					//if(_extractingContigs){
						cout << readMinimizers.size() << " " << kminmersInfos.size() << endl;
					//}
				}
			}
			*/



			//cout << readIndex << endl;
			//vector<u_int16_t> minimizerPosOffset(readMinimizers.size());
			//for(size_t i=0; i<minimizerPosOffset.size(); i++){
			//	minimizerPosOffset[i] = _minimizerSpacingMean;
			//}

				//if(!_isFirstPass && !_parsingContigs){
				//	u_int32_t size = readMinimizers.size();
				//	_readFile.write((const char*)&size, sizeof(size));
				//	_readFile.write((const char*)&readMinimizers[0], size*sizeof(u_int64_t));
					//_readFile.write((const char*)&minimizerPosOffset[0], size*sizeof(u_int16_t));
				//}

				//cout << readIndex << " " << kminmersInfos.size() << endl;
				for(size_t i=0; i<kminmersInfos.size(); i++){
					

					//const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];

					const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];
					const KmerVec& vec = kminmerInfo._vec;
					//if(vec.isPalindrome()) continue;

					float kminmerAbundance = -1;
					if(!_isFirstPass){
						kminmerAbundance = getAbundance(readMinimizers, kminmerInfo);
						if(kminmerAbundance == 1){//} && ab < abundanceCutoff){
							//cout << "Out: " << ab << " " << abundanceCutoff << endl;
							continue;
						}
					}


					bool exist = true;

					if(_graph._useBloomFilter){

						exist = false;


						#pragma omp critical
						{
							/*
							if(_graph._bloomFilterExact.find(vec) != _graph._bloomFilterExact.end()){
								exist = true;
							}
							else{
								_graph._bloomFilterExact.insert(vec);
							}
							*/

							
							if(_graph._bloomFilter->contains(vec.h())){
								exist = true;
							}
							else{
								exist = false;
								_graph._bloomFilter->insert(vec.h());
							}
							
						}
					}


					if(exist){
						
						bool isNewKey = false;



						_graph._mdbgNodesLight.lazy_emplace_l(vec.hash128(), 
						[this, &kminmerInfo](MdbgNodeMapLight::value_type& v) { // key exist
							//v.second._quality += kminmerInfo._quality;
							if(_isFirstPass){
								v.second._abundance += 1;
							}
						},           
						[&vec, this, &kminmerInfo, &readMinimizers, &readIndex, &isNewKey, &kminmerAbundance](const MdbgNodeMapLight::constructor& ctor) { // key inserted
							
							
							isNewKey = true;


							u_int32_t nodeName = -1;
							DbgNodeLight node = {false, nodeName, 2};

							if(_isFirstPass){
								if(_graph._useBloomFilter){
									node._abundance = 2;
								}
								else{
									node._abundance = 1;
								}
							}
							else{


								node._abundance = kminmerAbundance; 
							}

							//node._quality = 0;
							//node._quality += kminmerInfo._quality;

							ctor(vec.hash128(), node); 

							//cout << "jinsere" << endl;
						}); // construct value_type in place when key not present

						
						if(isNewKey){

							u_int32_t nodeName;

							#pragma omp critical
							{

								MDBG::writeKminmer(vec._kmers, _graph._kminmerFile);

								
								nodeName = _graph._node_id;
								_graph._node_id += 1;
								

							}

							auto set_value = [&nodeName](MdbgNodeMapLight::value_type& v) { v.second._index = nodeName; };
							_graph._mdbgNodesLight.modify_if(vec.hash128(), set_value);
						}


					}


			}


		}
	};
	


	class FilterKminmerFunctor {

		public:

		CreateMdbg& _graph;
		bool _isFirstPass;
		bool _parsingContigs;
		ofstream& _readFile;
		//MDBG* _mdbg;
		KminmerAbundanceMap& _kminmerAbundances;
		ofstream& _kminmerFile;
		float _minimizerSpacingMean;
		float _kminmerLengthMean;
		float _kminmerOverlapMean;
		size_t _minimizerSize;
		size_t _kminmerSize;
		size_t _kminmerSizeFirst;
		//BloomCacheCoherent<u_int64_t>* _bloomFilter;

		FilterKminmerFunctor(CreateMdbg& graph) : _graph(graph), _readFile(graph._readFile), _kminmerAbundances(graph._kminmerAbundances), _kminmerFile(graph._kminmerFile){
			_isFirstPass = graph._isFirstPass;
			_parsingContigs = graph._parsingContigs;
			//_mdbg = graph._mdbg;
			//_mdbgInit = graph._mdbgInit;
			_minimizerSpacingMean = graph._minimizerSpacingMean;
			_kminmerLengthMean = graph._kminmerLengthMean;
			_kminmerOverlapMean = graph._kminmerOverlapMean;
			_minimizerSize = graph._minimizerSize;
			_kminmerSize = graph._kminmerSize;
			_kminmerSizeFirst = graph._kminmerSizeFirst;
			//_bloomFilter = graph._bloomFilter;
		}

		FilterKminmerFunctor(const FilterKminmerFunctor& copy) : _graph(copy._graph), _readFile(copy._readFile), _kminmerAbundances(copy._kminmerAbundances), _kminmerFile(copy._kminmerFile){
			_isFirstPass = copy._isFirstPass;
			_parsingContigs = copy._parsingContigs;
			//_mdbg = copy._mdbg;
			//_mdbgInit = copy._mdbgInit;
			_minimizerSpacingMean = copy._minimizerSpacingMean;
			_kminmerLengthMean = copy._kminmerLengthMean;
			_kminmerOverlapMean = copy._kminmerOverlapMean;
			_minimizerSize = copy._minimizerSize;
			_kminmerSize = copy._kminmerSize;
			_kminmerSizeFirst = copy._kminmerSizeFirst;
			//_bloomFilter = copy._bloomFilter;
		}

		~FilterKminmerFunctor(){
			//delete _minimizerParser;
		}

	

		void operator () (const KminmerList& kminmerList) {


			u_int64_t readIndex = kminmerList._readIndex;
			const vector<MinimizerType>& readMinimizers = kminmerList._readMinimizers;
			//const vector<KmerVec>& kminmers = kminmerList._kminmers;
			const vector<ReadKminmerComplete>& kminmersInfos = kminmerList._kminmersInfo;


			vector<u_int32_t> abundances;
			//vector<u_int32_t> qualities;
			for(size_t i=0; i<kminmersInfos.size(); i++){
				
				const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];
				const KmerVec& vec = kminmerInfo._vec;
				u_int128_t hash = vec.hash128();
				
				//if(_kminmerAbundances.find(vec) == _kminmerAbundances.end()){
				if(_graph._mdbgNodesLight.find(hash) == _graph._mdbgNodesLight.end()){
					abundances.push_back(1);
				}
				else{
					abundances.push_back(_graph._mdbgNodesLight[hash]._abundance);
				}

				//qualities.push_back(kminmerInfo._quality);
			}

			u_int32_t median = Utils::compute_median(abundances);
			double cutoff = median * 0.1f;

			/*
			#pragma omp critical
			{
				cout << "----- " << readIndex << endl;
				for(u_int32_t ab : abundances){
					cout << ab << " ";
				}
				cout << endl;
				for(u_int8_t qual : qualities){
					cout << to_string(qual) << " ";
				}
				cout << endl;
				cout << cutoff << endl;
				getchar();
			}
			*/
			

			
			if(cutoff > 1) return;

			//cout << readIndex << " " << kminmersInfos.size() << endl;
			for(size_t i=0; i<kminmersInfos.size(); i++){
				
				const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];
				const KmerVec& vec = kminmerInfo._vec;
				u_int128_t hash = vec.hash128();
				//if(kminmerInfo._quality < 5) continue;
				//cout << (int) kminmerInfo._quality << endl;

				if(_graph._mdbgNodesLight.find(hash) != _graph._mdbgNodesLight.end()) continue;

				
				bool isNewKey = false;

				_graph._mdbgNodesLightUnique.lazy_emplace_l(hash, 
				[this](MdbgNodeMapLight::value_type& v) { // key exist
					//if(_isFirstPass){
					//	v.second._abundance += 1;
					//}
				},           
				[&hash, this, &kminmerInfo, &readMinimizers, &readIndex, &isNewKey](const MdbgNodeMapLight::constructor& ctor) { // key inserted
					
					
					isNewKey = true;

					u_int32_t nodeName = -1;
					DbgNodeLight node = {false, nodeName, 2};

					node._abundance = 1;
					//node._quality = kminmerInfo._quality;
					//node._quality += kminmerInfo._quality;

					ctor(hash, node); 

					//cout << "jinsere" << endl;
				}); // construct value_type in place when key not present

				if(isNewKey){

					u_int32_t nodeName;

					#pragma omp critical
					{
						
						MDBG::writeKminmer(vec._kmers, _graph._kminmerFile);

						nodeName = _graph._node_id;
						_graph._node_id += 1;

						//cout << "Indexed: " << nodeName << endl;
					}

					auto set_value = [&nodeName](MdbgNodeMapLight::value_type& v) { v.second._index = nodeName; };
					_graph._mdbgNodesLightUnique.modify_if(hash, set_value);
				}
				
			}

			//getchar();

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

			//KmerVec lol;
			//lol._kmers = {77646226911515249, 90273684072333150, 83465560366503996};

			bool isReversedPrefix;
			bool isReversedSuffix;
			KmerVec prefix = vec.prefix().normalize(isReversedPrefix);
			KmerVec suffix = vec.suffix().normalize(isReversedSuffix);
			MinimizerType prefix_minimizer = vec._kmers[vec._kmers.size()-1];
			MinimizerType suffix_minimizer = vec._kmers[0];
			
			for(KminmerEdge3& edge : _graph._mdbgEdges3[suffix.hash128()]){
				if(edge._isReversed != isReversedSuffix) continue;
				if(edge._isPrefix) continue;

				if(edge._minimizer == suffix_minimizer){
					//if(suffix == lol){
					//	cout << "omg suffix " << suffix.toString() << " " << isReversedSuffix << " " << suffix_minimizer << endl;
					//	getchar();
					//}
					edge._isUnitigged = true;
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
					break;
				}
			}
		}

		bool isNodeUnitigged(const KmerVec& vec){

			bool isReversedPrefix;
			bool isReversedSuffix;
			KmerVec prefix = vec.prefix().normalize(isReversedPrefix);
			KmerVec suffix = vec.suffix().normalize(isReversedSuffix);
			MinimizerType prefix_minimizer = vec._kmers[vec._kmers.size()-1];
			MinimizerType suffix_minimizer = vec._kmers[0];
			
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
			KmerVec v = suffix;
			if(edge._isReversed){
				v = v.reverse();
			}
			if(edge._isPrefix){
				v._kmers.push_back(edge._minimizer);
			}
			else{
				v._kmers.insert(v._kmers.begin(), edge._minimizer);
			}
			*/


			//Good
			//63003427831326793 54943957852837162 37090168052080844 45142317365530159


			//Bad
			//37090168052080844 54943957852837162 63003427831326793 45142317365530159

			//KmerVec lala;
			//lala._kmers = {37090168052080844, 54943957852837162, 63003427831326793};

			//if(prefix == lala){
				//cout << "omg: " << prefix_minimizer << endl;
				//getchar();
			//}

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
			

			const KmerVec& sourceNodeRC = sourceNode.reverse();

			if(!isNodeExists(sourceNode)) return;
		
		
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



			if(exist) return;

			vector<MinimizerType> unitig;
			KmerVec startNode = sourceNode;
			KmerVec endNode = sourceNode;

			u_int8_t isCircular = CONTIG_LINEAR;
			
			unitig = sourceNode._kmers;


			//forward
			while(true){


				if(!_graph.hasSingleSuccessor(endNode)) break;

				vector<KmerVec> successors;
				_graph.getSuccessors(endNode, successors);

				const KmerVec& successor = successors[0];

				if(!_graph.hasSinglePredecessor(successor)) break;


				if(successor == sourceNode){ //Circular
					endNode = successor;
					startNode = successor;
					isCircular = CONTIG_CIRCULAR;
					break; 
				}


				endNode = successor;

				unitig.push_back(endNode._kmers[endNode._kmers.size()-1]);
			}


			if(!isCircular){

				vector<MinimizerType> unitigBackward;

				startNode = sourceNode;

				//backward
				while(true){
					

					if(!_graph.hasSinglePredecessor(startNode)) break;

					vector<KmerVec> predecessors;
					_graph.getPredecessors(startNode, predecessors);
					
					KmerVec predecessor = predecessors[0];

					if(!_graph.hasSingleSuccessor(predecessor)) break;
					
					vector<KmerVec> successors;
					_graph.getSuccessors(predecessor, successors);

					startNode = predecessor;

					unitigBackward.push_back(startNode._kmers[0]);

				}

				std::reverse(unitigBackward.begin(), unitigBackward.end());

				unitig.insert(unitig.begin(), unitigBackward.begin(), unitigBackward.end());

				
			}
			

			KmerVec startNodeRC = endNode.reverse();
			KmerVec endNodeRC = startNode.reverse();

			
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



					_graph.dumpUnitigNode(_graph._unitigIndex, unitig);

					_graph._unitigIndex += 2;
					
				}

			}
			
			
		}

	};

};	


#endif 

