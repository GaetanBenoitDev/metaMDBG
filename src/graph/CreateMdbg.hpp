
#ifndef MDBG_METAG_CREATEMDBG
#define MDBG_METAG_CREATEMDBG

#include "Commons.hpp"

//#include "graph/Graph.hpp"
#include "GfaParser.hpp"
#include "GraphSimplify.hpp"
#include "../utils/BloomFilter.hpp"
//#include "./utils/ntHashIterator.hpp"

//#include <boost/graph/adjacency_list.hpp>
//#include <boost/graph/graph_utility.hpp>
//#include <boost/property_map/property_map.hpp>
//#include <boost/graph/connected_components.hpp>
//#include <boost/graph/dijkstra_shortest_paths.hpp>




struct KminmerAbundanceMapNode{
	u_int32_t _abundance;
	u_int32_t _contigIndex;
};

typedef phmap::parallel_flat_hash_map<KmerVec, KminmerAbundanceMapNode> KminmerAbundanceMap;




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
	string _filename_noKminmerReads;
	//ofstream _file_noKminmerReads;
	int _nbCores;
	bool _useBloomFilter;

	BloomCacheCoherent<u_int64_t>* _bloomFilter;
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

	unordered_set<u_int32_t> writtenNodeNames;
	unordered_set<KmerVec> _kminmerExist;
	ofstream _kminmerFile;
	//ofstream _kminmerFile2;
	ofstream _readFile;
	ofstream _fileSmallContigs;
	//string _filename_filteredMinimizers;
	//string _filename_readCompositions;

	//vector<u_int32_t> _evaluation_readToDataset;
	bool _parseReads;
	MDBG* _mdbg;
	MDBG* _mdbgSaved;
	//MDBG* _mdbgInit;
	KminmerAbundanceMap _kminmerAbundances;
	//MdbgEdgeMap _mdbgEdges;
	MdbgEdgeMap2 _mdbgEdges2;
	MDBG* _mdbgNoFilter;
	bool _parsingContigs;
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
	ofstream _outputFileGfa;

    CreateMdbg ();
    void execute ();
	void createMDBG();
	//void createMDBG_read(kseq_t* read, u_int64_t readIndex);
	void createMDBG_collectKminmers(const vector<u_int64_t>& minimizers, const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex);
	void createMDBG_collectKminmers_contig(const vector<u_int64_t>& minimizers, const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex);
	void removeErroneousKminmers(const vector<u_int64_t>& minimizers, const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex);
	void loadSolidKminmers();

	void createMDBG_collectKminmers_read(kseq_t* read, u_int64_t readIndex);
	void extractKminmerSequence(const char* sequenceOriginal, const ReadKminmer& kminmerInfo, string& sequence);
	void computeDeterministicNodeNames();

	void createMDBG_index(const vector<u_int64_t>& minimizers, const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex);
	void createGfa();
	void parseArgs(int argc, char* argv[]);
	void parseContigs();
	float computeKmerVecAbundance(const vector<u_int64_t>& minimizers, bool isContig);
	void computeContigAbundance();
	void indexEdge(const KmerVec& vec, u_int32_t nodeName);
	void computeEdge(const KmerVec& vec, u_int32_t id);
	void dumpEdge(u_int32_t nodeNameFrom, u_int8_t nodeNameFromOri, u_int32_t nodeNameTo, u_int8_t nodeNameToOri);
	void indexEdges();
	void computeEdges();
	//void createMDBG_collectKminmers_minspace_read(const vector<u_int64_t>& readMinimizers, const vector<ReadKminmerComplete>& kminmersInfos, u_int64_t readIndex);

	//void createSimilarityGraph(GraphInfo* graphInfo);
	//void execute_binning();
	//void execute_binning_cleanGraph();
	//void extract_kminmers();

	u_int64_t _debug_nbMinimizers;
	unordered_map<u_int64_t, u_int64_t> _minimizerCounts;
	//unordered_map<KmerVec, u_int32_t> kminmerCounts;
	unordered_map<KmerVec, KminmerData> _kminmersData;
	gzFile _file_readData;
	MinimizerParser* _minimizerParser;
	unordered_map<KmerVec, vector<u_int32_t>> _contigIndex;
	unordered_map<u_int32_t, u_int32_t> _contigAbundances;
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

	
	class IndexKminmerFunctor {

		public:

		CreateMdbg& _graph;
		bool _isFirstPass;
		bool _parsingContigs;
		ofstream& _readFile;
		unordered_set<KmerVec>& _kminmerExist;
		MDBG* _mdbg;
		KminmerAbundanceMap& _kminmerAbundances;
		ofstream& _kminmerFile;
		double _minimizerSpacingMean;
		double _kminmerLengthMean;
		double _kminmerOverlapMean;
		size_t _minimizerSize;
		size_t _kminmerSize;
		size_t _kminmerSizeFirst;
		BloomCacheCoherent<u_int64_t>* _bloomFilter;
		bool _extractingContigs;

		IndexKminmerFunctor(CreateMdbg& graph, bool extractingContigs) : _graph(graph), _readFile(graph._readFile), _kminmerExist(graph._kminmerExist), _kminmerAbundances(graph._kminmerAbundances), _kminmerFile(graph._kminmerFile){
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
			_bloomFilter = graph._bloomFilter;
			_extractingContigs = extractingContigs;
		}

		IndexKminmerFunctor(const IndexKminmerFunctor& copy) : _graph(copy._graph), _readFile(copy._readFile), _kminmerExist(copy._kminmerExist), _kminmerAbundances(copy._kminmerAbundances), _kminmerFile(copy._kminmerFile){
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
			_bloomFilter = copy._bloomFilter;
			_extractingContigs = copy._extractingContigs;
		}

		~IndexKminmerFunctor(){
			//delete _minimizerParser;
		}

		
		float getAbundance(const vector<u_int64_t>& readMinimizers, const ReadKminmerComplete& kminmerInfo){

			const KmerVec& vec = kminmerInfo._vec;

			vector<u_int64_t> minimizerSeq;
			for(size_t i=kminmerInfo._read_pos_start; i<=kminmerInfo._read_pos_end; i++){
				minimizerSeq.push_back(readMinimizers[i]);
			}
			if(kminmerInfo._isReversed){
				std::reverse(minimizerSeq.begin(), minimizerSeq.end());
			}

			return getAbundance(minimizerSeq);
		}
		
		float getAbundance(const vector<u_int64_t>& readMinimizers){

			vector<u_int64_t> rlePositions;
			vector<u_int64_t> minimizers_pos;//(minimizers.size());
			vector<KmerVec> kminmers; 
			vector<ReadKminmer> kminmersInfo;
			MDBG::getKminmers(_minimizerSize, _graph._kminmerSizePrev, readMinimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, 0, false);


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
					minAbundance = 1;
					break;
				}
			}

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

		void operator () (const KminmerList& kminmerList) {


			u_int64_t readIndex = kminmerList._readIndex;
			const vector<u_int64_t>& readMinimizers = kminmerList._readMinimizers;
			//const vector<KmerVec>& kminmers = kminmerList._kminmers;
			const vector<ReadKminmerComplete>& kminmersInfos = kminmerList._kminmersInfo;


			if(_extractingContigs && _kminmerSize > 8 && kminmersInfos.size() == 0 && getAbundance(readMinimizers) > 1){

				#pragma omp critical
				{
					u_int32_t contigSize = readMinimizers.size();
					_graph._fileSmallContigs.write((const char*)&contigSize, sizeof(contigSize));

					bool isCircular = false;
					_graph._fileSmallContigs.write((const char*)&isCircular, sizeof(isCircular));
					_graph._fileSmallContigs.write((const char*)&readMinimizers[0], contigSize*sizeof(u_int64_t));
					//cout << "small contig" << endl;
					//getchar();
				}
			}

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
							if(_bloomFilter->contains(vec.h())){
								exist = true;
							}
							else{
								exist = false;
								_bloomFilter->insert(vec.h());
							}
						}
					}


					if(exist){
						
						bool isNewKey = false;

						_mdbg->_dbg_nodes.lazy_emplace_l(vec, 
						[this, &kminmerInfo](MdbgNodeMap::value_type& v) { // key exist
							v.second._quality += kminmerInfo._quality;
							if(_isFirstPass){
								v.second._abundance += 1;
							}
						},           
						[&vec, this, &kminmerInfo, &readMinimizers, &readIndex, &isNewKey, &kminmerAbundance](const MdbgNodeMap::constructor& ctor) { // key inserted
							
							
							isNewKey = true;


							u_int32_t nodeName = -1;
							DbgNode node = {nodeName, 2};

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

							node._quality = 0;
							node._quality += kminmerInfo._quality;

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


		}
	};
	


	class FilterKminmerFunctor {

		public:

		CreateMdbg& _graph;
		bool _isFirstPass;
		bool _parsingContigs;
		ofstream& _readFile;
		unordered_set<KmerVec>& _kminmerExist;
		MDBG* _mdbg;
		KminmerAbundanceMap& _kminmerAbundances;
		ofstream& _kminmerFile;
		double _minimizerSpacingMean;
		double _kminmerLengthMean;
		double _kminmerOverlapMean;
		size_t _minimizerSize;
		size_t _kminmerSize;
		size_t _kminmerSizeFirst;
		BloomCacheCoherent<u_int64_t>* _bloomFilter;

		FilterKminmerFunctor(CreateMdbg& graph) : _graph(graph), _readFile(graph._readFile), _kminmerExist(graph._kminmerExist), _kminmerAbundances(graph._kminmerAbundances), _kminmerFile(graph._kminmerFile){
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
			_bloomFilter = graph._bloomFilter;
		}

		FilterKminmerFunctor(const FilterKminmerFunctor& copy) : _graph(copy._graph), _readFile(copy._readFile), _kminmerExist(copy._kminmerExist), _kminmerAbundances(copy._kminmerAbundances), _kminmerFile(copy._kminmerFile){
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
			_bloomFilter = copy._bloomFilter;
		}

		~FilterKminmerFunctor(){
			//delete _minimizerParser;
		}

	

		void operator () (const KminmerList& kminmerList) {


			u_int64_t readIndex = kminmerList._readIndex;
			const vector<u_int64_t>& readMinimizers = kminmerList._readMinimizers;
			//const vector<KmerVec>& kminmers = kminmerList._kminmers;
			const vector<ReadKminmerComplete>& kminmersInfos = kminmerList._kminmersInfo;


			vector<u_int32_t> abundances;
			//vector<u_int32_t> qualities;
			for(size_t i=0; i<kminmersInfos.size(); i++){
				
				const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];
				const KmerVec& vec = kminmerInfo._vec;

				
				//if(_kminmerAbundances.find(vec) == _kminmerAbundances.end()){
				if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()){
					abundances.push_back(1);
				}
				else{
					abundances.push_back(_mdbg->_dbg_nodes[vec]._abundance);
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

				if(_mdbg->_dbg_nodes.find(vec) != _mdbg->_dbg_nodes.end()) continue;

				
				bool isNewKey = false;

				_graph._mdbgSaved->_dbg_nodes.lazy_emplace_l(vec, 
				[this](MdbgNodeMap::value_type& v) { // key exist
					//if(_isFirstPass){
					//	v.second._abundance += 1;
					//}
				},           
				[&vec, this, &kminmerInfo, &readMinimizers, &readIndex, &isNewKey](const MdbgNodeMap::constructor& ctor) { // key inserted
					
					
					isNewKey = true;

					u_int32_t nodeName = -1;
					DbgNode node = {nodeName, 2};

					node._abundance = 1;
					node._quality = kminmerInfo._quality;
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

						//cout << "Indexed: " << nodeName << endl;
					}

					auto set_value = [&nodeName](MdbgNodeMap::value_type& v) { v.second._index = nodeName; };
					_graph._mdbgSaved->_dbg_nodes.modify_if(vec, set_value);
				}
				
			}

			//getchar();

		}
		

	};

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
		double _minimizerSpacingMean;
		double _kminmerLengthMean;
		double _kminmerOverlapMean;
		size_t _minimizerSize;
		size_t _kminmerSize;
		size_t _kminmerSizeFirst;
		BloomCacheCoherent<u_int64_t>* _bloomFilter;

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
			_bloomFilter = graph._bloomFilter;
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
			_bloomFilter = copy._bloomFilter;
		}

		~FilterKminmerFunctor2(){
			//delete _minimizerParser;
		}

	

		void operator () (const KminmerList& kminmerList) {

			/*
			if(_extractingContigs && _kminmerSize > 8 && kminmersInfos.size() == 0 && getAbundance(readMinimizers) > 1){

				#pragma omp critical
				{
					u_int32_t contigSize = readMinimizers.size();
					_graph._fileSmallContigs.write((const char*)&contigSize, sizeof(contigSize));

					bool isCircular = false;
					_graph._fileSmallContigs.write((const char*)&isCircular, sizeof(isCircular));
					_graph._fileSmallContigs.write((const char*)&readMinimizers[0], contigSize*sizeof(u_int64_t));
					//cout << "small contig" << endl;
					//getchar();
				}
			}
			*/

			u_int64_t readIndex = kminmerList._readIndex;
			const vector<u_int64_t>& readMinimizers = kminmerList._readMinimizers;
			//const vector<KmerVec>& kminmers = kminmerList._kminmers;
			const vector<ReadKminmerComplete>& kminmersInfos = kminmerList._kminmersInfo;

			/*
			vector<u_int32_t> abundances;
			//vector<u_int32_t> qualities;
			for(size_t i=0; i<kminmersInfos.size(); i++){
				
				const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];
				const KmerVec& vec = kminmerInfo._vec;


				if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()){
					abundances.push_back(1);
				}
				else{
					abundances.push_back(_mdbg->_dbg_nodes[vec]._abundance);
				}

				//qualities.push_back(kminmerInfo._quality);
			}

			u_int32_t median = Utils::compute_median(abundances);
			double cutoff = median * 0.1f;
			*/
		
			//float readAbundance = getReadAbundance(readMinimizers);

			int nbKminmers = 0;
			float readAbundance = getReadAbundance(readMinimizers, nbKminmers);
			double cutoff = readAbundance * 0.1f;
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
				/*
				float kminmerAbundance = -1;
				if(!_isFirstPass){
					kminmerAbundance = getAbundance(readMinimizers, kminmerInfo);
					if(kminmerAbundance == 1){//} && ab < abundanceCutoff){
						//cout << "Out: " << ab << " " << abundanceCutoff << endl;
						continue;
					}
				}
				*/

				bool isNewKey = false;

				_mdbg->_dbg_nodes.lazy_emplace_l(vec, 
				[this, &kminmerInfo](MdbgNodeMap::value_type& v) { // key exist
					v.second._quality += kminmerInfo._quality;
					//if(_isFirstPass){
					//	v.second._abundance += 1;
					//}
				},           
				[&vec, this, &kminmerInfo, &readMinimizers, &readIndex, &isNewKey, &kminmerAbundance](const MdbgNodeMap::constructor& ctor) { // key inserted
					
					
					isNewKey = true;


					u_int32_t nodeName = -1;
					DbgNode node = {nodeName, 2};

					node._abundance = kminmerAbundance; 

					node._quality = 0;
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
		
		float getReadAbundance(const vector<u_int64_t>& readMinimizers, int& nbKminmers){

			vector<u_int64_t> rlePositions;
			vector<u_int64_t> minimizers_pos;//(minimizers.size());
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

		float getKminmerAbundance(const vector<u_int64_t>& readMinimizers, const ReadKminmerComplete& kminmerInfo){

			const KmerVec& vec = kminmerInfo._vec;

			vector<u_int64_t> minimizerSeq;
			for(size_t i=kminmerInfo._read_pos_start; i<=kminmerInfo._read_pos_end; i++){
				minimizerSeq.push_back(readMinimizers[i]);
			}
			if(kminmerInfo._isReversed){
				std::reverse(minimizerSeq.begin(), minimizerSeq.end());
			}

			return getAbundance(minimizerSeq);
		}

		float getAbundance(const vector<u_int64_t>& readMinimizers){

			vector<u_int64_t> rlePositions;
			vector<u_int64_t> minimizers_pos;//(minimizers.size());
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

};	


#endif 

