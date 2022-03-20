
//Real data
//./bin/Bloocoo  -i ~/workspace/data/overlap_test/input.txt -l 21 -k 3 -d 0.005 -o ~/workspace/run/overlap_test/ -ihifiasm ~/workspace/run/hifiasm_meta/AD_components/big/component_3.fasta -idir ~/workspace/run/overlap_test/

//Simulation
///bin/Bloocoo  -i ~/workspace/data/overlap_test/input.txt -l 16 -k 3 -d 0.005 -o ~/workspace/run/overlap_test_3/ -ihifiasm ~/workspace/data/overlap_test/genome_2371_20x/truth_input.txt -idir ~/workspace/run/overlap_test_3/

#ifndef _BLOOCOO_HPP_
#define _BLOOCOO_HPP_

#include "Commons.hpp"

//#include "graph/Graph.hpp"
#include "graph/GfaParser.hpp"
#include "graph/GraphSimplify.hpp"
#include "utils/BloomFilter.hpp"
//#include "./utils/ntHashIterator.hpp"

//#include <boost/graph/adjacency_list.hpp>
//#include <boost/graph/graph_utility.hpp>
//#include <boost/property_map/property_map.hpp>
//#include <boost/graph/connected_components.hpp>
//#include <boost/graph/dijkstra_shortest_paths.hpp>












class Bloocoo : public Tool{
    
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
	string _inputFilename;
	string _outputDir;
	float _minimizerDensity;
    size_t _minimizerSize;
    size_t _kminmerSize;
    size_t _kminmerSizeFirst;
	string _filename_noKminmerReads;
	ofstream _file_noKminmerReads;
	int _nbCores;
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
	//string _filename_filteredMinimizers;
	//string _filename_readCompositions;

	//vector<u_int32_t> _evaluation_readToDataset;
	bool _parseReads;
	MDBG* _mdbg;
	MDBG* _mdbgInit;
	MDBG* _mdbgNoFilter;
	bool _parsingContigs;
	//MinimizerPairMap* _minimizerPairMap;

	double _minimizerSpacingMean;
	double _kminmerLengthMean;
	double _kminmerOverlapMean;

	double _minimizerSpacing_sum;
	double _minimizerSpacing_n;
	double _kminmerLength_sum;
	double _kminmerLength_n;

    Bloocoo ();
    void execute ();
	void createMDBG();
	//void createMDBG_read(kseq_t* read, u_int64_t readIndex);
	void createMDBG_collectKminmers(const vector<u_int64_t>& minimizers, const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex);
	void createMDBG_collectKminmers_contig(const vector<u_int64_t>& minimizers, const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex);
	void removeErroneousKminmers(const vector<u_int64_t>& minimizers, const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex);
	void loadSolidKminmers();

	void createMDBG_collectKminmers_read(kseq_t* read, u_int64_t readIndex);
	void extractKminmerSequence(const char* sequenceOriginal, const ReadKminmer& kminmerInfo, string& sequence);


	void createMDBG_index(const vector<u_int64_t>& minimizers, const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex);
	void createGfa();
	void parseArgs(int argc, char* argv[]);
	void parseContigs();
	float computeKmerVecAbundance(const vector<u_int64_t>& minimizers, bool isContig);
	void computeContigAbundance();
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




	class IndexKminmerFunctor {

		public:

		Bloocoo& _graph;
		bool _isFirstPass;
		bool _parsingContigs;
		ofstream& _readFile;
		unordered_set<KmerVec>& _kminmerExist;
		MDBG* _mdbg;
		MDBG* _mdbgInit;
		ofstream& _kminmerFile;
		double _minimizerSpacingMean;
		double _kminmerLengthMean;
		double _kminmerOverlapMean;
		size_t _minimizerSize;
		size_t _kminmerSize;
		size_t _kminmerSizeFirst;

		IndexKminmerFunctor(Bloocoo& graph) : _graph(graph), _readFile(graph._readFile), _kminmerExist(graph._kminmerExist), _kminmerFile(graph._kminmerFile){
			_isFirstPass = graph._isFirstPass;
			_parsingContigs = graph._parsingContigs;
			_mdbg = graph._mdbg;
			_mdbgInit = graph._mdbgInit;
			_minimizerSpacingMean = graph._minimizerSpacingMean;
			_kminmerLengthMean = graph._kminmerLengthMean;
			_kminmerOverlapMean = graph._kminmerOverlapMean;
			_minimizerSize = graph._minimizerSize;
			_kminmerSize = graph._kminmerSize;
			_kminmerSizeFirst = graph._kminmerSizeFirst;
		}

		IndexKminmerFunctor(const IndexKminmerFunctor& copy) : _graph(copy._graph), _readFile(copy._readFile), _kminmerExist(copy._kminmerExist), _kminmerFile(copy._kminmerFile){
			_isFirstPass = copy._isFirstPass;
			_parsingContigs = copy._parsingContigs;
			_mdbg = copy._mdbg;
			_mdbgInit = copy._mdbgInit;
			_minimizerSpacingMean = copy._minimizerSpacingMean;
			_kminmerLengthMean = copy._kminmerLengthMean;
			_kminmerOverlapMean = copy._kminmerOverlapMean;
			_minimizerSize = copy._minimizerSize;
			_kminmerSize = copy._kminmerSize;
			_kminmerSizeFirst = copy._kminmerSizeFirst;
		}

		~IndexKminmerFunctor(){
			//delete _minimizerParser;
		}

		void operator () (const KminmerList& kminmerList) {


			u_int64_t readIndex = kminmerList._readIndex;
			const vector<u_int64_t>& readMinimizers = kminmerList._readMinimizers;
			//const vector<KmerVec>& kminmers = kminmerList._kminmers;
			const vector<ReadKminmerComplete>& kminmersInfos = kminmerList._kminmersInfo;

			//vector<u_int16_t> minimizerPosOffset(readMinimizers.size());
			//for(size_t i=0; i<minimizerPosOffset.size(); i++){
			//	minimizerPosOffset[i] = _minimizerSpacingMean;
			//}

			#pragma omp critical
			{
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

					if(_kminmerExist.find(vec) != _kminmerExist.end() || _parsingContigs){

						if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()){

							_mdbg->addNode(vec, _kminmerLengthMean, _minimizerSpacingMean, _minimizerSpacingMean, kminmerInfo._isReversed);

							vector<u_int64_t> minimizerSeq;
							for(size_t i=kminmerInfo._read_pos_start; i<=kminmerInfo._read_pos_end; i++){
								minimizerSeq.push_back(readMinimizers[i]);
							}
							if(kminmerInfo._isReversed){
								std::reverse(minimizerSeq.begin(), minimizerSeq.end());
							}

							u_int16_t size = minimizerSeq.size();
							_kminmerFile.write((const char*)&size, sizeof(size));
							_kminmerFile.write((const char*)&minimizerSeq[0], size*sizeof(uint64_t));

							u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
							u_int32_t length = kminmerInfo._read_pos_end - kminmerInfo._read_pos_start;
							u_int32_t lengthStart = kminmerInfo._seq_length_start;
							u_int32_t lengthEnd = kminmerInfo._seq_length_end;
							//bool isReversed = kminmerInfo._isReversed;

							_kminmerFile.write((const char*)&nodeName, sizeof(nodeName));
							_kminmerFile.write((const char*)&length, sizeof(length));
							_kminmerFile.write((const char*)&lengthStart, sizeof(lengthStart));
							_kminmerFile.write((const char*)&lengthEnd, sizeof(lengthEnd));
							//_kminmerFile.write((const char*)&isReversed, sizeof(isReversed));

							//if(_mdbg->_dbg_nodes[vec]._index <= 0){
							//	cout << nodeName << " " << lengthStart << " " << lengthEnd << " " << isReversed << " " << kminmerSequence.size() << endl;
							//	cout << kminmerSequence << endl;
							//}

							if(_isFirstPass){
								_mdbg->_dbg_nodes[vec]._abundance = 2;
							}
							else{
								vector<u_int64_t> rlePositions;
								vector<u_int64_t> minimizers_pos;//(minimizers.size());
								vector<KmerVec> kminmers; 
								vector<ReadKminmer> kminmersInfo;
								MDBG::getKminmers(_minimizerSize, _kminmerSizeFirst, minimizerSeq, minimizers_pos, kminmers, kminmersInfo, rlePositions, readIndex, false);


								u_int32_t minAbundance = -1;
								float sum = 0;
								float n = 0;
								vector<float> abundances;
								//cout << kminmers.size() << endl;
								for(const KmerVec& vec : kminmers){
									//cout << (_mdbgInit->_dbg_nodes.find(vec) != _mdbgInit->_dbg_nodes.end()) << endl;
									if(_mdbgInit->_dbg_nodes.find(vec) != _mdbgInit->_dbg_nodes.end()){

										u_int32_t abundance = _mdbgInit->_dbg_nodes[vec]._abundance;
										if(abundance < minAbundance){
											minAbundance = abundance;
										}
										//cout << abundance << " ";
										abundances.push_back(abundance);

										sum += abundance;
										n += 1;
									}
								}
								//cout << endl;

								//cout << (sum / 2) << " " << Utils::compute_median_float(abundances) << endl;

								_mdbg->_dbg_nodes[vec]._abundance = minAbundance;
								//cout << minAbundance << endl;
							}

						}
						else{
							if(_isFirstPass){
								_mdbg->addNode(vec, _kminmerLengthMean, _minimizerSpacingMean, _minimizerSpacingMean, kminmerInfo._isReversed);
							}
						}

					}
					else{
						_kminmerExist.insert(vec);
					}


				}
			}

		}
	};
	

};	


#endif 

