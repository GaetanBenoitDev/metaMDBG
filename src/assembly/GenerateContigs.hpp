

#ifndef MDBG_METAG_GENERATECONTIGS
#define MDBG_METAG_GENERATECONTIGS

//#define PRINT_DEBUG_PREVRANK

#include "Commons.hpp"
#include "graph/ProgressiveAbundanceFilter.hpp"

//#include <string>
//#include <sstream>
//#include <unordered_set>
//#include <unordered_map>
//#include <regex>
//#include <algorithm>
//#include <libgen.h>
//#include <set>
//#include "graph/Graph.hpp"
//#include "graph/GfaParser.hpp"
#include "graph/GraphSimplify.hpp"
//#include "eval/ContigStatistics.hpp"
//#include "toBasespace/ToBasespaceOnTheFly.hpp"
//#include "contigFeatures/ContigFeature.hpp"










class GenerateContigs : public Tool{

public:

	string _inputDir;
	string _truthInputFilename;
	string _outputFilename;
	string _outputFilename_complete;
	bool _debug;
	//bool _isFinalAssembly;
	string _inputFilename_unitigNt;
	string _inputFilename_unitigCluster;
	string _filename_abundance;

	float _minimizerDensity;
    size_t _minimizerSize;
    size_t _kminmerSize;
    size_t _kminmerSizeFirst;
	float _minimizerSpacingMean;
	float _kminmerLengthMean;
	float _kminmerOverlapMean;

	string _partitionDir;
	string _filename_solidKminmers;
	string _filename_inputContigs;
	string _filename_outputContigs;
	string _filename_readMinimizers;
	string _filename_hifiasmGroundtruth;
	unordered_map<KmerVec, u_int16_t> _evaluation_hifiasmGroundTruth;
	unordered_map<KmerVec, u_int32_t> _evaluation_hifiasmGroundTruth_position;
	unordered_map<u_int32_t, u_int32_t> _evaluation_hifiasmGroundTruth_nodeNamePosition;
	gzFile _outputContigFile;
	gzFile _outputContigFile_complete;
	MDBG* _mdbg;
	//GraphSimplify* _graph;
	//ToBasespaceOnTheFly _toBasespace;
	//ContigFeature _contigFeature;
	string _gfaFilename;
	size_t _nbCores;

    vector<u_int32_t> _nodeName_to_abundance;

	GenerateContigs(): Tool (){

	}

	~GenerateContigs(){

	}


	void parseArgs(int argc, char* argv[]){

		args::ArgumentParser parser("contig", ""); //"This is a test program.", "This goes after the options."
		args::Positional<std::string> arg_outputDir(parser, "outputDir", "", args::Options::Required);
		//args::Positional<std::string> arg_contigs(parser, "contigs", "", args::Options::Required);
		//args::PositionalList<std::string> arg_readFilenames(parser, "reads", "Input filename(s) (separated by space)", args::Options::Required);
		//args::ValueFlag<int> arg_l(parser, "", "Minimizer length", {ARG_MINIMIZER_LENGTH2}, 13);
		//args::ValueFlag<float> arg_d(parser, "", "Minimizer density", {ARG_MINIMIZER_DENSITY2}, 0.005f);
		//args::ValueFlag<std::string> arg_contigs(parser, "", "", {ARG_INPUT_FILENAME_CONTIG}, "");
		args::ValueFlag<int> arg_nbCores(parser, "", "Number of cores", {ARG_NB_CORES2}, NB_CORES_DEFAULT_INT);
		//args::Flag arg_isFinalAssembly(parser, "", "Is final multi-k pass", {ARG_FINAL});
		args::Flag arg_firstPass(parser, "", "Is first pass of multi-k", {ARG_FIRST_PASS});
		args::Flag arg_help(parser, "", "", {'h', "help"}, args::Options::Hidden);
		//args::HelpFlag help(parser, "help", "Display this help menu", {'h'});
		//args::CompletionFlag completion(parser, {"complete"});

		//(ARG_INPUT_FILENAME_TRUTH, "", cxxopts::value<string>()->default_value(""))
		//(ARG_MINIMIZER_LENGTH, "", cxxopts::value<int>()->default_value("13"))
		//(ARG_MINIMIZER_DENSITY, "", cxxopts::value<float>()->default_value("0.005"))
		//(ARG_NB_CORES, "", cxxopts::value<int>()->default_value(NB_CORES_DEFAULT))
		//(ARG_BLOOM_FILTER, "", cxxopts::value<bool>()->default_value("false"));

		try
		{
			parser.ParseCLI(argc, argv);
		}
		catch (const std::exception& e)
		{
			cerr << parser;
			cerr << e.what() << endl;
			exit(0);
		}

		if(arg_help){
			cerr << parser;
			exit(0);
		}

		_inputDir = args::get(arg_outputDir);
		//_filename_inputContigs = args::get(arg_contigs);

		_nbCores = args::get(arg_nbCores);

		//_isFinalAssembly = false;
		//if(arg_isFinalAssembly){
		//	_isFinalAssembly = true;
		//}

		//_truthInputFilename = result[ARG_INPUT_FILENAME_TRUTH].as<string>();
		//_inputFilename_unitigNt = result[ARG_INPUT_FILENAME_UNITIG_NT].as<string>();
		//_inputFilename_unitigCluster = result[ARG_INPUT_FILENAME_UNITIG_CLUSTER].as<string>();
		//_debug = result[ARG_DEBUG].as<bool>();
		//_isFinalAssembly = result[ARG_FINAL].as<bool>();
		//_nbCores = result[ARG_NB_CORES].as<int>();
		//_nbCores = 1;
		//_filename_abundance = result[ARG_INPUT_FILENAME_ABUNDANCE].as<string>();


		/*
		cxxopts::Options options("Assembly", "");
		options.add_options()
		(ARG_OUTPUT_DIR, "", cxxopts::value<string>())
		(ARG_INPUT_FILENAME_CONTIG, "", cxxopts::value<string>()->default_value(""))
		(ARG_INPUT_FILENAME_TRUTH, "", cxxopts::value<string>()->default_value(""))
		//(ARG_DEBUG, "", cxxopts::value<bool>()->default_value("false"))
		//(ARG_INPUT_FILENAME_UNITIG_NT, "", cxxopts::value<string>()->default_value(""))
		//(ARG_INPUT_FILENAME_UNITIG_CLUSTER, "", cxxopts::value<string>()->default_value(""))
		(ARG_FINAL, "", cxxopts::value<bool>()->default_value("false"))
		(ARG_NB_CORES, "", cxxopts::value<int>()->default_value(NB_CORES_DEFAULT));
		//(ARG_INPUT_FILENAME_ABUNDANCE, "", cxxopts::value<string>()->default_value(""));



		if(argc <= 1){
			cout << options.help() << endl;
			exit(0);
		}

		cxxopts::ParseResult result;

		try{
			result = options.parse(argc, argv);

			_inputDir = result[ARG_OUTPUT_DIR].as<string>();
			_filename_inputContigs = result[ARG_INPUT_FILENAME_CONTIG].as<string>();
			_truthInputFilename = result[ARG_INPUT_FILENAME_TRUTH].as<string>();
			//_inputFilename_unitigNt = result[ARG_INPUT_FILENAME_UNITIG_NT].as<string>();
			//_inputFilename_unitigCluster = result[ARG_INPUT_FILENAME_UNITIG_CLUSTER].as<string>();
			//_debug = result[ARG_DEBUG].as<bool>();
			_isFinalAssembly = result[ARG_FINAL].as<bool>();
			_nbCores = result[ARG_NB_CORES].as<int>();
			//_nbCores = 1;
			//_filename_abundance = result[ARG_INPUT_FILENAME_ABUNDANCE].as<string>();
		}
		catch (const std::exception& e){
			std::cout << options.help() << std::endl;
			std::cerr << e.what() << std::endl;
			std::exit(EXIT_FAILURE);
		}
		*/


		string filename_parameters = _inputDir + "/parameters.gz";
		gzFile file_parameters = gzopen(filename_parameters.c_str(),"rb");
		gzread(file_parameters, (char*)&_minimizerSize, sizeof(_minimizerSize));
		gzread(file_parameters, (char*)&_kminmerSize, sizeof(_kminmerSize));
		gzread(file_parameters, (char*)&_minimizerDensity, sizeof(_minimizerDensity));
		gzread(file_parameters, (char*)&_kminmerSizeFirst, sizeof(_kminmerSizeFirst));
		gzread(file_parameters, (char*)&_minimizerSpacingMean, sizeof(_minimizerSpacingMean));
		gzread(file_parameters, (char*)&_kminmerLengthMean, sizeof(_kminmerLengthMean));
		gzread(file_parameters, (char*)&_kminmerOverlapMean, sizeof(_kminmerOverlapMean));
		gzclose(file_parameters);

		openLogFile(_inputDir);

		Logger::get().debug() << "";
		Logger::get().debug() << "Input dir: " << _inputDir;
		//cout << "Output filename: " << _outputFilename << endl;
		Logger::get().debug() << "Minimizer length: " << _minimizerSize;
		Logger::get().debug() << "Kminmer length: " << _kminmerSize;
		Logger::get().debug() << "Density: " << _minimizerDensity;
		Logger::get().debug() << "";

		_outputFilename = _inputDir + "/minimizer_contigs.gz";
		_outputFilename_complete = _inputDir + "/minimizer_contigs_complete.gz";
		_filename_readMinimizers = _inputDir + "/read_data.txt";
		_filename_hifiasmGroundtruth = _inputDir + "/hifiasmGroundtruth.gz";
		_filename_outputContigs = _inputDir + "/contigs.min.gz";
		_filename_solidKminmers = _inputDir + "/solid.min.gz";

	}

	unordered_map<u_int64_t, vector<u_int32_t>> _debug_readPath;
    vector<bool> _isBubble;

	class X
	{
		public: 
		void operator()(float cutoff) {}
	};

	void execute (){

		//cout << "Loading abundant nodes" << endl;
		//loadAbundantNodes();

		loadGraph();
		//exit(1);
		//generateUnitigs();
		//generateContigs2(_inputDir + "/contigs.nodepath");
		generateContigs3();

		//exit(1);
		
		//if(_kminmerSize == 4){
			_mdbg = new MDBG(_kminmerSize);
			_mdbg->load(_inputDir + "/kminmerData_min.txt", false);
			for(auto& it : _mdbg->_dbg_nodes){
				if(_nodeNameAbundances.find(it.second._index) == _nodeNameAbundances.end()) continue;
				const NodeAb& nodeAb = _nodeNameAbundances[it.second._index];
				it.second._abundance = ceil(nodeAb._abundance);

				//"reprise: essayer d'enelever pas ça pas sur que c'est correct en fait"
				if(nodeAb._nbNodes > _kminmerSize) it.second._abundance = max(it.second._abundance, (u_int32_t)2);

				//"generate contig finale: besoin de l'ancien systeme de cleaning avec prise en comtpe du coverage"
				//it.second._unitigNbNodes = nodeAb._nbNodes;
				//cout << nodeAb._abundance << " " << nodeAb._nbNodes << endl;
			}
			_mdbg->dump(_inputDir + "/kminmerData_min.txt");
		//}

		//closeLogFile();
	}

	unordered_set<u_int32_t> _abundantNodeNames;

	void loadAbundantNodes(){
		
		ifstream file(_inputDir + "/abundantNodes_backup.bin");


		while (true) {

			u_int32_t nodeName1;
			u_int32_t nodeName2;

			file.read((char*)&nodeName1, sizeof(nodeName1));

			if(file.eof()) break;

			file.read((char*)&nodeName2, sizeof(nodeName2));

			_abundantNodeNames.insert(nodeName1);
			_abundantNodeNames.insert(nodeName2);

		}

		file.close();

	}

	void loadGraph(){

		_gfaFilename = _inputDir + "/minimizer_graph.gfa";
		string gfa_filename_noUnsupportedEdges = _inputDir + "/minimizer_graph_noUnsupportedEdges.gfa";
		//string gfa_filename_unitigs = _inputDir + "/minimizer_graph_unitigs.gfa";
		string mdbg_filename = _inputDir + "/mdbg_nodes.gz";




		
		//cout << _gfaFilename << endl;
		//_mdbg = new MDBG(_kminmerSize);
		//_mdbg->load(mdbg_filename);
		//cout << "Nb nodes: " <<  _mdbg->_dbg_nodes.size() << endl;

		//extractContigKminmers2(_inputDir + "/contigs.fasta.gz");

		//if(_truthInputFilename != ""){
		//	extract_truth_kminmers();
		//}


		//file_groundTruth = ofstream(_inputDir + "/binning_results.csv");
		//file_groundTruth << "Name,Colour" << endl;

		//file_groundTruth_hifiasmContigs = ofstream(_inputDir + "/binning_results_hifiasm_" + to_string(_kminmerSize) + ".csv");
		//file_groundTruth_hifiasmContigs << "Name,Colour" << endl;

		//if(_debug){
            //gfa_filename = _inputDir + "/minimizer_graph_debug.gfa";
		//}
		
		BiGraph* graphSuccessors = GfaParser::createBiGraph_binary(_gfaFilename, true, 0, _kminmerLengthMean, _kminmerOverlapMean);
		GraphSimplify* graphSimplify = new GraphSimplify(graphSuccessors, _inputDir, _nbCores);
		//_graph = graphSimplify;
		
		u_int64_t nbNodes = graphSuccessors->_nodes.size();


		//_graph->clear(0);
		//_graph->compact(false, _unitigDatas);
		//_graph->removeErrors_4(_kminmerSize, _unitigDatas);

		//Generate unitigs
		//cout << "Indexing reads" << endl;
		//_unitigDatas.resize(_mdbg->_dbg_nodes.size());
		//_graph->clear(0);
		//_graph->compact(false, _unitigDatas);
		//removeUnsupportedEdges(_gfaFilename, gfa_filename_noUnsupportedEdges, _graph);
		

		//cout << "done" << endl;
	
    //void debug_writeGfaErrorfree(unitigDatas, crushBubble, smallBubbleOnly, detectRoundabout, insertBubble, saveAllState, doesSaveUnitigGraph, MDBG* mdbg, size_t minimizerSize, size_t nbCores, bool useLocalAbundanceFilter, bool removeLongTips){


		//if(_kminmerSize > 4){
		//GfaParser::binaryGraph_to_gfa(_gfaFilename, _kminmerLengthMean, _kminmerOverlapMean, _gfaFilename+".gfa", _graph->_graphSuccessors->_nodeDatas);
		//if(_kminmerSize > 45){
		//	getchar();
		//}
		//	cout << "lala" << endl;
		//	getchar();
		//}

	  //_graph->debug_writeGfaErrorfree(0, 0, -1, _kminmerSize, false, true, false, _unitigDatas, true, false, true, false, true, true, _mdbg, _minimizerSize, _nbCores, true, false);
		//_graph->debug_writeGfaErrorfree(0, 0, -1, _kminmerSize, false, true, false, _unitigDatas, true, false, false, false, false, true, _mdbg, _minimizerSize, _nbCores, false, false);
		
		//indexReads();


		graphSimplify->debug_writeGfaErrorfree();
		//removeUnsupportedEdges();
		//_graph->_unitigGraph->save("/home/gats/workspace/run/minimizer_graph_u.gfa");
		//getchar();
		        
		
		delete graphSimplify;
		delete graphSuccessors;



		Logger::get().debug() << "Load node data";
		ifstream kminmerFile(_inputDir + "/kminmerData_min.txt");
        _nodeName_to_abundance.resize(nbNodes/2, 0);


		while (true) {

			vector<MinimizerType> minimizerSeq;
			minimizerSeq.resize(_kminmerSize);
			kminmerFile.read((char*)&minimizerSeq[0], minimizerSeq.size()*sizeof(MinimizerType));

			if(kminmerFile.eof()) break;

			u_int32_t nodeName;
			u_int32_t abundance;
			//u_int32_t quality;
			//bool isReversed = false;

			kminmerFile.read((char*)&nodeName, sizeof(nodeName));
			kminmerFile.read((char*)&abundance, sizeof(abundance));
			//kminmerFile.read((char*)&quality, sizeof(quality));


			//cout << nodeName << " " << abSum << endl;

			//cout << nodeName << " " << abundance << " " << quality << endl;
			//if(quality == 0){
			//	cout << quality << endl;
			//	getchar();
			//}
			//if(abundance == 1) continue;
			//KmerVec vec;
			//vec._kmers = minimizerSeq;

			_nodeName_to_abundance[nodeName] = abundance;
			//_graph->_graphSuccessors->_nodeLengths[nodeName] = _kminmerLengthMean;
			//_kminmerAbundances[vec] = abundance;
		}

		//cout << abSum << " " << qualSum << endl;
		//getchar();

		kminmerFile.close();


		Logger::get().debug() << "Loading unitig graph";
		UnitigGraph* unitigGraph = new UnitigGraph(_kminmerLengthMean, _kminmerOverlapMean, _kminmerLengthMean-_kminmerOverlapMean);
		unitigGraph->load(_inputDir + "/unitigGraph.nodes.bin", _inputDir + "/unitigGraph.edges.successors.bin", _inputDir + "/unitigGraph.edges.predecessors.bin", _nodeName_to_abundance);

        Logger::get().debug() << "Nb unitigs: " << unitigGraph->_nodes.size();
        //_logFile << "Nb edges: " << nbEdgeSuccessors<< endl;

		Logger::get().debug() << "Simplifying";

        _progressiveAbundanceFilter = new ProgressiveAbundanceFilter(unitigGraph, _inputDir, _kminmerSize, true, _kminmerSize==_kminmerSizeFirst);
		X dummy;

		//delete _graph->_graphSuccessors;

        _progressiveAbundanceFilter->execute(dummy);

		//cout << _graph->_unitigGraph->computeChecksum_successor() << endl;
		//if(_kminmerSize == 20) exit(1);


        //_unitigGraph->save("/home/gats/workspace/run/unitigGraph_cleaned.gfa");

		
		//_graph->debug_selectUnitigIndex();
		//_graph->debug_writeGfaErrorfree(2000, 2000, -1, _kminmerSize, false, true, false, _unitigDatas, true, false, false, false, false, false, _mdbg, _minimizerSize, _nbCores, false, true);
			
		//delete _mdbg;
		
		_processedNodeNames = vector<bool>(nbNodes/2, false);
	}

	ProgressiveAbundanceFilter* _progressiveAbundanceFilter;

	/*
	u_int64_t computeSharedReads(UnitigGraph::Node* unitig1, UnitigGraph::Node* unitig2){

		const vector<u_int32_t>& source_readIndexes = _nodeName_to_readIndexes[BiGraph::nodeIndex_to_nodeName(unitig1->endNode())];
		const vector<u_int32_t>& readIndex = _nodeName_to_readIndexes[BiGraph::nodeIndex_to_nodeName(unitig2->startNode())];


		size_t i=0;
		size_t j=0;
		u_int64_t nbShared = 0;

		while(i < source_readIndexes.size() && j < readIndex.size()){

			if(source_readIndexes[i] == readIndex[j]){
				nbShared += 1;
				i += 1;
				j += 1;
			}
			else if(source_readIndexes[i] < readIndex[j]){
				i += 1;
			}
			else{
				j += 1;
			}

		}

		return nbShared;
	}

	//vector<vector<ReadNodePos>> _nodeName_to_readIndexes;
	vector<vector<u_int32_t>> _nodeName_to_readIndexes;
	unordered_set<u_int32_t> _startEndNodeNames;

	void removeUnsupportedEdges(){

		cout << "Indexing reads" << endl;


		for(UnitigGraph::Node* node : _graph->_unitigGraph->_nodes){
			if(node->_unitigIndex == -1) continue;
			_startEndNodeNames.insert(BiGraph::nodeIndex_to_nodeName(node->startNode()));
			_startEndNodeNames.insert(BiGraph::nodeIndex_to_nodeName(node->endNode()));
		}

		_nodeName_to_readIndexes.resize(_graph->_graphSuccessors->_nbNodes/2);

		string mdbg_filename = _inputDir + "/kminmerData_min.txt";
		_mdbg = new MDBG(_kminmerSize);
		_mdbg->load(mdbg_filename, false);

		KminmerParserParallel parser(_filename_readMinimizers, _minimizerSize, _kminmerSize, false, false, _nbCores);
		parser.parse(IndexReadsFunctor(*this));

		delete _mdbg;
		_startEndNodeNames.clear();


		for(UnitigGraph::Node* node : _graph->_unitigGraph->_nodes){
			if(node->_unitigIndex == -1) continue;


			vector<UnitigGraph::Node*> preds = node->_predecessors;

			for(UnitigGraph::Node* predecessor : preds){
				
				if(predecessor->_unitigIndex == -1) continue;
		
				u_int64_t nbSharedReads = computeSharedReads(predecessor, node);
				if(nbSharedReads > 0) continue;

				_graph->_unitigGraph->removePredecessor(node, predecessor);
				//recompactNodes.insert(predecessor);
				_graph->_unitigGraph->recompact(predecessor);
				cout << "Remove unsupported edge" << endl;
				//getchar();

			}

		}

		cout << "done" << endl;
		_nodeName_to_readIndexes.clear();
	}

	unordered_set<u_int32_t> _validReadIndex;
	unordered_set<u_int32_t> _isNodeNamePalindrome;

	class IndexReadsFunctor {

		public:

		GenerateContigs& _parent;

		IndexReadsFunctor(GenerateContigs& parent) : _parent(parent){
		}

		IndexReadsFunctor(const IndexReadsFunctor& copy) : _parent(copy._parent){
		}

		~IndexReadsFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {


			u_int32_t readIndex = kminmerList._readIndex;


			for(u_int16_t i=0; i<kminmerList._kminmersInfo.size(); i++){
			
				
				const ReadKminmerComplete& kminmerInfo = kminmerList._kminmersInfo[i];

				KmerVec vec = kminmerInfo._vec;
				
				if(_parent._mdbg->_dbg_nodes.find(vec) == _parent._mdbg->_dbg_nodes.end()){
					continue;
				}
				u_int32_t nodeName = _parent._mdbg->_dbg_nodes[vec]._index;
				
				if(_parent._startEndNodeNames.find(nodeName) == _parent._startEndNodeNames.end()) continue;

				//u_int32_t nodeIndex = BiGraph::nodeName_to_nodeIndex(nodeName, !kminmerInfo._isReversed);
				
				#pragma omp critical(indexReadsFunctor)
				{

					_parent._nodeName_to_readIndexes[nodeName].push_back(readIndex);
					//cout << nodeName << " " << _parent._nodeName_to_readIndexes[nodeName].size() << endl;
				}

				

			}

		}
	};
	*/

	/*
	phmap::parallel_flat_hash_set<vector<u_int32_t>> _nodeName_to_readIndexes;

	void indexReads(){

		for(u_int32_t nodeIndex=0; nodeIndex < _graph->_graphSuccessors->_nbNodes; nodeIndex++){
			const vector<AdjNode>& successors = _graph->_graphSuccessors->_nodes[nodeIndex];
			if(successors.size() > 1){
				_nodeName_to_readIndexes[BiGraph::nodeIndex_to_nodeName(nodeIndex)] = {};
				for(const AdjNode& node : successors){
					_nodeName_to_readIndexes[BiGraph::nodeIndex_to_nodeName(node._index)] = {};
				}
			}

		}

		//_nodeName_to_readIndexes.resize(_graph->_graphSuccessors->_nbNodes/2);

		_logFile << "Indexing reads" << endl;
		string mdbg_filename = _inputDir + "/kminmerData_min.txt";
		_mdbg = new MDBG(_kminmerSize);
		_mdbg->load(mdbg_filename, false);

		KminmerParserParallel parser(_filename_readMinimizers, _minimizerSize, _kminmerSize, false, false, _nbCores);
		parser.parse(IndexReadsFunctor(*this));

		delete _mdbg;

		for(u_int32_t nodeIndex=0; nodeIndex < _graph->_graphSuccessors->_nbNodes; nodeIndex++){
			const vector<AdjNode>& successors = _graph->_graphSuccessors->_nodes[nodeIndex];
			if(successors.size() > 1){
				for(const AdjNode& node : successors){
					u_int32_t nodeIndex_nn = node._index;
					if(Utils::shareAny(_nodeName_to_readIndexes[BiGraph::nodeIndex_to_nodeName(nodeIndex)], _nodeName_to_readIndexes[BiGraph::nodeIndex_to_nodeName(nodeIndex_nn)])){
						continue;
					}
					cout << "unsupported edge: " << BiGraph::nodeIndex_to_nodeName(nodeIndex) << " " << BiGraph::nodeIndex_to_nodeName(nodeIndex_nn) << endl;
					getchar();
				}
			}


		}
	}


	class IndexReadsFunctor {

		public:

		GenerateContigs& _generateContigs;

		IndexReadsFunctor(GenerateContigs& generateContigs) : _generateContigs(generateContigs){
		}

		IndexReadsFunctor(const IndexReadsFunctor& copy) : _generateContigs(copy._generateContigs){
		}

		~IndexReadsFunctor(){
		}

		//void collectBestSupportingReads_read(const vector<u_int64_t>& readMinimizers, const vector<ReadKminmerComplete>& kminmersInfos, u_int8_t isCircular, u_int64_t readIndex){

		void operator () (const KminmerList& kminmerList) {

			u_int64_t readIndex = kminmerList._readIndex;

			for(size_t i=0; i<kminmerList._kminmersInfo.size(); i++){
			
				
				//_logFile << readIndex << " " << i << endl;
				const ReadKminmerComplete& kminmerInfo = kminmerList._kminmersInfo[i];

				KmerVec vec = kminmerInfo._vec;
				
				if(_generateContigs._mdbg->_dbg_nodes.find(vec) == _generateContigs._mdbg->_dbg_nodes.end()){
					continue;
				}
				
				u_int32_t nodeName = _generateContigs._mdbg->_dbg_nodes[vec]._index;

				if(_nodeName_to_readIndexes.find(nodeName) == _nodeName_to_readIndexes.end()) continue;

				#pragma omp critical(indexReadsFunctor)
				{
					_generateContigs._nodeName_to_readIndexes[nodeName].push_back(readIndex);
				}

				

			}

		}
	};
	*/
	

	/*
	void indexReads_read(const vector<u_int64_t>& minimizers, const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex){//}, const vector<KmerVec>& kminmers_k3, const vector<ReadKminmer>& kminmersInfos_k3){

		//vector<ReadIndexType> unitigIndexex;

		for(const KmerVec& vec : kminmers){
			//if(mdbg->_dbg_nodes[vec]._index == 55479) cout << "AAAAA" << endl;

			//cout << mdbg->_dbg_nodes[vec]._index << endl;
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()) continue;

			//if(vec.isPalindrome()) cout << _mdbg->_dbg_nodes[vec]._index << endl;
			//getchar();


			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
			u_int32_t nodeIndex = BiGraph::nodeName_to_nodeIndex(nodeName, true);
			//if(_graph->_isNodeValid2.find(nodeIndex) == _graph->_isNodeValid2.end()) continue;
			//if(_nodeData.find(nodeName) == _nodeData.end()) continue;

			//UnitigData& unitigData = _nodeData[nodeName];
			UnitigData& unitigData = _unitigDatas[nodeName];
			//if(std::find(unitigData._readIndexes.begin(), unitigData._readIndexes.end(), readIndex) != unitigData._readIndexes.end()) continue;
			unitigData._readIndexes.push_back(readIndex);
			//cout << "indexing : " << readIndex << endl;
		}



	}

	void indexReads_contig(const vector<u_int64_t>& minimizers, const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex){//}, const vector<KmerVec>& kminmers_k3, const vector<ReadKminmer>& kminmersInfos_k3){

		readIndex += 2000000000;
		//vector<ReadIndexType> unitigIndexex;

		for(const KmerVec& vec : kminmers){
			//if(mdbg->_dbg_nodes[vec]._index == 55479) cout << "AAAAA" << endl;

			//cout << mdbg->_dbg_nodes[vec]._index << endl;
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()) continue;

			//if(vec.isPalindrome()) cout << _mdbg->_dbg_nodes[vec]._index << endl;
			//getchar();


			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
			u_int32_t nodeIndex = BiGraph::nodeName_to_nodeIndex(nodeName, true);
			//if(_graph->_isNodeValid2.find(nodeIndex) == _graph->_isNodeValid2.end()) continue;
			//if(_nodeData.find(nodeName) == _nodeData.end()) continue;

			//UnitigData& unitigData = _nodeData[nodeName];
			UnitigData& unitigData = _unitigDatas[nodeName];
			//if(std::find(unitigData._readIndexes.begin(), unitigData._readIndexes.end(), readIndex) != unitigData._readIndexes.end()) continue;
			unitigData._readIndexes.push_back(readIndex);
			//cout << "indexing : " << readIndex << endl;
		}



	}
	
	void removeUnsupportedEdges(const string& gfaFilename, const string& gfa_filename_noUnsupportedEdges, GraphSimplify* graph){

		KminmerParser parser(_filename_readMinimizers, _minimizerSize, _kminmerSize, false);
		//auto fp = std::bind(&Assembly::indexReads_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
		auto fp = std::bind(&GenerateContigs::indexReads_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
		parser.parse(fp);
		
		//if(_filename_inputContigs != ""){
		//	KminmerParser parserContig(_filename_inputContigs, _minimizerSize, _kminmerSize, true);
		//	auto fpContig = std::bind(&Assembly3::indexReads_contig, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
		//	parserContig.parse(fpContig);
		//}

		if(_kminmerSize > 3) return;

		vector<DbgEdge> removedEdges;

        for(const Unitig& unitig: _graph->_unitigs){

			//bool nodeName_ori;
			u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(unitig._endNode);

			vector<u_int32_t> successors;
			graph->getSuccessors(unitig._endNode, 0, successors);

			for(u_int32_t successor : successors){
				
				//bool nodeName_succ_ori;
				u_int32_t nodeName_succ = BiGraph::nodeIndex_to_nodeName(successor);

				if(Utils::shareAnyRead(_unitigDatas[nodeName], _unitigDatas[nodeName_succ])){
					continue;
				}

				Unitig& unitig_succ = graph->_unitigs[graph->nodeIndex_to_unitigIndex(successor)];

				//graph->_graphSuccessors->removeEdge(nodeName, nodeName_ori, nodeName_succ, nodeName_succ_ori);
				removedEdges.push_back({unitig._endNode, successor});
				
				
			}



		}


		cout << "nb edges: " << graph->_graphSuccessors->_nbEdges << endl;

		for(const DbgEdge& edge: removedEdges){

			bool edgeExist = graph->_graphSuccessors->removeEdge(edge._from, edge._to);
			//if(!edgeExist) cout << "nop 1" << endl;
			edgeExist = graph->_graphSuccessors->removeEdge(GraphSimplify::nodeIndex_toReverseDirection(edge._to), GraphSimplify::nodeIndex_toReverseDirection(edge._from));
			//if(!edgeExist) cout << "nop 2" << endl;
			
			//graph->_isEdgeUnsupported.insert(edge);
		}
		
		cout << "nb edges: " << graph->_graphSuccessors->_nbEdges << endl;
		
		cout << "Nb unsupported edges: " << graph->_isEdgeUnsupported.size() << endl;


	}



	size_t _iter;
	//ofstream file_groundTruth;
	//ofstream file_groundTruth_hifiasmContigs;

	ofstream file_kminmersContigs;
	


	MinimizerParser* _minimizerParser;
	u_int32_t _extract_truth_kminmers_read_position;
	unordered_map<u_int32_t, vector<string>> _evaluation_hifiasmGroundTruth_nodeName_to_unitigName;
	vector<u_int32_t> _evaluation_hifiasmGroundTruth_path;
	unordered_set<u_int32_t> _hifiasm_startingNodenames;
	ofstream _file_groundTruth_hifiasm_position;
	EncoderRLE _encoderRLE;

	void extract_truth_kminmers_read(const Read& read){
		
		u_int64_t readIndex = read._index;
		//ottalSize += strlen(read->seq.s);


		string rleSequence;
		vector<u_int64_t> rlePositions;
		_encoderRLE.execute(read._seq.c_str(), read._seq.size(), rleSequence, rlePositions);

		vector<u_int64_t> minimizers;
		vector<u_int64_t> minimizers_pos;
		_minimizerParser->parse(rleSequence, minimizers, minimizers_pos);

		vector<KmerVec> kminmers; 
		vector<ReadKminmer> kminmersInfo;
		MDBG::getKminmers(_minimizerSize, _kminmerSize, minimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, readIndex, false);

		for(size_t i=0; i<kminmers.size(); i++){

			KmerVec& vec = kminmers[i];
			//if(_mdbg->_dbg_nodes.find(kminmers[i]) == _mdbg->_dbg_nodes.end()) continue;

			//u_int32_t nodeName = _mdbg->_dbg_nodes[kminmers[i]]._index;
			if(_mdbg->_dbg_nodes.find(vec) != _mdbg->_dbg_nodes.end()){

				u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
				//2262408
				//if("utg009434l" == string(read->name.s)){
				//	cout << nodeName << endl;
				//	_hifiasm_startingNodenames.insert(nodeName);
				//}

				_evaluation_hifiasmGroundTruth_nodeName_to_unitigName[_mdbg->_dbg_nodes[vec]._index].push_back(read._header);
				_evaluation_hifiasmGroundTruth_path.push_back(_mdbg->_dbg_nodes[vec]._index);

				if(_evaluation_hifiasmGroundTruth_nodeNamePosition.find(_mdbg->_dbg_nodes[vec]._index) == _evaluation_hifiasmGroundTruth_nodeNamePosition.end()){
					_evaluation_hifiasmGroundTruth_nodeNamePosition[_mdbg->_dbg_nodes[vec]._index] = _extract_truth_kminmers_read_position;
				}
				//cout << mdbg->_dbg_nodes[vec]._index << " " << sequence.getComment() << endl;
				
				_file_groundTruth_hifiasm_position << nodeName << "," << _extract_truth_kminmers_read_position << endl;
			}
			else{
				//cout << "Not found position: " << position << endl;
				_evaluation_hifiasmGroundTruth_path.push_back(-1);
			}


			if(_evaluation_hifiasmGroundTruth.find(vec) != _evaluation_hifiasmGroundTruth.end()){
				//cout << position << endl;
				continue;
			}


			//_evaluation_hifiasmGroundTruth[vec] = datasetID;
			_evaluation_hifiasmGroundTruth_position[vec] = _extract_truth_kminmers_read_position;
			_extract_truth_kminmers_read_position += 1;

			
			//cout << _extract_truth_kminmers_read_position << endl;
		}


	}


	void extract_truth_kminmers(){

		_file_groundTruth_hifiasm_position.open(_inputDir + "/groundtruth_hifiasm_position.csv");
		_file_groundTruth_hifiasm_position << "Name,Position" << endl;

		_extract_truth_kminmers_read_position = 0;
		_minimizerParser = new MinimizerParser(_minimizerSize, _minimizerDensity);
		
		auto fp = std::bind(&GenerateContigs::extract_truth_kminmers_read, this, std::placeholders::_1);
		ReadParser readParser(_truthInputFilename, true, false);
		readParser.parse(fp);

		_file_groundTruth_hifiasm_position.close();

		delete _minimizerParser;
	}

	*/


	bool isContigAssembled(const vector<u_int32_t>& nodePath){

		//unordered_set<u_int32_t> distinctContigIndex;
		//u_int64_t nbBinnedNodes = 0;

		for(u_int32_t nodeIndex : nodePath){
			u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
			if(_processedNodeNames[nodeName]){
				return true;
				//distinctContigIndex.insert(processedNodeNames[nodeName]);
			}
			//else{
			//	return false;
			//}
		}

		return false;
		//float sharedRate = nbBinnedNodes / ((float) nodePath.size());

		//return sharedRate > 0.98;
		//return distinctContigIndex.size() == 1;
	}

/*
                ofstream outputFileContigColor(_outputDir + "/" + "contig_color.csv");
                outputFileContigColor << "Name,Color" << endl;

                for(const auto& it : unitigIndex_to_contigIndex){
                    u_int32_t unitigIndex = it.first;
                    const unordered_set<u_int32_t>& contigIndexes = it.second;

                    if(contigIndexes.size() > 1){
                        outputFileContigColor << unitigIndex << "," << "red" << endl;
                    }
                    else{
                        for(u_int32_t contigIndex : contigIndexes){
                            outputFileContigColor << unitigIndex << "," << "green" << endl; //contigIndex
                        }
                    }
                }

                outputFileContigColor.close();
				*/

/*
	void generateUnitigs(){
		//if(_kminmerSize < 10) return;

		ofstream outputFileContigColor(_inputDir + "/" + "contig_color.csv");
		outputFileContigColor << "Name,Color" << endl;

		const string& outputFilename = _inputDir + "/unitigs.nodepath.gz";

		cout << "Generating contigs" << endl;
		
		//const string& outputFilename = _inputDir + "/contigs.nodepath.gz";
		gzFile outputContigFile_min = gzopen(outputFilename.c_str(),"wb");

		unordered_set<u_int32_t> writtenUnitigs;

		//vector<float> readpathAbudance_values;
        for(const Unitig& u: _graph->_unitigs){

			
			//if(u._index == 880){
			//	cout << u._nodes.size() << endl;
			//	getchar();
 			//}



			//if(u._nbNodes < _kminmerSize*2) continue;

			if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._startNode)) != writtenUnitigs.end()) continue;
			if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._endNode)) != writtenUnitigs.end()) continue;

			writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._startNode));
			writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._endNode));

			vector<u_int32_t> nodepath = u._nodes;


			//if(u._index == 14) cout << u._nodes.size() << endl;
			//if(u._index == 64) cout << u._nodes.size() << endl;


			for(u_int32_t nodeIndex : u._nodes){
				u_int32_t unitigIndex = _graph->debug_nodeName_toSelectedUnitigIndex(BiGraph::nodeIndex_to_nodeName(nodeIndex));
				outputFileContigColor << unitigIndex << "," << "blue" << endl; //contigIndex
			}


				if(u._nbNodes < _kminmerSize) continue;
				
				if(u._nbNodes < _kminmerSize*2){

					double abundanceSum = 0;
					double abundanceN = 0;

					float minAbundance = std::numeric_limits<float>::max();

					vector<u_int32_t> successors;
					_graph->getSuccessors_unitig(u._index, 0, successors);
					vector<u_int32_t> predecessors;
					_graph->getPredecessors_unitig(u._index, 0, predecessors);

					for(u_int32_t unitigIndex : successors){
						abundanceSum += _graph->_unitigs[unitigIndex]._abundance * _graph->_unitigs[unitigIndex]._nbNodes;
						abundanceN += _graph->_unitigs[unitigIndex]._nbNodes;
						//if(_graph->_unitigs[unitigIndex]._abundance < minAbundance){
						//	minAbundance = _graph->_unitigs[unitigIndex]._abundance;
						//}
					}
					for(u_int32_t unitigIndex : predecessors){
						abundanceSum += _graph->_unitigs[unitigIndex]._abundance * _graph->_unitigs[unitigIndex]._nbNodes;
						abundanceN += _graph->_unitigs[unitigIndex]._nbNodes;
						//if(_graph->_unitigs[unitigIndex]._abundance < minAbundance){
						//	minAbundance = _graph->_unitigs[unitigIndex]._abundance;
						//}
					}

					if(abundanceN > 0){
						double mean = abundanceSum / abundanceN;
						if(u._abundance < mean*0.35){
							continue;
						}
					}
					//if(u._abundance < minAbundance){
					//	continue;
					//}
				}
				
				
			//}

			u_int64_t size = nodepath.size();

			//if(size < _kminmerSize*2) continue;
			for(u_int32_t nodeIndex : u._nodes){
				u_int32_t unitigIndex = _graph->debug_nodeName_toSelectedUnitigIndex(BiGraph::nodeIndex_to_nodeName(nodeIndex));
				outputFileContigColor << unitigIndex << "," << "red" << endl; //contigIndex
			}
			
			//cout << BiGraph::nodeIndex_to_nodeName(u._startNode) << " " << u._nbNodes << endl;
			gzwrite(outputContigFile_min, (const char*)&size, sizeof(size));
			gzwrite(outputContigFile_min, (const char*)&nodepath[0], size * sizeof(u_int32_t));
		}
		
	


		outputFileContigColor.close();
		gzclose(outputContigFile_min);
	}
*/
	/*

	
	struct Contig{
		vector<u_int32_t> _nodePath;
		vector<u_int32_t> _nodePath_sorted;
	};

	ofstream _fileTestLala;

	void extractContigKminmers (const string& outputFilename_fasta){
		
		_fileTestLala = ofstream(outputFilename_fasta + ".nodes.csv");
		_fileTestLala << "Name,Color" << endl;

		_mdbg = new MDBG(_kminmerSize);
		_mdbg->load(_inputDir + "/mdbg_nodes.gz");

		ReadParser parser(outputFilename_fasta, true, _minimizerSize, _kminmerSize, _minimizerDensity);
		auto fp = std::bind(&Assembly3::extractContigKminmers_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6);
		parser.parseKminmers(fp);

		for(auto& it: _contigFeature._nodeNameDuplicate){
			if(it.second == 1){
				_fileTestLala << it.first << "," << "blue" << endl;
			}
			else{
				_fileTestLala << it.first << "," << "red" << endl;
			}
		}

		delete _mdbg;
		_fileTestLala.close();
	}


	void extractContigKminmers_read(const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex, u_int64_t datasetIndex, const string& header, const string& seq){

		for(size_t i=0; i<kminmersInfos.size(); i++){
			

			KmerVec vec = kminmers[i];
			//for(KmerVec& vec : kminmers){
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()) continue;

			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
			_contigFeature._nodeNameDuplicate[nodeName] += 1;
			//_fileTestLala << nodeName << "," << "0" << endl;
			//_kminmerCounts[nodeName] = _countsInit;
			
		}

	}

	ofstream _file_contigToNode;

	void extractContigKminmers2 (const string& outputFilename_fasta){

		_file_contigToNode = ofstream(_inputDir + "/nodeToContig.csv");
		_file_contigToNode << "Name,Color" << endl;
		_contigFeature.loadAbundanceFile(_filename_abundance);

		ReadParser parser(outputFilename_fasta, true, _minimizerSize, _kminmerSize, _minimizerDensity);
		auto fp = std::bind(&Assembly3::extractContigKminmers_read2, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6);
		parser.parseKminmers(fp);

		_file_contigToNode.close();
		//vector<u_int32_t> nodePath = {933376, 1651014, 1762772, 732878};
		//vector<float> lala1;
		//vector<float> lala2;
		//_contigFeature.sequenceToAbundance(nodePath, lala1, lala2);

		//exit(1);
	}


	void extractContigKminmers_read2(const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex, u_int64_t datasetIndex, const string& header, const string& seq){

		if(seq.size() < 15000) return;
		//cout << seq.size() << endl;
		//if(readIndex == 20010) getchar(); 

		vector<float> composition;
		_contigFeature.sequenceToComposition(seq, composition);
		_contigFeature._contigCompositions[readIndex] = composition;

		unordered_set<u_int32_t> nodeNames;

		for(size_t i=0; i<kminmersInfos.size(); i++){
			

			KmerVec vec = kminmers[i];
			//for(KmerVec& vec : kminmers){
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()) continue;
			
			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
			//_file_contigToNode << nodeName << "," << readIndex << endl;
			//if(nodeName == 933376) cout << "HAHAHAHA" << endl;

			//if(_contigFeature._nodeName_to_contigIndex.find(nodeName) == _contigFeature._nodeName_to_contigIndex.end()) continue;

			_contigFeature._nodeNameDuplicate[nodeName] += 1;
			_contigFeature._nodeName_to_contigIndex[nodeName] = readIndex;
			//_fileTestLala << nodeName << "," << "0" << endl;
			//_kminmerCounts[nodeName] = _countsInit;
			
			nodeNames.insert(nodeName);
		}
		
		//std::sort(nodeNames.begin(), nodeNames.end());
		_contigFeature._contigNodes[readIndex] = nodeNames;

	}
	*/

	/*
	static bool ContigComparator_ByLength(const Contig &a, const Contig &b){
		return a._nodePath.size() > b._nodePath.size();
	}

	static bool ContigComparator_ByLength_Reverse(const Contig &a, const Contig &b){
		return a._nodePath.size() > b._nodePath.size();
	}


	void dereplicateContig2(vector<Contig>& contigs, vector<u_int32_t> nodePath, unordered_set<u_int32_t>& processedNodeNames, unordered_map<u_int32_t, u_int32_t>& nodeName_to_contigIndex, u_int64_t& contigIndex){


		if(nodePath.size() < _kminmerSize*2) return;

		vector<u_int32_t> component;
		for(u_int32_t nodeIndex : nodePath){
			component.push_back(BiGraph::nodeIndex_to_nodeName(nodeIndex));
		}
		//for(size_t i=_kminmerSize-1; i<nodePath.size()-(_kminmerSize-1); i++){
		//	component.push_back(BiGraph::nodeIndex_to_nodeName(nodePath[i]));
		//}
		std::sort(component.begin(), component.end());

		bool isValid = true;

		std::sort(contigs.begin(), contigs.end(), ContigComparator_ByLength);

		for(Contig& contig : contigs){
			if(contig._nodePath.size() == 0) continue;

			u_int64_t nbShared = Utils::computeSharedElements(component, contig._nodePath_sorted);
			if(nbShared == 0) continue;

			if(contig._nodePath.size() > component.size()){
				return;
			}
			else{
				contig._nodePath.clear();
				contig._nodePath_sorted.clear();
			}
		}
		
		addContig(nodePath, contigs, processedNodeNames, nodeName_to_contigIndex, contigIndex);
	}

	void addContig(vector<u_int32_t>& nodePath, vector<Contig>& contigs, unordered_set<u_int32_t>& processedNodeNames, unordered_map<u_int32_t, u_int32_t>& nodeName_to_contigIndex, u_int64_t& contigIndex){

		if(nodePath.size() <= _kminmerSize*2) return;

		vector<u_int32_t> nodePath_sorted;
		for(u_int32_t nodeIndex : nodePath){
			nodePath_sorted.push_back(BiGraph::nodeIndex_to_nodeName(nodeIndex));
		}

		//for(size_t i=_kminmerSize-1; i<nodePath.size()-(_kminmerSize-1); i++){
		//	nodePath_sorted.push_back(BiGraph::nodeIndex_to_nodeName(nodePath[i]));
		//}

		std::sort(nodePath_sorted.begin(), nodePath_sorted.end());
		
		contigs.push_back({nodePath, nodePath_sorted});



		//bool isLala = false;
		for(u_int32_t nodeIndex : nodePath){
			u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
			nodeName_to_contigIndex[nodeName] = contigIndex;
			//if(nodeName == 289994){
			//	isLala = true;
			//}
			processedNodeNames.insert(nodeName);
		}

		contigIndex += 1;




	}
	*/


	struct Contig{
		u_int64_t _readIndex;
		vector<u_int32_t> _nodepath;
	};

	//static bool UnitigComparator_ByLength2(const UnitigLength &a, const UnitigLength &b){
	//	return a._length > b._length;
	//}

	float _minUnitigAbundance;

	//struct Contig{
	//	vector<u_int32_t> _nodepath;
	//	vector<u_int32_t> _nodepath_sorted;
	//};

	static bool ContigComparator_ByLength(const Contig &a, const Contig &b){
		return a._nodepath.size() > b._nodepath.size();
	}
	
	static bool ContigComparator_ByLength2(const Contig &a, const Contig &b){

		if(a._nodepath.size() == b._nodepath.size()){
			for(size_t i=0; i<a._nodepath.size() && i<b._nodepath.size(); i++){
				if(a._nodepath[i] == b._nodepath[i]){
					continue;
				}
				else{
					return a._nodepath[i] > b._nodepath[i];
				}
			}
		}


		return a._nodepath.size() > b._nodepath.size();
	}

	struct NodeAb{
		float _abundance;
		u_int32_t _nbNodes;
	};

	vector<bool> _processedNodeNames;
	unordered_map<u_int32_t, NodeAb> _nodeNameAbundances;


	
	/*
	void generateContigPathFile(){

		_graph->saveUnitigGraph(_inputDir + "/minimizer_graph_u_cleaned.gfa", nullptr, _minimizerSize, _nbCores, true);

		_logFile << "Generating contig path file" << endl;
		ifstream contigFile(_inputDir + "/contigs.nodepath");
		ofstream outputFile(_inputDir + "/contigs_path.csv");

		//outputFile << "Name,Color" << endl;

		u_int64_t contigIndex = 0;
		vector<u_int32_t> unitigPath;

		while(true){

			vector<u_int32_t> nodePath;
			u_int64_t size;
			contigFile.read((char*)&size, sizeof(size));
			

			if(contigFile.eof()) break;

			u_int8_t isCircular;
			contigFile.read((char*)&isCircular, sizeof(isCircular));

			nodePath.resize(size);
			contigFile.read((char*)&nodePath[0], size * sizeof(u_int32_t));

			string contigName = "ctg" + to_string(contigIndex);
			outputFile << contigName;

            u_int32_t prevUnitigIndex = -1;
			
			for(size_t i=0; i<nodePath.size(); i++){
				u_int32_t nodeIndex = nodePath[i];
				bool orientation;
				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex, orientation);

                u_int32_t unitigIndex = _graph->nodeName_toSelectedUnitigIndex(nodeName, _graph->_selectedUnitigIndexTmp);

                if(unitigIndex != prevUnitigIndex){
					unitigPath.push_back(unitigIndex);
				    //nodeNames.push_back(unitigIndex);
                    prevUnitigIndex = unitigIndex;

					//outputFile << unitigIndex << "," << contigName << endl;
					outputFile << ";" << unitigIndex;
                }
			}

			//for(u_int32_t unitigIndex : unitigPath){
			//	outputFile << ";" << unitigIndex;
			//}

			outputFile << endl;

			contigIndex += 1;
		}

		contigFile.close();
		outputFile.close();

		//_logFile << _nodeName_entire.size() << " " << _nodeName_left.size() << " " << _nodeName_right.size() << endl;
		//_logFile << "Nb contigs: " << nbContigs << endl;
		
	}
	*/
/*
Nb contigs: 384
Nb nodes (checksum): 258228925

Nb contigs: 236
Nb nodes (checksum): 165166644
*/
	/*
	void getSupportingReads(const vector<u_int32_t>& pathNodes, vector<u_int64_t>& supportingReads){

		//cout << "lala" << endl;
		supportingReads.resize(pathNodes.size());
		return;

		supportingReads.clear();
		vector<u_int32_t> prevNodes;

		for(u_int32_t nodeIndex : pathNodes){

			//cout << nodeIndex << " " << prevNodes.size() <<  endl;
			u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
			UnitigData& unitigData = _unitigDatas[nodeName];

			if(prevNodes.size() == 0){
				
				cout << unitigData._readIndexes.size() << endl;

				prevNodes.push_back(nodeIndex);
				supportingReads.push_back(unitigData._readIndexes[0]);
				continue;
			}

			u_int32_t prevRank = 0;
			u_int64_t supportingRead = -1;

			while(true){
				
				int prevIndex = prevNodes.size() - prevRank - 1;
				//cout << "lalalala     " << prevIndex << endl;
				if(prevIndex < 0 ) break;
				
				u_int32_t prev_nodeIndex = prevNodes[prevIndex];
				u_int32_t prev_nodeName = BiGraph::nodeIndex_to_nodeName(prev_nodeIndex);
				UnitigData& prev_unitigData = _unitigDatas[prev_nodeName];
				
				vector<u_int64_t> sharedReads;
				Utils::collectSharedReads(unitigData, prev_unitigData, sharedReads);
				//cout << "sdfsdfsd     " << sharedReads.size() << endl;
				if(sharedReads.size() > 0){
					prevRank += 1;
					supportingRead = sharedReads[0];

					//cout << "lala: " << supportingRead << endl;
				}
				else{
					break;
				}
			}


			prevNodes.push_back(nodeIndex);
			supportingReads.push_back(supportingRead);

		}
		
		cout << "loulou" << endl;
	}
	*/




	/*
	void generateContigs(){
		if(_kminmerSize < 10) return;

		cout << "Generating contigs" << endl;
		
		const string& outputFilename = _inputDir + "/contigs.nodepath.gz";
		gzFile outputContigFile_min = gzopen(outputFilename.c_str(),"wb");

		unordered_set<u_int32_t> writtenUnitigs;

		vector<float> readpathAbudance_values;
		for(const Unitig& u : _graph->_unitigs){
			//if(u._nbNodes < _kminmerSize*2) continue;

			if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._startNode)) != writtenUnitigs.end()) continue;
			if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._endNode)) != writtenUnitigs.end()) continue;

			writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._startNode));
			writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._endNode));

			u_int64_t size = u._nodes.size();
			gzwrite(outputContigFile_min, (const char*)&size, sizeof(size));
			gzwrite(outputContigFile_min, (const char*)&u._nodes[0], size * sizeof(u_int32_t));
		}
		
		gzclose(outputContigFile_min);
	}
	*/

	void generateContigs3(){

		Logger::get().debug() << "Generating contigs";
		size_t nextMultikStep = Commons::getMultikStep(_kminmerSize);


		ofstream outputContigFile(_inputDir + "/contigs.nodepath");

		float prevCutoff = -1;
		u_int64_t contigIndex = 0;


		u_int64_t checksumTotal = 0;
		//unordered_set<u_int32_t> processedNodeNames;

		unordered_set<u_int32_t> actualAbundantNodes;
		unordered_set<u_int32_t> seenAbundantNodes;

		for(long i=_progressiveAbundanceFilter->_cutoffIndexes.size()-1; i>=0; i--){

			ProgressiveAbundanceFilter::CutoffIndex cutoffIndex = _progressiveAbundanceFilter->_cutoffIndexes[i];
			float cutoff = cutoffIndex._cutoffValue;
			Logger::get().debug() << cutoff;

        	ifstream nodepathFile(_inputDir + "/filter/unitigs_" + to_string(cutoffIndex._cutoffIndex) + ".bin");

			_minUnitigAbundance = cutoff / 0.5;

			while(true){

				u_int8_t isCircular = CONTIG_LINEAR;
				bool isRepeatSide = false;
				float contigAbundance = -1;
				vector<u_int32_t> nodePath;
				u_int32_t size;
				nodepathFile.read((char*)&size, sizeof(size));
				
				if(nodepathFile.eof()) break;


				nodepathFile.read((char*)&isCircular, sizeof(isCircular));
				nodepathFile.read((char*)&isRepeatSide, sizeof(isRepeatSide));
				nodepathFile.read((char*)&contigAbundance, sizeof(contigAbundance));

				nodePath.resize(size);
				nodepathFile.read((char*)&nodePath[0], size * sizeof(u_int32_t));

				if(isRepeatSide) continue;
				if(contigAbundance < _minUnitigAbundance) continue;
				if(isContigAssembled(nodePath)) continue;

				u_int64_t checksum = 0;
				for(u_int32_t nodeIndex : nodePath){
					checksum += BiGraph::nodeIndex_to_nodeName(nodeIndex);
				}
				checksum *= nodePath.size();
				checksumTotal += checksum;


				if(isCircular && nodePath.size() > 1){
					nodePath.push_back(nodePath[0]);
				}

				size = nodePath.size();
				//cout << size << "  " << nodePath.size() << " " << contigAbundance << endl;
				outputContigFile.write((const char*)&size, sizeof(size));
				outputContigFile.write((const char*)&isCircular, sizeof(isCircular));
				outputContigFile.write((const char*)&nodePath[0], size * sizeof(u_int32_t));

				//if(size > 100)
				//	cout << "contig: " << size << endl;

				
				//u_int32_t nodeName1 = BiGraph::nodeIndex_to_nodeName(nodePath[0]);
				//u_int32_t nodeName2 = BiGraph::nodeIndex_to_nodeName(nodePath[nodePath.size()-1]);

				for(u_int32_t nodeIndex : nodePath){
					u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
					//cout << nodeName << endl;
					_processedNodeNames[nodeName] = true;
					_nodeNameAbundances[nodeName] = {contigAbundance, (u_int32_t)nodePath.size()};

					/*
					if(_abundantNodeNames.find(nodeName) != _abundantNodeNames.end()){
						seenAbundantNodes.insert(nodeName);

						if(nodeName1 != nodeName && nodeName2 != nodeName){
						}
						else{
							cout << "omg: " << contigAbundance << endl;
						}
					}
					*/
				}

				/*
				if(_abundantNodeNames.find(nodeName1) != _abundantNodeNames.end()){
					actualAbundantNodes.insert(nodeName1);
				}
				if(_abundantNodeNames.find(nodeName2) != _abundantNodeNames.end()){
					actualAbundantNodes.insert(nodeName2);
				}

				cout << actualAbundantNodes.size() << " / " << seenAbundantNodes.size() << " / " << _abundantNodeNames.size() << endl;
				*/

				contigIndex += 1;
			}

			nodepathFile.close();
		}


		outputContigFile.close();
				//if(_kminmerSize > 18){
				//	getchar();
				//}
		//cout << "Checksum: " << checksumTotal << endl;
		//getchar();
	}


	
};





#endif
