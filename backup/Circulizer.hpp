

#ifndef MDBG_METAG_CIRCULIZER
#define MDBG_METAG_CIRCULIZER

#include "../Commons.hpp"

class Circulizer : public Tool{
    
public:

	string _mdbgDir;
	string _tmpDir;
	string _referenceFilename;
	string _outputFilename;
	size_t _nbCores;

	MDBG* _mdbg;
	ofstream _outputFile;
	GraphSimplify* _graph;
	string _filename_readMinimizers;

	float _minimizerDensity;
    size_t _minimizerSize;
    size_t _kminmerSize;
    size_t _kminmerSizeFirst;
	float _minimizerSpacingMean;
	float _kminmerLengthMean;
	float _kminmerOverlapMean;
	vector<UnitigData> _unitigDatas;
	vector<vector<u_int32_t>> _nodeName_to_readIndexes;


	Circulizer(): Tool (){

	}


	void parseArgs(int argc, char* argv[]){

		/*
		//_kminmerSize = 4;
		//_minimizerSize = 21;
		//_minimizerDensity = 0.05;

		cxxopts::Options options("ToBasespace", "");
		options.add_options()
		("mdbgDir", "", cxxopts::value<string>())
		("binDir", "", cxxopts::value<string>())
		("outputFilename", "", cxxopts::value<string>());

		options.parse_positional({"mdbgDir", "binDir", "outputFilename"});
		options.positional_help("mdbgDir binDir outputFilename");



		if(argc <= 1){
			cout << options.help() << endl;
			exit(0);
		}

		cxxopts::ParseResult result;

		try{
			result = options.parse(argc, argv);

			_mdbgDir = result["mdbgDir"].as<string>();
			_binDir = result["binDir"].as<string>();
			_outputFilename = result["outputFilename"].as<string>();
		}
		catch (const std::exception& e){
			std::cout << options.help() << std::endl;
			std::cerr << e.what() << std::endl;
			std::exit(EXIT_FAILURE);
		}
		*/
		args::ArgumentParser parser("circ", ""); //"This is a test program.", "This goes after the options."
		args::Positional<std::string> arg_tmpDir(parser, "tmpDir", "", args::Options::Required);
		args::Positional<std::string> arg_mdbgDir(parser, "mdbgDir", "", args::Options::Required);
		//args::Positional<std::string> arg_binDir(parser, "referenceFilename", "", args::Options::Required);
		//args::Positional<std::string> arg_outputFilename(parser, "outputFilename", "", args::Options::Required);
		args::ValueFlag<int> arg_nbCores(parser, "", "Number of cores", {ARG_NB_CORES2}, NB_CORES_DEFAULT_INT);
		//args::Positional<std::string> arg_contigs(parser, "contigs", "", args::Options::Required);
		//args::PositionalList<std::string> arg_readFilenames(parser, "reads", "Input filename(s) (separated by space)", args::Options::Required);
		//args::ValueFlag<int> arg_l(parser, "", "Minimizer length", {ARG_MINIMIZER_LENGTH2}, 13);
		//args::ValueFlag<float> arg_d(parser, "", "Minimizer density", {ARG_MINIMIZER_DENSITY2}, 0.005f);
		//args::ValueFlag<std::string> arg_contigs(parser, "", "", {ARG_INPUT_FILENAME_CONTIG}, "");
		//args::ValueFlag<int> arg_nbCores(parser, "", "Number of cores", {ARG_NB_CORES2}, NB_CORES_DEFAULT_INT);
		//args::Flag arg_cutInternal(parser, "", "", {ARG_CUT_INTERNAL});
		//args::Flag arg_noDump(parser, "", "", {ARG_NO_DUMP});
		//args::Flag arg_isFinalAssembly(parser, "", "Is final multi-k pass", {ARG_FINAL});
		//args::Flag arg_firstPass(parser, "", "Is first pass of multi-k", {ARG_FIRST_PASS});
		args::Flag arg_help(parser, "", "", {'h', "help"}, args::Options::Hidden);
		//args::HelpFlag help(parser, "help", "Display this help menu", {'h'});
		//args::CompletionFlag completion(parser, {"complete"});

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

		_tmpDir = args::get(arg_tmpDir);
		_mdbgDir = args::get(arg_mdbgDir);
		_nbCores = args::get(arg_nbCores);
		//_referenceFilename = args::get(arg_binDir);
		//_outputFilename = args::get(arg_outputFilename);



		string filename_parameters = _mdbgDir + "/parameters.gz";
		gzFile file_parameters = gzopen(filename_parameters.c_str(),"rb");
		gzread(file_parameters, (char*)&_minimizerSize, sizeof(_minimizerSize));
		gzread(file_parameters, (char*)&_kminmerSize, sizeof(_kminmerSize));
		gzread(file_parameters, (char*)&_minimizerDensity, sizeof(_minimizerDensity));
		gzread(file_parameters, (char*)&_kminmerSizeFirst, sizeof(_kminmerSizeFirst));
		gzread(file_parameters, (char*)&_minimizerSpacingMean, sizeof(_minimizerSpacingMean));
		gzread(file_parameters, (char*)&_kminmerLengthMean, sizeof(_kminmerLengthMean));
		gzread(file_parameters, (char*)&_kminmerOverlapMean, sizeof(_kminmerOverlapMean));
		gzclose(file_parameters);

		openLogFile(_tmpDir);

		//ifstream file_data(_mdbgDir + "/data.txt");
		//file_data.read((char*)&_nbReads, sizeof(_nbReads));
		//file_data.close();

		cout << endl;
		cout << "Input dir: " << _mdbgDir << endl;
		//cout << "Output filename: " << _outputFilename << endl;
		cout << "Minimizer length: " << _minimizerSize << endl;
		cout << "Kminmer length: " << _kminmerSize << endl;
		cout << "Density: " << _minimizerDensity << endl;
		//cout << "Nb reads: " << _nbReads << endl;
		cout << endl;

		//_kminmerSizeFirst = 4;

		_filename_readMinimizers = _tmpDir + "/read_data.txt";
		//_minimizerParser = new MinimizerParser(_minimizerSize, _minimizerDensity);
	}


    void execute (){



		loadGraph();
		indexReads();
		generateContigs2(_tmpDir + "/contigs.nodepath");

		closeLogFile();
		_outputContigFile.close();
	}


	void loadGraph(){






		string _gfaFilename = _mdbgDir + "/minimizer_graph.gfa";
		//string gfa_filename_noUnsupportedEdges = _inputDir + "/minimizer_graph_noUnsupportedEdges.gfa";
		//string gfa_filename_unitigs = _inputDir + "/minimizer_graph_unitigs.gfa";
		//string mdbg_filename = _inputDir + "/mdbg_nodes.gz";




		
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
		
		GraphSimplify* graphSimplify = new GraphSimplify(_gfaFilename, _mdbgDir, 0, _kminmerSize, _nbCores, _kminmerLengthMean, _kminmerOverlapMean, _logFile);
		_graph = graphSimplify;
		

		ifstream kminmerFile(_mdbgDir + "/kminmerData_min.txt");
        _graph->_graphSuccessors->_nodeDatas.resize(_graph->_graphSuccessors->_nbNodes/2, {0, 0, 0});

		u_int64_t abSum = 0;
		u_int64_t qualSum = 0;

		while (true) {

			vector<u_int64_t> minimizerSeq;
			minimizerSeq.resize(_kminmerSize);
			kminmerFile.read((char*)&minimizerSeq[0], minimizerSeq.size()*sizeof(u_int64_t));

			if(kminmerFile.eof()) break;

			u_int32_t nodeName;
			u_int32_t abundance;
			u_int32_t quality;
			//bool isReversed = false;

			kminmerFile.read((char*)&nodeName, sizeof(nodeName));
			kminmerFile.read((char*)&abundance, sizeof(abundance));
			kminmerFile.read((char*)&quality, sizeof(quality));

			abSum += abundance;
			qualSum += quality;

			//cout << nodeName << " " << abSum << endl;

			//cout << nodeName << " " << abundance << " " << quality << endl;
			//if(quality == 0){
			//	cout << quality << endl;
			//	getchar();
			//}
			//if(abundance == 1) continue;
			//KmerVec vec;
			//vec._kmers = minimizerSeq;

			_graph->_graphSuccessors->_nodeDatas[nodeName] = {abundance, (u_int32_t)_kminmerLengthMean, quality};
			//_graph->_graphSuccessors->_nodeLengths[nodeName] = _kminmerLengthMean;
			//_kminmerAbundances[vec] = abundance;
		}

		//cout << abSum << " " << qualSum << endl;
		//getchar();

		kminmerFile.close();

		//GfaParser::binaryGraph_to_gfa(_gfaFilename, _kminmerLengthMean, _kminmerOverlapMean,  _gfaFilename + "_debug.gfa", _graph->_graphSuccessors->_nodeDatas);
		
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

	  //_graph->debug_writeGfaErrorfree(0, 0, -1, _kminmerSize, false, true, false, _unitigDatas, true, false, true, false, true, true, _mdbg, _minimizerSize, _nbCores, true, false);
		//_graph->debug_writeGfaErrorfree(0, 0, -1, _kminmerSize, false, true, false, _unitigDatas, true, false, false, false, false, true, _mdbg, _minimizerSize, _nbCores, false, false);
		_graph->debug_writeGfaErrorfree(2000, 2000, -1, _kminmerSize, false, true, false, _unitigDatas, true, false, false, false, false, true, _mdbg, _minimizerSize, _nbCores, false, false);
		
		
		//_graph->debug_selectUnitigIndex();
		//_graph->debug_writeGfaErrorfree(2000, 2000, -1, _kminmerSize, false, true, false, _unitigDatas, true, false, false, false, false, false, _mdbg, _minimizerSize, _nbCores, false, true);
			
		//delete _mdbg;
		
		//cout << "to remove graph binary to gfa" << endl;
    	//GfaParser::binaryGraph_to_gfa(_gfaFilename, _kminmerLengthMean, _kminmerOverlapMean, _gfaFilename+".gfa", _graph->_graphSuccessors->_nodeDatas);


		//"reprise: extend from multiple nodes"
	}

	void indexReads(){

		_nodeName_to_readIndexes.resize(_graph->_graphSuccessors->_nbNodes/2);

		_logFile << "Indexing reads" << endl;
		string mdbg_filename = _mdbgDir + "/kminmerData_min.txt";
		_mdbg = new MDBG(_kminmerSize);
		_mdbg->load(mdbg_filename, false);

		KminmerParserParallel parser(_filename_readMinimizers, _minimizerSize, _kminmerSize, false, false, _nbCores);
		parser.parse(IndexReadsFunctor(*this));

		//delete _mdbg;
	}


	class IndexReadsFunctor {

		public:

		Circulizer& _generateContigs;

		IndexReadsFunctor(Circulizer& generateContigs) : _generateContigs(generateContigs){
		}

		IndexReadsFunctor(const IndexReadsFunctor& copy) : _generateContigs(copy._generateContigs){
		}

		~IndexReadsFunctor(){
		}

		//void collectBestSupportingReads_read(const vector<u_int64_t>& readMinimizers, const vector<ReadKminmerComplete>& kminmersInfos, bool isCircular, u_int64_t readIndex){

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

				#pragma omp critical(indexReadsFunctor)
				{
					_generateContigs._nodeName_to_readIndexes[nodeName].push_back(readIndex);
				}

			}

		}
	};

	vector<bool> _processedNodeNames;
	float _minUnitigAbundance;
	//unordered_map<u_int32_t, NodeAb> _nodeNameAbundances;

	ofstream _outputContigFile;

	void generateContigs2(const string& outputFilename){

		_outputContigFile = ofstream(outputFilename);


		_processedNodeNames = vector<bool>(_graph->_isNodeRemoved.size()/2, false);
		u_int64_t nbNodesCheckSum = 0;
		//string clusterDir = _inputDir + "/" + "binGreedy";
		//fs::path path(clusterDir);
		//if(!fs::exists (path)){
		//	fs::create_directory(path);
		//} 

		//string filename_binStats = clusterDir + "/binStats.txt";
		//ofstream file_binStats(filename_binStats);

		//ofstream fileHifiasmAll(_inputDir + "/binning_results_hifiasm.csv");
		//fileHifiasmAll << "Name,Colour" << endl;

		//gzFile outputContigFile_min = gzopen(outputFilename.c_str(),"wb");
		//gzFile outputContigFile_fasta = gzopen(outputFilename_fasta.c_str(),"wb");

		//ofstream file_asmResult = ofstream(_inputDir + "/binning_results.csv");
		//file_asmResult << "Name,Colour" << endl;

		//Assembly assembly(_unitigDatas, _contigFeature);

		float prevCutoff = -1;
		u_int64_t contigIndex = 0;


		//unordered_set<u_int32_t> processedNodeNames;
		u_int64_t processedUnitigs = 0;

		//u_int64_t contigIndex = 0;
		u_int64_t __loadState2_index = 0;


		
		vector<float> allCutoffs;

		for(const SaveState2& saveState : _graph->_cachedGraphStates){
			allCutoffs.push_back(saveState._abundanceCutoff_min);
			//cout << saveState._abundanceCutoff_min << endl;
		}
		
		std::reverse(allCutoffs.begin(), allCutoffs.end());




        //file_debug = ofstream("/home/gats/workspace/run//debug_graph.csv");
		//file_debug << "Name,Colour" << endl;

        //auto t1 = std::chrono::high_resolution_clock::now();
		
        _graph->clear(0);

        for(const SaveState2& saveState : _graph->_cachedGraphStates){
            //if(saveState._abundanceCutoff_min > abundanceCutoff_min) break;

            for(const u_int32_t nodeName : saveState._nodeNameRemoved){

                u_int32_t nodeIndex1 = BiGraph::nodeName_to_nodeIndex(nodeName, true);
                u_int32_t nodeIndex2 = BiGraph::nodeName_to_nodeIndex(nodeName, false);
				_graph->_isNodeRemoved[nodeIndex1] = true;
				_graph->_isNodeRemoved[nodeIndex2] = true;
                //_graph->_isNodeValid2.erase(nodeIndex1);
                //_graph->_isNodeValid2.erase(nodeIndex2);
            }

            for(const DbgEdge& dbgEdge : saveState._isEdgeRemoved){
				_graph->_graphSuccessors->setEdgeRemoved(dbgEdge._from, dbgEdge._to, true);
				_graph->_graphSuccessors->setEdgeRemoved(GraphSimplify::nodeIndex_toReverseDirection(dbgEdge._to), GraphSimplify::nodeIndex_toReverseDirection(dbgEdge._from), true);
                //_graph->_isEdgeRemoved.insert(dbgEdge);
            }
        }



		u_int64_t contigIndexLala = 0;
		//unordered_map<u_int32_t, u_int32_t> nodeName_to_contigIndex;
		//vector<Contig> contigs;
		u_int64_t checksum_global = 0;
		double checksum_abundance = 0;
		u_int64_t checksum_nbNodes = 0;
		u_int64_t checkSum = 0;

		for(size_t i=0; i<allCutoffs.size(); i++){

			float cutoff = allCutoffs[i];

			if(i > 0){

				//cout << i << endl;
				//if(cutoff == 0) continue;
				const vector<u_int32_t>& nodeNameRemoved = _graph->_cachedGraphStates[_graph->_cachedGraphStates.size()-1-i+1]._nodeNameRemoved;
				const vector<DbgEdge>& isEdgeRemoved = _graph->_cachedGraphStates[_graph->_cachedGraphStates.size()-1-i+1]._isEdgeRemoved;

				//cout << nodeNameRemoved.size() << endl;
				//cout << isEdgeRemoved.size() << endl;
				//cout << "Unload state: " << _graph->_cachedGraphStates[_graph->_cachedGraphStates.size()-1-i+1]._abundanceCutoff_min << endl;

				for(const u_int32_t nodeName : nodeNameRemoved){

					u_int32_t nodeIndex1 = BiGraph::nodeName_to_nodeIndex(nodeName, true);
					u_int32_t nodeIndex2 = BiGraph::nodeName_to_nodeIndex(nodeName, false);
					_graph->_isNodeRemoved[nodeIndex1] = false;
					_graph->_isNodeRemoved[nodeIndex2] = false;
					//_graph->_isNodeValid2.insert(nodeIndex1);
					//_graph->_isNodeValid2.insert(nodeIndex2);
				}


				for(const DbgEdge& dbgEdge : isEdgeRemoved){
					_graph->_graphSuccessors->setEdgeRemoved(dbgEdge._from, dbgEdge._to, false);
					_graph->_graphSuccessors->setEdgeRemoved(GraphSimplify::nodeIndex_toReverseDirection(dbgEdge._to), GraphSimplify::nodeIndex_toReverseDirection(dbgEdge._from), false);
					//_graph->_isEdgeRemoved.erase(dbgEdge);
				}
				//cout << _graph->_isEdgeRemoved.size() << endl;
				//cout << _graph->_isEdgeRemoved.size() << endl;

			}


        	_graph->compact(false, _unitigDatas);
			//cout << "Cutoff: " << cutoff << " " << _graph->getChecksumGlobal() << endl;
			checksum_global += _graph->getChecksumGlobal_utg();
			checksum_abundance += _graph->getChecksumGlobal_abundanceUtg();


			_minUnitigAbundance = cutoff / 0.5;

			cout << "Cutoff: " << cutoff << endl;

			const string& filename_contigs = _tmpDir + "/contig_data.txt";
			KminmerParserParallel parser3(filename_contigs, _minimizerSize, _kminmerSize, false, false, _nbCores);
			parser3.parse(ContigFunctor(*this));
			/*
			for(size_t i=0; i<_graph->_unitigs.size(); i+=2){

				const Unitig& u1 = _graph->_unitigs[i];
				const Unitig& u2 = _graph->_unitigs[i+1];

				Unitig u;

				if(u1._startNode < u2._startNode){
					u = u1;
				}
				else{
					u = u2;
				}
				//cout << unitig._length << " " << unitig._abundance << endl;
				//if(unitig._index % 2 == 1) continue;

				if(u._abundance < _minUnitigAbundance) continue;
				if(isContigAssembled(u._nodes)) continue;

				processedUnitigs += 1;

				vector<u_int32_t> nodePath = u._nodes;


                vector<u_int32_t> successors;
                _graph->getSuccessors_unitig(u._index, 0, successors);

                vector<u_int32_t> predecessors;
                _graph->getPredecessors_unitig(u._index, 0, predecessors);

				if(successors.size() == 0 && predecessors.size() == 0 && u._abundance == 1) continue;


				const string& filename_contigs = _tmpDir + "/unitig_data.txt";
				KminmerParserParallel parser3(filename_contigs, _minimizerSize, _kminmerSize, false, false, _nbCores);
				parser3.parse(ContigFunctor(*this));


				bool isCircular = u._startNode == u._endNode;



				u_int64_t size = nodePath.size();

				//gzwrite(outputContigFile_min, (const char*)&size, sizeof(size));
				//gzwrite(outputContigFile_min, (const char*)&nodePath[0], size * sizeof(u_int32_t));
				outputContigFile.write((const char*)&size, sizeof(size));
				outputContigFile.write((const char*)&isCircular, sizeof(isCircular));
				outputContigFile.write((const char*)&nodePath[0], size * sizeof(u_int32_t));

				for(u_int32_t nodeIndex : nodePath){
					u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
					if(u._startNode != u._endNode) nbNodesCheckSum += nodeName;
					checksum_nbNodes += 1;
					_processedNodeNames[nodeName] = true;

				}



				contigIndex += 1;
			}
			*/





		}

	}

	bool isContigAssembled(const vector<u_int32_t>& nodePath){

		for(u_int32_t nodeIndex : nodePath){
			u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
			if(_processedNodeNames[nodeName]){
				return true;
			}
		}

		return false;
	}

	void writeContig(const vector<u_int32_t>& nodePath){
		
		bool isCircular = nodePath[0] == nodePath[nodePath.size()-1];



		u_int64_t size = nodePath.size();

		//gzwrite(outputContigFile_min, (const char*)&size, sizeof(size));
		//gzwrite(outputContigFile_min, (const char*)&nodePath[0], size * sizeof(u_int32_t));
		_outputContigFile.write((const char*)&size, sizeof(size));
		_outputContigFile.write((const char*)&isCircular, sizeof(isCircular));
		_outputContigFile.write((const char*)&nodePath[0], size * sizeof(u_int32_t));
	}

	class ContigFunctor {

		public:

		Circulizer& _circulizer;

		ContigFunctor(Circulizer& circulizer) : _circulizer(circulizer){
		}

		ContigFunctor(const ContigFunctor& copy) : _circulizer(copy._circulizer){
		}

		~ContigFunctor(){
		}
		
		bool allContigNodeInGraph(const vector<u_int32_t>& nodePath){

			
			/*
			//return !_circulizer._graph->_isNodeRemoved[nodePath[0]] && !_circulizer._graph->_isNodeRemoved[nodePath[nodePath.size()-1]];
			for(u_int32_t nodeIndex : nodePath){
				if(_circulizer._graph->_isNodeRemoved[nodeIndex]){
					return false;
				}
				//u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
				//if(_processedNodeNames[nodeName]){
				//	return true;
				//}
			}

			return true;
			*/
		}

		void operator () (const KminmerList& kminmerList) {



			ofstream colorFile(_circulizer._mdbgDir + "/color.csv");
			colorFile << "Name,Color" << endl;

			u_int64_t readIndex = kminmerList._readIndex;
			const vector<u_int64_t>& readMinimizers = kminmerList._readMinimizers;
			//const vector<KmerVec>& kminmers = kminmerList._kminmers;
			const vector<ReadKminmerComplete>& kminmersInfos = kminmerList._kminmersInfo;

			if(kminmerList._isCircular) return;

			//float minNbMinimizers = 1000000 / _kminmerOverlapMean;
			//if(kminmerList._kminmersInfo.size() > )

			vector<u_int32_t> nodepath;
			vector<float> abundances;

			cout << kminmersInfos.size() << endl;

			for(size_t i=0; i<kminmersInfos.size(); i++){
				

				//const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];

				const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];
				const KmerVec& vec = kminmerInfo._vec;

				if(_circulizer._mdbg->_dbg_nodes.find(vec) == _circulizer._mdbg->_dbg_nodes.end()) continue;

				u_int32_t nodeName = _circulizer._mdbg->_dbg_nodes[vec]._index;
				u_int32_t nodeIndex = BiGraph::nodeName_to_nodeIndex(nodeName, !kminmerInfo._isReversed);

				nodepath.push_back(nodeIndex);

				if(!_circulizer._graph->_isNodeRemoved[nodeIndex]){
					abundances.push_back(_circulizer._graph->nodeIndex_to_unitig(nodeIndex)._abundance);
				}

				//if(i == 0 || i == kminmersInfos.size()-1){
				//	colorFile << nodeName << "," << "green" << endl;
				//}
				//else{
				//	colorFile << nodeName << "," << "red" << endl;
				//}
			}

			float medianAbundance = Utils::compute_median_float(abundances);
			if(medianAbundance < 30) return;

			cout << medianAbundance << " " << _circulizer._minUnitigAbundance << endl;

			if(medianAbundance < _circulizer._minUnitigAbundance) return;

			//if(!allContigNodeInGraph(nodepath)) return;

			createJoints(nodepath, medianAbundance);

			colorFile.close();
			getchar();
		}

		
		void createJoints(const vector<u_int32_t>& contigNodepath, u_int32_t sourceAbundance){

			cout << "------------" << endl;
			u_int32_t minSupportingReads = 1; //max((u_int32_t)2, (u_int32_t)(sourceAbundance*0.1));
			cout << "Source abundance: " << sourceAbundance << endl;
			cout << "Min supporint reads: " << minSupportingReads << endl;

			createJoint(contigNodepath, minSupportingReads, true, sourceAbundance);
			createJoint(contigNodepath, minSupportingReads, false, sourceAbundance);
		}

		void createJoint(const vector<u_int32_t>& contigNodepath, u_int32_t minSupportingReads, bool extendRight, float sourceAbundance){
			
			vector<u_int32_t> nodePathJoint;
			//u_int32_t prevUnitigIndex = -1;
			//vector<u_int32_t> unitigPath;
			//vector<u_int32_t> nodepath;

			cout << endl << "\nStart extend" << endl;

			//ofstream colorFile(_inputDir + "/color.csv");
			//colorFile << "Name,Color" << endl;
			//for(u_int32_t nodeName : contigNodepath){
			//	colorFile << BiGraph::nodeIndex_to_nodeName(nodeName) << ",red" << endl;
			//}

			u_int64_t maxLength = 20000;

			if(contigNodepath.size() < 4) return;

			u_int32_t sourceNodeIndex;
			u_int32_t prevNodeIndex;
			bool useSuccessors = false;
			u_int32_t destNodeIndex;

			if(extendRight){
				sourceNodeIndex = contigNodepath[contigNodepath.size()-1];
				prevNodeIndex = contigNodepath[contigNodepath.size()-2];
				destNodeIndex = contigNodepath[0];
			}
			else{
				sourceNodeIndex = contigNodepath[0];
				prevNodeIndex = contigNodepath[1];
				destNodeIndex = contigNodepath[contigNodepath.size()-1];
			}

			//if(!isNotRepeat(sourceNodeIndex, sourceAbundance)) return;
			
			//u_int32_t sourceUnitigIndex = _circulizer._graph->nodeIndex_to_unitigIndex(sourceNodeIndex);
			//cout << sourceUnitigIndex << endl;
			//unitigPath.push_back(sourceUnitigIndex);
			//u_int32_t prevUnitigIndex = sourceUnitigIndex;

			vector<u_int32_t> successors;
			_circulizer._graph->getSuccessors(sourceNodeIndex, 0, successors);

			if(std::find(successors.begin(), successors.end(), prevNodeIndex) == successors.end()){
				useSuccessors = true;
			}
			else{
				useSuccessors = false;
			}


			u_int32_t sourceNodeName = BiGraph::nodeIndex_to_nodeName(sourceNodeIndex);


			u_int64_t length = 0;
			
			u_int32_t currentNodeIndex = sourceNodeIndex;

			unordered_set<u_int32_t> visitedNodeNames;
			//u_int32_t lastSolvedUnitigIndex = -1;
			u_int32_t nbPathSolved = 0;

			while(true){
				
				//if(length > maxLength) break;

				vector<NodeSuccessor> nodeSuccessors;
				getSuccessors(currentNodeIndex, sourceNodeIndex, nodeSuccessors, useSuccessors, visitedNodeNames);


				if(nodeSuccessors.size() == 0){
					cout << "No succ" << endl;
					break;
				}

				std::sort(nodeSuccessors.begin(), nodeSuccessors.end(), NodeSuccessorComparator);

				cout << endl << "\t\t---------------------" << endl;
				cout << "\t\tCurrent node: " << BiGraph::nodeIndex_to_nodeName(currentNodeIndex) << "(length: " << length << ")" << endl;
				for(const NodeSuccessor& nodeSuccessor : nodeSuccessors){
					cout << "\t\t\t" << BiGraph::nodeIndex_to_nodeName(nodeSuccessor._nodeIndex) << " (" << nodeSuccessor._nbSupportingReads << ")" << endl;
				}

				const NodeSuccessor& bestSuccessorNode = nodeSuccessors[0];

				if(bestSuccessorNode._nbSupportingReads < minSupportingReads){
					cout << "No supp" << endl;
					break;
				}

				bool hasSolved = false;

				if(nodeSuccessors.size() > 1){
					float supportRatio = ((float)nodeSuccessors[1]._nbSupportingReads) / bestSuccessorNode._nbSupportingReads;
					//cout << supportRatio << endl;
					if(supportRatio > 0.4){
						cout << "Ambigous" << endl;
						break;
					}
					//if(nodeSuccessors[1]._nbSupportingReads > 0) break;
					hasSolved = true;
					nbPathSolved += 1;
				}

				u_int32_t bestSuccessorNodeName = BiGraph::nodeIndex_to_nodeName(bestSuccessorNode._nodeIndex);
				//if(visitedNodeNames.find(bestSuccessorNodeName) != visitedNodeNames.end()){
					//if(nbPathSolved > 1) break;
				//}

				//visitedNodeNames.insert(bestSuccessorNodeName);

				//break;
				currentNodeIndex = bestSuccessorNode._nodeIndex;
				length += (_circulizer._kminmerLengthMean - _circulizer._kminmerOverlapMean);
				
				//colorFile << BiGraph::nodeIndex_to_nodeName(currentNodeIndex) << ",green" << endl;

				u_int32_t unitigIndex = _circulizer._graph->nodeIndex_to_unitigIndex(bestSuccessorNode._nodeIndex);
				


				//if(unitigIndex != prevUnitigIndex){
					//if(hasSolved) lastSolvedUnitigIndex = unitigIndex;
					//prevUnitigIndex = unitigIndex;
					//unitigPath.push_back(unitigIndex);

					//if(unitigPath.size() == 5) break;
				//}

				if(isNotRepeat(currentNodeIndex, sourceAbundance)){
					sourceNodeIndex = currentNodeIndex;
				}

				nodePathJoint.push_back(currentNodeIndex);

				if(BiGraph::nodeIndex_to_nodeName(currentNodeIndex) == BiGraph::nodeIndex_to_nodeName(destNodeIndex)){

					cout << "Found circular path!" << endl;

					vector<u_int32_t> nodePathComplete = contigNodepath;
					nodePathComplete.insert(nodePathComplete.end(), nodePathJoint.begin(), nodePathJoint.end());
					
					_circulizer.writeContig(nodePathComplete);
					
					break;
				}

				if(length > 100000) break;
			}

			getchar();
			//colorFile.close();

			//getchar();

			/*
			if(unitigPath.size() == 5){
				for(size_t i=0; i<unitigPath.size(); i++){
					cout << unitigPath[i] << " ";
				}
				cout << endl;
				//getchar();
			}
			if(unitigPath.size() == 5 && (unitigPath[1] == unitigPath[3]) && (unitigPath[0] != unitigPath[2]) && (unitigPath[4] != unitigPath[2])) {// && lastSolvedUnitigIndex != -1){


				vector<u_int32_t> nodepathComplete;
				size_t n = 100;

				if(extendRight){
					//nodepathComplete = _graph->_unitigs[unitigPath[0]]._nodes;
				}
				else{
					
					std::reverse(unitigPath.begin(), unitigPath.end());
					
				}

				vector<u_int32_t> nodepathStarting = _circulizer._graph->_unitigs[unitigPath[0]]._nodes;
				vector<u_int32_t> nodepathStartingN(nodepathStarting.end() - std::min(nodepathStarting.size(), n), nodepathStarting.end());

				nodepathComplete = nodepathStartingN;

				for(size_t i=1; i<unitigPath.size()-1; i++){
					u_int32_t unitigIndex = unitigPath[i];
					//cout << unitigIndex << endl;
					const vector<u_int32_t>& nodepath = _circulizer._graph->_unitigs[unitigIndex]._nodes;
					nodepathComplete.insert(nodepathComplete.end(), nodepath.begin(), nodepath.end());

					//if(unitigIndex == lastSolvedUnitigIndex) break;
				}

				
				const Unitig& lastUnitig = _circulizer._graph->_unitigs[unitigPath[unitigPath.size()-1]];
				//if(lastUnitig._index == lastSolvedUnitigIndex){
				vector<u_int32_t> nodepathEnd = lastUnitig._nodes;
				vector<u_int32_t> nodepathEndN(nodepathEnd.begin(), nodepathEnd.begin() + std::min(nodepathEnd.size(), n));
				nodepathComplete.insert(nodepathComplete.end(), nodepathEndN.begin(), nodepathEndN.end());
				//}
				

				u_int64_t size = nodepathComplete.size();
				bool isCircular = false;

				_file_joints.write((const char*)&size, sizeof(size));
				_file_joints.write((const char*)&isCircular, sizeof(isCircular));
				_file_joints.write((const char*)&nodepathComplete[0], size * sizeof(u_int32_t));

				for(size_t i=0; i<nodepathComplete.size()-1; i++){
					u_int32_t nodeIndex = nodepathComplete[i];
					
					vector<u_int32_t> successors;
					_graph->getSuccessors(nodeIndex, 0, successors);

					if(std::find(successors.begin(), successors.end(), nodepathComplete[i+1]) == successors.end()){
						cout << "invalid path" << endl;
						getchar();
					}

				}

				cout << "Found repeat path" << endl;
				getchar();
			}
			*/
		}

		struct NodeSuccessor{
			u_int32_t _nodeIndex;
			u_int32_t _nbSupportingReads;
		};

		static bool NodeSuccessorComparator(const NodeSuccessor &a, const NodeSuccessor &b){
			return a._nbSupportingReads > b._nbSupportingReads;
		}

		void getSuccessors(u_int32_t nodeIndex, u_int32_t sourceNodeIndex, vector<NodeSuccessor>& nodeSuccessors, bool useSuccessors, unordered_set<u_int32_t>& visitedNodeNames){

			const vector<u_int32_t>& sourceReadIndexes = _circulizer._nodeName_to_readIndexes[BiGraph::nodeIndex_to_nodeName(sourceNodeIndex)];

			nodeSuccessors.clear();

			vector<u_int32_t> successors;

			if(useSuccessors){
				_circulizer._graph->getSuccessors(nodeIndex, 0, successors);
			}
			else{
				_circulizer._graph->getPredecessors(nodeIndex, 0, successors);
			}

			for(u_int32_t nodeIndex : successors){

				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
				if(visitedNodeNames.find(nodeName) != visitedNodeNames.end()) continue;

				const vector<u_int32_t>& readIndexes = _circulizer._nodeName_to_readIndexes[nodeName];
				u_int32_t nbSupportingReads = Utils::getNbSharedElements(readIndexes, sourceReadIndexes);

				nodeSuccessors.push_back({nodeIndex, nbSupportingReads});
			}
		}
		

		bool isNotRepeat(u_int32_t nodeIndex, u_int32_t sourceAbundance){
			float abundance = _circulizer._graph->nodeIndex_to_unitig(nodeIndex)._abundance; //_circulizer._graph->getNodeUnitigAbundance(nodeIndex);
			float variance = sourceAbundance*0.2;
			//cout << abundance << " " << sourceAbundance << endl;
			if(abundance > sourceAbundance-variance && abundance < sourceAbundance+variance){
				return true;
			}

			return false;
		}

	};

};	


#endif 
