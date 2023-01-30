

#ifndef MDBG_METAG_CIRCULIZER
#define MDBG_METAG_CIRCULIZER


#include "Commons.hpp"
#include "graph/GfaParser.hpp"
#include "graph/GraphSimplify.hpp"










class Circulizer : public Tool{

public:

	string _inputDir;
	string _passDir;
	string _truthInputFilename;
	string _outputFilename;
	string _outputFilename_complete;
	bool _debug;
	bool _isFinalAssembly;
	string _inputFilename_unitigNt;
	string _inputFilename_unitigCluster;
	string _filename_abundance;
	MDBG* _mdbg;

	float _minimizerDensity;
    size_t _minimizerSize;
    size_t _kminmerSize;
    size_t _kminmerSizeFirst;
	float _minimizerSpacingMean;
	float _kminmerLengthMean;
	float _kminmerOverlapMean;
	vector<UnitigData> _unitigDatas;

	string _partitionDir;
	string _filename_solidKminmers;
	string _filename_inputContigs;
	string _filename_outputContigs;
	string _filename_readMinimizers;
	string _filename_hifiasmGroundtruth;
	string _filename_contigs;
	string _referenceFilename;
	unordered_map<KmerVec, u_int16_t> _evaluation_hifiasmGroundTruth;
	unordered_map<KmerVec, u_int32_t> _evaluation_hifiasmGroundTruth_position;
	unordered_map<u_int32_t, u_int32_t> _evaluation_hifiasmGroundTruth_nodeNamePosition;
	gzFile _outputContigFile;
	gzFile _outputContigFile_complete;
	GraphSimplify* _graph;
	//ToBasespaceOnTheFly _toBasespace;
	//ContigFeature _contigFeature;
	string _gfaFilename;
	size_t _nbCores;

	unordered_set<u_int32_t> _initialRepeatedNodeNames;

	Circulizer(): Tool (){

	}

	~Circulizer(){

	}


	void parseArgs(int argc, char* argv[]){

		args::ArgumentParser parser("circ", ""); //"This is a test program.", "This goes after the options."
		args::Positional<std::string> arg_outputDir(parser, "outputDir", "", args::Options::Required);
		//args::Positional<std::string> arg_contigs(parser, "contigs", "", args::Options::Required);
		//args::PositionalList<std::string> arg_readFilenames(parser, "reads", "Input filename(s) (separated by space)", args::Options::Required);
		//args::ValueFlag<int> arg_l(parser, "", "Minimizer length", {ARG_MINIMIZER_LENGTH2}, 13);
		//args::ValueFlag<float> arg_d(parser, "", "Minimizer density", {ARG_MINIMIZER_DENSITY2}, 0.005f);
		//args::ValueFlag<std::string> arg_contigs(parser, "", "", {ARG_INPUT_FILENAME_CONTIG}, "");
		args::ValueFlag<string> arg_referenceFilename(parser, "", "reference filename", {"ref"}, "");
		args::ValueFlag<int> arg_nbCores(parser, "", "Number of cores", {ARG_NB_CORES2}, NB_CORES_DEFAULT_INT);
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
		_passDir = _inputDir + "/pass_k10/";
		//_passDir = _inputDir + "/pass_k4/";
		//_filename_inputContigs = args::get(arg_contigs);

		_nbCores = args::get(arg_nbCores);

		_filename_contigs = _inputDir + "/contig_data.txt";

		_referenceFilename = "";
		if(arg_referenceFilename){
			_referenceFilename = args::get(arg_referenceFilename);
		}

		string filename_parameters = _passDir + "/parameters.gz";
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

		//_kminmerSize = _kminmerSizeFirst;

		_logFile << endl;
		_logFile << "Input dir: " << _inputDir << endl;
		//cout << "Output filename: " << _outputFilename << endl;
		_logFile << "Minimizer length: " << _minimizerSize << endl;
		_logFile << "Kminmer length: " << _kminmerSize << endl;
		_logFile << "Density: " << _minimizerDensity << endl;
		_logFile << endl;

		_outputFilename = _inputDir + "/minimizer_contigs.gz";
		_outputFilename_complete = _inputDir + "/minimizer_contigs_complete.gz";
		_filename_readMinimizers = _inputDir + "/read_data.txt";
		_filename_hifiasmGroundtruth = _inputDir + "/hifiasmGroundtruth.gz";
		_filename_outputContigs = _inputDir + "/contigs.min.gz";
		_filename_solidKminmers = _inputDir + "/solid.min.gz";

	}

	void execute (){

		_nbPathSolved = 0;
		_nbLinearSolved = 0;

		collectStartingContigs();
		loadGraph();
		//exit(1);
		//generateUnitigs();
		//generateContigs2(_inputDir + "/contigs.nodepath");
		//generateContigs3();

		cout << "Nb path solved: " << _nbPathSolved << endl;
		closeLogFile();
	}

	struct ReferenceContig{
		vector<u_int32_t> _nodePath;
	};

	vector<ReferenceContig> _referenceContigs;

	ofstream _referencePathFile;

	void loadReferences(){

		if(_referenceFilename.empty()) return;

		auto fp = std::bind(&Circulizer::loadReferences_read, this, std::placeholders::_1);
		ReadParser readParser(_referenceFilename, true, false, _logFile);
		readParser.parse(fp);

	}

	void loadReferences_read(const Read& read){

		EncoderRLE _encoderRLE;
		MinimizerParser* minimizerParser = new MinimizerParser(_minimizerSize, _minimizerDensity);

		u_int64_t readIndex = read._index;

		string rleSequence;
		vector<u_int64_t> rlePositions;
		_encoderRLE.execute(read._seq.c_str(), read._seq.size(), rleSequence, rlePositions);

		vector<u_int64_t> minimizers;
		vector<u_int64_t> minimizers_pos;
		minimizerParser->parse(rleSequence, minimizers, minimizers_pos);

		vector<KmerVec> kminmers; 
		vector<ReadKminmer> kminmersInfos;
		MDBG::getKminmers(_minimizerSize, _kminmerSize, minimizers, minimizers_pos, kminmers, kminmersInfos, rlePositions, readIndex, false);

		vector<u_int32_t> nodePath;
		
		for(size_t i=0; i<kminmers.size(); i++){

			const KmerVec& vec = kminmers[i];

			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()){
				continue;
			}


			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
			u_int32_t nodeIndex = BiGraph::nodeName_to_nodeIndex(nodeName, !kminmersInfos[i]._isReversed);

			nodePath.push_back(nodeIndex);
		}

		_referenceContigs.push_back({nodePath});

	}

	void computeReferencePath(){
		
		_referencePathFile = ofstream(_passDir + "/refPath.csv");
		_referencePathFile << "Name,Pos" << endl;

		for(const ReferenceContig& refContig : _referenceContigs){

			unordered_map<u_int32_t, vector<u_int32_t>> unitigPos;
			u_int32_t pos = 0;
			u_int32_t currentUnitig = -1;

			for(u_int32_t nodeIndex : refContig._nodePath){
				u_int32_t unitigIndex = _nodeIndex_to_unitigIndex[nodeIndex];
				if(unitigIndex != currentUnitig){
					unitigPos[unitigIndex].push_back(pos);
					currentUnitig = unitigIndex;
					pos += 1;
				}
			}

			for(auto& it : unitigPos){

				string posStr = "";
				for(u_int32_t pos : it.second){
					posStr += to_string(pos) + "-";
				}
				posStr.pop_back();

				_referencePathFile << it.first << "," << posStr << endl;
				_referencePathFile << _progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(it.first) << "," << posStr << endl;
			}

		}
		
		_referencePathFile.close();
	}

	void loadGraph(){



		cout << "load graph" << endl;
		_gfaFilename = _passDir + "/minimizer_graph.gfa";
		//string gfa_filename_unitigs = _inputDir + "/minimizer_graph_unitigs.gfa";

		GraphSimplify* graphSimplify = new GraphSimplify(_gfaFilename, _inputDir, 0, _kminmerSize, _nbCores, _kminmerLengthMean, _kminmerOverlapMean, _logFile);
		_graph = graphSimplify;
		
		cout << _graph->_graphSuccessors->_nbNodes << endl;
		ifstream kminmerFile(_passDir + "/kminmerData_min.txt");
        _graph->_graphSuccessors->_nodeDatas.resize(_graph->_graphSuccessors->_nbNodes/2, {0, 0});

		u_int64_t abSum = 0;
		u_int64_t qualSum = 0;

		while (true) {

			vector<u_int64_t> minimizerSeq;
			minimizerSeq.resize(_kminmerSize);
			kminmerFile.read((char*)&minimizerSeq[0], minimizerSeq.size()*sizeof(u_int64_t));

			if(kminmerFile.eof()) break;

			u_int32_t nodeName;
			u_int32_t abundance;
			//u_int32_t quality;
			//bool isReversed = false;

			kminmerFile.read((char*)&nodeName, sizeof(nodeName));
			kminmerFile.read((char*)&abundance, sizeof(abundance));
			//kminmerFile.read((char*)&quality, sizeof(quality));


			abSum += abundance;

/*

		13684:
		13685:
1124985 875363
875362 1124984

39842572317101447 85402923321509858 15321455750989951 24680273669351159 29439733697870850 73953375187651903 52329015992061288 57251523862126129 87646769359969727 85469972183324879 
32794338962631351 33590639795828423 73953375187651903 29439733697870850 24680273669351159 15321455750989951 85402923321509858 39842572317101447 62213214356044024 60023953785580213 
39842572317101447 85402923321509858 15321455750989951 24680273669351159 29439733697870850 73953375187651903 33590639795828423 90943467114515961 32794338962631351 69948679337545569 
17703429498529502 42223744577734489 81822362348807848 36558248597177661 11030763198038594 52648394164445011 33384815089368208 48616834895537866 88483241758022389 66372980800638455

1124984
1124985
1770939
1770938
*/

			if(BiGraph::nodeName_to_nodeIndex(nodeName, true) == 1770939){
				for(u_int64_t m : minimizerSeq){
					cout << m << " ";
				}
				cout << endl;
			}
			if(BiGraph::nodeName_to_nodeIndex(nodeName, false) == 1770939){
				for(u_int64_t m : minimizerSeq){
					cout << m << " ";
				}
				cout << endl;
			}

			//
			// 
			//qualSum += quality;

			//cout << nodeName << " " << abSum << endl;

			//cout << nodeName << " " << abundance << " " << quality << endl;
			//if(quality == 0){
			//	cout << quality << endl;
			//	getchar();
			//}
			//if(abundance == 1) continue;
			//KmerVec vec;
			//vec._kmers = minimizerSeq;

			_graph->_graphSuccessors->_nodeDatas[nodeName] = {abundance, (u_int32_t)_kminmerLengthMean};
			//_graph->_graphSuccessors->_nodeLengths[nodeName] = _kminmerLengthMean;
			//_kminmerAbundances[vec] = abundance;
		}

		//cout << abSum << " " << qualSum << endl;
		//getchar();

		kminmerFile.close();



		//if(_kminmerSize > 4){
		//GfaParser::binaryGraph_to_gfa(_gfaFilename, _kminmerLengthMean, _kminmerOverlapMean, _gfaFilename+".gfa", _graph->_graphSuccessors->_nodeDatas);
		//if(_kminmerSize > 45){
		//	getchar();
		//}
		//	cout << "lala" << endl;
		//	getchar();
		//}


		_graph->debug_writeGfaErrorfree(2000, 2000, -1, _kminmerSize, false, true, false, _unitigDatas, true, false, false, false, false, true, _mdbg, _minimizerSize, _nbCores, false, false);
		
		//collectPalindrome();
		
		//indexReads();
		/*
		_graph->_unitigGraph->save(_inputDir + "/minimizer_graph_u.gfa", 0);

		ofstream palindromeColor(_inputDir + "/palindromeColor.csv");
		palindromeColor << "Name,Color" << endl;
		for(UnitigGraph::Node* node : _graph->_unitigGraph->_nodes){
			if(node->_unitigIndex == -1) continue;
			for(u_int32_t nodeIndex : node->_nodes){
				if(_isNodeNamePalindrome.find(BiGraph::nodeIndex_to_nodeName(nodeIndex)) != _isNodeNamePalindrome.end()){
					palindromeColor << node->_unitigIndex << ",red" << endl; 
					palindromeColor << _graph->_unitigGraph->unitigIndex_toReverseDirection(node)->_unitigIndex << ",red" << endl; 
				}
			}
		}
		palindromeColor.close();

		cout << "done" << endl;
		getchar();
		*/

		for(UnitigGraph::Node* node : _graph->_unitigGraph->_nodes){
			if(node->_unitigIndex == -1) continue;
			for(u_int32_t nodeIndex : node->_nodes){
				if(_isNodeNamePalindrome.find(BiGraph::nodeIndex_to_nodeName(nodeIndex)) != _isNodeNamePalindrome.end()){
					node->_isPalindrome = true;
					cout << node->_unitigIndex << endl;
				}
			}
		}
		


		GraphCleanedFunctor functor(*this);

        _progressiveAbundanceFilter = new ProgressiveAbundanceFilter(_graph->_unitigGraph, _inputDir, _kminmerSize, true);
        _progressiveAbundanceFilter->execute(functor);

		//cout << _graph->_unitigGraph->computeChecksum_successor() << endl;
		//if(_kminmerSize == 20) exit(1);


        //_unitigGraph->save("/home/gats/workspace/run/unitigGraph_cleaned.gfa");

		
		//_graph->debug_selectUnitigIndex();
		//_graph->debug_writeGfaErrorfree(2000, 2000, -1, _kminmerSize, false, true, false, _unitigDatas, true, false, false, false, false, false, _mdbg, _minimizerSize, _nbCores, false, true);
			
		//delete _mdbg;
		
	}

	struct ReadNodePos{
		u_int32_t _readIndex;
		u_int16_t _pos;
	};

	vector<vector<ReadNodePos>> _nodeName_to_readIndexes;
	unordered_set<u_int32_t> _startEndNodeNames;

	void indexReads(){

		cout << "Indexing reads" << endl;

		//UnitigGraph::Node* n1 = _progressiveAbundanceFilter->_unitigGraph->_nodes[13684];
		//UnitigGraph::Node* n2 = _progressiveAbundanceFilter->_unitigGraph->_nodes[4580];
		//cout << n1->startNode() << " " << n1->endNode() << endl;
		//cout << n2->startNode() << " " << n2->endNode() << endl;

		for(UnitigGraph::Node* node : _graph->_unitigGraph->_nodes){
			if(node->_unitigIndex == -1) continue;
			_startEndNodeNames.insert(BiGraph::nodeIndex_to_nodeName(node->startNode()));
			_startEndNodeNames.insert(BiGraph::nodeIndex_to_nodeName(node->endNode()));
		}

		_nodeName_to_readIndexes.resize(_graph->_graphSuccessors->_nbNodes);

		string mdbg_filename = _passDir + "/kminmerData_min.txt";
		_mdbg = new MDBG(_kminmerSize);
		_mdbg->load(mdbg_filename, false);

		KminmerParserParallel parser(_filename_readMinimizers, _minimizerSize, _kminmerSize, false, false, _nbCores);
		parser.parse(IndexReadsFunctor(*this, true));

		cout << "Selected reads: " << _validReadIndex.size() << endl;
		//_nodeName_to_readIndexes.clear();
		//_nodeName_to_readIndexes.resize(_graph->_graphSuccessors->_nbNodes/2);

		//KminmerParserParallel parser2(_filename_readMinimizers, _minimizerSize, _kminmerSize, false, false, _nbCores);
		//parser2.parse(IndexReadsFunctor(*this, false));

		_validReadIndex.clear();

		loadReferences();

		delete _mdbg;
		_startEndNodeNames.clear();

		/*
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
		*/
	}

	unordered_set<u_int32_t> _validReadIndex;
	unordered_set<u_int32_t> _isNodeNamePalindrome;

	class IndexReadsFunctor {

		public:

		Circulizer& _parent;
		bool _collectValidReads;

		IndexReadsFunctor(Circulizer& parent, bool collectValidReads) : _parent(parent){
			_collectValidReads = collectValidReads;
		}

		IndexReadsFunctor(const IndexReadsFunctor& copy) : _parent(copy._parent){
			_collectValidReads = copy._collectValidReads;
		}

		~IndexReadsFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {

			//if(kminmerList._kminmersInfo.size() < 40) return;
			unordered_set<u_int32_t> nodeNameExists;

			u_int32_t readIndex = kminmerList._readIndex;
			//if(readIndex == 21574){
			//	cout << "----- " << readIndex << endl;
			//}

			//if(!_collectValidReads){
			//	if(_parent._validReadIndex.find(readIndex) == _parent._validReadIndex.end()) return;
			//}

/*
		Read index: 8755
		Read index: 733364
		Read index: 1079023
		Read index: 2168951
		Read index: 3798840
		Read index: 71737
		Read index: 815250
		Read index: 1660626
		Read index: 2327457
		Read index: 2362219
		Read index: 3839166
		Read index: 4370311
		Read index: 4442808
		Read index: 4505767
		Read index: 4556509
		Read index: 5145133
		Read index: 5247811
		Read index: 5435902
		407416 18
		Read index: 8755
		Read index: 5247811
		Read index: 5435902
		*/
		/*
		unordered_set<u_int64_t> nodes;
		  
		nodes.insert(32794338962631351);
		nodes.insert(33590639795828423);
		nodes.insert(73953375187651903);

		nodes.insert(17703429498529502);
		nodes.insert(42223744577734489);
		nodes.insert(81822362348807848);

		unordered_set<u_int64_t> lala;
		for(size_t i=0; i<kminmerList._readMinimizers.size(); i++){
			u_int64_t m = kminmerList._readMinimizers[i];
			
			if(nodes.find(m) != nodes.end()){
				lala.insert(m);
			}
		}

		if(lala.size() == nodes.size()){
			cout << "RRR: " << readIndex << " " << kminmerList._kminmersInfo.size() << endl;
		}
		*/


			if(readIndex == 1677522){
				
				#pragma omp critical
				{		
					/*	
					for(u_int16_t i=0; i<kminmerList._kminmersInfo.size(); i++){
				
					
						//_logFile << readIndex << " " << i << endl;
						const ReadKminmerComplete& kminmerInfo = kminmerList._kminmersInfo[i];

						KmerVec vec = kminmerInfo._vec;
						
						if(_parent._mdbg->_dbg_nodes.find(vec) == _parent._mdbg->_dbg_nodes.end()){
							continue;
						}
						
						u_int32_t nodeName = _parent._mdbg->_dbg_nodes[vec]._index;
						u_int32_t nodeIndex = BiGraph::nodeName_to_nodeIndex(nodeName, kminmerInfo._isReversed);
						cout << nodeIndex << " ";
					}
					cout << endl;
					*/

					for(u_int64_t m : kminmerList._readMinimizers){
						cout << m << " ";
					}
					cout << endl;
				}
			}
			
			

			for(u_int16_t i=0; i<kminmerList._kminmersInfo.size(); i++){
			
				
				//_logFile << readIndex << " " << i << endl;
				const ReadKminmerComplete& kminmerInfo = kminmerList._kminmersInfo[i];

				KmerVec vec = kminmerInfo._vec;
				
				if(_parent._mdbg->_dbg_nodes.find(vec) == _parent._mdbg->_dbg_nodes.end()){
					continue;
				}
				u_int32_t nodeName = _parent._mdbg->_dbg_nodes[vec]._index;
				

				//if(_parent._startEndNodeNames.find(nodeName) == _parent._startEndNodeNames.end()) continue;

				
				u_int32_t nodeIndex = BiGraph::nodeName_to_nodeIndex(nodeName, !kminmerInfo._isReversed);
				
				/*
				if(nodeNameExists.find(nodeName) == nodeNameExists.end()){
					nodeNameExists.insert(nodeName);
				}
				else{
					#pragma omp critical(indexReadsFunctor_repeated)
					{
						_parent._initialRepeatedNodeNames.insert(nodeName);

						if(nodeName == 7144){
							cout << "-" << nodeName << endl;
							for(u_int16_t i=0; i<kminmerList._kminmersInfo.size(); i++){
						
							
								//_logFile << readIndex << " " << i << endl;
								const ReadKminmerComplete& kminmerInfo = kminmerList._kminmersInfo[i];

								KmerVec vec = kminmerInfo._vec;
								
								if(_parent._mdbg->_dbg_nodes.find(vec) == _parent._mdbg->_dbg_nodes.end()){
									continue;
								}
								
								u_int32_t nodeName = _parent._mdbg->_dbg_nodes[vec]._index;
								cout << nodeName << " ";
							}
							cout << endl;
							getchar();
						}

					}
				}
				//if(_parent._nodeName_to_readIndexes.find(nodeName) == _parent._nodeName_to_readIndexes.end()) continue;
				*/

				#pragma omp critical(indexReadsFunctor)
				{

					//if(_parent._nodeName_to_readIndexes[nodeIndex].size() < 150){
						_parent._nodeName_to_readIndexes[nodeIndex].push_back({readIndex, i});
					//}




					/*
					if(_collectValidReads){
						if(_parent._nodeName_to_readIndexes[nodeName].size() < 10){
							_parent._nodeName_to_readIndexes[nodeName].push_back(readIndex);
							_parent._validReadIndex.insert(readIndex);
						}

					}
					else{
						_parent._nodeName_to_readIndexes[nodeName].push_back(readIndex);
					}
					*/
				}

				

			}

		}
	};

	/*
	void collectPalindrome(){

		cout << "Collecting palindrome" << endl;

		string mdbg_filename = _inputDir + "/kminmerData_min_init.txt";
		_mdbg = new MDBG(_kminmerSize);
		_mdbg->load(mdbg_filename, false);

		KminmerParserParallel parser(_filename_readMinimizers, _minimizerSize, _kminmerSize, false, false, _nbCores);
		parser.parse(CollectPalindromeFunctor(*this));


		delete _mdbg;

	}

	class CollectPalindromeFunctor {

		public:

		Circulizer& _parent;

		CollectPalindromeFunctor(Circulizer& parent) : _parent(parent){
		}

		CollectPalindromeFunctor(const CollectPalindromeFunctor& copy) : _parent(copy._parent){
		}

		~CollectPalindromeFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {

			unordered_set<u_int32_t> nodeNameExists;

			u_int32_t readIndex = kminmerList._readIndex;

			u_int32_t prevNodeName;
			KmerVec prevVec;

			for(u_int16_t i=0; i<kminmerList._kminmersInfo.size(); i++){
			
				
				//_logFile << readIndex << " " << i << endl;
				const ReadKminmerComplete& kminmerInfo = kminmerList._kminmersInfo[i];

				KmerVec vec = kminmerInfo._vec;
				
				if(_parent._mdbg->_dbg_nodes.find(vec) == _parent._mdbg->_dbg_nodes.end()){
					continue;
				}
				u_int32_t nodeName = _parent._mdbg->_dbg_nodes[vec]._index;
				
				//if((vec.isPalindrome() || (i > 0 && vec.normalize() == prevVec.normalize()))){ 
				if(i > 0 && vec.normalize() == prevVec.normalize()){
					
					#pragma omp critical
					{			
						_parent._isNodeNamePalindrome.insert(nodeName);
						_parent._isNodeNamePalindrome.insert(prevNodeName);
					}
				}

				prevVec = vec;
				prevNodeName = nodeName;

			}

		}
	};
	*/

	struct StartingContig{
		u_int32_t _index;
		vector<u_int32_t> _nodePath;
		float _abundance;
	};

	vector<StartingContig> _startingContigs;

	void collectStartingContigs(){
		
		
		cout << "load starting contigs " << endl;
		
		_mdbg = new MDBG(_kminmerSize);
		_mdbg->load(_passDir + "/kminmerData_min.txt", false);


		KminmerParserParallel parser(_filename_contigs, _minimizerSize, _kminmerSize, false, false, 1);
		parser.parse(CollectStartingContigsFunctor(*this));

		delete _mdbg;

		cout << "done" << endl;
	}


	class CollectStartingContigsFunctor {

		public:

		Circulizer& _graph;

		CollectStartingContigsFunctor(Circulizer& graph) : _graph(graph){

		}

		CollectStartingContigsFunctor(const CollectStartingContigsFunctor& copy) : _graph(copy._graph){

		}

		~CollectStartingContigsFunctor(){
		}


		void operator () (const KminmerList& kminmerList) {


			u_int64_t readIndex = kminmerList._readIndex;
			const vector<u_int64_t>& readMinimizers = kminmerList._readMinimizers;
			const vector<ReadKminmerComplete>& kminmersInfos = kminmerList._kminmersInfo;

			//cout << kminmersInfos.size() << endl;
			if(kminmerList._isCircular) return;
			if(kminmersInfos.size() < 5000) return;

			vector<u_int32_t> nodePath;
			vector<float> abundances;

			for(size_t i=0; i<kminmersInfos.size(); i++){

				const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];
				const KmerVec& vec = kminmerInfo._vec;

				if(_graph._mdbg->_dbg_nodes.find(vec) == _graph._mdbg->_dbg_nodes.end()){
					//_logFile << "Not found kminmer" << endl;
					//getchar();
					continue;
				}



				u_int32_t nodeName = _graph._mdbg->_dbg_nodes[vec]._index;
				u_int32_t nodeIndex = BiGraph::nodeName_to_nodeIndex(nodeName, !kminmerInfo._isReversed);

				abundances.push_back(_graph._mdbg->_dbg_nodes[vec]._abundance);
				nodePath.push_back(nodeIndex);
			}

			float contigAbundance = Utils::compute_median_float(abundances);

			cout << kminmersInfos.size() << " " << contigAbundance << endl;

			_graph._startingContigs.push_back({readIndex, nodePath, contigAbundance});
		}
	};

	unordered_set<u_int32_t> _processedStartingContigs;
	ProgressiveAbundanceFilter* _progressiveAbundanceFilter;
	float _minUnitigAbundance;
	unordered_map<u_int32_t, u_int32_t> _nodeIndex_to_unitigIndex;
	bool _areReadsIndexed;
	unordered_set<u_int32_t> _isStartingUnitigs;
	u_int64_t _nbPathSolved;
	u_int64_t _nbLinearSolved;
	//unordered_set<u_int32_t> _solvedNodeNames;
	
	class GraphCleanedFunctor {

		public:

		u_int32_t _minSupportingReads;
		u_int32_t _minSupportingReads_2;
		float _lastCutoff;

		static bool UnitigComparator(UnitigGraph::Node* p1, UnitigGraph::Node* p2){
			if(p1->_nodes.size() == p2->_nodes.size()){
				if(p1->_abundance == p2->_abundance){
					return p1->startNode() > p2->startNode();
				}
				return p1->_abundance > p2->_abundance;
			}
			return p1->_nodes.size() < p2->_nodes.size();
		}

		Circulizer& _parent;

		GraphCleanedFunctor(Circulizer& parent) : _parent(parent){
			_lastCutoff = -1;
		}

		GraphCleanedFunctor(const GraphCleanedFunctor& copy) : _parent(copy._parent){

		}

		~GraphCleanedFunctor(){
		}


		void operator () (float cutoff) {
			
            float nextCutoff = cutoff;
			nextCutoff = nextCutoff * (1+0.1);

			float minUnitigAbundance = cutoff / 0.5;
			float maxUnitigAbundance = nextCutoff / 0.5;

			//if(cutoff < 100){
			//	cout << "skipping cutoff " << cutoff << endl;
			//	return;
			//}

			
			for(const StartingContig& startingContig : _parent._startingContigs){
				if(startingContig._abundance < 20){
					cout << "skipping starting unitig" << endl;
					continue;
				}

				if(_parent._processedStartingContigs.find(startingContig._index) != _parent._processedStartingContigs.end()) continue;

				/*
				bool isSolved = false;
				for(u_int32_t nodeIndex : startingContig._nodePath){
					if(_parent._solvedNodeNames.find(BiGraph::nodeIndex_to_nodeName(nodeIndex)) != _parent._solvedNodeNames.end()){
						isSolved = true;
						break;
					}
				}
				if(isSolved) continue;
				*/

				if(maxUnitigAbundance > startingContig._abundance){

					//_minSupportingReads = 1;//max((u_int32_t)(startingContig._abundance * 0.01), (u_int32_t)2);
					_minSupportingReads = max((u_int32_t)(startingContig._abundance * 0.01), (u_int32_t)2);

					if(!_parent._areReadsIndexed){
						_parent._areReadsIndexed = true;
						_parent.indexReads();
					}

					ofstream startingContigColor(_parent._passDir + "/startingContigColor.csv");
					startingContigColor << "Name,Color" << endl;

					_parent._nodeIndex_to_unitigIndex.clear();
					_parent._isStartingUnitigs.clear();
					for(UnitigGraph::Node* node : _parent._progressiveAbundanceFilter->_unitigGraph->_nodes){
						if(node->_unitigIndex == -1) continue;
						for(u_int32_t nodeIndex : node->_nodes){
							_parent._nodeIndex_to_unitigIndex[nodeIndex] = node->_unitigIndex;
						}
					}
					
					//"reprise: essayer d'output tous les linear solved pour commencer et check leur quality, ou les comparer au contigs de hifiasmmeta"
					u_int32_t firstUnitig = -1;
					u_int32_t lastUnitig = -1;

					for(size_t i=0; i<startingContig._nodePath.size(); i++){
						u_int32_t nodeIndex = startingContig._nodePath[i];

						_parent._isStartingUnitigs.insert(_parent._nodeIndex_to_unitigIndex[nodeIndex]);
						_parent._isStartingUnitigs.insert(_parent._progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(_parent._nodeIndex_to_unitigIndex[nodeIndex]));

						if(firstUnitig == -1){
							firstUnitig = _parent._nodeIndex_to_unitigIndex[nodeIndex];
						}
						lastUnitig = _parent._nodeIndex_to_unitigIndex[nodeIndex];
						
						startingContigColor << _parent._nodeIndex_to_unitigIndex[nodeIndex] << ",blue" << endl;
						startingContigColor << _parent._nodeIndex_to_unitigIndex[GraphSimplify::nodeIndex_toReverseDirection(nodeIndex)] << ",blue" << endl;
						
					}

					if(firstUnitig == lastUnitig || firstUnitig == _parent._progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(lastUnitig)){
						_parent._nbLinearSolved += 1;
					}
					/*
					for(u_int32_t unitigIndex : _parent._isStartingUnitigs){
						UnitigGraph::Node* unitig = _parent._progressiveAbundanceFilter->_unitigGraph->_nodes[unitigIndex];
						if(isPalindromeTip(unitig)){
							cout << "Palindrome tip: " << unitig->_unitigIndex << endl;
							for(UnitigGraph::Node* nn : _parent._progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(unitig)->_successors){
								cout << "\t" << nn->_unitigIndex << " " << computeSharedElements(unitig, nn) << endl;
							}
						}
					}
					*/

					startingContigColor << firstUnitig << ",green" << endl;
					startingContigColor << _parent._progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(firstUnitig) << ",green" << endl;
					startingContigColor << lastUnitig << ",red" << endl;
					startingContigColor << _parent._progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(lastUnitig) << ",red" << endl;

					startingContigColor.close();

					cout << "Abundance: " << startingContig._abundance << " " << cutoff << endl;

					if(cutoff != _lastCutoff){
						_lastCutoff = cutoff;
						_parent._progressiveAbundanceFilter->_unitigGraph->save(_parent._passDir + "/circGraph.gfa", startingContig._abundance);

						_parent.computeReferencePath();
						determineRepeatedUnitigs();
					}

					solvePath(startingContig);
					cout << "done" << endl;

					_parent._processedStartingContigs.insert(startingContig._index);

					cout << "Nb path solved: " << _parent._nbPathSolved << endl;
					cout << "Nb linear solved: " << _parent._nbLinearSolved << endl;
					//getchar();
				}

			}
			


		}

		vector<bool> _isUnitigUnique;
		
		struct UnitSuccessor{
			UnitigGraph::Node* _node;
			u_int64_t _nbSupportingReads;
		};



		void determineRepeatedUnitigs(){

			_isUnitigUnique.clear();
			_isUnitigUnique.resize(_parent._progressiveAbundanceFilter->_unitigGraph->_nodes.size(), false);

			ofstream ouputFile(_parent._passDir + "/debug.txt");
			vector<UnitigGraph::Node*> unitigs;

			for(UnitigGraph::Node* node : _parent._progressiveAbundanceFilter->_unitigGraph->_nodes){
				if(node->_unitigIndex == -1) continue;
				if(node->_unitigIndex % 2 == 1) continue;

				//if(_parent._nodeName_to_readIndexes[BiGraph::nodeIndex_to_nodeName(node->endNode())].size() > 1000){
				//	_isUnitigUnique[node->_unitigIndex] = false;
				//	continue;
				//}

				unitigs.push_back(node);
				

				
				_isUnitigUnique[node->_unitigIndex] = true;
				_isUnitigUnique[_parent._progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(node)->_unitigIndex] = true;
				/*
				if(node->_successors.size() > 1 || node->_predecessors.size() > 1){

					if(!isRepeatDirect(node)){
            			_isUnitigUnique[node->_unitigIndex] = true;
            			_isUnitigUnique[_parent._progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(node)->_unitigIndex] = true;
					}
					
				}
				else{
					_isUnitigUnique[node->_unitigIndex] = true;
					_isUnitigUnique[_parent._progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(node)->_unitigIndex] = true;
				}
				*/
			}

			/*
			cout << "Repeated: " << endl;
			for(u_int32_t nodeName : _parent._initialRepeatedNodeNames){
				cout << nodeName << endl;
				_isUnitigUnique[_parent._nodeIndex_to_unitigIndex[BiGraph::nodeName_to_nodeIndex(nodeName, true)]] = false;
				_isUnitigUnique[_parent._nodeIndex_to_unitigIndex[BiGraph::nodeName_to_nodeIndex(nodeName, false)]] = false;
			}

			cout << "----" << endl;
			for(u_int32_t nodeIndex : _parent._progressiveAbundanceFilter->_unitigGraph->_nodes[8124]->_nodes){
				//cout << BiGraph::nodeIndex_to_nodeName(nodeIndex) << endl;
				for(ReadNodePos r : _parent._nodeName_to_readIndexes[nodeIndex]){
					cout << "\t" << r._readIndex << endl;
				}
			}
			*/
			
			/*
			UnitigGraph::Node* n1 = _parent._progressiveAbundanceFilter->_unitigGraph->_nodes[90328];
			for(UnitigGraph::Node* nn : n1->_successors){
				cout << n1->_unitigIndex << " -> " << nn->_unitigIndex << " " << computeSharedElements(n1, nn) << endl;
			}
			for(UnitigGraph::Node* nn : n1->_predecessors){
				cout << n1->_unitigIndex << " -> " << nn->_unitigIndex << " " << computeSharedElements(nn, n1) << endl;
			}
			*/

			/*
			UnitigGraph::Node* n1 = _parent._progressiveAbundanceFilter->_unitigGraph->_nodes[6882];
			UnitigGraph::Node* n2 = _parent._progressiveAbundanceFilter->_unitigGraph->_nodes[303540];

			unordered_set<u_int32_t> allReads;
			for(size_t i=0; i<n1->_nodes.size(); i++){
				u_int32_t nodeIndex = n1->_nodes[i];
				
				for(const ReadNodePos& p: _parent._nodeName_to_readIndexes[nodeIndex]){
					allReads.insert(p._readIndex);
				}
				for(const ReadNodePos& p: _parent._nodeName_to_readIndexes[GraphSimplify::nodeIndex_toReverseDirection(nodeIndex)]){
					allReads.insert(p._readIndex);
				}

				
			}

			cout << "-----" << endl;

			//for(size_t i=0; i<n2->_nodes.size(); i++){
			//	cout << i  << ": " << unitig2->_nodes[i] << endl;
			//}
			
			for(size_t i=0; i<n2->_nodes.size(); i++){
				u_int32_t nodeIndex = n2->_nodes[i];

				for(const ReadNodePos& p: _parent._nodeName_to_readIndexes[nodeIndex]){
					if(allReads.find(p._readIndex) != allReads.end()){

						cout << p._readIndex << endl;
						cout << "omg: " << i << endl;
						break;
					}
				}

			}
			*/

			/*
			UnitigGraph::Node* n2 = _parent._progressiveAbundanceFilter->_unitigGraph->_nodes[90329];
			for(UnitigGraph::Node* nn : n2->_successors){
				cout << n2->_unitigIndex << " -> " << nn->_unitigIndex << " " << computeSharedElements(n2, nn) << endl;
			}
			for(UnitigGraph::Node* nn : n2->_predecessors){
				cout << n2->_unitigIndex << " -> " << nn->_unitigIndex << " " << computeSharedElements(n2, nn) << endl;
			}
			*/

			
			
			

			UnitigGraph::Node* n1 = _parent._progressiveAbundanceFilter->_unitigGraph->_nodes[493806];
			UnitigGraph::Node* n2 = _parent._progressiveAbundanceFilter->_unitigGraph->_nodes[74126];
			UnitigGraph::Node* n1_rc = _parent._progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(n1);
			UnitigGraph::Node* n2_rc = _parent._progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(n2);
			//UnitigGraph::Node* n3_rc = _parent._progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(n3);


			//cout << collectSharedElementsOld(n1, n2) << endl;
			//cout << collectSharedElementsOld(n1, n2_rc) << endl;
			//cout << collectSharedElementsOld(n1_rc, n2) << endl;
			//cout << collectSharedElementsOld(n1_rc, n2_rc) << endl;
			/*
			cout << computeSharedElements(n1, n2) << endl;
			cout << computeSharedElements(n1, n2_rc) << endl;
			cout << computeSharedElements(n1_rc, n2) << endl;
			cout << computeSharedElements(n1_rc, n2_rc) << endl;
			getchar();
			*/
			/*
			UnitigGraph::Node* n1 = _parent._progressiveAbundanceFilter->_unitigGraph->_nodes[13685];
			cout << "sdfsdfsdfsdf" << endl;
			cout << n1->endNode() << endl;
			cout << _parent._progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(n1)->startNode() << endl;
			
			
			UnitigGraph::Node* n2 = _parent._progressiveAbundanceFilter->_unitigGraph->_nodes[4580];
			cout << n2->startNode() << endl;
			cout << _parent._progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(n2)->endNode() << endl;
			*/

			//getchar();
			
			//cout << n1->startNode() << " " << n1->endNode() << endl;
			//cout << n2->startNode() << " " << n2->endNode() << endl;
			

			/*
			for(u_int32_t nodeIndex : n1->_nodes){
				cout << nodeIndex << " ";
			}
			cout << endl;
			cout << endl;
			for(u_int32_t nodeIndex : n2_rc->_nodes){
				cout << nodeIndex << " ";
			}
			cout << endl;
			cout << endl;
			for(u_int32_t nodeIndex : n3_rc->_nodes){
				cout << nodeIndex << " ";
			}
			cout << endl;
			cout << endl;


			cout << collectSharedElementsOld(n1, n2) << endl;
			cout << collectSharedElementsOld(n1, n2_rc) << endl;
			cout << collectSharedElementsOld(n1_rc, n2) << endl;
			cout << collectSharedElementsOld(n1_rc, n2_rc) << endl;

			cout << computeSharedElements(n1, n2) << endl;
			cout << computeSharedElements(n1, n2_rc) << endl;
			cout << computeSharedElements(n1_rc, n2) << endl;
			cout << computeSharedElements(n1_rc, n2_rc) << endl;
			
			//getchar();

			
			cout << "----" << endl;
			for(ReadNodePos r : _parent._nodeName_to_readIndexes[n1->startNode()]){
				cout << "\t" << r._readIndex << endl;
			}
			for(ReadNodePos r : _parent._nodeName_to_readIndexes[GraphSimplify::nodeIndex_toReverseDirection(n1->startNode())]){
				cout << "\t" << r._readIndex << endl;
			}
			cout << "----" << endl;
			for(ReadNodePos r : _parent._nodeName_to_readIndexes[n1->endNode()]){
				cout << "\t" << r._readIndex << endl;
			}
			for(ReadNodePos r : _parent._nodeName_to_readIndexes[GraphSimplify::nodeIndex_toReverseDirection(n1->endNode())]){
				cout << "\t" << r._readIndex << endl;
			}
			cout << "----" << endl;
			cout << "----" << endl;
			cout << "----" << endl;
			for(ReadNodePos r : _parent._nodeName_to_readIndexes[n2->startNode()]){
				cout << "\t" << r._readIndex << endl;
			}
			for(ReadNodePos r : _parent._nodeName_to_readIndexes[GraphSimplify::nodeIndex_toReverseDirection(n2->startNode())]){
				cout << "\t" << r._readIndex << endl;
			}
			cout << "----" << endl;
			cout << "----" << endl;
			cout << "----" << endl;
			for(ReadNodePos r : _parent._nodeName_to_readIndexes[n2->endNode()]){
				cout << "\t" << r._readIndex << endl;
			}
			for(ReadNodePos r : _parent._nodeName_to_readIndexes[GraphSimplify::nodeIndex_toReverseDirection(n2->endNode())]){
				cout << "\t" << r._readIndex << endl;
			}
			*/

			cout << "direct done" << endl;
			
			std::sort(unitigs.begin(), unitigs.end(), UnitigComparator);

			for(size_t i=0; i<5; i++){

				cout << "---------------------- " << i << endl;

				for(UnitigGraph::Node* unitig : unitigs){

					if(!_isUnitigUnique[unitig->_unitigIndex]) continue;
					if(isPalindromeTip(unitig)) continue;
					if(unitig->_length > 20000) continue;
					//if(_parent._nodeName_to_readIndexes[BiGraph::nodeIndex_to_nodeName(unitig->endNode())].size() > 1000){
					//	_isUnitigUnique[unitig->_unitigIndex] = false;
					//	continue;
					//}

					UnitigGraph::Node* unitig_rc = _parent._progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(unitig);

					vector<UnitigGraph::Node*> uniqueSuccessors;
					getNextUniqueSuccessors(unitig, uniqueSuccessors);

					ouputFile << unitig->_unitigIndex << endl;
					//if(uniqueSuccessors.size() > 1){
            		//	_isUnitigUnique[unitig->_unitigIndex] = false;
            		//	_isUnitigUnique[_parent._progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(unitig)->_unitigIndex] = false;
					//}

					vector<UnitigGraph::Node*> uniquePredecessors;
					getNextUniqueSuccessors(unitig_rc, uniquePredecessors);


					vector<UnitigGraph::Node*> uniqueSuccessorsDerep;
					dereplicateUniqueSuccessors(unitig, uniqueSuccessors, uniqueSuccessorsDerep);
					
					vector<UnitigGraph::Node*> uniquePredecessorsDerep;
					dereplicateUniqueSuccessors(unitig_rc, uniquePredecessors, uniquePredecessorsDerep);

					//if(uniquePredecessors.size() > 1){
            		//	_isUnitigUnique[unitig->_unitigIndex] = false;
            		//	_isUnitigUnique[_parent._progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(unitig)->_unitigIndex] = false;
					//}
					
					
					for(UnitigGraph::Node* node : uniqueSuccessorsDerep){
						
						//unordered_set<u_int32_t> sharedReads;
						//u_int64_t nbSharedReads = collectSharedElements(_parent._nodeName_to_readIndexes[BiGraph::nodeIndex_to_nodeName(unitig->endNode())], _parent._nodeName_to_readIndexes[BiGraph::nodeIndex_to_nodeName(node->startNode())], sharedReads);

						ouputFile << "\tS: " << node->_unitigIndex << endl; //<< " " << nbSharedReads << endl;
						//for(u_int32_t readIndex : sharedReads){
						//	ouputFile << "\t\t" << readIndex << endl;
						//}
					}
					for(UnitigGraph::Node* node : uniquePredecessorsDerep){

						
						//unordered_set<u_int32_t> sharedReads;
						//u_int64_t nbSharedReads = collectSharedElements(_parent._nodeName_to_readIndexes[BiGraph::nodeIndex_to_nodeName(unitig->endNode())], _parent._nodeName_to_readIndexes[BiGraph::nodeIndex_to_nodeName(node->startNode())], sharedReads);


						ouputFile << "\tP: " << node->_unitigIndex << endl;
						//for(u_int32_t readIndex : sharedReads){
						//	ouputFile << "\t\t" << readIndex << endl;
						//}
					}
					
					if(uniqueSuccessorsDerep.size() == 1 && uniquePredecessorsDerep.size() == 1){
					}
					else{
            			_isUnitigUnique[unitig->_unitigIndex] = false;
            			_isUnitigUnique[_parent._progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(unitig)->_unitigIndex] = false;

						if(_parent._isStartingUnitigs.find(unitig->_unitigIndex) != _parent._isStartingUnitigs.end()){
							cout << "Repeated: " << unitig->_unitigIndex << " " << unitig->_length << endl;
							//getchar();
						}
					}


				}
			}
			
			

			ofstream colorFile(_parent._passDir + "/repeatColor.csv");
			colorFile << "Name,Color" << endl; 
			for(UnitigGraph::Node* node : _parent._progressiveAbundanceFilter->_unitigGraph->_nodes){
				if(node->_unitigIndex == -1) continue;
				if(_isUnitigUnique[node->_unitigIndex]){
					colorFile << node->_unitigIndex << ",blue" << endl;
				}
				else{
					colorFile << node->_unitigIndex << ",red" << endl;
				}
			}

			colorFile.close();
			ouputFile.close();
			
		}

		bool isRepeatDirect(UnitigGraph::Node* node){
			
			UnitigGraph::Node* node_rc = _parent._progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(node);

			vector<UnitigGraph::Node*> uniqueSuccessors;
			for(UnitigGraph::Node* node_nn : node->_successors){
				if(node_nn == node || node_nn == node_rc) continue;
				uniqueSuccessors.push_back(node_nn);
			}


			vector<UnitigGraph::Node*> uniquePredecessors;
			for(UnitigGraph::Node* node_nn :  node_rc->_successors){
				if(node_nn == node || node_nn == node_rc) continue;
				uniquePredecessors.push_back(node_nn);
			}

			vector<UnitigGraph::Node*> uniqueSuccessorsDerep;
			dereplicateUniqueSuccessors(node, uniqueSuccessors, uniqueSuccessorsDerep);
			
			vector<UnitigGraph::Node*> uniquePredecessorsDerep;
			dereplicateUniqueSuccessors(node_rc, uniquePredecessors, uniquePredecessorsDerep);

			if(uniqueSuccessorsDerep.size() == 1 && uniquePredecessorsDerep.size() == 1) return false;

			return true;
		}

		
		void getNextUniqueSuccessors(UnitigGraph::Node* unitigSource, vector<UnitigGraph::Node*>& uniqueSuccessors){

			//u_int32_t source_nodeIndex = unitigSource->endNode();
			//vector<u_int32_t> source_readIndexes = _parent._nodeName_to_readIndexes[BiGraph::nodeIndex_to_nodeName(source_nodeIndex)]; //unitigDatas[BiGraph::nodeIndex_to_nodeName(source_nodeIndex)]._readIndexes;

			//if(unitigSource->_unitigIndex == 3484 || unitigSource->_unitigIndex == 3485){
			//	cout << "---- collect" << endl;
			//}

			uniqueSuccessors.clear();

			unordered_set<UnitigGraph::Node*> isVisited;
			queue<UnitigGraph::Node*> queue;

			queue.push(unitigSource);
			
			while (!queue.empty()){

				UnitigGraph::Node* unitig = queue.front();
				queue.pop();

				//if(unitigSource->_unitigIndex == 420008 || unitigSource->_unitigIndex == 420009){
				//	cout << "\t" << unitig->_unitigIndex << endl;
				//}

				if(isVisited.find(unitig) != isVisited.end()) continue;
				isVisited.insert(unitig);
				//isVisited.insert(_parent._progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(unitig));


				if(unitig != unitigSource){// && unitig != _parent._progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(unitigSource)){

					//u_int64_t minSupportingReads = _minSupportingReads;
					//if(_isUnitigUnique[unitig->_unitigIndex]){
					//minSupportingReads = 1;
					//}
					//const vector<u_int32_t>& readIndexes = _parent._nodeName_to_readIndexes[BiGraph::nodeIndex_to_nodeName(unitig->startNode())];
					//u_int64_t nbSharedReads = Utils::computeSharedElements(source_readIndexes, readIndexes);
				
					//unordered_set<u_int32_t> sharedReads;
					//u_int64_t nbSharedReads = Utils::collectSharedElements(source_readIndexes, readIndexes, sharedReads);

					//if(source_unitigIndex == unitigIndexLala || unitigIndex_toReverseDirection(source_unitigIndex) == unitigIndexLala){
					//    cout << BiGraph::nodeIndex_to_nodeName(source_nodeIndex) << " " << BiGraph::nodeIndex_to_nodeName(_unitigs[unitigIndex]._startNode) << " " <<  nbSharedReads << endl;
					//    //getchar();
					//}
					//281313

					u_int64_t nbSharedReads = computeSharedElements(unitigSource, unitig);

					if(nbSharedReads < _minSupportingReads) continue;

					if(unitigSource->_unitigIndex == 13685){
						cout << "\t\t" << unitig->_unitigIndex << " "<< nbSharedReads << endl;
					}
					if(_isUnitigUnique[unitig->_unitigIndex]){
						
						//for(u_int32_t readIndex : sharedReads){
        				//	source_readIndexes.erase(std::remove(source_readIndexes.begin(), source_readIndexes.end(), readIndex), source_readIndexes.end());
							
						//}
						uniqueSuccessors.push_back(unitig);
						continue;
					}

				}

				//if(unitigSource->_unitigIndex == 420008 || unitigSource->_unitigIndex == 420009){
				//	cout << "\t" << unitig->_successors.size() << endl;
				//}

				for(UnitigGraph::Node* nodeS : unitig->_successors){
					queue.push(nodeS);
				}

				


			}

			//file_scc.close();

		}
		
		bool isSuccessor(UnitigGraph::Node* unitig1, UnitigGraph::Node* unitig2, UnitigGraph::Node* unitigCurrent){
			
			//if((unitig1->_unitigIndex == 416192 && unitig2->_unitigIndex == 416193)){
			//	cout << "---- isSuccessor" << endl;
			//}


			//u_int32_t source_nodeIndex = unitigSource->endNode();

			unordered_set<UnitigGraph::Node*> isVisited;
			isVisited.insert(unitigCurrent);

			queue<UnitigGraph::Node*> queue;

			queue.push(unitig1);
			
			while (!queue.empty()){

				UnitigGraph::Node* unitig = queue.front();
				queue.pop();



				if (isVisited.find(unitig) != isVisited.end()) continue;
				isVisited.insert(unitig);
				//isVisited.insert(_parent._progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(unitig));


				if(unitig != unitig1){

					u_int64_t nbSharedReads = computeSharedElements(unitig1, unitig);

					//if((unitig1->_unitigIndex == 416192 && unitig2->_unitigIndex == 416193)){
					//	cout << "\t" << unitig->_unitigIndex << " " << nbSharedReads << endl;
					//}

					//if((unitig1->_unitigIndex == 18967 || unitig1->_unitigIndex == 18966) && (unitig2->_unitigIndex == 7987 || unitig2->_unitigIndex == 7986)){
					//	cout << " \t -" << unitig->_unitigIndex << " " << nbSharedReads << endl;
					//}

					if(nbSharedReads < _minSupportingReads) continue;

					if(unitig == unitig2) return true;
					//if(_isUnitigUnique[unitig->_unitigIndex]){
						
						//for(u_int32_t readIndex : sharedReads){
        				//	source_readIndexes.erase(std::remove(source_readIndexes.begin(), source_readIndexes.end(), readIndex), source_readIndexes.end());
							
						//}
					//	uniqueSuccessors.push_back(unitig);
					//	continue;
					//}

				}

				for(UnitigGraph::Node* nodeS : unitig->_successors){
					queue.push(nodeS);
				}
			}

			return false;
		}

		static bool UniqueSuccessorComparator(const UnitSuccessor& p1, const UnitSuccessor& p2){
			if(p1._nbSupportingReads == p2._nbSupportingReads){
				return p1._node->_nodes.size() < p2._node->_nodes.size();
			}
			return p1._nbSupportingReads > p2._nbSupportingReads;
		}




		void dereplicateUniqueSuccessors(UnitigGraph::Node* unitigSource, const vector<UnitigGraph::Node*>& uniqueSuccessors, vector<UnitigGraph::Node*>& uniqueSuccessorsDerep){

			//u_int64_t minSupportingReads = 1;

			
			uniqueSuccessorsDerep.clear();
			vector<UnitSuccessor> uniqueSuccessorsSup;

			//vector<ReadNodePos> source_readIndexes;
			//vector<ReadNodePos> source_readIndexes_rev;
			//source_readIndexes = _parent._nodeName_to_readIndexes[unitigSource->endNode()];
			//source_readIndexes_rev = _parent._nodeName_to_readIndexes[_parent._progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(unitigSource)->startNode()];
			//"reprise: isSuccessor ne marche pas pour le repeat loop side"
			
			if(unitigSource->_unitigIndex == 23062 || unitigSource->_unitigIndex == 23063){
				cout << "-----" << endl;
				cout << unitigSource->_unitigIndex << endl;
				//cout << unitigSource->endNode() << " " << _parent._progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(unitigSource)->startNode() << endl;
				//for(const ReadNodePos& r : source_readIndexes){
				//	cout << "\t" << r._readIndex << endl;
				//}
				//for(const ReadNodePos& r : source_readIndexes_rev){
				//	cout << "\t" << r._readIndex << endl;
				//}
			}
			
			
			for(UnitigGraph::Node* uniqueSuccessor : uniqueSuccessors){

				/*
				const vector<ReadNodePos>& readIndexes = _parent._nodeName_to_readIndexes[uniqueSuccessor->startNode()];
				const vector<ReadNodePos>& readIndexes_rev = _parent._nodeName_to_readIndexes[_parent._progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(uniqueSuccessor)->endNode()];
			
			
				vector<u_int32_t> sharedReads;
				u_int64_t nbSharedReads = collectSharedElements(source_readIndexes, readIndexes, true, sharedReads);



				//u_int64_t nbSharedReads = computeSharedElements(source_readIndexes, readIndexes, true);
				//uniqueSuccessorsSup.push_back({uniqueSuccessor, nbSharedReads});

				//nbSharedReads += computeSharedElements(source_readIndexes_rev, readIndexes_rev, false);
				//uniqueSuccessorsSup.push_back({uniqueSuccessor, nbSharedReads});

				vector<u_int32_t> sharedReads_rev;
				nbSharedReads += collectSharedElements(source_readIndexes_rev, readIndexes_rev, false, sharedReads_rev);

				//if(unitigSource->_unitigIndex == 456998 || unitigSource->_unitigIndex == 456999){
				//	cout << "loulou " << nbSharedReads << endl;
				//}

	
				if(unitigSource->_unitigIndex == 4800 || unitigSource->_unitigIndex == 4801){
					cout << "\t" << uniqueSuccessor->_unitigIndex << " " << nbSharedReads << endl;
					for(u_int32_t readIndex : sharedReads){
						cout << "\t\t" << readIndex << endl;
					}
					for(u_int32_t readIndex : sharedReads_rev){
						cout << "\t\t" << readIndex << endl;
					}
				}
				*/

				u_int64_t nbSharedReads = computeSharedElements(unitigSource, uniqueSuccessor);


				if(unitigSource->_unitigIndex == 23062 || unitigSource->_unitigIndex == 23063){
					cout << "\t" << uniqueSuccessor->_unitigIndex << " " << nbSharedReads << endl;
				}

				if(nbSharedReads < _minSupportingReads) continue;
				uniqueSuccessorsSup.push_back({uniqueSuccessor, nbSharedReads});

			}

			std::sort(uniqueSuccessorsSup.begin(), uniqueSuccessorsSup.end(), UniqueSuccessorComparator);


			vector<bool> isSuccessors(uniqueSuccessorsSup.size(), false);

			for(size_t i=0; i<uniqueSuccessorsSup.size(); i++){
				if(isSuccessors[i]) continue;
				for(size_t j=0; j<uniqueSuccessorsSup.size(); j++){
					if(i==j) continue;
					if(isSuccessors[j]) continue;
					if(isSuccessor(uniqueSuccessorsSup[i]._node, uniqueSuccessorsSup[j]._node, unitigSource)){
						
						if(unitigSource->_unitigIndex == 23062 || unitigSource->_unitigIndex == 23063){
							cout << "\tSuccessor: " << uniqueSuccessorsSup[i]._node->_unitigIndex << " -> " << uniqueSuccessorsSup[j]._node->_unitigIndex << endl;
						}

						isSuccessors[j] = true;
					}
				}
			}

			
			vector<UnitSuccessor> uniqueSuccessorsSup_noSucc;

			for(size_t i=0; i<uniqueSuccessorsSup.size(); i++){
				if(isSuccessors[i]) continue;
				
				if(unitigSource->_unitigIndex == 23062 || unitigSource->_unitigIndex == 23063){
					cout << "\tFinal Successor: " << uniqueSuccessorsSup[i]._node->_unitigIndex << endl;
				}

				uniqueSuccessorsDerep.push_back(uniqueSuccessorsSup[i]._node);
				//uniqueSuccessorsSup_noSucc.push_back(uniqueSuccessorsSup[i]);
			}

			/*
			if(uniqueSuccessorsSup_noSucc.size() == 1){
				uniqueSuccessorsDerep.push_back(uniqueSuccessorsSup_noSucc[0]._node);
				return;
			}

			vector<UnitSuccessor> uniqueSuccessorsSup_noSucc_cutoff;

			if(uniqueSuccessorsSup_noSucc.size() > 1){

				for(size_t i=0; i<uniqueSuccessorsSup_noSucc.size(); i++){
					if(uniqueSuccessorsSup_noSucc[i]._nbSupportingReads < _minSupportingReads_2){
						uniqueSuccessorsSup_noSucc_cutoff.push_back(uniqueSuccessorsSup_noSucc[i]);
					}
				}
			}

			if(uniqueSuccessorsSup_noSucc_cutoff.size() == 1){
				uniqueSuccessorsDerep.push_back(uniqueSuccessorsSup_noSucc_cutoff[0]._node);
				return;
			}
			*/
		}

		/*
		u_int64_t computeSharedElements(const vector<ReadNodePos>& reads1, const vector<ReadNodePos>& reads2, bool useSuccessors){

			size_t i=0;
			size_t j=0;
			u_int64_t nbShared = 0;

			while(i < reads1.size() && j < reads2.size()){

				if(reads1[i]._readIndex == reads2[j]._readIndex){
					if(useSuccessors){
						if(reads2[j]._pos > reads1[i]._pos){
							nbShared += 1;
						}
					}
					else{
						if(reads2[j]._pos < reads1[i]._pos){
							nbShared += 1;
						}
					}
					//nbShared += 1;
					i += 1;
					j += 1;
				}
				else if(reads1[i]._readIndex < reads2[j]._readIndex){
					i += 1;
				}
				else{
					j += 1;
				}

			}

			return nbShared;
		}
		*/

		u_int64_t computeSharedElements(UnitigGraph::Node* unitig1, UnitigGraph::Node* unitig2){
		

			const vector<ReadNodePos>& source_readIndexes = _parent._nodeName_to_readIndexes[unitig1->endNode()];
			const vector<ReadNodePos>& source_readIndexes_rev = _parent._nodeName_to_readIndexes[_parent._progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(unitig1)->startNode()];

			const vector<ReadNodePos>& readIndexes = _parent._nodeName_to_readIndexes[unitig2->startNode()];
			const vector<ReadNodePos>& readIndexes_rev = _parent._nodeName_to_readIndexes[_parent._progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(unitig2)->endNode()];
			

			vector<u_int32_t> sharedReads;
			u_int64_t nbSharedReads = collectSharedElements(source_readIndexes, readIndexes, true, sharedReads);

			vector<u_int32_t> sharedReads_rev;
			nbSharedReads += collectSharedElements(source_readIndexes_rev, readIndexes_rev, false, sharedReads_rev);

			if(unitig1->_unitigIndex == 13685){
				for(u_int32_t readIndex : sharedReads){
					cout << "\t\tRead index: " << readIndex << endl;
				}
				for(u_int32_t readIndex : sharedReads_rev){
					cout << "\t\tRead index: " << readIndex << endl;
				}
			}

			return nbSharedReads;
		}

		
		u_int64_t collectSharedElements(const vector<ReadNodePos>& reads1, const vector<ReadNodePos>& reads2, bool useSuccessors, vector<u_int32_t>& sharedElement){

			sharedElement.clear();

			size_t i=0;
			size_t j=0;
			u_int64_t nbShared = 0;

			while(i < reads1.size() && j < reads2.size()){

				if(reads1[i]._readIndex == reads2[j]._readIndex){

					//cout << "lala" << endl;
					if(useSuccessors){
						if(reads2[j]._pos > reads1[i]._pos){
							nbShared += 1;
							sharedElement.push_back(reads1[i]._readIndex);
						}
					}
					else{
						if(reads2[j]._pos < reads1[i]._pos){
							nbShared += 1;
							sharedElement.push_back(reads1[i]._readIndex);
						}
					}
					//nbShared += 1;
					i += 1;
					j += 1;
				}
				else if(reads1[i]._readIndex < reads2[j]._readIndex){
					i += 1;
				}
				else{
					j += 1;
				}

			}

			return nbShared;
		}
		
		u_int64_t collectSharedElementsOld(UnitigGraph::Node* unitig1, UnitigGraph::Node* unitig2){

			//const vector<ReadNodePos>& source_readIndexes = ;
			//const vector<ReadNodePos>& source_readIndexes_rev = _parent._nodeName_to_readIndexes[_parent._progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(unitig1)->startNode()];

			//const vector<ReadNodePos>& readIndexes = _parent._nodeName_to_readIndexes[unitig2->startNode()];
			//const vector<ReadNodePos>& readIndexes_rev = _parent._nodeName_to_readIndexes[_parent._progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(unitig2)->endNode()];
			
			vector<u_int32_t> allReadIndex1;
			for(const ReadNodePos& r : _parent._nodeName_to_readIndexes[unitig1->endNode()]){
				allReadIndex1.push_back(r._readIndex);
			}
			for(const ReadNodePos& r : _parent._nodeName_to_readIndexes[GraphSimplify::nodeIndex_toReverseDirection(unitig1->endNode())]){
				allReadIndex1.push_back(r._readIndex);
			}

			vector<u_int32_t> allReadIndex2;
			for(const ReadNodePos& r : _parent._nodeName_to_readIndexes[unitig2->startNode()]){
				allReadIndex2.push_back(r._readIndex);
			}
			for(const ReadNodePos& r : _parent._nodeName_to_readIndexes[GraphSimplify::nodeIndex_toReverseDirection(unitig2->startNode())]){
				allReadIndex2.push_back(r._readIndex);
			}

			std::sort(allReadIndex1.begin(), allReadIndex1.end());
			std::sort(allReadIndex2.begin(), allReadIndex2.end());

			size_t i=0;
			size_t j=0;
			u_int64_t nbShared = 0;

			while(i < allReadIndex1.size() && j < allReadIndex2.size()){

				if(allReadIndex1[i] == allReadIndex2[j]){
					nbShared += 1;
					i += 1;
					j += 1;
				}
				else if(allReadIndex1[i] < allReadIndex2[j]){
					i += 1;
				}
				else{
					j += 1;
				}

			}

			return nbShared;
		}

        bool isPalindromeTip(const UnitigGraph::Node* node){

            //if(node == nullptr) return false;
            if(node->_unitigIndex == -1) return false;
            
            //if(node->_length > _maxLength) return false;
            
            //if(node->_unitigIndex % 2 == 1) return false;
            if(node->_successors.size() > 0) return false;
            //if(node->_predecessors.size() != 1) return false;
            if(node->_predecessors.size() == 0) return false;

            //if(node->_isPalindrome) return true;

            return node->_isPalindrome;
        }


		void solvePath(const StartingContig& startingContig){

			cout << "Solving path" << endl;

			vector<UnitigGraph::Node*> unitigPath;

			UnitigGraph::Node* nodeSource = getStartingUnitig(startingContig);

			if(nodeSource == nullptr){
				cout << "No starting unitig" << endl;
				return;
			}

			UnitigGraph::Node* nodeCurrent = nodeSource;
			
			while(true){

				unitigPath.push_back(nodeCurrent);

				vector<UnitigGraph::Node*> uniqueSuccessors;
				getNextUniqueSuccessors(nodeCurrent, uniqueSuccessors);

				vector<UnitigGraph::Node*> uniqueSuccessorsDerep;
				dereplicateUniqueSuccessors(nodeCurrent, uniqueSuccessors, uniqueSuccessorsDerep);

				cout << "Current unitig: " << nodeCurrent->_unitigIndex << endl;
				cout << "\tSuccessors:" << endl;
				for(UnitigGraph::Node* nodeSuccessor : uniqueSuccessorsDerep){
					cout << "\t\t" <<  nodeSuccessor->_unitigIndex << " " << computeSharedElements(nodeCurrent, nodeSuccessor) << endl;
				}

				if(uniqueSuccessorsDerep.size() == 0){
					cout << "No successors" << endl;
					break;
				}
				else if(uniqueSuccessorsDerep.size() > 1){
					cout << "Multiple successors" << endl;
					break;
				}

				nodeCurrent = uniqueSuccessorsDerep[0];
				//getchar();

				/*
				if(nodeCurrent == nodeSource){
					cout << "Path solved!" << endl;
					_parent._nbPathSolved += 1;

					for(UnitigGraph::Node* node: unitigPath){
						for(u_int32_t nodeIndex : node->_nodes){
							_parent._solvedNodeNames.insert(BiGraph::nodeIndex_to_nodeName(nodeIndex));
						}
					}
					break;
					//getchar();
				}
				*/
				if(_parent._isStartingUnitigs.find(nodeCurrent->_unitigIndex) != _parent._isStartingUnitigs.end()){
					
					cout << "Path solved!" << endl;
					_parent._nbPathSolved += 1;

					//for(UnitigGraph::Node* node: unitigPath){
					//	for(u_int32_t nodeIndex : node->_nodes){
					//		_parent._solvedNodeNames.insert(BiGraph::nodeIndex_to_nodeName(nodeIndex));
					//	}
					//}
					break;
					//getchar();

				}
				
			}

		}

		UnitigGraph::Node* getStartingUnitig(const StartingContig& startingContig){

			UnitigGraph::Node* nodeSource = nullptr;
			u_int64_t longerUniqueUnitig = 0;

			for(size_t i=0; i<startingContig._nodePath.size(); i++){
				u_int32_t nodeIndex = startingContig._nodePath[i];
				u_int32_t unitigIndex = _parent._nodeIndex_to_unitigIndex[nodeIndex];
				UnitigGraph::Node* unitig = _parent._progressiveAbundanceFilter->_unitigGraph->_nodes[unitigIndex];

				if(_isUnitigUnique[unitigIndex]){
					nodeSource = unitig;
					break;
					//if(unitig->_length > longerUniqueUnitig){
					//	longerUniqueUnitig = unitig->_length;
					//	nodeSource = unitig;
					//}
				}
				
			}

			if(nodeSource == nullptr) return nodeSource;

			vector<UnitigGraph::Node*> uniqueSuccessors;
			getNextUniqueSuccessors(nodeSource, uniqueSuccessors);

			vector<UnitigGraph::Node*> uniqueSuccessorsDerep;
			dereplicateUniqueSuccessors(nodeSource, uniqueSuccessors, uniqueSuccessorsDerep);

			if(uniqueSuccessorsDerep.size() != 1) return nullptr;

			if(_parent._isStartingUnitigs.find(uniqueSuccessorsDerep[0]->_unitigIndex) != _parent._isStartingUnitigs.end()){
				return _parent._progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(nodeSource);
			}
			else{
				return nodeSource;
			}
			
		}

	};

};





#endif
