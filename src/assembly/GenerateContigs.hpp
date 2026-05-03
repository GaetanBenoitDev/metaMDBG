

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
	//string _contigGraphDir;
	//string _truthInputFilename;
	//string _outputFilename;
	//string _outputFilename_complete;
	//bool _debug;
	//bool _isFinalAssembly;
	//string _inputFilename_unitigNt;
	//string _inputFilename_unitigCluster;
	//string _filename_abundance;

	Parameters _params;
	/*
	float _minimizerDensity;
    size_t _minimizerSize;
    size_t _kminmerSize;
    size_t _kminmerSizeFirst;
	float _minimizerSpacingMean;
	float _kminmerLengthMean;
	float _kminmerOverlapMean;
	*/

	//string _partitionDir;
	//string _filename_solidKminmers;
	//string _filename_inputContigs;
	//string _filename_outputContigs;
	//string _filename_readMinimizers;
	//gzFile _outputContigFile;
	//gzFile _outputContigFile_complete;
	//MDBG* _mdbg;
	bool _genGraph;
	//GraphSimplify* _graph;
	//ToBasespaceOnTheFly _toBasespace;
	//ContigFeature _contigFeature;
	//string _gfaFilename;
	size_t _nbCores;

    //vector<u_int32_t> _nodeName_to_abundance;
	int _maxBubbleLength;
	int _maxTipLength;
	bool _isFinalAssembly;

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
		args::ValueFlag<int> arg_maxBubbleLength(parser, "", "Max bubble popping length", {ARG_MAX_BUBBLE_LENGTH}, 50000);
		args::ValueFlag<int> arg_maxTipLength(parser, "", "Max tip clipping length", {ARG_MAX_TIP_LENGTH}, 50000);
		args::ValueFlag<int> arg_nbCores(parser, "", "Number of cores", {ARG_NB_CORES2}, NB_CORES_DEFAULT_INT);
		
		//args::Flag arg_isFinalAssembly(parser, "", "Is final multi-k pass", {ARG_FINAL});
		//args::Flag arg_firstPass(parser, "", "Is first pass of multi-k", {ARG_FIRST_PASS});
		args::Flag arg_genGraph(parser, "", "Generate temporary file for genetarating assembly graphs", {"gen-graph"});
		args::Flag arg_finalAssembly(parser, "", "Solve final graph", {"final"});
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


		_genGraph = false;
		if(arg_genGraph){
			_genGraph = true;
		}

		//_contigGraphDir = _inputDir;

		_isFinalAssembly = false;
		if(arg_finalAssembly){
			_isFinalAssembly = true;
			//_contigGraphDir = _inputDir + "/contigGraph/";
		}

		_maxBubbleLength = args::get(arg_maxBubbleLength);
		_maxTipLength = args::get(arg_maxTipLength);


		if(_isFinalAssembly){
			const string contigGraphDir = _inputDir + "/contigGraph/";
			fs::copy(contigGraphDir + "/parameters.gz", _inputDir + "/parameters.gz", fs::copy_options::overwrite_existing);
			fs::copy(contigGraphDir + "/unitigGraph.edges.successors.bin", _inputDir + "/unitigGraph.edges.successors.bin", fs::copy_options::overwrite_existing);
			fs::copy(contigGraphDir + "/unitigGraph.nodes.abundances.bin", _inputDir + "/unitigGraph.nodes.abundances.bin", fs::copy_options::overwrite_existing);
			fs::copy(contigGraphDir + "/unitigGraph.nodes.bin", _inputDir + "/unitigGraph.nodes.bin", fs::copy_options::overwrite_existing);
			fs::copy(contigGraphDir + "/unitigGraph.stats.bin", _inputDir + "/unitigGraph.stats.bin", fs::copy_options::overwrite_existing);
			fs::copy(contigGraphDir + "/contig_data_final.bin", _inputDir + "/contig_data_final.bin", fs::copy_options::overwrite_existing);

		}

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


		_params.load(_inputDir + "/parameters.gz");
		/*
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
		*/

		openLogFile(_inputDir);

		/*
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
		_filename_outputContigs = _inputDir + "/contigs.min.gz";
		_filename_solidKminmers = _inputDir + "/solid.min.gz";
		*/

	}

	//unordered_map<u_int64_t, vector<u_int32_t>> _debug_readPath;
    //vector<bool> _isBubble;

	class X
	{
		public: 
		void operator()(float cutoff) {}
	};

	void execute(){


		/*
		ifstream inputFile2(_inputDir + "/kminmerData_min.txt");

		while (true) {

			u_int32_t nodeName;
			u_int32_t abundance;
			vector<MinimizerType> minimizers;

			bool iseof = MDBG::readKminmer(nodeName, abundance, minimizers, inputFile2, _kminmerSize);

			if(iseof) break;

			KmerVec vec;
			vec._kmers = minimizers;
			vec = vec.normalize();

			_kmerVec_to_nodeName[vec] = nodeName;
		}

		inputFile2.close();
		*/
		//cout << "Loading abundant nodes" << endl;
		//loadAbundantNodes();
		
		loadGraph();

		//if(_params._kminmerSize == 81){
		//	mapReadsToGraph();
		//}
		//exit(1);
		//generateUnitigs();
		//generateContigs2(_inputDir + "/contigs.nodepath");



		if(_isFinalAssembly){

			generateContigs3();
			_progressiveAbundanceFilter->_repeatSolver->writeUnitigSequences();
		}
		else{
			delete _progressiveAbundanceFilter->_unitigGraph2;

			generateContigs3();

			dumpUnitigAbundances();

		}
		//exit(1);
		
		/*
		u_int64_t checksum = 0;
		u_int64_t nbNodes = 0;
		
		unordered_map<KmerVec, float> mdbg;

		ifstream inputFile(_inputDir + "/kminmerData_min.txt");

		while (true) {

			u_int32_t nodeName;
			u_int32_t abundance;
			vector<MinimizerType> minimizers;

			bool iseof = MDBG::readKminmer(nodeName, abundance, minimizers, inputFile, _kminmerSize);

			if(iseof) break;

			KmerVec vec;
			vec._kmers = minimizers;
			vec = vec.normalize();

			checksum += abundance;

			mdbg[vec] = abundance;
			//if(_nodeNameAbundances.find(vec) == _nodeNameAbundances.end()){
			//	_nodeNameAbundances[vec] = {abundance, 1};
			//}
		}

		inputFile.close();
		
		//out << "Checkum init: " << checksum << endl;

		checksum = 0;

		ofstream kminmerFile = ofstream(_inputDir + "/kminmerData_min.txt");

		for(auto& it : mdbg){

			const KmerVec& vec = it.first;

			if(_nodeNameAbundances.find(vec) != _nodeNameAbundances.end()){
				const NodeAb& nodeAb = _nodeNameAbundances[vec];
				it.second = ceil(nodeAb._abundance);

				//"reprise: essayer d'enelever pas ça pas sur que c'est correct en fait"
				if(nodeAb._nbNodes-_kminmerSize+1 > _kminmerSize) it.second = max(it.second, (float)2);
			}

			MDBG::writeKminmer(0, it.second, vec._kmers, kminmerFile);

		}
		
		
		kminmerFile.close();
		
		
		for(auto& it : mdbg){
			checksum += it.second;
			nbNodes += 1;
		}
		cout << "Checksum abundance: " << nbNodes << " " << checksum << endl;
		*/
		//getchar();

		/*
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
		*/

		//closeLogFile();
	}



	void loadGraph(){
		/*
		KmerVec vec;
		vec._kmers = {706036053, 218005452, 609864957};
		printf("unsinged int 128 bit: %llu\n", (unsigned long long)vec.hash128());
		printf("unsinged int 128 bit: %llu\n", (unsigned long long)vec.normalize().hash128());
		vec._kmers = {218005452, 609864957, 218005452};
		printf("unsinged int 128 bit: %llu\n", (unsigned long long)vec.hash128());
		printf("unsinged int 128 bit: %llu\n", (unsigned long long)vec.normalize().hash128());

		
		vec._kmers = {494401525, 456413655, 220428157};
		printf("unsinged int 128 bit: %llu\n", (unsigned long long)vec.hash128());
		printf("unsinged int 128 bit: %llu\n", (unsigned long long)vec.normalize().hash128());
		vec._kmers = {456413655, 220428157, 460823664};
		printf("unsinged int 128 bit: %llu\n", (unsigned long long)vec.hash128());
		printf("unsinged int 128 bit: %llu\n", (unsigned long long)vec.normalize().hash128());
		vec._kmers = {220428157, 460823664, 899439573};
		printf("unsinged int 128 bit: %llu\n", (unsigned long long)vec.hash128());
		printf("unsinged int 128 bit: %llu\n", (unsigned long long)vec.normalize().hash128());
		vec._kmers = {460823664, 899439573, 225835357};
		printf("unsinged int 128 bit: %llu\n", (unsigned long long)vec.hash128());
		printf("unsinged int 128 bit: %llu\n", (unsigned long long)vec.normalize().hash128());
		*/


		//_gfaFilename = _inputDir + "/minimizer_graph.gfa";
		//string gfa_filename_noUnsupportedEdges = _inputDir + "/minimizer_graph_noUnsupportedEdges.gfa";
		//string gfa_filename_unitigs = _inputDir + "/minimizer_graph_unitigs.gfa";
		//string mdbg_filename = _inputDir + "/mdbg_nodes.gz";



		Logger::get().debug() << "Loading unitig graph";
		
		UnitigGraph2* unitigGraph2 = new UnitigGraph2(_params._kminmerSize, _params._kminmerLengthMean, _params._kminmerOverlapMean, _params._kminmerLengthMean-_params._kminmerOverlapMean, _params._minimizerSpacingMean, _nbCores);
		unitigGraph2->load(_inputDir + "/unitigGraph.stats.bin", _inputDir + "/unitigGraph.nodes.bin", _inputDir + "/unitigGraph.nodes.abundances.bin", _inputDir + "/unitigGraph.edges.successors.bin");


		//Logger::get().debug() << "Loading unitig graph";
		//UnitigGraph* unitigGraph = nullptr; //new UnitigGraph(_kminmerSize, _kminmerLengthMean, _kminmerOverlapMean, _kminmerLengthMean-_kminmerOverlapMean, _kmerVec_to_nodeName);
		//unitigGraph->load(_inputDir + "/unitigGraph.nodes.bin", _inputDir + "/unitigGraph.nodes.abundances.bin", _inputDir + "/unitigGraph.edges.successors.bin", _inputDir + "/unitigGraph.edges.predecessors.bin", _nodeName_to_abundance);

        //Logger::get().debug() << "Nb unitigs: " << unitigGraph->_nodes.size();
        //_logFile << "Nb edges: " << nbEdgeSuccessors<< endl;

		Logger::get().debug() << "Simplifying";

        _progressiveAbundanceFilter = new ProgressiveAbundanceFilter(unitigGraph2, _inputDir, _params._kminmerSize, true, _params._kminmerSize==_params._kminmerSizeFirst, _genGraph, _maxBubbleLength, _maxTipLength, _nbCores, _isFinalAssembly, _params);
		X dummy;

		//delete _graph->_graphSuccessors;

        _progressiveAbundanceFilter->execute(dummy);


		//cout << _graph->_unitigGraph->computeChecksum_successor() << endl;
		//if(_kminmerSize == 20) exit(1);


        //_unitigGraph->save("/home/gats/workspace/run/unitigGraph_cleaned.gfa");

		
		//_graph->debug_selectUnitigIndex();
		//_graph->debug_writeGfaErrorfree(2000, 2000, -1, _kminmerSize, false, true, false, _unitigDatas, true, false, false, false, false, false, _mdbg, _minimizerSize, _nbCores, false, true);
			
		//delete _mdbg;
		
		//_processedNodeNames = vector<bool>(nbNodes/2, false);
		
	}

	ProgressiveAbundanceFilter* _progressiveAbundanceFilter;



	bool isContigAssembled(const vector<UnitigType>& nodePath){

		for(const UnitigType& unitigIndex : nodePath){

			const UnitigType& unitigName = UnitigGraph2::unitigIndex_to_unitigName(unitigIndex);

			if(_processedNodeNames.find(unitigName) != _processedNodeNames.end()) return true;
		}

		return false;

	}




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

	unordered_set<UnitigType> _processedNodeNames;
	//vector<bool> _processedNodeNames;
	unordered_map<UnitigType, NodeAb> _nodeNameAbundances;



	void generateContigs3(){


		u_int64_t checkumAbundance = 0;
		Logger::get().debug() << "Generating contigs";
		//size_t nextMultikStep = Commons::getMultikStep(_kminmerSize);


		ofstream outputContigFile(_inputDir + "/contigs.nodepath");

		float prevCutoff = -1;
		u_int64_t contigIndex = 0;


		u_int64_t checksumTotal = 0;
		//unordered_set<u_int32_t> processedNodeNames;


		for(long i=_progressiveAbundanceFilter->_cutoffIndexes.size()-1; i>=0; i--){

			ProgressiveAbundanceFilter::CutoffIndex cutoffIndex = _progressiveAbundanceFilter->_cutoffIndexes[i];
			float cutoff = cutoffIndex._cutoffValue;
			Logger::get().debug() << cutoff;

        	ifstream nodepathFile(_inputDir + "/filter/unitigs_" + to_string(cutoffIndex._cutoffIndex) + ".bin");

			_minUnitigAbundance = cutoff / 0.5;
			//_minUnitigAbundance = min(_minUnitigAbundance, cutoff+30); //cuttof + 30x

			while(true){

				u_int8_t isCircular = CONTIG_LINEAR;
				bool isRepeatSide = false;
				float contigAbundance = -1;
				u_int32_t nbMinimizers = 0;
				vector<UnitigType> nodePath;
				u_int32_t size;
				nodepathFile.read((char*)&size, sizeof(size));
				
				if(nodepathFile.eof()) break;


				nodepathFile.read((char*)&isCircular, sizeof(isCircular));
				nodepathFile.read((char*)&isRepeatSide, sizeof(isRepeatSide));
				nodepathFile.read((char*)&contigAbundance, sizeof(contigAbundance));
				nodepathFile.read((char*)&nbMinimizers, sizeof(nbMinimizers));

				nodePath.resize(size);
				nodepathFile.read((char*)&nodePath[0], size * sizeof(UnitigType));

				//cout << isRepeatSide << " " << (contigAbundance < _minUnitigAbundance) << " " << isContigAssembled(nodePath) << endl;

				if(isCircular){
					//cout << "Discard circular 2" << endl;
					//continue;
				}




				//if(isRepeatSide) continue;
				if(contigAbundance < _minUnitigAbundance) continue;
				/*
				for(size_t i=0; i<nodePath.size(); i++){
					const UnitigType unitigName = UnitigGraph2::unitigIndex_to_unitigName(nodePath[i]);
					if(unitigName == 291){
						cout << "omg: " << contigAbundance << " " << nbMinimizers << " " << nodePath.size() << endl;
						
						for(const UnitigType& unitigIndex : nodePath){

							const UnitigType& unitigName = UnitigGraph2::unitigIndex_to_unitigName(unitigIndex);
							cout << "utg" << unitigName << "\t" << (_processedNodeNames.find(unitigName) != _processedNodeNames.end()) << "\t" << (_progressiveAbundanceFilter->_repeatSolver->_repeatedUnitigNames.find(unitigName) != _progressiveAbundanceFilter->_repeatSolver->_repeatedUnitigNames.end()) << endl;
						}
						getchar();
					}
				}
				*/

				if(isContigAssembled(nodePath)) continue;

				checkumAbundance += (contigAbundance*(nbMinimizers-_params._kminmerSize+1));
				//cout << contigIndex << " " << (int)isCircular << " " << (nbMinimizers-_kminmerSize+1) << " " << checkumAbundance << endl;
				//u_int64_t checksum = 0;
				//for(MinimizerType minimizer : nodePath){
				//	checksum += minimizer;
				//}
				//checksum *= nodePath.size();
				//checksumTotal += checksum;


				//if(isCircular){
				//	cout << "lol" << endl;
				//}
				if(isCircular && nbMinimizers-_params._kminmerSize+1 > 1){ //&& (nodePath.size()-_kminmerSize+1) > 1

					nbMinimizers += 1;


					//for(MinimizerType m : nodePath){
					//	cout << m << endl;
					//}

					//getchar();
					//for(size_t i=0; i<_kminmerSize; i++){
					//	nodePath.push_back(nodePath[i]);
					//}
					//nodePath.push_back(nodePath[0]);

					//cout << "plop: " << contigAbundance << endl;
					//for(u_int32_t unitigIndex : nodePath){
					//	cout << UnitigGraph2::unitigIndexToString(unitigIndex) << endl;
					//}
					//getchar();

				}

				size = nodePath.size();
				//cout << size << "  " << nodePath.size() << " " << contigAbundance << endl;
				outputContigFile.write((const char*)&size, sizeof(size));
				outputContigFile.write((const char*)&isCircular, sizeof(isCircular));
				outputContigFile.write((const char*)&nodePath[0], size * sizeof(UnitigType));



				/*
				for(size_t i=0; i<nodePath.size(); i++){
					const UnitigType unitigName = UnitigGraph2::unitigIndex_to_unitigName(nodePath[i]);
					if(unitigName == 291){
						cout << "omg: " << contigAbundance << " " << nbMinimizers << " " << nodePath.size() << endl;
						
						for(const UnitigType& unitigIndex : nodePath){

							const UnitigType& unitigName = UnitigGraph2::unitigIndex_to_unitigName(unitigIndex);
							cout << "utg" << unitigName << "\t" << (_processedNodeNames.find(unitigName) != _processedNodeNames.end()) << "\t" << (_progressiveAbundanceFilter->_repeatSolver->_repeatedUnitigNames.find(unitigName) != _progressiveAbundanceFilter->_repeatSolver->_repeatedUnitigNames.end()) << endl;
						}
					}
				}

				if(contigAbundance > 50){
					cout << "ctg" << contigIndex << endl;
					for(size_t i=0; i<nodePath.size(); i++){
						if(_progressiveAbundanceFilter->_repeatSolver->_unitigName_to_referenceName_init.find(nodePath[i]) != _progressiveAbundanceFilter->_repeatSolver->_unitigName_to_referenceName_init.end()){
							if(_progressiveAbundanceFilter->_repeatSolver->_unitigName_to_referenceName_init[nodePath[i]][0]._referenceName == "Prevotella_corporis.final.genome2"){
								cout << "\t" << i << "\t" << "utg" << UnitigGraph2::unitigIndex_to_unitigName(nodePath[i]) << "\t" << _progressiveAbundanceFilter->_repeatSolver->_unitigName_to_referenceName_init[nodePath[i]][0]._referenceName << endl;
							}
						}
					}
				}
				*/
				//cout << cutoff << " " << _minUnitigAbundance << " " << nbMinimizers << endl;
				//vector<MinimizerType> minimizers = _progressiveAbundanceFilter->_unitigGraph2->unitigSequenceToMinimizerSequence(nodePath);
				//cout << minimizers.size() << endl;
				
				//cout << "Output: " << _minUnitigAbundance << " " << contigAbundance << endl;
				//for(UnitigType unitigIndex : nodePath){
				//	cout << "\t" << UnitigGraph2::unitigIndexToString(unitigIndex) << endl;
				//}
				//if(size > 100)
				//	cout << "contig: " << size << endl;

				
				//u_int32_t nodeName1 = BiGraph::nodeIndex_to_nodeName(nodePath[0]);
				//u_int32_t nodeName2 = BiGraph::nodeIndex_to_nodeName(nodePath[nodePath.size()-1]);
				/*
				for(const KmerVec& kmerVec : Utils::minimizersToKminmers(nodePath, _kminmerSize)){
				//for(MinimizerType minimizer : nodePath){
					//u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
					//cout << nodeName << endl;
					_processedNodeNames.insert(kmerVec.normalize());
					_nodeNameAbundances[kmerVec.normalize()] = {contigAbundance, (u_int32_t)nodePath.size()};


				}
				*/


				for(const UnitigType& unitigIndex : nodePath){

					const UnitigType& unitigName = UnitigGraph2::unitigIndex_to_unitigName(unitigIndex);

					if(_isFinalAssembly && _progressiveAbundanceFilter->_repeatSolver->_repeatedUnitigNames.find(unitigName) != _progressiveAbundanceFilter->_repeatSolver->_repeatedUnitigNames.end()) continue; //Repeated unitigs from solved repeats can be present in multiple contigs

					_processedNodeNames.insert(unitigName);
					_nodeNameAbundances[unitigName] = {contigAbundance, nbMinimizers};
				}

				//if(contigIndex == 30) getchar();

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

		//cout << "Checksum 2: " << checkumAbundance << endl;
		//getchar();
		
		//if(_params._kminmerSize == 21) getchar();
	}

	void dumpUnitigAbundances(){


		//if(_params._kminmerSize == _params._kminmerSizeFirst){
			
		ofstream outputFile(_inputDir + "/unitigGraph.nodes.refined_abundances.bin");

		for(const auto& it : _nodeNameAbundances){

			const UnitigType& unitigName = it.first;
			const NodeAb& nodeAb = it.second;

			u_int32_t abundance = ceil(nodeAb._abundance);
			if(nodeAb._nbNodes-_params._kminmerSize+1 > _params._kminmerSize) abundance = max(abundance, (u_int32_t)2);
			//float abundance = nodeAb._abundance;

			outputFile.write((const char*)&unitigName, sizeof(unitigName));
			outputFile.write((const char*)&abundance, sizeof(abundance));
		}


		outputFile.close();

		//}

		fs::copy(_inputDir + "/unitigGraph.nodes.bin", _inputDir + "/unitigGraph_prev.nodes.bin", fs::copy_options::overwrite_existing);
		fs::copy(_inputDir + "/kminmerData_abundance.txt", _inputDir + "/kminmerData_abundance_prev.txt", fs::copy_options::overwrite_existing);


		if(_params._kminmerSize > _params._kminmerSizeFirst){ //Backup full unitig graph as soon as we start multiplex steps
			fs::copy(_inputDir + "/unitigGraph.edges.successors.bin", _inputDir + "/unitigGraph_prev.edges.successors.bin", fs::copy_options::overwrite_existing);
			fs::copy(_inputDir + "/unitigGraph.nodes.abundances.bin", _inputDir + "/unitigGraph_prev.nodes.abundances.bin", fs::copy_options::overwrite_existing);
			fs::copy(_inputDir + "/unitigGraph.stats.bin", _inputDir + "/unitigGraph_prev.stats.bin", fs::copy_options::overwrite_existing);
		}
		
		if(_params._kminmerSize == 21){ //contig graph

			string dir = _inputDir + "/contigGraph/";

			if(!fs::exists (dir)){
				fs::create_directories(dir); 
			} 

			fs::copy(_inputDir + "/parameters.gz", dir + "/parameters.gz", fs::copy_options::overwrite_existing);
			fs::copy(_inputDir + "/unitigGraph.nodes.bin", dir + "/unitigGraph.nodes.bin", fs::copy_options::overwrite_existing);
			fs::copy(_inputDir + "/unitigGraph.edges.successors.bin", dir + "/unitigGraph.edges.successors.bin", fs::copy_options::overwrite_existing);
			fs::copy(_inputDir + "/unitigGraph.nodes.abundances.bin", dir + "/unitigGraph.nodes.abundances.bin", fs::copy_options::overwrite_existing);
			fs::copy(_inputDir + "/unitigGraph.stats.bin", dir + "/unitigGraph.stats.bin", fs::copy_options::overwrite_existing);
		
		}
		
		

		//if(_kminmerSize == 4){
			//cout << "GenerateContigs: copy initial unitig abundance" << endl;
			
		//	fs::copy(_inputDir + "/unitigGraph.nodes.refined_abundances.bin", _inputDir + "/unitigGraph.nodes.refined_abundances.bin.init", fs::copy_options::overwrite_existing);
		//}
		/*
		_mdbg = new MDBG(_kminmerSize);
		_mdbg->load(_inputDir + "/kminmerData_min.txt", false);
		for(auto& it : _mdbg->_dbg_nodes){
			if(_nodeNameAbundances.find(it.second._index) == _nodeNameAbundances.end()) continue;
			const NodeAb& nodeAb = _nodeNameAbundances[it.second._index];
			it.second._abundance = ceil(nodeAb._abundance);

			if(nodeAb._nbNodes > _kminmerSize) it.second._abundance = max(it.second._abundance, (u_int32_t)2);

		}
		_mdbg->dump(_inputDir + "/kminmerData_min.txt");
		*/

		/*
		ifstream abundanceFile(abundanceFilename);

		while(true){


			u_int32_t unitigIndex;
			abundanceFile.read((char*)&unitigIndex, sizeof(unitigIndex));

			u_int32_t nbAbundances;
			abundanceFile.read((char*)&nbAbundances, sizeof(nbAbundances));
			

			if(abundanceFile.eof()) break;

			vector<float> abundances;
			abundances.resize(nbAbundances);
			abundanceFile.read((char*)&abundances[0], nbAbundances * sizeof(float));

            _nodes[unitigIndex]->_abundances = abundances;
            _nodes[unitigIndex]->_abundance = _nodes[unitigIndex]->computeMedianAbundance();

            u_int32_t unitigIndexRC = unitigIndex_toReverseDirection(unitigIndex);
            _nodes[unitigIndexRC]->_abundances = abundances;
            _nodes[unitigIndexRC]->_abundance = _nodes[unitigIndexRC]->computeMedianAbundance();

        }

        abundanceFile.close();









		_kminmerAbundances.clear();
		ifstream kminmerFile(_outputDir + "/kminmerData_min.txt");

		bool isEOF = false;
		KmerVec vec;
		u_int32_t nodeName;


		while (true) {


			vector<MinimizerType> minimizerSeq;
			minimizerSeq.resize(_kminmerSize);
			kminmerFile.read((char*)&minimizerSeq[0], minimizerSeq.size()*sizeof(MinimizerType));

			isEOF = kminmerFile.eof();
			if(isEOF) break;

			
			u_int32_t abundance;
			//u_int32_t quality;

			kminmerFile.read((char*)&nodeName, sizeof(nodeName));
			kminmerFile.read((char*)&abundance, sizeof(abundance));
			//kminmerFile.read((char*)&quality, sizeof(quality));
			vec._kmers = minimizerSeq;

			_kminmerAbundances[vec.hash128()] = abundance;

		}

		kminmerFile.close();


		*/
	}

	/*
	size_t _kminmerSizeMapping;
	ankerl::unordered_dense::map<u_int128_t, vector<u_int32_t>> _kminmer_to_contigIndexes;
	ankerl::unordered_dense::map<u_int128_t, u_int32_t> _unitigName_to_nbKminmers;

	void mapReadsToGraph(){

		_kminmerSizeMapping = 21;
		loadContigs();

		cout << "single core here" << endl;
		KminmerParserParallel parser(_inputDir + "/read_data_corrected.txt", 0, _kminmerSizeMapping, false, false, 1);
		//KminmerParserParallel parser(_inputDir + "/read_data_init.txt", 0, _kminmerSizeMapping, false, true, 1);
		parser.parse(ReadMapperFunctor(*this));
	}

	void loadContigs(){

		ifstream nodeFile(_inputDir + "/unitigGraph.nodes.bin");

		while(true){

			u_int32_t size;
			nodeFile.read((char*)&size, sizeof(size));
			

			if(nodeFile.eof()) break;

			vector<MinimizerType> minimizers;
			minimizers.resize(size);
			nodeFile.read((char*)&minimizers[0], size * sizeof(MinimizerType));

			UnitigType unitigIndex;
			nodeFile.read((char*)&unitigIndex, sizeof(unitigIndex));

			const UnitigType unitigName = unitigIndex/2;

			const vector<KmerVec>& kminmers = MDBG::minimizersToKminmers(minimizers, _kminmerSizeMapping);

			for(const KmerVec& kminmer : kminmers){

				const u_int128_t vecHash = kminmer.normalize().hash128();

				if(_kminmer_to_contigIndexes.find(vecHash) == _kminmer_to_contigIndexes.end()){
					_kminmer_to_contigIndexes[vecHash].push_back(unitigName);
				}
				else{
					auto& it = _kminmer_to_contigIndexes[vecHash];
					if(std::find(it.begin(), it.end(), unitigName) == it.end()){
						it.push_back(unitigName);
					}
				}
			}

			_unitigName_to_nbKminmers[unitigName] = kminmers.size();
			
			if(unitigName == 525 || unitigName == 1613){
				cout << "utg" << unitigName << endl;
				for(const auto& m : minimizers){
					cout << "\t" << m << endl;
				}
			}
            
        }

        nodeFile.close();
		

	}


	class ReadMapperFunctor {

		public:

		GenerateContigs& _parent;

		ReadMapperFunctor(GenerateContigs& parent) : _parent(parent){
		}

		ReadMapperFunctor(const ReadMapperFunctor& copy) : _parent(copy._parent){
			
		}

		~ReadMapperFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {

			u_int64_t readIndex = kminmerList._readIndex;
			const vector<MinimizerType>& readMinimizers = kminmerList._readMinimizers;
			//const vector<KmerVec>& kminmers = kminmerList._kminmers;
			const vector<ReadKminmerComplete>& kminmersInfos = kminmerList._kminmersInfo;

			unordered_map<UnitigType, u_int32_t> unitigName_to_hits;

			for(size_t i=0; i<kminmersInfos.size(); i++){
				
				const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];
				const KmerVec& vec = kminmerInfo._vec;
				
				const u_int128_t vecHash = vec.normalize().hash128();

				if(_parent._kminmer_to_contigIndexes.find(vecHash) == _parent._kminmer_to_contigIndexes.end()) continue;

				for(const UnitigType& unitigName : _parent._kminmer_to_contigIndexes[vecHash]){
					unitigName_to_hits[unitigName] += 1;
				}

			}


			cout << "----" << endl;
			cout << "read " << kminmerList._readIndex << endl;
			for(const auto& m : readMinimizers){
				cout << "\t" << m << endl;
			}

			cout << "Nb kminmers: " << kminmersInfos.size() << endl;

			ofstream colorFile = ofstream("/pasteur/appa/homes/gbenoit/zeus/tmp/readHits.csv");
			colorFile << "Name,Hits" << endl;

			for(const auto& it : unitigName_to_hits){

				float score = it.second; // / ((long double) _parent._unitigName_to_nbKminmers[it.first]);

				colorFile << "utg" << it.first << "," << score << endl;
				cout << "\t" << "utg" << it.first << "\t" << it.second << endl;
			}

			colorFile.close();
			cout << "Mapped read: " << kminmerList._readIndex << endl;
			getchar();

		}
	};
	*/
};





#endif
