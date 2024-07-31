

#ifndef MDBG_METAG_RELOCATIONREMOVER
#define MDBG_METAG_RELOCATIONREMOVER


#include "Commons.hpp"

//#include "graph/Graph.hpp"
#include "graph/GfaParser.hpp"
#include "graph/GraphSimplify.hpp"
//#include "eval/ContigStatistics.hpp"
//#include "toBasespace/ToBasespaceOnTheFly.hpp"
//#include "contigFeatures/ContigFeature.hpp"










class RelocationRemover : public Tool{

public:

	string _inputDir;
	string _truthInputFilename;
	string _outputFilename;
	string _outputFilename_complete;
	bool _debug;
	bool _isFinalAssembly;
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
	vector<UnitigData> _unitigDatas;

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
	GraphSimplify* _graph;
	string _inputFilenameContig;
	//ToBasespaceOnTheFly _toBasespace;
	//ContigFeature _contigFeature;
	string _gfaFilename;
	size_t _nbCores;

	RelocationRemover(): Tool (){

	}

	~RelocationRemover(){

	}


	void parseArgs(int argc, char* argv[]){

		args::ArgumentParser parser("reloc", ""); //"This is a test program.", "This goes after the options."
		args::Positional<std::string> arg_inputContigFilename(parser, "inputContigFilename", "", args::Options::Required);
		args::Positional<std::string> arg_outputDir(parser, "outputDir", "", args::Options::Required);
		//args::Positional<std::string> arg_contigs(parser, "contigs", "", args::Options::Required);
		//args::PositionalList<std::string> arg_readFilenames(parser, "reads", "Input filename(s) (separated by space)", args::Options::Required);
		//args::ValueFlag<int> arg_l(parser, "", "Minimizer length", {ARG_MINIMIZER_LENGTH2}, 13);
		//args::ValueFlag<float> arg_d(parser, "", "Minimizer density", {ARG_MINIMIZER_DENSITY2}, 0.005f);
		//args::ValueFlag<std::string> arg_contigs(parser, "", "", {ARG_INPUT_FILENAME_CONTIG}, "");
		args::ValueFlag<int> arg_nbCores(parser, "", "Number of cores", {ARG_NB_CORES2}, NB_CORES_DEFAULT_INT);
		//args::Flag arg_isFinalAssembly(parser, "", "Is final multi-k pass", {ARG_FINAL});
		//args::Flag arg_firstPass(parser, "", "Is first pass of multi-k", {ARG_FIRST_PASS});
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

		_inputFilenameContig = args::get(arg_inputContigFilename);


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

	unordered_map<u_int64_t, vector<u_int32_t>> _debug_readPath;
    vector<bool> _isBubble;

	void execute (){


		loadGraph();
		processContigs();

		closeLogFile();
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
		
		GraphSimplify* graphSimplify = new GraphSimplify(_gfaFilename, _inputDir, 0, _kminmerSize, _nbCores, _kminmerLengthMean, _kminmerOverlapMean, _logFile);
		_graph = graphSimplify;
		

		ifstream kminmerFile(_inputDir + "/kminmerData_min.txt");
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


		//if(_kminmerSize > 45){
		//	getchar();
		//}
		//	cout << "lala" << endl;
		//	getchar();
		//}

	  //_graph->debug_writeGfaErrorfree(0, 0, -1, _kminmerSize, false, true, false, _unitigDatas, true, false, true, false, true, true, _mdbg, _minimizerSize, _nbCores, true, false);
		//_graph->debug_writeGfaErrorfree(0, 0, -1, _kminmerSize, false, true, false, _unitigDatas, true, false, false, false, false, true, _mdbg, _minimizerSize, _nbCores, false, false);
		
		//_graph->debug_writeGfaErrorfree(2000, 2000, -1, _kminmerSize, false, true, false, _unitigDatas, true, false, false, false, false, true, _mdbg, _minimizerSize, _nbCores, false, false);
		
		//_graph->_unitigGraph->save("/home/gats/workspace/run/minimizer_graph_u.gfa");
		//getchar();
		
        //_progressiveAbundanceFilter = new ProgressiveAbundanceFilter(_graph->_unitigGraph, _inputDir, _kminmerSize);
        //_progressiveAbundanceFilter->execute();

		//cout << _graph->_unitigGraph->computeChecksum_successor() << endl;
		//if(_kminmerSize == 20) exit(1);


        //_unitigGraph->save("/home/gats/workspace/run/unitigGraph_cleaned.gfa");

		
		//_graph->debug_selectUnitigIndex();
		//_graph->debug_writeGfaErrorfree(2000, 2000, -1, _kminmerSize, false, true, false, _unitigDatas, true, false, false, false, false, false, _mdbg, _minimizerSize, _nbCores, false, true);
			
		//delete _mdbg;
		//_graph->clear();
		//_graph->compact();
	}

	void processContigs(){

		KminmerParserParallel parser(contigFilename, _minimizerSize, _kminmerSize, false, false, _nbCores);
		parser.parse(CollectBestSupportingReadsFunctor(*this));
	}
	
	
};





#endif
