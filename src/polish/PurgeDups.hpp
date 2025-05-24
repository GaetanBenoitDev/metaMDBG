

#ifndef MDBG_METAG_PURGEDUPS
#define MDBG_METAG_PURGEDUPS

#include "../Commons.hpp"
#include "../utils/edlib.h"
#include "../utils/DnaBitset.hpp"
#include "ContigDerep.hpp"



class PurgeDups : public Tool{
    
public:

	//string _inputFilename_reads;
	string _inputFilename_contigs;
	int _nbCores;
	//string _mapperOutputExeFilename;
	string _outputDir;
	string _tmpDir;
	//bool _cut_contigEnds;
	//bool _removeContainedOnly;
	//bool _cut_contigInternal;
	//bool _dontOuputContigs;
	float _minIdentity;
	//u_int64_t _minCircularLength;
	//u_int64_t _minLinearLength;
	
	string _outputFilename_contigs;
	string _outputFilename_mapping;
	//u_int64_t _maxBases;

	//u_int64_t _minDuplicationLength_ends;
	//float _minDuplicationIdentity_ends;
	//u_int64_t _minDuplicationLength_internal;
	//float _minDuplicationIdentity_internal;
	//ofstream _file_duplicatedContigs;
	//u_int64_t _maxMemory;

	struct Alignment{
		string _contigName;
		//u_int64_t _readIndex;
		bool _strand;
		u_int32_t _readStart;
		u_int32_t _readEnd;
		u_int32_t _contigStart;
		u_int32_t _contigEnd;
		//float _score;
		//u_int64_t _length;
		
		u_int64_t length() const{
			return std::max((u_int64_t)(_readEnd - _readStart), (u_int64_t)(_contigEnd - _contigStart));
		}
		
	};
	/*
	struct Alignment{
		u_int32_t _contigIndex;
		//u_int64_t _readIndex;
		bool _strand;
		u_int32_t _readStart;
		u_int32_t _readEnd;
		u_int32_t _contigStart;
		u_int32_t _contigEnd;
		//float _score;
		//u_int64_t _length;
		
		u_int64_t length() const{
			return std::max((u_int64_t)(_readEnd - _readStart), (u_int64_t)(_contigEnd - _contigStart));
		}
		
	};
	*/

	PurgeDups(): Tool (){

	}


	void parseArgs(int argc, char* argv[]){




		string ARG_CUT_INTERNAL = "cutinternal";
		string ARG_CONTAINED_ONLY = "containedOnly";
		string ARG_NO_DUMP= "nodump";

		
		args::ArgumentParser parser("derep", ""); //"This is a test program.", "This goes after the options."
		args::Positional<std::string> arg_contigs(parser, "contigs", "Input contig filename", args::Options::Required);
		args::Positional<std::string> arg_outputFilename(parser, "outputFilename", "Output contig filename", args::Options::Required);
		args::Positional<std::string> arg_outputDir(parser, "outputDir", "Output dir for temporary files", args::Options::Required);
		//args::Positional<std::string> arg_contigs(parser, "contigs", "", args::Options::Required);
		//args::PositionalList<std::string> arg_readFilenames(parser, "reads", "Input filename(s) (separated by space)", args::Options::Required);
		//args::ValueFlag<int> arg_l(parser, "", "Minimizer length", {ARG_MINIMIZER_LENGTH2}, 13);
		//args::ValueFlag<float> arg_d(parser, "", "Minimizer density", {ARG_MINIMIZER_DENSITY2}, 0.005f);
		//args::ValueFlag<std::string> arg_contigs(parser, "", "", {ARG_INPUT_FILENAME_CONTIG}, "");
		//args::ValueFlag<int> arg_length(parser, "", "Use contigs with length > l as reference for strain duplication purging", {ARG_LINEAR_LENGTH}, 1000000);
		args::ValueFlag<float> arg_minIdentity(parser, "", "Minimum identity for strain purging (0-1)", {ARG_MIN_IDENTITY}, 0.99);
		args::ValueFlag<int> arg_nbCores(parser, "", "Number of cores", {ARG_NB_CORES2}, NB_CORES_DEFAULT_INT);
		//args::Flag arg_containedOnly(parser, "", "", {ARG_CONTAINED_ONLY});
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

		_outputDir = args::get(arg_outputDir);
		_inputFilename_contigs = args::get(arg_contigs);
		_outputFilename_contigs = args::get(arg_outputFilename);
		_nbCores = args::get(arg_nbCores);

		//_minLinearLength = args::get(arg_length);
		_minIdentity = args::get(arg_minIdentity);
		//_minIdentity *= 100;

		//_dontOuputContigs = false;
		//if(arg_noDump){
		//	_dontOuputContigs = true;
		//}

		//_cut_contigInternal = false;
		//if(arg_cutInternal){
		//	_cut_contigInternal = true;
		//}

		//_removeContainedOnly = false;
		//if(arg_containedOnly){
		//	_removeContainedOnly = true;
		//}

		if (_outputFilename_contigs.find(".gz") == std::string::npos) {
			_outputFilename_contigs += ".gz";
		}


		if(_inputFilename_contigs == _outputFilename_contigs){
			cerr << "Output filename == input filename" << endl;
			exit(0);
		}

		/*
		//string ARG_CUT_ENDS = "cutends";
		string ARG_CUT_INTERNAL = "cutinternal";
		string ARG_NO_DUMP= "nodump";

		//string filenameExe = argv[0];
		//_logFile << filenameExe << endl;

		//fs::path pa(filenameExe);
		//_mapperOutputExeFilename = pa.parent_path().string() + "/mapper";
		//_logFile << _mapperOutputExeFilename << endl;
		//exit(1);

		cxxopts::Options options("ToBasespace", "");
		options.add_options()
		("contigs", "", cxxopts::value<string>())
		("outputFilename", "", cxxopts::value<string>())
		//("reads", "", cxxopts::value<string>())
		("tmpDir", "", cxxopts::value<string>())
		//(ARG_CUT_ENDS, "", cxxopts::value<bool>()->default_value("false"))
		(ARG_CUT_INTERNAL, "", cxxopts::value<bool>()->default_value("false"))
		(ARG_NO_DUMP, "", cxxopts::value<bool>()->default_value("false"))
		(ARG_NB_CORES, "", cxxopts::value<int>()->default_value(NB_CORES_DEFAULT));


		options.parse_positional({"contigs", "outputFilename", "tmpDir"});
		options.positional_help("contigs outputFilename tmpDir");


		//("k,kminmerSize", "File name", cxxopts::value<std::string>())
		//("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
		//;

		if(argc <= 1){
			_logFile << options.help() << endl;
			exit(0);
		}

		cxxopts::ParseResult result;

		try{
			result = options.parse(argc, argv);

			//_inputFilename_reads = result["reads"].as<string>();
			_inputFilename_contigs = result["contigs"].as<string>();
			_outputFilename_contigs = result["outputFilename"].as<string>();
			_outputDir = result["tmpDir"].as<string>();
			//_cut_contigEnds = result[ARG_CUT_ENDS].as<bool>();
			_cut_contigInternal = result[ARG_CUT_INTERNAL].as<bool>();
			_dontOuputContigs = result[ARG_NO_DUMP].as<bool>();
			_nbCores = result[ARG_NB_CORES].as<int>();
		}
		catch (const std::exception& e){
			std::_logFile << options.help() << std::endl;
			std::cerr << e.what() << std::endl;
			std::exit(EXIT_FAILURE);
		}
		*/




		//fs::path p(_inputFilename_contigs);
		//while(p.has_extension()){
		//	p.replace_extension("");
		//}

		_tmpDir = _outputDir;// + "/tmp/";
		if(!fs::exists(_tmpDir)){
			fs::create_directories(_tmpDir);
		}

		fs::path outputContigPath(_outputFilename_contigs);
		string contigDir = outputContigPath.parent_path();
		if(!fs::exists(contigDir)){
			fs::create_directories(contigDir);
		}

		//_outputFilename_contigs = p.string() + "_derep.fasta.gz";
		//_maxBases = 200000000ull;
		//_minDuplicationLength_ends = 1000;
		//_minDuplicationIdentity_ends = 0.95;
		//_minDuplicationLength_internal = 1000;
		//_minDuplicationIdentity_internal = 0.95;
		//_outputFilename_contigs = p.string() + "_corrected.fasta.gz";
		_outputFilename_mapping = _tmpDir + "/_tmp_mapping_derep__.paf.gz";
		//_maxMemory = 4000000000ull;

		//_file_duplicatedContigs = ofstream(_outputDir + "/duplicatedContigs.txt");
		
		
		openLogFile(_outputDir);

		Logger::get().debug() << "";
		Logger::get().debug() << "Contigs: " << _inputFilename_contigs;
		//_logFile << "Cut contigs ends: " << _cut_contigEnds << endl;
		//_logFile << "Cut contigs internal: " << _cut_contigInternal << endl;
		Logger::get().debug() << "Output filename: " << _outputFilename_contigs;
		Logger::get().debug() << "";

	}

	gzFile _outputContigFile;
	unordered_map<string, string> _contigSequences;
	unordered_map<string, vector<Alignment>> _alignments;

	/*
	//unordered_map<string, string> _debug_contigIndex_to_contigName;

	void debug_indexContigName(){
		
		auto fp = std::bind(&PurgeDups::indexContigName_read, this, std::placeholders::_1);
		ReadParser readParser(_inputFilename_contigs, true, false);
		readParser.parse(fp);
	}

	void indexContigName_read(const Read& read){
		_debug_contigIndex_to_contigName[read._index] = read._header;
	}
	*/

    void execute (){

		//const string derepTmp = _tmpDir + "/derep.fasta.gz";
		//{
			ContigDerep contigDerep(_inputFilename_contigs, _outputFilename_contigs, "4", _tmpDir, _tmpDir, 0.9, _nbCores);
			contigDerep.execute();
		//}

		//{
		//	ContigDerep contigDerep(derepTmp, _outputFilename_contigs, "4", _tmpDir, _tmpDir, 0.9, _nbCores);
		//	contigDerep.execute();
		//}

	}


};	

#endif 


