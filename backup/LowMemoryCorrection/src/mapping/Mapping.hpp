

#ifndef MDBG_METAG_MAPPING
#define MDBG_METAG_MAPPING

#include "../Commons.hpp"

class Mapping : public Tool{
    
public:

	string _contigFilename;
	string _binDir;
	string _outputFilename;

	MinimizerParser* _minimizerParser;
	size_t _minimizerSize;
	size_t _kminmerSize;
	float _minimizerDensity;

	Mapping(): Tool (){

	}


	void parseArgs(int argc, char* argv[]){

		_kminmerSize = 4;
		_minimizerSize = 21;
		_minimizerDensity = 0.05;

		/*
		cxxopts::Options options("ToBasespace", "");
		options.add_options()
		//(ARG_INPUT_FILENAME, "", cxxopts::value<string>())
		//(ARG_INPUT_FILENAME_CONTIG, "", cxxopts::value<string>())
		//(ARG_INPUT_FILENAME_TRUTH, "", cxxopts::value<string>()->default_value(""))
		//(ARG_OUTPUT_FILENAME, "", cxxopts::value<string>())
		//(ARG_FIRST_PASS, "", cxxopts::value<bool>()->default_value("false"))
		//(ARG_FASTA, "", cxxopts::value<bool>()->default_value("false"))
		//(ARG_OUTPUT_DIR, "", cxxopts::value<string>())
		//(ARG_NB_CORES, "", cxxopts::value<int>()->default_value(NB_CORES_DEFAULT));
		("contig", "", cxxopts::value<string>())
		("binDir", "", cxxopts::value<string>())
		("outputFilename", "", cxxopts::value<string>());
		//("outputDir", "", cxxopts::value<string>())
		//(ARG_INPUT_FILENAME_TRUTH, "", cxxopts::value<string>()->default_value(""))
		//(ARG_EVAL, "", cxxopts::value<bool>()->default_value("false"))
		//(ARG_INPUT_FILENAME_ABUNDANCE, "", cxxopts::value<string>()->default_value(""))
		//(ARG_NB_CORES, "", cxxopts::value<int>()->default_value("8"));
		//(ARG_KMINMER_LENGTH, "", cxxopts::value<int>()->default_value("3"))
		//(ARG_MINIMIZER_LENGTH, "", cxxopts::value<int>()->default_value("21"))
		//(ARG_MINIMIZER_DENSITY, "", cxxopts::value<float>()->default_value("0.005"));
		//("k,kminmerSize", "File name", cxxopts::value<std::string>())
		//("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
		//;
		options.parse_positional({"contig", "binDir", "outputFilename"});
		options.positional_help("contigs binDir outputFilename");


		//("k,kminmerSize", "File name", cxxopts::value<std::string>())
		//("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
		//;

		if(argc <= 1){
			cout << options.help() << endl;
			exit(0);
		}

		cxxopts::ParseResult result;

		try{
			result = options.parse(argc, argv);

			_contigFilename = result["contig"].as<string>();
			_binDir = result["binDir"].as<string>();
			_outputFilename = result["outputFilename"].as<string>();
			
		}
		catch (const std::exception& e){
			std::cout << options.help() << std::endl;
			std::cerr << e.what() << std::endl;
			std::exit(EXIT_FAILURE);
		}
		*/
		/*
		_inputFilename = getInput()->getStr(STR_INPUT);
		_inputDir = getInput()->getStr(STR_INPUT_DIR);
		*/
		/*
		string filename_parameters = _inputDir + "/parameters.gz";
		gzFile file_parameters = gzopen(filename_parameters.c_str(),"rb");
		gzread(file_parameters, (char*)&_minimizerSize, sizeof(_minimizerSize));
		gzread(file_parameters, (char*)&_kminmerSize, sizeof(_kminmerSize));
		gzread(file_parameters, (char*)&_minimizerDensity, sizeof(_minimizerDensity));
		gzread(file_parameters, (char*)&_kminmerSizeFirst, sizeof(_kminmerSizeFirst));
		gzclose(file_parameters);

		_kminmerSize = _kminmerSizeFirst;

		cout << endl;
		cout << "Input dir: " << _inputDir << endl;
		//cout << "Output filename: " << _outputFilename << endl;
		cout << "Minimizer length: " << _minimizerSize << endl;
		cout << "Kminmer length: " << _kminmerSize << endl;
		cout << "Density: " << _minimizerDensity << endl;
		cout << endl;

		//_filename_outputContigs = _inputDir + "/contigs.fasta.gz";
		//_filename_outputContigs = _inputFilenameContig + ".fasta.gz"; //_inputDir + "/tmpContigs.fasta.gz";
		*/

		args::ArgumentParser parser("mapBin", ""); //"This is a test program.", "This goes after the options."
		args::Positional<std::string> arg_mdbgDir(parser, "contig", "", args::Options::Required);
		args::Positional<std::string> arg_binDir(parser, "binDir", "", args::Options::Required);
		args::Positional<std::string> arg_outputFilename(parser, "outputFilename", "", args::Options::Required);
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

		_contigFilename = args::get(arg_mdbgDir);
		_binDir = args::get(arg_binDir);
		_outputFilename = args::get(arg_outputFilename);


		_minimizerParser = new MinimizerParser(_minimizerSize, _minimizerDensity);
		
		//openLogFile(_tmpDir);
	}


    void execute (){

		_outputFile = ofstream(_outputFilename);
		_outputFile << "Name,Color" << endl;

		extract_truth_kminmers();
		map();

		//closeLogFile();
	}



	unordered_map<KmerVec, vector<string>> _kmervec_to_unitigName;

	//MinimizerParser* _minimizerParser;
	//u_int32_t _extract_truth_kminmers_read_position;
	//unordered_map<u_int32_t, vector<string>> _evaluation_hifiasmGroundTruth_nodeName_to_unitigName;
	//vector<u_int32_t> _evaluation_hifiasmGroundTruth_path;
	//unordered_set<u_int32_t> _hifiasm_startingNodenames;
	//ofstream _file_groundTruth_hifiasm_position;
	EncoderRLE _encoderRLE;

	void extract_truth_kminmers(){

		cout << "Extract contig kminmers" << endl;
		//_file_groundTruth_hifiasm_position.open(_inputDir + "/groundtruth_hifiasm_position.csv");
		//_file_groundTruth_hifiasm_position << "Name,Position" << endl;

		//_extract_truth_kminmers_read_position = 0;
		_minimizerParser = new MinimizerParser(_minimizerSize, _minimizerDensity);
		
		auto fp = std::bind(&Mapping::extract_truth_kminmers_read, this, std::placeholders::_1);
		ReadParser readParser(_contigFilename, true, false, _logFile);
		readParser.parse(fp);

		//_file_groundTruth_hifiasm_position.close();

		//delete _minimizerParser;
	}

	void extract_truth_kminmers_read(const Read& read){
		//ottalSize += strlen(read->seq.s);


		_outputFile << read._header << ",0" << endl;

		u_int64_t readIndex = read._index;

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
			_kmervec_to_unitigName[vec].push_back(read._header);
			/*
			if(_mdbg->_dbg_nodes.find(vec) != _mdbg->_dbg_nodes.end()){

				u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
				

				_evaluation_hifiasmGroundTruth_nodeName_to_unitigName[_mdbg->_dbg_nodes[vec]._index].push_back(read._header);
				_evaluation_hifiasmGroundTruth_path.push_back(_mdbg->_dbg_nodes[vec]._index);


			}
			else{
				_evaluation_hifiasmGroundTruth_path.push_back(-1);
			}
			*/


		}


	}

	ofstream _outputFile;
	//u_int32_t _currentBinIndex;
	string _currentBinName;
	unordered_set<string> _writtenUnitigNames;
	
	void map(){

		cout << "Mapping bins" << endl;



		for (const auto & p : fs::directory_iterator(_binDir)){
			string ext = p.path().extension();
			cout << p.path() << endl;


			cout << ext << endl;
			if(ext == ".fa" || ext == ".fasta" || ext == ".fna"){
				string filename = p.path();
				cout << filename << endl;

				string binName = p.path().stem();
				//binName.erase(binName.find("bin."), 4);
				//binName.erase(binName.find(ext), ext.size());
				cout << binName << endl;
				//_currentBinIndex = stoull(binName);
				_currentBinName = binName;


				mapBin(filename);
			}

		}

		_outputFile.close();
	}


	void mapBin(const string& binFilename){

		auto fp = std::bind(&Mapping::mapBin_read, this, std::placeholders::_1);
		ReadParser readParser(binFilename, true, false, _logFile);
		readParser.parse(fp);
	}

	void mapBin_read(const Read& read){
		//ottalSize += strlen(read->seq.s);


		u_int64_t readIndex = read._index;

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

			if(_kmervec_to_unitigName.find(vec) == _kmervec_to_unitigName.end()) continue;

			const vector<string>& unitigNames = _kmervec_to_unitigName[vec];

			for(const string& unitigName : unitigNames){
				if(_writtenUnitigNames.find(unitigName) != _writtenUnitigNames.end()) continue;

				_outputFile << unitigName << "," << _currentBinName << endl;
				_writtenUnitigNames.insert(unitigName);
			}
			//_kmervec_to_unitigName[vec] = read._header;
			/*
			if(_mdbg->_dbg_nodes.find(vec) != _mdbg->_dbg_nodes.end()){

				u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
				

				_evaluation_hifiasmGroundTruth_nodeName_to_unitigName[_mdbg->_dbg_nodes[vec]._index].push_back(read._header);
				_evaluation_hifiasmGroundTruth_path.push_back(_mdbg->_dbg_nodes[vec]._index);


			}
			else{
				_evaluation_hifiasmGroundTruth_path.push_back(-1);
			}
			*/


		}


	}

};	


#endif 


