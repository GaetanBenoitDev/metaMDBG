

#ifndef MDBG_METAG_BINNINGPIPELINE
#define MDBG_METAG_BINNINGPIPELINE

#include "../Commons.hpp"


class BinningPipeline : public Tool{
    
public:

	string _inputFilename;
	string _inputDir;
	float _minimizerDensity;
    size_t _minimizerSize;
    size_t _kminmerSize;
	string _filename_readMinimizers;
	string _filename_exe;
	string _truthInputFilename;
	int _nbCores;
	string _filename_inputContigs;
	bool _computeBinStats;
	string _outputDir_binning;

	string _inputFilenameComplete;
	string _filename_abundance;

	BinningPipeline(): Tool (){
	}

    void execute (){
		execute_pipeline();
	}

	void parseArgs(int argc, char* argv[]){

		cxxopts::Options options("AssemblyPipeline", "");
		options.add_options()
		//("d,debug", "Enable debugging") // a bool parameter
		//(ARG_INPUT_FILENAME, "", cxxopts::value<string>())
		("contig", "", cxxopts::value<string>())
		("asmDir", "", cxxopts::value<string>())
		("outputDir", "", cxxopts::value<string>())
		(ARG_INPUT_FILENAME_TRUTH, "", cxxopts::value<string>()->default_value(""))
		(ARG_EVAL, "", cxxopts::value<bool>()->default_value("false"))
		(ARG_INPUT_FILENAME_ABUNDANCE, "", cxxopts::value<string>()->default_value(""))
		(ARG_NB_CORES, "", cxxopts::value<int>()->default_value("8"));
		//(ARG_KMINMER_LENGTH, "", cxxopts::value<int>()->default_value("3"))
		//(ARG_MINIMIZER_LENGTH, "", cxxopts::value<int>()->default_value("21"))
		//(ARG_MINIMIZER_DENSITY, "", cxxopts::value<float>()->default_value("0.005"));
		//("k,kminmerSize", "File name", cxxopts::value<std::string>())
		//("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
		//;
		options.parse_positional({"contig", "asmDir", "outputDir"});
		options.positional_help("contigs asmDir outputDir");
		if(argc <= 1){
			std::cout << options.help() << std::endl;
			exit(0);
		}


		cxxopts::ParseResult result;

		try{
			result = options.parse(argc, argv);

			//_kminmerSize = result[ARG_KMINMER_LENGTH].as<int>(); //getInput()->getInt(STR_KMINMER_SIZE);
			//_inputFilename = result[ARG_INPUT_FILENAME].as<string>(); //getInput()->getStr(STR_INPUT);
			_filename_inputContigs = result["contig"].as<string>();
			_inputDir = result["asmDir"].as<string>();
			_outputDir_binning = result["outputDir"].as<string>() + "/bins_/";
			//_minimizerSize = result[ARG_MINIMIZER_LENGTH].as<int>(); //getInput()->getInt(STR_MINIM_SIZE);
			//_inputDir = result[ARG_OUTPUT_DIR].as<string>(); //getInput()->getStr(STR_OUTPUT);
			_truthInputFilename = result[ARG_INPUT_FILENAME_TRUTH].as<string>();
			_filename_abundance = result[ARG_INPUT_FILENAME_ABUNDANCE].as<string>();
			_nbCores = result[ARG_NB_CORES].as<int>();
			_computeBinStats = result[ARG_EVAL].as<bool>();
			//_minimizerDensity = result[ARG_MINIMIZER_DENSITY].as<float>(); //getInput()->getDouble(STR_DENSITY);

		}
		catch (const std::exception& e){
			std::cout << options.help() << std::endl;
			std::cerr << e.what() << std::endl;
			std::exit(EXIT_FAILURE);
		}

        //fs::path path(_inputDir);
	    //if(!fs::exists (path)) fs::create_directory(path); 


	    if(fs::exists (_outputDir_binning)){
			fs::remove_all(_outputDir_binning);
		}
		fs::create_directories(_outputDir_binning); 
		//T exception{text};




		//_inputFilename = result[STR_INPUT].as<string>();

		//_compositionManager = new CompositionManager(4);

		//int w = 80;
		


		//if(!System::file().doesExist (_outputDir)) System::file().mkdir(_outputDir, -1);
		cout << _filename_inputContigs << " " << _inputDir << " " << _outputDir_binning << endl;

		cout << endl;
		cout << "Input: " << _filename_inputContigs << endl;
		cout << "Dir: " << _inputDir << endl;
		//cout << "Minimizer length: " << _minimizerSize << endl;
		//cout << "Kminmer length: " << _kminmerSize << endl;
		//cout << "Density: " << _minimizerDensity << endl;
		cout << endl;

		_filename_exe = argv[0];
		/*
		//_inputDir = getInput()->get(STR_INPUT_DIR) ? getInput()->getStr(STR_INPUT_DIR) : "";
		//_input_extractKminmers= getInput()->get(STR_INPUT_EXTRACT_KMINMERS) ? getInput()->getStr(STR_INPUT_EXTRACT_KMINMERS) : "";

		_filename_readMinimizers = _outputDir + "/read_data.gz";
		//_filename_readCompositions = _outputDir + "/read_compositions.gz";
		//_filename_filteredMinimizers = _outputDir + "/filteredMinimizers.gz";
		//_filename_hifiasmGroundtruth = _outputDir + "/hifiasmGroundtruth.gz";


		string filename_parameters = _outputDir + "/parameters.gz";
		gzFile file_parameters = gzopen(filename_parameters.c_str(),"wb");
		gzwrite(file_parameters, (const char*)&_minimizerSize, sizeof(_minimizerSize));
		gzwrite(file_parameters, (const char*)&_kminmerSize, sizeof(_kminmerSize));
		gzwrite(file_parameters, (const char*)&_minimizerDensity, sizeof(_minimizerDensity));
		gzclose(file_parameters);
		*/

		
	}

    void execute_pipeline(){
		
		long firstK = 4;
		string command = "";

		u_int64_t pass = 0;
		string binningFilename;
		long prevK = -1;

		string lastBinFilename = "";

		for(long k=firstK; k>=4; k-=10){

			string ouputDir = _inputDir + "/pass_k" + to_string(k);
			const string& binningFilename_input = _inputDir + "/pass_k" + to_string(prevK) + "/contigToBin.bin";
			const string& binningFilename_output = _inputDir + "/pass_k" + to_string(k) + "/contigToBin.bin";
			lastBinFilename = binningFilename_output;

			command = _filename_exe + " binPass -o " + ouputDir + " -c " + _filename_inputContigs + " --bo " + binningFilename_output;
			if(pass > 0) command += " --bi " + binningFilename_input;
			if(pass == 0) command += " --firstpass ";
			if(!_filename_abundance.empty()) command += " -a " + _filename_abundance;
			if(_computeBinStats || k == 4) command += " --eval ";
			executeCommand(command);
			
			prevK = k;
			pass += 1;
		}

		dumpBins(_filename_inputContigs, lastBinFilename);
    }

	/*
	void savePassData(u_int64_t k){

		const string& dir = _inputDir + "/pass_k" + to_string(k);

		if(fs::exists(dir)){
			fs::remove_all(dir);
		}

		fs::create_directory(dir);

		//const auto copyOptions = fs::copy_options::overwrite_existing;
		fs::copy(_inputDir + "/read_data.txt", dir + "/read_data.txt");
		fs::copy(_inputDir + "/minimizer_graph.gfa", dir + "/minimizer_graph.gfa");
		fs::copy(_inputDir + "/parameters.gz", dir + "/parameters.gz");
		fs::copy(_inputDir + "/mdbg_nodes.gz", dir + "/mdbg_nodes.gz");
		//fs::copy(_inputDir + "/mdbg_nodes_init.gz", dir + "/mdbg_nodes_init.gz", copyOptions);
	}
	*/

	ContigFeature _contigFeature; 

	void dumpBins(const string& contigFilename, const string binFilename){

		_contigFeature.loadContigBins(binFilename);
		loadContigs(contigFilename);

		u_int64_t binIndex = 0;
		
		

		for(const auto& it : _contigFeature._binIndex_to_contigIndex){
			bool isValid = dumpBin(binIndex, it.second);
			if(isValid) binIndex += 1;
		}

		cout << "Nb bins: " << binIndex << endl;
		cout << "Result dir: " << _outputDir_binning << endl;
	}
	
	
	void loadContigs(const string& filename){
		auto fp = std::bind(&BinningPipeline::loadContigs_read, this, std::placeholders::_1);
		ReadParser readParser(filename, true, false);
		readParser.parse(fp);
	}

	void loadContigs_read(const Read& read){
		if(_contigFeature._contigIndex_to_binIndex.find(read._index) == _contigFeature._contigIndex_to_binIndex.end()) return;
		_contigFeature._contigSequences[read._index] = new DnaBitset(read._seq);
	}


	struct ContigSequence{
		u_int64_t _contigIndex;
		string _sequence;
	};

	bool dumpBin(const u_int64_t binIndex, const vector<u_int32_t>& binContigIndexes){

		vector<ContigSequence> binSequences;
		for(u_int32_t contigIndex : binContigIndexes){

			char* seq = _contigFeature._contigSequences[contigIndex]->to_string();
			binSequences.push_back({contigIndex, string(seq)});
			free(seq);
		}


		u_int64_t lengthTotal = 0;
		for(const ContigSequence& contig : binSequences){
			lengthTotal += contig._sequence.size();
		}

		if(lengthTotal < 20000) return false;

		const string& filename = _outputDir_binning + "/bin_" + to_string(binIndex) + ".fasta";
		ofstream file = ofstream(filename);

		for(const ContigSequence& contig : binSequences){
			string header = ">ctg" + to_string(contig._contigIndex);
			file << header << endl;
			file << contig._sequence << endl;
			
		}

		file.close();

		return true;
	}

	void executeCommand(const string& command){
		cout << command << endl;
		int ret = system(command.c_str());
		if(ret != 0){
			cerr << "Command failed: " << ret << endl;
			exit(ret);
		}
	}

	void writeParameters(size_t minimizerSize, size_t k, float density, size_t firstK){

		string filename_parameters = _inputDir + "/parameters.gz";
		gzFile file_parameters = gzopen(filename_parameters.c_str(),"wb");
		gzwrite(file_parameters, (const char*)&minimizerSize, sizeof(minimizerSize));
		gzwrite(file_parameters, (const char*)&k, sizeof(k));
		gzwrite(file_parameters, (const char*)&density, sizeof(density));
		gzwrite(file_parameters, (const char*)&firstK, sizeof(firstK));
		gzclose(file_parameters);
	}

};	


#endif 


