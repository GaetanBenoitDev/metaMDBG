

#ifndef MDBG_METAG_ASSEMBLYPIPELINE
#define MDBG_METAG_ASSEMBLYPIPELINE

#include "../Commons.hpp"


class AssemblyPipeline : public Tool{
    
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
	bool _useBloomFilter;

	size_t _firstK;
	size_t _lastK;

	string _inputFilenameComplete;
	size_t _meanReadLength;

	AssemblyPipeline(): Tool (){
	}

    void execute (){
		//for(size_t i=0; i<1000; i++)
		execute_pipeline();
	}

	void parseArgs(int argc, char* argv[]){

		args::ArgumentParser parser("", ""); //"This is a test program.", "This goes after the options."
		args::Positional<std::string> arg_outputDir(parser, "outputDir", "Output dir", args::Options::Required);
		args::PositionalList<std::string> arg_readFilenames(parser, "reads", "Input filename(s) (separated by space)", args::Options::Required);
		args::ValueFlag<int> arg_l(parser, "", "Minimizer length", {ARG_MINIMIZER_LENGTH2}, 13);
		args::ValueFlag<float> arg_d(parser, "", "Minimizer density", {ARG_MINIMIZER_DENSITY2}, 0.005f);
		args::ValueFlag<int> arg_nbCores(parser, "", "Number of cores", {ARG_NB_CORES2}, 3);
		args::Flag arg_bf(parser, "", "Filter unique erroneous k-minmers", {ARG_BLOOM_FILTER});
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
		catch (const args::Help&)
		{
			std::cout << parser;
			exit(0);
		}
		catch (const std::exception& e)
		{
			std::cout << parser;
			//cout << endl;
			std::cout << e.what() << endl;
			exit(0);
		}

		if(arg_help){
			std::cout << parser;
			exit(0);
		}

		_inputDir = args::get(arg_outputDir);
		_minimizerSize = args::get(arg_l);
		_minimizerDensity = args::get(arg_d);
		_nbCores = args::get(arg_nbCores);

		_useBloomFilter = false;
		if(arg_bf){
			_useBloomFilter = true;
		}

		createInputFile(args::get(arg_readFilenames));
		//if (arg_l) { 
		//}

		/*
		catch (const args::Help&)
		{
			std::cout << parser;
			exit(0);
		}
		catch (const args::ParseError& e)
		{
			//std::cerr << e.what() << std::endl;
			//std::cerr << parser;
			std::cout << parser;
			std::exit(EXIT_FAILURE);
		}
		*/
		//return 0;

		/*
		cxxopts::Options options("AssemblyPipeline", "");
		options.add_options()
		//("d,debug", "Enable debugging") // a bool parameter
		("reads", "", cxxopts::value<string>())
		//("asmDir", "", cxxopts::value<string>())
		("outputDir", "", cxxopts::value<string>())
		//(ARG_INPUT_FILENAME, "", cxxopts::value<string>())
		//(ARG_OUTPUT_DIR, "", cxxopts::value<string>())
		(ARG_INPUT_FILENAME_TRUTH, "", cxxopts::value<string>()->default_value(""))
		(ARG_MINIMIZER_LENGTH, "", cxxopts::value<int>()->default_value("13"))
		(ARG_MINIMIZER_DENSITY, "", cxxopts::value<float>()->default_value("0.005"))
		(ARG_NB_CORES, "", cxxopts::value<int>()->default_value(NB_CORES_DEFAULT))
		(ARG_BLOOM_FILTER, "", cxxopts::value<bool>()->default_value("false"));
		//(ARG_KMINMER_LENGTH, "", cxxopts::value<int>()->default_value("3"))
		//("k,kminmerSize", "File name", cxxopts::value<std::string>())
		//("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
		//;

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
		options.parse_positional({"reads", "outputDir"});
		options.positional_help("reads outputDir");

		if(argc <= 1){
			std::cout << options.help() << std::endl;
			exit(0);
		}


		cxxopts::ParseResult result;

		try{
			result = options.parse(argc, argv);


			_inputFilename = result["reads"].as<string>();
			//_inputDir = result["asmDir"].as<string>();
			//_outputDir_binning = result["outputDir"].as<string>() + "/bins_/";

			//_kminmerSize = result[ARG_KMINMER_LENGTH].as<int>(); //getInput()->getInt(STR_KMINMER_SIZE);
			//_inputFilename = result[ARG_INPUT_FILENAME].as<string>(); //getInput()->getStr(STR_INPUT);
			_inputDir = result["outputDir"].as<string>(); //getInput()->getStr(STR_OUTPUT);
			_truthInputFilename = result[ARG_INPUT_FILENAME_TRUTH].as<string>();
			_nbCores = result[ARG_NB_CORES].as<int>();
			_minimizerSize = result[ARG_MINIMIZER_LENGTH].as<int>(); //getInput()->getInt(STR_MINIM_SIZE);
			_minimizerDensity = result[ARG_MINIMIZER_DENSITY].as<float>(); //getInput()->getDouble(STR_DENSITY);
			_useBloomFilter = result[ARG_BLOOM_FILTER].as<bool>();

		}
		catch (const std::exception& e){
			std::cout << options.help() << std::endl;
			std::cerr << e.what() << std::endl;
			std::exit(EXIT_FAILURE);
		}
		*/



        fs::path path(_inputDir);
	    if(!fs::exists (path)) fs::create_directories(path); 

		//T exception{text};




		//_inputFilename = result[STR_INPUT].as<string>();

		//_compositionManager = new CompositionManager(4);

		//int w = 80;
		


		//if(!System::file().doesExist (_outputDir)) System::file().mkdir(_outputDir, -1);
		//createInputFile();

		cout << endl;
		cout << "Input: " << _inputFilename << endl;
		cout << "Dir: " << _inputDir << endl;
		cout << "Minimizer length: " << _minimizerSize << endl;
		//cout << "Kminmer length: " << _kminmerSize << endl;
		cout << "Density: " << _minimizerDensity << endl;
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

		string filename = _inputDir + "/" + FILENAME_NO_KMINMER_READS;
		if(fs::exists(filename)){
			fs::remove(filename);
		}
		
	}

    void execute_pipeline(){
		//float density = 0.005;
		//u_int16_t minimizerSize = 21;

		_firstK = 4;

		u_int64_t meanReadLength = computeMeanReadLength(_inputFilename);
		_lastK = meanReadLength*_minimizerDensity*1.2f; //*0.95
		_meanReadLength = meanReadLength;
		//_lastK = 41;
		
		cout << "Mean read length: " << meanReadLength << endl;
		cout << "Min k: " << _firstK << endl;
		cout << "Max k: " << _lastK << endl;
		cout << endl;


		string command = "";

		ofstream fileSmallContigs(_inputDir + "/small_contigs.bin");
		fileSmallContigs.close();

		writeParameters(_minimizerSize, _firstK, _minimizerDensity, _firstK, _firstK, _lastK);
		//createInputFile(false);

		//Read selection
		command = _filename_exe + " readSelection -i " + _inputFilename + " -o " + _inputDir + " -f " + _inputDir + "/read_data_init.txt" + " -t " + to_string(_nbCores);
		Utils::executeCommand(command);
		

		u_int64_t pass = 0;
		u_int32_t prevK = -1;

		for(size_t k=_firstK; k<_lastK; k+=1){

			//cout << "Start asm: " << k << endl;


			//if(pass > 0) createInputFile(true);

			//if(pass <= 0){
			//	prevK = k;
			//	pass += 1;
			//	continue;
			//}

			//string inputFilename = "";
			//if(pass == 0){
			//	inputFilename = _inputFilename;
			//}
			//else{
			//	inputFilename = createInputFile(prevK);
			//}

			//Read selection
			//command = _filename_exe + " readSelection -i " + inputFilename + " -o " + _inputDir + " -f " + _inputDir + "/read_data.gz";
			//if(pass == 0) command += " --firstpass";
			//executeCommand(command);

			executePass(k, prevK, pass);

		
			//./bin/mdbgAsmMeta bin -o ~/workspace/run/overlap_test_multik_201/pass_0 -c ~/workspace/run/overlap_test_multik_201/contigs_10.fasta.gz

			//}
			
			//time ./bin/mdbgAsmMeta countKmer -i ~/workspace/data/AD/shortreads/input_10.txt -c /home/gats/workspace/run/overlap_test_multik_AD/contigs.nodepath.gz.fasta.gz.fasta.gz  -o ~/workspace/run/overlap_test_multik_AD/contig_coverages.tsv
			//!!!Binning: assembly2, permettre de choisir le nom du fichier de contig via -c, actuellement c'est contigs.fasta.gz
			//./bin/mdbgAsmMeta bin -o ~/workspace/run/overlap_test_multik_AD/ -a ~/workspace/run/overlap_test_multik_AD/contig_coverages.tsv

			/*
			//Contig selection
			const string& inputFilenameContig =  _inputDir + "/tmpInputContig.txt";
			ofstream inputFileContig(inputFilenameContig);
			inputFileContig << _inputDir + "/tmpContigs.fasta.gz" << endl;
			inputFileContig.close();
			
			command = _filename_exe + " readSelection -i " + inputFilenameContig + " -f " + _inputDir + "/contig_data.gz" + " -o " + _inputDir;
			executeCommand(command);
			*/



			prevK = k;
			//getchar();
			pass += 1;
			//if(pass == 2) break;
			//if(k >= 21) break;
			//exit(1);
			//break;
			cout << "pass done" << endl;
			//if(generatedContigs) getchar();
			//getchar();
			//if(pass >= 2) break;
			//if(k > 30) getchar();
		}


		executePass(_lastK, prevK, pass);
		cout << "pass done" << endl;

    }

	void executePass(size_t k, size_t prevK, size_t pass){

		bool isFinalPass = k == _lastK;

		writeParameters(_minimizerSize, k, _minimizerDensity, _firstK, prevK, _lastK);

		string command = "";

		if(pass == 0){
			command = _filename_exe + " graph -i " + _inputFilename + " -o " + _inputDir + " -t " + to_string(_nbCores);
			if(_useBloomFilter){
				command += " --bf";
			}
		}
		else{
			command = _filename_exe + " graph -i " + _inputDir + "/read_data.txt.corrected.txt " + " -o " + _inputDir + " -t " + to_string(_nbCores);	
		}
		if(pass == 0) command += " --firstpass";
		Utils::executeCommand(command);
		//getchar();

		//command = _filename_exe + " multik -o " + _inputDir + " -t " + to_string(_nbCores);
		//if(!_truthInputFilename.empty()) command += " --itruth " + _truthInputFilename;
		//if(pass > 0) command += " -c " +  _inputDir + "/contig_data.gz";
		//executeCommand(command);

		//Generate contigs


		//getchar();
		//getchar();


		/*
		command = _filename_exe + " toMinspace " + " -o " + _inputDir + " -c " + _inputDir + "/unitigs.nodepath.gz" + " -f " + _inputDir + "/unitig_data.txt";
		executeCommand(command);
		*/


		//bool generatedContigs = false;
		//if(k == 5 || k == 10 || k == 16 || k == 21 || k == 26 || k == 31){
		if(isFinalPass){

			
			//Generate contigs
			command = _filename_exe + " contig " + " -o " + _inputDir + " --final " + " -t " + to_string(_nbCores);
			if(!_truthInputFilename.empty()) command += " --itruth " + _truthInputFilename;
			//if(pass == 0) command += " --firstpass";
			Utils::executeCommand(command);

			//getchar();

			command = _filename_exe + " toMinspace " + " -o " + _inputDir + " -c " + _inputDir + "/contigs.nodepath" + " -f " + _inputDir + "/contig_data.txt";
			Utils::executeCommand(command);
			
			appendSmallContigs();

			//getchar();
			command = _filename_exe + " toBasespaceFast " + " -o " + _inputDir + " -i " + _inputFilename + " -c " + _inputDir + "/contig_data.txt " + " -f " + _inputDir + "/contigs_" + to_string(k) + ".fasta.gz " + " --fasta"  + " -t " + to_string(_nbCores);
			if(pass == 0) command += " --firstpass";
			Utils::executeCommand(command);

			//getchar();
			
			command = _filename_exe + " toBasespace " + " -o " + _inputDir + " -i " + _inputFilename + " -c " + _inputDir + "/contig_data.txt " + " -f " + _inputDir + "/contigs_" + to_string(k) + ".fasta.gz " + " --fasta"  + " -t " + to_string(_nbCores);
			if(pass == 0) command += " --firstpass";
			Utils::executeCommand(command);

			//getchar();
			command = _filename_exe + " polish " + _inputDir + "/contigs_" + to_string(k) + ".fasta.gz " + _inputFilename + " " + _inputDir + " -t " + to_string(_nbCores) + " --qual ";
			Utils::executeCommand(command);
			//generatedContigs = true;
		}
		else{

			command = _filename_exe + " contig " + " -o " + _inputDir + " -t " + to_string(_nbCores);;
			if(!_truthInputFilename.empty()) command += " --itruth " + _truthInputFilename;
			//if(pass == 0) command += " --firstpass";
			Utils::executeCommand(command);
			//getchar();
			

			command = _filename_exe + " toMinspace " + " -o " + _inputDir + " -c " + _inputDir + "/contigs.nodepath" + " -f " + _inputDir + "/unitig_data.txt";
			Utils::executeCommand(command);

			//getchar();
		}	
			
		savePassData(k);

	}

	ofstream _fileContigsAppend;

	void appendSmallContigs(){
		cout << "Append small contigs" << endl;

		string contigFilename = _inputDir + "/contig_data.txt";

		_fileContigsAppend = ofstream(contigFilename, std::ios_base::app);

		KminmerParser parser(_inputDir + "/small_contigs.bin", _minimizerSize, _kminmerSize, false, false);
		auto fp = std::bind(&AssemblyPipeline::appendSmallContigs_read, this, std::placeholders::_1, std::placeholders::_2);
		parser.parseSequences(fp);

		_fileContigsAppend.close();

	}

	void appendSmallContigs_read(const vector<u_int64_t>& readMinimizers, u_int64_t readIndex){
		
		u_int32_t contigSize = readMinimizers.size();
		_fileContigsAppend.write((const char*)&contigSize, sizeof(contigSize));
		_fileContigsAppend.write((const char*)&readMinimizers[0], contigSize*sizeof(u_int64_t));
	}

	void savePassData(u_int64_t k){

		const string& dir = _inputDir + "/pass_k" + to_string(k);

		if(fs::exists(dir)){
			fs::remove_all(dir);
		}

		fs::create_directories(dir);

		//const auto copyOptions = fs::copy_options::overwrite_existing;
		//fs::copy(_inputDir + "/read_data.txt", dir + "/read_data.txt");


		//if(fs::exists(_inputDir + "/read_index.txt")) fs::copy(_inputDir + "/read_index.txt", dir + "/read_index.txt");
		if(fs::exists(_inputDir + "/minimizer_graph_u.gfa")) fs::copy(_inputDir + "/minimizer_graph_u.gfa", dir + "/minimizer_graph_u.gfa");
		fs::copy(_inputDir + "/minimizer_graph_u_cleaned.gfa", dir + "/minimizer_graph_u_cleaned.gfa");
		//fs::copy(_inputDir + "/minimizer_graph_cleaned.gfa", dir + "/minimizer_graph_cleaned.gfa");
		fs::copy(_inputDir + "/parameters.gz", dir + "/parameters.gz");

		if(k == _firstK || k == _lastK){
			fs::copy(_inputDir + "/minimizer_graph.gfa", dir + "/minimizer_graph.gfa");
			fs::copy(_inputDir + "/kminmerData_min.txt", dir + "/kminmerData_min.txt");
			if(fs::exists(_inputDir + "/nodeName_to_unitigIndex.bin")) fs::copy(_inputDir + "/nodeName_to_unitigIndex.bin", dir + "/nodeName_to_unitigIndex.bin");
			if(fs::exists(_inputDir + "/groundtruth_position.csv")) fs::copy(_inputDir + "/groundtruth_position.csv", dir + "/groundtruth_position.csv");
			if(fs::exists(_inputDir + "/read_path.txt")) fs::copy(_inputDir + "/read_path.txt", dir + "/read_path.txt");
			if(fs::exists(_inputDir + "/read_path_cleaned.txt")) fs::copy(_inputDir + "/read_path_cleaned.txt", dir + "/read_path_cleaned.txt");
		}
		//fs::copy(_inputDir + "/mdbg_nodes_init.gz", dir + "/mdbg_nodes_init.gz", copyOptions);

		//fs::copy(_inputDir + "/unitig_data.txt", dir + "/unitig_data.txt"); //Debug
	}



	void writeParameters(size_t minimizerSize, size_t k, float density, size_t firstK, size_t prevK, size_t lastK){


		float minimizerSpacingMean = 1 / density;
		float kminmerLengthMean = minimizerSpacingMean * (k-1);
		float kminmerOverlapMean = kminmerLengthMean - minimizerSpacingMean;

		string filename_parameters = _inputDir + "/parameters.gz";
		gzFile file_parameters = gzopen(filename_parameters.c_str(),"wb");
		gzwrite(file_parameters, (const char*)&minimizerSize, sizeof(minimizerSize));
		gzwrite(file_parameters, (const char*)&k, sizeof(k));
		gzwrite(file_parameters, (const char*)&density, sizeof(density));
		gzwrite(file_parameters, (const char*)&firstK, sizeof(firstK));
		gzwrite(file_parameters, (const char*)&minimizerSpacingMean, sizeof(minimizerSpacingMean));
		gzwrite(file_parameters, (const char*)&kminmerLengthMean, sizeof(kminmerLengthMean));
		gzwrite(file_parameters, (const char*)&kminmerOverlapMean, sizeof(kminmerOverlapMean));
		gzwrite(file_parameters, (const char*)&prevK, sizeof(prevK));
		gzwrite(file_parameters, (const char*)&lastK, sizeof(lastK));
		gzwrite(file_parameters, (const char*)&_meanReadLength, sizeof(_meanReadLength));
		gzclose(file_parameters);
	}

	void createInputFile(auto& paths){

		//string inputFilenames = _inputFilename;
		_inputFilename = _inputDir + "/input.txt";
		ofstream inputFile(_inputFilename);

		//cout << inputFilenames << endl;
		//vector<string>* fields = new vector<string>();
		//GfaParser::tokenize(inputFilenames, fields, ',');

		for(auto &&path : paths){
			//string filenameAbs = fs::absolute(filename);
			cout << path << endl; //" " << fs::canonical(filename) << " " << fs::weakly_canonical(filename) << endl;
			inputFile << path << endl;
		}

		//delete fields;
		inputFile.close();


	}
	/*
	void createInputFile(bool useContigs){
		_inputFilenameComplete = _inputDir + "/input.txt";
		ofstream inputFile(_inputFilenameComplete);

		ReadParser readParser(_inputFilename, false);

		for(const string& filename : readParser._filenames){
			inputFile << filename << endl;
		}

		if(useContigs){
			inputFile << _inputDir + "/tmpContigs.fasta.gz"  << endl;
		}

		inputFile.close();
	}
	*/

	string createInputFile(u_int32_t k){
		string inputFilename = _inputDir + "/input.txt";
		ofstream inputFile(inputFilename);

		//const string& filename = _inputDir + "/correctedReads_" + to_string(k) + ".min.gz.bitset";
		const string& filename = _inputDir + "/read_data.txt.corrected.txt";
		inputFile << filename << endl;
		/*
		ReadParser readParser(_inputFilename, false);

		for(const string& filename : readParser._filenames){
			inputFile << filename << endl;
		}

		if(useContigs){
			inputFile << _inputDir + "/tmpContigs.fasta.gz"  << endl;
		}
		*/

		inputFile.close();

		return inputFilename;
	}

	u_int64_t _readLengthSum;
	u_int64_t _readLengthN;

	u_int64_t computeMeanReadLength(const string& filename){

		_readLengthSum = 0;
		_readLengthN = 0;

		auto fp = std::bind(&AssemblyPipeline::computeMeanReadLength_read, this, std::placeholders::_1);
		ReadParser readParser(filename, false, false);
		readParser._maxReads = 10000;
		readParser.parse(fp);

		return _readLengthSum / _readLengthN;
	}
	
	void computeMeanReadLength_read(const Read& read){
		_readLengthSum += read._seq.size();
		_readLengthN += 1;
	}

};	


#endif 



