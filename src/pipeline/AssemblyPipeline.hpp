

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

	string _inputFilenameComplete;

	AssemblyPipeline(): Tool (){
	}

    void execute (){
		//for(size_t i=0; i<1000; i++)
		execute_pipeline();
	}

	void parseArgs(int argc, char* argv[]){

		cxxopts::Options options("AssemblyPipeline", "");
		options.add_options()
		//("d,debug", "Enable debugging") // a bool parameter
		("reads", "", cxxopts::value<string>())
		//("asmDir", "", cxxopts::value<string>())
		("outputDir", "", cxxopts::value<string>())
		//(ARG_INPUT_FILENAME, "", cxxopts::value<string>())
		//(ARG_OUTPUT_DIR, "", cxxopts::value<string>())
		(ARG_INPUT_FILENAME_TRUTH, "", cxxopts::value<string>()->default_value(""))
		(ARG_MINIMIZER_LENGTH, "", cxxopts::value<int>()->default_value("21"))
		(ARG_MINIMIZER_DENSITY, "", cxxopts::value<float>()->default_value("0.005"))
		(ARG_NB_CORES, "", cxxopts::value<int>()->default_value("8"));
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

		}
		catch (const std::exception& e){
			std::cout << options.help() << std::endl;
			std::cerr << e.what() << std::endl;
			std::exit(EXIT_FAILURE);
		}

        fs::path path(_inputDir);
	    if(!fs::exists (path)) fs::create_directories(path); 

		//T exception{text};




		//_inputFilename = result[STR_INPUT].as<string>();

		//_compositionManager = new CompositionManager(4);

		//int w = 80;
		


		//if(!System::file().doesExist (_outputDir)) System::file().mkdir(_outputDir, -1);
		

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
		size_t firstK = 4;

		string command = "";

		ofstream fileSmallContigs(_inputDir + "/small_contigs.bin");

		writeParameters(_minimizerSize, firstK, _minimizerDensity, firstK);
		//createInputFile(false);

		//Read selection
		command = _filename_exe + " readSelection -i " + _inputFilename + " -o " + _inputDir + " -f " + _inputDir + "/read_data_init.txt" + " -t " + to_string(_nbCores);
		executeCommand(command);
		

		u_int64_t pass = 0;
		u_int32_t prevK = -1;

		for(size_t k=firstK; k<122; k+=1){

			//cout << "Start asm: " << k << endl;


			writeParameters(_minimizerSize, k, _minimizerDensity, firstK);
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

			if(pass == 0){
				command = _filename_exe + " graph -i " + _inputFilename + " -o " + _inputDir + " -t " + to_string(_nbCores);
			}
			else{
				command = _filename_exe + " graph -i " + _inputDir + "/read_data.txt.corrected.txt " + " -o " + _inputDir + " -t " + to_string(_nbCores);	
			}
			if(pass == 0) command += " --firstpass";
			executeCommand(command);
			//getchar();

			command = _filename_exe + " multik -o " + _inputDir + " -t " + to_string(_nbCores);
			if(!_truthInputFilename.empty()) command += " --itruth " + _truthInputFilename;
			//if(pass > 0) command += " -c " +  _inputDir + "/contig_data.gz";
			executeCommand(command);

			command = _filename_exe + " toMinspace " + " -o " + _inputDir + " -c " + _inputDir + "/unitigs.nodepath.gz" + " -f " + _inputDir + "/unitig_data.txt";
			executeCommand(command);

			//command = ./bin/mdbgAsmMeta toMinspace -o ~/workspace/run/overlap_test_multik_AD/ -c ~/workspace/run/overlap_test_multik_AD/contigs.nodepath.gz

			//command = _filename_exe + " toMinspace " + " -o " + _inputDir + " -c " + _inputDir + "/correctedReads_" + to_string(k) + ".min.gz" + " -f "  + _inputDir + "/read_data.txt.corrected.txt";
			//executeCommand(command);



			//command = _filename_exe + " toBasespace -i " + _inputFilename +  " -o " + _inputDir;
			//if(pass > 0) command += " -c " +  _inputDir + "/contigs.min.gz";
			//executeCommand(command);

			//command = _filename_exe + " toBasespaceFast " +  " -o " + _inputDir + " -i " + _inputFilenameComplete;
			//command = _filename_exe + " toBasespaceFast " + " -o " + _inputDir + " -i " + inputFilename + " -c " + _inputDir + "/correctedReads_" + to_string(k) + ".min.gz";
			//if(pass == 0) command += " --firstpass";
			//executeCommand(command);

			//command = _filename_exe + " toBasespaceFast " + " -o " + _inputDir + " -i " + inputFilename + " -c " + _inputDir + "/correctedReads_" + to_string(k) + ".min.gz";
			//if(pass == 0) command += " --firstpass";
			//executeCommand(command);

			
			//if(pass > 0){

			




			//command = _filename_exe + " toBasespaceFast " + " -o " + _inputDir + " -i " + _inputFilename + " -c " + _inputDir + "/contig_data.txt " + " -f " + _inputDir + "/contigs.fasta.gz " +  " --fasta";
			//if(pass == 0) command += " --firstpass";
			//executeCommand(command);	


			bool generatedContigs = false;
			//if(k == 5 || k == 10 || k == 16 || k == 21 || k == 26 || k == 31){
			if(k == 41 || k == 51 || k == 61 || k == 71 || k == 81 || k == 91 || k == 101 || k == 111 || k == 121){


				//Generate contigs
				command = _filename_exe + " contig " + " -o " + _inputDir;
				if(!_truthInputFilename.empty()) command += " --itruth " + _truthInputFilename;
				//if(pass == 0) command += " --firstpass";
				executeCommand(command);

				command = _filename_exe + " toMinspace " + " -o " + _inputDir + " -c " + _inputDir + "/contigs.nodepath.gz" + " -f " + _inputDir + "/contig_data.txt";
				executeCommand(command);

				command = _filename_exe + " toBasespace " + " -o " + _inputDir + " -i " + _inputFilename + " -c " + _inputDir + "/contig_data.txt " + " -f " + _inputDir + "/contigs_" + to_string(k) + ".fasta.gz " + " --fasta"  + " -t " + to_string(_nbCores);
				if(pass == 0) command += " --firstpass";
				executeCommand(command);

				generatedContigs = true;
			}	
		
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

			savePassData(k);

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

    }

	void savePassData(u_int64_t k){

		const string& dir = _inputDir + "/pass_k" + to_string(k);

		if(fs::exists(dir)){
			fs::remove_all(dir);
		}

		fs::create_directories(dir);

		//const auto copyOptions = fs::copy_options::overwrite_existing;
		//fs::copy(_inputDir + "/read_data.txt", dir + "/read_data.txt");

		
		if(fs::exists(_inputDir + "/groundtruth_position.csv")) fs::copy(_inputDir + "/groundtruth_position.csv", dir + "/groundtruth_position.csv");
		if(fs::exists(_inputDir + "/read_path.txt")) fs::copy(_inputDir + "/read_path.txt", dir + "/read_path.txt");
		if(fs::exists(_inputDir + "/read_path_cleaned.txt")) fs::copy(_inputDir + "/read_path_cleaned.txt", dir + "/read_path_cleaned.txt");
		//if(fs::exists(_inputDir + "/read_index.txt")) fs::copy(_inputDir + "/read_index.txt", dir + "/read_index.txt");
		fs::copy(_inputDir + "/minimizer_graph.gfa", dir + "/minimizer_graph.gfa");
		fs::copy(_inputDir + "/minimizer_graph_u.gfa", dir + "/minimizer_graph_u.gfa");
		fs::copy(_inputDir + "/minimizer_graph_u_cleaned.gfa", dir + "/minimizer_graph_u_cleaned.gfa");
		fs::copy(_inputDir + "/minimizer_graph_cleaned.gfa", dir + "/minimizer_graph_cleaned.gfa");
		fs::copy(_inputDir + "/parameters.gz", dir + "/parameters.gz");
		fs::copy(_inputDir + "/mdbg_nodes.gz", dir + "/mdbg_nodes.gz");
		//fs::copy(_inputDir + "/mdbg_nodes_init.gz", dir + "/mdbg_nodes_init.gz", copyOptions);
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

};	


#endif 



