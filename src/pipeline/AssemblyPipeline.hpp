

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


	AssemblyPipeline(): Tool (){
	}

    void execute (){
		execute_pipeline();
	}

	void parseArgs(int argc, char* argv[]){

		cxxopts::Options options("AssemblyPipeline", "");
		options.add_options()
		//("d,debug", "Enable debugging") // a bool parameter
		(ARG_INPUT_FILENAME, "", cxxopts::value<string>())
		(ARG_OUTPUT_DIR, "", cxxopts::value<string>())
		(ARG_INPUT_FILENAME_TRUTH, "", cxxopts::value<string>()->default_value(""));
		//(ARG_KMINMER_LENGTH, "", cxxopts::value<int>()->default_value("3"))
		//(ARG_MINIMIZER_LENGTH, "", cxxopts::value<int>()->default_value("21"))
		//(ARG_MINIMIZER_DENSITY, "", cxxopts::value<float>()->default_value("0.005"));
		//("k,kminmerSize", "File name", cxxopts::value<std::string>())
		//("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
		//;

		if(argc <= 1){
			std::cout << options.help() << std::endl;
			exit(0);
		}


		cxxopts::ParseResult result;

		try{
			result = options.parse(argc, argv);

			//_kminmerSize = result[ARG_KMINMER_LENGTH].as<int>(); //getInput()->getInt(STR_KMINMER_SIZE);
			_inputFilename = result[ARG_INPUT_FILENAME].as<string>(); //getInput()->getStr(STR_INPUT);
			//_minimizerSize = result[ARG_MINIMIZER_LENGTH].as<int>(); //getInput()->getInt(STR_MINIM_SIZE);
			_inputDir = result[ARG_OUTPUT_DIR].as<string>(); //getInput()->getStr(STR_OUTPUT);
			_truthInputFilename = result[ARG_INPUT_FILENAME_TRUTH].as<string>();
			//_minimizerDensity = result[ARG_MINIMIZER_DENSITY].as<float>(); //getInput()->getDouble(STR_DENSITY);

		}
		catch (const std::exception& e){
			std::cout << options.help() << std::endl;
			std::cerr << e.what() << std::endl;
			std::exit(EXIT_FAILURE);
		}


        fs::path path(_inputDir);
	    if(!fs::exists (path)) fs::create_directory(path); 

		//T exception{text};




		//_inputFilename = result[STR_INPUT].as<string>();

		//_compositionManager = new CompositionManager(4);

		//int w = 80;
		


		//if(!System::file().doesExist (_outputDir)) System::file().mkdir(_outputDir, -1);
		

		cout << endl;
		cout << "Input: " << _inputFilename << endl;
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
		float density = 0.005;
		u_int16_t minimizerSize = 21;

		string command = "";

		writeParameters(minimizerSize, 4, density);

		//Read selection
		command = _filename_exe + " readSelection -i " + _inputFilename + " -o " + _inputDir;
		executeCommand(command);




		u_int64_t pass = 0;

		for(size_t k=4; k<5; k+=1){
			//cout << "Start asm: " << k << endl;

			writeParameters(minimizerSize, k, density);


			command = _filename_exe + " graph -o " + _inputDir;
			if(pass > 0) command += " -c " +  _inputDir + "/contigs.min.gz";
			executeCommand(command);
			//getchar();

			command = _filename_exe + " toBasespace -i " + _inputFilename +  " -o " + _inputDir;
			if(pass > 0) command += " -c " +  _inputDir + "/contigs.min.gz";
			executeCommand(command);

			command = _filename_exe + " asm -o " + _inputDir;
			if(_truthInputFilename != "") command += " --itruth " + _truthInputFilename;
			if(pass > 0) command += " -c " +  _inputDir + "/contigs.min.gz";
			executeCommand(command);


			//getchar();
			pass += 1;
			//if(pass == 2) break;
			//if(k >= 21) break;
			//exit(1);
		}

    }

	void executeCommand(const string& command){
		cout << command << endl;

		int ret = system(command.c_str());
		if(ret != 0){
			cerr << "Command failed: " << ret << endl;
			exit(ret);
		}
	}

	void writeParameters(size_t minimizerSize, size_t k, float density){

		string filename_parameters = _inputDir + "/parameters.gz";
		gzFile file_parameters = gzopen(filename_parameters.c_str(),"wb");
		gzwrite(file_parameters, (const char*)&minimizerSize, sizeof(minimizerSize));
		gzwrite(file_parameters, (const char*)&k, sizeof(k));
		gzwrite(file_parameters, (const char*)&density, sizeof(density));
		gzclose(file_parameters);
	}

};	


#endif 



