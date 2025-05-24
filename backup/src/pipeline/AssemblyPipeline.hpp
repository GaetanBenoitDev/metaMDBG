

#ifndef MDBG_METAG_ASSEMBLYPIPELINE
#define MDBG_METAG_ASSEMBLYPIPELINE

#include "../Commons.hpp"


class AssemblyPipeline : public Tool{
    
public:

	string _inputFilename;
	string _outputDir;
	//float _minimizerDensity_correction;
	//float _minimizerDensity_assembly;
    size_t _minimizerSize;
    size_t _kminmerSize;
	string _filename_readMinimizers;
	string _filename_exe;
	string _truthInputFilename;
	int _nbCores;
	//bool _useInitialKminmerFilter;

	size_t _firstK;
	size_t _lastK;
	u_int64_t _maxK;

	string _inputFilenameComplete;
	size_t _meanReadLength;
	string _tmpDir;

	//float _minReadQuality;
	//u_int64_t _contigPolishing_nbReadFragments;
	//u_int64_t _strainPurging_minContigLength;
	//double _strainPurging_minIdentity;
	string _usedInputArgument;

	enum DataType{
		HiFi,
		Nanopore,
	};

	struct Params{
		DataType _dataType;
		float _minReadQuality;
		u_int32_t _readCorrectionMinOverlapLength;
		float _readCorrectionMinIdentity;
		float _contigDereplicationIdentityThreshold;
		float _minimizerDensityCorrection;
		float _minimizerDensityAssembly;
		u_int16_t _usedCoverageForContigPolishing;
		bool _useHomopolymerCompression;
		bool _useReadCorrection;
		string _minimap2Params;
	};

	struct ContigStats{
		u_int64_t _assemblySize;
		u_int32_t _n50;
		vector<u_int32_t> _contigLengths;
		u_int64_t _nbContigs;
		u_int64_t _nbContigsLonger1M;
		u_int64_t _nbContigsLonger1MCircular;

		ContigStats(){
			_assemblySize = 0;
			_n50 = 0;
			_nbContigs = 0;
			_nbContigsLonger1M = 0;
			_nbContigsLonger1MCircular = 0;
		}
	};

	Params _params;

	AssemblyPipeline(): Tool (){
	}

	void parseArgs(int argc, char* argv[]){

		//cout << argc << " " <<  string(*argv) << endl;
		//getchar();

		//cout << "Attention: dereplication set to 0.97" << endl;

		args::ArgumentParser parser("asm", "Example 1: " + string(argv[0]) +" asm --out-dir ./outputDir/ --in-hihi reads.fastq.gz --threads 4 \nExample 2: " + string(argv[0]) + " asm --out-dir ./outputDir/ --in-ont reads_A.fastq.gz reads_B.fastq.gz reads_C.fastq.gz --threads 4"); //"This is a test program.", "This goes after the options."
		args::Group groupInputOutput(parser, "Basic options:");
		args::Group groupAssembly(parser, "Assembly options:");
		args::Group groupCorrection(parser, "Correction options:");

		
		args::ValueFlag<std::string> arg_outputDir(groupInputOutput, "", "Output dir for contigs and temporary files", {ARG_OUTPUT_DIR2});
		args::NargsValueFlag<std::string> arg_readFilenames_hifi(groupInputOutput, "", "PacBio HiFi read filename(s) (separated by space)", {ARG_INPUT_HIFI}, 2);
		args::NargsValueFlag<std::string> arg_readFilenames_nanopore(groupInputOutput, "", "Nanopore R10.4+ read filename(s) (separated by space)", {ARG_INPUT_NANOPORE}, 2);
		args::ValueFlag<int> arg_nbCores(groupInputOutput, "", "Number of cores", {ARG_NB_CORES2}, NB_CORES_DEFAULT_INT);

		args::ValueFlag<int> arg_minimizerSize(groupAssembly, "", "k-mer size", {ARG_MINIMIZER_LENGTH2}, 13);
		args::ValueFlag<float> arg_densityAssembly(groupAssembly, "", "Fraction of total k-mers used for assembly", {ARG_MINIMIZER_DENSITY_ASSEMBLY}, 0.005f);
		args::ValueFlag<int> arg_maxK(groupAssembly, "", "Stop assembly after k iterations", {ARG_MAXK}, 0);
		
		args::ValueFlag<float> arg_minReadQuality(groupCorrection, "", "Minimum read average quality", {ARG_MIN_READ_QUALITY}, 0.0);
		args::ValueFlag<float> arg_densityCorrection(groupCorrection, "", "Fraction of total k-mers used for correction", {ARG_MINIMIZER_DENSITY_CORRECTION}, 0.025f);
		args::ValueFlag<float> arg_readCorrectionMinIdentity(groupCorrection, "", "Min read identity", {"min-read-identity"}, 0.96);
		args::ValueFlag<int> arg_readCorrectionMinOverlapLength(groupCorrection, "", "Min read overlap length", {"min-read-overlap"}, 1000);
		//args::ValueFlag<int> arg_nbWindows(groupAdvanced, "", "Maximum read coverage used for contig correction (increase for better correction)", {ARG_NB_WINDOWS}, 100);
		//args::ValueFlag<int> arg_length(groupAdvanced, "", "Use contigs with length > l as references for strain duplication purging", {ARG_LINEAR_LENGTH}, 1000000);
		//args::ValueFlag<float> arg_minIdentity(groupAdvanced, "", "Minimum identity for strain purging (0-1)", {ARG_MIN_IDENTITY}, 0.99);
		//args::Flag arg_bf(groupAdvanced, "", "Disable unique kminmer filter prior to graph construction", {ARG_BLOOM_FILTER});
		//args::Flag arg_homopolymerCompression(parser, "", "Activate read homopolymer compression", {ARG_HOMOPOLYMER_COMPRESSION});
		//args::Flag arg_correction(parser, "", "Activate read correction", {ARG_CORRECTION});
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
			std::cerr << parser;
			exit(0);
		}
		catch (const std::exception& e)
		{
			std::cerr << parser;
			//cout << endl;
			std::cerr << e.what() << endl;
			exit(0);
		}

		if(arg_help){
			std::cerr << parser;
			exit(0);
		}

		if(!arg_outputDir){
			std::cerr << parser;
			cerr << " Argument --out-dir is required" << endl;
			exit(0);
		}


		if(!arg_readFilenames_hifi && !arg_readFilenames_nanopore){
			std::cerr << parser;
			cerr << " One of the arguments ";
			for(const string& arg : possibleInputArguments){
				cerr << arg + " ";
			}
			cerr << "is required" << endl;
			exit(0);
		}

		if(arg_readFilenames_hifi && arg_readFilenames_nanopore){
			std::cerr << parser;
			cerr << " Choose only one of the arguments ";
			for(const string& arg : possibleInputArguments){
				cerr << arg + " ";
			}
			cerr << endl;
			exit(0);
		}

		_outputDir = args::get(arg_outputDir);
		_tmpDir = _outputDir + "/tmp/";
		_minimizerSize = args::get(arg_minimizerSize);
		//_minimizerDensity = args::get(arg_d);
		_nbCores = args::get(arg_nbCores);

		_maxK = 0;
		if(arg_maxK){
			_maxK = args::get(arg_maxK);
		}

		//_strainPurging_minContigLength = args::get(arg_length);
		//_strainPurging_minIdentity = args::get(arg_minIdentity);
		//_strainPurging_minIdentity = min(_strainPurging_minIdentity, 0.999); //Some duplication is required to remove erroneous contigs

		//_useInitialKminmerFilter = true;
		//if(arg_bf){
		//	_useInitialKminmerFilter = false;
		//}



	    if(!fs::exists (_outputDir)) fs::create_directories(_outputDir); 
	    if(!fs::exists (_tmpDir)) fs::create_directories(_tmpDir); 
	    if(!fs::exists (_tmpDir + "/filter")) fs::create_directories(_tmpDir + "/filter"); 
		
		
		//string bannedDir = _tmpDir + "/bannedKminmers/";
	    //if(fs::exists (bannedDir)){
		//	fs::remove_all(bannedDir);
		//}
		//fs::create_directories(bannedDir); 

		_inputFilename = _tmpDir + "/input.txt";

		if(arg_readFilenames_hifi){
			_usedInputArgument = ARG_INPUT_HIFI;
			Commons::createInputFile(args::get(arg_readFilenames_hifi), _inputFilename);
			_params = getParamsHifi();
		}
		else if(arg_readFilenames_nanopore){
			_usedInputArgument = ARG_INPUT_NANOPORE;
			Commons::createInputFile(args::get(arg_readFilenames_nanopore), _inputFilename);
			_params = getParamsNanopore();
		}

		_params._minimizerDensityCorrection = args::get(arg_densityCorrection);
		_params._minimizerDensityAssembly = args::get(arg_densityAssembly);
		_params._readCorrectionMinOverlapLength = args::get(arg_readCorrectionMinOverlapLength);
		_params._readCorrectionMinIdentity = args::get(arg_readCorrectionMinIdentity);
		_params._minReadQuality = args::get(arg_minReadQuality);
		//_contigPolishing_nbReadFragments = args::get(arg_nbWindows);




		_filename_exe = argv[0];

		
	}

	Params getParamsHifi(){

		Params params;

		params._dataType = DataType::HiFi;
		params._readCorrectionMinIdentity = 0.99;
		params._readCorrectionMinOverlapLength = 1000;
		params._minReadQuality = 0;
		params._contigDereplicationIdentityThreshold = 0.99;
		params._useHomopolymerCompression = true;
		params._useReadCorrection = false;
		params._usedCoverageForContigPolishing = 50;
		params._minimap2Params = "-x map-hifi";

		return params;
	}

	Params getParamsNanopore(){

		Params params;
		
		params._dataType = DataType::Nanopore;
		params._readCorrectionMinIdentity = 0.96;
		params._readCorrectionMinOverlapLength = 1000;
		params._minReadQuality = 0;
		params._contigDereplicationIdentityThreshold = 0.99;
		params._useHomopolymerCompression = false;
		params._useReadCorrection = true;
		params._usedCoverageForContigPolishing = 100;
		params._minimap2Params = "-x map-ont";

		return params;

	}

    void execute (){

		
		if(fs::exists(_outputDir + "/metaMDBG.log")){
			fs::remove(_outputDir + "/metaMDBG.log");
		}

		//const string& failedFilename = _tmpDir + "/failed.txt";
		//if(fs::exists(failedFilename)){
		//	fs::remove(_tmpDir + "/failed.txt");
		//}

		openLogFile(_tmpDir);

		checkDependencies();

		Logger::get().info() << "";
		//Logger::get().debug() << "Input: " << _inputFilename;
		Logger::get().info() << "Output dir:              " << _outputDir;
		Logger::get().info() << "Minimizer length:        " << _minimizerSize;
		Logger::get().info() << "Density correction:      " << _params._minimizerDensityCorrection;
		Logger::get().info() << "Density assembly:        " << _params._minimizerDensityAssembly;
		Logger::get().info() << "Min read quality:        " << _params._minReadQuality;
		Logger::get().info() << "Min read overlap length: " << _params._readCorrectionMinOverlapLength;
		Logger::get().info() << "Min read identity:       " << _params._readCorrectionMinIdentity;
		Logger::get().info() << "Contig Derep threshold:  " << _params._contigDereplicationIdentityThreshold;
		Logger::get().info() << "Homopolymer compression: " << _params._useHomopolymerCompression;
		Logger::get().info() << "Read correction:         " << _params._useReadCorrection;
		Logger::get().info() << "";
		
		auto start = high_resolution_clock::now();

		//for(size_t i=0; i<1000; i++)
		execute_pipeline();

		auto stop = high_resolution_clock::now();


		ContigStats contigStats = computeContigStats(_outputDir + "/contigs.fasta.gz");

		//cout << "Duration: " << duration_cast<seconds>(stop - start).count() << "s" << endl;
		//ofstream outfile(_tmpDir + "/perf.txt", std::ios_base::app);
		//outfile << "Run time (hh:mm:ss): " << Utils::formatTimeToHumanReadable(duration_cast<seconds>(stop - start)) << endl;
		//outfile.close();

		ifstream perfFile(_tmpDir + "/perf.txt");
		std::ostringstream ss;
		ss << perfFile.rdbuf();
		std::string maxMemory = ss.str();

		cleanTmpFolder();
		
		Logger::get().info() << "";
		Logger::get().info() << "\tRun time:                   " << Utils::formatTimeToHumanReadable(duration_cast<seconds>(stop - start)); //<< " (hh:mm:ss)";
		Logger::get().info() << "\tPeak memory:                " << maxMemory << " GB";
		Logger::get().info() << "\tAssembly length:            " << contigStats._assemblySize;
		Logger::get().info() << "\tContigs N50:                " << contigStats._n50;
		Logger::get().info() << "\tNb contigs:                 " << contigStats._nbContigs;
		Logger::get().info() << "\tNb Contigs (>1Mb):          " << contigStats._nbContigsLonger1M;
		Logger::get().info() << "\tNb circular contigs (>1Mb): " << contigStats._nbContigsLonger1MCircular;
		Logger::get().info() << "";
		Logger::get().info() << "Contig filename: " << _outputDir + "/contigs.fasta.gz";
		Logger::get().info() << "Done!";

		perfFile.close();
		
	}

	float getPeakMemory(){

		try
		{
			ifstream perfFile(_tmpDir + "/perf.txt");
			std::ostringstream ss;
			ss << perfFile.rdbuf();
			std::string maxMemory = ss.str();

			return stof(maxMemory);
		}
		catch (const args::Help&)
		{
			return 0;
		}


	}

	void cleanTmpFolder(){
		
		fs::remove_all(_tmpDir + "/filter");
		fs::remove(_tmpDir + "/assembly_graph.gfa.unitigs.nodepath");
		fs::remove(_tmpDir + "/contigs.nodepath");
		fs::remove(_tmpDir + "/kminmerData_min_init.txt");
		fs::remove(_tmpDir + "/kminmerData_min.txt");
		fs::remove(_tmpDir + "/minimizer_graph.gfa");
		fs::remove(_tmpDir + "/read_data_corrected.txt");
		//fs::remove(_tmpDir + "/read_data_init.txt");
		//fs::remove(_tmpDir + "/read_stats.txt");
		fs::remove(_tmpDir + "/unitig_data.txt");
		fs::remove(_tmpDir + "/unitigGraph.edges.predecessors.bin");
		fs::remove(_tmpDir + "/unitigGraph.edges.successors.bin");
		fs::remove(_tmpDir + "/unitigGraph.nodes.bin");
		fs::remove(_tmpDir + "/small_contigs.bin");
		fs::remove(_tmpDir + "/perf.txt");
		fs::remove(_tmpDir + "/time.txt");
		//fs::remove(_tmpDir + "/data.txt");
		fs::remove(_tmpDir + "/contig_data_backup.txt");
		//fs::remove(_tmpDir + "/contig_data.txt");
		
	}



    void execute_pipeline(){

		_firstK = 4;

		string command = "";

		if(_params._useReadCorrection){
			writeParameters(_minimizerSize, _firstK, _params._minimizerDensityCorrection, _firstK, _firstK, _lastK, _params._minimizerDensityCorrection, _params._useHomopolymerCompression, _params._dataType);
		}
		else{
			writeParameters(_minimizerSize, _firstK, _params._minimizerDensityAssembly, _firstK, _firstK, _lastK, _params._minimizerDensityCorrection, _params._useHomopolymerCompression, _params._dataType);
		}
		//createInputFile(false);

		Logger::get().info() << "Converting reads to minimizers";

		command = _filename_exe + " readSelection " + _tmpDir + " " + _tmpDir + "/read_data_init.txt" + " " + _inputFilename + " --threads " + to_string(_nbCores) + " --min-read-quality " + to_string(_params._minReadQuality);
		if(_params._useHomopolymerCompression) command += " --homopolymer-compression";
		if(_params._useReadCorrection) command += " --output-quality";
		executeCommand(command);
		


		//cout << "AssemblyPipeline: Setting density to 0.005" << endl;
		//_minimizerDensity = 0.005;
		//writeParameters(_minimizerSize, _firstK, _minimizerDensity, _firstK, _firstK, _lastK);

		//"reprise: il faut selectionner les small contigs, par abondance, si c'est un tip ou non etc"

		//float density = 0.005;
		//u_int16_t minimizerSize = 21;



		ifstream file_readStats(_tmpDir + "/read_stats.txt");
	
		u_int64_t nbReads;
		u_int32_t n50ReadLength;
		float actualDensity;
		u_int64_t nbBases;
		file_readStats.read((char*)&nbReads, sizeof(nbReads));
		file_readStats.read((char*)&n50ReadLength, sizeof(n50ReadLength));
		file_readStats.read((char*)&actualDensity, sizeof(actualDensity));
		file_readStats.read((char*)&nbBases, sizeof(nbBases));

		file_readStats.close();


		//u_int64_t meanReadLength = computeMeanReadLength(_inputFilename);
		_lastK = n50ReadLength*_params._minimizerDensityAssembly*2.0f; //1.2f; //2.0f; //*0.95
		_meanReadLength = n50ReadLength;
		//_lastK = 10;

		if(_maxK > 0){
			//int minimizerSize = _minimizerSize;
			//int maxAssemblyLength = _maxAssemblyLength;
			//int maxK = ((maxAssemblyLength-minimizerSize) * _params._minimizerDensityAssembly) + 1;;
			_lastK = _maxK; //max((u_int64_t)0, _maxK);
		}

		_lastK = max(_lastK, _firstK+1);

		//_lastK = 5;
		/*
		sort(_readLengths.begin(), _readLengths.end());
		_logFile << _readLengths[_readLengths.size() * 0.1] << endl;
		_logFile << _readLengths[_readLengths.size() * 0.2] << endl;
		_logFile << _readLengths[_readLengths.size() * 0.3] << endl;
		_logFile << _readLengths[_readLengths.size() * 0.4] << endl;
		_logFile << _readLengths[_readLengths.size() * 0.5] << endl;
		_logFile << _readLengths[_readLengths.size() * 0.6] << endl;
		_logFile << _readLengths[_readLengths.size() * 0.7] << endl;
		_logFile << _readLengths[_readLengths.size() * 0.8] << endl;
		_logFile << _readLengths[_readLengths.size() * 0.9] << endl;
		_readLengths.clear();
		//return scores[size * 0.1];
		*/

		/*
		//cerr << "Mean read length: " << meanReadLength << endl;
		cerr << "N50 read length: " << n50ReadLength << endl;
		cerr << "Min k: " << _firstK << endl;
		cerr << "Max k: " << _lastK << endl;
		cerr << endl;
		_logFile << "N50 read length: " << n50ReadLength << endl;
		_logFile << "Min k: " << _firstK << endl;
		_logFile << "Max k: " << _lastK << endl;
		_logFile << endl;
		*/
		Logger::get().info() << "";
		Logger::get().info() << "Total read bps:  " << nbBases;
		Logger::get().info() << "N50 read length: " << n50ReadLength;
		Logger::get().info() << "Start k:         " << _firstK;
		Logger::get().info() << "End k:           " << _lastK;
		Logger::get().info() << "";

		//if(params._dataType == DataType::HiFi){
		//	params._readCorrectionMinOverlapLength = n50ReadLength*0.8;
		//}
		//params._readCorrectionMinOverlapLength = max((u_int32_t)1000, n50ReadLength / 2);



		if(_params._useReadCorrection){
			Logger::get().info() << "Correcting reads";
			writeParameters(_minimizerSize, _firstK, _params._minimizerDensityAssembly, _firstK, _firstK, _lastK, _params._minimizerDensityCorrection, _params._useHomopolymerCompression, _params._dataType);
		
			command = _filename_exe + " readCorrection " + _tmpDir + " --min-identity " + to_string(_params._readCorrectionMinIdentity) + " --min-overlap-length " + to_string(_params._readCorrectionMinOverlapLength) + " --threads " + to_string(_nbCores);
			executeCommand(command);
		}


		ofstream fileSmallContigs(_tmpDir + "/small_contigs.bin");
		fileSmallContigs.close();
		//ofstream fileJoints(_tmpDir + "/joint_data.txt");
		//fileJoints.close();


		

		u_int64_t pass = 0;
		u_int32_t prevK = -1;

		//for(size_t k=_firstK; k<_lastK; k+=1){

		size_t k = _firstK;

		
		while(k < _lastK){


			//cerr << "Multi-k pass: " << k << "/" << _lastK << endl;
			//_logFile << "Multi-k pass: " << k << "/" << _lastK << endl;
			Logger::get().info() << "Multi-k pass: " << k << "/" << _lastK;
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
			//command = _filename_exe + " readSelection -i " + inputFilename + " -o " + _tmpDir + " -f " + _tmpDir + "/read_data.gz";
			//if(pass == 0) command += " --firstpass";
			//executeCommand(command);

			executePass(k, prevK, pass);

		
			//./bin/mdbgAsmMeta bin -o ~/workspace/run/overlap_test_multik_201/pass_0 -c ~/workspace/run/overlap_test_multik_201/contigs_10.fasta.gz

			//}
			
			//time ./bin/mdbgAsmMeta countKmer -i ~/workspace/data/AD/shortreads/input_10.txt -c /home/gats/workspace/run/overlap_test_multik_AD/contigs.nodepath.gz.fasta.gz.fasta.gz  -o ~/workspace/run/overlap_test_multik_AD/contig_coverages.tsv
			//!!!Binning: assembly2, permettre de choisir le nom du fichier de contig via -c, actuellement c'est contigs.fasta.gz
			//./bin/mdbgAsmMeta bin -o ~/workspace/run/overlap_test_multik_AD/ -a ~/workspace/run/overlap_test_multik_AD/contig_coverages.tsv



			prevK = k;
			//getchar();
			pass += 1;
			//if(pass == 2) break;
			//if(k >= 21) break;
			//exit(1);
			//break;
			//cout << "pass done" << endl;
			//if(generatedContigs) getchar();
			//getchar();
			//if(pass >= 2) break;
			//if(k > 30) getchar();

			k += Commons::getMultikStep(k);
		}

		Logger::get().info() << "Multi-k pass: " << _lastK << "/" << _lastK;
		//cerr << "Multi-k pass: " << _lastK << "/" << _lastK << endl;
		//_logFile << "Multi-k pass: " << _lastK << "/" << _lastK << endl;

		executePass(_lastK, prevK, pass);
		
		

		//execTest();
		Logger::get().debug() << "pass done";

    }

	/*
	void execTest(){

		writeParameters(_minimizerSize, _firstK, _minimizerDensity_assembly, _firstK, _firstK, _lastK);
		
			string contigFilenameDummy = _tmpDir + "/dummy.fasta.gz ";
			string contigFilenameCompressed = _tmpDir + "/contigs_uncorrected_underep.fasta.gz";

			//string command = _filename_exe + " toBasespaceFast " + " " + _tmpDir + " " + _tmpDir + "/contig_data.txt " + " " + contigFilenameDummy + " " + _inputFilename + " -t " + to_string(_nbCores);
			//executeCommand(command);

			string command = _filename_exe + " toBasespace " + " " + _tmpDir + " " + _tmpDir + "/contig_data.txt " + " " + contigFilenameCompressed + " " + _inputFilename  + " -t " + to_string(_nbCores);
			//executeCommand(command);

			command = _filename_exe + " derep " + contigFilenameCompressed + " " + contigFilenameDummy + " " + _tmpDir + " -t " + to_string(_nbCores) + " --nodump";
			//executeCommand(command);

			const string contigFilename_uncorrected = _tmpDir + "/contigs_uncorrected.fasta.gz";
			//dereplicate(contigFilenameCompressed, contigFilename_uncorrected);
			

			//command = _filename_exe + " circ " + _tmpDir + " -t " + to_string(_nbCores);
			//executeCommand(command);


			//cerr << "Constructing base-space contigs..." << endl;
			//_logFile << "Constructing base-space contigs..." << endl;

			//const string contigFilename_uncorrected = _tmpDir + "/contigs_uncorrected.fasta.gz";
			//command = _filename_exe + " toBasespace " + " " + _tmpDir + " " + _tmpDir + "/contig_data.txt " + " " + contigFilename_uncorrected + " " + _inputFilename  + " -t " + to_string(_nbCores);
			//executeCommand(command);

			cerr << "Polishing contigs..." << endl;
			_logFile << "Polishing contigs..." << endl;
			//./bin/metaMDBG polish ~/workspace/run/overlap_test_201/contigs_uncorrected.fasta.gz ~/workspace/run/overlap_test_201/ ~/workspace/data/overlap_test/genome_201_50x/simulatedReads_0.fastq.gz ~/workspace/data/overlap_test/genome_201_50x/simulatedReads_0.fastq.gz -t 15 --qual
			//getchar();
			command = _filename_exe + " polish " + contigFilename_uncorrected + " " + _tmpDir + " " + " -t " + to_string(_nbCores) + " -n " + to_string(_contigPolishing_nbReadFragments) + " --metaMDBG"; //--circ
			executeCommand(command);
			//generatedContigs = true;

			cerr << "Purging strain duplication..." << endl;
			_logFile << "Polishing contigs..." << endl;
			string polishedContigFilename = _tmpDir + "/contigs_polished.fasta.gz";
			string polishedContigFilenamederep = _outputDir + "/contigs.fasta.gz ";
			command = _filename_exe + " derep " + polishedContigFilename + " " + polishedContigFilenamederep + " " + _tmpDir + " -t " + to_string(_nbCores) + " -l " + to_string(_strainPurging_minContigLength) + " -i " + to_string(_strainPurging_minIdentity);
			executeCommand(command);
	}
	*/

	void executePass(size_t k, size_t prevK, size_t pass){

		bool isFinalPass = k == _lastK;
		bool isFirstPass = k == _firstK;

		writeParameters(_minimizerSize, k, _params._minimizerDensityAssembly, _firstK, prevK, _lastK, _params._minimizerDensityCorrection, _params._useHomopolymerCompression, _params._dataType);

		string command = "";

		if(pass == 0){
			command = _filename_exe + " graph " + _tmpDir + " --threads " + to_string(_nbCores);
			//if(!_useInitialKminmerFilter){
			//	command += " --nofilter ";
			//}
		}
		else{
			command = _filename_exe + " graph " + _tmpDir + " --threads " + to_string(_nbCores);	
		}
		if(pass == 0) command += " --firstpass ";
		if(_params._useReadCorrection) command += " --corrected-read ";
		executeCommand(command);
		//getchar();

		//command = _filename_exe + " multik -o " + _tmpDir + " -t " + to_string(_nbCores);
		//if(!_truthInputFilename.empty()) command += " --itruth " + _truthInputFilename;
		//if(pass > 0) command += " -c " +  _tmpDir + "/contig_data.gz";
		//executeCommand(command);

		//Generate contigs


		//getchar();
		//getchar();


		/*
		command = _filename_exe + " toMinspace " + " -o " + _tmpDir + " -c " + _tmpDir + "/unitigs.nodepath.gz" + " -f " + _tmpDir + "/unitig_data.txt";
		executeCommand(command);
		*/

		//bool generatedContigs = false;
		//if(k == 5 || k == 10 || k == 16 || k == 21 || k == 26 || k == 31){
		if(isFinalPass){

			
			//Generate contigs
			command = _filename_exe + " contig " + " " + _tmpDir + " " + " --threads " + to_string(_nbCores);
			if(!_truthInputFilename.empty()) command += " --itruth " + _truthInputFilename;
			//if(pass == 0) command += " --firstpass";
			executeCommand(command);

			//getchar();

			command = _filename_exe + " toMinspace " + " " + _tmpDir + " " + _tmpDir + "/contigs.nodepath" + " " + _tmpDir + "/contig_data.txt --threads "+ to_string(_nbCores);
			executeCommand(command);

			

			//command = _filename_exe + " toMinspace " + " " + _tmpDir + " " + _tmpDir + "/contigs.nodepath" + " " + _tmpDir + "/contig_data.txt";
			//executeCommand(command);

			//cerr << "Removing overlaps and duplication..." << endl;
			//_logFile << "Removing overlaps and duplication..." << endl;
			Logger::get().info() << "Removing overlaps and duplication";

			//cout << "AssemblyPipeline: Disabled small contigs append" << endl;
			appendSmallContigs();
			
			fs::copy(_tmpDir + "/contig_data.txt", _tmpDir + "/contig_data_backup.txt", fs::copy_options::overwrite_existing);

			string contigFilenameDummy = _tmpDir + "/dummy.fasta.gz ";
			string contigFilenameCompressed = _tmpDir + "/contigs_uncorrected.fasta.gz";

			string command = _filename_exe + " toBasespaceFast " + " " + _tmpDir + " " + _tmpDir + "/contig_data.txt " + " " + contigFilenameDummy + " " + _inputFilename + " --threads " + to_string(_nbCores);
			executeCommand(command);


			Logger::get().info() << "Constructing base-space contigs";

			if(_params._dataType == DataType::HiFi){
				command = _filename_exe + " toBasespace_hifi " + " " + _tmpDir + " " + _tmpDir + "/contig_data.txt " + " " + contigFilenameCompressed + " " + _inputFilename  + " --threads " + to_string(_nbCores);
				//if(_params._useHomopolymerCompression) command += " --homopolymer-compression";
				executeCommand(command);
			}
			else if(_params._dataType == DataType::Nanopore){
				command = _filename_exe + " toBasespace_ont " + " " + _tmpDir + " " + _tmpDir + "/contig_data.txt " + " " + contigFilenameCompressed + " " + _inputFilename  + " --threads " + to_string(_nbCores);
				//if(_params._useHomopolymerCompression) command += " --homopolymer-compression";
				executeCommand(command);
			} 


			//command = _filename_exe + " derep " + contigFilenameCompressed + " " + contigFilenameDummy + " " + _tmpDir + " -t " + to_string(_nbCores) + " --nodump";
			//executeCommand(command);

			//const string contigFilename_uncorrected = _tmpDir + "/contigs_uncorrected.fasta.gz";
			//dereplicate(contigFilenameCompressed, contigFilename_uncorrected);
			

			//command = _filename_exe + " circ " + _tmpDir + " -t " + to_string(_nbCores);
			//executeCommand(command);

			/*

		
			ofstream contigInputFile(_tmpDir + "/input_contig.txt");
			contigInputFile << contigFilenameCompressed_derep << endl;
			contigInputFile.close();
			command = _filename_exe + " readSelection -i " + _tmpDir + "/input_contig.txt" + " -o " + _tmpDir + " -f " + _tmpDir + "/contig_data.txt" + " -t " + to_string(_nbCores) + " --contig";
			executeCommand(command);
			*/
		
			//cerr << "Constructing base-space contigs..." << endl;
			//_logFile << "Constructing base-space contigs..." << endl;

			//const string contigFilename_uncorrected = _tmpDir + "/contigs_uncorrected.fasta.gz";
			//command = _filename_exe + " toBasespace " + " " + _tmpDir + " " + _tmpDir + "/contig_data.txt " + " " + contigFilename_uncorrected + " " + _inputFilename  + " -t " + to_string(_nbCores);
			//executeCommand(command);


			//cerr << "Polishing contigs..." << endl;
			//_logFile << "Polishing contigs..." << endl;
			//./bin/metaMDBG polish ~/workspace/run/overlap_test_201/contigs_uncorrected.fasta.gz ~/workspace/run/overlap_test_201/ ~/workspace/data/overlap_test/genome_201_50x/simulatedReads_0.fastq.gz ~/workspace/data/overlap_test/genome_201_50x/simulatedReads_0.fastq.gz -t 15 --qual
			//getchar();
			

			//if(maxMemory < 8){
			//	maxMemory = 8;
			//}

			Logger::get().info() << "Mapping reads vs contigs";

			string readPartitionDir = _tmpDir + "/_polish_readPartitions/";
			if(!fs::exists(readPartitionDir)){
				fs::create_directories(readPartitionDir);
			}

			//string readFilenames = "";
			//ReadParser readParser(_inputFilename_reads, false, false, _logFile);

			//for(const string& filename : readParser._filenames){
			//	readFilenames += filename + " ";
			//}

			int minimapBatchSize = 1;
			float peakMemory = getPeakMemory();
			if(peakMemory < 8 || peakMemory > 1000000){
				minimapBatchSize = 1;
			}
			else{
				minimapBatchSize = peakMemory / 8;
			}
			if(minimapBatchSize < 0){
				minimapBatchSize = 1;
			}
			if(minimapBatchSize > 100){
				minimapBatchSize = 100;
			}

			
			command = "minimap2 -I " + to_string(minimapBatchSize) + "G -t " + to_string(_nbCores) + " " + _params._minimap2Params + " " + contigFilenameCompressed + " " + Commons::inputFileToFilenames(_tmpDir + "/input.txt");
			command += " | gzip -c - > " + readPartitionDir + "/polish_mapping.paf.gz";
			Utils::executeCommand(command, _tmpDir);


			Logger::get().info() << "Polishing contigs";

			command = _filename_exe + " polish " + contigFilenameCompressed + " " + _tmpDir + " " + " --threads " + to_string(_nbCores) + " -n " + to_string(_params._usedCoverageForContigPolishing) + " --metaMDBG"; //--circ
			executeCommand(command);
			//generatedContigs = true;



			Logger::get().info() << "Mapping contigs vs contigs";
			string polishedContigFilename = _tmpDir + "/contigs_polished.fasta.gz";

			command = "minimap2 -X -I " + to_string(minimapBatchSize) + "G -t " + to_string(_nbCores) + " " + polishedContigFilename + " " + polishedContigFilename;
			command += " | gzip -c - > " + _tmpDir + "/_tmp_mapping_derep__.paf.gz";
			Utils::executeCommand(command, _tmpDir);


			Logger::get().info() << "Purging strain duplication";
			//cerr << "Purging strain duplication..." << endl;
			//_logFile << "Polishing contigs..." << endl;
			string polishedContigFilenamederep = _outputDir + "/contigs.fasta.gz ";
			command = _filename_exe + " derep " + polishedContigFilename + " " + polishedContigFilenamederep + " " + _tmpDir + " --threads " + to_string(_nbCores) + " -i " + to_string(_params._contigDereplicationIdentityThreshold); // + " -l " + to_string(_strainPurging_minContigLength) + " -i " + to_string(_strainPurging_minIdentity);
			executeCommand(command);

		}
		else{

			command = _filename_exe + " contig " + " " + _tmpDir + " --threads " + to_string(_nbCores);;
			if(!_truthInputFilename.empty()) command += " --itruth " + _truthInputFilename;
			//if(pass == 0) command += " --firstpass";
			executeCommand(command);
			
			//cout << "done" << endl;
			//getchar();
			

			//if(k > _firstK){
			//	const auto copyOptions = fs::copy_options::overwrite_existing;
			//	fs::copy(_tmpDir + "/unitig_data.txt", _tmpDir + "/unitig_data_prev.txt", copyOptions);
			//}

			command = _filename_exe + " toMinspace " + " " + _tmpDir + " " + _tmpDir + "/contigs.nodepath" + " " + _tmpDir + "/unitig_data.txt --threads " + to_string(_nbCores);
			executeCommand(command);

			if(!isFirstPass){
				command = _filename_exe + " toMinspace " + " " + _tmpDir + " " + _tmpDir + "/assembly_graph.gfa.unitigs.nodepath" + " " + _tmpDir + "/assembly_graph.gfa.unitigs --threads " + to_string(_nbCores);
				executeCommand(command);
			}

			
			//if(pass == 0){
			//	fs::copy(_tmpDir + "/unitig_data.txt", _tmpDir + "/unitig_data_first.txt", fs::copy_options::overwrite_existing);
			//}

			//const string contigFilename_uncorrected = _tmpDir + "/contigs_uncorrected.fasta.gz";
			//command = _filename_exe + " toBasespace " + " " + _tmpDir + " " + _tmpDir + "/assembly_graph.gfa.unitigs " + " " + _tmpDir + "/unitigs.fasta.gz " + _inputFilename  + " -t " + to_string(_nbCores);
			//if(pass == 0) command += " --firstpass";
			//executeCommand(command);

			//if(pass > 0) exit(1);
		}	
			
		savePassData(k);
		//cout << "done" << endl;
		//getchar();
	}

	unordered_set<u_int32_t> _duplicatedContigIndex;
	ofstream _outputContigFileDerep;
	gzFile _outputContigFile;

	void dereplicate(const string& inputContigFilename, const string& outputContigFilename){

		ifstream infile(_tmpDir + "/duplicatedContigs.txt");

		string line;

		while (infile >> line){
			//_logFile << "Duplicate: " << line << endl;
			_duplicatedContigIndex.insert(stoull(line));
		}

		infile.close();


		string contigFilename = _tmpDir + "/contig_data_derep.txt";
		_outputContigFileDerep = ofstream(contigFilename);

		KminmerParser parser(_tmpDir + "/contig_data.txt", _minimizerSize, _kminmerSize, false, false);
		auto fp = std::bind(&AssemblyPipeline::dumpDereplicatedContigs_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
		parser.parseSequences(fp);

		_outputContigFileDerep.close();

		fs::remove(_tmpDir + "/contig_data.txt");
		fs::rename(_tmpDir + "/contig_data_derep.txt", _tmpDir + "/contig_data.txt");

		_outputContigFile = gzopen(outputContigFilename.c_str(),"wb");

		ReadParserParallel readParser(inputContigFilename, true, false, 1);
		readParser.parse(ReadPartitionFunctor(*this));

		gzclose(_outputContigFile);

		_duplicatedContigIndex.clear();
	}


	void dumpDereplicatedContigs_read(const vector<MinimizerType>& readMinimizers, u_int8_t isCircular, u_int64_t readIndex){
		

		//_logFile << "contig: " << readIndex << " " << readMinimizers.size() << endl;
		if(readMinimizers.size() == 0){
			cout << "Read with no minimizer after derep" << endl;
		}
		if(_duplicatedContigIndex.find(readIndex) != _duplicatedContigIndex.end()) return;

		//_logFile << "dump contigs: " << readIndex << " " << readMinimizers.size() << endl;
		u_int32_t contigSize = readMinimizers.size();
		_outputContigFileDerep.write((const char*)&contigSize, sizeof(contigSize));
		_outputContigFileDerep.write((const char*)&isCircular, sizeof(isCircular));
		_outputContigFileDerep.write((const char*)&readMinimizers[0], contigSize*sizeof(MinimizerType));

		//_nbContigsPost += 1;
		//cout << "Dump: " << readIndex << " " << readMinimizers.size() << endl;

		//u_int64_t s = 0;
		//for(size_t i=0; i<readMinimizers.size(); i++){
		//	s += (readMinimizers[i]);
		//}

		//_checksum += (s*readMinimizers.size());

	}


	class ReadPartitionFunctor {

		public:

		AssemblyPipeline& _parent;


		ReadPartitionFunctor(AssemblyPipeline& parent) : _parent(parent){
		}

		ReadPartitionFunctor(const ReadPartitionFunctor& copy) : _parent(copy._parent){
		}

		~ReadPartitionFunctor(){
		}

		void operator () (const Read& read) {

			if(_parent._duplicatedContigIndex.find(read._index) != _parent._duplicatedContigIndex.end()) return;

			string header = ">" + read._header + '\n';
			gzwrite(_parent._outputContigFile, (const char*)&header[0], header.size());
			string seq = read._seq + '\n';
			gzwrite(_parent._outputContigFile, (const char*)&seq[0], seq.size());
			
		}
	};



	ofstream _fileContigsAppend;
	u_int64_t _nbSmallContigs;

	void appendSmallContigs(){

		_nbSmallContigs = 0;
		Logger::get().debug() << "Append small contigs";

		string contigFilename = _tmpDir + "/contig_data.txt";

		_fileContigsAppend = ofstream(contigFilename, std::ios_base::app);

		KminmerParser parser(_tmpDir + "/small_contigs.bin", _minimizerSize, _kminmerSize, false, false);
		auto fp = std::bind(&AssemblyPipeline::appendSmallContigs_read, this, std::placeholders::_1, std::placeholders::_2);
		parser.parseSequences(fp);

		_fileContigsAppend.close();

		//cout << "Nb small contigs: " << _nbSmallContigs << endl;
	}

	void appendSmallContigs_read(const vector<MinimizerType>& readMinimizers, u_int64_t readIndex){
		
		u_int32_t contigSize = readMinimizers.size();
		_fileContigsAppend.write((const char*)&contigSize, sizeof(contigSize));
		
		u_int8_t isCircular = CONTIG_LINEAR;
		_fileContigsAppend.write((const char*)&isCircular, sizeof(isCircular));
		_fileContigsAppend.write((const char*)&readMinimizers[0], contigSize*sizeof(MinimizerType));

		_nbSmallContigs += 1;
	}

	void savePassData(u_int64_t k){


		bool isFirstPass = k == _firstK;
		if(isFirstPass) return;

		const string& dir = _tmpDir + "/pass_k" + to_string(k);

		if(fs::exists(dir)){
			fs::remove_all(dir);
		}

		fs::create_directories(dir);

		//const auto copyOptions = fs::copy_options::overwrite_existing;
		//fs::copy(_tmpDir + "/read_data.txt", dir + "/read_data.txt");


		//if(fs::exists(_tmpDir + "/read_index.txt")) fs::copy(_tmpDir + "/read_index.txt", dir + "/read_index.txt");
		//if(fs::exists(_tmpDir + "/minimizer_graph.gfa.gfa")) fs::copy(_tmpDir + "/minimizer_graph.gfa.gfa", dir + "/minimizer_graph.gfa.gfa");
		//if(fs::exists(_tmpDir + "/minimizer_graph_u.gfa")) fs::copy(_tmpDir + "/minimizer_graph_u.gfa", dir + "/minimizer_graph_u.gfa");
		if(fs::exists(_tmpDir + "/assembly_graph.gfa")) fs::rename(_tmpDir + "/assembly_graph.gfa", dir + "/assembly_graph.gfa");
		if(fs::exists(_tmpDir + "/assembly_graph.gfa.unitigs")) fs::rename(_tmpDir + "/assembly_graph.gfa.unitigs", dir + "/assembly_graph.gfa.unitigs");
		//if(fs::exists(_tmpDir + "/assembly_graph.gfa.unitigs.index")) fs::copy(_tmpDir + "/assembly_graph.gfa.unitigs.index", dir + "/assembly_graph.gfa.unitigs.index");
		//fs::copy(_tmpDir + "/contigs_path.csv", dir + "/contigs_path.csv");
		//fs::copy(_tmpDir + "/minimizer_graph_cleaned.gfa", dir + "/minimizer_graph_cleaned.gfa");
		fs::copy(_tmpDir + "/parameters.gz", dir + "/parameters.gz");

		//if(k == _firstK || k == _lastK){
			//fs::copy(_tmpDir + "/minimizer_graph.gfa", dir + "/minimizer_graph.gfa");
			//fs::copy(_tmpDir + "/kminmerData_min.txt", dir + "/kminmerData_min.txt");
			//if(fs::exists(_tmpDir + "/nodeName_to_unitigIndex.bin")) fs::copy(_tmpDir + "/nodeName_to_unitigIndex.bin", dir + "/nodeName_to_unitigIndex.bin");
			//if(fs::exists(_tmpDir + "/groundtruth_position.csv")) fs::copy(_tmpDir + "/groundtruth_position.csv", dir + "/groundtruth_position.csv");
			//if(fs::exists(_tmpDir + "/read_path.txt")) fs::copy(_tmpDir + "/read_path.txt", dir + "/read_path.txt");
			//if(fs::exists(_tmpDir + "/read_path_cleaned.txt")) fs::copy(_tmpDir + "/read_path_cleaned.txt", dir + "/read_path_cleaned.txt");
		//}
		//fs::copy(_tmpDir + "/mdbg_nodes_init.gz", dir + "/mdbg_nodes_init.gz", copyOptions);

		//fs::copy(_tmpDir + "/unitig_data.txt", dir + "/unitig_data.txt"); //Debug
	}



	void writeParameters(size_t minimizerSize, size_t k, float assemblyDensity, size_t firstK, size_t prevK, size_t lastK, float correctionDensity, bool useHomopolymerCompression, DataType dataType){


		float minimizerSpacingMean = 1 / assemblyDensity;
		float kminmerLengthMean = minimizerSpacingMean * (k-1);
		float kminmerOverlapMean = kminmerLengthMean - minimizerSpacingMean;

		string filename_parameters = _tmpDir + "/parameters.gz";
		gzFile file_parameters = gzopen(filename_parameters.c_str(),"wb");
		gzwrite(file_parameters, (const char*)&minimizerSize, sizeof(minimizerSize));
		gzwrite(file_parameters, (const char*)&k, sizeof(k));
		gzwrite(file_parameters, (const char*)&assemblyDensity, sizeof(assemblyDensity));
		gzwrite(file_parameters, (const char*)&firstK, sizeof(firstK));
		gzwrite(file_parameters, (const char*)&minimizerSpacingMean, sizeof(minimizerSpacingMean));
		gzwrite(file_parameters, (const char*)&kminmerLengthMean, sizeof(kminmerLengthMean));
		gzwrite(file_parameters, (const char*)&kminmerOverlapMean, sizeof(kminmerOverlapMean));
		gzwrite(file_parameters, (const char*)&prevK, sizeof(prevK));
		gzwrite(file_parameters, (const char*)&lastK, sizeof(lastK));
		gzwrite(file_parameters, (const char*)&_meanReadLength, sizeof(_meanReadLength));
		gzwrite(file_parameters, (const char*)&correctionDensity, sizeof(correctionDensity));
		gzwrite(file_parameters, (const char*)&useHomopolymerCompression, sizeof(useHomopolymerCompression));

		int dataTypeInt = 0;
		if(dataType == DataType::HiFi){
			dataTypeInt = 0;
		}
		else if(dataType == DataType::Nanopore){
			dataTypeInt = 1;
		}

		gzwrite(file_parameters, (const char*)&dataTypeInt, sizeof(dataTypeInt));

		gzclose(file_parameters);
	}


	/*
	void createInputFile(bool useContigs){
		_inputFilenameComplete = _tmpDir + "/input.txt";
		ofstream inputFile(_inputFilenameComplete);

		ReadParser readParser(_inputFilename, false);

		for(const string& filename : readParser._filenames){
			inputFile << filename << endl;
		}

		if(useContigs){
			inputFile << _tmpDir + "/tmpContigs.fasta.gz"  << endl;
		}

		inputFile.close();
	}

	string createInputFile(u_int32_t k){
		string inputFilename = _tmpDir + "/input.txt";
		ofstream inputFile(inputFilename);

		//const string& filename = _tmpDir + "/correctedReads_" + to_string(k) + ".min.gz.bitset";
		const string& filename = _tmpDir + "/read_data.txt.corrected.txt";
		inputFile << filename << endl;


		inputFile.close();

		return inputFilename;
	}
	*/

	/*
	vector<u_int64_t> _readLengths;
	u_int64_t _readLengthSum;
	u_int64_t _readLengthN;

	u_int64_t computeMeanReadLength(const string& filename){




		_readLengthSum = 0;
		_readLengthN = 0;

		auto fp = std::bind(&AssemblyPipeline::computeMeanReadLength_read, this, std::placeholders::_1);
		ReadParser readParser(filename, false, false, _logFile);
		readParser._maxReads = 100000;
		readParser.parse(fp);

		if(_readLengthN == 0){
			cerr << "No reads provided (check input filenames)" << endl;
			exit(1);
		}

		return _readLengthSum / _readLengthN;
	}
	
	void computeMeanReadLength_read(const Read& read){
		_readLengthSum += read._seq.size();
		_readLengthN += 1;
		_readLengths.push_back(read._seq.size());
	}
	*/

	void executeCommand(const string& command){

		//cout << command << endl;
		//string command2 = "time -v \"" + command + "\" 2>&1 " + _tmpDir + "/time.txt";

		Utils::executeCommand(command, _tmpDir);
		//getchar();
	}

	void checkDependencies(){

		string read = "";
		for(size_t i=0; i<10000; i++){
			read += "A";
		}

		string filename = _tmpDir + "/testDependencies.fasta";
		//string outputFilename_mapping = _tmpDir + "/testDependencies_align.paf";
		//string filenameBzip = _tmpDir + "/testDependencies_bgzip.fasta.gz";
		string outputFilename_mapping_minimap = _tmpDir + "/testDependencies_minimap.paf";

		if(fs::exists(filename)) fs::remove(filename);
		//if(fs::exists(filenameBzip)) fs::remove(filenameBzip);
		//if(fs::exists(filenameBzip + ".fai")) fs::remove(filenameBzip + ".fai");
		if(fs::exists(outputFilename_mapping_minimap)) fs::remove(outputFilename_mapping_minimap);

		ofstream file(filename);

		file << ">read1" << endl;
		file << read << endl;
		file << ">read2" << endl;
		file << read << endl;
		
		file.close();


		/*
		//cerr << "Checking dependencies: bgzip (wfmash package) ";
		//_logFile << "Checking dependencies: bgzip (wfmash package) ";
		string command = "cat " + filename + " | bgzip --threads " + to_string(_nbCores) + " -c > " + filenameBzip;
		Utils::executeCommand(command, _tmpDir, _logFile);

		if(!fs::exists(filenameBzip)){
			cerr << endl << "Dependency is missing: bgzip (wfmash)" << endl;
			exit(1);
		}
		else{
			cerr << "ok" << endl;
		}

		cerr << "Checking dependencies: samtools ";
		_logFile << "Checking dependencies: samtools ";
		command = "samtools faidx " + filenameBzip;
		Utils::executeCommand(command, _tmpDir, _logFile);

		if(!fs::exists(filenameBzip + ".fai")){
			cerr << endl << "Dependency is missing: samtools" << endl;
			exit(1);
		}
		else{
			cerr << "ok" << endl;
		}

		cerr << "Checking dependencies: wfmash ";
		_logFile << "Checking dependencies: wfmash ";
		command = "wfmash " + filenameBzip + " -t " + to_string(_nbCores) + " > " + outputFilename_mapping; //-l 5000 -p 80
		Utils::executeCommand(command, _tmpDir, _logFile);

		if(!fs::exists(outputFilename_mapping)){
			cerr << endl << "Dependency is missing: wfmash" << endl;
			exit(1);
		}
		else{
			cerr << "ok" << endl;
		}
		*/

		Logger::get().info() << "";
		Logger::get().info() << "Checking dependencies: minimap2";

		//cerr << "Checking dependencies: minimap2 ";
		//_logFile << "Checking dependencies: minimap2 ";
		string command = "minimap2 -x map-hifi " + filename + " " + filename + " -o " + outputFilename_mapping_minimap;
		Utils::executeCommand(command, _tmpDir);

		if(!fs::exists(outputFilename_mapping_minimap)){
			Logger::get().info() << "\n" << "Dependency is missing: minimap2.24+";
			exit(1);
		}
		else{
			Logger::get().info() << "ok";
		}

		//Logger::get().info() << "";

		if(fs::exists(filename)) fs::remove(filename);
		//if(fs::exists(filenameBzip)) fs::remove(filenameBzip);
		//if(fs::exists(filenameBzip + ".fai")) fs::remove(filenameBzip + ".fai");
		if(fs::exists(outputFilename_mapping_minimap)) fs::remove(outputFilename_mapping_minimap);

	}


	ContigStats _contigStats;
	//u_int64_t _assemblySize;
	//vector<u_int32_t> _contigLengths;

	ContigStats computeContigStats(const string& contigFilename){

		_contigStats = ContigStats{};
		ReadParserParallel readParser(contigFilename, true, false, 1);
		readParser.parse(ComputeContigStatsFunctor(*this));

		_contigStats._n50 = Utils::computeN50(_contigStats._contigLengths);

		return _contigStats;
	}

	class ComputeContigStatsFunctor{

		public:

		AssemblyPipeline& _parent;

		ComputeContigStatsFunctor(AssemblyPipeline& parent) : _parent(parent){
		}

		ComputeContigStatsFunctor(const ComputeContigStatsFunctor& copy) : _parent(copy._parent){
		}

		~ComputeContigStatsFunctor(){
		}

		void operator () (const Read& read) {

			//cout << _parent._contigStats._nbContigs << " " << read._header << endl;
			_parent._contigStats._assemblySize += read._seq.size();
			_parent._contigStats._contigLengths.push_back(read._seq.size());
			_parent._contigStats._nbContigs += 1;
			
			if(read._seq.size() > 1000000){
				_parent._contigStats._nbContigsLonger1M += 1;

				if(Utils::isContigCircular(read._header)){
					_parent._contigStats._nbContigsLonger1MCircular += 1;
				}
			}
		}
	};

};	



#endif 



