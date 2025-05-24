

#ifndef MDBG_METAG_REASSEMBLEPIPELINE
#define MDBG_METAG_REASSEMBLEPIPELINE

#include "../Commons.hpp"


class ReassemblePipeline : public Tool{
    
public:

	string _readFilename;
	string _checkmQualityReportFilename;
	string _outputDir;
	string _tmpDir;
	string _filename_exe;
	string _binDir;
	int _nbCores;

	vector<string>* _fields;
	vector<string>* _fields_optional;

	ReassemblePipeline(): Tool (){
	}

	void parseArgs(int argc, char* argv[]){

		args::ArgumentParser parser("asm", ""); //"This is a test program.", "This goes after the options."
		args::Positional<std::string> arg_outputDir(parser, "outputDir", "Output dir for results and temporary files", args::Options::Required);
		args::Positional<std::string> arg_binDir(parser, "binDir", "Directory containing bins", args::Options::Required);
		args::Positional<std::string> arg_checkm2QualityReportFilename(parser, "checkm2Report", "checkm2 report file", args::Options::Required);
		args::PositionalList<std::string> arg_readFilenames(parser, "reads", "Read filename(s) (separated by space)", args::Options::Required);
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

		_outputDir = args::get(arg_outputDir);
		_tmpDir = _outputDir + "/tmp/";
		_nbCores = args::get(arg_nbCores);
		_checkmQualityReportFilename = args::get(arg_checkm2QualityReportFilename);
		_binDir = args::get(arg_binDir);

	    if(!fs::exists (_outputDir)) fs::create_directories(_outputDir); 
	    if(!fs::exists (_tmpDir)) fs::create_directories(_tmpDir); 
		
		_readFilename = _tmpDir + "/input.txt";
		Commons::createInputFile(args::get(arg_readFilenames), _readFilename);

		_filename_exe = argv[0];
		
		_fields = new vector<string>();
		_fields_optional = new vector<string>();
		
	}

    void execute (){
		cout << "todo: we can remove tmp files (mapped reads, alignement...)" << endl;
		
		ofstream fileLogs(_tmpDir + "/logs.txt");
		fileLogs.close();


		openLogFile(_tmpDir);
		checkDependencies();

		execute_pipeline();

		cerr << "done!" << endl;

	}

    void execute_pipeline(){


		loadQualityReport(_checkmQualityReportFilename);

		cerr << "Nb bins to process: " << _binToProcess.size() << endl;
		_logFile << "Nb bins to process: " << _binToProcess.size() << endl;

		string command = "";

		for(const string& binName : _binToProcess){

			processBin(binName);
		}

		//command = _filename_exe + " readSelection " + _tmpDir + " " + _tmpDir + "/read_data_init.txt" + " " + _inputFilename + " -t " + to_string(_nbCores);
		//executeCommand(command);



	}

	vector<string> _binToProcess;

	void loadQualityReport(const string& filename){

		cout << filename << endl;

		ifstream inputFile(filename);

		std::string line;
		getline(inputFile, line); //skip header

		while(!inputFile.eof()) {

			line.clear();
			std::getline( inputFile, line);
			std::stringstream buffer(line);
			
			cout << line << endl;

			if(line.empty()) continue;
			
			std::string temp;
			std::vector<string> fields;

			while( getline( buffer, temp, '\t') ) {
				fields.push_back(temp);
			}

			string binName = fields[0];
			float completeness = stof(fields[1]);
			float contamination = stof(fields[2]);

			//cout << binName << " " << completeness << " " << contamination << endl;
			
			if(completeness > 90) {//} && contamination > 5){
				_binToProcess.push_back(binName);
			}
		}

		inputFile.close();
	
	}

	void processBin(const string& binName){

		const string& binFilename = _binDir + "/" + binName + ".fa";
		cerr << "Processing: " << binFilename << endl;
		
		if(!fs::exists(binFilename)){
			cerr << "\tBin filename not found" << endl;
			return;
		}

		string subDir = _tmpDir + "/" + binName + "/";
		if(!fs::exists (subDir)){
			fs::create_directories(subDir); 
		}

		const string& pafFilename = subDir + "/readMapping.paf";
		const string& readFilename = subDir + "/reads.fastq";
		const string& asmDir = subDir + "/asm/";
		const string& evalDir = subDir + "/checkm/";


		cerr << "\tMapping reads" << endl;
		mapReads(binFilename, pafFilename);


		cerr << "\tLoading read mapping" << endl;
		loadReadMapping(pafFilename);

		cerr << "\tNb mapped reads: " << _mappedReadNames.size() << endl;

		cerr << "\tExtracting mapped reads" << endl;
		extractMappedReads(readFilename);

		cerr << "\tReassembling reads" << endl;
		//assemble(readFilename, asmDir);

		//assembleHifiasm(asmDir, readFilename, evalDir);
		//assembleMetaMDBG(asmDir, readFilename, evalDir);
		//exit(1);

		//if(fs::exists(pafFilename)) fs::remove(pafFilename);
	}

	void assembleMetaMDBG(const string& asmDir, const string& readFilename, const string& evalDir){

		//if(fs::exists(asmDir  + "/contigs.fasta.gz")){
		//	cout << "\tcontig file exists, skip" << endl;

		//}
		//else{
			string command = "/pasteur/appa/homes/gbenoit/zeus/tools/metaMDBG_correction/build/bin/metaMDBG asm -t " + to_string(_nbCores) + " " + asmDir + " " + readFilename + " -d 0.02";
			Utils::executeCommand(command, _tmpDir, _logFile);
		//}

		cerr << "\tEvaluating assembly" << endl;
		evaluateQuality(evalDir, asmDir  + "/contigs.fasta.gz");

	}

	void assembleHifiasm(const string& asmDir, const string& readFilename, const string& evalDir){

		const string& contigFilename = asmDir + "/asm.p_ctg.fasta";
		const string& contigFilenameHifiasm = asmDir + "/asm.p_ctg.gfa";

		if (fs::exists(contigFilenameHifiasm)){
			cout << "\tAssembly exists, skip" << endl;
			return;
		}

		string command = "conda run -n hifiasm-meta hifiasm_meta -o " + asmDir + "/asm -t " + to_string(_nbCores) + " " + readFilename;
		Utils::executeCommand(command, _tmpDir, _logFile);
		

		command = "conda run -n biopython python3 ~/zeus/scripts/misc/gfaToFasta.py " + contigFilenameHifiasm + " " + contigFilename;
		Utils::executeCommand(command, _tmpDir, _logFile);

		command = "pigz -f " + contigFilename;
		Utils::executeCommand(command, _tmpDir, _logFile);

		cerr << "\tEvaluating assembly" << endl;
		evaluateQuality(evalDir, contigFilename + ".gz");
	}

	void mapReads(const string& binFilename, const string& pafFilename){

		if(fs::exists(pafFilename)){
			cout << "\talign file exists, skip" << endl;
			return;
		}

		string command = "minimap2 -t " + to_string(_nbCores) + " -x map-ont " + binFilename + " " + Commons::inputFileToFilenames(_readFilename);
		command += " | gzip -c - > " + pafFilename;
		Utils::executeCommand(command, _tmpDir, _logFile);

	}

	void loadReadMapping(const string& pafFilename){

		_mappedReadNames.clear();

		PafParser pafParser(pafFilename);
		auto fp = std::bind(&ReassemblePipeline::parseAlignmentsGz_read, this, std::placeholders::_1);
		pafParser.parse(fp);

	}

	ofstream _readFile;
	unordered_set<string> _mappedReadNames;

	void parseAlignmentsGz_read(const string& line){

		//cout << line << endl;
		GfaParser::tokenize(line, _fields, '\t');

		const string& readName = Utils::shortenHeader((*_fields)[0]);

		double alignLength = stoull((*_fields)[10]);
		double readLength = stoull((*_fields)[1]);

		//if(alignLength < 500) return;
		if(alignLength / readLength < 0.75) return;

		_mappedReadNames.insert(readName);
	}

	void extractMappedReads(const string& readFilename){

		if(fs::exists(readFilename)){
			cout << "\tread file exists, skip" << endl;
			return;
		}

		_readFile = ofstream(readFilename);

		ReadParserParallel readParser(_tmpDir + "/input.txt", false, false, 1, _logFile);
		readParser.parse(ExtractReadFunctor(*this));

		_readFile.close();
	}

	class ExtractReadFunctor {

		public:

		ReassemblePipeline& _parent;

		ExtractReadFunctor(ReassemblePipeline& parent) : _parent(parent){
		}

		ExtractReadFunctor(const ExtractReadFunctor& copy) : _parent(copy._parent){
		}

		~ExtractReadFunctor(){
		}

		void operator () (const Read& read) {

			const string& readName = Utils::shortenHeader(read._header);

			if(_parent._mappedReadNames.find(readName) == _parent._mappedReadNames.end()) return;

			_parent._readFile << "@" << read._header << endl;
			_parent._readFile << read._seq << endl;
			_parent._readFile << "+" << endl;
			_parent._readFile << read._qual << endl;
		}

	};



	void evaluateQuality(const string& resultDir, const string& contigFilename){

		string command = "conda run -n biopython python3 ~/zeus/scripts/paper/run_singleContigs.py " + resultDir + " " + contigFilename + " " + contigFilename + " mdbg " + to_string(_nbCores);
		Utils::executeCommand(command, _tmpDir, _logFile);

	}

	void executeCommand(const string& command){
		Utils::executeCommand(command, _tmpDir, _logFile);
	}

	void checkDependencies(){

		string read = "";
		for(size_t i=0; i<10000; i++){
			read += "A";
		}
		read += "\n";

		string filename = _tmpDir + "/testDependencies.fasta";
		string outputFilename_mapping = _tmpDir + "/testDependencies_align.paf";
		string filenameBzip = _tmpDir + "/testDependencies_bgzip.fasta.gz";
		string outputFilename_mapping_minimap = _tmpDir + "/testDependencies_minimap.paf";

		if(fs::exists(filename)) fs::remove(filename);
		if(fs::exists(filenameBzip)) fs::remove(filenameBzip);
		if(fs::exists(filenameBzip + ".fai")) fs::remove(filenameBzip + ".fai");
		if(fs::exists(outputFilename_mapping)) fs::remove(outputFilename_mapping);

		ofstream file(filename);

		file << ">read1" << endl;
		file << read << endl;
		file << ">read2" << endl;
		file << read << endl;
		
		file.close();


		cerr << "Checking dependencies: minimap2 ";
		_logFile << "Checking dependencies: minimap2 ";
		string command = "minimap2 -x map-hifi " + filename + " " + filename + " > " + outputFilename_mapping_minimap;
		Utils::executeCommand(command, _tmpDir, _logFile);

		if(!fs::exists(outputFilename_mapping_minimap)){
			cerr << endl << "Dependency is missing: minimap2.24+" << endl;
			exit(1);
		}
		else{
			cerr << "ok" << endl;
		}

		cerr << endl;
	}


};	



#endif 


