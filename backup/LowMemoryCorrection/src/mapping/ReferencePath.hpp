

#ifndef MDBG_METAG_ReferencePath
#define MDBG_METAG_ReferencePath

#include "../Commons.hpp"

class ReferencePath : public Tool{
    
public:

	string _mdbgDir;
	string _referenceFilename;
	string _outputFilename;

	MinimizerParser* _minimizerParser;
	size_t _minimizerSize;
	size_t _kminmerSize;
	size_t _kminmerSizeFirst;
	float _minimizerDensity;

	MDBG* _mdbg;
	ofstream _outputFile;

	ReferencePath(): Tool (){

	}


	void parseArgs(int argc, char* argv[]){

		/*
		//_kminmerSize = 4;
		//_minimizerSize = 21;
		//_minimizerDensity = 0.05;

		cxxopts::Options options("ToBasespace", "");
		options.add_options()
		("mdbgDir", "", cxxopts::value<string>())
		("binDir", "", cxxopts::value<string>())
		("outputFilename", "", cxxopts::value<string>());

		options.parse_positional({"mdbgDir", "binDir", "outputFilename"});
		options.positional_help("mdbgDir binDir outputFilename");



		if(argc <= 1){
			cout << options.help() << endl;
			exit(0);
		}

		cxxopts::ParseResult result;

		try{
			result = options.parse(argc, argv);

			_mdbgDir = result["mdbgDir"].as<string>();
			_binDir = result["binDir"].as<string>();
			_outputFilename = result["outputFilename"].as<string>();
		}
		catch (const std::exception& e){
			std::cout << options.help() << std::endl;
			std::cerr << e.what() << std::endl;
			std::exit(EXIT_FAILURE);
		}
		*/
		args::ArgumentParser parser("derep", ""); //"This is a test program.", "This goes after the options."
		args::Positional<std::string> arg_mdbgDir(parser, "mdbgDir", "", args::Options::Required);
		args::Positional<std::string> arg_binDir(parser, "referenceFilename", "", args::Options::Required);
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

		_mdbgDir = args::get(arg_mdbgDir);
		_referenceFilename = args::get(arg_binDir);
		_outputFilename = args::get(arg_outputFilename);



		string filename_parameters = _mdbgDir + "/parameters.gz";
		gzFile file_parameters = gzopen(filename_parameters.c_str(),"rb");
		gzread(file_parameters, (char*)&_minimizerSize, sizeof(_minimizerSize));
		gzread(file_parameters, (char*)&_kminmerSize, sizeof(_kminmerSize));
		gzread(file_parameters, (char*)&_minimizerDensity, sizeof(_minimizerDensity));
		gzclose(file_parameters);

		//ifstream file_data(_mdbgDir + "/data.txt");
		//file_data.read((char*)&_nbReads, sizeof(_nbReads));
		//file_data.close();

		cout << endl;
		cout << "Input dir: " << _mdbgDir << endl;
		//cout << "Output filename: " << _outputFilename << endl;
		cout << "Minimizer length: " << _minimizerSize << endl;
		cout << "Kminmer length: " << _kminmerSize << endl;
		cout << "Density: " << _minimizerDensity << endl;
		//cout << "Nb reads: " << _nbReads << endl;
		cout << endl;

		_kminmerSizeFirst = 4;

		_minimizerParser = new MinimizerParser(_minimizerSize, _minimizerDensity);
	}


    void execute (){


		_outputFile = ofstream(_outputFilename);
		_outputFile << "Name,Position" << endl;

		string mdbg_filename = _mdbgDir + "/kminmerData_min.txt";

		//cout << _gfaFilename << endl;
		_mdbg = new MDBG(_kminmerSize);
		_mdbg->load(mdbg_filename, false);
		cout << "Nb nodes: " <<  _mdbg->_dbg_nodes.size() << endl;

		extract_truth_kminmers_bin();
		//mapToMdbg();
		//assingBinToUnitigs();
		//binReads();

		_outputFile.close();
	}

	unordered_map<u_int32_t, vector<u_int32_t>> _nodeName_to_referencePosition;
	//unordered_map<KmerVec, vector<string>> _kmerVec_to_binIndexes;
	//unordered_map<u_int32_t, vector<u_int32_t>> _nodeName_to_unitigIndexes;
	//unordered_map<KmerVec, string> _kmervec_to_binIndex;

	EncoderRLE _encoderRLE;
	//u_int32_t _currentBinIndex;

	/*
	void indexNodeNames(){

		ifstream file(_mdbgDir + "/nodeName_to_unitigIndex.bin");

		while(true){
			
			u_int32_t nodeName;
			u_int32_t unitigIndex;

			file.read((char*)&nodeName, sizeof(nodeName));
			if(file.eof()) break;
			file.read((char*)&unitigIndex, sizeof(unitigIndex));

			_nodeName_to_unitigIndexes[nodeName].push_back(unitigIndex);
		}

		file.close();
	}*/



	void extract_truth_kminmers_bin(){

		auto fp = std::bind(&ReferencePath::extract_truth_kminmers_bin_read, this, std::placeholders::_1);
		ReadParser readParser(_referenceFilename, true, false, _logFile);
		readParser.parse(fp);
	}

	void extract_truth_kminmers_bin_read(const Read& read){
		//ottalSize += strlen(read->seq.s);


		u_int64_t readIndex = read._index;

		//cout << read._seq.size() << endl;
		string rleSequence;
		vector<u_int64_t> rlePositions;
		_encoderRLE.execute(read._seq.c_str(), read._seq.size(), rleSequence, rlePositions);

		vector<u_int64_t> minimizers;
		vector<u_int64_t> minimizers_pos;
		_minimizerParser->parse(rleSequence, minimizers, minimizers_pos);

		vector<KmerVec> kminmers; 
		vector<ReadKminmer> kminmersInfo;
		MDBG::getKminmers(_minimizerSize, _kminmerSize, minimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, readIndex, false);

		cout << read._seq.size() << " " << kminmers.size() << endl;

		for(size_t i=0; i<kminmers.size(); i++){


			//cout << i << endl;
			KmerVec& vec = kminmers[i];
			cout << i << " " << (_mdbg->_dbg_nodes.find(vec) != _mdbg->_dbg_nodes.end()) << endl;
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()) continue;

			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
			
			_outputFile << nodeName << "," << i << endl;
			//_nodeName_to_referencePosition[nodeName] = i;

			//_kmerVec_to_binIndexes[vec].push_back(_currentBinName);
			//cout << _currentBinIndex << endl;
			/*
			cout << "oue?" << endl;
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()) continue;

			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
			cout << nodeName << endl;

			if(_nodeName_to_unitigIndexes.find(nodeName) == _nodeName_to_unitigIndexes.end()) continue;
			*/

			//if(_kmervec_to_unitigName.find(vec) == _kmervec_to_unitigName.end()) continue;

			//string unitigName = _kmervec_to_unitigName[vec];

			//if(_writtenUnitigNames.find(unitigName) != _writtenUnitigNames.end()) continue;

			//_outputFile << unitigName << "," << _currentBinIndex << endl;
			//_writtenUnitigNames.insert(unitigName);
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

