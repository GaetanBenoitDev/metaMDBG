

#ifndef MDBG_METAG_MAPPINGBINTOMDBG
#define MDBG_METAG_MAPPINGBINTOMDBG

#include "../Commons.hpp"

class Mapping_BinToMDBG : public Tool{
    
public:

	string _tmpDir;
	string _mdbgDir;
	string _binDir;
	string _outputFilename;

	MinimizerParser* _minimizerParser;
	size_t _minimizerSize;
	size_t _kminmerSize;
	//size_t _kminmerSizeFirst;
	float _minimizerDensity;

	MDBG* _mdbg;

	Mapping_BinToMDBG(): Tool (){

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
		args::Positional<std::string> arg_tmpDir(parser, "tmpDir", "", args::Options::Required);
		args::Positional<std::string> arg_mdbgDir(parser, "mdbgDir", "", args::Options::Required);
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

		_tmpDir = args::get(arg_tmpDir);
		_mdbgDir = args::get(arg_mdbgDir);
		_binDir = args::get(arg_binDir);
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

		//_kminmerSizeFirst = 4;

		_minimizerParser = new MinimizerParser(_minimizerSize, _minimizerDensity);
	}


    void execute (){

		string mdbg_filename = _mdbgDir + "/kminmerData_min.txt";

		//cout << _gfaFilename << endl;
		_mdbg = new MDBG(_kminmerSize);
		_mdbg->load(mdbg_filename, false);
		cout << "Nb nodes: " <<  _mdbg->_dbg_nodes.size() << endl;

		indexNodeNames();
		extract_truth_kminmers();
		mapToMdbg();
		assingBinToUnitigs();
		//binReads();
	}



	unordered_map<KmerVec, vector<string>> _kmerVec_to_binIndexes;
	unordered_map<u_int32_t, vector<u_int32_t>> _nodeName_to_unitigIndexes;
	//unordered_map<KmerVec, string> _kmervec_to_binIndex;

	EncoderRLE _encoderRLE;
	//u_int32_t _currentBinIndex;

	void indexNodeNames(){

		ifstream file(_mdbgDir + "/nodeName_to_unitigIndex.bin");

		while(true){
			
			u_int32_t nodeName;
			u_int32_t unitigIndex;

			file.read((char*)&nodeName, sizeof(nodeName));
			if(file.eof()) break;
			file.read((char*)&unitigIndex, sizeof(unitigIndex));

			//cout << unitigIndex << endl;
			_nodeName_to_unitigIndexes[nodeName].push_back(unitigIndex);

			//if(unitigIndex == 4052){
			//	cout << nodeName << endl;
			//}
		}

		file.close();
	}

	string _currentBinName;

	void extract_truth_kminmers(){

		//_currentBinIndex = 0;
		//_outputFile = ofstream(_outputFilename);
		//_outputFile << "Name,Color" << endl;

		for (const auto & p : fs::directory_iterator(_binDir)){
			string ext = p.path().extension();
			//cout << p.path() << endl;

			//cout << ext << endl;
			if(ext == ".fa" || ext == ".fasta" || ext == ".fna"){
				string filename = p.path();
				cout << filename << endl;

				_currentBinName = p.path().filename();
				//string binName = p.path().filename();
				//binName.erase(binName.find("bin."), 4);
				//binName.erase(binName.find(ext), ext.size());
				//cout << binName << endl;
				//_currentBinIndex = stoull(binName);

				_nbContigs = 0;
				countBinNbContigs(filename);

				if(_nbContigs > 50) continue;

				extract_truth_kminmers_bin(filename);

				//_currentBinIndex += 1;
			}

		}

		//_outputFile.close();
	}

	size_t _nbContigs;
	
	void countBinNbContigs(const string& binFilename){

		auto fp = std::bind(&Mapping_BinToMDBG::countBinNbContigs_read, this, std::placeholders::_1);
		ReadParser readParser(binFilename, true, false, _logFile);
		readParser.parse(fp);
	}

	void countBinNbContigs_read(const Read& read){
		_nbContigs += 1;
	}


	void extract_truth_kminmers_bin(const string& binFilename){

		auto fp = std::bind(&Mapping_BinToMDBG::extract_truth_kminmers_bin_read, this, std::placeholders::_1);
		ReadParser readParser(binFilename, true, false, _logFile);
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

		_contigKminmers.clear();
		_higherNbKminmers = 0;

		for(size_t i=0; i<kminmers.size(); i++){
			KmerVec& vec = kminmers[i];
			_contigKminmers.insert(vec);
		}

		extractContigSequence();
		//cout <<"seq:" <<  _originalContigKminmers.size() << endl;
		for(size_t i=0; i<_originalContigKminmers.size(); i++){


			//cout << i << endl;
			KmerVec& vec = _originalContigKminmers[i];
			_kmerVec_to_binIndexes[vec].push_back(_currentBinName);

		}


	}

	int _higherNbKminmers;
	unordered_set<KmerVec> _contigKminmers;
	vector<KmerVec> _originalContigKminmers;

	void extractContigSequence(){

		const string& filename_contigs = _tmpDir + "/contig_data.txt";
		KminmerParser parser(filename_contigs, -1, _kminmerSize, false, false);
		auto fp = std::bind(&Mapping_BinToMDBG::extractContigSequence_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
		parser.parse(fp);

	}

	void extractContigSequence_read(const vector<u_int64_t>& readMinimizers, const vector<KmerVec>& vecs, const vector<ReadKminmer>& kminmersInfos, u_int8_t isCircular, u_int32_t readIndex){
		

		//cout << readMinimizers.size() << endl;
		
		int nbKminmers = 0;

		for(u_int32_t i=0; i<vecs.size(); i++){
			
			KmerVec vec = vecs[i];

			if(_contigKminmers.find(vec) != _contigKminmers.end()){
				nbKminmers += 1;
			}
		}

		//cout << readMinimizers.size() << " " << nbKminmers << endl;

		if(nbKminmers > _higherNbKminmers){
			_higherNbKminmers = nbKminmers;
			_originalContigKminmers = vecs;
		}
	}

	void mapToMdbg(){


		for(const auto& it : _mdbg->_dbg_nodes){

			const KmerVec& vec = it.first;
			u_int32_t nodeName = it.second._index;

			//cout << nodeName << endl;
			//if(_nodeName_to_unitigIndexes.find(nodeName) == _nodeName_to_unitigIndexes.end()) continue;

			const vector<string>& binIndexes = _kmerVec_to_binIndexes[vec];

			//cout << nodeName << " " << binIndexes.size() << endl;

			for(u_int32_t unitigIndex : _nodeName_to_unitigIndexes[nodeName]){
				for(const string& binName : binIndexes){
					_unitigIndex_to_binIndexes[unitigIndex].push_back(binName);
					//cout << unitigIndex << " " << binName << endl;
				}
			}

			/*
			vector<u_int64_t> minimizers = vec._kmers;
			//vector<u_int64_t> minimizers = vec._kmers;
			//if(it.second._isReversed){
			//	std::reverse(minimizers.begin(), minimizers.end());
			//}
			vector<u_int64_t> rlePositions;
			vector<u_int64_t> minimizers_pos;//(minimizers.size());
			vector<KmerVec> kminmers; 
			vector<ReadKminmer> kminmersInfo;
			MDBG::getKminmers(_minimizerSize, _kminmerSize, minimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, 0, false);

			for(const KmerVec& vec : kminmers){

				//cout << "?" << endl;
				if(_kmerVec_to_binIndexes.find(vec) == _kmerVec_to_binIndexes.end()) continue;

				const vector<string>& binIndexes = _kmerVec_to_binIndexes[vec];

				for(u_int32_t unitigIndex : _nodeName_to_unitigIndexes[nodeName]){
					for(const string& binName : binIndexes){
						_unitigIndex_to_binIndexes[unitigIndex].push_back(binName);
					}
				}
			}
			*/
		
		}
	}
	
	unordered_map<u_int32_t, vector<string>> _unitigIndex_to_binIndexes;
	//unordered_map<KmerVec, vector<u_int32_t>> _kmerVec_to_binIndexes;
	//unordered_map<u_int32_t, vector<u_int32_t>> _nodeName_to_unitigIndexes;

	void assingBinToUnitigs(){

		ofstream file(_outputFilename);
		file << "Name,Color" << endl;

		for(const auto& it : _unitigIndex_to_binIndexes){
			u_int32_t unitigIndex = it.first;
			const vector<string>& binIndexes = it.second;

			unordered_map<string, u_int32_t> binCounts;

			for(const string& binIndex : binIndexes){
				binCounts[binIndex] += 1;
			}

			string maxBinIdnex = "";
			u_int32_t maxBinCount = 0;
			for(const auto& it : binCounts){
				const string& binIndex = it.first;
				u_int32_t binCount = it.second;

				if(binCount > maxBinCount){
					maxBinCount = binCount;
					maxBinIdnex = binIndex;
				}
			}

			file << unitigIndex << "," << maxBinIdnex << endl;
		}

		file.close();
	}
};	


#endif 


