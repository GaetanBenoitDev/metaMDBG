

#ifndef MDBG_METAG_KMINMERCOUNTER
#define MDBG_METAG_KMINMERCOUNTER

#include "../Commons.hpp"

typedef phmap::parallel_flat_hash_map<KmerVec, vector<u_int32_t>, phmap::priv::hash_default_hash<KmerVec>, phmap::priv::hash_default_eq<KmerVec>, std::allocator<std::pair<KmerVec, vector<u_int32_t>>>, 4, std::mutex> KminmerCountMap;

class KminmerCounter : public Tool{
    
public:

	string _readFilename;
	string _contigFilename;
	string _outputDir;
	int _nbCores;
	float _minimizerDensity;
    size_t _minimizerSize;
    size_t _kminmerSize;

	KminmerCountMap _kminmerCounts;
	size_t _nbDatasets;
	vector<u_int32_t> _countsInit;

	KminmerCounter(): Tool (){

	}

	void parseArgs(int argc, char* argv[]){

		args::ArgumentParser parser("count", "");
		args::Positional<std::string> arg_contigs(parser, "contigs", "", args::Options::Required);
		args::Positional<std::string> arg_outputDir(parser, "outputDir", "", args::Options::Required);
		//args::PositionalList<std::string> arg_readFilenames(parser, "reads", "Read filename(s) (separated by space)", args::Options::Required);
		//args::ValueFlag<int> arg_l(parser, "", "Minimizer length", {ARG_MINIMIZER_LENGTH2}, 13);
		//args::ValueFlag<float> arg_d(parser, "", "Minimizer density", {ARG_MINIMIZER_DENSITY2}, 0.005f);
		args::ValueFlag<int> arg_nbCores(parser, "", "Number of cores", {ARG_NB_CORES2}, NB_CORES_DEFAULT_INT);
		args::Flag arg_help(parser, "", "", {'h', "help"}, args::Options::Hidden);

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

		_contigFilename = args::get(arg_contigs);
		_outputDir = args::get(arg_outputDir);
		_nbCores = args::get(arg_nbCores);

		if(!fs::exists(_outputDir)){
			cerr << "Output dir must be asm dir" << endl;
			exit(1);
			//fs::create_directories(_outputDir);
		}

		_readFilename = _outputDir + "/input.txt";
		_kminmerSize = 4;
		_minimizerSize = 13;
		_minimizerDensity = 0.005;

		
		openLogFile(_outputDir);
	}

    void execute (){
    
		ReadParser parserDatasets(_readFilename, false, false, _logFile);
		_nbDatasets = parserDatasets._nbDatasets;
		cout << "Nb datasets: " << _nbDatasets << endl;

		for(u_int64_t i=0; i<_nbDatasets; i++) _countsInit.push_back(0);

		extractContigKminmers ();

		cout << "Nb contig kminmers: " << _kminmerCounts.size() << endl;

		countKminmers();
		dumpKminmerCounts();

		/*
		cout << "Loading mdbg" << endl;
		string mdbg_filename = _inputDir + "/mdbg_nodes.gz";
		_mdbg = new MDBG(_kminmerSize);
		_mdbg->load(mdbg_filename, false);
		cout << "MDBG nodes: " << _mdbg->_dbg_nodes.size() << endl;

		ReadParser parserDatasets(_inputFilename, false, _minimizerSize, _kminmerSize, _minimizerDensity);
		_nbDatasets = parserDatasets._nbDatasets;
		cout << "Nb datasets: " << _nbDatasets << endl;

		for(u_int64_t i=0; i<_nbDatasets; i++) _countsInit.push_back(0);

		cout << _useContigs << endl;
		if(_useContigs) extractContigKminmers();
		countKminmers();
		if(_useContigs) computeContigCoverage();
		dumpKminmerCount();

		delete _mdbg;
		*/

		closeLogFile();
	}




	void extractContigKminmers (){

		cout << "Extracting contig kminmers" << endl;

		ReadParserParallel readParser(_contigFilename, true, false, _nbCores, _logFile);
		readParser.parse(ExtractContigKminmersFunctor(*this));

	}

	class ExtractContigKminmersFunctor {

		public:

		KminmerCounter& _parent;
		MinimizerParser* _minimizerParser;
		EncoderRLE _encoderRLE;


		ExtractContigKminmersFunctor(KminmerCounter& parent) : _parent(parent){
			_minimizerParser = new MinimizerParser(_parent._minimizerSize, _parent._minimizerDensity);
		}

		ExtractContigKminmersFunctor(const ExtractContigKminmersFunctor& copy) : _parent(copy._parent){
			_minimizerParser = new MinimizerParser(_parent._minimizerSize, _parent._minimizerDensity);
		}

		~ExtractContigKminmersFunctor(){
			delete _minimizerParser;
		}

		void operator () (const Read& read) {

			vector<u_int32_t>& counts = _parent._countsInit;

			string rleSequence;
			vector<u_int64_t> rlePositions;
			_encoderRLE.execute(read._seq.c_str(), read._seq.size(), rleSequence, rlePositions);

			vector<u_int64_t> minimizers;
			vector<u_int64_t> minimizers_pos;
			_minimizerParser->parse(rleSequence, minimizers, minimizers_pos);

			vector<KmerVec> kminmers; 
			vector<ReadKminmer> kminmersInfos;
			MDBG::getKminmers(_parent._minimizerSize, _parent._kminmerSize, minimizers, minimizers_pos, kminmers, kminmersInfos, rlePositions, read._index, false);

			for(size_t i=0; i<kminmersInfos.size(); i++){
				
				KmerVec vec = kminmers[i];
				//for(KmerVec& vec : kminmers){
				//if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()) continue;

				//u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
				//_parent._kminmerCounts[nodeName] = _countsInit;
				_parent._kminmerCounts.lazy_emplace_l(vec, 
				[](KminmerCountMap::value_type& v) { // key exist
				},           
				[&vec, &counts](const KminmerCountMap::constructor& ctor) { // key inserted
					//KminmerSequenceEntire seq = {{}, false};
					ctor(vec, counts); 

				});

			}
		}

	};


    void countKminmers(){

		cout << "Counting kminmers" << endl;

		ReadParserParallel readParser(_readFilename, false, false, _nbCores, _logFile);
		readParser.parse(CountKminmersFunctor(*this));

    }


	class CountKminmersFunctor {

		public:

		KminmerCounter& _parent;
		MinimizerParser* _minimizerParser;
		EncoderRLE _encoderRLE;

		CountKminmersFunctor(KminmerCounter& parent) : _parent(parent){
			_minimizerParser = new MinimizerParser(_parent._minimizerSize, _parent._minimizerDensity);
		}

		CountKminmersFunctor(const CountKminmersFunctor& copy) : _parent(copy._parent){
			_minimizerParser = new MinimizerParser(_parent._minimizerSize, _parent._minimizerDensity);
		}

		~CountKminmersFunctor(){
			delete _minimizerParser;
		}

		void operator () (const Read& read) {

			u_int64_t datasetIndex = read._datasetIndex;

			string rleSequence;
			vector<u_int64_t> rlePositions;
			_encoderRLE.execute(read._seq.c_str(), read._seq.size(), rleSequence, rlePositions);

			vector<u_int64_t> minimizers;
			vector<u_int64_t> minimizers_pos;
			_minimizerParser->parse(rleSequence, minimizers, minimizers_pos);

			vector<KmerVec> kminmers; 
			vector<ReadKminmer> kminmersInfos;
			MDBG::getKminmers(_parent._minimizerSize, _parent._kminmerSize, minimizers, minimizers_pos, kminmers, kminmersInfos, rlePositions, read._index, false);

			for(size_t i=0; i<kminmersInfos.size(); i++){
				

				KmerVec vec = kminmers[i];

				bool lala = _parent._kminmerCounts.modify_if(vec, 
				[&datasetIndex](KminmerCountMap::value_type& v) { // key exist

					v.second[datasetIndex] += 1;
				});

			}
		}

	};

	void dumpKminmerCounts(){

		cout << "Dumping kminmer counts" << endl;

		ofstream outputFile(_outputDir + "/kminmerCounts.bin");

		outputFile.write((const char*)&_minimizerSize, sizeof(_minimizerSize));
		outputFile.write((const char*)&_minimizerDensity, sizeof(_minimizerDensity));
		outputFile.write((const char*)&_kminmerSize, sizeof(_kminmerSize));
		outputFile.write((const char*)&_nbDatasets, sizeof(_nbDatasets));

		for(auto& it : _kminmerCounts){

			vector<u_int64_t> minimizerSeq = it.first._kmers;
			vector<u_int32_t> counts = it.second;

			outputFile.write((const char*)&minimizerSeq[0], _kminmerSize*sizeof(uint64_t));
			outputFile.write((const char*)&counts[0], _nbDatasets*sizeof(uint32_t));

			//for(u_int32_t count : counts){
				//cout << count << " ";
			//}
			//cout << endl;

		}

		outputFile.close();
	}

	/*
	void countKminmers (){

		ReadParser parser(_inputFilename, false, _minimizerSize, _kminmerSize, _minimizerDensity);
		auto fp = std::bind(&KminmerCounter::countKminmers_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6);
		parser.parseKminmers(fp);

	}


	void countKminmers_read(const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex, u_int64_t datasetIndex, const string& header, const string& seq){

		if(readIndex % 100000 == 0){
			cout << datasetIndex << " " << readIndex << endl;
		}
		//cout << datasetIndex << " " << readIndex << endl;
		//cout << "------" << endl;

		//cout << readIndex << endl;
		for(size_t i=0; i<kminmersInfos.size(); i++){
			
			//const ReadKminmer& kminmerInfo = kminmersInfos[i];

			KmerVec vec = kminmers[i];
			//for(KmerVec& vec : kminmers){
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()) continue;

			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;


			if(_kminmerCounts.find(nodeName) == _kminmerCounts.end()){
				if(_useContigs){
					continue;
				}
				else{
					_kminmerCounts[nodeName] = _countsInit;
				}
			}


			_kminmerCounts[nodeName][datasetIndex] += 1;
		}

	}



	void computeContigCoverage (){

		ReadParser parser(_inputFilename_contig, true, _minimizerSize, _kminmerSize, _minimizerDensity);
		auto fp = std::bind(&KminmerCounter::computeContigCoverage_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6);
		parser.parseKminmers(fp);
	}


	void computeContigCoverage_read(const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex, u_int64_t datasetIndex, const string& header, const string& seq){


		vector<vector<u_int32_t>> abundances;
		for(u_int64_t i=0; i<_nbDatasets; i++) abundances.push_back({});
		
		//cout << datasetIndex << " " << readIndex << endl;
		//cout << "------" << endl;

		//cout << readIndex << endl;
		for(size_t i=0; i<kminmersInfos.size(); i++){
			
			//const ReadKminmer& kminmerInfo = kminmersInfos[i];

			KmerVec vec = kminmers[i];
			//for(KmerVec& vec : kminmers){
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()) continue;

			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;

			if(_kminmerCounts.find(nodeName) == _kminmerCounts.end()) continue;

			const vector<u_int32_t>& counts = _kminmerCounts[nodeName];

			for(u_int64_t i=0; i<counts.size(); i++){
				abundances[i].push_back(counts[i]);
			}
		}

		cout << header << ": ";
		for(u_int64_t i=0; i<_nbDatasets; i++){
			//cout << abundances[i].size() << endl;
			float median = Utils::compute_median(abundances[i]);
			cout << median << " ";
		}
		cout << endl;

	}

	void dumpKminmerCount(){

		ofstream outputFile(_outputFilename_kminmerCouts);
		
		string header = "NodeName";
		for(size_t i=0; i<_nbDatasets; i++){
			header += "\tS" + to_string(i);
		}
		outputFile << header << endl;

		for(const auto& it : _kminmerCounts){
			u_int32_t nodeName = it.first;
			const vector<u_int32_t>& counts = it.second;
			
			string line = to_string(nodeName);
			for(u_int32_t count : counts){
				line += "\t" + to_string(count);
			}

			outputFile << line << endl;

		}

		outputFile.close();
	}

	*/

};	


#endif 



