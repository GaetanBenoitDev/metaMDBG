

#ifndef MDBG_METAG_CONTIGLINKER
#define MDBG_METAG_CONTIGLINKER


#include "Commons.hpp"
#include "graph/GfaParser.hpp"
#include "graph/GraphSimplify.hpp"










class Circulizer2 : public Tool{

public:

	string _inputDir;
	string _passDir;
	string _truthInputFilename;
	string _outputFilename;
	string _outputFilename_complete;
	bool _debug;
	bool _isFinalAssembly;
	string _inputFilename_unitigNt;
	string _inputFilename_unitigCluster;
	string _filename_abundance;
	MDBG* _mdbg;

	float _minimizerDensity;
    size_t _minimizerSize;
    size_t _kminmerSize;
    size_t _kminmerSizeFirst;
	float _minimizerSpacingMean;
	float _kminmerLengthMean;
	float _kminmerOverlapMean;
	vector<UnitigData> _unitigDatas;

	string _partitionDir;
	string _filename_solidKminmers;
	string _filename_inputContigs;
	string _filename_outputContigs;
	string _filename_readMinimizers;
	string _filename_hifiasmGroundtruth;
	string _filename_contigs;
	string _referenceFilename;
	unordered_map<KmerVec, u_int16_t> _evaluation_hifiasmGroundTruth;
	unordered_map<KmerVec, u_int32_t> _evaluation_hifiasmGroundTruth_position;
	unordered_map<u_int32_t, u_int32_t> _evaluation_hifiasmGroundTruth_nodeNamePosition;
	gzFile _outputContigFile;
	gzFile _outputContigFile_complete;
	GraphSimplify* _graph;
	//ToBasespaceOnTheFly _toBasespace;
	//ContigFeature _contigFeature;
	string _gfaFilename;
	size_t _nbCores;

	u_int64_t _nbPathSolved;
	unordered_map<KmerVec, vector<u_int64_t>> _kmerVec_to_readIndexes;
	unordered_set<u_int32_t> _usedNodeNames;


	struct StartingContig{
		u_int64_t _index;
		u_int64_t _length;
		vector<KmerVec> _startKmerVecs;
		vector<KmerVec> _endKmerVecs;
		bool _isCircular;
	};

	unordered_map<u_int64_t, StartingContig> _startingContigs;

	Circulizer2(): Tool (){

	}

	~Circulizer2(){

	}


	void parseArgs(int argc, char* argv[]){

		args::ArgumentParser parser("circ", ""); //"This is a test program.", "This goes after the options."
		args::Positional<std::string> arg_outputDir(parser, "outputDir", "", args::Options::Required);
		//args::Positional<std::string> arg_contigs(parser, "contigs", "", args::Options::Required);
		//args::PositionalList<std::string> arg_readFilenames(parser, "reads", "Input filename(s) (separated by space)", args::Options::Required);
		//args::ValueFlag<int> arg_l(parser, "", "Minimizer length", {ARG_MINIMIZER_LENGTH2}, 13);
		//args::ValueFlag<float> arg_d(parser, "", "Minimizer density", {ARG_MINIMIZER_DENSITY2}, 0.005f);
		//args::ValueFlag<std::string> arg_contigs(parser, "", "", {ARG_INPUT_FILENAME_CONTIG}, "");
		args::ValueFlag<string> arg_referenceFilename(parser, "", "reference filename", {"ref"}, "");
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

		_inputDir = args::get(arg_outputDir);
		_passDir = _inputDir + "/pass_k4/";
		//_passDir = _inputDir + "/pass_k4/";
		//_filename_inputContigs = args::get(arg_contigs);

		_nbCores = args::get(arg_nbCores);

		_filename_contigs = _inputDir + "/contig_data.txt";

		_referenceFilename = "";
		if(arg_referenceFilename){
			_referenceFilename = args::get(arg_referenceFilename);
		}

		string filename_parameters = _passDir + "/parameters.gz";
		gzFile file_parameters = gzopen(filename_parameters.c_str(),"rb");
		gzread(file_parameters, (char*)&_minimizerSize, sizeof(_minimizerSize));
		gzread(file_parameters, (char*)&_kminmerSize, sizeof(_kminmerSize));
		gzread(file_parameters, (char*)&_minimizerDensity, sizeof(_minimizerDensity));
		gzread(file_parameters, (char*)&_kminmerSizeFirst, sizeof(_kminmerSizeFirst));
		gzread(file_parameters, (char*)&_minimizerSpacingMean, sizeof(_minimizerSpacingMean));
		gzread(file_parameters, (char*)&_kminmerLengthMean, sizeof(_kminmerLengthMean));
		gzread(file_parameters, (char*)&_kminmerOverlapMean, sizeof(_kminmerOverlapMean));
		gzclose(file_parameters);

		openLogFile(_inputDir);

		//_kminmerSize = _kminmerSizeFirst;

		_logFile << endl;
		_logFile << "Input dir: " << _inputDir << endl;
		//cout << "Output filename: " << _outputFilename << endl;
		_logFile << "Minimizer length: " << _minimizerSize << endl;
		_logFile << "Kminmer length: " << _kminmerSize << endl;
		_logFile << "Density: " << _minimizerDensity << endl;
		_logFile << endl;

		_outputFilename = _inputDir + "/minimizer_contigs.gz";
		_outputFilename_complete = _inputDir + "/minimizer_contigs_complete.gz";
		_filename_readMinimizers = _inputDir + "/read_data_init.txt";
		_filename_hifiasmGroundtruth = _inputDir + "/hifiasmGroundtruth.gz";
		_filename_outputContigs = _inputDir + "/contigs.min.gz";
		_filename_solidKminmers = _inputDir + "/solid.min.gz";

	}

	void execute (){

		_nbPathSolved = 0;

		collectStartingContigs();
		indexReads();
		circularize();
		dumpContigs();

		_logFile << "Nb path solved: " << _nbPathSolved << endl;
		closeLogFile();
	}



	void collectStartingContigs(){
		
		
		_logFile << "load starting contigs " << endl;
		

		KminmerParserParallel parser(_filename_contigs, _minimizerSize, _kminmerSize, false, false, 1);
		parser.parse(CollectStartingContigsFunctor(*this));


		_logFile << "Nb candidates: " << _startingContigs.size() << endl;
	} 


	class CollectStartingContigsFunctor {

		public:

		Circulizer2& _graph;

		CollectStartingContigsFunctor(Circulizer2& graph) : _graph(graph){

		}

		CollectStartingContigsFunctor(const CollectStartingContigsFunctor& copy) : _graph(copy._graph){

		}

		~CollectStartingContigsFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {

			u_int64_t readIndex = kminmerList._readIndex;
			const vector<u_int64_t>& readMinimizers = kminmerList._readMinimizers;
			const vector<ReadKminmerComplete>& kminmersInfos = kminmerList._kminmersInfo;

			//cout << kminmersInfos.size() << endl;
			if(kminmerList._isCircular) return;
			if(kminmersInfos.size() < 5000) return;

			vector<KmerVec> startNodeNames;
			vector<KmerVec> endNodeNames;

			for(size_t i=0; i<20; i++){

				const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];
				const KmerVec& vec = kminmerInfo._vec;

				startNodeNames.push_back(vec);
				_graph._kmerVec_to_readIndexes[vec] = {};
			}

			for(size_t i=kminmersInfos.size()-20; i<kminmersInfos.size(); i++){

				const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];
				const KmerVec& vec = kminmerInfo._vec;

				endNodeNames.push_back(vec);
				_graph._kmerVec_to_readIndexes[vec] = {};
			}

			_graph._startingContigs[readIndex] = {readIndex, kminmersInfos.size()*_graph._minimizerSpacingMean, startNodeNames, endNodeNames, false};
		}
	};



	void indexReads(){

		_logFile << "Indexing reads" << endl;

		KminmerParserParallel parser(_filename_readMinimizers, _minimizerSize, _kminmerSize, false, true, _nbCores);
		parser.parse(IndexReadsFunctor(*this));

	}

	class IndexReadsFunctor {

		public:

		Circulizer2& _parent;
		
		IndexReadsFunctor(Circulizer2& parent) : _parent(parent){
		}

		IndexReadsFunctor(const IndexReadsFunctor& copy) : _parent(copy._parent){
		}

		~IndexReadsFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {

			u_int32_t readIndex = kminmerList._readIndex;

			for(u_int16_t i=0; i<kminmerList._kminmersInfo.size(); i++){
			
				const ReadKminmerComplete& kminmerInfo = kminmerList._kminmersInfo[i];

				KmerVec vec = kminmerInfo._vec;
				if(_parent._kmerVec_to_readIndexes.find(vec) == _parent._kmerVec_to_readIndexes.end()) continue;

				#pragma omp critical
				{
					_parent._kmerVec_to_readIndexes[vec].push_back(readIndex);
				}

			}

		}
	};


	void circularize(){
		
		
		_logFile << "Circularizing " << endl;
		

		KminmerParserParallel parser(_filename_contigs, _minimizerSize, _kminmerSize, false, false, 1);
		parser.parse(CircularizeFunctor(*this));

	}


	class CircularizeFunctor {

		public:

		Circulizer2& _graph;

		CircularizeFunctor(Circulizer2& graph) : _graph(graph){

		}

		CircularizeFunctor(const CircularizeFunctor& copy) : _graph(copy._graph){

		}

		~CircularizeFunctor(){
		}


		void operator () (const KminmerList& kminmerList) {

			//cout << "---------------" << endl;

			u_int64_t readIndex = kminmerList._readIndex;
			const vector<u_int64_t>& readMinimizers = kminmerList._readMinimizers;
			const vector<ReadKminmerComplete>& kminmersInfos = kminmerList._kminmersInfo;

			if(_graph._startingContigs.find(readIndex) == _graph._startingContigs.end()) return;

			unordered_set<KmerVec> startKmerVecs;
			unordered_set<u_int64_t> startReadIndexes;
			unordered_set<u_int64_t> sharedReadIndex;

			unordered_map<u_int64_t, u_int32_t> startReadMapCount;
			unordered_map<u_int64_t, u_int32_t> endReadMapCount;

			for(size_t i=0; i<20; i++){

				const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];
				const KmerVec& vec = kminmerInfo._vec;

				if(_graph._kmerVec_to_readIndexes.find(vec) == _graph._kmerVec_to_readIndexes.end()) continue;

				for(u_int64_t readIndex : _graph._kmerVec_to_readIndexes[vec]){
					startReadIndexes.insert(readIndex);
					startReadMapCount[readIndex] += 1;
				}

				startKmerVecs.insert(vec);
			}

			for(size_t i=kminmersInfos.size()-20; i<kminmersInfos.size(); i++){

				const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];
				const KmerVec& vec = kminmerInfo._vec;

				if(startKmerVecs.find(vec) != startKmerVecs.end()) continue;
				if(_graph._kmerVec_to_readIndexes.find(vec) == _graph._kmerVec_to_readIndexes.end()) continue;

				for(u_int64_t readIndex : _graph._kmerVec_to_readIndexes[vec]){
					endReadMapCount[readIndex] += 1;

					if(sharedReadIndex.find(readIndex) != sharedReadIndex.end()) continue;
					if(startReadIndexes.find(readIndex) == startReadIndexes.end()) continue;

					sharedReadIndex.insert(readIndex);
					//nbSharedReads += 1;

				}
			}

			int nbSharedReads = 0;
			for(u_int64_t readIndex : sharedReadIndex){
				if(startReadMapCount[readIndex] > 4 && endReadMapCount[readIndex] > 4){
					nbSharedReads += 1;

					//cout << "Overlappin read: " << readIndex << endl;
				}
				//cout << "\t" << startReadMapCount[readIndex] << " " << endReadMapCount[readIndex] << endl;
			}


			if(nbSharedReads > 1){
				#pragma omp critical
				{

					_graph._nbPathSolved += 1;
					_graph._logFile << "Nb shared reads: " << _graph._startingContigs[readIndex]._length << " " << nbSharedReads << endl;
					for(u_int64_t readIndex : sharedReadIndex){
						//cout << "\t" << startReadMapCount[readIndex] << " " << endReadMapCount[readIndex] << endl;
					}

					_graph._startingContigs[readIndex]._isCircular = true;
				}
			}
		}
	};


	ofstream _outputContigFileTmp;

	void dumpContigs(){

		_logFile << "Dumping contigs" << endl;

		string contigFilenameTmp = _inputDir + "/contig_data_derep.txt";
		_outputContigFileTmp = ofstream(contigFilenameTmp);

		KminmerParser parser(_filename_contigs, _minimizerSize, _kminmerSize, false, false);
		auto fp = std::bind(&Circulizer2::dumpDereplicatedContigs_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
		parser.parseSequences(fp);

		_outputContigFileTmp.close();

		fs::remove(_filename_contigs);
		fs::rename(contigFilenameTmp, _filename_contigs);

	}

	void dumpDereplicatedContigs_read(const vector<u_int64_t>& readMinimizers, u_int8_t isCircular, u_int64_t readIndex){

		u_int32_t contigSize = readMinimizers.size();
		_outputContigFileTmp.write((const char*)&contigSize, sizeof(contigSize));

		if(_startingContigs.find(readIndex) == _startingContigs.end()){
			if(isCircular){
				//cout << "1" << endl;
				_outputContigFileTmp.write((const char*)&CONTIG_CIRCULAR, sizeof(isCircular));
			}
			else{
				//cout << "2" << endl;
				_outputContigFileTmp.write((const char*)&CONTIG_LINEAR, sizeof(isCircular));
			}
		}
		else{
			isCircular = _startingContigs[readIndex]._isCircular;
			if(isCircular){
				//cout << "3" << endl;
				_outputContigFileTmp.write((const char*)&CONTIG_CIRCULAR_RESCUED, sizeof(isCircular));
			}
			else{
				//cout << "4" << endl;
				_outputContigFileTmp.write((const char*)&CONTIG_LINEAR, sizeof(isCircular));
			}
			//_outputContigFileTmp.write((const char*)&isCircular, sizeof(isCircular));
		}
		_outputContigFileTmp.write((const char*)&readMinimizers[0], contigSize*sizeof(u_int64_t));

	}
};





#endif
