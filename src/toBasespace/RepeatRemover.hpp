



#ifndef MDBG_METAG_RepeatRemover
#define MDBG_METAG_RepeatRemover

#include "../Commons.hpp"
#include "ReadVsContigMapper.hpp"

class RepeatRemover : public Tool{
  

public:

	u_int32_t _debugFindFinalContigIndex = -1;
	unordered_set<MinimizerType> _debugFindFinalContigMinimizers;
	//typedef phmap::parallel_flat_hash_map<u_int128_t, vector<ReadType>, phmap::priv::hash_default_hash<u_int128_t>, phmap::priv::hash_default_eq<u_int128_t>, std::allocator<std::pair<u_int128_t, vector<ReadType>>>, 4, std::mutex> KminmerPosMap;
	//typedef phmap::parallel_flat_hash_map<u_int128_t, u_int32_t, phmap::priv::hash_default_hash<u_int128_t>, phmap::priv::hash_default_eq<u_int128_t>, std::allocator<std::pair<u_int128_t, u_int32_t>>, 4, std::mutex> KminmerMap;
	//typedef phmap::parallel_flat_hash_map<ReadType, u_int32_t, phmap::priv::hash_default_hash<ReadType>, phmap::priv::hash_default_eq<ReadType>, std::allocator<std::pair<ReadType, u_int32_t>>, 4, std::mutex> ReadLengthMap;
	
	string _inputFilenameContig;
	string _outputFilenameContig;
	string _inputDir;
	int _nbCores;

	Parameters _params;
	
	/*
	float _minimizerDensity;
    size_t _minimizerSize;
    size_t _kminmerSize;
	size_t _kminmerSizeFirst;

	float _minimizerSpacingMean;
	float _kminmerLengthMean;
	float _kminmerOverlapMean;
    size_t _kminmerSizePrev;
    size_t _kminmerSizeLast;
    size_t _meanReadLength;
	*/

	
	ankerl::unordered_dense::map<u_int128_t, u_int32_t> _kminmer_to_unitigIndex;
	ankerl::unordered_dense::map<u_int128_t, u_int32_t> _kminmer_to_abundance;
	//KminmerMap _kminmer_to_unitigIndex;
	//KminmerMap _kminmer_to_abundance;
	//KminmerPosMap _kminmer_to_readIndexes;
	//ofstream _outputFile;
	//unordered_map<ReadType, u_int32_t> _readIndex_to_readLength;

	struct ReadAbundance{
		u_int32_t _readLength;
		float _readAbundance;
	};

	struct ContigsWriter{
		u_int32_t _readIndex;
		vector<vector<MinimizerType>> _contigs;
		u_int8_t _isCircular;
	};

    struct ContigsWriter_Comparator {
        bool operator()(ContigsWriter const& p1, ContigsWriter const& p2){
            return p1._readIndex > p2._readIndex;
        }
    };

	struct FragmentCoverageDisk{
		u_int32_t _startPos;
		u_int32_t _endPos;
		float _coverage;
	};

	struct FragmentCoverageDiskWriter{
		u_int32_t _readIndex;
		vector<MinimizerType> _minimizers;
		u_int8_t _isCircular;
		vector<FragmentCoverageDisk> _fragments;
	};
	
	struct FragmentCoverageDiskWriter_Comparator {
		bool operator()(FragmentCoverageDiskWriter const& p1, FragmentCoverageDiskWriter const& p2){
			return p1._readIndex > p2._readIndex;
		}
	};

	struct Fragment{
		u_int32_t _fragmentIndex;
		u_int32_t _startPos;
		u_int32_t _endPos;
		u_int32_t _length;
		long double _coverage;
		int32_t _finalContigIndex;
		vector<ReadType> _readIndexes;
		unordered_map<u_int32_t, u_int32_t> _nbBridgingReads;
	};

	struct FragmentPath{
		u_int32_t _fragmentIndexStart;
		u_int32_t _fragmentIndexEnd;

		u_int32_t getPathLength() const{
			return _fragmentIndexEnd - _fragmentIndexStart;
		}
	};

	u_int64_t _checksum;

	//ankerl::unordered_dense::set<u_int32_t> _isContigIndexed;

	/*
	RepeatRemover(const string& inputDir, const string& inputFilenameContig, const string& outputFilenameContig, size_t kminmerSize, float minimizerDensity, bool hasQuality, int nbCores, float minimizerDensity_assembly){
		_inputDir = inputDir;
		_inputFilenameContig = inputFilenameContig;
		_outputFilenameContig = outputFilenameContig;
		_kminmerSize = kminmerSize;
		_minimizerDensity = minimizerDensity;
		_hasQuality = hasQuality;
		_nbCores = nbCores;
		
		_averageDistanceBetweenMinimizers = 1.0 / minimizerDensity_assembly;
	}
	*/

	RepeatRemover(): Tool (){


	}
	

	void parseArgs(int argc, char* argv[]){

		args::ArgumentParser parser("removeRepeats", "");
		args::Positional<std::string> arg_outputDir(parser, "outputDir", "", args::Options::Required);
		//args::Positional<std::string> arg_inputContigFilename(parser, "inputContigFilename", "", args::Options::Required);
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

		_inputDir = args::get(arg_outputDir);
		_nbCores = args::get(arg_nbCores);
		
		_inputFilenameContig = _inputDir + "/contig_data_init_small.txt.nooverlaps";
		_outputFilenameContig = _inputDir + "/contig_data_init_small.txt.norepeats";

		_params.load(_inputDir + "/parameters.gz");
		/*
		string filename_parameters = _inputDir + "/parameters.gz";
		gzFile file_parameters = gzopen(filename_parameters.c_str(),"rb");
		gzread(file_parameters, (char*)&_minimizerSize, sizeof(_minimizerSize));
		gzread(file_parameters, (char*)&_kminmerSize, sizeof(_kminmerSize));
		gzread(file_parameters, (char*)&_minimizerDensity, sizeof(_minimizerDensity));
		gzread(file_parameters, (char*)&_kminmerSizeFirst, sizeof(_kminmerSizeFirst));
		gzread(file_parameters, (char*)&_minimizerSpacingMean, sizeof(_minimizerSpacingMean));
		gzread(file_parameters, (char*)&_kminmerLengthMean, sizeof(_kminmerLengthMean));
		gzread(file_parameters, (char*)&_kminmerOverlapMean, sizeof(_kminmerOverlapMean));
		gzread(file_parameters, (char*)&_kminmerSizePrev, sizeof(_kminmerSizePrev));
		gzread(file_parameters, (char*)&_kminmerSizeLast, sizeof(_kminmerSizeLast));
		gzread(file_parameters, (char*)&_meanReadLength, sizeof(_meanReadLength));
		gzclose(file_parameters);
		*/

		_params._kminmerSize = _params._kminmerSizeFirst+1;

		openLogFile(_inputDir);

		/*
		Logger::get().debug() << "";
		Logger::get().debug() << "Input dir: " << _inputDir;
		//_logFile << "Output filename: " << _outputFilename << endl;
		Logger::get().debug() << "Minimizer length: " << _minimizerSize;
		Logger::get().debug() << "Kminmer length: " << _kminmerSize;
		Logger::get().debug() << "Density: " << _minimizerDensity;
		Logger::get().debug() << "";
		*/
	}


	class DebugFindContigFunctor {

		public:

		RepeatRemover& _parent;

		DebugFindContigFunctor(RepeatRemover& parent) : _parent(parent){
		}

		DebugFindContigFunctor(const DebugFindContigFunctor& copy) : _parent(copy._parent){
		}

		~DebugFindContigFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {
			if(kminmerList._readIndex == _parent._debugFindFinalContigIndex){
				for(const auto& m : kminmerList._readMinimizers){
					_parent._debugFindFinalContigMinimizers.insert(m);
				}
				cout << "Found debug contig: " << _parent._debugFindFinalContigMinimizers.size() << endl;
			}
		}
	};


	void execute(){

		if(_debugFindFinalContigIndex != -1){
			KminmerParserParallel parser(_inputDir + "/contig_data_final.bin", 0, 0, false, false, _nbCores);
			parser.parseSequences(DebugFindContigFunctor(*this));
		}

		alignReads();
		fragmentContigs();
		computeFragmentsCoverage();
		breakUnbridgedRepeats();

		Logger::get().debug() << "RepeatRemover checksum: " << _checksum; 

		fs::remove(_inputDir + "/contig_data_init_small.txt.nooverlaps.fragments");
		fs::remove(_inputDir + "/contig_data_init_small.txt.nooverlaps.fragments.coverage");
		fs::remove(_inputDir + "/readsVsContigsAlignments.bin");
		/*
		_nextReadIndexWriter = 0;
		_checksum = 0;
		_outputFile.open(_outputFilenameContig);

		processContigs();
		
		_outputFile.close();

		cout << "RepeatRemover checksum: " << _checksum << endl; 
		cout << _readWriterQueue.size() << endl;
		//cout << "Sorting read index" << endl;
		//for(auto& it : _kminmer_to_readIndexes){
		//	std::sort(it.second.begin(), it.second.end());
		//} 
		




		//Logger::get().debug() << "Break unbridged repeats";
		//breakUnbridgedRepeats();


		//fs::remove(_inputFilenameContig);
		//fs::rename(_inputFilenameContig + ".norepeats", _inputFilenameContig);

		//KminmerParser parser(_inputFilenameContig, -1, _kminmerSizeRepeat, false, false);
		//auto fp = std::bind(&OverlapRemover::detectWeakRepeats_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
		//parser.parse(fp);

		//_kminmerToIndex.clear();
		//_kminmerAbundances.clear();
		//_kminmer_to_unitigIndex.clear();
		//_kminmer_to_readIndexes.clear();
		//cout << "Done" << endl;
		*/
	}
	
	void alignReads(){

		Logger::get().debug() << "\tAlign reads vs contigs";
		auto start = high_resolution_clock::now();


		ReadVsContigMapper mapper(_inputDir + "/read_data_init.txt", _inputFilenameContig, _inputDir + "/readsVsContigsAlignments.bin", 1.0 / _params._minimizerDensity_assembly, _nbCores);
		mapper.execute();

		Logger::get().debug() << "\t\tDone: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";

	}

	ankerl::unordered_dense::map<u_int32_t, vector<pair<u_int32_t, u_int32_t>>> _contigIndex_to_alignments;

	void loadAlignments(){

		//_contigIndex_to_alignments.clear();
		/*
		for(const auto& it : _readIndex_to_bestContigMap){
			const ReadMapping2& al = it.second;
			_contigIndex_to_alignments[al._contigIndex].push_back(al);
		}

		_readIndex_to_bestContigMap.clear();
		*/

		ifstream alignmentFile(_inputDir + "/readsVsContigsAlignments.bin");

		while(true){

			ReadMapping2 alignment;
			bool isEof = alignment.read(alignmentFile);

			if(isEof) break;

			//if(_isContigIndexed.find(alignment._contigIndex) == _isContigIndexed.end()) continue; //Contig not in current chunk
			//if(_mReads.find(alignment._readIndex) == _mReads.end()) continue;

			//if(_contigIndex_to_alignments.find(alignment._contigIndex) == _contigIndex_to_alignments.end()){
			//	_contigIndex_to_alignments[alignment._contigIndex] = ContigData();
			//}

			_contigIndex_to_alignments[alignment._contigIndex].push_back({alignment._contigStart, alignment._contigEnd});
			//cout << "load al: " << alignment._contigIndex << endl;
			//_contigIndex_to_alignments[alignment._contigIndex]._alignments.push_back(alignment);

		}

		alignmentFile.close();

	}




	/*
	u_int64_t _nbContigs;


	void processContigs(){

		_nbContigs = 0;

		//u_int64_t maxBps = 1000000000;
		//u_int64_t maxBps = 3000000000;
		u_int64_t maxNbMinimizers = 15000000; //maxBps / (1/_minimizerDensity);
		//cout << "Max nb minimizers: " << maxNbMinimizers << endl;

		KminmerParserParallel parser(_inputFilenameContig, -1, _kminmerSize, false, false, 1);
		parser._maxChunkSize = maxNbMinimizers;
		parser.parseChunk(ContigFunctor(*this), ChunkFunctor(*this));

		//cout << "Nb contigs: " << _nbContigs << endl;
	}


	vector<KminmerList> _contigs;


	class ContigFunctor {

		public:

		RepeatRemover& _parent;

		ContigFunctor(RepeatRemover& parent) : _parent(parent){
		}

		ContigFunctor(const ContigFunctor& copy) : _parent(copy._parent){
		}

		~ContigFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {
			//cout << "Load contig: " << kminmerList._readIndex << endl;
			_parent._contigs.push_back(kminmerList);
		}
	};

	class ChunkFunctor {

		public:

		RepeatRemover& _parent;

		ChunkFunctor(RepeatRemover& parent) : _parent(parent){
		}

		ChunkFunctor(const ChunkFunctor& copy) : _parent(copy._parent){
		}

		~ChunkFunctor(){
		}

		void operator () () const {
			_parent.processChunk();
			_parent._contigs.clear();
		}

	};

	void processChunk(){


		//cout << "Process chunk: " << _contigs.size() << endl;

		//_nbContigs += _contigs.size();

		//_kminmer_to_unitigIndex.clear();
		//_kminmer_to_abundance.clear();
		//_kminmer_to_readIndexes.clear();
		//_readIndex_to_readLength.clear();
		
		vector<MinimizerPosition>().swap(_allMinimizerPositions); //_allMinimizerPositions.clear();
		ankerl::unordered_dense::map<MinimizerPairType, pair<u_int64_t, u_int64_t>>().swap(_minimizerLookupTable); //_minimizerLookupTable.clear();
		ankerl::unordered_dense::map<u_int128_t, u_int32_t>().swap(_kminmer_to_abundance);
		ankerl::unordered_dense::set<u_int32_t>().swap(_isContigIndexed);
		ankerl::unordered_dense::map<u_int32_t, vector<pair<u_int32_t, u_int32_t>>>().swap(_contigIndex_to_alignments);

		Logger::get().debug() << "\tLoad contigs";
		auto start = high_resolution_clock::now();
		//processChunkParallell(IndexContigKminmerFunctor(*this));
		for(const KminmerList& kminmerList : _contigs){
			//cout << "Index contig: " << kminmerList._readIndex << endl;
			_isContigIndexed.insert(kminmerList._readIndex);
			
			for(size_t i=0; i<kminmerList._kminmersInfo.size(); i++){
			
				const ReadKminmerComplete& kminmerInfo = kminmerList._kminmersInfo[i];

				const KmerVec& vec = kminmerInfo._vec;
				_kminmer_to_unitigIndex[vec.hash128()] = -1;

			}

		}


		Logger::get().debug() << "\t\tDone: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";
		

		//Logger::get().debug() << "\tLoading alignments";
		//loadAlignments();
		//Logger::get().debug() << "\tIndexing reads";
		//start = high_resolution_clock::now();
		//indexReads();
		//Logger::get().debug() << "\t\tDone: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";


		Logger::get().debug() << "\tBreak unbridged repeats";
		start = high_resolution_clock::now();
		processChunkParallell(BreakUnbridgedRepeatsFunctor(*this));
		Logger::get().debug() << "\t\tDone: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";


		//cout << _kminmer_to_unitigIndex.size() << endl;
		//cout << _kminmer_to_abundance.size() << endl;
		//cout << _kminmer_to_readIndexes.size() << endl;
		//cout << _readIndex_to_readLength.size() << endl;
		//cout << _nbContigs << endl;
		//getchar();
	}


	template<typename Functor>
	void processChunkParallell(const Functor& functor){

		//vector<ReadType> readIndexes;
		//for(size_t i=0; i<_contigs.size(); i++){
		//	readIndexes.push_back(i);
		//}

		//cout << "single core here" << endl;
		u_int64_t i = 0;

		#pragma omp parallel num_threads(_nbCores) //_nbCores
		{

			Functor functorSub(functor);

			KminmerList contig;
			//ReadType readIndex;
			//ReadWorker& read;
			bool isDone = false;


			while(true){

				
				#pragma omp critical(processWorkerParallell)
				{

					if(i >= _contigs.size()){
						isDone = true;
					}
					else{
						contig = _contigs[i];
						i += 1;
					}
				}

				if(isDone) break;

				functorSub(contig);

			}

			
		}
	}
	*/
	/*
	class IndexContigKminmerFunctor {

		public:

		RepeatRemover& _parent;

		IndexContigKminmerFunctor(RepeatRemover& parent) : _parent(parent){
		}

		IndexContigKminmerFunctor(const IndexContigKminmerFunctor& copy) : _parent(copy._parent){
		}

		~IndexContigKminmerFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {

			ReadType readIndex = kminmerList._readIndex;
			
			for(u_int32_t i=0; i<kminmerList._kminmersInfo.size(); i++){
			
				
				const ReadKminmerComplete& kminmerInfo = kminmerList._kminmersInfo[i];

				const KmerVec& vec = kminmerInfo._vec;
				
				_parent._kminmer_to_readIndexes.lazy_emplace_l(vec.hash128(),
				[&readIndex](KminmerPosMap::value_type& v) { // key exist
				},           
				[&vec, &readIndex](const KminmerPosMap::constructor& ctor) { // key inserted
					
					vector<ReadType> readIndexes = {}; //inital count of this kminmer

					ctor(vec.hash128(), readIndexes); 

				}); // construct value_type in place when key not present

			}
			
			
		}
		
	};
	*/

	
	void fragmentContigs(){

		Logger::get().debug() << "\tIndexing initial unitigs";
		auto start = high_resolution_clock::now();
		loadUnitigIndex();
		Logger::get().debug() << "\t\tDone: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";



		Logger::get().debug() << "\tFragment contigs";

		priority_queue<FragmentCoverageDiskWriter, vector<FragmentCoverageDiskWriter> , FragmentCoverageDiskWriter_Comparator> readWriterQueue;
		u_int64_t nextReadIndexWriter = 0;
		ofstream outputFile(_inputFilenameContig + ".fragments");

		KminmerParserParallel parser(_inputFilenameContig, 0, _params._kminmerSize, false, false, _nbCores);
		parser.parse(FragmentFunctor(*this, outputFile, readWriterQueue, nextReadIndexWriter));

		outputFile.close();
		ankerl::unordered_dense::map<u_int128_t, u_int32_t>().swap(_kminmer_to_unitigIndex);
				
		Logger::get().debug() << "\t\tDone: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";

	}

	void loadUnitigIndex(){

		//const string& unitigFilename = _inputDir + "/unitig_data.txt.init";
		const string& unitigFilename = _inputDir + "/unitig_data.txt.init.k" + to_string(_params._kminmerSize);

		KminmerParserParallel parser(unitigFilename, 0, _params._kminmerSize, false, false, _nbCores);
		parser.parse(IndexUnitigFunctor(*this));
	}


	class IndexUnitigFunctor {

		public:

		RepeatRemover& _parent;

		IndexUnitigFunctor(RepeatRemover& parent) : _parent(parent){

		}

		IndexUnitigFunctor(const IndexUnitigFunctor& copy) : _parent(copy._parent){
			
		}

		~IndexUnitigFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {

			u_int32_t readIndex = kminmerList._readIndex;

			const vector<MinimizerType>& readMinimizers = kminmerList._readMinimizers;
			const vector<ReadKminmerComplete>& kminmersInfos = kminmerList._kminmersInfo;

			for(size_t i=0; i<kminmersInfos.size(); i++){

				const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];
				const KmerVec& vec = kminmerInfo._vec;
				const u_int128_t vecHash = vec.hash128();

				#pragma omp critical(IndexUnitigFunctor)
				{
					_parent._kminmer_to_unitigIndex[vecHash] = readIndex;
				}
				
			}

		}
	};

	class FragmentFunctor {

		public:

		RepeatRemover& _parent;
		ofstream& _outputFile;
		priority_queue<FragmentCoverageDiskWriter, vector<FragmentCoverageDiskWriter> , FragmentCoverageDiskWriter_Comparator>& _readWriterQueue;
		u_int64_t& _nextReadIndexWriter;

		FragmentFunctor(RepeatRemover& parent, ofstream& outputFile, auto& readWriterQueue, auto& nextReadIndexWriter) : _parent(parent), _outputFile(outputFile), _readWriterQueue(readWriterQueue), _nextReadIndexWriter(nextReadIndexWriter){
		}

		FragmentFunctor(const FragmentFunctor& copy) : _parent(copy._parent), _outputFile(copy._outputFile), _readWriterQueue(copy._readWriterQueue), _nextReadIndexWriter(copy._nextReadIndexWriter){
		}

		~FragmentFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {

			ReadType readIndex = kminmerList._readIndex;
			const vector<MinimizerType>& readMinimizers = kminmerList._readMinimizers;
			const vector<ReadKminmerComplete>& kminmersInfos = kminmerList._kminmersInfo;
			u_int8_t isCircular = kminmerList._isCircular;

			if(isCircular){
				writeContigFragments(readIndex, readMinimizers, isCircular, {});
				return;
			}


			//u_int32_t fragmentIndex = 0;
			u_int32_t lastUnitigIndex = -1;
			u_int32_t fragmentStartPos = 0;
			vector<FragmentCoverageDisk> fragments;

			for(size_t i=0; i<kminmersInfos.size(); i++){

				const KmerVec& vec = kminmersInfos[i]._vec;
				const u_int128_t& vecHash = vec.hash128();

				u_int32_t unitigIndex = -1;
				if(_parent._kminmer_to_unitigIndex.find(vecHash) != _parent._kminmer_to_unitigIndex.end()){
					unitigIndex = _parent._kminmer_to_unitigIndex[vecHash];
				}

				//const u_int32_t& unitigIndex = _parent._kminmer_to_unitigIndex[vecHash];


				
				if(unitigIndex != lastUnitigIndex || i == kminmersInfos.size()-1){

					lastUnitigIndex = unitigIndex;

					if(i == 0){
						continue;
					}

					
					u_int32_t fragmentEndPos = i-1; 
					if(i == kminmersInfos.size()-1){
						fragmentEndPos = kminmersInfos.size()-1;
					}

					u_int32_t fragmentLength = fragmentEndPos - fragmentStartPos + 1;

					FragmentCoverageDisk fragment = {fragmentStartPos, fragmentEndPos, 0};
					fragments.push_back(fragment);
					fragmentStartPos = i;
					
					
				}
			}

			writeContigFragments(readIndex, readMinimizers, isCircular, fragments);
		}


		void writeContigFragments(const u_int32_t readIndex, const vector<MinimizerType>& minimizers, const u_int8_t isCircular, const vector<FragmentCoverageDisk>& fragments){

			#pragma omp critical(writeRead)
			{
				_readWriterQueue.push({readIndex, minimizers, isCircular, fragments});
				//cout << "Add read index: " << readIndex << " " << _nextReadIndexWriter << " " << _readWriterQueue.size() << endl;


				while(!_readWriterQueue.empty()){

					const FragmentCoverageDiskWriter& readWriter = _readWriterQueue.top();

					if(readWriter._readIndex == _nextReadIndexWriter){

						u_int32_t contigSize = readWriter._minimizers.size();
						_outputFile.write((const char*)&contigSize, sizeof(contigSize));
						_outputFile.write((const char*)&readWriter._isCircular, sizeof(readWriter._isCircular));
						_outputFile.write((const char*)&readWriter._minimizers[0], contigSize*sizeof(MinimizerType));
						
						u_int32_t nbFragments = readWriter._fragments.size();
						_outputFile.write((const char*)&nbFragments, sizeof(nbFragments));

						for(size_t i=0; i<readWriter._fragments.size(); i++){
							_outputFile.write((const char*)&readWriter._fragments[i]._startPos, sizeof(readWriter._fragments[i]._startPos));
							_outputFile.write((const char*)&readWriter._fragments[i]._endPos, sizeof(readWriter._fragments[i]._endPos));
							_outputFile.write((const char*)&readWriter._fragments[i]._coverage, sizeof(readWriter._fragments[i]._coverage));
						}

						//for(size_t i=0; i<minimizers.size(); i++){
						//	_checksum += minimizers.size() * readWriter._readIndex * minimizers[i];
						//}

						//_nbContigs += 1;


						//cout << "\tConsume: " << readWriter._readIndex << endl;
						//getchar();

						_readWriterQueue.pop();
						_nextReadIndexWriter += 1;
					}
					else{
						break;
					}
				}
			}

		}
		
	};

	void computeFragmentsCoverage(){

		Logger::get().debug() << "\tIndexing kminmer abundance";
		auto start = high_resolution_clock::now();
		loadKminmerAbundance();
		Logger::get().debug() << "\t\tDone: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";

		

		Logger::get().debug() << "\tCompute fragment coverage";

		priority_queue<FragmentCoverageDiskWriter, vector<FragmentCoverageDiskWriter> , FragmentCoverageDiskWriter_Comparator> readWriterQueue;
		u_int64_t nextReadIndexWriter = 0;
		ofstream outputFile(_inputFilenameContig + ".fragments.coverage");

		//cout << "single core here" << endl;
		FragmentParserParallel parser(_inputFilenameContig + ".fragments", _params._kminmerSize, _nbCores);
		parser.parse(CoverageFunctor(*this, outputFile, readWriterQueue, nextReadIndexWriter));

		outputFile.close();
		ankerl::unordered_dense::map<u_int128_t, u_int32_t>().swap(_kminmer_to_abundance);
				
		Logger::get().debug() << "\t\tDone: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";

	}

	void loadKminmerAbundance(){

		ifstream kminmerAbundanceFile(_inputDir + "/kminmerData_abundance_init_k" + to_string(_params._kminmerSize) + ".txt");
		//ifstream kminmerAbundanceFile(_inputDir + "/kminmerData_abundance_init.txt");

		while (true) {

			u_int128_t vecHash;
			u_int32_t abundance;

			kminmerAbundanceFile.read((char*)&vecHash, sizeof(vecHash));

			if(kminmerAbundanceFile.eof()) break;

			kminmerAbundanceFile.read((char*)&abundance, sizeof(abundance));
			
			//return false;

			//bool iseof = MDBG::readKminmerAbundance(vecHash, abundance, kminmerAbundanceFile);

			//if(iseof) break;

			if(abundance <= 1) continue;
			//if(_kminmer_to_unitigIndex.find(vecHash) == _kminmer_to_unitigIndex.end()) continue;

			_kminmer_to_abundance[vecHash] = abundance;

			//cout << "lala" << " " << abundance << endl;

		}

		kminmerAbundanceFile.close();
	}


	class CoverageFunctor {

		public:

		RepeatRemover& _parent;
		ofstream& _outputFile;
		priority_queue<FragmentCoverageDiskWriter, vector<FragmentCoverageDiskWriter> , FragmentCoverageDiskWriter_Comparator>& _readWriterQueue;
		u_int64_t& _nextReadIndexWriter;

		CoverageFunctor(RepeatRemover& parent, ofstream& outputFile, auto& readWriterQueue, auto& nextReadIndexWriter) : _parent(parent), _outputFile(outputFile), _readWriterQueue(readWriterQueue), _nextReadIndexWriter(nextReadIndexWriter){
		}

		CoverageFunctor(const CoverageFunctor& copy) : _parent(copy._parent), _outputFile(copy._outputFile), _readWriterQueue(copy._readWriterQueue), _nextReadIndexWriter(copy._nextReadIndexWriter){
		}

		~CoverageFunctor(){
		}

		void operator () (const KminmerList& kminmerList, const vector<FragmentCoverageDisk>& fragments) {

			ReadType readIndex = kminmerList._readIndex;
			const vector<MinimizerType>& readMinimizers = kminmerList._readMinimizers;
			const vector<ReadKminmerComplete>& kminmersInfos = kminmerList._kminmersInfo;
			u_int8_t isCircular = kminmerList._isCircular;

			//cout << readIndex << "\t" << readMinimizers.size() << "\t" << kminmersInfos.size() << "\t" << fragments.size() << endl;

			if(isCircular){
				writeContigFragments(readIndex, readMinimizers, isCircular, {});
				return;
			}
			
			vector<FragmentCoverageDisk> fragmentCoverages;

			for(size_t i=0; i<fragments.size(); i++){

				//cout << "\t" << fragments[i]._startPos << "\t" << fragments[i]._endPos << endl;
				long double sum = 0;
				long double n = 0;

				for(size_t j=fragments[i]._startPos; j<=fragments[i]._endPos; j++){
					
					//if(j >= kminmersInfos.size()){
					//	cout << "derp" << endl;
					//}
					const KmerVec& vec = kminmersInfos[j]._vec;
					const u_int128_t& vecHash = vec.hash128();

					if(_parent._kminmer_to_abundance.find(vecHash) == _parent._kminmer_to_abundance.end()){
						sum += 1;
					}
					else{
						sum += _parent._kminmer_to_abundance[vecHash];
					}

					n += 1;
				}


				long double fragmentCoverage = 0;
				if(n > 0){
					fragmentCoverage = sum / n;
				}

				fragmentCoverages.push_back({fragments[i]._startPos, fragments[i]._endPos, (float)fragmentCoverage});
			
			}


			writeContigFragments(readIndex, readMinimizers, isCircular, fragmentCoverages);

		}




		void writeContigFragments(const u_int32_t readIndex, const vector<MinimizerType>& minimizers, const u_int8_t isCircular, const vector<FragmentCoverageDisk>& fragments){

			#pragma omp critical(writeRead)
			{
				_readWriterQueue.push({readIndex, minimizers, isCircular, fragments});
				//cout << "Add read index: " << readIndex << " " << _nextReadIndexWriter << " " << _readWriterQueue.size() << endl;


				while(!_readWriterQueue.empty()){

					const FragmentCoverageDiskWriter& readWriter = _readWriterQueue.top();

					if(readWriter._readIndex == _nextReadIndexWriter){

						u_int32_t contigSize = readWriter._minimizers.size();
						_outputFile.write((const char*)&contigSize, sizeof(contigSize));
						_outputFile.write((const char*)&readWriter._isCircular, sizeof(readWriter._isCircular));
						_outputFile.write((const char*)&readWriter._minimizers[0], contigSize*sizeof(MinimizerType));
						
						u_int32_t nbFragments = readWriter._fragments.size();
						_outputFile.write((const char*)&nbFragments, sizeof(nbFragments));

						for(size_t i=0; i<readWriter._fragments.size(); i++){
							_outputFile.write((const char*)&readWriter._fragments[i]._startPos, sizeof(readWriter._fragments[i]._startPos));
							_outputFile.write((const char*)&readWriter._fragments[i]._endPos, sizeof(readWriter._fragments[i]._endPos));
							_outputFile.write((const char*)&readWriter._fragments[i]._coverage, sizeof(readWriter._fragments[i]._coverage));
						}

						//for(size_t i=0; i<minimizers.size(); i++){
						//	_checksum += minimizers.size() * readWriter._readIndex * minimizers[i];
						//}

						//_nbContigs += 1;


						//cout << "\tConsume: " << readWriter._readIndex << endl;
						//getchar();

						_readWriterQueue.pop();
						_nextReadIndexWriter += 1;
					}
					else{
						break;
					}
				}
			}

		}
	
	};


	void breakUnbridgedRepeats(){


		Logger::get().debug() << "\tLoading alignments";
		loadAlignments();
		auto start = high_resolution_clock::now();
		Logger::get().debug() << "\t\tDone: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";


		Logger::get().debug() << "\tBreak unbridged repeats";

		priority_queue<ContigsWriter, vector<ContigsWriter> , ContigsWriter_Comparator> readWriterQueue;
		u_int64_t nextReadIndexWriter = 0;
		ofstream outputFile(_outputFilenameContig);

		//cout << "single core here" << endl;
		FragmentParserParallel parser(_inputFilenameContig + ".fragments.coverage", _params._kminmerSize, _nbCores);
		parser.parse(BreakUnbridgedRepeatsFunctor(*this, outputFile, readWriterQueue, nextReadIndexWriter));

		outputFile.close();
		ankerl::unordered_dense::map<u_int128_t, u_int32_t>().swap(_kminmer_to_abundance);
				
		Logger::get().debug() << "\t\tDone: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";


	}

	//void loadContigs(){

	//}
	/*
	void indexReads(){

		//cout << "todo: pour ont utiliser read_data_corrected" << endl;

		//Logger::get().debug() << "Indexing reads 1";
		//KminmerParserParallel parser1(_inputFilenameContig, -1, _kminmerSize, false, false, _nbCores);
		//parser._densityThreshold = _minimizerDensity;
		//parser1.parse(IndexContigKminmerFunctor(*this));
		
		//if(_hasQuality){
		//	KminmerParserParallel parser2(_inputDir + "/read_data_corrected.txt", -1, _kminmerSize, false, false, _nbCores);
		//	parser2._densityThreshold = _minimizerDensity;
		//	parser2.parse(IndexReadsFunctor(*this));
		//}
		//else{
			KminmerParserParallel parser2(_inputDir + "/read_data_init.txt", -1, _kminmerSize, false, true, _nbCores);
			parser2._densityThreshold = _minimizerDensity;
			parser2.parse(IndexReadsFunctor(*this));
		//}
		//Logger::get().debug() << "Indexing reads 2";


	}




	class IndexReadsFunctor {

		public:

		RepeatRemover& _parent;

		IndexReadsFunctor(RepeatRemover& parent) : _parent(parent){
		}

		IndexReadsFunctor(const IndexReadsFunctor& copy) : _parent(copy._parent){
		}

		~IndexReadsFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {

			ReadType readIndex = kminmerList._readIndex;
			//if(readIndex % 100000 == 0)  Logger::get().debug() << "\tIndexing reads: " << readIndex;

			vector<float> abundances;
			unordered_set<u_int32_t> distinctUnitigIndex;

			for(size_t i=0; i<kminmerList._kminmersInfo.size(); i++){

				const ReadKminmerComplete& kminmerInfo = kminmerList._kminmersInfo[i];
				const KmerVec& vec = kminmerInfo._vec;
				const u_int128_t& vecHash = vec.hash128();

				if(_parent._kminmer_to_unitigIndex.find(vecHash) == _parent._kminmer_to_unitigIndex.end()) continue;

				distinctUnitigIndex.insert(_parent._kminmer_to_unitigIndex[vecHash]);

				if(_parent._kminmer_to_abundance.find(vecHash) != _parent._kminmer_to_abundance.end()){
					abundances.push_back(_parent._kminmer_to_abundance[vecHash]);
				}
			}

			

			if(distinctUnitigIndex.size() <= 1) return;

			for(size_t i=0; i<kminmerList._kminmersInfo.size(); i++){
			
				
				const ReadKminmerComplete& kminmerInfo = kminmerList._kminmersInfo[i];

				const KmerVec& vec = kminmerInfo._vec;
				const u_int128_t& vecHash = vec.hash128();
				
				_parent._kminmer_to_readIndexes.modify_if(vecHash, 
					[&readIndex](KminmerPosMap::value_type& v) { 
				
					v.second.push_back(readIndex);
						
				});

			}
			

			#pragma omp critical(IndexReadsFunctor)
			{
				_parent._readIndex_to_readLength[readIndex] = kminmerList._kminmersInfo.size();
			}

			
		}
		
	};

	*/

	
	class BreakUnbridgedRepeatsFunctor {

		public:

		RepeatRemover& _parent;
		ofstream& _outputFile;
		priority_queue<ContigsWriter, vector<ContigsWriter> , ContigsWriter_Comparator>& _readWriterQueue;
		u_int64_t& _nextReadIndexWriter;

		BreakUnbridgedRepeatsFunctor(RepeatRemover& parent, ofstream& outputFile, auto& readWriterQueue, auto& nextReadIndexWriter) : _parent(parent), _outputFile(outputFile), _readWriterQueue(readWriterQueue), _nextReadIndexWriter(nextReadIndexWriter){
		}

		BreakUnbridgedRepeatsFunctor(const BreakUnbridgedRepeatsFunctor& copy) : _parent(copy._parent), _outputFile(copy._outputFile), _readWriterQueue(copy._readWriterQueue), _nextReadIndexWriter(copy._nextReadIndexWriter){
		}

		~BreakUnbridgedRepeatsFunctor(){
		}

		void operator () (const KminmerList& kminmerList, const vector<FragmentCoverageDisk>& inputFragments) {

			ReadType readIndex = kminmerList._readIndex;
			const vector<MinimizerType>& readMinimizers = kminmerList._readMinimizers;
			const vector<ReadKminmerComplete>& kminmersInfos = kminmerList._kminmersInfo;
			u_int8_t isCircular = kminmerList._isCircular;

			if(isCircular){
				writeContig(readIndex, {readMinimizers}, isCircular);
				return;
			}

			vector<Fragment> fragments;

			for(u_int32_t i=0; i<inputFragments.size(); i++){

				const FragmentCoverageDisk& fragment = inputFragments[i];
				u_int32_t fragmentLength = fragment._endPos - fragment._startPos + 1;

				Fragment fragment2 = {i, fragment._startPos, fragment._endPos, fragmentLength, fragment._coverage, -1, {}};
				fragments.push_back(fragment2);
			}


			if(fragments.size() == 0){
				writeContig(readIndex, {readMinimizers}, isCircular);
				return;
			}

			/*
			//fragments[0]._startPos = 0;
			//fragments[fragments.size()-1]._endPos = kminmersInfos.size()-1; //contigSequence.size();

			for(size_t i=0; i<fragments.size(); i++){

				Fragment& fragment = fragments[i];
				
				for(size_t i=fragment._startPos; i <= fragment._endPos; i++){
					const KmerVec& vec = kminmersInfos[i]._vec;
					const u_int128_t& vecHash = vec.hash128();

					if(_parent._kminmer_to_unitigIndex.find(vecHash) == _parent._kminmer_to_unitigIndex.end()) continue;
					if(_parent._kminmer_to_readIndexes.find(vecHash) == _parent._kminmer_to_readIndexes.end()) continue;

					for(const ReadType& readIndex : _parent._kminmer_to_readIndexes[vecHash]){
						fragment._readIndexes.push_back(readIndex);
					}

					std::sort(fragment._readIndexes.begin(), fragment._readIndexes.end());
				}

			}
			*/

			if(_parent._contigIndex_to_alignments.find(readIndex) == _parent._contigIndex_to_alignments.end()){
				//cout << "contig has no alignments" << endl;
				writeContig(readIndex, {readMinimizers}, isCircular);
				return;
			}

			//cout << "a" << endl;
			//cout << (_parent._contigIndex_to_alignments.find(readIndex) != _parent._contigIndex_to_alignments.end()) << endl;

			const vector<pair<u_int32_t, u_int32_t>>& alignments = _parent._contigIndex_to_alignments[readIndex];
			if(alignments.size() == 0){
				//cout << "contig has no alignments" << endl;
				writeContig(readIndex, {readMinimizers}, isCircular);
				return;
			}

			//cout << "b" << endl;
			computeBridgingReads(readIndex, fragments, alignments);

			bool print = false;

			if(_parent._debugFindFinalContigIndex != -1 && readMinimizers.size() > _parent._debugFindFinalContigMinimizers.size()/2){
				float shared = 0;
				for(const auto& m : readMinimizers){
					if(_parent._debugFindFinalContigMinimizers.find(m) != _parent._debugFindFinalContigMinimizers.end()) shared += 1;
				}
				float sharedFraction = shared / readMinimizers.size();
				if(sharedFraction > 0.5){
					cout << sharedFraction << " " << readMinimizers.size() << "    ctg" << readIndex << " " << ((int)kminmerList._isCircular) << endl;

					//for(size_t i=0; i<readMinimizers.size(); i++){
					//	cout << i << "\t" << readMinimizers[i] << endl;
					//}
					//getchar();
					print = true;
				}
			}
			//cout << "c" << endl;
			
			if(print){
				for(size_t i=0; i<fragments.size(); i++){
					const Fragment& fragment = fragments[i];
					cout << fragment._fragmentIndex << "\t" << fragment._startPos << "\t" << fragment._endPos << "\t" << fragment._length << "\t" << (int) fragment._coverage << "\t" << fragment._readIndexes.size() << "\t" << fragment._nbBridgingReads.size() << endl; //<< "\t" << fragment._readIndexes.size() << endl; //<< "\t" << getCoverageInBounds(repeats, fragmentStartPos, fragmentEndPos) << endl;
				
				}
			}
			

			//getchar();
			vector<FragmentPath> paths;

			for(const Fragment& fragment : fragments){
				if(fragment._length * 1/_parent._params._minimizerDensity_assembly < 10000) continue;

				//cout << "\tSource fragment: " << fragment._fragmentIndex << " " << fragment._startPos << " " << fragment._coverage << endl;

				const FragmentPath& path = getCovPath(fragment, fragments);
				paths.push_back(path);

			}

			std::sort(paths.begin(), paths.end(), [](const FragmentPath& a, const FragmentPath& b){
				return a.getPathLength() < b.getPathLength();
			});


			//cout << endl;
			for(size_t i=0; i<paths.size(); i++){
				//cout << i << " " << paths[i]._fragmentIndexStart << " " << paths[i]._fragmentIndexEnd << endl;

				for(size_t j = paths[i]._fragmentIndexStart; j <= paths[i]._fragmentIndexEnd; j++){
					if(fragments[j]._finalContigIndex != -1) continue;
					fragments[j]._finalContigIndex = i;
				}
			}

			//cout << endl;



			u_int32_t currentContigIndex = fragments[0]._finalContigIndex;

			Fragment finalDummyFragment = {fragments.size(), 0, 0, 0, 0, -2};
			fragments.push_back(finalDummyFragment);

			//cout << endl;
			//cout << readIndex << " " << readMinimizers.size() << endl;
			//for(size_t i=0; i<fragments.size(); i++){
			//	const Fragment& fragment = fragments[i];
			//	cout << "\t" << fragment._fragmentIndex << "\t" << fragment._startPos << "\t" << fragment._endPos << "\t" << fragment._length << "\t" << (int) fragment._coverage << "\t" << fragment._readIndexes.size() << "\t" << fragment._nbBridgingReads.size() << "\t" << fragment._finalContigIndex << endl; //<< "\t" << fragment._readIndexes.size() << endl; //<< "\t" << getCoverageInBounds(repeats, fragmentStartPos, fragmentEndPos) << endl;
			//}

			int32_t nbContigsFinal = 0;

			for(size_t i=0; i<fragments.size(); i++){

				const Fragment& fragment = fragments[i];

				if(fragments[i]._finalContigIndex != currentContigIndex){
					currentContigIndex = fragments[i]._finalContigIndex;
					nbContigsFinal += 1;
				}
			}

			if(nbContigsFinal > 1){
				isCircular = 0;
				//cout << "Split: " << readIndex << endl;
			}

			//cout << "\tNb final contigs: " << nbContigsFinal << endl;
			vector<vector<MinimizerType>> contigs;

			int64_t startPos = 0;
			currentContigIndex = fragments[0]._finalContigIndex;

			for(size_t i=0; i<fragments.size(); i++){

				const Fragment& fragment = fragments[i];

				if(fragments[i]._finalContigIndex != currentContigIndex){
					
					//const string& originalContigSeq =  _contigSequences[contigIndex];

					int64_t endPos = fragments[i-1]._endPos;
					//if(i == fragments.size()-1) endPos = originalContigSeq.size();

					//if(endPos-startPos >= _kminmerSizeRepeat){

						//const string& contigSeq = originalContigSeq.substr(startPos, endPos-startPos);

						//startPos = 0;
						//endPos = 0;


						//cout << "\tWrite contig: " <<  startPos << "\t" << endPos << "\t" << endPos-startPos << "\t" << fragments[i-1]._coverage << endl; 
						//cout << "Write contig: " << endPos-startPos << "\t" << kminmersInfos.size() << endl; 

						vector<MinimizerType>::const_iterator first = readMinimizers.begin() + startPos;

						//int endPos2 = endPos+1;
						//if(fragments[i]._finalContigIndex == -2){ //final fragment
						//	endPos2 = endPos + _parent._kminmerSize;
						//}

						vector<MinimizerType>::const_iterator last = readMinimizers.begin() + endPos + _parent._params._kminmerSize;
						vector<MinimizerType> contigMinimizers(first, last);


						//for(size_t i=0; i<contigMinimizers.size(); i++){
						//	cout << i <<" " << contigMinimizers[i] << endl;
						//}

						//const vector<MinimizerType>& contigMinimizers = readMinimizers.substr(startPos, endPos-startPos+_kminmzerSize-1);
						//cout << contigMinimizers.size() << endl;
						//writeContig(readIndex, contigMinimizers, isCircular);
						contigs.push_back(contigMinimizers);
						//writeRepeatContig(splitIndex, _contigHeaders[contigIndex], contigSeq, originalContigSeq);
						//splitIndex += 1;

					//}

					startPos = fragments[i]._startPos;
					currentContigIndex = fragments[i]._finalContigIndex;

				}

				//cout << fragment._fragmentIndex << "\t" << fragment._startPos << "\t" << fragment._endPos << "\t" << fragment._length << "\t" << fragment._coverage << "\t" << fragment._finalContigIndex << endl; //<< "\t" << getCoverageInBounds(repeats, fragmentStartPos, fragmentEndPos) << endl;

			}
			
			
			writeContig(readIndex, contigs, isCircular);
			//getchar();
		}


		void computeBridgingReads(const u_int32_t readIndex, vector<Fragment>& fragments, const vector<pair<u_int32_t, u_int32_t>>& alignments){

			//vector<u_int32_t> coverages(5000000, 0);

			for(const auto& alignment : alignments){

				//for(size_t i=alignment.first; i<alignment.second; i++){
				//	coverages[i] += 1;
				//}

				vector<u_int32_t> mappedFragments;

				for(const Fragment& fragment : fragments){
					if(alignment.first < fragment._startPos && alignment.second > fragment._endPos){ //fragment contained
						mappedFragments.push_back(fragment._fragmentIndex);
					}
					else if(alignment.first > fragment._startPos && alignment.first < fragment._endPos){ //fragment overlap on the left
						mappedFragments.push_back(fragment._fragmentIndex);
					}
					else if(alignment.second > fragment._startPos && alignment.second < fragment._endPos){ //fragment overlap on the right
						mappedFragments.push_back(fragment._fragmentIndex);
					}
				}

				if(mappedFragments.size() <= 1) continue;

				//cout << "---" << endl;
				for(size_t i=0; i<mappedFragments.size(); i++){
					
					//cout << mappedFragments[i] << endl;
					Fragment& fragment1 = fragments[mappedFragments[i]];

					for(size_t j=i+1; j<mappedFragments.size(); j++){
						
						Fragment& fragment2 = fragments[mappedFragments[j]];

						fragment1._nbBridgingReads[fragment2._fragmentIndex] += 1;
						fragment2._nbBridgingReads[fragment1._fragmentIndex] += 1;
					}
				}
			}

		}


		
		FragmentPath getCovPath(const Fragment& sourceFragment, vector<Fragment>& fragments){

			long double sourceCoverage = sourceFragment._coverage;
			//float minRepeatCoverage = sourceCoverage * 2.5;

			int64_t maxFragmentIndex = 0;
			int64_t minFragmentIndex = 0;
			long double currentSourceCoverage = sourceCoverage;

			while(true){

				
				//cout << "\tCurrent cov: " << currentSourceCoverage << endl;

				long double loopCov = currentSourceCoverage;

				maxFragmentIndex = getCovPath_direction(sourceFragment, fragments, currentSourceCoverage, sourceCoverage, true);
				minFragmentIndex = getCovPath_direction(sourceFragment, fragments, currentSourceCoverage, sourceCoverage, false);


				if(currentSourceCoverage == loopCov) break;

				//cout << "\tCheck: " << currentSourceCoverage << endl;
				//getchar();
				
			}
			//int64_t minFragmentIndex = 0;

			//cout << "\t\t" << sourceFragment._fragmentIndex << ": " << minFragmentIndex << " " << maxFragmentIndex << endl;

			FragmentPath path = {minFragmentIndex, maxFragmentIndex};

			return path;
		}

		int64_t getCovPath_direction(const Fragment& sourceFragment, vector<Fragment>& fragments, long double& sourceCoverage, const long double& sourceCoverageInitial, const bool& useSuccessor){

			vector<u_int32_t> specificFragmentIndexes;
			specificFragmentIndexes.push_back(sourceFragment._fragmentIndex);
			//int64_t currentFragmentIndex =  sourceFragment._fragmentIndex;
			//int64_t maxI = sourceFragment._fragmentIndex + 1;

			//int i =0;

			while(true){

				

				const int64_t& nextFragmentIndex = getNextSpecificFragmentIndex(fragments, specificFragmentIndexes, sourceCoverage, useSuccessor);
				
				//cout << "\t\t" << i << " " << nextFragmentIndex << endl;
				//i += 1;

				if(nextFragmentIndex == -1) break;

				if(fragments[nextFragmentIndex]._coverage > sourceCoverage && fragments[nextFragmentIndex]._coverage < sourceCoverageInitial*1.5f){
					sourceCoverage = fragments[nextFragmentIndex]._coverage;
					//cout << "\tCov changed: " << sourceCoverage << endl;
					return -1;
				}

				specificFragmentIndexes.push_back(nextFragmentIndex);
			}

			return specificFragmentIndexes[specificFragmentIndexes.size()-1];
			/*
			while(true){

				if(useSuccessor){
					if(currentFragmentIndex == fragments.size()-1) break;
				}
				else{
					if(currentFragmentIndex == 0) break;
				}

				int64_t nextFragmentIndex = getNextSpecificFragmentIndex(fragments, minRepeatCoverage, currentFragmentIndex, useSuccessor);
				//cout << "\t\tRepeat: " << currentFragmentIndex << " " << nextFragmentIndex << endl;
				if(nextFragmentIndex == -1) break;

				//int nbBridgingReads = getNbBridgingReads(fragments[currentFragmentIndex], fragments[nextFragmentIndex]);
				//if(nbBridgingReads < 2) break;

				currentFragmentIndex = nextFragmentIndex;
			}

			return currentFragmentIndex;
			*/
		}

		int64_t getNextSpecificFragmentIndex(vector<Fragment>& fragments, vector<u_int32_t>& specificFragmentIndexes, long double& sourceCoverage, const bool& useSuccessor){
		
			float minRepeatCoverage = sourceCoverage * 2.0;

			//for(const u_int32_t& specificFragmentIndex : specificFragmentIndexes){
			for(int64_t ii=specificFragmentIndexes.size()-1; ii >= 0; ii--){

				const u_int32_t& specificFragmentIndex = specificFragmentIndexes[ii];


				Fragment& sourceFragment = fragments[specificFragmentIndex];
				const int64_t& latestSpecificIndex = specificFragmentIndexes[specificFragmentIndexes.size()-1];

				if(useSuccessor){


					for(int64_t i=latestSpecificIndex+1; i<fragments.size(); i++){
						
						const Fragment& fragment = fragments[i];

						if(fragment._coverage >= minRepeatCoverage) continue; //Repeat fragment
						if(sourceFragment._fragmentIndex+1 == fragment._fragmentIndex) return i; //Just two successive fragments without a repeat in between
						
						const u_int32_t& nbBridingReads = sourceFragment._nbBridgingReads[fragment._fragmentIndex];
						
						if(nbBridingReads == 0) continue;

						//if(bridgeLength * 1/_parent._minimizerDensity > 100000) continue;
						//if(std::find(specificFragmentIndexes.begin(), specificFragmentIndexes.end(), fragment._fragmentIndex) != specificFragmentIndexes.end()) continue;
						
						//cout <<"\t\t\tFound bridge: " <<  sourceFragment._fragmentIndex << " -> " << i << " " << fragment._coverage << endl; 
						return i;
						//if(fragment._length <= 1) continue;
						//if(fragment._coverage < minRepeatCoverage) return i;
						//if(fragment._coverage < minRepeatCoverage){
						//	int nbBridgingReads = getNbBridgingReads(fragments[sourceFragmentIndex], fragment);
							//cout << "\tRepeat: " << currentFragmentIndex << " " << nextFragmentIndex << " " << nbBridgingReads << endl;
						//	if(nbBridgingReads >= 1) return i;
						//}


					}

				}
				else{

					for(int64_t i=latestSpecificIndex-1; i >= 0; i--){
						
						const Fragment& fragment = fragments[i];

						if(fragment._coverage >= minRepeatCoverage) continue; //Repeat fragment
						if(sourceFragment._fragmentIndex == fragment._fragmentIndex+1) return i; //Just two successive fragments without a repeat in between

						const u_int32_t& nbBridingReads = sourceFragment._nbBridgingReads[fragment._fragmentIndex];
						
						if(nbBridingReads == 0) continue;

						//if(bridgeLength * 1/_parent._minimizerDensity > 100000) continue;
						//if(std::find(specificFragmentIndexes.begin(), specificFragmentIndexes.end(), fragment._fragmentIndex) != specificFragmentIndexes.end()) continue;
						
						//cout <<"\t\t\tFound bridge: " <<  sourceFragment._fragmentIndex << " -> " << i << " " << fragment._coverage << endl; 
						return i;
						//if(fragment._length <= 1) continue;
						//if(fragment._coverage < minRepeatCoverage) return i;
						//if(fragment._coverage < minRepeatCoverage){
						//	int nbBridgingReads = getNbBridgingReads(fragments[sourceFragmentIndex], fragment);
							//cout << "\tRepeat: " << currentFragmentIndex << " " << nextFragmentIndex << " " << nbBridgingReads << endl;
						//	if(nbBridgingReads >= 1) return i;
						//}


					}

				}

			}

			return -1;
		}
		
		/*
		int64_t getNextSpecificFragmentIndex(const vector<Fragment>& fragments, const float& minRepeatCoverage, const int64_t& sourceFragmentIndex, const bool& useSuccessor){

			if(useSuccessor){
				for(int64_t i=sourceFragmentIndex+1; i<fragments.size(); i++){
					const Fragment&  fragment = fragments[i];
					//if(fragment._length <= 1) continue;
					//if(fragment._coverage < minRepeatCoverage) return i;
					if(fragment._coverage < minRepeatCoverage){
						int nbBridgingReads = getNbBridgingReads(fragments[sourceFragmentIndex], fragment);
						//cout << "\tRepeat: " << currentFragmentIndex << " " << nextFragmentIndex << " " << nbBridgingReads << endl;
						if(nbBridgingReads >= 1) return i;
					}


				}

			}
			else{

				for(int64_t i=sourceFragmentIndex-1; i >= 0; i--){
					const Fragment&  fragment = fragments[i];
					//if(fragment._length <= 1) continue;
					//if(fragment._coverage < minRepeatCoverage) return i;
					if(fragment._coverage < minRepeatCoverage){
						int nbBridgingReads = getNbBridgingReads(fragments[sourceFragmentIndex], fragment);
						//cout << "\tRepeat: " << currentFragmentIndex << " " << nextFragmentIndex << " " << nbBridgingReads << endl;
						if(nbBridgingReads >= 1) return i;
					}
				}

			}

			return -1;
		}
		*/

		u_int64_t getNbBridgingReads(const Fragment& specificFragment1, const Fragment& specificFragment2){


			size_t i=0;
			size_t j=0;
			u_int64_t nbShared = 0;

			while(i < specificFragment1._readIndexes.size() && j < specificFragment2._readIndexes.size()){

				if(specificFragment1._readIndexes[i] == specificFragment2._readIndexes[j]){
					nbShared += 1;
					i += 1;
					j += 1;
					
				}
				else if(specificFragment1._readIndexes[i] < specificFragment2._readIndexes[j]){
					i += 1;
				}
				else{
					j += 1;
				}

			}

			return nbShared;

		}
		
		void writeContig(const u_int64_t readIndex, const vector<vector<MinimizerType>>& contigs, const u_int8_t& isCircular){
	
			writeRead(readIndex, contigs, isCircular);
			/*
			#pragma omp critical(writeContig)
			{
				for(const vector<MinimizerType>& minimizers : contigs){
					_parent._nbContigs += 1;

					u_int32_t contigSize = minimizers.size();
					_parent._outputFile.write((const char*)&contigSize, sizeof(contigSize));
					_parent._outputFile.write((const char*)&isCircular, sizeof(isCircular));
					_parent._outputFile.write((const char*)&minimizers[0], contigSize*sizeof(MinimizerType));
				}
			}
			*/
			
		}

		void writeRead(u_int64_t readIndex, const vector<vector<MinimizerType>>& contigSequences, u_int8_t isCircular){

			#pragma omp critical(writeRead)
			{
				_readWriterQueue.push({readIndex, contigSequences, isCircular});
				//cout << "Add read index: " << readIndex << " " << _nextReadIndexWriter << " " << _readWriterQueue.size() << endl;


				while(!_readWriterQueue.empty()){

					const ContigsWriter& readWriter = _readWriterQueue.top();

					if(readWriter._readIndex == _nextReadIndexWriter){

						for(const vector<MinimizerType>& minimizers : readWriter._contigs){
							if(minimizers.size() > 0){
								u_int32_t contigSize = minimizers.size();
								_outputFile.write((const char*)&contigSize, sizeof(contigSize));
								_outputFile.write((const char*)&readWriter._isCircular, sizeof(readWriter._isCircular));
								_outputFile.write((const char*)&minimizers[0], contigSize*sizeof(MinimizerType));

								for(size_t i=0; i<minimizers.size(); i++){
									_parent._checksum += minimizers.size() * readWriter._readIndex * minimizers[i];
								}

								//_nbContigs += 1;

							}
						}

						//cout << "\tConsume: " << readWriter._readIndex << endl;
						//getchar();

						_readWriterQueue.pop();
						_nextReadIndexWriter += 1;
					}
					else{
						break;
					}
				}
			}

		}


	};


	class FragmentParserParallel{

		public:

		string _inputFilename;
		size_t _k;
		int _nbCores;

		FragmentParserParallel(){
		}

		FragmentParserParallel(const string& inputFilename, size_t k, int nbCores){

			if(!fs::exists(inputFilename)){
				cout << "File not found: " << inputFilename << endl;
				exit(1);
			}

			_inputFilename = inputFilename;
			_k = k;
			_nbCores = nbCores;
		}

		template<typename Functor>
		void parse(const Functor& functor){

			ifstream file_readData(_inputFilename);

			u_int64_t readIndex = -1;

			#pragma omp parallel num_threads(_nbCores)
			{

				bool isEOF = false;
				Functor functorSub(functor);
				vector<MinimizerType> minimizers;
				vector<u_int32_t> minimizerPos;
				vector<u_int8_t> minimizerQualities;
				vector<u_int8_t> minimizerDirections;
				u_int32_t size;
				KminmerList kminmerList;
				u_int8_t isCircular;
				u_int32_t nbFragments;
				vector<FragmentCoverageDisk> fragments;

				while(true){
					

					#pragma omp critical(KminmerParserParallel_parse)
					{

						
						file_readData.read((char*)&size, sizeof(size));

						if(file_readData.eof()) isEOF = true;

						if(!isEOF){

							readIndex += 1;

							kminmerList = {readIndex};

							minimizers.resize(size);
							minimizerPos.resize(size);
							minimizerQualities.resize(size);
							minimizerDirections.resize(size);
							
							file_readData.read((char*)&isCircular, sizeof(isCircular));

							file_readData.read((char*)&minimizers[0], size*sizeof(MinimizerType));


							file_readData.read((char*)&nbFragments, sizeof(nbFragments));

							fragments.resize(nbFragments);

							for(size_t i=0; i<nbFragments; i++){

								u_int32_t startPos;
								file_readData.read((char*)&startPos, sizeof(startPos));
								
								u_int32_t endPos;
								file_readData.read((char*)&endPos, sizeof(endPos));

								float coverage;
								file_readData.read((char*)&coverage, sizeof(coverage));

								fragments[i] = {startPos, endPos, coverage};
							}

						}

					}


					if(isEOF) break;


					vector<ReadKminmerComplete> kminmersInfo;
					MDBG::getKminmers_complete(_k, minimizers, minimizerPos, kminmersInfo, readIndex, minimizerQualities);
					
					kminmerList._readMinimizers = minimizers;
					kminmerList._minimizerPos = minimizerPos;
					kminmerList._readMinimizerDirections = minimizerDirections;
					kminmerList._readQualities = minimizerQualities;
					kminmerList._kminmersInfo = kminmersInfo;
					kminmerList._isCircular = isCircular;
					functorSub(kminmerList, fragments);
				}
			}

			file_readData.close();


		}
	};

    //struct ReadWriter{
    //    u_int64_t _readIndex;
	//	vector<MinimizerType> _contigSequence;
	//	bool _isCircular;
    //};




	
};	


#endif 


