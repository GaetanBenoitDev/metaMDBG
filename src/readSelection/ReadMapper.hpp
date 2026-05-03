
#ifndef MDBG_METAG_READMAPPER
#define MDBG_METAG_READMAPPER

#include "../Commons.hpp"
#include "../readSelection/MinimizerChainer.hpp"


class ReadMapper{

	public:

	string _inputDir;
	string _alignmentFilename;
	size_t _totalNbReads;
	u_int64_t _maxChainingBand;
	u_int64_t _usedCoverageForCorrection;
	u_int64_t _nbMinimizersPerChunk;
	int _nbCores;

	size_t _kminmerSize;
	ofstream _alignmentFile;
	u_int64_t _nbReads;
	u_int64_t _alignmentCheckSum;
	int _nbPartitions;
	vector<BGZF*> _partitionFiles;
	string _partitionDir;


	struct MinimizerPosition2{
		MinimizerPairType _minimizerPair;
		ReadType _readIndex;
		u_int32_t _position;
		u_int32_t _positionIndex;
		bool _isReversed;

		friend bool operator<(const MinimizerPosition2 &a, const MinimizerPosition2 &b){
			return a._minimizerPair < b._minimizerPair;
		}
	};

	vector<MinimizerPosition2> _allMinimizerPositions;
	ankerl::unordered_dense::map<MinimizerPairType, pair<u_int64_t, u_int64_t>> _minimizerLookupTable;
	ofstream _debugFile;
	u_int64_t _nbReadsProcessed;
	u_int64_t _nbReadsInChunk;
	//phmap::flat_hash_map<MinimizerPairType, pair<u_int64_t, u_int64_t>> _minimizerLookupTable;

	struct Anchor2{

		ReadType _readIndex;
		int64_t _referencePosition;
		u_int32_t _referencePositionIndex;
		int64_t _queryPosition;
		u_int32_t _queryPositionIndex;
		bool _isReversed;
		//int32_t _referencePositionIndex;
		//int32_t _queryPositionIndex;

		Anchor2(const ReadType readIndex, const int64_t referencePosition, const int64_t referencePositionIndex, const int64_t queryPosition, const u_int32_t queryPositionIndex, const bool isReversed){
			_readIndex = readIndex;
			_referencePosition = referencePosition;
			_referencePositionIndex = referencePositionIndex;
			_queryPosition = queryPosition;
			_queryPositionIndex = queryPositionIndex;
			_isReversed = isReversed;
			//_referencePositionIndex = referencePositionIndex;
			//_queryPositionIndex = queryPositionIndex;
		}
	};

	struct AlignmentScore{
		ReadType _queryReadIndex;
		int32_t _score;
		//u_int32_t queryStart;
		//u_int32_t queryEnd;
		//bool _isQueryReversed;
	};

	struct AlignmentScoreMatches{
		ReadType _queryReadIndex;
		//int32_t _score;
		u_int32_t _nbMatches;
		vector<unsigned char> _compressedMatchPositions;
		//u_int32_t queryStart;
		//u_int32_t queryEnd;
		//bool _isQueryReversed;
	};

	struct AlignmentScoreComparator {
		bool operator()(AlignmentScore const& a, AlignmentScore const& b){
			if(a._score == b._score){
				return a._queryReadIndex < b._queryReadIndex;
			}
			return a._score > b._score;
		}
	};

	typedef priority_queue<AlignmentScore, vector<AlignmentScore> , AlignmentScoreComparator> AlignmentScoreQueue;


	struct Chain2{
		float _chainingScore;
		//vector<size_t> _anchorInterval;
		vector<u_int32_t> _queryAnchorPositions;
		vector<u_int32_t> _referenceAnchorPositions;
		int64_t _nbMatches;
		int64_t _nbDifferencesQuery;
		int64_t _nbDifferencesReference;
		int64_t _queryStart;
		int64_t _queryEnd;
		int64_t _referenceStart;
		int64_t _referenceEnd;
		bool _isQueryReversed;
	};

	struct AlignmentResult{
		ReadType _readIndex;
		float _score;
		Chain2 _chain;
	};

	struct QueryAlignmentResult{
		ReadType _queryReadIndex;
		u_int32_t _start;
		u_int32_t _end;
		float _score;
	};

	struct AlignmentMatches{
		ReadType _readIndex;
		u_int16_t _nbMatches;
		vector<unsigned char> _compressedMatchPositions;
	};


	//ankerl::unordered_dense::map<ReadType, ReferenceRead> _referenceReads;

	ReadMapper(const string& inputDir, const string& alignmentFilename, const size_t totalNbReads, const u_int64_t maxChainingBand, const u_int64_t usedCoverageForCorrection, const u_int64_t nbMinimizersPerChunk, const u_int64_t nbBases, const int nbCores){
		_inputDir = inputDir;
		_alignmentFilename = alignmentFilename;
		_totalNbReads = totalNbReads;
		_maxChainingBand = maxChainingBand;
		_usedCoverageForCorrection = usedCoverageForCorrection;
		_nbMinimizersPerChunk = nbMinimizersPerChunk;
		_nbCores = nbCores;

		_nbPartitions = nbBases / 20000000000ull;
		_nbPartitions = max(_nbPartitions, _nbCores);
		_nbPartitions = max(_nbPartitions, 1) ;//, 200);
		_nbPartitions = min(_nbPartitions, 5000);

		Logger::get().debug() << "\tNb alignment partitions: " << _nbPartitions << "\n";

		_partitionDir = inputDir + "/alignments_tmp/";
		if(!fs::exists(_partitionDir)){
			fs::create_directories(_partitionDir); 
		}
	}
	

	void execute(){

		//_debugFile = ofstream(_inputDir + "/lala_2.txt");


		auto start = high_resolution_clock::now();



		_nbReads = 0;
		_kminmerSize = 2;
		_alignmentCheckSum = 0;
		//_firstReadIndexInChunk = -1;
		_nbReadsInChunk = 0;
		_nbReadsProcessed = 0;


		_partitionFiles.resize(_nbPartitions);

		for(size_t i=0; i<_partitionFiles.size(); i++){
			_partitionFiles[i] = bgzf_open(getPartitionFilename(i).c_str(), "w1");  
		}

		//u_int64_t maxReads = 1000;
		//u_int64_t maxBps = 1000000000;
		//u_int64_t maxBps = 3000000000;
		//u_int64_t maxNbMinimizers = 15000000; //maxBps / (1/_minimizerDensity);
		//cout << "Max nb minimizers: " << maxNbMinimizers << endl;

		KminmerParserParallel parser(_inputDir + "/read_data_init.txt", 0, _kminmerSize, false, true, 1);
		parser._maxChunkSize = _nbMinimizersPerChunk;
		parser.parseChunk(ReadFunctor(*this), ChunkFunctor(*this));

		//cout << "Nb contigs: " << _nbContigs << endl;

		for(BGZF* outputFile : _partitionFiles){
			bgzf_close(outputFile);
		}


		_alignmentFile = ofstream(_alignmentFilename);

		mergeAlignmentPartitions();

		_alignmentFile.close();
		_debugFile.close();

		//cout << "\tAlignment checksum: " << _alignmentCheckSum << endl;
		Logger::get().debug() << "Nb reads processed: " << _nbReadsProcessed;
		Logger::get().debug() << "Alignment checksum: " << _alignmentCheckSum;
		Logger::get().debug() << "Done: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";

		if(!DEBUG) fs::remove_all(_partitionDir);
		
	}

	void mergeAlignmentPartitions(){

		Logger::get().debug() << "\tMerging alignment partitions";

		for(size_t i=0; i<_nbPartitions; i++){
			mergeAlignmentPartition(i);
		}

	}

	void mergeAlignmentPartition(const int partition){

		Logger::get().debug() << "\t\tMerging alignment partition: " << partition;

		const string partitionFilename = getPartitionFilename(partition);

		const string partitionFilenameDecompressed = partitionFilename + ".decomp";
		Utils::bgzip_decompress(partitionFilename, partitionFilenameDecompressed, _nbCores);

		ifstream partitionFile(partitionFilenameDecompressed);

		
		ankerl::unordered_dense::map<ReadType, vector<AlignmentScoreMatches>> referenceReadIndex_to_alignments;
		//ankerl::unordered_dense::map<ReadType, vector<AlignmentScore>> referenceReadIndex_to_size;


		while(true){

			ReadType referenceReadIndex;
			//u_int32_t referenceReadSize;
			u_int32_t nbAlignments;
			vector<AlignmentScoreMatches> alignments;

			partitionFile.read((char*)&referenceReadIndex, sizeof(referenceReadIndex));
			

			if(partitionFile.eof()) break;

			//partitionFile.read((char*)&referenceReadSize, sizeof(referenceReadSize));
			partitionFile.read((char*)&nbAlignments, sizeof(nbAlignments));
			

			for(size_t i=0; i<nbAlignments; i++){

				ReadType queryReadIndex;
				//int32_t score;
				u_int16_t nbMatches;
				u_int16_t compressedSize;
				vector<unsigned char> compressedPositionMatches;

				partitionFile.read((char*)&queryReadIndex, sizeof(queryReadIndex));
				//partitionFile.read((char*)&score, sizeof(score));

				partitionFile.read((char*)&nbMatches, sizeof(nbMatches));
				partitionFile.read((char*)&compressedSize, sizeof(compressedSize));
				compressedPositionMatches.resize(compressedSize);
				partitionFile.read((char*)&compressedPositionMatches[0], compressedPositionMatches.size());

				alignments.push_back({queryReadIndex, nbMatches, compressedPositionMatches});
			}



			/*
			partitionFile.read((char*)&alignments[0], alignments.size() * sizeof(AlignmentScore));

			//cout << referenceReadIndex << endl;
			//for(const AlignmentScore& alignment : alignments){
			//	cout << "\t" << alignment._queryReadIndex << endl;
			//}
			*/

			if(referenceReadIndex_to_alignments.find(referenceReadIndex) == referenceReadIndex_to_alignments.end()){
				referenceReadIndex_to_alignments[referenceReadIndex] = {};
				//cout << "init: " << referenceReadIndex << " " << referenceReadSize << " " << referenceReadIndex_to_queue[referenceReadIndex].size() << endl;
			}

			vector<AlignmentScoreMatches>& referenceAlignments = referenceReadIndex_to_alignments[referenceReadIndex];

			referenceAlignments.insert(referenceAlignments.end(), alignments.begin(), alignments.end());



			

        }

		partitionFile.close();
		fs::remove(partitionFilenameDecompressed);


		vector<ReadType> referenceReadIndexes;
		for(const auto& it : referenceReadIndex_to_alignments){
			referenceReadIndexes.push_back(it.first);
		}

		std::sort(referenceReadIndexes.begin(), referenceReadIndexes.end());

		#pragma omp parallel for num_threads(_nbCores)
		for(size_t i=0; i<referenceReadIndexes.size(); i++){

			const ReadType& referenceReadIndex = referenceReadIndexes[i];
			const vector<AlignmentScoreMatches>& alignments = referenceReadIndex_to_alignments[referenceReadIndex];

			//cout << referenceReadIndex << " " << alignments.size() << endl;
			vector<AlignmentScoreQueue> alignmentQueues;

			for(const AlignmentScoreMatches& alignment : alignments){
				mergeAlignmentScore(referenceReadIndex, alignment, alignmentQueues);
			}

			ankerl::unordered_dense::map<ReadType, AlignmentScore> selectedAlignments;
			//unordered_set<ReadType> selectedQueryReadIndex;
			vector<AlignmentResult2> bestAlignmentsLowDensity;
			
			for(AlignmentScoreQueue& queue : alignmentQueues){

				int i =0;
				while(queue.size() > 0){
					const AlignmentScore& alignmentScore = queue.top();
					selectedAlignments[alignmentScore._queryReadIndex] = alignmentScore;
					//selectedQueryReadIndex.insert(alignmentScore._queryReadIndex);
					queue.pop();
					i += 1;
				}

			}

			vector<AlignmentScore> selectedAlignmentsVec;

			for(const auto& it : selectedAlignments){
				selectedAlignmentsVec.push_back(it.second);
			}


			std::sort(selectedAlignmentsVec.begin(), selectedAlignmentsVec.end(), [](const AlignmentScore& a, const AlignmentScore& b){
				return a._queryReadIndex < b._queryReadIndex;
			});


			writeAlignments(referenceReadIndex, selectedAlignmentsVec);
		}


		Logger::get().debug() << "\t\tDone: " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";
	}

	void mergeAlignmentScore(const ReadType& referenceReadIndex, const AlignmentScoreMatches& currentAlignmentScore, vector<AlignmentScoreQueue>& alignmentQueues){
		

		vector<u_int32_t> queryMatchPositions(currentAlignmentScore._nbMatches, 0);
		p4nd1dec32((unsigned char*) currentAlignmentScore._compressedMatchPositions.data(), queryMatchPositions.size(), queryMatchPositions.data());

		


		//cout << "score: " << currentAlignmentScore._score << endl;

		int32_t score = 0;
		for(size_t i=0; i<((long)queryMatchPositions.size())-1; i++){
			score += 1;
			int32_t nbErrors = queryMatchPositions[i+1] - queryMatchPositions[i] - 1;
			score -= nbErrors;
		}
		score += 1;

		const AlignmentScore currentAlignmentScore2 = {currentAlignmentScore._queryReadIndex, score};
		//cout << "\tnew score: " << score << endl;

		//if(score != currentAlignmentScore._score) getchar();

		//if(referenceReadIndex == 0){
		//	cout << queryRead._readIndex << "\t" << chain._nbMatches << endl;
		//}

		//float score = chain._nbMatches - chain._nbDifferencesQuery;// - chain._nbDifferencesReference;
		//const AlignmentScore currentAlignmentScore = {referenceReadIndex, score, chain._queryStart, chain._queryEnd, chain._isQueryReversed};

		//if(queryRead._readIndex == 79 && referenceReadIndex == 92){
		//	cout << "lala: " << chain._nbMatches << " " << chain._nbDifferencesQuery << " " << score << endl;
		//}

		//if(referenceReadIndex == 0 && queryRead._readIndex == 157){
		//	cout << "lala: " << chain._nbMatches << " " << chain._nbDifferencesQuery << endl;
		//}


		//cout << chain._nbMatches << " " << chain._nbDifferencesQuery << endl;
		//cout << queryRead._readIndex << "\t" << referenceReadIndex << "\t" << chain._nbMatches << "\t" << (int32_t)score << endl;


		for(const u_int32_t& positionIndex : queryMatchPositions){
		
			while(positionIndex >= alignmentQueues.size()){
				alignmentQueues.push_back(AlignmentScoreQueue());
			}
			//cout << "\t" << i << "\t" << alignmentQueues.size() << endl;
			//const u_int32_t referencePosition = it.first;
			//if(referencePosition == -1) continue;
			//for(MinimizerType minimizer : queryRead._minimizers){
			//if(minimizer_to_alignmentScoreQueue.find(minimizer) == minimizer_to_alignmentScoreQueue.end()) continue;

			AlignmentScoreQueue& queue = alignmentQueues[positionIndex];

			if(queue.size() < _usedCoverageForCorrection){
				//if(referenceReadIndex == 5) cout << "Push:\ta\t" << i << "\t" << currentAlignmentScore._queryReadIndex << "\t" << currentAlignmentScore._score << endl;
				queue.push(currentAlignmentScore2);
				continue;
			}

			const AlignmentScore& worseAlignmentScore = queue.top();

			if(currentAlignmentScore2._score < worseAlignmentScore._score) continue;

			if(currentAlignmentScore2._score == worseAlignmentScore._score){
				if(currentAlignmentScore2._queryReadIndex > worseAlignmentScore._queryReadIndex) continue;
			}
			
			queue.pop();
			queue.push(currentAlignmentScore2);
			//if(referenceReadIndex == 5) cout << "Pop:\tb\t" << i << "\t" << worseAlignmentScore._queryReadIndex << endl;
			//if(referenceReadIndex == 5) cout << "Push:\tb\t" << i << "\t" << currentAlignmentScore._queryReadIndex << "\t" << currentAlignmentScore._score << endl;
		}
		

	}


	string getPartitionFilename(int partition){
		return _partitionDir + "/part_" + to_string(partition) + ".bin.gz";
	}

	class ReadFunctor {

		public:

		ReadMapper& _parent;

		ReadFunctor(ReadMapper& parent) : _parent(parent){
		}

		ReadFunctor(const ReadFunctor& copy) : _parent(copy._parent){
		}

		~ReadFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {

			ReadType readIndex = kminmerList._readIndex;
			if(readIndex % 100000 == 0) Logger::get().debug() << "\tLoad minimizers: " << readIndex;

			//bool isTooShort = Utils::isReadTooShort(kminmerList._readMinimizers.size());

			//cout << readIndex << " " << kminmerList._readMinimizers.size() << " " << kminmerList._kminmersInfo.size() << " " << _parent._allMinimizerPositions.size() << endl;

			//if(!isTooShort){
				for(u_int32_t i=0; i<kminmerList._kminmersInfo.size(); i++){
				
					
					//cout << "------- " << i << endl;

					//_logFile << readIndex << " " << i << endl;
					const ReadKminmerComplete& kminmerInfo = kminmerList._kminmersInfo[i];

					const KmerVec& vec = kminmerInfo._vec;

					//for(u_int32_t i=0; i<kminmerList._readMinimizers.size(); i++){
					
					//const MinimizerType minimizer = kminmerList._readMinimizers[i];
					//u_int32_t minimizerPosition = kminmerList._minimizerPos[i];
					//bool isReversed = kminmerList._readMinimizerDirections[i];

					const MinimizerPairType minimizerPair = vec.packPair();

					//const MinimizerPairType minimizerPair = vec.hash128();
					bool isReversed = kminmerInfo._isReversed;
					u_int32_t minimizerPosition = (kminmerList._minimizerPos[i] + kminmerList._minimizerPos[i+1]) / 2;

					_parent._allMinimizerPositions.push_back({minimizerPair, readIndex, minimizerPosition, i, isReversed});
					//_parent._allMinimizerPositions.push_back({minimizerPair, readIndex, minimizerPosition, isReversed});
				}
			//}

			_parent._nbReadsProcessed += 1;
			_parent._nbReadsInChunk += 1;

		}
	};

	class ChunkFunctor {

		public:

		ReadMapper& _parent;

		ChunkFunctor(ReadMapper& parent) : _parent(parent){
		}

		ChunkFunctor(const ChunkFunctor& copy) : _parent(copy._parent){
		}

		~ChunkFunctor(){
		}

		void operator () () const {


			auto start = high_resolution_clock::now();

			Logger::get().debug() << "\tProcessing chunk: " << _parent._nbReadsInChunk << " reads " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";
			//Logger::get().debug() << "\tDone: " << (peakrss() / 1024.0 / 1024.0 / 1024.0);

			_parent.processChunk();


			//Logger::get().debug() << "\tAlignment checksum: " << _parent._alignmentCheckSum;
			Logger::get().debug() << "\tDone: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB\n";
			//cout << "\tDone: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << endl;
			//phmap::flat_hash_map<MinimizerPairType, pair<u_int64_t, u_int64_t>>().swap(_parent._minimizerLookupTable);
			ankerl::unordered_dense::map<MinimizerPairType, pair<u_int64_t, u_int64_t>>().swap(_parent._minimizerLookupTable);
			vector<MinimizerPosition2>().swap(_parent._allMinimizerPositions);
			//ankerl::unordered_dense::map<ReadType, ReferenceRead>().swap(_parent._referenceReads);
			_parent._nbReadsInChunk = 0;

		}

	};


	//vector<KminmerList> _contigs;

	void processChunk(){

		//Logger::get().debug() << "Load minimizers";
		//loadReadMinimizers();
		//Logger::get().debug() << "\tDone: " << (peakrss() / 1024.0 / 1024.0 / 1024.0);
		
		Logger::get().debug() << "\t\tNb minimizers: " << _allMinimizerPositions.size();

		Logger::get().debug() << "\t\tSorting minimizers";
		Commons::sortParallel(_allMinimizerPositions, _allMinimizerPositions.size(), _nbCores);
		//sortReadMinimizers();
		Logger::get().debug() << "\t\tDone: " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";
		
		Logger::get().debug() << "\t\tCreating minimizer lookup table";
		createLookupTable();
		Logger::get().debug() << "\t\tDone: " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";

		//_file_readData = ofstream(_inputDir + "/read_data_corrected.txt");

		Logger::get().debug() << "\t\tPerform low density alignments";
		alignReadsLowDensity();
		
	}


	void createLookupTable(){
		u_int64_t index = -1;
		MinimizerPairType minimizer = -1;
		MinimizerPairType lastInsertedMinimizer = -1;

		for(size_t i=0; i<_allMinimizerPositions.size(); i++){

			if(minimizer != _allMinimizerPositions[i]._minimizerPair){

				if(index != -1){

					u_int64_t startIndex = index;
					u_int64_t endIndex = i-1;

					//if(startIndex != endIndex){
						_minimizerLookupTable[minimizer] = {startIndex, endIndex};
					//}

					lastInsertedMinimizer = minimizer;
					//cout << "Found chunk: " << startIndex << " -> " << endIndex << endl;
				}

				minimizer = _allMinimizerPositions[i]._minimizerPair;
				index = i;
			}
		}

		if(lastInsertedMinimizer != minimizer){
			u_int64_t startIndex = index;
			u_int64_t endIndex = _allMinimizerPositions.size()-1;

			//if(startIndex != endIndex){
				_minimizerLookupTable[minimizer] = {startIndex, endIndex};
			//}
			//cout << "Found chunk last: " << startIndex << " -> " << endIndex << endl;
		}

		/*
		u_int64_t nbPositions = 0;

		for(const auto& it : _minimizerLookupTable){
			for(size_t i=it.second.first; i <= it.second.second; i++){
				nbPositions += 1;
			}
		}

		cout << "Lookup table check: " << _allMinimizerPositions.size() << " " << nbPositions << endl;
		if(nbPositions != _allMinimizerPositions.size()){
			cout << "Lookup table issue" << endl;
			getchar();
		}
		*/
	}




	void alignReadsLowDensity(){


		//cout << "single core here" << endl;
		KminmerParserParallel parser(_inputDir + "/read_data_init.txt", 0, _kminmerSize, false, true, _nbCores);
		//parser._densityThreshold = _minimizerDensity_assembly;
		parser.parse(AlignReadsLowDensityFunctor(*this));

	}


	class AlignReadsLowDensityFunctor {

		public:

		ReadMapper& _parent;
		//MinimizerAligner* _minimizerAligner;
		//MinimizerChainer* _minimizerChainer;
		
		
		AlignReadsLowDensityFunctor(ReadMapper& parent) : _parent(parent){
			//_minimizerAligner = nullptr;
			//_minimizerChainer = nullptr;
		}

		AlignReadsLowDensityFunctor(const AlignReadsLowDensityFunctor& copy) : _parent(copy._parent){
			//_minimizerAligner = new MinimizerAligner(3, -1, -1, -1);
			//_minimizerChainer = new MinimizerChainer(_parent._minimizerSize);
		}

		~AlignReadsLowDensityFunctor(){
			//if(_minimizerAligner != nullptr) delete _minimizerAligner;
			//if(_minimizerChainer != nullptr) delete _minimizerChainer;
		}


		void operator () (const KminmerList& kminmerList) {

			//cout << "Align: " << kminmerList._readIndex << endl;


			if(kminmerList._readIndex % 100000 == 0)  Logger::get().debug() << "\t\tAligning reads: " << Utils::getProgress(kminmerList._readIndex, _parent._totalNbReads);

			MinimizerRead read = {kminmerList._readIndex, kminmerList._readMinimizers, kminmerList._minimizerPos, kminmerList._readQualities, kminmerList._readMinimizerDirections, kminmerList._readLength};
			//MinimizerRead readLowDensity = Utils::getLowDensityMinimizerRead(referenceRead, _parent._minimizerDensity_assembly);
			
			//if(readLowDensity._minimizers.size() < _parent._minReadLength){
			if(Utils::isReadTooShort(read)){
				//_parent.writeRead(referenceRead._readIndex, referenceRead._minimizers, referenceRead._qualities);
				//_parent.writeAlignments(read._readIndex, {});
				return;
			}





			ReadType readIndex = read._readIndex;

			//if(readIndex % 10000 == 0) Logger::get().debug() << readIndex;
			//cout << readIndex << endl;

			vector<MinimizerType> readMinimizers = read._minimizers;

			vector<Anchor2> anchors;

			for(u_int32_t i=0; i<kminmerList._kminmersInfo.size(); i++){
			
				
				//cout << "------- " << i << endl;

				//_logFile << readIndex << " " << i << endl;
				const ReadKminmerComplete& kminmerInfo = kminmerList._kminmersInfo[i];

				const KmerVec& vec = kminmerInfo._vec;

				//for(u_int32_t i=0; i<kminmerList._readMinimizers.size(); i++){
				
				//const MinimizerType minimizer = kminmerList._readMinimizers[i];
				//u_int32_t minimizerPosition = kminmerList._minimizerPos[i];
				//bool isReversed = kminmerList._readMinimizerDirections[i];

				//const MinimizerPairType minimizerPair = vec.hash128();
				const MinimizerPairType minimizerPair = vec.packPair();
				bool queryIsReversed = kminmerInfo._isReversed;
				u_int32_t queryPosition = (kminmerList._minimizerPos[i] + kminmerList._minimizerPos[i+1]) / 2;

			//for(u_int32_t i=0; i<read._minimizers.size(); i++){


				//const MinimizerType minimizer = read._minimizers[i];
				//const u_int32_t queryPosition = read._minimizersPos[i];
				//const bool queryIsReversed = read._readMinimizerDirections[i];

				if(_parent._minimizerLookupTable.find(minimizerPair) == _parent._minimizerLookupTable.end()) continue;

				const auto& range = _parent._minimizerLookupTable[minimizerPair];

				for(size_t j=range.first; j<=range.second; j++){
					const MinimizerPosition2& mPos = _parent._allMinimizerPositions[j];

					if(kminmerList._readIndex == mPos._readIndex) continue; //query == reference

					//Anchor2(const ReadType readIndex, const int64_t referencePosition, const int64_t queryPosition, const bool isReversed){
					//anchors.push_back({referenceMinimizerPosition._position, queryPosition, referenceMinimizerPosition._isReversed != queryIsReversed, referenceMinimizerPosition._positionIndex, i});
					//anchors.push_back({mPos._readIndex, mPos._position, queryPosition, i, mPos._isReversed != queryIsReversed});
					anchors.push_back({mPos._readIndex, mPos._position, mPos._positionIndex, queryPosition, i, mPos._isReversed != queryIsReversed});
					//_loul += 1;

					//if(mPos._position == 2600) print = true;
				}
			}

			std::sort(anchors.begin(), anchors.end(), [](const Anchor2& a, const Anchor2& b){
				if(a._readIndex == b._readIndex){

					if(a._referencePosition == b._referencePosition){
						return a._queryPosition < b._queryPosition;
					}
					return a._referencePosition < b._referencePosition;
				}
				
				return a._readIndex < b._readIndex;

			});

			/*
			if(print){
				cout << "Read: " << readIndex << endl;
				//unordered_set<u_int32_t> contigIndexes;


				for(size_t i=0; i<anchors.size(); i++){
					const Anchor& anchor = anchors[i];
					cout << i << "\t" << anchor._referenceIndex << "\t" << anchor._referencePosition << "\t" << anchor._queryPosition << "\t" << read._minimizersPos[anchor._queryPosition] << endl;
					//contigIndexes.insert(anchor._referenceIndex);
				}
			}
			*/

			vector<AlignmentMatches> alignmentMatches;
			vector<AlignmentScoreQueue> alignmentQueues(read._minimizers.size());

			u_int32_t nbChainings = 0;
			u_int32_t referenceIndex = -1;
			vector<Anchor2> subAnchors;

			for(size_t i=0; i<anchors.size(); i++){

				if(referenceIndex != anchors[i]._readIndex){

					if(referenceIndex != -1){
						processAnchors(subAnchors, read, referenceIndex, alignmentQueues, alignmentMatches);
						subAnchors.clear();
						nbChainings += 1;
					}

					referenceIndex = anchors[i]._readIndex;
				}

				subAnchors.push_back(anchors[i]);


			}

			if(subAnchors.size() > 0){
				processAnchors(subAnchors, read, referenceIndex, alignmentQueues, alignmentMatches);
				nbChainings += 1;
			}



			
			ankerl::unordered_dense::map<ReadType, AlignmentScoreMatches> selectedAlignments;
			//unordered_set<ReadType> selectedQueryReadIndex;
			//vector<AlignmentResult2> bestAlignmentsLowDensity;
			
			for(AlignmentScoreQueue& queue : alignmentQueues){

				int i =0;
				while(queue.size() > 0){
					const AlignmentScore& alignmentScore = queue.top();
					selectedAlignments[alignmentScore._queryReadIndex] = {alignmentScore._queryReadIndex, 0, {}};
					//selectedQueryReadIndex.insert(alignmentScore._queryReadIndex);
					queue.pop();
					i += 1;
				}

			}

			for(const AlignmentMatches& matches : alignmentMatches){
				if(selectedAlignments.find(matches._readIndex) == selectedAlignments.end()) continue;

				selectedAlignments[matches._readIndex]._nbMatches = matches._nbMatches;
				selectedAlignments[matches._readIndex]._compressedMatchPositions = matches._compressedMatchPositions;
			}

			vector<AlignmentScoreMatches> selectedAlignmentsVec;

			for(const auto& it : selectedAlignments){
				selectedAlignmentsVec.push_back(it.second);
			}


			std::sort(selectedAlignmentsVec.begin(), selectedAlignmentsVec.end(), [](const AlignmentScoreMatches& a, const AlignmentScoreMatches& b){
				return a._queryReadIndex < b._queryReadIndex;
			});


			_parent.writeAlignmentsPartition(read._readIndex, read._minimizers.size(), selectedAlignmentsVec);


			
		}

		
		void processAnchors(const vector<Anchor2>& anchors, const MinimizerRead& queryRead, const u_int32_t referenceReadIndex, vector<AlignmentScoreQueue>& alignmentQueues, vector<AlignmentMatches>& alignmentMatches){
			
			if(anchors.size() < 3) return;

			//cout << "perform align" << endl;
			/*
			if(contigIndex == 4735 && queryRead._readIndex == 576){
				cout << endl;
				for(size_t i=0; i<anchors.size(); i++){
					const Anchor2& anchor = anchors[i];
					cout << i << "\t" << anchor._readIndex << "\t" << anchor._referencePosition << "\t" << anchor._queryPosition << "\t" << anchor._queryPositionIndex << endl;

				}
			}
			*/
			
			

			Chain2 chain = chainAnchors(anchors, queryRead, _parent._maxChainingBand, referenceReadIndex);

			/*
			//bool print = queryRead._readIndex == 2479658;
			//if(print){
			cout << endl;
			for(size_t i=0; i<chain._anchorInterval.size(); i++){
				cout << chain._anchorInterval[i] << endl;
			}
			//}
			*/

			if(chain._chainingScore == 0) return;

			//cout << chain._chainingScore << " " << chain._queryAnchorPositions.size() << endl;
			addAlignmentScore(chain, queryRead, referenceReadIndex, alignmentQueues, alignmentMatches);
			//writeAlignment(chain, queryRead, contigIndex, isReversed);
			//getchar();
		}


		Chain2 chainAnchors(const vector<Anchor2>& anchors, const MinimizerRead& queryRead, const int64_t& maxChainingBand, const ReadType referenceReadIndex) {

			bool print = false;//referenceReadIndex == 4735 && queryRead._readIndex == 576;
			//size_t nbAnchors = _anchorIndex;
			vector<Chain> chainIntervals;

			//cout << "---- " << points.size() << endl;
			float w = 20;//_minimizerSize;// 20;// _minimizerSize; //_minimizerSize*5; //_parent._minimizerSize;
			vector<float> scores(anchors.size(), 0);
			vector<size_t> parents(anchors.size(), 0);

			//if(nbAnchors > _scores.size()){
			//	_scores.resize(nbAnchors);
			//	_parents.resize(nbAnchors);
			//}

			//for(size_t i=0; i<nbAnchors; i++){
			//	_scores[i] = 0;
			//	_parents[i] = 0;
			//}

			//for (var i in global_list_of_point) {
			for (size_t i=0; i<anchors.size(); i++) {

				//if(_parent._print_debug) cout << "\t--- " << i << endl;
				
				//cout << i << " " << points.size() << endl;
				//size_t j = i;
				//var j = parseInt(i)

				float bestScore = 0;
				size_t bestPrevIndex = i;

				//cout << "a" << endl;
				argmaxPosition(anchors, scores, i, w, bestScore, bestPrevIndex, maxChainingBand, false);

				
				if(bestPrevIndex != i) {
					scores[i] = bestScore;
					parents[i] = bestPrevIndex;
					//cout << "Add chain part: " << i << " " << bestPrevIndex << endl;
					//chain_part.union(i, best_prev_index);
				}
				else{
					scores[i] = w;
					parents[i] = -1;
				}
				/*
				if(bestScore < w){
					scores[i] = w;
					parents[i] = -1;
				}
				else{
					scores[i] = bestScore;
					parents[i] = bestPrevIndex;
				}
				*/

			}


			//for(size_t i=0; i<scores.size(); i++) {
			//	cout << i << "\t" << scores[i] << "\t" << parents[i] << endl;
			//}

			//if(_print_debug){
			//	cout << endl;
			//	for(size_t i=0; i<scores.size(); i++) {
			//		cout << "\tlala: " << i << " " << scores[i] << " " << getRoot(i, parents) << endl;
			//		F.push_back({i, scores[i], getRoot(i, parents)});
			//	}
			//}
			
			float maxScore = 0;
			size_t bestIndex = -1;

			for(size_t i=0; i<anchors.size(); i++) {
				if(scores[i] > maxScore){
					maxScore = scores[i];
					bestIndex = i;
				}
				//cout << "\tlala: " << i << " " << scores[i] << " " << getRoot(i, parents) << endl;
				//F.push_back({i, scores[i], getRoot(i, parents)});
			}

			if(print){
				cout << "\tParents:" << endl;
				for(size_t i=0; i<anchors.size(); i++){
					cout << "\t\t" << i << "\t" << anchors[i]._isReversed << "\t" << parents[i] << "\t" << scores[i] << endl;
				}
			}

			//cout << bestIndex << endl;

			vector<size_t> interval;

			size_t nullVal = -1;
			size_t index = bestIndex;
			while (index != nullVal) {
				interval.push_back(index);
				index = parents[index];
				//num_anchors += 1;
			}

			//cout << bestIndex << " " << index << endl;
			if(print){
				cout << "\tBest interval:" << endl;
				for(size_t i=0; i<interval.size(); i++){
					cout << "\t\t" << interval[i] << endl;
				}
			}

			std::reverse(interval.begin(), interval.end());
			chainIntervals.push_back({maxScore, interval});
			//std::sort(F.begin(), F.end(), [](const Point2& a, const Point2& b){
			//	if(a._score == b._score){
			//		return a._fromIndex > b._fromIndex;
			//	}
			//	return a._score > b._score;
			//});



			Chain2 chain;
			chain._chainingScore = 0;

			if(interval.size() < 3) return chain;

			std::reverse(interval.begin(), interval.end());

			for(size_t i=0; i<interval.size(); i++){
				chain._queryAnchorPositions.push_back(anchors[interval[i]]._queryPositionIndex);
			}

			chain._chainingScore = maxScore;
			//chain._anchorInterval = interval;
			chain._nbMatches = interval.size();

			const Anchor2& firstAnchor = anchors[interval[0]];
			const Anchor2& lastAnchor = anchors[interval[interval.size()-1]];

			bool isQueryReversed = firstAnchor._queryPositionIndex > lastAnchor._queryPositionIndex;

			if(isQueryReversed){
				chain._nbDifferencesQuery = (firstAnchor._queryPositionIndex - lastAnchor._queryPositionIndex + 1) - chain._nbMatches;
				chain._queryStart = lastAnchor._queryPositionIndex;
				chain._queryEnd = firstAnchor._queryPositionIndex;
				std::reverse(chain._queryAnchorPositions.begin(), chain._queryAnchorPositions.end());
			}
			else{
				chain._nbDifferencesQuery = (lastAnchor._queryPositionIndex - firstAnchor._queryPositionIndex + 1) - chain._nbMatches;
				chain._queryStart = firstAnchor._queryPositionIndex;
				chain._queryEnd = lastAnchor._queryPositionIndex;
			}
			//chain._nbDifferencesReference = 0;

			//cout << firstAnchor._referencePositionIndex << " " << lastAnchor._referencePositionIndex << endl;
			chain._nbDifferencesReference = (firstAnchor._referencePositionIndex - lastAnchor._referencePositionIndex + 1) - chain._nbMatches;
			chain._referenceStart = firstAnchor._referencePosition;
			chain._referenceEnd = lastAnchor._referencePosition;
			chain._isQueryReversed = isQueryReversed;
			//cout << contigIndex << "\t" <<firstAnchor._queryPositionIndex << "\t" << lastAnchor._queryPositionIndex << "\t" << chain._nbMatches << "\t" << chain._nbDifferencesQuery << endl;
			


			/*
			chain._nbMatches = 0;
			chain._nbDifferencesQuery = 0;

			for(size_t i=0; i<interval.size()-1; i++){

				const Anchor2& currentAnchor = anchors[interval[i]];
				const Anchor2& nextAnchor = anchors[interval[i+1]];

				//int referenceGapIndexSize = nextAnchor._referencePositionIndex - currentAnchor._referencePositionIndex - 1;
				int queryGapIndexSize = 0;//
				
				if(isQueryReversed){
					queryGapIndexSize = currentAnchor._queryPositionIndex - nextAnchor._queryPositionIndex - 1;
				}
				else{
					queryGapIndexSize = nextAnchor._queryPositionIndex - currentAnchor._queryPositionIndex - 1;
				}

				if(print) cout << i << "\t" << interval[i] << "\t" << queryGapIndexSize << endl;

				chain._nbDifferencesQuery += queryGapIndexSize;

				chain._nbMatches += 1;
			}

			chain._nbMatches += 1;
			*/


			if(print){
				cout << chain._queryStart << " " << chain._queryEnd << " " << chain._nbMatches << " " << chain._nbDifferencesQuery << " " << chain._nbDifferencesReference << endl;
			}

			return chain; 
		}

		int getChainMissmatches(const size_t& bestIndex, const vector<size_t>& parents, const vector<Anchor2>& anchors){

			vector<size_t> interval;

			size_t nullVal = -1;
			size_t index = bestIndex;
			while (index != nullVal) {
				interval.push_back(index);
				index = parents[index];
				//num_anchors += 1;
			}

			std::reverse(interval.begin(), interval.end());
			//chainIntervals.push_back({maxScore, interval});


			//if(interval.size() < 3) return chain;

			std::reverse(interval.begin(), interval.end());

			const Anchor2& firstAnchor = anchors[interval[0]];
			const Anchor2& lastAnchor = anchors[interval[interval.size()-1]];

			bool isQueryReversed = firstAnchor._queryPositionIndex > lastAnchor._queryPositionIndex;



			if(isQueryReversed){
				//chain._nbDifferencesReference = (firstAnchor._referencePositionIndex - lastAnchor._referencePositionIndex + 1) - chain._nbMatches;
				return (firstAnchor._queryPositionIndex - lastAnchor._queryPositionIndex + 1);
			}
			else{
				//chain._nbDifferencesReference= (lastAnchor._referencePositionIndex - firstAnchor._referencePositionIndex + 1) - chain._nbMatches;
				return (lastAnchor._queryPositionIndex - firstAnchor._queryPositionIndex + 1);
			}


		}

		//size_t _maxChainingBand;

		void argmaxPosition(const vector<Anchor2>& anchors, vector<float>& scores, size_t i, float w, float& bestScore, size_t& bestPrevIndex, const int64_t& maxChainingBand, const bool print) {
			
			//cout << "---" << i << endl;
			bestScore = 0;
			bestPrevIndex = i;

			//cout << "\tglurp" << endl;
			//size_t bestJ = 0;
			//float v = 0;
			const Anchor2& xi = anchors[i];
			//console.log("\t", xi);
			for (int64_t j = i-1; j >= 0; j--) {
			//for (int64_t k =0; k < i; k++) {

				//cout << "\t\t" << j << endl;

				//cout << "\t\t" << i << " " << k << " " << i-k << " " << maxBand << endl;
				if(i-j > maxChainingBand) break;
				const Anchor2& xj = anchors[j];


				//if(xi._queryPosition - xj._queryPosition > 2500) continue;
				if(xi._isReversed != xj._isReversed) continue;
				if(xi._referencePosition == xj._referencePosition || xi._queryPosition == xj._queryPosition) continue;
				//cout << "\t" << i << " " << k << " " << xi._isReversed << " " << xj._isReversed << " " << xi._referencePosition << "-" << xi._queryPosition << " " << xj._referencePosition << "-" << xj._queryPosition << endl;




				int32_t d_q;// = abs(xi._queryPosition - xj._queryPosition);
				if(xi._isReversed){
					d_q = xj._queryPosition - xi._queryPosition;
				}
				else{
					d_q = xi._queryPosition - xj._queryPosition;
				}
				int32_t d_r = xi._referencePosition - xj._referencePosition;

				//cout << "\t" << xi._referencePosition << " " << xj._referencePosition << " " << d_r << endl;
				//if(d_r > 40) continue;
				//if(d_q > 40) continue;

				//if(xi._isReversed) {
				//	d_r = aprpf64 - acrpf64;
				//} else {
				//	d_r = ;
				//}
				//if(_parent._print_debug) cout << "\t\t" << j << "\t" << d_q << "\t" << d_r << "\t" << abs(d_r - d_q) << "\t" << referenceRead._minimizers[xi._referencePositionIndex] << "\t" << referenceRead._minimizers[xj._referencePositionIndex] << endl;
				//cout << "\t" << d_r << " " << d_q << endl;

				//if d_q > D_MAX_LIN_LENGTH || d_r > D_MAX_LIN_LENGTH {
				if(d_q > 5000 || d_r > 5000) {
					continue;
				}

				if(d_r <= 0) {
					continue;
				}

				int32_t gap = abs(d_r - d_q);
				//cout << "\t" << xi._referencePosition << "\t" << xj._referencePosition << "\t" << gap << endl;

				if(gap > 100) {
					continue;
				}
				
				if(xi._isReversed){
					if(xi._queryPosition > xj._queryPosition) continue;
				}
				else{
					if(xi._queryPosition < xj._queryPosition) continue;
				}
				
				//cout << w << " " << min(smallerDistance(xi,xj,xi._isReversed),w) << " " << gamma(xi,xj,w,a,b) << endl;
				//float anchor_score = min(smallerDistance(xi,xj,xi._isReversed),w) - gamma(xi,xj,w,a,b,xi._isReversed); //w - gap;
				float anchor_score = w - gap;

				float new_score = scores[j] + anchor_score;
				if (new_score > bestScore) {
					bestScore = new_score;
					bestPrevIndex = j;
					//if(_parent._print_debug) cout << "\t\tBest score:\t" << bestScore << "\t" << anchor_score << "\t" << j << "\t" << referenceRead._minimizers[xj._referencePositionIndex] << endl;
				}

				//map_params.anchor_score - gap


				//if (xi._referencePosition-xj._referencePosition < maxGapLength && xi._queryPosition-xj._queryPosition < maxGapLength && xi._queryPosition > xj._queryPosition) {
				//float newv = F[j] + 20 - gap;//+ 20 - gapLenth(xi,xj);//F[k] + min(smallerDistance(xi,xj),w) - gamma(xi,xj,w,a,b);
				//if (v < newv) {
				//	v = newv;
				//	bestJ = j;
				//}
				//}
			}
			//cout << "\tBest: " << bestPrevIndex << " " << bestScore << endl;
			//cout << "\tglurpppp" << endl;
			//jReturn = bestJ;
			//vReturn = v;
			
		}


		void addAlignmentScore(Chain2& chain, const MinimizerRead& queryRead, const ReadType referenceReadIndex, vector<AlignmentScoreQueue>& alignmentQueues, vector<AlignmentMatches>& alignmentMatches){

			//if(referenceReadIndex == 0){
			//	cout << queryRead._readIndex << "\t" << chain._nbMatches << endl;
			//}

			//cout << chain._nbMatches << " " << chain._nbDifferencesQuery << " " << chain._nbDifferencesReference << endl;
			bool isAlignmentAdded = false;
			int32_t score = chain._nbMatches - chain._nbDifferencesQuery; // - chain._nbDifferencesReference;// - chain._nbDifferencesReference;

			//if(score <= 0) return; 

			const AlignmentScore currentAlignmentScore = {referenceReadIndex, score};

			//if(queryRead._readIndex == 79 && referenceReadIndex == 92){
			//	cout << "lala: " << chain._nbMatches << " " << chain._nbDifferencesQuery << " " << score << endl;
			//}

			//if(referenceReadIndex == 0 && queryRead._readIndex == 157){
			//	cout << "lala: " << chain._nbMatches << " " << chain._nbDifferencesQuery << endl;
			//}


			//cout << chain._nbMatches << " " << chain._nbDifferencesQuery << endl;
			//cout << queryRead._readIndex << "\t" << referenceReadIndex << "\t" << chain._nbMatches << "\t" << (int32_t)score << endl;

			for(const u_int32_t& positionIndex : chain._queryAnchorPositions){
				
				//const u_int32_t referencePosition = it.first;
				//if(referencePosition == -1) continue;
				//for(MinimizerType minimizer : queryRead._minimizers){
				//if(minimizer_to_alignmentScoreQueue.find(minimizer) == minimizer_to_alignmentScoreQueue.end()) continue;

				AlignmentScoreQueue& queue = alignmentQueues[positionIndex];

				if(queue.size() < _parent._usedCoverageForCorrection){

					if(!isAlignmentAdded){


						//#pragma omp critical
						//{
						//	cout << referenceReadIndex << "\t" << compressMatches(chain._queryAnchorPositions).size() << endl;
						//}

						isAlignmentAdded = true;
						alignmentMatches.push_back({referenceReadIndex, (u_int16_t)chain._queryAnchorPositions.size(), compressMatches(chain._queryAnchorPositions)});
					}

					//if(queryRead._readIndex == 5) cout << "Push:\ta\t" << positionIndex << "\t" << currentAlignmentScore._queryReadIndex << "\t" << currentAlignmentScore._score << endl;
					queue.push(currentAlignmentScore);
					continue;
				}

				const AlignmentScore& worseAlignmentScore = queue.top();

				if(currentAlignmentScore._score < worseAlignmentScore._score) continue;

				if(currentAlignmentScore._score == worseAlignmentScore._score){
					if(currentAlignmentScore._queryReadIndex > worseAlignmentScore._queryReadIndex) continue;
				}
				
				queue.pop();
				queue.push(currentAlignmentScore);

				if(!isAlignmentAdded){
					isAlignmentAdded = true;
					alignmentMatches.push_back({referenceReadIndex, (u_int16_t)chain._queryAnchorPositions.size(), compressMatches(chain._queryAnchorPositions)});

					//#pragma omp critical
					//{
					//	cout << referenceReadIndex << "\t" << compressMatches(chain._queryAnchorPositions).size() << endl;
					//}
				}

				//if(queryRead._readIndex == 5) cout << "Pop:\tb\t" << positionIndex << "\t" << worseAlignmentScore._queryReadIndex << endl;
				//if(queryRead._readIndex == 5) cout << "Push:\tb\t" << positionIndex << "\t" << currentAlignmentScore._queryReadIndex << "\t" << currentAlignmentScore._score << endl;
			}
				

		}

		vector<unsigned char> compressMatches(vector<u_int32_t>& matchPositions){
			vector<unsigned char> buffer(matchPositions.size()*32, 0);
			u_int32_t compressed_vector_size = p4nd1enc32(matchPositions.data(), matchPositions.size() , buffer.data());
			buffer.resize(compressed_vector_size);
			return buffer;
		}

	};


    struct AlignmentWriter{
        ReadType _referenceReadIndex;
        vector<ReadType> _alignedQueryReads;
    };

    struct AlignmentWriter_Comparator {
        bool operator()(const AlignmentWriter& p1, const AlignmentWriter& p2){
            return p1._referenceReadIndex > p2._referenceReadIndex;
        }
    };

	//priority_queue<AlignmentWriter, vector<AlignmentWriter> , AlignmentWriter_Comparator> _alignmentWriterQueue;
	//ReadType _nextAlignmentReadIndexWriter;

	void writeAlignmentsPartition(const ReadType referenceReadIndex, const u_int32_t referenceReadSize, const vector<AlignmentScoreMatches>& alignments){
		
		if(alignments.size() == 0) return;
		
		int partition = referenceReadIndex % _nbPartitions;
		BGZF* partitionFile = _partitionFiles[partition];

		#pragma omp critical(writeAlignments)
		{
			
			const u_int32_t nbAlignedReads = alignments.size();

			//vector<ReadType> alignedQueryReads;

			//for(const AlignmentScore& alignment : alignments){
				//_alignmentCheckSum += referenceReadIndex * nbAlignedReads * alignment._queryReadIndex;
				//alignedQueryReads.push_back(alignment._queryReadIndex);
			//}

			bgzf_write(partitionFile, (const char*)&referenceReadIndex, sizeof(referenceReadIndex));
			//partitionFile.write((const char*)&referenceReadSize, sizeof(referenceReadSize));
			bgzf_write(partitionFile, (const char*)&nbAlignedReads, sizeof(nbAlignedReads));

			for(const AlignmentScoreMatches& al : alignments){

				bgzf_write(partitionFile, (const char*)&al._queryReadIndex, sizeof(al._queryReadIndex));
				//bgzf_write(partitionFile, (const char*)&al._score, sizeof(al._score));

				u_int16_t nbMatches = al._nbMatches;
				bgzf_write(partitionFile, (const char*)&nbMatches, sizeof(nbMatches));

				u_int16_t compressedSize = al._compressedMatchPositions.size();
				bgzf_write(partitionFile, (const char*)&compressedSize, sizeof(compressedSize));

				bgzf_write(partitionFile, (const char*)&al._compressedMatchPositions[0], compressedSize);

			}
			

			//if(referenceReadIndex % 10000 == 0){
			//	Logger::get().debug() << "\tAlignment checksum: " << referenceReadIndex << " " << _alignmentCheckSum;
				//if(readWriter._referenceReadIndex > 1000) getchar();
			//}



		}


	}


	void writeAlignments(const ReadType referenceReadIndex, const vector<AlignmentScore>& alignments){
		
		if(alignments.size() == 0) return;

		#pragma omp critical(writeAlignments)
		{
			
			const u_int32_t nbAlignedReads = alignments.size();

			vector<ReadType> alignedQueryReads;

			for(const AlignmentScore& alignment : alignments){
				_alignmentCheckSum += referenceReadIndex * nbAlignedReads * alignment._queryReadIndex;
				alignedQueryReads.push_back(alignment._queryReadIndex);
			}


			_alignmentFile.write((const char*)&referenceReadIndex, sizeof(referenceReadIndex));
			_alignmentFile.write((const char*)&nbAlignedReads, sizeof(nbAlignedReads));
			_alignmentFile.write((const char*)&alignedQueryReads[0], nbAlignedReads*sizeof(ReadType));
			

			//if(referenceReadIndex % 10000 == 0){
			//	Logger::get().debug() << "\tAlignment checksum: " << referenceReadIndex << " " << _alignmentCheckSum;
				//if(readWriter._referenceReadIndex > 1000) getchar();
			//}

			//_debugFile << referenceReadIndex << endl;
			//for(const AlignmentScore& alignment : alignments){
			//	_debugFile << "\t" << alignment._queryReadIndex << endl;
			//}

		}


	}

};

#endif