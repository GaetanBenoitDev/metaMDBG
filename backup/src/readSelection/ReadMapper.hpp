
#ifndef MDBG_METAG_READMAPPER
#define MDBG_METAG_READMAPPER

#include "../Commons.hpp"
#include "../readSelection/MinimizerChainer.hpp"


class ReadMapper{
    
public:

	class ReadWorkerLock{

		public:

		vector<omp_lock_t> _locks;

		ReadWorkerLock(){

		}

		void init(int nbLocks){
			_locks.resize(nbLocks);

			for(size_t i=0; i<_locks.size(); i++){
				omp_init_lock(&_locks[i]);
			}
			
		
		}

		~ReadWorkerLock(){
			for(size_t i=0; i<_locks.size(); i++){
				omp_destroy_lock(&_locks[i]);
			}
		}

		void lock(const ReadType& readIndex){

			int lockIndex = readIndex % _locks.size();

			omp_set_lock(&_locks[lockIndex]);
		
		}

		void unlock(const ReadType& readIndex){

			int lockIndex = readIndex % _locks.size();

			omp_unset_lock(&_locks[lockIndex]);
		}

	};

	struct AlignmentScore{
		ReadType _queryReadIndex;
		float _score;
		float _identity;
	};

	struct MinimizerPairPosition{
		ReadType _readIndex;
		u_int16_t _positionIndex;
	};

	struct MinimizerPosition{
		int32_t _position;
		int32_t _positionIndex;
		bool _isReversed;
	};

	typedef phmap::parallel_flat_hash_map<KmerVec, vector<MinimizerPairPosition>, phmap::priv::hash_default_hash<KmerVec>, phmap::priv::hash_default_eq<KmerVec>, std::allocator<std::pair<KmerVec, vector<MinimizerPairPosition>>>, 4, std::mutex> KminmerPosMap;
	typedef unordered_map<MinimizerType, vector<MinimizerPosition>> ReadMinimizerPositionMap;

	class ReadWorker{


		public:

		struct AlignmentScore{
			ReadType _queryReadIndex;
			float _score;
			//float _identity;
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


		MinimizerRead _referenceRead;
		unordered_map<MinimizerType, AlignmentScoreQueue> _minimizer_to_alignmentScoreQueue;
		//vector<ReadType> _mappedReads;
		//unordered_map<MinimizerType, vector<AlignmentScore>> _minimizer_to_alignmentResultsNew;


		ReadWorker(){
		}

		ReadWorker(const MinimizerRead& referenceRead){
			_referenceRead = referenceRead;
			for(const MinimizerType& minimizer : referenceRead._minimizers){
				_minimizer_to_alignmentScoreQueue[minimizer] = {};
			}
		}

		void addAlignmentScore(const MinimizerRead& queryRead, const AlignmentResult2& alignment, ReadWorkerLock& readWorkerLock){

			readWorkerLock.lock(_referenceRead._readIndex);


			float score = alignment._nbMatches - alignment._nbMissmatches - alignment._nbInsertions - alignment._nbDeletions;
			const AlignmentScore currentAlignmentScore = {alignment._queryReadIndex, score};

			for(MinimizerType minimizer : queryRead._minimizers){
				if(_minimizer_to_alignmentScoreQueue.find(minimizer) == _minimizer_to_alignmentScoreQueue.end()) continue;

				AlignmentScoreQueue& queue = _minimizer_to_alignmentScoreQueue[minimizer];

				if(queue.size() < 20){ //_parent._usedCoverageForCorrection
					queue.push(currentAlignmentScore);
				}

				const AlignmentScore& worseAlignmentScore = queue.top();

				if(currentAlignmentScore._score < worseAlignmentScore._score) continue;

				if(currentAlignmentScore._score == worseAlignmentScore._score){
					if(currentAlignmentScore._queryReadIndex > worseAlignmentScore._queryReadIndex) continue;
				}
				
				queue.pop();
				queue.push(currentAlignmentScore);
			}
			
			
			readWorkerLock.unlock(_referenceRead._readIndex);

		}

		vector<ReadType> extractBestAlignment(){

			unordered_set<ReadType> selectedQueryReadIndex;

			for(auto& it : _minimizer_to_alignmentScoreQueue){

				const MinimizerType& minimizer = it.first;
				AlignmentScoreQueue& queue = it.second;

				while(queue.size() > 0){
					const AlignmentScore& alignmentScore = queue.top();
					selectedQueryReadIndex.insert(alignmentScore._queryReadIndex);
					queue.pop();
				}

			}


			vector<ReadType> selectedQueryReadIndexVec;

			for(const auto& it : selectedQueryReadIndex){
				selectedQueryReadIndexVec.push_back(it);
			}

			return selectedQueryReadIndexVec;
		}
		/*
		void addMappedRead(const MinimizerRead& queryRead, const AlignmentResult2& alignment){
			#pragma omp critical(ReadWorkerAddMappedRead)
			{


				float score = alignment._nbMatches - alignment._nbMissmatches - alignment._nbInsertions - alignment._nbDeletions;
				
				for(MinimizerType minimizer : queryRead._minimizers){
					if(_minimizer_to_alignmentResultsNew.find(minimizer) == _minimizer_to_alignmentResultsNew.end()) continue;

					_minimizer_to_alignmentResultsNew[minimizer].push_back({alignment._queryReadIndex, score, alignment._identity});
				}
					

			}
		}
		*/
		
	};

	ReadWorkerLock _readWorkerLock;
	unordered_map<ReadType, ReadWorker> _readWorkers;


	string _readFilename;
	size_t _minimizerSize;
	size_t _kminmerSize;
	int _nbCores;
	float _minimizerDensity_assembly;
	float _minimizerDensity_correction;

	KminmerPosMap _kminmer_to_readIndex;
	u_int64_t _alignmentCheckSum;
	string _alignmentFilename;
	ofstream _alignmentFile;
	u_int64_t _maxChainingBand_lowDensity;
	size_t _usedCoverageForCorrection;

	ReadMapper(const string& readFilename, const string& alignmentFilename, size_t minimizerSize, float minimizerDensity_assembly, float minimizerDensity_correction, size_t usedCoverageForCorrection, int nbCores){

		_readFilename = readFilename;
		_alignmentFilename = alignmentFilename;
		_minimizerSize = minimizerSize;
		_kminmerSize = 2;
		_minimizerDensity_assembly = minimizerDensity_assembly;
		_minimizerDensity_correction = minimizerDensity_correction;
		_usedCoverageForCorrection = usedCoverageForCorrection;
		_nbCores = nbCores;

		_nextAlignmentReadIndexWriter = 0;
		_alignmentCheckSum = 0;
		_maxChainingBand_lowDensity = (u_int64_t) 2500 * _minimizerDensity_correction;

		_readWorkerLock.init(1000);
	}

	void execute(){



		//cout << "Indexing genomic kminmers" << endl;
		//indexGenomicKminmers();

		_alignmentFile = ofstream(_alignmentFilename);

		MinimizerReadParserParallel parser(_readFilename, _kminmerSize, false, true, 1000000, _nbCores);
		//parser._densityThreshold = _minimizerDensity_assembly;
		parser.execute(ProcessReadChunkFunctor(*this));

		_alignmentFile.close();

		cout << "Alignment check sum: " << _alignmentCheckSum << endl;
	}

	/*
	BloomCacheCoherent<u_int64_t>* _isKminmerVisited;
	BloomCacheCoherent<u_int64_t>* _isKminmerGenomic;

	void indexGenomicKminmers(){
                                                               
		_isKminmerVisited = new BloomCacheCoherent<u_int64_t>(64000000000ull);
		_isKminmerGenomic = new BloomCacheCoherent<u_int64_t>(64000000000ull);

		KminmerParserParallel parser(_readFilename, _minimizerSize, _kminmerSize, false, true, 1);
		parser._densityThreshold = _minimizerDensity_assembly;
		parser.parseSequences(IndexGenomicMinimizersFunctor(*this));

		delete _isKminmerVisited;

	}


	class IndexGenomicMinimizersFunctor {

		public:

		ReadMapper& _parent;

		IndexGenomicMinimizersFunctor(ReadMapper& parent) : _parent(parent){
		}

		IndexGenomicMinimizersFunctor(const IndexGenomicMinimizersFunctor& copy) : _parent(copy._parent){
		}

		~IndexGenomicMinimizersFunctor(){
		}


		void operator () (const KminmerList& kminmerList) {

			u_int32_t readIndex = kminmerList._readIndex;
			if(readIndex % 100000 == 0) cout << "\t\tIndexing genomic kminmers: " << readIndex << endl;

			MinimizerRead read = {kminmerList._readIndex, kminmerList._readMinimizers, kminmerList._minimizerPos, kminmerList._readQualities, kminmerList._readMinimizerDirections, kminmerList._readLength};
			//MinimizerRead readLowDensity = Utils::getLowDensityMinimizerRead(read, _parent._minimizerDensity_assembly);
			
			if(Utils::isReadTooShort(read)) return;

			vector<ReadKminmerComplete> kminmersInfos;
			MDBG::getKminmers_complete(_parent._kminmerSize, kminmerList._readMinimizers, kminmerList._minimizerPos, kminmersInfos, kminmerList._readIndex, kminmerList._readQualities);
		

			for(size_t i=0; i<kminmersInfos.size(); i++){

				const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];
				const KmerVec& vec = kminmerInfo._vec;
				u_int64_t vec_hash = vec.h();

				if(!_parent._isKminmerVisited->contains(vec_hash)){
					_parent._isKminmerVisited->insert(vec_hash);
					continue;
				}

				_parent._isKminmerGenomic->insert(vec_hash);

			}	

			
		}

	};
	*/

	class ProcessReadChunkFunctor {

		public:

		ReadMapper& _parent;

		ProcessReadChunkFunctor(ReadMapper& parent) : _parent(parent){
		}

		ProcessReadChunkFunctor(const ProcessReadChunkFunctor& copy) : _parent(copy._parent){
		}

		~ProcessReadChunkFunctor(){
		}


		void operator () (vector<MinimizerRead>& reads) const {
			cout << "Loaded read chunk: " << reads.size() << endl;

			_parent.processChunk(reads);
			_parent.clearChunk();
			//_parent._totalReadProcessed += _parent._reads.size();
			//_parent._reads.clear();

			//_parent.mapReadChunk();

			//_parent._kminmer_to_readIndex.clear();
			//_parent._refReads.clear();
		}
	};

	void processChunk(vector<MinimizerRead>& reads){

		//_queryReadIndex_to_referenceReadIndexes.clear();
		_readWorkers.clear();
		for(const MinimizerRead& read : reads){
			_readWorkers[read._readIndex] = ReadWorker(read);
		}

		cout << "Read worker size: " << _readWorkers.size() << endl;

		reads.clear();
		_kminmer_to_readIndex.clear();

		cout << "\tIndexing reads" << endl;
		processWorkerParallell(IndexReadsFunctor(*this), _readWorkers);

		//Remove unique seeds (wont match two reads)
		for (auto it = _kminmer_to_readIndex.begin(); it != _kminmer_to_readIndex.end();) {
			if((*it).second.size() <= 1) {
				it = _kminmer_to_readIndex.erase(it);
			}
			else{
				it++;
			}
		}

		cout << "\tMapping reads" << endl;
		mapReads();

		cout << "\tWriting alignments" << endl;
		writeAlignments();

		cout << "\tDone" << endl;

	}

	void clearChunk(){

		/*
		ReadWorker& readWorker = _readWorkers[86];
		while(!readWorker._bestReadIndexes.empty()){

			const ReadWorker::ReadScore& readScore = readWorker._bestReadIndexes.top();
			cout << readScore._readIndex << " " << readScore._score << endl;

			readWorker._bestReadIndexes.pop();
		}

		getchar();
		*/

		//delete _bloomFilter;
		_kminmer_to_readIndex.clear();
		_readWorkers.clear();
		//_queryReadIndex_to_referenceReadIndexes.clear();
	}


	ReadType _readIndex;
	vector<MinimizerType> _minimizers;
	vector<u_int32_t> _minimizersPos;
	vector<u_int8_t> _qualities;
	vector<u_int8_t> _readMinimizerDirections;
	float _debugNbMatches;



	class IndexReadsFunctor {

		public:

		ReadMapper& _parent;

		IndexReadsFunctor(ReadMapper& parent) : _parent(parent){
		}

		IndexReadsFunctor(const IndexReadsFunctor& copy) : _parent(copy._parent){
		}

		~IndexReadsFunctor(){
		}


		void operator () (const ReadType referenceReadIndex, ReadWorker& readWorker) {

			u_int32_t readIndex = referenceReadIndex;
			if(readIndex % 100000 == 0) cout << "\t\tIndexing read: " << readIndex << endl;

			vector<ReadKminmerComplete> kminmersInfos;
			MDBG::getKminmers_complete(_parent._kminmerSize, readWorker._referenceRead._minimizers, readWorker._referenceRead._minimizersPos, kminmersInfos, readWorker._referenceRead._readIndex, readWorker._referenceRead._qualities);
		

			for(size_t i=0; i<kminmersInfos.size(); i++){

				u_int16_t positionIndex = i;

				const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];
				const KmerVec& vec = kminmerInfo._vec;
				
				//if(!_parent._isKminmerGenomic->contains(vec.h())) continue;

				_parent._kminmer_to_readIndex.lazy_emplace_l(vec,
				[&readIndex, &positionIndex](KminmerPosMap::value_type& v) { // key exist
					v.second.push_back({readIndex, positionIndex});
					//v.second[0] += 1; //Increment kminmer count
					//if(std::find(v.second.begin(), v.second.end(), readIndex) == v.second.end()){
					//	v.second.push_back(readIndex);
					//}
				},           
				[&vec, &readIndex, &positionIndex](const KminmerPosMap::constructor& ctor) { // key inserted
					
					//vector<u_int32_t> readIndexes = {2}; //inital count of this kminmer
					vector<MinimizerPairPosition> readIndexes = {{readIndex, positionIndex}}; //inital count of this kminmer

					ctor(vec, readIndexes); 

				}); // construct value_type in place when key not present
			}
		
			
		}
	};



	template<typename Functor>
	void processReadParallell(const Functor& functor, const vector<MinimizerRead>& reads){

		u_int64_t i = 0;

		#pragma omp parallel num_threads(_nbCores)
		{

			Functor functorSub(functor);

			MinimizerRead read;
			bool isDone = false;


			while(true){

				
				#pragma omp critical(indexRefReads)
				{

					if(i >= reads.size()){
						isDone = true;
					}
					else{
						read = reads[i];
						i += 1;
					}
				}

				if(isDone) break;

				functorSub(read);

			}

			
		}


	}

	template<typename Functor>
	void processWorkerParallell(const Functor& functor, unordered_map<ReadType, ReadWorker>& readWorkers){

		vector<ReadType> readIndexes;
		for(const auto& it : readWorkers){
			readIndexes.push_back(it.first);
		}
		std::sort(readIndexes.begin(), readIndexes.end());

		u_int64_t i = 0;

		#pragma omp parallel num_threads(_nbCores)
		{

			Functor functorSub(functor);

			ReadType readIndex;
			//ReadWorker& read;
			bool isDone = false;


			while(true){

				
				#pragma omp critical(processWorkerParallell)
				{

					if(i >= readWorkers.size()){
						isDone = true;
					}
					else{
						readIndex = readIndexes[i];
						i += 1;
					}
				}

				if(isDone) break;

				functorSub(readIndex, readWorkers[readIndex]);

			}

			
		}
	}

	struct KminmerCount{
		KmerVec _vec;
		u_int32_t _count;
	};

	bool _print_debug;

	void mapReads(){

		_print_debug = false;

		int nbCores = _nbCores;

		//if(_print_debug){
			//nbCores = 1; cout << "align 1 cores" << endl;
		//}

		KminmerParserParallel parser(_readFilename, _minimizerSize, _kminmerSize, false, true, nbCores);
		//parser._densityThreshold = _minimizerDensity_assembly;
		parser.parseSequences(MapReadsFunctor(*this));
	}


	class MapReadsFunctor {

		public:

		ReadMapper& _parent;
		MinimizerAligner* _minimizerAligner;
		MinimizerChainer* _minimizerChainer;
		
		
		MapReadsFunctor(ReadMapper& parent) : _parent(parent){
			_minimizerAligner = nullptr;
			_minimizerChainer = nullptr;
		}

		MapReadsFunctor(const MapReadsFunctor& copy) : _parent(copy._parent){
			_minimizerAligner = new MinimizerAligner(3, -1, -1, -1);
			_minimizerChainer = new MinimizerChainer(_parent._minimizerSize);
		}

		~MapReadsFunctor(){
			if(_minimizerAligner != nullptr) delete _minimizerAligner;
			if(_minimizerChainer != nullptr) delete _minimizerChainer;
		}


		void operator () (const KminmerList& kminmerList) {

			if(kminmerList._readIndex % 10000 == 0) cout << "Aligning reads (low density): " << kminmerList._readIndex << endl;

			MinimizerRead queryRead = {kminmerList._readIndex, kminmerList._readMinimizers, kminmerList._minimizerPos, kminmerList._readQualities, kminmerList._readMinimizerDirections, kminmerList._readLength};
			//MinimizerRead queryReadLowDensity = Utils::getLowDensityMinimizerRead(queryRead, _parent._minimizerDensity_assembly);
			
			//if(readLowDensity._minimizers.size() < _parent._minReadLength){
			if(Utils::isReadTooShort(queryRead)){
				//_parent.writeRead(referenceRead._readIndex, referenceRead._minimizers, referenceRead._qualities);
				//_parent.writeAlignments(referenceRead._readIndex, {});
				return;
			}

			ReadMinimizerPositionMap queryReadMinimizerPositionMap;
			//unordered_set<MinimizerType> minimizerSet;
			//unordered_map<MinimizerType, u_int8_t> readMinimizer_to_direction;

			

			for(u_int32_t i=0; i<queryRead._minimizers.size(); i++){

				MinimizerType minimizer = queryRead._minimizers[i];
				u_int32_t position = queryRead._minimizersPos[i];
				bool isReversed = queryRead._readMinimizerDirections[i];

				MinimizerPosition minmizerPosition = {position, i, isReversed};
				queryReadMinimizerPositionMap[minimizer].push_back(minmizerPosition);
			}


			vector<ReadKminmerComplete> kminmersInfos;
			MDBG::getKminmers_complete(_parent._kminmerSize, queryRead._minimizers, queryRead._minimizersPos, kminmersInfos, queryRead._readIndex, queryRead._qualities);

			
			unordered_set<KmerVec> uniqueReferenceKminmers;
			for(size_t i=0; i<kminmersInfos.size(); i++){

				const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];
				const KmerVec& vec = kminmerInfo._vec;
				uniqueReferenceKminmers.insert(vec);
			}
			

			//vector<AlignmentResult2> alignments;

			unordered_map<ReadType, ReadMatchBound> referenceReadIndex_to_anchors2;
			//unordered_map<ReadType, std::pair<u_int16_t, u_int16_t>> queryReadIndex_to_anchorBounds;

			vector<KminmerCount> kminmerCounts;
			for(const KmerVec& vec : uniqueReferenceKminmers){

				if(_parent._kminmer_to_readIndex.find(vec) == _parent._kminmer_to_readIndex.end()) continue;

				u_int32_t count = _parent._kminmer_to_readIndex[vec].size();
				kminmerCounts.push_back({vec, count});

			}

			std::sort(kminmerCounts.begin(), kminmerCounts.end(), [](const KminmerCount& a, const KminmerCount& b){
				return a._count > b._count;
			});

			/*
			bool writeLol = false;
			if(kminmerCounts.size() > 0 && kminmerCounts[0]._count > 5000){
				writeLol = true;
				for(size_t i=0; i<kminmerCounts.size(); i++){
					cout << i << "\t" << kminmerCounts[i]._vec._kmers[0] << "\t" << kminmerCounts[i]._vec._kmers[1] << "\t" << kminmerCounts[i]._count << endl;
				}

				getchar();
			}
			*/


			size_t skippedIndex = 0;//3;
			if(kminmerCounts.size() <= 3){
				skippedIndex = 0;
			}

			//for(const KmerVec& vec : uniqueReferenceKminmers){
			for(size_t i=skippedIndex; i<kminmerCounts.size(); i++){

				//if(writeLol){
				//	cout << i << " " << kminmerCounts[i]._count << endl;
				//}
				
				const KmerVec& vec = kminmerCounts[i]._vec;

				if(_parent._kminmer_to_readIndex.find(vec) == _parent._kminmer_to_readIndex.end()) continue;

				for(const MinimizerPairPosition& referenceReadPosition : _parent._kminmer_to_readIndex[vec]){
					if(queryRead._readIndex == referenceReadPosition._readIndex) continue;


					if(referenceReadIndex_to_anchors2.find(referenceReadPosition._readIndex) == referenceReadIndex_to_anchors2.end()){
						referenceReadIndex_to_anchors2[referenceReadPosition._readIndex] = ReadMatchBound();
					}

					ReadMatchBound& bound = referenceReadIndex_to_anchors2[referenceReadPosition._readIndex];
					
					if(referenceReadPosition._positionIndex < bound._minIndex){
						bound._minIndex = referenceReadPosition._positionIndex;
					}
					if(referenceReadPosition._positionIndex > bound._maxIndex){
						bound._maxIndex = referenceReadPosition._positionIndex;
					}

					bound._nbMatches += 1;

				}

			}

			
			for(const auto& it : referenceReadIndex_to_anchors2){

				ReadType referenceReadIndex = it.first;

				const ReadMatchBound& readMatchBound = it.second;

				if(readMatchBound._nbMatches < 1) continue;

				const MinimizerRead& referenceRead = _parent._readWorkers[referenceReadIndex]._referenceRead;// = _parent._mReads[queryReadIndex];

				//cout << referenceRead._minimizers.size() << " " << queryRead._minimizers.size() << endl;
				int64_t startQueryPositionIndex = getStartQueryPositionIndex(referenceRead, readMatchBound._minIndex, 5000, queryReadMinimizerPositionMap);
				int64_t endQueryPositionIndex = getEndQueryPositionIndex(referenceRead, readMatchBound._maxIndex, 5000, queryReadMinimizerPositionMap);
				
				//startQueryPositionIndex = 0;
				//endQueryPositionIndex = queryRead._minimizers.size()-1;

				
				vector<Anchor> anchors;

				//cout << startQueryPositionIndex << " " << endQueryPositionIndex << endl;
				for(size_t i=startQueryPositionIndex; i<=endQueryPositionIndex; i++){

					MinimizerType referenceMinimizer = referenceRead._minimizers[i];

					if(queryReadMinimizerPositionMap.find(referenceMinimizer) == queryReadMinimizerPositionMap.end()) continue;

					u_int32_t referencePosition = referenceRead._minimizersPos[i];
					bool referenceIsReversed = referenceRead._readMinimizerDirections[i];

					const vector<MinimizerPosition>& queryMinimizerPositions = queryReadMinimizerPositionMap[referenceMinimizer];
					
					for(const MinimizerPosition& queryMinimizerPosition : queryMinimizerPositions){
						anchors.push_back({referencePosition, queryMinimizerPosition._position, queryMinimizerPosition._isReversed != referenceIsReversed, i, queryMinimizerPosition._positionIndex});
					}



				}



				//cout << "Anchors: " << anchors.size() << endl;



				AlignmentResult2 chainingAlignment = _minimizerChainer->computeChainingAlignment(anchors, referenceRead, queryRead, _minimizerAligner, _parent._maxChainingBand_lowDensity);
				
				//cout << referenceRead._readIndex << " " << queryRead._readIndex << " " << chainingAlignment._overHangStart << " " << chainingAlignment._overHangEnd << " " << chainingAlignment._nbMatches << " " << chainingAlignment._nbMissmatches << endl;
				//getchar();
				//AlignmentResult2 chainingAlignment = _parent.computeAlignment(anchors, referenceRead, queryRead, _minimizerAligner);
				
				//cout << mashDistance << " " << chainingAlignment._divergence << endl;
				//if(chainingAlignment._nbMatches - chainingAlignment._nbMissmatches - chainingAlignment._nbInsertions - chainingAlignment._nbDeletions < 0) continue;
				//if(chainingAlignment._nbMatches - chainingAlignment._nbMissmatches < 5) continue;
				//if(chainingAlignment._nbMatches < 1000*_minimizerDensity_correction) continue;
				//if(chainingAlignment._overHangStart > 5000) continue;
				//if(chainingAlignment._overHangEnd > 5000) continue;
				//if(chainingAlignment._nbMatches - chainingAlignment._nbMissmatches < 0) continue;
				//if(chainingAlignment._divergence > 0.04) continue;

				if(chainingAlignment._alignments.empty()) continue;

				//alignments.push_back(chainingAlignment);
				
				_parent._readWorkers[referenceReadIndex].addAlignmentScore(queryRead, chainingAlignment, _parent._readWorkerLock);
			}
			
			/*
			if(_eval_correction){
				
				if(_simulatedReadMetadata[read._readIndex]._identity < 98) return read;
				if(_simulatedReadMetadata[read._readIndex]._isReversed) return read;
				
				cout << endl << "Read: " <<  read._readIndex << " " << read._minimizers.size() << " " << _mReads[read._readIndex]._meanReadQuality << endl;
				std::sort(alignments.begin(), alignments.end(), [](const AlignmentResult2& a, const AlignmentResult2& b){
					return a.getScore() > b.getScore();
				});

				cout << "Truth info:" << endl;
				const SimulatedReadMetadata& refMetadata = _simulatedReadMetadata[read._readIndex];
				cout << refMetadata.toString() << endl;

				cout << "Alignment info:" << endl;
				for(size_t i=0; i<alignments.size() && i < 50; i++){
					const AlignmentResult2& alignment = alignments[i];
					const SimulatedReadMetadata& metadata = _simulatedReadMetadata[alignment._queryReadIndex];
					cout << "\t" << i << ":\t" << metadata.toString() << "\t" << alignment._nbMatches << "\t" << alignment._nbMissmatches << "\t" << alignment._nbInsertions << "\t" << alignment._nbDeletions << "\t" << metadata.distanceFrom(refMetadata) << "\t" << alignment._divergence << "\t" << alignment.getSimilarity()<< endl;
				}

				getchar();
				
				//if(refMetadata._identity > 98.8) getchar();
				
				

			}
			*/
			/*
			//cout << referenceRead._readIndex << " " << alignments.size() << endl;
			//if(alignments.size() == 0){
			//	_parent.writeAlignments(referenceRead._readIndex, {});
			//	return;
			//}
			//cout << alignments.size() << endl;
			//cout << "copy here" << endl;
			vector<AlignmentResult2> bestAlignmentsLowDensity = _parent.selectBestAlignments(referenceRead, alignments, referenceReadMinimizerPositionMap, _parent._usedCoverageLowDensity);
			
			//const MinimizerRead& correctedRead = performPoaCorrection4(referenceRead, bestAlignmentsLowDensity, referenceReadMinimizerPositionMap);
			
			//_parent.writeRead(correctedRead._readIndex, correctedRead._minimizers, correctedRead._qualities);

			vector<ReadType> alignedQueryReads;
			for(const AlignmentResult2& alignment : bestAlignmentsLowDensity){
				alignedQueryReads.push_back(alignment._queryReadIndex);
			}

			//cout << referenceRead._readIndex << " " << alignedQueryReads.size() << endl;

			_parent.writeAlignments(referenceRead._readIndex, alignedQueryReads);
		
			*/
			//cout << "todo: select best reads, writealignment" << endl;
		

		}

		int64_t getStartQueryPositionIndex(const MinimizerRead& queryRead, int64_t startIndex, int64_t maxGapLength, ReadMinimizerPositionMap& referenceReadMinimizerPositionMap){
			
			//maxGapLength *= 2;
			int64_t lastMacthingIndex = startIndex;
			int64_t i = startIndex;
			int64_t currentGapLength = 0;
			int64_t prevPosition = queryRead._minimizersPos[i];
			i -= 1;

			while(true){
				if(i < 0) break;
				if(currentGapLength > maxGapLength) break;

				const MinimizerType currentMinimizer = queryRead._minimizers[i];

				if(referenceReadMinimizerPositionMap.find(currentMinimizer) == referenceReadMinimizerPositionMap.end()){
					currentGapLength += (prevPosition-queryRead._minimizersPos[i]);
				}
				else{
					currentGapLength = 0;
					lastMacthingIndex = i;
				}


				prevPosition = queryRead._minimizersPos[i];
				i -= 1;
			}

			return lastMacthingIndex;
		}

		int64_t getEndQueryPositionIndex(const MinimizerRead& queryRead, int64_t startIndex, int64_t maxGapLength, ReadMinimizerPositionMap& referenceReadMinimizerPositionMap){
			
			int64_t lastMacthingIndex = startIndex;
			int64_t i = startIndex;
			int64_t currentGapLength = 0;
			int64_t prevPosition = queryRead._minimizersPos[i];
			i += 1;

			while(true){
				if(i >= queryRead._minimizers.size()) break;
				if(currentGapLength > maxGapLength) break;

				const MinimizerType currentMinimizer = queryRead._minimizers[i];

				if(referenceReadMinimizerPositionMap.find(currentMinimizer) == referenceReadMinimizerPositionMap.end()){
					currentGapLength += (queryRead._minimizersPos[i]-prevPosition);
				}
				else{
					currentGapLength = 0;
					lastMacthingIndex = i;
				}


				prevPosition = queryRead._minimizersPos[i];
				i += 1;
			}

			return lastMacthingIndex;
		}
	
		
	};





	void writeAlignments(){

		vector<ReadType> readIndexes;
		for(auto& it : _readWorkers){
			readIndexes.push_back(it.first);
		}

		std::sort(readIndexes.begin(), readIndexes.end());

		for(const ReadType& referenceReadIndex : readIndexes){

			//const ReadType& referenceReadIndex = it.first;
			ReadWorker& readWorker = _readWorkers[referenceReadIndex];
			
			const vector<ReadType> alignedQueryReads = readWorker.extractBestAlignment();
			const u_int32_t nbAlignedReads = alignedQueryReads.size();
			//cout << nbAlignedReads << endl;

			if(nbAlignedReads > 0){

				for(ReadType queryReadIndex : alignedQueryReads){
					_alignmentCheckSum += referenceReadIndex * nbAlignedReads * queryReadIndex;
				}

				_alignmentFile.write((const char*)&referenceReadIndex, sizeof(referenceReadIndex));
				_alignmentFile.write((const char*)&nbAlignedReads, sizeof(nbAlignedReads));
				_alignmentFile.write((const char*)&alignedQueryReads[0], nbAlignedReads*sizeof(ReadType));
				

				if(referenceReadIndex % 10000 == 0){
					Logger::get().debug() << "\tAlignment checksum: " << referenceReadIndex << " " << _alignmentCheckSum;
					//if(readWriter._referenceReadIndex > 1000) getchar();
				}
				
			}

		}




	}


	/*
	vector<AlignmentResult2> selectBestAlignments(const MinimizerRead& referenceRead, const vector<AlignmentResult2>& alignments, ReadMapper::ReadMinimizerPositionMap& referenceReadMinimizerPositionMap, const size_t maxCoverage){
		

		vector<AlignmentResult2> selectedAlignments;
		
		


		unordered_set<ReadType> selectedReads;

		for(const auto& it : minimizer_to_alignmentResultsNew){
			u_int64_t minimizer = it.first;
			const vector<AlignmentScore>& alignmentResults = it.second;

			const vector<AlignmentScore>& bestAlignmentResults  = getBestAlignments(alignmentResults, maxCoverage);

			for(const AlignmentScore& alignmentResult : bestAlignmentResults){
				//if(alignmentResult.divergence() > 0.02) continue;
				//cout << alignmentResult._readIndex << " " << alignmentResult.score() << " " << readIndex_to_matchScore[alignmentResult._readIndex] << endl;
				selectedReads.insert(alignmentResult._queryReadIndex);
			}

			//cout << "check" << endl;
		}



		for(const AlignmentResult2& alignment : alignments){

			if(selectedReads.find(alignment._queryReadIndex) == selectedReads.end()) continue;

			selectedAlignments.push_back(alignment);

		}

		//cout << "lala: " << selectedReads.size() << endl;
		return selectedAlignments;

	}
	

	

	vector<AlignmentScore> getBestAlignments(vector<AlignmentScore> alignmentResults, int n){

		vector<AlignmentScore> bestAlignments;

		std::sort(alignmentResults.begin(), alignmentResults.end(), [](const AlignmentScore& a, const AlignmentScore& b){
			if(a._score == b._score){
				if(a._identity == b._identity){
					return a._queryReadIndex > b._queryReadIndex;
				}
				return a._identity > b._identity;
			}
			return a._score > b._score;
		});

		for(int i=0; i<n && i<alignmentResults.size(); i++){

			bool isHere = false;
			for(const AlignmentScore& al : bestAlignments){
				if(al._queryReadIndex == alignmentResults[i]._queryReadIndex){ //The same alignment can be present several time if a minimizer is repeated in a read
					isHere = true;
					break;
				}
			}

			if(isHere) continue;

			bestAlignments.push_back(alignmentResults[i]);
		}

		return bestAlignments;
	}
	*/

    struct AlignmentWriter{
        ReadType _referenceReadIndex;
        vector<ReadType> _alignedQueryReads;
    };

    struct AlignmentWriter_Comparator {
        bool operator()(const AlignmentWriter& p1, const AlignmentWriter& p2){
            return p1._referenceReadIndex > p2._referenceReadIndex;
        }
    };

	priority_queue<AlignmentWriter, vector<AlignmentWriter> , AlignmentWriter_Comparator> _alignmentWriterQueue;
	ReadType _nextAlignmentReadIndexWriter;

	void writeAlignments(ReadType referenceReadIndex, const vector<ReadType>& alignedQueryReads){

		#pragma omp critical
		{
			
			_alignmentWriterQueue.push({referenceReadIndex, alignedQueryReads});

			while(!_alignmentWriterQueue.empty()){


				const AlignmentWriter& readWriter = _alignmentWriterQueue.top();

				if(readWriter._referenceReadIndex == _nextAlignmentReadIndexWriter){


					const ReadType referenceReadIndexCurrent = readWriter._referenceReadIndex;
					const vector<ReadType> alignedQueryReadsCurrent = readWriter._alignedQueryReads;
					const u_int16_t nbAlignedReads = alignedQueryReadsCurrent.size();

					if(nbAlignedReads > 0){

						for(ReadType queryReadIndex : alignedQueryReadsCurrent){
							_alignmentCheckSum += referenceReadIndexCurrent * nbAlignedReads * queryReadIndex;
						}

						_alignmentFile.write((const char*)&referenceReadIndexCurrent, sizeof(referenceReadIndexCurrent));
						_alignmentFile.write((const char*)&nbAlignedReads, sizeof(nbAlignedReads));
						_alignmentFile.write((const char*)&alignedQueryReadsCurrent[0], nbAlignedReads*sizeof(ReadType));
						

						if(readWriter._referenceReadIndex % 10000 == 0){
							cout << readWriter._referenceReadIndex << " " << _alignmentCheckSum << endl;
							//if(readWriter._referenceReadIndex > 1000) getchar();
						}
						
					}


					
					_alignmentWriterQueue.pop();
					_nextAlignmentReadIndexWriter += 1;
				}
				else{
					break;
				}
			}
			
		}

	}



};

#endif