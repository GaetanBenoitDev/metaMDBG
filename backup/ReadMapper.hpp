
#ifndef MDBG_METAG_READMAPPER
#define MDBG_METAG_READMAPPER

#include "../Commons.hpp"



class ReadMapper{
    
public:

	typedef phmap::parallel_flat_hash_map<KmerVec, vector<u_int32_t>, phmap::priv::hash_default_hash<KmerVec>, phmap::priv::hash_default_eq<KmerVec>, std::allocator<std::pair<KmerVec, vector<u_int32_t>>>, 4, std::mutex> KminmerReadMap;

	string _refFilename;
	string _queryFilename;
	size_t _minReadSize;
	size_t _nbCores;

	size_t _kminmerSize;

	//vector<MinimizerRead> _refReads;
	vector<MinimizerRead> _queryReads;
	KminmerReadMap _kminmer_to_readIndex;

	class ReadWorker{

		const static u_int64_t _maxReads = 100;

		public:

		struct ReadScore{
			u_int64_t _readIndex;
			float _score;
		};

		vector<ReadScore> _readScores;
		//struct ReadComparator{
		//	bool operator() (const ReadScore& a, const ReadScore& b) { return a._score > b._score; }
		//};

		//std::priority_queue<ReadScore, vector<ReadScore>, ReadComparator> _bestReadIndexes;

		void addRead(u_int64_t readIndex, float score){

			#pragma omp critical(ReadWorker_addRead)
			{

				_readScores.push_back({readIndex, score});
				/*
				if(_bestReadIndexes.size() < _maxReads){
					_bestReadIndexes.push({readIndex, score});
				}
				else{
					if(score > _bestReadIndexes.top()._score){
						_bestReadIndexes.pop();
						_bestReadIndexes.push({readIndex, score});
					}
				}
				*/
			}

		}

	};

	unordered_map<u_int64_t, ReadWorker> _readWorkers;

	ReadMapper(const string& refFilename, const string& queryFilename, size_t minReadSize, int nbCores){

		_refFilename = refFilename;
		_queryFilename = queryFilename;
		_minReadSize = minReadSize;
		_nbCores = nbCores;

		_kminmerSize = 2;
	}

	void execute(){


		cout << "Loading minimizer reads" << endl;
		loadReads(_queryFilename);

		cout << "Mapping reads" << endl;
		mapReads();

	}


	void loadReads(const string& queryFilename){


		KminmerParserParallel parser(queryFilename, 13, _kminmerSize, false, true, 1);
		parser.parseSequences(LoadMinimizerReadsFunctor(*this));
	}


	class LoadMinimizerReadsFunctor {

		public:

		ReadMapper& _parent;

		LoadMinimizerReadsFunctor(ReadMapper& parent) : _parent(parent){
		}

		LoadMinimizerReadsFunctor(const LoadMinimizerReadsFunctor& copy) : _parent(copy._parent){
		}

		~LoadMinimizerReadsFunctor(){
		}


		void operator () (const KminmerList& kminmerList) {

			u_int32_t readIndex = kminmerList._readIndex;
			if(readIndex % 100000 == 0) cout << "\tLoading query reads: " << readIndex << endl;


			MinimizerRead read = {kminmerList._readIndex, kminmerList._readMinimizers, kminmerList._minimizerPos, kminmerList._readQualities, kminmerList._readMinimizerDirections};
			
			if(kminmerList._readMinimizers.size() < _parent._minReadSize){
				_parent._queryReads.push_back({});
			}
			else{
				_parent._queryReads.push_back(read);
			}
			

		}
	};


	void mapReads(){


		MinimizerReadParserParallel parser(_refFilename, _kminmerSize, false, true, 1000000, _nbCores);
		parser.execute(ProcessReadChunkFunctor(*this));

	}


	class ProcessReadChunkFunctor {

		public:

		ReadMapper& _parent;

		ProcessReadChunkFunctor(ReadMapper& parent) : _parent(parent){
		}

		ProcessReadChunkFunctor(const ProcessReadChunkFunctor& copy) : _parent(copy._parent){
		}

		~ProcessReadChunkFunctor(){
		}


		void operator () (const vector<MinimizerRead>& reads) const {
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

	void processChunk(const vector<MinimizerRead>& reads){

		cout << "Indexing reads" << endl;
		processReadParallell(IndexReadsFunctor(*this), reads);

		for(const MinimizerRead& read : reads){
			_readWorkers[read._readIndex] = ReadWorker();
		}
		//_readWorkers.resize(reads.size());

		cout << "Mapping reads" << endl;
		processReadParallell(MapReadFunctor(*this), _queryReads);
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

		_kminmer_to_readIndex.clear();
		_readWorkers.clear();
	}

	class IndexReadsFunctor {

		public:

		ReadMapper& _parent;

		IndexReadsFunctor(ReadMapper& parent) : _parent(parent){
		}

		IndexReadsFunctor(const IndexReadsFunctor& copy) : _parent(copy._parent){
		}

		~IndexReadsFunctor(){
		}


		void operator () (const MinimizerRead& read) {


			if(read._minimizers.size() < _parent._minReadSize) return;

			vector<u_int64_t> rlePositions;
			vector<ReadKminmerComplete> kminmersInfo;
			//MDBG::getKminmers(_l, _k, minimizers, minimizersPos, kminmers, kminmersInfo, rlePositions, 0, false);
			MDBG::getKminmers_complete(_parent._kminmerSize, read._minimizers, read._minimizersPos, kminmersInfo, read._readIndex, read._qualities);
				
			u_int32_t readIndex = read._readIndex;
			//if(readIndex % 100000 == 0) cout << "\tIndexing read: " << readIndex << endl;

			for(u_int32_t i=0; i<kminmersInfo.size(); i++){
			
				
				const ReadKminmerComplete& kminmerInfo = kminmersInfo[i];

				const KmerVec& vec = kminmerInfo._vec;
			

			
				_parent._kminmer_to_readIndex.lazy_emplace_l(vec,
				[&readIndex](KminmerReadMap::value_type& v) { // key exist
					if(std::find(v.second.begin(), v.second.end(), readIndex) == v.second.end()){
						v.second.push_back(readIndex);
					}
				},           
				[&vec, &readIndex](const KminmerReadMap::constructor& ctor) { // key inserted
					
					vector<u_int32_t> readIndexes = {readIndex};

					ctor(vec, readIndexes); 

				}); // construct value_type in place when key not present



			}

			
		}
	};

	class MapReadFunctor {

		public:

		ReadMapper& _parent;

		MapReadFunctor(ReadMapper& parent) : _parent(parent){
		}

		MapReadFunctor(const MapReadFunctor& copy) : _parent(copy._parent){
		}

		~MapReadFunctor(){
		}


		void operator () (const MinimizerRead& queryRead) {

			if(queryRead._minimizers.size() == 0) return; //Small read

			//"reprise: bug here"
			//cout << "lala: " <<queryRead._readIndex << " " << queryRead._minimizers.size() << endl;
			//getchar();
			if(queryRead._readIndex % 100000 == 0) cout << "\tMapping reads: " << queryRead._readIndex << endl;

			if(queryRead._minimizers.size() < _parent._minReadSize) return;

			vector<u_int64_t> rlePositions;
			vector<ReadKminmerComplete> kminmersInfos;
			//MDBG::getKminmers(_l, _k, minimizers, minimizersPos, kminmers, kminmersInfo, rlePositions, 0, false);
			MDBG::getKminmers_complete(_parent._kminmerSize, queryRead._minimizers, queryRead._minimizersPos, kminmersInfos, queryRead._readIndex, queryRead._qualities);
				
			//u_int32_t readIndex = read._readIndex;

			unordered_map<u_int32_t, u_int32_t> readIndex_to_matchCount;
			vector<KmerVec> readKminmers;

			for(const ReadKminmerComplete& kminmerInfo : kminmersInfos){

				const KmerVec& vec = kminmerInfo._vec;
				readKminmers.push_back(vec);

				if(_parent._kminmer_to_readIndex.find(vec) == _parent._kminmer_to_readIndex.end()) continue;

				for(u_int32_t referenceReadIndex : _parent._kminmer_to_readIndex[vec]){
					if(queryRead._readIndex == referenceReadIndex) continue; //Currently corrected read
					readIndex_to_matchCount[referenceReadIndex] += 1;
				}

			}

		

			for(const auto& it : readIndex_to_matchCount){

				u_int64_t referenceReadIndex = it.first;
				int64_t nbMatches = it.second;

				if(nbMatches < 2) continue;

				if(_parent._readWorkers.find(referenceReadIndex) == _parent._readWorkers.end()){
					cout << "error" << endl;
					exit(1);
				}

				_parent._readWorkers[referenceReadIndex].addRead(queryRead._readIndex, nbMatches);

				//if(referenceReadIndex == 86){

				//	cout << "Add read: " << queryRead._readIndex << " " << nbMatches << endl;
				//}
				//int64_t alignScoreFast = getAlignmentScoreFast(_mReads[readIndex], minimizerSet, minimizerPosition, readMinimizer_to_direction);
			

				//matchingReadIndexes.push_back({readIndex, nbMatches});
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

};

#endif