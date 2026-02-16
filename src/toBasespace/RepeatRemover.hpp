



#ifndef MDBG_METAG_RepeatRemover
#define MDBG_METAG_RepeatRemover

#include "../Commons.hpp"


class RepeatRemover {
  

public:

	typedef phmap::parallel_flat_hash_map<u_int128_t, vector<ReadType>, phmap::priv::hash_default_hash<u_int128_t>, phmap::priv::hash_default_eq<u_int128_t>, std::allocator<std::pair<u_int128_t, vector<ReadType>>>, 4, std::mutex> KminmerPosMap;
	typedef phmap::parallel_flat_hash_map<u_int128_t, u_int32_t, phmap::priv::hash_default_hash<u_int128_t>, phmap::priv::hash_default_eq<u_int128_t>, std::allocator<std::pair<u_int128_t, u_int32_t>>, 4, std::mutex> KminmerMap;
	//typedef phmap::parallel_flat_hash_map<ReadType, u_int32_t, phmap::priv::hash_default_hash<ReadType>, phmap::priv::hash_default_eq<ReadType>, std::allocator<std::pair<ReadType, u_int32_t>>, 4, std::mutex> ReadLengthMap;

	string _inputDir;
	string _inputFilenameContig;
	string _outputFilenameContig;
	size_t _kminmerSize;
	float _minimizerDensity;
	int _nbCores;
	bool _hasQuality;
	
	KminmerMap _kminmer_to_unitigIndex;
	KminmerMap _kminmer_to_abundance;
	KminmerPosMap _kminmer_to_readIndexes;
	ofstream _outputFile;
	unordered_map<ReadType, u_int32_t> _readIndex_to_readLength;

	struct ReadAbundance{
		u_int32_t _readLength;
		float _readAbundance;
	};

	//unordered_map<u_int128_t, ReadAbundance> _kminmer_to_readAbundance;

	RepeatRemover(const string& inputDir, const string& inputFilenameContig, const string& outputFilenameContig, size_t kminmerSize, float minimizerDensity, bool hasQuality, int nbCores){
		_inputDir = inputDir;
		_inputFilenameContig = inputFilenameContig;
		_outputFilenameContig = outputFilenameContig;
		_kminmerSize = kminmerSize;
		_minimizerDensity = minimizerDensity;
		_hasQuality = hasQuality;
		_nbCores = nbCores;
	}
	

	void execute(){

		_outputFile.open(_outputFilenameContig);

		processContigs();
		
		_outputFile.close();


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
	}
	
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

		_kminmer_to_unitigIndex.clear();
		_kminmer_to_abundance.clear();
		_kminmer_to_readIndexes.clear();
		_readIndex_to_readLength.clear();
		
		Logger::get().debug() << "\tLoad contigs";
		processChunkParallell(IndexContigKminmerFunctor(*this));

		Logger::get().debug() << "\tIndexing initial unitigs";
		loadUnitigIndex();

		Logger::get().debug() << "\tIndexing kminmer abundance";
		loadKminmerAbundance();

		Logger::get().debug() << "\tIndexing reads";
		indexReads();

		Logger::get().debug() << "\tBreak unbridged repeats";
		processChunkParallell(BreakUnbridgedRepeatsFunctor(*this));

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

	void loadUnitigIndex(){

		//const string& unitigFilename = _inputDir + "/unitig_data.txt.init";
		const string& unitigFilename = _inputDir + "/unitig_data.txt.init.k" + to_string(_kminmerSize);

		KminmerParserParallel parser(unitigFilename, -1, _kminmerSize, false, false, _nbCores);
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

			u_int64_t readIndex = kminmerList._readIndex;

			const vector<MinimizerType>& readMinimizers = kminmerList._readMinimizers;
			const vector<ReadKminmerComplete>& kminmersInfos = kminmerList._kminmersInfo;



			for(size_t i=0; i<kminmersInfos.size(); i++){

				const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];
				const KmerVec& vec = kminmerInfo._vec;
				u_int128_t vecHash = vec.hash128();

				if(_parent._kminmer_to_readIndexes.find(vecHash) == _parent._kminmer_to_readIndexes.end()) continue;

				_parent._kminmer_to_unitigIndex.lazy_emplace_l(vecHash, 
				[this](KminmerMap::value_type& v) { // key exist
				},           
				[&vecHash, &readIndex](const KminmerMap::constructor& ctor) { // key inserted
					
					ctor(vecHash, readIndex); 

				}); // construct value_type in place when key not present
				
			}


			/*
			#pragma omp critical(IndexUnitigFunctor)
			{

				for(size_t i=0; i<kminmersInfos.size(); i++){

					const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];
					const KmerVec& vec = kminmerInfo._vec;

					_parent._kminmer_to_unitigIndex[vec] = readIndex;


				}
			}
			*/
		}
	};
	
	//void loadContigs(){

	//}
	
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

			/*
			if(abundances.size() > 20){

				float readAbundance = Utils::compute_median_float(abundances);

				#pragma omp critical(IndexReadsFunctor)
				{

					for(size_t i=0; i<kminmerList._kminmersInfo.size(); i++){
						
						const ReadKminmerComplete& kminmerInfo = kminmerList._kminmersInfo[i];
						const KmerVec& vec = kminmerInfo._vec;
						const u_int128_t& vecHash = vec.hash128();
						
						if(_parent._kminmer_to_readAbundance.find(vecHash) == _parent._kminmer_to_readAbundance.end()){
							_parent._kminmer_to_readAbundance[vecHash] = {(u_int32_t) kminmerList._kminmersInfo.size(), readAbundance};
						}
						else{
							
							ReadAbundance existingReadAbundance = _parent._kminmer_to_readAbundance[vecHash];

							float abundance = min(existingReadAbundance._readAbundance, readAbundance);
							_parent._kminmer_to_readAbundance[vecHash] = {(u_int32_t) kminmerList._kminmersInfo.size(), abundance};
							//if(kminmerList._kminmersInfo.size() > existingReadAbundance._readLength){
							//	_parent._kminmer_to_readAbundance[vecHash] = {(u_int32_t) kminmerList._kminmersInfo.size(), readAbundance};
							//}
						}
					}
				}
			}
			*/

			

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


	void loadKminmerAbundance(){

		ifstream kminmerAbundanceFile(_inputDir + "/kminmerData_abundance_init_k" + to_string(_kminmerSize) + ".txt");
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
			if(_kminmer_to_readIndexes.find(vecHash) == _kminmer_to_readIndexes.end()) continue;

			_kminmer_to_abundance[vecHash] = abundance;

			//cout << "lala" << " " << abundance << endl;

		}

		kminmerAbundanceFile.close();
	}
	
	struct Fragment{
		u_int32_t _fragmentIndex;
		u_int32_t _startPos;
		u_int32_t _endPos;
		u_int32_t _length;
		long double _coverage;
		int32_t _finalContigIndex;
		vector<ReadType> _readIndexes;
		unordered_map<u_int32_t, u_int32_t> _nbBridgingReads;
		u_int32_t _maxReadLength;
	};

	struct FragmentPath{
		u_int32_t _fragmentIndexStart;
		u_int32_t _fragmentIndexEnd;

		u_int32_t getPathLength() const{
			return _fragmentIndexEnd - _fragmentIndexStart;
		}
	};




	class BreakUnbridgedRepeatsFunctor {

		public:

		RepeatRemover& _parent;

		BreakUnbridgedRepeatsFunctor(RepeatRemover& parent) : _parent(parent){
		}

		BreakUnbridgedRepeatsFunctor(const BreakUnbridgedRepeatsFunctor& copy) : _parent(copy._parent){
		}

		~BreakUnbridgedRepeatsFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {

			ReadType readIndex = kminmerList._readIndex;
			const vector<MinimizerType>& readMinimizers = kminmerList._readMinimizers;
			const vector<ReadKminmerComplete>& kminmersInfos = kminmerList._kminmersInfo;
			u_int8_t isCircular = kminmerList._isCircular;

			if(isCircular){
				writeContig(readMinimizers, isCircular);
				return;
			}

			//if(readIndex != 30261 && readIndex != 30260 && readIndex != 30262) return;
			//if(readMinimizers.size() < 12000) return;

			//cout << endl << endl << readIndex << " " << readMinimizers.size() << " " << kminmersInfos.size() << endl;

			u_int32_t fragmentIndex = 0;
			u_int32_t lastUnitigIndex = -1;
			u_int32_t fragmentStartPos = 0;
			vector<float> fragmentAbundances;
			vector<Fragment> fragments;

			for(size_t i=0; i<kminmersInfos.size(); i++){
				
				const KmerVec& vec = kminmersInfos[i]._vec;
				const u_int128_t& vecHash = vec.hash128();

				if(_parent._kminmer_to_abundance.find(vecHash) != _parent._kminmer_to_abundance.end()){
					fragmentAbundances.push_back(_parent._kminmer_to_abundance[vecHash]);
				}
				else{
					fragmentAbundances.push_back(1);
				}
			}

			for(size_t i=0; i<kminmersInfos.size(); i++){

				const KmerVec& vec = kminmersInfos[i]._vec;
				const u_int128_t& vecHash = vec.hash128();

				//cout << i << "\t" << _parent._kminmer_to_abundance[vecHash] << "\t" << _parent._kminmer_to_readAbundance[vecHash]._readAbundance << endl;//<< " " << _parent._kminmer_to_readAbundance[vecHash]._readLength << endl;
				if(_parent._kminmer_to_unitigIndex.find(vecHash) == _parent._kminmer_to_unitigIndex.end()) continue;


				const u_int32_t& unitigIndex = _parent._kminmer_to_unitigIndex[vecHash];


				//if(readIndex == 36203){
				//	cout << i << " " << unitigIndex << endl;
				//}

				//cout << i << " " << unitigIndex << " " << _parent._kminmer_to_abundance[vecHash] << " " << i*200 << endl;
				//getchar();
				
				if(unitigIndex != lastUnitigIndex || i == kminmersInfos.size()-1){

					lastUnitigIndex = unitigIndex;

					if(i == 0){
						//fragmentUnitigIndex = unitigIndex;
						continue;
					}

					
					u_int32_t fragmentEndPos = i-1; //rlePositions[minimizerPos[i]];
					if(i == kminmersInfos.size()-1){
						fragmentEndPos = kminmersInfos.size()-1;
					}

					u_int32_t fragmentLength = fragmentEndPos - fragmentStartPos + 1;


					
					long double sum = 0;
					long double n = 0;

					for(size_t i=fragmentStartPos; i <= fragmentEndPos; i++){
						sum += fragmentAbundances[i];
						n += 1;
					}

					long double fragmentCoverage = 0;
					if(n > 0){
						fragmentCoverage = sum / n;
					}

					/*
					long double sumSmooth = 0;
					long double nSmooth = 0;

					for(size_t i=fragmentStartPos; i <= fragmentEndPos; i++){

						const KmerVec& vec = kminmersInfos[i]._vec;
						const u_int128_t& vecHash = vec.hash128();

						if(_parent._kminmer_to_readAbundance.find(vecHash) != _parent._kminmer_to_readAbundance.end()){
							const ReadAbundance& readAbundance = _parent._kminmer_to_readAbundance[vecHash];
							if(readAbundance._readLength > fragmentLength){
								sumSmooth += readAbundance._readAbundance;
								nSmooth += 1;
							}
						}

					}
						
					if(nSmooth > 0){
						float fragmentCoverageSmoothed = sumSmooth / nSmooth;
						fragmentCoverage = max(fragmentCoverage, fragmentCoverageSmoothed);
					}
					*/
					
					//cout << fragmentIndex << "\t" << fragmentStartPos << "\t" << fragmentEndPos << "\t" << fragmentLength << "\t" << fragmentCoverage << endl; //<< "\t" << getCoverageInBounds(repeats, fragmentStartPos, fragmentEndPos) << endl;

					//for(size_t i=fragmentStartPos; i <= fragmentEndPos; i++){

					//}

					Fragment fragment = {fragmentIndex, fragmentStartPos, fragmentEndPos, fragmentLength, fragmentCoverage, -1, {}};
					fragments.push_back(fragment);
					fragmentStartPos = i;//rlePositions[minimizerPos[i]];
					fragmentIndex += 1;

					//fragmentAbundances.clear();
					
				}
			}




			if(fragments.size() == 0){
				writeContig(readMinimizers, isCircular);
				return;
			}

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
		
			computeBridgingReads(readIndex, fragments);

			/*
			if(readIndex == 36203){
				for(size_t i=0; i<fragments.size(); i++){
					const Fragment& fragment = fragments[i];
					cout << fragment._fragmentIndex << "\t" << fragment._startPos << "\t" << fragment._endPos << "\t" << fragment._length << "\t" << (int) fragment._coverage << "\t" << fragment._readIndexes.size() << "\t" << fragment._nbBridgingReads.size() << endl; //<< "\t" << fragment._readIndexes.size() << endl; //<< "\t" << getCoverageInBounds(repeats, fragmentStartPos, fragmentEndPos) << endl;
				
				}
			}
			*/

			//getchar();
			vector<FragmentPath> paths;

			for(const Fragment& fragment : fragments){
				if(fragment._length * 1/_parent._minimizerDensity < 10000) continue;

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

			/*
			if(readIndex == 36203){
				for(size_t i=0; i<fragments.size(); i++){
					const Fragment& fragment = fragments[i];
					cout << fragment._fragmentIndex << "\t" << fragment._startPos << "\t" << fragment._endPos << "\t" << fragment._length << "\t" << (int) fragment._coverage << "\t" << fragment._readIndexes.size() << "\t" << fragment._nbBridgingReads.size() << "\t" << fragment._finalContigIndex << endl; //<< "\t" << fragment._readIndexes.size() << endl; //<< "\t" << getCoverageInBounds(repeats, fragmentStartPos, fragmentEndPos) << endl;
				}
			}
			*/

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

						vector<MinimizerType>::const_iterator last = readMinimizers.begin() + endPos + _parent._kminmerSize;
						vector<MinimizerType> contigMinimizers(first, last);


						//for(size_t i=0; i<contigMinimizers.size(); i++){
						//	cout << i <<" " << contigMinimizers[i] << endl;
						//}

						//const vector<MinimizerType>& contigMinimizers = readMinimizers.substr(startPos, endPos-startPos+_kminmzerSize-1);
						//cout << contigMinimizers.size() << endl;
						writeContig(contigMinimizers, isCircular);

						//writeRepeatContig(splitIndex, _contigHeaders[contigIndex], contigSeq, originalContigSeq);
						//splitIndex += 1;

					//}

					startPos = fragments[i]._startPos;
					currentContigIndex = fragments[i]._finalContigIndex;

				}

				//cout << fragment._fragmentIndex << "\t" << fragment._startPos << "\t" << fragment._endPos << "\t" << fragment._length << "\t" << fragment._coverage << "\t" << fragment._finalContigIndex << endl; //<< "\t" << getCoverageInBounds(repeats, fragmentStartPos, fragmentEndPos) << endl;

			}
			
			
			//getchar();
		}

		//unordered_map<std::pair<u_int32_t, u_int32_t>, ReadType> _fragmentPair_to_nbBridgingReads;

		void computeBridgingReads(const u_int32_t readIndex, vector<Fragment>& fragments){

			for(u_int32_t i=0; i<fragments.size(); i++){
				
				Fragment& fragment = fragments[i];

				//u_int32_t maxReadLength = 0;
				fragment._maxReadLength = 0;

				for(const ReadType& readIndex : fragment._readIndexes){

					const u_int32_t& readLength = _parent._readIndex_to_readLength[readIndex];

					if(readLength > fragment._maxReadLength){
						fragment._maxReadLength = readLength;	
					}
				} 

				//cout << fragment._fragmentIndex << ": " << fragment._maxReadLength << endl; 
			}

			for(u_int32_t i=0; i<fragments.size(); i++){

				Fragment& fragment1 = fragments[i];

				for(u_int32_t j=i+1; j<fragments.size(); j++){
					
					Fragment& fragment2 = fragments[j];

					u_int32_t bridgeLength = fragment2._startPos - fragment1._endPos;
					if(bridgeLength > fragment1._maxReadLength) break;

					int nbBridgingReads = getNbBridgingReads(fragment1, fragment2);

					fragment1._nbBridgingReads[j] = nbBridgingReads;
					fragment2._nbBridgingReads[i] = nbBridgingReads;
					
					//cout << i << " -> " << j << " " << bridgeLength << " " << nbBridgingReads << endl;
					//getchar();
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
		
		void writeContig(const vector<MinimizerType>& minimizers, const u_int8_t& isCircular){
	
			#pragma omp critical(writeContig)
			{

				_parent._nbContigs += 1;

				u_int32_t contigSize = minimizers.size();
				_parent._outputFile.write((const char*)&contigSize, sizeof(contigSize));
				_parent._outputFile.write((const char*)&isCircular, sizeof(isCircular));
				_parent._outputFile.write((const char*)&minimizers[0], contigSize*sizeof(MinimizerType));
			}
		}
		
	};


	
};	


#endif 


