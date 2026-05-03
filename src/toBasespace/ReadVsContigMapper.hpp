

#ifndef MDBG_METAG_READVSCONTIGMAPPER
#define MDBG_METAG_READVSCONTIGMAPPER

#include "../Commons.hpp"

class ReadVsContigMapper{
    
public:

	string _contigFilename;
	string _readFilename;
	string _outputFilename;
	float _averageDistanceBetweenMinimizers;
	int _nbCores;

	ofstream _alignmentOutputFile;
	//unordered_map<u_int128_t, ReadAbundance> _kminmer_to_readAbundance;

	struct Anchor{

		u_int32_t _referenceIndex;
		u_int32_t _referencePosition;
		u_int32_t _queryPosition;
		bool _isReversed;
		//int32_t _referencePositionIndex;
		//int32_t _queryPositionIndex;

		Anchor(const u_int32_t referenceIndex, const u_int32_t referencePosition, const u_int32_t queryPosition, const bool isReversed){
			_referenceIndex = referenceIndex;
			_referencePosition = referencePosition;
			_queryPosition = queryPosition;
			_isReversed = isReversed;
			//_referencePositionIndex = referencePositionIndex;
			//_queryPositionIndex = queryPositionIndex;
		}
	};


	struct Chain{
		float _chainingScore;
		vector<size_t> _anchorInterval;
		int64_t _nbMatches;
		int64_t _nbDifferencesQuery;
		int64_t _nbDifferencesReference;
		int64_t _queryStart;
		int64_t _queryEnd;
		int64_t _referenceStart;
		int64_t _referenceEnd;
		bool _isReversed;
	};


	struct MinimizerPosition{
		MinimizerPairType _minimizerPair;
		ReadType _readIndex;
		u_int32_t _position;
		bool _isReversed;

		friend bool operator<(const MinimizerPosition &a, const MinimizerPosition &b){
			return a._minimizerPair < b._minimizerPair;
		}
	};


	vector<MinimizerPosition> _allMinimizerPositions;
	ankerl::unordered_dense::map<MinimizerPairType, pair<u_int64_t, u_int64_t>> _minimizerLookupTable;


	ReadVsContigMapper(const string readFilename, const string contigFilename, const string outputFilename, const float averageDistanceBetweenMinimizers, const int nbCores) {
		_readFilename = readFilename;
		_contigFilename = contigFilename;
		_outputFilename = outputFilename;
		_averageDistanceBetweenMinimizers = averageDistanceBetweenMinimizers;
		_nbCores = nbCores;
	}

	void execute(){


		Logger::get().debug() << "\tLoading minimizers";
		loadContigs();
		//loadReadsTest();
		Logger::get().debug() << "\t\tDone: " << (peakrss() / 1024.0 / 1024.0 / 1024.0);

		//cout << _allMinimizerPositions.size() << endl;
		Logger::get().debug() << "\tNb minimizers: " << _allMinimizerPositions.size();


		//std::sort(_allMinimizerPositions.begin(), _allMinimizerPositions.end(), [](const MinimizerPosition& a, const MinimizerPosition& b){
		//	return a._minimizer < b._minimizer;
		//});
		Logger::get().debug() << "\tSorting minimizers";
		Commons::sortParallel(_allMinimizerPositions, _allMinimizerPositions.size(), _nbCores);
		//sortIndex();
		Logger::get().debug() << "\t\tDone: " << (peakrss() / 1024.0 / 1024.0 / 1024.0);

		//cout << _allMinimizerPositions.size() << endl;
		//cout << "todo: enelver les minimizer de la table des positions: _allMinimizerPositions" << endl;
		//cout << "Attention: if(anchors.size() < 30) return; a reduire a 3" << endl;

		Logger::get().debug() << "\tCreating minimizer lookup table";
		createLookupTable();
		Logger::get().debug() << "\t\tDone: " << (peakrss() / 1024.0 / 1024.0 / 1024.0);


		//cout << "Peak memory: " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << endl;

		Logger::get().debug() << "\tAligning reads";
		alignReads();
		Logger::get().debug() << "\tDone: " << (peakrss() / 1024.0 / 1024.0 / 1024.0);





		vector<MinimizerPosition>().swap(_allMinimizerPositions); //_allMinimizerPositions.clear();
		ankerl::unordered_dense::map<MinimizerPairType, pair<u_int64_t, u_int64_t>>().swap(_minimizerLookupTable); //_minimizerLookupTable.clear();
		


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

	}
	
	void loadContigs(){


		KminmerParserParallel parser(_contigFilename, 0, 2, false, false, 1);
		parser.parse(LoadContigsFunctor(*this));

		//KminmerParserParallel parser(_inputFilenameContig, _minimizerSize, _kminmerSize, false, false, 1);
		//parser.parseSequences(LoadContigsFunctor(*this));
		//exit(1);
	}

	class LoadContigsFunctor {

		public:

		ReadVsContigMapper& _parent;

		LoadContigsFunctor(ReadVsContigMapper& parent) : _parent(parent){
		}

		LoadContigsFunctor(const LoadContigsFunctor& copy) : _parent(copy._parent){
		}

		~LoadContigsFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {

			ReadType readIndex = kminmerList._readIndex;
			if(readIndex % 100000 == 0) Logger::get().debug() << "\t" << readIndex;

			//for(size_t i=0; i<kminmerList._readMinimizers.size(); i++){
			//	cout << i << "\t" << kminmerList._readMinimizers[i] << endl;
			//}

			/*
			//if(readIndex != 24937) return;

			for(u_int32_t i=0; i<kminmerList._readMinimizers.size(); i++){
				
				_parent._allMinimizerPositions.push_back({kminmerList._readMinimizers[i], readIndex, i});

				u_int64_t checksum = readIndex * kminmerList._readMinimizers.size() * kminmerList._readMinimizers[i] * (i+1);

				#pragma omp atomic
				_parent._checksum_contigs += checksum;
			}

			//_parent._contigSizes.push_back(kminmerList._readMinimizers.size());
			
			//cout << readIndex << " " << kminmerList._readMinimizers.size() << endl;
			*/

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
				bool isReversed = kminmerInfo._isReversed;
				u_int32_t minimizerPosition =  i;// (kminmerList._minimizerPos[i] + kminmerList._minimizerPos[i+1]) / 2;

				_parent._allMinimizerPositions.push_back({minimizerPair, readIndex, minimizerPosition, isReversed});
			}

		}
	};

	/*
	void sortIndex() {

		//cout << "Sort merge" << endl;
		const u_int64_t data_count = _allMinimizerPositions.size();
		auto get_block_edge = [data_count](int i, int n) {
			return data_count * i / n;
		};


		//cout << "sorting" << endl;

		u_int64_t blocks = _nbCores;
		#pragma omp parallel num_threads(_nbCores)
		{
			//blocks = omp_get_num_threads();
			u_int64_t block = omp_get_thread_num();
			u_int64_t start = get_block_edge(block, blocks);
			u_int64_t finish = get_block_edge(block + 1, blocks);
			std::sort(std::begin(_allMinimizerPositions) + start, std::begin(_allMinimizerPositions) + finish, [](const MinimizerPosition& a, const MinimizerPosition& b){
				return a._minimizerPair < b._minimizerPair;
			});
		}

		//cout << "merging" << endl;
		for (int merge_step = 1; merge_step < blocks; merge_step *= 2) {
			#pragma omp parallel for num_threads(_nbCores)
			for (u_int64_t i = 0; i < blocks; i += 2 * merge_step) {
				u_int64_t start = get_block_edge(i, blocks);
				u_int64_t mid = std::min(get_block_edge(i + merge_step, blocks), data_count);
				u_int64_t finish = std::min(get_block_edge(i + 2 * merge_step, blocks), data_count);
				if (mid < finish)
					std::inplace_merge(std::begin(_allMinimizerPositions) + start, std::begin(_allMinimizerPositions) + mid, std::begin(_allMinimizerPositions) + finish, [](const MinimizerPosition& a, const MinimizerPosition& b){
						return a._minimizerPair < b._minimizerPair;
					});
			}
		}
		
		//cout << "done" << endl;
	}
	*/

	void alignReads(){

		auto start = high_resolution_clock::now();
		
		_alignmentOutputFile = ofstream(_outputFilename);

		//cout << "otodo: skip correction mode" << endl;
		//bool hasQuality = true;
		//if(_useHomopolymerCompression){ //hifi
		//	hasQuality = false;
		//}
		/*
		if(_skipCorrection){
			KminmerParserParallel parser(_inputDir + "/read_data_init.txt", _minimizerSize, _kminmerSize, false, true, _nbCores);
			parser._densityThreshold = _minimizerDensity_assembly;
			//parser._densityThreshold = _minimizerDensity_assembly;
			parser.parseSequences(AlignReadFunctor(*this));
			//KminmerParserParallel parser(_inputDir + "/read_data_init.txt", _minimizerSize, _kminmerSize, false, false, _nbCores);
			//parser.parseSequences(AlignReadFunctor(*this));
		}
		else if(_dataType == SequencingPlatformType::HiFi){
			KminmerParserParallel parser(_inputDir + "/read_data_init.txt", _minimizerSize, _kminmerSize, false, true, _nbCores);
			parser._densityThreshold = _minimizerDensity_assembly;
			//parser._densityThreshold = _minimizerDensity_assembly;
			parser.parseSequences(AlignReadFunctor(*this));
		}
		else{
			//KminmerParserParallel parser(_inputDir + "/read_data_init.txt", _minimizerSize, _kminmerSize, false, true, _nbCores);
			//cout << "nb snpmer: " << _isSnpmer.size() << " " << "attention hard copy ici"<< endl;
			//parser._densityThreshold = _minimizerDensity_assembly;
			//parser._isSnpmer = _isSnpmer;
			//parser.parseSequences(AlignReadFunctor(*this));
			//cout << "single core here" << endl;
			KminmerParserParallel parser(_inputDir + "/read_data_init.txt", _minimizerSize, _kminmerSize, false, true, _nbCores);
			parser._densityThreshold = _minimizerDensity_assembly;
			parser.parseSequences(AlignReadFunctor(*this));

		}
		*/

		//cout << "single core here" << endl;
		//loadAlignmentsDebug();
		KminmerParserParallel parser(_readFilename, 0, 2, false, true, _nbCores);
		//parser._densityThreshold = _minimizerDensity_assembly;
		parser.parse(AlignReadFunctor(*this));

		_alignmentOutputFile.close();
		
		Logger::get().debug() << "\tDone " << " (" << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s) " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << "GB";
	}

	class AlignReadFunctor {

		public:

		ReadVsContigMapper& _parent;
		MinimizerAligner* _minimizerAligner;
		ReadMapping2 _bestReadMapping;
		//std::unique_ptr<spoa64::AlignmentEngine> _alignmentEngine;// = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -4);
		//std::unique_ptr<spoa64::AlignmentEngine> _alignmentEngineTrimming;

		AlignReadFunctor(ReadVsContigMapper& parent) : _parent(parent){
			_minimizerAligner = nullptr;
		}

		AlignReadFunctor(const AlignReadFunctor& copy) : _parent(copy._parent){
			//_alignmentEngine = spoa64::AlignmentEngine::Create(spoa64::AlignmentType::kNW, 3, -5, -4);
			//_alignmentEngineTrimming = spoa64::AlignmentEngine::Create(spoa64::AlignmentType::kSW, 3, -5, -4);
			//_alignmentEngine = spoa64::AlignmentEngine::Create(spoa64::AlignmentType::kSW, 3, -2, -1);
			_minimizerAligner = new MinimizerAligner(3, -1, -1, -1);
		}

		~AlignReadFunctor(){
			if(_minimizerAligner != nullptr) delete _minimizerAligner;
		}



		void operator () (const KminmerList& kminmerList) {
			
			if(kminmerList._readIndex % 100000 == 0){
				Logger::get().debug() << "\t\tAligning reads: " << kminmerList._readIndex;
			}

			MinimizerRead read = {kminmerList._readIndex, kminmerList._readMinimizers, kminmerList._minimizerPos, kminmerList._readQualities, kminmerList._readMinimizerDirections, kminmerList._readLength};
			
			_bestReadMapping._contigIndex = -1;
			//const MinimizerRead readHighDensity = read2.toMinimizerRead(_parent._minimizerSize, _parent._snpmerSize);
			//const MinimizerRead read = Utils::getLowDensityMinimizerRead(readHighDensity, _parent._minimizerDensity_assembly);
			
			/*
			if(read._readIndex == 914057){

				char* dnaStr1 = read2->to_string();
				string read1 = string(dnaStr1, strlen(dnaStr1));
				free(dnaStr1);

				cout << "lala: " << read1.size() << " " << read._minimizers.size() << endl;


			}
			return;

			read2.destroy();
			*/
			/*
			//_parent.alignRead(read, _minimizerAligner);
			mapRead(read, false);

			const MinimizerRead& read_rc = read.toReverseComplement();
			mapRead(read_rc, true);

			writeAlignment2();
			

			*/
		

			mapRead2(kminmerList, read);

			//const ReadMapping2& bestReadMapping = _parent._debug_readIndex_to_bestAlignment[read._readIndex];
			//cout << bestReadMapping._contigIndex << "\t" << bestReadMapping._contigStart << "\t" << bestReadMapping._contigEnd << "        " << bestReadMapping._readStart << "\t" << bestReadMapping._readEnd << "\t" << bestReadMapping._isReversed << endl;

			writeAlignment2();
			
		}

		void writeAlignment2(){


			if(_bestReadMapping._contigIndex == -1){
				//cout << "\tNone" << endl;
				//getchar();
				return;
			}


			//if(_bestReadMapping._readIndex == 18731 || _bestReadMapping._readIndex == 21736){
			//	cout << _bestReadMapping._readIndex << "\t" << _bestReadMapping._contigIndex << "\t" << _bestReadMapping._contigStart << "\t" << _bestReadMapping._contigEnd << "        " << _bestReadMapping._readStart << "\t" << _bestReadMapping._readEnd << "\t" << _bestReadMapping._isReversed << "\t" << _bestReadMapping._readStartReal << "\t" << _bestReadMapping._readEndReal << endl;

			//}

			//cout << _bestReadMapping._contigIndex << "\t" << _bestReadMapping._contigStart << "\t" << _bestReadMapping._contigEnd << "        " << _bestReadMapping._readStart << "\t" << _bestReadMapping._readEnd << "\t" << _bestReadMapping._isReversed << endl;

			#pragma omp critical(indexMatchingPosition)
			{
				_bestReadMapping.write(_parent._alignmentOutputFile);
			}

			//getchar();
		}


		void mapRead2(const KminmerList& kminmerList, const MinimizerRead& read){

			//cout << "Align: " << kminmerList._readIndex << endl;

			vector<Anchor> anchors;

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
				bool queryIsReversed = kminmerInfo._isReversed;
				u_int32_t queryPosition = i; //(kminmerList._minimizerPos[i] + kminmerList._minimizerPos[i+1]) / 2;

				//for(u_int32_t i=0; i<read._minimizers.size(); i++){


				//const MinimizerType minimizer = read._minimizers[i];
				//const u_int32_t queryPosition = read._minimizersPos[i];
				//const bool queryIsReversed = read._readMinimizerDirections[i];

				if(_parent._minimizerLookupTable.find(minimizerPair) == _parent._minimizerLookupTable.end()) continue;

				const auto& range = _parent._minimizerLookupTable[minimizerPair];

				for(size_t j=range.first; j<=range.second; j++){
					//const MinimizerPosition2& mPos = _parent._allMinimizerPositions[j];

					const MinimizerPosition& mPos = _parent._allMinimizerPositions[j];
					//anchors.push_back({mPos._contigIndex, mPos._position, i});

					//if(kminmerList._readIndex == mPos._readIndex) continue; //query == reference

					//Anchor2(const ReadType readIndex, const int64_t referencePosition, const int64_t queryPosition, const bool isReversed){
					//anchors.push_back({referenceMinimizerPosition._position, queryPosition, referenceMinimizerPosition._isReversed != queryIsReversed, referenceMinimizerPosition._positionIndex, i});
					anchors.push_back({mPos._readIndex, mPos._position, queryPosition, mPos._isReversed != queryIsReversed});
					//_loul += 1;

					//if(mPos._position == 2600) print = true;
				}
			}

			std::sort(anchors.begin(), anchors.end(), [](const Anchor& a, const Anchor& b){
				if(a._referenceIndex == b._referenceIndex){

					if(a._referencePosition == b._referencePosition){
						return a._queryPosition < b._queryPosition;
					}
					return a._referencePosition < b._referencePosition;
				}
				
				return a._referenceIndex < b._referenceIndex;

			});

			/*
			//if(print){
			cout << endl << endl << "Read: " << read._readIndex << endl;
			//unordered_set<u_int32_t> contigIndexes;

			for(size_t i=0; i<read._minimizers.size(); i++){
				cout << i << "\t" << read._minimizers[i] << endl;
			}

			for(size_t i=0; i<anchors.size(); i++){
				const Anchor& anchor = anchors[i];
				cout << i << "\t" << anchor._referenceIndex << "\t" << anchor._referencePosition << "\t" << anchor._queryPosition << "\t" << anchor._isReversed << "\t" << read._minimizersPos[anchor._queryPosition] << endl;
				//contigIndexes.insert(anchor._referenceIndex);
			}
			//}
			*/
			

			u_int32_t nbChainings = 0;
			u_int32_t referenceIndex = -1;
			vector<Anchor> subAnchors;

			for(size_t i=0; i<anchors.size(); i++){

				if(referenceIndex != anchors[i]._referenceIndex){

					if(referenceIndex != -1){
						processAnchors(subAnchors, read, referenceIndex);
						subAnchors.clear();
						nbChainings += 1;
					}

					referenceIndex = anchors[i]._referenceIndex;
				}

				subAnchors.push_back(anchors[i]);


			}

			if(subAnchors.size() > 0){
				processAnchors(subAnchors, read, referenceIndex);
				nbChainings += 1;
			}



		}
		
		/*
		void mapRead(const MinimizerRead& read, const bool isReversed){

			//bool print = read._readIndex == 2479658;


			
			ReadType readIndex = read._readIndex;

			//if(readIndex % 10000 == 0) Logger::get().debug() << readIndex;
			//cout << readIndex << endl;

			vector<MinimizerType> readMinimizers = read._minimizers;

			vector<Anchor> anchors;

			for(u_int32_t i=0; i<read._minimizers.size(); i++){
				const MinimizerType minimizer = read._minimizers[i];

				if(_parent._minimizerLookupTable.find(minimizer) == _parent._minimizerLookupTable.end()) continue;

				const auto& range = _parent._minimizerLookupTable[minimizer];

				for(size_t j=range.first; j<=range.second; j++){
					const MinimizerPosition& mPos = _parent._allMinimizerPositions[j];
					anchors.push_back({mPos._contigIndex, mPos._position, i});
					//_loul += 1;

					//if(mPos._position == 2600) print = true;
				}
			}

			std::sort(anchors.begin(), anchors.end(), [](const Anchor& a, const Anchor& b){
				if(a._referenceIndex == b._referenceIndex){

					if(a._referencePosition == b._referencePosition){
						return a._queryPosition < b._queryPosition;
					}
					return a._referencePosition < b._referencePosition;
				}
				
				return a._referenceIndex < b._referenceIndex;

			});


			u_int32_t nbChainings = 0;
			u_int32_t referenceIndex = -1;
			vector<Anchor> subAnchors;

			for(size_t i=0; i<anchors.size(); i++){

				if(referenceIndex != anchors[i]._referenceIndex){

					if(referenceIndex != -1){
						processAnchors(subAnchors, read, referenceIndex, isReversed);
						subAnchors.clear();
						nbChainings += 1;
					}

					referenceIndex = anchors[i]._referenceIndex;
				}

				subAnchors.push_back(anchors[i]);


			}

			if(subAnchors.size() > 0){
				processAnchors(subAnchors, read, referenceIndex, isReversed);
				nbChainings += 1;
			}

			//cout << nbChainings << " " << contigIndexes.size() << endl;
			//if(nbChainings != contigIndexes.size()){
			//	cout << "derp" << endl;
			//	getchar();
			//}
		}
		*/

		void processAnchors(const vector<Anchor>& anchors, const MinimizerRead& queryRead, const u_int32_t contigIndex){
			if(anchors.size() < 2) return;

			
			
			//if(_bestReadMapping._readIndex == 18731 || _bestReadMapping._readIndex == 21736){
			//	cout << _bestReadMapping._readIndex << "\t" << _bestReadMapping._contigIndex << "\t" << _bestReadMapping._contigStart << "\t" << _bestReadMapping._contigEnd << "        " << _bestReadMapping._readStart << "\t" << _bestReadMapping._readEnd << "\t" << _bestReadMapping._isReversed << "\t" << _bestReadMapping._readStartReal << "\t" << _bestReadMapping._readEndReal << endl;

			//}
			/*
			if(queryRead._readIndex == 18731){
				cout << endl;
				for(size_t i=0; i<anchors.size(); i++){
					const Anchor& anchor = anchors[i];
					cout << i << "\t" << anchor._referenceIndex << "\t" << anchor._referencePosition << "\t" << anchor._queryPosition << "\t" << queryRead._minimizersPos[anchor._queryPosition] << endl;

				}
			}
			*/

			const Chain& chain = chainAnchors(anchors, 10, 20, queryRead);

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

			writeAlignment(chain, queryRead, contigIndex);
			//getchar();
		}



		Chain chainAnchors(const vector<Anchor>& anchors, const int64_t& maxChainingBand, const float w, const MinimizerRead& queryRead) {

			bool print = false;// queryRead._readIndex == 8195; //r1 = 3762
			//vector<Chain> chainIntervals;

			if(print){

				cout << endl;
				for(size_t i=0; i<queryRead._minimizers.size(); i++){
					cout << i << "\t" << queryRead._minimizers[i] << "\t" << queryRead._minimizersPos[i] << endl;
				}
				cout << endl;
				cout << queryRead._readIndex << endl;
				for(size_t i=0; i<anchors.size(); i++){
					const Anchor& anchor = anchors[i];
					cout << i << "\t" << anchor._referenceIndex << "\t" << anchor._referencePosition << "\t" << anchor._queryPosition << "\t" << queryRead._minimizersPos[anchor._queryPosition] << endl;

				}
			}

			vector<float> scores(anchors.size(), 0);
			vector<size_t> parents(anchors.size(), 0);

			for (size_t i=0; i<anchors.size(); i++) {


				float bestScore = 0;
				size_t bestPrevIndex = i;

				argmaxPosition(anchors, scores, i, w, bestScore, bestPrevIndex, maxChainingBand, queryRead);

				
				if(bestPrevIndex != i) {
					scores[i] = bestScore;
					parents[i] = bestPrevIndex;
				}
				else{
					scores[i] = w;
					parents[i] = -1;
				}
			

			}


			
			float maxScore = 0;
			size_t bestIndex = -1;

			for(size_t i=0; i<anchors.size(); i++) {
				if(scores[i] > maxScore){
					maxScore = scores[i];
					bestIndex = i;
				}
			}

			if(print){
				cout << "\tParents:" << endl;
				for(size_t i=0; i<anchors.size(); i++){
					cout << "\t\t" << i << "\t" << parents[i] << "\t" << scores[i] << endl;
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

			Chain chain;
			chain._chainingScore = 0;

			if(interval.size() < 2) return chain;

			std::reverse(interval.begin(), interval.end());

			chain._chainingScore = maxScore;
			chain._anchorInterval = interval;
			chain._nbMatches = interval.size();

			const Anchor& firstAnchor = anchors[interval[0]];
			const Anchor& lastAnchor = anchors[interval[interval.size()-1]];

			chain._isReversed = firstAnchor._queryPosition > lastAnchor._queryPosition;
			
			//cout << firstAnchor._referenceIndex << " " << maxScore << " " << chain._isReversed << endl;

			if(chain._isReversed){
				chain._nbDifferencesQuery = (firstAnchor._queryPosition - lastAnchor._queryPosition + 1) - chain._nbMatches;
				chain._queryStart = lastAnchor._queryPosition;
				chain._queryEnd = firstAnchor._queryPosition + 1; //+1 because we use pair for indexing;
			}
			else{
				chain._nbDifferencesQuery = (lastAnchor._queryPosition - firstAnchor._queryPosition + 1) - chain._nbMatches;
				chain._queryStart = firstAnchor._queryPosition;
				chain._queryEnd = lastAnchor._queryPosition + 1; //+1 because we use pair for indexing;
			}


			chain._nbDifferencesReference = (lastAnchor._referencePosition - firstAnchor._referencePosition + 1) - chain._nbMatches;
			chain._referenceStart = firstAnchor._referencePosition;
			chain._referenceEnd = lastAnchor._referencePosition + 1; //+1 because we use pair for indexing
			


			return chain;
			
		}


		void argmaxPosition(const vector<Anchor>& anchors, vector<float>& scores, size_t i, float w, float& bestScore, size_t& bestPrevIndex, const int64_t& maxChainingBand, const MinimizerRead& queryRead) {
			
			bestScore = 0;
			bestPrevIndex = i;

			const Anchor& xi = anchors[i];
			//console.log("\t", xi);
			for (int64_t j = i-1; j >= 0; j--) {
		
				if(i-j > maxChainingBand) break;
				const Anchor& xj = anchors[j];

				
				if(xi._referencePosition == xj._referencePosition || xi._queryPosition == xj._queryPosition) continue;
				if(xi._isReversed != xj._isReversed) continue;
				//cout << "\t" << i << " " << k << " " << xi._isReversed << " " << xj._isReversed << " " << xi._referencePosition << "-" << xi._queryPosition << " " << xj._referencePosition << "-" << xj._queryPosition << endl;

				//int32_t anchorSpacing = queryRead._minimizersPos[xi._queryPosition] - queryRead._minimizersPos[xj._queryPosition];
				//if(anchorSpacing > 5000) continue;

				int32_t d_q;// = abs(xi._queryPosition - xj._queryPosition);
				if(xi._isReversed){
					d_q = xj._queryPosition - xi._queryPosition;
				}
				else{
					d_q = xi._queryPosition - xj._queryPosition;
				}
				int32_t d_r = xi._referencePosition - xj._referencePosition;

				int32_t referenceSpacingApprox = (xi._referencePosition*_parent._averageDistanceBetweenMinimizers) - (xj._referencePosition*_parent._averageDistanceBetweenMinimizers);
				if(referenceSpacingApprox > 5000) continue;
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
				//if(d_q > 5000 || d_r > 5000) {
				//	continue;
				//}

				if(d_r <= 0) {
					continue;
				}

				int32_t gap = abs(d_r - d_q);
				//cout << "\t" << xi._referencePosition << "\t" << xj._referencePosition << "\t" << gap << endl;

				if(gap > 100) {
					continue;
				}
				
				if(xi._isReversed){
					if(queryRead._minimizersPos[xj._queryPosition] - queryRead._minimizersPos[xi._queryPosition] > 5000) continue;
					if(xi._queryPosition > xj._queryPosition) continue;
				}
				else{
					if(queryRead._minimizersPos[xi._queryPosition] - queryRead._minimizersPos[xj._queryPosition] > 5000) continue;
					if(xi._queryPosition < xj._queryPosition) continue;
				}

				//if(xi._isReversed){
				//	if(xi._queryPosition > xj._queryPosition) continue;
				//}
				//else{
				//if(xi._queryPosition < xj._queryPosition) continue;
				//}
				
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

		void writeAlignment(const Chain& chain, const MinimizerRead& queryRead, const u_int32_t contigIndex){

			//if(isReversed) return;
			//const Anchor& firstAnchor = chain
			//u_int32_t nbDifferencesQuery = 
			ReadType readIndex = queryRead._readIndex;
			u_int32_t readStart = chain._queryStart;
			u_int32_t readEnd = chain._queryEnd;
			u_int32_t contigStart = chain._referenceStart;
			u_int32_t contigEnd = chain._referenceEnd;
			u_int32_t nbMatches = chain._nbMatches;
			u_int32_t nbDifferences = chain._nbDifferencesQuery + chain._nbDifferencesReference;
			u_int32_t readLengthBp = queryRead._readLength;
			//u_int32_t readSize = queryRead._minimizers.size();
			//u_int32_t contigSize = _parent._contigSizes[contigIndex];

			int32_t overhangStart = queryRead._minimizersPos[readStart] / _parent._averageDistanceBetweenMinimizers;
			int32_t overhangEnd = (readLengthBp-queryRead._minimizersPos[readEnd]) / _parent._averageDistanceBetweenMinimizers;

			int32_t matchScore = nbMatches  - overhangStart - overhangEnd; //- nbDifferences
			bool isReversed = chain._isReversed;

			//if(matchScore <= 0) return;
			//cout << contigIndex << " " << matchScore << endl;
			//bool print = queryRead._readIndex == 2479658;
			//if(print){
			//	cout << matchScore << endl;
			//	exit(1);
			//}



			//if(matchScore <= 0) return;

			ReadMapping2 readMapping = {readIndex, contigIndex, readStart, readEnd, contigStart, contigEnd, isReversed, matchScore, queryRead._minimizersPos[readStart], queryRead._minimizersPos[readEnd], queryRead._readLength};
			
			//if(isReversed){
				//return;
			//	std::swap(readMapping._contigStart, readMapping._contigEnd);
			//}

			//if(print ){ //190058
			//	cout << contigIndex << "     " << readStart << "-" << readEnd << "    " << contigStart << "-" << contigEnd << "    " << readSize << "    " << nbMatches << "-" << nbDifferences << " " << queryRead._minimizersPos[readStart] << " " << minimizersPos[readEnd] << " " << overhangStart << " " << overhangEnd << endl;
			//}

			//if(isReversed){
			//	cout << contigStart << " " << contigEnd << endl;
			//}

			//getchar();
			//if(readMapping._readStartReal > 1000) return;
			//if(readLengthBp-readMapping._readEndReal > 1000) return;
			/*
			if(readEnd+2 < readMinimizers.size()){

				int n = 0;
				for(size_t i=0; i<readMinimizers.size(); i++){
					for(size_t j=0; j<contigMinimizers.size(); j++){
						if(readMinimizers[i] == contigMinimizers[j]){
							cout << "Match: " << i << " " << j << endl;
							n += 1;
						}
					}
				}

				cout << isReversed << endl;
				cout << contigEndLocal << endl;
				cout << readMinimizers[readEnd] << " " << contigMinimizers[contigEndLocal] << endl;
				cout << readMinimizers[readEnd+1] << " " << contigMinimizers[contigEndLocal+1] << endl;
				cout << readMinimizers[readEnd+2] << " " << contigMinimizers[contigEndLocal+2] << endl;

				//for(size_t i=0; i<contigMinimizers.size(); i++){
				//	cout << i << "\t" << contigMinimizers[i] << endl;
				//}
				cout << readStart << " " << readEnd << " " << readSize << " " << n << endl;
				cout << readMapping.isMaximalMapping(2) << endl;
				getchar();
			}
			*/

			//if(readIndex == 667){
			//	cout << readMapping._readStartReal << " " << readMapping._readEndReal << " " << readLengthBp << endl;
			//	getchar();
			//}

			//if(hasLargeIndel){
			//	cout << readIndex << "\t" << readMapping._contigStart << "\t" << readMapping._contigEnd << endl;
			//}

			//if(hasLargeIndel) return;
			
			//if(readIndex == 2724062){
			//	cout << "lala: " << readIndex << " " << readMapping._readStart << "-" << readMapping._readEnd << "\t" << readMapping._contigStart << "-" << readMapping._contigEnd << "\t" << nbMatches << " " << nbDifferences << " " << overhangStart << " " << overhangEnd << " " << queryRead._minimizersPos[readStart] << " " << (readLengthBp-queryRead._minimizersPos[readEnd]) << endl;
			//}
			//if(readIndex == 2701596){
			//	cout << "lala: " << readIndex << " " << readMapping._readStart << "-" << readMapping._readEnd << "\t" << readMapping._contigStart << "-" << readMapping._contigEnd << "\t" << nbMatches << " " << nbDifferences << " " << overhangStart << " " << overhangEnd << " " << queryRead._minimizersPos[readStart] << " " << (readLengthBp-queryRead._minimizersPos[readEnd]) << endl;
			//}

			//if(!readMapping.isMaximalMapping(2)) return;

			/*
			#pragma omp critical(indexMatchingPosition)
			{
				
				if(_parent._readIndex_to_bestContigMap.find(readIndex) == _parent._readIndex_to_bestContigMap.end()){
					//cout << "Add 1" << endl;
					//cout << readMapping._readStart << "\t" << readMapping._readEnd << "\t" << readMapping._contigStart << "\t" << readMapping._contigEnd << "\t" << readMapping._isReversed << endl;
					_parent._readIndex_to_bestContigMap[readIndex] = readMapping;
				}
				else{
					
					//cout << "Add 2" << endl;
					//exit(1);
					const ReadMapping2& existingReadMapping = _parent._readIndex_to_bestContigMap[readIndex];

					if(readMapping._matchScore > existingReadMapping._matchScore){
						_parent._readIndex_to_bestContigMap[readIndex] = readMapping;
					}
				}
			}
			*/

			if(_bestReadMapping._contigIndex == -1){
				_bestReadMapping = readMapping;	
			}
			else{
				if(readMapping._matchScore == _bestReadMapping._matchScore){ //Arbitrary resolve tie
					if(readMapping._readIndex > _bestReadMapping._readIndex){
						_bestReadMapping = readMapping;
					}
					else if(readMapping._contigStart < _bestReadMapping._contigStart){
						_bestReadMapping = readMapping;
					}
				}
				else if(readMapping._matchScore > _bestReadMapping._matchScore){
					_bestReadMapping = readMapping;
				}
			}
			
			//if(readIndex == 51750){
			//	cout << readMapping._readIndex << "\t" << readMapping._contigIndex << "\t" << readMapping._matchScore << endl;
			//}
		}

	};


};	

#endif 


