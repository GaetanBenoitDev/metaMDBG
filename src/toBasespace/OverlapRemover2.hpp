



#ifndef MDBG_METAG_OverlapRemover2
#define MDBG_METAG_OverlapRemover2

#include "../Commons.hpp"


class OverlapRemover2 {
  

public:


	struct ContigMinimizerPosition{

		u_int32_t _contigIndex;
		u_int32_t _positionIndex;
	};

	struct Anchor{
		u_int32_t _referencePositionIndex;
		u_int32_t _queryPositionIndex;
	};

	struct Chain{
		float _chainingScore;
		vector<size_t> _anchorInterval;
	};

	typedef phmap::parallel_flat_hash_map<KmerVecHashType, vector<ContigMinimizerPosition>, phmap::priv::hash_default_hash<KmerVecHashType>, phmap::priv::hash_default_eq<KmerVecHashType>, std::allocator<std::pair<KmerVecHashType, vector<ContigMinimizerPosition>>>, 10, std::mutex> KminmerPosMap;
	//typedef phmap::parallel_flat_hash_map<u_int32_t, u_int32_t, phmap::priv::hash_default_hash<u_int32_t>, phmap::priv::hash_default_eq<u_int32_t>, std::allocator<std::pair<u_int32_t, u_int32_t>>, 10, std::mutex> ContigSizeMap;


	KminmerPosMap _minimizerIndexer;
	//ContigSizeMap _contigLengths;

	string _inputDir;
	string _inputFilenameContig;
	string _outputFilenameContig;
	size_t _kminmerSize;
	int _nbCores;
	
	ofstream _outputFile;
	//u_int64_t _nbDuplicatedContigs;
	//u_int64_t _nbContigs;
	//unordered_set<u_int32_t> _isContigDuplicated;

	//struct ContigData{
	//	u_int32_t _size;
	//	u_int32_t _overlapLeft;
	//	u_int32_t _overlapRight;
	//};

	//vector<ContigData> _contigOverlaps;
	vector<u_int32_t> _contigSizes;

	OverlapRemover2(const string& inputDir, const string& inputFilenameContig, const string& outputFilenameContig, size_t kminmerSize, const int nbCores){
		_inputDir = inputDir;
		_inputFilenameContig = inputFilenameContig;
		_outputFilenameContig = outputFilenameContig;
		_kminmerSize = kminmerSize-1;
		_nbCores = nbCores;
	}

	

	void execute(){

		//cout << "todo: remove overlap self, attention circular" << endl;
		//_nbContigs = 0;

		Logger::get().debug() << "Count contigs";
		countContigs();

		//cout << "Nb contigs: " << _nbContigs << endl;
		//_contigOverlaps.resize(_nbContigs, {});


		Logger::get().debug() << "Indexing contigs";
		indexContigs();


		Logger::get().debug() << "Map contigs";
		mapContigs();
		
		//Logger::get().debug() << "Write contigs";
		//writeContigs();

		//cout << "Nb duplicated contigs: " << _nbDuplicatedContigs << endl;
	}
	
	void countContigs(){

		KminmerParserParallel parser(_inputFilenameContig, 0, 0, false, false, 1);
		parser.parseSequences(CountContigsFunctor(*this));
	}


	class CountContigsFunctor {

		public:

		OverlapRemover2& _parent;

		CountContigsFunctor(OverlapRemover2& parent) : _parent(parent){
		}

		CountContigsFunctor(const CountContigsFunctor& copy) : _parent(copy._parent){
		}

		~CountContigsFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {

			_parent._contigSizes.push_back(kminmerList._readMinimizers.size());
			//#pragma atomic
			//_parent._nbContigs += 1;

		}
	};
	

	void indexContigs(){

		KminmerParserParallel parser(_inputFilenameContig, 0, _kminmerSize, false, false, _nbCores);
		parser.parse(IndexContigsFunctor(*this));
	}

	class IndexContigsFunctor {

		public:

		OverlapRemover2& _parent;

		IndexContigsFunctor(OverlapRemover2& parent) : _parent(parent){
		}

		IndexContigsFunctor(const IndexContigsFunctor& copy) : _parent(copy._parent){
		}

		~IndexContigsFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {


			u_int32_t contigIndex = kminmerList._readIndex;
			

			if(contigIndex % 100000 == 0) Logger::get().debug() << "Index contigs: " << contigIndex;


			u_int32_t contigSize = kminmerList._kminmersInfo.size(); //kminmerList._readMinimizers.size();

			//for(size_t i=0; i<kminmerList._readMinimizers.size(); i++){

			for(u_int32_t i=0; i<kminmerList._kminmersInfo.size(); i++){
			
				
				//cout << "------- " << i << endl;

				//_logFile << readIndex << " " << i << endl;
				const ReadKminmerComplete& kminmerInfo = kminmerList._kminmersInfo[i];

				const KmerVec& vec = kminmerInfo._vec;
				KmerVecHashType vecHash = vec.hash128();
				//const MinimizerType minimizer = kminmerList._readMinimizers[i];
			
				
				//const ReadKminmerComplete& kminmerInfo = kminmerList._kminmersInfo[i];

				//const KmerVec& vec = kminmerInfo._vec;

				u_int32_t minimizerPositionIndex = i;
				ContigMinimizerPosition contigPos = {contigIndex, minimizerPositionIndex};

				_parent._minimizerIndexer.lazy_emplace_l(vecHash,
				[&contigPos](KminmerPosMap::value_type& v) { // key exist
					v.second.push_back(contigPos);
					//v.second[0] += 1; //Increment kminmer count
					//if(std::find(v.second.begin(), v.second.end(), readIndex) == v.second.end()){
					//	v.second.push_back(readIndex);
					//}
				},           
				[&vecHash, &contigPos](const KminmerPosMap::constructor& ctor) { // key inserted
					
					vector<ContigMinimizerPosition> readIndexes = {contigPos}; 

					ctor(vecHash, readIndexes); 

				}); // construct value_type in place when key not present
				



			}

			//_parent._contigOverlaps[contigIndex] = {contigSize, 0, 0};

		}
	};
	
	void mapContigs(){

		_outputFile.open(_outputFilenameContig);

		//cout << "single core here" << endl;
		KminmerParserParallel parser(_inputFilenameContig, 0, _kminmerSize, false, false, _nbCores);
		parser.parseSequences(MapContigsFunctor(*this));
		
		_outputFile.close();
	}

	class MapContigsFunctor {

		public:

		OverlapRemover2& _parent;

		MapContigsFunctor(OverlapRemover2& parent) : _parent(parent){
		}

		MapContigsFunctor(const MapContigsFunctor& copy) : _parent(copy._parent){
		}

		~MapContigsFunctor(){
		}



		void operator () (const KminmerList& kminmerList) {

			if(kminmerList._readIndex % 100000 == 0) Logger::get().debug() << "Map contigs: " << kminmerList._readIndex;

			vector<MinimizerType> minimizers = kminmerList._readMinimizers;
			//pair<u_int32_t, u_int32_t> contigOverlapResult = {0, 0};// = _parent._contigOverlaps[kminmerList._readIndex];

			//cout << kminmerList._readIndex << " " << minimizers.size() << endl;

			while(true){

				const pair<u_int32_t, u_int32_t> contigOverlap = computeOverlaps(kminmerList._readIndex, minimizers);
				if(contigOverlap.first == 0 && contigOverlap.second == 0) break;

				//cout << contigOverlap.first << " " << contigOverlap.second << endl;

				u_int32_t overlapLeft = 0;
				if(contigOverlap.first > 0){
					overlapLeft = contigOverlap.first + _parent._kminmerSize - 1;
				}

				u_int32_t overlapRight = 0;
				if(contigOverlap.second > 0){
					overlapRight = contigOverlap.second + _parent._kminmerSize - 1;
				}

				int64_t indexEnd = ((int64_t)minimizers.size()) - ((int64_t)overlapRight);

				if(overlapLeft + overlapRight >= minimizers.size()) return;
				if(overlapLeft >= indexEnd) return;
				
				vector<MinimizerType> newMinimizers(minimizers.begin() + overlapLeft, minimizers.begin() + indexEnd);

				if(newMinimizers.size() <= _parent._kminmerSize+1) return;

				minimizers = newMinimizers;

				//if(kminmerList._readIndex == 92635){
				//	for(size_t i=0; i<minimizers.size(); i++){
				//		cout << i << "\t" << minimizers[i] << endl;
				//	}
				//}
				//cout << "check: " << minimizers.size() << endl;
				//getchar();

			}

			//cout << "Write: " << minimizers.size() << endl;
			writeContig(minimizers, kminmerList._isCircular);

		}

		pair<u_int32_t, u_int32_t> computeOverlaps(const u_int32_t referenceContigIndex, const vector<MinimizerType>& minimizers){

			//const u_int32_t referenceContigIndex = kminmerList._readIndex;
			const u_int32_t referenceContigLength = minimizers.size();
			//ContigData& contigOverlap = _parent._contigOverlaps[referenceContigIndex];
			pair<u_int32_t, u_int32_t> contigOverlap;

			//cout << endl << "Contig index: " << referenceContigIndex << "\t" << kminmerList._kminmersInfo.size() << endl;

			unordered_map<u_int32_t, vector<Anchor>> _contigIndex_to_anchors;
			//unordered_map<u_int32_t, u_int32_t> _readIndex_to_nbMatches;

			vector<u_int64_t> rlePositions;
			vector<ReadKminmerComplete> kminmersInfo;
			vector<u_int32_t> minimizerPos(minimizers.size(), 0);
			vector<u_int8_t> minimizerQualities(minimizers.size(), 0);
			//MDBG::getKminmers(_l, _k, minimizers, minimizersPos, kminmers, kminmersInfo, rlePositions, 0, false);
			MDBG::getKminmers_complete(_parent._kminmerSize, minimizers, minimizerPos, kminmersInfo, referenceContigIndex, minimizerQualities);


			//for(size_t i=0; i<kminmerList._readMinimizers.size(); i++){
			for(u_int32_t i=0; i<kminmersInfo.size(); i++){
			
				
				//cout << "------- " << i << endl;

				//_logFile << readIndex << " " << i << endl;
				const ReadKminmerComplete& kminmerInfo = kminmersInfo[i];

				const KmerVec& vec = kminmerInfo._vec;
				KmerVecHashType vecHash = vec.hash128();

				//const MinimizerType minimizer = kminmerList._readMinimizers[i];
			
				
				//const ReadKminmerComplete& kminmerInfo = kminmerList._kminmersInfo[i];

				//const KmerVec& vec = kminmerInfo._vec;


				u_int32_t referenceMinimizerPositionIndex = i;
				//ContigMinimizerPosition contigPos = {contigIndex, minimizerPositionIndex};

				if(_parent._minimizerIndexer.find(vecHash) == _parent._minimizerIndexer.end()) continue;

				for(const ContigMinimizerPosition& contigPos : _parent._minimizerIndexer[vecHash]){

					const u_int32_t queryContigIndex = contigPos._contigIndex;
					const u_int32_t queryContigLength = _parent._contigSizes[queryContigIndex];

					if(queryContigIndex == referenceContigIndex) continue;
					if(queryContigLength < referenceContigLength) continue;

					const Anchor anchor = {referenceMinimizerPositionIndex, contigPos._positionIndex};
					_contigIndex_to_anchors[contigPos._contigIndex].push_back(anchor);
					//_readIndex_to_nbMatches[contigPos._contigIndex] += 1;
				}
				



			}

			//int max = 0;
			//for(const auto& it : _readIndex_to_nbMatches){
			//	if(it.second > max){
			//		max = it.second;
			//	}
			//}

			//cout << max << endl;
			
			//cout << _contigIndex_to_anchors.size() << endl;
			for(auto& it : _contigIndex_to_anchors){

				const u_int32_t queryContigIndex = it.first;
				vector<Anchor> anchors = it.second;

				if(anchors.size() == 0) continue;

				//if(referenceContigIndex == 234){
				//	cout << "\t" << queryContigIndex << endl;
				//	for(size_t i=0; i<anchors.size(); i++){
				//		cout << "\t\t" << i << "\t" << anchors[i]._queryPositionIndex << "\t" << anchors[i]._referencePositionIndex << endl;
				//	}
				//}
				//cout << "\t" << queryContigIndex << "\t" << anchors.size() << endl;
				//if(referenceContigIndex == 4 && queryContigIndex == 2) {

				//}
				//else{
				//	continue;
				//}

				std::sort(anchors.begin(), anchors.end(), [](const Anchor& a, const Anchor& b){
					if(a._referencePositionIndex == b._referencePositionIndex){
						return a._queryPositionIndex < b._queryPositionIndex;
					}
					return a._referencePositionIndex < b._referencePositionIndex;
				});

				//cout << "\tQuery: " << queryContigIndex << "\t" << _parent._contigOverlaps[queryContigIndex]._size << endl;

				//for(size_t i=0; i<anchors.size(); i++){
				//	cout << "\t\t" << i << "\t" << anchors[i]._referencePositionIndex << "\t" << anchors[i]._queryPositionIndex << endl;
				//}

				u_int32_t overlapLeft = getMaxOverlapLeft(anchors);
				u_int32_t overlapRight = getMaxOverlapRight(anchors, referenceContigLength);
				//cout << "\t" << overlapLeft << " " << overlapRight << endl;
				//cout << "\t" << overlapLeft << " " << overlapRight << endl;
				contigOverlap.first = max(contigOverlap.first, overlapLeft);
				contigOverlap.second = max(contigOverlap.second, overlapRight);
				//getchar();
				/*
				vector<Chain> chains = chainAnchors(anchors);

				cout << "\t\tNb chains: " << chains.size() << endl;

				for(size_t i=0; i<chains.size(); i++){
					
					cout << "\t\t" << i << "\t" << chains[i]._chainingScore << endl;
					for(size_t j=0; j<chains[i]._anchorInterval.size(); j++){
						cout << "\t\t\t" << j << "\t" << chains[i]._anchorInterval[j] << endl;
					}

					const Chain& chain = chains[i];
					const vector<size_t>& chainInterval = chain._anchorInterval;
					if(chainInterval.size() <= 0) continue;;

					const Anchor& firstAnchor = anchors[chainInterval[0]];
					const Anchor& lastAnchor = anchors[chainInterval[chainInterval.size()-1]];

					cout << firstAnchor._referencePositionIndex << " " << lastAnchor._referencePositionIndex << endl;
					
					if(firstAnchor._referencePositionIndex == 0){
						contigOverlap._overlapLeft = max(contigOverlap._overlapLeft, (u_int32_t)(chainInterval.size()-1+_parent._kminmerSize));
						cout << "Overlap left: " << contigOverlap._overlapLeft << endl;
					}
					//const Anchor& firstAnchor = anchors
					//if()
					//cout << "\t" << chains[i]._anchorInterval.size() << endl;

					//if(chains[i]._anchorInterval.size() == referenceContigLength){
						//#pragma atomic
						//_parent._nbDuplicatedContigs += 1;

						//#pragma omp critical
						//{
						//	_parent._isContigDuplicated.insert(referenceContigIndex);
						//}

						//return;
					//}

				}

				//if(referenceContigIndex == 4 && queryContigIndex == 2) getchar();
				*/
			}

			
			//getchar();
			return contigOverlap;
			
		}

		u_int32_t getMaxOverlapLeft(const vector<Anchor>& anchors){

			//cout << "allo" << endl;
			u_int32_t maxOverlapFinal = 0;

			for(int64_t i=0; i<anchors.size(); i++){

				const Anchor& xi = anchors[i];
				//cout << xi._referencePositionIndex << endl;
				if(xi._referencePositionIndex > 0) break;


				int64_t currentReferenceIndex = xi._referencePositionIndex;
				int64_t currentQueryIndex = xi._queryPositionIndex;
				u_int32_t maxOverlap = 1;

				for(int64_t j=i+1; j<anchors.size(); j++){

					const Anchor& xj = anchors[j];

					u_int32_t refGap = xj._referencePositionIndex - currentReferenceIndex;
					if(refGap > 1) break;

					if(xj._referencePositionIndex == currentReferenceIndex+1 && xj._queryPositionIndex == currentQueryIndex+1){
						maxOverlap += 1;
						currentReferenceIndex = xj._referencePositionIndex;
						currentQueryIndex = xj._queryPositionIndex;
					}
					else if(xj._referencePositionIndex == currentReferenceIndex+1 && xj._queryPositionIndex == currentQueryIndex-1){
						maxOverlap += 1;
						currentReferenceIndex = xj._referencePositionIndex;
						currentQueryIndex = xj._queryPositionIndex;
					}
				}

				maxOverlapFinal = max(maxOverlapFinal, maxOverlap);
				//cout << i << "\t" << maxOverlapFinal << "\t" << currentReferenceIndex << "\t" << currentQueryIndex << endl;
				//getchar();
			}

			//cout << "allo2" << endl;
			return maxOverlapFinal;
		}



		u_int32_t getMaxOverlapRight(const vector<Anchor>& anchors, const u_int32_t& referenceSize){

			//cout << "allo" << endl;
			u_int32_t maxOverlapFinal = 0;

			for(int64_t i=((int64_t)anchors.size())-1; i>= 0; i--){

				const Anchor& xi = anchors[i];
				//cout << "aaa " << xi._referencePositionIndex << endl;
				if(xi._referencePositionIndex != referenceSize-1-_parent._kminmerSize+1) break;


				int64_t currentReferenceIndex = xi._referencePositionIndex;
				int64_t currentQueryIndex = xi._queryPositionIndex;
				u_int32_t maxOverlap = 1;

				for(int64_t j=i-1; j >= 0; j--){

					const Anchor& xj = anchors[j];
					//cout << "\t" << j << " " << xj._referencePositionIndex << " " << xj._queryPositionIndex << endl;

					u_int32_t refGap = currentReferenceIndex - xj._referencePositionIndex;
					//cout << "\t" << refGap << endl;

					if(refGap > 1) break;

					if(xj._referencePositionIndex == currentReferenceIndex-1 && xj._queryPositionIndex == currentQueryIndex+1){
						maxOverlap += 1;
						currentReferenceIndex = xj._referencePositionIndex;
						currentQueryIndex = xj._queryPositionIndex;
						//cout << "\t" << "mioum1 " << currentReferenceIndex << endl;
					}
					else if(xj._referencePositionIndex == currentReferenceIndex-1 && xj._queryPositionIndex == currentQueryIndex-1){
						maxOverlap += 1;
						currentReferenceIndex = xj._referencePositionIndex;
						currentQueryIndex = xj._queryPositionIndex;
						//cout << "\t" << "mioum2 " << currentReferenceIndex << endl;
					}
				}

				maxOverlapFinal = max(maxOverlapFinal, maxOverlap);
				//cout << i << "\t" << maxOverlapFinal << "\t" << currentReferenceIndex << "\t" << currentQueryIndex << endl;
				//getchar();
			}

			//cout << "allo2" << endl;
			return maxOverlapFinal;
		}

		void writeContig(vector<MinimizerType> minimizers, const u_int8_t isCircular){

			removeOverlapsSelf(minimizers);

			#pragma omp critical(writeContig)
			{
				u_int32_t contigSize = minimizers.size();
				_parent._outputFile.write((const char*)&contigSize, sizeof(contigSize));
				_parent._outputFile.write((const char*)&isCircular, sizeof(isCircular));
				_parent._outputFile.write((const char*)&minimizers[0], contigSize*sizeof(MinimizerType));
			}

		}



		void removeOverlapsSelf(vector<MinimizerType>& minimizers){


			if(minimizers.size() == 0) return;

			long n = longestPrefixSuffix(minimizers)-1;

			if(n <= 0) return;
			
			minimizers.resize(minimizers.size()-n);
		}


		long longestPrefixSuffix(const vector<MinimizerType>& minimizers)
		{
			int n = minimizers.size();
		
			int lps[n];
			lps[0] = 0; // lps[0] is always 0
		
			// length of the previous
			// longest prefix suffix
			int len = 0;
		
			// the loop calculates lps[i]
			// for i = 1 to n-1
			int i = 1;
			while (i < n)
			{
				if (minimizers[i] == minimizers[len])
				{
					len++;
					lps[i] = len;
					i++;
				}
				else // (pat[i] != pat[len])
				{
					// This is tricky. Consider
					// the example. AAACAAAA
					// and i = 7. The idea is
					// similar to search step.
					if (len != 0)
					{
						len = lps[len-1];
		
						// Also, note that we do
						// not increment i here
					}
					else // if (len == 0)
					{
						lps[i] = 0;
						i++;
					}
				}
			}
		
			int res = lps[n-1];
		
			// Since we are looking for
			// non overlapping parts.
			return (res > n/2)? res/2 : res;
		}

		/*
		vector<Chain> chainAnchors(const vector<Anchor>& anchors) {

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
				argmaxPosition(anchors, scores, i, w, bestScore, bestPrevIndex);

				
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
			
			for(size_t i=0; i<anchors.size(); i++) {
				cout << "\t\t" << i << "\t" << scores[i] << endl;
			}


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

			//if(_print_debug){
			//	cout << "\tParents:" << endl;
			//	for(size_t i=0; i<anchors.size(); i++){
			//		cout << "\t\t" << i << "\t" << anchors[i]._isReversed << "\t" << parents[i] << "\t" << scores[i] << endl;
			//	}
			//}

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
			//if(_print_debug){
			//	cout << "\tBest interval:" << endl;
			//	for(size_t i=0; i<interval.size(); i++){
			//		cout << "\t\t" << interval[i] << endl;
			//	}
			//}

			std::reverse(interval.begin(), interval.end());
			chainIntervals.push_back({maxScore, interval});
			//std::sort(F.begin(), F.end(), [](const Point2& a, const Point2& b){
			//	if(a._score == b._score){
			//		return a._fromIndex > b._fromIndex;
			//	}
			//	return a._score > b._score;
			//});



			return chainIntervals; //F[0]._score;
		}


		//size_t _maxChainingBand;

		void argmaxPosition(const vector<Anchor>& anchors, vector<float>& scores, size_t i, float w, float& bestScore, size_t& bestPrevIndex) {
			
			//cout << "---" << i << endl;
			bestScore = 0;
			bestPrevIndex = i;

			//cout << "\tglurp" << endl;
			//size_t bestJ = 0;
			//float v = 0;
			const Anchor& xi = anchors[i];
			//console.log("\t", xi);
			//for (int64_t j = i-1; j >= 0; j--) {
			for (int64_t j =0; j < i; j++) {

				//cout << "\t\t" << j << endl;

				//cout << "\t\t" << i << " " << k << " " << i-k << " " << maxBand << endl;
				//if(i-j > maxChainingBand) break;
				const Anchor& xj = anchors[j];


				//if(xi._queryPosition - xj._queryPosition > 2500) continue;
				//if(xi._isReversed != xj._isReversed) continue;
				if(xi._referencePositionIndex == xj._referencePositionIndex || xi._queryPositionIndex == xj._queryPositionIndex) continue;
				//cout << "\t" << i << " " << k << " " << xi._isReversed << " " << xj._isReversed << " " << xi._referencePosition << "-" << xi._queryPosition << " " << xj._referencePosition << "-" << xj._queryPosition << endl;



				
				int32_t d_q;// = abs(xi._queryPosition - xj._queryPosition);
				//if(xi._isReversed){
				//	d_q = xj._queryPosition - xi._queryPosition;
				//}
				//else{
					d_q = xi._queryPositionIndex - xj._queryPositionIndex;
				//}
				int32_t d_r = xi._referencePositionIndex - xj._referencePositionIndex;
				
				if((d_q == 1 && d_r == 1) || (d_q == -1 && d_r == 1)){
				}
				else{
					continue;
				}
				
				//cout << xi._referencePositionIndex << "-" << xi._queryPositionIndex << "    " << xj._referencePositionIndex << "-" << xj._queryPositionIndex << "    " <<  d_r << " " << d_q << endl;
				
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

				//if(d_r <= 0) {
				//	continue;
				//}

				//int32_t gap = abs(d_r - d_q);
				//cout << "\t" << xi._referencePosition << "\t" << xj._referencePosition << "\t" << gap << endl;

				//if(gap > 100) {
				//	continue;
				//}
				
				//if(xi._isReversed){
				//	if(xi._queryPosition > xj._queryPosition) continue;
				//}
				//else{
				//if(xi._queryPositionIndex < xj._queryPositionIndex) continue;
				//}
				
				//cout << w << " " << min(smallerDistance(xi,xj,xi._isReversed),w) << " " << gamma(xi,xj,w,a,b) << endl;
				//float anchor_score = min(smallerDistance(xi,xj,xi._isReversed),w) - gamma(xi,xj,w,a,b,xi._isReversed); //w - gap;
				float anchor_score = 1;//w - gap;

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
		*/
	};

	/*
	void writeContigs(){

		_outputFile.open(_outputFilenameContig);


		KminmerParserParallel parser(_inputFilenameContig, 0, 0, false, false, 1);
		parser.parseSequences(WriteContigsFunctor(*this));

		_outputFile.close();

	}



	class WriteContigsFunctor {

		public:

		OverlapRemover2& _parent;

		WriteContigsFunctor(OverlapRemover2& parent) : _parent(parent){
		}

		WriteContigsFunctor(const WriteContigsFunctor& copy) : _parent(copy._parent){
		}

		~WriteContigsFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {


			u_int32_t contigIndex = kminmerList._readIndex;
			//cout << contigIndex << endl;
			//if(contigIndex % 10000 == 0) cout << "Index contigs: " << contigIndex << endl;

			const ContigData& contigOverlap = _parent._contigOverlaps[contigIndex];

			u_int32_t overlapLeft = 0;
			if(contigOverlap._overlapLeft > 0){
				overlapLeft = contigOverlap._overlapLeft + _parent._kminmerSize - 1;
			}

			u_int32_t overlapRight = 0;
			if(contigOverlap._overlapRight > 0){
				overlapRight = contigOverlap._overlapRight + _parent._kminmerSize - 1;
			}

			int64_t indexEnd = ((int64_t)kminmerList._readMinimizers.size()) - ((int64_t)overlapRight);
			//int64_t size = ((int64_t)kminmerList._readMinimizers.size()) - ((int64_t)overlapRight) - ((int64_t)overlapLeft);

			//cout << contigIndex << "\t" << kminmerList._readMinimizers.size() << "\t" << overlapLeft << "\t" << overlapRight << "\t" << indexEnd << endl;
			

			if(overlapLeft >= indexEnd) return;

			//vector<MinimizerType> minimizers = kminmerList._readMinimizers;
			
			vector<MinimizerType> minimizers(kminmerList._readMinimizers.begin() + overlapLeft, kminmerList._readMinimizers.begin() + indexEnd);

			if(minimizers.size() <= _parent._kminmerSize) return;

			u_int32_t contigSize = minimizers.size();
			_parent._outputFile.write((const char*)&contigSize, sizeof(contigSize));
			_parent._outputFile.write((const char*)&kminmerList._isCircular, sizeof(kminmerList._isCircular));
			_parent._outputFile.write((const char*)&minimizers[0], contigSize*sizeof(MinimizerType));

			//getchar();
		}
	};
	*/

};	


#endif 


