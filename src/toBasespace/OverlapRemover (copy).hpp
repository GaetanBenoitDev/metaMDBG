

#ifndef MDBG_METAG_OverlapRemover
#define MDBG_METAG_OverlapRemover

#include "../Commons.hpp"


class OverlapRemover {
  

public:

	string _inputDir;
	string _inputFilenameContig;
	size_t _kminmerSize;

	OverlapRemover(const string& inputDir, const string& inputFilenameContig, size_t kminmerSize){
		_inputDir = inputDir;
		_inputFilenameContig = inputFilenameContig;
		_kminmerSize = kminmerSize;
	}

	/*
	struct ContigOverlap{
		u_int64_t _contigIndex;
		vector<u_int64_t> _minimizers;
		vector<u_int32_t> _nodepath;
		vector<u_int32_t> _nodepath_sorted;
	};


	static bool ContigOverlapComparator_ByLength(const ContigOverlap &a, const ContigOverlap &b){

		if(a._nodepath.size() == b._nodepath.size()){
			for(size_t i=0; i<a._nodepath.size() && i<b._nodepath.size(); i++){
				if(BiGraph::nodeIndex_to_nodeName(a._nodepath[i]) == BiGraph::nodeIndex_to_nodeName(b._nodepath[i])){
					continue;
				}
				else{
					return BiGraph::nodeIndex_to_nodeName(a._nodepath[i]) > BiGraph::nodeIndex_to_nodeName(b._nodepath[i]);
				}
			}
		}


		return a._nodepath.size() > b._nodepath.size();
	}


	unordered_map<KmerVec, u_int32_t> _edgesToIndex;
	u_int32_t _edgeIndex;
	*/

	struct Contig{
		u_int32_t _contigIndex;
		vector<u_int64_t> _minimizers;
	};

	static bool ContigComparator_ByLength(const Contig &a, const Contig &b){

		if(a._minimizers.size() == b._minimizers.size()){
			for(size_t i=0; i<a._minimizers.size() && i<b._minimizers.size(); i++){
				if(a._minimizers[i] == b._minimizers[i]){
					continue;
				}
				else{
					return a._minimizers[i] > b._minimizers[i];
				}
			}
		}


		return a._minimizers.size() < b._minimizers.size();
	}




	vector<Contig> _contigs;
	//vector<bool> isContigRemoved;

	struct MinimizerIndex{
		u_int32_t _contigIndex;
		u_int32_t _pos;
	};

	unordered_map<u_int64_t, vector<MinimizerIndex>> _contigMinimizerIndex;
	unordered_map<u_int64_t, vector<MinimizerIndex>> _contigMinimizerIndex_rev;

	static bool MinimizerIndexComparator(const MinimizerIndex &a, const MinimizerIndex &b){
		return a._contigIndex < b._contigIndex;
	}

	void execute(){


		indexContigMinimizers();
		bool isModification = false;

		while(true){
			isModification = removeOverlaps();

			/*
			for(size_t i=0; i<_contigs.size(); i++){
				if(_contigs[i]._minimizers.size() == 0) continue;
				for(u_int64_t m : _contigs[i]._minimizers){
					cout << m << " ";
				}
				cout << endl;
				//cout << _contigs[i]._contigIndex << " " << _contigs[i]._minimizers.size() << endl;
			}
			*/

			//exit(1);




			//getchar();


			if(!isModification) break;

			vector<Contig> contigsTmp = _contigs;
			_contigs.clear();
			_contigMinimizerIndex.clear();

			u_int64_t contigIndex = 0;

			for(size_t i=0; i<contigsTmp.size(); i++){
				if(contigsTmp[i]._minimizers.size() == 0) continue;

				indexContigMinimizers_read(contigsTmp[i]._minimizers, contigIndex);
				contigIndex += 1;
			}

			//cout << "Nb contigs: " << contigIndex << endl;
			//getchar();
		}

		ofstream outputFile(_inputFilenameContig + ".nooverlaps");
		u_int64_t nbContigs = 0;

		for(size_t i=0; i<_contigs.size(); i++){
			if(_contigs[i]._minimizers.size() == 0) continue;
			
			u_int32_t contigSize = _contigs[i]._minimizers.size();
			outputFile.write((const char*)&contigSize, sizeof(contigSize));
			outputFile.write((const char*)&_contigs[i]._minimizers[0], contigSize*sizeof(u_int64_t));

			cout << contigSize << endl;
			nbContigs += 1;
		}
		outputFile.close();

		fs::remove(_inputFilenameContig);
		fs::rename(_inputFilenameContig + ".nooverlaps", _inputFilenameContig);

		cout << nbContigs << endl;
		getchar();


		//_inputFilenameContig = _inputFilenameContig + ".nooverlaps";

		/*
		_edgeIndex = 0;
		indexEdges();
		indexContigs();
		detectOverlaps();

		ofstream outputFile(_inputFilenameContig + ".nooverlaps");
		
		for(size_t i=0; i<_overContigs.size(); i++){
			if(_overContigs[i]._nodepath.size() == 0) continue;
			
			u_int32_t contigSize = _overContigs[i]._minimizers.size();
			outputFile.write((const char*)&contigSize, sizeof(contigSize));
			outputFile.write((const char*)&_overContigs[i]._minimizers[0], contigSize*sizeof(u_int64_t));
		}
		outputFile.close();

		_edgesToIndex.clear();
		_overContigs.clear();

		_inputFilenameContig = _inputFilenameContig + ".nooverlaps";
		*/
	}


	void indexContigMinimizers(){

		
		indexContigMinimizers_read({1, 2, 3, 4, 5, 6, 7, 8, 9}, 0);
		indexContigMinimizers_read({9, 8, 7, 6, 11, 12, 13, 14}, 1);
		//indexContigMinimizers_read({6, 7, 8, 9, 10, 11, 12, 13, 14}, 1);
		//indexContigMinimizers_read({0, 1, 2, 3, 4}, 2);
		//indexContigMinimizers_read({2, 3, 4, 5, 6}, 3);

		return;

		

		cout << "Indexing contigs" << endl;

		KminmerParser parser(_inputFilenameContig, -1, -1, false, false);
		auto fp = std::bind(&OverlapRemover::indexContigMinimizers_read, this, std::placeholders::_1, std::placeholders::_2);
		parser.parseSequences(fp);
	}

	void indexContigMinimizers_read(const vector<u_int64_t>& readMinimizers, u_int32_t readIndex){

		for(u_int32_t i=0; i<readMinimizers.size(); i++){
			u_int64_t m = readMinimizers[i];
			_contigMinimizerIndex[m].push_back({readIndex, i});

			u_int32_t i_rev = readMinimizers.size()-1-i;
			u_int64_t m_rev = readMinimizers[i_rev];
			_contigMinimizerIndex_rev[m_rev].push_back({readIndex, i_rev});
		}

		/*
		vector<u_int64_t> readMinimizersRC = readMinimizers;
		std::reverse(readMinimizersRC.begin(), readMinimizersRC.end());

		for(u_int32_t i=0; i<readMinimizersRC.size(); i++){
			u_int64_t m = readMinimizersRC[i];
			_contigMinimizerIndex[m].push_back({readIndex, i});
		}
		*/

		_contigs.push_back({readIndex, readMinimizers});
	}


	bool removeOverlaps(){

		//isContigRemoved.resize(_contigs.size(), false);
		cout << "detecting overlaps" << endl;

		//while(true){

			
			//cout << "loop" << endl;

		std::sort(_contigs.begin(), _contigs.end(), ContigComparator_ByLength);
		bool isModification = false;

		/*
		for(long i=_contigs.size()-1; i>=0; i--){

			Contig& contig = _contigs[i];

			if(contig._minimizers.size() == 0) continue;

			unordered_map<u_int32_t, u_int32_t> nbSharedElements;

			for(u_int64_t m : contig._minimizers){
				if(_contigMinimizerIndex.find(m) == _contigMinimizerIndex.end()) continue;

				for(const MinimizerIndex& mIndex : _contigMinimizerIndex[m]){
					if(isContigRemoved[mIndex._contigIndex]) continue;
					if(mIndex._contigIndex == contig._contigIndex) continue;
					nbSharedElements[mIndex._contigIndex] += 1;
				}
			}

			vector<u_int64_t> ms;
			for(const auto& it: nbSharedElements){
				u_int32_t contigIndex = it.first;
				double count = it.second;

				double sharedRate = count / contig._minimizers.size();
				//cout << contig._contigIndex << " " << contig._minimizers.size() << " " << count << " "  << sharedRate << endl;
				//getchar();

				if(sharedRate > 0.75){
					isModification = true;
					ms = contig._minimizers;
					contig._minimizers.clear();
					isContigRemoved[contig._contigIndex] = true;
					break;
				}
			}

			if(nbSharedElements.size() > 0){
				cout << "------" << endl;
				for(u_int64_t m : ms){
					if(_contigMinimizerIndex.find(m) == _contigMinimizerIndex.end()) continue;

					cout << "|" << m << endl;
					for(const MinimizerIndex& mIndex : _contigMinimizerIndex[m]){
						//if(isContigRemoved[mIndex._contigIndex]) continue;
						//if(mIndex._contigIndex == contig._contigIndex) continue;
						//nbSharedElements[mIndex._contigIndex] += 1;
						cout << mIndex._contigIndex << " " << mIndex << endl;
					}
				}
			}
		}
		*/

		
		for(long i=0; i<_contigs.size(); i++){

			Contig& contig = _contigs[i];
			//cout << "---------------------" << endl;
			//cout << i << " " << contig._minimizers.size() << endl;

			if(contig._minimizers.size() == 0) continue;

			//size_t contigSize = contig._minimizers.size();
			u_int64_t overlapSizeLeft = computeOverlapSize_left(contig);
			//cout << "done -----" << endl;

			u_int64_t overlapSizeRight = computeOverlapSize_right(contig);
			cout << overlapSizeLeft << " " << overlapSizeRight << endl;

			if(overlapSizeLeft + overlapSizeRight >= contig._minimizers.size()){
				isModification = true;
				for(size_t i=0; i<contig._minimizers.size(); i++){
					u_int64_t m = contig._minimizers[i];
					for(MinimizerIndex& mIndex : _contigMinimizerIndex[m]){
						if(mIndex._contigIndex == contig._contigIndex && mIndex._pos == i){
							//cout << "Removed: " << mIndex._contigIndex << " " << mIndex._pos << endl;
							mIndex._contigIndex = -1;
							mIndex._pos = -1;
							break;
						}
					}
				}
				contig._minimizers.clear();
			}
			else{
				if(overlapSizeLeft >= _kminmerSize-1){
					isModification = true;
					for(size_t i=0; i<overlapSizeLeft; i++){
						u_int64_t m = contig._minimizers[i];
						for(MinimizerIndex& mIndex : _contigMinimizerIndex[m]){
							if(mIndex._contigIndex == contig._contigIndex && mIndex._pos == i){
								//cout << "Removed: " << mIndex._contigIndex << " " << mIndex._pos << endl;
								mIndex._contigIndex = -1;
								mIndex._pos = -1;
								break;
							}
						}
					}
				}


				if(overlapSizeRight >= _kminmerSize-1){
					isModification = true;
					for(size_t i=0; i<overlapSizeRight; i++){
						size_t ii = contig._minimizers.size()-1-i;
						u_int64_t m = contig._minimizers[ii];
						//cout << contig._minimizers.size()-1-i << " " << m << endl;
						for(MinimizerIndex& mIndex : _contigMinimizerIndex[m]){
							if(mIndex._contigIndex == contig._contigIndex && mIndex._pos == ii){
								//cout << "Removed: " << mIndex._contigIndex << " " << mIndex._pos << endl;
								mIndex._contigIndex = -1;
								mIndex._pos = -1;
								break;
							}
						}
					}

					//contig._minimizers.erase(contig._minimizers.begin()+contigSize-overlapSizeRight, contig._minimizers.begin()+contigSize);
					//for(u_int64_t m : contig._minimizers){
					//	cout << m << " ";
					//}
					//cout << endl;
				}

				
				if(overlapSizeLeft >= _kminmerSize-1){
					contig._minimizers.erase(contig._minimizers.begin(), contig._minimizers.begin() + overlapSizeLeft);
				}
				if(overlapSizeRight >= _kminmerSize-1){
					contig._minimizers.resize(contig._minimizers.size()-overlapSizeRight);
					//contig._minimizers.erase(contig._minimizers.begin()+contig._minimizers.size()-1-overlapSizeRight, contig._minimizers.begin()+contig._minimizers.size()-1);
				}
				//cout << "lala: " << contig._minimizers.size() << endl;
				//for(u_int64_t m : contig._minimizers){
				//	cout << m << " ";
				//}
				//cout << endl;
			}



			if(contig._minimizers.size() <= _kminmerSize){
				isModification = true;
				contig._minimizers.clear();
			}


			/*
			if(overlapSize >= _kminmerSize-1){

				contig._minimizers.erase(contig._minimizers.begin(), contig._minimizers.begin() + overlapSize);

				for(size_t i=0; i<overlapSize; i++){
					u_int64_t m = contig._minimizers[i];
					for(MinimizerIndex& mIndex : _contigMinimizerIndex[m]){
						if(mIndex._contigIndex == contig._contigIndex && mIndex._pos == i){
							mIndex._contigIndex = -1;
							mIndex._pos = -1;
						}
					}
				}
			}
			*/

			/*
			//cout << computeOverlapSize(contig, false) << endl;
			//cout << computeOverlapSize(contig, true) << endl;

			std::reverse(contig._minimizers.begin(), contig._minimizers.end());
			
			overlapSize = computeOverlapSize(contig);
			cout << overlapSize << endl;
			if(overlapSize >= _kminmerSize-1){
				isModification = true;
				contig._minimizers.erase(contig._minimizers.begin(), contig._minimizers.begin() + overlapSize);
			}


			std::reverse(contig._minimizers.begin(), contig._minimizers.end());


			if(contig._minimizers.size() <= _kminmerSize){
				contig._minimizers.clear();
				isModification = true;
			}
			*/
			

			//}
			

			/*
			for(size_t i=0; i<_contigs.size(); i++){
				
				if(_contigs[i]._nodepath.size() == 0) continue;
				//if(_invalidContigIndex.find(_contigs[i]._readIndex) != _invalidContigIndex.end()) continue;

				for(long j=_contigs.size()-1; j>=i+1; j--){

					if(_contigs[j]._nodepath.size() == 0) continue;

					
					double nbShared = Utils::computeSharedElements(_overContigs[i]._nodepath_sorted, _overContigs[j]._nodepath_sorted);
					double sharedRate_1 = nbShared / _overContigs[i]._nodepath_sorted.size();
					double sharedRate_2 = nbShared / _overContigs[j]._nodepath_sorted.size();


					if(sharedRate_2 > 0.75){

						//cout << _overContigs[i]._nodepath_sorted.size() << " " << _overContigs[j]._nodepath_sorted.size() << " " << sharedRate_2 << endl;
						//_invalidContigIndex.insert(_contigs[j]._readIndex);
						//break;
						isModification = true;
						_overContigs[j]._minimizers.clear();
						_overContigs[j]._nodepath.clear();
						_overContigs[j]._nodepath_sorted.clear();

					}

				}
			}
			*/

			/*
			for(size_t i=0; i<_contigs.size(); i++){
				
				if(_overContigs[i]._nodepath.size() == 0) continue;
				//if(_invalidContigIndex.find(_contigs[i]._readIndex) != _invalidContigIndex.end()) continue;

				for(long j=_overContigs.size()-1; j>=i+1; j--){

					if(_overContigs[j]._nodepath.size() == 0) continue;
					
					
					//if(_invalidContigIndex.find(_contigs[j]._readIndex) != _invalidContigIndex.end()) continue;

					unordered_set<u_int32_t> sharedElements;
					Utils::collectSharedElements(_overContigs[i]._nodepath_sorted, _overContigs[j]._nodepath_sorted, sharedElements);

					if(sharedElements.size() == 0) continue;

					if(removeOverlap(_overContigs[i]._nodepath, _overContigs[j]._nodepath, sharedElements, true, _overContigs[j])){
						isModification = true;
					}

					if(removeOverlap(_overContigs[i]._nodepath, _overContigs[j]._nodepath, sharedElements, false, _overContigs[j])){
						isModification = true;
					}

					//if(_overContigs[j]._contigIndex == 1097) getchar();

					//getchar();
				}
			}

			if(!isModification) break;
			*/
		}

		return isModification;
	}

	u_int64_t computeOverlapSize_left(Contig& contig){



		//unordered_set<u_int32_t> validContigIndex;
		vector<MinimizerIndex> currentContigIndex;

		u_int64_t overlapSize = 0;

		for(size_t p=0; p<contig._minimizers.size(); p++){
			
			u_int64_t m = contig._minimizers[p];
			vector<MinimizerIndex> nextContigIndex;

			for(const MinimizerIndex& mIndex : _contigMinimizerIndex[m]){
				if(mIndex._contigIndex == contig._contigIndex) continue;

				if(p == 0){
					//validContigIndex.insert(mIndex._contigIndex);
					currentContigIndex.push_back(mIndex);
				}
				else{
					nextContigIndex.push_back(mIndex);
				}
			}

			//cout << "pos: " << p << endl;

			
			if(p > 0){
				vector<MinimizerIndex> sharedContigIndexValid;

				std::sort(currentContigIndex.begin(), currentContigIndex.end(), MinimizerIndexComparator);
				std::sort(nextContigIndex.begin(), nextContigIndex.end(), MinimizerIndexComparator);

				//for(const MinimizerIndex& mIndex : currentContigIndex){
				//	cout << mIndex._contigIndex << " " << mIndex._pos << endl;
				//}
				//for(const MinimizerIndex& mIndex : nextContigIndex){
				//	cout << mIndex._contigIndex << " " << mIndex._pos << endl;
				//}

				size_t i=0;
				size_t j=0;
				while(i < currentContigIndex.size() && j < nextContigIndex.size()){

					//cout << p << " " << currentContigIndex[i]._contigIndex << " " << nextContigIndex[j]._contigIndex << endl;

					if(currentContigIndex[i]._contigIndex == nextContigIndex[j]._contigIndex){

						if(nextContigIndex[j]._pos == currentContigIndex[i]._pos+1){
							sharedContigIndexValid.push_back(nextContigIndex[j]);
						}

						i += 1;
						j += 1;
					}
					else if(currentContigIndex[i]._contigIndex < nextContigIndex[j]._contigIndex){
						i += 1;
					}
					else{
						j += 1;
					}

				}

				currentContigIndex = sharedContigIndexValid;
			}

			if(currentContigIndex.size() == 0) break;
			overlapSize += 1;
		}

		return overlapSize;
	}

	u_int64_t computeOverlapSize_right(Contig& contig){

		long firstPos = contig._minimizers.size()-1;
		vector<MinimizerIndex> currentContigIndex;

		u_int64_t overlapSize = 0;

		for(long p=contig._minimizers.size()-1; p>=0; p--){
			
			u_int64_t m = contig._minimizers[p];
			vector<MinimizerIndex> nextContigIndex;

			for(const MinimizerIndex& mIndex : _contigMinimizerIndex[m]){
				if(mIndex._contigIndex == contig._contigIndex) continue;

				if(p == firstPos){
					//validContigIndex.insert(mIndex._contigIndex);
					currentContigIndex.push_back(mIndex);
				}
				else{
					nextContigIndex.push_back(mIndex);
				}
			}

			//cout << "pos: " << p << endl;

			
			if(p < firstPos){
				vector<MinimizerIndex> sharedContigIndexValid;

				std::sort(currentContigIndex.begin(), currentContigIndex.end(), MinimizerIndexComparator);
				std::sort(nextContigIndex.begin(), nextContigIndex.end(), MinimizerIndexComparator);

				//for(const MinimizerIndex& mIndex : currentContigIndex){
				//	cout << mIndex._contigIndex << " " << mIndex._pos << endl;
				//}
				//for(const MinimizerIndex& mIndex : nextContigIndex){
				//	cout << mIndex._contigIndex << " " << mIndex._pos << endl;
				//}

				size_t i=0;
				size_t j=0;
				while(i < currentContigIndex.size() && j < nextContigIndex.size()){

					//cout << p << " " << currentContigIndex[i]._contigIndex << " " << nextContigIndex[j]._contigIndex << endl;

					if(currentContigIndex[i]._contigIndex == nextContigIndex[j]._contigIndex){

						if(nextContigIndex[j]._pos == currentContigIndex[i]._pos-1){
							sharedContigIndexValid.push_back(nextContigIndex[j]);
						}

						i += 1;
						j += 1;
					}
					else if(currentContigIndex[i]._contigIndex < nextContigIndex[j]._contigIndex){
						i += 1;
					}
					else{
						j += 1;
					}

				}

				currentContigIndex = sharedContigIndexValid;
			}

			if(currentContigIndex.size() == 0) break;
			overlapSize += 1;
		}

		return overlapSize;
	}

	/*
	void indexContigs(){

		cout << "Loading min contigs" << endl;

		KminmerParser parser(_inputFilenameContig, _minimizerSize, _kminmerSize-1, false, false);
		auto fp = std::bind(&ToBasespaceNoCorrection::indexContigs_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
		parser.parseMinspace(fp);
	}

	void indexContigs_read(const vector<u_int64_t>& readMinimizers, const vector<ReadKminmerComplete>& kminmersInfos, u_int64_t readIndex){

		vector<u_int32_t> nodepath;

		for(size_t i=0; i<kminmersInfos.size(); i++){
			
			const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];
			KmerVec vec = kminmerInfo._vec;

			nodepath.push_back(_edgesToIndex[vec]);
		}


		vector<u_int32_t> nodepath_sorted = nodepath;
		std::sort(nodepath_sorted.begin(), nodepath_sorted.end());
		_overContigs.push_back({readIndex, readMinimizers, nodepath, nodepath_sorted});

	}
	*/
	/*
	void detectOverlaps(){

		cout << "detecting overlaps" << endl;

		while(true){

			
			cout << "loop" << endl;

			std::sort(_overContigs.begin(), _overContigs.end(), ContigOverlapComparator_ByLength);
			bool isModification = false;


			for(size_t i=0; i<_overContigs.size(); i++){
				
				if(_overContigs[i]._nodepath.size() == 0) continue;
				//if(_invalidContigIndex.find(_contigs[i]._readIndex) != _invalidContigIndex.end()) continue;

				for(long j=_overContigs.size()-1; j>=i+1; j--){

					if(_overContigs[j]._nodepath.size() == 0) continue;

					double nbShared = Utils::computeSharedElements(_overContigs[i]._nodepath_sorted, _overContigs[j]._nodepath_sorted);
					double sharedRate_1 = nbShared / _overContigs[i]._nodepath_sorted.size();
					double sharedRate_2 = nbShared / _overContigs[j]._nodepath_sorted.size();


					if(sharedRate_2 > 0.75){

						//cout << _overContigs[i]._nodepath_sorted.size() << " " << _overContigs[j]._nodepath_sorted.size() << " " << sharedRate_2 << endl;
						//_invalidContigIndex.insert(_contigs[j]._readIndex);
						//break;
						isModification = true;
						_overContigs[j]._minimizers.clear();
						_overContigs[j]._nodepath.clear();
						_overContigs[j]._nodepath_sorted.clear();

					}

				}
			}

			for(size_t i=0; i<_overContigs.size(); i++){
				
				if(_overContigs[i]._nodepath.size() == 0) continue;
				//if(_invalidContigIndex.find(_contigs[i]._readIndex) != _invalidContigIndex.end()) continue;

				for(long j=_overContigs.size()-1; j>=i+1; j--){

					if(_overContigs[j]._nodepath.size() == 0) continue;
					
					
					//if(_invalidContigIndex.find(_contigs[j]._readIndex) != _invalidContigIndex.end()) continue;

					unordered_set<u_int32_t> sharedElements;
					Utils::collectSharedElements(_overContigs[i]._nodepath_sorted, _overContigs[j]._nodepath_sorted, sharedElements);

					if(sharedElements.size() == 0) continue;

					if(removeOverlap(_overContigs[i]._nodepath, _overContigs[j]._nodepath, sharedElements, true, _overContigs[j])){
						isModification = true;
					}

					if(removeOverlap(_overContigs[i]._nodepath, _overContigs[j]._nodepath, sharedElements, false, _overContigs[j])){
						isModification = true;
					}

					//if(_overContigs[j]._contigIndex == 1097) getchar();

					//getchar();
				}
			}

			if(!isModification) break;
		}

	}

	bool removeOverlap(const vector<u_int32_t>& nodePath, vector<u_int32_t>& nodePath_shorter, unordered_set<u_int32_t>& sharedElements, bool right, ContigOverlap& contig){

		bool isModification = false;

		if(right){
			std::reverse(nodePath_shorter.begin(), nodePath_shorter.end());
		}
		
		u_int32_t overlapSize = computeOverlapSize(nodePath, nodePath_shorter, sharedElements);
		//cout << overlapSize << " " << contig._minimizers.size() << " " << contig._nodepath.size() << " " << contig._nodepath_sorted.size() << endl;
		//if(overlapSize > 0) overlapSize += 1;
		//if(right && overlapSize > 0) 

		if(overlapSize > 0){
			isModification = true;
			overlapSize += (_kminmerSize-1-1);
			if(overlapSize >= nodePath_shorter.size()){
				contig._minimizers.clear();
				contig._nodepath.clear();
				contig._nodepath_sorted.clear();
			}
			else{
				nodePath_shorter.erase(nodePath_shorter.begin(), nodePath_shorter.begin() + overlapSize);

				vector<u_int64_t> minimizers = contig._minimizers;

				if(right){
					std::reverse(minimizers.begin(), minimizers.end());
				}

				minimizers.erase(minimizers.begin(), minimizers.begin() + overlapSize);

				if(right){
					std::reverse(minimizers.begin(), minimizers.end());
					std::reverse(nodePath_shorter.begin(), nodePath_shorter.end());
				}

				contig._minimizers = minimizers;
				contig._nodepath = nodePath_shorter;
				contig._nodepath_sorted.clear();
				for(u_int32_t index : contig._nodepath){
					contig._nodepath_sorted.push_back(index);
				}
				std::sort(contig._nodepath_sorted.begin(), contig._nodepath_sorted.end());


			}

		}
		else{

			if(right){
				std::reverse(nodePath_shorter.begin(), nodePath_shorter.end());
			}

		}

		if(contig._minimizers.size() <= _kminmerSize){
			contig._minimizers.clear();
			contig._nodepath.clear();
			contig._nodepath_sorted.clear();
			isModification = true;
		}


		return isModification;

	}

	u_int32_t computeOverlapSize(const vector<u_int32_t>& nodePath, const vector<u_int32_t>& nodePath_shorter, unordered_set<u_int32_t>& sharedElements){

		if(sharedElements.size() == 0) return 0;

		size_t uniqueRunLength = 0;
		bool isUnique = false;

		u_int32_t nonUniquePos = 0;

		for(size_t i=0; i<nodePath_shorter.size(); i++){

			if(sharedElements.find(nodePath_shorter[i]) == sharedElements.end()){
				if(isUnique) uniqueRunLength += 1;
				isUnique = true;
				break;
			}
			else{
				nonUniquePos = (i+1);
				isUnique = false;
				uniqueRunLength = 0;
			}

			if(uniqueRunLength > 0) break;

		}

		//if(nonUniquePos > 0) nonUniquePos += 1;
		return nonUniquePos;
	}






	class OverlapRemover(){

	};




	gzFile _queryContigFile;
	unordered_set<u_int32_t> _duplicatedContigIndex;

	void removeDuplicatePost(){

		_nbContigsPost = 0;

		string outputMappingFilename = _filename_outputContigs + ".map";
		const string& filenameQuery = _filename_outputContigs + ".query";
		_queryContigFile = gzopen(filenameQuery.c_str(),"wb");
		

		auto fp = std::bind(&ToBasespaceNoCorrection::dumpSmallContigs_read, this, std::placeholders::_1);
		ReadParser readParser(_filename_outputContigs, true, false);
		readParser.parse(fp);

		gzclose(_queryContigFile);

		string command = "minimap2 -x map-hifi " + _filename_outputContigs + " " + filenameQuery + " > " + outputMappingFilename;
		Utils::executeCommand(command);



		ifstream mappingFile(outputMappingFilename);
        vector<string>* fields = new vector<string>();
        vector<string>* fields_optional = new vector<string>();

		string line;
		while (getline(mappingFile, line)) {


            GfaParser::tokenize(line, fields, '\t');

			const string& readName = (*fields)[0];
			const string& contigName = (*fields)[5];
			if(readName == contigName) continue;

			//u_int64_t queryLength = stoull((*fields)[1]);
			u_int64_t targetLength = stoull((*fields)[6]);
			double queryLength = stoull((*fields)[1]);

			if(targetLength < queryLength) continue;

			u_int64_t nbMatches = stoull((*fields)[9]);
			double alignLength = stoull((*fields)[10]);

			if(alignLength / queryLength < 0.95) continue;

			cout << (nbMatches / alignLength) << " " << (nbMatches / queryLength) << endl;
			
			for(size_t i=12; i<fields->size(); i++){

				//cout << (*fields)[i] << endl;

				GfaParser::tokenize((*fields)[i], fields_optional, ':');

				if((*fields_optional)[0] == "dv"){
					float divergence = std::stof((*fields_optional)[2]);

					//cout << (*fields_optional)[2] << endl;
					//cout << contigName << " " << readName << " " << (alignLength/queryLength*100) << " " << (divergence*100) << "     " << queryLength << " " << targetLength << endl;
					if(divergence < 0.05){
						string name = readName;
						size_t pos = name.find("ctg");
						name.erase(pos, 3);
						u_int32_t contigIndex = stoull(name);
						//cout << "Duplicate: " << contigIndex << endl;

						_duplicatedContigIndex.insert(contigIndex);
					}
				}

			}

			//getchar();
		}

		mappingFile.close();

		dumpDereplicatedContigs();

		fs::remove(_inputDir + "/contig_data.txt");
		fs::rename(_inputDir + "/contig_data_derep.txt", _inputDir + "/contig_data.txt");
		_duplicatedContigIndex.clear();
	}

	u_int64_t _nbContigsPost;

	void dumpSmallContigs_read(const Read& read){
		
		//cout << read._seq.size() << endl;
		if(read._seq.size() > _meanReadLength*3) return;

		string header = ">" + read._header + '\n';
		gzwrite(_queryContigFile, (const char*)&header[0], header.size());
		string contigSequence = read._seq + '\n';
		gzwrite(_queryContigFile, (const char*)&contigSequence[0], contigSequence.size());
		
	}

	ofstream _outputContigFileDerep;

	void dumpDereplicatedContigs(){

		string contigFilename = _inputDir + "/contig_data_derep.txt";
		_outputContigFileDerep = ofstream(contigFilename);

		KminmerParser parser(_inputFilenameContig, _minimizerSize, _kminmerSize, false, false);
		auto fp = std::bind(&ToBasespaceNoCorrection::dumpDereplicatedContigs_read, this, std::placeholders::_1, std::placeholders::_2);
		parser.parseSequences(fp);

		_outputContigFileDerep.close();

	}

	void dumpDereplicatedContigs_read(const vector<u_int64_t>& readMinimizers, u_int64_t readIndex){
		
		if(_duplicatedContigIndex.find(readIndex) != _duplicatedContigIndex.end()) return;

		u_int32_t contigSize = readMinimizers.size();
		_outputContigFileDerep.write((const char*)&contigSize, sizeof(contigSize));
		_outputContigFileDerep.write((const char*)&readMinimizers[0], contigSize*sizeof(u_int64_t));

		_nbContigsPost += 1;
		//cout << "Dump: " << readIndex << " " << readMinimizers.size() << endl;
	}
	*/
};	


#endif 


