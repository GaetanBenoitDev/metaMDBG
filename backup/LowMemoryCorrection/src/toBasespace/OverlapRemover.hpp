

#ifndef MDBG_METAG_OverlapRemover
#define MDBG_METAG_OverlapRemover

#include "../Commons.hpp"


class OverlapRemover {
  

public:

	string _inputDir;
	string _inputFilenameContig;
	size_t _kminmerSize;
	ofstream& _logFile;

	OverlapRemover(const string& inputDir, const string& inputFilenameContig, size_t kminmerSize, ofstream& logFile) : _logFile(logFile){
		_inputDir = inputDir;
		_inputFilenameContig = inputFilenameContig;
		_kminmerSize = kminmerSize-1;
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
		vector<MinimizerType> _minimizers;
		vector<u_int32_t> _kminmers;
		u_int8_t isCircular;
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

	struct KminmerIndex{
		u_int32_t _contigIndex;
		//u_int32_t _pos;
	};

	unordered_map<u_int32_t, vector<KminmerIndex>> _kminmerIndex;

	static bool KminmerIndexComparator(const KminmerIndex &a, const KminmerIndex &b){
		return a._contigIndex < b._contigIndex;
	}


	void execute(){

		indexKminmers();
		indexContigs();
		bool isModification = false;

		while(true){
			isModification = removeOverlaps();

			/*
			for(size_t i=0; i<_contigs.size(); i++){
				if(_contigs[i]._minimizers.size() == 0) continue;
				for(u_int64_t m : _contigs[i]._minimizers){
					_logFile << m << " ";
				}
				_logFile << endl;
				//_logFile << _contigs[i]._contigIndex << " " << _contigs[i]._minimizers.size() << endl;
			}
			*/

			//exit(1);




			//getchar();


			if(!isModification) break;

			vector<Contig> contigsTmp = _contigs;
			_contigs.clear();
			_kminmerIndex.clear();

			u_int32_t contigIndex = 0;

			u_int64_t checksum = 0;

			for(size_t i=0; i<contigsTmp.size(); i++){
				if(contigsTmp[i]._minimizers.size() == 0) continue;

				vector<u_int32_t> minimizersPos; 
				vector<u_int64_t> rlePositions; 
				vector<KmerVec> kminmers; 
				vector<ReadKminmer> kminmersInfo;
				MDBG::getKminmers(-1, _kminmerSize, contigsTmp[i]._minimizers, minimizersPos, kminmers, kminmersInfo, rlePositions, 0, false);

				indexContigs_read(contigsTmp[i]._minimizers, kminmers, kminmersInfo, contigsTmp[i].isCircular, contigIndex);
				contigIndex += 1;

				for(u_int64_t m : contigsTmp[i]._minimizers){
					checksum += m*contigsTmp[i]._minimizers.size();
				}
			}

			_logFile << "Nb contigs: " << contigIndex << endl;
			cout << "OverlapRemover checksum: " << checksum << endl;
		
			//getchar();
		}
		
		removeOverlapsSelf();

		ofstream outputFile(_inputFilenameContig + ".nooverlaps");
		u_int64_t nbContigs = 0;

		for(size_t i=0; i<_contigs.size(); i++){
			if(_contigs[i]._minimizers.size() == 0) continue;
			
			u_int32_t contigSize = _contigs[i]._minimizers.size();
			outputFile.write((const char*)&contigSize, sizeof(contigSize));
			outputFile.write((const char*)&_contigs[i].isCircular, sizeof(_contigs[i].isCircular));
			outputFile.write((const char*)&_contigs[i]._minimizers[0], contigSize*sizeof(MinimizerType));

			//_logFile << contigSize << endl;
			nbContigs += 1;
		}
		outputFile.close();

		fs::remove(_inputFilenameContig);
		fs::rename(_inputFilenameContig + ".nooverlaps", _inputFilenameContig);

		_logFile << nbContigs << endl;
		//getchar();


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


	unordered_map<KmerVec, u_int32_t> _kminmerToIndex;
	u_int32_t _kminmerID;

	void indexKminmers(){

		_logFile << "Indexing kminmers" << endl;

		_kminmerID = 0;

		KminmerParser parser(_inputFilenameContig, -1, _kminmerSize, false, false);
		auto fp = std::bind(&OverlapRemover::indexKminmers_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
		parser.parse(fp);
	}

	void indexKminmers_read(const vector<MinimizerType>& readMinimizers, const vector<KmerVec>& vecs, const vector<ReadKminmer>& kminmersInfos, u_int8_t isCircular, u_int32_t readIndex){
		
		for(u_int32_t i=0; i<vecs.size(); i++){
			
			KmerVec vec = vecs[i];

			if(_kminmerToIndex.find(vec) == _kminmerToIndex.end()){
				_kminmerToIndex[vec] = _kminmerID;
				_kminmerID += 1;
			}
		}
	}

	void indexContigs(){

		_logFile << "Indexing contigs" << endl;

		KminmerParser parser(_inputFilenameContig, -1, _kminmerSize, false, false);
		auto fp = std::bind(&OverlapRemover::indexContigs_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
		parser.parse(fp);

		//_kminmerToIndex.clear();
	}

	void indexContigs_read(const vector<MinimizerType>& readMinimizers, const vector<KmerVec>& vecs, const vector<ReadKminmer>& kminmersInfos, u_int8_t isCircular, u_int32_t readIndex){

		unordered_set<u_int32_t> indexedKminmer;
		vector<u_int32_t> nodepath;
		
		for(u_int32_t i=0; i<vecs.size(); i++){
			
			//const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];
			KmerVec vec = vecs[i];
			u_int32_t kminmerID = _kminmerToIndex[vec];

			if(indexedKminmer.find(kminmerID) == indexedKminmer.end()){
				_kminmerIndex[kminmerID].push_back({readIndex});
				indexedKminmer.insert(kminmerID);
			}

			nodepath.push_back(kminmerID);
			
		}
		
		_contigs.push_back({readIndex, readMinimizers, nodepath, isCircular});
	}



	void removeOverlapsSelf(){

		//isContigRemoved.resize(_contigs.size(), false);
		_logFile << "Removing self overlaps: " << _contigs.size()<< endl;

		//while(true){

			
			//_logFile << "loop" << endl;

		//std::sort(_contigs.begin(), _contigs.end(), ContigComparator_ByLength);
		//bool isModification = false;

		
		for(long i=0; i<_contigs.size(); i++){

			Contig& contig = _contigs[i];
			//_logFile << "---------------------" << endl;
			//_logFile << i << " " << contig._minimizers.size() << endl;

			if(contig._minimizers.size() == 0) continue;

			long n = longestPrefixSuffix(contig)-1;

			if(n <= 0) continue;
			

			//_logFile << "lala1: " << n << endl;
			
			contig._kminmers.resize(contig._kminmers.size()-n);
			//_logFile << contig._minimizers.size() << endl;
			contig._minimizers.resize(contig._minimizers.size()-n);

			//_logFile << "lala2: " << longestPrefixSuffix(contig) << endl;
		}
	}

	bool removeOverlaps(){

		//isContigRemoved.resize(_contigs.size(), false);
		_logFile << "detecting overlaps" << endl;

		//while(true){

			
			//_logFile << "loop" << endl;

		std::sort(_contigs.begin(), _contigs.end(), ContigComparator_ByLength);
		bool isModification = false;

		
		for(long i=0; i<_contigs.size(); i++){

			Contig& contig = _contigs[i];
			//_logFile << "---------------------" << endl;
			//_logFile << i << " " << contig._minimizers.size() << endl;

			if(contig._minimizers.size() == 0) continue;



			//_logFile << "-------" << endl;
			//for(u_int32_t nodeName : contig._kminmers){
			//	_logFile << nodeName << " ";
			//}
			//_logFile << endl;
			//_logFile << contig._kminmers.size() << endl;

			//size_t contigSize = contig._minimizers.size();
			u_int64_t overlapSizeLeft = computeOverlapSize_left(contig);
			//_logFile << "done -----" << endl;

			
			u_int64_t overlapSizeRight = computeOverlapSize_right(contig);
			//_logFile << overlapSizeLeft << " " << overlapSizeRight << endl;

			if(contig._minimizers.size() > 1000) continue;
			//u_int64_t overlapTotalMin = 0;
			if(overlapSizeLeft > 0){
				overlapSizeLeft += _kminmerSize -1;
				//overlapTotalMin += overlapSizeLeft + _kminmerSize -1;
			}
			if(overlapSizeRight > 0){
				overlapSizeRight += _kminmerSize -1;
				//overlapTotalMin += overlapSizeRight + _kminmerSize -1;
			}

			if(overlapSizeLeft + overlapSizeRight == 0){
				continue;
			}
			else if(overlapSizeLeft + overlapSizeRight >= contig._kminmers.size()){//} || overlapTotalMin >= contig._minimizers.size()){

				//_logFile << "remove total" << endl;

				isModification = true;
				for(size_t i=0; i<contig._kminmers.size(); i++){
					
					for(KminmerIndex& mIndex : _kminmerIndex[contig._kminmers[i]]){
						if(mIndex._contigIndex == contig._contigIndex){//} && mIndex._pos == i){
							//_logFile << "Removed: " << mIndex._contigIndex << " " << mIndex._pos << endl;
							mIndex._contigIndex = -1;
							//mIndex._pos = -1;
							break;
						}
					}
				}
				contig._minimizers.clear();
				contig._kminmers.clear();
			}
			else{
				
				//_logFile << "remove left and right" << endl;

				isModification = true;
				for(size_t i=0; i<overlapSizeLeft; i++){
					
					for(KminmerIndex& mIndex : _kminmerIndex[contig._kminmers[i]]){
						if(mIndex._contigIndex == contig._contigIndex){//} && mIndex._pos == i){
							//_logFile << "Removed: " << mIndex._contigIndex << " " << mIndex._pos << endl;
							mIndex._contigIndex = -1;
							//mIndex._pos = -1;
							break;
						}
					}
				}


				for(size_t i=0; i<overlapSizeRight; i++){
					size_t ii = contig._kminmers.size()-1-i;
					
					//_logFile << contig._minimizers.size()-1-i << " " << m << endl;
					for(KminmerIndex& mIndex : _kminmerIndex[contig._kminmers[ii]]){
						if(mIndex._contigIndex == contig._contigIndex){//} && mIndex._pos == i){
							//_logFile << "Removed: " << mIndex._contigIndex << " " << mIndex._pos << endl;
							mIndex._contigIndex = -1;
							//mIndex._pos = -1;
							break;
						}
					}
				}

				//contig._minimizers.erase(contig._minimizers.begin()+contigSize-overlapSizeRight, contig._minimizers.begin()+contigSize);
				//for(u_int64_t m : contig._minimizers){
				//	_logFile << m << " ";
				//}
				//_logFile << endl;

				

				if(overlapSizeLeft > 0){
					//overlapSizeLeft += (_kminmerSize-1);
					contig._kminmers.erase(contig._kminmers.begin(), contig._kminmers.begin() + overlapSizeLeft);
					contig._minimizers.erase(contig._minimizers.begin(), contig._minimizers.begin() + overlapSizeLeft);
				}
				if(overlapSizeRight > 0){
					//overlapSizeRight += (_kminmerSize-1);
					contig._kminmers.resize(contig._kminmers.size()-overlapSizeRight);
					//_logFile << contig._minimizers.size() << endl;
					contig._minimizers.resize(contig._minimizers.size()-overlapSizeRight);
					//_logFile << contig._minimizers.size() << endl;
					//contig._minimizers.erase(contig._minimizers.begin()+contig._minimizers.size()-1-overlapSizeRight, contig._minimizers.begin()+contig._minimizers.size()-1);
				}
				//_logFile << "lala: " << contig._minimizers.size() << endl;
				//for(u_int64_t m : contig._minimizers){
				//	_logFile << m << " ";
				//}
				//_logFile << endl;
			}



			if(contig._minimizers.size() <= _kminmerSize+1){
				isModification = true;
				contig._minimizers.clear();
				contig._kminmers.clear();
			}
			
			//for(u_int64_t m : contig._minimizers){
			//	_logFile << m << " ";
			//}
			//_logFile << endl;
			//for(u_int32_t m : contig._kminmers){
			//	_logFile << m << " ";
			//}
			//_logFile << endl;

			//if(contig._kminmers.size() > 0)	getchar();

		}

		return isModification;
	}

	u_int32_t computeOverlapSize_left(const Contig& contig){

		vector<KminmerIndex> currentContigIndex;
		u_int64_t overlapSize = 0;

		for(size_t p=0; p<contig._kminmers.size(); p++){

			vector<KminmerIndex> nextContigIndex;

			for(const KminmerIndex& mIndex : _kminmerIndex[contig._kminmers[p]]){
				if(mIndex._contigIndex == -1) continue;
				if(mIndex._contigIndex == contig._contigIndex) continue;

				if(p == 0){
					//validContigIndex.insert(mIndex._contigIndex);
					currentContigIndex.push_back(mIndex);
				}
				else{
					nextContigIndex.push_back(mIndex);
				}
			}

			//_logFile << "pos: " << p << endl;

			
			if(p > 0){
				vector<KminmerIndex> sharedContigIndexValid;

				std::sort(currentContigIndex.begin(), currentContigIndex.end(), KminmerIndexComparator);
				std::sort(nextContigIndex.begin(), nextContigIndex.end(), KminmerIndexComparator);

				//for(const MinimizerIndex& mIndex : currentContigIndex){
				//	_logFile << mIndex._contigIndex << " " << mIndex._pos << endl;
				//}
				//for(const MinimizerIndex& mIndex : nextContigIndex){
				//	_logFile << mIndex._contigIndex << " " << mIndex._pos << endl;
				//}

				size_t i=0;
				size_t j=0;
				while(i < currentContigIndex.size() && j < nextContigIndex.size()){

					//_logFile << p << " " << currentContigIndex[i]._contigIndex << " " << nextContigIndex[j]._contigIndex << endl;

					if(currentContigIndex[i]._contigIndex == nextContigIndex[j]._contigIndex){

						//if(contig._minimizers.size() > 1000){
						//	_logFile << p << " " << nextContigIndex[j]._contigIndex << endl;
						//}

						sharedContigIndexValid.push_back(nextContigIndex[j]);
						//if(nextContigIndex[j]._pos == currentContigIndex[i]._pos+1){
						//	sharedContigIndexValid.push_back(nextContigIndex[j]);
						//}

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

	u_int32_t computeOverlapSize_right(const Contig& contig){

		long firstPos = contig._kminmers.size()-1;
		vector<KminmerIndex> currentContigIndex;
		u_int64_t overlapSize = 0;

		for(long p=contig._kminmers.size()-1; p>=0; p--){

			vector<KminmerIndex> nextContigIndex;

			for(const KminmerIndex& mIndex : _kminmerIndex[contig._kminmers[p]]){
				if(mIndex._contigIndex == contig._contigIndex) continue;
				if(mIndex._contigIndex == -1) continue;

				if(p == firstPos){
					//validContigIndex.insert(mIndex._contigIndex);
					currentContigIndex.push_back(mIndex);
				}
				else{
					nextContigIndex.push_back(mIndex);
				}
			}

			//_logFile << "pos: " << p << endl;

			
			if(p < firstPos){
				vector<KminmerIndex> sharedContigIndexValid;

				std::sort(currentContigIndex.begin(), currentContigIndex.end(), KminmerIndexComparator);
				std::sort(nextContigIndex.begin(), nextContigIndex.end(), KminmerIndexComparator);

				//for(const MinimizerIndex& mIndex : currentContigIndex){
				//	_logFile << mIndex._contigIndex << " " << mIndex._pos << endl;
				//}
				//for(const MinimizerIndex& mIndex : nextContigIndex){
				//	_logFile << mIndex._contigIndex << " " << mIndex._pos << endl;
				//}

				size_t i=0;
				size_t j=0;
				while(i < currentContigIndex.size() && j < nextContigIndex.size()){

					//_logFile << p << " " << currentContigIndex[i]._contigIndex << " " << nextContigIndex[j]._contigIndex << endl;

					if(currentContigIndex[i]._contigIndex == nextContigIndex[j]._contigIndex){

						sharedContigIndexValid.push_back(nextContigIndex[j]);
						//if(nextContigIndex[j]._pos == currentContigIndex[i]._pos+1){
						//	sharedContigIndexValid.push_back(nextContigIndex[j]);
						//}

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


	long longestPrefixSuffix(const Contig& contig)
	{
		int n = contig._minimizers.size();
	
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
			if (contig._minimizers[i] == contig._minimizers[len])
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



};	


#endif 


