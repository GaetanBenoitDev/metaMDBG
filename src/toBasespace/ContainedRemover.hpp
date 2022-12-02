

#ifndef MDBG_METAG_ContainedRemover
#define MDBG_METAG_ContainedRemover

#include "../Commons.hpp"


class ContainedRemover {
  

public:

	string _inputDir;
	string _inputFilenameContig;
	size_t _kminmerSize;
	ofstream& _logFile;
	int _nbCores;

	ContainedRemover(const string& inputDir, const string& inputFilenameContig, size_t kminmerSize, ofstream& logFile, int nbCores) : _logFile(logFile){
		_inputDir = inputDir;
		_inputFilenameContig = inputFilenameContig;
		_kminmerSize = kminmerSize;
		_nbCores = nbCores;
	}

	void execute(){

		cerr << "Removing small contained contigs" << endl;

		indexSmallContigsNew();
		findDuplicateContigs();


		/*

		//_mdbgFirst = new MDBG(_kminmerSizeFirst);
		//_mdbgFirst->load(_outputDir + "/kminmerData_min_init.txt", false);

		//const string& smallContigFilename = _outputDir + "/small_contigs.bin";
		//KminmerParserParallel parser(_filename_output, _minimizerSize, _kminmerSizeFirst, false, false);
		//auto fp = std::bind(&ToMinspace::loadContigs_min_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
		//parser.parseMinspace(fp);

		//KminmerParserParallel parser(smallContigFilename, _minimizerSize, _kminmerSizeFirst, false, false, _nbCores);
		//parser.parse(LoadContigKminmerFunctor(*this));
		KminmerParser parser2(_filename_smallContigs, _minimizerSize, _kminmerSizeFirst, false, false);
		auto fp2 = std::bind(&CreateMdbg::loadSmallContigs_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
		parser2.parseSequences(fp2);
		//KminmerParser parser(_filename_smallContigsNew, _minimizerSize, _kminmerSizeFirst, false, false);
		//auto fp = std::bind(&CreateMdbg::loadSmallContigs_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
		//parser.parseSequences(fp);

		std::sort(_contigs.begin(), _contigs.end(), ContigComparator_ByLength);

		for(size_t i=0; i<_contigs.size(); i++){	
			for(size_t j=0; j<_contigs.size(); j++){	
				if(i == j) continue;
				bool isDuplicated = false;
				//if(_invalidContigIndex.find(_contigs[j]._readIndex) != _invalidContigIndex.end()) continue;

				if(isSubArray(_contigs[i]._minimizers, _contigs[j]._minimizers)){
					isDuplicated = true;
					//cout << "lala" << isSubArray(_contigs[i]._minimizers, _contigs[j]._minimizers) << endl;
					//cout << "Is sub array: " << _contigs[i]._minimizers.size() << " " << _contigs[j]._minimizers.size() << endl;

					_invalidContigIndex.insert(_contigs[j]._readIndex);
					_invalidContigIndexDebugPrev[_contigs[j]._readIndex] = _contigs[j]._minimizers;
					_invalidContigIndexDebugNew[_contigs[j]._readIndex] = _contigs[i]._minimizers;


				}

				
				std::reverse(_contigs[j]._minimizers.begin(), _contigs[j]._minimizers.end());

				if(isSubArray(_contigs[i]._minimizers, _contigs[j]._minimizers)){
					isDuplicated = true;
					//cout << "lala" << isSubArray(_contigs[i]._minimizers, _contigs[j]._minimizers) << endl;
					//cout << "Is sub array rev: " << _contigs[i]._minimizers.size() << " " << _contigs[j]._minimizers.size() << endl;

					_invalidContigIndex.insert(_contigs[j]._readIndex);
					_invalidContigIndexDebugPrev[_contigs[j]._readIndex] = _contigs[j]._minimizers;
					_invalidContigIndexDebugNew[_contigs[j]._readIndex] = _contigs[i]._minimizers;

				}
				std::reverse(_contigs[j]._minimizers.begin(), _contigs[j]._minimizers.end());


			}
		}	
		*/

		//_logFile << "Truth: " << _invalidContigIndex.size() << endl;
		//cout << "Nb small contigs: " << _nbSmallContigs << endl;
		cout << "Nb contained contigs: " << _duplicatedSmallContigs.size() << endl;

		/*
		for(u_int32_t contigIndex : _invalidContigIndex){


			if(_duplicatedSmallContigs.find(contigIndex) == _duplicatedSmallContigs.end()){
				cout << "woot" << endl;

				for(u_int64_t m : _invalidContigIndexDebugNew[contigIndex]){
					cout << m << " ";
				}
				cout << endl;
				for(u_int64_t m : _invalidContigIndexDebugPrev[contigIndex]){
					cout << m << " ";
				}
				cout << endl;
				getchar();
			}
			
		}

		
		for(u_int32_t contigIndex : _duplicatedSmallContigs){

			
			if(_invalidContigIndex.find(contigIndex) == _invalidContigIndex.end()){
				cout << "woot" << endl;

				for(u_int64_t m : _duplicatedSmallContigsDebug[contigIndex]){
					cout << m << " ";
				}
				cout << endl;
				for(u_int64_t m : _duplicatedSmallContigsDebugPrev[contigIndex]){
					cout << m << " ";
				}
				cout << endl;
				getchar();
			}
			
		}
		*/
		
		



		dumpDereplicatedSmallContigs();



		_logFile << "done" << endl;
		//getchar();
		/*
		delete _mdbgFirst;

		for(auto& it : _nodeName_to_contigs){
			std::sort(it.second.begin(), it.second.end());
		}

		#pragma omp parallel num_threads(_nbCores)
		{

			#pragma omp for
			for(size_t i=0; i<_contigs.size(); i++){

				const Contig& contig = _contigs[i];

				//if()
				//u_int64_t length = _kminmerLengthMean + ((contig._nodepath.size()-1) * (_kminmerLengthMean-_kminmerOverlapMean));
				//if(length > 200000) continue; //todo use length


				bool sharing = true;

				vector<u_int32_t> sharingContigIndexes = _nodeName_to_contigs[contig._nodepath[0]];

				for(size_t j=0; j<contig._nodepath.size(); j++){

					u_int32_t nodeName = contig._nodepath[j];

					vector<u_int32_t> sharingContigIndexesTmp = _nodeName_to_contigs[nodeName];
					//cout << sharingContigIndexesTmp.size() << endl;
					vector<u_int32_t> intersection;
					//std::sort(v1.begin(), v1.end());
					//std::sort(v2.begin(), v2.end());

					std::set_intersection(sharingContigIndexes.begin(),sharingContigIndexes.end(),
										sharingContigIndexesTmp.begin(),sharingContigIndexesTmp.end(),
										back_inserter(intersection));

					if(intersection.size() <= 1){
						sharing = false;
						break;
					}

					sharingContigIndexes = intersection;
				}

				if(sharing){
					cout << "TTT duplicate: " << contig._nodepath.size() << endl;
					
					#pragma omp critical(checkDuplication)
					{
						_invalidContigIndex.insert(contig._readIndex);
					}
				}

			}
		}
	
		rewriteContigs();
		*/
	}




	struct Contig{
		u_int64_t _readIndex;
		vector<u_int64_t> _minimizers;
		//vector<u_int32_t> _nodepath;
		//vector<u_int32_t> _nodepath_sorted;
	};

	vector<Contig> _contigs;
	//vector<vector<u_int64_t>> _contigsDebug;
	

	struct ContigMinimizerIndex{
		u_int32_t _contigIndex;
		u_int32_t _pos;
	};

	//phmap::parallel_flat_hash_map<u_int64_t, vector<ContigMinimizerIndex>> _contigMinimizerIndex;
	typedef phmap::parallel_flat_hash_map<u_int64_t, vector<ContigMinimizerIndex>, phmap::priv::hash_default_hash<u_int64_t>, phmap::priv::hash_default_eq<u_int64_t>, std::allocator<std::pair<u_int64_t, vector<ContigMinimizerIndex>>>, 4, std::mutex> ContigMinimizerMap;

	ContigMinimizerMap _contigMinimizerIndex;
	
	void indexSmallContigsNew(){
		//KminmerParserParallel parser(_filename_smallContigs, _minimizerSize, _kminmerSizeFirst, false, false);
		//auto fp = std::bind(&CreateMdbg::indexSmallContigsNew_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
		//parser.parseSequences(fp);
		
		_logFile << "Removing contained contigs" << endl;
		KminmerParserParallel parser(_inputFilenameContig, -1, _kminmerSize, false, false, _nbCores);
		parser.parseSequences(IndexSmallContigsNewFunctor(*this));
	}

	class IndexSmallContigsNewFunctor {

		public:

		ContainedRemover& _graph;

		IndexSmallContigsNewFunctor(ContainedRemover& graph) : _graph(graph){

		}

		IndexSmallContigsNewFunctor(const IndexSmallContigsNewFunctor& copy) : _graph(copy._graph){

		}

		~IndexSmallContigsNewFunctor(){
			//delete _minimizerParser;
		}

		void operator () (const KminmerList& kminmerList) {


			u_int64_t readIndex = kminmerList._readIndex;
			vector<u_int64_t> readMinimizers = kminmerList._readMinimizers;

			for(size_t i=0; i<readMinimizers.size(); i++){

				u_int64_t m = readMinimizers[i];

				_graph._contigMinimizerIndex.lazy_emplace_l(m, 
				[&readIndex, &i](ContigMinimizerMap::value_type& v) { // key exist
					v.second.push_back({readIndex, i});
				},           
				[&m, &readIndex, &i](const ContigMinimizerMap::constructor& ctor) { // key inserted
					
					vector<ContigMinimizerIndex> lala;
					lala.push_back({readIndex, i});
					ctor(m, lala); 

				}); // construct value_type in place when key not present

				//_contigMinimizerIndex[readMinimizers[i]].push_back({readIndex, i});
			}
			//_contigsDebug.push_back(readMinimizers);
		}
	};

	void findDuplicateContigs(){
		
		_logFile << "Determine duplicated small contigs" << endl;
		KminmerParserParallel parser(_inputFilenameContig, 1, _kminmerSize, false, false, _nbCores);
		parser.parseSequences(FindDuplicateContigsFunctor(*this));
	}


	phmap::parallel_flat_hash_set<u_int32_t> _duplicatedSmallContigs;
	phmap::parallel_flat_hash_map<u_int32_t, vector<u_int64_t>> _duplicatedSmallContigsDebug;
	phmap::parallel_flat_hash_map<u_int32_t, vector<u_int64_t>> _duplicatedSmallContigsDebugPrev;

	class FindDuplicateContigsFunctor {

		public:

		ContainedRemover& _graph;

		struct SubseqIndex{
			u_int64_t _contigIndex;
			u_int64_t _pos;
			u_int64_t _matchLength;
			bool _isValid;
		};

		FindDuplicateContigsFunctor(ContainedRemover& graph) : _graph(graph){

		}

		FindDuplicateContigsFunctor(const FindDuplicateContigsFunctor& copy) : _graph(copy._graph){

		}

		~FindDuplicateContigsFunctor(){
			//delete _minimizerParser;
		}




		void operator () (const KminmerList& kminmerList) {

			u_int64_t readIndex = kminmerList._readIndex;
			vector<u_int64_t> readMinimizers = kminmerList._readMinimizers;
			/*
			const vector<u_int64_t> target = {69637124315931369, 84477587394475608, 78143134765766670, 3439892788480131, 445750117931955, 445750117931955, 87646769359969727};
			if(readMinimizers == target){
				cout << "lala" << endl;
				cout << isDuplicated(readMinimizers) << endl;
				getchar();
			}
			*/
			
			if(isDuplicated(readIndex, readMinimizers, 1)){
				#pragma omp critical(checkDuplication)
				{
					_graph._duplicatedSmallContigs.insert(readIndex);
				}
				return;
			}

			//std::reverse(readMinimizers.begin(), readMinimizers.end());

			if(isDuplicated(readIndex, readMinimizers, -1)){
				#pragma omp critical(checkDuplication)
				{
					_graph._duplicatedSmallContigs.insert(readIndex);
				}
				return;
			}
			
		}

		
		bool isDuplicated(u_int64_t contigIndex, const vector<u_int64_t>& readMinimizers, int dir){
			
			//u_int64_t minOverlap = readMinimizers.size() - _graph._kminmerSizeF +1;

			//cout << "-----" << readMinimizers.size() << endl;
			bool isTarget = false;
			vector<SubseqIndex> matches;
			//vector<SubseqIndex> matches;
			//unordered_set<u_int32_t> matched

			for(size_t i=0; i<readMinimizers.size(); i++){
				u_int64_t m = readMinimizers[i];


				const vector<ContigMinimizerIndex>& currentMatches = _graph._contigMinimizerIndex[m];


				if(i==0){
					
					for(const ContigMinimizerIndex& currentMatch : currentMatches){
						if(currentMatch._contigIndex == contigIndex) continue;
						matches.push_back({currentMatch._contigIndex, currentMatch._pos, 1});
					}

					continue;
				}


				for(SubseqIndex& match : matches){
					match._isValid = false;
				}

				for(SubseqIndex& match : matches){
					for(const ContigMinimizerIndex& currentMatch : currentMatches){
						if(currentMatch._contigIndex == match._contigIndex){
							if(match._pos+dir == currentMatch._pos){
								match._pos += dir;
								match._matchLength += 1;
								match._isValid = true;

								if(match._matchLength >= readMinimizers.size()){
									/*
									#pragma omp critical(checkDuplication)
									{
										_graph._duplicatedSmallContigsDebug[contigIndex] = _graph._contigsDebug[match._contigIndex];
										_graph._duplicatedSmallContigsDebugPrev[contigIndex] = readMinimizers;


										//if(isTarget){
										//	for(u_int64_t m : _graph._contigsDebug[match._contigIndex]){
										//		cout << m << " ";
										//	}
										//	cout << endl;
										//}
									}
									*/

									return true;
								}

								break;
							}
						}
						
						//matches.push_back({currentMatch._contigIndex, currentMatch._pos, 1});
					}
				}

				vector<SubseqIndex> matchesTmp = matches;
				matches.clear();
				for(const SubseqIndex& match : matchesTmp){
					if(!match._isValid) continue;
					matches.push_back(match);
				}

				//cout << matches.size() << endl;
				if(matches.size() == 0) return false;

			}

			/*
			vector<SubseqIndex> successorContigIndexes;
			vector<SubseqIndex> nonSuccessorContigIndexes;

			//const vector<u_int64_t> target = {72097720224774221, 30942406196105332, 844230355215342, 83307865189328548, 72097720224774221};
			const vector<u_int64_t> target = {54786101128205186, 49710951260278902, 50123046530649727, 82839499731899160, 85030874348372966};
			bool isTarget = readMinimizers == target;

			if(isTarget){
				cout << "lala " << contigIndex << " " << dir << endl;
			}
			phmap::parallel_flat_hash_map<u_int32_t, vector<SubseqIndex>> matches;

			for(size_t i=0; i<readMinimizers.size(); i++){
				u_int64_t m = readMinimizers[i];

				if(isTarget){
					cout << m << endl;
				}

				const vector<ContigMinimizerIndex>& currentMatches = _graph._contigMinimizerIndex[m];

				bool isSuccessor = false;
				successorContigIndexes.clear();
				nonSuccessorContigIndexes.clear();

				for(const ContigMinimizerIndex& currentMatch : currentMatches){

					if(i == 0){

						if(currentMatch._contigIndex == contigIndex) continue;
						if(isTarget){

							vector<u_int64_t> lala = readMinimizers;
							std::reverse(lala.begin(), lala.end());
							if(_graph.isSubArray(_graph._contigsDebug[currentMatch._contigIndex], readMinimizers)){
								cout << "\t" << currentMatch._contigIndex << " " << currentMatch._pos << endl;

								cout << "\t";
								for(u_int64_t m : _graph._contigsDebug[currentMatch._contigIndex]){
									cout << m << " ";
								}
								cout << endl;
							}
							else if(_graph.isSubArray(_graph._contigsDebug[currentMatch._contigIndex], lala)){
								cout << "\t" << currentMatch._contigIndex << " " << currentMatch._pos << endl;

								cout << "\t";
								for(u_int64_t m : _graph._contigsDebug[currentMatch._contigIndex]){
									cout << m << " ";
								}
								cout << endl;
							} 
							


						}

						//if(dir && matches.find(currentMatch._contigIndex) != matches.end()) continue; //if sub array from left to right we keep the first pos

						matches[currentMatch._contigIndex].push_back({currentMatch._contigIndex, currentMatch._pos, 1});

						continue;
					}

					if(matches.find(currentMatch._contigIndex) == matches.end()){
						//matches.erase(currentMatch._contigIndex);
						//matches[currentMatch._contigIndex] = {currentMatch._contigIndex, i, 1};
						continue;
					}
					

					vector<SubseqIndex>& existingMatches = matches[currentMatch._contigIndex];

					for(SubseqIndex& existingMatch : existingMatches){
							
						if(isTarget){
							cout << "\t" << existingMatch._contigIndex << " " << existingMatch._pos << " " << existingMatch._matchLength << "    " << currentMatch._pos << " " << (existingMatch._pos+dir == currentMatch._pos) << endl;
						}

						//if(existingMatch._pos+1 == i){

						
						if(existingMatch._pos+dir == currentMatch._pos){
							
							existingMatch._pos += dir;
							existingMatch._matchLength += 1;

							successorContigIndexes.push_back(existingMatch);
							if(isTarget){
								//cout << "\t match" << endl;
							}

							//existingMatch._pos = i;
							//existingMatch._matchLength += 1;

							
							if(existingMatch._matchLength == readMinimizers.size()){

								#pragma omp critical(checkDuplication)
								{
									_graph._duplicatedSmallContigsDebug[contigIndex] = _graph._contigsDebug[existingMatch._contigIndex];
									_graph._duplicatedSmallContigsDebugPrev[contigIndex] = readMinimizers;

									if(isTarget){
										for(u_int64_t m : _graph._contigsDebug[existingMatch._contigIndex]){
											cout << m << " ";
										}
										cout << endl;
									}
								}

								return true;
							}
							
						}
						else{
							//nonSuccessorContigIndexes.push_back(existingMatch._contigIndex);

							//existingMatch._pos = i;
							//existingMatch._matchLength = 1;
						}
					}

				}

				
				//for(u_int32_t successorContigIndex : successorContigIndexes){
				for(SubseqIndex& existingMatchSuccessor : successorContigIndexes){
					
					//vector<SubseqIndex>& matchess = ;

					for(SubseqIndex& existingMatch : matches[existingMatchSuccessor._contigIndex]){

						if(existingMatch._pos != existingMatchSuccessor._pos) continue;

						//SubseqIndex& existingMatch = matches[successorContigIndex];

						existingMatch._pos += dir;
						existingMatch._matchLength += 1;

						if(existingMatch._matchLength == readMinimizers.size()){

							#pragma omp critical(checkDuplication)
							{
								_graph._duplicatedSmallContigsDebug[contigIndex] = _graph._contigsDebug[existingMatch._contigIndex];
								_graph._duplicatedSmallContigsDebugPrev[contigIndex] = readMinimizers;

								if(isTarget){
									for(u_int64_t m : _graph._contigsDebug[existingMatch._contigIndex]){
										cout << m << " ";
									}
									cout << endl;
								}
							}
							if(isTarget)
							cout << "yay " << existingMatch._contigIndex << " " << existingMatch._matchLength << endl;
							return true;
						}
					}
				}
				
				//handle repeated minimizer in readMinimizers
				for(u_int32_t successorContigIndex : nonSuccessorContigIndexes){
					
					SubseqIndex& existingMatch = matches[successorContigIndex];
					
					if(std::find(successorContigIndexes.begin(), successorContigIndexes.end(), successorContigIndex) == successorContigIndexes.end()){
						//existingMatch._pos = -1;
						//existingMatch._matchLength = 0;
						matches.erase(successorContigIndex);
					}
				}
				

			}
			*/

			return false;
		}

	};


	phmap::parallel_flat_hash_map<u_int32_t, vector<u_int32_t>> _nodeName_to_contigs;

	phmap::parallel_flat_hash_set<u_int32_t> _invalidContigIndex;
	phmap::parallel_flat_hash_map<u_int32_t, vector<u_int64_t>> _invalidContigIndexDebugPrev;
	phmap::parallel_flat_hash_map<u_int32_t, vector<u_int64_t>> _invalidContigIndexDebugNew;

	bool isSubArray(const vector<u_int64_t>& A, const vector<u_int64_t>& B){

		size_t n = A.size();
		size_t m = B.size();

		// Two pointers to traverse the arrays
		int i = 0, j = 0;

		// Traverse both arrays simultaneously
		while (i < n && j < m) {

			// If element matches
			// increment both pointers
			if (A[i] == B[j]) {

				i++;
				j++;

				// If array B is completely
				// traversed
				if (j == m)
					return true;
			}
			// If not,
			// increment i and reset j
			else {
				i = i - j + 1;
				j = 0;
			}
		}

		return false;
	}

	static bool ContigComparator_ByLength(const Contig &a, const Contig &b){
		return a._minimizers.size() > b._minimizers.size();
	}

	bool _laodingPrev;
	vector<Contig> _contigsPrev;

	unordered_map<u_int32_t, vector<u_int64_t>> _dupLala;


	void loadSmallContigs_read(const vector<u_int64_t>& readMinimizers, bool isCircular, u_int64_t readIndex){
		_contigs.push_back({readIndex, readMinimizers});

		if(_laodingPrev){
			_contigsPrev.push_back({readIndex, readMinimizers});
		}
	}

	/*
	class LoadContigKminmerFunctor {

		public:

		CreateMdbg& _graph;

		LoadContigKminmerFunctor(CreateMdbg& graph) : _graph(graph){

		}

		LoadContigKminmerFunctor(const LoadContigKminmerFunctor& copy) : _graph(copy._graph){

		}

		~LoadContigKminmerFunctor(){
			//delete _minimizerParser;
		}



		void operator () (const KminmerList& kminmerList) {


			u_int64_t readIndex = kminmerList._readIndex;
			const vector<u_int64_t>& readMinimizers = kminmerList._readMinimizers;
			//const vector<KmerVec>& kminmers = kminmerList._kminmers;
			const vector<ReadKminmerComplete>& kminmersInfos = kminmerList._kminmersInfo;


			//_nbContigs += 1;

			#pragma omp critical
			{
				//vector<u_int32_t> nodepath_sorted = nodepath;
				//std::sort(nodepath_sorted.begin(), nodepath_sorted.end());
				_graph._contigs.push_back({readIndex, readMinimizers});
				//cout << "load contig: " << nodepath.size() << endl;
			}
			
		}
	};
	*/
	
	ofstream _tmpOutputFileSmallContigs;

	void dumpDereplicatedSmallContigs(){

		/*
		_fileSmallContigs = ofstream(_filename_smallContigsNew, std::ios_base::app);

		KminmerParser parser(_filename_smallContigsPrev, _minimizerSize, _kminmerSizeFirst, false, false);
		auto fp = std::bind(&CreateMdbg::rewriteContigs_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
		parser.parseSequences(fp);

		_fileSmallContigs.close();
		*/

		
		const string& tmpFilename = _inputFilenameContig + ".tmp";
		_tmpOutputFileSmallContigs = ofstream(tmpFilename);

		KminmerParser parser(_inputFilenameContig, 1, _kminmerSize, false, false);
		auto fp = std::bind(&ContainedRemover::rewriteContigs_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
		parser.parseSequences(fp);

		_tmpOutputFileSmallContigs.close();
		
		
		fs::remove(_inputFilenameContig);
		fs::rename(tmpFilename, _inputFilenameContig);
		
	}

	void rewriteContigs_read(const vector<u_int64_t>& readMinimizers, bool isCircular, u_int64_t readIndex){

		if(_duplicatedSmallContigs.find(readIndex) != _duplicatedSmallContigs.end()) return;

		
		u_int32_t contigSize = readMinimizers.size();
		//_logFile << "Write: " << contigSize << endl;
		//getchar();
		_tmpOutputFileSmallContigs.write((const char*)&contigSize, sizeof(contigSize));
		_tmpOutputFileSmallContigs.write((const char*)&isCircular, sizeof(isCircular));
		_tmpOutputFileSmallContigs.write((const char*)&readMinimizers[0], contigSize*sizeof(u_int64_t));

	}
	
};	


#endif 


