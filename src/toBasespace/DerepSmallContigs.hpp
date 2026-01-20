




#ifndef MDBG_METAG_DEREPSMALLCONTIGS
#define MDBG_METAG_DEREPSMALLCONTIGS

#include "../Commons.hpp"


class DerepSmallContigs : public Tool{
    
public:

	struct MinimizerPosition{
		MinimizerType _minimizer;
		u_int32_t _contigIndex;
		u_int32_t _position;
	};

	vector<MinimizerPosition> _allMinimizerPositions;
	ankerl::unordered_dense::map<MinimizerType, pair<u_int64_t, u_int64_t>> _minimizerLookupTable;
	

	typedef phmap::parallel_flat_hash_map<KmerVec, vector<u_int32_t>, phmap::priv::hash_default_hash<KmerVec>, phmap::priv::hash_default_eq<KmerVec>, std::allocator<std::pair<KmerVec, vector<u_int32_t>>>, 4, std::mutex> KminmerReadMap;
	typedef phmap::parallel_flat_hash_map<ReadType, DnaBitset2*, phmap::priv::hash_default_hash<ReadType>, phmap::priv::hash_default_eq<ReadType>, std::allocator<std::pair<ReadType, vector<DnaBitset2*>>>, 4, std::mutex> BinaryReadMap;


	struct MinimizerPairPosition{

		ReadType _readIndex;
		u_int32_t _positionIndex;
	};

	typedef phmap::parallel_flat_hash_map<KmerVec, vector<MinimizerPairPosition>, phmap::priv::hash_default_hash<KmerVec>, phmap::priv::hash_default_eq<KmerVec>, std::allocator<std::pair<KmerVec, vector<MinimizerPairPosition>>>, 4, std::mutex> KminmerPosMap;

	struct Anchor{

		u_int32_t _referenceIndex;
		u_int32_t _referencePosition;
		u_int32_t _queryPosition;
		//bool _isReversed;
		//int32_t _referencePositionIndex;
		//int32_t _queryPositionIndex;

		Anchor(const u_int32_t referenceIndex, const u_int32_t referencePosition, const u_int32_t queryPosition){
			_referenceIndex = referenceIndex;
			_referencePosition = referencePosition;
			_queryPosition = queryPosition;
			//_isReversed = isReversed;
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
	};

	string _inputDir;
	int _nbCores;
	string _inputFilenameContig;
	string _outputFilename;
	int _startK;
	int _endK;

	ofstream _outputFile;

	DerepSmallContigs(): Tool (){

	}


	void parseArgs(int argc, char* argv[]){

		args::ArgumentParser parser("toBasespace", ""); //"This is a test program.", "This goes after the options."
		args::Positional<std::string> arg_outputDir(parser, "outputDir", "", args::Options::Required);
		args::ValueFlag<int> arg_startK(parser, "", "", {"first-k"}, args::Options::Required);
		args::ValueFlag<int> arg_endK(parser, "", "", {"last-k"}, args::Options::Required);
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
		_startK = args::get(arg_startK);
		_endK = args::get(arg_endK);
		_nbCores = args::get(arg_nbCores);
		
		_inputFilenameContig = _inputDir + "/contig_data_init.txt";
		_outputFilename = _inputDir + "/contig_data_init_small.txt";

		openLogFile(_inputDir);

	}



    void execute (){


		Logger::get().debug() << "Loading contigs";
		loadContigs();
		//loadReadsTest();

		//cout << _allMinimizerPositions.size() << endl;
		//cout << "Nb minimizers: " << _allMinimizerPositions.size() << endl;


		//std::sort(_allMinimizerPositions.begin(), _allMinimizerPositions.end(), [](const MinimizerPosition& a, const MinimizerPosition& b){
		//	return a._minimizer < b._minimizer;
		//});
		Logger::get().debug() << "Sort index";
		sortIndex();

		//cout << _allMinimizerPositions.size() << endl;
		//cout << "todo: enelver les minimizer de la table des positions: _allMinimizerPositions" << endl;
		//cout << "Attention: if(anchors.size() < 30) return; a reduire a 3" << endl;

		//cout << "Creating minimizer lookup table" << endl;
		Logger::get().debug() << "Create lookup table";
		createLookupTable();

		_outputFile = ofstream(_outputFilename);

		Logger::get().debug() << "Aligning reads";
		alignSmallContigs();
		

		Logger::get().debug() << "Append long contigs";
		appendLongContigs();
		
		_outputFile.close();
	}

	void appendLongContigs(){

		KminmerParserParallel parser(_inputFilenameContig, 0, 0, false, false, 1);
		parser.parseSequences(AppendLongContigsFunctor(*this));

	}


	class AppendLongContigsFunctor {

		public:

		DerepSmallContigs& _parent;

		AppendLongContigsFunctor(DerepSmallContigs& parent) : _parent(parent){
		}

		AppendLongContigsFunctor(const AppendLongContigsFunctor& copy) : _parent(copy._parent){
		}

		~AppendLongContigsFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {


			u_int32_t contigSize = kminmerList._readMinimizers.size();
			_parent._outputFile.write((const char*)&contigSize, sizeof(contigSize));
			
			u_int8_t isCircular = kminmerList._isCircular;
			_parent._outputFile.write((const char*)&isCircular, sizeof(isCircular));
			_parent._outputFile.write((const char*)&kminmerList._readMinimizers[0], contigSize*sizeof(MinimizerType));

		}
	};

	void createLookupTable(){
		u_int64_t index = -1;
		MinimizerType minimizer = -1;
		MinimizerType lastInsertedMinimizer = -1;

		for(size_t i=0; i<_allMinimizerPositions.size(); i++){

			if(minimizer != _allMinimizerPositions[i]._minimizer){

				if(index != -1){

					u_int64_t startIndex = index;
					u_int64_t endIndex = i-1;
					_minimizerLookupTable[minimizer] = {startIndex, endIndex};
					lastInsertedMinimizer = minimizer;
					//cout << "Found chunk: " << startIndex << " -> " << endIndex << endl;
				}

				minimizer = _allMinimizerPositions[i]._minimizer;
				index = i;
			}
		}

		if(lastInsertedMinimizer != minimizer){
			u_int64_t startIndex = index;
			u_int64_t endIndex = _allMinimizerPositions.size()-1;
			_minimizerLookupTable[minimizer] = {startIndex, endIndex};
			
			//cout << "Found chunk last: " << startIndex << " -> " << endIndex << endl;
		}

		u_int64_t nbPositions = 0;

		for(const auto& it : _minimizerLookupTable){
			for(size_t i=it.second.first; i <= it.second.second; i++){
				nbPositions += 1;
			}
		}

		//cout << "Lookup table check: " << _allMinimizerPositions.size() << " " << nbPositions << endl;
		if(nbPositions != _allMinimizerPositions.size()){
			cout << "Lookup table issue" << endl;
			getchar();
		}
	}

	void loadContigs(){
		KminmerParserParallel parser(_inputFilenameContig, 0, 0, false, false, 1);
		parser.parseSequences(LoadContigsFunctor(*this));
	}

	class LoadContigsFunctor {

		public:

		DerepSmallContigs& _parent;

		LoadContigsFunctor(DerepSmallContigs& parent) : _parent(parent){
		}

		LoadContigsFunctor(const LoadContigsFunctor& copy) : _parent(copy._parent){
		}

		~LoadContigsFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {

			u_int32_t readIndex = kminmerList._readIndex;
			if(readIndex % 100000 == 0) Logger::get().debug() << readIndex;


			for(u_int32_t i=0; i<kminmerList._readMinimizers.size(); i++){
				
				_parent._allMinimizerPositions.push_back({kminmerList._readMinimizers[i], readIndex, i});
			}


		}
	};


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
				return a._minimizer < b._minimizer;
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
						return a._minimizer < b._minimizer;
					});
			}
		}
		
		//cout << "done" << endl;
	}

	void alignSmallContigs(){


		for(size_t i=_startK; i<=_endK; i++){
			const string smallContigFilename = _inputDir + "/smallContigs/smallContigs_k" + to_string(i) + ".bin";
			//cout << smallContigFilename << " " << fs::exists(smallContigFilename) << endl;
			if(!fs::exists(smallContigFilename)) continue;

			KminmerParserParallel parser(smallContigFilename, 0, 0, false, false, _nbCores);
			parser.parseSequences(AlignSmallContigsFunctor(*this));
		}

	}

	
	class AlignSmallContigsFunctor {

		public:

		DerepSmallContigs& _parent;
		ReadMapping2 _bestReadMapping;

		AlignSmallContigsFunctor(DerepSmallContigs& parent) : _parent(parent){
		}

		AlignSmallContigsFunctor(const AlignSmallContigsFunctor& copy) : _parent(copy._parent){
		}

		~AlignSmallContigsFunctor(){
		}



		void operator () (const KminmerList& kminmerList) {
			

	
			//if(kminmerList._readIndex % 100000 == 0){
			//	cout << "Align small contigs: " << kminmerList._readIndex << endl;
			//}

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

			//_parent.alignRead(read, _minimizerAligner);
			mapRead(read, false);

			const MinimizerRead& read_rc = read.toReverseComplement();
			mapRead(read_rc, true);

			writeSmallContig(kminmerList._readIndex, kminmerList._readMinimizers);
			//writeAlignment2();
			
			//if(_bestReadMapping._contigIndex != -1){
				//cout << "Contig: " << kminmerList._readIndex << "\t" << kminmerList._readMinimizers.size() << "\t" << _bestReadMapping._matchScore << endl;	
				
			//}
		}

		void writeSmallContig(const ReadType readIndex, const vector<MinimizerType>& minimizers){

			//cout << "Contig: " << readIndex << "\t" << minimizers.size() << "\t" << _bestReadMapping._matchScore << endl;	
				
			if(_bestReadMapping._contigIndex != -1 && ((minimizers.size() - _bestReadMapping._matchScore) <= 2)) return;
			
			
			#pragma omp critical(writeSmallContig)
			{
				//cout << "Write! " << endl;
				u_int32_t contigSize = minimizers.size();
				_parent._outputFile.write((const char*)&contigSize, sizeof(contigSize));
				
				u_int8_t isCircular = CONTIG_LINEAR;
				_parent._outputFile.write((const char*)&isCircular, sizeof(isCircular));
				_parent._outputFile.write((const char*)&minimizers[0], contigSize*sizeof(MinimizerType));

			}
		}



		void mapRead(const MinimizerRead& read, const bool isReversed){

			//bool print = read._readIndex == 2479658;


			
			ReadType readIndex = read._readIndex;

			if(readIndex % 10000 == 0) Logger::get().debug() << readIndex;
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

		void processAnchors(const vector<Anchor>& anchors, const MinimizerRead& queryRead, const u_int32_t contigIndex, const bool isReversed){
			if(anchors.size() < 3) return;

			/*
			cout << endl;
			for(size_t i=0; i<anchors.size(); i++){
				const Anchor& anchor = anchors[i];
				cout << i << "\t" << anchor._referenceIndex << "\t" << anchor._referencePosition << "\t" << anchor._queryPosition << "\t" << queryRead._minimizersPos[anchor._queryPosition] << endl;

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

			writeAlignment(chain, queryRead, contigIndex, isReversed);
			//getchar();
		}



		Chain chainAnchors(const vector<Anchor>& anchors, const int64_t& maxChainingBand, const float w, const MinimizerRead& queryRead) {

			bool print = false; //queryRead._readIndex == 1118122;
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

			if(interval.size() < 3) return chain;

			std::reverse(interval.begin(), interval.end());

			chain._chainingScore = maxScore;
			chain._anchorInterval = interval;
			chain._nbMatches = interval.size();

			const Anchor& firstAnchor = anchors[interval[0]];
			const Anchor& lastAnchor = anchors[interval[interval.size()-1]];

			chain._nbDifferencesQuery = (lastAnchor._queryPosition - firstAnchor._queryPosition + 1) - chain._nbMatches;
			chain._nbDifferencesReference = (lastAnchor._referencePosition - firstAnchor._referencePosition + 1) - chain._nbMatches;
			chain._queryStart = firstAnchor._queryPosition;
			chain._queryEnd = lastAnchor._queryPosition;
			chain._referenceStart = firstAnchor._referencePosition;
			chain._referenceEnd = lastAnchor._referencePosition;

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
				//cout << "\t" << i << " " << k << " " << xi._isReversed << " " << xj._isReversed << " " << xi._referencePosition << "-" << xi._queryPosition << " " << xj._referencePosition << "-" << xj._queryPosition << endl;


				int32_t anchorSpacing = queryRead._minimizersPos[xi._queryPosition] - queryRead._minimizersPos[xj._queryPosition];
				if(anchorSpacing > 5000) continue;

				int32_t d_q;// = abs(xi._queryPosition - xj._queryPosition);
				//if(xi._isReversed){
				//	d_q = xj._queryPosition - xi._queryPosition;
				//}
				//else{
				d_q = xi._queryPosition - xj._queryPosition;
				//}
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
				
				//if(xi._isReversed){
				//	if(xi._queryPosition > xj._queryPosition) continue;
				//}
				//else{
				if(xi._queryPosition < xj._queryPosition) continue;
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

		void writeAlignment(const Chain& chain, const MinimizerRead& queryRead, const u_int32_t contigIndex, const bool isReversed){

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

			int32_t overhangStart = 0; //queryRead._minimizersPos[readStart] / _parent._averageDistanceBetweenMinimizers;
			int32_t overhangEnd = 0; //(readLengthBp-queryRead._minimizersPos[readEnd]) / _parent._averageDistanceBetweenMinimizers;

			int32_t matchScore = nbMatches; //  - overhangStart - overhangEnd; //- nbDifferences



			ReadMapping2 readMapping = {readIndex, contigIndex, readStart, readEnd, contigStart, contigEnd, isReversed, matchScore, queryRead._minimizersPos[readStart], queryRead._minimizersPos[readEnd], queryRead._readLength};
			


			if(_bestReadMapping._contigIndex == -1){
				_bestReadMapping = readMapping;	
			}
			else{
				if(readMapping._matchScore > _bestReadMapping._matchScore){
					_bestReadMapping = readMapping;
				}
			}
		}

	};





};	


#endif 

