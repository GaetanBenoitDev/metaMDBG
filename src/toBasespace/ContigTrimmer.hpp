

#ifndef MDBG_METAG_CONTIGTRIMMER
#define MDBG_METAG_CONTIGTRIMMER

#include "../Commons.hpp"
#include "../graph/GfaParser.hpp"



class ContigTrimmer{
    
public:

	struct AlignmentShort{
		u_int32_t _contigIndex;
		u_int32_t _contigStart;
		u_int32_t _contigEnd;
		float _score;
	};

	string _inputContigFilename;
	string _usedReadFilename;
	string _outputContigFilename;
	string _tmpDir;
	string _readPartitionDir;
	int _nbCores;
	int _minimapBatchSize;
	string _alignFilename;
	u_int32_t _minContigLength;
	string _minimap2Preset_map;

	gzFile _outputContigFile;
	phmap::parallel_flat_hash_map<ReadType, vector<Alignment>> _allAlignments;
	phmap::parallel_flat_hash_map<ReadType, vector<Alignment>> _allAlignmentsTrue;
	phmap::parallel_flat_hash_map<u_int32_t, vector<pair<u_int32_t, u_int32_t>>> _alignments_readsVsContigs;
	phmap::parallel_flat_hash_map<ReadType, vector<AlignmentShort>> _allAlignments2;
	



	
	ContigTrimmer(const string& inputContigFilename, const string& usedReadFilename, const string& outputContigFilename, const int& minimapBatchSize, const string& tmpDir, const string& readPartitionDir, const u_int32_t& minContigLength, const int& nbCores, const string minimap2Preset_map){
		_inputContigFilename = inputContigFilename;
		_usedReadFilename = usedReadFilename;
		_outputContigFilename = outputContigFilename;
		_minimapBatchSize = minimapBatchSize;
		_tmpDir = tmpDir;
		_readPartitionDir = readPartitionDir;
		_minContigLength = minContigLength;
		_nbCores = nbCores;
		_minimap2Preset_map = minimap2Preset_map;

		_alignFilename = _readPartitionDir + "/align_usedReadsVsContigs.paf.gz";

	}


    void execute (){

		_nbRepeatitiveContigs = 0;
		Logger::get().debug() << "\tAligning reads vs contigs";
		alignReads();
		Logger::get().debug() << "\tDone: " << (peakrss() / 1024.0 / 1024.0 / 1024.0);

		//Logger::get().debug() << "\tLoading contig coverages";
		//loadContigCoverages();

		//loadAlignmentsTrue();

		Logger::get().debug() << "\tLoading alignments";
		loadAlignments();
		Logger::get().debug() << "\tDone: " << (peakrss() / 1024.0 / 1024.0 / 1024.0);
		
		//exit(1);
		Logger::get().debug() << "\tTrimming contigs";
		trimContigs();
		Logger::get().debug() << "\tDone: " << (peakrss() / 1024.0 / 1024.0 / 1024.0);
		
	}

	void alignReads(){

		float minimapBatchSizeC = max(0.5, _minimapBatchSize / 4.0);
		//    let args = ["-t", &opts.nb_threads.to_string(), "-c", "-xasm20", "-DP", "--dual=no", "--no-long-join", "-r100", "-z200", "-g2k", fasta_path_str, fasta_path_str];
		//string command = "minimap2 -v 0 -p 1 -c -I " + to_string(minimapBatchSizeC) + "G -K 0.02G -t " + to_string(_nbCores) + " -x " + _minimap2Preset_map + " " + _inputContigFilename + " " + _usedReadFilename;
		string command = "minimap2 -v 0 -p 1 -I " + to_string(minimapBatchSizeC) + "G -t " + to_string(_nbCores) + " -x " + _minimap2Preset_map + " " + _inputContigFilename + " " + _usedReadFilename;
		//cout << command << endl;
		//string command = "minimap2 -v 0 -X -m 500 -x asm20 -I " + to_string(minimapBatchSizeC) + "G -t " + to_string(_nbCores) + " " + _inputContigFilename + " " + _inputContigFilename;
		//Utils::executeMinimap2(command, _alignFilename);
		//cout << command << endl;
		//command += " | gzip -c - > " + _alignFilename;
		//Utils::executeCommand(command, _tmpDir);

		//cout << _minimapBatchSize << " " << minimapBatchSizeC << endl;

		mm_idxopt_t iopt;
		mm_mapopt_t mopt;
		//int n_threads = nbCores;

		mm_verbose = 2; // disable message output to stderr
		mm_set_opt(0, &iopt, &mopt);
		mm_set_opt(_minimap2Preset_map.c_str(), &iopt, &mopt); 
		//mopt.flag |= MM_F_CIGAR; // perform alignment
		iopt.batch_size = minimapBatchSizeC*1000000000ULL; //0x7fffffffffffffffL; //always build a uni-part index

		//mopt.min_chain_score = 500; //-m 500
		mopt.pri_ratio = 1; //-p 1

		//cout << mopt.zdrop << " " << mopt.bw << " " << mopt.max_gap << endl;
		//mopt.zdrop = mopt.zdrop_inv = 200; //-z
		//mopt.bw = 100; //-r
		//mopt.max_gap = 15*2; //-g

		// open query file for reading; you may use your favorite FASTA/Q parser
		//gzFile f = gzopen(_inputContigFilename.c_str(), "r");
		//assert(f);
		//kseq_t *ks = kseq_init(f);

		// open index reader
		mm_idx_reader_t *r = mm_idx_reader_open(_inputContigFilename.c_str(), &iopt, 0);

		mm_idx_t *mi;
		while ((mi = mm_idx_reader_read(r, _nbCores)) != 0) {

			Logger::get().debug() << "\tBuild index: " << (peakrss() / 1024.0 / 1024.0 / 1024.0);
			//cout << "Index pass" << endl;
			mm_mapopt_update(&mopt, mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!
			//mopt.mid_occ = 5;

			ReadParserParallel readParser(_usedReadFilename, true, false, _nbCores);
			readParser.parse(MapReadsFunctor(*this, mi, mopt));
			
			mm_idx_destroy(mi);
			Logger::get().debug() << "\tPass done: " << (peakrss() / 1024.0 / 1024.0 / 1024.0);
		}


		mm_idx_reader_close(r); // close the index reader

	}


	

	class MapReadsFunctor {

		public:

		ContigTrimmer& _parent;
		mm_idx_t* _mi;
		mm_tbuf_t*_tbuf;
		mm_mapopt_t& _mopt;

		MapReadsFunctor(ContigTrimmer& parent, mm_idx_t* mi, mm_mapopt_t& mopt) : _parent(parent), _mopt(mopt){
			_tbuf = mm_tbuf_init();
			_mi = mi;
		}

		MapReadsFunctor(const MapReadsFunctor& copy) : _parent(copy._parent), _mopt(copy._mopt){
			_tbuf = mm_tbuf_init();
			_mi = copy._mi;
		}

		~MapReadsFunctor(){
			mm_tbuf_destroy(_tbuf);
		}

		void operator () (const Read& read) {
			
			if(read._index % 100000 == 0) Logger::get().debug() << "\t\tAlign reads " << read._index;
			
			const string queryName = Utils::shortenHeader(read._header);
			//u_int32_t readIndex = stoull(read._header);

			mm_reg1_t *reg;
			int j, i, n_reg;
			reg = mm_map(_mi, read._seq.size(), read._seq.c_str(), &n_reg, _tbuf, &_mopt, 0); // get all hits for the query
			
			#pragma omp critical(loadAlignments)
			{
				for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
					mm_reg1_t *r = &reg[j];

					AlignmentBounds alignment;
					alignment._referenceLength = _mi->seq[r->rid].len;
					alignment._queryLength = read._seq.size();
					alignment._queryStart = r->qs;
					alignment._queryEnd = r->qe;
					alignment._referenceStart = r->rs;
					alignment._referenceEnd = r->re;
					alignment._isReversed = r->rev;
					alignment._identity = ((double)r->mlen) / r->blen;// r->div;
					alignment._nbMatches = r->mlen;

					_parent.loadAlignments_read(queryName, string(_mi->seq[r->rid].name), alignment, r->score);
					
					//_parent._checksum += 1;
					//cout << r->qs << " " << r->qe << endl;
					//assert(r->p); // with MM_F_CIGAR, this should not be NULL
					//printf("%s\t%d\t%d\t%d\t%c\t", ks->name.s, ks->seq.l, r->qs, r->qe, "+-"[r->rev]);
					//printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\tcg:Z:", mi->seq[r->rid].name, mi->seq[r->rid].len, r->rs, r->re, r->mlen, r->blen, r->mapq);
					//for (i = 0; i < r->p->n_cigar; ++i) // IMPORTANT: this gives the CIGAR in the aligned regions. NO soft/hard clippings!
					//	printf("%d%c", r->p->cigar[i]>>4, MM_CIGAR_STR[r->p->cigar[i]&0xf]);
					//putchar('\n');
					//free(r->p);
				}
			}

			free(reg);
			

			
			/*
			//if(_parent._allAlignments.find(readIndex) == _parent._allAlignments.end()) return;
			if(_parent._alignments.find(readIndex) == _parent._alignments.end()) return;

			//const vector<Alignment>& als = _parent._allAlignments[readIndex];

			//for(const Alignment& al : als){
			const Alignment& al = _parent._alignments[readIndex];
			//for(const Alignment& al : _alignments[readIndex]){
			u_int32_t contigIndex = al._contigIndex;

			if(_parent._contigSequences.find(contigIndex) == _parent._contigSequences.end()) return;
			*/
		}

	};

	/*

	void loadAlignmentsTrue(){

		PafParser pafParser(_readPartitionDir + "/true_align.paf.gz");
		auto fp = std::bind(&ContigTrimmer::loadAlignments_readTrue, this, std::placeholders::_1);
		pafParser.parse(fp);
		
	}


	void loadAlignments_readTrue(const string& line){

		const vector<string>& _fields = Utils::split(line, '\t');


		//u_int32_t queryContigIndex = contigNameToContigIndex(_fields[0]);
		u_int32_t contigIndex = contigNameToContigIndex(_fields[5]);
		ReadType readIndex = readNameToReadIndex(_fields[0]);


		//if(queryContigIndex == targetContigIndex) return;

		u_int64_t readLength = stoull((_fields)[1]);
		u_int32_t readStart = stoull((_fields)[2]);
		u_int32_t readEnd = stoull((_fields)[3]);
		u_int32_t contigLength = stoull((_fields)[6]);
		u_int32_t contigStart = stoull((_fields)[7]);
		u_int32_t contigEnd = stoull((_fields)[8]);

		u_int64_t nbMatches = stoull((_fields)[9]);
		u_int64_t alignLength = stoull((_fields)[10]);

		bool strand = (_fields)[4] == "-";

		float identity =  ((float)nbMatches) / alignLength;
	

		//if((float)identity < (float)_minIdentity) return;

		//if(queryLength < targetLength){
			//if(queryLength > 30000) return;
		//}
		//else{
			//if(targetLength > 30000) return;
		//}

		
		//float queryFractionCovered = (queryEnd-queryStart) / ((long double) queryLength);
		//float referenceFractionCovered = (targetEnd-targetStart) / ((long double) targetLength);

		Alignment alignment = {contigIndex, readIndex, strand, readStart, readEnd, contigStart, contigEnd, identity, readLength, contigLength}; //, score

		indexReadAlignmentTrue(readIndex, alignment);
	}

	
	void indexReadAlignmentTrue(const ReadType& readIndex, const Alignment& alignment){
		

		if(_allAlignmentsTrue.find(readIndex) == _allAlignmentsTrue.end()){
			_allAlignmentsTrue[readIndex].push_back(alignment);
			return;
		}

		vector<Alignment>& existingAlignments = _allAlignmentsTrue[readIndex];

		vector<Alignment>::iterator it = existingAlignments.begin();
		bool isBetterAlignment = false;
		bool hasOverlap = false;
		bool overlapWithBetterAlignment = false;
		//bool canAddAlignment = true;


		for(const Alignment& existingAlignment : existingAlignments){

			if(alignmentOverlapExistingAlignment(alignment, existingAlignment)){

				if(alignment.score() < existingAlignment.score()){
					overlapWithBetterAlignment = true;
				}

				hasOverlap = true;
			}

		}

		if(overlapWithBetterAlignment){
			return;
		}

		while(it != existingAlignments.end()) {

			const Alignment& existingAlignment = *it;



			if(alignmentOverlapExistingAlignment(alignment, existingAlignment) && alignment.score() > existingAlignment.score()) {
				it = existingAlignments.erase(it);
				isBetterAlignment = true;
			}
			else{
				++it;
			}
		}



		if(isBetterAlignment){
			_allAlignmentsTrue[readIndex].push_back(alignment);
		}

		if(!isBetterAlignment && !hasOverlap){
			_allAlignmentsTrue[readIndex].push_back(alignment);
		}
	}


	*/












	
	void loadAlignments(){

		//PafParser pafParser(_alignFilename);
		//PafParser pafParser(_readPartitionDir + "/true_align.paf.gz");
		//auto fp = std::bind(&ContigTrimmer::loadAlignments_read, this, std::placeholders::_1);
		//pafParser.parse(fp);
		
		/*
		float minIdentity = 2;
		for(const auto& it : _allAlignmentsTrue){
			for(const Alignment& al : it.second){
				if(al._identity < minIdentity){
					minIdentity = al._identity;
				}
			}

			if(it.second.size() > 1){
				cout << it.second.size() << " " << it.first << endl;
			}
		}

		exit(1);

		for(const auto& it : _allAlignments){
			const ReadType readIndex = it.first;

			const vector<Alignment>& alignmentsTrue = _allAlignmentsTrue[readIndex];

			for(const Alignment& al : it.second){

				bool alignExist = false;

				for(size_t i=0; i<alignmentsTrue.size(); i++){
					const Alignment& al2 = alignmentsTrue[i];
					if(al._readIndex == al2._readIndex && al._contigStart == al2._contigStart && al._contigEnd == al2._contigEnd && al._identity == al2._identity){
						alignExist = true;
					}
				}

				if(!alignExist){
					cout << "Align diff: " << al._readIndex << "\t" << al._contigIndex << "\t" << al.score() << endl;
					for(size_t i=0; i<alignmentsTrue.size(); i++){
						const Alignment& al2 = alignmentsTrue[i];
						cout << "\t" << al2._contigIndex << "\t" << al2.score() << endl;
					}
				}
				//_alignments_readsVsContigs[al._contigIndex].push_back({al._contigStart, al._contigEnd});
			}
		}

		cout << "Min identity: " << minIdentity << endl;
		
		for(const auto& it : _allAlignments){
			const ReadType readIndex = it.first;
		
			if(readIndex == 86917){
				for(const Alignment& al : it.second){
					cout << "lul: " << al._contigIndex << "\t" << al._readStart << "\t" << al._readEnd << endl;
				}

			}
		}



		*/

		for(const auto& it : _allAlignments2){
			const ReadType readIndex = it.first;

			unordered_map<u_int32_t, float> contigIndex_to_score;

			for(const AlignmentShort& al : it.second){

				if(contigIndex_to_score.find(al._contigIndex) == contigIndex_to_score.end()){
					contigIndex_to_score[al._contigIndex] = 0;
				}

				contigIndex_to_score[al._contigIndex] += al._score;
			}

			float maxScore = -1;
			u_int32_t maxContigIndex = -1;

			for(const auto& it : contigIndex_to_score){

				if(it.second > maxScore){
					maxScore = it.second;
					maxContigIndex = it.first;
				}
			}

			if(maxContigIndex == -1) continue;


			u_int32_t contigStart = -1;
			u_int32_t contigEnd = 0;

			for(const AlignmentShort& al : it.second){
				if(al._contigIndex == maxContigIndex){
					if(al._contigStart < contigStart) contigStart = al._contigStart;
					if(al._contigEnd > contigEnd) contigEnd = al._contigEnd;
				}
			}
			//if(it.second._contigIndex == 4375){
			//	cout << readIndex << endl;
			//}
			//for(const Alignment& al : it.second){

			//	_alignments_readsVsContigs[al._contigIndex].push_back({al._contigStart, al._contigEnd});
			//}
			_alignments_readsVsContigs[maxContigIndex].push_back({contigStart, contigEnd});
		}

		//cout << _allAlignments.size() << " " << _alignments_readsVsContigs.size() << endl;
		_allAlignments2.clear();
	
		//exit(1);
	}


	void loadAlignments_read(const string& queryName, const string& targetName, const AlignmentBounds& bounds, const float chainScore){

		//cout << queryName << " " << targetName << " " << chainScore << endl;
		//getchar();
		//const vector<string>& _fields = Utils::split(line, '\t');


		//u_int32_t queryContigIndex = contigNameToContigIndex(_fields[0]);
		u_int32_t contigIndex = contigNameToContigIndex(targetName);
		ReadType readIndex = readNameToReadIndex(queryName);
		/*
		const string& queryName = Utils::shortenHeader((_fields)[0]);
		const string& targetName = Utils::shortenHeader((_fields)[5]);
		

		string contigIndexField = fields[0];
		contigIndexField.erase(0,3); //"remove "ctg" letters

		contigHeader._contigIndex = stoull(contigIndexField);

		queryName.erase(0,3);

		targetName.erase(0,3);
		u_int32_t targetContigIndex = stoull(targetName);
		*/

		//if(queryContigIndex == targetContigIndex) return;

		u_int64_t readLength = bounds._queryLength; // stoull((_fields)[1]);
		u_int32_t readStart = bounds._queryStart;// stoull((_fields)[2]);
		u_int32_t readEnd = bounds._queryEnd; // stoull((_fields)[3]);
		u_int32_t contigLength = bounds._referenceLength; // stoull((_fields)[6]);
		u_int32_t contigStart = bounds._referenceStart; // stoull((_fields)[7]);
		u_int32_t contigEnd = bounds._referenceEnd; // stoull((_fields)[8]);

		//u_int64_t nbMatches = stoull((_fields)[9]);
		//u_int64_t alignLength = stoull((_fields)[10]);

		bool strand = bounds._isReversed; // (_fields)[4] == "-";

		float identity = bounds._identity; // ((float)nbMatches) / alignLength;
	
		/*
		float chainScore = 0;
		
		for(size_t i=12; i<_fields.size(); i++){


			const vector<string>& fields2 = Utils::split(_fields[i], ':');

			//if((*_fields_optional)[0] == "gi"){
			if(fields2[0] == "s1"){
				chainScore = std::stoi(fields2[2]);
				break;
			}

		}
		*/
		//cout << chainScore << endl;

		//s1:i:15632

		//if((float)identity < (float)_minIdentity) return;

		//if(queryLength < targetLength){
			//if(queryLength > 30000) return;
		//}
		//else{
			//if(targetLength > 30000) return;
		//}

		
		//float queryFractionCovered = (queryEnd-queryStart) / ((long double) queryLength);
		//float referenceFractionCovered = (targetEnd-targetStart) / ((long double) targetLength);

		//Alignment alignment = {contigIndex, readIndex, strand, readStart, readEnd, contigStart, contigEnd, chainScore, readLength, contigLength}; //, score
		
		AlignmentShort alignment = {contigIndex, contigStart, contigEnd, chainScore}; //, score
		
		//if(readIndex == 86917){
		//	cout << contigIndex << " " << identity << " " << alignment.score() << endl;
		//	getchar();
		//}
		//if(identity < 0.95) return;
		//if(!alignment.isMaximalMapping(300)) return;


		//indexReadAlignment(readIndex, alignment);
		/*
		if(_allAlignments2.find(readIndex) == _allAlignments2.end()){
			_allAlignments2[readIndex].push_back(alignment);
		}
		else{
			//if(chainScore > _allAlignments2[readIndex]._identity){
			//	_allAlignments2[readIndex] = alignment;
			//}
		}*/

		
		_allAlignments2[readIndex].push_back(alignment);
	}


	
	void indexReadAlignment(const ReadType& readIndex, const Alignment& alignment){
		

		if(_allAlignments.find(readIndex) == _allAlignments.end()){
			_allAlignments[readIndex].push_back(alignment);
			return;
		}

		vector<Alignment>& existingAlignments = _allAlignments[readIndex];

		vector<Alignment>::iterator it = existingAlignments.begin();
		bool isBetterAlignment = false;
		bool hasOverlap = false;
		bool overlapWithBetterAlignment = false;
		//bool canAddAlignment = true;


		for(const Alignment& existingAlignment : existingAlignments){

			if(alignmentOverlapExistingAlignment(alignment, existingAlignment)){

				if(alignment.score() < existingAlignment.score()){
					overlapWithBetterAlignment = true;
				}

				hasOverlap = true;
			}

		}

		if(overlapWithBetterAlignment){
			return;
		}

		while(it != existingAlignments.end()) {

			const Alignment& existingAlignment = *it;



			if(alignmentOverlapExistingAlignment(alignment, existingAlignment) && alignment.score() > existingAlignment.score()) {
				it = existingAlignments.erase(it);
				isBetterAlignment = true;
			}
			else{
				++it;
			}
		}



		if(isBetterAlignment){
			_allAlignments[readIndex].push_back(alignment);
		}

		if(!isBetterAlignment && !hasOverlap){
			_allAlignments[readIndex].push_back(alignment);
		}
	}

	bool alignmentOverlapExistingAlignment(const Alignment& alignment, const Alignment& existingAlignment){

		float allowedOverlap = 0;

		if(alignment._readStart >= existingAlignment._readStart && alignment._readEnd <= existingAlignment._readEnd) return true; //alignment contained
		if(alignment._readStart <= existingAlignment._readStart && alignment._readEnd >= existingAlignment._readEnd) return true; //existing alignment contained
		//if(alignment._readStart >= existingAlignment._readStart && (alignment._readStart+allowedOverlap) <= existingAlignment._readEnd) return true;
		//if(alignment._readEnd >= (existingAlignment._readStart+allowedOverlap) && alignment._readEnd <= existingAlignment._readEnd) return true;

		if(alignment._readStart >= existingAlignment._readStart){
			if(existingAlignment._readEnd - alignment._readStart > allowedOverlap) return true;
		}

		if(alignment._readEnd <= existingAlignment._readEnd){
			if(alignment._readEnd - existingAlignment._readStart > allowedOverlap) return true;
		}
		//int64_t startPos = max(alignment._readStart, existingAlignment._readStart);
		//int64_t endPos = min(alignment._readEnd, existingAlignment._readEnd);

		//int64_t overlapLength = endPos - startPos;

		//if(overlapLength > allowedOverlap) return true;

		return false;
	}


	u_int32_t contigNameToContigIndex(string contigName){

		string name = Utils::shortenHeader(contigName);
		name.erase(0,3); //"remove "ctg" letters
		return stoull(name);

	}

	u_int32_t readNameToReadIndex(string readName){

		string name = Utils::shortenHeader(readName);
		name.erase(0,5); //"remove "read_" letters
		return stoull(name);

	}

	void trimContigs(){

		_outputContigFile = gzopen(_outputContigFilename.c_str(),"wb");

		//cout << "single core here" << endl;
		ReadParserParallel readParser(_inputContigFilename, true, false, _nbCores);
		readParser.parse(ContigTrimmerFunctor(*this));

		gzclose(_outputContigFile);

	}



	u_int64_t _nbRepeatitiveContigs;

	class ContigTrimmerFunctor {

		public:

		ContigTrimmer& _parent;
		mm_tbuf_t* _tbuf;


		ContigTrimmerFunctor(ContigTrimmer& parent) : _parent(parent){
			_tbuf = mm_tbuf_init();
		}

		ContigTrimmerFunctor(const ContigTrimmerFunctor& copy) : _parent(copy._parent){
			_tbuf = mm_tbuf_init();
		}

		~ContigTrimmerFunctor(){
			mm_tbuf_destroy(_tbuf);
		}

		void operator () (const Read& read) {

			//if(read._seq.size() < 1000000) return;

			string header = read._header;

			Utils::ContigHeader contigHeader = Utils::extractContigHeader(header);

			u_int32_t contigIndex = contigHeader._contigIndex;
			if(read._index % 10000 == 0) Logger::get().debug() << "\t" << read._index;
			
			//cout << contigIndex << " " <<  << endl;
			
			const string& contigSequence = read._seq;
			vector<bool> isContigCovered(contigSequence.size(), false);

			if(_parent._alignments_readsVsContigs.find(contigIndex) == _parent._alignments_readsVsContigs.end()) return;

			auto& alignments = _parent._alignments_readsVsContigs[contigIndex];

			for(const auto& it : alignments){

				const u_int32_t contigStart = it.first;
				const u_int32_t contigEnd = it.second;
				//const ReadType& readIndex = it.first;
				//const Alignment& alignment = it.second;

				for(size_t i=contigStart; i < contigEnd; i++){
					isContigCovered[i] = true;
				}
			}

			u_int32_t startLengthToRemove = 0;
			for(int64_t i=0; i<isContigCovered.size(); i++){
				startLengthToRemove = i;
				if(isContigCovered[i]) break;
			}

			u_int32_t endLengthToRemove = 0;
			for(int64_t i = isContigCovered.size()-1; i >= 0; i--){
				endLengthToRemove = (isContigCovered.size() - i - 1);
				if(isContigCovered[i]) break;
			}

			if(startLengthToRemove < 50) startLengthToRemove = 0;
			if(endLengthToRemove < 50) endLengthToRemove = 0;
			//for(int64_t i=0; i<isContigCovered.size(); i++){
			//	cout << isContigCovered[i] << " ";
			//}
			//cout << endl;

			//cout << startLengthToRemove << " " << endLengthToRemove << endl;

			if(startLengthToRemove + endLengthToRemove >= contigSequence.size()) return;

			string contigSequenceTrimmed = contigSequence;

			if(startLengthToRemove > 0){
				contigSequenceTrimmed = contigSequenceTrimmed.substr(startLengthToRemove, contigSequenceTrimmed.size()-startLengthToRemove);
			}

			if(endLengthToRemove > 0){
				contigSequenceTrimmed = contigSequenceTrimmed.substr(0, contigSequenceTrimmed.size()-endLengthToRemove);
			}

			if(contigSequenceTrimmed.size() < _parent._minContigLength) return;



			//cout << nbIters << endl;

			//if(mostAbundantRepeat > 20 ){
			//	cout << header << endl;
			//	cout << contigSequenceTrimmed << endl;
				//getchar();
			//	return;
			//}

			//cout << contigSequenceTrimmed << endl;
			//cout << header << endl;
			//cout << computeSequenceComplexity(read._seq, 64, 32) << endl;

			if(contigHeader._isCircular){
				u_int32_t selfOverlapLength = computeSelfOverlap(contigSequenceTrimmed, _tbuf, _parent._minimap2Preset_map);
				//cout << "Self overlap: " << selfOverlapLength << " " << read._seq.size() << endl;

				//if(!contigHeader._isCircular && selfOverlapLength > 0){
				//cout << "Self overlap: " << selfOverlapLength << " " << read._seq.size() << endl;
				//}

				if(selfOverlapLength > 0){
					contigSequenceTrimmed = contigSequenceTrimmed.substr(0, contigSequenceTrimmed.size()-selfOverlapLength);
				}

				if(contigSequenceTrimmed.size() < _parent._minContigLength) return;
					
				
			}


			


			header = Utils::createContigHeader(contigHeader._contigIndex, contigSequenceTrimmed.size(), contigHeader._coverage, contigHeader._isCircular);
			
			#pragma omp critical(writeContigs)
			{
				string headerFasta = ">" + header + '\n';
				gzwrite(_parent._outputContigFile, (const char*)&headerFasta[0], headerFasta.size());
				contigSequenceTrimmed +=  '\n';
				gzwrite(_parent._outputContigFile, (const char*)&contigSequenceTrimmed[0], contigSequenceTrimmed.size());
			}
			
		}


		u_int32_t computeSelfOverlap(const string& sequence, mm_tbuf_t* tbuf, const string preset){

			string fakeName = "target";
			//mm_tbuf_t* tbuf = mm_tbuf_init();

			//allAlignments.clear();
			//cout << "Target size: " << reference.size() << endl;
			//cout << "Query size: " << query.size() << endl;
			
			//AlignmentBounds alignment;
			//alignment._referenceLength = reference.size();
			//alignment._queryLength = query.size();

			//string preset = _minimap2Preset_ava
			mm_idxopt_t iopt;
			mm_mapopt_t mopt;
			mm_set_opt(0, &iopt, &mopt); //"ava-ont"
			mm_set_opt(preset.c_str(), &iopt, &mopt); //"ava-ont"
			iopt.batch_size = 0x7fffffffffffffffL; //always build a uni-part index

			//if(performBaseLevelAlignment){
			mopt.flag |= MM_F_CIGAR; // perform alignment 
			//}
			

			mopt.flag |= MM_F_ALL_CHAINS | MM_F_NO_DIAG | MM_F_NO_DUAL | MM_F_NO_LJOIN; // -D -P --no-long-join --dual=no

			mm_idx_t* mi = minimap2index(iopt.w, iopt.k, iopt.flag&1, iopt.bucket_bits, sequence);
			//mm_idx_t* mi = minimap2_indexRead(sequence);
			//mopt.min_chain_score = minChainScore; //-m 500
			//iopt.bucket_bits = 14;

			//cout << iopt.w << endl;
			//cout << iopt.k << endl;
			//cout << iopt.bucket_bits << endl;
			//cout << (iopt.flag&1) << endl;
			//cout << mopt.min_chain_score << endl;

			//mm_idxopt_t idx_opt;
			//mm_idx_t* mi = minimap2index(iopt.w, iopt.k, iopt.flag&1, iopt.bucket_bits, reference);
			//mm_idx_t* mi = minimap2index(iopt.w, iopt.k, iopt.flag&1, iopt.bucket_bits, reference);
			mm_mapopt_update(&mopt, mi);
			//mopt.mid_occ = 10;
			//mopt.mid_occ = 1000; // don't filter high-occ seeds

			mm_reg1_t *reg;
			//mm_tbuf_t *tbuf = mm_tbuf_init();
			int j, i, n_reg;
			reg = mm_map(mi, sequence.size(), sequence.c_str(), &n_reg, tbuf, &mopt, fakeName.c_str()); // get all hits for the query

			//int32_t expectedOverlapLength = getEstimatedOverlapLength(alignment1, alignment2);
			//int32_t expectedQueryStart = alignment2._readStart;
			//int32_t expectedQueryEnd = expectedQueryStart + expectedOverlapLength;

			u_int32_t maxSelfOverlapLength = 0;

			//if(n_reg > 1) cout << "----" << endl;
			for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
			
				mm_reg1_t *r = &reg[j];


				//cout << string(mi->seq[r->rid].name) << " " << r->qs << "\t" << r->qe << "\t" << r->rs << "\t" << r->re << "\t" << r->blen << endl;

				free(r->p);
				
				if(r->rev) continue; 
				if(r->qs > 50) continue; //Looking for alignment at the contig start
				if(sequence.size() - r->re > 50) continue; //Looking for alignment at the contig end

				//cout << r->qs << "\t" << r->qe << "\t" << r->rs << "\t" << r->re << endl;



				//cout << "Found self overlap:\t" << r->qs << "\t" << r->qe << "\t" << r->rs << "\t" << r->re << endl;

				u_int32_t overlapStart = r->qe;
				u_int32_t overlapEnd = sequence.size() - r->rs;

				u_int32_t selfOverlapLength = max(overlapStart, overlapEnd);
				if(selfOverlapLength >= sequence.size()) continue;

				if(selfOverlapLength > maxSelfOverlapLength){
					maxSelfOverlapLength = selfOverlapLength;
				}

			}
			
		
			free(reg);
			mm_idx_destroy(mi);
			

			return maxSelfOverlapLength;
		}

		mm_idx_t* minimap2index(int w, int k, int is_hpc, int bucket_bits, const string& sequence){
			const char *seq = sequence.c_str();
			int len = sequence.size();
			const char *fake_name = "target";
			char *s;
			mm_idx_t *mi;
			s = (char*)calloc(len + 1, 1);
			memcpy(s, seq, len);
			mi = mm_idx_str(w, k, is_hpc, bucket_bits, 1, (const char**)&s, (const char**)&fake_name);
			free(s);
			return mi;
		}


		
	};

};	

#endif 
