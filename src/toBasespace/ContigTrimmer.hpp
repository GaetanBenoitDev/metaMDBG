

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

		Logger::get().debug() << "\tAligning reads vs contigs";
		alignReads();

		//Logger::get().debug() << "\tLoading contig coverages";
		//loadContigCoverages();

		//loadAlignmentsTrue();

		Logger::get().debug() << "\tLoading alignments";
		loadAlignments();
		
		//exit(1);
		Logger::get().debug() << "\tTrimming contigs";
		trimContigs();
		
	}

	void alignReads(){

		float minimapBatchSizeC = max(1.0, _minimapBatchSize / 4.0);
		//    let args = ["-t", &opts.nb_threads.to_string(), "-c", "-xasm20", "-DP", "--dual=no", "--no-long-join", "-r100", "-z200", "-g2k", fasta_path_str, fasta_path_str];
		//string command = "minimap2 -v 0 -p 1 -c -I " + to_string(minimapBatchSizeC) + "G -K 0.02G -t " + to_string(_nbCores) + " -x " + _minimap2Preset_map + " " + _inputContigFilename + " " + _usedReadFilename;
		string command = "minimap2 -v 0 -p 1 -I " + to_string(minimapBatchSizeC) + "G -t " + to_string(_nbCores) + " -x " + _minimap2Preset_map + " " + _inputContigFilename + " " + _usedReadFilename;
		//string command = "minimap2 -v 0 -X -m 500 -x asm20 -I " + to_string(minimapBatchSizeC) + "G -t " + to_string(_nbCores) + " " + _inputContigFilename + " " + _inputContigFilename;
		Utils::executeMinimap2(command, _alignFilename);
		//cout << command << endl;
		//command += " | gzip -c - > " + _alignFilename;
		//Utils::executeCommand(command, _tmpDir);
		

	}




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














	
	void loadAlignments(){

		PafParser pafParser(_alignFilename);
		//PafParser pafParser(_readPartitionDir + "/true_align.paf.gz");
		auto fp = std::bind(&ContigTrimmer::loadAlignments_read, this, std::placeholders::_1);
		pafParser.parse(fp);
		
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


	void loadAlignments_read(const string& line){

		const vector<string>& _fields = Utils::split(line, '\t');


		//u_int32_t queryContigIndex = contigNameToContigIndex(_fields[0]);
		u_int32_t contigIndex = contigNameToContigIndex(_fields[5]);
		ReadType readIndex = readNameToReadIndex(_fields[0]);
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
	
		float chainScore = 0;
		
		for(size_t i=12; i<_fields.size(); i++){


			const vector<string>& fields2 = Utils::split(_fields[i], ':');

			//if((*_fields_optional)[0] == "gi"){
			if(fields2[0] == "s1"){
				chainScore = std::stoi(fields2[2]);
				break;
			}

		}

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

		ReadParserParallel readParser(_inputContigFilename, true, false, 1); //single core
		readParser.parse(ContigTrimmerFunctor(*this));

		gzclose(_outputContigFile);

	}

	class ContigTrimmerFunctor {

		public:

		ContigTrimmer& _parent;


		ContigTrimmerFunctor(ContigTrimmer& parent) : _parent(parent){
		}

		ContigTrimmerFunctor(const ContigTrimmerFunctor& copy) : _parent(copy._parent){
		}

		~ContigTrimmerFunctor(){
		}

		void operator () (const Read& read) {

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
			
			header = Utils::createContigHeader(contigHeader._contigIndex, contigSequenceTrimmed.size(), contigHeader._coverage, contigHeader._isCircular);
		

			string headerFasta = ">" + header + '\n';
			gzwrite(_parent._outputContigFile, (const char*)&headerFasta[0], headerFasta.size());
			contigSequenceTrimmed +=  '\n';
			gzwrite(_parent._outputContigFile, (const char*)&contigSequenceTrimmed[0], contigSequenceTrimmed.size());

		}
	};

};	

#endif 
