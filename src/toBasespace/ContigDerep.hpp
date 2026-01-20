

#ifndef MDBG_METAG_CONTIGDEREP
#define MDBG_METAG_CONTIGDEREP

#include "../Commons.hpp"
#include "../graph/GfaParser.hpp"



class ContigDerep{
    
public:

	string _inputContigFilename;
	string _outputContigFilename;
	string _tmpDir;
	string _tmpDirPolishing;
	float _minIdentity;
	int _nbCores;
	int _minimapBatchSize;
	string _alignFilename;
	u_int32_t _minContigLength;
	unordered_map<u_int32_t, float> _contigCoverages;

	typedef pair<int64_t, int64_t> Overlap;
	gzFile _outputContigFile;
	unordered_map<string, Overlap> _contigName_to_alignmentBounds;

	struct ContigOverlap{
		u_int32_t _contigIndex;
		u_int32_t _contigStart;
		u_int32_t _contigEnd;
	};

	//unordered_map<u_int32_t, float> _contigCoverages;
	unordered_map<u_int32_t, vector<ContigOverlap>> _contigOverlaps;

	ContigDerep(const string& inputContigFilename, const string& outputContigFilename, const int& minimapBatchSize, const string& tmpDir, const string& tmpDirPolishing, const float& minIdentity, const u_int32_t& minContigLength, const int& nbCores){
		_inputContigFilename = inputContigFilename;
		_outputContigFilename = outputContigFilename;
		_minimapBatchSize = minimapBatchSize;
		_tmpDir = tmpDir;
		_tmpDirPolishing = tmpDirPolishing;
		_minIdentity = minIdentity;
		_minContigLength = minContigLength;
		_nbCores = nbCores;

		_alignFilename = _tmpDirPolishing + "/align_contigsVsContigs.paf.gz";
		_fields = new vector<string>();
		_fields_optional = new vector<string>();

	}


    void execute (){

		Logger::get().debug() << "\tAligning contigs vs contigs";
		alignContigs();

		Logger::get().debug() << "\tLoading contig coverages";
		loadContigCoverages();

		Logger::get().debug() << "\tLoading contig overlaps";
		loadAlignments();
		
		Logger::get().debug() << "\tDereplicating contigs";
		derepContigs();
		
	}

	void alignContigs(){

		float minimapBatchSizeC = max(0.2, _minimapBatchSize / 4.0);
		//    let args = ["-t", &opts.nb_threads.to_string(), "-c", "-xasm20", "-DP", "--dual=no", "--no-long-join", "-r100", "-z200", "-g2k", fasta_path_str, fasta_path_str];
		string command = "minimap2 -v 0 -c -m 500 -x asm20 -I " + to_string(minimapBatchSizeC) + "G -t " + to_string(_nbCores) + " -DP --dual=no --no-long-join -r100 -z200 -g2k " + _inputContigFilename + " " + _inputContigFilename;
		//string command = "minimap2 -v 0 -X -m 500 -x asm20 -I " + to_string(minimapBatchSizeC) + "G -t " + to_string(_nbCores) + " " + _inputContigFilename + " " + _inputContigFilename;
		Utils::executeMinimap2(command, _alignFilename);
		//cout << command << endl;
		//command += " | gzip -c - > " + _alignFilename;
		//Utils::executeCommand(command, _tmpDir);

	}

	
	void loadAlignments(){
		Logger::get().debug() << "\tLoading alignments";

		PafParser pafParser(_alignFilename);
		auto fp = std::bind(&ContigDerep::loadAlignments_read, this, std::placeholders::_1);
		pafParser.parse(fp);
		
	}

	void loadContigCoverages(){
		

		auto fp = std::bind(&ContigDerep::loadContigCoverages_read, this, std::placeholders::_1);
		ReadParser readParser(_inputContigFilename, true, false);
		readParser.parse(fp);
		
	}

	void loadContigCoverages_read(const Read& read){

		Utils::ContigHeader contigHeader = Utils::extractContigHeader(read._header);

		//const u_int32_t contigIndex = contigNameToContigIndex(read._header);
		//const string& contigName = Utils::shortenHeader(read._header);
		
		//Utils::ContigHeader contigHeader = Utils::extractContigHeader(read._header);
		_contigCoverages[contigHeader._contigIndex] = contigHeader._coverage;

		//cout << contigName << " " << _contigCoverages[contigName] << endl;
	}

	enum MappingType {
		None,
		DovetailSuffix,
		QueryPrefix,
		QuerySuffix,
		ReferencePrefix,
		ReferenceSuffix,
		Internal,
		QueryContained,
		ReferenceContained,
		DovetailPrefix
	};

	vector<string>* _fields;
	vector<string>* _fields_optional;
	//unordered_set<string> _isDuplication;

	void loadAlignments_read(const string& line){

		const vector<string>& _fields = Utils::split(line, '\t');


		u_int32_t queryContigIndex = contigNameToContigIndex(_fields[0]);
		u_int32_t targetContigIndex = contigNameToContigIndex(_fields[5]);

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

		if(queryContigIndex == targetContigIndex) return;

		u_int32_t queryLength = stoull((_fields)[1]);
		u_int32_t queryStart = stoull((_fields)[2]);
		u_int32_t queryEnd = stoull((_fields)[3]);
		u_int32_t targetLength = stoull((_fields)[6]);
		u_int32_t targetStart = stoull((_fields)[7]);
		u_int32_t targetEnd = stoull((_fields)[8]);

		u_int32_t nbMatches = stoull((_fields)[9]);
		u_int32_t alignLength = stoull((_fields)[10]);
		//u_int64_t queryLength = stoull((*_fields)[1]);

		bool isReversed = (_fields)[4] == "-";

		float identity =  ((float)nbMatches) / alignLength;
	

		if((float)identity < (float)_minIdentity) return;

		//if(queryLength < targetLength){
			//if(queryLength > 30000) return;
		//}
		//else{
			//if(targetLength > 30000) return;
		//}

		float queryCoverage = _contigCoverages[queryContigIndex];
		float targetCoverage = _contigCoverages[targetContigIndex];
		
		//float queryFractionCovered = (queryEnd-queryStart) / ((long double) queryLength);
		//float referenceFractionCovered = (targetEnd-targetStart) / ((long double) targetLength);


		if(targetLength > queryLength){
			if(queryLength > 60000) return;
			if(queryCoverage > targetCoverage/2.0) return;
			_contigOverlaps[queryContigIndex].push_back({targetContigIndex, queryStart, queryEnd});

		}
		else{
			if(targetLength > 60000) return;
			if(targetCoverage > queryCoverage/2.0) return;
			_contigOverlaps[targetContigIndex].push_back({queryContigIndex, targetStart, targetEnd});

		}
	}

	u_int32_t contigNameToContigIndex(string contigName){

		string name = Utils::shortenHeader(contigName);
		name.erase(0,3); //"remove "ctg" letters
		return stoull(name);

	}

	void setContained(const string& contigName, const int64_t& contigLength){
		Overlap& existingBounds = _contigName_to_alignmentBounds[contigName];
		existingBounds.first = contigLength;
		existingBounds.second = 0;
	}

	void setStartBounds(const string& contigName, const int64_t& pos){
		Overlap& existingBounds = _contigName_to_alignmentBounds[contigName];
		existingBounds.first = max(existingBounds.first, pos);
	}

	void setEndBounds(const string& contigName, const int64_t& pos){
		Overlap& existingBounds = _contigName_to_alignmentBounds[contigName];
		existingBounds.second = min(existingBounds.second, pos);
	}

	void printBounds(const string& contigName, const int64_t& contigLength){
		cout << "0    " << _contigName_to_alignmentBounds[contigName].first << " " << _contigName_to_alignmentBounds[contigName].second << "    " << contigLength<< endl;
	}

	void initOverlap(const string& contigName, const int64_t& contigLength){

		if(_contigName_to_alignmentBounds.find(contigName) == _contigName_to_alignmentBounds.end()){
			_contigName_to_alignmentBounds[contigName] = {0, contigLength};
		}

	}
	/*
	void addOverlap(const string& contigName, const int32_t& contigStart, const int32_t& contigEnd){

		if(_contigName_to_alignmentBounds.find(contigName) == _contigName_to_alignmentBounds.end()){
			_contigName_to_alignmentBounds[contigName] = {contigStart, contigEnd};
			return;
		}

		if(contigStart > _contigName_to_alignmentBounds[contigName].first){
			_contigName_to_alignmentBounds[contigName].first = contigStart;
		}

		if(contigStart > _contigName_to_alignmentBounds[contigName].first){
			_contigName_to_alignmentBounds[contigName].first = contigStart;
		}
	}
	*/

	MappingType getMappingType(const int64_t& b1, const int64_t& e1, const int64_t& l1, const int64_t& b2, const int64_t& e2, const int64_t& l2){
		
		static float overhang = 1000;
		static float r = 0.8;

		int64_t left_overhang = min(b1,b2);
		int64_t right_overhang = min(l1-e1,l2-e2);
		int64_t maplen = max(e1-b1,e2-b2);
		int64_t oh_threshold = min(overhang, ((maplen)*r));

		if( b2 <= b1 && b2 <= oh_threshold && right_overhang > oh_threshold) {
			return MappingType::ReferencePrefix;
		} 
		else if (b1 <= b2 && b1 <= oh_threshold && right_overhang > oh_threshold) {
			return MappingType::QueryPrefix;
		} 
		else if (left_overhang > oh_threshold && l2-e2 <= oh_threshold && l2-e2 <= l1-e1) {
			return MappingType::ReferenceSuffix;
		} 
		else if (left_overhang > oh_threshold && l1-e1 <= oh_threshold && l1-e1 <= l2-e2) {
			return MappingType::QuerySuffix;
		} 
		else if (left_overhang > oh_threshold || right_overhang > oh_threshold) {
			return MappingType::Internal;
		} 
		else if (b1 >= b2 && l1-e1 >= l2-e2) {
			return MappingType::ReferenceContained;
		} 
		else if (b2 >= b1 && l2-e2 >= l1-e1) {
			return MappingType::QueryContained;
		} 
		else if (b1 <= b2) {
			return MappingType::DovetailPrefix;
		} 
		else {
			return MappingType::DovetailSuffix;
		}

		return MappingType::None;
	}

	void derepContigs(){

		_outputContigFile = gzopen(_outputContigFilename.c_str(),"wb");


		auto fp = std::bind(&ContigDerep::derepContigs_read, this, std::placeholders::_1);
		ReadParser readParser(_inputContigFilename, true, false);
		readParser.parse(fp);

		gzclose(_outputContigFile);
	}

	void derepContigs_read(const Read& read){


		Utils::ContigHeader contigHeader = Utils::extractContigHeader(read._header);
		//u_int32_t queryContigIndex = contigNameToContigIndex(read._header);
		//const string& contigName = Utils::shortenHeader(read._header);

		const pair<int64_t, int64_t> overlap = getOverlaps(contigHeader._contigIndex, read._seq.size());

		if(overlap.first == 0 && overlap.second == read._seq.size()){
			writeContig(read._header, read._seq);
		}
		else{

			//const Overlap& bounds = _contigName_to_alignmentBounds[contigName];
			if(overlap.first > overlap.second){ //Contig is contained in another contig or completely filled with overlaps
				return;
			}

			int64_t contigLength = overlap.second - overlap.first;
			if(contigLength < _minContigLength) return;

			string contigSequence = read._seq.substr(overlap.first, contigLength);
			writeContig(read._header, contigSequence);
		}
	}


	pair<int64_t, int64_t> getOverlaps(const u_int32_t queryContigIndex, const int64_t contigLength){

		const int64_t maxHang = 300;

		pair<int64_t, int64_t> overlapsResult = {0, contigLength};
		if(_contigOverlaps.find(queryContigIndex) == _contigOverlaps.end()) return overlapsResult;

		unordered_map<u_int32_t, vector<ContigOverlap>> contigIndex_to_overlaps;

		for(const ContigOverlap& overlap : _contigOverlaps[queryContigIndex]){
			contigIndex_to_overlaps[overlap._contigIndex].push_back(overlap);
		}

		float maxCoverage = -1;

		for(const auto& it : contigIndex_to_overlaps){
			const u_int32_t referenceContigIndex = it.first;

			vector<bool> coverages(contigLength, false);

			for(const ContigOverlap& overlap : it.second){
				for(size_t i=overlap._contigStart; i<overlap._contigEnd; i++){
					coverages[i] = true;
				}
			}

			pair<int64_t, int64_t> overlaps = {0, contigLength};
			vector<CoverageRegion> coverageRegions = collectCoveredFragments(coverages);
			
			for(int64_t i=0; i<coverageRegions.size(); i++){

				const CoverageRegion& fragment = coverageRegions[i];
				if(!fragment._isCovered && fragment.getLength() > maxHang) break;

				overlaps.first += fragment.getLength();

				//cout << i << " " << overlaps.first << " " << fragment.getLength() << endl;
			}

			for(int64_t i=((int64_t)coverageRegions.size())-1; i >= 0; i--){

				const CoverageRegion& fragment = coverageRegions[i];
				if(!fragment._isCovered && fragment.getLength() > maxHang) break;

				overlaps.second -= fragment.getLength();
			}


			//cout << overlaps.first << " " << overlaps.second << endl;

			overlapsResult.first = max(overlapsResult.first, overlaps.first);
			overlapsResult.second = min(overlapsResult.second, overlaps.second);

			//cout << overlapsResult.first << " " << overlapsResult.second << endl;

			//if(contigLength < 9000){
			//	cout << "lul" << endl;
			//	getchar();
			//}

		}
	
		//if(contigLength < 9000){
		//	cout << "luuuuul" << endl;
		//	cout << overlapsResult.first << " " << overlapsResult.second << endl;
		//	getchar();
		//}

		return overlapsResult;
	}

	struct CoverageRegion{


		size_t _startIndex;
		size_t _endIndex;
		bool _isCovered;

		int64_t getLength() const{
			return _endIndex - _startIndex + 1;
		}
	};


	vector<CoverageRegion> collectCoveredFragments(const vector<bool>& coverages){
		
		//for(size_t i=0; i<coverages.size(); i++){
		//	cout << i << "\t" << coverages[i] << endl;
		//}

		vector<CoverageRegion> regions;

		bool lastChar = coverages[0];
		u_int64_t lastPos = 0;

		for(size_t i=0; i<coverages.size(); i++){
			
			bool c = coverages[i];

			if(c == lastChar) continue;

			regions.push_back({lastPos, i-1, coverages[lastPos]});
			lastPos = i;
			lastChar = c;
		}

		regions.push_back({lastPos, coverages.size()-1, coverages[lastPos]});
		//cout << lastChar << endl;
		//rleSequence += lastChar;
		//rlePositions.push_back(lastPos);
		//rlePositions.push_back(length);


		//for(int64_t i=0; i<regions.size(); i++){
		//	cout << i << "\t" << regions[i]._startIndex << "\t" << regions[i]._endIndex << "\t" << regions[i]._isCovered << endl;
		//}



		return regions;
	}

	void writeContig(string header, const string& sequence){

		Utils::ContigHeader contigHeader = Utils::extractContigHeader(header);
		header = Utils::createContigHeader(contigHeader._contigIndex, sequence.size(), contigHeader._coverage, contigHeader._isCircular);
		

		string headerFasta = ">" + header + '\n';
		gzwrite(_outputContigFile, (const char*)&headerFasta[0], headerFasta.size());

		string sequenceFasta = sequence + '\n';
		gzwrite(_outputContigFile, (const char*)&sequenceFasta[0], sequenceFasta.size());
	}

};	

#endif 
