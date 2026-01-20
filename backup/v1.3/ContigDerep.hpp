

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
	bool _useMetamdbgHeaderStyle;

	typedef pair<int64_t, int64_t> Overlap;
	gzFile _outputContigFile;
	//unordered_map<string, Overlap> _contigName_to_alignmentBounds;

	//unordered_set<string> _isContigContained;

	//unordered_map<string, float> _mapToReferenceCoverages;

	struct ContigOverlap{
		u_int32_t _contigIndex;
		u_int32_t _contigStart;
		u_int32_t _contigEnd;
	};

	unordered_map<u_int32_t, float> _contigCoverages;
	unordered_map<u_int32_t, vector<ContigOverlap>> _contigOverlaps;

	//struct ReferenceMapping{
	//	float _referenceCoverage;
	//	u_int32_t _lagestAlignment;
	//};

	//unordered_map<string, float> _referenceMappingCoverage;

	ContigDerep(const string& inputContigFilename, const string& outputContigFilename, bool useMetamdbgHeaderStyle, const int& minimapBatchSize, const string& tmpDir, const string& tmpDirPolishing, const float& minIdentity, const u_int32_t& minContigLength, const int& nbCores){
		_inputContigFilename = inputContigFilename;
		_outputContigFilename = outputContigFilename;
		_useMetamdbgHeaderStyle = useMetamdbgHeaderStyle;
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

		loadContigCoverages();

		Logger::get().debug() << "\tLoading contig overlaps";
		loadAlignments();
		
		Logger::get().debug() << "\tDereplicating contigs";
		derepContigs();
		
	}
	
	u_int32_t contigNameToContigIndex(string contigName){

		string name = Utils::shortenHeader(contigName);
		name.erase(0,3); //"remove "ctg" letters
		return stoull(name);

	}

	void loadContigCoverages(){
		

		auto fp = std::bind(&ContigDerep::loadContigCoverages_read, this, std::placeholders::_1);
		ReadParser readParser(_inputContigFilename, true, false);
		readParser.parse(fp);
		
	}

	void loadContigCoverages_read(const Read& read){

		const u_int32_t contigIndex = contigNameToContigIndex(read._header);
		//const string& contigName = Utils::shortenHeader(read._header);
		
		Utils::ContigHeader contigHeader = Utils::extractContigHeader(read._header);
		_contigCoverages[contigIndex] = contigHeader._coverage;

		//cout << contigName << " " << _contigCoverages[contigName] << endl;
	}
	

	void alignContigs(){

		int minimapBatchSizeC = max(1, _minimapBatchSize / 4);
		//    let args = ["-t", &opts.nb_threads.to_string(), "-c", "-xasm20", "-DP", "--dual=no", "--no-long-join", "-r100", "-z200", "-g2k", fasta_path_str, fasta_path_str];

		//string command = "minimap2 -v 0 -c -m 500 -x asm20 -I " + to_string(minimapBatchSizeC) + "G -t " + to_string(_nbCores) + " -DP --dual=no --no-long-join -r100 -z200 -g2k " + _inputContigFilename + " " + _inputContigFilename;
		string command = "minimap2 -v 0 -X -m 500 -x asm20 -I " + to_string(minimapBatchSizeC) + "G -t " + to_string(_nbCores) + " " + _inputContigFilename + " " + _inputContigFilename;
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
	

		//if((float)identity < (float)_minIdentity) return;

		if(queryLength < targetLength){
			//if(queryLength > 30000) return;
		}
		else{
			//if(targetLength > 30000) return;
		}


		
		//float queryFractionCovered = (queryEnd-queryStart) / ((long double) queryLength);
		//float referenceFractionCovered = (targetEnd-targetStart) / ((long double) targetLength);

		
		if(targetLength > queryLength){
			_contigOverlaps[queryContigIndex].push_back({targetContigIndex, queryStart, queryEnd});
			/*
			if(queryFractionCovered > 0.5){

				if(_referenceMappingCoverage.find(queryName) ==_referenceMappingCoverage.end()){
					_referenceMappingCoverage[queryName] = _contigCoverages[targetName];
				}
				else{
					if(_contigCoverages[targetName] > _referenceMappingCoverage[queryName]){
						_referenceMappingCoverage[queryName] = _contigCoverages[targetName];
					}
				}
				//_isContigContained.insert(queryName);
			}
			*/
		}
		else{
			_contigOverlaps[targetContigIndex].push_back({queryContigIndex, targetStart, targetEnd});
			/*
			if(referenceFractionCovered > 0.5){
				//_isContigContained.insert(targetName);

				if(_referenceMappingCoverage.find(targetName) == _referenceMappingCoverage.end()){
					_referenceMappingCoverage[targetName] = _contigCoverages[queryName];
				}
				else{
					if(_contigCoverages[queryName] > _referenceMappingCoverage[targetName]){
						_referenceMappingCoverage[targetName] = _contigCoverages[queryName];
					}
				}
			}
			*/
		}


	}


	void derepContigs(){

		_outputContigFile = gzopen(_outputContigFilename.c_str(),"wb");


		auto fp = std::bind(&ContigDerep::derepContigs_read, this, std::placeholders::_1);
		ReadParser readParser(_inputContigFilename, true, false);
		readParser.parse(fp);

		gzclose(_outputContigFile);
	}


	void derepContigs_read(const Read& read){

		u_int32_t queryContigIndex = contigNameToContigIndex(read._header);

		//const string& queryName = Utils::shortenHeader(read._header);
		//queryName.erase(0,3);
		//u_int32_t queryContigIndex = stoull(queryName);

		float referenceCoverage = getReferenceCoverage(queryContigIndex, read._seq.size());

		if(referenceCoverage == -1){
			writeContig(read._header, read._seq);
			return;
		}


		//Utils::ContigHeader contigHeader = Utils::extractContigHeader(read._header);

		float contigCoverage = _contigCoverages[queryContigIndex];

		//cout << contigName << " " << contigCoverage << " " << referenceCoverage << endl;
		if(contigCoverage < referenceCoverage/2.0){
			//cout << "Removed:\t" << read._header << "\t" << read._seq.size() << "\t" << contigCoverage << "\t" << referenceCoverage << "\t" << contigCoverage/referenceCoverage << endl;
			return;
		}

		//cout << contigCoverage << " " << referenceCoverage << endl;
		//if(_isContigContained.find(contigName) != _isContigContained.end()) return;
		
		writeContig(read._header, read._seq);

		/*
		const string& contigName = Utils::shortenHeader(read._header);
		
		if(_referenceMappingCoverage.find(contigName) == _referenceMappingCoverage.end()){
			writeContig(read._header, read._seq);
			return;
		}


		Utils::ContigHeader contigHeader = Utils::extractContigHeader(read._header);

		float contigCoverage = _contigCoverages[contigName];
		float referenceCoverage = _referenceMappingCoverage[contigName];

		//cout << contigName << " " << contigCoverage << " " << referenceCoverage << endl;
		if(contigCoverage < referenceCoverage/2.0) return;

		//cout << contigCoverage << " " << referenceCoverage << endl;
		//if(_isContigContained.find(contigName) != _isContigContained.end()) return;
		
		writeContig(read._header, read._seq);
		*/
		
	}

	float getReferenceCoverage(const u_int32_t queryContigIndex, const u_int32_t contigLength){

		if(_contigOverlaps.find(queryContigIndex) == _contigOverlaps.end()) return -1;

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

			long double nb = 0;

			for(size_t i=0; i<coverages.size(); i++){
				if(coverages[i]){
					nb += 1;
				}
			}

			float fractionCovered = nb / coverages.size();

			if(fractionCovered > 0.95){
				if(_contigCoverages[referenceContigIndex] > maxCoverage){
					maxCoverage = _contigCoverages[referenceContigIndex];
				}
			}
		}
	
		return maxCoverage;
	}

	void writeContig(string header, const string& sequence){

		//if(_useMetamdbgHeaderStyle){
		//	Utils::ContigHeader contigHeader = Utils::extractContigHeader(header);
		//	header = Utils::createContigHeader(contigHeader._contigIndex, sequence.size(), contigHeader._coverage, contigHeader._isCircular);
		//}

		string headerFasta = ">" + header + '\n';
		gzwrite(_outputContigFile, (const char*)&headerFasta[0], headerFasta.size());

		string sequenceFasta = sequence + '\n';
		gzwrite(_outputContigFile, (const char*)&sequenceFasta[0], sequenceFasta.size());
	}

};	

#endif 


