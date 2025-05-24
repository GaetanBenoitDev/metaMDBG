

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
	string _maxMemory;
	string _alignFilename;

	typedef pair<int64_t, int64_t> Overlap;
	gzFile _outputContigFile;
	unordered_map<string, Overlap> _contigName_to_alignmentBounds;

	ContigDerep(const string& inputContigFilename, const string& outputContigFilename, const string& maxMemory, const string& tmpDir, const string& tmpDirPolishing, const float& minIdentity, const int& nbCores){
		_inputContigFilename = inputContigFilename;
		_outputContigFilename = outputContigFilename;
		_maxMemory = maxMemory;
		_tmpDir = tmpDir;
		_tmpDirPolishing = tmpDirPolishing;
		_minIdentity = minIdentity;
		_nbCores = nbCores;

		_alignFilename = _tmpDirPolishing + "/align_contigsVsContigs.paf.gz";
		_fields = new vector<string>();
		_fields_optional = new vector<string>();

	}


    void execute (){

		cout << "Aligning contigs vs contigs" << endl;
		alignContigs();

		cout << "Loading contig overlaps" << endl;
		loadAlignments();
		
		cout << "Dereplicating contigs" << endl;
		derepContigs();
		
	}

	void alignContigs(){

		//    let args = ["-t", &opts.nb_threads.to_string(), "-c", "-xasm20", "-DP", "--dual=no", "--no-long-join", "-r100", "-z200", "-g2k", fasta_path_str, fasta_path_str];

		string command = "minimap2 -c -m 500 -x asm20 -I " + _maxMemory + "G -t " + to_string(_nbCores) + " -DP --dual=no --no-long-join -r100 -z200 -g2k " + _inputContigFilename + " " + _inputContigFilename;
		Utils::executeMinimap2(command, _alignFilename);
		//cout << command << endl;
		//command += " | gzip -c - > " + _alignFilename;
		//Utils::executeCommand(command, _tmpDir);

	}

	
	void loadAlignments(){
		Logger::get().debug() << "Loading alignments";

		PafParser pafParser(_alignFilename);
		auto fp = std::bind(&ContigDerep::loadAlignments_read, this, std::placeholders::_1);
		pafParser.parse(fp);
		
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


		const string& queryName = Utils::shortenHeader((_fields)[0]);
		const string& targetName = Utils::shortenHeader((_fields)[5]);
		
		if(queryName == targetName) return;

		int64_t queryLength = stoull((_fields)[1]);
		int64_t queryStart = stoull((_fields)[2]);
		int64_t queryEnd = stoull((_fields)[3]);
		int64_t targetLength = stoull((_fields)[6]);
		int64_t targetStart = stoull((_fields)[7]);
		int64_t targetEnd = stoull((_fields)[8]);

		int64_t nbMatches = stoull((_fields)[9]);
		int64_t alignLength = stoull((_fields)[10]);
		//u_int64_t queryLength = stoull((*_fields)[1]);

		bool isReversed = (_fields)[4] == "-";

		float identity =  ((float)nbMatches) / alignLength;
	

		if((float)identity < (float)_minIdentity) return;

		MappingType mappingType;

		if(isReversed){
			mappingType = getMappingType(queryLength-queryEnd, queryLength-queryStart, queryLength, targetStart, targetEnd, targetLength);
		}
		else{
			mappingType = getMappingType(queryStart, queryEnd, queryLength, targetStart, targetEnd, targetLength);
		}

		if(mappingType == MappingType::QueryContained){

			if(queryLength > targetLength) return;
			
			initOverlap(queryName, queryLength);
			setContained(queryName, queryLength);

			//cout << "Query contained: " << queryName << endl;
			//printBounds(queryName, queryLength);
			//getchar();

		}
		else if(mappingType == MappingType::QueryPrefix){

			if(queryLength > targetLength) return;
			initOverlap(queryName, queryLength);

			if(isReversed){
				setEndBounds(queryName, queryStart);
			}
			else{
				setStartBounds(queryName, queryEnd);
			}

			//cout << "Prefix overlap: " << queryName << " " << overlap.first << " " << overlap.second << endl;
			//printBounds(queryName, queryLength);
			//getchar();
			
		}
		else if(mappingType == MappingType::QuerySuffix){
			
			if(queryLength > targetLength) return;
			initOverlap(queryName, queryLength);

			if(isReversed){
				setStartBounds(queryName, queryEnd);
			}
			else{
				setEndBounds(queryName, queryStart);
			}

			//cout << "Suffix overlap: " << queryName << " " << overlap.first << " " << overlap.second << endl;
			//printBounds(queryName, queryLength);
			//getchar();
			
		}
		else if(mappingType == MappingType::ReferenceContained){
			
			if(targetLength > queryLength) return;

			initOverlap(targetName, targetLength);
			setContained(targetName, targetLength);

			//cout << "Reference contained: " << targetName << endl;
			//printBounds(targetName, targetLength);
			//getchar();
			
		}
		else if(mappingType == MappingType::ReferencePrefix){
			
			if(targetLength > queryLength) return;

			initOverlap(targetName, targetLength);
			setStartBounds(targetName, targetEnd);

			//cout << "Reference prefix: " << targetName << " " << 0 << " " << targetEnd << endl;
			//printBounds(targetName, targetLength);
			//getchar();
			
		}
		else if(mappingType == MappingType::ReferenceSuffix){
			
			if(targetLength > queryLength) return;

			initOverlap(targetName, targetLength);
			setEndBounds(targetName, targetStart);

			//cout << "Reference suffix: " << targetName << " " << targetStart << " " << targetLength << endl;
			//printBounds(targetName, targetLength);
			//getchar();
			
		}
		else if(mappingType == MappingType::DovetailPrefix){
			

			if(queryLength < targetLength){

				initOverlap(queryName, queryLength);

				if(isReversed){
					setEndBounds(queryName, queryStart);
				}
				else{
					setStartBounds(queryName, queryEnd);
				}
			}
			else{
				initOverlap(targetName, targetLength);
				setEndBounds(targetName, targetStart);
			}

			//cout << "Dovetail prefix: " << endl;
			//printBounds(queryName, queryLength);
			//printBounds(targetName, targetLength);
			//getchar();
			
			
		}
		else if(mappingType == MappingType::DovetailSuffix){
			

			if(queryLength < targetLength){
				initOverlap(queryName, queryLength);
				if(isReversed){
					setStartBounds(queryName, queryEnd);
				}
				else{
					setEndBounds(queryName, queryStart);
				}
			}
			else{
				initOverlap(targetName, targetLength);
				setStartBounds(targetName, targetEnd);
			}

			//cout << "Dovetail suffix: " << endl;
			//printBounds(queryName, queryLength);
			//printBounds(targetName, targetLength);
			//getchar();
			
		}


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
		
		Logger::get().debug() << "Dereplicating contigs" ;

		_outputContigFile = gzopen(_outputContigFilename.c_str(),"wb");


		auto fp = std::bind(&ContigDerep::derepContigs_read, this, std::placeholders::_1);
		ReadParser readParser(_inputContigFilename, true, false);
		readParser.parse(fp);

		gzclose(_outputContigFile);
	}

	void derepContigs_read(const Read& read){

		const string& contigName = Utils::shortenHeader(read._header);

		if(_contigName_to_alignmentBounds.find(contigName) == _contigName_to_alignmentBounds.end()){
			writeContig(read._header, read._seq);
		}
		else{

			const Overlap& bounds = _contigName_to_alignmentBounds[contigName];
			if(bounds.first > bounds.second){ //Contig is contained in another contig or completely filled with overlaps
				return;
			}

			int64_t contigLength = bounds.second - bounds.first;
			if(contigLength < 1000) return;

			string contigSequence = read._seq.substr(bounds.first, contigLength);
			writeContig(read._header, contigSequence);
		}
	}

	void writeContig(string header, const string& sequence){

		//Utils::ContigHeader contigHeader = Utils::extractContigHeader(header);
		//header = Utils::createContigHeader(contigHeader._contigIndex, sequence.size(), contigHeader._coverage, contigHeader._isCircular);

		string headerFasta = ">" + header + '\n';
		gzwrite(_outputContigFile, (const char*)&headerFasta[0], headerFasta.size());

		string sequenceFasta = sequence + '\n';
		gzwrite(_outputContigFile, (const char*)&sequenceFasta[0], sequenceFasta.size());
	}

};	

#endif 


