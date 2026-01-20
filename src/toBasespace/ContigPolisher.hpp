

#ifndef MDBG_METAG_CONTIGPOLISHER
#define MDBG_METAG_CONTIGPOLISHER

#include "../Commons.hpp"
#include "../utils/edlib.h"
#include "../utils/spoa/include/spoa/spoa.hpp"
#include "ContigDerep.hpp"

class ContigPolisher{
    
public:


	constexpr static const int64_t _minContigLength = 50;
	constexpr static const float _minContigCoverage = 1.5;
	constexpr static const int64_t _maximalMappingOffset = 300;


	string _inputFilenameRead;
	string _outputContigFilename;
	int _nbCores;
	size_t _windowLength;
	size_t _windowPositionOffset;
	size_t _windowLengthVariance;
	size_t _maxWindowCopies;
	//string _mapperOutputExeFilename;
	bool _useQual;
	string _outputDir;
	string _tmpDir;
	string _readPartitionDir;
	//u_int64_t _minContigLength;
	
	string _outputFilename_mapping;
	string _outputMappingFilename_contigsVsUsedReads;
	int _minimapBatchSize;
	double _qualityThreshold;
	bool _useMetamdbgHeaderStyle;
	int _nbPartitions;

	ofstream _file_contigHeaders;
	

	unordered_map<u_int32_t, u_int32_t> _contigLength;
	unordered_map<u_int32_t, float> _contigCoverages;
	string _minimap2Preset_map;
	bool _useHpc;
	u_int64_t _checksumWrittenReads;
	u_int64_t _outputContigIndex;
	
	struct Window{
		DnaBitset2* _sequence;
		string _quality;
		u_int32_t _posStart;
		u_int32_t _posEnd;
		float _score;
	};

	ContigPolisher(const string tmpDir, const string readPartitionDir, const string outputContigFilename, const int sequencingTechnology, const int minimapBatchSize, const bool useQual, const int nbCores, const int nbPartitions, const string minimap2Preset_map){
		_tmpDir = tmpDir;
		_readPartitionDir = readPartitionDir;
		_minimapBatchSize = minimapBatchSize;
		_useQual = useQual;
		_nbCores = nbCores;
		_nbPartitions = nbPartitions;
		_minimap2Preset_map = minimap2Preset_map;

		_inputFilenameRead = _readPartitionDir + "/input.txt";
		_outputContigFilename = outputContigFilename; //_tmpDir + "/contigs.fasta.gz";
		_useMetamdbgHeaderStyle = true;
		_windowLength = 500;
		_windowLengthVariance = _windowLength*0.01;
		_maxWindowCopies = 100;
		_qualityThreshold = 10.0;

		/*
		if(sequencingTechnology == DataType::HiFi){
			cout << "Data type is HiFi" << endl;
			_minimap2Preset_map = "map-hifi";
			_minimap2Preset_ava = "ava-pb";
			_useHpc = true;
		}
		else if(sequencingTechnology == DataType::Nanopore){
			cout << "Data type is ONT" << endl;
			_minimap2Preset_map = "map-ont";
			_minimap2Preset_ava = "ava-ont";
			_useHpc = false;
		}
		*/
	}


    void execute (){



		
		Logger::get().debug() << "";
		Logger::get().debug() << "Polishing contigs (pass 1)";
		_outputContigIndex = 0;
		

		for(size_t i=0; i<_nbPartitions; i++){


			string inputContigFilename = _readPartitionDir + "/" + to_string(i) + "_contigs_uncorrected.fasta.gz";
			string outputContigFilename = _readPartitionDir + "/" + to_string(i) + "_contigs_polished.fasta.gz";
			gzFile outputContigFile = gzopen(outputContigFilename.c_str(),"wb");

			polishPartition(i, inputContigFilename, outputContigFile, 0);

			gzclose(outputContigFile);
		}
	
		//Logger::get().debug() << "";
		//Logger::get().debug() << "Polishing contigs (pass 2)";
		//const string& outputContigFilename_polished = _readPartitionDir + "/contigs_polished.fasta.gz";
		//pileupContigs(_nbPartitions);
		//cout << "done" << endl;

		
		Logger::get().debug() << "";
		Logger::get().debug() << "Polishing contigs (pass 2)";
		_outputContigIndex = 0;

		gzFile outputContigFile_polished = gzopen(_outputContigFilename.c_str(),"wb");
		_file_contigHeaders = ofstream(_tmpDir + "/contigHeaders.txt");

		for(size_t i=0; i<_nbPartitions; i++){


			//if(i != 30) continue;

			string inputContigFilename = _readPartitionDir + "/" + to_string(i) + "_contigs_polished.fasta.gz";

			polishPartition(i, inputContigFilename, outputContigFile_polished, 1);

		}

		_file_contigHeaders.close();
		gzclose(outputContigFile_polished);
		
		


		
		//Logger::get().debug() << "";
		//Logger::get().debug() << "Polished contigs: " << _outputContigFilename;
		Logger::get().debug() << "done";
		//closeLogFile();
	}




	void polishPartition(u_int32_t partition, const string& contigFilename, gzFile& outputContigFile, const int& pass){
		_currentPartition = partition;
		//if(_contigSequences.size() == 0) return;
		
		Logger::get().debug() << "";
		Logger::get().debug() << "\tPolishing partition: " << _currentPartition << "/" << _nbPartitions;
		//if(_debug_contigName != "") cout << "\tPolishing partition: " << _currentPartition << "/" << _nbPartitions << endl;
		auto start = high_resolution_clock::now();

		//if(_partitionNbReads[partition] == 0) return;

		clearPass();

		//if(_debug_contigName != ""){
		//	debugFindContigs(contigFilename);
		//	if(!_debug_foundContig) return;

		//	if(!fs::exists(contigFilename)) return;
		//}

		//string contigFilename = _readPartitionDir + "/part_" + to_string(partition) + "_contigs.gz";
		//string readFilename = _readPartitionDir + "/part_" + to_string(partition) + ".gz";


		string readFilename = _readPartitionDir + "/" + to_string(partition) + "_reads.gz";
		//string contigFilename = _readPartitionDir + "/" + to_string(partition) + "_contigsRepeats.gz";
		
		//string outputContigFilename = _readPartitionDir + "/" + to_string(partition) + "_contigs.gz";

		
		Logger::get().debug() << "\tMap reads to curated contigs";

		const string& alignFilename = _readPartitionDir + "/align.paf.gz";
		string command = "minimap2 -I 100G -v 0 -m 500 -t " + to_string(_nbCores) + " -x " + _minimap2Preset_map + " " + contigFilename + " " + readFilename;
		Utils::executeMinimap2(command, alignFilename);
		//command += " | gzip -c - > " + alignFilename;
		//cout << command << endl;
		//Utils::executeCommand(command, _tmpDir);

		indexContigName(contigFilename);
		loadContigs(contigFilename);
		//loadAlignments(alignFilename, true, true);
		//computeContigCoverages(contigFilename);
		
		loadAllAlignments(alignFilename, true);
		computeContigCoveragesAll(contigFilename);


		
		//cout << "Load contigs" << endl;
		//loadContigs(_inputFilename_contigs, false);

		Logger::get().debug() << "\tCollecting window sequences";

		ReadParserParallel readParser(readFilename, true, false, _nbCores);
		readParser.parse(CollectWindowSequencesFunctor(*this));

		//cout << "miuom" << endl;
		//getchar();
		
		Logger::get().debug() << "\tPerform correction";
		performCorrection(outputContigFile, pass);


		Logger::get().debug() << "\tPartition " << partition << " done " << " (" << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s)";


		clearPass();

	}


	void computeContigCoveragesAll(const string& contigFilename){
		

		_contigLength.clear();
		_contigCoverages.clear();

		Logger::get().debug() << "\tComputing contig coverages";

		unordered_map<u_int32_t, vector<pair<u_int32_t, u_int32_t>>> contigHits;
		auto fp = std::bind(&ContigPolisher::computeContigCoverages_setup_read, this, std::placeholders::_1);
		ReadParser readParser(contigFilename, true, false);
		readParser.parse(fp);


		for(const auto& it : _allAlignments){

			const vector<Alignment>& alignments = it.second;

			for(const Alignment& alignment : alignments){
				//if(_contigLength.find(alignment._contigIndex) == _contigLength.end()) continue;
				contigHits[alignment._contigIndex].push_back({alignment._contigStart, alignment._contigEnd});
			}

		}

		for(auto& it : contigHits){

			//if(_contigLength.find(it.first) == _contigLength.end()){
			//	cout << "ups: " << it.first << endl;
			//}

			//cout << _contigLength[it.first] << endl;
			vector<u_int32_t> coverages(_contigLength[it.first], 0);

			for(auto& interval : it.second){

				//cout << "\t" << interval.first << " " << interval.second << endl;
				for(size_t i=interval.first; i<interval.second; i++){
					coverages[i] += 1;
				}
			}
			if (coverages.size() < 80 *2){
				_contigCoverages[it.first] = 1;
			}
			else{
				u_int64_t sum = 0;
				for(long i=75; i<((long)coverages.size())-75; i++){
					sum += coverages[i];
				}

				float coverage = sum / ((double)coverages.size());  
				_contigCoverages[it.first] = coverage;
			}


		}

		_contigLength.clear();
	}

	void computeContigCoverages_setup_read(const Read& read){

		//cout << read._index << "\t" << read._header << "\t" << read._seq.size() << endl;
		_contigLength[read._index] = read._seq.size();
		//int nbCounts = read._seq.size() / _contigCoverageWindow;
		//cout << Utils::shortenHeader(read._header) << " " << nbCounts << endl;
		//_contigHitPos[read._index].resize(nbCounts, 0);
	}

	void trimContigs(const string& inputContigFilename, const string& outputContigFilename, const int& minimapBatchSize){
		
		//ContigCurator contigCurator(_readPartitionDir, _minimap2Preset_map, _minimap2Preset_ava, 0, _minimizerSize, ContigPolisher::_minContigLength, _useMetamdbgHeaderStyle, _nbCores);
		//contigCurator.executeContigTrimming(inputContigFilename, outputContigFilename, minimapBatchSize);
	}

	void pileupContigs(int nbPartitions){

		//auto start = high_resolution_clock::now();

		//ContigPileup contigPileup(_tmpDir, _readPartitionDir, _minimap2Preset_map, _minimap2Preset_ava, nbPartitions, _useMetamdbgHeaderStyle, _minimizerSize, _minContigLength, _maximalMappingOffset, _nbCores);
		//contigPileup.execute();
		
		//Logger::get().debug() << "\tDone " << " (" << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s)";

	}
	


	
	

	u_int32_t _currentPartition;
	
	/*
	void debugFindContigs(const string& contigFilename){

		_debug_foundContig = false;

		auto fp = std::bind(&ContigPolisher::debugFindContigs_read, this, std::placeholders::_1);
		ReadParser readParser(contigFilename, true, false);
		readParser.parse(fp);
		
	}
	
	void debugFindContigs_read(const Read& read){
		if (read._header.find(_debug_contigName) != std::string::npos) {
			_debug_foundContig = true;
		}
	}
	*/


	void loadContigs(const string& contigFilename){

		_contigSequences.clear();
		_contigHeaders.clear();
		_contigWindowSequences.clear();
		
		auto fp = std::bind(&ContigPolisher::loadContigs_read, this, std::placeholders::_1);
		ReadParser readParser(contigFilename, true, false);
		readParser.parse(fp);
		
	}
	
	void loadContigs_read(const Read& read){

		//if(_debug_contigName != ""){
		//	if (read._header.find(_debug_contigName) != std::string::npos) {
		//	}
		//	else{
		//		return;
		//	}
		//}

		u_int32_t contigIndex = read._index;
		//const string& contigName = Utils::shortenHeader(read._header);

		//if(_currentPartition == 0) _logFile << "Loading contig in partition 0: " << contigIndex << endl;
		_contigSequences[contigIndex] = read._seq;

		size_t nbWindows = ceil(((double)read._seq.size()) / (double)_windowLength);

		//cout << "loadContigs_read: Nb windows + 1 si offsrt polishing" << endl;
		vector<vector<Window>> windows(nbWindows);

		_contigWindowSequences[contigIndex] = windows;
		_contigHeaders[contigIndex] = read._header;

		//cout << "Loaded contig: " << contigIndex << " " << read._header << " " << read._seq.size() << endl;
	}


	void clearPass(){

		//clearOverlap();
		_allAlignments.clear();
		_contigSequences.clear();
		//_validContigIndexes.clear();
		_contigWindowSequences.clear();
		_contigHeaders.clear();
		//_alignments_isReversed.clear();
		//_usedReadIndexes.clear();
		_contigName_to_contigIndex.clear();
		//_alignments.clear();
	}


	void indexContigName(const string& contigFilename){
		 
		//cout << "Indexing contig names: " << contigFilename << endl;
		_contigName_to_contigIndex.clear();

		auto fp = std::bind(&ContigPolisher::indexContigName_read, this, std::placeholders::_1);
		ReadParser readParser(contigFilename, true, false);
		readParser.parse(fp);
	}

	void indexContigName_read(const Read& read){

		//cout << "Index contig: " << Utils::shortenHeader(read._header) << " " << read._index << endl;
		_contigName_to_contigIndex[Utils::shortenHeader(read._header)] = read._index;
	}

	phmap::parallel_flat_hash_map<ReadType, vector<Alignment>> _allAlignments;
	//phmap::parallel_flat_hash_map<ReadType, Alignment> _alignments;
	phmap::parallel_flat_hash_map<u_int32_t, string> _contigSequences;
	phmap::parallel_flat_hash_map<u_int32_t, vector<vector<Window>>> _contigWindowSequences;
	phmap::parallel_flat_hash_map<u_int32_t, string> _contigHeaders;
	//unordered_map<ContigRead, u_int32_t, ContigRead_hash> _alignmentCounts;
	//u_int64_t _correctedContigIndex;

	phmap::parallel_flat_hash_map<string, u_int32_t> _contigName_to_contigIndex;
	phmap::parallel_flat_hash_map<string, ReadType> _readName_to_readIndex;

	//bool _indexPerContig;
	bool _useErrorFilter;


	void loadAllAlignments(const string& alignFilename, bool useIndexedName){

		Logger::get().debug() << "\tLoading all alignments";

		_allAlignments.clear();

		PafParser pafParser(alignFilename);
		auto fp = std::bind(&ContigPolisher::loadAllAlignments_read, this, std::placeholders::_1);
		pafParser.parse(fp);

		/*
		cout << "Nb alignment: " << _allAlignments.size() << endl;

		for(const auto& it : _allAlignments){
			if(it.second.size() > 1){
				cout << it.first << " " << it.second[0]._readLength << endl;
				for(const Alignment& al : it.second){
					cout << "\t" <<  al._readStart << " " << al._readEnd << "    " << al._contigIndex << "\t" << al._contigStart << "\t" << al._contigEnd << endl;
				}
			}
		}

		cout << "S2" << endl;
		*/
		//getchar();
		
	}

	void loadAllAlignments_read(const string& line){

		vector<string> _fields = Utils::split(line, '\t');

		const string& readName = Utils::shortenHeader((_fields)[0]);
		const string& contigName = Utils::shortenHeader((_fields)[5]);
		
		u_int64_t readLength = stoull((_fields)[1]);
		u_int32_t readStart = stoull((_fields)[2]);
		u_int32_t readEnd = stoull((_fields)[3]);
		u_int32_t contigLength = stoull((_fields)[6]);
		u_int32_t contigStart = stoull((_fields)[7]);
		u_int32_t contigEnd = stoull((_fields)[8]);

		u_int64_t nbMatches = stoull((_fields)[9]);
		u_int64_t alignLength = stoull((_fields)[10]);

		bool strand = (_fields)[4] == "-";

		AlignmentBounds bounds;
		bounds._queryStart = readStart;
		bounds._queryEnd = readEnd;
		bounds._referenceStart = contigStart;
		bounds._referenceEnd = contigEnd;
		bounds._queryLength = readLength;
		bounds._referenceLength = contigLength;
		bounds._isReversed = strand;


		u_int32_t readSizeMappable = bounds.getMappableLength(); //A read size that do not consider alignment outside contig bounds (start and end of the contigs)
		
		float identity = ((long double) nbMatches) / ((long double) readSizeMappable);

		//u_int32_t contigIndex = _contigName_to_contigIndex[contigName];
		//ReadType readIndex = _readName_to_readIndex[readName];

		u_int32_t readIndex = stoull(readName);;

		if(strand) return;
		if(_contigName_to_contigIndex.find(contigName) == _contigName_to_contigIndex.end()) return;

		//int64_t overhang = readSizeMappable - (bounds._queryEnd - bounds._queryStart);

		//if(overhang > 100) return;
		//float identity2 = ((long double) bounds._queryEnd-bounds._queryStart) / ((long double) readSizeMappable);


		//if(identity2 < 0.99) return;

		

		u_int32_t contigIndex = _contigName_to_contigIndex[contigName];

		Alignment alignment = {contigIndex, readIndex, strand, readStart, readEnd, contigStart, contigEnd, identity, readLength, contigLength}; //, score

		if(!alignment.isMaximalMapping(_maximalMappingOffset)) return;

		indexReadAlignment(readIndex, alignment);

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

		float allowedOverlap = 500;

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



	class CollectWindowSequencesFunctor {

		public:

		ContigPolisher& _parent;


		CollectWindowSequencesFunctor(ContigPolisher& parent) : _parent(parent){
		}

		CollectWindowSequencesFunctor(const CollectWindowSequencesFunctor& copy) : _parent(copy._parent){
			
		}

		~CollectWindowSequencesFunctor(){
		}

		void operator () (const Read& read) {

			u_int32_t readIndex = stoull(read._header);

			if(_parent._allAlignments.find(readIndex) == _parent._allAlignments.end()) return;

			const vector<Alignment>& als = _parent._allAlignments[readIndex];

			for(const Alignment& al : als){
				//const Alignment& al = _alignments[readIndex];
				//for(const Alignment& al : _alignments[readIndex]){
				u_int32_t contigIndex = al._contigIndex;

				if(_parent._contigSequences.find(contigIndex) == _parent._contigSequences.end()) continue;

				//_logFile << read._seq.size() << " " << read._qual.size() << " " << _contigSequences[contigIndex].size() << " " << al._readStart << " " << al._readEnd << " " << al._contigStart << " " << al._contigEnd << endl;
				string readSeq = read._seq;
				string qualSeq = read._qual;
				string readSequence = readSeq.substr(al._readStart, al._readEnd-al._readStart);
				string contigSequence = _parent._contigSequences[contigIndex].substr(al._contigStart, al._contigEnd-al._contigStart);



				//if(al._strand){
				//	Utils::toReverseComplement(readSequence);
				//	Utils::toReverseComplement(readSeq);
				//	std::reverse(qualSeq.begin(), qualSeq.end());
				//}
				

				AlignmentBounds bounds;
				bounds._queryStart = al._readStart;
				bounds._queryEnd = al._readEnd;
				bounds._referenceStart = al._contigStart;
				bounds._referenceEnd = al._contigEnd;
				bounds._queryLength = readSeq.size();
				bounds._referenceLength = _parent._contigSequences[contigIndex].size();
				bounds._isReversed = al._strand;


				u_int32_t readSizeMappable = bounds.getMappableLength(); //A read size that do not consider alignment outside contig bounds (start and end of the contigs)
				
				


				//if(_contigPolisher.getMaxhang(bounds) > ContigPolisher::maxHang*2) return;

				//u_int32_t readSizeMappable = _contigPolisher.getMappableLength(bounds); //A read size that do not consider alignment outside contig bounds (start and end of the contigs)
				

				//int64_t overhang = readSeq.size() - readSizeMappable;
				//cout << overhang << endl;
				//if(overhang > 10) return;

				//cout << overhang << endl;


				//_logFile << readSequence << endl;
				//_logFile << contigSequence << endl;

				//_logFile << contigSequence.size() << " "<< readSequence.size() << endl;
				static EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0);


				EdlibAlignResult result = edlibAlign(readSequence.c_str(), readSequence.size(), contigSequence.c_str(), contigSequence.size(), config);


				char* cigar;

				if (result.status == EDLIB_STATUS_OK) {
					cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
				} else {
					Logger::get().error() << "Invalid edlib results";
					exit(1);
				}

				//_logFile << cigar << endl;

				edlibFreeAlignResult(result);
				
				int64_t nbMatches = Commons::getCigarNbMatches(cigar, readSequence, contigSequence);
				float identity = ((long double) nbMatches) / ((long double) readSizeMappable);
				if(identity < 0.9) continue;
				//if(al._contigEnd < 20000){
					//cout << "-----" << endl;
					//cout << al._contigStart << " " << al._contigEnd << "     " << al._readStart << " " << al._readEnd << " " << al._readLength << endl;
					//cout << _contigPolisher.getMaxhang(bounds) << endl;
					//u_int64_t nbMatches = getCigarNbMatches(al, cigar, readSequence, contigSequence);
					//getchar();

				//}

				find_breaking_points_from_cigar(_parent._windowLength, al, readSeq.size(), cigar, readSeq, qualSeq, identity);


				free(cigar);

			}
		}

		void find_breaking_points_from_cigar(uint32_t window_length, const Alignment& al, u_int64_t readLength, char* cigar_, const string& readSequence, const string& qualSequence, const float& identity)
		{
			
			vector<std::pair<uint32_t, uint32_t>> breaking_points_;
			u_int64_t t_begin_ = al._contigStart;
			u_int64_t t_end_ = al._contigEnd;
			u_int64_t q_begin_ = al._readStart;
			u_int64_t q_end_ = al._readEnd;
			bool strand_ = al._strand;
			u_int64_t q_length_ = readLength;

			// find breaking points from cigar
			std::vector<int32_t> window_ends;
			for (uint32_t i = 0; i < t_end_; i += window_length) {
				if (i > t_begin_) {
					window_ends.emplace_back(i - 1);
				}
			}
			window_ends.emplace_back(t_end_ - 1);

			uint32_t w = 0;
			bool found_first_match = false;
			std::pair<uint32_t, uint32_t> first_match = {0, 0}, last_match = {0, 0};

			int32_t q_ptr = (strand_ ? (q_length_ - q_end_) : q_begin_) - 1;
			int32_t t_ptr = t_begin_ - 1;

			for (uint32_t i = 0, j = 0; i < strlen(cigar_); ++i) {
				if (cigar_[i] == 'M' || cigar_[i] == '=' || cigar_[i] == 'X') {
					uint32_t k = 0, num_bases = atoi(&cigar_[j]);
					j = i + 1;
					while (k < num_bases) {
						++q_ptr;
						++t_ptr;

						if (!found_first_match) {
							found_first_match = true;
							first_match.first = t_ptr;
							first_match.second = q_ptr;
						}
						last_match.first = t_ptr + 1;
						last_match.second = q_ptr + 1;
						if (t_ptr == window_ends[w]) {
							if (found_first_match) {
								breaking_points_.emplace_back(first_match);
								breaking_points_.emplace_back(last_match);
							}
							found_first_match = false;
							++w;
						}


						++k;
					}
				} else if (cigar_[i] == 'I') {
					q_ptr += atoi(&cigar_[j]);
					j = i + 1;
				} else if (cigar_[i] == 'D' || cigar_[i] == 'N') {
					uint32_t k = 0, num_bases = atoi(&cigar_[j]);
					j = i + 1;
					while (k < num_bases) {
						++t_ptr;
						if (t_ptr == window_ends[w]) {
							if (found_first_match) {
								breaking_points_.emplace_back(first_match);
								breaking_points_.emplace_back(last_match);
							}
							found_first_match = false;
							++w;
						}
						++k;
					}
				} else if (cigar_[i] == 'S' || cigar_[i] == 'H' || cigar_[i] == 'P') {
					j = i + 1;
				}
			}

			//if(breaking_points_.size() > 0) breaking_points_.emplace_back(last_match);
				
			for (uint32_t j = 0; j < breaking_points_.size(); j += 2) {

				if(breaking_points_[j].second >= readSequence.size()) return;
				if(breaking_points_[j + 1].second >= readSequence.size()) return;

				if (breaking_points_[j + 1].second - breaking_points_[j].second < 0.02 * _parent._windowLength) {
					continue;
				}

				
				if (qualSequence.size() > 0) {

					//const auto& quality = overlaps[i]->strand() ? sequence->reverse_quality() : sequence->quality();
					double average_quality = 0;
					for (uint32_t k = breaking_points_[j].second; k < breaking_points_[j + 1].second; ++k) {
						average_quality += static_cast<uint32_t>(qualSequence[k]) - 33;
					}
					average_quality /= breaking_points_[j + 1].second - breaking_points_[j].second;

					if (average_quality < _parent._qualityThreshold) {
						continue;
					}
				}
				

				//uint64_t window_id = id_to_first_window_id[overlaps[i]->t_id()] +
				uint64_t window_id = breaking_points_[j].first / _parent._windowLength;
				uint32_t window_start = (breaking_points_[j].first / _parent._windowLength) * _parent._windowLength;

				//const char* data = overlaps[i]->strand() ?
				//	&(sequence->reverse_complement()[breaking_points[j].second]) :
				//	&(sequence->data()[breaking_points[j].second]);
				const char* data = &readSequence[breaking_points_[j].second];
				uint32_t data_length = breaking_points_[j + 1].second - breaking_points_[j].second;

				string sequence = string(data, data_length);

				/*
				if(al._contigIndex == 1 && window_id > 1 && window_id < 1000){
					if(window_id == 192 || window_id == 193){
						_logFile << readSequence << endl;
						_logFile << sequence << endl;
						if(breaking_points_[j].second > 1) _logFile << readSequence[breaking_points_[j].second-1] << " " << readSequence[breaking_points_[j].second] << endl;
						if(readSequence.size() > readSequence[breaking_points_[j].second+data_length-1]) _logFile << readSequence[breaking_points_[j].second+data_length-1] << " " << readSequence[breaking_points_[j].second+data_length] << endl;
						_logFile << window_id << endl;
						//getchar();
					}
					

				}
				
				long pos = breaking_points_[j].second+data_length;
				while(true){
					if(pos >= readSequence.size()) break;
					char prevChar = readSequence[pos-1];
					char c = readSequence[pos];
					if(prevChar != c) break;

					//sequence += c;
					breaking_points_[j + 1].second += 1;
					//data_length += 1;

					pos += 1;
				}

				pos = ((long)breaking_points_[j].second)-1;
				while(true){
					if(pos < 0) break;
					char prevChar = readSequence[pos+1];
					char c = readSequence[pos];
					if(prevChar != c) break;

					//sequence += c;
					breaking_points_[j].second -= 1;
					//data_length += 1;

					pos -= 1;
				}
				*/

				data = &readSequence[breaking_points_[j].second];
				data_length = breaking_points_[j + 1].second - breaking_points_[j].second;
				sequence = string(data, data_length);

				/*
				if(al._contigIndex == 1 && window_id > 1 && window_id < 1000){
					if(window_id == 192 || window_id == 193){
						_logFile << readSequence << endl;
						_logFile << sequence << endl;
						if(breaking_points_[j].second > 1){
							_logFile << readSequence[breaking_points_[j].second-1] << " " << readSequence[breaking_points_[j].second] << endl;
							if(readSequence[breaking_points_[j].second-1] == readSequence[breaking_points_[j].second]) getchar();
						}
						if(readSequence.size() > readSequence[breaking_points_[j].second+data_length-1]) _logFile << readSequence[breaking_points_[j].second+data_length-1] << " " << readSequence[breaking_points_[j].second+data_length] << endl;
						_logFile << window_id << endl;
						//getchar();
					}
					

				}
				*/

                u_int32_t posStart = breaking_points_[j].first - window_start;
                u_int32_t posEnd =  breaking_points_[j + 1].first - window_start - 1;

				string quality = "";
				if(_parent._useQual && qualSequence.size() > 0){
					quality = string(&qualSequence[breaking_points_[j].second], data_length);
				}

				//if(al._contigIndex == 0 && window_id == 161){
					
					//_logFile << al._contigIndex << " " << al._readIndex << endl;
					//_logFile << (breaking_points_[j + 1].second - breaking_points_[j].second) << " " << (0.02 * _windowLength) << endl;
					//_logFile << sequence << endl;
					//_logFile << quality << endl;
					//getchar();
				//}

				//if(sequence.size() < 490) continue;
				indexWindow(al, window_id, posStart, posEnd, sequence, quality, identity);

				//if(window_id == 24){
				//	cout << sequence << endl;
				//	getchar();
				//}

				//_logFile << window_id << " " << posStart << " " << posEnd << endl;
				//_logFile << sequence << endl;
				//getchar();
				/*
				const char* quality = overlaps[i]->strand() ?
					(sequence->reverse_quality().empty() ?
						nullptr : &(sequence->reverse_quality()[breaking_points[j].second]))
					:
					(sequence->quality().empty() ?
						nullptr : &(sequence->quality()[breaking_points[j].second]));
				uint32_t quality_length = quality == nullptr ? 0 : data_length;
				*/
			}

			/*
			for(size_t i=0; i<breaking_points_.size()-1; i++){
				u_int64_t contigWindowStart = breaking_points_[i].first;
				u_int64_t contigWindowEnd = breaking_points_[i+1].first;
				u_int64_t readWindowStart = breaking_points_[i].second;
				u_int64_t readWindowEnd = breaking_points_[i+1].second;

				//_logFile << contigWindowStart << " " << contigWindowEnd  << "      " << readWindowStart << " " << readWindowEnd << endl;

				//if(readWindowEnd-readWindowStart < _windowLength) continue; //window sides
				//string windowSequence = 

				string qualSeq = "";
				if(qualSequence.size() > 0){
					qualSeq = qualSequence.substr(readWindowStart, readWindowEnd-readWindowStart);
				}

				indexWindow(al, readWindowStart, readWindowEnd, contigWindowStart, contigWindowEnd, readSequence.substr(readWindowStart, readWindowEnd-readWindowStart), qualSeq);
			}
			*/
			

			//for(const auto& breakPoint : breaking_points_){
			//	_logFile << breakPoint.first << " " << breakPoint.second << endl;
			//}
		}

		//void indexWindow(const Alignment& al, size_t readWindowStart, size_t readWindowEnd, size_t contigWindowStart, size_t contigWindowEnd, const string& windowSequence, const string& windowQualities){
		void indexWindow(const Alignment& al, u_int64_t windowIndex, u_int32_t posStart, u_int32_t posEnd, const string& windowSequence, const string& windowQualities, const float& identity){

			#pragma omp critical(indexWindow)
			{
				
				//size_t contigWindowIndex = contigWindowStart / _windowLength;
				vector<Window>& windows = _parent._contigWindowSequences[al._contigIndex][windowIndex];
				
				/*
				if(contigWindowIndex == 0){
					//_logFile << contigWindowStart << endl;
					_logFile << windowSequence << endl;
					//_logFile << windows.size() << endl;
					_logFile << readWindowStart << " " << readWindowEnd << endl;
					getchar();
					//_logFile << windowQualities << endl;
				}
				*/
				
				bool interrupt = false;
				if(_parent._maxWindowCopies == 0 || windows.size() < (_parent._maxWindowCopies-1)){


					windows.push_back({new DnaBitset2(windowSequence), windowQualities, posStart, posEnd, identity}); //al.score()

					/*
					if(al._contigIndex == 1 && windowIndex == 2){
						_logFile << "1111111111111111" << endl;
						_logFile  << windowSequence << endl;
					}
					if(al._contigIndex == 1 && windowIndex == 3){
						_logFile << "2222222222222222" << endl;
						_logFile  << windowSequence << endl;
					}
					*/

					/*
					if(windowSequences.size() > 1){
						static EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0);

						char* dnaStrModel = windowSequences[0]->to_string();
						EdlibAlignResult result = edlibAlign(dnaStrModel, strlen(dnaStrModel), windowSequence.c_str(), windowSequence.size(), config);

						char* cigar;

						if (result.status == EDLIB_STATUS_OK) {
							cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
							_logFile << strlen(cigar) << endl;
							free(cigar);
						} else {
							_logFile << "Invalid edlib results" << endl;
							exit(1);
						}

						//_logFile << cigar << endl;

						edlibFreeAlignResult(result);
						free(dnaStrModel);

					}
					*/


					interrupt = true;
				}

				
				if(!interrupt){

					float score = identity; //al.score();

					u_int64_t incompleteWindowIndex = -1;
					//u_int64_t minWindowSize = -1;

					u_int64_t largerDistanceWindow = 0;

					for(size_t i=0; i<windows.size(); i++){

						const Window& window = windows[i];
						u_int64_t distance = abs(((long)window._sequence->m_len) - ((long)_parent._windowLength));

						if(distance > _parent._windowLengthVariance){
							if(distance > largerDistanceWindow){
								largerDistanceWindow = distance;
								incompleteWindowIndex = i;
							}
						}
					}


					//u_int64_t distance = abs(((long)windowSequence.size()) - ((long)_windowLength));

					if(incompleteWindowIndex != -1){
						Window& window = windows[incompleteWindowIndex];
						delete window._sequence;
						windows[incompleteWindowIndex] = {new DnaBitset2(windowSequence), windowQualities, posStart, posEnd, score};
					}
					else{
						//u_int64_t largestAligmenet = 0;
						//u_int64_t largestAligmenetIndex = 0;

						/*
						size_t largerWindowIndex = 0;
						u_int64_t largerDistanceWindow = 0;

						for(size_t i=0; i<windows.size(); i++){

							const Window& window = windows[i];
							u_int64_t distance = abs(((long)window._sequence->m_len) - ((long)_windowLength));

							if(distance > largerDistanceWindow){
								largerDistanceWindow = distance;
								largerWindowIndex = i;
							}
						}


						u_int64_t distance = abs(((long)windowSequence.size()) - ((long)_windowLength));

						if(distance < largerDistanceWindow){
							Window& window = windows[largerWindowIndex];
							delete window._sequence;
							windows[largerWindowIndex] = {new DnaBitset2(windowSequence), windowQualities, posStart, posEnd};
						}
						*/

						

        				static float maxVal = std::numeric_limits<float>::max();
						size_t largerWindowIndex = 0;
						//u_int64_t largerDistanceWindow = 0;
						float lowestScore = maxVal;
						
						for(size_t i=0; i<windows.size(); i++){

							const Window& window = windows[i];
							//u_int64_t distance = abs(((long)window._sequence->m_len) - ((long)_windowLength));
							//float score =

							if(window._score < lowestScore){
								lowestScore = window._score;
								largerWindowIndex = i;
							}
						}

						//u_int64_t distance = abs(((long)windowSequence.size()) - ((long)_windowLength));

						if(score > lowestScore){
							Window& window = windows[largerWindowIndex];
							delete window._sequence;
							windows[largerWindowIndex] = {new DnaBitset2(windowSequence), windowQualities, posStart, posEnd, score};
						}

					}



				} 




			}
			//_logFile << contigWindowIndex << endl;


			

			/*
			if(contigWindowIndex == _contigWindowSequences[al._contigIndex].size()-1){

				for(size_t i=0; i<windowSequences.size(); i++){
					_logFile << windowSequences[i]->m_len << " ";
				}
				_logFile << endl;
				//_logFile << contigWindowEnd-contigWindowStart << endl;
				//_logFile << windowSequence << endl;
			}
			*/
		}


	};

	struct CorrectedWindow{
		size_t _windowIndex;
		DnaBitset2* _correctedSequence;
		bool _success;
	};

	static bool CorrectedWindowComparator (const CorrectedWindow& p1, const CorrectedWindow& p2){
		return p1._windowIndex < p2._windowIndex;
	}


	unordered_map<u_int32_t, vector<CorrectedWindow>> _currentContigs;



	void performCorrection(gzFile& outputContigFile, const int& pass){
		
		u_int64_t checksum = 0;

		vector<u_int32_t> contigIndexes;
		for(auto& it : _contigWindowSequences){
			contigIndexes.push_back(it.first);
		}
		//vector<u_int32_t> contigIndexes;
		//for(auto& it : _contigWindowSequences){
		//	contigIndexes.push_back(it.first);
		//}


		size_t i = 0;
		size_t windowIndex = 0;
		//unordered_map<string, vector<vector<Window>>> _contigWindowSequences;

		//cout << "Perform correction single core" << endl;
		
		#pragma omp parallel num_threads(_nbCores)
		{

			bool isEOF = false;
			size_t contigIndexLocal;
			size_t windowIndexLocal;
			std::unique_ptr<spoa::AlignmentEngine> alignmentEngine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -4);
			

			while(true){


				#pragma omp critical(lala)
				{

					if(i >= _contigWindowSequences.size()){
						isEOF = true;
					}

					if(!isEOF){

						u_int32_t contigIndex = contigIndexes[i];
						//cout << "Contig: " << contigIndex << endl;
						//cout << windowIndex << " " << _contigWindowSequences[contigIndex].size() << endl;
						if(windowIndex >= _contigWindowSequences[contigIndex].size()){
							i += 1;
							windowIndex = 0;
							//cout << "inc contig index" << endl;
						}

						if(i >= _contigWindowSequences.size()){
							isEOF = true;
						}


						if(!isEOF){
							contigIndex = contigIndexes[i];
							contigIndexLocal = contigIndex;
							windowIndexLocal = windowIndex;
							windowIndex += 1;
						}
					}
					
				}

				if(isEOF) break;
				//functorSub(read);

				//#pragma omp critical(lala)
				//{
				//	cout << contigIndexLocal << ": " << windowIndexLocal << "/" << _contigWindowSequences[contigIndexLocal].size() << endl;
				//	getchar();
				//}

				//const string& contigName = contigIndexes[contigIndexLocal];

				//cout << _contigWindowSequences.size() << " " << contigIndexLocal << endl;
				//cout << "lala: " << contigIndexLocal << " " << endl;
				
				//if(_contigWindowSequences.find(contigIndexLocal) == _contigWindowSequences.end()){
				//	cout << "omg 1 " << endl;
				//}

				//if(windowIndexLocal >= _contigWindowSequences[contigIndexLocal].size()){
				//	cout << "omg 2" << endl;
				//	cout << _contigHeaders[contigIndexLocal] << endl;
				//	cout << _contigSequences[contigIndexLocal].size() << endl;
				//	cout << contigIndexLocal << " " << windowIndexLocal << " " << _contigWindowSequences[contigIndexLocal].size() << endl;
				//}

				if(windowIndexLocal >= _contigWindowSequences[contigIndexLocal].size()){
					continue;
					//addCorrectedWindow(false, new DnaBitset2(contigOriginalSequence), contigIndexLocal, windowIndexLocal, outputContigFile, pass);	
					//correctedWindows[w] = new DnaBitset2(contigOriginalSequence);
					//_logFile << "No sequences for window" << endl;
					//continue;
				}

				vector<string> sequenceStrs;
				vector<Window>& sequences = _contigWindowSequences[contigIndexLocal][windowIndexLocal];


				//cout << windowIndexLocal << " " << sequences.size() << endl;
				//u_int64_t nbCorrectedWindows = 0;
				//vector<DnaBitset2*> correctedWindows(windows.size());
			
			
					
				
				u_int64_t wStart = windowIndexLocal*_windowLength;
				u_int64_t wEnd = min(_contigSequences[contigIndexLocal].size(), (size_t)(wStart+_windowLength));
				string contigOriginalSequence = _contigSequences[contigIndexLocal].substr(wStart, wEnd-wStart);
				bool isLastWindow = (windowIndexLocal == _contigWindowSequences[contigIndexLocal].size()-1);


				if(sequences.size() < 2){

					for(size_t i=0; i<sequences.size(); i++){ 
						delete sequences[i]._sequence;
					}

					addCorrectedWindow(false, new DnaBitset2(contigOriginalSequence), contigIndexLocal, windowIndexLocal, outputContigFile, pass);	
					//correctedWindows[w] = new DnaBitset2(contigOriginalSequence);
					//_logFile << "No sequences for window" << endl;
					continue;
				}
				


				std::sort(sequences.begin(), sequences.end(), [&](const Window& lhs, const Window& rhs) {
					return lhs._posStart < rhs._posStart; });


				//return;

				//_logFile << sequences.size() << endl;
				//_logFile << "1" << endl;
				spoa::Graph graph{};

				string backboneQuality = "";
				for(size_t i=0; i<contigOriginalSequence.size(); i++){
					backboneQuality += '!';
				}

				graph.AddAlignment(
					spoa::Alignment(),
					contigOriginalSequence.c_str(), contigOriginalSequence.size(),
					backboneQuality.c_str(), backboneQuality.size()
				);

				u_int32_t offset = 0.01 * contigOriginalSequence.size();


				for(size_t i=0; i<sequences.size(); i++){ 

					//size_t i = order[ii];
					const Window& window = sequences[i];
					//const DnaBitset2* dna = variant._sequence; //sequenceCopies[s._sequenceIndex];
					char* dnaStr = window._sequence->to_string();

					sequenceStrs.push_back(string(dnaStr));
					//if(windowIndexLocal == 24){
					//	cout << string(dnaStr) << endl;
					//}

					spoa::Alignment alignment;
					if (window._posStart < offset && window._posEnd > contigOriginalSequence.size() - offset) {
						alignment = alignmentEngine->Align(
							dnaStr, strlen(dnaStr),
							graph);
					} else {
						std::vector<const spoa::Graph::Node*> mapping;
						auto subgraph = graph.Subgraph(
							window._posStart,
							window._posEnd,
							&mapping);
						alignment = alignmentEngine->Align(
							dnaStr, strlen(dnaStr),
							subgraph);
						subgraph.UpdateAlignment(mapping, &alignment);
					}
					
					if (window._quality.size() == 0) {
						graph.AddAlignment(
							alignment,
							dnaStr, strlen(dnaStr));
					} else {
						graph.AddAlignment(
							alignment,
							dnaStr, strlen(dnaStr),
							window._quality.c_str(), window._quality.size());
					}
					
					free(dnaStr);
					delete window._sequence;
				}


				vector<u_int32_t> coverages;
				string correctedSequence = graph.GenerateConsensus(&coverages);

				correctedSequence = trimConsensus(correctedSequence, coverages, sequences.size(), isLastWindow);

				//for(char letter : correctedSequence){
				//	checksum += letter;
				//}

				//ofstream file("/pasteur/appa/homes/gbenoit/appa/run/correction/test_deterministic/test_humanO1_4/loulou.fasta");
				//file << ">lala" << endl;
				//file << correctedSequence << endl;
				//file.close();
				//cout << windowIndexLocal << endl;
				//getchar();
				/*
				#pragma omp critical(debug)
				{
					
					
					if(correctedSequence.size() < 450){
						cout << endl;
						cout << "-----" << endl;
						cout << "Original sequence: " << contigOriginalSequence << endl;
						cout << correctedSequence.size() << " " << contigIndexLocal << " " << windowIndexLocal << endl;
						for(size_t i=0; i<sequenceStrs.size(); i++){
							cout << i << " " << sequenceStrs[i].size() << " " << sequenceStrs[i] << endl;
						}
						getchar();
					}

				}
				*/


				addCorrectedWindow(true, new DnaBitset2(correctedSequence), contigIndexLocal, windowIndexLocal, outputContigFile, pass);

				//correctedWindows[w] = new DnaBitset2(correctedSequence);

			}
		}
		
	}

	string trimConsensus(const string& correctedSequence, const vector<u_int32_t>& coverages, const int& nbSequences, const bool& isLastWindow){

		string trimmedSequence = "";
		uint32_t average_coverage = nbSequences / 2;

		while(true){


			int32_t begin = 0, end = correctedSequence.size() - 1;
			for (; begin < static_cast<int32_t>(correctedSequence.size()); ++begin) {
				if (coverages[begin] >= average_coverage) {
					break;
				}
			}
			for (; end >= 0; --end) {
				if (coverages[end] >= average_coverage) {
					break;
				}
			}

			if (begin >= end) {
				//fprintf(stderr, "[racon::Window::generate_consensus] warning: "
				//	"contig %lu might be chimeric in window %u!\n", id_, rank_);
			} else {
				trimmedSequence = correctedSequence.substr(begin, end - begin + 1);
			}

			if(isLastWindow) break;
			if(trimmedSequence.size() > _windowLength*0.8) break;
			
			average_coverage += 1;

			if(average_coverage > nbSequences) return correctedSequence;
		}


		return trimmedSequence;
	}


	void addCorrectedWindow(bool success, DnaBitset2* seq, size_t contigIndexLocal, size_t windowIndexLocal, gzFile& outputContigFile, const int& pass){

		
		#pragma omp critical(addCorrectedWindow)
		{
			
			//cout << "Add corrected window: " << windowIndexLocal << " " << seq->m_len << endl;

			_currentContigs[contigIndexLocal].push_back({windowIndexLocal, seq, success});

			if(_currentContigs[contigIndexLocal].size() == _contigWindowSequences[contigIndexLocal].size()){
				dumpCorrectedContig(contigIndexLocal, outputContigFile, pass);
			}
		}
		
	}

	void dumpCorrectedContig(const u_int32_t& contigIndex, gzFile& outputContigFile, const int& pass){


		u_int64_t nbCorrectedWindows = 0;
		vector<CorrectedWindow>& correctedWindows = _currentContigs[contigIndex];

		for(const CorrectedWindow& correctedWindow : correctedWindows){
			if(correctedWindow._success){
				nbCorrectedWindows += 1;
			}
		}
		

		if(nbCorrectedWindows > 0){

			std::sort(correctedWindows.begin(), correctedWindows.end(), CorrectedWindowComparator);

			string contigSequence = "";
			//for(size_t w=0; w<correctedWindows.size(); w++){
			for(const CorrectedWindow& correctedWindow : correctedWindows){
				if(correctedWindow._correctedSequence == nullptr) continue;
				
				//if(correctedWindows[w] == nullptr) continue;
				char* seq = correctedWindow._correctedSequence->to_string();
				
				//if(_debug_contigName != "") cout << "Contig length: " << contigSequence.size() << "   added: " << string(seq).size() << endl;

				contigSequence += string(seq);
				free(seq);

				//cout << contigSequence << endl;
				//cout << contigSequence.size() << endl;
				//getchar();
				//delete correctedWindows[w];
			}


			u_int64_t length = contigSequence.size();
			bool isValid = true;

				
			if(_contigCoverages[contigIndex] <= 1){
				isValid = false;
			}
			else if(length < _minContigLength){
				isValid = false;
			}
			else if(length < 7500 && _contigCoverages[contigIndex] < 4){
				isValid = false;
			}

			//else if(length < 7500 && _contigCoverages[contigIndex] < 3){
			//	isValid = false;
			//}
			


			if(isValid){

				string header = _contigHeaders[contigIndex];
				
				if(_useMetamdbgHeaderStyle && pass > 0){
					
					string originalHeader = header;
					header = Utils::split(header, '_')[0]; //Remove curator subpath suffix

					bool isCircular = false;
					//if(header.find("rc") != string::npos){
					//	string h = header.substr(0, header.size()-2);
					//	header = h + "_" + to_string(_contigCoverages[contigIndex]) + "x_rc";
					//}
					//else 
					if(header[header.size()-1] == 'c'){
						isCircular = true;
						//string h = header.substr(0, header.size()-1);
						//header = h + " length=" + to_string(contigSequence.size()) + " coverage=" + to_string(_contigCoverages[contigIndex]) + " circular=yes";
					}

					header.pop_back(); //remove circular indicator
					header.erase(0, 3); //remove "ctg"
					//u_int32_t contigIndex = stoull(header);

					//else{
					//	string h = header.substr(0, header.size()-1);
					//	header = h + " length=" + to_string(contigSequence.size()) + " coverage=" + to_string(_contigCoverages[contigIndex]) + " circular=no";
					//}
					header = Utils::createContigHeader(_outputContigIndex, contigSequence.size(), _contigCoverages[contigIndex], isCircular);
					
					_file_contigHeaders << originalHeader << "\t" << _outputContigIndex << endl;
				}
				

				//cout << "Polish contig: " << contigSequence.size() << endl;
				header = ">" + header + '\n';// ">ctg" + to_string(contigIndex) + '\n';
				//header += '\n';
				gzwrite(outputContigFile, (const char*)&header[0], header.size());
				contigSequence +=  '\n';
				gzwrite(outputContigFile, (const char*)&contigSequence[0], contigSequence.size());
				//_logFile << _contigHeaders[contigIndex] << " " << contigSequence.size() << endl;
				
				_outputContigIndex += 1;

			}



		}

		for(CorrectedWindow& correctedWindow : correctedWindows){
			if(correctedWindow._correctedSequence == nullptr) continue;
			delete correctedWindow._correctedSequence;
		}
		_currentContigs.erase(contigIndex);
	}


};	

#endif 


