

#ifndef MDBG_METAG_PURGEDUPS
#define MDBG_METAG_PURGEDUPS

#include "../Commons.hpp"
#include "../utils/edlib.h"
#include "../utils/DnaBitset.hpp"



class PurgeDups : public Tool{
    
public:

	//string _inputFilename_reads;
	string _inputFilename_contigs;
	int _nbCores;
	//string _mapperOutputExeFilename;
	string _outputDir;
	string _tmpDir;
	//bool _cut_contigEnds;
	bool _cut_contigInternal;
	bool _dontOuputContigs;
	
	string _outputFilename_contigs;
	string _outputFilename_mapping;
	u_int64_t _maxBases;

	u_int64_t _minDuplicationLength_ends;
	float _minDuplicationIdentity_ends;
	u_int64_t _minDuplicationLength_internal;
	float _minDuplicationIdentity_internal;
	ofstream _file_duplicatedContigs;
	//u_int64_t _maxMemory;

	struct Alignment{
		string _contigName;
		//u_int64_t _readIndex;
		bool _strand;
		u_int32_t _readStart;
		u_int32_t _readEnd;
		u_int32_t _contigStart;
		u_int32_t _contigEnd;
		//float _score;
		//u_int64_t _length;
		
		u_int64_t length() const{
			return std::max((u_int64_t)(_readEnd - _readStart), (u_int64_t)(_contigEnd - _contigStart));
		}
		
	};
	/*
	struct Alignment{
		u_int32_t _contigIndex;
		//u_int64_t _readIndex;
		bool _strand;
		u_int32_t _readStart;
		u_int32_t _readEnd;
		u_int32_t _contigStart;
		u_int32_t _contigEnd;
		//float _score;
		//u_int64_t _length;
		
		u_int64_t length() const{
			return std::max((u_int64_t)(_readEnd - _readStart), (u_int64_t)(_contigEnd - _contigStart));
		}
		
	};
	*/

	PurgeDups(): Tool (){

	}


	void parseArgs(int argc, char* argv[]){


		string ARG_CUT_INTERNAL = "cutinternal";
		string ARG_NO_DUMP= "nodump";

		
		args::ArgumentParser parser("derep", ""); //"This is a test program.", "This goes after the options."
		args::Positional<std::string> arg_contigs(parser, "contigs", "Input contig filename", args::Options::Required);
		args::Positional<std::string> arg_outputFilename(parser, "outputFilename", "Output contig filename", args::Options::Required);
		args::Positional<std::string> arg_outputDir(parser, "outputDir", "Output dir for temporary files", args::Options::Required);
		//args::Positional<std::string> arg_contigs(parser, "contigs", "", args::Options::Required);
		//args::PositionalList<std::string> arg_readFilenames(parser, "reads", "Input filename(s) (separated by space)", args::Options::Required);
		//args::ValueFlag<int> arg_l(parser, "", "Minimizer length", {ARG_MINIMIZER_LENGTH2}, 13);
		//args::ValueFlag<float> arg_d(parser, "", "Minimizer density", {ARG_MINIMIZER_DENSITY2}, 0.005f);
		//args::ValueFlag<std::string> arg_contigs(parser, "", "", {ARG_INPUT_FILENAME_CONTIG}, "");
		args::ValueFlag<int> arg_nbCores(parser, "", "Number of cores", {ARG_NB_CORES2}, NB_CORES_DEFAULT_INT);
		args::Flag arg_cutInternal(parser, "", "", {ARG_CUT_INTERNAL});
		args::Flag arg_noDump(parser, "", "", {ARG_NO_DUMP});
		//args::Flag arg_isFinalAssembly(parser, "", "Is final multi-k pass", {ARG_FINAL});
		//args::Flag arg_firstPass(parser, "", "Is first pass of multi-k", {ARG_FIRST_PASS});
		args::Flag arg_help(parser, "", "", {'h', "help"}, args::Options::Hidden);
		//args::HelpFlag help(parser, "help", "Display this help menu", {'h'});
		//args::CompletionFlag completion(parser, {"complete"});




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

		_outputDir = args::get(arg_outputDir);
		_inputFilename_contigs = args::get(arg_contigs);
		_outputFilename_contigs = args::get(arg_outputFilename);
		_nbCores = args::get(arg_nbCores);

		_dontOuputContigs = false;
		if(arg_noDump){
			_dontOuputContigs = true;
		}

		_cut_contigInternal = false;
		if(arg_cutInternal){
			_cut_contigInternal = true;
		}

		if (_outputFilename_contigs.find(".gz") == std::string::npos) {
			_outputFilename_contigs += ".gz";
		}


		if(_inputFilename_contigs == _outputFilename_contigs){
			cerr << "Output filename == input filename" << endl;
			exit(0);
		}

		/*
		//string ARG_CUT_ENDS = "cutends";
		string ARG_CUT_INTERNAL = "cutinternal";
		string ARG_NO_DUMP= "nodump";

		//string filenameExe = argv[0];
		//_logFile << filenameExe << endl;

		//fs::path pa(filenameExe);
		//_mapperOutputExeFilename = pa.parent_path().string() + "/mapper";
		//_logFile << _mapperOutputExeFilename << endl;
		//exit(1);

		cxxopts::Options options("ToBasespace", "");
		options.add_options()
		("contigs", "", cxxopts::value<string>())
		("outputFilename", "", cxxopts::value<string>())
		//("reads", "", cxxopts::value<string>())
		("tmpDir", "", cxxopts::value<string>())
		//(ARG_CUT_ENDS, "", cxxopts::value<bool>()->default_value("false"))
		(ARG_CUT_INTERNAL, "", cxxopts::value<bool>()->default_value("false"))
		(ARG_NO_DUMP, "", cxxopts::value<bool>()->default_value("false"))
		(ARG_NB_CORES, "", cxxopts::value<int>()->default_value(NB_CORES_DEFAULT));


		options.parse_positional({"contigs", "outputFilename", "tmpDir"});
		options.positional_help("contigs outputFilename tmpDir");


		//("k,kminmerSize", "File name", cxxopts::value<std::string>())
		//("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
		//;

		if(argc <= 1){
			_logFile << options.help() << endl;
			exit(0);
		}

		cxxopts::ParseResult result;

		try{
			result = options.parse(argc, argv);

			//_inputFilename_reads = result["reads"].as<string>();
			_inputFilename_contigs = result["contigs"].as<string>();
			_outputFilename_contigs = result["outputFilename"].as<string>();
			_outputDir = result["tmpDir"].as<string>();
			//_cut_contigEnds = result[ARG_CUT_ENDS].as<bool>();
			_cut_contigInternal = result[ARG_CUT_INTERNAL].as<bool>();
			_dontOuputContigs = result[ARG_NO_DUMP].as<bool>();
			_nbCores = result[ARG_NB_CORES].as<int>();
		}
		catch (const std::exception& e){
			std::_logFile << options.help() << std::endl;
			std::cerr << e.what() << std::endl;
			std::exit(EXIT_FAILURE);
		}
		*/




		//fs::path p(_inputFilename_contigs);
		//while(p.has_extension()){
		//	p.replace_extension("");
		//}

		_tmpDir = _outputDir + "/tmp/";
		if(!fs::exists(_tmpDir)){
			fs::create_directories(_tmpDir);
		}

		fs::path outputContigPath(_outputFilename_contigs);
		string contigDir = outputContigPath.parent_path();
		if(!fs::exists(contigDir)){
			fs::create_directories(contigDir);
		}

		//_outputFilename_contigs = p.string() + "_derep.fasta.gz";
		_maxBases = 200000000ull;
		_minDuplicationLength_ends = 1000;
		_minDuplicationIdentity_ends = 0.95;
		_minDuplicationLength_internal = 10000;
		_minDuplicationIdentity_internal = 0.95;
		//_outputFilename_contigs = p.string() + "_corrected.fasta.gz";
		_outputFilename_mapping = _tmpDir + "/_tmp_mapping_derep__.paf.gz";
		//_maxMemory = 4000000000ull;

		_file_duplicatedContigs = ofstream(_outputDir + "/duplicatedContigs.txt");
		
		
		openLogFile(_tmpDir);

		_logFile << "Contigs: " << _inputFilename_contigs << endl;
		//_logFile << "Cut contigs ends: " << _cut_contigEnds << endl;
		_logFile << "Cut contigs internal: " << _cut_contigInternal << endl;
		_logFile << "Output filename: " << _outputFilename_contigs << endl;
		
		_fields = new vector<string>();
		_fields_optional = new vector<string>();
	}

	gzFile _outputContigFile;
	unordered_map<string, string> _contigSequences;
	unordered_map<string, vector<Alignment>> _alignments;

	/*
	//unordered_map<string, string> _debug_contigIndex_to_contigName;

	void debug_indexContigName(){
		
		auto fp = std::bind(&PurgeDups::indexContigName_read, this, std::placeholders::_1);
		ReadParser readParser(_inputFilename_contigs, true, false);
		readParser.parse(fp);
	}

	void indexContigName_read(const Read& read){
		_debug_contigIndex_to_contigName[read._index] = read._header;
	}
	*/

    void execute (){

		//debug_indexContigName();
		mapReads();
		processContigs();
		
		dumpDereplicatedContigs();
		//fs::remove_all(_tmpDir);

		closeLogFile();
	}

	void mapReads(){
		
		cerr << "Mapping reads" << endl;

		string inputContigsFilename = _tmpDir + "/input_contigs.txt";
		ofstream input(inputContigsFilename);
		input << _inputFilename_contigs << endl;
		input.close();

		string command = "minimap2 -x asm20 -H -DP --dual=no -I 500M -K 500M -t " + to_string(_nbCores) + " " + _inputFilename_contigs + " " + _inputFilename_contigs;
		command += " | gzip -c - > " + _outputFilename_mapping;
		//command += " > " + _outputFilename_mapping;
		//command += " | " + _mapperOutputExeFilename + " " + _inputFilename_contigs + " " + inputContigsFilename + " " + _outputFilename_mapping;
		Utils::executeCommand(command, _outputDir, _logFile);
	}

	u_int64_t _currentLoadedBases;
	unordered_map<string, vector<DbgEdge>> _duplicationEnds;
	unordered_map<string, vector<DbgEdge>> _duplicationInternal;
	//unordered_set<DbgEdge, hash_pair> _performedPairs;

	void processContigs(){

		auto fp = std::bind(&PurgeDups::processContigs_read, this, std::placeholders::_1);
		ReadParser readParser(_inputFilename_contigs, true, false, _logFile);
		readParser.parse(fp);

		if(_currentLoadedBases > 0){
			processPass();
		}
	}

	
	void processContigs_read(const Read& read){

		const string& contigName = Utils::shortenHeader(read._header);
		//u_int32_t contigIndex = read._index;

		_contigSequences[contigName] = read._seq;
		_currentLoadedBases += read._seq.size();

		if(_currentLoadedBases > _maxBases){
			processPass();
		}

	}

	void processPass(){
		cerr << "Processing " <<  _contigSequences.size() << " contigs" << endl;
		indexAlignments();
		detectDuplication();
		clearPass();
		cerr << "Pass done" << endl;
	}

	void clearPass(){
		_alignments.clear();
		_contigSequences.clear();
		_currentLoadedBases = 0;
	}

	void indexAlignments(){
		cerr << "\tIndexing alignments" << endl;

		PafParser pafParser(_outputFilename_mapping);
		auto fp = std::bind(&PurgeDups::indexAlignments_read, this, std::placeholders::_1);
		pafParser.parse(fp);
	}

	vector<string>* _fields;
	vector<string>* _fields_optional;

	void indexAlignments_read(const string& line){


		double errorThreshold = 0.3;

		//_logFile << lineInput << endl;
		//getchar();


		GfaParser::tokenize(line, _fields, '\t');

		//_logFile << line << endl;

		const string& readName = Utils::shortenHeader((*_fields)[0]);
		const string& contigName = Utils::shortenHeader((*_fields)[5]);
		

		u_int32_t readStart = stoull((*_fields)[2]);
		u_int32_t readEnd = stoull((*_fields)[3]);
		u_int32_t contigLength = stoull((*_fields)[6]);
		u_int32_t contigStart = stoull((*_fields)[7]);
		u_int32_t contigEnd = stoull((*_fields)[8]);

		u_int64_t nbMatches = stoull((*_fields)[9]);
		u_int64_t alignLength = stoull((*_fields)[10]);
		u_int64_t queryLength = stoull((*_fields)[1]);

		bool strand = (*_fields)[4] == "-";

		if(readName == contigName) return;
		if(_contigSequences.find(contigName) == _contigSequences.end() && _contigSequences.find(readName) == _contigSequences.end()) return;

		u_int32_t length = std::max((u_int64_t)(readEnd - readStart), (u_int64_t)(contigEnd - contigStart));
		if(length < 1000) return;

		Alignment align = {contigName, strand, readStart, readEnd, contigStart, contigEnd}; //, score
		_alignments[readName].push_back(align);
		/*
		u_int32_t length = std::max((u_int64_t)(readEnd - readStart), (u_int64_t)(contigEnd - contigStart));

		if(_alignments.find(readIndex) == _alignments.end()){
			_alignments[readIndex] = align;
		}
		else{

			if(length > align.length()){
				_alignments[readIndex] = align;
			}
			
		}
		*/

		
	}

	void detectDuplication(){
		_logFile << "\tDetecting duplication" << endl;

		//auto fp = std::bind(&ContigPolisher::collectWindowCopies_read, this, std::placeholders::_1);
		//ReadParser readParser(_inputFilename_reads, false, false);
		//readParser.parse(fp);


		//const string& partitionFilename = _tmpDir + "/part_" + to_string(partition) + ".gz";
		ReadParserParallel readParser(_inputFilename_contigs, true, false, _nbCores, _logFile);
		readParser.parse(ContigAlignerFunctor(*this));
	}



	class ContigAlignerFunctor {

		public:

		PurgeDups& _purgeDups;
		//unordered_map<string, u_int32_t>& _contigName_to_contigIndex;
		//unordered_map<string, u_int64_t>& _readName_to_readIndex;
		//unordered_map<u_int64_t, vector<Alignment>>& _alignments;
		unordered_map<string, vector<Alignment>>& _alignments;
		unordered_map<string, string>& _contigSequences;
		//unordered_map<u_int32_t, vector<vector<Window>>>& _contigWindowSequences;
		//size_t _windowLength;


		ContigAlignerFunctor(PurgeDups& purgeDups) : _purgeDups(purgeDups), _alignments(purgeDups._alignments), _contigSequences(purgeDups._contigSequences){
		}

		ContigAlignerFunctor(const ContigAlignerFunctor& copy) : _purgeDups(copy._purgeDups), _alignments(copy._alignments), _contigSequences(copy._contigSequences){
			
		}

		~ContigAlignerFunctor(){
		}

		void operator () (const Read& read) {

			
			const string& readName = Utils::shortenHeader(read._header);

			//u_int64_t readIndex = stoull(read._header);

			//#pragma omp critical
			//{
				//_logFile << read._index << endl;
			//}

			if(_alignments.find(readName) == _alignments.end()) return;

			//_logFile << "\t" << _alignments[read._index].size() << endl;
			for(const Alignment& al : _alignments[readName]){

				//if(read._index == al._contigIndex){

					//_logFile << "Self AL not normal" << endl;
				//	continue;
				//}
				
				if(_contigSequences.find(al._contigName) == _contigSequences.end()) continue;

				//_logFile << al._readStart << " " << al._readEnd << endl;
				processAlignment(read, al);
				//_logFile << read._index << " " << al._contigIndex << endl;
			}
			//_logFile << "\tdone" << endl;
			/*
			//if(_contigPolisher._currentPartition == 0) _logFile << readIndex << " " << (_alignments.find(readIndex) != _alignments.end()) << endl;
			
			//if(readIndex % 100000 == 0) _logFile << "\t" << readIndex << endl;

			if(_alignments.find(readIndex) == _alignments.end()) return;

			//const vector<Alignment>& als = _alignments[readIndex];
			const Alignment& al = _alignments[readIndex];
			//for(const Alignment& al : _alignments[readIndex]){
			u_int64_t contigIndex = al._contigIndex;

			if(_contigSequences.find(contigIndex) == _contigSequences.end()) return;

			//_logFile << read._seq.size() << " " << read._qual.size() << " " << _contigSequences[contigIndex].size() << " " << al._readStart << " " << al._readEnd << " " << al._contigStart << " " << al._contigEnd << endl;
			string readSeq = read._seq;
			string qualSeq = read._qual;
			string readSequence = readSeq.substr(al._readStart, al._readEnd-al._readStart);
			string contigSequence = _contigSequences[contigIndex].substr(al._contigStart, al._contigEnd-al._contigStart);



			if(al._strand){
				Utils::toReverseComplement(readSequence);
				Utils::toReverseComplement(readSeq);
				std::reverse(qualSeq.begin(), qualSeq.end());
			}
			

			//_logFile << readSequence << endl;
			//_logFile << contigSequence << endl;

			//_logFile << contigSequence.size() << " "<< readSequence.size() << endl;
			static EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0);


			EdlibAlignResult result = edlibAlign(readSequence.c_str(), readSequence.size(), contigSequence.c_str(), contigSequence.size(), config);


			char* cigar;

			if (result.status == EDLIB_STATUS_OK) {
				cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
			} else {
				_logFile << "Invalid edlib results" << endl;
				exit(1);
			}

			//_logFile << cigar << endl;

			edlibFreeAlignResult(result);
			
			find_breaking_points_from_cigar(_windowLength, al, readSeq.size(), cigar, readSeq, qualSeq);
			free(cigar);

			//getchar();

			//}
			*/
		}

		u_int64_t _maxHang = 1000;

		void processAlignment(const Read& read, const Alignment& al){

			const string& readIndex = Utils::shortenHeader(read._header);
			const string& contigIndex = al._contigName;

			//u_int32_t readIndex = read._index;
			//u_int32_t contigIndex = al._contigIndex;


			const string& contigSequenceComplete = _contigSequences[contigIndex];


			//#pragma omp critical
			//{
			//	_logFile << readName << " " << read._seq.size() << " " << al._readStart << " " << al._readEnd << "    " << contigIndex << " " << contigSequenceComplete.size() << " " << al._contigStart << " " << al._contigEnd << endl;
			//}

			if(al.length() > 400000) return;

			string readSeq = read._seq;
			string readSequence = readSeq.substr(al._readStart, al._readEnd-al._readStart);
			string contigSequence = contigSequenceComplete.substr(al._contigStart, al._contigEnd-al._contigStart);

			//_logFile << readSequence.size() << " " << contigSequence.size() << endl;


			if(al._strand){
				Utils::toReverseComplement(readSequence);
				//Utils::toReverseComplement(readSeq);
			}

			//_logFile << contigSequence << endl;
			//_logFile << readSequence << endl;
			//exit(1);



			//_logFile << cigar << endl;
			//_logFile << contigSeq_al << endl;
			//_logFile << readSeq_al << endl;
			//_logFile << contigSeq_al.size() << endl;
			//_logFile << readSeq_al.size() << endl;

			bool isOverlap = false;

			if(contigSequenceComplete.size() < read._seq.size()){

				vector<u_int8_t> isMatches_read;
				vector<u_int8_t> isMatches_contig;
				performAlignement(readSequence, contigSequence, isMatches_read, isMatches_contig);

				bool allowOverlapLeft = true;
				bool allowOverlapRight = true;
				isOverlap = checkOverlap(contigIndex, al._contigStart, al._contigEnd, contigSequenceComplete.size(), isMatches_contig, allowOverlapLeft, allowOverlapRight);
				//if(!isOverlap) checkOverlap(readIndex, al._readStart, al._readEnd, read._seq.size(), isMatches_read, allowOverlapLeft, allowOverlapRight);
			}
			else{

				vector<u_int8_t> isMatches_read;
				vector<u_int8_t> isMatches_contig;
				performAlignement(readSequence, contigSequence, isMatches_read, isMatches_contig);

				bool allowOverlapLeft = true;
				bool allowOverlapRight = true;
				isOverlap = checkOverlap(readIndex, al._readStart, al._readEnd, read._seq.size(), isMatches_read, allowOverlapLeft, allowOverlapRight);
				//if(!isOverlap) checkOverlap(contigIndex, al._contigStart, al._contigEnd, contigSequenceComplete.size(), isMatches_contig, allowOverlapLeft, allowOverlapRight);
			}

			if(_purgeDups._cut_contigInternal){
				if(!isOverlap){
					if(contigSequenceComplete.size() < read._seq.size()){
						
						vector<u_int8_t> isMatches_read;
						vector<u_int8_t> isMatches_contig;
						performAlignement(readSequence, contigSequence, isMatches_read, isMatches_contig);

						//if(contigSequenceComplete.size() < 30000){
							checkInternalDuplication(contigIndex, al._contigStart, al._contigEnd, isMatches_contig);
						//}
					}
					else{
						
						vector<u_int8_t> isMatches_read;
						vector<u_int8_t> isMatches_contig;
						performAlignement(readSequence, contigSequence, isMatches_read, isMatches_contig);

						//if(readSeq.size() < 30000){
							checkInternalDuplication(readIndex, al._readStart, al._readEnd, isMatches_read);
						//}
					}
				}
			}


			//_logFile << cigar << endl;



			//getchar();
			//_logFile << cigar << endl;
			//exit(1);

		}

		void performAlignement(const string& readSequence, const string& contigSequence, vector<u_int8_t>& isMatches_read, vector<u_int8_t>& isMatches_contig){

			static EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0);


			//_logFile << "Comparing: " << readSequence.size() << " " << contigSequence.size() << endl;

			EdlibAlignResult result = edlibAlign(readSequence.c_str(), readSequence.size(), contigSequence.c_str(), contigSequence.size(), config);

			//_logFile << result.alignment << endl;
			//_logFile << result.alignmentLength << endl;

			char* cigar; 

			if (result.status == EDLIB_STATUS_OK) {
				cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
			} else {
				_purgeDups._logFile << "Invalid edlib results" << endl;
				exit(1);
			}






			//string contigSeq_al = "";
			//string readSeq_al = "";
			int contigPos = 0;
			int readPos = 0;

			for (uint32_t i = 0, j = 0; i < strlen(cigar); ++i) {
				if (cigar[i] == 'M' || cigar[i] == '=' || cigar[i] == 'X') {
					uint32_t k = 0, num_bases = atoi(&cigar[j]);
					j = i + 1;

					for(size_t i=0; i<num_bases; i++){

						if(contigSequence[contigPos] == readSequence[readPos]){
							isMatches_read.push_back(0);
							isMatches_contig.push_back(0);
						}
						else{
							isMatches_read.push_back(1);
							isMatches_contig.push_back(1);
						}
						//isMatches.push_back(true);
						//contigSeq_al += contigSequence[contigPos];
						//readSeq_al += readSequence[readPos];
						contigPos += 1;
						readPos += 1;
					}

					//_logFile << "Match: " << num_bases << endl;
				} else if (cigar[i] == 'I') {
					uint32_t k = 0, num_bases = atoi(&cigar[j]);
					j = i + 1;

					//for(size_t i=0; i<num_bases; i++){
					//	seq2 += readSequence[readPos];
					//	seq1 += "-";
					//}

				

					//if(countInsert){
					for(size_t i=0; i<num_bases; i++){
						
						//readSeq_al += readSequence[readPos];
						//contigSeq_al += "-";
						readPos += 1;
						isMatches_read.push_back(2);
						
						//readPos += 1;
					}
					//}

					/*
					for(size_t i=0; i<num_bases; i++){
						isMatches.push_back(true);
						seq1 += contigSequence[contigPos];
						seq2 += readSequence[readPos];
						contigPos += 1;
						readPos += 1;
					}
					*/

					//_logFile << "Insert: " << num_bases << endl;

				} else if (cigar[i] == 'D' || cigar[i] == 'N') {
					uint32_t k = 0, num_bases = atoi(&cigar[j]);
					j = i + 1;
					//for(size_t i=0; i<num_bases; i++){
					//	seq1 += contigSequence[readPos];
					//	seq2 += "-";
					//}
					//if(!countInsert){

					for(size_t i=0; i<num_bases; i++){

						//contigSeq_al += contigSequence[readPos];
						//readSeq_al += "-";

						contigPos += 1;
						isMatches_contig.push_back(2);
					}
					//}
					//_logFile << "Del: " << num_bases << endl;

				} else if (cigar[i] == 'S' || cigar[i] == 'H' || cigar[i] == 'P') {
					//uint32_t num_bases = atoi(&cigar_[j]);
					j = i + 1;
					//isMatches.push_back(num_bases);
				}
				//else{
				//	_logFile << cigar_[i] << endl;
				//	_logFile << "mdr" << endl;
				//}
			}

			//_logFile << "lala: " << isMatches_read.size() << " " << isMatches_contig.size() << endl;
			edlibFreeAlignResult(result);
			free(cigar);
		}

		void checkInternalDuplication(const string& targetIndex, u_int32_t alignStart, u_int32_t alignEnd, const vector<u_int8_t>& isMatches){

			size_t alignLength = alignEnd - alignStart;
			if(alignLength < _purgeDups._minDuplicationLength_internal) return;

			//_logFile << "check internal: " << _purgeDups._debug_contigIndex_to_contigName[targetIndex] << " " << alignStart << " " << alignEnd << " " << alignLength << endl;

			vector<bool> isDuplicated(isMatches.size(), false);

			updateInternalDuplication(0, alignLength, isMatches, isDuplicated);
			/*
			size_t w = 1000;
			size_t posStart = 0;
			while(true){

				size_t posEnd = posStart + _purgeDups._minDuplicationLength_internal;
				if(posEnd >= isMatches.size()) break;

				updateInternalDuplication(posStart, posEnd, isMatches, isDuplicated);

				posStart += w;
			}

			updateInternalDuplication(isMatches.size()-_purgeDups._minDuplicationLength_internal, isMatches.size(), isMatches, isDuplicated);
			*/

			//for(size_t i=0; i<isDuplicated.size(); i++){
				//_logFile << isDuplicated[i];
			//}
			//_logFile << endl;




			//_logFile << "lala: " << isMatches.size() << endl;

			bool isDuplicatedArea = false;
			u_int32_t startPos = 0;
			//size_t subSeqIndex = 0;
			u_int32_t endPos = 0;

			for(long i=0; i<isDuplicated.size(); i++){
				
				endPos = i;

				if(!isDuplicated[i]){
					if(isDuplicatedArea){
						//long length = endPos-startPos;
						//_logFile << length << endl;
						#pragma omp critical(dup)
						{
							_purgeDups._logFile << "Internal duplication: " << alignStart+startPos << " " << alignStart+endPos << endl;
							_purgeDups._duplicationInternal[targetIndex].push_back({alignStart+startPos, alignStart+endPos});
						}

					}
					isDuplicatedArea = false;
				}
				else{
					if(!isDuplicatedArea){
						startPos = i;
					}
					isDuplicatedArea = true;
				}
			}

			if(isDuplicatedArea){
				//long length = endPos-startPos;
				//_logFile << length << endl;

				#pragma omp critical(dup)
				{
					_purgeDups._logFile << "Internal duplication: " << alignStart+startPos << " " << alignStart+endPos << endl;
					_purgeDups._duplicationInternal[targetIndex].push_back({alignStart+startPos, alignStart+endPos});
				}
			}


		}

		void updateInternalDuplication(size_t posStart, size_t posEnd, const vector<u_int8_t>& isMatches, vector<bool>& isDuplicated){

			float identity = computeIdentity(isMatches);
			//_logFile << identity << endl;

			if(identity >= _purgeDups._minDuplicationIdentity_internal){
				for(size_t p=posStart; p<posEnd; p++){
					isDuplicated[p] = true;
				}
			}
		}

		float computeIdentity(const vector<u_int8_t>& isMatches){

			double nbMatches = 0;
			double length = 0;
			bool isGap = false;

			for(size_t i=0; i<isMatches.size(); i++){
				if(isMatches[i] == 0){
					isGap = false;
					nbMatches += 1;
					length += 1;
				}
				else if(isMatches[i] == 1){
					isGap = false;
					length += 1;
				}
				else if(isMatches[i] == 2){
					if(!isGap){
						length += 1;
					}
					isGap = true;
				}
				//length += 1;
			}

			double identity = nbMatches / length; //(posEnd-posStart);
			return identity;
		}

		bool checkOverlap(const string& targetIndex, u_int32_t alignStart, u_int32_t alignEnd, u_int32_t seqLength, const vector<u_int8_t>& isMatches, bool& allowOverlapLeft, bool& allowOverlapRight){

			u_int64_t hangLeft = alignStart;
			u_int64_t hangRight = seqLength - alignEnd;

			//_logFile << hangLeft << " " << hangRight << endl;
			//DbgEdge edge = {targetIndex, queryIndex};
			//edge = edge.normalize();
			bool isOverlap = false;

			if((hangLeft < _maxHang) ){
				u_int64_t overlapLength = determineOverlapLength(isMatches, true);
				//if(_performedPairs.find(edge) == _performedPairs.end()){
				//	_duplicationBounds[targetIndex].push_back({targetStart, targetEnd});
				//	isOverlap = true;
				//}
				if(overlapLength > _purgeDups._minDuplicationLength_ends){
					allowOverlapLeft = false;
					isOverlap = true;

					#pragma omp critical(dup)
					{
						//_logFile << "\tleft overlap: " << overlapLength << " " << alignStart << " " << alignEnd << endl;
						_purgeDups._logFile << "\tleft overlap: " << targetIndex << " " << alignStart << " " << alignEnd << endl;
						_purgeDups._duplicationEnds[targetIndex].push_back({alignStart, alignEnd});
					}
				}
			}

			if((hangRight < _maxHang)){
				u_int64_t overlapLength = determineOverlapLength(isMatches, false);

				if(overlapLength > _purgeDups._minDuplicationLength_ends){
					allowOverlapRight = false;
					isOverlap = true;
					#pragma omp critical(dup)
					{
						//_logFile << "\tright overlap: " << overlapLength << " " << alignStart << " " << alignEnd << endl;
						_purgeDups._logFile << "\tright overlap: " << targetIndex << " " << alignStart << " " << alignEnd << endl;
						_purgeDups._duplicationEnds[targetIndex].push_back({alignStart, alignEnd});
					}
				}

				//_logFile << "right overlap: " << overlapLength << endl;
				//if(_performedPairs.find(edge) == _performedPairs.end()){
				//	_duplicationBounds[targetIndex].push_back({targetStart, targetEnd});
				//	isOverlap = true;
				//}
			}



			/*
			u_int64_t targetIndex = Utils::contigName_to_contigIndex(contigName);
			u_int64_t queryIndex = Utils::contigName_to_contigIndex(readName);
			DbgEdge edge = {targetIndex, queryIndex};
			edge = edge.normalize();

			bool isOverlap = false;

			if(hangLeft < maxHang){
				if(_performedPairs.find(edge) == _performedPairs.end()){
					_duplicationBounds[targetIndex].push_back({targetStart, targetEnd});
					isOverlap = true;
				}
			}

			if(hangRight < maxHang){
				if(_performedPairs.find(edge) == _performedPairs.end()){
					_duplicationBounds[targetIndex].push_back({targetStart, targetEnd});
					isOverlap = true;
				}
			}

			if(isOverlap){
				_performedPairs.insert(edge);
			}
			*/


			//find_breaking_points_from_cigar(al, readSeq.size(), cigar, readSeq);

			return isOverlap;
		}



		long determineOverlapLength(const vector<u_int8_t>& isMatches, bool isLeftOverlap){
			
			if(isMatches.size() < _purgeDups._minDuplicationLength_ends) return 0;

			float identity = computeIdentity(isMatches);

			if(identity >= _purgeDups._minDuplicationIdentity_ends){
				return isMatches.size();
			}

			return 0;

			/*
			double nbMatches = 0;
			//double length = 0;

			for(size_t i=0; i<isMatches.size(); i++){
				if(isMatches[i]) nbMatches += 1;
				//length += 1;
			}

			double identity = nbMatches / (isMatches.size());
			//_logFile << identity << endl;

			if(identity >= _purgeDups._minDuplicationIdentity_ends){
				return isMatches.size();
				//for(size_t p=posStart; p<posEnd; p++){
				//	isDuplicated[p] = true;
				//}
			}

			return 0;
			*/

			//if(!isLeftOverlap){
			//	std::reverse(alTypes.begin(), alTypes.end());
			//}
			//_logFile << nbMatches / al.length() << endl;
			
			/*
			double nbMatches = 0;
			double currentLength = 0;

			if(isLeftOverlap){
				for(size_t i=0; i<isMatches.size(); i++){


					if(currentLength == _purgeDups._minDuplicationLength_ends){
						double identity = nbMatches / currentLength;
						if(identity < _purgeDups._minDuplicationIdentity_ends){
							return 0;
						}
					}
					else if(currentLength > _purgeDups._minDuplicationLength_ends){
						double identity = nbMatches / currentLength;
						//_logFile << currentLength << " " << identity << endl;
						if(identity < _purgeDups._minDuplicationIdentity_ends){
							return currentLength;
						}
					}

					if(isMatches[i]){
						nbMatches += 1;
					}

					currentLength += 1;
				}
			}
			else{
				for(long i=isMatches.size(); i>= 0; i--){


					if(currentLength == _purgeDups._minDuplicationLength_ends){
						double identity = nbMatches / currentLength;
						if(identity < _purgeDups._minDuplicationIdentity_ends){
							return 0;
						}
					}
					else if(currentLength > _purgeDups._minDuplicationLength_ends){
						double identity = nbMatches / currentLength;
						//_logFile << currentLength << " " << identity << endl;
						if(identity < _purgeDups._minDuplicationIdentity_ends){
							return currentLength;
						}
					}

					if(isMatches[i]){
						nbMatches += 1;
					}

					currentLength += 1;
				}
			}

			if(currentLength >= _purgeDups._minDuplicationLength_ends){
				return currentLength;
			}

			return 0;
			*/
		}

	};

	void dumpDereplicatedContigs(){
		
		_logFile << "Writing dereplicated contigs" << endl;

		_outputContigFile = gzopen(_outputFilename_contigs.c_str(),"wb");


		auto fp = std::bind(&PurgeDups::dumpDereplicatedContigs_read, this, std::placeholders::_1);
		ReadParser readParser(_inputFilename_contigs, true, false, _logFile);
		readParser.parse(fp);

		gzclose(_outputContigFile);
		_file_duplicatedContigs.close();
	}


	void dumpDereplicatedContigs_read(const Read& read){

		const string& contigName = Utils::shortenHeader(read._header);
		//u_int64_t contigIndex = read._index; //Utils::contigName_to_contigIndex(read._header);

		string seq = read._seq;

		//if(_duplicationBounds.find(contigIndex) != _duplicationBounds.end()){

			//if(read._header == "ctg89567") 
			//_logFile << "-----" << endl;
		vector<bool> isDuplicated(seq.size(), false);

		for(const DbgEdge& duplicatedBound : _duplicationEnds[contigName]){
			for(size_t i=duplicatedBound._from; i<duplicatedBound._to; i++){
				isDuplicated[i] = true;
			}
		}
		for(const DbgEdge& duplicatedBound : _duplicationInternal[contigName]){
			for(size_t i=duplicatedBound._from; i<duplicatedBound._to; i++){
				isDuplicated[i] = true;
			}
		}
		
		double nbDuplicatedPos = 0;
		for(size_t i=0; i<isDuplicated.size(); i++){
			if(isDuplicated[i]){
				nbDuplicatedPos += 1;
			}
		}

		_logFile << read._header << " " << read._seq.size() << " " << (nbDuplicatedPos / isDuplicated.size()) << endl;

		float duplicatuionRate = nbDuplicatedPos / isDuplicated.size();
		if(duplicatuionRate > 0.95){
			_file_duplicatedContigs << Utils::contigName_to_contigIndex(read._header) << endl;
			return;
		}
		if(_dontOuputContigs) return;










		isDuplicated.clear();

		for(const DbgEdge& duplicatedBound : _duplicationEnds[contigName]){
			for(size_t i=duplicatedBound._from; i<duplicatedBound._to; i++){
				isDuplicated[i] = true;
			}
		}
		if(_cut_contigInternal){
			for(const DbgEdge& duplicatedBound : _duplicationInternal[contigName]){
				for(size_t i=duplicatedBound._from; i<duplicatedBound._to; i++){
					isDuplicated[i] = true;
				}
			}
		}


		
		bool isDuplicatedArea = true;
		long startPos = 0;
		size_t subSeqIndex = 0;
		long endPos = 0;

		for(long i=0; i<seq.size(); i++){
			
			endPos = i;

			if(isDuplicated[i]){
				if(!isDuplicatedArea){
					long length = endPos-startPos-1;
					if(length >= 500){
						string header = read._header;


						char lastChar = header[header.size()-1];
						header.pop_back();

						header += "_" + to_string(subSeqIndex) + lastChar;
						
						string subSeq = seq.substr(startPos, length);

						header = ">" + header + '\n';
						gzwrite(_outputContigFile, (const char*)&header[0], header.size());
						string contigSequence = subSeq + '\n';
						gzwrite(_outputContigFile, (const char*)&contigSequence[0], contigSequence.size());


						//if(read._header == "ctg89567")
						//_logFile << "Dump area: " << startPos << " " << startPos+length << endl;

						subSeqIndex += 1;
					}


				}
				isDuplicatedArea = true;
			}
			else{
				if(isDuplicatedArea){
					startPos = i;
				}
				isDuplicatedArea = false;
			}
		}

		if(!isDuplicatedArea){
			long length = endPos-startPos-1;
			if(length > 500){
				string header = read._header;
				char lastChar = header[header.size()-1];
				header.pop_back();

				header += "_" + to_string(subSeqIndex) + lastChar;
				
				string subSeq = seq.substr(startPos, length);

				header = ">" + header + '\n';
				gzwrite(_outputContigFile, (const char*)&header[0], header.size());
				string contigSequence = subSeq + '\n';
				gzwrite(_outputContigFile, (const char*)&contigSequence[0], contigSequence.size());
				
				//_logFile << "Dump area: " << startPos << " " << startPos+length << endl;
			}


		}
			

		//}
		//else{
		//	string header = ">" + read._header + '\n';
		//	gzwrite(_outputContigFile, (const char*)&header[0], header.size());
		//	string contigSequence = seq + '\n';
		//	gzwrite(_outputContigFile, (const char*)&contigSequence[0], contigSequence.size());
		//}





	}

};	

#endif 


