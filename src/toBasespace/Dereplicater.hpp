

#ifndef MDBG_METAG_DEREPLICATER
#define MDBG_METAG_DEREPLICATER

#include "../Commons.hpp"

class Dereplicater : public Tool{
    
public:

	string _inputFilename_contigs;
	string _outputFilename_contigs;
	string _outputDir;
	string _tmpDir;
	int _nbCores;
	float _minIdentity;
	u_int64_t _minCircularLength;
	u_int64_t _minLinearLength;

	gzFile _outputContigFile;

	Dereplicater(): Tool (){


	}

	void parseArgs(int argc, char* argv[]){

		char ARG_MIN_IDENTITY = 'i';
		char ARG_CIRCULAR_LENGTH = 'c';
		char ARG_LINEAR_LENGTH = 'l';
		
		args::ArgumentParser parser("derepOld", ""); //"This is a test program.", "This goes after the options."
		args::Positional<std::string> arg_contigs(parser, "contigs", "Input contig filename", args::Options::Required);
		args::Positional<std::string> arg_outputFilename(parser, "outputFilename", "Output contig filename", args::Options::Required);
		args::Positional<std::string> arg_outputDir(parser, "outputDir", "Output dir for temporary files", args::Options::Required);
		args::ValueFlag<float> arg_minIdentity(parser, "", "Minimum identity for strains (0-1)", {ARG_MIN_IDENTITY}, 0.99);
		args::ValueFlag<int> arg_circular(parser, "", "Leave untouched circular contigs with length > minLength (0 = disable)", {ARG_CIRCULAR_LENGTH}, 0);
		args::ValueFlag<int> arg_linear(parser, "", "Leave untouched linear contigs with length > minLength (0 = disable)", {ARG_LINEAR_LENGTH}, 0);
		args::ValueFlag<int> arg_nbCores(parser, "", "Number of cores", {ARG_NB_CORES2}, NB_CORES_DEFAULT_INT);
		args::Flag arg_help(parser, "", "", {'h', "help"}, args::Options::Hidden);

		try
		{
			parser.ParseCLI(argc, argv);
		}
		catch (const args::Help&)
		{
			cerr << parser;
			exit(0);
		}
		catch (const std::exception& e)
		{
			cerr << parser;
			//_logFile << endl;
			cerr << e.what() << endl;
			exit(0);
		}

		if(arg_help){
			cerr << parser;
			exit(0);
		}



		_inputFilename_contigs = args::get(arg_contigs);
		_outputFilename_contigs = args::get(arg_outputFilename);
		_outputDir = args::get(arg_outputDir);
		_nbCores = args::get(arg_nbCores);
		_minCircularLength = args::get(arg_circular);
		_minLinearLength = args::get(arg_linear);
		_minIdentity = args::get(arg_minIdentity);
		_minIdentity *= 100;

		if (_outputFilename_contigs.find(".gz") == std::string::npos) {
			_outputFilename_contigs += ".gz";
		}


		if(_inputFilename_contigs == _outputFilename_contigs){
			cerr << "Output filename == input filename" << endl;
			exit(0);
		}
		//_tmpDir = _outputDir + "/__tmp_derepOld/";



		/*
		cxxopts::Options options("ToBasespace", "");
		options.add_options()
		("contigs", "", cxxopts::value<string>())
		("outputFilenme", "", cxxopts::value<string>())
		("tmpDir", "", cxxopts::value<string>())
		(ARG_NB_CORES, "", cxxopts::value<int>()->default_value(NB_CORES_DEFAULT));

		options.parse_positional({"contigs", "outputFilenme", "tmpDir"});
		options.positional_help("contigs outputFilenme tmpDir");

		if(argc <= 1){
			_logFile << options.help() << endl;
			exit(0);
		}

		cxxopts::ParseResult result;

		try{
			result = options.parse(argc, argv);

			_inputFilename_contigs = result["contigs"].as<string>();
			_outputFilename_contigs = result["outputFilenme"].as<string>();
			_tmpDir = result["tmpDir"].as<string>();;
			_nbCores = result[ARG_NB_CORES].as<int>();
		}
		catch (const std::exception& e){
			cerr << options.help() << std::endl;
			cerr << e.what() << std::endl;
			std::exit(EXIT_FAILURE);
		}
		*/

		//if (_outputFilename_contigs.find(".gz") == std::string::npos) {
		//	_outputFilename_contigs += ".gz";
		//}

		//_outputFilename_contigs = _outputDir + "/contigs_derepOld.fasta.gz";

		_tmpDir = _outputDir + "/tmp/";
		if(!fs::exists(_tmpDir)){
			fs::create_directories(_tmpDir);
		}

		fs::path outputContigPath(_outputFilename_contigs);
		string contigDir = outputContigPath.parent_path();
		if(!fs::exists(contigDir)){
			fs::create_directories(contigDir);
		}

		openLogFile(_tmpDir);

		_logFile << endl;
		_logFile << "Contigs: " << _inputFilename_contigs << endl;
		_logFile << "Output contigs: " << _outputFilename_contigs << endl;
		_logFile << "Tmp dir: " << _tmpDir << endl;
		_logFile << endl;

	}

	struct ContigOverlap{
		u_int64_t _leftPos;
		u_int64_t _rightPos;

		ContigOverlap(){
			_leftPos = 0;
			_rightPos = -1;
		}
	};

	unordered_map<u_int64_t, ContigOverlap> _contigOverlaps;

    void execute (){
		
		string outputMappingFilename = _tmpDir + "/derep_align.paf";

		mapContigs(outputMappingFilename);
		detectDuplicatedContigs(outputMappingFilename);
		dumpDereplicatedContigs();
		

		//fs::remove(outputMappingFilename);
		closeLogFile();
	}

	//zFile _queryContigFile;
	unordered_set<u_int32_t> _duplicatedContigIndex;
	unordered_map<u_int32_t, vector<DbgEdge>> _duplicationBounds;

	void mapContigs(const string& mapFilename){

		//fs::path p(_inputFilename_contigs);
		//string contigFilenameBgzip = p.parent_path() + "/" ;
		//cout << contigFilenameBgzip << endl;
		//cout << p.stem() << endl;
		//exit(1);

		// = _inputFilename_contigs + "_tmpBgZip.fasta.gz";

		string contigFilenameBgzip = _tmpDir + "/__tmp_contigs_bgzip.fasta.gz";

		string command = "zcat " + _inputFilename_contigs + " | bgzip --threads " + to_string(_nbCores) + " -c > " + contigFilenameBgzip;
		Utils::executeCommand(command, _tmpDir, _logFile);

		command = "samtools faidx " + contigFilenameBgzip;
		Utils::executeCommand(command, _tmpDir, _logFile);
		//-N 2 -p 0 --secondary=yes -I 2GB -x map-hifi -c -x asm20
		//string command = "minimap2 -m 500 --dual=no -H -DP -c -I 100M -t " + to_string(_nbCores) + " " + _inputFilename_contigs + " " + _inputFilename_contigs + " > " + mapFilename;
		//Utils::executeCommand(command, _tmpDir, _logFile);

		command = "wfmash " + contigFilenameBgzip + " -t " + to_string(_nbCores) + " > " + mapFilename; //-l 5000 -p 80
		Utils::executeCommand(command, _tmpDir, _logFile);

		fs::remove(contigFilenameBgzip);
	}


	unordered_set<DbgEdge, hash_pair> _performedPairs;

	void detectDuplicatedContigs(const string& mapFilename){

		ifstream mappingFile(mapFilename);
        vector<string>* fields = new vector<string>();
        vector<string>* fields_optional = new vector<string>();

		string line;
		while (getline(mappingFile, line)) {


            GfaParser::tokenize(line, fields, '\t');

			const string& queryName = (*fields)[0];
			const string& targetName = (*fields)[5];
			//if(readName == contigName) continue;

			bool isQueryCircular = queryName[queryName.size()-1] == 'c';
			bool isTargetCircular = targetName[targetName.size()-1] == 'c';

			//string search1 ="ctg88963";
			//string search2 ="ctg89401";
			//string search1 ="ctg89096";
			//string search2 ="ctg90131";
			//string search1 ="ctg87283";
			//string search2 ="ctg90059";
			//string search1 ="ctg90131";
			//string search2 ="ctg88674";
			string search1 ="ctg88963";
			string search2 ="ctg89401";
			



			//continue;
			//if((readName == search1 || contigName == search1) && (readName == search2 || contigName == search2)){
				//_logFile << line << endl;
			//}
			//if((readName == "ctg89567" || contigName == "ctg89567")){
				//_logFile << line << endl;
			//}


			//_logFile << line << endl;

			//continue;

			//u_int64_t queryLength = stoull((*fields)[1]);
			u_int64_t targetLength = stoull((*fields)[6]);
			double queryLength = stoull((*fields)[1]);

			//if(queryLength < targetLength) continue;
			//if(targetLength > queryLength) continue;


			u_int64_t queryStart = stoull((*fields)[2]);
			u_int64_t queryEnd = stoull((*fields)[3]);
			u_int64_t targetStart = stoull((*fields)[7]);
			u_int64_t targetEnd = stoull((*fields)[8]);
			bool strand = false;
			if((*fields)[4] == "-") strand = true;

			//u_int64_t ql = queryLength;
			//u_int64_t tl = targetLength;
			//u_int64_t qs = queryStart;
			//u_int64_t qe = queryEnd;
			//u_int64_t ts = targetStart;
			//u_int64_t te = targetEnd;

			u_int64_t nbMatches = stoull((*fields)[9]);
			double alignLength = stoull((*fields)[10]);

			float identity = 0;
			//_logFile << error << endl;
			//getchar();
			//_logFile << contigIndex << " " << readIndex << endl;

			for(size_t i=12; i<fields->size(); i++){

				//_logFile << (*fields)[i] << endl;

				GfaParser::tokenize((*fields)[i], fields_optional, ':');

				if((*fields_optional)[0] == "gi"){
					identity = std::stof((*fields_optional)[2]);
				}

			}

			if((float)identity < (float)_minIdentity) continue;
			cout << identity << " " << _minIdentity << " " << (identity < _minIdentity) << " " << _minLinearLength << " " << _minCircularLength << endl;
			
			
			
			if(targetLength < queryLength){
				if(_minCircularLength > 0 && isTargetCircular && targetLength >= _minCircularLength) continue;
				if(_minLinearLength > 0 && targetLength >= _minLinearLength) continue;
				//if(isTargetLinear && targetLength >= _minCircularLength) continue;
			}
			else{
				if(_minCircularLength > 0 && isQueryCircular && queryLength >= _minCircularLength) continue;
				if(_minLinearLength > 0 && queryLength >= _minLinearLength) continue;
			}
			//continue;
			//u_int64_t ml = nbMatches;
			//u_int64_t bl = alignLength;

			//_logFile << nbMatches / alignLength << " " << alignLength << endl;
			//if(nbMatches / alignLength < 0.8) continue;
			//if(alignLength < 1000) continue;

			u_int64_t maxHang = 100;

			

			u_int64_t queryIndex = Utils::contigName_to_contigIndex(queryName);
			u_int64_t targetIndex = Utils::contigName_to_contigIndex(targetName);
			//DbgEdge edge = {targetIndex, queryIndex};
			//edge = edge.normalize();

			bool isOverlap = false;

			if(targetLength < queryLength){

				u_int64_t hangLeft = targetStart;
				u_int64_t hangRight = targetLength - targetEnd;

				if(hangLeft < maxHang){
					_duplicationBounds[targetIndex].push_back({targetStart, targetEnd});
					_logFile << "overlap left: " << targetName << " " << targetStart << " " << targetEnd << endl;
					isOverlap = true;
				}

				if(hangRight < maxHang){
					_duplicationBounds[targetIndex].push_back({targetStart, targetEnd});
					_logFile << "overlap right: " << targetName << " " << targetStart << " " << targetEnd << endl;
					isOverlap = true;
				}
			}
			else{
				u_int64_t hangLeft = queryStart;
				u_int64_t hangRight = queryLength - queryEnd;

				if(hangLeft < maxHang){
					_duplicationBounds[queryIndex].push_back({queryStart, queryEnd});
					_logFile << "overlap left: " << queryName << " " << queryStart << " " << queryEnd << endl;
					isOverlap = true;
				}

				if(hangRight < maxHang){
					_duplicationBounds[queryIndex].push_back({queryStart, queryEnd});
					_logFile << "overlap right: " << queryName << " " << queryStart << " " << queryEnd << endl;
					isOverlap = true;
				}
			}

			if(isOverlap) continue;
			//if(isOverlap){
			//	_performedPairs.insert(edge);
			//}


			//if(nbMatches / alignLength < 0.95) continue;
			//if(alignLength < 10000) continue;
			//if(targetLength > queryLength) continue;


			if(targetLength < queryLength){
				_duplicationBounds[targetIndex].push_back({targetStart, targetEnd});
				_logFile << "Add internal " << targetName << " " << targetStart << " " << targetEnd << endl;
			}
			else{
				_duplicationBounds[queryIndex].push_back({queryStart, queryEnd});
				_logFile << "Add internal " << queryName <<  " " << queryStart << " " << queryEnd << endl;
			}
		}

		mappingFile.close();

	}
	
	void dumpDereplicatedContigs(){
		
		_logFile << "Writing dereplicated contigs" << endl;

		_outputContigFile = gzopen(_outputFilename_contigs.c_str(),"wb");

		string s1 = "/home/gats/workspace/tmp/file1.fa";
		string s2 = "/home/gats/workspace/tmp/file2.fa";
		string s3 = "/home/gats/workspace/tmp/file3.fa";
		string s4 = "/home/gats/workspace/tmp/file4.fa";
		_file1 = gzopen(s1.c_str(),"wb");
		_file2 = gzopen(s2.c_str(),"wb");
		_file3 = gzopen(s3.c_str(),"wb");
		_file4 = gzopen(s4.c_str(),"wb");

		auto fp = std::bind(&Dereplicater::dumpDereplicatedContigs_read, this, std::placeholders::_1);
		ReadParser readParser(_inputFilename_contigs, true, false, _logFile);
		readParser.parse(fp);

		gzclose(_outputContigFile);
		gzclose(_file1);
		gzclose(_file2);
		gzclose(_file3);
		gzclose(_file4);
	}

	gzFile _file1;
	gzFile _file2;
	gzFile _file3;
	gzFile _file4;

	void dumpDereplicatedContigs_read(const Read& read){

		/*
		if(read._header == "ctg88963"){
			string header = ">" + read._header + '\n';
			gzwrite(_file1, (const char*)&header[0], header.size());
			string contigSequence = read._seq + '\n';
			gzwrite(_file1, (const char*)&contigSequence[0], contigSequence.size());


		}
		else if(read._header == "ctg89401"){
			string header = ">" + read._header + '\n';
			gzwrite(_file2, (const char*)&header[0], header.size());
			string contigSequence = read._seq + '\n';
			gzwrite(_file2, (const char*)&contigSequence[0], contigSequence.size());

			string seq1 = read._seq;
			seq1 = seq1.substr(0, 31159);
			string seq2 = read._seq;
			seq2 = seq2.substr(31159, seq2.size());

			
			header = ">" + read._header + "_1" + '\n';
			gzwrite(_file3, (const char*)&header[0], header.size());
			contigSequence = seq1 + '\n';
			gzwrite(_file3, (const char*)&contigSequence[0], contigSequence.size());
			
			header = ">" + read._header + "_2" + '\n';
			gzwrite(_file4, (const char*)&header[0], header.size());
			contigSequence = seq2 + '\n';
			gzwrite(_file4, (const char*)&contigSequence[0], contigSequence.size());
		}
		u_int64_t contigIndex = Utils::contigName_to_contigIndex(read._header);
		//if(_duplicatedContigIndex.find(contigIndex) != _duplicatedContigIndex.end()){
		//	_logFile << "Discard: " << read._seq.size() << endl;
		//	return;
		//}

		string seq = read._seq;


		if(_contigOverlaps.find(contigIndex) != _contigOverlaps.end()){

			const ContigOverlap& ov = _contigOverlaps[contigIndex];

			_logFile << contigIndex << " " << ov._leftPos << " " << ov._rightPos << endl;
			if(ov._leftPos != 0 && ov._rightPos != -1){
				if(ov._leftPos >= ov._rightPos){
					_logFile << "\tSlicing complete: " << read._seq.size() << endl;
					return;
				}
				else{
					//u_int64_t rightLength = seq.size()-ov._rightPos;
					seq = seq.substr(ov._leftPos, ov._rightPos-ov._leftPos);
					_logFile << "\tSlicing both: " << read._seq.size() << endl;
				}
			}
			else if(ov._leftPos != 0){
				seq = seq.substr(ov._leftPos, seq.size()-ov._leftPos);
				_logFile << "\tSlicing left: " << read._seq.size() << endl;
			}
			else{
				//u_int64_t rightLength = seq.size()-ov._rightPos;
				seq = seq.substr(ov._leftPos, ov._rightPos-ov._leftPos);
				_logFile << "\tSlicing right: " << read._seq.size() << endl;
			}
		}

		*/

		u_int64_t contigIndex = Utils::contigName_to_contigIndex(read._header);

		string seq = read._seq;

		if(_duplicationBounds.find(contigIndex) != _duplicationBounds.end()){

			//if(read._header == "ctg89567") 
			_logFile << "-----" << endl;
			vector<bool> isDuplicated(seq.size(), false);
			double nbDuplis = 0;

			for(const DbgEdge& duplicatedBound : _duplicationBounds[contigIndex]){
				//if(read._header == "ctg89567") 
				_logFile << duplicatedBound._from << " " << duplicatedBound._to << endl;
				for(size_t i=duplicatedBound._from; i<duplicatedBound._to; i++){
					isDuplicated[i] = true;
				}
			}

			for(size_t i=0; i<isDuplicated.size(); i++){
				if(isDuplicated[i]) nbDuplis += 1;
			}

			cout << read._header << " " << read._seq.size() << " " << (nbDuplis/isDuplicated.size()) << endl;
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
							_logFile << "Dump area: " << startPos << " " << startPos+length << endl;

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
					
					_logFile << "Dump area: " << startPos << " " << startPos+length << endl;
				}


			}

		}
		else{
			string header = ">" + read._header + '\n';
			gzwrite(_outputContigFile, (const char*)&header[0], header.size());
			string contigSequence = seq + '\n';
			gzwrite(_outputContigFile, (const char*)&contigSequence[0], contigSequence.size());
		}



	}


};	



#endif 


