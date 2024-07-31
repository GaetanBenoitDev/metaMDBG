/*

./bin/metaMDBG setupCorr ~/appa/run/correction/test_deterministic/test_14/asm/tmp/ ~/appa/run/metaflye/humanGut_hac_duplex_hq/asm/assembly.fasta.gz ~/appa/run/metaflye/humanGut_hac_duplex_hq/asm/assembly_info.txt -t 32

./bin/metaMDBG setupCorr ~/appa/run/correction/nanopore_AD_10/asm/tmp/ ~/appa/run/metaflye/nanopore_AD/asm/assembly.fasta.gz ~/appa/run/metaflye/nanopore_AD/asm/assembly_info.txt -t 32

./bin/metaMDBG setupCorr ~/appa/run/correction/nanopore_AD_circ1_asm/tmp/ ~/appa/data/nanopore/subreads/circ1_metaflye/assembly.fasta ~/appa/data/nanopore/subreads/circ1_metaflye/assembly_info.txt -t 8
*/

#ifndef MDBG_METAG_SETUPCORRECTIONEVALUATION
#define MDBG_METAG_SETUPCORRECTIONEVALUATION

#include "../Commons.hpp"



class SetupCorrectionEvaluation : public Tool{
    
public:

	string _filename_exe;
	string _inputDir;
	int _nbCores;
	size_t _minimizerSize;
	size_t _kminmerSize;
	size_t _minimizerDensity;
	string _contigFilename;
	string _contigInfoFilename;
	string _inputFilename;

	vector<string>* _fields;
	vector<string>* _fields_optionnal;
	string _alignFilename;
	string _alignFilename_readsVsReads;
	unordered_map<string, u_int32_t> _contigName_to_contigCoverage;

	SetupCorrectionEvaluation(): Tool (){
	}


	void parseArgs(int argc, char* argv[]){


		_filename_exe = argv[0];

		args::ArgumentParser parser("lala", "");
		args::ValueFlag<std::string> arg_outputDir(parser, "", "outputDir", {ARG_OUTPUT_DIR2});
		//args::ValueFlag<std::string> arg_contigFilename(parser, "", "contigFilename", {"cc"});
		//args::ValueFlag<std::string> arg_contigInfoFilename(parser, "", "contigInfoFilename", {"ci"});
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
		//_contigFilename = args::get(arg_contigFilename);
		//_contigInfoFilename = args::get(arg_contigInfoFilename);
		_nbCores = args::get(arg_nbCores);
		_inputFilename = _inputDir + "/input.txt";
		_alignFilename = _inputDir + "/align.paf";
		_alignFilename_readsVsReads = _inputDir + "/align_readsVsReads.paf";
		
		string filename_parameters = _inputDir + "/parameters.gz";
		gzFile file_parameters = gzopen(filename_parameters.c_str(),"rb");
		gzread(file_parameters, (char*)&_minimizerSize, sizeof(_minimizerSize));
		gzread(file_parameters, (char*)&_kminmerSize, sizeof(_kminmerSize));
		gzread(file_parameters, (char*)&_minimizerDensity, sizeof(_minimizerDensity));
		gzclose(file_parameters);

		openLogFile(_inputDir);

		//_logFile << endl;
		//_logFile << "Input dir: " << _inputDir << endl;
		//_logFile << "Minimizer length: " << _minimizerSize << endl;
		//_logFile << "Kminmer length: " << _kminmerSize << endl;
		//_logFile << "Density: " << _minimizerDensity << endl;
		//_logFile << endl;

		//_filename_readMinimizers = _outputFilename; //_inputDir + "/read_data.gz";
		_kminmerSize = 2;

		_fields = new vector<string>();
		_fields_optionnal = new vector<string>();
	}

    void execute (){

		//cout << "Loading assembly info" << endl;
		//loadAssemblyInfoFile();

		//cout << "Converting contigs to minimizers-space" << endl;
		//convertContigToMinimizerSpace();

		//cout << "Mapping reads vs contigs" << endl;
		//mapReadsVsContigs();

		//cout << "Mapping reads vs reads" << endl;
		mapReadsVsReads();
		
		cout << "Indexing read name" << endl;
		indexReadName(_inputFilename, false);

		cout << "Indexing contig name" << endl;
		//indexReadName(_contigFilename, true);

		
		//cout << "Parsing alignments contigs" << endl;
		//PafParser pafParser(_alignFilename);
		//auto fp = std::bind(&SetupCorrectionEvaluation::parseAlignmentsGz_read, this, std::placeholders::_1);
		//pafParser.parse(fp);

		//cout << "Dumping best alignments contigs" << endl;
		//dumpBestMatches();
		

		//_readName_to_readIndex.clear();
		
		cout << "Parsing alignments reads" << endl;
		PafParser pafParser2(_alignFilename_readsVsReads);
		auto fp2 = std::bind(&SetupCorrectionEvaluation::parseAlignmentsGz_readVsRead, this, std::placeholders::_1);
		pafParser2.parse(fp2);

		cout << "Dumping read vs read matches" << endl;
		dumpReadvsReadMatches();
		

		//closeLogFile();
	}


	void loadAssemblyInfoFile(){

		std::string str;

		ifstream inputFile(_contigInfoFilename);
		std::getline(inputFile, str); //skip header

		while(!inputFile.eof()) {
			std::getline( inputFile, str);
			std::stringstream buffer(str);
			std::string temp;
			std::vector<string> fields;

			while( getline( buffer, temp, '\t') ) {
				fields.push_back(temp);
			}

			if(fields.size() == 0) break;

			//cout << fields[2] << endl;
			//cout << fields.size() << endl;
			string contigName = fields[0];
			u_int32_t contigCoverage = stoull(fields[2]);


			//cout << contigName << endl;

			_contigName_to_contigCoverage[contigName] = contigCoverage;

			
		}

		inputFile.close();
	}

	void convertContigToMinimizerSpace(){

		string contigDataFilename = _inputDir + "/contig_data_init.txt";

		if(fs::exists(contigDataFilename)){
			cout << "Skip, Contig minimizer file exists: " << contigDataFilename << endl;
			return;
		}
		

		ofstream file(_inputDir + "/inputContigs.txt");
		file << _contigFilename << endl;
		file.close();

		string command = _filename_exe + " readSelection " + _inputDir + " " + contigDataFilename + " " + _inputDir + "/inputContigs.txt" + " -t " + to_string(_nbCores);
		Utils::executeCommand(command, _inputDir);
		
	}

	void mapReadsVsContigs(){

		if(fs::exists(_alignFilename)){
			cout << "Skip, align file exists: " << _alignFilename << endl;
			return;
		}

		string command = "minimap2 -x map-ont -c -t " + to_string(_nbCores) + " " + _contigFilename + " " + Commons::inputFileToFilenames(_inputFilename);
		command += " | gzip -c - > " + _alignFilename;

		cout << command << endl;

		Utils::executeCommand(command, _inputDir);
	}

	void mapReadsVsReads(){

		if(fs::exists(_alignFilename_readsVsReads)){
			cout << "Skip, align file exists: " << _alignFilename_readsVsReads << endl;
			return;
		}


		//string command = "minimap2 -X -x map-ont -t " + to_string(_nbCores) + " " + Commons::inputFileToFilenames(_inputFilename) + " " + Commons::inputFileToFilenames(_inputFilename);
		string command = "minimap2 -x ava-ont -t " + to_string(_nbCores) + " " + Commons::inputFileToFilenames(_inputFilename) + " " + Commons::inputFileToFilenames(_inputFilename);
		command += " | gzip -c - > " + _alignFilename_readsVsReads;

		cout << command << endl;

		Utils::executeCommand(command, _inputDir);

	}

	void indexReadName(const string& filename, bool isFile){
		auto fp = std::bind(&SetupCorrectionEvaluation::indexReadName_read, this, std::placeholders::_1);
		ReadParser readParser(filename, isFile, false);
		readParser.parse(fp);
	}

	void indexReadName_read(const Read& read){
		if(read._index % 100000 == 0) cout << "\t" << read._index << endl;
		_readName_to_readIndex[Utils::shortenHeader(read._header)] = read._index;
	}

	struct ContigMatch{
		u_int32_t _contigIndex;
		bool _isReversed;
		u_int32_t _coverage;
		float _score;
		u_int32_t _contigStartPosition;
		float _identity;
	};

	unordered_map<string, u_int32_t> _readName_to_readIndex;
	unordered_map<u_int32_t, ContigMatch> _readIndex_to_bestContigMatch;

	void parseAlignmentsGz_read(const string& line){

		//cout << line << endl;
		GfaParser::tokenize(line, _fields, '\t');

		const string& readName = Utils::shortenHeader((*_fields)[0]);
		const string& contigName = Utils::shortenHeader((*_fields)[5]);

		//cout << (_readName_to_readIndex.find(readNameNanopore) != _readName_to_readIndex.end()) << endl;
		//cout << (_readName_to_readIndex.find(readNameHiFi) != _readName_to_readIndex.end()) << endl;

		//if(_readName_to_readIndex.find(readNameHiFi) == _readName_to_readIndex.end()) return;

		u_int32_t readIndex = _readName_to_readIndex[readName];
		u_int32_t contigIndex = _readName_to_readIndex[contigName];

		//if(readIndex1 == readIndex2) return;

		//u_int64_t alignLength = stoull((*_fields)[10]);

		u_int32_t readLength = stoull((*_fields)[1]);
		u_int32_t readStart = stoull((*_fields)[2]);
		u_int32_t readEnd = stoull((*_fields)[3]);
		
		//float readLength = stoull((*_fields)[1]);
		u_int32_t contigStart = stoull((*_fields)[7]);
		u_int32_t contigEnd = stoull((*_fields)[8]);

		//float score = (readEnd-readStart) / readLength;
		//u_int32_t contigLength = stoull((*_fields)[6]);
		//u_int32_t contigStart = stoull((*_fields)[7]);
		//u_int32_t contigEnd = stoull((*_fields)[8]);
		//double length = max((double)(contigEnd - contigStart), (double)(readEnd - readStart));
		//double error = 1 - min((double)(contigEnd - contigStart), (double)(readEnd - readStart)) / length;

		u_int64_t nbMatches = stoull((*_fields)[9]);
		u_int64_t alignLength = stoull((*_fields)[10]);

		float alignFraction = (readEnd-readStart) / ((float) readLength);
		float identity =  ((float)nbMatches) / alignLength;

		if(alignFraction < 0.95) return;
		if(identity < 0.95) return;

		float score = (contigEnd - contigStart) * identity;
		//if(error > 0.01) return;
		
		bool isReversed = (*_fields)[4] == "-";

		//u_int64_t score = stoull((*_fields)[10]);
		u_int32_t contigCoverage = _contigName_to_contigCoverage[contigName];

		ContigMatch match = {contigIndex, isReversed, contigCoverage, score, contigStart, identity};

		if(_readIndex_to_bestContigMatch.find(readIndex) == _readIndex_to_bestContigMatch.end()){
			_readIndex_to_bestContigMatch[readIndex] = match;
		}
		else{
			u_int64_t currentScore = _readIndex_to_bestContigMatch[readIndex]._score;
			if(score > currentScore){
				_readIndex_to_bestContigMatch[readIndex] = match;
				//cout << readIndex << " " << match._contigIndex << " " << match._score << endl;
			}
		}

	}


	unordered_map<u_int32_t, vector<ReadMatch>> _readIndex_to_bestReadMatch;

	void parseAlignmentsGz_readVsRead(const string& line){

		//cout << line << endl;
		GfaParser::tokenize(line, _fields, '\t');

		const string& queryName = Utils::shortenHeader((*_fields)[0]);
		const string& targetName = Utils::shortenHeader((*_fields)[5]);


		u_int32_t queryReadIndex = _readName_to_readIndex[queryName];
		u_int32_t referenceReadIndex = _readName_to_readIndex[targetName];

		//if(readIndex1 == readIndex2) return;

		//u_int64_t alignLength = stoull((*_fields)[10]);

		u_int32_t queryLength = stoull((*_fields)[1]);
		u_int32_t queryStart = stoull((*_fields)[2]);
		u_int32_t queryEnd = stoull((*_fields)[3]);
		
		//float readLength = stoull((*_fields)[1]);
		u_int32_t referenceLength = stoull((*_fields)[6]);
		u_int32_t referenceStart = stoull((*_fields)[7]);
		u_int32_t referenceEnd = stoull((*_fields)[8]);

		bool isReversed = (*_fields)[4] == "-";


		u_int64_t tl5 = 0;
		u_int64_t tl3 = 0;

		if(isReversed){
			tl5 = referenceLength - referenceEnd;
			tl3 = referenceStart;
		}
		else{
			tl5 = referenceStart;
			tl3 = referenceLength - referenceEnd;
		}

		u_int64_t hangLeft = tl5;
		u_int64_t hangRight = tl3;

		if(queryStart < tl5){
			hangLeft = queryStart;
		}

		if(queryLength - queryEnd < tl3){
			hangRight = queryLength - queryEnd;
		}


		//float score = (readEnd-readStart) / readLength;
		//u_int32_t contigLength = stoull((*_fields)[6]);
		//u_int32_t contigStart = stoull((*_fields)[7]);
		//u_int32_t contigEnd = stoull((*_fields)[8]);
		//double length = max((double)(contigEnd - contigStart), (double)(readEnd - readStart));
		//double error = 1 - min((double)(contigEnd - contigStart), (double)(readEnd - readStart)) / length;

		u_int32_t nbMatches = stoull((*_fields)[9]);
		u_int32_t alignLength = stoull((*_fields)[10]);

		//float alignFraction = (queryEnd-queryStart) / ((float) queryLength);
		//float identity =  ((float)nbMatches) / alignLength;

		if(alignLength < 1000) return;
		//if(identity < 0.95) return;

		//float score = (contigEnd - contigStart) * identity;
		//if(error > 0.01) return;
		

		float divergence = 0;
		//_logFile << error << endl;
		//getchar();
		//_logFile << contigIndex << " " << readIndex << endl;

		for(size_t i=12; i<_fields->size(); i++){

			//_logFile << (*fields)[i] << endl;

			GfaParser::tokenize((*_fields)[i], _fields_optionnal, ':');

			if((*_fields_optionnal)[0] == "dv"){
				divergence = std::stof((*_fields_optionnal)[2]);
			}
			else if((*_fields_optionnal)[0] == "id"){
				divergence = 1.0 - std::stof((*_fields_optionnal)[2]);
			}

		}

		//cout << divergence << endl;

		//if(divergence > 0.04) return;
		//if(isReversed){
		//	contigStart = contigLength - contigStart;
		//	contigEnd = contigLength - contigEnd;
		//}
		ReadMatch match = {queryReadIndex, isReversed, referenceStart, referenceEnd, queryStart, queryEnd, hangLeft, hangRight, nbMatches, alignLength, divergence};
		_readIndex_to_bestReadMatch[referenceReadIndex].push_back(match);

	}

	void dumpBestMatches(){

		ofstream file(_inputDir + "/readIndex_to_contigIndex.bin");

		for(const auto& it : _readIndex_to_bestContigMatch){
			
			u_int32_t readIndex = it.first;
			u_int32_t contigIndex = it.second._contigIndex;
			u_int32_t contigCoverage = it.second._coverage;
			u_int32_t contigStartPosition = it.second._contigStartPosition;
			bool isReversed = it.second._isReversed;
			float identity = it.second._identity;
			

			file.write((const char*)&readIndex, sizeof(readIndex));
			file.write((const char*)&contigIndex, sizeof(contigIndex));
			file.write((const char*)&isReversed, sizeof(isReversed));
			file.write((const char*)&contigCoverage, sizeof(contigCoverage));
			file.write((const char*)&contigStartPosition, sizeof(contigStartPosition));
			file.write((const char*)&identity, sizeof(identity));

			//cout << contigCoverage << endl;
		}

		file.close();
	}

	void dumpReadvsReadMatches(){

		ofstream file(_inputDir + "/readIndex_to_readIndex.bin");

		for(const auto& it : _readIndex_to_bestReadMatch){
			
			u_int32_t readIndex = it.first;
			const vector<ReadMatch>& matches = it.second;

			u_int32_t nbMatches = matches.size();
			file.write((const char*)&readIndex, sizeof(readIndex));
			file.write((const char*)&nbMatches, sizeof(nbMatches));
			file.write((const char*)&matches[0], matches.size()*sizeof(ReadMatch));
			/*
			for(const ReadMatch& readMatch : it.second){

				u_int32_t contigIndex = readMatch._readIndex;
				bool isReversed = readMatch._isReversed;
				

				file.write((const char*)&readIndex, sizeof(readIndex));
				file.write((const char*)&contigIndex, sizeof(contigIndex));
				file.write((const char*)&isReversed, sizeof(isReversed));
				file.write((const char*)&readMatch._alignStart, sizeof(readMatch._alignStart));
				file.write((const char*)&readMatch._alignEnd, sizeof(readMatch._alignEnd));
				file.write((const char*)&readMatch._alignNbMatches, sizeof(readMatch._alignNbMatches));
				file.write((const char*)&readMatch._alignLength, sizeof(readMatch._alignLength));

			}
			*/

		}

		file.close();
	}


};	


#endif 


