/*

./bin/metaMDBG setupCorr ~/appa/run/correction/test_deterministic/test_14/asm/tmp/ ~/appa/run/metaflye/humanGut_hac_duplex_hq/asm/assembly.fasta.gz ~/appa/run/metaflye/humanGut_hac_duplex_hq/asm/assembly_info.txt -t 32

./bin/metaMDBG setupCorr ~/appa/run/correction/nanopore_AD_10/asm/tmp/ ~/appa/run/metaflye/nanopore_AD/asm/assembly.fasta.gz ~/appa/run/metaflye/nanopore_AD/asm/assembly_info.txt -t 32

./bin/metaMDBG setupCorr ~/appa/run/correction/nanopore_AD_circ1_asm/tmp/ ~/appa/data/nanopore/subreads/circ1_metaflye/assembly.fasta ~/appa/data/nanopore/subreads/circ1_metaflye/assembly_info.txt -t 8
*/

#ifndef MDBG_METAG_SETUPCORRECTIONEVALUATIONSIMULATION
#define MDBG_METAG_SETUPCORRECTIONEVALUATIONSIMULATION

#include "../Commons.hpp"



class SetupCorrectionEvaluationSimulation : public Tool{
    
public:

	string _filename_exe;
	string _inputDir;
	int _nbCores;
	string _referenceFilename;
	string _simulatedReadMetadataFilename;
	//string _inputFilename;

	ofstream _file_readData;

	float _minimizerDensity_assembly;
    size_t _minimizerSize;
    size_t _kminmerSize;		
	float _minimizerSpacingMean;
	float _kminmerLengthMean;
	float _kminmerOverlapMean;
	size_t _kminmerSizeFirst;
	size_t _kminmerSizePrev;
	size_t _kminmerSizeLast;
	size_t _meanReadLength;
	float _minimizerDensity_correction;

	SetupCorrectionEvaluationSimulation(): Tool (){
	}


	void parseArgs(int argc, char* argv[]){


		_filename_exe = argv[0];

		args::ArgumentParser parser("lala", "");
		args::ValueFlag<std::string> arg_outputDir(parser, "", "outputDir", {ARG_OUTPUT_DIR2});
		args::ValueFlag<std::string> arg_referenceFilename(parser, "", "referenceFilename", {"ref"});
		args::ValueFlag<std::string> arg_metadataFilename(parser, "", "metadataFilename", {"metadata"});
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
		_referenceFilename = args::get(arg_referenceFilename);
		_simulatedReadMetadataFilename = args::get(arg_metadataFilename);
		_nbCores = args::get(arg_nbCores);
		//_inputFilename = _inputDir + "/input.txt";
		


		string filename_parameters = _inputDir + "/parameters.gz";
		gzFile file_parameters = gzopen(filename_parameters.c_str(),"rb");	
		gzread(file_parameters, (char*)&_minimizerSize, sizeof(_minimizerSize));
		gzread(file_parameters, (char*)&_kminmerSize, sizeof(_kminmerSize));
		gzread(file_parameters, (char*)&_minimizerDensity_assembly, sizeof(_minimizerDensity_assembly));
		gzread(file_parameters, (char*)&_kminmerSizeFirst, sizeof(_kminmerSizeFirst));
		gzread(file_parameters, (char*)&_minimizerSpacingMean, sizeof(_minimizerSpacingMean));
		gzread(file_parameters, (char*)&_kminmerLengthMean, sizeof(_kminmerLengthMean));
		gzread(file_parameters, (char*)&_kminmerOverlapMean, sizeof(_kminmerOverlapMean));
		gzread(file_parameters, (char*)&_kminmerSizePrev, sizeof(_kminmerSizePrev));
		gzread(file_parameters, (char*)&_kminmerSizeLast, sizeof(_kminmerSizeLast));
		gzread(file_parameters, (char*)&_meanReadLength, sizeof(_meanReadLength));
		gzread(file_parameters, (char*)&_minimizerDensity_correction, sizeof(_minimizerDensity_correction));

		gzclose(file_parameters);

		openLogFile(_inputDir);

		cout << endl;
		cout << "Input dir: " << _inputDir << endl;
		cout << "Minimizer length: " << _minimizerSize << endl;
		cout << "Kminmer length: " << _kminmerSize << endl;
		cout << "Density: " << _minimizerDensity_correction << endl;
		cout << endl;

	}

    void execute (){

		_file_readData = ofstream(_inputDir + "/read_data_simulated.txt");

		cout << "Loading assembly info" << endl;
		loadSimulatedReadMetadataFile();


		ReadParserParallel readParser(_referenceFilename, false, false, 1, _logFile);
		readParser.parse(ReadSelectionFunctor(*this, _minimizerSize, _minimizerDensity_correction));

		_file_readData.close();
		//closeLogFile();
	}

	struct SimulatedReadMetadata{
		ReadType _referenceIndex;
		ReadType _readIndex;
		bool _isReversed;
		u_int32_t _refCoordStart;
		u_int32_t _refCoordEnd;
		float _identity;
	};

	vector<SimulatedReadMetadata> _simulatedReadMetadata;

	void loadSimulatedReadMetadataFile(){

		std::string str;

		ifstream inputFile(_simulatedReadMetadataFilename);
		//std::getline(inputFile, str); //skip header

		while(!inputFile.eof()) {
			std::getline( inputFile, str);
			std::stringstream buffer(str);
			std::string temp;
			std::vector<string> fields;

			while( getline( buffer, temp, '\t') ) {
				fields.push_back(temp);
			}

			if(fields.size() == 0) break;

			//cout << str << endl;

			//cout << fields[2] << endl;
			//cout << fields.size() << endl;
			ReadType referenceIndex = stoull(fields[0]);
			ReadType readIndex = stoull(fields[1]);
			string strand = fields[2];
			u_int32_t refCoordStart = stoull(fields[3]);
			u_int32_t refCoordEnd = stoull(fields[4]);
			string identityStr = fields[5];
			float identity = 0;
			if(refCoordStart == -1) {
				
			}
			else{
				identity = stof(identityStr);
			}
			
			//cout << referenceIndex << endl;
			//cout << readIndex << endl;
			//cout << strand << endl;
			//cout << refCoordStart << endl;
			//cout << refCoordEnd << endl;
			//cout << identity << endl;

			bool isReversed = false;
			if(strand == "-") isReversed = true;
			
			SimulatedReadMetadata metadata = {referenceIndex, readIndex, isReversed, refCoordStart, refCoordEnd, identity};


			_simulatedReadMetadata.push_back(metadata);
			//_contigName_to_contigCoverage[contigName] = contigCoverage;
			//getchar();
			
		}

		inputFile.close();

	}



	class ReadSelectionFunctor {

		public:

		SetupCorrectionEvaluationSimulation& _readSelection;
		size_t _minimizerSize;
		float _minimizerDensity;
		MinimizerParser* _minimizerParser;
		EncoderRLE _encoderRLE;

		ReadSelectionFunctor(SetupCorrectionEvaluationSimulation& readSelection, size_t minimizerSize, float minimizerDensity) : _readSelection(readSelection){
			_minimizerSize = minimizerSize;
			_minimizerDensity = minimizerDensity;
			_minimizerParser = new MinimizerParser(minimizerSize, minimizerDensity);
		}

		ReadSelectionFunctor(const ReadSelectionFunctor& copy) : _readSelection(copy._readSelection){
			_minimizerSize = copy._minimizerSize;
			_minimizerDensity = copy._minimizerDensity;
			_minimizerParser = new MinimizerParser(_minimizerSize, _minimizerDensity);
		}

		~ReadSelectionFunctor(){
			delete _minimizerParser;
		}
		//1574237525



		void operator () (const Read& read) {

			u_int64_t readIndex = read._index;
			if(readIndex % 100000 == 0) _readSelection._logFile << readIndex << endl;

			//if(readIndex != 585611) return;
			//if(readIndex < 4000) return;

			//#pragma omp critical
			//{
			//cout << readIndex << " " << read._seq.size() << endl;
			//}
			string seq = read._seq; //.substr(0, 100);

			string rleSequence;
			vector<u_int64_t> rlePositions;
			_encoderRLE.execute(seq.c_str(), seq.size(), rleSequence, rlePositions);

			vector<MinimizerType> minimizers;
			vector<u_int32_t> minimizerPos;
			vector<u_int8_t> minimizerDirections;
			_minimizerParser->parse(rleSequence, minimizers, minimizerPos, minimizerDirections);

			//if(readIndex != 0) return;
			if(minimizers.size() < 5000) return;
			cout << read._datasetIndex << " " << read._header << " " << _minimizerDensity << " " << minimizers.size() << endl;
			//getchar();


			for(const SimulatedReadMetadata& readMetadata : _readSelection._simulatedReadMetadata){

				if(readMetadata._referenceIndex != read._datasetIndex) continue;


				vector<MinimizerType> readMinimizers;

				for(size_t i=0; i<minimizers.size(); i++){

					u_int32_t pos = minimizerPos[i];

					if(pos >= readMetadata._refCoordStart && pos <= readMetadata._refCoordEnd){
						readMinimizers.push_back(minimizers[i]);
					}

					if(pos > readMetadata._refCoordEnd) break;
				}

				if(readMetadata._isReversed){
					std::reverse(readMinimizers.begin(), readMinimizers.end());
				}

				//cout << readMetadata._readIndex << " " << readMinimizers.size() << endl;

				if(minimizers.size() == 0){
					cout << "issue" << endl;
					getchar();
				}

				u_int32_t size = readMinimizers.size();

				_readSelection._file_readData.write((const char*)&size, sizeof(size));

				u_int8_t isCircular = CONTIG_LINEAR;
				_readSelection._file_readData.write((const char*)&isCircular, sizeof(isCircular));

				_readSelection._file_readData.write((const char*)&readMinimizers[0], size*sizeof(MinimizerType));

			}
			
		}


	};
};	


#endif 


