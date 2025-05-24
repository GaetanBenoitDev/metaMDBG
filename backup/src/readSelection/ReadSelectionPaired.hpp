

#ifndef MDBG_METAG_ReadSelectionPaired
#define MDBG_METAG_ReadSelectionPaired

#include "../Commons.hpp"






class ReadSelectionPaired : public Tool{
    
public:

	string _inputFilename;
	string _inputDir;
	string _outputFilename;
	float _minimizerDensity;
    size_t _minimizerSize;
    size_t _kminmerSize;
	string _filename_readMinimizers;
	//bool _isFirstPass;
	int _nbCores;
	

    struct ReadWriter{
        u_int64_t _readIndex;
        vector<u_int64_t> _minimizers;
		vector<u_int8_t> _minimizerQualities;
    };

    struct ReadWriter_Comparator {
        bool operator()(ReadWriter const& p1, ReadWriter const& p2){
            return p1._readIndex > p2._readIndex;
        }
    };

	priority_queue<ReadWriter, vector<ReadWriter> , ReadWriter_Comparator> _readWriterQueue;
	u_int64_t _nextReadIndexWriter;

	ReadSelectionPaired(): Tool (){
	}

	void parseArgs(int argc, char* argv[]){

		args::ArgumentParser parser("readSelectionPaired", ""); //"This is a test program.", "This goes after the options."
		args::Positional<std::string> arg_outputDir(parser, "outputDir", "", args::Options::Required);
		args::Positional<std::string> arg_outputFilename(parser, "outputFilename", "", args::Options::Required);
		args::Positional<std::string> arg_inputReadFilename(parser, "inputReadFilename", "", args::Options::Required);
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
		_outputFilename = args::get(arg_outputFilename);
		_inputFilename = args::get(arg_inputReadFilename);
		_nbCores = args::get(arg_nbCores);

		string filename_parameters = _inputDir + "/parameters.gz";
		gzFile file_parameters = gzopen(filename_parameters.c_str(),"rb");
		gzread(file_parameters, (char*)&_minimizerSize, sizeof(_minimizerSize));
		gzread(file_parameters, (char*)&_kminmerSize, sizeof(_kminmerSize));
		gzread(file_parameters, (char*)&_minimizerDensity, sizeof(_minimizerDensity));
		gzclose(file_parameters);

		openLogFile(_inputDir);

		_logFile << endl;
		_logFile << "Input filename: " << _inputFilename << endl;
		_logFile << "Input dir: " << _inputDir << endl;
		_logFile << "Minimizer length: " << _minimizerSize << endl;
		_logFile << "Kminmer length: " << _kminmerSize << endl;
		_logFile << "Density: " << _minimizerDensity << endl;
		_logFile << endl;

		_filename_readMinimizers = _outputFilename; //_inputDir + "/read_data.gz";
	}

	ofstream _file_readData;


    void execute (){
		
		_nextReadIndexWriter = 0;
		_file_readData = ofstream(_filename_readMinimizers);


		ReadParser parser(_inputFilename, false, false, _logFile);

		for(size_t i=0; i<parser._filenames.size(); i+=2){
			
			string filename1 = parser._filenames[i];
			string filename2 = parser._filenames[i+1];

			_logFile << "Extracting minimizers: " << filename1 << " " << filename2 << endl;

			ReadParserParallelPaired readParser(filename1, filename2, _nbCores, _logFile);
			readParser.parse(ReadSelectionFunctor(*this, _minimizerSize, _minimizerDensity));
		}

		_file_readData.close();
		//closeLogFile();
	}



	class ReadSelectionFunctor {

		public:

		ReadSelectionPaired& _readSelection;
		size_t _minimizerSize;
		float _minimizerDensity;
		MinimizerParser* _minimizerParser;
		EncoderRLE _encoderRLE;

		ReadSelectionFunctor(ReadSelectionPaired& readSelection, size_t minimizerSize, float minimizerDensity) : _readSelection(readSelection){
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

		void operator () (const Read& read1, const Read& read2) {


			u_int64_t readIndex = read1._index;
			if(readIndex % 100000 == 0) _readSelection._logFile << readIndex << endl;

			vector<u_int64_t> minimizers1 = getMinimizers(read1);
			vector<u_int64_t> minimizers2 = getMinimizers(read2);

			minimizers1.insert( minimizers1.end(), minimizers2.begin(), minimizers2.end() );

			vector<u_int8_t> minimizerQualities;

			_readSelection.writeRead(read1, minimizers1, minimizerQualities);
		
		}

		vector<u_int64_t> getMinimizers(const Read& read){

			string rleSequence;
			vector<u_int64_t> rlePositions;
			_encoderRLE.execute(read._seq.c_str(), read._seq.size(), rleSequence, rlePositions);


			vector<u_int64_t> minimizers;
			vector<u_int64_t> minimizers_pos;
			_minimizerParser->parse(rleSequence, minimizers, minimizers_pos);

			return minimizers;
		}

	};



	void writeRead(const Read& read, const vector<u_int64_t>& minimizers, const vector<u_int8_t>& minimizerQualities){


		//#pragma omp critical(dataupdate)
		#pragma omp critical
		{
			_readWriterQueue.push({read._index, minimizers, minimizerQualities});
			//_logFile << _readWriterQueue.size() << " " << read._index << " " << _nextReadIndexWriter << endl;

			while(!_readWriterQueue.empty()){

				const ReadWriter& readWriter = _readWriterQueue.top();

				if(readWriter._readIndex == _nextReadIndexWriter){

					u_int32_t size = readWriter._minimizers.size();

					_file_readData.write((const char*)&size, sizeof(size));

					u_int8_t isCircular = CONTIG_LINEAR;
					_file_readData.write((const char*)&isCircular, sizeof(isCircular));

					_file_readData.write((const char*)&readWriter._minimizers[0], size*sizeof(u_int64_t));
					//_file_readData.write((const char*)&readWriter._minimizerQualities[0], size*sizeof(u_int8_t));


					_readWriterQueue.pop();
					_nextReadIndexWriter += 1;
				}
				else{
					break;
				}
			}
			
			//_logFile << readIndex << endl;
			//_file_readData.write((const char*)&minimizerPosOffset[0], size*sizeof(u_int16_t));
		}

	}

};	


#endif 



