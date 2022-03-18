

#ifndef MDBG_METAG_READSELECTION
#define MDBG_METAG_READSELECTION

#include "../Commons.hpp"






class ReadSelection : public Tool{
    
public:

	string _inputFilename;
	string _inputDir;
	string _outputFilename;
	float _minimizerDensity;
    size_t _minimizerSize;
    size_t _kminmerSize;
	string _filename_readMinimizers;
	bool _isFirstPass;
	int _nbCores;
	

	u_int64_t _debug_nbMinimizers;
	//unordered_map<u_int64_t, u_int64_t> _minimizerCounts;
	//unordered_map<KmerVec, KminmerData> _kminmersData;
	//gzFile _file_minimizerPos;

	ReadSelection(): Tool (){
	}

    void execute (){
		readSelection();
	}

	void parseArgs(int argc, char* argv[]){

		cxxopts::Options options("Assembly", "");
		options.add_options()
		(ARG_INPUT_FILENAME, "", cxxopts::value<string>())
		(ARG_OUTPUT_DIR, "", cxxopts::value<string>())
		(ARG_FIRST_PASS, "", cxxopts::value<bool>()->default_value("false"))
		(ARG_OUTPUT_FILENAME, "", cxxopts::value<string>())
		(ARG_NB_CORES, "", cxxopts::value<int>()->default_value(NB_CORES_DEFAULT));

		//("k,kminmerSize", "File name", cxxopts::value<std::string>())
		//("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
		//;

		if(argc <= 1){
			cout << options.help() << endl;
			exit(0);
		}

		cxxopts::ParseResult result;

		try{
			result = options.parse(argc, argv);

			_inputFilename = result[ARG_INPUT_FILENAME].as<string>();
			_inputDir = result[ARG_OUTPUT_DIR].as<string>();
			_outputFilename = result[ARG_OUTPUT_FILENAME].as<string>();
			_isFirstPass = result[ARG_FIRST_PASS].as<bool>();
			_nbCores = result[ARG_NB_CORES].as<int>();

		}
		catch (const std::exception& e){
			std::cout << options.help() << std::endl;
			std::cerr << e.what() << std::endl;
			std::exit(EXIT_FAILURE);
		}


		string filename_parameters = _inputDir + "/parameters.gz";
		gzFile file_parameters = gzopen(filename_parameters.c_str(),"rb");
		gzread(file_parameters, (char*)&_minimizerSize, sizeof(_minimizerSize));
		gzread(file_parameters, (char*)&_kminmerSize, sizeof(_kminmerSize));
		gzread(file_parameters, (char*)&_minimizerDensity, sizeof(_minimizerDensity));
		gzclose(file_parameters);

		cout << endl;
		cout << "Input filename: " << _inputFilename << endl;
		cout << "Input dir: " << _inputDir << endl;
		cout << "Minimizer length: " << _minimizerSize << endl;
		cout << "Kminmer length: " << _kminmerSize << endl;
		cout << "Density: " << _minimizerDensity << endl;
		cout << endl;

		_filename_readMinimizers = _outputFilename; //_inputDir + "/read_data.gz";
	}


    void readSelection(){
		_debug_nbMinimizers = 0;
		ofstream file_readData = ofstream(_filename_readMinimizers);
		//_file_minimizerPos = gzopen(_filename_readMinimizers.c_str(),"wb");
		
		//auto fp = std::bind(&ReadSelection::readSelection_read, this, std::placeholders::_1);
		ReadParserParallel readParser(_inputFilename, false, false, _nbCores);
		readParser.parse(ReadSelectionFunctor(_minimizerSize, _minimizerDensity, file_readData));


		file_readData.close();
		//delete _minimizerParser;
    }

	unordered_map<u_int32_t, u_int32_t> _lala;

	//void readSelection_read(){
		


	//}


	class ReadSelectionFunctor {

		public:

		size_t _minimizerSize;
		float _minimizerDensity;
		MinimizerParser* _minimizerParser;
		EncoderRLE _encoderRLE;
		ofstream& _file_readData;

		ReadSelectionFunctor(size_t minimizerSize, float minimizerDensity, ofstream& file_readData) : _file_readData(file_readData){
			_minimizerSize = minimizerSize;
			_minimizerDensity = minimizerDensity;
			_minimizerParser = new MinimizerParser(minimizerSize, minimizerDensity);
		}

		ReadSelectionFunctor(const ReadSelectionFunctor& copy) : _file_readData(copy._file_readData){
			_minimizerSize = copy._minimizerSize;
			_minimizerDensity = copy._minimizerDensity;
			_minimizerParser = new MinimizerParser(_minimizerSize, _minimizerDensity);
		}

		~ReadSelectionFunctor(){
			delete _minimizerParser;
		}

		void operator () (const Read& read) {

			u_int64_t readIndex = read._index;
			if(readIndex % 100000 == 0) cout << readIndex << endl;

			string rleSequence;
			vector<u_int64_t> rlePositions;
			_encoderRLE.execute(read._seq.c_str(), read._seq.size(), rleSequence, rlePositions);


			vector<u_int64_t> minimizers;
			vector<u_int64_t> minimizers_pos;
			_minimizerParser->parse(rleSequence, minimizers, minimizers_pos);
			//_debug_nbMinimizers += minimizers.size();

			//cout << strlen(read->seq.s) << " " << rleSequence.size() << endl;
			//cout << minimizers.size() << " " << minimizers_pos.size() << endl;
			
			//for(u_int64_t minimizer : minimizers){
			//	_minimizerCounts[minimizer] += 1;
			//}
			//for(size_t i=0; i<rlePositions.size(); i++){
			//	rlePositions[i] = i;
			//}

			//vector<KmerVec> kminmers; 
			//vector<ReadKminmer> kminmersInfo;
			//MDBG::getKminmers(_minimizerSize, _kminmerSize, minimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, readIndex);

			//for(size_t i=0; i<kminmers.size(); i++){

				//if(_kminmersData.find(kminmers[i]) == _kminmersData.end()){
				//	_kminmersData[kminmers[i]] = {0, kminmersInfo[i]._length - _minimizerSize, kminmersInfo[i]._seq_length_start, kminmersInfo[i]._seq_length_end, kminmersInfo[i]._isReversed};
				//}

				//_kminmersData[kminmers[i]]._count += 1;

			//}

			//DnaBitset* dnaBitset = new DnaBitset(string(read->seq.s));

			//u_int32_t size = dnaBitset->m_len;
			//gzwrite(_file_readData, (const char*)&size, sizeof(size));
			//gzwrite(_file_readData, dnaBitset->m_data, size);

			//u_int32_t size = strlen(read->seq.s);
			//gzwrite(_file_readData, (const char*)&size, sizeof(size));
			//gzwrite(_file_readData, read->seq.s, size);

			//cout << "----" << endl;
			
			vector<u_int16_t> minimizerPosOffset;

			if(minimizers.size() > 0){
				u_int16_t pos = minimizers_pos[0];
				//cout << pos << endl;
				minimizerPosOffset.push_back(pos);
				
				for(size_t i=1; i<minimizers_pos.size(); i++){
					u_int16_t posOffset = minimizers_pos[i] - pos;
					minimizerPosOffset.push_back(posOffset);
					pos = minimizers_pos[i];
					//cout << pos << " " << posOffset << endl;
				}
			}
	
			//#pragma omp critical(dataupdate)
			#pragma omp critical
			{
				//cout << readIndex << endl;
				u_int32_t size = minimizers.size();
				_file_readData.write((const char*)&size, sizeof(size));
				_file_readData.write((const char*)&minimizers[0], size*sizeof(u_int64_t));
				_file_readData.write((const char*)&minimizerPosOffset[0], size*sizeof(u_int16_t));
			}
			
		
		}

	};


};	


#endif 



