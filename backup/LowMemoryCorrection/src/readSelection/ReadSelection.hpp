

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
	//bool _isFirstPass;
	float _minReadQuality;
	int _nbCores;
	
	vector<u_int32_t> _allReadSizes;
	u_int64_t _debug_nbMinimizers;
	u_int64_t _nbLowQualityReads;

    struct ReadWriter{
        u_int64_t _readIndex;
        vector<MinimizerType> _minimizers;
		vector<u_int32_t> _minimizerPos;
        vector<u_int8_t> _minimizerDirections;
		vector<u_int8_t> _minimizerQualities;
		u_int32_t _readSize;
        //u_int32_t _prevNodeIndex;
    };

    struct ReadWriter_Comparator {
        bool operator()(ReadWriter const& p1, ReadWriter const& p2){
            return p1._readIndex > p2._readIndex;
        }
    };

	priority_queue<ReadWriter, vector<ReadWriter> , ReadWriter_Comparator> _readWriterQueue;
	u_int64_t _nextReadIndexWriter;
	ofstream _file_readStats;
	//unordered_map<u_int64_t, u_int64_t> _minimizerCounts;
	//unordered_map<KmerVec, KminmerData> _kminmersData;
	//gzFile _file_minimizerPos;

	ReadSelection(): Tool (){
	}

    void execute (){
		_nbLowQualityReads = 0;

		readSelection();

		cout << "Nb low quality reads: " << _nbLowQualityReads << endl;
		closeLogFile();
	}

	void parseArgs(int argc, char* argv[]){

		/*
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
			cerr << options.help() << endl;
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
			cerr << options.help() << std::endl;
			cerr << e.what() << std::endl;
			std::exit(EXIT_FAILURE);
		}

		*/

		args::ArgumentParser parser("readSelection", ""); //"This is a test program.", "This goes after the options."
		args::Positional<std::string> arg_outputDir(parser, "outputDir", "", args::Options::Required);
		args::Positional<std::string> arg_outputFilename(parser, "outputFilename", "", args::Options::Required);
		args::Positional<std::string> arg_inputReadFilename(parser, "inputReadFilename", "", args::Options::Required);
		//args::Positional<std::string> arg_contigs(parser, "contigs", "", args::Options::Required);
		//args::PositionalList<std::string> arg_readFilenames(parser, "reads", "Input filename(s) (separated by space)", args::Options::Required);
		//args::ValueFlag<int> arg_l(parser, "", "Minimizer length", {ARG_MINIMIZER_LENGTH2}, 13);
		//args::ValueFlag<float> arg_d(parser, "", "Minimizer density", {ARG_MINIMIZER_DENSITY2}, 0.005f);
		//args::ValueFlag<std::string> arg_contigs(parser, "", "", {ARG_INPUT_FILENAME_CONTIG}, "");
		args::ValueFlag<float> arg_minReadQuality(parser, "", "Minimum read average quality", {ARG_MIN_READ_QUALITY}, 0.0);
		args::ValueFlag<int> arg_nbCores(parser, "", "Number of cores", {ARG_NB_CORES2}, NB_CORES_DEFAULT_INT);
		//args::Flag arg_firstPass(parser, "", "Is first pass of multi-k", {ARG_FIRST_PASS});
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

		_inputDir = args::get(arg_outputDir);
		_outputFilename = args::get(arg_outputFilename);
		_inputFilename = args::get(arg_inputReadFilename);
		_minReadQuality = args::get(arg_minReadQuality);
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
		_logFile << "Min read quality: " << _minReadQuality << endl;
		_logFile << "Minimizer length: " << _minimizerSize << endl;
		_logFile << "Kminmer length: " << _kminmerSize << endl;
		_logFile << "Density: " << _minimizerDensity << endl;
		_logFile << endl;

		_filename_readMinimizers = _outputFilename; //_inputDir + "/read_data.gz";
	}

	ofstream _file_readData;

    void readSelection(){

		cout << "read selection: quality filter activated" << endl;

		_nextReadIndexWriter = 0;
		_debug_nbMinimizers = 0;
		_file_readData = ofstream(_filename_readMinimizers);
		//_file_minimizerPos = gzopen(_filename_readMinimizers.c_str(),"wb");
		
		//auto fp = std::bind(&ReadSelection::readSelection_read, this, std::placeholders::_1);
		ReadParserParallel readParser(_inputFilename, false, false, _nbCores, _logFile);
		readParser.parse(ReadSelectionFunctor(*this, _minimizerSize, _minimizerDensity));

		/*
		_logFile << _readWriterQueue.size() << endl;
		while(!_readWriterQueue.empty()){

			const ReadWriter& readWriter = _readWriterQueue.top();

			if(readWriter._readIndex == _nextReadIndexWriter){


				_logFile << "Writing read (end): " << _nextReadIndexWriter << endl;
				u_int32_t size = readWriter._minimizers.size();
				_file_readData.write((const char*)&size, sizeof(size));
				_file_readData.write((const char*)&readWriter._minimizers[0], size*sizeof(u_int64_t));

				_readWriterQueue.pop();
				_nextReadIndexWriter += 1;
			}

		}
		*/

		_file_readData.close();



	
		computeReadStats();

		//delete _minimizerParser;
    }

	void computeReadStats(){

		_file_readStats = ofstream(_inputDir + "/read_stats.txt");

		
		std::sort(_allReadSizes.begin(),  _allReadSizes.end(), std::greater<u_int32_t>());

		//u_int32_t n=_allReadSizes.size();
		//u_int32_t max=_allReadSizes[0];                 	
		//u_int32_t  sum = accumulate(_allReadSizes.begin(), _allReadSizes.end(), 0.0);
		vector<u_int64_t> readLengthCumuls;
		u_int64_t cumul = 0;

		for(size_t i=0; i<_allReadSizes.size(); i++){
			cumul += _allReadSizes[i];
			readLengthCumuls.push_back(cumul);
		}

		std::reverse(_allReadSizes.begin(), _allReadSizes.end());
		std::reverse(readLengthCumuls.begin(), readLengthCumuls.end());

		u_int32_t n50 = 0;
		u_int64_t halfsize = readLengthCumuls[0]/2;

		for(size_t i=0; i<_allReadSizes.size(); i++){
			//if(i < 10){
			//cout << _allReadSizes[i] << " " << readLengthCumuls[i] << " " << halfsize << endl;
			//}
			if(readLengthCumuls[i] < halfsize){
				n50 = _allReadSizes[i];
				break;
			}
		}

		//cout << "N50: " << n50 << endl;
		/*
		float mean = bases*1. / n;

		u_int32_t n50=0,l50=0;
		u_int32_t done=0;
		u_int32_t t50=0;
		u_int32_t ii=0;
		while(done<1){
			t50+=_allReadSizes[ii];
			if(t50 > bases*0.5) 
				done=1;
			ii++;
		}

		n50=ii;
		l50=_allReadSizes[n50];  //counting from 0
		*/
		//std::cout << std::fixed << std::setprecision(0) <<  "Bases= " << bases << " contigs= "<< n << " mean_length= " 
		//<< mean << " longest= " << max << " N50= "<< l50 << " n= " << n50   //counting from 1
		//<< std::endl;  

		//cout << n50 << " " << l50 << endl;
		

		//u_int32_t median = Utils::compute_median(_allReadSizes);
		_file_readStats.write((const char*)&n50, sizeof(n50));

		_file_readStats.close();
		
	}

	void writeRead(const Read& read, const vector<MinimizerType>& minimizers, const vector<u_int32_t>& minimizerPos, const vector<u_int8_t>& minimizerDirections, const vector<u_int8_t>& minimizerQualities){


		//#pragma omp critical(dataupdate)
		#pragma omp critical
		{
			_readWriterQueue.push({read._index, minimizers, minimizerPos, minimizerDirections, minimizerQualities, (u_int32_t) read._seq.size()});
			//_logFile << _readWriterQueue.size() << " " << read._index << " " << _nextReadIndexWriter << endl;

			while(!_readWriterQueue.empty()){

				const ReadWriter& readWriter = _readWriterQueue.top();

				if(readWriter._readIndex == _nextReadIndexWriter){

					/*
					for(auto qual : minimizerQualities){
						if(qual > 100){
							_logFile << qual << endl;
							getchar();
						}
					}
					*/
					//_logFile << "Writing read: " << _nextReadIndexWriter << endl;
					u_int32_t size = readWriter._minimizers.size();

					_file_readData.write((const char*)&size, sizeof(size));

					u_int8_t isCircular = CONTIG_LINEAR;
					_file_readData.write((const char*)&isCircular, sizeof(isCircular));

					_file_readData.write((const char*)&readWriter._minimizers[0], size*sizeof(MinimizerType));
					_file_readData.write((const char*)&readWriter._minimizerPos[0], (size+1)*sizeof(u_int32_t));
					_file_readData.write((const char*)&readWriter._minimizerDirections[0], size*sizeof(u_int8_t));
					_file_readData.write((const char*)&readWriter._minimizerQualities[0], size*sizeof(u_int8_t));

					_allReadSizes.push_back(readWriter._readSize);

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

	unordered_map<u_int32_t, u_int32_t> _lala;

	//void readSelection_read(){
		


	//}


	class ReadSelectionFunctor {

		public:

		ReadSelection& _readSelection;
		size_t _minimizerSize;
		float _minimizerDensity;
		MinimizerParser* _minimizerParser;
		EncoderRLE _encoderRLE;

		ReadSelectionFunctor(ReadSelection& readSelection, size_t minimizerSize, float minimizerDensity) : _readSelection(readSelection){
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

		void operator () (const Read& read) {

			u_int64_t readIndex = read._index;
			if(readIndex % 100000 == 0) _readSelection._logFile << readIndex << endl;

			string seq = read._seq; //.substr(0, 100);

			string rleSequence;
			vector<u_int64_t> rlePositions;
			_encoderRLE.execute(seq.c_str(), seq.size(), rleSequence, rlePositions);


			//cout << rleSequence << endl;

			vector<MinimizerType> minimizers;
			vector<u_int32_t> minimizerPos;
			vector<u_int8_t> minimizerDirections;
			_minimizerParser->parse(rleSequence, minimizers, minimizerPos, minimizerDirections);
			minimizerPos.push_back(rlePositions.size());

			long double qualitySum = 0;
			for(u_int8_t quality : read._qual){
				u_int8_t q = quality - 33;
				qualitySum += q;
			}

			float meanReadQuality = qualitySum / read._qual.size();
			//cout << meanReadQuality << endl;

			if(meanReadQuality < _readSelection._minReadQuality){

				#pragma omp atomic
				_readSelection._nbLowQualityReads += 1;

				minimizers.clear();
				minimizerPos.clear();
				minimizerDirections.clear();
				minimizerPos.push_back(rlePositions.size());
			}

			/*
			cout << seq << endl;
			cout << endl;
			cout << rleSequence << " " << rleSequence.size()  << endl;
			cout << endl;
			for(u_int64_t m : minimizers){
				cout << m << endl;
			}

			cout << endl;
			cout << endl;
			cout << endl;
			cout << endl;
			cout << endl;
			
			string rleSequenceRev;
			vector<u_int64_t> rlePositionsRev;
			string seqRev = seq;
			Utils::toReverseComplement(seqRev);
			_encoderRLE.execute(seqRev.c_str(), seqRev.size(), rleSequenceRev, rlePositionsRev);

			vector<u_int64_t> minimizersRev;
			vector<u_int64_t> minimizers_posRev;
			vector<u_int8_t> minimizerDirectionsRev;
			_minimizerParser->parse(rleSequenceRev, minimizersRev, minimizers_posRev, minimizerDirectionsRev);

			cout << seqRev << endl;
			cout << endl;
			cout << rleSequenceRev << " " << rleSequenceRev.size() << endl;
			cout << endl;
			for(u_int64_t m : minimizersRev){
				cout << m << endl;
			}

			Utils::toReverseComplement(rleSequenceRev);

			cout << (rleSequence == rleSequenceRev) << endl;
			getchar();
			*/
			//for(size_t i=0; i<minimizers.size(); i++){
			//	cout << i << ": " << rleSequence.substr(minimizers_pos[i], _minimizerSize) << " " << ((int)minimizers_direction[i]) << endl;
			//}

			//getchar();
			//size_t kmerSize = 201;
			//size_t smerSize = _minimizerSize;
			//size_t windowSize = kmerSize-30;

			//size_t windowSize = kmerSize - _minimizerSize + 1;
			//_minimizerParser->parseSyncmers(rleSequence, minimizers, minimizers_pos);
			//cout << endl << endl;
			//for(u_int64_t m : minimizers){
			//	cout << m << endl;
			//}
			
			//string rleSequenceRev = rleSequence;
			//Utils::toReverseComplement(rleSequenceRev);
			//_minimizerParser->parseSyncmers(rleSequenceRev, minimizers, minimizers_pos);
			//cout << endl << endl;
			//for(u_int64_t m : minimizers){
			//	cout << m << endl;
			//}

			//getchar();
			//cout << minimizers.size() << " " << minimizers_pos.size() << endl;
			//getchar();
			//iterateKmers(const vector<u_int16_t>& seq, size_t kmerSize, size_t smerSize, size_t windowSize)
			//_debug_nbMinimizers += minimizers.size();

			//_logFile << strlen(read->seq.s) << " " << rleSequence.size() << endl;
			//_logFile << minimizers.size() << " " << minimizers_pos.size() << endl;
			
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

			//_logFile << "----" << endl;
			
			/*
			vector<u_int16_t> minimizerPosOffset;

			if(minimizers.size() > 0){
				u_int16_t pos = minimizers_pos[0];
				//_logFile << pos << endl;
				minimizerPosOffset.push_back(pos);
				
				for(size_t i=1; i<minimizers_pos.size(); i++){
					u_int16_t posOffset = minimizers_pos[i] - pos;
					minimizerPosOffset.push_back(posOffset);
					pos = minimizers_pos[i];
					//_logFile << pos << " " << posOffset << endl;
				}
			}
			*/

			//_logFile << read._qual << endl;
			vector<u_int8_t> minimizerQualities;
			if(read._qual.empty()){
				for(size_t i=0; i<minimizers.size(); i++){
					minimizerQualities.push_back(1);
				}
			}
			else{

				

				for(size_t i=0; i<minimizers.size(); i++){

					MinimizerType minimizer = minimizers[i];
					u_int64_t pos = minimizerPos[i];

					/*
					if(minimizer == 58447242162608639){

						#pragma omp critical
						{
							cout << endl;
							cout << getMinimizerSequence(read._seq, read._qual, rlePositions[pos], rlePositions[pos+_readSelection._minimizerSize]) << endl;
							
							for(char q : getQualitySequence(read._seq, read._qual, rlePositions[pos], rlePositions[pos+_readSelection._minimizerSize])){
								cout << ((int)q)-33 << " ";
							}
							cout << endl;

							u_int8_t minQualityDebug = getMinQualityDebug(read._seq, read._qual, rlePositions[pos], rlePositions[pos+_readSelection._minimizerSize]);

							getchar();
						}


					}
					*/
					//_logFile << pos << endl;
					//double averageQuals_sum = 0;
					//double averageQuals_n = 0;
					/*
					u_int8_t minQuality = -1;

					char prevChar = '1';
					int maxQual = -1;
					
					u_int8_t minQuality = 
					//_logFile << read._qual.size() << pos+_readSelection._minimizerSize << endl;
					for(size_t i=rlePositions[pos]; i<rlePositions[pos+_readSelection._minimizerSize]; i++){
						u_int8_t quality = static_cast<u_int8_t>(read._qual[i]) - 33;

						if(read._seq[i] == prevChar){
							if(quality > maxQual){
								maxQual = quality;
							}
						}
						else{

							if(maxQual == -1){

							}

							prevChar = read._seq[i];
							maxQual = -1;


							if(quality < minQuality){
								minQuality = quality;
							}
						}

						//_logFile << i << " " << read._qual.size() << endl;
						//averageQuals_sum += quality;
						//averageQuals_n += 1;
						//_logFile << read._qual[i] << " " << to_string(quality) << endl;
						if(quality < minQuality){
							minQuality = quality;
						}
						//_logFile << quality << endl;
						//getchar();

					}
					*/

					//u_int8_t meanQuality = averageQuals_sum / averageQuals_n;


					u_int8_t minQuality = getMinQuality(read._seq, read._qual, rlePositions[pos], rlePositions[pos+_readSelection._minimizerSize]);
					minimizerQualities.push_back(minQuality);

				}

				//getchar();
			}

			/*
			if(read._index == 32){
				for(size_t i=0; i<read._qual.size(); i++){
					u_int8_t quality = static_cast<u_int8_t>(read._qual[i]) - 33;
					_logFile << to_string(quality) << " ";
				}
				_logFile << endl;
				//_logFile << minimizers.size() << endl;
				//_logFile << read._qual << endl;
				getchar();
			}
			*/

			_readSelection.writeRead(read, minimizers, minimizerPos, minimizerDirections, minimizerQualities);
		
		}

		string getMinimizerSequence(const string& seq, const string& qual, int startPos, int endPos){
			return seq.substr(startPos, endPos-startPos);
		}

		string getQualitySequence(const string& seq, const string& qual, int startPos, int endPos){
			return qual.substr(startPos, endPos-startPos);
		}

		/*
		u_int8_t getMinQualityDebug(const string& seq, const string& qual, int startPos, int endPos){

			vector<u_int8_t> qualities;

			char lastChar = '#'; //seq[0];
			u_int8_t lastQuality = 0;//qual[0];
			//u_int64_t lastPos = 0;

			for(size_t i=startPos; i<endPos; i++){
				char c = seq[i];
				cout << "\t" << i << " " << c << endl;

				if(c == lastChar){

					if(qual[i] > lastQuality){
						lastQuality = static_cast<u_int8_t>(qual[i]) - 33;
					}
					continue;
				}

				if(lastChar != '#'){

					cout << "\t- " << ((int)lastQuality) << endl;
					qualities.push_back(lastQuality);

					//runQuality = -1;
					//rleSequence += lastChar;
					//rlePositions.push_back(lastPos);
					//lastPos = i;
					//cout << lastChar << endl;
				}

				lastQuality = static_cast<u_int8_t>(qual[i]) - 33;
				lastChar = c;
			}

			qualities.push_back(lastQuality);
			//cout << lastChar << endl;
			//rleSequence += lastChar;

			u_int8_t minQuality = 200;
			for(u_int8_t quality : qualities){
				if(quality < minQuality){
					minQuality = quality;
				}
			}

			return minQuality;

		}
		*/

		u_int8_t getMinQuality(const string& seq, const string& qual, int startPos, int endPos){

			
			u_int8_t minQuality = -1;
			//cout << endl << (int) minQuality << endl;
			for(size_t i=startPos; i<endPos; i++){

				u_int8_t quality = qual[i];
				//cout << "\t" << (int) quality << endl;
				quality -= 33;

				if(quality < minQuality){
					minQuality = quality;
				}
			}

			//cout << (int) minQuality << endl;
			//getchar();
			return minQuality;

			/*
			vector<u_int8_t> qualities;

			char lastChar = '#';
			u_int64_t lastPos = 0;
			u_int8_t lastQuality = 0;

			for(size_t i=startPos; i<endPos; i++){
				//cout << i << " " << length << endl;
				char c = seq[i];

				if(c == lastChar){

					if(qual[i] > lastQuality){
						lastQuality = static_cast<u_int8_t>(qual[i]) - 33;
					}
					continue;
				}

				if(lastChar != '#'){

					qualities.push_back(lastQuality);

					//runQuality = -1;
					//rleSequence += lastChar;
					//rlePositions.push_back(lastPos);
					lastPos = i;
					//cout << lastChar << endl;
				}

				lastQuality = static_cast<u_int8_t>(qual[i]) - 33;
				lastChar = c;
			}

			qualities.push_back(lastQuality);
			//cout << lastChar << endl;
			//rleSequence += lastChar;

			u_int8_t minQuality = 20000;
			for(u_int8_t quality : qualities){
				if(quality < minQuality){
					minQuality = quality;
				}
			}

			return minQuality;
			*/
		}
	};


};	


#endif 



