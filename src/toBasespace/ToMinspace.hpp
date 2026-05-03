

#ifndef MDBG_METAG_TOMINSPACE
#define MDBG_METAG_TOMINSPACE

#include "../Commons.hpp"



class ToMinspace : public Tool{
    
public:


	string _inputDir;
	string _filename_inputContigs;
	string _filename_output;
	string _filename_rawUnitigs;
	int _nbCores;

	Parameters _params;
	u_int64_t _checksum;
	vector<vector<MinimizerType>> _unitigName_to_minimizers;
	vector<vector<u_int32_t>> _unitigName_to_abundances;
	ofstream _outputFile;
	string _abundanceFilename;

	/*
	//typedef phmap::parallel_flat_hash_map<UnitigType, vector<UnitigType>, phmap::priv::hash_default_hash<UnitigType>, phmap::priv::hash_default_eq<UnitigType>, std::allocator<std::pair<UnitigType, vector<UnitigType>>>, 4, std::mutex> UnitigSequenceMap;

	//string _inputFilename;
	float _minimizerDensity;
    size_t _minimizerSize;
    size_t _kminmerSize;
    size_t _kminmerSizeFirst;

	//string _filename_readMinimizers;
	//string _filename_outputContigs;
	string _filename_kminmerSequences;
	//string _filename_nodeContigs;
	//MDBG* _mdbg;
	
	float _minimizerSpacingMean;
	float _kminmerLengthMean;
	float _kminmerOverlapMean;
    size_t _kminmerSizePrev;


	//UnitigSequenceMap _unitigName_to_unitigSequence;
	ofstream _outputFileDebug;
	*/

	ToMinspace(): Tool (){

	}

	void parseArgs(int argc, char* argv[]){


		args::ArgumentParser parser("toMinspace", ""); //"This is a test program.", "This goes after the options."
		args::Positional<std::string> arg_outputDir(parser, "outputDir", "", args::Options::Required);
		args::Positional<std::string> arg_inputContigFilename(parser, "inputContigFilename", "", args::Options::Required);
		args::Positional<std::string> arg_outputContigFilename(parser, "outputContigFilename", "", args::Options::Required);
		args::Positional<std::string> arg_rawUnitigFilename(parser, "rawUnitigFilename", "", args::Options::Required);
		args::ValueFlag<std::string> arg_abundanceFilename(parser, "", "abundanceFilename", {"abundances"}, "");
		//args::Positional<std::string> arg_contigs(parser, "contigs", "", args::Options::Required);
		//args::PositionalList<std::string> arg_readFilenames(parser, "reads", "Input filename(s) (separated by space)", args::Options::Required);
		//args::ValueFlag<int> arg_l(parser, "", "Minimizer length", {ARG_MINIMIZER_LENGTH2}, 13);
		//args::ValueFlag<float> arg_d(parser, "", "Minimizer density", {ARG_MINIMIZER_DENSITY2}, 0.005f);
		//args::ValueFlag<std::string> arg_contigs(parser, "", "", {ARG_INPUT_FILENAME_CONTIG}, "");
		args::ValueFlag<int> arg_nbCores(parser, "", "Number of cores", {ARG_NB_CORES2}, NB_CORES_DEFAULT_INT);
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
		_filename_inputContigs = args::get(arg_inputContigFilename);
		_filename_output = args::get(arg_outputContigFilename);
		_filename_rawUnitigs = args::get(arg_rawUnitigFilename);
		
		_abundanceFilename = "";
		if(arg_abundanceFilename){
			_abundanceFilename = args::get(arg_abundanceFilename);
		}

		//cout << "lala: " << _abundanceFilename << endl;

		_nbCores = args::get(arg_nbCores);

		
		/*
		string filename_parameters = _inputDir + "/parameters.gz";
		gzFile file_parameters = gzopen(filename_parameters.c_str(),"rb");
		gzread(file_parameters, (char*)&_minimizerSize, sizeof(_minimizerSize));
		gzread(file_parameters, (char*)&_kminmerSize, sizeof(_kminmerSize));
		gzread(file_parameters, (char*)&_minimizerDensity, sizeof(_minimizerDensity));
		gzread(file_parameters, (char*)&_kminmerSizeFirst, sizeof(_kminmerSizeFirst));
		gzread(file_parameters, (char*)&_minimizerSpacingMean, sizeof(_minimizerSpacingMean));
		gzread(file_parameters, (char*)&_kminmerLengthMean, sizeof(_kminmerLengthMean));
		gzread(file_parameters, (char*)&_kminmerOverlapMean, sizeof(_kminmerOverlapMean));
		gzread(file_parameters, (char*)&_kminmerSizePrev, sizeof(_kminmerSizePrev));
		gzclose(file_parameters);
		*/

		_params.load(_inputDir + "/parameters.gz");

		openLogFile(_inputDir);

		/*
		Logger::get().debug() << "";
		Logger::get().debug() << "Input dir: " << _inputDir;
		//_logFile << "Output filename: " << _outputFilename << endl;
		Logger::get().debug() << "Minimizer length: " << _minimizerSize;
		Logger::get().debug() << "Kminmer length: " << _kminmerSize;
		Logger::get().debug() << "Density: " << _minimizerDensity;
		Logger::get().debug() << "";

		//_filename_outputContigs = _inputDir + "/contigs.min.gz";
		//_filename_output = _filename_inputContigs + ".min";
		//_filename_readMinimizers = _inputDir + "/readData_" + to_string(_kminmerSize) + ".txt";
		//_filename_readMinimizers = _inputDir + "/read_data.txt";
		//_filename_readMinimizers_output = _filename_readMinimizers + ".corrected.txt";
		//_filename_nodeContigs = _inputDir + "/contigs.nodepath.gz";
		*/
	}


    void execute (){
    
		//cout << "ToMinspace: on check la taille des overlap pour verifier l'histoire des palindrome dans les edges" << endl;

		loadUnitigSequences();
		createMinimizerContigs();

		Logger::get().debug() << "Contig filename: " << _filename_output;
		//closeLogFile();

		/*
		u_int64_t nbMinimizers = 0;
		u_int64_t checksum = 0;

		ofstream debugFile("/pasteur/appa/homes/gbenoit/zeus/tmp/lala_1.txt");
		for(const auto& it : _minimizerCount){
			debugFile << it.first << " " << it.second << endl;
			nbMinimizers += it.second;
			checksum += it.first * it.second;
		}
		debugFile.clear();
		cout << "Nb minimizers: " << nbMinimizers << " " << checksum << endl;
		*/

		Logger::get().debug() << "ToMinspace unitig checksum: " << _checksum;
		//exit(1);

		//if(_params._kminmerSize == _params._kminmerSizeFirst){
			//cout << "ToMinpsace: copy initial contigs" << endl;
			
		//	fs::copy(_filename_output, _filename_output + ".init", fs::copy_options::overwrite_existing);
		//}

		if(_params._kminmerSize == _params._kminmerSizeFirst+1){
			//cout << "ToMinpsace: copy initial contigs" << endl;
			
			fs::copy(_filename_output, _filename_output + ".init.k" + to_string(_params._kminmerSizeFirst+1), fs::copy_options::overwrite_existing);
		}
		
	}



	void loadUnitigSequences(){

		Logger::get().debug() << "Loading unitig sequences";
		
		//UnitigNodeParser parser(_inputDir + "/unitigGraph.nodes.bin", _nbCores);
		//parser.parse(UnitigNodeParserFunctor(*this));
		
		//ifstream nodeFile(_inputDir + "/unitigGraph.nodes.bin");

		ifstream nodeFile(_filename_rawUnitigs);

		while(true){

			u_int32_t size;
			nodeFile.read((char*)&size, sizeof(size));
			

			if(nodeFile.eof()) break;

			vector<MinimizerType> minimizers;
			minimizers.resize(size);
			nodeFile.read((char*)&minimizers[0], size * sizeof(MinimizerType));

			UnitigType unitigIndex;
			nodeFile.read((char*)&unitigIndex, sizeof(unitigIndex));

			//u_int32_t nodeName;
			//nodeFile.read((char*)&nodeName, sizeof(nodeName));

            addNode(unitigIndex/2, minimizers);

        }

        nodeFile.close();




		if(_abundanceFilename != ""){


			ifstream abundanceFile(_abundanceFilename);

			while(true){


				UnitigType unitigIndex;
				abundanceFile.read((char*)&unitigIndex, sizeof(unitigIndex));

				u_int32_t nbAbundances;
				abundanceFile.read((char*)&nbAbundances, sizeof(nbAbundances));
				

				if(abundanceFile.eof()) break;

				vector<u_int32_t> abundances;
				abundances.resize(nbAbundances);
				abundanceFile.read((char*)&abundances[0], nbAbundances * sizeof(u_int32_t));


				//_unitigs[unitigName]->_abundances = abundances;
				//_unitigs[unitigName]->_abundance = _unitigs[unitigName]->computeMedianAbundance();

				//u_int32_t unitigIndexRC = unitigIndex_toReverseDirection(unitigIndex);
				//_nodes[unitigIndexRC]->_abundances = abundances;
				//_nodes[unitigIndexRC]->_abundance = _nodes[unitigIndexRC]->computeMedianAbundance();

            	addAbundance(unitigIndex/2, abundances);
			}

			abundanceFile.close();

		}
	}

	void addNode(const UnitigType& unitigName, const vector<MinimizerType>& minimizers){

		//cout << "utg" << unitigName << " " << minimizers.size() << endl;

        while(_unitigName_to_minimizers.size() <= unitigName){
            _unitigName_to_minimizers.push_back({});
        }
		
        _unitigName_to_minimizers[unitigName] = minimizers;


	}

	void addAbundance(const UnitigType& unitigName, const vector<u_int32_t>& abundance){

		//cout << "utg" << unitigName << " " << minimizers.size() << endl;

        while(_unitigName_to_abundances.size() <= unitigName){
            _unitigName_to_abundances.push_back({});
        }
		
        _unitigName_to_abundances[unitigName] = abundance;


	}

	void createMinimizerContigs(){

		_nextReadIndexWriter = 0;

		_checksum = 0;

		Logger::get().debug() << "";
		Logger::get().debug() << "Creating mContigs";
		Logger::get().debug() << "";


		//_outputFileDebug = ofstream(_filename_output + ".debug.txt");
		_outputFile = ofstream(_filename_output);

		//cout << "single core here" << endl;
		NodePathParserParallel parser(_filename_inputContigs, _nbCores);
		parser.parse(CreateMinimizerContigsFunctor(*this));

		_outputFile.close();
		//_outputFileDebug.close();
	}


	class CreateMinimizerContigsFunctor { 

		public:

		ToMinspace& _parent;

		CreateMinimizerContigsFunctor(ToMinspace& parent) : _parent(parent){
		}

		CreateMinimizerContigsFunctor(const CreateMinimizerContigsFunctor& copy) : _parent(copy._parent){
		}

		~CreateMinimizerContigsFunctor(){
		}
		
		void operator () (const NodePath& nodePathObject) {
			
			/*
			vector<MinimizerType> contigSequence;

			long checkmSumLocal = 0;

			for(size_t i=0; i<nodePathObject._nodePath.size(); i++){
				
				u_int32_t nodeIndex = nodePathObject._nodePath[i];

				checkmSumLocal += nodeIndex;

				bool orientation;
				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex, orientation);

				if(i == 0){
					
					
					if(orientation){ //+

						vector<MinimizerType> minimizers = _graph._nodeName_entire[nodeName]._minimizers;
						for(MinimizerType minimizer : minimizers) contigSequence.push_back(minimizer);
					}
					else{
						vector<MinimizerType> minimizers = _graph._nodeName_entire[nodeName]._minimizers;
						std::reverse(minimizers.begin(), minimizers.end());
						for(MinimizerType minimizer : minimizers) contigSequence.push_back(minimizer);

					}
				}
				else {
					if(orientation){
						if(_graph._nodeName_right.find(nodeName) == _graph._nodeName_right.end()){
							 Logger::get().debug() << "omg: " << nodeName;
							 //getchar();
						}
						vector<MinimizerType> minimizers = _graph._nodeName_right[nodeName]._minimizers;
						for(MinimizerType minimizer : minimizers) contigSequence.push_back(minimizer);
					}
					else{
						if(_graph._nodeName_left.find(nodeName) == _graph._nodeName_left.end()){
							 Logger::get().debug() << "omg: " << nodeName;
							 //getchar();
						}
						vector<MinimizerType> minimizers = _graph._nodeName_left[nodeName]._minimizers;
						for(MinimizerType minimizer : minimizers) contigSequence.push_back(minimizer);
					}
				}

				//_checksum += checkmSumLocal*nodePath.size();
			}
			*/
			/*
			#pragma omp critical(createMinimizerContigs)
			{
				u_int32_t contigSize = contigSequence.size();
				//_logFile << "Write: " << contigSize << endl;
				//getchar();
				_graph._outputFile.write((const char*)&contigSize, sizeof(contigSize));
				_graph._outputFile.write((const char*)&nodePathObject._isCircular, sizeof(nodePathObject._isCircular));
				_graph._outputFile.write((const char*)&contigSequence[0], contigSize*sizeof(u_int64_t));

				//cout << contigSequence.size() << endl;

			}
			*/


			vector<MinimizerType> minimizers;
			vector<u_int32_t> abundances;
			unitigSequenceToMinimizerSequence(nodePathObject._nodePath, minimizers, abundances);

			if(nodePathObject._isCircular && minimizers.size() > _parent._params._kminmerSize){

				//cout <<  endl;
				//for(auto m : minimizers){
				//	cout << m << endl;
				//}
				//cout << "OMG" << endl;
				minimizers.push_back(minimizers[_parent._params._kminmerSize-1]);
				if(abundances.size() > 0) abundances.push_back(abundances[0]);
				//getchar();
			}

			_parent.writeRead(nodePathObject._readIndex, minimizers, nodePathObject._isCircular, abundances);
		}

		void unitigSequenceToMinimizerSequence(const vector<UnitigType>& unitigs, vector<MinimizerType>& minimizersSequence, vector<u_int32_t>& abundancesSequence){

			bool print = false;


			/*
			for(size_t i=0; i<unitigs.size(); i++){

				UnitigType unitigIndex = unitigs[i];
				bool isReversed;
				const UnitigType& unitigName = UnitigGraph2::unitigIndex_to_unitigName(unitigIndex, isReversed);
				vector<MinimizerType> minimizers = _parent._unitigName_to_minimizers[unitigName];
				if(minimizers.size() == _parent._kminmerSize-1) print = true;
			}
			*/

			//cout << "---" << endl;
			//cout << unitigs.size() << endl;
			minimizersSequence.clear();
			abundancesSequence.clear();

			int prevUnitigSize = 0;
			UnitigType prevUnitigIndex = -1;
			vector<MinimizerType> prevMinimizers;

			for(size_t i=0; i<unitigs.size(); i++){

				UnitigType unitigIndex = unitigs[i];
				
				//cout << "\t" << unitigIndex << endl;
				bool isReversed;
				const UnitigType& unitigName = UnitigGraph2::unitigIndex_to_unitigName(unitigIndex, isReversed);

				//cout << "\t" << unitigName << " " << isReversed << endl;
				//cout << "\t" << unitigName << " " << _parent._unitigName_to_minimizers.size() << " " << endl;
				//cout << "\t" << _parent._unitigName_to_minimizers[unitigName].size() << " " << endl;

				vector<MinimizerType> minimizers = _parent._unitigName_to_minimizers[unitigName];

				vector<u_int32_t> abundances;
				if(_parent._unitigName_to_abundances.size() > 0){
					abundances = _parent._unitigName_to_abundances[unitigName];
				}


				if(isReversed){
					std::reverse(minimizers.begin(), minimizers.end());
					if(abundances.size() > 0 ) std::reverse(abundances.begin(), abundances.end());
				}

				if(i==0){

					minimizersSequence = minimizers;
					abundancesSequence = abundances;

					//cout << "\tInit: " << minimizers.size() << "\t" << abundances.size() << endl;
				}
				else{

					/*
					int overlapSize = 0;

					if(minimizers.size() >= _parent._kminmerSize){
						overlapSize = _parent._kminmerSize-1;
					}
					else if(prevUnitigSize == _parent._kminmerSize-1 && minimizers.size() == _parent._kminmerSize-1){
						
						bool isReversedPrev;
						const UnitigType unitigNamePrev = UnitigGraph2::unitigIndex_to_unitigName(prevUnitigIndex, isReversedPrev);

						vector<MinimizerType> prevMinimizers = _parent._unitigName_to_minimizers[unitigNamePrev];
						if(isReversed){
							std::reverse(prevMinimizers.begin(), prevMinimizers.end());
						}

						if(prevMinimizers == minimizers){
							overlapSize = _parent._kminmerSize-1;
						}
						else{
							overlapSize = _parent._kminmerSize-1-1;
						}
					}
					else if(prevUnitigSize >= _parent._kminmerSize && minimizers.size() == _parent._kminmerSize-1){
						overlapSize = _parent._kminmerSize-1; //Skip
					}

					if(print) cout << "\tlul: " <<  overlapSize << endl;
					*/

					//cout << endl;
					//cout << UnitigGraph2::unitigIndexToString(prevUnitigIndex) << endl;
					//for(size_t i=0; i<prevMinimizers.size(); i++){
					//	cout << prevMinimizers[i] << endl;
					//}
					//cout << UnitigGraph2::unitigIndexToString(unitigIndex) << endl;
					//for(size_t i=0; i<minimizers.size(); i++){
					//	cout << minimizers[i] << endl;
					//}

					/*
					int overlapSize = 0;
					if(hasOverlapSize(prevMinimizers, minimizers, _parent._kminmerSize - 1)){
						overlapSize = _parent._kminmerSize - 1;
					}
					else if(hasOverlapSize(prevMinimizers, minimizers, _parent._kminmerSize - 2)){
						overlapSize = _parent._kminmerSize - 2;
					}
					else{
						cout << "No valid overlap" << endl;
						exit(1);
					}
					*/
					//int overlapSize = longestOverlap(prevMinimizers, minimizers);
					//if(overlapSize <= 2){
					//	cout << "Invalid overlap size: " << overlapSize << endl;
					//	exit(1);
					//}

					//cout << _parent._kminmerSize << " " << longestOverlap(prevMinimizers, minimizers) << " " << longestOverlap2(prevMinimizers, minimizers) << endl;
					
					int overlapSize = longestOverlap2(prevMinimizers, minimizers);
					
					//if(longestOverlap2(prevMinimizers, minimizers) < _parent._kminmerSize-1){
					//	cout << "derp" << endl;
					//	getchar();
					//}

					//if(overlapSize < _parent._kminmerSize){
					//	cout << "Invalid overlap size 2: " << overlapSize << endl;
					//	exit(1);
					//}
					//cout << overlapSize << endl;
					//if(overlapSize != _parent._kminmerSize-1 && overlapSize != _parent._kminmerSize) getchar();
					
					minimizersSequence.insert( minimizersSequence.end(), minimizers.begin()+overlapSize, minimizers.end());
					
					if(abundances.size() > 0) abundancesSequence.insert( abundancesSequence.end(), abundances.begin(), abundances.end());
					
					//cout << "\tAdd: " << minimizers.size() << "\t" << abundances.size() << endl;
					//for(size_t i=_kminmerSize-1; i<minimizers.size(); i++){
					//    minimizersSequence.push_back(minimizers[i]);
					//}
					/*
					int overlapSizeX = longestOverlap(prevMinimizers, minimizers);
					if(overlapSizeX < _parent._kminmerSize - 1){
						cout << endl;
						cout << "Issue ToMinspace overlap size: " << overlapSizeX << endl;
						cout << UnitigGraph2::unitigIndexToString(prevUnitigIndex) << endl;
						for(size_t i=0; i<prevMinimizers.size(); i++){
							cout << prevMinimizers[i] << endl;
						}
						cout << UnitigGraph2::unitigIndexToString(unitigIndex) << endl;
						for(size_t i=0; i<minimizers.size(); i++){
							cout << minimizers[i] << endl;
						}
					}
					*/

				}



				if(print) cout << "\t" << minimizers.size() << endl;
				if(print){
					for(MinimizerType m : minimizers){
						cout << "\t\t" << m << endl;
					}
				}

				prevUnitigSize = minimizers.size();
				prevUnitigIndex = unitigIndex;
				prevMinimizers = minimizers;
			}

			if(print) cout << "\t\t\t" << minimizersSequence.size() << endl;
			if(print){
				for(MinimizerType m : minimizersSequence){
					cout << "\t\t\t" << m << endl;
				}
			}


		}


		int longestOverlap(const vector<MinimizerType>& list1, const vector<MinimizerType>& list2) {

			
			size_t n = list1.size();
			size_t m = list2.size();
			size_t maxk = min(n, m);

			for (size_t k = maxk; k >= 1; --k) {
				bool ok = true;
				for (size_t i = 0; i < k; ++i) {
					if (list1[n - k + i] != list2[i]) {
						ok = false;
						break;
					}
				}
				if (ok) return k;
			}
			return 0;
			
		}
		
		int longestOverlap2(const vector<MinimizerType>& list1, const vector<MinimizerType>& list2) {

			if(list1.size() == _parent._params._kminmerSize && list2.size() == _parent._params._kminmerSize){
				if(list1 == list2) return _parent._params._kminmerSize;
			}

			return _parent._params._kminmerSize - 1;
		
		}


	};



	
    struct ReadWriter{
        u_int64_t _readIndex;
		vector<MinimizerType> _contigSequence;
		bool _isCircular;
		vector<u_int32_t> _abundances;
    };

    struct ReadWriter_Comparator {
        bool operator()(ReadWriter const& p1, ReadWriter const& p2){
            return p1._readIndex > p2._readIndex;
        }
    };

	priority_queue<ReadWriter, vector<ReadWriter> , ReadWriter_Comparator> _readWriterQueue;
	u_int64_t _nextReadIndexWriter;
	unordered_map<MinimizerType, u_int64_t> _minimizerCount;

	void writeRead(u_int64_t readIndex, const vector<MinimizerType>& contigSequence, u_int8_t isCircular, const vector<u_int32_t>& abundances){

		//#pragma omp critical(dataupdate)
		#pragma omp critical
		{
			_readWriterQueue.push({readIndex, contigSequence, isCircular, abundances});
			//_logFile << _readWriterQueue.size() << " " << read._index << " " << _nextReadIndexWriter << endl;

			while(!_readWriterQueue.empty()){

				const ReadWriter& readWriter = _readWriterQueue.top();

				if(readWriter._readIndex == _nextReadIndexWriter){

					
					//if(readWriter._isCircular){
					//	cout << "Discard circular" << endl;
					//	_readWriterQueue.pop();
					//	_nextReadIndexWriter += 1;
					//	continue;
						//for(MinimizerType m : readWriter._contigSequence){
						//	cout << "\t" << m << endl;
						//}

						//getchar();
					//}
					

					//if(!readWriter._isCircular){
					u_int32_t contigSize = readWriter._contigSequence.size();
					_outputFile.write((const char*)&contigSize, sizeof(contigSize));
					_outputFile.write((const char*)&readWriter._isCircular, sizeof(readWriter._isCircular));
					_outputFile.write((const char*)&readWriter._contigSequence[0], contigSize*sizeof(MinimizerType));

					//cout << readWriter._contigSequence.size() << " " << readWriter._abundances.size() << endl;
					if(abundances.size() > 0){
						u_int32_t abundanceSize = readWriter._abundances.size();
						_outputFile.write((const char*)&abundanceSize, sizeof(abundanceSize));
						_outputFile.write((const char*)&readWriter._abundances[0], abundanceSize*sizeof(u_int32_t));

						//cout << "Write unitig: " << readWriter._contigSequence.size() << "\t" << readWriter._abundances.size() << endl;
					}
					//_checksum += (readWriter._readIndex*contigSize);

					for(MinimizerType m : readWriter._contigSequence){
						_checksum += (m*contigSize);
						//_minimizerCount[m] += 1;
						//cout << "\t" << m << endl;
					}


					//_outputFileDebug << readWriter._readIndex << " " << readWriter._contigSequence.size() << endl;

					//for(MinimizerType m : readWriter._contigSequence){
					//	_outputFileDebug << m << " ";
					//}
					//_outputFileDebug << endl;
					//getchar();

					//cout << readWriter._readIndex << " " << readWriter._contigSequence.size() << " " << _checksum << endl;
					//}
					/*
					//if(readWriter._isCircular){

					//	cout << "circ" << endl;

					//}

					if(!readWriter._isCircular){
						



					}
					else{
						cout << "circ" << endl;
						cout << readWriter._readIndex << " " << readWriter._contigSequence.size() << " " << _checksum << endl;

						for(MinimizerType m : readWriter._contigSequence){
							cout << "\t" << m << endl;
						}
					}
					*/
					//cout << readWriter._readIndex << " " << readWriter._contigSequence.size() << " " << _checksum << endl;
					//getchar();
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
	

	/*
	struct Contig{
		u_int64_t _readIndex;
		vector<u_int32_t> _nodepath;
		//vector<u_int32_t> _nodepath_sorted;
	};

	vector<Contig> _contigs;



	phmap::parallel_flat_hash_map<u_int32_t, vector<u_int32_t>> _nodeName_to_contigs;

	phmap::parallel_flat_hash_set<u_int32_t> _invalidContigIndex;

	void checkDuplication(){

		cout << "Removing duplication" << endl;
		
		//KminmerParserParallel parser(_filename_output, _minimizerSize, _kminmerSizeFirst, false, false);
		//auto fp = std::bind(&ToMinspace::loadContigs_min_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
		//parser.parseMinspace(fp);

		KminmerParserParallel parser(_filename_output, _minimizerSize, _kminmerSizeFirst, false, false, _nbCores);
		parser.parse(LoadContigKminmerFunctor(*this));

		for(auto& it : _nodeName_to_contigs){
			std::sort(it.second.begin(), it.second.end());
		}

		#pragma omp parallel num_threads(_nbCores)
		{

			#pragma omp for
			for(size_t i=0; i<_contigs.size(); i++){

				const Contig& contig = _contigs[i];

				u_int64_t length = _kminmerLengthMean + ((contig._nodepath.size()-1) * (_kminmerLengthMean-_kminmerOverlapMean));
				if(length > 200000) continue; //todo use length


				bool sharing = true;

				vector<u_int32_t> sharingContigIndexes = _nodeName_to_contigs[contig._nodepath[0]];

				for(size_t j=0; j<contig._nodepath.size(); j++){

					u_int32_t nodeName = contig._nodepath[j];

					vector<u_int32_t> sharingContigIndexesTmp = _nodeName_to_contigs[nodeName];
					//cout << sharingContigIndexesTmp.size() << endl;
					vector<u_int32_t> intersection;
					//std::sort(v1.begin(), v1.end());
					//std::sort(v2.begin(), v2.end());

					std::set_intersection(sharingContigIndexes.begin(),sharingContigIndexes.end(),
										sharingContigIndexesTmp.begin(),sharingContigIndexesTmp.end(),
										back_inserter(intersection));

					if(intersection.size() <= 1){
						sharing = false;
						break;
					}

					sharingContigIndexes = intersection;
				}

				if(sharing){
					cout << "TTT duplicate: " << contig._nodepath.size() << endl;
					
					#pragma omp critical(checkDuplication)
					{
						_invalidContigIndex.insert(contig._readIndex);
					}
				}

			}
		}
		
	
		rewriteContigs();
	}
	

	class LoadContigKminmerFunctor {

		public:

		ToMinspace& _graph;

		LoadContigKminmerFunctor(ToMinspace& graph) : _graph(graph){

		}

		LoadContigKminmerFunctor(const LoadContigKminmerFunctor& copy) : _graph(copy._graph){

		}

		~LoadContigKminmerFunctor(){
			//delete _minimizerParser;
		}



		void operator () (const KminmerList& kminmerList) {


			u_int64_t readIndex = kminmerList._readIndex;
			const vector<u_int64_t>& readMinimizers = kminmerList._readMinimizers;
			//const vector<KmerVec>& kminmers = kminmerList._kminmers;
			const vector<ReadKminmerComplete>& kminmersInfos = kminmerList._kminmersInfo;

			vector<u_int32_t> nodepath;

			for(size_t i=0; i<kminmersInfos.size(); i++){
				
				const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];

				KmerVec vec = kminmerInfo._vec;
				
				if(_graph._mdbg->_dbg_nodes.find(vec) == _graph._mdbg->_dbg_nodes.end()){
					_graph._logFile << "Not found kminmer" << endl;
					//getchar();
					continue;
				}



				u_int32_t nodeName = _graph._mdbg->_dbg_nodes[vec]._index;
				nodepath.push_back(nodeName);
				
				#pragma omp critical
				{
					
					vector<u_int32_t>& contigIndexes = _graph._nodeName_to_contigs[nodeName]; 
					if(std::find(contigIndexes.begin(), contigIndexes.end(), readIndex) == contigIndexes.end()){
						contigIndexes.push_back(readIndex);
					}
				}
				

			}

			//_nbContigs += 1;

			#pragma omp critical
			{
				//vector<u_int32_t> nodepath_sorted = nodepath;
				//std::sort(nodepath_sorted.begin(), nodepath_sorted.end());
				_graph._contigs.push_back({readIndex, nodepath});
				//cout << "load contig: " << nodepath.size() << endl;
			}
			
		}
	};
	
	ofstream _tmpOutputFile;

	void rewriteContigs(){

		const string& tmpFilename = _filename_output + ".tmp";
		_tmpOutputFile = ofstream(tmpFilename);

		KminmerParser parser(_filename_output, _minimizerSize, _kminmerSizeFirst, false, false);
		auto fp = std::bind(&ToMinspace::rewriteContigs_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
		parser.parseSequences(fp);

		_tmpOutputFile.close();
		
		
		fs::remove(_inputDir + "/unitig_data.txt");
		fs::rename(tmpFilename, _inputDir + "/unitig_data.txt");
	}

	void rewriteContigs_read(const vector<u_int64_t>& readMinimizers, u_int8_t isCircular, u_int64_t readIndex){

		if(_invalidContigIndex.find(readIndex) != _invalidContigIndex.end()) return;

		
		u_int32_t contigSize = readMinimizers.size();
		//_logFile << "Write: " << contigSize << endl;
		//getchar();
		_tmpOutputFile.write((const char*)&contigSize, sizeof(contigSize));
		_tmpOutputFile.write((const char*)&isCircular, sizeof(isCircular));
		_tmpOutputFile.write((const char*)&readMinimizers[0], contigSize*sizeof(u_int64_t));

	}

	*/

	

	

};	


#endif 



