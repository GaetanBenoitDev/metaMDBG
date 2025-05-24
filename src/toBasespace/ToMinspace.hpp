

#ifndef MDBG_METAG_TOMINSPACE
#define MDBG_METAG_TOMINSPACE

#include "../Commons.hpp"



class ToMinspace : public Tool{
    
public:

	//typedef phmap::parallel_flat_hash_map<UnitigType, vector<UnitigType>, phmap::priv::hash_default_hash<UnitigType>, phmap::priv::hash_default_eq<UnitigType>, std::allocator<std::pair<UnitigType, vector<UnitigType>>>, 4, std::mutex> UnitigSequenceMap;

	//string _inputFilename;
	string _inputDir;
	float _minimizerDensity;
    size_t _minimizerSize;
    size_t _kminmerSize;
    size_t _kminmerSizeFirst;

	string _filename_inputContigs;
	//string _filename_readMinimizers;
	//string _filename_outputContigs;
	string _filename_kminmerSequences;
	string _filename_output;
	//string _filename_nodeContigs;
	//MDBG* _mdbg;
	
	float _minimizerSpacingMean;
	float _kminmerLengthMean;
	float _kminmerOverlapMean;
    size_t _kminmerSizePrev;

	u_int64_t _checksum;
	int _nbCores;

	//UnitigSequenceMap _unitigName_to_unitigSequence;
	vector<vector<MinimizerType>> _unitigName_to_minimizers;
	ofstream _outputFile;

	ToMinspace(): Tool (){

	}

	void parseArgs(int argc, char* argv[]){


		args::ArgumentParser parser("toMinspace", ""); //"This is a test program.", "This goes after the options."
		args::Positional<std::string> arg_outputDir(parser, "outputDir", "", args::Options::Required);
		args::Positional<std::string> arg_inputContigFilename(parser, "inputContigFilename", "", args::Options::Required);
		args::Positional<std::string> arg_outputContigFilename(parser, "outputContigFilename", "", args::Options::Required);
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

		_nbCores = args::get(arg_nbCores);

		
		
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

		openLogFile(_inputDir);

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
	}


    void execute (){
    

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

		//cout << "Checksum: " << _checksum << endl;
		//exit(1);

		if(_kminmerSize == _kminmerSizeFirst){
			//cout << "ToMinpsace: copy initial contigs" << endl;
			
			fs::copy(_filename_output, _filename_output + ".init", fs::copy_options::overwrite_existing);
		}

		if(_kminmerSize == _kminmerSizeFirst+1){
			//cout << "ToMinpsace: copy initial contigs" << endl;
			
			fs::copy(_filename_output, _filename_output + ".init.k" + to_string(_kminmerSizeFirst+1), fs::copy_options::overwrite_existing);
		}
	}



	void loadUnitigSequences(){

		Logger::get().debug() << "Loading unitig sequences";
		
		//UnitigNodeParser parser(_inputDir + "/unitigGraph.nodes.bin", _nbCores);
		//parser.parse(UnitigNodeParserFunctor(*this));
		
		ifstream nodeFile(_inputDir + "/unitigGraph.nodes.bin");

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

	}

	void addNode(const UnitigType& unitigName, const vector<MinimizerType>& minimizers){

        while(_unitigName_to_minimizers.size() <= unitigName){
            _unitigName_to_minimizers.push_back({});
        }
		
        _unitigName_to_minimizers[unitigName] = minimizers;


	}


	void createMinimizerContigs(){

		_nextReadIndexWriter = 0;

		_checksum = 0;

		Logger::get().debug() << "";
		Logger::get().debug() << "Creating mContigs";
		Logger::get().debug() << "";


		_outputFile = ofstream(_filename_output);

		NodePathParserParallel parser(_filename_inputContigs, _nbCores);
		parser.parse(CreateMinimizerContigsFunctor(*this));

		_outputFile.close();
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
			unitigSequenceToMinimizerSequence(nodePathObject._nodePath, minimizers);

			if(nodePathObject._isCircular && minimizers.size() > _parent._kminmerSize){
				minimizers.push_back(minimizers[_parent._kminmerSize-1]);
			}

			_parent.writeRead(nodePathObject._readIndex, minimizers, nodePathObject._isCircular);
		}

		void unitigSequenceToMinimizerSequence(const vector<UnitigType>& unitigs, vector<MinimizerType>& minimizersSequence){

			minimizersSequence.clear();


			for(size_t i=0; i<unitigs.size(); i++){

				UnitigType unitigIndex = unitigs[i];
				//cout << UnitigGraph2::unitigIndexToString(unitigIndex) << endl;

				bool isReversed;
				const UnitigType& unitigName = UnitigGraph2::unitigIndex_to_unitigName(unitigIndex, isReversed);

				vector<MinimizerType> minimizers = _parent._unitigName_to_minimizers[unitigName];

				if(isReversed){
					std::reverse(minimizers.begin(), minimizers.end());
				}

				if(i==0){

					minimizersSequence = minimizers;
				}
				else{

					minimizersSequence.insert( minimizersSequence.end(), minimizers.begin()+_parent._kminmerSize-1, minimizers.end());
					//for(size_t i=_kminmerSize-1; i<minimizers.size(); i++){
					//    minimizersSequence.push_back(minimizers[i]);
					//}
				}
			}



		}

	};


	
    struct ReadWriter{
        u_int64_t _readIndex;
		vector<MinimizerType> _contigSequence;
		bool _isCircular;
    };

    struct ReadWriter_Comparator {
        bool operator()(ReadWriter const& p1, ReadWriter const& p2){
            return p1._readIndex > p2._readIndex;
        }
    };

	priority_queue<ReadWriter, vector<ReadWriter> , ReadWriter_Comparator> _readWriterQueue;
	u_int64_t _nextReadIndexWriter;
	unordered_map<MinimizerType, u_int64_t> _minimizerCount;

	void writeRead(u_int64_t readIndex, const vector<MinimizerType>& contigSequence, u_int8_t isCircular){

		//#pragma omp critical(dataupdate)
		#pragma omp critical
		{
			_readWriterQueue.push({readIndex, contigSequence, isCircular});
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

					//_checksum += (readWriter._readIndex*contigSize);

					for(MinimizerType m : readWriter._contigSequence){
						_checksum += (m*contigSize);
						//_minimizerCount[m] += 1;
						//cout << "\t" << m << endl;
					}

					//cout << readWriter._readIndex << " " << readWriter._contigSequence.size() << " " << _checksum << endl;

					for(MinimizerType m : readWriter._contigSequence){
					//	cout << "\t" << m << endl;
					}
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



