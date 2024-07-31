

#ifndef MDBG_METAG_TOMINSPACE
#define MDBG_METAG_TOMINSPACE

#include "../Commons.hpp"

/*
struct KminmerSequence{
	u_int32_t _readIndex;
	vector<u_int64_t> _minimizers;
};*/

struct KminmerSequenceEntire{
	vector<MinimizerType> _minimizers;
	bool _loaded;
};

struct KminmerSequenceLeftRight{
	vector<MinimizerType> _minimizer;
	bool _loaded;
};

typedef phmap::parallel_flat_hash_map<u_int32_t, KminmerSequenceEntire, phmap::priv::hash_default_hash<u_int32_t>, phmap::priv::hash_default_eq<u_int32_t>, std::allocator<std::pair<u_int32_t, KminmerSequenceEntire>>, 4, std::mutex> KminmerSequenceMap;




class ToMinspace : public Tool{
    
public:

	//string _inputFilename;
	string _inputDir;
	float _minimizerDensity;
    size_t _minimizerSize;
    size_t _kminmerSize;
    size_t _kminmerSizeFirst;

	string _filename_inputContigs;
	string _filename_readMinimizers;
	//string _filename_outputContigs;
	string _filename_kminmerSequences;
	string _filename_output;
	//string _filename_nodeContigs;
	MDBG* _mdbg;
	
	float _minimizerSpacingMean;
	float _kminmerLengthMean;
	float _kminmerOverlapMean;
    size_t _kminmerSizePrev;

	KminmerSequenceMap _nodeName_entire;
	//unordered_map<u_int32_t, KminmerSequenceEntire> _nodeName_entireRC;
	KminmerSequenceMap _nodeName_right;
	KminmerSequenceMap _nodeName_left;

	u_int64_t _checksum;
	int _nbCores;

	//unordered_map<u_int32_t, vector<u_int64_t>> _debugMinimizers;
	//unordered_map<u_int32_t, bool> _debugMinimizers_isReversed;
	//unordered_map<u_int32_t, bool> _debugMinimizers;

	ToMinspace(): Tool (){

		/*
		//getParser()->push_back (new OptionOneParam (STR_INPUT, "input file", true));
		getParser()->push_back (new OptionOneParam (STR_INPUT_DIR, "input dir", true));
		//getParser()->push_back (new OptionOneParam (STR_OUTPUT, "output contig filename in basespace", true));
		*/
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
		//_filename_inputContigs = args::get(arg_contigs);

		_nbCores = args::get(arg_nbCores);

		/*
		cxxopts::Options options("ToMinspace", "");
		options.add_options()
		//(ARG_INPUT_FILENAME, "", cxxopts::value<string>())
		(ARG_OUTPUT_DIR, "", cxxopts::value<string>())
		(ARG_OUTPUT_FILENAME, "", cxxopts::value<string>())
		(ARG_INPUT_FILENAME_CONTIG, "", cxxopts::value<string>()->default_value(""));

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

			//_inputFilename = result[ARG_INPUT_FILENAME].as<string>();
			_inputDir = result[ARG_OUTPUT_DIR].as<string>();
			_filename_inputContigs = result[ARG_INPUT_FILENAME_CONTIG].as<string>();
			_filename_output = result[ARG_OUTPUT_FILENAME].as<string>();
			
		}
		catch (const std::exception& e){
			std::_logFile << options.help() << std::endl;
			std::cerr << e.what() << std::endl;
			std::exit(EXIT_FAILURE);
		}
		*/
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

		_logFile << endl;
		_logFile << "Input dir: " << _inputDir << endl;
		//_logFile << "Output filename: " << _outputFilename << endl;
		_logFile << "Minimizer length: " << _minimizerSize << endl;
		_logFile << "Kminmer length: " << _kminmerSize << endl;
		_logFile << "Density: " << _minimizerDensity << endl;
		_logFile << endl;

		//_filename_outputContigs = _inputDir + "/contigs.min.gz";
		//_filename_output = _filename_inputContigs + ".min";
		//_filename_readMinimizers = _inputDir + "/readData_" + to_string(_kminmerSize) + ".txt";
		_filename_readMinimizers = _inputDir + "/read_data.txt";
		//_filename_readMinimizers_output = _filename_readMinimizers + ".corrected.txt";
		//_filename_nodeContigs = _inputDir + "/contigs.nodepath.gz";
	}


    void execute (){
    
		loadContigs();

		_logFile << "Loading mdbg" << endl;
		string mdbg_filename = _inputDir + "/mdbg_nodes.gz";
		//_mdbg = new MDBG(_kminmerSize);
		//_mdbg->load(mdbg_filename, false);
		//_logFile << "MDBG nodes: " << _mdbg->_dbg_nodes.size() << endl;

		extractKminmerSequences();
		//_logFile << "delete mdbg a remetre !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
		//delete _mdbg;

		createMinimizerContigs();
		//dereplicateContigs();
		//checkDuplication();

		_logFile << endl << "Contig filename: " << _filename_output << endl;
		closeLogFile();
	}



	void loadContigs(){

		_logFile << "Loading contigs" << endl;
		
		NodePathParserParallel parser2(_filename_inputContigs, _nbCores);
		parser2.parse(LoadContigFunctor(*this));
		/*
		ifstream contigFile(_filename_inputContigs);
		u_int64_t nbContigs = 0;
		
		while(true){

			vector<u_int32_t> nodePath;
			//vector<u_int64_t> supportingReads;
			u_int64_t size;
			contigFile.read((char*)&size, sizeof(size));
			

			if(contigFile.eof()) break;

			u_int8_t isCircular;
			contigFile.read((char*)&isCircular, sizeof(isCircular));

			nodePath.resize(size);
			//supportingReads.resize(size);
			contigFile.read((char*)&nodePath[0], size * sizeof(u_int32_t));
			//gzread(contigFile, (char*)&supportingReads[0], size * sizeof(u_int64_t));

			//for(u_int32_t nodeIndex : nodePath){
			for(size_t i=0; i<nodePath.size(); i++){
				u_int32_t nodeIndex = nodePath[i];
				bool orientation;
				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex, orientation);
				//u_int32_t readIndex = supportingReads[i];


				if(i == 0){
					_nodeName_entire[nodeName] = {{}, false};
					//_logFile << nodeName << endl;
					//if(orientation){ //+
					//	_nodeName_entire[nodeName] = {{}, false};
					//}
					//else{ //-
					//	_nodeName_entireRC[nodeName] = {{}, false};
					//}
				}
				else {
					//if(i==1) _logFile << orientation << endl;

					if(orientation){ //+
						_nodeName_right[nodeName] = {{}, false};
					}
					else{ //-
						_nodeName_left[nodeName] = {{}, false};
					}
				}

			}

			nbContigs += 1;
			//_logFile << nodePath.size() << endl;
		}

		contigFile.close();

		_logFile << _nodeName_entire.size() << " " << _nodeName_left.size() << " " << _nodeName_right.size() << endl;
		_logFile << "Nb contigs: " << nbContigs << endl;
		*/
		
	}

	class LoadContigFunctor { 

		public:

		ToMinspace& _graph;

		LoadContigFunctor(ToMinspace& graph) : _graph(graph){
		}

		LoadContigFunctor(const LoadContigFunctor& copy) : _graph(copy._graph){
		}

		~LoadContigFunctor(){
		}
		
		void operator () (const NodePath& nodePathObject) {
			
			for(size_t i=0; i<nodePathObject._nodePath.size(); i++){
				u_int32_t nodeIndex = nodePathObject._nodePath[i];
				bool orientation;
				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex, orientation);
				//u_int32_t readIndex = supportingReads[i];


				if(i == 0){


					_graph._nodeName_entire.lazy_emplace_l(nodeName, 
					[](KminmerSequenceMap::value_type& v) { // key exist
					},           
					[&nodeName](const KminmerSequenceMap::constructor& ctor) { // key inserted
						KminmerSequenceEntire seq = {{}, false};
						ctor(nodeName, seq); 

					});

					//_graph._nodeName_entire[nodeName] = {{}, false};
					//_logFile << nodeName << endl;
					//if(orientation){ //+
					//	_nodeName_entire[nodeName] = {{}, false};
					//}
					//else{ //-
					//	_nodeName_entireRC[nodeName] = {{}, false};
					//}
				}
				else {
					//if(i==1) _logFile << orientation << endl;

					if(orientation){ //+
						_graph._nodeName_right.lazy_emplace_l(nodeName, 
						[](KminmerSequenceMap::value_type& v) { // key exist
						},           
						[&nodeName](const KminmerSequenceMap::constructor& ctor) { // key inserted
							KminmerSequenceEntire seq = {{}, false};
							ctor(nodeName, seq); 

						});
						//_graph._nodeName_right[nodeName] = {{}, false};
					}
					else{ //-
						_graph._nodeName_left.lazy_emplace_l(nodeName, 
						[](KminmerSequenceMap::value_type& v) { // key exist
						},           
						[&nodeName](const KminmerSequenceMap::constructor& ctor) { // key inserted
							KminmerSequenceEntire seq = {{}, false};
							ctor(nodeName, seq); 

						});
						//_graph._nodeName_left[nodeName] = {{}, false};
					}
				}

			}

		}

	};

	/*
	void extractKminmerSequences_read(const vector<u_int64_t>& readMinimizers, const vector<ReadKminmerComplete>& kminmersInfos, u_int64_t readIndex){

		//_logFile << "----" << endl;
		//_logFile << readIndex << " " << kminmersInfos.size() << endl;

		for(size_t i=0; i<kminmersInfos.size(); i++){
			
			const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];

			KmerVec vec = kminmerInfo._vec;
			//for(KmerVec& vec : kminmers){
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()) continue;

			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;

			//_logFile << nodeName << endl;
			//_logFile << nodeName << " " << kminmerInfo._length << endl;
			
			vector<u_int64_t> minimizerSeq;

			//_logFile << kminmerInfo._read_pos_start << " " << kminmerInfo._read_pos_end << endl;

			for(size_t i=kminmerInfo._read_pos_start; i<=kminmerInfo._read_pos_end; i++){
				minimizerSeq.push_back(readMinimizers[i]);
			}
			//if(kminmerInfo._length != 3){
			//	getchar();
			//}


			if(kminmerInfo._isReversed){
				std::reverse(minimizerSeq.begin(), minimizerSeq.end());
			}



			//_debugMinimizers_isReversed[nodeName] = kminmerInfo._isReversed;
			//_debugMinimizers[nodeName] = vec._kmers;
			
			auto found = _nodeName_entire.find(nodeName);
			if(found != _nodeName_entire.end() && !found->second._loaded){

				//_logFile << "entire" << endl;
				
				u_int32_t startPosition = 0; //kminmerInfo._read_pos_start;
				u_int32_t len = kminmerInfo._read_pos_end - kminmerInfo._read_pos_start; //startPosition + kminmerInfo._read_pos_end - kminmerInfo._read_pos_start;
				

				//_logFile << startPosition << " " << len << " "  << minimizerSeq.size() << endl;
				vector<u_int64_t> minimizerSeq_sub; 
				for(size_t i=startPosition; i <= len; i++){
					//_logFile << i << " " << " " << minimizerSeq.size() << endl;
					minimizerSeq_sub.push_back(minimizerSeq[i]);
					//_logFile << minimizerSeq[i] << " ";
				}
			
				//_logFile << startPosition << " " << len << " " << minimizerSeq_sub.size() << endl;

				//if(kminmerInfo._isReversed){
				//	std::reverse(vec._kmers.begin(), vec._kmers.end());
				//}
				if(kminmerInfo._isReversed){
					//std::reverse(minimizerSeq.begin(), minimizerSeq.end());
				}
				_nodeName_entire[nodeName] = {minimizerSeq, true};

				//getchar();
			}
			



			auto found2 = _nodeName_left.find(nodeName);
			if(_nodeName_left.find(nodeName) != _nodeName_left.end() && !found2->second._loaded){
				
				//_logFile << "left" << endl;
				//std::reverse(vec._kmers.begin(), vec._kmers.end());
				//if(12424 == nodeName){
				//	_logFile << "left: " <<  vec._kmers[0] << " " << vec._kmers[vec._kmers.size()-1] << endl;
				//}


				vector<u_int64_t> minimizerSeq_sub;

				u_int32_t startPosition = 0; //kminmerInfo._read_pos_start;
				u_int32_t len = startPosition+kminmerInfo._seq_length_start;

				//_logFile << "left: " << startPosition << " " << len << endl;
				
				for(size_t i=startPosition; i<len; i++){
					//_logFile << i << " " << " " << minimizerSeq.size() << endl;
					minimizerSeq_sub.push_back(minimizerSeq[i]);
					//_logFile << minimizerSeq[i] << " ";
				}
				//_logFile << endl;

				_nodeName_left[nodeName] = {minimizerSeq_sub, true};
				
				//getchar();
			}
			

			found2 = _nodeName_right.find(nodeName);
			if(_nodeName_right.find(nodeName) != _nodeName_right.end() && !found2->second._loaded){

				//_logFile << "right" << endl;

			
				vector<u_int64_t> minimizerSeq_sub;

				u_int32_t startPosition = minimizerSeq.size() - kminmerInfo._seq_length_end;  //kminmerInfo._position_of_second_to_last_minimizer_seq;
				u_int32_t len = startPosition+kminmerInfo._seq_length_end; //kminmerInfo._read_pos_end

				//_logFile << "right: " << startPosition << " " << len << endl;

				for(size_t i=startPosition; i<len; i++){
					//_logFile << i << " " << " " << minimizerSeq.size() << endl;

					minimizerSeq_sub.push_back(minimizerSeq[i]);
					//_logFile << minimizerSeq[i] << " ";
				}
				//_logFile << endl;

				_nodeName_right[nodeName] = {minimizerSeq_sub, true};
				//getchar();

			}

			//if(12424 == nodeName){
				//exit(1);
			//}

		}

		//_logFile << "done" << endl;
	}
	*/
	
	void extractKminmerSequences (){

		u_int16_t size = _kminmerSize;
		_logFile << "Extracting kminmer sequences" << endl;
		
		
		ifstream kminmerFile(_inputDir + "/kminmerData_min.txt");
			
		#pragma omp parallel num_threads(_nbCores)
		{

			bool isEOF = false;
			vector<MinimizerType> minimizerSeq;
			u_int32_t nodeName;
			u_int32_t length;
			u_int32_t lengthStart;
			u_int32_t lengthEnd;
			u_int32_t abundance;
			//u_int32_t quality;
			
			while (true) {


				#pragma omp critical
				{
					
					minimizerSeq.resize(size);
					kminmerFile.read((char*)&minimizerSeq[0], size*sizeof(MinimizerType));

					//int result = kseq_read(seq);
					isEOF = kminmerFile.eof();

					if(!isEOF){

						kminmerFile.read((char*)&nodeName, sizeof(nodeName));
						kminmerFile.read((char*)&abundance, sizeof(abundance));
						//kminmerFile.read((char*)&quality, sizeof(quality));
						//kminmerFile.read((char*)&length, sizeof(length));
						//kminmerFile.read((char*)&lengthStart, sizeof(lengthStart));
						//kminmerFile.read((char*)&lengthEnd, sizeof(lengthEnd));

						length = _kminmerSize-1;
						lengthStart = 1;
						lengthEnd = 1;

					}

				}

				if(isEOF) break;


				//bool isReversed = false;

				//kminmerFile.read((char*)&length, sizeof(length));
				//kminmerFile.read((char*)&lengthStart, sizeof(lengthStart));
				//kminmerFile.read((char*)&lengthEnd, sizeof(lengthEnd));

				length = _kminmerSize-1;
				lengthStart = 1;
				lengthEnd = 1;
				//length = _kminmerSize;
				//lengthStart = 1;
				//lengthEnd = 1;
				//_logFile << size << " " << nodeName << endl;

				//_logFile << nodeName << endl;
				//auto found = _nodeName_entire.find(nodeName);
				//if(found != _nodeName_entire.end() && !found->second._loaded){
				//	_nodeName_entire[nodeName] = {minimizerSeq, true};
				//}

				bool lala = _nodeName_entire.modify_if(nodeName, 
				[&minimizerSeq](KminmerSequenceMap::value_type& v) { // key exist
					if(!v.second._loaded){
						v.second._loaded = true;
						v.second._minimizers = minimizerSeq;
					}
				});



				bool loulou = _nodeName_left.modify_if(nodeName, 
				[&lengthStart, &minimizerSeq](KminmerSequenceMap::value_type& v) { // key exist
					if(!v.second._loaded){

						vector<MinimizerType> minimizerSeq_sub;

						u_int32_t startPosition = 0; //kminmerInfo._read_pos_start;
						u_int32_t len = startPosition+lengthStart;

						for(size_t i=startPosition; i<len; i++){
							minimizerSeq_sub.push_back(minimizerSeq[i]);
						}

						v.second._loaded = true;
						v.second._minimizers = minimizerSeq_sub;
					}
				});


				_nodeName_right.modify_if(nodeName, 
				[&lengthEnd, &minimizerSeq](KminmerSequenceMap::value_type& v) { // key exist
					if(!v.second._loaded){
			
						vector<MinimizerType> minimizerSeq_sub;

						u_int32_t startPosition = minimizerSeq.size() - lengthEnd;  //kminmerInfo._position_of_second_to_last_minimizer_seq;
						u_int32_t len = startPosition+lengthEnd; //kminmerInfo._read_pos_end

						//_logFile << "right: " << startPosition << " " << len << endl;

						for(size_t i=startPosition; i<len; i++){
							//_logFile << i << " " << " " << minimizerSeq.size() << endl;

							minimizerSeq_sub.push_back(minimizerSeq[i]);
							//_logFile << minimizerSeq[i] << " ";
						}

						v.second._loaded = true;
						v.second._minimizers = minimizerSeq_sub;
					}
				});



			}
		}

		kminmerFile.close();
	}
	

	/*
	void extractKminmerSequences (){


		_logFile << "Extracting kminmers (reads)" << endl;
		KminmerParser parser(_filename_readMinimizers, _minimizerSize, _kminmerSize, true);
		auto fp = std::bind(&ToMinspace::extractKminmerSequences_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
		parser.parseMinspace(fp);
	}
	*/
	

	ofstream _outputFile;

	void createMinimizerContigs(){

		_nextReadIndexWriter = 0;

		_checksum = 0;
		size_t firstK = _kminmerSizeFirst;

		_logFile << endl;
		_logFile << "Creating mContigs" << endl;
		_logFile << endl;

		//_mdbg = new MDBG(firstK);
		//_mdbg->load(_inputDir + "/kminmerData_min_init.txt", false);

		//ofstream testCsv(_inputDir + "/minimizerContigsDebug.csv");
		//testCsv << "Name,Colour" << endl;

		//gzFile outputContigFile = gzopen(_filename_outputContigs.c_str(),"wb");

		_outputFile = ofstream(_filename_output);

		NodePathParserParallel parser2(_filename_inputContigs, _nbCores);
		parser2.parse(CreateMinimizerContigsFunctor(*this));


		_outputFile.close();

		/*
		//gzFile contigFile = gzopen(_filename_inputContigs.c_str(), "rb");
		ifstream contigFile(_filename_inputContigs);

		u_int64_t contig_index = 0;

		int nbFailed = 0;
		while(true){

			vector<u_int32_t> nodePath;
			//vector<u_int64_t> supportingReads;
			u_int64_t size;
			contigFile.read((char*)&size, sizeof(size));
			

			if(contigFile.eof()) break;


			u_int8_t isCircular;
			contigFile.read((char*)&isCircular, sizeof(isCircular));
			//u_int8_t isCircular;
			//gzread(contigFile, (char*)&isCircular, sizeof(isCircular));

			nodePath.resize(size);
			//supportingReads.resize(size);
			contigFile.read((char*)&nodePath[0], size * sizeof(u_int32_t));
			//gzread(contigFile, (char*)&supportingReads[0], size * sizeof(u_int64_t));


			vector<u_int64_t> contigSequence;

			long checkmSumLocal = 0;

			for(size_t i=0; i<nodePath.size(); i++){
				
				u_int32_t nodeIndex = nodePath[i];

				checkmSumLocal += nodeIndex;

				bool orientation;
				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex, orientation);

				if(i == 0){
					
					
					if(orientation){ //+

						vector<u_int64_t> minimizers = _nodeName_entire[nodeName]._minimizers;
						for(u_int64_t minimizer : minimizers) contigSequence.push_back(minimizer);
					}
					else{
						vector<u_int64_t> minimizers = _nodeName_entire[nodeName]._minimizers;
						std::reverse(minimizers.begin(), minimizers.end());
						for(u_int64_t minimizer : minimizers) contigSequence.push_back(minimizer);

					}
				}
				else {
					if(orientation){
						if(_nodeName_right.find(nodeName) == _nodeName_right.end()){
							 _logFile << "omg: " << nodeName << endl;
							 //getchar();
						}
						vector<u_int64_t> minimizers = _nodeName_right[nodeName]._minimizer;
						for(u_int64_t minimizer : minimizers) contigSequence.push_back(minimizer);
					}
					else{
						if(_nodeName_left.find(nodeName) == _nodeName_left.end()){
							 _logFile << "omg: " << nodeName << endl;
							 //getchar();
						}
						vector<u_int64_t> minimizers = _nodeName_left[nodeName]._minimizer;
						for(u_int64_t minimizer : minimizers) contigSequence.push_back(minimizer);
					}
				}

				_checksum += checkmSumLocal*nodePath.size();
			}



			u_int32_t contigSize = contigSequence.size();
			//_logFile << "Write: " << contigSize << endl;
			//getchar();
			outputFile.write((const char*)&contigSize, sizeof(contigSize));
			outputFile.write((const char*)&isCircular, sizeof(isCircular));
			outputFile.write((const char*)&contigSequence[0], contigSize*sizeof(u_int64_t));
			*/


			
			/*
			vector<KmerVec> kminmers; 
			vector<ReadKminmer> kminmersInfo;
			vector<u_int64_t> minimizersPos; 
			vector<u_int64_t> rlePositions;
			MDBG::getKminmers(_minimizerSize, firstK, contigSequence, minimizersPos, kminmers, kminmersInfo, rlePositions, -1, false);

			u_int64_t nbFoundMinimizers = 0;
			for(size_t i=0; i<kminmers.size(); i++){
				KmerVec& vec = kminmers[i];
				

				if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()){
					_logFile << "Not good: " << i << endl;
					continue;	
				}

				u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;

				nbFoundMinimizers += 1;

			}

			//_logFile << "Is circular: " << ((int)isCircular) << endl;
			//_logFile << "Found nodes: " << nbFoundMinimizers << endl;
			if(nbFoundMinimizers != kminmers.size()){
				_logFile << "Nb kminmers: " << kminmers.size() << endl;
				_logFile << "Found nodes: " << nbFoundMinimizers << endl;
				_logFile << nodePath.size() << " " << contigSequence.size() << endl;
				_logFile << "issue minspace" << endl;
				//getchar();
			}
			//if(contigSequence.size() > 10000)
			//getchar();
			

			contig_index += 1;
		}


		*/
		//_logFile << nbFailed << " " << contig_index << endl;
		//contigFile.close();
		//outputFile.close();
		
		//_logFile << "Checksum: " << _checksum << endl;
		//gzclose(outputContigFile);
		//testCsv.close();

		//if(_kminmerSize == 21) exit(1);


	}


	class CreateMinimizerContigsFunctor { 

		public:

		ToMinspace& _graph;

		CreateMinimizerContigsFunctor(ToMinspace& graph) : _graph(graph){
		}

		CreateMinimizerContigsFunctor(const CreateMinimizerContigsFunctor& copy) : _graph(copy._graph){
		}

		~CreateMinimizerContigsFunctor(){
		}
		
		void operator () (const NodePath& nodePathObject) {
			
				
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
							 _graph._logFile << "omg: " << nodeName << endl;
							 //getchar();
						}
						vector<MinimizerType> minimizers = _graph._nodeName_right[nodeName]._minimizers;
						for(MinimizerType minimizer : minimizers) contigSequence.push_back(minimizer);
					}
					else{
						if(_graph._nodeName_left.find(nodeName) == _graph._nodeName_left.end()){
							 _graph._logFile << "omg: " << nodeName << endl;
							 //getchar();
						}
						vector<MinimizerType> minimizers = _graph._nodeName_left[nodeName]._minimizers;
						for(MinimizerType minimizer : minimizers) contigSequence.push_back(minimizer);
					}
				}

				//_checksum += checkmSumLocal*nodePath.size();
			}

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
			_graph.writeRead(nodePathObject._readIndex, contigSequence, nodePathObject._isCircular);
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

	void writeRead(u_int64_t readIndex, const vector<MinimizerType>& contigSequence, u_int8_t isCircular){

		//#pragma omp critical(dataupdate)
		#pragma omp critical
		{
			_readWriterQueue.push({readIndex, contigSequence, isCircular});
			//_logFile << _readWriterQueue.size() << " " << read._index << " " << _nextReadIndexWriter << endl;

			while(!_readWriterQueue.empty()){

				const ReadWriter& readWriter = _readWriterQueue.top();

				if(readWriter._readIndex == _nextReadIndexWriter){


					u_int32_t contigSize = readWriter._contigSequence.size();
					_outputFile.write((const char*)&contigSize, sizeof(contigSize));
					_outputFile.write((const char*)&readWriter._isCircular, sizeof(readWriter._isCircular));
					_outputFile.write((const char*)&readWriter._contigSequence[0], contigSize*sizeof(MinimizerType));

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



