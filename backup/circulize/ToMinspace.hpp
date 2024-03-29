

#ifndef MDBG_METAG_TOMINSPACE
#define MDBG_METAG_TOMINSPACE

#include "../Commons.hpp"

/*
struct KminmerSequence{
	u_int32_t _readIndex;
	vector<u_int64_t> _minimizers;
};*/


struct KminmerSequenceEntire{
	vector<u_int64_t> _minimizers;
	bool _loaded;
};

struct KminmerSequenceLeftRight{
	vector<u_int64_t> _minimizer;
	bool _loaded;
};

class ToMinspace : public Tool{
    
public:

	//string _inputFilename;
	string _inputDir;
	float _minimizerDensity;
    size_t _minimizerSize;
    size_t _kminmerSize;
    size_t _kminmerSizeFirst;

	string _filename_joints;
	string _filename_inputContigs;
	string _filename_readMinimizers;
	//string _filename_outputContigs;
	string _filename_kminmerSequences;
	string _filename_output;
	//string _filename_nodeContigs;
	MDBG* _mdbg;
	
	unordered_map<u_int32_t, KminmerSequenceEntire> _nodeName_entire;
	//unordered_map<u_int32_t, KminmerSequenceEntire> _nodeName_entireRC;
	unordered_map<u_int32_t, KminmerSequenceLeftRight> _nodeName_right;
	unordered_map<u_int32_t, KminmerSequenceLeftRight> _nodeName_left;

	u_int64_t _checksum;

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

		//_nbCores = args::get(arg_nbCores);

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
		_filename_joints = _inputDir + "/joint_data.txt";
		_filename_readMinimizers = _inputDir + "/read_data.txt";
		//_filename_readMinimizers_output = _filename_readMinimizers + ".corrected.txt";
		//_filename_nodeContigs = _inputDir + "/contigs.nodepath.gz";
	}


    void execute (){
    
		loadContigs(_filename_inputContigs);
		loadContigs(_filename_joints);

		_logFile << "Loading mdbg" << endl;
		string mdbg_filename = _inputDir + "/mdbg_nodes.gz";
		//_mdbg = new MDBG(_kminmerSize);
		//_mdbg->load(mdbg_filename, false);
		//_logFile << "MDBG nodes: " << _mdbg->_dbg_nodes.size() << endl;

		extractKminmerSequences();
		//_logFile << "delete mdbg a remetre !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
		//delete _mdbg;

		ofstream outputFile(_filename_output);
		createMinimizerContigs(_filename_inputContigs, outputFile, false);

		
		ofstream fileReadMinimizers = ofstream(_filename_readMinimizers, std::ofstream::app);
		createMinimizerContigs(_filename_joints, fileReadMinimizers, true);
		//dereplicateContigs();

		_logFile << endl << "Contig filename: " << _filename_output << endl;
		closeLogFile();
	}



	void loadContigs(const string& filename){

		_logFile << "Loading contigs" << endl;
		ifstream contigFile(filename);
		u_int64_t nbContigs = 0;
		
		while(true){

			vector<u_int32_t> nodePath;
			//vector<u_int64_t> supportingReads;
			u_int64_t size;
			contigFile.read((char*)&size, sizeof(size));
			

			if(contigFile.eof()) break;

			bool isCircular;
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
		
	}

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

		_logFile << "Extracting kminmer sequences" << endl;
		
		
		ifstream kminmerFile(_inputDir + "/kminmerData_min.txt");

		while (true) {

			u_int16_t size = _kminmerSize;
			//kminmerFile.read((char*)&size, sizeof(size));


			vector<u_int64_t> minimizerSeq;
			minimizerSeq.resize(size);
			kminmerFile.read((char*)&minimizerSeq[0], size*sizeof(u_int64_t));

			if(kminmerFile.eof())break;

			u_int32_t nodeName;
			u_int32_t length;
			u_int32_t lengthStart;
			u_int32_t lengthEnd;
			u_int32_t abundance;
			u_int32_t quality;
			//bool isReversed = false;

			kminmerFile.read((char*)&nodeName, sizeof(nodeName));
			kminmerFile.read((char*)&abundance, sizeof(abundance));
			kminmerFile.read((char*)&quality, sizeof(quality));
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
			auto found = _nodeName_entire.find(nodeName);
			if(found != _nodeName_entire.end() && !found->second._loaded){

				//_logFile << "entire" << endl;
				
				//u_int32_t startPosition = 0; //kminmerInfo._read_pos_start;
				//u_int32_t len = length; //startPosition + kminmerInfo._read_pos_end - kminmerInfo._read_pos_start;
				

				//_logFile << startPosition << " " << len << " "  << minimizerSeq.size() << endl;
				//vector<u_int64_t> minimizerSeq_sub; 
				//for(size_t i=startPosition; i <= len; i++){
					//_logFile << i << " " << " " << minimizerSeq.size() << endl;
					//minimizerSeq_sub.push_back(minimizerSeq[i]);
					//_logFile << minimizerSeq[i] << " ";
				//}
			
				//_logFile << startPosition << " " << len << " " << minimizerSeq_sub.size() << endl;

				//if(kminmerInfo._isReversed){
				//	std::reverse(vec._kmers.begin(), vec._kmers.end());
				//}
				//if(kminmerInfo._isReversed){
					//std::reverse(minimizerSeq.begin(), minimizerSeq.end());
				//}
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
				u_int32_t len = startPosition+lengthStart;

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

				u_int32_t startPosition = minimizerSeq.size() - lengthEnd;  //kminmerInfo._position_of_second_to_last_minimizer_seq;
				u_int32_t len = startPosition+lengthEnd; //kminmerInfo._read_pos_end

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
	

	void createMinimizerContigs(const string& filename, ofstream& outputFile, bool isJoint){

		_checksum = 0;
		size_t firstK = _kminmerSizeFirst;

		_logFile << endl;
		_logFile << "Creating mContigs" << endl;
		_logFile << endl;

		_mdbg = new MDBG(firstK);
		_mdbg->load(_inputDir + "/kminmerData_min_init.txt", false);

		//ofstream testCsv(_inputDir + "/minimizerContigsDebug.csv");
		//testCsv << "Name,Colour" << endl;

		//gzFile outputContigFile = gzopen(_filename_outputContigs.c_str(),"wb");



		//gzFile contigFile = gzopen(_filename_inputContigs.c_str(), "rb");
		ifstream contigFile(filename);

		u_int64_t contig_index = 0;

		int nbFailed = 0;
		while(true){

			vector<u_int32_t> nodePath;
			//vector<u_int64_t> supportingReads;
			u_int64_t size;
			contigFile.read((char*)&size, sizeof(size));
			

			if(contigFile.eof()) break;


			bool isCircular;
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
				//_checksum += nodeIndex; //BiGraph::nodeIndex_to_nodeName(nodeIndex);

				bool orientation;
				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex, orientation);
				//u_int64_t readIndex = supportingReads[i];
				//ContigNode contigNode = {nodeName, readIndex};

				//if(12524 == nodeName){
				//	_logFile << BiGraph::nodeIndex_to_nodeName(nodePath[i+1], orientation) << endl;
				//	getchar();
				//}
				if(i == 0){
					
					//_logFile << endl << endl << endl << endl << "Entire" << endl;
					
					if(orientation){ //+
						//_logFile << endl << endl << endl << endl << "Entire" << endl;

						vector<u_int64_t> minimizers = _nodeName_entire[nodeName]._minimizers;
						for(u_int64_t minimizer : minimizers) contigSequence.push_back(minimizer);
					}
					else{
						//_logFile << endl << endl << endl << endl << "Entire RC" << endl;
						//_logFile << _debugMinimizers_isReversed[nodeName] << endl;
						vector<u_int64_t> minimizers = _nodeName_entire[nodeName]._minimizers;
						std::reverse(minimizers.begin(), minimizers.end());
						for(u_int64_t minimizer : minimizers) contigSequence.push_back(minimizer);

					}
				}
				else {
					if(orientation){
						//if(i > 490 && i < 510){
						//	_logFile << "HIHIHI: " << i << " " << _nodeName_right[nodeName]._minimizer << endl;
						//}
						if(_nodeName_right.find(nodeName) == _nodeName_right.end()){
							 _logFile << "omg: " << nodeName << endl;
							 //getchar();
						}
						vector<u_int64_t> minimizers = _nodeName_right[nodeName]._minimizer;
						for(u_int64_t minimizer : minimizers) contigSequence.push_back(minimizer);
						//contigSequence.push_back(minimizer);

						//if(minimizer == 0){ _logFile <<"HAAAAAAAAAAAAAAAAAAAA" << endl; exit(1);}
						//_logFile << "Add minimizer (right): " << minimizer << endl;
					}
					else{
						//if(i > 490 && i < 510){
						//	_logFile << "HIHIHI: " << i << " " << _nodeName_left[nodeName]._minimizer << endl;
						//}
						if(_nodeName_left.find(nodeName) == _nodeName_left.end()){
							 _logFile << "omg: " << nodeName << endl;
							 //getchar();
						}
						vector<u_int64_t> minimizers = _nodeName_left[nodeName]._minimizer;
						for(u_int64_t minimizer : minimizers) contigSequence.push_back(minimizer);
						//if(minimizer == 0){ _logFile <<"HAAAAAAAAAAAAAAAAAAAA" << endl; exit(1);}
						//contigSequence.push_back(minimizer);
						//_logFile << "Add minimizer (left): " << minimizer << endl;
					}
				}

				/*
				if(i == 0){

					_logFile << "Next node: " << nodeName << endl;
				
					if(_debugMinimizers.find(nodeName) != _debugMinimizers.end()){
						_logFile << nodeName << endl;
						for(u_int64_t minimizer : _debugMinimizers[nodeName]){
							_logFile << minimizer << " ";
						}
						_logFile << endl;
					}
				}
				*/

				/*
				_logFile << "----" << endl;
				for(size_t i=0; i<contigSequence.size(); i++){
					_logFile << contigSequence[i] << " ";
				}
				_logFile << endl;

				for(auto& it : mdbg->_dbg_nodes){
					if(it.second._index == nodeName){
						_logFile << it.first._kmers[0] << " " << it.first._kmers[1] << " " << it.first._kmers[2] << endl;
					}
				}

				if(i > 5) break;
				*/
				
				/*
				vector<KmerVec> kminmers; 
				vector<ReadKminmer> kminmersInfo;
				vector<u_int64_t> minimizersPos; 
				vector<u_int64_t> rlePositions;
				MDBG::getKminmers(21, 3, contigSequence, minimizersPos, kminmers, kminmersInfo, rlePositions, -1);

				_logFile << "-------------------" << kminmers.size() << endl;
				u_int64_t nbFoundMinimizers = 0;
				for(size_t i=0; i<kminmers.size(); i++){
					KmerVec& vec = kminmers[i];
					

					if(mdbg->_dbg_nodes.find(vec) == mdbg->_dbg_nodes.end()){
						_logFile << vec._kmers[0] << " " << vec._kmers[1] << " " << vec._kmers[2] << endl;

						
						for(auto& it : mdbg->_dbg_nodes){
							if(it.second._index == 6247){
								_logFile << it.first._kmers[0] << " " << it.first._kmers[1] << " " << it.first._kmers[2] << endl;
							}
						}
						
						getchar();
					}

					u_int32_t nodeName = mdbg->_dbg_nodes[vec]._index;
					_logFile << nodeName << endl;
					//_logFile << vec._kmers[0] << " " << vec._kmers[1] << " " << vec._kmers[2] << "     " << nodeName << endl;



				}*/


				_checksum += checkmSumLocal*nodePath.size();
			}

			//if(isCircular){
				//_logFile << "blurp: " << contigSequence[_kminmerSize-1] << endl;
				//contigSequence.pop_back();
				//contigSequence.push_back(contigSequence[_kminmerSize-1]); //k += 1
				//contigSequence.push_back(contigSequence[_kminmerSize]); //k += 2
			//}
			//exit(1);
			//_logFile << contigSequence.substr(362170, 100) << endl;
			//_logFile << endl << endl;

			//if(contigSequence.size() > 3000){
			//_logFile << nodePath.size() << " " << contigSequence.size() << endl;
			//}
			//_logFile << contigSequence << endl;

			if(isJoint) cout << contigSequence.size() << endl;
			u_int32_t contigSize = contigSequence.size();
			//_logFile << "Write: " << contigSize << endl;
			//getchar();
			outputFile.write((const char*)&contigSize, sizeof(contigSize));
			outputFile.write((const char*)&isCircular, sizeof(isCircular));
			outputFile.write((const char*)&contigSequence[0], contigSize*sizeof(u_int64_t));

			//_logFile << contigSize << " " << isCircular << endl;
			//_logFile << contigSize << endl;
			//if(contigSize <= 0){
			//	_logFile << "Empty contig: " << nodePath.size() << endl;
			//	getchar();
			//}

			//vector<u_int16_t> minimizerPosOffset(contigSequence.size(), 0);
			//outputFile.write((const char*)&minimizerPosOffset[0], size*sizeof(u_int16_t));

			//gzwrite(outputContigFile, (const char*)&contigSize, sizeof(contigSize));
			//gzwrite(outputContigFile, (const char*)&contigSequence[0], contigSize * sizeof(u_int64_t));

			/*
			if(contigSequence.size() > 4000){
				_logFile << _kminmerSize << endl;
				_logFile << "---" << endl;
				for(size_t i=0; i<_kminmerSize; i++){
					_logFile << contigSequence[i] << " ";
				}
				_logFile << endl;
				_logFile << "---" << endl;
				for(size_t i=contigSequence.size()-_kminmerSize-1; i<contigSequence.size(); i++){
					_logFile << contigSequence[i] << " ";
				}
				_logFile << endl;
				//getchar();
			}*/
			
			
			vector<KmerVec> kminmers; 
			vector<ReadKminmer> kminmersInfo;
			vector<u_int64_t> minimizersPos; 
			vector<u_int64_t> rlePositions;
			MDBG::getKminmers(_minimizerSize, firstK, contigSequence, minimizersPos, kminmers, kminmersInfo, rlePositions, -1, false);

			u_int64_t nbFoundMinimizers = 0;
			for(size_t i=0; i<kminmers.size(); i++){
				KmerVec& vec = kminmers[i];
				
				//if(vec._kmers[0] == 66918863945726617 && vec._kmers[1] == 46622693399843280 && vec._kmers[2] == 91239340561015544){
				//if(vec._kmers[0] == 4901087538459952 && vec._kmers[1] == 14618695441502469 && vec._kmers[2] == 79737352576542472){
				//	_logFile << "HAHAHAHA: " << kminmers.size() << endl;
					//getchar();
				//}

				/*
				if(i  == 0){
					_logFile << "---" << endl;
					for(size_t i=0; i<vec._kmers.size(); i++){
						_logFile << vec._kmers[i] << " ";
					}
					_logFile << endl;
				}
				else if(i == kminmers.size()-1){
					_logFile << "---" << endl;
					for(size_t i=0; i<vec._kmers.size(); i++){
						_logFile << vec._kmers[i] << " ";
					}
					_logFile << endl;
				}*/
				//if(i==0){
				//	_logFile << vec._kmers[0] << " " << vec._kmers[1] << " " << vec._kmers[2] << endl;
				//}
				//_logFile << "----" << endl;
				//_logFile << vec._kmers[0] << " " << vec._kmers[1] << " " << vec._kmers[2] << endl;
	
				//if(i>5) break;
				//_logFile << _kminmerSize << endl;
				//if(i > 490 && i < 500){
				//	_logFile << vec._kmers[0] << " " << vec._kmers[1] << " " << vec._kmers[2] << endl;
				//}

				if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()){
					//if(i==2){ nbFailed += 1; }
					_logFile << "Not good: " << i << endl;
					//_logFile << vec._kmers[0] << " " << vec._kmers[1] << " " << vec._kmers[2] << " " << vec._kmers[3] << endl;
					
					//if(mdbg->_dbg_nodes.find(kminmers[i-1]) != mdbg->_dbg_nodes.end()){
					//	_logFile << mdbg->_dbg_nodes[kminmers[i-1]]._index << endl;
					//	_logFile << kminmers[i-1]._kmers[0] << " " << kminmers[i-1]._kmers[1] << " " << kminmers[i-1]._kmers[2] << endl;
					//}
					//if(mdbg->_dbg_nodes.find(kminmers[i+1]) != mdbg->_dbg_nodes.end()){
					//	_logFile << mdbg->_dbg_nodes[kminmers[i+1]]._index << endl;
					//	_logFile << kminmers[i-1]._kmers[0] << " " << kminmers[i-1]._kmers[1] << " " << kminmers[i-1]._kmers[2] << endl;
					//}

					//getchar();
					continue;	
				}

				u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
				//_logFile << nodeName << endl;
				//_logFile << vec._kmers[0] << " " << vec._kmers[1] << " " << vec._kmers[2] << "     " << nodeName << endl;
				//_logFile << nodeName << endl;

				nbFoundMinimizers += 1;
				//ReadKminmer& kminmerInfo = kminmersInfo[i];

				//testCsv << nodeName << "," << contig_index << endl;

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


		//_logFile << nbFailed << " " << contig_index << endl;
		contigFile.close();
		outputFile.close();
		
		_logFile << "Checksum: " << _checksum << endl;
		//gzclose(outputContigFile);
		//testCsv.close();

		//if(_kminmerSize == 21) exit(1);


	}


	/*
	struct Contig{
		u_int64_t _readIndex;
		vector<u_int32_t> _nodepath;
		vector<u_int32_t> _nodepath_sorted;
	};

	vector<Contig> _contigs;
	unordered_map<u_int32_t, u_int32_t> _kminmerCounts;
	static bool ContigComparator_ByLength(const Contig &a, const Contig &b){

		if(a._nodepath.size() == b._nodepath.size()){
			for(size_t i=0; i<a._nodepath.size() && i<b._nodepath.size(); i++){
				if(BiGraph::nodeIndex_to_nodeName(a._nodepath[i]) == BiGraph::nodeIndex_to_nodeName(b._nodepath[i])){
					continue;
				}
				else{
					return BiGraph::nodeIndex_to_nodeName(a._nodepath[i]) > BiGraph::nodeIndex_to_nodeName(b._nodepath[i]);
				}
			}
		}


		return a._nodepath.size() > b._nodepath.size();
	}
	
	unordered_set<u_int32_t> _invalidContigIndex;
	u_int64_t _nbContigs;
	ofstream _outputFileDerep;


	void dereplicateContigs(){

		_logFile << "Dereplicating contigs" << endl;


		_outputFileDerep.open(_filename_output); 
		_nbContigs = 0;

		KminmerParser parser(_filename_output + ".tmp", _minimizerSize, _kminmerSizeFirst, false, false);
		auto fp = std::bind(&ToMinspace::loadContigs_min_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
		parser.parseMinspace(fp);

		_logFile << "Nb contigs: " << _nbContigs << endl;

		double nbRepeatedKminmers = 0;
		for(auto& it : _kminmerCounts){
			if(it.second > 1){
				nbRepeatedKminmers += 1;
			}
		}
		_logFile << "Nb kminmer: " << _kminmerCounts.size() << endl;
		_logFile << "Repeated kminmer: " << nbRepeatedKminmers << endl;
		_logFile << "Repeated kminmer rate: " << (nbRepeatedKminmers / _kminmerCounts.size()) << endl;

		
		std::sort(_contigs.begin(), _contigs.end(), ContigComparator_ByLength);

		for(size_t i=0; i<_contigs.size(); i++){
			
			if(_invalidContigIndex.find(_contigs[i]._readIndex) != _invalidContigIndex.end()) continue;

			for(long j=_contigs.size()-1; j>=i+1; j--){
			
				if(_invalidContigIndex.find(_contigs[j]._readIndex) != _invalidContigIndex.end()) continue;

				double nbShared = Utils::computeSharedElements(_contigs[i]._nodepath_sorted, _contigs[j]._nodepath_sorted);
				double sharedRate_1 = nbShared / _contigs[i]._nodepath_sorted.size();
				double sharedRate_2 = nbShared / _contigs[j]._nodepath_sorted.size();

				if(sharedRate_1 > 0.9 || sharedRate_2 > 0.9){

					_logFile << _contigs[j]._nodepath_sorted.size() << " " << sharedRate_1 << " " << sharedRate_2 << endl;
					_invalidContigIndex.insert(_contigs[j]._readIndex);
					//break;

				}

				//if(_contigs[j]._nodepath.size() < 100) _invalidContigIndex.insert(_contigs[j]._readIndex);
			}
		}
		


		_kminmerCounts.clear();

		_logFile << _nbContigs << " " << _invalidContigIndex.size() << endl;
 		_logFile << "Nb contigs (no duplicate): " << (_nbContigs-_invalidContigIndex.size()) << endl;

		KminmerParser parser2(_filename_output + ".tmp", _minimizerSize, _kminmerSizeFirst, false, false);
		auto fp2 = std::bind(&ToMinspace::loadContigs_min_read2, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
		parser2.parseMinspace(fp2);

		nbRepeatedKminmers = 0;
		for(auto& it : _kminmerCounts){
			if(it.second > 1){
				nbRepeatedKminmers += 1;
			}
		}
		_logFile << "Nb kminmer: " << _kminmerCounts.size() << endl;
		_logFile << "Repeated kminmer: " << nbRepeatedKminmers << endl;
		_logFile << "Repeated kminmer rate: " << (nbRepeatedKminmers / _kminmerCounts.size()) << endl;

		_contigs.clear();
		_outputFileDerep.close();
	}





	void loadContigs_min_read(const vector<u_int64_t>& readMinimizers, const vector<ReadKminmerComplete>& kminmersInfos, u_int64_t readIndex){

		vector<u_int32_t> nodepath;
		
		for(size_t i=0; i<kminmersInfos.size(); i++){
			
			const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];

			KmerVec vec = kminmerInfo._vec;
			
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()){
				_logFile << "Not found kminmer" << endl;
				//getchar();
				continue;
			}


			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
			nodepath.push_back(nodeName);
			
			_kminmerCounts[nodeName] += 1;

		}

		_nbContigs += 1;

		vector<u_int32_t> nodepath_sorted = nodepath;
		std::sort(nodepath_sorted.begin(), nodepath_sorted.end());
		_contigs.push_back({readIndex, nodepath, nodepath_sorted});

	}

	void loadContigs_min_read2(const vector<u_int64_t>& readMinimizers, const vector<ReadKminmerComplete>& kminmersInfos, u_int64_t readIndex){

		if(_invalidContigIndex.find(readIndex) != _invalidContigIndex.end()) return;
		
		u_int32_t contigSize = readMinimizers.size();
		//_logFile << "Write: " << contigSize << endl;
		//getchar();
		_outputFileDerep.write((const char*)&contigSize, sizeof(contigSize));
		_outputFileDerep.write((const char*)&readMinimizers[0], contigSize*sizeof(u_int64_t));

		for(size_t i=0; i<kminmersInfos.size(); i++){
			
			const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];

			KmerVec vec = kminmerInfo._vec;
			
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()){
				continue;
			}

			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
			
			_kminmerCounts[nodeName] += 1;

		}
	}
	*/

};	


#endif 



