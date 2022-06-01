

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
		/*
		//_inputFilename = getInput()->getStr(STR_INPUT);
		_inputDir = getInput()->getStr(STR_INPUT_DIR);
		*/
		
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
			cout << options.help() << endl;
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
			std::cout << options.help() << std::endl;
			std::cerr << e.what() << std::endl;
			std::exit(EXIT_FAILURE);
		}

		string filename_parameters = _inputDir + "/parameters.gz";
		gzFile file_parameters = gzopen(filename_parameters.c_str(),"rb");
		gzread(file_parameters, (char*)&_minimizerSize, sizeof(_minimizerSize));
		gzread(file_parameters, (char*)&_kminmerSize, sizeof(_kminmerSize));
		gzread(file_parameters, (char*)&_minimizerDensity, sizeof(_minimizerDensity));
		gzread(file_parameters, (char*)&_kminmerSizeFirst, sizeof(_kminmerSizeFirst));
		gzclose(file_parameters);


		cout << endl;
		cout << "Input dir: " << _inputDir << endl;
		//cout << "Output filename: " << _outputFilename << endl;
		cout << "Minimizer length: " << _minimizerSize << endl;
		cout << "Kminmer length: " << _kminmerSize << endl;
		cout << "Density: " << _minimizerDensity << endl;
		cout << endl;

		//_filename_outputContigs = _inputDir + "/contigs.min.gz";
		//_filename_output = _filename_inputContigs + ".min";
		//_filename_readMinimizers = _inputDir + "/readData_" + to_string(_kminmerSize) + ".txt";
		_filename_readMinimizers = _inputDir + "/read_data.txt";
		//_filename_readMinimizers_output = _filename_readMinimizers + ".corrected.txt";
		//_filename_nodeContigs = _inputDir + "/contigs.nodepath.gz";
	}


    void execute (){
    
		loadContigs();

		cout << "Loading mdbg" << endl;
		string mdbg_filename = _inputDir + "/mdbg_nodes.gz";
		//_mdbg = new MDBG(_kminmerSize);
		//_mdbg->load(mdbg_filename);
		//cout << "MDBG nodes: " << _mdbg->_dbg_nodes.size() << endl;

		extractKminmerSequences();
		//cout << "delete mdbg a remetre !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
		//delete _mdbg;

		createMinimizerContigs();

		cout << endl << "Contig filename: " << _filename_output << endl;
	}



	void loadContigs(){

		cout << "Loading contigs" << endl;
		ifstream contigFile(_filename_inputContigs);
		u_int64_t nbContigs = 0;
		
		while(true){

			vector<u_int32_t> nodePath;
			//vector<u_int64_t> supportingReads;
			u_int64_t size;
			contigFile.read((char*)&size, sizeof(size));
			

			if(contigFile.eof()) break;

			//u_int8_t isCircular;
			//gzread(contigFile, (char*)&isCircular, sizeof(isCircular));

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
					//cout << nodeName << endl;
					//if(orientation){ //+
					//	_nodeName_entire[nodeName] = {{}, false};
					//}
					//else{ //-
					//	_nodeName_entireRC[nodeName] = {{}, false};
					//}
				}
				else {
					//if(i==1) cout << orientation << endl;

					if(orientation){ //+
						_nodeName_right[nodeName] = {{}, false};
					}
					else{ //-
						_nodeName_left[nodeName] = {{}, false};
					}
				}

			}

			nbContigs += 1;
			//cout << nodePath.size() << endl;
		}

		contigFile.close();

		cout << _nodeName_entire.size() << " " << _nodeName_left.size() << " " << _nodeName_right.size() << endl;
		cout << "Nb contigs: " << nbContigs << endl;
		
	}

	/*
	void extractKminmerSequences_read(const vector<u_int64_t>& readMinimizers, const vector<ReadKminmerComplete>& kminmersInfos, u_int64_t readIndex){

		//cout << "----" << endl;
		//cout << readIndex << " " << kminmersInfos.size() << endl;

		for(size_t i=0; i<kminmersInfos.size(); i++){
			
			const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];

			KmerVec vec = kminmerInfo._vec;
			//for(KmerVec& vec : kminmers){
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()) continue;

			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;

			//cout << nodeName << endl;
			//cout << nodeName << " " << kminmerInfo._length << endl;
			
			vector<u_int64_t> minimizerSeq;

			//cout << kminmerInfo._read_pos_start << " " << kminmerInfo._read_pos_end << endl;

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

				//cout << "entire" << endl;
				
				u_int32_t startPosition = 0; //kminmerInfo._read_pos_start;
				u_int32_t len = kminmerInfo._read_pos_end - kminmerInfo._read_pos_start; //startPosition + kminmerInfo._read_pos_end - kminmerInfo._read_pos_start;
				

				//cout << startPosition << " " << len << " "  << minimizerSeq.size() << endl;
				vector<u_int64_t> minimizerSeq_sub; 
				for(size_t i=startPosition; i <= len; i++){
					//cout << i << " " << " " << minimizerSeq.size() << endl;
					minimizerSeq_sub.push_back(minimizerSeq[i]);
					//cout << minimizerSeq[i] << " ";
				}
			
				//cout << startPosition << " " << len << " " << minimizerSeq_sub.size() << endl;

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
				
				//cout << "left" << endl;
				//std::reverse(vec._kmers.begin(), vec._kmers.end());
				//if(12424 == nodeName){
				//	cout << "left: " <<  vec._kmers[0] << " " << vec._kmers[vec._kmers.size()-1] << endl;
				//}


				vector<u_int64_t> minimizerSeq_sub;

				u_int32_t startPosition = 0; //kminmerInfo._read_pos_start;
				u_int32_t len = startPosition+kminmerInfo._seq_length_start;

				//cout << "left: " << startPosition << " " << len << endl;
				
				for(size_t i=startPosition; i<len; i++){
					//cout << i << " " << " " << minimizerSeq.size() << endl;
					minimizerSeq_sub.push_back(minimizerSeq[i]);
					//cout << minimizerSeq[i] << " ";
				}
				//cout << endl;

				_nodeName_left[nodeName] = {minimizerSeq_sub, true};
				
				//getchar();
			}
			

			found2 = _nodeName_right.find(nodeName);
			if(_nodeName_right.find(nodeName) != _nodeName_right.end() && !found2->second._loaded){

				//cout << "right" << endl;

			
				vector<u_int64_t> minimizerSeq_sub;

				u_int32_t startPosition = minimizerSeq.size() - kminmerInfo._seq_length_end;  //kminmerInfo._position_of_second_to_last_minimizer_seq;
				u_int32_t len = startPosition+kminmerInfo._seq_length_end; //kminmerInfo._read_pos_end

				//cout << "right: " << startPosition << " " << len << endl;

				for(size_t i=startPosition; i<len; i++){
					//cout << i << " " << " " << minimizerSeq.size() << endl;

					minimizerSeq_sub.push_back(minimizerSeq[i]);
					//cout << minimizerSeq[i] << " ";
				}
				//cout << endl;

				_nodeName_right[nodeName] = {minimizerSeq_sub, true};
				//getchar();

			}

			//if(12424 == nodeName){
				//exit(1);
			//}

		}

		//cout << "done" << endl;
	}
	*/
	
	void extractKminmerSequences (){

		cout << "Extracting kminmer sequences" << endl;
		
		
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
			//bool isReversed = false;

			kminmerFile.read((char*)&nodeName, sizeof(nodeName));
			kminmerFile.read((char*)&abundance, sizeof(abundance));
			//kminmerFile.read((char*)&length, sizeof(length));
			//kminmerFile.read((char*)&lengthStart, sizeof(lengthStart));
			//kminmerFile.read((char*)&lengthEnd, sizeof(lengthEnd));

			length = _kminmerSize-1;
			lengthStart = 1;
			lengthEnd = 1;
			//length = _kminmerSize;
			//lengthStart = 1;
			//lengthEnd = 1;
			//cout << size << " " << nodeName << endl;

			//cout << nodeName << endl;
			auto found = _nodeName_entire.find(nodeName);
			if(found != _nodeName_entire.end() && !found->second._loaded){

				//cout << "entire" << endl;
				
				//u_int32_t startPosition = 0; //kminmerInfo._read_pos_start;
				//u_int32_t len = length; //startPosition + kminmerInfo._read_pos_end - kminmerInfo._read_pos_start;
				

				//cout << startPosition << " " << len << " "  << minimizerSeq.size() << endl;
				//vector<u_int64_t> minimizerSeq_sub; 
				//for(size_t i=startPosition; i <= len; i++){
					//cout << i << " " << " " << minimizerSeq.size() << endl;
					//minimizerSeq_sub.push_back(minimizerSeq[i]);
					//cout << minimizerSeq[i] << " ";
				//}
			
				//cout << startPosition << " " << len << " " << minimizerSeq_sub.size() << endl;

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
				
				//cout << "left" << endl;
				//std::reverse(vec._kmers.begin(), vec._kmers.end());
				//if(12424 == nodeName){
				//	cout << "left: " <<  vec._kmers[0] << " " << vec._kmers[vec._kmers.size()-1] << endl;
				//}


				vector<u_int64_t> minimizerSeq_sub;

				u_int32_t startPosition = 0; //kminmerInfo._read_pos_start;
				u_int32_t len = startPosition+lengthStart;

				//cout << "left: " << startPosition << " " << len << endl;
				
				for(size_t i=startPosition; i<len; i++){
					//cout << i << " " << " " << minimizerSeq.size() << endl;
					minimizerSeq_sub.push_back(minimizerSeq[i]);
					//cout << minimizerSeq[i] << " ";
				}
				//cout << endl;

				_nodeName_left[nodeName] = {minimizerSeq_sub, true};
				
				//getchar();
			}
			

			found2 = _nodeName_right.find(nodeName);
			if(_nodeName_right.find(nodeName) != _nodeName_right.end() && !found2->second._loaded){

				//cout << "right" << endl;

			
				vector<u_int64_t> minimizerSeq_sub;

				u_int32_t startPosition = minimizerSeq.size() - lengthEnd;  //kminmerInfo._position_of_second_to_last_minimizer_seq;
				u_int32_t len = startPosition+lengthEnd; //kminmerInfo._read_pos_end

				//cout << "right: " << startPosition << " " << len << endl;

				for(size_t i=startPosition; i<len; i++){
					//cout << i << " " << " " << minimizerSeq.size() << endl;

					minimizerSeq_sub.push_back(minimizerSeq[i]);
					//cout << minimizerSeq[i] << " ";
				}
				//cout << endl;

				_nodeName_right[nodeName] = {minimizerSeq_sub, true};
				//getchar();

			}


		}

		kminmerFile.close();
	}
	

	/*
	void extractKminmerSequences (){


		cout << "Extracting kminmers (reads)" << endl;
		KminmerParser parser(_filename_readMinimizers, _minimizerSize, _kminmerSize, true);
		auto fp = std::bind(&ToMinspace::extractKminmerSequences_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
		parser.parseMinspace(fp);
	}
	*/
	

	void createMinimizerContigs(){

		size_t firstK = _kminmerSizeFirst;

		cout << endl;
		cout << "Creating mContigs" << endl;
		cout << endl;

		string mdbg_filename =  _inputDir + "/kminmerData_min_init.txt";
		MDBG* mdbg = new MDBG(firstK);
		mdbg->load(mdbg_filename);

		//ofstream testCsv(_inputDir + "/minimizerContigsDebug.csv");
		//testCsv << "Name,Colour" << endl;

		//gzFile outputContigFile = gzopen(_filename_outputContigs.c_str(),"wb");
		ofstream outputFile(_filename_output, std::ios::binary);



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


			//u_int8_t isCircular;
			//gzread(contigFile, (char*)&isCircular, sizeof(isCircular));

			nodePath.resize(size);
			//supportingReads.resize(size);
			contigFile.read((char*)&nodePath[0], size * sizeof(u_int32_t));
			//gzread(contigFile, (char*)&supportingReads[0], size * sizeof(u_int64_t));


			vector<u_int64_t> contigSequence;

			for(size_t i=0; i<nodePath.size(); i++){
				
				u_int32_t nodeIndex = nodePath[i];
				bool orientation;
				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex, orientation);
				//u_int64_t readIndex = supportingReads[i];
				//ContigNode contigNode = {nodeName, readIndex};

				//if(12524 == nodeName){
				//	cout << BiGraph::nodeIndex_to_nodeName(nodePath[i+1], orientation) << endl;
				//	getchar();
				//}
				if(i == 0){
					
					//cout << endl << endl << endl << endl << "Entire" << endl;
					
					if(orientation){ //+
						//cout << endl << endl << endl << endl << "Entire" << endl;

						vector<u_int64_t> minimizers = _nodeName_entire[nodeName]._minimizers;
						for(u_int64_t minimizer : minimizers) contigSequence.push_back(minimizer);
					}
					else{
						//cout << endl << endl << endl << endl << "Entire RC" << endl;
						//cout << _debugMinimizers_isReversed[nodeName] << endl;
						vector<u_int64_t> minimizers = _nodeName_entire[nodeName]._minimizers;
						std::reverse(minimizers.begin(), minimizers.end());
						for(u_int64_t minimizer : minimizers) contigSequence.push_back(minimizer);

					}
				}
				else {
					if(orientation){
						//if(i > 490 && i < 510){
						//	cout << "HIHIHI: " << i << " " << _nodeName_right[nodeName]._minimizer << endl;
						//}
						if(_nodeName_right.find(nodeName) == _nodeName_right.end()){
							 cout << "omg: " << nodeName << endl;
							 //getchar();
						}
						vector<u_int64_t> minimizers = _nodeName_right[nodeName]._minimizer;
						for(u_int64_t minimizer : minimizers) contigSequence.push_back(minimizer);
						//contigSequence.push_back(minimizer);

						//if(minimizer == 0){ cout <<"HAAAAAAAAAAAAAAAAAAAA" << endl; exit(1);}
						//cout << "Add minimizer (right): " << minimizer << endl;
					}
					else{
						//if(i > 490 && i < 510){
						//	cout << "HIHIHI: " << i << " " << _nodeName_left[nodeName]._minimizer << endl;
						//}
						if(_nodeName_left.find(nodeName) == _nodeName_left.end()){
							 cout << "omg: " << nodeName << endl;
							 //getchar();
						}
						vector<u_int64_t> minimizers = _nodeName_left[nodeName]._minimizer;
						for(u_int64_t minimizer : minimizers) contigSequence.push_back(minimizer);
						//if(minimizer == 0){ cout <<"HAAAAAAAAAAAAAAAAAAAA" << endl; exit(1);}
						//contigSequence.push_back(minimizer);
						//cout << "Add minimizer (left): " << minimizer << endl;
					}
				}

				/*
				if(i == 0){

					cout << "Next node: " << nodeName << endl;
				
					if(_debugMinimizers.find(nodeName) != _debugMinimizers.end()){
						cout << nodeName << endl;
						for(u_int64_t minimizer : _debugMinimizers[nodeName]){
							cout << minimizer << " ";
						}
						cout << endl;
					}
				}
				*/

				/*
				cout << "----" << endl;
				for(size_t i=0; i<contigSequence.size(); i++){
					cout << contigSequence[i] << " ";
				}
				cout << endl;

				for(auto& it : mdbg->_dbg_nodes){
					if(it.second._index == nodeName){
						cout << it.first._kmers[0] << " " << it.first._kmers[1] << " " << it.first._kmers[2] << endl;
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

				cout << "-------------------" << kminmers.size() << endl;
				u_int64_t nbFoundMinimizers = 0;
				for(size_t i=0; i<kminmers.size(); i++){
					KmerVec& vec = kminmers[i];
					

					if(mdbg->_dbg_nodes.find(vec) == mdbg->_dbg_nodes.end()){
						cout << vec._kmers[0] << " " << vec._kmers[1] << " " << vec._kmers[2] << endl;

						
						for(auto& it : mdbg->_dbg_nodes){
							if(it.second._index == 6247){
								cout << it.first._kmers[0] << " " << it.first._kmers[1] << " " << it.first._kmers[2] << endl;
							}
						}
						
						getchar();
					}

					u_int32_t nodeName = mdbg->_dbg_nodes[vec]._index;
					cout << nodeName << endl;
					//cout << vec._kmers[0] << " " << vec._kmers[1] << " " << vec._kmers[2] << "     " << nodeName << endl;



				}*/


			}

			//if(isCircular){
				//cout << "blurp: " << contigSequence[_kminmerSize-1] << endl;
				//contigSequence.pop_back();
				//contigSequence.push_back(contigSequence[_kminmerSize-1]); //k += 1
				//contigSequence.push_back(contigSequence[_kminmerSize]); //k += 2
			//}
			//exit(1);
			//cout << contigSequence.substr(362170, 100) << endl;
			//cout << endl << endl;

			//if(contigSequence.size() > 3000){
			//cout << nodePath.size() << " " << contigSequence.size() << endl;
			//}
			//cout << contigSequence << endl;


			u_int32_t contigSize = contigSequence.size();
			//cout << "Write: " << contigSize << endl;
			//getchar();
			outputFile.write((const char*)&contigSize, sizeof(contigSize));
			outputFile.write((const char*)&contigSequence[0], contigSize*sizeof(u_int64_t));

			//cout << contigSize << endl;
			//if(contigSize <= 0){
			//	cout << "Empty contig: " << nodePath.size() << endl;
			//	getchar();
			//}

			//vector<u_int16_t> minimizerPosOffset(contigSequence.size(), 0);
			//outputFile.write((const char*)&minimizerPosOffset[0], size*sizeof(u_int16_t));

			//gzwrite(outputContigFile, (const char*)&contigSize, sizeof(contigSize));
			//gzwrite(outputContigFile, (const char*)&contigSequence[0], contigSize * sizeof(u_int64_t));

			/*
			if(contigSequence.size() > 4000){
				cout << _kminmerSize << endl;
				cout << "---" << endl;
				for(size_t i=0; i<_kminmerSize; i++){
					cout << contigSequence[i] << " ";
				}
				cout << endl;
				cout << "---" << endl;
				for(size_t i=contigSequence.size()-_kminmerSize-1; i<contigSequence.size(); i++){
					cout << contigSequence[i] << " ";
				}
				cout << endl;
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
				//	cout << "HAHAHAHA: " << kminmers.size() << endl;
					//getchar();
				//}

				/*
				if(i  == 0){
					cout << "---" << endl;
					for(size_t i=0; i<vec._kmers.size(); i++){
						cout << vec._kmers[i] << " ";
					}
					cout << endl;
				}
				else if(i == kminmers.size()-1){
					cout << "---" << endl;
					for(size_t i=0; i<vec._kmers.size(); i++){
						cout << vec._kmers[i] << " ";
					}
					cout << endl;
				}*/
				//if(i==0){
				//	cout << vec._kmers[0] << " " << vec._kmers[1] << " " << vec._kmers[2] << endl;
				//}
				//cout << "----" << endl;
				//cout << vec._kmers[0] << " " << vec._kmers[1] << " " << vec._kmers[2] << endl;
	
				//if(i>5) break;
				//cout << _kminmerSize << endl;
				//if(i > 490 && i < 500){
				//	cout << vec._kmers[0] << " " << vec._kmers[1] << " " << vec._kmers[2] << endl;
				//}

				if(mdbg->_dbg_nodes.find(vec) == mdbg->_dbg_nodes.end()){
					//if(i==2){ nbFailed += 1; }
					cout << "Not good: " << i << endl;
					//cout << vec._kmers[0] << " " << vec._kmers[1] << " " << vec._kmers[2] << " " << vec._kmers[3] << endl;
					
					//if(mdbg->_dbg_nodes.find(kminmers[i-1]) != mdbg->_dbg_nodes.end()){
					//	cout << mdbg->_dbg_nodes[kminmers[i-1]]._index << endl;
					//	cout << kminmers[i-1]._kmers[0] << " " << kminmers[i-1]._kmers[1] << " " << kminmers[i-1]._kmers[2] << endl;
					//}
					//if(mdbg->_dbg_nodes.find(kminmers[i+1]) != mdbg->_dbg_nodes.end()){
					//	cout << mdbg->_dbg_nodes[kminmers[i+1]]._index << endl;
					//	cout << kminmers[i-1]._kmers[0] << " " << kminmers[i-1]._kmers[1] << " " << kminmers[i-1]._kmers[2] << endl;
					//}

					//getchar();
					continue;	
				}

				u_int32_t nodeName = mdbg->_dbg_nodes[vec]._index;
				//cout << nodeName << endl;
				//cout << vec._kmers[0] << " " << vec._kmers[1] << " " << vec._kmers[2] << "     " << nodeName << endl;
				//cout << nodeName << endl;

				nbFoundMinimizers += 1;
				//ReadKminmer& kminmerInfo = kminmersInfo[i];

				//testCsv << nodeName << "," << contig_index << endl;

			}

			//cout << "Is circular: " << ((int)isCircular) << endl;
			//cout << "Found nodes: " << nbFoundMinimizers << endl;
			if(nbFoundMinimizers != kminmers.size()){
				cout << "Nb kminmers: " << kminmers.size() << endl;
				cout << "Found nodes: " << nbFoundMinimizers << endl;
				cout << nodePath.size() << " " << contigSequence.size() << endl;
				cout << "issue minspace" << endl;
				//getchar();
			}
			//if(contigSequence.size() > 10000)
			//getchar();
			

			contig_index += 1;
		}


		//cout << nbFailed << " " << contig_index << endl;
		contigFile.close();
		outputFile.close();
		//gzclose(outputContigFile);
		//testCsv.close();

		//if(_kminmerSize == 21) exit(1);


	}


};	


#endif 



