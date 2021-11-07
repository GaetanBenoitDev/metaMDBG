

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
	u_int64_t _minimizer;
	bool _loaded;
};

class ToMinspace : public Tool{
    
public:

	//string _inputFilename;
	string _inputDir;
	float _minimizerDensity;
    size_t _minimizerSize;
    size_t _kminmerSize;

	string _filename_readMinimizers;
	string _filename_outputContigs;
	string _filename_kminmerSequences;
	string _filename_nodeContigs;
	MDBG* _mdbg;
	
	unordered_map<u_int32_t, KminmerSequenceEntire> _nodeName_entire;
	unordered_map<u_int32_t, KminmerSequenceLeftRight> _nodeName_right;
	unordered_map<u_int32_t, KminmerSequenceLeftRight> _nodeName_left;

	//unordered_map<u_int32_t, vector<u_int64_t>> _debugMinimizers;

	ToMinspace(): Tool ("toMinspace"){

		//getParser()->push_back (new OptionOneParam (STR_INPUT, "input file", true));
		getParser()->push_back (new OptionOneParam (STR_INPUT_DIR, "input dir", true));
		//getParser()->push_back (new OptionOneParam (STR_OUTPUT, "output contig filename in basespace", true));

	}

    void execute (){
		parseArgs();
		loadContigs();

		cout << "Loading mdbg" << endl;
		string mdbg_filename = _inputDir + "/mdbg_nodes.gz";
		_mdbg = new MDBG(_kminmerSize);
		_mdbg->load(mdbg_filename);
		cout << "MDBG nodes: " << _mdbg->_dbg_nodes.size() << endl;


		extractKminmerSequences();
		//cout << "delete mdbg a remetre !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
		delete _mdbg;

		createMinimizerContigs();

		cout << endl << "Contig filename: " << _filename_outputContigs << endl;
	}

	void parseArgs(){
		//_inputFilename = getInput()->getStr(STR_INPUT);
		_inputDir = getInput()->getStr(STR_INPUT_DIR);

		string filename_parameters = _inputDir + "/parameters.gz";
		gzFile file_parameters = gzopen(filename_parameters.c_str(),"rb");
		gzread(file_parameters, (char*)&_minimizerSize, sizeof(_minimizerSize));
		gzread(file_parameters, (char*)&_kminmerSize, sizeof(_kminmerSize));
		gzread(file_parameters, (char*)&_minimizerDensity, sizeof(_minimizerDensity));
		gzclose(file_parameters);


		cout << endl;
		cout << "Input dir: " << _inputDir << endl;
		//cout << "Output filename: " << _outputFilename << endl;
		cout << "Minimizer length: " << _minimizerSize << endl;
		cout << "Kminmer length: " << _kminmerSize << endl;
		cout << "Density: " << _minimizerDensity << endl;
		cout << endl;

		_filename_outputContigs = _inputDir + "/contigs.min.gz";
		_filename_readMinimizers = _inputDir + "/read_data.gz";
		_filename_nodeContigs = _inputDir + "/minimizer_contigs.gz";
	}

	void loadContigs(){

		gzFile contigFile = gzopen(_filename_nodeContigs.c_str(),"rb");
		u_int64_t nbContigs = 0;
		
		while(true){

			vector<u_int32_t> nodePath;
			vector<u_int64_t> supportingReads;
			u_int64_t size;
			gzread(contigFile, (char*)&size, sizeof(size));
			

			if(gzeof(contigFile)) break;

			nodePath.resize(size);
			supportingReads.resize(size);
			gzread(contigFile, (char*)&nodePath[0], size * sizeof(u_int32_t));
			gzread(contigFile, (char*)&supportingReads[0], size * sizeof(u_int64_t));

			//for(u_int32_t nodeIndex : nodePath){
			for(size_t i=0; i<nodePath.size(); i++){
				u_int32_t nodeIndex = nodePath[i];
				bool orientation;
				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex, orientation);
				u_int32_t readIndex = supportingReads[i];


				if(i == 0){
					_nodeName_entire[nodeName] = {{}, false};
				}
				else {
					if(orientation){ //+
						_nodeName_right[nodeName] = {0, false};
					}
					else{ //-
						_nodeName_left[nodeName] = {0, false};
					}
				}

			}

			nbContigs += 1;
			//cout << nodePath.size() << endl;
		}

		gzclose(contigFile);

		cout << "Nb contigs: " << nbContigs << endl;
	}


	void extractKminmerSequences (){




		gzFile file_readData = gzopen(_filename_readMinimizers.c_str(),"rb");
		ReadIndexType readIndex = 0;

		while(true){
			
			//cout << readIndex << endl;
			//"reprise: essayer avec gfatools unitigs"
			u_int16_t size;
			vector<u_int64_t> minimizers;
			gzread(file_readData, (char*)&size, sizeof(size));

			if(gzeof(file_readData)) break;
			
			minimizers.resize(size);
			gzread(file_readData, (char*)&minimizers[0], size * sizeof(u_int64_t));


			vector<KmerVec> kminmers; 
			vector<ReadKminmer> kminmersInfo;
			vector<u_int64_t> minimizersPos; 
			vector<u_int64_t> rlePositions;
			MDBG::getKminmers(_minimizerSize, _kminmerSize, minimizers, minimizersPos, kminmers, kminmersInfo, rlePositions, readIndex);



			for(size_t i=0; i<kminmers.size(); i++){
				KmerVec& vec = kminmers[i];
				//for(KmerVec& vec : kminmers){
				if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()) continue;

				u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
				ReadKminmer& kminmerInfo = kminmersInfo[i];

				
				if(12424 == nodeName){
					for(u_int64_t minimizer : vec._kmers){
						cout << minimizer << " ";
					}
					cout << endl;
				}
				

				if(kminmerInfo._isReversed){
					std::reverse(vec._kmers.begin(), vec._kmers.end());
				}

				
				if(12424 == nodeName){
					for(u_int64_t minimizer : vec._kmers){
						cout << minimizer << " ";
					}
					cout << endl;
					//exit(1);
				}
				


				//_debugMinimizers[nodeName] = vec._kmers;
				
				auto found = _nodeName_entire.find(nodeName);
				if(found != _nodeName_entire.end() && !found->second._loaded){
					_nodeName_entire[nodeName] = {vec._kmers, true};
				}
				
				auto found2 = _nodeName_left.find(nodeName);
				if(_nodeName_left.find(nodeName) != _nodeName_left.end() && !found2->second._loaded){
					//std::reverse(vec._kmers.begin(), vec._kmers.end());
					//if(12424 == nodeName){
					//	cout << "left: " <<  vec._kmers[0] << " " << vec._kmers[vec._kmers.size()-1] << endl;
					//}
					if(kminmerInfo._isReversed){
						_nodeName_left[nodeName] = {vec._kmers[vec._kmers.size()-1], true};
					}
					else{
						_nodeName_left[nodeName] = {vec._kmers[0], true};
					}
				}
				
				found2 = _nodeName_right.find(nodeName);
				if(_nodeName_right.find(nodeName) != _nodeName_right.end() && !found2->second._loaded){
					if(kminmerInfo._isReversed){
						_nodeName_right[nodeName] = {vec._kmers[0], true};
					}
					else{
						_nodeName_right[nodeName] = {vec._kmers[vec._kmers.size()-1], true};
					}
				}

				if(12424 == nodeName){
					//exit(1);
				}

			}


			
			readIndex += 1;
		}



	}

	void createMinimizerContigs(){


		string mdbg_filename = "/home/gats/workspace/run/overlap_test_562_k3/mdbg_nodes.gz";
		MDBG* mdbg = new MDBG(3);
		mdbg->load(mdbg_filename);

		ofstream testCsv(_inputDir + "/minimizerContigsDebug.csv");
		testCsv << "Name,Colour" << endl;

		gzFile outputContigFile = gzopen(_filename_outputContigs.c_str(),"wb");


		cout << endl;
		cout << "Creating basespace contigs" << endl;
		cout << endl;

		gzFile contigFile = gzopen(_filename_nodeContigs.c_str(), "rb");

		u_int64_t contig_index = 0;

		while(true){

			vector<u_int32_t> nodePath;
			vector<u_int64_t> supportingReads;
			u_int64_t size;
			gzread(contigFile, (char*)&size, sizeof(size));
			

			if(gzeof(contigFile)) break;

			nodePath.resize(size);
			supportingReads.resize(size);
			gzread(contigFile, (char*)&nodePath[0], size * sizeof(u_int32_t));
			gzread(contigFile, (char*)&supportingReads[0], size * sizeof(u_int64_t));


			vector<u_int64_t> contigSequence;

			for(size_t i=0; i<nodePath.size(); i++){
				
				u_int32_t nodeIndex = nodePath[i];
				bool orientation;
				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex, orientation);
				u_int64_t readIndex = supportingReads[i];
				ContigNode contigNode = {nodeName, readIndex};

				/*
				cout << "Next node: " << nodeName << endl;
			
				if(_debugMinimizers.find(nodeName) != _debugMinimizers.end()){
					cout << nodeName << endl;
					for(u_int64_t minimizer : _debugMinimizers[nodeName]){
						cout << minimizer << " ";
					}
					cout << endl;
				}*/
				

				if(i == 0){
					
					if(orientation){ //+
						cout << "Entire" << endl;

						vector<u_int64_t> minimizers = _nodeName_entire[nodeName]._minimizers;
						for(u_int64_t minimizer : minimizers) contigSequence.push_back(minimizer);
					}
					else{
						cout << "Entire RC" << endl;
						vector<u_int64_t> minimizers = _nodeName_entire[nodeName]._minimizers;
						std::reverse(minimizers.begin(), minimizers.end());
						for(u_int64_t minimizer : minimizers) contigSequence.push_back(minimizer);
					}
				}
				else {
					if(orientation){
						u_int64_t minimizer = _nodeName_right[nodeName]._minimizer;
						contigSequence.push_back(minimizer);
						//cout << "Add minimizer (right): " << minimizer << endl;
					}
					else{
						u_int64_t minimizer = _nodeName_left[nodeName]._minimizer;
						contigSequence.push_back(minimizer);
						//cout << "Add minimizer (left): " << minimizer << endl;
					}
				}

				
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

			//exit(1);
			//cout << contigSequence.substr(362170, 100) << endl;
			cout << endl << endl;
			cout << nodePath.size() << endl;
			cout << contigSequence.size() << endl;
			//cout << contigSequence << endl;


			u_int64_t contigSize = contigSequence.size();
			gzwrite(outputContigFile, (const char*)&contigSize, sizeof(contigSize));
			gzwrite(outputContigFile, (const char*)&contigSequence[0], contigSize * sizeof(u_int64_t));




			vector<KmerVec> kminmers; 
			vector<ReadKminmer> kminmersInfo;
			vector<u_int64_t> minimizersPos; 
			vector<u_int64_t> rlePositions;
			MDBG::getKminmers(21, 3, contigSequence, minimizersPos, kminmers, kminmersInfo, rlePositions, -1);

			cout << "Nb kminmers: " << kminmers.size() << endl;
			u_int64_t nbFoundMinimizers = 0;
			for(size_t i=0; i<kminmers.size(); i++){
				KmerVec& vec = kminmers[i];
				
				if(mdbg->_dbg_nodes.find(vec) == mdbg->_dbg_nodes.end()) continue;

				u_int32_t nodeName = mdbg->_dbg_nodes[vec]._index;
				//cout << vec._kmers[0] << " " << vec._kmers[1] << " " << vec._kmers[2] << "     " << nodeName << endl;
				//cout << nodeName << endl;

				nbFoundMinimizers += 1;
				//ReadKminmer& kminmerInfo = kminmersInfo[i];

				testCsv << nodeName << "," << contig_index << endl;

			}

			cout << "Found nodes: " << nbFoundMinimizers << endl;










			contig_index += 1;
		}

		gzclose(contigFile);
		gzclose(outputContigFile);
		testCsv.close();

	}
};	


#endif 



