

#ifndef MDBG_METAG_MAPPINGCONTIGTOMDBG
#define MDBG_METAG_MAPPINGCONTIGTOMDBG

#include "../Commons.hpp"

class Mapping_ContigToMDBG : public Tool{
    
public:

	string _mdbgDir;
	string _binDir;
	string _outputFilename;

	MinimizerParser* _minimizerParser;
	size_t _minimizerSize;
	size_t _kminmerSize;
	size_t _kminmerSizeFirst;
	float _minimizerDensity;

	MDBG* _mdbg;

	Mapping_ContigToMDBG(): Tool (){

	}


	void parseArgs(int argc, char* argv[]){

		//_kminmerSize = 4;
		//_minimizerSize = 21;
		//_minimizerDensity = 0.05;

		cxxopts::Options options("ToBasespace", "");
		options.add_options()
		("mdbgDir", "", cxxopts::value<string>())
		("binDir", "", cxxopts::value<string>())
		("outputFilename", "", cxxopts::value<string>());

		options.parse_positional({"mdbgDir", "binDir", "outputFilename"});
		options.positional_help("mdbgDir binDir outputFilename");



		if(argc <= 1){
			cout << options.help() << endl;
			exit(0);
		}

		cxxopts::ParseResult result;

		try{
			result = options.parse(argc, argv);

			_mdbgDir = result["mdbgDir"].as<string>();
			_binDir = result["binDir"].as<string>();
			_outputFilename = result["outputFilename"].as<string>();
		}
		catch (const std::exception& e){
			std::cout << options.help() << std::endl;
			std::cerr << e.what() << std::endl;
			std::exit(EXIT_FAILURE);
		}

		string filename_parameters = _mdbgDir + "/parameters.gz";
		gzFile file_parameters = gzopen(filename_parameters.c_str(),"rb");
		gzread(file_parameters, (char*)&_minimizerSize, sizeof(_minimizerSize));
		gzread(file_parameters, (char*)&_kminmerSize, sizeof(_kminmerSize));
		gzread(file_parameters, (char*)&_minimizerDensity, sizeof(_minimizerDensity));
		gzclose(file_parameters);

		//ifstream file_data(_mdbgDir + "/data.txt");
		//file_data.read((char*)&_nbReads, sizeof(_nbReads));
		//file_data.close();

		cout << endl;
		cout << "Input dir: " << _mdbgDir << endl;
		//cout << "Output filename: " << _outputFilename << endl;
		cout << "Minimizer length: " << _minimizerSize << endl;
		cout << "Kminmer length: " << _kminmerSize << endl;
		cout << "Density: " << _minimizerDensity << endl;
		//cout << "Nb reads: " << _nbReads << endl;
		cout << endl;

		_minimizerSize = 17;
		_minimizerDensity = 0.01;
		_kminmerSizeFirst = 4;

		_minimizerParser = new MinimizerParser(_minimizerSize, _minimizerDensity);
	}


    void execute (){

		string mdbg_filename = _mdbgDir + "/mdbg_nodes.gz";

		//cout << _gfaFilename << endl;
		_mdbg = new MDBG(_kminmerSize);
		_mdbg->load(mdbg_filename);
		cout << "Nb nodes: " <<  _mdbg->_dbg_nodes.size() << endl;

		indexNodeNames();
		extract_truth_kminmers();
		mapToMdbg();
		assingBinToUnitigs();
		//binReads();
	}

	unordered_map<KmerVec, vector<u_int32_t>> _kmerVec_to_binIndexes;
	unordered_map<u_int32_t, vector<u_int32_t>> _nodeName_to_unitigIndexes;
	//unordered_map<KmerVec, string> _kmervec_to_binIndex;

	EncoderRLE _encoderRLE;
	u_int32_t _currentBinIndex;

	void indexNodeNames(){

		ifstream file(_mdbgDir + "/nodeName_to_unitigIndex.bin");

		while(true){
			
			u_int32_t nodeName;
			u_int32_t unitigIndex;

			file.read((char*)&nodeName, sizeof(nodeName));
			if(file.eof()) break;
			file.read((char*)&unitigIndex, sizeof(unitigIndex));

			_nodeName_to_unitigIndexes[nodeName].push_back(unitigIndex);
		}

		file.close();
	}

	void extract_truth_kminmers(){

		_currentBinIndex = 0;
		//_outputFile = ofstream(_outputFilename);
		//_outputFile << "Name,Color" << endl;

		for (const auto & p : fs::directory_iterator(_binDir)){
			string ext = p.path().extension();
			//cout << p.path() << endl;

			//cout << ext << endl;
			if(ext == ".fa" || ext == ".fasta" || ext == ".fna"){
				string filename = p.path();
				cout << filename << endl;

				//string binName = p.path().filename();
				//binName.erase(binName.find("bin."), 4);
				//binName.erase(binName.find(ext), ext.size());
				//cout << binName << endl;
				//_currentBinIndex = stoull(binName);


				extract_truth_kminmers_bin(filename);

				_currentBinIndex += 1;
			}

		}

		//_outputFile.close();
	}


	void extract_truth_kminmers_bin(const string& binFilename){

		auto fp = std::bind(&Mapping_ContigToMDBG::extract_truth_kminmers_bin_read, this, std::placeholders::_1);
		ReadParser readParser(binFilename, true, false);
		readParser.parse(fp);
	}

	void extract_truth_kminmers_bin_read(const Read& read){
		//ottalSize += strlen(read->seq.s);


		u_int64_t readIndex = read._index;

		//cout << read._seq.size() << endl;
		string rleSequence;
		vector<u_int64_t> rlePositions;
		_encoderRLE.execute(read._seq.c_str(), read._seq.size(), rleSequence, rlePositions);

		vector<u_int64_t> minimizers;
		vector<u_int64_t> minimizers_pos;
		_minimizerParser->parse(rleSequence, minimizers, minimizers_pos);

		vector<KmerVec> kminmers; 
		vector<ReadKminmer> kminmersInfo;
		MDBG::getKminmers(_minimizerSize, _kminmerSizeFirst, minimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, readIndex, false);

		for(size_t i=0; i<kminmers.size(); i++){


			//cout << i << endl;
			KmerVec& vec = kminmers[i];
			_kmerVec_to_binIndexes[vec].push_back(_currentBinIndex);
			//cout << _currentBinIndex << endl;
			/*
			cout << "oue?" << endl;
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()) continue;

			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
			cout << nodeName << endl;

			if(_nodeName_to_unitigIndexes.find(nodeName) == _nodeName_to_unitigIndexes.end()) continue;
			*/

			//if(_kmervec_to_unitigName.find(vec) == _kmervec_to_unitigName.end()) continue;

			//string unitigName = _kmervec_to_unitigName[vec];

			//if(_writtenUnitigNames.find(unitigName) != _writtenUnitigNames.end()) continue;

			//_outputFile << unitigName << "," << _currentBinIndex << endl;
			//_writtenUnitigNames.insert(unitigName);
			//_kmervec_to_unitigName[vec] = read._header;
			/*
			if(_mdbg->_dbg_nodes.find(vec) != _mdbg->_dbg_nodes.end()){

				u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
				

				_evaluation_hifiasmGroundTruth_nodeName_to_unitigName[_mdbg->_dbg_nodes[vec]._index].push_back(read._header);
				_evaluation_hifiasmGroundTruth_path.push_back(_mdbg->_dbg_nodes[vec]._index);


			}
			else{
				_evaluation_hifiasmGroundTruth_path.push_back(-1);
			}
			*/


		}


	}

	void mapToMdbg(){

		for(const auto& it : _mdbg->_dbg_nodes){

			const KmerVec& vec = it.first;
			u_int32_t nodeName = it.second._index;

			//cout << nodeName << endl;
			if(_nodeName_to_unitigIndexes.find(nodeName) == _nodeName_to_unitigIndexes.end()) continue;


			vector<u_int64_t> rlePositions;
			vector<u_int64_t> minimizers_pos;//(minimizers.size());
			vector<KmerVec> kminmers; 
			vector<ReadKminmer> kminmersInfo;
			MDBG::getKminmers(_minimizerSize, _kminmerSizeFirst, vec._kmers, minimizers_pos, kminmers, kminmersInfo, rlePositions, 0, false);

			for(const KmerVec& vec : kminmers){

				//cout << "?" << endl;
				if(_kmerVec_to_binIndexes.find(vec) == _kmerVec_to_binIndexes.end()) continue;

				const vector<u_int32_t>& binIndexes = _kmerVec_to_binIndexes[vec];

				for(u_int32_t unitigIndex : _nodeName_to_unitigIndexes[nodeName]){
					for(u_int32_t binIndex : binIndexes){
						_unitigIndex_to_binIndexes[unitigIndex].push_back(binIndex);
					}
				}
			}
			
		}
	}
	
	unordered_map<u_int32_t, vector<u_int32_t>> _unitigIndex_to_binIndexes;
	//unordered_map<KmerVec, vector<u_int32_t>> _kmerVec_to_binIndexes;
	//unordered_map<u_int32_t, vector<u_int32_t>> _nodeName_to_unitigIndexes;

	void assingBinToUnitigs(){

		ofstream file(_outputFilename);
		file << "Name,Color" << endl;

		for(const auto& it : _unitigIndex_to_binIndexes){
			u_int32_t unitigIndex = it.first;
			const vector<u_int32_t>& binIndexes = it.second;

			unordered_map<u_int32_t, u_int32_t> binCounts;

			for(u_int32_t binIndex : binIndexes){
				binCounts[binIndex] += 1;
			}

			u_int32_t maxBinIdnex = -1;
			u_int32_t maxBinCount = 0;
			for(const auto& it : binCounts){
				u_int32_t binIndex = it.first;
				u_int32_t binCount = it.second;

				if(binCount > maxBinCount){
					maxBinCount = binCount;
					maxBinIdnex = binIndex;
				}
			}

			file << unitigIndex << "," << maxBinIdnex << endl;
		}

		file.close();
	}
};	


#endif 


