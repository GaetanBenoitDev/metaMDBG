

#ifndef MDBG_METAG_MAPPINGBINREADS
#define MDBG_METAG_MAPPINGBINREADS

#include "../Commons.hpp"

class Mapping_BinReads : public Tool{
    
public:

	string _contigFilename;
	string _readFilename;
	string _outputDir;

	MinimizerParser* _minimizerParser;
	size_t _minimizerSize;
	size_t _kminmerSize;
	float _minimizerDensity;

	Mapping_BinReads(): Tool (){

	}


	void parseArgs(int argc, char* argv[]){

		_kminmerSize = 4;
		_minimizerSize = 21;
		_minimizerDensity = 0.05;

		cxxopts::Options options("ToBasespace", "");
		options.add_options()
		//(ARG_INPUT_FILENAME, "", cxxopts::value<string>())
		//(ARG_INPUT_FILENAME_CONTIG, "", cxxopts::value<string>())
		//(ARG_INPUT_FILENAME_TRUTH, "", cxxopts::value<string>()->default_value(""))
		//(ARG_OUTPUT_FILENAME, "", cxxopts::value<string>())
		//(ARG_FIRST_PASS, "", cxxopts::value<bool>()->default_value("false"))
		//(ARG_FASTA, "", cxxopts::value<bool>()->default_value("false"))
		//(ARG_OUTPUT_DIR, "", cxxopts::value<string>())
		//(ARG_NB_CORES, "", cxxopts::value<int>()->default_value(NB_CORES_DEFAULT));
		("contig", "", cxxopts::value<string>())
		("read", "", cxxopts::value<string>())
		("outputDir", "", cxxopts::value<string>());
		//("outputDir", "", cxxopts::value<string>())
		//(ARG_INPUT_FILENAME_TRUTH, "", cxxopts::value<string>()->default_value(""))
		//(ARG_EVAL, "", cxxopts::value<bool>()->default_value("false"))
		//(ARG_INPUT_FILENAME_ABUNDANCE, "", cxxopts::value<string>()->default_value(""))
		//(ARG_NB_CORES, "", cxxopts::value<int>()->default_value("8"));
		//(ARG_KMINMER_LENGTH, "", cxxopts::value<int>()->default_value("3"))
		//(ARG_MINIMIZER_LENGTH, "", cxxopts::value<int>()->default_value("21"))
		//(ARG_MINIMIZER_DENSITY, "", cxxopts::value<float>()->default_value("0.005"));
		//("k,kminmerSize", "File name", cxxopts::value<std::string>())
		//("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
		//;
		options.parse_positional({"contig", "read", "outputDir"});
		options.positional_help("contigs readFilename outputDir");


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

			_contigFilename = result["contig"].as<string>();
			_readFilename = result["read"].as<string>();
			_outputDir = result["outputDir"].as<string>();
			
		}
		catch (const std::exception& e){
			std::cout << options.help() << std::endl;
			std::cerr << e.what() << std::endl;
			std::exit(EXIT_FAILURE);
		}

		/*
		_inputFilename = getInput()->getStr(STR_INPUT);
		_inputDir = getInput()->getStr(STR_INPUT_DIR);
		*/
		/*
		string filename_parameters = _inputDir + "/parameters.gz";
		gzFile file_parameters = gzopen(filename_parameters.c_str(),"rb");
		gzread(file_parameters, (char*)&_minimizerSize, sizeof(_minimizerSize));
		gzread(file_parameters, (char*)&_kminmerSize, sizeof(_kminmerSize));
		gzread(file_parameters, (char*)&_minimizerDensity, sizeof(_minimizerDensity));
		gzread(file_parameters, (char*)&_kminmerSizeFirst, sizeof(_kminmerSizeFirst));
		gzclose(file_parameters);

		_kminmerSize = _kminmerSizeFirst;

		cout << endl;
		cout << "Input dir: " << _inputDir << endl;
		//cout << "Output filename: " << _outputFilename << endl;
		cout << "Minimizer length: " << _minimizerSize << endl;
		cout << "Kminmer length: " << _kminmerSize << endl;
		cout << "Density: " << _minimizerDensity << endl;
		cout << endl;

		//_filename_outputContigs = _inputDir + "/contigs.fasta.gz";
		//_filename_outputContigs = _inputFilenameContig + ".fasta.gz"; //_inputDir + "/tmpContigs.fasta.gz";
		*/
		_minimizerParser = new MinimizerParser(_minimizerSize, _minimizerDensity);
	}


    void execute (){
		extract_truth_kminmers();
		//binReads();
	}


	unordered_map<KmerVec, vector<u_int32_t>> _kmervec_to_contigIndex;
	unordered_map<u_int32_t, ofstream> _binFiles;
	u_int32_t _binIndex;

	//MinimizerParser* _minimizerParser;
	//u_int32_t _extract_truth_kminmers_read_position;
	//unordered_map<u_int32_t, vector<string>> _evaluation_hifiasmGroundTruth_nodeName_to_unitigName;
	//vector<u_int32_t> _evaluation_hifiasmGroundTruth_path;
	//unordered_set<u_int32_t> _hifiasm_startingNodenames;
	//ofstream _file_groundTruth_hifiasm_position;
	EncoderRLE _encoderRLE;

	void extract_truth_kminmers(){

		//_file_groundTruth_hifiasm_position.open(_inputDir + "/groundtruth_hifiasm_position.csv");
		//_file_groundTruth_hifiasm_position << "Name,Position" << endl;

		//_extract_truth_kminmers_read_position = 0;
		_minimizerParser = new MinimizerParser(_minimizerSize, _minimizerDensity);
		_binIndex = 0;

		auto fp = std::bind(&Mapping_BinReads::extract_truth_kminmers_read, this, std::placeholders::_1);
		ReadParser readParser(_contigFilename, true, false);
		readParser.parse(fp);


		//_file_groundTruth_hifiasm_position.close();

		//delete _minimizerParser;
	}

	void extract_truth_kminmers_read(const Read& read){
		//ottalSize += strlen(read->seq.s);

		if(read._seq.size() < 2000000) return;

		

		u_int64_t readIndex = read._index;

		string rleSequence;
		vector<u_int64_t> rlePositions;
		_encoderRLE.execute(read._seq.c_str(), read._seq.size(), rleSequence, rlePositions);

		vector<u_int64_t> minimizers;
		vector<u_int64_t> minimizers_pos;
		_minimizerParser->parse(rleSequence, minimizers, minimizers_pos);

		vector<KmerVec> kminmers; 
		vector<ReadKminmer> kminmersInfo;
		MDBG::getKminmers(_minimizerSize, _kminmerSize, minimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, readIndex, false);

		for(size_t i=0; i<kminmers.size(); i++){

			KmerVec& vec = kminmers[i];
			_kmervec_to_contigIndex[vec].push_back(_binIndex);
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

		string referenceFilename = _outputDir + "/reads_" + to_string(_binIndex) + "_reference.fasta.gz"; 
		gzFile referenceFile = gzopen(referenceFilename.c_str(),"wb");

		string header = ">ctg" + to_string(_binIndex) + '\n';
		gzwrite(referenceFile, (const char*)&header[0], header.size());
		string seq = read._seq + '\n';
		gzwrite(referenceFile, (const char*)&seq[0], seq.size());
		gzclose(referenceFile);
		
		string binFilename = _outputDir + "/reads_" + to_string(_binIndex) + ".fastq"; 
		_binFiles[_binIndex] = ofstream(binFilename);

		string inputFilename = _outputDir + "/input_reads_" + to_string(_binIndex) + ".txt"; 
		ofstream inputFile(inputFilename);
		inputFile << binFilename << endl;
		inputFile.close();

		_binIndex += 1;
	}

	void binReads(){

		auto fp = std::bind(&Mapping_BinReads::binReads_read, this, std::placeholders::_1);
		ReadParser readParser(_readFilename, true, false);
		readParser.parse(fp);

		for(auto& it : _binFiles){
			it.second.close();
		}
	}

	void binReads_read(const Read& read){
		//ottalSize += strlen(read->seq.s);


		u_int64_t readIndex = read._index;
		//cout << readIndex << " " << _kmervec_to_contigIndex.size() << endl;

		string rleSequence;
		vector<u_int64_t> rlePositions;
		_encoderRLE.execute(read._seq.c_str(), read._seq.size(), rleSequence, rlePositions);

		vector<u_int64_t> minimizers;
		vector<u_int64_t> minimizers_pos;
		_minimizerParser->parse(rleSequence, minimizers, minimizers_pos);

		vector<KmerVec> kminmers; 
		vector<ReadKminmer> kminmersInfo;
		MDBG::getKminmers(_minimizerSize, _kminmerSize, minimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, readIndex, false);

		unordered_set<u_int32_t> isWrittenBinIndex;

		for(size_t i=0; i<kminmers.size(); i++){

			KmerVec& vec = kminmers[i];
			if(_kmervec_to_contigIndex.find(vec) == _kmervec_to_contigIndex.end()) continue;

			for(u_int32_t binIndex : _kmervec_to_contigIndex[vec]){
				
				if(isWrittenBinIndex.find(binIndex) != isWrittenBinIndex.end()) continue;
				isWrittenBinIndex.insert(binIndex);

				//cout << "lala: " << binIndex << endl;
				ofstream& file = _binFiles[binIndex];
				//cout << file.is_open() << " " << read._header << endl;
				file << "@" << read._header << endl;
				file << read._seq << endl;
				file << "+" << endl;
				file << read._qual << endl;

				//file.flush();
			}

			/*
			if(_kmervec_to_unitigName.find(vec) == _kmervec_to_unitigName.end()) continue;

			string unitigName = _kmervec_to_unitigName[vec];

			if(_writtenUnitigNames.find(unitigName) != _writtenUnitigNames.end()) continue;

			_outputFile << unitigName << "," << _currentBinIndex << endl;
			_writtenUnitigNames.insert(unitigName);
			*/
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

};	


#endif 


