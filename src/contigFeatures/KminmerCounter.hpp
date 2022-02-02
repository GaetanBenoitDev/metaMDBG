

#ifndef MDBG_METAG_KMINMERCOUNTER
#define MDBG_METAG_KMINMERCOUNTER

#include "../Commons.hpp"


class KminmerCounter : public Tool{
    
public:

	string _inputFilename;
	string _inputDir;
	string _inputFilename_contig;
	float _minimizerDensity;
    size_t _minimizerSize;
    size_t _kminmerSize;

	MDBG* _mdbg;
	u_int64_t _nbDatasets;
	unordered_map<u_int32_t, vector<u_int32_t>> _kminmerCounts;

	KminmerCounter(): Tool (){

	}

    void execute (){
    
		cout << "Loading mdbg" << endl;
		string mdbg_filename = _inputDir + "/mdbg_nodes.gz";
		_mdbg = new MDBG(_kminmerSize);
		_mdbg->load(mdbg_filename);
		cout << "MDBG nodes: " << _mdbg->_dbg_nodes.size() << endl;

		extractContigKminmers();
		countKminmers();
		computeContigCoverage();
		
		delete _mdbg;

	}

	void parseArgs(int argc, char* argv[]){


		cxxopts::Options options("ToMinspace", "");
		options.add_options()
		//(ARG_INPUT_FILENAME, "", cxxopts::value<string>())
		(ARG_OUTPUT_DIR, "", cxxopts::value<string>())
		(ARG_INPUT_FILENAME, "", cxxopts::value<string>()->default_value(""))
		(ARG_INPUT_FILENAME_CONTIG, "", cxxopts::value<string>()->default_value(""));


		if(argc <= 1){
			cout << options.help() << endl;
			exit(0);
		}

		cxxopts::ParseResult result;

		try{
			result = options.parse(argc, argv);

			_inputFilename = result[ARG_INPUT_FILENAME].as<string>();
			_inputDir = result[ARG_OUTPUT_DIR].as<string>();
			_inputFilename_contig = result[ARG_INPUT_FILENAME_CONTIG].as<string>();
			
			
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
		cout << "Input dir: " << _inputDir << endl;
		//cout << "Output filename: " << _outputFilename << endl;
		cout << "Minimizer length: " << _minimizerSize << endl;
		cout << "Kminmer length: " << _kminmerSize << endl;
		cout << "Density: " << _minimizerDensity << endl;
		cout << endl;

	}

	void extractContigKminmers (){

		ReadParser parserDatasets(_inputFilename, false, _minimizerSize, _kminmerSize, _minimizerDensity);
		_nbDatasets = parserDatasets.getNbDatasets();
		cout << "Nb datasets: " << _nbDatasets << endl;

		ReadParser parser(_inputFilename_contig, true, _minimizerSize, _kminmerSize, _minimizerDensity);
		auto fp = std::bind(&KminmerCounter::extractContigKminmers_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
		parser.parseKminmers(fp);

	}


	void extractContigKminmers_read(const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex, u_int64_t datasetIndex){

		vector<u_int32_t> counts;
		for(u_int64_t i=0; i<_nbDatasets; i++) counts.push_back(0);

		for(size_t i=0; i<kminmersInfos.size(); i++){
			

			KmerVec vec = kminmers[i];
			//for(KmerVec& vec : kminmers){
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()) continue;

			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
			_kminmerCounts[nodeName] = counts;
			
		}

	}


	void countKminmers (){

		ReadParser parser(_inputFilename, false, _minimizerSize, _kminmerSize, _minimizerDensity);
		auto fp = std::bind(&KminmerCounter::countKminmers_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
		parser.parseKminmers(fp);

		/*
		for(const auto& it : _kminmerCounts){

			int sum = 0;
			for(u_int32_t count : _kminmerCounts[it.first]){
				sum += count;
			}

			if(sum == 0) continue;

			cout << it.first << ": ";
			for(u_int32_t count : _kminmerCounts[it.first]){
				cout << count << " ";
			}
			cout << endl;
		}
		*/
	}


	void countKminmers_read(const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex, u_int64_t datasetIndex){

		if(readIndex % 100000 == 0){
			cout << datasetIndex << " " << readIndex << endl;
		}
		//cout << datasetIndex << " " << readIndex << endl;
		//cout << "------" << endl;

		//cout << readIndex << endl;
		for(size_t i=0; i<kminmersInfos.size(); i++){
			
			//const ReadKminmer& kminmerInfo = kminmersInfos[i];

			KmerVec vec = kminmers[i];
			//for(KmerVec& vec : kminmers){
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()) continue;

			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;

			if(_kminmerCounts.find(nodeName) == _kminmerCounts.end()) continue;

			_kminmerCounts[nodeName][datasetIndex] += 1;
		}

	}



	void computeContigCoverage (){

		ReadParser parser(_inputFilename_contig, true, _minimizerSize, _kminmerSize, _minimizerDensity);
		auto fp = std::bind(&KminmerCounter::computeContigCoverage_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
		parser.parseKminmers(fp);
	}


	void computeContigCoverage_read(const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex, u_int64_t datasetIndex){


		vector<vector<u_int32_t>> abundances;
		for(u_int64_t i=0; i<_nbDatasets; i++) abundances.push_back({});
		
		//cout << datasetIndex << " " << readIndex << endl;
		//cout << "------" << endl;

		//cout << readIndex << endl;
		for(size_t i=0; i<kminmersInfos.size(); i++){
			
			//const ReadKminmer& kminmerInfo = kminmersInfos[i];

			KmerVec vec = kminmers[i];
			//for(KmerVec& vec : kminmers){
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()) continue;

			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;

			if(_kminmerCounts.find(nodeName) == _kminmerCounts.end()) continue;

			const vector<u_int32_t>& counts = _kminmerCounts[nodeName];

			for(u_int64_t i=0; i<counts.size(); i++){
				abundances[i].push_back(counts[i]);
			}
		}

		cout << "Ctg" << readIndex << ": ";
		for(u_int64_t i=0; i<_nbDatasets; i++){
			//cout << abundances[i].size() << endl;
			float median = Utils::compute_median(abundances[i]);
			cout << median << " ";
		}
		cout << endl;

	}




};	


#endif 



