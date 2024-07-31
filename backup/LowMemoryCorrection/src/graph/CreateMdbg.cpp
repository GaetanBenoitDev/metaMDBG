

#include "CreateMdbg.hpp"

CreateMdbg::CreateMdbg () : Tool()
{

}


void CreateMdbg::parseArgs(int argc, char* argv[]){

	
	args::ArgumentParser parser("", ""); //"This is a test program.", "This goes after the options."
	args::Positional<std::string> arg_outputDir(parser, "outputDir", "Output dir", args::Options::Required);
	//args::PositionalList<std::string> arg_readFilenames(parser, "reads", "Input filename(s) (separated by space)", args::Options::Required);
	//args::ValueFlag<int> arg_l(parser, "", "Minimizer length", {ARG_MINIMIZER_LENGTH2}, 13);
	//args::ValueFlag<float> arg_d(parser, "", "Minimizer density", {ARG_MINIMIZER_DENSITY2}, 0.005f);
	args::ValueFlag<int> arg_nbCores(parser, "", "Number of cores", {ARG_NB_CORES2}, NB_CORES_DEFAULT_INT);
	args::Flag arg_bf(parser, "", "Filter unique erroneous k-minmers", {ARG_BLOOM_FILTER});
	args::Flag arg_firstPass(parser, "", "Is first pass of multi-k", {ARG_FIRST_PASS});
	args::Flag arg_help(parser, "", "", {'h', "help"}, args::Options::Hidden);
	//args::HelpFlag help(parser, "help", "Display this help menu", {'h'});
	//args::CompletionFlag completion(parser, {"complete"});

	//(ARG_INPUT_FILENAME_TRUTH, "", cxxopts::value<string>()->default_value(""))
	//(ARG_MINIMIZER_LENGTH, "", cxxopts::value<int>()->default_value("13"))
	//(ARG_MINIMIZER_DENSITY, "", cxxopts::value<float>()->default_value("0.005"))
	//(ARG_NB_CORES, "", cxxopts::value<int>()->default_value(NB_CORES_DEFAULT))
	//(ARG_BLOOM_FILTER, "", cxxopts::value<bool>()->default_value("false"));

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

	_outputDir = args::get(arg_outputDir);
	_nbCores = args::get(arg_nbCores);

	_useBloomFilter = true;
	if(arg_bf){
		_useBloomFilter = false;
	}

	_isFirstPass = false;
	if(arg_firstPass){
		_isFirstPass = true;
	}

	if(!_isFirstPass){
		_useBloomFilter = false;	
	}

	/*
	cxxopts::Options options("Graph", "Create MDBG");
	options.add_options()
	//("d,debug", "Enable debugging") // a bool parameter
	(ARG_INPUT_FILENAME, "", cxxopts::value<string>())
	(ARG_OUTPUT_DIR, "", cxxopts::value<string>())
	(ARG_FIRST_PASS, "", cxxopts::value<bool>()->default_value("false"))
	(ARG_INPUT_FILENAME_CONTIG, "", cxxopts::value<string>()->default_value(""))
	(ARG_NB_CORES, "", cxxopts::value<int>()->default_value(NB_CORES_DEFAULT))
	(ARG_BLOOM_FILTER, "", cxxopts::value<bool>()->default_value("false"));
	//(ARG_KMINMER_LENGTH, "", cxxopts::value<int>()->default_value("3"))
	//(ARG_MINIMIZER_LENGTH, "", cxxopts::value<int>()->default_value("21"))
	//(ARG_MINIMIZER_DENSITY, "", cxxopts::value<float>()->default_value("0.005"));
	//("k,kminmerSize", "File name", cxxopts::value<std::string>())
	//("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
	//;

	if(argc <= 1){
		std::cout << options.help() << std::endl;
		exit(0);
	}


	cxxopts::ParseResult result;

    try{
		result = options.parse(argc, argv);

		_inputFilename = result[ARG_INPUT_FILENAME].as<string>();
		_outputDir = result[ARG_OUTPUT_DIR].as<string>(); 
		_filename_inputContigs = result[ARG_INPUT_FILENAME_CONTIG].as<string>();
		_isFirstPass = result[ARG_FIRST_PASS].as<bool>();
		_nbCores = result[ARG_NB_CORES].as<int>();
		_useBloomFilter = result[ARG_BLOOM_FILTER].as<bool>();

    }
    catch (const std::exception& e){
		std::cout << options.help() << std::endl;
    	std::cerr << e.what() << std::endl;
    	std::exit(EXIT_FAILURE);
    }
	*/


	
	string filename_parameters = _outputDir + "/parameters.gz";
	gzFile file_parameters = gzopen(filename_parameters.c_str(),"rb");
	gzread(file_parameters, (char*)&_minimizerSize, sizeof(_minimizerSize));
	gzread(file_parameters, (char*)&_kminmerSize, sizeof(_kminmerSize));
	gzread(file_parameters, (char*)&_minimizerDensity, sizeof(_minimizerDensity));
	gzread(file_parameters, (char*)&_kminmerSizeFirst, sizeof(_kminmerSizeFirst));
	gzread(file_parameters, (char*)&_minimizerSpacingMean, sizeof(_minimizerSpacingMean));
	gzread(file_parameters, (char*)&_kminmerLengthMean, sizeof(_kminmerLengthMean));
	gzread(file_parameters, (char*)&_kminmerOverlapMean, sizeof(_kminmerOverlapMean));
	gzread(file_parameters, (char*)&_kminmerSizePrev, sizeof(_kminmerSizePrev));
	gzread(file_parameters, (char*)&_kminmerSizeLast, sizeof(_kminmerSizeLast));
	gzclose(file_parameters);

	openLogFile(_outputDir);


	_logFile << endl;
	_logFile << "Contig filename: " << _filename_inputContigs << endl;
	_logFile << "Output dir: " << _outputDir << endl;
	_logFile << "Minimizer length: " << _minimizerSize << endl;
	_logFile << "Kminmer length: " << _kminmerSize << endl;
	_logFile << "Density: " << _minimizerDensity << endl;
	_logFile << endl;

	//_inputDir = getInput()->get(STR_INPUT_DIR) ? getInput()->getStr(STR_INPUT_DIR) : "";
	//_input_extractKminmers= getInput()->get(STR_INPUT_EXTRACT_KMINMERS) ? getInput()->getStr(STR_INPUT_EXTRACT_KMINMERS) : "";

	_filename_readMinimizers = _outputDir + "/read_data.gz";
	_filename_contigMinimizers = _outputDir + "/contig_data.gz";
	//_filename_outputContigs = _inputDir + "/contigs.min.gz";
	//_filename_readCompositions = _outputDir + "/read_compositions.gz";
	//_filename_filteredMinimizers = _outputDir + "/filteredMinimizers.gz";
	//_filename_hifiasmGroundtruth = _outputDir + "/hifiasmGroundtruth.gz";



	//_filename_noKminmerReads = _outputDir + "/" + FILENAME_NO_KMINMER_READS;
	//cout << _minimizerSpacingMean << endl;
	//cout << _kminmerLengthMean << endl;
	//cout << _kminmerOverlapMean << endl;
	//getchar();
}

void CreateMdbg::execute (){

	//_file_noKminmerReads = ofstream(_filename_noKminmerReads, std::ofstream::app);
	_node_id = 0;
	//parseArgs();
	createMDBG();
	createGfa();


	if(_isFirstPass){
		const auto copyOptions = fs::copy_options::overwrite_existing;
		fs::copy(_outputDir + "/kminmerData_min.txt", _outputDir + "/kminmerData_min_init.txt", copyOptions);
		//_mdbg->dump(_outputDir + "/mdbg_nodes_init.gz");

		/*
		_mdbg = new MDBG(_kminmerSize);
		_mdbg->load(_outputDir + "/kminmerData_min_init.txt", false);

		KmerVec vec;
		vec._kmers = {7892019983131394, 62446187524336769, 47730431539434127, 90993437045220547};
		cout << _mdbg->_dbg_nodes[vec]._index << endl;
		getchar();
		*/

	}
	


	if(_isFirstPass){
		ofstream file_data(_outputDir + "/data.txt");

		file_data.write((const char*)&_nbReads, sizeof(_nbReads));

		file_data.close();
	}

	//_file_noKminmerReads.close();

	closeLogFile();
}


void CreateMdbg::createMDBG (){

	


	if(_useBloomFilter){
		//_bloomFilter = new BloomCacheCoherent<u_int64_t>(16000000000ull);
	}

	_filename_smallContigs = _outputDir + "/small_contigs.bin";
	_fileSmallContigs = ofstream(_filename_smallContigs, std::ios_base::app);
	_nbSmallContigs = 0;
	/*
	double nbFP = 0;

	#pragma omp parallel for
	for(size_t i=0; i<100000000ull; i++){
		if(i%2==0) continue;
		_bloomFilter->insert(i);
	}

	cout << "lala" << endl;

	#pragma omp parallel for
	for(size_t i=0; i<100000000ull; i++){
		if(i%2==0){
			if(_bloomFilter->contains(i)){
				nbFP += 1;
			}
		}
		
	}
	cout << "lala" << endl;
	cout << nbFP << endl;
	exit(1);
	*/
	
	//cout << _bloomFilter->contains(5) << endl;
	//_bloomFilter->insert(5);
	//cout << _bloomFilter->contains(5) << endl;
	//while(true){};

	//_kminmerExist.clear();
	_mdbg = new MDBG(_kminmerSize);
	_mdbgNoFilter = new MDBG(_kminmerSize);
	_mdbgSaved = new MDBG(_kminmerSize);



	//string inputFilename = _inputFilename;

	if(!_isFirstPass){
		
		//_mdbgInit = new MDBG(_kminmerSizePrev);
		
		ifstream kminmerFile(_outputDir + "/kminmerData_min.txt");

		while (true) {

			vector<MinimizerType> minimizerSeq;
			minimizerSeq.resize(_kminmerSizePrev);
			kminmerFile.read((char*)&minimizerSeq[0], minimizerSeq.size()*sizeof(MinimizerType));

			if(kminmerFile.eof())break;

			u_int32_t nodeName;
			u_int32_t abundance;
			//u_int32_t quality;
			//bool isReversed = false;

			kminmerFile.read((char*)&nodeName, sizeof(nodeName));
			kminmerFile.read((char*)&abundance, sizeof(abundance));
			//kminmerFile.read((char*)&quality, sizeof(quality));


			if(abundance == 1) continue;
			KmerVec vec;
			vec._kmers = minimizerSeq;

			_kminmerAbundances[vec] = {abundance, 0};
			//_kmerVec_to_nodeName[vec] = nodeName;
		}

		kminmerFile.close();
		


	}

	/*
	if(_kminmerSize >= 6){
		_mdbgFirst = new MDBG(_kminmerSizeFirst);
		_mdbgFirst->load(_outputDir + "/kminmerData_min_init.txt", false);
		checkDuplication();
		delete _mdbgFirst;
		//loadBannedKminmers();
	}
	*/

	string inputFilename_min;
	bool usePos = false;
	
	//_kminmerFile = ofstream(_outputDir + "/kminmerData_min.txt");

	_logFile << "Building kminmer graph" << endl;

	_parsingContigs = false;
	
	auto start = high_resolution_clock::now();

	if(_isFirstPass){
		
		usePos = false;
		inputFilename_min = _outputDir + "/read_data_corrected.txt";

		_logFile << "Filling bloom filter" << endl;
		_nbSolidKminmers = 0;
		//KminmerParserParallel parser(inputFilename_min, _minimizerSize, _kminmerSize, usePos, _nbCores);
		//parser.parse(FillBloomFilter(*this));
		_logFile << "Solid kminmers: " << _nbSolidKminmers << endl;
		//getchar();
		
		_logFile << "Building mdbg" << endl;
		KminmerParserParallel parser2(inputFilename_min, _minimizerSize, _kminmerSize, usePos, true, _nbCores);
		parser2.parse(IndexKminmerFunctor(*this, false));
		_logFile << "Indexed: " << _mdbg->_dbg_nodes.size() << endl;;
		//auto fp = std::bind(&CreateMdbg::createMDBG_collectKminmers_minspace_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
		//parser.parseMinspace(fp);

		if(_useBloomFilter){
			cout << "Attention: bloom filter exact" << endl;
			_bloomFilterExact.clear();
			//delete _bloomFilter;
				
			//cout << "!! Bloom filter disabled!!!!!" << endl;
			KminmerParserParallel parser2(inputFilename_min, _minimizerSize, _kminmerSize, usePos, true, _nbCores);
			parser2.parse(FilterKminmerFunctor(*this));
			_logFile << "Indexed: " << _mdbgSaved->_dbg_nodes.size() << endl;;


			for(auto& it : _mdbgSaved->_dbg_nodes){
				KmerVec vec = it.first;

				u_int32_t nodeName = it.second._index;
				u_int32_t abundance = it.second._abundance;
				//u_int32_t quality = it.second._quality;

				//cout << nodeName << " " << abundance << endl;

				_mdbg->_dbg_nodes[vec] = {nodeName, abundance};
			}

			/*
    		auto it = _mdbg->_dbg_nodes.begin();
			while (it != _mdbg->_dbg_nodes.end()) {
				// Check if key's first character is F
				if (it->second._abundance < 500) {
					// erase() function returns the iterator of the next
					// to last deleted element.
					it = _mdbg->_dbg_nodes.erase(it);
				} else
					it++;
			}
			*/

			delete _mdbgSaved;
		}


	}

	if(!_isFirstPass){
		//const string& filename_uncorrectedReads = _outputDir + "/read_uncorrected.txt";
		const string& filename_contigs = _outputDir + "/unitig_data.txt";
		const string& filename_uncorrectedReads = _outputDir + "/read_data_corrected.txt";

		//KminmerParserParallel parser4(filename_contigs, _minimizerSize, _kminmerSizePrev, false, false, 1);
		//parser4.parse(IndexContigFunctor(*this));

		//_logFile << "Filling bloom filter" << endl;
		//KminmerParserParallel parser(filename_uncorrectedReads, _minimizerSize, _kminmerSize, false, _nbCores);
		//parser.parse(FillBloomFilter(*this));

		_logFile << "Building mdbg" << endl;
		KminmerParserParallel parser2(filename_uncorrectedReads, _minimizerSize, _kminmerSize, false, true, _nbCores);
		//parser2.parse(FilterKminmerFunctor2(*this));
		parser2.parse(IndexKminmerFunctor(*this, false));


		_parsingContigs = true;
		_savingSmallContigs = false;
		
		KminmerParserParallel parser3(filename_contigs, _minimizerSize, _kminmerSize, false, false, _nbCores);
		parser3.parse(IndexKminmerFunctor(*this, true));
		
		_savingSmallContigs = true;
		KminmerParserParallel parser4(filename_contigs, _minimizerSize, _kminmerSize, false, false, _nbCores);
		parser4.parse(IndexKminmerFunctor(*this, true));
	}

	auto stop = high_resolution_clock::now();
	_logFile << "Duration: " << duration_cast<seconds>(stop - start).count() << endl;
	//cout << _mdbg->_dbg_nodes.size() << endl;
	
	
	if(!_isFirstPass){
		_kminmerAbundances.clear();  
		//delete _mdbgInit;
		//_readFile.close();
		
		const auto copyOptions = fs::copy_options::overwrite_existing;
		//fs::copy(_outputDir + "/read_uncorrected.txt", _outputDir + "/read_data.txt", copyOptions);
	}


	writtenNodeNames.clear();

	if(_isFirstPass){
		//if(fs::exists(_outputDir + "/read_data_init.txt")){
		//	fs::remove(_outputDir + "/read_data_init.txt");
		//}
		const auto copyOptions = fs::copy_options::overwrite_existing;
		//fs::copy(_outputDir + "/read_data_init.txt", _outputDir + "/read_data.txt", copyOptions);
		//fs::copy(_outputDir + "/read_data_init.txt", _outputDir + "/read_uncorrected.txt", copyOptions); //disable if correction is enabled
	}


	_logFile << "Nb kminmers (reads): " << _kminmerExist.size() << endl;
	_logFile << "Nb solid kminmers: " << _mdbg->_dbg_nodes.size() << endl;
	
	_kminmerExist.clear();
	_minimizerCounts.clear();
	_kminmersData.clear();
	_contigIndex.clear();
	_contigAbundances.clear();
	_fileSmallContigs.close();

	_logFile << "Nb small contigs written: " << _nbSmallContigs << endl;
	//if(_kminmerSize == _kminmerSizeLast){
	//	removeDuplicatedSmallContigs();
	//}
}


//void CreateMdbg::createMDBG_collectKminmers_minspace_read(const vector<u_int64_t>& readMinimizers, const vector<ReadKminmerComplete>& kminmersInfos, u_int64_t readIndex){



//}


struct KmerVecSorterData{
	vector<MinimizerType> _minimizers;
	u_int32_t _abundance;
	//u_int32_t _quality;
};

static bool KmerVecComparator(const KmerVecSorterData &a, const KmerVecSorterData &b){

	for(size_t i =0; i<a._minimizers.size(); i++){
		if(a._minimizers[i] == b._minimizers[i]){
			continue;
		}
		else{
			return a._minimizers[i] > b._minimizers[i];
		}
	}
	
	return a._abundance > b._abundance;
	//return a._kmerVec._kmers[0] > b._kmerVec._kmers[0];
	//if(a._startNodeIndex == b._startNodeIndex){
	//    return a._length < b._length;
	//}
	//return a._startNodeIndex < b._startNodeIndex;
	//if(a._length == b._length){
	//    return a._startNodeIndex < b._startNodeIndex;
	//}
	//return a._length < b._length;
}


void CreateMdbg::createGfa(){

	u_int32_t nbNodes = _mdbg->_dbg_nodes.size();
	_logFile << "Dumping mdbg nodes" << endl;
	_mdbg->dump(_outputDir + "/kminmerData_min.txt");



	//cout << "A RMETTRE" << endl;
	delete _mdbg;




	//_logFile << "Cleaning repeats..." << endl;

	/*
	vector<KmerVecSorterData> kmerVecs;


	for(const auto& it : _mdbg->_dbg_nodes){
		KmerVec vec = it.first;
		u_int32_t nodeName = it.second._index;

		kmerVecs.push_back({nodeName, vec});
	}

	std::sort(kmerVecs.begin(), kmerVecs.end(), KmerVecComparator);


	unordered_map<u_int32_t, u_int32_t> _nodeName_to_deterministicNodeName;
	
	
	for(size_t i=0; i<kmerVecs.size(); i++){
		_nodeName_to_deterministicNodeName[kmerVecs[i]._nodeName] = i;
	}
	
	for(auto& it : _mdbg->_dbg_nodes){
		it.second._index = _nodeName_to_deterministicNodeName[it.second._index];
	}
	*/

	/*
	for(auto& it : _mdbg->_dbg_nodes){
		KmerVec vec = it.first;

		it.second._index = _nodeName_to_deterministicNodeName[it.second._index];
		
		//KmerVecSorterData& d = kmerVecs[i];
		u_int32_t nodeName = it.second._index;
		u_int32_t abundance = it.second._abundance;

		//bool isReversed;
		//d._kmerVec.normalize(isReversed);

		//vector<u_int64_t> minimizerSeq = d._kmerVec.normalize()._kmers;

		vector<u_int64_t> minimizerSeq = vec._kmers;
		//if(kminmerInfo._isReversed){
		//	std::reverse(minimizerSeq.begin(), minimizerSeq.end());
		//}

		u_int16_t size = minimizerSeq.size();
		//_kminmerFile.write((const char*)&size, sizeof(size));
		_kminmerFile.write((const char*)&minimizerSeq[0], size*sizeof(uint64_t));

		_kminmerFile.write((const char*)&nodeName, sizeof(nodeName));
		_kminmerFile.write((const char*)&abundance, sizeof(abundance));


	}
	*/
	
	
	//"reprise: générer des contigs sans correction"

	//cout << kmerVecs.size() << endl;
	//cout << kmerVecs[0]._kmerVec._kmers.size() << endl;
	//getchar();
	/*
	for(size_t i=0; i<kmerVecs.size(); i++){

		KmerVecSorterData& d = kmerVecs[i];
		u_int32_t nodeName = i;

		vector<u_int64_t> minimizerSeq = d._kmerVec.normalize()._kmers;


		u_int16_t size = minimizerSeq.size();
		_kminmerFile.write((const char*)&size, sizeof(size));
		_kminmerFile.write((const char*)&minimizerSeq, size*sizeof(uint64_t));

		//u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
		u_int32_t length = d._kmerVec._kmers.size();
		u_int32_t lengthStart = 0;
		u_int32_t lengthEnd = length;
		//bool isReversed = kminmerInfo._isReversed;

		_kminmerFile.write((const char*)&nodeName, sizeof(nodeName));
		_kminmerFile.write((const char*)&length, sizeof(length));
		_kminmerFile.write((const char*)&lengthStart, sizeof(lengthStart));
		_kminmerFile.write((const char*)&lengthEnd, sizeof(lengthEnd));
	}
	*/

	//_kminmerFile.close();

	//delete mdbg_repeatFree;
	//cout << _dbg_edges.size() << endl;

	_nbEdges = 0;

	_outputFileGfa.open(_outputDir + "/minimizer_graph.gfa");
	ofstream output_file_gfa_debug(_outputDir + "/minimizer_graph_debug.gfa");

	
	_outputFileGfa.write((const char*)&nbNodes, sizeof(nbNodes));
	//auto start = high_resolution_clock::now();

	/*
	for(const auto& vec_id : _mdbg->_dbg_nodes){

		KmerVec vec = vec_id.first;
		KmerVec vec_rev = vec_id.first.reverse();
		u_int32_t id = vec_id.second._index;
		//u_int32_t id = vec_id.second._index;

		//output_file_gfa << "S" << "\t" << id << "\t" << "*" << "\t" << "LN:i:" << _kminmerLengthMean << "\t" << "dp:i:" << vec_id.second._abundance << endl;
		
		//cout << mdbg->_dbg_edges[vec.prefix().normalize()].size() << endl;
		for(KmerVec& v : _mdbg->_dbg_edges[vec.prefix().normalize()]){
			if(v==vec) continue;
			KmerVec v_rev = v.reverse();
			u_int32_t id2 = _mdbg->_dbg_nodes[v]._index;

			if (vec.suffix() == v.prefix()) {
				nbEdges += 1;
				//u_int16_t overlapLength =  _mdbg->_dbg_nodes[v]._length - _mdbg->_dbg_nodes[v]._overlapLength_end; //   min(_mdbg->_dbg_nodes[v]._overlapLength_end, _mdbg->_dbg_nodes[vec]._overlapLength_end);
				//output_file_gfa << "L" << "\t" << id << "\t" << "+" << "\t" << id2 << "\t" << "+" << "\t" << _kminmerOverlapMean << "M" << endl;
				//vec_add_edge("+", "+");
			}
			if (vec.suffix() == v_rev.prefix()) {
				nbEdges += 1;
				//u_int16_t overlapLength = _mdbg->_dbg_nodes[v]._length - _mdbg->_dbg_nodes[v]._overlapLength_start; //min(_mdbg->_dbg_nodes[v]._overlapLength_end, _mdbg->_dbg_nodes[vec]._overlapLength_end);
				//output_file_gfa << "L" << "\t" << id << "\t" << "+" << "\t" << id2 << "\t" << "-" << "\t" << _kminmerOverlapMean << "M" << endl;
				//vec_add_edge("+", "-");
			}
			if (vec_rev.suffix() == v.prefix()) {
				nbEdges += 1;
				//u_int16_t overlapLength = _mdbg->_dbg_nodes[v]._length - _mdbg->_dbg_nodes[v]._overlapLength_end; //min(_mdbg->_dbg_nodes[v]._overlapLength_start, _mdbg->_dbg_nodes[vec]._overlapLength_start);
				//output_file_gfa << "L" << "\t" << id << "\t" << "-" << "\t" << id2 << "\t" << "+" << "\t" << _kminmerOverlapMean << "M" << endl;
				//vec_add_edge("-", "+");
			}
			if (vec_rev.suffix() == v_rev.prefix()) {
				nbEdges += 1;
				//u_int16_t overlapLength = _mdbg->_dbg_nodes[v]._length - _mdbg->_dbg_nodes[v]._overlapLength_start;; //min(_mdbg->_dbg_nodes[v]._overlapLength_start, _mdbg->_dbg_nodes[vec]._overlapLength_start);
				//output_file_gfa << "L" << "\t" << id << "\t" << "-" << "\t" << id2 << "\t" << "-" << "\t" << _kminmerOverlapMean << "M" << endl;
				//vec_add_edge("-", "-");
			}

		}
		for(KmerVec& v : _mdbg->_dbg_edges[vec.suffix().normalize()]){
			if(v==vec) continue;
			KmerVec v_rev = v.reverse();
			u_int32_t id2 = _mdbg->_dbg_nodes[v]._index;

			if (vec.suffix() == v.prefix()) {
				nbEdges += 1;
				//u_int16_t overlapLength = _mdbg->_dbg_nodes[v]._length - _mdbg->_dbg_nodes[v]._overlapLength_end; //min(_mdbg->_dbg_nodes[v]._overlapLength_end, _mdbg->_dbg_nodes[vec]._overlapLength_end);
				//cout << overlapLength << " " << _mdbg->_dbg_nodes[v]._length<< endl;
				//output_file_gfa << "L" << "\t" << id << "\t" << "+" << "\t" << id2 << "\t" << "+" << "\t" << _kminmerOverlapMean << "M" << endl;
				//vec_add_edge("+", "+");
			}
			if (vec.suffix() == v_rev.prefix()) {
				nbEdges += 1;
				//u_int16_t overlapLength = _mdbg->_dbg_nodes[v]._length - _mdbg->_dbg_nodes[v]._overlapLength_start; //min(_mdbg->_dbg_nodes[v]._overlapLength_end, _mdbg->_dbg_nodes[vec]._overlapLength_end);
				//cout << overlapLength << " " << _mdbg->_dbg_nodes[v]._length<< endl;
				//output_file_gfa << "L" << "\t" << id << "\t" << "+" << "\t" << id2 << "\t" << "-" << "\t" << _kminmerOverlapMean << "M" << endl;
				//vec_add_edge("+", "-");
			}
			if (vec_rev.suffix() == v.prefix()) {
				nbEdges += 1;
				//u_int16_t overlapLength = _mdbg->_dbg_nodes[v]._length - _mdbg->_dbg_nodes[v]._overlapLength_end; //min(_mdbg->_dbg_nodes[v]._overlapLength_start, _mdbg->_dbg_nodes[vec]._overlapLength_start);
				//cout << overlapLength << " " << _mdbg->_dbg_nodes[v]._length << endl;
				//output_file_gfa << "L" << "\t" << id << "\t" << "-" << "\t" << id2 << "\t" << "+" << "\t" << _kminmerOverlapMean << "M" << endl;
				//vec_add_edge("-", "+");
			}
			if (vec_rev.suffix() == v_rev.prefix()) {
				nbEdges += 1;
				//u_int16_t overlapLength = _mdbg->_dbg_nodes[v]._length - _mdbg->_dbg_nodes[v]._overlapLength_start; //min(_mdbg->_dbg_nodes[v]._overlapLength_start, _mdbg->_dbg_nodes[vec]._overlapLength_start);
				//cout << overlapLength << " " << _mdbg->_dbg_nodes[v]._length<< endl;
				//output_file_gfa << "L" << "\t" << id << "\t" << "-" << "\t" << id2 << "\t" << "-" << "\t" << _kminmerOverlapMean << "M" << endl;
				//vec_add_edge("-", "-");
			}
		}
	}
	*/



	
	//vector<KmerVec> keys;
	//for(const auto& it : _mdbg->_dbg_nodes){
	//	keys.push_back(it.first);
	//}
	_logFile << "Computing deterministic node index" << endl;
	computeDeterministicNodeNames();


	_logFile << "Indexing edges" << endl;

	indexEdges();
	
	_logFile << "Computing edges" << endl;

	computeEdges();



	//#pragma omp parallel num_threads(_nbCores)
	//{

		//#pragma omp for
		//for(size_t i=0; i<keys.size(); i++){

			//foo(element->first, element->second);
		//}
		
		//for(KmerVec& v : _dbg_edges_prefix[vec.suffix().normalize()]){
		//	output_file_gfa << "L" << "\t" << vec_id.second << "\t" << "+" << "\t" << _dbg_nodes[v] << "\t" << "+" << "\t" << "0M" << endl;
		//}
	//}

	//auto stop = high_resolution_clock::now();
	//cout << "Duration: " << duration_cast<seconds>(stop - start).count() << endl;

	_outputFileGfa.close();
	output_file_gfa_debug.close();
	
	_logFile << "Nb nodes: " << nbNodes  << endl;
	_logFile << "Nb edges: " << _nbEdges << endl;


	
	//_mdbgNoFilter->dump(_outputDir + "/mdbg_nodes_noFilter.gz");
	//_mdbg->dump(_outputDir + "/mdbg_nodes.gz");
	//_mdbg->_dbg_nodes.clear();
	//_mdbg->_dbg_edges.clear();



}

void CreateMdbg::computeDeterministicNodeNames(){

	ifstream kminmerFile(_outputDir + "/kminmerData_min.txt");

	bool isEOF = false;
	KmerVec vec;
	u_int32_t nodeName;

	vector<KmerVecSorterData> kmerVecs;

	while (true) {


		vector<MinimizerType> minimizerSeq;
		minimizerSeq.resize(_kminmerSize);
		kminmerFile.read((char*)&minimizerSeq[0], minimizerSeq.size()*sizeof(MinimizerType));

		isEOF = kminmerFile.eof();
		if(isEOF) break;

		
		u_int32_t abundance;
		//u_int32_t quality;

		kminmerFile.read((char*)&nodeName, sizeof(nodeName));
		kminmerFile.read((char*)&abundance, sizeof(abundance));
		//kminmerFile.read((char*)&quality, sizeof(quality));
		//vec._kmers = minimizerSeq;


		kmerVecs.push_back({minimizerSeq, abundance});

	}

	kminmerFile.close();

	std::sort(kmerVecs.begin(), kmerVecs.end(), KmerVecComparator);
	nodeName = 0;


	ofstream kminmerFileOut = ofstream(_outputDir + "/kminmerData_min.txt");

	for(auto& entry : kmerVecs){
		//KmerVec vec = it.first;

		//it.second._index = _nodeName_to_deterministicNodeName[it.second._index];
		
		//KmerVecSorterData& d = kmerVecs[i];
		//u_int32_t nodeName = it.second._index;
		u_int32_t abundance = entry._abundance;
		//u_int32_t quality = entry._quality;

		//if(quality==0){
		//	cout << "omg " << nodeName << endl;
		//	getchar();
		//}
		//bool isReversed;
		//d._kmerVec.normalize(isReversed);

		//vector<u_int64_t> minimizerSeq = d._kmerVec.normalize()._kmers;

		vector<MinimizerType> minimizerSeq = entry._minimizers;
		//if(kminmerInfo._isReversed){
		//	std::reverse(minimizerSeq.begin(), minimizerSeq.end());
		//}

		u_int16_t size = minimizerSeq.size();
		//_kminmerFile.write((const char*)&size, sizeof(size));
		kminmerFileOut.write((const char*)&minimizerSeq[0], size*sizeof(MinimizerType));


		kminmerFileOut.write((const char*)&nodeName, sizeof(nodeName));
		kminmerFileOut.write((const char*)&abundance, sizeof(abundance));
		//kminmerFileOut.write((const char*)&quality, sizeof(quality));

		//cout << nodeName << " " <<  abundance << endl;
		
		//if(minimizerSeq[0] == 7892019983131394 && minimizerSeq[1] == 62446187524336769 && minimizerSeq[2] == 47730431539434127 && minimizerSeq[3] == 90993437045220547){
		//	cout << nodeName << endl;
		//	getchar();
		//}
		//cout << nodeName << endl;
		//for(u_int64_t m : minimizerSeq) cout << m << " ";
		//cout << abundance << " " << quality << endl;
		//cout << endl;
		//getchar();
		
		nodeName += 1;
	}


	kminmerFileOut.close();


}

void CreateMdbg::indexEdges(){

	ifstream kminmerFile(_outputDir + "/kminmerData_min.txt");

	#pragma omp parallel num_threads(_nbCores)
	{

		bool isEOF = false;
		KmerVec vec;
		u_int32_t nodeName;

		while (true) {

			#pragma omp critical(indexEdge)
			{

				vector<MinimizerType> minimizerSeq;
				minimizerSeq.resize(_kminmerSize);
				kminmerFile.read((char*)&minimizerSeq[0], minimizerSeq.size()*sizeof(MinimizerType));

				isEOF = kminmerFile.eof();

				
				u_int32_t abundance;
				//u_int32_t quality;

				if(!isEOF){
					kminmerFile.read((char*)&nodeName, sizeof(nodeName));
					kminmerFile.read((char*)&abundance, sizeof(abundance));
					//kminmerFile.read((char*)&quality, sizeof(quality));
					vec._kmers = minimizerSeq;
				}

				//cout << nodeName << endl;
			}

			if(isEOF) break;
			
			indexEdge(vec, nodeName);

		}
	}

	kminmerFile.close();
}

void CreateMdbg::indexEdge(const KmerVec& vec, u_int32_t nodeName){

	
	bool needprint = false;

	bool isReversedPrefix;
	bool isReversedSuffix;
	KmerVec prefix = vec.prefix().normalize(isReversedPrefix);
	KmerVec suffix = vec.suffix().normalize(isReversedSuffix);
	MinimizerType prefix_minimizer = vec._kmers[vec._kmers.size()-1];
	MinimizerType suffix_minimizer = vec._kmers[0];
	
	/*
	//5759311437506739 24292623053282490 90668338080075164 46010529315713468
	if(std::find(vec._kmers.begin(), vec._kmers.end(), 5759311437506739) != vec._kmers.end() && std::find(vec._kmers.begin(), vec._kmers.end(), 24292623053282490) != vec._kmers.end() && std::find(vec._kmers.begin(), vec._kmers.end(), 90668338080075164) != vec._kmers.end()){
		cout << "--------------" << endl;
		cout << isReversedPrefix << endl;
		for(u_int64_t m : vec._kmers){
			cout << m << " ";
		}
		cout << endl;

		cout << prefix_minimizer << endl;
		needprint = true;

		getchar();
	}
	*/
	

	  

	/*
	cout << "--------------" << endl;
	cout << isReversedPrefix << endl;
	for(u_int64_t m : vec._kmers){
		cout << m << " ";
	}
	cout << endl;

	for(u_int64_t m : prefix._kmers){
		cout << m << " ";
	}
	cout << endl;

	cout << prefix_minimizer << endl;

	KmerVec vecNew = prefix;
	if(isReversedPrefix){
		vecNew = vecNew.reverse();
	}
	vecNew._kmers.push_back(prefix_minimizer);
	for(u_int64_t m : vecNew._kmers){
		cout << m << " ";
	}
	cout << endl;
	getchar();
	*/
	/*
	cout << "--------------" << endl;
	cout << isReversedSuffix << endl;
	for(u_int64_t m : vec._kmers){
		cout << m << " ";
	}
	cout << endl;

	for(u_int64_t m : suffix._kmers){
		cout << m << " ";
	}
	cout << endl;

	cout << suffix_minimizer << endl;

	KmerVec vecNew = suffix;
	if(isReversedSuffix){
		vecNew = vecNew.reverse();
	}
	vecNew._kmers.insert(vecNew._kmers.begin(), suffix_minimizer);
	for(u_int64_t m : vecNew._kmers){
		cout << m << " ";
	}
	cout << endl;
	getchar();
	*/

	/*
	_mdbgEdgesNew.lazy_emplace_l(prefix, 
	[&prefix_minimizer, &nodeName, &isReversedPrefix](MdbgEdgeMap2::value_type& v) { // key exist
		v.second.push_back({nodeName, prefix_minimizer, isReversedPrefix, true});
	},           
	[&prefix, &prefix_minimizer, &nodeName, &isReversedPrefix](const MdbgEdgeMap2::constructor& ctor) { // key inserted
		vector<KminmerEdge2> nodes;
		nodes.push_back({nodeName, prefix_minimizer, isReversedPrefix, true});
		ctor(prefix, nodes); 
	});

	_mdbgEdgesNew.lazy_emplace_l(suffix, 
	[&suffix_minimizer, &nodeName, &isReversedSuffix](MdbgEdgeMap2::value_type& v) { // key exist
		v.second.push_back({nodeName, suffix_minimizer, isReversedSuffix, false});
	},           
	[&suffix, &suffix_minimizer, &nodeName, &isReversedSuffix](const MdbgEdgeMap2::constructor& ctor) { // key inserted
		vector<KminmerEdge2> nodes;
		nodes.push_back({nodeName, suffix_minimizer, isReversedSuffix, false});
		ctor(suffix, nodes); 
	});
	*/

	/*
	KmerVec edge1 = vec.prefix().normalize();
	KmerVec edge2 = vec.suffix().normalize();

	_mdbgEdges.lazy_emplace_l(edge1, 
	[&vec, &nodeName](MdbgEdgeMap::value_type& v) { // key exist
		v.second.push_back({nodeName, vec});
	},           
	[&edge1, &vec, &nodeName](const MdbgEdgeMap::constructor& ctor) { // key inserted
		vector<KminmerEdge> nodes;
		nodes.push_back({nodeName, vec});
		ctor(edge1, nodes); 
	});

	_mdbgEdges.lazy_emplace_l(edge2, 
	[&vec, &nodeName](MdbgEdgeMap::value_type& v) { // key exist
		v.second.push_back({nodeName, vec});
	},           
	[&edge2, &vec, &nodeName](const MdbgEdgeMap::constructor& ctor) { // key inserted
		vector<KminmerEdge> nodes;
		nodes.push_back({nodeName, vec});
		ctor(edge2, nodes); 
	});
	*/

	
	_mdbgEdges2.lazy_emplace_l(prefix, 
	[&prefix_minimizer, &nodeName, &isReversedPrefix](MdbgEdgeMap2::value_type& v) { // key exist
		v.second.push_back({nodeName, prefix_minimizer, isReversedPrefix, true});
	},           
	[&prefix, &prefix_minimizer, &nodeName, &isReversedPrefix](const MdbgEdgeMap2::constructor& ctor) { // key inserted
		vector<KminmerEdge2> nodes;
		nodes.push_back({nodeName, prefix_minimizer, isReversedPrefix, true});
		ctor(prefix, nodes); 
	});

	_mdbgEdges2.lazy_emplace_l(suffix, 
	[&suffix_minimizer, &nodeName, &isReversedSuffix](MdbgEdgeMap2::value_type& v) { // key exist
		v.second.push_back({nodeName, suffix_minimizer, isReversedSuffix, false});
	},           
	[&suffix, &suffix_minimizer, &nodeName, &isReversedSuffix](const MdbgEdgeMap2::constructor& ctor) { // key inserted
		vector<KminmerEdge2> nodes;
		nodes.push_back({nodeName, suffix_minimizer, isReversedSuffix, false});
		ctor(suffix, nodes); 
	});
	
	//cout << _mdbgEdges.size() << endl;
}

void CreateMdbg::computeEdges(){
	
	ifstream kminmerFile(_outputDir + "/kminmerData_min.txt");

	#pragma omp parallel num_threads(_nbCores)
	{

		bool isEOF = false;
		KmerVec vec;
		u_int32_t nodeName;

		while (true) {

			#pragma omp critical(indexEdge)
			{

				vector<MinimizerType> minimizerSeq;
				minimizerSeq.resize(_kminmerSize);
				kminmerFile.read((char*)&minimizerSeq[0], minimizerSeq.size()*sizeof(MinimizerType));

				isEOF = kminmerFile.eof();

				
				u_int32_t abundance;
				//u_int32_t quality;

				if(!isEOF){
					kminmerFile.read((char*)&nodeName, sizeof(nodeName));
					kminmerFile.read((char*)&abundance, sizeof(abundance));
					//kminmerFile.read((char*)&quality, sizeof(quality));
					vec._kmers = minimizerSeq;
				}

			}

			if(isEOF) break;
			
			computeEdge(vec, nodeName);

		}
	}

	kminmerFile.close();
}

void CreateMdbg::computeEdge(const KmerVec& vec, u_int32_t id){


	static u_int8_t isS = 0;
	static u_int8_t isL = 1;
	static u_int8_t edgePlus = 0;
	static u_int8_t edgeMinus = 1;

	//KmerVec vec = keys[i];
	KmerVec vec_rev = vec.reverse();

	const KmerVec prefix = vec.prefix().normalize();


	for(KminmerEdge2& edge : _mdbgEdges2[prefix]){

			
		//KmerVec v = edge._vec;
		
		KmerVec v = prefix;
		if(edge._isReversed){
			v = v.reverse();
			//std::reverse(v._kmers.begin(), v._kmers.end());
		}
		//v._kmers.push_back(edge._minimizer);
		
		
		if(edge._isPrefix){
			v._kmers.push_back(edge._minimizer);
		}
		else{
			v._kmers.insert(v._kmers.begin(), edge._minimizer);
		}
		
		//v = v.normalize();

		/*
		if(_mdbg->_dbg_nodes.find(v) == _mdbg->_dbg_nodes.end()){

			for(u_int64_t m : v._kmers){
				cout << m << " ";
			}
			cout << endl;
			for(u_int64_t m : prefix._kmers){
				cout << m << " ";
			}
			cout << endl;
			
			cout << "prefix" << endl;
			getchar();
		}
		*/
		

		if(v==vec) continue;

		KmerVec v_rev = v.reverse();
		u_int32_t id2 = edge._nodeName;

		if (vec.suffix() == v.prefix()) {
			dumpEdge(id, edgePlus, id2, edgePlus);
		}
		if (vec.suffix() == v_rev.prefix()) {
			dumpEdge(id, edgePlus, id2, edgeMinus);
		}
		if (vec_rev.suffix() == v.prefix()) {
			dumpEdge(id, edgeMinus, id2, edgePlus);
		}
		if (vec_rev.suffix() == v_rev.prefix()) {
			dumpEdge(id, edgeMinus, id2, edgeMinus);
		}

	}

	
	KmerVec suffix = vec.suffix().normalize();
	for(KminmerEdge2& edge : _mdbgEdges2[suffix]){
		
		//KmerVec v = edge._vec;

		
		KmerVec v = suffix;
		if(edge._isReversed){
			v = v.reverse();
		}
		if(edge._isPrefix){
			v._kmers.push_back(edge._minimizer);
		}
		else{
			v._kmers.insert(v._kmers.begin(), edge._minimizer);
		}
		//v = v.normalize();
		//v._kmers.insert(v._kmers.begin(), edge._minimizer);

		/*
		if(_mdbg->_dbg_nodes.find(v) == _mdbg->_dbg_nodes.end()){
			cout << "suffix" << endl;
			getchar();
		}
		*/
		

		if(v==vec) continue;

		KmerVec v_rev = v.reverse();
		u_int32_t id2 = edge._nodeName;

		if (vec.suffix() == v.prefix()) {
			dumpEdge(id, edgePlus, id2, edgePlus);
		}
		if (vec.suffix() == v_rev.prefix()) {
			dumpEdge(id, edgePlus, id2, edgeMinus);
		}
		if (vec_rev.suffix() == v.prefix()) {
			dumpEdge(id, edgeMinus, id2, edgePlus);
		}
		if (vec_rev.suffix() == v_rev.prefix()) {
			dumpEdge(id, edgeMinus, id2, edgeMinus);
		}
	}

}

void CreateMdbg::dumpEdge(u_int32_t nodeNameFrom, u_int8_t nodeNameFromOri, u_int32_t nodeNameTo, u_int8_t nodeNameToOri){
	
	#pragma omp critical(gfa)
	{
		
		_nbEdges += 1;
		_outputFileGfa.write((const char*)&nodeNameFrom, sizeof(nodeNameFrom));
		_outputFileGfa.write((const char*)&nodeNameFromOri, sizeof(nodeNameFromOri));
		_outputFileGfa.write((const char*)&nodeNameTo, sizeof(nodeNameTo));
		_outputFileGfa.write((const char*)&nodeNameToOri, sizeof(nodeNameToOri));
	}
}

