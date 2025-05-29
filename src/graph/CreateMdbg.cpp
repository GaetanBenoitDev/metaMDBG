

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
	args::ValueFlag<int> arg_minAbundance(parser, "", "Minimum k-mer-mer abundance", {ARG_MIN_KMINMER_ABUNDANCE}, 0);
	args::ValueFlag<int> arg_nbCores(parser, "", "Number of cores", {ARG_NB_CORES2}, NB_CORES_DEFAULT_INT);
	args::Flag arg_useCorrectedRead(parser, "", "Use corrected reads", {"corrected-read"});
	//args::Flag arg_bf(parser, "", "Filter unique erroneous k-minmers", {ARG_BLOOM_FILTER});
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
	_minAbundance = args::get(arg_minAbundance);

	_useBloomFilter = true;
	//if(arg_bf){
	//	_useBloomFilter = false;
	//}

	_isFirstPass = false;
	if(arg_firstPass){
		_isFirstPass = true;
	}

	if(!_isFirstPass){
		_useBloomFilter = false;	
	}

	_useCorrectedRead = false;
	if(arg_useCorrectedRead){
		_useCorrectedRead = true;
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


	Logger::get().debug() << "";
	Logger::get().debug() << "Contig filename: " << _filename_inputContigs;
	Logger::get().debug() << "Output dir: " << _outputDir;
	Logger::get().debug() << "Minimizer length: " << _minimizerSize;
	Logger::get().debug() << "Kminmer length: " << _kminmerSize;
	Logger::get().debug() << "Density: " << _minimizerDensity;
	Logger::get().debug() << "Min abundance: " << _minAbundance;
	Logger::get().debug() << "";

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


	_lalaLol = 0;

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

	//closeLogFile();

	
}



void CreateMdbg::createMDBG (){

	

	_kminmerFile.open(_outputDir + "/kminmerData_min.txt");

	if(_useBloomFilter){
		_bloomFilter = new BloomCacheCoherent<u_int64_t>(32000000000ull);
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
	//_mdbg = new MDBG(_kminmerSize);
	//_mdbgNoFilter = new MDBG(_kminmerSize);
	//_mdbgSaved = new MDBG(_kminmerSize);



	//string inputFilename = _inputFilename;

	if(!_isFirstPass){
		
		Logger::get().debug() << "Loading refined abundances";
		auto start = high_resolution_clock::now();
		loadRefinedAbundances();
		Logger::get().debug() << "\tDuration: " << duration_cast<seconds>(high_resolution_clock::now() - start).count();
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
	//bool usePos = false;
	
	//_kminmerFile = ofstream(_outputDir + "/kminmerData_min.txt");

	Logger::get().debug() << "Collecting kminmers";
	auto start = high_resolution_clock::now();

	_parsingContigs = false;
	

	if(_isFirstPass){
		
		//usePos = false;
		inputFilename_min = _outputDir + "/read_data_init.txt";
		if(_useCorrectedRead){
			inputFilename_min = _outputDir + "/read_data_corrected.txt";
		}

		Logger::get().debug() << "Filling bloom filter";
		_nbSolidKminmers = 0;
		//KminmerParserParallel parser(inputFilename_min, _minimizerSize, _kminmerSize, usePos, _nbCores);
		//parser.parse(FillBloomFilter(*this));
		Logger::get().debug() << "Solid kminmers: " << _nbSolidKminmers;
		//getchar();
		
		Logger::get().debug() << "Building mdbg";
		KminmerParserParallel parser2(inputFilename_min, _minimizerSize, _kminmerSize, false, false, _nbCores);
		parser2.parse(IndexKminmerFunctor(*this, false));
		Logger::get().debug() << "Indexed: " << _mdbgNodesLight.size();
		//auto fp = std::bind(&CreateMdbg::createMDBG_collectKminmers_minspace_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
		//parser.parseMinspace(fp);

		if(_useBloomFilter){
			//cout << "Attention: bloom filter exact" << endl;
			//_bloomFilterExact.clear();
			delete _bloomFilter;
				
			if(_minAbundance <= 1){
				//cout << "!! Bloom filter disabled!!!!!" << endl;
				KminmerParserParallel parser2(inputFilename_min, _minimizerSize, _kminmerSize, false, false, _nbCores);
				parser2.parse(FilterKminmerFunctor(*this));
				Logger::get().debug() << "Indexed: " << _mdbgNodesLightUnique.size();


				for(auto& it : _mdbgNodesLightUnique){
					u_int128_t hash = it.first;

					u_int32_t nodeName = it.second._index;
					u_int32_t abundance = it.second._abundance;

					//u_int32_t quality = it.second._quality;

					//cout << nodeName << " " << abundance << endl;

					_mdbgNodesLight[hash] = {false, nodeName, abundance};
				}


			
			}

			_mdbgNodesLightUnique.clear();
			//delete _mdbgSaved;

				
			if(_minAbundance > 1){
				
				auto it = _mdbgNodesLight.begin();

				while (it != _mdbgNodesLight.end()) {
					
					if (it->second._abundance < _minAbundance) {
						it = _mdbgNodesLight.erase(it);
					}
					else{
						it++;
					}
				}
			}
				
		}


	}

	if(!_isFirstPass){
		

		//const string& filename_uncorrectedReads = _outputDir + "/read_uncorrected.txt";
		const string& filename_contigs = _outputDir + "/unitig_data.txt";
		string filename_uncorrectedReads = _outputDir + "/read_data_init.txt";
		if(_useCorrectedRead){
			filename_uncorrectedReads = _outputDir + "/read_data_corrected.txt";
		}

		//KminmerParserParallel parser4(filename_contigs, _minimizerSize, _kminmerSizePrev, false, false, 1);
		//parser4.parse(IndexContigFunctor(*this));

		//_logFile << "Filling bloom filter" << endl;
		//KminmerParserParallel parser(filename_uncorrectedReads, _minimizerSize, _kminmerSize, false, _nbCores);
		//parser.parse(FillBloomFilter(*this));

		Logger::get().debug() << "Building mdbg";
		
		KminmerParserParallel parser2(filename_uncorrectedReads, _minimizerSize, _kminmerSize, false, false, _nbCores);
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

	Logger::get().debug() << "Duration: " << duration_cast<seconds>(high_resolution_clock::now() - start).count();
	//cout << _mdbg->_dbg_nodes.size() << endl;
	
	
	//if(!_isFirstPass){
	_kminmerAbundances.clear();  
		//delete _mdbgInit;
		//_readFile.close();
		
		//const auto copyOptions = fs::copy_options::overwrite_existing;
		//fs::copy(_outputDir + "/read_uncorrected.txt", _outputDir + "/read_data.txt", copyOptions);
	//}


	//writtenNodeNames.clear();

	//if(_isFirstPass){
		//if(fs::exists(_outputDir + "/read_data_init.txt")){
		//	fs::remove(_outputDir + "/read_data_init.txt");
		//}
		//const auto copyOptions = fs::copy_options::overwrite_existing;
		//fs::copy(_outputDir + "/read_data_init.txt", _outputDir + "/read_data.txt", copyOptions);
		//fs::copy(_outputDir + "/read_data_init.txt", _outputDir + "/read_uncorrected.txt", copyOptions); //disable if correction is enabled
	//}


	Logger::get().debug() << "Nb solid kminmers: " << _mdbgNodesLight.size();
	
	_fileSmallContigs.close();

	Logger::get().debug() << "Nb small contigs written: " << _nbSmallContigs;


	Logger::get().debug() << "Dump kminmer abundances";
	start = high_resolution_clock::now();
	ofstream abundanceFile(_outputDir + "/kminmerData_abundance.txt");
	for(const auto& it : _mdbgNodesLight){

		u_int128_t vecHash = it.first;

		//u_int32_t nodeName = it.second._index;
		u_int32_t abundance = it.second._abundance;

		MDBG::writeKminmerAbundance(vecHash, abundance, abundanceFile);
	}
	Logger::get().debug() << "\tDuration: " << duration_cast<seconds>(high_resolution_clock::now() - start).count();
	abundanceFile.close();

	if(_isFirstPass){
		fs::copy(_outputDir + "/kminmerData_abundance.txt", _outputDir + "/kminmerData_abundance_init.txt", fs::copy_options::overwrite_existing);
	}


	if(_kminmerSize == _kminmerSizeFirst+1){
		fs::copy(_outputDir + "/kminmerData_abundance.txt", _outputDir + "/kminmerData_abundance_init_k" + to_string(_kminmerSizeFirst+1) + ".txt", fs::copy_options::overwrite_existing);
	}

	//if(_kminmerSize == _kminmerSizeLast){
	//	removeDuplicatedSmallContigs();
	//}

	//u_int32_t nbNodes = _mdbgNodesLight.size();
	//
	//dumpNodes();
	//cout << "MDBG nb nodes: " << _mdbgNodesLight.size() << endl;
	//Logger::get().debug() << "Duration: " << duration_cast<seconds>(high_resolution_clock::now() - start).count();
	_mdbgNodesLight.clear();
	//_mdbg->dump(_outputDir + "/kminmerData_min.txt");

	_kminmerFile.close();

	//cout << "A RMETTRE" << endl;
	//delete _mdbg;

}





struct KmerVecSorterData{
	UnitigType _unitigIndex;
	vector<MinimizerType> _minimizers;
	//u_int32_t _abundance;
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
	
	return false;
};



/*
void CreateMdbg::dumpNodes(){
	

	ofstream kminmerFile = ofstream(_outputDir + "/kminmerData_min.txt");


	string readFilename = _outputDir + "/read_data_init.txt";
	if(_useCorrectedRead){
		readFilename = _outputDir + "/read_data_corrected.txt";
	}

	
	KminmerParserParallel parser(readFilename, _minimizerSize, _kminmerSize, false, false, _nbCores);
	parser.parse(DumpKminmerFunctor(*this, kminmerFile));

	if(!_isFirstPass){
		//const string& filename_uncorrectedReads = _outputDir + "/read_uncorrected.txt";
		const string& contigFilename = _outputDir + "/unitig_data.txt";


		KminmerParserParallel parser2(contigFilename, _minimizerSize, _kminmerSize, false, false, _nbCores);
		parser2.parse(DumpKminmerFunctor(*this, kminmerFile));
		
	}

	kminmerFile.close();

}
*/

void CreateMdbg::createGfa(){






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

	//ofstream output_file_gfa_debug(_outputDir + "/minimizer_graph_debug.gfa");

	
	//_outputFileGfa.write((const char*)&nbNodes, sizeof(nbNodes));
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
	//Logger::get().debug() << "Computing deterministic node index";
	//computeDeterministicNodeNames();
	/*
	_kminmerAbundances.clear();
	ifstream kminmerFile(_outputDir + "/kminmerData_min.txt");

	bool isEOF = false;
	KmerVec vec;
	u_int32_t nodeName;


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
		vec._kmers = minimizerSeq;

		_kminmerAbundances[vec.hash128()] = abundance;

	}

	kminmerFile.close();
	*/
	

	Logger::get().debug() << "Indexing edges";
	auto start = high_resolution_clock::now();
	indexEdges();
	Logger::get().debug() << "\tDuration: " << duration_cast<seconds>(high_resolution_clock::now() - start).count();
	

	//Logger::get().debug() << "Computing edges";
	//computeEdgesOld();

	//cout << "gfa nb edges: " << _nbEdges << endl;
	//getchar();




	Logger::get().debug() << "Compute unitig nodes";
	start = high_resolution_clock::now();
	_unitigGraphFile_nodes = ofstream(_outputDir + "/unitigGraph.nodes.bin");
	computeUnitigNodes();
	Logger::get().debug() << "\tDuration: " << duration_cast<seconds>(high_resolution_clock::now() - start).count();
	_unitigGraphFile_nodes.close();
	_mdbgEdges3.clear();


	//Logger::get().debug() << "Compute deterministic start end nodes";
	//computeUnitiStartEndNode();



	Logger::get().debug() << "Index unitig edges";
	start = high_resolution_clock::now();
	indexUnitigEdges();
	Logger::get().debug() << "\tDuration: " << duration_cast<seconds>(high_resolution_clock::now() - start).count();



	Logger::get().debug() << "Compute unitig edges";
	start = high_resolution_clock::now();

	_unitigGraphFile_edges_successors = ofstream(_outputDir + "/unitigGraph.edges.successors.bin");
	_unitigGraphFile_edges_predecessors = ofstream(_outputDir + "/unitigGraph.edges.predecessors.bin");

	computeUnitigEdges();

	Logger::get().debug() << "\tDuration: " << duration_cast<seconds>(high_resolution_clock::now() - start).count();
	_unitigGraphFile_edges_successors.close();
	_unitigGraphFile_edges_predecessors.close();
	_startNode_to_unitigIndex.clear();
	_endNode_to_unitigIndex.clear();

	_mdbgEdges4.clear();
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

	//output_file_gfa_debug.close();
	
	

	Logger::get().debug() << "Dump unitig abundances";
	start = high_resolution_clock::now();
	_unitigGraphFile_nodes_abundances = ofstream(_outputDir + "/unitigGraph.nodes.abundances.bin");
	dumpUnitigAbundances();
	Logger::get().debug() << "\tDuration: " << duration_cast<seconds>(high_resolution_clock::now() - start).count();
	_unitigGraphFile_nodes_abundances.close();


	Logger::get().debug() << "Nb edges:   " << _nbEdges;


	
	//_mdbgNoFilter->dump(_outputDir + "/mdbg_nodes_noFilter.gz");
	//_mdbg->dump(_outputDir + "/mdbg_nodes.gz");
	//_mdbg->_dbg_nodes.clear();
	//_mdbg->_dbg_edges.clear();



}


void CreateMdbg::computeDeterministicNodeNames(){


	ifstream kminmerFile(_outputDir + "/kminmerData_min.txt");

	bool isEOF = false;
	//KmerVec vec;
	u_int32_t nodeName;

	vector<KmerVecSorterData> kmerVecs;
	vector<KmerVecSorterData> kmerVecs2;
	//vector<vector<MinimizerType>> kmerVecs;

	while (true) {


		vector<MinimizerType> minimizers;
		isEOF = MDBG::readKminmer(minimizers, kminmerFile, _kminmerSize);
		
		isEOF = kminmerFile.eof();
		if(isEOF) break;

		//vec._kmers = minimizers;

		//vector<MinimizerType> minimizerSeq;
		//minimizerSeq.resize(_kminmerSize);
		//kminmerFile.read((char*)&minimizerSeq[0], minimizerSeq.size()*sizeof(MinimizerType));


		
		//u_int32_t abundance;
		//u_int32_t quality;

		//kminmerFile.read((char*)&nodeName, sizeof(nodeName));
		//kminmerFile.read((char*)&abundance, sizeof(abundance));
		//kminmerFile.read((char*)&quality, sizeof(quality));
		//vec._kmers = minimizerSeq;


		//kmerVecs.push_back({minimizerSeq, abundance});
		kmerVecs.push_back({0, minimizers});


		kmerVecs2.push_back({0, minimizers});
		std::reverse(minimizers.begin(), minimizers.end());
		kmerVecs2.push_back({0, minimizers});
	}

	kminmerFile.close();

	std::sort(kmerVecs.begin(), kmerVecs.end(), KmerVecComparator);
	std::sort(kmerVecs2.begin(), kmerVecs2.end(), KmerVecComparator);
	nodeName = 0;


	ofstream kminmerFileOut = ofstream(_outputDir + "/kminmerData_min.txt");

	for(const KmerVecSorterData& data : kmerVecs){
		//KmerVec vec = it.first;

		//it.second._index = _nodeName_to_deterministicNodeName[it.second._index];
		
		//KmerVecSorterData& d = kmerVecs[i];
		//u_int32_t nodeName = it.second._index;
		//u_int32_t abundance = entry._abundance;
		//u_int32_t quality = entry._quality;

		//if(quality==0){
		//	cout << "omg " << nodeName << endl;
		//	getchar();
		//}
		//bool isReversed;
		//d._kmerVec.normalize(isReversed);

		//vector<u_int64_t> minimizerSeq = d._kmerVec.normalize()._kmers;

		//vector<MinimizerType> minimizerSeq = entry._minimizers;
		
		MDBG::writeKminmer(data._minimizers, kminmerFileOut);

		KmerVec vec;
		vec._kmers = data._minimizers;

		_kmervec_to_nodeName[vec] = nodeName;
		
		nodeName += 1;
		//if(kminmerInfo._isReversed){
		//	std::reverse(minimizerSeq.begin(), minimizerSeq.end());
		//}

		//u_int16_t size = minimizerSeq.size();
		//_kminmerFile.write((const char*)&size, sizeof(size));
		//kminmerFileOut.write((const char*)&minimizerSeq[0], size*sizeof(MinimizerType));


		//kminmerFileOut.write((const char*)&nodeName, sizeof(nodeName));
		//kminmerFileOut.write((const char*)&abundance, sizeof(abundance));
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
		
	}


	kminmerFileOut.close();

	//cout << "gfa nb nodes: " << kmerVecs.size() << endl;


}


void CreateMdbg::indexEdges(){

	//cout << "About to index edges" << endl;
	//getchar();
	
	_edgeIndexElements = 0;

	ifstream kminmerFile(_outputDir + "/kminmerData_min.txt");

	#pragma omp parallel num_threads(_nbCores)
	{

		bool isEOF = false;
		KmerVec vec;

		while (true) {

			#pragma omp critical(indexEdge)
			{

				vector<MinimizerType> minimizers;
				isEOF = MDBG::readKminmer(minimizers, kminmerFile, _kminmerSize);
				
				vec._kmers = minimizers;
				//kminmerFile.read((char*)&minimizerSeq[0], minimizerSeq.size()*sizeof(MinimizerType));

				//isEOF = kminmerFile.eof();

				
				//u_int32_t abundance;
				//u_int32_t quality;

				//if(!isEOF){
				//	kminmerFile.read((char*)&nodeName, sizeof(nodeName));
				//	kminmerFile.read((char*)&abundance, sizeof(abundance));
					//kminmerFile.read((char*)&quality, sizeof(quality));
				//}

				//cout << nodeName << endl;
			}

			if(isEOF) break;
			
			indexEdge(vec, 0);
			//indexEdge(vec.reverse(), nodeName);

		}
	}

	kminmerFile.close();

	Logger::get().debug() << "Edge index size: " << _mdbgEdges3.size();
	Logger::get().debug() << "Edge index elements: " << _edgeIndexElements;
	//cout << _mdbgEdges.size() << endl;
	//getchar();

}

void CreateMdbg::indexEdge(const KmerVec& vec, u_int32_t nodeName){

	u_int64_t& lala = _edgeIndexElements;
	bool needprint = false;

	bool isReversedPrefix;
	bool isReversedSuffix;
	KmerVec prefix = vec.prefix().normalize(isReversedPrefix);
	KmerVec suffix = vec.suffix().normalize(isReversedSuffix);
	MinimizerType prefix_minimizer = vec._kmers[vec._kmers.size()-1];
	MinimizerType suffix_minimizer = vec._kmers[0];
	
	/*
	_mdbgEdges2.lazy_emplace_l(prefix.hash128(), 
	[&prefix_minimizer, &nodeName, &isReversedPrefix](MdbgEdgeMap2::value_type& v) { // key exist
		v.second.push_back({nodeName, prefix_minimizer, isReversedPrefix, true});
	},           
	[&prefix, &prefix_minimizer, &nodeName, &isReversedPrefix](const MdbgEdgeMap2::constructor& ctor) { // key inserted
		vector<KminmerEdge2> nodes;
		nodes.push_back({nodeName, prefix_minimizer, isReversedPrefix, true});
		ctor(prefix.hash128(), nodes); 
	});

	_mdbgEdges2.lazy_emplace_l(suffix.hash128(), 
	[&suffix_minimizer, &nodeName, &isReversedSuffix](MdbgEdgeMap2::value_type& v) { // key exist
		v.second.push_back({nodeName, suffix_minimizer, isReversedSuffix, false});
	},           
	[&suffix, &suffix_minimizer, &nodeName, &isReversedSuffix](const MdbgEdgeMap2::constructor& ctor) { // key inserted
		vector<KminmerEdge2> nodes;
		nodes.push_back({nodeName, suffix_minimizer, isReversedSuffix, false});
		ctor(suffix.hash128(), nodes); 
	});
	*/

	_mdbgEdges3.lazy_emplace_l(prefix.hash128(), 
	[&prefix_minimizer, &isReversedPrefix, &lala](MdbgEdgeMap3::value_type& v) { // key exist
		lala += 1;
		v.second.push_back({prefix_minimizer, isReversedPrefix, true, false});
	},           
	[&prefix, &prefix_minimizer, &isReversedPrefix, &lala](const MdbgEdgeMap3::constructor& ctor) { // key inserted
		lala += 1;
		vector<KminmerEdge3> nodes;
		nodes.push_back({prefix_minimizer, isReversedPrefix, true, false});
		ctor(prefix.hash128(), nodes); 
	});

	_mdbgEdges3.lazy_emplace_l(suffix.hash128(), 
	[&suffix_minimizer, &isReversedSuffix, &lala](MdbgEdgeMap3::value_type& v) { // key exist
		lala += 1;
		v.second.push_back({suffix_minimizer, isReversedSuffix, false, false});
	},           
	[&suffix, &suffix_minimizer, &isReversedSuffix, &lala](const MdbgEdgeMap3::constructor& ctor) { // key inserted
		lala += 1;
		vector<KminmerEdge3> nodes;
		nodes.push_back({suffix_minimizer, isReversedSuffix, false, false});
		ctor(suffix.hash128(), nodes); 
	});
	


}

/*

void CreateMdbg::computeUnitigNodes(){
	

	cout << "computeUnitigNodes: besoin d'activer le if(_graph->_isNodeIndexIndexed[nodeIndex]) exist = true (requiert d'indexer tous les kmervec des unitigs)" << endl;
	cout << "Parallelisation unitig creation: bien faire en sorte que les nodeSource soit randomiser, surtout ne pas paralleliser a partir des node d'un read successif qui vont probablement etre dans le meme unitig" << endl;
	cout << "nb cores 1" << endl;
	cout << "Opti: getSuccessor, getPredecessor: on peut arreter la collect des que ya plus de 1 successor/predecessor pendant l'unitigage" << endl;
	cout << "Todo: computeUnitigNodes: a partir de snodes pas des reads, besoin de randomiser les nodes" << endl;

	
	string readFilename = _outputDir + "/read_data_init.txt";
	if(_useCorrectedRead){
		readFilename = _outputDir + "/read_data_corrected.txt";
	}

	
	KminmerParserParallel parser(readFilename, _minimizerSize, _kminmerSize, false, false, _nbCores);
	parser.parse(ComputeUnitigFunctor(*this));

	if(!_isFirstPass){
		//const string& filename_uncorrectedReads = _outputDir + "/read_uncorrected.txt";
		const string& contigFilename = _outputDir + "/unitig_data.txt";


		KminmerParserParallel parser2(contigFilename, _minimizerSize, _kminmerSize, false, false, _nbCores);
		parser2.parse(ComputeUnitigFunctor(*this));
		
	}


}
*/


void CreateMdbg::computeUnitigNodes(){

	ComputeUnitigFunctor computeUnitigFunctor(*this);


	//cout << "computeUnitigNodes: besoin d'activer le if(_graph->_isNodeIndexIndexed[nodeIndex]) exist = true (requiert d'indexer tous les kmervec des unitigs)" << endl;
	//cout << "Parallelisation unitig creation: bien faire en sorte que les nodeSource soit randomiser, surtout ne pas paralleliser a partir des node d'un read successif qui vont probablement etre dans le meme unitig" << endl;
	//cout << "nb cores 1" << endl;
	


	_unitigIndex = 0;

	ifstream kminmerFile(_outputDir + "/kminmerData_min.txt");

	//cout << "ComputeUnitig nb cores 1 " << endl;

	#pragma omp parallel num_threads(_nbCores)
	{

		bool isEOF = false;
		KmerVec vec;

		while (true) {

			#pragma omp critical(indexEdge)
			{


				vector<MinimizerType> minimizers;
				isEOF = MDBG::readKminmer(minimizers, kminmerFile, _kminmerSize);
				
				vec._kmers = minimizers;

				//vector<MinimizerType> minimizerSeq;
				//minimizerSeq.resize(_kminmerSize);
				//kminmerFile.read((char*)&minimizerSeq[0], minimizerSeq.size()*sizeof(MinimizerType));

				//isEOF = kminmerFile.eof();

				
				//u_int32_t abundance;
				//u_int32_t quality;

				//if(!isEOF){
				//	kminmerFile.read((char*)&nodeName, sizeof(nodeName));
				//	kminmerFile.read((char*)&abundance, sizeof(abundance));
					//kminmerFile.read((char*)&quality, sizeof(quality));
				//	vec._kmers = minimizerSeq;
					
					//cout << endl << "Start: " << vec.toString() << " " << nodeName << endl;
				//}

			}

			if(isEOF) break;
			
			//cout << "1" << endl;
			//computeUnitigFunctor.computeUnitigNode(vec);
			//cout << "2" << endl;
			computeUnitigFunctor.computeUnitigNode2(vec);
			
			//computeUnitigNode(vec);

		}
	}

	kminmerFile.close();


}





void CreateMdbg::getSuccessors(const KmerVec& vec, vector<KmerVec>& successors){

	successors.clear();
	//KmerVec vec = keys[i];
	KmerVec vec_rev = vec.reverse();

	KmerVec vec_rev_suffix = vec_rev.suffix();
	KmerVec vec_suffix = vec.suffix();
	
	KmerVec suffix = vec.suffix().normalize();

	if(_mdbgEdges3.find(suffix.hash128()) == _mdbgEdges3.end()){
		cout << "derp succ" << endl;
		return;
	}
		
	//cout << "\tRofl succ: " << _mdbgEdges3[suffix.hash128()].size() << endl;

	
	for(KminmerEdge3& edge : _mdbgEdges3[suffix.hash128()]){
		
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
		
		if(v_rev==vec) continue;

		//u_int32_t id2 = edge._nodeName;

		if (vec_suffix == v.prefix()) {
			successors.push_back(v);

			//dumpEdge(id, edgePlus, id2, edgePlus);
		}
		if (vec_suffix == v_rev.prefix()) {
			successors.push_back(v_rev);
			//dumpEdge(id, edgePlus, id2, edgeMinus);
		}
		if (vec_rev_suffix == v.prefix()) {
			successors.push_back(v_rev);
			//dumpEdge(id, edgeMinus, id2, edgePlus);
		}
		if (vec_rev_suffix == v_rev.prefix()) {
			successors.push_back(v);
			//dumpEdge(id, edgeMinus, id2, edgeMinus);
		}
	}

}

void CreateMdbg::getPredecessors(const KmerVec& vec, vector<KmerVec>& predecessors){


	predecessors.clear();
	//KmerVec vec = keys[i];
	KmerVec vec_rev = vec.reverse();

	KmerVec vec_rev_suffix = vec_rev.suffix();
	KmerVec vec_suffix = vec.suffix();
	
	const KmerVec prefix = vec.prefix().normalize();

	
	//cout << "\tRofl preds: " << _mdbgEdges3[prefix.hash128()].size() << endl;

	for(KminmerEdge3& edge : _mdbgEdges3[prefix.hash128()]){

		//cout << "\tPred: " << edge._isReversed << " " << edge._minimizer << endl;

		KmerVec v = prefix;
		if(edge._isReversed){
			v = v.reverse();
		}
		
		
		if(edge._isPrefix){
			v._kmers.push_back(edge._minimizer);
		}
		else{
			v._kmers.insert(v._kmers.begin(), edge._minimizer);
		}
		

		if(v==vec) continue;

		KmerVec v_rev = v.reverse();
		
		if(v_rev==vec) continue;

		//cout << "\tALlo? " << v.toString() << endl;
		//cout << "\tALlo? " << v_rev.toString() << endl;

		//u_int32_t id2 = edge._nodeName;

		if (vec_suffix == v.prefix()) {
			//cout << "rofl 1 " << " " << v.toString() << endl;
			predecessors.push_back(v);
			//dumpEdge(id, edgePlus, id2, edgePlus);
		}
		if (vec_suffix == v_rev.prefix()) {
			//cout << "rofl 2 " << v_rev.toString() << endl;
			predecessors.push_back(v_rev);
		}
		if (vec_rev_suffix == v.prefix()) {
			//cout << "rofl 3 " << v_rev.toString() << endl;
			predecessors.push_back(v_rev);
			//dumpEdge(id, edgeMinus, id2, edgePlus);
		}
		if (vec_rev_suffix == v_rev.prefix()) {
			//cout << "rofl 4 " << v.toString() << endl;
			predecessors.push_back(v);
			//dumpEdge(id, edgeMinus, id2, edgeMinus);
			//predecessors.push_back(v_rev);
		}

	}


}

bool CreateMdbg::hasSingleSuccessor(const KmerVec& vec){

	int nbSuccessors = 0;

	KmerVec vec_rev = vec.reverse();
	KmerVec vec_rev_suffix = vec_rev.suffix();
	KmerVec vec_suffix = vec.suffix();

	KmerVec suffix = vec.suffix().normalize();
	
	for(KminmerEdge3& edge : _mdbgEdges3[suffix.hash128()]){
		
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
		

		if(v==vec) continue;

		KmerVec v_rev = v.reverse();
		
		if(v_rev==vec) continue;

		if (vec_suffix == v.prefix()) {
			nbSuccessors += 1;
		}
		if (vec_suffix == v_rev.prefix()) {
			nbSuccessors += 1;
		}
		if (vec_rev_suffix == v.prefix()) {
			nbSuccessors += 1;
		}
		if (vec_rev_suffix == v_rev.prefix()) {
			nbSuccessors += 1;
		}

		if(nbSuccessors > 1) return false;
	}

	return nbSuccessors == 1;

}


bool CreateMdbg::hasSinglePredecessor(const KmerVec& vec){


	int nbPredecessors = 0;

	KmerVec vec_rev = vec.reverse();

	
	const KmerVec prefix = vec.prefix().normalize();

	KmerVec vec_rev_suffix = vec_rev.suffix();
	KmerVec vec_suffix = vec.suffix();
	
	for(KminmerEdge3& edge : _mdbgEdges3[prefix.hash128()]){

		KmerVec v = prefix;
		if(edge._isReversed){
			v = v.reverse();
		}
		
		
		if(edge._isPrefix){
			v._kmers.push_back(edge._minimizer);
		}
		else{
			v._kmers.insert(v._kmers.begin(), edge._minimizer);
		}
		

		if(v==vec) continue;

		KmerVec v_rev = v.reverse();
		
		if(v_rev==vec) continue;


		if (vec_suffix == v.prefix()) {
			nbPredecessors += 1;
		}
		if (vec_suffix == v_rev.prefix()) {
			nbPredecessors += 1;
		}
		if (vec_rev_suffix == v.prefix()) {
			nbPredecessors += 1;
		}
		if (vec_rev_suffix == v_rev.prefix()) {
			nbPredecessors += 1;
		}

		if(nbPredecessors > 1) return false;
	}

	return nbPredecessors == 1;

}


void CreateMdbg::getSuccessors_unitig(const KmerVec& vec, UnitigType unitigIndexFrom, vector<UnitigEdge>& successors){

	//bool isPalindomre = vec.isPalindrome();
	//if(isPalindomre){
	//	cout << "palouf" << endl;
	//	cout << vec.toString() << endl;
	//}

	//cout << "\tgetSuccessors_unitig" << endl;
	successors.clear();
	//KmerVec vec = keys[i];
	KmerVec vec_rev = vec.reverse();

	static u_int8_t isS = 0;
	static u_int8_t isL = 1;
	static u_int8_t edgePlus = 0;
	static u_int8_t edgeMinus = 1;
	
	KmerVec suffix = vec.suffix().normalize();

	if(_mdbgEdges4.find(suffix.hash128()) == _mdbgEdges4.end()){
		cout << "derp succ" << endl;
		return;
	}
		
	//cout << "\tRofl succ: " << _mdbgEdges3[suffix.hash128()].size() << endl;

	
	for(KminmerEdge4& edge : _mdbgEdges4[suffix.hash128()]){
		
		//cout << "\t\t" << suffix.toString() << " " << edge._isPrefix << " " << edge._isReversed << " " << edge._minimizer << " " << edge._unitigIndex << endl;
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
			
		//if(isPalindomre){
		//	cout << "palouf" << endl;
			//cout << suffix.toString() << "    " << v.toString() << endl;
		//}


		if(v==vec) continue;

		KmerVec v_rev = v.reverse();
		
		if(v_rev==vec) continue;

		//u_int32_t id2 = edge._nodeName;

		if (vec.suffix() == v.prefix()) {
			//if(isPalindomre){
			//	cout << "suf 1 " << edge._unitigIndex << endl;	
			//}
			//cout << "\t\tSucc1: " << edge._isReversed << " " << edge._unitigIndex << endl;
			successors.push_back({unitigIndexFrom, edge._unitigIndex});

			//dumpEdge(id, edgePlus, id2, edgePlus);
		}
		if (vec.suffix() == v_rev.prefix()) {
			//if(isPalindomre){
			//	cout << "suf 2 " << edge._unitigIndex << endl;	
			//	cout << vec.toString() << "         " <<  v_rev.toString() << endl;
			//}
			//cout << "\t\tSucc2: " << edge._isReversed << " " << edge._unitigIndex << endl;
			successors.push_back({unitigIndexFrom, Utils::unitigIndexToReverseDirection(edge._unitigIndex)});
			//if(edge._unitigIndex % 2 == 0){
			//	successors.push_back(edge._unitigIndex+1);
			//}
			//else{
			//	successors.push_back(edge._unitigIndex-1);
			//}
			//dumpEdge(id, edgePlus, id2, edgeMinus);
		}
		if (vec_rev.suffix() == v.prefix()) {
			//if(isPalindomre){
			//	cout << "suf 3 " << edge._unitigIndex << endl;	
			//}
			//cout << "\t\tSucc3: " << edge._isReversed << " " << edge._unitigIndex << endl;
			successors.push_back({Utils::unitigIndexToReverseDirection(unitigIndexFrom), edge._unitigIndex});
			//successors.push_back(edge._unitigIndex);
			//if(edge._unitigIndex % 2 == 0){
			//	successors.push_back(edge._unitigIndex+1);
			//}
			//else{
			//	successors.push_back(edge._unitigIndex-1);
			//}
			//successors.push_back(v_rev);
			//dumpEdge(id, edgeMinus, id2, edgePlus);
		}
		if (vec_rev.suffix() == v_rev.prefix()) {
			//if(isPalindomre){
			//	cout << "suf 4 " << edge._unitigIndex << endl;	
			//	cout << vec_rev.toString() << "         " <<  v_rev.toString() << endl;
			//}
			//cout << "\t\tSucc4: " << edge._isReversed << " " << edge._unitigIndex << endl;
			successors.push_back({Utils::unitigIndexToReverseDirection(unitigIndexFrom), Utils::unitigIndexToReverseDirection(edge._unitigIndex)});
			//dumpEdge(id, edgeMinus, id2, edgeMinus);
		}
	}

}

void CreateMdbg::getPredecessors_unitig(const KmerVec& vec, UnitigType unitigIndexFrom, vector<UnitigEdge>& predecessors){

	//bool isPalindomre = vec.isPalindrome();
	//if(isPalindomre){
	//	cout << "palouf" << endl;
	//}

	//cout << "\tgetPredecessors_unitig" << endl;

	predecessors.clear();
	//KmerVec vec = keys[i];
	KmerVec vec_rev = vec.reverse();

	static u_int8_t isS = 0;
	static u_int8_t isL = 1;
	static u_int8_t edgePlus = 0;
	static u_int8_t edgeMinus = 1;
	
	const KmerVec prefix = vec.prefix().normalize();
	//cout << "\t" << prefix.toString() << endl;

	if(_mdbgEdges4.find(prefix.hash128()) == _mdbgEdges4.end()){
		//cout << "derp pred" << endl;
		return;
	}
	
	//cout << "\tRofl preds: " << _mdbgEdges3[prefix.hash128()].size() << endl;

	for(KminmerEdge4& edge : _mdbgEdges4[prefix.hash128()]){

		//cout << "\t\t" << prefix.toString() << " " << edge._isPrefix << " " << edge._isReversed << " " << edge._minimizer << " " << edge._unitigIndex << endl;
		//cout << "\tPred: " << edge._isReversed << " " << edge._minimizer << endl;

		KmerVec v = prefix;
		if(edge._isReversed){
			v = v.reverse();
		}
		
		
		if(edge._isPrefix){
			v._kmers.push_back(edge._minimizer);
		}
		else{
			v._kmers.insert(v._kmers.begin(), edge._minimizer);
		}
		

		if(v==vec) continue;

		KmerVec v_rev = v.reverse();
		
		if(v_rev==vec) continue;

		//cout << "\tALlo? " << v.toString() << endl;
		//cout << "\tALlo? " << v_rev.toString() << endl;

		//u_int32_t id2 = edge._nodeName;

		if (vec.suffix() == v.prefix()) {
			//if(isPalindomre){
			//	cout << "pre 1" << endl;	
			//}
			//cout << "\t\tPred1: " << edge._isReversed << " " << edge._unitigIndex << endl;
			//cout << "rofl 1 " << " " << v.toString() << endl;
			predecessors.push_back({unitigIndexFrom, edge._unitigIndex});
			//dumpEdge(id, edgePlus, id2, edgePlus);
		}
		if (vec.suffix() == v_rev.prefix()) {
			//if(isPalindomre){
			//	cout << "pre 2" << endl;	
			//}
			//cout << "\t\tPred2: " << edge._isReversed << " " << edge._unitigIndex << endl;
			//predecessors.push_back(edge._unitigIndex);
			predecessors.push_back({unitigIndexFrom, Utils::unitigIndexToReverseDirection(edge._unitigIndex)});
			//if(edge._unitigIndex % 2 == 0){
			//	predecessors.push_back(edge._unitigIndex+1);
			//}
			//else{
			//	predecessors.push_back(edge._unitigIndex-1);
			//}
			//cout << "rofl 2 " << v_rev.toString() << endl;
			//predecessors.push_back(v_rev);
		}
		if (vec_rev.suffix() == v.prefix()) {
			//if(isPalindomre){
			//	cout << "pre 3" << endl;	
			//}
			//cout << "\t\tPred3: " << edge._isReversed << " " << edge._unitigIndex << endl;
			predecessors.push_back({Utils::unitigIndexToReverseDirection(unitigIndexFrom), edge._unitigIndex});
			//if(edge._unitigIndex % 2 == 0){
			//	predecessors.push_back(edge._unitigIndex+1);
			//}
			//else{
			//	predecessors.push_back(edge._unitigIndex-1);
			//}
			//cout << "rofl 3 " << v_rev.toString() << endl;
			//predecessors.push_back(v_rev);
			//dumpEdge(id, edgeMinus, id2, edgePlus);
		}
		if (vec_rev.suffix() == v_rev.prefix()) {
			//if(isPalindomre){
			//	cout << "pre 4" << endl;	
			//}
			//cout << "\t\tPred4: " << edge._isReversed << " " << edge._unitigIndex << endl;
			//cout << "rofl 4 " << v.toString() << endl;
			predecessors.push_back({Utils::unitigIndexToReverseDirection(unitigIndexFrom), Utils::unitigIndexToReverseDirection(edge._unitigIndex)});
			//dumpEdge(id, edgeMinus, id2, edgeMinus);
			//predecessors.push_back(v_rev);
		}

	}


}

void CreateMdbg::dumpUnitigNode(const UnitigType& unitigIndex, const vector<MinimizerType>& unitig){
	/*
	if(_lalaLol % 2 == 0){
		_lalaLol2 = unitig;
	}
	else{


		cout << endl << "----" << endl;
		for(MinimizerType minimizer : _lalaLol2){
			cout << "\t" << minimizer << endl;
		}

		cout << endl ;
		for(MinimizerType minimizer : unitig){
			cout << "\t" << minimizer << endl;
		}


		if(_lalaLol2 == unitig){

		}
		else {
			cout << "derp" << endl;
			getchar();
		}
	}

	_lalaLol += 1;
	//cout << "dump node: " << fromUnitigIndex << " " << toUnitigIndex << endl;
	//cout << "Dump unitig: utg" << unitigIndex << " " <<  (unitig.size()-3) << endl;


		
	//if(unitigIndex == 4540){
	//	getchar();
	//}
	*/

	u_int32_t size = unitig.size();
	
	_unitigGraphFile_nodes.write((const char*)&size, sizeof(size));
	_unitigGraphFile_nodes.write((const char*)&unitig[0], size * sizeof(MinimizerType));
	_unitigGraphFile_nodes.write((const char*)&unitigIndex, sizeof(unitigIndex));
	//_unitigGraphFile_nodes.write((const char*)&nodeName, sizeof(nodeName));



}

void CreateMdbg::dumpUnitigEdge(u_int32_t fromUnitigIndex, u_int32_t toUnitigIndex, bool isSuccessor){


	#pragma omp critical(dumpUnitigEdge)
	{

		//if(fromUnitigIndex == 1916 && toUnitigIndex == 4914){
		//	cout << "Add edge: " << fromUnitigIndex << " -> " << toUnitigIndex << endl;
		//	cout << "derp" << endl;
		//	getchar(); 
		//}
        //if(fromUnitigIndex == 1916 || toUnitigIndex == 1916 || fromUnitigIndex == 1917 || toUnitigIndex == 1917){
        //    cout << "Add edge: " << fromUnitigIndex << " -> " << toUnitigIndex << endl;
        //}

		//getchar();

		if(isSuccessor){
			_unitigGraphFile_edges_successors.write((const char*)&fromUnitigIndex, sizeof(fromUnitigIndex));
			_unitigGraphFile_edges_successors.write((const char*)&toUnitigIndex, sizeof(toUnitigIndex));
		}
		else{
			_unitigGraphFile_edges_predecessors.write((const char*)&fromUnitigIndex, sizeof(fromUnitigIndex));
			_unitigGraphFile_edges_predecessors.write((const char*)&toUnitigIndex, sizeof(toUnitigIndex));
		}
	}

}

void CreateMdbg::computeEdgesOld(){
	
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
			
			computeEdgeOld(vec, nodeName);

		}
	}

	kminmerFile.close();
}

void CreateMdbg::computeEdgeOld(const KmerVec& vec, u_int32_t id){


	static u_int8_t isS = 0;
	static u_int8_t isL = 1;
	static u_int8_t edgePlus = 0;
	static u_int8_t edgeMinus = 1;

	//KmerVec vec = keys[i];
	KmerVec vec_rev = vec.reverse();

	const KmerVec prefix = vec.prefix().normalize();


	for(KminmerEdge3& edge : _mdbgEdges3[prefix.hash128()]){

			
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


		if(v==vec) continue;

		KmerVec v_rev = v.reverse();
		
		if(v_rev==vec) continue;

		u_int32_t id2;// = edge._nodeName;

		if (vec.suffix() == v.prefix()) {
			dumpEdgeOld(id, edgePlus, id2, edgePlus);
		}
		if (vec.suffix() == v_rev.prefix()) {
			dumpEdgeOld(id, edgePlus, id2, edgeMinus);
		}
		if (vec_rev.suffix() == v.prefix()) {
			dumpEdgeOld(id, edgeMinus, id2, edgePlus);
		}
		if (vec_rev.suffix() == v_rev.prefix()) {
			dumpEdgeOld(id, edgeMinus, id2, edgeMinus);
		}

	}

	
	KmerVec suffix = vec.suffix().normalize();
	for(KminmerEdge3& edge : _mdbgEdges3[suffix.hash128()]){
		
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


		if(v==vec) continue;

		KmerVec v_rev = v.reverse();
		
		if(v_rev==vec) continue;

		u_int32_t id2;// = edge._nodeName;

		if (vec.suffix() == v.prefix()) {
			dumpEdgeOld(id, edgePlus, id2, edgePlus);
		}
		if (vec.suffix() == v_rev.prefix()) {
			dumpEdgeOld(id, edgePlus, id2, edgeMinus);
		}
		if (vec_rev.suffix() == v.prefix()) {
			dumpEdgeOld(id, edgeMinus, id2, edgePlus);
		}
		if (vec_rev.suffix() == v_rev.prefix()) {
			dumpEdgeOld(id, edgeMinus, id2, edgeMinus);
		}
	}

}

void CreateMdbg::dumpEdgeOld(u_int32_t nodeNameFrom, u_int8_t nodeNameFromOri, u_int32_t nodeNameTo, u_int8_t nodeNameToOri){
	
	#pragma omp critical(gfa)
	{
		
		_nbEdges += 1;
		//_outputFileGfa.write((const char*)&nodeNameFrom, sizeof(nodeNameFrom));
		//_outputFileGfa.write((const char*)&nodeNameFromOri, sizeof(nodeNameFromOri));
		//_outputFileGfa.write((const char*)&nodeNameTo, sizeof(nodeNameTo));
		//_outputFileGfa.write((const char*)&nodeNameToOri, sizeof(nodeNameToOri));
	}
}




void CreateMdbg::computeUnitiStartEndNode(){

	//unordered_map<UnitigType, u_int32_t> unitigIndex_to_deterministicStartNode;
	vector<KmerVecSorterData> kmerVecs;


	ofstream startEndFile(_outputDir + "/unitigGraph.startEnd.bin");
	ifstream nodeFile(_outputDir + "/unitigGraph.nodes.bin");

	while(true){

		u_int32_t size;
		nodeFile.read((char*)&size, sizeof(size));
		

		if(nodeFile.eof()) break;

		vector<MinimizerType> minimizers;
		minimizers.resize(size);
		nodeFile.read((char*)&minimizers[0], size * sizeof(MinimizerType));


		UnitigType unitigIndex;
		nodeFile.read((char*)&unitigIndex, sizeof(unitigIndex));
		

		const vector<KmerVec>& nodes = Utils::minimizersToKminmers(minimizers, _kminmerSize);

		KmerVec startNode = nodes[0];
		KmerVec endNode = nodes[nodes.size()-1];
		KmerVec startNodeRC = endNode.reverse();
		KmerVec endNodeRC = startNode.reverse();
		
		kmerVecs.push_back({unitigIndex, startNode._kmers});
		kmerVecs.push_back({unitigIndex+1, startNodeRC._kmers});
		//indexUnitigEdge(unitigIndex, unitig);

		bool isReversed;
		startNode = startNode.normalize(isReversed);
		u_int32_t nodeName = _kmervec_to_nodeName[startNode];
		if(isReversed){
			nodeName = nodeName*2 + 1;
		}
		else{
			nodeName = nodeName*2;
		}
		

		startEndFile.write((const char*)&unitigIndex, sizeof(unitigIndex));
		startEndFile.write((const char*)&nodeName, sizeof(nodeName));

		unitigIndex += 1;

		startNodeRC = startNodeRC.normalize(isReversed);
		nodeName = _kmervec_to_nodeName[startNodeRC];
		if(isReversed){
			nodeName = nodeName*2 + 1;
		}
		else{
			nodeName = nodeName*2;
		}
		//cout << nodeName << " " << nodeNamerc << endl;

		startEndFile.write((const char*)&unitigIndex, sizeof(unitigIndex));
		startEndFile.write((const char*)&nodeName, sizeof(nodeName));


	}

	nodeFile.close();

	startEndFile.close();
	/*

	std::sort(kmerVecs.begin(), kmerVecs.end(), KmerVecComparator);
	u_int32_t startNode = 0;
	
	ofstream startEndFile(_outputDir + "/unitigGraph.startEnd.bin");

	for(const KmerVecSorterData& data : kmerVecs){

		
		startEndFile.write((const char*)&data._unitigIndex, sizeof(data._unitigIndex));
		startEndFile.write((const char*)&startNode, sizeof(startNode));

		startNode += 1;
	}




	startEndFile.close();
	*/

}

void CreateMdbg::indexUnitigEdges(){



	ifstream nodeFile(_outputDir + "/unitigGraph.nodes.bin");

	#pragma omp parallel num_threads(1)
	{

		bool isEOF = false;
		vector<MinimizerType> unitig;
		UnitigType unitigIndex;
		//u_int32_t nodeName;

		while (true) {

			#pragma omp critical(indexEdge)
			{

				u_int32_t size;
				nodeFile.read((char*)&size, sizeof(size));

				isEOF = nodeFile.eof();

				

				if(!isEOF){
					unitig.resize(size);
					nodeFile.read((char*)&unitig[0], size * sizeof(MinimizerType));


					nodeFile.read((char*)&unitigIndex, sizeof(unitigIndex));
					//nodeFile.read((char*)&nodeName, sizeof(nodeName));
				}

				//cout << nodeName << endl;
			}

			if(isEOF) break;
			
			indexUnitigEdge(unitigIndex, unitig);

		}
	}

	nodeFile.close();


}
	

void CreateMdbg::indexUnitigEdge(const u_int32_t unitigIndex, const vector<MinimizerType>& minimizers){

	//cout << "Index unitig edge: utg" << unitigIndex << endl;
	//for(u_int64_t minimizer : minimizers){
	//	cout << "\t" << minimizer << endl;
	//}

	u_int32_t unitigIndexRC = unitigIndex + 1;
	const vector<KmerVec>& nodes = Utils::minimizersToKminmers(minimizers, _kminmerSize);

	const KmerVec& startNode = nodes[0];
	const KmerVec& endNode = nodes[nodes.size()-1];
	const KmerVec& startNodeRC = endNode.reverse();
	const KmerVec& endNodeRC = startNode.reverse();

	


	
	//cout << "\tutg" << unitigIndex << " " << minimizers.size() << endl;
	//return;
	

	bool isReversed;
	const KmerVec& startNodeNorm = startNode.normalize(isReversed);
	//cout << "\tStart node: " << startNodeNorm.toString() << " " << isReversed << endl;

	//const KmerVec& endNodeNorm = endNode.normalize(isReversed);

	//indexEdgeUnitig(startNodeNorm, unitigIndex);
	//indexEdgeUnitig(endNodeNorm, unitigIndex);
	
	
	if(isReversed){
		indexEdgeUnitig(startNodeNorm, unitigIndexRC);
	}
	else{
		indexEdgeUnitig(startNodeNorm, unitigIndex);
	}

	if(startNode == endNode) return;
	
	const KmerVec& endNodeNorm = endNode.normalize(isReversed);

	//cout << "\tEnd node: " << endNodeNorm.toString() << " " << isReversed << endl;

	if(isReversed){
		indexEdgeUnitig(endNodeNorm, unitigIndexRC);
	}
	else{
		indexEdgeUnitig(endNodeNorm, unitigIndex);
	}
	
	//getchar();

	/*
	indexEdge(startNode, 0);

	if(startNode == endNode){
	}
	else{
		indexEdge(endNode, 0);
	}

	*/

	//_startNode_to_unitigIndex[startNode] = unitigIndex;
	//_endNode_to_unitigIndex[endNode] = unitigIndex;
	//_startNode_to_unitigIndex[startNodeRC] = unitigIndexRC;
	//_endNode_to_unitigIndex[endNodeRC] = unitigIndexRC;


	/*
	_startNode_to_unitigIndex.lazy_emplace_l(startNode, 
	[&unitigIndex](UnitigEdgeMap::value_type& v) { // key exist
		v.second = unitigIndex;
		//v.second.push_back({prefix_minimizer, isReversedPrefix, true, false});
	},           
	[&startNode, &unitigIndex](const UnitigEdgeMap::constructor& ctor) { // key inserted
		ctor(startNode, unitigIndex); 
	});

	_startNode_to_unitigIndex.lazy_emplace_l(startNodeRC, 
	[&unitigIndexRC](UnitigEdgeMap::value_type& v) { // key exist
		v.second = unitigIndexRC;
		//v.second.push_back({prefix_minimizer, isReversedPrefix, true, false});
	},           
	[&startNodeRC, &unitigIndexRC](const UnitigEdgeMap::constructor& ctor) { // key inserted
		ctor(startNodeRC, unitigIndexRC); 
	});

	_endNode_to_unitigIndex.lazy_emplace_l(endNode, 
	[&unitigIndex](UnitigEdgeMap::value_type& v) { // key exist
		v.second = unitigIndex;
		//v.second.push_back({prefix_minimizer, isReversedPrefix, true, false});
	},           
	[&endNode, &unitigIndex](const UnitigEdgeMap::constructor& ctor) { // key inserted
		ctor(endNode, unitigIndex); 
	});

	_endNode_to_unitigIndex.lazy_emplace_l(endNodeRC, 
	[&unitigIndexRC](UnitigEdgeMap::value_type& v) { // key exist
		v.second = unitigIndexRC;
		//v.second.push_back({prefix_minimizer, isReversedPrefix, true, false});
	},           
	[&endNodeRC, &unitigIndexRC](const UnitigEdgeMap::constructor& ctor) { // key inserted
		ctor(endNodeRC, unitigIndexRC); 
	});
	*/
	//cout << "index " << _startNode_to_unitigIndex.size() << " " << _endNode_to_unitigIndex.size() << endl;
	//_startNode_to_unitigIndex[startNode] = _graph._unitigIndex;
	//_endNode_to_unitigIndex[endNode] = _graph._unitigIndex;
	//_startNode_to_unitigIndex[startNodeRC] = _graph._unitigIndex+1;
	//_endNode_to_unitigIndex[endNodeRC] = _graph._unitigIndex+1;

}
	

void CreateMdbg::indexEdgeUnitig(const KmerVec& vec, u_int32_t unitigIndex){


	bool isReversedPrefix;
	bool isReversedSuffix;
	KmerVec prefix = vec.prefix().normalize(isReversedPrefix);
	KmerVec suffix = vec.suffix().normalize(isReversedSuffix);
	MinimizerType prefix_minimizer = vec._kmers[vec._kmers.size()-1];
	MinimizerType suffix_minimizer = vec._kmers[0];
	
	//cout << isReversedPrefix << " " << isReversedSuffix << endl;
	//getchar();
	_mdbgEdges4.lazy_emplace_l(prefix.hash128(), 
	[&prefix_minimizer, &isReversedPrefix, &unitigIndex](MdbgEdgeMap4::value_type& v) { // key exist
		v.second.push_back({prefix_minimizer, isReversedPrefix, true, unitigIndex});
	},           
	[&prefix, &prefix_minimizer, &isReversedPrefix, &unitigIndex](const MdbgEdgeMap4::constructor& ctor) { // key inserted
		vector<KminmerEdge4> nodes;
		nodes.push_back({prefix_minimizer, isReversedPrefix, true, unitigIndex});
		ctor(prefix.hash128(), nodes); 
	});

	_mdbgEdges4.lazy_emplace_l(suffix.hash128(), 
	[&suffix_minimizer, &isReversedSuffix, &unitigIndex](MdbgEdgeMap4::value_type& v) { // key exist
		v.second.push_back({suffix_minimizer, isReversedSuffix, false, unitigIndex});
	},           
	[&suffix, &suffix_minimizer, &isReversedSuffix, &unitigIndex](const MdbgEdgeMap4::constructor& ctor) { // key inserted
		vector<KminmerEdge4> nodes;
		nodes.push_back({suffix_minimizer, isReversedSuffix, false, unitigIndex});
		ctor(suffix.hash128(), nodes); 
	});
	
}



void CreateMdbg::computeUnitigEdges(){


	
	ifstream nodeFile(_outputDir + "/unitigGraph.nodes.bin");

	#pragma omp parallel num_threads(_nbCores)
	{

		bool isEOF = false;
		vector<MinimizerType> unitig;
		UnitigType unitigIndex;
		//u_int32_t nodeName;

		while (true) {

			#pragma omp critical(indexEdge)
			{

				u_int32_t size;
				nodeFile.read((char*)&size, sizeof(size));

				isEOF = nodeFile.eof();

				

				if(!isEOF){
					unitig.resize(size);
					nodeFile.read((char*)&unitig[0], size * sizeof(MinimizerType));


					nodeFile.read((char*)&unitigIndex, sizeof(unitigIndex));
					//nodeFile.read((char*)&nodeName, sizeof(nodeName));
				}

				//cout << nodeName << endl;
			}

			if(isEOF) break;
			
			computeUnitigEdge(unitigIndex, unitig);

		}
	}

	nodeFile.close();
	

	/*
	int nbEdges  =0 ;
	for(const auto& it : _endNode_to_unitigIndex){
		
		const KmerVec& endNode = it.first;
		u_int32_t unitigIndex = it.second;


		vector<KmerVec> successors;
		getSuccessors(endNode, successors);
			

		for(KmerVec node : successors){
			u_int32_t toUnitigIndex = _startNode_to_unitigIndex[node];
			dumpUnitigEdge(unitigIndex, toUnitigIndex, true);
			//cout << "Edge: utg" << unitigIndex << " utg" << toUnitigIndex << "    " << endNode.toString() << "    " << node.toString() << endl;
			nbEdges += 1;
		}
		
	}

	cout << "Nb unitig edges: " << nbEdges << endl;
	*/

}


void CreateMdbg::computeUnitigEdge(const u_int32_t unitigIndex, const vector<MinimizerType>& minimizers){

	
	u_int32_t unitigIndexRC = unitigIndex + 1;
	const vector<KmerVec>& nodes = Utils::minimizersToKminmers(minimizers, _kminmerSize);

	
	//cout << "\tutg" << unitigIndex << " " << minimizers.size() << endl;
	//return;
	
	const KmerVec& startNode = nodes[0];
	const KmerVec& endNode = nodes[nodes.size()-1];
	//const KmerVec& startNodeRC = endNode.reverse();
	//const KmerVec& endNodeRC = startNode.reverse();

	//const KmerVec& startNodeNorm = startNode.normalize(isReversed);

	bool isReversed;
	const KmerVec& endNodeNorm = endNode; //.normalize(isReversed);
	
	//cout << "computeUnitigEdge_successors " << unitigIndex << endl;
	vector<UnitigEdge> successors;
	getSuccessors_unitig(endNodeNorm, unitigIndex, successors);

	for(const UnitigEdge& edge : successors){

		//cout << " \t" << toUnitigIndex << endl;
		//u_int32_t unitigIndexFrom = unitigIndex;
		/*

		if(isReversed){
			if(unitigIndex % 2 == 0){
				unitigIndexFrom = unitigIndex + 1;
			}
			else{
				unitigIndexFrom = unitigIndex - 1;
			}
		}
		*/
		dumpUnitigEdge(edge._unitigIndexFrom, edge._unitigIndexTo, true);
	}



	//cout << "computeUnitigEdge_predecessors " << unitigIndex << endl;
	const KmerVec& startNodeNorm = startNode;//.normalize(isReversed);

	vector<UnitigEdge> predecessors;
	getPredecessors_unitig(startNodeNorm, unitigIndex, predecessors);

	for(const UnitigEdge& edge : predecessors){

		//cout << " \t" << toUnitigIndex << endl;
		//u_int32_t unitigIndexFrom = unitigIndex;


		dumpUnitigEdge(edge._unitigIndexFrom, edge._unitigIndexTo, true);
		//dumpUnitigEdge(Utils::unitigIndexToReverseDirection(unitigIndexFrom), Utils::unitigIndexToReverseDirection(toUnitigIndex), true);
		/*
		if(isReversed){
			if(unitigIndex % 2 == 0){
				unitigIndexFrom = unitigIndex + 1;
			}
			else{
				unitigIndexFrom = unitigIndex - 1;
			}
		}
		*/

		//dumpUnitigEdge(toUnitigIndex, unitigIndexFrom, true);
	}

	//if(unitigIndex == 10002) getchar();
	/*
	if(isReversed){
		indexEdgeUnitig(startNodeNorm, unitigIndexRC);
	}
	else{
		indexEdgeUnitig(startNodeNorm, unitigIndex);
	}

	const KmerVec& endNodeNorm = endNode.normalize(isReversed);

	if(isReversed){
		indexEdgeUnitig(endNodeNorm, unitigIndexRC);
	}
	else{
		indexEdgeUnitig(endNodeNorm, unitigIndex);
	}
	*/
}

void CreateMdbg::dumpUnitigAbundances(){

	_kminmerAbundances.clear();
	ifstream kminmerAbundanceFile(_outputDir + "/kminmerData_abundance.txt");

	//bool isEOF = false;
	//KmerVec vec;
	//u_int32_t nodeName;


	while (true) {

		u_int128_t vecHash;
		AbundanceType abundance;

		bool iseof = MDBG::readKminmerAbundance(vecHash, abundance, kminmerAbundanceFile);

		if(iseof) break;

		if(abundance == 1) continue;

		_kminmerAbundances[vecHash] = abundance;

	}

	kminmerAbundanceFile.close();


	ifstream nodeFile(_outputDir + "/unitigGraph.nodes.bin");

	while(true){

		u_int32_t size;
		nodeFile.read((char*)&size, sizeof(size));
		

		if(nodeFile.eof()) break;

		vector<MinimizerType> unitig;
		unitig.resize(size);
		nodeFile.read((char*)&unitig[0], size * sizeof(MinimizerType));


		UnitigType unitigIndex;
		u_int32_t nodeName;
		nodeFile.read((char*)&unitigIndex, sizeof(unitigIndex));
		//nodeFile.read((char*)&nodeName, sizeof(nodeName));

		vector<float> abundances;
		const vector<KmerVec>& kminmers = Utils::minimizersToKminmers(unitig, _kminmerSize);

		for(const KmerVec& kminmer : kminmers){

			const u_int128_t vecHash = kminmer.normalize().hash128();

			if(_kminmerAbundances.find(vecHash) == _kminmerAbundances.end()){
				abundances.push_back(1);
			}
			else{
				abundances.push_back(_kminmerAbundances[vecHash]);
			}
		}


		u_int32_t nbAbundances = abundances.size();
		_unitigGraphFile_nodes_abundances.write((const char*)&unitigIndex, sizeof(unitigIndex));
		_unitigGraphFile_nodes_abundances.write((const char*)&nbAbundances, sizeof(nbAbundances));
		_unitigGraphFile_nodes_abundances.write((const char*)&abundances[0], nbAbundances * sizeof(float));



	}

	nodeFile.close();
	_kminmerAbundances.clear();
}

void CreateMdbg::loadRefinedAbundances(){

	/*
	ifstream kminmerAbundanceFile(_outputDir + "/kminmerData_abundance_init.txt");

	while (true) {

		u_int128_t vecHash;
		AbundanceType abundance;

		bool iseof = MDBG::readKminmerAbundance(vecHash, abundance, kminmerAbundanceFile);

		if(iseof) break;

		if(abundance == 1) continue;
		_kminmerAbundances[vecHash] = abundance;


	}
	*/
	
	u_int64_t nbNodesOriginal = 0;
	u_int64_t sumAbundanceOriginal = 0;


	ifstream kminmerAbundanceFile(_outputDir + "/kminmerData_abundance.txt");
	//u_int64_t checksum = 0;

	while (true) {

		u_int128_t vecHash;
		AbundanceType abundance;

		bool iseof = MDBG::readKminmerAbundance(vecHash, abundance, kminmerAbundanceFile);

		if(iseof) break;

		//nbNodesOriginal += 1;
		sumAbundanceOriginal += abundance;

		//if(abundance == 1) continue;

		if(abundance == 1) continue;
		_kminmerAbundances[vecHash] = abundance;


	}

	kminmerAbundanceFile.close();
	
	//cout << "Nb nodes original: " << _kminmerAbundances.size() << endl;



	unordered_map<UnitigType, u_int32_t> unitigName_to_abundance;

	ifstream abundanceFile(_outputDir + "/unitigGraph.nodes.refined_abundances.bin");

	while (true) {

		UnitigType unitigName;
		u_int32_t abundance;

		abundanceFile.read((char*)&unitigName, sizeof(unitigName));
		
		if(abundanceFile.eof()) break;

		abundanceFile.read((char*)&abundance, sizeof(abundance));

		unitigName_to_abundance[unitigName] = abundance;
	}


	abundanceFile.close();




	ifstream nodeFile(_outputDir + "/unitigGraph.nodes.bin");

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


		const UnitigType& unitigName = unitigIndex/2;

		if(unitigName_to_abundance.find(unitigName) == unitigName_to_abundance.end()) continue;
		
		const u_int32_t unitigAbundance = unitigName_to_abundance[unitigName];
		//if(unitigAbundance == 1) continue;

		for(const KmerVec& kminmer : Utils::minimizersToKminmers(minimizers, _kminmerSizePrev)){

			const u_int128_t hash = kminmer.normalize().hash128();

			if(unitigAbundance == 1){
				if(_kminmerAbundances.find(hash) != _kminmerAbundances.end()){
					_kminmerAbundances.erase(hash);
				}
			}
			else{
				_kminmerAbundances[hash] = unitigAbundance;
			}
			//if(_kminmerAbundances.find(hash) == _kminmerAbundances.end()) continue;

			//cout << unitigAbundance << endl;
			
		}
	}

	nodeFile.close();
	
	/*
	u_int64_t nbNodes = 0;
	u_int64_t nbAbundantNodes = 0;
	u_int64_t sumAbundance = 0;
	u_int64_t sumAbundantAbundance = 0;



	for(const auto& it : _kminmerAbundances){

		nbNodes += 1;
		sumAbundance += it.second;

		if(it.second > 1){
			nbAbundantNodes += 1;
			sumAbundantAbundance += it.second;
		}
	}


	//cout << endl;
	//cout << "Nb nodes original: " << nbNodesOriginal << endl;
	cout << "Abundance original: " << sumAbundanceOriginal << endl;
	cout << endl;
	cout << "Nb nodes: " << nbNodes << endl;
	cout << "Nb nodes abundant: " << nbAbundantNodes << endl;
	cout << "Abundance: " << sumAbundance << endl;
	cout << "Abundance abundant: " << sumAbundantAbundance << endl;

	//getchar();

	//cout << "Nb nodes abundant: " << nbAbundantNodes << endl;
	//cout << "nb kminmers: " << nbKminmers << endl;
	//cout << "sum abundance: " << sumAbundance << endl;
	//getchar();

	*/
}
