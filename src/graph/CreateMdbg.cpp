

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
	//args::Flag arg_useCorrectedRead(parser, "", "Use corrected reads", {"corrected-read"});
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

	//_useBloomFilter = true;
	//if(arg_bf){
	//	_useBloomFilter = false;
	//}

	_isFirstPass = false;
	if(arg_firstPass){
		_isFirstPass = true;
	}

	//if(!_isFirstPass){
	//	_useBloomFilter = false;	
	//}

	//_useCorrectedRead = false;
	//if(arg_useCorrectedRead){
	//	_useCorrectedRead = true;
	//}

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


	/*
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
	*/

	_params.load(_outputDir + "/parameters.gz");
	_kminmerSize = _params._kminmerSize;
	_kminmerSizeFirst = _params._kminmerSizeFirst;
	_kminmerSizePrev = _params._kminmerSizePrev;

	openLogFile(_outputDir);

	/*
	Logger::get().debug() << "";
	Logger::get().debug() << "Contig filename: " << _filename_inputContigs;
	Logger::get().debug() << "Output dir: " << _outputDir;
	Logger::get().debug() << "Minimizer length: " << _minimizerSize;
	Logger::get().debug() << "Kminmer length: " << _kminmerSize;
	Logger::get().debug() << "Density: " << _minimizerDensity;
	Logger::get().debug() << "Min abundance: " << _minAbundance;
	Logger::get().debug() << "";
	*/
	//_inputDir = getInput()->get(STR_INPUT_DIR) ? getInput()->getStr(STR_INPUT_DIR) : "";
	//_input_extractKminmers= getInput()->get(STR_INPUT_EXTRACT_KMINMERS) ? getInput()->getStr(STR_INPUT_EXTRACT_KMINMERS) : "";

	//_filename_readMinimizers = _outputDir + "/read_data.gz";
	//_filename_contigMinimizers = _outputDir + "/contig_data.gz";
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



	_checksum_unitigNodes = 0;
	_checksum_unitigEdges = 0;
	_checksum_unitigAbundances = 0;
	//_lalaLol = 0;

	_nbUnitigEdges = 0;
	_nbUnitigNodes = 0;
	//_file_noKminmerReads = ofstream(_filename_noKminmerReads, std::ofstream::app);
	//_node_id = 0;
	//parseArgs();

	_mutexes.resize(1000);

	for(size_t i=0; i<_mutexes.size(); i++){
		omp_init_lock(&_mutexes[i]);
	}
	
	createMDBG();
	
	for(size_t i=0; i<_mutexes.size(); i++){
		omp_destroy_lock(&_mutexes[i]);
	}
	
}



void CreateMdbg::createMDBG (){


	_readStats.load(_outputDir + "/read_stats.txt");
	/*
	ifstream file_readStats(_outputDir + "/read_stats.txt");

	size_t totalNbReads;
	float minimizerDensity;
	u_int32_t n50ReadLength;
	float averageQuality;
	u_int64_t nbBases;
	u_int32_t meanReadLength;

	file_readStats.read((char*)&totalNbReads, sizeof(totalNbReads));
	file_readStats.read((char*)&n50ReadLength, sizeof(n50ReadLength));
	file_readStats.read((char*)&minimizerDensity, sizeof(minimizerDensity));
	file_readStats.read((char*)&nbBases, sizeof(nbBases));
	file_readStats.read((char*)&averageQuality, sizeof(averageQuality));
	file_readStats.read((char*)&meanReadLength, sizeof(meanReadLength));

	file_readStats.close();
	*/

	_nbPartitions = _readStats._nbBases / 20000000000ull;
	_nbPartitions = max(_nbPartitions, _nbCores);
	_nbPartitions = max(_nbPartitions, 1) ;//, 200);
	_nbPartitions = min(_nbPartitions, 5000);
	
	Logger::get().debug() << "Nb partitions: " << _nbPartitions;



	/*
	cout << "Checksum unitig nodes:   " << _checksum_unitigNodes << endl;
	cout << "Checksum unitig edges:   " << _checksum_unitigEdges << endl;
	cout << "Checksum unitig abundance:   " << _checksum_unitigAbundances << endl;

	createGfa();


	cout << "Checksum unitig nodes:   " << _checksum_unitigNodes << endl;
	cout << "Checksum unitig edges:   " << _checksum_unitigEdges << endl;
	cout << "Checksum unitig abundance:   " << _checksum_unitigAbundances << endl;

	exit(1);
	*/





	//_savingSmallContigs = false;


	//if(_useBloomFilter){
		//_bloomFilter = new BloomCacheCoherent<u_int64_t>(32000000000ull);
	//}

	_filename_smallContigs = _outputDir + "/smallContigs/smallContigs_k" + to_string(_kminmerSize) + ".bin";
	_fileSmallContigs = ofstream(_filename_smallContigs, std::ios_base::binary);
	_nbSmallContigs = 0;
	_nbKminmersTotal = 0;

	
	//cout << _bloomFilter->contains(5) << endl;
	//_bloomFilter->insert(5);
	//cout << _bloomFilter->contains(5) << endl;
	//while(true){};

	//_kminmerExist.clear();
	//_mdbg = new MDBG(_kminmerSize);
	//_mdbgNoFilter = new MDBG(_kminmerSize);
	//_mdbgSaved = new MDBG(_kminmerSize);



	//string inputFilename = _inputFilename;


	

	string inputFilename_min;
	//bool usePos = false;

	Logger::get().debug() << "Collecting kminmers";
	//auto start = high_resolution_clock::now();

	//_parsingContigs = false;
	

	if(_isFirstPass){

		//cout << "Counting kminmers" << endl;
		auto start = high_resolution_clock::now();
		
		_kminmerFile.open(_outputDir + "/kminmerData_min.txt");
		_kminmerAbundanceFile.open(_outputDir + "/kminmerData_abundance.txt");

		KminmerCounter kminmerCounter(*this);
		kminmerCounter.execute();
		_nbKminmersTotal += kminmerCounter._nbSolidKminmers;
		Logger::get().debug() << "\tDuration: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";
		//Logger::get().debug() << "\tNb kminmers: " << kminmerCounter._nbKminmers;
		Logger::get().debug() << "\tNb solid kminmers: " << kminmerCounter._nbSolidKminmers;

		_kminmerFile.close();
		_kminmerAbundanceFile.close();
		



		_kminmerFile.open(_outputDir + "/kminmerData_min.txt", std::ios_base::app);
		_kminmerAbundanceFile.open(_outputDir + "/kminmerData_abundance.txt", std::ios_base::app);

		Logger::get().debug() << "Rescuing kminmers";
		start = high_resolution_clock::now();
		
		if(_minAbundance <= 1){
			rescueKminmers();
		}

		_kminmerFile.close();
		_kminmerAbundanceFile.close();
		Logger::get().debug() << "\tDuration: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";
		Logger::get().debug() << "\tNb rescued kminmers: " << _nbRescuedKminmers;

		_nbKminmersTotal += _nbRescuedKminmers;

		Logger::get().debug() << "\tNb kminmers: " << _nbKminmersTotal;

		/*
		//exit(1);

		if(_useBloomFilter){
			//cout << "Attention: bloom filter exact" << endl;
			//_bloomFilterExact.clear();
			delete _bloomFilter;
				
			if(_minAbundance <= 1){
				//cout << "!! Bloom filter disabled!!!!!" << endl;
				KminmerParserParallel parser2(inputFilename_min, _minimizerSize, _kminmerSize, false, hasQuality, _nbCores);
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

			MdbgNodeMapLight().swap(_mdbgNodesLightUnique); //_mdbgNodesLightUnique.clear();
			//
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
		*/

	}
	else{
		
			
		Logger::get().debug() << "Loading refined abundances";
		auto start = high_resolution_clock::now();
		loadRefinedAbundances();
		Logger::get().debug() << "\tDuration: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";	


		if(_kminmerSize == _kminmerSizeFirst+1){
			



			auto start = high_resolution_clock::now();
			
			_kminmerFile.open(_outputDir + "/kminmerData_min.txt");
			_kminmerAbundanceFile.open(_outputDir + "/kminmerData_abundance.txt");

			KminmerCounter kminmerCounter(*this);
			kminmerCounter.execute();
			_nbKminmersTotal += kminmerCounter._nbSolidKminmers;
			Logger::get().debug() << "\tDuration: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";
			//cout << "\tNb kminmers: " << kminmerCounter._nbKminmers << endl;
			Logger::get().debug() << "\tNb kminmers: " << kminmerCounter._nbSolidKminmers;



			_kminmerFile.close();
			_kminmerAbundanceFile.close();
		}
		else{

			auto start = high_resolution_clock::now();

			/*
			if(_kminmerSize == _kminmerSizeFirst+2){ //First multiplex pass
				Logger::get().debug() << "Loading refined abundances";
				auto start = high_resolution_clock::now();
				loadRefinedAbundances();
				Logger::get().debug() << "\tDuration: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";
	

				rewriteMinimizerReadForMultiplex();

				KminmerAbundanceMap().swap(_kminmerAbundances);
			}
			*/

			//cout << "single core here" << endl;
			KminmerParserParallel parser2(_outputDir + "/read_data_corrected.txt", 0, _kminmerSize, false, false, _nbCores);
			parser2.parse(IndexKminmerFunctor(*this, false));


			//_parsingContigs = true;
			_nextReadIndexWriter = 0;
			
			//cout << "single core here" << endl;
			KminmerParserParallel parser3(_outputDir + "/unitig_data.txt", 0, _kminmerSize, false, false, _nbCores);
			parser3.parse(IndexKminmerFunctor(*this, true));
			
			Logger::get().debug() << "\tDuration: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";
	
			Logger::get().debug() << "\tNb kminmers: " << _mdbgNodesLight.size();
			
			Logger::get().debug() << "Dump kminmer abundances";
			start = high_resolution_clock::now();
			_kminmerAbundanceFile.open(_outputDir + "/kminmerData_abundance.txt");
			for(const auto& it : _mdbgNodesLight){

				u_int128_t vecHash = it.first;

				//u_int32_t nodeName = it.second._index;
				u_int32_t abundance = it.second;

				MDBG::writeKminmerAbundance(vecHash, abundance, _kminmerAbundanceFile);
			}
			Logger::get().debug() << "\tDuration: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";
			_kminmerAbundanceFile.close();
			
	
			//getchar();
		}







	}


	//Logger::get().debug() << "Duration: " << duration_cast<seconds>(high_resolution_clock::now() - start).count();
	//cout << _mdbg->_dbg_nodes.size() << endl;
	
	
	//if(!_isFirstPass){
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

	//u_int64_t nbKminmers = _mdbgNodesLight.size();
	//Logger::get().debug() << "Nb solid kminmers: " << _mdbgNodesLight.size();
	
	_fileSmallContigs.close();

	Logger::get().debug() << "Nb small contigs written:\t" << _kminmerSize << "\t" << _nbSmallContigs;

	_kminmerFile.close();
	_kminmerAbundanceFile.close();



	if(_isFirstPass){
		fs::copy(_outputDir + "/kminmerData_abundance.txt", _outputDir + "/kminmerData_abundance_init.txt", fs::copy_options::overwrite_existing);
	}


	if(_kminmerSize == _kminmerSizeFirst+1){
		fs::copy(_outputDir + "/kminmerData_abundance.txt", _outputDir + "/kminmerData_abundance_init_k" + to_string(_kminmerSizeFirst+1) + ".txt", fs::copy_options::overwrite_existing);
	}

	KminmerAbundanceMap().swap(_kminmerAbundances); //_kminmerAbundances.clear();  
	

	if(_kminmerSize > _kminmerSizeFirst+1){
		//cout << "Compute k+1 unitig graph" << endl;
		computeNextUnitigGraph();
		return;
	}

	//if(_kminmerSize == _kminmerSizeLast){
	//	removeDuplicatedSmallContigs();
	//}

	//u_int32_t nbNodes = _mdbgNodesLight.size();
	//
	//dumpNodes();
	//cout << "MDBG nb nodes: " << _mdbgNodesLight.size() << endl;
	//Logger::get().debug() << "Duration: " << duration_cast<seconds>(high_resolution_clock::now() - start).count();
	KminmerAbundanceMap().swap(_mdbgNodesLight); //_mdbgNodesLight.clear();
	
	//_mdbg->dump(_outputDir + "/kminmerData_min.txt");

	_kminmerFile.close();
	_kminmerAbundanceFile.close();
	

	//cout << "A RMETTRE" << endl;
	//delete _mdbg;
	
	createGfa();


	//if(_isFirstPass){
		//const auto copyOptions = fs::copy_options::overwrite_existing;
		//fs::copy(_outputDir + "/kminmerData_min.txt", _outputDir + "/kminmerData_min_init.txt", copyOptions);
		//_mdbg->dump(_outputDir + "/mdbg_nodes_init.gz");

		/*
		_mdbg = new MDBG(_kminmerSize);
		_mdbg->load(_outputDir + "/kminmerData_min_init.txt", false);

		KmerVec vec;
		vec._kmers = {7892019983131394, 62446187524336769, 47730431539434127, 90993437045220547};
		cout << _mdbg->_dbg_nodes[vec]._index << endl;
		getchar();
		*/

	//}

	
	Logger::get().debug() << "Checksum unitig nodes:   " << _checksum_unitigNodes;
	Logger::get().debug() << "Checksum unitig edges:   " << _checksum_unitigEdges;
	Logger::get().debug() << "Checksum unitig abundance:   " << _checksum_unitigAbundances;

	//exit(1);
	//if(_isFirstPass){
	//	ofstream file_data(_outputDir + "/data.txt");

	//	file_data.write((const char*)&_nbReads, sizeof(_nbReads));

	//	file_data.close();
	//}

	//_file_noKminmerReads.close();

	//closeLogFile();
	//exit(1);

	//if(_kminmerSize == _kminmerSizeFirst+1){
	//	cout << "kill k5" << endl;
	//	exit(1);
	//}


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
	MdbgEdgeMap3().swap(_mdbgEdges3); //_mdbgEdges3.clear();
	//MdbgEdgeMap6().swap(_mdbgEdges5); //_mdbgEdges5.clear();
	_mdbgEdges10.clear();


	Logger::get().debug() << "Deterministic unitigs";
	start = high_resolution_clock::now();
	computeDeterministicUnitigs();
	Logger::get().debug() << "\tDuration: " << duration_cast<seconds>(high_resolution_clock::now() - start).count();


	Logger::get().debug() << "Index unitig edges";
	start = high_resolution_clock::now();
	indexUnitigEdges();
	Logger::get().debug() << "\tDuration: " << duration_cast<seconds>(high_resolution_clock::now() - start).count();



	Logger::get().debug() << "Compute unitig edges";
	start = high_resolution_clock::now();

	_unitigGraphFile_edges_successors = ofstream(_outputDir + "/unitigGraph.edges.successors.bin");

	computeUnitigEdges();

	Logger::get().debug() << "\tDuration: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s " << (peakrss() / 1024.0 / 1024.0 / 1024.0);
	_unitigGraphFile_edges_successors.close();
	//UnitigEdgeMap().swap(_startNode_to_unitigIndex); //_startNode_to_unitigIndex.clear();
	//UnitigEdgeMap().swap(_endNode_to_unitigIndex); //_endNode_to_unitigIndex.clear();
	MdbgEdgeMap4().swap(_mdbgEdges4); //_mdbgEdges4.clear();
	_unitigEdgeMap.clear();
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

	writeUnitigGraphStat(_nbUnitigNodes, _nbUnitigEdges);
	//Logger::get().debug() << "Nb edges:   " << _nbEdges;


	
	//_mdbgNoFilter->dump(_outputDir + "/mdbg_nodes_noFilter.gz");
	//_mdbg->dump(_outputDir + "/mdbg_nodes.gz");
	//_mdbg->_dbg_nodes.clear();
	//_mdbg->_dbg_edges.clear();



}

void CreateMdbg::writeUnitigGraphStat(const u_int64_t nbUnitigNodes, const u_int64_t nbUnitigEdges){
	
	ofstream unitigGraphFile_stats(_outputDir + "/unitigGraph.stats.bin");
	unitigGraphFile_stats.write((const char*)&nbUnitigNodes, sizeof(nbUnitigNodes));
	unitigGraphFile_stats.write((const char*)&nbUnitigEdges, sizeof(nbUnitigEdges));
	unitigGraphFile_stats.close();

}

struct DeterministicUnitigNode{
	u_int128_t _hash;
	vector<MinimizerType> _minimizers;
};

void CreateMdbg::computeDeterministicUnitigs(){

	vector<DeterministicUnitigNode> unitigs;
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
		

		KmerVec vec;
		vec._kmers = minimizers;
		vec = vec.normalize();

		u_int128_t vecHash = vec.hash128();

		unitigs.push_back({vecHash, vec._kmers});
	}

	nodeFile.close();

	_nbUnitigNodes = unitigs.size();
	Logger::get().debug() << "Nb unitigs: " << unitigs.size();
	
	std::sort(unitigs.begin(), unitigs.end(), [](const DeterministicUnitigNode & a, const DeterministicUnitigNode & b){
		return a._hash < b._hash;
	});

	u_int32_t unitigIndex = 0;
	ofstream outputFile(_outputDir + "/unitigGraph.nodes.bin");

	for(const DeterministicUnitigNode& unitig : unitigs){
		
		dumpUnitigNode(unitigIndex, unitig._minimizers, outputFile);
		unitigIndex += 2;
	}

	outputFile.close();
}

/*

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
*/

void CreateMdbg::indexEdges(){

	Logger::get().debug() << "\tDereplicating edges";
	auto start = high_resolution_clock::now();
	EdgeIndexer edgeIndexer(*this);
	edgeIndexer.execute();
	Logger::get().debug() << "\t\tDone: " << edgeIndexer._nbEdges << " " << edgeIndexer._checksum << " " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";

	const u_int64_t nbEdges = edgeIndexer._nbEdges;
	const string edgeFilename = edgeIndexer.getOutputFilename();

	Logger::get().debug() << "\tBuilding mphf: " << nbEdges;
	start = high_resolution_clock::now();
	boomphf::file_binary<u_int128_t> inputFile(edgeFilename.c_str());

	double gammaFactor = 1.0; // lowest bit/elem is achieved with gamma=1, higher values lead to larger mphf but faster construction/query
	_mdbgEdges10._keys = new boomphf::mphf<u_int128_t, hasher_t>(nbEdges, inputFile, _nbCores, gammaFactor, true, false);
	Logger::get().debug() << "\t\tDone: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";
	
	Logger::get().debug() << "\tCreating edge values";
	start = high_resolution_clock::now();
	_mdbgEdges10._values.resize(nbEdges);
	Logger::get().debug() << "\t\tDone: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";
	


	Logger::get().debug() << "\tIndexing edges";
	start = high_resolution_clock::now();

	//_mdbgEdges5.setup();
	size_t i = 0;
	//cout << "indexEdges single core" << endl;
	//cout << "About to index edges" << endl;
	//getchar();
	
	//_edgeIndexElements = 0;

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
				
				i+= 1;
				//if(i % 1000000 == 0) cout << "\t" << i << " " << _mdbgEdges3.size() << "    " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << endl;
			}

			if(isEOF) break;
			
			indexEdge(vec, 0);
			//indexEdge(vec.reverse(), nodeName);

		}
	}

	kminmerFile.close();
	fs::remove(edgeFilename);

	Logger::get().debug() << "\t\tDone: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";
	
	//cout << "done" << endl;
	//Logger::get().debug() << "Edge index elements: " << _edgeIndexElements;
	//cout << _mdbgEdges.size() << endl;
	//getchar();
	//_mdbgEdges5.clear();


}

void CreateMdbg::indexEdge(const KmerVec& vec, u_int32_t nodeName){

	//u_int64_t& lala = _edgeIndexElements;
	bool needprint = false;

	bool isReversedPrefix;
	bool isReversedSuffix;
	KmerVec prefix = vec.prefix().normalize(isReversedPrefix);
	KmerVec suffix = vec.suffix().normalize(isReversedSuffix);
	MinimizerType prefix_minimizer = vec._kmers[vec._kmers.size()-1];
	MinimizerType suffix_minimizer = vec._kmers[0];
	

	const u_int128_t suffixHash = suffix.hash128();
	const u_int128_t suffixKey = _mdbgEdges10._keys->lookup(suffixHash);
	int partition_suffix = suffixHash % _mutexes.size();
	






	omp_set_lock(&_mutexes[partition_suffix]);

	bool exists1 = successorExists(suffix, suffixHash, isReversedSuffix, false);

	if(!exists1){
		
		KminmerEdge33& edge = _mdbgEdges10._values[suffixKey];

		if(edge._minimizer1 == -1){

			//KminmerEdge33 edge;

			edge._minimizer1 = suffix_minimizer;
			edge._isReversed1 = isReversedSuffix;
			edge._isPrefix1 = false;
			edge._isUnitigged1 = false;
			edge._hasMultipleSuccessors1 = false;

			//_mdbgEdges5[suffixHash] = edge;

			//if(((u_int64_t) suffixHash) == 564419951353098858){
			//	cout << "inserting a: " << suffix_minimizer << endl;
			//}

		}
		else{
			//KminmerEdge33& edge = _mdbgEdges5[suffixHash];

			edge._minimizer2 = suffix_minimizer;
			edge._isReversed2 = isReversedSuffix;
			edge._isPrefix2 = false;
			edge._isUnitigged2 = false;
			edge._hasMultipleSuccessors2 = false;

			//if(((u_int64_t) suffixHash) == 564419951353098858){
			//	cout << "inserting b: " << suffix_minimizer << endl;
			//}

		}
		//_mdbgEdgesTest2[suffixHash] = {};
		//_mdbgEdges5.insert(suffix, suffixHash, {suffix_minimizer, isReversedSuffix, false, false, false});
		//_mdbgEdges5[suffixHash].push_back({suffix_minimizer, isReversedSuffix, false, false, false});
	}

	omp_unset_lock(&_mutexes[partition_suffix]);









	const u_int128_t prefixHash = prefix.hash128();
	const u_int128_t prefixKey = _mdbgEdges10._keys->lookup(prefixHash);
	int partition_prefix = prefixHash % _mutexes.size();


	omp_set_lock(&_mutexes[partition_prefix]);

	bool exists2 = successorExists(prefix, prefixHash, isReversedPrefix, true);

	if(!exists2){
		//_mdbgEdgesTest2[prefixHash] = {};
		//_mdbgEdges5[prefixHash].push_back({prefix_minimizer, isReversedPrefix, true, false, false});
		//_mdbgEdges5.insert(prefix, prefixHash, {prefix_minimizer, isReversedPrefix, true, false, false});

		KminmerEdge33& edge = _mdbgEdges10._values[prefixKey];


		//if(_mdbgEdges5.find(prefixHash) == _mdbgEdges5.end()){
		if(edge._minimizer1 == -1){

			//KminmerEdge33 edge;

			edge._minimizer1 = prefix_minimizer;
			edge._isReversed1 = isReversedPrefix;
			edge._isPrefix1 = true;
			edge._isUnitigged1 = false;
			edge._hasMultipleSuccessors1 = false;

			//_mdbgEdges5[prefixHash] = edge;

			//if(((u_int64_t) prefixHash) == 564419951353098858){
			//	cout << "inserting c: " << prefix_minimizer << endl;
			//}

		}
		else{
			//KminmerEdge33& edge = _mdbgEdges5[prefixHash];

			edge._minimizer2 = prefix_minimizer;
			edge._isReversed2 = isReversedPrefix;
			edge._isPrefix2 = true;
			edge._isUnitigged2 = false;
			edge._hasMultipleSuccessors2 = false;

			//if(((u_int64_t) prefixHash) == 564419951353098858){
			//	cout << "inserting d: " << prefix_minimizer << endl;
			//}

		}

	}
	
	omp_unset_lock(&_mutexes[partition_prefix]);

	/*
	#pragma omp critical
	{
		
		
		//const vector<MinimizerType> ms = {33824, 75773092, 33759410, 42640641, 451446017};
		//bool print = (vec._kmers == ms);

		//if(print){
		//	cout << vec.toString() << endl;
		//	cout << "Add prefix: " << prefix.toString() << "    " << prefix_minimizer << endl;
		//	cout << "Add suffix: " << suffix.toString() << "    " << suffix_minimizer << endl;
		//}

		_mdbgEdges3[suffixHash].push_back({suffix_minimizer, isReversedSuffix, false, false, false});
		_mdbgEdges3[prefixHash].push_back({prefix_minimizer, isReversedPrefix, true, false, false});
	}
	*/
}


bool CreateMdbg::successorExists(const KmerVec& edgeVec, const u_int128_t& edgeVecHash, const bool isReversed, const bool isPrefix){


	const u_int128_t key = _mdbgEdges10._keys->lookup(edgeVecHash);
	KminmerEdge33& edge = _mdbgEdges10._values[key];

	if(edge._minimizer1 == -1) return false;

	//return false;
	//if(_mdbgEdges5.find(key) == _mdbgEdges5.end()) return false;

	bool isPalindrome = edgeVec.isPalindrome();

	//KminmerEdge33& edge = _mdbgEdges5[key];


	if(edge._minimizer1 != -1){
		
		//if(((u_int64_t) key) == 564419951353098858){
		//	cout << "\tblurp: " << edge._minimizer1 << endl;
		//}


		if(isPalindrome){
			edge._hasMultipleSuccessors1 = true;
			return true;
		}
		
		if(edge._isReversed1 == isReversed && edge._isPrefix1 == isPrefix){
			edge._hasMultipleSuccessors1 = true;
			return true;
		}

		if(edge._isReversed1 == !isReversed && edge._isPrefix1 == !isPrefix){
			edge._hasMultipleSuccessors1 = true;
			return true;
		}

	}

	if(edge._minimizer2 != -1){
		
		//if(((u_int64_t) key) == 564419951353098858){
		//	cout << "\tblurp:    " << edge._minimizer2 << " " << isReversed << "    " << isPrefix << " " << edge._isReversed2 << " " << edge._isPrefix2 << endl;
		//}


		if(isPalindrome){
			edge._hasMultipleSuccessors2 = true;
			return true;
		}

		if(edge._isReversed2 == isReversed && edge._isPrefix2 == isPrefix){
			edge._hasMultipleSuccessors2 = true;
			return true;
		}

		if(edge._isReversed2 == !isReversed && edge._isPrefix2 == !isPrefix){
			edge._hasMultipleSuccessors2 = true;
			return true;
		}
	}

	/*
	//for(KminmerEdge3& edge : _mdbgEdges3[key]){
	for(KminmerEdge3& edge : vecs){
		
		if(isPalindrome){
			edge._hasMultipleSuccessors = true;
			return true;
		}

		if(edge._isReversed == isReversed && edge._isPrefix == isPrefix){
			edge._hasMultipleSuccessors = true;
			return true;
		}

		if(edge._isReversed == !isReversed && edge._isPrefix == !isPrefix){
			edge._hasMultipleSuccessors = true;
			return true;
		}

	}
	*/

	//if(((u_int64_t) key) == 564419951353098858){
	//	cout << "blurp" << endl;
	//}

	return false;


}


void CreateMdbg::computeUnitigNodes(){

	//_debugFile = ofstream(_outputDir + "/debug.txt");

	ComputeUnitigFunctor computeUnitigFunctor(*this);


	//cout << "computeUnitigNodes: besoin d'activer le if(_graph->_isNodeIndexIndexed[nodeIndex]) exist = true (requiert d'indexer tous les kmervec des unitigs)" << endl;
	//cout << "Parallelisation unitig creation: bien faire en sorte que les nodeSource soit randomiser, surtout ne pas paralleliser a partir des node d'un read successif qui vont probablement etre dans le meme unitig" << endl;
	//cout << "nb cores 1" << endl;
	


	_unitigIndex = 0;

	ifstream kminmerFile(_outputDir + "/kminmerData_min.txt");

	//cout << "ComputeUnitig nb cores 1 " << endl;

	//cout << "single core here" << endl;
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

	//_debugFile.close();
	//exit(1);
}




void CreateMdbg::getSuccessors(const KmerVec& sourceVec, vector<KmerVec>& successors){

	successors.clear();
	//KmerVec vec = keys[i];
	//KmerVec vec_rev = vec.reverse();

	//KmerVec vec_rev_suffix = vec_rev.suffix();
	KmerVec vec_suffix = sourceVec.suffix();
	
	KmerVec suffix = vec_suffix.normalize();

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

			if (vec_suffix == v.prefix()) {
				successors.push_back(v);
			}

		}
		else{
			v._kmers.insert(v._kmers.begin(), edge._minimizer);

			KmerVec v_rev = v.reverse();

			if (vec_suffix == v_rev.prefix()) {
				successors.push_back(v_rev);
				//dumpEdge(id, edgePlus, id2, edgeMinus);
			}

		}
		//v = v.normalize();
		//v._kmers.insert(v._kmers.begin(), edge._minimizer);

		/*
		if(_mdbg->_dbg_nodes.find(v) == _mdbg->_dbg_nodes.end()){
			cout << "suffix" << endl;
			getchar();
		}
		*/
		

		//if(v==vec) continue;

		//KmerVec v_rev = v.reverse();
		
		//if(v_rev==vec) continue;

		//u_int32_t id2 = edge._nodeName;

		//if (vec_suffix == v.prefix()) {
		//	successors.push_back(v);

			//dumpEdge(id, edgePlus, id2, edgePlus);
		//}
		//if (vec_suffix == v_rev.prefix()) {
		//	successors.push_back(v_rev);
			//dumpEdge(id, edgePlus, id2, edgeMinus);
		//}
		//if (vec_rev_suffix == v.prefix()) {
		//	successors.push_back(v_rev);
		//	//dumpEdge(id, edgeMinus, id2, edgePlus);
		//}
		//if (vec_rev_suffix == v_rev.prefix()) {
		//	successors.push_back(v);
		//	//dumpEdge(id, edgeMinus, id2, edgeMinus);
		//}
	}

}

void CreateMdbg::getPredecessors(const KmerVec& vec, vector<KmerVec>& predecessors){


	predecessors.clear();
	//KmerVec vec = keys[i];
	KmerVec vec_rev = vec.reverse();

	KmerVec vec_rev_suffix = vec_rev.suffix();
	//KmerVec vec_suffix = vec.suffix();
	
	const KmerVec prefix = vec.prefix().normalize();

	
	//cout << "\tRofl preds: " << _mdbgEdges3[prefix.hash128()].size() << endl;

	for(KminmerEdge3& edge : _mdbgEdges3[prefix.hash128()]){

		//cout << "\tPred: " << edge._isReversed << " " << edge._minimizer << endl;

		KmerVec v = prefix;
		if(edge._isReversed){
			v = v.reverse();
		}
		
		
		if(edge._isPrefix){

			if (vec_rev_suffix == v) {
				
				v._kmers.push_back(edge._minimizer);
				KmerVec v_rev = v.reverse();
				//cout << "rofl 3 " << v_rev.toString() << endl;
				predecessors.push_back(v_rev);
				//dumpEdge(id, edgeMinus, id2, edgePlus);
			}

		}
		else{
			KmerVec v_rev = v.reverse();

			if (vec_rev_suffix == v_rev) {
				
				v._kmers.insert(v._kmers.begin(), edge._minimizer);
				//cout << "rofl 4 " << v.toString() << endl;
				predecessors.push_back(v);
				//dumpEdge(id, edgeMinus, id2, edgeMinus);
				//predecessors.push_back(v_rev);
			}
			
		}
		

		//if(v==vec) continue;

		
		//if(v_rev==vec) continue;

		//cout << "\tALlo? " << v.toString() << endl;
		//cout << "\tALlo? " << v_rev.toString() << endl;

		//u_int32_t id2 = edge._nodeName;

		//if (vec_suffix == v.prefix()) {
		//	//cout << "rofl 1 " << " " << v.toString() << endl;
		//	predecessors.push_back(v);
		//	//dumpEdge(id, edgePlus, id2, edgePlus);
		//}
		//if (vec_suffix == v_rev.prefix()) {
		//	//cout << "rofl 2 " << v_rev.toString() << endl;
		//	predecessors.push_back(v_rev);
		//}



	}


}





bool CreateMdbg::hasSingleSuccessor(const KmerVec& vec){


	int nbSuccessors = 0;

	//KmerVec vec_rev = vec.reverse();
	//KmerVec vec_rev_suffix = vec_rev.suffix();
	KmerVec vec_suffix = vec.suffix();

	KmerVec suffix = vec_suffix.normalize();
	
	//bool print = _mdbgEdges3[suffix.hash128()].size() > 4;

	//if(print) cout << "----" << endl;

	//cout << (u_int64_t) suffix.hash128() << endl;

	for(KminmerEdge3& edge : _mdbgEdges3[suffix.hash128()]){
		
		//cout << "\tedge" << endl;

		KmerVec v = suffix;

		if(edge._isReversed){
			//cout << "\t\tIs reversed" << endl;
			v = v.reverse();
		}

		if(edge._isPrefix){
			//cout << "\t\tIs prefix" << endl;
			//v._kmers.push_back(edge._minimizer);

			
			if (vec_suffix == v) {
				//cout << "\ta: " << edge._minimizer << " " << edge._isReversed << " " <<  edge._isPrefix << " " << endl; //edge._isPrefix toujours 1
				//if(print) cout << "a: " << edge._isReversed << endl;
				//if(!edge._isReversed){
				nbSuccessors += 1;
			}

		}
		else{
			//cout << "\t\tNo prefix" << endl;
			//v._kmers.insert(v._kmers.begin(), edge._minimizer);


			KmerVec v_rev = v.reverse();

			//if (vec_suffix == v_rev.prefix()) {
			if (vec_suffix == v_rev) {
				//if(print) cout << "b: " << edge._isReversed << endl;
				//cout << "\tb: " << edge._minimizer << " " << edge._isReversed << " " <<  edge._isPrefix << " " << endl; //edge._isPrefix toujours 0
				//if(edge._isReversed){
				nbSuccessors += 1;
			}
			
		}
		

		/*
		if(v==vec){
			//cout << "lul1: " << edge._isReversed << " " << edge._isPrefix << endl; //edge._isPrefix toujours 0

			if (vec_suffix == v.prefix()) {
				cout << "derp a" << endl;
				getchar();
			}
			if (vec_suffix == v_rev.prefix()) {
				cout << vec.toString() << endl;
				cout << "derp b" << endl;
				getchar();
			}

			continue;
		}
		
		if(v_rev==vec){
			//cout << "lul2: " << (vec_suffix == v.prefix()) << " " << (vec_suffix == v_rev.prefix()) << endl;

			if (vec_suffix == v.prefix()) {
				cout << "derp c" << endl;
				getchar();
			}
			if (vec_suffix == v_rev.prefix()) {
				cout << "derp d" << endl;
				getchar();
			}
			
			//cout << "lul2: " << edge._isReversed << " " << edge._isPrefix << endl; //edge._isPrefix toujours 1
			continue;
		}
		*/

		//if(edge._isPrefix){
		//	nbSuccessorsNew += 1;
		//}

		//if (vec_rev_suffix == v.prefix()) {
		//	//cout << "c: " << edge._isReversed << " " <<  edge._isPrefix << " " << endl;
		//	nbSuccessors += 1;
		//	cout << "derp" << endl;
		//	getchar();
		//}
		//if (vec_rev_suffix == v_rev.prefix()) {
		//	//cout << "d: " << edge._isReversed << " " <<  edge._isPrefix << " " << endl;
		//	nbSuccessors += 1;
		//	cout << "derp" << endl;
		//	getchar();
		//}

		//if(nbSuccessors > 1) return false;
	}

	//cout << "\tNb succs: " << nbSuccessors << endl;
	//cout << nbSuccessors << " " << _mdbgEdges3[suffix.hash128()].size() << endl;

	//if(nbSuccessors != _mdbgEdges3[suffix.hash128()].size()){
	//	getchar();
	//}

	return nbSuccessors == 1;

}






bool CreateMdbg::getNbSuccessors(const KmerVec& vec, KmerVec& singleSuccessor){

	bool print = false;
	bool hasMultiSuccessors = false;

	int nbSuccessors = 0;

	//KmerVec vec_rev = vec.reverse();
	//KmerVec vec_rev_suffix = vec_rev.suffix();
	KmerVec vec_suffix = vec.suffix();

	KmerVec suffix = vec_suffix.normalize();
	
	//if(suffix.isPalindrome()){
		//print = true;
	//}
	//bool print = _mdbgEdges3[suffix.hash128()].size() > 4;

	if(print){
		cout << endl;
		cout << "----" << endl;
		cout << suffix.toString() << endl;
	}

	/*
	KminmerEdge33& edge = _mdbgEdges5[suffix.hash128()];

	if(edge._minimizer1 != -1){
		KmerVec v = suffix;

		if(edge._isReversed1){
			//cout << "\t\tIs reversed" << endl;
			v = v.reverse();
		}

		if(edge._isPrefix1){
			//cout << "\t\tIs prefix" << endl;
			//v._kmers.push_back(edge._minimizer);

			
			if (vec_suffix == v) {
				if(edge._hasMultipleSuccessors1){
					hasMultiSuccessors = true;
					return false;
				}
				
				v._kmers.push_back(edge._minimizer1);

				singleSuccessor = v;
				//if(print) cout << "\ta: " << edge._minimizer1 << " " <<  edge._isReversed1 << " " <<  edge._isPrefix1 << " " << edge._hasMultipleSuccessors << endl; //edge._isPrefix toujours 1
				
				nbSuccessors += 1;
			}

		}
		else{
			//cout << "\t\tNo prefix" << endl;


			KmerVec v_rev = v.reverse();

			//if (vec_suffix == v_rev.prefix()) {
			if (vec_suffix == v_rev) {
				
				if(edge._hasMultipleSuccessors1){
					hasMultiSuccessors = true;
					return false;
				}
				
				v._kmers.insert(v._kmers.begin(), edge._minimizer1);
				KmerVec v_rev = v.reverse();

				singleSuccessor = v_rev;
				//if(print) cout << "\tb: " << edge._minimizer1 << " " << edge._isReversed1 << " " <<  edge._isPrefix1 << " " << edge._hasMultipleSuccessors1 << endl; //edge._isPrefix toujours 0
				
				nbSuccessors += 1;
			}
			
		}
		
	}

	if(edge._minimizer2 != -1){

		KmerVec v = suffix;

		if(edge._isReversed2){
			//cout << "\t\tIs reversed" << endl;
			v = v.reverse();
		}

		if(edge._isPrefix2){
			//cout << "\t\tIs prefix" << endl;

			
			if (vec_suffix == v) {
				if(edge._hasMultipleSuccessors2){
					hasMultiSuccessors = true;
					return false;
				}
				
				v._kmers.push_back(edge._minimizer2);

				singleSuccessor = v;
				//if(print) cout << "\ta: " << edge._minimizer2 << " " <<  edge._isReversed2 << " " <<  edge._isPrefix2 << " " << edge._hasMultipleSuccessors1 << endl; //edge._isPrefix toujours 1
				
				nbSuccessors += 1;
			}

		}
		else{
			//cout << "\t\tNo prefix" << endl;


			KmerVec v_rev = v.reverse();

			//if (vec_suffix == v_rev.prefix()) {
			if (vec_suffix == v_rev) {
				
				if(edge._hasMultipleSuccessors2){
					hasMultiSuccessors = true;
					return false;
				}
				
				v._kmers.insert(v._kmers.begin(), edge._minimizer2);
				KmerVec v_rev = v.reverse();

				singleSuccessor = v_rev;
				//if(print) cout << "\tb: " << edge._minimizer2 << " " << edge._isReversed2 << " " <<  edge._isPrefix2 << " " << edge._hasMultipleSuccessors2 << endl; //edge._isPrefix toujours 0
				
				nbSuccessors += 1;
			}
			
		}
	}
	*/

	const u_int128_t key = _mdbgEdges10._keys->lookup(suffix.hash128());
	const KminmerEdge33& e = _mdbgEdges10._values[key];



	//const KminmerEdge33& e = _mdbgEdges5[suffix.hash128()];
	vector<KminmerEdge3> edges;
	if(e._minimizer1 != -1){
		//cout << (u_int64_t) suffix.hash128() << endl;
		//cout << "\t" << e._minimizer1 << endl;
		edges.push_back({e._minimizer1, e._isReversed1, e._isPrefix1, e._isUnitigged1, e._hasMultipleSuccessors1});
	}

	if(e._minimizer2 != -1){
		//cout << "\t" << e._minimizer2 << endl;
		edges.push_back({e._minimizer2, e._isReversed2, e._isPrefix2, e._isUnitigged2, e._hasMultipleSuccessors2});
	}

	//cout << edges.size() << endl;

	//for(KminmerEdge3& edge : _mdbgEdges3[suffix.hash128()]){
	for(KminmerEdge3& edge : edges){
		
		//cout << "\tedge" << endl;

		KmerVec v = suffix;

		if(edge._isReversed){
			//cout << "\t\tIs reversed" << endl;
			v = v.reverse();
		}


		if(edge._isPrefix){
			//cout << "\t\tIs prefix" << endl;
			//v._kmers.push_back(edge._minimizer);

			
			if (vec_suffix == v) {
				if(edge._hasMultipleSuccessors){
					hasMultiSuccessors = true;
					//return false;
				}
				
				v._kmers.push_back(edge._minimizer);

				singleSuccessor = v;
				//cout << "\taa: " << edge._minimizer << " " <<  edge._isReversed << " " <<  edge._isPrefix << " " << edge._hasMultipleSuccessors << endl; //edge._isPrefix toujours 1
				//if(!edge._isReversed){
				nbSuccessors += 1;
			}

		}
		else{
			//cout << "\t\tNo prefix" << endl;
			//v._kmers.insert(v._kmers.begin(), edge._minimizer);


			KmerVec v_rev = v.reverse();

			//if (vec_suffix == v_rev.prefix()) {
			if (vec_suffix == v_rev) {
				
				if(edge._hasMultipleSuccessors){
					hasMultiSuccessors = true;
					//return false;
				}
				
				v._kmers.insert(v._kmers.begin(), edge._minimizer);
				KmerVec v_rev = v.reverse();

				singleSuccessor = v_rev;
				//cout << "\tbb: " << edge._minimizer << " " << edge._isReversed << " " <<  edge._isPrefix << " " << edge._hasMultipleSuccessors << endl; //edge._isPrefix toujours 0
				//if(edge._isReversed){
				nbSuccessors += 1;
			}
			
		}


		
	}
	

	//vector<MinimizerType> ms = {644249177, 644245349, 644249177};
	//bool print2 = (suffix._kmers == ms);

	//if(print2) getchar();

	if(print){
		cout << "\tNb successors: " << nbSuccessors  << " " << hasMultiSuccessors << endl;
		for(KminmerEdge3& edge : _mdbgEdges3[suffix.hash128()]){
			cout << "\t\t" << edge._minimizer << " " << edge._isReversed << " " << edge._isPrefix << endl;
		}
	}

	if(nbSuccessors == 1 && hasMultiSuccessors){
		//cout << "derp: " << suffix.toString() << endl;
		//getchar();
	}

	if(nbSuccessors > 1 && !hasMultiSuccessors){
		//cout << "derp: " << suffix.toString() << endl;
		//getchar();
	}

	//return hasMultiSuccessors;
	if(hasMultiSuccessors) return false;
	if(nbSuccessors == 0) return false;

	return true;
}


bool CreateMdbg::getNbPredecessors(const KmerVec& vec, KmerVec& singlePredecessor){

	bool print = false;
	//cout << "huhu: " << vec.toString() << endl;
	bool hasMultiPredecessors = false;

	int nbPredecessors = 0;

	KmerVec vec_rev = vec.reverse();

	
	const KmerVec prefix = vec.prefix().normalize();

	KmerVec vec_rev_suffix = vec_rev.suffix();
	//KmerVec vec_suffix = vec.suffix();
	
	//if(prefix.isPalindrome()){
		//print = true;
	//}
	//bool print = _mdbgEdges3[suffix.hash128()].size() > 4;

	if(print){
		cout << endl;
		cout << "----" << endl;
		cout << prefix.toString() << endl;
	}

	//cout << endl;
	//cout << "----" << endl;
	//cout << prefix.toString() << endl;



	/*
	if(edge._minimizer1 != -1){

		KmerVec v = prefix;
		if(edge._isReversed1){
			v = v.reverse();
		}
		
		
		if(edge._isPrefix1){

			if (vec_rev_suffix == v) {

				if(edge._hasMultipleSuccessors1){
					hasMultiPredecessors = true;
					return false;
				}

				v._kmers.push_back(edge._minimizer1);
				KmerVec v_rev = v.reverse();
				singlePredecessor = v_rev;

				if(print) cout << "\tcc: " << edge._isReversed1 << " " <<  edge._isPrefix1 << " " << edge._hasMultipleSuccessors1 << endl; //edge._isPrefix toujours 1
				nbPredecessors += 1;
			}

		}
		else{

			KmerVec v_rev = v.reverse();

			if (vec_rev_suffix == v_rev) {
				
				if(edge._hasMultipleSuccessors1){
					hasMultiPredecessors = true;
					return false;
				}

				v._kmers.insert(v._kmers.begin(), edge._minimizer1);
				singlePredecessor = v;

				if(print) cout << "\tdd: " << edge._isReversed1 << " " <<  edge._isPrefix1 << " " << edge._hasMultipleSuccessors1 << endl; //edge._isPrefix toujours 0
				nbPredecessors += 1;
			}



		}
	}

	if(edge._minimizer2 != -1){

		KmerVec v = prefix;
		if(edge._isReversed2){
			v = v.reverse();
		}
		
		
		if(edge._isPrefix2){

			if (vec_rev_suffix == v) {

				if(edge._hasMultipleSuccessors2){
					hasMultiPredecessors = true;
					return false;
				}

				v._kmers.push_back(edge._minimizer2);
				KmerVec v_rev = v.reverse();
				singlePredecessor = v_rev;

				if(print) cout << "\tcc: " << edge._isReversed2 << " " <<  edge._isPrefix2 << " " << edge._hasMultipleSuccessors2 << endl; //edge._isPrefix toujours 1
				nbPredecessors += 1;
			}

		}
		else{

			KmerVec v_rev = v.reverse();

			if (vec_rev_suffix == v_rev) {
				
				if(edge._hasMultipleSuccessors2){
					hasMultiPredecessors = true;
					return false;
				}

				v._kmers.insert(v._kmers.begin(), edge._minimizer2);
				singlePredecessor = v;

				if(print) cout << "\tdd: " << edge._isReversed2 << " " <<  edge._isPrefix2 << " " << edge._hasMultipleSuccessors2 << endl; //edge._isPrefix toujours 0
				nbPredecessors += 1;
			}



		}
	}
	*/
	

	const u_int128_t key = _mdbgEdges10._keys->lookup(prefix.hash128());
	const KminmerEdge33& e = _mdbgEdges10._values[key];

	//KminmerEdge33& e = _mdbgEdges5[prefix.hash128()];
	vector<KminmerEdge3> edges;
	if(e._minimizer1 != -1){
		edges.push_back({e._minimizer1, e._isReversed1, e._isPrefix1, e._isUnitigged1, e._hasMultipleSuccessors1});
	}

	if(e._minimizer2 != -1){
		edges.push_back({e._minimizer2, e._isReversed2, e._isPrefix2, e._isUnitigged2, e._hasMultipleSuccessors2});
	}

	//for(KminmerEdge3& edge : _mdbgEdges3[prefix.hash128()]){
	for(KminmerEdge3& edge : edges){

		KmerVec v = prefix;
		if(edge._isReversed){
			v = v.reverse();
		}
		
		
		if(edge._isPrefix){

			if (vec_rev_suffix == v) {

				if(edge._hasMultipleSuccessors){
					hasMultiPredecessors = true;
					return false;
				}

				v._kmers.push_back(edge._minimizer);
				KmerVec v_rev = v.reverse();
				singlePredecessor = v_rev;

				if(print) cout << "\tcc: " << edge._isReversed << " " <<  edge._isPrefix << " " << edge._hasMultipleSuccessors << endl; //edge._isPrefix toujours 1
				nbPredecessors += 1;
			}

		}
		else{

			KmerVec v_rev = v.reverse();

			if (vec_rev_suffix == v_rev) {
				
				if(edge._hasMultipleSuccessors){
					hasMultiPredecessors = true;
					return false;
				}

				v._kmers.insert(v._kmers.begin(), edge._minimizer);
				singlePredecessor = v;

				if(print) cout << "\tdd: " << edge._isReversed << " " <<  edge._isPrefix << " " << edge._hasMultipleSuccessors << endl; //edge._isPrefix toujours 0
				nbPredecessors += 1;
			}



		}
		
	}
	

	if(print){
		cout << "\tNb predecessors: " << nbPredecessors << " " << hasMultiPredecessors << endl;
		for(KminmerEdge3& edge : _mdbgEdges3[prefix.hash128()]){
			cout << "\t\t" << edge._minimizer << " " << edge._isReversed << " " << edge._isPrefix << endl;
		}
	}

	if(nbPredecessors == 1 && hasMultiPredecessors){
		//cout << "derp: " << prefix.toString() << endl;
		//getchar();
	}

	if(nbPredecessors > 1 && !hasMultiPredecessors){
		//cout << "derp: " << prefix.toString() << endl;
		//getchar();
	}

	//return nbPredecessors;

	if(hasMultiPredecessors) return false;
	if(nbPredecessors == 0) return false;

	return true;
}

bool CreateMdbg::hasSinglePredecessor(const KmerVec& vec){

	//cout << "hihi: " << vec.toString() << endl;

	int nbPredecessors = 0;

	KmerVec vec_rev = vec.reverse();

	
	const KmerVec prefix = vec.prefix().normalize();

	KmerVec vec_rev_suffix = vec_rev.suffix();
	//KmerVec vec_suffix = vec.suffix();
	
	for(KminmerEdge3& edge : _mdbgEdges3[prefix.hash128()]){

		KmerVec v = prefix;
		if(edge._isReversed){
			v = v.reverse();
		}
		
		
		if(edge._isPrefix){
			//v._kmers.push_back(edge._minimizer);

			if (vec_rev_suffix == v) {
				//cout << "\tc: " << edge._isReversed << " " <<  edge._isPrefix << " " << v.toString() << endl; //edge._isPrefix toujours 1
				nbPredecessors += 1;
			}


		}
		else{
			//v._kmers.insert(v._kmers.begin(), edge._minimizer);

			KmerVec v_rev = v.reverse();

			if (vec_rev_suffix == v_rev) {
				//cout << "\td: " << edge._isReversed << " " <<  edge._isPrefix << " " << v.toString() << endl; //edge._isPrefix toujours 0
				nbPredecessors += 1;
			}

		}
		

		//if(v==vec) continue;

		
		//if(v_rev==vec) continue;


		//if (vec_suffix == v.prefix()) {
		//	//cout << "a: " << edge._isReversed << " " <<  edge._isPrefix << " " << endl;
		//	nbPredecessors += 1;
		//	cout << "derp" << endl;
		//	getchar();
		//}
		//if (vec_suffix == v_rev.prefix()) {
		//	//cout << "b: " << edge._isReversed << " " <<  edge._isPrefix << " " << endl;
		//	nbPredecessors += 1;
		//	cout << "derp" << endl;
		//	getchar();
		//}



		//if(nbPredecessors > 1) return false;
	}

	//cout << "\tNb preds: " << nbPredecessors << endl;
	
	return nbPredecessors == 1;

}


void CreateMdbg::getSuccessors_unitig(const KmerVec& vec, const UnitigType& unitigIndexFrom, vector<UnitigType>& successors){

	//cout << "----" << endl;
	//cout << vec.toString() << endl;

	successors.clear();

	KmerVec vec_suffix = vec.suffix();
	KmerVec suffix = vec_suffix.normalize();

	KmerVec vec_rev = vec.reverse();
	


	const u_int128_t key = _unitigEdgeMap._keys->lookup(suffix.hash128());
	const vector<KminmerEdgeU>& edges = _unitigEdgeMap._values[key];


	//if(_mdbgEdges4.find(suffix.hash128()) == _mdbgEdges4.end()){
	//	cout << "derp succ" << endl;
	//	return;
	//}
		
	//cout << "\tRofl succ: " << _mdbgEdges3[suffix.hash128()].size() << endl;

	//vector<UnitigEdge> successors2;

	for(const KminmerEdgeU& edge : edges){
		

		//if(unitigIndexFrom == edge._unitigIndex) continue;
		//if(unitigIndexFrom == Utils::unitigIndexToReverseDirection(edge._unitigIndex)) continue;


		KmerVec v = suffix;

		if(edge._isReversed){
			//cout << "\t\tIs reversed" << endl;
			v = v.reverse();
		}


		if(edge._isPrefix){
			
			if (vec_suffix == v) {
				//if(unitigIndexFrom == edge._unitigIndex && (vec.suffix().isPalindrome() || vec.prefix().isPalindrome())) continue;
				if(unitigIndexFrom == Utils::unitigIndexToReverseDirection(edge._unitigIndex)) continue;
				//cout << "a: " << unitigIndexFrom << " -> " << edge._unitigIndex << endl;
				successors.push_back(edge._unitigIndex);
			}

		}
		else{
			
			KmerVec v_rev = v.reverse();

			//if (vec_suffix == v_rev.prefix()) {
			if (vec_suffix == v_rev) {
				//if(unitigIndexFrom == Utils::unitigIndexToReverseDirection(edge._unitigIndex) && (vec.suffix().isPalindrome() || vec.prefix().isPalindrome())) continue;
				if(unitigIndexFrom == edge._unitigIndex) continue;
				//cout << "b: " << unitigIndexFrom << " -> " << Utils::unitigIndexToReverseDirection(edge._unitigIndex) << endl;
				successors.push_back(Utils::unitigIndexToReverseDirection(edge._unitigIndex));
			}
			
		}

		
	}

	/*
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

			
		//if(isPalindomre){
		//	cout << "palouf" << endl;
			//cout << suffix.toString() << "    " << v.toString() << endl;
		//}


		if(v==vec){
			//cout << "mioum 1" << endl;
			continue;
		}

		KmerVec v_rev = v.reverse();
		
		if(v_rev==vec){
			//cout << "mioum 2" << endl;
			continue;
		}

		//u_int32_t id2 = edge._nodeName;

		if (vec.suffix() == v.prefix()) {

			//cout << "aa: " << unitigIndexFrom << " -> " << edge._unitigIndex << endl;
			//if(isPalindomre){
			//	cout << "suf 1 " << edge._unitigIndex << endl;	
			//}
			//cout << "\t\tSucc1: " << edge._isReversed << " " << edge._unitigIndex << endl;
			successors.push_back({unitigIndexFrom, edge._unitigIndex});

			//dumpEdge(id, edgePlus, id2, edgePlus);
		}
		if (vec.suffix() == v_rev.prefix()) {
			//cout << "bb: " << unitigIndexFrom << " -> " << Utils::unitigIndexToReverseDirection(edge._unitigIndex) << endl;
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
			cout << "c" << endl;
			getchar();
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
			cout << "d" << endl;
			getchar();
			//if(isPalindomre){
			//	cout << "suf 4 " << edge._unitigIndex << endl;	
			//	cout << vec_rev.toString() << "         " <<  v_rev.toString() << endl;
			//}
			//cout << "\t\tSucc4: " << edge._isReversed << " " << edge._unitigIndex << endl;
			successors.push_back({Utils::unitigIndexToReverseDirection(unitigIndexFrom), Utils::unitigIndexToReverseDirection(edge._unitigIndex)});
			//dumpEdge(id, edgeMinus, id2, edgeMinus);
		}
		
	}

	if(successors.size() != successors2.size()){
		cout << "derp succ: " << successors.size() << " " << successors2.size() << endl;
		//getchar();
	}
	*/
}

void CreateMdbg::getPredecessors_unitig(const KmerVec& vec, const UnitigType& unitigIndexFrom, vector<UnitigType>& predecessors){

	bool isPalouf = false;
	//cout << "----" << endl;
	//cout << vec.toString() << endl;

	predecessors.clear();
	//KmerVec vec = keys[i];
	KmerVec vec_rev = vec.reverse();
	KmerVec vec_rev_suffix = vec_rev.suffix();

	
	const KmerVec prefix = vec.prefix().normalize();
	//cout << "\t" << prefix.toString() << endl;


	const u_int128_t key = _unitigEdgeMap._keys->lookup(prefix.hash128());
	const vector<KminmerEdgeU>& edges = _unitigEdgeMap._values[key];

	//if(_mdbgEdges4.find(prefix.hash128()) == _mdbgEdges4.end()){
	//	//cout << "derp pred" << endl;
	//	return;
	//}
	
	//cout << "\tRofl preds: " << _mdbgEdges3[prefix.hash128()].size() << endl;


	//vector<UnitigEdge> predecessors2;


	for(const KminmerEdgeU& edge : edges){


		KmerVec v = prefix;
		if(edge._isReversed){
			v = v.reverse();
		}
		
		
		if(edge._isPrefix){

			if (vec_rev_suffix == v) {
				if(unitigIndexFrom == edge._unitigIndex) continue;
				//if(Utils::unitigIndexToReverseDirection(unitigIndexFrom) == edge._unitigIndex && (vec.suffix().isPalindrome() || vec.prefix().isPalindrome())) continue;
				//cout << "a: " << Utils::unitigIndexToReverseDirection(unitigIndexFrom) << " -> " << edge._unitigIndex << endl;
				predecessors.push_back(edge._unitigIndex);
			}

		}
		else{

			KmerVec v_rev = v.reverse();

			if (vec_rev_suffix == v_rev) {
				if(unitigIndexFrom == Utils::unitigIndexToReverseDirection(edge._unitigIndex)) continue;
				//if(unitigIndexFrom == edge._unitigIndex && (vec.suffix().isPalindrome() || vec.prefix().isPalindrome())) continue;
				//cout << "b: " << Utils::unitigIndexToReverseDirection(unitigIndexFrom) << " -> " << Utils::unitigIndexToReverseDirection(edge._unitigIndex) << endl;
				predecessors.push_back(Utils::unitigIndexToReverseDirection(edge._unitigIndex));
			}



		}

	}
	
	/*
	for(const KminmerEdge4& edge : _mdbgEdges4[prefix.hash128()]){

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
			//cout << "aa: " << Utils::unitigIndexToReverseDirection(unitigIndexFrom) << " -> " << edge._unitigIndex << endl;
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
			//cout << "bb: " << Utils::unitigIndexToReverseDirection(unitigIndexFrom) << " -> " << Utils::unitigIndexToReverseDirection(edge._unitigIndex) << endl;
			predecessors.push_back({Utils::unitigIndexToReverseDirection(unitigIndexFrom), Utils::unitigIndexToReverseDirection(edge._unitigIndex)});
			//dumpEdge(id, edgeMinus, id2, edgeMinus);
			//predecessors.push_back(v_rev);
		}

	}
	

	if(predecessors.size() != predecessors2.size()){
		cout << "derp pred: " << predecessors.size() << " " << predecessors2.size() << endl;
		//if(!isPalouf) getchar();
	}
	*/

}


void CreateMdbg::dumpUnitigNode(const UnitigType& unitigIndex, const vector<MinimizerType>& unitig, ofstream& outputFile){

	
	//cout << "\tDump unitig: " << unitigIndex << " " << unitig.size() << endl;

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
	
	outputFile.write((const char*)&size, sizeof(size));
	outputFile.write((const char*)&unitig[0], size * sizeof(MinimizerType));
	outputFile.write((const char*)&unitigIndex, sizeof(unitigIndex));
	//_unitigGraphFile_nodes.write((const char*)&nodeName, sizeof(nodeName));

	//for(size_t i=0; i<unitig.size(); i++){
	//	const MinimizerType m = unitig[i];
	//	_checksum_unitigNodes += m*unitig.size()*unitigIndex;
	//}

	//cout << "write unitig: " << unitigIndex << " " << unitig.size() << endl;
}

//void CreateMdbg::dumpUnitigEdge(u_int32_t fromUnitigIndex, u_int32_t toUnitigIndex, bool isSuccessor){
void CreateMdbg::dumpUnitigEdge(const UnitigType& unitigIndexFrom, const vector<UnitigType>& successors, const vector<UnitigType>& predecessors){

	//vector<UnitigType> edges;
	//for(const UnitigEdge& edge : successors){
	//	edges.push_back(edge._unitigIndexTo);
	//}
	//for(const UnitigEdge& edge : predecessors){
	//	edges.push_back(edge._unitigIndexTo);
	//}

	UnitigType unitigIndexFrom_rev = Utils::unitigIndexToReverseDirection(unitigIndexFrom);

	#pragma omp critical(dumpUnitigEdge)
	{
		/*
		cout << "---" << endl;
		for(const UnitigEdge& edge : successors){
			cout << "a " << edge._unitigIndexFrom << " " << edge._unitigIndexTo << endl;
		}
		for(const UnitigEdge& edge : predecessors){
			cout << "b " << edge._unitigIndexFrom << " " << edge._unitigIndexTo << endl;
		}
		*/

		u_int32_t nbSuccessors = successors.size();
		u_int32_t nbPredecessors = predecessors.size();
		//if(fromUnitigIndex == 1916 && toUnitigIndex == 4914){
		//	cout << "Add edge: " << fromUnitigIndex << " -> " << toUnitigIndex << endl;
		//	cout << "derp" << endl;
		//	getchar(); 
		//}
        //if(fromUnitigIndex == 1916 || toUnitigIndex == 1916 || fromUnitigIndex == 1917 || toUnitigIndex == 1917){
        //    cout << "Add edge: " << fromUnitigIndex << " -> " << toUnitigIndex << endl;
        //}

		//getchar();

		//_unitigGraphFile_edges_successors.write((const char*)&fromUnitigIndex, sizeof(fromUnitigIndex));
		//_unitigGraphFile_edges_successors.write((const char*)&toUnitigIndex, sizeof(toUnitigIndex));

		_unitigGraphFile_edges_successors.write((const char*)&unitigIndexFrom, sizeof(unitigIndexFrom));
		_unitigGraphFile_edges_successors.write((const char*)&nbSuccessors, sizeof(nbSuccessors));
		_unitigGraphFile_edges_successors.write((const char*)&successors[0], successors.size() * sizeof(UnitigType));
		_unitigGraphFile_edges_successors.write((const char*)&nbPredecessors, sizeof(nbPredecessors));
		_unitigGraphFile_edges_successors.write((const char*)&predecessors[0], predecessors.size() * sizeof(UnitigType));

		for(const UnitigType& unitigIndexTo : successors){
			_nbUnitigEdges += 1;
			_checksum_unitigEdges += unitigIndexFrom * unitigIndexTo;
		}
		for(const UnitigType& unitigIndexTo : predecessors){
			_nbUnitigEdges += 1;
			_checksum_unitigEdges += unitigIndexFrom_rev * unitigIndexTo;
		}

	}

}



void CreateMdbg::indexUnitigEdges(){

	
	Logger::get().debug() << "\tDereplicating unitig edges";
	auto start = high_resolution_clock::now();
	UnitigEdgeIndexer edgeIndexer(*this);
	edgeIndexer.execute();
	Logger::get().debug() << "\t\tDone: " << edgeIndexer._nbEdges << " " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";


	const u_int64_t nbEdges = edgeIndexer._nbEdges;
	const string edgeFilename = edgeIndexer.getOutputFilename();

	Logger::get().debug() << "\tBuilding mphf: " << nbEdges;
	start = high_resolution_clock::now();
	boomphf::file_binary<u_int128_t> inputFile(edgeFilename.c_str());

	double gammaFactor = 1.0; // lowest bit/elem is achieved with gamma=1, higher values lead to larger mphf but faster construction/query
	_unitigEdgeMap._keys = new boomphf::mphf<u_int128_t, hasher_t>(nbEdges, inputFile, _nbCores, gammaFactor, true, false);
	Logger::get().debug() << "\t\tDone: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";

	Logger::get().debug() << "\tCreating edge values";
	_unitigEdgeMap._values.resize(nbEdges);
	Logger::get().debug() << "\t\tDone: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";

	Logger::get().debug() << "\tIndexing edges";

	ifstream nodeFile(_outputDir + "/unitigGraph.nodes.bin");

	
	size_t i = 0;

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

				i+= 1;
				//if(i % 1000000 == 0) cout << "\t" << i << " " << _mdbgEdges4.size() << "    " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << endl;
				//cout << nodeName << endl;
			}

			if(isEOF) break;
			
			indexUnitigEdge(unitigIndex, unitig);
			

		}
	}

	nodeFile.close();
	fs::remove(edgeFilename);

	Logger::get().debug() << "\t\tDone: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";

}
	

void CreateMdbg::indexUnitigEdge(const UnitigType unitigIndex, const vector<MinimizerType>& minimizers){

	//cout << "Index unitig edge: utg" << unitigIndex << endl;
	//for(u_int64_t minimizer : minimizers){
	//	cout << "\t" << minimizer << endl;
	//}

	u_int32_t unitigIndexRC = unitigIndex + 1;
	const vector<KmerVec>& nodes = MDBG::minimizersToKminmers(minimizers, _kminmerSize);

	const KmerVec& startNode = nodes[0];
	const KmerVec& endNode = nodes[nodes.size()-1];
	const KmerVec& startNodeRC = endNode.reverse();
	const KmerVec& endNodeRC = startNode.reverse();

	

	

	bool isReversed;
	const KmerVec& startNodeNorm = startNode.normalize(isReversed);
	
	
	if(isReversed){
		indexEdgeUnitig(startNodeNorm, unitigIndexRC, false);
	}
	else{
		indexEdgeUnitig(startNodeNorm, unitigIndex, true);
	}

	if(startNode == endNode) return;
	
	const KmerVec& endNodeNorm = endNode.normalize(isReversed);

	//cout << "\tEnd node: " << endNodeNorm.toString() << " " << isReversed << endl;

	if(isReversed){
		indexEdgeUnitig(endNodeNorm, unitigIndexRC, true);
	}
	else{
		indexEdgeUnitig(endNodeNorm, unitigIndex, false);
	}


}
	

void CreateMdbg::indexEdgeUnitig(const KmerVec& vec, UnitigType unitigIndex, bool indexPrefix){


	bool isReversedPrefix;
	bool isReversedSuffix;
	KmerVec prefix = vec.prefix().normalize(isReversedPrefix);
	KmerVec suffix = vec.suffix().normalize(isReversedSuffix);
	MinimizerType prefix_minimizer = vec._kmers[vec._kmers.size()-1];
	MinimizerType suffix_minimizer = vec._kmers[0];
	
	const u_int128_t prefixHash = prefix.hash128();
	const u_int128_t suffixHash = suffix.hash128();



	//cout << isReversedPrefix << " " << isReversedSuffix << endl;
	//getchar();
	/*
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
	*/


	//if(indexPrefix){
	u_int128_t key_prefix = _unitigEdgeMap._keys->lookup(prefixHash);
	int partition_prefix = prefixHash % _mutexes.size();
	
	omp_set_lock(&_mutexes[partition_prefix]);

	vector<KminmerEdgeU>& edges_prefix = _unitigEdgeMap._values[key_prefix];
	edges_prefix.push_back({unitigIndex, isReversedPrefix, true});
	
	omp_unset_lock(&_mutexes[partition_prefix]);
	//}
	//else{
	u_int128_t key_suffix = _unitigEdgeMap._keys->lookup(suffixHash);
	int partition_suffix = suffixHash % _mutexes.size();

	omp_set_lock(&_mutexes[partition_suffix]);

	vector<KminmerEdgeU>& edges_suffix = _unitigEdgeMap._values[key_suffix];
	edges_suffix.push_back({unitigIndex, isReversedSuffix, false});
	
	omp_unset_lock(&_mutexes[partition_suffix]);
	//}





}



void CreateMdbg::computeUnitigEdges(){

	_nbUnitigEdges = 0;
	
	ifstream nodeFile(_outputDir + "/unitigGraph.nodes.bin");

	//cout << "single core here" << endl;

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
	
	Logger::get().debug() << "Nb unitig edges: " << _nbUnitigEdges;
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


void CreateMdbg::computeUnitigEdge(const UnitigType& unitigIndex, const vector<MinimizerType>& minimizers){

	
	u_int32_t unitigIndexRC = unitigIndex + 1;
	const vector<KmerVec>& nodes = MDBG::minimizersToKminmers(minimizers, _kminmerSize);

	
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
	vector<UnitigType> successors;
	getSuccessors_unitig(endNodeNorm, unitigIndex, successors);

	//for(const UnitigEdge& edge : successors){

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
		//dumpUnitigEdge(edge._unitigIndexFrom, edge._unitigIndexTo, true);
	//}



	//cout << "computeUnitigEdge_predecessors " << unitigIndex << endl;
	const KmerVec& startNodeNorm = startNode;//.normalize(isReversed);

	vector<UnitigType> predecessors;
	getPredecessors_unitig(startNodeNorm, unitigIndex, predecessors);


	dumpUnitigEdge(unitigIndex, successors, predecessors);
	//for(const UnitigEdge& edge : predecessors){

		//cout << " \t" << toUnitigIndex << endl;
		//u_int32_t unitigIndexFrom = unitigIndex;


		//dumpUnitigEdge(edge._unitigIndexFrom, edge._unitigIndexTo, true);
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
	//}

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
	/*
	cout << "TEST: " << _nbKminmersTotal << endl;
	u_int64_t nbKminmers2 = 0;
	file_binary_kminmerAbundance inputFile((_outputDir + "/kminmerData_abundance.txt").c_str());
	boophf_t* lala = new boomphf::mphf<u_int128_t, hasher_t>(_nbKminmersTotal, inputFile, _nbCores, 1.0);
	cout << "TEST" << endl;
	*/

	u_int64_t nbKminmers = 0;
	u_int64_t nbSolidKminmers = 0;
	u_int64_t checksumKminmerAbundance = 0;
	u_int64_t nbUnitigs = 0;

	KminmerAbundanceMap().swap(_kminmerAbundances); //_kminmerAbundances.clear();
	
	ifstream kminmerAbundanceFile(_outputDir + "/kminmerData_abundance.txt");

	//bool isEOF = false;
	//KmerVec vec;
	//u_int32_t nodeName;


	while (true) {

		u_int128_t vecHash;
		AbundanceType abundance;

		bool iseof = MDBG::readKminmerAbundance(vecHash, abundance, kminmerAbundanceFile);

		if(iseof) break;

		checksumKminmerAbundance += abundance * vecHash;
		nbKminmers += 1;

		if(abundance == 1) continue;

		_kminmerAbundances[vecHash] = abundance;
		nbSolidKminmers += 1;

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

		vector<u_int32_t> abundances;
		const vector<KmerVec>& kminmers = MDBG::minimizersToKminmers(unitig, _kminmerSize);

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
		_unitigGraphFile_nodes_abundances.write((const char*)&abundances[0], nbAbundances * sizeof(u_int32_t));


		//cout << "write unitig: " << unitigIndex << " " << abundances.size() << endl;

		for(size_t i=0; i<unitig.size(); i++){
			const MinimizerType m = unitig[i];
			_checksum_unitigNodes += m*unitig.size()*unitigIndex;
		}

		for(size_t i=0; i<abundances.size(); i++){
			_checksum_unitigAbundances += ((u_int64_t)abundances[i])* ((u_int64_t)abundances.size());
		}

		nbUnitigs += 1;
	}

	nodeFile.close();
	
	KminmerAbundanceMap().swap(_kminmerAbundances); //_kminmerAbundances.clear();


	//cout << "Nb kminmers: " << nbKminmers << endl;
	//cout << "Nb solid kminmers: " << nbSolidKminmers << endl;
	Logger::get().debug() << "Checksum kminmer abundance: " << checksumKminmerAbundance;
	//cout << "Nb unitigs: " << nbUnitigs << endl;
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

	//KmerVec v;
	//v._kmers = {34886160748729436, 85562229786978109, 51808891639906056, 88809709545039515, 10937167467189174};
	/*
	if(_kminmerSize > _kminmerSizeFirst+1){
		
		
		ifstream kminmerAbundanceFile(_outputDir + "/kminmerData_abundance_init_k" + to_string(_kminmerSizeFirst+1) + ".txt");
		//ifstream kminmerAbundanceFile(_inputDir + "/kminmerData_abundance_init.txt");

		while (true) {

			u_int128_t vecHash;
			u_int32_t abundance;

			kminmerAbundanceFile.read((char*)&vecHash, sizeof(vecHash));

			if(kminmerAbundanceFile.eof()) break;

			kminmerAbundanceFile.read((char*)&abundance, sizeof(abundance));
			
			//return false;

			//bool iseof = MDBG::readKminmerAbundance(vecHash, abundance, kminmerAbundanceFile);

			//if(iseof) break;

			if(abundance <= 1) continue;
			//if(_kminmer_to_unitigIndex.find(vecHash) == _kminmer_to_unitigIndex.end()) continue;

			_kminmerAbundances[vecHash] = abundance;



		}

		kminmerAbundanceFile.close();

		return;
		
	}
	*/

	u_int64_t nbNodesOriginal = 0;
	u_int64_t sumAbundanceOriginal = 0;


	ifstream kminmerAbundanceFile(_outputDir + "/kminmerData_abundance_prev.txt");
	//u_int64_t checksum = 0;

	while (true) {

		u_int128_t vecHash;
		AbundanceType abundance;

		bool iseof = MDBG::readKminmerAbundance(vecHash, abundance, kminmerAbundanceFile);

		if(iseof) break;

		//nbNodesOriginal += 1;
		sumAbundanceOriginal += abundance;

		//if(abundance == 1) continue;

		//if(vecHash == v.normalize().hash128()){
		//	cout << "oyé" << endl;
		//	exit(1);
		//}

		if(abundance == 1) continue;
		_kminmerAbundances[vecHash] = abundance;


	}

	kminmerAbundanceFile.close();
	
	//cout << "skip refined" << endl;
	//return;
	Logger::get().debug() << "\tDone: " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";
	//cout << "Nb nodes original: " << _kminmerAbundances.size() << endl;
	//if(_kminmerSize > _kminmerSizeFirst+1) return;

	ankerl::unordered_dense::map<UnitigType, u_int32_t> unitigName_to_abundance;

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

	Logger::get().debug() << "\tDone: " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";


	ifstream nodeFile(_outputDir + "/unitigGraph_prev.nodes.bin");

	#pragma omp parallel num_threads(_nbCores)
	{

		bool isEOF = false;
		UnitigType unitigIndex;
		vector<MinimizerType> minimizers;
		//KmerVec vec;

		while (true) {

			#pragma omp critical(loadRefinedAbundances)
			{



				u_int32_t size;
				nodeFile.read((char*)&size, sizeof(size));
				
				isEOF = nodeFile.eof();
				//if(nodeFile.eof()) break;

				if(!isEOF){
					minimizers.resize(size);
					nodeFile.read((char*)&minimizers[0], size * sizeof(MinimizerType));


					nodeFile.read((char*)&unitigIndex, sizeof(unitigIndex));
				}

			}

			if(isEOF) break;
			
			const UnitigType& unitigName = unitigIndex/2;

			if(unitigName_to_abundance.find(unitigName) == unitigName_to_abundance.end()) continue;
			
			const u_int32_t unitigAbundance = unitigName_to_abundance[unitigName];
			//if(unitigAbundance == 1) continue;

			for(const KmerVec& kminmer : MDBG::minimizersToKminmers(minimizers, _kminmerSizePrev)){

				const u_int128_t hash = kminmer.normalize().hash128();

				//if(hash == v.normalize().hash128()){
				//	cout << "oyé" << endl;
				//	exit(1);
				//}


				if(unitigAbundance == 1){
					//if(_kminmerAbundances.find(hash) != _kminmerAbundances.end()){
					//	_kminmerAbundances[hash] = 0;
					//}
					_kminmerAbundances.modify_if(hash, 
						[this](KminmerAbundanceMap::value_type& v) { 
						v.second = 0;	
					});

				}
				else{

					_kminmerAbundances.lazy_emplace_l(hash, 
					[this, &unitigAbundance](KminmerAbundanceMap::value_type& v) { // key exist
						v.second = unitigAbundance;
					},           
					[&hash, this, &unitigAbundance](const KminmerAbundanceMap::constructor& ctor) { // key inserted
						
						//DbgNodeLight node = {};

						ctor(hash, unitigAbundance); 
					});

					//_kminmerAbundances[hash] = unitigAbundance;
				}
				
			}

		}


	}

	/*
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

		for(const KmerVec& kminmer : MDBG::minimizersToKminmers(minimizers, _kminmerSizePrev)){

			const u_int128_t hash = kminmer.normalize().hash128();

			//if(hash == v.normalize().hash128()){
			//	cout << "oyé" << endl;
			//	exit(1);
			//}


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
	*/

	nodeFile.close();
	
	Logger::get().debug() << "\tDone: " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";
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


void CreateMdbg::computeNextUnitigGraph(){

	ifstream nodeFile(_outputDir + "/unitigGraph_prev.nodes.bin");

	while(true){

		u_int32_t size;
		nodeFile.read((char*)&size, sizeof(size));
		
		if(nodeFile.eof()) break;

		vector<MinimizerType> minimizers;
		minimizers.resize(size);
		nodeFile.read((char*)&minimizers[0], size * sizeof(MinimizerType));

		UnitigType unitigIndex;
		nodeFile.read((char*)&unitigIndex, sizeof(unitigIndex));

		const UnitigType unitigName = unitigIndex/2;

        while(_unitigName_to_minimizers.size() <= unitigName){
            _unitigName_to_minimizers.push_back({});
        }

		_unitigName_to_minimizers[unitigName] = minimizers;
	}

	nodeFile.close();


	Logger::get().debug() << "Loading unitig graph";
	UnitigGraph2* unitigGraph = new UnitigGraph2(_kminmerSizePrev, _kminmerLengthMean, _kminmerOverlapMean, _kminmerLengthMean-_kminmerOverlapMean, _minimizerSpacingMean, _nbCores);
	unitigGraph->load(_outputDir + "/unitigGraph_prev.stats.bin", _outputDir + "/unitigGraph_prev.nodes.bin", _outputDir + "/unitigGraph_prev.nodes.abundances.bin", _outputDir + "/unitigGraph_prev.edges.successors.bin");

	//cout << unitigGraph->_kminmerLength << " " <<  unitigGraph->_kminmerSize << " " <<  unitigGraph->_kminmerLengthNonOverlap << endl;

	/*
	for(size_t i=0; i<_unitigName_to_minimizers.size(); i++){
		
		int n = 0;
		for(auto m : _unitigName_to_minimizers[i]){
			if(m == 945004921997704){
				n += 1;
			}
			if(m == 54561592015907157){
				n += 1;
			}
			if(m == 8063797441767696){
				n += 1;
			}
			if(m == 15467155777103272){
				n += 1;
			}
		}

		if(n == 4){
			cout << "hahaha: " << i << endl;
			getchar();
		}
	}
	*/




	//loadReference(_kminmerSizePrev);
	//loadReference(_kminmerSize);
	//unitigGraph->save("/pasteur/appa/homes/gbenoit/zeus/tmp/unitigs/assembly_graph.gfa", _outputDir);
	//writeReferencePositions(unitigGraph, "/pasteur/appa/homes/gbenoit/zeus/tmp/unitigs/assembly_graph.gfa.referencePos.csv");
	//getchar();

	/*
	cout << "utg39" << endl; 
	for(MinimizerType m : _unitigName_to_minimizers[39]){
		cout << m << endl;
	}


	for(const UnitigGraph2::UnitigNode* unitig : unitigGraph->_unitigs){
		
		if(unitig == nullptr) continue;

		UnitigType unitigIndex = UnitigGraph2::unitigName_to_unitigIndex(unitig->_unitigName, false);

		const vector<MinimizerType>& minimizers = _unitigName_to_minimizers[unitig->_unitigName];
		const vector<KmerVec>& kminmers = Utils::minimizersToKminmers(minimizers, _kminmerSize);

		vector<float> abundances;

		for(const KmerVec& kminmer : kminmers){

			const u_int128_t vecHash = kminmer.normalize().hash128();

			if(unitig->_unitigName == 39){
				cout << kminmer.normalize().toString() << " " << (_mdbgNodesLight.find(vecHash) != _mdbgNodesLight.end()) << " " << (_debug_vecHash_to_referencePosition.find(vecHash) != _debug_vecHash_to_referencePosition.end()) << endl;
				//printf("%llx\n", (unsigned long long)(vecHash & 0xFFFFFFFFFFFFFFFF));
				printf("unsinged int 128 bit: %llu\n", (unsigned long long)vecHash);
			}

			//cout << "lul: " << _mdbgNodesLight[vecHash]._abundance << endl;
			if(_mdbgNodesLight.find(vecHash) == _mdbgNodesLight.end()){
				abundances.push_back(1);
			}
			else{
				abundances.push_back(_mdbgNodesLight[vecHash]._abundance);
			}
		}

		const vector<KmerVec>& kminmersPrev = Utils::minimizersToKminmers(minimizers, _kminmerSizePrev);
		for(const KmerVec& kminmer : kminmersPrev){

			const u_int128_t vecHash = kminmer.normalize().hash128();

			if(unitig->_unitigName == 39){
				cout << kminmer.normalize().toString() << " " << (_debug_vecHash_to_referencePosition.find(vecHash) != _debug_vecHash_to_referencePosition.end()) << endl;
			}
		}
	}

	

	for(const KmerVec& kminmer : Utils::minimizersToKminmers(_unitigName_to_minimizers[6], _kminmerSizePrev)){

		cout << kminmer.toString() << endl;
		const u_int128_t hash = kminmer.normalize().hash128();

		if(_debug_vecHash_to_referencePosition.find(hash) == _debug_vecHash_to_referencePosition.end()){
			cout << "Found invalid kminmer" << endl; 
		}
	}
	*/

	
	//cout << "utg344" << endl; 
	//for(MinimizerType m : _unitigName_to_minimizers[344]){
	//	cout << m << endl;
	//}


	//cout << "Nb unitigs: " << unitigGraph->_unitigs.size() << endl;
	/*
	for(UnitigGraph2::UnitigNode* node : unitigGraph->_unitigs){

		if(node == nullptr) continue;

		if(node->_abundance < 10){
			cout << "Removing low cov unitig: utg" << node->_unitigName << endl;
			unitigGraph->removeNode(node);
		}
	}
	*/

	//extendUnitigs(unitigGraph);
	//if(_kminmerSize == _kminmerSizeFirst+1){
	//}
	
	Logger::get().debug() << "Solve edges";
	auto start = high_resolution_clock::now();
	solveEdges(unitigGraph);
	Logger::get().debug() << "\tDuration: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";
	
	
	Logger::get().debug() << "Remove unsupported unitigs";
	start = high_resolution_clock::now();
	removeUnsupportedUnitigs(unitigGraph);
	Logger::get().debug() << "\tDuration: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";
	//unitigGraph->save("/pasteur/appa/homes/gbenoit/zeus/tmp/unitigs/assembly_graph_1.gfa", _outputDir);
	

	Logger::get().debug() << "Solve small unitigs";
	start = high_resolution_clock::now();
	solveSmallUnitigs(unitigGraph);
	Logger::get().debug() << "\tDuration: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";
	
	//cout << "Remove unsupported edge disabled" << endl;

	//exit(1);

	


	Logger::get().debug() << "Write unitigs";
	start = high_resolution_clock::now();
	writeUnitigs(unitigGraph);
	Logger::get().debug() << "\tDuration: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";
	

}



void CreateMdbg::solveEdges(UnitigGraph2* unitigGraph){

	vector<UnitigGraph2::UnitigNode*> nodes;

	for(UnitigGraph2::UnitigNode* node : unitigGraph->_unitigs){

		if(node == nullptr) continue;
		if(node->_nbMinimizers == _kminmerSizePrev) continue; //Small unitig that going to be removed

		nodes.push_back(node);
		//cout << node->_unitigName << endl;
	

	}

	unordered_set<UnitigType> processedUnitigNames;
	
	for(UnitigGraph2::UnitigNode* node : nodes){
		const UnitigType& unitigIndex1 = UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, false);
		solveEdge(unitigGraph, node->_unitigName, unitigIndex1, processedUnitigNames);

		const UnitigType& unitigIndex2 = UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, true);
		solveEdge(unitigGraph, node->_unitigName, unitigIndex2, processedUnitigNames);

		processedUnitigNames.insert(node->_unitigName);
	}
}


void CreateMdbg::solveEdge(UnitigGraph2* unitigGraph, const UnitigType unitigName, const UnitigType unitigIndex, unordered_set<UnitigType>& processedUnitigNames){

	//cout << "Solve edge: " << unitigName << endl; 

	//int expectedOverlapSize = _kminmerSizePrev-1;

	//vector<MinimizerType> minimizers = UnitigGraph2::getMinimizers(unitigIndex, _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(unitigIndex)]);

	//cout << "\tutg" << UnitigGraph2::unitigIndexToString(unitigIndex) << " " << minimizers.size() << endl;

	vector<UnitigType> successors;
	unitigGraph->getSuccessors(unitigIndex, successors);

	//if(unitigName == 88){
	//	cout << "\tutg" << UnitigGraph2::unitigIndexToString(unitigIndex) << endl;
	//	for(MinimizerType m : minimizers){
	//		cout << "\t\t" << m << endl;
	//	}
	//}
	
	//vector<MinimizerType> sourceEndNode = UnitigGraph2::getEndNode(unitigIndex, _unitigName_to_minimizers[unitigName], _kminmerSizePrev);

	for(size_t i=0; i<successors.size(); i++){
		UnitigType unitigIndexSucc = successors[i];
		UnitigType unitigNameSucc = UnitigGraph2::unitigIndex_to_unitigName(unitigIndexSucc);

		UnitigGraph2::UnitigNode* nodeSuccessor = unitigGraph->_unitigs[unitigNameSucc];
		if(nodeSuccessor->_nbMinimizers == _kminmerSizePrev) continue; //Successor is a small unitig that going to be removed
		if(processedUnitigNames.find(unitigNameSucc) != processedUnitigNames.end()) continue;

		//cout << "\t" << unitigNameSucc << endl;
		//UnitigType unitigNameSucc = UnitigGraph2::unitigIndex_to_unitigName(unitigIndexSucc);
		//vector<MinimizerType> startNode = UnitigGraph2::getStartNode(unitigIndexSucc, _unitigName_to_minimizers[unitigNameSucc], _kminmerSizePrev);
	

		vector<MinimizerType> doublet;
		getDoublet2(unitigGraph, unitigIndex, unitigIndexSucc, doublet);

		KmerVec kminmer;
		kminmer._kmers = doublet;
		kminmer = kminmer.normalize();
		u_int128_t vecHash = kminmer.hash128();

		//vector<MinimizerType> succMinimizers = UnitigGraph2::getMinimizers(unitigIndexSucc, _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(unitigIndexSucc)]);
		
		//if(unitigName == 88){
		//	cout << "\t\tutg" << UnitigGraph2::unitigIndexToString(unitigIndexSucc) << endl;
		//	for(MinimizerType m : succMinimizers){
		//		cout << "\t\t\t" << m << endl;
		//	}
		//}

		//cout << "\t\tutg" << UnitigGraph2::unitigIndexToString(unitigIndexSucc) << " " << succMinimizers.size() << endl;
		//if(unitigName == 88){
		//	cout << "\t\tCheck edge: utg" << UnitigGraph2::unitigIndexToString(unitigIndex) << " -> " << UnitigGraph2::unitigIndexToString(unitigIndexSucc) << endl;
		//}
		

		if(isEdgeSupported(vecHash, _mdbgNodesLight)){

			createDoubletNode(unitigGraph, unitigIndex, unitigIndexSucc, doublet, processedUnitigNames);
			//cout << "\t\t\tEdge unsupported: utg" << UnitigGraph2::unitigIndexToString(unitigIndex) << " -> utg" << UnitigGraph2::unitigIndexToString(unitigIndexSucc) << endl;
		
		}

		unitigGraph->removeSuccessor(unitigIndex, unitigIndexSucc);
		unitigGraph->removeSuccessor(UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndexSucc), UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndex));


		/*
		int overlapSize = longestOverlap(minimizers, succMinimizers);

		if(overlapSize != expectedOverlapSize) continue;

		vector<MinimizerType> edge;

		//cout << "---" << endl;
		//cout << _kminmerSize << " " << expectedOverlapSize << " " << minimizers.size()-expectedOverlapSize << " " << minimizers.size() << endl;
		for(size_t i=minimizers.size()-expectedOverlapSize-1; i<minimizers.size(); i++){
			edge.push_back(minimizers[i]);
		}
		//cout << edge.size() << endl;
		edge.push_back(succMinimizers[expectedOverlapSize]);

		KmerVec kminmer;
		kminmer._kmers = edge;
		kminmer = kminmer.normalize();

		//cout << kminmer.toString() << endl;

		//if(_mdbgNodesLight.find(kminmer.hash128()) == _mdbgNodesLight.end()){
		//	cout << "Actual: Unsupported: utg" << UnitigGraph2::unitigIndexToString(unitigIndex) << " -> utg" << UnitigGraph2::unitigIndexToString(unitigIndexSucc) << endl;
		//}
		*/
		//getchar();
		/*
		int overlapSize = longestOverlap(minimizers, succMinimizers);

		if(overlapSize != expectedOverlapSize) continue;

		vector<MinimizerType> edge;

		//cout << "---" << endl;
		//cout << _kminmerSize << " " << expectedOverlapSize << " " << minimizers.size()-expectedOverlapSize << " " << minimizers.size() << endl;
		for(size_t i=minimizers.size()-expectedOverlapSize-1; i<minimizers.size(); i++){
			edge.push_back(minimizers[i]);
		}
		//cout << edge.size() << endl;
		edge.push_back(succMinimizers[expectedOverlapSize]);


		//cout << edge.size() << endl;
		//vector<MinimizerType> edge = sourceEndNode;
		//edge.push_back(startNode[startNode.size()-1]);

		KmerVec kminmer;
		kminmer._kmers = edge;
		kminmer = kminmer.normalize();

		if(_mdbgNodesLight.find(kminmer.hash128()) == _mdbgNodesLight.end()){
			unitigGraph->removeSuccessor(unitigIndex, unitigIndexSucc);
			unitigGraph->removeSuccessor(UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndexSucc), UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndex));

			//cout << "Unsupported: utg" << UnitigGraph2::unitigIndexToString(unitigIndex) << " -> utg" << UnitigGraph2::unitigIndexToString(unitigIndexSucc) << endl;
		}
		*/
	}

}


void CreateMdbg::createDoubletNode(UnitigGraph2* unitigGraph, const UnitigType unitigIndex, const UnitigType unitigIndexSuccessor, const vector<MinimizerType>& minimizers, unordered_set<UnitigType>& processedUnitigNames){


	pair<UnitigType, UnitigType> unitigPair = {unitigIndex, unitigIndexSuccessor};

	//if(edgeNodes.find(unitigPair) != edgeNodes.end()){
	//	return;// edgeNodes[unitigPair];
	//}

	//UnitigGraph2::UnitigNode* node = unitigGraph->_unitigs[UnitigGraph2::unitigIndex_to_unitigName(unitigIndex)];
	UnitigGraph2::UnitigNode* successor = unitigGraph->_unitigs[UnitigGraph2::unitigIndex_to_unitigName(unitigIndexSuccessor)];
	//vector<MinimizerType> predMinimizers = UnitigGraph2::getMinimizers(unitigIndex, _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(unitigIndex)]);
	//vector<MinimizerType> succMinimizers = UnitigGraph2::getMinimizers(unitigIndexSuccessor, _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(unitigIndexSuccessor)]);
	//vector<MinimizerType> minimizers = UnitigGraph2::getMinimizers(unitigIndex, _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(unitigIndex)]);

	//const vector<MinimizerType> minimizers = getDoublet(unitigGraph, unitigIndex, predMinimizers, unitigIndexSuccessor, succMinimizers, _kminmerSize);
	//int overlapSize = longestOverlap(minimizers, succMinimizers);

	//bool connectTwoSmallNodes = successor->_nbMinimizers == _kminmerSizePrev;

	//vector<MinimizerType> minimizers = _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(unitigIndex)];
	//if(!connectTwoSmallNodes){
	//minimizers.push_back(succMinimizers[overlapSize]);
	//}

	//if(isUnsupportedUnitig(minimizers)) return;

	UnitigGraph2::UnitigNode* edgeNode = createEdgeNode(unitigGraph, minimizers);
	processedUnitigNames.insert(edgeNode->_unitigName);
	//UnitigGraph2::UnitigNode* edgeNode = unitigGraph->addNode(unitigGraph->_unitigs.size(), minimizers);
	//edgeNodes[unitigPair] = edgeNode;
	//edgeNode->_abundance = 40;
	//_unitigName_to_minimizers[edgeNode->_unitigName] = minimizers;

	const UnitigType edgeNodeUnitigIndex = UnitigGraph2::unitigName_to_unitigIndex(edgeNode->_unitigName, false);
	
	//if(edgeNode->_unitigName == 200){
	//	cout << "lala" << endl;
	//	getchar();
	//}
	
	unitigGraph->addSuccessor(unitigIndex, edgeNodeUnitigIndex);
	unitigGraph->addSuccessor(UnitigGraph2::unitigIndex_to_reverseDirection(edgeNodeUnitigIndex), UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndex));
	
	unitigGraph->addSuccessor(edgeNodeUnitigIndex, unitigIndexSuccessor);
	unitigGraph->addSuccessor(UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndexSuccessor), UnitigGraph2::unitigIndex_to_reverseDirection(edgeNodeUnitigIndex));

	//cout << "\tCreate edge node S: " << UnitigGraph2::unitigIndexToString(unitigIndex) << " -> " << UnitigGraph2::unitigIndexToString(unitigIndexSuccessor) << "    utg" << edgeNode->_unitigName << endl;
	//for(MinimizerType m : minimizers){
	//	cout << "\t\t" << m << endl;
	//}


	//vector<UnitigType> successors;
	//unitigGraph->getSuccessors(edgeNodeUnitigIndex, successors);

	//for(UnitigType unitigIndex : successors){
	//	cout << "\t\tutg" << UnitigGraph2::unitigIndex_to_unitigName(unitigIndex) << endl;
	//	for(MinimizerType m : UnitigGraph2::getMinimizers(unitigIndex, _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(unitigIndex)])){
	//		cout << "\t\t\t" << m << endl;
	//	}
	//}

	/*
	if(edgeNode->_unitigName == 195){
		cout << "lul" << endl;
		for(MinimizerType m : _unitigName_to_minimizers[edgeNode->_unitigName]){
			cout << "\t\t\t" << m << endl;
		}
		getchar();
	}
	*/

}

void CreateMdbg::removeUnsupportedUnitigs(UnitigGraph2* unitigGraph){

	vector<UnitigGraph2::UnitigNode*> nodesToRemove;



	#pragma omp parallel for num_threads(_nbCores)
	for(size_t i=0; i<unitigGraph->_unitigs.size(); i++){

		UnitigGraph2::UnitigNode* node = unitigGraph->_unitigs[i];
		if(node == nullptr) continue;

		const vector<MinimizerType>& minimizers = _unitigName_to_minimizers[node->_unitigName];

		for(const KmerVec& kminmer : MDBG::minimizersToKminmers(minimizers, _kminmerSize)){

			const u_int128_t hash = kminmer.normalize().hash128();
			
			if(_mdbgNodesLight.find(hash) != _mdbgNodesLight.end()){
				//if(_kminmerSize == _kminmerSizeFirst+1){
				//	if (_mdbgNodesLight[hash]._abundance > 10 && _mdbgNodesLight[hash]._index <= 1) {
				//		cout << "Found chimeric kminmer" << endl;
				//		return true;
				//	}
				//}
			}
			else{
				//cout << "Found unsupported unitig" << endl;

				#pragma omp critical(removeUnsupportedUnitigs)
				{
					nodesToRemove.push_back(node);
				}
				break;
				//return true;
			}
		}


	}

	
	for(UnitigGraph2::UnitigNode* node : nodesToRemove){
		//cout << "Remove node: " << node->_nbMinimizers << " " << node->_abundance << endl;
		unitigGraph->removeNode(node);
	}

	/*
	for(UnitigGraph2::UnitigNode* node : nodesToRemove){
		
		cout << "\tRemove unsupported unitig: utg" << node->_unitigName << endl;

		const UnitigType& unitigIndex = UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, false);

		vector<UnitigType> predecessors;
		unitigGraph->getPredecessors(unitigIndex, predecessors);

		vector<UnitigType> successors;
		unitigGraph->getSuccessors(unitigIndex, successors);

		const vector<MinimizerType>& minimizers = _unitigName_to_minimizers[node->_unitigName];
		const vector<KmerVec> kminmers = MDBG::minimizersToKminmers(minimizers, _kminmerSize);
		
		bool isConnectedToPredecessors = false;
		bool isConnectedToSuccessors = false;
		
		vector<vector<KmerVec>> unitigs;

		vector<KmerVec> unitig;

		unitigGraph->removeNode(node);
		cout << "\t\tNb kminmers: " << kminmers.size() << endl;

		for(size_t i=0; i<kminmers.size(); i++){
			
			const KmerVec& vec = kminmers[i];
			const u_int128_t hash = vec.normalize().hash128();
			
			if(_mdbgNodesLight.find(hash) != _mdbgNodesLight.end()){
				cout << "\t\t\t1" << endl;
				if(i==0){
					isConnectedToPredecessors = true;
				}
				else if(i == kminmers.size()-1){
					isConnectedToSuccessors = true;
				}

				unitig.push_back(vec);
			}
			else{
				cout << "\t\t\t0" << endl;
				if(unitig.size() > 0){
					unitigs.push_back(unitig);
					unitig.clear();
				}
			}

		}

		if(unitig.size() > 0){
			unitigs.push_back(unitig);
			unitig.clear();
		}

		for(size_t i=0; i<unitigs.size(); i++){
			const vector<MinimizerType>& unitigMinimizers = kminmersToMinimizers(unitigs[i]);
			
			UnitigGraph2::UnitigNode* node = createEdgeNode(unitigGraph, unitigMinimizers);

			//UnitigGraph2::UnitigNode* node = unitigGraph->addNode(unitigGraph->_unitigs.size(), unitigMinimizers);
			const UnitigType nodeUnitigIndex = UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, false);

			cout << "\t\tCreate subnode: utg" << node->_unitigName << " " << node->_nbMinimizers << endl;  

			if(i == 0 && isConnectedToPredecessors){
				for(const UnitigType unitigIndexPred : predecessors){
					if(unitigGraph->_unitigs[UnitigGraph2::unitigIndex_to_unitigName(unitigIndexPred)] == nullptr) continue;
					unitigGraph->addSuccessor(unitigIndexPred, nodeUnitigIndex);
					unitigGraph->addSuccessor(UnitigGraph2::unitigIndex_to_reverseDirection(nodeUnitigIndex), UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndexPred));
				}
			}

			if(i == unitigs.size()-1 && isConnectedToSuccessors){
				for(const UnitigType unitigIndexSucc : successors){
					if(unitigGraph->_unitigs[UnitigGraph2::unitigIndex_to_unitigName(unitigIndexSucc)] == nullptr) continue;
					unitigGraph->addSuccessor(nodeUnitigIndex, unitigIndexSucc);
					unitigGraph->addSuccessor(UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndexSucc), UnitigGraph2::unitigIndex_to_reverseDirection(nodeUnitigIndex));
				}
			}

		}
	}
	*/
	Logger::get().debug() << "\tRemove unsupported unitigs: " << nodesToRemove.size();
}


void CreateMdbg::solveSmallUnitigs(UnitigGraph2* unitigGraph){

	//unordered_set<UnitigType> isUnitigNameToRemove;
	vector<UnitigGraph2::UnitigNode*> nodeToSolve;

	for(UnitigGraph2::UnitigNode* node : unitigGraph->_unitigs){

		if(node == nullptr) continue;
		if(node->_nbMinimizers != _kminmerSizePrev) continue;

		//solveSmallUnitigsSub(unitigGraph, node);
		//const UnitigType& unitigIndex1 = UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, false);
		
		//if(unitigGraph->nbSuccessors(unitigIndex1) == 1 || unitigGraph->nbPredecessors(unitigIndex1) == 1) continue;

		//solveSmallUnitigsSub(unitigGraph, node);
		nodeToSolve.push_back(node);
		//isUnitigNameToRemove.insert(node->_unitigName);

	}

	//cout << "Nb small unitigs: " << nodesToRemove.size() << endl;

	for(size_t i=0; i<nodeToSolve.size(); i++){
		solveSmallUnitigsSub2(unitigGraph, nodeToSolve[i]);
	}

}

/*
void CreateMdbg::solveSmallUnitigsSub(UnitigGraph2* unitigGraph, UnitigGraph2::UnitigNode* node){
	
	unordered_set<UnitigType> isNodeExtended;
	
	//if(node->_unitigName == 24 || node->_unitigName == 25 || node->_unitigName == 44 || node->_unitigName == 56){
	//if(node->_unitigName == 24 || node->_unitigName == 25 || node->_unitigName == 56){
	//if(node->_unitigName == 25){
	//if(node->_unitigName == 26){

	//}
	//else{
	//	return;
	//}
	
	//if(node->_unitigName == 29){

	//}
	//else{
	//	return;
	//}
	

	vector<vector<UnitigType>> supportedTriplets;
	phmap::parallel_flat_hash_map<pair<UnitigType, UnitigType>, UnitigGraph2::UnitigNode*> edgeNodes;

	//cout << node->_unitigName << endl;

	const UnitigType& unitigIndex1 = UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, false);

	vector<UnitigType> predecessors;
	unitigGraph->getPredecessors(unitigIndex1, predecessors);

	vector<UnitigType> successors;
	unitigGraph->getSuccessors(unitigIndex1, successors);


	const vector<MinimizerType> minimizers = _unitigName_to_minimizers[node->_unitigName];
	const UnitigType unitigName = node->_unitigName;

	//cout << endl << endl;
	//cout << "Remove node: utg" << node->_unitigName << endl;
	unitigGraph->removeNode(node);

	//for(MinimizerType m : _unitigName_to_minimizers[unitigName]){
	//	cout << "\t" << m << endl;
	//}

	//cout << "\tutg" << UnitigGraph2::unitigIndexToString(unitigIndex1) << endl;
	//for(MinimizerType m : minimizers){
	//	cout << "\t\t" << m << endl;
	//}

	for(size_t i=0; i<predecessors.size(); i++){

		//cout << "p " << i << " " << predecessors.size() << endl;

		vector<MinimizerType> predMinimizers = UnitigGraph2::getMinimizers(predecessors[i], _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(predecessors[i])]);

		//vector<MinimizerType> endNode = UnitigGraph2::getEndNode(predecessors[i], _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(predecessors[i])], _kminmerSizePrev);

		//cout << "\tutg" << UnitigGraph2::unitigIndexToString(predecessors[i]) << endl;
		//for(MinimizerType m : predMinimizers) {
		//	cout << "\t\t\t" << m << endl;
		//}

		int predecessorOverlapSize = longestOverlap(predMinimizers, minimizers);

		for(size_t j=0; j<successors.size(); j++){

			//cout << "\ts " << j << " " << successors.size() << endl;

			if(predecessors[i] == successors[j]){
				cout << "Skip circular triplet" << endl;
				continue;
			}
			if(predecessors[i] == unitigIndex1){
				cout << "Skip self loop 1" << endl;
				continue;
			}
			if(successors[j] == unitigIndex1){
				cout << "Skip self loop 2" << endl;
				continue;
			}

			vector<MinimizerType> succMinimizers = UnitigGraph2::getMinimizers(successors[j], _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(successors[j])]);
		


			int successorOverlapSize = longestOverlap(minimizers, succMinimizers);

			//cout << predMinimizers.size() << " " << minimizers.size() << " " << succMinimizers.size() << " " << predecessorOverlapSize << " " << successorOverlapSize << endl;
			//if(predecessorOverlapSize != successorOverlapSize){
				//cout << "Overlap predecessor != overlap successor: " << predecessorOverlapSize << " " << successorOverlapSize << endl;
			//}
			
			//vector<MinimizerType> startNode = UnitigGraph2::getStartNode(successors[j], _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(successors[j])], _kminmerSizePrev);
			//cout << "\t\tutg" << UnitigGraph2::unitigIndexToString(successors[j]) << endl;
			//for(MinimizerType m : succMinimizers){
			//	cout << "\t\t\t\t" << m << endl;
			//}


			vector<MinimizerType> triplet;
			//triplet.push_back(endNode[0]);
			triplet.push_back(predMinimizers[predMinimizers.size()-predecessorOverlapSize-1]);
			for(const MinimizerType m : minimizers){
				triplet.push_back(m);
			}
			triplet.push_back(succMinimizers[successorOverlapSize]);
			//triplet.push_back(startNode[startNode.size()-1]);

			


			KmerVec kminmer;
			kminmer._kmers = triplet;
			kminmer = kminmer.normalize();

			if(_solidTriplets.find(kminmer.hash128()) == _solidTriplets.end()) continue;
			if(isUnsupportedUnitig(triplet)) continue;

			
			cout << "\t\tResolvable: " << "utg" << UnitigGraph2::unitigIndexToString(predecessors[i]) << " -> " << "utg" << UnitigGraph2::unitigIndexToString(unitigIndex1) << " -> " << "utg" << UnitigGraph2::unitigIndexToString(successors[j]) << endl;

			cout << "\t\tTriplet:" << endl;
			for(const MinimizerType m : triplet){
				cout << "\t\t\t" << m << endl;
			}

			cout << "\t\t\tutg" << UnitigGraph2::unitigIndexToString(predecessors[i]) << endl;
			for(MinimizerType m : predMinimizers) {
				cout << "\t\t\t\t" << m << endl;
			}
			cout << "\t\t\tutg" << UnitigGraph2::unitigIndexToString(unitigIndex1) << endl;
			for(MinimizerType m : minimizers){
				cout << "\t\t\t\t" << m << endl;
			}
			cout << "\t\t\tutg" << UnitigGraph2::unitigIndexToString(successors[j]) << endl;
			for(MinimizerType m : succMinimizers){
				cout << "\t\t\t\t" << m << endl;
			}
			

			createEdgeNodePredecessor(unitigGraph, unitigIndex1, predecessors[i], edgeNodes);
			createEdgeNodeSuccessor(unitigGraph, unitigIndex1, successors[j], edgeNodes);

			supportedTriplets.push_back({predecessors[i], unitigIndex1, successors[j]});


		}
	}

	//cout << "yo a" << endl;

	for(const vector<UnitigType>& triplet : supportedTriplets){
		pair<UnitigType, UnitigType> unitigPairPredecessor = {triplet[0], triplet[1]};
		if(edgeNodes.find(unitigPairPredecessor) == edgeNodes.end()) continue;

		UnitigGraph2::UnitigNode* edgeNodePredecessor = edgeNodes[unitigPairPredecessor];

		pair<UnitigType, UnitigType> unitigPairSuccessor = {triplet[1], triplet[2]};
		if(edgeNodes.find(unitigPairSuccessor) == edgeNodes.end()) continue;

		UnitigGraph2::UnitigNode* edgeNodeSuccessor = edgeNodes[unitigPairSuccessor];

		const UnitigType unitigIndexPredecessor = UnitigGraph2::unitigName_to_unitigIndex(edgeNodePredecessor->_unitigName, false);
		const UnitigType unitigIndexSuccessor = UnitigGraph2::unitigName_to_unitigIndex(edgeNodeSuccessor->_unitigName, false);

		unitigGraph->addSuccessor(unitigIndexPredecessor, unitigIndexSuccessor);
		unitigGraph->addSuccessor(UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndexSuccessor), UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndexPredecessor));

	}

	//cout << "yo b" << endl;
	//if(unitigName == 189){
	//	cout << "lulul" << endl;
	//	getchar();
	//}

}

*/


void CreateMdbg::solveSmallUnitigsSub2(UnitigGraph2* unitigGraph, UnitigGraph2::UnitigNode* node){
	
	vector<UnitigType> supportedPredecessors;
	vector<UnitigType> supportedSuccessors;
	phmap::parallel_flat_hash_map<pair<UnitigType, UnitigType>, UnitigGraph2::UnitigNode*> edgeNodes;

	//cout << node->_unitigName << endl;

	const UnitigType& unitigIndex = UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, false);

	vector<UnitigType> predecessors;
	unitigGraph->getPredecessors(unitigIndex, predecessors);

	vector<UnitigType> successors;
	unitigGraph->getSuccessors(unitigIndex, successors);


	const vector<MinimizerType> minimizers = _unitigName_to_minimizers[node->_unitigName];
	const UnitigType unitigName = node->_unitigName;

	//cout << endl << endl;
	//cout << "Remove node: utg" << node->_unitigName << endl;

	//for(MinimizerType m : _unitigName_to_minimizers[unitigName]){
	//	cout << "\t" << m << endl;
	//}

	//cout << "\tutg" << UnitigGraph2::unitigIndexToString(unitigIndex1) << endl;
	//for(MinimizerType m : minimizers){
	//	cout << "\t\t" << m << endl;
	//}

	for(size_t i=0; i<predecessors.size(); i++){

		if(predecessors[i] == unitigIndex){
			//cout << "Skip circular triplet" << endl;
			continue;
		}

		UnitigGraph2::UnitigNode* predecessor = unitigGraph->_unitigs[UnitigGraph2::unitigIndex_to_unitigName(predecessors[i])];
		vector<MinimizerType> predMinimizers = UnitigGraph2::getMinimizers(predecessors[i], _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(predecessors[i])]);

		int predecessorOverlapSize = longestOverlap(predecessors[i], predMinimizers, predecessor->_isEdgeNode, unitigIndex, minimizers, node->_isEdgeNode);
		//cout << predecessorOverlapSize << endl;

		vector<MinimizerType> triplet;
		
		triplet.push_back(predMinimizers[predMinimizers.size()-predecessorOverlapSize-1]);
		for(const MinimizerType m : minimizers){
			triplet.push_back(m);
		}

		//cout << triplet.size() << endl;

		KmerVec kminmer;
		kminmer._kmers = triplet;
		kminmer = kminmer.normalize();

		if(_mdbgNodesLight.find(kminmer.hash128()) != _mdbgNodesLight.end()){
			supportedPredecessors.push_back(predecessors[i]);
		}

	}


	for(size_t i=0; i<successors.size(); i++){

		if(successors[i] == unitigIndex){
			//cout << "Skip circular triplet" << endl;
			continue;
		}

		UnitigGraph2::UnitigNode* successor = unitigGraph->_unitigs[UnitigGraph2::unitigIndex_to_unitigName(successors[i])];
		vector<MinimizerType> succMinimizers = UnitigGraph2::getMinimizers(successors[i], _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(successors[i])]);
		
		int successorOverlapSize = longestOverlap(unitigIndex, minimizers, node->_isEdgeNode, successors[i], succMinimizers, successor->_isEdgeNode);

		vector<MinimizerType> triplet;
		
		for(const MinimizerType m : minimizers){
			triplet.push_back(m);
		}
		triplet.push_back(succMinimizers[successorOverlapSize]);

		//cout << triplet.size() << endl;

		KmerVec kminmer;
		kminmer._kmers = triplet;
		kminmer = kminmer.normalize();

		if(_mdbgNodesLight.find(kminmer.hash128()) != _mdbgNodesLight.end()){
			supportedSuccessors.push_back(successors[i]);
		}

	}

	/*
		for(size_t j=0; j<successors.size(); j++){

			//cout << "\ts " << j << " " << successors.size() << endl;

			if(predecessors[i] == successors[j]){
				cout << "Skip circular triplet" << endl;
				continue;
			}
			if(predecessors[i] == unitigIndex1){
				cout << "Skip self loop 1" << endl;
				continue;
			}
			if(successors[j] == unitigIndex1){
				cout << "Skip self loop 2" << endl;
				continue;
			}

			vector<MinimizerType> succMinimizers = UnitigGraph2::getMinimizers(successors[j], _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(successors[j])]);
		


			int successorOverlapSize = longestOverlap(minimizers, succMinimizers);

			//cout << predMinimizers.size() << " " << minimizers.size() << " " << succMinimizers.size() << " " << predecessorOverlapSize << " " << successorOverlapSize << endl;
			//if(predecessorOverlapSize != successorOverlapSize){
				//cout << "Overlap predecessor != overlap successor: " << predecessorOverlapSize << " " << successorOverlapSize << endl;
			//}
			
			//vector<MinimizerType> startNode = UnitigGraph2::getStartNode(successors[j], _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(successors[j])], _kminmerSizePrev);
			//cout << "\t\tutg" << UnitigGraph2::unitigIndexToString(successors[j]) << endl;
			//for(MinimizerType m : succMinimizers){
			//	cout << "\t\t\t\t" << m << endl;
			//}


			vector<MinimizerType> triplet;
			//triplet.push_back(endNode[0]);
			triplet.push_back(predMinimizers[predMinimizers.size()-predecessorOverlapSize-1]);
			for(const MinimizerType m : minimizers){
				triplet.push_back(m);
			}
			triplet.push_back(succMinimizers[successorOverlapSize]);
			//triplet.push_back(startNode[startNode.size()-1]);

			


			KmerVec kminmer;
			kminmer._kmers = triplet;
			kminmer = kminmer.normalize();

			if(_solidTriplets.find(kminmer.hash128()) == _solidTriplets.end()) continue;
			if(isUnsupportedUnitig(triplet)) continue;

			

			(unitigGraph, unitigIndex1, successors[j], edgeNodes);

			supportedTriplets.push_back({predecessors[i], unitigIndex1, successors[j]});


		}
	}
	*/

	for(const UnitigType unitigIndexPred : supportedPredecessors){
		createEdgeNodePredecessor(unitigGraph, unitigIndex, unitigIndexPred, edgeNodes);
	}

	for(const UnitigType unitigIndexSucc : supportedSuccessors){
		createEdgeNodeSuccessor(unitigGraph, unitigIndex, unitigIndexSucc, edgeNodes);
	}

	
	for(const UnitigType unitigIndexPred : supportedPredecessors){
		
		pair<UnitigType, UnitigType> unitigPairPredecessor = {unitigIndexPred, unitigIndex};
		if(edgeNodes.find(unitigPairPredecessor) == edgeNodes.end()) continue;
		
		UnitigGraph2::UnitigNode* edgeNodePredecessor = edgeNodes[unitigPairPredecessor];

		for(const UnitigType unitigIndexSucc : supportedSuccessors){

			pair<UnitigType, UnitigType> unitigPairSuccessor = {unitigIndex, unitigIndexSucc};
			if(edgeNodes.find(unitigPairSuccessor) == edgeNodes.end()) continue;

			UnitigGraph2::UnitigNode* edgeNodeSuccessor = edgeNodes[unitigPairSuccessor];

			const UnitigType unitigIndexPredecessor = UnitigGraph2::unitigName_to_unitigIndex(edgeNodePredecessor->_unitigName, false);
			const UnitigType unitigIndexSuccessor = UnitigGraph2::unitigName_to_unitigIndex(edgeNodeSuccessor->_unitigName, false);

			//vector<MinimizerType> minimizers = UnitigGraph2::getMinimizers(unitigIndexPredecessor, _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(unitigIndexPredecessor)]);
			//vector<MinimizerType> succMinimizers = UnitigGraph2::getMinimizers(unitigIndexSuccessor, _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(unitigIndexSuccessor)]);
		

			//if(!isEdgeSupported(unitigIndexPredecessor, minimizers, unitigIndexSuccessor, succMinimizers, _kminmerSizePrev, _kminmerAbundances)) continue;

			unitigGraph->addSuccessor(unitigIndexPredecessor, unitigIndexSuccessor);
			unitigGraph->addSuccessor(UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndexSuccessor), UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndexPredecessor));

		}
	}


	unitigGraph->removeNode(node);

	//cout << "yo b" << endl;
	//if(unitigName == 189){
	//	cout << "lulul" << endl;
	//	getchar();
	//}

	/*
	const UnitigType& unitigIndex1 = UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, false);

	vector<UnitigType> predecessors;
	unitigGraph->getPredecessors(unitigIndex1, predecessors);

	vector<UnitigType> successors;
	unitigGraph->getSuccessors(unitigIndex1, successors);
	

	//for(size_t i=0; i<predecessors.size(); i++){
		//if(unitigGraph->_unitigs[UnitigGraph2::unitigIndex_to_unitigName(predecessors[i])]->_nbMinimizers == _kminmerSizePrev) return;
	//}
		
	//for(size_t i=0; i<successors.size(); i++){
		//if(unitigGraph->_unitigs[UnitigGraph2::unitigIndex_to_unitigName(successors[i])]->_nbMinimizers == _kminmerSizePrev) return;
	//}


	if(node->_unitigName == 44 || node->_unitigName == 25 || node->_unitigName == 56){

		//if(node->_unitigName == 25) return;
	}
	else{
		cout << "Remove node: " << node->_unitigName << " " << successors.size() << " " << predecessors.size() << endl;

		unitigGraph->removeNode(node);
	}


	for(size_t i=0; i<predecessors.size(); i++){

		if(isUnitigNameToRemove.find(UnitigGraph2::unitigIndex_to_unitigName(predecessors[i])) != isUnitigNameToRemove.end()) continue;
		//if(node->_unitigName == 56){
		//	cout << "lul: " << UnitigGraph2::unitigIndex_to_unitigName(predecessors[i]) << " " << unitigGraph->_unitigs[UnitigGraph2::unitigIndex_to_unitigName(predecessors[i])]->_nbMinimizers << endl;
		//}

		//if(unitigGraph->_unitigs[UnitigGraph2::unitigIndex_to_unitigName(predecessors[i])]->_nbMinimizers == _kminmerSizePrev) continue;

		for(size_t j=0; j<successors.size(); j++){

			//if(unitigGraph->_unitigs[UnitigGraph2::unitigIndex_to_unitigName(successors[j])]->_nbMinimizers == _kminmerSizePrev) continue;
			if(isUnitigNameToRemove.find(UnitigGraph2::unitigIndex_to_unitigName(successors[j])) != isUnitigNameToRemove.end()) continue;

			unitigGraph->addSuccessor(predecessors[i], successors[j]);
			unitigGraph->addSuccessor(UnitigGraph2::unitigIndex_to_reverseDirection(successors[j]), UnitigGraph2::unitigIndex_to_reverseDirection(predecessors[i]));

			//addEdgeIfSupported(unitigGraph, predecessors[i], successors[j]);
			//addEdgeIfSupported(unitigGraph, UnitigGraph2::unitigIndex_to_reverseDirection(successors[j]), UnitigGraph2::unitigIndex_to_reverseDirection(predecessors[i]));

		}
	}
	*/

}

void CreateMdbg::createEdgeNodePredecessor(UnitigGraph2* unitigGraph, const UnitigType unitigIndex, const UnitigType unitigIndexPredecessor, phmap::parallel_flat_hash_map<pair<UnitigType, UnitigType>, UnitigGraph2::UnitigNode*>& edgeNodes){

	pair<UnitigType, UnitigType> unitigPair = {unitigIndexPredecessor, unitigIndex};

	if(edgeNodes.find(unitigPair) != edgeNodes.end()){
		return;// edgeNodes[unitigPair];
	}

	UnitigGraph2::UnitigNode* node = unitigGraph->_unitigs[UnitigGraph2::unitigIndex_to_unitigName(unitigIndex)];
	//UnitigGraph2::UnitigNode* node = unitigGraph->_unitigs[UnitigGraph2::unitigIndex_to_unitigName(unitigIndex)];
	UnitigGraph2::UnitigNode* predecessor = unitigGraph->_unitigs[UnitigGraph2::unitigIndex_to_unitigName(unitigIndexPredecessor)];
	vector<MinimizerType> predMinimizers = UnitigGraph2::getMinimizers(unitigIndexPredecessor, _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(unitigIndexPredecessor)]);
	vector<MinimizerType> minimizers = UnitigGraph2::getMinimizers(unitigIndex, _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(unitigIndex)]);

	//vector<MinimizerType> endNode = UnitigGraph2::getEndNode(unitigIndexPredecessor, _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(unitigIndexPredecessor)], _kminmerSizePrev);

	int overlapSize = longestOverlap(unitigIndexPredecessor, predMinimizers, predecessor->_isEdgeNode, unitigIndex, minimizers, node->_isEdgeNode);

	//bool connectTwoSmallNodes = predecessor->_nbMinimizers == _kminmerSizePrev;

	//vector<MinimizerType> minimizers = _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(unitigIndex)];
	//if(!connectTwoSmallNodes){
	minimizers.insert(minimizers.begin(), predMinimizers[predMinimizers.size()-overlapSize-1]);
	//}

	//if(isUnsupportedUnitig(minimizers)) return;

	UnitigGraph2::UnitigNode* edgeNode = createEdgeNode(unitigGraph, minimizers);
	edgeNodes[unitigPair] = edgeNode;

	const UnitigType edgeNodeUnitigIndex = UnitigGraph2::unitigName_to_unitigIndex(edgeNode->_unitigName, false);
	unitigGraph->addSuccessor(unitigIndexPredecessor, edgeNodeUnitigIndex);
	unitigGraph->addSuccessor(UnitigGraph2::unitigIndex_to_reverseDirection(edgeNodeUnitigIndex), UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndexPredecessor));

	//for(MinimizerType m : minimizers){
	//	cout << "\t\t" << m << endl;
	//}

	//cout << "\tCreate edge node P: " << UnitigGraph2::unitigIndexToString(unitigIndexPredecessor) << " -> " << UnitigGraph2::unitigIndexToString(unitigIndex) << "    utg" << edgeNode->_unitigName << endl;
	
	//vector<UnitigType> predecessors;
	//unitigGraph->getPredecessors(edgeNodeUnitigIndex, predecessors);

	//for(UnitigType unitigIndex : predecessors){
	//	cout << "\t\tutg" << UnitigGraph2::unitigIndex_to_unitigName(unitigIndex) << endl;
	//	for(MinimizerType m : UnitigGraph2::getMinimizers(unitigIndex, _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(unitigIndex)])){
	//		cout << "\t\t\t" << m << endl;
	//	}
	//}

	/*
	if(edgeNode->_unitigName == 195){
		cout << "lul" << endl;
		for(MinimizerType m : _unitigName_to_minimizers[edgeNode->_unitigName]){
			cout << "\t\t\t" << m << endl;
		}
		getchar();
	}
	*/

}

void CreateMdbg::createEdgeNodeSuccessor(UnitigGraph2* unitigGraph, const UnitigType unitigIndex, const UnitigType unitigIndexSuccessor, phmap::parallel_flat_hash_map<pair<UnitigType, UnitigType>, UnitigGraph2::UnitigNode*>& edgeNodes){


	pair<UnitigType, UnitigType> unitigPair = {unitigIndex, unitigIndexSuccessor};

	if(edgeNodes.find(unitigPair) != edgeNodes.end()){
		return;// edgeNodes[unitigPair];
	}

	//UnitigGraph2::UnitigNode* node = unitigGraph->_unitigs[UnitigGraph2::unitigIndex_to_unitigName(unitigIndex)];
	UnitigGraph2::UnitigNode* node = unitigGraph->_unitigs[UnitigGraph2::unitigIndex_to_unitigName(unitigIndex)];
	UnitigGraph2::UnitigNode* successor = unitigGraph->_unitigs[UnitigGraph2::unitigIndex_to_unitigName(unitigIndexSuccessor)];
	vector<MinimizerType> succMinimizers = UnitigGraph2::getMinimizers(unitigIndexSuccessor, _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(unitigIndexSuccessor)]);
	vector<MinimizerType> minimizers = UnitigGraph2::getMinimizers(unitigIndex, _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(unitigIndex)]);

	int overlapSize = longestOverlap(unitigIndex, minimizers, node->_isEdgeNode, unitigIndexSuccessor, succMinimizers, successor->_isEdgeNode);

	//bool connectTwoSmallNodes = successor->_nbMinimizers == _kminmerSizePrev;

	//vector<MinimizerType> minimizers = _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(unitigIndex)];
	//if(!connectTwoSmallNodes){
	minimizers.push_back(succMinimizers[overlapSize]);
	//}

	//if(isUnsupportedUnitig(minimizers)) return;

	UnitigGraph2::UnitigNode* edgeNode = createEdgeNode(unitigGraph, minimizers);
	//UnitigGraph2::UnitigNode* edgeNode = unitigGraph->addNode(unitigGraph->_unitigs.size(), minimizers);
	edgeNodes[unitigPair] = edgeNode;
	//edgeNode->_abundance = 40;
	//_unitigName_to_minimizers[edgeNode->_unitigName] = minimizers;

	const UnitigType edgeNodeUnitigIndex = UnitigGraph2::unitigName_to_unitigIndex(edgeNode->_unitigName, false);
	unitigGraph->addSuccessor(edgeNodeUnitigIndex, unitigIndexSuccessor);
	unitigGraph->addSuccessor(UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndexSuccessor), UnitigGraph2::unitigIndex_to_reverseDirection(edgeNodeUnitigIndex));

	//cout << "\tCreate edge node S: " << UnitigGraph2::unitigIndexToString(unitigIndex) << " -> " << UnitigGraph2::unitigIndexToString(unitigIndexSuccessor) << "    utg" << edgeNode->_unitigName << endl;
	//for(MinimizerType m : minimizers){
	//	cout << "\t\t" << m << endl;
	//}


	//vector<UnitigType> successors;
	//unitigGraph->getSuccessors(edgeNodeUnitigIndex, successors);

	//for(UnitigType unitigIndex : successors){
	//	cout << "\t\tutg" << UnitigGraph2::unitigIndex_to_unitigName(unitigIndex) << endl;
	//	for(MinimizerType m : UnitigGraph2::getMinimizers(unitigIndex, _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(unitigIndex)])){
	//		cout << "\t\t\t" << m << endl;
	//	}
	//}

	/*
	if(edgeNode->_unitigName == 195){
		cout << "lul" << endl;
		for(MinimizerType m : _unitigName_to_minimizers[edgeNode->_unitigName]){
			cout << "\t\t\t" << m << endl;
		}
		getchar();
	}
	*/

}

/*
bool CreateMdbg::isUnsupportedUnitig(const vector<MinimizerType>& minimizers){
	
	//cout << _kminmerSize << " " << minimizers.size() << endl;
	

	//cout << minimizers.size() << endl;
	//const vector<MinimizerType>& minimizers = _unitigName_to_minimizers[node->_unitigName];

	for(const KmerVec& kminmer : MDBG::minimizersToKminmers(minimizers, _kminmerSize)){

		const u_int128_t hash = kminmer.normalize().hash128();
		
		if(_mdbgNodesLight.find(hash) != _mdbgNodesLight.end()){
			if(_kminmerSize == _kminmerSizeFirst+1){
				if (_mdbgNodesLight[hash]._abundance > 10 && _mdbgNodesLight[hash]._index <= 1) {
					cout << "Found chimeric kminmer" << endl;
					return true;
				}
			}
		}
		else{
			cout << "Unitig has kminmers not in set" << endl;
			return true;
		}
	}

	return false;

}
*/
UnitigGraph2::UnitigNode* CreateMdbg::createEdgeNode(UnitigGraph2* unitigGraph, const vector<MinimizerType>& minimizers){

	UnitigGraph2::UnitigNode* edgeNode = unitigGraph->addNode(unitigGraph->_unitigs.size(), minimizers.size(), true);
	//edgeNode->_abundance = 40;
	//_unitigName_to_minimizers[edgeNode->_unitigName] = minimizers;
	_unitigName_to_minimizers.push_back(minimizers);

	u_int32_t minAbundance = -1;
	/*
	for(const KmerVec& kminmer : Utils::minimizersToKminmers(minimizers, _kminmerSizePrev)){

		const u_int128_t hash = kminmer.normalize().hash128();



		//cout << minimizers.size() << " " << kminmer._kmers.size() << " " << (_kminmerAbundances.find(hash) != _kminmerAbundances.end()) << endl;

		if(_kminmerAbundances.find(hash) != _kminmerAbundances.end()){

			const u_int32_t abundance = _kminmerAbundances[hash];
			edgeNode->_abundances.push_back(abundance);
			
			if(abundance < minAbundance){
				minAbundance = abundance;
			}

		}
		else{
			edgeNode->_abundances.push_back(1);
			minAbundance = 1;
		}


	}
	*/
	/*
	for(const KmerVec& kminmer : MDBG::minimizersToKminmers(minimizers, _kminmerSize)){

		const u_int128_t hash = kminmer.normalize().hash128();



		//cout << minimizers.size() << " " << kminmer._kmers.size() << " " << (_kminmerAbundances.find(hash) != _kminmerAbundances.end()) << endl;

		if(_mdbgNodesLight.find(hash) != _mdbgNodesLight.end()){

			const u_int32_t abundance = _mdbgNodesLight[hash]._abundance;
			edgeNode->_abundances.push_back(abundance);
			
			if(abundance < minAbundance){
				minAbundance = abundance;
			}

		}
		else{
			edgeNode->_abundances.push_back(1);
			minAbundance = 1;
		}


	}
	*/
	/*
	bool isError = false;

	for(const KmerVec& kminmer : Utils::minimizersToKminmers(minimizers, _debug_kminmerSize)){

		const u_int128_t hash = kminmer.normalize().hash128();

		if(_debug_vecHash_to_referencePosition.find(hash) == _debug_vecHash_to_referencePosition.end()){
			isError = true;
		}
	}

	if(isError){
		cout << "Erroneous: utg" <<  edgeNode->_unitigName << endl;
		for(const KmerVec& kminmer : Utils::minimizersToKminmers(minimizers, _kminmerSizePrev)){

			const u_int128_t hash = kminmer.normalize().hash128();

			if(_kminmerAbundances.find(hash) != _kminmerAbundances.end()){
				cout << "\t" << kminmer.toString() << " " << _kminmerAbundances[hash] << endl;
			}
			else{
				cout << "\t" << kminmer.toString() << " 1" << endl;
			}

			
		}

	}
	*/
	if(minimizers.size() < _kminmerSize){
		cout << "derp" << endl;
		getchar();
	}


	KmerVec kminmer;
	kminmer._kmers = minimizers;
	const u_int128_t hash = kminmer.normalize().hash128();

	if(_mdbgNodesLight.find(hash) != _mdbgNodesLight.end()){

		const u_int32_t abundance = _mdbgNodesLight[hash];
		edgeNode->_abundances.push_back(abundance);

	}
	else{
		edgeNode->_abundances.push_back(1);
	}

	//for(size_t i=0; i<edgeNode->_abundances.size(); i++){
	//	edgeNode->_abundances[i] = minAbundance;
	//}

	std::sort(edgeNode->_abundances.begin(), edgeNode->_abundances.end());
	edgeNode->_abundance = edgeNode->computeMedianAbundance();

	//cout << "utg" << edgeNode->_unitigName << endl;
	//for(MinimizerType m : minimizers){
	//	cout << m << endl;
	//}

	//if(edgeNode->_unitigName == 311){
		//cout << "derp" << endl;
		//getchar();
	//}
	
	//cout << "Add edge node: " << edgeNode->_unitigName << endl;
	//for(size_t i=0; i<minimizers.size(); i++){
	//	cout << "\t" << minimizers[i] << endl;
	//}

	return edgeNode;
}


void CreateMdbg::writeUnitigs(UnitigGraph2* unitigGraph){
	/*
	cout << "lalala:" << endl;
	cout << endl;
	cout << "utg160" << endl; 
	for(MinimizerType m : _unitigName_to_minimizers[160]){
		cout << m << endl;
	}
	cout << "utg175" << endl; 
	for(MinimizerType m : _unitigName_to_minimizers[175]){
		cout << m << endl;
	}
	cout << "utg156" << endl; 
	for(MinimizerType m : _unitigName_to_minimizers[156]){
		cout << m << endl;
	}

	cout << "utg160" << endl; 
	for(size_t i=0; i<unitigGraph->_unitigs[160]->_unitigMerge.size(); i++){
		cout << UnitigGraph2::unitigIndexToString(unitigGraph->_unitigs[160]->_unitigMerge[i]) << endl;
	}
	
	
	cout << "utg156" << endl; 
	for(size_t i=0; i<unitigGraph->_unitigs[156]->_unitigMerge.size(); i++){
		cout << UnitigGraph2::unitigIndexToString(unitigGraph->_unitigs[156]->_unitigMerge[i]) << endl;
	}
	getchar();
	cout << "lalala:" << endl;
	*/

	//computeUnitigSequences(unitigGraph, "");


	/*
	for(size_t i=0; i<unitigGraph->_unitigs.size(); i++){

		UnitigGraph2::UnitigNode* node = unitigGraph->_unitigs[i];

		if(node == nullptr) continue;

		extendNode(unitigGraph, UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, false), false);
		extendNode(unitigGraph, UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, true), false);
	}



	for(size_t i=0; i<unitigGraph->_unitigs.size(); i++){

		UnitigGraph2::UnitigNode* node = unitigGraph->_unitigs[i];

		if(node == nullptr) continue;

		extendNodeBranching(unitigGraph, UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, false), false);
		extendNodeBranching(unitigGraph, UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, true), false);
	}


	cout << endl << endl;

	for(size_t i=0; i<unitigGraph->_unitigs.size(); i++){

		UnitigGraph2::UnitigNode* node = unitigGraph->_unitigs[i];

		if(node == nullptr) continue;

		extendNode(unitigGraph, UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, false), true);
		extendNode(unitigGraph, UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, true), true);
	}

	*/

	//loadReference(_kminmerSize);
	//unitigGraph->save("/pasteur/appa/homes/gbenoit/zeus/tmp/unitigs/assembly_graph_2.gfa", _outputDir);
	//writeReferencePositions(unitigGraph, "/pasteur/appa/homes/gbenoit/zeus/tmp/unitigs/assembly_graph_2.gfa.referencePos.csv");


	/*

	for(const UnitigGraph2::UnitigNode* unitig : unitigGraph->_unitigs){
		
		if(unitig == nullptr) continue;

		UnitigType unitigIndex = UnitigGraph2::unitigName_to_unitigIndex(unitig->_unitigName, false);

		const vector<MinimizerType>& minimizers = _unitigName_to_minimizers[unitig->_unitigName];
		const vector<KmerVec>& kminmers = Utils::minimizersToKminmers(minimizers, _kminmerSize);

		vector<float> abundances;

		for(const KmerVec& kminmer : kminmers){

			const u_int128_t vecHash = kminmer.normalize().hash128();

			if(unitig->_unitigName == 38){
				cout << kminmer.normalize().toString() << " " << (_mdbgNodesLight.find(vecHash) != _mdbgNodesLight.end()) << endl;
			}
			if(unitig->_unitigName == 337){
				cout << kminmer.normalize().toString() << " " << (_mdbgNodesLight.find(vecHash) != _mdbgNodesLight.end()) << endl;
			}

			//cout << "lul: " << _mdbgNodesLight[vecHash]._abundance << endl;
			if(_mdbgNodesLight.find(vecHash) == _mdbgNodesLight.end()){
				abundances.push_back(1);
			}
			else{
				abundances.push_back(_mdbgNodesLight[vecHash]._abundance);
			}
		}

	}
	*/


	for(UnitigGraph2::UnitigNode* node : unitigGraph->_unitigs){

		if(node == nullptr) continue;

		unitigGraph->recompact(node);
	}

	for(UnitigGraph2::UnitigNode* node : unitigGraph->_unitigs){

		if(node == nullptr) continue;

		if(node->_unitigMerge.size() == 0 && node->_nbMinimizers == _kminmerSizePrev){
			cout << "Removing remaining small unitig: utg" << node->_unitigName << endl;
			unitigGraph->removeNode(node);
		}
	}



	/*
	vector<UnitigType> successors;
	unitigGraph->getSuccessors(UnitigGraph2::unitigName_to_unitigIndex(428, false), successors);

	vector<UnitigType> predecessors;
	unitigGraph->getPredecessors(UnitigGraph2::unitigName_to_unitigIndex(428, false), predecessors);


	for(size_t i=0; i<predecessors.size(); i++){
		cout << "\tutg" << UnitigGraph2::unitigIndexToString(predecessors[i]) << endl;
		for(MinimizerType m : UnitigGraph2::getMinimizers(predecessors[i], _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(predecessors[i])])){
			cout << "\t\t" << m << endl;
		}

	}

	cout << "\tutg" << UnitigGraph2::unitigIndexToString(UnitigGraph2::unitigName_to_unitigIndex(428, false)) << endl;
	for(MinimizerType m : UnitigGraph2::getMinimizers(UnitigGraph2::unitigName_to_unitigIndex(428, false), _unitigName_to_minimizers[428])){
		cout << "\t\t" << m << endl;
	}

	for(size_t i=0; i<successors.size(); i++){
		cout << "\tutg" << UnitigGraph2::unitigIndexToString(successors[i]) << endl;
		for(MinimizerType m : UnitigGraph2::getMinimizers(successors[i], _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(successors[i])])){
			cout << "\t\t" << m << endl;
		}

	}
	*/


	//loadReference(_kminmerSize);
	//unitigGraph->save("/pasteur/appa/homes/gbenoit/zeus/tmp/unitigs/assembly_graph_3.gfa", _outputDir);
	//writeReferencePositions(unitigGraph, "/pasteur/appa/homes/gbenoit/zeus/tmp/unitigs/assembly_graph_3.gfa.referencePos.csv");

	/*
	for(auto unitigIndex : unitigGraph->_unitigs[311]->_unitigMerge){
		cout << UnitigGraph2::unitigIndexToString(unitigIndex) << endl;
	}

	for(UnitigGraph2::UnitigNode* node : unitigGraph->_unitigs){

		if(node == nullptr) continue;

		if(node->_unitigName == 326){
			cout << "utg" << node->_unitigName << endl;
			for(size_t i=0; i<node->_abundances.size(); i++){
				cout << "\t" << node->_abundances[i] << endl;
			}
		}
		
	}

	exit(1);
	*/








	//exit(1);
	//cout << "oploplop" << endl;
	//if(_kminmerSize == 5) exit(1);
	//exit(1);
	
	/*
	cout << "utg195" << endl; 
	for(MinimizerType m : _unitigName_to_minimizers[195]){
		cout << m << endl;
	}

	cout << "utg239" << endl; 
	for(MinimizerType m : _unitigName_to_minimizers[239]){
		cout << m << endl;
	}


	for(size_t i=0; i<unitigGraph->_unitigs.size(); i++){

		if(unitigGraph->_unitigs[i] == nullptr) continue;

		if(_unitigName_to_minimizers[unitigGraph->_unitigs[i]->_unitigName].size() == 4){
			cout << "wtf: utg" << unitigGraph->_unitigs[i]->_unitigName << endl;
			getchar();
		}
	}
	*/

	//if(_kminmerSize == 6) exit(1);

	//computeUnitigSequences(unitigGraph, _outputDir + "/unitigGraph.nodes.bin");

	UnitigType unitigIndexNew = 0;
	vector<UnitigType> currentUnitigName_to_newUnitigName(unitigGraph->_unitigs.size(), 0);

	for(size_t i=0; i<unitigGraph->_unitigs.size(); i++){

		UnitigGraph2::UnitigNode* node = unitigGraph->_unitigs[i];

		if(node == nullptr) continue;

		currentUnitigName_to_newUnitigName[node->_unitigName] = unitigIndexNew;

		unitigIndexNew += 1;
	}



	_nbUnitigNodes = 0;

	Logger::get().debug() << "Write unitig nodes";
	auto start = high_resolution_clock::now();
	writeUnitigSequences(unitigGraph, currentUnitigName_to_newUnitigName);
	Logger::get().debug() << "\tDuration: " << duration_cast<seconds>(high_resolution_clock::now() - start).count();






	/*
	u_int32_t nbSolidUnitigs = 0;
	u_int32_t nbUnitigsToRemove = 0;
	u_int32_t nbUnitigsHybrid = 0;

	vector<UnitigGraph2::UnitigNode*> nodesToRemove;

	for(UnitigGraph2::UnitigNode* node : unitigGraph->_unitigs){

		if(node == nullptr) continue;

		const vector<MinimizerType>& minimizers = _unitigName_to_minimizers[node->_unitigName];

		int nbSolidKminmers = 0;

		const vector<KmerVec>& kminmers = MDBG::minimizersToKminmers(minimizers, _kminmerSize);

		for(const KmerVec& kminmer : kminmers){

			const u_int128_t hash = kminmer.normalize().hash128();
			
			if(_mdbgNodesLight.find(hash) != _mdbgNodesLight.end()){
				nbSolidKminmers += 1;
			}
		}

		if(nbSolidKminmers == kminmers.size()){
			nbSolidUnitigs += 1;
		}
		else if(nbSolidKminmers > 0){
			nbUnitigsHybrid += 1;
		}
		else{
			nbUnitigsToRemove += 1;
		}


	}

	cout << "Unitig stats: " << nbSolidUnitigs << " " << nbUnitigsToRemove << " " << nbUnitigsHybrid << endl;
	*/


	Logger::get().debug() << "Write unitig edges";
	start = high_resolution_clock::now();

	int nbEdgesWritten = 0;
	_nbUnitigEdges = 0;

	ofstream unitigGraphFile_edges(_outputDir + "/unitigGraph.edges.successors.bin");

	unordered_set<UnitigType> isVisited;

	for(const UnitigGraph2::UnitigNode* node : unitigGraph->_unitigs){
		
		if(node == nullptr) continue;

		UnitigType sourceUnitigIndex = UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, false);

		if(isVisited.find(node->_unitigName) != isVisited.end()) continue;

		queue<UnitigType> queue;
		queue.push(sourceUnitigIndex);

		while(queue.size() > 0){

			UnitigType unitigIndex = queue.front();
			queue.pop();

			bool isReversed;
			UnitigType unitigName = UnitigGraph2::unitigIndex_to_unitigName(unitigIndex, isReversed);
			

			if(isVisited.find(unitigName) != isVisited.end()) continue;
			isVisited.insert(unitigName);

			//UnitigGraph2::UnitigNode* unitig = _unitigs[unitigName];

			vector<UnitigType> successors;
			unitigGraph->getSuccessors(unitigIndex, successors);

			vector<UnitigType> predecessors;
			unitigGraph->getPredecessors(unitigIndex, predecessors);



/*


	UnitigType unitigIndexFrom_rev = Utils::unitigIndexToReverseDirection(unitigIndexFrom);


		u_int32_t nbSuccessors = successors.size();
		u_int32_t nbPredecessors = predecessors.size();


		//_unitigGraphFile_edges_successors.write((const char*)&fromUnitigIndex, sizeof(fromUnitigIndex));
		//_unitigGraphFile_edges_successors.write((const char*)&toUnitigIndex, sizeof(toUnitigIndex));

		_unitigGraphFile_edges_successors.write((const char*)&unitigIndexFrom, sizeof(unitigIndexFrom));
		_unitigGraphFile_edges_successors.write((const char*)&nbSuccessors, sizeof(nbSuccessors));
		_unitigGraphFile_edges_successors.write((const char*)&successors[0], successors.size() * sizeof(UnitigType));
		_unitigGraphFile_edges_successors.write((const char*)&nbPredecessors, sizeof(nbPredecessors));
		_unitigGraphFile_edges_successors.write((const char*)&predecessors[0], predecessors.size() * sizeof(UnitigType));

		for(const UnitigType& unitigIndexTo : successors){
			_checksum_unitigEdges += unitigIndexFrom * unitigIndexTo;
		}
		for(const UnitigType& unitigIndexTo : predecessors){
			_checksum_unitigEdges += unitigIndexFrom_rev * unitigIndexTo;
		}
*/


			const UnitigType newUnitigName = currentUnitigName_to_newUnitigName[unitigName];
			unitigIndex = UnitigGraph2::unitigName_to_unitigIndex(newUnitigName, isReversed);

			//cout << "----" << endl;
			//cout << unitigIndex << endl;

			//UnitigType unitigIndexFrom;
			vector<UnitigType> successors2;

			for(UnitigType unitigIndexSuccessor : successors){
				
				//cout << "---" << endl;
				//bool isReversedSuccessor;
				//const UnitigType& unitigNameSuccessor = UnitigGraph2::unitigIndex_to_unitigName(unitigIndexSuccessor, isReversedSuccessor);
				
				bool isReversedSucc;
				UnitigType unitigNameSucc = UnitigGraph2::unitigIndex_to_unitigName(unitigIndexSuccessor, isReversedSucc);
				queue.push(unitigIndexSuccessor);

				const UnitigType newUnitigNameSucc = currentUnitigName_to_newUnitigName[unitigNameSucc];
				unitigIndexSuccessor = UnitigGraph2::unitigName_to_unitigIndex(newUnitigNameSucc, isReversedSucc);

				//cout << "utg" << unitigIndex << " -> " << "utg" << unitigIndexSuccessor << endl;
				//unitigGraphFile_edges.write((const char*)&unitigIndex, sizeof(unitigIndex));
				//unitigGraphFile_edges.write((const char*)&unitigIndexSuccessor, sizeof(unitigIndexSuccessor));

				nbEdgesWritten += 1;
				//UnitigType unitigIndexRev = UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndex);
				//UnitigType unitigIndexSuccessorRev = UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndexSuccessor);

				//cout << "utg" << unitigIndexSuccessorRev << " -> " << "utg" << unitigIndexRev << endl;

				//unitigGraphFile_edges.write((const char*)&unitigIndexSuccessorRev, sizeof(unitigIndexSuccessorRev));
				//unitigGraphFile_edges.write((const char*)&unitigIndexRev, sizeof(unitigIndexRev));
				//cout << "a " << unitigIndex << " " << unitigIndexSuccessor << endl;
				//unitigIndexFrom = unitigIndex;
				successors2.push_back(unitigIndexSuccessor);
			}

			vector<UnitigType> predecessors2;

			for(UnitigType unitigIndexPredecessor : predecessors){

				bool isReversedPred;
				UnitigType unitigNamePred = UnitigGraph2::unitigIndex_to_unitigName(unitigIndexPredecessor, isReversedPred);
				queue.push(unitigIndexPredecessor);

				const UnitigType newUnitigNamePred = currentUnitigName_to_newUnitigName[unitigNamePred];
				unitigIndexPredecessor = UnitigGraph2::unitigName_to_unitigIndex(newUnitigNamePred, isReversedPred);

				nbEdgesWritten += 1;

				
				//unitigGraphFile_edges.write((const char*)&unitigIndexPredecessor, sizeof(unitigIndexPredecessor));
				//unitigGraphFile_edges.write((const char*)&unitigIndex, sizeof(unitigIndex));
				
				UnitigType unitigIndexRev = UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndex);
				UnitigType unitigIndexPredecessorRev = UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndexPredecessor);

				//unitigGraphFile_edges.write((const char*)&unitigIndexRev, sizeof(unitigIndexRev));
				//unitigGraphFile_edges.write((const char*)&unitigIndexPredecessorRev, sizeof(unitigIndexPredecessorRev));


				//bool isReversedPredecessor;
				//const UnitigType& unitigNamePredecessor = UnitigGraph2::unitigIndex_to_unitigName(unitigIndexPredecessor, isReversedPredecessor);
				
				//cout << "b " << unitigIndexRev << " " << unitigIndexPredecessorRev << endl;
				
				predecessors2.push_back(unitigIndexPredecessorRev);
			}
			

			const u_int32_t nbSuccessors = successors2.size();
			const u_int32_t nbPredecessors = predecessors2.size();

			unitigGraphFile_edges.write((const char*)&unitigIndex, sizeof(unitigIndex));
			unitigGraphFile_edges.write((const char*)&nbSuccessors, sizeof(nbSuccessors));
			unitigGraphFile_edges.write((const char*)&successors2[0], successors2.size() * sizeof(UnitigType));
			unitigGraphFile_edges.write((const char*)&nbPredecessors, sizeof(nbPredecessors));
			unitigGraphFile_edges.write((const char*)&predecessors2[0], predecessors2.size() * sizeof(UnitigType));
			
			_nbUnitigEdges += nbSuccessors + nbPredecessors;
		}

	}

	unitigGraphFile_edges.close();
	Logger::get().debug() << "\tDuration: " << duration_cast<seconds>(high_resolution_clock::now() - start).count();





	Logger::get().debug() << "Write unitig abundances";
	start = high_resolution_clock::now();

	ofstream unitigGraphFile_nodes_abundances(_outputDir + "/unitigGraph.nodes.abundances.bin");

	/*
	#pragma omp parallel for num_threads(_nbCores)
	for(size_t i=0; i<unitigGraph->_unitigs.size(); i++){

		UnitigGraph2::UnitigNode* node = unitigGraph->_unitigs[i];
		if(node == nullptr) continue;
		
		UnitigType unitigName = currentUnitigName_to_newUnitigName[node->_unitigName];
		UnitigType unitigIndex = unitigName*2; //UnitigGraph2::unitigName_to_unitigIndex(unitigName, false);
		

		const vector<MinimizerType>& minimizers = _unitigName_to_minimizers[node->_unitigName];
		const vector<KmerVec>& kminmers = MDBG::minimizersToKminmers(minimizers, _kminmerSize);

		vector<u_int32_t> abundances;

		for(const KmerVec& kminmer : kminmers){

			const u_int128_t vecHash = kminmer.normalize().hash128();

			//if(unitig->_unitigName == 39){
			//	cout << kminmer.normalize().toString() << " " << (_mdbgNodesLight.find(vecHash) != _mdbgNodesLight.end()) << " " << (_debug_vecHash_to_referencePosition.find(vecHash) != _debug_vecHash_to_referencePosition.end()) << endl;
				//printf("%llx\n", (unsigned long long)(vecHash & 0xFFFFFFFFFFFFFFFF));
			//	printf("unsinged int 128 bit: %llu\n", (unsigned long long)vecHash);
			//}

			//if(unitig->_unitigName == 38){
			//	cout << kminmer.toString() << " " << (_mdbgNodesLight.find(vecHash) != _mdbgNodesLight.end()) << endl;
			//}
			//cout << "lul: " << _mdbgNodesLight[vecHash]._abundance << endl;
			if(_mdbgNodesLight.find(vecHash) == _mdbgNodesLight.end()){
				abundances.push_back(1);
			}
			else{
				abundances.push_back(_mdbgNodesLight[vecHash]);
			}
		}
		

		//const vector<float>& abundances = node->_abundances;


		#pragma omp critical(writeAbundances)
		{
			u_int32_t nbAbundances = abundances.size();
			unitigGraphFile_nodes_abundances.write((const char*)&unitigIndex, sizeof(unitigIndex));
			unitigGraphFile_nodes_abundances.write((const char*)&nbAbundances, sizeof(nbAbundances));
			unitigGraphFile_nodes_abundances.write((const char*)&abundances[0], nbAbundances * sizeof(u_int32_t));
			
			cout << "write abundance: " << unitigIndex << " " << abundances.size() << " " << minimizers.size() << endl;
			if(unitigIndex == 88) getchar();

		}
	}
	*/


	ifstream nodeFile(_outputDir + "/unitigGraph.nodes.bin");

	#pragma omp parallel num_threads(_nbCores)
	{

		bool isEOF = false;
		UnitigType unitigIndex;
		vector<MinimizerType> minimizers;
		//KmerVec vec;

		while (true) {

			#pragma omp critical(loadRefinedAbundances)
			{



				u_int32_t size;
				nodeFile.read((char*)&size, sizeof(size));
				
				isEOF = nodeFile.eof();
				//if(nodeFile.eof()) break;

				if(!isEOF){
					minimizers.resize(size);
					nodeFile.read((char*)&minimizers[0], size * sizeof(MinimizerType));


					nodeFile.read((char*)&unitigIndex, sizeof(unitigIndex));
				}

			}

			if(isEOF) break;
			
			//const UnitigType& unitigName = unitigIndex/2;

			const vector<KmerVec>& kminmers = MDBG::minimizersToKminmers(minimizers, _kminmerSize);

			vector<u_int32_t> abundances;

			for(const KmerVec& kminmer : kminmers){

				const u_int128_t vecHash = kminmer.normalize().hash128();

				//if(unitig->_unitigName == 39){
				//	cout << kminmer.normalize().toString() << " " << (_mdbgNodesLight.find(vecHash) != _mdbgNodesLight.end()) << " " << (_debug_vecHash_to_referencePosition.find(vecHash) != _debug_vecHash_to_referencePosition.end()) << endl;
					//printf("%llx\n", (unsigned long long)(vecHash & 0xFFFFFFFFFFFFFFFF));
				//	printf("unsinged int 128 bit: %llu\n", (unsigned long long)vecHash);
				//}

				//if(unitig->_unitigName == 38){
				//	cout << kminmer.toString() << " " << (_mdbgNodesLight.find(vecHash) != _mdbgNodesLight.end()) << endl;
				//}
				//cout << "lul: " << _mdbgNodesLight[vecHash]._abundance << endl;
				if(_mdbgNodesLight.find(vecHash) == _mdbgNodesLight.end()){
					abundances.push_back(1);
				}
				else{
					abundances.push_back(_mdbgNodesLight[vecHash]);
				}
			}
			

			//const vector<float>& abundances = node->_abundances;


			#pragma omp critical(writeAbundances)
			{
				u_int32_t nbAbundances = abundances.size();
				unitigGraphFile_nodes_abundances.write((const char*)&unitigIndex, sizeof(unitigIndex));
				unitigGraphFile_nodes_abundances.write((const char*)&nbAbundances, sizeof(nbAbundances));
				unitigGraphFile_nodes_abundances.write((const char*)&abundances[0], nbAbundances * sizeof(u_int32_t));
				
				//cout << "write abundance: " << unitigIndex << " " << abundances.size() << " " << minimizers.size() << endl;
				//if(unitigIndex == 88) getchar();

			}


		}


	}


	unitigGraphFile_nodes_abundances.close();

	
	writeUnitigGraphStat(_nbUnitigNodes, _nbUnitigEdges);
	
	Logger::get().debug() << "\tDuration: " << duration_cast<seconds>(high_resolution_clock::now() - start).count();


	//cout << "Nb written edges:   " << nbEdgesWritten << endl;

	/*
	Logger::get().debug() << "Loading unitig graph";
	UnitigGraph2* unitigGraphNew = new UnitigGraph2(_kminmerSizePrev, _kminmerLengthMean, _kminmerOverlapMean, _kminmerLengthMean-_kminmerOverlapMean, _minimizerSpacingMean, _kmerVec_to_nodeName, _nbCores);
	unitigGraphNew->load(_outputDir);

	//loadReference(_kminmerSize);
	//unitigGraphNew->save("/pasteur/appa/homes/gbenoit/zeus/tmp/unitigs/assembly_graph_final.gfa", _outputDir);
	//writeReferencePositions(unitigGraphNew, "/pasteur/appa/homes/gbenoit/zeus/tmp/unitigs/assembly_graph_final.gfa.referencePos.csv");

	for(UnitigGraph2::UnitigNode* node : unitigGraphNew->_unitigs){

		if(node == nullptr) continue;

		
		if(node->_unitigName == 326){
			cout << "utg" << node->_unitigName << endl;
			for(size_t i=0; i<node->_abundances.size(); i++){
				cout << "\t" << node->_abundances[i] << endl;
			}
		}
		
		if(_unitigName_to_minimizers[node->_unitigName].size() == _kminmerSizePrev){
			cout << "wtf: utg" << node->_unitigName << endl;
			for(UnitigType unitigIndex : node->_unitigMerge){
				cout << "\tutg" << UnitigGraph2::unitigIndexToString(unitigIndex) << endl;
			}

			for(auto m : _unitigName_to_minimizers[node->_unitigName]){
				cout << m << endl;
			}
			//getchar();
		}
	}
	*/

	//if(_kminmerSize == 6){
		//cout << "lul" << endl;
		//	getchar();
	//	exit(1);
	//}
	//cout << "Skip graph " << _kminmerSize << endl;
	//getchar();

	//exit(0);
}

/*
void CreateMdbg::extendNode(UnitigGraph2* unitigGraph, const UnitigType unitigIndex, bool checkOnly){

	if(unitigGraph->isCompactable(unitigIndex)) return;

	bool isReversed;
	const UnitigType& unitigName = UnitigGraph2::unitigIndex_to_unitigName(unitigIndex, isReversed);

	vector<MinimizerType>& minimizers = _unitigName_to_minimizers[unitigName];

	if(unitigName == 160){
		cout << endl;
		for(MinimizerType m : minimizers){
			cout << m << endl;
		}
	}



	vector<MinimizerType> endNode = UnitigGraph2::getEndNode(unitigIndex, _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(unitigIndex)], _kminmerSizePrev);

	vector<UnitigType> successors;
	unitigGraph->getSuccessors(unitigIndex, successors);

	if(!checkOnly && successors.size() != 1) return;


	for(size_t i=0; i<successors.size(); i++){

		vector<MinimizerType> startNode = UnitigGraph2::getStartNode(successors[i], _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(successors[i])], _kminmerSizePrev);
		
		if(endNode == startNode){
			cout << "\tutg" << UnitigGraph2::unitigIndexToString(unitigIndex) << " -> " << UnitigGraph2::unitigIndexToString(successors[i]) << " is ok" << endl;
		}
		else{
			cout << "\tutg" << UnitigGraph2::unitigIndexToString(unitigIndex) << " -> " << UnitigGraph2::unitigIndexToString(successors[i]) << " need to be extended" << endl;
			
			if(!checkOnly){
				if(isReversed){
					cout << "lul 1 " << endl;
					minimizers.insert(minimizers.begin(), startNode[startNode.size()-1]);
				}
				else{
					cout << "lul 2 " << endl;
					minimizers.push_back(startNode[startNode.size()-1]);
				}
			}

		}

	}
	
	if(unitigName == 160){
		//cout << endl;
		//for(MinimizerType m : minimizers){
		//	cout << m << endl;
		//}
		//cout << endl;
		//for(MinimizerType m : endNode){
		//	cout << m << endl;
		//}
		cout << endl;


		for(MinimizerType m : UnitigGraph2::getMinimizers(unitigIndex, _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(unitigIndex)])){
			cout << "\t\t\t" << m << endl;
		}
		cout << endl;

		for(size_t i=0; i<successors.size(); i++){

			
			cout << endl;
			for(MinimizerType m : UnitigGraph2::getMinimizers(successors[i], _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(successors[i])])){
				cout << "\t\t\t\t" << m << endl;
			}

			//vector<MinimizerType> startNode = UnitigGraph2::getStartNode(successors[i], _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(successors[i])], _kminmerSizePrev);
			
			//cout << endl;
			//for(auto m : startNode){
			//	cout << "\t" << m << endl;
			//}

		}
	}
}


void CreateMdbg::extendNodeBranching(UnitigGraph2* unitigGraph, const UnitigType unitigIndex, bool checkOnly){
	
	if(unitigGraph->isCompactable(unitigIndex)) return;

	bool isReversed;
	const UnitigType& unitigName = UnitigGraph2::unitigIndex_to_unitigName(unitigIndex, isReversed);

	//if(unitigName != 175) return;

	vector<MinimizerType>& minimizers = _unitigName_to_minimizers[unitigName];

	vector<MinimizerType> endNode = UnitigGraph2::getEndNode(unitigIndex, _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(unitigIndex)], _kminmerSizePrev);

	vector<UnitigType> successors;
	unitigGraph->getSuccessors(unitigIndex, successors);

	if(!checkOnly && successors.size() <= 1) return;


	for(size_t i=0; i<successors.size(); i++){

		vector<MinimizerType> startNode = UnitigGraph2::getStartNode(successors[i], _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(successors[i])], _kminmerSizePrev);
		
		if(endNode == startNode){
			cout << "\tutg" << UnitigGraph2::unitigIndexToString(unitigIndex) << " -> " << UnitigGraph2::unitigIndexToString(successors[i]) << " is ok" << endl;
		}
		else{
			cout << "\tutg" << UnitigGraph2::unitigIndexToString(unitigIndex) << " -> " << UnitigGraph2::unitigIndexToString(successors[i]) << " need to be extended" << endl;
			
			if(!checkOnly){
				
				//cout << "Removing: " << UnitigGraph2::unitigIndexToString(unitigIndex) << " -> " << UnitigGraph2::unitigIndexToString(successors[i]) << endl;
				unitigGraph->removeSuccessor(unitigIndex, successors[i]);
				unitigGraph->removeSuccessor(UnitigGraph2::unitigIndex_to_reverseDirection(successors[i]), UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndex));

				vector<MinimizerType> minimizers = endNode;
				minimizers.push_back(startNode[startNode.size()-1]);
				//if(!connectTwoSmallNodes){
				//	minimizers.insert(minimizers.begin(), endNode[0]);
				//}

				UnitigGraph2::UnitigNode* edgeNode = unitigGraph->addNode(unitigGraph->_unitigs.size(), minimizers);
				//edgeNodes[unitigPair] = edgeNode;
				edgeNode->_abundance = 40;
				_unitigName_to_minimizers[edgeNode->_unitigName] = minimizers;

				const UnitigType edgeNodeUnitigIndex = UnitigGraph2::unitigName_to_unitigIndex(edgeNode->_unitigName, false);
				unitigGraph->addSuccessor(unitigIndex, edgeNodeUnitigIndex);
				unitigGraph->addSuccessor(UnitigGraph2::unitigIndex_to_reverseDirection(edgeNodeUnitigIndex), UnitigGraph2::unitigIndex_to_reverseDirection(unitigIndex));

				unitigGraph->addSuccessor(edgeNodeUnitigIndex, successors[i]);
				unitigGraph->addSuccessor(UnitigGraph2::unitigIndex_to_reverseDirection(successors[i]), UnitigGraph2::unitigIndex_to_reverseDirection(edgeNodeUnitigIndex));

			}

		}

	}

	if(unitigName == 175){
		//cout << endl;
		//for(MinimizerType m : minimizers){
		//	cout << m << endl;
		//}
		//cout << endl;
		//for(MinimizerType m : endNode){
		//	cout << m << endl;
		//}
		cout << endl;


		for(MinimizerType m : UnitigGraph2::getMinimizers(unitigIndex, _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(unitigIndex)])){
			cout << "\t\t\t" << m << endl;
		}
		cout << endl;

		for(size_t i=0; i<successors.size(); i++){

			
			cout << endl;
			for(MinimizerType m : UnitigGraph2::getMinimizers(successors[i], _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(successors[i])])){
				cout << "\t\t\t\t" << m << endl;
			}

			//vector<MinimizerType> startNode = UnitigGraph2::getStartNode(successors[i], _unitigName_to_minimizers[UnitigGraph2::unitigIndex_to_unitigName(successors[i])], _kminmerSizePrev);
			
			//cout << endl;
			//for(auto m : startNode){
			//	cout << "\t" << m << endl;
			//}

		}
	}
}
*/
/*
void CreateMdbg::computeUnitigSequences(UnitigGraph2* unitigGraph, string outputFilename){

	bool writeUnitigs = outputFilename != "";

	if(!writeUnitigs){
		outputFilename = _outputDir + "/unitigGraph.nodes.bin.dummy";
	}

	UnitigType unitigIndexNew = 0;
	vector<UnitigType> unitigNameMap;

	for(size_t i=0; i<unitigGraph->_unitigs.size(); i++){

		UnitigGraph2::UnitigNode* node = unitigGraph->_unitigs[i];

		if(node == nullptr) continue;

		unitigNameMap.push_back(node->_unitigName);

	}
	





	int nbUnitigsWritten = 0;

	ofstream outputContigFile(_outputDir + "/contigs.nodepath");

	for(size_t i=0; i<unitigGraph->_unitigs.size(); i++){

		UnitigGraph2::UnitigNode* node = unitigGraph->_unitigs[i];

		if(node == nullptr) continue;

		vector<UnitigType> unitigs;
		if(node->_unitigMerge.size() == 0){
			unitigs.push_back(UnitigGraph2::unitigName_to_unitigIndex(node->_unitigName, false));    
		}
		else{
			unitigs = node->_unitigMerge;
		}

		if(node->_isReversed){
            vector<UnitigType> unitigsRC;
            UnitigGraph2::reverseComplementUnitigs(unitigs, unitigsRC);

			unitigs = unitigsRC;
		}

		bool isCircular = false;
		u_int32_t size = unitigs.size();

		outputContigFile.write((const char*)&size, sizeof(size));
		outputContigFile.write((const char*)&isCircular, sizeof(isCircular));
		outputContigFile.write((const char*)&unitigs[0], size * sizeof(UnitigType));

		nbUnitigsWritten += 1;

	}

	cout << "Nb written unitigs: " << nbUnitigsWritten << endl;

	outputContigFile.close();








	ofstream unitigGraphFile_nodes_tmp(_outputDir + "/unitigGraph.nodes.bin.tmp");

	for(const auto& it : _unitigName_to_minimizers){

		const UnitigType unitigName = it.first * 2;

		//if(unitigGraph->_unitigs[unitigName] == nullptr) continue;

		const vector<MinimizerType>& minimizers = it.second;
		u_int32_t size = minimizers.size();

		//if(minimizers.size() == 0){
		//	cout << "derp 1: " << unitigName << " " << minimizers.size() << endl;
		//}
		unitigGraphFile_nodes_tmp.write((const char*)&size, sizeof(size));
		unitigGraphFile_nodes_tmp.write((const char*)&minimizers[0], size * sizeof(MinimizerType));
		unitigGraphFile_nodes_tmp.write((const char*)&unitigName, sizeof(unitigName));

	}

	unitigGraphFile_nodes_tmp.close();





	_unitigName_to_minimizers.clear();

	string command = "/pasteur/appa/homes/gbenoit/zeus/tools/merging/metaMDBG_1.4/metaMDBG/build/bin/metaMDBG toMinspace " + _outputDir + " " + _outputDir + "/contigs.nodepath " + _outputDir + "/unitigs.bin " + _outputDir + "/unitigGraph.nodes.bin.tmp " + " --threads " + to_string(_nbCores);
	Utils::executeCommand(command, _outputDir);


	
	ofstream unitigGraphFile_nodes_dummy(outputFilename);

	KminmerParserParallel parser(_outputDir + "/unitigs.bin", 0, 0, false, false, 1);
	parser.parseSequences(RemapUnitigFunctor(*this, unitigGraphFile_nodes_dummy, unitigNameMap, writeUnitigs));

	unitigGraphFile_nodes_dummy.close();


	
}
*/
