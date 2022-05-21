

#include <Bloocoo.hpp>

Bloocoo::Bloocoo () : Tool()
{

}


void Bloocoo::parseArgs(int argc, char* argv[]){



	cxxopts::Options options("Graph", "Create MDBG");
	options.add_options()
	//("d,debug", "Enable debugging") // a bool parameter
	(ARG_INPUT_FILENAME, "", cxxopts::value<string>())
	(ARG_OUTPUT_DIR, "", cxxopts::value<string>())
	(ARG_FIRST_PASS, "", cxxopts::value<bool>()->default_value("false"))
	(ARG_INPUT_FILENAME_CONTIG, "", cxxopts::value<string>()->default_value(""))
		(ARG_NB_CORES, "", cxxopts::value<int>()->default_value(NB_CORES_DEFAULT));
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

    }
    catch (const std::exception& e){
		std::cout << options.help() << std::endl;
    	std::cerr << e.what() << std::endl;
    	std::exit(EXIT_FAILURE);
    }


	
	string filename_parameters = _outputDir + "/parameters.gz";
	gzFile file_parameters = gzopen(filename_parameters.c_str(),"rb");
	gzread(file_parameters, (char*)&_minimizerSize, sizeof(_minimizerSize));
	gzread(file_parameters, (char*)&_kminmerSize, sizeof(_kminmerSize));
	gzread(file_parameters, (char*)&_minimizerDensity, sizeof(_minimizerDensity));
	gzread(file_parameters, (char*)&_kminmerSizeFirst, sizeof(_kminmerSizeFirst));
	gzclose(file_parameters);


	cout << endl;
	cout << "Contig filename: " << _filename_inputContigs << endl;
	cout << "Output dir: " << _outputDir << endl;
	cout << "Minimizer length: " << _minimizerSize << endl;
	cout << "Kminmer length: " << _kminmerSize << endl;
	cout << "Density: " << _minimizerDensity << endl;
	cout << endl;

	//_inputDir = getInput()->get(STR_INPUT_DIR) ? getInput()->getStr(STR_INPUT_DIR) : "";
	//_input_extractKminmers= getInput()->get(STR_INPUT_EXTRACT_KMINMERS) ? getInput()->getStr(STR_INPUT_EXTRACT_KMINMERS) : "";

	_filename_readMinimizers = _outputDir + "/read_data.gz";
	_filename_contigMinimizers = _outputDir + "/contig_data.gz";
	//_filename_outputContigs = _inputDir + "/contigs.min.gz";
	//_filename_readCompositions = _outputDir + "/read_compositions.gz";
	//_filename_filteredMinimizers = _outputDir + "/filteredMinimizers.gz";
	//_filename_hifiasmGroundtruth = _outputDir + "/hifiasmGroundtruth.gz";


	_minimizerSpacingMean = 1 / _minimizerDensity;
	_kminmerLengthMean = _minimizerSpacingMean * (_kminmerSize-1);
	_kminmerOverlapMean = _kminmerLengthMean - _minimizerSpacingMean;

	_filename_noKminmerReads = _outputDir + "/" + FILENAME_NO_KMINMER_READS;
	//cout << _minimizerSpacingMean << endl;
	//cout << _kminmerLengthMean << endl;
	//cout << _kminmerOverlapMean << endl;
	//getchar();
}

void Bloocoo::execute (){

	_file_noKminmerReads = ofstream(_filename_noKminmerReads, std::ofstream::app);
	_node_id = 0;
	//parseArgs();
	createMDBG();
	createGfa();

	if(_isFirstPass){
		ofstream file_data(_outputDir + "/data.txt");

		file_data.write((const char*)&_nbReads, sizeof(_nbReads));

		file_data.close();
	}

	_file_noKminmerReads.close();


}


void Bloocoo::createMDBG (){

	
	//_bloomFilter = new BloomCacheCoherent<u_int64_t>(16000000000ull);
	//_bloomFilter = new BloomCacheCoherent<u_int64_t>(16000000000ull);

	_fileSmallContigs = ofstream(_outputDir + "/small_contigs.bin", std::ios_base::app);
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

	string inputFilename = _inputFilename;

	if(!_isFirstPass){
		
		_mdbgInit = new MDBG(_kminmerSize-1);
		_mdbgInit->load(_outputDir + "/mdbg_nodes.gz");

		_readFile = ofstream(_outputDir + "/read_data.txt");


	}

	string inputFilename_min;
	bool usePos = false;
	
	_kminmerFile = ofstream(_outputDir + "/kminmerData_min.txt");

	cout << "Extracting kminmers (read)" << endl;

	_parsingContigs = false;
	
	if(_isFirstPass){
		
		usePos = false;
		inputFilename_min = _outputDir + "/read_data_init.txt";

		cout << "Filling bloom filter" << endl;
		_nbSolidKminmers = 0;
		//KminmerParserParallel parser(inputFilename_min, _minimizerSize, _kminmerSize, usePos, _nbCores);
		//parser.parse(FillBloomFilter(*this));
		cout << "Solid kminmers: " << _nbSolidKminmers << endl;
		//getchar();
		
		cout << "Building mdbg" << endl;
		KminmerParserParallel parser2(inputFilename_min, _minimizerSize, _kminmerSize, usePos, _nbCores);
		parser2.parse(IndexKminmerFunctor(*this, false));
		//auto fp = std::bind(&Bloocoo::createMDBG_collectKminmers_minspace_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
		//parser.parseMinspace(fp);
	}

	if(!_isFirstPass){
		const string& filename_uncorrectedReads = _outputDir + "/read_uncorrected.txt";

		cout << "Filling bloom filter" << endl;
		//KminmerParserParallel parser(filename_uncorrectedReads, _minimizerSize, _kminmerSize, false, _nbCores);
		//parser.parse(FillBloomFilter(*this));

		cout << "Building mdbg" << endl;
		KminmerParserParallel parser2(filename_uncorrectedReads, _minimizerSize, _kminmerSize, false, _nbCores);
		parser2.parse(IndexKminmerFunctor(*this, false));

		_parsingContigs = true;
		const string& filename_contigs = _outputDir + "/unitig_data.txt";

		KminmerParserParallel parser3(filename_contigs, _minimizerSize, _kminmerSize, false, _nbCores);
		parser3.parse(IndexKminmerFunctor(*this, true));
		
	}

	//cout << _mdbg->_dbg_nodes.size() << endl;
	
	
	if(!_isFirstPass){
		delete _mdbgInit;
		_readFile.close();
		
		const auto copyOptions = fs::copy_options::overwrite_existing;
		fs::copy(_outputDir + "/read_uncorrected.txt", _outputDir + "/read_data.txt", copyOptions);
	}


	writtenNodeNames.clear();

	if(_isFirstPass){
		//if(fs::exists(_outputDir + "/read_data_init.txt")){
		//	fs::remove(_outputDir + "/read_data_init.txt");
		//}
		const auto copyOptions = fs::copy_options::overwrite_existing;
		fs::copy(_outputDir + "/read_data_init.txt", _outputDir + "/read_data.txt", copyOptions);
		fs::copy(_outputDir + "/read_data_init.txt", _outputDir + "/read_uncorrected.txt", copyOptions); //disable if correction is enabled
	}


	cout << "Nb kminmers (reads): " << _kminmerExist.size() << endl;
	cout << "Nb solid kminmers: " << _mdbg->_dbg_nodes.size() << endl;
	
	_kminmerExist.clear();
	_minimizerCounts.clear();
	_kminmersData.clear();
	_contigIndex.clear();
	_contigAbundances.clear();
	_fileSmallContigs.close();

}


//void Bloocoo::createMDBG_collectKminmers_minspace_read(const vector<u_int64_t>& readMinimizers, const vector<ReadKminmerComplete>& kminmersInfos, u_int64_t readIndex){



//}


struct KmerVecSorterData{
	u_int32_t _nodeName;
	KmerVec _kmerVec;
};

static bool KmerVecComparator(const KmerVecSorterData &a, const KmerVecSorterData &b){

	for(size_t i =0; i<a._kmerVec._kmers.size(); i++){
		if(a._kmerVec._kmers[i] == b._kmerVec._kmers[i]){
			continue;
		}
		else{
			return a._kmerVec._kmers[i] > b._kmerVec._kmers[i];
		}
	}
	
	return a._kmerVec._kmers[0] > b._kmerVec._kmers[0];
	//if(a._startNodeIndex == b._startNodeIndex){
	//    return a._length < b._length;
	//}
	//return a._startNodeIndex < b._startNodeIndex;
	//if(a._length == b._length){
	//    return a._startNodeIndex < b._startNodeIndex;
	//}
	//return a._length < b._length;
}


void Bloocoo::createGfa(){

	cout << "Writing gfa..." << endl;
	//cout << "Cleaning repeats..." << endl;

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
		KmerVec vec = it.first;

		it.second._index = _nodeName_to_deterministicNodeName[it.second._index];
		
		//KmerVecSorterData& d = kmerVecs[i];
		u_int32_t nodeName = it.second._index;


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


	}
	
	


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

	_kminmerFile.close();

	//delete mdbg_repeatFree;
	//cout << _dbg_edges.size() << endl;

	u_int64_t nbEdges = 0;

	string gfa_filename = _outputDir + "/minimizer_graph.gfa";
	ofstream output_file_gfa(gfa_filename);

	for(const auto& vec_id : _mdbg->_dbg_nodes){

		KmerVec vec = vec_id.first;
		KmerVec vec_rev = vec_id.first.reverse();
		u_int32_t id = vec_id.second._index;
		//u_int32_t id = vec_id.second._index;

		output_file_gfa << "S" << "\t" << id << "\t" << "*" << "\t" << "LN:i:" << _kminmerLengthMean << "\t" << "dp:i:" << vec_id.second._abundance << endl;

		//cout << mdbg->_dbg_edges[vec.prefix().normalize()].size() << endl;
		for(KmerVec& v : _mdbg->_dbg_edges[vec.prefix().normalize()]){
			if(v==vec) continue;
			KmerVec v_rev = v.reverse();
			u_int32_t id2 = _mdbg->_dbg_nodes[v]._index;

			if (vec.suffix() == v.prefix()) {
				nbEdges += 1;
				//u_int16_t overlapLength =  _mdbg->_dbg_nodes[v]._length - _mdbg->_dbg_nodes[v]._overlapLength_end; //   min(_mdbg->_dbg_nodes[v]._overlapLength_end, _mdbg->_dbg_nodes[vec]._overlapLength_end);
				output_file_gfa << "L" << "\t" << id << "\t" << "+" << "\t" << id2 << "\t" << "+" << "\t" << _kminmerOverlapMean << "M" << endl;
				//vec_add_edge("+", "+");
			}
			if (vec.suffix() == v_rev.prefix()) {
				nbEdges += 1;
				//u_int16_t overlapLength = _mdbg->_dbg_nodes[v]._length - _mdbg->_dbg_nodes[v]._overlapLength_start; //min(_mdbg->_dbg_nodes[v]._overlapLength_end, _mdbg->_dbg_nodes[vec]._overlapLength_end);
				output_file_gfa << "L" << "\t" << id << "\t" << "+" << "\t" << id2 << "\t" << "-" << "\t" << _kminmerOverlapMean << "M" << endl;
				//vec_add_edge("+", "-");
			}
			if (vec_rev.suffix() == v.prefix()) {
				nbEdges += 1;
				//u_int16_t overlapLength = _mdbg->_dbg_nodes[v]._length - _mdbg->_dbg_nodes[v]._overlapLength_end; //min(_mdbg->_dbg_nodes[v]._overlapLength_start, _mdbg->_dbg_nodes[vec]._overlapLength_start);
				output_file_gfa << "L" << "\t" << id << "\t" << "-" << "\t" << id2 << "\t" << "+" << "\t" << _kminmerOverlapMean << "M" << endl;
				//vec_add_edge("-", "+");
			}
			if (vec_rev.suffix() == v_rev.prefix()) {
				nbEdges += 1;
				//u_int16_t overlapLength = _mdbg->_dbg_nodes[v]._length - _mdbg->_dbg_nodes[v]._overlapLength_start;; //min(_mdbg->_dbg_nodes[v]._overlapLength_start, _mdbg->_dbg_nodes[vec]._overlapLength_start);
				output_file_gfa << "L" << "\t" << id << "\t" << "-" << "\t" << id2 << "\t" << "-" << "\t" << _kminmerOverlapMean << "M" << endl;
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
				output_file_gfa << "L" << "\t" << id << "\t" << "+" << "\t" << id2 << "\t" << "+" << "\t" << _kminmerOverlapMean << "M" << endl;
				//vec_add_edge("+", "+");
			}
			if (vec.suffix() == v_rev.prefix()) {
				nbEdges += 1;
				//u_int16_t overlapLength = _mdbg->_dbg_nodes[v]._length - _mdbg->_dbg_nodes[v]._overlapLength_start; //min(_mdbg->_dbg_nodes[v]._overlapLength_end, _mdbg->_dbg_nodes[vec]._overlapLength_end);
				//cout << overlapLength << " " << _mdbg->_dbg_nodes[v]._length<< endl;
				output_file_gfa << "L" << "\t" << id << "\t" << "+" << "\t" << id2 << "\t" << "-" << "\t" << _kminmerOverlapMean << "M" << endl;
				//vec_add_edge("+", "-");
			}
			if (vec_rev.suffix() == v.prefix()) {
				nbEdges += 1;
				//u_int16_t overlapLength = _mdbg->_dbg_nodes[v]._length - _mdbg->_dbg_nodes[v]._overlapLength_end; //min(_mdbg->_dbg_nodes[v]._overlapLength_start, _mdbg->_dbg_nodes[vec]._overlapLength_start);
				//cout << overlapLength << " " << _mdbg->_dbg_nodes[v]._length << endl;
				output_file_gfa << "L" << "\t" << id << "\t" << "-" << "\t" << id2 << "\t" << "+" << "\t" << _kminmerOverlapMean << "M" << endl;
				//vec_add_edge("-", "+");
			}
			if (vec_rev.suffix() == v_rev.prefix()) {
				nbEdges += 1;
				//u_int16_t overlapLength = _mdbg->_dbg_nodes[v]._length - _mdbg->_dbg_nodes[v]._overlapLength_start; //min(_mdbg->_dbg_nodes[v]._overlapLength_start, _mdbg->_dbg_nodes[vec]._overlapLength_start);
				//cout << overlapLength << " " << _mdbg->_dbg_nodes[v]._length<< endl;
				output_file_gfa << "L" << "\t" << id << "\t" << "-" << "\t" << id2 << "\t" << "-" << "\t" << _kminmerOverlapMean << "M" << endl;
				//vec_add_edge("-", "-");
			}
		}

		/*
			#pragma omp parallel num_threads(_nbCores)
	{

		#pragma omp for
		for(size_t i=0; i<keys.size(); i++){
			KmerVec vec = keys[i];
			KmerVec vec_rev = vec.reverse();
			u_int32_t id = _mdbg->_dbg_nodes[vec]._index;
			//u_int32_t id = vec_id.second._index;

			#pragma omp critical(gfa)
			{
				output_file_gfa << "S" << "\t" << id << "\t" << "*" << "\t" << "LN:i:" << _kminmerLengthMean << "\t" << "dp:i:" << _mdbg->_dbg_nodes[vec]._abundance << endl;
			}
			
			//cout << mdbg->_dbg_edges[vec.prefix().normalize()].size() << endl;
			for(KmerVec& v : _mdbg->_dbg_edges[vec.prefix().normalize()]){
				if(v==vec) continue;
				KmerVec v_rev = v.reverse();
				u_int32_t id2 = _mdbg->_dbg_nodes[v]._index;

				if (vec.suffix() == v.prefix()) {
					
					#pragma omp critical(gfa)
					{
					nbEdges += 1;
					//u_int16_t overlapLength =  _mdbg->_dbg_nodes[v]._length - _mdbg->_dbg_nodes[v]._overlapLength_end; //   min(_mdbg->_dbg_nodes[v]._overlapLength_end, _mdbg->_dbg_nodes[vec]._overlapLength_end);
					output_file_gfa << "L" << "\t" << id << "\t" << "+" << "\t" << id2 << "\t" << "+" << "\t" << _kminmerOverlapMean << "M" << endl;
					//vec_add_edge("+", "+");
					}
				}
				if (vec.suffix() == v_rev.prefix()) {
					
					#pragma omp critical(gfa)
					{
					nbEdges += 1;
					//u_int16_t overlapLength = _mdbg->_dbg_nodes[v]._length - _mdbg->_dbg_nodes[v]._overlapLength_start; //min(_mdbg->_dbg_nodes[v]._overlapLength_end, _mdbg->_dbg_nodes[vec]._overlapLength_end);
					output_file_gfa << "L" << "\t" << id << "\t" << "+" << "\t" << id2 << "\t" << "-" << "\t" << _kminmerOverlapMean << "M" << endl;
					//vec_add_edge("+", "-");
					}
				}
				if (vec_rev.suffix() == v.prefix()) {
					
					#pragma omp critical(gfa)
					{
					nbEdges += 1;
					//u_int16_t overlapLength = _mdbg->_dbg_nodes[v]._length - _mdbg->_dbg_nodes[v]._overlapLength_end; //min(_mdbg->_dbg_nodes[v]._overlapLength_start, _mdbg->_dbg_nodes[vec]._overlapLength_start);
					output_file_gfa << "L" << "\t" << id << "\t" << "-" << "\t" << id2 << "\t" << "+" << "\t" << _kminmerOverlapMean << "M" << endl;
					//vec_add_edge("-", "+");
					}
				}
				if (vec_rev.suffix() == v_rev.prefix()) {
					
					#pragma omp critical(gfa)
					{
					nbEdges += 1;
					//u_int16_t overlapLength = _mdbg->_dbg_nodes[v]._length - _mdbg->_dbg_nodes[v]._overlapLength_start;; //min(_mdbg->_dbg_nodes[v]._overlapLength_start, _mdbg->_dbg_nodes[vec]._overlapLength_start);
					output_file_gfa << "L" << "\t" << id << "\t" << "-" << "\t" << id2 << "\t" << "-" << "\t" << _kminmerOverlapMean << "M" << endl;
					//vec_add_edge("-", "-");
					}
				}

			}
			for(KmerVec& v : _mdbg->_dbg_edges[vec.suffix().normalize()]){
				if(v==vec) continue;
				KmerVec v_rev = v.reverse();
				u_int32_t id2 = _mdbg->_dbg_nodes[v]._index;

				if (vec.suffix() == v.prefix()) {
					
					#pragma omp critical(gfa)
					{
					nbEdges += 1;
					//u_int16_t overlapLength = _mdbg->_dbg_nodes[v]._length - _mdbg->_dbg_nodes[v]._overlapLength_end; //min(_mdbg->_dbg_nodes[v]._overlapLength_end, _mdbg->_dbg_nodes[vec]._overlapLength_end);
					//cout << overlapLength << " " << _mdbg->_dbg_nodes[v]._length<< endl;
					output_file_gfa << "L" << "\t" << id << "\t" << "+" << "\t" << id2 << "\t" << "+" << "\t" << _kminmerOverlapMean << "M" << endl;
					//vec_add_edge("+", "+");
					}
				}
				if (vec.suffix() == v_rev.prefix()) {
					
					#pragma omp critical(gfa)
					{
					nbEdges += 1;
					//u_int16_t overlapLength = _mdbg->_dbg_nodes[v]._length - _mdbg->_dbg_nodes[v]._overlapLength_start; //min(_mdbg->_dbg_nodes[v]._overlapLength_end, _mdbg->_dbg_nodes[vec]._overlapLength_end);
					//cout << overlapLength << " " << _mdbg->_dbg_nodes[v]._length<< endl;
					output_file_gfa << "L" << "\t" << id << "\t" << "+" << "\t" << id2 << "\t" << "-" << "\t" << _kminmerOverlapMean << "M" << endl;
					//vec_add_edge("+", "-");
					}
				}
				if (vec_rev.suffix() == v.prefix()) {
					
					#pragma omp critical(gfa)
					{
					nbEdges += 1;
					//u_int16_t overlapLength = _mdbg->_dbg_nodes[v]._length - _mdbg->_dbg_nodes[v]._overlapLength_end; //min(_mdbg->_dbg_nodes[v]._overlapLength_start, _mdbg->_dbg_nodes[vec]._overlapLength_start);
					//cout << overlapLength << " " << _mdbg->_dbg_nodes[v]._length << endl;
					output_file_gfa << "L" << "\t" << id << "\t" << "-" << "\t" << id2 << "\t" << "+" << "\t" << _kminmerOverlapMean << "M" << endl;
					//vec_add_edge("-", "+");
					}
				}
				if (vec_rev.suffix() == v_rev.prefix()) {
					
					#pragma omp critical(gfa)
					{
					nbEdges += 1;
					//u_int16_t overlapLength = _mdbg->_dbg_nodes[v]._length - _mdbg->_dbg_nodes[v]._overlapLength_start; //min(_mdbg->_dbg_nodes[v]._overlapLength_start, _mdbg->_dbg_nodes[vec]._overlapLength_start);
					//cout << overlapLength << " " << _mdbg->_dbg_nodes[v]._length<< endl;
					output_file_gfa << "L" << "\t" << id << "\t" << "-" << "\t" << id2 << "\t" << "-" << "\t" << _kminmerOverlapMean << "M" << endl;
					//vec_add_edge("-", "-");
					}
				}
			}
			//foo(element->first, element->second);
		}
		*/
		//for(KmerVec& v : _dbg_edges_prefix[vec.suffix().normalize()]){
		//	output_file_gfa << "L" << "\t" << vec_id.second << "\t" << "+" << "\t" << _dbg_nodes[v] << "\t" << "+" << "\t" << "0M" << endl;
		//}
	}

	output_file_gfa.close();
	
	cout << "Nb nodes: " << _mdbg->_dbg_nodes.size() << endl;
	cout << "Nb edges: " << nbEdges << endl;


	if(_isFirstPass){
		_mdbg->dump(_outputDir + "/mdbg_nodes_init.gz");
	}
	
	_mdbgNoFilter->dump(_outputDir + "/mdbg_nodes_noFilter.gz");
	_mdbg->dump(_outputDir + "/mdbg_nodes.gz");
	//_mdbg->_dbg_nodes.clear();
	//_mdbg->_dbg_edges.clear();



}


