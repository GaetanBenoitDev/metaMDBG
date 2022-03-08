

#include <Bloocoo.hpp>

Bloocoo::Bloocoo () : Tool()
{


	/*
	getParser()->push_back (new OptionOneParam (STR_INPUT, "input file", true));
	getParser()->push_back (new OptionOneParam (STR_OUTPUT, "output dir", true));
	getParser()->push_back (new OptionOneParam (STR_MINIM_SIZE, "minimizer length", false, "21"));
	getParser()->push_back (new OptionOneParam (STR_KMINMER_SIZE, "k-min-mer length", false, "3"));
	getParser()->push_back (new OptionOneParam (STR_DENSITY, "density of minimizers", false, "0.005"));*/
	//getParser()->push_back (new OptionOneParam (STR_INPUT_DIR, "input dir", false, ""));
	//getParser()->push_back (new OptionOneParam (STR_INPUT_EXTRACT_KMINMERS, "", false, ""));

}





void Bloocoo::execute (){
	//parseArgs();
	createMDBG();
	createGfa();
}

void Bloocoo::parseArgs(int argc, char* argv[]){



	cxxopts::Options options("Graph", "Create MDBG");
	options.add_options()
	//("d,debug", "Enable debugging") // a bool parameter
	(ARG_INPUT_FILENAME, "", cxxopts::value<string>())
	(ARG_OUTPUT_DIR, "", cxxopts::value<string>())
	(ARG_FIRST_PASS, "", cxxopts::value<bool>()->default_value("false"))
	(ARG_INPUT_FILENAME_CONTIG, "", cxxopts::value<string>()->default_value(""));
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

    }
    catch (const std::exception& e){
		std::cout << options.help() << std::endl;
    	std::cerr << e.what() << std::endl;
    	std::exit(EXIT_FAILURE);
    }




    //T exception{text};




	//_inputFilename = result[STR_INPUT].as<string>();

	//_compositionManager = new CompositionManager(4);

	//int w = 80;
	


	//if(!System::file().doesExist (_outputDir)) System::file().mkdir(_outputDir, -1);
	
	string filename_parameters = _outputDir + "/parameters.gz";
	gzFile file_parameters = gzopen(filename_parameters.c_str(),"rb");
	gzread(file_parameters, (char*)&_minimizerSize, sizeof(_minimizerSize));
	gzread(file_parameters, (char*)&_kminmerSize, sizeof(_kminmerSize));
	gzread(file_parameters, (char*)&_minimizerDensity, sizeof(_minimizerDensity));
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




}

/*
void Bloocoo::createMDBG_read(kseq_t* read, u_int64_t readIndex){
	
	if(readIndex % 100000 == 0) cout << readIndex << " " << _debug_nbMinimizers << endl;

	string rleSequence;
	vector<u_int64_t> rlePositions;
	Encoder::encode_rle(read->seq.s, strlen(read->seq.s), rleSequence, rlePositions);

	vector<u_int64_t> minimizers;
	vector<u_int64_t> minimizers_pos;
	_minimizerParser->parse(rleSequence, minimizers, minimizers_pos);
	_debug_nbMinimizers += minimizers.size();

	//for(u_int64_t minimizer : minimizers){
	//	_minimizerCounts[minimizer] += 1;
	//}
	for(size_t i=0; i<rlePositions.size(); i++){
		rlePositions[i] = i;
	}

	vector<KmerVec> kminmers; 
	vector<ReadKminmer> kminmersInfo;
	MDBG::getKminmers(_minimizerSize, _kminmerSize, minimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, readIndex);

	for(size_t i=0; i<kminmers.size(); i++){

		if(_kminmersData.find(kminmers[i]) == _kminmersData.end()){
			_kminmersData[kminmers[i]] = {0, kminmersInfo[i]._length - _minimizerSize, kminmersInfo[i]._seq_length_start, kminmersInfo[i]._seq_length_end, kminmersInfo[i]._isReversed};
		}

		_kminmersData[kminmers[i]]._count += 1;

	}

	u_int16_t size = minimizers.size();
	gzwrite(_file_readData, (const char*)&size, sizeof(size));
	gzwrite(_file_readData, (const char*)&minimizers[0], size * sizeof(u_int64_t));

}
*/





/*
void Bloocoo::createMDBG_collectKminmers_contig(const vector<u_int64_t>& minimizers, const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex){
	


	for(size_t i=0; i<kminmers.size(); i++){

		//if(kminmers[i].isPalindrome()){
		//	cout << kminmers[i]._kmers[0] << " " << kminmers[i]._kmers[1] << " " << kminmers[i]._kmers[2] << " " << kminmers[i]._kmers[3] << endl;
 			//getchar();
		//}

		if(_kminmersData.find(kminmers[i]) == _kminmersData.end()){
			_kminmersData[kminmers[i]] = {0, kminmersInfos[i]._length - _minimizerSize, kminmersInfos[i]._seq_length_start, kminmersInfos[i]._seq_length_end, kminmersInfos[i]._isReversed};
			_kminmersData[kminmers[i]]._count = 2;
		}




	}
}
*/

/*
void Bloocoo::createMDBG_index(const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex){

	
	if(readIndex == 0){
		cout << _lala._kmers.size() << endl;
		if(_lala._kmers.size() != 0){
			for(u_int64_t minimizer : _lala._kmers) cout << minimizer << " ";
			cout << endl;
		}
		//cout << _lala._kmers[0] << " " << _lala._kmers[1] << " " << _lala._kmers[2] << endl;
		getchar(); 
	}
	vector<u_int32_t> abundances;

	for(const KmerVec& vec : kminmers){
		abundances.push_back(_kminmersData[vec]._count);
	}

	float minAbundance_cutoff = compute_median(abundances) / 2;
	u_int32_t nodeNameLala = -1;

	for(size_t i=0; i<kminmers.size(); i++){

		const KmerVec& vec = kminmers[i];
		if (vec.h() == _lala.h()){
			nodeNameLala = _mdbg->_dbg_nodes[vec]._index;
			cout << _mdbg->_dbg_nodes[vec]._index << endl;
			getchar();
		}

		if(_kminmersData[vec]._count == 1) continue; //ATTENTION, ON ENELEVE TOUS LES KMINMER VU UNE FOIS
		//if(kminmerCounts[vec] > 1000) cout << kminmerCounts[vec] << endl;
		if(_kminmersData[vec]._count <= 2){
			//if(_kminmersData[vec]._count <= minAbundance_cutoff) continue;
		}


		//cout << kminmersData[vec] << endl;

		KminmerData& kminmerData = _kminmersData[vec];
		//cout << kminmerData._length << endl;
		_mdbg->addNode(vec, kminmerData._length, kminmerData._overlapLength_start, kminmerData._overlapLength_end, kminmerData._isReversed);
	
		if (vec.h() == _lala.h()){
			nodeNameLala = _mdbg->_dbg_nodes[vec]._index;
			cout << _mdbg->_dbg_nodes[vec]._index << endl;
			getchar();
		}
	}

	if(nodeNameLala != -1){
		cout << nodeNameLala << endl;
		getchar();
		nodeNameLala = -1;
	}

}
*/

/*
void Bloocoo::createMDBG_collectKminmers_contigs(const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex){

	u_int32_t minContigAbundance = -1;

	vector<u_int32_t> abundances;


	bool requireIndexing = kminmers.size() > (10000*_minimizerDensity*2);
	u_int64_t nbExistingKminmers = 0;

	for(size_t i=0; i<kminmers.size(); i++){

		const KmerVec& vec = kminmers[i];

		if(requireIndexing){
			_contigIndex[vec].push_back(readIndex);
		}

		if(_kminmersData.find(vec) == _kminmersData.end()){
			_kminmersData[vec] = {0, kminmersInfos[i]._length - _minimizerSize, kminmersInfos[i]._seq_length_start, kminmersInfos[i]._seq_length_end, kminmersInfos[i]._isReversed};
			_kminmersData[vec]._count += 1; //bypass abundance filter
			
		}
		else{
			nbExistingKminmers += 1;
		}
		

		_kminmersData[vec]._count += 1;
		
		if(_kminmersData[vec]._count < minContigAbundance){
			minContigAbundance = _kminmersData[vec]._count;
		}
		abundances.push_back(_kminmersData[vec]._count);
	}
		
	if(requireIndexing){
		u_int32_t abundance = compute_median(abundances);
		cout << "Contig abundance: " << minContigAbundance << " (nb kminmers: " << kminmers.size() << ")" << endl;
		_contigAbundances[readIndex] = minContigAbundance;
	}
	
}

void Bloocoo::removeErroneousKminmers(const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex){
	
	vector<u_int32_t> contigAbundances;

	//cout << "----" << endl;
	vector<u_int32_t> mappingContigs;
	u_int32_t minContigAbundance = -1;

	for(size_t i=0; i<kminmers.size(); i++){
		const KmerVec& vec = kminmers[i];

		if(_contigIndex.find(vec) != _contigIndex.end()){
			vector<u_int32_t> contigIndexes = _contigIndex[vec];

			for(u_int32_t contigIndex : contigIndexes){

				u_int32_t abundance = _contigAbundances[contigIndex];
				if(abundance < minContigAbundance){
					minContigAbundance = abundance;
				}
				//if(std::find(mappingContigs.begin(), mappingContigs.end(), contigIndex) == mappingContigs.end()){
				//	mappingContigs.push_back(contigIndex);
				//	contigAbundances.push_back(_contigAbundances[contigIndex]);
				//}
			}
			//cout << "Kminmers filter: " << kminmerData._count << " (Nb contigs:  " << contigIndexes.size() << ")    ";
			//if(contigIndexes.size() == 1){
				//mappingContigs.push_back(contigIndexes[0] - 2000000000 + 1);
				//cout << _contigAbundances[contigIndexes[0]];
			//}
			mappingContigs.push_back(contigIndexes[0] - 2000000000 + 1);

			//cout << endl;

			//float cutoff = ((float) _contigAbundances[contigIndexes[0]]) / 2.0;
			//if(kminmerData._count < cutoff) continue;
		}
		else{
			mappingContigs.push_back(0);
		}
	}


	//cout << "Contig abundance: " << minContigAbundance << endl;


	float cutoff = ((float) minContigAbundance) /2.0;

	for(size_t i=0; i<kminmers.size(); i++){
		const KmerVec& vec = kminmers[i];
		
		if(_kminmersData.find(vec) == _kminmersData.end()) continue; //Already removed


		const KminmerData& kminmerData = _kminmersData[vec];

		//cout << kminmerData._count << endl;
		
		if(kminmerData._count == 1){//ATTENTION, ON ENELEVE TOUS LES KMINMER VU UNE FOIS
			_kminmersData.erase(vec);
			//cout << "Erased" << endl;
			continue; 
		}


	}
}
*/

/*
void Bloocoo::createMDBG_collectKminmers(const vector<u_int64_t>& minimizers, const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex){
	


	for(size_t i=0; i<kminmers.size(); i++){

		const KmerVec& vec = kminmers[i];

		if(_kminmerExist.find(vec)){
			_mdbg->addNode(vec, kminmersInfos[i]._length - _minimizerSize, kminmersInfos[i]._seq_length_start, kminmersInfos[i]._seq_length_end, kminmersInfos[i]._isReversed);
			//_mdbg->_dbg_nodes[vec]._abundance += 1;


			
			DnaBitset* contigSequenceBitset = new DnaBitset(contigSequence);
			u_int32_t sizeData = contigSequenceBitset->_bitsetSize;
			u_int32_t sizeSeq = contigSequenceBitset->m_len;
			uint8_t* m_data = contigSequenceBitset->m_data;
			contigFile_bitset.write((const char*)&sizeData, sizeof(sizeData));
			contigFile_bitset.write((const char*)&sizeSeq, sizeof(sizeSeq));
			contigFile_bitset.write((const char*)&m_data[0], sizeData*sizeof(uint8_t));
			delete contigSequenceBitset;
			

			"reprise: besoin de la sequence du kminmer ici"
			_kminmerFile.write();
		}
		else{
			_kminmerExist.insert(vec);
		}
		//if(kminmers[i].isPalindrome()){
		//	cout << kminmers[i]._kmers[0] << " " << kminmers[i]._kmers[1] << " " << kminmers[i]._kmers[2] << " " << kminmers[i]._kmers[3] << endl;
 			//getchar();
		//}

		//if(_kminmersData.find(kminmers[i]) == _kminmersData.end()){
		//	_kminmersData[kminmers[i]] = {0, kminmersInfos[i]._length - _minimizerSize, kminmersInfos[i]._seq_length_start, kminmersInfos[i]._seq_length_end, kminmersInfos[i]._isReversed};
		//}

		//_kminmersData[kminmers[i]]._count += 1;



	}
}
*/


void Bloocoo::createMDBG_collectKminmers_read(kseq_t* read, u_int64_t readIndex){
				
	if(readIndex % 100000 == 0) cout << "Processed reads: " << readIndex << endl;
	
	string kminmerSequence;
	char* sequenceOriginal = read->seq.s;

	string rleSequence;
	vector<u_int64_t> rlePositions;
	Encoder::encode_rle(read->seq.s, strlen(read->seq.s), rleSequence, rlePositions);

	vector<u_int64_t> minimizers;
	vector<u_int64_t> minimizers_pos;
	_minimizerParser->parse(rleSequence, minimizers, minimizers_pos);

	vector<KmerVec> kminmers; 
	vector<ReadKminmer> kminmersInfos;
	MDBG::getKminmers(_minimizerSize, _kminmerSize, minimizers, minimizers_pos, kminmers, kminmersInfos, rlePositions, readIndex, false);


	vector<u_int16_t> minimizerPosOffset;

	if(minimizers.size() > 0){
		u_int16_t pos = minimizers_pos[0];
		//cout << pos << endl;
		minimizerPosOffset.push_back(pos);
		
		for(size_t i=1; i<minimizers_pos.size(); i++){
			u_int16_t posOffset = minimizers_pos[i] - pos;
			minimizerPosOffset.push_back(posOffset);
			pos = minimizers_pos[i];
			//cout << pos << " " << posOffset << endl;
		}
	}

	u_int32_t size = minimizers.size();
	_readFile.write((const char*)&size, sizeof(size));
	_readFile.write((const char*)&minimizers[0], size*sizeof(u_int64_t));
	_readFile.write((const char*)&minimizerPosOffset[0], size*sizeof(u_int16_t));

	/*
	for(auto& it : _mdbg->_dbg_nodes){
		cout << it.second._index << endl;
		if(it.second._index == 21594){
			cout << "wesh ? " << endl;
			getchar();
		}
	}
	*/

	for(size_t i=0; i<kminmers.size(); i++){


		const KmerVec& vec = kminmers[i];
		const ReadKminmer& kminmerInfo = kminmersInfos[i];


		if(_kminmerExist.find(vec) != _kminmerExist.end()){
			_mdbg->addNode(vec, kminmerInfo._length - _minimizerSize, kminmerInfo._seq_length_start, kminmerInfo._seq_length_end, kminmersInfos[i]._isReversed);


			if(_mdbg->_dbg_nodes[vec]._abundance == 2){
				
				extractKminmerSequence(sequenceOriginal, kminmerInfo, kminmerSequence);
				DnaBitset* sequenceBitset = new DnaBitset(kminmerSequence);

				//DnaBitset* contigSequenceBitset = new DnaBitset(contigSequence);
				u_int32_t sizeData = sequenceBitset->_bitsetSize;
				u_int32_t sizeSeq = sequenceBitset->m_len;
				uint8_t* m_data = sequenceBitset->m_data;
				_kminmerFile.write((const char*)&sizeData, sizeof(sizeData));
				_kminmerFile.write((const char*)&sizeSeq, sizeof(sizeSeq));
				_kminmerFile.write((const char*)&m_data[0], sizeData*sizeof(uint8_t));
				delete sequenceBitset;

				u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
				u_int32_t lengthStart = kminmerInfo._seq_length_start;
				u_int32_t lengthEnd = kminmerInfo._seq_length_end;
				//bool isReversed = kminmerInfo._isReversed;

				_kminmerFile.write((const char*)&nodeName, sizeof(nodeName));
				_kminmerFile.write((const char*)&lengthStart, sizeof(lengthStart));
				_kminmerFile.write((const char*)&lengthEnd, sizeof(lengthEnd));
				//_kminmerFile.write((const char*)&isReversed, sizeof(isReversed));

				//if(_mdbg->_dbg_nodes[vec]._index <= 0){
				//	cout << nodeName << " " << lengthStart << " " << lengthEnd << " " << isReversed << " " << kminmerSequence.size() << endl;
				//	cout << kminmerSequence << endl;
				//}
			}
		}
		else{
			_kminmerExist.insert(vec);
		}

			
	}

}



void Bloocoo::extractKminmerSequence(const char* sequenceOriginal, const ReadKminmer& kminmerInfo, string& sequence){

	u_int32_t startPosition = 0;
	u_int32_t len = 0;

	startPosition = kminmerInfo._read_pos_start;
	len = kminmerInfo._read_pos_end - kminmerInfo._read_pos_start;

	char subbuff[len+1];
	memcpy( subbuff, &sequenceOriginal[startPosition], len);
	subbuff[len] = '\0';
	sequence = string(subbuff);

	if(kminmerInfo._isReversed){
		Utils::revcomp(sequence);
	}

}



void Bloocoo::createMDBG_collectKminmers_minspace_read(const vector<u_int64_t>& readMinimizers, const vector<ReadKminmerComplete>& kminmersInfos, u_int64_t readIndex){

	//cout << readMinimizers.size() << endl;
	vector<u_int16_t> minimizerPosOffset(readMinimizers.size());
	for(size_t i=0; i<minimizerPosOffset.size(); i++){
		minimizerPosOffset[i] = 200;
	}

	u_int32_t size = readMinimizers.size();
	_readFile.write((const char*)&size, sizeof(size));
	_readFile.write((const char*)&readMinimizers[0], size*sizeof(u_int64_t));
	_readFile.write((const char*)&minimizerPosOffset[0], size*sizeof(u_int16_t));


	//cout << readIndex << " " << kminmersInfos.size() << endl;
	for(size_t i=0; i<kminmersInfos.size(); i++){
		
		const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];

		KmerVec vec = kminmerInfo._vec;
		//vec = vec.normalize();
		//if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()){
		//	cout << "Not found" << endl;
		//	continue;
		//}

		//u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
		
		if(_kminmerExist.find(vec) != _kminmerExist.end()){
			_mdbg->addNode(vec, 200, 50, 50, kminmerInfo._isReversed);
		}
		else{
			_kminmerExist.insert(vec);
		}

	}

}



void Bloocoo::createMDBG (){



	//_kminmerExist.clear();
	_mdbg = new MDBG(_kminmerSize);

	string inputFilename = _inputFilename;
	_readFile = ofstream(_outputDir + "/read_data.txt", std::ios::binary);

	if(_isFirstPass){
		_minimizerParser = new MinimizerParser(_minimizerSize, _minimizerDensity);
		_kminmerFile = ofstream(_outputDir + "/kminmerData.txt", std::ios::binary);

		cout << "Extracting kminmers (read)" << endl;
		//KminmerParser parser(_filename_readMinimizers, _minimizerSize, _kminmerSize, true);
		ReadParser readParser(_inputFilename, false, !_isFirstPass);
		auto fp = std::bind(&Bloocoo::createMDBG_collectKminmers_read, this, std::placeholders::_1, std::placeholders::_2);
		readParser.parse(fp);

		_kminmerFile.close();

		delete _minimizerParser;
		//inputFilename = _outputDir + "/read_data.txt";
	}
	else{
		
		/*
		ifstream file_uncorrectedReads(_inputDir + "/read_uncorrected.txt");

		while(true){
			
			u_int32_t size;
			vector<u_int64_t> minimizers;
			
			file_uncorrectedReads.read((char*)&size, sizeof(size));

			if(file_uncorrectedReads.eof())break;

			minimizers.resize(size);
			file_uncorrectedReads.read((char*)&minimizers[0], size*sizeof(u_int64_t));

			
			outputFile.write((const char*)&size, sizeof(size));
			outputFile.write((const char*)&minimizers[0], size*sizeof(u_int64_t));
		}

		file_uncorrectedReads.close();
		*/

		//cout << _inputFilename << endl;

		cout << "Extracting kminmers (read)" << endl;
		KminmerParser parser(_inputFilename, _minimizerSize, _kminmerSize, false);
		auto fp = std::bind(&Bloocoo::createMDBG_collectKminmers_minspace_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
		parser.parseMinspace(fp);
		
		const string& filename_uncorrectedReads = _outputDir + "/read_uncorrected.txt";

		
		KminmerParser parser2(filename_uncorrectedReads, _minimizerSize, _kminmerSize, false);
		auto fp2 = std::bind(&Bloocoo::createMDBG_collectKminmers_minspace_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
		parser2.parseMinspace(fp2);
		
	}

	_readFile.close();



	if(_isFirstPass){
		//if(fs::exists(_outputDir + "/read_data_init.txt")){
		//	fs::remove(_outputDir + "/read_data_init.txt");
		//}
		//const auto copyOptions = fs::copy_options::overwrite_existing;
		//fs::copy(_outputDir + "/read_data.txt", _outputDir + "/read_data_init.txt", copyOptions);
	}


	cout << "Nb kminmers (reads): " << _kminmerExist.size() << endl;

	/*
	if(_filename_inputContigs != ""){
		cout << "Extracting kminmers (contigs)" << endl;
		//KminmerParser contigParser(_filename_inputContigs, _minimizerSize, _kminmerSize);
		//auto fp2 = std::bind(&Bloocoo::createMDBG_collectKminmers_contigs, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
		//contigParser.parse_mContigs(fp2);

		KminmerParser parser(_filename_contigMinimizers, _minimizerSize, _kminmerSize, true);
		auto fp = std::bind(&Bloocoo::createMDBG_collectKminmers_contig, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
		parser.parse(fp);

		cout << "Nb kminmers (Contigs): " << _kminmersData.size() << endl;

		//cout << "Removing erroneous kminmers" << endl;
		//KminmerParser readParser(_filename_readMinimizers, _minimizerSize, _kminmerSize);
		//auto fp3 = std::bind(&Bloocoo::removeErroneousKminmers, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
		//readParser.parse(fp3);
	}
	*/
	


	//unordered_map<KmerVec, KminmerData> _kminmersData;
	//auto fp2 = std::bind(&Bloocoo::createMDBG_index, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
	//parser.parse(fp2);

	/*
	for(const auto& it : _kminmersData){

		const KmerVec& vec = it.first;
		const KminmerData& kminmerData = it.second;

		//cout << "A ENLEVER" << endl;
		//if(_kminmerSize == 4 && kminmerData._count == 1) continue;
		if(kminmerData._count == 1) continue;

		//if(kminmerCounts[vec] > 1000) cout << kminmerCounts[vec] << endl;
		//if(kminmerData._count <= 2){
			//if(_kminmersData[vec]._count <= minAbundance_cutoff) continue;
		//}

		//if(vec.isPalindrome()){
		//	cout << "Palindrome:" << _mdbg->_dbg_nodes[vec]._index << endl;
		//}

		//cout << kminmersData[vec] << endl;

		//KminmerData& kminmerData = _kminmersData[vec];
		//cout << kminmerData._length << endl;
		_mdbg->addNode(vec, kminmerData._length, kminmerData._overlapLength_start, kminmerData._overlapLength_end, kminmerData._isReversed);
		_mdbg->_dbg_nodes[vec]._abundance = kminmerData._count;

	}
	*/



	cout << "Nb solid kminmers: " << _mdbg->_dbg_nodes.size() << endl;
	
	_kminmerExist.clear();
	_minimizerCounts.clear();
	_kminmersData.clear();
	_contigIndex.clear();
	_contigAbundances.clear();

	//exit(1);
	/*
	_debug_nbMinimizers = 0;
	_minimizerParser = new MinimizerParser(_minimizerSize, _minimizerDensity);
	_file_readData = gzopen(_filename_readMinimizers.c_str(),"wb");
	
	auto fp = std::bind(&Bloocoo::createMDBG_read, this, std::placeholders::_1, std::placeholders::_2);
	ReadParser readParser(_inputFilename, false);
	readParser.parse(fp);

	*/

	/*
	if(_input_extractKminmers != ""){
		extract_kminmers();
		//return;
	}
	if(_inputDir != ""){
		//execute_binning_cleanGraph();
		execute_binning();
		return;
	}
	*/

	//unordered_map<MinimizerPair, vector<u_int32_t>> minimizerPairs_to_reads;

	//vector<string> read_headers;

	//float minimizerDensity = 0.005; // 0.004; 0.0008
	//u_int64_t maxHashValue = -1;
	//cout << maxHashValue << endl;
	//u_int64_t nbMinimizers = 0;

	//setDispatcher (new SerialDispatcher());



	//cout << _inputFilename << endl;

	//u_int32_t readIndex = 0;
	//u_int32_t datasetID = 0;

	//gzFile file_readComposition = gzopen(_filename_readCompositions.c_str(),"wb");
	
	//TCACATACTATCTGCACTGAC
	//KmerModel kmerModel(_minimizerSize);
	//MinimizerParser minimizerParser(_minimizerSize, _minimizerDensity);



	/*
	std::ifstream infile(_inputFilename.c_str());
	std::string line;
	int total_kminmers = 0;

	//double minimizerBound = maxHashValue * ((double)_minimizerDensity);
	u_int32_t nbMinimizers = 0;

	while (std::getline(infile, line))
	{
		cout << line << endl;

		gzFile fp;
		kseq_t *seq;
		int slen = 0, qlen = 0;
		fp = gzopen(line.c_str(), "r");
		seq = kseq_init(fp);

		while (kseq_read(seq) >= 0){
			
			if(readIndex % 100000 == 0) cout << readIndex << " " << nbMinimizers << endl;
 			//cout << readIndex << " " << strlen(seq->seq.s) << " " << nbMinimizers << endl;

			//cout << "--------------" << endl;
			//string sequence_str = "";
			//char lastChar = '0';

			//cout << strlen(seq->seq.s) << endl;

			u_int64_t pos = 0;
			//if(readIndex > 100000){
			//	cout << "HAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" << endl;
			//	break;
			//}

			//cout << "-" << endl;
			//Sequence& sequence = itSeq->item();

			//cout << sequence.toString() << endl;
			string rleSequence;
			vector<u_int64_t> rlePositions;

			//u_int64_t lastMinimizer = -1;
			vector<u_int64_t> minimizers;
			vector<u_int64_t> minimizers_pos;
			u_int64_t nbMinizersRead = 0;

			size_t nbMinimizersPerRead = 0;
			size_t nbHashsLala = 0;

			Encoder::encode_rle(seq->seq.s, strlen(seq->seq.s), rleSequence, rlePositions);
			//cout << rleSequence << endl;

			minimizerParser.parse(rleSequence, minimizers, minimizers_pos);
			nbMinimizers += minimizers.size();
			
			

			for(u_int64_t minimizer : minimizers){
				minimizerCounts[minimizer] += 1;
			}
			for(size_t i=0; i<rlePositions.size(); i++){
				rlePositions[i] = i;
			}

			vector<KmerVec> kminmers; 
			vector<ReadKminmer> kminmersInfo;
			MDBG::getKminmers(_minimizerSize, _kminmerSize, minimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, readIndex);

			bool isLala = false;
			bool isLoulou = false;

			for(size_t i=0; i<kminmers.size(); i++){

				if(kminmersData.find(kminmers[i]) == kminmersData.end()){
					kminmersData[kminmers[i]] = {0, kminmersInfo[i]._length - _minimizerSize, kminmersInfo[i]._seq_length_start, kminmersInfo[i]._seq_length_end, kminmersInfo[i]._isReversed};
				}

				kminmersData[kminmers[i]]._count += 1;

			}

			u_int16_t size = minimizers.size();
			gzwrite(file_readData, (const char*)&size, sizeof(size));
			gzwrite(file_readData, (const char*)&minimizers[0], size * sizeof(u_int64_t));
			

			readIndex += 1;
		}

		datasetID += 1;

	}
	*/
	
	//gzclose(_file_readData);
	
	//exit(1);
	//unordered_map<u_int64_t, u_int64_t> minimizerCounts;
	//unordered_map<KmerVec, KminmerData> kminmersData;


	/*
	IBank* inbank = Bank::open(_inputFilename);

	
	Iterator<Sequence>* itSeq = createIterator<Sequence> (
														inbank->iterator(),
														inbank->estimateNbItems(),
														"Parsing reads"
														);

	LOCAL (itSeq);
		
	std::vector<Iterator<Sequence>*> itBanks =  itSeq->getComposition();

	ModelCanonical model (_minimizerSize);
	ModelCanonical::Iterator itKmer (model);
	//vector<u_int32_t> neighbors;

	hash<KmerVec> h;
	//_node_id = 0;

	ReadData lala1;
	ReadData lala2;
	vector<u_int64_t> minim1;
	vector<u_int64_t> minim2;
	vector<ReadData> comps;

	u_int32_t nbMinimizers = 0;
	int lol = 0;

	for (size_t i=0; i<itBanks.size(); i++)
	{
		itSeq = createIterator<Sequence> (itBanks[i], inbank->estimateNbItemsBanki(i), "lala");

		for (itSeq->first(); !itSeq->isDone(); itSeq->next()){

			if(readIndex % 100000 == 0) cout << readIndex << " " << nbMinimizers << endl;
			u_int64_t pos = 0;
			//if(readIndex > 100000){
			//	cout << "HAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" << endl;
			//	break;
			//}

			//cout << "-" << endl;
			Sequence& sequence = itSeq->item();

			//cout << sequence.toString() << endl;
			string rleSequence;
			vector<u_int64_t> rlePositions;



			//u_int64_t lastMinimizer = -1;
			vector<u_int64_t> minimizers;
			vector<u_int64_t> minimizers_pos;
			u_int64_t nbMinizersRead = 0;

			size_t nbMinimizersPerRead = 0;
			size_t nbHashsLala = 0;

			Encoder::encode_rle(sequence.getDataBuffer(), sequence.getDataSize(), rleSequence, rlePositions);

			

			
			//for(size_t i=0; i<rlePositions.size(); i++) {
			//	cout << i << ": " << rlePositions[i] << " " << endl;
			//}

			//_readData.push_back(readData);

			//sequence.setData();
        	//char* sequence_rle = (char*)  MALLOC (bs->read->max);
			//cout << sequence.toString() << endl;
			//cout << sequence_str << endl;
			//cout << readseq[0] << endl;
			//cout << readseq[1] << endl;
			//cout << readseq[2] << endl;
			//cout << ((readseq[0]>>1)&3) << endl;
			//cout << ((readseq[1]>>1)&3) << endl;
			//cout << ((readseq[2]>>1)&3) << endl;

			//exit(1);
			
			//u_int64_t lastMinimizer = -1;
			Data buf((char*)rleSequence.c_str());
			//u_int32_t sequenceLength = rleSequence.size();

			//vector<float> composition;
			//_compositionManager->readToComposition(buf, sequenceLength, composition);

			 

			//gzwrite(file_readComposition, (const char*)&sequenceLength, sizeof(sequenceLength));
			//gzwrite(file_readComposition, (const char*)&composition[0], _compositionManager->_compositionVectorSize * sizeof(float));
			

			//read_headers.push_back(sequence.getComment()); //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! to add

			//cout << sequence.getComment() << endl;
			//cout << sequence.getDataSize() << endl;
			itKmer.setData (buf);

			//u_int64_t lastMinimizer = -1;
			//vector<u_int64_t> minimizers;
			//vector<u_int64_t> minimizers_pos;
			//u_int64_t nbMinizersRead = 0;

			//vector<MinimizerPair> minimizerPairs;
			
			//cout << "---------------" << endl;

			//int lala = 0;
			//cout << rleSequence << endl;
			u_int32_t lastMinimizerPos = -1;
			pos = 0;
			for (itKmer.first(); !itKmer.isDone(); itKmer.next()){

				kmer_type kmerMin = itKmer->value();

				//cout << itKmer->value().getVal() << endl;
				//cout << itKmer->value() << endl;
				//cout << pos << " " << (rleSequence.size()-_minimizerSize) << endl;

				if(pos == 0){
					//cout << "lala1" << endl;
					pos += 1;
					continue;
				}
				else if(pos == rleSequence.size()-_minimizerSize){
					//cout << "lala2" << endl;
					continue;
				}

				//if(!itKmer->value().isValid()) continue;
				//kmer_type kmerMin = min(itKmer->value(), revcomp(itKmer->value(), _kmerSize));
				//if(lala < 100 ) cout << model.toString(itKmer->value()) << endl;
				//lala += 1;
				u_int64_t kmerValue = kmerMin.getVal();
				u_int64_t minimizer;
				MurmurHash3_x64_128 ((const char*)&kmerValue, sizeof(kmerValue), _seed, &_hash_otpt);
				minimizer = _hash_otpt[0];



				//if(minimizerCounts[minimizer] > 1000) cout << minimizer << endl;
				double kmerHashed_norm = ((double) minimizer) / maxHashValue;
				if(kmerHashed_norm < _minimizerDensity){


					minimizers.push_back(minimizer);
					minimizers_pos.push_back(pos);

					
					minimizerCounts[minimizer] += 1;
					nbMinimizers += 1;

				}

				//cout << kmerHashed << endl;
				pos += 1;
			}
			

			//Small hack to disable rle positions mapping, the mdbg work on compressed positions
			for(size_t i=0; i<rlePositions.size(); i++){
				rlePositions[i] = i;
			}

			vector<KmerVec> kminmers; 
			vector<ReadKminmer> kminmersInfo;
			MDBG::getKminmers(_minimizerSize, _kminmerSize, minimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, readIndex);

			bool isLala = false;
			bool isLoulou = false;

			for(size_t i=0; i<kminmers.size(); i++){

				if(kminmersData.find(kminmers[i]) == kminmersData.end()){
					kminmersData[kminmers[i]] = {0, kminmersInfo[i]._length - _minimizerSize, kminmersInfo[i]._seq_length_start, kminmersInfo[i]._seq_length_end, kminmersInfo[i]._isReversed};
				}

				kminmersData[kminmers[i]]._count += 1;
				//kminmersData[kminmers[i]] = kminmersInfo[i]._length;

				//if(kminmers[i]._kmers[0] == 81904984297519754 && kminmers[i]._kmers[1] == 37284240752535023 && kminmers[i]._kmers[2] == 82851683898915704){
				//	isLala = true;
				//}
				//if(kminmers[i]._kmers[0] == 45618398143196555 && kminmers[i]._kmers[1] == 55406240415237664 && kminmers[i]._kmers[2] == 33230919625359447){
				//	isLoulou = true;
				//}




			}

			//if(isLala && isLoulou){
			//	lol += 1;
			//	cout << lol << endl;
			//}


			//cout << minimizers.size() << endl;
			u_int16_t size = minimizers.size();
			gzwrite(file_readData, (const char*)&size, sizeof(size));
			gzwrite(file_readData, (const char*)&minimizers[0], size * sizeof(u_int64_t));
			
			//comps.push_back({sequenceLength, composition});
		
			readIndex += 1;
		}


		datasetID += 1;
	}
	*/
	/*
	for(size_t i=0; i<comps.size(); i++){
		for(size_t j=i+1; j<comps.size(); j++){
			cout << computeDistanceTNF(comps[i], comps[j]) << endl;
		}
	}
	*/
	//gzclose(file_readComposition);


	/*
	file_readData = gzopen(_filename_readMinimizers.c_str(),"rb");




	//unordered_set<u_int64_t> filteredMinimizers;



	while(true){
		
		u_int16_t size;
		vector<u_int64_t> minimizers;
		gzread(file_readData, (char*)&size, sizeof(size));

		if(gzeof(file_readData)) break;
		
		minimizers.resize(size);
		gzread(file_readData, (char*)&minimizers[0], size * sizeof(u_int64_t));

		vector<u_int16_t> abundances_minimizers;
		for(u_int64_t minimizer : minimizers){
			abundances_minimizers.push_back(minimizerCounts[minimizer]);
			//cout << minimizerCounts[minimizers[i]] << endl;
		}

		//float median_abundance = compute_median(abundances_minimizers);
		//size_t minAbundance_cutoff = median_abundance / 10;
		//size_t cutoff = median_abundance * 8;

		//for(u_int64_t minimizer : minimizers){
		//	if(minimizerCounts[minimizer] > cutoff){
		//		filteredMinimizers.insert(minimizer);
		//	}
		//}

	}

	gzclose(file_readData);
	*/
	//cout << "Nb repeated minimizers: " << filteredMinimizers.size() << endl;

	/*
	gzFile file_filteredMinimiers = gzopen(_filename_filteredMinimizers.c_str(),"wb");
	for(u_int64_t minimizer : filteredMinimizers){
		gzwrite(file_filteredMinimiers, (char*)&minimizer, sizeof(minimizer));
		//cout << minimizer << endl;
	}
	gzclose(file_filteredMinimiers);
	*/


	/*
	file_readData = gzopen(_filename_readMinimizers.c_str(),"rb");
	unordered_map<KmerVec, u_int16_t> kminmerCounts;
	//unordered_map<KmerVec, u_int16_t> kminmerPreSuf_Counts;

	while(true){
		
		u_int16_t size;
		vector<u_int64_t> minimizers;
		gzread(file_readData, (char*)&size, sizeof(size));

		if(gzeof(file_readData)) break;
		
		minimizers.resize(size);
		gzread(file_readData, (char*)&minimizers[0], size * sizeof(u_int64_t));

		vector<KmerVec> kminmers; 
		MDBG::getKminmers(_kminmerSize, minimizers, kminmers);

		//cout << kminmers.size() << endl;
		for(KmerVec& vec : kminmers){
			kminmerCounts[vec] += 1;

			//kminmerPreSuf_Counts[vec.prefix()] += 1;
			//kminmerPreSuf_Counts[vec.suffix()] += 1;

			//cout << endl;
			//cout << vec.toString() << " " << kminmerCounts[vec] << endl;
		}

	}
	*/

	//gzclose(file_readData);


	/*
	cout << "Nb kminmers: " << _kminmersData.size() << endl;

	_file_readData = gzopen(_filename_readMinimizers.c_str(),"rb");
	//cout << "Nb edges: " << nbEdges << endl;

	//cout << "Cleaning errors..." << endl;

	//delete mdbg;
	
	_mdbg = new MDBG(_kminmerSize);

	//for(size_t readIndex=0; readIndex<_readData.size(); readIndex++) {
	//	ReadData& readData = _readData[readIndex];
	//	vector<u_int64_t>& minimizers = readData._minimizers;

	
	while(true){
		
		u_int16_t size;
		vector<u_int64_t> minimizers;
		gzread(_file_readData, (char*)&size, sizeof(size));

		if(gzeof(_file_readData)) break;
		
		minimizers.resize(size);
		gzread(_file_readData, (char*)&minimizers[0], size * sizeof(u_int64_t));


		vector<KmerVec> kminmers; 
		vector<ReadKminmer> kminmersInfo;
		vector<u_int64_t> minimizersPos; 
		vector<u_int64_t> rlePositions;
		MDBG::getKminmers(_minimizerSize, _kminmerSize, minimizers, minimizersPos, kminmers, kminmersInfo, rlePositions, 0);
		//getKminmers_filterRepeatedEdge(minimizers, filteredMinimizers, kminmers, kminmerCounts);


		vector<u_int32_t> abundances;

		for(KmerVec& vec : kminmers){
			abundances.push_back(_kminmersData[vec]._count);
		}

		float minAbundance_cutoff = compute_median(abundances) / 2;

		//cout << "---------------" << endl;
		//for(u_int64_t minimizer : minimizers){
		//	cout << minimizerCounts[minimizer] << endl;
		//}
		//for(KmerVec& vec : kminmers){
		//	cout << kminmerCounts[vec] << endl;
		//}

		for(KmerVec& vec : kminmers){

			if(_kminmersData[vec]._count == 1) continue; //ATTENTION, ON ENELEVE TOUS LES KMINMER VU UNE FOIS
			//if(kminmerCounts[vec] > 1000) cout << kminmerCounts[vec] << endl;
			if(_kminmersData[vec]._count <= 2){
				if(_kminmersData[vec]._count <= minAbundance_cutoff) continue;
			}


			//cout << kminmersData[vec] << endl;

			KminmerData& kminmerData = _kminmersData[vec];

			

			_mdbg->addNode(vec, kminmerData._length, kminmerData._overlapLength_start, kminmerData._overlapLength_end, kminmerData._isReversed);
			*/
			/*
			if(_mdbg->_dbg_nodes[vec]._index == 16494){
				cout << "--------" << endl;
				cout << vec._kmers[0] << endl;
				cout << vec._kmers[1] << endl;
				cout << vec._kmers[2] << endl;
				cout << _mdbg->_dbg_nodes[vec]._index << ": " << kminmerData._length << " " << kminmerData._overlapLength_start << " " << kminmerData._overlapLength_end << " " << kminmerData._isReversed << endl;
				//exit(1);
			}

			if(_mdbg->_dbg_nodes[vec]._index == 8581){
				cout << "--------" << endl;
				cout << vec._kmers[0] << endl;
				cout << vec._kmers[1] << endl;
				cout << vec._kmers[2] << endl;
				cout << _mdbg->_dbg_nodes[vec]._index << ": " << kminmerData._length << " " << kminmerData._overlapLength_start << " " << kminmerData._overlapLength_end << " " << kminmerData._isReversed << endl;
				//exit(1);
			}
			*/
			/*
			if(_mdbg->_dbg_nodes[vec]._index == 7837 || _mdbg->_dbg_nodes[vec]._index == 7836){
				cout << "--------" << endl;
				cout << vec._kmers[0] << endl;
				cout << vec._kmers[1] << endl;
				cout << vec._kmers[2] << endl;
			}
			*/

			/*
			if(vec.isPalindrome() && _mdbg->_dbg_nodes[vec]._abundance == 1){
				cout << "--------" << endl;
				cout << _mdbg->_dbg_nodes[vec]._index << endl;
				cout << vec._kmers[0] << endl;
				cout << vec._kmers[1] << endl;
				cout << vec._kmers[2] << endl;
			}
			*/
			/*
			if(mdbg->_dbg_nodes[vec]._index == 8814){
				cout << "8814:" << endl;
				cout << mdbg->_dbg_nodes[vec]._abundance << endl;
				cout << kminmerPreSuf_Counts[vec.prefix()] << endl;
				cout << kminmerPreSuf_Counts[vec.suffix()] << endl;
			}
			else if(mdbg->_dbg_nodes[vec]._index == 4961){
				cout << "4961:" << endl;
				cout << mdbg->_dbg_nodes[vec]._abundance << endl;
				cout << kminmerPreSuf_Counts[vec.prefix()] << endl;
				cout << kminmerPreSuf_Counts[vec.suffix()] << endl;
			}*/

			//mdbg->_dbg_nodes[vec]._abundance = kminmerCounts[vec]._abundance;
		/*
		}
*/
		/*
		vector<u_int16_t> abundances;

		//double abundance_mean = 0;
		for(KmerVec& kminmer : kminmers){
			//if(mdbg->_dbg_nodes.find(kminmer) == mdbg->_dbg_nodes.end()) cout << "haaaa" << endl;
			//abundance_mean += mdbg->_dbg_nodes[kminmer]._abundance;
			abundances.push_back(mdbg->_dbg_nodes[kminmer]._abundance);
		}
		//abundance_mean /= kminmers.size();

		float minAbundance_cutoff = median(abundances) / 10;
		*/

		

		//if(abundance_mean < 10){ 
		//cout << "----------------" << endl;
		//cout << "Mean: " << abundance_mean << endl;
		//cout << "Median: " << median(abundances) << endl;
		//for(KmerVec& kminmer : kminmers){
		//	cout << mdbg->_dbg_nodes[kminmer]._abundance << endl;
		//}

		//}
		//cout << "----------------" << endl;
		//cout << "Mean: " << abundance_mean << endl;
		
		/*
		bool lala = false;
		for(KmerVec& kminmer : kminmers){
			if(mdbg->_dbg_nodes[kminmer]._abundance > 500) lala = true;
		}

		//if(lala){
			//cout << "haaa " << minimizers.size() << endl;
			cout << "--------------" << endl;
			
			cout << "Cutoff min: " << minAbundance_cutoff << endl;
			cout << "Cutoff max: " << maxAbundance_cutoff << endl;

			for(size_t i=0; i<minimizers.size(); i++){
				cout << minimizerCounts[minimizers[i]] << endl;
			}
			cout << "-" << endl;
			for(KmerVec& kminmer : kminmers){
				cout << mdbg->_dbg_nodes[kminmer]._abundance << endl;
			}
		//}
		*/

		/*
		for(KmerVec& kminmer : kminmers){
			//if(mdbg->_dbg_nodes[kminmer]._abundance <= 2){
			if(mdbg->_dbg_nodes[kminmer]._abundance <= minAbundance_cutoff) continue;
			//}
			//if(mdbg->_dbg_nodes[kminmer]._abundance > 2000) continue;

			mdbg_errorFree->addNode(kminmer, mdbg->_dbg_nodes[kminmer]._length);
			mdbg_errorFree->_dbg_nodes[kminmer]._abundance = mdbg->_dbg_nodes[kminmer]._abundance;

		}
		*/

	/*
	}

	cout << "Nb solid kminmers: " << _mdbg->_dbg_nodes.size() << endl;
	
	_minimizerCounts.clear();
	_kminmersData.clear();*/
}

void Bloocoo::createGfa(){
	cout << "Writing gfa..." << endl;
	//cout << "Cleaning repeats..." << endl;

	//delete mdbg_repeatFree;
	//cout << _dbg_edges.size() << endl;

	u_int64_t nbEdges = 0;

	string gfa_filename = _outputDir + "/minimizer_graph.gfa";
	ofstream output_file_gfa(gfa_filename);

	for(auto vec_id : _mdbg->_dbg_nodes){

		KmerVec vec = vec_id.first;
		KmerVec vec_rev = vec_id.first.reverse();
		u_int32_t id = vec_id.second._index;

		output_file_gfa << "S" << "\t" << id << "\t" << "*" << "\t" << "LN:i:" << vec_id.second._length << "\t" << "dp:i:" << vec_id.second._abundance << endl;

		//cout << mdbg->_dbg_edges[vec.prefix().normalize()].size() << endl;
		for(KmerVec& v : _mdbg->_dbg_edges[vec.prefix().normalize()]){
			if(v==vec) continue;
			KmerVec v_rev = v.reverse();

			if (vec.suffix() == v.prefix()) {
				nbEdges += 1;
				u_int16_t overlapLength =  _mdbg->_dbg_nodes[v]._length - _mdbg->_dbg_nodes[v]._overlapLength_end; //   min(_mdbg->_dbg_nodes[v]._overlapLength_end, _mdbg->_dbg_nodes[vec]._overlapLength_end);
				output_file_gfa << "L" << "\t" << id << "\t" << "+" << "\t" << _mdbg->_dbg_nodes[v]._index << "\t" << "+" << "\t" << overlapLength << "M" << endl;
				//vec_add_edge("+", "+");
			}
			if (vec.suffix() == v_rev.prefix()) {
				nbEdges += 1;
				u_int16_t overlapLength = _mdbg->_dbg_nodes[v]._length - _mdbg->_dbg_nodes[v]._overlapLength_start; //min(_mdbg->_dbg_nodes[v]._overlapLength_end, _mdbg->_dbg_nodes[vec]._overlapLength_end);
				output_file_gfa << "L" << "\t" << id << "\t" << "+" << "\t" << _mdbg->_dbg_nodes[v]._index << "\t" << "-" << "\t" << overlapLength << "M" << endl;
				//vec_add_edge("+", "-");
			}
			if (vec_rev.suffix() == v.prefix()) {
				nbEdges += 1;
				u_int16_t overlapLength = _mdbg->_dbg_nodes[v]._length - _mdbg->_dbg_nodes[v]._overlapLength_end; //min(_mdbg->_dbg_nodes[v]._overlapLength_start, _mdbg->_dbg_nodes[vec]._overlapLength_start);
				output_file_gfa << "L" << "\t" << id << "\t" << "-" << "\t" << _mdbg->_dbg_nodes[v]._index << "\t" << "+" << "\t" << overlapLength << "M" << endl;
				//vec_add_edge("-", "+");
			}
			if (vec_rev.suffix() == v_rev.prefix()) {
				nbEdges += 1;
				u_int16_t overlapLength = _mdbg->_dbg_nodes[v]._length - _mdbg->_dbg_nodes[v]._overlapLength_start;; //min(_mdbg->_dbg_nodes[v]._overlapLength_start, _mdbg->_dbg_nodes[vec]._overlapLength_start);
				output_file_gfa << "L" << "\t" << id << "\t" << "-" << "\t" << _mdbg->_dbg_nodes[v]._index << "\t" << "-" << "\t" << overlapLength << "M" << endl;
				//vec_add_edge("-", "-");
			}

		}
		for(KmerVec& v : _mdbg->_dbg_edges[vec.suffix().normalize()]){
			if(v==vec) continue;
			KmerVec v_rev = v.reverse();

			if (vec.suffix() == v.prefix()) {
				nbEdges += 1;
				u_int16_t overlapLength = _mdbg->_dbg_nodes[v]._length - _mdbg->_dbg_nodes[v]._overlapLength_end; //min(_mdbg->_dbg_nodes[v]._overlapLength_end, _mdbg->_dbg_nodes[vec]._overlapLength_end);
				//cout << overlapLength << " " << _mdbg->_dbg_nodes[v]._length<< endl;
				output_file_gfa << "L" << "\t" << id << "\t" << "+" << "\t" << _mdbg->_dbg_nodes[v]._index << "\t" << "+" << "\t" << overlapLength << "M" << endl;
				//vec_add_edge("+", "+");
			}
			if (vec.suffix() == v_rev.prefix()) {
				nbEdges += 1;
				u_int16_t overlapLength = _mdbg->_dbg_nodes[v]._length - _mdbg->_dbg_nodes[v]._overlapLength_start; //min(_mdbg->_dbg_nodes[v]._overlapLength_end, _mdbg->_dbg_nodes[vec]._overlapLength_end);
				//cout << overlapLength << " " << _mdbg->_dbg_nodes[v]._length<< endl;
				output_file_gfa << "L" << "\t" << id << "\t" << "+" << "\t" << _mdbg->_dbg_nodes[v]._index << "\t" << "-" << "\t" << overlapLength << "M" << endl;
				//vec_add_edge("+", "-");
			}
			if (vec_rev.suffix() == v.prefix()) {
				nbEdges += 1;
				u_int16_t overlapLength = _mdbg->_dbg_nodes[v]._length - _mdbg->_dbg_nodes[v]._overlapLength_end; //min(_mdbg->_dbg_nodes[v]._overlapLength_start, _mdbg->_dbg_nodes[vec]._overlapLength_start);
				//cout << overlapLength << " " << _mdbg->_dbg_nodes[v]._length << endl;
				output_file_gfa << "L" << "\t" << id << "\t" << "-" << "\t" << _mdbg->_dbg_nodes[v]._index << "\t" << "+" << "\t" << overlapLength << "M" << endl;
				//vec_add_edge("-", "+");
			}
			if (vec_rev.suffix() == v_rev.prefix()) {
				nbEdges += 1;
				u_int16_t overlapLength = _mdbg->_dbg_nodes[v]._length - _mdbg->_dbg_nodes[v]._overlapLength_start; //min(_mdbg->_dbg_nodes[v]._overlapLength_start, _mdbg->_dbg_nodes[vec]._overlapLength_start);
				//cout << overlapLength << " " << _mdbg->_dbg_nodes[v]._length<< endl;
				output_file_gfa << "L" << "\t" << id << "\t" << "-" << "\t" << _mdbg->_dbg_nodes[v]._index << "\t" << "-" << "\t" << overlapLength << "M" << endl;
				//vec_add_edge("-", "-");
			}
		}
		//for(KmerVec& v : _dbg_edges_prefix[vec.suffix().normalize()]){
		//	output_file_gfa << "L" << "\t" << vec_id.second << "\t" << "+" << "\t" << _dbg_nodes[v] << "\t" << "+" << "\t" << "0M" << endl;
		//}
	}

	output_file_gfa.close();
	
	cout << "Nb nodes: " << _mdbg->_dbg_nodes.size() << endl;
	cout << "Nb edges: " << nbEdges << endl;
	











	/*
	GfaParser gfaParser;
	AdjGraph* graph = gfaParser.createGraph_lol(gfa_filename);
	cout << "haaaaa: " << graph->_nbNodes << endl;

	unordered_set<u_int32_t> abundant_kminmers;

	int lala = 0;
	for(auto it : mdbg->_dbg_nodes){
		if(it.second._abundance > 300){

			abundant_kminmers.insert(it.second._index);
			//cout << "------------ " << it.second._index << endl;
			//cout << it.second._abundance << endl;
			
			vector<u_int32_t> neighbors;
			graph->collectNeighbors(it.second._index, 5, neighbors);
			for(u_int32_t nn : neighbors){
				//cout << nn << endl;
				abundant_kminmers.insert(nn);
			}
			//cout << "n: " << neighbors.size() << endl;
			//cout << n << endl;
			lala += 1;
		}
	}
	cout << "Nb nodes abundant: " << lala << endl;
	gfaParser.rewriteGfa_withNodes(gfa_filename, gfa_filename + "_abundant.gfa", abundant_kminmers);

	*/





	if(_isFirstPass){
		_mdbg->dump(_outputDir + "/mdbg_nodes_init.gz");
	}
	
	_mdbg->dump(_outputDir + "/mdbg_nodes.gz");
	//_mdbg->_dbg_nodes.clear();
	//_mdbg->_dbg_edges.clear();



	//createGroundTruth();
}



void Bloocoo::createGroundTruth(){
}

/*
void Bloocoo::createGroundTruth(){


	string gfa_filename = _outputDir + "/minimizer_graph_notips_nobubbles.gfa";
	cout << gfa_filename << endl;
	gzFile file_readData = gzopen(_filename_readMinimizers.c_str(),"rb");

	ofstream file_groundTruth(_outputDir + "/groundtruth.csv");
	file_groundTruth << "Name,Colour" << endl;

	MDBG* mdbg = new MDBG(_kminmerSize);
	mdbg->load(_outputDir + "/mdbg_nodes.gz");

	vector<u_int32_t> unitigLengths;
	vector<int32_t> node_to_unitig(mdbg->_dbg_nodes.size(), 0);
	GfaParser gfaParser;
	AdjGraph* graph = gfaParser.createGraph(gfa_filename, node_to_unitig, unitigLengths);


	ReadIndexType readIndex = 0;

	unordered_set<u_int32_t> processedUnitigs;



	unordered_set<u_int64_t> filteredMinimizers;


	while(true){

		u_int16_t size;
		vector<u_int64_t> minimizers;
		gzread(file_readData, (char*)&size, sizeof(size));

		if(gzeof(file_readData)) break;
		
		minimizers.resize(size);
		gzread(file_readData, (char*)&minimizers[0], size * sizeof(u_int64_t));

		
		vector<KmerVec> kminmers; 
		getKminmers(minimizers, filteredMinimizers, kminmers);

		vector<u_int32_t> datasets;

		for(KmerVec& vec : kminmers){
			if(mdbg->_dbg_nodes.find(vec) == mdbg->_dbg_nodes.end()) continue;
			
			kminmers.push_back(vec);
			//kminmers.push_back(vec.normalize());

			u_int32_t dataset = _evaluation_readToDataset[readIndex];
			if(std::find(datasets.begin(), datasets.end(), dataset) != datasets.end()) continue;
			
			datasets.push_back(dataset);
		}


		for(KmerVec& vec : kminmers){

			u_int32_t kminmer_index = mdbg->_dbg_nodes[vec]._index;
			u_int32_t unitigIndex = node_to_unitig[kminmer_index];

			if(processedUnitigs.find(unitigIndex) != processedUnitigs.end()) continue;

			processedUnitigs.insert(unitigIndex);

			//if(unitigIndex == 14) cout << "allo" << endl;
			//cout << kminmer_index << " " << unitigIndex << endl;

			string unitigName = "utg";
			string unitig_name_id = to_string(unitigIndex+1);
			size_t nbZeros = 7 - unitig_name_id.size();
			//cout << unitigIndex << " " << nbZeros << endl;
			for(size_t i=0; i<nbZeros; i++){
				unitigName += "0";
			}
			unitigName += unitig_name_id + "l";

			if(datasets.size() == 0){
				file_groundTruth << unitigName << "," << -2 << endl;
			}
			else if(datasets.size() == 1){
				file_groundTruth << unitigName << "," << datasets[0] << endl;
			}
			else{
				file_groundTruth << unitigName << "," << -1 << endl;
			}

			//size_t nbZeros = utg0000001l
		}

		
		readIndex += 1;
	}

	file_groundTruth.close();
	delete mdbg;
}
*/


	/*
	vector<vector<u_int32_t>> components;
	
	graph->computeConnectedComponents(components);
	
	cout << endl << "Nb connected components: " << components.size() << endl;
	for(size_t i=0; i<components.size(); i++){

		u_int64_t component_size_nt = 0;
		for(u_int32_t utg : components[i]){
			component_size_nt += unitigLengths[utg];
		}
		if(component_size_nt > 1000000){
			cout << i <<": " << component_size_nt << endl;
		}
	}
	cout << endl;
	
	unordered_set<u_int64_t> filteredMinimizers;
	gzFile file_filteredMinimiers = gzopen(_filename_filteredMinimizers.c_str(),"rb");

	while(true){
		u_int64_t minimizer;
		gzread(file_filteredMinimiers, (char*)&minimizer, sizeof(minimizer));
		if(gzeof(file_filteredMinimiers)) break;
		filteredMinimizers.insert(minimizer);
		//cout << minimizer << endl;
	}

	gzclose(file_filteredMinimiers);
	cout << "Nb filtered minimizers: " << filteredMinimizers.size() << endl;


	

	gzFile file_readData = gzopen(_filename_readMinimizers.c_str(),"rb");

	ReadIndexType readIndex = 0;
	vector<UnitigData> unitigDatas;
	for(size_t i=0; i<graph->_nbNodes; i++){
		unitigDatas.push_back({i, {}});
	}

	while(true){
		
		//cout << readIndex << endl;

		u_int16_t size;
		vector<u_int64_t> minimizers;
		gzread(file_readData, (char*)&size, sizeof(size));

		if(gzeof(file_readData)) break;
		
		minimizers.resize(size);
		gzread(file_readData, (char*)&minimizers[0], size * sizeof(u_int64_t));


		vector<KmerVec> kminmers; 
		getKminmers(minimizers, filteredMinimizers, kminmers);



		vector<ReadIndexType> unitigIndexex;

		for(KmerVec& vec : kminmers){
			if(mdbg->_dbg_nodes.find(vec) == mdbg->_dbg_nodes.end()) continue;

			u_int32_t kminmer_index = mdbg->_dbg_nodes[vec]._index;
			u_int32_t unitigIndex = node_to_unitig[kminmer_index];
			if(unitigIndex == -1) continue;

			//cout << kminmer_index << " " << unitigIndex << endl;
			if(std::find(unitigIndexex.begin(), unitigIndexex.end(), unitigIndex) != unitigIndexex.end()) continue;

			unitigIndexex.push_back(unitigIndex);
		}

		if(unitigIndexex.size() <= 1) continue;

		for(KmerVec& vec : kminmers){
			if(mdbg->_dbg_nodes.find(vec) == mdbg->_dbg_nodes.end()) continue;

			u_int32_t kminmer_index = mdbg->_dbg_nodes[vec]._index;
			u_int32_t unitigIndex = node_to_unitig[kminmer_index];
			if(unitigIndex == -1) continue;

			UnitigData& unitigData = unitigDatas[unitigIndex];
			
			if(unitigData._readIndexes_exists.find(readIndex) != unitigData._readIndexes_exists.end()) continue;
			//cout << unitigData._readIndexes.size() << endl;
			//if(std::find(unitigData._readIndexes.begin(), unitigData._readIndexes.end(), readIndex) != unitigData._readIndexes.end()) continue;
				
			unitigData._readIndexes_exists.insert(readIndex);
			unitigData._readIndexes.push_back(readIndex);

			//cout << kminmer_index << " " << unitigIndex << " " << unitigData._readIndexes.size() << " " << readIndex << endl;
			//cout << unitigData._readIndexes.size() << " " << unitigLengths[unitigIndex] << endl;
			//if(std::find(unitigIndexex.begin(), unitigIndexex.end(), unitigIndex) != unitigIndexex.end()) continue
			//unitigIndexex.push_back(unitigIndex);

		}
		
		readIndex += 1;
	}






	
	unordered_set<DbgEdge, hash_pair> unsupportedEdges;

	for(size_t utg=0; utg<graph->_nbNodes; utg++){

		//cout << utg << " " << graph->_nbNodes << endl;

		adjNode* node = graph->_nodes[utg];

        while (node != nullptr) {
			
			ReadIndexType utg_n = node->val;

			//cout << utg << " " << utg_n << " " << unitigDatas[utg]._readIndexes.size() << " " << unitigDatas[utg_n]._readIndexes.size() << endl;

			if(unsupportedEdges.find({utg, utg_n}) != unsupportedEdges.end() || unsupportedEdges.find({utg_n, utg}) != unsupportedEdges.end()) {	
				node = node->next;
				continue;
			}

			if(shareAnyRead(unitigDatas[utg], unitigDatas[utg_n])){
				node = node->next;
				continue;
			}

			//if(shareAnyRead(unitigDatas[utg], unitigDatas[utg_n])){
			unsupportedEdges.insert({utg, utg_n});
			unsupportedEdges.insert({utg_n, utg});
			//}

			node = node->next;
        }

	}

	cout << "Nb unsupported edges: " << unsupportedEdges.size() << endl;

	gfaParser.rewriteGfa_withoutEdges(gfa_filename, gfa_filename +"_2.gfa", unsupportedEdges);



	delete graph;
	node_to_unitig.clear();
	unitigLengths.clear();
	components.clear();
	graph = gfaParser.createGraph(gfa_filename +"_2.gfa", node_to_unitig, unitigLengths);
	
	cout << "Nb nodes: " << graph->_nbNodes << endl;
	cout << "Nb edges: " << graph->_nbEdges << endl;

	graph->computeConnectedComponents(components);

	cout << endl << "Nb connected components: " << components.size() << endl;
	for(size_t i=0; i<components.size(); i++){

		u_int64_t component_size_nt = 0;
		for(u_int32_t utg : components[i]){
			component_size_nt += unitigLengths[utg];
		}
		if(component_size_nt > 1000000){
			cout << i <<": " << component_size_nt << endl;
		}
	}
	cout << endl;

	*/
	/*
	ofstream file_groundTruth(_outputDir + "/binning_results.csv");
	file_groundTruth << "Name,Colour" << endl;

	int iter = 0;
	//utg0000256l
	ReadIndexType utg = 33; //255
	stack<ReadIndexType> stack;
	stack.push(utg);

	unordered_set<DbgEdge, hash_pair> isEdgeVisited;

	while (!stack.empty() && iter < 500){
		
		utg = stack.top();
		stack.pop();


		string unitigName = "utg";
		string unitig_name_id = to_string(utg+1);
		size_t nbZeros = 7 - unitig_name_id.size();
		//cout << unitigIndex << " " << nbZeros << endl;
		for(size_t i=0; i<nbZeros; i++){
			unitigName += "0";
		}
		unitigName += unitig_name_id + "l";

		cout << "Visit node: " << unitigName << endl;

		//cout << unitigName << endl;
		file_groundTruth << unitigName << "," << "red" << endl;




		adjNode* node = graph->_nodes[utg];
		vector<UnitigEdgeScore> unitigScores;

        while (node != nullptr) {
			
			ReadIndexType utg_n = node->val;
			if(isEdgeVisited.find({utg, utg_n}) != isEdgeVisited.end() || isEdgeVisited.find({utg_n, utg}) != isEdgeVisited.end()) {	
				node = node->next;
				continue;
			}

			unitigScores.push_back({utg, utg_n, computeSharedReads(unitigDatas[utg], unitigDatas[utg_n])});

			node = node->next;
        }

		if(unitigScores.size() == 0) continue;

		std::sort(unitigScores.begin(), unitigScores.end(), UnitigEdgeScoreComparator);

		UnitigEdgeScore& bestEdge = unitigScores[0];

		isEdgeVisited.insert({bestEdge._from, bestEdge._to});
		isEdgeVisited.insert({bestEdge._to, bestEdge._from});

		for(size_t i=0; i<unitigScores.size(); i++) {
			cout << unitigScores[i]._from << " " << unitigScores[i]._to << " " <<  unitigScores[i]._score << endl;
			if(unitigScores[i]._score == 0) continue;
			stack.push(unitigScores[i]._to);
		}

		//if(bestEdge._score == 0) continue;



		iter += 1;


	}
	

	file_groundTruth.close();
	*/
	/*
	for(UnitigData& unitigData : unitigDatas){
		cout << "------------- " << unitigData._index << endl;
		for(ReadIndexType readIndex : unitigData._readIndexes){
			cout << readIndex << " ";
		}
		cout << endl;
	}
	*/
//}
