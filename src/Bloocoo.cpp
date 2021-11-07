

#include <Bloocoo.hpp>


Bloocoo::Bloocoo () : Tool("bloocoo")
{


	getParser()->push_back (new OptionOneParam (STR_INPUT, "input file", true));
	getParser()->push_back (new OptionOneParam (STR_OUTPUT, "output dir", true));
	getParser()->push_back (new OptionOneParam (STR_MINIM_SIZE, "minimizer length", false, "21"));
	getParser()->push_back (new OptionOneParam (STR_KMINMER_SIZE, "k-min-mer length", false, "3"));
	getParser()->push_back (new OptionOneParam (STR_DENSITY, "density of minimizers", false, "0.005"));
	//getParser()->push_back (new OptionOneParam (STR_INPUT_DIR, "input dir", false, ""));
	//getParser()->push_back (new OptionOneParam (STR_INPUT_EXTRACT_KMINMERS, "", false, ""));

}





void Bloocoo::execute (){
	parseArgs();
	createMDBG();
	createGfa();
}

void Bloocoo::parseArgs(){

	//_compositionManager = new CompositionManager(4);

	//int w = 80;
	_kminmerSize = getInput()->getInt(STR_KMINMER_SIZE);
	_inputFilename = getInput()->getStr(STR_INPUT);
	_minimizerSize = getInput()->getInt(STR_MINIM_SIZE);
	_outputDir = getInput()->getStr(STR_OUTPUT);
	_minimizerDensity = getInput()->getDouble(STR_DENSITY);

	if(!System::file().doesExist (_outputDir)) System::file().mkdir(_outputDir, -1);

	cout << endl;
	cout << "Input: " << _inputFilename << endl;
	cout << "Output dir: " << _outputDir << endl;
	cout << "Minimizer length: " << _minimizerSize << endl;
	cout << "Kminmer length: " << _kminmerSize << endl;
	cout << "Density: " << _minimizerDensity << endl;
	cout << endl;

	//_inputDir = getInput()->get(STR_INPUT_DIR) ? getInput()->getStr(STR_INPUT_DIR) : "";
	//_input_extractKminmers= getInput()->get(STR_INPUT_EXTRACT_KMINMERS) ? getInput()->getStr(STR_INPUT_EXTRACT_KMINMERS) : "";

	_filename_readMinimizers = _outputDir + "/read_data.gz";
	//_filename_readCompositions = _outputDir + "/read_compositions.gz";
	//_filename_filteredMinimizers = _outputDir + "/filteredMinimizers.gz";
	//_filename_hifiasmGroundtruth = _outputDir + "/hifiasmGroundtruth.gz";


	string filename_parameters = _outputDir + "/parameters.gz";
	gzFile file_parameters = gzopen(filename_parameters.c_str(),"wb");
	gzwrite(file_parameters, (const char*)&_minimizerSize, sizeof(_minimizerSize));
	gzwrite(file_parameters, (const char*)&_kminmerSize, sizeof(_kminmerSize));
	gzwrite(file_parameters, (const char*)&_minimizerDensity, sizeof(_minimizerDensity));
	gzclose(file_parameters);

}

void Bloocoo::createMDBG (){



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

	unordered_map<MinimizerPair, vector<u_int32_t>> minimizerPairs_to_reads;

	vector<string> read_headers;

	//float minimizerDensity = 0.005; // 0.004; 0.0008
	u_int64_t maxHashValue = -1;
	//cout << maxHashValue << endl;
	//u_int64_t nbMinimizers = 0;

	u_int64_t _hash_otpt[2];
	int _seed = 42;
	setDispatcher (new SerialDispatcher());



	cout << _inputFilename << endl;
	IBank* inbank = Bank::open(_inputFilename);

	
	Iterator<Sequence>* itSeq = createIterator<Sequence> (
														inbank->iterator(),
														inbank->estimateNbItems(),
														"Parsing reads"
														);

	LOCAL (itSeq);
		
	std::vector<Iterator<Sequence>*> itBanks =  itSeq->getComposition();
	u_int32_t readIndex = 0;
	u_int32_t datasetID = 0;

	gzFile file_readData = gzopen(_filename_readMinimizers.c_str(),"wb");
	//gzFile file_readComposition = gzopen(_filename_readCompositions.c_str(),"wb");
	
	/*
	ModelCanonical model (_kmerSize);
	ModelCanonical::Iterator itKmer (model);

	//ModelMinimizer model (63, 14);
    //const ModelCanonical& modelMinimizer = model.getMmersModel();
	vector<u_int64_t> kmers;
	vector<u_int64_t> minimizers_0;
	vector<u_int64_t> minimizers_1;
	vector<u_int64_t> minimizers_2;

	//unordered_map<u_int64_t, u_int64_t> minimCount;
	for (size_t i=0; i<itBanks.size(); i++)
	{
		itSeq = createIterator<Sequence> (itBanks[i], inbank->estimateNbItemsBanki(i), "lala");

		for (itSeq->first(); !itSeq->isDone(); itSeq->next()){


			Sequence& sequence = itSeq->item();

			ReadData readData = {sequence.getDataSize(), {}, NULL, false};
			_compositionManager->readToComposition(sequence, readData._composition);
			_readData.push_back(readData);

			//set<u_int64_t> minimizers_0_set;
			kmers.clear();
			minimizers_0.clear();
			minimizers_1.clear();
			minimizers_2.clear();
			//cout << "------------" << endl;
			read_headers.push_back(sequence.getComment());

			char* readseq = sequence.getDataBuffer();
			string sequence_str;

			char lastChar = '0';
			for(size_t i=0; i<sequence.getDataSize(); i++){
				if(readseq[i] == lastChar) continue;
				sequence_str += readseq[i];
				lastChar = readseq[i];
			}



			//exit(1);
			size_t nbMinimizersPerRead = 0;

			//u_int64_t lastMinimizer = -1;
			Data buf((char*)sequence_str.c_str());
			

			itKmer.setData (buf);

			vector<MinimizerPair> minimizerPairs;

			//cout << "----------" << endl;
			u_int64_t pos = 0;
			for (itKmer.first(); !itKmer.isDone(); itKmer.next()){

				kmer_type kmerMin = itKmer->value();
				//if(!itKmer->value().isValid()) continue;
				//kmer_type kmerMin = min(itKmer->value(), revcomp(itKmer->value(), _kmerSize));
				//if(lala < 100 ) cout << model.toString(itKmer->value()) << endl;
				//lala += 1;
				u_int64_t kmerValue = kmerMin.getVal();
				u_int64_t kmerHashed;
				MurmurHash3_x64_128 ((const char*)&kmerValue, sizeof(kmerValue), _seed, &_hash_otpt);
				kmerHashed = _hash_otpt[0];

				//cout << kmerHashed << endl;
				kmers.push_back(kmerHashed);
				//minimizers_0_set.insert(kmerHashed);

				if(pos >= w){
					u_int64_t minimizer = -1; //(*minimizers_0_set.begin());
					for(u_int64_t kmer : kmers){
						if(kmer < minimizer){
							minimizer = kmer;
						}
					}
					//cout << kmers[0] << endl;
					//auto it = minimizers_0_set.find(kmers[0]);
					//minimizers_0_set.erase(it);

					kmers.erase(kmers.begin());

					if(minimizers_0.size() == 0 || minimizer != minimizers_0[minimizers_0.size()-1]){
						minimizers_0.push_back(minimizer);
					}
					//cout << "Min: " << (*minimizers_0_set.begin()) << endl;
				}

				pos += 1;
			}

			//cout << minimizers_0.size() << endl;
			for(size_t w=0; w<minimizers_0.size()-windowSize; w++){
				u_int64_t m = -1;
				for(size_t i=0; i<windowSize; i++){
					if(minimizers_0[w+i] < m){
						m = minimizers_0[w+i];
					}
				}
				if(std::find(minimizers_1.begin(), minimizers_1.end(), m) == minimizers_1.end()){
					minimizers_1.push_back(m);
				}
			}

			//cout << minimizers_1.size() << endl;

			for(size_t w=0; w<minimizers_1.size()-windowSize; w++){
				u_int64_t m = -1;
				for(size_t i=0; i<windowSize; i++){
					if(minimizers_1[w+i] < m){
						m = minimizers_1[w+i];
					}
				}
				if(std::find(minimizers_2.begin(), minimizers_2.end(), m) == minimizers_2.end()){
					minimizers_2.push_back(m);
				}
			}

			cout << minimizers_2.size() << endl;

			u_int64_t lastMinimizer = -1;
			for(u_int64_t minimizer : minimizers_2){



				if(lastMinimizer != -1){

					//nbMinizersRead += 1;

					MinimizerPair minimizerPair1 = {lastMinimizer, minimizer};
					MinimizerPair minimizerPair2 = {minimizer, lastMinimizer};

					bool pair1_exists = minimizerPairs_to_reads.find(minimizerPair1) != minimizerPairs_to_reads.end();
					bool pair2_exists = minimizerPairs_to_reads.find(minimizerPair2) != minimizerPairs_to_reads.end();
					if(pair1_exists){
						minimizerPairs.push_back(minimizerPair1);
						//minimizerPairs_to_reads[minimizerPair1].push_back(readIndex);
					}
					else if(pair2_exists){
						minimizerPairs.push_back(minimizerPair2);
						//minimizerPairs_to_reads[minimizerPair2].push_back(readIndex);
					}
					else{
						minimizerPairs.push_back(minimizerPair1);
						//minimizerPairs_to_reads[minimizerPair1] = vector<u_int32_t>();
						//minimizerPairs_to_reads[minimizerPair1].push_back(readIndex);
					}

					//if(minimCount.find(minimizer) == minimCount.end()){
					//	minimCount[minimizer] = 0;
					//}

					//minimCount[minimizer] += 1;

					//if(minimCount[minimizer] > 1000){
						//cout << model.toString(itKmer->value()) << endl;
					//}



				}
				
				lastMinimizer = minimizer;
			}






			//---------------------------------------------------------------------------------
			//cout << minimizerPairs.size() << endl;
			unordered_map<u_int64_t, u_int64_t> minimizerPairCount;
			for(MinimizerPair& minimizerPair : minimizerPairs){
				if(minimizerPairs_to_reads.find(minimizerPair) == minimizerPairs_to_reads.end()) continue;
				
				
				for(u_int64_t r : minimizerPairs_to_reads[minimizerPair]) {

					if(minimizerPairCount.find(r) == minimizerPairCount.end()){
						minimizerPairCount[r] = 0;
					}
					minimizerPairCount[r] += 1;
				}
			}

			//cout << minimizerPairCount.size() << endl;

			//cout << "----------------" << endl;
			//for (auto& it: minimizerPairCount) {
			//	cout << it.first << " " << it.second << endl;
			//}

			std::vector<std::pair<u_int64_t, u_int64_t>> elems(minimizerPairCount.begin(), minimizerPairCount.end());
			std::sort(elems.begin(), elems.end(), Ralalalala);
			
			//cout << "-" << endl;

			u_int64_t readIndex_merged = -1;
			for (size_t i=0; i<elems.size(); i++) {
				

				u_int64_t read = elems[i].first;
				float dist = computeDistanceTNF(_readData[readIndex], _readData[read]);

				
				//cout << elems[i].first << " " << elems[i].second << " " << dist << endl;

				if(dist < 0.002){ //0.05
					//readIndex_merged = read; //! disabled
					//break;
				}
			}

			if(readIndex_merged == -1){
				for(MinimizerPair& minimizerPair : minimizerPairs){
					minimizerPairs_to_reads[minimizerPair].push_back(readIndex);
				}
			}
			else{
				for(MinimizerPair& minimizerPair : minimizerPairs){

					if(minimizerPairs_to_reads.find(minimizerPair) == minimizerPairs_to_reads.end()){
						minimizerPairs_to_reads[minimizerPair].push_back(readIndex_merged);
					}
					else{
						vector<u_int32_t>& minimizerPairReads = minimizerPairs_to_reads[minimizerPair];
						if(std::find(minimizerPairReads.begin(), minimizerPairReads.end(), readIndex_merged) == minimizerPairReads.end()){
							minimizerPairReads.push_back(readIndex_merged);
						}
					}


				}
				_readData[readIndex]._composition.clear();
			}

			//cout << minimizerPairs.size() << endl;
			//if(read_headers[readIndex] == "G3_0 S1_1"){
				//cout << minimizerPairs.size() << " " << _readData[readIndex]._length << endl;
			//}



			//---------------------------------------------------------------------------------



			read_to_dataset.push_back(datasetID);
			readIndex += 1;
		}
		
		datasetID += 1;
	}
	*/



	//_minimizerPairMap = new MinimizerPairMap();
	//_overlapGraph = new AdjGraph();
	
	//unordered_map<u_int64_t, u_int64_t> minimCount;
	//set<MinimizerPair> _minimizerPairs;
	//vector<MinimizerPair_Edge> _minimizerPairs_edges;

	//unordered_map<MinimizerPair, u_int16_t> minimizerCounts;
	unordered_map<u_int64_t, u_int64_t> minimizerCounts;
	//unordered_map<KmerVec, u_int32_t> kminmerCounts;
	unordered_map<KmerVec, KminmerData> kminmersData;



	ModelCanonical model (_minimizerSize);
	ModelCanonical::Iterator itKmer (model);
	vector<u_int32_t> neighbors;

	hash<KmerVec> h;
	//_node_id = 0;

	ReadData lala1;
	ReadData lala2;
	vector<u_int64_t> minim1;
	vector<u_int64_t> minim2;
	vector<ReadData> comps;

	int lol = 0;

	for (size_t i=0; i<itBanks.size(); i++)
	{
		itSeq = createIterator<Sequence> (itBanks[i], inbank->estimateNbItemsBanki(i), "lala");

		for (itSeq->first(); !itSeq->isDone(); itSeq->next()){

			//if(readIndex > 100000){
			//	cout << "HAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" << endl;
			//	break;
			//}

			//cout << "-" << endl;
			Sequence& sequence = itSeq->item();

			//cout << sequence.toString() << endl;
			string rleSequence;
			vector<u_int64_t> rlePositions;

			/*
			string testseq = "AAAGGTCGTTTTA";
			encode_rle(testseq.c_str(), testseq.size(), rleSequence, rlePositions);
			cout << testseq << endl;
			cout << rleSequence << endl; 
			for(size_t i=0; i<rlePositions.size(); i++){
				cout << rlePositions[i] << " ";
			}
			cout << endl;
			*/

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
			size_t nbMinimizersPerRead = 0;

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

			u_int64_t lastMinimizer = -1;
			vector<u_int64_t> minimizers;
			vector<u_int64_t> minimizers_pos;
			u_int64_t nbMinizersRead = 0;

			vector<MinimizerPair> minimizerPairs;
			
			//cout << "---------------" << endl;

			//int lala = 0;
			u_int64_t pos = 0;
			u_int32_t lastMinimizerPos = -1;
			for (itKmer.first(); !itKmer.isDone(); itKmer.next()){

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

				kmer_type kmerMin = itKmer->value();
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

					/*
					if(lastMinimizerPos == -1){
						lastMinimizerPos = pos;
					}
					else{
						u_int16_t distance = pos - lastMinimizerPos;
						if(distance < 100) continue;
						lastMinimizerPos = pos;
					}*/
					minimizers.push_back(minimizer);
					minimizers_pos.push_back(pos);

					//cout << pos << endl;
					//cout << rlePositions[pos] << endl; 
					
					minimizerCounts[minimizer] += 1;
					
					//cout << pos << endl;
					//_readData[readIndex]._minimizers.push_back(minimizer);
					//readData << minimizer << endl;

					//if(readIndex == 0){
						//cout << minimizer << endl;
						//gzwrite(file_readData, (char*)&minimizer, 8);
						//cout << () << endl;
						//cout << readData.str().size() << endl;
					//}
					//nbMinimizers += 1;
					//minimizers.push_back(kmerHashed);
					/*
					if(lastMinimizer != -1){
						nbMinizersRead += 1;

						MinimizerPair minimizerPair1 = {lastMinimizer, minimizer};
						MinimizerPair minimizerPair2 = {minimizer, lastMinimizer};

						bool pair1_exists = minimizerPairs_to_reads.find(minimizerPair1) != minimizerPairs_to_reads.end();
						bool pair2_exists = minimizerPairs_to_reads.find(minimizerPair2) != minimizerPairs_to_reads.end();
						if(pair1_exists){

							//u_int16_t minimizerCount = minimizerCounts[minimizerPair1];
							//if(minimizerCount > 10000){
							//	continue;
							//}
							//else{
							//	minimizerCounts[minimizerPair1] += 1;
							//}

							minimizerPairs.push_back(minimizerPair1);
							//minimizerPairs_to_reads[minimizerPair1].push_back(readIndex);
						}
						else if(pair2_exists){

							//u_int16_t minimizerCount = minimizerCounts[minimizerPair2];
							//if(minimizerCount > 10000){
							//	continue;
							//}
							//else{
							//	minimizerCounts[minimizerPair2] += 1;
							//}

							minimizerPairs.push_back(minimizerPair2);
							//minimizerPairs_to_reads[minimizerPair2].push_back(readIndex);
						}
						else{

							//u_int16_t minimizerCount = minimizerCounts[minimizerPair1];
							//if(minimizerCount > 10000){
							//	continue;
							//}
							//else{
							//	minimizerCounts[minimizerPair2] += 1;
							//}

							minimizerPairs.push_back(minimizerPair1);
							//minimizerPairs_to_reads[minimizerPair1] = vector<u_int32_t>();
							//minimizerPairs_to_reads[minimizerPair1].push_back(readIndex);
						}

						//if(minimCount.find(minimizer) == minimCount.end()){
						//	minimCount[minimizer] = 0;
						//}

						//minimCount[minimizer] += 1;

						//if(minimCount[minimizer] > 1000){
							//cout << model.toString(itKmer->value()) << endl;
						//}
					}

					lastMinimizer = minimizer;
					*/
				}

				//cout << kmerHashed << endl;
				pos += 1;
			}

			//cout << endl;
			//for(u_int64_t m : minimizers){
			//	cout << m << endl;
			//}
			//cout << endl;
			/*
			vector<bool> bannedPositions(minimizers.size(), false);

			//cout << "-------" << endl;
			while(true){

				//cout << "w" << endl;

				bool hasPalindrome = false;

				int i_max = ((int)minimizers.size()) - (int)_kminmerSize + 1;
				for(int i=0; i<i_max; i++){
					if(bannedPositions[i]) continue;

					//cout << i << "/" << minimizers.size() << endl;
					KmerVec vec;

					vector<u_int32_t> currentMinimizerIndex;
					int j=i;
					while(true){
						
						if(j >= minimizers.size()){
							break;
						}
						if(bannedPositions[j]){
							//cout << "    is banned " << j << endl;
							j += 1;

							continue;
						}

						u_int64_t minimizer = minimizers[j];
						//if(bannedMinimizers.find(minimizer) != bannedMinimizers.end()) continue;

						currentMinimizerIndex.push_back(j);
						vec._kmers.push_back(minimizer);

						if(vec._kmers.size() == _kminmerSize){
							if(vec.isPalindrome()){

								for(size_t m=0; m<_kminmerSize-1; m++){
									bannedPositions[currentMinimizerIndex[m]] = true;
									//cout << "Banned: " << currentMinimizerIndex[m] << endl;
								}

								//exit(1);

								//cout << j << endl;
								//cout << j-1 << endl;
								//cout << j-2 << endl;
								//bannedPositions[j] = true;
								//bannedPositions[j-1] = true;
								//bannedPositions[j-2] = true;
								hasPalindrome = true;
								break;
							}
							else{
								//kminmers.push_back(vec.normalize());
								//cout << i << " " << j << endl;
								u_int32_t length = minimizers_pos[j] - minimizers_pos[i]; //minimizers_pos[i+windowSize-1] - minimizers_pos[i];
								kminmersData[vec.normalize()] = length;
								//cout << vec.toString() << endl;
								break;
							}
						}

						j += 1;
					}

					if(hasPalindrome) break;
				}

				
				if(!hasPalindrome) break;
			}
			*/


			/*
			vector<KmerVec> kminmers; 
			getKminmers(minimizers, filteredMinimizers, kminmers);

			foreach(KmerVec& vec : kminmers){

			}
			

			int i_max = ((int)minimizers.size()) - (int)_kminmerSize + 1;
			for(int i=0; i<i_max; i++){

				KmerVec vec;

				for(int j=i; j<i+_kminmerSize; j++){
					u_int64_t minimizer = minimizers[j];

					vec._kmers.push_back(minimizer);
				}

				u_int32_t length = minimizers_pos[i+_kminmerSize-1] - minimizers_pos[i]; //minimizers_pos[i+windowSize-1] - minimizers_pos[i];
				//kminmers.push_back(vec.normalize());

				//cout << "\t" << length << endl;
				kminmersData[vec.normalize()] = length;
				//mdbg_repeatFree->addNode(vec.normalize());
				//kminmers.push_back(vec.normalize());
			}
			*/

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



			u_int16_t size = minimizers.size();
			gzwrite(file_readData, (const char*)&size, sizeof(size));
			gzwrite(file_readData, (const char*)&minimizers[0], size * sizeof(u_int64_t));
			
			//comps.push_back({sequenceLength, composition});
			/*
			if(readIndex == 2229){
				cout << sequence_str.size() << endl;
				lala1 = {sequenceLength, composition};
				cout << "-------------------------- " << size << endl;
				for(u_int64_t m : minimizers){
					cout << m << endl;
				}
			}
			else if(readIndex == 3426){
				cout << sequence_str.size() << endl;
				lala2 = {sequenceLength, composition};
				cout << "-------------------------- " << size << endl;
				for(u_int64_t m : minimizers){
					cout << m << endl;
				}
				cout << "-------------------------- " << size << endl;
				cout << computeDistanceTNF(lala1, lala2) << endl;;
				cout << euclidianDistance(lala1._composition, lala2._composition) << endl;;
				cout << lala1._length << " " << lala2._length << endl;;
			}
			*/
			/*
			//cout << "-------------------------- " << size << endl;

			//gzwrite(file_readData, "\n", 1);
			//readData << endl;
			//string readData_str = readData.str();
			//gzwrite(file_readData, readData_str.c_str(), readData_str.size());

			//cout << minimizers.size() << endl;
			//if(minimizers.size() >= windowSize)
			int i_max = ((int)minimizers.size()) - (int)_kminmerSize + 1;
			//cout << minimizers.size() << " " << i_max << endl;
			for(int i=0; i<i_max; i++){

				//cout << i << " " << i_max << endl;
				KmerVec vec;

				for(int j=i; j<i+_kminmerSize; j++){
					//cout << minimizers.size() << " " << (j) << endl;
					vec._kmers.push_back(minimizers[j]);
				}
				
				//cout << vec._kmers.size() << endl;
				u_int32_t length = minimizers_pos[i+_kminmerSize-1] - minimizers_pos[i]; //minimizers_pos[i+windowSize-1] - minimizers_pos[i];
				//cout << "L: " << length << endl;
				//mdbg->addNode(vec.normalize(), length);



			}*/

			/*
			for(MinimizerPair& minimizerPair : minimizerPairs){
				
				minimizerPairs_to_reads[minimizerPair].push_back(readIndex);

				minimizerCounts[minimizerPair] += 1;

				_minimizerPairMap->addNode(minimizerPair);
				bool added = _overlapGraph->addNode(_minimizerPairMap->pair_to_id(minimizerPair));

				if(added){
					//output_file_gfa << "S" << "\t" << _minimizerPairMap->pair_to_id(minimizerPair) << "\t" << "*" << "\t" << "LN:i:" << 0 << endl;
				}
			}

			if(minimizerPairs.size() >= 2){
				
				MinimizerPair lastPair = minimizerPairs[0];
				for(size_t i=1; i<minimizerPairs.size(); i++){
					bool added = _overlapGraph->addEdge_checkDuplicate(_minimizerPairMap->pair_to_id(lastPair), _minimizerPairMap->pair_to_id(minimizerPairs[i]), 0);
					
					if(added){
						//output_file_gfa << "L" << "\t" << _minimizerPairMap->pair_to_id(lastPair) << "\t" << "+" << "\t" << _minimizerPairMap->pair_to_id(minimizerPairs[i]) << "\t" << "+" << "\t" << "0M" << endl;
					}
					
					//_overlapGraph->addEdge_checkDuplicate(_minimizerPairMap->pair_to_id(minimizerPairs[i]), _minimizerPairMap->pair_to_id(lastPair), 0);

					lastPair = minimizerPairs[i];
				}
				
			}
			*/

			/*
			//---------------------------------------------------------------------------------
			//cout << minimizerPairs.size() << endl;
			unordered_map<u_int64_t, u_int64_t> minimizerPairCount;
			for(MinimizerPair& minimizerPair : minimizerPairs){
				if(minimizerPairs_to_reads.find(minimizerPair) == minimizerPairs_to_reads.end()) continue;
				
				
				for(u_int64_t r : minimizerPairs_to_reads[minimizerPair]) {

					if(minimizerPairCount.find(r) == minimizerPairCount.end()){
						minimizerPairCount[r] = 0;
					}
					minimizerPairCount[r] += 1;
				}
			}

			//cout << minimizerPairCount.size() << endl;

			//cout << "----------------" << endl;
			//for (auto& it: minimizerPairCount) {
			//	cout << it.first << " " << it.second << endl;
			//}

			std::vector<std::pair<u_int64_t, u_int64_t>> elems(minimizerPairCount.begin(), minimizerPairCount.end());
			std::sort(elems.begin(), elems.end(), Ralalalala);
			
			//if(elems.size() > 1000){
				//cout << elems.size() << endl;
			//}
			//cout << "-" << endl;
			//cout << elems.size() << endl;
			u_int64_t readIndex_merged = -1;
			for (size_t i=0; i<elems.size(); i++) {
				

				u_int64_t read = elems[i].first;
				float dist = computeDistanceTNF(_readData[readIndex], _readData[read]);

				
				//cout << elems[i].first << " " << elems[i].second << " " << dist << endl;

				if(dist < 0.0005){ //0.05 0.002 0.0005
					readIndex_merged = read; //! disabled
					break;
				}
			}

			if(readIndex_merged == -1){
				for(MinimizerPair& minimizerPair : minimizerPairs){
					minimizerPairs_to_reads[minimizerPair].push_back(readIndex);
				}
			}
			else{
				for(MinimizerPair& minimizerPair : minimizerPairs){

					if(minimizerPairs_to_reads.find(minimizerPair) == minimizerPairs_to_reads.end()){
						minimizerPairs_to_reads[minimizerPair].push_back(readIndex_merged);
					}
					else{
						vector<u_int32_t>& minimizerPairReads = minimizerPairs_to_reads[minimizerPair];
						if(std::find(minimizerPairReads.begin(), minimizerPairReads.end(), readIndex_merged) == minimizerPairReads.end()){
							minimizerPairReads.push_back(readIndex_merged);
						}
					}


				}
				_readData[readIndex]._composition.clear();
			}
			*/
			//cout << minimizerPairs.size() << endl;
			//if(read_headers[readIndex] == "G3_0 S1_1"){
				//cout << minimizerPairs.size() << " " << _readData[readIndex]._length << endl;
			//}

			/*
			//S1_2132 S1_4386
			string str2 = sequence.getComment();
			//cout << str2 << endl;
			if (str2.find("S1_2132") != std::string::npos || str2.find("S1_4386") != std::string::npos) {
				cout << str2 << endl;
				cout << sequence.getDataSize() << endl;

				for(MinimizerPair& minimizerPair : minimizerPairs){
					cout << minimizerPair._first << " " << minimizerPair._second << endl;
				}
			}
			*/

			//---------------------------------------------------------------------------------


			//cout << "Nb minimizers: " << nbMinizersRead << endl;
			//_evaluation_readToDataset.push_back(datasetID);
			readIndex += 1;
		}


		datasetID += 1;
	}
	/*
	for(size_t i=0; i<comps.size(); i++){
		for(size_t j=i+1; j<comps.size(); j++){
			cout << computeDistanceTNF(comps[i], comps[j]) << endl;
		}
	}
	*/
	gzclose(file_readData);
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

	cout << "Nb kminmers: " << kminmersData.size() << endl;

	file_readData = gzopen(_filename_readMinimizers.c_str(),"rb");
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
		gzread(file_readData, (char*)&size, sizeof(size));

		if(gzeof(file_readData)) break;
		
		minimizers.resize(size);
		gzread(file_readData, (char*)&minimizers[0], size * sizeof(u_int64_t));


		vector<KmerVec> kminmers; 
		vector<ReadKminmer> kminmersInfo;
		vector<u_int64_t> minimizersPos; 
		vector<u_int64_t> rlePositions;
		MDBG::getKminmers(_minimizerSize, _kminmerSize, minimizers, minimizersPos, kminmers, kminmersInfo, rlePositions, 0);
		//getKminmers_filterRepeatedEdge(minimizers, filteredMinimizers, kminmers, kminmerCounts);


		vector<u_int32_t> abundances;

		for(KmerVec& vec : kminmers){
			abundances.push_back(kminmersData[vec]._count);
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
			//if(kminmerCounts[vec] > 1000) cout << kminmerCounts[vec] << endl;
			if(kminmersData[vec]._count <= 2){
				if(kminmersData[vec]._count <= minAbundance_cutoff) continue;
			}


			//cout << kminmersData[vec] << endl;

			KminmerData& kminmerData = kminmersData[vec];

			

			_mdbg->addNode(vec, kminmerData._length, kminmerData._overlapLength_start, kminmerData._overlapLength_end, kminmerData._isReversed);
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
		}

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

	}

	cout << "Nb solid kminmers: " << _mdbg->_dbg_nodes.size() << endl;


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
				u_int16_t overlapLength = min(_mdbg->_dbg_nodes[v]._overlapLength_end, _mdbg->_dbg_nodes[vec]._overlapLength_end);
				output_file_gfa << "L" << "\t" << id << "\t" << "+" << "\t" << _mdbg->_dbg_nodes[v]._index << "\t" << "+" << "\t" << overlapLength << "M" << endl;
				//vec_add_edge("+", "+");
			}
			if (vec.suffix() == v_rev.prefix()) {
				nbEdges += 1;
				u_int16_t overlapLength = min(_mdbg->_dbg_nodes[v]._overlapLength_end, _mdbg->_dbg_nodes[vec]._overlapLength_end);
				output_file_gfa << "L" << "\t" << id << "\t" << "+" << "\t" << _mdbg->_dbg_nodes[v]._index << "\t" << "-" << "\t" << overlapLength << "M" << endl;
				//vec_add_edge("+", "-");
			}
			if (vec_rev.suffix() == v.prefix()) {
				nbEdges += 1;
				u_int16_t overlapLength = min(_mdbg->_dbg_nodes[v]._overlapLength_start, _mdbg->_dbg_nodes[vec]._overlapLength_start);
				output_file_gfa << "L" << "\t" << id << "\t" << "-" << "\t" << _mdbg->_dbg_nodes[v]._index << "\t" << "+" << "\t" << overlapLength << "M" << endl;
				//vec_add_edge("-", "+");
			}
			if (vec_rev.suffix() == v_rev.prefix()) {
				nbEdges += 1;
				u_int16_t overlapLength = min(_mdbg->_dbg_nodes[v]._overlapLength_start, _mdbg->_dbg_nodes[vec]._overlapLength_start);
				output_file_gfa << "L" << "\t" << id << "\t" << "-" << "\t" << _mdbg->_dbg_nodes[v]._index << "\t" << "-" << "\t" << overlapLength << "M" << endl;
				//vec_add_edge("-", "-");
			}

		}
		for(KmerVec& v : _mdbg->_dbg_edges[vec.suffix().normalize()]){
			if(v==vec) continue;
			KmerVec v_rev = v.reverse();

			if (vec.suffix() == v.prefix()) {
				nbEdges += 1;
				u_int16_t overlapLength = min(_mdbg->_dbg_nodes[v]._overlapLength_end, _mdbg->_dbg_nodes[vec]._overlapLength_end);
				output_file_gfa << "L" << "\t" << id << "\t" << "+" << "\t" << _mdbg->_dbg_nodes[v]._index << "\t" << "+" << "\t" << overlapLength << "M" << endl;
				//vec_add_edge("+", "+");
			}
			if (vec.suffix() == v_rev.prefix()) {
				nbEdges += 1;
				u_int16_t overlapLength = min(_mdbg->_dbg_nodes[v]._overlapLength_end, _mdbg->_dbg_nodes[vec]._overlapLength_end);
				output_file_gfa << "L" << "\t" << id << "\t" << "+" << "\t" << _mdbg->_dbg_nodes[v]._index << "\t" << "-" << "\t" << overlapLength << "M" << endl;
				//vec_add_edge("+", "-");
			}
			if (vec_rev.suffix() == v.prefix()) {
				nbEdges += 1;
				u_int16_t overlapLength = min(_mdbg->_dbg_nodes[v]._overlapLength_start, _mdbg->_dbg_nodes[vec]._overlapLength_start);
				output_file_gfa << "L" << "\t" << id << "\t" << "-" << "\t" << _mdbg->_dbg_nodes[v]._index << "\t" << "+" << "\t" << overlapLength << "M" << endl;
				//vec_add_edge("-", "+");
			}
			if (vec_rev.suffix() == v_rev.prefix()) {
				nbEdges += 1;
				u_int16_t overlapLength = min(_mdbg->_dbg_nodes[v]._overlapLength_start, _mdbg->_dbg_nodes[vec]._overlapLength_start);
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
