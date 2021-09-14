

#include <Bloocoo.hpp>


Bloocoo::Bloocoo () : Tool("bloocoo")
{

	getParser()->push_back (new OptionOneParam (STR_INPUT, "input file", true));
	getParser()->push_back (new OptionOneParam (STR_OUTPUT, "output dir", true));
	getParser()->push_back (new OptionOneParam (STR_MINIM_SIZE, "minimizer length", false, "16"));
	getParser()->push_back (new OptionOneParam (STR_KMINMER_SIZE, "k-min-mer length", false, "3"));
	getParser()->push_back (new OptionOneParam (STR_DENSITY, "density of minimizers", false, "0.005"));
	getParser()->push_back (new OptionOneParam (STR_INPUT_DIR, "input dir", false, ""));
	getParser()->push_back (new OptionOneParam (STR_INPUT_EXTRACT_KMINMERS, "", false, ""));

}





void Bloocoo::execute ()
{

	_compositionManager = new CompositionManager(4);

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
	cout << "Kminmer length: " << _kminmerSize << endl;
	cout << "Minimizer length: " << _minimizerSize << endl;
	cout << "Density: " << _minimizerDensity << endl;
	cout << endl;

	_inputDir = getInput()->get(STR_INPUT_DIR) ? getInput()->getStr(STR_INPUT_DIR) : "";
	_input_extractKminmers= getInput()->get(STR_INPUT_EXTRACT_KMINMERS) ? getInput()->getStr(STR_INPUT_EXTRACT_KMINMERS) : "";

	_filename_readMinimizers = _outputDir + "/read_data.gz";
	_filename_readCompositions = _outputDir + "/read_compositions.gz";
	_filename_filteredMinimizers = _outputDir + "/filteredMinimizers.gz";
	_filename_hifiasmGroundtruth = _outputDir + "/hifiasmGroundtruth.gz";

	if(_input_extractKminmers != ""){
		extract_kminmers();
		//return;
	}
	if(_inputDir != ""){
		//execute_binning_cleanGraph();
		execute_binning();
		return;
	}

	unordered_map<MinimizerPair, vector<u_int32_t>> minimizerPairs_to_reads;

	vector<string> read_headers;

	//float minimizerDensity = 0.005; // 0.004; 0.0008
	u_int64_t maxHashValue = -1;
	//cout << maxHashValue << endl;
	//u_int64_t nbMinimizers = 0;

	u_int64_t _hash_otpt[2];
	int _seed = 42;
	setDispatcher (new SerialDispatcher());




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
	gzFile file_readComposition = gzopen(_filename_readCompositions.c_str(),"wb");
	
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



	_minimizerPairMap = new MinimizerPairMap();
	_overlapGraph = new AdjGraph();
	
	//unordered_map<u_int64_t, u_int64_t> minimCount;
	//set<MinimizerPair> _minimizerPairs;
	//vector<MinimizerPair_Edge> _minimizerPairs_edges;

	//unordered_map<MinimizerPair, u_int16_t> minimizerCounts;
	unordered_map<u_int64_t, u_int64_t> minimizerCounts;
	unordered_map<KmerVec, u_int16_t> kminmersData;



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



			char* readseq = sequence.getDataBuffer();
			string sequence_str;

			char lastChar = '0';
			for(size_t i=0; i<sequence.getDataSize(); i++){
				if(readseq[i] == lastChar) continue;
				sequence_str += readseq[i];
				lastChar = readseq[i];
			}


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
			Data buf((char*)sequence_str.c_str());
			u_int32_t sequenceLength = sequence_str.size();

			vector<float> composition;
			_compositionManager->readToComposition(buf, sequenceLength, composition);

			 

			gzwrite(file_readComposition, (const char*)&sequenceLength, sizeof(sequenceLength));
			gzwrite(file_readComposition, (const char*)&composition[0], _compositionManager->_compositionVectorSize * sizeof(float));
			

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
								//cout << sequence.toString() << endl;
								//exit(1);
								//cout << "PALOUF" << endl;
								//vec._kmers.pop_back();
								//for(size_t p=j-_kminmerSize+1; p<=j-1; p++){
								//	cout << "Banned: " << p << endl;
								//	bannedPositions[p] = true;
									//cout << "banned " << p << endl;
								//}

								/*
								for(size_t p=j-_kminmerSize+1; p<=j-1; p++){
									cout << "Banned: " << p << endl;
									bannedPositions[p] = true;
									//cout << "banned " << p << endl;
								}
								*/

								//cout << endl;
								//int loulalsdf = 0;
								//for(u_int64_t m : minimizers){
								//	cout << loulalsdf << ": " << m << endl;
								//	loulalsdf += 1;
								//}
								//cout << endl;

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
			_evaluation_readToDataset.push_back(datasetID);
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
	gzclose(file_readComposition);


	file_readData = gzopen(_filename_readMinimizers.c_str(),"rb");




	unordered_set<u_int64_t> filteredMinimizers;




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

		float median_abundance = compute_median(abundances_minimizers);
		//size_t minAbundance_cutoff = median_abundance / 10;
		size_t cutoff = median_abundance * 8;

		for(u_int64_t minimizer : minimizers){
			if(minimizerCounts[minimizer] > cutoff){
				filteredMinimizers.insert(minimizer);
			}
		}

	}

	gzclose(file_readData);
	cout << "Nb repeated minimizers: " << filteredMinimizers.size() << endl;

	gzFile file_filteredMinimiers = gzopen(_filename_filteredMinimizers.c_str(),"wb");
	for(u_int64_t minimizer : filteredMinimizers){
		gzwrite(file_filteredMinimiers, (char*)&minimizer, sizeof(minimizer));
		//cout << minimizer << endl;
	}
	gzclose(file_filteredMinimiers);



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
		getKminmers(minimizers, filteredMinimizers, kminmers);

		//cout << kminmers.size() << endl;
		for(KmerVec& vec : kminmers){
			kminmerCounts[vec] += 1;

			//kminmerPreSuf_Counts[vec.prefix()] += 1;
			//kminmerPreSuf_Counts[vec.suffix()] += 1;

			//cout << endl;
			//cout << vec.toString() << " " << kminmerCounts[vec] << endl;
		}

	}


	gzclose(file_readData);
	cout << "Nb kminmers: " << kminmerCounts.size() << endl;


	file_readData = gzopen(_filename_readMinimizers.c_str(),"rb");
	//cout << "Nb edges: " << nbEdges << endl;

	//cout << "Cleaning errors..." << endl;

	//delete mdbg;
	

	MDBG* mdbg = new MDBG(_kminmerSize);

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
		getKminmers(minimizers, filteredMinimizers, kminmers);
		//getKminmers_filterRepeatedEdge(minimizers, filteredMinimizers, kminmers, kminmerCounts);


		vector<u_int16_t> abundances;

		for(KmerVec& vec : kminmers){
			abundances.push_back(kminmerCounts[vec]);
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
			if(kminmerCounts[vec] <= 2){
				if(kminmerCounts[vec] <= minAbundance_cutoff) continue;
			}

			//cout << kminmersData[vec] << endl;
			mdbg->addNode(vec, kminmersData[vec]); //TODO length
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







	cout << "Writing gfa..." << endl;
	//cout << "Cleaning repeats..." << endl;

	//delete mdbg_repeatFree;
	//cout << _dbg_edges.size() << endl;

	u_int64_t nbEdges = 0;

	string gfa_filename = _outputDir + "/minimizer_graph.gfa";
	ofstream output_file_gfa(gfa_filename);

	for(auto vec_id : mdbg->_dbg_nodes){

		KmerVec vec = vec_id.first;
		KmerVec vec_rev = vec_id.first.reverse();
		u_int32_t id = vec_id.second._index;

		output_file_gfa << "S" << "\t" << id << "\t" << "*" << "\t" << "LN:i:" << vec_id.second._length << "\t" << "dp:i:" << vec_id.second._abundance << endl;

		//cout << mdbg->_dbg_edges[vec.prefix().normalize()].size() << endl;
		for(KmerVec& v : mdbg->_dbg_edges[vec.prefix().normalize()]){
			if(v==vec) continue;
			KmerVec v_rev = v.reverse();

			if (vec.suffix() == v.prefix()) {
				nbEdges += 1;
				output_file_gfa << "L" << "\t" << id << "\t" << "+" << "\t" << mdbg->_dbg_nodes[v]._index << "\t" << "+" << "\t" << "0M" << endl;
				//vec_add_edge("+", "+");
			}
			if (vec.suffix() == v_rev.prefix()) {
				nbEdges += 1;
				output_file_gfa << "L" << "\t" << id << "\t" << "+" << "\t" << mdbg->_dbg_nodes[v]._index << "\t" << "-" << "\t" << "0M" << endl;
				//vec_add_edge("+", "-");
			}
			if (vec_rev.suffix() == v.prefix()) {
				nbEdges += 1;
				output_file_gfa << "L" << "\t" << id << "\t" << "-" << "\t" << mdbg->_dbg_nodes[v]._index << "\t" << "+" << "\t" << "0M" << endl;
				//vec_add_edge("-", "+");
			}
			if (vec_rev.suffix() == v_rev.prefix()) {
				nbEdges += 1;
				output_file_gfa << "L" << "\t" << id << "\t" << "-" << "\t" << mdbg->_dbg_nodes[v]._index << "\t" << "-" << "\t" << "0M" << endl;
				//vec_add_edge("-", "-");
			}

		}
		for(KmerVec& v : mdbg->_dbg_edges[vec.suffix().normalize()]){
			if(v==vec) continue;
			KmerVec v_rev = v.reverse();

			if (vec.suffix() == v.prefix()) {
				nbEdges += 1;
				output_file_gfa << "L" << "\t" << id << "\t" << "+" << "\t" << mdbg->_dbg_nodes[v]._index << "\t" << "+" << "\t" << "0M" << endl;
				//vec_add_edge("+", "+");
			}
			if (vec.suffix() == v_rev.prefix()) {
				nbEdges += 1;
				output_file_gfa << "L" << "\t" << id << "\t" << "+" << "\t" << mdbg->_dbg_nodes[v]._index << "\t" << "-" << "\t" << "0M" << endl;
				//vec_add_edge("+", "-");
			}
			if (vec_rev.suffix() == v.prefix()) {
				nbEdges += 1;
				output_file_gfa << "L" << "\t" << id << "\t" << "-" << "\t" << mdbg->_dbg_nodes[v]._index << "\t" << "+" << "\t" << "0M" << endl;
				//vec_add_edge("-", "+");
			}
			if (vec_rev.suffix() == v_rev.prefix()) {
				nbEdges += 1;
				output_file_gfa << "L" << "\t" << id << "\t" << "-" << "\t" << mdbg->_dbg_nodes[v]._index << "\t" << "-" << "\t" << "0M" << endl;
				//vec_add_edge("-", "-");
			}
		}
		//for(KmerVec& v : _dbg_edges_prefix[vec.suffix().normalize()]){
		//	output_file_gfa << "L" << "\t" << vec_id.second << "\t" << "+" << "\t" << _dbg_nodes[v] << "\t" << "+" << "\t" << "0M" << endl;
		//}
	}

	output_file_gfa.close();
	
	cout << "Nb nodes: " << mdbg->_dbg_nodes.size() << endl;
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





	mdbg->dump(_outputDir + "/mdbg_nodes.gz");
	mdbg->_dbg_nodes.clear();
	mdbg->_dbg_edges.clear();



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

void Bloocoo::execute_binning(){

	//size_t globalAbundanceCutoff_min = 3;

	string gfa_filename = _inputDir + "/minimizer_graph.gfa";
	string gfa_filename_noUnsupportedEdges = _inputDir + "/minimizer_graph_noUnsupportedEdges.gfa";
	string gfa_filename_unitigs = _inputDir + "/minimizer_graph_unitigs.gfa";
	string mdbg_filename = _inputDir + "/mdbg_nodes.gz";

	cout << gfa_filename << endl;
	MDBG* mdbg = new MDBG(_kminmerSize);
	mdbg->load(mdbg_filename);

	cout << mdbg->_dbg_nodes.size() << endl;

	//vector<u_int32_t> unitigLengths;
	//vector<int32_t> node_to_unitig(mdbg->_dbg_nodes.size(), -1);
	GfaParser gfaParser;
	//BiGraph* graph = gfaParser.createBiGraph_lol(gfa_filename);
	AdjGraph* graph = gfaParser.createGraph_lol(gfa_filename);
	
	cout << "Nb nodes: " << graph->_nbNodes << endl;
	cout << "Nb edges: " << graph->_nbEdges << endl;








	unordered_set<u_int64_t> filteredMinimizers;



	//load_read_compositions();


	gzFile file_readData = gzopen(_filename_readMinimizers.c_str(),"rb");

	ReadIndexType readIndex = 0;
	//vector<UnitigData> unitigDatas;
	_unitigDatas.resize(graph->_nbNodes);

	//for(auto it : mdbg->_dbg_nodes){

		//const KmerVec& vec = it.first;

		//u_int32_t kminmer_index = mdbg->_dbg_nodes[vec]._index;
		//u_int32_t unitigIndex = node_to_unitig[kminmer_index];
		//if(unitigIndex == -1) continue;

		//unitigDatas[unitigIndex]._nbKminmers += 1;
		//unitigDatas[unitigIndex]._meanAbundance += mdbg->_dbg_nodes[vec]._abundance;
		//_unitigDatas[kminmer_index]._meanAbundance = mdbg->_dbg_nodes[vec]._abundance;
	//}
	/*
	for(size_t i=0; i<graph->_nbNodes; i++){
		//cout << i << endl;
		//unitigDatas.push_back({i, {}});
		//unitigDatas.push_back({i, {}, {}, 0});
		unitigDatas[i]._compositionMean.resize(_compositionManager->_compositionVectorSize);
	}*/

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

		/*
		if(readIndex == 75286){
			cout << "-------------------------- " << size << endl;
			for(u_int64_t m : minimizers){
				cout << m << endl;
			}
		}
		else if(readIndex == 75888){
			cout << "-------------------------- " << size << endl;
			for(u_int64_t m : minimizers){
				cout << m << endl;
			}
		}*/

		vector<ReadIndexType> unitigIndexex;
	
		for(KmerVec& vec : kminmers){
			//if(mdbg->_dbg_nodes[vec]._index == 55479) cout << "AAAAA" << endl;

			//cout << mdbg->_dbg_nodes[vec]._index << endl;
			if(mdbg->_dbg_nodes.find(vec) == mdbg->_dbg_nodes.end()) continue;

			u_int32_t kminmer_index = mdbg->_dbg_nodes[vec]._index;
			u_int32_t unitigIndex = kminmer_index;//node_to_unitig[kminmer_index];
			if(unitigIndex == -1) continue;

			/*
			UnitigData& unitigData = unitigDatas[unitigIndex];
			ReadData& readData = _readDatas[readIndex];
			for(size_t i=0; i<readData._composition.size(); i++){
				unitigData._compositionMean[i] += readData._composition[i];
			}
			unitigData._compositionNb += 1;
			*/

			//if(kminmer_index == 55479) cout << "AAAAA" << endl;
			//cout << kminmer_index << " " << unitigIndex << endl;
			if(std::find(unitigIndexex.begin(), unitigIndexex.end(), unitigIndex) != unitigIndexex.end()) continue;

			unitigIndexex.push_back(unitigIndex);
			if(unitigIndexex.size() >= 2) break;
		}

		if(unitigIndexex.size() >= 2){
			for(KmerVec& vec : kminmers){


				if(mdbg->_dbg_nodes.find(vec) == mdbg->_dbg_nodes.end()) continue;

				u_int32_t kminmer_index = mdbg->_dbg_nodes[vec]._index;
				u_int32_t unitigIndex = kminmer_index; //node_to_unitig[kminmer_index];
				//if(unitigIndex == -1) continue;

				UnitigData& unitigData = _unitigDatas[unitigIndex];
				
				//if(unitigData._readIndexes_exists.find(readIndex) != unitigData._readIndexes_exists.end()) continue;
				//cout << unitigData._readIndexes.size() << endl;
				//if(std::find(unitigData._readIndexes.begin(), unitigData._readIndexes.end(), readIndex) != unitigData._readIndexes.end()) continue;
					
				//unitigData._readIndexes_exists.insert(readIndex);
				unitigData._readIndexes.push_back(readIndex);
				/*
				cout << unitigDatas[0]._readIndexes.size() << endl;
				if(kminmer_index == 0){
					cout << (mdbg->_dbg_nodes.find(vec) == mdbg->_dbg_nodes.end()) << endl;
					cout << mdbg->_dbg_nodes[vec]._abundance << endl;
					cout << unitigDatas[0]._readIndexes.size() << endl;
					cout << kminmer_index << " " << unitigData._readIndexes.size() << " " << readIndex << endl;
				}
					cout << kminmer_index << " " << unitigData._readIndexes.size() << " " << readIndex << endl;*/
				//cout << unitigData._readIndexes.size() << " " << unitigLengths[unitigIndex] << endl;
				//if(std::find(unitigIndexex.begin(), unitigIndexex.end(), unitigIndex) != unitigIndexex.end()) continue
				//unitigIndexex.push_back(unitigIndex);

			}
		}


		
		readIndex += 1;
	}





	unordered_set<DbgEdge, hash_pair> unsupportedEdges;

	for(size_t utg=0; utg<graph->_nbNodes; utg++){

		//cout << utg << " " << graph->_nbNodes << endl;

		adjNode* node = graph->_nodes[utg];

        while (node != nullptr) {
			
			ReadIndexType utg_n = node->val;

			if(unsupportedEdges.find({utg, utg_n}) != unsupportedEdges.end() || unsupportedEdges.find({utg_n, utg}) != unsupportedEdges.end()) {	
				node = node->next;
				continue;
			}

			if(shareAnyRead(_unitigDatas[utg], _unitigDatas[utg_n])){
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
	gfaParser.rewriteGfa_withoutEdges(gfa_filename, gfa_filename_noUnsupportedEdges, unsupportedEdges);
	unsupportedEdges.clear();

	//delete graph; //NEED THIS FOR DEBUG TO EXTRACT HIFIASM SUB GRAPH
	//BiGraph* graphBi_successors = gfaParser.createBiGraph_lol(gfa_filename_noUnsupportedEdges, true);
	//BiGraph* graphBi_predecessors = gfaParser.createBiGraph_lol(gfa_filename_noUnsupportedEdges, false);
	//cout << "Nb nodes: " << (graphBi_successors->_nbNodes/2) << endl;
	//cout << "Nb edges: " << graphBi_successors->_nbEdges << endl;



	/*

	string command = "python3 ~/workspace/scripts/assembly/simplify_gfa.py " + gfa_filename_noUnsupportedEdges + " -out " + gfa_filename_unitigs;
	cout << command << endl;
	int ret = system(command.c_str());
	if(ret != 0){
		cout << "ERROR IN GFA TOOLS" << endl;
		exit(ret);
	}

	
	//vector<u_int32_t> unitigLengths;
	_node_to_unitig.resize(mdbg->_dbg_nodes.size(), -1);
	GfaParser::getNodeToUnitig(gfa_filename_unitigs, _node_to_unitig);
	
	*/

	//vector<int32_t> node_to_unitig(mdbg->_dbg_nodes.size(), -1);
	//string gfa_filename_unitigs = _inputDir + "/minimizer_graph_unitigs.gfa";
	//UnitigGraph* graph_unitig = gfaParser.createGraph(gfa_filename_unitigs, _node_to_unitig, unitigLengths);
	//delete graph_unitig;



	/*
	cout << "---------------------" << endl;
	for(u_int32_t r1 : unitigDatas[857640]._readIndexes){
		for(u_int32_t r2 : unitigDatas[857640]._readIndexes){
			//float tnf_dist = computeDistanceTNF(_readDatas[r1], _readDatas[r2]);
			//cout << tnf_dist << " " << readIndex_1 << " " << readIndex_2 << " " << utg << " " << utg_n << endl;
			cout << euclidianDistance(_readDatas[r1]._composition, _readDatas[r2]._composition) << " " << r1 << " " << r2 << " " << endl;
			//cout << euclidianDistance(_readDatas[readIndex_1]._composition, _readDatas[readIndex_2]._composition) << " " << utg << " " << utg_n << endl;
		}
	}
	cout << "---------------------" << endl;
	*/


	/*
	for(UnitigData& data : unitigDatas){
		for(size_t i=0; i<data._compositionMean.size(); i++){
			data._compositionMean[i] /= data._compositionNb;
		}

	}
	*/


	
	//Simulation
	ofstream file_groundTruth_hifiasm_position(_outputDir + "/groundtruth_hifiasm_position.csv");
	ofstream file_groundTruth_hifiasm(_outputDir + "/groundtruth_hifiasm.csv");
	file_groundTruth_hifiasm << "Name,Colour" << endl;
	file_groundTruth_hifiasm_position << "Name,Order" << endl;

	//unordered_set<u_int32_t> groundTruth_kminmers;
	int founded = 0;
	for(auto it : mdbg->_dbg_nodes){
		const KmerVec& vec = it.first;

		//vec = vec.normalize();
		if(_evaluation_hifiasmGroundTruth.find(vec) == _evaluation_hifiasmGroundTruth.end()) continue;

		founded += 1;
		//groundTruth_kminmers.insert(it.second._index);

		file_groundTruth_hifiasm << it.second._index << "," << _evaluation_hifiasmGroundTruth[vec] << endl;
		file_groundTruth_hifiasm_position << it.second._index << "," << _evaluation_hifiasmGroundTruth_position[vec] << endl;

		//vector<u_int32_t> neighbors;
		//graph->collectNeighbors(it.second._index, 100, neighbors, 100, visitedNodes);
		//for(u_int32_t nn : neighbors){
			//cout << nn << endl;
			//groundTruth_kminmers.insert(nn);
		//}
		//cout << "n: " << neighbors.size() << endl;
		//cout << n << endl;


		//cout << groundTruth_kminmers.size() << endl;
	}
	//cout << "Nb nodes abundant: " << groundTruth_kminmers.size() << endl;
	cout << "Founded: " << founded << endl;
	//gfaParser.rewriteGfa_withNodes(gfa_filename, gfa_filename + "_groundTruth_hifiasm.gfa", groundTruth_kminmers);
	file_groundTruth_hifiasm.close();
	file_groundTruth_hifiasm_position.close();
	

	
	/*
	cout << "Collecting truth kminmers" << endl;
	//ofstream file_groundTruth_hifiasm_position(_outputDir + "/groundtruth_hifiasm_position.csv");
	ofstream file_groundTruth_hifiasm(_outputDir + "/groundtruth_hifiasm.csv");
	file_groundTruth_hifiasm << "Name,Colour" << endl;
	//file_groundTruth_hifiasm_position << "Name,Colour" << endl;

	unordered_set<u_int32_t> visitedNodes;
	for(auto it : mdbg->_dbg_nodes){
		if(_evaluation_hifiasmGroundTruth.find(it.first) == _evaluation_hifiasmGroundTruth.end()) continue;
		visitedNodes.insert(it.second._index);
	}

	unordered_set<u_int32_t> groundTruth_kminmers;
	int founded = 0;
	for(auto it : mdbg->_dbg_nodes){

		const KmerVec& vec = it.first;

		//vec = vec.normalize();
		if(_evaluation_hifiasmGroundTruth.find(vec) == _evaluation_hifiasmGroundTruth.end()) continue;

		founded += 1;
		groundTruth_kminmers.insert(it.second._index);

		file_groundTruth_hifiasm << it.second._index << "," << _evaluation_hifiasmGroundTruth[vec] << endl;
		//file_groundTruth_hifiasm_position << it.second._index << "," << _evaluation_hifiasmGroundTruth_position[vec] << endl;

		vector<u_int32_t> neighbors;
		graph->collectNeighbors(it.second._index, 100, neighbors, 100, visitedNodes);
		for(u_int32_t nn : neighbors){
			//cout << nn << endl;
			groundTruth_kminmers.insert(nn);
		}
		//cout << "n: " << neighbors.size() << endl;
		//cout << n << endl;


		cout << groundTruth_kminmers.size() << endl;
	}
	cout << "Nb nodes abundant: " << groundTruth_kminmers.size() << endl;
	cout << "Founded: " << founded << endl;
	gfaParser.rewriteGfa_withNodes(gfa_filename_noUnsupportedEdges, gfa_filename + "_groundTruth_hifiasm.gfa", groundTruth_kminmers);
	file_groundTruth_hifiasm.close();
	//file_groundTruth_hifiasm_position.close();
	*/
	


	/*
	computeNodeAbundance(mdbg, gfa_filename_noUnsupportedEdges);
	getUnitigLengths(gfa_filename_unitigs);
	*/

	cout << "Simplifying graph" << endl;
	GraphSimplify* graphSimplify = new GraphSimplify(gfa_filename_noUnsupportedEdges, _outputDir);
	graphSimplify->execute(1);
	//graphSimplify->compact();
	
	
	exit(1);
	
	file_groundTruth = ofstream(_outputDir + "/binning_results.csv");
	file_groundTruth << "Name,Colour" << endl;

	//562 (ecoli)
	solveBin(graphSimplify->_graphSuccessors->nodeName_to_nodeIndex(9847, true), graphSimplify->nodeIndex_to_unitig(9847)._abundance, graphSimplify, 0);
	//solveBin(graphSimplify->_graphSuccessors->nodeName_to_nodeIndex(10645, false), graphSimplify->nodeIndex_to_unitig(10645)._abundance, graphSimplify, 0);


	exit(1);



	//file_groundTruth.close();

	unordered_set<u_int32_t> binNodes;
	
	
	unordered_set<u_int32_t> visitedNodes;
	vector<u_int32_t> startingNodesIndex;

	vector<UnitigLength> unitigLengths;
	for(Unitig& unitig : graphSimplify->_unitigs){
		unitigLengths.push_back({unitig._index, unitig._length, unitig._startNode});
	}
	std::sort(unitigLengths.begin(), unitigLengths.end(), UnitigComparator_ByLength);

	bool orient_dummy;
	//size_t binIndex=0;
	for(UnitigLength& unitigLength : unitigLengths){
		cout << "Unitig length: " << unitigLength._length << " " << unitigLength._index << endl;
		
		if(unitigLength._length < 10000) continue;
		//Unitig& unitig = graphSimplify->_unitigs[unitigLength._index];
		//if(graphSimplify->_nodeToUnitig[unitig._startNode] == -1) continue;

		vector<u_int32_t> nodes;
		Unitig& unitig = graphSimplify->_unitigs[unitigLength._index];

		unitigLength._abundance = unitig._abundance;

		float abundanceCutoff_min = computeAbundanceCutoff_min(unitigLength._abundance);
		//if(abundanceCutoff_min < 30) continue;


		
		graphSimplify->getUnitigNodes(unitig, nodes);

		u_int32_t nodeIndex = nodes[rand() % nodes.size()];
		//u_int32_t nodeName = graphSimplify->_graphSuccessors->nodeIndex_to_nodeName(nodeIndex, orient_dummy);
		//if(_binnedNodes.find(nodeName) != _binnedNodes.end()) continue;



		unitigLength._startNodeIndex = nodeIndex;
		//startingNodesIndex.push_back(nodeIndex);
	}

	//for(size_t n=0; n<graphSimplify->_nodeToUnitig.size(); n++){
	//	u_int32_t unitigIndex = graphSimplify->_nodeToUnitig[n];
	//	graphSimplify->_nodeAbundances[n] = graphSimplify->_unitigs[unitigIndex]._abundance;
	//}

	//genome3
	//solveBin(graphSimplify->_graphSuccessors->nodeName_to_nodeIndex(1869, true), graphSimplify, 0);
	
	//genome_201_50x
	//solveBin(graphSimplify->_graphSuccessors->nodeName_to_nodeIndex(115862, false), graphSimplify, 0);
	
	//exit(1);

	size_t binIndex=0;
	//for(u_int32_t nodeIndex : startingNodesIndex){
	for(UnitigLength& unitigLength : unitigLengths){
		if(unitigLength._index % 2 == 1) continue;
		cout << "Unitig length: " << unitigLength._length << " " << unitigLength._index << endl;

		if(unitigLength._length < 10000) continue;
		float abundanceCutoff_min = computeAbundanceCutoff_min(unitigLength._abundance);
		//if(abundanceCutoff_min < 30) continue;
		if(visitedNodes.find(unitigLength._index) != visitedNodes.end()) continue;

		u_int32_t nodeIndex = unitigLength._startNodeIndex;
		
		u_int32_t nodeName = graphSimplify->_graphSuccessors->nodeIndex_to_nodeName(nodeIndex, orient_dummy);
		if(_binnedNodes.find(nodeName) != _binnedNodes.end()) continue;

		visitedNodes.insert(unitigLength._index);
		//vector<u_int32_t> nodes;
		//Unitig& unitig = graphSimplify->_unitigs[unitigLength._index];
		//graphSimplify->getUnitigNodes(unitig, nodes);
		//for(u_int32_t node : nodes){
		//	visitedNodes.insert(node);
		//}
		/*
		cout << "Unitig length: " << unitigLength._length << " " << unitigLength._index << endl;
		
		if(unitigLength._length < 10000) continue;

		//Unitig& unitig = graphSimplify->_unitigs[unitigLength._index];
		//if(graphSimplify->_nodeToUnitig[unitig._startNode] == -1) continue;

		vector<u_int32_t> nodes;
		Unitig& unitig = graphSimplify->_unitigs[unitigLength._index];
		graphSimplify->getUnitigNodes(unitig, nodes);

		u_int32_t nodeIndex = nodes[rand() % nodes.size()];
		u_int32_t nodeName = graphSimplify->_graphSuccessors->nodeIndex_to_nodeName(nodeIndex, orient_dummy);
		if(_binnedNodes.find(nodeName) != _binnedNodes.end()) continue;

		float abundanceCutoff_min = computeAbundanceCutoff_min(nodeIndex, graphSimplify);
		if(abundanceCutoff_min < 30) continue;
		//u_int32_t nodeName = graphSimplify->_graphSuccessors->nodeIndex_to_nodeName(unitigLength.startNode, orient_dummy);
		//if(_binnedNodes.find(nodeName) != _binnedNodes.end()) continue;
			
		*/

		//for(size_t n=0; n<graphSimplify->_nodeToUnitig.size(); n++){
		//	if(n%2 == 1) continue;


			//if(graphSimplify->_nodeToUnitig[n] == unitigLength._index){

				
				//file_groundTruth = ofstream(_outputDir + "/binning_results.csv");
	
				solveBin(nodeIndex, unitigLength._abundance, graphSimplify, binIndex);
				//file_groundTruth.close();

				//exit(1);
				binIndex += 1;
				//break;
			//}
		//}
		
	}

	cout << "Nb bins: " << binIndex << endl;
	//solveBin(611, graphBi);

	//solveBin(16057, graphBi, 0); //1
	//solveBin(18318, graphBi, 1); //715

	//solveBin(11274, graphBi); //2154
	//solveBin(15152, graphBi); //2490
	//solveBin(10075, graphBi); //3717
	//solveBin(17845, graphBi); //4156

	


	//solveBin(4229, graph, unitigDatas, file_groundTruth, 1, binNodes);
	//solveBin(16049, graph);


	//solveBin(6751, graph, unitigDatas, file_groundTruth, 1, binNodes);


	//home/gats/workspace/run/hifiasm_meta/AD_components/component_5.fasta
	//solveBin(90422, graphBi);
	
	//home/gats/workspace/run/hifiasm_meta/AD_components/component_6.fasta
	//solveBin(2144353, graph, unitigDatas, file_groundTruth, 1);
	
	//home/gats/workspace/run/hifiasm_meta/AD_components/component_7.fasta
	//solveBin(137362, graph, unitigDatas, file_groundTruth, 1);

	//home/gats/workspace/run/hifiasm_meta/AD_components/component_8.fasta
	//solveBin(446417, graph, unitigDatas, file_groundTruth, 1);



	//solveBin(721346, graph, unitigDatas, file_groundTruth, 1);

	//solveBin(1057286, graph, unitigDatas, file_groundTruth, 1);

	/*
	solveBin(4106, graph, unitigDatas, file_groundTruth, 1);
	solveBin(51899, graph, unitigDatas, file_groundTruth, 2);
	solveBin(95193, graph, unitigDatas, file_groundTruth, 3);
	solveBin(10985, graph, unitigDatas, file_groundTruth, 4);

	solveBin(78693, graph, unitigDatas, file_groundTruth, 5);
	solveBin(92561, graph, unitigDatas, file_groundTruth, 6);
	solveBin(109366, graph, unitigDatas, file_groundTruth, 7);
	solveBin(234177, graph, unitigDatas, file_groundTruth, 8);
	solveBin(143694, graph, unitigDatas, file_groundTruth, 9);
	solveBin(25798, graph, unitigDatas, file_groundTruth, 10);
	*/

	file_groundTruth.close();

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
}


void Bloocoo::execute_binning_cleanGraph(){
	string gfa_filename = _inputDir + "/minimizer_graph_unitigs.gfa";
	string mdbg_filename = _inputDir + "/mdbg_nodes.gz";

	cout << gfa_filename << endl;
	MDBG* mdbg = new MDBG(_kminmerSize);
	mdbg->load(mdbg_filename);

	vector<u_int32_t> unitigLengths;
	vector<int32_t> node_to_unitig(mdbg->_dbg_nodes.size(), -1);
	GfaParser gfaParser;
	UnitigGraph* graph = gfaParser.createGraph(gfa_filename, node_to_unitig, unitigLengths);
	
	cout << "Nb nodes: " << (graph->_nbNodes/2) << endl;
	cout << "Nb edges: " << graph->_nbEdges << endl;




	BiGraph* graph_bi = gfaParser.createBiGraph_lol(_inputDir + "/minimizer_graph_notips_nobubbles.gfa", true);
	














	
	load_read_compositions();
	//load_filtered_minimizers();



	
	
	gzFile file_readData = gzopen(_filename_readMinimizers.c_str(),"rb");

	ReadIndexType readIndex = 0;
	_unitigDatas.resize(graph->_nbNodes);
	//vector<UnitigData> unitigDatas;
	//for(size_t i=0; i<graph->_nbNodes; i++){
	//	unitigDatas.push_back({i, {}});
		//unitigDatas[i]._compositionMean.resize(_compositionManager->_compositionVectorSize);
	//}



	/*
	cout << "Collecting truth kminmers" << endl;
	//ofstream file_groundTruth_hifiasm_position(_outputDir + "/groundtruth_hifiasm_position.csv");
	ofstream file_groundTruth_hifiasm(_outputDir + "/groundtruth_hifiasm.csv");
	file_groundTruth_hifiasm << "Name,Colour" << endl;
	//file_groundTruth_hifiasm_position << "Name,Colour" << endl;

	unordered_set<u_int32_t> groundTruth_kminmers;
	unordered_set<u_int32_t> groundTruth_unitig;
	int founded = 0;
	for(auto it : mdbg->_dbg_nodes){
		//for(auto it : _evaluation_hifiasmGroundTruth){

		const KmerVec& vec = it.first;

		u_int32_t kminmer_index = mdbg->_dbg_nodes[vec]._index;
		u_int32_t unitigIndex = node_to_unitig[kminmer_index];
		if(unitigIndex == -1) continue;

		unitigDatas[unitigIndex]._nbKminmers += 1;
		unitigDatas[unitigIndex]._meanAbundance += mdbg->_dbg_nodes[vec]._abundance;



		//vec = vec.normalize();
		if(_evaluation_hifiasmGroundTruth.find(vec) == _evaluation_hifiasmGroundTruth.end()) continue;
		//if(mdbg->_dbg_nodes.find(vec) == mdbg->_dbg_nodes.end()) continue;


		groundTruth_kminmers.insert(kminmer_index);
		founded += 1;


		groundTruth_unitig.insert(unitigIndex);
		file_groundTruth_hifiasm << gfaParser.unitigIndex_to_unitigName(unitigIndex) << "," << _evaluation_hifiasmGroundTruth[vec] << endl;
		//file_groundTruth_hifiasm_position << unitigIndex << "," << _evaluation_hifiasmGroundTruth_position[vec] << endl;

		vector<u_int32_t> neighbors;
		graph->collectNeighbors(unitigIndex, 100, neighbors, 100);
		for(u_int32_t nn : neighbors){
			//cout << nn << endl;
			groundTruth_kminmers.insert(nn);
			groundTruth_unitig.insert(nn);
			//unitigDatas[nn]._nbKminmers += 1;
			//unitigDatas[nn]._meanAbundance += mdbg->_dbg_nodes[vec]._abundance;
			//u_int32_t unitigIndex = node_to_unitig[nn];
			//if(unitigIndex == -1) continue;
			//groundTruth_unitig.insert(unitigIndex);
		}

		
		//cout << founded << endl;
		//cout << "n: " << neighbors.size() << endl;
		//cout << n << endl;


		cout << groundTruth_kminmers.size() << endl;
	}

	for(UnitigData& data : unitigDatas){
		if(data._nbKminmers != 0){
			data._meanAbundance /= data._nbKminmers;
		}
	}
	cout << "Nb nodes abundant: " << groundTruth_kminmers.size() << endl;
	cout << "Founded: " << founded << endl;
	gfaParser.rewriteGfa_withUnitigs(gfa_filename, gfa_filename + "_groundTruth_hifiasm.gfa", groundTruth_unitig, unitigDatas);
	cout << "endo" << endl;
	file_groundTruth_hifiasm.close();
	//file_groundTruth_hifiasm_position.close();
	*/


	while(true){
		
		//cout << readIndex << endl;

		u_int16_t size;
		vector<u_int64_t> minimizers;
		gzread(file_readData, (char*)&size, sizeof(size));

		if(gzeof(file_readData)) break;
		
		minimizers.resize(size);
		gzread(file_readData, (char*)&minimizers[0], size * sizeof(u_int64_t));


		vector<KmerVec> kminmers; 
		getKminmers(minimizers, _filteredMinimizers, kminmers);

		vector<ReadIndexType> unitigIndexex;

		for(KmerVec& vec : kminmers){
			if(mdbg->_dbg_nodes.find(vec) == mdbg->_dbg_nodes.end()) continue;

			u_int32_t kminmer_index = mdbg->_dbg_nodes[vec]._index;
			u_int32_t unitigIndex = node_to_unitig[kminmer_index];
			if(unitigIndex == -1) continue;

			/*
			UnitigData& unitigData = unitigDatas[unitigIndex];
			ReadData& readData = _readDatas[readIndex];
			for(size_t i=0; i<readData._composition.size(); i++){
				unitigData._compositionMean[i] += readData._composition[i];
			}
			unitigData._compositionNb += 1;
			*/

			//if(groundTruth_unitig.find(unitigIndex) == groundTruth_unitig.end()) continue;

			//cout << kminmer_index << " " << unitigIndex << endl;
			if(std::find(unitigIndexex.begin(), unitigIndexex.end(), unitigIndex) != unitigIndexex.end()) continue;

			unitigIndexex.push_back(unitigIndex);
			if(unitigIndexex.size() >= 2) break;
		}

		if(unitigIndexex.size() >= 2){
			for(KmerVec& vec : kminmers){
				if(mdbg->_dbg_nodes.find(vec) == mdbg->_dbg_nodes.end()) continue;

				u_int32_t kminmer_index = mdbg->_dbg_nodes[vec]._index;
				u_int32_t unitigIndex = node_to_unitig[kminmer_index];
				if(unitigIndex == -1) continue;

				//if(groundTruth_unitig.find(unitigIndex) == groundTruth_unitig.end()) continue;
				//if(groundTruth_unitig.find(unitigIndex) == groundTruth_unitig.end()) continue;

				UnitigData& unitigData = _unitigDatas[unitigIndex];
				
				//if(unitigData._readIndexes_exists.find(readIndex) != unitigData._readIndexes_exists.end()) continue;
				//cout << unitigData._readIndexes.size() << endl;
				//if(std::find(unitigData._readIndexes.begin(), unitigData._readIndexes.end(), readIndex) != unitigData._readIndexes.end()) continue;
					
				//unitigData._readIndexes_exists.insert(readIndex);
				//cout << unitigIndex << " " << unitigData._readIndexes.size() << endl;
				//cout << unitigData._readIndexes.size() << endl;
				unitigData._readIndexes.push_back(readIndex);

				//cout << kminmer_index << " " << unitigIndex << " " << unitigData._readIndexes.size() << " " << readIndex << endl;
				//cout << unitigData._readIndexes.size() << " " << unitigLengths[unitigIndex] << endl;
				//if(std::find(unitigIndexex.begin(), unitigIndexex.end(), unitigIndex) != unitigIndexex.end()) continue
				//unitigIndexex.push_back(unitigIndex);

			}
		}


		
		readIndex += 1;
	}

	/*
	for(UnitigData& data : unitigDatas){
		for(size_t i=0; i<data._compositionMean.size(); i++){
			data._compositionMean[i] /= data._compositionNb;
		}

	}
	*/




	//solveBin(GfaParser::unitigName_to_id("utg0000804l"), graph);
	

	//solveBin(GfaParser::unitigName_to_id("utg0000805l"), graph, unitigDatas, file_groundTruth, 1);

	///home/gats/workspace/run/hifiasm_meta/AD_components/component_5.fasta
	//solveBin(GfaParser::unitigName_to_id("utg0437330l"), graph, unitigDatas, file_groundTruth, 1);
	
	///home/gats/workspace/run/hifiasm_meta/AD_components/component_6.fasta
	//solveBin(GfaParser::unitigName_to_id("utg0459582l"), graph, unitigDatas, file_groundTruth, 1);


	///home/gats/workspace/run/hifiasm_meta/AD_components/component_5.fasta
	//solveBin(721346, graph, unitigDatas, file_groundTruth, 1);

	///home/gats/workspace/run/hifiasm_meta/AD_components/component_6.fasta
	//solveBin(1057286, graph, unitigDatas, file_groundTruth, 1);

	///home/gats/workspace/run/hifiasm_meta/AD_components/component_8.fasta
	//solveBin(GfaParser::unitigName_to_id("utg0198809l"), graph, unitigDatas, file_groundTruth, 1, binNodes);

	///home/gats/workspace/run/hifiasm_meta/AD_components/component_9.fasta
	//solveBin(GfaParser::unitigName_to_id("utg0435813l"), graph, unitigDatas, file_groundTruth, 1, binNodes);



	/*
	solveBin(4106, graph, unitigDatas, file_groundTruth, 1);
	solveBin(51899, graph, unitigDatas, file_groundTruth, 2);
	solveBin(95193, graph, unitigDatas, file_groundTruth, 3);
	solveBin(10985, graph, unitigDatas, file_groundTruth, 4);

	solveBin(78693, graph, unitigDatas, file_groundTruth, 5);
	solveBin(92561, graph, unitigDatas, file_groundTruth, 6);
	solveBin(109366, graph, unitigDatas, file_groundTruth, 7);
	solveBin(234177, graph, unitigDatas, file_groundTruth, 8);
	solveBin(143694, graph, unitigDatas, file_groundTruth, 9);
	solveBin(25798, graph, unitigDatas, file_groundTruth, 10);
	*/


	//solveBin(GfaParser::unitigName_to_id("utg0000799l"), graph);

	ofstream file_groundTruth(_outputDir + "/binning_results.csv");
	file_groundTruth << "Name,Colour" << endl;
	unordered_set<u_int32_t> binNodes;
	//TODO
	file_groundTruth.close();

	
	/*
	cout << "Writing gfa" << endl;

	unordered_set<u_int32_t> unitigs;
	for(auto it : mdbg->_dbg_nodes){

		const KmerVec& vec = it.first;

		u_int32_t kminmer_index = mdbg->_dbg_nodes[vec]._index;
		u_int32_t unitigIndex = node_to_unitig[kminmer_index];
		if(unitigIndex == -1) continue;


		if(binNodes.find(unitigIndex) == binNodes.end()) continue;
		unitigs.insert(unitigIndex);

		vector<u_int32_t> neighbors;
		graph->collectNeighbors(unitigIndex, 100, neighbors, 100);
		for(u_int32_t nn : neighbors){
			unitigs.insert(nn);
		}

	}
	gfaParser.rewriteGfa_withUnitigs(gfa_filename, gfa_filename + "_binResults.gfa", unitigs, unitigDatas);
	*/
	//file_groundTruth_hifiasm_position.close();



	/*
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
}


	/*
	vector<u_int32_t> node_to_unitig(mdbg_errorFree->_dbg_nodes.size(), 0);
	GfaParser gfaParser;
	AdjGraph* graph = gfaParser.createGraph("/home/gats/workspace/data/overlap_test/read_overlaps_notips_nobubbles.gfa", node_to_unitig);
	
	vector<vector<u_int32_t>> components;
	graph->computeConnectedComponents(components);
	
	cout << "Nb connected components: " << components.size() << endl;
	for(size_t i=0; i<components.size(); i++){
		if(components[i].size() > 50){
			cout << i <<": " << components[i].size() << endl;
		}
		//if(components[i].size() == 1) cout << components[i][0] << endl;
	}
	*/


	/*
	cout << _overlapGraph->_nbNodes << endl;
	cout << _overlapGraph->_nbEdges << endl;


	
	AdjGraph* minimizerGraph_cleaned = new AdjGraph();
	MinimizerPairMap* _minimizerPairMap_cleaned = new MinimizerPairMap();

	//for (auto& it: minimizerCounts) {
	//	if(it.second <= 1) continue;
	//	if(it.second <= 10) continue;
	//	cout << it.second << endl;
	//}

	for(size_t n=0; n<_overlapGraph->_nbNodes; n++){

		MinimizerPair& pair = _minimizerPairMap->id_to_pair(n);



		if(!checkRemoveNode(n, _overlapGraph, _minimizerPairMap, minimizerCounts)){
			//_overlapGraph->collectNeighbors(n, 1, neighbors);
			//if(minimizerCounts[pair] > 1 && minimizerCounts[pair] < 35){
			_minimizerPairMap_cleaned->addNode(pair);
			bool added = minimizerGraph_cleaned->addNode(_minimizerPairMap_cleaned->pair_to_id(pair));
			if(added){
				//output_file_gfa << "S" << "\t" << _minimizerPairMap->pair_to_id(pair) << "\t" << "*" << "\t" << "LN:i:" << 0 << "\t" << "dp:i:" << minimizerCounts[pair] << endl;
			}
		}


		adjNode* node = _overlapGraph->_nodes[n];

        while (node != nullptr) {
			
			u_int32_t nn = node->val;
			MinimizerPair& pair_nn = _minimizerPairMap->id_to_pair(nn);

			if(!checkRemoveNode(nn, _overlapGraph, _minimizerPairMap, minimizerCounts)){
				//if(minimizerCounts[pair_nn] > 1 && minimizerCounts[pair_nn] < 35){
				_minimizerPairMap_cleaned->addNode(pair_nn);
				bool added = minimizerGraph_cleaned->addNode(_minimizerPairMap_cleaned->pair_to_id(pair_nn));
				if(added){
					//output_file_gfa << "S" << "\t" << _minimizerPairMap->pair_to_id(pair_nn) << "\t" << "*" << "\t" << "LN:i:" << 0 << "\t" << "dp:i:" << minimizerCounts[pair_nn] << endl;
				}

				added = minimizerGraph_cleaned->addEdge_checkDuplicate(_minimizerPairMap_cleaned->pair_to_id(pair), _minimizerPairMap_cleaned->pair_to_id(pair_nn), 0);
				if(added){
					//output_file_gfa << "L" << "\t" << _minimizerPairMap->pair_to_id(pair) << "\t" << "+" << "\t" << _minimizerPairMap->pair_to_id(pair_nn) << "\t" << "+" << "\t" << "0M" << endl;
				}

			}

			node = node->next;
        }

	}

	cout << minimizerGraph_cleaned->_nbNodes << endl;
	cout << minimizerGraph_cleaned->_nbEdges << endl;

	output_file_gfa.close();
	*/

	/*
	vector<float> dists;
	for(size_t n=0; n<_readData.size(); n++){
		

		for(size_t nn=n+1; nn<_readData.size(); nn++){
			float dist = computeDistanceTNF(_readData[n], _readData[nn]);

			//cout << dist << endl;
			
			dists.push_back(dist);

		}
	}

	print_stats(dists);
	*/
	
	//cout << nbMinimizers << endl;

	/*
	u_int64_t processed = 0;

	vector<vector<u_int16_t>> readCoverages(readIndex);
	unordered_map<MinimizerPair, u_int16_t> readPairCount;

	cout << minimizerPairs_to_reads.size() << endl;

	for (auto& it: minimizerPairs_to_reads) {

		processed += 1;
		if(processed % 100000 == 0){
			float progress = ((float) processed) / minimizerPairs_to_reads.size();
			//cout << progress << endl;
		} 


		MinimizerPair minimizerPair = it.first;
		vector<u_int32_t>& readIds = it.second;



		//u_int32_t m1 = (u_int32_t)((minimizerPair & 0xFFFFFFFF00000000LL) >> 32);
		//u_int32_t m2 = (u_int32_t)(minimizerPair & 0xFFFFFFFFLL);
		//cout << "Read pair: " << m1 << " " << m2 << " " << minimizerPair << endl;

		//cout << readIds.size() << endl;
		//cout << "------------" << endl;
		//for(size_t k=0; k<readIds.size(); k++) {
		//	cout << readIds[k] << endl;
		//}



		if(readIds.size() == 1) continue;

		//cout << readIds.size() << endl;
		if(readIds.size() > 10000) continue;

		//cout << readIds.size() << endl;
		for(size_t i=0; i<readIds.size(); i++) {
			
				//cout << i << " " << readIds.size() << endl;

			u_int32_t read1 = readIds[i];

			//if(read1 == 9363) cout << "Reads: " << readIds.size() << endl; 

			for(size_t j=i+1; j<readIds.size(); j++) {

				u_int32_t read2 = readIds[j];
				//if(read2 == 9363) cout << "Reads: " << readIds.size() << endl; 


				//if ((read1 == 7721 && read2 == 9442) || read2 == 7721 && read1 == 9442) {
					//cout << "lala" << endl;
				//}

				//cout << read1 << " " << read2 << endl;

				if(read1 == read2){
					//cout << "------------" << endl;
					//for(size_t k=0; k<readIds.size(); k++) {
						//cout << readIds[k] << endl;
					//}
					continue;
				}

				float dist = computeDistanceTNF(_readData[read1], _readData[read2]);
				//float t = 0.01f + 
				if(dist > 0.6f){
					//cout << read1 << " " << read2 << " " << dist << endl;
					continue;
				}

				

				


				//if(_evaluation_readToDataset[read1] != _evaluation_readToDataset[read2]){
				//	cout << "Diff dataset: " << computeDistanceTNF(_readData[read1], _readData[read2]) << endl;
				//}
				//else{
				//	cout << "Same dataset: " << computeDistanceTNF(_readData[read1], _readData[read2]) << endl;
				//}

				MinimizerPair readPair = {read1, read2};
				//u_int64_t readPair = (u_int64_t) read1 << 32 | read2;

				if(readPairCount.find(readPair) == readPairCount.end()){
					//cout << "Overlap: " << read1 << " " << read2 << endl;
					readPairCount[readPair] = 0;
					//cout << readPairCount.size() << endl;
				}

				readPairCount[readPair] += 1;


				//if ((readPair._first == 7721 && read2 == 9442) || read2 == 7721 && read1 == 9442) {
					//cout << readPairCount[readPair] << endl;
				//}

				//cout << j << endl;

			}
		}

		
		//cout << "Reads: " << readIds.size() << endl;


	}


	*/
	/*
	for(size_t i=0; i<readCoverages.size(); i++){
		if(read_headers[i] == "G3_0 S1_1"){

			for(size_t cov : readCoverages[i]){
				cout << cov << endl;
			}

			float mean = 0;
			
			for(size_t cov : readCoverages[i]){
				mean += cov;
			}
			mean /= readCoverages[i].size();

			float var = 0;
			for(size_t cov : readCoverages[i]){
				var += pow((cov - mean), 2);
			}
			var /= (readCoverages[i].size() - 1);
			float sd = SQRT(var);
			
			cout << "Mean: " << mean << endl;
			cout << "Sd: " << sd << endl;
			cout << "Var: " << var << endl;

		}
	}
	*/

	/*
	set<u_int32_t> nodes;


	//vector<u_int32_t> nodes;

	//ofstream output_file("/home/gats/workspace/data/overlap_test/read_overlaps.txt");
	ofstream output_file_gfa("/home/gats/workspace/data/overlap_test/read_overlaps.gfa");

	u_int64_t nbOverlaps = 0;

	u_int64_t nbOverlaps_interGenomes = 0;

	
	//u_int64_t nbNodes = 0;
	for (auto& it: readPairCount) {
		u_int16_t count = it.second;
		//if(count <= 1) continue;


		MinimizerPair readPair = it.first;
		u_int32_t read1 = readPair._first;//(u_int32_t)((readPair & 0xFFFFFFFF00000000LL) >> 32);
		u_int32_t read2 = readPair._second;//(u_int32_t)(readPair & 0xFFFFFFFFLL);

		nodes.insert(read1);
		nodes.insert(read2);
		//nbNodes += 1;
	}
	
	GraphInfo* graphInfo = new GraphInfo();

	cout << "Nb nodes: " << nodes.size() << endl;
	//cout << readIndex << " " << nbNodes << endl;
	//u_int64_t nbNodes = readIndex;
	_overlapGraph = new AdjGraph(nodes.size());

	vector<Overlap> overlaps;

	int nbOverlapsLala = 0;
	for (auto& it: readPairCount) {
		
		u_int16_t count = it.second;
		if(count <= 1) continue;


		MinimizerPair readPair = it.first;

		u_int32_t read1 = readPair._first;//(u_int32_t)((readPair & 0xFFFFFFFF00000000LL) >> 32);
		u_int32_t read2 = readPair._second;//(u_int32_t)(readPair & 0xFFFFFFFFLL);

		overlaps.push_back({read1, read2, count});


		//if(read1 == read2) cout << "rofl " << read1 << endl;


		if(read1 == 9363 || read2 == 9363){
			nbOverlapsLala += 1;
			//cout << "ReadPairCount: " << nbOverlapsLala << " " << count << endl; 
		}


		//float dist = computeDistanceTNF(_readData[read1], _readData[read2]);
		//float t = 0.01f + count * 0.005f;
		//if(dist > t) continue; //To remove I think

		if(_evaluation_readToDataset[read1] != _evaluation_readToDataset[read2]){
			nbOverlaps_interGenomes += 1;
			//float dist = computeDistanceTNF(_readData[read1], _readData[read2]);
			//cout << "Bad overlaps: " <<  read1 << " " << read2 << " " << count << " " << dist << endl;
		}

		//if ((read1 == 7721 && read2 == 9442) || read2 == 7721 && read1 == 9442) {
			//cout << read_headers[read1] << " " << read_headers[read2] << endl;
		//}


		//if(read1 == 101) cout << _readData[read1]._graphVertex << endl;
		//if(read2 == 101) cout << _readData[read2]._graphVertex << endl;


		//output_file << read_headers[read1] << ";" << read_headers[read2] << ";" << count << endl; //!

		//cout << read1 << " " << read2 << " " << count << " " << (_evaluation_readToDataset[read1] == _evaluation_readToDataset[read2]) << endl;
		//if(!_readData[read1]._vertexCreated){
		//	_readData[read1]._vertexCreated = true;
			//nodes.push_back(read1);
			//boost::adjacency_list<>::vertex_descriptor v1 = boost::add_vertex(_overlapGraph);
			//_readData[read1]._graphVertex = v1;
			//_overlapGraph[v1]._readIndex = read1;

			//cout << "Create node: " << v1 << " " << read1 << endl;
			//if(read1 == 101) cout << "Added " << v1 << endl;

		//}
		//if(!_readData[read2]._vertexCreated){
		//	_readData[read2]._vertexCreated = true;
			//nodes.push_back(read2);
			//boost::adjacency_list<>::vertex_descriptor v1 = boost::add_vertex(_overlapGraph);
			//_readData[read2]._graphVertex = v1;
			//_overlapGraph[v1]._readIndex = read2;
			
			//cout << "Create node: " << v1 << " " << read2 << endl;
			//if(read2 == 101) cout << "Added " << v1 << endl;
		//}

		//if(read1 > nbNodes || read2 > nbNodes) cout << "allo" << endl;

		graphInfo->addNode(read1);
		graphInfo->addNode(read2);

		_overlapGraph->addEdge(graphInfo->readIndex_to_id(read1), graphInfo->readIndex_to_id(read2), 0);
		//cout << read1 << " " << read2 << " " << nbNodes << endl;

		//boost::add_edge(_readData[read1]._graphVertex, _readData[read2]._graphVertex, _overlapGraph);

		//output_file_gfa << "L" << "\t" << read1 << "\t" << "+" << "\t" << read2 << "\t" << "+" << "\t" << "0M" << endl;

		nbOverlaps += 1;

	}

	//vector<bool> hasMatch(readIndex, false);
	unordered_map<ReadIndexType, vector<ReadIndexType>> bestOverlaps;
	unordered_map<ReadIndexType, vector<ReadIndexType>> bestOverlaps_reverse;
	//cout << "----" << endl;
	std::sort(overlaps.begin(), overlaps.end(), OverlapComparator);
	for(Overlap& overlap : overlaps){

		ReadIndexType read1 = overlap._r1;
		ReadIndexType read2 = overlap._r2;


		if(bestOverlaps.find(read1) == bestOverlaps.end() && bestOverlaps_reverse.find(read2) == bestOverlaps_reverse.end()){

			bool isCycle = false;
			ReadIndexType currentRead = read1;
			while(true){

				if(bestOverlaps_reverse.find(currentRead) == bestOverlaps_reverse.end()) break;

				currentRead = bestOverlaps_reverse[currentRead];
				if(currentRead == read2){
					isCycle = true;
					break;
				}
			}

			if(!isCycle){
				bestOverlaps[read1] = read2;
				bestOverlaps_reverse[read2] = read1;
				continue;
			}
			
		}
		
		if(bestOverlaps.find(read2) == bestOverlaps.end() && bestOverlaps_reverse.find(read1) == bestOverlaps_reverse.end()){

			bool isCycle = false;
			ReadIndexType currentRead = read2;
			while(true){

				if(bestOverlaps_reverse.find(currentRead) == bestOverlaps_reverse.end()) break;
				currentRead = bestOverlaps_reverse[currentRead];
				if(currentRead == read1){
					isCycle = true;
					break;
				}
			}

			if(!isCycle){
				bestOverlaps[read2] = read1;
				bestOverlaps_reverse[read1] = read2;
			}
			
		}
		
		



		//else if(bestOverlaps.find(read2) == bestOverlaps.end()){
		//	bestOverlaps[read2] = read1;
		//	hasMatch[read1] = true;
		//}


		//cout << "Overlap: " << overlap._nbMinimizers << endl;
	}

	for (auto& it: bestOverlaps) {
		
		u_int32_t read1 = it.first;//(u_int32_t)((readPair & 0xFFFFFFFF00000000LL) >> 32);
		u_int32_t read2 = it.second;//(u_int32_t)(readPair & 0xFFFFFFFFLL);
		output_file_gfa << "L" << "\t" << read1 << "\t" << "+" << "\t" << read2 << "\t" << "+" << "\t" << "0M" << endl;
	}
	//for(size_t i=0; i<_readData.size(); i++){
	//	if(_readData[i]._composition.size() == 0) continue;
	//	output_file_gfa << "L" << "\t" << i << "\t" << "+" << "\t" << _readData[i]._readBestMatch << "\t" << "+" << "\t" << "0M" << endl;
	//}
	
	for(u_int32_t node : nodes){
		output_file_gfa << "S" << "\t" << node << "\t" << "*" << "\t" << "LN:i:" << _readData[node]._length << endl;
	}

	output_file_gfa.close();
	nodes.clear();

	cout << "Nb overlaps: " <<  nbOverlaps << endl;
	cout << "Nb overlaps inter genomes: " <<  nbOverlaps_interGenomes << endl;

	*/
	/*
	//_overlapGraph->display_AdjList(0, graphInfo);
	cout << "--------" << endl;
	vector<u_int32_t> neighbors;
	_overlapGraph->collectNeighbors(0, 2, neighbors);
	for(u_int32_t n : neighbors){
		cout << n << endl;
	}

	set<u_int32_t> lala;
	cout << "--------" << endl;
	_overlapGraph->collectNeighbors(0, 1, neighbors);
	for(u_int32_t n : neighbors){
		lala.insert(n);
		vector<u_int32_t> neighbors2;
		_overlapGraph->collectNeighbors(n, 1, neighbors2);
		for(u_int32_t nn : neighbors2){
			lala.insert(nn);
		}
	}

	for(auto lalala : lala){
		cout << lalala << endl;
	}

	//output_file.close();

	//_overlapGraph->display_AdjList(0);
	*/
	/*
	cout << "Nb nodes: " << boost::num_vertices(_overlapGraph) << endl;

	std::vector<int> component (boost::num_vertices (_overlapGraph));
	size_t num_components = boost::connected_components (_overlapGraph, &component[0]);
	cout << "Nb components: " << num_components << endl;

	vector<vector<ReadIndexType>> readPerCompoenents(num_components);
	vector<u_int64_t> nbNodePerComponents(num_components, 0);

	for (size_t i = 0; i < boost::num_vertices (_overlapGraph); i++){
		nbNodePerComponents[component[i]] += 1;

		readPerCompoenents[component[i]].push_back(_overlapGraph[i]._readIndex);
		//cout <<  "Node: " << i << " Read: " << graph[i]._readIndex << endl;
	}
	
	for(size_t i=0; i<num_components; i++){
		if(nbNodePerComponents[i] < 100) continue;
		cout << "Compoenents " << i << ": " << nbNodePerComponents[i] << endl;
	}

	cout << "Nb components: " << num_components << endl;
	*/
	/*
	for(size_t i=0; i<readPerCompoenents.size(); i++){
		cout << "---------------- " << i << endl;
		for(ReadIndexType ri : readPerCompoenents[i]){
			cout << _evaluation_readToDataset[ri] << endl;
		}
	}
	*/

	//createSimilarityGraph(graphInfo);


/*
static void writeGfa(AdjGraph* graph, const string& outputFilename){

        cout << "Dumping graph: " << outputFilename << endl;

        //cout << graph->_nbNodes << endl;
	    vector<u_int64_t> neighbors;

	    ofstream outputFile(outputFilename);

        for(size_t n=0; n<graph->_nbNodes; n++){
            //if(graph->_nodes[n]->isBidirection) continue;

            //cout << n << " " << graph->_graphInfo->_unitigs_length.size() << endl;
            outputFile << "S" << "\t" << n << "\t" << "*" << "\t" << "LN:i:" << graph->_graphInfo->_unitigs_length[n] << endl;


            adjNode* node = graph->_nodes[n];
            while (node != nullptr) {

                if(node->isBidirection){
                    node = node->next;
                    continue;
                }

                string from = "";
                if(node->directionFrom) from = "+"; else from = "-";

                string to = "";
                if(node->directionTo) to = "+"; else to = "-";

                outputFile << "L" << "\t" << n << "\t" << from << "\t" << node->val << "\t" << to << "\t" << "0M" << endl;

                node = node->next;

            }


        }

        outputFile.close();
    }

*/



/*
void Bloocoo::createSimilarityGraph(GraphInfo* graphInfo){


	
	ofstream stats_file("/home/gats/workspace/run/histos/binner_stats.csv");
	ofstream graph_output_file("/home/gats/workspace/run/histos/similarityGraph.txt");


	u_int64_t nbEdgeIntra = 0;
	u_int64_t nbEdgeInter = 0;

	vector<float> distTNF_intra;
	vector<float> distTNF_inter;
	vector<u_int32_t> neighbors;

	u_int64_t nbEdges = 0;

	for(size_t n=0; n<_overlapGraph->_nbNodes; n++){

		ReadIndexType read = graphInfo->id_to_readIndex(n);
		//ReadIndexType read = _overlapGraph[n]._readIndex;
		u_int32_t read_datasetId = _evaluation_readToDataset[read];

		//cout << n << " " << read << endl;
		_overlapGraph->collectNeighbors(n, 10, neighbors);

		//cout << "------" << endl;
		//auto neighbours = boost::adjacent_vertices(n, _overlapGraph);
		//for (auto nn : make_iterator_range(neighbours)){

		for (u_int64_t nn : neighbors){

			//cout << nn << endl;
			//ReadIndexType read_neighbor = _overlapGraph[nn]._readIndex;
			ReadIndexType read_neighbor = graphInfo->id_to_readIndex(nn);
			u_int32_t read_neighbor_datasetId = _evaluation_readToDataset[read_neighbor];

			

			float dist = computeDistanceTNF(_readData[read], _readData[read_neighbor]);

			//if(n == 0){
				//cout << dist << endl;
			//}
			if(read_datasetId == read_neighbor_datasetId){
				nbEdgeIntra += 1;
				distTNF_intra.push_back(dist);
				//cout << dist << endl;
			}
			else{
				nbEdgeInter += 1;
				distTNF_inter.push_back(dist);
				//cout << dist << endl;
			}
			
			if(dist > 0.05) continue;
			

			float weight = 1.0 - dist;
			graph_output_file << read << " " << read_neighbor << " " << weight << endl;
			nbEdges += 1;

		}
		
	}
	
	cout << "Nb edges total: " << nbEdges << endl;
	cout << "Nb edges intra: " << nbEdgeIntra << endl;
	cout << "Nb edges inter: " << nbEdgeInter << endl;
	cout << "Intra dists:" << endl;
	print_stats(distTNF_intra);
	cout << "Inter dists:" << endl;
	print_stats(distTNF_inter);

	graph_output_file.close();


	stats_file.close();
}
*/

void Bloocoo::extract_kminmers(){


	u_int64_t maxHashValue = -1;

	u_int64_t _hash_otpt[2];
	int _seed = 42;
	setDispatcher (new SerialDispatcher());


	IBank* inbank = Bank::open(_input_extractKminmers);

	
	Iterator<Sequence>* itSeq = createIterator<Sequence> (
														inbank->iterator(),
														inbank->estimateNbItems(),
														"Parsing reads"
														);

	LOCAL (itSeq);
		
	std::vector<Iterator<Sequence>*> itBanks =  itSeq->getComposition();
	u_int32_t readIndex = 0;
	u_int32_t datasetID = 0;

	//gzFile file = gzopen(_filename_hifiasmGroundtruth.c_str(),"wb");
	

	ModelCanonical model (_minimizerSize);
	ModelCanonical::Iterator itKmer (model);

	hash<KmerVec> h;

	for (size_t i=0; i<itBanks.size(); i++)
	{
		itSeq = createIterator<Sequence> (itBanks[i], inbank->estimateNbItemsBanki(i), "lala");

		u_int32_t position = 0;

		for (itSeq->first(); !itSeq->isDone(); itSeq->next()){


			Sequence& sequence = itSeq->item();


			char* readseq = sequence.getDataBuffer();
			string sequence_str;

			char lastChar = '0';
			for(size_t i=0; i<sequence.getDataSize(); i++){
				if(readseq[i] == lastChar) continue;
				sequence_str += readseq[i];
				lastChar = readseq[i];
			}


			size_t nbMinimizersPerRead = 0;

			Data buf((char*)sequence_str.c_str());



			itKmer.setData (buf);

			//u_int64_t lastMinimizer = -1;
			vector<u_int64_t> minimizers;
			//vector<u_int64_t> minimizers_pos;
			//u_int64_t nbMinizersRead = 0;

			//vector<MinimizerPair> minimizerPairs;
			

			//u_int64_t pos = 0;
			//u_int32_t lastMinimizerPos = -1;
			for (itKmer.first(); !itKmer.isDone(); itKmer.next()){

				kmer_type kmerMin = itKmer->value();
				u_int64_t kmerValue = kmerMin.getVal();
				u_int64_t minimizer;
				MurmurHash3_x64_128 ((const char*)&kmerValue, sizeof(kmerValue), _seed, &_hash_otpt);
				minimizer = _hash_otpt[0];



				//if(minimizerCounts[minimizer] > 1000) cout << minimizer << endl;
				double kmerHashed_norm = ((double) minimizer) / maxHashValue;
				if(kmerHashed_norm < _minimizerDensity){


					minimizers.push_back(minimizer);
					//minimizers_pos.push_back(pos);

					//cout << pos << endl;

					//minimizerCounts[minimizer] += 1;
					

				}

				//cout << kmerHashed << endl;
				//pos += 1;
			}

			
			int i_max = ((int)minimizers.size()) - (int)_kminmerSize + 1;
			for(int i=0; i<i_max; i++){

				KmerVec vec;

				for(int j=i; j<i+_kminmerSize; j++){
					u_int64_t minimizer = minimizers[j];

					vec._kmers.push_back(minimizer);
				}

				vec = vec.normalize();

				if(_evaluation_hifiasmGroundTruth.find(vec) != _evaluation_hifiasmGroundTruth.end()) continue;
				_evaluation_hifiasmGroundTruth[vec] = datasetID;
				_evaluation_hifiasmGroundTruth_position[vec] = position;
				position += 1;
				//cout << position << endl;
				//gzwrite(file, (const char*)&vec._kmers[0], _kminmerSize * sizeof(u_int64_t));
				//gzwrite(file, (const char*)&datasetID, sizeof(datasetID));
				
			}

			readIndex += 1;
		}


		datasetID += 1;
	}
	
	cout << "Nb minimizers groundtruth: " << _evaluation_hifiasmGroundTruth.size() << endl;
	//gzclose(file);

}