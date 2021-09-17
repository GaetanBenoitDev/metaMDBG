

#ifndef MDBG_METAG_TOBASESPACE
#define MDBG_METAG_TOBASESPACE

#include "../Commons.hpp"




class ToBasespace : public Tool{
    
public:

	string _inputFilename;
	string _inputDir;
	float _minimizerDensity;
    size_t _minimizerSize;
    size_t _kminmerSize;

	string _filename_outputContigs;
	string _filename_kminmerSequences;
	MDBG* _mdbg;
	
	unordered_map<u_int32_t, bool> _nodeName_entire;
	unordered_map<u_int32_t, bool> _nodeName_entire_rc;
	unordered_map<u_int32_t, bool> _nodeName_right;
	unordered_map<u_int32_t, bool> _nodeName_left;
	unordered_map<u_int32_t, bool> _nodeName_rightLast;
	unordered_map<u_int32_t, bool> _nodeName_leftLast;

	unordered_map<u_int32_t, string> _kminmerSequences_entire;
	unordered_map<u_int32_t, string> _kminmerSequences_left;
	unordered_map<u_int32_t, string> _kminmerSequences_right;

	enum LoadType{
        Entire,
        EntireRc,
        Left,
        Right,
        LeftLast,
        RightLast
	};

	ToBasespace(): Tool ("toBasepace"){

		getParser()->push_back (new OptionOneParam (STR_INPUT, "input file", true));
		getParser()->push_back (new OptionOneParam (STR_INPUT_DIR, "input dir", true));
		//getParser()->push_back (new OptionOneParam (STR_OUTPUT, "output contig filename in basespace", true));

	}

    void execute (){
		parseArgs();
		loadContigs();



		cout << "Loading mdbg" << endl;
		string mdbg_filename = _inputDir + "/mdbg_nodes.gz";
		_mdbg = new MDBG(_kminmerSize);
		_mdbg->load(mdbg_filename);
		cout << "MDBG nodes: " << _mdbg->_dbg_nodes.size() << endl;


		extractKminmerSequences();
		delete _mdbg;

		//exit(1);
		loadKminmerSequences();
		createBaseContigs();

		cout << endl << "Contig filename: " << _filename_outputContigs << endl;
	}

	void parseArgs(){
		_inputFilename = getInput()->getStr(STR_INPUT);
		_inputDir = getInput()->getStr(STR_INPUT_DIR);

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

		_filename_outputContigs = _inputDir + "/contigs.fasta.gz";
	}

	void loadContigs(){
		string contigFilename = _inputDir + "/minimizer_contigs.gz";

		gzFile contigFile = gzopen(contigFilename.c_str(),"rb");

		while(true){

			vector<u_int32_t> nodePath;
			u_int64_t size;
			gzread(contigFile, (char*)&size, sizeof(size));
			

			if(gzeof(contigFile)) break;

			nodePath.resize(size);
			gzread(contigFile, (char*)&nodePath[0], size * sizeof(u_int32_t));

			//for(u_int32_t nodeIndex : nodePath){
			for(size_t i=0; i<nodePath.size(); i++){
				u_int32_t nodeIndex = nodePath[i];
				bool orientation;
				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex, orientation);

				//cout << nodeName << endl;
				//_usedNodeNames.insert(nodeName);


				if(i == 0){
					_nodeName_entire[nodeName] = false;

					//if(orientation){ //+
					//	_nodeName_entire[nodeName] = false;
					//}
					//else{
					//	_nodeName_entire_rc[nodeName] = false;
					//}
				}
				else {
					if(orientation){
						_nodeName_right[nodeName] = false;
						//if(i == nodePath.size()-1){
						//	_nodeName_rightLast[nodeName] = false;
						//}
						//else{
						//	_nodeName_right[nodeName] = false;
						//}
					}
					else{
						_nodeName_left[nodeName] = false;
						//if(i == nodePath.size()-1){
						//	_nodeName_leftLast[nodeName] = false;
						//}
						//else{
						//	_nodeName_left[nodeName] = false;
						//}
					}
				}

			}
			//cout << nodePath.size() << endl;
		}

		gzclose(contigFile);

	}

	void extractKminmerSequences (){

		_filename_kminmerSequences = _inputDir + "/kminmerSequences";

		gzFile outputFile_entire = gzopen((_filename_kminmerSequences + "_entire.gz").c_str(),"wb");
		gzFile outputFile_entire_rc = gzopen((_filename_kminmerSequences + "_entire_rc.gz").c_str(),"wb");
		gzFile outputFile_right = gzopen((_filename_kminmerSequences + "_right.gz").c_str(),"wb");
		gzFile outputFile_rightLast = gzopen((_filename_kminmerSequences + "_rightLast.gz").c_str(),"wb");
		gzFile outputFile_left = gzopen((_filename_kminmerSequences + "_left.gz").c_str(),"wb");
		gzFile outputFile_leftLast = gzopen((_filename_kminmerSequences + "_leftLast.gz").c_str(),"wb");

		u_int64_t maxHashValue = -1;
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

		ModelCanonical model (_minimizerSize);
		ModelCanonical::Iterator itKmer (model);

		string kminmerSequence;

		for (size_t i=0; i<itBanks.size(); i++)
		{
			itSeq = createIterator<Sequence> (itBanks[i], inbank->estimateNbItemsBanki(i), "lala");

			for (itSeq->first(); !itSeq->isDone(); itSeq->next()){


				cout << readIndex << endl;




				Sequence& sequence = itSeq->item();

				string rleSequence;
				vector<u_int64_t> rlePositions;

				char* sequenceOriginal = sequence.getDataBuffer();
				Encoder::encode_rle(sequenceOriginal, sequence.getDataSize(), rleSequence, rlePositions);



				size_t nbMinimizersPerRead = 0;

				Data buf((char*)rleSequence.c_str());

				itKmer.setData (buf);

				u_int64_t lastMinimizer = -1;
				vector<u_int64_t> minimizers;
				vector<u_int64_t> minimizers_pos;
				u_int64_t nbMinizersRead = 0;

				vector<MinimizerPair> minimizerPairs;
				
				u_int64_t pos = 0;
				u_int32_t lastMinimizerPos = -1;
				for (itKmer.first(); !itKmer.isDone(); itKmer.next()){

					kmer_type kmerMin = itKmer->value();
					u_int64_t kmerValue = kmerMin.getVal();
					u_int64_t minimizer;
					MurmurHash3_x64_128 ((const char*)&kmerValue, sizeof(kmerValue), _seed, &_hash_otpt);
					minimizer = _hash_otpt[0];



					double kmerHashed_norm = ((double) minimizer) / maxHashValue;
					if(kmerHashed_norm < _minimizerDensity){


						minimizers.push_back(minimizer);
						minimizers_pos.push_back(pos); //rlePositions[pos]

						
						//minimizerCounts[minimizer] += 1;

					}

					pos += 1;
				}

				/*
				cout << endl << endl;
				for(size_t i=0; i<minimizers_pos.size(); i++){
					cout << minimizers_pos[i] << endl;
				}
				cout << endl << endl;
				*/

				//for(size_t i=0; i<rlePositions.size(); i++) {
				//	cout << i << ": " << rlePositions[i] << " " << endl;
				//}
				//cout << endl;

				vector<KmerVec> kminmers; 
				vector<ReadKminmer> kminmersInfo;
				MDBG::getKminmers(_minimizerSize, _kminmerSize, minimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions);

			/*	
		cout << endl;
		cout << endl;
		cout << endl;
		cout << endl;
		cout << endl;
		cout << endl;
		cout << endl;
				cout << endl;
				for(size_t i=0; i<kminmers.size(); i++){
					
					u_int32_t nodeName = _mdbg->_dbg_nodes[kminmers[i]]._index;
					ReadKminmer& kminmerInfo = kminmersInfo[i];
					
					extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Entire, kminmerSequence);
					//writeKminmerSequence(nodeName, kminmerSequence, outputFile_entire);

					//cout << kminmerInfo._isReversed << endl;
					//cout << kminmerSequence << endl;
				}
				
		cout << endl;
		cout << endl;
		cout << endl;
		cout << endl;
		cout << endl;
		cout << endl;
		cout << endl;
			*/
				//exit(1);
				
				for(size_t i=0; i<kminmers.size(); i++){

					//let read_offsets = (read_obj.minimizers_pos[i] as usize, (read_obj.minimizers_pos[i+k-1] as usize + l), (read_obj.minimizers_pos[i+k-1] + 1 - read_obj.minimizers_pos[i] + 1));

					u_int32_t nodeName = _mdbg->_dbg_nodes[kminmers[i]]._index;
					ReadKminmer& kminmerInfo = kminmersInfo[i];

					

					//if(_usedNodeNames.find(nodeName) == _usedNodeNames.end()) continue;

					auto found = _nodeName_entire.find(nodeName);
					if(found != _nodeName_entire.end() && !found->second){

						_nodeName_entire[nodeName] = true;
						extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Entire, kminmerSequence);
						writeKminmerSequence(nodeName, kminmerSequence, outputFile_entire);

						//cout << "sdfsdfs" << endl;
						//cout << nodeName << " " << kminmerSequence << endl;
						//cout << endl;
						//extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Entire, kminmerSequence);
						//cout << kminmerSequence << endl;
						//extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Left, kminmerSequence);
						//cout << kminmerSequence << endl;
						//extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Right, kminmerSequence);
						//cout << kminmerSequence << endl;


						//if(kminmerInfo.is)
						//outputFile_entire
						/*
						cout << endl;
						cout << endl;
						cout << nodeName << endl;
						cout << kminmers[i]._kmers[0] << endl;
						cout << kminmers[i]._kmers[1] << endl;
						cout << kminmers[i]._kmers[2] << endl;
						cout << endl;
						cout << endl;
					
						string lala =  string(sequenceOriginal); //string(sequenceOriginal); //rleSequence;//
						cout << lala.size() << " " << kminmerInfo._read_pos_start << " " << len << endl;
						cout << lala.substr(kminmerInfo._read_pos_start, len) << endl;

						
						string rleSequence2;
						vector<u_int64_t> rlePositions2;
						string loulou(lala.substr(kminmerInfo._read_pos_start, len));
						Encoder::encode_rle(loulou.c_str(), loulou.size(), rleSequence2, rlePositions2);
						cout << rleSequence2 << endl;

						for(size_t i=0; i<rlePositions2.size(); i++){
							cout << rlePositions2[i] << " ";
						}
						cout << endl;
						
						//exit(1);

						*/
						

						
						//cout << subbuff << " " << kminmerInfo._isReversed << endl;
						
						//cout << subbuff << endl << endl;

						/*
						string rleSequence2;
						vector<u_int64_t> rlePositions2;
						string loulou(subbuff);
						Encoder::encode_rle(loulou.c_str(), loulou.size(), rleSequence2, rlePositions2);
						cout << rleSequence2 << endl;
						*/
						//string kminmer_seq = sequenceOriginal.substr (kminmerInfo._read_pos_start, len);
						//cout << kminmerInfo._read_pos_start << " " << kminmer_seq << endl;
						//string sequence = sequenceOriginal[]
						//gzwrite(outputFile_entire, (const char*)&nodeName, sizeof(nodeName));
						//u_int16_t sequenceLength = ;
					}
					
					found = _nodeName_left.find(nodeName);
					if(found != _nodeName_left.end() && !found->second){
						
						_nodeName_left[nodeName] = true;
						extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Left, kminmerSequence);
						writeKminmerSequence(nodeName, kminmerSequence, outputFile_left);

						//extractKminmerSequence(sequenceOriginal, kminmerInfo, kminmerSequence);
						//u_int16_t len = 
						//string left_seq = kminmerSequence.substr(0, len);
						//writeKminmerSequence(nodeName, kminmerSequence, outputFile_entire);
						//let left_seq = utils::revcomp(&v[2][0..minim_pos[0] as usize]);
					}
					
					found = _nodeName_right.find(nodeName);
					if(found != _nodeName_right.end() && !found->second){
						
						_nodeName_right[nodeName] = true;
						extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Right, kminmerSequence);
						writeKminmerSequence(nodeName, kminmerSequence, outputFile_right);
					}



					//let seq = match seq {Some(read) => read, None => &read_seq.unwrap()[read_offsets.unwrap().0..read_offsets.unwrap().1]};
					//let seq = if *seq_reversed {utils::revcomp(&seq)} else {seq.to_string()};
					//let seq_line = format!("{}\t{}\t{}\t{}\t{}\t{:?}",cur_node_index, node.print_as_string(), seq, "*", origin, shift);

					/*
					for i in 0..(read_obj.transformed.len() - k + 1) {
						let mut node : Kmer = Kmer::make_from(&read_obj.transformed[i..i+k]);
						let mut seq_reversed = false;
						if REVCOMP_AWARE { 
							let (node_norm, reversed) = node.normalize(); 
							node = node_norm;
							seq_reversed = reversed;
						} 
						let origin = "*".to_string(); // uncomment the line below to track where the kmer is coming from (but not needed in production)
						let minimizers_pos = &read_obj.minimizers_pos;
						let position_of_second_minimizer = match seq_reversed {
							true => minimizers_pos[i+k-1] - minimizers_pos[i+k-2],
							false => minimizers_pos[i+1] - minimizers_pos[i]
						};
						let position_of_second_to_last_minimizer = match seq_reversed {
							true => minimizers_pos[i+1] - minimizers_pos[i],
							false => minimizers_pos[i+k-1] - minimizers_pos[i+k-2]
						};
						let shift = (position_of_second_minimizer, position_of_second_to_last_minimizer);
						let read_offsets = (read_obj.minimizers_pos[i] as usize, (read_obj.minimizers_pos[i+k-1] as usize + l), (read_obj.minimizers_pos[i+k-1] + 1 - read_obj.minimizers_pos[i] + 1));
					*/

					//u_int32_t nodeName = minimizers.size();
					//gzwrite(outputFile, (const char*)&minimizers[0], size * sizeof(u_int64_t));
					//kminmerCounts[kminmers[i]] += 1;
					//kminmersData[kminmers[i]] = kminmersLength[i];
				}

				//cout << "ENELEVR CE BREAK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
				//break;
				readIndex += 1;
			}


			datasetID += 1;
		}

		gzclose(outputFile_entire);
		gzclose(outputFile_entire_rc);
		gzclose(outputFile_right);
		gzclose(outputFile_rightLast);
		gzclose(outputFile_left);
		gzclose(outputFile_leftLast);

	}

	void extractKminmerSequence(const char* sequenceOriginal, const ReadKminmer& kminmerInfo, LoadType loadType, string& sequence){

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

		cout << "-------------------" << endl;
		cout << "-------------------" << endl;
		cout << "-------------------" << endl;
		cout << kminmerInfo._isReversed << endl;
		cout << sequence << endl;


		if(loadType == LoadType::Entire){
			return;
			//startPosition = kminmerInfo._read_pos_start;
			//len = kminmerInfo._read_pos_end - kminmerInfo._read_pos_start;
		}
		else if(loadType == LoadType::Left){
			startPosition = 0; //kminmerInfo._read_pos_start;
			len = kminmerInfo._seq_length_start; //kminmerInfo._position_of_second_minimizer_seq - kminmerInfo._read_pos_start;
			//cout << kminmerInfo._read_pos_start << " " << kminmerInfo._position_of_second_minimizer << endl;
		}
		else if(loadType == LoadType::Right){
			//return;
			startPosition = sequence.size() - kminmerInfo._seq_length_end;  //kminmerInfo._position_of_second_to_last_minimizer_seq;
			len = kminmerInfo._seq_length_end; //kminmerInfo._read_pos_end - kminmerInfo._position_of_second_to_last_minimizer_seq;
			//return;
			//cout << kminmerInfo._read_pos_end << " " << kminmerInfo._position_of_second_to_last_minimizer << endl;
		}
		
		//char* seq = sequence[0];
		//char subbuff2[len+1];
		//memcpy( subbuff, &seq[startPosition], len);
		//subbuff2[len] = '\0';
		//sequence = string(subbuff2);
		cout << sequence.size() << " " << startPosition << " " << len  << endl;
		sequence = sequence.substr(startPosition, len);


		//if(loadType == LoadType::Left){
			//Utils::revcomp(sequence);
			//exit(1);
		//}
		//if(loadType == LoadType::EntireRc){
		//	Utils::revcomp(sequence);
		//}

		//return len;



		cout << startPosition << " " << len << endl;
		cout << "allo" << endl;
		cout << endl << endl;
		//cout << 
		cout << sequence << endl;
		//Utils::revcomp(sequence);
		//cout << sequence << endl;

	}

	void writeKminmerSequence(u_int32_t nodeName, const string& sequence, const gzFile& file){
		
		u_int16_t length = sequence.size();

		cout << "Writing: " << endl;
		cout << nodeName << " " << length << " " << sequence << endl;
		gzwrite(file, (const char*)&nodeName, sizeof(nodeName));
		gzwrite(file, (const char*)&length, sizeof(length));
		gzwrite(file, (const char*)&sequence[0], length);
	}

	void loadKminmerSequences(){
		loadKminmerSequences_aux(_filename_kminmerSequences + "_entire.gz", _kminmerSequences_entire);
		loadKminmerSequences_aux(_filename_kminmerSequences + "_left.gz", _kminmerSequences_left);
		loadKminmerSequences_aux(_filename_kminmerSequences + "_right.gz", _kminmerSequences_right);
	}

	void loadKminmerSequences_aux(const string& filename, unordered_map<u_int32_t, string>& nodeSequences){
		
		gzFile file = gzopen(filename.c_str(), "rb");

		while(true){

			u_int32_t nodeName;
			u_int16_t sequenceLength;
			string sequence;
			gzread(file, (char*)&nodeName, sizeof(nodeName));
			
			if(gzeof(file)) break;

			gzread(file, (char*)&sequenceLength, sizeof(sequenceLength));
			sequence.resize(sequenceLength);
			gzread(file, (char*)&sequence[0], sequenceLength);

			//cout << nodeName << " " << sequenceLength << " " << sequence << endl;

			nodeSequences[nodeName] = sequence;
		}

		gzclose(file);
	}

	void createBaseContigs(){


		gzFile basespaceContigFile = gzopen(_filename_outputContigs.c_str(),"wb");


		cout << endl;
		cout << "Creating basespace contigs" << endl;
		cout << endl;

		string contigFilename = _inputDir + "/minimizer_contigs.gz";

		gzFile contigFile = gzopen(contigFilename.c_str(),"rb");

		u_int64_t contig_index = 0;

		while(true){

			vector<u_int32_t> nodePath;
			u_int64_t size;
			gzread(contigFile, (char*)&size, sizeof(size));
			

			if(gzeof(contigFile)) break;

			nodePath.resize(size);
			gzread(contigFile, (char*)&nodePath[0], size * sizeof(u_int32_t));


			string contigSequence = "";

			//for(u_int32_t nodeIndex : nodePath){
			for(size_t i=0; i<nodePath.size(); i++){
				cout << "|||||||||||||||||||||||||||||||||||||||||||||||" << endl;
				cout << endl << "Step: " << i << endl;
				u_int32_t nodeIndex = nodePath[i];
				bool orientation;
				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex, orientation);

				//_usedNodeNames.insert(nodeName);

				if(i == 0){
					//cout << nodeName << " " << orientation << endl;
					if(orientation){ //+
						cout << "Entire" << endl;
						contigSequence += _kminmerSequences_entire[nodeName];
						cout << contigSequence << endl;
					}
					else{
						cout << "Entire RC" << endl;
						string seq = _kminmerSequences_entire[nodeName];
						Utils::revcomp(seq);
						contigSequence += seq;
						cout << contigSequence << endl;
					}
				}
				else {
					if(orientation){
						cout << "right" << endl;
						contigSequence += _kminmerSequences_right[nodeName];
						cout << nodeName << endl;
						cout << _kminmerSequences_right[nodeName] << endl;
						
						if(i == nodePath.size()-1){
							//_nodeName_rightLast[nodeName] = false;
						}
						else{
							//_nodeName_right[nodeName] = false;
						}
						//cout << contigSequence << endl;
					}
					else{
						cout << "left" << endl;
						string seq = _kminmerSequences_left[nodeName];
						Utils::revcomp(seq);
						contigSequence += seq;

						//cout << "todo: revcomp je crois" << endl;
						cout << nodeName << endl;
						cout << seq << endl;
						//exit(1);
						if(i == nodePath.size()-1){
							//_nodeName_leftLast[nodeName] = false;
						}
						else{
							//_nodeName_left[nodeName] = false;
						}
					}
				}
				
				//cout << contigSequence << endl;

			}

			cout << endl << endl;
			cout << nodePath.size() << endl;
			cout << contigSequence.size() << endl;
			//cout << contigSequence << endl;

			string header = ">ctg" + to_string(contig_index) + '\n';
			gzwrite(basespaceContigFile, (const char*)&header[0], header.size());
			contigSequence +=  '\n';
			gzwrite(basespaceContigFile, (const char*)&contigSequence[0], contigSequence.size());

			contig_index += 1;
		}

		gzclose(contigFile);
		gzclose(basespaceContigFile);

	}

};	


#endif 



