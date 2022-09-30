

#ifndef MDBG_METAG_TOBASESPACE
#define MDBG_METAG_TOBASESPACE

#include "../Commons.hpp"
#include "../utils/edlib.h"
#include "../utils/spoa/include/spoa/spoa.hpp"
//#include <seqan/align.h>
//#include <seqan/graph_msa.h>
//#include <cstring>

class ToBasespace : public Tool{
    
public:

	string _inputFilename;
	string _inputDir;
	float _minimizerDensity;
    size_t _minimizerSize;
    size_t _kminmerSize;
	bool _isFirstPass;

	string _filename_outputContigs;
	string _filename_kminmerSequences;
	MDBG* _mdbg;
	MinimizerParser* _minimizerParser;
	
	unordered_map<ContigNode, string> _debug_node_sequences;
	unordered_map<u_int32_t, DnaBitset*> _nodeName_entire;
	unordered_map<u_int32_t, DnaBitset*> _nodeName_right;
	unordered_map<u_int32_t, DnaBitset*> _nodeName_left;

	//unordered_map<u_int32_t, DnaBitset*> _nodeName_all_right;
	//unordered_map<u_int32_t, DnaBitset*> _nodeName_all_left;
	//unordered_map<ContigNode, DnaBitset*> _kminmerSequences_entire;
	//unordered_map<ContigNode, DnaBitset*> _kminmerSequences_left;
	//unordered_map<ContigNode, DnaBitset*> _kminmerSequences_right;

	unordered_map<ContigNode, vector<DnaBitset*>> _kminmerSequenceCopies_entire;
	unordered_map<ContigNode, vector<DnaBitset*>> _kminmerSequenceCopies_left;
	unordered_map<ContigNode, vector<DnaBitset*>> _kminmerSequenceCopies_right;
	//unordered_map<u_int32_t, vector<DnaBitset*>> _kminmerSequenceCopies_left;
	//unordered_map<u_int32_t, vector<DnaBitset*>> _kminmerSequenceCopies_right;




	enum LoadType{
        Entire,
        EntireRc,
        Left,
        Right,
        LeftLast,
        RightLast
	};

	ToBasespace(): Tool (){

		/*
		getParser()->push_back (new OptionOneParam (STR_INPUT, "input file", true));
		getParser()->push_back (new OptionOneParam (STR_INPUT_DIR, "input dir", true));
		//getParser()->push_back (new OptionOneParam (STR_OUTPUT, "output contig filename in basespace", true));
		*/

	}

	gzFile _outputFile_left;
	gzFile _outputFile_right;
	//std::unique_ptr<spoa::AlignmentEngine> _alignment_engine;
	spoa::Graph _graph;

    void execute (){
		
		cout << "Loading mdbg" << endl;
		string mdbg_filename = _inputDir + "/mdbg_nodes.gz";
		_mdbg = new MDBG(_kminmerSize);
		_mdbg->load(mdbg_filename);
		cout << "MDBG nodes: " << _mdbg->_dbg_nodes.size() << endl;

		loadContigs(_inputDir + "/contigs.nodepath.gz");
		//loadContigs(_inputDir + "/minimizer_contigs_complete.gz");
		//loadContigs(_inputDir + "/eval/composition/22/debug_longUnitigs.gz");
		//loadContigs(_inputDir + "/eval/composition//3/debug_longUnitigsNeighbors.gz");

		cout << "TODO: remove ContigNode and _requiredCopiers_entire, instead use unorderedmap nodeName => vector<(ReadIndex, string)>" << endl;




		//extractKminmerSequences();
		lalalala();
		delete _mdbg;

		//exit(1);
		//loadKminmerSequences();
		createBaseContigs(_inputDir + "/contigs.nodepath.gz", _filename_outputContigs.c_str());
		//createBaseContigs(_inputDir + "/minimizer_contigs_complete.gz", _filename_outputContigs.c_str());
		//createBaseContigs(_inputDir + "/eval/composition//22//debug_longUnitigs.gz", _inputDir + "/eval/composition//22//debug_longUnitigs.fasta.gz");
		//createBaseContigs(_inputDir + "/eval/composition//3//debug_longUnitigsNeighbors.gz", _inputDir + "/eval/composition//3/debug_longUnitigsNeighbors.fasta.gz");

		cout << endl << "Contig filename: " << _filename_outputContigs << endl;
	}

	void removeVariants(){

	}

	/*
	void writeKminmerSequence(u_int32_t nodeName, unordered_map<u_int32_t, VariantQueue>& variants){

		_isDoneCorrection.insert()
		string correctedSequence;
		performErrorCorrection(nodeName, dnaSeq, _kminmerSequenceCopies_all_left[nodeName], correctedSequence, alignment_engine, graph);

		writeKminmerSequence_all(nodeName, correctedSequence, outputFile_left);


	}*/

	void parseArgs(int argc, char* argv[]){


		cxxopts::Options options("ToBasespace", "");
		options.add_options()
		(ARG_INPUT_FILENAME, "", cxxopts::value<string>())
		(ARG_FIRST_PASS, "", cxxopts::value<bool>()->default_value("false"))
		(ARG_OUTPUT_DIR, "", cxxopts::value<string>());

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

			_inputFilename = result[ARG_INPUT_FILENAME].as<string>();
			_inputDir = result[ARG_OUTPUT_DIR].as<string>();
			_isFirstPass = result[ARG_FIRST_PASS].as<bool>();
			
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
		_minimizerParser = new MinimizerParser(_minimizerSize, _minimizerDensity);
	}

	void loadContigs(const string& contigFilename){

		cout << "Loading mContigs: " << contigFilename << endl;
		gzFile contigFile = gzopen(contigFilename.c_str(),"rb");
		u_int64_t nbContigs = 0;
		
		while(true){

			vector<u_int32_t> nodePath;
			vector<u_int64_t> supportingReads;
			u_int64_t size;
			gzread(contigFile, (char*)&size, sizeof(size));
			

			if(gzeof(contigFile)) break;

			nodePath.resize(size);
			supportingReads.resize(size);
			gzread(contigFile, (char*)&nodePath[0], size * sizeof(u_int32_t));
			gzread(contigFile, (char*)&supportingReads[0], size * sizeof(u_int64_t));

			//for(u_int32_t nodeIndex : nodePath){
			for(size_t i=0; i<nodePath.size(); i++){
				u_int32_t nodeIndex = nodePath[i];
				bool orientation;
				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex, orientation);
				u_int32_t readIndex = supportingReads[i];

				//ContigNode contigNode = {nodeName, supportingReads[i]};

				//cout << nodeName << " " << supportingReads[i] << endl;

				//cout << nodeName << endl;
				//_usedNodeNames.insert(nodeName);


				if(i == 0){
					initKminmerSequence(nodeName, readIndex, _nodeName_entire, _nodeName_entire_multi, _requiredCopiers_entire);

					//_nodeName_entire[nodeName] = {0, nullptr};
					//_requiredCopiers_entire.insert(nodeName);

					//if(orientation){ //+
					//	_nodeName_entire[nodeName] = false;
					//}
					//else{
					//	_nodeName_entire_rc[nodeName] = false;
					//}
				}
				else {
					if(orientation){ //+
						initKminmerSequence(nodeName, readIndex, _nodeName_right, _nodeName_right_multi, _requiredCopiers_right);
						//_nodeName_right[nodeName] = {0, nullptr};
						//_requiredCopiers_right.insert(nodeName);
						//if(i == nodePath.size()-1){
						//	_nodeName_rightLast[nodeName] = false;
						//}
						//else{
						//	_nodeName_right[nodeName] = false;
						//}
					}
					else{ //-
						initKminmerSequence(nodeName, readIndex, _nodeName_left, _nodeName_left_multi, _requiredCopiers_left);
						//_nodeName_left[nodeName] = {0, nullptr};
						//_requiredCopiers_left.insert(nodeName);
						//if(i == nodePath.size()-1){
						//	_nodeName_leftLast[nodeName] = false;
						//}
						//else{
						//	_nodeName_left[nodeName] = false;
						//}
					}
				}

			}

			nbContigs += 1;
			//cout << nodePath.size() << endl;
		}

		gzclose(contigFile);

		cout << "Nb contigs: " << nbContigs << endl;
	}

	void initKminmerSequence(u_int32_t nodeName, u_int32_t readIndex, auto& dictSimple, auto& dictMulti, auto& requiredCopies){


		if(dictMulti.find(nodeName) != dictMulti.end()){
			dictMulti[nodeName].push_back({readIndex, nullptr});

			//cout << "----" << endl;
			//for(KminmerSequence& seq : dictMulti[nodeName]){
			//	cout << seq._readIndex << endl;
			//}
			return;
		}

		if(dictSimple.find(nodeName) != dictSimple.end()){
			dictMulti[nodeName].push_back({dictSimple[nodeName]._readIndex, nullptr});
			dictSimple.erase(nodeName);
			dictMulti[nodeName].push_back({readIndex, nullptr});
			return;
		}

		dictSimple[nodeName] = {readIndex, nullptr};
		requiredCopies.insert(nodeName);
	}


	void extractKminmerSequences_read(kseq_t* read, u_int64_t readIndex){
		//cout << readIndex << endl;
		//ottalSize += strlen(read->seq.s);
					
		string kminmerSequence;
		char* sequenceOriginal = read->seq.s;

		string rleSequence;
		vector<u_int64_t> rlePositions;
		Encoder::encode_rle(read->seq.s, strlen(read->seq.s), rleSequence, rlePositions);

		vector<u_int64_t> minimizers;
		vector<u_int64_t> minimizers_pos;
		_minimizerParser->parse(rleSequence, minimizers, minimizers_pos);

		vector<KmerVec> kminmers; 
		vector<ReadKminmer> kminmersInfo;
		MDBG::getKminmers(_minimizerSize, _kminmerSize, minimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, readIndex, false);

		for(size_t i=0; i<kminmers.size(); i++){
			if(_mdbg->_dbg_nodes.find(kminmers[i]) == _mdbg->_dbg_nodes.end()) continue;

			u_int32_t nodeName = _mdbg->_dbg_nodes[kminmers[i]]._index;
			ReadKminmer& kminmerInfo = kminmersInfo[i];


			//ContigNode contigNode = {nodeName, readIndex};
	

			auto found = _nodeName_entire.find(nodeName);
			if(found != _nodeName_entire.end() && found->second._readIndex == readIndex){


				extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Entire, kminmerSequence);
				_nodeName_entire[nodeName] = {readIndex, new DnaBitset(kminmerSequence)};

			}
			else{
				if(_nodeName_entire_multi.find(nodeName) != _nodeName_entire_multi.end()){
					for(KminmerSequence& seq : _nodeName_entire_multi[nodeName]){
						if(seq._readIndex == readIndex){
							extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Entire, kminmerSequence);
							seq._sequence = new DnaBitset(kminmerSequence);
						}
					}
				}
			}
			



			found = _nodeName_left.find(nodeName);
			if(found != _nodeName_left.end() && found->second._readIndex == readIndex){
				
				extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Left, kminmerSequence);
				_nodeName_left[nodeName] = {readIndex, new DnaBitset(kminmerSequence)};

			}
			else{
				if(_nodeName_left_multi.find(nodeName) != _nodeName_left_multi.end()){
					for(KminmerSequence& seq : _nodeName_left_multi[nodeName]){
						if(seq._readIndex == readIndex){
							extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Left, kminmerSequence);
							seq._sequence = new DnaBitset(kminmerSequence);
						}
					}
				}
			}
			

			found = _nodeName_right.find(nodeName);
			if(found != _nodeName_right.end() && found->second._readIndex == readIndex){
				
				extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Right, kminmerSequence);
				_nodeName_right[nodeName] = {readIndex, new DnaBitset(kminmerSequence)};
			}
			else{
				if(_nodeName_right_multi.find(nodeName) != _nodeName_right_multi.end()){
					for(KminmerSequence& seq : _nodeName_right_multi[nodeName]){
						if(seq._readIndex == readIndex){
							extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Right, kminmerSequence);
							seq._sequence = new DnaBitset(kminmerSequence);
						}
					}
				}
			}


				
		}
	}

	void extractKminmerSequences (){

		auto fp = std::bind(&ToBasespace::extractKminmerSequences_read, this, std::placeholders::_1, std::placeholders::_2);
		ReadParser readParser(_inputFilename, false, !_isFirstPass);
		readParser.parse(fp);


		/*
		//_filename_kminmerSequences = _inputDir + "/kminmerSequences";

		//gzFile outputFile_entire = gzopen((_filename_kminmerSequences + "_entire.gz").c_str(),"wb");
		//gzFile outputFile_entire_rc = gzopen((_filename_kminmerSequences + "_entire_rc.gz").c_str(),"wb");
		//gzFile outputFile_right = gzopen((_filename_kminmerSequences + "_right.gz").c_str(),"wb");
		//gzFile outputFile_rightLast = gzopen((_filename_kminmerSequences + "_rightLast.gz").c_str(),"wb");
		//gzFile outputFile_left = gzopen((_filename_kminmerSequences + "_left.gz").c_str(),"wb");
		//gzFile outputFile_leftLast = gzopen((_filename_kminmerSequences + "_leftLast.gz").c_str(),"wb");

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
		u_int64_t readIndex = 0;
		u_int32_t datasetID = 0;

		ModelCanonical model (_minimizerSize);
		ModelCanonical::Iterator itKmer (model);

		string kminmerSequence;

		for (size_t i=0; i<itBanks.size(); i++)
		{
			itSeq = createIterator<Sequence> (itBanks[i], inbank->estimateNbItemsBanki(i), "lala");

					bool lala = false;
					bool loulou = false;
			for (itSeq->first(); !itSeq->isDone(); itSeq->next()){

				//cout << readIndex << endl;





				Sequence& sequence = itSeq->item();
				//if(readIndex == 593){
				//	cout << sequence.toString() << endl;
				//}


				string rleSequence;
				vector<u_int64_t> rlePositions;

				char* sequenceOriginal = sequence.getDataBuffer();
				Encoder::encode_rle(sequenceOriginal, sequence.getDataSize(), rleSequence, rlePositions);




				vector<u_int64_t> minimizers;
				vector<u_int64_t> minimizers_pos;
				u_int64_t pos = 0;

				ntHashIterator ntHashIt(rleSequence, 1, _minimizerSize);

				while (ntHashIt != ntHashIt.end()) {

					if(pos == 0){
						pos += 1;
						++ntHashIt;
						continue;
					}
					else if(pos == rleSequence.size()-_minimizerSize){
						++ntHashIt;
						continue;
					}

					u_int64_t minimizer = (*ntHashIt)[0];

					double kmerHashed_norm = ((double) minimizer) / maxHashValue;
					if(kmerHashed_norm < (_minimizerDensity*0.5)){
						minimizers.push_back(minimizer);
						minimizers_pos.push_back(pos);

						//cout << pos << endl;
						//cout << rlePositions[pos] << endl; 
						
						//minimizerCounts[minimizer] += 1;
					}

					//bloom.insert(*itr);
					++ntHashIt;
					pos += 1;
				}
				
			

				vector<KmerVec> kminmers; 
				vector<ReadKminmer> kminmersInfo;
				MDBG::getKminmers(_minimizerSize, _kminmerSize, minimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, readIndex);

			
				//exit(1);
				for(size_t i=0; i<kminmers.size(); i++){

					if(_mdbg->_dbg_nodes.find(kminmers[i]) == _mdbg->_dbg_nodes.end()) continue;
					//let read_offsets = (read_obj.minimizers_pos[i] as usize, (read_obj.minimizers_pos[i+k-1] as usize + l), (read_obj.minimizers_pos[i+k-1] + 1 - read_obj.minimizers_pos[i] + 1));

					u_int32_t nodeName = _mdbg->_dbg_nodes[kminmers[i]]._index;
					ReadKminmer& kminmerInfo = kminmersInfo[i];




					ContigNode contigNode = {nodeName, readIndex};
			
				

					auto found = _nodeName_entire.find(nodeName);
					if(found != _nodeName_entire.end() && found->second._readIndex == readIndex){


						extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Entire, kminmerSequence);
						_nodeName_entire[nodeName] = {readIndex, new DnaBitset(kminmerSequence)};

					}
					else{
						if(_nodeName_entire_multi.find(nodeName) != _nodeName_entire_multi.end()){
							for(KminmerSequence& seq : _nodeName_entire_multi[nodeName]){
								if(seq._readIndex == readIndex){
									extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Entire, kminmerSequence);
									seq._sequence = new DnaBitset(kminmerSequence);
								}
							}
						}
					}
					



					found = _nodeName_left.find(nodeName);
					if(found != _nodeName_left.end() && found->second._readIndex == readIndex){
						
						extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Left, kminmerSequence);
						_nodeName_left[nodeName] = {readIndex, new DnaBitset(kminmerSequence)};

					}
					else{
						if(_nodeName_left_multi.find(nodeName) != _nodeName_left_multi.end()){
							for(KminmerSequence& seq : _nodeName_left_multi[nodeName]){
								if(seq._readIndex == readIndex){
									extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Left, kminmerSequence);
									seq._sequence = new DnaBitset(kminmerSequence);
								}
							}
						}
					}
					

					found = _nodeName_right.find(nodeName);
					if(found != _nodeName_right.end() && found->second._readIndex == readIndex){
						
						extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Right, kminmerSequence);
						_nodeName_right[nodeName] = {readIndex, new DnaBitset(kminmerSequence)};
					}
					else{
						if(_nodeName_right_multi.find(nodeName) != _nodeName_right_multi.end()){
							for(KminmerSequence& seq : _nodeName_right_multi[nodeName]){
								if(seq._readIndex == readIndex){
									extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Right, kminmerSequence);
									seq._sequence = new DnaBitset(kminmerSequence);
								}
							}
						}
					}


				}

				readIndex += 1;
			}


			datasetID += 1;
		}

		//gzclose(outputFile_entire);
		//gzclose(outputFile_entire_rc);
		//gzclose(outputFile_right);
		//gzclose(outputFile_rightLast);
		//gzclose(outputFile_left);
		//gzclose(outputFile_leftLast);
		*/
	}



	void lalalala_read(kseq_t* read, u_int64_t readIndex){
		//ottalSize += strlen(read->seq.s);
					
		string kminmerSequence;
		char* sequenceOriginal = read->seq.s;

		string rleSequence;
		vector<u_int64_t> rlePositions;
		Encoder::encode_rle(read->seq.s, strlen(read->seq.s), rleSequence, rlePositions);

		vector<u_int64_t> minimizers;
		vector<u_int64_t> minimizers_pos;
		_minimizerParser->parse(rleSequence, minimizers, minimizers_pos);

		vector<KmerVec> kminmers; 
		vector<ReadKminmer> kminmersInfo;
		MDBG::getKminmers(_minimizerSize, _kminmerSize, minimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, readIndex, false);

		for(size_t i=0; i<kminmers.size(); i++){
			if(_mdbg->_dbg_nodes.find(kminmers[i]) == _mdbg->_dbg_nodes.end()) continue;

			u_int32_t nodeName = _mdbg->_dbg_nodes[kminmers[i]]._index;
			ReadKminmer& kminmerInfo = kminmersInfo[i];


			if(_requiredCopiers_entire.find(nodeName) != _requiredCopiers_entire.end()){
				extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Entire, kminmerSequence);
				addKminmerSequenceVariant(nodeName, readIndex, _nodeName_entire, _nodeName_entire_multi, _kminmerSequenceCopies_entire, kminmerSequence);
				//_kminmerSequenceCopies_entire[nodeName].push_back(new DnaBitset(kminmerSequence));
			}
			if(_requiredCopiers_left.find(nodeName) != _requiredCopiers_left.end()){//} && _kminmerSequenceCopies_left[nodeName].size() < 50){
				extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Left, kminmerSequence);
				addKminmerSequenceVariant(nodeName, readIndex, _nodeName_left, _nodeName_left_multi, _kminmerSequenceCopies_left, kminmerSequence);
				//_kminmerSequenceCopies_left[nodeName].push_back(new DnaBitset(kminmerSequence));
			}
			if(_requiredCopiers_right.find(nodeName) != _requiredCopiers_right.end()){// && _kminmerSequenceCopies_right[nodeName].size() < 50){
				extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Right, kminmerSequence);
				addKminmerSequenceVariant(nodeName, readIndex, _nodeName_right, _nodeName_right_multi, _kminmerSequenceCopies_right, kminmerSequence);
				//_kminmerSequenceCopies_right[nodeName].push_back(new DnaBitset(kminmerSequence));
			}


				
		}
	}

	void lalalala (){

		auto fp = std::bind(&ToBasespace::lalalala_read, this, std::placeholders::_1, std::placeholders::_2);
		ReadParser readParser(_inputFilename, false, !_isFirstPass);
		readParser.parse(fp);
		/*
		//_filename_kminmerSequences = _inputDir + "/kminmerSequences";


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
		u_int64_t readIndex = 0;
		u_int32_t datasetID = 0;

		ModelCanonical model (_minimizerSize);
		ModelCanonical::Iterator itKmer (model);

		string kminmerSequence;

		for (size_t i=0; i<itBanks.size(); i++)
		{
			itSeq = createIterator<Sequence> (itBanks[i], inbank->estimateNbItemsBanki(i), "lala");

			for (itSeq->first(); !itSeq->isDone(); itSeq->next()){

				Sequence& sequence = itSeq->item();

				string rleSequence;
				vector<u_int64_t> rlePositions;

				char* sequenceOriginal = sequence.getDataBuffer();
				Encoder::encode_rle(sequenceOriginal, sequence.getDataSize(), rleSequence, rlePositions);


				vector<u_int64_t> minimizers;
				vector<u_int64_t> minimizers_pos;
				u_int64_t pos = 0;

				ntHashIterator ntHashIt(rleSequence, 1, _minimizerSize);

				while (ntHashIt != ntHashIt.end()) {

					if(pos == 0){
						pos += 1;
						++ntHashIt;
						continue;
					}
					else if(pos == rleSequence.size()-_minimizerSize){
						++ntHashIt;
						continue;
					}

					u_int64_t minimizer = (*ntHashIt)[0];

					double kmerHashed_norm = ((double) minimizer) / maxHashValue;
					if(kmerHashed_norm < (_minimizerDensity*0.5)){
						minimizers.push_back(minimizer);
						minimizers_pos.push_back(pos);

						//cout << pos << endl;
						//cout << rlePositions[pos] << endl; 
						
						//minimizerCounts[minimizer] += 1;
					}

					//bloom.insert(*itr);
					++ntHashIt;
					pos += 1;
				}

		

				vector<KmerVec> kminmers; 
				vector<ReadKminmer> kminmersInfo;
				MDBG::getKminmers(_minimizerSize, _kminmerSize, minimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, readIndex);

				for(size_t i=0; i<kminmers.size(); i++){

					if(_mdbg->_dbg_nodes.find(kminmers[i]) == _mdbg->_dbg_nodes.end()) continue;
					//let read_offsets = (read_obj.minimizers_pos[i] as usize, (read_obj.minimizers_pos[i+k-1] as usize + l), (read_obj.minimizers_pos[i+k-1] + 1 - read_obj.minimizers_pos[i] + 1));

					u_int32_t nodeName = _mdbg->_dbg_nodes[kminmers[i]]._index;
					ReadKminmer& kminmerInfo = kminmersInfo[i];




					ContigNode contigNode = {nodeName, readIndex};

					if(_requiredCopiers_entire.find(nodeName) != _requiredCopiers_entire.end()){
						extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Entire, kminmerSequence);
						addKminmerSequenceVariant(nodeName, readIndex, _nodeName_entire, _nodeName_entire_multi, _kminmerSequenceCopies_entire, kminmerSequence);
						//_kminmerSequenceCopies_entire[nodeName].push_back(new DnaBitset(kminmerSequence));
					}
					if(_requiredCopiers_left.find(nodeName) != _requiredCopiers_left.end()){//} && _kminmerSequenceCopies_left[nodeName].size() < 50){
						extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Left, kminmerSequence);
						addKminmerSequenceVariant(nodeName, readIndex, _nodeName_left, _nodeName_left_multi, _kminmerSequenceCopies_left, kminmerSequence);
						//_kminmerSequenceCopies_left[nodeName].push_back(new DnaBitset(kminmerSequence));
					}
					if(_requiredCopiers_right.find(nodeName) != _requiredCopiers_right.end()){// && _kminmerSequenceCopies_right[nodeName].size() < 50){
						extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Right, kminmerSequence);
						addKminmerSequenceVariant(nodeName, readIndex, _nodeName_right, _nodeName_right_multi, _kminmerSequenceCopies_right, kminmerSequence);
						//_kminmerSequenceCopies_right[nodeName].push_back(new DnaBitset(kminmerSequence));
					}
				}

				readIndex += 1;
			}


			datasetID += 1;
		}

		*/
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

		/*
		cout << "-------------------" << endl;
		cout << "-------------------" << endl;
		cout << "-------------------" << endl;
		cout << kminmerInfo._isReversed << endl;
		cout << sequence << endl;
		*/

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
		/*
		cout << sequence.size() << " " << startPosition << " " << len  << endl;
		*/
		sequence = sequence.substr(startPosition, len);


		//if(loadType == LoadType::Left){
			//Utils::revcomp(sequence);
			//exit(1);
		//}
		//if(loadType == LoadType::EntireRc){
		//	Utils::revcomp(sequence);
		//}

		//return len;


		/*
		cout << startPosition << " " << len << endl;
		cout << "allo" << endl;
		cout << endl << endl;
		//cout << 
		cout << sequence << endl;
		//Utils::revcomp(sequence);
		//cout << sequence << endl;
		*/
	}

	void addKminmerSequenceVariant(u_int32_t nodeName, u_int32_t readIndex, auto& dictSimple, auto& dictMulti, auto& variants, const string& sequence){
		
		
		if(dictSimple.find(nodeName) != dictSimple.end()){
			if(dictSimple[nodeName]._readIndex == readIndex) return; //is model
			addKminmerSequenceVariant_add(nodeName, dictSimple[nodeName]._readIndex, variants, dictSimple[nodeName]._sequence, sequence);
		}
		else{
			if(dictMulti.find(nodeName) != dictMulti.end()){
				for(KminmerSequence& seq : dictMulti[nodeName]){
					if(seq._readIndex == readIndex) continue; //is model
					addKminmerSequenceVariant_add(nodeName, seq._readIndex, variants, seq._sequence, sequence);
				}
			}
		}

		//bool isSimple;
		//DnaBitset* model = getKminmerSequence(nodeName, readIndex, dictSimple, dictMulti, isSimple);
		//if(model == nullptr) return;

		//cout << (models.find(node) == models.end()) << " " << (models[node] == nullptr) << endl;
		//if(models.find(node) == models.end() || models[node] == nullptr) return;


		//if(variants.find(node) == variants.end()){
			//variants[node] = new VariantQueue();
		//}
		//14256 3725
		//cout << nodeName << " " << readIndex << endl;


		
	}

	void addKminmerSequenceVariant_add(u_int32_t nodeName, u_int32_t readIndex, auto& variants, DnaBitset* sequenceModel, const string& sequence){
		if(sequenceModel == nullptr){
			cout << "pas normal" << endl;
			return; //model not found yet
		}


		ContigNode contigNode = {nodeName, readIndex};
		VariantQueue& queue = variants[contigNode];

		
		static EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0);


		char* sequenceModelStr = sequenceModel->to_string();

		EdlibAlignResult result = edlibAlign(sequenceModelStr, sequenceModel->m_len, sequence.c_str(), sequence.size(), config);
		free(sequenceModelStr);

		if (result.status != EDLIB_STATUS_OK){
			edlibFreeAlignResult(result);
			return;
		}
		
		//queue.push({result.editDistance, new DnaBitset(sequence)});


		//edlibFreeAlignResult(result);
		/*
		if(nodeName == 17326){
			cout << "------------" << endl;
			cout << result.editDistance << endl;

			//cout << sequenceModelStr << endl;
			//cout << sequence << endl;
		}
		*/

		//queue.push({result.editDistance, new DnaBitset(sequence)});

		
		
		if(queue.size() < 20){
			queue.push({result.editDistance, new DnaBitset(sequence)});
		}
		else{
			if(result.editDistance < queue.top()._editDistance){

				const KminmerSequenceVariant& variant = queue.top();
				delete variant._sequence;

				queue.pop();
				queue.push({result.editDistance, new DnaBitset(sequence)});
			}
		}
		
		//cout << queue.size() << endl;
		edlibFreeAlignResult(result);
		//if(nodeName == 17326){
			//cout << queue.top()._editDistance << endl;
		//}
		
	}


	/*
	void loadKminmerSequences_aux(const string& filename, unordered_map<ContigNode, string>& nodeSequences){
		
		gzFile file = gzopen(filename.c_str(), "rb");

		while(true){

			u_int32_t nodeName;
			u_int64_t readIndex;
			u_int16_t sequenceLength;
			string sequence;
			gzread(file, (char*)&nodeName, sizeof(nodeName));
			
			if(gzeof(file)) break;

			gzread(file, (char*)&readIndex, sizeof(readIndex));
			gzread(file, (char*)&sequenceLength, sizeof(sequenceLength));
			sequence.resize(sequenceLength);
			gzread(file, (char*)&sequence[0], sequenceLength);

			//cout << nodeName << " " << sequenceLength << " " << sequence << endl;

			nodeSequences[{nodeName, readIndex}] = sequence;
		}

		gzclose(file);
	}
	*/

	DnaBitset* getKminmerSequence(u_int32_t nodeName, u_int32_t readIndex, auto& dictSimple, auto& dictMulti){
		bool dummy;
		return getKminmerSequence(nodeName, readIndex, dictSimple, dictMulti, dummy);
	}

	DnaBitset* getKminmerSequence(u_int32_t nodeName, u_int32_t readIndex, auto& dictSimple, auto& dictMulti, bool& isSimple){
		if(dictSimple.find(nodeName) != dictSimple.end()){
			isSimple = true;
			return dictSimple[nodeName]._sequence;
		}

		for(KminmerSequence& seq : dictMulti[nodeName]){
			if(seq._readIndex == readIndex){
				isSimple = false;
				return seq._sequence;
			}
		}

		return nullptr;
	}

	void createBaseContigs(const string& contigFilename, const string& outputFilename){

		cout << "Creating basespace contigs: " << contigFilename << " " << outputFilename << endl;

		auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -3);  // linear gaps
		spoa::Graph graph{};


		gzFile basespaceContigFile = gzopen(outputFilename.c_str(),"wb");


		cout << endl;
		cout << "Creating basespace contigs" << endl;
		cout << endl;


		gzFile contigFile = gzopen(contigFilename.c_str(),"rb");

		u_int64_t contig_index = 0;

		while(true){

			vector<u_int32_t> nodePath;
			vector<u_int64_t> supportingReads;
			u_int64_t size;
			gzread(contigFile, (char*)&size, sizeof(size));
			

			if(gzeof(contigFile)) break;

			nodePath.resize(size);
			supportingReads.resize(size);
			gzread(contigFile, (char*)&nodePath[0], size * sizeof(u_int32_t));
			gzread(contigFile, (char*)&supportingReads[0], size * sizeof(u_int64_t));


			string contigSequence = "";

			//for(u_int32_t nodeIndex : nodePath){
			for(size_t i=0; i<nodePath.size(); i++){
				//cout << "|||||||||||||||||||||||||||||||||||||||||||||||" << endl;
				//cout << endl << "Step: " << i << endl;
				u_int32_t nodeIndex = nodePath[i];
				bool orientation;
				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex, orientation);
				u_int64_t readIndex = supportingReads[i];
				ContigNode contigNode = {nodeName, readIndex};

				//_usedNodeNames.insert(nodeName);

				if(i == 0){
					//cout << nodeName << " " << orientation << endl;
					if(orientation){ //+
						cout << "Entire" << endl;

						//string seq = _nodeName_entire[contigNode];
						string correctedSequence;
						performErrorCorrection(nodeName, getKminmerSequence(nodeName, readIndex, _nodeName_entire, _nodeName_entire_multi), _kminmerSequenceCopies_entire[contigNode], correctedSequence, alignment_engine, graph);
						

						contigSequence += correctedSequence;
						//cout << contigSequence << endl;
					}
					else{
						cout << "Entire RC" << endl;
						//string seq = ;
						string correctedSequence;
						performErrorCorrection(nodeName, getKminmerSequence(nodeName, readIndex, _nodeName_entire, _nodeName_entire_multi), _kminmerSequenceCopies_entire[contigNode], correctedSequence, alignment_engine, graph);
						

						Utils::revcomp(correctedSequence);
						contigSequence += correctedSequence;
						//cout << contigSequence << endl;
					}
				}
				else {
					if(orientation){
						//string seq = ;
						//contigSequence += seq;

						//cout << _kminmerSequenceCopies_right[nodeName].size() << endl;
						
						string correctedSequence;
						performErrorCorrection(nodeName, getKminmerSequence(nodeName, readIndex, _nodeName_right, _nodeName_right_multi), _kminmerSequenceCopies_right[contigNode], correctedSequence, alignment_engine, graph);
						contigSequence += correctedSequence;
						
						//cout << contigSequence.size() << endl;
						//cout << nodeName << endl;
						//cout << _kminmerSequences_right[contigNode] << endl;
						
				/*
				if(contigSequence.size() >= 362100 && contigSequence.size() < 363100){
						cout << "right" << endl;
						cout << nodeName << endl;
					cout << seq << endl;
					//cout << nodeName << endl;
					//exit(1);
				}
										if(i == nodePath.size()-1){
							//_nodeName_rightLast[nodeName] = false;
						}
						else{
							//_nodeName_right[nodeName] = false;
						}*/
						//cout << contigSequence << endl;
					}
					else{
						//string seq = _nodeName_left[contigNode];
						//Utils::revcomp(seq);
						//contigSequence += seq;
						
						//cout << _kminmerSequenceCopies_left[nodeName].size() << endl;
						
						string correctedSequence;
						performErrorCorrection(nodeName, getKminmerSequence(nodeName, readIndex, _nodeName_left, _nodeName_left_multi), _kminmerSequenceCopies_left[contigNode], correctedSequence, alignment_engine, graph);
						
						Utils::revcomp(correctedSequence);
						contigSequence += correctedSequence;
						//cout << contigSequence.size() << endl;
/*
				if(contigSequence.size() >= 362100 && contigSequence.size() < 363100){
						cout << "left" << endl;
						cout << nodeName << endl;
					cout << seq << endl;
					//cout << nodeName << endl;
					//exit(1);
				}

						//cout << "todo: revcomp je crois" << endl;
						//cout << nodeName << endl;
						//cout << seq << endl;
						//exit(1);
						if(i == nodePath.size()-1){
							//_nodeName_leftLast[nodeName] = false;
						}
						else{
							//_nodeName_left[nodeName] = false;
						}*/
					}
				}
				
				

				//cout << contigSequence.size() << endl;
				
				//cout << contigSequence << endl;

			}

			//cout << contigSequence.substr(362170, 100) << endl;
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

	struct AlignScore{
		u_int32_t _sequenceIndex;
		u_int32_t _editDistance;
	};

	struct AlignScore_Comparator {
		bool operator()(AlignScore const& p1, AlignScore const& p2)
		{
			return p1._editDistance > p2._editDistance;
		}
	};

	//unordered_set<u_int32_t> _isDoneCorrection;

	void performErrorCorrection(u_int32_t nodeName, const DnaBitset* sequenceModel, VariantQueue& sequenceCopies, string& correctedSequence, const auto& alignment_engine, auto& graph){
		
		if(sequenceModel == nullptr){
			cout << "pas normal" << endl;
			correctedSequence = "";
			return;
		}

		if(sequenceCopies.size() == 0){
			char* seq = sequenceModel->to_string();
			correctedSequence = string(seq);
			free(seq);
			return;
		}

		
		//cout << (sequenceModel == nullptr) << endl;
		//cout << sequenceModel->to_string() << endl;
		//correctedSequence = string(sequenceModel->to_string());
		//return;

		//correctedSequence = string(sequenceModel->to_string(), sequenceModel->m_len);
		//cout << correctedSequence << endl;
		
		static EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0);

		/*
		char* sequenceModelStr = sequenceModel->to_string();

		//cout << "------" << endl;
		priority_queue<AlignScore, vector<AlignScore>, AlignScore_Comparator> queue;
		for(size_t i=0; i<sequenceCopies.size(); i++){
			//const string& sequence = sequenceCopies[i];
			const DnaBitset* dna = sequenceCopies[i];
			char* dnaStr = dna->to_string();

			EdlibAlignResult result = edlibAlign(sequenceModelStr, sequenceModel->m_len, dnaStr, dna->m_len, config);
			if (result.status == EDLIB_STATUS_OK) {
				queue.push({i, result.editDistance});
			}
	
			cout << endl << sequenceModelStr << endl;
			cout << dnaStr << endl;

			cout << result.editDistance << " " << result.alignmentLength << endl;
			//cout << result.alignment << endl;
			//printf("%s", result.alignment);
			//cout << string(result.alignment, result.alignmentLength) << endl;
			//char* cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
			//printf("%s", cigar);
			//cout << strlen(cigar) << endl;
			//free(cigar);
			edlibFreeAlignResult(result);

		}
		*/
		



		graph.Clear();
		//auto alignment = alignment_engine->Align(sequenceModel, graph);

		//char* sequenceModelStr = sequenceModel->to_string();
		//auto alignment = alignment_engine->Align(sequenceModelStr, sequenceModel->m_len, graph);
		//graph.AddAlignment(alignment, sequenceModelStr, sequenceModel->m_len);

		u_int32_t processedSequences = 0;

		/*
		size_t sketchSize = sequenceCopies.size();
		vector<KminmerSequenceVariant> variant_reversed(sketchSize);

		for(size_t i=0; i<sketchSize; i++){

			const KminmerSequenceVariant& variant = sequenceCopies.top();
			if(nodeName == 17326){
				cout << variant._editDistance << endl;
			}

			variant_reversed[i] = variant; //variant_reversed.size()-1-i
			sequenceCopies.pop();

		}
		*/

		/*
		vector<KminmerSequenceVariant> variant_reversed;
		while(!sequenceCopies.empty() && processedSequences < 5){
			const KminmerSequenceVariant& variant = sequenceCopies.top();
			
			if(nodeName == 17326){
				cout << variant._editDistance << endl;
			}

			sequenceCopies.pop();
			variant_reversed.push_back(variant);
		}

		std::reverse(variant_reversed.begin(), variant_reversed.end());
		*/

		vector<DnaBitset*> sequences;

		while(!sequenceCopies.empty()){
		//for(size_t i=0; i<variant_reversed.size() && processedSequences < 5; i++){
			//const KminmerSequenceVariant& variant = variant_reversed[i]; //sequenceCopies.top();
			const KminmerSequenceVariant& variant = sequenceCopies.top();
			
			//if(nodeName == 17326){
			//	cout << variant._editDistance << endl;
				//cout << variant._sequence << endl;
			//}

			
			//if(nodeName == 17326){
			//	cout << variant._sequence << endl;
			//}
			sequences.push_back(variant._sequence);
			
			sequenceCopies.pop();
			//processedSequences += 1;
		}

		//while(!sequenceCopies.empty() && processedSequences < 5){
		//for(size_t i=0; i<variant_reversed.size() && processedSequences < 5; i++){
			//const KminmerSequenceVariant& variant = variant_reversed[i]; //sequenceCopies.top();
			//const KminmerSequenceVariant& variant = sequenceCopies.top();
			
			//if(nodeName == 17326){
			//	cout << variant._editDistance << endl;
			//}
			//sequenceCopies.pop();
			

			//const AlignScore& s = queue.top();
			//queue.pop();

		//::reverse(sequences.begin(), sequences.end());
		std::reverse(sequences.begin(), sequences.end());

		for(size_t i=0; i<sequences.size(); i++){ //&& processedSequences < 5

			DnaBitset* dna = sequences[i];
			//const DnaBitset* dna = variant._sequence; //sequenceCopies[s._sequenceIndex];
			char* dnaStr = dna->to_string();

			//cout << s._sequenceIndex << " " << s._editDistance << endl;

			auto alignment = alignment_engine->Align(dnaStr, dna->m_len, graph);
			graph.AddAlignment(alignment, dnaStr, dna->m_len);

			free(dnaStr);

			processedSequences += 1;
		}


		//cout << endl;
		const vector<string>& msa = graph.GenerateMultipleSequenceAlignment();
		vector<vector<u_int32_t>> counts(msa[0].size(), vector<u_int32_t>(4, 0));

		for(const string& seq : msa){
			for(size_t i=0; i<seq.size(); i++){
				if(seq[i] == 'A'){
					counts[i][0] += 1;
				}
				else if(seq[i] == 'C'){
					counts[i][1] += 1;
				}
				else if(seq[i] == 'G'){
					counts[i][2] += 1;
				}
				else if(seq[i] == 'T'){
					counts[i][3] += 1;
				}
			}
		}

		float t = msa.size() * 0.45;

		correctedSequence.clear();

		for(size_t i=0; i<counts.size(); i++){
			for(size_t j=0; j<4; j++){
				if(counts[i][j] > t){

					if(j == 0){
						correctedSequence += 'A';
					}
					else if(j == 1){
						correctedSequence += 'C';
					}
					else if(j == 2){
						correctedSequence += 'G';
					}
					else if(j == 3){
						correctedSequence += 'T';
					}
					
					break;
				}
			}
		}

		/*
		if(nodeName == 17326){
			cout << "T: " << t << endl;
			for (const auto& it : msa) {
				std::cerr << it << std::endl;
			}
			for(size_t i=0; i<4; i++){
				for(size_t j=0; j<counts.size(); j++){
					cout << counts[j][i];
				}
				cout << endl;
			}
			cout << correctedSequence << endl;
		}
		*/

		//correctedSequence.clear();
		//correctedSequence = consensus;
		/*

		*/

		//correctedSequence = graph.GenerateConsensus();
		//cout << correctedSequence << endl;




		//if(sequenceCopies.size() > 60){
		//	correctedSequence = sequenceModel;
		//	return;
		//}

		//static EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0);
		/*
		for(const string& sequence : sequenceCopies){
			EdlibAlignResult result = edlibAlign(sequenceModel.c_str(), sequenceModel.size(), sequence.c_str(), sequence.size(), edlibDefaultAlignConfig());
			if (result.status == EDLIB_STATUS_OK) {
				//cout << sequenceCopies.size() << " " << result.editDistance << endl;
				//printf("edit_distance('hello', 'world!') = %d\n", result.editDistance);
			}
			edlibFreeAlignResult(result);

			//if(sequenceModel.size() > 50 && sequenceModel.size() < 70 && sequenceCopies.size() > 100){
			//	cout << sequenceModel << endl;
			//	cout << sequence << endl;
			//}
		}*/

		/*
		std::vector<std::string> sequences = {
		"CATAAAAGAACGTAGGTCGCCCGTCCGTAACCTGTCGGATCACCGGAAAGGACCCGTAAAGTGATAATGAT",
		"ATAAAGGCAGTCGCTCTGTAAGCTGTCGATTCACCGGAAAGATGGCGTTACCACGTAAAGTGATAATGATTAT",
		"ATCAAAGAACGTGTAGCCTGTCCGTAATCTAGCGCATTTCACACGAGACCCGCGTAATGGG",
		"CGTAAATAGGTAATGATTATCATTACATATCACAACTAGGGCCGTATTAATCATGATATCATCA",
		"GTCGCTAGAGGCATCGTGAGTCGCTTCCGTACCGCAAGGATGACGAGTCACTTAAAGTGATAAT",
		"CCGTAACCTTCATCGGATCACCGGAAAGGACCCGTAAATAGACCTGATTATCATCTACAT"
		};




		*/


		/*
		std::cerr << ">Consensus LN:i:" << consensus.size() << std::endl << consensus << std::endl;

		auto msa = graph.GenerateMultipleSequenceAlignment();

		for (const auto& it : msa) {
			std::cerr << it << std::endl;
		}
		*/


		//if(sequenceModel.size() > 50 && sequenceModel.size() < 70  && sequenceCopies.size() > 100){
		//	exit(1);
		//}

		/*
		if(sequenceCopies.size() > 60){
			correctedSequence = sequenceModel;
			return;
		}

		//correctedSequence = sequenceModel;
		//return;
		
		seqan::Align<seqan::String<seqan::Dna> > align;
		seqan::resize(rows(align), 10);
		for (int i = 0; i < 10; ++i){
			seqan::assignSource(row(align, i), sequenceCopies[i].c_str());
		}


		seqan::globalMsaAlignment(align, seqan::Blosum62(-1, -11));
    	//cout << align << "\n";

		seqan::String<seqan::ProfileChar<seqan::Dna> > profile;
		seqan::resize(profile, length(row(align, 0)));
		for (unsigned rowNo = 0; rowNo < 4u; ++rowNo)
			for (unsigned i = 0; i < length(row(align, rowNo)); ++i)
				profile[i].count[ordValue(getValue(row(align, rowNo), i))] += 1;

		// call consensus from this string
		seqan::DnaString consensus;
		for (unsigned i = 0; i < length(profile); ++i)
		{
			int idx = getMaxIndex(profile[i]);
			if (idx < 4)  // is not gap
				appendValue(consensus, seqan::Dna(getMaxIndex(profile[i])));
		}

		correctedSequence.clear();
		correctedSequence.assign(begin(consensus), end(consensus));

		cout << correctedSequence << endl;

		//const char* lala = consensus.cstr();
		//std::basic_string lala(consensus.cstr());
		//string lala = string());
		//CString cs("Hello");
		//std::string s((LPCTSTR) toCString(consensus));

		
		//exit(1);
		*/

	}

	/*
	void extractKminmerSequences_all (){

		cout << "Extracting kminmer sequence model" << endl;

		auto fp = std::bind(&ToBasespace::extractKminmerSequences_all_read, this, std::placeholders::_1, std::placeholders::_2);
		ReadParser readParser(_inputFilename, false);
		readParser.parse(fp);
	}

	void extractKminmerSequences_all_read(kseq_t* read, u_int64_t readIndex){
		

		//cout << readIndex << endl;
		//ottalSize += strlen(read->seq.s);
					
		string kminmerSequence;
		char* sequenceOriginal = read->seq.s;

		string rleSequence;
		vector<u_int64_t> rlePositions;
		Encoder::encode_rle(read->seq.s, strlen(read->seq.s), rleSequence, rlePositions);

		vector<u_int64_t> minimizers;
		vector<u_int64_t> minimizers_pos;
		_minimizerParser->parse(rleSequence, minimizers, minimizers_pos);

		vector<KmerVec> kminmers; 
		vector<ReadKminmer> kminmersInfo;
		MDBG::getKminmers(_minimizerSize, _kminmerSize, minimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, readIndex, false);

		for(size_t i=0; i<kminmers.size(); i++){
			if(_mdbg->_dbg_nodes.find(kminmers[i]) == _mdbg->_dbg_nodes.end()) continue;

			u_int32_t nodeName = _mdbg->_dbg_nodes[kminmers[i]]._index;
			ReadKminmer& kminmerInfo = kminmersInfo[i];

			if(_nodeName_all_left.find(nodeName) == _nodeName_all_left.end()){
				extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Left, kminmerSequence);
				_nodeName_all_left[nodeName] = new DnaBitset(kminmerSequence);
			}	
			if(_nodeName_all_right.find(nodeName) == _nodeName_all_right.end()){
				extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Right, kminmerSequence);
				_nodeName_all_right[nodeName] = new DnaBitset(kminmerSequence);
			}
		}

	}
	*/

	void extractKminmerSequences_allVariants (){

		cout << "Extracting kminmer sequence variants" << endl;

		auto fp = std::bind(&ToBasespace::extractKminmerSequences_allVariants_read, this, std::placeholders::_1, std::placeholders::_2);
		ReadParser readParser(_inputFilename, false, !_isFirstPass);
		readParser.parse(fp);

	}

	unordered_set<u_int32_t> _isDoneNodeName_left;
	unordered_set<u_int32_t> _isDoneNodeName_right;

	void extractKminmerSequences_allVariants_read(kseq_t* read, u_int64_t readIndex){
		//ottalSize += strlen(read->seq.s);
					
		string kminmerSequence;
		char* sequenceOriginal = read->seq.s;

		string rleSequence;
		vector<u_int64_t> rlePositions;
		Encoder::encode_rle(read->seq.s, strlen(read->seq.s), rleSequence, rlePositions);

		vector<u_int64_t> minimizers;
		vector<u_int64_t> minimizers_pos;
		_minimizerParser->parse(rleSequence, minimizers, minimizers_pos);

		vector<KmerVec> kminmers; 
		vector<ReadKminmer> kminmersInfo;
		MDBG::getKminmers(_minimizerSize, _kminmerSize, minimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, readIndex, false);

		for(size_t i=0; i<kminmers.size(); i++){
			if(_mdbg->_dbg_nodes.find(kminmers[i]) == _mdbg->_dbg_nodes.end()) continue;

			u_int32_t nodeName = _mdbg->_dbg_nodes[kminmers[i]]._index;
			ReadKminmer& kminmerInfo = kminmersInfo[i];

			//cout << nodeName << endl;
			if(_isDoneNodeName_left.find(nodeName) == _isDoneNodeName_left.end()){
				//DnaBitset* sequenceModel = _nodeName_all_left[nodeName];
				extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Left, kminmerSequence);
				addKminmerSequenceVariant_add_all(nodeName, _kminmerSequenceCopies_all_left, kminmerSequence, _outputFile_left, _isDoneNodeName_left);
			}
			if(_isDoneNodeName_right.find(nodeName) == _isDoneNodeName_right.end()){
				//DnaBitset* sequenceModel = _nodeName_all_right[nodeName];
				extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Right, kminmerSequence);
				addKminmerSequenceVariant_add_all(nodeName, _kminmerSequenceCopies_all_right, kminmerSequence, _outputFile_right, _isDoneNodeName_right);
			}

			/*
			if(_requiredCopiers_entire.find(nodeName) != _requiredCopiers_entire.end()){
				extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Entire, kminmerSequence);
				addKminmerSequenceVariant_add(nodeName, seq._readIndex, variants, seq._sequence, sequence);
				//addKminmerSequenceVariant(nodeName, readIndex, _nodeName_entire, _nodeName_entire_multi, _kminmerSequenceCopies_entire, kminmerSequence);
				//_kminmerSequenceCopies_entire[nodeName].push_back(new DnaBitset(kminmerSequence));
			}
			if(_requiredCopiers_left.find(nodeName) != _requiredCopiers_left.end()){//} && _kminmerSequenceCopies_left[nodeName].size() < 50){
				extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Left, kminmerSequence);
				addKminmerSequenceVariant(nodeName, readIndex, _nodeName_left, _nodeName_left_multi, _kminmerSequenceCopies_left, kminmerSequence);
				//_kminmerSequenceCopies_left[nodeName].push_back(new DnaBitset(kminmerSequence));
			}
			if(_requiredCopiers_right.find(nodeName) != _requiredCopiers_right.end()){// && _kminmerSequenceCopies_right[nodeName].size() < 50){
				extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Right, kminmerSequence);
				addKminmerSequenceVariant(nodeName, readIndex, _nodeName_right, _nodeName_right_multi, _kminmerSequenceCopies_right, kminmerSequence);
				//_kminmerSequenceCopies_right[nodeName].push_back(new DnaBitset(kminmerSequence));
			}
			*/

				
		}
	}


	void addKminmerSequenceVariant_add_all(u_int32_t nodeName, auto& variants, const string& sequence, gzFile& file, auto& isDoneNodeName){
		//if(sequenceModel == nullptr){
		//	cout << "pas normal" << endl;
		//	return; //model not found yet
		//}

		//if(_isDoneCorrection.find(nodeName) != _isDoneCorrection.end()) return;

		vector<DnaBitset*>& queue = variants[nodeName];
		queue.push_back(new DnaBitset(sequence));

		if(queue.size() < 20) return;

		//cout << "correcting" << endl;
		//DnaBitset* sequenceModel = models[nodeName];

		string correctedSequence;
		performErrorCorrection_all(nodeName, queue, correctedSequence);

		writeKminmerSequence_all(nodeName, correctedSequence, file);

		for(DnaBitset* dna : queue){
			delete dna;
		}
		
		variants.erase(nodeName);

		isDoneNodeName.insert(nodeName);
		//variants.erase(nodeName);

		//if(queue.size() < 20){
		//	queue.push({0, new DnaBitset(sequence)});
		//}

		/*
		return;
		
		static EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0);


		char* sequenceModelStr = sequenceModel->to_string();

		EdlibAlignResult result = edlibAlign(sequenceModelStr, sequenceModel->m_len, sequence.c_str(), sequence.size(), config);
		free(sequenceModelStr);

		if (result.status != EDLIB_STATUS_OK){
			edlibFreeAlignResult(result);
			return;
		}
		
		
		if(queue.size() < 20){
			queue.push({result.editDistance, new DnaBitset(sequence)});
		}
		else{
			if(result.editDistance < queue.top()._editDistance){

				const KminmerSequenceVariant& variant = queue.top();
				delete variant._sequence;

				queue.pop();
				queue.push({result.editDistance, new DnaBitset(sequence)});
			}
		}
		
		edlibFreeAlignResult(result);
		*/
	}

	void performErrorCorrection_all(u_int32_t nodeName, const vector<DnaBitset*>& sequences, string& correctedSequence){
		
		
		if(sequences.size() == 0){
			cout << "pas normal" << endl;
			correctedSequence = "";
			return;
		}

		/*
		if(sequences.size() == 1){
			char* seq = sequences[0]->to_string();
			correctedSequence = string(seq);
			free(seq);
			return;
		}*/



		
		//static EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0);
		static auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -3);  // linear gaps





		_graph.Clear();

		/*
		u_int32_t processedSequences = 0;

		vector<DnaBitset*> sequences;

		while(!sequenceCopies.empty()){
		//for(size_t i=0; i<variant_reversed.size() && processedSequences < 5; i++){
			//const KminmerSequenceVariant& variant = variant_reversed[i]; //sequenceCopies.top();
			const KminmerSequenceVariant& variant = sequenceCopies.top();
			
			//if(nodeName == 17326){
			//	cout << variant._editDistance << endl;
				//cout << variant._sequence << endl;
			//}

			
			//if(nodeName == 17326){
			//	cout << variant._sequence << endl;
			//}
			sequences.push_back(variant._sequence);
			
			sequenceCopies.pop();
			//processedSequences += 1;
		}
		*/


		//std::reverse(sequences.begin(), sequences.end());

		for(size_t i=0; i<sequences.size(); i++){ //&& processedSequences < 5

			DnaBitset* dna = sequences[i];
			//const DnaBitset* dna = variant._sequence; //sequenceCopies[s._sequenceIndex];
			char* dnaStr = dna->to_string();

			//cout << s._sequenceIndex << " " << s._editDistance << endl;

			auto alignment = alignment_engine->Align(dnaStr, dna->m_len, _graph);
			_graph.AddAlignment(alignment, dnaStr, dna->m_len);

			free(dnaStr);

			//processedSequences += 1;
		}


		//cout << endl;
		const vector<string>& msa = _graph.GenerateMultipleSequenceAlignment();
		vector<vector<u_int32_t>> counts(msa[0].size(), vector<u_int32_t>(4, 0));

		for(const string& seq : msa){
			for(size_t i=0; i<seq.size(); i++){
				if(seq[i] == 'A'){
					counts[i][0] += 1;
				}
				else if(seq[i] == 'C'){
					counts[i][1] += 1;
				}
				else if(seq[i] == 'G'){
					counts[i][2] += 1;
				}
				else if(seq[i] == 'T'){
					counts[i][3] += 1;
				}
			}
		}

		float t = msa.size() * 0.45;

		correctedSequence.clear();

		for(size_t i=0; i<counts.size(); i++){
			for(size_t j=0; j<4; j++){
				if(counts[i][j] > t){

					if(j == 0){
						correctedSequence += 'A';
					}
					else if(j == 1){
						correctedSequence += 'C';
					}
					else if(j == 2){
						correctedSequence += 'G';
					}
					else if(j == 3){
						correctedSequence += 'T';
					}
					
					break;
				}
			}
		}


	}


	void writeKminmerSequence_all(u_int32_t nodeName, const string& sequence, const gzFile& file){
		
		//char* dnaStr = dna->to_string();

		//const string& sequence = 
		u_int16_t length = sequence.size();

		//cout << "Writing: " << endl;
		//cout << nodeName << " " << length << " " << sequence << endl;
		gzwrite(file, (const char*)&nodeName, sizeof(nodeName));
		gzwrite(file, (const char*)&length, sizeof(length));
		gzwrite(file, (const char*)&sequence[0], length);
	}

	/*
	void writeKminmerSequence(u_int32_t nodeName, u_int64_t readIndex, const string& sequence, const gzFile& file){
		
		u_int16_t length = sequence.size();

		//cout << "Writing: " << endl;
		//cout << nodeName << " " << length << " " << sequence << endl;
		gzwrite(file, (const char*)&nodeName, sizeof(nodeName));
		gzwrite(file, (const char*)&readIndex, sizeof(readIndex));
		gzwrite(file, (const char*)&length, sizeof(length));
		gzwrite(file, (const char*)&sequence[0], length);
	}

	
	void loadKminmerSequences(){
		loadKminmerSequences_aux(_filename_kminmerSequences + "_entire.gz", _kminmerSequences_entire);
		loadKminmerSequences_aux(_filename_kminmerSequences + "_left.gz", _kminmerSequences_left);
		loadKminmerSequences_aux(_filename_kminmerSequences + "_right.gz", _kminmerSequences_right);
	}*/

};	


#endif 


