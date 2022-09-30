

#ifndef MDBG_METAG_TOBASESPACENOCORRECTION
#define MDBG_METAG_TOBASESPACENOCORRECTION

#include "../Commons.hpp"
//#include "../utils/edlib.h"
//#include "../utils/spoa/include/spoa/spoa.hpp"
#include "../utils/DnaBitset.hpp"
//#include <seqan/align.h>
//#include <seqan/graph_msa.h>
//#include <cstring>

struct KminmerSequence{
	u_int32_t _readIndex;
	DnaBitset* _sequence;
};

class ToBasespaceNoCorrection : public Tool{
    
public:

	string _inputFilename;
	string _inputFilenameContig;
	//string _inputFilenameContig_fasta;
	string _inputDir;
	float _minimizerDensity;
    size_t _minimizerSize;
    size_t _kminmerSize;
	bool _isFirstPass;
	bool _isOutputFasta;

	string _filename_outputContigs;
	string _filename_kminmerSequences;
	MDBG* _mdbg;
	MinimizerParser* _minimizerParser;
	
	unordered_map<u_int32_t, KminmerSequence> _nodeName_entire;
	unordered_map<u_int32_t, KminmerSequence> _nodeName_right;
	unordered_map<u_int32_t, KminmerSequence> _nodeName_left;

	enum LoadType{
        Entire,
        EntireRc,
        Left,
        Right,
        LeftLast,
        RightLast
	};

	ToBasespaceNoCorrection(): Tool (){


	}

	void parseArgs(int argc, char* argv[]){


		cxxopts::Options options("ToBasespace", "");
		options.add_options()
		(ARG_INPUT_FILENAME, "", cxxopts::value<string>())
		(ARG_INPUT_FILENAME_CONTIG, "", cxxopts::value<string>())
		(ARG_OUTPUT_FILENAME, "", cxxopts::value<string>())
		//(ARG_INPUT_FILENAME_CONTIG_FASTA, "", cxxopts::value<string>()->default_value(""))
		(ARG_FIRST_PASS, "", cxxopts::value<bool>()->default_value("false"))
		(ARG_FASTA, "", cxxopts::value<bool>()->default_value("false"))
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
			_inputFilenameContig = result[ARG_INPUT_FILENAME_CONTIG].as<string>();
			//_inputFilenameContig_fasta = result[ARG_INPUT_FILENAME_CONTIG_FASTA].as<string>();
			_isFirstPass = result[ARG_FIRST_PASS].as<bool>();
			_isOutputFasta = result[ARG_FASTA].as<bool>();
			_filename_outputContigs = result[ARG_OUTPUT_FILENAME].as<string>();
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

		_kminmerSize = 4;

		cout << endl;
		cout << "Input dir: " << _inputDir << endl;
		//cout << "Output filename: " << _outputFilename << endl;
		cout << "Minimizer length: " << _minimizerSize << endl;
		cout << "Kminmer length: " << _kminmerSize << endl;
		cout << "Density: " << _minimizerDensity << endl;
		cout << endl;

		//_filename_outputContigs = _inputFilenameContig + ".fasta.gz"; //_inputDir + "/tmpContigs.fasta.gz";
		_minimizerParser = new MinimizerParser(_minimizerSize, _minimizerDensity);
	}


    void execute (){
		
		cout << "Loading mdbg" << endl;
		string mdbg_filename = _inputDir + "/mdbg_nodes_init.gz";
		_mdbg = new MDBG(_kminmerSize);
		_mdbg->load(mdbg_filename);
		cout << "MDBG nodes: " << _mdbg->_dbg_nodes.size() << endl;

		//cout << _inputFilenameContig_fasta << endl;
		loadContigs_min(_inputFilenameContig);

		extractKminmerSequences();
		//lalalala();

		//exit(1);
		//loadKminmerSequences();
		createBaseContigs(_inputFilenameContig, _filename_outputContigs.c_str());
		//createBaseContigs(_inputDir + "/minimizer_contigs_complete.gz", _filename_outputContigs.c_str());
		//createBaseContigs(_inputDir + "/eval/composition//22//debug_longUnitigs.gz", _inputDir + "/eval/composition//22//debug_longUnitigs.fasta.gz");
		//createBaseContigs(_inputDir + "/eval/composition//3//debug_longUnitigsNeighbors.gz", _inputDir + "/eval/composition//3/debug_longUnitigsNeighbors.fasta.gz");

		cout << endl << "Contig filename: " << _filename_outputContigs << endl;
		delete _mdbg;
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




	void loadContigs_min(const string& contigFilename){

		cout << "Extracting kminmers (contigs)" << endl;
		KminmerParser parser(contigFilename, _minimizerSize, _kminmerSize, false);
		auto fp = std::bind(&ToBasespaceNoCorrection::loadContigs_min_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
		parser.parseMinspace(fp);

		/*
		cout << "Loading mContigs: " << contigFilename << endl;
		gzFile contigFile = gzopen(contigFilename.c_str(),"rb");
		u_int64_t nbContigs = 0;
		
		while(true){

			vector<u_int32_t> minimizers;
			//vector<u_int64_t> supportingReads;
			u_int64_t size;
			gzread(contigFile, (char*)&size, sizeof(size));
			

			if(gzeof(contigFile)) break;

			minimizers.resize(size);
			//supportingReads.resize(size);
			gzread(contigFile, (char*)&minimizers[0], size * sizeof(u_int32_t));
			//gzread(contigFile, (char*)&supportingReads[0], size * sizeof(u_int64_t));
			


			
			//for(u_int32_t nodeIndex : nodePath){
			for(size_t i=0; i<minimizers.size(); i++){
				u_int32_t nodeIndex = nodePath[i];
				bool orientation;
				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex, orientation);
				//u_int32_t readIndex = supportingReads[i];

				//if(nbContigs == 0) cout << nodeName << endl;
				//if(nodeName == 21594){
				//	cout << "wesh1 ? " << endl;
				//	getchar();
				//}

				//ContigNode contigNode = {nodeName, supportingReads[i]};

				//cout << nodeName << " " << supportingReads[i] << endl;

				//cout << nodeName << endl;
				//_usedNodeNames.insert(nodeName);


				if(i == 0){
					initKminmerSequence(nodeName, readIndex, _nodeName_entire);

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
						initKminmerSequence(nodeName, readIndex, _nodeName_right);
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
						initKminmerSequence(nodeName, readIndex, _nodeName_left);
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
		*/
	}

	void loadContigs_min_read(const vector<u_int64_t>& readMinimizers, const vector<ReadKminmerComplete>& kminmersInfos, u_int64_t readIndex){

		//cout << readIndex << " " << kminmersInfos.size() << endl;
		for(size_t i=0; i<kminmersInfos.size(); i++){
			
			const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];

			KmerVec vec = kminmerInfo._vec;
			
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()){
				cout << "Not found" << endl;
				continue;
			}

			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
			
			//vector<u_int64_t> minimizerSeq;
			
			//for(size_t i=kminmerInfo._read_pos_start; i<=kminmerInfo._read_pos_end; i++){
			//	minimizerSeq.push_back(readMinimizers[i]);
			//}
			

			//if(kminmerInfo._isReversed){
			//	std::reverse(minimizerSeq.begin(), minimizerSeq.end());
			//}

			bool orientation = !kminmerInfo._isReversed;

			if(i == 0){
				initKminmerSequence(nodeName, readIndex, _nodeName_entire);
			}
			else {
				if(orientation){ //+
					initKminmerSequence(nodeName, readIndex, _nodeName_right);
				}
				else{ //-
					initKminmerSequence(nodeName, readIndex, _nodeName_left);
				}
			}

		}

	}

	void initKminmerSequence(u_int32_t nodeName, u_int32_t readIndex, auto& dictSimple){


		dictSimple[nodeName] = {readIndex, nullptr};
		/*
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
		*/
	}




	void extractKminmerSequences (){

		cout << "Extracting kminmer sequences" << endl;
		//auto fp = std::bind(&ToBasespaceNoCorrection::extractKminmerSequences_read, this, std::placeholders::_1, std::placeholders::_2);
		//ReadParser readParser(_inputFilename, false, !_isFirstPass);
		//readParser.parse(fp);
		string kminmerSequence;
		ifstream kminmerFile(_inputDir + "/kminmerData.txt");

		while (true) {

			u_int32_t sizeData = -1;
			kminmerFile.read((char*)&sizeData, sizeof(sizeData));

			if(kminmerFile.eof())break;

			u_int32_t sizeSequence = -1;
			kminmerFile.read((char*)&sizeSequence, sizeof(sizeSequence));

			uint8_t* m_data = new uint8_t[sizeData];
			kminmerFile.read((char*)&m_data[0], sizeData*sizeof(uint8_t));

			DnaBitset* dnaBitset = new DnaBitset(m_data, sizeData, sizeSequence);
			char* kminmerSequenceComplete = dnaBitset->to_string();

			u_int32_t nodeName = -1;;
			u_int32_t lengthStart = -1;
			u_int32_t lengthEnd = -1;
			//bool isReversed = false;

			kminmerFile.read((char*)&nodeName, sizeof(nodeName));
			kminmerFile.read((char*)&lengthStart, sizeof(lengthStart));
			kminmerFile.read((char*)&lengthEnd, sizeof(lengthEnd));

			//kminmerFile.read((char*)&isReversed, sizeof(isReversed));

			//cout << nodeName << " " << lengthStart << " " << lengthEnd << " " << isReversed << " " << strlen(kminmerSequenceComplete) << endl;
			//cout << kminmerSequenceComplete << endl;
			//getchar();
			//u_int32_t nodeName = _mdbg->_dbg_nodes[kminmers[i]]._index;
			//ReadKminmer& kminmerInfo = kminmersInfo[i];

			//cout << i << endl;
			//ContigNode contigNode = {nodeName, readIndex};
			//void extractKminmerSequence(const char* kminmerSequenceComplete, bool isReversed, u_int32_t lengthStart, u_int32_t lengthEnd, u_int32_tLoadType loadType, string& sequence){



			if(_nodeName_entire.find(nodeName) != _nodeName_entire.end() && _nodeName_entire[nodeName]._sequence == nullptr){
				//cout << "entire" << endl;
				extractKminmerSequence(kminmerSequenceComplete, lengthStart, lengthEnd, LoadType::Entire, kminmerSequence);
				_nodeName_entire[nodeName] = {0, new DnaBitset(kminmerSequence)};
			}
			

			if(_nodeName_left.find(nodeName) != _nodeName_left.end() && _nodeName_left[nodeName]._sequence == nullptr){
				//cout << "left" << endl;
				extractKminmerSequence(kminmerSequenceComplete, lengthStart, lengthEnd, LoadType::Left, kminmerSequence);
				_nodeName_left[nodeName] = {0, new DnaBitset(kminmerSequence)};

			}
			

			if(_nodeName_right.find(nodeName) != _nodeName_right.end() && _nodeName_right[nodeName]._sequence == nullptr){
				//cout << "right" << endl;
				extractKminmerSequence(kminmerSequenceComplete, lengthStart, lengthEnd, LoadType::Right, kminmerSequence);
				_nodeName_right[nodeName] = {0, new DnaBitset(kminmerSequence)};
			}

			//if(nodeName == 2108){
			//	cout << "lala" << endl;
			//	getchar();
			//}

			free(kminmerSequenceComplete);
			delete dnaBitset;

		}

		kminmerFile.close();
	}

	/*
			if(_mdbg->_dbg_nodes[vec]._abundance == 2){
				//DnaBitset* contigSequenceBitset = new DnaBitset(contigSequence);
				u_int32_t sizeData = sequenceBitset->_bitsetSize;
				u_int32_t sizeSeq = sequenceBitset->m_len;
				uint8_t* m_data = sequenceBitset->m_data;
				_kminmerFile.write((const char*)&sizeData, sizeof(sizeData));
				_kminmerFile.write((const char*)&sizeSeq, sizeof(sizeSeq));
				_kminmerFile.write((const char*)&m_data[0], sizeData*sizeof(uint8_t));
				delete sequenceBitset;

				u_int32_t lengthStart = kminmerInfo._seq_length_start;
				u_int32_t lengthEnd = kminmerInfo._seq_length_end;
				bool isReversed = kminmerInfo._isReversed;

				_kminmerFile.write((const char*)&lengthStart, sizeof(lengthStart));
				_kminmerFile.write((const char*)&lengthEnd, sizeof(lengthEnd));
				_kminmerFile.write((const char*)&isReversed, sizeof(isReversed));
			}
	*/

	/*
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

		//for(size_t i=0; i<minimizers.size(); i++){
		//	cout << i << ": " << minimizers[i] << " " << minimizers_pos[i] << endl;
		//}

		vector<KmerVec> kminmers; 
		vector<ReadKminmer> kminmersInfo;
		MDBG::getKminmers(_minimizerSize, _kminmerSize, minimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, readIndex, false);



		for(size_t i=0; i<kminmers.size(); i++){
			if(_mdbg->_dbg_nodes.find(kminmers[i]) == _mdbg->_dbg_nodes.end()) continue;

			u_int32_t nodeName = _mdbg->_dbg_nodes[kminmers[i]]._index;
			ReadKminmer& kminmerInfo = kminmersInfo[i];

			//cout << i << endl;
			//ContigNode contigNode = {nodeName, readIndex};
	
			if(_nodeName_entire.find(nodeName) != _nodeName_entire.end() && _nodeName_entire[nodeName]._sequence == nullptr){
				extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Entire, kminmerSequence);
				_nodeName_entire[nodeName] = {readIndex, new DnaBitset(kminmerSequence)};
			}
			

			if(_nodeName_left.find(nodeName) != _nodeName_left.end() && _nodeName_left[nodeName]._sequence == nullptr){
				
				extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Left, kminmerSequence);
				_nodeName_left[nodeName] = {readIndex, new DnaBitset(kminmerSequence)};

			}
			

			if(_nodeName_right.find(nodeName) != _nodeName_right.end() && _nodeName_right[nodeName]._sequence == nullptr){
				extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Right, kminmerSequence);
				_nodeName_right[nodeName] = {readIndex, new DnaBitset(kminmerSequence)};
			}



				
		}

		//cout << readIndex << " " << _nodeName_left.size() << " " << _nodeName_right.size() << endl;
	}
	*/


	void extractKminmerSequence(const char* kminmerSequenceComplete, u_int32_t lengthStart, u_int32_t lengthEnd, LoadType loadType, string& sequence){

		u_int32_t startPosition = 0;
		u_int32_t len = 0;

		/*
		//cout << "lala    " << string(sequenceOriginal) << endl;


		startPosition = kminmerInfo._read_pos_start;
		len = kminmerInfo._read_pos_end - kminmerInfo._read_pos_start;

		char subbuff[len+1];
		memcpy( subbuff, &sequenceOriginal[startPosition], len);
		subbuff[len] = '\0';
		sequence = string(subbuff);
		*/

		sequence = string(kminmerSequenceComplete);

		//if(isReversed){
		//	Utils::revcomp(sequence);
		//}


		if(loadType == LoadType::Entire){
			return;
			//startPosition = kminmerInfo._read_pos_start;
			//len = kminmerInfo._read_pos_end - kminmerInfo._read_pos_start;
		}
		else if(loadType == LoadType::Left){
			startPosition = 0; //kminmerInfo._read_pos_start;
			len = lengthStart; //kminmerInfo._position_of_second_minimizer_seq - kminmerInfo._read_pos_start;
			//cout << kminmerInfo._read_pos_start << " " << kminmerInfo._position_of_second_minimizer << endl;
		}
		else if(loadType == LoadType::Right){
			//return;
			startPosition = sequence.size() - lengthEnd;  //kminmerInfo._position_of_second_to_last_minimizer_seq;
			len = lengthEnd; //kminmerInfo._read_pos_end - kminmerInfo._position_of_second_to_last_minimizer_seq;
			//return;
			//cout << kminmerInfo._read_pos_end << " " << kminmerInfo._position_of_second_to_last_minimizer << endl;
		}
		
		/*
		if(startPosition > 99999999){
			//cout << readIndex << endl;
			//cout << sequenceOriginal << endl << endl;
			cout << sequence << endl;
			cout << sequence.size() << endl;
			cout << kminmerInfo._read_pos_start << endl;
			cout << kminmerInfo._read_pos_end << endl;
			cout << kminmerInfo._seq_length_start << endl;
			cout << kminmerInfo._seq_length_end << endl;
			cout << len << endl;
			cout << kminmerInfo._isReversed << endl;
			cout << loadType << endl;
			cout << startPosition << endl;
			getchar();
		}
		*/


		sequence = sequence.substr(startPosition, len);


	}

	/*
	void extractKminmerSequence(const char* sequenceOriginal, const ReadKminmer& kminmerInfo, LoadType loadType, string& sequence){

		//cout << "lala    " << string(sequenceOriginal) << endl;
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
		



		sequence = sequence.substr(startPosition, len);


	}
	*/

	void createBaseContigs(const string& contigFilename, const string& outputFilename){

		cout << "Creating basespace contigs: " << contigFilename << " " << outputFilename << endl;

		if(_isOutputFasta){
			_basespaceContigFile = gzopen(outputFilename.c_str(),"wb");
		}
		else{
			_contigFile_bitset = ofstream(contigFilename + ".bitset");
		}

		KminmerParser parser(contigFilename, _minimizerSize, _kminmerSize, false);
		auto fp = std::bind(&ToBasespaceNoCorrection::createBaseContigs_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
		parser.parseMinspace(fp);

		if(_isOutputFasta){
			gzclose(_basespaceContigFile);
		}
		else{
			_contigFile_bitset.close();
		}
		

	}

	gzFile _basespaceContigFile;
	ofstream _contigFile_bitset;

	void createBaseContigs_read(const vector<u_int64_t>& readMinimizers, const vector<ReadKminmerComplete>& kminmersInfos, u_int64_t readIndex){


		cout << readIndex << " " << kminmersInfos.size() << endl;
		
		string contigSequence = "";

		for(size_t i=0; i<kminmersInfos.size(); i++){
			
			//cout << i << endl;
			const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];

			KmerVec vec = kminmerInfo._vec;
			
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()){
				cout << "Not found" << endl;
				continue;
			}

			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
			
			//vector<u_int64_t> minimizerSeq;
			
			//for(size_t i=kminmerInfo._read_pos_start; i<=kminmerInfo._read_pos_end; i++){
			//	minimizerSeq.push_back(readMinimizers[i]);
			//}
			

			//if(kminmerInfo._isReversed){
			//	std::reverse(minimizerSeq.begin(), minimizerSeq.end());
			//}

			bool orientation = !kminmerInfo._isReversed;

			if(i == 0){
				//cout << nodeName << " " << orientation << endl;
				if(orientation){ //+
					//cout << "Entire" << endl;

					//string seq = _nodeName_entire[contigNode];
					//string correctedSequence;
					if(_nodeName_entire.find(nodeName) == _nodeName_entire.end() || _nodeName_entire[nodeName]._sequence == nullptr){
						cout << "not found entire " << nodeName << endl;
					}
					
					char* seq = _nodeName_entire[nodeName]._sequence->to_string();
					string kminmerSequence = string(seq);
					free(seq);


					//performErrorCorrection(nodeName, getKminmerSequence(nodeName, readIndex, _nodeName_entire, _nodeName_entire_multi), _kminmerSequenceCopies_entire[contigNode], correctedSequence, alignment_engine, graph);
					

					contigSequence += kminmerSequence;
					//cout << contigSequence << endl;
				}
				else{
					//cout << "Entire RC" << endl;
					//string seq = ;
					//string correctedSequence;
					//performErrorCorrection(nodeName, getKminmerSequence(nodeName, readIndex, _nodeName_entire, _nodeName_entire_multi), _kminmerSequenceCopies_entire[contigNode], correctedSequence, alignment_engine, graph);
					
					if(_nodeName_entire.find(nodeName) == _nodeName_entire.end() || _nodeName_entire[nodeName]._sequence == nullptr){
						cout << "not found entire RC " << nodeName << endl;
					}

					char* seq = _nodeName_entire[nodeName]._sequence->to_string();
					string kminmerSequence = string(seq);
					free(seq);

					Utils::revcomp(kminmerSequence);
					contigSequence += kminmerSequence;
					//cout << contigSequence << endl;
				}
			}
			else {
				if(orientation){
					
					if(_nodeName_right.find(nodeName) == _nodeName_right.end() || _nodeName_right[nodeName]._sequence == nullptr){
						cout << "not found right " << nodeName << endl;
					}
					//cout << (_nodeName_right.find(nodeName) != _nodeName_right.end()) << endl;
					//cout << (_nodeName_right[nodeName]._sequence == nullptr) << endl;
					//cout << nodeName << endl;

					char* seq = _nodeName_right[nodeName]._sequence->to_string();
					string kminmerSequence = string(seq);
					free(seq);

					//string kminmerSequence = _nodeName_right[nodeName];
					contigSequence += kminmerSequence;
					

				}
				else{
					
					if(_nodeName_left.find(nodeName) == _nodeName_left.end() || _nodeName_left[nodeName]._sequence == nullptr){
						cout << "not found left " << nodeName << endl;
					}
					//cout << (_nodeName_left.find(nodeName) != _nodeName_left.end()) << endl;
					//cout << (_nodeName_left[nodeName]._sequence == nullptr) << endl;
					//cout << nodeName << endl;

					char* seq = _nodeName_left[nodeName]._sequence->to_string();
					string kminmerSequence = string(seq);
					free(seq);

					//string correctedSequence;
					//performErrorCorrection(nodeName, getKminmerSequence(nodeName, readIndex, _nodeName_left, _nodeName_left_multi), _kminmerSequenceCopies_left[contigNode], correctedSequence, alignment_engine, graph);
					
					//string kminmerSequence = _nodeName_left[nodeName];

					Utils::revcomp(kminmerSequence);
					contigSequence += kminmerSequence;


				}
			}


		}

		
		if(_isOutputFasta){
			string header = ">ctg" + to_string(readIndex) + '\n';
			gzwrite(_basespaceContigFile, (const char*)&header[0], header.size());
			contigSequence +=  '\n';
			gzwrite(_basespaceContigFile, (const char*)&contigSequence[0], contigSequence.size());
		}
		else{
			DnaBitset* contigSequenceBitset = new DnaBitset(contigSequence);
			u_int32_t sizeData = contigSequenceBitset->_bitsetSize;
			u_int32_t sizeSeq = contigSequenceBitset->m_len;
			uint8_t* m_data = contigSequenceBitset->m_data;
			_contigFile_bitset.write((const char*)&sizeData, sizeof(sizeData));
			_contigFile_bitset.write((const char*)&sizeSeq, sizeof(sizeSeq));
			_contigFile_bitset.write((const char*)&m_data[0], sizeData*sizeof(uint8_t));
			delete contigSequenceBitset;
		}


	}


};	


#endif 



