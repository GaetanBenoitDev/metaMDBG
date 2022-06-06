

#ifndef MDBG_METAG_TOBASESPACENOCORRECTION
#define MDBG_METAG_TOBASESPACENOCORRECTION

#include "../Commons.hpp"
//#include "../utils/edlib.h"
//#include "../utils/spoa/include/spoa/spoa.hpp"
#include "../utils/DnaBitset.hpp"
//#include <seqan/align.h>
//#include <seqan/graph_msa.h>
//#include <cstring>

class ToBasespaceNoCorrection : public Tool{
    
public:

	string _inputFilename;
	string _inputFilenameContig;
	//string _inputFilenameContig_fasta;
	string _inputDir;
	float _minimizerDensity;
    size_t _minimizerSize;
    size_t _kminmerSize;
	size_t _kminmerSizeFirst;
	bool _isFirstPass;
	bool _isOutputFasta;
	int _nbCores;

	string _filename_outputContigs;
	string _filename_kminmerSequences;
	MDBG* _mdbg;
	//MinimizerParser* _minimizerParser;
	
	unordered_map<u_int32_t, DnaBitset*> _nodeName_entire;
	unordered_map<u_int32_t, DnaBitset*> _nodeName_right;
	unordered_map<u_int32_t, DnaBitset*> _nodeName_left;

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
		(ARG_OUTPUT_DIR, "", cxxopts::value<string>())
		(ARG_NB_CORES, "", cxxopts::value<int>()->default_value(NB_CORES_DEFAULT));

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
			_nbCores = result[ARG_NB_CORES].as<int>();
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
		gzread(file_parameters, (char*)&_kminmerSizeFirst, sizeof(_kminmerSizeFirst));
		gzclose(file_parameters);

		_kminmerSize = _kminmerSizeFirst;

		cout << endl;
		cout << "Input dir: " << _inputDir << endl;
		//cout << "Output filename: " << _outputFilename << endl;
		cout << "Minimizer length: " << _minimizerSize << endl;
		cout << "Kminmer length: " << _kminmerSize << endl;
		cout << "Density: " << _minimizerDensity << endl;
		cout << endl;

		//_filename_outputContigs = _inputFilenameContig + ".fasta.gz"; //_inputDir + "/tmpContigs.fasta.gz";
		//_minimizerParser = new MinimizerParser(_minimizerSize, _minimizerDensity);
	}


    void execute (){
		
		cout << "Loading mdbg" << endl;
		string mdbg_filename = _inputDir + "/kminmerData_min_init.txt";
		_mdbg = new MDBG(_kminmerSize);
		_mdbg->load(mdbg_filename, false);
		cout << "MDBG nodes: " << _mdbg->_dbg_nodes.size() << endl;

		//cout << _inputFilenameContig_fasta << endl;
		loadContigs_min(_inputFilenameContig);

		extractKminmerSequences();
		//lalalala();

		//exit(1);
		//loadKminmerSequences();
		createBaseContigs(_inputFilenameContig, _filename_outputContigs);
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



	struct Contig{
		u_int64_t _readIndex;
		vector<u_int32_t> _nodepath;
		vector<u_int32_t> _nodepath_sorted;
	};

	vector<Contig> _contigs;
	static bool ContigComparator_ByLength(const Contig &a, const Contig &b){

		if(a._nodepath.size() == b._nodepath.size()){
			for(size_t i=0; i<a._nodepath.size() && i<b._nodepath.size(); i++){
				if(BiGraph::nodeIndex_to_nodeName(a._nodepath[i]) == BiGraph::nodeIndex_to_nodeName(b._nodepath[i])){
					continue;
				}
				else{
					return BiGraph::nodeIndex_to_nodeName(a._nodepath[i]) > BiGraph::nodeIndex_to_nodeName(b._nodepath[i]);
				}
			}
		}


		return a._nodepath.size() > b._nodepath.size();
	}

	unordered_set<u_int32_t> _invalidContigIndex;
	u_int64_t _nbContigs;
	unordered_map<u_int32_t, u_int32_t> _kminmerCounts;

	void loadContigs_min(const string& contigFilename){

		_nbContigs = 0;
		
		cout << "Extracting kminmers: " << contigFilename << endl;
		KminmerParser parser(contigFilename, _minimizerSize, _kminmerSize, false, false);
		auto fp = std::bind(&ToBasespaceNoCorrection::loadContigs_min_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
		parser.parseMinspace(fp);

		cout << "Nb contigs: " << _nbContigs << endl;

		double nbRepeatedKminmers = 0;
		for(auto& it : _kminmerCounts){
			if(it.second > 1){
				nbRepeatedKminmers += 1;
			}
		}
		cout << "Nb kminmer: " << _kminmerCounts.size() << endl;
		cout << "Repeated kminmer: " << nbRepeatedKminmers << endl;
		cout << "Repeated kminmer rate: " << (nbRepeatedKminmers / _kminmerCounts.size()) << endl;

		
		std::sort(_contigs.begin(), _contigs.end(), ContigComparator_ByLength);

		for(size_t i=0; i<_contigs.size(); i++){
			
			if(_invalidContigIndex.find(_contigs[i]._readIndex) != _invalidContigIndex.end()) continue;

			for(long j=_contigs.size()-1; j>=i+1; j--){
			
				if(_invalidContigIndex.find(_contigs[j]._readIndex) != _invalidContigIndex.end()) continue;

				double nbShared = Utils::computeSharedElements(_contigs[i]._nodepath_sorted, _contigs[j]._nodepath_sorted);
				double sharedRate_1 = nbShared / _contigs[i]._nodepath_sorted.size();
				double sharedRate_2 = nbShared / _contigs[j]._nodepath_sorted.size();

				if(sharedRate_1 > 0.33 || sharedRate_2 > 0.33){

					cout << _contigs[j]._nodepath_sorted.size() << " " << sharedRate_1 << " " << sharedRate_2 << endl;
					_invalidContigIndex.insert(_contigs[j]._readIndex);
					//break;

				}

			}
		}
		


		_kminmerCounts.clear();

		cout << _nbContigs << " " << _invalidContigIndex.size() << endl;
 		cout << "Nb contigs (no duplicate): " << (_nbContigs-_invalidContigIndex.size()) << endl;

		cout << "Extracting kminmers: " << contigFilename << endl;
		KminmerParser parser2(contigFilename, _minimizerSize, _kminmerSize, false, false);
		auto fp2 = std::bind(&ToBasespaceNoCorrection::loadContigs_min_read2, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
		parser2.parseMinspace(fp2);

		nbRepeatedKminmers = 0;
		for(auto& it : _kminmerCounts){
			if(it.second > 1){
				nbRepeatedKminmers += 1;
			}
		}
		cout << "Nb kminmer: " << _kminmerCounts.size() << endl;
		cout << "Repeated kminmer: " << nbRepeatedKminmers << endl;
		cout << "Repeated kminmer rate: " << (nbRepeatedKminmers / _kminmerCounts.size()) << endl;

		_contigs.clear();
		_kminmerCounts.clear();

	}


	void loadContigs_min_read(const vector<u_int64_t>& readMinimizers, const vector<ReadKminmerComplete>& kminmersInfos, u_int64_t readIndex){

		vector<u_int32_t> nodepath;

		for(size_t i=0; i<kminmersInfos.size(); i++){
			
			const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];

			KmerVec vec = kminmerInfo._vec;
			
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()){
				cout << "Not found kminmer" << endl;
				//getchar();
				continue;
			}



			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
			nodepath.push_back(nodeName);
			
			_kminmerCounts[nodeName] += 1;

		}

		_nbContigs += 1;

		vector<u_int32_t> nodepath_sorted = nodepath;
		std::sort(nodepath_sorted.begin(), nodepath_sorted.end());
		_contigs.push_back({readIndex, nodepath, nodepath_sorted});

	}

	void loadContigs_min_read2(const vector<u_int64_t>& readMinimizers, const vector<ReadKminmerComplete>& kminmersInfos, u_int64_t readIndex){

		if(_invalidContigIndex.find(readIndex) != _invalidContigIndex.end()) return;

		for(size_t i=0; i<kminmersInfos.size(); i++){
			
			const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];

			KmerVec vec = kminmerInfo._vec;
			
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()){
				cout << "Not found kminmer" << endl;
				//getchar();
				continue;
			}



			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
			
			_kminmerCounts[nodeName] += 1;


			bool orientation = !kminmerInfo._isReversed;

			if(i == 0){
				_nodeName_entire[nodeName] = nullptr;
			}
			else {
				if(orientation){ //+
					_nodeName_right[nodeName] = nullptr;
				}
				else{ //-
					_nodeName_left[nodeName] = nullptr;
				}
			}
			

		}


	}




	void extractKminmerSequences (){

		cout << "Extracting kminmer sequences" << endl;
		ReadParserParallel readParser(_inputFilename, false, false, _nbCores);
		readParser.parse(ExtractKminmerSequenceFunctor(_minimizerSize, _minimizerDensity, *this));
	}


	class ExtractKminmerSequenceFunctor{

		public:

		ToBasespaceNoCorrection& _toBasespace;
		EncoderRLE _encoderRLE;
		MinimizerParser* _minimizerParser;
		size_t _minimizerSize;
		float _minimizerDensity;
		

		ExtractKminmerSequenceFunctor(size_t minimizerSize, float minimizerDensity, ToBasespaceNoCorrection& toBasespace) : _toBasespace(toBasespace){
			_minimizerSize = minimizerSize;
			_minimizerDensity = minimizerDensity;
			_minimizerParser = new MinimizerParser(minimizerSize, minimizerDensity);
		}

		ExtractKminmerSequenceFunctor(const ExtractKminmerSequenceFunctor& copy) : _toBasespace(copy._toBasespace){
			_minimizerSize = copy._minimizerSize;
			_minimizerDensity = copy._minimizerDensity;
			_minimizerParser = new MinimizerParser(_minimizerSize, _minimizerDensity);
		}

		~ExtractKminmerSequenceFunctor(){
			delete _minimizerParser;
		}

		void operator () (const Read& read) {
			u_int64_t readIndex = read._index;

			//if(readIndex % 1000 == 0){
			//	cout << "Correcting kminmer " << readIndex << endl;
			//}
			//ottalSize += strlen(read->seq.s);
						
			string kminmerSequence;
			const char* sequenceOriginal = read._seq.c_str();

			string rleSequence;
			vector<u_int64_t> rlePositions;
			_encoderRLE.execute(read._seq.c_str(), read._seq.size(), rleSequence, rlePositions);

			vector<u_int64_t> minimizers;
			vector<u_int64_t> minimizers_pos;
			_minimizerParser->parse(rleSequence, minimizers, minimizers_pos);

			vector<KmerVec> kminmers; 
			vector<ReadKminmer> kminmersInfo;
			MDBG::getKminmers(_minimizerSize, _toBasespace._kminmerSize, minimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, readIndex, false);

			for(size_t i=0; i<kminmers.size(); i++){
				if(_toBasespace._mdbg->_dbg_nodes.find(kminmers[i]) == _toBasespace._mdbg->_dbg_nodes.end()) continue;

				u_int32_t nodeName = _toBasespace._mdbg->_dbg_nodes[kminmers[i]]._index;
				ReadKminmer& kminmerInfo = kminmersInfo[i];


				#pragma omp critical(indexKminmer)
				{
					if(_toBasespace._nodeName_entire.find(nodeName) != _toBasespace._nodeName_entire.end() && _toBasespace._nodeName_entire[nodeName] == nullptr){
						
						extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Entire, kminmerSequence);
						_toBasespace._nodeName_entire[nodeName] = new DnaBitset(kminmerSequence);
						//addKminmerSequenceVariant_add_all(nodeName, readIndex, _toBasespace._kminmerSequenceCopies_all_entire, kminmerSequence, _toBasespace._isDoneNodeName_entire, _toBasespace._kminmerSequence_entire, _toBasespace._kminmerSequence_entire_multi);
					}
					
					if(_toBasespace._nodeName_left.find(nodeName) != _toBasespace._nodeName_left.end() && _toBasespace._nodeName_left[nodeName] == nullptr){
						extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Left, kminmerSequence);
						_toBasespace._nodeName_left[nodeName] = new DnaBitset(kminmerSequence);
						//addKminmerSequenceVariant_add_all(nodeName, readIndex, _toBasespace._kminmerSequenceCopies_all_left, kminmerSequence, _toBasespace._isDoneNodeName_left, _toBasespace._kminmerSequence_left, _toBasespace._kminmerSequence_left_multi);
					}
					if(_toBasespace._nodeName_right.find(nodeName) != _toBasespace._nodeName_right.end() && _toBasespace._nodeName_right[nodeName] == nullptr){
						extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Right, kminmerSequence);
						_toBasespace._nodeName_right[nodeName] = new DnaBitset(kminmerSequence);
						//addKminmerSequenceVariant_add_all(nodeName, readIndex, _toBasespace._kminmerSequenceCopies_all_right, kminmerSequence, _toBasespace._isDoneNodeName_right, _toBasespace._kminmerSequence_right, _toBasespace._kminmerSequence_right_multi);
					}
				}


			}


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
			}
			

			sequence = sequence.substr(startPosition, len);

		}

	};

	void createBaseContigs(const string& contigFilename, const string& outputFilename){

		cout << "Creating basespace contigs: " << contigFilename << " " << outputFilename << endl;

		_contigIndex = 0;
		if(_isOutputFasta){
			_basespaceContigFile = gzopen(outputFilename.c_str(),"wb");
		}
		else{
			_contigFile_bitset = ofstream(contigFilename + ".bitset");
		}

		KminmerParser parser(contigFilename, _minimizerSize, _kminmerSize, false, false);
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
	u_int64_t _contigIndex;

	void createBaseContigs_read(const vector<u_int64_t>& readMinimizers, const vector<ReadKminmerComplete>& kminmersInfos, u_int64_t readIndex){


		if(_invalidContigIndex.find(readIndex) != _invalidContigIndex.end()) return;
		//cout << readIndex << " " << kminmersInfos.size() << endl;
		
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
					if(_nodeName_entire.find(nodeName) == _nodeName_entire.end() || _nodeName_entire[nodeName] == nullptr){
						cout << "not found entire " << nodeName << endl;
						continue;
					}
					
					char* seq = _nodeName_entire[nodeName]->to_string();
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
					
					if(_nodeName_entire.find(nodeName) == _nodeName_entire.end() || _nodeName_entire[nodeName] == nullptr){
						cout << "not found entire RC " << nodeName << endl;
						continue;
					}

					char* seq = _nodeName_entire[nodeName]->to_string();
					string kminmerSequence = string(seq);
					free(seq);

					Utils::revcomp(kminmerSequence);
					contigSequence += kminmerSequence;
					//cout << contigSequence << endl;
				}
			}
			else {
				if(orientation){
					
					if(_nodeName_right.find(nodeName) == _nodeName_right.end() || _nodeName_right[nodeName] == nullptr){
						cout << "not found right " << nodeName << endl;
						continue;
					}
					//cout << (_nodeName_right.find(nodeName) != _nodeName_right.end()) << endl;
					//cout << (_nodeName_right[nodeName]._sequence == nullptr) << endl;
					//cout << nodeName << endl;

					char* seq = _nodeName_right[nodeName]->to_string();
					string kminmerSequence = string(seq);
					free(seq);

					//string kminmerSequence = _nodeName_right[nodeName];
					contigSequence += kminmerSequence;
					

				}
				else{
					
					if(_nodeName_left.find(nodeName) == _nodeName_left.end() || _nodeName_left[nodeName] == nullptr){
						cout << "not found left " << nodeName << endl;
						continue;
					}
					//cout << (_nodeName_left.find(nodeName) != _nodeName_left.end()) << endl;
					//cout << (_nodeName_left[nodeName]._sequence == nullptr) << endl;
					//cout << nodeName << endl;

					char* seq = _nodeName_left[nodeName]->to_string();
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

		
		if(contigSequence.size() == 0){
			cout << "empty contig" << endl;
			return;
		}

		if(_isOutputFasta){
			string header = ">ctg" + to_string(_contigIndex) + '\n';
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

		_contigIndex += 1;

	}


};	


#endif 


