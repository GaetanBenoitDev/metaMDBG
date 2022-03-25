

/*
Tester la vitesse de racon pour verifier que je suis sur la bonne voix
1) Compter les kminmer dans les contigs
2) Si un kminmers est répété, on veut récolté autant de miode
3) Recuperation des modeles: Indexer les readIndex, selectioner le readIndex qui suporte le plus de noeuds du contig
4) Utiliser edlib pour récupérer les 20 copies les plus proches
*/

#ifndef MDBG_METAG_TOBASESPACE
#define MDBG_METAG_TOBASESPACE

#include "../Commons.hpp"
#include "../utils/edlib.h"
#include "../utils/spoa/include/spoa/spoa.hpp"
#include "../utils/DnaBitset.hpp"
//#include <seqan/align.h>
//#include <seqan/graph_msa.h>
//#include <cstring>

struct KminmerSequenceVariant{
	u_int16_t _editDistance;
	string _sequence;
};

struct KminmerSequenceVariant_Comparator {
	bool operator()(KminmerSequenceVariant const& p1, KminmerSequenceVariant const& p2)
	{
		return p1._editDistance < p2._editDistance;
	}
};

typedef priority_queue<KminmerSequenceVariant, vector<KminmerSequenceVariant>, KminmerSequenceVariant_Comparator> VariantQueue;

class ToBasespace : public Tool{
    
public:

	string _inputFilename;
	string _inputFilenameContig;
	string _inputDir;
	float _minimizerDensity;
    size_t _minimizerSize;
    size_t _kminmerSize;
    size_t _kminmerSizeFirst;
	bool _isFirstPass;
	bool _isOutputFasta;

	string _filename_outputContigs;
	string _filename_kminmerSequences;
	MDBG* _mdbg;
	MinimizerParser* _minimizerParser;
	string _truthInputFilename;
	int _nbCores;
	
	//unordered_map<ContigNode, string> _debug_node_sequences;
	//unordered_map<u_int32_t, KminmerSequence> _nodeName_entire;
	//unordered_map<u_int32_t, KminmerSequence> _nodeName_right;
	//unordered_map<u_int32_t, KminmerSequence> _nodeName_left;
	//unordered_map<u_int32_t, vector<KminmerSequence>> _nodeName_entire_multi;
	//unordered_map<u_int32_t, vector<KminmerSequence>> _nodeName_right_multi;
	//unordered_map<u_int32_t, vector<KminmerSequence>> _nodeName_left_multi;

	//unordered_map<u_int32_t, DnaBitset*> _nodeName_all_right;
	//unordered_map<u_int32_t, DnaBitset*> _nodeName_all_left;
	//unordered_map<ContigNode, DnaBitset*> _kminmerSequences_entire;
	//unordered_map<ContigNode, DnaBitset*> _kminmerSequences_left;
	//unordered_map<ContigNode, DnaBitset*> _kminmerSequences_right;

	//unordered_set<u_int32_t> _requiredCopiers_entire;
	//unordered_set<u_int32_t> _requiredCopiers_left;
	//unordered_set<u_int32_t> _requiredCopiers_right;

	unordered_map<u_int32_t, vector<DnaBitset*>> _kminmerSequenceCopies_all_entire;
	unordered_map<u_int32_t, vector<DnaBitset*>> _kminmerSequenceCopies_all_left;
	unordered_map<u_int32_t, vector<DnaBitset*>> _kminmerSequenceCopies_all_right;
	
	unordered_set<u_int32_t> isKminmerRepeated;

	unordered_map<ReadNodeName, DnaBitset*> _repeatedKminmerSequence_entire;
	unordered_map<ReadNodeName, DnaBitset*> _repeatedKminmerSequence_left;
	unordered_map<ReadNodeName, DnaBitset*> _repeatedKminmerSequence_right;

	unordered_map<u_int32_t, DnaBitset*> _kminmerSequence_entire;
	unordered_map<u_int32_t, DnaBitset*> _kminmerSequence_left;
	unordered_map<u_int32_t, DnaBitset*> _kminmerSequence_right;
	//unordered_map<u_int32_t, vector<DnaBitset*>> _kminmerSequenceCopies_left;
	//unordered_map<u_int32_t, vector<DnaBitset*>> _kminmerSequenceCopies_right;


	//unordered_map<u_int32_t, vector<DnaBitset*>> _kminmerSequenceCopies_all_left;
	//unordered_map<u_int32_t, vector<DnaBitset*>> _kminmerSequenceCopies_all_right;


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

	void parseArgs(int argc, char* argv[]){


		cxxopts::Options options("ToBasespace", "");
		options.add_options()
		(ARG_INPUT_FILENAME, "", cxxopts::value<string>())
		(ARG_INPUT_FILENAME_CONTIG, "", cxxopts::value<string>())
		(ARG_INPUT_FILENAME_TRUTH, "", cxxopts::value<string>()->default_value(""))
		(ARG_OUTPUT_FILENAME, "", cxxopts::value<string>())
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
			_isFirstPass = result[ARG_FIRST_PASS].as<bool>();
			_isOutputFasta = result[ARG_FASTA].as<bool>();
			_filename_outputContigs = result[ARG_OUTPUT_FILENAME].as<string>();
			_truthInputFilename = result[ARG_INPUT_FILENAME_TRUTH].as<string>();
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

		//_filename_outputContigs = _inputDir + "/contigs.fasta.gz";
		//_filename_outputContigs = _inputFilenameContig + ".fasta.gz"; //_inputDir + "/tmpContigs.fasta.gz";
		_minimizerParser = new MinimizerParser(_minimizerSize, _minimizerDensity);
	}

	vector<vector<u_int64_t>> _unitigDatas;
	ofstream _contigFileSupported;
	ifstream _contigFileSupported_input;

    void execute (){
		




		//_kminmerSize = 4;

		cout << "Loading mdbg" << endl;
		string mdbg_filename = _inputDir + "/mdbg_nodes_init.gz";
		_mdbg = new MDBG(_kminmerSize);
		_mdbg->load(mdbg_filename);
		cout << "MDBG nodes: " << _mdbg->_dbg_nodes.size() << endl;

		if(_truthInputFilename != ""){
			extract_truth_kminmers();
		}

		cout << "Indexing reads" << endl;
		_unitigDatas.resize(_mdbg->_dbg_nodes.size());
		indexReads();

		//cout << "Loading original mdbg" << endl;
		//string mdbg_filename = _inputDir + "/mdbg_nodes_init.gz";
		//_mdbgInit = new MDBG(_originalKminmerSize);
		//_mdbgInit->load(mdbg_filename);
		//cout << "MDBG nodes: " << _mdbgInit->_dbg_nodes.size() << endl;

		loadContigs_min(_inputFilenameContig);
		collectBestSupportingReads(_inputFilenameContig);
		//extractKminmerSequences_all();
		_unitigDatas.clear();

		//cout << isKminmerRepeated.size() << endl;
		//getchar();

		extractKminmerSequences_allVariants();

		//cout << _nodeName_all_left.size() << " " << _nodeName_all_right.size() << endl;

		cout << "Correcting kminmer sequences" << endl;

		//auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -3);  // linear gaps
		//spoa::Graph graph{};
		//_alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -3);  // linear gaps
		//spoa::Graph graph{};
		

		endCorrection(_kminmerSequenceCopies_all_entire, _isDoneNodeName_entire, _kminmerSequence_entire, _repeatedKminmerSequence_entire);
		endCorrection(_kminmerSequenceCopies_all_left, _isDoneNodeName_left, _kminmerSequence_left, _repeatedKminmerSequence_left);
		endCorrection(_kminmerSequenceCopies_all_right, _isDoneNodeName_right, _kminmerSequence_right, _repeatedKminmerSequence_right);

		isKminmerRepeated.clear();
		/*
		for(auto& it : _kminmerSequenceCopies_all_left){
			u_int32_t nodeName = it.first;
			const vector<DnaBitset*>& dnaSeq = it.second;

			if(_isDoneNodeName_left.find(nodeName) == _isDoneNodeName_left.end()) continue;

			//if(_isDoneNodeName_left.find(nodeName) == _isDoneNodeName_left.end()){
			string correctedSequence;
			performErrorCorrection_all(nodeName, dnaSeq, correctedSequence);
			_kminmerSequence_left[nodeName] = new DnaBitset(correctedSequence);
			//writeKminmerSequence_all(nodeName, correctedSequence, _outputFile_left);

			//}
			for(DnaBitset* dna : dnaSeq){
				delete dna;
			}
			it.second.clear();
		}

		
		for(auto& it : _kminmerSequenceCopies_all_right){
			u_int32_t nodeName = it.first;
			const vector<DnaBitset*>& dnaSeq = it.second;

			if(_isDoneNodeName_right.find(nodeName) == _isDoneNodeName_right.end()) continue;

			//if(_isDoneNodeName_right.find(nodeName) == _isDoneNodeName_right.end()){
			string correctedSequence;
			performErrorCorrection_all(nodeName, dnaSeq, correctedSequence);
			_kminmerSequence_right[nodeName] = new DnaBitset(correctedSequence);
			//writeKminmerSequence_all(nodeName, correctedSequence, _outputFile_right);
			//}
			for(DnaBitset* dna : dnaSeq){
				delete dna;
			}
			it.second.clear();
		}
		*/
		
		createBaseContigs(_inputFilenameContig, _filename_outputContigs.c_str());

		//delete _mdbg;
		/*
		gzclose(_outputFile_left);
		gzclose(_outputFile_right);


		return;

		//loadContigs(_inputDir + "/minimizer_contigs_complete.gz");
		//loadContigs(_inputDir + "/eval/composition/22/debug_longUnitigs.gz");
		//loadContigs(_inputDir + "/eval/composition//3/debug_longUnitigsNeighbors.gz");

		cout << "TODO: remove ContigNode and _requiredCopiers_entire, instead use unorderedmap nodeName => vector<(ReadIndex, string)>" << endl;




		extractKminmerSequences();
		lalalala();
		delete _mdbg;

		//exit(1);
		//loadKminmerSequences();
		createBaseContigs(_inputDir + "/minimizer_contigs.gz", _filename_outputContigs.c_str());
		//createBaseContigs(_inputDir + "/minimizer_contigs_complete.gz", _filename_outputContigs.c_str());
		//createBaseContigs(_inputDir + "/eval/composition//22//debug_longUnitigs.gz", _inputDir + "/eval/composition//22//debug_longUnitigs.fasta.gz");
		//createBaseContigs(_inputDir + "/eval/composition//3//debug_longUnitigsNeighbors.gz", _inputDir + "/eval/composition//3/debug_longUnitigsNeighbors.fasta.gz");

		cout << endl << "Contig filename: " << _filename_outputContigs << endl;
		*/
	}

	void endCorrection(auto& copies, auto& isDoneNodeName, auto& correctedSequences, auto& repeatedKminmers){

		/*
		cout << "-----" << endl;
		for(const auto& it : repeatedKminmers){
			cout << (it.second != nullptr) << endl;
			if(it.second == nullptr){
				cout << it.first._nodeIndex << " " << it.first._supportingReadIndex << endl;
				getchar();
			}

			if(it.first._nodeIndex == 602){
				cout << (it.second != nullptr) << endl;
				getchar();
			}
		}
		*/

		vector<u_int32_t> nodeNames;
		for(auto& it : copies){
			if(isDoneNodeName.find(it.first) != isDoneNodeName.end()) continue;
			nodeNames.push_back(it.first);
		}
		
		cout << nodeNames.size() << endl;
		#pragma omp parallel num_threads(_nbCores)
		{

			ExtractKminmerSequenceFunctor functor(_minimizerSize, _minimizerDensity, *this);

			#pragma omp for
			for(size_t i=0; i<nodeNames.size(); i++){
				
				u_int32_t nodeName = nodeNames[i];
				if(isKminmerRepeated.find(nodeName) != isKminmerRepeated.end()) continue;

				vector<DnaBitset*>* dnaSeq;

				#pragma omp critical
				dnaSeq = &copies[nodeName];


				string correctedSequence;
				functor.performErrorCorrection_all(nodeName, *dnaSeq, correctedSequence);
				//writeKminmerSequence_all(nodeName, correctedSequence, _outputFile_left);

				//}
				
				#pragma omp critical
				{
					
					correctedSequences[nodeName] = new DnaBitset(correctedSequence);

					for(DnaBitset* dna : *dnaSeq){
						delete dna;
					}
					//dnaSeq.clear();

					//if(nodeName == 10042){
					//	cout << "omg" << endl;
					//	getchar();
					//}
				}
			}
		}

		cout << correctedSequences.size() << endl;




		vector<ReadNodeName> readNodeNames;
		for(const auto& it : repeatedKminmers){
			readNodeNames.push_back(it.first);
		}
		
		cout << "Correcting repeated kminmers: " << readNodeNames.size() << endl;
		#pragma omp parallel num_threads(1)
		{

			ExtractKminmerSequenceFunctor functor(_minimizerSize, _minimizerDensity, *this);

			#pragma omp for
			for(size_t i=0; i<readNodeNames.size(); i++){
				
				ReadNodeName readNodeName = readNodeNames[i];
				
				//cout << readNodeName._nodeIndex << endl;

				DnaBitset* sequenceModel = nullptr;

				#pragma omp critical
				sequenceModel = repeatedKminmers[readNodeName];


				//cout << (sequenceModel != nullptr) << endl;
				char* sequenceModelStr = sequenceModel->to_string();

				VariantQueue queue;
				const vector<DnaBitset*>& kminmerCopies = copies[readNodeName._nodeIndex];

				for(DnaBitset* dnaSeq : kminmerCopies){

					static EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0);

					//cout << (dnaSeq != nullptr) << " " << dnaSeq->m_len << endl;
					char* dnaSeqStr = dnaSeq->to_string();
					string sequence = string(dnaSeqStr, dnaSeq->m_len);
					//cout << string(sequenceModelStr).size() << " " << sequence.size() << endl;

					EdlibAlignResult result = edlibAlign(sequenceModelStr, sequenceModel->m_len, sequence.c_str(), sequence.size(), config);
					free(dnaSeqStr);

					if (result.status != EDLIB_STATUS_OK){
						edlibFreeAlignResult(result);
						cout << "Invalid edlib status" << endl;
						continue;
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
						queue.push({result.editDistance, sequence});
					}
					else{
						if(result.editDistance < queue.top()._editDistance){

							const KminmerSequenceVariant& variant = queue.top();
							//delete variant._sequence;

							queue.pop();
							queue.push({result.editDistance, sequence});
						}
					}
					
					edlibFreeAlignResult(result);

				}

				free(sequenceModelStr);

				vector<DnaBitset*> sequences;

				while(!queue.empty()){
					KminmerSequenceVariant variant = queue.top();
					sequences.push_back(new DnaBitset(variant._sequence));
					queue.pop();
				}

				string correctedSequence;
				functor.performErrorCorrection_all(readNodeName._nodeIndex, sequences, correctedSequence);
				
				#pragma omp critical
				{
					delete repeatedKminmers[readNodeName];
					repeatedKminmers[readNodeName] = new DnaBitset(correctedSequence);
				}

				for(DnaBitset* dnaSeq : sequences){
					delete dnaSeq;
				}
				sequences.clear();

				/*
				vector<DnaBitset*>* dnaSeq;

				#pragma omp critical
				dnaSeq = &copies[nodeName];


				string correctedSequence;
				functor.performErrorCorrection_all(nodeName, *dnaSeq, correctedSequence);
				//writeKminmerSequence_all(nodeName, correctedSequence, _outputFile_left);

				//}
				#pragma omp critical
				{
					
					correctedSequences[nodeName] = new DnaBitset(correctedSequence);

					for(DnaBitset* dna : *dnaSeq){
						delete dna;
					}
					//dnaSeq.clear();
				}
				*/
			}
		}

	}

	/*
	void writeKminmerSequence(u_int32_t nodeName, unordered_map<u_int32_t, VariantQueue>& variants){

		_isDoneCorrection.insert()
		string correctedSequence;
		performErrorCorrection(nodeName, dnaSeq, _kminmerSequenceCopies_all_left[nodeName], correctedSequence, alignment_engine, graph);

		writeKminmerSequence_all(nodeName, correctedSequence, outputFile_left);


	}*/


	/*
	void loadContigs(const string& contigFilename){

		cout << "Loading mContigs: " << contigFilename << endl;
		gzFile contigFile = gzopen(contigFilename.c_str(),"rb");
		u_int64_t nbContigs = 0;
		
		while(true){

			vector<u_int32_t> nodePath;
			//vector<u_int64_t> supportingReads;
			u_int64_t size;
			gzread(contigFile, (char*)&size, sizeof(size));
			

			if(gzeof(contigFile)) break;

			nodePath.resize(size);
			//supportingReads.resize(size);
			gzread(contigFile, (char*)&nodePath[0], size * sizeof(u_int32_t));
			//gzread(contigFile, (char*)&supportingReads[0], size * sizeof(u_int64_t));

			//for(u_int32_t nodeIndex : nodePath){
			for(size_t i=0; i<nodePath.size(); i++){
				u_int32_t nodeIndex = nodePath[i];
				bool orientation;
				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex, orientation);

				if(i == 0){
					_kminmerSequenceCopies_all_entire[nodeName] = {};
				}
				else {
					if(orientation){ //+
						_kminmerSequenceCopies_all_right[nodeName] = {};
					}
					else{ //-
						_kminmerSequenceCopies_all_left[nodeName] = {};
					}
				}

			}

			nbContigs += 1;
			//cout << nodePath.size() << endl;
		}

		gzclose(contigFile);

		cout << "Nb contigs: " << nbContigs << endl;
	}
	*/
	//MDBG* _mdbgInit;
	//size_t _originalKminmerSize;
	u_int64_t _nbContigs;
	unordered_map<u_int32_t, u_int32_t> _kminmerCounts;

	/*
	void loadContigs_fasta (const string& contigFilename){


		_nbContigs = 0;

		auto fp = std::bind(&ToBasespace::loadContigs_fasta_read, this, std::placeholders::_1, std::placeholders::_2);
		ReadParser readParser(contigFilename, true, false);
		readParser.parse(fp);

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

		_kminmerCounts.clear();
	}
	*/





	void indexReads(){

		KminmerParser parser(_inputDir + "/read_data_init.txt", _minimizerSize, _kminmerSize, false);
		auto fp = std::bind(&ToBasespace::indexReads_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
		parser.parse(fp);
	}

	void indexReads_read(const vector<u_int64_t>& minimizers, const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex){//}, const vector<KmerVec>& kminmers_k3, const vector<ReadKminmer>& kminmersInfos_k3){

		for(const KmerVec& vec : kminmers){
			
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()) continue;

			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
			_unitigDatas[nodeName].push_back(readIndex);
		}
	}

	void loadContigs_min(const string& contigFilename){

		_nbContigs = 0;

		cout << "Extracting kminmers (contigs)" << endl;
		KminmerParser parser(contigFilename, _minimizerSize, _kminmerSize, false);
		auto fp = std::bind(&ToBasespace::loadContigs_min_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
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

	}

	
	void loadContigs_min_read(const vector<u_int64_t>& readMinimizers, const vector<ReadKminmerComplete>& kminmersInfos, u_int64_t readIndex){

		//cout << readIndex << " " << kminmersInfos.size() << endl;
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

			//if(_kminmerCounts[nodeName] > 1){
			//	isKminmerRepeated.insert(nodeName);
			//}
			//vector<u_int64_t> minimizerSeq;
			
			//for(size_t i=kminmerInfo._read_pos_start; i<=kminmerInfo._read_pos_end; i++){
			//	minimizerSeq.push_back(readMinimizers[i]);
			//}
			

			//if(kminmerInfo._isReversed){
			//	std::reverse(minimizerSeq.begin(), minimizerSeq.end());
			//}

			bool orientation = !kminmerInfo._isReversed;

			if(i == 0){
				_kminmerSequenceCopies_all_entire[nodeName] = {};
			}
			else {
				if(orientation){ //+
					_kminmerSequenceCopies_all_right[nodeName] = {};
				}
				else{ //-
					_kminmerSequenceCopies_all_left[nodeName] = {};
				}
			}

		}

		_nbContigs += 1;

	}


	void collectBestSupportingReads(const string& contigFilename){

		_contigFileSupported = ofstream(contigFilename + ".tmp");

		cout << "Collecting best supporting reads" << endl;
		KminmerParser parser(contigFilename, _minimizerSize, _kminmerSize, false);
		auto fp = std::bind(&ToBasespace::collectBestSupportingReads_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
		parser.parseMinspace(fp);

		_kminmerCounts.clear();

		_contigFileSupported.close();

	}

	
	void collectBestSupportingReads_read(const vector<u_int64_t>& readMinimizers, const vector<ReadKminmerComplete>& kminmersInfos, u_int64_t readIndex){

		//cout << "----------------------" << endl;
		vector<u_int64_t> supportingReads;

		//cout << readIndex << " " << kminmersInfos.size() << endl;
		for(size_t i=0; i<kminmersInfos.size(); i++){
			
			const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];

			KmerVec vec = kminmerInfo._vec;
			
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()){
				cout << "Not found kminmer" << endl;
				supportingReads.push_back(-1);
				//getchar();
				continue;
			}

			bool orientation = !kminmerInfo._isReversed;

			
			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;

			//_kminmerCounts.clear();
			if(_kminmerCounts[nodeName] > 1){
				isKminmerRepeated.insert(nodeName);
				u_int64_t readIndex = getBestSupportingRead(nodeName, i, kminmersInfos);
				supportingReads.push_back(readIndex);
				
				ReadNodeName readNodeName = {nodeName, readIndex};

				if(i == 0){
					_repeatedKminmerSequence_entire[readNodeName] = nullptr;
				}
				else {
					if(orientation){ //+
						//cout << "right: " << nodeName << " " << readIndex << endl;
						_repeatedKminmerSequence_right[readNodeName] = nullptr;
					}
					else{ //-
						//cout << "left: " << nodeName << " " << readIndex << endl;
						_repeatedKminmerSequence_left[readNodeName] = nullptr;
					}
				}

				//if(nodeName==293 && readIndex==8320){
				//	cout << "omg" << endl;
				//	getchar();
				//}

				//getchar();
				//cout << _repeatedKminmerSequence_entire.size() << " " << _repeatedKminmerSequence_right.size() << " " << _repeatedKminmerSequence_left.size() << endl;
			}
			else{
				supportingReads.push_back(-1);
			}
			//vector<u_int64_t> minimizerSeq;
			
			//for(size_t i=kminmerInfo._read_pos_start; i<=kminmerInfo._read_pos_end; i++){
			//	minimizerSeq.push_back(readMinimizers[i]);
			//}
			

			//if(kminmerInfo._isReversed){
			//	std::reverse(minimizerSeq.begin(), minimizerSeq.end());
			//}


			if(i == 0){
				_kminmerSequenceCopies_all_entire[nodeName] = {};
			}
			else {
				if(orientation){ //+
					_kminmerSequenceCopies_all_right[nodeName] = {};
				}
				else{ //-
					_kminmerSequenceCopies_all_left[nodeName] = {};
				}
			}

		}

		_nbContigs += 1;

		u_int32_t size = supportingReads.size();
		_contigFileSupported.write((const char*)&size, sizeof(size));
		_contigFileSupported.write((const char*)&supportingReads[0], size*sizeof(u_int64_t));

	}

	u_int64_t getBestSupportingRead(u_int32_t nodeName, size_t nodeNamePosition, const vector<ReadKminmerComplete>& kminmersInfos){

		//cout << "-----" << endl;
		//cout << "Size: " << kminmersInfos.size() << endl;

		//unordered_map<u_int64_t, u_int32_t> readSupports;

		u_int64_t minSharedReads = -1;

		u_int64_t readIndexMin = -1;
		
		bool isValidRight;
		bool isValidLeft;
		readIndexMin = getBestSupportingRead_direction(nodeName, nodeNamePosition+1, kminmersInfos, 1, isValidRight, minSharedReads);
		u_int64_t nodeNamePosition_left = getBestSupportingRead_direction(nodeName, nodeNamePosition-1, kminmersInfos, -1, isValidLeft, minSharedReads);

		if(nodeNamePosition_left != -1) readIndexMin = nodeNamePosition_left;
		//long dist_right = nodeNamePosition_right - nodeNamePosition;
		//long dist_left = nodeNamePosition - nodeNamePosition_left;
		//cout << nodeNamePosition << " " << nodeNamePosition_left << " " << nodeNamePosition_right << "     " << dist_left << " " << dist_right << endl;


		//u_int64_t longuestSupportingReadIndex = -1;

		/*
		if((!isValidLeft && !isValidRight) || (isValidLeft && isValidRight)){

			if(dist_right > dist_left){

				u_int32_t nodeName_right = _mdbg->_dbg_nodes[kminmersInfos[nodeNamePosition_right]._vec]._index;

				vector<u_int64_t> sharedElements;
				Utils::collectSharedElements(_unitigDatas[nodeName], _unitigDatas[nodeName_right], sharedElements);
				longuestSupportingReadIndex = sharedElements[0];
				cout << sharedElements.size() << endl;
				if(sharedElements.size() > 5) getchar();
			}
			else{
				
				u_int32_t nodeName_left = _mdbg->_dbg_nodes[kminmersInfos[nodeNamePosition_left]._vec]._index;

				vector<u_int64_t> sharedElements;
				Utils::collectSharedElements(_unitigDatas[nodeName], _unitigDatas[nodeName_left], sharedElements);
				longuestSupportingReadIndex = sharedElements[0];
				cout << sharedElements.size() << endl;
				if(sharedElements.size() > 5) getchar();
			}
		}
		else if(isValidRight){
			u_int32_t nodeName_right = _mdbg->_dbg_nodes[kminmersInfos[nodeNamePosition_right]._vec]._index;

			vector<u_int64_t> sharedElements;
			Utils::collectSharedElements(_unitigDatas[nodeName], _unitigDatas[nodeName_right], sharedElements);
			longuestSupportingReadIndex = sharedElements[0];
			cout << sharedElements.size() << endl;
			if(sharedElements.size() > 5) getchar();
		}
		else{
			u_int32_t nodeName_left = _mdbg->_dbg_nodes[kminmersInfos[nodeNamePosition_left]._vec]._index;

			vector<u_int64_t> sharedElements;
			Utils::collectSharedElements(_unitigDatas[nodeName], _unitigDatas[nodeName_left], sharedElements);
			longuestSupportingReadIndex = sharedElements[0];
			cout << sharedElements.size() << endl;
			if(sharedElements.size() > 5) getchar();
		}
		*/

		return readIndexMin;
		//cout << kminmersInfos.size() << " " << nodeNamePosition_left << " " << nodeNamePosition_right << endl;


		//cout << nodeName_left << " " << nodeName_right << endl;
		//cout << sharedElements.size() << endl;
		//return sharedElements[0];
		/*
		u_int64_t maxSupportRead = -1;

		if(readSupports.size() == 0){
			return _unitigDatas[nodeName][0];
		}
		else{

			u_int64_t maxSupport = 0;
			for(const auto& it : readSupports){
				if(it.second > maxSupport){
					maxSupport = it.second;
					maxSupportRead = it.first;
				}
			}


		}


		return maxSupportRead;
		*/

	}

	u_int64_t getBestSupportingRead_direction(u_int32_t nodeName, size_t nodeNamePosition, const vector<ReadKminmerComplete>& kminmersInfos, int inc, bool& isValid, u_int64_t& minSharedReads){
		
		u_int64_t readIndexResult = -1;

		if(inc == 1){
			//cout << "RIGHT" << endl;
		}
		else{
			//cout << "LEFT" << endl;
		}

		long i = nodeNamePosition;

		while(true){
			if(i < 0 || i >= kminmersInfos.size()){
				//cout << "invalid" << endl;
				isValid = false;
				break;
			}

			const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];

			KmerVec vec = kminmerInfo._vec;

			if(_mdbg->_dbg_nodes.find(vec) != _mdbg->_dbg_nodes.end()){

				u_int32_t nodeName2 = _mdbg->_dbg_nodes[vec]._index;
				
				vector<u_int64_t> sharedElements;
				Utils::collectSharedElements(_unitigDatas[nodeName], _unitigDatas[nodeName2], sharedElements);

				if(sharedElements.size() == 0) break;

				if(sharedElements.size() < minSharedReads){
					minSharedReads = sharedElements.size();
					readIndexResult = sharedElements[0];
				}

				/*
				if(!Utils::shareAny(_unitigDatas[nodeName], _unitigDatas[nodeName2])){
					cout << "valid" << endl;
					isValid = true;
					//isValid = false;
					break;
				}
				*/

				//if(nodeName2 == nodeName) break;
				//if(sharedElements.size() == 0) break;

				//for(u_int64_t readIndex : sharedElements){
				//	readSupports[readIndex] += 1;
				//}


				//cout << i << " " << sharedElements.size() << endl;

			}

			i += inc;
		}

		//if(isValid) return i;
		return readIndexResult;
		/*
		if(inc == 1){
			return i-1;
		}
		else{
			return i+1;
		}
		*/

	}

	/*
	void loadContigs_fasta_read(kseq_t* read, u_int64_t readIndex){
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


		
		//if(readIndex == 31){
			//minimizers.erase(minimizers.begin()+920);
			//minimizers_pos.erase(minimizers_pos.begin()+920);
		//}
		//for(size_t i=0; i<minimizers.size(); i++){
		//	cout << i << ": " << minimizers[i] << " " << minimizers_pos[i] << endl;
		//}

		



		vector<KmerVec> kminmers; 
		vector<ReadKminmer> kminmersInfo;
		MDBG::getKminmers(_minimizerSize, _kminmerSize, minimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, readIndex, false);
		

		for(size_t i=0; i<kminmers.size(); i++){
			if(_mdbg->_dbg_nodes.find(kminmers[i]) == _mdbg->_dbg_nodes.end()){
				cout << "Unknown original kminmer" << endl;
				cout << readIndex << " " << strlen(sequenceOriginal) << " " << minimizers.size() << endl;
				cout << i << endl;
				//exit(1);
				continue;
			}

			u_int32_t nodeName = _mdbg->_dbg_nodes[kminmers[i]]._index;
			//cout << nodeName << endl;
			_kminmerCounts[nodeName] += 1;

			bool orientation = true;
			if(kminmersInfo[i]._isReversed){
				orientation = false;
			}

			if(i == 0){
				_kminmerSequenceCopies_all_entire[nodeName] = {};
			}
			else {
				if(orientation){ //+
					_kminmerSequenceCopies_all_right[nodeName] = {};
				}
				else{ //-
					_kminmerSequenceCopies_all_left[nodeName] = {};
				}
			}
				
		}

		//cout << readIndex << " " << _nodeName_left.size() << " " << _nodeName_right.size() << endl;

		_nbContigs += 1;
	}
	*/



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


	void extractKminmerSequences_allVariants (){

		cout << "Extracting kminmer sequence variants" << endl;

		//auto fp = std::bind(&ToBasespace::extractKminmerSequences_allVariants_read, this, std::placeholders::_1);
		ReadParserParallel readParser(_inputFilename, false, false, _nbCores);
		readParser.parse(ExtractKminmerSequenceFunctor(_minimizerSize, _minimizerDensity, *this));

	}

	class ExtractKminmerSequenceFunctor{

		public:

		ToBasespace& _toBasespace;
		EncoderRLE _encoderRLE;
		MinimizerParser* _minimizerParser;
		size_t _minimizerSize;
		float _minimizerDensity;
		
		spoa::Graph _graph;
		std::unique_ptr<spoa::AlignmentEngine> _alignementEngine;

		ExtractKminmerSequenceFunctor(size_t minimizerSize, float minimizerDensity, ToBasespace& toBasespace) : _toBasespace(toBasespace){
			_minimizerSize = minimizerSize;
			_minimizerDensity = minimizerDensity;
			_minimizerParser = new MinimizerParser(minimizerSize, minimizerDensity);
			
			_alignementEngine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -3);  // linear gaps
		}

		ExtractKminmerSequenceFunctor(const ExtractKminmerSequenceFunctor& copy) : _toBasespace(copy._toBasespace){
			_minimizerSize = copy._minimizerSize;
			_minimizerDensity = copy._minimizerDensity;
			_minimizerParser = new MinimizerParser(_minimizerSize, _minimizerDensity);

			_alignementEngine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -3);  // linear gaps


		}

		~ExtractKminmerSequenceFunctor(){
			delete _minimizerParser;
		}

		void operator () (const Read& read) {
			u_int64_t readIndex = read._index;

			if(readIndex % 1000 == 0){
				cout << "Correcting kminmer " << readIndex << endl;
			}
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


				//3908 7040
				//298 8320
				//if(readIndex == 8320){
				//	cout << nodeName << endl;
				//}

				//cout << _kminmerSequenceCopies_all_entire.size() << " " << _kminmerSequenceCopies_all_left.size() << " " << _kminmerSequenceCopies_all_right.size() << " " << _kminmerSequence_entire.size() << " " << _kminmerSequence_left.size() << " " << _kminmerSequence_right.size() << endl;

				bool isEntire;
				bool isLeft;
				bool isRight;

				//
				//{

				//}

				#pragma omp critical
				{
					isEntire = _toBasespace._kminmerSequenceCopies_all_entire.find(nodeName) != _toBasespace._kminmerSequenceCopies_all_entire.end() && _toBasespace._isDoneNodeName_entire.find(nodeName) == _toBasespace._isDoneNodeName_entire.end();
					isLeft = _toBasespace._kminmerSequenceCopies_all_left.find(nodeName) != _toBasespace._kminmerSequenceCopies_all_left.end() && _toBasespace._isDoneNodeName_left.find(nodeName) == _toBasespace._isDoneNodeName_left.end();
					isRight = _toBasespace._kminmerSequenceCopies_all_right.find(nodeName) != _toBasespace._kminmerSequenceCopies_all_right.end() && _toBasespace._isDoneNodeName_right.find(nodeName) == _toBasespace._isDoneNodeName_right.end();

					if(!isEntire){
						ReadNodeName readNodeName = {nodeName, readIndex};
						isEntire = _toBasespace._repeatedKminmerSequence_entire.find(readNodeName) != _toBasespace._repeatedKminmerSequence_entire.end();
					}

					if(!isLeft){
						ReadNodeName readNodeName = {nodeName, readIndex};
						isLeft = _toBasespace._repeatedKminmerSequence_left.find(readNodeName) != _toBasespace._repeatedKminmerSequence_left.end();
					}

					if(!isRight){
						ReadNodeName readNodeName = {nodeName, readIndex};
						isRight = _toBasespace._repeatedKminmerSequence_right.find(readNodeName) != _toBasespace._repeatedKminmerSequence_right.end();
					}
				}

				//if(nodeName==293 && readIndex==8320){
				//	ReadNodeName readNodeName = {nodeName, readIndex};
				//	cout << (_toBasespace._repeatedKminmerSequence_entire.find(readNodeName) != _toBasespace._repeatedKminmerSequence_entire.end()) << endl;
				//	cout << "lol" << endl;
				//	cout << isEntire << endl;
				//	getchar();
				//}

				if(isEntire){
					_toBasespace.extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Entire, kminmerSequence);
					addKminmerSequenceVariant_add_all(nodeName, readIndex, _toBasespace._kminmerSequenceCopies_all_entire, kminmerSequence, _toBasespace._isDoneNodeName_entire, _toBasespace._kminmerSequence_entire, _toBasespace._repeatedKminmerSequence_entire);
				}
				
				if(isLeft){
					_toBasespace.extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Left, kminmerSequence);
					addKminmerSequenceVariant_add_all(nodeName, readIndex, _toBasespace._kminmerSequenceCopies_all_left, kminmerSequence, _toBasespace._isDoneNodeName_left, _toBasespace._kminmerSequence_left, _toBasespace._repeatedKminmerSequence_left);
				}
				if(isRight){
					_toBasespace.extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Right, kminmerSequence);
					addKminmerSequenceVariant_add_all(nodeName, readIndex, _toBasespace._kminmerSequenceCopies_all_right, kminmerSequence, _toBasespace._isDoneNodeName_right, _toBasespace._kminmerSequence_right, _toBasespace._repeatedKminmerSequence_right);
				}

				//if(nodeName==293 && readIndex==8320){
				//	ReadNodeName readNodeName = {nodeName, readIndex};
				//	//cout << (_toBasespace._repeatedKminmerSequence_left.find(readNodeName) != _toBasespace._repeatedKminmerSequence_left.end()) << " " << (_toBasespace._repeatedKminmerSequence_right.find(readNodeName) != _toBasespace._repeatedKminmerSequence_right.end()) << endl;
				//	cout << "lol" << endl;
				//	getchar();
				//}

			}

			//if(readIndex == 7210){
			//	cout << "lal" << endl;
			//	getchar();
			//}

		}


		void addKminmerSequenceVariant_add_all(u_int32_t nodeName, u_int64_t readIndex, auto& variants, const string& sequence, auto& isDoneNodeName, auto& correctedSequences, auto& repeatedKminmers){
			
			//cout << (isDoneNodeName.find(nodeName) != isDoneNodeName.end()) << endl;
			//if(sequenceModel == nullptr){
			//	cout << "pas normal" << endl;
			//	return; //model not found yet
			//}

			//cout << nodeName << " " << sequence.size() << endl;
						

			vector<DnaBitset*>* queue;

			bool cancel = false;
			//if(_isDoneCorrection.find(nodeName) != _isDoneCorrection.end()) return;
			#pragma omp critical
			{
				
				if(_toBasespace.isKminmerRepeated.find(nodeName) == _toBasespace.isKminmerRepeated.end() && isDoneNodeName.find(nodeName) != isDoneNodeName.end()) cancel = true;

				if(!cancel){
					queue = &variants[nodeName];
					queue->push_back(new DnaBitset(sequence));

					if(_toBasespace.isKminmerRepeated.find(nodeName) != _toBasespace.isKminmerRepeated.end()) cancel = true;

					//Collect all sequences for repeated kminmer
					//if(_toBasespace.isKminmerRepeated.find(nodeName) != _toBasespace.isKminmerRepeated.end()){ 
					ReadNodeName readNodeName = {nodeName, readIndex};
					if(repeatedKminmers.find(readNodeName) != repeatedKminmers.end()){
						//cout << "lala: " << nodeName << " " << readIndex << endl;
						repeatedKminmers[readNodeName] = new DnaBitset(sequence);
						//if(nodeName==293 && readIndex==8320){
							//cout << "omg!!" << endl;
							//getchar();
						//}
						cancel = true;
					}
					//}
					if(queue->size() < 20) cancel = true;
					
					//cout << nodeName << " " << queue->size() << endl;
					if(!cancel) isDoneNodeName.insert(nodeName);
				}
			}

			//return;
			if(cancel) return;


			//cout << "correcting" << endl;
			//DnaBitset* sequenceModel = models[nodeName];

			string correctedSequence;
			performErrorCorrection_all(nodeName, *queue, correctedSequence);

			//writeKminmerSequence_all(nodeName, correctedSequence, file);
			#pragma omp critical
			{
				

				correctedSequences[nodeName] = new DnaBitset(correctedSequence);
				//variants.erase(nodeName);

				//for(DnaBitset* dna : *queue){
					//delete dna;
				//}
			}



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
				//getchar();
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

			/*
			#pragma omp critical
			{
				cout << "-----" << endl;
				for(size_t i=0; i<sequences.size(); i++){ //&& processedSequences < 5

					DnaBitset* dna = sequences[i];
					//const DnaBitset* dna = variant._sequence; //sequenceCopies[s._sequenceIndex];
					char* dnaStr = dna->to_string();

					cout << dnaStr << endl;

					free(dnaStr);

					//processedSequences += 1;
				}
				getchar();
			}
			*/

			//std::reverse(sequences.begin(), sequences.end());

			for(size_t i=0; i<sequences.size(); i++){ //&& processedSequences < 5

				DnaBitset* dna = sequences[i];
				//const DnaBitset* dna = variant._sequence; //sequenceCopies[s._sequenceIndex];
				char* dnaStr = dna->to_string();

				//cout << s._sequenceIndex << " " << s._editDistance << endl;

				auto alignment = _alignementEngine->Align(dnaStr, dna->m_len, _graph);
				_graph.AddAlignment(alignment, dnaStr, dna->m_len);

				free(dnaStr);

				//processedSequences += 1;
			}


			//string consensus = _graph.GenerateConsensus();
			//correctedSequence = consensus;
			//return;

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

			float t = msa.size() * 0.5;

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

	};

	unordered_set<u_int32_t> _isDoneNodeName_entire;
	unordered_set<u_int32_t> _isDoneNodeName_left;
	unordered_set<u_int32_t> _isDoneNodeName_right;

	//void extractKminmerSequences_allVariants_read(const Read& read){


	//}






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


	gzFile _basespaceContigFile;

	/*
	void createBaseContigs (const string& contigFilename, const string& outputFilename){

		cout << "Creating basespace contigs: " << contigFilename << " " << outputFilename << endl;

		_basespaceContigFile = gzopen(outputFilename.c_str(),"wb");

		auto fp = std::bind(&ToBasespace::createBaseContigs_read, this, std::placeholders::_1, std::placeholders::_2);
		ReadParser readParser(contigFilename, true, false);
		readParser.parse(fp);

		
		gzclose(_basespaceContigFile);
	}
	*/
	
	ofstream _fileHifiasmAll;
	u_int64_t _hifiasmContigIndex;
	unordered_set<string> _hifiasmWrittenNodeNames;

	void createBaseContigs(const string& contigFilename, const string& outputFilename){

		_contigFileSupported_input = ifstream(contigFilename + ".tmp");

		_hifiasmContigIndex = 0;
		_fileHifiasmAll = ofstream(_inputDir + "/binning_results_hifiasm.csv");
		_fileHifiasmAll << "Name,Colour" << endl;

		cout << "Creating basespace contigs: " << contigFilename << " " << outputFilename << endl;

		_basespaceContigFile = gzopen(outputFilename.c_str(),"wb");

		KminmerParser parser(contigFilename, _minimizerSize, _kminmerSize, false);
		auto fp = std::bind(&ToBasespace::createBaseContigs_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
		parser.parseMinspace(fp);

		gzclose(_basespaceContigFile);
		_fileHifiasmAll.close();
		_contigFileSupported_input.close();
	}

	void createBaseContigs_read(const vector<u_int64_t>& readMinimizers, const vector<ReadKminmerComplete>& kminmersInfos, u_int64_t readIndex){


		vector<u_int64_t> supportingReads;
		u_int32_t size;
		_contigFileSupported_input.read((char*)&size, sizeof(size));
		supportingReads.resize(size);
		_contigFileSupported_input.read((char*)&supportingReads[0], size*sizeof(u_int64_t));

		//cout << readIndex << " " << kminmersInfos.size() << endl;
		
		string contigSequence = "";

		for(size_t i=0; i<kminmersInfos.size(); i++){
			
			u_int64_t readIndex = supportingReads[i];
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
				if(orientation){ //+

					char* seq;


					if(_kminmerSequence_entire.find(nodeName) != _kminmerSequence_entire.end()){
						if(_kminmerSequence_entire[nodeName] == nullptr){
							cout << "No sequence for kminmer" << endl;
							continue;
						}
						seq = _kminmerSequence_entire[nodeName]->to_string();
					}
					else{
						ReadNodeName readNodeName = {nodeName, readIndex};
						if(_repeatedKminmerSequence_entire.find(readNodeName) != _repeatedKminmerSequence_entire.end()){
							if(_repeatedKminmerSequence_entire[readNodeName] == nullptr){
								cout << "No sequence for kminmer" << endl;
								continue;
							}
						}
						seq = _repeatedKminmerSequence_entire[readNodeName]->to_string();
					}
					//else if
					//continue;
					//char* seq = _kminmerSequence_entire[nodeName]->to_string();
					string kminmerSequence = string(seq);
					free(seq);

					contigSequence += kminmerSequence;
					//cout << contigSequence << endl;
				}
				else{
					
					
					char* seq;


					if(_kminmerSequence_entire.find(nodeName) != _kminmerSequence_entire.end()){
						if(_kminmerSequence_entire[nodeName] == nullptr){
							cout << "No sequence for kminmer" << endl;
							continue;
						}
						seq = _kminmerSequence_entire[nodeName]->to_string();
					}
					else{
						ReadNodeName readNodeName = {nodeName, readIndex};
						if(_repeatedKminmerSequence_entire.find(readNodeName) != _repeatedKminmerSequence_entire.end()){
							if(_repeatedKminmerSequence_entire[readNodeName] == nullptr){
								cout << "No sequence for kminmer" << endl;
								continue;
							}
						}
						seq = _repeatedKminmerSequence_entire[readNodeName]->to_string();
					}

					string kminmerSequence = string(seq);
					free(seq);

					Utils::revcomp(kminmerSequence);
					contigSequence += kminmerSequence;
				}
			}
			else {
				if(orientation){

					char* seq;


					if(_kminmerSequence_right.find(nodeName) != _kminmerSequence_right.end()){
						if(_kminmerSequence_right[nodeName] == nullptr){
							cout << "No sequence for kminmer" << endl;
							continue;
						}
						seq = _kminmerSequence_right[nodeName]->to_string();
					}
					else{
						ReadNodeName readNodeName = {nodeName, readIndex};
						if(_repeatedKminmerSequence_right.find(readNodeName) != _repeatedKminmerSequence_right.end()){
							if(_repeatedKminmerSequence_right[readNodeName] == nullptr){
								cout << "No sequence for kminmer" << endl;
								continue;
							}
						}
						seq = _repeatedKminmerSequence_right[readNodeName]->to_string();
					}

					string kminmerSequence = string(seq);
					free(seq);
					
					contigSequence += kminmerSequence;
					

				}
				else{
					
					char* seq;


					if(_kminmerSequence_left.find(nodeName) != _kminmerSequence_left.end()){
						if(_kminmerSequence_left[nodeName] == nullptr){
							cout << "No sequence for kminmer" << endl;
							continue;
						}
						seq = _kminmerSequence_left[nodeName]->to_string();
					}
					else{
						ReadNodeName readNodeName = {nodeName, readIndex};
						if(_repeatedKminmerSequence_left.find(readNodeName) != _repeatedKminmerSequence_left.end()){
							if(_repeatedKminmerSequence_left[readNodeName] == nullptr){
								cout << "No sequence for kminmer" << endl;
								continue;
							}
						}
						seq = _repeatedKminmerSequence_left[readNodeName]->to_string();
					}

					string kminmerSequence = string(seq);
					free(seq);

					Utils::revcomp(kminmerSequence);
					contigSequence += kminmerSequence;


				}
			}
		}

		
		string header = ">ctg" + to_string(readIndex) + '\n';
		gzwrite(_basespaceContigFile, (const char*)&header[0], header.size());
		contigSequence +=  '\n';
		gzwrite(_basespaceContigFile, (const char*)&contigSequence[0], contigSequence.size());

		if(contigSequence.size() > 100000){

			if(_truthInputFilename != ""){
				
				for(size_t i=0; i<kminmersInfos.size(); i++){
					const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];
					KmerVec vec = kminmerInfo._vec;
				
					if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()){
						continue;
					}

					u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;

					if(_evaluation_hifiasmGroundTruth_nodeName_to_unitigName.find(nodeName) != _evaluation_hifiasmGroundTruth_nodeName_to_unitigName.end()){
						for(string& unitigName : _evaluation_hifiasmGroundTruth_nodeName_to_unitigName[nodeName]){
							if(_hifiasmWrittenNodeNames.find(unitigName) != _hifiasmWrittenNodeNames.end()) continue;
							_fileHifiasmAll << unitigName << "," << _hifiasmContigIndex << endl;
							_hifiasmWrittenNodeNames.insert(unitigName);
						}
					}
				}

				_hifiasmContigIndex += 1;

			}

		}

	}

	/*
	void createBaseContigs_read(kseq_t* read, u_int64_t readIndex){

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
		
		
		string contigSequence = "";


		for(size_t i=0; i<kminmers.size(); i++){
			if(_mdbg->_dbg_nodes.find(kminmers[i]) == _mdbg->_dbg_nodes.end()){
				cout << "Unknown original kminmer" << endl;
				//exit(1);
				continue;
			}

			u_int32_t nodeName = _mdbg->_dbg_nodes[kminmers[i]]._index;
			//cout << nodeName << endl;

			bool orientation = true;
			if(kminmersInfo[i]._isReversed){
				orientation = false;
			}

			if(i == 0){
				if(orientation){ //+
					if(_kminmerSequence_entire.find(nodeName) == _kminmerSequence_entire.end() || _kminmerSequence_entire[nodeName] == nullptr) continue;
					char* seq = _kminmerSequence_entire[nodeName]->to_string();
					string kminmerSequence = string(seq);
					free(seq);

					contigSequence += kminmerSequence;
					//cout << contigSequence << endl;
				}
				else{
					
					if(_kminmerSequence_entire.find(nodeName) == _kminmerSequence_entire.end() || _kminmerSequence_entire[nodeName] == nullptr) continue;
					char* seq = _kminmerSequence_entire[nodeName]->to_string();
					string kminmerSequence = string(seq);
					free(seq);

					Utils::revcomp(kminmerSequence);
					contigSequence += kminmerSequence;
				}
			}
			else {
				if(orientation){

					if(_kminmerSequence_right.find(nodeName) == _kminmerSequence_right.end() || _kminmerSequence_right[nodeName] == nullptr) continue;
					char* seq = _kminmerSequence_right[nodeName]->to_string();
					string kminmerSequence = string(seq);
					free(seq);
					
					contigSequence += kminmerSequence;
					

				}
				else{
					
					if(_kminmerSequence_left.find(nodeName) == _kminmerSequence_left.end() || _kminmerSequence_left[nodeName] == nullptr) continue;
					char* seq = _kminmerSequence_left[nodeName]->to_string();
					string kminmerSequence = string(seq);
					free(seq);

					Utils::revcomp(kminmerSequence);
					contigSequence += kminmerSequence;


				}
			}

				
		}


		string header = ">ctg" + to_string(readIndex) + '\n';
		gzwrite(_basespaceContigFile, (const char*)&header[0], header.size());
		contigSequence +=  '\n';
		gzwrite(_basespaceContigFile, (const char*)&contigSequence[0], contigSequence.size());


		//cout << readIndex << " " << _nodeName_left.size() << " " << _nodeName_right.size() << endl;
	}
	*/

	/*
	void createBaseContigs(const string& contigFilename, const string& outputFilename){

		gzFile contigFile = gzopen(contigFilename.c_str(),"rb");
		ofstream contigFile_bitset(contigFilename + ".bitset");

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
			//gzread(contigFile, (char*)&supportingReads[0], size * sizeof(u_int64_t));



			//for(u_int32_t nodeIndex : nodePath){
			for(size_t i=0; i<nodePath.size(); i++){
				
				u_int32_t nodeIndex = nodePath[i];

				

			}

			//cout << contigSequence.substr(362170, 100) << endl;
			if(contig_index % 100000 == 0){
				cout << endl << endl;
				cout << nodePath.size() << endl;
				cout << contigSequence.size() << endl;
			}
			//cout << contigSequence << endl;


			if(_isOutputFasta){
			}
			else{
				DnaBitset* contigSequenceBitset = new DnaBitset(contigSequence);
				u_int32_t sizeData = contigSequenceBitset->_bitsetSize;
				u_int32_t sizeSeq = contigSequenceBitset->m_len;
				uint8_t* m_data = contigSequenceBitset->m_data;
				contigFile_bitset.write((const char*)&sizeData, sizeof(sizeData));
				contigFile_bitset.write((const char*)&sizeSeq, sizeof(sizeSeq));
				contigFile_bitset.write((const char*)&m_data[0], sizeData*sizeof(uint8_t));
				delete contigSequenceBitset;
			}


		}

		gzclose(contigFile);
		contigFile_bitset.close();

	}
	*/




	//MinimizerParser* _minimizerParser;
	u_int32_t _extract_truth_kminmers_read_position;
	unordered_map<u_int32_t, vector<string>> _evaluation_hifiasmGroundTruth_nodeName_to_unitigName;
	vector<u_int32_t> _evaluation_hifiasmGroundTruth_path;
	unordered_set<u_int32_t> _hifiasm_startingNodenames;
	ofstream _file_groundTruth_hifiasm_position;
	EncoderRLE _encoderRLE;

	void extract_truth_kminmers_read(const Read& read){
		//ottalSize += strlen(read->seq.s);


		u_int64_t readIndex = read._index;

		string rleSequence;
		vector<u_int64_t> rlePositions;
		_encoderRLE.execute(read._seq.c_str(), read._seq.size(), rleSequence, rlePositions);

		vector<u_int64_t> minimizers;
		vector<u_int64_t> minimizers_pos;
		_minimizerParser->parse(rleSequence, minimizers, minimizers_pos);

		vector<KmerVec> kminmers; 
		vector<ReadKminmer> kminmersInfo;
		MDBG::getKminmers(_minimizerSize, _kminmerSize, minimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, readIndex, false);

		for(size_t i=0; i<kminmers.size(); i++){

			KmerVec& vec = kminmers[i];
			//if(_mdbg->_dbg_nodes.find(kminmers[i]) == _mdbg->_dbg_nodes.end()) continue;

			//u_int32_t nodeName = _mdbg->_dbg_nodes[kminmers[i]]._index;
			if(_mdbg->_dbg_nodes.find(vec) != _mdbg->_dbg_nodes.end()){

				u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
				//2262408
				//if("utg009434l" == string(read->name.s)){
				//	cout << nodeName << endl;
				//	_hifiasm_startingNodenames.insert(nodeName);
				//}
				//cout << _evaluation_hifiasmGroundTruth_nodeName_to_unitigName.size() << endl;

				_evaluation_hifiasmGroundTruth_nodeName_to_unitigName[_mdbg->_dbg_nodes[vec]._index].push_back(read._header);
				_evaluation_hifiasmGroundTruth_path.push_back(_mdbg->_dbg_nodes[vec]._index);

				//if(_evaluation_hifiasmGroundTruth_nodeNamePosition.find(_mdbg->_dbg_nodes[vec]._index) == _evaluation_hifiasmGroundTruth_nodeNamePosition.end()){
				//	_evaluation_hifiasmGroundTruth_nodeNamePosition[_mdbg->_dbg_nodes[vec]._index] = _extract_truth_kminmers_read_position;
				//}
				//cout << mdbg->_dbg_nodes[vec]._index << " " << sequence.getComment() << endl;
				
				//_file_groundTruth_hifiasm_position << nodeName << "," << _extract_truth_kminmers_read_position << endl;
			}
			else{
				//cout << "Not found position: " << position << endl;
				_evaluation_hifiasmGroundTruth_path.push_back(-1);
			}


			//if(_evaluation_hifiasmGroundTruth.find(vec) != _evaluation_hifiasmGroundTruth.end()){
				//cout << position << endl;
			//	continue;
			//}


			//_evaluation_hifiasmGroundTruth[vec] = datasetID;
			//_evaluation_hifiasmGroundTruth_position[vec] = _extract_truth_kminmers_read_position;
			//_extract_truth_kminmers_read_position += 1;

			
			//cout << _extract_truth_kminmers_read_position << endl;
		}


	}


	void extract_truth_kminmers(){

		_file_groundTruth_hifiasm_position.open(_inputDir + "/groundtruth_hifiasm_position.csv");
		_file_groundTruth_hifiasm_position << "Name,Position" << endl;

		_extract_truth_kminmers_read_position = 0;
		_minimizerParser = new MinimizerParser(_minimizerSize, _minimizerDensity);
		
		auto fp = std::bind(&ToBasespace::extract_truth_kminmers_read, this, std::placeholders::_1);
		ReadParser readParser(_truthInputFilename, true, false);
		readParser.parse(fp);

		_file_groundTruth_hifiasm_position.close();

		//delete _minimizerParser;
	}

};	


#endif 


