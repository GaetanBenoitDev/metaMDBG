

/*
Tester la vitesse de racon pour verifier que je suis sur la bonne voix
1) Compter les kminmer dans les contigs
2) Si un kminmers est répété, on veut récolté autant de miode
3) Recuperation des modeles: Indexer les readIndex, selectioner le readIndex qui suporte le plus de noeuds du contig
4) Utiliser edlib pour récupérer les 20 copies les plus proches
*/

#ifndef MDBG_METAG_TOBASESPACE2
#define MDBG_METAG_TOBASESPACE2

#include "../Commons.hpp"
//#include "../utils/edlib.h"
//#include "../utils/spoa/include/spoa/spoa.hpp"
#include "../utils/DnaBitset.hpp"

//#include "../utils/spoa64/include/spoa64/spoa.hpp"
//#include "../utils/abPOA2/include/abpoa.h"
//#include <seqan/align.h>
//#include <seqan/graph_msa.h>
//#include <cstring>

/*
struct KminmerSequenceVariant{
	u_int16_t _editDistance;
	DnaBitset* _sequence;
};

struct KminmerSequenceVariant_Comparator {
	bool operator()(KminmerSequenceVariant const& p1, KminmerSequenceVariant const& p2)
	{
		return p1._editDistance < p2._editDistance;
	}
};

typedef priority_queue<KminmerSequenceVariant, vector<KminmerSequenceVariant>, KminmerSequenceVariant_Comparator> VariantQueue;

struct ReadSequence{
	u_int64_t _readIndex;
	DnaBitset* _sequence;
	VariantQueue _variants;
};
*/
class ToBasespace2 : public Tool{
    
public:

	typedef phmap::parallel_flat_hash_map<KmerVec, vector<u_int32_t>, phmap::priv::hash_default_hash<KmerVec>, phmap::priv::hash_default_eq<KmerVec>, std::allocator<std::pair<KmerVec, vector<u_int32_t>>>, 4, std::mutex> KminmerReadMap;



	struct MinimizerPairPosition{

		ReadType _readIndex;
		u_int32_t _positionIndex;
	};

	typedef phmap::parallel_flat_hash_map<KmerVec, vector<MinimizerPairPosition>, phmap::priv::hash_default_hash<KmerVec>, phmap::priv::hash_default_eq<KmerVec>, std::allocator<std::pair<KmerVec, vector<MinimizerPairPosition>>>, 4, std::mutex> KminmerPosMap;

	

	string _inputFilename;
	string _inputFilenameContig;
	string _inputDir;
	float _minimizerDensity;
    size_t _minimizerSize;
    size_t _kminmerSize;
    size_t _kminmerSizeFirst;
	//bool _useHomopolymerCompression;
	//bool _isFirstPass;
	//bool _isOutputFasta;

	string _filename_outputContigs;
	string _filename_kminmerSequences;
	//MDBG* _mdbg;
	MinimizerParser* _minimizerParser;
	string _truthInputFilename;
	int _nbCores;
	long double _averageDistanceBetweenMinimizers;
	bool _skipCorrection;

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

	phmap::parallel_flat_hash_map<u_int32_t, vector<DnaBitset*>> _kminmerSequenceCopies_all_entire;
	phmap::parallel_flat_hash_map<u_int32_t, vector<DnaBitset*>> _kminmerSequenceCopies_all_left;
	phmap::parallel_flat_hash_map<u_int32_t, vector<DnaBitset*>> _kminmerSequenceCopies_all_right;
	
	//unordered_set<u_int32_t> isKminmerRepeated;

	phmap::parallel_flat_hash_map<ReadNodeName, DnaBitset*> _repeatedKminmerSequence_entire;
	phmap::parallel_flat_hash_map<ReadNodeName, DnaBitset*> _repeatedKminmerSequence_left;
	phmap::parallel_flat_hash_map<ReadNodeName, DnaBitset*> _repeatedKminmerSequence_right;

	//unordered_map<ReadNodeName, VariantQueue> _repeatedKminmerSequence_copies_entire;
	//unordered_map<ReadNodeName, VariantQueue> _repeatedKminmerSequence_copies_left;
	//unordered_map<ReadNodeName, VariantQueue> _repeatedKminmerSequence_copies_right;

	phmap::parallel_flat_hash_map<u_int32_t, DnaBitset*> _kminmerSequence_entire;
	phmap::parallel_flat_hash_map<u_int32_t, DnaBitset*> _kminmerSequence_left;
	phmap::parallel_flat_hash_map<u_int32_t, DnaBitset*> _kminmerSequence_right;
	phmap::parallel_flat_hash_map<u_int32_t, vector<ReadSequence>> _kminmerSequence_entire_multi;
	phmap::parallel_flat_hash_map<u_int32_t, vector<ReadSequence>> _kminmerSequence_left_multi;
	phmap::parallel_flat_hash_map<u_int32_t, vector<ReadSequence>> _kminmerSequence_right_multi;
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

	ToBasespace2(): Tool (){

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

		args::ArgumentParser parser("toBasespace", ""); //"This is a test program.", "This goes after the options."
		args::Positional<std::string> arg_outputDir(parser, "outputDir", "", args::Options::Required);
		args::Positional<std::string> arg_inputContigFilename(parser, "inputContigFilename", "", args::Options::Required);
		args::Positional<std::string> arg_outputContigFilename(parser, "outputContigFilename", "", args::Options::Required);
		args::Positional<std::string> arg_inputReadFilename(parser, "inputReadFilename", "", args::Options::Required);
		args::Flag arg_skipCorrection(parser, "", "Skip read correction", {"skip-correction"});
		//args::Positional<std::string> arg_contigs(parser, "contigs", "", args::Options::Required);
		//args::PositionalList<std::string> arg_readFilenames(parser, "reads", "Input filename(s) (separated by space)", args::Options::Required);
		//args::ValueFlag<int> arg_l(parser, "", "Minimizer length", {ARG_MINIMIZER_LENGTH2}, 13);
		//args::ValueFlag<float> arg_d(parser, "", "Minimizer density", {ARG_MINIMIZER_DENSITY2}, 0.005f);
		//args::ValueFlag<std::string> arg_contigs(parser, "", "", {ARG_INPUT_FILENAME_CONTIG}, "");
		//args::Flag arg_homopolymerCompression(parser, "", "Activate homopolymer compression", {ARG_HOMOPOLYMER_COMPRESSION});
		args::ValueFlag<int> arg_nbCores(parser, "", "Number of cores", {ARG_NB_CORES2}, NB_CORES_DEFAULT_INT);
		//args::Flag arg_firstPass(parser, "", "Is first pass of multi-k", {ARG_FIRST_PASS});
		//args::Flag arg_isFinalAssembly(parser, "", "Is final multi-k pass", {ARG_FINAL});
		//args::Flag arg_firstPass(parser, "", "Is first pass of multi-k", {ARG_FIRST_PASS});
		args::Flag arg_help(parser, "", "", {'h', "help"}, args::Options::Hidden);
		//args::HelpFlag help(parser, "help", "Display this help menu", {'h'});
		//args::CompletionFlag completion(parser, {"complete"});

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

		_inputDir = args::get(arg_outputDir);
		_inputFilenameContig = args::get(arg_inputContigFilename);
		_filename_outputContigs = args::get(arg_outputContigFilename);
		_inputFilename = args::get(arg_inputReadFilename);
		_nbCores = args::get(arg_nbCores);

		//_useHomopolymerCompression = false;
		//if(arg_homopolymerCompression){
		//	_useHomopolymerCompression = true;
		//}

		/*
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
			_logFile << options.help() << endl;
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
			std::_logFile << options.help() << std::endl;
			std::cerr << e.what() << std::endl;
			std::exit(EXIT_FAILURE);
		}
		*/
		/*
		_inputFilename = getInput()->getStr(STR_INPUT);
		_inputDir = getInput()->getStr(STR_INPUT_DIR);
		*/

		_skipCorrection = false;
		if(arg_skipCorrection){
			_skipCorrection = true;
		}

		string filename_parameters = _inputDir + "/parameters.gz";
		gzFile file_parameters = gzopen(filename_parameters.c_str(),"rb");
		gzread(file_parameters, (char*)&_minimizerSize, sizeof(_minimizerSize));
		gzread(file_parameters, (char*)&_kminmerSize, sizeof(_kminmerSize));
		gzread(file_parameters, (char*)&_minimizerDensity, sizeof(_minimizerDensity));
		gzread(file_parameters, (char*)&_kminmerSizeFirst, sizeof(_kminmerSizeFirst));
		gzclose(file_parameters);

		_kminmerSize = 2;
		_averageDistanceBetweenMinimizers = 1.0 / _minimizerDensity;
		
		openLogFile(_inputDir);

		Logger::get().debug() << "";
		Logger::get().debug() << "Input dir: " << _inputDir;
		//_logFile << "Output filename: " << _outputFilename << endl;
		Logger::get().debug() << "Minimizer length: " << _minimizerSize;
		Logger::get().debug() << "Kminmer length: " << _kminmerSize;
		Logger::get().debug() << "Density: " << _minimizerDensity;
		Logger::get().debug() << "";

		//_filename_outputContigs = _inputDir + "/contigs.fasta.gz";
		//_filename_outputContigs = _inputFilenameContig + ".fasta.gz"; //_inputDir + "/tmpContigs.fasta.gz";
		_minimizerParser = new MinimizerParser(_minimizerSize, _minimizerDensity);
	}

	vector<vector<u_int32_t>> _unitigDatas;
	ofstream _contigFileSupported;
	ifstream _contigFileSupported_input;



    void execute (){
		

		_checksum = 0;


		//_kminmerSize = 4;

		//_logFile << "Loading mdbg" << endl;
		//string mdbg_filename = _inputDir + "/kminmerData_min_init.txt";
		//_mdbg = new MDBG(_kminmerSize);
		//_mdbg->load(mdbg_filename, false);
		//_logFile << "MDBG nodes: " << _mdbg->_dbg_nodes.size() << endl;

		Logger::get().debug() << "Indexing contigs";
		indexContigs();

		Logger::get().debug() << "Loading contigs";
		loadContigs();

		Logger::get().debug() << "Aligning reads";
		alignReads();

		//exit(1);
		//indexContigs();
		//indexKminmers();
		//indexReads();
		setupKminmerSequenceExtraction();

		//_logFile << "Loading original mdbg" << endl;
		//string mdbg_filename = _inputDir + "/mdbg_nodes_init.gz";
		//_mdbgInit = new MDBG(_originalKminmerSize);
		//_mdbgInit->load(mdbg_filename);
		//_logFile << "MDBG nodes: " << _mdbgInit->_dbg_nodes.size() << endl;

		//string inputFilenameContigSmall = _inputDir + "/small_contigs.bin";

		//loadContigs_min(_inputFilenameContig);
		//loadContigs_min(inputFilenameContigSmall);
		//collectBestSupportingReads(_inputFilenameContig);
		//collectBestSupportingReads(inputFilenameContigSmall);
		//extractKminmerSequences_all();
		_kminmerCounts.clear();
		_unitigDatas.clear();

		//isKminmerRepeated.clear();
		//_logFile << isKminmerRepeated.size() << endl;
		//getchar();

		extractKminmerSequences_model();

		//cout << _readPosition_entire.size() << " " << _readPosition_left.size() << " " << _readPosition_right.size() << endl;
		//extractKminmerSequences_allVariants();

		//_logFile << _nodeName_all_left.size() << " " << _nodeName_all_right.size() << endl;

		//_logFile << "Correcting kminmer sequences" << endl;
		//exit(1);

		//auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -3);  // linear gaps
		//spoa::Graph graph{};
		//_alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -3);  // linear gaps
		//spoa::Graph graph{};
		
		/*
		endCorrection(_kminmerSequenceCopies_all_entire, _isDoneNodeName_entire, _kminmerSequence_entire, _repeatedKminmerSequence_entire, _kminmerSequence_entire_multi);
		endCorrection(_kminmerSequenceCopies_all_left, _isDoneNodeName_left, _kminmerSequence_left, _repeatedKminmerSequence_left, _kminmerSequence_left_multi);
		endCorrection(_kminmerSequenceCopies_all_right, _isDoneNodeName_right, _kminmerSequence_right, _repeatedKminmerSequence_right, _kminmerSequence_right_multi);
		*/
	
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
		
		_contigIndex = 0;
		_basespaceContigFile = gzopen(_filename_outputContigs.c_str(),"wb");

		createBaseContigs(_inputFilenameContig);
		//createBaseContigs(inputFilenameContigSmall);

		gzclose(_basespaceContigFile);
		//delete _mdbg;
		/*
		gzclose(_outputFile_left);
		gzclose(_outputFile_right);


		return;

		//loadContigs(_inputDir + "/minimizer_contigs_complete.gz");
		//loadContigs(_inputDir + "/eval/composition/22/debug_longUnitigs.gz");
		//loadContigs(_inputDir + "/eval/composition//3/debug_longUnitigsNeighbors.gz");

		_logFile << "TODO: remove ContigNode and _requiredCopiers_entire, instead use unorderedmap nodeName => vector<(ReadIndex, string)>" << endl;




		extractKminmerSequences();
		lalalala();
		delete _mdbg;

		//exit(1);
		//loadKminmerSequences();
		createBaseContigs(_inputDir + "/minimizer_contigs.gz", _filename_outputContigs.c_str());
		//createBaseContigs(_inputDir + "/minimizer_contigs_complete.gz", _filename_outputContigs.c_str());
		//createBaseContigs(_inputDir + "/eval/composition//22//debug_longUnitigs.gz", _inputDir + "/eval/composition//22//debug_longUnitigs.fasta.gz");
		//createBaseContigs(_inputDir + "/eval/composition//3//debug_longUnitigsNeighbors.gz", _inputDir + "/eval/composition//3/debug_longUnitigsNeighbors.fasta.gz");

		_logFile << endl << "Contig filename: " << _filename_outputContigs << endl;
		*/

		//closeLogFile();
	}

	struct KminmerPos{
		u_int32_t _contigIndex;
		u_int32_t _pos;
	};

	struct ReadPositionMatch{
		u_int32_t _readIndex;
		u_int32_t _readPosition;
		u_int32_t _matchLength;
	};

	u_int32_t _kminmerNodeName;
	phmap::parallel_flat_hash_map<KmerVec, u_int32_t> _kmerVec_to_nodeName;
	phmap::parallel_flat_hash_map<KmerVec, vector<KminmerPos>> _kmerVecIndex;
	phmap::parallel_flat_hash_map<u_int32_t, KminmerPos> _bestReadMatch;
	phmap::parallel_flat_hash_map<ContigPosition, ReadPositionMatch> _contigPosition_to_readPosition;
	phmap::parallel_flat_hash_map<ReadPosition, DnaBitset*> _readPosition_entire;
	phmap::parallel_flat_hash_map<ReadPosition, DnaBitset*> _readPosition_left;
	phmap::parallel_flat_hash_map<ReadPosition, DnaBitset*> _readPosition_right;

	KminmerPosMap _kminmer_to_readIndex;

	void indexContigs(){

		KminmerParserParallel parser(_inputFilenameContig, _minimizerSize, _kminmerSize, false, false, _nbCores);
		parser.parse(IndexContigsFunctor(*this));
	}

	class IndexContigsFunctor {

		public:

		ToBasespace2& _parent;

		IndexContigsFunctor(ToBasespace2& parent) : _parent(parent){
		}

		IndexContigsFunctor(const IndexContigsFunctor& copy) : _parent(copy._parent){
		}

		~IndexContigsFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {


			u_int32_t readIndex = kminmerList._readIndex;
			if(readIndex % 100000 == 0) Logger::get().debug() << readIndex;


			for(u_int32_t i=0; i<kminmerList._kminmersInfo.size(); i++){
			
				
				const ReadKminmerComplete& kminmerInfo = kminmerList._kminmersInfo[i];

				const KmerVec& vec = kminmerInfo._vec;
				
				u_int32_t positionIndex = i;

				_parent._kminmer_to_readIndex.lazy_emplace_l(vec,
				[&readIndex, &positionIndex](KminmerPosMap::value_type& v) { // key exist
					v.second.push_back({readIndex, positionIndex});
				},           
				[&vec, &readIndex, &positionIndex](const KminmerPosMap::constructor& ctor) { // key inserted
					
					//vector<u_int32_t> readIndexes = {2}; //inital count of this kminmer
					vector<MinimizerPairPosition> readIndexes = {{readIndex, positionIndex}}; //inital count of this kminmer

					ctor(vec, readIndexes); 

				}); // construct value_type in place when key not present

			}
			
			
		}
		
	};


	//struct MinimizerRead{
	//	u_int64_t _readIndex;
	//	vector<MinimizerType> _minimizers;
	//	vector<u_int8_t> _qualities;
	//	vector<u_int8_t> _readMinimizerDirections;
	//};

	struct MinimizerContig{
		u_int32_t _contigIndex;
		vector<MinimizerType> _contigMinimizers;
	};

	vector<MinimizerContig> _mContigs;

	void loadContigs(){
		KminmerParserParallel parser(_inputFilenameContig, _minimizerSize, _kminmerSize, false, false, 1);
		parser.parseSequences(LoadContigsFunctor(*this));
	}

	class LoadContigsFunctor {

		public:

		ToBasespace2& _parent;

		LoadContigsFunctor(ToBasespace2& parent) : _parent(parent){
		}

		LoadContigsFunctor(const LoadContigsFunctor& copy) : _parent(copy._parent){
		}

		~LoadContigsFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {

			u_int32_t readIndex = kminmerList._readIndex;
			if(readIndex % 100000 == 0) Logger::get().debug() << readIndex;

			_parent._mContigs.push_back({readIndex, kminmerList._readMinimizers});
			

		}
	};


	void alignReads(){

		//bool hasQuality = true;
		//if(_useHomopolymerCompression){ //hifi
		//	hasQuality = false;
		//}

		if(_skipCorrection){
			KminmerParserParallel parser(_inputDir + "/read_data_init.txt", _minimizerSize, _kminmerSize, false, false, _nbCores);
			parser.parseSequences(AlignReadFunctor(*this));
		}
		else{
			KminmerParserParallel parser(_inputDir + "/read_data_init.txt", _minimizerSize, _kminmerSize, false, true, _nbCores);
			parser._densityThreshold = _minimizerDensity;
			parser.parseSequences(AlignReadFunctor(*this));
		}

	}

	
	class AlignReadFunctor {

		public:

		ToBasespace2& _parent;
		MinimizerAligner* _minimizerAligner;
		//std::unique_ptr<spoa64::AlignmentEngine> _alignmentEngine;// = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -4);
		//std::unique_ptr<spoa64::AlignmentEngine> _alignmentEngineTrimming;

		AlignReadFunctor(ToBasespace2& parent) : _parent(parent){
			_minimizerAligner = nullptr;
		}

		AlignReadFunctor(const AlignReadFunctor& copy) : _parent(copy._parent){
			//_alignmentEngine = spoa64::AlignmentEngine::Create(spoa64::AlignmentType::kNW, 3, -5, -4);
			//_alignmentEngineTrimming = spoa64::AlignmentEngine::Create(spoa64::AlignmentType::kSW, 3, -5, -4);
			//_alignmentEngine = spoa64::AlignmentEngine::Create(spoa64::AlignmentType::kSW, 3, -2, -1);
			_minimizerAligner = new MinimizerAligner(3, -1, -1, -1);
		}

		~AlignReadFunctor(){
			if(_minimizerAligner != nullptr) delete _minimizerAligner;
		}


		void operator () (const KminmerList& kminmerList) {
			_parent.alignRead(kminmerList, _minimizerAligner);
		}

	};



	void alignRead(const KminmerList& kminmerList, MinimizerAligner* minimizerAligner){



		u_int64_t readIndex = kminmerList._readIndex;

		if(readIndex % 10000 == 0) Logger::get().debug() << readIndex;
		//cout << readIndex << endl;

		vector<MinimizerType> readMinimizers = kminmerList._readMinimizers;


		const MinimizerRead queryRead = {kminmerList._readIndex, kminmerList._readMinimizers, kminmerList._minimizerPos, kminmerList._readQualities, kminmerList._readMinimizerDirections};
		
		unordered_map<MinimizerType, vector<u_int32_t>> queryReadMinimizerPositionMap;

		
		for(u_int32_t i=0; i<queryRead._minimizers.size(); i++){

			MinimizerType minimizer = queryRead._minimizers[i];
			//u_int32_t position = read._positions[i];
			//bool isReversed = read._readMinimizerDirections[i];

			queryReadMinimizerPositionMap[minimizer].push_back(i);
		}

		/*
		vector<MinimizerType> readMinimizers;// = kminmerList._readMinimizers;

		static MinimizerType maxHashValue = -1;
		static MinimizerType minimizerBound = 0.005 * maxHashValue;


		for(size_t i=0; i<kminmerList._readMinimizers.size(); i++){
			u_int64_t m = kminmerList._readMinimizers[i];
			if(m > minimizerBound) continue;

			readMinimizers.push_back(kminmerList._readMinimizers[i]);
			//minimizerPosDummy.push_back(minimizers[i]);

		}
		*/

		unordered_set<MinimizerType> minimizerSet;

		for(size_t i=0; i<readMinimizers.size(); i++){
			MinimizerType m = readMinimizers[i];

			minimizerSet.insert(m);
		}

		//Kminmers !!
		vector<ReadKminmerComplete> kminmersInfos;
		vector<u_int32_t> minimizersPos(readMinimizers.size(), 0);
		vector<u_int8_t> minimizerQualities(readMinimizers.size(), 0);
		MDBG::getKminmers_complete(_kminmerSize, readMinimizers, minimizersPos, kminmersInfos, readIndex, minimizerQualities);
		
		/*
		unordered_map<u_int32_t, u_int32_t> readIndex_to_matchCount;
		vector<KmerVec> readKminmers;

		for(const ReadKminmerComplete& kminmerInfo : kminmersInfos){

			const KmerVec& vec = kminmerInfo._vec;
			readKminmers.push_back(vec);

			if(_kminmer_to_readIndex.find(vec) == _kminmer_to_readIndex.end()) continue;

			for(const MinimizerPairPosition& mPos : _kminmer_to_readIndex[vec]){
				if(readIndex == mPos._readIndex) continue; //Currently corrected read
				readIndex_to_matchCount[mPos._readIndex] += 1;
			}

		}
		*/
		unordered_map<ReadType, ReadMatchBound> referenceReadIndex_to_anchors2;

		for(size_t i=0; i<kminmersInfos.size(); i++){

			const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];
			const KmerVec& vec = kminmerInfo._vec;

			if(_kminmer_to_readIndex.find(vec) == _kminmer_to_readIndex.end()) continue;

			for(const MinimizerPairPosition& contigPosition : _kminmer_to_readIndex[vec]){
				//if(read._readIndex == queryReadPosition._readIndex) continue;

				if(referenceReadIndex_to_anchors2.find(contigPosition._readIndex) == referenceReadIndex_to_anchors2.end()){
					referenceReadIndex_to_anchors2[contigPosition._readIndex] = ReadMatchBound();
				}

				ReadMatchBound& bound = referenceReadIndex_to_anchors2[contigPosition._readIndex];

				bound._nbMatches += 1;

				if(contigPosition._positionIndex < bound._minIndex){
					bound._minIndex = contigPosition._positionIndex;
				}
				if(contigPosition._positionIndex > bound._maxIndex){
					bound._maxIndex = contigPosition._positionIndex;
				}

			}

		}

		//spoa64::Graph graph{};
		

		//vector<u_int64_t> weights(readMinimizers.size(), 1);
		//graph.AddAlignment(spoa64::Alignment(), Utils::vec32_to_vec64(readMinimizers), readMinimizers.size(), weights);
		


		vector<MinimizerContig> contigSequences;

		for(const auto& it : referenceReadIndex_to_anchors2){

			ReadType contigIndex = it.first;
			const ReadMatchBound& readMatchBound = it.second;

			if(readMatchBound._nbMatches < 2) continue;


			const MinimizerContig& mContig = _mContigs[contigIndex];
			//u_int32_t contigIndex = mContig._contigIndex;

			const vector<MinimizerType>& contigMinimizers = mContig._contigMinimizers;

			int64_t startQueryPositionIndex = getStartQueryPositionIndex(contigMinimizers, readMatchBound._minIndex, 10000, queryReadMinimizerPositionMap);
			int64_t endQueryPositionIndex = getEndQueryPositionIndex(contigMinimizers, readMatchBound._maxIndex, 10000, queryReadMinimizerPositionMap);

			int64_t overHangOffset = 10000*_minimizerDensity;
			startQueryPositionIndex = max((int64_t)0, startQueryPositionIndex-overHangOffset);
			endQueryPositionIndex = min((int64_t)contigMinimizers.size(), endQueryPositionIndex+overHangOffset);

			//startQueryPositionIndex = 0;
			//endQueryPositionIndex = contigMinimizers.size();

			const vector<MinimizerType> subContigMinimizers_forward(contigMinimizers.begin() + startQueryPositionIndex, contigMinimizers.begin() + endQueryPositionIndex);
			vector<MinimizerType> subContigMinimizers_reverse = subContigMinimizers_forward;
			std::reverse(subContigMinimizers_reverse.begin(), subContigMinimizers_reverse.end());

			//vector<MinimizerType> minimizerSequence_forward = mContig._contigMinimizers;
			//vector<MinimizerType> minimizerSequence_reverse = mContig._contigMinimizers;
			//std::reverse(minimizerSequence_reverse.begin(), minimizerSequence_reverse.end());

			//vector<u_int8_t> qualities(minimizerSequence_forward.size(), 1);

			//if(_print_debug){
			//	cout << endl << "\tAlign sequence: " << i << " " << minimizerSequence.size() << " " << graph.nodes_.size() << endl; //<< " " << graphPOA->_graph->_nodes.size() << endl;
			//	for(MinimizerType m : minimizerSequence){
			//		cout << "\t" << m << endl; 
			//	}
			
			//}

			int64_t nbMatches_forward = 0;
			int64_t alignScore_forward = 0;
			MinimizerAligner::Alignment alignment_forward;
			getAlignmentScore(readMinimizers, subContigMinimizers_forward, minimizerAligner, alignment_forward, nbMatches_forward, alignScore_forward);
			
			int64_t nbMatches_reverse = 0;
			int64_t alignScore_reverse = 0;
			MinimizerAligner::Alignment alignment_reverse;
			getAlignmentScore(readMinimizers, subContigMinimizers_reverse, minimizerAligner, alignment_reverse, nbMatches_reverse, alignScore_reverse);

			if(nbMatches_forward < 2 && nbMatches_reverse < 2) continue;

			if(nbMatches_forward > nbMatches_reverse){
				extractAlignment(alignment_forward, readIndex, readMinimizers, contigIndex, subContigMinimizers_forward, alignScore_forward, false, startQueryPositionIndex);
			}
			else{
				extractAlignment(alignment_reverse, readIndex, readMinimizers, contigIndex, subContigMinimizers_reverse, alignScore_reverse, true, startQueryPositionIndex);
			}
			//cout << nbMatches_forward << " " << nbMatches_reverse << endl;

			//double alignmentScore = nbMatches / (double) readMinimizers.size();// / (double) minimizerSequence.size();
			
			//if(_print_debug){
			//	cout << "\tAl score: " << alignmentScore << endl; //<< " " << getAlignmentScore(minimizerSequence, readMinimizers, al) << endl;
			//}

			//if(nbMatches < 5) continue;
			//if(alignmentScore < 0.2) continue;
			//if(alignmentScore < 5) continue;

			//if(isReversed){
			//	std::reverse(qualities.begin(), qualities.end());
			//}

			
			//mappedReads.push_back({mOverlap._readIndex, minimizerSequence, qualities});
		}

	}

	int64_t getStartQueryPositionIndex(const vector<MinimizerType>& contigMinimizers, int64_t startIndex, int64_t maxGapLength, unordered_map<MinimizerType, vector<u_int32_t>>& queryReadMinimizerPositionMap){
		
		int64_t lastMacthingIndex = startIndex;
		int64_t i = startIndex;
		int64_t currentGapLength = 0;
		int64_t prevPosition = i * _averageDistanceBetweenMinimizers; //queryRead._minimizersPos[i];
		i -= 1;

		while(true){
			if(i < 0) break;
			if(currentGapLength > maxGapLength) break;

			const MinimizerType currentMinimizer = contigMinimizers[i];

			if(queryReadMinimizerPositionMap.find(currentMinimizer) == queryReadMinimizerPositionMap.end()){
				currentGapLength += (prevPosition-(i*_averageDistanceBetweenMinimizers));
			}
			else{
				currentGapLength = 0;
				lastMacthingIndex = i;
			}


			prevPosition = i*_averageDistanceBetweenMinimizers;
			i -= 1;
		}

		return lastMacthingIndex;
	}

	int64_t getEndQueryPositionIndex(const vector<MinimizerType>& contigMinimizers, int64_t startIndex, int64_t maxGapLength, unordered_map<MinimizerType, vector<u_int32_t>>& queryReadMinimizerPositionMap){

		int64_t lastMacthingIndex = startIndex;
		int64_t i = startIndex;
		int64_t currentGapLength = 0;
		int64_t prevPosition = i * _averageDistanceBetweenMinimizers; //queryRead._minimizersPos[i];
		i += 1;

		while(true){
			if(i >= contigMinimizers.size()) break;
			if(currentGapLength > maxGapLength) break;

			const MinimizerType currentMinimizer = contigMinimizers[i];

			if(queryReadMinimizerPositionMap.find(currentMinimizer) == queryReadMinimizerPositionMap.end()){
				currentGapLength += ((i * _averageDistanceBetweenMinimizers)-prevPosition);
			}
			else{
				currentGapLength = 0;
				lastMacthingIndex = i;
			}


			prevPosition = (i * _averageDistanceBetweenMinimizers);
			i += 1;
		}

		return lastMacthingIndex;
	}

	void getAlignmentScore(const vector<MinimizerType>& s1, const vector<MinimizerType>& s2, MinimizerAligner* minimizerAligner, MinimizerAligner::Alignment& alignmentResult, int64_t& nbMatches, int64_t& alignScore){
		

		MinimizerAligner::Alignment alignment = minimizerAligner->performAlignment(s1, s2); //al->Align(Utils::vec32_to_vec64(s2), s2.size(), graph);
		alignmentResult = alignment;

		nbMatches = 0;
		alignScore = 0;


		
		for (size_t i=0; i<alignment.size(); i++) {

			int64_t v1 = alignment[i].first;
			int64_t v2 = alignment[i].second;
			
			if(v1 == -1){ //insert in
				//alignScore -= 1;
			}
			else if(v2 == -1){ //insert in
				//alignScore -= 1;
			}
			else if(s1[v1] == s2[v2]){
				alignScore += 1;
				nbMatches += 1;
			}
			else{
				alignScore -= 1;
			}
		}


	}

	void extractAlignment(const MinimizerAligner::Alignment& alignment, u_int32_t readIndex, const vector<MinimizerType>& readMinimizers, u_int32_t contigIndex, const vector<MinimizerType>& contigMinimizers, u_int32_t nbMatches, bool isReversed, const int64_t startQueryPositionIndex){
		
		int64_t prevReadPosition = -1;
		int64_t prevContigPosition = -1;
		
		for (size_t i=0; i<alignment.size()-1; i++) {

			int64_t readPosition = alignment[i].first;
			int64_t contigPosition = alignment[i].second;
			
			if(readPosition == -1){ //insert in
				prevReadPosition = -1;
				prevContigPosition = -1;
			}
			else if(contigPosition == -1){ //insert in
				prevReadPosition = -1;
				prevContigPosition = -1;
			}
			else if(readMinimizers[readPosition] == contigMinimizers[contigPosition]){

				if(prevReadPosition != -1 && prevContigPosition != -1){
					u_int32_t prevContigPositionNorm = prevContigPosition;
					if(isReversed){
						prevContigPositionNorm = contigMinimizers.size() - contigPosition - 1;
						//return;
					}
					prevContigPositionNorm += startQueryPositionIndex;
					indexMatchingPosition(readMinimizers[prevReadPosition], readMinimizers[readPosition], {readIndex, (u_int32_t)prevReadPosition, nbMatches}, {contigIndex, prevContigPositionNorm});
				}

				prevReadPosition = readPosition;
				prevContigPosition = contigPosition;
			}
			else{
				prevReadPosition = -1;
				prevContigPosition = -1;
			}

		}

	}

	void indexMatchingPosition(MinimizerType minimizer1, MinimizerType minimizer2, const ReadPositionMatch& readPosition, const ContigPosition& contigPosition){


		
		#pragma omp critical(indexMatchingPosition)
		{
			if(_contigPosition_to_readPosition.find(contigPosition) == _contigPosition_to_readPosition.end()){
				_contigPosition_to_readPosition[contigPosition] = readPosition;
			}
			else{
				
				const ReadPositionMatch& readPositionExisting = _contigPosition_to_readPosition[contigPosition];


				if(readPosition._matchLength == readPositionExisting._matchLength){
					if(readPosition._readIndex < readPositionExisting._readIndex){
						_contigPosition_to_readPosition[contigPosition] = readPosition;
					}

				}
				else if(readPosition._matchLength > readPositionExisting._matchLength){
					_contigPosition_to_readPosition[contigPosition] = readPosition;

				}
			}


		}
		

	}


	/*
	struct MinimizerAlignment{
		vector<MinimizerType> _minimizerSequence;
		vector<MinimizerType> _minimizerQualities;
		MinimizerAligner::Alignment _alignment;
		double _alignmentScore;
		int _lastMatchPos;
	};

	void indexReads(){


		Logger::get().debug() << "Indexing reads";
		_unitigDatas.resize(_kmerVec_to_nodeName.size());
		
		KminmerParserParallel parser(_inputDir + "/read_data_init.txt", _minimizerSize, _kminmerSize, false, true, _nbCores);
		//parser2.parse(FilterKminmerFunctor2(*this));
		parser.parse(IndexReadsFunctor(*this));

		//KminmerParser parser(_inputDir + "/read_data_init.txt", _minimizerSize, _kminmerSize, false, true);
		//auto fp = std::bind(&ToBasespace::indexReads_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
		//parser.parse(fp);
	}

	struct MatchPosition{
		u_int32_t _contigIndex;
		u_int32_t _contigPos;
		u_int32_t _contigMatchPos_right;
		u_int32_t _contigMatchPos_left;
		u_int32_t _matchLength;
		bool _isValid_right;
		bool _isValid_left;
	};

	class IndexReadsFunctor {

		public:

		ToBasespace2& _parent;

		IndexReadsFunctor(ToBasespace2& parent) : _parent(parent){
		}

		IndexReadsFunctor(const IndexReadsFunctor& copy) : _parent(copy._parent){
		}

		~IndexReadsFunctor(){
		}

		void operator () (KminmerList& kminmerList) {

			processKminmers(kminmerList, false);

			std::reverse(kminmerList._kminmersInfo.begin(), kminmerList._kminmersInfo.end());

			processKminmers(kminmerList, true);
		}

			
		void processKminmers(const KminmerList& kminmerList, bool isReversed){

			for(u_int32_t i=0; i<kminmerList._kminmersInfo.size(); i++){
			
				
				//cout << "------- " << i << endl;

				//_logFile << readIndex << " " << i << endl;
				const ReadKminmerComplete& kminmerInfo = kminmerList._kminmersInfo[i];

				const KmerVec& vec = kminmerInfo._vec;
				
				if(_parent._kmerVecIndex.find(vec) == _parent._kmerVecIndex.end()) continue;

				vector<MatchPosition> matches;

				const vector<KminmerPos>& initialMatches = _parent._kmerVecIndex[vec];
				for(const KminmerPos& initialMatch : initialMatches){
					matches.push_back({initialMatch._contigIndex, initialMatch._pos, initialMatch._pos, initialMatch._pos, 1, true, true});
				}

				findBestMatches(kminmerList, matches, i);

				for(const MatchPosition& match : matches){

					ContigPosition contigPosition =  {match._contigIndex, match._contigPos};
					ReadPositionMatch readPosition = {kminmerList._readIndex, i, match._matchLength};

					if(isReversed){
						readPosition._readPosition = kminmerList._kminmersInfo.size()-1-readPosition._readPosition;
					}

					#pragma omp critical(processKminmers)
					{
						if(_parent._contigPosition_to_readPosition.find(contigPosition) == _parent._contigPosition_to_readPosition.end()){
							_parent._contigPosition_to_readPosition[contigPosition] = readPosition;
							//cout << isReversed << " " << readPosition._matchLength << endl;
						}
						else{
							
							const ReadPositionMatch& readPositionExisting = _parent._contigPosition_to_readPosition[contigPosition];

							if(readPosition._matchLength > readPositionExisting._matchLength){
								_parent._contigPosition_to_readPosition[contigPosition] = readPosition;
								//cout << isReversed << " " << readPosition._matchLength << endl;

								//if(contigPosition._contigIndex == 6 && contigPosition._contigPosition == 2401)
								//cout << contigPosition._contigIndex << " " << contigPosition._contigPosition << " " << readPosition._matchLength << endl;
							}
						}


						//cout << match._contigIndex << " " << match._contigPos << " " << match._matchLength << endl;
						//getchar();
					}
				}

			}

		}

		void findBestMatches(const KminmerList& kminmerList, vector<MatchPosition>& initialMatches, long readStartPos){

			for(u_int32_t i=readStartPos+1; i<kminmerList._kminmersInfo.size(); i++){
			
				//cout << i << " " << kminmerList._kminmersInfo.size() << endl;
				//cout << "\t" << i << endl;
				
				//_logFile << readIndex << " " << i << endl;
				const ReadKminmerComplete& kminmerInfo = kminmerList._kminmersInfo[i];

				const KmerVec& vec = kminmerInfo._vec;
				if(_parent._kmerVecIndex.find(vec) == _parent._kmerVecIndex.end()){
					//for(const MatchPosition& initialMatch : initialMatches){
					//	initialMatch._isValid_left = false;
					//	break;
					//}
					break;
				} 

				const vector<KminmerPos>& currentMatches = _parent._kmerVecIndex[vec];
				
				//if(_parent._kmerVecIndex.find(vec) == _parent._kmerVecIndex.end()) continue;

				for(MatchPosition& initialMatch : initialMatches){

					//cout << "\t" <<  initialMatch._contigIndex << " " << initialMatch._contigPos << " " << initialMatch._contigMatchPos << endl;

					if(!initialMatch._isValid_right) continue;

					bool isValid = false;

					for(const KminmerPos& currentMatch : currentMatches){


						if(initialMatch._contigIndex == currentMatch._contigIndex){
							if(initialMatch._contigMatchPos_right+1 == currentMatch._pos){
								initialMatch._contigMatchPos_right += 1;
								initialMatch._matchLength += 1;
								isValid = true;
								break;
							}
						}
					}

					initialMatch._isValid_right = isValid;
				}

				bool allMatchInvalid = true;
				for(const MatchPosition& initialMatch : initialMatches){
					if(initialMatch._isValid_right){
						allMatchInvalid = false;
						break;
					}
				}

				if(allMatchInvalid) break;



			}

			

			for(long i=readStartPos-1; i >= 0; i--){
			
				//cout << "\t" << i << endl;
				
				//_logFile << readIndex << " " << i << endl;
				const ReadKminmerComplete& kminmerInfo = kminmerList._kminmersInfo[i];

				const KmerVec& vec = kminmerInfo._vec;

				if(_parent._kmerVecIndex.find(vec) == _parent._kmerVecIndex.end()){
					//for(MatchPosition& initialMatch : initialMatches){
					//	initialMatch._isValid_left = false;
					//	break;
					//}
					break;
				} 

				const vector<KminmerPos>& currentMatches = _parent._kmerVecIndex[vec];
				
				//if(_parent._kmerVecIndex.find(vec) == _parent._kmerVecIndex.end()) continue;

				for(MatchPosition& initialMatch : initialMatches){

					//cout << "\t" <<  initialMatch._contigIndex << " " << initialMatch._contigPos << " " << initialMatch._contigMatchPos << endl;

					if(!initialMatch._isValid_left) continue;

					bool isValid = false;

					for(const KminmerPos& currentMatch : currentMatches){


						if(initialMatch._contigIndex == currentMatch._contigIndex){
							if(initialMatch._contigMatchPos_left-1 == currentMatch._pos){
								initialMatch._contigMatchPos_left -= 1;
								initialMatch._matchLength += 1;
								isValid = true;
								break;
							}
						}
					}

					initialMatch._isValid_left = isValid;
				}

				bool allMatchInvalid = true;
				for(const MatchPosition& initialMatch : initialMatches){
					if(initialMatch._isValid_left){
						allMatchInvalid = false;
						break;
					}
				}

				if(allMatchInvalid) break;


			}


		}

	};
	*/


	void setupKminmerSequenceExtraction(){


		Logger::get().debug() << "Parsing contigs";
		KminmerParserParallel parser(_inputFilenameContig, _minimizerSize, _kminmerSize, false, false, 1);
		parser.parse(SetupKminmerSequenceExtractionFunctor(*this));

	}


	class SetupKminmerSequenceExtractionFunctor {

		public:

		ToBasespace2& _toBasespace;

		SetupKminmerSequenceExtractionFunctor(ToBasespace2& toBasespace) : _toBasespace(toBasespace){
		}

		SetupKminmerSequenceExtractionFunctor(const SetupKminmerSequenceExtractionFunctor& copy) : _toBasespace(copy._toBasespace){
		}

		~SetupKminmerSequenceExtractionFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {

			
			#pragma omp critical
			{

				for(u_int32_t i=0; i<kminmerList._kminmersInfo.size(); i++){
					
					//_logFile << readIndex << " " << i << endl;
					const ReadKminmerComplete& kminmerInfo = kminmerList._kminmersInfo[i];

					KmerVec vec = kminmerInfo._vec;

					bool orientation = !kminmerInfo._isReversed;

					ContigPosition contigPosition =  {(u_int32_t) kminmerList._readIndex, i};
					if(_toBasespace._contigPosition_to_readPosition.find(contigPosition) == _toBasespace._contigPosition_to_readPosition.end()){
						//Logger::get().debug() << "No read position for contig position";
						continue;
					}

					const ReadPositionMatch& readPositionMatch = _toBasespace._contigPosition_to_readPosition[contigPosition];
					ReadPosition readPosition = {readPositionMatch._readIndex, readPositionMatch._readPosition};
					//_toBasespace._kmerVecIndex[vec].push_back({(u_int32_t) kminmerList._readIndex, (u_int32_t) i});

					if(i == 0){
						//ReadSequence rs = {readIndex, nullptr, {}};
						_toBasespace._readPosition_entire[readPosition] = nullptr;
					}
					else {
						if(orientation){ //+
							//_logFile << "right: " << nodeName << " " << readIndex << endl;
							//ReadSequence rs = {readIndex, nullptr, {}};
							_toBasespace._readPosition_right[readPosition] = nullptr;
						}
						else{ //-
							//_logFile << "left: " << nodeName << " " << readIndex << endl;
							//ReadSequence rs = {readIndex, nullptr, {}};
							_toBasespace._readPosition_left[readPosition] = nullptr;
						}
					}

				}
			}
		}
	};

	/*
	void indexKminmers(){

		cerr << "Indexing reads" << endl;

		_kminmerNodeName = 0;

		KminmerParserParallel parser(_inputDir + "/read_data_init.txt", _minimizerSize, _kminmerSize, false, true, _nbCores);
		parser.parse(IndexKminmerFunctor(*this));

	}


	class IndexKminmerFunctor {

		public:

		ToBasespace& _parent;

		IndexKminmerFunctor(ToBasespace& parent) : _parent(parent){
		}

		IndexKminmerFunctor(const IndexKminmerFunctor& copy) : _parent(copy._parent){
		}

		~IndexKminmerFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {

			#pragma omp critical(indexKminmerFunctor)
			{

				for(u_int16_t i=0; i<kminmerList._kminmersInfo.size(); i++){
				
					
					const ReadKminmerComplete& kminmerInfo = kminmerList._kminmersInfo[i];

					KmerVec vec = kminmerInfo._vec;
					if(_parent._kmerVec_to_nodeName.find(vec) == _parent._kmerVec_to_nodeName.end()){
						_parent._kminmerNodeName += 1;
						_parent._kmerVec_to_nodeName[vec] = _parent._kminmerNodeName;
					}

				}
			}


		}
	};


	void indexReads(){


		_logFile << "Indexing reads" << endl;
		_unitigDatas.resize(_kmerVec_to_nodeName.size());
		
		KminmerParserParallel parser(_inputDir + "/read_data_init.txt", _minimizerSize, _kminmerSize, false, true, _nbCores);
		//parser2.parse(FilterKminmerFunctor2(*this));
		parser.parse(IndexReadsFunctor(*this));

		//KminmerParser parser(_inputDir + "/read_data_init.txt", _minimizerSize, _kminmerSize, false, true);
		//auto fp = std::bind(&ToBasespace::indexReads_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
		//parser.parse(fp);
	}

	class IndexReadsFunctor {

		public:

		ToBasespace& _parent;

		IndexReadsFunctor(ToBasespace& parent) : _parent(parent){
		}

		IndexReadsFunctor(const IndexReadsFunctor& copy) : _parent(copy._parent){
		}

		~IndexReadsFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {

			u_int64_t readIndex = kminmerList._readIndex;

			for(size_t i=0; i<kminmerList._kminmersInfo.size(); i++){
			
				
				//_logFile << readIndex << " " << i << endl;
				const ReadKminmerComplete& kminmerInfo = kminmerList._kminmersInfo[i];

				KmerVec vec = kminmerInfo._vec;
				
				if(_parent._kmerVec_to_nodeName.find(vec) == _parent._kmerVec_to_nodeName.end()){
					continue;
				}
				
				u_int32_t nodeName = _parent._kmerVec_to_nodeName[vec];

				#pragma omp critical(indexReadsFunctor)
				{
					if(_parent._unitigDatas[nodeName].size() < 50){
						_parent._unitigDatas[nodeName].push_back(readIndex);
					}
				}

				

			}

		}
	};

	*/

	void endCorrection(auto& copies, auto& isDoneNodeName, auto& correctedSequences, auto& repeatedKminmers_model, auto& repeatedKminmers_variants){

		/*
		_logFile << "-----" << endl;
		for(const auto& it : repeatedKminmers){
			_logFile << (it.second != nullptr) << endl;
			if(it.second == nullptr){
				_logFile << it.first._nodeIndex << " " << it.first._supportingReadIndex << endl;
				getchar();
			}

			if(it.first._nodeIndex == 602){
				_logFile << (it.second != nullptr) << endl;
				getchar();
			}
		}
		*/

		vector<u_int32_t> nodeNames;
		for(auto& it : copies){
			if(isDoneNodeName.find(it.first) != isDoneNodeName.end()) continue;
			nodeNames.push_back(it.first);
		}
		
		/*
		_logFile << nodeNames.size() << endl;
		#pragma omp parallel num_threads(_nbCores)
		{

			ExtractKminmerSequenceFunctor functor(_minimizerSize, _minimizerDensity, *this, false);

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
					//	_logFile << "omg" << endl;
					//	getchar();
					//}
				}
			}
		}

		_logFile << correctedSequences.size() << endl;
		*/



		
		vector<u_int32_t> readNodeNames;
		for(const auto& it : repeatedKminmers_variants){
			readNodeNames.push_back(it.first);
		}
		
		//"reprise: verifier que cette partie est bien paralleliser, ça a tourner a 100% de cpu la derniere run"
		Logger::get().debug() << "Correcting repeated kminmers: " << readNodeNames.size();
		#pragma omp parallel num_threads(_nbCores)
		{

			ExtractKminmerSequenceFunctor functor(_minimizerSize, _minimizerDensity, *this, false);

			#pragma omp for
			for(size_t i=0; i<readNodeNames.size(); i++){
				
				u_int32_t nodeName = readNodeNames[i];
				
				for(ReadSequence& readSequence : repeatedKminmers_variants[nodeName]){

					if(readSequence._sequence == nullptr){
						//_logFile << "No model sequence" << endl;
						continue;
					}

					VariantQueue& queue = readSequence._variants;

					//#pragma omp critical
					//queue = &repeatedKminmers_variants[readNodeName];

					//const vector<DnaBitset*>& kminmerCopies = copies[readNodeName._nodeIndex];


					vector<DnaBitset*> sequences;

					while(!queue.empty()){
						sequences.push_back(queue.top()._sequence);
						queue.pop();
					}

					sequences.push_back(readSequence._sequence);

					string correctedSequence;
					functor.performErrorCorrection_all(nodeName, sequences, correctedSequence);
					
					#pragma omp critical
					{
						//delete repeatedKminmers_model[nodeName];
						repeatedKminmers_model[{nodeName, readSequence._readIndex}] = new DnaBitset(correctedSequence);
					}

					for(DnaBitset* dnaSeq : sequences){
						delete dnaSeq;
					}
					sequences.clear();

				}
				//_logFile << readNodeName._nodeIndex << endl;




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

		_logFile << "Loading mContigs: " << contigFilename << endl;
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
			//_logFile << nodePath.size() << endl;
		}

		gzclose(contigFile);

		_logFile << "Nb contigs: " << nbContigs << endl;
	}
	*/
	//MDBG* _mdbgInit;
	//size_t _originalKminmerSize;
	u_int64_t _nbContigs;
	phmap::parallel_flat_hash_map<u_int32_t, u_int32_t> _kminmerCounts;

	/*
	void loadContigs_fasta (const string& contigFilename){


		_nbContigs = 0;

		auto fp = std::bind(&ToBasespace::loadContigs_fasta_read, this, std::placeholders::_1, std::placeholders::_2);
		ReadParser readParser(contigFilename, true, false);
		readParser.parse(fp);

		_logFile << "Nb contigs: " << _nbContigs << endl;

		double nbRepeatedKminmers = 0;
		for(auto& it : _kminmerCounts){
			if(it.second > 1){
				nbRepeatedKminmers += 1;
			}
		}
		_logFile << "Nb kminmer: " << _kminmerCounts.size() << endl;
		_logFile << "Repeated kminmer: " << nbRepeatedKminmers << endl;
		_logFile << "Repeated kminmer rate: " << (nbRepeatedKminmers / _kminmerCounts.size()) << endl;

		_kminmerCounts.clear();
	}
	*/






	/*
	void indexReads_read(const vector<u_int64_t>& minimizers, const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int8_t isCircular, u_int64_t readIndex){//}, const vector<KmerVec>& kminmers_k3, const vector<ReadKminmer>& kminmersInfos_k3){

		for(const KmerVec& vec : kminmers){
			
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()) continue;

			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
			
			if(_unitigDatas[nodeName].size() > 50) continue;

			//_logFile << nodeName << " " << readIndex << endl;
			_unitigDatas[nodeName].push_back(readIndex);
		}
	}
	*/

	/*
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
	*/

	//unordered_set<u_int32_t> _invalidContigIndex;
	/*
	void loadContigs_min(const string& contigFilename){

		_logFile << "Extracting kminmers: " << contigFilename << endl;
		KminmerParser parser2(contigFilename, _minimizerSize, _kminmerSize, false, false);
		auto fp2 = std::bind(&ToBasespace::loadContigs_min_read2, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
		parser2.parseMinspace(fp2);

	}

	void loadContigs_min_read2(const vector<u_int64_t>& readMinimizers, const vector<ReadKminmerComplete>& kminmersInfos, u_int8_t isCircular, u_int64_t readIndex){

		//_logFile << kminmersInfos.size() << " " << (_invalidContigIndex.find(readIndex) != _invalidContigIndex.end()) << endl;
		//if(_invalidContigIndex.find(readIndex) != _invalidContigIndex.end()) return;
		double s = 0;
		for(size_t i=0; i<kminmersInfos.size(); i++){
			
			const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];

			KmerVec vec = kminmerInfo._vec;
			
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()){
				_logFile << "Not found kminmer" << endl;
				//getchar();
				continue;
			}


			//if(kminmersInfos.size() <= 1) continue;

			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
			//s += nodeName;
			//nodepath.push_back(nodeName);
			//if(nodeName == 38376){

			//	_logFile << "lala" << endl;
			//	getchar();
			//}
			//_logFile << kminmersInfos.size() << " " << nodeName << endl;

			
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
			u_int32_t nodeIndex = BiGraph::nodeName_to_nodeIndex(nodeName, orientation);
			s += nodeIndex;

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

		_checksum += s*kminmersInfos.size();


	}
	*/


	/*
	void collectBestSupportingReads(const string& contigFilename){

		_nextReadIndexWriter = 0;
		_contigFileSupported = ofstream(contigFilename + ".tmp");

		_logFile << "Collecting best supporting reads" << endl;
		KminmerParserParallel parser(contigFilename, _minimizerSize, _kminmerSize, false, false, _nbCores);
		parser.parse(CollectBestSupportingReadsFunctor(*this));
		//KminmerParser parser(contigFilename, _minimizerSize, _kminmerSize, false, false);
		//auto fp = std::bind(&ToBasespace::collectBestSupportingReads_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
		//parser.parse(fp);


		_contigFileSupported.close();

	}

	
    struct CollectBestSupportingReadsFunctor_Writer{
        u_int64_t _readIndex;
        vector<u_int64_t> _supportingReads;
        //u_int32_t _prevNodeIndex;
    };

    struct CollectBestSupportingReadsFunctor_WriterComparator {
        bool operator()(CollectBestSupportingReadsFunctor_Writer const& p1, CollectBestSupportingReadsFunctor_Writer const& p2){
            return p1._readIndex > p2._readIndex;
        }
    };

	priority_queue<CollectBestSupportingReadsFunctor_Writer, vector<CollectBestSupportingReadsFunctor_Writer> , CollectBestSupportingReadsFunctor_WriterComparator> _readWriterQueue;
	u_int64_t _nextReadIndexWriter;

	class CollectBestSupportingReadsFunctor {

		public:

		ToBasespace& _toBasespace;

		CollectBestSupportingReadsFunctor(ToBasespace& toBasespace) : _toBasespace(toBasespace){
		}

		CollectBestSupportingReadsFunctor(const CollectBestSupportingReadsFunctor& copy) : _toBasespace(copy._toBasespace){
		}

		~CollectBestSupportingReadsFunctor(){
		}

		//void collectBestSupportingReads_read(const vector<u_int64_t>& readMinimizers, const vector<ReadKminmerComplete>& kminmersInfos, u_int8_t isCircular, u_int64_t readIndex){

		void operator () (const KminmerList& kminmerList) {
			//if(_invalidContigIndex.find(readIndex) != _invalidContigIndex.end()) return;

			//_logFile << "----------------------" << endl;
			vector<u_int64_t> supportingReads;

			//cout << kminmerList._readMinimizers.size() << endl;
			//_logFile << readIndex << " " << kminmersInfos.size() << endl;
			for(size_t i=0; i<kminmerList._kminmersInfo.size(); i++){
				
				//_logFile << readIndex << " " << i << endl;
				const ReadKminmerComplete& kminmerInfo = kminmerList._kminmersInfo[i];

				KmerVec vec = kminmerInfo._vec;
				
				if(_toBasespace._kmerVec_to_nodeName.find(vec) == _toBasespace._kmerVec_to_nodeName.end()){
					_toBasespace._logFile << "Not found kminmer" << endl;
					supportingReads.push_back(-1);
					//getchar();
					continue;
				}

				bool orientation = !kminmerInfo._isReversed;

				
				u_int32_t nodeName = _toBasespace._kmerVec_to_nodeName[vec];

				//_kminmerCounts.clear();
				//if(_kminmerCounts[nodeName] > 1){
					//_logFile << nodeName << endl;
					//isKminmerRepeated.insert(nodeName);
					u_int64_t readIndex = getBestSupportingRead(nodeName, i, kminmerList._kminmersInfo);
					supportingReads.push_back(readIndex);
					
					ReadNodeName readNodeName = {nodeName, readIndex};

					#pragma omp critical(bestSupport)
					{
						if(i == 0){
							ReadSequence rs = {readIndex, nullptr, {}};
							_toBasespace._kminmerSequence_entire_multi[nodeName].push_back(rs);
						}
						else {
							if(orientation){ //+
								//_logFile << "right: " << nodeName << " " << readIndex << endl;
								ReadSequence rs = {readIndex, nullptr, {}};
								_toBasespace._kminmerSequence_right_multi[nodeName].push_back(rs);
							}
							else{ //-
								//_logFile << "left: " << nodeName << " " << readIndex << endl;
								ReadSequence rs = {readIndex, nullptr, {}};
								_toBasespace._kminmerSequence_left_multi[nodeName].push_back(rs);
							}
						}
					}

			}

			//_toBasespace._nbContigs += 1;

			#pragma omp critical(bestSupport)
			{
				_toBasespace._readWriterQueue.push({kminmerList._readIndex, supportingReads});
				//_logFile << _readWriterQueue.size() << " " << read._index << " " << _nextReadIndexWriter << endl;

				while(!_toBasespace._readWriterQueue.empty()){

					const CollectBestSupportingReadsFunctor_Writer& readWriter = _toBasespace._readWriterQueue.top();

					if(readWriter._readIndex == _toBasespace._nextReadIndexWriter){

						for(u_int64_t readIndex : readWriter._supportingReads){
							_toBasespace._bestSupportChecksum += readIndex;
						}

						u_int32_t size = readWriter._supportingReads.size();
						_toBasespace._contigFileSupported.write((const char*)&size, sizeof(size));
						_toBasespace._contigFileSupported.write((const char*)&readWriter._supportingReads[0], size*sizeof(u_int64_t));

						_toBasespace._readWriterQueue.pop();
						_toBasespace._nextReadIndexWriter += 1;
					}
					else{
						break;
					}
				}
				
			}



		}

		u_int64_t getBestSupportingRead(u_int32_t nodeName, size_t nodeNamePosition, const vector<ReadKminmerComplete>& kminmersInfos){

			const vector<u_int32_t>& reads = _toBasespace._unitigDatas[nodeName];

			u_int64_t maxNbmatches = 0;
			u_int64_t bestReadIndex = -1;

			for(u_int64_t readIndex : reads){
				u_int64_t nbMatchesRight = getMatchDirection(nodeName, nodeNamePosition+1, kminmersInfos, 1, readIndex);
				u_int64_t nbMatchesLeft = getMatchDirection(nodeName, nodeNamePosition-1, kminmersInfos, -1, readIndex);
				u_int64_t nbMatches = nbMatchesRight + nbMatchesLeft + 1; //+1 seed match not counted
				if(nbMatches > maxNbmatches){
					maxNbmatches = nbMatches;
					bestReadIndex = readIndex;
				}
			}

			return bestReadIndex;


		}

		u_int64_t getBestSupportingRead_direction(u_int32_t nodeName, size_t nodeNamePosition, const vector<ReadKminmerComplete>& kminmersInfos, int inc){
			
			u_int64_t readIndex = -1;
			vector<u_int32_t> sharedElements;

			long i = nodeNamePosition;

			while(true){
				if(i < 0 || i >= kminmersInfos.size()){
					break;
				}

				const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];

				KmerVec vec = kminmerInfo._vec;

				if(_toBasespace._kmerVec_to_nodeName.find(vec) != _toBasespace._kmerVec_to_nodeName.end()){

					u_int32_t nodeName2 = _toBasespace._kmerVec_to_nodeName[vec];
					
					if(sharedElements.size() == 0){
						Utils::collectSharedElements(_toBasespace._unitigDatas[nodeName], _toBasespace._unitigDatas[nodeName2], sharedElements);
					}
					else{
						vector<u_int32_t> sharedElementTmp;
						Utils::collectSharedElements(sharedElements, _toBasespace._unitigDatas[nodeName2], sharedElementTmp);
						sharedElements = sharedElementTmp;
					}

					if(sharedElements.size() == 0) break;


					readIndex = sharedElements[0];



				}

				i += inc;

			}

			return readIndex;


		}

		u_int64_t getMatchDirection(u_int32_t nodeName, size_t nodeNamePosition, const vector<ReadKminmerComplete>& kminmersInfos, int inc, u_int64_t readIndex){
			
			//vector<u_int64_t> readIndexSource = {readIndex};
			u_int64_t nbMatches = 0;

			long i = nodeNamePosition;

			while(true){
				if(i < 0 || i >= kminmersInfos.size()){
					break;
				}

				const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];

				KmerVec vec = kminmerInfo._vec;

				if(_toBasespace._kmerVec_to_nodeName.find(vec) != _toBasespace._kmerVec_to_nodeName.end()){

					u_int32_t nodeName2 = _toBasespace._kmerVec_to_nodeName[vec];

					vector<u_int32_t>& readIndexes = _toBasespace._unitigDatas[nodeName2];
					if(std::find(readIndexes.begin(), readIndexes.end(), readIndex) == readIndexes.end()){
						break; 
					}

				}
				else{
					break;
				}

				nbMatches += 1;
				i += inc;

			}

			return nbMatches;


		}

	};

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
		_logFile << "-------------------" << endl;
		_logFile << "-------------------" << endl;
		_logFile << "-------------------" << endl;
		_logFile << kminmerInfo._isReversed << endl;
		_logFile << sequence << endl;
		*/

		if(loadType == LoadType::Entire){
			return;
			//startPosition = kminmerInfo._read_pos_start;
			//len = kminmerInfo._read_pos_end - kminmerInfo._read_pos_start;
		}
		else if(loadType == LoadType::Left){
			startPosition = 0; //kminmerInfo._read_pos_start;
			len = kminmerInfo._seq_length_start; //kminmerInfo._position_of_second_minimizer_seq - kminmerInfo._read_pos_start;
			//_logFile << kminmerInfo._read_pos_start << " " << kminmerInfo._position_of_second_minimizer << endl;
		}
		else if(loadType == LoadType::Right){
			//return;
			startPosition = sequence.size() - kminmerInfo._seq_length_end;  //kminmerInfo._position_of_second_to_last_minimizer_seq;
			len = kminmerInfo._seq_length_end; //kminmerInfo._read_pos_end - kminmerInfo._position_of_second_to_last_minimizer_seq;
			//return;
			//_logFile << kminmerInfo._read_pos_end << " " << kminmerInfo._position_of_second_to_last_minimizer << endl;
		}
		
		//char* seq = sequence[0];
		//char subbuff2[len+1];
		//memcpy( subbuff, &seq[startPosition], len);
		//subbuff2[len] = '\0';
		//sequence = string(subbuff2);
		/*
		_logFile << sequence.size() << " " << startPosition << " " << len  << endl;
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
		_logFile << startPosition << " " << len << endl;
		_logFile << "allo" << endl;
		_logFile << endl << endl;
		//_logFile << 
		_logFile << sequence << endl;
		//Utils::revcomp(sequence);
		//_logFile << sequence << endl;
		*/
	}

	void extractKminmerSequences_model (){

		Logger::get().debug() << "Extracting kminmer models" ;
		ReadParserParallel readParser(_inputFilename, false, false, _nbCores);
		readParser.parse(ExtractKminmerSequenceFunctor(_minimizerSize, _minimizerDensity, *this, true));
	}

	void extractKminmerSequences_allVariants (){

		Logger::get().debug() << "Extracting kminmer sequence variants";

		//auto fp = std::bind(&ToBasespace::extractKminmerSequences_allVariants_read, this, std::placeholders::_1);
		ReadParserParallel readParser(_inputFilename, false, false, _nbCores);
		readParser.parse(ExtractKminmerSequenceFunctor(_minimizerSize, _minimizerDensity, *this, false));

	}

	class ExtractKminmerSequenceFunctor{

		public:

		bool _collectModel;
		ToBasespace2& _toBasespace;
		EncoderRLE _encoderRLE;
		MinimizerParser* _minimizerParser;
		size_t _minimizerSize;
		float _minimizerDensity;

		ExtractKminmerSequenceFunctor(size_t minimizerSize, float minimizerDensity, ToBasespace2& toBasespace, bool collectModel) : _toBasespace(toBasespace){
			_minimizerSize = minimizerSize;
			_minimizerDensity = minimizerDensity;
			_minimizerParser = new MinimizerParser(minimizerSize, minimizerDensity);
			
			_collectModel = collectModel;
		}

		ExtractKminmerSequenceFunctor(const ExtractKminmerSequenceFunctor& copy) : _toBasespace(copy._toBasespace){
			_minimizerSize = copy._minimizerSize;
			_minimizerDensity = copy._minimizerDensity;
			_minimizerParser = new MinimizerParser(_minimizerSize, _minimizerDensity);

			_collectModel = copy._collectModel;


		}

		~ExtractKminmerSequenceFunctor(){
			delete _minimizerParser;
		}

		void operator () (const Read& read) {
			u_int32_t readIndex = read._index;

			if(readIndex % 100000 == 0){
				Logger::get().debug() << "Correcting kminmer " << readIndex;
			}
			//ottalSize += strlen(read->seq.s);
						
			string kminmerSequence;
			const char* sequenceOriginal = read._seq.c_str();

			string rleSequence;
			vector<u_int64_t> rlePositions;
			_encoderRLE.execute(read._seq.c_str(), read._seq.size(), rleSequence, rlePositions, false);

			vector<MinimizerType> minimizers;
			vector<u_int32_t> minimizers_pos;
			vector<u_int8_t> minimizers_direction;
			_minimizerParser->parse(rleSequence, minimizers, minimizers_pos, minimizers_direction);

			vector<KmerVec> kminmers; 
			vector<ReadKminmer> kminmersInfo;
			MDBG::getKminmers(_minimizerSize, _toBasespace._kminmerSize, minimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, readIndex, false);

			for(u_int32_t i=0; i<kminmers.size(); i++){
				
				
				ReadPosition readPosition = {readIndex, i};
				ReadKminmer& kminmerInfo = kminmersInfo[i];

				bool isEntire = false;
				bool isLeft = false;
				bool isRight = false;

				if(_toBasespace._readPosition_entire.find(readPosition) != _toBasespace._readPosition_entire.end()){
					isEntire = true;
				}
				if(_toBasespace._readPosition_left.find(readPosition) != _toBasespace._readPosition_left.end()){
					isLeft = true;
				}
				if(_toBasespace._readPosition_right.find(readPosition) != _toBasespace._readPosition_right.end()){
					isRight = true;
				}


				if(isEntire){
					_toBasespace.extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Entire, kminmerSequence);
					_toBasespace._readPosition_entire[readPosition] = new DnaBitset(kminmerSequence);
					//addKminmerSequenceVariant_add_all(nodeName, readIndex, _toBasespace._kminmerSequenceCopies_all_entire, kminmerSequence, _toBasespace._isDoneNodeName_entire, _toBasespace._kminmerSequence_entire, _toBasespace._kminmerSequence_entire_multi);
				}
				
				if(isLeft){
					_toBasespace.extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Left, kminmerSequence);
					_toBasespace._readPosition_left[readPosition] = new DnaBitset(kminmerSequence);
					//addKminmerSequenceVariant_add_all(nodeName, readIndex, _toBasespace._kminmerSequenceCopies_all_left, kminmerSequence, _toBasespace._isDoneNodeName_left, _toBasespace._kminmerSequence_left, _toBasespace._kminmerSequence_left_multi);
				}
				if(isRight){
					_toBasespace.extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Right, kminmerSequence);
					_toBasespace._readPosition_right[readPosition] = new DnaBitset(kminmerSequence);
					//addKminmerSequenceVariant_add_all(nodeName, readIndex, _toBasespace._kminmerSequenceCopies_all_right, kminmerSequence, _toBasespace._isDoneNodeName_right, _toBasespace._kminmerSequence_right, _toBasespace._kminmerSequence_right_multi);
				}


				/*
				if(_toBasespace._kmerVec_to_nodeName.find(kminmers[i]) == _toBasespace._kmerVec_to_nodeName.end()) continue;

				u_int32_t nodeName = _toBasespace._kmerVec_to_nodeName[kminmers[i]];
				ReadKminmer& kminmerInfo = kminmersInfo[i];

				//if(nodeName == 3114) _logFile << ">>>>>>>>>>> " << nodeName << endl;
				if(_collectModel){
					//if(_toBasespace.isKminmerRepeated.find(nodeName) == _toBasespace.isKminmerRepeated.end()) continue;
				}
				*/
				
				//3908 7040
				//298 8320
				//if(readIndex == 8320){
				//	_logFile << nodeName << endl;
				//}

				//_logFile << _kminmerSequenceCopies_all_entire.size() << " " << _kminmerSequenceCopies_all_left.size() << " " << _kminmerSequenceCopies_all_right.size() << " " << _kminmerSequence_entire.size() << " " << _kminmerSequence_left.size() << " " << _kminmerSequence_right.size() << endl;

				/*
				bool isEntire = false;
				bool isLeft = false;
				bool isRight = false;

				#pragma omp critical(indexKminmer)
				{
					if(_toBasespace._kminmerSequence_entire_multi.find(nodeName) != _toBasespace._kminmerSequence_entire_multi.end()){
						for(const ReadSequence& rs : _toBasespace._kminmerSequence_entire_multi[nodeName]){
							if(rs._readIndex == readIndex){
								isEntire = true;
							}
						}
					}

					if(_toBasespace._kminmerSequence_left_multi.find(nodeName) != _toBasespace._kminmerSequence_left_multi.end()){
						for(const ReadSequence& rs : _toBasespace._kminmerSequence_left_multi[nodeName]){
							if(rs._readIndex == readIndex){
								isLeft = true;
							}
						}
					}

					if(_toBasespace._kminmerSequence_right_multi.find(nodeName) != _toBasespace._kminmerSequence_right_multi.end()){
						for(const ReadSequence& rs : _toBasespace._kminmerSequence_right_multi[nodeName]){
							if(rs._readIndex == readIndex){
								isRight = true;
							}
						}
					}
				}
				*/

				//#pragma omp critical
				//{
					/*
					if(_toBasespace.isKminmerRepeated.find(nodeName) == _toBasespace.isKminmerRepeated.end()){
						isEntire = _toBasespace._kminmerSequenceCopies_all_entire.find(nodeName) != _toBasespace._kminmerSequenceCopies_all_entire.end(); // && _toBasespace._isDoneNodeName_entire.find(nodeName) == _toBasespace._isDoneNodeName_entire.end();
						isLeft = _toBasespace._kminmerSequenceCopies_all_left.find(nodeName) != _toBasespace._kminmerSequenceCopies_all_left.end(); // && _toBasespace._isDoneNodeName_left.find(nodeName) == _toBasespace._isDoneNodeName_left.end();
						isRight = _toBasespace._kminmerSequenceCopies_all_right.find(nodeName) != _toBasespace._kminmerSequenceCopies_all_right.end(); // && _toBasespace._isDoneNodeName_right.find(nodeName) == _toBasespace._isDoneNodeName_right.end();
					}
					else{

						if(_collectModel){
							if(_toBasespace._kminmerSequence_entire_multi.find(nodeName) != _toBasespace._kminmerSequence_entire_multi.end()){
								for(const ReadSequence& rs : _toBasespace._kminmerSequence_entire_multi[nodeName]){
									if(rs._readIndex == readIndex){
										isEntire = true;
									}
								}
							}

							if(_toBasespace._kminmerSequence_left_multi.find(nodeName) != _toBasespace._kminmerSequence_left_multi.end()){
								for(const ReadSequence& rs : _toBasespace._kminmerSequence_left_multi[nodeName]){
									if(rs._readIndex == readIndex){
										isLeft = true;
									}
								}
							}

							if(_toBasespace._kminmerSequence_right_multi.find(nodeName) != _toBasespace._kminmerSequence_right_multi.end()){
								for(const ReadSequence& rs : _toBasespace._kminmerSequence_right_multi[nodeName]){
									if(rs._readIndex == readIndex){
										isRight = true;
									}
								}
							}
						}
						else{
							isEntire = _toBasespace._kminmerSequence_entire_multi.find(nodeName) != _toBasespace._kminmerSequence_entire_multi.end();
							isLeft = _toBasespace._kminmerSequence_left_multi.find(nodeName) != _toBasespace._kminmerSequence_left_multi.end();
							isRight = _toBasespace._kminmerSequence_right_multi.find(nodeName) != _toBasespace._kminmerSequence_right_multi.end();
						}


						//if(nodeName == 3114) _logFile << isEntire << " " << isLeft << " " << isRight << endl;
						//ReadNodeName readNodeName = {nodeName, readIndex};
						//isEntire = _toBasespace._repeatedKminmerSequence_entire.find(readNodeName) != _toBasespace._repeatedKminmerSequence_entire.end();
						//isLeft = _toBasespace._repeatedKminmerSequence_left.find(readNodeName) != _toBasespace._repeatedKminmerSequence_left.end();
						//isRight = _toBasespace._repeatedKminmerSequence_right.find(readNodeName) != _toBasespace._repeatedKminmerSequence_right.end();
					}
					*/


				//}

				//if(nodeName==293 && readIndex==8320){
				//	ReadNodeName readNodeName = {nodeName, readIndex};
				//	_logFile << (_toBasespace._repeatedKminmerSequence_entire.find(readNodeName) != _toBasespace._repeatedKminmerSequence_entire.end()) << endl;
				//	_logFile << "lol" << endl;
				//	_logFile << isEntire << endl;
				//	getchar();
				//}
				/*
				if(isEntire){
					_toBasespace.extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Entire, kminmerSequence);
					addKminmerSequenceVariant_add_all(nodeName, readIndex, _toBasespace._kminmerSequenceCopies_all_entire, kminmerSequence, _toBasespace._isDoneNodeName_entire, _toBasespace._kminmerSequence_entire, _toBasespace._kminmerSequence_entire_multi);
				}
				
				if(isLeft){
					_toBasespace.extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Left, kminmerSequence);
					addKminmerSequenceVariant_add_all(nodeName, readIndex, _toBasespace._kminmerSequenceCopies_all_left, kminmerSequence, _toBasespace._isDoneNodeName_left, _toBasespace._kminmerSequence_left, _toBasespace._kminmerSequence_left_multi);
				}
				if(isRight){
					_toBasespace.extractKminmerSequence(sequenceOriginal, kminmerInfo, LoadType::Right, kminmerSequence);
					addKminmerSequenceVariant_add_all(nodeName, readIndex, _toBasespace._kminmerSequenceCopies_all_right, kminmerSequence, _toBasespace._isDoneNodeName_right, _toBasespace._kminmerSequence_right, _toBasespace._kminmerSequence_right_multi);
				}
				*/

				//if(nodeName==293 && readIndex==8320){
				//	ReadNodeName readNodeName = {nodeName, readIndex};
				//	//_logFile << (_toBasespace._repeatedKminmerSequence_left.find(readNodeName) != _toBasespace._repeatedKminmerSequence_left.end()) << " " << (_toBasespace._repeatedKminmerSequence_right.find(readNodeName) != _toBasespace._repeatedKminmerSequence_right.end()) << endl;
				//	_logFile << "lol" << endl;
				//	getchar();
				//}

			}

			//if(readIndex == 7210){
			//	_logFile << "lal" << endl;
			//	getchar();
			//}

		}


		void addKminmerSequenceVariant_add_all(u_int32_t nodeName, u_int64_t readIndex, auto& variants, const string& sequence, auto& isDoneNodeName, auto& correctedSequences, auto& repeatedKminmers){
			
			//static EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0);

			//_logFile << (isDoneNodeName.find(nodeName) != isDoneNodeName.end()) << endl;
			//if(sequenceModel == nullptr){
			//	_logFile << "pas normal" << endl;
			//	return; //model not found yet
			//}

			//_logFile << nodeName << " " << sequence.size() << endl;
						

			vector<DnaBitset*>* queue;
			VariantQueue* variantQueue;
			DnaBitset* dnaSeq_model;

			vector<ReadSequence>* readSequences = nullptr;

			bool isRepeated = false;
			bool cancel = false;
			
			#pragma omp critical(indexKminmer)
			{
				
				if(_collectModel){

					for(ReadSequence& rs : repeatedKminmers[nodeName]){
						if(rs._readIndex == readIndex){
							rs._sequence = new DnaBitset(sequence);
						}
					}

					//ReadNodeName readNodeName = {nodeName, readIndex};
					//repeatedKminmers_model[readNodeName] = new DnaBitset(sequence);
					cancel = true;

					//_logFile << repeatedKminmers_model.size() << " " << _toBasespace.isKminmerRepeated.size() << endl;
					//_logFile << nodeName << " " << readIndex << endl;
					//if(nodeName == 24 && readIndex == 6){
					//	getchar();
					//}
				}
				else{
					/*
					if(_toBasespace.isKminmerRepeated.find(nodeName) == _toBasespace.isKminmerRepeated.end()){
						if(isDoneNodeName.find(nodeName) != isDoneNodeName.end()) cancel = true;

						if(!cancel){
							queue = &variants[nodeName];
							queue->push_back(new DnaBitset(sequence));

							//if(_toBasespace.isKminmerRepeated.find(nodeName) != _toBasespace.isKminmerRepeated.end()) cancel = true;

							if(queue->size() < 20) cancel = true;
							
							if(!cancel) isDoneNodeName.insert(nodeName);
						}

						isRepeated = false;
					}
					else{*/

						//_logFile << "1" << endl;
						cancel = true;
						isRepeated = true;

						//if(repeatedKminmers.find(nodeName) != repeatedKminmers.end()){
						readSequences = &repeatedKminmers[nodeName];
						//}
						/*
						ReadNodeName readNodeName = {nodeName, readIndex};
						dnaSeq_model = repeatedKminmers_model[readNodeName];
						//_logFile << nodeName << " " << readIndex << " " << (dnaSeq_model != nullptr) << " " << (_toBasespace.isKminmerRepeated.find(nodeName) != _toBasespace.isKminmerRepeated.end()) << endl;
						variantQueue = &repeatedKminmers_variants[readNodeName];

						for(const ReadSequence& rs : _toBasespace._kminmerSequence_right_multi){
							if(rs._readIndex == readIndex){
								isRight = true;
							}
						}
						*/

					//}
				}




			}


			if(cancel) return;

			string correctedSequence;
			performErrorCorrection_all(nodeName, *queue, correctedSequence);

			#pragma omp critical
			{
				correctedSequences[nodeName] = new DnaBitset(correctedSequence);
			}


		}


		void performErrorCorrection_all(u_int32_t nodeName, const vector<DnaBitset*>& sequences, string& correctedSequence){
			
		}

	};

	phmap::parallel_flat_hash_set<u_int32_t> _isDoneNodeName_entire;
	phmap::parallel_flat_hash_set<u_int32_t> _isDoneNodeName_left;
	phmap::parallel_flat_hash_set<u_int32_t> _isDoneNodeName_right;

	//void extractKminmerSequences_allVariants_read(const Read& read){


	//}



	u_int64_t _nbBps;
	u_int64_t _checksum;


	void writeKminmerSequence_all(u_int32_t nodeName, const string& sequence, const gzFile& file){
		
		//char* dnaStr = dna->to_string();

		//const string& sequence = 
		u_int16_t length = sequence.size();

		//_logFile << "Writing: " << endl;
		//_logFile << nodeName << " " << length << " " << sequence << endl;
		gzwrite(file, (const char*)&nodeName, sizeof(nodeName));
		gzwrite(file, (const char*)&length, sizeof(length));
		gzwrite(file, (const char*)&sequence[0], length);
	}


	gzFile _basespaceContigFile;
	u_int64_t _contigIndex;

	/*
	void createBaseContigs (const string& contigFilename, const string& outputFilename){

		_logFile << "Creating basespace contigs: " << contigFilename << " " << outputFilename << endl;

		_basespaceContigFile = gzopen(outputFilename.c_str(),"wb");

		auto fp = std::bind(&ToBasespace::createBaseContigs_read, this, std::placeholders::_1, std::placeholders::_2);
		ReadParser readParser(contigFilename, true, false);
		readParser.parse(fp);

		
		gzclose(_basespaceContigFile);
	}
	*/
	
	void createBaseContigs(const string& contigFilename){

		_nbBps = 0;

		_contigFileSupported_input = ifstream(contigFilename + ".tmp");

		Logger::get().debug() << "Creating basespace contigs: " << contigFilename;


		//KminmerParserParallel parser(contigFilename, _minimizerSize, _kminmerSize, false, false);
		KminmerParserParallel parser(_inputFilenameContig, _minimizerSize, _kminmerSize, false, false, 1);
		parser.parse(CreateBaseContigsFunctor(*this));
		//auto fp = std::bind(&ToBasespace::createBaseContigs_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
		//parser.parseMinspace(fp);

		_contigFileSupported_input.close();


		Logger::get().debug() << "Nb contigs: " << (_contigIndex);
		Logger::get().debug() << "Nb bps: " << _nbBps;
		Logger::get().debug() << "Checksum: " << _checksum;
		//Logger::get().debug() << "Checksum (best support): " << _bestSupportChecksum;

	}



	class CreateBaseContigsFunctor {

		public:

		ToBasespace2& _toBasespace;

		CreateBaseContigsFunctor(ToBasespace2& toBasespace) : _toBasespace(toBasespace){
		}

		CreateBaseContigsFunctor(const CreateBaseContigsFunctor& copy) : _toBasespace(copy._toBasespace){
		}

		~CreateBaseContigsFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {

			ReadType readIndex = kminmerList._readIndex;
			//MinimizerReadPacked read = {kminmerList._readIndex, kminmerList._readMinimizers, kminmerList._minimizerPos, kminmerList._meanReadQuality, kminmerList._readLength};
			//MinimizerReadPacked readLowDensity = Utils::getLowDensityMinimizerRead(read, _parent._minimizerDensity_assembly);
			

			string contigSequence = "";
			bool isCircular = kminmerList._isCircular;

			for(u_int32_t i=0; i<kminmerList._kminmersInfo.size(); i++){
			
				
				//cout << "------- " << i << endl;

				//_logFile << readIndex << " " << i << endl;
				const ReadKminmerComplete& kminmerInfo = kminmerList._kminmersInfo[i];

				const KmerVec& vec = kminmerInfo._vec;
				
				
				
				ContigPosition contigPosition =  {(u_int32_t) readIndex, i};
				const ReadPositionMatch& readPositionMatch = _toBasespace._contigPosition_to_readPosition[contigPosition];
				const ReadPosition& readPosition = {readPositionMatch._readIndex, readPositionMatch._readPosition};

				//_logFile << i << " " << supportingReads.size() << endl;
				//u_int64_t readIndex = supportingReads[i];
				//_logFile << readIndex << endl;
				//if(readIndex == -1){
				//	_logFile << "No sequence for kminmer: " << readIndex << endl;
				//	continue;
				//}
				
				//_logFile << i << endl;
				//const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];

				//KmerVec vec = kminmerInfo._vec;
				
				//if(_kmerVec_to_nodeName.find(vec) == _kmerVec_to_nodeName.end()){
				//	_logFile << "Not found node name" << endl;
				//	continue;
				//}

				//u_int32_t nodeName = _kmerVec_to_nodeName[vec];
				//if(nodeName == 38376){
				//	_logFile << "lala" << endl;
				//	getchar();
				//}
				//_logFile << nodeName << endl;
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

						char* seq = nullptr;

						if(_toBasespace._readPosition_entire.find(readPosition) != _toBasespace._readPosition_entire.end() && _toBasespace._readPosition_entire[readPosition] != nullptr){
							seq = _toBasespace._readPosition_entire[readPosition]->to_string();
						}

						
						if(seq == nullptr){
							//cerr << "No sequence for kminmer" << endl;
							continue;
						}

						/*
						char* seq = nullptr;

						if(_kminmerSequence_entire_multi.find(nodeName) != _kminmerSequence_entire_multi.end()){
							for(const ReadSequence& rs : _kminmerSequence_entire_multi[nodeName]){
								if(rs._readIndex == readIndex){
									//_logFile << (rs._sequence != nullptr) << endl;
									seq = rs._sequence->to_string();
									break;
								}
							}
						}

						if(seq == nullptr){
							_logFile << "No sequence for kminmer" << endl;
							continue;
						}
						*/

						string kminmerSequence = string(seq);
						free(seq);

						contigSequence += kminmerSequence;
						//_logFile << contigSequence << endl;
					}
					else{
						
						char* seq = nullptr;

						if(_toBasespace._readPosition_entire.find(readPosition) != _toBasespace._readPosition_entire.end() && _toBasespace._readPosition_entire[readPosition] != nullptr){
							seq = _toBasespace._readPosition_entire[readPosition]->to_string();
						}
						
						if(seq == nullptr){
							//cerr << "No sequence for kminmer" << endl;
							continue;
						}

						/*
						char* seq = nullptr;

						if(_kminmerSequence_entire_multi.find(nodeName) != _kminmerSequence_entire_multi.end()){
							for(const ReadSequence& rs : _kminmerSequence_entire_multi[nodeName]){
								if(rs._readIndex == readIndex){

									//_logFile << (rs._sequence != nullptr) << endl;

									seq = rs._sequence->to_string();
									break;
								}
							}
						}

						if(seq == nullptr){
							_logFile << "No sequence for kminmer" << endl;
							continue;
						}
						*/

						/*
						if(_kminmerSequence_entire.find(nodeName) != _kminmerSequence_entire.end()){
							if(_kminmerSequence_entire[nodeName] == nullptr){
								_logFile << "No sequence for kminmer" << endl;
								continue;
							}
							seq = _kminmerSequence_entire[nodeName]->to_string();
						}
						else{
							ReadNodeName readNodeName = {nodeName, readIndex};
							if(_kminmerSequence_entire_multi.find(readNodeName) == _kminmerSequence_entire_multi.end() || _kminmerSequence_entire_multi[readNodeName] == nullptr){
								_logFile << "No sequence for kminmer" << endl;
								continue;
							}
							seq = _kminmerSequence_entire_multi[readNodeName]->to_string();
						//}
						*/

						string kminmerSequence = string(seq);
						free(seq);

						Utils::revcomp(kminmerSequence);
						contigSequence += kminmerSequence;
					}
				}
				else {
					if(orientation){

						char* seq = nullptr;

						if(_toBasespace._readPosition_right.find(readPosition) != _toBasespace._readPosition_right.end() && _toBasespace._readPosition_right[readPosition] != nullptr){
							seq = _toBasespace._readPosition_right[readPosition]->to_string();
						}

						if(seq == nullptr){
							//cerr << "No sequence for kminmer" << endl;
							continue;
						}

						/*
						char* seq = nullptr;

						if(_kminmerSequence_right_multi.find(nodeName) != _kminmerSequence_right_multi.end()){
							for(const ReadSequence& rs : _kminmerSequence_right_multi[nodeName]){
								if(rs._readIndex == readIndex){
									seq = rs._sequence->to_string();
									break;
								}
							}
						}

						if(seq == nullptr){
							_logFile << "No sequence for kminmer" << endl;
							continue;
						}
						*/
						/*
						if(_kminmerSequence_right.find(nodeName) != _kminmerSequence_right.end()){
							if(_kminmerSequence_right[nodeName] == nullptr){
								_logFile << "No sequence for kminmer" << endl;
								continue;
							}
							seq = _kminmerSequence_right[nodeName]->to_string();
						}
						else{
							ReadNodeName readNodeName = {nodeName, readIndex};
							if(_kminmerSequence_entire_right.find(readNodeName) == _kminmerSequence_entire_right.end() || _kminmerSequence_entire_right[readNodeName] == nullptr){
								_logFile << "No sequence for kminmer" << endl;
								continue;
							}
							seq = _kminmerSequence_entire_right[readNodeName]->to_string();
						//}
						*/

						string kminmerSequence = string(seq);
						free(seq);
						
						contigSequence += kminmerSequence;
						

					}
					else{
						
						char* seq = nullptr;

						if(_toBasespace._readPosition_left.find(readPosition) != _toBasespace._readPosition_left.end() && _toBasespace._readPosition_left[readPosition] != nullptr){
							seq = _toBasespace._readPosition_left[readPosition]->to_string();
						}
						
						if(seq == nullptr){
							//cerr << "No sequence for kminmer" << endl;
							continue;
						}

						/*
						char* seq = nullptr;

						if(_kminmerSequence_left_multi.find(nodeName) != _kminmerSequence_left_multi.end()){
							for(const ReadSequence& rs : _kminmerSequence_left_multi[nodeName]){
								if(rs._readIndex == readIndex){
									seq = rs._sequence->to_string();
									break;
								}
							}
						}

						if(seq == nullptr){
							_logFile << "No sequence for kminmer" << endl;
							continue;
						}
						*/
						/*
						if(_kminmerSequence_left.find(nodeName) != _kminmerSequence_left.end()){
							if(_kminmerSequence_left[nodeName] == nullptr){
								_logFile << "No sequence for kminmer" << endl;
								continue;
							}
							seq = _kminmerSequence_left[nodeName]->to_string();
						}
						else{
							ReadNodeName readNodeName = {nodeName, readIndex};
							if(_kminmerSequence_entire_left.find(readNodeName) == _kminmerSequence_entire_left.end() || _kminmerSequence_entire_left[readNodeName] == nullptr){
								_logFile << "No sequence for kminmer" << endl;
								continue;
							}
							seq = _kminmerSequence_entire_left[readNodeName]->to_string();
						//}
						*/

						string kminmerSequence = string(seq);
						free(seq);

						Utils::revcomp(kminmerSequence);
						contigSequence += kminmerSequence;


					}
				}
			}

			if(contigSequence.size() == 0){
				//Logger::get().debug() << "Empty contig " << kminmerList._kminmersInfo.size();
				//return;
			}
			else{
				//cout << "lala:" << ((u_int32_t) isCircular) << endl;
				string linearOrCircular;
				if(isCircular == CONTIG_LINEAR){
					linearOrCircular = "l";
				}
				else if(isCircular == CONTIG_CIRCULAR){
					linearOrCircular = "c";
				}
				else if(isCircular == CONTIG_CIRCULAR_RESCUED){
					linearOrCircular = "rc";
				}

				//_logFile << contigSequence.size() << endl;

				string header = ">ctg" + to_string(readIndex) + linearOrCircular + '\n';
				gzwrite(_toBasespace._basespaceContigFile, (const char*)&header[0], header.size());
				contigSequence +=  '\n';
				gzwrite(_toBasespace._basespaceContigFile, (const char*)&contigSequence[0], contigSequence.size());
			}




			_toBasespace._nbBps += contigSequence.size();

			_toBasespace._contigIndex += 1;

		}
	};
	

};	


#endif 

