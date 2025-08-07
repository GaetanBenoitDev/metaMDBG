
/*

- Collecting repetitive minimizers: desactiver completement en mode hifi (le sumbsample)
- avant de soumettre v2: supprimer le polishPartitionDir

- preparer le benchmark qui integre:
	- clipping
	- zero cov
	- contamination
	- besoin d'opti un peu le script clipping en mémoire

- si on fait le minimap2 polishing en -c, on peut utiliser le cigar fournit par minimap2
	- trop couteux a stocker en mémoire, autre alternative: au lieu de edlib on utilise un minimap2align

- checkpoint system: ajouter des checkpint dans le polisher

- todo circ check:
	- alignement sur le debut et la fin d'un contig
	- pas d'alignement sur d'autre contigs
	- on doit avoir pas mal de contig circulaire detruit vu qu'on enleve les repeats sur les bords des contigs dans le repeatremover sans checker si un read la bridge

- spoa polishing:
	- max hang ultra important
	- pour bien detecter le maxHang en ont, faire un alignement minimap2align (on peut aligner une fenetre du contig et du read un peu plus grande que le mapping)
	- attention le trimming des windows peut generer une courte sequence si on a selectionner une majorité de window incomplete par exemple, il faudrait faire evoluer le seuile de trimming tant que le consensus est trop court
	- on pourrait aussi empecher des window trop courte d'etre utiliser pendant la correction

- subsample for repetitive minimizers: essayer de predire un peu cb de reads il faut parser en fonction de la taille du jeu d'entrée

- quand on charche le successeur d'un read dans le graph, il est possible qu'on ne choississe pas le meilleur successeur (genre le read d'une autre strain), si quand on effecture un le minimap2 alignment, on se rend compte qu'il y a un petit overhang entre les deux reads (< 200), on pourrait essayer d'autrees successeurs

- read contigStart/End, fuat-il rajouter la length non alignée, pour tirer les alignement par ordre du contig aussi

- rafinage du bord des contigs

- mode --genomic:
	- discard all kminmer since once in first pass
	- desactiver le repeat remover
	- tip removing?


- ont read correction: moins bon resultat peut venir de la correction final via le chaining plutot que l'alignement (a test sur ZymoFecal)
- en mode ont on eneleve vraiment bcp de minimizer reptitive, pas sur de ça (check sur le soil aussi)

- RepeatRemover: baisser le 50 à 10 sourceUnitig
- RepeatRemover: si on a aucun source unitig (nbMinimizer > 50), il faudrait au moins prendre le plus grand par defaut, pour que chaque contig soit checké
- RepeatRemover: comment gerer l'endroit ou ont split les contig? (Actuellement on laisse un overlap mais il faudrait enlever l'overlap pour le fragment le plus petit), on peut aussi limité la taille des contigs a derepliquer (max 30000 par exemple, le but est d'enlever les erreurs)
- isCircular checker: readMapping + rafinement des bord

*/

/*

- transitive reudction: possible de l'accelerer en reactivant le code de filtrage basé sur la length
- sur le soil va falloir check la conso mémoire, check les potentiel memory leak
- _contigSequences: on peut le clear des qu'on a refait les alignment reads vs contigs (avant graphOverlapPath)

- si on valide k=15, on pourra ptet changer les minimizer pair en minimizer classique dans la module de read correction
- CreateMDBG: cout << "Opti: getSuccessor, getPredecessor: on peut arreter la collect des que ya plus de 1 successor/predecessor pendant l'unitigage" << endl;

./bin/metaMDBG polish --polish-target ~/appa/run/overlap_test_201/tmp/contigs_uncorrected.fasta.gz --in-hifi /pasteur/appa/scratch/gbenoit/data/overlap_test/genome_201_50x/simulatedReads_0.fastq.gz --out-dir ~/appa/run/overlap_test_201/ --threads 32 -n 50
- assemblage: paralleliser bubble, enlever tout le code pour les start/end node, nouveau parametre? k = 15, 0.005 fixer dans le readSelection a check

~/zeus/tools/merging/metaMDBG_7/metaMDBG/build/bin/metaMDBG polish asm/tmp/contigs_uncorrected.fasta.gz asm/contigs_2.fasta.gz asm/tmp/  --threads 32 -n 100 --metaMDBG --minimap2-options "-x map-ont"
python3 ~/zeus/scripts/nanopore/evaluation/computeReferenceCompleteness.py ~/appa/data/nanopore/mock/Zymo_D6331/references_filenames.txt asm/contigs.fasta.gz asm/contigs.fasta.gz metalol ~/appa/tmp/refComp_13 0.99 32
python3 ~/zeus/scripts/nanopore/evaluation/run_metaquast.py asm/contigs.fasta.gz ~/appa/data/nanopore/mock/Zymo_D6331/references_filenames.txt ~/zeus/tmp/metaquast_1/ 32
python3 ~/zeus/scripts/nanopore/eval_clipping.py asm/contigs_2.fasta.gz /pasteur/appa/scratch/gbenoit/data/nanopore/mock/Zymo_D6331/SRR17913199_Q20/SRR17913199.1.fastq.gz 32

/pasteur/appa/homes/gbenoit/zeus/tools/merging/metaMDBG_7/metaMDBG/build/bin/metaMDBG polish --polish-target ctg40713l_uncorrected.fasta --out-dir ./ctg40713l_polish/ --in-ont ctg40713l_reads.fastq --max-memory 17 --threads 32 -n 100 --metaMDBG

/pasteur/appa/homes/gbenoit/zeus/tools/merging/metaMDBG/build/bin/metaMDBG polish --polish-target ctg5583c_uncorrected.fasta --out-dir ./ctg5583c_polish/ --in-ont ctg5583c_reads.fastq --max-memory 17 --threads 32 -n 100 --metaMDBG
*/


#ifndef MDBG_METAG_CONTIGPOLISHER
#define MDBG_METAG_CONTIGPOLISHER

#include "../Commons.hpp"
#include "../utils/edlib.h"
#include "../utils/spoa/include/spoa/spoa.hpp"
#include "../utils/DnaBitset.hpp"
//#include "../readSelection/BaseAligner.hpp"
#include "ContigDerep.hpp"

class ContigPolisher : public Tool{
    
public:

	constexpr static const size_t _minimizerSize = 15;
	constexpr static const float _minimizerDensity = 0.1;
	constexpr static const int32_t _refineAlignEndsSize = 1000;
	constexpr static const float intFrac = 0.8;
	constexpr static const int64_t maxHang = 50;
	constexpr static const int64_t minOverlap = 500;
	constexpr static const int32_t _useAlignedLength = 200;
	constexpr static const u_int32_t _minContigLength = 2500;
	//constexpr static const bool _testAllVsAllAlignments = false;


	struct AlignmentBounds{

		public:

		int32_t _queryStart;
		int32_t _queryEnd;
		int32_t _referenceStart;
		int32_t _referenceEnd;
		int32_t _queryLength;
		int32_t _referenceLength;
		bool _isReversed;
		int32_t _originalChainingAlignmentLength;

		AlignmentBounds(){

			_queryStart = -1;
			_queryEnd = -1;
			_referenceStart = -1;
			_referenceEnd = -1;
			_queryLength = -1;
			_referenceLength = -1;
			_isReversed = false;
			_originalChainingAlignmentLength = 0;
		}

		void printAll() const{
			cout << _queryStart << " " << _queryEnd << "    " << _referenceStart << " " << _referenceEnd << " " << _isReversed << endl;
		}

	};


	struct OverlapGraphPath{
		u_int32_t _contigIndex;
		u_int32_t _pathIndex;
		u_int32_t _referenceStart;
		u_int32_t _referenceEnd;
		vector<u_int32_t> _path;
	};

	class OverlapGraph {
		
	public:

		unordered_map<ReadType, u_int32_t> _readIndex_to_nodeIndex;
		vector<ReadType> _nodeIndex_to_readIndex;
		
		struct Edge{
			u_int32_t _nodeIndex;
			//u_int32_t _nodeName;
			bool _isRemoved;
			//bool _reduce;
			//int32_t _overlapLength;
		};

		u_int32_t _contigIndex;
		vector<vector<Edge>> _successors;
		vector<vector<Edge>> _predecessors;

		OverlapGraph(const u_int32_t& contigIndex){
			_contigIndex = contigIndex;
		}

		void addEdge(const u_int32_t& fromNodeName, const u_int32_t& toNodeName){


			//cout << "\tAdd edge1: " << fromNodeName  << " -> " << toNodeName << endl;
			u_int32_t fromNodeIndex = getNodeIndex(fromNodeName);
			u_int32_t toNodeIndex = getNodeIndex(toNodeName);
			//cout << "\tAdd edge2: " << fromNodeIndex  << " -> " << toNodeIndex << endl;

			_successors[fromNodeIndex].push_back({toNodeIndex, false});
			_predecessors[toNodeIndex].push_back({fromNodeIndex, false});
			//cout << "\tAdd edge3" << endl;
		}

		u_int32_t getNodeIndex(const u_int32_t& readIndex){

			if(_readIndex_to_nodeIndex.find(readIndex) != _readIndex_to_nodeIndex.end()) return _readIndex_to_nodeIndex[readIndex];

			u_int32_t nodeIndex = _successors.size();
			_readIndex_to_nodeIndex[readIndex] = nodeIndex;
			_nodeIndex_to_readIndex.push_back(readIndex);

			_successors.push_back({});
			_predecessors.push_back({});

			//cout << "Add node: " << nodeIndex << " " << nodeName << endl;

			return nodeIndex;
		}

		void removeEdge(const u_int32_t& nodeIndexFrom, const u_int32_t& nodeIndexTo){

			for(Edge& successor : _successors[nodeIndexFrom]){
				if(successor._nodeIndex == nodeIndexTo){
					successor._isRemoved = true;
					break;
				}
			}

			for(Edge& predecessor : _predecessors[nodeIndexTo]){
				if(predecessor._nodeIndex == nodeIndexFrom){
					predecessor._isRemoved = true;
					break;
				}
			}

		}

	};

	struct Alignment{
		u_int32_t _contigIndex;
		u_int32_t _readIndex;
		bool _strand;
		u_int32_t _readStart;
		u_int32_t _readEnd;
		u_int32_t _contigStart;
		u_int32_t _contigEnd;
		float _identity;
		u_int32_t _readLength;
		u_int32_t _contigLength;
		//float _score;
		//u_int64_t _length;
		
		u_int64_t length() const{
			return std::max((u_int64_t)(_readEnd - _readStart), (u_int64_t)(_contigEnd - _contigStart));
		}

		float score() const{
			return (_contigEnd - _contigStart) * _identity;
		}

		float alignLength() const{
			return min(_readEnd - _readStart, _contigEnd - _contigStart);
		}

		void print() const{
			cout << _readIndex << " " << _identity << " " << _contigStart << " " << _contigEnd << " " << _strand << endl;
		}

		void printAll() const{
			cout << _readIndex << " " << _identity << "    " << _readStart << " " << _readEnd << "    " << _contigStart << " " << _contigEnd << " " << _strand << endl;
		}

	};

	/*
	struct Alignment2{
		u_int32_t _contigIndex;
		u_int32_t _readIndex;
		bool _strand;
		u_int32_t _readStart;
		u_int32_t _readEnd;
		u_int32_t _contigStart;
		u_int32_t _contigEnd;
		float _identity;
		//float _score;
		//u_int64_t _length;
		
		u_int64_t length() const{
			return std::max((u_int64_t)(_readEnd - _readStart), (u_int64_t)(_contigEnd - _contigStart));
		}

		float score() const{
			return (_contigEnd - _contigStart) * _identity;
		}

		float alignLength() const{
			return min(_readEnd - _readStart, _contigEnd - _contigStart);
		}

		void print() const{
			cout << _readIndex << " " << _identity << " " << _contigStart << " " << _contigEnd << " " << _strand << endl;
		}
	};
	*/


	string _inputFilename_reads;
	string _inputFilename_contigs;
	string _outputContigFilename;
	int _nbCores;
	size_t _windowLength;
	size_t _windowPositionOffset;
	size_t _windowLengthVariance;
	size_t _maxWindowCopies;
	//string _mapperOutputExeFilename;
	bool _useQual;
	string _outputDir;
	string _tmpDir;
	string _readPartitionDir;
	bool _circularize;
	//u_int64_t _minContigLength;
	
	string _outputFilename_mapping;
	string _outputMappingFilename_contigsVsUsedReads;
	int _maxMemoryGB;
	double _qualityThreshold;
	bool _useMetamdbgHeaderStyle;

	//gzFile _outputContigFile_polished;
	gzFile _outputContigFile_trimmed;
	bool _print_debug;
	//abpoa_para_t *abpt;
	


	//unordered_set<string> _isContigCircular;
	//unordered_set<u_int32_t> _isReadCircular;
	//unordered_map<u_int64_t, u_int64_t> _contigLength;
	//unordered_map<string, vector<string>> _contigMapLeft;
	//unordered_map<string, vector<string>> _contigMapRight;
	//unordered_map<u_int32_t, vector<u_int32_t>> _contigHitPos;
	unordered_map<u_int32_t, u_int32_t> _contigLength;
	unordered_map<u_int32_t, float> _contigCoverages;
	string _minimap2Preset_map;
	string _minimap2Preset_ava;
	bool _useHpc;
	u_int64_t _checksumWrittenReads;
	u_int64_t _outputContigIndex;

	//unordered_set<u_int32_t> _debugClipReads;

	struct ContigRead{
		u_int32_t _contigIndex;
		u_int64_t _readIndex;

		bool operator==(const ContigRead &other) const{
			return _contigIndex == other._contigIndex && _readIndex == other._readIndex;
		}

		size_t hash() const{
			std::size_t seed = 2;
			seed ^= _contigIndex + 0x9e3779b9 + (seed << 6) + (seed >> 2);
			seed ^= _readIndex + 0x9e3779b9 + (seed << 6) + (seed >> 2);
			return seed;
		}

	};

	struct ContigRead_hash {
		size_t operator()(const ContigRead& p) const
		{
			return p.hash();
		}
	};
	
	struct Window{
		DnaBitset2* _sequence;
		string _quality;
		u_int32_t _posStart;
		u_int32_t _posEnd;
		float _score;
	};

	ContigPolisher(): Tool (){

	}


	void parseArgs(int argc, char* argv[]){


		_print_debug = false;
		string ARG_NO_QUAL = "noqual";
		string ARG_IS_METAMDBG = "metaMDBG";
		string filenameExe = argv[0];
		//_logFile << filenameExe << endl;

		//fs::path pa(filenameExe);
		//_mapperOutputExeFilename = pa.parent_path().string() + "/mapper";
		//_logFile << _mapperOutputExeFilename << endl;
		//exit(1);


		/*
		args::ArgumentParser parser("polish", ""); //"This is a test program.", "This goes after the options."
		
		//args::ValueFlag<std::string> arg_outputDir(parser, "", "Output dir for contigs and temporary files", {ARG_OUTPUT_DIR2});
		//args::NargsValueFlag<std::string> arg_readFilenames_hifi(parser, "", "PacBio HiFi read filename(s) (separated by space)", {ARG_INPUT_HIFI}, 2);
		//args::NargsValueFlag<std::string> arg_readFilenames_nanopore(parser, "", "nanopore R10.4+ read filename(s) (separated by space)", {ARG_INPUT_NANOPORE}, 2);
		//args::ValueFlag<int> arg_nbCores(groupInputOutput, "", "Number of cores", {ARG_NB_CORES2}, NB_CORES_DEFAULT_INT);

		args::Positional<std::string> arg_contigs(parser, "contigs", "Contig filename to be corrected", args::Options::Required);
		args::Positional<std::string> arg_outputContigFilename(parser, "outputContigFilename", "Output contig filename", args::Options::Required);
		args::Positional<std::string> arg_outputDir(parser, "outputDir", "Output dir for contigs and temporary files", args::Options::Required);
		//args::PositionalList<std::string> arg_readFilenames(parser, "reads", "Read filename(s) used for correction (separated by space)", args::Options::Required);
		args::ValueFlag<int> arg_nbWindows(parser, "", "Maximum read coverage used for contig correction (increase for better correction)", {ARG_NB_WINDOWS}, 0);
		args::ValueFlag<string> arg_minimap2params(parser, "", "Set any minimap2 options (e.g. \"-x map-ont\")", {"minimap2-options"}, "");
		args::Flag arg_noQual(parser, "", "Do not use qualities during correction", {ARG_NO_QUAL});
		args::Flag arg_isMetaMDBG(parser, "", "Do not use qualities during correction", {ARG_IS_METAMDBG}, args::Options::Hidden);
		args::ValueFlag<int> arg_nbCores(parser, "", "Number of cores", {ARG_NB_CORES2}, NB_CORES_DEFAULT_INT);
		args::ValueFlag<int> arg_maxMemory(parser, "", "Maximum memory usage", {ARG_MAX_MEMORY}, 8);
		//args::Flag arg_useCirculize(parser, "", "Check if contigs are circular and add a flag in contig header (l: linear, c: circular)", {ARG_CIRCULARIZE});
		args::Flag arg_help(parser, "", "", {'h', "help"}, args::Options::Hidden);
		//args::HelpFlag help(parser, "help", "Display this help menu", {'h'});
		//args::CompletionFlag completion(parser, {"complete"});

		//(ARG_INPUT_FILENAME_TRUTH, "", cxxopts::value<string>()->default_value(""))
		//(ARG_MINIMIZER_LENGTH, "", cxxopts::value<int>()->default_value("13"))
		//(ARG_MINIMIZER_DENSITY, "", cxxopts::value<float>()->default_value("0.005"))
		//(ARG_NB_CORES, "", cxxopts::value<int>()->default_value(NB_CORES_DEFAULT))
		//(ARG_BLOOM_FILTER, "", cxxopts::value<bool>()->default_value("false"));
		*/



		args::ArgumentParser parser("polish", "Example: " + string(argv[0]) + " polish --polish-target ./contigs.fasta.gz --in-ont reads.fastq.gz --out-dir ./outputDir/ --threads 4"); //"This is a test program.", "This goes after the options."
		args::Group groupInputOutput(parser, "Basic options:");
		args::Group groupOther(parser, "Other options:");

		
		args::ValueFlag<std::string> arg_outputDir(groupInputOutput, "", "Output dir for contigs and temporary files", {ARG_OUTPUT_DIR2});
		args::NargsValueFlag<std::string> arg_readFilenames_hifi(groupInputOutput, "", "PacBio HiFi read filename(s) (separated by space)", {ARG_INPUT_HIFI}, 2);
		args::NargsValueFlag<std::string> arg_readFilenames_nanopore(groupInputOutput, "", "Nanopore R10.4+ read filename(s) (separated by space)", {ARG_INPUT_NANOPORE}, 2);
		args::ValueFlag<std::string> arg_polishTarget(groupInputOutput, "", "Contigs to polish", {ARG_POLISH_TARGET});
		args::ValueFlag<int> arg_nbCores(groupInputOutput, "", "Number of cores", {ARG_NB_CORES2}, NB_CORES_DEFAULT_INT);

		args::ValueFlag<int> arg_nbWindows(groupOther, "", "Maximum read coverage used for contig correction (increase for better correction)", {ARG_NB_WINDOWS}, 0);
		args::Flag arg_isMetaMDBG(groupOther, "", "Use metaMDBG style headers", {ARG_IS_METAMDBG}, args::Options::Hidden);
		args::ValueFlag<int> arg_maxMemory(groupOther, "", "Maximum memory usage for read mapping", {ARG_MAX_MEMORY}, 8);

		//args::ValueFlag<int> arg_minimizerSize(groupAssembly, "", "k-mer size", {ARG_MINIMIZER_LENGTH2}, 13);
		//args::ValueFlag<float> arg_densityAssembly(groupAssembly, "", "Fraction of total k-mers used for assembly", {ARG_MINIMIZER_DENSITY_ASSEMBLY}, 0.005f);
		//args::ValueFlag<int> arg_maxK(groupAssembly, "", "Stop assembly after k iterations", {ARG_MAXK}, 0);
		
		args::Flag arg_help(parser, "", "", {'h', "help"}, args::Options::Hidden);


		try
		{
			parser.ParseCLI(argc, argv);
		}
		catch (const args::Help&)
		{
			cerr << parser;
			exit(0);
		}
		catch (const std::exception& e)
		{
			cerr << parser;
			//_logFile << endl;
			cerr << e.what() << endl;
			exit(0);
		}

		if(arg_help){
			cerr << parser;
			exit(0);
		}

		Commons::checkRequiredArgs(parser, arg_outputDir, arg_readFilenames_hifi, arg_readFilenames_nanopore);

		_inputFilename_contigs = args::get(arg_polishTarget);
		//_outputContigFilename = args::get(arg_outputContigFilename);
		_outputDir = args::get(arg_outputDir);
		_nbCores = args::get(arg_nbCores);
		//_minimap2Params = args::get(arg_minimap2params);
		
		_maxMemoryGB = args::get(arg_maxMemory);
		_maxMemoryGB = max(4, _maxMemoryGB);
		//_maxMemory = maxMemoryGB * 1000000000ull;

		//cout << _maxMemoryGB << endl;

		_useQual = true;
		//if(arg_noQual){
		//	_useQual = false;
		//}

		_useMetamdbgHeaderStyle = false;
		if(arg_isMetaMDBG){
			_useMetamdbgHeaderStyle = true;
		}

		_circularize = false;
		//if(arg_useCirculize){
		//	_circularize = true;
		//}



		_windowLength = 500;
		//_windowPositionOffset = 0;
		_windowLengthVariance = _windowLength*0.01;
		_maxWindowCopies = args::get(arg_nbWindows);
		_qualityThreshold = 10.0;
		//_minContigLength = 0;

		_tmpDir = _outputDir;// + "/tmp/";
		//if(!fs::exists(_tmpDir)){
		//	fs::create_directories(_tmpDir);
		//}

		_readPartitionDir = _tmpDir + "/_polish_readPartitions/";
		_inputFilename_reads = _readPartitionDir + "/input.txt";
		_outputContigFilename = _tmpDir + "/contigs.fasta.gz";
		_usedReadFilename = _readPartitionDir + "/usedReads.fasta.gz";

		if(!fs::exists(_tmpDir)){
			fs::create_directories(_tmpDir);
		}

		if(!fs::exists(_readPartitionDir)){
			fs::create_directories(_readPartitionDir);
		}

		if(arg_readFilenames_hifi){
			Commons::createInputFile(args::get(arg_readFilenames_hifi), _inputFilename_reads);
			_minimap2Preset_map = "map-hifi";
			_minimap2Preset_ava = "ava-pb";
			_useHpc = true;
			//_minimap2Params = "-x map-hifi";
			//_params = getParamsHifi();
		}
		else if(arg_readFilenames_nanopore){
			Commons::createInputFile(args::get(arg_readFilenames_nanopore), _inputFilename_reads);
			_minimap2Preset_map = "map-ont";
			_minimap2Preset_ava = "ava-ont";
			_useHpc = false;
			//_params = getParamsNanopore();
		}

		/*
		if(arg_readFilenames_hifi){
			Commons::createInputFile(args::get(arg_readFilenames_hifi), _inputFilename_reads);
		}
		else if(arg_readFilenames_nanopore){
			Commons::createInputFile(args::get(arg_readFilenames_nanopore), _inputFilename_reads);
		}
		*/

		//Commons::createInputFile(args::get(arg_readFilenames), _inputFilename_reads);

		openLogFile(_tmpDir);


		Logger::get().debug() << "";
		Logger::get().debug() << "Contigs: " << _inputFilename_contigs;
		Logger::get().debug() << "Reads: " << _inputFilename_reads;
		Logger::get().debug() << "Use quality: " << _useQual ;
		Logger::get().debug() << "Max window variants: " << _maxWindowCopies ;
		Logger::get().debug() << "";

		//fs::path p(_inputFilename_contigs);
		//while(p.has_extension()){
		//	p.replace_extension("");
		//}



		//_outputFilename_contigs = p.string() + "_corrected.fasta.gz";
		//_outputFilename_mapping = p.string() + "_tmp_mapping__.paf";
		_outputFilename_mapping = _readPartitionDir + "/polish_mapping.paf.gz";
		//_outputFilename_mapping = "/mnt/gpfs/gaetan/tmp/debug_polish/racon_align.paf";

		//_maxMemory = 4000000000ull;
		//_maxMemory = 2000000ull;
		//_nbCores = 2;


		_fields = new vector<string>();
		_fields_optional = new vector<string>();


	}


	u_int64_t _currentContigSize;
	//u_int64_t _windowByteSize;
	//unordered_set<u_int32_t> _validContigIndexes;
	
	string _usedReadFilename;
	gzFile _usedReadFile;
	gzFile _outputContigFilePartition;
	gzFile _outputContigFileRepeat;
	ofstream _file_contigHeaders;

	/*
	abpt = abpoa_init_para();
	abpt->out_msa = 0; // generate Row-Column multiple sequence alignment(RC-MSA), set 0 to disable
	abpt->out_cons = 1; // generate consensus sequence, set 0 to disable
	abpt->w = 6, abpt->k = 9; abpt->min_w = 10; // minimizer-based seeding and partition
	abpt->progressive_poa = 1;
	abpt->max_n_cons = 1;

	abpoa_post_set_para(abpt);
	*/
		


    void execute (){

		//for(size_t i=0; i<40; i++){
			//_partitionNbReads[i] = 10000;
		//	processPartition(i);
		//}
		//exit(1);
		
		//for(size_t i=2; i<50; i++){
		//	_partitionNbReads[i] = 10000;
		//	detectWeakRepeats2(i);
		//}
		//exit(1);
		//detectWeakRepeats2(9);
		//exit(1);
		//cout << "Load kminmer abundance" << endl;
		//loadKminmerAbundance();

		//cout << "Loading unitig index" << endl;
		//loadUnitigIndex();

		//cout << "Loading contig data" << endl;
		//loadContigData();
		/*
		cout << "todo: remove clip file debug" << endl;
		ifstream infile("/pasteur/appa/homes/gbenoit/appa/run/correction/test_deterministic/test_zymo_1/clipReads.txt");
		std::string line;
		while (std::getline(infile, line))
		{
			std::istringstream iss(line);
			int a;
			if (!(iss >> a)) { break; } // error

			cout << a << endl;
			_debugClipReads.insert(a);
			// process pair (a,b)
		}
		*/


		int minimapBatchSize = 1;
		//float peakMemory = getPeakMemory();
		if(_maxMemoryGB < 8 || _maxMemoryGB > 1000000){
			minimapBatchSize = 1;
		}
		else{
			minimapBatchSize = _maxMemoryGB / 8;
		}
		if(minimapBatchSize < 0){
			minimapBatchSize = 1;
		}
		if(minimapBatchSize > 100){
			minimapBatchSize = 100;
		}



		//trimContigs(minimapBatchSize);
		//exit(1);

		//const string& inputContigFilenameDerep = _readPartitionDir + "/contigs_uncorrected_derep.fasta.gz";
		//derepContigs(_inputFilename_contigs, inputContigFilenameDerep, minimapBatchSize);
		
		//_inputFilename_contigs = inputContigFilenameDerep;
		
		
		Logger::get().debug() << "Map reads to all contigs";
		string command = "minimap2 -v 0 -m 500 -I " + to_string(minimapBatchSize) + "G -t " + to_string(_nbCores) + " -x " + _minimap2Preset_map + " " + _inputFilename_contigs + " " + Commons::inputFileToFilenames(_inputFilename_reads);
		Utils::executeMinimap2(command, _outputFilename_mapping);
		

		//int error = Utils::pipeCommands(command1, command2, _readPartitionDir + "/polish_mapping.paf.gz");
		//cout << "Return: " << error << endl;
		
		//mapReads();
		
		
		Logger::get().debug() << "";
		Logger::get().debug() << "Partitionning reads on the disk";
		auto start = high_resolution_clock::now();
		
		indexContigName(_inputFilename_contigs);
		indexReadName();
		//if(_circularize) executeCircularize();
		
		loadAlignments(_outputFilename_mapping, false, false, false);

		_contigName_to_contigIndex.clear();
		_readName_to_readIndex.clear();

		computeContigCoverages(_inputFilename_contigs);
		//writeAlignmentBestHits();
		//parseAlignments(true, false);
		partitionReads();
		Logger::get().debug() << "\tDone " << " (" << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s)";
		
		
		//_nbPartitions = _nbCores;



		Logger::get().debug() << "";
		Logger::get().debug() << "Curating contigs";

		//_outputContigFile = gzopen(_outputFilename_contigs.c_str(),"wb");
		_usedReadFile = gzopen(_usedReadFilename.c_str(),"wb");

		//_partitionNbReads[30] = 10000;
		//processPartition(30);
		//exit(1);

		

		for(size_t i=0; i<_nbPartitions; i++){
			//_partitionNbReads[i] = 10000;
			
			//if(i != 30) continue;

			//processPartition(i);
		}
		

		gzclose(_usedReadFile);
		

		//for(size_t i=0; i<_nbPartitions; i++){
		//	_partitionNbReads[i] = 10000;
			//detectWeakRepeats2(i);
		//}



		Logger::get().debug() << "";
		Logger::get().debug() << "Polishing contigs (pass 1)";
		_outputContigIndex = 0;
		
		//_nbPartitions = 50;

		for(size_t i=0; i<_nbPartitions; i++){

			//if(i != 30) continue;

			//_partitionNbReads[i] = 10000;

			//string inputContigFilename = _readPartitionDir + "/" + to_string(i) + "_contigsCurated.gz";
			string inputContigFilename = _readPartitionDir + "/" + to_string(i) + "_contigs.gz";
			string outputContigFilename = _readPartitionDir + "/" + to_string(i) + "_contigsPolished.gz";
			gzFile outputContigFile = gzopen(outputContigFilename.c_str(),"wb");

			//_partitionNbReads[i] = 10000;
			polishPartition(i, inputContigFilename, outputContigFile, 0);

			gzclose(outputContigFile);
		}

		//cout << "done" << endl;



		Logger::get().debug() << "";
		Logger::get().debug() << "Polishing contigs (pass 2)";
		_outputContigIndex = 0;

		const string& outputContigFilename_polished = _readPartitionDir + "/contigs_polished.fasta.gz";
		gzFile outputContigFile_polished = gzopen(outputContigFilename_polished.c_str(),"wb");
		_file_contigHeaders = ofstream(_tmpDir + "/contigHeaders.txt");

		for(size_t i=0; i<_nbPartitions; i++){


			//if(i != 30) continue;

			string inputContigFilename = _readPartitionDir + "/" + to_string(i) + "_contigsPolished.gz";

			polishPartition(i, inputContigFilename, outputContigFile_polished, 1);

		}

		_file_contigHeaders.close();
		gzclose(outputContigFile_polished);
		
		//exit(1);
		//indexContigName();
		//indexReadName();

		/*
		u_int64_t totalByteSize = 0;

		u_int64_t contigIndex = 0;
		gzFile fp;
		kseq_t *seq;
		int slen = 0, qlen = 0;
		fp = gzopen(_inputFilename_contigs.c_str(), "r");
		seq = kseq_init(fp);

		while (kseq_read(seq) >= 0){

			totalByteSize += loadContig(contigIndex, string(seq->seq.s));
			//_logFile << totalByteSize << " " << _windowByteSize << endl;
			if(totalByteSize >= _maxMemory){
				processPass();
				clearPass();
				totalByteSize = 0;
			}

			contigIndex += 1;
		}
		*/
			
		//processPass();

		/* !!!!!!!!!!!!
		kseq_destroy(seq);
		gzclose(fp);
		*/
		

		Logger::get().debug() << "";
		Logger::get().debug() << "Dereplicating contigs";
		//const string& derepContigFilename1 = _tmpDir + "/contigs_derep_1.fasta.gz";
		const string& outputContigFilename_derep = _readPartitionDir + "/contigs_derep.fasta.gz";
		derepContigs(outputContigFilename_polished, outputContigFilename_derep, minimapBatchSize);

		Logger::get().debug() << "";
		Logger::get().debug() << "Dereplicating contigs";
		//const string& derepContigFilename1 = _tmpDir + "/contigs_derep_1.fasta.gz";
		const string& outputContigFilename_derep_2 = _readPartitionDir + "/contigs_derep_2.fasta.gz";
		derepContigs(outputContigFilename_derep, outputContigFilename_derep_2, minimapBatchSize);

		Logger::get().debug() << "";
		Logger::get().debug() << "Trimming contigs";
		const string& outputContigFilename_trimmed = _readPartitionDir + "/contigs_trimmed.fasta.gz";
		//trimContigs(outputContigFilename_derep, outputContigFilename_trimmed, minimapBatchSize);


		Logger::get().debug() << "";
		Logger::get().debug() << "Moving final contigs to destination";
		//fs::rename(outputContigFilename_trimmed, _outputContigFilename);
		fs::rename(outputContigFilename_derep_2, _outputContigFilename);


		cout << "A REMETTRE: delete tmp dir for polishing" << endl;
		fs::remove_all(_readPartitionDir);
		
		Logger::get().debug() << "";
		Logger::get().debug() << "Polished contigs: " << _outputContigFilename;
		Logger::get().debug() << "done";
		//closeLogFile();
	}


	void trimContigs(const string& inputContigFilename, const string& outputContigFilename, const int& minimapBatchSize){
		
		int minimapBatchSizeC = max(1, minimapBatchSize / 4);

		auto start = high_resolution_clock::now();
		Logger::get().debug() << "\tMap used reads to final contigs";
		
		const string& alignFilename = _readPartitionDir + "/align_contigsVsUsedReads.paf.gz";
		string command1 = "minimap2 -v 0 -p 1 -c -I " + to_string(minimapBatchSizeC) + "G -t " + to_string(_nbCores) + " -x " + _minimap2Preset_map + " " + inputContigFilename + " " + _usedReadFilename;
		//cout << _outputFilename_mapping << endl;
		Utils::executeMinimap2(command1, alignFilename);

		indexContigName(inputContigFilename);
		loadAlignments(alignFilename, true, true, false);

		Logger::get().debug() << "\tTrimming contigs";
		_outputContigFile_trimmed = gzopen(outputContigFilename.c_str(),"wb");

		ReadParserParallel readParser(inputContigFilename, true, false, 1); //single core
		readParser.parse(ContigTrimmerFunctor(*this));

		gzclose(_outputContigFile_trimmed);

		Logger::get().debug() << "\tDone " << " (" << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s)";
	}


	class ContigTrimmerFunctor {

		public:

		ContigPolisher& _parent;


		ContigTrimmerFunctor(ContigPolisher& parent) : _parent(parent){
		}

		ContigTrimmerFunctor(const ContigTrimmerFunctor& copy) : _parent(copy._parent){
		}

		~ContigTrimmerFunctor(){
		}

		void operator () (const Read& read) {

			u_int32_t contigIndex = read._index;
			if(read._index % 10000 == 0) Logger::get().debug() << "\t" << read._index;
			
			const string& contigSequence = read._seq;
			vector<bool> isContigCovered(contigSequence.size(), false);

			if(_parent._alignments_readsVsContigs.find(contigIndex) == _parent._alignments_readsVsContigs.end()) return;

			auto& alignments = _parent._alignments_readsVsContigs[contigIndex];

			for(const auto& it : alignments){

				const ReadType& readIndex = it.first;
				const Alignment& alignment = it.second;

				for(size_t i=alignment._contigStart; i <= alignment._contigEnd; i++){
					isContigCovered[i] = true;
				}
			}

			u_int32_t startLengthToRemove = 0;
			for(int64_t i=0; i<isContigCovered.size(); i++){
				startLengthToRemove = i;
				if(isContigCovered[i]) break;
			}

			u_int32_t endLengthToRemove = 0;
			for(int64_t i = isContigCovered.size()-1; i >= 0; i--){
				endLengthToRemove = (isContigCovered.size() - i - 1);
				if(isContigCovered[i]) break;
			}

			//for(int64_t i=0; i<isContigCovered.size(); i++){
			//	cout << isContigCovered[i] << " ";
			//}
			//cout << endl;

			//cout << startLengthToRemove << " " << endLengthToRemove << endl;

			if(startLengthToRemove + endLengthToRemove >= contigSequence.size()) return;

			string contigSequenceTrimmed = contigSequence;

			if(startLengthToRemove > 0){
				contigSequenceTrimmed = contigSequenceTrimmed.substr(startLengthToRemove, contigSequenceTrimmed.size()-startLengthToRemove);
			}

			if(endLengthToRemove > 0){
				contigSequenceTrimmed = contigSequenceTrimmed.substr(0, contigSequenceTrimmed.size()-endLengthToRemove);
			}

			if(contigSequenceTrimmed.size() < ContigPolisher::_minContigLength) return;

			string header = read._header;

			if(_parent._useMetamdbgHeaderStyle){
				Utils::ContigHeader contigHeader = Utils::extractContigHeader(header);
				header = Utils::createContigHeader(contigHeader._contigIndex, contigSequenceTrimmed.size(), contigHeader._coverage, contigHeader._isCircular);
			}

			string headerFasta = ">" + header + '\n';
			gzwrite(_parent._outputContigFile_trimmed, (const char*)&headerFasta[0], headerFasta.size());
			contigSequenceTrimmed +=  '\n';
			gzwrite(_parent._outputContigFile_trimmed, (const char*)&contigSequenceTrimmed[0], contigSequenceTrimmed.size());

		}
	};

	
	void derepContigs(const string& inputContigFilename, const string& outputContigFilename, const int& minimapBatchSize){

		auto start = high_resolution_clock::now();

		ContigDerep contigDerep(inputContigFilename, outputContigFilename, _useMetamdbgHeaderStyle, minimapBatchSize, _tmpDir, _readPartitionDir, 0.9, _minContigLength, _nbCores);
		contigDerep.execute();
		
		Logger::get().debug() << "\tDone " << " (" << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s)";
		/*
		return;
		const string derepTmp = _readPartitionDir + "/derep.fasta.gz";
		{
			ContigDerep contigDerep(inputContigFilename, derepTmp, to_string(minimapBatchSize), _tmpDir, _readPartitionDir, 0.9, _nbCores);
			contigDerep.execute();
		}

		{
			ContigDerep contigDerep(derepTmp, outputContigFilename, to_string(minimapBatchSize), _tmpDir, _readPartitionDir, 0.9, _nbCores);
			contigDerep.execute();
		}
		*/
	}
	/*
	void executeCircularize(){

		cerr << "Detecting circular contigs..." << endl;

		indexCircularContigs();
		indexMappingOnContigEnds();
		detectCircularContigs();

		_contigMapRight.clear();
		_contigMapLeft.clear();
		//_contigLength.clear();
		//_isContigCircular.clear();
	}

	void indexCircularContigs(){

		auto fp = std::bind(&ContigPolisher::indexCircularContigs_read, this, std::placeholders::_1);
		ReadParser readParser(_inputFilename_contigs, true, false, _logFile);
		readParser.parse(fp);
    
	}

	void indexCircularContigs_read(const Read& read){

		//u_int64_t contigIndex = Utils::contigName_to_contigIndex(read._header);

		//_logFile << "Contig: " << read._index << " " << read._seq.size() << endl;

		if (read._seq.size() < _minContigLength) return;

		//u_int32_t contigIndex = Utils::contigName_to_contigIndex(read._header);
		//_contigLength[read._index] = read._seq.size();

		char lastChar = read._header[read._header.size()-1];
		if(lastChar == 'c'){
			_isContigCircular.insert(Utils::shortenHeader(read._header));
		}

	}

	void indexMappingOnContigEnds(){

		cerr << "\tParsing alignments" << endl;

		_alignments.clear();
		//_indexPartitionOnly = indexPartitionOnly;

		PafParser pafParser(_outputFilename_mapping);
		auto fp = std::bind(&ContigPolisher::indexMappingOnContigEnds_read, this, std::placeholders::_1);
		pafParser.parse(fp);

	}

	void indexMappingOnContigEnds_read(const string& line){

    	long maxHang = 100;

		double errorThreshold = 0.3;

		GfaParser::tokenize(line, _fields, '\t');

		//_logFile << line << endl;

		const string& readName = Utils::shortenHeader((*_fields)[0]);
		const string& contigName = Utils::shortenHeader((*_fields)[5]);
		

		u_int32_t readStart = stoull((*_fields)[2]);
		u_int32_t readEnd = stoull((*_fields)[3]);
		u_int32_t contigLength = stoull((*_fields)[6]);
		u_int32_t contigStart = stoull((*_fields)[7]);
		u_int32_t contigEnd = stoull((*_fields)[8]);

		u_int64_t nbMatches = stoull((*_fields)[9]);
		u_int64_t alignLength = stoull((*_fields)[10]);
		u_int64_t queryLength = stoull((*_fields)[1]);

		bool strand = (*_fields)[4] == "-";
		//float score = (double) nbMatches / (double) queryLength;
		//float score2 = (double) nbMatches / (double) alignLength;

		//u_int32_t contigIndex = mappingIndex._contigName_to_contigIndex[contigName];
		//u_int64_t readIndex = mappingIndex._readName_to_readIndex[readName];
		//u_int64_t length = std::max(readEnd - readStart, contigEnd - contigStart);

		if(alignLength < 3000) return;
		if (contigLength < _minContigLength) return;
		//_logFile << contigIndex << " " << readIndex << endl;

		double length = max((double)(contigEnd - contigStart), (double)(readEnd - readStart));
		double error = 1 - min((double)(contigEnd - contigStart), (double)(readEnd - readStart)) / length;
		//_logFile << error << " " << errorThreshold << endl;
		
		if(error > errorThreshold) return;




		long hangLeft = contigStart;
		long hangRight = contigLength - contigEnd;

		if (hangLeft < maxHang){
			_contigMapLeft[contigName].push_back(readName);
			//_logFile << lineInput << endl;
		}
		if (hangRight < maxHang){
			_contigMapRight[contigName].push_back(readName);
			//_logFile << lineInput << endl;
		}

	}

	
	void detectCircularContigs(){

		auto fp = std::bind(&ContigPolisher::detectCircularContigs_read, this, std::placeholders::_1);
		ReadParser readParser(_inputFilename_contigs, true, false, _logFile);
		readParser.parse(fp);
    
	}

	void detectCircularContigs_read(const Read& read){

		const string& contigName = Utils::shortenHeader(read._header);

		//u_int64_t contigIndex = read._index; //Utils::contigName_to_contigIndex(read._header);

		u_int8_t isCircular = CONTIG_LINEAR
		u_int64_t nbSupportingReads = countCircularReads(contigName);

		//char lastChar = read._header[read._header.size()-1];




		//string header = read._header;
		//header.pop_back();
		
		//if(isCircular){
		//	header += "c";
			//_isContigCircular[read._index] = true;
		//}
		//else{
			//_isContigCircular[read._index] = false;
		//	header += "l";
		//}

		//header = ">" + header + '\n';
		//gzwrite(_outputContigFile, (const char*)&header[0], header.size());
		//string contigSequence = read._seq + '\n';
		//gzwrite(_outputContigFile, (const char*)&contigSequence[0], contigSequence.size());

		if(read._seq.size() >= _minContigLength){
			_logFile << read._header << " " << (_isContigCircular.find(contigName) != _isContigCircular.end()) << " " << nbSupportingReads << endl;
		}

		if(_isContigCircular.find(contigName) != _isContigCircular.end()){
			//isCircular = true;
		}
		else{
			if(nbSupportingReads >= 3){
				//isCircular = true;
				_isContigCircular.insert(contigName);
			}
		}

	}
	



	u_int64_t countCircularReads(const string& contigIndex){

		if((_contigMapLeft.find(contigIndex) == _contigMapLeft.end()) || (_contigMapRight.find(contigIndex) == _contigMapRight.end())) return 0;

		//_logFile << "---" << endl;
		u_int64_t nbSupportingReads = 0;


		const vector<string>& leftReads = _contigMapLeft[contigIndex];
		const vector<string>& rightReads = _contigMapRight[contigIndex];

		//_logFile << leftReads.size() << endl;
		//_logFile << rightReads.size() << endl;

		for(const string& readIndex : leftReads){
			//_logFile << readIndex << endl;
			if(std::find(rightReads.begin(), rightReads.end(), readIndex) != rightReads.end()){
				nbSupportingReads += 1;
				//_logFile << readIndex << endl;
				//_isReadCircular.insert(readIndex);
			}
		}

		return nbSupportingReads;


	}
	*/

	int _contigCoverageWindow = 100;

	void computeContigCoverages(const string& contigFilename){
		
		_contigLength.clear();
		_contigCoverages.clear();

		Logger::get().debug() << "\tComputing contig coverages";

		unordered_map<u_int32_t, vector<pair<u_int32_t, u_int32_t>>> contigHits;
		auto fp = std::bind(&ContigPolisher::computeContigCoverages_setup_read, this, std::placeholders::_1);
		ReadParser readParser(contigFilename, true, false);
		readParser.parse(fp);


		for(const auto& it : _alignments){

			//const string& contigName = it.second._contigName;
			u_int32_t contigIndex = it.second._contigIndex;

			//if(contigHits.find(contigIndex) == contigHits.end()){
			//	contigHits[contigIndex] = {};
			//}

			contigHits[contigIndex].push_back({it.second._contigStart, it.second._contigEnd});

			//for(size_t i=it.second._contigStart; i<it.second._contigEnd; i++){
			//	if(i%_contigCoverageWindow != 0) continue;
				//cout << i/_contigCoverageWindow << " " << contigHitPos.size() << endl;
			//	if(i/_contigCoverageWindow >= contigHitPos.size()) continue;
			//	contigHitPos[i/_contigCoverageWindow] += 1;
				//cout << i/_contigCoverageWindow << " " << contigHitPos.size() << endl;
			//}

		}

		for(auto& it : contigHits){

			//cout << it.first << " " << _contigLength[it.first] << endl;
			vector<u_int32_t> coverages(_contigLength[it.first], 0);

			for(auto& interval : it.second){

				for(size_t i=interval.first; i<interval.second; i++){
					coverages[i] += 1;
				}
			}
			//cout << it.first << endl;
			//for(u_int32_t count : it.second){
			//	cout << count << " ";
			//}
			//cout << endl;

			if (coverages.size() < 80 *2){
				_contigCoverages[it.first] = 1;
			}
			else{
				u_int64_t sum = 0;
				for(long i=75; i<((long)coverages.size())-75; i++){
					sum += coverages[i];
				}

				float coverage = sum / ((double)coverages.size());  //Utils::compute_median(coverages);
				_contigCoverages[it.first] = coverage;
			}


			//Logger::get().debug() << "Coverage " << it.first << ": " << coverage << endl;
		}

		//_contigHitPos.clear();
		_contigLength.clear();
	}

	
	void computeContigCoverages_setup_read(const Read& read){

		_contigLength[read._index] = read._seq.size();
		//int nbCounts = read._seq.size() / _contigCoverageWindow;
		//cout << Utils::shortenHeader(read._header) << " " << nbCounts << endl;
		//_contigHitPos[read._index].resize(nbCounts, 0);
	}

	/*
	void computeContigCoverages(){
		
		Logger::get().debug() << "Computing contig coverages...";

		auto fp = std::bind(&ContigPolisher::computeContigCoverages_setup_read, this, std::placeholders::_1);
		ReadParser readParser(_inputFilename_contigs, true, false);
		readParser.parse(fp);


		for(const auto& it : _alignments){

			//const string& contigName = it.second._contigName;
			u_int32_t contigIndex = it.second._contigIndex;
			vector<u_int32_t>& contigHitPos = _contigHitPos[contigIndex];

			for(size_t i=it.second._contigStart; i<it.second._contigEnd; i++){
				if(i%_contigCoverageWindow != 0) continue;
				//cout << i/_contigCoverageWindow << " " << contigHitPos.size() << endl;
				if(i/_contigCoverageWindow >= contigHitPos.size()) continue;
				contigHitPos[i/_contigCoverageWindow] += 1;
				//cout << i/_contigCoverageWindow << " " << contigHitPos.size() << endl;
			}

		}

		for(auto& it : _contigHitPos){
			//cout << it.first << endl;
			//for(u_int32_t count : it.second){
			//	cout << count << " ";
			//}
			//cout << endl;

			
			float coverage = Utils::compute_median(it.second);
			_contigCoverages[it.first] = coverage;

			//Logger::get().debug() << "Coverage " << it.first << ": " << coverage << endl;
		}

		_contigHitPos.clear();
	}


	void computeContigCoverages_setup_read(const Read& read){

		int nbCounts = read._seq.size() / _contigCoverageWindow;
		//cout << Utils::shortenHeader(read._header) << " " << nbCounts << endl;
		_contigHitPos[read._index].resize(nbCounts, 0);
	}
	*/
	/*
	u_int64_t loadContig(u_int64_t contigIndex, const string& seq){

		//_validContigIndexes.insert(contigIndex);

		u_int64_t totalByteSize = 0;

		_contigSequences[contigIndex] = seq;

		//_logFile << "load contig: " << contigIndex << " " << seq.size() << endl;
		size_t nbWindows = ceil((double)seq.size() / (double)_windowLength);
		vector<vector<Window>> windows(nbWindows);
		//_logFile << "Nb windows: " << nbWindows << endl;

		_contigWindowSequences[contigIndex] = windows;

		//_logFile << nbWindows << " " << nbWindows * _windowByteSize << endl;
		return nbWindows * _windowByteSize;
	}
	*/
	
	class PartitionFile{

		public:

		gzFile _file;
		omp_lock_t _mutex;

		PartitionFile(u_int32_t partition, const string& partitionFilename){
			_file = gzopen(partitionFilename.c_str(), "wb");

			omp_init_lock(&_mutex);
		}

		~PartitionFile(){
			gzclose(_file);
			omp_destroy_lock(&_mutex);
		}
	};

	unordered_map<u_int32_t, u_int32_t> _contigToPartition;
	vector<PartitionFile*> _partitionFiles;
	unordered_map<u_int32_t, u_int64_t> _partitionNbReads;

	
	struct ContigStats{
		u_int32_t _contigIndex;
		u_int32_t _length;
	};

	vector<ContigStats> _contigStats;
	u_int32_t _nbPartitions;


	void partitionReads(){

		collectContigStats();

		int nbPartitionsInit = max((u_int32_t)1, (u_int32_t)_nbCores);
		u_int64_t totalByteSize = 0;

		vector<u_int64_t> memoryPerPartitions(nbPartitionsInit, 0);

		u_int64_t maxMemory = 4000000000ull;

		for(const ContigStats& contigStat : _contigStats){

			int partitionIndex = 0;
			u_int64_t minMemory = -1;

			for(size_t i=0; i<memoryPerPartitions.size(); i++){
				if(memoryPerPartitions[i] < minMemory){
					minMemory = memoryPerPartitions[i];
					partitionIndex = i;
				}
			}

			u_int64_t currentParitionMemory = memoryPerPartitions[partitionIndex];
			u_int64_t contigMemory = computeContigMemory(contigStat._contigIndex, contigStat._length);

			//cout << "lala " << contigMemory << endl;
			
			if(currentParitionMemory > 0 && currentParitionMemory+contigMemory > maxMemory){
				memoryPerPartitions.push_back(0);
				partitionIndex = memoryPerPartitions.size()-1;
			}
			
			memoryPerPartitions[partitionIndex] += contigMemory;
			_contigToPartition[contigStat._contigIndex] = partitionIndex;
			

		}



		_nbPartitions = 0;
		for(size_t i=0; i<memoryPerPartitions.size(); i++){
			if(memoryPerPartitions[i] > 0){
				_nbPartitions += 1;
			}
		}
		

		writeReadPartitions();
		writeContigPartitions();

		_contigStats.clear();
		_contigCoverages.clear();
		//cout << "_contigCoverages.clear(); a remettre" << endl;
	}

	void writeReadPartitions(){

		for(u_int32_t i=0; i<_nbPartitions; i++){
			_partitionNbReads[i] = 0;
			_partitionFiles.push_back(new PartitionFile(i, _readPartitionDir + "/" + to_string(i) + "_reads.gz"));
		}


		//_contigCoverages.clear();

		Logger::get().debug() << "\tNb partitions: " << _nbPartitions;
		
		_checksumWrittenReads = 0;
		ReadParserParallel readParser(_inputFilename_reads, false, false, _nbCores);
		readParser.parse(ReadPartitionFunctor(*this));

		for(PartitionFile* partitionFile : _partitionFiles){
			delete partitionFile;
		}
		_partitionFiles.clear();

	}

	void writeContigPartitions(){

		for(u_int32_t i=0; i<_nbPartitions; i++){
			//_partitionNbReads[i] = 0;
			//_partitionFiles.push_back(new PartitionFile(i, _readPartitionDir));
			_partitionFiles.push_back(new PartitionFile(i, _readPartitionDir + "/" + to_string(i) + "_contigs.gz"));
		}

		ReadParserParallel readParser(_inputFilename_contigs, true, false, _nbCores);
		readParser.parse(ContigPartitionFunctor(*this));

		for(PartitionFile* partitionFile : _partitionFiles){
			delete partitionFile;
		}
		_partitionFiles.clear();

	}

	void collectContigStats(){
		auto fp = std::bind(&ContigPolisher::collectContigStats_read, this, std::placeholders::_1);
		ReadParser readParser(_inputFilename_contigs, true, false);
		readParser.parse(fp);
	}
	
	void collectContigStats_read(const Read& read){
		_contigStats.push_back({(u_int32_t)read._index, (u_int32_t)read._seq.size()});
	}

	u_int64_t computeContigMemory(u_int32_t contigIndex, size_t contigLength){
		
		u_int64_t windowByteSize = 0;
		u_int32_t contigCoverage = _contigCoverages[contigIndex];
		contigCoverage = max((u_int32_t)1, contigCoverage);

		//cout << contigCoverage << endl;
		if(_useQual){
			windowByteSize = contigCoverage  * ((_windowLength/4) + _windowLength);
		}
		else{
			windowByteSize = contigCoverage  * (_windowLength/4);
		}

		size_t nbWindows = ceil((double)contigLength / (double)_windowLength);
		return nbWindows * windowByteSize;
	}


	class ReadPartitionFunctor {

		public:

		ContigPolisher& _contigPolisher;


		ReadPartitionFunctor(ContigPolisher& contigPolisher) : _contigPolisher(contigPolisher){
		}

		ReadPartitionFunctor(const ReadPartitionFunctor& copy) : _contigPolisher(copy._contigPolisher){
		}

		~ReadPartitionFunctor(){
		}

		void operator () (const Read& read) {

			u_int32_t readIndex = read._index;
			//const string& readName = Utils::shortenHeader(read._header);

			if(read._index % 10000 == 0) Logger::get().debug() << "\t" << read._index;
			
			if(_contigPolisher._alignments.find(readIndex) == _contigPolisher._alignments.end()) return;


			//unordered_set<u_int32_t> writtenPartitions;

			//for(const Alignment& al : _contigPolisher._alignments[readIndex]){
			const Alignment& al = _contigPolisher._alignments[readIndex];
			u_int32_t contigIndex = al._contigIndex;
			//if(contigName != "ctg90297c") return;
			//u_int32_t contigIndex = al._contigIndex; //_contigPolisher._alignments[readIndex]._contigIndex;
			//_logFile << contigIndex << " " << (_contigPolisher._contigToPartition.find(contigIndex) != _contigPolisher._contigToPartition.end()) << endl;
			u_int32_t partition = _contigPolisher._contigToPartition[contigIndex];

			//if(writtenPartitions.find(partition) != writtenPartitions.end()) return;
			//writtenPartitions.insert(partition);

			//_logFile << partition << endl;
			PartitionFile* partitionFile = _contigPolisher._partitionFiles[partition];
			bool isFastq = read._qual.size() > 0 && _contigPolisher._useQual;

			string readSeq = read._seq;
			string qualSeq = read._qual;
			if(al._strand){
				Utils::toReverseComplement(readSeq);
				if(isFastq) std::reverse(qualSeq.begin(), qualSeq.end());
			}

			omp_set_lock(&partitionFile->_mutex);
			
			//if(partition == 0){
				//_logFile << "Read: " << readIndex << endl;
			//}


			_contigPolisher._partitionNbReads[partition] += 1;

			u_int32_t readSize = readSeq.size();

			if(isFastq){
				//string header = '@' + to_string(read._index) + '\n';
				string header = '@' + to_string(readIndex) + '\n';
				string seq = readSeq + '\n';
				gzwrite(partitionFile->_file, (const char*)&header[0], header.size());
				gzwrite(partitionFile->_file, (const char*)&seq[0], seq.size());
				static string strPlus = "+\n";
				gzwrite(partitionFile->_file, (const char*)&strPlus[0], strPlus.size());
				string qual = qualSeq + '\n';
				gzwrite(partitionFile->_file, (const char*)&qual[0], qual.size());
			}
			else{
				//string header = '>' + to_string(read._index) + '\n';
				string header = '>' + to_string(readIndex) + '\n';
				string seq = readSeq + '\n';
				gzwrite(partitionFile->_file, (const char*)&header[0], header.size());
				gzwrite(partitionFile->_file, (const char*)&seq[0], seq.size());
			}

			//for(char letter : readName){
			//	_contigPolisher._checksumWrittenReads += letter;
			//}
			//_logFile << readName << " " << _contigPolisher._checksumWrittenReads << endl;

			//gzwrite(partitionFile->_file, (const char*)&readIndex, sizeof(readIndex));
			//gzwrite(partitionFile->_file, (const char*)&readSize, sizeof(readSize));
			//gzwrite(partitionFile->_file, read._seq.c_str(), read._seq.size());
			//gzwrite(partitionFile->_file, read._qual.c_str(), read._qual.size());
			omp_unset_lock(&partitionFile->_mutex);
			//}
		}
	};

	class ContigPartitionFunctor {

		public:

		ContigPolisher& _parent;


		ContigPartitionFunctor(ContigPolisher& parent) : _parent(parent){
		}

		ContigPartitionFunctor(const ContigPartitionFunctor& copy) : _parent(copy._parent){
		}

		~ContigPartitionFunctor(){
		}

		void operator () (const Read& read) {

			u_int32_t contigIndex = read._index;

			if(read._index % 10000 == 0) Logger::get().debug() << "\t" << read._index;
			
			if(_parent._contigToPartition.find(contigIndex) == _parent._contigToPartition.end()) return;

			u_int32_t partition = _parent._contigToPartition[contigIndex];

			PartitionFile* partitionFile = _parent._partitionFiles[partition];

			omp_set_lock(&partitionFile->_mutex);
			
			string header = '>' + read._header + '\n';
			string seq = read._seq + '\n';
			gzwrite(partitionFile->_file, (const char*)&header[0], header.size());
			gzwrite(partitionFile->_file, (const char*)&seq[0], seq.size());
			

			omp_unset_lock(&partitionFile->_mutex);
		}
	};

	struct MinimizerPosition{
		MinimizerType _minimizer;
		int32_t _position;
		u_int8_t _isReversed;
	};

	u_int32_t _currentPartition;
	phmap::parallel_flat_hash_map<u_int32_t, phmap::parallel_flat_hash_map<u_int32_t, Alignment>> _alignments_readsVsContigs;
	//phmap::parallel_flat_hash_map<ReadType, bool> _alignments_isReversed;
	phmap::parallel_flat_hash_map<std::pair<u_int32_t, u_int32_t>, Alignment> _alignments_readsVsReads;
	vector<OverlapGraphPath> _contigOverlapPaths;
	phmap::parallel_flat_hash_map<u_int32_t, DnaBitset2*> _readIndex_to_sequence;
	phmap::parallel_flat_hash_map<u_int32_t, vector<MinimizerPosition>> _readIndex_to_minimizers;
	unordered_set<ReadType> _usedReadIndexes;

	
	bool _loadAllContigs;

	void loadContigs(const string& contigFilename, bool loadAllContigs){

		_contigSequences.clear();
		_contigHeaders.clear();
		_contigWindowSequences.clear();
		
		_loadAllContigs = loadAllContigs;
		auto fp = std::bind(&ContigPolisher::loadContigs_read, this, std::placeholders::_1);
		ReadParser readParser(contigFilename, true, false);
		readParser.parse(fp);
		
	}
	
	void loadContigs_read(const Read& read){

		
		u_int32_t contigIndex = read._index;
		//const string& contigName = Utils::shortenHeader(read._header);
		u_int32_t contigPartition = _contigToPartition[contigIndex];

		if(!_loadAllContigs && contigPartition != _currentPartition) return;
		

		//if(_currentPartition == 0) _logFile << "Loading contig in partition 0: " << contigIndex << endl;
		_contigSequences[contigIndex] = read._seq;

		size_t nbWindows = ceil(((double)read._seq.size()) / (double)_windowLength);

		//cout << "loadContigs_read: Nb windows + 1 si offsrt polishing" << endl;
		vector<vector<Window>> windows(nbWindows);

		_contigWindowSequences[contigIndex] = windows;
		_contigHeaders[contigIndex] = read._header;

		//cout << "Loaded contig: " << contigIndex << " " << read._header << " " << read._seq.size() << endl;
	}

	/*
	void mapReadsVsReads(const string& readFilename, const string& alignFilename){
		
    	string command = "minimap2 -v 0 -m 500 -t " + to_string(_nbCores) + " -x ava-ont -c " + readFilename + " " + readFilename;
		command += " | gzip -c - > " + alignFilename;

		Utils::executeCommand(command, _tmpDir);

		//executeCommand(command);
    	//os.system(command)
	}

	void loadReadsVsReadsAlignments(string alignFilename){
				
		cout << "Load alignment reads vs reads: " << alignFilename << endl;
		
		PafParser pafParser(alignFilename);
		auto fp = std::bind(&ContigPolisher::loadReadsVsReadsAlignments_read, this, std::placeholders::_1);
		pafParser.parse(fp);

	}

	void loadReadsVsReadsAlignments_read(const string& line){

		vector<string> _fields = Utils::split(line, '\t');
		//GfaParser::tokenize(line, _fields, '\t');

		//_logFile << line << endl;

		
		//const string& queryName = Utils::shortenHeader((*_fields)[0]);
		u_int32_t queryIndex = stoull((_fields)[0]);
		int64_t queryLength = stoull((_fields)[1]);
		int64_t queryStart = stoull((_fields)[2]);
		int64_t queryEnd = stoull((_fields)[3]);
		bool isReversed = (_fields)[4] == "-";
		//const string& targetName = Utils::shortenHeader((*_fields)[5]);
		u_int32_t targetIndex = stoull((_fields)[5]);
		int64_t targetLength = stoull((_fields)[6]);
		int64_t targetStart = stoull((_fields)[7]);
		int64_t targetEnd = stoull((_fields)[8]);

		int64_t nbMatches = stoull((_fields)[9]);
		int64_t alignLength = stoull((_fields)[10]);

		

		//double length = max((double)(contigEnd - contigStart), (double)(readEnd - readStart));
		//double error = 1 - min((double)(contigEnd - contigStart), (double)(readEnd - readStart)) / length;
		//_logFile << error << " " << errorThreshold << endl;
		
		//if(error > errorThreshold) return;

		float identity =  ((float)nbMatches) / alignLength;

		//u_int32_t queryIndex = _readName_to_readIndex[readName];
		//u_int32_t targetIndex = _readName_to_readIndex[contigName];


		int64_t tl3 = 0;
		int64_t tl5 = 0;
		int64_t ext3 = 0;
		int64_t ext5 = 0;

        if(isReversed){
            tl5 = targetLength - targetEnd;
            tl3 = targetStart;
		}
        else{
            tl5 = targetStart;
            tl3 = targetLength - targetEnd;
		}

        if(queryStart < tl5){
            ext5 = queryStart;
		}
        else{
            ext5 = tl5;
		}

        if (queryLength - queryEnd < tl3){
            ext3 = queryLength - queryEnd;
		}
        else{
            ext3 = tl3;
		}

        if (ext5 > maxHang || ext3 > maxHang || queryEnd - queryStart < (queryEnd - queryStart + ext5 + ext3) * intFrac) return;

        if (queryStart <= tl5 && queryLength - queryEnd <= tl3) return; //query contained
        else if (queryStart >= tl5 && queryLength - queryEnd >= tl3) return; //target contained

        if (queryEnd - queryStart + ext5 + ext3 < minOverlap) return; //Short overlap
        if (targetEnd - targetStart + ext5 + ext3 < minOverlap) return; //Short overlap

		std::pair<u_int32_t, u_int32_t> key;

        if (queryIndex < targetIndex){
			key.first = queryIndex;
			key.second = targetIndex;
		}
        else{
			key.first = targetIndex;
			key.second = queryIndex;
		}

		Alignment align = {targetIndex, queryIndex, isReversed, queryStart, queryEnd, targetStart, targetEnd, identity};

		if(_alignments_readsVsReads.find(key) != _alignments_readsVsReads.end()){
			if(align.alignLength() > _alignments_readsVsReads[key].alignLength()){
				_alignments_readsVsReads[key] = align;
			}
		}
		else{
			_alignments_readsVsReads[key] = align;
		}
       
	}
	*/


	//_validContigIndexes.insert(contigIndex);

	//u_int64_t totalByteSize = 0;


	//_logFile << "load contig: " << contigIndex << " " << seq.size() << endl;
	//_logFile << "Nb windows: " << nbWindows << endl;

	

	//_logFile << nbWindows << " " << nbWindows * _windowByteSize << endl;
	//return nbWindows * _windowByteSize;

	int32_t getReadLength(const ReadType& readIndex){
		return _readIndex_to_sequence[readIndex]->m_len;
	}

	void clearPass(){

		clearOverlap();
		_contigSequences.clear();
		//_validContigIndexes.clear();
		_contigWindowSequences.clear();
		_contigHeaders.clear();
		//_alignments_isReversed.clear();
		_usedReadIndexes.clear();
		//_alignments.clear();
	}

	void clearOverlap(){

		_readIndex_to_minimizers.clear();

		for(const auto& it: _readIndex_to_sequence){
			delete it.second;
		}
		_readIndex_to_sequence.clear();

		_contigOverlapPaths.clear();
		_alignments_readsVsContigs.clear();
		_alignments_readsVsReads.clear();
		//_alignments_readsVsReads_cached.clear();
	}


	void indexContigName(const string& contigFilename){
		 
		//cout << "Indexing contig names: " << contigFilename << endl;
		_contigName_to_contigIndex.clear();

		auto fp = std::bind(&ContigPolisher::indexContigName_read, this, std::placeholders::_1);
		ReadParser readParser(contigFilename, true, false);
		readParser.parse(fp);
	}

	void indexContigName_read(const Read& read){

		//cout << "Index contig: " << Utils::shortenHeader(read._header) << " " << read._index << endl;
		_contigName_to_contigIndex[Utils::shortenHeader(read._header)] = read._index;
	}

	void indexReadName(){

		//cout << "Indexing read names" << endl;
		//cout << _inputFilename_reads << endl;
		auto fp = std::bind(&ContigPolisher::indexReadName_read, this, std::placeholders::_1);
		ReadParser readParser(_inputFilename_reads, false, false);
		readParser.parse(fp);

	}

	void indexReadName_read(const Read& read){
		//if(read._index % 100000 == 0) cout << read._index << endl;
		_readName_to_readIndex[Utils::shortenHeader(read._header)] = read._index;
	}


	struct AlignmentPartitionning{
		u_int32_t _contigIndex;
		u_int32_t _length;
	};

	//unordered_map<string, u_int32_t> _contigName_to_contigIndex;
	//unordered_map<string, u_int64_t> _readName_to_readIndex;
	//unordered_map<u_int64_t, AlignmentPartitionning> _readToContigIndex;
	//unordered_map<u_int64_t, vector<Alignment>> _alignments;
	phmap::parallel_flat_hash_map<u_int32_t, Alignment> _alignments;
	phmap::parallel_flat_hash_map<u_int32_t, string> _contigSequences;
	phmap::parallel_flat_hash_map<u_int32_t, vector<vector<Window>>> _contigWindowSequences;
	phmap::parallel_flat_hash_map<u_int32_t, string> _contigHeaders;
	//unordered_map<ContigRead, u_int32_t, ContigRead_hash> _alignmentCounts;
	//u_int64_t _correctedContigIndex;

	phmap::parallel_flat_hash_map<string, u_int32_t> _contigName_to_contigIndex;
	phmap::parallel_flat_hash_map<string, u_int32_t> _readName_to_readIndex;

	/*
	void indexContigName(){
		
		_logFile << "Indexing contig names" << endl;

		auto fp = std::bind(&ContigPolisher::indexContigName_read, this, std::placeholders::_1);
		ReadParser readParser(_inputFilename_contigs, true, false);
		readParser.parse(fp);
	}
	
	void indexContigName_read(const Read& read){

		string minimap_name;

		auto find = read._header.find(' ');
		if(find == std::string::npos){
			minimap_name = read._header;
		}
		else{
		 	minimap_name = read._header.substr(0, find);
		}

		//_logFile << minimap_name << endl;
		_contigName_to_contigIndex[minimap_name] = read._index;
	}

	void indexReadName(){

		_logFile << "Indexing read names" << endl;
		auto fp = std::bind(&ContigPolisher::indexReadName_read, this, std::placeholders::_1);
		ReadParser readParser(_inputFilename_reads, false, false);
		readParser.parse(fp);

	}
	
	void indexReadName_read(const Read& read){
		//_logFile << read._index << endl;

		string minimap_name;

		auto find = read._header.find(' ');
		if(find == std::string::npos){
			minimap_name = read._header;
		}
		else{
		 	minimap_name = read._header.substr(0, find);
		}

		//string minimap_name = read._header.substr(0, read._header.find(' '));
		//_logFile << minimap_name << endl;
		_readName_to_readIndex[minimap_name] = read._index;
	}
	*/


	vector<string>* _fields;
	vector<string>* _fields_optional;
	bool _useIndexedName;
	bool _indexPerContig;
	bool _useErrorFilter;

	void loadAlignments(const string& alignFilename, bool useIndexedName, bool indexPerContig, bool useErrorFilter){

		Logger::get().debug() << "\tLoading alignments";

		_alignments.clear();
		_alignments_readsVsContigs.clear();

		_useIndexedName = useIndexedName;
		_indexPerContig = indexPerContig;
		_useErrorFilter = useErrorFilter;
		//_indexPartitionOnly = indexPartitionOnly;

		PafParser pafParser(alignFilename);
		auto fp = std::bind(&ContigPolisher::loadAlignments_read, this, std::placeholders::_1);
		pafParser.parse(fp);

	}

	void loadAlignments_read(const string& line){

		//cout << line << endl;
		double errorThreshold = 0.3;
		//_logFile << lineInput << endl;
		//getchar();


		vector<string> _fields = Utils::split(line, '\t');
		//GfaParser::tokenize(line, _fields, '\t');

		//_logFile << line << endl;

		const string& readName = Utils::shortenHeader((_fields)[0]);
		const string& contigName = Utils::shortenHeader((_fields)[5]);
		

		u_int64_t readLength = stoull((_fields)[1]);
		u_int32_t readStart = stoull((_fields)[2]);
		u_int32_t readEnd = stoull((_fields)[3]);
		u_int32_t contigLength = stoull((_fields)[6]);
		u_int32_t contigStart = stoull((_fields)[7]);
		u_int32_t contigEnd = stoull((_fields)[8]);

		u_int64_t nbMatches = stoull((_fields)[9]);
		u_int64_t alignLength = stoull((_fields)[10]);

		bool strand = (_fields)[4] == "-";
		//float score = (double) nbMatches / (double) queryLength;
		//float score2 = (double) nbMatches / (double) alignLength;

		//u_int32_t contigIndex = mappingIndex._contigName_to_contigIndex[contigName];
		//u_int64_t readIndex = mappingIndex._readName_to_readIndex[readName];
		//u_int64_t length = std::max(readEnd - readStart, contigEnd - contigStart);


		//_logFile << contigIndex << " " << readIndex << endl;

		double length = max((double)(contigEnd - contigStart), (double)(readEnd - readStart));
		double error = 1 - min((double)(contigEnd - contigStart), (double)(readEnd - readStart)) / length;
		//_logFile << error << " " << errorThreshold << endl;
		
		if(_useErrorFilter && error > errorThreshold) return;

		//float identity =  ((float)nbMatches) / alignLength;


		u_int32_t readIndex;
		u_int32_t contigIndex = _contigName_to_contigIndex[contigName];

		if(_useIndexedName){
			readIndex = stoull(readName);
		}
		else{
			readIndex = _readName_to_readIndex[readName];
		}
		
		//cout << readIndex<< " " << contigIndex << " " << contigName << endl;
		
		//if(indexPartitionOnly && _contigSequences.find(contigName) == _contigSequences.end()) continue;


		//u_int32_t length = std::max((u_int64_t)(readEnd - readStart), (u_int64_t)(contigEnd - contigStart));

		AlignmentBounds bounds;
		bounds._queryStart = readStart;
		bounds._queryEnd = readEnd;
		bounds._referenceStart = contigStart;
		bounds._referenceEnd = contigEnd;
		bounds._queryLength = readLength;
		bounds._referenceLength = contigLength;
		bounds._isReversed = strand;

		//_contigHits[contigIndex].push_back({contigStart, contigEnd});

		u_int32_t readSizeMappable = getMappableLength(bounds); //A read size that do not consider alignment outside contig bounds (start and end of the contigs)
		
		float identity = ((long double) nbMatches) / ((long double) readSizeMappable);

		Alignment alignment = {contigIndex, readIndex, strand, readStart, readEnd, contigStart, contigEnd, identity, readLength, contigLength}; //, score

		if(_indexPerContig){

			if(_alignments_readsVsContigs.find(contigIndex) == _alignments_readsVsContigs.end()){
				_alignments_readsVsContigs[contigIndex][readIndex] = alignment;
				//_alignments_isReversed[readIndex] = alignment._strand;
			}
			else{
				auto& map = _alignments_readsVsContigs[contigIndex];

				if(map.find(readIndex) == map.end()){
					map[readIndex] = alignment;
					//_alignments_isReversed[readIndex] = alignment._strand;
				}
				else{

					if(alignment.alignLength() > map[readIndex].alignLength()){
						map[readIndex] = alignment;
						//_alignments_isReversed[readIndex] = alignment._strand;
					}
				}
			}

			//cout << "Load align: " << readIndex << " " << contigIndex << endl;
		}
		else{
			

			
			if(_alignments.find(readIndex) == _alignments.end()){
				_alignments[readIndex] = alignment;
			}
			else{

				if(alignment.score() >= _alignments[readIndex].score()){
					_alignments[readIndex] = alignment;
				}
				
			}


		}



	}


	
	class ReadLoaderFunctor {

		public:

		ContigPolisher& _parent;
		//MinimizerParser* _minimizerParser;


		ReadLoaderFunctor(ContigPolisher& parent) : _parent(parent){
			//_minimizerParser = new MinimizerParser(ContigPolisher::_minimizerSize, ContigPolisher::_minimizerDensity);
		}

		ReadLoaderFunctor(const ReadLoaderFunctor& copy) : _parent(copy._parent){
			//_minimizerParser = new MinimizerParser(ContigPolisher::_minimizerSize, ContigPolisher::_minimizerDensity);
		}

		~ReadLoaderFunctor(){
			//delete _minimizerParser;
		}

		void operator () (const Read& read) {

			u_int32_t readIndex = stoull(read._header);

			//if(_parent._readIndex_to_sequence.find(readIndex) ==_parent._readIndex_to_sequence.end()) return;

			//if(_parent._alignments.find(readIndex) == _parent._alignments.end()) return;

			//const Alignment& readAlignment = _parent._alignments[readIndex];
			/*
			string readSequence = read._seq;
			//if(_parent._alignments_isReversed[readIndex]){
			if(_parent._alignments[readIndex]._strand){
				Utils::toReverseComplement(readSequence);
			}

			vector<MinimizerType> minimizers;
			vector<u_int32_t> minimizerPos;
			vector<u_int8_t> minimizerDirections;
			_minimizerParser->parse(readSequence, minimizers, minimizerPos, minimizerDirections);


			vector<MinimizerPosition> minimizerPositions;

			for(size_t i=0; i<minimizers.size(); i++){

				MinimizerType minimizer = minimizers[i];
				u_int32_t position = minimizerPos[i];
				bool isReversed = minimizerDirections[i];

				MinimizerPosition minmizerPosition = {minimizer, position, isReversed};

				minimizerPositions.push_back(minmizerPosition);
			}

			std::sort(minimizerPositions.begin(), minimizerPositions.end(), [](const MinimizerPosition& a, const MinimizerPosition& b){

				if(a._minimizer == b._minimizer){
					return a._position < b._position;
				}

				return a._minimizer < b._minimizer;
			});
			*/
			#pragma omp critical(ReadLoaderFunctor)
			{
				_parent._readIndex_to_sequence[readIndex] = new DnaBitset2(read._seq);
				//_parent._readIndex_to_minimizers[readIndex] = minimizerPositions;
			}

		}
	};

	/*
	void computeReadsToContigsAlignments(u_int32_t partition){
		
		

		//auto fp = std::bind(&ContigPolisher::collectWindowCopies_read, this, std::placeholders::_1);
		//ReadParser readParser(_inputFilename_reads, false, false);
		//readParser.parse(fp);

		//cout << "single core here" << endl;

		const string& partitionFilename = _readPartitionDir + "/part_" + to_string(partition) + ".gz";
		ReadParserParallel readParser(partitionFilename, true, false, _nbCores);
		readParser.parse(ReadsToContigsAlignFunctor(*this));
	}
	*/

	class ReadsToContigsAlignFunctor {

		public:

		ContigPolisher& _contigPolisher;
		phmap::parallel_flat_hash_map<u_int32_t, Alignment>& _alignments;
		phmap::parallel_flat_hash_map<u_int32_t, string>& _contigSequences;
		phmap::parallel_flat_hash_map<u_int32_t, phmap::parallel_flat_hash_map<u_int32_t, Alignment>>& _alignments_readsVsContigs;

		ReadsToContigsAlignFunctor(ContigPolisher& contigPolisher) : _contigPolisher(contigPolisher), _alignments(contigPolisher._alignments), _contigSequences(contigPolisher._contigSequences), _alignments_readsVsContigs(_contigPolisher._alignments_readsVsContigs){
		}

		ReadsToContigsAlignFunctor(const ReadsToContigsAlignFunctor& copy) : _contigPolisher(copy._contigPolisher), _alignments(copy._alignments), _contigSequences(copy._contigSequences), _alignments_readsVsContigs(copy._alignments_readsVsContigs){
			
		}

		~ReadsToContigsAlignFunctor(){
		}

		void operator () (const Read& read) {

			//const string& readName = Utils::shortenHeader(read._header);
			u_int32_t readIndex = stoull(read._header);

			//if(_contigPolisher._currentPartition == 0) _logFile << readIndex << " " << (_alignments.find(readIndex) != _alignments.end()) << endl;
			
			//if(readIndex % 100000 == 0) _logFile << "\t" << readIndex << endl;

			if(_alignments.find(readIndex) == _alignments.end()) return;

			//cout << readIndex << endl;

			//const vector<Alignment>& als = _alignments[readIndex];
			Alignment& al = _alignments[readIndex];
			//for(const Alignment& al : _alignments[readIndex]){
			u_int32_t contigIndex = al._contigIndex;

			if(_contigSequences.find(contigIndex) == _contigSequences.end()) return;

			//_logFile << read._seq.size() << " " << read._qual.size() << " " << _contigSequences[contigIndex].size() << " " << al._readStart << " " << al._readEnd << " " << al._contigStart << " " << al._contigEnd << endl;
			string readSeq = read._seq;
			string qualSeq = read._qual;

			if(al._strand){
				Utils::toReverseComplement(readSeq);
				std::reverse(qualSeq.begin(), qualSeq.end());
			}

			string readSequence = readSeq.substr(al._readStart, al._readEnd-al._readStart);
			string contigSequence = _contigSequences[contigIndex].substr(al._contigStart, al._contigEnd-al._contigStart);




			if(al._strand){
				Utils::toReverseComplement(readSequence);
				Utils::toReverseComplement(readSeq);
				std::reverse(qualSeq.begin(), qualSeq.end());
			}
			

			AlignmentBounds bounds;
			bounds._queryStart = al._readStart;
			bounds._queryEnd = al._readEnd;
			bounds._referenceStart = al._contigStart;
			bounds._referenceEnd = al._contigEnd;
			bounds._queryLength = readSeq.size();
			bounds._referenceLength = _contigSequences[contigIndex].size();
			bounds._isReversed = al._strand;


			u_int32_t readSizeMappable = _contigPolisher.getMappableLength(bounds); //A read size that do not consider alignment outside contig bounds (start and end of the contigs)
			

			int64_t overhang = readSeq.size() - readSizeMappable;
			//cout << overhang << endl;
			//if(overhang > ContigPolisher::maxHang) return;





			static EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0);


			EdlibAlignResult result = edlibAlign(readSequence.c_str(), readSequence.size(), contigSequence.c_str(), contigSequence.size(), config);


			char* cigar;

			if (result.status == EDLIB_STATUS_OK) {
				cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
			} else {
				Logger::get().error() << "Invalid edlib results";
				exit(1);
			}

			edlibFreeAlignResult(result);
			
			u_int64_t nbMatches = getCigarNbMatches(al, cigar, readSequence, contigSequence);
			free(cigar);





			float identity = ((long double) nbMatches) / ((long double) readSizeMappable);
			Alignment alignment = {al._contigIndex, readIndex, al._strand, al._readStart, al._readEnd, al._contigStart, al._contigEnd, identity, al._readLength, al._contigLength};
			

			//cout << identity << endl;

			//if(!_contigPolisher.isValidAlignment(bounds)) return;

			//	al.printAll();


			//cout << contigIndex << " " << readIndex << " " << _contigSequences[contigIndex].size() << endl;

			#pragma omp critical(ReadsToContigsAlignFunctor)
			{

				
				if(_alignments_readsVsContigs.find(contigIndex) == _alignments_readsVsContigs.end()){
					_alignments_readsVsContigs[contigIndex][readIndex] = alignment;
				}
				else{
					auto& map = _alignments_readsVsContigs[contigIndex];

					if(map.find(readIndex) == map.end()){
						map[readIndex] = alignment;
					}
					else{

						if(alignment.alignLength() > map[readIndex].alignLength()){
							map[readIndex] = alignment;
						}
					}
				}
			}

			//Alignment align = {contigIndex, strand, readStart, readEnd, contigStart, contigEnd, identity}; //, score



			//getchar();

			//}

		}

		u_int64_t getCigarNbMatches(const Alignment& al, char* cigar_, const string& readSequence, const string& contigSequence){
			
			u_int64_t nbMatches = 0;

			//u_int64_t t_begin_ = al._contigStart;
			//u_int64_t t_end_ = al._contigEnd;
			//u_int64_t q_begin_ = al._readStart;
			//u_int64_t q_end_ = al._readEnd;
			bool strand_ = al._strand;
			//u_int64_t q_length_ = readSequence.size();

			int32_t q_ptr = -1; //(strand_ ? (q_length_ - q_end_) : q_begin_) - 1;
			int32_t t_ptr = -1; //t_begin_ - 1;
			

			for (uint32_t i = 0, j = 0; i < strlen(cigar_); ++i) {


				if (cigar_[i] == 'M' || cigar_[i] == '=' || cigar_[i] == 'X') {
					
					//cout << "M: " << atoi(&cigar_[j]) << endl;

					uint32_t k = 0, num_bases = atoi(&cigar_[j]);
					j = i + 1;
					while (k < num_bases) {
						++q_ptr;
						++t_ptr;

						//cout << q_ptr << " " << readSequence.size() << "    " << t_ptr << " " << contigSequence.size() << endl;
						if(readSequence[q_ptr] == contigSequence[t_ptr]){
							nbMatches += 1;
						}

						++k;
					}
				} else if (cigar_[i] == 'I') {
					//cout << "I: " << atoi(&cigar_[j]) << endl;
					q_ptr += atoi(&cigar_[j]);
					j = i + 1;
				} else if (cigar_[i] == 'D' || cigar_[i] == 'N') {
					//cout << "D: " << atoi(&cigar_[j]) << endl;
					//uint32_t k = 0, num_bases = atoi(&cigar_[j]);
					t_ptr += atoi(&cigar_[j]);
					j = i + 1;
					//while (k < num_bases) {
					//	++t_ptr;
					//	++k;
					//}
				} else if (cigar_[i] == 'S' || cigar_[i] == 'H' || cigar_[i] == 'P') {
					//cout << "SSSS: " << atoi(&cigar_[j]) << endl;
					j = i + 1;
				}
			}

			//cout << contigSequence << endl;
			//cout << readSequence << endl;
			//cout << nbMatches << endl;
			//getchar();

			return nbMatches;
		}



	};


	void polishPartition(u_int32_t partition, const string& contigFilename, gzFile& outputContigFile, const int& pass){
		_currentPartition = partition;
		//if(_contigSequences.size() == 0) return;
		
		Logger::get().debug() << "";
		Logger::get().debug() << "\tPolishing partition: " << _currentPartition << "/" << _nbPartitions;
		auto start = high_resolution_clock::now();

		if(_partitionNbReads[partition] == 0) return;

		clearPass();
		//string contigFilename = _readPartitionDir + "/part_" + to_string(partition) + "_contigs.gz";
		//string readFilename = _readPartitionDir + "/part_" + to_string(partition) + ".gz";


		string readFilename = _readPartitionDir + "/" + to_string(partition) + "_reads.gz";
		//string contigFilename = _readPartitionDir + "/" + to_string(partition) + "_contigsRepeats.gz";
		
		//string outputContigFilename = _readPartitionDir + "/" + to_string(partition) + "_contigs.gz";

		
		Logger::get().debug() << "\tMap reads to curated contigs";

		const string& alignFilename = _readPartitionDir + "/align.paf.gz";
		string command = "minimap2 -v 0 -m 500 -t " + to_string(_nbCores) + " -x " + _minimap2Preset_map + " " + contigFilename + " " + readFilename;
		Utils::executeMinimap2(command, alignFilename);
		//command += " | gzip -c - > " + alignFilename;
		//cout << command << endl;
		//Utils::executeCommand(command, _tmpDir);

		indexContigName(contigFilename);
		loadContigs(contigFilename, true);
		loadAlignments(alignFilename, true, false, true);
		computeContigCoverages(contigFilename);
		


		
		//cout << "Load contigs" << endl;
		//loadContigs(_inputFilename_contigs, false);

		Logger::get().debug() << "\tCollecting window sequences";

		ReadParserParallel readParser(readFilename, true, false, _nbCores);
		readParser.parse(CollectWindowSequencesFunctor(*this));

		//cout << "miuom" << endl;
		//getchar();
		
		Logger::get().debug() << "\tPerform correction";
		performCorrection(outputContigFile, pass);


		Logger::get().debug() << "\tPartition " << partition << " done " << " (" << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s)";



	}

	class CollectWindowSequencesFunctor {

		public:

		ContigPolisher& _contigPolisher;
		//unordered_map<string, u_int32_t>& _contigName_to_contigIndex;
		//unordered_map<string, u_int64_t>& _readName_to_readIndex;
		//unordered_map<u_int64_t, vector<Alignment>>& _alignments;
		phmap::parallel_flat_hash_map<u_int32_t, Alignment>& _alignments;
		phmap::parallel_flat_hash_map<u_int32_t, string>& _contigSequences;
		phmap::parallel_flat_hash_map<u_int32_t, vector<vector<Window>>>& _contigWindowSequences;
		size_t _windowLength;


		CollectWindowSequencesFunctor(ContigPolisher& contigPolisher) : _contigPolisher(contigPolisher), _alignments(contigPolisher._alignments), _contigSequences(contigPolisher._contigSequences), _contigWindowSequences(contigPolisher._contigWindowSequences), _windowLength(contigPolisher._windowLength){
		}

		CollectWindowSequencesFunctor(const CollectWindowSequencesFunctor& copy) : _contigPolisher(copy._contigPolisher), _alignments(copy._alignments), _contigSequences(copy._contigSequences), _contigWindowSequences(copy._contigWindowSequences), _windowLength(copy._windowLength){
			
		}

		~CollectWindowSequencesFunctor(){
		}

		void operator () (const Read& read) {

			//const string& readName = Utils::shortenHeader(read._header);
			u_int32_t readIndex = stoull(read._header);

			//if(_contigPolisher._currentPartition == 0) _logFile << readIndex << " " << (_alignments.find(readIndex) != _alignments.end()) << endl;
			
			//if(readIndex % 100000 == 0) _logFile << "\t" << readIndex << endl;

			if(_alignments.find(readIndex) == _alignments.end()) return;

			//const vector<Alignment>& als = _alignments[readIndex];
			const Alignment& al = _alignments[readIndex];
			//for(const Alignment& al : _alignments[readIndex]){
			u_int32_t contigIndex = al._contigIndex;

			if(_contigSequences.find(contigIndex) == _contigSequences.end()) return;

			//_logFile << read._seq.size() << " " << read._qual.size() << " " << _contigSequences[contigIndex].size() << " " << al._readStart << " " << al._readEnd << " " << al._contigStart << " " << al._contigEnd << endl;
			string readSeq = read._seq;
			string qualSeq = read._qual;
			string readSequence = readSeq.substr(al._readStart, al._readEnd-al._readStart);
			string contigSequence = _contigSequences[contigIndex].substr(al._contigStart, al._contigEnd-al._contigStart);



			if(al._strand){
				Utils::toReverseComplement(readSequence);
				Utils::toReverseComplement(readSeq);
				std::reverse(qualSeq.begin(), qualSeq.end());
			}
			

			AlignmentBounds bounds;
			bounds._queryStart = al._readStart;
			bounds._queryEnd = al._readEnd;
			bounds._referenceStart = al._contigStart;
			bounds._referenceEnd = al._contigEnd;
			bounds._queryLength = readSeq.size();
			bounds._referenceLength = _contigSequences[contigIndex].size();
			bounds._isReversed = al._strand;

			//if(_contigPolisher.getMaxhang(bounds) > ContigPolisher::maxHang*2) return;

			//u_int32_t readSizeMappable = _contigPolisher.getMappableLength(bounds); //A read size that do not consider alignment outside contig bounds (start and end of the contigs)
			

			//int64_t overhang = readSeq.size() - readSizeMappable;
			//cout << overhang << endl;
			//if(overhang > 10) return;

			//cout << overhang << endl;


			//_logFile << readSequence << endl;
			//_logFile << contigSequence << endl;

			//_logFile << contigSequence.size() << " "<< readSequence.size() << endl;
			static EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0);


			EdlibAlignResult result = edlibAlign(readSequence.c_str(), readSequence.size(), contigSequence.c_str(), contigSequence.size(), config);


			char* cigar;

			if (result.status == EDLIB_STATUS_OK) {
				cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
			} else {
				Logger::get().error() << "Invalid edlib results";
				exit(1);
			}

			//_logFile << cigar << endl;

			edlibFreeAlignResult(result);
			
			//if(al._contigEnd < 20000){
				//cout << "-----" << endl;
				//cout << al._contigStart << " " << al._contigEnd << "     " << al._readStart << " " << al._readEnd << " " << al._readLength << endl;
				//cout << _contigPolisher.getMaxhang(bounds) << endl;
				//u_int64_t nbMatches = getCigarNbMatches(al, cigar, readSequence, contigSequence);
				//getchar();

			//}

			find_breaking_points_from_cigar(_windowLength, al, readSeq.size(), cigar, readSeq, qualSeq);


			free(cigar);

			//getchar();

			//}

		}

		/*
		u_int64_t getCigarNbMatches(const Alignment& al, char* cigar_, const string& readSequence, const string& contigSequence){
			
			u_int64_t nbMatches = 0;

			//u_int64_t t_begin_ = al._contigStart;
			//u_int64_t t_end_ = al._contigEnd;
			//u_int64_t q_begin_ = al._readStart;
			//u_int64_t q_end_ = al._readEnd;
			bool strand_ = al._strand;
			//u_int64_t q_length_ = readSequence.size();

			int32_t q_ptr = -1; //(strand_ ? (q_length_ - q_end_) : q_begin_) - 1;
			int32_t t_ptr = -1; //t_begin_ - 1;
			

			for (uint32_t i = 0, j = 0; i < strlen(cigar_); ++i) {


				if (cigar_[i] == 'M' || cigar_[i] == '=' || cigar_[i] == 'X') {
					
					//cout << "M: " << atoi(&cigar_[j]) << endl;

					uint32_t k = 0, num_bases = atoi(&cigar_[j]);
					j = i + 1;
					while (k < num_bases) {
						++q_ptr;
						++t_ptr;

						//cout << q_ptr << " " << readSequence.size() << "    " << t_ptr << " " << contigSequence.size() << endl;
						if(readSequence[q_ptr] == contigSequence[t_ptr]){
							nbMatches += 1;
						}

						++k;
					}
				} else if (cigar_[i] == 'I') {
					cout << "I: " << atoi(&cigar_[j]) << endl;
					q_ptr += atoi(&cigar_[j]);
					j = i + 1;
				} else if (cigar_[i] == 'D' || cigar_[i] == 'N') {
					cout << "D: " << atoi(&cigar_[j]) << endl;
					//uint32_t k = 0, num_bases = atoi(&cigar_[j]);
					t_ptr += atoi(&cigar_[j]);
					j = i + 1;
					//while (k < num_bases) {
					//	++t_ptr;
					//	++k;
					//}
				} else if (cigar_[i] == 'S' || cigar_[i] == 'H' || cigar_[i] == 'P') {
					cout << "SSSS: " << atoi(&cigar_[j]) << endl;
					j = i + 1;
				}
			}

			//cout << contigSequence << endl;
			//cout << readSequence << endl;
			//cout << nbMatches << endl;
			//getchar();

			return nbMatches;
		}
		*/

		void find_breaking_points_from_cigar(uint32_t window_length, const Alignment& al, u_int64_t readLength, char* cigar_, const string& readSequence, const string& qualSequence)
		{
			
			vector<std::pair<uint32_t, uint32_t>> breaking_points_;
			u_int64_t t_begin_ = al._contigStart;
			u_int64_t t_end_ = al._contigEnd;
			u_int64_t q_begin_ = al._readStart;
			u_int64_t q_end_ = al._readEnd;
			bool strand_ = al._strand;
			u_int64_t q_length_ = readLength;

			// find breaking points from cigar
			std::vector<int32_t> window_ends;
			for (uint32_t i = 0; i < t_end_; i += window_length) {
				if (i > t_begin_) {
					window_ends.emplace_back(i - 1);
				}
			}
			window_ends.emplace_back(t_end_ - 1);

			uint32_t w = 0;
			bool found_first_match = false;
			std::pair<uint32_t, uint32_t> first_match = {0, 0}, last_match = {0, 0};

			int32_t q_ptr = (strand_ ? (q_length_ - q_end_) : q_begin_) - 1;
			int32_t t_ptr = t_begin_ - 1;

			for (uint32_t i = 0, j = 0; i < strlen(cigar_); ++i) {
				if (cigar_[i] == 'M' || cigar_[i] == '=' || cigar_[i] == 'X') {
					uint32_t k = 0, num_bases = atoi(&cigar_[j]);
					j = i + 1;
					while (k < num_bases) {
						++q_ptr;
						++t_ptr;

						if (!found_first_match) {
							found_first_match = true;
							first_match.first = t_ptr;
							first_match.second = q_ptr;
						}
						last_match.first = t_ptr + 1;
						last_match.second = q_ptr + 1;
						if (t_ptr == window_ends[w]) {
							if (found_first_match) {
								breaking_points_.emplace_back(first_match);
								breaking_points_.emplace_back(last_match);
							}
							found_first_match = false;
							++w;
						}


						++k;
					}
				} else if (cigar_[i] == 'I') {
					q_ptr += atoi(&cigar_[j]);
					j = i + 1;
				} else if (cigar_[i] == 'D' || cigar_[i] == 'N') {
					uint32_t k = 0, num_bases = atoi(&cigar_[j]);
					j = i + 1;
					while (k < num_bases) {
						++t_ptr;
						if (t_ptr == window_ends[w]) {
							if (found_first_match) {
								breaking_points_.emplace_back(first_match);
								breaking_points_.emplace_back(last_match);
							}
							found_first_match = false;
							++w;
						}
						++k;
					}
				} else if (cigar_[i] == 'S' || cigar_[i] == 'H' || cigar_[i] == 'P') {
					j = i + 1;
				}
			}

			//if(breaking_points_.size() > 0) breaking_points_.emplace_back(last_match);
				
			for (uint32_t j = 0; j < breaking_points_.size(); j += 2) {

				if(breaking_points_[j].second >= readSequence.size()) return;
				if(breaking_points_[j + 1].second >= readSequence.size()) return;

				if (breaking_points_[j + 1].second - breaking_points_[j].second < 0.02 * _windowLength) {
					continue;
				}

				
				if (qualSequence.size() > 0) {

					//const auto& quality = overlaps[i]->strand() ? sequence->reverse_quality() : sequence->quality();
					double average_quality = 0;
					for (uint32_t k = breaking_points_[j].second; k < breaking_points_[j + 1].second; ++k) {
						average_quality += static_cast<uint32_t>(qualSequence[k]) - 33;
					}
					average_quality /= breaking_points_[j + 1].second - breaking_points_[j].second;

					if (average_quality < _contigPolisher._qualityThreshold) {
						continue;
					}
				}
				

				//uint64_t window_id = id_to_first_window_id[overlaps[i]->t_id()] +
				uint64_t window_id = breaking_points_[j].first / _windowLength;
				uint32_t window_start = (breaking_points_[j].first / _windowLength) * _windowLength;

				//const char* data = overlaps[i]->strand() ?
				//	&(sequence->reverse_complement()[breaking_points[j].second]) :
				//	&(sequence->data()[breaking_points[j].second]);
				const char* data = &readSequence[breaking_points_[j].second];
				uint32_t data_length = breaking_points_[j + 1].second - breaking_points_[j].second;

				string sequence = string(data, data_length);

				/*
				if(al._contigIndex == 1 && window_id > 1 && window_id < 1000){
					if(window_id == 192 || window_id == 193){
						_logFile << readSequence << endl;
						_logFile << sequence << endl;
						if(breaking_points_[j].second > 1) _logFile << readSequence[breaking_points_[j].second-1] << " " << readSequence[breaking_points_[j].second] << endl;
						if(readSequence.size() > readSequence[breaking_points_[j].second+data_length-1]) _logFile << readSequence[breaking_points_[j].second+data_length-1] << " " << readSequence[breaking_points_[j].second+data_length] << endl;
						_logFile << window_id << endl;
						//getchar();
					}
					

				}
				
				long pos = breaking_points_[j].second+data_length;
				while(true){
					if(pos >= readSequence.size()) break;
					char prevChar = readSequence[pos-1];
					char c = readSequence[pos];
					if(prevChar != c) break;

					//sequence += c;
					breaking_points_[j + 1].second += 1;
					//data_length += 1;

					pos += 1;
				}

				pos = ((long)breaking_points_[j].second)-1;
				while(true){
					if(pos < 0) break;
					char prevChar = readSequence[pos+1];
					char c = readSequence[pos];
					if(prevChar != c) break;

					//sequence += c;
					breaking_points_[j].second -= 1;
					//data_length += 1;

					pos -= 1;
				}
				*/

				data = &readSequence[breaking_points_[j].second];
				data_length = breaking_points_[j + 1].second - breaking_points_[j].second;
				sequence = string(data, data_length);

				/*
				if(al._contigIndex == 1 && window_id > 1 && window_id < 1000){
					if(window_id == 192 || window_id == 193){
						_logFile << readSequence << endl;
						_logFile << sequence << endl;
						if(breaking_points_[j].second > 1){
							_logFile << readSequence[breaking_points_[j].second-1] << " " << readSequence[breaking_points_[j].second] << endl;
							if(readSequence[breaking_points_[j].second-1] == readSequence[breaking_points_[j].second]) getchar();
						}
						if(readSequence.size() > readSequence[breaking_points_[j].second+data_length-1]) _logFile << readSequence[breaking_points_[j].second+data_length-1] << " " << readSequence[breaking_points_[j].second+data_length] << endl;
						_logFile << window_id << endl;
						//getchar();
					}
					

				}
				*/

                u_int32_t posStart = breaking_points_[j].first - window_start;
                u_int32_t posEnd =  breaking_points_[j + 1].first - window_start - 1;

				string quality = "";
				if(_contigPolisher._useQual && qualSequence.size() > 0){
					quality = string(&qualSequence[breaking_points_[j].second], data_length);
				}

				//if(al._contigIndex == 0 && window_id == 161){
					
					//_logFile << al._contigIndex << " " << al._readIndex << endl;
					//_logFile << (breaking_points_[j + 1].second - breaking_points_[j].second) << " " << (0.02 * _windowLength) << endl;
					//_logFile << sequence << endl;
					//_logFile << quality << endl;
					//getchar();
				//}

				//if(sequence.size() < 490) continue;
				indexWindow(al, window_id, posStart, posEnd, sequence, quality);

				//if(window_id == 24){
				//	cout << sequence << endl;
				//	getchar();
				//}

				//_logFile << window_id << " " << posStart << " " << posEnd << endl;
				//_logFile << sequence << endl;
				//getchar();
				/*
				const char* quality = overlaps[i]->strand() ?
					(sequence->reverse_quality().empty() ?
						nullptr : &(sequence->reverse_quality()[breaking_points[j].second]))
					:
					(sequence->quality().empty() ?
						nullptr : &(sequence->quality()[breaking_points[j].second]));
				uint32_t quality_length = quality == nullptr ? 0 : data_length;
				*/
			}

			/*
			for(size_t i=0; i<breaking_points_.size()-1; i++){
				u_int64_t contigWindowStart = breaking_points_[i].first;
				u_int64_t contigWindowEnd = breaking_points_[i+1].first;
				u_int64_t readWindowStart = breaking_points_[i].second;
				u_int64_t readWindowEnd = breaking_points_[i+1].second;

				//_logFile << contigWindowStart << " " << contigWindowEnd  << "      " << readWindowStart << " " << readWindowEnd << endl;

				//if(readWindowEnd-readWindowStart < _windowLength) continue; //window sides
				//string windowSequence = 

				string qualSeq = "";
				if(qualSequence.size() > 0){
					qualSeq = qualSequence.substr(readWindowStart, readWindowEnd-readWindowStart);
				}

				indexWindow(al, readWindowStart, readWindowEnd, contigWindowStart, contigWindowEnd, readSequence.substr(readWindowStart, readWindowEnd-readWindowStart), qualSeq);
			}
			*/
			

			//for(const auto& breakPoint : breaking_points_){
			//	_logFile << breakPoint.first << " " << breakPoint.second << endl;
			//}
		}

		//void indexWindow(const Alignment& al, size_t readWindowStart, size_t readWindowEnd, size_t contigWindowStart, size_t contigWindowEnd, const string& windowSequence, const string& windowQualities){
		void indexWindow(const Alignment& al, u_int64_t windowIndex, u_int32_t posStart, u_int32_t posEnd, const string& windowSequence, const string& windowQualities){

			#pragma omp critical(indexWindow)
			{
				
				//size_t contigWindowIndex = contigWindowStart / _windowLength;
				vector<Window>& windows = _contigWindowSequences[al._contigIndex][windowIndex];
				
				/*
				if(contigWindowIndex == 0){
					//_logFile << contigWindowStart << endl;
					_logFile << windowSequence << endl;
					//_logFile << windows.size() << endl;
					_logFile << readWindowStart << " " << readWindowEnd << endl;
					getchar();
					//_logFile << windowQualities << endl;
				}
				*/
				
				bool interrupt = false;
				if(_contigPolisher._maxWindowCopies == 0 || windows.size() < (_contigPolisher._maxWindowCopies-1)){


					windows.push_back({new DnaBitset2(windowSequence), windowQualities, posStart, posEnd, al.score()});

					/*
					if(al._contigIndex == 1 && windowIndex == 2){
						_logFile << "1111111111111111" << endl;
						_logFile  << windowSequence << endl;
					}
					if(al._contigIndex == 1 && windowIndex == 3){
						_logFile << "2222222222222222" << endl;
						_logFile  << windowSequence << endl;
					}
					*/

					/*
					if(windowSequences.size() > 1){
						static EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0);

						char* dnaStrModel = windowSequences[0]->to_string();
						EdlibAlignResult result = edlibAlign(dnaStrModel, strlen(dnaStrModel), windowSequence.c_str(), windowSequence.size(), config);

						char* cigar;

						if (result.status == EDLIB_STATUS_OK) {
							cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
							_logFile << strlen(cigar) << endl;
							free(cigar);
						} else {
							_logFile << "Invalid edlib results" << endl;
							exit(1);
						}

						//_logFile << cigar << endl;

						edlibFreeAlignResult(result);
						free(dnaStrModel);

					}
					*/


					interrupt = true;
				}

				
				if(!interrupt){

					float score = al.score();

					u_int64_t incompleteWindowIndex = -1;
					//u_int64_t minWindowSize = -1;

					u_int64_t largerDistanceWindow = 0;

					for(size_t i=0; i<windows.size(); i++){

						const Window& window = windows[i];
						u_int64_t distance = abs(((long)window._sequence->m_len) - ((long)_windowLength));

						if(distance > _contigPolisher._windowLengthVariance){
							if(distance > largerDistanceWindow){
								largerDistanceWindow = distance;
								incompleteWindowIndex = i;
							}
						}
					}


					//u_int64_t distance = abs(((long)windowSequence.size()) - ((long)_windowLength));

					if(incompleteWindowIndex != -1){
						Window& window = windows[incompleteWindowIndex];
						delete window._sequence;
						windows[incompleteWindowIndex] = {new DnaBitset2(windowSequence), windowQualities, posStart, posEnd, score};
					}
					else{
						//u_int64_t largestAligmenet = 0;
						//u_int64_t largestAligmenetIndex = 0;

						/*
						size_t largerWindowIndex = 0;
						u_int64_t largerDistanceWindow = 0;

						for(size_t i=0; i<windows.size(); i++){

							const Window& window = windows[i];
							u_int64_t distance = abs(((long)window._sequence->m_len) - ((long)_windowLength));

							if(distance > largerDistanceWindow){
								largerDistanceWindow = distance;
								largerWindowIndex = i;
							}
						}


						u_int64_t distance = abs(((long)windowSequence.size()) - ((long)_windowLength));

						if(distance < largerDistanceWindow){
							Window& window = windows[largerWindowIndex];
							delete window._sequence;
							windows[largerWindowIndex] = {new DnaBitset2(windowSequence), windowQualities, posStart, posEnd};
						}
						*/

						

        				static float maxVal = std::numeric_limits<float>::max();
						size_t largerWindowIndex = 0;
						//u_int64_t largerDistanceWindow = 0;
						float lowestScore = maxVal;
						
						for(size_t i=0; i<windows.size(); i++){

							const Window& window = windows[i];
							//u_int64_t distance = abs(((long)window._sequence->m_len) - ((long)_windowLength));
							//float score =

							if(window._score < lowestScore){
								lowestScore = window._score;
								largerWindowIndex = i;
							}
						}

						//u_int64_t distance = abs(((long)windowSequence.size()) - ((long)_windowLength));

						if(score > lowestScore){
							Window& window = windows[largerWindowIndex];
							delete window._sequence;
							windows[largerWindowIndex] = {new DnaBitset2(windowSequence), windowQualities, posStart, posEnd, score};
						}

					}



				} 




			}
			//_logFile << contigWindowIndex << endl;


			

			/*
			if(contigWindowIndex == _contigWindowSequences[al._contigIndex].size()-1){

				for(size_t i=0; i<windowSequences.size(); i++){
					_logFile << windowSequences[i]->m_len << " ";
				}
				_logFile << endl;
				//_logFile << contigWindowEnd-contigWindowStart << endl;
				//_logFile << windowSequence << endl;
			}
			*/
		}


	};

	struct CorrectedWindow{
		size_t _windowIndex;
		DnaBitset2* _correctedSequence;
		bool _success;
	};

	static bool CorrectedWindowComparator (const CorrectedWindow& p1, const CorrectedWindow& p2){
		return p1._windowIndex < p2._windowIndex;
	}


	unordered_map<u_int32_t, vector<CorrectedWindow>> _currentContigs;



	void performCorrection(gzFile& outputContigFile, const int& pass){
		
		u_int64_t checksum = 0;

		vector<u_int32_t> contigIndexes;
		for(auto& it : _contigWindowSequences){
			contigIndexes.push_back(it.first);
		}
		//vector<u_int32_t> contigIndexes;
		//for(auto& it : _contigWindowSequences){
		//	contigIndexes.push_back(it.first);
		//}


		size_t i = 0;
		size_t windowIndex = 0;
		//unordered_map<string, vector<vector<Window>>> _contigWindowSequences;

		#pragma omp parallel num_threads(_nbCores)
		{

			bool isEOF = false;
			size_t contigIndexLocal;
			size_t windowIndexLocal;
			std::unique_ptr<spoa::AlignmentEngine> alignmentEngine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -4);
			

			while(true){


				#pragma omp critical(lala)
				{

					if(i >= _contigWindowSequences.size()){
						isEOF = true;
					}

					if(!isEOF){

						u_int32_t contigIndex = contigIndexes[i];
						//cout << "Contig: " << contigIndex << endl;
						//cout << windowIndex << " " << _contigWindowSequences[contigIndex].size() << endl;
						if(windowIndex >= _contigWindowSequences[contigIndex].size()){
							i += 1;
							windowIndex = 0;
							//cout << "inc contig index" << endl;
						}

						if(i >= _contigWindowSequences.size()){
							isEOF = true;
						}


						if(!isEOF){
							contigIndex = contigIndexes[i];
							contigIndexLocal = contigIndex;
							windowIndexLocal = windowIndex;
							windowIndex += 1;
						}
					}
					
				}

				if(isEOF) break;
				//functorSub(read);

				//#pragma omp critical(lala)
				//{
				//	cout << contigIndexLocal << ": " << windowIndexLocal << "/" << _contigWindowSequences[contigIndexLocal].size() << endl;
				//	getchar();
				//}

				//const string& contigName = contigIndexes[contigIndexLocal];

				//cout << _contigWindowSequences.size() << " " << contigIndexLocal << endl;
				//cout << "lala: " << contigIndexLocal << " " << endl;
				

				vector<Window>& sequences = _contigWindowSequences[contigIndexLocal][windowIndexLocal];

				//cout << windowIndexLocal << " " << sequences.size() << endl;
				//u_int64_t nbCorrectedWindows = 0;
				//vector<DnaBitset2*> correctedWindows(windows.size());
			
					
				
				u_int64_t wStart = windowIndexLocal*_windowLength;
				u_int64_t wEnd = min(_contigSequences[contigIndexLocal].size(), (size_t)(wStart+_windowLength));
				string contigOriginalSequence = _contigSequences[contigIndexLocal].substr(wStart, wEnd-wStart);
				bool isLastWindow = (windowIndexLocal == _contigWindowSequences[contigIndexLocal].size()-1);

				if(sequences.size() < 2){

					for(size_t i=0; i<sequences.size(); i++){ 
						delete sequences[i]._sequence;
					}

					addCorrectedWindow(false, new DnaBitset2(contigOriginalSequence), contigIndexLocal, windowIndexLocal, outputContigFile, pass);	
					//correctedWindows[w] = new DnaBitset2(contigOriginalSequence);
					//_logFile << "No sequences for window" << endl;
					continue;
				}
				


				std::sort(sequences.begin(), sequences.end(), [&](const Window& lhs, const Window& rhs) {
					return lhs._posStart < rhs._posStart; });

				//_logFile << sequences.size() << endl;
				//_logFile << "1" << endl;
				spoa::Graph graph{};

				string backboneQuality = "";
				for(size_t i=0; i<contigOriginalSequence.size(); i++){
					backboneQuality += '!';
				}

				graph.AddAlignment(
					spoa::Alignment(),
					contigOriginalSequence.c_str(), contigOriginalSequence.size(),
					backboneQuality.c_str(), backboneQuality.size()
				);

				u_int32_t offset = 0.01 * contigOriginalSequence.size();


				for(size_t i=0; i<sequences.size(); i++){ 

					//size_t i = order[ii];
					const Window& window = sequences[i];
					//const DnaBitset2* dna = variant._sequence; //sequenceCopies[s._sequenceIndex];
					char* dnaStr = window._sequence->to_string();

					//if(windowIndexLocal == 24){
					//	cout << string(dnaStr) << endl;
					//}

					spoa::Alignment alignment;
					if (window._posStart < offset && window._posEnd > contigOriginalSequence.size() - offset) {
						alignment = alignmentEngine->Align(
							dnaStr, strlen(dnaStr),
							graph);
					} else {
						std::vector<const spoa::Graph::Node*> mapping;
						auto subgraph = graph.Subgraph(
							window._posStart,
							window._posEnd,
							&mapping);
						alignment = alignmentEngine->Align(
							dnaStr, strlen(dnaStr),
							subgraph);
						subgraph.UpdateAlignment(mapping, &alignment);
					}
					
					if (window._quality.size() == 0) {
						graph.AddAlignment(
							alignment,
							dnaStr, strlen(dnaStr));
					} else {
						graph.AddAlignment(
							alignment,
							dnaStr, strlen(dnaStr),
							window._quality.c_str(), window._quality.size());
					}
					
					free(dnaStr);
					delete window._sequence;
				}


				vector<u_int32_t> coverages;
				string correctedSequence = graph.GenerateConsensus(&coverages);

				correctedSequence = trimConsensus(correctedSequence, coverages, sequences.size(), isLastWindow);

				//for(char letter : correctedSequence){
				//	checksum += letter;
				//}

				//ofstream file("/pasteur/appa/homes/gbenoit/appa/run/correction/test_deterministic/test_humanO1_4/loulou.fasta");
				//file << ">lala" << endl;
				//file << correctedSequence << endl;
				//file.close();
				//cout << windowIndexLocal << endl;
				//getchar();

				addCorrectedWindow(true, new DnaBitset2(correctedSequence), contigIndexLocal, windowIndexLocal, outputContigFile, pass);

				//correctedWindows[w] = new DnaBitset2(correctedSequence);

			}
		}
		
	}

	string trimConsensus(const string& correctedSequence, const vector<u_int32_t>& coverages, const int& nbSequences, const bool& isLastWindow){

		string trimmedSequence = "";
		uint32_t average_coverage = nbSequences / 2;

		while(true){


			int32_t begin = 0, end = correctedSequence.size() - 1;
			for (; begin < static_cast<int32_t>(correctedSequence.size()); ++begin) {
				if (coverages[begin] >= average_coverage) {
					break;
				}
			}
			for (; end >= 0; --end) {
				if (coverages[end] >= average_coverage) {
					break;
				}
			}

			if (begin >= end) {
				//fprintf(stderr, "[racon::Window::generate_consensus] warning: "
				//	"contig %lu might be chimeric in window %u!\n", id_, rank_);
			} else {
				trimmedSequence = correctedSequence.substr(begin, end - begin + 1);
			}

			if(isLastWindow) break;
			if(trimmedSequence.size() > _windowLength*0.8) break;
			
			average_coverage += 1;

			if(average_coverage > nbSequences) return correctedSequence;
		}


		return trimmedSequence;
	}


	void addCorrectedWindow(bool success, DnaBitset2* seq, size_t contigIndexLocal, size_t windowIndexLocal, gzFile& outputContigFile, const int& pass){

		
		#pragma omp critical(addCorrectedWindow)
		{
			
			//cout << "Add corrected window: " << windowIndexLocal << " " << seq->m_len << endl;

			_currentContigs[contigIndexLocal].push_back({windowIndexLocal, seq, success});

			if(_currentContigs[contigIndexLocal].size() == _contigWindowSequences[contigIndexLocal].size()){
				dumpCorrectedContig(contigIndexLocal, outputContigFile, pass);
			}
		}
		
	}

	void dumpCorrectedContig(const u_int32_t& contigIndex, gzFile& outputContigFile, const int& pass){


		u_int64_t nbCorrectedWindows = 0;
		vector<CorrectedWindow>& correctedWindows = _currentContigs[contigIndex];

		for(const CorrectedWindow& correctedWindow : correctedWindows){
			if(correctedWindow._success){
				nbCorrectedWindows += 1;
			}
		}
		

		if(nbCorrectedWindows > 0){

			std::sort(correctedWindows.begin(), correctedWindows.end(), CorrectedWindowComparator);

			string contigSequence = "";
			//for(size_t w=0; w<correctedWindows.size(); w++){
			for(const CorrectedWindow& correctedWindow : correctedWindows){
				if(correctedWindow._correctedSequence == nullptr) continue;
				//if(correctedWindows[w] == nullptr) continue;
				char* seq = correctedWindow._correctedSequence->to_string();
				contigSequence += string(seq);
				free(seq);

				//cout << contigSequence << endl;
				//cout << contigSequence.size() << endl;
				//getchar();
				//delete correctedWindows[w];
			}


			u_int64_t length = contigSequence.size();
			bool isValid = true;

				
			if(_contigCoverages[contigIndex] <= 1){
				isValid = false;
			}
			else if(length < _minContigLength){
				isValid = false;
			}
			else if(length < 7500 && _contigCoverages[contigIndex] < 4){
				isValid = false;
			}

			//else if(length < 7500 && _contigCoverages[contigIndex] < 3){
			//	isValid = false;
			//}
			


			if(isValid){

				string header = _contigHeaders[contigIndex];
				
				if(_useMetamdbgHeaderStyle && pass > 0){
					
					string originalHeader = header;
					header = Utils::split(header, '_')[0]; //Remove curator subpath suffix

					bool isCircular = false;
					//if(header.find("rc") != string::npos){
					//	string h = header.substr(0, header.size()-2);
					//	header = h + "_" + to_string(_contigCoverages[contigIndex]) + "x_rc";
					//}
					//else 
					if(header[header.size()-1] == 'c'){
						isCircular = true;
						//string h = header.substr(0, header.size()-1);
						//header = h + " length=" + to_string(contigSequence.size()) + " coverage=" + to_string(_contigCoverages[contigIndex]) + " circular=yes";
					}

					header.pop_back(); //remove circular indicator
					header.erase(0, 3); //remove "ctg"
					//u_int32_t contigIndex = stoull(header);

					//else{
					//	string h = header.substr(0, header.size()-1);
					//	header = h + " length=" + to_string(contigSequence.size()) + " coverage=" + to_string(_contigCoverages[contigIndex]) + " circular=no";
					//}
					header = Utils::createContigHeader(_outputContigIndex, contigSequence.size(), _contigCoverages[contigIndex], isCircular);
					
					_file_contigHeaders << originalHeader << "\t" << _outputContigIndex << endl;
				}
				

				//cout << "Polish contig: " << contigSequence.size() << endl;
				header = ">" + header + '\n';// ">ctg" + to_string(contigIndex) + '\n';
				//header += '\n';
				gzwrite(outputContigFile, (const char*)&header[0], header.size());
				contigSequence +=  '\n';
				gzwrite(outputContigFile, (const char*)&contigSequence[0], contigSequence.size());
				//_logFile << _contigHeaders[contigIndex] << " " << contigSequence.size() << endl;
				
				_outputContigIndex += 1;

			}



		}

		for(CorrectedWindow& correctedWindow : correctedWindows){
			if(correctedWindow._correctedSequence == nullptr) continue;
			delete correctedWindow._correctedSequence;
		}
		_currentContigs.erase(contigIndex);
	}



	int32_t getMappableLength(const AlignmentBounds& alignmentBounds){

		int64_t queryLength = alignmentBounds._queryLength;
		int64_t queryStart = alignmentBounds._queryStart;
		int64_t queryEnd = alignmentBounds._queryEnd;
		
		int64_t targetLength = alignmentBounds._referenceLength;
		int64_t targetStart = alignmentBounds._referenceStart;
		int64_t targetEnd = alignmentBounds._referenceEnd;
		bool isReversed = alignmentBounds._isReversed;


		int64_t tl3 = 0;
		int64_t tl5 = 0;
		int64_t ext3 = 0;
		int64_t ext5 = 0;

        if(isReversed){
            tl5 = targetLength - targetEnd;
            tl3 = targetStart;
		}
        else{
            tl5 = targetStart;
            tl3 = targetLength - targetEnd;
		}

        if(queryStart < tl5){
            ext5 = queryStart;
		}
        else{
            ext5 = tl5;
		}

        if (queryLength - queryEnd < tl3){
            ext3 = queryLength - queryEnd;
		}
        else{
            ext3 = tl3;
		}

		//cout << "Ext5: " << ext5 << endl;
		//cout << "Ext3: " << ext3 << endl;
		return queryEnd - queryStart + ext5 + ext3;
	}


};	

#endif 


