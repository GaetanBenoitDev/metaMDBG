
/*

- flag isMetaMDBG, l'utiliser partout ou on créé des header metaMDBG (trimCOntig, COntigDerep)

- enelever le bordel de _windowPositionOffset
- contig header desactiver pendant le debuggage
- quand on charche le successeur d'un read dans le graph, il est possible qu'on ne choississe pas le meilleur successeur (genre le read d'une autre strain), si quand on effecture un le minimap2 alignment, on se rend compte qu'il y a un petit overhang entre les deux reads (< 200), on pourrait essayer d'autrees successeurs

- read contigStart/End, fuat-il rajouter la length non alignée, pour tirer les alignement par ordre du contig aussi

- on pourrait rafiner les bord des contigs en ajoutant les reads qui mappent au bout puis une dereplication massive des overlaps, mais risqué ptet
- contig derep: gros peak de mémoire vu qu'on utilise le meme -I 8G que pour le mapping en minimizer-space, essayer de fix le -I XG de 1 à 4

- mode --genomic: discard all kminmer since once in first pass
- checkpoint system

- ont read correction: moins bon resultat peut venir de la correction final via le chaining plutot que l'alignement (a test sur ZymoFecal)
- en mode ont on eneleve vraiment bcp de minimizer reptitive, pas sur de ça (check sur le soil aussi)

- RepeatRemover: baisser le 50 à 10 sourceUnitig
- RepeatRemover: si on a aucun source unitig (nbMinimizer > 50), il faudrait au moins prendre le plus grand par defaut, pour que chaque contig soit checké
- RepeatRemover: comment gerer l'endroit ou ont split les contig? (Actuellement on laisse un overlap mais il faudrait enlever l'overlap pour le fragment le plus petit), on peut aussi limité la taille des contigs a derepliquer (max 30000 par exemple, le but est d'enlever les erreurs)
- isCircular checker: readMapping + rafinement des bord

- maintenant qu'on dereplique bien, plus besoin d'enlever des contigs juste par leur taille (3000 ou 2500)
*/

/*
- en test:
	- is circular: verifier qu'on enleve le circular flag si on a decouper un contig
	- ReadCorrection ont: performPoaCorrection5: remettre le chaining alignment pour le poa en low-density

- getBestAlignment: acutellement on peut selection un tous petit overlap avec 0.99% identity, il faudrait faire une selection en deux temps, tout d'abord, minimum 5k overlap, puis si on a rien tout overlap length autorisé
- computeOverlapPath: plus que 1000?
- transitive reudction: possible de l'accelerer en reactivant le code de filtrage basé sur la length
- parallelisation
- sur le soil va falloir check la conso mémoire, check les potentiel memory leak
- _contigSequences: on peut le clear des qu'on a refait les alignment reads vs contigs (avant graphOverlapPath)

- Polishing:
	- racon: selection des N window: doivent provenir de la strain dominante en priorité (donc plus grand alignment)
	- racon: verifier erreurs sur le bord des window?
	- integrer abpoa ou virer du code de spoa qui n'est pas utiliser a la compil
	- mm2-fast, version rapide de minimap2: https://github.com/bwa-mem2/mm2-fast

- repeat post-process:
	- reessayer de detcter les repeat dans le polisher via de l'alignement
	- si le coverage splitting est validé, il faudra bien delimité les bounds des fragments avant de split
	- estimation coverage: peut etre mauvais pour les reads qui mappe sur plusieurs bord de contig (vu qu'on garde que le meilleur alignement)
	- fix le k et density si on conserve le decoupage en unitig

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

/*

- A Evaluer:
	- init contig length (on ne doit pas utiliser tout le premier read mais seulement la partir qui map sur le contig), mais attention quand on aligne le read suivant du coup (cas particulier)
	- Contig de taille 0 dans le fichier de sortie (on les enleve actuellement mais il faudrait check d'ou il viennent a la base)
	- cas particulier: petit contig dont le overlap path est de taille 1 (donc un read qui chevauche le petit contig), on a fait un cas particulier dans le convert path to sequence
	- attention: gestion de l'erreur si minimap2 crash + evaluation peak memory

*/


#ifndef MDBG_METAG_CONTIGPOLISHER
#define MDBG_METAG_CONTIGPOLISHER

#include "../Commons.hpp"
#include "../utils/edlib.h"
#include "../utils/spoa/include/spoa/spoa.hpp"
#include "../utils/DnaBitset.hpp"
//#include "../readSelection/BaseAligner.hpp"
#include "ContigDerep.hpp"
#include "../utils/minimap2/minimap.h"
#include "../utils/minimap2/mmpriv.h"

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
	bool _isMetaMDBG;

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
		args::Flag arg_isMetaMDBG(groupOther, "", "Do not use qualities during correction", {ARG_IS_METAMDBG}, args::Options::Hidden);
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

		_isMetaMDBG = false;
		if(arg_isMetaMDBG){
			_isMetaMDBG = true;
		}

		_circularize = false;
		//if(arg_useCirculize){
		//	_circularize = true;
		//}



		_windowLength = 500;
		_windowPositionOffset = 0;
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

		_outputContigIndex = 0;
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

		//_partitionNbReads[28] = 10000;
		//processPartition(28);
		//exit(1);

		for(size_t i=0; i<_nbPartitions; i++){
			_partitionNbReads[i] = 10000;
			processPartition(i);
		}

		gzclose(_usedReadFile);

		//for(size_t i=0; i<_nbPartitions; i++){
		//	_partitionNbReads[i] = 10000;
			//detectWeakRepeats2(i);
		//}

		Logger::get().debug() << "";
		Logger::get().debug() << "Polishing contigs (pass 1)";
		_windowPositionOffset = 0;

		for(size_t i=0; i<_nbPartitions; i++){

			
			string inputContigFilename = _readPartitionDir + "/" + to_string(i) + "_contigsCurated.gz";
			string outputContigFilename = _readPartitionDir + "/" + to_string(i) + "_contigsPolished.gz";
			gzFile outputContigFile = gzopen(outputContigFilename.c_str(),"wb");

			//_partitionNbReads[i] = 10000;
			polishPartition(i, inputContigFilename, outputContigFile);

			gzclose(outputContigFile);
		}


		Logger::get().debug() << "";
		Logger::get().debug() << "Polishing contigs (pass 2)";
		_windowPositionOffset = 0; //_windowLength/2;

		const string& outputContigFilename_polished = _readPartitionDir + "/contigs_polished.fasta.gz";
		gzFile outputContigFile_polished = gzopen(outputContigFilename_polished.c_str(),"wb");

		for(size_t i=0; i<_nbPartitions; i++){


			string inputContigFilename = _readPartitionDir + "/" + to_string(i) + "_contigsPolished.gz";

			polishPartition(i, inputContigFilename, outputContigFile_polished);

		}

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
		Logger::get().debug() << "Trimming contigs";
		const string& outputContigFilename_trimmed = _readPartitionDir + "/contigs_trimmed.fasta.gz";
		trimContigs(outputContigFilename_derep, outputContigFilename_trimmed, minimapBatchSize);

		cout << "A REMETTRE: delete tmp dir for polishing" << endl;
		//fs::remove_all(_readPartitionDir);

		Logger::get().debug() << "";
		Logger::get().debug() << "Moving final contigs to destination";
		fs::rename(outputContigFilename_trimmed, _outputContigFilename);

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
			Utils::ContigHeader contigHeader = Utils::extractContigHeader(header);
			header = Utils::createContigHeader(contigHeader._contigIndex, contigSequenceTrimmed.size(), contigHeader._coverage, contigHeader._isCircular);

			string headerFasta = ">" + header + '\n';
			gzwrite(_parent._outputContigFile_trimmed, (const char*)&headerFasta[0], headerFasta.size());
			contigSequenceTrimmed +=  '\n';
			gzwrite(_parent._outputContigFile_trimmed, (const char*)&contigSequenceTrimmed[0], contigSequenceTrimmed.size());

		}
	};

	
	void derepContigs(const string& inputContigFilename, const string& outputContigFilename, const int& minimapBatchSize){

		auto start = high_resolution_clock::now();

		ContigDerep contigDerep(inputContigFilename, outputContigFilename, minimapBatchSize, _tmpDir, _readPartitionDir, 0.9, _minContigLength, _nbCores);
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
		//writeContigPartitions();

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
			_partitionNbReads[i] = 0;
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



	
	//typedef phmap::parallel_flat_hash_map<MinimizerType, vector<MinimizerPosition>, phmap::priv::hash_default_hash<MinimizerType>, phmap::priv::hash_default_eq<MinimizerType>, std::allocator<std::pair<MinimizerType, vector<MinimizerPosition>>>, 4, std::mutex> MinimizerPosMap;
	//MinimizerPosMap _minimizerIndex;


	void processPartition(u_int32_t partition){


		_currentPartition = partition;
		//if(_contigSequences.size() == 0) return;

		Logger::get().debug() << "";
		Logger::get().debug() << "\tCurating partition: " << _currentPartition << "/" << _nbPartitions;
		auto start = high_resolution_clock::now();

		if(_partitionNbReads[partition] == 0) return;


		//cout << "--------------" << endl;

		string readFilename = _readPartitionDir + "/" + to_string(partition) + "_reads.gz";
		//string contigFilename = _readPartitionDir + "/" + to_string(partition) + "_contigs.gz";
		string outputContigFilename = _readPartitionDir + "/" + to_string(partition) + "_contigsCurated.gz";

		//indexContigs();
		clearPass();
		
		
		_outputContigFilePartition = gzopen(outputContigFilename.c_str(),"wb");
		
		//loadContigs(_inputFilename_contigs, false);
		//parseAlignmentsGz(true);

		
		//cout << "Load contigs" << endl;
		Logger::get().debug() << "\tLoading contigs";
		loadContigs(_inputFilename_contigs, false);
		

		Logger::get().debug() << "\tAligning reads vs contigs";
		auto start2 = high_resolution_clock::now();
		//const string& partitionFilename = _readPartitionDir + "/part_" + to_string(partition) + ".gz";
		ReadParserParallel readParser2(readFilename, true, false, _nbCores);
		readParser2.parse(ReadsToContigsAlignFunctor(*this));
		Logger::get().debug() << "\tDone " << " (" << duration_cast<seconds>(high_resolution_clock::now() - start2).count() << "s)";





		Logger::get().debug() << "\tIndexing reads"; 
		auto start3 = high_resolution_clock::now();
		ReadParserParallel readParser(readFilename, true, false, _nbCores);
		readParser.parse(ReadLoaderFunctor(*this));
		Logger::get().debug() << "\tDone " << " (" << duration_cast<seconds>(high_resolution_clock::now() - start3).count() << "s)";
		//indexReads(readFilename);


		//_alignments_isReversed.clear();
		//cout << "Actuellement, les minimizers sont stockés dans le sens aligné sur les contigs, mais on ne stock pas les seuqnece reads dans le sens aligné sur le contig mais dans le sens original: _parent._readIndex_to_sequence[readIndex] = new DnaBitset2(read._seq);" << endl;


		Logger::get().debug() << "\tTransforming contig sequence to overlap paths";
		auto start4 = high_resolution_clock::now();
		computeContigOverlapPaths();
		Logger::get().debug() << "\tDone " << " (" << duration_cast<seconds>(high_resolution_clock::now() - start4).count() << "s)";

		Logger::get().debug() << "\tConvert overlap paths to sequences";
		auto start5 = high_resolution_clock::now();
		convertOverlapPathToSequences();
		Logger::get().debug() << "\tDone " << " (" << duration_cast<seconds>(high_resolution_clock::now() - start5).count() << "s)";

		//if(fs::exists(_outputFilename_mapping)) fs::remove(_outputFilename_mapping);

		gzclose(_outputContigFilePartition);

		//cout << endl << endl;
		//cout << "Initial contigs:" << endl;
		//for(const auto& it : _contigSequences){
		//	cout << "\t" << it.first << ": " << it.second.size() << endl;
		//}

		Logger::get().debug() << "\tWriting used reads";
		auto start6 = high_resolution_clock::now();
		writeUsedReads();
		Logger::get().debug() << "\tDone " << " (" << duration_cast<seconds>(high_resolution_clock::now() - start6).count() << "s)";

		mm_idx_destroy(_minimap2Index);
		

		//collectOverlappingReads(outputContigFilename, readFilename);

		Logger::get().debug() << "\tPartition " << partition << " done " << " (" << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s)";
		//cout << "done" << endl;
		//getchar();
	}

	
	void collectOverlappingReads(const string& contigFilename, const string& readFilename){

		const string& alignFilename = _readPartitionDir + "/align.paf.gz";
		string command = "minimap2 -v 0 -m 500 -t " + to_string(_nbCores) + " -x " + _minimap2Preset_map + " " + contigFilename + " " + readFilename;
		cout << command << endl;
		Utils::executeMinimap2(command, alignFilename);

		indexContigName(contigFilename);
		loadContigs(contigFilename, true);
		loadAlignments(alignFilename, true, true, false);

		for(const auto& it : _alignments_readsVsContigs){

			const u_int32_t& contigIndex = it.first;
			const auto& alignments = it.second;

			for(const auto& it2: alignments){

				const ReadType& readIndex = it2.first;
				const Alignment& alignment = it2.second;

				AlignmentBounds alignmentBounds;
				alignmentBounds._queryStart = alignment._readStart;
				alignmentBounds._queryEnd = alignment._readEnd;
				alignmentBounds._referenceStart = alignment._contigStart;
				alignmentBounds._referenceEnd = alignment._contigEnd;
				alignmentBounds._queryLength = alignment._readLength;
				alignmentBounds._referenceLength = alignment._contigLength;
				alignmentBounds._isReversed = alignment._strand;

				if(isValidAlignment(alignmentBounds, 0)){
					cout << _contigHeaders[contigIndex] << "\t" << alignmentBounds._referenceStart << "\t" << alignmentBounds._referenceEnd << "\t" << _contigSequences[contigIndex].size() << "\t\t" << readIndex << " " << alignmentBounds._isReversed << endl;
					//alignmentBounds.printAll();
				}


				//if(_contigHeaders[contigIndex] == "ctg52409l_1" && alignmentBounds._referenceStart > 20000){
				//	cout << "b " << _contigHeaders[contigIndex] << "\t" << alignmentBounds._referenceStart << "\t" << alignmentBounds._referenceEnd << "\t" << _contigSequences[contigIndex].size() << "\t\t" << readIndex << "\t" << alignmentBounds._isReversed << "\t" << alignment._readLength << endl;
					
				//}
			}
		}

		getchar();
	}
	

	void writeUsedReads(){

		for(const ReadType& readIndex : _usedReadIndexes){

			char* dnaStr = _readIndex_to_sequence[readIndex]->to_string();
			string readSeq = string(dnaStr, strlen(dnaStr));
			free(dnaStr);

			string header = ">" + to_string(readIndex) + '\n';
			//header += '\n';
			gzwrite(_usedReadFile, (const char*)&header[0], header.size());
			string readSequence = readSeq +  '\n';
			gzwrite(_usedReadFile, (const char*)&readSequence[0], readSequence.size());

		}
		
	}

	mm_idx_t* _minimap2Index;

	/*
	void indexReads(const string& readFilename){


		string preset = "ava-ont";
		mm_idxopt_t iopt;
		mm_mapopt_t mopt;
		mm_set_opt(0, &iopt, &mopt); //"ava-ont"
		mm_set_opt(preset.c_str(), &iopt, &mopt); //"ava-ont"
		iopt.batch_size = 0x7fffffffffffffffL; //always build a uni-part index
		
		mm_idx_reader_t *r = mm_idx_reader_open(readFilename.c_str(), &iopt, 0);
		_minimap2Index = mm_idx_reader_read(r, _nbCores);
		//mm_idx_t* mi = minimap2index(iopt.w, iopt.k, iopt.flag&1, iopt.bucket_bits, target);
		mm_mapopt_update(&mopt, _minimap2Index);
		//mopt.mid_occ = 1000; // don't filter high-occ seeds

		mm_idx_reader_close(r);

	}

	string getReadSequence(const ReadType& readIndex){

		static string nucletoides = "ACGTN";

		int length = _minimap2Index->seq[readIndex].len;
		
		int en = length;
		int st = 0;
		char *s;
		s = (char*)malloc(en - st + 1);

		int len = mm_idx_getseq(_minimap2Index, readIndex, st, en, (uint8_t*)s);
		for (size_t i = 0; i < len; ++i){
			s[i] = nucletoides[(uint8_t)s[i]];
		}
		s[len] = 0;

		string seq(s);
		//cout << s << endl;
		free(s);

		return seq;
	}

	*/

	unordered_map<u_int32_t, u_int32_t> _contigIndex_to_nbPaths;

	void convertOverlapPathToSequences(){

		//vector<OverlapGraphPath> paths;
		//vector<u_int32_t> contigIndexes;
		for(const OverlapGraphPath& path : _contigOverlapPaths){
			_contigIndex_to_nbPaths[path._contigIndex] += 1;
			//paths.push_back(path);
		}


		std::sort(_contigOverlapPaths.begin(), _contigOverlapPaths.end(), [](const OverlapGraphPath& a, const OverlapGraphPath& b){
			u_int32_t length1 = a._referenceEnd - a._referenceStart;
			u_int32_t length2 = b._referenceEnd - b._referenceStart;
			return length1 > length2;
		});

		size_t i = 0;


		#pragma omp parallel num_threads(_nbCores)
		{

			ConvertOverlapPathToSequenceFunctor functorSub(*this, _print_debug);//(functor);
			bool isEOF = false;
			OverlapGraphPath path;
			//u_int32_t contigIndex = -1;

			while(true){


				#pragma omp critical(lala)
				{

					if(i >= _contigOverlapPaths.size()){
						isEOF = true;
					}

					if(!isEOF){

						//contigIndex = contigLengths[i]._contigIndex;
						path = _contigOverlapPaths[i];
						i += 1;
						//cout << "Contig: " << contigIndex << endl;
						//cout << windowIndex << " " << _contigWindowSequences[contigIndex].size() << endl;
						//if(windowIndex >= _contigWindowSequences[contigIndex].size()){
						//	i += 1;
						//	windowIndex = 0;
							//cout << "inc contig index" << endl;
						//}

						//if(i >= _contigWindowSequences.size()){
						//	isEOF = true;
						//}


						//if(!isEOF){
						//	contigIndex = contigIndexes[i];
						//	contigIndexLocal = contigIndex;
						//	windowIndexLocal = windowIndex;
						//	windowIndex += 1;
						//}
					}
					
				}

				if(isEOF) break;

				functorSub(path);
				//addCorrectedWindow(true, new DnaBitset2(correctedSequence), contigIndexLocal, windowIndexLocal);


			}
		}

		/*
		//_contigIndex_to_overlapPaths.clear();

		#pragma omp parallel num_threads(_nbCores)
		{

			ConvertOverlapPathToSequenceFunctor functorSub(*this, _print_debug);//(functor);

            #pragma omp for
            for(size_t i=0; i < _contigOverlapPaths.size(); i++){


				//cout << "lala: " << contigIndexes[i] << endl;
                //functorSub(_unitigs[i]);
				functorSub(_contigOverlapPaths[i]);
            }
			
		}
		*/
	}


	class ConvertOverlapPathToSequenceFunctor {

		public:

		ContigPolisher& _parent;
		bool _print_debug;
		mm_tbuf_t* _tbuf;

		ConvertOverlapPathToSequenceFunctor(ContigPolisher& parent, bool print_debug) : _parent(parent){
			_print_debug = print_debug;
			_tbuf = mm_tbuf_init();
		}

		//ComputeContigOverlapPathFunctor(const ComputeContigOverlapPathFunctor& copy) : _parent(copy._parent){
		//}

		~ConvertOverlapPathToSequenceFunctor(){
			mm_tbuf_destroy(_tbuf);
		}

		void operator () (const OverlapGraphPath& overlapPath) {


			string contigSeq = "";

			if(overlapPath._path.size() == 0){
			}
			else if(overlapPath._path.size() == 1){

				u_int32_t readIndex = overlapPath._path[0];
				const Alignment& al = _parent._alignments_readsVsContigs[overlapPath._contigIndex][readIndex];

				char* dnaStr = _parent._readIndex_to_sequence[readIndex]->to_string();
				string readSeq = string(dnaStr, strlen(dnaStr));
				free(dnaStr);

				readSeq = readSeq.substr(al._readStart, al._readEnd-al._readStart);

				//if(al._strand){
				//	Utils::toReverseComplement(readSeq);
				//}

				contigSeq = readSeq;
			}
			else{
					
				int32_t firstReadLength = 0;
				int32_t lastReadLength = 0;

				for(size_t i=0; i<overlapPath._path.size()-1; i++){
					u_int32_t readIndex1 = overlapPath._path[i];
					u_int32_t readIndex2 = overlapPath._path[i+1];

					const Alignment& alignment1 = _parent._alignments_readsVsContigs[overlapPath._contigIndex][readIndex1];
					const Alignment& alignment2 = _parent._alignments_readsVsContigs[overlapPath._contigIndex][readIndex2];
					if(_print_debug) cout << endl << i << endl;// ": " << readIndex1 << " " << _alignments_readsVsContigs[readIndex1] << " -> "
					if(_print_debug) alignment1.printAll();
					if(_print_debug) alignment2.printAll();
					
					//bool isReversed1 = alignment1._strand;
					//bool isReversed2 = alignment2._strand;

					char* dnaStr1 = _parent._readIndex_to_sequence[readIndex1]->to_string();
					string readSeq1 = string(dnaStr1, strlen(dnaStr1));
					free(dnaStr1);

					char* dnaStr2 = _parent._readIndex_to_sequence[readIndex2]->to_string();
					string readSeq2 = string(dnaStr2, strlen(dnaStr2));
					free(dnaStr2);

					AlignmentBounds alignmentBounds = _parent.computeReadOverlap(alignment1, alignment2, _tbuf);
					
					string endSeq2 = readSeq2.substr(alignmentBounds._queryEnd, readSeq2.size()-alignmentBounds._queryEnd);

					if(i == 0){
						firstReadLength = readSeq1.size();
						//cout << readSeq1_tmp << endl;
					}
					lastReadLength = readSeq2.size();


					//print(str(i) + ": " + nodeName1 + " " + str(name_to_alignment[nodeName1]) + " -> " + str(name_to_alignment[nodeName2]))

					/*
					const Alignment& readOverlap = getOverlap(readIndex1, readIndex2);
					if(_print_debug) cout << "Read overlap: " << endl;
					if(_print_debug) readOverlap.print();


					//cout << (_readIndex_to_sequence.find(readIndex1) != _readIndex_to_sequence.end()) << endl;
					//cout << (_readIndex_to_sequence.find(readIndex2) != _readIndex_to_sequence.end()) << endl;




					u_int32_t readSeq1_start;
					u_int32_t readSeq1_end;
					u_int32_t readSeq2_start;
					u_int32_t readSeq2_end;

					if (alignment1._readIndex == readOverlap._readIndex){
						readSeq1_start = readOverlap._readStart;
						readSeq1_end = readOverlap._readEnd;
						readSeq2_start = readOverlap._contigStart;
						readSeq2_end = readOverlap._contigEnd;
					}
					else{
						readSeq2_start = readOverlap._readStart;
						readSeq2_end = readOverlap._readEnd;
						readSeq1_start = readOverlap._contigStart;
						readSeq1_end = readOverlap._contigEnd;
					}

					//cout << "Expected: " << readSeq1_start << " " << readSeq1_end << " " << readSeq2_start << " " << readSeq2_end << endl;


					if (!isReversed1 && !isReversed2){
						//endSeq2 = readSeq2[readSeq2_end:]
						endSeq2 = readSeq2.substr(readSeq2_end, readSeq2.size()-readSeq2_end);
					}
					else if (!isReversed1 && isReversed2){
						endSeq2 = readSeq2.substr(0, readSeq2_start+1);
						Utils::toReverseComplement(endSeq2);
						//endSeq2 = readSeq2[:readSeq2_start]
						//endSeq2 = str(Seq(endSeq2).reverse_complement())
					}
					else if (isReversed1 && !isReversed2){
						//endSeq2 = readSeq2[readSeq2_end:]
						endSeq2 = readSeq2.substr(readSeq2_end, readSeq2.size()-readSeq2_end);
					}
					else if (isReversed1 && isReversed2){
						endSeq2 = readSeq2.substr(0, readSeq2_start+1);
						Utils::toReverseComplement(endSeq2);
						//endSeq2 = readSeq2[0:readSeq2_start]
						//endSeq2 = str(Seq(endSeq2).reverse_complement())
					}
					*/
					/*
					string readSeq1_tmp = readSeq1;
					if(isReversed1){
						Utils::toReverseComplement(readSeq1_tmp);
					}
					string readSeq2_tmp = readSeq2;
					if(isReversed2){
						Utils::toReverseComplement(readSeq2_tmp);
					}



					cout << "todo" << endl;
					//cout << isReversed1 << " " << isReversed2 << endl;
					AlignmentBounds alignmentBounds;// = _parent.computeReadOverlap(alignment1, alignment2, false, _parent.computeReadOverlapBounds(alignment1._readIndex, alignment2._readIndex, _minimizerChainer), _aligner);
					if(_print_debug) alignmentBounds.printAll();
					//cout << "Actual: " << alignmentBounds._queryStart << " " << alignmentBounds._queryEnd << " " << alignmentBounds._referenceStart << " " << alignmentBounds._referenceEnd << endl;
					//cout << isReversed1 << " " << isReversed2 << endl;
					//getchar();

					//cout << "Expected:" << endl;
					//cout << endSeq2.size() << endl;
					//cout << endSeq2 << endl;

					endSeq2 = readSeq2_tmp.substr(alignmentBounds._queryEnd, readSeq2_tmp.size()-alignmentBounds._queryEnd);
					
					//cout << "Actual:" << endl;
					//cout << endSeq2.size() << endl;
					//cout << endSeq2 << endl;
					
					//getchar();

					*/


					if (i == 0){
						//if (isReversed1){
						//	string readSeq1copy = readSeq1;
							//readSeq1copy = readSeq1copy.substr(alignment1._readStart, alignment1._readEnd-alignment1._readStart);
						//	Utils::toReverseComplement(readSeq1copy);
						//	contigSeq = readSeq1copy; //str(Seq(readSeq1).reverse_complement());
						//}
						//else{
						string readSeq1copy = readSeq1;
						//readSeq1copy = readSeq1copy.substr(alignment1._readStart, alignment1._readEnd-alignment1._readStart);
						contigSeq = readSeq1copy;
						//}

						const Alignment& firstAlignment = _parent._alignments_readsVsContigs[overlapPath._contigIndex][overlapPath._path[0]];
						
						int32_t beginingLengthToRemove = 0;
						if(firstAlignment._strand){
							beginingLengthToRemove = firstReadLength - firstAlignment._readEnd;
						}
						else{
							beginingLengthToRemove = firstAlignment._readStart;
						}
						if(beginingLengthToRemove > 0) contigSeq = contigSeq.substr(beginingLengthToRemove, contigSeq.size()-beginingLengthToRemove);

						if(_print_debug) cout << alignment1._readStart << " " << alignment1._readEnd << " " << firstReadLength << endl;
						if(_print_debug) cout << "Init contig seq: " << contigSeq.size() << endl;

					}
					
					u_int32_t contigOverHang = readSeq1.size() - alignmentBounds._referenceEnd;
					if(contigOverHang > 0){
						contigSeq = contigSeq.substr(0, contigSeq.size()-contigOverHang);
					}

					contigSeq += endSeq2;
					if(_print_debug) cout << "Contig length: " << contigSeq.size() << "   added: " << endSeq2.size() << endl;
					
					//getchar();
				}

				//cout << contigSeq.substr(0, 1000) << endl; 

				//cout << "lala" << endl;
				//cout << contigSeq.size() << endl;


				
				const Alignment& lastAlignment = _parent._alignments_readsVsContigs[overlapPath._contigIndex][overlapPath._path[overlapPath._path.size()-1]];

				int32_t endingLengthToRemove = 0;
				if(lastAlignment._strand){
					endingLengthToRemove = lastAlignment._readStart;
				}
				else{
					endingLengthToRemove = lastReadLength - lastAlignment._readEnd;
				}
				if(endingLengthToRemove > 0) contigSeq = contigSeq.substr(0, contigSeq.size()-endingLengthToRemove);

				//firstAlignment.printAll();
				//cout << beginingLengthToRemove << endl;

				//lastAlignment.printAll();
				//cout << endingLengthToRemove << endl;


				//cout << contigSeq.size() << endl;
				//cout << contigSeq.substr(0, 1000) << endl; 
			}



			if(contigSeq.size() <= _parent._windowLength) return;
			//return contigSeq

			#pragma omp critical(convertOverlapPathToSequence)
			{				
				
				for(size_t i=0; i<overlapPath._path.size(); i++){
					_parent._usedReadIndexes.insert(overlapPath._path[i]);
				}

				//contigSeq = contigSeq.substr(0, contigSeq.size()-20000);
				//contigSeq = contigSeq.substr(20000, contigSeq.size()-20000);
				if(_print_debug) cout << "Write contig: " << contigSeq.size() << " " << _parent._contigSequences[overlapPath._contigIndex].size() << endl;

				string header = _parent._contigHeaders[overlapPath._contigIndex];

				if(_parent._contigIndex_to_nbPaths[overlapPath._contigIndex] > 1){
					header.pop_back();
					header += 'l';
				}

				header = ">" + header + "_" + to_string(overlapPath._pathIndex) + '\n';// ">ctg" + to_string(contigIndex) + '\n';
				//header += '\n';
				gzwrite(_parent._outputContigFilePartition, (const char*)&header[0], header.size());
				string contigSequence = contigSeq +  '\n';
				gzwrite(_parent._outputContigFilePartition, (const char*)&contigSequence[0], contigSequence.size());
				//_logFile << _contigHeaders[contigIndex] << " " << contigSequence.size() << endl;
				//exit(1);
			}

		
		}



	};

	Alignment getOverlap(const ReadType& prevReadIndex, const ReadType& nextReadIndex){

		std::pair<ReadType, ReadType> key;

		if (prevReadIndex < nextReadIndex){
			key.first = prevReadIndex;
			key.second = nextReadIndex;
		}
		else{
			key.first = nextReadIndex;
			key.second = prevReadIndex;
		}

		return _alignments_readsVsReads[key];
	}
	
	struct ContigLength{
		u_int32_t _contigIndex;
		u_int32_t _contigLength;
	};

	void computeContigOverlapPaths(){
		
		vector<ContigLength> contigLengths;
		//cout << "single core here" << endl;
		//int nbCores = _nbCores;
		//if(_testAllVsAllAlignments) nbCores = 1;

		//vector<u_int32_t> contigIndexes;
		for(auto& it : _alignments_readsVsContigs){

			const u_int32_t& contigIndex = it.first;
			const u_int32_t& contigLength = _contigSequences[it.first].size();

			//if(_print_debug){
			//	if(_contigHeaders[contigIndex] != "ctg53061l") continue;
			//}
			
			contigLengths.push_back({contigIndex, contigLength});
		}

		std::sort(contigLengths.begin(), contigLengths.end(), [](const ContigLength& a, const ContigLength& b){
			return a._contigLength > b._contigLength;
		});


		size_t i = 0;

		int nbCores = _nbCores;
		if(_print_debug) nbCores = 1;

		#pragma omp parallel num_threads(nbCores)
		{

			ComputeContigOverlapPathFunctor functorSub(*this, _print_debug);//(functor);
			bool isEOF = false;
			u_int32_t contigIndex = -1;

			while(true){


				#pragma omp critical(lala)
				{

					if(i >= contigLengths.size()){
						isEOF = true;
					}

					if(!isEOF){

						contigIndex = contigLengths[i]._contigIndex;
						i += 1;
						//cout << "Contig: " << contigIndex << endl;
						//cout << windowIndex << " " << _contigWindowSequences[contigIndex].size() << endl;
						//if(windowIndex >= _contigWindowSequences[contigIndex].size()){
						//	i += 1;
						//	windowIndex = 0;
							//cout << "inc contig index" << endl;
						//}

						//if(i >= _contigWindowSequences.size()){
						//	isEOF = true;
						//}


						//if(!isEOF){
						//	contigIndex = contigIndexes[i];
						//	contigIndexLocal = contigIndex;
						//	windowIndexLocal = windowIndex;
						//	windowIndex += 1;
						//}
					}
					
				}

				if(isEOF) break;

				functorSub(contigIndex);
				//addCorrectedWindow(true, new DnaBitset2(correctedSequence), contigIndexLocal, windowIndexLocal);


			}
		}


		/*
		#pragma omp parallel num_threads(nbCores)
		{

			ComputeContigOverlapPathFunctor functorSub(*this, _print_debug);//(functor);

            //#pragma omp for
            for(size_t i=0; i < contigLengths.size(); i++){


				//cout << "lala: " << contigIndexes[i] << endl;
                //functorSub(_unitigs[i]);
				functorSub(contigIndexes[i]);
            }
			
		}
		*/
	}

	class ComputeContigOverlapPathFunctor {

		public:

		ContigPolisher& _parent;
		phmap::parallel_flat_hash_map<std::pair<u_int32_t, u_int32_t>, bool> _alignments_readsVsReads_cached;
		mm_tbuf_t* _tbuf;
		bool _print_debug;
		
		ComputeContigOverlapPathFunctor(ContigPolisher& parent, bool print_debug) : _parent(parent){
			_print_debug = print_debug;
			_tbuf = mm_tbuf_init();
		}

		//ComputeContigOverlapPathFunctor(const ComputeContigOverlapPathFunctor& copy) : _parent(copy._parent){
		//}

		~ComputeContigOverlapPathFunctor(){
			mm_tbuf_destroy(_tbuf);
		}

		void operator () (const u_int32_t& contigIndex) {

			_alignments_readsVsReads_cached.clear();

			//if(_parent._contigSequences[contigIndex].size() < 5100000) return;
			//if(_parent._contigHeaders[contigIndex] != "ctg5582l") return;

			if(_print_debug) cout << "Curating contig: " << _parent._contigSequences[contigIndex].size() << endl;

			vector<Alignment> alignments;

			for(auto& it : _parent._alignments_readsVsContigs[contigIndex]){

				const ReadType& readIndex = it.first;
				const Alignment& al = it.second;

				//Alignment2 al2 = {al._contigIndex, readIndex, al._strand, al._readStart, al._readEnd, al._contigStart, al._contigEnd, al._identity};
				alignments.push_back(al);
			}

			if(_print_debug) cout << "Sorting alignment: " << alignments.size() << endl;

			std::sort(alignments.begin(), alignments.end(), [](const Alignment& a, const Alignment& b){

				if(a._contigStart == b._contigStart){
					return a._contigEnd < b._contigEnd;
				}

				return a._contigStart < b._contigStart;
			});

			if(_print_debug) cout << "Create overlap graph" << endl;
			OverlapGraph* overlapGraph = createOverlapGraph(contigIndex, alignments);

			if(_print_debug) cout << "Transitive reduction" << endl;
			transitiveReduction2(overlapGraph);

			if(_print_debug) cout << "Repeat overlap cleaning" << endl;
			//cleanRepeatOverlaps(overlapGraph);

			//exportToGfa(overlapGraph, "/pasteur/appa/homes/gbenoit/zeus/tmp//assembly_graph.gfa");
			//exit(1);
			//cout << "Connected compojnents" << endl;
			//connectedComponents(overlapGraph);
			//for(const Alignment2 al : alignments){
			//	cout << al._contigStart << " " << al._contigEnd << endl;
				//getchar();
			//}

			if(_print_debug) cout << "Compute overlap paths" << endl;

			vector<OverlapGraphPath> paths;
			computeOverlapPaths(overlapGraph, alignments, paths);


			if(_print_debug) cout << "Most supported paths:" << endl;
			for (const OverlapGraphPath& path : paths){
				if(_print_debug) cout << "\t" << path._referenceStart << " " << path._referenceEnd << endl;
				//printPathSupport(overlapGraph, path);
				//getchar();
			}

			//order overlap path per reference position
			std::sort(paths.begin(), paths.end(), [](const OverlapGraphPath& a, const OverlapGraphPath& b){
				return a._referenceStart < b._referenceStart;
			});

			for(size_t i=0; i<paths.size(); i++){
				paths[i]._pathIndex = i;
			}

			//getchar();

			#pragma omp critical(computeContigOverlapPath)
			{

				

				for (const OverlapGraphPath& path : paths){
					_parent._contigOverlapPaths.push_back(path);
				}
			}
		
			delete overlapGraph;

		}



		struct SuccessorOverlapSorter{
			ReadType _readIndex;
			int32_t _estimatedOverlapLength;
		};

		OverlapGraph* createOverlapGraph(const u_int32_t& contigIndex, const vector<Alignment>& alignments){

			u_int64_t nbEdges = 0;

			OverlapGraph* graph = new OverlapGraph(contigIndex);
			
			for(size_t i=0; i<alignments.size(); i++){



				//getchar();

				const Alignment& alignment1 = alignments[i];

				//if(alignment1._contigStart > 412000){
				//	cout << "-------" << endl;
				//	cout << i << endl;
				//	alignment1.printAll();
				//	getchar();
				//}

				vector<SuccessorOverlapSorter> successorAlignemnts;

				for(size_t j=i+1; j<alignments.size(); j++){
					
					//cout << "\t" << j << endl;

					const Alignment& alignment2 = alignments[j];


				
					if(alignment2._contigStart > alignment1._contigEnd) break;
					if(!overlapOnTheReference(alignment1, alignment2)) continue;

					//if(alignment1._contigStart > 412000){
					//	cout << j << endl;
					//	alignment2.printAll();
					//}

					/*
					if(_testAllVsAllAlignments && alignment1._contigStart > 50000){
 

						//if(alignment1._readIndex != 4451) continue;
						//if(alignment2._readIndex != 1650) continue;
						 
						//if(alignment1._strand) continue;
						//if(alignment2._strand) continue;

						char* dnaStr1 = _parent._readIndex_to_sequence[alignment1._readIndex]->to_string();
						string seq1 = string(dnaStr1, strlen(dnaStr1));
						free(dnaStr1);

						char* dnaStr2 = _parent._readIndex_to_sequence[alignment2._readIndex]->to_string();
						string seq2 = string(dnaStr2, strlen(dnaStr2));
						free(dnaStr2);

						if(alignment1._strand){
							Utils::toReverseComplement(seq1);
						}
						
						if(alignment2._strand){
							Utils::toReverseComplement(seq2);
						}

						cout << endl;
						cout << endl;
						cout << endl;
						cout << endl;
						cout << endl;
						cout << alignment1._readIndex << " " << alignment2._readIndex << endl;
						alignment1.printAll();
						alignment2.printAll();


						Alignment minimap2Alignment = minimap2Align(seq2, seq1);
						cout << minimap2Alignment._readStart << " " << minimap2Alignment._readEnd << "    " << minimap2Alignment._contigStart << " " << minimap2Alignment._contigEnd << endl;
						//getchar();
						//Alignment al = _parent.getOverlap(alignment1._readIndex, alignment2._readIndex);
						//cout << al._readIndex << " " << al._contigIndex << endl;
						//cout << al._readStart << " " << al._readEnd << "    " << al._contigStart << " " << al._contigEnd << "     " << al._strand << endl;

					}
					*/

					//alignment1.print();
					//alignment2.print();

					//cout << "Expected: " << al._readEnd-al._readStart << " " << al._contigEnd-al._contigStart << endl;
					//cout << "Estimated: " << estimatedOverlapLength(alignment1, alignment2)  << endl;
					
					//if(alignment1._contigStart >  100000)				getchar();
					//bool lala1 = isOverlap(alignment1._readIndex, alignment2._readIndex);
					//bool lala2 = isOverlap2(alignment1, alignment2);
					//cout << lala1 << " " << lala2 << endl;

					//if(lala1 != lala2){
					//	Alignment2 al = getOverlap(alignment1._readIndex, alignment2._readIndex);
					//	cout << "Expected: " << al._readStart << " " << al._readEnd << " " << al._contigStart << " " << al._contigEnd << endl;
					//	getchar();
					//}

					//if(!isOverlap(alignment1._readIndex, alignment2._readIndex)) continue;

					//cout << "\tAdd edge: " << endl;
					//cout << "\tdone" << endl;
					int32_t overlapLength = _parent.getEstimatedOverlapLength(alignment1, alignment2);
					if(overlapLength < minOverlap) continue;

					successorAlignemnts.push_back({alignment2._readIndex, overlapLength});
					//graph->addEdge(alignment1._readIndex, alignment2._readIndex);
				}

				std::sort(successorAlignemnts.begin(), successorAlignemnts.end(), [](const SuccessorOverlapSorter& a, const SuccessorOverlapSorter& b){
					return a._estimatedOverlapLength > b._estimatedOverlapLength;
				});

				for(size_t i=0; i<successorAlignemnts.size() && i< 500; i++){
					graph->addEdge(alignment1._readIndex, successorAlignemnts[i]._readIndex);
					nbEdges += 1;
				}
			}

			

			//cout << "doneeee" << endl;

			return graph;
		}
		

		void transitiveReduction2(OverlapGraph* graph){

			int nbReduction = 0;

			const int FUZZ = 10;
			const int VACANT = 0;
			const int INPLAY = 1;
			const int ELIMINATED = 2;

			unordered_map<u_int32_t, int> mark;

			//for(int i=0; i<10; i++){

			for(const auto& it : graph->_readIndex_to_nodeIndex){


				const ReadType& readIndex = it.first;
				const u_int32_t& nodeIndex = it.second;

				mark[nodeIndex] = VACANT;

				//for(OverlapGraph::Edge& successor : graph->_successors[nodeIndex]){
				//	successor._reduce = false;
				//}

			}


			for(const auto& it : graph->_readIndex_to_nodeIndex){

				const ReadType& readIndex = it.first;
				const u_int32_t& nodeIndex = it.second;

				//cout << "-----" << endl;
				//cout << nodeIndex << " " << readIndex << endl;

				if(graph->_successors[nodeIndex].size() <= 1) continue;

				for(const OverlapGraph::Edge& successor : graph->_successors[nodeIndex]){
					if(successor._isRemoved) continue;
					mark[successor._nodeIndex] = INPLAY;
				}

				//int32_t longest = graph->_successors[nodeIndex][0]._overlapLength;
				//for(const OverlapGraph::Edge& successor : graph->_successors[nodeIndex]){
				//	cout << "\tSuccessor: " << successor._nodeIndex << " " << successor._overlapLength << endl;
				//}
				//cout << "\tLongest: " << longest << endl;
				//getchar();
				for(const OverlapGraph::Edge& successor : graph->_successors[nodeIndex]){
					if(successor._isRemoved) continue;
					if (mark[successor._nodeIndex] != INPLAY) continue;
					
					for(const OverlapGraph::Edge& successor2 : graph->_successors[successor._nodeIndex]){
						if(successor2._isRemoved) continue;
						//if(!overlapOnthe reference)
						if (mark[successor2._nodeIndex] != INPLAY) continue;
						//if(successor2._overlapLength + successor._overlapLength <= longest){
							mark[successor2._nodeIndex] = ELIMINATED;
						//}
					}
				}

				/*
				for(OverlapGraph::Edge& successor : graph->_successors[nodeIndex]){
					if(successor._isRemoved) continue;

					for(size_t i=0; i<graph->_successors[successor._nodeIndex].size(); i++){
						const OverlapGraph::Edge& successor2 = graph->_successors[successor._nodeIndex][i];
						if(successor2._isRemoved) continue;

						if(successor2._overlapLength < FUZZ || i == graph->_successors[successor._nodeIndex].size()-1){
							if (mark[successor2._nodeIndex] != INPLAY) continue;
							mark[successor2._nodeIndex] = ELIMINATED;
						}
					}
				}
				*/

				//cout << "-----" << endl;
				//for(OverlapGraph::Edge& successor : graph->_successors[nodeIndex]){
				//	cout << "\t" << successor._nodeIndex << " " << successor._overlapLength << endl;
				//}

				for(OverlapGraph::Edge& successor : graph->_successors[nodeIndex]){
					//if(successor._isRemoved) continue;

					if (mark[successor._nodeIndex] == ELIMINATED){
						graph->removeEdge(nodeIndex, successor._nodeIndex);
						//cout << "\tEliminating: " << successor._nodeIndex << " " << successor._overlapLength << endl;
						//successor._isRemoved = true;
						nbReduction += 1;
					}

					mark[successor._nodeIndex] = VACANT;


				}

				//getchar();
			}

			if(_print_debug) cout << "nb reduction: " << nbReduction << endl;
			//}

			//getchar();
		}

		void cleanRepeatOverlaps(OverlapGraph* graph){

			u_int64_t nbReperatedOverlaps = 0;

			for(const auto& it : graph->_readIndex_to_nodeIndex){

				const ReadType& readIndex = it.first;
				const u_int32_t& nodeIndex = it.second;

				if(graph->_successors[nodeIndex].size() <= 1) continue;

				const Alignment& currentAlignment = _parent._alignments_readsVsContigs[graph->_contigIndex][readIndex];
				const OverlapGraph::Edge& firstSuccessor = graph->_successors[nodeIndex][0];
				const ReadType& firstSuccessorReadIndex = graph->_nodeIndex_to_readIndex[firstSuccessor._nodeIndex];
				int32_t longestOverlapLength = _parent.getEstimatedOverlapLength(currentAlignment, _parent._alignments_readsVsContigs[graph->_contigIndex][firstSuccessorReadIndex]);
				
				//cout << "Longest overlap: " << longestOverlapLength << endl;
				float minOverlapLength = longestOverlapLength * 0.5;

				for(OverlapGraph::Edge& successor : graph->_successors[nodeIndex]){
					if(successor._isRemoved) continue;

					const ReadType& successorReadIndex = graph->_nodeIndex_to_readIndex[successor._nodeIndex];
					int32_t overlapLength = _parent.getEstimatedOverlapLength(currentAlignment, _parent._alignments_readsVsContigs[graph->_contigIndex][successorReadIndex]);
					//cout << successor._nodeName << " " << overlapLength << endl;
					if(overlapLength < minOverlapLength){
						graph->removeEdge(nodeIndex, successor._nodeIndex);
						nbReperatedOverlaps += 1;
					}
				}

				//getchar();
			}

			if(_print_debug) cout << "Nb repeated overlaps: " << nbReperatedOverlaps << endl;
		}

		/*
		void exportToGfa(OverlapGraph* graph, const string& outputFilename){

			ofstream outputFile(outputFilename);
			ofstream posFile(outputFilename + "_pos.csv");
			posFile << "Name,Pos" << endl;
			
			vector<ReadType> nodes;
			for(const auto& it : graph->_readIndex_to_nodeIndex){
				nodes.push_back(it.first);
			}

			int i=0;

			#pragma omp parallel num_threads(_parent._nbCores)
			{

				#pragma omp for
				for(size_t i=0; i<nodes.size(); i++){


					const ReadType& readIndex = nodes[i];
					const u_int32_t& nodeIndex = graph->_readIndex_to_nodeIndex[readIndex];

					const Alignment& currentAlignment = _parent._alignments_readsVsContigs[graph->_contigIndex][readIndex];

					char* dnaStr = _parent._readIndex_to_sequence[readIndex]->to_string();
					string readSeq = string(dnaStr, strlen(dnaStr));
					free(dnaStr);

					#pragma omp critical(exportToGfa)
					{
						cout << i << " " << graph->_readIndex_to_nodeIndex.size() << endl;
						i += 1;
						outputFile << "S" << "\t" << readIndex << "\t" << readSeq << "\t" << "LN:i:" << readSeq.size() << endl; //<< "\t" << "dp:i:" << _nodes[unitigIndex]->_abundance << endl;
					
						posFile << readIndex << "," << currentAlignment._contigStart << endl;
					}

					for(const OverlapGraph::Edge& successor : graph->_successors[nodeIndex]){
						if(successor._isRemoved) continue;

						const ReadType& successorReadIndex = graph->_nodeIndex_to_readIndex[successor._nodeIndex];
						const Alignment& successorAlignment = _parent._alignments_readsVsContigs[graph->_contigIndex][successorReadIndex];

						if(_parent.getEstimatedOverlapLength(currentAlignment, successorAlignment) < 1000) continue;

						AlignmentBounds chaingingALignment = _parent.computeReadOverlapBounds(currentAlignment, successorAlignment, _tbuf);


						if(!isValidAlignment(chaingingALignment, _parent.getEstimatedOverlapLength(currentAlignment, successorAlignment))) continue;

						u_int32_t overlap = chaingingALignment._queryEnd - chaingingALignment._queryStart;

						#pragma omp critical(exportToGfa)
						{	
							//cout << overlap << " " << _parent.getEstimatedOverlapLength(currentAlignment, successorAlignment) << endl;
							outputFile << "L" << "\t" << currentAlignment._readIndex << "\t" << "+" << "\t" << successorAlignment._readIndex << "\t" << "+" << "\t" << overlap << "M" << endl;

							//getchar();
						}

					}
				}

			}

			outputFile.close();
			posFile.close();




		}
		*/

		void computeOverlapPaths(OverlapGraph* graph, const vector<Alignment>& alignments, vector<OverlapGraphPath>& paths){


			Alignment dummyAlignment;
			dummyAlignment._readIndex = -1;

			if(_print_debug) cout << "a" << endl;

			paths.clear();
			
			while(true){

				Alignment sourceAlignment;
				bool isAlignment = getBestAlignment(dummyAlignment, alignments, paths, sourceAlignment, graph);

				if(!isAlignment) break;

				u_int32_t sourceNodeName = sourceAlignment._readIndex;

				vector<u_int32_t> mostSupportedPath_successors;
				computeMostSupportedPath_direction(graph, sourceNodeName, true, paths, mostSupportedPath_successors);

				vector<u_int32_t> mostSupportedPath_predecessors;
				computeMostSupportedPath_direction(graph, sourceNodeName, false, paths, mostSupportedPath_predecessors);
				
				vector<u_int32_t> path;

				for(size_t i=0; i<mostSupportedPath_predecessors.size()-1; i++){ //-1: Discard source node here
					path.push_back(mostSupportedPath_predecessors[i]);
				}

				for(size_t i=0; i<mostSupportedPath_successors.size(); i++){
					path.push_back(mostSupportedPath_successors[i]);
				}
				
				const Alignment& startAlignment = _parent._alignments_readsVsContigs[graph->_contigIndex][path[0]];
				const Alignment& endAlignment = _parent._alignments_readsVsContigs[graph->_contigIndex][path[path.size()-1]];

				u_int32_t pathIndex = paths.size();
				paths.push_back({graph->_contigIndex, pathIndex, startAlignment._contigStart, endAlignment._contigEnd, path});
				
				if(_print_debug) cout << "Existing paths:" << endl;

				for(const OverlapGraphPath& path : paths){
					if(_print_debug) cout << "\t" << path._referenceStart << " " << path._referenceEnd << endl;
					/*
					Alignment currentAlignment = _parent._alignments_readsVsContigs[graph->_contigIndex][path._path[path._path.size()-1]];
					currentAlignment.printAll();

					for(const OverlapGraph::Edge& successor : graph->_successors[graph->_nodeName_to_nodeIndex[path._path[path._path.size()-1]]]){
						if(successor._isRemoved) continue;

						Alignment al = _parent._alignments_readsVsContigs[graph->_contigIndex][successor._nodeName];
						al.printAll();

						AlignmentBounds alignment = _parent.computeReadOverlap(currentAlignment, al, _tbuf);
						alignment.printAll();
						//alignmentBounds.printAll();
						cout << isValidAlignment(alignment, _parent.getEstimatedOverlapLength(currentAlignment, al)) << endl;

					}
					*/	
				}
				
				//getchar();
			}

			if(_print_debug) cout << "b" << endl;

		}
			

		void computeMostSupportedPath_direction(OverlapGraph* graph, u_int32_t sourceNodeName, bool useSuccessors, const vector<OverlapGraphPath>& existingPaths, vector<u_int32_t>& longestPath){

			//cout << "aa" << endl;

			int32_t nbAttempts = 1000;

			u_int64_t longestPathLength = 0;
			longestPath.clear();
			unordered_set<u_int32_t> removedNodes;

			while(true){
			
				vector<u_int32_t> mostSupportedPath;
				computeMostSupportedPath(graph, sourceNodeName, removedNodes, useSuccessors, existingPaths, mostSupportedPath);

				if(!useSuccessors){
					std::reverse(mostSupportedPath.begin(), mostSupportedPath.end());
				}

				//if(useSuccessors){
				//	cout << nbAttempts << endl;
				//	cout << getPathLength(graph->_contigIndex, longestPath) << endl;
				//	_parent._alignments_readsVsContigs[graph->_contigIndex][mostSupportedPath[mostSupportedPath.size()-1]].printAll();
					
				//}
				if(useSuccessors){
					removedNodes.insert(mostSupportedPath[mostSupportedPath.size()-1]);
				}
				else{
					removedNodes.insert(mostSupportedPath[0]);
				}

				if (longestPath.size() == 0){
					longestPath = mostSupportedPath;
					longestPathLength = getPathLength(graph->_contigIndex, longestPath);
				}
				else{
					int64_t pathLength = getPathLength(graph->_contigIndex, mostSupportedPath);
					
					if (pathLength > longestPathLength){
						longestPath = mostSupportedPath;
						longestPathLength = pathLength;
					}
				}


				if (mostSupportedPath.size() <= 1) break;

				nbAttempts -= 1;
				if (nbAttempts <= 0) break;
				
				//cout << removedNodes.size() << " " << nbAttempts << " " << longestPathLength << endl;
			}
			
			//cout << "bb" << endl;
		}

		/*
		void printPathSupport(OverlapGraph* graph, const OverlapGraphPath& path){

			for(int64_t i=0; i <= path._path.size(); i++){

				const ReadType& readIndex = path._path[i];
				int nbSuccessors = graph->_successors[graph->_readIndex_to_nodeIndex[readIndex]].size();

				cout << i << " " << nbSuccessors << endl;

			}

			for(int64_t i=path._path.size()-1; i>=0; i--){
				
				const ReadType& readIndex = path._path[i];

				int nbSuccessors = graph->_predecessors[graph->_readIndex_to_nodeIndex[readIndex]].size();

				cout << i << " " << nbSuccessors << endl;
			}

		}
		*/

		int64_t getPathLength(const u_int32_t& contigIndex, const vector<u_int32_t>& path){
			const Alignment& startAlignment = _parent._alignments_readsVsContigs[contigIndex][path[0]];
			const Alignment& endAlignment = _parent._alignments_readsVsContigs[contigIndex][path[path.size()-1]];

			return endAlignment._contigEnd - startAlignment._contigStart;
		}

		void computeMostSupportedPath(OverlapGraph* graph, u_int32_t currentNodeName, unordered_set<u_int32_t>& removedNodes, bool useSuccessors, const vector<OverlapGraphPath>& existingPaths, vector<u_int32_t>& path){

			//cout << "aaa" << endl;
			Alignment bestSuccessorAlignment;
			unordered_set<u_int32_t> isVisited;
			path.clear();

			while(true){


				u_int32_t maxNodeName = -1;

				path.push_back(currentNodeName);

				//cout << currentNodeName << " " << path.size() << endl;

				if(isVisited.find(currentNodeName) != isVisited.end()){
					cout << "Cycle in the graph!" << endl;
					break;
				}
				

				isVisited.insert(currentNodeName);

				const Alignment& currentAlignment = _parent._alignments_readsVsContigs[graph->_contigIndex][currentNodeName];


				//cout << "ccc: " << currentNodeName << endl;

				bool isBestAlignment = false;

				while(true){

					vector<Alignment> successorAlignments;
					/*
					const vector<OverlapGraph::Edge>& successors = getSuccessors(graph, currentNodeName, useSuccessors);
					if(successors.size() == 0) break;

					//int32_t maxEstimatedOverlapLength = getEstimatedOverlapLength(currentAlignment, _alignments_readsVsContigs[graph->_contigIndex][successors[0]._nodeName]);
					//int32_t minOverlapLength = maxEstimatedOverlapLength * 0.8;
					//cout << "\tSucc: " << currentNodeName << " " << (graph->_nodeName_to_nodeIndex.find(currentNodeName) != graph->_nodeName_to_nodeIndex.end()) << endl;

					for(const OverlapGraph::Edge& successor : successors){
						if(successor._isRemoved) continue;
						if(removedNodes.find(successor._nodeName) != removedNodes.end()) continue;
						if(isNotOverlap(currentNodeName, successor._nodeName)) continue;
						
						//int32_t estimatedOverlapLength = getEstimatedOverlapLength(currentAlignment, _alignments_readsVsContigs[graph->_contigIndex][successor._nodeName]);
						//if(estimatedOverlapLength < minOverlapLength) continue;

						successorAlignments.push_back(_alignments_readsVsContigs[graph->_contigIndex][successor._nodeName]);
					}
					*/
					
					if (useSuccessors){
						/*
						cout << "----" << endl;

						for(const OverlapGraph::Edge& successor : graph->_successors[graph->_nodeName_to_nodeIndex[currentNodeName]]){
							if(successor._isRemoved) continue;
							Alignment al = _parent._alignments_readsVsContigs[graph->_contigIndex][successor._nodeName];
							//cout << successor._nodeName << "\t" << al._identity << "\t" << al._contigEnd - al._contigStart << "\t" << getEstimatedOverlapLength(currentAlignment, al) << endl;
							cout << "\t" << successor._nodeName << "\t" << al._identity << "\t" << al._contigStart << "\t" << al._contigEnd  << "\t" << getEstimatedOverlapLength(currentAlignment, al) << "\t" << isPathFrom(graph, graph->_nodeName_to_nodeIndex[currentNodeName], successor._nodeIndex, currentAlignment) << endl;

							for(const OverlapGraph::Edge& successor2 : graph->_successors[successor._nodeIndex]){
								if(successor2._isRemoved) continue;
								Alignment al2 = _parent._alignments_readsVsContigs[graph->_contigIndex][successor2._nodeName];
								//cout << successor._nodeName << "\t" << al._identity << "\t" << al._contigEnd - al._contigStart << "\t" << getEstimatedOverlapLength(currentAlignment, al) << endl;
								cout << "\t\t" << successor2._nodeName << "\t" << al2._identity << "\t" << al2._contigStart << "\t" << al2._contigEnd  << endl;//<< "\t" << getEstimatedOverlapLength(currentAlignment, al) << endl;

							}
						}
						
						getchar();
						*/
					
						//const vector<OverlapGraph::Edge>& successors = graph->_successors[graph->_nodeName_to_nodeIndex[currentNodeName]];
						//if(successors.size() == 0) break;

						//int32_t maxEstimatedOverlapLength = getEstimatedOverlapLength(currentAlignment, _alignments_readsVsContigs[graph->_contigIndex][successors[0]._nodeName]);
						//int32_t minOverlapLength = maxEstimatedOverlapLength * 0.8;
						//cout << "\tSucc: " << currentNodeName << " " << (graph->_nodeName_to_nodeIndex.find(currentNodeName) != graph->_nodeName_to_nodeIndex.end()) << endl;

						for(const OverlapGraph::Edge& successor : graph->_successors[graph->_readIndex_to_nodeIndex[currentNodeName]]){
							if(successor._isRemoved) continue;
							const ReadType& successorReadIndex = graph->_nodeIndex_to_readIndex[successor._nodeIndex];
							if(removedNodes.find(successorReadIndex) != removedNodes.end()) continue;
							if(isNotOverlap(currentNodeName, successorReadIndex)) continue;
							
							//int32_t estimatedOverlapLength = getEstimatedOverlapLength(currentAlignment, _alignments_readsVsContigs[graph->_contigIndex][successor._nodeName]);
							//if(estimatedOverlapLength < minOverlapLength) continue;

							successorAlignments.push_back(_parent._alignments_readsVsContigs[graph->_contigIndex][successorReadIndex]);
						}
					}
					else{
						for(const OverlapGraph::Edge& predecessor : graph->_predecessors[graph->_readIndex_to_nodeIndex[currentNodeName]]){
							if(predecessor._isRemoved) continue;
							const ReadType& predecessorReadIndex = graph->_nodeIndex_to_readIndex[predecessor._nodeIndex];
							if(removedNodes.find(predecessorReadIndex) != removedNodes.end()) continue;
							if(isNotOverlap(predecessorReadIndex, currentNodeName)) continue;
							successorAlignments.push_back(_parent._alignments_readsVsContigs[graph->_contigIndex][predecessorReadIndex]);
						}
					}
					
					

					

					//cout << "2" << endl;
					isBestAlignment = getBestOverlap(currentAlignment, successorAlignments, existingPaths, bestSuccessorAlignment);
					


					//cout << "3" << endl;
					if(!isBestAlignment) break;

					//cout << "4" << endl;
					if (useSuccessors){

						//cout << currentAlignment._readIndex << " " << _alignments_readsVsContigs[graph->_contigIndex][bestSuccessorAlignment._readIndex]._readIndex << endl;
						if(isOverlap2_cached(currentAlignment, _parent._alignments_readsVsContigs[graph->_contigIndex][bestSuccessorAlignment._readIndex])){
							break;
						}
						//cout << "lol: " << isOverlap2_cached(currentAlignment, _alignments_readsVsContigs[graph->_contigIndex][bestSuccessorAlignment._readIndex]) << endl;
						//cout << "lol: " << isNotOverlap(currentNodeName, bestSuccessorAlignment._readIndex) << endl;
					}
					else{
						if(isOverlap2_cached(_parent._alignments_readsVsContigs[graph->_contigIndex][bestSuccessorAlignment._readIndex], currentAlignment)){
							break;
						}
						//cout << "lol: " << isNotOverlap(bestSuccessorAlignment._readIndex, currentNodeName) << endl;
					}



					//cout << "5" << endl;
				}


				//cout << "ddd: " << currentNodeName << endl;


				if(isBestAlignment){
					maxNodeName = bestSuccessorAlignment._readIndex;
				}
				else{
					break;
				}

				//cout << "fff: " << currentNodeName << endl;

				currentNodeName = maxNodeName;
				
				//cout << "ggg: " << currentNodeName << endl;
			}

			//cout << "bbb" << endl;
		}

		bool getBestAlignment(const Alignment& prevAlignment, const vector<Alignment>& alignments, const vector<OverlapGraphPath>& existingPaths, Alignment& resultAlignment, OverlapGraph* graph){

			///static vector<u_int32_t> overlapLengths = {10000, 9000, 8000, 7000, 6000, 5000, 4000, 3000, 2000, 1000, minOverlap};
			static vector<u_int32_t> lengths1 = {100000, 50000, 25000, 10000, 5000};
			static vector<float> identities1 = {0.99, 0.98};

			if(getBestAlignmentSub(prevAlignment, alignments, existingPaths, resultAlignment, graph, lengths1, identities1)) return true;			

			static vector<u_int32_t> lengths2 = {100000, 50000, 25000, 10000, 5000, 0};
			static vector<float> identities2 = {0.99, 0.98, 0.97, 0.96, 0.95, 0.94, 0.93, 0.92, 0.91, 0.9, 0.8, 0.7, 0.6, 0.5};

			if(getBestAlignmentSub(prevAlignment, alignments, existingPaths, resultAlignment, graph, lengths2, identities2)) return true;	

			return false;
		}
		
		bool getBestAlignmentSub(const Alignment& prevAlignment, const vector<Alignment>& alignments, const vector<OverlapGraphPath>& existingPaths, Alignment& resultAlignment, OverlapGraph* graph, const vector<u_int32_t>& lengths, const vector<float>& identities){

			for(float minIdentity : identities){
				for(u_int32_t minLength : lengths){

					float highestIdentity = 0;

					for(const Alignment& alignment : alignments){

						if(alignment._identity < minIdentity) continue;
						if(alignment._contigEnd - alignment._contigStart < minLength) continue;

						
						if(graph != nullptr){
							if(graph->_readIndex_to_nodeIndex.find(alignment._readIndex) == graph->_readIndex_to_nodeIndex.end()) continue;
						}

						if(alignmentIsContainsInExistingPath(alignment, existingPaths)) continue;

						if(alignment._identity > highestIdentity){
							highestIdentity = alignment._identity;
							resultAlignment = alignment;
						}
					}

					if(highestIdentity > 0){
						return true;
					}
						
				}
			}
						

			return false;
		}

		bool getBestOverlap(const Alignment& prevAlignment, const vector<Alignment>& alignments, const vector<OverlapGraphPath>& existingPaths, Alignment& resultAlignment){
			/*
			static vector<u_int32_t> overlapLengths = {10000, 9000, 8000, 7000, 6000, 5000, 4000, 3000, 2000, 1000, minOverlap};
			//static vector<u_int32_t> lengths = {100000, 50000, 25000, 10000, 5000, 0};
			static vector<float> identities = {0.97, 0.96, 0.95, 0.94, 0.93, 0.92, 0.91, 0.9, 0.8, 0.7, 0.6, 0.5};
			
			for(float minIdentity : identities){
				for(u_int32_t minLength : overlapLengths){

					float highestIdentity = 0;

					for(const Alignment& alignment : alignments){

						if(alignment._identity < minIdentity) continue;

						int32_t overlapLength = _parent.getEstimatedOverlapLength(prevAlignment, alignment);
						if(overlapLength < minLength) continue;

						if(alignmentIsContainsInExistingPath(alignment, existingPaths)) continue;


						if(alignment._identity > highestIdentity){
							highestIdentity = alignment._identity;
							resultAlignment = alignment;
						}
					}

					if(highestIdentity > 0){
						return true;
					}
						
				}
			}
						

			return false;
			*/
			///static vector<u_int32_t> overlapLengths = {10000, 9000, 8000, 7000, 6000, 5000, 4000, 3000, 2000, 1000, minOverlap};
			static vector<u_int32_t> lengths1 = {10000, 9000, 8000, 7000, 6000, 5000};
			static vector<float> identities1 = {0.99, 0.98};

			if(getBestOverlapSub(prevAlignment, alignments, existingPaths, resultAlignment, lengths1, identities1)) return true;			

			static vector<u_int32_t> lengths2 = {10000, 9000, 8000, 7000, 6000, 5000, 4000, 3000, 2000, 1000, minOverlap};
			static vector<float> identities2 = {0.99, 0.98, 0.97, 0.96, 0.95, 0.94, 0.93, 0.92, 0.91, 0.9, 0.8, 0.7, 0.6, 0.5};

			if(getBestOverlapSub(prevAlignment, alignments, existingPaths, resultAlignment, lengths2, identities2)) return true;	

			return false;

		}

		bool getBestOverlapSub(const Alignment& prevAlignment, const vector<Alignment>& alignments, const vector<OverlapGraphPath>& existingPaths, Alignment& resultAlignment, const vector<u_int32_t>& overlapLengths, const vector<float>& identities){

			//static vector<u_int32_t> overlapLengths = {10000, 9000, 8000, 7000, 6000, 5000, 4000, 3000, 2000, 1000, minOverlap};
			//static vector<u_int32_t> lengths = {100000, 50000, 25000, 10000, 5000, 0};
			//static vector<float> identities = {0.97, 0.96, 0.95, 0.94, 0.93, 0.92, 0.91, 0.9, 0.8, 0.7, 0.6, 0.5};
			
			for(float minIdentity : identities){
				for(u_int32_t minLength : overlapLengths){

					float highestIdentity = 0;

					for(const Alignment& alignment : alignments){

						if(alignment._identity < minIdentity) continue;

						int32_t overlapLength = _parent.getEstimatedOverlapLength(prevAlignment, alignment);
						if(overlapLength < minLength) continue;

						if(alignmentIsContainsInExistingPath(alignment, existingPaths)) continue;


						if(alignment._identity > highestIdentity){
							highestIdentity = alignment._identity;
							resultAlignment = alignment;
						}
					}

					if(highestIdentity > 0){
						return true;
					}
						
				}
			}
						

			return false;
		}

		bool isOverlap(const ReadType& prevReadIndex, const ReadType& nextReadIndex){


			std::pair<ReadType, ReadType> key;

			if (prevReadIndex < nextReadIndex){
				key.first = prevReadIndex;
				key.second = nextReadIndex;
			}
			else{
				key.first = nextReadIndex;
				key.second = prevReadIndex;
			}

			return _parent._alignments_readsVsReads.find(key) != _parent._alignments_readsVsReads.end();
			
		}
		
		bool isNotOverlap(const ReadType& prevReadIndex, const ReadType& nextReadIndex){
			
			std::pair<ReadType, ReadType> key;

			if (prevReadIndex < nextReadIndex){
				key.first = prevReadIndex;
				key.second = nextReadIndex;
			}
			else{
				key.first = nextReadIndex;
				key.second = prevReadIndex;
			}

			if(_alignments_readsVsReads_cached.find(key) == _alignments_readsVsReads_cached.end()) return false;

			return !_alignments_readsVsReads_cached[key];
		}

		bool isOverlap2_cached(const Alignment& prevAlignment, const Alignment& nextAlignment){

			std::pair<ReadType, ReadType> key;

			if (prevAlignment._readIndex < nextAlignment._readIndex){
				key.first = prevAlignment._readIndex;
				key.second = nextAlignment._readIndex;
			}
			else{
				key.first = nextAlignment._readIndex;
				key.second = prevAlignment._readIndex;
			}

			if(_alignments_readsVsReads_cached.find(key) != _alignments_readsVsReads_cached.end()){
				return _alignments_readsVsReads_cached[key];
			}

			_alignments_readsVsReads_cached[key] = isOverlap2(prevAlignment, nextAlignment);

			return _alignments_readsVsReads_cached[key];
		}

		bool isOverlap2(const Alignment& prevAlignment, const Alignment& nextAlignment){

			AlignmentBounds chaingingALignment = _parent.computeReadOverlapBounds(prevAlignment, nextAlignment, _tbuf);

			

			//chaingingALignment.printAll();
			
			//if(prevAlignment._readIndex == 5735 && nextAlignment._readIndex == 37504){
			//	isValidAlignment(chaingingALignment, _parent.getEstimatedOverlapLength(prevAlignment, nextAlignment));
			//	cout << "allo??" << endl;
			//	getchar();
			//}
			//cout << isValidAlignment(chaingingALignment, _parent.getEstimatedOverlapLength(prevAlignment, nextAlignment)) << endl;

			//if(chaingingALignment._queryStart != -1){
			//	if(isValidAlignment(chaingingALignment, _parent.getEstimatedOverlapLength(prevAlignment, nextAlignment))){
			//		return true;
			//	}
			//}

			//cout << prevAlignment._readIndex << " " << nextAlignment._readIndex << endl;
			//getchar();
			return _parent.isValidAlignment(chaingingALignment, _parent.getEstimatedOverlapLength(prevAlignment, nextAlignment));
			
			//cout << "aaaaaa" << endl;

			//cout << (_readIndex_to_sequence.find(prevAlignment._readIndex) != _readIndex_to_sequence.end()) << endl;
			//cout << (_readIndex_to_sequence.find(nextAlignment._readIndex) != _readIndex_to_sequence.end()) << endl;


			/*
			char* dnaStr1 = _readIndex_to_sequence[prevAlignment._readIndex]->to_string();
			string readSeq1 = string(dnaStr1, strlen(dnaStr1));
			free(dnaStr1);

			char* dnaStr2 = _readIndex_to_sequence[nextAlignment._readIndex]->to_string();
			string readSeq2 = string(dnaStr2, strlen(dnaStr2));
			free(dnaStr2);

			if(prevAlignment._strand){
				Utils::toReverseComplement(readSeq1);
			}
			
			if(nextAlignment._strand){
				Utils::toReverseComplement(readSeq2);
			}
			*/

			/*
			//cout << "bbbbbb" << endl;
			//cout << "Strand: " << prevAlignment._strand << " " << nextAlignment._strand << endl;
			//cout << isReversed1 << " " << isReversed2 << endl;
			//AlignmentBounds alignmentBounds = _parent.computeReadOverlap(prevAlignment, nextAlignment, true, chaingingALignment, _aligner);
			AlignmentBounds alignmentBounds = _parent.computeReadOverlap(prevAlignment, nextAlignment, _tbuf);
			//cout << alignmentBounds._queryStart << endl;

			if(alignmentBounds._queryStart == -1){
				//getchar();
				return false;
			}



			//cout << "Expected: " << expectedOverlapLength << endl;
			//cout << "Actual:   " << alignmentBounds._queryEnd - alignmentBounds._queryStart << endl;
			//cout << "Error:    " << error << endl;

			return isValidAlignment(alignmentBounds, _parent.getEstimatedOverlapLength(prevAlignment, nextAlignment));
			
			//AlignmentBounds alignmentBounds = computeReadOverlapBounds(prevAlignment._readIndex, nextAlignment._readIndex);
			
			//cout << "Bounds: " << alignmentBounds._queryStart << " " << alignmentBounds._queryEnd << " " << alignmentBounds._referenceStart << " " << alignmentBounds._referenceEnd << endl;
			

			//cout << "cccccc" << endl;
			*/

		}

		bool alignmentIsContainsInExistingPath(const Alignment& alignment, const vector<OverlapGraphPath>& existingPaths){

			for(const OverlapGraphPath& path : existingPaths){
				if (alignment._contigStart >= path._referenceStart && alignment._contigStart <= path._referenceEnd) return true;
				if (alignment._contigEnd >= path._referenceStart && alignment._contigEnd <= path._referenceEnd) return true;
			}

			return false;
		}
		



		bool overlapOnTheReference(const Alignment& prevAlignment, const Alignment& nextAlignment){

			
			int32_t prevLeftOverhang = 0;
			int32_t prevRightOverhang = 0;

			if(prevAlignment._strand){
				prevLeftOverhang = prevAlignment._readLength - prevAlignment._readEnd;
				prevRightOverhang = prevAlignment._readStart;
			}
			else{
				prevLeftOverhang = prevAlignment._readStart;
				prevRightOverhang = prevAlignment._readLength - prevAlignment._readEnd;
			}

			int32_t nextLeftOverhang = 0;
			int32_t nextRightOverhang = 0;

			if(nextAlignment._strand){
				nextLeftOverhang = nextAlignment._readLength - nextAlignment._readEnd;
				nextRightOverhang = nextAlignment._readStart;
			}
			else{
				nextLeftOverhang = nextAlignment._readStart;
				nextRightOverhang = nextAlignment._readLength - nextAlignment._readEnd;
			}

			//prevLeftOverhang = 0;
			//prevRightOverhang = 0;
			//nextLeftOverhang = 0;
			//nextRightOverhang = 0;

			int32_t prevStart = ((int32_t) prevAlignment._contigStart) - prevLeftOverhang;
			int32_t prevEnd = ((int32_t) prevAlignment._contigEnd) + prevRightOverhang;

			int32_t nextStart = ((int32_t) nextAlignment._contigStart) - nextLeftOverhang;
			int32_t nextEnd = ((int32_t) nextAlignment._contigEnd) + nextRightOverhang;


			u_int32_t alignmentOffset = minOverlap;

			//u_int32_t refStart = prevAlignment._contigStart;
			//u_int32_t refEnd = prevAlignment._contigEnd;

			/*
			cout << "----" << endl;
			cout << prevAlignment._strand << " " << nextAlignment._strand << endl;
			cout << prevAlignment._contigStart << "\t" << prevAlignment._contigEnd << endl;
			cout << nextAlignment._contigStart << "\t" << nextAlignment._contigEnd << endl;
			cout << endl;
			cout << prevLeftOverhang << "\t" << prevRightOverhang << endl;
			cout << prevStart << "\t" << prevEnd << endl;
			cout << nextStart << "\t" << nextEnd << endl;
			getchar();
			*/

			if (nextStart > prevStart+alignmentOffset && nextStart < prevEnd-alignmentOffset && nextEnd > prevEnd+alignmentOffset) return true;

			return false;
		}

		



	};




	bool isValidAlignment(const AlignmentBounds& alignmentBounds, const int32_t& expectedOverlapLength){

		if(alignmentBounds._queryStart == -1) return false;
		if(alignmentBounds._isReversed) return false; 

		int32_t actualOverlapLength = alignmentBounds._queryEnd - alignmentBounds._queryStart;
		int32_t error = abs(expectedOverlapLength - actualOverlapLength);
		
		//cout << error << endl;

		//if(error > 2000) return false;


		int64_t queryLength = alignmentBounds._queryLength;
		int64_t queryStart = alignmentBounds._queryStart;
		int64_t queryEnd = alignmentBounds._queryEnd;
		
		int64_t targetLength = alignmentBounds._referenceLength;
		int64_t targetStart = alignmentBounds._referenceStart;
		int64_t targetEnd = alignmentBounds._referenceEnd;
		bool isReversed = alignmentBounds._isReversed;
		//int64_t nbMatches = stoull((*_fields)[9]);
		//int64_t alignLength = stoull((*_fields)[10]);

		
		if(targetStart < queryStart) return false;
		//double length = max((double)(contigEnd - contigStart), (double)(readEnd - readStart));
		//double error = 1 - min((double)(contigEnd - contigStart), (double)(readEnd - readStart)) / length;
		//_logFile << error << " " << errorThreshold << endl;
		
		//if(error > errorThreshold) return;

		//float identity =  ((float)nbMatches) / alignLength;

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

		//cout << tl3 << " " << tl5 << " " << ext3 << " " << ext5 << " " << queryLength << " " << targetLength << endl;
		
		//cout << queryStart << " " << queryEnd << " " << ext5 << " " << ext3 << endl;
		//getchar();

		if (ext5 > maxHang || ext3 > maxHang || queryEnd - queryStart < (queryEnd - queryStart + ext5 + ext3) * intFrac) return false;


		if (queryStart <= tl5 && queryLength - queryEnd <= tl3) return false; //query contained
		else if (queryStart >= tl5 && queryLength - queryEnd >= tl3) return false; //target contained

		if (queryEnd - queryStart + ext5 + ext3 < minOverlap) return false; //Short overlap
		if (targetEnd - targetStart + ext5 + ext3 < minOverlap) return false; //Short overlap


		return true;

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

		size_t nbWindows = ceil(((double)read._seq.size()-_windowPositionOffset) / (double)_windowLength);
		
		if(_windowPositionOffset > 0){
			if(read._seq.size() > _windowPositionOffset) nbWindows += 1;
		}

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


		_contigIndex_to_nbPaths.clear();
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


	void polishPartition(u_int32_t partition, const string& contigFilename, gzFile& outputContigFile){
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
		performCorrection(outputContigFile);


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

			//cout << readIndex << " " << (_alignments.find(readIndex) != _alignments.end()) << endl;
			//if(_contigPolisher._currentPartition == 0) _logFile << readIndex << " " << (_alignments.find(readIndex) != _alignments.end()) << endl;
			
			//if(readIndex % 100000 == 0) _logFile << "\t" << readIndex << endl;

			if(_alignments.find(readIndex) == _alignments.end()) return;

			//const vector<Alignment>& als = _alignments[readIndex];
			const Alignment& al = _alignments[readIndex];
			//for(const Alignment& al : _alignments[readIndex]){
			u_int32_t contigIndex = al._contigIndex;

			//cout << "a " << readIndex << " " << contigIndex << " " << _contigSequences.size() << " " << (_contigSequences.find(contigIndex) != _contigSequences.end()) << endl;

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
			
			find_breaking_points_from_cigar(_windowLength, al, readSeq.size(), cigar, readSeq, qualSeq);
			free(cigar);

			//getchar();

			//}

		}

		void find_breaking_points_from_cigar(uint32_t window_length, const Alignment& al, u_int64_t readLength, char* cigar_, const string& readSequence, const string& qualSequence)
		{
			
			//cout << "Alignment: " << al._contigStart << " " << al._contigEnd << endl;
			vector<std::pair<int32_t, int32_t>> breaking_points_;
			int32_t t_begin_ = al._contigStart;
			int32_t t_end_ = al._contigEnd;
			int32_t q_begin_ = al._readStart;
			int32_t q_end_ = al._readEnd;
			bool strand_ = al._strand;
			int32_t q_length_ = readLength;

			// find breaking points from cigar
			std::vector<int32_t> window_ends;

			for (int32_t i = -_contigPolisher._windowPositionOffset; i < t_end_; i += window_length) {
				
				if (i > t_begin_) {
					window_ends.emplace_back(i - 1);
					//cout << "\tWindow end: " << (i-1) << endl;
				}
			}
			window_ends.emplace_back(t_end_ - 1);
			//cout << "\tWindow end: " << (t_end_ - 1) << endl;

			//getchar();
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

				//cout << "---" << endl;
				//cout << "Break point: " << j << " " << breaking_points_[j].first << " " << breaking_points_[j + 1].first << endl;
				//cout << "Break point: " << j << " " << breaking_points_[j].second << " " << breaking_points_[j + 1].second << endl;

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
				uint64_t window_id = (breaking_points_[j].first+_contigPolisher._windowPositionOffset) / _windowLength;
				int32_t window_start = window_id * _windowLength - _contigPolisher._windowPositionOffset;
				window_start = max(0, window_start);

				//cout << "Window id: " << window_id << endl;
				//cout << "Window start: " << window_start << endl;
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

                int32_t posStart = (breaking_points_[j].first) - window_start;
				posStart = max(0, posStart);
                int32_t posEnd =  (breaking_points_[j + 1].first) - window_start - 1;

				//cout << (breaking_points_[j].first) << " " << window_start << " " << posStart << endl;
				//cout << (breaking_points_[j + 1].first) << " " << window_start - 1 << " " << posEnd << endl;

				//getchar();

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

			//getchar();
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
				
				//cout << windowIndex << " " << _contigWindowSequences[al._contigIndex].size() << " " << posStart << " " << posEnd << " " << windowSequence << endl;

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


	void performCorrection(gzFile& outputContigFile){
		
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

				//cout << contigIndexLocal << " " << windowIndexLocal << " " << sequences.size() << endl;
				//u_int64_t nbCorrectedWindows = 0;
				//vector<DnaBitset2*> correctedWindows(windows.size());
			
				std::unique_ptr<spoa::AlignmentEngine> alignmentEngine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -4);
					
				
				int64_t wStart = (int32_t)windowIndexLocal * (int32_t) _windowLength - _windowPositionOffset;
				int64_t wEnd = min((int64_t)_contigSequences[contigIndexLocal].size(), (int64_t)(wStart+_windowLength)) ;
				wStart = max((int64_t)0, wStart);
				//cout << "a" << " " << _contigSequences[contigIndexLocal].size() << " " << wStart << " " << wEnd-wStart << endl;
				string contigOriginalSequence = _contigSequences[contigIndexLocal].substr(wStart, wEnd-wStart);
				//cout << "b" << endl;
				//cout << contigOriginalSequence << endl;

				if(sequences.size() < 2){

					//cout << "c" << endl;
					for(size_t i=0; i<sequences.size(); i++){ 
						delete sequences[i]._sequence;
					}

					//cout << "d" << endl;
					addCorrectedWindow(false, new DnaBitset2(contigOriginalSequence), contigIndexLocal, windowIndexLocal, outputContigFile);	
					//correctedWindows[w] = new DnaBitset2(contigOriginalSequence);
					//_logFile << "No sequences for window" << endl;
					continue;
				}
				

				//cout << "e" << endl;

				std::sort(sequences.begin(), sequences.end(), [&](const Window& lhs, const Window& rhs) {
					return lhs._posStart < rhs._posStart; });

				//_logFile << sequences.size() << endl;
				//_logFile << "1" << endl;
				spoa::Graph graph{};

				string backboneQuality = "";
				for(size_t i=0; i<contigOriginalSequence.size(); i++){
					backboneQuality += '!';
				}

				//cout << "f" << endl;
				graph.AddAlignment(
					spoa::Alignment(),
					contigOriginalSequence.c_str(), contigOriginalSequence.size(),
					backboneQuality.c_str(), backboneQuality.size()
				);

				u_int32_t offset = 0.01 * contigOriginalSequence.size();


				//cout << "g" << endl;
				for(size_t i=0; i<sequences.size(); i++){ 

					//size_t i = order[ii];
					const Window& window = sequences[i];
					//const DnaBitset2* dna = variant._sequence; //sequenceCopies[s._sequenceIndex];
					char* dnaStr = window._sequence->to_string();

					//cout << i << " " << string(dnaStr, strlen(dnaStr)) << " " << window._posStart << " " << window._posEnd << endl;
						
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
					
					//cout << "aa" << endl;
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
					
					//cout << "bb" << endl;
					free(dnaStr);
					delete window._sequence;
					//cout << "cc" << endl;
				}

				//cout << "h" << endl;

				std::vector<uint32_t> coverages;
				string correctedSequence = graph.GenerateConsensus(&coverages);

				uint32_t average_coverage = (sequences.size()) / 2;

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
					correctedSequence = correctedSequence.substr(begin, end - begin + 1);
				}

				//for(char letter : correctedSequence){
				//	checksum += letter;
				//}

				addCorrectedWindow(true, new DnaBitset2(correctedSequence), contigIndexLocal, windowIndexLocal, outputContigFile);

				//correctedWindows[w] = new DnaBitset2(correctedSequence);

			}
		}
		
		/*
		for(auto& it : _contigWindowSequences){

			const string& contigName = it.first;

			vector<vector<Window>>& windows = it.second;
			u_int64_t nbCorrectedWindows = 0;
			vector<DnaBitset2*> correctedWindows(windows.size());
		

			#pragma omp parallel num_threads(_nbCores)
			{

				std::unique_ptr<spoa::AlignmentEngine> alignmentEngine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -4);
				
				#pragma omp for
				for(size_t w=0; w<windows.size(); w++){

					vector<Window>& sequences = windows[w];
					
					u_int64_t wStart = w*_windowLength;
					u_int64_t wEnd = min(_contigSequences[contigName].size(), wStart+_windowLength);
					string contigOriginalSequence = _contigSequences[contigName].substr(wStart, wEnd-wStart);

					if(sequences.size() < 2){
							
						
						correctedWindows[w] = new DnaBitset2(contigOriginalSequence);
						//_logFile << "No sequences for window" << endl;
						continue;
					}
					
					#pragma omp atomic
					nbCorrectedWindows += 1;


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


    				std::vector<uint32_t> coverages;
    				string correctedSequence = graph.GenerateConsensus(&coverages);

					uint32_t average_coverage = (sequences.size()) / 2;

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
						correctedSequence = correctedSequence.substr(begin, end - begin + 1);
					}

					for(char letter : correctedSequence){
						checksum += letter;
					}

					correctedWindows[w] = new DnaBitset2(correctedSequence);
				}
			}

			if(nbCorrectedWindows == 0) continue;


			string contigSequence = "";
			for(size_t w=0; w<correctedWindows.size(); w++){
				if(correctedWindows[w] == nullptr) continue;
				char* seq = correctedWindows[w]->to_string();
				contigSequence += string(seq);
				free(seq);
				delete correctedWindows[w];
			}


			string header = _contigHeaders[contigName];
			if(_circularize){
				
				char lastChar = header[header.size()-1];
				if(lastChar == 'c' || lastChar == 'l'){
					header.pop_back();
				}

				if(_isContigCircular.find(contigName) == _isContigCircular.end()){
					header += 'l';
				}
				else{
					header += 'c';
				}
			}


			header = ">" + header + '\n';// ">ctg" + to_string(contigIndex) + '\n';
			//header += '\n';
			gzwrite(_outputContigFile, (const char*)&header[0], header.size());
			contigSequence +=  '\n';
			gzwrite(_outputContigFile, (const char*)&contigSequence[0], contigSequence.size());
			//_logFile << _contigHeaders[contigIndex] << " " << contigSequence.size() << endl;
		}

		_logFile << "Checksum: " << checksum << endl;
		*/
	}


	void addCorrectedWindow(bool success, DnaBitset2* seq, size_t contigIndexLocal, size_t windowIndexLocal, gzFile& outputContigFile){

		
		#pragma omp critical(addCorrectedWindow)
		{
			_currentContigs[contigIndexLocal].push_back({windowIndexLocal, seq, success});

			if(_currentContigs[contigIndexLocal].size() == _contigWindowSequences[contigIndexLocal].size()){
				dumpCorrectedContig(contigIndexLocal, outputContigFile);
			}
		}
		
	}

	void dumpCorrectedContig(const u_int32_t& contigIndex, gzFile& outputContigFile){


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
			//else if(length < 7500 && _contigCoverages[contigIndex] < 3){
			//	isValid = false;
			//}
			


			if(isValid){

				string header = _contigHeaders[contigIndex];
				
				
				
				header = Utils::split(header, '_')[0]; //Remove curator subpath suffix

				if(_isMetaMDBG){
					
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



	AlignmentBounds computeReadOverlapBounds(const Alignment& alignment1, const Alignment& alignment2, mm_tbuf_t* tbuf){



		const AlignmentBounds& alignment = minimap2Align(alignment1, alignment2, tbuf);
		
		return alignment;
	
	}


	AlignmentBounds computeReadOverlap(const Alignment& alignment1, const Alignment& alignment2, mm_tbuf_t* tbuf){
		
		return computeReadOverlapBounds(alignment1, alignment2, tbuf);
	}
	

	mm_idx_t* minimap2index(int w, int k, int is_hpc, int bucket_bits, const string& sequence){
		const char *seq = sequence.c_str();
		int len = sequence.size();
		const char *fake_name = "target";
		char *s;
		mm_idx_t *mi;
		s = (char*)calloc(len + 1, 1);
		memcpy(s, seq, len);
		mi = mm_idx_str(w, k, is_hpc, bucket_bits, 1, (const char**)&s, (const char**)&fake_name);
		free(s);
		return mi;
	}


	int32_t getEstimatedOverlapLength(const Alignment& prevAlignment, const Alignment& nextAlignment){

		
		int32_t prevLeftOverhang = 0;
		int32_t prevRightOverhang = 0;

		if(prevAlignment._strand){
			prevLeftOverhang = prevAlignment._readLength - prevAlignment._readEnd;
			prevRightOverhang = prevAlignment._readStart;
		}
		else{
			prevLeftOverhang = prevAlignment._readStart;
			prevRightOverhang = prevAlignment._readLength - prevAlignment._readEnd;
		}

		int32_t nextLeftOverhang = 0;
		int32_t nextRightOverhang = 0;

		if(nextAlignment._strand){
			nextLeftOverhang = nextAlignment._readLength - nextAlignment._readEnd;
			nextRightOverhang = nextAlignment._readStart;
		}
		else{
			nextLeftOverhang = nextAlignment._readStart;
			nextRightOverhang = nextAlignment._readLength - nextAlignment._readEnd;
		}

		int32_t prevStart = ((int32_t) prevAlignment._contigStart) - prevLeftOverhang;
		int32_t prevEnd = ((int32_t) prevAlignment._contigEnd) + prevRightOverhang;

		int32_t nextStart = ((int32_t) nextAlignment._contigStart) - nextLeftOverhang;
		int32_t nextEnd = ((int32_t) nextAlignment._contigEnd) + nextRightOverhang;


		return (int32_t)prevEnd - (int32_t)nextStart;
		
		return (int32_t)prevAlignment._contigEnd - (int32_t)nextAlignment._contigStart;
	}

	AlignmentBounds minimap2Align(const Alignment& alignment1, const Alignment& alignment2, mm_tbuf_t* tbuf){

		const ReadType& readIndex1 = alignment1._readIndex;
		const ReadType& readIndex2 = alignment2._readIndex;

		//cout << readIndex1 << " " << readIndex2 << endl;


		char* dnaStr1 = _readIndex_to_sequence[readIndex1]->to_string();
		string reference = string(dnaStr1, strlen(dnaStr1));
		free(dnaStr1);

		char* dnaStr2 = _readIndex_to_sequence[readIndex2]->to_string();
		string query = string(dnaStr2, strlen(dnaStr2));
		free(dnaStr2);
		
		AlignmentBounds alignment;
		alignment._referenceLength = reference.size();
		alignment._queryLength = query.size();

		//string preset = _minimap2Preset_ava
		mm_idxopt_t iopt;
		mm_mapopt_t mopt;
		mm_set_opt(0, &iopt, &mopt); //"ava-ont"
		mm_set_opt(_minimap2Preset_ava.c_str(), &iopt, &mopt); //"ava-ont"
		iopt.batch_size = 0x7fffffffffffffffL; //always build a uni-part index

		//" reprise: tester chacinign only, on a eu une seg fault la"
		mopt.flag |= MM_F_CIGAR; // perform alignment 
		mopt.min_chain_score = 500; //-m 500
		//iopt.bucket_bits = 14;

		//cout << iopt.w << endl;
		//cout << iopt.k << endl;
		//cout << iopt.bucket_bits << endl;
		//cout << (iopt.flag&1) << endl;
		//cout << mopt.min_chain_score << endl;

		//mm_idxopt_t idx_opt;
		mm_idx_t* mi = minimap2index(iopt.w, iopt.k, iopt.flag&1, iopt.bucket_bits, reference);
		mm_mapopt_update(&mopt, mi);
		//mopt.mid_occ = 1000; // don't filter high-occ seeds

		mm_reg1_t *reg;
		//mm_tbuf_t *tbuf = mm_tbuf_init();
		int j, i, n_reg;
		reg = mm_map(mi, query.size(), query.c_str(), &n_reg, tbuf, &mopt, 0); // get all hits for the query

		int32_t expectedOverlapLength = getEstimatedOverlapLength(alignment1, alignment2);
		int32_t expectedQueryStart = alignment2._readStart;
		int32_t expectedQueryEnd = expectedQueryStart + expectedOverlapLength;


		u_int32_t maxError = -1;
		u_int32_t minLength = 0;

		/*
		if(n_reg > 1){

			int nbValidAlignments = 0;
			cout << alignment1._readIndex << " " << alignment2._readIndex << endl;

			for (j = 0; j < n_reg; ++j) { // traverse hits and print them out


				mm_reg1_t *r = &reg[j];

				int32_t actualOverlapLength = r->qe - r->qs;

				int error = abs(expectedOverlapLength - actualOverlapLength) + abs(expectedQueryStart - r->qs);

				AlignmentBounds alignment2;
				alignment2._referenceLength = reference.size();
				alignment2._queryLength = query.size();
				alignment2._queryStart = r->qs;
				alignment2._queryEnd = r->qe;
				alignment2._referenceStart = r->rs;
				alignment2._referenceEnd = r->re;
				alignment2._isReversed = r->rev;
				

				cout << "Minimap2 align: " << r->qs << "\t" << r->qe << "\t\t" << r->rs << "\t" << r->re << "\t\t" << error << "\t" << isValidAlignment(alignment2, 0) << endl;

				if(isValidAlignment(alignment2, 0)){
					nbValidAlignments += 1;
				}
			}

			if(nbValidAlignments > 1) getchar();
		}
		*/

		for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
		
			mm_reg1_t *r = &reg[j];

			AlignmentBounds alignment2;
			alignment2._referenceLength = reference.size();
			alignment2._queryLength = query.size();
			alignment2._queryStart = r->qs;
			alignment2._queryEnd = r->qe;
			alignment2._referenceStart = r->rs;
			alignment2._referenceEnd = r->re;
			alignment2._isReversed = r->rev;
			
			if(!isValidAlignment(alignment2, 0)) continue;

			int32_t actualOverlapLength = r->qe - r->qs;

			int error = abs(expectedOverlapLength - actualOverlapLength) + abs(expectedQueryStart - r->qs);

			u_int32_t alignLength = min(r->qe - r->qs, r->re - r->rs);
			//cout << "Minimap2 align: " << r->qs << " " << r->qe << "    " << r->rs << " " << r->re << "    " << error << endl;
			

			//if(error < maxError){
			if(alignLength > minLength){
				minLength = alignLength;
				//maxError = error;

				alignment._queryStart = r->qs;
				alignment._queryEnd = r->qe;
				alignment._referenceStart = r->rs;
				alignment._referenceEnd = r->re;
				alignment._isReversed = r->rev;
			}
			/*
			qs, qe, rs, re
			assert(r->p); // with MM_F_CIGAR, this should not be NULL
			printf("%s\t%d\t%d\t%d\t%c\t", "query", query.size(), r->qs, r->qe, "+-"[r->rev]);
			printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\tcg:Z:", mi->seq[r->rid].name, mi->seq[r->rid].len, r->rs, r->re, r->mlen, r->blen, r->mapq);
			for (i = 0; i < r->p->n_cigar; ++i) // IMPORTANT: this gives the CIGAR in the aligned regions. NO soft/hard clippings!
				printf("%d%c", r->p->cigar[i]>>4, MM_CIGAR_STR[r->p->cigar[i]&0xf]);
			putchar('\n');
			*/
			free(r->p);
		}
		
		//cout << "Expected:       " << expectedQueryStart << " " << expectedQueryEnd << endl;
		//alignment.printAll();

		free(reg);
		//mm_tbuf_destroy(tbuf);
		mm_idx_destroy(mi);

		//if(alignment._isReversed){
			//cout << "omg" << endl;
			//getchar();
		//}
		//if(n_reg > 1) getchar();

		return alignment;
	}


	
	void loadUnitigIndex(){
		const string& unitigFilename = _tmpDir + "/unitig_data.txt.init";


		KminmerParserParallel parser(unitigFilename, _minimizerSize, 4, false, false, _nbCores);
		parser.parse(IndexUnitigFunctor(*this));
	}

	unordered_map<KmerVec, u_int32_t> _kminmer_to_unitigIndex;

	class IndexUnitigFunctor {

		public:

		ContigPolisher& _parent;

		IndexUnitigFunctor(ContigPolisher& parent) : _parent(parent){

		}

		IndexUnitigFunctor(const IndexUnitigFunctor& copy) : _parent(copy._parent){
			
		}

		~IndexUnitigFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {

			u_int64_t readIndex = kminmerList._readIndex;
			const vector<MinimizerType>& readMinimizers = kminmerList._readMinimizers;
			const vector<ReadKminmerComplete>& kminmersInfos = kminmerList._kminmersInfo;

			#pragma omp critical(IndexUnitigFunctor)
			{

				for(size_t i=0; i<kminmersInfos.size(); i++){

					const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];
					const KmerVec& vec = kminmerInfo._vec;

					_parent._kminmer_to_unitigIndex[vec] = readIndex;


				}
			}
		}
	};
	

	unordered_map<u_int128_t, u_int32_t> _kminmerAbundances;

	void loadKminmerAbundance(){


		ifstream kminmerAbundanceFile(_tmpDir + "/kminmerData_abundance_init.txt");

		while (true) {

			u_int128_t vecHash;
			u_int32_t abundance;

			kminmerAbundanceFile.read((char*)&vecHash, sizeof(vecHash));

			if(kminmerAbundanceFile.eof()) break;

			kminmerAbundanceFile.read((char*)&abundance, sizeof(abundance));
			
			//return false;

			//bool iseof = MDBG::readKminmerAbundance(vecHash, abundance, kminmerAbundanceFile);

			//if(iseof) break;

			if(abundance == 1) continue;
			_kminmerAbundances[vecHash] = abundance;

			//cout << "lala" << " " << abundance << endl;

		}

		kminmerAbundanceFile.close();
	}
	
	u_int128_t hash128(const KmerVec& vec) const{
		u_int128_t seed = vec._kmers.size();
		for(auto& i : vec._kmers) {
			seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		}
		return seed;
	}

	/*
	unordered_map<u_int32_t, u_int32_t> _unitigIndex_to_abundance;

	void loadUnitigAbundance(){


		ifstream abundanceFile(_tmpDir + "/unitigGraph.nodes.refined_abundances.bin.init");

		while (true) {

			u_int32_t unitigName;
			float abundance;

			abundanceFile.read((char*)&unitigName, sizeof(unitigName));
			
			if(abundanceFile.eof()) break;

			abundanceFile.read((char*)&abundance, sizeof(abundance));

			//cout << unitigName << " " << abundance << endl;
			_unitigIndex_to_abundance[unitigName*2] = abundance;
		}


		abundanceFile.close();

	}




	void loadContigData(){
		const string& unitigFilename = _tmpDir + "/contig_data.txt";


		KminmerParserParallel parser(unitigFilename, _minimizerSize, 4, false, false, _nbCores);
		parser.parse(ContigDataFunctor(*this));

		exit(1);
	}

	class ContigDataFunctor {

		public:

		ContigPolisher& _parent;

		ContigDataFunctor(ContigPolisher& parent) : _parent(parent){

		}

		ContigDataFunctor(const ContigDataFunctor& copy) : _parent(copy._parent){
			
		}

		~ContigDataFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {

			u_int64_t readIndex = kminmerList._readIndex;
			const vector<MinimizerType>& readMinimizers = kminmerList._readMinimizers;
			const vector<ReadKminmerComplete>& kminmersInfos = kminmerList._kminmersInfo;

			if(readIndex != 40703) return;

			for(size_t i=0; i<readMinimizers.size(); i++){
				cout << i << " " << readMinimizers[i] << endl;
			}

			return;

			#pragma omp critical(IndexUnitigFunctor)
			{

				cout << readIndex << " " << kminmersInfos.size() << endl;

				for(size_t i=0; i<kminmersInfos.size(); i++){

					const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];
					const KmerVec& vec = kminmerInfo._vec;

					cout << i << " " << _parent._kminmer_to_unitigIndex[vec] << endl;
					//cout << i << " " << _parent._kminmerAbundances.size() << " " << _parent._kminmerAbundances[hash128(vec)] << endl;
					//getchar();
					//_parent._kminmer_to_unitigIndex[vec] = readIndex;
				}
			}
		}

		u_int128_t hash128(const KmerVec& vec) const{
			u_int128_t seed = vec._kmers.size();
			for(auto& i : vec._kmers) {
				seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
			}
			return seed;
		}

	};
	*/

};	

#endif 


