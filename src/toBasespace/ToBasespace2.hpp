



/*

overlapOnTheReference: on pourrait utiliser un offset negatif pour essayer de mapper en ava-ont les reads qui map au bord d'une missassembly critique (comme le overlapOnTheReference du contigCurator initial)
	- on a un +100 hardcoded maintenant a ajuster pour aller aligner plus loin dans les overlap sur la reference

- continuer d'intégrer les modifs de la v1.3
	- read correction
	- fix 50000 bubble tip cleaning
	- finir le fast ToBasespaceFast
	- MinimizerType = 32 bits
	- HashType = 64 bits
	- verifier les repeated minimizers (surement a enlever ou ajuster)
	- read correction (actuellement faites en low density), on est repasser en high density mais pas sur que ça soit utile (il faudra essayer low density avec les minimizer de taille 13)
	- read correction: reactiver les truc de detection des chimeres
	- read correction alignLowDensity: est-ce qu'on pourrait enlever la queue qui force a ecrire les alignment dans le bon ordre? SInon ça stack en memoire 
	- contig derep: utiliser la methode actuel, mais en plus on ajouter la limite de ref coverage (ajuster la taille minimum des contig a derep à 200 000, enlever la limite par identity)

- TODO: HiFi mode et skip correction mode dans le toBasespace2
- verifier que le datatype converti en int fonctionne bien avec l'enum _dataType == DataType::Ont...
- alignReads: ecrire les alignment sur le disque, puis load

- renforcer le decision making pour tag isCircular
- rafiner les bord des contigs pour les contig circular + qui ont des mapping sur les deux bouts

- peut etre trop de start contig détecté comme érroné a cause du isErroneous qui n'a pas assez de reads mapping au debut du read? (<=2 cov maintenant)
- enlever les length filter dans ContigPolisher
- utiliser le Covered fraction pour tracker les derniers contigs mal convertis
- méthode additionel: on peut coller les differents fragments d'un contig en remappant les reads dessus (minimap2), et en verifiant qu'un read fait bien la jonction
- verifier qu'on est bien deterministe (le readIndex doit trancher si egalité dans le ssorting)
- ReadMapping2: potentiellement des variable a clean (contigSize par exemple)

- toujours des misassemblies, surement dans le isErronous besoin d'un seuil de couverture, le <= 2 doit etre trop faible pour les organisme tres abondant

- supprimer le readPartitionDir dans la version finale
- usedReads: pas besoin de les output dans la version finale

- read correction avec nouvelle indexation comme dans le toBasespace2 avec du sorting
	- assez lent sur un full soil sample, mais peut etre qu'on peut encore plus trier l'index (actellement on tri juste par minimizers) pour ne pas avoir a trier pendant le alignReads
	- rappel: une fois qu'on a trier l'index, on peut en théorie enlever les minimizers de l'index (avant de creer la lookup table)

./bin/metaMDBG toBasespace_ont ~/appa/run/correction/test_deterministic/test_zymo_7/asm/tmp/ ~/appa/run/correction/test_deterministic/test_zymo_7/asm/tmp/contig_data_init_small.txt.norepeats ~/appa/run/correction/test_deterministic/test_zymo_7/asm/tmp/contigs_uncorrected.fasta.gz ~/appa/run/correction/test_deterministic/test_zymo_7/asm/tmp/input.txt --threads 32
conda run -n quast quast ~/appa/run/correction/test_deterministic/test_zymo_7/asm/tmp/contigs_uncorrected.fasta.gz  -o /pasteur/appa/scratch/gbenoit/run/overlap_test_201/quast/ -r Roseburia_hominis.fa
python3 ~/zeus/scripts/nanopore/evaluation/run_metaquast.py ~/appa/run/correction/test_deterministic/test_zymo_7/asm/contigs.fasta.gz ~/appa/data/nanopore/mock/Zymo_D6331/references_filenames.txt ~/zeus/tmp/metaquast/ 32
grep -H "# misassemblies" /pasteur/appa/homes/gbenoit/zeus/tmp/metaquast_12/runs_per_reference//report.tsv

*/

#ifndef MDBG_METAG_TOBASESPACE2
#define MDBG_METAG_TOBASESPACE2

#include "../Commons.hpp"
//#include "../utils/edlib.h"
//#include "../utils/spoa/include/spoa/spoa.hpp"
#include "../utils/DnaBitset.hpp"
#include "ReadPartitionner.hpp"
#include "ContigPolisher.hpp"
#include "ContigDerep.hpp"
#include "ContigTrimmer.hpp"
//#include "ContigRepeatRemover.hpp"
#include "ReadVsContigMapper.hpp"

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

	u_int32_t _debug_contigIndex = -1;

	constexpr static const int64_t _minOverlap = 500;
	constexpr static const float _intFrac = 0.8;
	constexpr static const int64_t _maxHang = 200;
	constexpr static const int64_t _startEndOffset = 50;
	constexpr static const float _minOverlapIdentity = 0.9;
	//constexpr static const float _minContigLength = 50;


	//struct MinimizerPosition{
	//	MinimizerType _minimizer;
	//	u_int32_t _contigIndex;
	//	u_int32_t _position;
	//};

	//vector<MinimizerPosition> _allMinimizerPositions2;
	//vector<u_int32_t> _contigSizes;
	unordered_map<ReadType, string> _debug_readHeaders;


	string _inputFilename;
	string _inputFilenameContig;
	string _inputDir;

	
	Parameters _params;
	string _readPartitionDir;
	u_int32_t _nbPartitions;

	BGZF* _usedReadFile;
	int _maxMemoryGB;
	string _outputContigFilename;
	//bool _isFirstPass;
	//bool _isOutputFasta;

	//string _filename_outputContigs;
	//MDBG* _mdbg;
	//MinimizerParser* _minimizerParser;
	int _nbCores;
	long double _averageDistanceBetweenMinimizers;
	bool _skipCorrection;
	ofstream _outputMinimizerContigFile;



	u_int64_t _nbContigs;
	u_int64_t _nbBps;
	//u_int64_t _checksum;
	//gzFile _basespaceContigFile;
	u_int64_t _contigIndex;
	bool _print_debug_main = false;
	bool _print_debug = false;
	ofstream _alignmentOutputFile;


	string _minimap2Preset_map;
	string _minimap2Preset_ava;
	int _minContigLength;
	int _minContigCoverage;
	bool _isFinalAssembly;
	ankerl::unordered_dense::set<UnitigType> _validContigIndexes;

	ankerl::unordered_dense::map<ReadType, string> _mReads;
	
	u_int64_t _checksum_curatedContigs;
	u_int64_t _checksum_curatedContigs_total;

	struct ContigData{
		float _coverage;
		vector<ReadMapping2> _alignments;
		ankerl::unordered_dense::map<ReadType, u_int32_t> _readIndex_to_i;
		ankerl::unordered_dense::map<ReadType, bool> _isErroneousReadCached;
		ankerl::unordered_dense::map<pair<ReadType, ReadType>, vector<AlignmentBounds>> _cachedAlignments;
	};

	enum CoverageRegionType{
		Low,
		Normal,
		High,
	};

	struct CoverageRegion{

		//const static int low = 0;
		//const static int normal = 1;
		//const static int high = 2;

		size_t _startIndex;
		size_t _endIndex;
		int _coverageType;
		//u_int64_t _coverage;
		//u_int64_t _length;

		int64_t getLength() const{
			return _endIndex - _startIndex;
		}

		bool isSupportedByRead(const vector<u_int64_t>& coveragesMapping) const{

			for(size_t i=_startIndex; i<_endIndex; i++){
				if(coveragesMapping[i] == 0) return false;
			}

			return true;
		}

	};

	//phmap::parallel_flat_hash_map<ReadType, ReadMapping2> _readIndex_to_bestContigMap;
	ankerl::unordered_dense::map<ReadType, ContigData> _contigIndex_to_alignments;

	ToBasespace2(): Tool (){

	}

	//gzFile _outputFile_left;
	//gzFile _outputFile_right;
	//std::unique_ptr<spoa::AlignmentEngine> _alignment_engine;

	void parseArgs(int argc, char* argv[]){

		args::ArgumentParser parser("toBasespace", ""); //"This is a test program.", "This goes after the options."
		args::Positional<std::string> arg_outputDir(parser, "outputDir", "", args::Options::Required);
		//args::Positional<std::string> arg_inputContigFilename(parser, "inputContigFilename", "", args::Options::Required);
		//args::Positional<std::string> arg_outputContigFilename(parser, "outputContigFilename", "", args::Options::Required);
		//args::Positional<std::string> arg_inputReadFilename(parser, "inputReadFilename", "", args::Options::Required);
		args::Flag arg_skipCorrection(parser, "", "Skip read correction", {"skip-correction"});
		args::ValueFlag<int> arg_maxMemory(parser, "", "Maximum memory usage", {ARG_MAX_MEMORY}, 8);
		args::ValueFlag<float> arg_minContigLength(parser, "", "Minimum contig length", {ARG_MIN_CONTIG_LENGTH}, 50);
		args::ValueFlag<float> arg_minContigCoverage(parser, "", "Minimum contig coverage", {ARG_MIN_CONTIG_COVERAGE}, 1);
		args::Flag arg_finalAssembly(parser, "", "Solve final graph", {"final"});
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
		//_inputFilenameContig = _inputDir + "/contig_data_init_small.txt.nooverlaps";  //args::get(arg_inputContigFilename);
		_inputFilenameContig = _inputDir + "/contig_data_init_small.txt.norepeats";  //args::get(arg_inputContigFilename);
		//_filename_outputContigs = _inputDir + "/contigs_uncorrected.fasta.gz"; //args::get(arg_outputContigFilename);
		_inputFilename = _inputDir + "/input.txt"; //args::get(arg_inputReadFilename);
		_nbCores = args::get(arg_nbCores);

		_maxMemoryGB = args::get(arg_maxMemory);
		_maxMemoryGB = max(4, _maxMemoryGB);

		_skipCorrection = false;
		if(arg_skipCorrection){
			_skipCorrection = true;
		}

		_isFinalAssembly = false;
		if(arg_finalAssembly){
			_isFinalAssembly = true;
		}

		_minContigLength = args::get(arg_minContigLength);
		_minContigCoverage = args::get(arg_minContigCoverage);


		/*
		string filename_parameters = _inputDir + "/parameters.gz";
		gzFile file_parameters = gzopen(filename_parameters.c_str(),"rb");
		gzread(file_parameters, (char*)&_minimizerSize, sizeof(_minimizerSize));
		gzread(file_parameters, (char*)&_kminmerSize, sizeof(_kminmerSize));
		gzread(file_parameters, (char*)&_minimizerDensity_assembly, sizeof(_minimizerDensity_assembly));
		gzread(file_parameters, (char*)&_kminmerSizeFirst, sizeof(_kminmerSizeFirst));
		gzread(file_parameters, (char*)&_minimizerSpacingMean, sizeof(_minimizerSpacingMean));
		gzread(file_parameters, (char*)&_kminmerLengthMean, sizeof(_kminmerLengthMean));
		gzread(file_parameters, (char*)&_kminmerOverlapMean, sizeof(_kminmerOverlapMean));
		gzread(file_parameters, (char*)&_kminmerSizePrev, sizeof(_kminmerSizePrev));
		gzread(file_parameters, (char*)&_kminmerSizeLast, sizeof(_kminmerSizeLast));
		gzread(file_parameters, (char*)&_meanReadLength, sizeof(_meanReadLength));
		gzread(file_parameters, (char*)&_minimizerDensity_correction, sizeof(_minimizerDensity_correction));
		gzread(file_parameters, (char*)&_useHomopolymerCompression, sizeof(_useHomopolymerCompression));
		gzread(file_parameters, (char*)&_dataType, sizeof(_dataType));
		gzread(file_parameters, (char*)&_snpmerSize, sizeof(_snpmerSize));
		gzclose(file_parameters);
		*/

		_params.load(_inputDir + "/parameters.gz");
		_params._kminmerSize = 2;
		_averageDistanceBetweenMinimizers = 1.0 / _params._minimizerDensity_assembly;
		_readPartitionDir = _inputDir + "/_polish_readPartitions/";
		_outputContigFilename = _inputDir + "/contigs.fasta.gz";
		
		if(!fs::exists(_readPartitionDir)){
			fs::create_directories(_readPartitionDir);
		}

		openLogFile(_inputDir);
		/*
		Logger::get().debug() << "";
		Logger::get().debug() << "Input dir: " << _inputDir;
		//_logFile << "Output filename: " << _outputFilename << endl;
		Logger::get().debug() << "Minimizer length: " << _minimizerSize;
		Logger::get().debug() << "Kminmer length: " << _kminmerSize;
		Logger::get().debug() << "Density: " << _minimizerDensity_assembly;
		Logger::get().debug() << "";

		//_filename_outputContigs = _inputDir + "/contigs.fasta.gz";
		//_filename_outputContigs = _inputFilenameContig + ".fasta.gz"; //_inputDir + "/tmpContigs.fasta.gz";
		//_minimizerParser = new MinimizerParser(_minimizerSize, _minimizerDensity);
		*/
	}



    void execute (){
		
		//cout << "Set minimap2 in no-kalloc mode" << endl;
		mm_dbg_flag |= MM_DBG_NO_KALLOC; // --no-kalloc

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

		if(_params._dataType == SequencingPlatformType::HiFi){
			//cout << "Data type is HiFi" << endl;
			_minimap2Preset_map = "map-hifi";
			_minimap2Preset_ava = "ava-pb";
			//_useHpc = true;
		}
		else if(_params._dataType == SequencingPlatformType::Nanopore){
			//cout << "Data type is ONT" << endl;
			_minimap2Preset_map = "map-ont";
			_minimap2Preset_ava = "ava-ont";
			//_useHpc = false;
		}
		
		//_nbPartitions = 30;
		//removeContigRepeats();
		//exit(1);
		//derepContigs(_readPartitionDir + "/contigs_polished.fasta.gz", _readPartitionDir + "/contigs_derep.fasta.gz", 10);
		//trimContigs(_readPartitionDir + "/contigs_derep.fasta.gz", _readPartitionDir + "/contigs_trimmed.fasta.gz", minimapBatchSize);
		//exit(1);
		

		
		//cout << "load read headers" << endl;
		//ReadParserParallel readParser(_inputDir + "/input.txt", false, false, 1);
		//readParser.parse(ReadParserFunctor(*this));
		
		//_checksum_contigs = 0;
		/*
		Logger::get().debug() << "Loading minimizers";
		loadContigs();
		//loadReadsTest();
		Logger::get().debug() << "\tDone: " << (peakrss() / 1024.0 / 1024.0 / 1024.0);

		//cout << _allMinimizerPositions.size() << endl;
		Logger::get().debug() << "Nb minimizers: " << _allMinimizerPositions.size();


		//std::sort(_allMinimizerPositions.begin(), _allMinimizerPositions.end(), [](const MinimizerPosition& a, const MinimizerPosition& b){
		//	return a._minimizer < b._minimizer;
		//});
		Logger::get().debug() << "Sorting minimizers";
		Commons::sortParallel(_allMinimizerPositions, _allMinimizerPositions.size(), _nbCores);
		//sortIndex();
		Logger::get().debug() << "\tDone: " << (peakrss() / 1024.0 / 1024.0 / 1024.0);

		//cout << _allMinimizerPositions.size() << endl;
		//cout << "todo: enelver les minimizer de la table des positions: _allMinimizerPositions" << endl;
		//cout << "Attention: if(anchors.size() < 30) return; a reduire a 3" << endl;

		Logger::get().debug() << "Creating minimizer lookup table";
		createLookupTable();
		Logger::get().debug() << "\tDone: " << (peakrss() / 1024.0 / 1024.0 / 1024.0);


		//cout << "Peak memory: " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << endl;

		Logger::get().debug() << "Aligning reads";
		alignReads();
		Logger::get().debug() << "\tDone: " << (peakrss() / 1024.0 / 1024.0 / 1024.0);
		//cout << "Nb queries: " << _loul << endl;

		//exit(1);

		//cout << "Peak memory: " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << endl;

		vector<MinimizerPosition>().swap(_allMinimizerPositions); //_allMinimizerPositions.clear();
		ankerl::unordered_dense::map<MinimizerPairType, pair<u_int64_t, u_int64_t>>().swap(_minimizerLookupTable); //_minimizerLookupTable.clear();
		
		*/

		
		alignReads();
		


		if(_debug_contigIndex != -1){
			cout << "Skip read partitionning" << endl;
			_nbPartitions = 0;
			for(int i=0; i<10000; i++){
				if(fs::exists(_readPartitionDir + "/" + to_string(i) + "_contigs.bin")){
					_nbPartitions += 1;
				}
				else{
					break;
				}
			}
			cout << "Nb partitions: " << _nbPartitions << endl;
		}
		else{
			Logger::get().debug() << "";
			Logger::get().debug() << "Partitionning reads";
			if(_print_debug_main) cout << "Partitionning reads" << endl;
			partitionReads();
		}


	
		
		
		const string& outputContigFilename_polished = _readPartitionDir + "/contigs_polished.fasta.gz";

		_contigIndex = 0;
		//_basespaceContigFile = gzopen(_filename_outputContigs.c_str(),"wb");
		
		_usedReadFile = bgzf_open((_readPartitionDir + "/usedReads.fasta.gz").c_str(), "w1"); 
		_outputMinimizerContigFile = ofstream(_inputDir + "/contig_data_final.bin");
		
		createBaseContigs(_inputFilenameContig, outputContigFilename_polished);
		//createBaseContigs(inputFilenameContigSmall);

		//exit(1);

		Logger::get().debug() << "Checksum nb contigs: " << _contigIndex;
		Logger::get().debug() << "Checksum total bases: " << _nbBps;
		Logger::get().debug() << "Checksum: " << _checksum_curatedContigs_total;
		Logger::get().debug() << "Peak memory: " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";

		//gzclose(_basespaceContigFile);
		_outputMinimizerContigFile.close();
		
		bgzf_close(_usedReadFile);
		
		//Logger::get().debug() << "";
		//Logger::get().debug() << "Remove repeats";
		//removeContigRepeats();

		


		
		//Logger::get().debug() << "";
		//Logger::get().debug() << "Polishing contigs";
		//const string& outputContigFilename_polished = _readPartitionDir + "/contigs_polished.fasta.gz";
		//polishContigs(outputContigFilename_polished, minimapBatchSize);


		//cout << "Peak memory: " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << endl;
		
	
		//exit(1);



		Logger::get().debug() << "";
		Logger::get().debug() << "Dereplicating contigs";
		//const string& derepContigFilename1 = _tmpDir + "/contigs_derep_1.fasta.gz";
		const string& outputContigFilename_derep = _readPartitionDir + "/contigs_derep.fasta.gz";
		derepContigs(outputContigFilename_polished, outputContigFilename_derep, minimapBatchSize);
		
		
	
		//cout << "Peak memory: " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << endl;
		
		Logger::get().debug() << "";
		Logger::get().debug() << "Trimming contigs";
		const string& outputContigFilename_trimmed = _readPartitionDir + "/contigs_trimmed.fasta.gz";
		trimContigs(outputContigFilename_derep, outputContigFilename_trimmed, minimapBatchSize);



		//cout << "Peak memory: " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << endl;

		Logger::get().debug() << "";
		Logger::get().debug() << "Moving final contigs to destination";
		fs::rename(outputContigFilename_trimmed, _outputContigFilename);
		//fs::rename(outputContigFilename_derep, _outputContigFilename);

		if(!_isFinalAssembly){
			rewriteFinalContigFileWithoutErrors();
			//fs::copy(_inputDir + "/contig_data_final.bin", _inputDir + "/contigGraph/contig_data_final.bin", fs::copy_options::overwrite_existing);
		}

		if(!DEBUG) fs::remove_all(_readPartitionDir);

	}

	void rewriteFinalContigFileWithoutErrors(){

		ofstream outputFile(_inputDir + "/contigGraph/contig_data_final.bin");

		indexValidContigs();

		KminmerParserParallel parser(_inputDir + "/contig_data_final.bin", 0, 0, false, false, 1);
		parser.parseSequences(RewriteValidContigFunctor(*this, outputFile));

		outputFile.close();
		_validContigIndexes.clear();
	}


	void indexValidContigs(){

		auto fp = std::bind(&ToBasespace2::indexValidContigs_read, this, std::placeholders::_1);
		ReadParser readParser(_outputContigFilename, true, false);
		readParser.parse(fp);
		
	}
	
	void indexValidContigs_read(const Read& read){

		string contigName = Utils::shortenHeader(read._header);
		contigName.erase(0, 3); //remove "ctg"
		u_int32_t contigIndex = stoull(contigName);

		_validContigIndexes.insert(contigIndex);
	}



	class RewriteValidContigFunctor {

		public:

		ToBasespace2& _parent;
		ofstream& _outputFile;

		RewriteValidContigFunctor(ToBasespace2& parent, ofstream& outputFile) : _parent(parent), _outputFile(outputFile){
		}

		RewriteValidContigFunctor(const RewriteValidContigFunctor& copy) : _parent(copy._parent), _outputFile(copy._outputFile){
		}

		~RewriteValidContigFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {

			bool isCircular = kminmerList._isCircular;
			const UnitigType contigIndex = kminmerList._readIndex;

			vector<MinimizerType> minimizers = kminmerList._readMinimizers;

			if(_parent._validContigIndexes.find(contigIndex) == _parent._validContigIndexes.end()){
				minimizers.clear();
			}

			u_int32_t contigSize = minimizers.size();
			_outputFile.write((const char*)&contigSize, sizeof(contigSize));
			_outputFile.write((const char*)&isCircular, sizeof(isCircular));
			_outputFile.write((const char*)&minimizers[0], contigSize*sizeof(MinimizerType));
		}
	};

	void alignReads(){

		Logger::get().debug() << "Align reads vs contigs";
		auto start = high_resolution_clock::now();


		ReadVsContigMapper mapper(_inputDir + "/read_data_init.txt", _inputFilenameContig, _readPartitionDir + "/readsVsContigsAlignments.bin", 1.0 / _params._minimizerDensity_assembly, _nbCores);
		mapper.execute();

		Logger::get().debug() << "\tDone: " << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";

	}
	

	class ReadParserFunctor {

		public:

		ToBasespace2& _parent;

		ReadParserFunctor(ToBasespace2& parent) : _parent(parent){
		}

		ReadParserFunctor(const ReadParserFunctor& copy) : _parent(copy._parent){
		}

		~ReadParserFunctor(){
		}

		
		void operator () (Read& read) {

			/*
			//if(read._index == 667){
			//	cout << read._header << endl;
			//	getchar();
			//}
			
			//if(Utils::shortenHeader(read._header) == "5f33be6b-b155-4e3f-8712-31a1b44189b4"){
			if(read._index == 2724062){
				cout << read._index << endl;
				cout << read._header << endl;
				cout << read._seq << endl;
				//exit(1);
			}

			//if(Utils::shortenHeader(read._header) == "5f33be6b-b155-4e3f-8712-31a1b44189b4"){
			if(read._index == 866271){
				cout << read._index << endl;
				cout << read._header << endl;
				cout << read._seq << endl;
				//exit(1);
			}
			*/

			_parent._debug_readHeaders[read._index] = read._header;
		}
	};


	void partitionReads(){
		
		auto start = high_resolution_clock::now();

		ReadPartitionner readPartitionner(_inputDir + "/input.txt", _inputFilenameContig, _readPartitionDir, _averageDistanceBetweenMinimizers, _maxMemoryGB, _nbCores);
		//readPartitionner.computeContigCoverages(_contigIndex_to_alignments, _contigSizes);
		readPartitionner.execute();

		_nbPartitions = readPartitionner._nbPartitions;
		
		Logger::get().debug() << "\tDone " << " (" << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s) " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << "GB";
	}

	void loadAlignments(){

		ankerl::unordered_dense::map<ReadType, ContigData>().swap(_contigIndex_to_alignments);
		//_contigIndex_to_alignments.clear();
		/*
		for(const auto& it : _readIndex_to_bestContigMap){
			const ReadMapping2& al = it.second;
			_contigIndex_to_alignments[al._contigIndex].push_back(al);
		}

		_readIndex_to_bestContigMap.clear();
		*/

		ifstream alignmentFile(_readPartitionDir + "/readsVsContigsAlignments.bin");

		while(true){

			ReadMapping2 alignment;
			bool isEof = alignment.read(alignmentFile);

			if(isEof) break;

			if(_mReads.find(alignment._readIndex) == _mReads.end()) continue;

			if(_contigIndex_to_alignments.find(alignment._contigIndex) == _contigIndex_to_alignments.end()){
				_contigIndex_to_alignments[alignment._contigIndex] = ContigData();
			}

			_contigIndex_to_alignments[alignment._contigIndex]._alignments.push_back(alignment);

		}

		alignmentFile.close();

	}

	unordered_map<ReadType, ReadMapping2> _debug_readIndex_to_bestAlignment;

	void loadAlignmentsDebug(){

		ankerl::unordered_dense::map<ReadType, ContigData>().swap(_contigIndex_to_alignments);
		//_contigIndex_to_alignments.clear();
		/*
		for(const auto& it : _readIndex_to_bestContigMap){
			const ReadMapping2& al = it.second;
			_contigIndex_to_alignments[al._contigIndex].push_back(al);
		}

		_readIndex_to_bestContigMap.clear();
		*/

		ifstream alignmentFile(_readPartitionDir + "/readsVsContigsAlignments.bin");

		while(true){

			ReadMapping2 alignment;
			bool isEof = alignment.read(alignmentFile);

			if(isEof) break;

			_debug_readIndex_to_bestAlignment[alignment._readIndex] = alignment;

		}

		alignmentFile.close();

	}

	void clearPass(){

		//for(const auto& it: _mReads){
		//	delete it.second;
		//}
		
		
		//_contigIndex_to_breakpoints.clear();

		ankerl::unordered_dense::map<ReadType, string>().swap(_mReads); //_mReads.clear();
		ankerl::unordered_dense::map<ReadType, ContigData>().swap(_contigIndex_to_alignments); //_contigIndex_to_alignments.clear();
	}

	void loadReads(const string& readFilename){

		//KminmerParserParallelOnt parser(_inputDir + "/read_data_init.txt.ont", _minimizerSize, _kminmerSize, false, true, _nbCores);
		//parser.parse(LoadMinimizerReadsFunctor(*this));

		ReadParserParallel readParser(readFilename, true, false, 1);
		readParser.parse(LoadReadsFunctor(*this));
	}

	//vector<string> _mReads4;

	class LoadReadsFunctor {

		public:

		ToBasespace2& _parent;

		LoadReadsFunctor(ToBasespace2& parent) : _parent(parent){
		}

		LoadReadsFunctor(const LoadReadsFunctor& copy) : _parent(copy._parent){
		}

		~LoadReadsFunctor(){
		}

		void operator () (Read& read) {

			ReadType readIndex = read._index;

			if(readIndex % 100000 == 0) Logger::get().debug() << "\t\t" << readIndex;

			readIndex = stoull(read._header);

			_parent._mReads[readIndex] = read._seq;


			//_parent._mReads4.push_back(read._seq);
			//if(_parent._readIndex_to_bestContigMap.find(readIndex) == _parent._readIndex_to_bestContigMap.end()) return;

			/*
			//_parent._mReads[readIndex] = seqBinary;
			BinaryRead binaryRead = {DnaBitset2(read._seq), read._qual};

			_parent._mReads.lazy_emplace_l(readIndex,
			[](BinaryReadMap::value_type& v) { // key exist
				//v.second.;
				//v.second[0] += 1; //Increment kminmer count
				//if(std::find(v.second.begin(), v.second.end(), readIndex) == v.second.end()){
				//	v.second.push_back(readIndex);
				//}
			},           
			[&readIndex, &binaryRead](const BinaryReadMap::constructor& ctor) { // key inserted
				
				ctor(readIndex, binaryRead); 

			}); // construct value_type in place when key not present
			*/
			

		}
	};



	BGZF* _outputContigFile;
	bool _foundDebugContig;
	//u_int32_t _currentPartition;
	//u_int64_t _nbPrecomputedAlignments;

	void createBaseContigs(const string contigFilename, const string outputContigFilenamePolished){

		
		BGZF* outputContigFile_polished = bgzf_open(outputContigFilenamePolished.c_str(),"w1");
		ofstream file_contigHeaders(_inputDir + "/contigHeaders.bin");

		auto startTotal = high_resolution_clock::now();

		_nbBps = 0;
		_checksum_curatedContigs_total = 0;

		Logger::get().debug() << "";
		Logger::get().debug() << "Creating basespace contigs: " << contigFilename;

		for(u_int32_t i=0; i<_nbPartitions; i++){

			//if(i != 6){
			//	cout << "skip" << endl;
			//	continue;
			//}

			auto start = high_resolution_clock::now();
			Logger::get().debug() << "";
			Logger::get().debug() << "";
			Logger::get().debug() << "Process partition: " << i << "/" << _nbPartitions;
			if(_print_debug_main) cout << "Process partition: " << i << "/" << _nbPartitions << endl;



			_checksum_curatedContigs = 0;
			clearPass();
			//_currentPartition = partition;
			
			const string readFilenameCompressed = _readPartitionDir + "/" + to_string(i) + "_reads.gz";
			const string readFilename = _readPartitionDir + "/" + to_string(i) + "_reads";
			const string contigFilename = _readPartitionDir + "/" + to_string(i) + "_contigs.bin";
			//const string outputContigFilename = _readPartitionDir + "/" + to_string(i) + "_contigs_uncorrected.fasta.gz";
			const string outputContigFilename = _readPartitionDir + "/contigs_uncorrected.fasta.gz";

			//if(fs::exists(outputContigFilename)) continue;

			if(_debug_contigIndex != -1){
				searchDebugContig(contigFilename);
				if(!_foundDebugContig) continue;
			}

			//cout << "a" << endl;
			//auto start2 = high_resolution_clock::now();
			//string command = "bgzip -d -f -k --threads " + to_string(_nbCores) + " " + readFilenameCompressed ;
			//system(command.c_str());
			//Logger::get().debug() << "\tDone " << " (" << duration_cast<seconds>(high_resolution_clock::now() - start2).count() << "s)";

			//cout << "b" << endl;
			
			Logger::get().debug() << "";
			Logger::get().debug() << "\tDecompressing read file";
			auto start2 = high_resolution_clock::now();
			Utils::bgzip_decompress(readFilenameCompressed, readFilename, _nbCores);
			Logger::get().debug() << "\tDone " << " (" << duration_cast<seconds>(high_resolution_clock::now() - start2).count() << "s) " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";

			
			Logger::get().debug() << "";
			Logger::get().debug() << "\tLoading reads";
			start2 = high_resolution_clock::now();
			loadReads(readFilename);
			Logger::get().debug() << "\tDone " << " (" << duration_cast<seconds>(high_resolution_clock::now() - start2).count() << "s) " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";

			Logger::get().debug() << "";
			Logger::get().debug() << "\tLoading alignments";
			start2 = high_resolution_clock::now();
			loadAlignments();
			Logger::get().debug() << "\tDone " << " (" << duration_cast<seconds>(high_resolution_clock::now() - start2).count() << "s) " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";

			

			Logger::get().debug() << "";
			Logger::get().debug() << "\tPrecompute reads vs reads alignments";
			start2 = high_resolution_clock::now();
			precomputeAlignments(contigFilename);
			Logger::get().debug() << "\tDone " << " (" << duration_cast<seconds>(high_resolution_clock::now() - start2).count() << "s) " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";
			
			
			Logger::get().debug() << "";
			Logger::get().debug() << "\tCreating contigs";
			start2 = high_resolution_clock::now();
			//cout << "Create contigs" << endl;
			
			_outputContigFile = bgzf_open(outputContigFilename.c_str(), "w1");

			//cout << "single core here" << endl;
			KminmerParserParallel parser(contigFilename, _params._minimizerSize, _params._kminmerSize, false, false, _nbCores);
			parser._hasContigIndex = true;
			parser.parseSequences(CreateBaseContigsFunctor(*this));
			
			bgzf_close(_outputContigFile);

			Logger::get().debug() << "\tDone " << " (" << duration_cast<seconds>(high_resolution_clock::now() - start2).count() << "s) " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";
			//cout << "\tDone " << " (" << duration_cast<seconds>(high_resolution_clock::now() - start2).count() << "s)" << endl;
			




			clearPass();
			
			ContigPolisher contigPolisher(_inputDir, _readPartitionDir, outputContigFilename, _params._dataType, true, _nbCores, _nbPartitions, _minimap2Preset_map, _minContigLength, _minContigCoverage);
			contigPolisher.execute2(i, outputContigFile_polished, file_contigHeaders);




			fs::remove(readFilename);

			Logger::get().debug() << "\tChecksum partition " << i <<  ": " << _checksum_curatedContigs;
			Logger::get().debug() << "\tDone " << " (" << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s) " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << "GB";
		}
		//cout << "todo: missing last minimizer" << endl;


		file_contigHeaders.close();
		bgzf_close(outputContigFile_polished);

		Logger::get().debug() << "Nb contigs: " << (_contigIndex);
		Logger::get().debug() << "Nb bps: " << _nbBps;
		//Logger::get().debug() << "Checksum: " << _checksum_contigs;
		Logger::get().debug() << "Done " << " (" << duration_cast<seconds>(high_resolution_clock::now() - startTotal).count() << "s) " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << "GB";

	}

	void precomputeAlignments(const string contigFilename){

		KminmerParserParallel parser2(contigFilename, _params._minimizerSize, _params._kminmerSize, false, false, _nbCores);
		parser2._hasContigIndex = true;
		parser2.parseSequences(AlignmentCacheFunctor(*this));

		//cout << "Nb cached alignments: " << _alignmentsToPrecompute.size() << endl;

		//_nbPrecomputedAlignments = 0;
		//start2 = high_resolution_clock::now();
		//vector<pair<ReadType, ReadType>> cachedAlignments;
		//for(const auto& it : _readsVsReadsAlignmentsCached){
		//	cachedAlignments.push_back(it.first);
		//}
	
		#pragma omp parallel num_threads(_nbCores)
		{

			mm_tbuf_t* tbuf = mm_tbuf_init();

			#pragma omp for
			for (u_int64_t i=0; i < _alignmentsToPrecompute.size(); i++) {

				const PrecomputeAlignmentData alignmentData = _alignmentsToPrecompute[i];
				//const pair<ReadType, ReadType> readPair = cachedAlignments[i];

				const string& read1 = _mReads[alignmentData._readIndex1];
				//string read1 = string(dnaStr1, strlen(dnaStr1));
				//free(dnaStr1);

				const string& read2 = _mReads[alignmentData._readIndex2];
				//string read2 = string(dnaStr2, strlen(dnaStr2));
				//free(dnaStr2);

				mm_idx_t* mi = minimap2_indexRead(read1);

				vector<AlignmentBounds> allAlignments;
				computeAlignment(read1, read2, tbuf, _minimap2Preset_map, true, ToBasespace2::_minOverlap, mi, allAlignments);
				//if(alignment._queryStart == -1) continue;
				
				mm_idx_destroy(mi);


				ContigData& contigData = _contigIndex_to_alignments[alignmentData._contigIndex];
				u_int32_t ii = contigData._readIndex_to_i[alignmentData._readIndex1];
				ankerl::unordered_dense::map<ReadType, bool> isErroneousReadCached;

				bool isErroneous = isErroneousRead(ii, contigData._alignments, isErroneousReadCached, contigData._coverage, tbuf);
				pair<ReadType, ReadType> readPair = {alignmentData._readIndex1, alignmentData._readIndex2};


				#pragma omp critical(precomputeAlignments)
				{
					contigData._isErroneousReadCached[alignmentData._readIndex1] = isErroneous;
					contigData._cachedAlignments[readPair] = allAlignments;
				}
			}

			mm_tbuf_destroy(tbuf);
		}

		vector<PrecomputeAlignmentData>().swap(_alignmentsToPrecompute); //_alignmentsToPrecompute.clear();
		

	}

	/*
	unordered_map<u_int32_t, vector<u_int32_t>> _contigIndex_to_breakpoints;

	void loadContigBreakpoints(const string& filename){

		ifstream inputFile(filename);

		while(true){

			u_int32_t contigIndex;
			u_int32_t position;

			inputFile.read((char*)&contigIndex, sizeof(contigIndex));

			if(inputFile.eof()) break;

			inputFile.read((char*)&position, sizeof(position));

			_contigIndex_to_breakpoints[contigIndex].push_back(position);

			cout << "Load breakpoint: " << contigIndex << " " << position << endl;
		}

		inputFile.close();

	}
	*/

	void searchDebugContig(const string& contigFilename){
		
		_foundDebugContig = false;

		KminmerParserParallel parser(contigFilename, _params._minimizerSize, _params._kminmerSize, false, false, 1);
		parser._hasContigIndex = true;
		parser.parseSequences(SearchDebugContigFunctor(*this));
	}

	class SearchDebugContigFunctor {

		public:

		ToBasespace2& _toBasespace;

		SearchDebugContigFunctor(ToBasespace2& toBasespace) : _toBasespace(toBasespace){
		}

		SearchDebugContigFunctor(const SearchDebugContigFunctor& copy) : _toBasespace(copy._toBasespace){
		}

		~SearchDebugContigFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {
			ReadType contigIndex = kminmerList._readIndex;
			if(_toBasespace._debug_contigIndex == contigIndex){
				_toBasespace._foundDebugContig = true;
			}
		}
	};

	struct ReadCount{
		ReadType _readIndex;
		u_int32_t _count;
	};

	struct ReadMapping{
		ReadType _readIndex;
		bool _isReversed;
	};

	struct Overlap{
		ReadType _readIndex;
		u_int32_t _overlapLength;
	};

	struct NextAligment{
		ReadMapping2 _mapping;
		size_t _j;
	};

	struct ReadPath{
		vector<ReadType> _readPath;
		int64_t _contigStart;
		int64_t _contigEnd;

		//int64_t getSize(){
		//	return _contigEnd - _contigStart;
		//}
	};

	struct PrecomputeAlignmentData{
		u_int32_t _contigIndex;
		ReadType _readIndex1;
		ReadType _readIndex2;
	};
	//phmap::parallel_flat_hash_map<pair<ReadType, ReadType>, vector<AlignmentBounds>> _readsVsReadsAlignmentsCached;
	vector<PrecomputeAlignmentData> _alignmentsToPrecompute;

	class AlignmentCacheFunctor {

		public:

		ToBasespace2& _toBasespace;
		mm_tbuf_t* _tbuf;

		AlignmentCacheFunctor(ToBasespace2& toBasespace) : _toBasespace(toBasespace){
			_tbuf = mm_tbuf_init();
		}

		AlignmentCacheFunctor(const AlignmentCacheFunctor& copy) : _toBasespace(copy._toBasespace){
			_tbuf = mm_tbuf_init();
		}

		~AlignmentCacheFunctor(){
			mm_tbuf_destroy(_tbuf);
		}


		void operator () (const KminmerList& kminmerList) {

			ReadType contigIndex = kminmerList._readIndex;

			if(_toBasespace._contigIndex_to_alignments.find(contigIndex) == _toBasespace._contigIndex_to_alignments.end()) return;

			vector<ReadMapping2> alignments;
			ankerl::unordered_dense::map<ReadType, u_int32_t> readIndex_to_i;
			phmap::parallel_flat_hash_map<pair<ReadType, ReadType>, vector<AlignmentBounds>> readsVsReadsAlignmentsCached;
			phmap::parallel_flat_hash_map<pair<ReadType, ReadType>, AlignmentBounds> readsVsReadsAlignmentsUsed;

			ReadType startingReadIndex = -1;
			u_int32_t minPosition = -1;

			unordered_map<ReadType, ReadMapping2> readIndex_to_alignment;

			vector<u_int32_t> coverages(kminmerList._readMinimizers.size(), 0);
			for(const ReadMapping2& al : _toBasespace._contigIndex_to_alignments[contigIndex]._alignments){

				for(size_t i=al._contigStart; i<al._contigEnd; i++){
					coverages[i] += 1;
				}
				
				alignments.push_back(al);

				readIndex_to_alignment[al._readIndex] = al;
			}

			long double sum = 0;
			long double n = 0;
			for(size_t i=0; i<coverages.size(); i++){
				sum += coverages[i];
				n += 1;
				//cout << i << "\t" << coverages[i] << endl;
				//if(i % 5000 == 0 && i != 0) getchar();
			}
			
			float contigCoverage = sum / n;
			
			if(_toBasespace._print_debug_main){
				#pragma omp critical(writeContig)
				{
					cout << "Start contig:\t" << contigIndex << "\t" << kminmerList._readMinimizers.size() << "\t" << ((u_int32_t)contigCoverage) << "x\t" << alignments.size() << endl;
				}
			}
			
			#pragma omp critical(writeContig)
			{
				_toBasespace._contigIndex_to_alignments[contigIndex]._coverage = contigCoverage;
			}
			

			if(_toBasespace._print_debug) cout << "Nb alignments: " << alignments.size() << endl;

			if(contigCoverage <= 1) return;
			if(alignments.size() == 0) return;

			std::sort(alignments.begin(), alignments.end(), [](const ReadMapping2& a, const ReadMapping2& b){

				if(a._contigStart == b._contigStart){
					if(a._contigEnd == b._contigEnd){
						return a._readIndex < b._readIndex;
					}

					return a._contigEnd < b._contigEnd;
				}

				return a._contigStart < b._contigStart;
			});

			for(size_t i=0; i<alignments.size(); i++){
				const ReadMapping2& alignment = alignments[i];
				readIndex_to_i[alignment._readIndex] = i;
			}

			#pragma omp critical(writeContig)
			{
				_toBasespace._contigIndex_to_alignments[contigIndex]._alignments = alignments;
				_toBasespace._contigIndex_to_alignments[contigIndex]._readIndex_to_i = readIndex_to_i;
			}


			vector<ReadPath> readPaths;

			while(true){

				bool foundStart = getPath(kminmerList, readPaths, alignments, readIndex_to_alignment, readIndex_to_i, readsVsReadsAlignmentsCached, readsVsReadsAlignmentsUsed, contigCoverage);
				//getchar();
				if(!foundStart) break;
			}

			#pragma omp critical(AlignmentCacheFunctor)
			{

				//for(long i=0; i<((long)alignments.size())-1; i++){
				//	pair<ReadType, ReadType> readPair = {alignments[i]._readIndex, alignments[i+1]._readIndex};
				//	_toBasespace._readsVsReadsAlignmentsCached[readPair] = {};
				//}
				
				
				for(const ReadPath& readPath : readPaths){
					for(long i=0; i<((long)readPath._readPath.size())-1; i++){
						//cout << "Collect:\t" << readPath._readPath[i] << "\t" << readPath._readPath[i+1] << endl;
						_toBasespace._alignmentsToPrecompute.push_back({contigIndex, readPath._readPath[i], readPath._readPath[i+1]});
						//pair<ReadType, ReadType> readPair = {readPath._readPath[i], readPath._readPath[i+1]};
						//_toBasespace._readsVsReadsAlignmentsCached[readPair] = {};
					}
					//cout << readPath._readPath.size() << endl;
				}
				
			}

		}


		bool getPath(const KminmerList& kminmerList, vector<ReadPath>& readPaths, const vector<ReadMapping2>& alignments, auto& readIndex_to_alignment, auto& readIndex_to_i, auto& readsVsReadsAlignmentsCached, auto& readsVsReadsAlignmentsUsed, const float contigCoverage){

			//bool isAggressive = false;
			//u_int32_t maxAggressiveContigEnd = 0;
			//cout << isErroneousRead(readIndex_to_i[2221053], alignments, contigCoverage) << endl;
			//exit(1);

			ReadPath currentReadPath;
			//unordered_set<ReadType> excludedReadIndexes;
			//unordered_map<ReadType, bool> isErroneousReadCached;
			
			int64_t startI = 0;
			ReadMapping2 bestStartingAlignment; // = alignments[0];
			int64_t maxScore = std::numeric_limits<int64_t>::min();
			
			bool foundStart = false;

			u_int64_t minContigStart = -1;

			for(size_t i=0; i<alignments.size(); i++){

				const ReadMapping2& alignment = alignments[i];
				//cout << i << "/" << alignments.size() << "\t" << alignment._contigStart << "\t" << alignment._contigEnd << "\t" << alignmentOverlapExistingReadPath(alignment, readPaths) << endl;
				if(_toBasespace.alignmentOverlapExistingReadPath(alignment, readPaths)) continue;
				//if(alignment._isReversed) continue;

				//if(isErroneousRead(readIndex_to_i[alignment._readIndex], alignments, isErroneousReadCached, contigCoverage)){
					//cout << "\tErroneous" << endl;
					//excludedReadIndexes.insert(readIndex2);
				//	continue;
				//}

				if(minContigStart == -1){
					minContigStart = alignment._contigStart;
				}

				//if(_toBasespace._print_debug) cout << alignment._readIndex << "\t" <<  alignment._contigStart << "\t" << alignment._contigEnd << "\t" << alignment._matchScore << endl;

				if(alignment._contigStart > minContigStart) break;

				if(alignment._matchScore > maxScore){
					bestStartingAlignment = alignment;
					maxScore = alignment._matchScore;
					startI = i;
					foundStart = true;
				}

			}

			if(!foundStart) return false;
			//cout << "d" << endl;
			//cout << bestStartingAlignment._readStartReal << " " << bestStartingAlignment._readEndReal << endl;

			//if(_toBasespace._print_debug) cout << bestStartingAlignment._readIndex << " " << bestStartingAlignment._isReversed << endl; //<< " " << (_toBasespace._mReads.find(bestStartingAlignment._readIndex) != _toBasespace._mReads.end()) << endl;

			//char* dnaStr1 = _toBasespace._mReads[bestStartingAlignment._readIndex]->to_string();
			//string read1 = string(dnaStr1, strlen(dnaStr1));
			//free(dnaStr1);

			//if(bestStartingAlignment._isReversed){
			//	Utils::toReverseComplement(read1);
			//}

			//u_int64_t nbFailed = 0;
			//u_int64_t failedContigEnd = 0;

			ReadMapping2 lastAlignment;
			//string lastRead;

			vector<ReadType> readPath;
			readPath.push_back(bestStartingAlignment._readIndex);
			//contigSequence += read1.substr(bestStartingAlignment._readStartReal, read1.size());

			//cout << "e" << endl;
			//if(_toBasespace._print_debug) cout << "todo startI a choisir en boucle aussi" << endl;

			for(size_t i=startI; i<alignments.size(); i++){

				//cout << "louloul: " << isErroneousRead(readIndex_to_i[2221053], alignments, contigCoverage) << endl;
				//if(isErroneousRead(readIndex_to_i[2221053], alignments, contigCoverage)) getchar();

				const ReadMapping2& alignment1 = alignments[i];
				//cout << i << endl;

				const ReadType readIndex1 = readPath[readPath.size()-1];
				//if(excludedReadIndexes.find(readIndex1) != excludedReadIndexes.end()) continue;

				//char* dnaStr1 = _toBasespace._mReads[readIndex1]->to_string();
				//string read1 = string(dnaStr1, strlen(dnaStr1));
				//free(dnaStr1);
				//if(alignment1._isReversed){
				//	Utils::toReverseComplement(read1);
				//}
				//mm_idx_t* mi = minimap2_indexRead(read1);


				//if(_toBasespace._print_debug) cout << "Start: " << i << "\t" << _toBasespace._debug_readHeaders[alignment1._readIndex] << "\t" << alignment1._readIndex << "\t" << alignment1._contigStart << "\t" << alignment1._contigEnd << "\t" << alignment1._isReversed << endl;
				//int nbAlignments = 0;

				//if(alignment1._contigStart > 1000) break;

				//u_int64_t maxAlignLength = 0;
				//bool foundSuccessor = false;

				ReadMapping2 lastAlignment = getBestSuccessor(alignments, i, readIndex1, alignment1, readIndex_to_i);
				
				//cout << "Lul:\t" << readIndex1 << "\t" << lastAlignment._readIndex << endl;

				//if(lastAlignment._contigIndex == -1){
				//	lastAlignment = getBestSuccessor(alignments, i, read1, mi, readIndex1, alignment1, contigCoverage, readsVsReadsAlignmentsCached, readsVsReadsAlignmentsUsed, readIndex_to_i, excludedReadIndexes, isErroneousReadCached, false, isAggressive);
					
					//if(lastAlignment._contigIndex == -1){
					//	lastAlignment = getBestSuccessor(alignments, i, read1, mi, readIndex1, alignment1, contigCoverage, readsVsReadsAlignmentsCached, readsVsReadsAlignmentsUsed, readIndex_to_i, excludedReadIndexes, isErroneousReadCached, true, true);
						
					//	if(lastAlignment._contigIndex == -1){
					//		lastAlignment = getBestSuccessor(alignments, i, read1, mi, readIndex1, alignment1, contigCoverage, readsVsReadsAlignmentsCached, readsVsReadsAlignmentsUsed, readIndex_to_i, excludedReadIndexes, isErroneousReadCached, false, true);
					//	}
					//}
					
				//}
				
				bool foundSuccessor = lastAlignment._contigIndex != -1;

				//mm_idx_destroy(mi);

				if(foundSuccessor){
					
					
					//if(lastAlignment._contigEnd > maxAggressiveContigEnd){
					//	isAggressive = false;
					//}
					
					//processedReadIndexes.insert(lastAlignment._readIndex);
					readPath.push_back(lastAlignment._readIndex);
					//if(_toBasespace._print_debug) cout << "\tBest successor: " << lastAlignment._readIndex << " " <<  lastAlignment._contigStart*200 << " " << lastAlignment._contigEnd*200 << endl;
				}

				//u_int32_t prevOverhang =  bestAlignment._referenceLength - bestAlignment._referenceEnd;
				//if(prevOverhang > 0) contigSequence = contigSequence.substr(0, contigSequence.size()-prevOverhang);

				//cout << prevOverhang << endl;
				if(!foundSuccessor){
					
					if(readPath.size() > 0){
						if(currentReadPath._readPath.size() == 0){ //new path
							const ReadMapping2& alignmentStart = readIndex_to_alignment[readPath[0]];
							const ReadMapping2& alignmentEnd = readIndex_to_alignment[readPath[readPath.size()-1]];
							//cout << "\tAdd read path (new):\t" <<  readPath.size() << "\t" << alignmentStart._contigStart << "\t" << alignmentEnd._contigEnd << endl;
							currentReadPath = {readPath, alignmentStart._contigStart, alignmentEnd._contigEnd};
							//getchar();
						}
						else{
							const ReadMapping2& alignmentStart = readIndex_to_alignment[readPath[0]];
							const ReadMapping2& alignmentEnd = readIndex_to_alignment[readPath[readPath.size()-1]];

							if(alignmentEnd._contigEnd > currentReadPath._contigEnd){
								//cout << "\tAdd read path (existing):\t" <<  readPath.size() << "\t" << alignmentStart._contigStart << "\t" << alignmentEnd._contigEnd << endl;
								currentReadPath = {readPath, alignmentStart._contigStart, alignmentEnd._contigEnd};
							}
						}
					}

					break;
					//if(lastAlignment._contigEnd+10 > kminmerList._readMinimizers.size()) break;
					
					//if(_toBasespace._print_debug) cout << lastAlignment._contigEnd << " " << kminmerList._readMinimizers.size() << endl;

					//ofstream readFile("/pasteur/appa/homes/gbenoit/zeus/tools/merging/metaMDBG/build/read.fasta");
					//readFile << ">read_" << lastAlignment._readIndex << endl;
					//readFile << lastRead << endl;
					//readFile.close();

					//if(_toBasespace._print_debug) cout << "derp: " << alignment1._contigEnd << endl;

					//excludedReadIndexes.insert(readPath[readPath.size()-1]);
					//readPath.pop_back();
					
					//if(alignment1._contigEnd > failedContigEnd){
					//	nbFailed = 0;
					//	failedContigEnd = alignment1._contigEnd;
					//}

					//if(readPath.size() == 0) break; //No successor after starting read

					/*
					//nbFailed += 1;
					//failedContigEnd = alignment1._contigEnd;

					if(nbFailed > 10){
						//cout << "Failed to find successor: " << kminmerList._readIndex << " " << kminmerList._readMinimizers.size() << endl;
						if(isAggressive){
							break;
						}
						else{

							if(maxAggressiveContigEnd == currentReadPath._contigEnd) break;

							//cout << "Start agressive" << endl;
							//cout << currentReadPath._readPath.size() << " " << currentReadPath._contigEnd << endl;
							isAggressive = true;
							maxAggressiveContigEnd = currentReadPath._contigEnd;
							nbFailed = 0;
							excludedReadIndexes.clear();
							readPath = currentReadPath._readPath;

							//getchar();

							
						}
					}

					//getchar();
					//getchar();
					//continue;
					//break;
					*/
				}





				//string readSeqEnd = bestRead.substr(bestAlignment._queryEnd, bestRead.size());
				//contigSequence += readSeqEnd;



				//cout << "\tNb aligmnents: " << nbAlignments << endl;
				//cout << "\tBest successor: " << lastAlignment._readIndex << endl;
				//cout << "\tContig size: " << contigSequence.size() << "  (Added: " << readSeqEnd.size() << ")" << endl;

				i = readIndex_to_i[readPath[readPath.size()-1]]-1;

				//getchar();
			}

			if(currentReadPath._readPath.size() > 0){
				//cout << "\tAdd read path: " << currentReadPath._contigStart << " " << currentReadPath._contigEnd << endl;
				//getchar();
				readPaths.push_back(currentReadPath);
			}

			return true;
		}


		ReadMapping2 getBestSuccessor(const vector<ReadMapping2>& alignments, const size_t i, const ReadType readIndex1, const ReadMapping2& alignment1, auto& readIndex_to_i){

			//cout << "getBestSuccessor:\t" << readIndex1 << "\t" << i << endl;

			vector<ReadMapping2> nextAlignments;
			for(size_t j=i+1; j<alignments.size(); j++){
				
				const ReadMapping2& alignment2 = alignments[j];
				
				//cout << "\t" << j << "\t" << alignment2._readIndex << endl;	

				//cout << alignment1._contigStart << " " << alignment1._contigEnd << "\t" << alignment2._contigStart << " " << alignment2._contigEnd << endl;
				//if(excludedReadIndexes.find(alignment2._readIndex) != excludedReadIndexes.end()) continue;
				//if(_toBasespace._print_debug) cout << "\tPossible:\t" << alignment2._readIndex << "\t" << _toBasespace._debug_readHeaders[alignment2._readIndex] << "\t" << alignment2._contigStart << "\t" << alignment2._contigEnd << "\t" << alignment2._isReversed << endl;

				//if(overlapOnTheReferenceOnly){
				if(alignment2._contigStart > alignment1._contigEnd) break;
				if(!_toBasespace.overlapOnTheReference(alignment1, alignment2)) continue;
				//}
				//else{
				//	if(alignment2._contigStart == alignment1._contigStart) continue;
				//	if(alignment2._contigEnd < alignment1._contigEnd) continue;

				//	if(alignment2._contigStart > alignment1._contigEnd+100) break;
				//}

				nextAlignments.push_back(alignment2);
			}

			//if(readIndex1 == 4535) getchar();
			std::sort(nextAlignments.begin(), nextAlignments.end(), [](const ReadMapping2& a, const ReadMapping2& b){
				if(a._matchScore == b._matchScore){
					return a._readIndex < b._readIndex;
				}
				return a._matchScore > b._matchScore;
			});


			//for(size_t j=i+1; j<alignments.size(); j++){
			for(const ReadMapping2& alignment2: nextAlignments){
				
				//const ReadMapping2& alignment2 = alignments[j];

				//cout << "\tTrying:\t" << alignment2._readIndex << "\t" << alignment2._contigStart << "\t" << alignment2._contigEnd << endl;

				//if(alignment2._contigStart > alignment1._contigEnd) break;
				//if(!overlapOnTheReference(alignment1, alignment2)) continue;

				//const ReadMapping2& alignment2 = alignment22._mapping;

				//const ReadType readIndex2 = alignment2._readIndex;
				//char* dnaStr2 = _toBasespace._mReads[readIndex2]->to_string();
				//string read2 = string(dnaStr2, strlen(dnaStr2));
				//free(dnaStr2);

				//if(alignment2._isReversed){
				//	Utils::toReverseComplement(read2);
				//}

				//int64_t maxContigSeq = min((int64_t) contigSequence.size(), (int64_t)(read2.size()+5000));
				//string contigSeq = contigSequence.substr(contigSequence.size()-maxContigSeq, maxContigSeq);

				//const pair<ReadType, ReadType> readPair = {readIndex1, readIndex2};
				//vector<AlignmentBounds> allAlignments;
				//AlignmentBounds alignment;
				
				/*
				if(readsVsReadsAlignmentsCached.find(readPair) == readsVsReadsAlignmentsCached.end()){
					computeAlignment(read1, read2, _tbuf, _toBasespace._minimap2Preset_ava, true, ToBasespace2::_minOverlap, mi, allAlignments);
					readsVsReadsAlignmentsCached[readPair] = allAlignments;
				}
				else{
					allAlignments = readsVsReadsAlignmentsCached[readPair];
				}

				AlignmentBounds alignment;
				u_int32_t minLength = 0;

				for(const AlignmentBounds& al : allAlignments){

					if(allowErroneousReads){
						if(!isValidOverlapAlignment(al, false)) continue;
					}
					else{
						if(!isValidOverlapAlignment(al, true)) continue;
					}

					u_int32_t alignLength = min(al._queryEnd - al._queryStart, al._referenceEnd - al._referenceStart);
					
					if(alignLength > minLength){
						minLength = alignLength;
						alignment = al;
					}
				}
				

				//if(_toBasespace._print_debug) cout << "\t" << _toBasespace._debug_readHeaders[alignment2._readIndex] << "\t" << readIndex1 << " -> " << readIndex2 << " " << alignment._queryStart << " " << alignment2._isReversed << endl;

				//cout << read1 << endl;
				//cout << read2 << endl;
				//getchar();
				//nbAlignments += 1;

				//cout << "\t\t" << j << " " << alignment2._readIndex << "     " << alignment._queryStart << endl;
			

				if(alignment._queryStart == -1){
					//cout << "\tNo valid overlap" << endl;
					continue;
				}


				//if(readIndex2 == 866271){
				//	cout << "haaaaaaaaaaaa" << endl;
				//	cout << isErroneousRead(readIndex_to_i[readIndex2], alignments) << endl;
				//	getchar();
				//}

				if(!allowErroneousReads && isErroneousRead(readIndex_to_i[readIndex2], alignments, isErroneousReadCached, contigCoverage)){
					//cout << "\tErroneous" << endl;
					//if(readIndex2 == 838925 || readIndex2 == 864136 || readIndex2 == 2221053){
					//	cout << "derp444: "<< readIndex2 << endl;
					//	getchar();
					//}
					//excludedReadIndexes.insert(readIndex2);
					continue;
				}
				
				//if(readIndex1 == 108086){
				//	cout << "Overlap on the ref only: " << overlapOnTheReferenceOnly << endl;
				//	cout << "Allow error: " << allowErroneousReads << endl;
				//	cout << readIndex2 << endl;
				//	for(const auto& it : allAlignments){
				//		cout << it._queryStart << "-" << it._queryEnd << "-" << it._isReversed << "   " << it._referenceStart << "-" << it._referenceEnd<< endl;
				//	}
				//	cout << "Used: " << alignment._referenceStart << "-" << alignment._referenceEnd << "-" << alignment._isReversed << endl;
				//}
				
				//cout << readIndex1 << " -> " << readIndex2 << " " << alignment._referenceStart << "-" << alignment._referenceEnd << endl;
				//getchar();

				//if(alignment._queryEnd - alignment._queryStart > maxAlignLength){
				//maxAlignLength = alignment._queryEnd - alignment._queryStart;
				//bestAlignment = alignment;
				//bestRead = read2;
				//i2 = alignment22._j;
				//lastAlignment = alignment2;
				//lastRead = read2;
				//foundSuccessor = true;

				//}
				readsVsReadsAlignmentsUsed[readPair] = alignment;
				*/

				//cout << alignment2._readIndex << endl;
				return alignment2;
				//cout << "\t" << j << " " << alignment2._readIndex << "     " << alignment1._contigStart << "-" << alignment1._contigEnd << "\t" << alignment2._contigStart << "-" << alignment2._contigEnd << endl;
			

			}

			//cout << "null" << endl;
			ReadMapping2 nullMapping;
			nullMapping._readIndex = -1;
			nullMapping._contigIndex = -1;
			return nullMapping;
		}


	};
	
	
	class CreateBaseContigsFunctor {

		public:

		ToBasespace2& _toBasespace;
		mm_tbuf_t* _tbuf;
		KmerModel* _modelRepetitive;
		KmerModel* _modelTri;


		struct KmerAbundance{
			u_int64_t _kmer;
			u_int64_t _abundance;
		};

		CreateBaseContigsFunctor(ToBasespace2& toBasespace) : _toBasespace(toBasespace){
			_tbuf = mm_tbuf_init();
			_modelRepetitive = new KmerModel(21);
			_modelTri = new KmerModel(3);
		}

		CreateBaseContigsFunctor(const CreateBaseContigsFunctor& copy) : _toBasespace(copy._toBasespace){
			_tbuf = mm_tbuf_init();
			_modelRepetitive = new KmerModel(21);
			_modelTri = new KmerModel(3);
		}

		~CreateBaseContigsFunctor(){
			mm_tbuf_destroy(_tbuf);
			delete _modelRepetitive;
			delete _modelTri;
		}

		void operator () (const KminmerList& kminmerList) {


			//if(kminmerList._readMinimizers.size() > 1000){
			//	cout << "SKIP" << endl;
			//	return;
			//}
			//#pragma atomic
			//_toBasespace.__nbLala += 1;

			//if(_toBasespace.__nbLala > 5000) return;
			//cout << "Start contig: " << kminmerList._readIndex << " " << kminmerList._readMinimizers.size() << endl;
			//bool isCircular = kminmerList._isCircular;

		
			ReadType contigIndex = kminmerList._readIndex;
			
			if(_toBasespace._debug_contigIndex != -1 && contigIndex != _toBasespace._debug_contigIndex) return;
			//cout << "Start contig: " << contigIndex << " " << kminmerList._readMinimizers.size() << endl;
			//if(contigIndex != 24817) return; //2851 2873
			//if(kminmerList._readMinimizers.size() < 4000) return;



			if(_toBasespace._contigIndex_to_alignments.find(contigIndex) == _toBasespace._contigIndex_to_alignments.end()) return;

			//cout << contigIndex << " " << kminmerList._readMinimizers.size() << " " << _toBasespace._contigIndex_to_alignments[contigIndex].size() << endl;
			//MinimizerReadPacked read = {kminmerList._readIndex, kminmerList._readMinimizers, kminmerList._minimizerPos, kminmerList._meanReadQuality, kminmerList._readLength};
			//MinimizerReadPacked readLowDensity = Utils::getLowDensityMinimizerRead(read, _parent._minimizerDensity_assembly);
			
			const ContigData& contigData = _toBasespace._contigIndex_to_alignments[contigIndex];

			vector<ReadMapping2> alignments = contigData._alignments;
			ankerl::unordered_dense::map<ReadType, u_int32_t> readIndex_to_i;
			//unordered_map<ReadType, ReadMapping2> readIndex_to_alignment;
			phmap::parallel_flat_hash_map<pair<ReadType, ReadType>, vector<AlignmentBounds>> readsVsReadsAlignmentsCached;
			phmap::parallel_flat_hash_map<pair<ReadType, ReadType>, AlignmentBounds> readsVsReadsAlignmentsUsed;
			ankerl::unordered_dense::map<ReadType, bool> isErroneousReadCached;

			ReadType startingReadIndex = -1;
			u_int32_t minPosition = -1;

			unordered_map<ReadType, ReadMapping2> readIndex_to_alignment;

			u_int64_t maxContigEnd = 0;

			//vector<u_int32_t> coverages(kminmerList._readMinimizers.size(), 0);
			for(const ReadMapping2& al : alignments) {//_toBasespace._contigIndex_to_alignments[contigIndex]._alignments){

				//for(size_t i=al._contigStart; i<al._contigEnd; i++){
				//	coverages[i] += 1;
				//}
				//if(al._readIndex == 108086){
				//	cout << "skip" << endl;
				//	continue;
				//}

				/*

				if(al._readIndex == 1188374){
					cout << "skip" << endl;
					continue;
				}

				if(al._readIndex == 1263006){
					cout << "skip" << endl;
					continue;
				}
				*/
				
				//alignments.push_back(al);

				readIndex_to_alignment[al._readIndex] = al;

				if(al._contigEnd > maxContigEnd){
					maxContigEnd = al._contigEnd;
				}
			}

			/*
			long double sum = 0;
			long double n = 0;
			for(size_t i=0; i<coverages.size(); i++){
				sum += coverages[i];
				n += 1;
				//cout << i << "\t" << coverages[i] << endl;
				//if(i % 5000 == 0 && i != 0) getchar();
			}
			
			float contigCoverage = sum / n;
			*/
			float contigCoverage = contigData._coverage;

			if(_toBasespace._print_debug_main){
				#pragma omp critical(writeContig)
				{
					cout << "Start contig:\t" << contigIndex << "\t" << kminmerList._readMinimizers.size() << "\t" << ((u_int32_t)contigCoverage) << "x\t" << alignments.size() << endl;
				}
			}

			if(_toBasespace._print_debug) cout << "Nb alignments: " << alignments.size() << endl;

			if(contigCoverage <= 1) return;
			if(alignments.size() == 0) return;
			/*
			std::sort(alignments.begin(), alignments.end(), [](const ReadMapping2& a, const ReadMapping2& b){

				if(a._contigStart == b._contigStart){
					if(a._contigEnd == b._contigEnd){
						return a._readIndex < b._readIndex;
					}

					return a._contigEnd < b._contigEnd;
				}

				return a._contigStart < b._contigStart;
			});


			for(size_t i=0; i<alignments.size(); i++){
				const ReadMapping2& alignment = alignments[i];
				readIndex_to_i[alignment._readIndex] = i;
				//readIndex_to_alignment[alignment._readIndex] = alignment;
				
				//cout << alignment._contigStart << "\t" << alignment._contigEnd << endl;

			}
			*/

			readIndex_to_i = contigData._readIndex_to_i;

			for(const auto& it : contigData._isErroneousReadCached){
				const ReadType readIndex = it.first;
				const bool isErroneous = it.second;
				isErroneousReadCached[readIndex] = isErroneous;
			}

			for(const auto& it : contigData._cachedAlignments){
				const pair<ReadType, ReadType>& readPair = it.first;
				const vector<AlignmentBounds>& alignment = it.second;
				readsVsReadsAlignmentsCached[readPair] = alignment;
			}

			//cout << isErroneousReadCached.size() << " " << readsVsReadsAlignmentsCached.size() << endl;
			/*
			if(_toBasespace._print_debug){
				ofstream allReadFile("/pasteur/appa/homes/gbenoit/zeus/tools/merging/metaMDBG/build/allReads.fasta");

				for(size_t i=0; i<alignments.size(); i++){
					const ReadMapping2& alignment = alignments[i];

					char* dnaStr1 = _toBasespace._mReads[alignment._readIndex]->to_string();
					string read1 = string(dnaStr1, strlen(dnaStr1));
					free(dnaStr1);
					if(alignment._isReversed){
						Utils::toReverseComplement(read1);
					}

					allReadFile << ">read_" << alignment._readIndex << endl;
					allReadFile << read1 << endl;

				}

				allReadFile.close();
			}
			*/

			vector<ReadPath> readPaths;

			while(true){

				bool foundStart = getPath(kminmerList, readPaths, alignments, readIndex_to_alignment, readIndex_to_i, readsVsReadsAlignmentsCached, readsVsReadsAlignmentsUsed, isErroneousReadCached, contigCoverage, maxContigEnd);

				//cout << foundStart << " " << readPaths[0]._contigStart << " " <<  readPaths[0]._contigEnd << endl;
				//getchar();
				if(!foundStart) break;
			}


			u_int64_t nbPrecomputedAlignments = 0;
			u_int64_t nbPrecomputedErroneous = 0;

			for(const ReadPath& readPath : readPaths){
				for(long i=0; i<((long)readPath._readPath.size())-1; i++){
					pair<ReadType, ReadType> readPair = {readPath._readPath[i], readPath._readPath[i+1]};
					if(contigData._cachedAlignments.find(readPair) != contigData._cachedAlignments.end()){
						nbPrecomputedAlignments += 1;
						//#pragma omp atomic
						//_toBasespace._nbPrecomputedAlignments += 1;
					}


					if(contigData._isErroneousReadCached.find(readPath._readPath[i]) != contigData._isErroneousReadCached.end()){
						nbPrecomputedErroneous += 1;
						//#pragma omp atomic
						//_toBasespace._nbPrecomputedAlignments += 1;
					}

				}
				//cout << readPath._readPath.size() << endl;


			}

			/*
			#pragma omp critical(AlignmentCacheFunctor)
			{
				cout << endl;
				cout << contigIndex << endl;
				cout << "omg: " << nbPrecomputedAlignments << " " << readsVsReadsAlignmentsCached.size() << endl;
				cout << "omg: " << nbPrecomputedErroneous << " " << isErroneousReadCached.size() << endl;
				//getchar();

				//if(contigIndex == 19476){
				//	getchar();
				//}
			}
			*/
			/*
			#pragma omp critical(AlignmentCacheFunctor)
			{

				//for(long i=0; i<((long)alignments.size())-1; i++){
				//	pair<ReadType, ReadType> readPair = {alignments[i]._readIndex, alignments[i+1]._readIndex};
				//	_toBasespace._readsVsReadsAlignmentsCached[readPair] = {};
				//}
				
				
				for(const ReadPath& readPath : readPaths){
					for(long i=0; i<((long)readPath._readPath.size())-1; i++){
						pair<ReadType, ReadType> readPair = {readPath._readPath[i], readPath._readPath[i+1]};
						if(_toBasespace._readsVsReadsAlignmentsCached.find(readPair) != _toBasespace._readsVsReadsAlignmentsCached.end()){
							_toBasespace._nbPrecomputedAlignments += 1;
						}
					}
					//cout << readPath._readPath.size() << endl;
				}

				cout << "omg: " << _toBasespace._nbPrecomputedAlignments << " " << _toBasespace._readsVsReadsAlignmentsCached.size() << endl;
				
			}

			//cout << "Nb read path: " << readPaths.size() << endl;
			//exit(1);
			*/

			readPathsToContigs(kminmerList, readPaths, readsVsReadsAlignmentsUsed, readIndex_to_alignment, contigCoverage);

			if(_toBasespace._print_debug_main){

				#pragma omp critical(writeContig)
				{
					int64_t coveredMinimizers = 0;

					for(size_t i=0; i<readPaths.size(); i++){

						const ReadPath& readPath = readPaths[i];
						coveredMinimizers += (readPath._contigEnd - readPath._contigStart);
					}

					float coveredFraction = ((double) coveredMinimizers) / kminmerList._readMinimizers.size();
					cout << "\tCovered fraction: " << coveredFraction << endl;

					if(coveredFraction < 0.5){
						cout << coveredFraction << endl;
						//getchar();
					}

				}
			}


			//cout << "Nb alignments performed: " << readsVsReadsAlignmentsCached.size() << " " << readsVsReadsAlignmentsUsed.size() << endl;
		}




		bool getPath(const KminmerList& kminmerList, vector<ReadPath>& readPaths, const vector<ReadMapping2>& alignments, auto& readIndex_to_alignment, auto& readIndex_to_i, auto& readsVsReadsAlignmentsCached, auto& readsVsReadsAlignmentsUsed, auto& isErroneousReadCached, const float contigCoverage, const u_int64_t maxContigEnd){


			bool isAggressive = false;
			u_int32_t maxAggressiveContigEnd = 0;
			//cout << isErroneousRead(readIndex_to_i[2221053], alignments, contigCoverage) << endl;
			//exit(1);

			ReadPath currentReadPath;
			unordered_set<ReadType> excludedReadIndexes;
			
			int64_t startI = 0;
			ReadMapping2 bestStartingAlignment; // = alignments[0];
			int64_t maxScore = std::numeric_limits<int64_t>::min();
			
			bool foundStart = false;

			u_int64_t minContigStart = -1;

			for(size_t i=0; i<alignments.size(); i++){

				const ReadMapping2& alignment = alignments[i];
				//cout << i << "/" << alignments.size() << "\t" << alignment._contigStart << "\t" << alignment._contigEnd << "\t" << alignmentOverlapExistingReadPath(alignment, readPaths) << endl;
				if(_toBasespace.alignmentOverlapExistingReadPath(alignment, readPaths)) continue;
				//if(alignment._isReversed) continue;

				if(_toBasespace.isErroneousRead(readIndex_to_i[alignment._readIndex], alignments, isErroneousReadCached, contigCoverage, _tbuf)){
					//cout << "\tErroneous" << endl;
					//excludedReadIndexes.insert(readIndex2);
					continue;
				}

				if(minContigStart == -1){
					minContigStart = alignment._contigStart;
				}

				if(_toBasespace._print_debug) cout << alignment._readIndex << "\t" <<  alignment._contigStart << "\t" << alignment._contigEnd << "\t" << alignment._matchScore << endl;

				if(alignment._contigStart > minContigStart) break;

				if(alignment._matchScore > maxScore){
					bestStartingAlignment = alignment;
					maxScore = alignment._matchScore;
					startI = i;
					foundStart = true;
				}

			}

			if(!foundStart) return false;
			//cout << "d" << endl;
			//cout << bestStartingAlignment._readStartReal << " " << bestStartingAlignment._readEndReal << endl;

			if(_toBasespace._print_debug) cout << bestStartingAlignment._readIndex << " " << bestStartingAlignment._isReversed << endl; //<< " " << (_toBasespace._mReads.find(bestStartingAlignment._readIndex) != _toBasespace._mReads.end()) << endl;

			//char* dnaStr1 = _toBasespace._mReads[bestStartingAlignment._readIndex]->to_string();
			//string read1 = string(dnaStr1, strlen(dnaStr1));
			//free(dnaStr1);

			//if(bestStartingAlignment._isReversed){
			//	Utils::toReverseComplement(read1);
			//}

			u_int64_t nbFailed = 0;
			u_int64_t failedContigEnd = 0;

			ReadMapping2 lastAlignment;
			//string lastRead;

			vector<ReadType> readPath;
			readPath.push_back(bestStartingAlignment._readIndex);
			//contigSequence += read1.substr(bestStartingAlignment._readStartReal, read1.size());

			//cout << "e" << endl;
			if(_toBasespace._print_debug) cout << "todo startI a choisir en boucle aussi" << endl;

			for(size_t i=startI; i<alignments.size(); i++){

				//cout << "louloul: " << isErroneousRead(readIndex_to_i[2221053], alignments, contigCoverage) << endl;
				//if(isErroneousRead(readIndex_to_i[2221053], alignments, contigCoverage)) getchar();

				const ReadMapping2& alignment1 = alignments[i];
				//cout << i << endl;

				const ReadType readIndex1 = readPath[readPath.size()-1];
				if(excludedReadIndexes.find(readIndex1) != excludedReadIndexes.end()) continue;

				const string& read1 = _toBasespace._mReads[readIndex1];
				//string read1 = string(dnaStr1, strlen(dnaStr1));
				//free(dnaStr1);
				//if(alignment1._isReversed){
				//	Utils::toReverseComplement(read1);
				//}
				mm_idx_t* mi = _toBasespace.minimap2_indexRead(read1);


				if(_toBasespace._print_debug) cout << "Start: " << i << "\t" << _toBasespace._debug_readHeaders[alignment1._readIndex] << "\t" << alignment1._readIndex << "\t" << alignment1._contigStart << "\t" << alignment1._contigEnd << "\t" << alignment1._isReversed << endl;
				//int nbAlignments = 0;

				//if(alignment1._contigStart > 1000) break;

				//u_int64_t maxAlignLength = 0;
				//bool foundSuccessor = false;
				ReadMapping2 lastAlignment = getBestSuccessor(alignments, i, read1, mi, readIndex1, alignment1, contigCoverage, readsVsReadsAlignmentsCached, readsVsReadsAlignmentsUsed, readIndex_to_i, excludedReadIndexes, isErroneousReadCached, true, isAggressive);
				
				//cout << lastAlignment._readIndex << "\t" <<  lastAlignment._contigEnd << "\t" << maxContigEnd << endl;

				if(lastAlignment._contigIndex == -1){
					lastAlignment = getBestSuccessor(alignments, i, read1, mi, readIndex1, alignment1, contigCoverage, readsVsReadsAlignmentsCached, readsVsReadsAlignmentsUsed, readIndex_to_i, excludedReadIndexes, isErroneousReadCached, false, isAggressive);
					
					//if(lastAlignment._contigIndex == -1){
					//	lastAlignment = getBestSuccessor(alignments, i, read1, mi, readIndex1, alignment1, contigCoverage, readsVsReadsAlignmentsCached, readsVsReadsAlignmentsUsed, readIndex_to_i, excludedReadIndexes, isErroneousReadCached, true, true);
						
					//	if(lastAlignment._contigIndex == -1){
					//		lastAlignment = getBestSuccessor(alignments, i, read1, mi, readIndex1, alignment1, contigCoverage, readsVsReadsAlignmentsCached, readsVsReadsAlignmentsUsed, readIndex_to_i, excludedReadIndexes, isErroneousReadCached, false, true);
					//	}
					//}
					
				}
				
				bool foundSuccessor = lastAlignment._contigIndex != -1;

				mm_idx_destroy(mi);

				if(foundSuccessor){

					if(lastAlignment._contigEnd > maxAggressiveContigEnd){
						isAggressive = false;
					}
					
					//processedReadIndexes.insert(lastAlignment._readIndex);
					readPath.push_back(lastAlignment._readIndex);
					if(_toBasespace._print_debug) cout << "\tBest successor: " << lastAlignment._readIndex << " " <<  lastAlignment._contigStart*200 << " " << lastAlignment._contigEnd*200 << endl;
				}

				if(foundSuccessor && lastAlignment._contigEnd >= maxContigEnd){ //Le contig ne peut pas etre étendu d'avantage
					foundSuccessor = false;
				}
				//u_int32_t prevOverhang =  bestAlignment._referenceLength - bestAlignment._referenceEnd;
				//if(prevOverhang > 0) contigSequence = contigSequence.substr(0, contigSequence.size()-prevOverhang);

				//cout << prevOverhang << endl;
				if(!foundSuccessor){
					
					if(readPath.size() > 0){
						if(currentReadPath._readPath.size() == 0){ //new path
							const ReadMapping2& alignmentStart = readIndex_to_alignment[readPath[0]];
							const ReadMapping2& alignmentEnd = readIndex_to_alignment[readPath[readPath.size()-1]];
							//cout << "\tAdd read path (new):\t" <<  readPath.size() << "\t" << alignmentStart._contigStart << "\t" << alignmentEnd._contigEnd << endl;
							currentReadPath = {readPath, alignmentStart._contigStart, alignmentEnd._contigEnd};
							//getchar();
						}
						else{
							const ReadMapping2& alignmentStart = readIndex_to_alignment[readPath[0]];
							const ReadMapping2& alignmentEnd = readIndex_to_alignment[readPath[readPath.size()-1]];

							if(alignmentEnd._contigEnd > currentReadPath._contigEnd){
								//cout << "\tAdd read path (existing):\t" <<  readPath.size() << "\t" << alignmentStart._contigStart << "\t" << alignmentEnd._contigEnd << endl;
								currentReadPath = {readPath, alignmentStart._contigStart, alignmentEnd._contigEnd};
							}
						}
					}

					if(lastAlignment._contigEnd >= maxContigEnd) break; //Le contig ne peut pas etre étendu d'avantage

					//if(lastAlignment._contigEnd+10 > kminmerList._readMinimizers.size()) break;
					
					//if(_toBasespace._print_debug) cout << lastAlignment._contigEnd << " " << kminmerList._readMinimizers.size() << endl;

					//ofstream readFile("/pasteur/appa/homes/gbenoit/zeus/tools/merging/metaMDBG/build/read.fasta");
					//readFile << ">read_" << lastAlignment._readIndex << endl;
					//readFile << lastRead << endl;
					//readFile.close();

					if(_toBasespace._print_debug) cout << "derp: " << alignment1._contigEnd << endl;

					excludedReadIndexes.insert(readPath[readPath.size()-1]);
					readPath.pop_back();
					
					if(alignment1._contigEnd > failedContigEnd){
						nbFailed = 0;
						failedContigEnd = alignment1._contigEnd;
					}

					if(readPath.size() == 0) break; //No successor after starting read

					nbFailed += 1;
					//failedContigEnd = alignment1._contigEnd;

					if(nbFailed > 10){
						//cout << "Failed to find successor: " << kminmerList._readIndex << " " << kminmerList._readMinimizers.size() << endl;
						if(isAggressive){
							break;
						}
						else{

							if(maxAggressiveContigEnd == currentReadPath._contigEnd) break;

							//cout << "Start agressive" << endl;
							//cout << currentReadPath._readPath.size() << " " << currentReadPath._contigEnd << endl;
							isAggressive = true;
							maxAggressiveContigEnd = currentReadPath._contigEnd;
							nbFailed = 0;
							excludedReadIndexes.clear();
							readPath = currentReadPath._readPath;

							//getchar();

							
						}
					}

					//getchar();
					//getchar();
					//continue;
					//break;
				}





				//string readSeqEnd = bestRead.substr(bestAlignment._queryEnd, bestRead.size());
				//contigSequence += readSeqEnd;



				//cout << "\tNb aligmnents: " << nbAlignments << endl;
				//cout << "\tBest successor: " << lastAlignment._readIndex << endl;
				//cout << "\tContig size: " << contigSequence.size() << "  (Added: " << readSeqEnd.size() << ")" << endl;

				i = readIndex_to_i[readPath[readPath.size()-1]]-1;

				//getchar();
			}

			if(currentReadPath._readPath.size() > 0){
				//cout << "\tAdd read path: " << currentReadPath._contigStart << " " << currentReadPath._contigEnd << endl;
				//getchar();
				readPaths.push_back(currentReadPath);
			}

			return true;
		}

		//"reprise: allow erroneous: i lfaut aussi augmenter le maxHang dans isValidOverlap?"

		ReadMapping2 getBestSuccessor(const vector<ReadMapping2>& alignments, const size_t i, const string& read1, mm_idx_t* mi, const ReadType readIndex1, const ReadMapping2& alignment1, const float contigCoverage, auto& readsVsReadsAlignmentsCached, auto& readsVsReadsAlignmentsUsed, auto& readIndex_to_i, auto& excludedReadIndexes, auto& isErroneousReadCached, const bool overlapOnTheReferenceOnly, const bool allowErroneousReads){

			//AlignmentBounds bestAlignment;
			//string bestRead = ""; 
			//vector<Overlap> overlaps;

			//size_t i2 = i;

			vector<ReadMapping2> nextAlignments;
			for(size_t j=i+1; j<alignments.size(); j++){
				
				const ReadMapping2& alignment2 = alignments[j];
				
				if(excludedReadIndexes.find(alignment2._readIndex) != excludedReadIndexes.end()) continue;
				//if(_toBasespace._print_debug) cout << "\tPossible:\t" << alignment2._readIndex << "\t" << _toBasespace._debug_readHeaders[alignment2._readIndex] << "\t" << alignment2._contigStart << "\t" << alignment2._contigEnd << "\t" << alignment2._isReversed << endl;

				if(overlapOnTheReferenceOnly){
					if(alignment2._contigStart > alignment1._contigEnd) break;
					if(!_toBasespace.overlapOnTheReference(alignment1, alignment2)) continue;
				}
				else{
					if(alignment2._contigStart == alignment1._contigStart) continue;
					if(alignment2._contigEnd < alignment1._contigEnd) continue;

					if(alignment2._contigStart > alignment1._contigEnd+100) break;
				}

				nextAlignments.push_back(alignment2);
			}


			std::sort(nextAlignments.begin(), nextAlignments.end(), [](const ReadMapping2& a, const ReadMapping2& b){
				if(a._matchScore == b._matchScore){
					return a._readIndex < b._readIndex;
				}
				return a._matchScore > b._matchScore;
			});

	
			//bool cachedExists = false;
			//for(const ReadMapping2& alignment2: nextAlignments){
			//	if(readsVsReadsAlignmentsCached.find({readIndex1, alignment2._readIndex}) != readsVsReadsAlignmentsCached.end()){
			//		cachedExists = true;
			//	}
			//}

			//for(const ReadMapping2& alignment2: nextAlignments){	
			for(size_t i=0; i<nextAlignments.size() ; i++){ //&& i < 10
				const ReadMapping2& alignment2 =  nextAlignments[i];
				//for(size_t j=i+1; j<alignments.size(); j++){
				//const ReadMapping2& alignment2 = alignments[j];

				//cout << "\tTrying:\t" << alignment2._readIndex << "\t" << alignment2._contigStart << "\t" << alignment2._contigEnd << endl;

				//if(alignment2._contigStart > alignment1._contigEnd) break;
				//if(!overlapOnTheReference(alignment1, alignment2)) continue;

				//const ReadMapping2& alignment2 = alignment22._mapping;

				const ReadType readIndex2 = alignment2._readIndex;
				const string& read2 = _toBasespace._mReads[readIndex2];
				//string read2 = string(dnaStr2, strlen(dnaStr2));
				//free(dnaStr2);

				//if(alignment2._isReversed){
				//	Utils::toReverseComplement(read2);
				//}

				//int64_t maxContigSeq = min((int64_t) contigSequence.size(), (int64_t)(read2.size()+5000));
				//string contigSeq = contigSequence.substr(contigSequence.size()-maxContigSeq, maxContigSeq);

				const pair<ReadType, ReadType> readPair = {readIndex1, readIndex2};
				vector<AlignmentBounds> allAlignments;
				//AlignmentBounds alignment;

				if(readsVsReadsAlignmentsCached.find(readPair) == readsVsReadsAlignmentsCached.end()){
					_toBasespace.computeAlignment(read1, read2, _tbuf, _toBasespace._minimap2Preset_ava, true, ToBasespace2::_minOverlap, mi, allAlignments);
					readsVsReadsAlignmentsCached[readPair] = allAlignments;
					//cout << "derp:\t" << readIndex1 << "\t" << readIndex2 << "\t" << cachedExists << endl;
				}
				else{
					allAlignments = readsVsReadsAlignmentsCached[readPair];
					//cout << "OK:\t" << readIndex1 << "\t" << readIndex2 << "\t" << cachedExists << endl;
				}

				AlignmentBounds alignment;
				u_int32_t minLength = 0;

				for(const AlignmentBounds& al : allAlignments){

					if(allowErroneousReads){
						if(!_toBasespace.isValidOverlapAlignment(al, false)) continue;
					}
					else{
						if(!_toBasespace.isValidOverlapAlignment(al, true)) continue;
					}

					u_int32_t alignLength = min(al._queryEnd - al._queryStart, al._referenceEnd - al._referenceStart);
					
					if(alignLength > minLength){
						minLength = alignLength;
						alignment = al;
					}
				}
				

				//if(_toBasespace._print_debug) cout << "\t" << _toBasespace._debug_readHeaders[alignment2._readIndex] << "\t" << readIndex1 << " -> " << readIndex2 << " " << alignment._queryStart << " " << alignment2._isReversed << endl;

				//cout << read1 << endl;
				//cout << read2 << endl;
				//getchar();
				//nbAlignments += 1;

				//cout << "\t\t" << j << " " << alignment2._readIndex << "     " << alignment._queryStart << endl;
			

				if(alignment._queryStart == -1){
					//cout << "\tNo valid overlap" << endl;
					continue;
				}


				//if(readIndex2 == 866271){
				//	cout << "haaaaaaaaaaaa" << endl;
				//	cout << isErroneousRead(readIndex_to_i[readIndex2], alignments) << endl;
				//	getchar();
				//}

				if(!allowErroneousReads && _toBasespace.isErroneousRead(readIndex_to_i[readIndex2], alignments, isErroneousReadCached, contigCoverage, _tbuf)){
					//cout << "\tErroneous" << endl;
					//if(readIndex2 == 838925 || readIndex2 == 864136 || readIndex2 == 2221053){
					//	cout << "derp444: "<< readIndex2 << endl;
					//	getchar();
					//}
					//excludedReadIndexes.insert(readIndex2);
					continue;
				}
				
				//if(readIndex1 == 108086){
				//	cout << "Overlap on the ref only: " << overlapOnTheReferenceOnly << endl;
				//	cout << "Allow error: " << allowErroneousReads << endl;
				//	cout << readIndex2 << endl;
				//	for(const auto& it : allAlignments){
				//		cout << it._queryStart << "-" << it._queryEnd << "-" << it._isReversed << "   " << it._referenceStart << "-" << it._referenceEnd<< endl;
				//	}
				//	cout << "Used: " << alignment._referenceStart << "-" << alignment._referenceEnd << "-" << alignment._isReversed << endl;
				//}
				
				//cout << readIndex1 << " -> " << readIndex2 << " " << alignment._referenceStart << "-" << alignment._referenceEnd << endl;
				//getchar();

				//if(alignment._queryEnd - alignment._queryStart > maxAlignLength){
				//maxAlignLength = alignment._queryEnd - alignment._queryStart;
				//bestAlignment = alignment;
				//bestRead = read2;
				//i2 = alignment22._j;
				//lastAlignment = alignment2;
				//lastRead = read2;
				//foundSuccessor = true;

				//}
				readsVsReadsAlignmentsUsed[readPair] = alignment;
				
				return alignment2;
				//cout << "\t" << j << " " << alignment2._readIndex << "     " << alignment1._contigStart << "-" << alignment1._contigEnd << "\t" << alignment2._contigStart << "-" << alignment2._contigEnd << endl;
			

			}

			ReadMapping2 nullMapping;
			nullMapping._contigIndex = -1;
			return nullMapping;
		}

		void readPathsToContigs(const KminmerList& kminmerList, const vector<ReadPath>& readPaths, auto& readsVsReadsAlignments, auto& readIndex_to_alignment, const float contigCoverage){

			//if(kminmerList._readMinimizers.size() > 10000){
			//	cout << "SKIP LONG CONTIG" << endl;
			//	return;
			//}


			u_int64_t outputBps = 0;

			u_int32_t contigIndex = kminmerList._readIndex;
			bool isCircular = kminmerList._isCircular;
			if(readPaths.size() > 1) isCircular = false;

			if(readPaths.size() == 0) return;

			for(size_t pathIndex=0; pathIndex < readPaths.size(); pathIndex++){

				const ReadPath& readPath = readPaths[pathIndex];
				
				
				string contigSequence = "";

				if(readPath._readPath.size() == 1){
					
					const ReadType readIndex = readPath._readPath[0];
					const ReadMapping2& mapping1 = readIndex_to_alignment[readIndex];

					const string& read1 = _toBasespace._mReads[readIndex];
					//string read1 = string(dnaStr1, strlen(dnaStr1));
					//free(dnaStr1);
					//if(mapping1._isReversed){
					//	Utils::toReverseComplement(read1);
					//}

					contigSequence = read1;

				}
				else{

					for(size_t i=0; i<readPath._readPath.size()-1; i++){
						
						const ReadType readIndex1 = readPath._readPath[i];
						const ReadType readIndex2 = readPath._readPath[i+1];
						
						if(_toBasespace._print_debug) cout << i << "\t" << readIndex1 << " -> " << readIndex2 << endl;

						const pair<ReadType, ReadType> readPair = {readIndex1, readIndex2};
						if(readsVsReadsAlignments.find(readPair) == readsVsReadsAlignments.end()){
							cout << "Alignment not found" << endl;
							//getchar();
						}
						const AlignmentBounds& aligment = readsVsReadsAlignments[readPair];
						const ReadMapping2& mapping1 = readIndex_to_alignment[readIndex1];
						const ReadMapping2& mapping2 = readIndex_to_alignment[readIndex2];

						const string& read1 = _toBasespace._mReads[readIndex1];
						//string read1 = string(dnaStr1, strlen(dnaStr1));
						//free(dnaStr1);
						//if(mapping1._isReversed){
						//	Utils::toReverseComplement(read1);
						//}

						const string& read2 = _toBasespace._mReads[readIndex2];
						//string read2 = string(dnaStr2, strlen(dnaStr2));
						//free(dnaStr2);
						//if(mapping2._isReversed){
						//	Utils::toReverseComplement(read2);
						//}

						if(i == 0){
							contigSequence += read1;
							//contigSequence += read1.substr(mapping1._readStartReal, read1.size());
							if(_toBasespace._print_debug) cout << "\tContig size init: " << contigSequence.size() << endl;
							
						}
							
						u_int32_t prevOverhang =  aligment._referenceLength - aligment._referenceEnd;
						if(prevOverhang > 0) contigSequence = contigSequence.substr(0, contigSequence.size()-prevOverhang);

						string readSeqEnd = read2.substr(aligment._queryEnd, read2.size());
						contigSequence += readSeqEnd;

						if(_toBasespace._print_debug) cout << aligment._referenceStart << "-" << aligment._referenceEnd << "-" << aligment._referenceLength << "\t" << aligment._queryStart << "-" << aligment._queryEnd << "-" << aligment._queryLength << endl;
						if(_toBasespace._print_debug) cout << "\tContig size: " << contigSequence.size() << "  (Added: " << readSeqEnd.size() << ")" << endl;
						
						//getchar();


					}

				}
				




				const ReadMapping2& alignmentStart = readIndex_to_alignment[readPath._readPath[0]];
				u_int32_t oversizeStart = alignmentStart._readStartReal;
				if(alignmentStart._isReversed){
					oversizeStart = _toBasespace._mReads[alignmentStart._readIndex].size() - alignmentStart._readEndReal;
				}

				const ReadMapping2& alignmentEnd = readIndex_to_alignment[readPath._readPath[readPath._readPath.size()-1]];
				u_int32_t oversizeEnd = _toBasespace._mReads[alignmentEnd._readIndex].size() - alignmentEnd._readEndReal;
				if(alignmentEnd._isReversed){
					oversizeEnd = alignmentEnd._readStartReal;
				}

				if(isCircular){ //We leave at most 1000 extra overlaping bps at contig ends which will be removed accurately in the contig trimmer

					if(oversizeStart > 1000){
						oversizeStart = oversizeStart-1000;
					}
					else{
						oversizeStart = 0;
					}

					if(oversizeEnd > 1000){
						oversizeEnd = oversizeEnd-1000;
					}
					else{
						oversizeEnd = 0;
					}

				}

				/*
				cout << "aaaa: " << alignmentStart._readStartReal << " " << alignmentStart._readEndReal << " " << _toBasespace._mReads[alignmentStart._readIndex].size() << endl;
				cout << "bbbb: " << alignmentEnd._readStartReal << " " << alignmentEnd._readEndReal << " " << _toBasespace._mReads[alignmentEnd._readIndex].size() << endl;

				cout << oversizeStart << endl;
				cout << oversizeEnd << endl;

				if(isCircular){
					cout << endl;
					cout << "CONTIG: " << kminmerList._readIndex << " " << kminmerList._readMinimizers.size() << endl;
					for(size_t i=0; i<200; i++){
						cout << kminmerList._readMinimizers[i] << endl;
					}
					cout << endl;
					for(size_t i=kminmerList._readMinimizers.size()-10; i<kminmerList._readMinimizers.size(); i++){
						cout << kminmerList._readMinimizers[i] << endl;
					}

					cout << alignmentStart._readIndex << " " << alignmentEnd._readIndex << endl;
					//cout << alignmentEnd._readStartReal << " " << alignmentEnd._readEndReal << endl;
					//cout << _toBasespace._mReads[alignmentEnd._readIndex].size() << " " << alignmentEnd._readLengthBp << endl;
					//cout << oversizeEnd << endl;
					
					//cout << contigSequence.size() << endl;

					cout << _toBasespace._mReads[alignmentStart._readIndex] << endl;
					cout << _toBasespace._mReads[alignmentEnd._readIndex] << endl;
					
					//cout << "aaa:" << alignmentStart._readStartReal << " " << alignmentStart._readEndReal << endl;
					cout << alignmentStart._isReversed << " " << alignmentEnd._isReversed << endl;
				}
				*/

				contigSequence = contigSequence.substr(oversizeStart, contigSequence.size()-oversizeStart-oversizeEnd);
				//contigSequence = contigSequence.substr(0, contigSequence.size()-oversize);

				/*
				if(isCircular){
					cout << oversizeStart << endl;
					cout << contigSequence.size() << endl;
					getchar();
				}
				*/

				if(contigSequence.size() < _toBasespace._minContigLength) continue;




				if(computeSequenceComplexity(contigSequence, 64, 32) > 8 && contigCoverage < 6 && contigSequence.size() < 50000){
					//cout << read._seq << endl;
					//cout << alignments.size() << endl;
					//getchar();
					continue;
				}

				
				//if(contigHeader._isCircular){

				//cout << header << endl;


				bool isInvalid = false;
				bool isRepetitive = false;
				int nbIters = 0;
				
				while(true){

					u_int32_t mostAbundantRepeat = isHighlyRepetitive(contigSequence);

					if(mostAbundantRepeat == -1 && contigCoverage < 10){ //Super highly repetitive and low depth are discarded
						isInvalid = true;
						break; 
					}
					//if(!contigHeader._isCircular){
					//	cout << mostAbundantRepeat << endl;
					//}

					//cout << mostAbundantRepeat << endl;
					//cout << contigSequenceTrimmed.size() << endl;

					if(mostAbundantRepeat < 20) break;
					if(contigSequence.size() < 1000) break;

					//cout << "lul" << endl;
					float lengthToRemove = contigSequence.size()*0.1;
					contigSequence = contigSequence.substr(0, contigSequence.size()-lengthToRemove);

					nbIters += 1;
					isRepetitive = true;

					if(nbIters > 1000) break;
				}

				if(isInvalid) continue;
				if(contigSequence.size() < _toBasespace._minContigLength) continue;

				
				//if(!isCircular){
				if(isRepetitive){
					u_int32_t selfOverlapLength = computeSelfOverlap(contigSequence, _tbuf, _toBasespace._minimap2Preset_map);

					if(selfOverlapLength > 0){
						contigSequence = contigSequence.substr(0, contigSequence.size()-selfOverlapLength);
					}
				}

				if(contigSequence.size() < _toBasespace._minContigLength) continue;

				//if(!isCircular && isHighlyRepetitive(contigSequence)) continue;
				//"reprise: voir pourquoi on peut avoir des overlap massif quand circular"
				//u_int32_t selfOverlapLength = _toBasespace.computeSelfOverlap(contigSequence, _tbuf, _toBasespace._minimap2Preset_map);
				//cout << "self overlap: " << selfOverlapLength << endl;
				//if(selfOverlapLength > 0){
					//contigSequence = contigSequence.substr(0, contigSequence.size()-selfOverlapLength);
				//}

				//if(contigSequence.size() < _minContigLength) return;

				#pragma omp critical(writeContig)
				{
					
					if(_toBasespace._print_debug_main) cout << "\tWrite contig: " << contigSequence.size() << endl;
					outputBps += contigSequence.size();
					//if(contigSequence.size() > 200000) getchar();
					//if(contigSequence.size() == 0){
						//Logger::get().debug() << "Empty contig " << kminmerList._kminmersInfo.size();
						//return;
					//}
					//else{
					//cout << "lala:" << ((u_int32_t) isCircular) << endl;
					string linearOrCircular;
					if(isCircular){
						linearOrCircular = "c";
					}
					else{
						linearOrCircular = "l";
					}
					
					u_int64_t checksum = 0;
					for(size_t i=0; i<contigSequence.size(); i++){
						checksum += (contigSequence[i]*contigSequence.size()*contigIndex);
					}

					#pragma omp atomic
					_toBasespace._checksum_curatedContigs += checksum;

					#pragma omp atomic
					_toBasespace._checksum_curatedContigs_total += checksum;

					//for(char nt : contigSequence){
					//	_toBasespace._checksum_contigs += (nt*contigSequence.size());
					//}

					//_logFile << contigSequence.size() << endl;

					//string header = ">ctg" + to_string(contigIndex) + linearOrCircular + "_" + to_string(pathIndex) + '\n';
					string header = ">ctg" + to_string(_toBasespace._contigIndex) + linearOrCircular + '\n';
					bgzf_write(_toBasespace._outputContigFile, (const char*)&header[0], header.size());
					contigSequence +=  '\n';
					bgzf_write(_toBasespace._outputContigFile, (const char*)&contigSequence[0], contigSequence.size());


					//cout << kminmerList._readIndex << "\t" << kminmerList._readMinimizers.size() << "\t" << readPath._contigStart << "\t" << readPath._contigEnd << endl;
					vector<MinimizerType> contigMinimizers;
					for(size_t i=readPath._contigStart; i<=readPath._contigEnd; i++){
						contigMinimizers.push_back(kminmerList._readMinimizers[i]);
					}

					u_int32_t contigSize = contigMinimizers.size();
					_toBasespace._outputMinimizerContigFile.write((const char*)&contigSize, sizeof(contigSize));
					_toBasespace._outputMinimizerContigFile.write((const char*)&isCircular, sizeof(isCircular));
					_toBasespace._outputMinimizerContigFile.write((const char*)&contigMinimizers[0], contigSize*sizeof(MinimizerType));
					//_toBasespace._outputMinimizerContigFile.write((const char*)&contigIndex, sizeof(contigIndex));


					for(const ReadType readIndex : readPath._readPath){

						const ReadMapping2& mapping1 = readIndex_to_alignment[readIndex];

						const string& read1 = _toBasespace._mReads[readIndex];
						//string read1 = string(dnaStr1, strlen(dnaStr1));
						//free(dnaStr1);
						//if(mapping1._isReversed){
						//	Utils::toReverseComplement(read1);
						//}

						string readHeader = ">read_" + to_string(readIndex) + '\n';
						bgzf_write(_toBasespace._usedReadFile, (const char*)&readHeader[0], readHeader.size());
						string readSequence = read1 + '\n';
						bgzf_write(_toBasespace._usedReadFile, (const char*)&readSequence[0], readSequence.size());

						//_toBasespace._usedReadFile << ">read_" << readIndex << endl;
						//_toBasespace._usedReadFile << read1 << endl;
					}

					_toBasespace._nbBps += contigSequence.size();
					_toBasespace._contigIndex += 1;
				}


			}

			if(_toBasespace._print_debug_main){
				#pragma omp critical(writeContig)
				{
					long double expectedBp = kminmerList._readMinimizers.size()*_toBasespace._averageDistanceBetweenMinimizers;
					cout << "\tOutput bp fraction: " << outputBps << " " << outputBps / expectedBp << endl;
					if(outputBps < 1000000){
						cout << "\tExpected bp: " << kminmerList._readMinimizers.size()*_toBasespace._averageDistanceBetweenMinimizers << endl;
						cout << "\tOutput bp: " << outputBps << endl;
						//getchar();
					}
				}
			}
		}


		u_int32_t isHighlyRepetitive(const string& readSeq){

			//KmerModel kmerModel(21);

			vector<u_int64_t> kmers;
			vector<u_int8_t> kmerDirections;
			_modelRepetitive->iterate(readSeq.c_str(), readSeq.size(), kmers, kmerDirections);

			unordered_map<u_int64_t, u_int32_t> mCounts;
			for(const auto& it : kmers){
				mCounts[it] += 1;
			}

			long double nbRepeatedKmers = 0;
			for(const auto& it : kmers){
				if(mCounts[it] > 1) nbRepeatedKmers += 1;
			}

			vector<KmerAbundance> kmerAbundances;
			for(const auto& it : mCounts){
				kmerAbundances.push_back({it.first, it.second});
				//cout << it.first << "\t" << it.second << endl;
			}

			std::sort(kmerAbundances.begin(), kmerAbundances.end(), [&](const KmerAbundance& a, const KmerAbundance& b) {
				return a._abundance > b._abundance;
			});

			//cout << nbRepeatedKmers / kmers.size() << endl;
			//for(size_t i=0; i<kmerAbundances.size() && i < 20; i++){
			//	cout << kmerAbundances[i]._kmer << "\t" << kmerAbundances[i]._abundance << endl;
			//}

			//cout << nbRepeatedKmers / kmers.size() << endl;
			if(nbRepeatedKmers / kmers.size() > 0.9) return -1;
			if(nbRepeatedKmers / kmers.size() > 0.4) return kmerAbundances[0]._abundance;

			return 0;
		}


		double computeSequenceComplexity(const string& seq, const u_int64_t& w, const u_int64_t& step){

			float maxScore = 0;

			//KmerModel kmerModel(3);

			double l = w-2; //The number of trinucletoide kmers in a window

			//cout << "----" << endl;
			//cout << "a" << endl;

			vector<u_int64_t> kmers;
			vector<u_int8_t> kmerDirections;
			_modelTri->iterate(seq.c_str(), seq.size(), kmers, kmerDirections);

			//cout << "b" << endl;
			//cout << seq.size() << " " << kmers.size() << endl; 

			double nbWindows = 0;
			double windowScoreSum = 0;

			for(size_t ii=0; ii<kmers.size(); ii += step){
				
				vector<double> kmerCounts(64, 0); //64 = 4^3 (number of possible trinucletoide kmers)
				vector<double> sa(64, 0);

				int nbKmers = 0;
				for(size_t i=ii; i<kmers.size(); i++){
					//cout << "\t" << i << " " << kmers[i] << endl;
					kmerCounts[kmers[i]] += 1;
					nbKmers += 1;
					if(nbKmers == w) break;
				}

				//cout << "\t" << ii << " " << nbKmers << " " << w << endl;

				if(nbKmers < w) continue; //ignore end of the read


				for(size_t i=0; i<kmerCounts.size(); i++){
					sa[i] = kmerCounts[i] * (kmerCounts[i]-1) / 2.0;
				}

				double score = 0;

				for(size_t i=0; i<kmerCounts.size(); i++){
					score += sa[i];
				}

				score /= (l-1);

				if(score > maxScore){
					maxScore = score;
				}
				//if(score > 5) return score;
				//cout << "\t" << ii << " " << score << endl;
				//getchar();

				nbWindows += 1;
				windowScoreSum += score;
			}


			return maxScore; //windowScoreSum / nbWindows;
		}


		u_int32_t computeSelfOverlap(const string& sequence, mm_tbuf_t* tbuf, const string preset){

			string fakeName = "target";
			
			mm_idxopt_t iopt;
			mm_mapopt_t mopt;
			mm_set_opt(0, &iopt, &mopt); //"ava-ont"
			mm_set_opt(preset.c_str(), &iopt, &mopt); //"ava-ont"
			iopt.batch_size = 0x7fffffffffffffffL; //always build a uni-part index

			//if(performBaseLevelAlignment){
			mopt.flag |= MM_F_CIGAR; // perform alignment 
			//}
			
			mopt.min_chain_score = 500; //-m 500

			mopt.flag |= MM_F_ALL_CHAINS | MM_F_NO_DIAG | MM_F_NO_DUAL | MM_F_NO_LJOIN; // -D -P --no-long-join --dual=no

			mm_idx_t* mi = _toBasespace.minimap2index(iopt.w, iopt.k, iopt.flag&1, iopt.bucket_bits, sequence);
			
			//mopt.mid_occ = 5;

			mm_mapopt_update(&mopt, mi);

			mm_reg1_t *reg;
			//mm_tbuf_t *tbuf = mm_tbuf_init();
			int j, i, n_reg;
			reg = mm_map(mi, sequence.size(), sequence.c_str(), &n_reg, tbuf, &mopt, fakeName.c_str()); // get all hits for the query

			u_int32_t maxSelfOverlapLength = 0;

			//if(n_reg > 1) cout << "----" << endl;
			for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
			
				mm_reg1_t *r = &reg[j];

				free(r->p);
				
				if(r->rev) continue; 
				if(r->qs > 50) continue; //Looking for alignment at the contig start
				if(sequence.size() - r->re > 50) continue; //Looking for alignment at the contig end

				//cout << r->qs << "\t" << r->qe << "\t" << r->rs << "\t" << r->re << endl;



				//cout << "Found self overlap:\t" << r->qs << "\t" << r->qe << "\t" << r->rs << "\t" << r->re << endl;

				u_int32_t overlapStart = r->qe;
				u_int32_t overlapEnd = sequence.size() - r->rs;

				u_int32_t selfOverlapLength = max(overlapStart, overlapEnd);
				if(selfOverlapLength >= sequence.size()) continue;

				if(selfOverlapLength > maxSelfOverlapLength){
					maxSelfOverlapLength = selfOverlapLength;
				}

			}
			
		
			free(reg);
			mm_idx_destroy(mi);
			

			return maxSelfOverlapLength;
		}

		/*
		bool isHighlyRepetitive(const string& readSeq){

			vector<u_int64_t> kmers;
			vector<u_int8_t> kmerDirections;
			_modelRepetitive->iterate(readSeq.c_str(), readSeq.size(), kmers, kmerDirections);

			unordered_map<u_int64_t, u_int32_t> mCounts;
			for(const auto& it : kmers){
				mCounts[it] += 1;
			}

			long double nbRepeatedKmers = 0;
			for(const auto& it : kmers){
				if(mCounts[it] > 1) nbRepeatedKmers += 1;
			}

			vector<KmerAbundance> kmerAbundances;
			for(const auto& it : mCounts){
				kmerAbundances.push_back({it.first, it.second});
				//cout << it.first << "\t" << it.second << endl;
			}

			std::sort(kmerAbundances.begin(), kmerAbundances.end(), [&](const KmerAbundance& a, const KmerAbundance& b) {
				return a._abundance > b._abundance;
			});

			//for(size_t i=0; i<kmerAbundances.size() && i < 20; i++){
			//	cout << kmerAbundances[i]._kmer << "\t" << kmerAbundances[i]._abundance << endl;
			//}

			
			if(nbRepeatedKmers / kmers.size() > 0.6) return true;

			return false;
		}
		*/
		/*
		bool overlapOnTheReference(const ReadMapping2& prevAlignment, const ReadMapping2& nextAlignment){


			int64_t prevLeftOverhang = 0;
			int64_t prevRightOverhang = 0;

			if(prevAlignment._isReversed){
				prevLeftOverhang = prevAlignment._readLengthBp - prevAlignment._readEndReal;
				prevRightOverhang = prevAlignment._readStartReal;
			}
			else{
				prevLeftOverhang = prevAlignment._readStartReal;
				prevRightOverhang = prevAlignment._readLengthBp - prevAlignment._readEndReal;
			}

			int64_t nextLeftOverhang = 0;
			int64_t nextRightOverhang = 0;

			if(nextAlignment._isReversed){
				nextLeftOverhang = nextAlignment._readLengthBp - nextAlignment._readEndReal;
				nextRightOverhang = nextAlignment._readStartReal;
			}
			else{
				nextLeftOverhang = nextAlignment._readStartReal;
				nextRightOverhang = nextAlignment._readLengthBp - nextAlignment._readEndReal;
			}



			int64_t prevStart = ((int64_t) prevAlignment._contigStart*_toBasespace._averageDistanceBetweenMinimizers) - prevLeftOverhang;
			int64_t prevEnd = ((int64_t) prevAlignment._contigEnd*_toBasespace._averageDistanceBetweenMinimizers) + prevRightOverhang;

			int64_t nextStart = ((int64_t) nextAlignment._contigStart*_toBasespace._averageDistanceBetweenMinimizers) - nextLeftOverhang;
			int64_t nextEnd = ((int64_t) nextAlignment._contigEnd*_toBasespace._averageDistanceBetweenMinimizers) + nextRightOverhang;


			int64_t alignmentOffset = _minOverlap;

			if (nextStart > prevStart+alignmentOffset && nextStart < prevEnd-alignmentOffset && nextEnd > prevEnd+alignmentOffset) return true;

			return false;
		}
		*/



	};



	bool isErroneousRead(const u_int32_t ii, const vector<ReadMapping2>& alignments, auto& isErroneousReadCached, const float contigCoverage, mm_tbuf_t* tbuf){

		int64_t usedCoverage = 10;
		u_int32_t largeDeletions = 0;
		u_int32_t largeInsertions = 0;

		const ReadMapping2& alignment1 = alignments[ii];
		const ReadType readIndex1 = alignment1._readIndex;

		if(isErroneousReadCached.find(readIndex1) != isErroneousReadCached.end()){
			return isErroneousReadCached[readIndex1];
		}

		//cout << readIndex1 << endl;
		const string& read1 = _mReads[readIndex1];
		//string read1 = string(dnaStr1, strlen(dnaStr1));
		//free(dnaStr1);

		//cout << read1 << endl;

		//bool print = alignment1._readIndex == 1751275;

		//if(!print) return false;

		mm_idx_t* mi = minimap2_indexRead(read1);

		vector<u_int64_t> coveragesMapping(read1.size(), 0);
		vector<u_int64_t> coverages(read1.size(), 0);
		//vector<u_int64_t> deletions(read1.size(), 0);
		//vector<u_int64_t> insertions(read1.size(), 0);
		
		vector<ReadMapping2> selectedAlignments;
		subsampleMappedReads(ii, alignments, selectedAlignments, usedCoverage);

		for(const ReadMapping2& alignment2 : selectedAlignments){
			
			//cout << i << " " << j << endl;
			//const ReadMapping2& alignment2 = alignments[j];
			
			if(alignment2._contigStart > alignment1._contigEnd) break;
			
			const ReadType readIndex2 = alignment2._readIndex;
			const string& read2 = _mReads[readIndex2];
			//string read2 = string(dnaStr2, strlen(dnaStr2));
			//free(dnaStr2);

			vector<AlignmentBounds> results;
			const AlignmentBounds& alignment = computeAlignment(read1, read2, tbuf, _minimap2Preset_map, true, ToBasespace2::_minOverlap, mi, results);
			if(alignment._queryStart == -1) continue;

			//for(size_t k=alignment._referenceStart; k<alignment._referenceEnd; k++){
			//	coverages[k] += 1;
			//}
			
			//if(print){
			
			u_int32_t referencePos = alignment._referenceStart;

			for(size_t k=0; k<alignment._cigarElements.size(); k++){
				const CigarElement& c = alignment._cigarElements[k];

				if(c._type == 'M'){
					for(size_t m=referencePos; m<referencePos+c._size; m++){
						coverages[m] += 1;
						coveragesMapping[m] += 1;
					}
					referencePos += c._size;
				}
				else if(c._type == 'I'){
					//if(c._size > 50){
					//	largeInsertions += 1;
					//	cout << "large insertion: " << c._size << " " << referencePos << endl;
					//}
				}
				else if(c._type == 'D'){
					//if(c._size > 50){
					//	largeDeletions += 1;
					//	cout << "large deletion: " << c._size << " " << referencePos << endl;
					//}
					//for(size_t m=referencePos; m<referencePos+c._size; m++){
						//if(k >= deletions.size()) cout << "derp: " << k << " " << deletions.size() << endl;
					//	deletions[m] += 1;
					//}
					for(size_t m=referencePos; m<referencePos+c._size; m++){
						coveragesMapping[m] += 1;
					}

					referencePos += c._size;
				}
			}
			
			//}
			

		}


		mm_idx_destroy(mi);

		/*
		if(alignment1._readIndex == 28263){

			cout << endl;
			for(size_t i=0; i<coverages.size(); i+=20){
				cout << i << "\t" << coverages[i] << endl;
			}
			
			cout << "is chimeric: " <<  isChimeric(coverages, contigCoverage, usedCoverage) << endl;
		}
		*/

		bool isError = isChimeric(coverages, coveragesMapping, contigCoverage, usedCoverage);
		isErroneousReadCached[readIndex1] = isError;

		return isError;


		/*

		int64_t deletionStart = -1;

		//cout << "lul" << endl;
		for(int64_t i=0; i<coverages.size(); i++){

			//cout << i << "\t" << coverages[i] << endl;
			//if(i != 0 && i%5000==0) getchar();
			if(coverages[i] == 0) continue;

			float deletionFraction = (double) deletions[i] / (double) coverages[i];

			if(deletionFraction > 0.9){
				//cout << "Delection:\t" << i << "\t" << deletionFraction << endl;
				if(deletionStart == -1){
					deletionStart = i;
				}
			}
			else{
				if(deletionStart != -1){
					int64_t deletionSize = i - deletionStart;
					if(deletionSize > 50) return true; 
					//cout << deletionSize << endl;
					deletionStart = -1;
				}
			}
			//cout << i << "\t" << coverages[i] << endl;
		}

		if(deletionStart != -1){
			int64_t deletionSize = ((int64_t)coverages.size()) - deletionStart;
			if(deletionSize > 50) return true; 
			//cout << deletionSize << endl;
		}

		//cout << largeInsertions << " " << largeDeletions << endl;
		//getchar();
		*/

		return false;
	}


	void subsampleMappedReads(const int64_t ii, const vector<ReadMapping2>& alignments, vector<ReadMapping2>& selectedAlignments, const int64_t usedCoverage){
		


		selectedAlignments.clear();

		unordered_set<ReadType> removedReadIndexes;

		const ReadMapping2& alignment1 = alignments[ii];
		const ReadType readIndex1 = alignment1._readIndex;

		//vector<u_int32_t> coverages(read1.size(), 0);
		vector<ReadMapping2> nextAlignments;
		u_int32_t contigStart = alignment1._contigStart;
		u_int32_t contigEnd = alignment1._contigEnd;


		//if(alignment1._readIndex == 2221053){
		//	cout << "lala: " << ii << endl;
		//	for(size_t i=0; i<alignments.size(); i++){
		//		cout << i << "\t" << alignments[i]._contigStart << "\t" << alignments[i]._contigEnd << endl; 
		//	}
		//}


		for(int64_t j=ii-1; j>=0; j--){
			
			//cout << i << " " << j << endl;
			const ReadMapping2& alignment2 = alignments[j];
			
			if(alignment2._contigEnd < alignment1._contigStart+3) continue;

			nextAlignments.push_back(alignment2);

			//cout << alignment2._readIndex << "\t" << alignment2._matchScore << "\t" << alignment2._contigStart << "\t" << alignment2._contigEnd << endl;
		}

		for(int64_t j=ii+1; j<alignments.size(); j++){
			
			//cout << i << " " << j << endl;
			const ReadMapping2& alignment2 = alignments[j];
			
			if(alignment2._contigStart+3 > alignment1._contigEnd) break;

			nextAlignments.push_back(alignment2);

			//cout << alignment2._readIndex << "\t" << alignment2._matchScore << "\t" << alignment2._contigStart << "\t" << alignment2._contigEnd << endl;
		}

		std::sort(nextAlignments.begin(), nextAlignments.end(), [](const ReadMapping2& a, const ReadMapping2& b){
			if(a._matchScore == b._matchScore){
				return a._readIndex < b._readIndex;
			}
			return a._matchScore < b._matchScore;
		});

		//if(alignment1._readIndex == 2221053){
		//	cout << "huuu: " << nextAlignments.size() << " " << alignment1._contigStart << " " << alignment1._contigEnd << endl;
			//getchar();
		//}

		vector<u_int32_t> coverages(contigEnd-contigStart, 0);

		for(const ReadMapping2& al : nextAlignments){
			
			//cout << "\t" << al._readIndex << "\t" << al._contigStart << "\t" << al._contigEnd << endl;

			for(size_t i=al._contigStart; i<al._contigEnd; i++){
				int64_t pos = i - contigStart;
				//cout << "\t\t" << pos << endl;
				if(pos < 0 || pos >= coverages.size()) continue;
				coverages[pos] += 1;
			}
		}

		//if(alignment1._readIndex == 2221053){
		//	cout << "huuu: " << nextAlignments.size() << " " << coverages.size() << endl;
		//	for(size_t i=0; i<coverages.size(); i++){
		//		cout << i << "\t" << coverages[i] << endl;
		//	}
		//}


		int64_t offset = 0;
		//int64_t minCoverage = 10;

		for(const ReadMapping2& alignment : nextAlignments){

			bool canRemove = true;
			bool isRemoveAllow = false;

			const ReadType readIndex = alignment._readIndex;
			int64_t start = alignment._contigStart;
			int64_t end = alignment._contigEnd;

			//cout << readIndex << " " << start << " " << end << endl;

			for(int64_t i=start+offset; i<end-offset; i++){
				int64_t pos = i - contigStart;

				if(pos < 0 || pos >= coverages.size()) continue;

				//cout << coverages[pos] << " " << usedCoverage << endl;
				isRemoveAllow = true;
				if(coverages[pos] <= usedCoverage){
					canRemove = false;
					break;
				}
				
				//if(pos >= coverages.size()-1) break;
			}

			//cout << isRemoveAllow << " " << canRemove << endl;
			if(isRemoveAllow && canRemove){
				removedReadIndexes.insert(readIndex);

				for(int64_t i=start+offset; i<end-offset; i++){
					
					int64_t pos = i - contigStart;
					if(pos < 0 || pos >= coverages.size()) continue;
					//if(pos >= coverages.size()) break;

					coverages[pos] -= 1;
				}
			}

		}

		vector<u_int32_t> coverages2(contigEnd-contigStart, 0);

		for(const ReadMapping2& al : nextAlignments){
			if(removedReadIndexes.find(al._readIndex) != removedReadIndexes.end()) continue;
			selectedAlignments.push_back(al);
			//cout << al._readIndex << "\t" << al._matchScore << "\t" << al._contigStart << "\t" << al._contigEnd << endl;


			for(int64_t i=al._contigStart; i<al._contigEnd; i++){
				
				int64_t pos = i - contigStart;
				if(pos >= coverages.size()) break;

				coverages2[pos] += 1;
			}

		}

		//if(alignment1._readIndex == 2221053){
			//cout << "endo: " << selectedAlignments.size() << endl;
			//getchar();
		//}

		//for(size_t i=0; i<coverages.size(); i++){
		//	cout << i << "\t" << coverages[i] << endl;
		//}
		//cout << endl;
		//for(size_t i=0; i<coverages2.size(); i++){
		//	cout << i << "\t" << coverages2[i] << endl;
		//}

		//cout << nextAlignments.size() << " " << selectedAlignments.size() << " " << removedReadIndexes.size() << endl;
		//getchar();
	}


	bool isChimeric(const vector<u_int64_t>& coverages, const vector<u_int64_t>& coveragesMapping, const float contigCoverage, const int64_t usedCoverage){
		//if(contigCoverage < 10) return false; 
		//if(contigCoverage < 10) return false; 
		//bool isChimeric = false;

		vector<MinimizerRead> reads;

		vector<CoverageRegion> coverageRegions = collectLowHighDepthRegions(coverages, contigCoverage, usedCoverage);


		//for(size_t i=0; i<referenceRead._minimizers.size(); i++){
		//	cout << "\t" << i << "\t" << coverages[i] << endl;
		//}

		//cout << endl;
		//for(size_t i=0; i<coverageRegions.size(); i++){
		//	cout << "\t" << i << "\t" << coverageRegions[i]._startIndex << "\t" << coverageRegions[i]._endIndex << "\t" << coverageRegions[i]._coverageType << endl;
		//}
		
		/*
		vector<MinimizerType> minimizers;
		vector<u_int32_t> pos;
		vector<u_int8_t> directions;
		vector<u_int8_t> qualities;

		for(int64_t i=0; i<coverageRegions.size(); i++){

			bool isChimeric = false;

			if(coverageRegions[i]._coverageType == CoverageRegion::low){
				if(i > 0){  
					if(coverageRegions[i-1]._coverageType == CoverageRegion::high && coverageRegions[i-1].getLength(referenceRead._minimizersPos) >= 200){
						//cout << "Chimeric position: " << coverageRegions[i-1]._endIndex << endl;
						isChimeric = true;
					}
				}

				if(i < ((int64_t)coverageRegions.size())-1){  
					if(coverageRegions[i+1]._coverageType == CoverageRegion::high && coverageRegions[i+1].getLength(referenceRead._minimizersPos) >= 200){
						//cout << "Chimeric position: " << coverageRegions[i+1]._startIndex << endl;
						isChimeric = true;
					}
				}

			}

			if(!isChimeric){
				//cout << "add: " << coverageRegions[i]._startIndex << " " << coverageRegions[i]._endIndex << endl;
				for(size_t j=coverageRegions[i]._startIndex; j<=coverageRegions[i]._endIndex; j++){
					minimizers.push_back(referenceRead._minimizers[j]);
					pos.push_back(referenceRead._minimizersPos[j]);
					directions.push_back(referenceRead._readMinimizerDirections[j]);
					qualities.push_back(referenceRead._qualities[j]);
				}
			}

			if(isChimeric){

				if(minimizers.size() > 0){
					//cout << "push: " << minimizers[0] << " " << minimizers[minimizers.size()-1];
					u_int32_t readLength = pos[pos.size()-1] - pos[0];
					MinimizerRead read = {referenceRead._readIndex, minimizers, pos, qualities, directions, readLength, referenceRead._meanReadQuality};
					reads.push_back(read);
				}

				minimizers.clear();
				pos.clear();
				directions.clear();
				qualities.clear();
			}

			//cout << "\t" << i << "\t" << coverageRegions[i]._startIndex << "\t" << coverageRegions[i]._endIndex << "\t" << coverageRegions[i]._coverageType << endl;
		}

		if(reads.size() > 0 && minimizers.size() > 0){ //output last minimizer read (if no chimeric region at all, we just output the origina referenceRead)
			u_int32_t readLength = pos[pos.size()-1] - pos[0];
			MinimizerRead read = {referenceRead._readIndex, minimizers, pos, qualities, directions, readLength, referenceRead._meanReadQuality};
			reads.push_back(read);
		}

		//for(size_t i=0; i<reads.size(); i++){
		//	cout << i << endl;
		//	for(size_t j=0; j<reads[i]._minimizers.size(); j++){
		//		cout << "\t" << j << "\t" << reads[i]._minimizers[j] << endl;
		//	}
		//}

		if(reads.size() == 0){
			return {referenceRead};
		}
		

		if(reads.size() > 1){
			//getchar();
		}

		return reads;
		*/


		for(int64_t i=0; i<coverageRegions.size(); i++){

			if(coverageRegions[i]._coverageType == CoverageRegionType::Low && coverageRegions[i].getLength() >= 200){

				if(contigCoverage < 10){
					if(coverageRegions[i].isSupportedByRead(coveragesMapping)){
						return true;
					}
					else{
						return false;
					}
				}

				return true;
			}
			/*
			bool isChimeric = false;


			if(coverageRegions[i]._coverageType == CoverageRegion::low){
				if(i > 0){  
					if(coverageRegions[i-1]._coverageType == CoverageRegion::high && coverageRegions[i-1].getLength() >= 200){
						//cout << "Chimeric position: " << coverageRegions[i-1]._endIndex << endl;
						isChimeric = true;
					}
				}

				if(i < ((int64_t)coverageRegions.size())-1){  
					if(coverageRegions[i+1]._coverageType == CoverageRegion::high && coverageRegions[i+1].getLength() >= 200){
						//cout << "Chimeric position: " << coverageRegions[i+1]._startIndex << endl;
						isChimeric = true;
					}
				}

			}

			if(isChimeric) return true;
			*/
		}


		return false;
	}

	vector<CoverageRegion> collectLowHighDepthRegions(const vector<u_int64_t>& coverages, const float contigCoverage, const float usedCoverage){
		
		float minCoverage = 0;
		if(contigCoverage > 30) minCoverage = 1;
		if(contigCoverage > 70) minCoverage = 2;
		if(contigCoverage > 200) minCoverage = 3;

		//float minCoverage = min(contigCoverage*0.2, usedCoverage*0.2);

		vector<int> coverageTypes;

		for(size_t i=0; i<coverages.size(); i++){

			u_int64_t coverage = coverages[i];

			//if(coverage <= minCoverage){
			if(coverage <= minCoverage){
				coverageTypes.push_back(CoverageRegionType::Low);
			}
			else{
				coverageTypes.push_back(CoverageRegionType::Normal);
			}
			//cout << i << "\t" << coverageBreadth[i] << endl;

			/*
			if(coverage <= 2){
				coverageTypes.push_back(CoverageRegion::low);
			}
			else if(coverage > 5*2){
				coverageTypes.push_back(CoverageRegion::high);
			}
			else{
				coverageTypes.push_back(CoverageRegion::normal);
			}
			*/
		}

		vector<CoverageRegion> regions;



		int lastChar = -1;
		u_int64_t lastPos = 0;

		for(size_t i=0; i<coverageTypes.size(); i++){
			
			int c = coverageTypes[i];

			if(c == lastChar) continue;

			if(lastChar != -1){

				regions.push_back({lastPos, i-1, coverageTypes[lastPos]});
				//rleSequence += lastChar;
				//rlePositions.push_back(lastPos);
				lastPos = i;
				//cout << lastChar << endl;
			}
			lastChar = c;
		}

		regions.push_back({lastPos, coverageTypes.size()-1, coverageTypes[lastPos]});
		//cout << lastChar << endl;
		//rleSequence += lastChar;
		//rlePositions.push_back(lastPos);
		//rlePositions.push_back(length);


		return regions;
	}


	mm_idx_t* minimap2_indexRead(const string& reference){

		//string preset = _minimap2Preset_ava
		mm_idxopt_t iopt;
		mm_mapopt_t mopt;
		mm_set_opt(0, &iopt, &mopt); //"ava-ont"
		//mm_set_opt(_minimap2Preset_ava.c_str(), &iopt, &mopt); //"ava-ont"
		iopt.batch_size = 0x7fffffffffffffffL; //always build a uni-part index

		//if(useBaseLevelAlignment){
		//	mopt.flag |= MM_F_CIGAR; // perform alignment 
		//}

		//mopt.min_chain_score = 500; //-m 500
		//iopt.bucket_bits = 14;

		//cout << iopt.w << endl;
		//cout << iopt.k << endl;
		//cout << iopt.bucket_bits << endl;
		//cout << (iopt.flag&1) << endl;
		//cout << mopt.min_chain_score << endl;

		//mm_idxopt_t idx_opt;
		mm_idx_t* mi = minimap2index(iopt.w, iopt.k, iopt.flag&1, iopt.bucket_bits, reference);
		return mi;
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

	AlignmentBounds computeAlignment(const string& reference, const string& query, mm_tbuf_t* tbuf, const string preset, bool performBaseLevelAlignment, const size_t minChainScore, mm_idx_t* mi, vector<AlignmentBounds>& allAlignments){

		//mm_tbuf_t* tbuf = mm_tbuf_init();

		allAlignments.clear();
		//cout << "Target size: " << reference.size() << endl;
		//cout << "Query size: " << query.size() << endl;
		
		AlignmentBounds alignment;
		alignment._referenceLength = reference.size();
		alignment._queryLength = query.size();

		//string preset = _minimap2Preset_ava
		mm_idxopt_t iopt;
		mm_mapopt_t mopt;
		mm_set_opt(0, &iopt, &mopt); //"ava-ont"
		mm_set_opt(preset.c_str(), &iopt, &mopt); //"ava-ont"
		iopt.batch_size = 0x7fffffffffffffffL; //always build a uni-part index

		if(performBaseLevelAlignment){
			mopt.flag |= MM_F_CIGAR; // perform alignment 
		}
		

		mopt.min_chain_score = minChainScore; //-m 500
		//iopt.bucket_bits = 14;

		//cout << iopt.w << endl;
		//cout << iopt.k << endl;
		//cout << iopt.bucket_bits << endl;
		//cout << (iopt.flag&1) << endl;
		//cout << mopt.min_chain_score << endl;

		//mm_idxopt_t idx_opt;
		//mm_idx_t* mi = minimap2index(iopt.w, iopt.k, iopt.flag&1, iopt.bucket_bits, reference);
		//mm_idx_t* mi = minimap2index(iopt.w, iopt.k, iopt.flag&1, iopt.bucket_bits, reference);
		mm_mapopt_update(&mopt, mi);
		mopt.mid_occ = 5;
		//mopt.mid_occ = 1000; // don't filter high-occ seeds

		mm_reg1_t *reg;
		//mm_tbuf_t *tbuf = mm_tbuf_init();
		int j, i, n_reg;
		reg = mm_map(mi, query.size(), query.c_str(), &n_reg, tbuf, &mopt, 0); // get all hits for the query

		//int32_t expectedOverlapLength = getEstimatedOverlapLength(alignment1, alignment2);
		//int32_t expectedQueryStart = alignment2._readStart;
		//int32_t expectedQueryEnd = expectedQueryStart + expectedOverlapLength;


		u_int32_t maxError = -1;
		u_int32_t minLength = 0;


		//if(n_reg > 1) cout << "----" << endl;
		for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
		
			mm_reg1_t *r = &reg[j];

			//string cigar = "";
			
			vector<CigarElement> cigar;
			if(performBaseLevelAlignment){
				for (i = 0; i < r->p->n_cigar; ++i){
					u_int32_t size = r->p->cigar[i]>>4;
					char type = MM_CIGAR_STR[r->p->cigar[i]&0xf];
					cigar.push_back({size, type});
					//cigar += to_string(size) + type;
					//cout << size << " " << type << endl;
				}
			}
			

			free(r->p);

			AlignmentBounds alignment2;
			alignment2._referenceLength = reference.size();
			alignment2._queryLength = query.size();
			alignment2._queryStart = r->qs;
			alignment2._queryEnd = r->qe;
			alignment2._referenceStart = r->rs;
			alignment2._referenceEnd = r->re;
			alignment2._isReversed = r->rev;
			alignment2._identity = ((double)r->mlen) / r->blen;// r->div;
			
			allAlignments.push_back(alignment2);
			//if(alignment2._isReversed){
				//free(reg);
				//return alignment;
				//continue;
			//}

			//cout << "\t\t" << alignment2._referenceStart << "-" << alignment2._referenceEnd << "-" << alignment2._referenceLength << "\t" << alignment2._queryStart << "-" << alignment2._queryEnd <<  "-" << alignment2._queryLength << endl;

			//if(isOverlap && !isValidOverlapAlignment(alignment2)){
				//free(reg);
				//return alignment;
			//	continue;
			//}
			//if(!isValidAlignment(alignment2, 0)) continue; //A remettre?

			int32_t actualOverlapLength = r->qe - r->qs;

			//int error = abs(expectedOverlapLength - actualOverlapLength) + abs(expectedQueryStart - r->qs);

			u_int32_t alignLength = min(r->qe - r->qs, r->re - r->rs);
			

			//if(error < maxError){
			if(alignLength > minLength){
				minLength = alignLength;
				//maxError = error;

				alignment._queryStart = r->qs;
				alignment._queryEnd = r->qe;
				alignment._referenceStart = r->rs;
				alignment._referenceEnd = r->re;
				alignment._isReversed = r->rev;
				alignment._identity = ((double)r->mlen) / r->blen;// r->div;
				alignment._nbMatches = r->mlen;
				alignment._cigarElements = cigar;
			}
			
		}
		
	
		free(reg);
		//mm_idx_destroy(mi);
		
		//mm_tbuf_destroy(tbuf);

		return alignment;
	}


	bool isValidOverlapAlignment(const AlignmentBounds& alignmentBounds, const bool checkMaxhang){

		if(alignmentBounds._queryStart == -1) return false;
		if(alignmentBounds._isReversed) return false; 

		//if(alignmentBounds._identity < ContigCurator::_minOverlapIdentity) return false;

		int32_t actualOverlapLength = alignmentBounds._queryEnd - alignmentBounds._queryStart;
		//int32_t error = abs(expectedOverlapLength - actualOverlapLength);
		
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

		if(checkMaxhang){
			if (ext5 > _maxHang || ext3 > _maxHang || queryEnd - queryStart < (queryEnd - queryStart + ext5 + ext3) * _intFrac) return false;
		}


		if (queryStart <= tl5 && queryLength - queryEnd <= tl3) return false; //query contained
		else if (queryStart >= tl5 && queryLength - queryEnd >= tl3) return false; //target contained

		if (queryEnd - queryStart + ext5 + ext3 < _minOverlap) return false; //Short overlap
		if (targetEnd - targetStart + ext5 + ext3 < _minOverlap) return false; //Short overlap


		//cout << queryLength << " " <<  targetLength  << " " << alignmentBounds._identity << " " << ContigCurator::_minOverlapIdentity << endl;

		return true;

	}


	bool alignmentOverlapExistingReadPath(const ReadMapping2& alignment, const vector<ReadPath>& readPaths){



		int64_t allowedOverlap = 0;

		for(const ReadPath& readPath : readPaths){

			if(alignment._contigStart >= readPath._contigStart && alignment._contigEnd <= readPath._contigEnd) return true; //alignment contained
			if(alignment._contigStart <= readPath._contigStart && alignment._contigEnd >= readPath._contigEnd) return true; //existing alignment contained

			if(alignment._contigStart >= readPath._contigStart){
				if(readPath._contigEnd - alignment._contigStart > allowedOverlap) return true;
			}

			if(alignment._contigEnd <= readPath._contigEnd){
				if(alignment._contigEnd - readPath._contigStart > allowedOverlap) return true;
			}

		}

		return false;
	}

	bool isAlignmentContainedInExistingReadPath(const ReadMapping2& alignment, const vector<ReadPath>& readPaths){
		for(const ReadPath& path : readPaths){
			if (alignment._contigStart >= path._contigStart && alignment._contigEnd <= path._contigEnd) return true;
		}

		return false;

	}

	
	bool overlapOnTheReference(const ReadMapping2& prevAlignment, const ReadMapping2& nextAlignment){


		int64_t alignmentOffset = 1;



		if (nextAlignment._contigStart > prevAlignment._contigStart+alignmentOffset && nextAlignment._contigStart < prevAlignment._contigEnd-alignmentOffset && nextAlignment._contigEnd > prevAlignment._contigEnd+alignmentOffset) return true;

		//if (nextAlignment._contigStart >= prevAlignment._contigStart && nextAlignment._contigStart <= prevAlignment._contigEnd && nextAlignment._contigEnd >= prevAlignment._contigEnd) return true;

		return false;
	}

	/*
	bool overlapOnTheReference(const ReadMapping2& prevAlignment, const ReadMapping2& nextAlignment){

		int64_t minOffset = 2500 / (1/_params._minimizerDensity_assembly);
		
		int64_t alignmentOffset = (prevAlignment._contigEnd-prevAlignment._contigStart)*0.5;
		alignmentOffset = max((int64_t)1, alignmentOffset);
		alignmentOffset = min((int64_t)minOffset, alignmentOffset);


		if (nextAlignment._contigStart > prevAlignment._contigStart+alignmentOffset && nextAlignment._contigStart < prevAlignment._contigEnd-alignmentOffset && nextAlignment._contigEnd > prevAlignment._contigEnd+alignmentOffset) return true;

		//if (nextAlignment._contigStart >= prevAlignment._contigStart && nextAlignment._contigStart <= prevAlignment._contigEnd && nextAlignment._contigEnd >= prevAlignment._contigEnd) return true;

		return false;
	}
	*/
	/*
	u_int32_t computeSelfOverlap(const string& sequence, mm_tbuf_t* tbuf, const string preset){

		string fakeName = "target";
		//mm_tbuf_t* tbuf = mm_tbuf_init();

		//allAlignments.clear();
		//cout << "Target size: " << reference.size() << endl;
		//cout << "Query size: " << query.size() << endl;
		
		//AlignmentBounds alignment;
		//alignment._referenceLength = reference.size();
		//alignment._queryLength = query.size();

		//string preset = _minimap2Preset_ava
		mm_idxopt_t iopt;
		mm_mapopt_t mopt;
		mm_set_opt(0, &iopt, &mopt); //"ava-ont"
		mm_set_opt(preset.c_str(), &iopt, &mopt); //"ava-ont"
		iopt.batch_size = 0x7fffffffffffffffL; //always build a uni-part index

		//if(performBaseLevelAlignment){
		mopt.flag |= MM_F_CIGAR; // perform alignment 
		//}
		

		mopt.flag |= MM_F_ALL_CHAINS | MM_F_NO_DIAG | MM_F_NO_DUAL | MM_F_NO_LJOIN; // -D -P --no-long-join --dual=no

		mm_idx_t* mi = minimap2index(iopt.w, iopt.k, iopt.flag&1, iopt.bucket_bits, sequence);
		//mm_idx_t* mi = minimap2_indexRead(sequence);
		//mopt.min_chain_score = minChainScore; //-m 500
		//iopt.bucket_bits = 14;

		//cout << iopt.w << endl;
		//cout << iopt.k << endl;
		//cout << iopt.bucket_bits << endl;
		//cout << (iopt.flag&1) << endl;
		//cout << mopt.min_chain_score << endl;

		//mm_idxopt_t idx_opt;
		//mm_idx_t* mi = minimap2index(iopt.w, iopt.k, iopt.flag&1, iopt.bucket_bits, reference);
		//mm_idx_t* mi = minimap2index(iopt.w, iopt.k, iopt.flag&1, iopt.bucket_bits, reference);
		mm_mapopt_update(&mopt, mi);
		//mopt.mid_occ = 5;
		//mopt.mid_occ = 1000; // don't filter high-occ seeds

		mm_reg1_t *reg;
		//mm_tbuf_t *tbuf = mm_tbuf_init();
		int j, i, n_reg;
		reg = mm_map(mi, sequence.size(), sequence.c_str(), &n_reg, tbuf, &mopt, fakeName.c_str()); // get all hits for the query

		//int32_t expectedOverlapLength = getEstimatedOverlapLength(alignment1, alignment2);
		//int32_t expectedQueryStart = alignment2._readStart;
		//int32_t expectedQueryEnd = expectedQueryStart + expectedOverlapLength;

		u_int32_t maxSelfOverlapLength = 0;

		//if(n_reg > 1) cout << "----" << endl;
		for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
		
			mm_reg1_t *r = &reg[j];


			//cout << string(mi->seq[r->rid].name) << " " << r->qs << "\t" << r->qe << "\t" << r->rs << "\t" << r->re << "\t" << r->blen << endl;

			free(r->p);
			
			if(r->rev) continue; 
			if(r->qs > 50) continue; //Looking for alignment at the contig start
			if(sequence.size() - r->re > 50) continue; //Looking for alignment at the contig end

			//cout << r->qs << "\t" << r->qe << "\t" << r->rs << "\t" << r->re << endl;



			cout << "Found self overlap:\t" << r->qs << "\t" << r->qe << "\t" << r->rs << "\t" << r->re << endl;

			u_int32_t overlapStart = r->qe;
			u_int32_t overlapEnd = sequence.size() - r->rs;

			u_int32_t selfOverlapLength = max(overlapStart, overlapEnd);
			if(selfOverlapLength >= sequence.size()) continue;

			if(selfOverlapLength > maxSelfOverlapLength){
				maxSelfOverlapLength = selfOverlapLength;
			}

		}
		
	
		free(reg);
		mm_idx_destroy(mi);
		

		return maxSelfOverlapLength;
	}
	*/
	/*
	void polishContigs(const string& outputContigFilename, const int minimapBatchSize){

		
		ContigPolisher contigPolisher(_inputDir, _readPartitionDir, outputContigFilename, _dataType, minimapBatchSize, true, _nbCores, _nbPartitions, _minimap2Preset_map);
		//readPartitionner.computeContigCoverages(_contigIndex_to_alignments, _contigSizes);
		contigPolisher.execute();

	}
	*/


	void derepContigs(const string& inputContigFilename, const string& outputContigFilename, const int minimapBatchSize){

		auto start = high_resolution_clock::now();

		ContigDerep contigDerep(inputContigFilename, outputContigFilename, minimapBatchSize, _inputDir, _readPartitionDir, 0.9, _minContigLength, _nbCores);
		contigDerep.execute();
		
		Logger::get().debug() << "\tDone " << " (" << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s) " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << "GB";

	}

	void trimContigs(const string& inputContigFilename, const string& outputContigFilename, const int minimapBatchSize){

		auto start = high_resolution_clock::now();

		ContigTrimmer contigTrimmer(inputContigFilename, _readPartitionDir + "/usedReads.fasta.gz", outputContigFilename, minimapBatchSize, _inputDir, _readPartitionDir, _minContigLength, _nbCores, _minimap2Preset_map);
		contigTrimmer.execute();
		
		Logger::get().debug() << "\tDone " << " (" << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s) " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << "GB";

	}

	/*
	void removeContigRepeats(){
		
		ContigRepeatRemover contigRepeatRemover(_inputDir, _nbPartitions, _minimizerSize, _minimizerDensity_assembly, _useHomopolymerCompression, _kminmerSizeFirst+1, _nbCores, _readPartitionDir, _minimap2Preset_map);
		
		contigRepeatRemover.execute();
	}
	*/
};	


#endif 

