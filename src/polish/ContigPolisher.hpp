
/*
- Polisher pas deterministe dans la selection des 20 copies
	- utiliser un meilleur critiere (qualité de l'alignement etc)

- memory usage: attention si on a trier les contigs par longueur, les premiere passes pourrait demander un max de mémoire
	- !! attention les contigs sont trier dans le sens inverse (les plus grand a la fin)
 	
-ToBasespace index reads:
	- pas besoin d'indexer les reads qui ne contiennent pas de minimizer dupliquer

- ToBasespaceNoCorrection a paralleliser (remove overlaps etc)

- Tester racon sur contigs 0.33 avec l'option -c (write cigar in ouput file), voir si ça ameliore les résulats de binning ou non

- read best hit: la maniere dont on le selection actuellement est pas forcement la meilleur car elle ne favorise pas les petit contigs (les grand contig vont attirer les reads meme si leur sililarité est plus faible que sur un petit contig)
- attention ne pas outupt empty contig
- ToBasespace: générer un assemblage grossier
- tester les scores de qualité ?
- quand un read a été process, on peut erase son entrée dans _alignments ?
- window sequence a selectionner en priorité: 
	- la distance est une mauvaise metrique car on ne sait pas si la sequence de reference est erroné ou non
	- un mapping tres long sur un contig a plus de valeur qu'un mapping court (on est plus sûr que ce read appartient au contig)+
	- NEW: on peut garder la window avec la quality la plus haute
- si "no sequencer fo window": utiliser la fenetre original du contig ?
- minimap2 output: écrire un petit programme pour compresser les résultats d'alignements on the fly
- version finale: remove le minimap2 align filename
- voir comment racon handle multimap
- minimap2 read name: lors d'une premiere passe sur les reads, ecrire tous les headers (minimap2 name) dans un fichier
- vu que le polisher est standalone, ne pas changer le nom des header des contigs (les load en memoire t restituer dans le fichier output)
- mettre minimap2 dans le process de contig subsampling (donc écrire un chunk de contig dans un fichier séparer, puis refaire le minimap2 a cahque pass etc)

"reprise: ajouter les dernieres filtres de qualité de racon (easy) et refaire des test, aussi potentiellement 20 copies est pas l'optimum, essayer entre 20 et 25, tester correction k=5, revoir comment calculer superbubble lineairement"
- Blocoo creat graph: load les edge apres que les nodes du graphe soit tous indexé, comme ça on peut delete mdbg_init (abundance calculation) avant d'indexé les edges
*/

#ifndef MDBG_METAG_CONTIGPOLISHER
#define MDBG_METAG_CONTIGPOLISHER

#include "../Commons.hpp"
#include "../utils/edlib.h"
#include "../utils/spoa/include/spoa/spoa.hpp"
#include "../utils/DnaBitset.hpp"
#include "../utils/abPOA2/include/abpoa.h"


/*
extern unsigned char nt4_table;
unsigned char nt4_table2[256] = {
       0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

// 65,97=>A, 67,99=>C, 71,103=>G, 84,85,116,117=>T, else=>N
const char nt256_table2[256] = {
       'A', 'C', 'G', 'T',  'N', '-', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', '-',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'A', 'N', 'C',  'N', 'N', 'N', 'G',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'T', 'T', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'A', 'N', 'C',  'N', 'N', 'N', 'G',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'T', 'T', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N'
};
*/

class ContigPolisher : public Tool{
    
public:

	string _inputFilename_reads;
	string _inputFilename_contigs;
	int _nbCores;
	size_t _windowLength;
	size_t _maxWindowCopies;
	string _mapperOutputExeFilename;
	bool _useQual;
	string _outputDir;
	string _tmpDir;
	
	string _outputFilename_contigs;
	string _outputFilename_mapping;
	u_int64_t _maxMemory;
	double _qualityThreshold;

	//abpoa_para_t *abpt;
	
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
	};

	ContigPolisher(): Tool (){

	}


	void parseArgs(int argc, char* argv[]){

		string ARG_USE_QUAL = "qual";
		string filenameExe = argv[0];
		//cout << filenameExe << endl;

		fs::path pa(filenameExe);
		_mapperOutputExeFilename = pa.parent_path().string() + "/mapper";
		//cout << _mapperOutputExeFilename << endl;
		//exit(1);

		cxxopts::Options options("ToBasespace", "");
		options.add_options()
		("contigs", "", cxxopts::value<string>())
		("reads", "", cxxopts::value<string>())
		("tmpDir", "", cxxopts::value<string>())
		(ARG_USE_QUAL, "", cxxopts::value<bool>()->default_value("false"))
		(ARG_NB_CORES, "", cxxopts::value<int>()->default_value("4"));

		options.parse_positional({"contigs", "reads", "tmpDir"});
		options.positional_help("contigs reads tmpDir");


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

			_inputFilename_reads = result["reads"].as<string>();
			_inputFilename_contigs = result["contigs"].as<string>();
			_outputDir = result["tmpDir"].as<string>();;
			_nbCores = result[ARG_NB_CORES].as<int>();
			_useQual = result[ARG_USE_QUAL].as<bool>();
			_windowLength = 500;
			_maxWindowCopies = 21; //21;
			_qualityThreshold = 10.0;
			
		}
		catch (const std::exception& e){
			std::cout << options.help() << std::endl;
			std::cerr << e.what() << std::endl;
			std::exit(EXIT_FAILURE);
		}

		cout << "Contigs: " << _inputFilename_contigs << endl;
		cout << "Reads: " << _inputFilename_reads << endl;
		cout << "Use quality: " << _useQual << endl;

		fs::path p(_inputFilename_contigs);
		while(p.has_extension()){
			p.replace_extension("");
		}

		_tmpDir = _outputDir + "/__tmp/";
		if(!fs::exists(_tmpDir)){
			fs::create_directories(_tmpDir);
		}

		_outputFilename_contigs = p.string() + "_corrected.fasta.gz";
		_outputFilename_mapping = p.string() + "_tmp_mapping__.paf";
		_maxMemory = 4000000000ull;

		if(_useQual){
			_windowByteSize = _maxWindowCopies  * ((_windowLength/4) + _windowLength);
		}
		else{
			_windowByteSize = _maxWindowCopies  * (_windowLength/4);
		}
	}


	u_int64_t _currentContigSize;
	u_int64_t _windowByteSize;
	//unordered_set<u_int32_t> _validContigIndexes;
	gzFile _outputContigFile;

    void execute (){


		_outputContigFile = gzopen(_outputFilename_contigs.c_str(),"wb");

		/*
		abpt = abpoa_init_para();
		abpt->out_msa = 0; // generate Row-Column multiple sequence alignment(RC-MSA), set 0 to disable
		abpt->out_cons = 1; // generate consensus sequence, set 0 to disable
		abpt->w = 6, abpt->k = 9; abpt->min_w = 10; // minimizer-based seeding and partition
		abpt->progressive_poa = 1;
		abpt->max_n_cons = 1;

		abpoa_post_set_para(abpt);
		*/
		
		/*
		mapReads();

		
		parseAlignments(false);
		writeAlignmentBestHits();
		//parseAlignments(true, false);
		partitionReads();
		*/

		for(size_t i=0; i<10000000ull; i++){
			_contigToPartition[i] = 0;
		}

		_nbPartitions = 1;
		_partitionNbReads[0] = 10000;

		for(size_t i=0; i<_nbPartitions; i++){
			processPartition(i);
		}

		exit(1);
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
			//cout << totalByteSize << " " << _windowByteSize << endl;
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
		gzclose(_outputContigFile);
		fs::remove_all(_tmpDir);

	}


	/*
	u_int64_t loadContig(u_int64_t contigIndex, const string& seq){

		//_validContigIndexes.insert(contigIndex);

		u_int64_t totalByteSize = 0;

		_contigSequences[contigIndex] = seq;

		//cout << "load contig: " << contigIndex << " " << seq.size() << endl;
		size_t nbWindows = ceil((double)seq.size() / (double)_windowLength);
		vector<vector<Window>> windows(nbWindows);
		//cout << "Nb windows: " << nbWindows << endl;

		_contigWindowSequences[contigIndex] = windows;

		//cout << nbWindows << " " << nbWindows * _windowByteSize << endl;
		return nbWindows * _windowByteSize;
	}
	*/
	
	class PartitionFile{

		public:

		gzFile _file;
		omp_lock_t _mutex;

		PartitionFile(u_int32_t partition, const string& tmpDir){
			const string& partitionFilename = tmpDir + "/part_" + to_string(partition) + ".gz";
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

		cout << "Partitionning reads on the disk" << endl;
		collectContigStats();

        srand(time(NULL));
        std::random_shuffle(_contigStats.begin(), _contigStats.end());
		
		
		u_int64_t totalByteSize = 0;
		_nbPartitions = 0;

		for(const ContigStats& contigStat : _contigStats){

			totalByteSize += computeContigMemory(contigStat._length);
			if(totalByteSize >= _maxMemory){
				_nbPartitions += 1;
				totalByteSize = 0;
			}
			//_contigToPartition[contigIndex] = partition;

		}
		if(totalByteSize >= _maxMemory){
			_nbPartitions += 1;
		}

		_nbPartitions = max(_nbPartitions, (u_int32_t)_nbCores);
		_nbPartitions = min(_nbPartitions, (u_int32_t)_contigStats.size());
		u_int64_t nbContigsPerPartition = _contigStats.size() / _nbPartitions;
		cout << "Nb partitions: " << _nbPartitions << endl;
		cout << "Contigs per partition: " << nbContigsPerPartition << endl;

		for(u_int32_t i=0; i<_nbPartitions; i++){
			_partitionNbReads[i] = 0;
			_partitionFiles.push_back(new PartitionFile(i, _tmpDir));
		}

		u_int32_t partition = 0;
		u_int64_t nbContigs = 0;
		//totalByteSize = 0;

		for(const ContigStats& contigStat : _contigStats){

			//totalByteSize += computeContigMemory(contigStat._length);
			_contigToPartition[contigStat._contigIndex] = partition;
			//cout << nbContigs << " " << partition << endl;
			nbContigs += 1;

			//if(partition == 0){
			//	cout << "Contig in partition 0: " << contigStat._contigIndex << endl;
			//}
			
			if(nbContigs >= nbContigsPerPartition){
				//totalByteSize = 0;
				partition += 1;
				partition = min(partition, _nbPartitions-1);
				nbContigs = 0;
			}

			//_contigToPartition[contigIndex] = partition;

		}
		//if(totalByteSize >= _maxMemory){
		//	_partitionFiles.push_back(new PartitionFile(partition, _tmpDir));
		//}

		_contigStats.clear();


		//vector<vector<u_int32_t>> contigPartitions;
		//vector<u_int32_t> contigPartition;

		/*
		u_int64_t totalByteSize = 0;
		u_int64_t contigIndex = 0;
		u_int32_t partition = 0;

		gzFile fp;
		kseq_t *seq;
		int slen = 0, qlen = 0;
		fp = gzopen(_inputFilename_contigs.c_str(), "r");
		seq = kseq_init(fp);

		//"todo: nbPartition en fonction du maxMemory mais aussi nb thread available + besoin de randomiser l'ordre des contigs (créer des partitions equilibré"
		while (kseq_read(seq) >= 0){

			totalByteSize += computeContigMemory(strlen(seq->seq.s));
			_contigToPartition[contigIndex] = partition;

			//cout << totalByteSize << " " << _windowByteSize << endl;
			if(totalByteSize >= _maxMemory){

				//omp_lock_t writelock;
				//omp_init_lock(&writelock);

				_partitionFiles.push_back(new PartitionFile(partition, _tmpDir));
				//contigPartitions.push_back(contigPartition);
				//contigPartition.clear();
				totalByteSize = 0;
				partition += 1;
			}

			contigIndex += 1;
		}

		if(totalByteSize > 0){
			_partitionFiles.push_back(new PartitionFile(partition, _tmpDir));
		}


		kseq_destroy(seq);
		gzclose(fp);
		*/



		writeReadPartitions();

		//readToContigIndex.clear();
		//_contigToPartition.clear();
		_alignments.clear();

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

	u_int64_t computeContigMemory(size_t contigLength){
		size_t nbWindows = ceil((double)contigLength / (double)_windowLength);
		return nbWindows * _windowByteSize;
	}


	void writeReadPartitions(){
		ReadParserParallel readParser(_inputFilename_reads, false, false, _nbCores);
		readParser.parse(ReadPartitionFunctor(*this));
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

			u_int64_t readIndex = read._index;

			if(readIndex % 10000 == 0) cout << "\t" << readIndex << endl;
			
			if(_contigPolisher._alignments.find(readIndex) == _contigPolisher._alignments.end()) return;

			unordered_set<u_int32_t> writtenPartitions;

			for(const Alignment& al : _contigPolisher._alignments[readIndex]){
				u_int32_t contigIndex = al._contigIndex; //_contigPolisher._alignments[readIndex]._contigIndex;
				//cout << contigIndex << " " << (_contigPolisher._contigToPartition.find(contigIndex) != _contigPolisher._contigToPartition.end()) << endl;
				u_int32_t partition = _contigPolisher._contigToPartition[contigIndex];

				if(writtenPartitions.find(partition) != writtenPartitions.end()) continue;
				writtenPartitions.insert(partition);

				//cout << partition << endl;
				PartitionFile* partitionFile = _contigPolisher._partitionFiles[partition];

				omp_set_lock(&partitionFile->_mutex);
				
				//if(partition == 0){
					//cout << "Read: " << readIndex << endl;
				//}

				bool isFastq = read._qual.size() > 0 && _contigPolisher._useQual;

				_contigPolisher._partitionNbReads[partition] += 1;

				u_int32_t readSize = read._seq.size();

				if(isFastq){
					string header = '@' + to_string(read._index) + '\n';
					string seq = read._seq + '\n';
					gzwrite(partitionFile->_file, (const char*)&header[0], header.size());
					gzwrite(partitionFile->_file, (const char*)&seq[0], seq.size());
					static string strPlus = "+\n";
					gzwrite(partitionFile->_file, (const char*)&strPlus[0], strPlus.size());
					string qual = read._qual + '\n';
					gzwrite(partitionFile->_file, (const char*)&qual[0], qual.size());
				}
				else{
					string header = '>' + to_string(read._index) + '\n';
					string seq = read._seq + '\n';
					gzwrite(partitionFile->_file, (const char*)&header[0], header.size());
					gzwrite(partitionFile->_file, (const char*)&seq[0], seq.size());
				}


				//gzwrite(partitionFile->_file, (const char*)&readIndex, sizeof(readIndex));
				//gzwrite(partitionFile->_file, (const char*)&readSize, sizeof(readSize));
				//gzwrite(partitionFile->_file, read._seq.c_str(), read._seq.size());
				//gzwrite(partitionFile->_file, read._qual.c_str(), read._qual.size());
				omp_unset_lock(&partitionFile->_mutex);
			}
		}
	};

	u_int32_t _currentPartition;

	void processPartition(u_int32_t partition){
		_currentPartition = partition;
		//if(_contigSequences.size() == 0) return;
		cout << "Processing partition: " << _currentPartition << endl;

		if(_partitionNbReads[partition] == 0) return;


		//indexContigs();
		clearPass();
		loadContigs();
		parseAlignments(true);
		//_contigName_to_contigIndex.clear();
		//_readName_to_readIndex.clear();

		collectWindowSequences(_currentPartition);
		performCorrection();
		//if(fs::exists(_outputFilename_mapping)) fs::remove(_outputFilename_mapping);
	}

	void loadContigs(){
		auto fp = std::bind(&ContigPolisher::loadContigs_read, this, std::placeholders::_1);
		ReadParser readParser(_inputFilename_contigs, true, false);
		readParser.parse(fp);
	}
	
	void loadContigs_read(const Read& read){

		u_int32_t contigIndex = read._index;
		u_int32_t contigPartition = _contigToPartition[contigIndex];

		if(contigPartition != _currentPartition) return;
		

		//if(_currentPartition == 0) cout << "Loading contig in partition 0: " << contigIndex << endl;
		_contigSequences[contigIndex] = read._seq;
		size_t nbWindows = ceil((double)read._seq.size() / (double)_windowLength);
		vector<vector<Window>> windows(nbWindows);

		_contigWindowSequences[contigIndex] = windows;
	}

	//_validContigIndexes.insert(contigIndex);

	//u_int64_t totalByteSize = 0;


	//cout << "load contig: " << contigIndex << " " << seq.size() << endl;
	//cout << "Nb windows: " << nbWindows << endl;

	

	//cout << nbWindows << " " << nbWindows * _windowByteSize << endl;
	//return nbWindows * _windowByteSize;


	void clearPass(){
		_contigSequences.clear();
		//_validContigIndexes.clear();
		_contigWindowSequences.clear();
		//_alignments.clear();
	}
	/*
	void loadContigs(){
		auto fp = std::bind(&ContigPolisher::loadContigs_read, this, std::placeholders::_1);
		ReadParser readParser(_inputFilename_contigs, true, false);
		readParser.parse(fp);
	}
	*/

	void mapReads(){
		
		string readFilenames = "";
		ReadParser readParser(_inputFilename_reads, false, false);

		for(const string& filename : readParser._filenames){
			readFilenames += filename + " ";
		}

		string command = "minimap2 -H -I 2GB -t " + to_string(_nbCores) + " -x map-hifi " + _inputFilename_contigs + " " + readFilenames;
		command += " | " + _mapperOutputExeFilename + " " + _inputFilename_contigs + " " + _inputFilename_reads + " " + _outputFilename_mapping;
		Utils::executeCommand(command, _outputDir);

		//minimap2 -x map-hifi ~/workspace/run/overlap_test_201/contigs_47.fasta.gz ~/workspace/data/overlap_test/genome_201_50x/simulatedReads_0.fastq.gz | ./bin/mapper ~/workspace/run/overlap_test_201/contigs_47.fasta.gz ~/workspace/data/overlap_test/genome_201_50x/input.txt ~/workspace/run/overlap_test_201/align.bin
	}

	struct Alignment{
		u_int32_t _contigIndex;
		//u_int64_t _readIndex;
		bool _strand;
		u_int32_t _readStart;
		u_int32_t _readEnd;
		u_int64_t _contigStart;
		u_int64_t _contigEnd;
		//float _score;
		//u_int64_t _length;
		
		u_int64_t length(){
			return std::max((u_int64_t)(_readEnd - _readStart), (u_int64_t)(_contigEnd - _contigStart));
		}
		
	};

	struct AlignmentPartitionning{
		u_int32_t _contigIndex;
		u_int32_t _length;
	};

	//unordered_map<string, u_int32_t> _contigName_to_contigIndex;
	//unordered_map<string, u_int64_t> _readName_to_readIndex;
	//unordered_map<u_int64_t, AlignmentPartitionning> _readToContigIndex;
	unordered_map<u_int64_t, vector<Alignment>> _alignments;
	unordered_map<u_int32_t, string> _contigSequences;
	unordered_map<u_int32_t, vector<vector<Window>>> _contigWindowSequences;
	unordered_map<ContigRead, u_int32_t, ContigRead_hash> _alignmentCounts;
	//u_int64_t _correctedContigIndex;


	/*
	void indexContigName(){
		
		cout << "Indexing contig names" << endl;

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

		//cout << minimap_name << endl;
		_contigName_to_contigIndex[minimap_name] = read._index;
	}

	void indexReadName(){

		cout << "Indexing read names" << endl;
		auto fp = std::bind(&ContigPolisher::indexReadName_read, this, std::placeholders::_1);
		ReadParser readParser(_inputFilename_reads, false, false);
		readParser.parse(fp);

	}
	
	void indexReadName_read(const Read& read){
		//cout << read._index << endl;

		string minimap_name;

		auto find = read._header.find(' ');
		if(find == std::string::npos){
			minimap_name = read._header;
		}
		else{
		 	minimap_name = read._header.substr(0, find);
		}

		//string minimap_name = read._header.substr(0, read._header.find(' '));
		//cout << minimap_name << endl;
		_readName_to_readIndex[minimap_name] = read._index;
	}
	*/

	void parseAlignments(bool indexPartitionOnly){

		cout << "\tIndexing read alignments" << endl;

        ifstream infile(_outputFilename_mapping);

        std::string line;
        vector<string>* fields = new vector<string>();
        //vector<string>* fields_optional = new vector<string>();

		while(true){
        	/*
			//while (std::getline(infile, line)){

            //GfaParser::tokenize(line, fields, '\t');

			//cout << line << endl;

			const string& readName = (*fields)[0];
			const string& contigName = (*fields)[5];

			u_int32_t readStart = stoull((*fields)[2]);
			u_int32_t readEnd = stoull((*fields)[3]);
			u_int64_t contigStart = stoull((*fields)[7]);
			u_int64_t contigEnd = stoull((*fields)[8]);

			u_int64_t nbMatches = stoull((*fields)[9]);
			u_int64_t alignLength = stoull((*fields)[10]);
        	u_int64_t queryLength = stoull((*fields)[1]);

			bool strand = (*fields)[4] == "-";
			float score = (double) nbMatches / (double) queryLength;
			float score2 = (double) nbMatches / (double) alignLength;

			u_int32_t contigIndex = _contigName_to_contigIndex[contigName];
			u_int64_t readIndex = _readName_to_readIndex[readName];
			*/

			u_int32_t contigIndex;
			u_int64_t readIndex;
			u_int32_t readStart;
			u_int32_t readEnd;
			u_int64_t contigStart;
			u_int64_t contigEnd;
			float score;
			bool strand;

			//u_int64_t nbMatches;
			//u_int64_t alignLength;
			//u_int64_t queryLength;


			infile.read((char*)&contigIndex, sizeof(contigIndex));
			if(infile.eof()) break;
			infile.read((char*)&contigStart, sizeof(contigStart));
			infile.read((char*)&contigEnd, sizeof(contigEnd));
			infile.read((char*)&readIndex, sizeof(readIndex));
			infile.read((char*)&readStart, sizeof(readStart));
			infile.read((char*)&readEnd, sizeof(readEnd));
			infile.read((char*)&strand, sizeof(strand));
			infile.read((char*)&score, sizeof(score));

			if(indexPartitionOnly && _contigSequences.find(contigIndex) == _contigSequences.end()) continue;


			u_int32_t length = std::max((u_int64_t)(readEnd - readStart), (u_int64_t)(contigEnd - contigStart));
			Alignment align = {contigIndex, strand, readStart, readEnd, contigStart, contigEnd}; //, score

			if(_alignments.find(readIndex) == _alignments.end()){
				_alignments[readIndex].push_back(align);
			}
			else{

				//bool isBetter = false;
				for(Alignment& al: _alignments[readIndex]){
					if(al._contigIndex != contigIndex) continue;
					
					if(length > al.length()){

						al._contigIndex = contigIndex;
						al._strand = strand;
						al._readStart = readStart;
						al._readEnd = readEnd;
						al._contigStart = contigStart;
						al._contigEnd = contigEnd;

						//isBetter = true;
						break;
						//_alignments[readIndex] = align;
					}
				}

				//if(length > _alignments[readIndex].length()){
				//	_alignments[readIndex] = align;
				//}
				
			}

			/*
			if(_alignments.find(readIndex) == _alignments.end()){
				_alignments[readIndex] = align;
			}
			else{
				if(length > _alignments[readIndex].length()){
					_alignments[readIndex] = align;
				}
				
			}
			*/
			//ContigRead alignKey = {_contigName_to_contigIndex[contigName], _readName_to_readIndex[readName]};
			/*
			if(indexOnlyContigIndex){

				if(_readToContigIndex.find(readIndex) == _readToContigIndex.end()){
					_readToContigIndex[readIndex] = {contigIndex, length};
				}
				else{
					if(length > _readToContigIndex[readIndex]._length){
						_readToContigIndex[readIndex] = {contigIndex, length};
					}
					
				}

			}
			else{
				if(_alignments.find(readIndex) == _alignments.end()){
					_alignments[readIndex] = align;
				}
				else{
					if(length > _alignments[readIndex].length()){
						_alignments[readIndex] = align;
					}
					
				}

			}
			*/
			

        }

		/*
		if(!isPartionning){
			for(auto& it : _alignments){

				auto& als = it.second;
				vector<u_int64_t> removedIndex;

				if(als.size() == 1) continue;
				
				for(size_t i=0; i<als.size(); i++){
					
					if(std::find(removedIndex.begin(), removedIndex.end(), i) != removedIndex.end()) continue;

					for(size_t j=i+1; j<als.size(); j++){

						//if(als[i]._contigIndex != als[j]._contigIndex) continue;

						if(std::find(removedIndex.begin(), removedIndex.end(), j) != removedIndex.end()) continue;

						//cout << als[i]._score << " " << als[j]._score << endl;
						if(als[i].length() > als[j].length()){
							//if(als[i]._score < als[j]._score){
							removedIndex.push_back(j);
							//alsFiltered.push_back(als[i]);
						}
						else{
							removedIndex.push_back(i);
							//alsFiltered.push_back(als[j]);
						}
						
		
						//if(als[i]._contigIndex == 0 && als[i]._readIndex == 705){
						//	cout << endl;
						//	cout << i << " " << als[i]._contigIndex << " " << als[i]._readIndex << " " << als[i]._length << endl;
						//	cout << j << " " << als[j]._contigIndex << " " << als[j]._readIndex << " " << als[j]._length << endl;
						//	//getchar();
						//}
					}
				}
				
				vector<Alignment> alsFiltered;
				for(size_t i=0; i<als.size(); i++){
					if(std::find(removedIndex.begin(), removedIndex.end(), i) != removedIndex.end()) continue;
					alsFiltered.push_back(als[i]);
				}

				als = alsFiltered;
			}
		}
		*/

	}

	
	void writeAlignmentBestHits(){

        ofstream outputFile(_outputFilename_mapping);

		for(const auto& it : _alignments){

			u_int64_t readIndex = it.first;
			
			for(const Alignment& al : it.second){
				float score = 0;

				outputFile.write((const char*)&al._contigIndex, sizeof(al._contigIndex));
				outputFile.write((const char*)&al._contigStart, sizeof(al._contigStart));
				outputFile.write((const char*)&al._contigEnd, sizeof(al._contigEnd));
				outputFile.write((const char*)&readIndex, sizeof(readIndex));
				outputFile.write((const char*)&al._readStart, sizeof(al._readStart));
				outputFile.write((const char*)&al._readEnd, sizeof(al._readEnd));
				outputFile.write((const char*)&al._strand, sizeof(al._strand));
				outputFile.write((const char*)&score, sizeof(score));
			}
		}

		outputFile.close();
	}


	void collectWindowSequences(u_int32_t partition){
		
		cout << "\tCollecting window sequences" << endl;

		//auto fp = std::bind(&ContigPolisher::collectWindowCopies_read, this, std::placeholders::_1);
		//ReadParser readParser(_inputFilename_reads, false, false);
		//readParser.parse(fp);


		const string& partitionFilename = _tmpDir + "/part_" + to_string(partition) + ".gz";
		ReadParserParallel readParser(partitionFilename, true, false, _nbCores);
		readParser.parse(CollectWindowSequencesFunctor(*this));
	}
	
	class CollectWindowSequencesFunctor {

		public:

		ContigPolisher& _contigPolisher;
		//unordered_map<string, u_int32_t>& _contigName_to_contigIndex;
		//unordered_map<string, u_int64_t>& _readName_to_readIndex;
		unordered_map<u_int64_t, vector<Alignment>>& _alignments;
		unordered_map<u_int32_t, string>& _contigSequences;
		unordered_map<u_int32_t, vector<vector<Window>>>& _contigWindowSequences;
		size_t _windowLength;


		CollectWindowSequencesFunctor(ContigPolisher& contigPolisher) : _contigPolisher(contigPolisher), _alignments(contigPolisher._alignments), _contigSequences(contigPolisher._contigSequences), _contigWindowSequences(contigPolisher._contigWindowSequences), _windowLength(contigPolisher._windowLength){
		}

		CollectWindowSequencesFunctor(const CollectWindowSequencesFunctor& copy) : _contigPolisher(copy._contigPolisher), _alignments(copy._alignments), _contigSequences(copy._contigSequences), _contigWindowSequences(copy._contigWindowSequences), _windowLength(copy._windowLength){
			
		}

		~CollectWindowSequencesFunctor(){
		}

		void operator () (const Read& read) {

			u_int64_t readIndex = stoull(read._header);

			//if(_contigPolisher._currentPartition == 0) cout << readIndex << " " << (_alignments.find(readIndex) != _alignments.end()) << endl;
			
			//if(readIndex % 100000 == 0) cout << "\t" << readIndex << endl;

			if(_alignments.find(readIndex) == _alignments.end()) return;

			//const vector<Alignment>& als = _alignments[readIndex];
			//const Alignment& al = _alignments[readIndex];
			for(const Alignment& al : _alignments[readIndex]){
				u_int64_t contigIndex = al._contigIndex;

				if(_contigSequences.find(contigIndex) == _contigSequences.end()) continue;

				//cout << read._seq.size() << " " << read._qual.size() << " " << _contigSequences[contigIndex].size() << " " << al._readStart << " " << al._readEnd << " " << al._contigStart << " " << al._contigEnd << endl;
				string readSeq = read._seq;
				string qualSeq = read._qual;
				string readSequence = readSeq.substr(al._readStart, al._readEnd-al._readStart);
				string contigSequence = _contigSequences[contigIndex].substr(al._contigStart, al._contigEnd-al._contigStart);



				if(al._strand){
					Utils::toReverseComplement(readSequence);
					Utils::toReverseComplement(readSeq);
					std::reverse(qualSeq.begin(), qualSeq.end());
				}
				

				//cout << readSequence << endl;
				//cout << contigSequence << endl;

				//cout << contigSequence.size() << " "<< readSequence.size() << endl;
				static EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0);


				EdlibAlignResult result = edlibAlign(readSequence.c_str(), readSequence.size(), contigSequence.c_str(), contigSequence.size(), config);


				char* cigar;

				if (result.status == EDLIB_STATUS_OK) {
					cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
				} else {
					cout << "Invalid edlib results" << endl;
					exit(1);
				}

				//cout << cigar << endl;

				edlibFreeAlignResult(result);
				
				find_breaking_points_from_cigar(_windowLength, al, readSeq.size(), cigar, readSeq, qualSeq);
				free(cigar);

				//getchar();

			}

		}

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
				uint32_t window_start = (breaking_points_[j].first / _windowLength) *
					_windowLength;

				//const char* data = overlaps[i]->strand() ?
				//	&(sequence->reverse_complement()[breaking_points[j].second]) :
				//	&(sequence->data()[breaking_points[j].second]);
				const char* data = &readSequence[breaking_points_[j].second];
				uint32_t data_length = breaking_points_[j + 1].second - breaking_points_[j].second;

				string sequence = string(data, data_length);

				/*
				if(al._contigIndex == 1 && window_id > 1 && window_id < 1000){
					if(window_id == 192 || window_id == 193){
						cout << readSequence << endl;
						cout << sequence << endl;
						if(breaking_points_[j].second > 1) cout << readSequence[breaking_points_[j].second-1] << " " << readSequence[breaking_points_[j].second] << endl;
						if(readSequence.size() > readSequence[breaking_points_[j].second+data_length-1]) cout << readSequence[breaking_points_[j].second+data_length-1] << " " << readSequence[breaking_points_[j].second+data_length] << endl;
						cout << window_id << endl;
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
						cout << readSequence << endl;
						cout << sequence << endl;
						if(breaking_points_[j].second > 1){
							cout << readSequence[breaking_points_[j].second-1] << " " << readSequence[breaking_points_[j].second] << endl;
							if(readSequence[breaking_points_[j].second-1] == readSequence[breaking_points_[j].second]) getchar();
						}
						if(readSequence.size() > readSequence[breaking_points_[j].second+data_length-1]) cout << readSequence[breaking_points_[j].second+data_length-1] << " " << readSequence[breaking_points_[j].second+data_length] << endl;
						cout << window_id << endl;
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
					
					//cout << al._contigIndex << " " << al._readIndex << endl;
					//cout << (breaking_points_[j + 1].second - breaking_points_[j].second) << " " << (0.02 * _windowLength) << endl;
					//cout << sequence << endl;
					//cout << quality << endl;
					//getchar();
				//}

				//if(sequence.size() < 490) continue;
				indexWindow(al, window_id, posStart, posEnd, sequence, quality);

				//cout << window_id << " " << posStart << " " << posEnd << endl;
				//cout << sequence << endl;
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

				//cout << contigWindowStart << " " << contigWindowEnd  << "      " << readWindowStart << " " << readWindowEnd << endl;

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
			//	cout << breakPoint.first << " " << breakPoint.second << endl;
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
					//cout << contigWindowStart << endl;
					cout << windowSequence << endl;
					//cout << windows.size() << endl;
					cout << readWindowStart << " " << readWindowEnd << endl;
					getchar();
					//cout << windowQualities << endl;
				}
				*/
				
				bool interrupt = false;
				if(windows.size() < (_contigPolisher._maxWindowCopies-1)){


					windows.push_back({new DnaBitset2(windowSequence), windowQualities, posStart, posEnd});

					/*
					if(al._contigIndex == 1 && windowIndex == 2){
						cout << "1111111111111111" << endl;
						cout  << windowSequence << endl;
					}
					if(al._contigIndex == 1 && windowIndex == 3){
						cout << "2222222222222222" << endl;
						cout  << windowSequence << endl;
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
							cout << strlen(cigar) << endl;
							free(cigar);
						} else {
							cout << "Invalid edlib results" << endl;
							exit(1);
						}

						//cout << cigar << endl;

						edlibFreeAlignResult(result);
						free(dnaStrModel);

					}
					*/


					interrupt = true;
				}

				
				if(!interrupt){

					//u_int64_t largestAligmenet = 0;
					//u_int64_t largestAligmenetIndex = 0;

					
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
					

				} 




			}
			//cout << contigWindowIndex << endl;


			

			/*
			if(contigWindowIndex == _contigWindowSequences[al._contigIndex].size()-1){

				for(size_t i=0; i<windowSequences.size(); i++){
					cout << windowSequences[i]->m_len << " ";
				}
				cout << endl;
				//cout << contigWindowEnd-contigWindowStart << endl;
				//cout << windowSequence << endl;
			}
			*/
		}


	};

	void performCorrection(){
		
		u_int64_t checksum = 0;

		cout << "\tPerform correction" << endl;


		vector<u_int32_t> contigIndexes;
		for(const auto& it : _contigWindowSequences){
			contigIndexes.push_back(it.first);
		}
		std::sort(contigIndexes.begin(), contigIndexes.end());

		for(u_int32_t contigIndex : contigIndexes){


			vector<vector<Window>>& windows = _contigWindowSequences[contigIndex];
			u_int64_t nbCorrectedWindows = 0;
			vector<DnaBitset2*> correctedWindows(windows.size());
		

			#pragma omp parallel num_threads(_nbCores)
			{

				std::unique_ptr<spoa::AlignmentEngine> alignmentEngine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -4);
				
				#pragma omp for
				for(size_t w=0; w<windows.size(); w++){

					vector<Window>& sequences = windows[w];
					
					u_int64_t wStart = w*_windowLength;
					u_int64_t wEnd = min(_contigSequences[contigIndex].size(), wStart+_windowLength);
					string contigOriginalSequence = _contigSequences[contigIndex].substr(wStart, wEnd-wStart);

					checksum += sequences.size()+1;
					//cout << w << ": " << (sequences.size()+1) << endl;
					//continue;
					//cout << contigOriginalSequence << endl;
					//getchar();
					if(sequences.size() < 2){
							
						
						correctedWindows[w] = new DnaBitset2(contigOriginalSequence);
						//cout << "No sequences for window" << endl;
						continue;
					}
					
					#pragma omp atomic
					nbCorrectedWindows += 1;

					//std::sort(sequences.begin(), sequences.end(), []
					//(Window& first, Window& second){
					//	return first._sequence->m_len > second._sequence->m_len;
					//});

					//sequences.insert(sequences.begin(), new DnaBitset2(contigOriginalSequence));

				
					//vector<u_int32_t> windowLengths;
					//for(size_t i=0; i<sequences.size(); i++){
					//	windowLengths.push_back(sequences[i]->m_len);
					//}
					//cout << Utils::compute_median(windowLengths) << endl;


					//abpoa_t *ab = abpoa_init();

					//vector<size_t> order;
					//for(size_t i=0; i<sequences.size(); i++){
					//	order.push_back(i);
					//}
					//srand(time(NULL));
					//std::random_shuffle(order.begin(), order.end());

					/*
					vector<string> seqSorted;
					for(size_t i=1; i<sequences.size(); i++){
						DnaBitset2* dna = sequences[i];
						//const DnaBitset2* dna = variant._sequence; //sequenceCopies[s._sequenceIndex];
						char* dnaStr = dna->to_string();
						seqSorted.push_back(string(dnaStr));
						free(dnaStr);
					}


					std::sort(seqSorted.begin(), seqSorted.end());

					DnaBitset2* dnaModel = sequences[0];
					//const DnaBitset2* dna = variant._sequence; //sequenceCopies[s._sequenceIndex];
					char* dnaStrModel = dnaModel->to_string();
					seqSorted.insert(seqSorted.begin(), string(dnaStrModel));
					free(dnaStrModel);
					*/

					/*
					//cout << "1" << endl;
					int n_seqs = sequences.size();
					int *seq_lens = (int*)malloc(sizeof(int) * n_seqs);
					uint8_t **bseqs = (uint8_t**)malloc(sizeof(uint8_t*) * n_seqs);
					
					cout << "-----" << endl;
					cout << w << endl;
					//cout << contigOriginalSequence << endl;

					for(size_t i=0; i<sequences.size(); i++){ 

						//size_t i = order[ii];
						DnaBitset2* dna = sequences[i];
						//const DnaBitset2* dna = variant._sequence; //sequenceCopies[s._sequenceIndex];
						char* dnaStr = dna->to_string();
						//cout << dnaStr << endl;

						seq_lens[i] = strlen(dnaStr);
						bseqs[i] = (uint8_t*)malloc(sizeof(uint8_t) * seq_lens[i]);

						for (int j = 0; j < seq_lens[i]; ++j){
							bseqs[i][j] = nt4_table2[(int)(dnaStr[j])];
							//cout << dnaStr[j] << " " << bseqs[i][j] << endl;
						}

						//cout << s._sequenceIndex << " " << s._editDistance << endl;

						//auto alignment = _alignementEngine->Align(dnaStr, dna->m_len, _graph);
						//_graph.AddAlignment(alignment, dnaStr, dna->m_len);
						//sequencesStr.push_back()
						free(dnaStr);

						//processedSequences += 1;
					}



					//getchar();
					abpoa_msa(ab, abpt, n_seqs, NULL, seq_lens, bseqs, NULL);
					abpoa_cons_t *abc = ab->abc;



					string correctedSequence = "";

					for (size_t i = 0; i < abc->n_cons; ++i) {
						for (size_t j = 0; j < abc->cons_len[i]; ++j){
							correctedSequence += nt256_table2[abc->cons_base[i][j]];
						}
					}


					for (int i = 0; i < n_seqs; ++i) free(bseqs[i]); free(bseqs); free(seq_lens);
					abpoa_free(ab);
					*/

					/*
					std::vector<uint32_t> rank;
					rank.reserve(sequences.size());
					for (uint32_t i = 0; i < sequences_.size(); ++i) {
						rank.emplace_back(i);
					}
					*/


					std::sort(sequences.begin(), sequences.end(), [&](const Window& lhs, const Window& rhs) {
						return lhs._posStart < rhs._posStart; });

					//cout << sequences.size() << endl;
					//cout << "1" << endl;
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


						/*
						std::vector<uint32_t> rank;
						rank.reserve(sequences.size());
						for (uint32_t i = 0; i < sequences_.size(); ++i) {
							rank.emplace_back(i);
						}

						std::sort(rank.begin() + 1, rank.end(), [&](uint32_t lhs, uint32_t rhs) {
							return positions_[lhs].first < positions_[rhs].first; });
						*/

						//uint32_t offset = 0.01 * sequences[0].size();
						//for (uint32_t i = 1; i < sequences.size(); ++i) {
							//uint32_t i = rank[j];

							
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
							
							//cout << dnaStr << endl;
							//cout << window._quality << endl;
							//cout << window._posStart << " " << window._posEnd << endl;
							//getchar();
							//auto alignment = alignmentEngine->Align(dnaStr, window._sequence->m_len, graph);

							//cout << dnaStr << endl;
							//getchar();
							//fprintf(stdout, "%s\n", std::string(sequences_[i].first, sequences_[i].second).c_str());
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
						//}
					}


					//cout << "2" << endl;



					/*

					vector<vector<u_int32_t>> counts(abc->msa_len, vector<u_int32_t>(4, 0));

					for (int i = 0; i < abc->n_seq; ++i) {
						for (int j = 0; j < abc->msa_len; ++j) {

							char c = nt256_table2[abc->msa_base[i][j]];
							if(c == 'A'){
								counts[j][0] += 1;
							}
							else if(c == 'C'){
								counts[j][1] += 1;
							}
							else if(c == 'G'){
								counts[j][2] += 1;
							}
							else if(c == 'T'){
								counts[j][3] += 1;
							}

							//cout << (int) nt256_table2[abc->msa_base[i][j]];
							//fprintf(stdout, "%c", nt256_table2[abc->msa_base[i][j]]);
						}
						//cout << endl;
						//fprintf(stdout, "\n");
					}

					float t = n_seqs * 0.5;

					string correctedSequence;

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

					for (int i = 0; i < n_seqs; ++i) free(bseqs[i]); free(bseqs); free(seq_lens);
					abpoa_free(ab);
					*/

					//if(correctedSequence.size() < _windowLength-20){
					//	cout << correctedSequence.size() << " " << sequences.size() << " " << contigOriginalSequence.size() << endl;
					//}


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

					//cout << "3" << endl;
					//cout << correctedSequence << endl;
					//getchar();
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

			//cout << contigSequence.size() << " " << nbCorrectedWindows << endl;
			
			string header = ">ctg" + to_string(contigIndex) + '\n';
			gzwrite(_outputContigFile, (const char*)&header[0], header.size());
			contigSequence +=  '\n';
			gzwrite(_outputContigFile, (const char*)&contigSequence[0], contigSequence.size());
			//cout << contigSequence.size() << endl;
		}

		cout << "Checksum: " << checksum << endl;
	}

};	


#endif 


