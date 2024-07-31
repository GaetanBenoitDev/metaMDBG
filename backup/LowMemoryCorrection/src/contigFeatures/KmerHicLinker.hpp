

/*

./bin/metaMDBG link ~/appa/run/metaMDBG/plankton_1dataset/singleContigs/__results/mag_singleContigs/__checkmCircularContigs/bins/bin.circ.83.fa ~/appa/data/plankton/hicAll_1.fq.gz ~/appa/data/plankton/hicAll_2.fq.gz ~/appa/tmp/linkedReads.fastq.gz ~/appa/data/plankton/Sample1.fastq -t 16
./bin/metaMDBG link ~/appa/data/genomes/micromonas/assembly.fasta ~/appa/data/plankton/hicAll_1.fq.gz ~/appa/data/plankton/hicAll_2.fq.gz ~/appa/tmp/linkedReads_micromonas.fastq ~/appa/data/plankton/Sample1.fastq.gz -t 64
*/

#ifndef MDBG_METAG_KmerHicLinker
#define MDBG_METAG_KmerHicLinker

#include "../Commons.hpp"

typedef u_int128_t KmerType;
typedef phmap::parallel_flat_hash_map<KmerType, vector<u_int32_t>, phmap::priv::hash_default_hash<KmerType>, phmap::priv::hash_default_eq<KmerType>, std::allocator<std::pair<KmerType, vector<u_int32_t>>>, 4, std::mutex> KmerCountMap5;
typedef phmap::parallel_flat_hash_map<KmerType, u_int32_t, phmap::priv::hash_default_hash<KmerType>, phmap::priv::hash_default_eq<KmerType>, std::allocator<std::pair<KmerType, u_int32_t>>, 4, std::mutex> KmerCountMap6;
typedef phmap::parallel_flat_hash_set<KmerType, phmap::priv::hash_default_hash<KmerType>, phmap::priv::hash_default_eq<KmerType>, std::allocator<KmerType>, 4, std::mutex> ReferenceKmerSet;



class KmerHicLinker : public Tool{
    
public:

	string _inputContigFilename;
	string _readHic1;
	string _readHic2;
	string _outputFilename;
	int _nbCores;
	size_t _kmerSize;
	float _kmerDensity;
	string _readFilename;
	ofstream _outputFile;

	KmerCountMap5 _readKmers;
	//KmerCountMap6 _kmerCounts2;
	unordered_map<u_int32_t, u_int32_t> _contigIndex_to_nbHits;
	//SparseMatrix::SparseMatrix<u_int32_t>* _matrix;
	ofstream _contigDataFile;

	ReferenceKmerSet _referenceKmers;

	KmerHicLinker(): Tool (){

		//_matrix = new SparseMatrix::SparseMatrix<u_int32_t>(100000); // 3×3 matrix of integers
	}

	void parseArgs(int argc, char* argv[]){
		args::ArgumentParser parser("toBasespace", ""); //"This is a test program.", "This goes after the options."
		args::Positional<std::string> arg_inputContigFilename(parser, "inputContigFilename", "", args::Options::Required);
		args::Positional<std::string> arg_readHic1(parser, "readHic1", "", args::Options::Required);
		args::Positional<std::string> arg_readHic2(parser, "readHic2", "", args::Options::Required);
		args::Positional<std::string> arg_outputFilename(parser, "outputFilename", "", args::Options::Required);
		args::PositionalList<std::string> arg_readFilenames(parser, "reads", "Read filename(s) (separated by space)", args::Options::Required);
		//args::Positional<std::string> arg_contigs(parser, "contigs", "", args::Options::Required);
		//args::PositionalList<std::string> arg_readFilenames(parser, "reads", "Input filename(s) (separated by space)", args::Options::Required);
		args::ValueFlag<int> arg_kmerSize(parser, "", "Minimizer length", {ARG_MINIMIZER_LENGTH2}, 31);
		args::ValueFlag<float> arg_kmerDensity(parser, "", "Minimizer density", {ARG_MINIMIZER_DENSITY2}, 0.1);
		//args::ValueFlag<std::string> arg_contigs(parser, "", "", {ARG_INPUT_FILENAME_CONTIG}, "");
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

		_inputContigFilename = args::get(arg_inputContigFilename);
		_readHic1 = args::get(arg_readHic1);
		_readHic2 = args::get(arg_readHic2);
		_outputFilename = args::get(arg_outputFilename);
		_kmerSize = args::get(arg_kmerSize);
		_kmerDensity = args::get(arg_kmerDensity);
		_nbCores = args::get(arg_nbCores);

		_minHicLink = 4;

		_readFilename = _outputFilename + ".input.tmp.txt";
		Commons::createInputFile(args::get(arg_readFilenames), _readFilename);

	}

    void execute (){
    

		cout << "Indexing contigs" << endl;
		extractContigKmers();

		cout << "Collecting Hic links" << endl;
		collectHicLinks();

		cout << "Extracting linked reads" << endl;
		extractLinkedReads();

	}

	void extractContigKmers (){

		//ReadParserParallel readParser2(_inputContigFilename, true, false, _nbCores, _logFile);
		//readParser._maxReads = 10000;
		//readParser2.parse(ContigFunctorCount(*this));
		
		ReadParserParallel readParser(_inputContigFilename, true, false, _nbCores, _logFile);
		//readParser._maxReads = 10000;
		readParser.parse(ContigFunctor(*this));

		ReadParserParallel readParser2(_readFilename, false, false, _nbCores, _logFile);
		//readParser._maxReads = 100000;
		readParser2.parse(ReadFunctor(*this));

		cout << "Nb kmers ref: " << _referenceKmers.size() << endl;
		cout << "Nb kmers reads: " << _readKmers.size() << endl;

	}

	class ReadFunctor {

		public:

		KmerHicLinker& _kmerCounter;
		MinimizerParser128* _minimizerParser;
		//KmerModel128* _kmerModel;
		u_int64_t _currentDatasetIndex;

		ReadFunctor(KmerHicLinker& kmerCounter) : _kmerCounter(kmerCounter){
			
		}

		ReadFunctor(const ReadFunctor& copy) : _kmerCounter(copy._kmerCounter){
			_minimizerParser = new MinimizerParser128(_kmerCounter._kmerSize, _kmerCounter._kmerDensity);
		}

		~ReadFunctor(){
		}

		void operator () (const Read& read) {


			u_int32_t readIndex = read._index;
			if(readIndex % 100000 == 0) cout << readIndex << endl;
			
			vector<KmerType> minimizers;
			vector<u_int64_t> minimizers_pos;
			_minimizerParser->parse(read._seq, minimizers, minimizers_pos);


			//#pragma omp critical
			//{
			for(KmerType minimizer : minimizers){
				//if(_kmerCounter._kmerCounts2[minimizer] <= 1) continue;
				
				/*
				_kmerCounter._kmerCounts.try_emplace_l(minimizer, 
				[](KmerCountMap::value_type& v) { // key exist
					//v.second += 1;
				}, 0);
				*/
				
				_kmerCounter._readKmers.lazy_emplace_l(minimizer, 
				[&readIndex](KmerCountMap5::value_type& v) { // key exist
					if(std::find(v.second.begin(), v.second.end(), readIndex) == v.second.end()){
						v.second.push_back(readIndex);
					}
				},           
				[&minimizer, &readIndex](const KmerCountMap5::constructor& ctor) { // key inserted
					vector<u_int32_t> v = {readIndex};
					ctor(minimizer, v);
				});
				

				//if(_kmerCounter._kmerCounts.find(minimizer) == _kmerCounter._kmerCounts.end()) _kmerCounter._kmerCounts[minimizer] = 0;
			}

			//#pragma omp critical
			//{
				//_kmerCounter._contigDataFile << (read._index+1) << "\t" << read._header << "\t" << read._seq.size() << "\t" << "50" << "\t" << "50" << "\t" << "-" << "\t" << "-" << endl;
			//}
		
		}

	};
	
	class ContigFunctor {

		public:

		KmerHicLinker& _kmerCounter;
		MinimizerParser128* _minimizerParser;
		//KmerModel128* _kmerModel;
		u_int64_t _currentDatasetIndex;

		ContigFunctor(KmerHicLinker& kmerCounter) : _kmerCounter(kmerCounter){
			
		}

		ContigFunctor(const ContigFunctor& copy) : _kmerCounter(copy._kmerCounter){
			_minimizerParser = new MinimizerParser128(_kmerCounter._kmerSize, _kmerCounter._kmerDensity);
		}

		~ContigFunctor(){
		}

		void operator () (const Read& read) {


			u_int32_t readIndex = read._index;
			if(readIndex % 1000 == 0) cout << readIndex << endl;
			
			vector<KmerType> minimizers;
			vector<u_int64_t> minimizers_pos;
			_minimizerParser->parse(read._seq, minimizers, minimizers_pos);


			//#pragma omp critical
			//{
			for(KmerType minimizer : minimizers){
				//if(_kmerCounter._kmerCounts2[minimizer] <= 1) continue;
				
				/*
				_kmerCounter._kmerCounts.try_emplace_l(minimizer, 
				[](KmerCountMap::value_type& v) { // key exist
					//v.second += 1;
				}, 0);
				*/
				
				_kmerCounter._referenceKmers.lazy_emplace_l(minimizer, 
				[&readIndex](ReferenceKmerSet::value_type& v) { // key exist
					//if(std::find(v.second.begin(), v.second.end(), readIndex) == v.second.end()){
					//	v.second.push_back(readIndex);
					//}
				},           
				[&minimizer, &readIndex](const ReferenceKmerSet::constructor& ctor) { // key inserted
					//vector<u_int32_t> v = {readIndex};
					ctor(minimizer);
				});
				

				//if(_kmerCounter._kmerCounts.find(minimizer) == _kmerCounter._kmerCounts.end()) _kmerCounter._kmerCounts[minimizer] = 0;
			}

			//#pragma omp critical
			//{
				//_kmerCounter._contigDataFile << (read._index+1) << "\t" << read._header << "\t" << read._seq.size() << "\t" << "50" << "\t" << "50" << "\t" << "-" << "\t" << "-" << endl;
			//}
		
		}

	};

	unordered_map<u_int32_t, u_int32_t> _readIndex_to_nbLinks;

	
	void collectHicLinks (){


		ReadParserParallelPaired readParser(_readHic1, _readHic2, _nbCores, _logFile);
		//readParser._maxReads = 1000000;
		readParser.parse(HicReadFunctor(*this));

	
	}

	class HicReadFunctor {

		public:

		KmerHicLinker& _parent;
		MinimizerParser128* _minimizerParser;
		//KmerModel128* _kmerModel;

		HicReadFunctor(KmerHicLinker& parent) : _parent(parent){
			
		}

		HicReadFunctor(const HicReadFunctor& copy) : _parent(copy._parent){
			//_kmerModel = new KmerModel128(_kmerCounter._kmerSize);
			_minimizerParser = new MinimizerParser128(_parent._kmerSize, _parent._kmerDensity);
		}

		~HicReadFunctor(){
		}

		bool isMatchingContig(const vector<KmerType>& hicReadMinimizers){

			u_int64_t nbMatches = 0;
			
			for(KmerType minimizer : hicReadMinimizers){

				if(_parent._referenceKmers.find(minimizer) != _parent._referenceKmers.end()){
					nbMatches += 1;
				}
			}

			return nbMatches > (hicReadMinimizers.size() * 0.5);
		}
		
		void operator () (const Read& read1, const Read& read2) {


			u_int64_t readIndex = read1._index;
			if(readIndex % 1000000 == 0) cout << readIndex << " " << _parent.getNbLinks() << endl;
			
			//string rleSequence;
			//vector<u_int64_t> rlePositions;
			//Encoder::encode_rle(read->seq.s, strlen(read->seq.s), rleSequence, rlePositions);

			vector<KmerType> minimizers1;
			vector<u_int64_t> minimizers_pos1;
			_minimizerParser->parse(read1._seq, minimizers1, minimizers_pos1);


			vector<KmerType> minimizers2;
			vector<u_int64_t> minimizers_pos2;
			_minimizerParser->parse(read2._seq, minimizers2, minimizers_pos2);
			

			bool isMatchContig = isMatchingContig(minimizers1) || isMatchingContig(minimizers2);

			if(!isMatchContig) return;

			


			unordered_map<u_int32_t, u_int32_t> nbKmerMatches1;
			unordered_map<u_int32_t, u_int32_t> nbKmerMatches2;

			for(KmerType minimizer : minimizers1){

				if(_parent._readKmers.find(minimizer) == _parent._readKmers.end()){
					continue;
				}

				for(u_int32_t readIndex : _parent._readKmers[minimizer]){
					nbKmerMatches1[readIndex] += 1;
				}


			}

			for(KmerType minimizer : minimizers2){


				if(_parent._readKmers.find(minimizer) == _parent._readKmers.end()){
					continue;
				}

				for(u_int32_t readIndex : _parent._readKmers[minimizer]){
					nbKmerMatches2[readIndex] += 1;
				}

			}

			unordered_set<u_int32_t> matchingReadIndexes;

			for(const auto& it : nbKmerMatches1){
				if(it.second >= minimizers1.size()*0.5){
					matchingReadIndexes.insert(it.first);
				}
			}
			for(const auto& it : nbKmerMatches2){
				if(it.second >= minimizers2.size()*0.5){
					matchingReadIndexes.insert(it.first);
				}
			}


			#pragma omp critical
			{
				
				for(u_int32_t readIndex : matchingReadIndexes){
					_parent._readIndex_to_nbLinks[readIndex] += 1;
				}
			}
		
		}

	};



	/*
	void dumpNetwork(){

		//ofstream networkFile(_outputFilename);
		ofstream networkFile_normalized(_outputFilename + ".norm.txt");

		for(auto& it : _nbContigLinks){

			//auto& it2 = it.second;
			//const DbgEdge& contigPair = it.first;
			//u_int32_t nbLinks = it.second;

			for(auto& it2: it.second){

				double nbHits1 = _contigIndex_to_nbHits[it.first];
				double nbHits2 = _contigIndex_to_nbHits[it2.first];

				//cout << it.first << " " << it2.first << " " << it2.second << endl;

				//cout << nbHits1 << " " << nbHits2 << " " << it2.second << endl;
				if(it2.second <= 4) continue;
				
				float factor = sqrt(nbHits1 * nbHits2);
				float weightNormalized = it2.second / factor;

				networkFile_normalized.write((const char*)&it.first, sizeof(it.first));
				networkFile_normalized.write((const char*)&it2.first, sizeof(it2.first));
				networkFile_normalized.write((const char*)&weightNormalized, sizeof(weightNormalized));
				
				//networkFile_normalized << (it.first+1) << "\t" << (it2.first+1) << "\t" << weightNormalized << endl;
				//networkFile_normalized << (it.first) << "\t" << (it2.first) << "\t" << weightNormalized << endl;
			}
		}

		//networkFile.close();
		networkFile_normalized.close();


	}
	*/

	u_int64_t getNbLinks(){

		u_int64_t nbLinkedReads = 0;

		for(const auto& it : _readIndex_to_nbLinks){
			if(it.second >= _minHicLink){
				nbLinkedReads += 1;
			}
		}

		return nbLinkedReads;
	}

	u_int64_t _minHicLink;

	void extractLinkedReads(){

		cout << "Nb linked reads: " << getNbLinks() << endl;

		_outputFile = ofstream(_outputFilename);

		ReadParserParallel readParser(_readFilename, false, false, 1, _logFile);
		readParser.parse(ReadExtractor(*this));

		_outputFile.close();

	}


	class ReadExtractor {

		public:

		KmerHicLinker& _parent;

		ReadExtractor(KmerHicLinker& parent) : _parent(parent){
			
		}

		ReadExtractor(const ReadExtractor& copy) : _parent(copy._parent){
		}

		~ReadExtractor(){
		}

		void operator () (const Read& read) {

			u_int32_t readIndex = read._index;
			if(readIndex % 100000 == 0) cout << readIndex << endl;
			
			if(_parent._readIndex_to_nbLinks.find(readIndex) == _parent._readIndex_to_nbLinks.end()) return;
			if(_parent._readIndex_to_nbLinks[readIndex] < _parent._minHicLink) return;

			_parent._outputFile << ">" << read._header << endl;
			_parent._outputFile << read._seq << endl;
			_parent._outputFile << "+" << endl;
			_parent._outputFile << read._qual << endl;
		}

	};

};	


#endif 



