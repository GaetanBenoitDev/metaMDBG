
//time ./bin/mdbgAsmMeta countKmer -t 8 -o ~/workspace/run/overlap_test_multik_AD/contigCoverages_k81_k61_all.tsv -i ~/workspace/data/AD/shortreads/input.txt -c ~/workspace/run/overlap_test_multik_AD/contigs_81.fasta.gz  -k 61


#ifndef MDBG_METAG_KMERCOUNTER
#define MDBG_METAG_KMERCOUNTER

#include "../Commons.hpp"

typedef u_int128_t KmerType;
typedef phmap::parallel_flat_hash_map<KmerType, u_int32_t, phmap::priv::hash_default_hash<KmerType>, phmap::priv::hash_default_eq<KmerType>, std::allocator<std::pair<KmerType, u_int32_t>>, 4, std::mutex> KmerCountMap;


class KmerCounter : public Tool{
    
public:

	string _inputFilename_contig;
	string _inputFilename_shortreads;
	string _outputFilename;
	float _trimmedOutliers;
	int _nbCores;
	size_t _kmerSize;
	float _kmerDensity;

	//string _tmpDir;
	//string _outputFilename_kmerCouts;
	vector<float> _countsInit;
	size_t _nbDatasets;

	//KmerModel* _kmerModel;
	KmerCountMap _kmerCounts;
	vector<vector<float>> _contigCoverages_mean;
	vector<vector<float>> _contigCoverages_var;
	u_int32_t _currentDatasetIndex;

	//MinimizerParser* _minimizerParser;
	u_int64_t _nbContigs;

	KmerCounter(): Tool (){

	}

	void parseArgs(int argc, char* argv[]){

		string ARG_TRIMMED_OUTLIERS = "to";

		cxxopts::Options options("ToMinspace", "");
		options.add_options()
		//(ARG_INPUT_FILENAME, "", cxxopts::value<string>())
		(ARG_OUTPUT_DIR, "", cxxopts::value<string>())
		(ARG_INPUT_FILENAME, "", cxxopts::value<string>()->default_value(""))
		(ARG_INPUT_FILENAME_CONTIG, "", cxxopts::value<string>()->default_value(""))
		(ARG_KMINMER_LENGTH, "", cxxopts::value<int>()->default_value("31"))
		(ARG_MINIMIZER_DENSITY, "", cxxopts::value<float>()->default_value("0.1"))
		(ARG_TRIMMED_OUTLIERS, "", cxxopts::value<float>()->default_value("0.05"))
		(ARG_NB_CORES, "", cxxopts::value<int>()->default_value(NB_CORES_DEFAULT));


		if(argc <= 1){
			cout << options.help() << endl;
			exit(0);
		}

		cxxopts::ParseResult result;

		try{
			result = options.parse(argc, argv);

			_inputFilename_shortreads = result[ARG_INPUT_FILENAME].as<string>();
			_outputFilename = result[ARG_OUTPUT_DIR].as<string>();
			_inputFilename_contig = result[ARG_INPUT_FILENAME_CONTIG].as<string>();
			_trimmedOutliers = result[ARG_TRIMMED_OUTLIERS].as<float>();
			_nbCores = result[ARG_NB_CORES].as<int>();
			_kmerDensity = result[ARG_MINIMIZER_DENSITY].as<float>(); //getInput()->getDouble(STR_DENSITY);
			_kmerSize = result[ARG_KMINMER_LENGTH].as<int>(); //getInput()->getInt(STR_MINIM_SIZE);
			
			cout << "Kmer size: " << _kmerSize << endl;
			cout << "Kmer density: " << _kmerDensity << endl;
			cout << "Trimmed val: " << _trimmedOutliers << endl;
		}
		catch (const std::exception& e){
			std::cout << options.help() << std::endl;
			std::cerr << e.what() << std::endl;
			std::exit(EXIT_FAILURE);
		}

		//_outputFilename_kmerCouts = _inputDir + "/" + "contigCoverages.tsv";
	}

    void execute (){
    
		_nbContigs = 0;
		//_kmerModel = new KmerModel(31);
		//_minimizerParser = new MinimizerParser(31, 0.1);

		ReadParser parserDatasets(_inputFilename_shortreads, false, false);
		_nbDatasets = parserDatasets._nbDatasets;
		_nbDatasets /= 2; //PAIRED
		cout << "Nb datasets: " << _nbDatasets << endl;

		for(u_int64_t i=0; i<_nbDatasets; i++) _countsInit.push_back(0);

		//_tmpDir = _outputFilename + ".tmp";
		//fs::path path(_tmpDir);
	    //if(!fs::exists (path)){
        //    fs::create_directory(path);
        //} 

		extractContigKmers();
		countShortreadsKmers();
		dumpContigCoverages();

		//countKminmers();
		//computeContigCoverage();
		//dumpKminmerCount();


	}

	void extractContigKmers (){

		ReadParserParallel readParser(_inputFilename_contig, true, false, _nbCores);
		readParser.parse(ContigFunctor(*this));

		for(size_t i=0; i<_nbContigs; i++){
			_contigCoverages_mean.push_back(_countsInit);
			_contigCoverages_var.push_back(_countsInit);
		}

		cout << "Nb kmers: " << _kmerCounts.size() << endl;
		//ReadParser parser(_inputFilename_contig, true, false);
		//auto fp = std::bind(&KmerCounter::extractContigKmers_read, this, std::placeholders::_1);
		//parser.parse(fp);

		//for(const auto& it : _kmerCounts){
		//	_sortedKmers.push_back(it.first);
		//}

		//std::sort(_sortedKmers.begin(), _sortedKmers.end());
	}

	/*
	void extractContigKmers_read(const Read& read){

		u_int64_t readIndex = read._index;

		_contigCoverages_mean.push_back(_countsInit);
		_contigCoverages_var.push_back(_countsInit);
		//string rleSequence;
		//vector<u_int64_t> rlePositions;
		//Encoder::encode_rle(read->seq.s, strlen(read->seq.s), rleSequence, rlePositions);

		vector<u_int64_t> minimizers;
		vector<u_int64_t> minimizers_pos;
		_minimizerParser->parse(read._seq, minimizers, minimizers_pos);

		cout << "----" << endl;
		cout << readIndex << endl;
		cout << _countsInit.size() << endl;
		cout << _kmerCounts.size() << endl;

		for(u_int64_t minimizer : minimizers){
			if(_kmerCounts.find(minimizer) == _kmerCounts.end()) _kmerCounts[minimizer] = 0;
		}


	}
	*/

	void countShortreadsKmers (){

		_currentDatasetIndex = 0;
		ReadParser parser(_inputFilename_shortreads, false, false);
		u_int64_t datasetIndex = 0;

		for(const string& filename : parser._filenames){
			
			cout << "Countign kmers: " << filename << endl;

			if(datasetIndex % 2 == 0){
				cout << "\tNew dataset" << endl;
				for(auto& it : _kmerCounts){
					_kmerCounts[it.first] = 0;
				}

			}
			else{
				cout << "\tProcessing pair" << endl;
			}



			ReadParserParallel readParser(filename, true, false, _nbCores);
			readParser.parse(KmerCounterFunctor(*this));
			//ReadParser parserDataset(filename, true, false);
			//auto fp = std::bind(&KmerCounter::countShortreadsKmers_read, this, std::placeholders::_1);
			//parserDataset.parse(fp);

			if(datasetIndex % 2 == 1){
				computeContigCoverage();
				_currentDatasetIndex += 1;
			}
			/*
			string kmerCountFilename = _tmpDir + "/" + "kmerCounts_" + to_string(datasetIndex) + ".gz";
			gzFile kmerCountFile = gzopen(kmerCountFilename.c_str(), "wb");

			for(u_int64_t kmer : _sortedKmers){
				u_int32_t count = _kmerCounts[kmer];
				gzwrite(kmerCountFile, (const char*)&count, sizeof(count));
			}

			gzclose(kmerCountFile);
			*/
			datasetIndex += 1;
		}
		


	}

	/*
	void countShortreadsKmers_read(const Read& read){

		u_int64_t readIndex = read._index;
		if(readIndex % 100000 == 0) cout << _currentDatasetIndex << " " << readIndex << endl;
		
		//string rleSequence;
		//vector<u_int64_t> rlePositions;
		//Encoder::encode_rle(read->seq.s, strlen(read->seq.s), rleSequence, rlePositions);

		vector<u_int64_t> minimizers;
		vector<u_int64_t> minimizers_pos;
		_minimizerParser->parse(read._seq, minimizers, minimizers_pos);


		for(u_int64_t minimizer : minimizers){
			if(_kmerCounts.find(minimizer) != _kmerCounts.end()) _kmerCounts[minimizer] += 1;
		}

		
	}
	*/


	void computeContigCoverage (){

		cout << "\tComputing contig coverage" << endl;
		
		ReadParserParallel readParser(_inputFilename_contig, true, false, _nbCores);
		readParser.parse(ContigCoverageFunctor(*this));

		//ReadParser parser(_inputFilename_contig, true, false);
		//auto fp = std::bind(&KmerCounter::computeContigCoverage_read, this, std::placeholders::_1);
		//parser.parse(fp);

	}

	/*
	void computeContigCoverage_read(const Read& read){

		u_int64_t readIndex = read._index;
		//string rleSequence;
		//vector<u_int64_t> rlePositions;
		//Encoder::encode_rle(read->seq.s, strlen(read->seq.s), rleSequence, rlePositions);

		vector<u_int32_t> count_values;

		vector<u_int64_t> minimizers;
		vector<u_int64_t> minimizers_pos;
		_minimizerParser->parse(read._seq, minimizers, minimizers_pos);


		for(u_int64_t minimizer : minimizers){
			count_values.push_back(_kmerCounts[minimizer]);
		}

		std::sort(count_values.begin(), count_values.end());

		double sum = 0;
		double n = 0;

		u_int64_t nbTrimmedValues = count_values.size() * _trimmedOutliers;

		for(size_t i=nbTrimmedValues; i<count_values.size()-nbTrimmedValues; i++){
			sum += count_values[i];
			n += 1;
		}

		if(n <= 1) return;

		float mean = sum / n;

		double var = 0;
		for(size_t i=nbTrimmedValues; i<count_values.size()-nbTrimmedValues; i++){
			double count = count_values[i];
			var += ((count - mean) * (count - mean));
		}

		var /= (n-1);

		cout << readIndex << " " << mean << " " << var << endl; 

		_contigCoverages_mean[readIndex][_currentDatasetIndex] = mean;
		_contigCoverages_var[readIndex][_currentDatasetIndex] = var;


	}
	*/

	void dumpContigCoverages(){

		ofstream outputFile(_outputFilename);

		outputFile << "contigName\tcontigLen\ttotalAvgDepth" << endl;

		for(size_t i=0; i<_contigCoverages_mean.size(); i++){
			//u_int32_t contigIndex = i;
			const vector<float>& means = _contigCoverages_mean[i];
			const vector<float>& vars = _contigCoverages_var[i];

			outputFile << "ctg" << i;
			for(size_t i=0; i<means.size(); i++){
				outputFile << "\t" << means[i] << "\t" << vars[i];
			}
			outputFile << endl;
		}

		outputFile.close();
	}


	class ContigFunctor {

		public:

		KmerCounter& _kmerCounter;
		MinimizerParser128* _minimizerParser;
		//KmerModel128* _kmerModel;
		u_int64_t _currentDatasetIndex;

		ContigFunctor(KmerCounter& kmerCounter) : _kmerCounter(kmerCounter){
			
		}

		ContigFunctor(const ContigFunctor& copy) : _kmerCounter(copy._kmerCounter){
			//_kmerModel = new KmerModel128(_kmerCounter._kmerSize);
			_minimizerParser = new MinimizerParser128(_kmerCounter._kmerSize, _kmerCounter._kmerDensity);
			_currentDatasetIndex = _kmerCounter._currentDatasetIndex;
		}

		~ContigFunctor(){
		}

		void operator () (const Read& read) {

			//#pragma omp critical
			//{
			//	cout << read._index << " " << read._seq.size() << endl;
			//}

			#pragma omp atomic
			_kmerCounter._nbContigs += 1;

			u_int64_t readIndex = read._index;
			if(readIndex % 1000 == 0) cout << readIndex << endl;
			
			//string rleSequence;
			//vector<u_int64_t> rlePositions;
			//Encoder::encode_rle(read->seq.s, strlen(read->seq.s), rleSequence, rlePositions);

			vector<KmerType> minimizers;
			vector<u_int64_t> minimizers_pos;
			_minimizerParser->parse(read._seq, minimizers, minimizers_pos);

			//cout << "----" << endl;
			//cout << readIndex << endl;
			//cout << _countsInit.size() << endl;
			//cout << _kmerCounts.size() << endl;

			//#pragma omp critical
			//{
			for(KmerType minimizer : minimizers){

					
				/*
				_kmerCounter._kmerCounts.try_emplace_l(minimizer, 
				[](KmerCountMap::value_type& v) { // key exist
					//v.second += 1;
				}, 0);
				*/
				/*
				_kmerCounter._kmerCounts.lazy_emplace_l(minimizer, 
				[](KmerCountMap::value_type& v) { // key exist
					//v.second += 1;
				},           
				[&minimizer](const KmerCountMap::constructor& ctor) { // key inserted
					ctor(minimizer, 0);
				});
				*/

				if(_kmerCounter._kmerCounts.find(minimizer) == _kmerCounter._kmerCounts.end()) _kmerCounter._kmerCounts[minimizer] = 0;
			}
			//}


		
		}

	};

	class KmerCounterFunctor {

		public:

		KmerCounter& _kmerCounter;
		MinimizerParser128* _minimizerParser;
		//KmerModel128* _kmerModel;
		u_int64_t _currentDatasetIndex;

		KmerCounterFunctor(KmerCounter& kmerCounter) : _kmerCounter(kmerCounter){
			
		}

		KmerCounterFunctor(const KmerCounterFunctor& copy) : _kmerCounter(copy._kmerCounter){
			//_kmerModel = new KmerModel128(_kmerCounter._kmerSize);
			_minimizerParser = new MinimizerParser128(_kmerCounter._kmerSize, _kmerCounter._kmerDensity);
			_currentDatasetIndex = _kmerCounter._currentDatasetIndex;
		}

		~KmerCounterFunctor(){
		}

		void operator () (const Read& read) {


			u_int64_t readIndex = read._index;
			if(readIndex % 100000 == 0) cout << _currentDatasetIndex << " " << readIndex << endl;
			
			//string rleSequence;
			//vector<u_int64_t> rlePositions;
			//Encoder::encode_rle(read->seq.s, strlen(read->seq.s), rleSequence, rlePositions);

			vector<KmerType> minimizers;
			vector<u_int64_t> minimizers_pos;
			_minimizerParser->parse(read._seq, minimizers, minimizers_pos);


			//#pragma omp critical
			//{
			for(KmerType minimizer : minimizers){

				auto set_value = [](KmerCountMap::value_type& v) { v.second += 1; };
    			_kmerCounter._kmerCounts.modify_if(minimizer, set_value);

				/*
				_kmerCounter._kmerCounts->_dbg_nodes.lazy_emplace_l(minimizer, 
				[](MdbgNodeMap::value_type& v) { // key exist
					v.second += 1;
				},           
				[&minimizer](const MdbgNodeMap::constructor& ctor) { // key inserted
					//ctor(minimizer, 1);
				});
				*/

				//if(_kmerCounter._kmerCounts.find(minimizer) != _kmerCounter._kmerCounts.end()) _kmerCounter._kmerCounts[minimizer] += 1;
			}
			//}



		
		}

	};


	class ContigCoverageFunctor {

		public:

		KmerCounter& _kmerCounter;
		MinimizerParser128* _minimizerParser;
		//KmerModel128* _kmerModel;
		u_int64_t _currentDatasetIndex;

		ContigCoverageFunctor(KmerCounter& kmerCounter) : _kmerCounter(kmerCounter){
			
		}

		ContigCoverageFunctor(const ContigCoverageFunctor& copy) : _kmerCounter(copy._kmerCounter){
			//_kmerModel = new KmerModel128(_kmerCounter._kmerSize);
			_minimizerParser = new MinimizerParser128(_kmerCounter._kmerSize, _kmerCounter._kmerDensity);
			_currentDatasetIndex = _kmerCounter._currentDatasetIndex;
		}

		~ContigCoverageFunctor(){
		}

		void operator () (const Read& read) {

			u_int64_t readIndex = read._index;
			if(readIndex % 1000 == 0) cout << readIndex << endl;
			
			//string rleSequence;
			//vector<u_int64_t> rlePositions;
			//Encoder::encode_rle(read->seq.s, strlen(read->seq.s), rleSequence, rlePositions);

			vector<KmerType> minimizers;
			vector<u_int64_t> minimizers_pos;
			_minimizerParser->parse(read._seq, minimizers, minimizers_pos);

			vector<u_int32_t> count_values;
			for(KmerType minimizer : minimizers){
				count_values.push_back(_kmerCounter._kmerCounts[minimizer]);
			}

			std::sort(count_values.begin(), count_values.end());

			double sum = 0;
			double n = 0;

			u_int64_t nbTrimmedValues = count_values.size() * _kmerCounter._trimmedOutliers;

			for(size_t i=nbTrimmedValues; i<count_values.size()-nbTrimmedValues; i++){
				sum += count_values[i];
				n += 1;
			}

			if(n <= 1) return;

			float mean = sum / n;

			double var = 0;
			for(size_t i=nbTrimmedValues; i<count_values.size()-nbTrimmedValues; i++){
				double count = count_values[i];
				var += ((count - mean) * (count - mean));
			}

			var /= (n-1);

			//

			_kmerCounter._contigCoverages_mean[readIndex][_currentDatasetIndex] = mean;
			_kmerCounter._contigCoverages_var[readIndex][_currentDatasetIndex] = var;

			//#pragma omp critical
			//{
			//	cout << readIndex << " " << mean << " " << var << endl; 
			//}


		
		}

		

	};

};	


#endif 



