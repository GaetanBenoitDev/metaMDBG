

#ifndef MDBG_METAG_KMERCOUNTER
#define MDBG_METAG_KMERCOUNTER

#include "../Commons.hpp"


class KmerCounter : public Tool{
    
public:

	string _inputFilename_contig;
	string _inputFilename_shortreads;
	string _outputFilename;
    size_t _kmerSize;

	//string _tmpDir;
	//string _outputFilename_kmerCouts;
	vector<float> _countsInit;
	size_t _nbDatasets;

	KmerModel* _kmerModel;
	unordered_map<u_int64_t, u_int32_t> _kmerCounts;
	vector<vector<float>> _contigCoverages_mean;
	vector<vector<float>> _contigCoverages_var;
	u_int32_t _currentDatasetIndex;

	MinimizerParser* _minimizerParser;

	KmerCounter(): Tool (){

	}

	void parseArgs(int argc, char* argv[]){


		cxxopts::Options options("ToMinspace", "");
		options.add_options()
		//(ARG_INPUT_FILENAME, "", cxxopts::value<string>())
		(ARG_OUTPUT_DIR, "", cxxopts::value<string>())
		(ARG_INPUT_FILENAME, "", cxxopts::value<string>()->default_value(""))
		(ARG_INPUT_FILENAME_CONTIG, "", cxxopts::value<string>()->default_value(""));


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
			
			
		}
		catch (const std::exception& e){
			std::cout << options.help() << std::endl;
			std::cerr << e.what() << std::endl;
			std::exit(EXIT_FAILURE);
		}

		//_outputFilename_kmerCouts = _inputDir + "/" + "contigCoverages.tsv";
	}

    void execute (){
    
		_kmerModel = new KmerModel(31);
		_minimizerParser = new MinimizerParser(31, 0.1);

		ReadParser parserDatasets(_inputFilename_shortreads, false);
		_nbDatasets = parserDatasets._nbDatasets;
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

		"reprise: enelever bord des contigs
		retirer outliers"
		//countKminmers();
		//computeContigCoverage();
		//dumpKminmerCount();


	}

	void extractContigKmers (){

		ReadParser parser(_inputFilename_contig, true);
		auto fp = std::bind(&KmerCounter::extractContigKmers_read, this, std::placeholders::_1, std::placeholders::_2);
		parser.parse(fp);

		//for(const auto& it : _kmerCounts){
		//	_sortedKmers.push_back(it.first);
		//}

		//std::sort(_sortedKmers.begin(), _sortedKmers.end());
	}


	void extractContigKmers_read(kseq_t* read, u_int64_t readIndex){

		_contigCoverages_mean.push_back(_countsInit);
		_contigCoverages_var.push_back(_countsInit);
		//string rleSequence;
		//vector<u_int64_t> rlePositions;
		//Encoder::encode_rle(read->seq.s, strlen(read->seq.s), rleSequence, rlePositions);

		vector<u_int64_t> minimizers;
		vector<u_int64_t> minimizers_pos;
		_minimizerParser->parse(string(read->seq.s), minimizers, minimizers_pos);

		cout << "----" << endl;
		cout << readIndex << endl;
		cout << _countsInit.size() << endl;
		cout << _kmerCounts.size() << endl;

		for(u_int64_t minimizer : minimizers){
			if(_kmerCounts.find(minimizer) == _kmerCounts.end()) _kmerCounts[minimizer] = 0;
		}

		/*
		vector<u_int64_t> kmers;
		_kmerModel->iterate(read->seq.s, strlen(read->seq.s), kmers);

		for(u_int64_t kmer : kmers){
			if(_kmerCounts.find(kmer) == _kmerCounts.end()) _kmerCounts[kmer] = _countsInit;
		}
		cout << _kmerCounts.size() << endl;
		*/
	}

	void countShortreadsKmers (){

		_currentDatasetIndex = 0;
		ReadParser parser(_inputFilename_shortreads, false);
		u_int64_t datasetIndex = 0;

		for(const string& filename : parser._filenames){
			
			cout << "Countign kmers: " << filename << endl;

			for(auto& it : _kmerCounts){
				_kmerCounts[it.first] = 0;
			}

			ReadParser parserDataset(filename, true);
			auto fp = std::bind(&KmerCounter::countShortreadsKmers_read, this, std::placeholders::_1, std::placeholders::_2);
			parserDataset.parse(fp);

			computeContigCoverage();
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
			_currentDatasetIndex += 1;
		}
		


	}


	void countShortreadsKmers_read(kseq_t* read, u_int64_t readIndex){

		if(readIndex % 100000 == 0) cout << _currentDatasetIndex << " " << readIndex << endl;
		
		//string rleSequence;
		//vector<u_int64_t> rlePositions;
		//Encoder::encode_rle(read->seq.s, strlen(read->seq.s), rleSequence, rlePositions);

		vector<u_int64_t> minimizers;
		vector<u_int64_t> minimizers_pos;
		_minimizerParser->parse(string(read->seq.s), minimizers, minimizers_pos);


		for(u_int64_t minimizer : minimizers){
			if(_kmerCounts.find(minimizer) != _kmerCounts.end()) _kmerCounts[minimizer] += 1;
		}

		
	}


	void computeContigCoverage (){

		ReadParser parser(_inputFilename_contig, true);
		auto fp = std::bind(&KmerCounter::computeContigCoverage_read, this, std::placeholders::_1, std::placeholders::_2);
		parser.parse(fp);

	}


	void computeContigCoverage_read(kseq_t* read, u_int64_t readIndex){

		//string rleSequence;
		//vector<u_int64_t> rlePositions;
		//Encoder::encode_rle(read->seq.s, strlen(read->seq.s), rleSequence, rlePositions);

		vector<u_int64_t> minimizers;
		vector<u_int64_t> minimizers_pos;
		_minimizerParser->parse(string(read->seq.s), minimizers, minimizers_pos);

		double sum = 0;
		//u_int64_t n = 0;

		for(u_int64_t minimizer : minimizers){
			sum += _kmerCounts[minimizer];
			//n += 1;
			//if(_kmerCounts.find(minimizer) == _kmerCounts.end()) _kmerCounts[minimizer] = 0;
		}

		float mean = sum / minimizers.size();

		double var = 0;
		for(u_int64_t minimizer : minimizers){
			double count = _kmerCounts[minimizer];
			var += ((count - mean) * (count - mean));
		}

		var /= (minimizers.size()-1);

		cout << readIndex << " " << mean << " " << var << endl; 

		_contigCoverages_mean[readIndex][_currentDatasetIndex] = mean;
		_contigCoverages_var[readIndex][_currentDatasetIndex] = var;
		/*
		vector<u_int64_t> kmers;
		_kmerModel->iterate(read->seq.s, strlen(read->seq.s), kmers);

		for(u_int64_t kmer : kmers){
			if(_kmerCounts.find(kmer) == _kmerCounts.end()) _kmerCounts[kmer] = _countsInit;
		}
		cout << _kmerCounts.size() << endl;
		*/
	}

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

};	


#endif 



