

#ifndef MDBG_METAG_PURGEGRAPH
#define MDBG_METAG_PURGEGRAPH

#include "../Commons.hpp"
#include "../utils/edlib.h"
#include "../utils/spoa/include/spoa/spoa.hpp"
#include "../utils/DnaBitset.hpp"
#include "../utils/abPOA2/include/abpoa.h"


class PurgeGraph : public Tool{
    
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
	bool _circularize;
	u_int64_t _minContigLength;
	
	string _outputFilename_contigs;
	string _outputFilename_mapping;


	PurgeGraph(): Tool (){

	}


	void parseArgs(int argc, char* argv[]){


		string ARG_USE_QUAL = "qual";
		string ARG_CIRCULARIZE = "circ";
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
		(ARG_CIRCULARIZE, "", cxxopts::value<bool>()->default_value("false"))
		(ARG_NB_CORES, "", cxxopts::value<int>()->default_value(NB_CORES_DEFAULT));


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
			_circularize = result[ARG_CIRCULARIZE].as<bool>();
			_windowLength = 500;
			_maxWindowCopies = 21; //21;
			_minContigLength = 1000000;
			
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

	}


	gzFile _outputContigFile;
	string _alignFilename;

    void execute (){

		_alignFilename = "/mnt/gpfs/gaetan/tmp/bin.1081_2.paf";
		//mapReads();
		indexMappingOnContigEnds();

		//gzclose(_outputContigFile);
		//fs::remove_all(_tmpDir);

	}

	void mapReads(){
		
		string readFilenames = "";
		ReadParser readParser(_inputFilename_reads, false, false);

		for(const string& filename : readParser._filenames){
			readFilenames += filename + " ";
		}

		string command = "minimap2 -H -I 2G -t " + to_string(_nbCores) + " -x map-hifi " + _inputFilename_contigs + " " + readFilenames;
		command += " | " + _mapperOutputExeFilename + " " + _inputFilename_contigs + " " + _inputFilename_reads + " " + _outputFilename_mapping;
		Utils::executeCommand(command, _outputDir);

		//minimap2 -x map-hifi ~/workspace/run/overlap_test_201/contigs_47.fasta.gz ~/workspace/data/overlap_test/genome_201_50x/simulatedReads_0.fastq.gz | ./bin/mapper ~/workspace/run/overlap_test_201/contigs_47.fasta.gz ~/workspace/data/overlap_test/genome_201_50x/input.txt ~/workspace/run/overlap_test_201/align.bin
	}

	void indexMappingOnContigEnds(){

		cout << "\tIndexing read alignments" << endl;

    	long maxHang = 300;

		vector<string>* fields = new vector<string>();
		vector<string>* fields_optional = new vector<string>();

		ifstream infile(_alignFilename);
		string line;

		while (getline(infile, line)) {
			//cout << lineInput << endl;
			//getchar();


			GfaParser::tokenize(line, fields, '\t');

			//cout << line << endl;

			const string& readName = (*fields)[0];
			const string& contigName = (*fields)[5];

			u_int32_t readStart = stoull((*fields)[2]);
			u_int32_t readEnd = stoull((*fields)[3]);
			u_int32_t contigLength = stoull((*fields)[6]);
			u_int32_t contigStart = stoull((*fields)[7]);
			u_int32_t contigEnd = stoull((*fields)[8]);

			u_int64_t nbMatches = stoull((*fields)[9]);
			u_int64_t alignLength = stoull((*fields)[10]);
			u_int64_t queryLength = stoull((*fields)[1]);

			bool strand = (*fields)[4] == "-";

			long hangLeft = contigStart;
			long hangRight = contigLength - contigEnd;

			if (hangLeft < maxHang){
				cout << line << endl;
				//_contigMapLeft[contigIndex].push_back(readIndex);
				//cout << readIndex << endl;
			}
			if (hangRight < maxHang){
				cout << line << endl;
				//_contigMapRight[contigIndex].push_back(readIndex);
				//cout << readIndex << endl;
			}

		}

		infile.close();

	}



};	

#endif 


