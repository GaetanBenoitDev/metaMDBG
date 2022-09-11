

#ifndef MDBG_METAG_PURGEDUPS
#define MDBG_METAG_PURGEDUPS

#include "../Commons.hpp"
#include "../utils/edlib.h"
#include "../utils/DnaBitset.hpp"



class PurgeDups : public Tool{
    
public:

	//string _inputFilename_reads;
	string _inputFilename_contigs;
	int _nbCores;
	string _mapperOutputExeFilename;
	string _outputDir;
	string _tmpDir;
	
	string _outputFilename_contigs;
	string _outputFilename_mapping;
	u_int64_t _maxBases;

	u_int64_t _minDuplicationLength_ends;
	float _minDuplicationIdentity_ends;
	u_int64_t _minDuplicationLength_internal;
	float _minDuplicationIdentity_internal;
	//u_int64_t _maxMemory;

	struct Alignment{
		u_int32_t _contigIndex;
		//u_int64_t _readIndex;
		bool _strand;
		u_int32_t _readStart;
		u_int32_t _readEnd;
		u_int32_t _contigStart;
		u_int32_t _contigEnd;
		//float _score;
		//u_int64_t _length;
		
		u_int64_t length(){
			return std::max((u_int64_t)(_readEnd - _readStart), (u_int64_t)(_contigEnd - _contigStart));
		}
		
	};

	PurgeDups(): Tool (){

	}


	void parseArgs(int argc, char* argv[]){

		string filenameExe = argv[0];
		//cout << filenameExe << endl;

		fs::path pa(filenameExe);
		_mapperOutputExeFilename = pa.parent_path().string() + "/mapper";
		//cout << _mapperOutputExeFilename << endl;
		//exit(1);

		cxxopts::Options options("ToBasespace", "");
		options.add_options()
		("contigs", "", cxxopts::value<string>())
		//("reads", "", cxxopts::value<string>())
		("tmpDir", "", cxxopts::value<string>())
		(ARG_NB_CORES, "", cxxopts::value<int>()->default_value(NB_CORES_DEFAULT));


		options.parse_positional({"contigs", "tmpDir"});
		options.positional_help("contigs tmpDir");


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

			//_inputFilename_reads = result["reads"].as<string>();
			_inputFilename_contigs = result["contigs"].as<string>();
			_outputDir = result["tmpDir"].as<string>();;
			_nbCores = result[ARG_NB_CORES].as<int>();
		}
		catch (const std::exception& e){
			std::cout << options.help() << std::endl;
			std::cerr << e.what() << std::endl;
			std::exit(EXIT_FAILURE);
		}

		cout << "Contigs: " << _inputFilename_contigs << endl;

		fs::path p(_inputFilename_contigs);
		while(p.has_extension()){
			p.replace_extension("");
		}

		_tmpDir = _outputDir + "/__tmp/";
		if(!fs::exists(_tmpDir)){
			fs::create_directories(_tmpDir);
		}

		_outputFilename_contigs = p.string() + "_derep.fasta.gz";
		_maxBases = 100000000ull;
		_minDuplicationLength_ends = 0.9;
		_minDuplicationIdentity_ends = 1000;
		_minDuplicationLength_internal = 10000;
		_minDuplicationIdentity_internal = 0.95;
		//_outputFilename_contigs = p.string() + "_corrected.fasta.gz";
		_outputFilename_mapping = _tmpDir + "/_tmp_mapping_derep__.paf";
		//_maxMemory = 4000000000ull;
		cout << _outputFilename_contigs << endl;
	}

	gzFile _outputContigFile;
	unordered_map<u_int32_t, string> _contigSequences;
	unordered_map<u_int32_t, vector<Alignment>> _alignments;

    void execute (){

		_outputContigFile = gzopen(_outputFilename_contigs.c_str(), "wb");
		
		mapReads();
		processContigs();
		
		gzclose(_outputContigFile);
		//fs::remove_all(_tmpDir);

	}

	void mapReads(){
		
		string inputContigsFilename = _tmpDir + "/input_contigs.txt";
		ofstream input(inputContigsFilename);
		input << _inputFilename_contigs << endl;
		input.close();

		string command = "minimap2 -m 1000 -H -X -I 2G -t " + to_string(_nbCores) + " -x map-hifi " + _inputFilename_contigs + " " + _inputFilename_contigs;
		command += " | " + _mapperOutputExeFilename + " " + _inputFilename_contigs + " " + inputContigsFilename + " " + _outputFilename_mapping;
		Utils::executeCommand(command, _outputDir);
	}

	u_int64_t _currentLoadedBases;

	void processContigs(){

		auto fp = std::bind(&PurgeDups::processContigs_read, this, std::placeholders::_1);
		ReadParser readParser(_inputFilename_contigs, true, false);
		readParser.parse(fp);

		if(_currentLoadedBases > 0){
			processPass();
		}
	}

	
	void processContigs_read(const Read& read){

		u_int32_t contigIndex = read._index;

		_contigSequences[contigIndex] = read._seq;
		_currentLoadedBases += read._seq.size();

		if(_currentLoadedBases > _maxBases){
			processPass();
		}

	}

	void processPass(){
		cout << "Processing " <<  _contigSequences.size() << " contigs" << endl;
		indexAlignments();
		detectDuplication();
		clearPass();
		cout << "Pass done" << endl;
	}

	void clearPass(){
		_alignments.clear();
		_contigSequences.clear();
		_currentLoadedBases = 0;
	}

	void indexAlignments(){

		cout << "\tIndexing alignments" << endl;

        ifstream infile(_outputFilename_mapping);

        std::string line;

		while(true){

			u_int32_t contigIndex;
			u_int32_t contigLength;
			u_int64_t readIndex;
			u_int32_t readStart;
			u_int32_t readEnd;
			u_int32_t contigStart;
			u_int32_t contigEnd;
			float score;
			bool strand;

			infile.read((char*)&contigIndex, sizeof(contigIndex));
			if(infile.eof()) break;
			infile.read((char*)&contigLength, sizeof(contigLength));
			infile.read((char*)&contigStart, sizeof(contigStart));
			infile.read((char*)&contigEnd, sizeof(contigEnd));
			infile.read((char*)&readIndex, sizeof(readIndex));
			infile.read((char*)&readStart, sizeof(readStart));
			infile.read((char*)&readEnd, sizeof(readEnd));
			infile.read((char*)&strand, sizeof(strand));
			infile.read((char*)&score, sizeof(score));


			if(_contigSequences.find(contigIndex) == _contigSequences.end() && _contigSequences.find(readIndex) == _contigSequences.end()) continue;

			u_int32_t length = std::max((u_int64_t)(readEnd - readStart), (u_int64_t)(contigEnd - contigStart));
			if(length < 1000) continue;

			Alignment align = {contigIndex, strand, readStart, readEnd, contigStart, contigEnd}; //, score
			_alignments[readIndex].push_back(align);
			/*
			u_int32_t length = std::max((u_int64_t)(readEnd - readStart), (u_int64_t)(contigEnd - contigStart));

			if(_alignments.find(readIndex) == _alignments.end()){
				_alignments[readIndex] = align;
			}
			else{

				if(length > align.length()){
					_alignments[readIndex] = align;
				}
				
			}
			*/
        }
	}

	void detectDuplication(){
		cout << "\tDetecting duplication" << endl;

		//auto fp = std::bind(&ContigPolisher::collectWindowCopies_read, this, std::placeholders::_1);
		//ReadParser readParser(_inputFilename_reads, false, false);
		//readParser.parse(fp);


		//const string& partitionFilename = _tmpDir + "/part_" + to_string(partition) + ".gz";
		ReadParserParallel readParser(_inputFilename_contigs, true, false, _nbCores);
		readParser.parse(ContigAlignerFunctor(*this));
	}



	class ContigAlignerFunctor {

		public:

		PurgeDups& _purgeDups;
		//unordered_map<string, u_int32_t>& _contigName_to_contigIndex;
		//unordered_map<string, u_int64_t>& _readName_to_readIndex;
		//unordered_map<u_int64_t, vector<Alignment>>& _alignments;
		unordered_map<u_int32_t, vector<Alignment>>& _alignments;
		unordered_map<u_int32_t, string>& _contigSequences;
		//unordered_map<u_int32_t, vector<vector<Window>>>& _contigWindowSequences;
		//size_t _windowLength;


		ContigAlignerFunctor(PurgeDups& purgeDups) : _purgeDups(purgeDups), _alignments(purgeDups._alignments), _contigSequences(purgeDups._contigSequences){
		}

		ContigAlignerFunctor(const ContigAlignerFunctor& copy) : _purgeDups(copy._purgeDups), _alignments(copy._alignments), _contigSequences(copy._contigSequences){
			
		}

		~ContigAlignerFunctor(){
		}

		void operator () (const Read& read) {

			
			//u_int64_t readIndex = stoull(read._header);

			//#pragma omp critical
			//{
				//cout << read._index << endl;
			//}

			if(_alignments.find(read._index) == _alignments.end()) return;

			//cout << "\t" << _alignments[read._index].size() << endl;
			for(const Alignment& al : _alignments[read._index]){

				if(read._index == al._contigIndex){

					//cout << "Self AL not normal" << endl;
					continue;
				}
				
				if(_contigSequences.find(al._contigIndex) == _contigSequences.end()) continue;

				//cout << al._readStart << " " << al._readEnd << endl;
				processAlignment(read, al);
				//cout << read._index << " " << al._contigIndex << endl;
			}
			//cout << "\tdone" << endl;
			/*
			//if(_contigPolisher._currentPartition == 0) cout << readIndex << " " << (_alignments.find(readIndex) != _alignments.end()) << endl;
			
			//if(readIndex % 100000 == 0) cout << "\t" << readIndex << endl;

			if(_alignments.find(readIndex) == _alignments.end()) return;

			//const vector<Alignment>& als = _alignments[readIndex];
			const Alignment& al = _alignments[readIndex];
			//for(const Alignment& al : _alignments[readIndex]){
			u_int64_t contigIndex = al._contigIndex;

			if(_contigSequences.find(contigIndex) == _contigSequences.end()) return;

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

			//}
			*/
		}

		void processAlignment(const Read& read, const Alignment& al){

			u_int32_t readIndex = read._index;
			u_int32_t contigIndex = al._contigIndex;

			string readSeq = read._seq;
			string qualSeq = read._qual;
			string readSequence = readSeq.substr(al._readStart, al._readEnd-al._readStart);
			string contigSequence = _contigSequences[contigIndex].substr(al._contigStart, al._contigEnd-al._contigStart);



			if(al._strand){
				Utils::toReverseComplement(readSequence);
				//Utils::toReverseComplement(readSeq);
			}

			//cout << contigSequence << endl;
			//cout << readSequence << endl;
			//exit(1);

			static EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0);


			//cout << "Comparing: " << readSequence.size() << " " << contigSequence.size() << endl;

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
			free(cigar);

			//cout << cigar << endl;
			//exit(1);

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

			/*
			//if(breaking_points_.size() > 0) breaking_points_.emplace_back(last_match);
				
			for (uint32_t j = 0; j < breaking_points_.size(); j += 2) {
				//if (breaking_points_[j + 1].second - breaking_points_[j].second < 0.02 * _windowLength) {
				//	continue;
				//}				

				//uint64_t window_id = id_to_first_window_id[overlaps[i]->t_id()] +


				//const char* data = overlaps[i]->strand() ?
				//	&(sequence->reverse_complement()[breaking_points[j].second]) :
				//	&(sequence->data()[breaking_points[j].second]);
				const char* data = &readSequence[breaking_points_[j].second];
				uint32_t data_length = breaking_points_[j + 1].second - breaking_points_[j].second;

				string sequence = string(data, data_length);


				data = &readSequence[breaking_points_[j].second];
				data_length = breaking_points_[j + 1].second - breaking_points_[j].second;
				sequence = string(data, data_length);


                u_int32_t posStart = breaking_points_[j].first - window_start;
                u_int32_t posEnd =  breaking_points_[j + 1].first - window_start - 1;


				//indexWindow(al, window_id, posStart, posEnd, sequence, quality);

			}
			*/
			
		}

	};

};	

#endif 


