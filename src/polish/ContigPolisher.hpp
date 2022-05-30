
/*
- attention ne pas outupt empty contig
- ToBasespace: générer un assemblage grossier
- tester les scores de qualité ?
- quand un read a été process, on peut erase son entrée dans _alignments ?
- window sequence a selectionner en priorité: 
	- la distance est une mauvaise metrique car on ne sait pas si la sequence de reference est erroné ou non
	- un mapping tres long sur un contig a plus de valeur qu'un mapping court (on est plus sûr que ce read appartient au contig)+
- si "no sequencer fo window": utiliser la fenetre original du contig ?
- minimap2 output: écrire un petit programme pour compresser les résultats d'alignements on the fly
- version finale: remove le minimap2 align filename
- voir comment racon handle multimap
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
	
	string _outputFilename_contigs;
	string _outputFilename_mapping;

	abpoa_para_t *abpt;
	
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


		cxxopts::Options options("ToBasespace", "");
		options.add_options()
		("contigs", "", cxxopts::value<string>())
		("reads", "", cxxopts::value<string>())
		(ARG_NB_CORES, "", cxxopts::value<int>()->default_value("4"));

		options.parse_positional({"contigs", "reads"});
		options.positional_help("contigs reads");


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
			_nbCores = result[ARG_NB_CORES].as<int>();
			_windowLength = 500;
			
		}
		catch (const std::exception& e){
			std::cout << options.help() << std::endl;
			std::cerr << e.what() << std::endl;
			std::exit(EXIT_FAILURE);
		}

		cout << "Contigs: " << _inputFilename_contigs << endl;
		cout << "Reads: " << _inputFilename_reads << endl;

		fs::path p(_inputFilename_contigs);
		while(p.has_extension()){
			p.replace_extension("");
		}

		_outputFilename_contigs = p.string() + "_corrected.fasta.gz";
		_outputFilename_mapping = p.string() + "_tmp_mapping__.paf";


	}



    void execute (){


		//mapReads();
		indexContigName();
		indexReadName();
		parseAlignments();

		_contigName_to_contigIndex.clear();
		_readName_to_readIndex.clear();

		loadContigs();
		collectWindowSequences();
		performCorrection();
		//if(fs::exists(_outputFilename_mapping)) fs::remove(_outputFilename_mapping);
	}

	void mapReads(){
		
		string readFilenames = "";
		ReadParser readParser(_inputFilename_reads, false, false);

		for(const string& filename : readParser._filenames){
			readFilenames += filename + " ";
		}

		string command = "minimap2 -t " + to_string(_nbCores) + " -x map-hifi " + _inputFilename_contigs + " " + readFilenames + " > " + _outputFilename_mapping;
		Utils::executeCommand(command);

	}

	struct Alignment{
		u_int32_t _contigIndex;
		u_int64_t _readIndex;
		bool _strand;
		u_int64_t _readStart;
		u_int64_t _readEnd;
		u_int64_t _contigStart;
		u_int64_t _contigEnd;
		float _score;
	};

	unordered_map<string, u_int32_t> _contigName_to_contigIndex;
	unordered_map<string, u_int64_t> _readName_to_readIndex;
	unordered_map<u_int64_t, Alignment> _alignments;
	vector<string> _contigSequences;
	vector<vector<vector<Window>>> _contigWindowSequences;
	unordered_map<ContigRead, u_int32_t, ContigRead_hash> _alignmentCounts;




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

	void parseAlignments(){

		cout << "Indexing read alignements" << endl;

        ifstream infile(_outputFilename_mapping);

        std::string line;
        vector<string>* fields = new vector<string>();
        //vector<string>* fields_optional = new vector<string>();


        while (std::getline(infile, line)){

            GfaParser::tokenize(line, fields, '\t');

			//cout << line << endl;

			const string& readName = (*fields)[0];
			const string& contigName = (*fields)[5];

			u_int64_t readStart = stoull((*fields)[2]);
			u_int64_t readEnd = stoull((*fields)[3]);
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
			Alignment align = {contigIndex, readIndex, strand, readStart, readEnd, contigStart, contigEnd, score};

			//ContigRead alignKey = {_contigName_to_contigIndex[contigName], _readName_to_readIndex[readName]};

			if(_alignments.find(readIndex) != _alignments.end()){
				const Alignment& existingAlignment = _alignments[readIndex];
				if(align._score > existingAlignment._score){
					_alignments[readIndex] = align;
				}
			}
			else{
				_alignments[readIndex] = align;
			}
			
			/*
			_alignmentCounts[{_contigName_to_contigIndex[contigName], _readName_to_readIndex[readName]}] += 1;

			if(_alignmentCounts[{_contigName_to_contigIndex[contigName], _readName_to_readIndex[readName]}] > 1){
				//cout << "multi map: " << contigName << " " << readName << " " << score << endl;
			}

			if(readName == "_92"){
				cout << contigName << " " << readName << " " << nbMatches << " " << alignLength << " " << contigStart << " " << contigEnd << " " << score << " " << score2<< endl;
			}
			if(readName == "_284"){
				cout << contigName << " " << readName << " " << nbMatches << " " << alignLength << " " << contigStart << " " << contigEnd << " " << score  << " " << score2<< endl;
			}
			if(readName == "_293"){
				cout << contigName << " " << readName << " " << nbMatches << " " << alignLength << " " << contigStart << " " << contigEnd << " " << score  << " " << score2<< endl;
			}
			if(readName == "_308"){
				cout << contigName << " " << readName << " " << nbMatches << " " << alignLength << " " << contigStart << " " << contigEnd << " " << score  << " " << score2<< endl;
			}
			if(readName == "_503"){
				cout << contigName << " " << readName << " " << nbMatches << " " << alignLength << " " << contigStart << " " << contigEnd << " " << score  << " " << score2<< endl;
			}
			if(readName == "_876"){
				cout << contigName << " " << readName << " " << nbMatches << " " << alignLength << " " << contigStart << " " << contigEnd << " " << score  << " " << score2<< endl;
			}
			*/
			//cout << readName << " " << contigName << " " << contigStart << " " << contigEnd << " " << readStart << " " << readEnd << " " << score << " " << strand << endl;
			//for(string field : (*fields)){
			//	cout << field << " ";
			//}
			//cout << endl;
        }
	}

	void loadContigs(){
		auto fp = std::bind(&ContigPolisher::loadContigs_read, this, std::placeholders::_1);
		ReadParser readParser(_inputFilename_contigs, true, false);
		readParser.parse(fp);
	}
	
	void loadContigs_read(const Read& read){
		_contigSequences.push_back(read._seq);

		size_t nbWindows = ceil((double)read._seq.size() / (double)_windowLength);
		vector<vector<Window>> windows(nbWindows);
		//cout << "Nb windows: " << nbWindows << endl;

		_contigWindowSequences.push_back(windows);
	}

	void collectWindowSequences(){
		
		cout << "Collecting window sequences" << endl;

		//auto fp = std::bind(&ContigPolisher::collectWindowCopies_read, this, std::placeholders::_1);
		//ReadParser readParser(_inputFilename_reads, false, false);
		//readParser.parse(fp);


		ReadParserParallel readParser(_inputFilename_reads, false, false, _nbCores);
		readParser.parse(CollectWindowSequencesFunctor(*this));
	}
	
	class CollectWindowSequencesFunctor {

		public:

		ContigPolisher& _contigPolisher;
		unordered_map<string, u_int32_t>& _contigName_to_contigIndex;
		unordered_map<string, u_int64_t>& _readName_to_readIndex;
		unordered_map<u_int64_t, Alignment>& _alignments;
		vector<string>& _contigSequences;
		vector<vector<vector<Window>>>& _contigWindowSequences;
		size_t _windowLength;


		CollectWindowSequencesFunctor(ContigPolisher& contigPolisher) : _contigPolisher(contigPolisher), _contigName_to_contigIndex(contigPolisher._contigName_to_contigIndex), _readName_to_readIndex(contigPolisher._readName_to_readIndex), _alignments(contigPolisher._alignments), _contigSequences(contigPolisher._contigSequences), _contigWindowSequences(contigPolisher._contigWindowSequences), _windowLength(contigPolisher._windowLength){
		}

		CollectWindowSequencesFunctor(const CollectWindowSequencesFunctor& copy) : _contigPolisher(copy._contigPolisher), _contigName_to_contigIndex(copy._contigName_to_contigIndex), _readName_to_readIndex(copy._readName_to_readIndex), _alignments(copy._alignments), _contigSequences(copy._contigSequences), _contigWindowSequences(copy._contigWindowSequences), _windowLength(copy._windowLength){
			
		}

		~CollectWindowSequencesFunctor(){
		}

		void operator () (const Read& read) {

			u_int64_t readIndex = read._index;
			
			if(readIndex % 100000 == 0) cout << readIndex << endl;

			if(_alignments.find(read._index) == _alignments.end()) return;

			const Alignment& al = _alignments[read._index];
			u_int64_t contigIndex = al._contigIndex;

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

				/*
				if (!sequence->quality().empty() ||
					!sequence->reverse_quality().empty()) {

					const auto& quality = overlaps[i]->strand() ?
						sequence->reverse_quality() : sequence->quality();
					double average_quality = 0;
					for (uint32_t k = breaking_points[j].second; k < breaking_points[j + 1].second; ++k) {
						average_quality += static_cast<uint32_t>(quality[k]) - 33;
					}
					average_quality /= breaking_points[j + 1].second - breaking_points[j].second;

					if (average_quality < quality_threshold_) {
						//continue;
					}
				}
				*/

				//uint64_t window_id = id_to_first_window_id[overlaps[i]->t_id()] +
				uint64_t window_id = breaking_points_[j].first / _windowLength;
				uint32_t window_start = (breaking_points_[j].first / _windowLength) *
					_windowLength;

				//const char* data = overlaps[i]->strand() ?
				//	&(sequence->reverse_complement()[breaking_points[j].second]) :
				//	&(sequence->data()[breaking_points[j].second]);
				const char* data = &readSequence[breaking_points_[j].second];
				uint32_t data_length = breaking_points_[j + 1].second -
					breaking_points_[j].second;

				string sequence = string(data, data_length);

                u_int32_t posStart = breaking_points_[j].first - window_start;
                u_int32_t posEnd =  breaking_points_[j + 1].first - window_start - 1;

				string quality = "";
				if(qualSequence.size() > 0){
					quality = string(&qualSequence[breaking_points_[j].second], data_length);
				}

				if(window_id == 0){
					//cout << sequence << endl;
					//cout << quality << endl;
					//getchar();
				}

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
				if(windows.size() < 20){


					windows.push_back({new DnaBitset2(windowSequence), windowQualities, posStart, posEnd});



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

		cout << "Perform correction" << endl;

		gzFile outputContigFile = gzopen(_outputFilename_contigs.c_str(),"wb");;

		abpt = abpoa_init_para();
		abpt->out_msa = 0; // generate Row-Column multiple sequence alignment(RC-MSA), set 0 to disable
		abpt->out_cons = 1; // generate consensus sequence, set 0 to disable
		abpt->w = 6, abpt->k = 9; abpt->min_w = 10; // minimizer-based seeding and partition
		abpt->progressive_poa = 1;
		abpt->max_n_cons = 1;

		abpoa_post_set_para(abpt);

		for(size_t contigIndex=0; contigIndex < _contigWindowSequences.size(); contigIndex++){

			vector<DnaBitset2*> correctedWindows(_contigWindowSequences[contigIndex].size());
		

			#pragma omp parallel num_threads(_nbCores)
			{

				std::unique_ptr<spoa::AlignmentEngine> alignmentEngine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -4);
				
				#pragma omp for
				for(size_t w=0; w<_contigWindowSequences[contigIndex].size(); w++){

					vector<Window>& sequences = _contigWindowSequences[contigIndex][w];
					
					u_int64_t wStart = w*_windowLength;
					u_int64_t wEnd = min(_contigSequences[contigIndex].size(), wStart+_windowLength);
					string contigOriginalSequence = _contigSequences[contigIndex].substr(wStart, wEnd-wStart);

					//cout << w << ": " << (sequences.size()+1) << endl;
					//continue;
					//cout << contigOriginalSequence << endl;
					//getchar();
					if(sequences.size() < 2){
							
						
						correctedWindows[w] = new DnaBitset2(contigOriginalSequence);
						//cout << "No sequences for window" << endl;
						continue;
					}
					
					std::sort(sequences.begin(), sequences.end(), []
					(Window& first, Window& second){
						return first._sequence->m_len > second._sequence->m_len;
					});

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

						//cout << dnaStr << endl;



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
						//}
					}





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

					//cout << correctedSequence << endl;
					//getchar();
					correctedWindows[w] = new DnaBitset2(correctedSequence);
				}
			}


			string contigSequence = "";
			for(size_t w=0; w<correctedWindows.size(); w++){
				if(correctedWindows[w] == nullptr) continue;
				char* seq = correctedWindows[w]->to_string();
				contigSequence += string(seq);
				free(seq);
				delete correctedWindows[w];
			}

			string header = ">ctg" + to_string(contigIndex) + '\n';
			gzwrite(outputContigFile, (const char*)&header[0], header.size());
			contigSequence +=  '\n';
			gzwrite(outputContigFile, (const char*)&contigSequence[0], contigSequence.size());
			//cout << contigSequence.size() << endl;
		}

		gzclose(outputContigFile);
	}

};	


#endif 


