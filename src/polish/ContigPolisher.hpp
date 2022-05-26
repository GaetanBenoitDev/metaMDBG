
/*
- ToBasespace: générer un assemblage grossier
- parallelisation
- stockage des windows via leur CIGAR et une seuqence de reference (verifier la taille des cigar sur données réelles)
- utiliser les scores de quality 
- quand un read a été process, on peut erase son entrée dans _alignments ?
- window sequence a selectionner en priorité: 
	- la distance est une mauvaise metrique car on ne sait pas si la sequence de reference est erroné ou non
	- un mapping tres long sur un contig a plus de valeur qu'un mapping court (on est plus sûr que ce read appartient au contig)
*/

#ifndef MDBG_METAG_CONTIGPOLISHER
#define MDBG_METAG_CONTIGPOLISHER

#include "../Commons.hpp"
#include "../utils/edlib.h"
//#include "../utils/spoa/include/spoa/spoa.hpp"
#include "../utils/DnaBitset.hpp"
#include "../utils/abPOA2/include/abpoa.h"



class ContigPolisher : public Tool{
    
public:

	string _inputFilename_reads;
	string _inputFilename_contigs;
	int _nbCores;
	size_t _windowLength;
	
	string _outputFilename_contigs;
	string _outputFilename_mapping;

	abpoa_para_t *abpt;
	/*
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
	*/


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
		_outputFilename_mapping= p.string() + "_tmp_mapping__.paf";


	}



    void execute (){

		mapReads();
		indexContigName();
		indexReadName();
		parseAlignments();

		_contigName_to_contigIndex.clear();
		_readName_to_readIndex.clear();

		loadContigs();
		collectWindowCopies();
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
	vector<vector<vector<DnaBitset*>>> _contigWindowSequences;
	//unordered_map<ContigRead, u_int32_t, ContigRead_hash> _alignmentCounts;




	void indexContigName(){
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
        	u_int64_t query_length = stoull((*fields)[1]);

			bool strand = (*fields)[4] == "-";
			float score = (double) nbMatches / (double) query_length;

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
			//_alignmentCounts[{_contigName_to_contigIndex[contigName], _readName_to_readIndex[readName]}] += 1;

			//if(_alignmentCounts[{_contigName_to_contigIndex[contigName], _readName_to_readIndex[readName]}] > 1){
			//	cout << "multi map: " << contigName << " " << readName << " " << score << endl;
			//}
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
		vector<vector<DnaBitset*>> windows(nbWindows);
		//cout << "Nb windows: " << nbWindows << endl;

		_contigWindowSequences.push_back(windows);
	}

	void collectWindowCopies(){
		
		cout << "Collecting window sequences" << endl;

		auto fp = std::bind(&ContigPolisher::collectWindowCopies_read, this, std::placeholders::_1);
		ReadParser readParser(_inputFilename_reads, false, false);
		readParser.parse(fp);
	}
	
	void collectWindowCopies_read(const Read& read){

		if(read._index % 100000 == 0) cout << read._index << endl;
		if(_alignments.find(read._index) == _alignments.end()) return;

		const Alignment& al = _alignments[read._index];
		u_int64_t contigIndex = al._contigIndex;

		string readSeq = read._seq;
		string readSequence = readSeq.substr(al._readStart, al._readEnd-al._readStart);
		string contigSequence = _contigSequences[contigIndex].substr(al._contigStart, al._contigEnd-al._contigStart);


		if(al._strand){
			Utils::toReverseComplement(readSequence);
			Utils::toReverseComplement(readSeq);
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
		
		find_breaking_points_from_cigar(_windowLength, al, readSeq.size(), cigar, readSeq);
		free(cigar);

		//getchar();
	}

	void find_breaking_points_from_cigar(uint32_t window_length, const Alignment& al, u_int64_t readLength, char* cigar_, const string& readSequence)
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
							//breaking_points_.emplace_back(last_match);
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
							//breaking_points_.emplace_back(last_match);
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

		if(breaking_points_.size() > 0) breaking_points_.emplace_back(last_match);
		
		for(size_t i=0; i<breaking_points_.size()-1; i++){
			u_int64_t contigWindowStart = breaking_points_[i].first;
			u_int64_t contigWindowEnd = breaking_points_[i+1].first;
			u_int64_t readWindowStart = breaking_points_[i].second;
			u_int64_t readWindowEnd = breaking_points_[i+1].second;

			//cout << contigWindowStart << " " << contigWindowEnd  << "      " << readWindowStart << " " << readWindowEnd << endl;

			//if(readWindowEnd-readWindowStart < _windowLength) continue; //window sides
			//string windowSequence = 
			indexWindow(al, contigWindowStart, contigWindowEnd, readSequence.substr(readWindowStart, readWindowEnd-readWindowStart));
		}
		

		//for(const auto& breakPoint : breaking_points_){
		//	cout << breakPoint.first << " " << breakPoint.second << endl;
		//}
	}

	void indexWindow(const Alignment& al, size_t contigWindowStart, size_t contigWindowEnd, const string& windowSequence){

		size_t contigWindowIndex = contigWindowStart / _windowLength;
		vector<DnaBitset*>& windowSequences = _contigWindowSequences[al._contigIndex][contigWindowIndex];
		
		if(windowSequences.size() < 20){
			windowSequences.push_back(new DnaBitset(windowSequence));

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


			return;
		}


		size_t largerWindowIndex = 0;
		u_int64_t largerDistanceWindow = 0;

		for(size_t i=0; i<windowSequences.size(); i++){

			DnaBitset* dnaSeq = windowSequences[i];
			u_int64_t distance = abs(((long)dnaSeq->m_len) - ((long)_windowLength));

			if(distance > largerDistanceWindow){
				largerDistanceWindow = distance;
				largerWindowIndex = i;
			}
		}


		u_int64_t distance = abs(((long)windowSequence.size()) - ((long)_windowLength));

		if(distance < largerDistanceWindow){
			DnaBitset* dnaSeq = windowSequences[largerWindowIndex];
			delete dnaSeq;
			windowSequences[largerWindowIndex] = new DnaBitset(windowSequence);
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

	void performCorrection(){

		cout << "Perform correction" << endl;

		gzFile outputContigFile = gzopen(_outputFilename_contigs.c_str(),"wb");;

		abpt = abpoa_init_para();
		abpt->out_msa = 1; // generate Row-Column multiple sequence alignment(RC-MSA), set 0 to disable
		abpt->out_cons = 0; // generate consensus sequence, set 0 to disable
		abpt->w = 6, abpt->k = 9; abpt->min_w = 10; // minimizer-based seeding and partition
		abpt->progressive_poa = 1;
		abpt->max_n_cons = 2;

		abpoa_post_set_para(abpt);

		for(size_t contigIndex=0; contigIndex < _contigWindowSequences.size(); contigIndex++){

			string contigSequence = "";

			for(size_t w=0; w<_contigWindowSequences[contigIndex].size(); w++){

				const vector<DnaBitset*>& sequences = _contigWindowSequences[contigIndex][w];
				if(sequences.size() == 0){
					cout << "No sequences for window" << endl;
					continue;
				}
				//vector<u_int32_t> windowLengths;
				//for(size_t i=0; i<sequences.size(); i++){
				//	windowLengths.push_back(sequences[i]->m_len);
				//}
				//cout << Utils::compute_median(windowLengths) << endl;


				abpoa_t *ab = abpoa_init();

				//vector<size_t> order;
				//for(size_t i=0; i<sequences.size(); i++){
				//	order.push_back(i);
				//}
				//srand(time(NULL));
				//std::random_shuffle(order.begin(), order.end());

				/*
				vector<string> seqSorted;
				for(size_t i=1; i<sequences.size(); i++){
					DnaBitset* dna = sequences[i];
					//const DnaBitset* dna = variant._sequence; //sequenceCopies[s._sequenceIndex];
					char* dnaStr = dna->to_string();
					seqSorted.push_back(string(dnaStr));
					free(dnaStr);
				}


				std::sort(seqSorted.begin(), seqSorted.end());

				DnaBitset* dnaModel = sequences[0];
				//const DnaBitset* dna = variant._sequence; //sequenceCopies[s._sequenceIndex];
				char* dnaStrModel = dnaModel->to_string();
				seqSorted.insert(seqSorted.begin(), string(dnaStrModel));
				free(dnaStrModel);
				*/

				//cout << "1" << endl;
				int n_seqs = sequences.size();
				int *seq_lens = (int*)malloc(sizeof(int) * n_seqs);
				uint8_t **bseqs = (uint8_t**)malloc(sizeof(uint8_t*) * n_seqs);
				
				for(size_t i=0; i<sequences.size(); i++){ 

					//size_t i = order[ii];
					DnaBitset* dna = sequences[i];
					//const DnaBitset* dna = variant._sequence; //sequenceCopies[s._sequenceIndex];
					char* dnaStr = dna->to_string();

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

				
				abpoa_msa(ab, abpt, n_seqs, NULL, seq_lens, bseqs, NULL);


				
				abpoa_cons_t *abc = ab->abc;
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

				contigSequence += correctedSequence;

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


