

#ifndef MDBG_METAG_DEREPLICATER
#define MDBG_METAG_DEREPLICATER

#include "../Commons.hpp"

class Dereplicater : public Tool{
    
public:

	string _inputFilename_contigs;
	string _outputFilename_contigs;
	string _tmpDir;
	int _nbCores;

	gzFile _outputContigFile;

	Dereplicater(): Tool (){


	}

	void parseArgs(int argc, char* argv[]){


		cxxopts::Options options("ToBasespace", "");
		options.add_options()
		("contigs", "", cxxopts::value<string>())
		("outputFilenme", "", cxxopts::value<string>())
		("tmpDir", "", cxxopts::value<string>())
		(ARG_NB_CORES, "", cxxopts::value<int>()->default_value(NB_CORES_DEFAULT));

		options.parse_positional({"contigs", "outputFilenme", "tmpDir"});
		options.positional_help("contigs outputFilenme tmpDir");

		if(argc <= 1){
			cout << options.help() << endl;
			exit(0);
		}

		cxxopts::ParseResult result;

		try{
			result = options.parse(argc, argv);

			_inputFilename_contigs = result["contigs"].as<string>();
			_outputFilename_contigs = result["outputFilenme"].as<string>();
			_tmpDir = result["tmpDir"].as<string>();;
			_nbCores = result[ARG_NB_CORES].as<int>();
		}
		catch (const std::exception& e){
			std::cout << options.help() << std::endl;
			std::cerr << e.what() << std::endl;
			std::exit(EXIT_FAILURE);
		}

		if (_outputFilename_contigs.find(".gz") == std::string::npos) {
			_outputFilename_contigs += ".gz";
		}

		_tmpDir = _tmpDir + "/__tmp/";
		if(!fs::exists(_tmpDir)){
			fs::create_directories(_tmpDir);
		}

		cout << endl;
		cout << "Contigs: " << _inputFilename_contigs << endl;
		cout << "Output contigs: " << _outputFilename_contigs << endl;
		cout << "Tmp dir: " << _tmpDir << endl;
		cout << endl;

	}

	struct ContigOverlap{
		u_int64_t _leftPos;
		u_int64_t _rightPos;

		ContigOverlap(){
			_leftPos = 0;
			_rightPos = -1;
		}
	};

	unordered_map<u_int64_t, ContigOverlap> _contigOverlaps;

    void execute (){
		
		string outputMappingFilename = _tmpDir + "/derep_align.paf";

		mapContigs(outputMappingFilename);
		detectDuplicatedContigs(outputMappingFilename);
		dumpDereplicatedContigs();
		

		//fs::remove(outputMappingFilename);
	}


	//zFile _queryContigFile;
	unordered_set<u_int32_t> _duplicatedContigIndex;

	void mapContigs(const string& mapFilename){



		string command = "minimap2 -N 2 -p 0 --secondary=yes -I 2GB -x map-hifi -t " + to_string(_nbCores) + " " + _inputFilename_contigs + " " + _inputFilename_contigs + " > " + mapFilename;
		Utils::executeCommand(command, _tmpDir);

	}

	void detectDuplicatedContigs(const string& mapFilename){

		ifstream mappingFile(mapFilename);
        vector<string>* fields = new vector<string>();
        vector<string>* fields_optional = new vector<string>();

		string line;
		while (getline(mappingFile, line)) {


            GfaParser::tokenize(line, fields, '\t');

			const string& readName = (*fields)[0];
			const string& contigName = (*fields)[5];
			if(readName == contigName) continue;

			//string search1 ="ctg88963";
			//string search2 ="ctg89401";
			//string search1 ="ctg89096";
			//string search2 ="ctg90131";
			//string search1 ="ctg87283";
			//string search2 ="ctg90059";
			//string search1 ="ctg90131";
			//string search2 ="ctg88674";
			string search1 ="ctg88963";
			string search2 ="ctg89401";
			



			//continue;
			if((readName == search1 || contigName == search1) && (readName == search2 || contigName == search2)){
				cout << line << endl;
			}



			//continue;

			//u_int64_t queryLength = stoull((*fields)[1]);
			u_int64_t targetLength = stoull((*fields)[6]);
			double queryLength = stoull((*fields)[1]);

			//if(queryLength < targetLength) continue;
			//if(targetLength > queryLength) continue;


			u_int64_t queryStart = stoull((*fields)[2]);
			u_int64_t queryEnd = stoull((*fields)[3]);
			u_int64_t targetStart = stoull((*fields)[7]);
			u_int64_t targetEnd = stoull((*fields)[8]);
			bool strand = false;
			if((*fields)[4] == "-") strand = true;

			//u_int64_t ql = queryLength;
			//u_int64_t tl = targetLength;
			//u_int64_t qs = queryStart;
			//u_int64_t qe = queryEnd;
			//u_int64_t ts = targetStart;
			//u_int64_t te = targetEnd;

			u_int64_t nbMatches = stoull((*fields)[9]);
			double alignLength = stoull((*fields)[10]);
			//u_int64_t ml = nbMatches;
			//u_int64_t bl = alignLength;


			if(nbMatches / alignLength < 0.2) continue;
			if(alignLength < 2000) continue;

			u_int64_t maxHang = 1000;
			u_int64_t hangLeft = targetStart;
			u_int64_t hangRight = targetLength - targetEnd;

			u_int64_t targetIndex = Utils::contigName_to_contigIndex(contigName);

			//if(readName == "ctg89401" || contigName == "ctg89401"){
				//cout << line << endl;
			//}
			
			if(hangLeft < maxHang){

				cout << "left overlap: " << targetIndex << "    " << line << endl;
				//cout << "left overlap: " << targetIndex << " " << contigOverlap._leftPos << endl;
				
				if(_contigOverlaps.find(targetIndex) == _contigOverlaps.end()){
					_contigOverlaps[targetIndex] = ContigOverlap();
				}
				ContigOverlap& contigOverlap = _contigOverlaps[targetIndex];

				//if(strand){
				//	contigOverlap._rightPos = min(contigOverlap._rightPos, targetLength-targetEnd);
					//contigOverlap._rightPos = max(contigOverlap._leftPos, targetStart);
				//}
				//else{
					contigOverlap._leftPos = max(contigOverlap._leftPos, targetEnd);
				//}

				//_contigOverlaps[targetIndex] = contigOverlap;
			}

			if(hangRight < maxHang){

				cout << "right overlap: " << targetIndex << "    " << line << endl;
				//cout << "right overlap: " << targetIndex << " " << contigOverlap._leftPos << endl;

				if(_contigOverlaps.find(targetIndex) == _contigOverlaps.end()){
					_contigOverlaps[targetIndex] = ContigOverlap();
				}
				ContigOverlap& contigOverlap = _contigOverlaps[targetIndex];

				//if(strand){
					//contigOverlap._leftPos = max(contigOverlap._leftPos, targetLength-targetStart);
					//contigOverlap._rightPos = min(contigOverlap._rightPos, targetEnd);
				//}
				//else{
					contigOverlap._rightPos = min(contigOverlap._rightPos, targetStart);
					//contigOverlap._rightPos = min(contigOverlap._rightPos, targetStart);
				//}
				
				//_contigOverlaps[targetIndex] = contigOverlap;
			}

			//cout << alignLength / queryLength << endl;



			/*
			if(alignLength / queryLength < 0.1) continue;

			//cout << (nbMatches / alignLength) << " " << (nbMatches / queryLength) << endl;
			
			float divergence = 1;

			for(size_t i=12; i<fields->size(); i++){

				//cout << (*fields)[i] << endl;

				GfaParser::tokenize((*fields)[i], fields_optional, ':');

				if((*fields_optional)[0] == "dv"){
					divergence = std::stof((*fields_optional)[2]);

				
					
					break;
				}

			}
			*/
			/*
			u_int64_t queryIndex = Utils::contigName_to_contigIndex(readName);
			u_int64_t targetIndex = Utils::contigName_to_contigIndex(contigName);

			u_int64_t min_span = 2000;
			u_int64_t min_match = 100;
			u_int64_t max_hang = 1000;
			float int_frac = 0.8;
			
			if(alignLength / queryLength > 1){
				cout << "lala" << endl;
				cout << line << endl;
			}

			u_int64_t maxHang = 10000;


			float queryRatio = (queryEnd-queryStart) / queryLength;
			//if(queryRatio > 0.5){
			//	cout << "overlap: " << line << endl; 
			//	_duplicatedContigIndex.insert(queryIndex);
			//	continue;
			//}

			//if(strand){
				
				u_int64_t leftHang = queryStart;
				u_int64_t rightHang_target = queryLength - queryEnd;
				u_int64_t rightHand_query = targetLength - targetEnd;
				u_int64_t rightHang = min(rightHang_target, rightHand_query);
				
				if(leftHang < maxHang && rightHang < maxHang){
					cout << "overlap: " << line << endl; 
					_duplicatedContigIndex.insert(queryIndex);
				}
			*/
			/*
			}
			else{

				u_int64_t leftHang = queryLength - queryEnd;
				u_int64_t rightHang_target = queryStart;
				u_int64_t rightHand_query = targetLength - targetEnd;
				u_int64_t rightHang = min(rightHang_target, rightHand_query);

				if(leftHang < maxHang && rightHang < maxHang){
					cout << "overlap: " << line << endl; 
					_duplicatedContigIndex.insert(queryIndex);
				}

			}
			*/

			/*
			int l5, l3;
			if (qe - qs < min_span || te - ts < min_span || ml < min_match){
				cout << "pre" << endl;
				continue;
			}
			l5 = (!strand)? tl - te : ts;
			l3 = (!strand)? ts : tl - te;
			if (ql>>1 > tl) {
				//cout << (l5 > max_hang>>2) << " " << (l3 > max_hang>>2) << " " << (te - ts < tl * int_frac) << endl;
				if (l5 > max_hang>>2 || l3 > max_hang>>2 || te - ts < tl * int_frac){
					cout << line << endl;
					cout << "internal 1" << endl;
					continue; // internal match
				}
				if ((int)qs - l5 > max_hang<<1 && (int)(ql - qe) - l3 > max_hang<<1){
					//_duplicatedContigIndex.insert();
					cout << "contained 1" << endl;
					//sd_put(d, r.tn, r.tl);
					//continue;
				}
			} else if (ql < tl>>1) {
				//cout << qs << " " << (max_hang>>2) << " " << (qs > max_hang>>2) << " " << (qe - qs < ql * int_frac) << endl;
				if (qs > max_hang>>2 || ql - qe > max_hang>>2 || qe - qs < ql * int_frac){
					cout << line << endl;
					cout << "internal 1" << endl;
					continue; // internal match
				}
				if (l5 - (int)qs > max_hang<<1 && l3 - (int)(ql - qe) > max_hang<<1){
					cout << "contained 1" << endl;
					//sd_put(d, r.qn, r.ql);
					//continue;
				} 
			}
			


			//uint64_t qns = sd_put(d, r.qn, r.ql)<<32 | r.qs
			
			int32_t tl5, tl3, ext5, ext3;//, qs = (int32_t)h->qns;
			uint32_t u, v, l; // u: query end; v: target end; l: length from u to v
			if (!strand) tl5 = tl - te, tl3 = ts; // tl5: 5'-end overhang (on the query strand); tl3: similar
			else tl5 = ts, tl3 = tl - te;
			ext5 = qs < tl5? qs : tl5;
			ext3 = ql - qe < tl3? ql - qe : tl3;
			if (ext5 > max_hang || ext3 > max_hang){// || qe - qs < (qe - qs + ext5 + ext3) * int_frac)
				cout << "max hang" << endl;
				continue;
			}
			if (qs <= tl5 && ql - qe <= tl3){
				_duplicatedContigIndex.insert(queryIndex);
				cout << "contained 2" << endl;
				continue; // query contained	
			} 
			else if (qs >= tl5 && ql - qe >= tl3){
				_duplicatedContigIndex.insert(targetIndex);
				cout << "contained 2" << endl;
				continue; // target contained
			}
			else if (qs > tl5) u = 0, v = !!(!strand), l = qs - tl5;
			else u = 1, v = !(!strand), l = (ql - qe) - tl3;
			//if (qe - qs + ext5 + ext3 < min_ovlp || h->te - h->ts + ext5 + ext3 < min_ovlp) return MA_HT_SHORT_OVLP; // short overlap
			//u |= h->qns>>32<<1, v |= h->tn<<1;
			//ul = (uint64_t)u<<32 | l, p->v = v, p->ol = ql - l, p->del = 0;


			_duplicatedContigIndex.insert(targetIndex);

			cout << line << endl;
			cout << "overlap" << endl;
			*/
		}
		
		mappingFile.close();

	}
	
	void dumpDereplicatedContigs(){
		
		cout << "Writing dereplicated contigs" << endl;

		_outputContigFile = gzopen(_outputFilename_contigs.c_str(),"wb");

		string s1 = "/home/gats/workspace/tmp/file1.fa";
		string s2 = "/home/gats/workspace/tmp/file2.fa";
		string s3 = "/home/gats/workspace/tmp/file3.fa";
		string s4 = "/home/gats/workspace/tmp/file4.fa";
		_file1 = gzopen(s1.c_str(),"wb");
		_file2 = gzopen(s2.c_str(),"wb");
		_file3 = gzopen(s3.c_str(),"wb");
		_file4 = gzopen(s4.c_str(),"wb");

		auto fp = std::bind(&Dereplicater::dumpDereplicatedContigs_read, this, std::placeholders::_1);
		ReadParser readParser(_inputFilename_contigs, true, false);
		readParser.parse(fp);

		gzclose(_outputContigFile);
		gzclose(_file1);
		gzclose(_file2);
		gzclose(_file3);
		gzclose(_file4);
	}

	gzFile _file1;
	gzFile _file2;
	gzFile _file3;
	gzFile _file4;

	void dumpDereplicatedContigs_read(const Read& read){

		if(read._header == "ctg88963"){
			string header = ">" + read._header + '\n';
			gzwrite(_file1, (const char*)&header[0], header.size());
			string contigSequence = read._seq + '\n';
			gzwrite(_file1, (const char*)&contigSequence[0], contigSequence.size());

			/*
			string seq1 = read._seq;
			seq1 = seq1.substr(0, 167682+30000);
			string seq2 = read._seq;
			seq2 = seq2.substr(167682+30000, seq2.size());

			
			header = ">" + read._header + "_1" + '\n';
			gzwrite(_file3, (const char*)&header[0], header.size());
			contigSequence = seq1 + '\n';
			gzwrite(_file3, (const char*)&contigSequence[0], contigSequence.size());
			
			header = ">" + read._header + "_2" + '\n';
			gzwrite(_file4, (const char*)&header[0], header.size());
			contigSequence = seq2 + '\n';
			gzwrite(_file4, (const char*)&contigSequence[0], contigSequence.size());
			*/
		}
		else if(read._header == "ctg89401"){
			string header = ">" + read._header + '\n';
			gzwrite(_file2, (const char*)&header[0], header.size());
			string contigSequence = read._seq + '\n';
			gzwrite(_file2, (const char*)&contigSequence[0], contigSequence.size());

			string seq1 = read._seq;
			seq1 = seq1.substr(0, 31159);
			string seq2 = read._seq;
			seq2 = seq2.substr(31159, seq2.size());

			
			header = ">" + read._header + "_1" + '\n';
			gzwrite(_file3, (const char*)&header[0], header.size());
			contigSequence = seq1 + '\n';
			gzwrite(_file3, (const char*)&contigSequence[0], contigSequence.size());
			
			header = ">" + read._header + "_2" + '\n';
			gzwrite(_file4, (const char*)&header[0], header.size());
			contigSequence = seq2 + '\n';
			gzwrite(_file4, (const char*)&contigSequence[0], contigSequence.size());
		}

		u_int64_t contigIndex = Utils::contigName_to_contigIndex(read._header);
		//if(_duplicatedContigIndex.find(contigIndex) != _duplicatedContigIndex.end()){
		//	cout << "Discard: " << read._seq.size() << endl;
		//	return;
		//}

		string seq = read._seq;


		if(_contigOverlaps.find(contigIndex) != _contigOverlaps.end()){

			const ContigOverlap& ov = _contigOverlaps[contigIndex];

			cout << contigIndex << " " << ov._leftPos << " " << ov._rightPos << endl;
			if(ov._leftPos != 0 && ov._rightPos != -1){
				if(ov._leftPos >= ov._rightPos){
					cout << "\tSlicing complete: " << read._seq.size() << endl;
					return;
				}
				else{
					//u_int64_t rightLength = seq.size()-ov._rightPos;
					seq = seq.substr(ov._leftPos, ov._rightPos-ov._leftPos);
					cout << "\tSlicing both: " << read._seq.size() << endl;
				}
			}
			else if(ov._leftPos != 0){
				seq = seq.substr(ov._leftPos, seq.size()-ov._leftPos);
				cout << "\tSlicing left: " << read._seq.size() << endl;
			}
			else{
				//u_int64_t rightLength = seq.size()-ov._rightPos;
				seq = seq.substr(ov._leftPos, ov._rightPos-ov._leftPos);
				cout << "\tSlicing right: " << read._seq.size() << endl;
			}
		}

		string header = ">" + read._header + '\n';
		gzwrite(_outputContigFile, (const char*)&header[0], header.size());
		string contigSequence = seq + '\n';
		gzwrite(_outputContigFile, (const char*)&contigSequence[0], contigSequence.size());

	}
};	







/*

90059 29513 18446744073709551615
	Slicing left: 847077
89096 53191 168468
	Slicing both: 295145
89401 3705 18446744073709551615
	Slicing left: 380917
89567 67242 18446744073709551615
	Slicing left: 448801
86936 124326 13
	Slicing complete: 124337
88674 0 185221
	Slicing right: 224475
83777 64387 59
	Slicing complete: 64390
85082 81257 56
	Slicing complete: 81344
87283 0 44197
	Slicing right: 134520
88963 0 83159
	Slicing right: 250826
87938 154041 18446744073709551615
	Slicing left: 156245
84539 74709 10
	Slicing complete: 74720
89441 369573 18446744073709551615
	Slicing left: 399063
61502 3324 9961
	Slicing both: 13285










bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	TIGR00488	ctg90059_919
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	TIGR00488	ctg89401_330

bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	TIGR01351	ctg90131_610
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	TIGR00244	ctg90059_491
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	TIGR00416	ctg89401_93
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	TIGR00691	ctg89401_73
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	TIGR03598	ctg89401_211
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	TIGR00182	ctg90059_766
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	TIGR01575	ctg89401_313
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	TIGR00460	ctg90059_905
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	TIGR03263	ctg89401_338
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	TIGR03263	ctg90059_911
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	TIGR00436	ctg89401_260
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	TIGR00154	ctg90131_816
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	TIGR00422	ctg90131_554
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	TIGR00967	ctg90131_611
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	TIGR01079	ctg90131_621
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	TIGR00406	ctg90059_497
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	TIGR00159	ctg89401_276
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	TIGR00360	ctg89401_311
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	TIGR00185	ctg90059_729
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	TIGR03596	ctg90059_750
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	TIGR03725	ctg89401_314
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	TIGR03725	ctg90059_944
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	TIGR01051	ctg90059_741
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	TIGR00595	ctg90059_908
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	TIGR02168	ctg90059_763
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	TIGR03974	ctg87283_65
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	TIGR03974	ctg90059_284
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	TIGR00634	ctg90059_459
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	TIGR00482	ctg90059_919
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	TIGR00482	ctg89401_330
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	TIGR03594	ctg87283_117
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	TIGR03594	ctg90059_361
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	TIGR01069	ctg87283_83
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	TIGR01510	ctg90059_779
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	TIGR00233	ctg89401_83
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF03602.10	ctg90059_780
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF04816.7	ctg90059_830
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00366.15	ctg90131_623
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF10458.4	ctg90131_554
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF01649.13	ctg90059_824
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF04998.12	ctg90131_870
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF03726.9	ctg90131_385
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF03947.13	ctg90131_629
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00181.18	ctg90131_629
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF05833.6	ctg90059_893
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF01715.12	ctg90059_353
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF01715.12	ctg87283_111
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF09285.6	ctg90131_359
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF08207.7	ctg90131_359
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF01132.15	ctg90131_359
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF01000.21	ctg90131_601
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF01193.19	ctg90131_601
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00861.17	ctg90131_615
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF01782.13	ctg90059_753
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00297.17	ctg90131_632
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00318.15	ctg90131_427
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00416.17	ctg90131_604
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00831.18	ctg90131_624
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF08459.6	ctg89401_132
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF01351.13	ctg90059_747
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF05190.13	ctg90059_356
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF05190.13	ctg87283_113
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF05188.12	ctg90059_356
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF05188.12	ctg87283_113
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF08676.6	ctg90059_354
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF08676.6	ctg87283_112
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF01119.14	ctg87283_112
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF01119.14	ctg90059_355
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00410.14	ctg90131_618
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF06135.7	ctg90131_431
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF06838.6	ctg90059_351
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF06838.6	ctg87283_109
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF13742.1	ctg90059_472
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF02021.12	ctg90059_746
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF02132.10	ctg89401_86
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF01245.15	ctg90059_751
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00562.23	ctg90131_872
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF04563.10	ctg90131_872
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF04561.9	ctg90131_872
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF04560.15	ctg90131_872
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF04565.11	ctg90131_872
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF10385.4	ctg90131_872
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF02542.11	ctg89401_241
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF02224.13	ctg87283_92
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF02224.13	ctg90059_330
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF01121.15	ctg90059_806
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF01746.16	ctg90059_752
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00453.13	ctg90131_340
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF06421.7	ctg90059_504
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF11987.3	ctg90131_395
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF06574.7	ctg90131_390
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF06574.7	ctg90131_302
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF06574.7	ctg90059_35&&ctg90059_36
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF01687.12	ctg90131_302
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF01687.12	ctg90059_36
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF02401.13	ctg90059_332
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF02401.13	ctg87283_94
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00338.17	ctg90131_633
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF01632.14	ctg90131_341
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF01765.14	ctg90131_422
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00889.14	ctg90131_426
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF12072.3	ctg87283_100
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF12072.3	ctg90059_340
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF13184.1	ctg90131_400
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF04997.7	ctg90131_871
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00623.15	ctg90131_871
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF04983.13	ctg90131_871
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF05000.12	ctg90131_871
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF05198.11	ctg90131_342
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00707.17	ctg90131_342
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF04079.11	ctg90059_449
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF01078.16	ctg90059_743
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF01078.16	ctg89401_209
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF05491.8	ctg87283_56
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF05491.8	ctg90059_275
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF03147.9	ctg90131_50
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF03147.9	ctg90059_195
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF03484.10	ctg90131_50
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF03484.10	ctg90059_195
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF04297.9	ctg90059_758
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00673.16	ctg90131_620
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00281.14	ctg90131_620
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF01628.16	ctg90059_502
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00189.15	ctg90131_626
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00380.14	ctg90059_727
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF02130.12	ctg89401_258
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF08529.6	ctg90131_402
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF03485.11	ctg90131_1203
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00154.16	ctg87283_99
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00154.16	ctg90059_339
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF05697.8	ctg89401_204
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00828.14	ctg90131_612
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF06071.8	ctg90131_1054
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00237.14	ctg90131_627
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF01205.14	ctg90131_1123
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF09186.6	ctg90131_1123
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF09269.6	ctg90059_922
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF09269.6	ctg89401_327
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF01018.17	ctg89401_327
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF01018.17	ctg90059_922&&ctg90059_923
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF04468.7	ctg90131_844
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF07499.8	ctg87283_55
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF07499.8	ctg90059_274
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF01250.12	ctg90131_1145
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF01084.15	ctg90131_1147
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF06969.11	ctg90059_503
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00118.19	ctg90131_641
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF12344.3	ctg90131_261
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF12344.3	ctg90059_77
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF03461.10	ctg89401_238
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF03668.10	ctg90059_426
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF02565.10	ctg89401_262
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF02609.11	ctg90059_470
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00411.14	ctg90131_603
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF02616.9	ctg90059_450
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF06144.8	ctg90059_949
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF12392.3	ctg90059_314
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF12392.3	ctg87283_82
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00347.18	ctg90131_616&&ctg90131_617
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00573.17	ctg90131_631
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF02912.13	ctg90059_194
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF02912.13	ctg90131_51
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00886.14	ctg90059_755
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF03652.10	ctg90131_430
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF08275.6	ctg90059_820
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF01176.14	ctg90131_606
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF02650.9	ctg90059_424
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF02603.11	ctg89401_134
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF01725.11	ctg90059_891
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF02631.11	ctg90059_785
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00750.14	ctg90131_1204
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00238.14	ctg90131_622
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00252.13	ctg90131_625
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF01196.14	ctg90131_599
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF02367.12	ctg90059_943
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF02367.12	ctg89401_315
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF01195.14	ctg89401_236
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF03719.10	ctg90131_614
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00333.15	ctg90131_614
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00276.15	ctg90131_630
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF02547.10	ctg89401_317
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF02547.10	ctg90059_941
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF02033.13	ctg90131_394
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF03948.9	ctg90131_1141
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF01281.14	ctg90131_1141
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF02575.11	ctg89401_85
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF02576.12	ctg90131_404
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF02601.10	ctg90059_471
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF02978.14	ctg90059_756
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00164.20	ctg90131_867
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00572.13	ctg90059_728
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF01509.13	ctg90131_391&&ctg90131_392
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF02873.11	ctg90059_427
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00203.16	ctg90131_628
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF01702.13	ctg89401_316
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF01702.13	ctg90059_942
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF01016.14	ctg90131_246
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF01016.14	ctg90059_86
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00312.17	ctg90131_389
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00829.16	ctg90131_244
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00829.16	ctg90059_88
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF00177.16	ctg90131_866
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF01709.15	ctg90131_216
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF05636.6	ctg90059_776
bin.1081.derep.derep.derep.derep.derep.derep.derep.derep	PF04296.8	ctg90131_399









bin.1081.derep	TIGR00615	ctg89401_85
bin.1081.derep	TIGR03723	ctg89567_174
bin.1081.derep	TIGR00855	ctg89567_165
bin.1081.derep	TIGR00250	ctg90131_426
bin.1081.derep	TIGR00250	ctg88674_165
bin.1081.derep	TIGR02432	ctg89096_19
bin.1081.derep	TIGR02432	ctg90131_1042
bin.1081.derep	TIGR02075	ctg88674_161
bin.1081.derep	TIGR00922	ctg89567_169
bin.1081.derep	TIGR00392	ctg89567_31
bin.1081.derep	TIGR03263	ctg90059_854
bin.1081.derep	TIGR00329	ctg89567_174
bin.1081.derep	TIGR00459	ctg89401_77
bin.1081.derep	TIGR00755	ctg89567_48
bin.1081.derep	TIGR01079	ctg90131_615
bin.1081.derep	TIGR00344	ctg89567_77
bin.1081.derep	TIGR00019	ctg89567_277
bin.1081.derep	TIGR03594	ctg90059_307
bin.1081.derep	TIGR00967	ctg90131_605
bin.1081.derep	TIGR00460	ctg90059_848
bin.1081.derep	TIGR00460	ctg89441_20
bin.1081.derep	TIGR00084	ctg90059_221
bin.1081.derep	PF00366.15	ctg90131_617
bin.1081.derep	PF04997.7	ctg90131_863
bin.1081.derep	PF04997.7	ctg89096_90&&ctg89096_91
bin.1081.derep	PF00623.15	ctg90131_863
bin.1081.derep	PF00623.15	ctg89096_91
bin.1081.derep	PF04983.13	ctg90131_863
bin.1081.derep	PF04983.13	ctg89096_91
bin.1081.derep	PF05000.12	ctg90131_863
bin.1081.derep	PF05000.12	ctg89096_91
bin.1081.derep	PF04998.12	ctg89096_91
bin.1081.derep	PF04998.12	ctg90131_862

bin.1081.derep	PF01795.14	ctg89401_214
bin.1081.derep	PF01795.14	ctg88963_1

bin.1081.derep	PF01016.14	ctg90131_243
bin.1081.derep	PF01016.14	ctg90059_38
bin.1081.derep	PF12344.3	ctg90131_258
bin.1081.derep	PF00886.14	ctg90059_699
bin.1081.derep	PF01250.12	ctg89567_354
bin.1081.derep	PF01250.12	ctg90131_1135
bin.1081.derep	PF08529.6	ctg90131_398
bin.1081.derep	PF08529.6	ctg88674_151
bin.1081.derep	PF13184.1	ctg88674_151
bin.1081.derep	PF13184.1	ctg90131_396
bin.1081.derep	PF06421.7	ctg90059_449
bin.1081.derep	PF11987.3	ctg88674_148
bin.1081.derep	PF11987.3	ctg90131_391
bin.1081.derep	PF00889.14	ctg90131_422
bin.1081.derep	PF00889.14	ctg88674_162
bin.1081.derep	PF00562.23	ctg90131_864
bin.1081.derep	PF00562.23	ctg89096_89
bin.1081.derep	PF04563.10	ctg90131_864
bin.1081.derep	PF04563.10	ctg89096_89
bin.1081.derep	PF04561.9	ctg90131_864
bin.1081.derep	PF04561.9	ctg89096_89
bin.1081.derep	PF04560.15	ctg90131_864
bin.1081.derep	PF04560.15	ctg89096_89
bin.1081.derep	PF04565.11	ctg90131_864
bin.1081.derep	PF04565.11	ctg89096_89
bin.1081.derep	PF10385.4	ctg90131_864
bin.1081.derep	PF10385.4	ctg89096_89
bin.1081.derep	PF01509.13	ctg88674_145
bin.1081.derep	PF01509.13	ctg90131_387&&ctg90131_388
bin.1081.derep	PF05491.8	ctg90059_222
bin.1081.derep	PF00203.16	ctg90131_622
bin.1081.derep	PF02978.14	ctg90059_700

bin.1081.derep	PF01018.17	ctg89401_326
bin.1081.derep	PF01018.17	ctg90059_865&&ctg90059_866

bin.1081.derep	PF06071.8	ctg89096_18
bin.1081.derep	PF06071.8	ctg90131_1044
bin.1081.derep	PF00573.17	ctg90131_625
bin.1081.derep	PF00416.17	ctg90131_598
bin.1081.derep	PF00411.14	ctg90131_597
bin.1081.derep	PF00861.17	ctg90131_609
bin.1081.derep	PF01765.14	ctg88674_160
bin.1081.derep	PF01765.14	ctg90131_418
bin.1081.derep	PF01649.13	ctg90059_767
bin.1081.derep	PF00189.15	ctg90131_620
bin.1081.derep	PF00380.14	ctg90059_672
bin.1081.derep	PF00177.16	ctg89096_94
bin.1081.derep	PF00177.16	ctg90131_858
bin.1081.derep	PF03719.10	ctg90131_608
bin.1081.derep	PF00333.15	ctg90131_608
bin.1081.derep	PF01632.14	ctg88674_114
bin.1081.derep	PF01632.14	ctg90131_338

bin.1081.derep	PF02367.12	ctg90059_886
bin.1081.derep	PF02367.12	ctg89401_314

bin.1081.derep	PF01409.15	ctg90059_144
bin.1081.derep	PF01409.15	ctg90131_51
bin.1081.derep	PF02912.13	ctg90059_144
bin.1081.derep	PF02912.13	ctg90131_51
bin.1081.derep	PF03948.9	ctg90131_1131
bin.1081.derep	PF03948.9	ctg89567_357
bin.1081.derep	PF01281.14	ctg90131_1131
bin.1081.derep	PF01281.14	ctg89567_357
bin.1081.derep	PF01195.14	ctg89401_235
bin.1081.derep	PF00276.15	ctg90131_624
bin.1081.derep	PF00252.13	ctg90131_619
bin.1081.derep	PF00466.15	ctg89567_166
bin.1081.derep	PF03947.13	ctg90131_623
bin.1081.derep	PF00181.18	ctg90131_623
bin.1081.derep	PF00453.13	ctg88674_113
bin.1081.derep	PF00453.13	ctg90131_337
bin.1081.derep	PF00410.14	ctg90131_612
bin.1081.derep	PF01746.16	ctg90059_696
bin.1081.derep	PF08459.6	ctg89401_131
bin.1081.derep	PF00338.17	ctg90131_627
bin.1081.derep	PF00687.16	ctg89567_167
bin.1081.derep	PF01196.14	ctg90131_593
bin.1081.derep	PF01000.21	ctg90131_595
bin.1081.derep	PF01193.19	ctg90131_595
bin.1081.derep	PF00164.20	ctg89096_93
bin.1081.derep	PF00164.20	ctg90131_859
bin.1081.derep	PF03484.10	ctg90131_50
bin.1081.derep	PF03484.10	ctg90059_145
bin.1081.derep	PF00238.14	ctg90131_616
bin.1081.derep	PF00297.17	ctg90131_626
bin.1081.derep	PF00237.14	ctg90131_621
bin.1081.derep	PF02130.12	ctg89401_257
bin.1081.derep	PF00312.17	ctg88674_143
bin.1081.derep	PF00312.17	ctg90131_385
bin.1081.derep	PF01121.15	ctg90059_749
bin.1081.derep	PF02033.13	ctg88674_147
bin.1081.derep	PF02033.13	ctg90131_390
bin.1081.derep	PF03946.9	ctg89567_168
bin.1081.derep	PF00298.14	ctg89567_168
bin.1081.derep	PF00831.18	ctg90131_618
bin.1081.derep	PF00572.13	ctg90059_673
bin.1081.derep	PF01245.15	ctg90059_695
bin.1081.derep	PF00347.18	ctg90131_610&&ctg90131_611
bin.1081.derep	PF00828.14	ctg90131_606

bin.1081.derep	PF05697.8	ctg88963_12
bin.1081.derep	PF05697.8	ctg89401_203

bin.1081.derep	PF00673.16	ctg90131_614
bin.1081.derep	PF00281.14	ctg90131_614
bin.1081.derep	PF00318.15	ctg88674_163
bin.1081.derep	PF00318.15	ctg90131_423
bin.1081.derep	PF00829.16	ctg90131_241
bin.1081.derep	PF00829.16	ctg90059_4





bin.1081.derep	TIGR02432	ctg89096_19
bin.1081.derep	TIGR02432	ctg90131_1042
bin.1081.derep	TIGR00755	ctg89567_48
bin.1081.derep	TIGR00019	ctg89567_277
bin.1081.derep	TIGR00344	ctg89567_77
bin.1081.derep	TIGR03263	ctg89401_337
bin.1081.derep	TIGR03263	ctg90059_902
bin.1081.derep	TIGR03723	ctg89567_174
bin.1081.derep	TIGR00392	ctg89567_31
bin.1081.derep	TIGR03594	ctg90059_354
bin.1081.derep	TIGR00250	ctg90131_426
bin.1081.derep	TIGR00250	ctg88674_165
bin.1081.derep	TIGR02075	ctg88674_161
bin.1081.derep	TIGR01079	ctg90131_615
bin.1081.derep	TIGR00329	ctg89567_174
bin.1081.derep	TIGR00922	ctg89567_169
bin.1081.derep	TIGR00460	ctg90059_896
bin.1081.derep	TIGR00460	ctg89441_20
bin.1081.derep	TIGR00459	ctg89401_77
bin.1081.derep	TIGR00855	ctg89567_165
bin.1081.derep	TIGR00084	ctg90059_268
bin.1081.derep	TIGR00615	ctg89401_85
bin.1081.derep	PF01409.15	ctg90059_191
bin.1081.derep	PF01409.15	ctg90131_51
bin.1081.derep	PF02912.13	ctg90059_191
bin.1081.derep	PF02912.13	ctg90131_51
bin.1081.derep	PF00828.14	ctg90131_606
bin.1081.derep	PF06071.8	ctg89096_18
bin.1081.derep	PF06071.8	ctg90131_1044
bin.1081.derep	PF04997.7	ctg90131_863
bin.1081.derep	PF04997.7	ctg89096_90&&ctg89096_91
bin.1081.derep	PF00623.15	ctg90131_863
bin.1081.derep	PF00623.15	ctg89096_91
bin.1081.derep	PF04983.13	ctg90131_863
bin.1081.derep	PF04983.13	ctg89096_91
bin.1081.derep	PF05000.12	ctg90131_863
bin.1081.derep	PF05000.12	ctg89096_91
bin.1081.derep	PF04998.12	ctg89096_91
bin.1081.derep	PF04998.12	ctg90131_862
bin.1081.derep	PF01016.14	ctg90131_243
bin.1081.derep	PF01016.14	ctg90059_85
bin.1081.derep	PF01245.15	ctg90059_743
bin.1081.derep	PF02367.12	ctg90059_934
bin.1081.derep	PF02367.12	ctg89401_314
bin.1081.derep	PF08529.6	ctg88674_151
bin.1081.derep	PF08529.6	ctg90131_398
bin.1081.derep	PF13184.1	ctg88674_151
bin.1081.derep	PF13184.1	ctg90131_396
bin.1081.derep	PF01765.14	ctg88674_160
bin.1081.derep	PF01765.14	ctg90131_418
bin.1081.derep	PF00673.16	ctg90131_614
bin.1081.derep	PF00281.14	ctg90131_614
bin.1081.derep	PF00410.14	ctg90131_612
bin.1081.derep	PF03946.9	ctg89567_168
bin.1081.derep	PF00298.14	ctg89567_168
bin.1081.derep	PF00562.23	ctg89096_89
bin.1081.derep	PF00562.23	ctg90131_864
bin.1081.derep	PF04563.10	ctg89096_89
bin.1081.derep	PF04563.10	ctg90131_864
bin.1081.derep	PF04561.9	ctg89096_89
bin.1081.derep	PF04561.9	ctg90131_864
bin.1081.derep	PF04560.15	ctg89096_89
bin.1081.derep	PF04560.15	ctg90131_864
bin.1081.derep	PF04565.11	ctg89096_89
bin.1081.derep	PF04565.11	ctg90131_864
bin.1081.derep	PF10385.4	ctg89096_89
bin.1081.derep	PF10385.4	ctg90131_864
bin.1081.derep	PF05491.8	ctg90059_269
bin.1081.derep	PF00318.15	ctg88674_163
bin.1081.derep	PF00318.15	ctg90131_423
bin.1081.derep	PF00297.17	ctg90131_626
bin.1081.derep	PF00237.14	ctg90131_621
bin.1081.derep	PF03719.10	ctg90131_608
bin.1081.derep	PF00333.15	ctg90131_608
bin.1081.derep	PF12344.3	ctg90131_258
bin.1081.derep	PF12344.3	ctg90059_76
bin.1081.derep	PF11987.3	ctg88674_148
bin.1081.derep	PF11987.3	ctg90131_391
bin.1081.derep	PF00164.20	ctg89096_93
bin.1081.derep	PF00164.20	ctg90131_859
bin.1081.derep	PF00338.17	ctg90131_627
bin.1081.derep	PF00453.13	ctg88674_113
bin.1081.derep	PF00453.13	ctg90131_337
bin.1081.derep	PF00831.18	ctg90131_618
bin.1081.derep	PF02978.14	ctg90059_748
bin.1081.derep	PF00573.17	ctg90131_625
bin.1081.derep	PF01746.16	ctg90059_744
bin.1081.derep	PF03484.10	ctg90059_192
bin.1081.derep	PF03484.10	ctg90131_50
bin.1081.derep	PF01018.17	ctg89401_326
bin.1081.derep	PF01018.17	ctg90059_913&&ctg90059_914
bin.1081.derep	PF00276.15	ctg90131_624
bin.1081.derep	PF00572.13	ctg90059_721
bin.1081.derep	PF00177.16	ctg89096_94
bin.1081.derep	PF00177.16	ctg90131_858
bin.1081.derep	PF00189.15	ctg90131_620
bin.1081.derep	PF01000.21	ctg90131_595
bin.1081.derep	PF01193.19	ctg90131_595
bin.1081.derep	PF01632.14	ctg88674_114
bin.1081.derep	PF01632.14	ctg90131_338
bin.1081.derep	PF02130.12	ctg89401_257
bin.1081.derep	PF00687.16	ctg89567_167
bin.1081.derep	PF03947.13	ctg90131_623
bin.1081.derep	PF00181.18	ctg90131_623
bin.1081.derep	PF06421.7	ctg90059_497
bin.1081.derep	PF00203.16	ctg90131_622
bin.1081.derep	PF02033.13	ctg88674_147
bin.1081.derep	PF02033.13	ctg90131_390
bin.1081.derep	PF00380.14	ctg90059_720
bin.1081.derep	PF05697.8	ctg88963_12
bin.1081.derep	PF05697.8	ctg89401_203
bin.1081.derep	PF00238.14	ctg90131_616
bin.1081.derep	PF01250.12	ctg89567_354
bin.1081.derep	PF01250.12	ctg90131_1135
bin.1081.derep	PF01195.14	ctg89401_235
bin.1081.derep	PF08459.6	ctg89401_131
bin.1081.derep	PF00366.15	ctg90131_617
bin.1081.derep	PF01649.13	ctg90059_815
bin.1081.derep	PF00312.17	ctg88674_143
bin.1081.derep	PF00312.17	ctg90131_385
bin.1081.derep	PF00347.18	ctg90131_610&&ctg90131_611
bin.1081.derep	PF00466.15	ctg89567_166
bin.1081.derep	PF03948.9	ctg90131_1131
bin.1081.derep	PF03948.9	ctg89567_357
bin.1081.derep	PF01281.14	ctg90131_1131
bin.1081.derep	PF01281.14	ctg89567_357
bin.1081.derep	PF00411.14	ctg90131_597
bin.1081.derep	PF00861.17	ctg90131_609
bin.1081.derep	PF00886.14	ctg90059_747
bin.1081.derep	PF01196.14	ctg90131_593
bin.1081.derep	PF01121.15	ctg90059_797
bin.1081.derep	PF01795.14	ctg89401_214
bin.1081.derep	PF01795.14	ctg88963_1
bin.1081.derep	PF00889.14	ctg90131_422
bin.1081.derep	PF00889.14	ctg88674_162
bin.1081.derep	PF00252.13	ctg90131_619
bin.1081.derep	PF00829.16	ctg90131_241
bin.1081.derep	PF00829.16	ctg90059_87
bin.1081.derep	PF01509.13	ctg88674_145
bin.1081.derep	PF01509.13	ctg90131_387&&ctg90131_388
bin.1081.derep	PF00416.17	ctg90131_598
















bin.1081.derep	TIGR00084	ctg87283_55
bin.1081.derep	TIGR00084	ctg90059_268

bin.1081.derep	TIGR03594	ctg87283_117
bin.1081.derep	TIGR03594	ctg90059_354

bin.1081.derep	PF05491.8	ctg87283_56
bin.1081.derep	PF05491.8	ctg90059_269

bin.1081.derep	TIGR00329	ctg89567_174
bin.1081.derep	TIGR00855	ctg89567_165
bin.1081.derep	TIGR00392	ctg89567_31
bin.1081.derep	TIGR02075	ctg88674_161
bin.1081.derep	TIGR00250	ctg90131_427
bin.1081.derep	TIGR00250	ctg88674_165
bin.1081.derep	TIGR01079	ctg90131_616

bin.1081.derep	TIGR00460	ctg90059_896
bin.1081.derep	TIGR00460	ctg89441_20

bin.1081.derep	TIGR03723	ctg89567_174
bin.1081.derep	TIGR00459	ctg89401_78
bin.1081.derep	TIGR00344	ctg89567_77
bin.1081.derep	TIGR00019	ctg89567_277
bin.1081.derep	TIGR02432	ctg89096_19
bin.1081.derep	TIGR02432	ctg90131_1044
bin.1081.derep	TIGR00922	ctg89567_169
bin.1081.derep	TIGR03263	ctg89401_338
bin.1081.derep	TIGR03263	ctg90059_902
bin.1081.derep	TIGR00615	ctg89401_86
bin.1081.derep	TIGR00967	ctg90131_606
bin.1081.derep	TIGR00755	ctg89567_48
bin.1081.derep	PF01409.15	ctg90059_191
bin.1081.derep	PF01409.15	ctg90131_51
bin.1081.derep	PF02912.13	ctg90059_191
bin.1081.derep	PF02912.13	ctg90131_51
bin.1081.derep	PF02367.12	ctg90059_934
bin.1081.derep	PF02367.12	ctg89401_315
bin.1081.derep	PF02130.12	ctg89401_258
bin.1081.derep	PF00687.16	ctg89567_167
bin.1081.derep	PF00829.16	ctg90131_242
bin.1081.derep	PF00829.16	ctg90059_87
bin.1081.derep	PF01509.13	ctg88674_145
bin.1081.derep	PF01509.13	ctg90131_388&&ctg90131_389
bin.1081.derep	PF01795.14	ctg89401_215
bin.1081.derep	PF01795.14	ctg88963_1
bin.1081.derep	PF04998.12	ctg89096_91
bin.1081.derep	PF04998.12	ctg90131_864
bin.1081.derep	PF04997.7	ctg90131_865
bin.1081.derep	PF04997.7	ctg89096_90&&ctg89096_91
bin.1081.derep	PF00623.15	ctg89096_91
bin.1081.derep	PF00623.15	ctg90131_865
bin.1081.derep	PF04983.13	ctg89096_91
bin.1081.derep	PF04983.13	ctg90131_865
bin.1081.derep	PF05000.12	ctg89096_91
bin.1081.derep	PF05000.12	ctg90131_865
bin.1081.derep	PF03484.10	ctg90059_192
bin.1081.derep	PF03484.10	ctg90131_50
bin.1081.derep	PF01196.14	ctg90131_594
bin.1081.derep	PF01000.21	ctg90131_596
bin.1081.derep	PF01193.19	ctg90131_596
bin.1081.derep	PF08459.6	ctg89401_132
bin.1081.derep	PF12344.3	ctg90131_259
bin.1081.derep	PF12344.3	ctg90059_76
bin.1081.derep	PF00411.14	ctg90131_598
bin.1081.derep	PF00861.17	ctg90131_610
bin.1081.derep	PF11987.3	ctg90131_392
bin.1081.derep	PF11987.3	ctg88674_148
bin.1081.derep	PF01632.14	ctg88674_114
bin.1081.derep	PF01632.14	ctg90131_339
bin.1081.derep	PF03719.10	ctg90131_609
bin.1081.derep	PF00333.15	ctg90131_609
bin.1081.derep	PF00203.16	ctg90131_623
bin.1081.derep	PF00238.14	ctg90131_617
bin.1081.derep	PF13184.1	ctg90131_397
bin.1081.derep	PF13184.1	ctg88674_151
bin.1081.derep	PF00312.17	ctg88674_143
bin.1081.derep	PF00312.17	ctg90131_386
bin.1081.derep	PF01018.17	ctg89401_327
bin.1081.derep	PF01018.17	ctg90059_913&&ctg90059_914
bin.1081.derep	PF06421.7	ctg90059_497
bin.1081.derep	PF00416.17	ctg90131_599
bin.1081.derep	PF03946.9	ctg89567_168
bin.1081.derep	PF00298.14	ctg89567_168
bin.1081.derep	PF03948.9	ctg89567_357
bin.1081.derep	PF03948.9	ctg90131_1133
bin.1081.derep	PF01281.14	ctg89567_357
bin.1081.derep	PF01281.14	ctg90131_1133
bin.1081.derep	PF08529.6	ctg90131_399
bin.1081.derep	PF08529.6	ctg88674_151
bin.1081.derep	PF00572.13	ctg90059_721
bin.1081.derep	PF00562.23	ctg89096_89
bin.1081.derep	PF00562.23	ctg90131_866
bin.1081.derep	PF04563.10	ctg89096_89
bin.1081.derep	PF04563.10	ctg90131_866
bin.1081.derep	PF04561.9	ctg89096_89
bin.1081.derep	PF04561.9	ctg90131_866
bin.1081.derep	PF04560.15	ctg89096_89
bin.1081.derep	PF04560.15	ctg90131_866
bin.1081.derep	PF04565.11	ctg89096_89
bin.1081.derep	PF04565.11	ctg90131_866
bin.1081.derep	PF10385.4	ctg89096_89
bin.1081.derep	PF10385.4	ctg90131_866
bin.1081.derep	PF00380.14	ctg90059_720
bin.1081.derep	PF00347.18	ctg90131_611&&ctg90131_612
bin.1081.derep	PF00573.17	ctg90131_626
bin.1081.derep	PF01016.14	ctg90131_244
bin.1081.derep	PF01016.14	ctg90059_85
bin.1081.derep	PF03947.13	ctg90131_624
bin.1081.derep	PF00181.18	ctg90131_624
bin.1081.derep	PF00177.16	ctg89096_94
bin.1081.derep	PF00177.16	ctg90131_860
bin.1081.derep	PF00673.16	ctg90131_615
bin.1081.derep	PF00281.14	ctg90131_615
bin.1081.derep	PF00828.14	ctg90131_607
bin.1081.derep	PF00410.14	ctg90131_613
bin.1081.derep	PF01649.13	ctg90059_815
bin.1081.derep	PF00831.18	ctg90131_619
bin.1081.derep	PF01121.15	ctg90059_797
bin.1081.derep	PF00276.15	ctg90131_625
bin.1081.derep	PF01195.14	ctg89401_236
bin.1081.derep	PF00318.15	ctg88674_163
bin.1081.derep	PF00318.15	ctg90131_424
bin.1081.derep	PF00466.15	ctg89567_166
bin.1081.derep	PF00189.15	ctg90131_621
bin.1081.derep	PF02033.13	ctg88674_147
bin.1081.derep	PF02033.13	ctg90131_391
bin.1081.derep	PF01746.16	ctg90059_744
bin.1081.derep	PF00453.13	ctg88674_113
bin.1081.derep	PF00453.13	ctg90131_338
bin.1081.derep	PF00252.13	ctg90131_620
bin.1081.derep	PF00886.14	ctg90059_747
bin.1081.derep	PF00366.15	ctg90131_618
bin.1081.derep	PF06071.8	ctg89096_18
bin.1081.derep	PF06071.8	ctg90131_1046
bin.1081.derep	PF01250.12	ctg89567_354
bin.1081.derep	PF01250.12	ctg90131_1137
bin.1081.derep	PF01245.15	ctg90059_743
bin.1081.derep	PF00889.14	ctg90131_423
bin.1081.derep	PF00889.14	ctg88674_162
bin.1081.derep	PF02978.14	ctg90059_748
bin.1081.derep	PF00237.14	ctg90131_622
bin.1081.derep	PF00297.17	ctg90131_627
bin.1081.derep	PF00338.17	ctg90131_628
bin.1081.derep	PF05697.8	ctg88963_12
bin.1081.derep	PF05697.8	ctg89401_204
bin.1081.derep	PF01765.14	ctg88674_160
bin.1081.derep	PF01765.14	ctg90131_419
bin.1081.derep	PF00164.20	ctg89096_93
bin.1081.derep	PF00164.20	ctg90131_861




















bin.1081.derep	TIGR00615	ctg88963_123
bin.1081.derep	TIGR00615	ctg89401_87 --

bin.1081.derep	TIGR01079	ctg89096_282 --
bin.1081.derep	TIGR01079	ctg90131_612

bin.1081.derep	TIGR03594	ctg87283_117
bin.1081.derep	TIGR03594	ctg90059_397 --

bin.1081.derep	TIGR00250	ctg90131_424
bin.1081.derep	TIGR00250	ctg88674_165

bin.1081.derep	TIGR00084	ctg87283_55
bin.1081.derep	TIGR00084	ctg90059_311

bin.1081.derep	TIGR00967	ctg89096_291
bin.1081.derep	TIGR00967	ctg90131_602

bin.1081.derep	TIGR00459	ctg88963_132
bin.1081.derep	TIGR00459	ctg89401_79 --

bin.1081.derep	TIGR03263	ctg89401_339 --
bin.1081.derep	TIGR03263	ctg90059_945 --

bin.1081.derep	TIGR02432	ctg89096_72
bin.1081.derep	TIGR02432	ctg90131_1038

bin.1081.derep	PF01018.17	ctg89401_328 --
bin.1081.derep	PF01018.17	ctg90059_956&&ctg90059_957

bin.1081.derep	PF06071.8	ctg89096_71
bin.1081.derep	PF06071.8	ctg90131_1040

bin.1081.derep	PF11987.3	ctg88674_148
bin.1081.derep	PF11987.3	ctg90131_389

bin.1081.derep	PF04998.12	ctg89096_144
bin.1081.derep	PF04998.12	ctg90131_858

bin.1081.derep	PF04997.7	ctg90131_859
bin.1081.derep	PF04997.7	ctg89096_143&&ctg89096_144

bin.1081.derep	PF00623.15	ctg89096_144
bin.1081.derep	PF00623.15	ctg90131_859

bin.1081.derep	PF04983.13	ctg89096_144
bin.1081.derep	PF04983.13	ctg90131_859

bin.1081.derep	PF05000.12	ctg89096_144
bin.1081.derep	PF05000.12	ctg90131_859

bin.1081.derep	PF00276.15	ctg89096_273
bin.1081.derep	PF00276.15	ctg90131_621

bin.1081.derep	PF00366.15	ctg89096_280
bin.1081.derep	PF00366.15	ctg90131_614

bin.1081.derep	PF05491.8	ctg87283_56
bin.1081.derep	PF05491.8	ctg90059_312 --

bin.1081.derep	PF00177.16	ctg89096_147
bin.1081.derep	PF00177.16	ctg90131_854

bin.1081.derep	PF00861.17	ctg90131_606
bin.1081.derep	PF00861.17	ctg89096_287

bin.1081.derep	PF01409.15	ctg90059_234
bin.1081.derep	PF01409.15	ctg90131_51

bin.1081.derep	PF02912.13	ctg90059_234
bin.1081.derep	PF02912.13	ctg90131_51

bin.1081.derep	PF03947.13	ctg90131_620
bin.1081.derep	PF03947.13	ctg89096_274

bin.1081.derep	PF00181.18	ctg90131_620
bin.1081.derep	PF00181.18	ctg89096_274

bin.1081.derep	PF00453.13	ctg88674_113
bin.1081.derep	PF00453.13	ctg90131_335

bin.1081.derep	PF00562.23	ctg90131_860
bin.1081.derep	PF00562.23	ctg89096_142

bin.1081.derep	PF04563.10	ctg90131_860
bin.1081.derep	PF04563.10	ctg89096_142

bin.1081.derep	PF04561.9	ctg90131_860
bin.1081.derep	PF04561.9	ctg89096_142

bin.1081.derep	PF04560.15	ctg90131_860
bin.1081.derep	PF04560.15	ctg89096_142

bin.1081.derep	PF04565.11	ctg90131_860
bin.1081.derep	PF04565.11	ctg89096_142

bin.1081.derep	PF10385.4	ctg90131_860
bin.1081.derep	PF10385.4	ctg89096_142

bin.1081.derep	PF03719.10	ctg89096_288
bin.1081.derep	PF03719.10	ctg90131_605

bin.1081.derep	PF00333.15	ctg89096_288
bin.1081.derep	PF00333.15	ctg90131_605

bin.1081.derep	PF00673.16	ctg89096_283
bin.1081.derep	PF00673.16	ctg90131_611

bin.1081.derep	PF00281.14	ctg89096_283
bin.1081.derep	PF00281.14	ctg90131_611

bin.1081.derep	PF00831.18	ctg89096_279
bin.1081.derep	PF00831.18	ctg90131_615

bin.1081.derep	PF00829.16	ctg90131_240
bin.1081.derep	PF00829.16	ctg90059_132

bin.1081.derep	PF00410.14	ctg89096_285
bin.1081.derep	PF00410.14	ctg90131_609

bin.1081.derep	PF00828.14	ctg89096_290
bin.1081.derep	PF00828.14	ctg90131_603

bin.1081.derep	PF00203.16	ctg89096_275
bin.1081.derep	PF00203.16	ctg90131_619
bin.1081.derep	PF00889.14	ctg90131_420
bin.1081.derep	PF00889.14	ctg88674_162
bin.1081.derep	PF00252.13	ctg89096_278
bin.1081.derep	PF00252.13	ctg90131_616

bin.1081.derep	PF02367.12	ctg90059_976
bin.1081.derep	PF02367.12	ctg89401_316 -

bin.1081.derep	PF03484.10	ctg90059_235
bin.1081.derep	PF03484.10	ctg90131_50
bin.1081.derep	PF00347.18	ctg89096_286
bin.1081.derep	PF00347.18	ctg90131_607&&ctg90131_608
bin.1081.derep	PF02978.14	ctg90059_790
bin.1081.derep	PF00380.14	ctg90059_762
bin.1081.derep	PF00237.14	ctg89096_276
bin.1081.derep	PF00237.14	ctg90131_618
bin.1081.derep	PF00886.14	ctg90059_789
bin.1081.derep	PF02033.13	ctg88674_147
bin.1081.derep	PF02033.13	ctg90131_388
bin.1081.derep	PF00416.17	ctg89096_297
bin.1081.derep	PF00416.17	ctg90131_595
bin.1081.derep	PF01746.16	ctg90059_786
bin.1081.derep	PF03948.9	ctg89567_415
bin.1081.derep	PF03948.9	ctg90131_1127
bin.1081.derep	PF01281.14	ctg89567_415
bin.1081.derep	PF01281.14	ctg90131_1127

bin.1081.derep	PF01795.14	ctg89401_216 -
bin.1081.derep	PF01795.14	ctg88963_1

bin.1081.derep	PF01509.13	ctg88674_145
bin.1081.derep	PF01509.13	ctg90131_385&&ctg90131_386
bin.1081.derep	PF00318.15	ctg88674_163
bin.1081.derep	PF00318.15	ctg90131_421

bin.1081.derep	PF05697.8	ctg88963_12
bin.1081.derep	PF05697.8	ctg89401_205 -

bin.1081.derep	PF00338.17	ctg89096_270
bin.1081.derep	PF00338.17	ctg90131_624
bin.1081.derep	PF08529.6	ctg88674_151
bin.1081.derep	PF08529.6	ctg90131_396
bin.1081.derep	PF13184.1	ctg88674_151
bin.1081.derep	PF13184.1	ctg90131_394

bin.1081.derep	PF08459.6	ctg88963_82
bin.1081.derep	PF08459.6	ctg89401_133 -

bin.1081.derep	PF01016.14	ctg90131_242
bin.1081.derep	PF01016.14	ctg90059_130
bin.1081.derep	PF01250.12	ctg89567_412
bin.1081.derep	PF01250.12	ctg90131_1131
bin.1081.derep	PF00189.15	ctg89096_277
bin.1081.derep	PF00189.15	ctg90131_617
bin.1081.derep	PF00297.17	ctg89096_271
bin.1081.derep	PF00297.17	ctg90131_623
bin.1081.derep	PF00238.14	ctg89096_281
bin.1081.derep	PF00238.14	ctg90131_613
bin.1081.derep	PF00573.17	ctg89096_272
bin.1081.derep	PF00573.17	ctg90131_622
bin.1081.derep	PF00312.17	ctg88674_143
bin.1081.derep	PF00312.17	ctg90131_383
bin.1081.derep	PF00164.20	ctg89096_146
bin.1081.derep	PF00164.20	ctg90131_855
bin.1081.derep	PF01765.14	ctg88674_160
bin.1081.derep	PF01765.14	ctg90131_416
bin.1081.derep	PF00572.13	ctg90059_763
bin.1081.derep	PF01196.14	ctg90131_590
bin.1081.derep	PF00411.14	ctg90131_594
bin.1081.derep	PF00411.14	ctg89096_298
bin.1081.derep	PF01632.14	ctg88674_114
bin.1081.derep	PF01632.14	ctg90131_336
*/
#endif 


