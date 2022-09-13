

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
		

		fs::remove(outputMappingFilename);
	}

	//zFile _queryContigFile;
	unordered_set<u_int32_t> _duplicatedContigIndex;
	unordered_map<u_int32_t, vector<DbgEdge>> _duplicationBounds;

	void mapContigs(const string& mapFilename){


		//-N 2 -p 0 --secondary=yes -I 2GB -x map-hifi -c -x asm20
		string command = "minimap2 -H -DP -c -I 100M -t " + to_string(_nbCores) + " " + _inputFilename_contigs + " " + _inputFilename_contigs + " > " + mapFilename;
		Utils::executeCommand(command, _tmpDir);

	}


	unordered_set<DbgEdge, hash_pair> _performedPairs;

	void detectDuplicatedContigs(const string& mapFilename){

		ifstream mappingFile(mapFilename);
        vector<string>* fields = new vector<string>();
        vector<string>* fields_optional = new vector<string>();

		string line;
		while (getline(mappingFile, line)) {


            GfaParser::tokenize(line, fields, '\t');

			const string& readName = (*fields)[0];
			const string& contigName = (*fields)[5];
			//if(readName == contigName) continue;



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
				//cout << line << endl;
			}
			if((readName == "ctg89567" || contigName == "ctg89567")){
				//cout << line << endl;
			}


			//cout << line << endl;

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

			//cout << nbMatches / alignLength << " " << alignLength << endl;
			if(nbMatches / alignLength < 0.8) continue;
			if(alignLength < 1000) continue;

			u_int64_t maxHang = 100;
			u_int64_t hangLeft = targetStart;
			u_int64_t hangRight = targetLength - targetEnd;

			

			u_int64_t targetIndex = Utils::contigName_to_contigIndex(contigName);
			u_int64_t queryIndex = Utils::contigName_to_contigIndex(readName);
			DbgEdge edge = {targetIndex, queryIndex};
			edge = edge.normalize();

			bool isOverlap = false;

			if(hangLeft < maxHang){
				if(_performedPairs.find(edge) == _performedPairs.end()){
					//cout << "Left overlap: " << line << endl;
					//cout << line << endl;
					_duplicationBounds[targetIndex].push_back({targetStart, targetEnd});
					cout << "overlap left: " << contigName << " " << targetStart << " " << targetEnd << endl;
					isOverlap = true;
				}
			}

			if(hangRight < maxHang){
				if(_performedPairs.find(edge) == _performedPairs.end()){
					//cout << "Right overlap: " << line << endl;
					//cout << line << endl;
					_duplicationBounds[targetIndex].push_back({targetStart, targetEnd});
					cout << "overlap right: " << contigName << " " << targetStart << " " << targetEnd << endl;
					isOverlap = true;
				}
			}

			if(isOverlap){
				_performedPairs.insert(edge);
			}


			if(nbMatches / alignLength < 0.95) continue;
			if(alignLength < 10000) continue;
			if(targetLength > queryLength) continue;


			//_duplicationBounds[targetIndex].push_back({targetStart, targetEnd});
			//cout << "Add internal " << targetStart << " " << targetEnd << endl;
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

		/*
		if(read._header == "ctg88963"){
			string header = ">" + read._header + '\n';
			gzwrite(_file1, (const char*)&header[0], header.size());
			string contigSequence = read._seq + '\n';
			gzwrite(_file1, (const char*)&contigSequence[0], contigSequence.size());


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

		*/

		u_int64_t contigIndex = Utils::contigName_to_contigIndex(read._header);

		string seq = read._seq;

		if(_duplicationBounds.find(contigIndex) != _duplicationBounds.end()){

			//if(read._header == "ctg89567") 
			cout << "-----" << endl;
			vector<bool> isDuplicated(seq.size(), false);

			for(const DbgEdge& duplicatedBound : _duplicationBounds[contigIndex]){
				//if(read._header == "ctg89567") 
				cout << duplicatedBound._from << " " << duplicatedBound._to << endl;
				for(size_t i=duplicatedBound._from; i<duplicatedBound._to; i++){
					isDuplicated[i] = true;
				}
			}

			bool isDuplicatedArea = true;
			long startPos = 0;
			size_t subSeqIndex = 0;
			long endPos = 0;

			for(long i=0; i<seq.size(); i++){
				
				endPos = i;

				if(isDuplicated[i]){
					if(!isDuplicatedArea){
						long length = endPos-startPos-1;
						if(length >= 500){
							string header = read._header;


							char lastChar = header[header.size()-1];
							header.pop_back();

							header += "_" + to_string(subSeqIndex) + lastChar;
							
							string subSeq = seq.substr(startPos, length);

							header = ">" + header + '\n';
							gzwrite(_outputContigFile, (const char*)&header[0], header.size());
							string contigSequence = subSeq + '\n';
							gzwrite(_outputContigFile, (const char*)&contigSequence[0], contigSequence.size());


							//if(read._header == "ctg89567")
							cout << "Dump area: " << startPos << " " << startPos+length << endl;

							subSeqIndex += 1;
						}


					}
					isDuplicatedArea = true;
				}
				else{
					if(isDuplicatedArea){
						startPos = i;
					}
					isDuplicatedArea = false;
				}
			}

			if(!isDuplicatedArea){
				long length = endPos-startPos-1;
				if(length > 500){
					string header = read._header;
					char lastChar = header[header.size()-1];
					header.pop_back();

					header += "_" + to_string(subSeqIndex) + lastChar;
					
					string subSeq = seq.substr(startPos, length);

					header = ">" + header + '\n';
					gzwrite(_outputContigFile, (const char*)&header[0], header.size());
					string contigSequence = subSeq + '\n';
					gzwrite(_outputContigFile, (const char*)&contigSequence[0], contigSequence.size());
					
					cout << "Dump area: " << startPos << " " << startPos+length << endl;
				}


			}

		}
		else{
			string header = ">" + read._header + '\n';
			gzwrite(_outputContigFile, (const char*)&header[0], header.size());
			string contigSequence = seq + '\n';
			gzwrite(_outputContigFile, (const char*)&contigSequence[0], contigSequence.size());
		}



	}


};	







/*

bin.1081.derep	TIGR03594	ctg87283_117
bin.1081.derep	TIGR03594	ctg90059_400
bin.1081.derep	TIGR00855	ctg89567_224
bin.1081.derep	TIGR01079	ctg90131_612
bin.1081.derep	TIGR02075	ctg88674_149
bin.1081.derep	TIGR00250	ctg90131_424
bin.1081.derep	TIGR00250	ctg88674_153
bin.1081.derep	TIGR00967	ctg90131_602
bin.1081.derep	TIGR00084	ctg87283_55
bin.1081.derep	TIGR00084	ctg90059_314
bin.1081.derep	TIGR00922	ctg89567_228
bin.1081.derep	TIGR03723	ctg89567_233
bin.1081.derep	TIGR00392	ctg89567_90
bin.1081.derep	TIGR00329	ctg89567_233
bin.1081.derep	TIGR00019	ctg89567_336
bin.1081.derep	TIGR00019	ctg86936_63
bin.1081.derep	TIGR03263	ctg89401_340
bin.1081.derep	TIGR00344	ctg89567_136
bin.1081.derep	TIGR00810	ctg85082_33
bin.1081.derep	TIGR00810	ctg89567_52
bin.1081.derep	TIGR00460	ctg90059_942
bin.1081.derep	TIGR00460	ctg89441_25
bin.1081.derep	TIGR00459	ctg89401_80
bin.1081.derep	TIGR02432	ctg89096_72
bin.1081.derep	TIGR02432	ctg90131_1040
bin.1081.derep	TIGR00755	ctg89567_107
bin.1081.derep	TIGR00615	ctg89401_88
bin.1081.derep	PF00318.15	ctg88674_151
bin.1081.derep	PF00318.15	ctg90131_421
bin.1081.derep	PF01746.16	ctg90059_790
bin.1081.derep	PF12344.3	ctg90059_122
bin.1081.derep	PF12344.3	ctg90131_255
bin.1081.derep	PF00572.13	ctg90059_767
bin.1081.derep	PF04998.12	ctg89096_144
bin.1081.derep	PF04998.12	ctg90131_860
bin.1081.derep	PF04997.7	ctg90131_861
bin.1081.derep	PF04997.7	ctg89096_143&&ctg89096_144
bin.1081.derep	PF00623.15	ctg89096_144
bin.1081.derep	PF00623.15	ctg90131_861
bin.1081.derep	PF04983.13	ctg89096_144
bin.1081.derep	PF04983.13	ctg90131_861
bin.1081.derep	PF05000.12	ctg89096_144
bin.1081.derep	PF05000.12	ctg90131_861
bin.1081.derep	PF11987.3	ctg88674_136
bin.1081.derep	PF11987.3	ctg90131_389
bin.1081.derep	PF00237.14	ctg90131_618
bin.1081.derep	PF03947.13	ctg90131_620
bin.1081.derep	PF00181.18	ctg90131_620
bin.1081.derep	PF01409.15	ctg90059_236
bin.1081.derep	PF01409.15	ctg90131_51
bin.1081.derep	PF02912.13	ctg90059_236
bin.1081.derep	PF02912.13	ctg90131_51
bin.1081.derep	PF02978.14	ctg90059_794
bin.1081.derep	PF03946.9	ctg89567_227
bin.1081.derep	PF00298.14	ctg89567_227
bin.1081.derep	PF01668.13	ctg89567_50
bin.1081.derep	PF01668.13	ctg85082_35
bin.1081.derep	PF05697.8	ctg89401_206
bin.1081.derep	PF06071.8	ctg89096_71
bin.1081.derep	PF06071.8	ctg90131_1042
bin.1081.derep	PF00252.13	ctg90131_616
bin.1081.derep	PF01632.14	ctg88674_102
bin.1081.derep	PF01632.14	ctg90131_335
bin.1081.derep	PF00886.14	ctg90059_793
bin.1081.derep	PF01018.17	ctg89401_329
bin.1081.derep	PF05491.8	ctg87283_56
bin.1081.derep	PF05491.8	ctg90059_315
bin.1081.derep	PF00276.15	ctg90131_621
bin.1081.derep	PF08529.6	ctg88674_139
bin.1081.derep	PF08529.6	ctg90131_396
bin.1081.derep	PF13184.1	ctg88674_139
bin.1081.derep	PF13184.1	ctg90131_394
bin.1081.derep	PF00411.14	ctg90131_594
bin.1081.derep	PF00861.17	ctg90131_606
bin.1081.derep	PF00562.23	ctg89096_142
bin.1081.derep	PF00562.23	ctg90131_862
bin.1081.derep	PF04563.10	ctg89096_142
bin.1081.derep	PF04563.10	ctg90131_862
bin.1081.derep	PF04561.9	ctg89096_142
bin.1081.derep	PF04561.9	ctg90131_862
bin.1081.derep	PF04560.15	ctg89096_142
bin.1081.derep	PF04560.15	ctg90131_862
bin.1081.derep	PF04565.11	ctg89096_142
bin.1081.derep	PF04565.11	ctg90131_862
bin.1081.derep	PF10385.4	ctg89096_142
bin.1081.derep	PF10385.4	ctg90131_862
bin.1081.derep	PF00673.16	ctg90131_611
bin.1081.derep	PF00281.14	ctg90131_611
bin.1081.derep	PF00889.14	ctg90131_420
bin.1081.derep	PF00889.14	ctg88674_150
bin.1081.derep	PF00312.17	ctg88674_131
bin.1081.derep	PF00312.17	ctg90131_383
bin.1081.derep	PF01245.15	ctg90059_789
bin.1081.derep	PF01016.14	ctg90131_241
bin.1081.derep	PF01016.14	ctg90059_131
bin.1081.derep	PF03484.10	ctg90059_237
bin.1081.derep	PF03484.10	ctg90131_50
bin.1081.derep	PF00347.18	ctg90131_607&&ctg90131_608
bin.1081.derep	PF00829.16	ctg90131_239
bin.1081.derep	PF00829.16	ctg90059_133
bin.1081.derep	PF01196.14	ctg90131_590
bin.1081.derep	PF00416.17	ctg90131_595
bin.1081.derep	PF06421.7	ctg90059_542
bin.1081.derep	PF01765.14	ctg88674_148
bin.1081.derep	PF01765.14	ctg90131_416
bin.1081.derep	PF00238.14	ctg90131_613
bin.1081.derep	PF01649.13	ctg90059_861
bin.1081.derep	PF00189.15	ctg90131_617
bin.1081.derep	PF02130.12	ctg89401_260
bin.1081.derep	PF01121.15	ctg90059_843
bin.1081.derep	PF00164.20	ctg89096_146
bin.1081.derep	PF00164.20	ctg90131_857
bin.1081.derep	PF00366.15	ctg90131_614
bin.1081.derep	PF00687.16	ctg89567_226
bin.1081.derep	PF00410.14	ctg90131_609
bin.1081.derep	PF01509.13	ctg88674_133
bin.1081.derep	PF01509.13	ctg90131_385&&ctg90131_386
bin.1081.derep	PF00162.14	ctg89567_55
bin.1081.derep	PF00162.14	ctg85082_30
bin.1081.derep	PF00828.14	ctg90131_603
bin.1081.derep	PF00453.13	ctg88674_101
bin.1081.derep	PF00453.13	ctg90131_334
bin.1081.derep	PF01000.21	ctg90131_592
bin.1081.derep	PF01193.19	ctg90131_592
bin.1081.derep	PF00573.17	ctg90131_622
bin.1081.derep	PF00338.17	ctg90131_624
bin.1081.derep	PF03719.10	ctg90131_605
bin.1081.derep	PF00333.15	ctg90131_605
bin.1081.derep	PF02033.13	ctg88674_135
bin.1081.derep	PF02033.13	ctg90131_388
bin.1081.derep	PF08459.6	ctg89401_134
bin.1081.derep	PF01250.12	ctg90131_1133
bin.1081.derep	PF01795.14	ctg89401_217
bin.1081.derep	PF03948.9	ctg90131_1129
bin.1081.derep	PF01281.14	ctg90131_1129
bin.1081.derep	PF00177.16	ctg89096_147
bin.1081.derep	PF00177.16	ctg90131_856
bin.1081.derep	PF02367.12	ctg89401_317
bin.1081.derep	PF00203.16	ctg90131_619
bin.1081.derep	PF00380.14	ctg90059_766
bin.1081.derep	PF00831.18	ctg90131_615
bin.1081.derep	PF01195.14	ctg89401_238
bin.1081.derep	PF00466.15	ctg89567_225
bin.1081.derep	PF00297.17	ctg90131_623

*/
#endif 


