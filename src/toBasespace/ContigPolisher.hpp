

#ifndef MDBG_METAG_CONTIGPOLISHER
#define MDBG_METAG_CONTIGPOLISHER

#include "../Commons.hpp"
//#include "../src/utils/abPOA/include/abpoa.h"
#include "ContigDerep.hpp"

class ContigPolisher{
    
public:


	//constexpr static const int64_t _minContigLength = 50;
	//constexpr static const float _minContigCoverage = 1.5;
	constexpr static const int64_t _maximalMappingOffset = 300;


	string _inputFilenameRead;
	string _outputContigFilename;
	int _nbCores;
	size_t _windowLength;
	size_t _windowPositionOffset;
	size_t _windowLengthVariance;
	size_t _maxWindowCopies;
	//string _mapperOutputExeFilename;
	bool _useQual;
	string _outputDir;
	string _tmpDir;
	string _readPartitionDir;
	int _minContigLength;
	int _minContigCoverage;
	
	string _outputFilename_mapping;
	string _outputMappingFilename_contigsVsUsedReads;
	//int _minimapBatchSize;
	double _qualityThreshold;
	bool _useMetamdbgHeaderStyle;
	int _nbPartitions;

	
	//BinaryReadMap& _mReads;
	string _minimap2Preset_map;
	bool _useHpc;
	u_int64_t _checksumWrittenReads;
	u_int64_t _checksum;
	u_int64_t _checksum_total;
	u_int64_t _checksum_inputContigs;

	struct Window{
		DnaBitset2 _sequence;
		string _quality;
		u_int32_t _posStart;
		u_int32_t _posEnd;
		float _score;

		u_int64_t hash() const{
			
			u_int64_t h = 0;

			string seq = _sequence.to_string();
			//string seq = string(dnaStr);
			//free(dnaStr);

			if(_quality.size() > 0){
				for(size_t i=0; i<seq.size(); i++){
					h += seq[i] * _quality[i];
				}
			}
			else{
				for(size_t i=0; i<seq.size(); i++){
					h += seq[i];
				}
			}

			return h;
		}
	};

	// AaCcGgTtNn ==> 0,1,2,3,4
	unsigned char nt4_table[256] = {
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
	};

	const char nt256_table[256] = {
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

	//u_int64_t& _outputContigIndex;

	ContigPolisher(const string tmpDir, const string readPartitionDir, const string outputContigFilename, const int sequencingTechnology, const bool useQual, const int nbCores, const int nbPartitions, const string minimap2Preset_map, const int minContigLength, const int minContigCoverage)  {
		_tmpDir = tmpDir;
		_readPartitionDir = readPartitionDir;
		//_minimapBatchSize = minimapBatchSize;
		_useQual = useQual;
		_nbCores = nbCores;
		_nbPartitions = nbPartitions;
		_minimap2Preset_map = minimap2Preset_map;

		_inputFilenameRead = _readPartitionDir + "/input.txt";
		_outputContigFilename = outputContigFilename; //_tmpDir + "/contigs.fasta.gz";
		_useMetamdbgHeaderStyle = true;
		_windowLength = 500;
		_windowLengthVariance = _windowLength*0.02;
		_maxWindowCopies = 100;
		_qualityThreshold = 10.0;

		_minContigLength = minContigLength;
		_minContigCoverage = minContigCoverage;

		//_maxWindowCopies = 100;//
		/*
		if(sequencingTechnology == DataType::HiFi){
			cout << "Data type is HiFi" << endl;
			_minimap2Preset_map = "map-hifi";
			_minimap2Preset_ava = "ava-pb";
			_useHpc = true;
		}
		else if(sequencingTechnology == DataType::Nanopore){
			cout << "Data type is ONT" << endl;
			_minimap2Preset_map = "map-ont";
			_minimap2Preset_ava = "ava-ont";
			_useHpc = false;
		}
		*/
	}

	/*
    void execute (){



		auto startTotal = high_resolution_clock::now();
		
		Logger::get().debug() << "";
		Logger::get().debug() << "Polishing contigs (pass 1)";
		_outputContigIndex = 0;
		_checksum_total = 0;
		_checksum_inputContigs = 0;
		
		
		for(size_t i=0; i<_nbPartitions; i++){

			//if(i != 2){
			//	cout << "skip: " << i << endl;
			//	continue;
			//}

			_checksum = 0;

			//string inputContigFilename = _readPartitionDir + "/" + to_string(i) + "_contigs_norepeats.fasta.gz";
			string inputContigFilename = _readPartitionDir + "/" + to_string(i) + "_contigs_uncorrected.fasta.gz";
			string outputContigFilename = _readPartitionDir + "/" + to_string(i) + "_contigs_polished.fasta.gz";
			gzFile outputContigFile = gzopen(outputContigFilename.c_str(),"wb");

			polishPartition(i, inputContigFilename, outputContigFile, 0);

			gzclose(outputContigFile);

			cout << "ContigPolisher checksum " << i << ": " << _checksum << endl;
			cout << "Peak memory: " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << endl;
		}
	
		cout << "ContigPolisher checksum total: " << _checksum_total << endl;
		Logger::get().debug() << "Done " << " (" << duration_cast<seconds>(high_resolution_clock::now() - startTotal).count() << "s) " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << "GB";
		//

		//exit(1);
		//Logger::get().debug() << "";
		//Logger::get().debug() << "Polishing contigs (pass 2)";
		//const string& outputContigFilename_polished = _readPartitionDir + "/contigs_polished.fasta.gz";
		//pileupContigs(_nbPartitions);
		//cout << "done" << endl;
		
		
		Logger::get().debug() << "";
		Logger::get().debug() << "Polishing contigs (pass 2)";
		_outputContigIndex = 0;
		_checksum_total = 0;
		startTotal = high_resolution_clock::now();

		gzFile outputContigFile_polished = gzopen(_outputContigFilename.c_str(),"wb");
		_file_contigHeaders = ofstream(_tmpDir + "/contigHeaders.txt");

		for(size_t i=0; i<_nbPartitions; i++){

				
			_checksum = 0;
			//if(i != 30) continue;

			string inputContigFilename = _readPartitionDir + "/" + to_string(i) + "_contigs_polished.fasta.gz";

			polishPartition(i, inputContigFilename, outputContigFile_polished, 1);
			
			cout << "ContigPolisher checksum " << i << ": " << _checksum << endl;
			cout << "Peak memory: " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << endl;


		}

		_file_contigHeaders.close();
		gzclose(outputContigFile_polished);
		
		cout << "ContigPolisher checksum total: " << _checksum_total << endl;
		Logger::get().debug() << "Done " << " (" << duration_cast<seconds>(high_resolution_clock::now() - startTotal).count() << "s) " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << "GB";
		
		//exit(1);

		
		//Logger::get().debug() << "";
		//Logger::get().debug() << "Polished contigs: " << _outputContigFilename;
		//Logger::get().debug() << "done";
		//closeLogFile();
	}
	*/


    void execute2 (const u_int32_t partition, BGZF* outputContigFile_polished, ofstream& file_contigHeaders){

		string inputContigFilename = _readPartitionDir + "/contigs_uncorrected.fasta.gz";
		//breakIntraRepeats(partition, inputContigFilename);
		
		//return;
		//string inputContigFilename = _readPartitionDir + "/" + to_string(partition) + "_contigs_uncorrected.fasta.gz";
		//string outputContigFilename = _readPartitionDir + "/" + to_string(partition) + "_contigsPolished_pass0.fasta.gz";
		string outputContigFilename = _readPartitionDir + "/contigsPolished_pass0.fasta.gz";
		BGZF* outputContigFile = bgzf_open(outputContigFilename.c_str(),"w1");

		polishPartition(partition, inputContigFilename, outputContigFile, file_contigHeaders, 0);

		bgzf_close(outputContigFile);
		


		inputContigFilename = outputContigFilename; //_readPartitionDir + "/" + to_string(partition) + "_contigs_uncorrected.fasta.gz";
		
		//outputContigFilename = _readPartitionDir + "/" + to_string(partition) + "_contigsPolished_pass2.fasta.gz";
		//outputContigFile = gzopen(outputContigFilename.c_str(),"wb");

		polishPartition(partition, inputContigFilename, outputContigFile_polished, file_contigHeaders, 1);

		//gzclose(outputContigFile);

		//fs::remove(outputContigFilename);


	}


	void polishPartition(u_int32_t partition, const string& contigFilename, BGZF* outputContigFile, ofstream& file_contigHeaders, const int& pass){
		_currentPartition = partition;
		//if(_contigSequences.size() == 0) return;

		Logger::get().debug() << "";
		Logger::get().debug() << "\tPolishing contigs pass " << pass;
		//if(_debug_contigName != "") cout << "\tPolishing partition: " << _currentPartition << "/" << _nbPartitions << endl;
		auto start = high_resolution_clock::now();

		//if(_partitionNbReads[partition] == 0) return;

		clearPass();

		//if(_debug_contigName != ""){
		//	debugFindContigs(contigFilename);
		//	if(!_debug_foundContig) return;

		//	if(!fs::exists(contigFilename)) return;
		//}

		//string contigFilename = _readPartitionDir + "/part_" + to_string(partition) + "_contigs.gz";
		//string readFilename = _readPartitionDir + "/part_" + to_string(partition) + ".gz";


		string readFilename = _readPartitionDir + "/" + to_string(partition) + "_reads";
		//string contigFilename = _readPartitionDir + "/" + to_string(partition) + "_contigsRepeats.gz";
		
		//string outputContigFilename = _readPartitionDir + "/" + to_string(partition) + "_contigs.gz";

		
		indexContigName(contigFilename);

		Logger::get().debug() << "\tMap reads to curated contigs";

		const string& alignFilename = _readPartitionDir + "/align.paf.gz";
		string command = "minimap2 -I 100G -v 0 -m 500 -t " + to_string(_nbCores) + " -x " + _minimap2Preset_map + " " + contigFilename + " " + readFilename;
		//Utils::executeMinimap2(command, alignFilename);
		auto start2 = high_resolution_clock::now();
		alignReads(contigFilename, readFilename, _minimap2Preset_map, _nbCores);
		Logger::get().debug() << "\tDone " << " (" << duration_cast<seconds>(high_resolution_clock::now() - start2).count() << "s) " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";
		//command += " | gzip -c - > " + alignFilename;
		//cout << command << endl;
		//Utils::executeCommand(command, _tmpDir);

		Logger::get().debug() << "\tLoad contigs";
		start2 = high_resolution_clock::now();
		loadContigs(contigFilename, pass);
		Logger::get().debug() << "\tDone " << " (" << duration_cast<seconds>(high_resolution_clock::now() - start2).count() << "s) " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";

		//loadAlignments(alignFilename, true, true);
		//computeContigCoverages(contigFilename);
		
		//loadAllAlignments(alignFilename, true);
		computeContigCoveragesAll();

		ankerl::unordered_dense::map<u_int32_t, vector<pair<u_int32_t, u_int32_t>>>().swap(_contigIndex_to_alignments);
		
		u_int64_t checksumTotal = 0;

		//for(const auto& it : _alignments){
		//	u_int64_t checksum = ((it.second._readIndex+1) * (it.second._readStart+1) * (it.second._readEnd+1) * (it.second._contigStart+1) * (it.second._contigEnd+1));
		//	checksumTotal += checksum;
			//luls.push_back({it.second._readIndex, checksum});
		//}
		
		for(const auto& it2 : _allAlignments){
			for(const auto& it : it2.second){
				u_int64_t checksum = ((it._readIndex+1) * (it._readStart+1) * (it._readEnd+1) * (it._contigStart+1) * (it._contigEnd+1));
				checksumTotal += checksum;
			}
			//luls.push_back({it.second._readIndex, checksum});
		}

		Logger::get().debug() << "\t\tChecksum alignment: " << checksumTotal;
		
		//cout << "Load contigs" << endl;
		//loadContigs(_inputFilename_contigs, false);

		Logger::get().debug() << "\tCollecting window sequences";
		start2 = high_resolution_clock::now();
		ReadParserParallel readParser(readFilename, true, false, _nbCores);
		readParser.parse(CollectWindowSequencesFunctor(*this));
		Logger::get().debug() << "\tDone " << " (" << duration_cast<seconds>(high_resolution_clock::now() - start2).count() << "s) " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";


		u_int64_t checksum = 0;
		u_int64_t checksum2 = 0;

		for(const auto& it : _contigWindowSequences){

			const u_int32_t contigIndex = it.first;
			const vector<vector<Window>>& contigWindows = it.second;
			u_int64_t checksumContig = 0;

			//cout << "Contig " << contigIndex << " " << contigWindows.size() << endl;

			for(size_t i=0; i<contigWindows.size(); i++){

				const vector<Window>& windows = contigWindows[i];

				u_int64_t checksumLocal = 0;
				//vector<Window> windowsDebug;

				for(size_t j=0; j<windows.size(); j++){
					const Window& window = windows[j];
					
					checksumLocal += (i+1) * window.hash();

					checksum += ((i+1) * window.hash()); 

					//if(contigIndex == 538 && i == 921){
					//	windowsDebug.push_back(window);
					//}
				}

				checksumContig += checksumLocal;
				
				//if(contigIndex == 538){
				//	cout << "\t" << i << "\t" << checksumLocal << endl;
				//}
				//checksum2 += checksumLocal;
				//
				/*
				if(contigIndex == 538 && i == 921){
					std::sort(windowsDebug.begin(), windowsDebug.end(), [&](const Window& a, const Window& b) {
						return a.hash() < b.hash();
					});

					for(size_t i=0; i<windowsDebug.size(); i++){
						
						char* dnaStr = windowsDebug[i]._sequence->to_string();
						string seq = string(dnaStr);
						free(dnaStr);

						cout << i << "\t" << windowsDebug[i].hash() << "\t" << seq << endl;
					}
				}
				*/
				
				//if(contigIndex == 538 && i == 921){
				//	cout << "\t" << i << "\t" << checksumLocal << endl;
				//	getchar();
				//}
				
			}

			//cout << "\tChecksum: " << checksumContig << endl;
			//cout << contigIndex << " " << checksum << endl;
		}

		Logger::get().debug() << "\t\tWindow checksum: " << checksum;
		//cout << "Window checksum2: " << checksum2 << endl;

		//cout << "miuom" << endl;
		//getchar();
		
		Logger::get().debug() << "\tPerform correction";
		start2 = high_resolution_clock::now();
		performCorrection(outputContigFile, file_contigHeaders, pass);
		Logger::get().debug() << "\tDone " << " (" << duration_cast<seconds>(high_resolution_clock::now() - start2).count() << "s) " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";


		Logger::get().debug() << "\tPolishing done " << " (" << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s)";
		

		clearPass();

	}


	void alignReads(const string contigFilename, const string readFilename, const string preset, int nbCores){
		
		//cout << nbCores << endl;
		//string command = "minimap2 -I 100G -v 0 -m 500 -t " + to_string(_nbCores) + " -x " + _minimap2Preset_map + " " + contigFilename + " " + readFilename;
		

		mm_idxopt_t iopt;
		mm_mapopt_t mopt;
		//int n_threads = nbCores;

		mm_verbose = 2; // disable message output to stderr
		mm_set_opt(0, &iopt, &mopt);
		mm_set_opt(preset.c_str(), &iopt, &mopt); //"ava-ont"
		//mopt.flag |= MM_F_CIGAR; // perform alignment
		iopt.batch_size = 0x7fffffffffffffffL; //always build a uni-part index

		mopt.min_chain_score = 500; //-m 500

		// open query file for reading; you may use your favorite FASTA/Q parser
		//gzFile f = gzopen(readFilename.c_str(), "r");
		//assert(f);
		//kseq_t *ks = kseq_init(f);

		// open index reader
		mm_idx_reader_t *r = mm_idx_reader_open(contigFilename.c_str(), &iopt, 0);
		mm_idx_t *mi = mm_idx_reader_read(r, nbCores);
		//cout << (mi != 0) << endl;
		//cout << mi << endl;
		mm_mapopt_update(&mopt, mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!
		//mopt.mid_occ = 5;

		ReadParserParallel readParser(readFilename, true, false, _nbCores);
		readParser.parse(MapReadsFunctor(*this, mi, mopt));
		/*
		while ((mi = mm_idx_reader_read(r, nbCores)) != 0) { // traverse each part of the index
			mm_mapopt_update(&mopt, mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!
			mopt.mid_occ = 5;
			mm_tbuf_t *tbuf = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread
			gzrewind(f);
			kseq_rewind(ks);
			while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence
				mm_reg1_t *reg;
				int j, i, n_reg;
				reg = mm_map(mi, ks->seq.l, ks->seq.s, &n_reg, tbuf, &mopt, 0); // get all hits for the query
				for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
					mm_reg1_t *r = &reg[j];

					//cout << r->qs << " " << r->qe << endl;
					//assert(r->p); // with MM_F_CIGAR, this should not be NULL
					//printf("%s\t%d\t%d\t%d\t%c\t", ks->name.s, ks->seq.l, r->qs, r->qe, "+-"[r->rev]);
					//printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\tcg:Z:", mi->seq[r->rid].name, mi->seq[r->rid].len, r->rs, r->re, r->mlen, r->blen, r->mapq);
					//for (i = 0; i < r->p->n_cigar; ++i) // IMPORTANT: this gives the CIGAR in the aligned regions. NO soft/hard clippings!
					//	printf("%d%c", r->p->cigar[i]>>4, MM_CIGAR_STR[r->p->cigar[i]&0xf]);
					//putchar('\n');
					//free(r->p);
				}
				free(reg);
			}
			mm_tbuf_destroy(tbuf);
			mm_idx_destroy(mi);
		}
		*/
		mm_idx_destroy(mi);
		mm_idx_reader_close(r); // close the index reader
		//kseq_destroy(ks); // close the query file
		//gzclose(f);

	}



	class MapReadsFunctor {

		public:

		ContigPolisher& _parent;
		mm_idx_t* _mi;
		mm_tbuf_t*_tbuf;
		mm_mapopt_t& _mopt;

		MapReadsFunctor(ContigPolisher& parent, mm_idx_t* mi, mm_mapopt_t& mopt) : _parent(parent), _mopt(mopt){
			_tbuf = mm_tbuf_init();
			_mi = mi;
		}

		MapReadsFunctor(const MapReadsFunctor& copy) : _parent(copy._parent), _mopt(copy._mopt){
			_tbuf = mm_tbuf_init();
			_mi = copy._mi;
		}

		~MapReadsFunctor(){
			mm_tbuf_destroy(_tbuf);
		}

		void operator () (const Read& read) {
			/*
			cout << read._index << endl;
			cout << read._seq << endl;
			cout << _mi << endl;

			string preset = "map-ont";
			mm_idxopt_t iopt;
			mm_mapopt_t mopt;
			//int n_threads = nbCores;
			mm_verbose = 2; // disable message output to stderr
			mm_set_opt(0, &iopt, &mopt);
			mm_set_opt(preset.c_str(), &iopt, &mopt); //"ava-ont"
			//mopt.flag |= MM_F_CIGAR; // perform alignment
			iopt.batch_size = 0x7fffffffffffffffL; //always build a uni-part index
			mopt.min_chain_score = 500; //-m 500
			mm_mapopt_update(&mopt, _mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!
			mopt.mid_occ = 5;
			*/
			
			u_int32_t readIndex = stoull(read._header);

			mm_reg1_t *reg;
			int j, i, n_reg;
			reg = mm_map(_mi, read._seq.size(), read._seq.c_str(), &n_reg, _tbuf, &_mopt, 0); // get all hits for the query
			
			#pragma omp critical(loadAlignments)
			{
				for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
					mm_reg1_t *r = &reg[j];

					AlignmentBounds alignment;
					alignment._referenceLength = _mi->seq[r->rid].len;
					alignment._queryLength = read._seq.size();
					alignment._queryStart = r->qs;
					alignment._queryEnd = r->qe;
					alignment._referenceStart = r->rs;
					alignment._referenceEnd = r->re;
					alignment._isReversed = r->rev;
					alignment._identity = ((double)r->mlen) / r->blen;// r->div;
					alignment._nbMatches = r->mlen;

					_parent.loadAllAlignments_read2(readIndex, string(_mi->seq[r->rid].name), alignment);
					//_parent._checksum += 1;
					//cout << r->qs << " " << r->qe << endl;
					//assert(r->p); // with MM_F_CIGAR, this should not be NULL
					//printf("%s\t%d\t%d\t%d\t%c\t", ks->name.s, ks->seq.l, r->qs, r->qe, "+-"[r->rev]);
					//printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\tcg:Z:", mi->seq[r->rid].name, mi->seq[r->rid].len, r->rs, r->re, r->mlen, r->blen, r->mapq);
					//for (i = 0; i < r->p->n_cigar; ++i) // IMPORTANT: this gives the CIGAR in the aligned regions. NO soft/hard clippings!
					//	printf("%d%c", r->p->cigar[i]>>4, MM_CIGAR_STR[r->p->cigar[i]&0xf]);
					//putchar('\n');
					//free(r->p);
				}
			}

			free(reg);
			

			
			/*
			//if(_parent._allAlignments.find(readIndex) == _parent._allAlignments.end()) return;
			if(_parent._alignments.find(readIndex) == _parent._alignments.end()) return;

			//const vector<Alignment>& als = _parent._allAlignments[readIndex];

			//for(const Alignment& al : als){
			const Alignment& al = _parent._alignments[readIndex];
			//for(const Alignment& al : _alignments[readIndex]){
			u_int32_t contigIndex = al._contigIndex;

			if(_parent._contigSequences.find(contigIndex) == _parent._contigSequences.end()) return;
			*/
		}
	};

	void computeContigCoveragesAll(){
		
		//_contigLength.clear();
		_contigCoverages.clear();

		Logger::get().debug() << "\tComputing contig coverages";

		/*
		ankerl::unordered_dense::map<u_int32_t, vector<pair<u_int32_t, u_int32_t>>> contigHits;
		//auto fp = std::bind(&ContigPolisher::computeContigCoverages_setup_read, this, std::placeholders::_1);
		//ReadParser readParser(contigFilename, true, false);
		//readParser.parse(fp);

		
		for(const auto& it : _allAlignments){
		//for(const auto& it : _alignments){
			
			const vector<Alignment>& alignments = it.second;


			for(const Alignment& alignment : alignments){
				//if(_contigLength.find(alignment._contigIndex) == _contigLength.end()) continue;
				contigHits[alignment._contigIndex].push_back({alignment._contigStart, alignment._contigEnd});
			}
			
			//contigHits[it.second._contigIndex].push_back({it.second._contigStart, it.second._contigEnd});

		}
		*/

		for(auto& it : _contigIndex_to_alignments){


			//if(_contigLength.find(it.first) == _contigLength.end()){
			//	cout << "ups: " << it.first << endl;
			//}

			//cout << _contigLength[it.first] << endl;
			vector<u_int32_t> coverages(_contigSequences[it.first].size(), 0);

			for(auto& interval : it.second){

				u_int32_t contigStart = interval.first;
				u_int32_t contigEnd = interval.second;

				
				if(contigStart >= coverages.size()) continue;  //Break intra repeats
				contigEnd = min((u_int32_t)contigEnd, (u_int32_t)coverages.size()); //Break intra repeats

				//cout << "\t" << interval.first << " " << interval.second << endl;
				for(size_t i=contigStart; i<contigEnd; i++){
					coverages[i] += 1;
				}
			}
			if (coverages.size() < 80 *2){
				_contigCoverages[it.first] = 1;
			}
			else{
				u_int64_t sum = 0;
				for(long i=75; i<((long)coverages.size())-75; i++){
					sum += coverages[i];
				}

				float coverage = sum / ((double)coverages.size());  
				_contigCoverages[it.first] = coverage;
			}


		}

		//_contigLength.clear();
	}

	//void computeContigCoverages_setup_read(const Read& read){

		//cout << read._index << "\t" << read._header << "\t" << read._seq.size() << endl;
	//	_contigLength[read._index] = read._seq.size();
		//int nbCounts = read._seq.size() / _contigCoverageWindow;
		//cout << Utils::shortenHeader(read._header) << " " << nbCounts << endl;
		//_contigHitPos[read._index].resize(nbCounts, 0);
	//}

	void trimContigs(const string& inputContigFilename, const string& outputContigFilename, const int& minimapBatchSize){
		
		//ContigCurator contigCurator(_readPartitionDir, _minimap2Preset_map, _minimap2Preset_ava, 0, _minimizerSize, ContigPolisher::_minContigLength, _useMetamdbgHeaderStyle, _nbCores);
		//contigCurator.executeContigTrimming(inputContigFilename, outputContigFilename, minimapBatchSize);
	}

	void pileupContigs(int nbPartitions){

		//auto start = high_resolution_clock::now();

		//ContigPileup contigPileup(_tmpDir, _readPartitionDir, _minimap2Preset_map, _minimap2Preset_ava, nbPartitions, _useMetamdbgHeaderStyle, _minimizerSize, _minContigLength, _maximalMappingOffset, _nbCores);
		//contigPileup.execute();
		
		//Logger::get().debug() << "\tDone " << " (" << duration_cast<seconds>(high_resolution_clock::now() - start).count() << "s)";

	}
	


	
	

	u_int32_t _currentPartition;
	
	/*
	void debugFindContigs(const string& contigFilename){

		_debug_foundContig = false;

		auto fp = std::bind(&ContigPolisher::debugFindContigs_read, this, std::placeholders::_1);
		ReadParser readParser(contigFilename, true, false);
		readParser.parse(fp);
		
	}
	
	void debugFindContigs_read(const Read& read){
		if (read._header.find(_debug_contigName) != std::string::npos) {
			_debug_foundContig = true;
		}
	}
	*/

	
	void loadContigs(const string& contigFilename, int pass){

		_contigSequences.clear();
		_contigHeaders.clear();
		_contigWindowSequences.clear();
		
		auto fp = std::bind(&ContigPolisher::loadContigs_read, this, std::placeholders::_1);
		ReadParser readParser(contigFilename, true, false);
		readParser.parse(fp);
		
		for(const auto& it: _allAlignments){
			for(const auto& alignment: it.second){
				_contigIndex_to_alignments[alignment._contigIndex].push_back({alignment._contigStart, alignment._contigEnd});
			}
		}

		/*
		//cout << "single core here" << endl;
		//ReadParserParallel readParser(contigFilename, true, false, _nbCores);
		//readParser.parse(BreakIntraRepeatFunctor(*this, pass));
		*/

		
	}
	
	void loadContigs_read(const Read& read){

		//if(_debug_contigName != ""){
		//	if (read._header.find(_debug_contigName) != std::string::npos) {
		//	}
		//	else{
		//		return;
		//	}
		//}

		u_int32_t contigIndex = read._index;
		//const string& contigName = Utils::shortenHeader(read._header);

		//if(_currentPartition == 0) _logFile << "Loading contig in partition 0: " << contigIndex << endl;
		_contigSequences[contigIndex] = read._seq;

		size_t nbWindows = ceil(((double)read._seq.size()) / (double)_windowLength);

		//cout << "loadContigs_read: Nb windows + 1 si offsrt polishing" << endl;
		vector<vector<Window>> windows(nbWindows);

		_contigWindowSequences[contigIndex] = windows;
		_contigHeaders[contigIndex] = read._header;

		//cout << "Loaded contig: " << contigIndex << " " << read._header << " " << read._seq.size() << endl;

		for(size_t i=0; i<read._seq.size(); i++){
			_checksum_inputContigs += (read._seq[i]*read._seq.size());
		}

		//cout << "load: " <<  read._header << " " << read._seq.size() << endl;
	}
	


	/*
	void breakIntraRepeats(u_int32_t partition, const string& contigFilename){
		clearPass();

		string readFilename = _readPartitionDir + "/" + to_string(partition) + "_reads";
		indexContigName(contigFilename);

		Logger::get().debug() << "\tMap reads to contigs";

		const string& alignFilename = _readPartitionDir + "/align.paf.gz";
		string command = "minimap2 -I 100G -v 0 -m 500 -t " + to_string(_nbCores) + " -x " + _minimap2Preset_map + " " + contigFilename + " " + readFilename;
		//Utils::executeMinimap2(command, alignFilename);
		auto start2 = high_resolution_clock::now();
		alignReads(contigFilename, readFilename, _minimap2Preset_map, _nbCores);
		Logger::get().debug() << "\tDone " << " (" << duration_cast<seconds>(high_resolution_clock::now() - start2).count() << "s) " << (peakrss() / 1024.0 / 1024.0 / 1024.0) << " GB";
		//command += " | gzip -c - > " + alignFilename;
		//cout << command << endl;
		//Utils::executeCommand(command, _tmpDir);
		
	}
	*/


	class BreakIntraRepeatFunctor {

		public:

		ContigPolisher& _parent;
		mm_tbuf_t*_tbuf;
		int _pass;

		BreakIntraRepeatFunctor(ContigPolisher& parent, int pass) : _parent(parent){
			_tbuf = mm_tbuf_init();
			_pass = pass;
		}

		BreakIntraRepeatFunctor(const BreakIntraRepeatFunctor& copy) : _parent(copy._parent){
			_tbuf = mm_tbuf_init();
			_pass = copy._pass;
		}

		~BreakIntraRepeatFunctor(){
			mm_tbuf_destroy(_tbuf);
		}

		void operator () (const Read& read) {
			
			bool isCircular = false;
			if(read._header[read._header.size()-1] == 'c'){
				isCircular = true;
			}

			const u_int32_t contigIndex2 = _parent._contigName_to_contigIndex[Utils::shortenHeader(read._header)];
			string contigSequence = read._seq;
			
			if(!isCircular && _pass == 0){
					

				while(true){
						
					vector<pair<u_int32_t, u_int32_t>> repeats = findIntraRepeats(contigSequence, _tbuf, _parent._minimap2Preset_map);

					if(repeats.size() > 0){
						for(size_t i=0; i<repeats.size(); i++){
							cout << "Intra repeat: " << read._header << " " << contigSequence.size() << " " << (repeats[i].second-repeats[i].first) << "    " << repeats[i].first << " " << repeats[i].second << endl;
						}
					}

					u_int32_t unbridgedRepeatStart = -1;

					for(size_t i=0; i<repeats.size(); i++){

						u_int32_t repeatStart = repeats[i].first;
						u_int32_t repeatEnd = repeats[i].second;

						u_int32_t repeatStartOffset = 500;
						if(repeatStart < 500) repeatStartOffset = 1; //Disable extra required read alignment if the repeat is at the begining of the contig (this is probably a spurious small repeat)

						bool isRepeatBridged = false;

						if(_parent._contigIndex_to_alignments.find(contigIndex2) != _parent._contigIndex_to_alignments.end()){

							for(const auto& alignment : _parent._contigIndex_to_alignments[contigIndex2]){
								if(alignment.first+repeatStartOffset < repeatStart && alignment.second > repeatEnd+500){
									cout << "\tBridged: " << alignment.first << " " << alignment.second << endl;
									isRepeatBridged = true;
									break;
								}
							}
						}

						if(!isRepeatBridged){
							if(repeatStart < unbridgedRepeatStart){
								unbridgedRepeatStart = repeatStart;
							}
						}


					}

					if(unbridgedRepeatStart != -1){
						cout << "\tYOOOOO Contig: " << read._header << " " << contigIndex2 << " " << contigSequence.size() << "    " << unbridgedRepeatStart << endl;
						cout << "\tRemoving: " << (contigSequence.size() - unbridgedRepeatStart) << endl;
						contigSequence = contigSequence.substr(0, unbridgedRepeatStart);
					}
					
					if(unbridgedRepeatStart == -1) break; //No unsupported repeat found
				}

			}

			#pragma omp critical(loadContig)
			{
				u_int32_t contigIndex = read._index;
				//const string& contigName = Utils::shortenHeader(read._header);

				//if(_currentPartition == 0) _logFile << "Loading contig in partition 0: " << contigIndex << endl;
				_parent._contigSequences[contigIndex] = contigSequence;

				size_t nbWindows = ceil(((double)contigSequence.size()) / (double)_parent._windowLength);

				//cout << "loadContigs_read: Nb windows + 1 si offsrt polishing" << endl;
				vector<vector<Window>> windows(nbWindows);

				_parent._contigWindowSequences[contigIndex] = windows;
				_parent._contigHeaders[contigIndex] = read._header;

				//cout << "Loaded contig: " << contigIndex << " " << read._header << " " << read._seq.size() << endl;

				for(size_t i=0; i<contigSequence.size(); i++){
					_parent._checksum_inputContigs += (contigSequence[i]*contigSequence.size());
				}
			}
		}




		vector<pair<u_int32_t, u_int32_t>> findIntraRepeats(const string& sequence, mm_tbuf_t* tbuf, const string preset){

			vector<pair<u_int32_t, u_int32_t>> repeats;
			string fakeName = "target";
			
			mm_idxopt_t iopt;
			mm_mapopt_t mopt;
			mm_set_opt(0, &iopt, &mopt); //"ava-ont"
			mm_set_opt(preset.c_str(), &iopt, &mopt); //"ava-ont"
			iopt.batch_size = 0x7fffffffffffffffL; //always build a uni-part index

			//if(performBaseLevelAlignment){
			mopt.flag |= MM_F_CIGAR; // perform alignment 
			//}
			
			//mopt.min_chain_score = 500; //-m 500

			mopt.flag |= MM_F_ALL_CHAINS | MM_F_NO_DIAG | MM_F_NO_DUAL | MM_F_NO_LJOIN; // -D -P --no-long-join --dual=no

			mm_idx_t* mi = minimap2index(iopt.w, iopt.k, iopt.flag&1, iopt.bucket_bits, sequence);
			
			mopt.mid_occ = 10;

			mm_mapopt_update(&mopt, mi);

			mm_reg1_t *reg;
			//mm_tbuf_t *tbuf = mm_tbuf_init();
			int j, i, n_reg;
			reg = mm_map(mi, sequence.size(), sequence.c_str(), &n_reg, tbuf, &mopt, fakeName.c_str()); // get all hits for the query

			//u_int32_t maxSelfOverlapLength = 0;

			//if(n_reg > 1) cout << "----" << endl;
			for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
			
				mm_reg1_t *r = &reg[j];

				free(r->p);
				
				if(r->rev) continue; 

				//cout << r->qs << " " << r->qe << "    " << r->rs << " " << r->re << endl;
				if(r->qs > 50) continue; //Looking for alignment at the contig start
				if(r->rs < 2000) continue; //Make sure that the repeat is not just the contig against itself sligtly shifted

				float identity = ((double)r->mlen) / r->blen;
				//cout << r->mlen << " " << r->blen << " " <<  identity << endl;

				if(identity < 0.9) continue;

				repeats.push_back({r->rs, r->re});
				/*
				//if(sequence.size() - r->re > 50) continue; //Looking for alignment at the contig end

				//cout << r->qs << "\t" << r->qe << "\t" << r->rs << "\t" << r->re << endl;



				//cout << "Found self overlap:\t" << r->qs << "\t" << r->qe << "\t" << r->rs << "\t" << r->re << endl;

				u_int32_t overlapStart = r->qe;
				u_int32_t overlapEnd = sequence.size() - r->rs;

				u_int32_t selfOverlapLength = max(overlapStart, overlapEnd);
				if(selfOverlapLength >= sequence.size()) continue;

				if(selfOverlapLength > maxSelfOverlapLength){
					maxSelfOverlapLength = selfOverlapLength;
				}
				*/
			}
			
		
			free(reg);
			mm_idx_destroy(mi);
			
			//cout << repeats.size() << endl;
			return repeats;
		}


		mm_idx_t* minimap2index(int w, int k, int is_hpc, int bucket_bits, const string& sequence){
			const char *seq = sequence.c_str();
			int len = sequence.size();
			const char *fake_name = "target";
			char *s;
			mm_idx_t *mi;
			s = (char*)calloc(len + 1, 1);
			memcpy(s, seq, len);
			mi = mm_idx_str(w, k, is_hpc, bucket_bits, 1, (const char**)&s, (const char**)&fake_name);
			free(s);
			return mi;
		}

	};


	void clearPass(){

		ankerl::unordered_dense::map<ReadType, vector<Alignment>>().swap(_allAlignments); //_allAlignments.clear();
		ankerl::unordered_dense::map<u_int32_t, string>().swap(_contigSequences); //_contigSequences.clear();
		ankerl::unordered_dense::map<u_int32_t, vector<vector<Window>>>().swap(_contigWindowSequences); //_contigWindowSequences.clear();
		ankerl::unordered_dense::map<u_int32_t, string>().swap(_contigHeaders); //_contigHeaders.clear();
		ankerl::unordered_dense::map<string, u_int32_t>().swap(_contigName_to_contigIndex);
		//ankerl::unordered_dense::map<u_int32_t, u_int32_t>().swap(_contigLength);
		ankerl::unordered_dense::map<u_int32_t, float>().swap(_contigCoverages);
		ankerl::unordered_dense::map<u_int32_t, vector<pair<u_int32_t, u_int32_t>>>().swap(_contigIndex_to_alignments);
		//_alignments.clear();
		//clearOverlap();
		
		
		//_validContigIndexes.clear();
		
		
		//_alignments_isReversed.clear();
		//_usedReadIndexes.clear();
		//_contigName_to_contigIndex.clear();
		//_alignments.clear();
	}


	void indexContigName(const string& contigFilename){
		 
		//cout << "Indexing contig names: " << contigFilename << endl;
		//_contigName_to_contigIndex.clear();

		auto fp = std::bind(&ContigPolisher::indexContigName_read, this, std::placeholders::_1);
		ReadParser readParser(contigFilename, true, false);
		readParser.parse(fp);
	}

	void indexContigName_read(const Read& read){

		//cout << "Index contig: " << Utils::shortenHeader(read._header) << " " << read._index << endl;
		_contigName_to_contigIndex[Utils::shortenHeader(read._header)] = read._index;
	}

	//ankerl::unordered_dense::map<u_int32_t, u_int32_t> _contigLength;
	ankerl::unordered_dense::map<u_int32_t, float> _contigCoverages;
	ankerl::unordered_dense::map<ReadType, vector<Alignment>> _allAlignments;
	ankerl::unordered_dense::map<u_int32_t, vector<pair<u_int32_t, u_int32_t>>> _contigIndex_to_alignments;
	//phmap::parallel_flat_hash_map<ReadType, Alignment> _alignments;
	ankerl::unordered_dense::map<u_int32_t, string> _contigSequences;
	ankerl::unordered_dense::map<u_int32_t, vector<vector<Window>>> _contigWindowSequences;
	ankerl::unordered_dense::map<u_int32_t, string> _contigHeaders;
	//unordered_map<ContigRead, u_int32_t, ContigRead_hash> _alignmentCounts;
	//u_int64_t _correctedContigIndex;

	ankerl::unordered_dense::map<string, u_int32_t> _contigName_to_contigIndex;
	//ankerl::unordered_dense::map<string, ReadType> _readName_to_readIndex;

	//bool _indexPerContig;
	bool _useErrorFilter;

	//struct Lul{
	//	ReadType _readIndex;
	//	u_int64_t _hash;
	//};

	void loadAllAlignments(const string& alignFilename, bool useIndexedName){
		/*
		Logger::get().debug() << "\tLoading all alignments";

		//_alignments.clear();
		_allAlignments.clear();

		PafParser pafParser(alignFilename);
		auto fp = std::bind(&ContigPolisher::loadAllAlignments_read, this, std::placeholders::_1);
		pafParser.parse(fp);


		//vector<Lul> luls;
		//ofstream outputFile("./lala.txt");
		u_int64_t checksumTotal = 0;

		for(const auto& it : _alignments){
			u_int64_t checksum = ((it.second._readIndex+1) * (it.second._readStart+1) * (it.second._readEnd+1) * (it.second._contigStart+1) * (it.second._contigEnd+1));
			checksumTotal += checksum;
			//luls.push_back({it.second._readIndex, checksum});
		}

		//std::sort(luls.begin(), luls.end(), [&](const Lul& a, const Lul& b) {
		//	return a._readIndex < b._readIndex;
		//});

		cout << "ContigPolisher checksum alignment: " << checksumTotal << endl;
		*/
		//for(const auto& it : luls){
		//	outputFile << it._readIndex << "\t" << it._hash << endl;
		//}
		//outputFile.close();
		//exit(1);
		
		/*
		cout << "Nb alignment: " << _allAlignments.size() << endl;

		for(const auto& it : _allAlignments){
			if(it.second.size() > 1){
				cout << it.first << " " << it.second[0]._readLength << endl;
				for(const Alignment& al : it.second){
					cout << "\t" <<  al._readStart << " " << al._readEnd << "    " << al._contigIndex << "\t" << al._contigStart << "\t" << al._contigEnd << endl;
				}
			}
		}

		cout << "S2" << endl;
		*/
		//getchar();
		
	}


	void loadAllAlignments_read2(const ReadType readIndex, const string& contigName, const AlignmentBounds& bounds){
		/*
		vector<string> _fields = Utils::split(line, '\t');

		const string& readName = Utils::shortenHeader((_fields)[0]);
		const string& contigName = Utils::shortenHeader((_fields)[5]);
		
		u_int64_t readLength = stoull((_fields)[1]);
		u_int32_t readStart = stoull((_fields)[2]);
		u_int32_t readEnd = stoull((_fields)[3]);
		u_int32_t contigLength = stoull((_fields)[6]);
		u_int32_t contigStart = stoull((_fields)[7]);
		u_int32_t contigEnd = stoull((_fields)[8]);

		u_int64_t nbMatches = stoull((_fields)[9]);
		u_int64_t alignLength = stoull((_fields)[10]);

		bool strand = (_fields)[4] == "-";

		AlignmentBounds bounds;
		bounds._queryStart = readStart;
		bounds._queryEnd = readEnd;
		bounds._referenceStart = contigStart;
		bounds._referenceEnd = contigEnd;
		bounds._queryLength = readLength;
		bounds._referenceLength = contigLength;
		bounds._isReversed = strand;
		*/

		u_int32_t readSizeMappable = bounds.getMappableLength(); //A read size that do not consider alignment outside contig bounds (start and end of the contigs)
		
		float identity = ((long double) bounds._nbMatches) / ((long double) readSizeMappable);

		//u_int32_t contigIndex = _contigName_to_contigIndex[contigName];
		//ReadType readIndex = _readName_to_readIndex[readName];

		//u_int32_t readIndex = stoull(readName);;

		if(bounds._isReversed) return;
		if(_contigName_to_contigIndex.find(contigName) == _contigName_to_contigIndex.end()) return;

		//int64_t overhang = readSizeMappable - (bounds._queryEnd - bounds._queryStart);

		//if(overhang > 100) return;
		//float identity2 = ((long double) bounds._queryEnd-bounds._queryStart) / ((long double) readSizeMappable);


		//if(identity2 < 0.99) return;

		

		u_int32_t contigIndex = _contigName_to_contigIndex[contigName];

		Alignment alignment = {contigIndex, readIndex, bounds._isReversed, bounds._queryStart, bounds._queryEnd, bounds._referenceStart, bounds._referenceEnd, identity, bounds._queryLength, bounds._referenceLength}; //, score

		if(!alignment.isMaximalMapping(_maximalMappingOffset)) return;

		indexReadAlignment(readIndex, alignment);

		//if(readIndex == 1568525){
		//	cout << readIndex << "\t" << contigIndex << "\t" << contigStart << "\t" << contigEnd << "\t" << readStart << "\t" << readEnd << endl;
		//}

		/*
		if(_alignments.find(readIndex) == _alignments.end()){
			_alignments[readIndex] = alignment;
		}
		else{

			if(alignment.score() > _alignments[readIndex].score()){
				_alignments[readIndex] = alignment;
			}
			else if(alignment.score() == _alignments[readIndex].score()){ //arbitrary checks to solve tie
				if(alignment._contigStart < _alignments[readIndex]._contigStart){
					_alignments[readIndex] = alignment;
				}
				else if(alignment._contigEnd < _alignments[readIndex]._contigEnd){
					_alignments[readIndex] = alignment;
				}
				else if(alignment._readStart < _alignments[readIndex]._readStart){
					_alignments[readIndex] = alignment;
				}
				else if(alignment._readEnd < _alignments[readIndex]._readEnd){
					_alignments[readIndex] = alignment;
				}
			}
			
		}
		*/

	}

	/*
	void loadAllAlignments_read(const string& line){

		vector<string> _fields = Utils::split(line, '\t');

		const string& readName = Utils::shortenHeader((_fields)[0]);
		const string& contigName = Utils::shortenHeader((_fields)[5]);
		
		u_int64_t readLength = stoull((_fields)[1]);
		u_int32_t readStart = stoull((_fields)[2]);
		u_int32_t readEnd = stoull((_fields)[3]);
		u_int32_t contigLength = stoull((_fields)[6]);
		u_int32_t contigStart = stoull((_fields)[7]);
		u_int32_t contigEnd = stoull((_fields)[8]);

		u_int64_t nbMatches = stoull((_fields)[9]);
		u_int64_t alignLength = stoull((_fields)[10]);

		bool strand = (_fields)[4] == "-";

		AlignmentBounds bounds;
		bounds._queryStart = readStart;
		bounds._queryEnd = readEnd;
		bounds._referenceStart = contigStart;
		bounds._referenceEnd = contigEnd;
		bounds._queryLength = readLength;
		bounds._referenceLength = contigLength;
		bounds._isReversed = strand;


		u_int32_t readSizeMappable = bounds.getMappableLength(); //A read size that do not consider alignment outside contig bounds (start and end of the contigs)
		
		float identity = ((long double) nbMatches) / ((long double) readSizeMappable);

		//u_int32_t contigIndex = _contigName_to_contigIndex[contigName];
		//ReadType readIndex = _readName_to_readIndex[readName];

		u_int32_t readIndex = stoull(readName);;

		if(strand) return;
		if(_contigName_to_contigIndex.find(contigName) == _contigName_to_contigIndex.end()) return;

		//int64_t overhang = readSizeMappable - (bounds._queryEnd - bounds._queryStart);

		//if(overhang > 100) return;
		//float identity2 = ((long double) bounds._queryEnd-bounds._queryStart) / ((long double) readSizeMappable);


		//if(identity2 < 0.99) return;

		

		u_int32_t contigIndex = _contigName_to_contigIndex[contigName];

		Alignment alignment = {contigIndex, readIndex, strand, readStart, readEnd, contigStart, contigEnd, identity, readLength, contigLength}; //, score

		if(!alignment.isMaximalMapping(_maximalMappingOffset)) return;

		//indexReadAlignment(readIndex, alignment);

		//if(readIndex == 1568525){
		//	cout << readIndex << "\t" << contigIndex << "\t" << contigStart << "\t" << contigEnd << "\t" << readStart << "\t" << readEnd << endl;
		//}

		if(_alignments.find(readIndex) == _alignments.end()){
			_alignments[readIndex] = alignment;
		}
		else{

			if(alignment.score() > _alignments[readIndex].score()){
				_alignments[readIndex] = alignment;
			}
			else if(alignment.score() == _alignments[readIndex].score()){ //arbitrary checks to solve tie
				if(alignment._contigStart < _alignments[readIndex]._contigStart){
					_alignments[readIndex] = alignment;
				}
				else if(alignment._contigEnd < _alignments[readIndex]._contigEnd){
					_alignments[readIndex] = alignment;
				}
				else if(alignment._readStart < _alignments[readIndex]._readStart){
					_alignments[readIndex] = alignment;
				}
				else if(alignment._readEnd < _alignments[readIndex]._readEnd){
					_alignments[readIndex] = alignment;
				}
			}
			
		}

	}
	*/
	
	
	void indexReadAlignment(const ReadType& readIndex, const Alignment& alignment){
		

		if(_allAlignments.find(readIndex) == _allAlignments.end()){
			_allAlignments[readIndex].push_back(alignment);
			return;
		}


		vector<Alignment>& existingAlignments = _allAlignments[readIndex];

		vector<Alignment>::iterator it = existingAlignments.begin();
		bool isBetterAlignment = false;
		bool hasOverlap = false;
		bool overlapWithBetterAlignment = false;
		//bool canAddAlignment = true;


		for(const Alignment& existingAlignment : existingAlignments){

			if(alignmentOverlapExistingAlignment(alignment, existingAlignment)){

				if(alignment.score() < existingAlignment.score()){
					overlapWithBetterAlignment = true;
				}

				hasOverlap = true;
			}

		}

		if(overlapWithBetterAlignment){
			return;
		}

		while(it != existingAlignments.end()) {

			const Alignment& existingAlignment = *it;



			if(alignmentOverlapExistingAlignment(alignment, existingAlignment) && ((alignment.score() > existingAlignment.score()) || (alignment.score() == existingAlignment.score() && alignment._readIndex > existingAlignment._readIndex))  ) {
				it = existingAlignments.erase(it);
				isBetterAlignment = true;
			}
			else{
				++it;
			}
		}



		if(isBetterAlignment){
			_allAlignments[readIndex].push_back(alignment);
		}

		if(!isBetterAlignment && !hasOverlap){
			_allAlignments[readIndex].push_back(alignment);
		}
	}
	
	bool alignmentOverlapExistingAlignment(const Alignment& alignment, const Alignment& existingAlignment){

		float allowedOverlap = 500;

		if(alignment._readStart >= existingAlignment._readStart && alignment._readEnd <= existingAlignment._readEnd) return true; //alignment contained
		if(alignment._readStart <= existingAlignment._readStart && alignment._readEnd >= existingAlignment._readEnd) return true; //existing alignment contained
		//if(alignment._readStart >= existingAlignment._readStart && (alignment._readStart+allowedOverlap) <= existingAlignment._readEnd) return true;
		//if(alignment._readEnd >= (existingAlignment._readStart+allowedOverlap) && alignment._readEnd <= existingAlignment._readEnd) return true;

		if(alignment._readStart >= existingAlignment._readStart){
			if(existingAlignment._readEnd - alignment._readStart > allowedOverlap) return true;
		}

		if(alignment._readEnd <= existingAlignment._readEnd){
			if(alignment._readEnd - existingAlignment._readStart > allowedOverlap) return true;
		}
		//int64_t startPos = max(alignment._readStart, existingAlignment._readStart);
		//int64_t endPos = min(alignment._readEnd, existingAlignment._readEnd);

		//int64_t overlapLength = endPos - startPos;

		//if(overlapLength > allowedOverlap) return true;

		return false;
	}
	


	class CollectWindowSequencesFunctor {

		public:

		ContigPolisher& _parent;


		CollectWindowSequencesFunctor(ContigPolisher& parent) : _parent(parent){
		}

		CollectWindowSequencesFunctor(const CollectWindowSequencesFunctor& copy) : _parent(copy._parent){
			
		}

		~CollectWindowSequencesFunctor(){
		}

		void operator () (const Read& read) {

			u_int32_t readIndex = stoull(read._header);

			if(_parent._allAlignments.find(readIndex) == _parent._allAlignments.end()) return;
			//if(_parent._alignments.find(readIndex) == _parent._alignments.end()) return;

			vector<Alignment>& als = _parent._allAlignments[readIndex];

			for(Alignment& al : als){
			//const Alignment& al = _parent._alignments[readIndex];
			//for(const Alignment& al : _alignments[readIndex]){
				u_int32_t contigIndex = al._contigIndex;

				if(_parent._contigSequences.find(contigIndex) == _parent._contigSequences.end()) continue;

				//_logFile << read._seq.size() << " " << read._qual.size() << " " << _contigSequences[contigIndex].size() << " " << al._readStart << " " << al._readEnd << " " << al._contigStart << " " << al._contigEnd << endl;
				string readSeq = read._seq;
				string qualSeq = read._qual;
				string readSequence = readSeq.substr(al._readStart, al._readEnd-al._readStart);


				//if(al._strand){
				//	Utils::toReverseComplement(readSequence);
				//	Utils::toReverseComplement(readSeq);
				//	std::reverse(qualSeq.begin(), qualSeq.end());
				//}
				

				AlignmentBounds bounds;
				bounds._queryStart = al._readStart;
				bounds._queryEnd = al._readEnd;
				bounds._referenceStart = al._contigStart;
				bounds._referenceEnd = al._contigEnd;
				bounds._queryLength = readSeq.size();
				bounds._referenceLength = _parent._contigSequences[contigIndex].size();
				bounds._isReversed = al._strand;


				if(al._contigStart >= bounds._referenceLength) continue;  //Break intra repeats
				al._contigEnd = min((u_int32_t)al._contigEnd, (u_int32_t)bounds._referenceLength); //Break intra repeats

				string contigSequence = _parent._contigSequences[contigIndex].substr(al._contigStart, al._contigEnd-al._contigStart);

				u_int32_t readSizeMappable = bounds.getMappableLength(); //A read size that do not consider alignment outside contig bounds (start and end of the contigs)
				
				


				//if(_contigPolisher.getMaxhang(bounds) > ContigPolisher::maxHang*2) return;

				//u_int32_t readSizeMappable = _contigPolisher.getMappableLength(bounds); //A read size that do not consider alignment outside contig bounds (start and end of the contigs)
				

				//int64_t overhang = readSeq.size() - readSizeMappable;
				//cout << overhang << endl;
				//if(overhang > 10) return;

				//cout << overhang << endl;


				//_logFile << readSequence << endl;
				//_logFile << contigSequence << endl;

				//_logFile << contigSequence.size() << " "<< readSequence.size() << endl;
				static EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0);


				EdlibAlignResult result = edlibAlign(readSequence.c_str(), readSequence.size(), contigSequence.c_str(), contigSequence.size(), config);


				char* cigar;

				if (result.status == EDLIB_STATUS_OK) {
					cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
				} else {
					Logger::get().error() << "Invalid edlib results";
					exit(1);
				}

				//_logFile << cigar << endl;

				edlibFreeAlignResult(result);
				
				int64_t nbMatches = Commons::getCigarNbMatches(cigar, readSequence, contigSequence);
				float identity = ((long double) nbMatches) / ((long double) readSizeMappable);
				if(identity < 0.9) continue;
				//if(al._contigEnd < 20000){
					//cout << "-----" << endl;
					//cout << al._contigStart << " " << al._contigEnd << "     " << al._readStart << " " << al._readEnd << " " << al._readLength << endl;
					//cout << _contigPolisher.getMaxhang(bounds) << endl;
					//u_int64_t nbMatches = getCigarNbMatches(al, cigar, readSequence, contigSequence);
					//getchar();

				//}

				find_breaking_points_from_cigar(_parent._windowLength, al, readSeq.size(), cigar, readSeq, qualSeq, identity);


				free(cigar);

			}
		}

		void find_breaking_points_from_cigar(uint32_t window_length, const Alignment& al, u_int64_t readLength, char* cigar_, const string& readSequence, const string& qualSequence, const float& identity)
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

				if(breaking_points_[j].second >= readSequence.size()) return;
				if(breaking_points_[j + 1].second >= readSequence.size()) return;

				if (breaking_points_[j + 1].second - breaking_points_[j].second < 0.02 * _parent._windowLength) {
					continue;
				}

				
				if (qualSequence.size() > 0) {

					//const auto& quality = overlaps[i]->strand() ? sequence->reverse_quality() : sequence->quality();
					double average_quality = 0;
					for (uint32_t k = breaking_points_[j].second; k < breaking_points_[j + 1].second; ++k) {
						average_quality += static_cast<uint32_t>(qualSequence[k]) - 33;
					}
					average_quality /= breaking_points_[j + 1].second - breaking_points_[j].second;

					if (average_quality < _parent._qualityThreshold) {
						continue;
					}
				}
				

				//uint64_t window_id = id_to_first_window_id[overlaps[i]->t_id()] +
				uint64_t window_id = breaking_points_[j].first / _parent._windowLength;
				uint32_t window_start = (breaking_points_[j].first / _parent._windowLength) * _parent._windowLength;

				//const char* data = overlaps[i]->strand() ?
				//	&(sequence->reverse_complement()[breaking_points[j].second]) :
				//	&(sequence->data()[breaking_points[j].second]);
				const char* data = &readSequence[breaking_points_[j].second];
				uint32_t data_length = breaking_points_[j + 1].second - breaking_points_[j].second;

				string sequence = string(data, data_length);

				/*
				if(al._contigIndex == 1 && window_id > 1 && window_id < 1000){
					if(window_id == 192 || window_id == 193){
						_logFile << readSequence << endl;
						_logFile << sequence << endl;
						if(breaking_points_[j].second > 1) _logFile << readSequence[breaking_points_[j].second-1] << " " << readSequence[breaking_points_[j].second] << endl;
						if(readSequence.size() > readSequence[breaking_points_[j].second+data_length-1]) _logFile << readSequence[breaking_points_[j].second+data_length-1] << " " << readSequence[breaking_points_[j].second+data_length] << endl;
						_logFile << window_id << endl;
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
						_logFile << readSequence << endl;
						_logFile << sequence << endl;
						if(breaking_points_[j].second > 1){
							_logFile << readSequence[breaking_points_[j].second-1] << " " << readSequence[breaking_points_[j].second] << endl;
							if(readSequence[breaking_points_[j].second-1] == readSequence[breaking_points_[j].second]) getchar();
						}
						if(readSequence.size() > readSequence[breaking_points_[j].second+data_length-1]) _logFile << readSequence[breaking_points_[j].second+data_length-1] << " " << readSequence[breaking_points_[j].second+data_length] << endl;
						_logFile << window_id << endl;
						//getchar();
					}
					

				}
				*/

                u_int32_t posStart = breaking_points_[j].first - window_start;
                u_int32_t posEnd =  breaking_points_[j + 1].first - window_start - 1;

				string quality = "";
				if(_parent._useQual && qualSequence.size() > 0){
					quality = string(&qualSequence[breaking_points_[j].second], data_length);
				}

				//if(al._contigIndex == 0 && window_id == 161){
					
					//_logFile << al._contigIndex << " " << al._readIndex << endl;
					//_logFile << (breaking_points_[j + 1].second - breaking_points_[j].second) << " " << (0.02 * _windowLength) << endl;
					//_logFile << sequence << endl;
					//_logFile << quality << endl;
					//getchar();
				//}

				//if(sequence.size() < 490) continue;
				indexWindow(al, window_id, posStart, posEnd, sequence, quality, identity);

				//if(window_id == 24){
				//	cout << sequence << endl;
				//	getchar();
				//}

				//_logFile << window_id << " " << posStart << " " << posEnd << endl;
				//_logFile << sequence << endl;
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

				//_logFile << contigWindowStart << " " << contigWindowEnd  << "      " << readWindowStart << " " << readWindowEnd << endl;

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
			//	_logFile << breakPoint.first << " " << breakPoint.second << endl;
			//}
		}

		//void indexWindow(const Alignment& al, size_t readWindowStart, size_t readWindowEnd, size_t contigWindowStart, size_t contigWindowEnd, const string& windowSequence, const string& windowQualities){
		void indexWindow(const Alignment& al, u_int64_t windowIndex, u_int32_t posStart, u_int32_t posEnd, const string& windowSequence, const string& windowQualities, const float& identity){

			//bool print2 = false;

			//if(windowSequence == "CATCGACCGAAGCTGTTCTGCAAATAAAATTGGACAGCGATAAAAAAGATCAGATATTGAATATTGTGTTTGGGGAACAGATTGCGGAATTGTTGAAATCGAAATGATTTCCTTTTGGCACGGAAGGAGAGGCTCGAACTCCCGACACCTGGTTTTGGAGACCAGTGCTCTACCAACTGAGCTACTTCCGTGTTTGCGGGTGCAAAGGTAAGAGAATATTCTATAATTCCAACTATAT"){
				
			//	Window newWindowDebug = {new DnaBitset2(windowSequence), windowQualities, posStart, posEnd, 0};
			//	cout << al._contigIndex << " " << al._readIndex << " " << windowIndex << endl;
			//	cout << windowQualities << endl;
			//	cout << "hash: " << newWindowDebug.hash() << endl;
			//	print2 = true;
			//}


			/*
			#pragma omp critical(indexWindow)
			{
				
				//size_t contigWindowIndex = contigWindowStart / _windowLength;
				vector<Window>& windows = _parent._contigWindowSequences[al._contigIndex][windowIndex];

				windows.push_back({new DnaBitset2(windowSequence), windowQualities, posStart, posEnd, identity}); //al.score()
			}
			*/

			
			#pragma omp critical(indexWindow)
			{
				//size_t contigWindowIndex = contigWindowStart / _windowLength;
				vector<Window>& windows = _parent._contigWindowSequences[al._contigIndex][windowIndex];

				
				bool interrupt = false;
				if(_parent._maxWindowCopies == 0 || windows.size() < (_parent._maxWindowCopies-1)){


					windows.push_back({DnaBitset2(windowSequence), windowQualities, posStart, posEnd, identity}); //al.score()



					interrupt = true;
				}

				
				
				if(!interrupt){

					bool isIncompleteWindow = false;
					float score = identity; //al.score();

					
					//Window newWindowDebug = {new DnaBitset2(windowSequence), windowQualities, posStart, posEnd, score};
					//bool print = false;
					
					/*
					if(newWindowDebug.hash() == 2819066){
						cout << "a" << endl;
						cout << newWindowDebug.hash() << endl;
						cout << windows.size() << endl;
						cout << al._contigIndex << " " << al._readIndex << " " << windowIndex << endl;
						print = true;
					}
					*/

					//if(newWindowDebug.hash() == 2836249){
					//	cout << "b" << endl;
					//	cout << newWindowDebug.hash() << endl;
					//	cout << windows.size() << endl;
					//	cout << al._contigIndex << " " << al._readIndex << endl;
					//	print = true;
					//}



					u_int64_t currentWindowDistance =  abs(((long)windowSequence.size()) - ((long)_parent._windowLength));
					if(currentWindowDistance > _parent._windowLengthVariance) isIncompleteWindow = true;

					//if(print){
					//	cout << "cc: " << currentWindowDistance << endl;
					//}

					u_int64_t incompleteWindowIndex = -1;
					//u_int64_t minWindowSize = -1;
					u_int64_t largerDistanceWindow = 0;
					

					for(size_t i=0; i<windows.size(); i++){

						const Window& window = windows[i];
						u_int64_t distance = abs(((long)window._sequence._m_len) - ((long)_parent._windowLength));

						if(distance < currentWindowDistance) continue;

						if(distance > _parent._windowLengthVariance){
							if(distance > largerDistanceWindow){
								largerDistanceWindow = distance;
								incompleteWindowIndex = i;

								//if(print2){
								//	cout << "XX: " << largerDistanceWindow << " " << window._sequence->m_len << endl;
								//}
							}
							else if(distance == largerDistanceWindow){
								if(window.hash() > windows[incompleteWindowIndex].hash()){
									largerDistanceWindow = distance;
									incompleteWindowIndex = i;
									
									//if(print2){
									//	cout << "ZZ: " << largerDistanceWindow << " " << window._sequence->m_len << endl;
									//}
									
								}
							}
						}
					}
					

					//if(print){
					//	cout << "dd: " << incompleteWindowIndex << " " << largerDistanceWindow  << endl;
					//}

					//if(print2){
					//	cout << windowSequence << "\t" << currentWindowDistance << "\t" << largerDistanceWindow << "\t" << score << "\t" << incompleteWindowIndex << "\t" << isIncompleteWindow << endl;
					//}

					//u_int64_t distance = abs(((long)windowSequence.size()) - ((long)_windowLength));

					if(incompleteWindowIndex != -1){
						
						//if(print2){
						//	cout << "a" << endl;
						//}

						if(largerDistanceWindow == currentWindowDistance){
							
							Window& window = windows[incompleteWindowIndex];
							Window newWindow = {DnaBitset2(windowSequence), windowQualities, posStart, posEnd, score};

							if(newWindow.hash() < window.hash()){

								//if(al._contigIndex == 538 && windowIndex == 921){
								//	cout << "c:\tAdd: " << newWindow._sequence->m_len << "\tRemove: " << window._sequence->m_len << endl;
								//}

								//if(print2){
								//	cout << "b" << endl;
								//	cout << "Remove: " << window._sequence->m_len << endl;
								//}


								//if(1336573 == windows[incompleteWindowIndex].hash()){
								//	cout << "derp a" << endl;
								//		getchar();
								//}

								//delete window._sequence;
								windows[incompleteWindowIndex] = newWindow;
							}
							else{
								
								//if(al._contigIndex == 538 && windowIndex == 921){
								//	cout << "b:\tIgnore: " << newWindow._sequence->m_len << endl;
								//}

								//if(print2){
								//	cout << "c" << endl;
								//	cout << "Remove: " << newWindow._sequence->m_len << endl;
								//}
								//delete newWindow._sequence;
							}
						}
						else{
							

							Window& window = windows[incompleteWindowIndex];

							//if(al._contigIndex == 538 && windowIndex == 921){
							//	cout << "a:\tAdd: " << windowSequence.size() << "\tRemove: " << window._sequence->m_len << " " << window.hash() << endl;
							//}

							//if(1336573 == windows[incompleteWindowIndex].hash()){
							//	cout << "derp b: " << windowSequence << endl;
							//		getchar();
							//}

							//if(print2){
							//	cout << "d" << endl;
							//	cout << "Remove: " << window._sequence->m_len << endl;
							//}

							//delete window._sequence;
							windows[incompleteWindowIndex] = {DnaBitset2(windowSequence), windowQualities, posStart, posEnd, score};
						}
					}
					else if(!isIncompleteWindow){
						
						//if(al._contigIndex == 538 && windowIndex == 921){
						//	cout << "e" << endl;
						//}
						//u_int64_t largestAligmenet = 0;
						//u_int64_t largestAligmenetIndex = 0;

						
						
        				static float maxVal = std::numeric_limits<float>::max();
						size_t largerWindowIndex = 0;
						//u_int64_t largerDistanceWindow = 0;
						float lowestScore = maxVal;
						
						for(size_t i=0; i<windows.size(); i++){

							const Window& window = windows[i];
							//u_int64_t distance = abs(((long)window._sequence->m_len) - ((long)_windowLength));
							//float score =

							if(window._score < lowestScore){
								lowestScore = window._score;
								largerWindowIndex = i;
							}
							else if(window._score == lowestScore){
								if(window.hash() > windows[largerWindowIndex].hash()){
									lowestScore = window._score;
									largerWindowIndex = i;
								}
							}
						}

						//u_int64_t distance = abs(((long)windowSequence.size()) - ((long)_windowLength));

						if(score == lowestScore){
							
							Window& window = windows[largerWindowIndex];
							Window newWindow = {DnaBitset2(windowSequence), windowQualities, posStart, posEnd, score};

							if(newWindow.hash() < window.hash()){
								
								//if(al._contigIndex == 538 && windowIndex == 921){
								//	cout << "f" << endl;
								//}
								//if(1336573 == windows[largerWindowIndex].hash()){
								//	cout << "derp c" << endl;
								//	getchar();
								//}
								//delete window._sequence;
								windows[largerWindowIndex] = newWindow;
							}
							else{
								
								//if(al._contigIndex == 538 && windowIndex == 921){
								//	cout << "g" << endl;
								//}
								//delete newWindow._sequence;
							}
						}
						else if(score > lowestScore){
							
							//if(1336573 == windows[largerWindowIndex].hash()){
							//	cout << "derp d" << endl;
							//		getchar();
							//}
							//if(al._contigIndex == 538 && windowIndex == 921){
							//	cout << "h" << endl;
							//}
							Window& window = windows[largerWindowIndex];
							//delete window._sequence;
							windows[largerWindowIndex] = {DnaBitset2(windowSequence), windowQualities, posStart, posEnd, score};
						}
						
					}



				} 




			}

			
			//_logFile << contigWindowIndex << endl;


			

			/*
			if(contigWindowIndex == _contigWindowSequences[al._contigIndex].size()-1){

				for(size_t i=0; i<windowSequences.size(); i++){
					_logFile << windowSequences[i]->m_len << " ";
				}
				_logFile << endl;
				//_logFile << contigWindowEnd-contigWindowStart << endl;
				//_logFile << windowSequence << endl;
			}
			*/
		}


	};

	struct CorrectedWindow{
		size_t _windowIndex;
		DnaBitset2 _correctedSequence;
		bool _success;
	};

	static bool CorrectedWindowComparator (const CorrectedWindow& p1, const CorrectedWindow& p2){
		return p1._windowIndex < p2._windowIndex;
	}


	ankerl::unordered_dense::map<u_int32_t, vector<CorrectedWindow>> _currentContigs;



	void performCorrection(BGZF* outputContigFile, ofstream& file_contigHeaders, const int& pass){
		
		u_int64_t checksum = 0;

		vector<u_int32_t> contigIndexes;
		for(auto& it : _contigWindowSequences){
			contigIndexes.push_back(it.first);
		}
		//vector<u_int32_t> contigIndexes;
		//for(auto& it : _contigWindowSequences){
		//	contigIndexes.push_back(it.first);
		//}


		size_t i = 0;
		size_t windowIndex = 0;
		//unordered_map<string, vector<vector<Window>>> _contigWindowSequences;

		//cout << "Perform correction single core" << endl;
		
		//cout << "single core here" << endl;
		#pragma omp parallel num_threads(_nbCores)
		{

			bool isEOF = false;
			size_t contigIndexLocal;
			size_t windowIndexLocal;
			std::unique_ptr<spoa::AlignmentEngine> alignmentEngine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -4);
			/*
			abpoa_para_t *abpt = abpoa_init_para();

			abpt->out_msa = 0; // generate Row-Column multiple sequence alignment(RC-MSA), set 0 to disable
			abpt->out_cons = 1; // generate consensus sequence, set 0 to disable
			abpt->w = 6, abpt->k = 9; abpt->min_w = 10; // minimizer-based seeding and partition
			abpt->progressive_poa = 0;
			abpt->max_n_cons = 1; // to generate 2 consensus sequences
			abpt->use_qv = 1; //Use weight ?

			abpoa_post_set_para(abpt);
			*/

			while(true){


				#pragma omp critical(performCorrection)
				{

					if(i >= _contigWindowSequences.size()){
						isEOF = true;
					}

					if(!isEOF){

						u_int32_t contigIndex = contigIndexes[i];
						//cout << "Contig: " << contigIndex << endl;
						//cout << windowIndex << " " << _contigWindowSequences[contigIndex].size() << endl;
						if(windowIndex >= _contigWindowSequences[contigIndex].size()){
							i += 1;
							windowIndex = 0;
							//cout << "inc contig index" << endl;
						}

						if(i >= _contigWindowSequences.size()){
							isEOF = true;
						}


						if(!isEOF){
							contigIndex = contigIndexes[i];
							contigIndexLocal = contigIndex;
							windowIndexLocal = windowIndex;
							windowIndex += 1;
						}
					}
					
				}

				if(isEOF) break;
				//functorSub(read);

				//#pragma omp critical(lala)
				//{
				//	cout << contigIndexLocal << ": " << windowIndexLocal << "/" << _contigWindowSequences[contigIndexLocal].size() << endl;
				//	getchar();
				//}

				//const string& contigName = contigIndexes[contigIndexLocal];

				//cout << _contigWindowSequences.size() << " " << contigIndexLocal << endl;
				//cout << "lala: " << contigIndexLocal << " " << endl;
				
				//if(_contigWindowSequences.find(contigIndexLocal) == _contigWindowSequences.end()){
				//	cout << "omg 1 " << endl;
				//}

				//if(windowIndexLocal >= _contigWindowSequences[contigIndexLocal].size()){
				//	cout << "omg 2" << endl;
				//	cout << _contigHeaders[contigIndexLocal] << endl;
				//	cout << _contigSequences[contigIndexLocal].size() << endl;
				//	cout << contigIndexLocal << " " << windowIndexLocal << " " << _contigWindowSequences[contigIndexLocal].size() << endl;
				//}

				
				if(windowIndexLocal >= _contigWindowSequences[contigIndexLocal].size()){
					continue;
					//addCorrectedWindow(false, new DnaBitset2(contigOriginalSequence), contigIndexLocal, windowIndexLocal, outputContigFile, pass);	
					//correctedWindows[w] = new DnaBitset2(contigOriginalSequence);
					//_logFile << "No sequences for window" << endl;
					//continue;
				}

				u_int64_t wStart = windowIndexLocal*_windowLength;
				u_int64_t wEnd = min(_contigSequences[contigIndexLocal].size(), (size_t)(wStart+_windowLength));
				string contigOriginalSequence = _contigSequences[contigIndexLocal].substr(wStart, wEnd-wStart);

				//cout << "---" << endl;
				//cout << contigOriginalSequence << endl;
				string correctedSequence1 = computeConsensus(contigIndexLocal, windowIndexLocal, contigOriginalSequence, alignmentEngine);
				//cout << correctedSequence1 << endl;
				//string correctedSequence2 = computeConsensus(contigIndexLocal, windowIndexLocal, abpt, correctedSequence1);
				//cout << correctedSequence2 << endl;

				//if(correctedSequence1 != correctedSequence2){
				//	cout << "yaya" << endl;
				//	getchar();
				//}

			
				addCorrectedWindow(true, DnaBitset2(correctedSequence1), contigIndexLocal, windowIndexLocal, outputContigFile, file_contigHeaders, pass);
				//getchar();
				//correctedWindows[w] = new DnaBitset2(correctedSequence);

			}

			
			//abpoa_free_para(abpt);
		}
		
	}

	/*
	string performCorrection2(const vector<Window>& windows){


        para_t* para = new para_t();
        para->match = 2; //match;
        para->mismatch = -4; //mismatch;
        para->gap_ext1 = -2; //gap;
        para->gap_open1 = 2 * para->gap_ext1;
        para->b = 25;
        para->f = 40;
        para->ab_band = false;
        para->result = 0;
        para->enable_seeding = false;
        para->thread = 1;
        // para->verbose = 1;
        initPara(para);
        //alignment_engines_.emplace_back(para);

		cout << "lul" << endl;
		string consensus = generate_consensus2(para, true, windows);
		cout << consensus << endl;
		//getchar();

        delete para;
        //alignment_engines_[i] = nullptr;

		return consensus;

	}

	string generate_consensus2(para_t* para, bool trim, const vector<Window>& windows) {

		vector<string> sequences;
		//std::vector<std::pair<const char*, uint32_t>> sequences_;
		std::vector<std::pair<const char*, uint32_t>> qualities_;
		std::vector<std::pair<uint32_t, uint32_t>> positions_;

		if (windows.size() < 3) {
			return "";
		}

		cout << "todo: verifier qu'on a le window contig en front() position" << endl;
		cout << "a" << endl;
		cout << sequences.size() << endl;
		for(size_t i=0; i<windows.size(); i++){ 

			//size_t i = order[ii];
			const Window& window = windows[i];
			//const DnaBitset2* dna = variant._sequence; //sequenceCopies[s._sequenceIndex];
			string dnaStr = window._sequence.to_string();

			//cout << dnaStr << endl;

			if(dnaStr.size() < 490) continue;

			sequences.push_back(dnaStr);
			qualities_.emplace_back(window._quality.c_str(), window._quality.size());
			positions_.emplace_back(window._posStart, window._posEnd);

			//cout << "Add: " << dnaStr << " " << window._quality << " " << window._posStart << " " << window._posEnd << endl;
		}

		std::sort(sequences.begin(), sequences.end(), [&](const string& a, const string& b) {
			return a.size() > b.size();
		});


		if(sequences.front().size() < 490) return "";
		if (sequences.size() < 3) {
			//return std::string(sequences_.front().first, sequences_.front().second);
			//return false;
			return "";
		}

		cout << "b" << endl;
		graph* DAG = new graph();
		DAG->init(para);
		aligned_buff_t* mpool = new aligned_buff_t;
		
		cout << sequences.front() << endl;
		cout << "c" << endl;
		std::vector<res_t> res = alignment(para, DAG, nullptr, 0, sequences.front().c_str(), sequences.front().size(), mpool);
		DAG->add_path(para->m, 0, res, 1);
		DAG->topsort(para, 0);

		cout << "d" << endl;
		std::vector<uint32_t> rank;
		rank.reserve(sequences.size());
		for (uint32_t i = 0; i < sequences.size(); ++i) {
			rank.emplace_back(i);
		}

		cout << "e" << endl;
		std::sort(rank.begin() + 1, rank.end(), [&](uint32_t lhs, uint32_t rhs) {
			return positions_[lhs].first < positions_[rhs].first; });
			
		cout << "f" << endl;
		// perform minipoa
		uint32_t offset = 0.02 * sequences.front().size();
		for (uint32_t j = 1; j < sequences.size(); ++j) {
			
			cout << "\t" << j << endl;

			uint32_t i = rank[j];

			// spoa::Alignment alignment;
			std::vector<res_t> res;
			int sink_id = -1;
			if (positions_[i].first < offset && positions_[i].second > sequences.front().size() - offset) {


				cout << "\taa: " << endl;
				cout << "\t" << sequences[i] << endl;

				res = poa(para, DAG, 0, 1, i, sequences[i].c_str(), sequences[i].size(), mpool, para->ab_band);
				cout << "\tbb" << endl;
				// alignment = alignment_engine->align(sequences_[i].first,
				//     sequences_[i].second, graph);
				sink_id = 1;
			}
			else {
				
				cout << "\tcc " << positions_[i].first << " " << positions_[i].second << endl;
				int beg_id = positions_[i].first + 2;
				int end_id = positions_[i].second + 1 == sequences.front().size() ? 1 : positions_[i].second + 1 - 2;
				//int end_id = positions_[i].second + 1 == sequences.front().size() ? 1 : positions_[i].second;

				end_id = min((int)end_id, (int) (((long)DAG->node.size()) -1));

				
				cout << beg_id << " " << end_id << endl;
				cout << sequences[i] << endl;
				//cout << string(sequences_[i].first + 1, sequences_[i].second - 1) << endl;

				res = poa(para, DAG, beg_id, end_id, i, sequences[i].c_str() + 1, sequences[i].size() - 1, mpool, para->ab_band);
				if(end_id == 1) sink_id = 1;
				// std::vector<int32_t> mapping;
				// auto subgraph = graph->subgraph(positions_[i].first,
				//     positions_[i].second, mapping);
				// alignment = alignment_engine->align(sequences_[i].first,
				//     sequences_[i].second, subgraph);
				// subgraph->update_alignment(alignment, mapping);
				cout << "\tdd" << endl;
				
			}
			
			cout << "\tee" << endl;
			// add res to graph
			DAG->add_path(para->m, i, res, sink_id);
			cout << "\tff" << endl;
			DAG->topsort(para, 0);
			cout << "\tgg" << endl;
			// if (qualities_[i].first == nullptr) {
			//     graph->add_alignment(alignment, sequences_[i].first,
			//         sequences_[i].second);
			// }
			// else {
			//     graph->add_alignment(alignment, sequences_[i].first,
			//         sequences_[i].second, qualities_[i].first,
			//         qualities_[i].second);
			// }
		}

		cout << "f" << endl;
		DAG->build_consensus(true);
		string consensus_ = DAG->cons;
		std::vector<int> coverages = DAG->coverages;

		//if (type_ == WindowType::kTGS && trim) {
		if (trim) {
			uint32_t average_coverage = (sequences.size() - 1) / 2;

			int32_t begin = 0, end = consensus_.size() - 1;
			for (; begin < static_cast<int32_t>(consensus_.size()); ++begin) {
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
			}
			else {
				consensus_ = consensus_.substr(begin, end - begin + 1);
			}
		}

		
		cout << "g" << endl;
		delete mpool;
		mpool = nullptr;  // 防止后续误用
		delete DAG;
		DAG = nullptr;  // 防止后续误用
		return consensus_;
	}
	*/

	string computeConsensus(const size_t contigIndexLocal, const size_t windowIndexLocal, const string contigOriginalSequence, auto& alignmentEngine){

		//cout << contigOriginalSequence << endl;
		//vector<string> sequenceStrs;
		vector<Window>& sequences = _contigWindowSequences[contigIndexLocal][windowIndexLocal];

		//cout << windowIndexLocal << " " << sequences.size() << endl;
		//u_int64_t nbCorrectedWindows = 0;
		//vector<DnaBitset2*> correctedWindows(windows.size());
	
	
		//cout << sequences.size() << endl;
			
		
		bool isLastWindow = (windowIndexLocal == _contigWindowSequences[contigIndexLocal].size()-1);


		if(sequences.size() < 2){

			//for(size_t i=0; i<sequences.size(); i++){ 
			//	delete sequences[i]._sequence;
			//}

			//addCorrectedWindow(false, DnaBitset2(contigOriginalSequence), contigIndexLocal, windowIndexLocal, outputContigFile, file_contigHeaders, pass);	
			//correctedWindows[w] = new DnaBitset2(contigOriginalSequence);
			//_logFile << "No sequences for window" << endl;
			return contigOriginalSequence;
		}
		


		std::sort(sequences.begin(), sequences.end(), [&](const Window& a, const Window& b) {
			if(a._posStart == b._posStart){
				return a.hash() < b.hash();
			}
			return a._posStart < b._posStart;
		});

		//cout << "a: " << contigOriginalSequence << endl;

		//std::unique_ptr<spoa::AlignmentEngine> alignmentEngine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -4);
		
		/*
		abpoa_t *ab = abpoa_init();


		string backboneQuality = "";
		for(size_t i=0; i<contigOriginalSequence.size(); i++){
			backboneQuality += '!';
		}

		sequences.insert(sequences.begin(), {DnaBitset2(contigOriginalSequence), backboneQuality, 0, contigOriginalSequence.size(), 1.0});
		

		int32_t n_seqs = sequences.size();
		// collect sequence length, trasform ACGT to 0123
		int *seq_lens = (int*)malloc(sizeof(int) * n_seqs);
		uint8_t **bseqs = (uint8_t**)malloc(sizeof(uint8_t*) * n_seqs);
		int **weights = (int**)malloc(sizeof(int*) * n_seqs);
		for (size_t i = 0; i < sequences.size(); ++i) {
			
			const string& seq = sequences[i]._sequence.to_string();
			const string& qual = sequences[i]._quality;

			seq_lens[i] = seq.size();
			bseqs[i] = (uint8_t*)malloc(sizeof(uint8_t) * seq_lens[i]);
			weights[i] = (int*)malloc(sizeof(int) * seq_lens[i]);
			for (size_t j = 0; j < seq_lens[i]; ++j) {
				bseqs[i][j] = nt4_table[(int)seq[j]];
				//weights[i][j] = 0;
				weights[i][j] = (int) qual[j];

				//if(i == 0){
				//	cout << "a: " << (char)(qual[j]) << " " << ((int) qual[j] )<< endl;
				//	//getchar();
				//}
				//else{
				//	cout << "b: " << (char)(qual[j]) << " " << ((int) qual[j] )<< endl;
				//	getchar();
				//}
				//if (j >= 12) weights[i][j] = 2;
				//else weights[i][j] = 0;
			}
		}
		
		// 1. directly output to stdout
		//fprintf(stdout, "=== output to stdout ===\n");
		// perform abpoa-msa
		// set weights as NULL if no quality score weights are used
		abpoa_msa(ab, abpt, n_seqs, NULL, seq_lens, bseqs, weights, NULL);


		//cout << (int)abpt->disable_seeding << endl;
		//cout << (int)abpt->progressive_poa << endl;
		//cout << (int)abpt->align_mode << endl;
		//cout << ((abpt->disable_seeding && abpt->progressive_poa==0) || abpt->align_mode != ABPOA_GLOBAL_MODE) << endl;

		string consensus = "";

		abpoa_cons_t *abc = ab->abc;

		//cout << abc->n_cons << endl;
		
		for (size_t i = 0; i < abc->n_cons; ++i) {

			for (size_t j = 0; j < abc->cons_len[i]; ++j){
				consensus += nt256_table[abc->cons_base[i][j]];
			}
		}
		
		//cout << "b: " << consensus << endl;
		
		string correctedSequence = consensus;
		
		// free seq-related variables
		for (size_t i = 0; i < n_seqs; ++i) { free(bseqs[i]); free(weights[i]); }
		free(bseqs); free(seq_lens); free(weights);

		// free abpoa-related variables
		abpoa_free(ab); 
		
		

		*/
		
		//return;

		//_logFile << sequences.size() << endl;
		//_logFile << "1" << endl;
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
			string dnaStr = window._sequence.to_string();

			//sequenceStrs.push_back(dnaStr);
			//if(windowIndexLocal == 24){
			//	cout << string(dnaStr) << endl;
			//}

			spoa::Alignment alignment;
			if (window._posStart < offset && window._posEnd > contigOriginalSequence.size() - offset) {
				alignment = alignmentEngine->Align(
					dnaStr.c_str(), dnaStr.size(),
					graph);
			} else {
				std::vector<const spoa::Graph::Node*> mapping;
				auto subgraph = graph.Subgraph(
					window._posStart,
					window._posEnd,
					&mapping);
				alignment = alignmentEngine->Align(
					dnaStr.c_str(), dnaStr.size(),
					subgraph);
				subgraph.UpdateAlignment(mapping, &alignment);
			}
			
			if (window._quality.size() == 0) {
				graph.AddAlignment(
					alignment,
					dnaStr.c_str(), dnaStr.size());
			} else {
				graph.AddAlignment(
					alignment,
					dnaStr.c_str(), dnaStr.size(),
					window._quality.c_str(), window._quality.size());
			}
			
			//free(dnaStr);
			//delete window._sequence;
		}


		vector<u_int32_t> coverages;
		string correctedSequence = graph.GenerateConsensus(&coverages);

		correctedSequence = trimConsensus(correctedSequence, coverages, sequences.size(), isLastWindow);
		

		//string correctedSequence = performCorrection2(sequences);
		//for(char letter : correctedSequence){
		//	checksum += letter;
		//}

		//ofstream file("/pasteur/appa/homes/gbenoit/appa/run/correction/test_deterministic/test_humanO1_4/loulou.fasta");
		//file << ">lala" << endl;
		//file << correctedSequence << endl;
		//file.close();
		//cout << windowIndexLocal << endl;
		//getchar();
		/*
		#pragma omp critical(debug)
		{
			
			
			if(correctedSequence.size() < 450){
				cout << endl;
				cout << "-----" << endl;
				cout << "Original sequence: " << contigOriginalSequence << endl;
				cout << correctedSequence.size() << " " << contigIndexLocal << " " << windowIndexLocal << endl;
				for(size_t i=0; i<sequenceStrs.size(); i++){
					cout << i << " " << sequenceStrs[i].size() << " " << sequenceStrs[i] << endl;
				}
				getchar();
			}

		}
		*/

		return correctedSequence;
	}

	string trimConsensus(const string& correctedSequence, const vector<u_int32_t>& coverages, const int& nbSequences, const bool& isLastWindow){

		string trimmedSequence = "";
		uint32_t average_coverage = nbSequences / 2;

		while(true){


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
				trimmedSequence = correctedSequence.substr(begin, end - begin + 1);
			}

			if(isLastWindow) break;
			if(trimmedSequence.size() > _windowLength*0.8) break;
			
			average_coverage += 1;

			if(average_coverage > nbSequences) return correctedSequence;
		}


		return trimmedSequence;
	}


	void addCorrectedWindow(bool success, const DnaBitset2& seq, size_t contigIndexLocal, size_t windowIndexLocal, BGZF* outputContigFile, ofstream& file_contigHeaders, const int& pass){

		
		#pragma omp critical(addCorrectedWindow)
		{
			
			//cout << "Add corrected window: " << windowIndexLocal << " " << seq->m_len << endl;

			_currentContigs[contigIndexLocal].push_back({windowIndexLocal, seq, success});

			if(_currentContigs[contigIndexLocal].size() == _contigWindowSequences[contigIndexLocal].size()){
				dumpCorrectedContig(contigIndexLocal, outputContigFile, file_contigHeaders, pass);
			}
		}
		
	}

	void dumpCorrectedContig(const u_int32_t& contigIndex, BGZF* outputContigFile, ofstream& file_contigHeaders, const int& pass){


		u_int64_t nbCorrectedWindows = 0;
		vector<CorrectedWindow>& correctedWindows = _currentContigs[contigIndex];

		for(const CorrectedWindow& correctedWindow : correctedWindows){
			if(correctedWindow._success){
				nbCorrectedWindows += 1;
			}
		}
		

		if(nbCorrectedWindows > 0){

			std::sort(correctedWindows.begin(), correctedWindows.end(), CorrectedWindowComparator);

			string contigSequence = "";
			//for(size_t w=0; w<correctedWindows.size(); w++){
			for(const CorrectedWindow& correctedWindow : correctedWindows){
				//if(correctedWindow._correctedSequence == nullptr) continue;
				
				//if(correctedWindows[w] == nullptr) continue;
				string seq = correctedWindow._correctedSequence.to_string();
				
				//if(_debug_contigName != "") cout << "Contig length: " << contigSequence.size() << "   added: " << string(seq).size() << endl;

				contigSequence += string(seq);
				//free(seq);

				//cout << contigSequence << endl;
				//cout << contigSequence.size() << endl;
				//getchar();
				//delete correctedWindows[w];
			}


			u_int64_t length = contigSequence.size();
			bool isValid = true;

			//cout << _minContigLength << " "<< _minContigCoverage << endl;
			if(_contigCoverages[contigIndex] <= _minContigCoverage){
				isValid = false;
			}
			else if(length < _minContigLength){
				isValid = false;
			}
			else if(length < 7500 && _contigCoverages[contigIndex] < 4){
				isValid = false;
			}

			//else if(length < 7500 && _contigCoverages[contigIndex] < 3){
			//	isValid = false;
			//}
			


			if(isValid){

				string header = _contigHeaders[contigIndex];
				
				if(_useMetamdbgHeaderStyle && pass > 0){
					
					float coverage = _contigCoverages[contigIndex];
					//string originalHeader = header;
					//header = Utils::split(header, '_')[0]; //Remove curator subpath suffix

					bool isCircular = false;
					//if(header.find("rc") != string::npos){
					//	string h = header.substr(0, header.size()-2);
					//	header = h + "_" + to_string(_contigCoverages[contigIndex]) + "x_rc";
					//}
					//else 
					if(header[header.size()-1] == 'c'){
						isCircular = true;
						//string h = header.substr(0, header.size()-1);
						//header = h + " length=" + to_string(contigSequence.size()) + " coverage=" + to_string(_contigCoverages[contigIndex]) + " circular=yes";
					}

					header.pop_back(); //remove circular indicator
					header.erase(0, 3); //remove "ctg"
					u_int32_t originalContigIndex = stoull(header);

					//else{
					//	string h = header.substr(0, header.size()-1);
					//	header = h + " length=" + to_string(contigSequence.size()) + " coverage=" + to_string(_contigCoverages[contigIndex]) + " circular=no";
					//}
					header = Utils::createContigHeader(originalContigIndex, contigSequence.size(), coverage, isCircular);
					
					//file_contigHeaders << originalHeader << "\t" << _outputContigIndex << endl;
				}
				

				u_int64_t checksum = 0;
				for(size_t i=0; i<contigSequence.size(); i++){
					checksum += (contigSequence[i]*contigSequence.size());
				}
				_checksum += checksum;
				_checksum_total += checksum;


				//cout << "Polish contig: " << contigSequence.size() << endl;
				header = ">" + header + '\n';// ">ctg" + to_string(contigIndex) + '\n';
				//header += '\n';
				bgzf_write(outputContigFile, (const char*)&header[0], header.size());
				contigSequence +=  '\n';
				bgzf_write(outputContigFile, (const char*)&contigSequence[0], contigSequence.size());
				//_logFile << _contigHeaders[contigIndex] << " " << contigSequence.size() << endl;
				


				//_outputContigIndex += 1;

			}



		}

		//for(CorrectedWindow& correctedWindow : correctedWindows){
		//	if(correctedWindow._correctedSequence == nullptr) continue;
		//	delete correctedWindow._correctedSequence;
		//}
		_currentContigs.erase(contigIndex);
	}


};	

#endif 


