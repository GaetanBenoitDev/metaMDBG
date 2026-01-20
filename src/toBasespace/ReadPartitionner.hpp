

#ifndef MDBG_METAG_READPARTITIONNER
#define MDBG_METAG_READPARTITIONNER

#include "../Commons.hpp"

class ReadPartitionner{
    
public:


	string _inputFilenameRead;
	string _inputFilenameContig;
	string _readPartitionDir;
	u_int32_t _nbPartitions;
	//phmap::parallel_flat_hash_map<ReadType, ReadMapping2>& _alignments;
	bool _useQual;
	long double _averageDistanceBetweenMinimizers;

	int _nbCores;

	class PartitionFile{

		public:

		gzFile _file;
		omp_lock_t _mutex;

		PartitionFile(u_int32_t partition, const string& partitionFilename){
			_file = gzopen(partitionFilename.c_str(), "wb");

			omp_init_lock(&_mutex);
		}

		~PartitionFile(){
			gzclose(_file);
			omp_destroy_lock(&_mutex);
		}
	};

	unordered_map<u_int32_t, u_int32_t> _contigToPartition;
	vector<PartitionFile*> _partitionFiles;
	vector<ofstream> _contigFiles;
	//unordered_map<u_int32_t, u_int64_t> _partitionNbReads;
	unordered_map<u_int32_t, float> _contigCoverages;

	struct AlignmentInfo{
		u_int32_t _contigIndex;
		bool _isReversed;
	};
	unordered_map<ReadType, AlignmentInfo> _readIndex_to_contigIndex;

	
	struct ContigStats{
		u_int32_t _contigIndex;
		u_int32_t _nbMinimizers;
	};

	vector<ContigStats> _contigStats;
	u_int64_t _checksumWrittenReads;

	ReadPartitionner(const string inputFilenameRead, const string inputFilenameContig, const string readPartitionDir, long double averageDistanceBetweenMinimizers, const int nbCores) {
		_inputFilenameRead = inputFilenameRead;
		_inputFilenameContig = inputFilenameContig;
		_readPartitionDir = readPartitionDir;
		_averageDistanceBetweenMinimizers = averageDistanceBetweenMinimizers;
		_nbCores = nbCores;
		_useQual = true;
	}

	void execute(){

		collectContigStats();
		computeContigCoverages();

		int nbPartitionsInit = max((u_int32_t)1, (u_int32_t)_nbCores);
		u_int64_t totalByteSize = 0;

		vector<u_int64_t> memoryPerPartitions(nbPartitionsInit, 0);

		u_int64_t maxMemory = 4000000000ull;

		for(const ContigStats& contigStat : _contigStats){

			int partitionIndex = 0;
			u_int64_t minMemory = -1;

			for(size_t i=0; i<memoryPerPartitions.size(); i++){
				if(memoryPerPartitions[i] < minMemory){
					minMemory = memoryPerPartitions[i];
					partitionIndex = i;
				}
			}

			u_int64_t currentParitionMemory = memoryPerPartitions[partitionIndex];
			u_int64_t contigMemory = computeContigMemory(contigStat._contigIndex, contigStat._nbMinimizers);

			//cout << "lala " << contigMemory << endl;
			
			if(currentParitionMemory > 0 && currentParitionMemory+contigMemory > maxMemory){
				memoryPerPartitions.push_back(0);
				partitionIndex = memoryPerPartitions.size()-1;
			}
			
			memoryPerPartitions[partitionIndex] += contigMemory;
			_contigToPartition[contigStat._contigIndex] = partitionIndex;
			

		}



		_nbPartitions = 0;
		for(size_t i=0; i<memoryPerPartitions.size(); i++){
			if(memoryPerPartitions[i] > 0){
				_nbPartitions += 1;
			}
		}
		
		Logger::get().debug() << "Nb partitions: " << _nbPartitions;


		Logger::get().debug() << "Write read partitions";
		writeReadPartitions();

		Logger::get().debug() << "Write contig partitions";
		writeContigPartitions();

		//_contigStats.clear();
		//_contigCoverages.clear();
		//cout << "_contigCoverages.clear(); a remettre" << endl;
	}


	void computeContigCoverages(){


		unordered_map<u_int32_t, vector<pair<u_int32_t, u_int32_t>>> contigHits;

		ifstream alignmentFile(_readPartitionDir + "/readsVsContigsAlignments.bin");

		while(true){

			ReadMapping2 alignment;
			bool isEof = alignment.read(alignmentFile);

			if(isEof) break;

			contigHits[alignment._contigIndex].push_back({alignment._contigStart, alignment._contigEnd});

		}

		alignmentFile.close();


		for(auto& it : contigHits){

			const u_int32_t contigIndex = it.first;

			//cout << it.first << " " << _contigLength[it.first] << endl;
			vector<u_int32_t> coverages(_contigStats[contigIndex]._nbMinimizers, 0);

			for(auto& interval : it.second){

				for(size_t i=interval.first; i<interval.second; i++){
					coverages[i] += 1;
				}
			}

			long double sum = 0;
			long double n = 0;
			for(size_t i=0; i<coverages.size(); i++){
				sum += coverages[i];
				n += 1;
			}
			
			float contigCoverage = sum / n;
			_contigCoverages[contigIndex] = contigCoverage;

		}

		contigHits.clear();
	}

	void writeReadPartitions(){

		for(u_int32_t i=0; i<_nbPartitions; i++){
			//_partitionNbReads[i] = 0;
			_partitionFiles.push_back(new PartitionFile(i, _readPartitionDir + "/" + to_string(i) + "_reads.gz"));
		}



		ifstream alignmentFile(_readPartitionDir + "/readsVsContigsAlignments.bin");

		while(true){

			ReadMapping2 alignment;
			bool isEof = alignment.read(alignmentFile);

			if(isEof) break;

			_readIndex_to_contigIndex[alignment._readIndex] = {alignment._contigIndex, alignment._isReversed};

		}

		alignmentFile.close();


		//_contigCoverages.clear();

		Logger::get().debug() << "\tNb partitions: " << _nbPartitions;
		
		_checksumWrittenReads = 0;
		ReadParserParallel readParser(_inputFilenameRead, false, false, _nbCores);
		readParser.parse(ReadPartitionFunctor(*this));

		for(PartitionFile* partitionFile : _partitionFiles){
			delete partitionFile;
		}
		_partitionFiles.clear();

	}


	void collectContigStats(){
		KminmerParserParallel parser(_inputFilenameContig, 0, 0, false, false, 1);
		parser.parseSequences(LoadContigsFunctor(*this));
	}

	class LoadContigsFunctor {

		public:

		ReadPartitionner& _parent;

		LoadContigsFunctor(ReadPartitionner& parent) : _parent(parent){
		}

		LoadContigsFunctor(const LoadContigsFunctor& copy) : _parent(copy._parent){
		}

		~LoadContigsFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {

			u_int32_t readIndex = kminmerList._readIndex;

			u_int32_t contigLength = kminmerList._readMinimizers.size();

			_parent._contigStats.push_back({readIndex, contigLength});
			//_parent._mContigs.push_back({readIndex, kminmerList._readMinimizers});
			
			//cout << readIndex << " " << kminmerList._readMinimizers.size() << endl;

		}
	};

	/*
	void collectContigStats(){
		auto fp = std::bind(&ReadPartitionner::collectContigStats_read, this, std::placeholders::_1);
		ReadParser readParser(_inputFilename_contigs, true, false);
		readParser.parse(fp);
	}
	
	void collectContigStats_read(const Read& read){
		_contigStats.push_back({(u_int32_t)read._index, (u_int32_t)read._seq.size()});
	}
	*/

	u_int64_t computeContigMemory(const u_int32_t contigIndex, const size_t nbMinimizers){
		
		u_int32_t contigLength = nbMinimizers * _averageDistanceBetweenMinimizers;
		u_int64_t windowByteSize = 0;
		u_int32_t contigCoverage = _contigCoverages[contigIndex];
		contigCoverage = max((u_int32_t)1, contigCoverage);

		//cout << contigCoverage << endl;
		//if(_useQual){
		//	windowByteSize = contigCoverage  * ((_windowLength/4) + _windowLength);
		//}
		//else{
		//	windowByteSize = contigCoverage  * (_windowLength/4);
		//}

		//size_t nbWindows = ceil((double)contigLength / (double)_windowLength);
		//return nbWindows * windowByteSize;

		if(_useQual){
			return ceil(contigCoverage * (contigLength + contigLength/4.0));
		}
		
		return ceil(contigCoverage * (contigLength/4.0));
	}


	class ReadPartitionFunctor {

		public:

		ReadPartitionner& _parent;


		ReadPartitionFunctor(ReadPartitionner& parent) : _parent(parent){
		}

		ReadPartitionFunctor(const ReadPartitionFunctor& copy) : _parent(copy._parent){
		}

		~ReadPartitionFunctor(){
		}

		void operator () (const Read& read) {


			ReadType readIndex = read._index;
			//cout << readIndex << endl;
			//const string& readName = Utils::shortenHeader(read._header);

			if(read._index % 100000 == 0) Logger::get().debug() << "\t" << read._index;
			
			if(_parent._readIndex_to_contigIndex.find(readIndex) == _parent._readIndex_to_contigIndex.end()) return;

			const AlignmentInfo& al = _parent._readIndex_to_contigIndex[readIndex];
			//unordered_set<u_int32_t> writtenPartitions;

			//for(const Alignment& al : _contigPolisher._alignments[readIndex]){
			const u_int32_t contigIndex = al._contigIndex;
			//if(contigName != "ctg90297c") return;
			//u_int32_t contigIndex = al._contigIndex; //_contigPolisher._alignments[readIndex]._contigIndex;
			//_logFile << contigIndex << " " << (_contigPolisher._contigToPartition.find(contigIndex) != _contigPolisher._contigToPartition.end()) << endl;
			if(_parent._contigToPartition.find(contigIndex) == _parent._contigToPartition.end()) return;

			const u_int32_t partition = _parent._contigToPartition[contigIndex];

			//if(writtenPartitions.find(partition) != writtenPartitions.end()) return;
			//writtenPartitions.insert(partition);

			//_logFile << partition << endl;
			PartitionFile* partitionFile = _parent._partitionFiles[partition];
			bool isFastq = read._qual.size() > 0 && _parent._useQual;

			string readSeq = read._seq;
			string qualSeq = read._qual;
			if(al._isReversed){
				Utils::toReverseComplement(readSeq);
				if(isFastq) std::reverse(qualSeq.begin(), qualSeq.end());
			}

			omp_set_lock(&partitionFile->_mutex);
			
			//if(partition == 0){
				//_logFile << "Read: " << readIndex << endl;
			//}


			//_parent._partitionNbReads[partition] += 1;

			u_int32_t readSize = readSeq.size();

			if(isFastq){
				//string header = '@' + to_string(read._index) + '\n';
				string header = '@' + to_string(readIndex) + '\n';
				string seq = readSeq + '\n';
				gzwrite(partitionFile->_file, (const char*)&header[0], header.size());
				gzwrite(partitionFile->_file, (const char*)&seq[0], seq.size());
				static string strPlus = "+\n";
				gzwrite(partitionFile->_file, (const char*)&strPlus[0], strPlus.size());
				string qual = qualSeq + '\n';
				gzwrite(partitionFile->_file, (const char*)&qual[0], qual.size());
			}
			else{
				//string header = '>' + to_string(read._index) + '\n';
				string header = '>' + to_string(readIndex) + '\n';
				string seq = readSeq + '\n';
				gzwrite(partitionFile->_file, (const char*)&header[0], header.size());
				gzwrite(partitionFile->_file, (const char*)&seq[0], seq.size());
			}

			//for(char letter : readName){
			//	_contigPolisher._checksumWrittenReads += letter;
			//}
			//_logFile << readName << " " << _contigPolisher._checksumWrittenReads << endl;

			//gzwrite(partitionFile->_file, (const char*)&readIndex, sizeof(readIndex));
			//gzwrite(partitionFile->_file, (const char*)&readSize, sizeof(readSize));
			//gzwrite(partitionFile->_file, read._seq.c_str(), read._seq.size());
			//gzwrite(partitionFile->_file, read._qual.c_str(), read._qual.size());
			omp_unset_lock(&partitionFile->_mutex);
			//}
		}
	};



	void writeContigPartitions(){

		for(u_int32_t i=0; i<_nbPartitions; i++){
			_contigFiles.push_back(ofstream(_readPartitionDir + "/" + to_string(i) + "_contigs.bin"));
		}

		//ReadParserParallel readParser(_inputFilenameContig, true, false, _nbCores);
		//readParser.parse(ContigPartitionFunctor(*this));

		KminmerParserParallel parser(_inputFilenameContig, 0, 0, false, false, 1);
		//parser._hasContigIndex = true;
		parser.parseSequences(ContigPartitionFunctor(*this));
		
		for(ofstream& file : _contigFiles){
			file.close();
		}

	}


	class ContigPartitionFunctor {

		public:

		ReadPartitionner& _parent;

		ContigPartitionFunctor(ReadPartitionner& parent) : _parent(parent){
		}

		ContigPartitionFunctor(const ContigPartitionFunctor& copy) : _parent(copy._parent){
		}

		~ContigPartitionFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {

			const u_int32_t contigIndex = kminmerList._readIndex;
			//cout << contigIndex << endl;
			const u_int8_t isCircular = kminmerList._isCircular;

			if(_parent._contigToPartition.find(contigIndex) == _parent._contigToPartition.end()) return;

			const u_int32_t partition = _parent._contigToPartition[contigIndex];

			ofstream& contigFile = _parent._contigFiles[partition];

			u_int32_t contigSize = kminmerList._readMinimizers.size();
			contigFile.write((const char*)&contigSize, sizeof(contigSize));
			contigFile.write((const char*)&isCircular, sizeof(isCircular));
			contigFile.write((const char*)&kminmerList._readMinimizers[0], contigSize*sizeof(MinimizerType));
			contigFile.write((const char*)&contigIndex, sizeof(contigIndex));


		}
	};


};	

#endif 


