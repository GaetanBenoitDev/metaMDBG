
/*

//"reprise: extraire un meilleur gfa qui ne contient que les nodes de la ref + du voisinage"

./bin/metaMDBG extend ~/appa/run/correction/nanopore_AD/asm/tmp/ --ref ~/appa/data/nanopore/contigLink/metaMDBG_bins/bin_166.fasta -t 8

./bin/metaMDBG toMinspace /pasteur/appa/homes/gbenoit/appa/run/correction/nanopore_AD/asm/tmp/ /pasteur/appa/homes/gbenoit/appa/run/correction/nanopore_AD/asm/tmp/_circ_pass/3.138429//assembly_graph.gfa.unitigs.nodepath /pasteur/appa/homes/gbenoit/appa/run/correction/nanopore_AD/asm/tmp/_circ_pass/3.138429//assembly_graph.gfa.unitigs -t 16
./bin/metaMDBG map ~/appa/run/correction/nanopore_AD/asm/tmp/_circ_pass/3.138429/ ~/appa/data/nanopore/metaflye/AD_nanopore_nocor/binning_c/checkm/__checkm/bins_/complete/bin.166.fa
*/

#ifndef MDBG_METAG_EXTEND
#define MDBG_METAG_EXTEND


#include "Commons.hpp"
#include "graph/GfaParser.hpp"
#include "graph/GraphSimplify.hpp"
#include "../utils/spoa64/include/spoa64/spoa.hpp"
 

#define PRINT_DEBUG_REPEAT




class Circulizer : public Tool{

public:

	string _inputDir;
	int _nbCores;
	string _referenceFilename;
	MDBG* _mdbg;

	float _minimizerDensity;
    size_t _minimizerSize;
    size_t _kminmerSize;
    size_t _kminmerSizeFirst;
	float _minimizerSpacingMean;
	float _kminmerLengthMean;
	float _kminmerOverlapMean;
    size_t _kminmerSizePrev;
	size_t _kminmerSizeLast;
	size_t _meanReadLength;
	u_int16_t _nbDatasets;

	string _tmpDir;


	Circulizer(): Tool (){

	}

	~Circulizer(){

	}


	void parseArgs(int argc, char* argv[]){

		args::ArgumentParser parser("circ", ""); //"This is a test program.", "This goes after the options."
		args::Positional<std::string> arg_outputDir(parser, "outputDir", "", args::Options::Required);
		args::ValueFlag<string> arg_referenceFilename(parser, "", "reference filename", {"ref"}, "");
		args::ValueFlag<int> arg_nbCores(parser, "", "Number of cores", {ARG_NB_CORES2}, NB_CORES_DEFAULT_INT);
		args::Flag arg_help(parser, "", "", {'h', "help"}, args::Options::Hidden);


		try
		{
			parser.ParseCLI(argc, argv);
		}
		catch (const std::exception& e)
		{
			cerr << parser;
			cerr << e.what() << endl;
			exit(0);
		}

		if(arg_help){
			cerr << parser;
			exit(0);
		}

		_inputDir = args::get(arg_outputDir);

		_nbCores = args::get(arg_nbCores);

		//_filename_contigs = _inputDir + "/contig_data.txt";

		_referenceFilename = "";
		if(arg_referenceFilename){
			_referenceFilename = args::get(arg_referenceFilename);
		}

		string filename_parameters = _inputDir + "/parameters.gz";
		gzFile file_parameters = gzopen(filename_parameters.c_str(),"rb");
		gzread(file_parameters, (char*)&_minimizerSize, sizeof(_minimizerSize));
		gzread(file_parameters, (char*)&_kminmerSize, sizeof(_kminmerSize));
		gzread(file_parameters, (char*)&_minimizerDensity, sizeof(_minimizerDensity));
		gzread(file_parameters, (char*)&_kminmerSizeFirst, sizeof(_kminmerSizeFirst));
		gzread(file_parameters, (char*)&_minimizerSpacingMean, sizeof(_minimizerSpacingMean));
		gzread(file_parameters, (char*)&_kminmerLengthMean, sizeof(_kminmerLengthMean));
		gzread(file_parameters, (char*)&_kminmerOverlapMean, sizeof(_kminmerOverlapMean));
		gzread(file_parameters, (char*)&_kminmerSizePrev, sizeof(_kminmerSizePrev));
		gzread(file_parameters, (char*)&_kminmerSizeLast, sizeof(_kminmerSizeLast));
		gzread(file_parameters, (char*)&_meanReadLength, sizeof(_meanReadLength));
		gzread(file_parameters, (char*)&_nbDatasets, sizeof(_nbDatasets));
		gzclose(file_parameters);

		openLogFile(_inputDir);

		//_kminmerSize = _kminmerSizeFirst;

		_logFile << endl;
		_logFile << "Input dir: " << _inputDir << endl;
		//cout << "Output filename: " << _outputFilename << endl;
		_logFile << "Minimizer length: " << _minimizerSize << endl;
		_logFile << "Kminmer length: " << _kminmerSize << endl;
		_logFile << "Density: " << _minimizerDensity << endl;
		_logFile << endl;
	}


	std::unique_ptr<spoa64::AlignmentEngine> _alignmentEngine;

	void execute (){

		_minPoaNodeCoverage = 2;
		_print_debug = false;
		_alignmentEngine = spoa64::AlignmentEngine::Create(spoa64::AlignmentType::kNW, 3, -5, -4);

		_tmpDir = _inputDir + "/_circ_pass/";

		loadReferences();
		loadGraph();

		closeLogFile();
	}

	unordered_map<u_int32_t, vector<u_int32_t>> _ref_nodeName_to_unitigIndexes;
	struct ReferenceContig{
		vector<u_int32_t> _nodePath;
	};

	unordered_map<u_int32_t, ReferenceContig> _referenceContigs;

	ofstream _referencePathFile;
	ofstream _referenceColorFile;
	ofstream _referenceSelectedColorFile;
    vector<u_int32_t> _nodeName_to_abundance;

	void loadReferences(){

		if(_referenceFilename.empty()) return;

		cout << "Loading reference" << endl;

		string mdbg_filename = _inputDir + "/kminmerData_min.txt";
		_mdbg = new MDBG(_kminmerSize);
		_mdbg->load(mdbg_filename, false);


		auto fp = std::bind(&Circulizer::loadReferences_read_index, this, std::placeholders::_1);
		ReadParser readParser(_referenceFilename, true, false, _logFile);
		readParser.parse(fp);

		cout << "Loading reference minimizer sequences" << endl;
		getRealReferenceSequence();

		//auto fp2 = std::bind(&BinnerHic::loadReferences_read, this, std::placeholders::_1);
		//ReadParser readParser2(_referenceFilename, true, false, _logFile);
		//readParser2.parse(fp2);

		delete _mdbg;
		_ref_nodeName_to_unitigIndexes.clear();
	}

	unordered_map<u_int32_t, u_int32_t> _utgIndex_to_referenceIndex;
	unordered_map<u_int32_t, u_int32_t> _referenceIndex_to_referenceLength;
	unordered_map<u_int32_t, string> _referenceIndex_to_referenceName;

	void loadReferences_read_index(const Read& read){

		vector<string>* fields = new vector<string>();
		GfaParser::tokenize(read._header, fields, '_');

		//for(const string& field: (*fields)){
		//	cout << field << endl;
		//}

		string index = (*fields)[0];
		index.erase(0, 3);

		u_int32_t refIndex = stoull(index);
		//cout << (*fields)[0] << " " << refIndex << endl;
		delete fields;

		//_referenceContigIndexes.push_back(refIndex);
		_utgIndex_to_referenceIndex[refIndex] = read._index;
		_referenceIndex_to_referenceLength[read._index] = read._seq.size();
		_referenceIndex_to_referenceName[read._index] = read._header;

	}
	
	void getRealReferenceSequence(){

		cerr << "Get real reference sequences" << endl;
		KminmerParserParallel parser(_inputDir + "/contig_data.txt", _minimizerSize, _kminmerSize, false, false, 1);
		parser.parse(SplitContigsFunctor(*this));

	}


	class SplitContigsFunctor {

		public:

		Circulizer& _parent;

		SplitContigsFunctor(Circulizer& parent) : _parent(parent){
		}

		SplitContigsFunctor(const SplitContigsFunctor& copy) : _parent(copy._parent){
		}

		~SplitContigsFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {

			if(_parent._utgIndex_to_referenceIndex.find(kminmerList._readIndex) == _parent._utgIndex_to_referenceIndex.end()) return;

			u_int32_t refIndex = _parent._utgIndex_to_referenceIndex[kminmerList._readIndex];
			//cout << _parent._utgIndex_to_referenceIndex[kminmerList._readIndex] << " " << kminmerList._kminmersInfo.size()*270 << " " << _parent._referenceIndex_to_referenceLength[_parent._utgIndex_to_referenceIndex[kminmerList._readIndex]] << endl;
			
			vector<u_int32_t> nodePath;
			for(const ReadKminmerComplete& kminmerInfo : kminmerList._kminmersInfo){

				const KmerVec& vec = kminmerInfo._vec;

				if(_parent._mdbg->_dbg_nodes.find(vec) == _parent._mdbg->_dbg_nodes.end()){
					continue;
				}


				u_int32_t nodeName = _parent._mdbg->_dbg_nodes[vec]._index;
				u_int32_t nodeIndex = BiGraph::nodeName_to_nodeIndex(nodeName, !kminmerInfo._isReversed);

				nodePath.push_back(nodeIndex);
			}

			_parent._referenceContigs[refIndex] = {nodePath};
		}
	};
	

	unordered_set<u_int32_t> _debug_referenceUnitigIndexes;
	vector<vector<u_int32_t>> _debug_referenceUnitigPaths;


	void computeReferencePath(){
		
		_debug_referenceUnitigPaths.clear();

		_referencePathFile = ofstream(_passDir + "/refPath.csv");
		_referencePathFile << "Name,Pos" << endl;

		_referenceColorFile = ofstream(_passDir + "/refColor.csv");
		_referenceColorFile << "Name,Color" << endl;

		ofstream referenceNameFile(_passDir + "/refName.csv");
		referenceNameFile << "Name,Text" << endl;

		_referenceSelectedColorFile = ofstream(_passDir + "/refSelected.csv");
		_referenceSelectedColorFile << "Name,Color" << endl;

		ofstream clusterFile(_passDir + "/metatorColor.csv");
		clusterFile << "Name,Color" << endl;

		u_int32_t nbMissingGenomicKminmers = 0;

		unordered_set<u_int32_t> unitigIndexes;


		for(const auto& it : _referenceContigs){


			u_int32_t refIndex = it.first;

			//if(refIndex != 3) continue;

			const ReferenceContig& refContig = it.second;

			vector<u_int32_t> unitigPath;
			unordered_map<u_int32_t, vector<u_int32_t>> unitigPos;
			u_int32_t pos = 0;
			u_int32_t currentUnitig = -1;
			u_int32_t currentUnitigPath = -1;

			for(u_int32_t nodeIndex : refContig._nodePath){

				if(_nodeIndex_to_unitigIndex.find(nodeIndex) == _nodeIndex_to_unitigIndex.end()){


					if(currentUnitigPath != -1){
						currentUnitigPath = -1;
						unitigPath.push_back(currentUnitigPath);
					}


					//cout << "missing reference kminmer" << endl;
					nbMissingGenomicKminmers += 1;
					continue;
				}

				u_int32_t unitigIndex = _nodeIndex_to_unitigIndex[nodeIndex];

				//_progressiveAbundanceFilter->_unitigGraph->_nodes[unitigIndex]->_readIndexes.push_back(refIndex);
				//_progressiveAbundanceFilter->_unitigGraph->_nodes[_progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(unitigIndex)]->_readIndexes.push_back(refIndex);



				if(unitigIndex != currentUnitigPath){
					currentUnitigPath = unitigIndex;
					unitigPath.push_back(unitigIndex);
				}

				if(unitigIndex != currentUnitig){
					currentUnitig = unitigIndex;
					unitigPos[unitigIndex].push_back(pos);
					//cout << "utg" << unitigIndex << " " << pos << endl;
					pos += 1;
					//getchar();

					_debug_referenceUnitigIndexes.insert(unitigIndex);
					_debug_referenceUnitigIndexes.insert(_progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(unitigIndex));
				}

				unitigIndexes.insert(unitigIndex);
			}

			_debug_referenceUnitigPaths.push_back(unitigPath);

			for(auto& it : unitigPos){

				string posStr = "";
				for(u_int32_t pos : it.second){
					posStr += to_string(pos) + "-";
				}
				posStr.pop_back();

				_referencePathFile << "utg" << it.first << "," << posStr << endl;
				_referencePathFile << "utg" << _progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(it.first) << "," << posStr << endl;
				
				referenceNameFile << "utg" << it.first << "," << _referenceIndex_to_referenceName[refIndex] << endl;
				referenceNameFile << "utg" << _progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(it.first) << "," << _referenceIndex_to_referenceName[refIndex] << endl;
				
				_referenceColorFile << "utg" << it.first << "," << _referenceIndex_to_referenceName[refIndex] << endl;
				_referenceColorFile << "utg" << _progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(it.first) << "," << _referenceIndex_to_referenceName[refIndex] << endl;
				
				_referenceSelectedColorFile << "utg" << it.first << "," << "red" << endl;
				_referenceSelectedColorFile << "utg" << _progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(it.first) << "," << "red" << endl;
			
			
			}

		}
		
		_referencePathFile.close();
		_referenceColorFile.close();
		_referenceSelectedColorFile.close();
		clusterFile.close();

		ofstream successorPathFile(_passDir + "/successorUnitigPath.csv");
		successorPathFile << "Name,Text" << endl;

		for(u_int32_t unitigIndex : unitigIndexes){
			
			vector<u_int32_t> path = getMostSupportedPath(_unitigGraph->_nodes[unitigIndex]);

			string s = "";
			for(u_int32_t unitigIndex : path){
				s += "utg" + to_string(unitigIndex) + " ";
			}

			vector<u_int32_t> pathRev = getMostSupportedPath(_unitigGraph->unitigIndex_toReverseDirection(_unitigGraph->_nodes[unitigIndex]));

			string s_rev = "";
			for(u_int32_t unitigIndex : pathRev){
				s_rev += "utg" + to_string(unitigIndex) + " ";
			}

			if(unitigIndex % 2 == 0){
				successorPathFile << "utg" << unitigIndex << "," << s << " ||| " << s_rev << endl;
			}
			else{
				successorPathFile << "utg" << (unitigIndex-1) << "," << s << " ||| " << s_rev << endl;
			}

		}

		successorPathFile.close();

		//cout << "Nb missing genomic kminmers: " << nbMissingGenomicKminmers << endl;
		//if(nbMissingGenomicKminmers > 0){
		//	getchar();
		//}
	}


	void loadGraph(){


		string gfaFilename = _inputDir + "/minimizer_graph.gfa";

		cout << "Loading gfa" << endl;
		GraphSimplify* graphSimplify;// = new GraphSimplify(gfaFilename, _inputDir, 0, _kminmerSize, _nbCores, _kminmerLengthMean, _kminmerOverlapMean, _logFile, _nodeName_to_abundance);
		//!

		ifstream kminmerFile(_inputDir + "/kminmerData_min.txt");
        _nodeName_to_abundance.resize(graphSimplify->_graphSuccessors->_nbNodes/2, 0);

		u_int64_t abSum = 0;
		u_int64_t qualSum = 0;

		while (true) {

			vector<MinimizerType> minimizerSeq;
			minimizerSeq.resize(_kminmerSize);
			kminmerFile.read((char*)&minimizerSeq[0], minimizerSeq.size()*sizeof(MinimizerType));

			if(kminmerFile.eof()) break;

			u_int32_t nodeName;
			u_int32_t abundance;
			//u_int32_t quality;
			//bool isReversed = false;

			kminmerFile.read((char*)&nodeName, sizeof(nodeName));
			kminmerFile.read((char*)&abundance, sizeof(abundance));
			//kminmerFile.read((char*)&quality, sizeof(quality));

			abSum += abundance;
			//qualSum += quality;

			//cout << nodeName << " " << abSum << endl;

			//cout << nodeName << " " << abundance << " " << quality << endl;
			//if(quality == 0){
			//	cout << quality << endl;
			//	getchar();
			//}
			//if(abundance == 1) continue;
			//KmerVec vec;
			//vec._kmers = minimizerSeq;

			_nodeName_to_abundance[nodeName] = abundance;
			//_graph->_graphSuccessors->_nodeLengths[nodeName] = _kminmerLengthMean;
			//_kminmerAbundances[vec] = abundance;
		}

		//cout << abSum << " " << qualSum << endl;
		//getchar();

		kminmerFile.close();



		//if(_kminmerSize > 4){
		//GfaParser::binaryGraph_to_gfa(_gfaFilename, _kminmerLengthMean, _kminmerOverlapMean, _gfaFilename+".gfa", _graph->_graphSuccessors->_nodeDatas);
		//if(_kminmerSize > 45){
		//	getchar();
		//}
		//	cout << "lala" << endl;
		//	getchar();
		//}

		cout << "Compating graph" << endl;

		vector<UnitigData> unitigDatas;
		graphSimplify->debug_writeGfaErrorfree(2000, 2000, -1, _kminmerSize, false, true, false, unitigDatas, true, false, false, false, false, true, _mdbg, _minimizerSize, 1, false, false);
		
		//collectPalindrome();
		
		//indexReads();



		cout << "Filering graph" << endl;
		_areReadsIndexed = false;
		//_unitigGraph = graphSimplify->_unitigGraph;
		//!

		GraphCleanedFunctor functor(*this);

		cout << "Superbubble disabled to speed up" << endl;
        //_progressiveAbundanceFilter = new ProgressiveAbundanceFilter(graphSimplify->_unitigGraph, _inputDir, _kminmerSize, true, _logFile);
		//!

		delete graphSimplify;

        _progressiveAbundanceFilter->execute(functor);

		
	}


	ProgressiveAbundanceFilter* _progressiveAbundanceFilter;
	UnitigGraph* _unitigGraph;
	string _passDir;
	bool _areReadsIndexed;
	unordered_map<u_int32_t, u_int32_t> _nodeIndex_to_unitigIndex;

	class GraphCleanedFunctor {

		public:

		Circulizer& _parent;
		bool _isDone;

		GraphCleanedFunctor(Circulizer& parent) : _parent(parent){
			_isDone = false;
			//_lastCutoff = -1;
		}

		GraphCleanedFunctor(const GraphCleanedFunctor& copy) : _parent(copy._parent){
			_isDone = false;
		}
		

		~GraphCleanedFunctor(){
		}



		void operator () (float cutoff) {

            float nextCutoff = cutoff;
			nextCutoff = nextCutoff * (1+0.1);

			//float minUnitigAbundance = cutoff / 0.5;
			float maxUnitigAbundance = cutoff / 0.5;
			cout << cutoff << " " << nextCutoff << " " << maxUnitigAbundance << endl;


			//if(cutoff < 5){
				//cout << "skipping cutoff " << cutoff << endl;
				//return;
			//}




			_isDone = true;

			_parent._passDir = _parent._tmpDir + "/" + to_string(cutoff) + "/";

			if(!fs::exists (_parent._passDir)){
				fs::create_directories(_parent._passDir); 
			}
			
			_parent._unitigGraph->save(_parent._passDir + "/assembly_graph.gfa", 0);

			cout << "Indexing unitigs" << endl;
			_parent.indexUnitigs();
						
			if(!_parent._areReadsIndexed){
				_parent._areReadsIndexed = true;
				_parent.indexReads(true, true);
			}
			
			cout << "Compute reference path" << endl;
			_parent.computeReferencePath();

			
			_parent._isUnitigVisited.clear();
			_parent._isUnitigVisitedUnique.clear();

			_parent._isUnitigVisited.resize(_parent._unitigGraph->_nodes.size(), false);
			_parent._isUnitigVisitedUnique.resize(_parent._unitigGraph->_nodes.size(), false);






			vector<u_int32_t> unitigPath;
			//bool isCircular = _parent.extendPath(_parent._unitigGraph->_nodes[1802], unitigPath);

			//"reprise: poa bubuluuuu detection"
			
			u_int32_t unitigIndex = 772;
			vector<u_int32_t> path = _parent.getMostSupportedPath(_parent._unitigGraph->_nodes[unitigIndex]);
			
			cout << endl;
			for(const vector<u_int32_t>& unitigPath : _parent._unitigGraph->_nodes[unitigIndex]->_successorUnitigPaths){
				for(u_int32_t unitigIndex : unitigPath){
					cout << "utg" << unitigIndex << " ";
				}
				cout << endl;
			}

			cout << "Most supported path:" << endl;
			for(u_int32_t unitigIndex: path){
				cout << "\tutg" << unitigIndex << endl;
			}
			

			exit(1);
		}





	};

	void indexUnitigs(){

		_nodeIndex_to_unitigIndex.clear();

		for(UnitigGraph::Node* node : _unitigGraph->_nodes){

			if(node->_unitigIndex == -1) continue;

			for(u_int32_t nodeIndex : node->_nodes){
				_nodeIndex_to_unitigIndex[nodeIndex] = node->_unitigIndex;
			}

		}
	}


	void indexReads(bool indexContigs, bool indexUnitigs){

		cout << "Indexing read" << endl;

		for(UnitigGraph::Node* node : _progressiveAbundanceFilter->_unitigGraph->_nodes){

			if(node->_unitigIndex == -1) continue;

			node->_readIndexes.clear();
			node->_successorUnitigPaths.clear();
		}

		
		cout << "Load node names" << endl;
		string mdbg_filename = _inputDir + "/kminmerData_min.txt";
		_mdbg = new MDBG(_kminmerSize);
		_mdbg->load(mdbg_filename, false);



		cout << "Indexing reads" << endl;

		KminmerParserParallel parser(_inputDir + "/read_data_corrected.txt", _minimizerSize, _kminmerSize, false, true, _nbCores);
		parser.parse(IndexReadsFunctor(*this, indexContigs, indexUnitigs));


		delete _mdbg;


	}


	class IndexReadsFunctor {

		public:

		Circulizer& _parent;
		bool _indexContigs;
		bool _indexUnitigs;

		IndexReadsFunctor(Circulizer& parent, bool indexContigs, bool indexUnitigs) : _parent(parent){
			_indexContigs = indexContigs;
			_indexUnitigs = indexUnitigs;
		}

		IndexReadsFunctor(const IndexReadsFunctor& copy) : _parent(copy._parent){
			_indexContigs = copy._indexContigs;
			_indexUnitigs = copy._indexUnitigs;
		}

		~IndexReadsFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {

			u_int32_t readIndex = kminmerList._readIndex;
			if(readIndex % 100000 == 0) cout << "Indexing reads: " << readIndex << endl;

			vector<u_int32_t> readNodeIndex;

			for(u_int16_t i=0; i<kminmerList._kminmersInfo.size(); i++){
			
				
				const ReadKminmerComplete& kminmerInfo = kminmerList._kminmersInfo[i];

				KmerVec vec = kminmerInfo._vec;


				if(_parent._mdbg->_dbg_nodes.find(vec) == _parent._mdbg->_dbg_nodes.end()){
					//cout << "XXX1" << " ";
				}


				/*
				if(vec.isPalindrome()){
					cout << endl;
					cout << "palindrome" << endl;
					for(u_int32_t m : vec._kmers){
						cout << m << " ";
					}
					cout << endl;
				}
				*/

				if(_parent._mdbg->_dbg_nodes.find(vec) == _parent._mdbg->_dbg_nodes.end()){
					continue;
				}


				u_int32_t nodeName = _parent._mdbg->_dbg_nodes[vec]._index;
				u_int32_t nodeIndex = BiGraph::nodeName_to_nodeIndex(nodeName, !kminmerInfo._isReversed);

				readNodeIndex.push_back(nodeIndex);

				if(_parent._nodeIndex_to_unitigIndex.find(nodeIndex) == _parent._nodeIndex_to_unitigIndex.end()){
					//cout << "XXX2" << " ";
				}
				else{
					//cout << "OOOO" << " ";
				}

			}

			/*
			#pragma omp critical(indexReadsFunctor)
			{
				bool isHere1 = false;
				bool isHere2 = false;

				for(u_int32_t nodeIndex : readNodeIndex){

					//cout << nodeIndex << " " << _parent._nodeIndex_to_unitigIndex[nodeIndex] << endl;

					if((_parent._nodeIndex_to_unitigIndex[nodeIndex] == 3140842 || _parent._nodeIndex_to_unitigIndex[nodeIndex] == 3140843)){
						isHere1 = true;
						//cout << "lol" << endl;
						//getchar();
					}

					if((_parent._nodeIndex_to_unitigIndex[nodeIndex] == 4067790 || _parent._nodeIndex_to_unitigIndex[nodeIndex] == 4067791)){

						isHere2 = true;
					}

				}

				if(isHere1 && isHere2){
					
					cout << endl;
					for(u_int32_t nodeIndex : readNodeIndex){
						cout << nodeIndex << " " << _parent._nodeIndex_to_unitigIndex[nodeIndex] << endl;
					}
					getchar();
				}
			}
			*/

			processKminmers(readIndex, readNodeIndex, kminmerList);

			std::reverse(readNodeIndex.begin(), readNodeIndex.end());
			for(size_t i=0; i<readNodeIndex.size(); i++){
				readNodeIndex[i] = UnitigGraph::nodeIndex_toReverseDirection(readNodeIndex[i]);
			}

			processKminmers(readIndex, readNodeIndex, kminmerList);


		}

		void processKminmers(u_int64_t readIndex, const vector<u_int32_t>& readNodeIndex, const KminmerList& kminmerList){


			vector<u_int32_t> unitigs;
			//vector<u_int32_t> contigs;

			u_int32_t prevUnitigIndex = -1;
			//u_int32_t prevContigIndex = -1;

			/*
			unordered_map<u_int32_t, u_int32_t> contigLastPositions;

			u_int32_t contigPos = 0;

			for(u_int32_t nodeIndex : readNodeIndex){
			
				if(_parent._nodeIndex_to_unitigIndex.find(nodeIndex) == _parent._nodeIndex_to_unitigIndex.end()) continue;
				
				u_int32_t unitigIndex = _parent._nodeIndex_to_unitigIndex[nodeIndex];



				u_int32_t contigIndex = _parent._progressiveAbundanceFilter->_unitigGraph->_nodes[unitigIndex]->_contigIndex;
				
				if(unitigIndex != prevUnitigIndex){
					prevUnitigIndex = unitigIndex;
					unitigs.push_back(unitigIndex);
				}
				
				if(contigIndex != prevContigIndex){
					prevContigIndex = contigIndex;
					contigs.push_back(contigIndex);
					contigLastPositions[contigIndex] = contigPos;
					contigPos += 1;
				}

				//if(contigIndex == 198760){
				if(unitigIndex == 136808){
					isHere = true;
				}
			}
			*/
			
			for(u_int32_t nodeIndex : readNodeIndex){
			
				if(_parent._nodeIndex_to_unitigIndex.find(nodeIndex) == _parent._nodeIndex_to_unitigIndex.end()) continue;
				
				u_int32_t unitigIndex = _parent._nodeIndex_to_unitigIndex[nodeIndex];
				
				if(unitigIndex != prevUnitigIndex){

					//if(unitigIndex == 505934 || unitigIndex == 505935){
					//	_parent._selectedReads1.insert(readIndex);
					//	cout << "allo1" << endl;
					//}
					//else if(unitigIndex == 265690 || unitigIndex == 265691){
					//	cout << "allo2" << endl;
					//	_parent._selectedReads2.insert(readIndex);
					//}

					prevUnitigIndex = unitigIndex;
					unitigs.push_back(unitigIndex);
				}

			}

			#pragma omp critical(indexReadsFunctor)
			{



				//for(u_int32_t unitigIndex : unitigs){
				//	_parent._progressiveAbundanceFilter->_unitigGraph->_nodes[unitigIndex]->_readIndexes.push_back(readIndex);
				//}
				/*
				if(_indexContigs){
						
					if(isHere){
						
						for(size_t i=0; i<contigs.size(); i++){
							
							vector<u_int32_t> contigPath;
							u_int32_t contigIndex = contigs[i];
							cout << contigIndex << " ";
						}
						cout << endl;
					}

					//unordered_set<u_int32_t> indexedContigs;
					//cout << "-----" << endl;
					//cout << unitigs.size() << endl;

					for(size_t i=0; i<contigs.size(); i++){
						
						vector<u_int32_t> contigPath;
						u_int32_t contigIndex = contigs[i];

						if(i == contigLastPositions[contigIndex]){
							
							for(size_t j=i+1; j<contigs.size(); j++){
								
								//UnitigGraph::Node* unitig2 = _parent._progressiveAbundanceFilter->_unitigGraph->_nodes[unitigs[j]];
								contigPath.push_back(contigs[j]);
							}
							

							if(contigPath.size() > 0){
								if(_parent._contigDatas[contigIndex]._isLongContig){
									_parent._contigDatas[contigIndex]._successorUnitigPaths.push_back(contigPath);
								}
							}
						}
					}


				}
				*/
				
				if(_indexUnitigs){

					/*
					bool hasUnitig = false;
					for(size_t i=0; i<unitigs.size(); i++){
						if(unitigs[i] == 3140842){
							hasUnitig = true;
							break;
						}
					}

					if(hasUnitig){
						for(size_t i=0; i<unitigs.size(); i++){
							cout << unitigs[i] << " ";
						}
						cout << endl;
					}
					*/
					
					for(size_t i=0; i<unitigs.size(); i++){

						
						vector<u_int32_t> unitigPath;
						
						//cout << unitigs[i] << endl;
						UnitigGraph::Node* unitig1 = _parent._progressiveAbundanceFilter->_unitigGraph->_nodes[unitigs[i]];
						//u_int32_t contigIndex = unitig1->_contigIndex;

						//if(unitig1->_unitigIndex == 310){
						//	cout << "lala" << " " << unitigs.size() << endl;
						//	getchar();
						//}
						
						for(size_t j=i+1; j<unitigs.size(); j++){
							
							UnitigGraph::Node* unitig2 = _parent._progressiveAbundanceFilter->_unitigGraph->_nodes[unitigs[j]];
							unitigPath.push_back(unitig2->_unitigIndex);

							//if(unitig1->_successorSupports.find(unitig2->_unitigIndex) == unitig1->_successorSupports.end()){
							//	unitig1->_successorSupports[unitig2->_unitigIndex] = 0;
							//}
							
							//unitig1->_successorSupports[unitig2->_unitigIndex] += 1;
						}
						//cout << unitig1->_successorSupports.size() << endl;


							
						if(unitigPath.size() > 0){
							unitig1->_successorUnitigPaths.push_back(unitigPath);



						}

					}
				}
				
				//getchar();
			}
		}

	};
	
	vector<bool> _isUnitigVisitedUnique;
	vector<bool> _isUnitigVisited;

	void addUnitigToPath(u_int32_t unitigIndex, vector<u_int32_t>& unitigPath){


		unitigPath.push_back(unitigIndex);
		/*
		writePoaGraphs(_parent._progressiveAbundanceFilter->_unitigGraph->_nodes[unitigIndex]);
		//for(u_int32_t unitigIndex : unitigPath){
		//}

		unitigPath.push_back(unitigIndex);
		
		int index = unitigPath.size()-1;

		if(_parent._debug_referenceUnitigPaths[_currentReferenceIndex][index] == unitigPath[index]){
			cout << "\tPath index " + to_string(index) + ": OK" << endl;
		}
		else{
			cout << "\tPath index " + to_string(index) + ": NOK" << "    utg" << unitigPath[index] << " utg" <<  _parent._debug_referenceUnitigPaths[_currentReferenceIndex][index] << endl;
		}
		*/

	}

	bool extendPath(UnitigGraph::Node* startingUnitig, vector<u_int32_t>& unitigPath){

		u_int64_t _maxBubbleLength = 75000;
		float minAbundance = startingUnitig->_abundance * 0.2;
		//_unitigPathLast.clear();
		//_unitigPathLast.insert(startingUnitig->_unitigIndex);

		_isUnitigVisitedUnique[startingUnitig->_unitigIndex] = true;
		_isUnitigVisitedUnique[_unitigGraph->unitigIndex_toReverseDirection(startingUnitig)->_unitigIndex] = true;


		UnitigGraph::Node* currentUnitig = startingUnitig;
		UnitigGraph::Node* currentUnitigRev = _unitigGraph->unitigIndex_toReverseDirection(currentUnitig);

		cout << endl << "------" << endl;

		cout << "Starting unitig: utg" <<  startingUnitig->_unitigIndex << endl;
		cout << "Abundance filter: " << startingUnitig->_abundance << " " << minAbundance << endl;

		addUnitigToPath(currentUnitig->_unitigIndex, unitigPath);
		//unitigPath.push_back(currentUnitig->_unitigIndex);

		while(true){

			cout << "Current unitig: utg" << currentUnitig->_unitigIndex << endl; 
			

			UnitigGraph::Node* selectedNextUnitig = nullptr;

			vector<UnitigGraph::Node*> validSuccessors;


			vector<u_int32_t> mostSupportedPath = getMostSupportedPath(currentUnitig);
			for(u_int32_t unitigIndex : mostSupportedPath){

				UnitigGraph::Node* nextUnitig = _unitigGraph->_nodes[unitigIndex];
				UnitigGraph::Node* nextUnitigRev = _unitigGraph->unitigIndex_toReverseDirection(nextUnitig);
	
				cout << "\t- utg" << nextUnitig->_unitigIndex << endl;
				if(isValidPath(nextUnitigRev, currentUnitigRev)){
					cout << "\t\tValid path: utg" << currentUnitig->_unitigIndex << " -> utg" << nextUnitig->_unitigIndex << endl;
					
					if(_isUnitigVisitedUnique[nextUnitig->_unitigIndex] && nextUnitig != startingUnitig){
						cout << "\t\tUnique already visited" << endl; //Visiting the starting unitig two time is allowed for circularity
					}
					else if(isPathUnique(nextUnitigRev, currentUnitigRev, _maxBubbleLength, minAbundance)){
						cout << "\t\tUnique" << endl;
						validSuccessors.push_back(nextUnitig);
					}
					else{
						cout << "\t\tAmbiguous" << endl;
					}
				}
			}

			if(validSuccessors.size() == 0){
				cout << "No valid successor" << endl;
				break;
			}

			selectedNextUnitig = validSuccessors[0];

			for(u_int32_t unitigIndex : mostSupportedPath){
				if(unitigIndex == selectedNextUnitig->_unitigIndex) break;
				addUnitigToPath(unitigIndex, unitigPath);
				//unitigPath.push_back(unitigIndex);
			}

			/*
			if(currentUnitig->_bubbleExit != nullptr && !_isUnitigVisited[currentUnitig->_bubbleExit->_unitigIndex]){

				cout << "\tFound bubble exit: utg" << currentUnitig->_bubbleExit->_unitigIndex << endl;

				vector<UnitigGraph::Node*> superbubbleNodes;
				_parent._progressiveAbundanceFilter->_unitigGraph->collectSuperbubbleNodes(currentUnitig, currentUnitig->_bubbleExit, superbubbleNodes);
			
				for(UnitigGraph::Node* superbubbleNode : superbubbleNodes){
					_isUnitigVisited[superbubbleNode->_unitigIndex] = true;
					_isUnitigVisited[_parent._progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(superbubbleNode)->_unitigIndex] = true;

					if(superbubbleNode == startingUnitig){
						return true;
					}
				}

				vector<u_int32_t> superbubblePath;
				_parent._progressiveAbundanceFilter->_unitigGraph->getSuperbubbleMostSupportedPath(currentUnitig, currentUnitig->_bubbleExit, superbubblePath);
				
				for(u_int32_t unitigIndex : superbubblePath){
					//unitigPath.push_back(unitigIndex);
					addUnitigToPath(unitigIndex, unitigPath);
					cout << "\t\tAdd bubble unitig: utg" << unitigIndex << endl;
				}

				selectedNextUnitig = currentUnitig->_bubbleExit;
				//getchar();
			}
			else{

				UnitigGraph::Node* unitigTangleExit = isTangle(currentUnitig);

				if(unitigTangleExit != nullptr && !_isUnitigVisited[unitigTangleExit->_unitigIndex]){
					
					cout << "\tFound bridged tangle exit: utg" << unitigTangleExit->_unitigIndex << endl;
					//getchar();

					for(u_int32_t unitigIndex : currentUnitig->_poaBestSupportingPath){

						if(unitigIndex == startingUnitig->_unitigIndex){
							return true;
						}
						else if(unitigIndex == unitigTangleExit->_unitigIndex){
							break;
						}
						else{
							_isUnitigVisited[unitigIndex] = true;
							_isUnitigVisited[_parent._progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(unitigIndex)] = true;
							//unitigPath.push_back(unitigIndex);
							addUnitigToPath(unitigIndex, unitigPath);
						}
					}
					
					selectedNextUnitig = unitigTangleExit;
				}
				else{

					vector<UnitigGraph::Node*> validSuccessors;

					for(u_int32_t unitigIndex : currentUnitig->_poaBestSupportingPath){

						UnitigGraph::Node* nextUnitig = _parent._progressiveAbundanceFilter->_unitigGraph->_nodes[unitigIndex];
						UnitigGraph::Node* nextUnitigRev = _parent._progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(nextUnitig);
			
						cout << "\t- utg" << nextUnitig->_unitigIndex << endl;
						if(isValidPath(nextUnitigRev, currentUnitigRev)){
							cout << "\t\tValid path: utg" << currentUnitig->_unitigIndex << " -> utg" << nextUnitig->_unitigIndex << endl;
							
							if(_isUnitigVisitedUnique[nextUnitig->_unitigIndex] && nextUnitig != startingUnitig){
								cout << "\t\tUnique already visited" << endl; //Visiting the starting unitig two time is allowed for circularity
							}
							else if(isPathUnique(nextUnitigRev, currentUnitigRev, _maxBubbleLength, minAbundance)){
								cout << "\t\tUnique" << endl;
								validSuccessors.push_back(nextUnitig);
							}
							else{
								cout << "\t\tAmbiguous" << endl;
							}
						}
					}

					if(validSuccessors.size() == 0){
						cout << "No valid successor" << endl;
						break;
					}

					selectedNextUnitig = validSuccessors[0];

					for(u_int32_t unitigIndex : currentUnitig->_poaBestSupportingPath){
						if(unitigIndex == selectedNextUnitig->_unitigIndex) break;
						addUnitigToPath(unitigIndex, unitigPath);
						//unitigPath.push_back(unitigIndex);
					}

				}


			}
			*/

			//"reprise: most supported path a ameliorer, on doit suivre le chemin le plus supporté mais aussi s'arreter si le chemin est ambigue, si le chemin est ambigue on test si c'est une bubble etc..."
			cout << "\tSelected successor: utg" << selectedNextUnitig->_unitigIndex << endl;
			currentUnitig = selectedNextUnitig;
			currentUnitigRev = _unitigGraph->unitigIndex_toReverseDirection(currentUnitig);

			if(currentUnitig == startingUnitig){
				return true;
			}
			else{
				addUnitigToPath(selectedNextUnitig->_unitigIndex, unitigPath);
				//unitigPath.push_back(selectedNextUnitig->_unitigIndex);
			}

			_isUnitigVisited[currentUnitig->_unitigIndex] = true;
			_isUnitigVisited[_unitigGraph->unitigIndex_toReverseDirection(currentUnitig)->_unitigIndex] = true;
			_isUnitigVisitedUnique[currentUnitig->_unitigIndex] = true;
			_isUnitigVisitedUnique[_unitigGraph->unitigIndex_toReverseDirection(currentUnitig)->_unitigIndex] = true;


			ofstream colorFile(_passDir + "/currentContigPath.csv");
			colorFile << "Name,Color" << endl;
			for(u_int32_t unitigIndex : unitigPath){
				colorFile << "utg" << unitigIndex << ",green" << endl;
				colorFile << "utg" << _unitigGraph->unitigIndex_toReverseDirection(unitigIndex) << ",green" << endl;
			}
			colorFile.close();

			//getchar();
		}

		return false;

	}

	bool isValidPath(UnitigGraph::Node* fromUnitig, UnitigGraph::Node* toUnitig){
						
		cout << "\tisPathValid:" << endl;
		for(u_int32_t unitigIndex : getMostSupportedPath(fromUnitig)){
			cout << "\t\tutg" << unitigIndex << endl;
			if(unitigIndex == toUnitig->_unitigIndex) return true;
		}
	
		return false;
	}

	bool isPathUnique(UnitigGraph::Node* fromUnitig, UnitigGraph::Node* toUnitig, u_int64_t maxLength, float minAbundance){

		vector<u_int32_t> mostSupportedPath = getMostSupportedPath(fromUnitig);

		for(u_int32_t unitigIndex : mostSupportedPath){
			if(unitigIndex == toUnitig->_unitigIndex) return true;
		}


		cout << "Most supported path: " << endl;
		for(u_int32_t unitigIndex : mostSupportedPath){
			cout << "\tutg" << unitigIndex << endl;
		}

		return false;
		/*
		GraphPOA::Node* node = nullptr;
		int maxWeight = 0;

		for(GraphPOA::Node* n : fromUnitig->_poaGraph->_nodes){
			if(n->_unitigIndex == fromUnitig->_unitigIndex){
				for(GraphPOA::Edge* edge : n->_successors){
					if(edge->_weight > maxWeight){
						maxWeight = edge->_weight;
						node = n;
					}
				}
			}
		}


		while(true){

			UnitigGraph::Node* unitig = _parent._progressiveAbundanceFilter->_unitigGraph->_nodes[node->_unitigIndex];
			unordered_set<u_int32_t> validUnitigSuccessors;

			for(UnitigGraph::Node* nn : unitig->_successors){
				validUnitigSuccessors.insert(nn->_unitigIndex);
			}

			maxWeight = 0;

			vector<GraphPOA::Edge*> validPoaSuccessors;
			for(GraphPOA::Edge* edge : node->_successors){
				if(validUnitigSuccessors.find(edge->_toNode->_unitigIndex) == validUnitigSuccessors.end()) continue;

				validPoaSuccessors.push_back(edge);

				if(edge->_weight > maxWeight){
					maxWeight = edge->_weight;
				}
			}

			int minSupport = ceil(maxWeight * _parent._chimericSuccessorDelta)+1;


			cout << "\t\tutg" << node->_unitigIndex << " " << validPoaSuccessors.size() << endl;
			cout << "\t\t\tMax weight: " << maxWeight << endl;
			cout << "\t\t\tChimeric threshold: " << minSupport << endl;
			for(GraphPOA::Edge* edge : validPoaSuccessors){
				cout << "\t\t\tutg" << edge->_toNode->_unitigIndex << " " << edge->_weight << endl;
			}


			vector<GraphPOA::Edge*> validPoaSuccessorsNoChimeric;
			for(GraphPOA::Edge* edge : validPoaSuccessors){
				if(edge->_weight < minSupport) continue; //Chimeric edge
				validPoaSuccessorsNoChimeric.push_back(edge);
			}

			if(validPoaSuccessorsNoChimeric.size() == 0){
				validPoaSuccessorsNoChimeric = validPoaSuccessors; //Disable chimeric filter if all sucessor filtered out (can happen if a repeat have super high value)
			}

			if(validPoaSuccessorsNoChimeric.size() == 0){ 
				cout << "\t\tNo successor" << endl;
				return false; //No successor
			}
			else{
				
				

				if(validPoaSuccessorsNoChimeric.size() == 1){
				}
				else if(validPoaSuccessorsNoChimeric.size() >= 2){
					cout << "\t\tMultiple successors" << endl;
					
					GraphPOA::Edge* uniqueBubbleSuccessor = isUniqueBubbleSuccessor(unitig, validPoaSuccessorsNoChimeric);
					
					if(uniqueBubbleSuccessor == nullptr){
						return false; //Current unitig successors is ambiguous (repeat)
					}

					cout << "\t\tMultiple successors are bubble" << endl;
					//return false; //Current unitig successors is ambiguous (repeat)
					validPoaSuccessorsNoChimeric.clear();
					validPoaSuccessorsNoChimeric.push_back(uniqueBubbleSuccessor);

					//getchar();
				}


				GraphPOA::Edge* mostSupportedSuccessor =  validPoaSuccessorsNoChimeric[0];
				if(mostSupportedSuccessor->_weight <= 0){
					cout << "\t\tNot enough support" << endl;
					return false; //Not enough support
				}

				node = mostSupportedSuccessor->_toNode;

				if(node->_unitigIndex == toUnitig->_unitigIndex) return true;
			}
		}
		
		return false;
		*/
	}




	vector<u_int32_t> getMostSupportedPath(UnitigGraph::Node* unitig){

		vector<u_int32_t> path;

		if(unitig->_successorUnitigPaths.size() == 0) return path;

		vector<vector<u_int32_t>> unitigPaths = unitig->_successorUnitigPaths;
		std::sort(unitigPaths.begin(), unitigPaths.end(), [](const vector<u_int32_t> & a, const vector<u_int32_t> & b){ return a.size() > b.size(); });

		spoa64::Graph graph{};
		
		int i = 0;
		for(const vector<u_int32_t>& unitigPath: unitigPaths){

			vector<MinimizerType> unitigPath64(unitigPath.size(), -1);
			for(size_t i=0; i<unitigPath.size(); i++){
				unitigPath64[i] = unitigPath[i];
			}
			unitigPath64.insert(unitigPath64.begin(), unitig->_unitigIndex);

			vector<MinimizerType> weights(unitigPath64.size(), 1);

			if(i == 0){
				graph.AddAlignment(spoa64::Alignment(), unitigPath64, unitigPath64.size(), weights);

			}
			else{


				spoa64::Alignment alignment = _alignmentEngine->Align(unitigPath64, unitigPath64.size(), graph);

				graph.AddAlignment(alignment, unitigPath64, unitigPath64.size(), weights);

			}

			i+= 1;
		}

		savePoaGraph(graph, _passDir + "/poaGraph.gfa");

		path = computePath(graph, unitig->_unitigIndex);

		return path;
	}
	

	void savePoaGraph(const spoa64::Graph& graph, const string& outputFilename){

		ofstream outputFile(outputFilename);

		ofstream colorFile(outputFilename + ".color.csv");
		colorFile << "Name,Color" << endl;
		
		ofstream edgeFile(outputFilename + ".edge.csv");
		edgeFile << "Name,Edge" << endl;

		//std::cout << "H\tVN:Z:1.0" << std::endl;
		for (const auto& it : graph.nodes()) {

			MinimizerType minimizer = graph.decoder(it->code);
			string id = to_string(it->id) + "-" + to_string(graph.decoder(it->code));

			//graph.decoder(it->code)
			outputFile << "S\t" << id << "\t" << "*" << "\t" << "LN:i:500" << "\t" << "dp:i:" << it->Coverage() << endl;
			//if (is_consensus_node[it->id]) {
			//std::cout << "\tic:Z:true";
			//}
			//std::cout << std::endl;

			string edgeStr = "";

			for (const auto& jt : it->outedges) {
				string edgeId = to_string(jt->head->id) + "-" + to_string(graph.decoder(jt->head->code));
				outputFile << "L\t" << id << "\t" << "+\t" << edgeId << "\t" << "+\t" << "1M" << endl; //\t" << "ew:f:" << jt->weight << endl;
				//if (is_consensus_node[it->id] &&
				//	is_consensus_node[jt->head->id]) {
				//	std::cout << "\tic:Z:true";
				//}
				//std::cout << std::endl;
				if(jt->weight > 2){
					edgeStr += "[" + to_string(jt->head->id) + " - " + to_string(jt->weight) + "] ";
				}
			}

			edgeFile << (id) << "," << edgeStr << endl; 
		}

		outputFile.close();
		colorFile.close();
		edgeFile.close();
	}

	bool _print_debug;
	u_int64_t _minPoaNodeCoverage;

	struct DereplicatedEdge{
		spoa64::Graph::Edge* _edge;
		u_int64_t _weight;
	};
	
	vector<u_int32_t> computePath(const spoa64::Graph& graph, u_int32_t startUnitigIndex){


		PathWeight bestPath = {{}, 0};

		for (const auto& it : graph.nodes()) {

			u_int32_t unitigIndex = graph.decoder(it->code);

			if(unitigIndex != startUnitigIndex) continue;

			if(it->Coverage() < _minPoaNodeCoverage) continue;
			//string id = to_string(it->id) + "-" + to_string(graph.decoder(it->code));

			PathWeight possiblePath = computePathSub(graph, startUnitigIndex, graph.nodes_[it->id].get());

			if(possiblePath._weight > bestPath._weight){
				bestPath = possiblePath;
			}
		}

		return bestPath._path;
	}

	struct PathWeight{
		vector<u_int32_t> _path;
		float _weight;
	};

	PathWeight computePathSub(const spoa64::Graph& graph, u_int32_t startUnitigIndex, spoa64::Graph::Node* startNode){

		if(_print_debug) cout << "Compute path: " << startNode->id << "-" << graph.decoder(startNode->code) << endl;
		
		PathWeight path = {{}, 0};
		//unordered_set<UnitigGraph::Node*> isVisited;
		//list<PathSuccessor> queue;
		//unordered_map<u_int64_t, u_int64_t> prev;
	
		//isVisited.insert(unitigSource);
		//isVisited.insert(_progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(unitigSource));

		//queue.push_back({fromNode, 0});
		//prev[fromNode] = -1;
		//bool found = false;
		spoa64::Graph::Node* currentNode = startNode;
		//path._path.push_back(graph.decoder(currentNode->code));
		
		while (true) {
	

			if(_print_debug) cout << "\tVisit: " << currentNode->id << "-" << graph.decoder(currentNode->code) << endl;

			vector<DereplicatedEdge> successors = getSolidSuccessors(graph, currentNode);

	


			for(const DereplicatedEdge& succ : successors){
				if(_print_debug) cout << "\t\t" << succ._edge->head->id << "-" << graph.decoder(succ._edge->head->code) << " " << succ._edge->weight << endl;

			}

			if(successors.size() == 0) break;
			//if(successors.size() == 0) return nullPath;
			if(successors.size() > 2 && successors[0]._weight <= 2) return path;

			currentNode = successors[0]._edge->head;

			if(graph.decoder(currentNode->code) == startUnitigIndex) break; //join a starting node

			path._weight += successors[0]._edge->weight;
			path._path.push_back(graph.decoder(currentNode->code));


			//if(path.size() > maxSize) return nullPath;
			
		}


		
		return path;
	}

	vector<DereplicatedEdge> getSolidSuccessors(const spoa64::Graph& graph, spoa64::Graph::Node* node){

		/*
		u_int32_t unitigIndex = graph.decoder(node->code);

		unordered_set<u_int32_t> validUnitigSuccessors;
		for(UnitigGraph::Node* nn : _unitigGraph->_nodes[unitigIndex]->_successors){
			validUnitigSuccessors.insert(nn->_unitigIndex);

			if(_print_debug) cout << "\tValid assembly successors: utg" << nn->_unitigIndex << endl;
		}
		*/

		vector<DereplicatedEdge> successors = getDereplicatedSuccessors(graph, node);

		std::sort(successors.begin(), successors.end(), [](const DereplicatedEdge a, const DereplicatedEdge b){
			return a._weight > b._weight;
		});

		float maxWeight = 0;

		for(const DereplicatedEdge& edge : successors){
			//if(edge._weight == 1 && edge._edge->head->Coverage() <= 1) continue;

			if(edge._weight > maxWeight){
				maxWeight = edge._weight;
			}
		}

		float minWeight = maxWeight * 0.5;
		if(_print_debug) cout << "\t\tMin weight: " << minWeight << endl;

		vector<DereplicatedEdge> solidEdges;

		for(const DereplicatedEdge& edge : successors){
			//if(edge._weight == 1 && edge._edge->head->Coverage() <= 1) continue;

			//if(edge._edge->head->Coverage() < _minPoaNodeCoverage) continue;
			if(edge._weight < _minPoaNodeCoverage) continue;
			if(edge._weight < minWeight) continue;
				
			u_int32_t unitigIndexSuccessor = graph.decoder(edge._edge->head->code);
			
			//cout << "utg" << unitigIndexSuccessor << endl;

			//if(validUnitigSuccessors.find(unitigIndexSuccessor) == validUnitigSuccessors.end()) continue;

			solidEdges.push_back(edge);
		}

		return solidEdges;
	}


	struct SuccessorCompletion{
		spoa64::Graph::Edge* _edge;
		u_int64_t _completion;
	};

	vector<DereplicatedEdge> getDereplicatedSuccessors(const spoa64::Graph& graph, spoa64::Graph::Node* node){


		for(spoa64::Graph::Edge* edge : node->outedges){
			//if(_print_debug) cout << "\t\tPossible successor: " << edge->head->id << " " << edge->weight << " " << computeEdgeCompletion(graph, edge, readMinimizers) << endl;
			//for(spoa64::Graph::Edge* edge2 : edge->head->outedges){
			//	cout << "\t\t\tLala: " << edge2->head->id << " " << edge2->weight << " " << computeEdgeCompletion(graph, edge2, readMinimizers) << endl;
			//}
		}

		unordered_map<MinimizerType, SuccessorCompletion> minimizer_to_edge;
		unordered_map<MinimizerType, u_int64_t> minimizer_to_weight;

		for(spoa64::Graph::Edge* edge : getPoaAssemblySuccessor(graph, node)){//->outedges){

			MinimizerType minimizer = edge->head->code;

			minimizer_to_weight[minimizer] += edge->weight;
			u_int64_t completion = computeEdgeCompletion(graph, edge);

			if(_print_debug) cout << "\t\tPossible successor: " << edge->head->id << " utg" << graph.decoder(edge->head->code) << " Completion: " << completion << endl;

			if(minimizer_to_edge.find(minimizer) != minimizer_to_edge.end()){
				if(completion > minimizer_to_edge[minimizer]._completion){
					minimizer_to_edge[minimizer] = {edge, completion};
				}
			}
			else{
				minimizer_to_edge[minimizer] = {edge, completion};
			}

		}


		vector<DereplicatedEdge> derepEdges;

		for(const auto& it : minimizer_to_edge){

			MinimizerType minimizer = it.first;
			spoa64::Graph::Edge* edge = it.second._edge;

			u_int64_t derepWeight = minimizer_to_weight[minimizer];

			//edge->weight = minimizer_to_weight[it.first];
			derepEdges.push_back({edge, derepWeight});
			

			if(_print_debug) cout << "\t\tPossible successor derep: " << edge->head->id << " utg" << graph.decoder(edge->head->code) << " " << derepWeight << endl;
		}

		return derepEdges;

	}


	u_int64_t computeEdgeCompletion(const spoa64::Graph& graph, spoa64::Graph::Edge* edge){


		list<spoa64::Graph::Node*> queue;

		queue.push_back(edge->head);
		unordered_set<MinimizerType> isVisited;
	
		u_int64_t nbVisitedReadNodes = edge->weight;

		while (!queue.empty()) {
	
			spoa64::Graph::Node* currentNode = queue.front();
			queue.pop_front();

			if(isVisited.find(currentNode->id) != isVisited.end()) continue;
			isVisited.insert(currentNode->id);


			for(spoa64::Graph::Edge* successor : getPoaAssemblySuccessor(graph, currentNode)){

				if(successor->head->Coverage() < _minPoaNodeCoverage) continue;
				//if(readMinimizers.find(graph.decoder(successor->head->code)) != readMinimizers.end()){
				//	nbVisitedReadNodes += successor->weight;
				//}
				nbVisitedReadNodes += successor->weight;

				queue.push_back(successor->head);
			}


		}

		return nbVisitedReadNodes;
	}

	vector<spoa64::Graph::Edge*> getPoaAssemblySuccessor(const spoa64::Graph& graph, spoa64::Graph::Node* node){

		vector<spoa64::Graph::Edge*> successors;

		u_int32_t unitigIndex = graph.decoder(node->code);

		unordered_set<u_int32_t> validUnitigSuccessors;
		for(UnitigGraph::Node* nn : _unitigGraph->_nodes[unitigIndex]->_successors){
			validUnitigSuccessors.insert(nn->_unitigIndex);
		}

		for(spoa64::Graph::Edge* successor : node->outedges){
			successors.push_back(successor);
		}

		return successors;
	}


};





#endif