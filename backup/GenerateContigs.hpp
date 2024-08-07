
#ifndef MDBG_METAG_GENERATECONTIGS
#define MDBG_METAG_GENERATECONTIGS

//#define PRINT_DEBUG_PREVRANK

#include "Commons.hpp"

#include <string>
#include <sstream>
#include <unordered_set>
#include <unordered_map>
#include <regex>
#include <algorithm>
#include <libgen.h>
#include <set>
//#include "graph/Graph.hpp"
#include "graph/GfaParser.hpp"
#include "graph/GraphSimplify.hpp"
#include "eval/ContigStatistics.hpp"
//#include "toBasespace/ToBasespaceOnTheFly.hpp"
//#include "contigFeatures/ContigFeature.hpp"










class GenerateContigs : public Tool{

public:

	string _inputDir;
	string _truthInputFilename;
	string _outputFilename;
	string _outputFilename_complete;
	bool _debug;
	string _inputFilename_unitigNt;
	string _inputFilename_unitigCluster;
	string _filename_abundance;

	float _minimizerDensity;
    size_t _minimizerSize;
    size_t _kminmerSize;
	vector<UnitigData> _unitigDatas;

	string _partitionDir;
	string _filename_solidKminmers;
	string _filename_inputContigs;
	string _filename_outputContigs;
	string _filename_readMinimizers;
	string _filename_hifiasmGroundtruth;
	unordered_map<KmerVec, u_int16_t> _evaluation_hifiasmGroundTruth;
	unordered_map<KmerVec, u_int32_t> _evaluation_hifiasmGroundTruth_position;
	unordered_map<u_int32_t, u_int32_t> _evaluation_hifiasmGroundTruth_nodeNamePosition;
	gzFile _outputContigFile;
	gzFile _outputContigFile_complete;
	MDBG* _mdbg;
	GraphSimplify* _graph;
	//ToBasespaceOnTheFly _toBasespace;
	//ContigFeature _contigFeature;
	string _gfaFilename;

	GenerateContigs(): Tool (){

	}

	~GenerateContigs(){

	}

	void execute (){

		loadGraph();
		generateContigs2(_inputDir + "/contigs.nodepath.gz", _inputDir + "/contigs.fasta.gz");

	}

	void parseArgs(int argc, char* argv[]){



		cxxopts::Options options("Assembly", "");
		options.add_options()
		(ARG_OUTPUT_DIR, "", cxxopts::value<string>())
		(ARG_INPUT_FILENAME_CONTIG, "", cxxopts::value<string>()->default_value(""))
		(ARG_INPUT_FILENAME_TRUTH, "", cxxopts::value<string>()->default_value(""))
		(ARG_DEBUG, "", cxxopts::value<bool>()->default_value("false"))
		(ARG_INPUT_FILENAME_UNITIG_NT, "", cxxopts::value<string>()->default_value(""))
		(ARG_INPUT_FILENAME_UNITIG_CLUSTER, "", cxxopts::value<string>()->default_value(""))
		(ARG_INPUT_FILENAME_ABUNDANCE, "", cxxopts::value<string>()->default_value(""));



		if(argc <= 1){
			cout << options.help() << endl;
			exit(0);
		}

		cxxopts::ParseResult result;

		try{
			result = options.parse(argc, argv);

			_inputDir = result[ARG_OUTPUT_DIR].as<string>();
			_filename_inputContigs = result[ARG_INPUT_FILENAME_CONTIG].as<string>();
			_truthInputFilename = result[ARG_INPUT_FILENAME_TRUTH].as<string>();
			_inputFilename_unitigNt = result[ARG_INPUT_FILENAME_UNITIG_NT].as<string>();
			_inputFilename_unitigCluster = result[ARG_INPUT_FILENAME_UNITIG_CLUSTER].as<string>();
			_debug = result[ARG_DEBUG].as<bool>();
			_filename_abundance = result[ARG_INPUT_FILENAME_ABUNDANCE].as<string>();
		}
		catch (const std::exception& e){
			std::cout << options.help() << std::endl;
			std::cerr << e.what() << std::endl;
			std::exit(EXIT_FAILURE);
		}


		string filename_parameters = _inputDir + "/parameters.gz";
		gzFile file_parameters = gzopen(filename_parameters.c_str(),"rb");
		gzread(file_parameters, (char*)&_minimizerSize, sizeof(_minimizerSize));
		gzread(file_parameters, (char*)&_kminmerSize, sizeof(_kminmerSize));
		gzread(file_parameters, (char*)&_minimizerDensity, sizeof(_minimizerDensity));
		gzclose(file_parameters);

		cout << endl;
		cout << "Input dir: " << _inputDir << endl;
		//cout << "Output filename: " << _outputFilename << endl;
		cout << "Minimizer length: " << _minimizerSize << endl;
		cout << "Kminmer length: " << _kminmerSize << endl;
		cout << "Density: " << _minimizerDensity << endl;
		cout << endl;

		_outputFilename = _inputDir + "/minimizer_contigs.gz";
		_outputFilename_complete = _inputDir + "/minimizer_contigs_complete.gz";
		_filename_readMinimizers = _inputDir + "/read_data.txt";
		_filename_hifiasmGroundtruth = _inputDir + "/hifiasmGroundtruth.gz";
		_filename_outputContigs = _inputDir + "/contigs.min.gz";
		_filename_solidKminmers = _inputDir + "/solid.min.gz";

	}

	unordered_map<u_int64_t, vector<u_int32_t>> _debug_readPath;
    vector<bool> _isBubble;

	void loadGraph(){

		_gfaFilename = _inputDir + "/minimizer_graph.gfa";
		string gfa_filename_noUnsupportedEdges = _inputDir + "/minimizer_graph_noUnsupportedEdges.gfa";
		//string gfa_filename_unitigs = _inputDir + "/minimizer_graph_unitigs.gfa";
		string mdbg_filename = _inputDir + "/mdbg_nodes.gz";




		
		cout << _gfaFilename << endl;
		_mdbg = new MDBG(_kminmerSize);
		_mdbg->load(mdbg_filename);
		cout << "Nb nodes: " <<  _mdbg->_dbg_nodes.size() << endl;

		//extractContigKminmers2(_inputDir + "/contigs.fasta.gz");

		if(_truthInputFilename != ""){
			extract_truth_kminmers();
		}


		//file_groundTruth = ofstream(_inputDir + "/binning_results.csv");
		//file_groundTruth << "Name,Colour" << endl;

		//file_groundTruth_hifiasmContigs = ofstream(_inputDir + "/binning_results_hifiasm_" + to_string(_kminmerSize) + ".csv");
		//file_groundTruth_hifiasmContigs << "Name,Colour" << endl;

		//if(_debug){
            //gfa_filename = _inputDir + "/minimizer_graph_debug.gfa";
		//}
		
		GraphSimplify* graphSimplify = new GraphSimplify(_gfaFilename, _inputDir, 0, _kminmerSize);
		_graph = graphSimplify;
		

		//Generate unitigs
		cout << "Indexing reads" << endl;
		_unitigDatas.resize(_mdbg->_dbg_nodes.size());
		_graph->clear(0);
		_graph->compact(false, _unitigDatas);
		removeUnsupportedEdges(_gfaFilename, gfa_filename_noUnsupportedEdges, _graph);
		
		delete _mdbg;

		cout << "done" << endl;
	

		_graph->debug_writeGfaErrorfree(500, 500, -1, _kminmerSize, false, true, false, _unitigDatas, true, false, false, false);
			
	}

	void indexReads_read(const vector<u_int64_t>& minimizers, const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex){//}, const vector<KmerVec>& kminmers_k3, const vector<ReadKminmer>& kminmersInfos_k3){

		//vector<ReadIndexType> unitigIndexex;

		for(const KmerVec& vec : kminmers){
			//if(mdbg->_dbg_nodes[vec]._index == 55479) cout << "AAAAA" << endl;

			//cout << mdbg->_dbg_nodes[vec]._index << endl;
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()) continue;

			//if(vec.isPalindrome()) cout << _mdbg->_dbg_nodes[vec]._index << endl;
			//getchar();


			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
			u_int32_t nodeIndex = BiGraph::nodeName_to_nodeIndex(nodeName, true);
			//if(_graph->_isNodeValid2.find(nodeIndex) == _graph->_isNodeValid2.end()) continue;
			//if(_nodeData.find(nodeName) == _nodeData.end()) continue;

			//UnitigData& unitigData = _nodeData[nodeName];
			UnitigData& unitigData = _unitigDatas[nodeName];
			//if(std::find(unitigData._readIndexes.begin(), unitigData._readIndexes.end(), readIndex) != unitigData._readIndexes.end()) continue;
			unitigData._readIndexes.push_back(readIndex);
			//cout << "indexing : " << readIndex << endl;
		}



	}

	void indexReads_contig(const vector<u_int64_t>& minimizers, const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex){//}, const vector<KmerVec>& kminmers_k3, const vector<ReadKminmer>& kminmersInfos_k3){

		readIndex += 2000000000;
		//vector<ReadIndexType> unitigIndexex;

		for(const KmerVec& vec : kminmers){
			//if(mdbg->_dbg_nodes[vec]._index == 55479) cout << "AAAAA" << endl;

			//cout << mdbg->_dbg_nodes[vec]._index << endl;
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()) continue;

			//if(vec.isPalindrome()) cout << _mdbg->_dbg_nodes[vec]._index << endl;
			//getchar();


			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
			u_int32_t nodeIndex = BiGraph::nodeName_to_nodeIndex(nodeName, true);
			//if(_graph->_isNodeValid2.find(nodeIndex) == _graph->_isNodeValid2.end()) continue;
			//if(_nodeData.find(nodeName) == _nodeData.end()) continue;

			//UnitigData& unitigData = _nodeData[nodeName];
			UnitigData& unitigData = _unitigDatas[nodeName];
			//if(std::find(unitigData._readIndexes.begin(), unitigData._readIndexes.end(), readIndex) != unitigData._readIndexes.end()) continue;
			unitigData._readIndexes.push_back(readIndex);
			//cout << "indexing : " << readIndex << endl;
		}



	}
	
	void removeUnsupportedEdges(const string& gfaFilename, const string& gfa_filename_noUnsupportedEdges, GraphSimplify* graph){

		KminmerParser parser(_filename_readMinimizers, _minimizerSize, _kminmerSize, true);
		//auto fp = std::bind(&Assembly::indexReads_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
		auto fp = std::bind(&GenerateContigs::indexReads_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
		parser.parse(fp);
		
		//if(_filename_inputContigs != ""){
		//	KminmerParser parserContig(_filename_inputContigs, _minimizerSize, _kminmerSize, true);
		//	auto fpContig = std::bind(&Assembly3::indexReads_contig, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
		//	parserContig.parse(fpContig);
		//}

		if(_kminmerSize > 3) return;

		vector<DbgEdge> removedEdges;

		for(Unitig& unitig : graph->_unitigs){

			//bool nodeName_ori;
			u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(unitig._endNode);

			vector<u_int32_t> successors;
			graph->getSuccessors(unitig._endNode, 0, successors);

			for(u_int32_t successor : successors){
				
				//bool nodeName_succ_ori;
				u_int32_t nodeName_succ = BiGraph::nodeIndex_to_nodeName(successor);

				if(Utils::shareAnyRead(_unitigDatas[nodeName], _unitigDatas[nodeName_succ])){
					continue;
				}

				Unitig& unitig_succ = graph->_unitigs[graph->nodeIndex_to_unitigIndex(successor)];

				//graph->_graphSuccessors->removeEdge(nodeName, nodeName_ori, nodeName_succ, nodeName_succ_ori);
				removedEdges.push_back({unitig._endNode, successor});
				
				
			}



		}


		cout << "nb edges: " << graph->_graphSuccessors->_nbEdges << endl;

		for(const DbgEdge& edge: removedEdges){

			bool edgeExist = graph->_graphSuccessors->removeEdge(edge._from, edge._to);
			//if(!edgeExist) cout << "nop 1" << endl;
			edgeExist = graph->_graphSuccessors->removeEdge(GraphSimplify::nodeIndex_toReverseDirection(edge._to), GraphSimplify::nodeIndex_toReverseDirection(edge._from));
			//if(!edgeExist) cout << "nop 2" << endl;
			
			//graph->_isEdgeUnsupported.insert(edge);
		}
		
		cout << "nb edges: " << graph->_graphSuccessors->_nbEdges << endl;
		
		cout << "Nb unsupported edges: " << graph->_isEdgeUnsupported.size() << endl;


	}



	size_t _iter;
	//ofstream file_groundTruth;
	//ofstream file_groundTruth_hifiasmContigs;

	ofstream file_kminmersContigs;
	


	MinimizerParser* _minimizerParser;
	u_int32_t _extract_truth_kminmers_read_position;
	unordered_map<u_int32_t, vector<string>> _evaluation_hifiasmGroundTruth_nodeName_to_unitigName;
	vector<u_int32_t> _evaluation_hifiasmGroundTruth_path;
	unordered_set<u_int32_t> _hifiasm_startingNodenames;
	ofstream _file_groundTruth_hifiasm_position;

	void extract_truth_kminmers_read(kseq_t* read, u_int64_t readIndex){
		//ottalSize += strlen(read->seq.s);


		string rleSequence;
		vector<u_int64_t> rlePositions;
		Encoder::encode_rle(read->seq.s, strlen(read->seq.s), rleSequence, rlePositions);

		vector<u_int64_t> minimizers;
		vector<u_int64_t> minimizers_pos;
		_minimizerParser->parse(rleSequence, minimizers, minimizers_pos);

		vector<KmerVec> kminmers; 
		vector<ReadKminmer> kminmersInfo;
		MDBG::getKminmers(_minimizerSize, _kminmerSize, minimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, readIndex, false);

		for(size_t i=0; i<kminmers.size(); i++){

			KmerVec& vec = kminmers[i];
			//if(_mdbg->_dbg_nodes.find(kminmers[i]) == _mdbg->_dbg_nodes.end()) continue;

			//u_int32_t nodeName = _mdbg->_dbg_nodes[kminmers[i]]._index;
			if(_mdbg->_dbg_nodes.find(vec) != _mdbg->_dbg_nodes.end()){

				u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
				//2262408
				if("utg009434l" == string(read->name.s)){
					cout << nodeName << endl;
					_hifiasm_startingNodenames.insert(nodeName);
				}

				_evaluation_hifiasmGroundTruth_nodeName_to_unitigName[_mdbg->_dbg_nodes[vec]._index].push_back(string(read->name.s));
				_evaluation_hifiasmGroundTruth_path.push_back(_mdbg->_dbg_nodes[vec]._index);

				if(_evaluation_hifiasmGroundTruth_nodeNamePosition.find(_mdbg->_dbg_nodes[vec]._index) == _evaluation_hifiasmGroundTruth_nodeNamePosition.end()){
					_evaluation_hifiasmGroundTruth_nodeNamePosition[_mdbg->_dbg_nodes[vec]._index] = _extract_truth_kminmers_read_position;
				}
				//cout << mdbg->_dbg_nodes[vec]._index << " " << sequence.getComment() << endl;
				
				_file_groundTruth_hifiasm_position << nodeName << "," << _extract_truth_kminmers_read_position << endl;
			}
			else{
				//cout << "Not found position: " << position << endl;
				_evaluation_hifiasmGroundTruth_path.push_back(-1);
			}


			if(_evaluation_hifiasmGroundTruth.find(vec) != _evaluation_hifiasmGroundTruth.end()){
				//cout << position << endl;
				continue;
			}


			//_evaluation_hifiasmGroundTruth[vec] = datasetID;
			_evaluation_hifiasmGroundTruth_position[vec] = _extract_truth_kminmers_read_position;
			_extract_truth_kminmers_read_position += 1;

			
			//cout << _extract_truth_kminmers_read_position << endl;
		}


	}


	void extract_truth_kminmers(){

		_file_groundTruth_hifiasm_position.open(_inputDir + "/groundtruth_hifiasm_position.csv");
		_file_groundTruth_hifiasm_position << "Name,Position" << endl;

		_extract_truth_kminmers_read_position = 0;
		_minimizerParser = new MinimizerParser(_minimizerSize, _minimizerDensity);
		
		auto fp = std::bind(&GenerateContigs::extract_truth_kminmers_read, this, std::placeholders::_1, std::placeholders::_2);
		ReadParser readParser(_truthInputFilename, true, false);
		readParser.parse(fp);

		_file_groundTruth_hifiasm_position.close();

		delete _minimizerParser;
	}


	struct Contig{
		vector<u_int32_t> _nodePath;
		vector<u_int32_t> _nodePath_sorted;
	};


	bool isContigAssembled(const vector<u_int32_t>& nodePath, unordered_map<u_int32_t, u_int32_t>& processedNodeNames){

		unordered_set<u_int32_t> distinctContigIndex;
		u_int64_t nbBinnedNodes = 0;

		for(u_int32_t nodeIndex : nodePath){
			u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
			if(processedNodeNames.find(nodeName) != processedNodeNames.end()){
				distinctContigIndex.insert(processedNodeNames[nodeName]);
			}
			else{
				return false;
			}
		}

		//float sharedRate = nbBinnedNodes / ((float) nodePath.size());

		//return sharedRate > 0.98;
		return distinctContigIndex.size() == 1;
	}

	/*
	ofstream _fileTestLala;

	void extractContigKminmers (const string& outputFilename_fasta){
		
		_fileTestLala = ofstream(outputFilename_fasta + ".nodes.csv");
		_fileTestLala << "Name,Color" << endl;

		_mdbg = new MDBG(_kminmerSize);
		_mdbg->load(_inputDir + "/mdbg_nodes.gz");

		ReadParser parser(outputFilename_fasta, true, _minimizerSize, _kminmerSize, _minimizerDensity);
		auto fp = std::bind(&Assembly3::extractContigKminmers_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6);
		parser.parseKminmers(fp);

		for(auto& it: _contigFeature._nodeNameDuplicate){
			if(it.second == 1){
				_fileTestLala << it.first << "," << "blue" << endl;
			}
			else{
				_fileTestLala << it.first << "," << "red" << endl;
			}
		}

		delete _mdbg;
		_fileTestLala.close();
	}


	void extractContigKminmers_read(const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex, u_int64_t datasetIndex, const string& header, const string& seq){

		for(size_t i=0; i<kminmersInfos.size(); i++){
			

			KmerVec vec = kminmers[i];
			//for(KmerVec& vec : kminmers){
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()) continue;

			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
			_contigFeature._nodeNameDuplicate[nodeName] += 1;
			//_fileTestLala << nodeName << "," << "0" << endl;
			//_kminmerCounts[nodeName] = _countsInit;
			
		}

	}

	ofstream _file_contigToNode;

	void extractContigKminmers2 (const string& outputFilename_fasta){

		_file_contigToNode = ofstream(_inputDir + "/nodeToContig.csv");
		_file_contigToNode << "Name,Color" << endl;
		_contigFeature.loadAbundanceFile(_filename_abundance);

		ReadParser parser(outputFilename_fasta, true, _minimizerSize, _kminmerSize, _minimizerDensity);
		auto fp = std::bind(&Assembly3::extractContigKminmers_read2, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6);
		parser.parseKminmers(fp);

		_file_contigToNode.close();
		//vector<u_int32_t> nodePath = {933376, 1651014, 1762772, 732878};
		//vector<float> lala1;
		//vector<float> lala2;
		//_contigFeature.sequenceToAbundance(nodePath, lala1, lala2);

		//exit(1);
	}


	void extractContigKminmers_read2(const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex, u_int64_t datasetIndex, const string& header, const string& seq){

		if(seq.size() < 15000) return;
		//cout << seq.size() << endl;
		//if(readIndex == 20010) getchar(); 

		vector<float> composition;
		_contigFeature.sequenceToComposition(seq, composition);
		_contigFeature._contigCompositions[readIndex] = composition;

		unordered_set<u_int32_t> nodeNames;

		for(size_t i=0; i<kminmersInfos.size(); i++){
			

			KmerVec vec = kminmers[i];
			//for(KmerVec& vec : kminmers){
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()) continue;
			
			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
			//_file_contigToNode << nodeName << "," << readIndex << endl;
			//if(nodeName == 933376) cout << "HAHAHAHA" << endl;

			//if(_contigFeature._nodeName_to_contigIndex.find(nodeName) == _contigFeature._nodeName_to_contigIndex.end()) continue;

			_contigFeature._nodeNameDuplicate[nodeName] += 1;
			_contigFeature._nodeName_to_contigIndex[nodeName] = readIndex;
			//_fileTestLala << nodeName << "," << "0" << endl;
			//_kminmerCounts[nodeName] = _countsInit;
			
			nodeNames.insert(nodeName);
		}
		
		//std::sort(nodeNames.begin(), nodeNames.end());
		_contigFeature._contigNodes[readIndex] = nodeNames;

	}
	*/

	static bool ContigComparator_ByLength(const Contig &a, const Contig &b){
		return a._nodePath.size() > b._nodePath.size();
	}

	static bool ContigComparator_ByLength_Reverse(const Contig &a, const Contig &b){
		return a._nodePath.size() > b._nodePath.size();
	}

	void dereplicateContig2(vector<Contig>& contigs, vector<u_int32_t> nodePath, unordered_set<u_int32_t>& processedNodeNames, unordered_map<u_int32_t, u_int32_t>& nodeName_to_contigIndex, u_int64_t& contigIndex){


		if(nodePath.size() < _kminmerSize*2) return;

		vector<u_int32_t> component;
		for(u_int32_t nodeIndex : nodePath){
			component.push_back(BiGraph::nodeIndex_to_nodeName(nodeIndex));
		}
		//for(size_t i=_kminmerSize-1; i<nodePath.size()-(_kminmerSize-1); i++){
		//	component.push_back(BiGraph::nodeIndex_to_nodeName(nodePath[i]));
		//}
		std::sort(component.begin(), component.end());

		bool isValid = true;

		std::sort(contigs.begin(), contigs.end(), ContigComparator_ByLength);

		for(Contig& contig : contigs){
			if(contig._nodePath.size() == 0) continue;

			u_int64_t nbShared = Utils::computeSharedElements(component, contig._nodePath_sorted);
			if(nbShared == 0) continue;

			if(contig._nodePath.size() > component.size()){
				return;
			}
			else{
				contig._nodePath.clear();
				contig._nodePath_sorted.clear();
			}
		}
		
		addContig(nodePath, contigs, processedNodeNames, nodeName_to_contigIndex, contigIndex);
	}

	void addContig(vector<u_int32_t>& nodePath, vector<Contig>& contigs, unordered_set<u_int32_t>& processedNodeNames, unordered_map<u_int32_t, u_int32_t>& nodeName_to_contigIndex, u_int64_t& contigIndex){

		if(nodePath.size() <= _kminmerSize*2) return;

		vector<u_int32_t> nodePath_sorted;
		for(u_int32_t nodeIndex : nodePath){
			nodePath_sorted.push_back(BiGraph::nodeIndex_to_nodeName(nodeIndex));
		}

		//for(size_t i=_kminmerSize-1; i<nodePath.size()-(_kminmerSize-1); i++){
		//	nodePath_sorted.push_back(BiGraph::nodeIndex_to_nodeName(nodePath[i]));
		//}

		std::sort(nodePath_sorted.begin(), nodePath_sorted.end());
		
		contigs.push_back({nodePath, nodePath_sorted});



		//bool isLala = false;
		for(u_int32_t nodeIndex : nodePath){
			u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
			nodeName_to_contigIndex[nodeName] = contigIndex;
			//if(nodeName == 289994){
			//	isLala = true;
			//}
			processedNodeNames.insert(nodeName);
		}

		contigIndex += 1;

		/*
		if(isLala){
			cout << "Size: " << nodePath.size() << endl;
			for(u_int32_t nodeIndex : nodePath){
				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
				cout << nodeName << endl;
			}


			unordered_set<u_int32_t> component;
			_graph->getConnectedComponent_unitig(_graph->nodeIndex_to_unitigIndex(BiGraph::nodeName_to_nodeIndex(289994, false)), component);

			unordered_set<u_int32_t> validNodes;
			for (u_int32_t unitigIndex : component){
				for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
					u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
					validNodes.insert(nodeName);
				}
			}
			string outputFilename = _inputDir + "/minimizer_graph_sub_2.gfa";
			GfaParser::rewriteGfa_withoutNodes(_gfaFilename, outputFilename, validNodes, _graph->_isEdgeRemoved, _graph->_graphSuccessors);
			

			getchar();
		}
		*/


	}


	static bool UnitigComparator_ByLength2(const UnitigLength &a, const UnitigLength &b){
		return a._length > b._length;
	}

	float _minUnitigAbundance;

	void generateContigs2(const string& outputFilename, const string& outputFilename_fasta){

		string clusterDir = _inputDir + "/" + "binGreedy";
		fs::path path(clusterDir);
		if(!fs::exists (path)){
			fs::create_directory(path);
		} 

		string filename_binStats = clusterDir + "/binStats.txt";
		ofstream file_binStats(filename_binStats);

		ofstream fileHifiasmAll(_inputDir + "/binning_results_hifiasm.csv");
		fileHifiasmAll << "Name,Colour" << endl;

		gzFile outputContigFile_min = gzopen(outputFilename.c_str(),"wb");
		gzFile outputContigFile_fasta = gzopen(outputFilename_fasta.c_str(),"wb");

		ofstream file_asmResult = ofstream(_inputDir + "/binning_results.csv");
		file_asmResult << "Name,Colour" << endl;

		//Assembly assembly(_unitigDatas, _contigFeature);

		float prevCutoff = -1;


		unordered_set<u_int32_t> processedNodeNames;
		u_int64_t processedUnitigs = 0;

		//u_int64_t contigIndex = 0;
		u_int64_t __loadState2_index = 0;

		vector<Contig> contigs;
		vector<vector<u_int32_t>> clusters;













		vector<float> allCutoffs;

		for(const SaveState2& saveState : _graph->_cachedGraphStates){
			allCutoffs.push_back(saveState._abundanceCutoff_min);
			cout << saveState._abundanceCutoff_min << endl;
		}
		
		std::reverse(allCutoffs.begin(), allCutoffs.end());



		u_int64_t contigIndexLala = 0;
		unordered_map<u_int32_t, u_int32_t> nodeName_to_contigIndex;

		/*
		vector<u_int32_t> unitigLength_cutoffs = {100000, 50000, 30000, 10000, 5000};
		//vector<vector<UnitigLength>> startingUnitigs;
		//startingUnitigs.resize(unitigLength_cutoffs.size());		

		for(size_t i=0; i<unitigLength_cutoffs.size(); i++){
				
			u_int32_t unitigLength_cutoff_min = unitigLength_cutoffs[i];
			u_int32_t unitigLength_cutoff_max = -1;
			if(i > 0){
				unitigLength_cutoff_max = unitigLength_cutoffs[i-1];
			}

			//vector<UnitigLength>& unitigLengths = _startingUnitigs[i];
			*/

			for(float cutoff : allCutoffs){

				vector<UnitigLength> startingUnitigs;

				_graph->loadState2(cutoff, -1, _unitigDatas);
				_minUnitigAbundance = cutoff / 0.2;


				for(const Unitig& unitig : _graph->_unitigs){
					//cout << unitig._length << " " << unitig._abundance << endl;
					if(unitig._index % 2 == 1) continue;
					/*
					if(unitig._length < unitigLength_cutoff_min) continue;
					if(unitig._length > unitigLength_cutoff_max) continue;
					*/
					if(unitig._abundance < _minUnitigAbundance) continue;
					if(unitig._nbNodes < _kminmerSize*2) continue;
					if(cutoff == 0 && unitig._length < 10000) continue;
					//if(cutoff == 0) continue;

					startingUnitigs.push_back({unitig._length, unitig._abundance, unitig._startNode});
				}

				std::sort(startingUnitigs.begin(), startingUnitigs.end(), UnitigComparator_ByLength2);

				//cout << "Cutoff: " << unitigLength_cutoff_min << "-" << unitigLength_cutoff_max << " " << cutoff << endl;
					
				for(const UnitigLength& unitigLength : startingUnitigs){

					u_int32_t source_nodeIndex = unitigLength._startNodeIndex;
					const Unitig& unitig = _graph->nodeIndex_to_unitig(source_nodeIndex);

					//cout << "Cutoff: " << unitigLength_cutoff_min << "-" << unitigLength_cutoff_max << " " << cutoff << " " << processedUnitigs << " " << startingUnitigs.size() << "     " << unitig._length << endl;
					processedUnitigs += 1;

					/*
					bool isLala = false;
					for(u_int32_t nodeIndex : unitig._nodes){
						u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
						if(nodeName == 289994){
							isLala = true;
						}
					}

					if(isLala){

						unordered_set<u_int32_t> component;
						_graph->getConnectedComponent_unitig(_graph->nodeIndex_to_unitigIndex(BiGraph::nodeName_to_nodeIndex(289994, false)), component);

						unordered_set<u_int32_t> validNodes;
						for (u_int32_t unitigIndex : component){
							for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
								u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
								validNodes.insert(nodeName);
							}
						}
						string outputFilename = _inputDir + "/minimizer_graph_sub_3.gfa";
						GfaParser::rewriteGfa_withoutNodes(_gfaFilename, outputFilename, validNodes, _graph->_isEdgeRemoved, _graph->_graphSuccessors);
			
						cout << component.size() << endl;
						getchar();
					}
					*/

					if(isContigAssembled(unitig._nodes, nodeName_to_contigIndex)) continue;


					u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(unitig._startNode);


					//cout << endl << "Starting unitig: " << BiGraph::nodeIndex_to_nodeName(unitig._startNode) << " " << unitig._length << endl;


					for(u_int32_t nodeIndex : unitig._nodes){
						processedNodeNames.insert(BiGraph::nodeIndex_to_nodeName(nodeIndex));
					}


					vector<u_int32_t> nodePath = unitig._nodes;
					vector<u_int64_t> nodePath_supportingReads;
					//assembly.solveBin2(unitig._startNode, unitig._abundance, _graph, 0, 0, false, nodePath, nodePath_supportingReads, 0);




					string unitigSequenceExpanded;
					//_toBasespace.createSequence(nodePath, unitigSequenceExpanded);
					//cout << "\tExpanded contigs: " << nodePath.size() << endl; //<< " " << unitigSequenceExpanded.size() << endl;

					dereplicateContig2(contigs, nodePath, processedNodeNames, nodeName_to_contigIndex, contigIndexLala);





					/*
					nodePath = contig._nodePath;
					for(size_t i=0; i<_kminmerSize-1 && nodePath.size() > 0; i++){
						nodePath.erase(nodePath.begin());
					}
					for(size_t i=0; i<_kminmerSize-1 && nodePath.size() > 0; i++){
						nodePath.pop_back();
					}
					*/

					/*
					//bool isValid = true;
					Contig contig;// = {nodePath, {}};
					bool isValid = dereplicateContig(contigs, nodePath, contig);


					
					if(isValid && nodePath.size() > _kminmerSize*2){

						contigs.push_back(contig);
						contigIndex += 1;

					}
					*/

					int validContigs = 0;
					for(Contig& contig : contigs){
						if(contig._nodePath.size() > 0){
							validContigs += 1;
						}
					}
					//cout << "\tNb valid contigs: " << validContigs << endl;
					
					/*
					if(validContigs % 100 == 0){

						for(size_t i=0; i<contigs.size(); i++){
							if(contigs[i]._nodePath.size() == 0) continue;
							string unitigSequenceLala;
							_toBasespace.createSequence(contigs[i]._nodePath, unitigSequenceLala);
							cout << "Contig length: " << unitigSequenceLala.size() << " " << contigs[i]._nodePath.size() << endl;
						}

						for(size_t i=0; i<contigs.size(); i++){
							for(size_t j=i+1; j<contigs.size(); j++){
								if(contigs[i]._nodePath.size() == 0 || contigs[j]._nodePath.size() == 0) continue;
								//cout << endl;
								//cout << contigs[i]._nodePath.size() << " " << contigs[i]._nodePath_sorted.size() << endl;
								//cout << contigs[j]._nodePath.size() << " " << contigs[j]._nodePath_sorted.size() << endl;
								u_int64_t nbShared = Utils::computeSharedElements(contigs[i]._nodePath_sorted, contigs[j]._nodePath_sorted);
								if(nbShared == 0) continue;
								cout << nbShared /((float)contigs[i]._nodePath.size()) << " " << nbShared /((float)contigs[j]._nodePath.size()) << endl;
							
								if(nbShared > 0){
								//if(nbShared /((float)contigs[i]._nodePath.size()) > 0.05){

									cout << "Contig size: " << contigs[i]._nodePath.size() << " " << contigs[j]._nodePath.size() << endl;
									unordered_set<u_int32_t> sharedElements;
									Utils::collectSharedElements(contigs[i]._nodePath_sorted, contigs[j]._nodePath_sorted, sharedElements);

									cout << endl;
									for(size_t p=0; p<contigs[i]._nodePath.size(); p++){
										if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(contigs[i]._nodePath[p])) == sharedElements.end()){
											cout << "0";
										}
										else{
											cout << "1";
										}
									}
									cout << endl;

									getchar();
								}
								if(nbShared > 0){
								//if(nbShared /((float)contigs[j]._nodePath.size()) > 0.05){
									
									cout << "Contig size: " << contigs[i]._nodePath.size() << " " << contigs[j]._nodePath.size() << endl;
									
									unordered_set<u_int32_t> sharedElements;
									Utils::collectSharedElements(contigs[i]._nodePath_sorted, contigs[j]._nodePath_sorted, sharedElements);

									cout << endl;
									for(size_t p=0; p<contigs[j]._nodePath.size(); p++){
										if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(contigs[j]._nodePath[p])) == sharedElements.end()){
											cout << "0";
										}
										else{
											cout << "1";
										}
									}
									cout << endl;

									getchar();
								}

							}
						}

					}
					*/

					//cout << "----" << endl;
					
				}







			}


			

		//}








		/*
		for(const Unitig& unitig : _graph->_unitigs){
			//if(unitig._nbNodes <= _kminmerSize*2) continue;
			if(unitig._length < 20000) continue;


			if(unitig._abundance < _minUnitigAbundance) continue;

			startingUnitigs.push_back({unitig._length, unitig._abundance, unitig._startNode});
		}
		*/



		/*
		cout << contigs.size() << endl;
		std::sort(contigs.begin(), contigs.end(), ContigComparator_ByLength);
		//vector<Contig> contigs_reverse = contigs;
		//std::reverse(contigs_reverse.begin(), contigs_reverse.end());

		for(long i=0; i<contigs.size(); i++){
			for(long j=contigs.size()-1; j>=0; j--){
				if(i == j) continue;

				Contig& contig1 = contigs[i];
				Contig& contig2 = contigs[j];

				if(contig1._nodePath.size() == 0) continue;
				if(contig2._nodePath.size() == 0) continue;

				//cout << i << " " << j << " " << contigs.size()-j-1 << endl;

				u_int64_t nbShared = Utils::computeSharedElements(contig1._nodePath_sorted, contig2._nodePath_sorted);
				if(nbShared == 0) continue;

				if(nbShared / ((float)contig2._nodePath.size()) > 0.98){
					contig2._nodePath.clear();
					contig2._nodePath_sorted.clear();
					continue;
				}


				unordered_set<u_int32_t> sharedElements;
				Utils::collectSharedElements(contig1._nodePath_sorted, contig2._nodePath_sorted, sharedElements);


				removeOverlap(contig1._nodePath, contig2._nodePath, sharedElements, false);
				removeOverlap(contig1._nodePath, contig2._nodePath, sharedElements, true);
				
				if(contig2._nodePath.size() > _kminmerSize*2){

					vector<u_int32_t> nodePath_sorted;
					for(u_int32_t nodeIndex : contig2._nodePath){
						nodePath_sorted.push_back(BiGraph::nodeIndex_to_nodeName(nodeIndex));
					}
					std::sort(nodePath_sorted.begin(), nodePath_sorted.end());
					contig2._nodePath_sorted = nodePath_sorted;
				}
				else{
					contig2._nodePath.clear();
					contig2._nodePath_sorted.clear();
				}

			}
		}
		*/

		/*
		cout << "check 1" << endl;
		
		for(size_t i=0; i<contigs.size(); i++){
			for(size_t j=i+1; j<contigs.size(); j++){
				if(contigs[i]._nodePath.size() == 0 || contigs[j]._nodePath.size() == 0) continue;
				//cout << endl;
				//cout << contigs[i]._nodePath.size() << " " << contigs[i]._nodePath_sorted.size() << endl;
				//cout << contigs[j]._nodePath.size() << " " << contigs[j]._nodePath_sorted.size() << endl;
				u_int64_t nbShared = Utils::computeSharedElements(contigs[i]._nodePath_sorted, contigs[j]._nodePath_sorted);
				if(nbShared == 0) continue;
				cout << nbShared /((float)contigs[i]._nodePath.size()) << " " << nbShared /((float)contigs[j]._nodePath.size()) << endl;
			
				//if(nbShared /((float)contigs[i]._nodePath.size()) > 0.05 || nbShared /((float)contigs[j]._nodePath.size()) > 0.05){
				if(nbShared > 0){

					unordered_set<u_int32_t> sharedElements;
					Utils::collectSharedElements(contigs[i]._nodePath_sorted, contigs[j]._nodePath_sorted, sharedElements);

					cout << endl;
					cout << endl;
					cout << endl;
					for(size_t p=0; p<contigs[i]._nodePath.size(); p++){
						if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(contigs[i]._nodePath[p])) == sharedElements.end()){
							cout << "0";
						}
						else{
							cout << "1";
						}
					}
					cout << endl;

					cout << endl;
					for(size_t p=0; p<contigs[j]._nodePath.size(); p++){
						if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(contigs[j]._nodePath[p])) == sharedElements.end()){
							cout << "0";
						}
						else{
							cout << "1";
						}
					}
					cout << endl;

					getchar();
				}

			}
		}
		


		//remove k-1 overlaps

		for(size_t i=0; i<contigs.size(); i++){

			Contig& c = contigs[i];
			if(c._nodePath.size() == 0) continue;

			c._nodePath_sorted.clear();
			for(u_int32_t nodeIndex : c._nodePath){
				c._nodePath_sorted.push_back(BiGraph::nodeIndex_to_nodeName(nodeIndex));
			}
			std::sort(c._nodePath_sorted.begin(), c._nodePath_sorted.end());
		}


		cout << "check 2" << endl;

		for(size_t i=0; i<contigs.size(); i++){
			for(size_t j=i+1; j<contigs.size(); j++){
				if(contigs[i]._nodePath.size() == 0 || contigs[j]._nodePath.size() == 0) continue;
				//cout << endl;
				//cout << contigs[i]._nodePath.size() << " " << contigs[i]._nodePath_sorted.size() << endl;
				//cout << contigs[j]._nodePath.size() << " " << contigs[j]._nodePath_sorted.size() << endl;
				u_int64_t nbShared = Utils::computeSharedElements(contigs[i]._nodePath_sorted, contigs[j]._nodePath_sorted);
				if(nbShared == 0) continue;
				cout << nbShared /((float)contigs[i]._nodePath.size()) << " " << nbShared /((float)contigs[j]._nodePath.size()) << endl;
			
				//if(nbShared /((float)contigs[i]._nodePath.size()) > 0.05 || nbShared /((float)contigs[j]._nodePath.size()) > 0.05){
				if(nbShared > 0){

					unordered_set<u_int32_t> sharedElements;
					Utils::collectSharedElements(contigs[i]._nodePath_sorted, contigs[j]._nodePath_sorted, sharedElements);

					cout << endl;
					cout << endl;
					cout << endl;
					for(size_t p=0; p<contigs[i]._nodePath.size(); p++){
						if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(contigs[i]._nodePath[p])) == sharedElements.end()){
							cout << "0";
						}
						else{
							cout << "1";
						}
					}
					cout << endl;

					cout << endl;
					for(size_t p=0; p<contigs[j]._nodePath.size(); p++){
						if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(contigs[j]._nodePath[p])) == sharedElements.end()){
							cout << "0";
						}
						else{
							cout << "1";
						}
					}
					cout << endl;

				}

			}
		}



		for(size_t i=0; i<contigs.size(); i++){

			Contig& c1 = contigs[i];
			if(c1._nodePath.size() == 0) continue;

			for(size_t j=i+1; j<contigs.size(); j++){
				
				Contig& c2 = contigs[j];
				if(c2._nodePath.size() == 0) continue;

				unordered_set<u_int32_t> sharedElements;
				Utils::collectSharedElements(c1._nodePath_sorted, c2._nodePath_sorted, sharedElements);
				
				u_int64_t overlapSizeLeft = 0;
				u_int64_t overlapSizeRight = 0;

				//left
				for(size_t p=0; p<_kminmerSize-1; p++){
					if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(c2._nodePath[p])) != sharedElements.end()){
						overlapSizeLeft += 1;
					}
					if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(c2._nodePath[c2._nodePath.size() - p - 1])) != sharedElements.end()){
						overlapSizeRight += 1;
					}
				}

				if(overlapSizeLeft == _kminmerSize-1){
					for(size_t p=0; p<_kminmerSize-1; p++){
						u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(c2._nodePath[p]);
						c2._nodePath_sorted.erase(std::remove(c2._nodePath_sorted.begin(), c2._nodePath_sorted.end(), nodeName), c2._nodePath_sorted.end());
					}
					c2._nodePath.erase(c2._nodePath.begin(), c2._nodePath.begin()+_kminmerSize-1);
				}

				if(overlapSizeRight == _kminmerSize-1){
					for(size_t p=0; p<_kminmerSize-1; p++){
						c2._nodePath.pop_back();
						u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(c2._nodePath[c2._nodePath.size() - p - 1]);
						c2._nodePath_sorted.erase(std::remove(c2._nodePath_sorted.begin(), c2._nodePath_sorted.end(), nodeName), c2._nodePath_sorted.end());
					}
				}

			}
		}




		cout << "check 3" << endl;

		for(size_t i=0; i<contigs.size(); i++){
			for(size_t j=i+1; j<contigs.size(); j++){
				if(contigs[i]._nodePath.size() == 0 || contigs[j]._nodePath.size() == 0) continue;
				//cout << endl;
				//cout << contigs[i]._nodePath.size() << " " << contigs[i]._nodePath_sorted.size() << endl;
				//cout << contigs[j]._nodePath.size() << " " << contigs[j]._nodePath_sorted.size() << endl;
				u_int64_t nbShared = Utils::computeSharedElements(contigs[i]._nodePath_sorted, contigs[j]._nodePath_sorted);
				if(nbShared == 0) continue;
				cout << nbShared /((float)contigs[i]._nodePath.size()) << " " << nbShared /((float)contigs[j]._nodePath.size()) << endl;
			
				//if(nbShared /((float)contigs[i]._nodePath.size()) > 0.05 || nbShared /((float)contigs[j]._nodePath.size()) > 0.05){
				if(nbShared > 0){

					unordered_set<u_int32_t> sharedElements;
					Utils::collectSharedElements(contigs[i]._nodePath_sorted, contigs[j]._nodePath_sorted, sharedElements);

					cout << endl;
					cout << endl;
					cout << endl;
					for(size_t p=0; p<contigs[i]._nodePath.size(); p++){
						if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(contigs[i]._nodePath[p])) == sharedElements.end()){
							cout << "0";
						}
						else{
							cout << "1";
						}
					}
					cout << endl;

					cout << endl;
					for(size_t p=0; p<contigs[j]._nodePath.size(); p++){
						if(sharedElements.find(BiGraph::nodeIndex_to_nodeName(contigs[j]._nodePath[p])) == sharedElements.end()){
							cout << "0";
						}
						else{
							cout << "1";
						}
					}
					cout << endl;

					getchar();
				}

			}
		}
		*/

		//ofstream file_contigToNode(outputFilename_fasta + ".nodes.csv");
		unordered_map<u_int32_t, u_int32_t> nodeCounts;

		u_int64_t contigIndex = 0;
		for(Contig& contig : contigs){
			if(contig._nodePath.size() > 0){

				for(u_int32_t nodeIndex : contig._nodePath){
					nodeCounts[BiGraph::nodeIndex_to_nodeName(nodeIndex)] += 1;
				}

				/*
				vector<u_int32_t> nodepath = contig._nodePath;
				if(nodepath.size() > 100){
				nodepath.erase(nodepath.begin(), nodepath.begin()+20);
				nodepath.resize(nodepath.size()-20);
				}
				*/
				
				u_int64_t size = contig._nodePath.size();
				gzwrite(outputContigFile_min, (const char*)&size, sizeof(size));
				gzwrite(outputContigFile_min, (const char*)&contig._nodePath[0], size * sizeof(u_int32_t));
				
				//vector<u_int64_t> supportingReads;
				//getSupportingReads(contig._nodePath, supportingReads);

				//gzwrite(outputContigFile_min, (const char*)&supportingReads[0], size * sizeof(u_int64_t));

				/*
				string unitigSequenceLala;
				_toBasespace.createSequence(contig._nodePath, unitigSequenceLala);

				
				string header = ">ctg" + to_string(contigIndex) + '\n';
				gzwrite(outputContigFile_fasta, (const char*)&header[0], header.size());
				unitigSequenceLala +=  '\n';
				gzwrite(outputContigFile_fasta, (const char*)&unitigSequenceLala[0], unitigSequenceLala.size());

				contigIndex += 1;
				*/
				/*
				continue;

				
				size_t splitSize = 20000;

				if(unitigSequenceLala.size() < splitSize) continue;


				for (size_t i = 0; i < unitigSequenceLala.length(); i += splitSize) {

					string unitigSequence_split = unitigSequenceLala.substr(i, splitSize);
					//cout << str.substr(i, 4) << endl;

					if(unitigSequence_split.size() < splitSize/2) break;

					string header = ">ctg" + to_string(contigIndex) + '\n';
					gzwrite(outputContigFile_fasta, (const char*)&header[0], header.size());
					unitigSequence_split +=  '\n';
					gzwrite(outputContigFile_fasta, (const char*)&unitigSequence_split[0], unitigSequence_split.size());


					//for(u_int32_t nodeIndex : contig._nodePath){
					//	file_contigToNode << BiGraph::nodeIndex_to_nodeName(nodeIndex) << ";" << contigIndex << endl;
					//}

					contigIndex += 1;

				}
				*/
				

			}
		}


		ofstream fileTestLala(outputFilename_fasta + ".nodes2.csv");
		fileTestLala << "Name,Color" << endl;

		for(auto& it: nodeCounts){
			if(it.second == 1){
				fileTestLala << it.first << "," << "blue" << endl;
			}
			else{
				fileTestLala << it.first << "," << "red" << endl;
			}
		}

		fileTestLala.close();


		gzclose(outputContigFile_min);
		gzclose(outputContigFile_fasta);
		fileHifiasmAll.close();
		//file_contigToNode.close();

		//extractContigKminmers(outputFilename_fasta);
		file_asmResult.close();
		
	}

	void getSupportingReads(const vector<u_int32_t>& pathNodes, vector<u_int64_t>& supportingReads){

		//cout << "lala" << endl;
		supportingReads.resize(pathNodes.size());
		return;

		supportingReads.clear();
		vector<u_int32_t> prevNodes;

		for(u_int32_t nodeIndex : pathNodes){

			//cout << nodeIndex << " " << prevNodes.size() <<  endl;
			u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
			UnitigData& unitigData = _unitigDatas[nodeName];

			if(prevNodes.size() == 0){
				
				cout << unitigData._readIndexes.size() << endl;

				prevNodes.push_back(nodeIndex);
				supportingReads.push_back(unitigData._readIndexes[0]);
				continue;
			}

			u_int32_t prevRank = 0;
			u_int64_t supportingRead = -1;

			while(true){
				
				int prevIndex = prevNodes.size() - prevRank - 1;
				//cout << "lalalala     " << prevIndex << endl;
				if(prevIndex < 0 ) break;
				
				u_int32_t prev_nodeIndex = prevNodes[prevIndex];
				u_int32_t prev_nodeName = BiGraph::nodeIndex_to_nodeName(prev_nodeIndex);
				UnitigData& prev_unitigData = _unitigDatas[prev_nodeName];
				
				vector<u_int64_t> sharedReads;
				Utils::collectSharedReads(unitigData, prev_unitigData, sharedReads);
				//cout << "sdfsdfsd     " << sharedReads.size() << endl;
				if(sharedReads.size() > 0){
					prevRank += 1;
					supportingRead = sharedReads[0];

					//cout << "lala: " << supportingRead << endl;
				}
				else{
					break;
				}
			}


			prevNodes.push_back(nodeIndex);
			supportingReads.push_back(supportingRead);

		}
		
		cout << "loulou" << endl;
	}





	/*
	void generateContigs(){
		if(_kminmerSize < 10) return;

		cout << "Generating contigs" << endl;
		
		const string& outputFilename = _inputDir + "/contigs.nodepath.gz";
		gzFile outputContigFile_min = gzopen(outputFilename.c_str(),"wb");

		unordered_set<u_int32_t> writtenUnitigs;

		vector<float> readpathAbudance_values;
		for(const Unitig& u : _graph->_unitigs){
			//if(u._nbNodes < _kminmerSize*2) continue;

			if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._startNode)) != writtenUnitigs.end()) continue;
			if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._endNode)) != writtenUnitigs.end()) continue;

			writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._startNode));
			writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._endNode));

			u_int64_t size = u._nodes.size();
			gzwrite(outputContigFile_min, (const char*)&size, sizeof(size));
			gzwrite(outputContigFile_min, (const char*)&u._nodes[0], size * sizeof(u_int32_t));
		}
		
		gzclose(outputContigFile_min);
	}
	*/

};





#endif
