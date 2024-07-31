

#ifndef MDBG_METAG_READCORRECTION
#define MDBG_METAG_READCORRECTION

#include "../Commons.hpp"
#include "../graph/GraphPOA.hpp"
#include "../utils/spoa64/include/spoa64/spoa.hpp"
#include "../graph/GfaParser.hpp"
/*
Perform Correction:
	- IsSuccessor comme a l'ancienne: parfois on peut avoir des litiges au niveau du successeurs, par exemple entre deux successors, mais souvent l'un de ces deux successeur est le successeur du second ou vis versa 
	- en fait on doit faire de la transitive edge reduction je pense, un node map sur plusieurs successeurs a cause des insertion/deletion

Overlappin read: 2070
Overlappin read: 8958
Overlappin read: 6060
Overlappin read: 8235

./bin/metaMDBG readCorrection ~/appa/run/correction/nanopore_AD_circ1_asm/tmp/ --paf /pasteur/appa/scratch/gbenoit/data/nanopore/subreads/circ1_hifi.fastq.paf.tmp.gz --hifi ~/appa/data/nanopore/subreads/circ1_hifi.fastq -t 32

Illumina:
./bin/metaMDBG readCorrection ~/appa/run/correction/nanopore_AD_circ1_asm/tmp/ -t 32 --illumina ~/appa/data/nanopore/subreads/circ1_illumina.fastq --paf /pasteur/appa/scratch/gbenoit/data/nanopore/subreads/circ1_hifi.fastq.paf.tmp.gz --hifi ~/appa/data/nanopore/subreads/circ1_hifi.fastq
*/


class ReadCorrection : public Tool{
    
public:



	class Edge;

	class Node{
		public:

		u_int64_t _minimizer;
		vector<Edge*> _successors;
		vector<Node*> _predecessors;
		u_int64_t _abundance;

		Node(u_int64_t minimizer){
			_minimizer = minimizer;
			_abundance = 0;
		}

		string toString(){
			return to_string(_minimizer);
		}

	};


	class Edge{

		public:

		Node* _tail;
		Node* _head;
		u_int32_t _weight;

		Edge(Node* tail, Node* head, u_int32_t weight){
			_tail = tail;
			_head = head;
			_weight = weight;
		}
		
	};

	class Graph{
		public:


		bool _needTopoSort;
		unordered_map<u_int64_t, Node*> _nodes;
		vector<Node*> _topoSortNodes;

		Graph(){
			_needTopoSort = true;
		}

		~Graph(){
			for(auto& it : _nodes){

				for(Edge* edge : it.second->_successors){
					delete edge;
				}
				delete it.second;
			}
			
		}

		void addSequence(const vector<u_int64_t>& minimizers){

			if(minimizers.size() < 2) return;

			for(size_t i=0; i<minimizers.size()-1; i++){

				u_int64_t m1 = minimizers[i];
				u_int64_t m2 = minimizers[i+1];

				Node* tail = addNode(m1);
				Node* head = addNode(m2);

				addEdge(tail, head);

			}

			for(size_t i=0; i<minimizers.size(); i++){
				Node* node = addNode(minimizers[i]);
				node->_abundance += 1;
			}
		}

		Node* addNode(u_int64_t minimizer){

			if(_nodes.find(minimizer) != _nodes.end()){
				return _nodes[minimizer];
			}
			//cout << "Add Node: " << unitigIndex << endl;

			Node* node = new Node(minimizer);
			_nodes[minimizer] = node;

			_needTopoSort = true;

			return node;
		}

		void addEdge(Node* tail, Node* head){


			if(tail == nullptr) return;
			if(tail == head) return;

			for(Edge* edge : tail->_successors){
				if(edge->_head == head){
					edge->_weight += 1;
					return;
				}
			}

			Edge* edge = new Edge(tail, head, 1);
			tail->_successors.push_back(edge);
			head->_predecessors.push_back(tail);

			_needTopoSort = true;

			//cout << "Add edge: " << fromNode->_nodeIndex << " -> " << toNode->_nodeIndex<< endl;
		}

		void save(const string& outputFilename, const vector<u_int64_t>& readMinimizers){

			ofstream outputFile(outputFilename);

			ofstream colorFile(outputFilename + ".color.csv");
			colorFile << "Name,Color" << endl;
			
			ofstream edgeFile(outputFilename + ".edge.csv");
			edgeFile << "Name,Edge" << endl;

			//std::cout << "H\tVN:Z:1.0" << std::endl;
			for (const auto& it : _nodes) {

				u_int64_t minimizer = it.first;
				Node* node = it.second;

				string id = node->toString();

				if(std::find(readMinimizers.begin(), readMinimizers.end(), minimizer) != readMinimizers.end()){
					colorFile << id << ",green" << endl; 
				}
				else{
					colorFile << id << ",grey" << endl; 
				}

				//graph.decoder(it->code)
				outputFile << "S\t" << id << "\t" << "*" << "\t" << "LN:i:500" << "\t" << "dp:i:" << node->_abundance << endl;
				//if (is_consensus_node[it->id]) {
				//std::cout << "\tic:Z:true";
				//}
				//std::cout << std::endl;

				string edgeStr = "";

				for (Edge* edge : node->_successors) {
					
					//cout << id << " -> " << edge->_head->toString() << " " << edge->_weight << endl;

					if(edge->_weight < 2) continue;

					outputFile << "L\t" << id << "\t" << "+\t" << edge->_head->toString() << "\t" << "+\t" << "1M" << endl; //\t" << "ew:f:" << jt->weight << endl;


					if(edge->_weight > 1){
						edgeStr += "[" + edge->_head->toString() + " - " + to_string(edge->_weight) + "] ";
					}
				}

				edgeFile << id << "," << edgeStr << endl; 
			}



			outputFile.close();
			colorFile.close();
			edgeFile.close();
		}

		vector<u_int64_t> computeConsensus(const vector<u_int64_t>& readMinimizers){

			vector<u_int64_t> consensus;
			unordered_map<u_int64_t, u_int64_t> nbVisitedTimes;

			float coverage = computeCoverage();
			//cout << "Coverage: " << coverage << endl;

			//cout << _nodes[46601840752821896]->_abundance << " " << max((u_int64_t)1, (u_int64_t)round(_nodes[46601840752821896]->_abundance / coverage)) << endl;
			float minAbundance = coverage * 0.25;

			Node* startNode = nullptr;
			for(u_int64_t m : readMinimizers){
				if(_nodes[m]->_abundance > minAbundance){
					startNode = _nodes[m];
					break;
				}
			}

			if(startNode == nullptr){
				cout << "No valid start node" << endl;
				return consensus;
			}


			consensus.push_back(startNode->_minimizer);

			Node* currentNode = startNode;

			while(true){

				vector<Edge*> successors = getSolidSuccessors(currentNode, nbVisitedTimes, coverage);
				if(successors.size() == 0) break;

				currentNode = successors[0]->_head;
				nbVisitedTimes[currentNode->_minimizer] += 1;

				consensus.push_back(currentNode->_minimizer);
			}

			return consensus;
		}

		vector<Edge*> getSolidSuccessors(Node* node, unordered_map<u_int64_t, u_int64_t>& nbVisitedTimes, float coverage){


			vector<Edge*> solidSuccessors;

			if(node->_successors.size() == 0) return solidSuccessors;

			vector<Edge*> successors;// = node->_successors;
			
			//cout << "\t\tVisit: " << node->_minimizer << endl;
			for(Edge* edge : node->_successors){

				u_int64_t maxVisitables = 1; //max((u_int64_t)1, (u_int64_t)(edge->_head->_abundance / coverage));
				
				//cout << "\t\t\tSuccessor: " << edge->_head->_minimizer << " " << edge->_head->_abundance << " " << maxVisitables << " " << nbVisitedTimes[edge->_head->_minimizer] << endl;
				u_int64_t currentVisitedTimes = nbVisitedTimes[edge->_head->_minimizer];

				//if(currentVisitedTimes)
				if(currentVisitedTimes >= maxVisitables){
					continue;
				}

				successors.push_back(edge);
			}
			



			std::sort(successors.begin(), successors.end(), [](Edge* a, Edge* b){
				return a->_weight > b->_weight;
			});

			float maxWeight = successors[0]->_weight;

			float minWeight = maxWeight * 0.5;

			for(Edge* edge : successors){
				if(edge->_weight < minWeight) continue;
					
				solidSuccessors.push_back(edge);
			}

			return solidSuccessors;
		}

		float computeCoverage(){
			
			vector<float> abundances;

			for (const auto& it : _nodes) {

				u_int64_t minimizer = it.first;
				Node* node = it.second;

				if(node->_abundance <= 2) continue;

				abundances.push_back(node->_abundance);
			}

			float medianAbundance = Utils::compute_median_float(abundances);

			return medianAbundance;
		}
		/*
		Edge* getEdge(Node* fromNode, Node* toNode){

			if(fromNode == nullptr) return nullptr;

			for(Edge* edge : fromNode->_successors){
				if(edge->_toNode == toNode){
					return edge;
				}
			}

			return nullptr;
		}
		*/
		/*
		vector<int> _nodeIDtoIndex;
		vector<int> _nodeIndexToID;
		//unordered_map<u_int32_t, int> _nodeIDtoIndex;
		//unordered_map<int, u_int32_t> _nodeIndexToID;
		
		void computeTopologicalSort(){
			if(!_needTopoSort) return;

			_topoSortNodes.clear();

			stack<int> Stack;
		
			// Mark all the vertices as not visited
			unordered_set<Node*> isVisited;

			//vector<bool> visited(_nodes.size(), false);
			//for (int i = 0; i < _nodes.size(); i++)
			//    visited[i] = false;
		
			// Call the recursive helper function
			// to store Topological
			// Sort starting from all
			// vertices one by one
			//for (int i = 0; i < _nodes.size(); i++){

			for(const auto& it : _nodes){
				Node* node = it.second;
				//if(_nodes[i]->_nodeIndex == -1) continue;
				if(isVisited.find(node) == isVisited.end()){
					topologicalSortUtil(node->_nodeIndex, isVisited, Stack);
				}
			}
		
			// Print contents of stack
			while (!Stack.empty()) {
				_topoSortNodes.push_back(_nodes[Stack.top()]);
				//cout << Stack.top() << " ";
				Stack.pop();
			}
			
			//delete [] visited;


			_nodeIDtoIndex.clear();
			_nodeIndexToID.clear();

			_nodeIDtoIndex.resize(_topoSortNodes.size(), 0);
			_nodeIndexToID.resize(_topoSortNodes.size(), 0);
			//_nodeIDtoIndex[-1] = std::numeric_limits<int>::min();

			for(int i=0; i<_topoSortNodes.size(); i++){
				Node* node = _topoSortNodes[i];
				_nodeIDtoIndex[node->_nodeIndex] = i;
				_nodeIndexToID[i] = node->_nodeIndex;
			}

		}


		void topologicalSortUtil(int v, unordered_set<Node*>& isVisited, stack<int>& Stack){
			// Mark the current node as visited.
			isVisited.insert(_nodes[v]);
		
			for(Edge* edge : _nodes[v]->_successors){
				if(isVisited.find(edge->_toNode) == isVisited.end()){
					topologicalSortUtil(edge->_toNode->_nodeIndex, isVisited, Stack);
				}
			}
			// Recur for all the vertices
			// adjacent to this vertex
			//list<int>::iterator i;
			//for (i = adj[v].begin(); i != adj[v].end(); ++i)
			//    if (!visited[*i])
			//        topologicalSortUtil(*i, visited, Stack);
		
			// Push current vertex to stack
			// which stores result
			Stack.push(v);
		}
		*/
		/*
		void save(const string& outputFilename){

			vector<Node*> nodes;
			for(const auto& it : _nodes){
				nodes.push_back(it.second);
			}

			std::sort(nodes.begin(), nodes.end(), [](Node* a, Node* b){
				return a->_nodeIndex < b->_nodeIndex;
				//return a.size() > b.size();
			});

			//if(fs::exists (outputFilename)) return;

			ofstream outputFile(outputFilename);

			//for(Node* node : _nodes){
			for(Node* node : nodes){
				
				for(Edge* edge : node->_successors){

					//if(getMaxInWeight(edge->_toNode) < 5) continue;
					if(edge->_weight <= 1) continue;
					outputFile << node->_unitigIndex << "(" << node->_nodeIndex  << ")" << "\t" << edge->_toNode->_unitigIndex << "(" << edge->_toNode->_nodeIndex  << ")" << "\t" << edge->_weight << endl; //<< "    " << getMaxInWeight(edge->_toNode) << endl;
				}
			}

			outputFile.close();
		}
		*/

	};

	typedef phmap::parallel_flat_hash_map<KmerVec, vector<u_int32_t>, phmap::priv::hash_default_hash<KmerVec>, phmap::priv::hash_default_eq<KmerVec>, std::allocator<std::pair<KmerVec, vector<u_int32_t>>>, 4, std::mutex> KminmerReadMap;
	typedef phmap::parallel_flat_hash_map<KmerVec, u_int32_t, phmap::priv::hash_default_hash<KmerVec>, phmap::priv::hash_default_eq<KmerVec>, std::allocator<std::pair<KmerVec, u_int32_t>>, 4, std::mutex> KminmerCountMap;


	string _inputFilename;
	string _inputDir;
	string _outputFilename;
	float _minimizerDensity;
    size_t _minimizerSize;
    size_t _kminmerSize;
	string _filename_readMinimizers;
	//bool _isFirstPass;
	int _nbCores;
	string _pafFilename;
	string _hifiFilename;
	string _illuminaFilename;
	string _nanoporeFilename;
	string _filename_exe;
	
	vector<u_int32_t> _allReadSizes;
	u_int64_t _debug_nbMinimizers;

	u_int64_t _minNbMinimizersQuery;
	u_int64_t _minNbMinimizersReference;

    struct ReadWriter{
        u_int64_t _readIndex;
        vector<u_int64_t> _minimizers;
    };

    struct ReadWriter_Comparator {
        bool operator()(ReadWriter const& p1, ReadWriter const& p2){
            return p1._readIndex > p2._readIndex;
        }
    };

	priority_queue<ReadWriter, vector<ReadWriter> , ReadWriter_Comparator> _readWriterQueue;
	u_int64_t _nextReadIndexWriter;
	ofstream _file_readData;
	//unordered_map<u_int64_t, u_int64_t> _minimizerCounts;
	//unordered_map<KmerVec, KminmerData> _kminmersData;
	//gzFile _file_minimizerPos;

	ReadCorrection(): Tool (){
	}


	void parseArgs(int argc, char* argv[]){


		_filename_exe = argv[0];

		args::ArgumentParser parser("readSelection", ""); //"This is a test program.", "This goes after the options."
		args::Positional<std::string> arg_outputDir(parser, "outputDir", "", args::Options::Required);
		args::ValueFlag<string> arg_pafFilename(parser, "", "paf filename", {"paf"}, "");
		args::ValueFlag<string> arg_hifiFilename(parser, "", "hifi filename", {"hifi"}, "");
		args::ValueFlag<string> arg_illuminaFilename(parser, "", "illumina filename", {"illumina"}, "");
		//args::Positional<std::string> arg_outputFilename(parser, "outputFilename", "", args::Options::Required);
		//args::Positional<std::string> arg_inputReadFilename(parser, "inputReadFilename", "", args::Options::Required);
		//args::Positional<std::string> arg_contigs(parser, "contigs", "", args::Options::Required);
		//args::PositionalList<std::string> arg_readFilenames(parser, "reads", "Input filename(s) (separated by space)", args::Options::Required);
		//args::ValueFlag<int> arg_l(parser, "", "Minimizer length", {ARG_MINIMIZER_LENGTH2}, 13);
		//args::ValueFlag<float> arg_d(parser, "", "Minimizer density", {ARG_MINIMIZER_DENSITY2}, 0.005f);
		//args::ValueFlag<std::string> arg_contigs(parser, "", "", {ARG_INPUT_FILENAME_CONTIG}, "");
		args::ValueFlag<int> arg_nbCores(parser, "", "Number of cores", {ARG_NB_CORES2}, NB_CORES_DEFAULT_INT);
		//args::Flag arg_firstPass(parser, "", "Is first pass of multi-k", {ARG_FIRST_PASS});
		//args::Flag arg_isFinalAssembly(parser, "", "Is final multi-k pass", {ARG_FINAL});
		//args::Flag arg_firstPass(parser, "", "Is first pass of multi-k", {ARG_FIRST_PASS});
		args::Flag arg_help(parser, "", "", {'h', "help"}, args::Options::Hidden);
		//args::HelpFlag help(parser, "help", "Display this help menu", {'h'});
		//args::CompletionFlag completion(parser, {"complete"});

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
		//_outputFilename = args::get(arg_outputFilename);
		//_inputFilename = args::get(arg_inputReadFilename);
		_nbCores = args::get(arg_nbCores);

		_pafFilename = "";
		if(arg_pafFilename){
			_pafFilename = args::get(arg_pafFilename);
		}

		_hifiFilename = "";
		if(arg_hifiFilename){
			_hifiFilename = args::get(arg_hifiFilename);
		}

		_illuminaFilename = "";
		if(arg_illuminaFilename){
			_illuminaFilename = args::get(arg_illuminaFilename);
		}

		string filename_parameters = _inputDir + "/parameters.gz";
		gzFile file_parameters = gzopen(filename_parameters.c_str(),"rb");
		gzread(file_parameters, (char*)&_minimizerSize, sizeof(_minimizerSize));
		gzread(file_parameters, (char*)&_kminmerSize, sizeof(_kminmerSize));
		gzread(file_parameters, (char*)&_minimizerDensity, sizeof(_minimizerDensity));
		gzclose(file_parameters);

		openLogFile(_inputDir);

		_logFile << endl;
		_logFile << "Input filename: " << _inputFilename << endl;
		_logFile << "Input dir: " << _inputDir << endl;
		_logFile << "Minimizer length: " << _minimizerSize << endl;
		_logFile << "Kminmer length: " << _kminmerSize << endl;
		_logFile << "Density: " << _minimizerDensity << endl;
		_logFile << endl;

		//_filename_readMinimizers = _outputFilename; //_inputDir + "/read_data.gz";
		_kminmerSize = 2;
		_minNbMinimizersReference = 5;
		_minNbMinimizersQuery = 5;

	
	}


    void execute (){

		
		_nextReadIndexWriter = 0;
		_file_readData = ofstream(_inputDir + "/read_data_corrected.txt");

		ReadParserParallel readParser(_inputDir + "/input.txt", false, false, _nbCores, _logFile);
		_nanoporeFilename = readParser._filenames[0];

		//cout << "ReadCorrection: Setting density to 0.005" << endl << endl<< endl;

		cout << "Attention: actuellement on ne prend que le premier read filename si plusieurs sont fourni" << endl;
		//cout << "Using read filename: " << _usedReadFilename << endl;
		cout << "Index reads: on skip les reads <= 5 minimizers" << endl;
		cout << "reprise: read index 26691: comparer les overlaps minimizer-sapce avec les overlap trouver avec minimap2" << endl;
		//cout << "Indexing read names" << endl;
		//indexReadName();

		cout << "todo: revoir des reads ont un getAlignmentScore vraiment bas alors que la sequene a l'air bonne" << endl;

		cout << "todo: reesayer quality maintenant qu'on les reverse" << endl;
		
		//cout << "Mapping reads" << endl;
		//mapReads();

		//cout << "Loading read overlaps" << endl;

		if(_pafFilename != ""){
			cout << "Loading read overlaps" << endl;
			loadReadOverlaps();
		}

		if(_nanoporeFilename != ""){
			cout << "Indexing Nanopore reads" << endl;

			if(fs::exists(_inputDir + "/read_data_init.txt")){
				cout << "Skip nanopore readSelection" << endl;
			}
			else{
				cout << "Read selection nanopore" << endl;
				string command = _filename_exe + " readSelection " + _inputDir + " " + _inputDir + "/read_data_init.txt" + " " + _inputDir + "/input.txt" + " -t " + to_string(_nbCores);
				Utils::executeCommand(command, _inputDir, _logFile);
			}

		}


		cout << "Indexing mReads" << endl;
		indexReads();

		if(_illuminaFilename != ""){
			cout << "Indexing illumina reads" << endl;

			if(fs::exists(_inputDir + "/read_data_illumina.txt")){
				cout << "Skip illumina readSelection" << endl;
			}
			else{

				ofstream file(_inputDir + "/inputIllumina.txt");
				file << _illuminaFilename << endl;
				file.close();

				cout << "Read selection illumina" << endl;
				string command = _filename_exe + " readSelection " + _inputDir + " " + _inputDir + "/read_data_illumina.txt" + " " + _inputDir + "/inputIllumina.txt" + " -t " + to_string(_nbCores);
				Utils::executeCommand(command, _inputDir, _logFile);
			}

			indexIlluminaReads();
		}
		
		cout << "Loading minimizer reads" << endl;
		loadMinimizerReads(_inputDir + "/read_data_init.txt", false);

		cout << "Correcting mReads" << endl;
		correctReads();

		_file_readData.close();
		
		closeLogFile();
	}


	vector<string>* _fields;
	
	void loadReadOverlaps(){

		if(fs::exists(_inputDir + "/read_data_hifi.txt")){
			cout << "Skip hifi readSelection" << endl;
		}
		else{

			ofstream file(_inputDir + "/inputHifi.txt");
			file << _hifiFilename << endl;
			file.close();

			cout << "Read selection hifi" << endl;
			string command = _filename_exe + " readSelection " + _inputDir + " " + _inputDir + "/read_data_hifi.txt" + " " + _inputDir + "/inputHifi.txt" + " -t " + to_string(_nbCores);
			Utils::executeCommand(command, _inputDir, _logFile);
		}

		cout << "Indexing read name hifi" << endl;
		indexReadName(_hifiFilename);

		cout << "Indexing read name nanopore" << endl;
		indexReadName(_nanoporeFilename);

		_fields = new vector<string>();

		cout << "Loading hifi alignments" << endl;
		loadMinimizerReads(_inputDir + "/read_data_hifi.txt", true);
		 
		PafParser pafParser(_pafFilename);
		auto fp = std::bind(&ReadCorrection::parseAlignmentsGz_read, this, std::placeholders::_1);
		pafParser.parse(fp);

	}


	struct HifiMinimizerReadMatch{
		u_int32_t _readIndex;
		vector<u_int64_t> _minimizers;
		bool _isReversed;
		u_int64_t _score;
	};

	phmap::parallel_flat_hash_map<string, u_int32_t> _readName_to_readIndex;
	unordered_map<u_int32_t, HifiMinimizerReadMatch> _nanoporeReadIndex_to_bestHifiMinimizerRead;

	void parseAlignmentsGz_read(const string& line){

		//cout << line << endl;
		GfaParser::tokenize(line, _fields, '\t');

		const string& readNameNanopore = Utils::shortenHeader((*_fields)[5]);
		const string& readNameHiFi = Utils::shortenHeader((*_fields)[0]);

		//cout << (_readName_to_readIndex.find(readNameNanopore) != _readName_to_readIndex.end()) << endl;
		//cout << (_readName_to_readIndex.find(readNameHiFi) != _readName_to_readIndex.end()) << endl;

		if(_readName_to_readIndex.find(readNameHiFi) == _readName_to_readIndex.end()) return;

		u_int32_t readIndexNanopore = _readName_to_readIndex[readNameNanopore];
		u_int64_t readIndexHiFi = _readName_to_readIndex[readNameHiFi];

		//if(readIndex1 == readIndex2) return;

		u_int64_t alignLength = stoull((*_fields)[10]);

		
		u_int32_t readStart = stoull((*_fields)[2]);
		u_int32_t readEnd = stoull((*_fields)[3]);
		u_int32_t contigLength = stoull((*_fields)[6]);
		u_int32_t contigStart = stoull((*_fields)[7]);
		u_int32_t contigEnd = stoull((*_fields)[8]);
		double length = max((double)(contigEnd - contigStart), (double)(readEnd - readStart));
		double error = 1 - min((double)(contigEnd - contigStart), (double)(readEnd - readStart)) / length;

		//if(error > 0.01) return;
		
		bool isReversed = (*_fields)[4] == "-";

		u_int64_t score = stoull((*_fields)[10]);

		vector<u_int64_t> minimizers = _mReadsHiFi[readIndexHiFi]._minimizers;
		HifiMinimizerReadMatch match = {readIndexHiFi, minimizers, isReversed, score};

		if(_nanoporeReadIndex_to_bestHifiMinimizerRead.find(readIndexNanopore) == _nanoporeReadIndex_to_bestHifiMinimizerRead.end()){
			_nanoporeReadIndex_to_bestHifiMinimizerRead[readIndexNanopore] = match;
		}
		else{
			u_int64_t currentScore = _nanoporeReadIndex_to_bestHifiMinimizerRead[readIndexNanopore]._score;
			if(score > currentScore){
				_nanoporeReadIndex_to_bestHifiMinimizerRead[readIndexNanopore] = match;
			}
		}
		//addReadOverlap(readIndex1, readIndex2, error, isReversed);
		//addReadOverlap(readIndex2, readIndex1, error, isReversed);

		//cout << readName1 << " " << readName2 << " " << readIndex1 << " " << readIndex2 << endl;
		//getchar();
		
		//cout << line << endl;
		//getchar();



	}

	void indexReadName(const string& readFilename){

		//cout << "Indexing read names" << endl;
		//cout << _inputFilename_reads << endl;
		auto fp = std::bind(&ReadCorrection::indexReadName_read, this, std::placeholders::_1);
		ReadParser readParser(readFilename, true, false, _logFile);
		readParser.parse(fp);

	}




	void indexReadName_read(const Read& read){
		_readName_to_readIndex[Utils::shortenHeader(read._header)] = read._index;
	}

	KminmerReadMap _kminmer_to_readIndex;
	KminmerCountMap _illuminaKminmers;
	
	/*
	string _readOverlapsFilename;
	string _usedReadFilename;

	void mapReads(){



		_readOverlapsFilename = _inputDir + "/readOverlaps.paf.gz";

		if(fs::exists(_readOverlapsFilename)){
			cout << "Read overlap file found: " << _readOverlapsFilename << endl;
			return;
		}

		string command = "minimap2 -x ava-ont --dual=yes " + _usedReadFilename + " " + _usedReadFilename + " -t " + to_string(_nbCores) + " | awk '$11>=500' ";
		command += " | gzip -c - > " + _readOverlapsFilename;

		cout << command << endl;
		Utils::executeCommand(command, _inputDir, _logFile);
	}

	struct ReadOverlap{
		u_int32_t _readIndex;
		float _error;
		bool _isReversed;
	};

	void loadReadOverlaps(){

		_fields = new vector<string>();
		_fields_optional = new vector<string>();

		PafParser pafParser(_readOverlapsFilename);
		auto fp = std::bind(&ReadCorrection::parseAlignmentsGz_read, this, std::placeholders::_1);
		pafParser.parse(fp);


		_readName_to_readIndex.clear();
	}

	vector<string>* _fields;
	vector<string>* _fields_optional;

	void indexReadName(){

		//cout << "Indexing read names" << endl;
		//cout << _inputFilename_reads << endl;
		auto fp = std::bind(&ReadCorrection::indexReadName_read, this, std::placeholders::_1);
		ReadParser readParser(_usedReadFilename, true, false, _logFile);
		readParser.parse(fp);

	}

	phmap::parallel_flat_hash_map<string, u_int32_t> _readName_to_readIndex;
	unordered_map<u_int32_t, vector<ReadOverlap>> _readIndex_toOverlappingReadIndexes;

	void indexReadName_read(const Read& read){
		//cout << read._index << " " << read._header << " " << Utils::shortenHeader(read._header) << endl;
		//if(read._index % 100000 == 0) cout << read._index << endl;
		_readName_to_readIndex[Utils::shortenHeader(read._header)] = read._index;

		//getchar();
	}

	void parseAlignmentsGz_read(const string& line){

		//cout << line << endl;
		GfaParser::tokenize(line, _fields, '\t');

		const string& readName1 = Utils::shortenHeader((*_fields)[0]);
		const string& readName2 = Utils::shortenHeader((*_fields)[5]);


		u_int32_t readIndex1 = _readName_to_readIndex[readName1];
		u_int64_t readIndex2 = _readName_to_readIndex[readName2];

		if(readIndex1 == readIndex2) return;

		u_int64_t alignLength = stoull((*_fields)[10]);

		
		u_int32_t readStart = stoull((*_fields)[2]);
		u_int32_t readEnd = stoull((*_fields)[3]);
		u_int32_t contigLength = stoull((*_fields)[6]);
		u_int32_t contigStart = stoull((*_fields)[7]);
		u_int32_t contigEnd = stoull((*_fields)[8]);
		double length = max((double)(contigEnd - contigStart), (double)(readEnd - readStart));
		double error = 1 - min((double)(contigEnd - contigStart), (double)(readEnd - readStart)) / length;

		//if(error > 0.01) return;
		
		bool isReversed = (*_fields)[4] == "-";
		addReadOverlap(readIndex1, readIndex2, error, isReversed);
		addReadOverlap(readIndex2, readIndex1, error, isReversed);

		//cout << readName1 << " " << readName2 << " " << readIndex1 << " " << readIndex2 << endl;
		//getchar();
		
		//cout << line << endl;
		//getchar();



	}

	void addReadOverlap(u_int32_t readIndex1, u_int32_t readIndex2, float error, bool isReversed){


		vector<ReadOverlap>& overlaps = _readIndex_toOverlappingReadIndexes[readIndex1];

		for(const ReadOverlap& overlap : _readIndex_toOverlappingReadIndexes[readIndex1]){
			if(overlap._readIndex == readIndex2) return;
		}

		//if(std::find(overlaps.begin(), overlaps.end(), readIndex2) != overlaps.end()) return;

		overlaps.push_back({readIndex2, error, isReversed});
	}
	*/

	struct MinimizerRead{
		vector<u_int64_t> _minimizers;
		vector<u_int8_t> _qualities;
	};

	//unordered_map<u_int32_t, vector<u_int64_t>> _mReads;
	vector<MinimizerRead> _mReads;
	vector<MinimizerRead> _mReadsHiFi;
	MinimizerReadMap _minimizer_to_readIndex;

	void loadMinimizerReads(const string& filename, bool loadHiFi){
		KminmerParserParallel parser(filename, _minimizerSize, _kminmerSize, false, true, 1);
		parser.parseSequences(LoadMinimizerReadsFunctor(*this, loadHiFi));
	}

	class LoadMinimizerReadsFunctor {

		public:

		ReadCorrection& _parent;
		bool _loadHifi;

		LoadMinimizerReadsFunctor(ReadCorrection& parent, bool loadHifi) : _parent(parent){
			_loadHifi = loadHifi;
		}

		LoadMinimizerReadsFunctor(const LoadMinimizerReadsFunctor& copy) : _parent(copy._parent){
			_loadHifi = copy._loadHifi;
		}

		~LoadMinimizerReadsFunctor(){
		}


		void operator () (const KminmerList& kminmerList) {

			u_int32_t readIndex = kminmerList._readIndex;
			if(readIndex % 100000 == 0) cout << readIndex << endl;

			if(_loadHifi){
				_parent._mReadsHiFi.push_back({kminmerList._readMinimizers, kminmerList._readQualities});
			}
			else{
				_parent._mReads.push_back({kminmerList._readMinimizers, kminmerList._readQualities});
			}

		}
	};

	/*
	void indexReads(){

		KminmerParserParallel parser(_inputDir + "/read_data_all.txt", _minimizerSize, _kminmerSize, false, false, _nbCores);
		parser.parseSequences(IndexReadsFunctor(*this));


	}
	
	class IndexReadsFunctor {

		public:

		ReadCorrection& _parent;

		IndexReadsFunctor(ReadCorrection& parent) : _parent(parent){
		}

		IndexReadsFunctor(const IndexReadsFunctor& copy) : _parent(copy._parent){
		}

		~IndexReadsFunctor(){
		}


		void operator () (const KminmerList& kminmerList) {

			u_int32_t readIndex = kminmerList._readIndex;
			if(readIndex % 100000 == 0) cout << readIndex << endl;
			
			if(kminmerList._readMinimizers.size() <= 4) return; //read too short


			for(u_int64_t minimizer : kminmerList._readMinimizers){
				

				_parent._minimizer_to_readIndex.lazy_emplace_l(minimizer,
				[&readIndex](MinimizerReadMap::value_type& v) { // key exist
					if(std::find(v.second.begin(), v.second.end(), readIndex) == v.second.end()){
						v.second.push_back(readIndex);
					}
				},           
				[&minimizer, &readIndex](const MinimizerReadMap::constructor& ctor) { // key inserted
					
					vector<u_int32_t> readIndexes = {readIndex};

					ctor(minimizer, readIndexes); 

				}); // construct value_type in place when key not present

			}

			
		}
	};
	*/
	
	void indexReads(){

		KminmerParserParallel parser(_inputDir + "/read_data_init.txt", _minimizerSize, _kminmerSize, false, true, _nbCores);
		parser.parse(IndexReadsFunctor(*this));


	}
	
	class IndexReadsFunctor {

		public:

		ReadCorrection& _parent;

		IndexReadsFunctor(ReadCorrection& parent) : _parent(parent){
		}

		IndexReadsFunctor(const IndexReadsFunctor& copy) : _parent(copy._parent){
		}

		~IndexReadsFunctor(){
		}


		void operator () (const KminmerList& kminmerList) {

			if(kminmerList._readMinimizers.size() < _parent._minNbMinimizersQuery) return;

			u_int32_t readIndex = kminmerList._readIndex;
			if(readIndex % 100000 == 0) cout << readIndex << endl;

			for(u_int32_t i=0; i<kminmerList._kminmersInfo.size(); i++){
			
				
				//cout << "------- " << i << endl;

				//_logFile << readIndex << " " << i << endl;
				const ReadKminmerComplete& kminmerInfo = kminmerList._kminmersInfo[i];

				const KmerVec& vec = kminmerInfo._vec;
			

			
				_parent._kminmer_to_readIndex.lazy_emplace_l(vec,
				[&readIndex](KminmerReadMap::value_type& v) { // key exist
					if(std::find(v.second.begin(), v.second.end(), readIndex) == v.second.end()){
						v.second.push_back(readIndex);
					}
				},           
				[&vec, &readIndex](const KminmerReadMap::constructor& ctor) { // key inserted
					
					vector<u_int32_t> readIndexes = {readIndex};

					ctor(vec, readIndexes); 

				}); // construct value_type in place when key not present



			}


			
		}
	};
	
	u_int64_t _illuminaLala;

	void indexIlluminaReads(){

		_illuminaLala = 0;

		KminmerParserParallel parser(_inputDir + "/read_data_illumina.txt", _minimizerSize, 2, false, true, _nbCores);
		parser.parse(IndexIlluminaReadsFunctor(*this));


	}
	

	class IndexIlluminaReadsFunctor {

		public:

		ReadCorrection& _parent;

		IndexIlluminaReadsFunctor(ReadCorrection& parent) : _parent(parent){
		}

		IndexIlluminaReadsFunctor(const IndexIlluminaReadsFunctor& copy) : _parent(copy._parent){
		}

		~IndexIlluminaReadsFunctor(){
		}


		void operator () (const KminmerList& kminmerList) {

			u_int32_t readIndex = kminmerList._readIndex;

			if(readIndex % 100000 == 0){
				cout << readIndex << endl;
				cout << _parent._illuminaLala << " " << readIndex << endl;
			}

			if(kminmerList._kminmersInfo.size() <= 1){


				#pragma omp critical
				{
					_parent._illuminaLala += 1;
				}

			}
			
			//cout << readIndex << " " << kminmerList._kminmersInfo.size() << endl;


			for(u_int32_t i=0; i<kminmerList._kminmersInfo.size(); i++){
			
				
				//cout << "------- " << i << endl;

				//_logFile << readIndex << " " << i << endl;
				const ReadKminmerComplete& kminmerInfo = kminmerList._kminmersInfo[i];

				const KmerVec& vec = kminmerInfo._vec;
			

			
				_parent._illuminaKminmers.lazy_emplace_l(vec,
				[&readIndex](KminmerCountMap::value_type& v) { // key exist
					v.second += 1;
				},           
				[&vec, &readIndex](const KminmerCountMap::constructor& ctor) { // key inserted
					
					ctor(vec, 1); 

				}); // construct value_type in place when key not present



			}


			
		}
	};


	void correctReads(){

		_nbReadWithDuplicatedMinimizers = 0;

		KminmerParserParallel parser2(_inputDir + "/read_data_init.txt", _minimizerSize, _kminmerSize, false, true, 1);
		//parser2.parse(FilterKminmerFunctor2(*this));
		parser2.parseSequences(ReadCorrectionFunctor(*this));

	}

	class ReadCorrectionFunctor {

		public:

		ReadCorrection& _parent;
		std::unique_ptr<spoa64::AlignmentEngine> _alignmentEngine;// = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -4);

		ReadCorrectionFunctor(ReadCorrection& parent) : _parent(parent){
		}

		ReadCorrectionFunctor(const ReadCorrectionFunctor& copy) : _parent(copy._parent){
			_alignmentEngine = spoa64::AlignmentEngine::Create(spoa64::AlignmentType::kNW, 3, -5, -4);
		}

		~ReadCorrectionFunctor(){
		}


		void operator () (const KminmerList& kminmerList) {
			_parent.correctRead(kminmerList, _alignmentEngine);
		}
	};
	
	struct MinimizerReadOverlap{
		vector<u_int64_t> _minimizers;
		float _error;

	};
	void correctRead(const KminmerList& kminmerList, const std::unique_ptr<spoa64::AlignmentEngine>& al){


		_print_debug = true;

		u_int64_t readIndex = kminmerList._readIndex;
		//if(readIndex != 169) return;
		

		//Overlappin read: 2070
		//Overlappin read: 8958
		//Overlappin read: 6060
		//Overlappin read: 8235

		//if(readIndex != 8235) return;
		
		//if(readIndex % 1000 == 0) cout << readIndex << endl;
		//if(readIndex < 100) return;
		

		if(readIndex % 10000 == 0) cout << readIndex << endl;
		if(readIndex < 30) return;

		//#pragma omp critical
		//{
		//	cout << readIndex << endl;
		//}

		





		if(kminmerList._readMinimizers.size() < _minNbMinimizersReference){ //|| kminmerList._readMinimizers.size() > 250
			//cout << "Skip short read" << endl;
			writeRead(readIndex, kminmerList._readMinimizers);
			return;
		}


		if(_print_debug){
			cout << endl << endl;
			cout << "-------------------------------------------------------" << endl;
			cout << "Correcting read: " << readIndex << " " << kminmerList._readMinimizers.size() << endl;
			for(u_int64_t m : kminmerList._readMinimizers){
				cout << m << endl;
			}
			cout << endl << endl;
			//
			//{
			//}
		}

		vector<u_int64_t> correctedReadMinimizers = correctWindow(0, readIndex, kminmerList._readMinimizers, kminmerList._readQualities, _print_debug, al);
		
		/*
		//getchar();
		//return;

		vector<u_int64_t> correctedReadMinimizers;
		std::vector<std::vector<u_int64_t>> windowSequences = {kminmerList._readMinimizers};
		std::vector<std::vector<u_int8_t>> windowQualities = {kminmerList._readQualities};

		int nbSplits = 2;

		while(windowSequences[0].size() > 20){
			windowSequences = splitVector<u_int64_t>(kminmerList._readMinimizers, nbSplits);
			windowQualities = splitVector<u_int8_t>(kminmerList._readQualities, nbSplits);
			nbSplits += 1;
		}
		
		//cout << kminmerList._readMinimizers.size() << " " << windowSequences.size() << endl;


		//for(const vector<u_int64_t>& windowSequence : windowSequences){
		//	cout << windowSequence.size() << endl;
		//}

		for(size_t i=0; i<windowSequences.size(); i++){

			if(print_debug){
				cout << "Correcting window: " << i << endl;
				for(u_int64_t m : windowSequences[i]){
					cout << m << endl;
				}
			}

			correctWindow(i, readIndex, windowSequences[i], windowQualities[i], print_debug, al);

			if(print_debug){
				cout << "done" << endl;
			}

			//getchar();
		}
		*/

		/*
		if(correctedReadMinimizers.size() == 0){ //No correction because not enough coverage
			writeRead(readIndex, readMinimizers);
		}
		else{
			writeRead(readIndex, correctedReadMinimizers);
		}
		*/

		writeRead(readIndex, correctedReadMinimizers);
		//getchar();
		
	}

	struct MinimizerAlignment{
		vector<u_int64_t> _minimizerSequence;
		vector<u_int64_t> _minimizerQualities;
		spoa64::Alignment _alignment;
		double _alignmentScore;
		int _lastMatchPos;
	};

	u_int64_t _nbReadWithDuplicatedMinimizers;

	vector<u_int64_t> correctWindow(size_t windowIndex, u_int64_t readIndex, const vector<u_int64_t>& readMinimizers, const vector<u_int8_t>& readQualities, bool print_debug, const std::unique_ptr<spoa64::AlignmentEngine>& al){

		bool isMinimizerRepeat = false;
		unordered_map<u_int64_t, u_int16_t> minimizerCount;

		unordered_set<u_int64_t> minimizerSet;
		for(u_int64_t m : readMinimizers){
			minimizerSet.insert(m);
			minimizerCount[m] += 1;
		}

		for(const auto& it: minimizerCount){
			if(it.second > 1){

				
				#pragma omp critical
				{
					_nbReadWithDuplicatedMinimizers += 1;

					cout << _nbReadWithDuplicatedMinimizers << " " << readIndex << endl;
					isMinimizerRepeat = true;
				}

				break;
			}
		}

		

		//Kminmers !!
		vector<ReadKminmerComplete> kminmersInfos;
		vector<u_int64_t> minimizersPos(readMinimizers.size(), 0);
		vector<u_int8_t> minimizerQualities(readMinimizers.size(), 0);
		MDBG::getKminmers_complete(_kminmerSize, readMinimizers, minimizersPos, kminmersInfos, readIndex, minimizerQualities);
		

		unordered_map<u_int32_t, u_int32_t> readIndex_to_matchCount;

		for(const ReadKminmerComplete& kminmerInfo : kminmersInfos){

			const KmerVec& vec = kminmerInfo._vec;

			if(_kminmer_to_readIndex.find(vec) == _kminmer_to_readIndex.end()) continue;

			for(u_int32_t rIndex : _kminmer_to_readIndex[vec]){
				if(readIndex == rIndex) continue; //Currently corrected read
				readIndex_to_matchCount[rIndex] += 1;
			}

		}

		vector<u_int32_t> matchingReadIndexes;

		for(const auto& it : readIndex_to_matchCount){

			u_int32_t nbMatches = it.second;
			if(nbMatches < 2) continue;

			//float sim = nbMatches / minimizers.size();
			
			//if(sim < 0.2) continue;

			//cout << "Read index: " << it.first << " " << it.second << endl; 

			matchingReadIndexes.push_back(it.first);
		}

		//cout << "\tNb matching reads: " << matchingReadIndexes.size() << endl;

		//return;

		vector<MinimizerRead> minimizerSequences;
		double nbAlignedMinimizers = 0;

		for(u_int32_t rIndex : matchingReadIndexes){

			const MinimizerRead& matchingReadMinimizers = _mReads[rIndex];

			vector<u_int64_t> minimizerSequence = matchingReadMinimizers._minimizers;

			bool isReversed = isSequenceReversed(readMinimizers, minimizerSequence, al);
			
			if(isReversed){
				std::reverse(minimizerSequence.begin(), minimizerSequence.end());
			}

			int lastMatchPos = 0;
			double alignmentScore = getAlignmentScore(readMinimizers, minimizerSequence, al, lastMatchPos) / (double) minimizerSequence.size();
			
			if(print_debug){
				cout << "\tAl score: " << alignmentScore << endl; //<< " " << getAlignmentScore(minimizerSequence, readMinimizers, al) << endl;
			}

			if(alignmentScore < 0.2) continue;


			//if(readMinimizers.size() > 100) continue;
			//if(computeSimilarity(minimizerSet, matchingReadMinimizers._minimizers) < 0.5) continue;
			//if(nbSharedMinimizers(minimizerSet, matchingReadMinimizers._minimizers) < 10) continue;
			minimizerSequences.push_back({minimizerSequence, {}});

			nbAlignedMinimizers += minimizerSequence.size();

		}

		if(print_debug) cout << "Nb valid alignments: " << minimizerSequences.size() << endl;

		//if(alignments.size() < 10) return readMinimizers;
		
		float alignCoverage = nbAlignedMinimizers / readMinimizers.size();
		if(alignCoverage < 10) return readMinimizers;

		//if(isMinimizerRepeat) return readMinimizers;


		Graph* dbgGraph = new Graph();
		dbgGraph->addSequence(readMinimizers);


		for(size_t i=0; i<minimizerSequences.size(); i++){

			const MinimizerRead& mOverlap = minimizerSequences[i];
			vector<u_int64_t> minimizerSequence = mOverlap._minimizers;

			bool isReversed = isSequenceReversed(readMinimizers, minimizerSequence, al);
			if(isReversed){
				std::reverse(minimizerSequence.begin(), minimizerSequence.end());
			}

			if(_print_debug){
				cout << "\tAdd sequence: " << minimizerSequence.size() << endl;
				for(u_int64_t m : minimizerSequence){
					//cout << "\t\t" << m << endl;
				}
			}

			dbgGraph->addSequence(minimizerSequence);

		}

		if(_print_debug){
			cout << "\tCoverage: " << dbgGraph->computeCoverage() << endl;
		}

		vector<u_int64_t> consensus = dbgGraph->computeConsensus(readMinimizers);

		if(_print_debug){
			dbgGraph->save(_inputDir + "/poaGraph.gfa", readMinimizers);
		}


		vector<u_int64_t> consensusTrimmed = trimCorrectedPath(readMinimizers, consensus, al);


		if(_print_debug){
			cout << endl << "\toriginal sequence: " << endl;
			for(u_int64_t m : readMinimizers){
				cout << "\t\t" << m << endl;
			}

			cout << endl << "\tCorrected sequence: " << endl;
			for(u_int64_t m : consensusTrimmed){
				cout << "\t\t" << m << endl;
			}

			cout << readMinimizers.size() << " " << consensusTrimmed.size() << endl;
			cout << "\tdone" << endl;
			if(isMinimizerRepeat) getchar();
		}

		delete dbgGraph;

		
		return consensusTrimmed;

		//if(isMinimizerRepeat){
			//cout << "Read with repeat" << endl;
			//getchar();
		//}

		//Original method
		/*

		
		unordered_map<u_int32_t, u_int32_t> readIndex_to_matchCount;


		for(u_int64_t minimizer : readMinimizers){

			if(_minimizer_to_readIndex.find(minimizer) == _minimizer_to_readIndex.end()) continue;

			for(u_int32_t rIndex : _minimizer_to_readIndex[minimizer]){
				if(readIndex == rIndex) continue; //Currently corrected read
				readIndex_to_matchCount[rIndex] += 1;
			}
		}

		//cout << readIndex_to_matchCount.size() << endl;
		//cout << "Nb minimizers: " << minimizers.size() << endl;
		//cout << "Total read matches: " << readIndex_to_matchCount.size() << endl;

		vector<u_int32_t> matchingReadIndexes;

		for(const auto& it : readIndex_to_matchCount){

			u_int32_t nbMatches = it.second;
			if(nbMatches < 5) continue;

			//float sim = nbMatches / minimizers.size();
			
			//if(sim < 0.2) continue;

			//cout << "Read index: " << it.first << " " << it.second << endl; 

			matchingReadIndexes.push_back(it.first);
		}

		//cout << "\tNb matching reads: " << matchingReadIndexes.size() << endl;

		//return;

		vector<MinimizerRead> minimizerSequences;

		for(u_int32_t rIndex : matchingReadIndexes){


			const MinimizerRead& matchingReadMinimizers = _mReads[rIndex];
			//if(readMinimizers.size() > 100) continue;
			if(computeSimilarity(minimizerSet, matchingReadMinimizers._minimizers) < 0.5) continue;

			minimizerSequences.push_back(matchingReadMinimizers);


		}
		*/
		

		/* Minimap2 !!!!
		unordered_map<u_int64_t, u_int32_t> minimizer_to_abundance;

		for(const ReadOverlap& overlap : _readIndex_toOverlappingReadIndexes[readIndex]){
			//cout << overlappingReadIndex << endl;

			vector<u_int64_t> readMinimizers = _mReads[overlap._readIndex];
			//cout << "lala: " << overlap._isReversed << endl;
			if(overlap._isReversed){
				//std::reverse(readMinimizers.begin(), readMinimizers.end());
			}
			if(readMinimizers.size() > 100) continue;
			if(computeSimilarity(minimizerSet, readMinimizers) < 0.5) continue;

			minimizerSequences.push_back({readMinimizers, overlap._error});

			for(u_int64_t minimizer : readMinimizers){
				minimizer_to_abundance[minimizer] += 1;
			}

		}
		*/

		/* //Original method
		//getchar();
		vector<float> readAbundances;
		for(u_int64_t minimizer : readMinimizers){
			readAbundances.push_back(minimizer_to_abundance[minimizer]);
		}

		float readAbundance = Utils::compute_median_float(readAbundances);
		if(print_debug) cout << "Read abundance: " << Utils::compute_median_float(readAbundances) << endl;


		if(print_debug) cout << "Nb reads used for correction: " << minimizerSequences.size() << endl;
		//return;

		//getchar();
		//return;

		unordered_map<u_int64_t, u_int32_t> minimizerPos;
		for(size_t i=0; i<readMinimizers.size(); i++){
			u_int64_t m = readMinimizers[i];
			minimizerPos[m] = i;
		}

		std::sort(minimizerSequences.begin(), minimizerSequences.end(), [this, &minimizerSet](const MinimizerReadOverlap & a, const MinimizerReadOverlap & b){
			return computeSimilarity(minimizerSet, a._minimizers) > computeSimilarity(minimizerSet, b._minimizers);
			//return a.size() > b.size();
		});

		//std::sort(minimizerSequences.begin(), minimizerSequences.end(), [this, &minimizerSet](const MinimizerReadOverlap & a, const MinimizerReadOverlap & b){
		//	return a._error < b._error;
			//return a.size() > b.size();
		//});
		*/

		/*
		GraphPOA::GraphPOA* graphPOA = new GraphPOA::GraphPOA(kminmerList._readMinimizers);

		int i = 0;
		for(vector<u_int64_t> minimizerSequence: minimizerSequences){

			bool isReversed = isSequenceReversed(kminmerList._readMinimizers, minimizerSequence);

			//cout << "\tAdd sequence: " << i << " " << minimizerSequence.size() << " " << graphPOA->_graph->_nodes.size() << endl;

			
			if(print_debug){
				cout << endl << endl;
				cout << "\tAdd sequence: " << i << " " << minimizerSequence.size() << " " << graphPOA->_graph->_nodes.size() << endl;


				for(u_int64_t m : minimizerSequence){
					cout << "\t" << m << endl;
				}
				
				cout << "\tIs reversed: " << isReversed << endl;
			}
			

			if(graphPOA->_graph->_nodes.size() > 5000){
				//cout << "Complex read skipped" << endl;
				//return;
			}

			if(isReversed){
				std::reverse(minimizerSequence.begin(), minimizerSequence.end());
			}
			//cout << "Similarity: " << computeSimilarity(minimizerSet, minimizerSequence);

			graphPOA->addSequence(minimizerSequence);
			//graphPOA->_graph->save(_outputDir + "/poaGraph.txt");

			i += 1;
			
			//cout << "\tNb nodes: " << graphPOA->_graph->_nodes.size() << endl;
			//getchar();
		}

		//graphPOA->_graph->save(_inputDir + "/poaGraph.txt");

		graphPOA->_graph->filterGraph();
		

		//graphPOA->_graph->save(_inputDir + "/poaGraph_filtered.txt");

		vector<u_int64_t> readMinimizersCorrected = graphPOA->performCorrection(kminmerList._readMinimizers);

		if(print_debug){
			for(u_int64_t minimizer : readMinimizersCorrected){
				cout << minimizer << endl;
			}

			cout << "done" << endl;
			//getchar();
		}
		//if(readIndex == 25) getchar();

		if(readMinimizersCorrected.size() == 0){ //No correction because not enough coverage
			writeRead(readIndex, kminmerList._readMinimizers);
		}
		else{
			writeRead(readIndex, readMinimizersCorrected);
		}
		
		delete graphPOA;

		//exit(1);
		*/

		/*
		if(print_debug) cout << "Nb reads used for correction: " << minimizerSequences.size() << endl;




		//std::unique_ptr<spoa::AlignmentEngine> alignmentEngine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -4);
		spoa64::Graph graph{};
		
		//string readQuality = 
		//if()
		//string backboneQuality = "";
		//for(size_t i=0; i<contigOriginalSequence.size(); i++){
		//	backboneQuality += '!';
		//}

		vector<u_int64_t> weights(readQualities.size(), 0);
		for(size_t i=0; i<readQualities.size(); i++){
			weights[i] = readQualities[i];
		}

		vector<u_int64_t> readNodes;
		graph.AddAlignment(spoa64::Alignment(), readMinimizers, readMinimizers.size(), weights);
		
		
		for (const auto& it : graph.nodes()) {
			readNodes.push_back(it->id);
		}




		//vector<u_int64_t> lala;// = readMinimizers;
		//lala.push_back(1);
		//lala.erase(lala.end()-1);
		//lala.erase(lala.end()-1);
		//lala.push_back(1);
		//lala.push_back(1);
		//getAlignmentScore(readMinimizers, lala, al);
		//getchar();

		//vector<u_int64_t> lala = kminmerList._readMinimizers;

		//lala.insert(lala.begin()+5, 5);
		//lala.erase(lala.begin()+10);
		//lala[15] = 18;
		//lala[15] = 18;
		//minimizerSequences.insert(minimizerSequences.begin(), lala);
		//minimizerSequences.insert(minimizerSequences.begin(), lala);
		//minimizerSequences.insert(minimizerSequences.begin(), lala);
		//minimizerSequences.insert(minimizerSequences.begin(), lala);
		//minimizerSequences.insert(minimizerSequences.begin(), lala);

		//getAlignmentScore(kminmerList._readMinimizers, lala);

		//exit(1);

		double nbAlignedMinimizers = 0;

		vector<MinimizerAlignment> alignments;

		for(size_t i=0; i<minimizerSequences.size(); i++){

			const MinimizerRead& mOverlap = minimizerSequences[i];

			vector<u_int64_t> minimizerSequence = mOverlap._minimizers;

			if(print_debug){
				cout << endl << "\tAlign sequence: " << i << " " << minimizerSequence.size() << " " << graph.nodes_.size() << endl; //<< " " << graphPOA->_graph->_nodes.size() << endl;
				for(u_int64_t m : minimizerSequence){
					//cout << "\t" << m << endl; 
				}
			
			}

			bool isReversed = isSequenceReversed(readMinimizers, minimizerSequence, al);
			//isReversed = !isReversed;

			//if(print_debug) cout << "\tIs reversed: " << isReversed << endl;

			//if(print_debug){
			//	double alignmentScore = getAlignmentScore(readMinimizers, minimizerSequence, al);
			//	cout << "\tAl score: " << alignmentScore << endl;
			//}
			
			if(isReversed){
				std::reverse(minimizerSequence.begin(), minimizerSequence.end());
			}

			//if(print_debug){
			//}

			int lastMatchPos = 0;
			double alignmentScore = getAlignmentScore(readMinimizers, minimizerSequence, al, lastMatchPos) / (double) minimizerSequence.size();
			
			if(print_debug){
				cout << "\tAl score: " << alignmentScore << endl; //<< " " << getAlignmentScore(minimizerSequence, readMinimizers, al) << endl;
			}

			if(alignmentScore < 0.2) continue;
			//for(u_int64_t m : minimizerSequence._minimizers){
			//	cout << "\t" << m << endl;
			//}

			vector<u_int64_t> qualities(mOverlap._qualities.size(), 0);
			for(size_t i=0; i<mOverlap._qualities.size(); i++){
				qualities[i] = mOverlap._qualities[i];
			}


			if(isReversed){
				std::reverse(qualities.begin(), qualities.end());
			}

			spoa64::Alignment alignment = al->Align(
							minimizerSequence,
							minimizerSequence.size(),
							graph);

			alignments.push_back({minimizerSequence, qualities, alignment, lastMatchPos});

			nbAlignedMinimizers += minimizerSequence.size();
		}


		if(print_debug) cout << "Nb valid alignments: " << alignments.size() << endl;

		//if(alignments.size() < 10) return readMinimizers;
		
		float alignCoverage = nbAlignedMinimizers / readMinimizers.size();
		if(alignCoverage < 10) return readMinimizers;

		std::sort(alignments.begin(), alignments.end(), [](const MinimizerAlignment & a, const MinimizerAlignment & b){
			return a._alignmentScore > b._alignmentScore;
			//return a._lastMatchPos < b._lastMatchPos;
		});


		for(size_t i=0; i<alignments.size(); i++){

			MinimizerAlignment& al = alignments[i];


			if(print_debug){
				cout << endl << "\tAdd alignment: " << i << " " << al._minimizerSequence.size() << " " << al._alignmentScore << " " << graph.nodes_.size() << endl; //<< " " << graphPOA->_graph->_nodes.size() << endl;
				for(u_int64_t m : al._minimizerSequence){
					//cout << "\t" << m << endl; 
				}
			
			}


			graph.AddAlignment(al._alignment, al._minimizerSequence, al._minimizerSequence.size(), al._minimizerQualities);

		}
		*/

		/*
		if(print_debug){

			//cout << endl;
			//for(u_int64_t nodeId : readNodes){
			//	cout << graph.nodes_[nodeId]->Coverage() << endl;
				//for(spoa64::Graph::Edge* edge : graph.nodes_[nodeId]->outedges){
				//	cout << "\t"
				//}
			//}

			cout << endl;
			for(size_t i=0; i<readMinimizers.size(); i++){
				cout << readMinimizers[i] << " " << graph.nodes_[i]->Coverage() << " "  << endl;
			}
		}
		*/

		//vector<u_int64_t> correctedReadMinimizers = performCorrection(graph, kminmerList._readMinimizers, al, readNodes);
		
		//cout << endl << "Corrected read: " << endl;
		//for(u_int64_t m : correctedReadMinimizers){ 
		//	cout << "\t" << m << endl;
		//}
		//cout << endl;

		//getchar();
		



		//if(print_debug) cout << "\tCleaning graph" << endl;

		//spoa64::Graph graphCleaned{};
		//cleanGraph(graph, graphCleaned);

		//if(print_debug) cout << "\tCleaning done" << endl;

		//if(graphCleaned.nodes_.size() < readMinimizers.size()) return readMinimizers;
		/*
 		std::vector<u_int64_t> correctedCoverages;
		//std::reverse(correctedSequence.begin(), correctedSequence.end());
		vector<u_int64_t> correctedReadMinimizers = performCorrection(graph, readMinimizers, al, correctedCoverages, minimizerSet);
		*/

		/*
		std::vector<u_int64_t> coverages;
		std::vector<u_int64_t> correctedReadMinimizers = graph.GenerateConsensusVec(coverages);

		double covSum = 0;// (nbAddedSequences.size()) / 2;
		double covN = 0;

		for(size_t i=0; i<coverages.size(); i++){
			//cout << correctedReadMinimizers[i] << " " << coverages[i] << endl;

			covSum += coverages[i];
			covN += 1;
		}

		double average_coverage = covSum / covN;
		average_coverage *= 0.75;

		int32_t begin = 0, end = correctedReadMinimizers.size() - 1;
		for (; begin < static_cast<int32_t>(correctedReadMinimizers.size()); ++begin) {
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
			correctedReadMinimizers = vector<u_int64_t>(&correctedReadMinimizers[begin], &correctedReadMinimizers[end]);
			//correctedSequence = correctedSequence.substr(begin, end - begin + 1);
		}
		*/

		//cout << endl << "Original reads:" << endl;
		//for(size_t i=0; i<readMinimizers.size(); i++){
		//	cout << readMinimizers[i] << endl;
		//}

		//cout << endl << "Corrected reads:" << endl;
		//for(size_t i=0; i<correctedReadMinimizers.size(); i++){
		//	cout << correctedReadMinimizers[i] << endl;
		//}
		/*
		if(print_debug){

			

			cout << endl << "\tOriginal sequence:" << endl;
			for(size_t i=0; i<readMinimizers.size(); i++){
				cout << "\t" << readMinimizers[i] << " " << ((int)readQualities[i]) << endl;
			}

			cout << endl << "\tCorrected sequence OLD:" << endl;
			for(size_t i=0; i<correctedReadMinimizers.size(); i++){
				cout << "\t" << correctedReadMinimizers[i] << endl;//<< " " << correctedCoverages[i] << endl;
			}

			if(_illuminaFilename != ""){

				cout << endl << "\tOriginal sequence:" << endl;
				for(size_t i=0; i<readMinimizers.size()-1; i++){

					u_int64_t m1 = readMinimizers[i];
					u_int64_t m2 = readMinimizers[i+1];
					KmerVec vec;
					vec._kmers = {m1, m2};
					vec = vec.normalize();

					if(_illuminaKminmers.find(vec) == _illuminaKminmers.end()){
						cout << m1 << " 0" << endl; 
					}
					else{
						cout << m1 << " " << _illuminaKminmers[vec] << endl; 
					}
				}

				cout << endl << "\tCorrected sequence OLD:" << endl;
				if(correctedReadMinimizers.size() > 0){
					for(size_t i=0; i<correctedReadMinimizers.size()-1; i++){

						u_int64_t m1 = correctedReadMinimizers[i];
						u_int64_t m2 = correctedReadMinimizers[i+1];
						KmerVec vec;
						vec._kmers = {m1, m2};
						vec = vec.normalize();

						if(_illuminaKminmers.find(vec) == _illuminaKminmers.end()){
							cout << m1 << " 0" << endl; 
						}
						else{
							cout << m1 << " " << _illuminaKminmers[vec] << endl; 
						}
					}
				}

			}
			
			if(_hifiFilename != ""){

				if(_nanoporeReadIndex_to_bestHifiMinimizerRead.find(readIndex) != _nanoporeReadIndex_to_bestHifiMinimizerRead.end()){
					cout << endl << "\tHifi sequence: " << _nanoporeReadIndex_to_bestHifiMinimizerRead[readIndex]._isReversed << endl;
					
					vector<u_int64_t> hifiMinimizers = _nanoporeReadIndex_to_bestHifiMinimizerRead[readIndex]._minimizers;
					if(_nanoporeReadIndex_to_bestHifiMinimizerRead[readIndex]._isReversed){
						std::reverse(hifiMinimizers.begin(), hifiMinimizers.end());
					}


					if(_illuminaFilename != ""){
						for(size_t i=0; i<hifiMinimizers.size()-1; i++){

							u_int64_t m1 = hifiMinimizers[i];
							u_int64_t m2 = hifiMinimizers[i+1];
							KmerVec vec;
							vec._kmers = {m1, m2};
							vec = vec.normalize();

							if(_illuminaKminmers.find(vec) == _illuminaKminmers.end()){
								cout << m1 << " 0" << endl; 
							}
							else{
								cout << m1 << " " << _illuminaKminmers[vec] << endl; 
							}
						}
					}
					else{
						for(u_int64_t m : hifiMinimizers){
							cout << "\t" << m << endl;
						}
					}

					cout << endl << "Score: " << endl;
					printAlignmentScore(hifiMinimizers, correctedReadMinimizers, al);

					cout << readIndex << " " << _nanoporeReadIndex_to_bestHifiMinimizerRead[readIndex]._readIndex << endl;

				
					
					savePoaGraph(graph, _inputDir + "/poaGraph.gfa", readNodes);
					//savePoaGraph(graphCleaned, _inputDir + "/poaGraph.gfa.cleaned", readNodes);
					cout << "Original sequence size: " << readMinimizers.size() << endl;
					cout << "Corrected sequence size: " << correctedReadMinimizers.size() << endl;
					cout << "Read index: " << readIndex << endl;
					
					getchar();
				}

			}
			

		}
		*/
		/*
		#pragma omp critical
		{
			cout << "Read index: " << readIndex << endl;
			cout << "Original sequence size: " << readMinimizers.size() << endl;
			cout << "Corrected sequence size: " << correctedReadMinimizers.size() << endl;
		}
		*/

		//getchar();
		/*
		cout << "Spoa path: " << endl;
		std::vector<u_int64_t> coverages;
		std::vector<u_int64_t> correctedSequence = graph.GenerateConsensusVec(coverages);

		cout << correctedSequence.size() << " " << coverages.size() << endl;
		for(size_t i=0; i<correctedSequence.size(); i++){ 
			u_int64_t minimizer = correctedSequence[i];
			cout << "\t" << minimizer << " " << coverages[i] << endl;
		}
		cout << endl;

		GraphPOA::GraphPOA* graphPOA_1 = new GraphPOA::GraphPOA(correctedSequence);
		graphPOA_1->getAlignment(kminmerList._readMinimizers);

		vector<u_int64_t> s = kminmerList._readMinimizers;
		std::reverse(s.begin(), s.end());
		GraphPOA::GraphPOA* graphPOA_rev = new GraphPOA::GraphPOA(correctedSequence);
		graphPOA_rev->getAlignment(s);
		*/


		
		//getchar();

		//return correctedReadMinimizers;
	}

	template<typename T>
	std::vector<std::vector<T>> splitVector(const std::vector<T>& vec, size_t n)
	{
		std::vector<std::vector<T>> outVec;

		size_t length = vec.size() / n;
		size_t remain = vec.size() % n;

		size_t begin = 0;
		size_t end = 0;

		for (size_t i = 0; i < std::min(n, vec.size()); ++i)
		{
			end += (remain > 0) ? (length + !!(remain--)) : length;

			outVec.push_back(std::vector<T>(vec.begin() + begin, vec.begin() + end));

			begin = end;
		}

		return outVec;
	}

	//void findMostSupportedPath(GraphPOA::GraphPOA* graphPOA, u_int64_t startMinimizer, u_int64_t endMinimizer){

		//unordered_set<u_int32_t> startingNodeIndex;
		//for()
	//}


	float computeSimilarity(const unordered_set<u_int64_t>& s1, const vector<u_int64_t>& s2){

		double nbSharedElements = 0;

		for(u_int64_t m : s2){
			if(s1.find(m) == s1.end()) continue;
			nbSharedElements += 1;
		}

		return nbSharedElements / s2.size();
	}

	float nbSharedMinimizers(const unordered_set<u_int64_t>& s1, const vector<u_int64_t>& s2){

		double nbSharedElements = 0;

		for(u_int64_t m : s2){
			if(s1.find(m) == s1.end()) continue;
			nbSharedElements += 1;
		}

		return nbSharedElements;
	}

	/*
	bool isSequenceReversed(unordered_map<u_int64_t, u_int32_t>& minimizerPos, const vector<u_int64_t>& s2){

		int nbCorrect = 0;
		int nbReversed = 0;

		for(long i=0; i<s2.size()-1; i++){

			u_int64_t m1 = s2[i];
			u_int64_t m2 = s2[i+1];

			if(minimizerPos.find(m1) == minimizerPos.end()) continue;
			if(minimizerPos.find(m2) == minimizerPos.end()) continue;

			if(minimizerPos[m1] < minimizerPos[m2]){
				nbCorrect += 1;
			}
			else{
				nbReversed += 1;
			}
		}

		return nbReversed > nbCorrect;

	}
	*/

	double getAlignmentScore(const vector<u_int64_t>& s1, const vector<u_int64_t>& s2, const std::unique_ptr<spoa64::AlignmentEngine>& al, int& lastMatch){
		
		vector<u_int64_t> weights(s1.size(), 1);
		//cout << "-" << endl;
		spoa64::Graph graph{};

		graph.AddAlignment(spoa64::Alignment(), s1, s1.size(), weights);
		spoa64::Alignment alignment = al->Align(s2, s2.size(), graph);

		double score = 0;

		//cout << "Size: " << s1.size() << " " << s2.size() << endl;
		//for (const auto& it : alignment) {
		//	cout << it.first << " " << it.second << endl;
		//}

		bool isOpening = true;
		int nbDeletions = 0;

		for (size_t i=0; i<alignment.size(); i++) {

			int64_t v1 = alignment[i].first;
			int64_t v2 = alignment[i].second;
			
			if(v1 == -1){ //insert in
				if(isOpening){
					//cout << "insertion Opening" << endl;
				}
				else{
					//cout << "insertion" << endl;
					//nbDeletions += 1;
				}

				//score -= nbDeletions;
				//nbDeletions = 0;
				//score -= 1;
			}
			else if(v2 == -1){ //insert in
				if(isOpening){
					//cout << "deletion Opening" << endl;
				}
				else{
					//cout << "deletion" << endl;
					//nbDeletions += 1;
				}
				//score -= 1;
			}
			else if(s1[v1] == s2[v2]){
				lastMatch = i;
				isOpening = false;
				//cout << "match" << endl;
				//score -= nbDeletions;
				//nbDeletions = 0;
				score += 1;
			}
			else{
				isOpening = false;
				//cout << "missmatch" << endl;
				//score -= nbDeletions;
				//nbDeletions = 0;
				score -= 1;
			}
		}

		//cout << "done" << endl;

		return score;
	}

	void printAlignmentScore(const vector<u_int64_t>& s1, const vector<u_int64_t>& s2, const std::unique_ptr<spoa64::AlignmentEngine>& al){
		
		vector<u_int64_t> weights(s1.size(), 1);
		//cout << "-" << endl;
		spoa64::Graph graph{};

		graph.AddAlignment(spoa64::Alignment(), s1, s1.size(), weights);
		spoa64::Alignment alignment = al->Align(s2, s2.size(), graph);

		for (const auto& it : alignment) {
			int64_t v1 = it.first;
			int64_t v2 = it.second;
			
			if(v1 == -1){ //insert in
				cout << "insertion" << endl;
			}
			else if(v2 == -1){ //insert in
				cout << "deletion" << endl;
			}
			else if(s1[v1] == s2[v2]){
				cout << "match" << endl;
			}
			else{
				cout << "missmatch" << endl;
			}
		}

	}

	/*
	double getAlignmentScoreNorm(const vector<u_int64_t>& s1, const vector<u_int64_t>& s2, const std::unique_ptr<spoa64::AlignmentEngine>& al){
		
		spoa64::Graph graph{};

		graph.AddAlignment(spoa64::Alignment(), s1, s1.size());
		spoa64::Alignment alignment = al->Align(s2, s2.size(), graph);

		double score = 0;

		//cout << "Size: " << s1.size() << " " << s2.size() << endl;
		//for (const auto& it : alignment) {
		//	cout << it.first << " " << it.second << endl;
		//}
		double n = 0;

		for (const auto& it : alignment) {
			int64_t v1 = it.first;
			int64_t v2 = it.second;
			
			if(v1 == -1){ //insert in
				//cout << "insertion" << endl;
				//score -= 1;
			}
			else if(v2 == -1){ //insert in
				//cout << "deletion" << endl;
				//score -= 1;
			}
			else if(s1[v1] == s2[v2]){
				//cout << "match" << endl;
				score += 1;
				n += 1;
			}
			else{
				//cout << "missmatch" << endl;
				n += 1;
				//score -= 1;
			}
		}

		//cout << "done" << endl;

		if(n == 0) return 0;
		return score / n;
	}
	*/

	bool isSequenceReversed(const vector<u_int64_t>& s1, const vector<u_int64_t>& s2, const std::unique_ptr<spoa64::AlignmentEngine>& al){

		int lastMatchpos = 0;
		double score_forward = getAlignmentScore(s1, s2, al, lastMatchpos);

		vector<u_int64_t> s2_rev = s2; 
		std::reverse(s2_rev.begin(), s2_rev.end());
		double score_reverse = getAlignmentScore(s1, s2_rev, al, lastMatchpos);

		//cout << score_forward << " " << score_reverse << endl;

		return score_reverse > score_forward;
		/*
		vector<u_int64_t> s2_rev = s2;
		std::reverse(s2_rev.begin(), s2_rev.end());

		std::unique_ptr<spoa::AlignmentEngine> alignmentEngine2 = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -4);
		
		spoa::Graph graph2{};
		graph2.AddAlignment(spoa::Alignment(), s1, s1.size());
		spoa::Alignment alignment2 = alignmentEngine2->Align(s2_rev, s2_rev.size(), graph2);

		GraphPOA::GraphPOA* graphPOA_1 = new GraphPOA::GraphPOA(s1);
		graphPOA_1->addSequence(s2);

		vector<u_int64_t> s2_rev = s2;
		std::reverse(s2_rev.begin(), s2_rev.end());

		GraphPOA::GraphPOA* graphPOA_2 = new GraphPOA::GraphPOA(s1);
		graphPOA_2->addSequence(s2_rev);

		u_int32_t size1 = graphPOA_1->_scores.get(graphPOA_1->_scores._nbRows-1, graphPOA_1->_scores._nbCols-1);
		u_int32_t size2 = graphPOA_2->_scores.get(graphPOA_2->_scores._nbRows-1, graphPOA_2->_scores._nbCols-1);

		//cout <<  << endl;
		//cout <<  << endl;
		//cout << size1 << " " << size2 << endl;
		delete graphPOA_1;
		delete graphPOA_2;

		return size2 > size1;
		*/
	}

	/*
	vector<u_int64_t> performCorrection(spoa64::Graph& graph, const vector<u_int64_t>& readMinimizers, const std::unique_ptr<spoa64::AlignmentEngine>& al, const vector<u_int64_t>& readNodes){
	
		bool print_debug = false;

		vector<float> coverages;

		for (u_int64_t nodeId : readNodes) {

			float coverage = graph.nodes_[nodeId]->Coverage();
			if(coverage <= 2) continue;

			coverages.push_back(coverage);
		}

		float readCoverage = Utils::compute_median_float(coverages);
		if(print_debug) cout << "Read coverage: " << readCoverage << endl;

		float minCoverage = readCoverage * 0.5;

		vector<u_int64_t> solidReadNodes;
		for (u_int64_t nodeId : readNodes) {

			float coverage = graph.nodes_[nodeId]->Coverage();
			if(coverage < minCoverage) continue;

			solidReadNodes.push_back(nodeId);
		}
		

		if(solidReadNodes.size() == 0){
			if(print_debug) cout << "No solid nodes" << endl;
			
			return {};
		}

		if(print_debug){	
			cout << endl << "Solid read nodes: " << endl;
			for (u_int64_t nodeId : solidReadNodes) {
				cout << nodeId << " " << graph.decoder(graph.nodes_[nodeId]->code) << " " << graph.nodes_[nodeId]->Coverage() << endl;
			}
		}

		vector<u_int64_t> nodePath;


		//Extend left
		vector<u_int64_t> nodePathLeft = extendNode(graph, graph.nodes_[solidReadNodes[0]].get(), readCoverage, true);
		std::reverse(nodePathLeft.begin(), nodePathLeft.end());

		for(u_int64_t nodeId : nodePathLeft){
			//cout << nodeId << " " << graph.decoder(graph.nodes_[nodeId]->code) << " " << graph.nodes_[nodeId]->Coverage() << endl;
			nodePath.push_back(nodeId);
		}

		//Extend between solid nodes
		for(size_t i=0; i<solidReadNodes.size()-1; i++){

			if(print_debug) cout << "Search path: " << solidReadNodes[i] << " -> " << solidReadNodes[i+1] << endl;

			vector<u_int64_t> path = computePath(graph, solidReadNodes[i], solidReadNodes[i+1], -1);
			if(print_debug) cout << "\t" << path.size() << endl;

			nodePath.push_back(solidReadNodes[i]);

			if(path.size() < 2) continue;

			for(size_t i=1; i<path.size()-1; i++){
				nodePath.push_back(path[i]);
			}
			//getchar();
		}
		nodePath.push_back(solidReadNodes[solidReadNodes.size()-1]);


		//Extend right
		vector<u_int64_t> nodePathRight = extendNode(graph, graph.nodes_[solidReadNodes[solidReadNodes.size()-1]].get(), readCoverage*0.75, false);
		for(u_int64_t nodeId : nodePathRight){
			//cout << nodeId << " " << graph.decoder(graph.nodes_[nodeId]->code) << " " << graph.nodes_[nodeId]->Coverage() << endl;
			nodePath.push_back(nodeId);
		}
		



		if(print_debug) cout << endl << "Corrected read: " << endl;
		vector<u_int64_t> correctReadMinimizers;
		for(u_int64_t nodeId : nodePath){
			if(print_debug) cout << nodeId << " " << graph.decoder(graph.nodes_[nodeId]->code) << " " << graph.nodes_[nodeId]->Coverage() << endl;
			correctReadMinimizers.push_back(graph.decoder(graph.nodes_[nodeId]->code));
		}

		//spoa64::Alignment alignment = al->Align(readMinimizers, readMinimizers.size(), graph);

		//graph.performCorrection(alignment, readMinimizers, readMinimizers.size(), 1);

		return correctReadMinimizers;

	}

	vector<u_int64_t> extendNode(spoa64::Graph& graph, spoa64::Graph::Node* startNode, float minAbundance, bool usePredecessors){

		spoa64::Graph::Node* currentNode = startNode;

		//cout << currentNode->id << endl;
		//cout << usePredecessors << endl;
		vector<u_int64_t> path;
		//path.push_back(currentNode->id);

		while(true){

			vector<spoa64::Graph::Node*> successors = getSolidSuccessors(currentNode, usePredecessors, minAbundance);
			if(successors.size() != 1) break;

			currentNode = successors[0];
			path.push_back(currentNode->id);


			//cout << currentNode->id << endl;
			//getchar();

		}

		return path;
	}

	*/
	
	vector<u_int64_t> performCorrection(spoa64::Graph& poaGraph, const vector<u_int64_t>& readMinimizers, const std::unique_ptr<spoa64::AlignmentEngine>& al, vector<u_int64_t>& correctedCoverages, unordered_set<u_int64_t>& minimizerSet){

		if(_illuminaFilename != ""){
			for (const auto& it : poaGraph.nodes()) {

				u_int64_t m1 = poaGraph.decoder(it->code);

				for (const auto& edge : it->outedges) {
					u_int64_t m2 = poaGraph.decoder(edge->head->code);
					
					KmerVec vec;
					vec._kmers = {m1, m2};
					vec = vec.normalize();

					if(_illuminaKminmers.find(vec) == _illuminaKminmers.end()){
					}
					else{
						edge->weight = _illuminaKminmers[vec];
					}

				}
			}
		}

		std::vector<u_int64_t> coverages;
		std::vector<u_int64_t> consensus = poaGraph.GenerateConsensusVec(coverages);

		/*
		//std::unique_ptr<spoa::AlignmentEngine> alignmentEngine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -4);
		
		spoa64::Graph graphAln{};

		vector<u_int64_t> weights(consensus.size(), 1);
		graphAln.AddAlignment(spoa64::Alignment(), consensus, consensus.size(), weights);

		spoa64::Alignment alignment = al->Align(readMinimizers, readMinimizers.size(), graphAln);


		u_int64_t startMatchPos = -1;
		u_int64_t endMatchPos = -1;

		for (const auto& it : alignment) {
			int64_t v1 = it.first;
			int64_t v2 = it.second;

			if(v1 == -1){ //insert in
				//cout << "insertion" << endl;
			}
			else if(v2 == -1){ //insert in
				//cout << "deletion" << endl;
			}
			else if(consensus[v1] == readMinimizers[v2]){

				if(startMatchPos == -1){
					startMatchPos = v1;
				}
				
				endMatchPos = v1;
			}

		}


		if(startMatchPos == -1 || startMatchPos == endMatchPos){
			vector<u_int64_t> correctReadMinimizers;
			return correctReadMinimizers;
		}
	
		//cout << startMatchPos << " " << endMatchPos << endl;

		vector<u_int64_t>::const_iterator first = consensus.begin() + startMatchPos;
		vector<u_int64_t>::const_iterator last = consensus.begin() + endMatchPos;
		vector<u_int64_t> correctReadMinimizers(first, last);

		vector<u_int64_t>::const_iterator first2 = coverages.begin() + startMatchPos;
		vector<u_int64_t>::const_iterator last2 = coverages.begin() + endMatchPos;
		vector<u_int64_t> correctCoverages(first2, last2);
		correctedCoverages = correctCoverages;

		//return correctReadMinimizers;
		*/

		spoa64::Alignment consensusAlignment = al->Align(consensus, consensus.size(), poaGraph);
		
		u_int64_t maxPathSize = readMinimizers.size() * 2;
		u_int64_t startNode;
		u_int64_t endNode;
		poaGraph.performCorrection(consensusAlignment, consensus, consensus.size(), startNode, endNode);

		//cout << (startNode ) << " " << (endNode ) << endl;
		vector<u_int64_t> correctReadMinimizers2 = computePath(poaGraph, poaGraph.nodes_[startNode].get(), poaGraph.nodes_[endNode].get(), maxPathSize, minimizerSet);

		//return correctReadMinimizers;
		//if(correctReadMinimizers2.size() == 0) return correctReadMinimizers;
		
		//cout << "Correction is valid" << endl;
		
		vector<u_int64_t> correctReadMinimizersTrimmed = trimCorrectedPath(readMinimizers, correctReadMinimizers2, al);
		
		return correctReadMinimizersTrimmed;
	}
	
	vector<u_int64_t> trimCorrectedPath(const vector<u_int64_t>& readMinimizers, const std::vector<u_int64_t>& correctedReadMinimizers, const std::unique_ptr<spoa64::AlignmentEngine>& al){

		spoa64::Graph graphAln{};

		vector<u_int64_t> weights(correctedReadMinimizers.size(), 1);
		graphAln.AddAlignment(spoa64::Alignment(), correctedReadMinimizers, correctedReadMinimizers.size(), weights);

		spoa64::Alignment alignment = al->Align(readMinimizers, readMinimizers.size(), graphAln);


		u_int64_t startMatchPos = -1;
		u_int64_t endMatchPos = -1;

		for (const auto& it : alignment) {
			int64_t v1 = it.first;
			int64_t v2 = it.second;

			if(v1 == -1){ //insert in
				//cout << "insertion" << endl;
			}
			else if(v2 == -1){ //insert in
				//cout << "deletion" << endl;
			}
			else if(correctedReadMinimizers[v1] == readMinimizers[v2]){

				//cout << "match" << endl;
				if(startMatchPos == -1){
					startMatchPos = v1;
				}
				
				endMatchPos = v1+1;
			}
			else{
				//cout << "missmatch" << endl;
			}

		}


		if(startMatchPos == -1 || startMatchPos == endMatchPos){
			vector<u_int64_t> correctReadMinimizers;
			return correctReadMinimizers;
		}
	
		//cout << startMatchPos << " " << endMatchPos << endl;

		vector<u_int64_t>::const_iterator first = correctedReadMinimizers.begin() + startMatchPos;
		vector<u_int64_t>::const_iterator last = correctedReadMinimizers.begin() + endMatchPos;
		vector<u_int64_t> correctedReadMinimizersTrimmed(first, last);

		return correctedReadMinimizersTrimmed;
	}

	
	struct PathSuccessor{
		u_int64_t _node;
		u_int32_t _distance;
	};


	struct DereplicatedEdge{
		spoa64::Graph::Edge* _edge;
		u_int64_t _weight;
	};
	
	vector<u_int64_t> computePath(const spoa64::Graph& graph, spoa64::Graph::Node* startNode, spoa64::Graph::Node* endNode, u_int64_t maxSize, unordered_set<u_int64_t>& readMinimizers){

		if(_print_debug) cout << "Compute path: " << startNode->id << "-" << graph.decoder(startNode->code) << " -> " << endNode->id << "-" << graph.decoder(endNode->code) << endl;
		
		vector<u_int64_t> path;
		//unordered_set<UnitigGraph::Node*> isVisited;
		//list<PathSuccessor> queue;
		//unordered_map<u_int64_t, u_int64_t> prev;
	
		//isVisited.insert(unitigSource);
		//isVisited.insert(_progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(unitigSource));

		//queue.push_back({fromNode, 0});
		//prev[fromNode] = -1;
		//bool found = false;
		spoa64::Graph::Node* currentNode = startNode;
		path.push_back(graph.decoder(currentNode->code));
		
		while (true) {
	

			if(_print_debug) cout << "\tVisit: " << currentNode->id << "-" << graph.decoder(currentNode->code) << endl;

			vector<DereplicatedEdge> successors = getSolidSuccessors(graph, currentNode, readMinimizers);

	


			for(const DereplicatedEdge& succ : successors){
				if(_print_debug) cout << "\t\t" << succ._edge->head->id << "-" << graph.decoder(succ._edge->head->code) << " " << succ._edge->weight << endl;

			}

			if(successors.size() == 0) break;
			//if(successors.size() == 0) return nullPath;
			if(successors.size() > 2 && successors[0]._weight <= 2) return path;

			currentNode = successors[0]._edge->head;
			path.push_back(graph.decoder(currentNode->code));

			if(currentNode == endNode) break;

			//if(path.size() > maxSize) return nullPath;
			
		}


		
		return path;
	}

	vector<DereplicatedEdge> getSolidSuccessors(const spoa64::Graph& graph, spoa64::Graph::Node* node, unordered_set<u_int64_t>& readMinimizers){

		vector<DereplicatedEdge> successors = getDereplicatedSuccessors(graph, node, readMinimizers);

		std::sort(successors.begin(), successors.end(), [](const DereplicatedEdge& a, const DereplicatedEdge& b){
			return a._weight > b._weight;
		});

		float maxWeight = 0;

		for(const DereplicatedEdge& edge : successors){
			if(edge._weight > maxWeight){
				maxWeight = edge._weight;
			}
		}

		float minWeight = maxWeight * 0.5;

		vector<DereplicatedEdge> solidEdges;

		for(const DereplicatedEdge& edge : successors){
			if(edge._weight < minWeight) continue;
				
			solidEdges.push_back(edge);
		}

		return solidEdges;
	}


	struct SuccessorCompletion{
		spoa64::Graph::Edge* _edge;
		u_int64_t _completion;
	};

	vector<DereplicatedEdge> getDereplicatedSuccessors(const spoa64::Graph& graph, spoa64::Graph::Node* node, unordered_set<u_int64_t>& readMinimizers){

		for(spoa64::Graph::Edge* edge : node->outedges){
			if(_print_debug) cout << "\t\tPossible successor: " << edge->head->id << " " << edge->weight << " " << computeEdgeCompletion(graph, edge, readMinimizers) << endl;
			//for(spoa64::Graph::Edge* edge2 : edge->head->outedges){
			//	cout << "\t\t\tLala: " << edge2->head->id << " " << edge2->weight << " " << computeEdgeCompletion(graph, edge2, readMinimizers) << endl;
			//}
		}

		unordered_map<u_int64_t, SuccessorCompletion> minimizer_to_edge;
		unordered_map<u_int64_t, u_int64_t> minimizer_to_weight;

		for(spoa64::Graph::Edge* edge : node->outedges){

			u_int64_t minimizer = edge->head->code;

			minimizer_to_weight[minimizer] += edge->weight;
			u_int64_t completion = computeEdgeCompletion(graph, edge, readMinimizers);

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

			u_int64_t minimizer = it.first;
			spoa64::Graph::Edge* edge = it.second._edge;

			u_int64_t derepWeight = minimizer_to_weight[minimizer];

			//edge->weight = minimizer_to_weight[it.first];
			derepEdges.push_back({edge, derepWeight});
			if(_print_debug) cout << "\t\tPossible successor derep: " << edge->head->id << " " << derepWeight << endl;
		}

		return derepEdges;

	}

	bool _print_debug;

	u_int64_t computeEdgeCompletion(const spoa64::Graph& graph, spoa64::Graph::Edge* edge, unordered_set<u_int64_t>& readMinimizers){

		u_int64_t nbVisitedReadNodes = 0;

		list<spoa64::Graph::Node*> queue;

		queue.push_back(edge->head);
		unordered_set<u_int64_t> isVisited;
	
		while (!queue.empty()) {
	
			spoa64::Graph::Node* currentNode = queue.front();
			queue.pop_front();

			if(isVisited.find(currentNode->id) != isVisited.end()) continue;
			isVisited.insert(currentNode->id);


			for(spoa64::Graph::Edge* successor : currentNode->outedges){

				if(readMinimizers.find(graph.decoder(successor->head->code)) != readMinimizers.end()){
					nbVisitedReadNodes += successor->weight;
				}

				queue.push_back(successor->head);
			}


		}

		return nbVisitedReadNodes;
	}

	/*
	struct PathSuccessor{
		u_int64_t _node;
		u_int32_t _distance;
	};

	vector<u_int64_t> computePath(const spoa64::Graph& graph, u_int64_t fromNode, u_int64_t toNode, float minAbundance){


		//unordered_set<UnitigGraph::Node*> isVisited;
		list<PathSuccessor> queue;
		unordered_map<u_int64_t, u_int64_t> prev;
	
		//isVisited.insert(unitigSource);
		//isVisited.insert(_progressiveAbundanceFilter->_unitigGraph->unitigIndex_toReverseDirection(unitigSource));

		queue.push_back({fromNode, 0});
		prev[fromNode] = -1;
		bool found = false;
	
		while (!queue.empty() && !found) {
	
			PathSuccessor pathSuccessor = queue.front();
			queue.pop_front();

			//cout << "\tVisit: " << pathSuccessor._node << endl;

			vector<spoa64::Graph::Node*> successors = getSolidSuccessors(graph.nodes_[pathSuccessor._node].get(), false, minAbundance);

			vector<spoa64::Graph::Node*> successors_noIsSucc = getSuccessors_noIsSuccessor(successors, false, minAbundance);

			for(spoa64::Graph::Node* node : successors_noIsSucc){
				//cout << "\t\t" << node->id << endl;

				prev[node->id] = pathSuccessor._node;

				if(node->id == toNode){
					found = true;
					break;
				}
				queue.push_back({node->id, 0});
			}


		}

		vector<u_int64_t> path;

		if(found){

			u_int64_t currentNode = toNode;
			while(currentNode != -1){
				path.push_back(currentNode);
				currentNode = prev[currentNode];
			}

			std::reverse(path.begin(), path.end());
		}
		//std::reverse(path.begin(), path.end());
		
		return path;
	}

	vector<spoa64::Graph::Node*> getSolidSuccessors(spoa64::Graph::Node* node, bool usePredecessors, float minAbundance){

		vector<spoa64::Graph::Edge*> successors = getDereplicatedSuccessors(node, usePredecessors, minAbundance);

		float maxWeight = 0;

		for(spoa64::Graph::Edge* edge : successors){
			if(edge->weight > maxWeight){
				maxWeight = edge->weight;
			}
		}

		float minWeight = maxWeight * 0.5;

		vector<spoa64::Graph::Node*> solidEdges;

		for(spoa64::Graph::Edge* edge : successors){
			if(minAbundance != -1 && edge->weight < minAbundance) continue;
			if(edge->weight < minWeight) continue;
				
			if(usePredecessors){
				solidEdges.push_back(edge->tail);
			}
			else{
				solidEdges.push_back(edge->head);
			}
		}

		return solidEdges;
	}


	vector<spoa64::Graph::Edge*> getDereplicatedSuccessors(spoa64::Graph::Node* node, bool usePredecessors, float minAbundance){

		unordered_map<u_int64_t, spoa64::Graph::Edge*> minimizer_to_nodeId;

		if(usePredecessors){
			for(spoa64::Graph::Edge* edge : node->inedges){

				u_int64_t minimizer = edge->tail->code;

				if(minimizer_to_nodeId.find(minimizer) != minimizer_to_nodeId.end()){
					if(edge->weight > minimizer_to_nodeId[minimizer]->weight){
						minimizer_to_nodeId[minimizer] = edge;
					}
				}
				else{
					minimizer_to_nodeId[minimizer] = edge;
				}
			}
		}
		else{
			for(spoa64::Graph::Edge* edge : node->outedges){

				u_int64_t minimizer = edge->head->code;

				if(minimizer_to_nodeId.find(minimizer) != minimizer_to_nodeId.end()){
					if(edge->weight > minimizer_to_nodeId[minimizer]->weight){
						minimizer_to_nodeId[minimizer] = edge;
					}
				}
				else{
					minimizer_to_nodeId[minimizer] = edge;
				}
			}
		}


		vector<spoa64::Graph::Edge*> derepEdges;

		for(const auto& it : minimizer_to_nodeId){
			derepEdges.push_back(it.second);
		}

		return derepEdges;

	}

	vector<spoa64::Graph::Node*> getSuccessors_noIsSuccessor(const vector<spoa64::Graph::Node*>& nodes, bool usePredecessors, float minAbundance){

		vector<spoa64::Graph::Node*> resultNodes;
		vector<spoa64::Graph::Node*> invalidNodes;

		for(size_t i=0; i<nodes.size(); i++){
			for(size_t j=0; j<nodes.size(); j++){
				
				if(isSuccessor(nodes[i], nodes[j], usePredecessors, minAbundance)){
					invalidNodes.push_back(nodes[j]);
				}
			}
		}

		for(spoa64::Graph::Node* node : nodes){
			if(std::find(invalidNodes.begin(), invalidNodes.end(), node) != invalidNodes.end()) continue;

			resultNodes.push_back(node);
		}

		return resultNodes;
	}

	bool isSuccessor(spoa64::Graph::Node* node1, spoa64::Graph::Node* node2, bool usePredecessors, float minAbundance){

		for(spoa64::Graph::Node* successor : getSolidSuccessors(node1, usePredecessors, minAbundance)){

			if(successor->code == node2->code) return true;
		}

		return false;

	}
	*/

	void writeRead(u_int32_t readIndex, const vector<u_int64_t>& minimizers){

		static u_int64_t maxHashValue = -1;
		static u_int64_t minimizerBound = 0.005 * maxHashValue;

		//#pragma omp critical(dataupdate)
		#pragma omp critical
		{
			

			//vector<u_int64_t> validMinimizers;
			//for(u_int64_t m : minimizers){
			//	if(m > minimizerBound) continue;
			//	validMinimizers.push_back(m);
			//}

			
			_readWriterQueue.push({readIndex, minimizers});

			while(!_readWriterQueue.empty()){


				const ReadWriter& readWriter = _readWriterQueue.top();

				//cout << readWriter._readIndex << " " << _nextReadIndexWriter << " " << _readWriterQueue.size() << endl;

				if(readWriter._readIndex == _nextReadIndexWriter){

					if(readWriter._minimizers.size() > 0){
						u_int32_t size = readWriter._minimizers.size();
						_file_readData.write((const char*)&size, sizeof(size));

						u_int8_t isCircular = CONTIG_LINEAR;
						_file_readData.write((const char*)&isCircular, sizeof(isCircular));

						_file_readData.write((const char*)&readWriter._minimizers[0], size*sizeof(u_int64_t));
					}

					_readWriterQueue.pop();
					_nextReadIndexWriter += 1;
				}
				else{
					break;
				}
			}
			
		}

	}

	void savePoaGraph(const spoa64::Graph& graph, const string& outputFilename, const vector<u_int64_t>& readNodes){

		ofstream outputFile(outputFilename);

		ofstream colorFile(outputFilename + ".color.csv");
		colorFile << "Name,Color" << endl;
		
		ofstream edgeFile(outputFilename + ".edge.csv");
		edgeFile << "Name,Edge" << endl;

		//std::cout << "H\tVN:Z:1.0" << std::endl;
		for (const auto& it : graph.nodes()) {

			string id = to_string(it->id) + "-" + to_string(graph.decoder(it->code));

			if(std::find(readNodes.begin(), readNodes.end(), it->id) != readNodes.end()){
				colorFile << (id) << ",green" << endl; 
			}
			else{
				colorFile << (id) << ",grey" << endl; 
			}

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

		/*
		for (std::uint32_t i = 0; i < graph.sequences().size(); ++i) {
			std::cout << "P\t" << headers[i] << "\t";

			std::vector<std::uint32_t> path;
			auto curr = graph.sequences()[i];
			while (true) {
			path.emplace_back(curr->id + 1);
			if (!(curr = curr->Successor(i))) {
				break;
			}
			}

			bool ir = !is_reversed.empty() && is_reversed[i];
			if (ir) {
			std::reverse(path.begin(), path.end());
			}
			for (std::uint32_t j = 0; j < path.size(); ++j) {
			if (j != 0) {
				std::cout << ",";
			}
			std::cout << path[j] << (ir ? "-" : "+");
			}
			std::cout << "\t*" << std::endl;
		}

		if (include_consensus) {
			std::cout << "P\tConsensus\t";
			for (std::uint32_t i = 0; i < graph.consensus().size(); ++i) {
			if (i != 0) {
				std::cout << ",";
			}
			std::cout << graph.consensus()[i]->id + 1 << "+";
			}
			std::cout << "\t*" << std::endl;
		}
		*/

		outputFile.close();
		colorFile.close();
		edgeFile.close();
	}


	/*
	void cleanGraph(spoa64::Graph& graph, spoa64::Graph& graphCleaned){


		graphCleaned.coder_ = graph.coder_;
		graphCleaned.decoder_ = graph.decoder_;
		
		unordered_map<u_int64_t, spoa64::Graph::Node*> id_to_node;
		unordered_map<u_int64_t, u_int64_t> id_to_coverage;


		for (const auto& it : graph.nodes()) {
			id_to_coverage[it->id] = it->Coverage();
		}

		for (const auto& it : graph.nodes()) {
			if(id_to_coverage[it->id] <= 1) continue;

			spoa64::Graph::Node* node = graphCleaned.AddNode(it->code);
			id_to_node[it->id] = node;
		}

		for (const auto& it : graph.nodes()) {

			if(id_to_node.find(it->id) == id_to_node.end()) continue;

			spoa64::Graph::Node* node = id_to_node[it->id];
			//string id = to_string(it->id) + "-" + to_string(graph.decoder(it->code));

			u_int64_t maxWeight = 0;

			for (const auto& edge : it->outedges) {
				if(id_to_node.find(edge->head->id) == id_to_node.end()) continue;

				if(edge->weight > maxWeight){
					maxWeight = edge->weight;
				}
			}

			double minWeight = maxWeight * 0.5;

			for (const auto& edge : it->outedges) {
				//if(edge->weight < minWeight) continue;
				if(id_to_node.find(edge->head->id) == id_to_node.end()) continue;

				spoa64::Graph::Node* toNode = id_to_node[edge->head->id];
      			graphCleaned.AddEdge(node, toNode, edge->weight);
			}
		}

		graphCleaned.TopologicalSort();
	}
	*/


};	


#endif 




