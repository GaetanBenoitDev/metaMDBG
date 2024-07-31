
//metaMDBG_testOv1

#ifndef MDBG_METAG_READCORRECTION
#define MDBG_METAG_READCORRECTION

#include "../Commons.hpp"
#include "../graph/GraphPOA.hpp"
#include "../utils/spoa64/include/spoa64/spoa.hpp"
#include "../graph/GfaParser.hpp"
#include "ReadMapper.hpp"

/*

./bin/metaMDBG asm ~/appa/run/correction/nanopore_mock_1M/asm/ ~/appa/data/nanopore/mock/zymo_hmw_r104_1M.fastq.gz -t 32

python3 ~/zeus/scripts/paper/computeReferenceCompleteness.py ~/appa/data/nanopore/mock/D6322.refseq/input.txt ~/appa/run/correction/nanopore_mock_1M/asm/contigs.fasta.gz ~/appa/run/correction/nanopore_mock_1M/asm/contigs.fasta.gz mdbg ~/appa/run/correction/nanopore_mock_1M/refComp/ 0.95 32

./bin/metaMDBG asm ~/appa/run/correction/nanopore_AD_circ1_asm/ ~/appa/data/nanopore/subreads/circ1.fastq -t 8
"determinisme: l'index kminmer -> readIndexes ne doit pas etre dans le meme ordre selon les run, il faut juste trier la liste des readIndex avant de collecter les mapepdReads"

"reprise: utilisé la mock nanopore pour comprendre prk l'assembly size est si differente entre k=11 k=13 firstK=3 firstk=4"

*/

class ReadCorrection : public Tool{
    
public:

	struct AlignmentResult{

		public:

		u_int64_t _readIndex;
		int64_t _nbMatches;
		int64_t _nbMissmatches;
		int64_t _nbInsertions;
		int64_t _nbDeletions;
		int64_t _alignLengthBps;

		int64_t score() const{
			int64_t nbErrors = _nbMissmatches + _nbInsertions + _nbDeletions;
			return _nbMatches - nbErrors;
		}
		
		//int64_t length() const{
		//	int64_t nbErrors = _nbMissmatches + _nbInsertions + _nbDeletions;
		//	return _nbMatches - nbErrors;
		//}

		float divergence() const{
			int64_t nbErrors = _nbMissmatches + _nbInsertions + _nbDeletions;
			return nbErrors / ((long double)(_nbMatches+nbErrors));
		}

	};



	
	struct MinimizerReadAlignment{
		u_int64_t _readIndex;
		vector<MinimizerType> _minimizers;
		vector<u_int8_t> _qualities;
		vector<u_int8_t> _readMinimizerDirections;
		AlignmentResult _alignmentResult;
	};


	typedef phmap::parallel_flat_hash_map<KmerVec, vector<u_int32_t>, phmap::priv::hash_default_hash<KmerVec>, phmap::priv::hash_default_eq<KmerVec>, std::allocator<std::pair<KmerVec, vector<u_int32_t>>>, 4, std::mutex> KminmerReadMap;
	typedef phmap::parallel_flat_hash_map<KmerVec, u_int32_t, phmap::priv::hash_default_hash<KmerVec>, phmap::priv::hash_default_eq<KmerVec>, std::allocator<std::pair<KmerVec, u_int32_t>>, 4, std::mutex> KminmerCountMap;

	typedef phmap::parallel_flat_hash_map<u_int64_t, vector<u_int64_t>, phmap::priv::hash_default_hash<u_int64_t>, phmap::priv::hash_default_eq<u_int64_t>, std::allocator<std::pair<u_int64_t, vector<u_int64_t>>>, 4, std::mutex> ReadToReadsMap;
	typedef phmap::parallel_flat_hash_map<u_int64_t, vector<MinimizerType>, phmap::priv::hash_default_hash<u_int64_t>, phmap::priv::hash_default_eq<u_int64_t>, std::allocator<std::pair<u_int64_t, vector<MinimizerType>>>, 4, std::mutex> ReadToMinimizersMap;

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
	u_int64_t _maxMappedReads;
	ReadToReadsMap _nanoporeIndex_to_illuminaIndexes;
	ReadToMinimizersMap _illuminaIndex_to_minimizers;

	u_int64_t _maxUsedReadOverlapForCorrection;
	bool _print_debug;

    struct ReadWriter{
        u_int64_t _readIndex;
        vector<MinimizerType> _minimizers;
        vector<u_int8_t> _qualities;
    };

    struct ReadWriter_Comparator {
        bool operator()(ReadWriter const& p1, ReadWriter const& p2){
            return p1._readIndex > p2._readIndex;
        }
    };

	priority_queue<ReadWriter, vector<ReadWriter> , ReadWriter_Comparator> _readWriterQueue;
	u_int64_t _nextReadIndexWriter;
	ofstream _file_readData;

	string _illuminaToNanoporeMappingFilename;
	//unordered_map<u_int64_t, u_int64_t> _minimizerCounts;
	//unordered_map<KmerVec, KminmerData> _kminmersData;
	//gzFile _file_minimizerPos;

	ofstream _debug_outputFile_nbMinimizerPairs;
	unordered_set<DbgEdge, hash_pair> _isReadPairWritten;



	/*
	class Edge;

	class Node{
		public:

		MinimizerType _minimizer;
		vector<Edge*> _successors;
		vector<Node*> _predecessors;
		u_int64_t _abundance;
		u_int64_t _quality;

		Node(const MinimizerType& minimizer){
			_minimizer = minimizer;
			_abundance = 0;
			_quality = 0;
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
		unordered_map<MinimizerType, Node*> _nodes;
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

		void addSequence(const vector<MinimizerType>& minimizers, const vector<u_int8_t>& qualities){


			if(minimizers.size() < 2) return;



			//for(size_t i=0; i<minimizers.size(); i++){
			//	cout << minimizers[i] << " " << (int)qualities[i] << endl;
			//}

			for(size_t i=0; i<minimizers.size()-1; i++){

				MinimizerType m1 = minimizers[i];
				MinimizerType m2 = minimizers[i+1];

				u_int8_t qual1 = qualities[i];
				u_int8_t qual2 = qualities[i+1];
				u_int8_t qual = min(qual1, qual2);

				Node* tail = addNode(m1);
				Node* head = addNode(m2);

				addEdge(tail, head, qual);

				//cout << (int)qual1 << " " << (int)qual2 << endl;
			}

			//getchar();

			for(size_t i=0; i<minimizers.size(); i++){
				Node* node = addNode(minimizers[i]);
				node->_abundance += 1;
				node->_quality += qualities[i];
			}
		}

		Node* addNode(const MinimizerType& minimizer){

			if(_nodes.find(minimizer) != _nodes.end()){
				return _nodes[minimizer];
			}
			//cout << "Add Node: " << unitigIndex << endl;

			Node* node = new Node(minimizer);
			_nodes[minimizer] = node;

			_needTopoSort = true;

			return node;
		}

		void addEdge(Node* tail, Node* head, u_int8_t weight){


			if(tail == nullptr) return;
			//if(tail == head) return;

			for(Edge* edge : tail->_successors){
				if(edge->_head == head){
					edge->_weight += weight;
					return;
				}
			}

			Edge* edge = new Edge(tail, head, weight);
			tail->_successors.push_back(edge);
			head->_predecessors.push_back(tail);

			_needTopoSort = true;

			//cout << "Add edge: " << fromNode->_nodeIndex << " -> " << toNode->_nodeIndex<< endl;
		}

		void save(const string& outputFilename, const vector<MinimizerType>& readMinimizers){

			cout << "save dbg graph" << endl;

			ofstream outputFile(outputFilename);

			ofstream colorFile(outputFilename + ".color.csv");
			colorFile << "Name,Color" << endl;
			
			ofstream edgeFile(outputFilename + ".edge.csv");
			edgeFile << "Name,Edge" << endl;

			//std::cout << "H\tVN:Z:1.0" << std::endl;
			for (const auto& it : _nodes) {

				MinimizerType minimizer = it.first;
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


		vector<u_int64_t> computeConsensus(const vector<MinimizerType>& readMinimizers, bool& isCycle){

			isCycle = false;

			unordered_set<MinimizerType> readKminmers;
			for(const MinimizerType& vec : readMinimizers){
				readKminmers.insert(vec);
			}

			PathWeight path = computePathWeight(_nodes[readMinimizers[0]], readKminmers);

			vector<u_int64_t> minimizerPath;

			//cout << endl << "Most supported path: " << endl;
			//int trimStart = 0;
			//int trimEnd = 0;

			for(size_t i =0; i<path._path.size(); i++){

				Node* node = path._path[i];

				minimizerPath.push_back(node->_minimizer);
				//cout << "\t" << node->toString() << " " << node->_isReversed << endl; 


				//if(readKminmers.find(node->_minimizer) != readKminmers.end()){
				//	trimEnd = minimizerPath.size();
				//}

			}


			//vector<u_int64_t>::const_iterator first = minimizerPath.begin() + trimStart;
			//vector<u_int64_t>::const_iterator last = minimizerPath.begin() + trimEnd;
			//vector<u_int64_t> minimizerPathTrimmed(first, last);

			isCycle = path._isCycle;

			return minimizerPath;
		}

		struct PathWeight{
			vector<Node*> _path;
			u_int64_t _weight;
			bool _isCycle;
		};


		PathWeight computePathWeight(Node* startNode, unordered_set<MinimizerType>& readKminmers){

			PathWeight pathWeight;


			//list<spoa64::Graph::Node*> queue;

			//queue.push_back(edge->head);
			Node* currentNode = startNode;
			unordered_set<Node*> isVisited;
			
			pathWeight._weight = 0;
			pathWeight._path.push_back(currentNode);
			pathWeight._isCycle = false;
		
			//cout << "\t\tStart node: " << startNode->toString() << endl;
			
			while (true) {
		

				//cout << "\t\tVisit node: " << currentNode->toString() << endl;
				if(isVisited.find(currentNode) != isVisited.end()){
					//
					
					//cout << "\t\tCycle detected allo ??: " << currentNode->toString() << endl;

					if(pathWeight._path.size() < readKminmers.size() +2){
						//cout << "\t\tCycle detected: " << currentNode->toString() << endl;
						pathWeight._isCycle = true;
					}

					break;
				}

				isVisited.insert(currentNode);

				u_int64_t maxWeight = 0;
				Edge* maxSuccessor = nullptr;

				for(Edge* successor : currentNode->_successors){
					
					if(successor->_weight <= 1) continue;

					if(successor->_weight > maxWeight){
						maxWeight = successor->_weight;
						maxSuccessor = successor;
					}
				}

				float repeatCutoff = maxWeight*0.9;
				int nbSolidSuccessors = 0;

				for(Edge* successor : currentNode->_successors){
					if(successor->_weight > repeatCutoff){
						nbSolidSuccessors += 1;
					}
				}

				if(nbSolidSuccessors > 1 && pathWeight._path.size() < readKminmers.size() +2){
					//cout << maxWeight << endl;
					//pathWeight._isCycle = true;
				}

				if(maxSuccessor == nullptr) break;

				//cout << "\t\tBest successor:: " << maxSuccessor->_head->toString() << endl;

				currentNode = maxSuccessor->_head;

				pathWeight._path.push_back(currentNode);

				//cout << readKminmers.size() << " " << (readKminmers.find(currentNode->_minimizer) != readKminmers.end()) << endl;

				//for(const auto& it : readKminmers){
				//	cout << it._kmers[0] << " " << it._kmers[1] << endl;
				//}
				if(readKminmers.find(currentNode->_minimizer) != readKminmers.end()){
					//cout << "Add weight: " << currentNode->_abundance << endl;
					pathWeight._weight += maxSuccessor->_weight;
				}

				//getchar();
			}

			//if(pathWeight._isCycle){
			//	cout << "loul" << endl;
			//	getchar();
			//}

			return pathWeight;
		}


	};
	*/

	class Edge;

	class Node{
		public:

		u_int32_t _nodeIndex;
		MinimizerType _minimizer;
		vector<Edge*> _successors;
		vector<Node*> _predecessors;
		u_int64_t _abundance;
		u_int64_t _quality;
		u_int64_t _maxQuality;

		vector<AlignmentResult> _alignmentResults;

		Node(u_int32_t nodeIndex, const MinimizerType& minimizer){
			_nodeIndex = nodeIndex;
			_minimizer = minimizer;
			_abundance = 0;
			_quality = 0;
			_maxQuality = 0;
		}

		string toString(){
			return to_string(_nodeIndex) + "-" + to_string(_minimizer);
		}

		void addQuality(u_int64_t quality){
			_quality += quality;
			_maxQuality = max(_maxQuality, quality);
		}

		vector<AlignmentResult> getBestAlignments(int n){

			vector<AlignmentResult> bestAlignments;

			std::sort(_alignmentResults.begin(), _alignmentResults.end(), [](const AlignmentResult& a, const AlignmentResult& b){
				return a.score() > b.score();
			});

			for(int i=0; i<n && i<_alignmentResults.size(); i++){
				bestAlignments.push_back(_alignmentResults[i]);
			}

			return bestAlignments;
		}
	};


	class Edge{

		public:

		Node* _tail;
		Node* _head;
		u_int32_t _weight;
		u_int32_t _support;

		Edge(Node* tail, Node* head, u_int32_t weight){
			_tail = tail;
			_head = head;
			_weight = weight;
			_support = 1;
		}
		
	};

	class Graph{
		public:

		bool _print;
		unordered_map<u_int32_t, Node*> _nodes;

		Graph(const vector<MinimizerType>& minimizers, const vector<u_int8_t>& qualities){

			_print = false;

			for(size_t i=0; i<minimizers.size(); i++){
				Node* node = addNode(i, minimizers[i]);
				node->_abundance += 1;
				node->addQuality(qualities[i]);
			}

			for(size_t i=0; i<minimizers.size()-1; i++){

				//MinimizerType m1 = minimizers[i];
				//MinimizerType m2 = minimizers[i+1];

				u_int8_t qual1 = qualities[i];
				u_int8_t qual2 = qualities[i+1];
				u_int8_t qual = min(qual1, qual2);

				//Node* tail = addNode(i, m1, qual1);
				//Node* head = addNode(i+1, m2, qual2);

				//cout << _nodes[i]->toString() << " " << _nodes[i+1]->toString() << " " << (int) qual1 << " " << (int)  qual2 << endl;
				//cout << (int)  edge->_weight << " " << (int) weight << endl;
				addEdge(_nodes[i], _nodes[i+1], qual);
			}
		}

		~Graph(){
			for(auto& it : _nodes){

				for(Edge* edge : it.second->_successors){
					delete edge;
				}
				delete it.second;
			}
			
		}

		void addAlignment(const spoa64::Alignment& alignment, const vector<MinimizerType>& referenceMinimizers, const vector<MinimizerType>& readMinimizers, const vector<MinimizerType>& readQualities,const AlignmentResult alignmentResult){

			if(_print) {
				cout << endl << "---" << endl;
				for(u_int64_t m : readMinimizers){
					cout << "\t" << m << endl;
				}
			}
			//for(Edge* edge : _nodes[7]->_successors){
			//	cout << "\t" << edge->_head->toString() << " " << (int) edge->_weight << endl;
			//}
			Node* prevNode = nullptr;
			Node* currentNode = nullptr;
			//if(readMinimizers.size() < 2) return;

			//spoa64::Alignment alignment = al->Align(s2, s2.size(), graph);

			for (const auto& it : alignment) {
				int64_t v1 = it.first;
				int64_t v2 = it.second;
				
				if(v1 == -1){ //insert in
					if(_print) cout << "Insertion " << readMinimizers[v2] << endl;

					currentNode = addNode(prevNode, readMinimizers[v2]);
					currentNode->_abundance += 1;
					currentNode->addQuality(readQualities[v2]);
					currentNode->_alignmentResults.push_back(alignmentResult);

					if(prevNode != nullptr){
						addEdge(prevNode, currentNode, readQualities[v2]);
					}

				}
				else if(v2 == -1){ //insert in
					if(_print) cout << "Deletion " << referenceMinimizers[v1] << endl;
				}
				else if(referenceMinimizers[v1] == readMinimizers[v2]){

					if(_print) cout << "Match " << referenceMinimizers[v1] << " " << readMinimizers[v2] << endl;
					//cout << _nodes.size() << endl;
					currentNode = _nodes[v1];
					currentNode->_abundance += 1;
					currentNode->addQuality(readQualities[v2]);
					currentNode->_alignmentResults.push_back(alignmentResult);

					//cout << (currentNode != nullptr) << " " << (prevNode!= nullptr) << endl;

					if(prevNode != nullptr){
						addEdge(prevNode, currentNode, readQualities[v2]);
					}

					
				}
				else{

					if(_print) cout << "Missmatch " << referenceMinimizers[v1] << " " << readMinimizers[v2] << endl;


					currentNode = addNode(prevNode, readMinimizers[v2]);
					currentNode->_abundance += 1;
					currentNode->addQuality(readQualities[v2]);
					currentNode->_alignmentResults.push_back(alignmentResult);

					if(prevNode != nullptr){
						addEdge(prevNode, currentNode, readQualities[v2]);
					}


				}

				//cout << "Current node: " << currentNode->toString() << endl;
				prevNode = currentNode;
			}

		}

			/*
			//for(size_t i=0; i<minimizers.size(); i++){
			//	cout << minimizers[i] << " " << (int)qualities[i] << endl;
			//}

			for(size_t i=0; i<minimizers.size()-1; i++){

				MinimizerType m1 = minimizers[i];
				MinimizerType m2 = minimizers[i+1];

				u_int8_t qual1 = qualities[i];
				u_int8_t qual2 = qualities[i+1];
				u_int8_t qual = min(qual1, qual2);

				Node* tail = addNode(m1);
				Node* head = addNode(m2);

				addEdge(tail, head, qual);

				//cout << (int)qual1 << " " << (int)qual2 << endl;
			}

			//getchar();

			for(size_t i=0; i<minimizers.size(); i++){
				Node* node = addNode(minimizers[i]);
				node->_abundance += 1;
				node->_quality += qualities[i];
			}
		}

		*/


		Node* addNode(u_int32_t nodeIndex, const MinimizerType& minimizer){

			if(_nodes.find(nodeIndex) != _nodes.end()){
				//_nodes[nodeIndex]->_abundance += 1;
				//_nodes[nodeIndex]->_quality += quality;
				return _nodes[nodeIndex];
			}

			Node* node = new Node(nodeIndex, minimizer);
			_nodes[nodeIndex] = node;

			if(_print) cout << "\tAdd node: " << nodeIndex << " " << minimizer << endl;

			return node;

		}

		Node* addNode(Node* prevNode, const MinimizerType& minimizer){

			if(prevNode != nullptr){
				for(Edge* nn : prevNode->_successors){
					if(nn->_head->_minimizer == minimizer){
						//nn->_head->_abundance += 1;
						//nn->_head->_quality += quality;
						return nn->_head;
					}
				}
			}

			u_int32_t nodeIndex = _nodes.size();

			Node * node = addNode(nodeIndex, minimizer);
			/*
			for(Node* nn :)
			//if(_nodes.find(minimizer) != _nodes.end()){
			//	return _nodes[minimizer];
			//}
			//cout << "Add Node: " << unitigIndex << endl;

			Node* node = new Node(minimizer);
			_nodes[minimizer] = node;


			return node;
			*/
			return node;
		}

		void addEdge(Node* tail, Node* head, u_int8_t weight){

			
			if(tail == nullptr) return;
			//if(tail == head) return;

			for(Edge* edge : tail->_successors){
				if(edge->_head == head){

					//if(edge->_head->toString() == "8-71093523068023552" && edge->_tail->toString() == "7-6626939888570") 
					//cout << (int)  edge->_weight << " " << (int) weight << endl;
					edge->_weight += weight;
					edge->_support += 1;
					return;
				}
			}

			Edge* edge = new Edge(tail, head, weight);
			tail->_successors.push_back(edge);
			head->_predecessors.push_back(tail);
			
			if(_print) cout << "\tAdd edge: " << tail->toString() << " -> " << head->toString() << endl;

			//cout << "Add edge: " << fromNode->_nodeIndex << " -> " << toNode->_nodeIndex<< endl;
		}

		struct EdgeWeight{
			Edge* _edge;
			u_int64_t _weight;
		};

		void save(const string& outputFilename, const vector<MinimizerType>& readMinimizers){

			cout << "save dbg graph" << endl;

			ofstream outputFile(outputFilename);

			ofstream colorFile(outputFilename + ".color.csv");
			colorFile << "Name,Color" << endl;
			
			ofstream edgeFile(outputFilename + ".edge.csv");
			edgeFile << "Name,Edge" << endl;

			ofstream alignResultFile(outputFilename + ".alignResult_mean.csv");
			alignResultFile << "Name,Edge" << endl;

			ofstream alignResultMaxFile(outputFilename + ".alignResult_max.csv");
			alignResultMaxFile << "Name,Edge" << endl;

			ofstream alignResultLengthFile(outputFilename + ".alignResult_length.csv");
			alignResultLengthFile << "Name,Edge" << endl;

			//std::cout << "H\tVN:Z:1.0" << std::endl;
			for (const auto& it : _nodes) {

				MinimizerType minimizer = it.first;
				Node* node = it.second;

				string id = node->toString();

				
				int64_t n = 0;
				int64_t sumMatches = 0;
				int64_t sumMissmatches = 0;
				int64_t sumInsertions = 0;
				int64_t sumDeletions = 0;
				int64_t sumLength = 0;


				int64_t maxMatches = 0;
				int64_t maxMissmatches = 0;
				int64_t maxInsertions = 0;
				int64_t maxDeletions = 0;
				int64_t maxLength = 0;

				//cout << "lala " << node->_alignmentResults.size() << endl;
				
				for(const AlignmentResult al : node->getBestAlignments(5)){

					maxMatches = max(maxMatches, al._nbMatches);
					maxMissmatches = max(maxMissmatches, al._nbMissmatches);
					maxInsertions = max(maxInsertions, al._nbInsertions);
					maxDeletions = max(maxDeletions, al._nbDeletions);
					maxLength = max(maxLength, al._alignLengthBps);

					//cout << "loulou " << al._nbMatches << endl;
					sumMatches += al._nbMatches;
					sumMissmatches += al._nbMissmatches;
					sumInsertions += al._nbInsertions;
					sumDeletions += al._nbDeletions;
					sumLength += al._alignLengthBps;
					n += 1;
				}

				if(n > 0){
					int64_t meanMatches = sumMatches / n;
					int64_t meanMissmatches = sumMissmatches / n;
					int64_t meanInsertions = sumInsertions / n;
					int64_t meanDeletions = sumDeletions / n;
					int64_t meanLength = sumLength / n;

					string alignmentResultStr = to_string(meanMatches) + " | " + to_string(meanMissmatches) + " | " +to_string(meanInsertions) + " | " +to_string(meanDeletions);
					alignResultFile << id << "," << alignmentResultStr << endl; 

					string alignmentResultMaxStr = to_string(maxMatches) + " | " + to_string(maxMissmatches) + " | " +to_string(maxInsertions) + " | " +to_string(maxDeletions);
					alignResultMaxFile << id << "," << alignmentResultMaxStr << endl; 

					string alignmentResultLengthStr = to_string(meanLength) + " | " + to_string(maxLength);
					alignResultLengthFile << id << "," << alignmentResultLengthStr << endl; 
				}
				

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

				vector<EdgeWeight> successors;

				for (Edge* edge : node->_successors) {
					successors.push_back({edge, edge->_weight});
				}

				std::sort(successors.begin(), successors.end(), [](const EdgeWeight& a, const EdgeWeight& b){
					return a._weight > b._weight;
				});

				string edgeStr = "";


				for (size_t i=0; i<successors.size() && i < 5; i++) {
					
					const EdgeWeight& successor = successors[i];
					//cout << id << " -> " << edge->_head->toString() << " " << edge->_weight << endl;

					//if(edge->_weight < 2) continue;

					outputFile << "L\t" << id << "\t" << "+\t" << successor._edge->_head->toString() << "\t" << "+\t" << "1M" << endl; //\t" << "ew:f:" << jt->weight << endl;


					//if(successor._edge->_weight > 1){
						edgeStr += "[" + successor._edge->_head->toString() + " - " + to_string(successor._edge->_weight) + "] ";
					//}
				}

				edgeFile << id << "," << edgeStr << endl; 



			}



			outputFile.close();
			colorFile.close();
			edgeFile.close();
			alignResultFile.close();
			alignResultMaxFile.close();
			alignResultLengthFile.close();

			cout << "done" << endl;
		}

		/*
		vector<u_int64_t> computeConsensus(const vector<MinimizerType>& readMinimizers, bool& isCycle){

			isCycle = false;

			unordered_set<MinimizerType> readKminmers;
			for(const MinimizerType& vec : readMinimizers){
				readKminmers.insert(vec);
			}



			PathWeight path = computePathWeight(_nodes[readMinimizers[0]], readKminmers);

			vector<u_int64_t> minimizerPath;

			//cout << endl << "Most supported path: " << endl;
			//int trimStart = 0;
			//int trimEnd = 0;

			for(size_t i =0; i<path._path.size(); i++){

				Node* node = path._path[i];

				minimizerPath.push_back(node->_minimizer);
				//cout << "\t" << node->toString() << " " << node->_isReversed << endl; 


				//if(readKminmers.find(node->_minimizer) != readKminmers.end()){
				//	trimEnd = minimizerPath.size();
				//}

			}


			//vector<u_int64_t>::const_iterator first = minimizerPath.begin() + trimStart;
			//vector<u_int64_t>::const_iterator last = minimizerPath.begin() + trimEnd;
			//vector<u_int64_t> minimizerPathTrimmed(first, last);

			isCycle = path._isCycle;

			return minimizerPath;
		}

		struct PathWeight{
			vector<Node*> _path;
			u_int64_t _weight;
			bool _isCycle;
		};


		PathWeight computePathWeight(Node* startNode, unordered_set<MinimizerType>& readKminmers){

			PathWeight pathWeight;


			//list<spoa64::Graph::Node*> queue;

			//queue.push_back(edge->head);
			Node* currentNode = startNode;
			unordered_set<Node*> isVisited;
			
			pathWeight._weight = 0;
			pathWeight._path.push_back(currentNode);
			pathWeight._isCycle = false;
		
			//cout << "\t\tStart node: " << startNode->toString() << endl;
			
			while (true) {
		

				//cout << "\t\tVisit node: " << currentNode->toString() << endl;
				if(isVisited.find(currentNode) != isVisited.end()){
					//
					
					//cout << "\t\tCycle detected allo ??: " << currentNode->toString() << endl;

					if(pathWeight._path.size() < readKminmers.size() +2){
						//cout << "\t\tCycle detected: " << currentNode->toString() << endl;
						pathWeight._isCycle = true;
					}

					break;
				}

				isVisited.insert(currentNode);

				u_int64_t maxWeight = 0;
				Edge* maxSuccessor = nullptr;

				for(Edge* successor : currentNode->_successors){
					
					if(successor->_weight <= 1) continue;

					if(successor->_weight > maxWeight){
						maxWeight = successor->_weight;
						maxSuccessor = successor;
					}
				}

				float repeatCutoff = maxWeight*0.9;
				int nbSolidSuccessors = 0;

				for(Edge* successor : currentNode->_successors){
					if(successor->_weight > repeatCutoff){
						nbSolidSuccessors += 1;
					}
				}

				if(nbSolidSuccessors > 1 && pathWeight._path.size() < readKminmers.size() +2){
					//cout << maxWeight << endl;
					//pathWeight._isCycle = true;
				}

				if(maxSuccessor == nullptr) break;

				//cout << "\t\tBest successor:: " << maxSuccessor->_head->toString() << endl;

				currentNode = maxSuccessor->_head;

				pathWeight._path.push_back(currentNode);

				//cout << readKminmers.size() << " " << (readKminmers.find(currentNode->_minimizer) != readKminmers.end()) << endl;

				//for(const auto& it : readKminmers){
				//	cout << it._kmers[0] << " " << it._kmers[1] << endl;
				//}
				if(readKminmers.find(currentNode->_minimizer) != readKminmers.end()){
					//cout << "Add weight: " << currentNode->_abundance << endl;
					pathWeight._weight += maxSuccessor->_weight;
				}

				//getchar();
			}

			//if(pathWeight._isCycle){
			//	cout << "loul" << endl;
			//	getchar();
			//}

			return pathWeight;
		}

		*/
	};


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

		//_illuminaFilename = "/pasteur/appa/homes/gbenoit/appa/data/nanopore/AD/input_paired.txt";
		//_illuminaFilename = "/pasteur/appa/homes/gbenoit/appa/data/nanopore/subreads/circ1_illumina.txt";
		
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
		_maxMappedReads = 100;
		_minNbMinimizersReference = 10;
		_minNbMinimizersQuery = 10;
		_maxUsedReadOverlapForCorrection = 100;
		_illuminaToNanoporeMappingFilename = _inputDir + "/illuminaMapping.paf"; //_illuminaFilename + ".paf.tmp";
		_nbRepeats = 0;
		_nbRepeatsDBG = 0;
		_nbRepeatsMissed = 0;
		_nbCorrectedRepeats = 0;
		_correctionCheckSum = 0;
		_eval_nbBases = 0;
		_eval_nbMatches = 0;

		_fields = new vector<string>();
	}
	u_int64_t _nbRepeats;
	u_int64_t _nbRepeatsDBG;
	u_int64_t _nbRepeatsMissed;
	u_int64_t _nbCorrectedRepeats;

	u_int64_t _eval_nbBases;
	u_int64_t _eval_nbMatches;
	vector<string>* _fields;


    void execute (){


		_nextReadIndexWriter = 0;

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
		cout << "pb: trimming: si il y a un minimizer répété, parfois le match final est trop loin car l'alignement accept mass insertion, plutot que prendre le minimizer repete le plus proche comme match" << endl;
		
		//cout << "Mapping reads" << endl;
		//mapReads();

		//cout << "Loading read overlaps" << endl;


		_print_debug = false;
		_eval_correction = false;

		
		if(_eval_correction){
			cout << "Loading read to contig mapping" << endl;
			loadReadToContigMapping();

			cout << "Loading contigs" << endl;
			loadContigs();



			//cout << "Loading read to read mapping" << endl;
			//loadReadToReadMapping();


		}

		//processReads();
		
		cout << "Loading minimizer reads" << endl;
		loadMinimizerReads(_inputDir + "/read_data_init.txt", false);


		cout << "Indexing mReads" << endl;
		indexReads();



		cout << "Correcting mReads" << endl;
		_file_readData = ofstream(_inputDir + "/read_data_corrected.txt");


		correctReads();

		_file_readData.close();
		


		cout << "Correction check sum: " << _correctionCheckSum << endl;
		

		closeLogFile();
	}



	/*
	void processReads(){

		ReadMapper readMapper(_inputDir + "/read_data_init.txt", _inputDir + "/read_data_init.txt", _minNbMinimizersQuery, _nbCores);
		readMapper.execute();

	}
	*/

	class LoadReadFunctor {

		public:

		ReadCorrection& _parent;

		LoadReadFunctor(ReadCorrection& parent) : _parent(parent){
		}

		LoadReadFunctor(const LoadReadFunctor& copy) : _parent(copy._parent){
		}

		~LoadReadFunctor(){
		}


		void operator () (const KminmerList& kminmerList) {

			u_int32_t readIndex = kminmerList._readIndex;
			if(readIndex % 100000 == 0) cout << readIndex << endl;

			MinimizerRead read = {kminmerList._readIndex, kminmerList._readMinimizers, kminmerList._minimizerPos, kminmerList._readQualities, kminmerList._readMinimizerDirections};
			
			_parent._mReads.push_back(read);
			
		}
	};

	class LoadReadChunkFunctor {

		public:

		ReadCorrection& _parent;

		LoadReadChunkFunctor(ReadCorrection& parent) : _parent(parent){
		}

		LoadReadChunkFunctor(const LoadReadChunkFunctor& copy) : _parent(copy._parent){
		}

		~LoadReadChunkFunctor(){
		}


		void operator () () const {
			cout << "Loaded read chunk: " << _parent._mReads.size() << endl;

			_parent._totalReadProcessed += _parent._mReads.size();
			_parent._mReads.clear();
		}
	};


	u_int64_t _totalReadProcessed;



















	KminmerReadMap _kminmer_to_readIndex;
	


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
				//_parent._mReadsHiFi.push_back({kminmerList._readIndex, kminmerList._readMinimizers, kminmerList._readQualities, kminmerList._readMinimizerDirections});
			}
			else{

				MinimizerRead read = {kminmerList._readIndex, kminmerList._readMinimizers, kminmerList._minimizerPos, kminmerList._readQualities, kminmerList._readMinimizerDirections};
				
				_parent._mReads.push_back(read);
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
	
	/*
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
	*/
	//"reprise: verfiier qu'on ne perd pas des reads, noeud de part erroné, mauvais aprcours de graphe, trimming etc faire le isReversed check avec le local alignment"

	//"reprise: evaluer le coverage d'un read, le restituer non corrigé si coverage trop bas, assess seulement le coverage apres correction? eng ros trim des partie basses couverture avant evaluation"
	//"min coverage: au moins 3 pour le first et last solid position ?"

	u_int64_t _correctionCheckSum;

	bool _eval_correction;
	void correctReads(){





		int nbCores = _nbCores;
		if(_print_debug || _eval_correction){  
			nbCores = 1;
		}

		//nbCores = 1;

		KminmerParserParallel parser2(_inputDir + "/read_data_init.txt", _minimizerSize, _kminmerSize, false, true, nbCores);
		//parser2.parse(FilterKminmerFunctor2(*this));
		parser2.parseSequences(ReadCorrectionFunctor(*this));

	}

	class ReadCorrectionFunctor {

		public:

		ReadCorrection& _parent;
		//std::unique_ptr<spoa64::AlignmentEngine> _alignmentEngine;// = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -4);
		//std::unique_ptr<spoa64::AlignmentEngine> _alignmentEngineTrimming;
		std::unique_ptr<spoa64::AlignmentEngine> _alignmentEngine_ov;
		
		ReadCorrectionFunctor(ReadCorrection& parent) : _parent(parent){
		}

		ReadCorrectionFunctor(const ReadCorrectionFunctor& copy) : _parent(copy._parent){
			//_alignmentEngine = spoa64::AlignmentEngine::Create(spoa64::AlignmentType::kNW, 3, -5, -4);
			//_alignmentEngineTrimming = spoa64::AlignmentEngine::Create(spoa64::AlignmentType::kSW, 3, -2, -2);
			_alignmentEngine_ov = spoa64::AlignmentEngine::Create(spoa64::AlignmentType::kOV, 3, -1, -1, -1);
		}

		~ReadCorrectionFunctor(){
		}


		void operator () (const KminmerList& kminmerList) {
			_parent.correctRead(kminmerList, _alignmentEngine_ov);
		}
	};
	
	struct MinimizerReadOverlap{
		vector<MinimizerType> _minimizers;
		float _error;

	};

	struct CorrectedRead{
		u_int64_t _readIndex;
		vector<MinimizerType> _minimizers;
		vector<u_int8_t> _qualities;
		vector<u_int8_t> _readMinimizerDirections;
		vector<MinimizerType> _originalMinimizers;
		vector<u_int8_t> _originalQualities;
	};


	struct MinimizerAlignment{
		vector<MinimizerType> _minimizerSequence;
		vector<MinimizerType> _minimizerQualities;
		spoa64::Alignment _alignment;
		AlignmentResult _alignmentResult;
	};

	void correctRead(const KminmerList& kminmerList, const std::unique_ptr<spoa64::AlignmentEngine>& alOverlap){

		u_int64_t readIndex = kminmerList._readIndex;
		if(readIndex % 10000 == 0) cout << kminmerList._readIndex << endl;

		if(_print_debug){
			cout << endl << endl;
			cout << "-------------------------------------------------------" << endl;
			cout << "Correcting read: " << "ReadIndex=" << kminmerList._readIndex << " Readsize=" << kminmerList._readMinimizers.size() << endl;
			for(MinimizerType m : kminmerList._readMinimizers){
				cout << m << endl;
			}
			cout << endl << endl;

		}

		#pragma omp atomic
		_nbCorrectedRepeats += 1;

		CorrectedRead correctedRead = {kminmerList._readIndex, kminmerList._readMinimizers, kminmerList._readQualities, kminmerList._readMinimizerDirections, kminmerList._readMinimizers, kminmerList._readQualities};

		for(size_t i=0; i<1; i++){

			if(correctedRead._minimizers.size() < _minNbMinimizersReference){ //|| kminmerList._readMinimizers.size() > 250
				break;
			}

			correctedRead = correctRead_iteration(correctedRead, alOverlap);

			//if(_print_debug){
			//	cout << "Iteration " << i << ": Readsize=" << correctedRead._minimizers.size() << endl;
			//	getchar();
			//}
		}

		writeRead(correctedRead._readIndex, correctedRead._minimizers, correctedRead._qualities);
	}

	CorrectedRead correctRead_iteration(const CorrectedRead& read, const std::unique_ptr<spoa64::AlignmentEngine>& alOverlap){



		//if(read._readIndex != 718) return read;
		

		//Overlappin read: 2070
		//Overlappin read: 8958
		//Overlappin read: 6060
		//Overlappin read: 8235

		//if(readIndex != 8235) return;
		
		//if(readIndex % 1000 == 0) cout << readIndex << endl;
		//if(readIndex < 100) return;
		

		//if(readIndex < 93136) return;
		//if(readIndex < 221572) return;

		//#pragma omp critical
		//{
		//	cout << readIndex << endl;
		//}



		//if(readIndex % 1000 == 0) cout << "tension read correction disabled" << endl;





		//#pragma omp critical
		//{
		//	cout << _nbRepeats << " / " << _nbRepeatsDBG << " / " << _nbRepeatsMissed << " / " << _nbCorrectedRepeats << endl;
		//}
		

		CorrectedRead readWindow = {read._readIndex, read._minimizers, read._qualities, read._readMinimizerDirections, read._originalMinimizers, read._originalQualities};
		const CorrectedRead& correctedRead = correctWindow(0, readWindow, alOverlap);
		
		return correctedRead;
		
		/*
		vector<MinimizerType> correctedReadMinimizersTrim = correctWindow(0, readIndex, kminmerList._readMinimizers, kminmerList._readMinimizerDirections, kminmerList._readQualities, _print_debug, al, alTrimming, false);
		
		cout << endl;
		cout << endl;
		cout << "Read index: " << readIndex << endl;

		cout << "lala" << endl;
		for(size_t i=0; i<kminmerList._readMinimizers.size(); i++){
			cout << kminmerList._readMinimizers[i] << " " << (int) kminmerList._readQualities[i] << endl;
		}
		cout << "---" << endl;
		for(u_int64_t m : correctedReadMinimizers){
			cout << m << endl;
		}
		cout << "---" << endl;
		for(u_int64_t m : correctedReadMinimizersTrim){
			cout << m << endl;
		}

		cout << kminmerList._readMinimizers.size() << " " << correctedReadMinimizers.size() << " " << correctedReadMinimizersTrim.size() << endl;
		
		size_t diff = abs(((int)correctedReadMinimizers.size()) - ((int)correctedReadMinimizersTrim.size()));
		if(diff > 5) getchar();
		//getchar();
		*/

		/*
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

		//writeRead(readIndex, correctedReadMinimizers);
		//getchar();
		
	}


	size_t _timeQueryingReads;
	size_t _timeMappingReads;
	size_t _timeCorrecting;

	struct ReadIndexMatch{
		u_int64_t _readIndex;
		int64_t _matchScore;
	};

	CorrectedRead correctWindow(size_t windowIndex, const CorrectedRead& read, const std::unique_ptr<spoa64::AlignmentEngine>& alOverlap){



		unordered_map<MinimizerType, u_int32_t> minimizerPosition;
		unordered_map<MinimizerType, u_int8_t> readMinimizer_to_direction;
		unordered_set<MinimizerType> minimizerSet;

		for(size_t i=0; i<read._minimizers.size(); i++){
			MinimizerType m = read._minimizers[i];

			readMinimizer_to_direction[m] = read._readMinimizerDirections[i];
			minimizerSet.insert(m);

			if(minimizerPosition.find(m) != minimizerPosition.end()){ //repeated minimizer
				minimizerPosition[m] = -1;
			}
			else{
				minimizerPosition[m] = i;
			}
		}
		
		auto time = high_resolution_clock::now();

		//Kminmers !!
		vector<ReadKminmerComplete> kminmersInfos;
		vector<u_int32_t> minimizersPos(read._minimizers.size(), 0);
		//vector<u_int8_t> minimizerQualities(readMinimizers.size(), 0);
		MDBG::getKminmers_complete(_kminmerSize, read._minimizers, minimizersPos, kminmersInfos, read._readIndex, read._qualities);
		

		unordered_map<u_int32_t, u_int32_t> readIndex_to_matchCount;
		vector<KmerVec> readKminmers;

		for(const ReadKminmerComplete& kminmerInfo : kminmersInfos){

			const KmerVec& vec = kminmerInfo._vec;
			readKminmers.push_back(vec);

			if(_kminmer_to_readIndex.find(vec) == _kminmer_to_readIndex.end()) continue;

			for(u_int32_t rIndex : _kminmer_to_readIndex[vec]){
				if(read._readIndex == rIndex) continue; //Currently corrected read
				readIndex_to_matchCount[rIndex] += 1;
			}

		}

		vector<ReadIndexMatch> matchingReadIndexes;

		for(const auto& it : readIndex_to_matchCount){

			u_int64_t readIndex = it.first;
			int64_t nbMatches = it.second;

			if(nbMatches < 2) continue;

			//int64_t alignScoreFast = getAlignmentScoreFast(_mReads[readIndex], minimizerSet, minimizerPosition, readMinimizer_to_direction);
		

			matchingReadIndexes.push_back({readIndex, nbMatches});
		}

		//return read;

		std::sort(matchingReadIndexes.begin(), matchingReadIndexes.end(), [](const ReadIndexMatch & a, const ReadIndexMatch & b){
			if(a._matchScore == b._matchScore){
				return a._readIndex > b._readIndex;
			}
			return a._matchScore > b._matchScore;
		});



		unordered_map<MinimizerType, vector<u_int64_t>> minimizer_to_bestReadIndexes;

		for(size_t i=0; i<matchingReadIndexes.size(); i++){

			u_int64_t readIndex = matchingReadIndexes[i]._readIndex;

			for(MinimizerType minimizer : _mReads[readIndex]._minimizers){
				if(minimizerSet.find(minimizer) == minimizerSet.end()) continue;
				if(minimizer_to_bestReadIndexes[minimizer].size() >= _maxMappedReads) continue;

				minimizer_to_bestReadIndexes[minimizer].push_back(readIndex);
			}
		}


		unordered_set<u_int64_t> selectedReads;

		for(const auto& it : minimizer_to_bestReadIndexes){
			u_int64_t minimizer = it.first;
			const vector<u_int64_t>& readIndexes = it.second;

			for(u_int64_t readIndex : readIndexes){
				selectedReads.insert(readIndex);
			}
		}





		vector<MinimizerRead> minimizerSequences;


		for(u_int64_t readIndex : selectedReads){

			const MinimizerRead& matchingReadMinimizers = _mReads[readIndex];

			minimizerSequences.push_back(matchingReadMinimizers);

		}
		
		
		//std::sort(minimizerSequences.begin(), minimizerSequences.end(), [this, &minimizerSet](const MinimizerRead & a, const MinimizerRead & b){
		//	if(a._minimizers.size() == b._minimizers.size()){
		//		return a._readIndex < b._readIndex;
		//	}
		//	return a._minimizers.size() > b._minimizers.size();
		//	//return a.size() > b.size();
		//});

		#pragma omp atomic
		_timeQueryingReads += duration_cast<microseconds>(high_resolution_clock::now() - time).count();
		//auto time = high_resolution_clock::now();
		//getchar();
		

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
		if(_print_debug) cout << "Nb reads used for correction: " << minimizerSequences.size() << endl;


		time = high_resolution_clock::now();

		vector<MinimizerReadAlignment> mappedReads;
		getMappedMinimizerReads(minimizerSequences, read, minimizerSet, minimizerPosition, readMinimizer_to_direction, mappedReads, alOverlap);

		#pragma omp atomic
		_timeMappingReads += duration_cast<microseconds>(high_resolution_clock::now() - time).count();

		if(_print_debug) cout << "Nb mapped mReads: " << mappedReads.size() << endl;
		
		/*
		#pragma omp critical
		{

			if(isArtifact(read)){
				cout << "Nb mapped mReads: " << mappedReads.size() << endl;
				getchar();
			}

		}
		*/


		//vector<MinimizerType> correctedReadMinimizers;

		//if(trimMappedReads){
		//	correctedReadMinimizers = performPoaCorrection(windowIndex, readIndex, readMinimizers, readMinimizerDirections, readQualities, al, alTrimming, mappedReads, minimizerSet);
		//}
		//else{


		time = high_resolution_clock::now();

		const CorrectedRead& correctedRead = performPoaCorrection2(windowIndex, read, alOverlap, mappedReads, minimizerSet);
		
		#pragma omp atomic
		_timeCorrecting += duration_cast<microseconds>(high_resolution_clock::now() - time).count();

		
		long double totalTime = _timeQueryingReads + _timeMappingReads + _timeCorrecting;
		//cout << "Time querying: " << _timeQueryingReads/totalTime << endl;
		//cout << "Time mapping: " << _timeMappingReads/totalTime << endl;
		//cout << "Time correcting: " << _timeCorrecting/totalTime << endl;
		//}

		/*
		if(isRepeatedMinimizer){
			isRepeat = true;


			
		}
		else{
			
			Graph* dbgGraph = new Graph();
			dbgGraph->addSequence(readMinimizers, readQualities);


			for(size_t i=0; i<mappedReads.size(); i++){

				const MinimizerRead& minimizerRead = mappedReads[i];


				dbgGraph->addSequence(minimizerRead._minimizers, minimizerRead._qualities);

			}

			vector<u_int64_t> consensusDBG = dbgGraph->computeConsensus(readMinimizers, isCycle);
			correctedReadMinimizersDBG = trimCorrectedPath(readMinimizers, consensusDBG, alTrimming);
			

			//if(_print_debug) dbgGraph->save(_inputDir + "/poaGraph_DBG.gfa", readMinimizers);

			delete dbgGraph;

			if(isCycle){
				correctedReadMinimizers = performPoaCorrection(windowIndex, readIndex, readMinimizers, readMinimizerDirections, readQualities, al, alTrimming, mappedReads, minimizerSet);
			
			}
			else{
				correctedReadMinimizers = correctedReadMinimizersDBG;
			}

			//correctedReadMinimizersLala = performPoaCorrection(windowIndex, readIndex, readMinimizers, readMinimizerDirections, readQualities, al, alTrimming, mappedReads, minimizerSet);

			unordered_set<u_int64_t> mSet;
			for(u_int64_t m :correctedReadMinimizersLala){

				if(mSet.find(m) == mSet.end()){
					mSet.insert(m);
				}
				else{
					repeatedMinimizers.insert(m);
					isRepeat = true;
				}
			}



			if(isCycle){
				#pragma omp atomic
				_nbRepeatsDBG += 1;

			}

			if(isRepeat && !isCycle){
				#pragma omp atomic
				_nbRepeatsMissed += 1;

			}

		}


	
		if(isRepeat){
			#pragma omp atomic
			_nbRepeats += 1;

		}
		*/

		//correctedReadMinimizers = trimByAbundance(correctedReadMinimizers, mappedReads);
		
		//vector<MinimizerType> correctedReadMinimizers = performPoaCorrection(windowIndex, readIndex, readMinimizers, readMinimizerDirections, readQualities, al, alTrimming, mappedReads, minimizerSet);

		if(_eval_correction){
			CorrectedRead fakeCorrectedRead = evaluateCorrection(read._readIndex, correctedRead, alOverlap);
			//return fakeCorrectedRead;
			//if(!isValid){
			//	return {correctedRead._readIndex, correctedRead._originalMinimizers, correctedRead._originalQualities, correctedRead._readMinimizerDirections, correctedRead._originalMinimizers, correctedRead._originalQualities};
	
			//}
		}

		
		if(_print_debug){

			

			cout << endl << "\tOriginal sequence:" << endl;
			for(size_t i=0; i<read._minimizers.size(); i++){
				cout << "\t" << read._minimizers[i] << " " << ((int)read._qualities[i]) << endl;
			}


			cout << endl << "\tCorrected sequence:" << endl;
			for(size_t i=0; i<correctedRead._minimizers.size(); i++){
				cout << "\t" << correctedRead._minimizers[i] << " " << ((int)correctedRead._qualities[i]) << endl;;
			}


			cout << read._minimizers.size() << " " << correctedRead._minimizers.size() << endl;


			getchar();
			
			
			




		}
		
		

		return correctedRead;
	}

	//vector<MinimizerType> trimSequenceFast(const MinimizerRead& read, unordered_set<MinimizerType>& minimizerSet){

	//	int64_t start = -1;
	//	int64_t end = -1;
	//}
	/*
	bool isArtifact(const CorrectedRead& read){

		unordered_map<u_int64_t, u_int64_t> minimizer_to_count;

		for(u_int64_t m : read._minimizers){
			minimizer_to_count[m] += 1;
		}

		long double nbRepeatedMinimizers = 0;

		for(const auto& it : minimizer_to_count){
			if(it.second > 1){
				nbRepeatedMinimizers += 1;
			}
		}

		double fractionRepeatedMinimizer = nbRepeatedMinimizers / minimizer_to_count.size();
		if(fractionRepeatedMinimizer > 0.5){
			cout << endl;
			cout << "Read index: " << read._readIndex << endl;
			cout << "Read size: " << read._minimizers.size() << endl;
			for(size_t i=0; i<read._minimizers.size(); i++){

				u_int64_t m = read._minimizers[i];
				u_int64_t qual = read._qualities[i];

				if(minimizer_to_count[m] > 1){
					cout << "\t" << m << " " << qual << " *" << endl;
				}
				else{
					cout << "\t" << m << " " << qual << endl;
				}
			}

			return true;
		}

		return false;
	}
	*/

	vector<u_int64_t> vec32_to_vec64(const vector<u_int32_t>& vec32){

		vector<u_int64_t> vec64(vec32.size(), 0);

		for(size_t i=0; i<vec32.size(); i++){
			vec64[i] = vec32[i];
		}

		return vec64;
	}

	void getMappedMinimizerReads(const vector<MinimizerRead>& minimizerSequences, const CorrectedRead& read, unordered_set<MinimizerType>& minimizerSet, unordered_map<MinimizerType, u_int32_t>& minimizerPosition, unordered_map<MinimizerType, u_int8_t>& readMinimizer_to_direction, vector<MinimizerReadAlignment>& mappedReads, const std::unique_ptr<spoa64::AlignmentEngine>& al){


		//cout << "Nb collected reads: " <<  minimizerSequences.size() << endl;
		mappedReads.clear();
		vector<MinimizerReadAlignment> mappedReadsTmp;
		spoa64::Graph graph{};
		

		vector<MinimizerType> weights(read._minimizers.size(), 1);
		graph.AddAlignment(spoa64::Alignment(), vec32_to_vec64(read._minimizers), read._minimizers.size(), vec32_to_vec64(weights));
		
		unordered_map<MinimizerType, vector<AlignmentResult>> minimizer_to_alignmentResults;
		

		for(size_t i=0; i<minimizerSequences.size(); i++){

			const MinimizerRead& mappedRead = minimizerSequences[i];

			vector<MinimizerType> minimizerSequence = mappedRead._minimizers;
			vector<u_int8_t> qualities = mappedRead._qualities;
			vector<u_int8_t> directions = mappedRead._readMinimizerDirections;

			//trimSimple(mOverlap._minimizers, mOverlap._qualities, mOverlap._readMinimizerDirections, minimizerSet, minimizerSequence, qualities, directions);




			//bool isReversed = isSequenceReversed(readMinimizers, minimizerSequence, al);
			bool isReversed = isSequenceReversed2(mappedRead, readMinimizer_to_direction);

			//cout << isReversed << " " << isReversed2 << endl;
			if(isReversed){
				std::reverse(minimizerSequence.begin(), minimizerSequence.end());
				std::reverse(qualities.begin(), qualities.end());
			}



			/*
			cout << endl << "---" << endl;
			for(u_int64_t m : readMinimizers){
				cout << m << endl;
			}
			cout << "---" << endl;
			for(u_int64_t m : minimizerSequence){
				cout << m << endl;
			}
			cout << "---" << endl;
			for(u_int64_t m : minimizerSequenceT){
				cout << m << endl;
			}
			getchar();
			*/

			//int nbMatches = 0;
			//int nbInsertions = 0;
			//int nbDeletions = 0;
			//int nbMissmatch = 0;
			//int nbErrors = 0;
			//int lastMatchPos = 0;
			AlignmentResult alignmentResult = getAlignmentScore(mappedRead._readIndex, graph, read._minimizers, minimizerSequence, mappedRead._minimizersPos, al);
			//double alignmentScore = nbMatches / (double) read._minimizers.size();// / (double) minimizerSequence.size();
			
			//if(_print_debug){
			//	cout << "\tAl score: " << alignmentScore << endl; //<< " " << getAlignmentScore(minimizerSequence, readMinimizers, al) << endl;
			//}

			int nbErrors = alignmentResult._nbMissmatches;// + alignmentResult._nbInsertions + alignmentResult._nbDeletions;
			
			//if(alignmentResult._nbMatches <= 5) continue;
			//if(alignmentResult.divergence() > 0.3) continue;

			//int64_t minNbMatches = 20;//min((int64_t)20, (int64_t)(read._minimizers.size() / 2));
			//if(alignmentResult._nbMatches < minNbMatches) continue;
			if(alignmentResult._nbMatches-nbErrors < 5) continue;

			//if(nbMatches-nbErrors < 10) continue;
			//if(alignmentScore < 0.2) continue;
			//if(alignmentScore < 5) continue;

			//if(isReversed){
			//	std::reverse(qualities.begin(), qualities.end());
			//}
	
			/*
			if(_print_debug){

				if(mappedRead._readIndex == 20527){
					
					//cout << endl << "\tAlign sequence: " << i << " " << minimizerSequence.size() << " " << graph.nodes_.size() << endl; //<< " " << graphPOA->_graph->_nodes.size() << endl;
					
					for(MinimizerType m : read._minimizers){
						cout << "\t" << m << endl; 
					}
					cout << endl;
					for(MinimizerType m : minimizerSequence){

						
						if(std::find(read._minimizers.begin(), read._minimizers.end(), m) != read._minimizers.end()){
							cout << "\t" << m << " *" << endl; 
						}
						else{
							cout << "\t" << m << endl; 
						}

					}
				}
			
			}
			*/

			for(MinimizerType minimizer : mappedRead._minimizers){
				if(minimizerSet.find(minimizer) == minimizerSet.end()) continue;

				minimizer_to_alignmentResults[minimizer].push_back(alignmentResult);
			}
			
			vector<MinimizerType> minimizerSequenceT;
			vector<u_int8_t> qualitiesT;
			//"reprise: si c'est bien le mode trimNew2 la, reessayer sans le +=5 minimizers, et le +completion dans computepath"
			trimCorrectedPathFast(graph, read._minimizers, minimizerSequence, qualities, al, minimizerSequenceT, qualitiesT);
			//if(trimMappedReads) trimSimple(mOverlap._minimizers, mOverlap._qualities, mOverlap._readMinimizerDirections, minimizerSet, minimizerSequence, qualities, directions);

			//cout << endl;
			//for(u_int64_t m : mOverlap._minimizers){
			//	cout << "\t" << m << endl;
			//}
			//cout << "---" << endl;
			//for(u_int64_t m : minimizerSequence){
			//	cout << "\t" << m << endl;
			//}

			
			mappedReadsTmp.push_back({mappedRead._readIndex, minimizerSequenceT, qualitiesT, {}, alignmentResult});
			
		}

		unordered_set<u_int64_t> selectedReads;

		for(const auto& it : minimizer_to_alignmentResults){
			u_int64_t minimizer = it.first;
			const vector<AlignmentResult>& alignmentResults = it.second;

			const vector<AlignmentResult>& bestAlignmentResults  = getBestAlignments(alignmentResults, 20);

			for(const AlignmentResult& alignmentResult : bestAlignmentResults){
				//cout << alignmentResult._readIndex << " " << alignmentResult.score() << " " << getAlignmentScoreFast(_mReads[alignmentResult._readIndex], minimizerSet, minimizerPosition, readMinimizer_to_direction) << endl;
				//cout << read._minimizers.size() << "\t" << alignmentResult._nbMatches << "\t" << alignmentResult._nbMissmatches << "\t" << alignmentResult._nbInsertions << "\t" << alignmentResult._nbDeletions << endl;//<< " " << getAlignmentScoreFast(_mReads[alignmentResult._readIndex], minimizerSet, minimizerPosition, readMinimizer_to_direction) << endl;
				selectedReads.insert(alignmentResult._readIndex);
			}

			//cout << "check" << endl;
			//getchar();
		}

		for(const MinimizerReadAlignment& mappedRead : mappedReadsTmp){

			if(selectedReads.find(mappedRead._readIndex) == selectedReads.end()) continue;

			mappedReads.push_back(mappedRead);
		}

	}

	
	int64_t getAlignmentScoreFast(const MinimizerRead& mappedRead, unordered_set<MinimizerType>& minimizerSet, unordered_map<MinimizerType, u_int32_t>& minimizerPosition, unordered_map<MinimizerType, u_int8_t>& readMinimizer_to_direction){
		

		bool isReversed = isSequenceReversed2(mappedRead, readMinimizer_to_direction);
		
		vector<MinimizerType> minimizers = mappedRead._minimizers;
		if(isReversed) std::reverse(minimizers.begin(), minimizers.end());

		int64_t startPositionOriginal = -1;
		int64_t endPositionOriginal = -1;
		int64_t startPositionMapped = -1;
		int64_t endPositionMapped = -1;

		int64_t startPositionOriginalBest = -1;
		int64_t endPositionOriginalBest = -1;
		int64_t startPositionMappedBest = -1;
		int64_t endPositionMappedBest = -1;

		int64_t bestScore = 0;
		int64_t currentScore = 0;
		int64_t missmatchChain = 0;
		bool alignStarted = false;

		for(size_t i=0; i<minimizers.size(); i++){

			MinimizerType m = minimizers[i];

			if(minimizerPosition.find(m) == minimizerPosition.end()){
				
				missmatchChain += 1;

				if(missmatchChain > 5){
					if(currentScore > bestScore){
						startPositionOriginalBest = startPositionOriginal;
						endPositionOriginalBest = endPositionOriginal;
						startPositionMappedBest = startPositionMapped;
						endPositionMappedBest = endPositionMapped;
						bestScore = currentScore;
					}
					alignStarted = false;
				}
			}
			else{
				if(!alignStarted){
					startPositionMapped = i;
					startPositionOriginal = minimizerPosition[m];
					alignStarted = true;
					missmatchChain = 0;
					currentScore = 0;
				}
				currentScore -= missmatchChain;
				missmatchChain = 0;
				currentScore += 1;
				endPositionOriginal = minimizerPosition[m];
				endPositionMapped = i;
			}
		}

		if(currentScore > bestScore){
			bestScore = currentScore;
			startPositionOriginalBest = startPositionOriginal;
			endPositionOriginalBest = endPositionOriginal;
			startPositionMappedBest = startPositionMapped;
			endPositionMappedBest = endPositionMapped;
		}

		int64_t hangLeft = min(startPositionOriginalBest, startPositionMappedBest);
		int64_t hangRight = min(minimizerSet.size()-endPositionOriginalBest, minimizers.size()-endPositionMappedBest);

		bestScore -= hangLeft;
		bestScore -= hangRight;

		/*
		if(mappedRead._readIndex == 20527){
			//cout << isReversed << endl;
			cout << startPositionOriginalBest << " " << endPositionOriginalBest << endl;
			cout << startPositionMappedBest << " " << endPositionMappedBest << endl;
			cout << minimizerSet.size() << " " << minimizers.size() << endl;
			cout << hangLeft << " " << hangRight << endl;
		}
		*/
	
		return bestScore;
	}
	

	vector<AlignmentResult> getBestAlignments(vector<AlignmentResult> alignmentResults, int n){

		vector<AlignmentResult> bestAlignments;

		std::sort(alignmentResults.begin(), alignmentResults.end(), [](const AlignmentResult& a, const AlignmentResult& b){
			return a.score() > b.score();
		});

		for(int i=0; bestAlignments.size() < n && i < alignmentResults.size(); i++){

			bool isHere = false;
			for(const AlignmentResult& al : bestAlignments){
				if(al._readIndex == alignmentResults[i]._readIndex){ //The same alignment can be present several time if a minimizer is repeated in a read
					isHere = true;
					break;
				}
			}

			if(isHere) continue;

			bestAlignments.push_back(alignmentResults[i]);
		}

		return bestAlignments;
	}

	/*
	vector<MinimizerType> performPoaCorrection(size_t windowIndex, u_int64_t readIndex, const vector<MinimizerType>& readMinimizers, const vector<u_int8_t>& readMinimizerDirections, const vector<u_int8_t>& readQualities, const std::unique_ptr<spoa64::AlignmentEngine>& al, const std::unique_ptr<spoa64::AlignmentEngine>& alTrimming, const vector<MinimizerRead>& mappedReads, unordered_set<MinimizerType>& minimizerSet){

		vector<ReadKminmerComplete> kminmersInfos;
		vector<u_int32_t> minimizersPos(readMinimizers.size(), 0);
		MDBG::getKminmers_complete(_kminmerSize, readMinimizers, minimizersPos, kminmersInfos, readIndex, readQualities);
		

		unordered_map<KmerVec, u_int32_t> kminmer_to_abundance;

		for(const ReadKminmerComplete& kminmerInfo : kminmersInfos){
			const KmerVec& vec = kminmerInfo._vec;
			kminmer_to_abundance[vec] = 1;
		}




		//cout << endl << endl;
		//for(u_int64_t m : readMinimizers){
		//	cout << m << endl;
		//}

		
		//Graph* dbgGraph = new Graph(readMinimizers, readQualities);

		
		spoa64::Graph graph{};

		
		MinimizerRead minimizerReadCorr = {0, readMinimizers, readQualities}; //applyIlluminaCorrection({0, readMinimizers, readQualities}, illuminaMinimizers);

		vector<MinimizerType> weights(minimizerReadCorr._qualities.size(), 1);
		for(size_t i=0; i<minimizerReadCorr._qualities.size(); i++){
			weights[i] = minimizerReadCorr._qualities[i];
		}

		//vector<MinimizerType> readNodes;
		graph.AddAlignment(spoa64::Alignment(), minimizerReadCorr._minimizers, minimizerReadCorr._minimizers.size(), weights);




		
		vector<MinimizerAlignment> alignments;

		for(size_t i=0; i<mappedReads.size(); i++){

			const MinimizerRead& minimizerRead = mappedReads[i];
			//MinimizerRead minimizerReadIlluminaCorrected = applyIlluminaCorrection(mappedReads[i], illuminaMinimizers);


			vector<u_int64_t> seq = minimizerRead._minimizers;
			//seq.insert(seq.begin()+4, 10);
			//seq.insert(seq.begin()+4, 11);

			spoa64::Alignment alignment = al->Align(
							seq,
							seq.size(),
							graph);

			vector<MinimizerType> qualities(minimizerRead._qualities.size(), 0);
			for(size_t i=0; i<minimizerRead._qualities.size(); i++){
				qualities[i] = minimizerRead._qualities[i];
			}


			alignments.push_back({seq, qualities, alignment, 0});
			
		}

		std::shuffle (alignments.begin(), alignments.end(), std::default_random_engine(42));


		//std::sort(alignments.begin(), alignments.end(), [](const MinimizerAlignment & a, const MinimizerAlignment & b){
		//	return a._alignmentScore > b._alignmentScore;
			//return a._lastMatchPos < b._lastMatchPos;
		//});


		u_int64_t nbUsedReads = 0;
		double nbMappedMinimizers;

		for(size_t i=0; i<alignments.size(); i++){

			MinimizerAlignment& alignment = alignments[i];

			vector<ReadKminmerComplete> kminmersInfos;
			vector<u_int32_t> minimizersPos(alignment._minimizerSequence.size(), 0);
			vector<u_int8_t> minimizersQualities(alignment._minimizerSequence.size(), 0);
			MDBG::getKminmers_complete(_kminmerSize, alignment._minimizerSequence, minimizersPos, kminmersInfos, readIndex, minimizersQualities);
			

			bool isCovered = false;

			for(const ReadKminmerComplete& kminmerInfo : kminmersInfos){

				const KmerVec& vec = kminmerInfo._vec;
				if(kminmer_to_abundance.find(vec) == kminmer_to_abundance.end()) continue;

				kminmer_to_abundance[vec] += 1;

				if(kminmer_to_abundance[vec] > 200) isCovered = true;
				
			}


			if(isCovered) break;





			if(_print_debug){
				cout << endl << "\tAdd alignment: " << i << " " << alignment._minimizerSequence.size() << " " << alignment._alignmentScore << " " << graph.nodes_.size() << endl; //<< " " << graphPOA->_graph->_nodes.size() << endl;
				for(MinimizerType m : alignment._minimizerSequence){
					cout << "\t" << m << endl; 
				}
			
			}

		
			//dbgGraph->addAlignment(alignment._alignment, readMinimizers, alignment._minimizerSequence, alignment._minimizerQualities);

			//printAlignmentScore(alignment._alignment, readMinimizers, alignment._minimizerSequence, al);

			graph.AddAlignment(alignment._alignment, alignment._minimizerSequence, alignment._minimizerSequence.size(), alignment._minimizerQualities);

			nbMappedMinimizers += alignment._alignmentScore;
			nbUsedReads += 1;

			//double coverage = nbMappedMinimizers / readMinimizers.size();
			//if(coverage > 50) break;

			//savePoaGraph(graph, _inputDir + "/poaGraph.gfa", readNodes);
			//getchar();


		}
		
		//cout << endl;
		//for(const auto& it: kminmer_to_abundance){
		//	cout << it.second << endl;
				
		//}

		//cout << "Size: " << readMinimizers.size() << endl;
		
		//cout << "Nb reads used in graph: " << nbUsedReads << endl;

		//if(readMinimizers.size() > 80) 	getchar();
		
 		std::vector<MinimizerType> correctedCoverages;
		//std::reverse(correctedSequence.begin(), correctedSequence.end());
		vector<MinimizerType> correctedReadMinimizers = performCorrection(graph, readMinimizers, al, alTrimming, correctedCoverages, minimizerSet);
		
		//if(_print_debug)
		//dbgGraph->save(_inputDir + "/poaGraph.gfa", readMinimizers);

		//vector<MinimizerType> correctedReadMinimizers = computePath2(dbgGraph);
		//delete dbgGraph;


		//if(correctedReadMinimizers.size() == 0) return readMinimizers;

		//vector<MinimizerType> correctReadMinimizersTrimmed = trimCorrectedPath(readMinimizers, correctedReadMinimizers, alTrimming);
		
		//
		//getchar();

		savePoaGraph(graph, _inputDir + "/poaGraph.gfa", readMinimizers);

		//vector<MinimizerType> correctedReadMinimizers;
		return correctedReadMinimizers;
	}
	*/
		
	Graph* _debugDbgGraph;
	/*
	CorrectedRead performPoaCorrection2(size_t windowIndex, const CorrectedRead& read, const std::unique_ptr<spoa64::AlignmentEngine>& alOverlap, const vector<MinimizerReadAlignment>& mappedReads, unordered_set<MinimizerType>& minimizerSet){

		
		Graph* dbgGraph = new Graph(read._minimizers, read._qualities);

		
		spoa64::Graph graph{};

		
		//MinimizerRead minimizerReadCorr = {0, readMinimizers, readQualities}; //applyIlluminaCorrection({0, readMinimizers, readQualities}, illuminaMinimizers);

		vector<MinimizerType> weights(read._qualities.size(), 1);
		for(size_t i=0; i<read._qualities.size(); i++){
			weights[i] = read._qualities[i];
		}

		//vector<MinimizerType> readNodes;
		graph.AddAlignment(spoa64::Alignment(), vec32_to_vec64(read._minimizers), read._minimizers.size(), vec32_to_vec64(weights));




		
		vector<MinimizerAlignment> alignments;

		for(size_t i=0; i<mappedReads.size(); i++){

			const MinimizerReadAlignment& minimizerRead = mappedReads[i];
			//MinimizerRead minimizerReadIlluminaCorrected = applyIlluminaCorrection(mappedReads[i], illuminaMinimizers);


			vector<u_int64_t> seq = minimizerRead._minimizers;
			//seq.insert(seq.begin()+4, 10);
			//seq.insert(seq.begin()+4, 11);

			spoa64::Alignment alignment = alOverlap->Align(seq, seq.size(), graph);

			vector<MinimizerType> qualities(minimizerRead._qualities.size(), 0);
			for(size_t i=0; i<minimizerRead._qualities.size(); i++){
				qualities[i] = minimizerRead._qualities[i];
			}


			alignments.push_back({seq, qualities, alignment, minimizerRead._alignmentResult});
			
		}

		//std::shuffle (alignments.begin(), alignments.end(), std::default_random_engine(42));


		//std::sort(alignments.begin(), alignments.end(), [](const MinimizerAlignment & a, const MinimizerAlignment & b){
		//	return a._alignmentScore > b._alignmentScore;
			//return a._lastMatchPos < b._lastMatchPos;
		//});


		u_int64_t nbUsedReads = 0;

		for(size_t i=0; i<alignments.size(); i++){

			MinimizerAlignment& alignment = alignments[i];




			//if(_print_debug){
			//	cout << endl << "\tAdd alignment: " << i << " " << alignment._minimizerSequence.size() << " " << alignment._alignmentScore << " " << graph.nodes_.size() << endl; //<< " " << graphPOA->_graph->_nodes.size() << endl;
			//	for(MinimizerType m : alignment._minimizerSequence){
			//		cout << "\t" << m << endl; 
			//	}
			
			//}

		
			dbgGraph->addAlignment(alignment._alignment, read._minimizers, alignment._minimizerSequence, alignment._minimizerQualities, alignment._alignmentResult);

			//printAlignmentScore(alignment._alignment, readMinimizers, alignment._minimizerSequence, al);

			//graph.AddAlignment(alignment._alignment, alignment._minimizerSequence, alignment._minimizerSequence.size(), alignment._minimizerQualities);

			nbUsedReads += 1;

			//double coverage = nbMappedMinimizers / readMinimizers.size();
			//if(coverage > 50) break;

			//savePoaGraph(graph, _inputDir + "/poaGraph.gfa", readNodes);



		}
		//7-773482283752690 6-10963043653607913
		//cout << endl;
		//for(const auto& it: kminmer_to_abundance){
		//	cout << it.second << endl;
				
		//}

		//cout << "Size: " << readMinimizers.size() << endl;
		
		//cout << "Nb reads used in graph: " << nbUsedReads << endl;

		//if(readMinimizers.size() > 80) 	getchar();
		
 		//std::vector<MinimizerType> correctedCoverages;
		//std::reverse(correctedSequence.begin(), correctedSequence.end());
		//vector<MinimizerType> correctedReadMinimizers = performCorrection(graph, readMinimizers, al, alTrimming, correctedCoverages, minimizerSet);
		
		//if(_print_debug)
		//dbgGraph->save(_inputDir + "/poaGraph_new.gfa", readMinimizers);

	
		const CorrectedRead& correctedRead = computePath2(read, dbgGraph, minimizerSet);
		//const CorrectedRead& correctedRead = computePath3(read, dbgGraph, minimizerSet);

		if(_eval_correction){
			_debugDbgGraph = dbgGraph;
		}
		else{
			delete dbgGraph;
		}


		//if(correctedReadMinimizers.size() == 0) return readMinimizers;

		const CorrectedRead& correctedReadTrimmed = trimCorrectedPath(read, correctedRead, alOverlap);
		
		//
		//getchar();


		//vector<MinimizerType> correctedReadMinimizers;
		return correctedReadTrimmed;
	}
	
	*/

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





	AlignmentResult getAlignmentScore(u_int64_t readIndex, const vector<MinimizerType>& s1, const vector<MinimizerType>& s2, const vector<u_int32_t>& s2_pos, const std::unique_ptr<spoa64::AlignmentEngine>& al){
		
		vector<MinimizerType> weights(s1.size(), 1);
		//cout << "-" << endl;
		spoa64::Graph graph{};

		graph.AddAlignment(spoa64::Alignment(), vec32_to_vec64(s1), s1.size(), vec32_to_vec64(weights));
		//spoa64::Alignment alignment = al->Align(s2, s2.size(), graph);

		return getAlignmentScore(readIndex, graph, s1, s2, s2_pos, al);

	}

	AlignmentResult getAlignmentScore(u_int64_t readIndex, const spoa64::Graph& graph, const vector<MinimizerType>& s1, const vector<MinimizerType>& s2, const vector<u_int32_t>& s2_pos, const std::unique_ptr<spoa64::AlignmentEngine>& al){
		

		u_int32_t posStartBps = -1;
		u_int32_t posEndBps = -1;

		AlignmentResult alignmentResult = {readIndex, 0, 0, 0, 0, 0};

		spoa64::Alignment alignment = al->Align(vec32_to_vec64(s2), s2.size(), graph);


		for (size_t i=0; i<alignment.size(); i++) {

			int64_t v1 = alignment[i].first;
			int64_t v2 = alignment[i].second;
			
			if(v1 == -1){ //insert in
				//cout << "Insertion: " << s2[v2] << endl;
				alignmentResult._nbInsertions += 1;
			}
			else if(v2 == -1){ //insert in
				//cout << "Deletion: " << s1[v1] << endl;
				alignmentResult._nbDeletions += 1;
			}
			else if(s1[v1] == s2[v2]){
				alignmentResult._nbMatches += 1;
				if(posStartBps == -1) posStartBps = s2_pos[v2];
				posEndBps = s2_pos[v2];
			}
			else{
				alignmentResult._nbMissmatches += 1;
			}
		}

		if(posStartBps != -1){
			alignmentResult._alignLengthBps = posEndBps - posStartBps;
		}
		
		return alignmentResult;
	}

	/*
	double computeGlobalSimilarity(const vector<MinimizerType>& s1, const vector<MinimizerType>& s2, const std::unique_ptr<spoa64::AlignmentEngine>& al){
		
		vector<MinimizerType> weights(s1.size(), 1);
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
				//score -= 1;
			}
		}

		//cout << "done" << endl;
		double length = max(s1.size(), s2.size());
		return score / length;
	}

	void printAlignmentScore(const vector<MinimizerType>& s1, const vector<MinimizerType>& s2, const std::unique_ptr<spoa64::AlignmentEngine>& al){
		
		vector<MinimizerType> weights(s1.size(), 1);
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

	void printAlignmentScore(const spoa64::Alignment& alignment, const vector<MinimizerType>& s1, const vector<MinimizerType>& s2, const std::unique_ptr<spoa64::AlignmentEngine>& al){
		

		//spoa64::Alignment alignment = al->Align(s2, s2.size(), graph);
		size_t i = 0;
		size_t j = 0;
		for (const auto& it : alignment) {
			int64_t v1 = it.first;
			int64_t v2 = it.second;
			
			if(v1 == -1){ //insert in
				cout << s2[j] << " Insertion" << endl;
				j += 1;
			}
			else if(v2 == -1){ //insert in
				cout << s1[i] << " Deletion" << endl;
				i += 1;
			}
			else if(s1[v1] == s2[v2]){
				cout << s1[i] << " " << s2[j] << " Match" << endl;
				i += 1;
				j += 1;
			}
			else{
				cout << s1[i] << " " << s2[j] << " Missmatch" << endl;
				i += 1;
				j += 1;
			}
		}

	}
	*/
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

	bool isSequenceReversed2(const MinimizerRead& mappedRead, unordered_map<MinimizerType, u_int8_t>& readMinimizer_to_direction){
	
		
		int nbSame = 0;
		int nbDiff = 0;

		for(size_t i=0; i<mappedRead._minimizers.size(); i++){

			MinimizerType minimizer = mappedRead._minimizers[i];
			if(readMinimizer_to_direction.find(minimizer) == readMinimizer_to_direction.end()) continue;

			if(readMinimizer_to_direction[minimizer] == mappedRead._readMinimizerDirections[i]){
				nbSame += 1;
			}
			else{
				nbDiff += 1;
			}

		}

		return nbDiff > nbSame;
	}

	/*
	bool isSequenceReversed(const vector<u_int64_t>& s1, const vector<u_int64_t>& s2, const std::unique_ptr<spoa64::AlignmentEngine>& al){

		int lastMatchpos = 0;
		double score_forward = getAlignmentScore(s1, s2, al, lastMatchpos);

		vector<u_int64_t> s2_rev = s2; 
		std::reverse(s2_rev.begin(), s2_rev.end());
		double score_reverse = getAlignmentScore(s1, s2_rev, al, lastMatchpos);

		//cout << score_forward << " " << score_reverse << endl;

		return score_reverse > score_forward;
		
	}
	*/
	/*
	vector<u_int64_t> performCorrection(spoa64::Graph& graph, const vector<u_int64_t>& readMinimizers, const std::unique_ptr<spoa64::AlignmentEngine>& al, const vector<u_int64_t>& readNodes){
	

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
	/*
	vector<MinimizerType> performCorrection(spoa64::Graph& poaGraph, const vector<MinimizerType>& readMinimizers, const std::unique_ptr<spoa64::AlignmentEngine>& al, const std::unique_ptr<spoa64::AlignmentEngine>& alTrimming, vector<MinimizerType>& correctedCoverages, unordered_set<MinimizerType>& minimizerSet){


		vector<MinimizerType> correctReadMinimizers2 = computePath(poaGraph, poaGraph.nodes_[0].get(), nullptr, maxPathSize, minimizerSet);


		vector<MinimizerType> correctReadMinimizersTrimmed = trimCorrectedPath(readMinimizers, correctReadMinimizers2, alTrimming);
		
		return correctReadMinimizersTrimmed;
	}
	*/
	/*
	void trimSimple(const vector<MinimizerType>& readMinimizers, const vector<u_int8_t>& readQualities, const vector<u_int8_t>& readMinimizerDirections, unordered_set<MinimizerType>& readMinimizersSet, vector<MinimizerType>& readMinimizersTrimmed, vector<u_int8_t>& readQualitiesTrimmed, vector<u_int8_t>& readMinimizerDirectionsTrimmed){

		readMinimizersTrimmed.clear();
		readQualitiesTrimmed.clear();
		readMinimizerDirectionsTrimmed.clear();

		int startPosition = -1;
		int endPosition = -1;

		for(size_t i=0; i<readMinimizers.size(); i++){

			MinimizerType minimizer = readMinimizers[i];

			if(readMinimizersSet.find(minimizer) != readMinimizersSet.end()){
				if(startPosition == -1) startPosition = i;
				endPosition = i;
			}
		}


		if(startPosition == -1) return;

		//startPosition = max(0, startPosition-5);
		//endPosition = min((int)(readMinimizers.size()-1), endPosition+5);

		for(size_t i=startPosition; i<=endPosition; i++){
			readMinimizersTrimmed.push_back(readMinimizers[i]);
			readQualitiesTrimmed.push_back(readQualities[i]);
			readMinimizerDirectionsTrimmed.push_back(readMinimizerDirections[i]);
		}

	}
	*/
	CorrectedRead trimCorrectedPath(const CorrectedRead& read, const CorrectedRead& correctedRead, const std::unique_ptr<spoa64::AlignmentEngine>& al){

		/*
		if(_print_debug){
			cout << "trimming" << endl;

			cout << endl << "Original sequence:" << endl;
			for(MinimizerType m : readMinimizers){
				cout << m << endl;
			}
			cout << endl << "Corrected sequence:" << endl;
			for(size_t i=0; i<correctedReadMinimizers.size(); i++){
				cout << i << ": " << correctedReadMinimizers[i] << endl;
			}
		}
		*/
		
		
		spoa64::Graph graphAln{};

		vector<MinimizerType> weights(correctedRead._minimizers.size(), 1);
		graphAln.AddAlignment(spoa64::Alignment(), correctedRead._minimizers, correctedRead._minimizers.size(), weights);

		spoa64::Alignment alignment = al->Align(read._originalMinimizers, read._originalMinimizers.size(), graphAln);


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
			else{

				//cout << correctedReadMinimizers[v1] << " " << readMinimizers[v2] << endl;

				if(correctedRead._minimizers[v1] == read._originalMinimizers[v2]){

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

		}


		if(startMatchPos == -1 || startMatchPos == endMatchPos){
			return {correctedRead._readIndex, {}, {}, {}};
		}
	
		//cout << startMatchPos << " " << endMatchPos << endl;
		//if(_print_debug){
		//	cout << "\tTrim:" << startMatchPos << " " << endMatchPos << endl;
		//}

		vector<MinimizerType> minimizersTrimmed;
		vector<u_int32_t> minimizersPosTrimmed;
		vector<u_int8_t> qualitiesTrimmed;

		for(size_t i=startMatchPos; i<endMatchPos; i++){
			minimizersTrimmed.push_back(correctedRead._minimizers[i]);
			qualitiesTrimmed.push_back(correctedRead._qualities[i]);
			//minimizersPosTrimmed.push_back(correctedRead._minimizersPos[i]);
		}
		//cout << "check si besoin de < ou <= ici" << endl;
		//vector<MinimizerType>::const_iterator first = correctedRead._minimizers.begin() + startMatchPos;
		//vector<MinimizerType>::const_iterator last = correctedRead._minimizers.begin() + endMatchPos;
		//vector<MinimizerType> correctedReadMinimizersTrimmed(first, last);

		//for(size_t i=0; i<minimizersTrimmed.size(); i++){
		//	cout << i << " " << minimizersTrimmed[i] << endl;
		//}
		//for(size_t i=0; i<correctedReadMinimizersTrimmed.size(); i++){
		//	cout << i << " " << correctedReadMinimizersTrimmed[i] << endl;
		//}
		//getchar();
		return {correctedRead._readIndex, minimizersTrimmed, qualitiesTrimmed, correctedRead._readMinimizerDirections, correctedRead._originalMinimizers, correctedRead._originalQualities};
	}

	void trimCorrectedPathFast(const spoa64::Graph& graph, const vector<MinimizerType>& readMinimizers, const std::vector<MinimizerType>& correctedReadMinimizers, const vector<u_int8_t>& readQualities, const std::unique_ptr<spoa64::AlignmentEngine>& al, vector<MinimizerType>& readMinimizersTrimmed, vector<u_int8_t>& readQualitiesTrimmed){

		
		readMinimizersTrimmed.clear();
		readQualitiesTrimmed.clear();
		//spoa64::Graph graphAln{};

		//vector<MinimizerType> weights(correctedReadMinimizers.size(), 1);
		//graphAln.AddAlignment(spoa64::Alignment(), correctedReadMinimizers, correctedReadMinimizers.size(), weights);

		spoa64::Alignment alignment = al->Align(correctedReadMinimizers, correctedReadMinimizers.size(), graph);


		int startMatchPos = -1;
		int endMatchPos = -1;

		for (const auto& it : alignment) {
			int64_t v1 = it.first;
			int64_t v2 = it.second;

			if(v1 == -1){ //insert in
				//cout << "insertion" << endl;
			}
			else if(v2 == -1){ //insert in
				//cout << "deletion" << endl;
			}
			else{

				//cout << correctedReadMinimizers[v1] << " " << readMinimizers[v2] << endl;

				if(readMinimizers[v1] == correctedReadMinimizers[v2]){

					//cout << "match" << endl;
					if(startMatchPos == -1){
						startMatchPos = v2;
					}
					
					endMatchPos = v2+1;
				}
				else{
					//cout << "missmatch" << endl;
				}

			}

		}


		if(startMatchPos == -1 || startMatchPos == endMatchPos){
			return;
		}
	
		startMatchPos = max(0, startMatchPos-1);
		endMatchPos = min((int)correctedReadMinimizers.size(), endMatchPos+1);
		//cout << startMatchPos << " " << endMatchPos << endl;

		vector<MinimizerType>::const_iterator first = correctedReadMinimizers.begin() + startMatchPos;
		vector<MinimizerType>::const_iterator last = correctedReadMinimizers.begin() + endMatchPos;
		vector<MinimizerType> correctedReadMinimizersTrimmed(first, last);

		vector<u_int8_t>::const_iterator first2 = readQualities.begin() + startMatchPos;
		vector<u_int8_t>::const_iterator last2 = readQualities.begin() + endMatchPos;
		vector<u_int8_t> quals(first2, last2);

		readMinimizersTrimmed = correctedReadMinimizersTrimmed;
		readQualitiesTrimmed = quals;

	}

	//struct PathSuccessor{
	//	u_int64_t _node;
	//	u_int32_t _distance;
	//};


	struct DereplicatedEdge{
		spoa64::Graph::Edge* _edge;
		u_int64_t _weight;
	};
	
	vector<MinimizerType> trimByAbundance(const vector<MinimizerType>& readMinimizers, const vector<MinimizerRead>& mappedReads){



		unordered_map<KmerVec, u_int32_t> kminmer_to_abundance;

		vector<ReadKminmerComplete> kminmersInfos;
		vector<u_int32_t> minimizersPos(readMinimizers.size(), 0);
		vector<u_int8_t> minimizersQualities(readMinimizers.size(), 0);
		MDBG::getKminmers_complete(_kminmerSize, readMinimizers, minimizersPos, kminmersInfos, 0, minimizersQualities);

		for(const ReadKminmerComplete& kminmerInfo : kminmersInfos){

			const KmerVec& vec = kminmerInfo._vec;
			kminmer_to_abundance[vec] = 1;

			
		}

		
		
		for(size_t i=0; i<mappedReads.size(); i++){

			const MinimizerRead& minimizerRead = mappedReads[i];

			vector<ReadKminmerComplete> kminmersInfos;
			vector<u_int32_t> minimizersPos(minimizerRead._minimizers.size(), 0);
			vector<u_int8_t> minimizersQualities(minimizerRead._minimizers.size(), 0);
			MDBG::getKminmers_complete(_kminmerSize, minimizerRead._minimizers, minimizersPos, kminmersInfos, 0, minimizersQualities);

			for(const ReadKminmerComplete& kminmerInfo : kminmersInfos){

				const KmerVec& vec = kminmerInfo._vec;

				if(kminmer_to_abundance.find(vec) == kminmer_to_abundance.end()) continue;

				kminmer_to_abundance[vec] += 1;

				
			}

		}

		vector<float> abundances;
		for(const auto& it : kminmer_to_abundance){
			abundances.push_back(it.second);
		}

		float median = Utils::compute_median_float(abundances);
		float minAbundance = median * 0.2;

		//cout << "\t----" << endl;
		//cout << "\tMedian: " << median << endl;

		size_t i=0;
		size_t startPosition = -1;
		size_t endPosition = -1;

		for(const ReadKminmerComplete& kminmerInfo : kminmersInfos){

			const KmerVec& vec = kminmerInfo._vec;

			if(kminmer_to_abundance[vec] >= minAbundance){
				if(startPosition == -1) startPosition = i;
				endPosition = i+1;
			}
			//cout << "\t\t" << kminmer_to_abundance[vec] << endl;

			i += 1;
		}

		if(startPosition == -1) return readMinimizers;

		vector<MinimizerType> readMinimizersTrimmed;

		for(size_t i=startPosition; i<=endPosition; i++){
			readMinimizersTrimmed.push_back(readMinimizers[i]);
		}

		return readMinimizersTrimmed;
	}
	





	CorrectedRead computePath3(const CorrectedRead& read, Graph* graph, unordered_set<MinimizerType>& readMinimizers){


		vector<Node*> correctedPath;
		//vector<u_int8_t> pathQualities;
		//vector<u_int32_t> pathSupport;

		unordered_set<Node*> isNodeVisited;
		Node* currentNode = graph->_nodes[0];

		isNodeVisited.insert(currentNode);

		correctedPath.push_back(currentNode);	
		//pathQualities.push_back(currentNode->_maxQuality);


		while (true) {
	
			const vector<Node*>& subPath = findPathFromNode(currentNode, readMinimizers);

			if(subPath.size() == 0){
				MinimizerType currentNodeIndex = currentNode->_nodeIndex;
				MinimizerType nextNodeIndex = currentNodeIndex+1;

				if(nextNodeIndex >= readMinimizers.size()) break;

				currentNode = graph->_nodes[nextNodeIndex];

				if(isNodeVisited.find(currentNode) != isNodeVisited.end()) break;
				isNodeVisited.insert(currentNode);

				correctedPath.push_back(currentNode);
			}
			else{

				bool isVisited = false;

				for(Node* node : subPath){

					if(isNodeVisited.find(node) != isNodeVisited.end()){
						isVisited = true;
						break;
					}
					isNodeVisited.insert(node);

					correctedPath.push_back(node);
				}

				if(isVisited) break;

				currentNode = subPath[subPath.size()-1];
			}



		}


		vector<MinimizerType> path;
		vector<u_int8_t> pathQualities;
		//vector<u_int32_t> pathSupport;

		for(Node* node : correctedPath){
			path.push_back(node->_minimizer);
			pathQualities.push_back(node->_maxQuality);
		}

		/*
		if(_print_debug){

			cout << "\tPath support: " << endl;
			for(size_t i=0; i<path.size(); i++){
				if(i < pathSupport.size()){
					cout << "\t\t" << path[i] << " " << pathSupport[i] << endl;
				}
				else{
					cout << "\t\t" << path[i] << endl;
				}
			}
			
		}
		*/
		
		//vector<MinimizerType> solidPath;

		//if(firstSolidPosition == -1){
		//	return solidPath;
		//}


		//for(size_t i=firstSolidPosition; i<=lastSolidPosition; i++){
		//	solidPath.push_back(path[i]);
		//}
		
		return {read._readIndex, path, pathQualities, read._readMinimizerDirections, read._originalMinimizers, read._originalQualities};
	}


	vector<Node*> findPathFromNode(Node* fromNode, unordered_set<MinimizerType>& readMinimizers){

		vector<Node*> path;
		//vector<u_int8_t> pathQualities;
		//vector<u_int32_t> pathSupport;

		unordered_set<Node*> isNodeVisited;
		Node* currentNode = fromNode;

		isNodeVisited.insert(currentNode);

		//path.push_back(currentNode);	
		//pathQualities.push_back(currentNode->_maxQuality);

		//size_t i=1;

		while (true) {
	
			int maxWeight = 0;

			for(Edge* successor : currentNode->_successors){
				if(successor->_head == currentNode) continue; //prevent self repeat node

				if(successor->_weight > maxWeight){
					maxWeight = successor->_weight;
				}
			}


			Edge* maxSuccessor = nullptr;
			float minWeight = maxWeight * 0.75;
			

			vector<Edge*> solidSuccessors;
			for(Edge* successor : currentNode->_successors){
				if(successor->_head == currentNode) continue; //prevent self repeat node

				if(successor->_weight >= minWeight){
					solidSuccessors.push_back(successor);
				}
			}

			if(solidSuccessors.size() == 0){
				break;
			}
			else if(solidSuccessors.size() == 1){
				maxSuccessor = solidSuccessors[0];
			}
			else{
				break;
				//maxSuccessor = computeBestSuccessor(solidSuccessors, readMinimizers);
			}
			//cout << "\tVisit: " << currentNode->id << "-" << graph.decoder(currentNode->code) << endl;

			//vector<DereplicatedEdge> successors = getSolidSuccessors(graph, currentNode, readMinimizers);

	


			//for(const DereplicatedEdge& succ : successors){
			//	cout << "\t\t" << succ._edge->head->id << "-" << graph.decoder(succ._edge->head->code) << " " << succ._edge->weight << endl;

			//}

			//if(maxSuccessor == nullptr) break;
			//if(successors.size() == 0) return nullPath;
			//if(successors.size() > 2 && successors[0]._weight <= 2) return path;

			currentNode = maxSuccessor->_head;

			if(isNodeVisited.find(currentNode) != isNodeVisited.end()) break;
			isNodeVisited.insert(currentNode);

			//if(currentNode->_abundance > 2){
			//if(maxSuccessor->_support > 2){
			//	if(firstSolidPosition == -1) firstSolidPosition = i-1;
			//	lastSolidPosition = i;
			//}
			//if(currentNode->_abundance > 2){
			//	if(firstSolidPosition == -1) firstSolidPosition = i;
			//	lastSolidPosition = i;
			//}

			path.push_back(currentNode);
			//pathQualities.push_back(currentNode->_maxQuality);
			//pathSupport.push_back(maxSuccessor->_support);

			if(readMinimizers.find(currentNode->_minimizer) != readMinimizers.end()){
				return path;
			}
			//nodePath.push_back(currentNode);
			//if(currentNode == endNode) break;
			//if(path.size() > maxSize) return nullPath;
			//cout << currentNode->toString() << endl;
			//getchar();

			//i += 1;
		}

		path.clear();
		return path;
	}

	CorrectedRead computePath2(const CorrectedRead& read, Graph* graph, unordered_set<MinimizerType>& minimizerSet){


		vector<MinimizerType> path;
		vector<u_int8_t> pathQualities;
		vector<u_int32_t> pathSupport;

		unordered_set<Node*> isNodeVisited;
		Node* currentNode = graph->_nodes[0];

		isNodeVisited.insert(currentNode);


		//int firstSolidPosition = -1;
		//int lastSolidPosition = -1;

		//if(currentNode->_abundance > 2){
		//	firstSolidPosition = 0;
		//}


		path.push_back(currentNode->_minimizer);	
		pathQualities.push_back(currentNode->_maxQuality);

		size_t i=1;

		while (true) {
	
			int maxWeight = 0;

			for(Edge* successor : currentNode->_successors){
				if(successor->_head == currentNode) continue; //prevent self repeat node

				if(successor->_weight > maxWeight){
					maxWeight = successor->_weight;
				}
			}


			Edge* maxSuccessor = nullptr;
			float minWeight = maxWeight * 0.75;
			

			vector<Edge*> solidSuccessors;
			for(Edge* successor : currentNode->_successors){
				if(successor->_head == currentNode) continue; //prevent self repeat node

				if(successor->_weight >= minWeight){
					solidSuccessors.push_back(successor);
				}
			}

			if(solidSuccessors.size() == 0){
				break;
			}
			else if(solidSuccessors.size() == 1){
				maxSuccessor = solidSuccessors[0];
			}
			else{
				maxSuccessor = computeBestSuccessor(solidSuccessors, minimizerSet);
			}
			//cout << "\tVisit: " << currentNode->id << "-" << graph.decoder(currentNode->code) << endl;

			//vector<DereplicatedEdge> successors = getSolidSuccessors(graph, currentNode, readMinimizers);

	


			//for(const DereplicatedEdge& succ : successors){
			//	cout << "\t\t" << succ._edge->head->id << "-" << graph.decoder(succ._edge->head->code) << " " << succ._edge->weight << endl;

			//}

			//if(maxSuccessor == nullptr) break;
			//if(successors.size() == 0) return nullPath;
			//if(successors.size() > 2 && successors[0]._weight <= 2) return path;

			currentNode = maxSuccessor->_head;

			if(isNodeVisited.find(currentNode) != isNodeVisited.end()) break;
			isNodeVisited.insert(currentNode);

			//if(currentNode->_abundance > 2){
			//if(maxSuccessor->_support > 2){
			//	if(firstSolidPosition == -1) firstSolidPosition = i-1;
			//	lastSolidPosition = i;
			//}
			//if(currentNode->_abundance > 2){
			//	if(firstSolidPosition == -1) firstSolidPosition = i;
			//	lastSolidPosition = i;
			//}

			path.push_back(currentNode->_minimizer);
			pathQualities.push_back(currentNode->_maxQuality);
			pathSupport.push_back(maxSuccessor->_support);
			//nodePath.push_back(currentNode);
			//if(currentNode == endNode) break;
			//if(path.size() > maxSize) return nullPath;
			//cout << currentNode->toString() << endl;
			//getchar();

			i += 1;
		}

		if(_print_debug){

			cout << "\tPath support: " << endl;
			for(size_t i=0; i<path.size(); i++){
				if(i < pathSupport.size()){
					cout << "\t\t" << path[i] << " " << pathSupport[i] << endl;
				}
				else{
					cout << "\t\t" << path[i] << endl;
				}
			}
			
		}
		
		//vector<MinimizerType> solidPath;

		//if(firstSolidPosition == -1){
		//	return solidPath;
		//}


		//for(size_t i=firstSolidPosition; i<=lastSolidPosition; i++){
		//	solidPath.push_back(path[i]);
		//}
		
		return {read._readIndex, path, pathQualities, read._readMinimizerDirections, read._originalMinimizers, read._originalQualities};
	}

	
	Edge* computeBestSuccessor(const vector<Edge*>& solidSuccessors, unordered_set<MinimizerType>& readMinimizers){

		u_int64_t maxCompletion = 0;
		Edge* maxSuccessor = nullptr;

		for(Edge* successor : solidSuccessors){
			u_int64_t completion = computeSuccessorCompletion(successor, readMinimizers);

			if(completion > maxCompletion){
				maxCompletion = completion;
				maxSuccessor = successor;
			}
		}

		return maxSuccessor;
	}

	u_int64_t computeSuccessorCompletion(Edge* successor, unordered_set<MinimizerType>& readMinimizers){


		list<Node*> queue;

		queue.push_back(successor->_head);
		unordered_set<Node*> isVisited;
	
		u_int64_t completion = successor->_weight;

		while (!queue.empty()) {
	
			Node* currentNode = queue.front();
			queue.pop_front();

			if(isVisited.find(currentNode) != isVisited.end()) continue;
			isVisited.insert(currentNode);


			for(Edge* nn : currentNode->_successors){

				if(readMinimizers.find(nn->_head->_minimizer) != readMinimizers.end()){
					completion += nn->_weight;
				}

				queue.push_back(nn->_head);
			}


		}

		return completion;
	}

	/*

	vector<MinimizerType> computePath(const spoa64::Graph& graph, spoa64::Graph::Node* startNode, spoa64::Graph::Node* endNode, u_int64_t maxSize, unordered_set<MinimizerType>& readMinimizers){

		if(_print_debug){
			//cout << "Compute path: " << startNode->id << "-" << graph.decoder(startNode->code) << endl; //<< " -> " << endNode->id << "-" << graph.decoder(endNode->code) << endl;
		}

		//vector<spoa64::Graph::Node*> nodePath;
		vector<MinimizerType> path;
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
		//nodePath.push_back(currentNode);
		
		while (true) {
	

			//cout << "\tVisit: " << currentNode->id << "-" << graph.decoder(currentNode->code) << endl;

			vector<DereplicatedEdge> successors = getSolidSuccessors(graph, currentNode, readMinimizers);

	


			//for(const DereplicatedEdge& succ : successors){
			//	cout << "\t\t" << succ._edge->head->id << "-" << graph.decoder(succ._edge->head->code) << " " << succ._edge->weight << endl;

			//}

			if(successors.size() == 0) break;
			//if(successors.size() == 0) return nullPath;
			//if(successors.size() > 2 && successors[0]._weight <= 2) return path;

			currentNode = successors[0]._edge->head;
			path.push_back(graph.decoder(currentNode->code));
			//nodePath.push_back(currentNode);
			//if(currentNode == endNode) break;

			//if(path.size() > maxSize) return nullPath;
			
		}

		//cout << "----- HAHAHA" << endl;
		//for(spoa64::Graph::Node* node : nodePath){
		//	cout << "\t " << graph.decoder(node->code) << " " << node->Coverage() << endl;
		//}

		
		return path;
	}

	vector<DereplicatedEdge> getSolidSuccessors(const spoa64::Graph& graph, spoa64::Graph::Node* node, unordered_set<MinimizerType>& readMinimizers){

		vector<DereplicatedEdge> successors = getDereplicatedSuccessors(graph, node, readMinimizers);

		std::sort(successors.begin(), successors.end(), [](const DereplicatedEdge& a, const DereplicatedEdge& b){
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

		vector<DereplicatedEdge> solidEdges;

		for(const DereplicatedEdge& edge : successors){
			//if(edge._weight == 1 && edge._edge->head->Coverage() <= 1) continue;

			if(edge._weight < minWeight) continue;
				
			solidEdges.push_back(edge);
		}

		return solidEdges;
	}


	struct SuccessorCompletion{
		spoa64::Graph::Edge* _edge;
		u_int64_t _completion;
	};

	vector<DereplicatedEdge> getDereplicatedSuccessors(const spoa64::Graph& graph, spoa64::Graph::Node* node, unordered_set<MinimizerType>& readMinimizers){

		for(spoa64::Graph::Edge* edge : node->outedges){
			//if(_print_debug) cout << "\t\tPossible successor: " << edge->head->id << " " << edge->weight << " " << computeEdgeCompletion(graph, edge, readMinimizers) << endl;
			//for(spoa64::Graph::Edge* edge2 : edge->head->outedges){
			//	cout << "\t\t\tLala: " << edge2->head->id << " " << edge2->weight << " " << computeEdgeCompletion(graph, edge2, readMinimizers) << endl;
			//}
		}

		unordered_map<MinimizerType, SuccessorCompletion> minimizer_to_edge;
		unordered_map<MinimizerType, u_int64_t> minimizer_to_weight;

		for(spoa64::Graph::Edge* edge : node->outedges){

			MinimizerType minimizer = edge->head->code;

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

			MinimizerType minimizer = it.first;
			spoa64::Graph::Edge* edge = it.second._edge;
			u_int64_t completion = it.second._completion;

			u_int64_t derepWeight = minimizer_to_weight[minimizer];

			//edge->weight = minimizer_to_weight[it.first];
			derepEdges.push_back({edge, derepWeight}); //+completion
			//cout << "\t\tPossible successor derep: " << edge->head->id << " " << derepWeight << " " << it.second._completion << endl;
		}

		return derepEdges;

	}


	u_int64_t computeEdgeCompletion(const spoa64::Graph& graph, spoa64::Graph::Edge* edge, unordered_set<MinimizerType>& readMinimizers){


		list<spoa64::Graph::Node*> queue;

		queue.push_back(edge->head);
		unordered_set<MinimizerType> isVisited;
	
		u_int64_t nbVisitedReadNodes = edge->weight;

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
	*/
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

	void writeRead(u_int32_t readIndex, const vector<MinimizerType>& minimizers, const vector<u_int8_t>& qualities){

		//static u_int64_t maxHashValue = -1;
		//static u_int64_t minimizerBound = 0.005 * maxHashValue;

		//#pragma omp critical(dataupdate)
		#pragma omp critical
		{
			

			//vector<u_int64_t> validMinimizers;
			//for(u_int64_t m : minimizers){
			//	if(m > minimizerBound) continue;
			//	validMinimizers.push_back(m);
			//}

			
			_readWriterQueue.push({readIndex, minimizers, qualities});

			while(!_readWriterQueue.empty()){


				const ReadWriter& readWriter = _readWriterQueue.top();

				//cout << readWriter._readIndex << " " << _nextReadIndexWriter << " " << _readWriterQueue.size() << endl;

				if(readWriter._readIndex == _nextReadIndexWriter){

					//if(readWriter._minimizers.size() > 0){

						for(u_int64_t m : readWriter._minimizers){
							_correctionCheckSum += m*readWriter._minimizers.size();
						}

						u_int32_t size = readWriter._minimizers.size();
						_file_readData.write((const char*)&size, sizeof(size));

						u_int8_t isCircular = CONTIG_LINEAR;
						_file_readData.write((const char*)&isCircular, sizeof(isCircular));

						_file_readData.write((const char*)&readWriter._minimizers[0], size*sizeof(MinimizerType));

						vector<u_int32_t> minimizerPosDummy(size+1, 0);
						vector<u_int8_t> minimizerDirectionDummy(size, 0);
						_file_readData.write((const char*)&minimizerPosDummy[0], (size+1)*sizeof(u_int32_t));
						_file_readData.write((const char*)&minimizerDirectionDummy[0], size*sizeof(u_int8_t));
						_file_readData.write((const char*)&readWriter._qualities[0], size*sizeof(u_int8_t));
					//}

					if(readWriter._readIndex % 10000 == 0){
						cout << readWriter._readIndex << " " << _correctionCheckSum << endl;
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

	void savePoaGraph(const spoa64::Graph& graph, const string& outputFilename, const vector<MinimizerType>& readNodes){

		cout << "save poa graph" << endl;

		ofstream outputFile(outputFilename);

		ofstream colorFile(outputFilename + ".color.csv");
		colorFile << "Name,Color" << endl;
		
		ofstream edgeFile(outputFilename + ".edge.csv");
		edgeFile << "Name,Edge" << endl;

		//std::cout << "H\tVN:Z:1.0" << std::endl;
		for (const auto& it : graph.nodes()) {

			MinimizerType minimizer = graph.decoder(it->code);
			string id = to_string(it->id) + "-" + to_string(graph.decoder(it->code));

			//if(std::find(readNodes.begin(), readNodes.end(), it->id) != readNodes.end()){
			//if(illuminaMinimizers.find(minimizer) != illuminaMinimizers.end()){
			//	colorFile << (id) << ",green" << endl; 
			//}
			//else{
			//	colorFile << (id) << ",grey" << endl; 
			//}

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


	struct ContigMatch{
		u_int32_t _contigIndex;
		bool _isReversed;
		u_int32_t _coverage;
		u_int32_t _contigStartPosition;
		float _identity;
	};

	unordered_map<u_int32_t, ContigMatch> _readIndex_to_bestContigMatch;

	void loadReadToContigMapping(){

		if(!fs::exists(_inputDir + "/readIndex_to_contigIndex.bin")) return;

		ifstream file(_inputDir + "/readIndex_to_contigIndex.bin");

		while (true) {

			u_int32_t readIndex;
			u_int32_t contigIndex;
			bool isReversed;
			u_int32_t contigCoverage;
			u_int32_t contigStartPosition;
			float identity;

			file.read((char*)&readIndex, sizeof(readIndex));

			if(file.eof()) break;

			file.read((char*)&contigIndex, sizeof(contigIndex));
			file.read((char*)&isReversed, sizeof(isReversed));
			file.read((char*)&contigCoverage, sizeof(contigCoverage));
			file.read((char*)&contigStartPosition, sizeof(contigStartPosition));
			file.read((char*)&identity, sizeof(identity));

			_readIndex_to_bestContigMatch[readIndex] = {contigIndex, isReversed, contigCoverage, contigStartPosition, identity};
			//cout << readIndex << " " << contigIndex << " " << isReversed << " " << contigCoverage << endl;
		}

		file.close();

	}

	vector<u_int64_t> collectReadWithSimilarContigPosition(u_int64_t readIndex){

		vector<u_int64_t> readIndexes;

		if(_readIndex_to_bestContigMatch.find(readIndex) == _readIndex_to_bestContigMatch.end()) return readIndexes;

		const ContigMatch& contigMatch = _readIndex_to_bestContigMatch[readIndex];

		for(const auto& it : _readIndex_to_bestContigMatch){
			u_int64_t otherReadIndex = it.first;
			const ContigMatch& otherContigMatch = it.second;

			if(otherContigMatch._contigIndex != contigMatch._contigIndex) continue;

			int32_t positionDistance = abs((int32_t)contigMatch._contigStartPosition - (int32_t)otherContigMatch._contigStartPosition);
			
			if(positionDistance < 2000){
				readIndexes.push_back(otherReadIndex);
			}
		}

		return readIndexes;
	}


	unordered_map<u_int32_t, vector<ReadMatch>> _readIndex_to_bestReadMatch;

	void loadReadToReadMapping(){

		if(!fs::exists(_inputDir + "/readIndex_to_readIndex.bin")) return;

		ifstream file(_inputDir + "/readIndex_to_readIndex.bin");

		while (true) {

			u_int32_t readIndex;
			u_int32_t nbMatches;
			vector<ReadMatch> matches;

			file.read((char*)&readIndex, sizeof(readIndex));

			if(file.eof()) break;

			file.read((char*)&nbMatches, sizeof(nbMatches));

			matches.resize(nbMatches);

			file.read((char*)&matches[0], matches.size()*sizeof(ReadMatch));

			_readIndex_to_bestReadMatch[readIndex] = matches;
		}

		file.close();

	}

	vector<vector<MinimizerType>> _contigs;

	void loadContigs(){
		if(!fs::exists(_inputDir + "/contig_data_init.txt")) return;
		KminmerParserParallel parser(_inputDir + "/contig_data_init.txt", _minimizerSize, _kminmerSize, false, true, 1);
		parser.parseSequences(LoadContigFunctor(*this));
	}

	class LoadContigFunctor {

		public:

		ReadCorrection& _parent;

		LoadContigFunctor(ReadCorrection& parent) : _parent(parent){
		}

		LoadContigFunctor(const LoadContigFunctor& copy) : _parent(copy._parent){
		}

		~LoadContigFunctor(){
		}


		void operator () (const KminmerList& kminmerList) {

			u_int32_t readIndex = kminmerList._readIndex;
			if(readIndex % 100000 == 0) cout << readIndex << endl;

			_parent._contigs.push_back({kminmerList._readMinimizers});


		}
	};


	CorrectedRead evaluateCorrection(u_int32_t readIndex, const CorrectedRead& correctedRead, const std::unique_ptr<spoa64::AlignmentEngine>& al){
		//u_int32_t readIndex = kminmerList._readIndex;


		CorrectedRead fakeCorrectedReads = {correctedRead._readIndex, correctedRead._originalMinimizers, correctedRead._originalQualities, correctedRead._readMinimizerDirections, correctedRead._originalMinimizers, correctedRead._originalQualities};

		//int readLength = kminmerList._readMinimizers1pos[kminmerList._readMinimizers1pos.size()-1];
		//if(readLength < MIN_NB_MINIMIZER_REFERENCE || readMinimizersOriginal.size() < 4) return;

		//if(_print_debug){
		//	cout << endl << endl;
		//	cout << "-------------------" << endl;
		//	cout << "Read index: " << readIndex << endl;
		//	cout << "Read size original:  " << readMinimizersOriginal.size() << endl;
		//	cout << "Read size corrected: " << readMinimizersCorrected.size() << endl;
		//}

		//return;

		if(_readIndex_to_bestContigMatch.find(readIndex) == _readIndex_to_bestContigMatch.end()) return fakeCorrectedReads;

		const ContigMatch& contig = _readIndex_to_bestContigMatch[readIndex];
		u_int32_t contigIndex = contig._contigIndex;
		bool isContigReversed = contig._isReversed;
		u_int32_t contigCoverage = contig._coverage;
		float identity = contig._identity;
		vector<MinimizerType> contigMinimizers = _contigs[contigIndex];
		if(contigCoverage < 20){
			cout << "Skip low coverage contig" << endl;
			return fakeCorrectedReads;
		}
		//if(identity < 0.98){
			//cout << "Skip low identity read" << endl;
		//	return false;
		//}

		//if(_print_debug){
		//	cout << "Contig index: " << contigIndex << " " << contigMinimizers.size() << " " << isContigReversed << endl;
		//	cout << "Contig coverage: " << contigCoverage << endl;
		//}

		if(isContigReversed){
			//std::reverse(readMinimizersOriginal.begin(), readMinimizersOriginal.end());
			//std::reverse(readMinimizersCorrected.begin(), readMinimizersCorrected.end());
			std::reverse(contigMinimizers.begin(), contigMinimizers.end());
		}

		int nbMatchesCorrected = 0;
		int nbErrorsCorrected = 0;
		int contigStartCorrected = 0;
		int contigEndCorrected = 0;
		int contigLengthCorrected = 0;

		int nbMatches = 0;
		int nbErrors = 0;
		int contigStart = 0;
		int contigEnd = 0;
		int contigLength = 0;
		//if(_print_debug) cout << endl << "Original alignment:" << endl;
		computeAccuracy(correctedRead._originalMinimizers, contigMinimizers, false, al, contigLength, contigStart, contigEnd, nbMatches, nbErrors);

		if(correctedRead._minimizers.size() > 0){
			//if(_print_debug) cout << endl << "Corrected alignment:" << endl;
			computeAccuracy(correctedRead._minimizers, contigMinimizers, true, al, contigLengthCorrected, contigStartCorrected, contigEndCorrected, nbMatchesCorrected, nbErrorsCorrected);
		}

		#pragma omp atomic
		_eval_nbBases += contigLength;

		#pragma omp atomic
		_eval_nbMatches += nbMatchesCorrected;

		//cout << "Correction accuracy: " << (long double) _eval_nbMatches / (long double) _eval_nbBases << endl;

		int correctedLengthDiff = abs((int64_t)contigLength - (int64_t)correctedRead._minimizers.size());


		
		if(nbErrorsCorrected > 5 || correctedLengthDiff > 5){


			cout << endl << endl;
			cout << "-------------------" << endl;
			cout << "Read index: " << readIndex << endl;
			cout << "Read size original:  " << correctedRead._originalMinimizers.size() << endl;
			cout << "Read size corrected: " << correctedRead._minimizers.size() << endl;
			cout << "Contig coverage: " << contigCoverage << endl;
			cout << "Read identity: " << identity << endl;
			cout << "Nb Errors original: " << nbErrors << endl;
			cout << "Nb Errors corrected: " << nbErrorsCorrected << endl;

			unordered_set<MinimizerType> subContigMinimizers;
			for(size_t i=contigStart; i<=contigEnd; i++){
				u_int64_t minimizer = contigMinimizers[i];
				subContigMinimizers.insert(minimizer);
			}
			
			//cout << "Contig index: " << contigIndex << " " << contigMinimizers.size() << " " << isContigReversed << endl;
			//cout << "Contig coverage: " << contigCoverage << endl;

			cout << endl << "\tRead original: " << endl;
			for(size_t i=0; i<correctedRead._originalMinimizers.size(); i++){

				if(subContigMinimizers.find(correctedRead._originalMinimizers[i]) == subContigMinimizers.end()){
					cout << "\t" << correctedRead._originalMinimizers[i] << " " << (int) correctedRead._originalQualities[i] << endl;
				}
				else{
					cout << "\t" << correctedRead._originalMinimizers[i] << " " << (int) correctedRead._originalQualities[i] << " *" << endl;
				}
			}

			cout << endl << "\tRead corrected: " << endl;
			for(size_t i=0; i<correctedRead._minimizers.size(); i++){

				if(subContigMinimizers.find(correctedRead._minimizers[i]) == subContigMinimizers.end()){
					cout << "\t" << correctedRead._minimizers[i] << endl;
				}
				else{
					cout << "\t" << correctedRead._minimizers[i] << " *" << endl;
				}
			}

			cout << endl << "\tContig sequences: " << endl;
			for(size_t i=contigStart; i<=contigEnd; i++){

				u_int64_t minimizer = contigMinimizers[i];
				
				if(std::find(correctedRead._minimizers.begin(), correctedRead._minimizers.end(), minimizer) == correctedRead._minimizers.end()){

					cout << "\t\t" << minimizer << endl;
				}
				else{
					cout << "\t\t" << minimizer << " *"<< endl;
				}
			}

			_debugDbgGraph->save("/pasteur/zeus/projets/p02/seqbio/gbenoit/workspace/home/tmp/poaGraph.gfa", correctedRead._originalMinimizers);

			cout << "Read index: " << readIndex << endl;
			cout << "Read size original:  " << correctedRead._originalMinimizers.size() << endl;
			cout << "Read size corrected: " << correctedRead._minimizers.size() << endl;
			cout << "Contig coverage: " << contigCoverage << endl;
			cout << "Read identity: " << identity << endl;
			cout << "Nb Errors original: " << nbErrors << endl;
			cout << "Nb Errors corrected: " << nbErrorsCorrected << endl;
			cout << "Correction accuracy: " << (long double) _eval_nbMatches / (long double) _eval_nbBases << endl;

			getchar();
		}
		

		//if(_print_debug){
		//	cout << "ReadIndex: " << readIndex << endl;
		//	cout << readMinimizersOriginal.size() << " " << readMinimizersCorrected.size() << endl;
		//}


		delete _debugDbgGraph;
		
		vector<MinimizerType> fakeMinimizers;
		vector<u_int8_t> fakeQualities;
		for(size_t i=contigStart; i<=contigEnd; i++){
			fakeMinimizers.push_back(contigMinimizers[i]);
			fakeQualities.push_back(40);
		}

		fakeCorrectedReads._minimizers = fakeMinimizers;
		fakeCorrectedReads._qualities = fakeQualities;

		return fakeCorrectedReads;
		//if(readIndex == 46) getchar();
	}


	void computeAccuracy(const vector<MinimizerType>& readMinimizers, const vector<MinimizerType>& contigMinimizers, bool isCorrectedRead, const std::unique_ptr<spoa64::AlignmentEngine>& al, int& contigLength, int& firstMatchPos, int& lastMatchPos, int& nbMatches, int& nbErrors){

		nbMatches = 0;
		nbErrors = 0;

		firstMatchPos = -1;
		lastMatchPos = -1;

		contigLength = 0;
		//int startPos = -1;
		//int endPos = -1;

		vector<MinimizerType> weights(readMinimizers.size(), 1);
		
		spoa64::Graph graph{};

		graph.AddAlignment(spoa64::Alignment(), readMinimizers, readMinimizers.size(), weights);
		spoa64::Alignment alignment = al->Align(contigMinimizers, contigMinimizers.size(), graph);



		for (const auto& it : alignment) {

			int64_t v1 = it.first;
			int64_t v2 = it.second;
			
			
			//if(startPos == -1) startPos = v1;
			//endPos = v1;
			//cout << "\t" << v1 << " " << v2 << endl;

			if(v1 == -1){ //insert in
				
				//nbMatches -= 1;
				nbErrors += 1;
				contigLength += 1;
				if(_print_debug) cout << "\tInsertion " << endl;

			}
			else if(v2 == -1){ //insert in

				nbErrors += 1;
				if(_print_debug) cout << "\tDeletion  "  << endl;
			}
			else if(readMinimizers[v1] == contigMinimizers[v2]){

				//cout << "Match: " << readMinimizers[v1] << " " << contigMinimizers[v2] << endl;
				if(firstMatchPos == -1) firstMatchPos = v2;
				lastMatchPos = v2;

				contigLength += 1;
				nbMatches += 1;
				if(_print_debug) cout << "\tMatch     " << readMinimizers[v1] << " " << contigMinimizers[v2] << endl;

			}
			else{

				//nbMatches -= 1;

				nbErrors += 1;
				contigLength += 1;
				if(_print_debug) cout << "\tMissmatch " << readMinimizers[v1] << " " << contigMinimizers[v2] << endl;

			}
		}

		//int nbStartingMissmatches = startPos;
		//int nbEndingMissmatches = readMinimizers.size() - 1 - endPos;

		//cout << nbStartingMissmatches << " " << nbEndingMissmatches << endl;

		//contigLength += nbStartingMissmatches;
		//contigLength += nbEndingMissmatches;

		//cout << contigLength << endl;
		//cout << "lala" << firstMatchPos << " " << lastMatchPos << endl;
		/*
		if(isCorrectedRead && nbErrors > 5){
			cout << "Nb errors: " << nbErrors << endl;
			//getchar();
			return true;
		}

		return false;
		*/
	}


};	


#endif 


