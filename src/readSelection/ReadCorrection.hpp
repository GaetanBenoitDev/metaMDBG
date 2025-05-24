
//metaMDBG_testOv1

#ifndef MDBG_METAG_READCORRECTION
#define MDBG_METAG_READCORRECTION

#include "../Commons.hpp"
//#include "../graph/GraphPOA.hpp"
//#include "../utils/spoa64/include/spoa64/spoa.hpp"
#include "../graph/GfaParser.hpp"
//#include "../utils/BloomFilter.hpp"
#include "../readSelection/MinimizerChainer.hpp"
//#include "../readSelection/ReadMapper.hpp"

/*


./bin/metaMDBG asm --out-dir ~/appa/run/correction/nanopore_AD_circ1_asm/ --in-ont ~/appa/data/nanopore/subreads/circ1.fastq --threads 32


./dorado correct ~/appa/data/nanopore/humanGut/hac_duplex/bgzip/DON_020D10_1.fastq.gz -t 12 -m .temp_dorado_model-a2c679592e85eee7/herro-v1/ > ~/appa/data/nanopore/humanGut/hac_duplex/bgzip/DON_020D10_1_herro.fasta



./bin/metaMDBG asm --out-dir ~/appa/run/correction/simulation/simulation_1/ --in-ont ~/appa/data/nanopore/simulation/simulation_1/reads.fastq --threads 32
./bin/metaMDBG setupCorr -o ~/appa/run/correction/simulation/simulation_1/tmp/ --ref ~/appa/data/nanopore/simulation/simulation_1/reference_input.txt --metadata ~/appa/data/nanopore/simulation/simulation_1/reads_metadata.tsv
./bin/metaMDBG readCorrection ~/appa/run/correction/simulation/simulation_1/tmp/ --ref ~/appa/data/nanopore/simulation/simulation_1/reference_input.txt --metadata ~/appa/data/nanopore/simulation/simulation_1/reads_metadata.tsv -t 32



./bin/metaMDBG setupCorr --out-dir ~/appa/run/correction/eval/tmp/ --threads 32
*/




class ReadCorrection : public Tool{
    
public:



	struct MinimizerPairPosition{

		ReadType _readIndex;
		u_int16_t _positionIndex;
	};

	typedef phmap::parallel_flat_hash_map<KmerVec, vector<MinimizerPairPosition>, phmap::priv::hash_default_hash<KmerVec>, phmap::priv::hash_default_eq<KmerVec>, std::allocator<std::pair<KmerVec, vector<MinimizerPairPosition>>>, 4, std::mutex> KminmerPosMap;

	struct MinimizerPos{
		u_int32_t _readIndex;
		u_int32_t _position1;
		u_int32_t _position2;
		bool _isReversed;
	};


	//typedef phmap::parallel_flat_hash_map<u_int32_t, u_int32_t, phmap::priv::hash_default_hash<u_int32_t>, phmap::priv::hash_default_eq<u_int32_t>, std::allocator<std::pair<u_int32_t, u_int32_t>>, 4, std::mutex> ReadLengthMap;

	typedef phmap::parallel_flat_hash_map<KmerVec, vector<MinimizerPos>, phmap::priv::hash_default_hash<KmerVec>, phmap::priv::hash_default_eq<KmerVec>, std::allocator<std::pair<KmerVec, vector<MinimizerPos>>>, 4, std::mutex> KminmerIndex;
	//typedef phmap::parallel_flat_hash_map<KmerVec, vector<u_int32_t>, phmap::priv::hash_default_hash<KmerVec>, phmap::priv::hash_default_eq<KmerVec>, std::allocator<std::pair<KmerVec, vector<u_int32_t>>>, 4, std::mutex> KminmerReadMap;
	//typedef phmap::parallel_flat_hash_map<KmerVec, u_int32_t, phmap::priv::hash_default_hash<KmerVec>, phmap::priv::hash_default_eq<KmerVec>, std::allocator<std::pair<KmerVec, u_int32_t>>, 4, std::mutex> KminmerCountMap;

	//typedef phmap::parallel_flat_hash_map<u_int64_t, vector<u_int64_t>, phmap::priv::hash_default_hash<u_int64_t>, phmap::priv::hash_default_eq<u_int64_t>, std::allocator<std::pair<u_int64_t, vector<u_int64_t>>>, 4, std::mutex> ReadToReadsMap;
	//typedef phmap::parallel_flat_hash_map<u_int64_t, vector<MinimizerType>, phmap::priv::hash_default_hash<u_int64_t>, phmap::priv::hash_default_eq<u_int64_t>, std::allocator<std::pair<u_int64_t, vector<MinimizerType>>>, 4, std::mutex> ReadToMinimizersMap;

	string _inputFilename;
	string _inputDir;
	string _outputFilename;
	float _minimizerDensity_assembly;
    size_t _minimizerSize;
    size_t _kminmerSize;		
	float _minimizerSpacingMean;
	float _kminmerLengthMean;
	float _kminmerOverlapMean;
	size_t _kminmerSizeFirst;
	size_t _kminmerSizePrev;
	size_t _kminmerSizeLast;
	size_t _meanReadLength;
	float _minimizerDensity_correction;
	float _minReadQuality;
	float _minIdentity;
	float _minOverlapLength;

	int _nbCores;
	string _filename_exe;
	

	u_int64_t _minReadLength;
	size_t _totalNbReads;
	u_int64_t _usedCoverageLowDensity;
	u_int64_t _usedCoverageForCorrection;

	u_int64_t _maxUsedReadOverlapForCorrection;

	string _alignmentFilename;
	ofstream _alignmentFile;
	string _partitionFilename;
	//ReadLengthMap _readIndex_to_readSize;



	vector<vector<MinimizerType>> _mReads_evaluation;

	struct MinimizerPosition{

		//ReadType _readIndex;
		int32_t _position;
		int32_t _positionIndex;
		bool _isReversed;
	};


	struct SimulatedReadMetadata{

		public:

		ReadType _referenceIndex;
		ReadType _readIndex;
		bool _isReversed;
		u_int32_t _refCoordStart;
		u_int32_t _refCoordEnd;
		float _identity;

		string toString() const{
			return to_string(_referenceIndex) + "\t" + to_string(_readIndex) + "\t" + to_string(_isReversed) + "\t" + to_string(_refCoordStart) + "\t" + to_string(_refCoordEnd) + "\t" + to_string(_identity) ;
		}

		int32_t distanceFrom(const SimulatedReadMetadata& other) const{
			return abs((int64_t)_refCoordStart - (int64_t)other._refCoordStart);
		}
	};

	vector<SimulatedReadMetadata> _simulatedReadMetadata;
	typedef unordered_map<MinimizerType, vector<MinimizerPosition>> ReadMinimizerPositionMap;

	typedef phmap::parallel_flat_hash_map<KmerVec, vector<ReadType>, phmap::priv::hash_default_hash<KmerVec>, phmap::priv::hash_default_eq<KmerVec>, std::allocator<std::pair<KmerVec, vector<ReadType>>>, 4, std::mutex> KminmerReadMap;
	typedef phmap::parallel_flat_hash_map<MinimizerType, vector<MinimizerPosition>, phmap::priv::hash_default_hash<MinimizerType>, phmap::priv::hash_default_eq<MinimizerType>, std::allocator<std::pair<MinimizerType, vector<MinimizerPosition>>>, 4, std::mutex> MinimizerReadMap;
	typedef phmap::parallel_flat_hash_map<ReadType, vector<ReadType>, phmap::priv::hash_default_hash<ReadType>, phmap::priv::hash_default_eq<ReadType>, std::allocator<std::pair<ReadType, vector<ReadType>>>, 4, std::mutex> ReadIndexMap;


    struct ReadWriter{
        ReadType _readIndex;
        vector<MinimizerType> _minimizers;
        vector<u_int8_t> _qualities;
    };

    struct ReadWriter_Comparator {
        bool operator()(ReadWriter const& p1, ReadWriter const& p2){
            return p1._readIndex > p2._readIndex;
        }
    };

	priority_queue<ReadWriter, vector<ReadWriter> , ReadWriter_Comparator> _readWriterQueue;
	ReadType _nextReadIndexWriter;
	ofstream _file_readData;

    struct KminmerAbundance{
        KmerVec _kminmer;
        u_int32_t _abundance;
    };

    struct KminmerAbundanceComparator {
        bool operator()(KminmerAbundance const& p1, KminmerAbundance const& p2){
            return p1._abundance < p2._abundance;
        }
    };


	//string _illuminaToNanoporeMappingFilename;
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

		//vector<AlignmentResult> _alignmentResults;

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
		/*
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
		*/
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
		//unordered_set<MinimizerType> _referenceMinimizers;
		unordered_map<u_int32_t, Node*> _nodes;
		vector<unordered_map<MinimizerType, vector<ReadType>>> _referencePosition_to_minimizerCounts;

		Graph(const vector<MinimizerType>& minimizers, const vector<u_int8_t>& qualities){

			_referencePosition_to_minimizerCounts.resize(minimizers.size());

			_print = false;

			for(size_t i=0; i<minimizers.size(); i++){
				//_referenceMinimizers.insert(minimizers[i]);
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

		/*
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
		*/


		void addAlignment2(const AlignmentResult2& alignment, const vector<MinimizerType>& referenceMinimizers, const vector<MinimizerType>& queryMinimizers, const vector<u_int8_t>& queryQualities){

			u_int32_t firstMatchPosition = -1;
			u_int32_t endMatchPosition = -1;


			for(const auto& it : alignment._alignments){

				u_int32_t referencePosition = it.first;
				u_int32_t queryPosition = it.second;

				if(referencePosition == -1){ //insertion
				}
				else if(queryPosition == -1){ //deletion

				}
				else if(referenceMinimizers[referencePosition] == queryMinimizers[queryPosition]){ //Match
					if(firstMatchPosition == -1) firstMatchPosition = referencePosition;
					endMatchPosition = referencePosition;
				}
				else { //missmatch
				}

			}
			

			Node* prevNode = nullptr;
			Node* currentNode = nullptr;
			//if(readMinimizers.size() < 2) return;

			for(const auto& it : alignment._alignments){

				u_int32_t referencePosition = it.first;
				u_int32_t queryPosition = it.second;

				//cout << endl;
				//cout << referencePosition << " " << referenceMinimizers.size() << endl;
				//cout << queryPosition << " " << queryMinimizers.size() << endl;
				//cout << endl;

				if(referencePosition == -1){ //insertion
					currentNode = addNode(prevNode, queryMinimizers[queryPosition]);
					currentNode->_abundance += 1;
					currentNode->addQuality(queryQualities[queryPosition]);

					if(prevNode != nullptr){
						addEdge(prevNode, currentNode, queryQualities[queryPosition]);
					}

					//cout << referencePosition << "\t" << referenceReadMinimizers[referencePosition] << endl;
					//cout << queryPosition << "\t" << queryReadMinimizers[queryPosition] << endl;


					prevNode = currentNode;
				}
				else if(queryPosition == -1){ //deletion

				}
				else if(referenceMinimizers[referencePosition] == queryMinimizers[queryPosition]){ //Match
					currentNode = _nodes[referencePosition];
					currentNode->_abundance += 1;
					currentNode->addQuality(queryQualities[queryPosition]);

					//cout << (currentNode != nullptr) << " " << (prevNode!= nullptr) << endl;

					if(prevNode != nullptr){
						addEdge(prevNode, currentNode, queryQualities[queryPosition]);
					}

					prevNode = currentNode;
				}
				else { //missmatch

					//if(referencePosition < firstMatchPosition) continue;
					//if(referencePosition > endMatchPosition) continue;

					//if(queryPosition > startMatchPos && queryPosition < endMatchPos){
					currentNode = addNode(prevNode, queryMinimizers[queryPosition]);
					currentNode->_abundance += 1;
					currentNode->addQuality(queryQualities[queryPosition]);

					if(prevNode != nullptr){
						addEdge(prevNode, currentNode, queryQualities[queryPosition]);
					}
					//}

					prevNode = currentNode;
				}

			}
			

		}

		
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
				/*
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
				*/

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
		u_int64_t computeReadAbundance(u_int64_t nbErroneousMinimizers, const vector<MinimizerType>& referenceMinimizers){

			vector<u_int32_t> minimizerAbundances;


			for(size_t i=0; i<_referencePosition_to_minimizerCounts.size(); i++){

				MinimizerType referenceMinimizer = referenceMinimizers[i];


				for(const auto& it2 : _referencePosition_to_minimizerCounts[i]){

					MinimizerType minimizer = it2.first;
					const vector<ReadType>& readIndexes = it2.second;
					u_int32_t count = readIndexes.size();

					if(minimizer == referenceMinimizer){
						minimizerAbundances.push_back(count);
					}
				}
			}

			//while(minimizerAbundances.size() < referenceMinimizers.size()){
			//	minimizerAbundances.push_back(1);
			//}

			std::sort(minimizerAbundances.begin(), minimizerAbundances.end());

			cout << "Read size: " << referenceMinimizers.size() << endl;
			//int error = 0;
			for(u_int32_t ab : minimizerAbundances){
				cout << ab << endl;
			}

			//getchar();

			if(minimizerAbundances.size() == 0) return 0;
			if(nbErroneousMinimizers >= minimizerAbundances.size()) return minimizerAbundances[minimizerAbundances.size()-1];

			return minimizerAbundances[nbErroneousMinimizers];
		}

		bool isLowCovered(const vector<MinimizerType>& referenceMinimizers){

			vector<u_int32_t> minimizerAbundances;


			for(size_t i=0; i<_referencePosition_to_minimizerCounts.size(); i++){

				MinimizerType referenceMinimizer = referenceMinimizers[i];


				for(const auto& it2 : _referencePosition_to_minimizerCounts[i]){

					MinimizerType minimizer = it2.first;
					const vector<ReadType>& readIndexes = it2.second;
					u_int32_t count = readIndexes.size();

					if(minimizer == referenceMinimizer){
						minimizerAbundances.push_back(count);
					}
				}
			}


			while(minimizerAbundances.size() < referenceMinimizers.size()){
				minimizerAbundances.push_back(1);
			}

			u_int64_t nbLowCoveragedMinimizers = 0;

			for(u_int32_t ab : minimizerAbundances){
				//cout << ab << endl;
				if(ab < 5) nbLowCoveragedMinimizers += 1;
			}

			//cout << nbLowCoveragedMinimizers << " " << referenceMinimizers.size() << endl;
			//getchar();
			return nbLowCoveragedMinimizers > (referenceMinimizers.size() / 2);
		}

		void clearPileup(){
			for(size_t i=0; i<_referencePosition_to_minimizerCounts.size(); i++){
				_referencePosition_to_minimizerCounts[i].clear();
			}
		}

		void pileupAlignment(const spoa64::Alignment& alignment, const vector<MinimizerType>& referenceMinimizers, u_int64_t readIndex, const vector<MinimizerType>& readMinimizers, const vector<MinimizerType>& readQualities,const AlignmentResult alignmentResult){


			for (const auto& it : alignment) {
				int64_t v1 = it.first;
				int64_t v2 = it.second;
				
				if(v1 == -1){ //insert in

				}
				else if(v2 == -1){ //insert in
				}
				else if(referenceMinimizers[v1] == readMinimizers[v2]){ //match

					_referencePosition_to_minimizerCounts[v1][readMinimizers[v2]].push_back(readIndex);
					
				}
				else{ //missmatch


					_referencePosition_to_minimizerCounts[v1][readMinimizers[v2]].push_back(readIndex);
				}

			}

		}

		void pileupAlignment2(const AlignmentResult2& alignment, const vector<MinimizerType>& referenceMinimizers, const vector<MinimizerType>& queryMinimizers){


			vector<Utils::CigarElement> cigarSequence = Utils::extractCigarSequence(alignment._cigar);

			size_t referencePosition = alignment._cigarReferenceStart;
			size_t queryPosition = alignment._cigarQueryStart;


			for(const Utils::CigarElement& cigarElement : cigarSequence){

				if(cigarElement._cigarType == Utils::CigarType::Match){

					for(size_t i=0; i<cigarElement._nbOccurences; i++){

						if(referenceMinimizers[referencePosition] == queryMinimizers[queryPosition]){ //Match
							_referencePosition_to_minimizerCounts[referencePosition][queryMinimizers[queryPosition]].push_back(alignment._queryReadIndex);
						}
						else{ //missmatch

							_referencePosition_to_minimizerCounts[referencePosition][queryMinimizers[queryPosition]].push_back(alignment._queryReadIndex);


						}

						referencePosition += 1;
						if(alignment._isQueryReversed){
							queryPosition -= 1;
						}	
						else{
							queryPosition += 1;
						}	
					}
				}
				else if(cigarElement._cigarType == Utils::CigarType::Insertion){

					for(size_t i=0; i<cigarElement._nbOccurences; i++){


						if(alignment._isQueryReversed){
							queryPosition -= 1;
						}	
						else{
							queryPosition += 1;
						}	
					}
					
				}
				else if(cigarElement._cigarType == Utils::CigarType::Deletion){

					for(size_t i=0; i<cigarElement._nbOccurences; i++){

						referencePosition += 1;
					}


				}

			}

		}

		unordered_set<ReadType> detectStrainReads(const vector<MinimizerType>& referenceMinimizers){

			bool print = false;

			unordered_set<ReadType> allStrainReadIndexes;
			if(print) cout << endl << endl;

			for(size_t i=0; i<_referencePosition_to_minimizerCounts.size(); i++){

				if(print) cout << "\tColumn: " << i << endl;
				MinimizerType referenceMinimizer = referenceMinimizers[i];

				bool isSolidPosition = false;
				vector<ReadType> strainReadIndexes; 

				for(const auto& it2 : _referencePosition_to_minimizerCounts[i]){

					MinimizerType minimizer = it2.first;
					const vector<ReadType>& readIndexes = it2.second;
					u_int32_t count = readIndexes.size();

					if(minimizer == referenceMinimizer){
						if(count >= 5){
							isSolidPosition = true;
						}
						if(print) cout << "\t\t" << minimizer << " " << count << " ****" << endl;
					}
					else{
						if(count >= 5){
							for(ReadType readIndex : readIndexes){
								strainReadIndexes.push_back(readIndex);
							}
							if(print) cout << "\t\t" << minimizer << " " << count << " |||" << endl;
						}
						else{
							if(print) cout << "\t\t" << minimizer << " " << count << endl;
						}
					}



				}

				if(isSolidPosition){
					for(ReadType readIndex : strainReadIndexes){
						allStrainReadIndexes.insert(readIndex);
					}
				}

			}

			//vector<u_int64_t> strainReadIndexesVec;
			//for(u_int64_t readIndex : allStrainReadIndexes){
			//	strainReadIndexesVec.push_back(readIndex);
			//}

			return allStrainReadIndexes;
		}
		*/
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

	string _simulatedReadMetadataFilename;

	ReadCorrection(): Tool (){
	}


	void parseArgs(int argc, char* argv[]){

		_filename_exe = argv[0];

		args::ArgumentParser parser("readSelection", ""); //"This is a test program.", "This goes after the options."
		args::Positional<std::string> arg_outputDir(parser, "outputDir", "", args::Options::Required);
		//args::ValueFlag<string> arg_pafFilename(parser, "", "paf filename", {"paf"}, "");
		//args::ValueFlag<string> arg_hifiFilename(parser, "", "hifi filename", {"hifi"}, "");
		//args::ValueFlag<string> arg_illuminaFilename(parser, "", "illumina filename", {"illumina"}, "");
		args::ValueFlag<std::string> arg_referenceFilename(parser, "", "referenceFilename", {"ref"});
		args::ValueFlag<std::string> arg_metadataFilename(parser, "", "metadataFilename", {"metadata"});
		//args::Positional<std::string> arg_outputFilename(parser, "outputFilename", "", args::Options::Required);
		//args::Positional<std::string> arg_inputReadFilename(parser, "inputReadFilename", "", args::Options::Required);
		//args::Positional<std::string> arg_contigs(parser, "contigs", "", args::Options::Required);
		//args::PositionalList<std::string> arg_readFilenames(parser, "reads", "Input filename(s) (separated by space)", args::Options::Required);
		//args::ValueFlag<int> arg_l(parser, "", "Minimizer length", {ARG_MINIMIZER_LENGTH2}, 13);
		//args::ValueFlag<float> arg_d(parser, "", "Minimizer density", {ARG_MINIMIZER_DENSITY2}, 0.005f);
		//args::ValueFlag<std::string> arg_contigs(parser, "", "", {ARG_INPUT_FILENAME_CONTIG}, "");
		args::ValueFlag<int> arg_minOverlapLength(parser, "", "Min overlap length", {"min-overlap-length"}, 1000);
		args::ValueFlag<float> arg_minIdentity(parser, "", "Min identity", {"min-identity"}, 0.96f);
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

		/*
		u_int32_t maxHashValue_32bit = -1;
		u_int64_t maxHashValue_64bit = -1;

		float density = 0.005;
		u_int32_t minimizerBound_32bit = density * maxHashValue_32bit;
		u_int64_t minimizerBound_64bit = density * maxHashValue_64bit;

		cout << std::fixed << minimizerBound_32bit << endl;
		cout << std::fixed  << minimizerBound_64bit << endl;
		exit(1);
		*/
		_inputDir = args::get(arg_outputDir);
		_simulatedReadMetadataFilename = args::get(arg_metadataFilename);
		//_outputFilename = args::get(arg_outputFilename);
		//_inputFilename = args::get(arg_inputReadFilename);
		_nbCores = args::get(arg_nbCores);
		_minIdentity = args::get(arg_minIdentity);
		_minOverlapLength = args::get(arg_minOverlapLength);

		//_illuminaFilename = "/pasteur/appa/homes/gbenoit/appa/data/nanopore/AD/input_paired.txt";
		//_illuminaFilename = "/pasteur/appa/homes/gbenoit/appa/data/nanopore/subreads/circ1_illumina.txt";
		
		string filename_parameters = _inputDir + "/parameters.gz";
		gzFile file_parameters = gzopen(filename_parameters.c_str(),"rb");	
		gzread(file_parameters, (char*)&_minimizerSize, sizeof(_minimizerSize));
		gzread(file_parameters, (char*)&_kminmerSize, sizeof(_kminmerSize));
		gzread(file_parameters, (char*)&_minimizerDensity_assembly, sizeof(_minimizerDensity_assembly));
		gzread(file_parameters, (char*)&_kminmerSizeFirst, sizeof(_kminmerSizeFirst));
		gzread(file_parameters, (char*)&_minimizerSpacingMean, sizeof(_minimizerSpacingMean));
		gzread(file_parameters, (char*)&_kminmerLengthMean, sizeof(_kminmerLengthMean));
		gzread(file_parameters, (char*)&_kminmerOverlapMean, sizeof(_kminmerOverlapMean));
		gzread(file_parameters, (char*)&_kminmerSizePrev, sizeof(_kminmerSizePrev));
		gzread(file_parameters, (char*)&_kminmerSizeLast, sizeof(_kminmerSizeLast));
		gzread(file_parameters, (char*)&_meanReadLength, sizeof(_meanReadLength));
		gzread(file_parameters, (char*)&_minimizerDensity_correction, sizeof(_minimizerDensity_correction));

		//gzread(file_parameters, (char*)&_minimizerSize, sizeof(_minimizerSize));
		//gzread(file_parameters, (char*)&_kminmerSize, sizeof(_kminmerSize));
		//gzread(file_parameters, (char*)&_minimizerDensity_assembly, sizeof(_minimizerDensity_assembly));
		gzclose(file_parameters);

		openLogFile(_inputDir);

		Logger::get().debug() << "";
		Logger::get().debug() << "Input filename: " << _inputFilename;
		Logger::get().debug() << "Input dir: " << _inputDir;
		Logger::get().debug() << "Minimizer length: " << _minimizerSize;
		Logger::get().debug() << "Kminmer length: " << _kminmerSize;
		Logger::get().debug() << "Assembly Density: " << _minimizerDensity_assembly;
		Logger::get().debug() << "Correction Density: " << _minimizerDensity_correction;
		Logger::get().debug() << "Min overlap length: " << _minOverlapLength;
		Logger::get().debug() << "Min identity: " << _minIdentity;
		Logger::get().debug() << "";

		_kminmerSize = 2;
		_minReadQuality = 10;
		_nbIgnoredMostAbundantMinimizers = 10000;
		_usedCoverageLowDensity = 20;
		_usedCoverageForCorrection = 20;
		_minReadLength = 10; //2000*_minimizerDensity_assembly;
		//_minNbMinimizersReference = 10;
		//_minNbMinimizersQuery = 10;
		_maxUsedReadOverlapForCorrection = 100;
		_correctionCheckSum = 0;
		_alignmentCheckSum = 0;
		_eval_nbBases = 0;
		_eval_nbMatches = 0;
		_fractionRepititiveMinimizersToFilterOut = 0.0002;
		_maxChainingBand_lowDensity = (u_int64_t) 2500 * _minimizerDensity_correction;
		_maxChainingBand_highDensity = (u_int64_t) 2500 * _minimizerDensity_correction;
		_alignmentFilename = _inputDir + "/readAlignmentsLowDensity.bin";
		_partitionFilename = _inputDir + "/readPartition.bin";

		_fields = new vector<string>();
	}

	u_int64_t _eval_nbBases;
	u_int64_t _eval_nbMatches;
	//u_int64_t _nbReads;
	u_int64_t _nbIgnoredMostAbundantMinimizers;
	float _fractionRepititiveMinimizersToFilterOut;
	u_int64_t _maxChainingBand_lowDensity;
	u_int64_t _maxChainingBand_highDensity;

    void execute (){

		ifstream file_readStats(_inputDir + "/read_stats.txt");

		float minimizerDensity;
		u_int32_t n50ReadLength;
		file_readStats.read((char*)&_totalNbReads, sizeof(_totalNbReads));
		file_readStats.read((char*)&n50ReadLength, sizeof(n50ReadLength));
		file_readStats.read((char*)&minimizerDensity, sizeof(minimizerDensity));

		file_readStats.close();

		//cout << "Nb reads: " << _totalNbReads << endl;
		//cout << "Minimizer density: " << minimizerDensity << endl;
		//cout << "Average distance between minimizers: " << 1 / minimizerDensity << " bps" << endl;
		//cout << "Assembly density: " << _minimizerDensity_assembly << endl;


		
		//cout << "a metter dans une fonction pour bien clear la mémoire" << endl;
		//{
			
		//	size_t minNbAnchors = 3; //max((size_t)2, (size_t)(_minMatchLength*_minimizerDensity_correction));
			//ReadMapper readMapper(_inputDir + "/read_data_init.txt", _inputDir + "/read_data_init.txt", _inputDir, _minimizerSize, _minimizerDensity_assembly, _minimizerDensity_correction, _totalNbReads, _nbReadsPerPartition, _minReadLength, _maxMappedReads, minNbAnchors, 0.04, _nbCores);
			//readMapper.execute();
			//exit(1);
		//}
		_nextReadIndexWriter = 0;
		_nextAlignmentReadIndexWriter = 0;


		_print_debug = false;
		_eval_correction = false;

		if(_eval_correction){

			//loadSimulatedReadMetadataFile();

			//cout << "Loading minimizer reads evaluation" << endl;
			//loadMinimizerReadsEvaluation(_inputDir + "/read_data_simulated.txt");

			//cout << "Loading read to contig mapping" << endl;
			//loadReadToContigMapping();

			//cout << "Loading contigs" << endl;
			//loadContigs();



			//cout << "Loading read to read mapping" << endl;
			//loadReadToReadMapping();


		}

		//cout << "Loading read to read mapping" << endl;
		//loadReadToReadMapping();
		//cout << "Loading minimizer reads" << endl;
		//loadMinimizerReads(_inputDir + "/read_data_init.txt");
		evalMapping();
		//return;


		
		Logger::get().debug() << "\nIndexing genomic minimizers";
		indexGenomicMinimizers();
		
		Logger::get().debug() << "\nIndexing reads";
		indexReads();

		
		Logger::get().debug() << "\nLoading minimizer reads";
		loadMinimizerReads(_inputDir + "/read_data_init.txt", true, false);


		//_file_readData = ofstream(_inputDir + "/read_data_corrected.txt");
		


		Logger::get().debug() << "\nPerform low density alignments";
		alignReadsLowDensity();
		
		//_file_readData.close();

		for(size_t i=0; i<_mReads.size(); i++){
			_mReads[i] = {};
		}

		//_mReads.clear();
		_kminmer_to_readIndex.clear();

		Logger::get().debug() << "\nLoading low density alignments";
		loadAlignments(false);

		Logger::get().debug() << "\nPartitioning reads";
		int nbPartitions = partitionReads();


		ifstream partitionFile(_partitionFilename);
		_file_readData = ofstream(_inputDir + "/read_data_corrected.txt");


		_isReadToLoad.resize(_totalNbReads, false);
		_isReadToCorrect.resize(_totalNbReads, false);

		int currentPartition = 0;

		while(true){

			Logger::get().debug() << "\nProcessing partition: " << currentPartition << " / " << nbPartitions;

			for(size_t i=0; i<_isReadToLoad.size(); i++){
				_isReadToLoad[i] = false;
				_isReadToCorrect[i] = false;
			}

			u_int64_t readToLoadSize;
			u_int64_t readToCorrectSize;
			ReadPartition partition;


			partitionFile.read((char*)&readToLoadSize, sizeof(readToLoadSize));

			if(partitionFile.eof()) break;

			partition._readsToLoad.resize(readToLoadSize);
			partitionFile.read((char*)&partition._readsToLoad[0], partition._readsToLoad.size() * sizeof(ReadType));

			partitionFile.read((char*)&readToCorrectSize, sizeof(readToCorrectSize));
			partition._readsToCorrect.resize(readToCorrectSize);
			partitionFile.read((char*)&partition._readsToCorrect[0], partition._readsToCorrect.size() * sizeof(ReadType));

			for(const ReadType& readIndex : partition._readsToLoad){
				_isReadToLoad[readIndex] = true;
			}

			for(const ReadType& readIndex : partition._readsToCorrect){
				_isReadToCorrect[readIndex] = true;
				
			}

			//cout << partition._readsToLoad.size() << " " << partition._readsToCorrect.size() << endl;
			//getchar();

			loadAlignments(true);

			Logger::get().debug() << "\tLoading minimizer reads (high density)";
			loadMinimizerReads(_inputDir + "/read_data_init.txt", false, true);

			Logger::get().debug() << "\tCorrecting reads";
			correctReads();

			Logger::get().debug() << "\tCorrection checksum: " << _correctionCheckSum;
			currentPartition += 1;
        }

		_file_readData.close();

		partitionFile.close();
		fs::remove(_inputDir + "/readAlignmentsLowDensity.bin");
		fs::remove(_inputDir + "/readPartition.bin");
		

		//unordered_set<ReadType> correctedReadsReal;
		//unordered_set<ReadType> correctedReads;

		/*
		for(size_t i=0; i<_readPartitions.size(); i++){

			const ReadPartition& partition = _readPartitions[i];

			_mReads.clear();
			_isReadToLoad.resize(_lowDensityAlignments.size(), false);
			_isReadToCorrect.resize(_lowDensityAlignments.size(), false);

			for(size_t j=0; j<partition._readsToLoad.size(); j++){
				_isReadToLoad[partition._readsToLoad[j]] = true;
			}
			for(size_t j=0; j<partition._readsToCorrect.size(); j++){
				//correctedReads.insert(partition._readsToCorrect[j]);
				_isReadToCorrect[partition._readsToCorrect[j]] = true;
			}


			for(size_t j=0; j<_isReadToCorrect.size(); j++){
				//_isReadToCorrect[j] = true;
			}

			for(size_t j=0; j<_isReadToLoad.size(); j++){
				//_isReadToLoad[j] = true;
				//correctedReadsReal.insert(j);
			}

			cout << "Loading minimizer reads (high density)" << endl;
			loadMinimizerReads(_inputDir + "/read_data_init.txt", false, true);

			cout << "Correcting mReads" << endl;
			correctReads();

			//break;
		}
		*/

		//cout << correctedReadsReal.size() << " " << correctedReads.size() << endl;
		//getchar();


		Logger::get().debug() << "\tAlignment check sum: " << _alignmentCheckSum;
		Logger::get().debug() << "\tCorrection check sum: " << _correctionCheckSum;

		//Experiment nb pairs vs mapping
		//_debug_outputFile_nbMinimizerPairs.close();
		
		//closeLogFile();
	}

	vector<bool> _isReadToLoad;
	vector<bool> _isReadToCorrect;
	unordered_set<KmerVec> _isRepetitiveKminmers;


	BloomCacheCoherent<u_int64_t>* _bloomFilter;
	//BloomCacheCoherent<u_int64_t>* _bloomGenomic;
	
	void indexGenomicMinimizers(){

		//cout << "Bloom filter en mode determinist! " << endl;
		_bloomFilter = new BloomCacheCoherent<u_int64_t>(64000000000ull);
		//_bloomGenomic = new BloomCacheCoherent<u_int64_t>(16000000000ull);

		KminmerParserParallel parser(_inputDir + "/read_data_init.txt", _minimizerSize, _kminmerSize, false, true, _nbCores);
		parser._densityThreshold = _minimizerDensity_assembly;
		parser.parse(IndexGenomicMinimizersFunctor(*this));

		delete _bloomFilter;
	}

	class IndexGenomicMinimizersFunctor {

		public:

		ReadCorrection& _parent;

		IndexGenomicMinimizersFunctor(ReadCorrection& parent) : _parent(parent){
		}

		IndexGenomicMinimizersFunctor(const IndexGenomicMinimizersFunctor& copy) : _parent(copy._parent){
		}

		~IndexGenomicMinimizersFunctor(){
		}

		void operator () (const KminmerList& kminmerList) {

			//if(kminmerList._minimizerPos.size() == 0 || kminmerList._minimizerPos[kminmerList._minimizerPos.size()-1] < _parent._minReadLength) return;
			
			//if(!Utils::isValidPositions(kminmerList._minimizerPos)) return;

			ReadType readIndex = kminmerList._readIndex;
			if(readIndex % 100000 == 0)  Logger::get().debug() << "\tIndexing genomic minimizers: " << Utils::getProgress(readIndex, _parent._totalNbReads);


			MinimizerRead read = {kminmerList._readIndex, kminmerList._readMinimizers, kminmerList._minimizerPos, kminmerList._readQualities, kminmerList._readMinimizerDirections, kminmerList._readLength};
			MinimizerRead readLowDensity = Utils::getLowDensityMinimizerRead(read, _parent._minimizerDensity_assembly);
			
			if(_parent.isReadTooShort(readLowDensity)) return;
			
			//if(kminmerList._readMinimizers.size() < _parent._minReadLength) return;
			//if(kminmerList._readLength < _parent._minOverlapLength) return;

			for(u_int32_t i=0; i<kminmerList._kminmersInfo.size(); i++){
			
				
				//cout << "------- " << i << endl;

				//_logFile << readIndex << " " << i << endl;
				const ReadKminmerComplete& kminmerInfo = kminmerList._kminmersInfo[i];

				const KmerVec& vec = kminmerInfo._vec;
			
				bool exist = false;


				//#pragma omp critical(IndexGenomicMinimizersFunctor)
				//{
					
					if(_parent._bloomFilter->contains(vec.h())){
						exist = true;
					}
					else{
						exist = false;
						_parent._bloomFilter->insert(vec.h());
					}
				//}


				if(!exist) continue;

				//_parent._bloomGenomic->insert(vec.h());
				
				_parent._kminmer_to_readIndex.lazy_emplace_l(vec,
				[&readIndex](KminmerPosMap::value_type& v) { // key exist
					//v.second[0] += 1; //Increment kminmer count
					//if(std::find(v.second.begin(), v.second.end(), readIndex) == v.second.end()){
					//	v.second.push_back(readIndex);
					//}
				},           
				[&vec, &readIndex](const KminmerPosMap::constructor& ctor) { // key inserted
					
					//vector<u_int32_t> readIndexes = {2}; //inital count of this kminmer
					vector<MinimizerPairPosition> readIndexes = {}; //inital count of this kminmer

					ctor(vec, readIndexes); 

				}); // construct value_type in place when key not present
				



			}


		}
	};
	
	vector<string>* _fields;



	//KminmerIndex _kminmer_to_readIndex;
	KminmerPosMap _kminmer_to_readIndex;

	vector<MinimizerRead> _mReadsDebug;
	vector<MinimizerRead> _mReads;
	//vector<MinimizerRead> _mReadsHiFi;
	//vector<MinimizerReadBinary> _mReadsBinary;
	MinimizerReadMap _minimizer_to_readIndex;

	void loadMinimizerReads(const string& filename, bool isLowDensity, bool useReadPartionning){

		if(_mReads.size() == 0){
			_mReads.resize(_totalNbReads, {});
		}
		//_mReads.clear();

		KminmerParserParallel parser(filename, _minimizerSize, _kminmerSize, false, true, _nbCores);

		if(isLowDensity){
			parser._densityThreshold = _minimizerDensity_assembly;
		}

		parser.parseSequences(LoadMinimizerReadsFunctor(*this, isLowDensity, useReadPartionning));
	}

	class LoadMinimizerReadsFunctor {

		public:

		ReadCorrection& _parent;
		bool _isLowDensity;
		bool _useReadPartionning;

		LoadMinimizerReadsFunctor(ReadCorrection& parent, bool isLowDensity, bool useReadPartionning) : _parent(parent){
			_isLowDensity = isLowDensity;
			_useReadPartionning = useReadPartionning;
		}

		LoadMinimizerReadsFunctor(const LoadMinimizerReadsFunctor& copy) : _parent(copy._parent){
			_isLowDensity = copy._isLowDensity;
			_useReadPartionning = copy._useReadPartionning;
		}

		~LoadMinimizerReadsFunctor(){
		}


		void operator () (const KminmerList& kminmerList) {

			ReadType readIndex = kminmerList._readIndex;			
			if(readIndex % 100000 == 0 && _isLowDensity)  Logger::get().debug() << "\tLoading minimizer reads: " << Utils::getProgress(readIndex, _parent._totalNbReads);


			MinimizerRead read = {kminmerList._readIndex, kminmerList._readMinimizers, kminmerList._minimizerPos, kminmerList._readQualities, kminmerList._readMinimizerDirections, kminmerList._readLength};
			MinimizerRead readLowDensity = Utils::getLowDensityMinimizerRead(read, _parent._minimizerDensity_assembly);
			
			if(_isLowDensity){
				read._qualities.clear();
				readLowDensity._qualities.clear();
			}

			if(_parent.isReadTooShort(readLowDensity)){
				//_parent._mReads[readIndex] 
				//_parent._mReads.push_back({});
			}
			else if(_useReadPartionning && !_parent._isReadToLoad[readIndex]){
				_parent._mReads[readIndex] = {};
				//_parent._mReads.push_back({});
			}
			else{
				if(_parent._mReads[readIndex]._minimizers.size() > 0) return;
				_parent._mReads[readIndex] = read;
				//_parent._mReads.push_back(read);
				//_parent._mReads.push_back(read);
			}

			

		}
	};

	void loadMinimizerReadsEvaluation(const string& filename){
		KminmerParserParallel parser(filename, _minimizerSize, _kminmerSize, false, false, 1);
		parser.parseSequences(LoadMinimizerReadsEvaluationFunctor(*this));
	}

	class LoadMinimizerReadsEvaluationFunctor {

		public:

		ReadCorrection& _parent;

		LoadMinimizerReadsEvaluationFunctor(ReadCorrection& parent) : _parent(parent){
		}

		LoadMinimizerReadsEvaluationFunctor(const LoadMinimizerReadsEvaluationFunctor& copy) : _parent(copy._parent){
		}

		~LoadMinimizerReadsEvaluationFunctor(){
		}


		void operator () (const KminmerList& kminmerList) {

			ReadType readIndex = kminmerList._readIndex;
			if(readIndex % 100000 == 0) cout << "Loading minimizer reads: " << readIndex << " / " << _parent._totalNbReads << endl;

			_parent._mReads_evaluation.push_back(kminmerList._readMinimizers);

		}
	};

	bool isReadTooShort(const MinimizerRead& minimizerReadLowDensity){
		if(minimizerReadLowDensity._minimizers.size() < 10) return true;
		//if(minimizerReadLowDensity._readLength < 500) return true;
		//if(minimizerReadLowDensity._readLength < _minOverlapLength) return true;

		return false;
	}

	vector<vector<MinimizerPosition>> _minimizer_to_readPositions;
	//MinimizerReadMap _minimizer_to_readPositions;

	void indexReads(){

		KminmerParserParallel parser(_inputDir + "/read_data_init.txt", _minimizerSize, _kminmerSize, false, true, _nbCores);
		parser._densityThreshold = _minimizerDensity_assembly;
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

			ReadType readIndex = kminmerList._readIndex;
			if(readIndex % 100000 == 0)  Logger::get().debug() << "\tIndexing reads: " << Utils::getProgress(readIndex, _parent._totalNbReads);

			MinimizerRead read = {kminmerList._readIndex, kminmerList._readMinimizers, kminmerList._minimizerPos, kminmerList._readQualities, kminmerList._readMinimizerDirections, kminmerList._readLength};
			MinimizerRead readLowDensity = Utils::getLowDensityMinimizerRead(read, _parent._minimizerDensity_assembly);
			
			if(_parent.isReadTooShort(readLowDensity)) return;
			//if(readLowDensity._minimizers.size() < _parent._minReadLength) return;
			//if(kminmerList._readLength < _parent._minOverlapLength) return;
			//if(read._minimizersPos.size() == 0 || read._minimizersPos[read._minimizersPos.size()-1] < _parent._minReadLength) return;
			//if(readWorker._referenceRead._minimizers.size() < _parent._minReadLength) return;
			//if(read._minimizers.size() > 5000) return;

			//const MinimizerRead& referenceRead = Utils::getLowDensityMinimizerRead(readWorker._referenceRead, _parent._minimizerDensityAssembly);
			//if(kminmerList._readMinimizers.size() < _parent._minReadLength) return;
			//if(!Utils::isValidPositions(kminmerList._minimizerPos)) return;

			for(u_int32_t i=0; i<kminmerList._kminmersInfo.size(); i++){
			
				
				const ReadKminmerComplete& kminmerInfo = kminmerList._kminmersInfo[i];

				const KmerVec& vec = kminmerInfo._vec;
				
				u_int16_t positionIndex = i;
				//const u_int32_t& readPosition1 = kminmerList._minimizerPos[kminmerInfo._read_pos_start];
				//const u_int32_t& readPosition2 = kminmerList._minimizerPos[kminmerInfo._read_pos_end];

				//MinimizerType minimizer = kminmerList._readMinimizers[i];
				//int32_t position_1 = kminmerList._minimizerPos[i];
				//int32_t position_index_1 = i;
				//bool isReversed_1 = kminmerList._readMinimizerDirections[i];

				//int32_t position_2 = kminmerList._minimizerPos[i+1];
				//int32_t position_index_2 = i+1;
				//bool isReversed_2 = kminmerList._readMinimizerDirections[i+1];

				//MinimizerPairPosition minmizerPosition = {kminmerList._readIndex, position_1, position_index_1, isReversed_1, position_2, position_index_2, isReversed_2};

				//_parent._minimizer_to_readPositions[minimizer].push_back(minmizerPosition);
				//_minimizer_to_readPositions[minimizer].push_back(minmizerPosition);
				_parent._kminmer_to_readIndex.modify_if(vec, 
					[&readIndex, &positionIndex](KminmerPosMap::value_type& v) { 
					
					//if(std::find(v.second.begin(), v.second.end(), readIndex) == v.second.end()){
					//	v.second.push_back(readIndex);
					//}

					v.second.push_back({readIndex, positionIndex});
						
				});

			}
			
			
		}
		
	};


	u_int64_t _correctionCheckSum;
	u_int64_t _alignmentCheckSum;

	bool _eval_correction;


	void alignReadsLowDensity(){

		_alignmentFile = ofstream(_alignmentFilename);

		KminmerParserParallel parser(_inputDir + "/read_data_init.txt", _minimizerSize, _kminmerSize, false, true, _nbCores);
		parser._densityThreshold = _minimizerDensity_assembly;
		parser.parseSequences(AlignReadsLowDensityFunctor(*this));

		_alignmentFile.close();
	}


	class AlignReadsLowDensityFunctor {

		public:


		struct AlignmentScore{
			ReadType _queryReadIndex;
			float _score;
			//float _identity;
		};

		struct AlignmentScoreComparator {
			bool operator()(AlignmentScore const& a, AlignmentScore const& b){
				if(a._score == b._score){
					return a._queryReadIndex < b._queryReadIndex;
				}
				return a._score > b._score;
			}
		};

		typedef priority_queue<AlignmentScore, vector<AlignmentScore> , AlignmentScoreComparator> AlignmentScoreQueue;

		ReadCorrection& _parent;
		MinimizerAligner* _minimizerAligner;
		MinimizerChainer* _minimizerChainer;
		
		
		AlignReadsLowDensityFunctor(ReadCorrection& parent) : _parent(parent){
			_minimizerAligner = nullptr;
			_minimizerChainer = nullptr;
		}

		AlignReadsLowDensityFunctor(const AlignReadsLowDensityFunctor& copy) : _parent(copy._parent){
			_minimizerAligner = new MinimizerAligner(3, -1, -1, -1);
			_minimizerChainer = new MinimizerChainer(_parent._minimizerSize);
		}

		~AlignReadsLowDensityFunctor(){
			if(_minimizerAligner != nullptr) delete _minimizerAligner;
			if(_minimizerChainer != nullptr) delete _minimizerChainer;
		}


		void operator () (const KminmerList& kminmerList) {

			if(kminmerList._readIndex % 10000 == 0)  Logger::get().debug() << "\tAligning reads (low density): " << Utils::getProgress(kminmerList._readIndex, _parent._totalNbReads);

			MinimizerRead referenceRead = {kminmerList._readIndex, kminmerList._readMinimizers, kminmerList._minimizerPos, kminmerList._readQualities, kminmerList._readMinimizerDirections, kminmerList._readLength};
			MinimizerRead readLowDensity = Utils::getLowDensityMinimizerRead(referenceRead, _parent._minimizerDensity_assembly);
			
			//if(readLowDensity._minimizers.size() < _parent._minReadLength){
			if(_parent.isReadTooShort(readLowDensity)){
				//_parent.writeRead(referenceRead._readIndex, referenceRead._minimizers, referenceRead._qualities);
				_parent.writeAlignments(referenceRead._readIndex, {});
				return;
			}

			ReadMinimizerPositionMap referenceReadMinimizerPositionMap;
			//unordered_set<MinimizerType> minimizerSet;
			//unordered_map<MinimizerType, u_int8_t> readMinimizer_to_direction;

			
			unordered_map<MinimizerType, AlignmentScoreQueue> minimizer_to_alignmentScoreQueue;

			for(u_int32_t i=0; i<referenceRead._minimizers.size(); i++){

				MinimizerType minimizer = referenceRead._minimizers[i];
				u_int32_t position = referenceRead._minimizersPos[i];
				bool isReversed = referenceRead._readMinimizerDirections[i];

				MinimizerPosition minmizerPosition = {position, i, isReversed};
				referenceReadMinimizerPositionMap[minimizer].push_back(minmizerPosition);

				minimizer_to_alignmentScoreQueue[minimizer] = AlignmentScoreQueue{};
			}


			vector<ReadKminmerComplete> kminmersInfos;
			MDBG::getKminmers_complete(_parent._kminmerSize, referenceRead._minimizers, referenceRead._minimizersPos, kminmersInfos, referenceRead._readIndex, referenceRead._qualities);

			
			unordered_set<KmerVec> uniqueReferenceKminmers;
			for(size_t i=0; i<kminmersInfos.size(); i++){

				const ReadKminmerComplete& kminmerInfo = kminmersInfos[i];
				const KmerVec& vec = kminmerInfo._vec;
				uniqueReferenceKminmers.insert(vec);
			}
			

			vector<AlignmentResult2> alignments;

			unordered_map<ReadType, ReadMatchBound> referenceReadIndex_to_anchors2;
			unordered_map<ReadType, std::pair<u_int16_t, u_int16_t>> queryReadIndex_to_anchorBounds;

			for(const KmerVec& vec : uniqueReferenceKminmers){


				if(_parent._kminmer_to_readIndex.find(vec) == _parent._kminmer_to_readIndex.end()) continue;

				for(const MinimizerPairPosition& queryReadPosition : _parent._kminmer_to_readIndex[vec]){
					if(referenceRead._readIndex == queryReadPosition._readIndex) continue;


					if(referenceReadIndex_to_anchors2.find(queryReadPosition._readIndex) == referenceReadIndex_to_anchors2.end()){
						referenceReadIndex_to_anchors2[queryReadPosition._readIndex] = ReadMatchBound();
					}

					ReadMatchBound& bound = referenceReadIndex_to_anchors2[queryReadPosition._readIndex];
					
					if(queryReadPosition._positionIndex < bound._minIndex){
						bound._minIndex = queryReadPosition._positionIndex;
					}
					if(queryReadPosition._positionIndex > bound._maxIndex){
						bound._maxIndex = queryReadPosition._positionIndex;
					}

					bound._nbMatches += 1;

				}

			}

			
			
			for(const auto& it : referenceReadIndex_to_anchors2){

				ReadType queryReadIndex = it.first;

				const ReadMatchBound& readMatchBound = it.second;

				if(readMatchBound._nbMatches < 1) continue;

				const MinimizerRead& queryRead = _parent._mReads[queryReadIndex];

				//cout << referenceRead._minimizers.size() << " " << queryRead._minimizers.size() << endl;
				int64_t startQueryPositionIndex = getStartQueryPositionIndex(queryRead, readMatchBound._minIndex, 5000, referenceReadMinimizerPositionMap);
				int64_t endQueryPositionIndex = getEndQueryPositionIndex(queryRead, readMatchBound._maxIndex, 5000, referenceReadMinimizerPositionMap);
				
				//int64_t bestPossibleAlignmentScore = getBestPossibleAlignmentScore(queryRead, startQueryPositionIndex, endQueryPositionIndex, referenceReadMinimizerPositionMap);

				//if(!canBeBetterAlignment(queryRead, bestPossibleAlignmentScore, minimizer_to_alignmentScoreQueue)){
					//cout << "Filtered alignment: " << bestPossibleAlignmentScore << endl;
				//	continue;
				//}

				//startQueryPositionIndex = 0;
				//endQueryPositionIndex = queryRead._minimizers.size()-1;

				//_minimizerChainer->_anchorIndex = 0;
				vector<Anchor> anchors;

				//cout << startQueryPositionIndex << " " << endQueryPositionIndex << endl;
				for(size_t i=startQueryPositionIndex; i<=endQueryPositionIndex; i++){

					MinimizerType queryMinimizer = queryRead._minimizers[i];

					if(referenceReadMinimizerPositionMap.find(queryMinimizer) == referenceReadMinimizerPositionMap.end()) continue;

					u_int32_t queryPosition = queryRead._minimizersPos[i];
					bool queryIsReversed = queryRead._readMinimizerDirections[i];

					const vector<MinimizerPosition>& referenceMinimizerPositions = referenceReadMinimizerPositionMap[queryMinimizer];
					
					for(const MinimizerPosition& referenceMinimizerPosition : referenceMinimizerPositions){
						anchors.push_back({referenceMinimizerPosition._position, queryPosition, referenceMinimizerPosition._isReversed != queryIsReversed, referenceMinimizerPosition._positionIndex, i});
					}



				}



				//cout << "Anchors: " << anchors.size() << endl;



				AlignmentResult2 chainingAlignment = _minimizerChainer->computeChainingAlignment(anchors, referenceRead, queryRead, _minimizerAligner, _parent._maxChainingBand_lowDensity);
				
				//cout << referenceRead._readIndex << " " << queryRead._readIndex << " " << chainingAlignment._overHangStart << " " << chainingAlignment._overHangEnd << " " << chainingAlignment._nbMatches << " " << chainingAlignment._nbMissmatches << endl;
				//getchar();
				//AlignmentResult2 chainingAlignment = _parent.computeAlignment(anchors, referenceRead, queryRead, _minimizerAligner);
				
				//cout << mashDistance << " " << chainingAlignment._divergence << endl;
				//if(chainingAlignment._nbMatches - chainingAlignment._nbMissmatches - chainingAlignment._nbInsertions - chainingAlignment._nbDeletions < 0) continue;
				//if(chainingAlignment._nbMatches - chainingAlignment._nbMissmatches < 5) continue;
				//if(chainingAlignment._nbMatches < 1000*_minimizerDensity_correction) continue;
				//if(chainingAlignment._overHangStart > 5000) continue;
				//if(chainingAlignment._overHangEnd > 5000) continue;
				//if(chainingAlignment._nbMatches - chainingAlignment._nbMissmatches < 0) continue;
				//if(chainingAlignment._divergence > 0.04) continue;

				if(chainingAlignment._alignments.empty()) continue;

				addAlignmentScore(chainingAlignment, queryRead, minimizer_to_alignmentScoreQueue);
				//alignments.push_back(chainingAlignment);
				
			}
			
			/*
			if(_eval_correction){
				
				if(_simulatedReadMetadata[read._readIndex]._identity < 98) return read;
				if(_simulatedReadMetadata[read._readIndex]._isReversed) return read;
				
				cout << endl << "Read: " <<  read._readIndex << " " << read._minimizers.size() << " " << _mReads[read._readIndex]._meanReadQuality << endl;
				std::sort(alignments.begin(), alignments.end(), [](const AlignmentResult2& a, const AlignmentResult2& b){
					return a.getScore() > b.getScore();
				});

				cout << "Truth info:" << endl;
				const SimulatedReadMetadata& refMetadata = _simulatedReadMetadata[read._readIndex];
				cout << refMetadata.toString() << endl;

				cout << "Alignment info:" << endl;
				for(size_t i=0; i<alignments.size() && i < 50; i++){
					const AlignmentResult2& alignment = alignments[i];
					const SimulatedReadMetadata& metadata = _simulatedReadMetadata[alignment._queryReadIndex];
					cout << "\t" << i << ":\t" << metadata.toString() << "\t" << alignment._nbMatches << "\t" << alignment._nbMissmatches << "\t" << alignment._nbInsertions << "\t" << alignment._nbDeletions << "\t" << metadata.distanceFrom(refMetadata) << "\t" << alignment._divergence << "\t" << alignment.getSimilarity()<< endl;
				}

				getchar();
				
				//if(refMetadata._identity > 98.8) getchar();
				
				

			}
			*/

			//cout << referenceRead._readIndex << " " << alignments.size() << endl;
			//if(alignments.size() == 0){
			//	_parent.writeAlignments(referenceRead._readIndex, {});
			//	return;
			//}
			//cout << alignments.size() << endl;
			//cout << "copy here" << endl;
			unordered_set<ReadType> selectedQueryReadIndex;
			vector<AlignmentResult2> bestAlignmentsLowDensity;// = selectBestAlignments(referenceRead, alignments, referenceReadMinimizerPositionMap, _parent._usedCoverageLowDensity);
			
			for(auto& it : minimizer_to_alignmentScoreQueue){

				const MinimizerType& minimizer = it.first;
				AlignmentScoreQueue& queue = it.second;

				//cout << "\t----" << endl;
				int i =0;
				while(queue.size() > 0){
					const AlignmentScore& alignmentScore = queue.top();
					selectedQueryReadIndex.insert(alignmentScore._queryReadIndex);
					queue.pop();
					// << "\t" << i << ": " << alignmentScore._queryReadIndex << " " << alignmentScore._score << endl;
					i += 1;
				}

			}
			//const MinimizerRead& correctedRead = performPoaCorrection4(referenceRead, bestAlignmentsLowDensity, referenceReadMinimizerPositionMap);
			
			//_parent.writeRead(correctedRead._readIndex, correctedRead._minimizers, correctedRead._qualities);

			vector<ReadType> alignedQueryReads;
			//for(const AlignmentResult2& alignment : bestAlignmentsLowDensity){
			//	alignedQueryReads.push_back(alignment._queryReadIndex);
			//}

			for(const ReadType& queryReadIndex : selectedQueryReadIndex){
				alignedQueryReads.push_back(queryReadIndex);
			}

			//cout << referenceRead._readIndex << " " << alignedQueryReads.size() << endl;

			_parent.writeAlignments(referenceRead._readIndex, alignedQueryReads);
		
		
		

		}

		int64_t getStartQueryPositionIndex(const MinimizerRead& queryRead, int64_t startIndex, int64_t maxGapLength, ReadMinimizerPositionMap& referenceReadMinimizerPositionMap){
			
			//maxGapLength *= 2;
			int64_t lastMacthingIndex = startIndex;
			int64_t i = startIndex;
			int64_t currentGapLength = 0;
			int64_t prevPosition = queryRead._minimizersPos[i];
			i -= 1;

			while(true){
				if(i < 0) break;
				if(currentGapLength > maxGapLength) break;

				const MinimizerType currentMinimizer = queryRead._minimizers[i];

				if(referenceReadMinimizerPositionMap.find(currentMinimizer) == referenceReadMinimizerPositionMap.end()){
					currentGapLength += (prevPosition-queryRead._minimizersPos[i]);
				}
				else{
					currentGapLength = 0;
					lastMacthingIndex = i;
				}


				prevPosition = queryRead._minimizersPos[i];
				i -= 1;
			}

			return lastMacthingIndex;
		}

		int64_t getEndQueryPositionIndex(const MinimizerRead& queryRead, int64_t startIndex, int64_t maxGapLength, ReadMinimizerPositionMap& referenceReadMinimizerPositionMap){
			
			int64_t lastMacthingIndex = startIndex;
			int64_t i = startIndex;
			int64_t currentGapLength = 0;
			int64_t prevPosition = queryRead._minimizersPos[i];
			i += 1;

			while(true){
				if(i >= queryRead._minimizers.size()) break;
				if(currentGapLength > maxGapLength) break;

				const MinimizerType currentMinimizer = queryRead._minimizers[i];

				if(referenceReadMinimizerPositionMap.find(currentMinimizer) == referenceReadMinimizerPositionMap.end()){
					currentGapLength += (queryRead._minimizersPos[i]-prevPosition);
				}
				else{
					currentGapLength = 0;
					lastMacthingIndex = i;
				}


				prevPosition = queryRead._minimizersPos[i];
				i += 1;
			}

			return lastMacthingIndex;
		}
	

		int64_t getBestPossibleAlignmentScore(const MinimizerRead& queryRead, int64_t startQueryPositionIndex, int64_t endQueryPositionIndex, ReadMinimizerPositionMap& referenceReadMinimizerPositionMap){

			int64_t nbMatches = 0;
			u_int64_t startPositionReference = -1;
			u_int64_t endPositionReference = 0;
			
			for(size_t i=startQueryPositionIndex; i<=endQueryPositionIndex; i++){

				MinimizerType queryMinimizer = queryRead._minimizers[i];

				if(referenceReadMinimizerPositionMap.find(queryMinimizer) == referenceReadMinimizerPositionMap.end()) continue;

				nbMatches += 1;
			}

			for(const MinimizerPosition& referenceMinimizerPosition : referenceReadMinimizerPositionMap[queryRead._minimizers[startQueryPositionIndex]]){
				u_int64_t referencePositionIndex = referenceMinimizerPosition._positionIndex;
				if(referencePositionIndex < startPositionReference){
					startPositionReference = referencePositionIndex;
				}
			}

			for(const MinimizerPosition& referenceMinimizerPosition : referenceReadMinimizerPositionMap[queryRead._minimizers[endQueryPositionIndex]]){
				u_int64_t referencePositionIndex = referenceMinimizerPosition._positionIndex;
				if(referencePositionIndex > endPositionReference){
					endPositionReference = referencePositionIndex;
				}
			}

			int64_t referenceSize = endPositionReference - startPositionReference;
			int64_t nbDeletions = referenceSize - nbMatches;


			int64_t querySize = endQueryPositionIndex - startQueryPositionIndex;
			int64_t nbInsertions = querySize - nbMatches;

			return nbMatches - nbDeletions - nbInsertions;
		}

		bool canBeBetterAlignment(const MinimizerRead& queryRead, const u_int32_t estimatedAlignmentScore, unordered_map<MinimizerType, AlignmentScoreQueue>& minimizer_to_alignmentScoreQueue){

			for(MinimizerType minimizer : queryRead._minimizers){
				if(minimizer_to_alignmentScoreQueue.find(minimizer) == minimizer_to_alignmentScoreQueue.end()) continue;

				AlignmentScoreQueue& queue = minimizer_to_alignmentScoreQueue[minimizer];

				if(queue.size() < _parent._usedCoverageForCorrection){
					return true;
				}

				const AlignmentScore& worseAlignmentScore = queue.top();

				if(estimatedAlignmentScore < worseAlignmentScore._score) continue;

				return true;
			}

			return false;
		}

		void addAlignmentScore(const AlignmentResult2& alignment, const MinimizerRead& queryRead, unordered_map<MinimizerType, AlignmentScoreQueue>& minimizer_to_alignmentScoreQueue){


			float score = alignment._nbMatches - alignment._nbMissmatches - alignment._nbInsertions - alignment._nbDeletions;
			const AlignmentScore currentAlignmentScore = {alignment._queryReadIndex, score};

			for(MinimizerType minimizer : queryRead._minimizers){
				if(minimizer_to_alignmentScoreQueue.find(minimizer) == minimizer_to_alignmentScoreQueue.end()) continue;

				AlignmentScoreQueue& queue = minimizer_to_alignmentScoreQueue[minimizer];

				if(queue.size() < _parent._usedCoverageForCorrection){
					queue.push(currentAlignmentScore);
					continue;
				}

				const AlignmentScore& worseAlignmentScore = queue.top();

				if(currentAlignmentScore._score < worseAlignmentScore._score) continue;

				if(currentAlignmentScore._score == worseAlignmentScore._score){
					if(currentAlignmentScore._queryReadIndex > worseAlignmentScore._queryReadIndex) continue;
				}
				
				queue.pop();
				queue.push(currentAlignmentScore);
			}
				

		}
		/*
		struct AlignmentScore{
			ReadType _queryReadIndex;
			float _score;
			float _identity;
		};
		
		vector<AlignmentResult2> selectBestAlignments(const MinimizerRead& referenceRead, const vector<AlignmentResult2>& alignments, ReadMinimizerPositionMap& referenceReadMinimizerPositionMap, const size_t maxCoverage){
			

			vector<AlignmentResult2> selectedAlignments;
			
			
			unordered_map<MinimizerType, vector<AlignmentScore>> minimizer_to_alignmentResultsNew;
			
			for(size_t i=0; i<alignments.size(); i++){

				const AlignmentResult2& alignment = alignments[i];

				float score = alignment._nbMatches - alignment._nbMissmatches - alignment._nbInsertions - alignment._nbDeletions;
				
				for(MinimizerType minimizer : _parent._mReads[alignment._queryReadIndex]._minimizers){
					if(referenceReadMinimizerPositionMap.find(minimizer) == referenceReadMinimizerPositionMap.end()) continue;

					minimizer_to_alignmentResultsNew[minimizer].push_back({alignment._queryReadIndex, score, alignment._identity});
				}
				
			}


			unordered_set<ReadType> selectedReads;

			for(const auto& it : minimizer_to_alignmentResultsNew){
				MinimizerType minimizer = it.first;
				const vector<AlignmentScore>& alignmentResults = it.second;

				const vector<AlignmentScore>& bestAlignmentResults  = getBestAlignments(alignmentResults, maxCoverage);

				for(const AlignmentScore& alignmentResult : bestAlignmentResults){
					//if(alignmentResult.divergence() > 0.02) continue;
					//cout << alignmentResult._readIndex << " " << alignmentResult.score() << " " << readIndex_to_matchScore[alignmentResult._readIndex] << endl;
					selectedReads.insert(alignmentResult._queryReadIndex);
				}

				//cout << "check" << endl;
			}



			for(const AlignmentResult2& alignment : alignments){

				if(selectedReads.find(alignment._queryReadIndex) == selectedReads.end()) continue;

				selectedAlignments.push_back(alignment);

			}

			//cout << "lala: " << selectedReads.size() << endl;
			return selectedAlignments;

		}
		*/
		/*
		vector<AlignmentResult2> selectBestAlignments2(const MinimizerRead& referenceRead, const vector<AlignmentResult2>& alignments, ReadMinimizerPositionMap& referenceReadMinimizerPositionMap, const size_t maxCoverage){
			

			vector<AlignmentResult2> selectedAlignments;
			
			vector<vector<AlignmentScore>> position_to_alignmentScores(referenceRead._minimizers.size());

			for(const AlignmentResult2& alignment : alignments){

				float score = alignment._nbMatches - alignment._nbMissmatches - alignment._nbInsertions - alignment._nbDeletions;
				const MinimizerRead& queryRead = _parent._mReads[alignment._queryReadIndex];

				for(const auto& it : alignment._alignments){

					u_int32_t referencePosition = it.first;
					u_int32_t queryPosition = it.second;

					if(referencePosition == -1){ //insertion
					}
					else if(queryPosition == -1){ //deletion

					}
					else if(referenceRead._minimizers[referencePosition] == queryRead._minimizers[queryPosition]){ //Match
						position_to_alignmentScores[referencePosition].push_back({alignment._queryReadIndex, score, alignment._identity});
					}
					else { //missmatch
					}

				}
			}


			unordered_set<ReadType> selectedReads;

			for(size_t i=0; i<position_to_alignmentScores.size(); i++){
				//MinimizerType minimizer = it.first;
				const vector<AlignmentScore>& alignmentResults = position_to_alignmentScores[i];

				const vector<AlignmentScore>& bestAlignmentResults  = getBestAlignments(alignmentResults, maxCoverage);

				for(const AlignmentScore& alignmentResult : bestAlignmentResults){
					//if(alignmentResult.divergence() > 0.02) continue;
					//cout << alignmentResult._readIndex << " " << alignmentResult.score() << " " << readIndex_to_matchScore[alignmentResult._readIndex] << endl;
					selectedReads.insert(alignmentResult._queryReadIndex);
				}

				//cout << "check" << endl;
			}



			for(const AlignmentResult2& alignment : alignments){

				if(selectedReads.find(alignment._queryReadIndex) == selectedReads.end()) continue;

				selectedAlignments.push_back(alignment);

			}

			//cout << "lala: " << selectedReads.size() << endl;
			return selectedAlignments;

		}
		*/
		/*
		vector<AlignmentScore> getBestAlignments(vector<AlignmentScore> alignmentResults, int n){

			vector<AlignmentScore> bestAlignments;

			std::sort(alignmentResults.begin(), alignmentResults.end(), [](const AlignmentScore& a, const AlignmentScore& b){
				if(a._score == b._score){
					//if(a._divergence == b._divergence){
					return a._queryReadIndex < b._queryReadIndex;
					//}
					//return a._divergence < b._divergence;
				}
				return a._score > b._score;
			});

			for(int i=0; i<n && i<alignmentResults.size(); i++){

				bool isHere = false;
				for(const AlignmentScore& al : bestAlignments){
					if(al._queryReadIndex == alignmentResults[i]._queryReadIndex){ //The same alignment can be present several time if a minimizer is repeated in a read
						isHere = true;
						break;
					}
				}

				if(isHere) continue;

				bestAlignments.push_back(alignmentResults[i]);
			}

			return bestAlignments;
		}
		*/

	};

    struct AlignmentWriter{
        ReadType _referenceReadIndex;
        vector<ReadType> _alignedQueryReads;
    };

    struct AlignmentWriter_Comparator {
        bool operator()(const AlignmentWriter& p1, const AlignmentWriter& p2){
            return p1._referenceReadIndex > p2._referenceReadIndex;
        }
    };

	priority_queue<AlignmentWriter, vector<AlignmentWriter> , AlignmentWriter_Comparator> _alignmentWriterQueue;
	ReadType _nextAlignmentReadIndexWriter;

	void writeAlignments(ReadType referenceReadIndex, const vector<ReadType>& alignedQueryReads){

		#pragma omp critical
		{
			
			_alignmentWriterQueue.push({referenceReadIndex, alignedQueryReads});

			while(!_alignmentWriterQueue.empty()){


				const AlignmentWriter& readWriter = _alignmentWriterQueue.top();

				if(readWriter._referenceReadIndex == _nextAlignmentReadIndexWriter){


					const ReadType referenceReadIndexCurrent = readWriter._referenceReadIndex;
					const vector<ReadType> alignedQueryReadsCurrent = readWriter._alignedQueryReads;
					const u_int16_t nbAlignedReads = alignedQueryReadsCurrent.size();

					if(nbAlignedReads > 0){

						for(ReadType queryReadIndex : alignedQueryReadsCurrent){
							_alignmentCheckSum += referenceReadIndexCurrent * nbAlignedReads * queryReadIndex;
						}

						_alignmentFile.write((const char*)&referenceReadIndexCurrent, sizeof(referenceReadIndexCurrent));
						_alignmentFile.write((const char*)&nbAlignedReads, sizeof(nbAlignedReads));
						_alignmentFile.write((const char*)&alignedQueryReadsCurrent[0], nbAlignedReads*sizeof(ReadType));
						

						if(readWriter._referenceReadIndex % 10000 == 0){
							Logger::get().debug() << "\tAlignment checksum: " << readWriter._referenceReadIndex << " " << _alignmentCheckSum;
							//if(readWriter._referenceReadIndex > 1000) getchar();
						}
						
					}


					
					_alignmentWriterQueue.pop();
					_nextAlignmentReadIndexWriter += 1;
				}
				else{
					break;
				}
			}
			
		}

	}













	vector<vector<ReadType>> _lowDensityAlignments;

	void loadAlignments(bool useReadPartionning){

		_lowDensityAlignments.clear();
		ifstream alignmentFile(_alignmentFilename);


		while(true){

			ReadType referenceReadIndex;
			u_int16_t nbAlignedReads;
			vector<ReadType> alignedQueryReads;

			alignmentFile.read((char*)&referenceReadIndex, sizeof(referenceReadIndex));
			

			if(alignmentFile.eof()) break;

			alignmentFile.read((char*)&nbAlignedReads, sizeof(nbAlignedReads));

			alignedQueryReads.resize(nbAlignedReads);
			alignmentFile.read((char*)&alignedQueryReads[0], alignedQueryReads.size() * sizeof(ReadType));

			while(_lowDensityAlignments.size() <= referenceReadIndex){
				_lowDensityAlignments.push_back({});
			}

			if(useReadPartionning && !_isReadToCorrect[referenceReadIndex]){
				//_lowDensityAlignments.push_back({});
			}
			else{
				_lowDensityAlignments[_lowDensityAlignments.size()-1] = alignedQueryReads;
			}

			//cout << referenceReadIndex << " " << alignedQueryReads.size() << endl;
        }

		alignmentFile.close();

		//cout << _lowDensityAlignments.size() << endl;
		//getchar();
	}

	struct ReadPartition{
		vector<ReadType> _readsToLoad;
		vector<ReadType> _readsToCorrect;
	};

	//vector<ReadPartition> _readPartitions;

	int partitionReads(){
		
		ofstream partitionFile(_partitionFilename);

		ReadPartition currentPartition;

		vector<bool> isCorrected(_lowDensityAlignments.size(), false);
		vector<bool> isVisited(_lowDensityAlignments.size(), false);
		size_t nbPartitions = 0;

		for(ReadType readIndex = 0; readIndex < _lowDensityAlignments.size(); readIndex++){

			//if(isVisited[readIndex]) continue;
			if(isCorrected[readIndex]) continue;

			currentPartition._readsToLoad.push_back(readIndex);
			isVisited[readIndex] = true;

			queue<ReadType> queue;
			queue.push(readIndex);

			while (!queue.empty()) {
				
				ReadType currentReadIndex = queue.front();
				queue.pop();

				if(isCorrected[currentReadIndex]) continue;

				//currentPartition._readsToLoad.push_back(currentReadIndex);
				currentPartition._readsToCorrect.push_back(currentReadIndex);
				isCorrected[currentReadIndex] = true;
				
				for(const ReadType readIndexNeighbor : _lowDensityAlignments[currentReadIndex]) {
					if(isVisited[readIndexNeighbor]) continue;
					//if(isCorrected[readIndexNeighbor]) continue;

					currentPartition._readsToLoad.push_back(readIndexNeighbor);
					isVisited[readIndexNeighbor] = true;
					queue.push(readIndexNeighbor);
				}
				
				if(currentPartition._readsToLoad.size() > 2000000) break;
			}


			if(currentPartition._readsToLoad.size() > 2000000){
				Logger::get().debug() << "\tAdd partition: " << currentPartition._readsToCorrect.size() << " " << currentPartition._readsToLoad.size();
				/*
				unordered_set<ReadType> realNbreadToLoad;
				for(ReadType readIndex : currentPartition._readsToCorrect){
					realNbreadToLoad.insert(readIndex);
					for(const ReadType readIndexNeighbor : _lowDensityAlignments[readIndex]) {
						realNbreadToLoad.insert(readIndexNeighbor);
					}
				}

				//cout << currentPartition._readsToLoad.size() << " " << realNbreadToLoad.size() << endl;

				if(currentPartition._readsToLoad.size() != realNbreadToLoad.size()){
					cout << "issue partitionning" << endl;
					exit(1);
				}
				*/

				//getchar();
				nbPartitions += 1;
				writePartition(partitionFile, currentPartition);
				//_readPartitions.push_back(currentPartition);
				currentPartition._readsToCorrect.clear();
				currentPartition._readsToLoad.clear();
				for(size_t i=0; i<isVisited.size(); i++){
					isVisited[i] = false;
				}
				
			}

		}

		if(currentPartition._readsToLoad.size() > 0){
			Logger::get().debug() << "\tAdd partition: " << currentPartition._readsToCorrect.size() << " " << currentPartition._readsToLoad.size();

			nbPartitions += 1;
			writePartition(partitionFile, currentPartition);
			//_readPartitions.push_back(currentPartition);
			currentPartition._readsToCorrect.clear();
			currentPartition._readsToLoad.clear();
			for(size_t i=0; i<isVisited.size(); i++){
				isVisited[i] = false;
			}
			
		}

		Logger::get().debug() << "\tNb partitions: " << nbPartitions; 
		//for(size_t i=0; i<_readPartitions.size(); i++){
		//	cout << "\t" << i << "\t" << _readPartitions[i]._readsToCorrect.size() << " " << _readPartitions[i]._readsToLoad.size() << endl;
		//}

		//getchar();
		partitionFile.close();

		return nbPartitions;
	}

	void writePartition(ofstream& partitionFile, const ReadPartition& partition){

		u_int64_t readToLoadSize = partition._readsToLoad.size();
		u_int64_t readToCorrectSize = partition._readsToCorrect.size();

		partitionFile.write((const char*)&readToLoadSize, sizeof(readToLoadSize));
		partitionFile.write((const char*)&partition._readsToLoad[0], readToLoadSize*sizeof(ReadType));

		partitionFile.write((const char*)&readToCorrectSize, sizeof(readToCorrectSize));
		partitionFile.write((const char*)&partition._readsToCorrect[0], readToCorrectSize*sizeof(ReadType));				

	}

	void correctReads(){

		

		//AlignmentFile* alignmentFile = new AlignmentFile(_inputDir, _nbReads, _nbReadsPerPartition, 'r');
		//alignmentFile->setupLinearReading();
		//alignmentFile->readAll();


		int nbCores = _nbCores;
		if(_print_debug || _eval_correction){  
			nbCores = 1;
		}


		//nbCores = 1; cout << "correction one cores" << endl;

		KminmerParserParallel parser(_inputDir + "/read_data_init.txt", _minimizerSize, _kminmerSize, false, true, nbCores);
		//parser2._alignmentFile = alignmentFile;
		//parser2.parse(FilterKminmerFunctor2(*this));
		parser.parseSequences(ReadCorrectionFunctor(*this));

		//alignmentFile->closeLinearReading();
		//delete alignmentFile;


	}

	class ReadCorrectionFunctor {

		public:

		ReadCorrection& _parent;
		//std::unique_ptr<spoa64::AlignmentEngine> _alignmentEngine;// = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 3, -5, -4);
		//std::unique_ptr<spoa64::AlignmentEngine> _alignmentEngineTrimming;
		//std::unique_ptr<spoa64::AlignmentEngine> _alignmentEngine_ov;
		MinimizerAligner* _minimizerAligner;
		MinimizerChainer* _minimizerChainer;
		
		
		ReadCorrectionFunctor(ReadCorrection& parent) : _parent(parent){
			_minimizerAligner = nullptr;
			_minimizerChainer = nullptr;
		}

		ReadCorrectionFunctor(const ReadCorrectionFunctor& copy) : _parent(copy._parent){
			//_alignmentEngine = spoa64::AlignmentEngine::Create(spoa64::AlignmentType::kNW, 3, -5, -4);
			//_alignmentEngineTrimming = spoa64::AlignmentEngine::Create(spoa64::AlignmentType::kSW, 3, -2, -2);
			//_alignmentEngine_ov = spoa64::AlignmentEngine::Create(spoa64::AlignmentType::kOV, 3, -1, -1, -1);
			_minimizerAligner = new MinimizerAligner(3, -1, -1, -1);
			_minimizerChainer = new MinimizerChainer(_parent._minimizerSize);
		}

		~ReadCorrectionFunctor(){
			if(_minimizerAligner != nullptr) delete _minimizerAligner;
			if(_minimizerChainer != nullptr) delete _minimizerChainer;
		}


		void operator () (const KminmerList& kminmerList) {
			


			//if(kminmerList._readIndex < 14995) return;

			//cout << kminmerList._readIndex << endl;
			//if(kminmerList._readIndex % 10000 == 0)  Logger::get().debug() << "\tCorrecting reads: " << Utils::getProgress(kminmerList._readIndex, _parent._totalNbReads);
			//if(kminmerList._readIndex % 10000 == 0) cout << "Correcting read: " << kminmerList._readIndex << endl;

			if(!_parent._isReadToCorrect[kminmerList._readIndex]) return;

			if(_parent._print_debug){
				cout << endl << endl;
				cout << "-------------------------------------------------------" << endl;
				cout << "Correcting read: " << "ReadIndex=" << kminmerList._readIndex << " Readsize=" << kminmerList._readMinimizers.size() << endl;
				for(MinimizerType m : kminmerList._readMinimizers){
					cout << m << endl;
				}
				cout << endl << endl;

			}


			MinimizerRead referenceRead = {kminmerList._readIndex, kminmerList._readMinimizers, kminmerList._minimizerPos, kminmerList._readQualities, kminmerList._readMinimizerDirections, kminmerList._readLength};
			MinimizerRead referenceReadLowDensity = Utils::getLowDensityMinimizerRead(referenceRead, _parent._minimizerDensity_assembly);
			
			//if(referenceReadLowDensity._minimizers.size() < _parent._minReadLength){
			if(_parent.isReadTooShort(referenceReadLowDensity)){
				//if(referenceRead._meanReadQuality >= _parent._minReadQuality) {
				_parent.writeRead(referenceRead._readIndex, referenceRead._minimizers, referenceRead._qualities);
				//}

			}
			else{
				MinimizerRead correctedRead = correctRead(referenceRead);
				_parent.writeRead(correctedRead._readIndex, correctedRead._minimizers, correctedRead._qualities);
			}


			
		}

		MinimizerRead correctRead(const MinimizerRead& referenceRead){
			

			
			const vector<ReadType>& alignedQueryReads = _parent._lowDensityAlignments[referenceRead._readIndex];
			if(alignedQueryReads.size() <= 0) return referenceRead;

			//if(alignedQueryReads.size() <= 2){
			//	if(referenceRead._meanReadQuality < _parent._minReadQuality) return {referenceRead._readIndex};
			//	return referenceRead;
			//}
			
			
			const vector<ReadType>& bestAlignments = filterAlignments(referenceRead, alignedQueryReads);
			if(bestAlignments.size() <= 0) return referenceRead;
			
			//const vector<AlignmentResult2>& bestAlignmentsLowDensity = convertToLowDensityAlignments(referenceReadLowDensity, alignedQueryReads);
			//float alignmentCoverage = estimateAlignmentCoverage(referenceRead, bestAlignments, referenceReadMinimizerPositionMapHighDensity);
			
			//if(alignmentCoverage <= 2){
				//if(referenceRead._meanReadQuality < _parent._minReadQuality) return {referenceRead._readIndex};
			//	return referenceRead;
			//}
			//#pragma omp critical
			//{
			//cout << alignmentCoverage << endl;
			//}

			const MinimizerRead& referenceReadLowDensity = Utils::getLowDensityMinimizerRead(referenceRead, _parent._minimizerDensity_assembly);
			const MinimizerRead& correctedRead = performPoaCorrection4(referenceReadLowDensity, bestAlignments);
			
			//const MinimizerRead& correctedReadLow = performPoaCorrection5(referenceRead, bestAlignments);
			

			if(_parent._eval_correction){
				//_parent.evaluateCorrection(correctedRead, alOverlap);
			}
			
			return correctedRead;
		}

		

		vector<ReadType> filterAlignments(const MinimizerRead& referenceRead, const vector<ReadType>& alignedQueryReads){

			ReadMinimizerPositionMap referenceReadMinimizerPositionMap;


			for(u_int32_t i=0; i<referenceRead._minimizers.size(); i++){

				MinimizerType minimizer = referenceRead._minimizers[i];
				u_int32_t position = referenceRead._minimizersPos[i];
				bool isReversed = referenceRead._readMinimizerDirections[i];

				MinimizerPosition minmizerPosition = {position, i, isReversed};
				referenceReadMinimizerPositionMap[minimizer].push_back(minmizerPosition);
			}

			//cout << referenceReadHighDensity._readIndex << " " << referenceReadHighDensity._minimizers.size() << endl;
			vector<ReadType> alignments;

			//const MinimizerRead& referenceRead = _mReads[referenceReadIndex];
			//const MinimizerRead& referenceRead_lowDensity = Utils::getLowDensityMinimizerRead(_mReads[read._readIndex], _minimizerDensity_assembly);
			//const MinimizerRead& referenceRead_highDensity = _mReads[read._readIndex];
			//cout << read._readIndex << endl;



			for(const ReadType& queryReadIndex : alignedQueryReads){

				//if(queryReadIndex == referenceReadHighDensity._readIndex) continue;

				const MinimizerRead& queryReadHighDensity = _parent._mReads[queryReadIndex];

				vector<Anchor> anchors;
				//_minimizerChainer->_anchorIndex = 0;

				for(size_t i=0; i<queryReadHighDensity._minimizers.size(); i++){

					MinimizerType queryMinimizer = queryReadHighDensity._minimizers[i];

					if(referenceReadMinimizerPositionMap.find(queryMinimizer) == referenceReadMinimizerPositionMap.end()) continue;

					u_int32_t queryPosition = queryReadHighDensity._minimizersPos[i];
					bool queryIsReversed = queryReadHighDensity._readMinimizerDirections[i];

					const vector<MinimizerPosition>& referenceMinimizerPositions = referenceReadMinimizerPositionMap[queryMinimizer];
					
					for(const MinimizerPosition& referenceMinimizerPosition : referenceMinimizerPositions){
						anchors.push_back({referenceMinimizerPosition._position, queryPosition, referenceMinimizerPosition._isReversed != queryIsReversed, referenceMinimizerPosition._positionIndex, i});
					}




				}

			
				AlignmentResult2 chainingAlignment = _minimizerChainer->computeChainingAlignment(anchors, referenceRead, queryReadHighDensity, _minimizerAligner, _parent._maxChainingBand_highDensity);
				
				
				//cout << mashDistance << " " << chainingAlignment._divergence << endl;
				//if(chainingAlignment._nbMatches - chainingAlignment._nbMissmatches - chainingAlignment._nbInsertions - chainingAlignment._nbDeletions < 0) continue;
				//if(chainingAlignment._nbMatches - chainingAlignment._nbMissmatches < 5) continue;
				//if(chainingAlignment._nbMatches < 1000*_minimizerDensity_correction) continue;
				//if(chainingAlignment._divergence > 0.04) continue;
				if(chainingAlignment._overHangStart > 1000) continue;
				if(chainingAlignment._overHangEnd > 1000) continue;
				if(chainingAlignment._alignLength < _parent._minOverlapLength) continue;
				if(chainingAlignment._identity < _parent._minIdentity) continue;
				//if(chainingAlignment._divergence > 0.04) continue;
				if(chainingAlignment._alignments.empty()) continue;

				alignments.push_back(queryReadIndex);

			}

		

			//return _parent.selectBestAlignments(referenceReadHighDensity, alignments, referenceReadMinimizerPositionMap, _parent._usedCoverageForCorrection);

			
			//for(size_t i=0; i<alignments.size(); i++){

			//	const AlignmentResult2& alignment = alignments[i];

			//	alignments[i] = _parent.computeAlignment(_parent._mReads[alignment._referenceReadIndex], _parent._mReads[alignment._queryReadIndex], _minimizerAligner, alignment._isQueryReversed);
				
			//}
			

			return alignments;

			//return alignments;
		}
	
		/*
		float estimateAlignmentCoverage(const MinimizerRead& referenceRead, const vector<AlignmentResult2>& alignments, ReadMinimizerPositionMap& referenceReadMinimizerPositionMap){
			
			
			unordered_map<MinimizerType, u_int32_t> minimizerCounts;
			
			for(size_t i=0; i<alignments.size(); i++){

				const AlignmentResult2& alignment = alignments[i];

				for(MinimizerType minimizer : _parent._mReads[alignment._queryReadIndex]._minimizers){
					if(referenceReadMinimizerPositionMap.find(minimizer) == referenceReadMinimizerPositionMap.end()) continue;

					minimizerCounts[minimizer] += 1;
				}
				
			}

			vector<float> minimizerAbundances;

			for(const auto& it : minimizerCounts){
				minimizerAbundances.push_back(it.second);
			}

			float median = Utils::compute_median_float(minimizerAbundances);

			return median;
		}
		*/

		Graph* _debugDbgGraph;
		
		MinimizerRead performPoaCorrection4(const MinimizerRead& referenceRead, const vector<ReadType>& alignedQueryReads){

			Graph* dbgGraph = new Graph(referenceRead._minimizers, referenceRead._qualities);

			
			//vector<u_int64_t> weights(referenceRead._qualities.size(), 1);
			//for(size_t i=0; i<referenceRead._qualities.size(); i++){
			//	weights[i] = referenceRead._qualities[i];
			//}


			


			/*
			for(size_t i=0; i<alignments.size(); i++){

				const AlignmentResult2& alignment = alignments[i];
			
				const MinimizerRead& queryRead = _parent._mReads[alignment._queryReadIndex]; //Utils::getLowDensityMinimizerRead(_mReads[alignment._queryReadIndex], _minimizerDensity_assembly);
			
				dbgGraph->pileupAlignment2(alignment, referenceRead._minimizers, queryRead._minimizers);

			}

			unordered_set<ReadType> strainReadIndexes = dbgGraph->detectStrainReads(referenceRead._minimizers);

			//cout << "Nb reads: " << alignments.size() << endl;
			//cout << "Nb strain reads: " << strainReadIndexes.size() << endl;
			//getchar();

			vector<AlignmentResult2> filteredAlignments;

			for(const AlignmentResult2& alignment : alignments){
				if(strainReadIndexes.find(alignment._queryReadIndex) != strainReadIndexes.end()) continue;
				filteredAlignments.push_back(alignment);
			}

			//if(alignments.size() == filteredAlignments.size()) break;
			//alignments = filteredAlignments;
			

			//if(!dbgGraph->isLowCovered(read._minimizers)){
			//	delete dbgGraph;
			//	return read;
			//}
			//u_int64_t nbErroneousMinimizers =  computeNbErroneousMinimizers(read);
			//cout << "Nb erroneous minimizers: " << nbErroneousMinimizers << endl;
			//u_int64_t readAbundance = dbgGraph->computeReadAbundance(nbErroneousMinimizers, read._minimizers);
			//cout << "Read abundance: " << readAbundance << endl;
			//getchar();



			dbgGraph->clearPileup();
			*/
			
			//std::shuffle (alignments.begin(), alignments.end(), std::default_random_engine(42));


			//std::sort(alignments.begin(), alignments.end(), [](const MinimizerAlignment & a, const MinimizerAlignment & b){
			//	return a._alignmentScore > b._alignmentScore;
				//return a._lastMatchPos < b._lastMatchPos;
			//});


			ReadMinimizerPositionMap referenceReadMinimizerPositionMap;


			for(u_int32_t i=0; i<referenceRead._minimizers.size(); i++){

				MinimizerType minimizer = referenceRead._minimizers[i];
				u_int32_t position = referenceRead._minimizersPos[i];
				bool isReversed = referenceRead._readMinimizerDirections[i];

				MinimizerPosition minmizerPosition = {position, i, isReversed};
				referenceReadMinimizerPositionMap[minimizer].push_back(minmizerPosition);
			}

			//cout << referenceReadHighDensity._readIndex << " " << referenceReadHighDensity._minimizers.size() << endl;
			vector<ReadType> alignments;

			//const MinimizerRead& referenceRead = _mReads[referenceReadIndex];
			//const MinimizerRead& referenceRead_lowDensity = Utils::getLowDensityMinimizerRead(_mReads[read._readIndex], _minimizerDensity_assembly);
			//const MinimizerRead& referenceRead_highDensity = _mReads[read._readIndex];
			//cout << read._readIndex << endl;



			for(const ReadType& queryReadIndex : alignedQueryReads){

				//if(queryReadIndex == referenceReadHighDensity._readIndex) continue;

				//const MinimizerRead& queryReadHighDensity = _parent._mReads[queryReadIndex];
				const MinimizerRead& queryRead = Utils::getLowDensityMinimizerRead(_parent._mReads[queryReadIndex], _parent._minimizerDensity_assembly);
				


				vector<Anchor> anchors;
				//_minimizerChainer->_anchorIndex = 0;

				for(size_t i=0; i<queryRead._minimizers.size(); i++){

					MinimizerType queryMinimizer = queryRead._minimizers[i];

					if(referenceReadMinimizerPositionMap.find(queryMinimizer) == referenceReadMinimizerPositionMap.end()) continue;

					u_int32_t queryPosition = queryRead._minimizersPos[i];
					bool queryIsReversed = queryRead._readMinimizerDirections[i];

					const vector<MinimizerPosition>& referenceMinimizerPositions = referenceReadMinimizerPositionMap[queryMinimizer];
					
					for(const MinimizerPosition& referenceMinimizerPosition : referenceMinimizerPositions){
						anchors.push_back({referenceMinimizerPosition._position, queryPosition, referenceMinimizerPosition._isReversed != queryIsReversed, referenceMinimizerPosition._positionIndex, i});
					}




				}

			
				AlignmentResult2 alignment = _minimizerChainer->computeChainingAlignment(anchors, referenceRead, queryRead, _minimizerAligner, _parent._maxChainingBand_highDensity);
				

				//const MinimizerRead& queryRead = _parent._mReads[alignment._queryReadIndex];

				dbgGraph->addAlignment2(alignment, referenceRead._minimizers, queryRead._minimizers, queryRead._qualities);


			}
			
			//dbgGraph->save("/pasteur/zeus/projets/p02/seqbio/gbenoit/workspace/home/tmp/poaGraph.gfa", referenceRead._originalMinimizers);
			//getchar();
		
			MinimizerRead correctedRead = computePath2(referenceRead, dbgGraph, referenceReadMinimizerPositionMap);
			//dbgGraph->save("/pasteur/zeus/projets/p02/seqbio/gbenoit/workspace/home/tmp/poaGraph.gfa", referenceRead._originalMinimizers);
			//getchar();
			delete dbgGraph;
			
			/*
			int32_t firstMatchPosition = -1;
			int32_t lastMatchPosition = -1;
			for(size_t i=0; i<correctedRead._minimizers.size(); i++){
				if(referenceReadMinimizerPositionMap.find(correctedRead._minimizers[i]) == referenceReadMinimizerPositionMap.end()) continue;
				if(firstMatchPosition == -1) firstMatchPosition = i;
				lastMatchPosition = i;
			}

			vector<MinimizerType> minimizers;
			for(size_t i=firstMatchPosition; i<=lastMatchPosition; i++){
				minimizers.push_back(correctedRead._minimizers[i]);
			}

			correctedRead._minimizers = minimizers;
			*/

			const MinimizerRead& correctedReadTrimmed = trimCorrectedPath(referenceRead, correctedRead);
			
			return correctedReadTrimmed;
		}

		/*
		MinimizerRead performPoaCorrection5(const MinimizerRead& referenceReadHigh, const vector<AlignmentResult2>& alignmentsHigh){


			if(referenceReadHigh._minimizers.size() < 5){
				return referenceReadHigh;
			}

			const MinimizerRead& referenceRead = Utils::getLowDensityMinimizerRead(referenceReadHigh, _parent._minimizerDensity_assembly);
			if(referenceRead._minimizers.size() < 5) return referenceRead;

			Graph* dbgGraph = new Graph(referenceRead._minimizers, referenceRead._qualities);

			ReadMinimizerPositionMap referenceReadMinimizerPositionMapHighDensity;


			for(u_int32_t i=0; i<referenceRead._minimizers.size(); i++){

				MinimizerType minimizer = referenceRead._minimizers[i];
				u_int32_t position = 0;//referenceRead._minimizersPos[i];
				bool isReversed = false;//referenceRead._readMinimizerDirections[i];

				MinimizerPosition minmizerPosition = {position, i, isReversed};
				referenceReadMinimizerPositionMapHighDensity[minimizer].push_back(minmizerPosition);
			}

			for(size_t i=0; i<alignmentsHigh.size(); i++){

				const AlignmentResult2& alignmentHigh = alignmentsHigh[i];

				const MinimizerRead& queryRead = Utils::getLowDensityMinimizerRead(_parent._mReads[alignmentHigh._queryReadIndex], _parent._minimizerDensity_assembly);
				

				//const vector<MinimizerType> subContigMinimizers_forward(contigMinimizers.begin() + startQueryPositionIndex, contigMinimizers.begin() + endQueryPositionIndex);
				vector<MinimizerType> queryRead_reverse = queryRead._minimizers;
				std::reverse(queryRead_reverse.begin(), queryRead_reverse.end());

				//vector<MinimizerType> minimizerSequence_forward = mContig._contigMinimizers;
				//vector<MinimizerType> minimizerSequence_reverse = mContig._contigMinimizers;
				//std::reverse(minimizerSequence_reverse.begin(), minimizerSequence_reverse.end());

				//vector<u_int8_t> qualities(minimizerSequence_forward.size(), 1);

				//if(_print_debug){
				//	cout << endl << "\tAlign sequence: " << i << " " << minimizerSequence.size() << " " << graph.nodes_.size() << endl; //<< " " << graphPOA->_graph->_nodes.size() << endl;
				//	for(MinimizerType m : minimizerSequence){
				//		cout << "\t" << m << endl; 
				//	}
				
				//}

				
				int64_t nbMatches_forward = 0;
				int64_t alignScore_forward = 0;
				MinimizerAligner::Alignment alignment_forward;
				getAlignmentScore(referenceRead._minimizers, queryRead._minimizers, _minimizerAligner, alignment_forward, nbMatches_forward, alignScore_forward);
				
				int64_t nbMatches_reverse = 0;
				int64_t alignScore_reverse = 0;
				MinimizerAligner::Alignment alignment_reverse;
				getAlignmentScore(referenceRead._minimizers, queryRead_reverse, _minimizerAligner, alignment_reverse, nbMatches_reverse, alignScore_reverse);

				if(nbMatches_forward < 2 && nbMatches_reverse < 2) continue;


				if(nbMatches_forward > nbMatches_reverse){

					
					AlignmentResult2 al;

					for(size_t i=0; i<alignment_forward.size(); i++){
						al._alignments.push_back({alignment_forward[i].first, alignment_forward[i].second});
					}
					
					dbgGraph->addAlignment2(al, referenceRead._minimizers, queryRead._minimizers, queryRead._qualities);

					//extractAlignment(alignment_forward, readIndex, readMinimizers, contigIndex, subContigMinimizers_forward, alignScore_forward, false, startQueryPositionIndex);
				}
				else{

					vector<u_int8_t> quals = queryRead._qualities;
					std::reverse(quals.begin(), quals.end());

					AlignmentResult2 al;
					for(size_t i=0; i<alignment_reverse.size(); i++){
						al._alignments.push_back({alignment_reverse[i].first, alignment_reverse[i].second});
					}
					
					dbgGraph->addAlignment2(al, referenceRead._minimizers, queryRead_reverse, quals);

					//extractAlignment(alignment_reverse, readIndex, readMinimizers, contigIndex, subContigMinimizers_reverse, alignScore_reverse, true, startQueryPositionIndex);
				}

				

			}
			
			//dbgGraph->save("/pasteur/zeus/projets/p02/seqbio/gbenoit/workspace/home/tmp/poaGraph.gfa", referenceRead._originalMinimizers);
			//getchar();
		
			MinimizerRead correctedRead = computePath2(referenceRead, dbgGraph, referenceReadMinimizerPositionMapHighDensity);
			//dbgGraph->save("/pasteur/zeus/projets/p02/seqbio/gbenoit/workspace/home/tmp/poaGraph.gfa", referenceRead._originalMinimizers);
			//getchar();
			delete dbgGraph;
			


			const MinimizerRead& correctedReadTrimmed = trimCorrectedPath(referenceRead, correctedRead);
			
			return correctedReadTrimmed;
		}
		*/

		void getAlignmentScore(const vector<MinimizerType>& s1, const vector<MinimizerType>& s2, MinimizerAligner* minimizerAligner, MinimizerAligner::Alignment& alignmentResult, int64_t& nbMatches, int64_t& alignScore){
			

			MinimizerAligner::Alignment alignment = minimizerAligner->performAlignment(s1, s2); //al->Align(Utils::vec32_to_vec64(s2), s2.size(), graph);
			alignmentResult = alignment;

			nbMatches = 0;
			alignScore = 0;


			
			for (size_t i=0; i<alignment.size(); i++) {

				int64_t v1 = alignment[i].first;
				int64_t v2 = alignment[i].second;
				
				if(v1 == -1){ //insert in
					//alignScore -= 1;
				}
				else if(v2 == -1){ //insert in
					//alignScore -= 1;
				}
				else if(s1[v1] == s2[v2]){
					alignScore += 1;
					nbMatches += 1;
				}
				else{
					alignScore -= 1;
				}
			}


		}

		MinimizerRead computePath2(const MinimizerRead& read, Graph* graph, ReadMinimizerPositionMap& referenceReadMinimizerPositionMap){


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
				float minWeight = maxWeight*0.75;
				

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
					maxSuccessor = computeBestSuccessor(solidSuccessors, referenceReadMinimizerPositionMap);
				}

				if(maxSuccessor == nullptr){ //All successor completion score is 0
					break;
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

			if(_parent._print_debug){

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
			/*
			float supportMedian = Utils::compute_median(pathSupport);
			float minSupport = supportMedian * 0.25;

			int64_t trimStart = -1;
			int64_t trimEnd = -1;
			for(size_t i=0; i<pathSupport.size(); i++){

				u_int32_t support = pathSupport[i];

				if(support > minSupport){
					if(trimStart == -1) trimStart = i;
					trimEnd = i+1;
				}
			}

			//cout << trimStart << " " << trimEnd << endl;
			//for(size_t i=0; i<path.size(); i++){
			//	cout << i << " " << path[i] << endl;
			//}
			

			vector<MinimizerType> solidMinimizers;
			vector<u_int8_t> solidQualities;

			if(trimStart == -1){
				return read;
			}


			for(size_t i=trimStart; i<=trimEnd; i++){
				solidMinimizers.push_back(path[i]);
				solidQualities.push_back(pathQualities[i]);
			}
			

			//vector<MinimizerType> solidPath;

			//if(firstSolidPosition == -1){
			//	return solidPath;
			//}


			//for(size_t i=firstSolidPosition; i<=lastSolidPosition; i++){
			//	solidPath.push_back(path[i]);
			//}
			*/
			MinimizerRead correctedRead = {read._readIndex, path, {}, pathQualities, {}, read._readLength};
			return correctedRead;
			//return {read._readIndex, path, {}, pathQualities, read._readMinimizerDirections, read._originalMinimizers, read._originalQualities};
		}

		
		Edge* computeBestSuccessor(const vector<Edge*>& solidSuccessors, ReadMinimizerPositionMap& referenceReadMinimizerPositionMap){

			u_int64_t maxCompletion = 0;
			Edge* maxSuccessor = nullptr;

			for(Edge* successor : solidSuccessors){
				u_int64_t completion = computeSuccessorCompletion(successor, referenceReadMinimizerPositionMap);

				if(completion > maxCompletion){
					maxCompletion = completion;
					maxSuccessor = successor;
				}
			}

			return maxSuccessor;
		}

		u_int64_t computeSuccessorCompletion(Edge* successor, ReadMinimizerPositionMap& referenceReadMinimizerPositionMap){


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

					if(referenceReadMinimizerPositionMap.find(nn->_head->_minimizer) != referenceReadMinimizerPositionMap.end()){
						completion += nn->_weight;
					}

					queue.push_back(nn->_head);
				}


			}

			return completion;
		}

		
		MinimizerRead trimCorrectedPath(const MinimizerRead& originalRead, const MinimizerRead& correctedRead){

			/*
			cout << "trimming" << endl;

			cout << endl << "Original sequence:" << endl;
			for(MinimizerType m : read._originalMinimizers){
				cout << m << endl;
			}
			cout << endl << "Corrected sequence:" << endl;
			for(size_t i=0; i<correctedRead._minimizers.size(); i++){
				cout << i << ": " << correctedRead._minimizers[i] << endl;
			}
			*/
			
			/*
			spoa64::Graph graphAln{};

			vector<u_int64_t> weights(correctedRead._minimizers.size(), 1);
			graphAln.AddAlignment(spoa64::Alignment(), Utils::vec32_to_vec64(correctedRead._minimizers), correctedRead._minimizers.size(), weights);

			spoa64::Alignment alignment = al->Align(Utils::vec32_to_vec64(read._originalMinimizers), read._originalMinimizers.size(), graphAln);
			*/


			//cout << "a" << endl;
			MinimizerAligner::Alignment alignment;

			//if(minimizerAligner == nullptr){
			//	spoa64::Graph graphAln{};

			//	vector<u_int64_t> weights(correctedRead._minimizers.size(), 1);
			//	graphAln.AddAlignment(spoa64::Alignment(), Utils::vec32_to_vec64(correctedRead._minimizers), correctedRead._minimizers.size(), weights);

			//	alignment = al->Align(Utils::vec32_to_vec64(read._originalMinimizers), read._originalMinimizers.size(), graphAln);
			//}
			//else{
			alignment = _minimizerAligner->performAlignment(correctedRead._minimizers, originalRead._minimizers);
			//}
			//; = minimizerAligner->performAlignment(correctedRead._minimizers, read._originalMinimizers);
			//cout << "b" << endl;
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

					if(correctedRead._minimizers[v1] == originalRead._minimizers[v2]){

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

			//cout << startMatchPos << " " << endMatchPos << endl;
			//getchar();

			if(startMatchPos == -1 || startMatchPos == endMatchPos){
				return {correctedRead._readIndex, {}, {}, {}};
			}
		
			//cout << startMatchPos << " " << endMatchPos << endl;
			//if(_print_debug){
			//	cout << "\tTrim:" << startMatchPos << " " << endMatchPos << endl;
			//}

			vector<MinimizerType> minimizersTrimmed;
			//vector<u_int32_t> minimizersPosTrimmed;
			vector<u_int8_t> qualitiesTrimmed;

			for(size_t i=startMatchPos; i<endMatchPos; i++){
				minimizersTrimmed.push_back(correctedRead._minimizers[i]);
				qualitiesTrimmed.push_back(correctedRead._qualities[i]);
				//minimizersPosTrimmed.push_back(correctedRead._positions[i]);
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

			MinimizerRead correctedReadTrimmed = {correctedRead._readIndex, minimizersTrimmed, {}, qualitiesTrimmed, {}, correctedRead._readLength};

			return correctedReadTrimmed;
			//return {correctedRead._readIndex, minimizersTrimmed, {}, qualitiesTrimmed, correctedRead._readMinimizerDirections, correctedRead._originalMinimizers, correctedRead._originalQualities};
		}


	};








	/*
	AlignmentResult2 computeAlignment(const MinimizerRead& referenceRead, const MinimizerRead& queryRead, MinimizerAligner* minimizerAligner, bool isQueryReversed){



		vector<MinimizerType> queryMinimizers = queryRead._minimizers;
		if(isQueryReversed){
			std::reverse(queryMinimizers.begin(), queryMinimizers.end());
		}


		AlignmentResult2 alignmentResult = {0, 0, 0, 0, "", false, 0, 0, 0, 0, 0, 0};
		
		//spoa64::Graph graph{};
		//vector<u_int64_t> weights(referenceRead._minimizers.size(), 1);
		//graph.AddAlignment(spoa64::Alignment(), Utils::vec32_to_vec64(referenceRead._minimizers), referenceRead._minimizers.size(), weights);
		//spoa64::Alignment alignment = _alignmentEngine->Align(Utils::vec32_to_vec64(queryMinimizers), queryMinimizers.size(), graph);

		MinimizerAligner::Alignment alignment = minimizerAligner->performAlignment(referenceRead._minimizers, queryMinimizers);
		if(alignment.size() == 0) return alignmentResult;


	
		u_int32_t cigarReferenceStart = alignment[0].first;
		u_int32_t cigarQueryStart = alignment[0].second;
		string cigar = "";

		int64_t nbMatches = 0;
		int64_t nbMissmatches = 0;
		int64_t nbInsertions = 0;
		int64_t nbDeletions = 0;

		for (size_t i=0; i<alignment.size(); i++) {

			int64_t v1 = alignment[i].first;
			int64_t v2 = alignment[i].second;
			
			if(v1 == -1){ //insert in
				cigar += "I";
				nbInsertions += 1;
				if(_print_debug) cout << "\tInsertion: " << v2 << " " << queryMinimizers[v2] << "   " << queryRead._minimizersPos[v2] << endl;
			}
			else if(v2 == -1){ //insert in
				cigar += "D";
				nbDeletions += 1;
				if(_print_debug) cout << "\tDeletion : " << v1 << " " << referenceRead._minimizers[v1] << "    " << referenceRead._minimizersPos[v1] << endl;
			}
			else if(referenceRead._minimizers[v1] == queryMinimizers[v2]){
				cigar += "M";
				nbMatches += 1;
				if(_print_debug) cout << "\tMatch    : " << v1 << " " << v2 << " " << referenceRead._minimizers[v1] << " " << queryMinimizers[v2] << "    " << referenceRead._minimizersPos[v1] << " " << queryRead._minimizersPos[v2] << endl;
				//if(posStartBps == -1) posStartBps = queryMinimizersPos[v2];
				//posEndBps = queryMinimizersPos[v2];
			}
			else{
				cigar += "M";
				nbMissmatches += 1;
				if(_print_debug) cout << "\tMissmatch: " << v1 << " " << v2 << " " << referenceRead._minimizers[v1] << " " << queryMinimizers[v2] << endl;
				
			}
		}

		
		u_int64_t referenceSize = alignmentResult._nbMatches + alignmentResult._nbMissmatches + alignmentResult._nbDeletions;
		u_int64_t querySize = alignmentResult._nbMatches + alignmentResult._nbMissmatches + alignmentResult._nbInsertions;
		
		long double nbSeeds = min(referenceSize, querySize);
		float divergence = 0;

		//cout << nbMatches << " " << nbSeeds << endl;
		//double nbSeeds = alignmentResult._nbMatches + alignmentResult._nbMissmatches + alignmentResult._nbInsertions + alignmentResult._nbDeletions;

		if(alignmentResult._nbMatches == nbSeeds){
			divergence = 0;
		}
		else if(alignmentResult._nbMatches == 0){
			divergence = 1;
		}
		else{
			divergence = 1.0 - pow((alignmentResult._nbMatches / nbSeeds), 1.0/_minimizerSize);
		}
		
		

		string cigarCompressed = compressCigar(cigar);
		
		alignmentResult._referenceReadIndex = referenceRead._readIndex;
		alignmentResult._queryReadIndex = queryRead._readIndex;
		alignmentResult._cigarReferenceStart = cigarReferenceStart;
		alignmentResult._cigarQueryStart = cigarQueryStart;
		alignmentResult._cigar = cigarCompressed;
		alignmentResult._isQueryReversed = isQueryReversed;
		//alignmentResult._chainingScore = chain._chainingScore;
		alignmentResult._identity = 1.0-divergence;

		
		return alignmentResult;
	}
	*/

	/*
	string compressCigar(const string& cigar){

		string cigarCompressed = "";
		if(cigar.size() == 0) return cigarCompressed;

		char lastChar = '#';
		int nbOccurences = 1;
		//u_int64_t lastPos = 0;

		for(size_t i=0; i<cigar.size(); i++){
			//cout << i << " " << length << endl;
			char c = cigar[i];
			if(c == lastChar){
				nbOccurences += 1;
				continue;
			}
			if(lastChar != '#'){
				cigarCompressed += to_string(nbOccurences);
				cigarCompressed += lastChar;
				nbOccurences = 1;
				//rlePositions.push_back(lastPos);
				//lastPos = i;
				//cout << lastChar << endl;
			}
			lastChar = c;
		}

		cigarCompressed += to_string(nbOccurences);
		cigarCompressed += lastChar;
		//cout << lastChar << endl;
		//rleSequence += lastChar;
		//rlePositions.push_back(lastPos);
		//rlePositions.push_back(length);
		return cigarCompressed;
	}

	*/




	/*
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

	*/


	/*
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
	*/



	bool _print_debug;



	void writeRead(ReadType readIndex, const vector<MinimizerType>& minimizersToAdd, const vector<u_int8_t>& qualitiesToAdd){

		//cout << _kminmerSizeFirst << endl;
		if(minimizersToAdd.size() < _kminmerSizeFirst) return;
		//static MinimizerType maxHashValue = -1;
		//static MinimizerType minimizerBound = 0.005 * maxHashValue;

		//#pragma omp critical(dataupdate)
		#pragma omp critical
		{
			

			//vector<u_int64_t> validMinimizers;
			//for(u_int64_t m : minimizers){
			//	if(m > minimizerBound) continue;
			//	validMinimizers.push_back(m);
			//}

			
			//_readWriterQueue.push({readIndex, minimizers, qualities});

			//while(!_readWriterQueue.empty()){


			//	const ReadWriter& readWriter = _readWriterQueue.top();

				//cout << readWriter._readIndex << " " << _nextReadIndexWriter << " " << _readWriterQueue.size() << endl;

			//	if(readWriter._readIndex == _nextReadIndexWriter){

					//if(readWriter._minimizers.size() > 0){

						for(MinimizerType m : minimizersToAdd){
							_correctionCheckSum += m*minimizersToAdd.size();
						}

						
						vector<MinimizerType> minimizers = minimizersToAdd;
						vector<u_int32_t> minimizerPos(minimizers.size(), 0); 
						vector<u_int8_t> minimizerDirections(minimizers.size(), 0); 
						vector<u_int8_t> minimizerQualities = qualitiesToAdd; 
						vector<MinimizerType> minimizersFiltered;
						vector<u_int32_t> minimizerPosFiltered; 
						vector<u_int8_t> minimizerDirectionsFiltered; 
						vector<u_int8_t> minimizerQualitiesFiltered; 

						Utils::applyDensityThreshold(_minimizerDensity_assembly, minimizers, minimizerPos, minimizerDirections, minimizerQualities, minimizersFiltered, minimizerPosFiltered, minimizerDirectionsFiltered, minimizerQualitiesFiltered);
						
						minimizers = minimizersFiltered;
						minimizerPos = minimizerPosFiltered;
						minimizerDirections = minimizerDirectionsFiltered; 
						minimizerQualities = minimizerQualitiesFiltered; 
						
						//vector<MinimizerType> validMinimizers;

						//for(size_t i=0; i<readWriter._minimizers.size(); i++){
						//	u_int64_t m = readWriter._minimizers[i];
						//	if(m > minimizerBound) continue;

						//	validMinimizers.push_back(readWriter._minimizers[i]);
							//minimizerPosDummy.push_back(minimizers[i]);

						//}

						u_int32_t size = minimizers.size();
						_file_readData.write((const char*)&size, sizeof(size));

						u_int8_t isCircular = CONTIG_LINEAR;
						_file_readData.write((const char*)&isCircular, sizeof(isCircular));

						_file_readData.write((const char*)&minimizers[0], size*sizeof(MinimizerType));

						//vector<u_int32_t> minimizerPosDummy(size+1, 0);
						//vector<u_int8_t> minimizerDirectionDummy(size, 0);
						//vector<u_int8_t> qualitiesDummy(size, 0);
						//_file_readData.write((const char*)&minimizerPosDummy[0], (size+1)*sizeof(u_int32_t));
						//_file_readData.write((const char*)&minimizerDirectionDummy[0], size*sizeof(u_int8_t));
						//_file_readData.write((const char*)&qualitiesDummy[0], size*sizeof(u_int8_t));
					//}

					//if(readIndex % 10000 == 0){
					//	Logger::get().debug() << "\tCorrection checksum: " << readIndex << " " << _correctionCheckSum;
					//}
					//if(readWriter._readIndex > 1000) exit(1);
 					//_readWriterQueue.pop();
					//_nextReadIndexWriter += 1;
				//}
				//else{
				//	break;
				//}
			//}
			
		}

	}


	/*
	void evaluateCorrection(const CorrectedRead& correctedRead, const std::unique_ptr<spoa64::AlignmentEngine>& al){

		const vector<MinimizerType>& truthMinimizers = _mReads_evaluation[correctedRead._readIndex];


		u_int64_t nbMatches = 0;
		u_int64_t nbErrors = 0;

		computeAccuracy(correctedRead._minimizers, truthMinimizers, al, nbMatches, nbErrors);


		u_int64_t nbMatches_reverse = 0;
		u_int64_t nbErrors_reverse = 0;

		vector<MinimizerType> truthMinimizersReverse = truthMinimizers;
		std::reverse(truthMinimizersReverse.begin(), truthMinimizersReverse.end());

		computeAccuracy(correctedRead._minimizers, truthMinimizersReverse, al, nbMatches_reverse, nbErrors_reverse);

		if(nbMatches > nbMatches_reverse){

			_eval_nbBases += nbMatches + nbErrors;
			_eval_nbMatches += nbMatches;
		}
		else{
			_eval_nbBases += nbMatches_reverse + nbErrors_reverse;
			_eval_nbMatches += nbMatches_reverse;
		}

		cout << _eval_nbMatches << "\t" << _eval_nbBases << "\t" << (long double) _eval_nbMatches / _eval_nbBases << endl;
		
	}


	void computeAccuracy(const vector<MinimizerType>& readMinimizers, const vector<MinimizerType>& contigMinimizers, const std::unique_ptr<spoa64::AlignmentEngine>& al, u_int64_t& nbMatches, u_int64_t& nbErrors){
		
		nbMatches = 0;
		nbErrors = 0;

		bool print = false;

		if(print) cout << endl;
		vector<u_int64_t> weights(contigMinimizers.size(), 1);
		
		spoa64::Graph graph{};

		graph.AddAlignment(spoa64::Alignment(), Utils::vec32_to_vec64(contigMinimizers), contigMinimizers.size(), weights);
		spoa64::Alignment alignment = al->Align(Utils::vec32_to_vec64(readMinimizers), readMinimizers.size(), graph);



		for (const auto& it : alignment) {

			int64_t v1 = it.first;
			int64_t v2 = it.second;
			

			if(v1 == -1){ //insert in
				nbErrors += 1;
				if(print) cout << "\tInsertion: " << readMinimizers[v2] << endl;

			}
			else if(v2 == -1){ //insert in

				nbErrors += 1;
				if(print) cout << "\tDeletion : " << contigMinimizers[v1] << endl;
			}
			else if(contigMinimizers[v1] == readMinimizers[v2]){

				nbMatches += 1;
				if(print) cout << "\tMatch    : " << contigMinimizers[v1] << " " << readMinimizers[v2] << endl;

			}
			else{

				nbErrors += 1;
				if(print) cout << "\tMissmatch: " << contigMinimizers[v1] << " " << readMinimizers[v2] << endl;

			}
		}
		

	}



	void loadSimulatedReadMetadataFile(){

		std::string str;

		ifstream inputFile(_simulatedReadMetadataFilename);
		//std::getline(inputFile, str); //skip header

		while(!inputFile.eof()) {
			std::getline( inputFile, str);
			std::stringstream buffer(str);
			std::string temp;
			std::vector<string> fields;

			while( getline( buffer, temp, '\t') ) {
				fields.push_back(temp);
			}

			if(fields.size() == 0) break;

			//cout << str << endl;

			//cout << fields[2] << endl;
			//cout << fields.size() << endl;
			ReadType referenceIndex = stoull(fields[0]);
			ReadType readIndex = stoull(fields[1]);
			string strand = fields[2];
			u_int32_t refCoordStart = stoull(fields[3]);
			u_int32_t refCoordEnd = stoull(fields[4]);
			string identityStr = fields[5];
			float identity = 0;
			if(refCoordStart == -1) {
				
			}
			else{
				identity = stof(identityStr);
			}
			
			//cout << referenceIndex << endl;
			//cout << readIndex << endl;
			//cout << strand << endl;
			//cout << refCoordStart << endl;
			//cout << refCoordEnd << endl;
			//cout << identity << endl;

			bool isReversed = false;
			if(strand == "-") isReversed = true;
			
			SimulatedReadMetadata metadata = {referenceIndex, readIndex, isReversed, refCoordStart, refCoordEnd, identity};


			_simulatedReadMetadata.push_back(metadata);
			//_contigName_to_contigCoverage[contigName] = contigCoverage;
			//getchar();
			
		}

		inputFile.close();

	}

	*/


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

			unordered_map<ReadType, ReadMatch> bestMatches;

			for(const ReadMatch& match: matches){
				if(readIndex == match._readIndex) continue;

				if(bestMatches.find(match._readIndex) == bestMatches.end()){
					bestMatches[match._readIndex] = match;
				}
				else{
					if(match._alignLength > bestMatches[match._readIndex]._alignLength){
						bestMatches[match._readIndex] = match;
					}
				}
				//counts[match._readIndex] += 1;
				//if(readIndex < 1000000) _readIndex_to_bestReadMatch[readIndex].push_back(match);

				//ReadMatch match2 = match;
				//match2._readIndex = readIndex;
				//if(match._readIndex < 100000) _readIndex_to_bestReadMatch[match._readIndex].push_back(match2);
			}

			for(const auto& it: bestMatches){
				_readIndex_to_bestReadMatch[readIndex].push_back(it.second);
			}

			//_readIndex_to_bestReadMatch[readIndex] = matches;
		}

		file.close();

	}


	void evalMapping(){

		if(!fs::exists(_inputDir + "/readIndex_to_readIndex.bin")) return;

		MinimizerChainer* minimizerChainer = new MinimizerChainer(_minimizerSize);

		cout << "load read mapping" << endl;
		loadReadToReadMapping();

		cout << "load reads" << endl;
		loadMinimizerReads(_inputDir + "/read_data_init.txt", false, false);

		cout << "eval" << endl;

		ofstream resultFile("/pasteur/zeus/projets/p02/seqbio/gbenoit/workspace/home/scripts/nanopore/evalMapping/data_evalMapping.tsv");
		resultFile << "WfmashDivergence\tNanoMdbgDivergence" << endl;

		for(const auto& it : _readIndex_to_bestReadMatch){

			u_int32_t referenceReadIndex = it.first;

			for(const ReadMatch& match: it.second){

				u_int32_t queryReadIndex = match._readIndex;
				const MinimizerRead referenceRead = _mReads[referenceReadIndex];
				const MinimizerRead queryRead = _mReads[queryReadIndex];

				ReadMinimizerPositionMap referenceReadMinimizerPositionMap;


				for(u_int32_t i=0; i<referenceRead._minimizers.size(); i++){

					MinimizerType minimizer = referenceRead._minimizers[i];
					u_int32_t position = referenceRead._minimizersPos[i];
					bool isReversed = referenceRead._readMinimizerDirections[i];

					MinimizerPosition minmizerPosition = {position, i, isReversed};
					referenceReadMinimizerPositionMap[minimizer].push_back(minmizerPosition);
				}

				vector<Anchor> anchors;
				//_minimizerChainer->_anchorIndex = 0;

				for(size_t i=0; i<queryRead._minimizers.size(); i++){

					MinimizerType queryMinimizer = queryRead._minimizers[i];

					if(referenceReadMinimizerPositionMap.find(queryMinimizer) == referenceReadMinimizerPositionMap.end()) continue;

					u_int32_t queryPosition = queryRead._minimizersPos[i];
					bool queryIsReversed = queryRead._readMinimizerDirections[i];

					const vector<MinimizerPosition>& referenceMinimizerPositions = referenceReadMinimizerPositionMap[queryMinimizer];
					
					for(const MinimizerPosition& referenceMinimizerPosition : referenceMinimizerPositions){
						anchors.push_back({referenceMinimizerPosition._position, queryPosition, referenceMinimizerPosition._isReversed != queryIsReversed, referenceMinimizerPosition._positionIndex, i});
					}




				}

				//cout << endl;
				//cout << "----" << endl;
				//cout << anchors.size() << endl;
				AlignmentResult2 chainingAlignment = minimizerChainer->computeChainingAlignment(anchors, referenceRead, queryRead, nullptr, _maxChainingBand_highDensity);
				

				//cout << match._isReversed << "\t" << chainingAlignment._isQueryReversed << endl;
				//cout << match._referenceStart << "\t" << match._referenceEnd << "\t" << match._queryStart << "\t" << match._queryEnd << endl;
				//cout << chainingAlignment._referenceStart << "\t" << chainingAlignment._referenceEnd << "\t" << chainingAlignment._queryStart << "\t" << chainingAlignment._queryEnd << endl;

				//cout << match._hangLeft << "\t" << match._hangRight << endl;
				//cout << chainingAlignment._overHangStart << "\t" << chainingAlignment._overHangEnd << endl;
				//cout << match._divergence << "\t" << chainingAlignment._divergence << endl;

				if(chainingAlignment._alignLength < 1000) continue;

				float alignIdentity = 1 - (match._alignNbMatches / (long double)(match._queryEnd - match._queryStart));
				//resultFile << match._divergence << "\t" << chainingAlignment._divergence << endl;
				//getchar();
				//if(chainingAlignment._divergence == 0) getchar();
			}

		}

		resultFile.close();
	}

};	


#endif 


