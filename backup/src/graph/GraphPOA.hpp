

#ifndef MDBG_GRAPHPOA
#define MDBG_GRAPHPOA

#include "Commons.hpp"

namespace GraphPOA{
    

/*
class Scoring{

    public:

    int _gap_open;
    int _gap_extend;
    //u_int32_t _match_fn;
    u_int32_t _match_scores;
    int _xclip_prefix;
    int _xclip_suffix;
    int _yclip_prefix;
    int _yclip_suffix;

    Scoring(int gap_open, int gap_extend){
        _gap_open = gap_open;
        _gap_extend = gap_extend;
    }

    int match_fn(u_int32_t a, u_int32_t b){
        
        if(a == b){
            return 1;
        }
        else{
            return -1;
        }

    }

};
*/



class Matrix{

    public:

    const static int NULL_VALUE = -1;

    u_int32_t _nbRows;
    u_int32_t _nbCols;

    vector<int> _elements;

    Matrix(){

    }

    void init(u_int32_t nbRows ,u_int32_t nbCols, int initialValue){

        _elements.clear();

        _nbRows = nbRows;
        _nbCols = nbCols;

        u_int64_t nbElements = _nbRows*_nbCols;
        _elements.resize(nbElements, initialValue);
        
    }

    void set(u_int32_t row, u_int32_t col, int value){
        u_int32_t index = col + row*_nbCols;
        _elements[index] = value;
    }

    int get(u_int32_t row, u_int32_t col){
        u_int32_t index = col + row*_nbCols;
        return _elements[index];
    }

    void getMaxCoord(int& maxRow, int& maxCol){

        maxRow = -1;
        maxCol = -1;

        int maxValue = std::numeric_limits<int>::min();

        for(u_int32_t row=0; row<_nbRows; row++){
            for(u_int32_t col=0; col<_nbCols; col++){


                int value = get(row, col);

                if(value == maxValue){
                    if(col > maxCol){
                        maxRow = row;
                        maxCol = col;
                    }
                }
                else if(value > maxValue){
                    maxValue = value;
                    maxRow = row;
                    maxCol = col;

                }

            }
        }
    }

    void save(const string& outputFilename){

        ofstream matrixFile(outputFilename);


        for(size_t row=0; row<_nbRows; row++){

            for(size_t col=0; col<_nbCols; col++){

                matrixFile << get(row, col) << "\t";
            }

            matrixFile << endl;
        }

        matrixFile.close();

    }

    void print(){

        for(size_t row=0; row<_nbRows; row++){

            cout << row << "\t";
            
            for(size_t col=0; col<_nbCols; col++){

                cout << get(row, col) << "\t";
            }

            cout << endl;
        }

        cout << endl;

    }

};

class Edge;

class Node{
    public:

    u_int32_t _nodeIndex;
    vector<Edge*> _successors;
    vector<Node*> _predecessors;
    u_int64_t _unitigIndex;
    u_int32_t _weight;
    vector<Node*> _alignedTo;

    Node(u_int32_t nodeIndex, u_int64_t unitigIndex){
        _nodeIndex = nodeIndex;
        _unitigIndex = unitigIndex;
        _weight = 0;
    }

    
    string toString(){
        return to_string(_unitigIndex) + "(" + to_string(_nodeIndex) + ")";
    }




    //void addSuccessor(Node* node){
    //    _successors.push_back(node);
    //    node->_predecessors.push_back(this);
    //}
};


class Edge{

    public:

    Node* _fromNode;
    Node* _toNode;
    u_int32_t _weight;

    Edge(Node* fromNode, Node* toNode, u_int32_t weight){
        _fromNode = fromNode;
        _toNode = toNode;
        _weight = weight;
    }
};

class Graph{
    public:


    bool _needTopoSort;
    unordered_map<u_int32_t, Node*> _nodes;
    u_int32_t _nodeIndex;
    vector<Node*> _topoSortNodes;

    Graph(){
        _nodeIndex = 0;
        _needTopoSort = true;
    }

    ~Graph(){
        for(auto& it : _nodes){

            for(Edge* edge : it.second->_successors){
                delete edge;
            }
            delete it.second;
        }
        //for(size_t i=0; i<_nodes.size(); i++){
        //    delete _nodes[i];
        //}
    }

    /*
    vector<vector<u_int32_t>> getSolidConnectedComponents(){


        vector<vector<u_int32_t>> components;
        unordered_set<u_int32_t> isVisited;

        for(Node* node : _nodes){

            if(isVisited.find(node->_nodeIndex) != isVisited.end()) continue;
            isVisited.insert(node->_nodeIndex);

            queue<u_int32_t> queue;
            queue.push(node->_nodeIndex);

            vector<u_int32_t> component = {node->_nodeIndex};

            while (!queue.empty()){

                u_int32_t currentNodeIndex = queue.front();

                queue.pop();

                //unordered_set<Node*> neighbors;
                //for(Node* nodeIndexNN : currentUnitig->_successors){
                //    neighbors.insert(nodeIndexNN);
                //    neighbors.insert(unitigIndex_toReverseDirection(nodeIndexNN));
                //}
                //for(UnitigGraph::Node* nodeIndexNN : currentUnitig->_predecessors){
                //    neighbors.insert(nodeIndexNN);
                //    neighbors.insert(unitigIndex_toReverseDirection(nodeIndexNN));
                //}

                //vector<UnitigGraph::Node*> neighbors = get_neighbours(currentNodeIndex);


                //for(UnitigGraph::Node* nodeIndexNN : neighbors){
                for(u_int32_t nodeIndexSuccessor : getSolidSuccessors(currentNodeIndex)){
                    
                    Node* nodeSuccessor = _nodes[nodeIndexSuccessor];

                    if (isVisited.find(nodeSuccessor->_nodeIndex) != isVisited.end()) continue;

                    isVisited.insert(nodeSuccessor->_nodeIndex);
                    component.push_back(nodeSuccessor->_nodeIndex);

                    queue.push(nodeSuccessor->_nodeIndex);

                }


            }

            components.push_back(component);

        }
        
        return components;
    }


    */

    Node* addNode(u_int64_t unitigIndex){

        //cout << "Add Node: " << unitigIndex << endl;

        Node* node = new Node(_nodeIndex, unitigIndex);
        _nodes[_nodeIndex] = node;
        //_nodes.push_back(node);
        _nodeIndex += 1;

        _needTopoSort = true;

        return node;
    }

    void addEdge(Node* fromNode, Node* toNode){


        if(fromNode == nullptr) return;
        if(fromNode == toNode) return;

        for(Edge* edge : fromNode->_successors){
            if(edge->_toNode == toNode){
                //edge->_weight += 1;
                return;
            }
        }

        Edge* edge = new Edge(fromNode, toNode, 1);
        fromNode->_successors.push_back(edge);
        toNode->_predecessors.push_back(fromNode);

        _needTopoSort = true;

        //cout << "Add edge: " << fromNode->_nodeIndex << " -> " << toNode->_nodeIndex<< endl;
    }

    //void removeEdge(Node* fromNode, Node* toNode){


    //    fromNode->_successors.erase(std::remove(fromNode->_successors.begin(), fromNode->_successors.end(), 8), fromNode->_successors.end());
    //    toNode->_predecessors.erase(std::remove(toNode->_predecessors.begin(), toNode->_predecessors.end(), 8), toNode->_predecessors.end());

    //}

    Edge* getEdge(Node* fromNode, Node* toNode){

        if(fromNode == nullptr) return nullptr;

        for(Edge* edge : fromNode->_successors){
            if(edge->_toNode == toNode){
                return edge;
            }
        }

        return nullptr;
    }

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
    
    /*
    void removeDuplicatedEdges(){
        //return;
        
        while(true){

            bool isDuplication = false;

            //cout << "a" << endl;
            for(Node* node : _nodes){

                unordered_map<u_int64_t, Edge*> unitigIndexes;

                //cout << "b" << endl;
                for(Edge* edge : node->_successors){
                    
                    if(edge->_weight == 0) continue; //Duplicated removed
                    if(edge->_toNode == nullptr) continue;
                    if(edge == nullptr) continue;
                    //cout << edge << endl;
                    //cout << edge->_toNode << endl;
                   // cout << edge->_toNode->_unitigIndex << endl;
                    u_int64_t unitigIndex = edge->_toNode->_unitigIndex;

                    //cout << "d" << endl;
                    if(unitigIndexes.find(unitigIndex) == unitigIndexes.end()){
                        unitigIndexes[unitigIndex] = edge;
                    }
                    else{

                        Edge* edgeKeep = unitigIndexes[unitigIndex];
                        Edge* edgeDuplicated = edge;

                        edgeKeep->_weight += edgeDuplicated->_weight;
                        edgeDuplicated->_weight = 0;
                        isDuplication = true;



                        for(Edge* edge : edgeDuplicated->_toNode->_successors){
                            //cout << (edge == nullptr) << endl;
                            edgeKeep->_toNode->_successors.push_back(edge);
                        }

                        edgeDuplicated->_toNode->_successors.clear();

                        //cout << "h" << endl;

                    }
                    //if(edge->_weight > 5)
                    //unitigIndexes[edge->_toNode->_unitigIndex] += 1;
                }



                

            }

            if(!isDuplication) break;
        }
        */

        /*
        for(Node* node : _nodes){
            
            unordered_map<u_int32_t, u_int32_t> unitigIndexesCount;

            for(Edge* edge : node->_successors){
                if(edge->_weight > 0)
                unitigIndexesCount[edge->_toNode->_unitigIndex] += 1;
            }

            
            for(const auto& it : unitigIndexesCount){
                if(it.second > 1){
                    cout << "duplicated edges: utg" << node->_unitigIndex << " utg" << it.first << endl;
                    fs::remove("/home/gats/workspace/tmp/poagraph.txt");
                    save("/home/gats/workspace/tmp/poagraph.txt");
                    getchar();
                }
            }


        }
        */
    //}

    /*
    struct Pseudonode{
        int _pnode_id;
        vector<int> _predecessors;
        vector<int> _successors;
        vector<int> _node_ids;
    };

    void simplified_graph_rep(vector<Pseudonode>& pseudonodes){

        pseudonodes.clear();

        unordered_map<int, int> node_to_pn;
        unordered_map<int, vector<int>> pn_to_nodes;

        int cur_pnid = 0;

        for(Node* node : _nodes){

            if(node_to_pn.find(node->_nodeIndex) == node_to_pn.end()){

                vector<int> node_ids;
                node_ids.push_back(node->_nodeIndex);
                for(Node* n : node->_alignedTo){
                    node_ids.push_back(n->_nodeIndex);
                }

                pn_to_nodes[cur_pnid] = node_ids;

                for(int n : node_ids){
                    node_to_pn[n] = cur_pnid;
                }

                cur_pnid += 1;
            }
        }

        for(int pnid=0; pnid<cur_pnid; pnid++){

            const vector<int>& nids = pn_to_nodes[pnid];
            vector<int> preds;
            vector<int> succs;

            for(int nodeIndex : nids){
                for(Node* nn : _nodes[nodeIndex]->_predecessors){
                    preds.push_back(node_to_pn[nn->_nodeIndex]);
                }
                for(Edge* edge : _nodes[nodeIndex]->_successors){
                    succs.push_back(node_to_pn[edge->_toNode->_nodeIndex]);
                }
            }

            pseudonodes.push_back({pnid, preds, succs, nids});
        }
    }
    
    void computeTopologicalSort(){

        if(!_needTopoSort) return;

        _needTopoSort = false;
        //cout << "---" << endl;
        _topoSortNodes.clear();

        vector<int> sortedlist;
        unordered_set<int> completed;

        vector<Pseudonode> pseudonodes;
        simplified_graph_rep(pseudonodes);

        //for(Pseudonode p: pseudonodes){
        //    cout << p._pnode_id << " " << p._successors.size() << " " << p._predecessors.size() << " "<< p._node_ids.size() << endl;
        //}
        while(sortedlist.size() < _nodes.size()){

            int found = -1;

            for(size_t pnid=0; pnid<pseudonodes.size(); pnid++){
                if(completed.find(pnid) == completed.end() && pseudonodes[pnid]._predecessors.size() == 0){
                    found = pnid;
                    break;
                }
            }

            dfs(found, completed, sortedlist, pseudonodes);
        }

        for(int id : sortedlist){
            _topoSortNodes.push_back(_nodes[id]);
        }

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
    

    void dfs(int start, unordered_set<int>& complete, vector<int>& sortedlist, const vector<Pseudonode>& pseudonodes){

        //cout << start << endl;

        stack<int> stack;
        stack.push(start);
        unordered_set<int> started;

        while(!stack.empty()){
            int pnodeID = stack.top();
            //cout << "\t" << pnodeID << endl;
            stack.pop();

            if(complete.find(pnodeID) != complete.end()) continue;

            if(started.find(pnodeID) != started.end()){
                complete.insert(pnodeID);

                for(int nid : pseudonodes[pnodeID]._node_ids){
                    sortedlist.insert(sortedlist.begin(), nid);
                }

                started.erase(pnodeID);
                continue;
            }

            const vector<int> successors = pseudonodes[pnodeID]._successors;
            started.insert(pnodeID);
            stack.push(pnodeID);
            //cout << "\t\tPush: " << pnodeID << endl;
            for(int succ : successors){
                stack.push(succ);
                //cout << "\t\tPush: " << succ << endl;
            }
            

        }

    }
    */

    u_int32_t getMaxInWeight(Node* n){

        u_int32_t maxWeight = 0;

        for(Node* node : n->_predecessors){
            for(Edge* edge : node->_successors){
                if(edge->_toNode->_nodeIndex == n->_nodeIndex){
                    float weight = edge->_weight;
                    if(weight > maxWeight){
                        maxWeight = weight;
                    }
                    break;
                }
            }
        }

        return maxWeight;
    }

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

    void filterGraph(){

        unordered_set<Node*> solidNodes;
        //for(Node* node : _nodes){
        for(const auto& it : _nodes){
            Node* node = it.second;
            for(Edge* edge : node->_successors){

                //if(getMaxInWeight(edge->_toNode) < 5) continue;
                if(edge->_weight <= 1) continue;
                solidNodes.insert(edge->_fromNode);
                solidNodes.insert(edge->_toNode);
                //outputFile << node->_unitigIndex << "(" << node->_nodeIndex  << ")" << "\t" << edge->_toNode->_unitigIndex << "(" << edge->_toNode->_nodeIndex  << ")" << "\t" << edge->_weight << "    " << getMaxInWeight(edge->_toNode) << endl;
            }
        }



        /*
        //for(Node* node : _nodes){
        for(const auto& it : _nodes){
            Node* node = it.second;
            for(Edge* edge : node->_successors){
                if(solidNodes.find(edge->_toNode) == solidNodes.end() ||  solidNodes.find(edge->_fromNode) == solidNodes.end()){
                    edgeToRemove.insert(edge);
                }
            }
        }
        */

        unordered_set<Node*> nodeToRemove;

        //for(Node* node : _nodes){
        for(const auto& it : _nodes){
            Node* node = it.second;
            if(solidNodes.find(node) != solidNodes.end()) continue;
            nodeToRemove.insert(node);
        }

        /*
        for(Edge* edge : edgeToRemove){

            Node* fromNode = edge->_fromNode;
            Node* toNode = edge->_toNode;

            
            //for(Edge* e : fromNode->_successors){

            //}
            fromNode->_successors.erase(std::remove(fromNode->_successors.begin(), fromNode->_successors.end(), edge), fromNode->_successors.end());
            

            for(Node* node : toNode->_predecessors){
                //edge->_toNode->_predecessors.erase(std::remove(edge->_toNode->_predecessors.begin(), edge->_toNode->_predecessors.end(), node), edge->_toNode->_predecessors.end());
                
            }

            //toNode->_predecessors.erase(std::remove(toNode->_predecessors.begin(), toNode->_predecessors.end(), 8), toNode->_predecessors.end());

            delete edge;
        }
        */
        
        
        //cout << "Size: " << _nodes.size() << endl;
        for(Node* node : nodeToRemove){

            //cout << "erase: " << node->_nodeIndex << endl;

            for(Edge* edge : node->_successors){
                edge->_toNode->_predecessors.erase(std::remove(edge->_toNode->_predecessors.begin(), edge->_toNode->_predecessors.end(), node), edge->_toNode->_predecessors.end());
                delete edge;
            }



            for(Node* n : node->_predecessors){

                Edge* edgeToRemove = nullptr;

                for(Edge* edge : n->_successors){
                    if(edge->_toNode == node){
                        edgeToRemove = edge;
                        break;
                    }
                }

                n->_successors.erase(std::remove(n->_successors.begin(), n->_successors.end(), edgeToRemove), n->_successors.end());
                delete edgeToRemove;
            }

            u_int32_t nodeIndex = node->_nodeIndex;

            delete _nodes[node->_nodeIndex];
            _nodes.erase(nodeIndex);

        }

        //cout << "Size: " << _nodes.size() << endl;
        
        //cout << "Key present: " << endl;
        //for(const auto& it : _nodes){
        //    cout << "\t" << it.first << " " << it.second->_nodeIndex << endl;
        //}
    }
};

/*
class CompactGraph{

public: 




    class UnitigEdge;

    class UnitigNode{
        public:

        u_int32_t _nodeIndex;
        vector<UnitigEdge*> _successors;
        vector<UnitigNode*> _predecessors;
        vector<Node*> _nodes;

        UnitigNode(u_int32_t nodeIndex, vector<Node*>& nodes){
            _nodeIndex = nodeIndex;
            _nodes = nodes;
        }

        string toString(){
            
            string s = "";
            for(Node* node : _nodes){
                s += to_string(node->_unitigIndex) + "-";
            }
            //s += "(" + to_string(weight) + ")";
            s.pop_back();
            return s;
            //return to_string(_unitigIndex) + "(" + to_string(_nodeIndex) + ")";
        }

    };


    class UnitigEdge{

        public:

        UnitigNode* _toNode;
        u_int32_t _weight;

        UnitigEdge(UnitigNode* toNode, u_int32_t weight){
            _toNode = toNode;
            _weight = weight;
        }
    };

    u_int32_t _nodeIndex;
    vector<UnitigNode*> _nodes;
    unordered_map<Node*, UnitigNode*> _nodeToUnitigNode;

    CompactGraph(Graph* graph){
        _nodeIndex = 0;

        initNodes(graph);
        initEdges();
    }

    ~CompactGraph(){
        for(size_t i=0; i<_nodes.size(); i++){
            delete _nodes[i];
        }
    }


    void initNodes(Graph* graph){

        //UnitigNode* prevUnitigNode = nullptr;

        vector<Node*> queue;
        queue.push_back(graph->_nodes[0]);
        unordered_set<Node*> isVisited;

        while(!queue.empty()){

            Node* currentNode = queue[queue.size()-1];
            queue.pop_back();

            if(isVisited.find(currentNode) != isVisited.end()) continue;
            isVisited.insert(currentNode);

            vector<Node*> nodes;
            nodes.push_back(currentNode);

            while(true){

                vector<Edge*>& successors = currentNode->_successors;

                if(successors.size() != 1) break;

                Node* successorNode = successors[0]->_toNode;
                vector<Node*>& predecessors = successorNode->_predecessors;
                
                if(predecessors.size() != 1) break;

                nodes.push_back(successorNode);
                currentNode = successorNode;

            }

            UnitigNode* unitigNode = addNode(nodes);
            cout << "Add unitig: ";
            for(Node * node: nodes){
                cout << node->_nodeIndex << " "; 
            }
            cout << endl;

            for(Node* node : nodes){
                _nodeToUnitigNode[node] = unitigNode;
            }

            //if(prevUnitigNode != nullptr){
            //    addEdge(prevUnitigNode, unitigNode, weight);
            //}

            //cout << "\tUnitig: " << endl;
            //for(Node* node : nodes){
            //    cout << "\t" << node->_unitigIndex << endl;
            //}

            for(Edge* successor : currentNode->_successors){
                queue.push_back(successor->_toNode);
            }


        }

    }

    void initEdges(){

        for(size_t i=0; i<_nodes.size(); i++){
            
            UnitigNode* node = _nodes[i];

            Node* endNode = node->_nodes[node->_nodes.size()-1];
            
            for(Edge* edge : endNode->_successors){
                addEdge(node, _nodeToUnitigNode[edge->_toNode], edge->_weight);
                cout << "Add edge: " << endNode->_nodeIndex << " -> " << edge->_toNode->_nodeIndex << endl;
            }

            //Node* startNode = node->_nodes[0];

            //for(Node* node : startNode->_predecessors){
                //u_int32_t toUnitigIndex = _unitigGraph->_nodeIndex_to_unitigIndex[nodeIndex];
                //_unitigGraph->addPredecessors(node->_unitigIndex, toUnitigIndex);
                //addEdge(node, _nodeToUnitigNode[edge->_toNode], edge->_weight);
            //}

        }


    }



    void detectBubbles(unordered_set<u_int32_t>& bubbleUnitigIndexes){

        bubbleUnitigIndexes.clear();

        simplify();

        for(UnitigNode* node : _nodes){

            UnitigNode* bubbleExit = isSuperbubble(node);
            if(bubbleExit != nullptr){

                vector<UnitigNode*> mostSupportedNodes;
                getSuperbubbleMostSupportedPath(node, bubbleExit, mostSupportedNodes);

                vector<UnitigNode*> bubbleNodes;
                collectSuperbubbleNodes(node, bubbleExit, bubbleNodes);

                unordered_set<UnitigNode*> mostSupportedNodesIndex;
                for(UnitigNode* node : mostSupportedNodes){
                    mostSupportedNodesIndex.insert(node);
                }

                for(UnitigNode* node : bubbleNodes){
                    if(mostSupportedNodesIndex.find(node) != mostSupportedNodesIndex.end()) continue;

                    for(Node* internalNode : node->_nodes){
                        bubbleUnitigIndexes.insert(internalNode->_unitigIndex);
                    }

                }

            }
        }

        printDebug();
    }

    void simplify(){

        while(true){

            bool isModification = false;
            

            bool isTipRemoved = removeTips();
            if(isTipRemoved){
                isModification = true;
            }


            if(!isModification) break;
        }

    }

    bool removeTips(){

        bool isModificationAll = false;

        while(true){

            bool isModification = false;

            for(UnitigNode* node : _nodes){

                if(isTip(node)){
                    //cout << node->toString() << " " << isTip(node) << endl;
                    isModification = true;
                    isModificationAll = true;
                    removeNode(node);
                }

            }

            if(!isModification) break;
        }

        if(isModificationAll){

            for(UnitigNode* node : _nodes){
                if(node->_nodeIndex == -1) continue;
                recompact(node);
            }
        }

        return isModificationAll;

    }

    bool isTip(const UnitigNode* node){

        if(node->_nodeIndex == -1) return false;
        
        if(node->_predecessors.size() == 0) return false; //root node

        if(node->_successors.size() > 0) return false;
        if(node->_predecessors.size() != 1) return false;
        //if(node->_predecessors.size() == 0) return false;

        return true;
    }


    UnitigNode* isSuperbubble(UnitigNode* nodeSource){

        unordered_set<UnitigNode*> isVisited;
        unordered_set<UnitigNode*> seen;
        vector<UnitigNode*> queue;

        queue.push_back(nodeSource);

        while(queue.size() > 0){

            UnitigNode* vNode = queue[queue.size()-1];
            //u_int32_t v = vNode->_unitigIndex;

            queue.pop_back();

            isVisited.insert(vNode);
            if(seen.find(vNode) != seen.end()){
                seen.erase(vNode);
            }

            if(vNode->_successors.size() == 0){
                return nullptr; //abort tip
            }

            for(UnitigEdge* uEdge : vNode->_successors){

                UnitigNode* u = uEdge->_toNode;
                
                if(u == nodeSource){
                    return nullptr; //cycle including s
                }

                if(isVisited.find(u) == isVisited.end()){
                    seen.insert(u);      
                }
                else{

                    return nullptr; //Cycle within superbubble
                }

            }


            for(UnitigEdge* uEdge : vNode->_successors){

                UnitigNode* u = uEdge->_toNode;
                
                bool allPredecessorsAreVisited = true;

                for(UnitigNode* pNode : u->_predecessors){
                    if(isVisited.find(pNode) == isVisited.end()){
                        allPredecessorsAreVisited = false;
                        break;
                    }
                }

                if(allPredecessorsAreVisited){
                    queue.push_back(u);
                }

                if(queue.size() == 1 && seen.size() == 1 && seen.find(queue[0]) != seen.end()){ //only one vertex t is left in S and no other vertex is seen 
                    
                    UnitigNode* t = queue[0];

                    return t;
                    //u_int32_t t = tNode->_unitigIndex;

                    //vector<UnitigGraph::Node*> tNodeSuccessors;
                    //getValidSuccessors(tNode, isUnitigIndexValid, tNodeSuccessors);

                    //for(UnitigEdge* edgeSuccessor : t->_successors){
                    //    if(edgeSuccessor->_toNode == nodeSource)
                    //}
                    //if(std::find(t->_successors.begin(), t->_successors.end(), nodeSource) == t->_successors.end()){
                    //    return t;
                    //}
                    //else{
                    //    return nullptr; // cycle including s
                    //}
                }
            }
        }

        return nullptr;

    }

    void collectSuperbubbleNodes(UnitigNode* sourceNode, UnitigNode* exitNode, vector<UnitigNode*>& nodes){

        nodes.clear();

        //unordered_set<UnitigNode*> isVisited;
        vector<UnitigNode*> queue;

        //isVisited.insert(sourceNode);
        //isVisited.insert(exitNode);
        queue.push_back(sourceNode);

        while (!queue.empty()){

            UnitigNode* currentNode = queue[queue.size()-1];
            queue.pop_back();

            for(UnitigEdge* edgeSuccessor : currentNode->_successors){

                if(edgeSuccessor->_toNode == exitNode) continue;
                //if (isVisited.find(nn->_unitigIndex) != isVisited.end()) continue;

                queue.push_back(edgeSuccessor->_toNode);
                //isVisited.insert(nn);
                nodes.push_back(edgeSuccessor->_toNode);

            }
        }

    }
        

    void getSuperbubbleMostSupportedPath(UnitigNode* sourceNode, UnitigNode* exitNode, vector<UnitigNode*>& mostSupportedNodes){

        mostSupportedNodes.clear();

        UnitigNode* nodeCurrent = sourceNode;
        
        while(true){

            vector<UnitigNode*> maxAbNodes;
            float maxAbundance = 0;
            u_int32_t maxV = -1;

            //bool isExitAllowed = true;
            //if(nodeCurrent->_successors.size() > 1) isExitAllowed = false; //Do not allow entering exit if there is still node to visit

            for(UnitigEdge* edgeSuccessor : nodeCurrent->_successors){
                
                //if(!isExitAllowed && nn == exitNode) continue;

                if(edgeSuccessor->_weight == maxAbundance){
                    maxAbNodes.push_back(edgeSuccessor->_toNode);
                }
                else if(edgeSuccessor->_weight > maxAbundance){
                    maxAbundance = edgeSuccessor->_weight;
                    maxAbNodes.clear();
                    maxAbNodes.push_back(edgeSuccessor->_toNode);
                }
            }

            if(maxAbNodes.size() == 1){
                nodeCurrent = maxAbNodes[0];
            }
            else{

                
                //vector<BubbleSide> bubbleSides;
                //for(UnitigGraph::Node* nn : maxAbNodes){
                //    bubbleSides.push_back({nn, nn->_nodes.size(), nn->startNode()}); //_unitigs[unitigIndex]._quality
                //}


                //std::sort(bubbleSides.begin(), bubbleSides.end(), BubbleSideComparator);

                nodeCurrent = maxAbNodes[0];

            }
            

            if(nodeCurrent == exitNode) break;

            mostSupportedNodes.push_back(nodeCurrent);
        }


        
    }
    
    void removeNode(UnitigNode* node){

        for(UnitigEdge* edge : node->_successors){
            edge->_toNode->_predecessors.erase(std::remove(edge->_toNode->_predecessors.begin(), edge->_toNode->_predecessors.end(), node), edge->_toNode->_predecessors.end());
        }

        for(UnitigNode* nn : node->_predecessors){

            UnitigEdge* edgeToRemove = nullptr;

            for(UnitigEdge* edge : nn->_successors){
                if(edge->_toNode == node){
                    edgeToRemove = edge;
                    break;
                }
            }

            if(edgeToRemove != nullptr){
                nn->_successors.erase(std::remove(nn->_successors.begin(), nn->_successors.end(), edgeToRemove), nn->_successors.end());
            }
        }

        node->_nodeIndex = -1;

    }

    void recompact(UnitigNode* node){

        while(true){

            bool isChanged = false;

            if(node->_successors.size() == 1){

                UnitigEdge* edgeSuccessor = node->_successors[0];

                if(edgeSuccessor->_toNode->_predecessors.size() == 1){

                    mergeNode(node, edgeSuccessor->_toNode);

                }
            }

            if(!isChanged) break;
        }
    }

    
    void mergeNode(UnitigNode* node1, UnitigNode* node2){


        node1->_nodes.insert( node1->_nodes.end(), node2->_nodes.begin(), node2->_nodes.end());
        
        node1->_successors = node2->_successors;

        for(UnitigEdge* edge : node2->_successors){

            for(size_t i=0; i<edge->_toNode->_predecessors.size(); i++){
                if(edge->_toNode->_predecessors[i] == node2){
                    edge->_toNode->_predecessors[i] = node1;
                    break;
                }
            }



        }

        node2->_nodeIndex = -1;

    }


    UnitigNode* addNode(vector<Node*>& nodes){

        UnitigNode* node = new UnitigNode(_nodeIndex, nodes);
        _nodes.push_back(node);
        _nodeIndex += 1;

        return node;
    }

    void addEdge(UnitigNode* fromNode, UnitigNode* toNode, int weight){

        UnitigEdge* edge = new UnitigEdge(toNode, weight);
        fromNode->_successors.push_back(edge);
        toNode->_predecessors.push_back(fromNode);

    }

    void save(const string& outputFilename){


        cerr << "\tSaving graph" << endl;


        
        ofstream outputFile(outputFilename);


        for(UnitigNode* node : _nodes){
            
            if(node->_nodeIndex == -1) continue;

            outputFile << "S" << "\t" << node->toString() << "\t" << "*" << "\t" << "LN:i:" << node->_nodes.size()*100 << "\t" << "dp:i:" << "0" << endl;

            for(UnitigEdge* edge : node->_successors){
                //outputFile << "L" << "\t" << node->toString() << "\t" << "+" << "\t" << edge->_toNode->toString() << "\t" << "+" << "\t" << "50" << "M" << endl;     
            }


        }

        outputFile.close();
        
        cerr << "\tdone" << endl;

    }

    void printDebug(){

        for(UnitigNode* node : _nodes){
            
            if(node->_nodeIndex == -1) continue;

            cout << node->toString() << endl;
            for(UnitigEdge* edge : node->_successors){
                cout << "\t" << edge->_toNode->toString() << " (" << edge->_weight << ")" << endl;
            }
        }


    }
};
*/
/*
enum AlignmentOperationType {
    Match,
    Del,
    Ins,
};

struct AlignmentOperation{
    
    AlignmentOperationType _type;
    u_int32_t _prevRow;
    u_int32_t _prevCol;

};

class Alignment {

    public:

    int _score;
    u_int32_t _yStart;
    vector<AlignmentOperation> _operations;

    Alignment(int score, uint32_t yStart, const vector<AlignmentOperation>& operations){
        _score = score;
        _yStart = yStart;
        _operations = operations;
    }
};

class TracebackMatrix{

    public:

    const static u_int32_t CELL_NO_PREV = -1;
    const static int CELL_MIN_SCORE = -858993459;

    u_int32_t _nbRows;
    u_int32_t _nbCols;
    Node* _lastNode;

    struct Cell{
        int _score;
        AlignmentOperation _op;
    };

    vector<Cell> _elements;

    /// * `m` - the number of nodes in the DAG
    /// * `n` - the length of the query sequence
    TracebackMatrix(u_int32_t m ,u_int32_t n){

        _nbRows = m;
        _nbCols = n;

        u_int64_t nbElements = (m+1)*(n+1);
        _elements.resize(nbElements, {0, {AlignmentOperationType::Match, CELL_NO_PREV, CELL_NO_PREV}});
        
        _lastNode = nullptr;
    }

    void set(u_int32_t row, u_int32_t col, const Cell& cell){
        u_int32_t index = col + row*(_nbCols+1);
        _elements[index] = cell;
    }

    Cell get(u_int32_t row, u_int32_t col){
        u_int32_t index = col + row*(_nbCols+1);
        return _elements[index];
    }
    
    Cell max(Cell a, Cell b){

        if(a._score == b._score){
            return a;
        }
        else if(a._score > b._score){
            return a;
        }
        else{
            return b;
        }
    }

};
*/

class GraphPOA{

    public:

    Graph* _graph;
    //int _gap_open;
    //int _gap_extend;
    //u_int32_t _match_fn;
    //u_int32_t _match_scores;
    //int _xclip_prefix;
    //int _xclip_suffix;
    //int _yclip_prefix;
    //int _yclip_suffix;

    //const static int MATCH_SCORE = 1;
    //const static int MISMATCH_SCORE = -1;
    //const static int GAP_SCORE = -1;
    //int NONE;
    int _matchScore;
    int _mismatchScore;
    int _gapScore;
    int _none;

    Node* _lastSequenceNode;

    GraphPOA(const vector<u_int64_t>& sequence){


        _print_debug = false;

        _time_computeTopologicalSort = 0;
        _time_initializeDynamicProgrammingData = 0;
        _time_alignSequenceToGraph = 0;
        _time_backtrack = 0;
        _time_addAlignmentToGraph = 0;
        _time_total = 0;

        _matchScore = 1;
        _mismatchScore = -1;
        _gapScore = -1;
        _none = std::numeric_limits<int>::min();

        if(_print_debug){
            cout << "Nb sequences: " << sequence.size() << endl;
            cout << "Init with sequence: ";
            for(u_int64_t unitigIndex : sequence){
                cout << unitigIndex << " ";
            }
            cout << endl;
        }

        //_gap_open = -1;
        //_gap_extend = -1;
        /*
        _graph = new Graph();
        Node* node1 = _graph->addNode();
        Node* node2 = _graph->addNode();
        Node* node3 = _graph->addNode();
        Node* node4 = _graph->addNode();
        Node* node5 = _graph->addNode();

        node5->addSuccessor(node4);
        node5->addSuccessor(node3);
        node4->addSuccessor(node2);
        node2->addSuccessor(node1);
        node3->addSuccessor(node1);

        
        vector<Node*> topoSortNodes;
        _graph->computeTopologicalSort(topoSortNodes);


        for(Node* node : topoSortNodes){
            cout << node->_nodeIndex << endl;
        }
        */


        _graph = new Graph();

        Node* prevNode = _graph->addNode(sequence[0]);
        
        for(size_t i=1; i<sequence.size(); i++){

            Node* currentNode = _graph->addNode(sequence[i]);
            _graph->addEdge(prevNode, currentNode);

            prevNode = currentNode;
            _lastSequenceNode = currentNode;
        }

    }


    ~GraphPOA(){
        delete _graph;
    }
    /*
    void addSequence(const vector<u_int32_t>& sequence){

        cout << "\n\n\nAdd sequence: ";
        for(u_int32_t unitigIndex : sequence){
            cout << unitigIndex << " ";
        }
        cout << endl;

        TracebackMatrix* traceback = alignSequenceToGraph(sequence);
        Alignment* alignment = getAlignmentOperations(traceback);

        addAlignmentToGraph(sequence, alignment);

        delete alignment;
        delete traceback;

    }


    TracebackMatrix* alignSequenceToGraph(const vector<u_int32_t>& sequence){


        u_int32_t nbNodesGraph = _graph->_nodes.size();
        u_int32_t sequenceLength = sequence.size();

        TracebackMatrix* traceback = new TracebackMatrix(nbNodesGraph, sequenceLength);
         
        
        for(int row=0; row < traceback->_nbRows+1; row++){
            traceback->set(row, 0, {0, {AlignmentOperationType::Del, TracebackMatrix::CELL_NO_PREV, TracebackMatrix::CELL_NO_PREV}});
        }
        for(int col=0; col < traceback->_nbCols+1; col++){
            traceback->set(0, col, {col*_gap_open, {AlignmentOperationType::Ins, TracebackMatrix::CELL_NO_PREV, TracebackMatrix::CELL_NO_PREV}});
        }


        traceback->set(0, 0, {0, {AlignmentOperationType::Match, TracebackMatrix::CELL_NO_PREV, TracebackMatrix::CELL_NO_PREV}});

        // construct the score matrix (O(n^2) space)
        //let mut topo = Topo::new(&self.graph);

        vector<Node*> topoSortNodes;
        _graph->computeTopologicalSort(topoSortNodes);

        //writeMatrix(sequence, traceback, topoSortNodes);
        //getchar();

        for(Node* node : topoSortNodes){


            //cout << endl << "\t" << node->_unitigIndex << "(" << node->_nodeIndex << ")" << endl;

            u_int32_t referenceUnitigIndex = node->_unitigIndex;
            u_int32_t i = node->_nodeIndex + 1;

            traceback->_lastNode = node;
            vector<Node*> predecessors = node->_predecessors;

            for(u_int32_t j_p=0; j_p<sequence.size(); j_p++){
                
                u_int32_t sequenceUnitigIndex = sequence[j_p];
                u_int32_t j = j_p + 1;
                TracebackMatrix::Cell maxCell;

                if(predecessors.size() == 0){
                    maxCell._score = traceback->get(0, j - 1)._score + match_fn(referenceUnitigIndex, sequenceUnitigIndex);
                    maxCell._op._type = AlignmentOperationType::Match;
                    maxCell._op._prevRow = TracebackMatrix::CELL_NO_PREV;
                    maxCell._op._prevCol = TracebackMatrix::CELL_NO_PREV;
                    cout << "\tNo pred score: " << maxCell._score << endl;
                }
                else{
                    maxCell._score = TracebackMatrix::CELL_MIN_SCORE;
                    maxCell._op._type = AlignmentOperationType::Match;
                    maxCell._op._prevRow = TracebackMatrix::CELL_NO_PREV;
                    maxCell._op._prevCol = TracebackMatrix::CELL_NO_PREV;

                    for(Node* prevNode : predecessors){

                        u_int32_t i_p = prevNode->_nodeIndex + 1;
                        int gapPenalty = determineGapPenalty(traceback->get(i_p, j), {AlignmentOperationType::Del, i_p-1, i});


                        TracebackMatrix::Cell cellMatch;
                        cellMatch._score = traceback->get(i_p, j - 1)._score + match_fn(referenceUnitigIndex, sequenceUnitigIndex);
                        cellMatch._op._type = AlignmentOperationType::Match;
                        cellMatch._op._prevRow = i_p - 1;
                        cellMatch._op._prevCol = i - 1;
                        
                        TracebackMatrix::Cell cellDel;
                        cellDel._score = traceback->get(i_p, j)._score + gapPenalty;
                        cellDel._op._type = AlignmentOperationType::Del;
                        cellDel._op._prevRow = i_p - 1;
                        cellDel._op._prevCol = i;

                        TracebackMatrix::Cell maxPrevCell = traceback->max(cellMatch, cellDel);

                        maxCell = traceback->max(maxCell, maxPrevCell);

                    }
                }

                int gapPenalty = determineGapPenalty(traceback->get(i, j-1), {AlignmentOperationType::Ins, i - 1, TracebackMatrix::CELL_NO_PREV});
                

                TracebackMatrix::Cell cellIns;
                cellIns._score = traceback->get(i, j - 1)._score + gapPenalty;
                cellIns._op._type = AlignmentOperationType::Ins;
                cellIns._op._prevRow = i - 1;
                cellIns._op._prevCol = TracebackMatrix::CELL_NO_PREV;

                TracebackMatrix::Cell cellScore;
                cellScore = traceback->max(maxCell, cellIns);

                traceback->set(i, j, cellScore);

                writeMatrix(sequence, traceback, topoSortNodes);
                cout << endl << "\t\t" << i << " " << j << " " << cellScore._score << " " << cellScore._op._type << " " << cellScore._op._prevRow << " " << cellScore._op._prevCol << endl;
                getchar();

            }  

        }



        return traceback;
    }

    void writeMatrix(const vector<u_int32_t>& sequence, TracebackMatrix* traceback, const vector<Node*>& topoSortNodes){

        ofstream matrixFile("/home/gats/workspace/tmp/alignMatrix.tsv");

        matrixFile << "\t";
        for(u_int32_t unitigIndex : sequence){
            matrixFile << unitigIndex << "\t";
        }
        matrixFile << endl;

        for(size_t row=0; row<traceback->_nbRows+1; row++){

            if(row < topoSortNodes.size()){
                matrixFile << topoSortNodes[row]->_unitigIndex << "(" << topoSortNodes[row]->_nodeIndex << ")" << "\t";
            }
            else{
                matrixFile << "\t";
            }

            for(size_t col=0; col<traceback->_nbCols+1; col++){

                matrixFile << traceback->get(row, col)._score << "\t";
            }

            matrixFile << endl;
        }

        matrixFile.close();

    }

    int match_fn(u_int32_t a, u_int32_t b){
        
        if(a == b){
            return 1;
        }
        else{
            return -1;
        }

    }

    int determineGapPenalty(const TracebackMatrix::Cell& cell, AlignmentOperation current_op){
        return _gap_extend;
    }

    struct MaxCandidate{
        u_int32_t _i;
        u_int32_t _j;
        int _score;
    };

    static bool MaxCandidateComparator(const MaxCandidate& a, const MaxCandidate& b){
        return a._score > b._score;
    }

    Alignment* getAlignmentOperations(TracebackMatrix* tracebackMatrix) {
        
        Alignment* alignment = nullptr;
        
        vector<Node*> startingNodes;
        for(Node* node : _graph->_nodes){
            //if(node->_successors.size() == 0){
                startingNodes.push_back(node);
            //}
        }

        //for(size_t i=0; i<tracebackMatrix->_nbRows; i++){
        //    cout << tracebackMatrix->get(i, tracebackMatrix->_nbCols-1)._score << endl;
        //}
        vector<MaxCandidate> maxCandidate;
        for(Node* node : startingNodes){
            //cout << node->_unitigIndex << " " << tracebackMatrix->get(node->_nodeIndex+1, tracebackMatrix->_nbCols-1)._score << endl;
            maxCandidate.push_back({node->_nodeIndex+1, tracebackMatrix->_nbCols, tracebackMatrix->get(node->_nodeIndex+1, tracebackMatrix->_nbCols)._score});
        }

        std::sort(maxCandidate.begin(), maxCandidate.end(), MaxCandidateComparator);

        int maxScore = maxCandidate[0]._score;

        for(const MaxCandidate& candidate: maxCandidate){
        
        cout << "\t" << candidate._score << endl;
        //cout <<  candidate._score << " " << maxScore << endl;;
        //if(candidate._score < maxScore) break;


        string debugAlignStr = "";
        vector<AlignmentOperation> ops;

        u_int32_t i0 = candidate._i; //maxCandidate[0]._i;
        u_int32_t j0 = candidate._j; //maxCandidate[0]._j;


        u_int32_t i = i0;
        u_int32_t j = j0;

        while(i > 0 && j > 0){
            
            //cout << "\t" << i << " " << j << endl;

            AlignmentOperation op = tracebackMatrix->get(i, j)._op;
            ops.push_back(op);

            if(op._prevRow == TracebackMatrix::CELL_NO_PREV && op._prevCol == TracebackMatrix::CELL_NO_PREV){
                if(op._type == AlignmentOperationType::Match){
                    debugAlignStr += "m";
                    //if debug{print!("m");}
                    j -= 1; //?
                    break;
                }
                else if(op._type == AlignmentOperationType::Del){
                    debugAlignStr += "d";
                    //if debug{print!("d");}
                    break;
                }
                else if(op._type == AlignmentOperationType::Ins){
                    debugAlignStr += "i";
                    //if debug{print!("i");}
                    i -= 1;
                    j -= 1; //?
                }
            }
            else{
                if(op._type == AlignmentOperationType::Match){
                    debugAlignStr += "M";
                    //if debug{print!("M");}
                    i = op._prevRow + 1;
                    j -= 1;
                }
                else if(op._type == AlignmentOperationType::Del){
                    debugAlignStr += "-";
                    //if debug{print!("-");}
                    i = op._prevRow + 1;
                }
                else if(op._type == AlignmentOperationType::Ins){
                    debugAlignStr += "I";
                    //if debug{print!("I");}
                    i = op._prevRow + 1;
                    j -= 1; 
                }
            }

        }  
        

        std::reverse(ops.begin(), ops.end());


        std::reverse(debugAlignStr.begin(), debugAlignStr.end());
        cout << debugAlignStr << endl;

        if(alignment == nullptr){
        alignment = new Alignment(tracebackMatrix->get(i0, j0)._score, j, ops);
        }
        }


        return alignment;
    }


    void addAlignmentToGraph(const vector<u_int32_t>& sequence, const Alignment* alignment){

        cout << "Add alignment" << endl;

        Node* prevNode = _graph->_nodes[0];
        u_int32_t prev_i = 0;
        u_int32_t i = alignment->_yStart;

        for(const AlignmentOperation op : alignment->_operations){

            //cout << i << " " << op._type << " " << op._prevRow << " " << op._prevCol << endl;

            if(op._type == AlignmentOperationType::Match){
                if(op._prevRow == TracebackMatrix::CELL_NO_PREV && op._prevCol == TracebackMatrix::CELL_NO_PREV){
                    cout << "\tMatch skipped" << endl;
                    i += 1;
                }
                else{

                    Node* node = _graph->_nodes[op._prevCol];

                    if(sequence[i] == node->_unitigIndex){

                        cout << "\tMatch" << endl;
                        Edge* edge = _graph->getEdge(prevNode , node);

                        //cout << "\t" << prevNode->_unitigIndex << " " << node->_unitigIndex << " " << (edge != nullptr) << endl;

                        if(edge == nullptr){
                            _graph->addEdge(prevNode, node);
                        }
                        else{
                            edge->_weight += 1;
                        }

                        prevNode = node;
                        prev_i = i;
                    }
                    else{

                        cout << "\tMissmatch" << endl;

                        Node* node = _graph->addNode(sequence[i]);
                        _graph->addEdge(prevNode, node);

                        prevNode = node;
                        prev_i = i;
                    }

                    i += 1;
                }

            }
            else if(op._type == AlignmentOperationType::Ins){
                if(op._prevRow == TracebackMatrix::CELL_NO_PREV && op._prevCol == TracebackMatrix::CELL_NO_PREV){
                    cout << "\tInsert skipped" << endl;
                    i += 1;
                }
                else{

                    cout << "\tInsert" << endl;
                    Node* node = _graph->addNode(sequence[i]);
                    _graph->addEdge(prevNode, node);

                    prevNode = node;
                    prev_i = i;

                    i += 1;
                }
            }
            else if(op._type == AlignmentOperationType::Del){
            }


        }

    }
    */









    enum AlignmentOperationType {
        Match,
        Del,
        Ins,
    };

    struct AlignmentOperation{
        int _score;
        int _prevRow;
        int _prevCol;
        AlignmentOperationType _operationType;
    };

    static bool AlignmentOperationComparator(const AlignmentOperation& a, const AlignmentOperation& b){
        return a._score > b._score;
    }


    int matchscore(u_int64_t c1, u_int64_t c2){
        if(c1 == c2)
            return _matchScore;
        else
            return _mismatchScore;
    }

    float _time_computeTopologicalSort;
    float _time_initializeDynamicProgrammingData;
    float _time_alignSequenceToGraph;
    float _time_backtrack;
    float _time_addAlignmentToGraph;
    float _time_total;

    bool _print_debug;
    void addSequence(const vector<u_int64_t>& sequence){


        if(_print_debug){
            cout << "\n\n\nAdd sequence: ";
            for(u_int64_t unitigIndex : sequence){
                cout << unitigIndex << " ";
            }
            cout << endl;
        }

        if(_print_debug) cout << "topo sort" << endl;
        
        auto startTotal = high_resolution_clock::now();

        auto start = high_resolution_clock::now();
        _graph->computeTopologicalSort();
        auto stop = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(stop - start);
        _time_computeTopologicalSort += duration.count();

        if(_print_debug) cout << "init dynamic programming" << endl;
        start = high_resolution_clock::now();
        initializeDynamicProgrammingData(sequence);
        stop = high_resolution_clock::now();
		duration = duration_cast<microseconds>(stop - start);
        _time_initializeDynamicProgrammingData += duration.count();

        if(_print_debug) cout << "align sequence to graph" << endl;
        start = high_resolution_clock::now();
        alignSequenceToGraph(sequence);
        stop = high_resolution_clock::now();
		duration = duration_cast<microseconds>(stop - start);
        _time_alignSequenceToGraph += duration.count();

        //_scores.save("/home/gats/workspace/tmp/alignMatrix.tsv");
        
        if(_print_debug) cout << "backtrack" << endl;
        start = high_resolution_clock::now();
        backtrack(sequence);
        stop = high_resolution_clock::now();
		duration = duration_cast<microseconds>(stop - start);
        _time_backtrack += duration.count();

        if(_print_debug) cout << "add alignment to graph" << endl;
        //cout << printAlignmentStrings(_matches, _strindexes, sequence) << endl;
        start = high_resolution_clock::now();
        addAlignmentToGraph(sequence);
        stop = high_resolution_clock::now();
		duration = duration_cast<microseconds>(stop - start);
        _time_addAlignmentToGraph += duration.count();


        stop = high_resolution_clock::now();
		duration = duration_cast<microseconds>(stop - startTotal);
        _time_total += duration.count();
        //for(Node* node : _graph->_topoSortNodes){
        //    cout << node->_unitigIndex << endl;
        //}
        //cout << "done" << endl;
    }

    int _lastj;
    void alignSequenceToGraph(const vector<u_int64_t>& sequence){

        for(int i=0; i<_graph->_topoSortNodes.size(); i++){

            Node* node = _graph->_topoSortNodes[i];
            u_int64_t pbase = node->_unitigIndex;
            
            if(node == _lastSequenceNode) _lastj = i;

            for(int j=0; j<sequence.size(); j++){

                //cout << i << " " << j << endl;
                u_int64_t sbase = sequence[j];

                vector<int> prevIndices;
                collectPrevIndices(node, prevIndices);

                vector<AlignmentOperation> candidates;
                candidates.push_back({_scores.get(i+1, j) + _gapScore, i+1, j, AlignmentOperationType::Ins});

                for(int predIndex : prevIndices){
                    candidates.push_back({_scores.get(predIndex+1, j+1) + _gapScore, predIndex+1, j+1, AlignmentOperationType::Del});
                    candidates.push_back({_scores.get(predIndex+1, j) + matchscore(sbase, pbase), predIndex+1, j, AlignmentOperationType::Match});
                }

                //cout << candidates.size() << endl;
                std::sort(candidates.begin(), candidates.end(), AlignmentOperationComparator);

                const AlignmentOperation& maxCandidate = candidates[0];
                _scores.set(i+1, j+1, maxCandidate._score);
                _backGrphIdx.set(i+1, j+1, maxCandidate._prevRow);
                _backStrIdx.set(i+1, j+1, maxCandidate._prevCol);
                AlignmentOperationType movetype = maxCandidate._operationType;

                //if(_scores.get(i+1, j+1) < 0){
                //    _scores.set(i+1, j+1, 0);
                //    _backGrphIdx.set(i+1, j+1, -1);
                //    _backStrIdx.set(i+1, j+1, -1);
                //}

            }
        }

    }

    void insertions(int i, int l2, vector<bool>& inserted){

        for(size_t i=0; i<inserted.size(); i++){
            inserted[i] = false;
        }

        for(int j=0; j<l2; j++){

            int insscore = _scores.get(i+1, j) + _gapScore;

            if(insscore >= _scores.get(i+1, j+1)){
                _scores.set(i+1, j+1, insscore);
                inserted[j] = true;
            }
        }

    }


    void matchscoreVec(u_int64_t c, const vector<u_int64_t>& sequence, vector<int>& scores){

        scores.clear();

        for(size_t i=0; i<sequence.size(); i++){
            if(sequence[i] == c){
                scores.push_back(_matchScore);
            }
            else{
                scores.push_back(_mismatchScore);
            }
        }
    }

	template<typename T>
    void printVector(const vector<T>& vec){
        for(T val : vec){
            cout << val << " ";
        }
        cout << endl;
    }

	template<typename T>
    void maximum(vector<T>& vec1, const vector<T>& vec2){
        for(size_t i=0; i<vec1.size(); i++){
            if(vec1[i] > vec2[i]){

            }
            else{
                vec1[i] = vec2[i];
            }
        }
    }

    void alignSequenceToGraphFast(const vector<u_int64_t>& sequence){
        
        int l2 = sequence.size();
        //seqvec = numpy.array(list(self.sequence))


        vector<int> arrange0;
        for(size_t i=0; i<l2; i++){
            arrange0.push_back(i);
        }

        vector<int> arrange1;
        for(size_t i=1; i<l2+1; i++){
            arrange1.push_back(i);
        }


        vector<bool> inserted(l2, false);

        for(int i=0; i<_graph->_topoSortNodes.size(); i++){

            Node* node = _graph->_topoSortNodes[i];
            u_int64_t gbase = node->_unitigIndex;
            
            vector<int> predecessors;
            collectPrevIndices(node, predecessors);

            vector<int> deletescore;
            for(size_t j=1; j<_scores._nbCols; j++){
                deletescore.push_back(_scores.get(predecessors[0]+1, j) + _gapScore);
            }

            vector<int> bestdelete(l2, predecessors[0]+1);

            //printVector<int>(deletescore);
            //printVector<int>(bestdelete);

            vector<int> matchpoints;
            matchscoreVec(gbase, sequence, matchpoints);

            vector<int> matchscore;
            for(size_t j=0; j<_scores._nbCols-1; j++){
                matchscore.push_back(_scores.get(predecessors[0]+1, j) + matchpoints[j]);
            }

            vector<int> bestmatch(l2, predecessors[0]+1);


            //printVector<int>(matchpoints);
            //printVector<int>(matchscore);
            //printVector<int>(bestmatch);

            for(size_t i=1; i<predecessors.size(); i++){

                int predecessor = predecessors[i];

                vector<int> newdeletescore;
                for(size_t j=1; j<_scores._nbCols; j++){
                    newdeletescore.push_back(_scores.get(predecessor+1, j) + _gapScore);
                }

                for(size_t i=0; i<newdeletescore.size(); i++){
                    if(newdeletescore[i] > deletescore[i]){
                        bestdelete[i] = predecessor+1;
                    }
                    //else{
                        //bestdelete[i] 
                    //}
                }

                maximum(deletescore, newdeletescore);

                gbase = _graph->_nodes[predecessor]->_unitigIndex;
                matchscoreVec(gbase, sequence, matchpoints);

                vector<int> newmatchscore;
                for(size_t j=0; j<_scores._nbCols-1; j++){
                    newmatchscore.push_back(_scores.get(predecessor+1, j) + matchpoints[j]);
                }

                for(size_t i=0; i<newmatchscore.size(); i++){
                    if(newmatchscore[i] > matchscore[i]){
                        bestmatch[i] = predecessor+1;
                    }
                }

                maximum(matchscore, newmatchscore);

            }

            vector<bool> deleted;
            for(size_t i=0; i<matchscore.size(); i++){
                deleted.push_back(deletescore[i] >= matchscore[i]);
            }

            for(size_t j=0; j<deleted.size(); j++){
                if(deleted[j]){
                    _backGrphIdx.set(i+1, j+1, bestdelete[j]);
                }
                else{
                    _backGrphIdx.set(i+1, j+1, bestmatch[j]);
                }
            }

            for(size_t j=0; j<deleted.size(); j++){
                if(deleted[j]){
                    _backStrIdx.set(i+1, j+1, arrange1[j]);
                }
                else{
                    _backStrIdx.set(i+1, j+1, arrange0[j]);
                }
            }


            for(size_t j=0; j<deleted.size(); j++){
                if(deleted[j]){
                    _scores.set(i+1, j+1, deletescore[j]);
                }
                else{
                    _scores.set(i+1, j+1, matchscore[j]);
                }
            }

            //printVector<bool>(deleted);
            //_backGrphIdx.print();
            //_backStrIdx.print();
            //_scores.print();


            insertions(i, l2, inserted);

            for(size_t j=0; j<inserted.size(); j++){
                if(inserted[j]){
                    _backGrphIdx.set(i+1, j+1, i+1);
                }
                //else{
                //    _backGrphIdx.set(i+1, j+1, _backGrphIdx.get(i+1, j+1));
                //}
            }

            for(size_t j=0; j<inserted.size(); j++){
                if(inserted[j]){
                    _backStrIdx.set(i+1, j+1, arrange0[j]);
                }
                //else{
                //    _backStrIdx.set(i+1, j+1, matchscore[j]);
                //}
            }

            //getchar();
        }
    }

    /*
            # choose best options available of match, delete
            deleted       = deletescore >= matchscore
            backGrphIdx[i+1, 1:] = numpy.where(deleted, bestdelete, bestmatch)
            backStrIdx [i+1, 1:] = numpy.where(deleted, numpy.arange(1, l2+1), numpy.arange(0, l2))
            scores[i+1, 1:] = numpy.where(deleted, deletescore, matchscore)

            # insertions: updated in place, don't depend on predecessors
            insertions(i, l2, scores, inserted)
            backGrphIdx[i+1, 1:] = numpy.where(inserted, i+1, backGrphIdx[i+1, 1:])
            backStrIdx[i+1, 1:] = numpy.where(inserted, numpy.arange(l2), backStrIdx[i+1, 1:])


        return self.backtrack(scores, backStrIdx, backGrphIdx, nodeIndexToID)
    */

    void initializeDynamicProgrammingData(const vector<u_int64_t>& sequence){


        u_int64_t l1 = _graph->_nodes.size();
        u_int64_t l2 = sequence.size();


        _scores.init(l1+1, l2+1, 0);
        _backStrIdx.init(l1+1, l2+1, 0);
        _backGrphIdx.init(l1+1, l2+1, 0);



        for(size_t i=0; i<_scores._nbCols; i++){
            _scores.set(0, i, _gapScore*i);
        }

        
        for(int i=0; i<_graph->_topoSortNodes.size(); i++){
            Node* node = _graph->_topoSortNodes[i];

            vector<int> prevIdxs;
            collectPrevIndices(node, prevIdxs);

            int best = _scores.get(prevIdxs[0]+1, 0);
            for(int prevIdx : prevIdxs){
                best = max(best, _scores.get(prevIdx+1, 0));
            }
            //_scores.set(i+1, 0, best + _gapScore);
        }
        


    }

    void collectPrevIndices(Node* node, vector<int>& indices){

        indices.clear();

        for(Node* predNode : node->_predecessors){
            //if(_nodeIDtoIndex.find(predNode->_nodeIndex) == _nodeIDtoIndex.end()){
            //    cout << "LALA" << endl;
            //    getchar();
            //}
            indices.push_back(_graph->_nodeIDtoIndex[predNode->_nodeIndex]);
        }

        if(indices.size() == 0){
            indices.push_back(-1);
        }
    }

    Matrix _scores;
    Matrix _backStrIdx;
    Matrix _backGrphIdx;
    vector<int> _matches;
    vector<int> _strindexes;

    void backtrack(const vector<u_int64_t>& sequence){

        //_scores.print();

        _matches.clear();
        _strindexes.clear();

        int besti = -1;
        int bestj = -1;


        //_scores.getMaxCoord(besti, bestj);

        
        //besti = _scores._nbRows-1;
        besti = 1 ;//_lastj+1;
        bestj = _scores._nbCols-1;
        

        
        int bestscore = _scores.get(besti, bestj);

        for(int i=1; i<_graph->_topoSortNodes.size(); i++){

            //Node* node = _graph->_topoSortNodes[index];

            int score = _scores.get(i+1, bestj);
            if(score > bestscore){
                bestscore = score;
                besti = i+1;
            }
            

        }
        
        //cout << _graph->_nodes.size() << endl;
        //cout << besti << " " << bestj << endl;
        //cout << endl << endl << endl;
        //getchar();


        /*
        int maxValue = std::numeric_limits<int>::min();

        for(int row=0; row<_scores._nbRows; row++){
        //for(int row=_scores._nbRows-1; row>=0; row--){

            int value = _scores.get(row, bestj);

            if(value > maxValue){
                maxValue = value;
                besti = row;
                //bestj = bestj;
            }
        }
        
        vector<u_int32_t> maxCandidates;

        for(int row=0; row<_scores._nbRows; row++){

            int value = _scores.get(row, bestj);

            if(value == maxValue){
                maxCandidates.push_back(row);
            }
        }

        cout << "Best score: " << maxValue << endl;
        cout << "Start coord: " << besti << " " << bestj << endl;   
        */

        /*
        vector<int> terminalIndices;
        int maxScore = std::numeric_limits<int>::min();

        for(int i=0; i<_graph->_topoSortNodes.size(); i++){

            Node* node = _graph->_topoSortNodes[i];
            if(node->_successors.size() == 0){
                terminalIndices.push_back(i);
            }
        }


        besti = terminalIndices[0] + 1;
        int bestscore = _scores.get(besti, bestj);

        for(size_t ii=1; ii<terminalIndices.size(); ii++){
            int i = terminalIndices[ii];

            int score = _scores.get(i+1, bestj);

            cout << i+1 << " " << score << endl;
            if (score > bestscore){
                bestscore = score;
                besti = i+1;
            }
        }

        cout << "Best score: " << bestscore << endl;
        for(size_t ii=0; ii<terminalIndices.size(); ii++){
            int i = terminalIndices[ii];

            int score = _scores.get(i+1, bestj);
            if (score == bestscore){
                cout << "Possible start coord: " << besti << " " << bestj << endl;    
            }
        }

        cout << "Start coord: " << besti << " " << bestj << endl;    
        */
        //while(_scores.get(besti, bestj) > 0){
        //int maxAln = 0;

        //cout << "Max candidates: " << maxCandidates.size() << endl;

        //for(u_int32_t row : maxCandidates){

        //    cout << row << endl;
            //besti = row;
            //bestj = _scores._nbCols-1;

            //cout << "haha: " << besti << " " << bestj << endl;
            vector<int> matches;
            vector<int> strindexes;

            while(besti != 0 && bestj != 0){

                //cout << besti << " " << bestj << endl;

                int nexti = _backGrphIdx.get(besti, bestj);
                int nextj = _backStrIdx.get(besti, bestj);
                int curstridx = bestj-1;
                int curnodeidx = _graph->_nodeIndexToID[besti-1];

                if(nextj != bestj){
                    strindexes.push_back(curstridx);
                }
                else{
                    strindexes.push_back(_none);
                }

                if(nexti != besti){
                    matches.push_back(curnodeidx);
                }
                else{
                    matches.push_back(_none);
                }

                //cout << besti << " " << bestj << " " << curstridx << " " << curnodeidx << endl;

                besti = nexti;
                bestj = nextj;
            }

            std::reverse(matches.begin(), matches.end());
            std::reverse(strindexes.begin(), strindexes.end());

            //string aln = printAlignmentStrings(matches, strindexes, sequence);
            
            //if(aln.size() > maxAln){
            //    maxAln = aln.size();
            //    _matches = matches;
            //    _strindexes = strindexes;
            //}
        //}
        /*
        cout << endl << endl << endl;
        for(size_t i=0; i<_matches.size(); i++){
            cout << _matches[i] << endl;
        }
        cout << endl;
        for(size_t i=0; i<_strindexes.size(); i++){
            cout << _strindexes[i] << endl;
        }
        */
        _matches = matches;
        _strindexes = strindexes;
    }
    
    string printAlignmentStrings(const vector<int>& matches, const vector<int>& strindexes, const vector<u_int64_t>& sequence){

        size_t ii=0; 
        size_t jj=0;

        vector<u_int64_t> aln1;
        vector<u_int64_t> aln2;
        string alignStr1 = "";
        
        for(int i : strindexes){
            if(i == _none){
                alignStr1 += "- ";
                aln1.push_back(-1);
            }
            else{
                alignStr1 += to_string(sequence[i]) + " ";
                aln1.push_back(sequence[i]);
            }
        }

        string alignStr2 = "";

        for(int j : matches){
            if(j == _none){
                alignStr2 += "- ";
                aln2.push_back(-1);
            }
            else{
                alignStr2 += to_string(_graph->_nodes[j]->_unitigIndex) + " ";
                aln2.push_back(_graph->_nodes[j]->_unitigIndex);
            }
        }

        
        string alignStr = "";

        for(size_t i=0; i<aln1.size(); i++){
            if(aln1[i] == -1){
                alignStr += "-";
            }
            else if(aln1[i] == aln2[i]){
                alignStr += "M";
            }
            else{
                alignStr += "U";
            }
            
        }

        //cout << alignStr1 << endl;
        //cout << alignStr2 << endl;
        //cout << alignStr << endl;

        return alignStr;

    }

    void addAlignmentToGraph(const vector<u_int64_t>& sequence){

        //cout << "---" << endl;
        Node * prevNode = nullptr; //_graph->_nodes[0];

        //for(Node* node : _graph->_nodes){
        //    node->_alignedTo.clear();
        //}

        for(size_t i=0; i<_matches.size(); i++){
            int sindex = _strindexes[i];
            int matchID = _matches[i];

            //cout << sindex << " " << matchID << endl;
            if(sindex == _none) continue;

            u_int64_t base = sequence[sindex];
            Node* node = nullptr;

            /*
            if(prevNode != nullptr){
                for(Edge* edge : prevNode->_successors){
                    if(edge->_toNode->_unitigIndex == base){
                        foundNode = edge->_toNode;
                        break;
                    }
                }
            }
            */

            //if(foundNode != nullptr){
            //    node = foundNode;
            //}
            //else 
            if(matchID == _none){
                //cout << "lala1" << endl;
                node = _graph->addNode(base);
            }
            else if(_graph->_nodes[matchID]->_unitigIndex == base){
                node = _graph->_nodes[matchID];
            }
            else{

                Node* foundNode = nullptr;
                vector<Node*> otherAligns = _graph->_nodes[matchID]->_alignedTo;

                //if(foundNode == nullptr){
                
                for(Node* otherNode : otherAligns){
                    if(otherNode->_unitigIndex == base){
                        foundNode = otherNode;
                        break;
                    }
                }
                //}

                

                if(foundNode == nullptr){
                    //cout << "lala2" << endl;
                    node = _graph->addNode(base);

                    node->_alignedTo.clear();
                    for(Node* otherNode : otherAligns){
                        node->_alignedTo.push_back(otherNode);
                    }
                    node->_alignedTo.push_back(_graph->_nodes[matchID]);

                    _graph->_nodes[matchID]->_alignedTo.push_back(node);
                    for(Node* otherNode : otherAligns){
                        otherNode->_alignedTo.push_back(node);
                    }
                }
                else{
                    node = foundNode;
                }

            }

            //cout << (prevNode!=nullptr) << " " << (node!=nullptr) << endl;
            Edge* edge = _graph->getEdge(prevNode, node);
            //cout << "\t" << (edge != nullptr) << endl;
            
            /*
            //bool isDuplicatedEdge = false;
            for(Edge* successor : prevNode->_successors){
                if(successor->_toNode->_unitigIndex == node->_unitigIndex){
                    //cout << "duplicated successor: " << prevNode->toString() << " -> " << node->toString() << endl;
                    //getchar();
                    //successor->_weight += 1;
                    //isDuplicatedEdge = true;
                    edge = successor;
                    node = successor->_toNode;
                    break;
                }
            }
            */

            if(edge == nullptr){
                //if(!isDuplicatedEdge){
                _graph->addEdge(prevNode, node);
                //}
            }
            else{
                edge->_weight += 1;
            }


            prevNode = node;

        }


        //for(Node* node : _graph->_nodes){
        //    cout << node->_nodeIndex << " " << node->_alignedTo.size() << endl;
        //}
    }


    vector<u_int32_t> getAlignment(const vector<u_int64_t>& sequence){

        vector<u_int32_t> alignment;
        if(_graph->_nodes.size() == 0) return alignment;
        
        _graph->computeTopologicalSort();
        initializeDynamicProgrammingData(sequence);
        alignSequenceToGraph(sequence);
        backtrack(sequence);
        cout << printAlignmentStrings(_matches, _strindexes, sequence) << endl;
    

        Node * prevNode = nullptr; //_graph->_nodes[0];

        for(size_t i=0; i<_matches.size(); i++){
            int sindex = _strindexes[i];
            int matchID = _matches[i];

            if(sindex == _none) continue;

            u_int64_t base = sequence[sindex];
            Node* node = nullptr;


            if(matchID == _none){
                //cout << "derp 1" << endl;
                continue;
                //cout << "lala1" << endl;
                //node = _graph->addNode(base);
            }
            else if(_graph->_nodes[matchID]->_unitigIndex == base){
                //cout << "yes 1" << endl; 
                node = _graph->_nodes[matchID];
            }
            else{

                Node* foundNode = nullptr;
                /*
                vector<Node*> otherAligns = _graph->_nodes[matchID]->_alignedTo;

                //if(foundNode == nullptr){
                
                for(Node* otherNode : otherAligns){
                    if(otherNode->_unitigIndex == base){
                        foundNode = otherNode;
                        break;
                    }
                }
                //}
                */
                

                if(foundNode == nullptr){
                    //cout << "derp 2" << endl;
                    continue;
                }
                else{
                    //cout << "yes 2" << endl; 
                    node = foundNode;
                }

            }

            Edge* edge = _graph->getEdge(prevNode, node);

            alignment.push_back(node->_nodeIndex);

            //cout << node->_nodeIndex << " " << node->_unitigIndex << endl;
            prevNode = node;

        }

        return alignment;

    }

    vector<u_int64_t> performCorrection(const vector<u_int64_t>& readMinimizers){

        vector<u_int64_t> readMinimizersCorrected;

        vector<u_int32_t> alignment = getAlignment(readMinimizers);
        if(alignment.size() <= 1) return readMinimizersCorrected;

        for(size_t i=0; i<alignment.size()-1; i++){
            Node* startNode = _graph->_nodes[alignment[i]];
            Node* endNode = _graph->_nodes[alignment[i+1]];

            //cout << "Find path: " << startNode->_nodeIndex << " -> " << endNode->_nodeIndex << endl;
            vector<u_int32_t> nodePath = getMostSupportedPath(startNode, endNode);
            //cout << "\t" << nodePath.size() << endl;

            if(nodePath.size() < 2){
                cout << "small path?" << endl;
                //getchar();
                continue;
            }

            for(size_t i=0; i<nodePath.size()-1; i++){
                readMinimizersCorrected.push_back(_graph->_nodes[nodePath[i]]->_unitigIndex);
            }
        }

        readMinimizersCorrected.push_back(_graph->_nodes[alignment[alignment.size()-1]]->_unitigIndex);
        

        return readMinimizersCorrected;
    }

    vector<u_int32_t> getMostSupportedPath(Node* startNode, Node* endNode){


		//u_int32_t unitigIndexSource = bin._unitigIndexSource;
		//vector<float> coveragesSource;// = getUnitigCoverages(_progressiveAbundanceFilter->_unitigGraph->_nodes[unitigIndexSource]);


		unordered_map<u_int32_t, u_int32_t> prev;
		list<u_int32_t> queue;
		//unordered_set<u_int32_t> isVisited;

		//isVisited.insert(startNode->_nodeIndex);
		queue.push_back(startNode->_nodeIndex);
        prev[startNode->_nodeIndex] = -1;
        //path.push_back(startNode->_nodeIndex);

        bool found = false;

		while(!queue.empty()){

            //cout << "lala" << endl;
			u_int32_t nodeIndex = queue.front();
			queue.pop_front();

            //if(startNode->_nodeIndex == 32 && endNode->_nodeIndex == 33){
            //    cout << nodeIndex << endl;
            //    getchar();
            //}
            //cout << queue.size() << " " << nodeIndex << endl;

            vector<u_int32_t> solidSuccessors = getSolidSuccessors(nodeIndex);

            for(u_int32_t nodeIndexSuccessor : solidSuccessors){

                //cout << "\t" << nodeIndexSuccessor << endl;

                if(nodeIndexSuccessor == endNode->_nodeIndex){
                    found = true;
                }

                prev[nodeIndexSuccessor] = nodeIndex;
                queue.push_back(nodeIndexSuccessor);

            }


            if(found) break;


		}

        vector<u_int32_t> path;

        if(found){
            u_int32_t nodeIndex = endNode->_nodeIndex;
            while(nodeIndex != -1){
                path.push_back(nodeIndex);
                nodeIndex = prev[nodeIndex];

            }
        }

        std::reverse(path.begin(), path.end());

		return path;

    }
    
    vector<u_int32_t> getSolidSuccessors(u_int32_t nodeIndex){

        vector<u_int32_t> solidSuccessors;

        double maxWeight = 0;

        for(Edge* successorEdge : _graph->_nodes[nodeIndex]->_successors){
            if(successorEdge->_weight > maxWeight){
                maxWeight = successorEdge->_weight;
            }
        }

        float minWeight = maxWeight * 0.25;

        for(Edge* successorEdge : _graph->_nodes[nodeIndex]->_successors){
            if(successorEdge->_weight >= minWeight){
                solidSuccessors.push_back(successorEdge->_toNode->_nodeIndex);
            }
        }

        return solidSuccessors;
    }

};

};


#endif