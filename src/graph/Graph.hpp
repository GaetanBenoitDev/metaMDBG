


#ifndef MDBG_METAG_GRAPH
#define MDBG_METAG_GRAPH


#include <iostream>
#include <unordered_map>
#include <vector>
#include <queue>

using namespace std;
typedef pair<int, int> iPair;


/*
static u_int32_t unitigName_to_id(string unitig_name){
    unitig_name.erase(unitig_name.begin());
    unitig_name.erase(unitig_name.begin());
    unitig_name.erase(unitig_name.begin());
    unitig_name.erase(unitig_name.end()-1);
    return std::stoull(unitig_name)-1;
}

static string unitigIndex_to_unitigName(u_int32_t unitigIndex){

    string unitigName = "utg";
    string unitig_name_id = to_string(unitigIndex+1);
    size_t nbZeros = 7 - unitig_name_id.size();
    //cout << unitigIndex << " " << nbZeros << endl;
    for(size_t i=0; i<nbZeros; i++){
        unitigName += "0";
    }
    unitigName += unitig_name_id + "l";

    return unitigName;
}
*/

struct adjNode {
    u_int32_t val;
    float weight;
    adjNode* next;
};

struct GraphEdge {
    u_int64_t start_ver;
    u_int64_t end_ver;
    float weight;
};

struct GraphNode {
    u_int64_t _nodeID, _length;
};

struct AdjNode {
    u_int32_t _index;
    //u_int16_t _overlap;
    bool _isRemoved;
};



/*


class GraphInfo{

public:

    GraphInfo(){
        _current_node_id = 0;
    }

    u_int64_t _current_node_id;
    unordered_map<u_int32_t, u_int32_t> _readIndex_to_id;
    vector<u_int32_t> _id_to_readIndex;


    //vector<u_int64_t> _unitigs_length;

    void addNode(u_int32_t readIndex){
        if (_readIndex_to_id.find(readIndex) != _readIndex_to_id.end()) return;

        //cout << "Add node: " << readIndex << " " << _current_node_id << endl;y
        _readIndex_to_id[readIndex] = _current_node_id;
        _id_to_readIndex.push_back(readIndex);
        _current_node_id += 1;
    }

    u_int64_t readIndex_to_id(u_int32_t readIndex){
        return _readIndex_to_id[readIndex];
    }

    u_int32_t id_to_readIndex(u_int64_t id){
        return _id_to_readIndex[id];
    }


};
*/








/*

class AdjGraph{
    
    adjNode* getAdjListNode(u_int32_t to, float weight, adjNode* head, bool isSuccessor)   {
        adjNode* newNode = new adjNode;
        newNode->val = to;
        newNode->weight = weight;

        //if(to == 1650 || to == 1652) cout << "lala " <<  isSuccessor << endl;
        //newNode->isSuccessor = isSuccessor;
        //newNode->to_direction = to_direction;
         
        newNode->next = head;
        return newNode;
    }

public:

    GraphInfo* _graphInfo;
    u_int64_t _nbNodes;
    u_int64_t _nbEdges;
    u_int32_t _nodeIndex;

    vector<adjNode*> _nodes;
    vector<bool> isVisited;
    vector<u_int32_t> distance;
    vector<u_int32_t> prev;

    //adjNode **_nodes;                //adjacency list as array of pointers
    // Constructor

    AdjGraph(){
        _nbNodes = 0;
        _nbEdges = 0;
    }

    AdjGraph(u_int32_t nbNodes){
        _nbNodes = nbNodes;
        _nodes.resize(_nbNodes, nullptr);

        prev.resize(_nbNodes, 0);
        isVisited.resize(_nbNodes, 0);
        distance.resize(_nbNodes, 0);

        _nodeIndex = 0;
    }

    bool addNode(u_int32_t id){
        //cout << "Add node: " << id << endl;
        if(_nodes.size() <= id){
            //cout << "\tPush: " << id << endl;
            _nodes.push_back(nullptr);

            prev.push_back(0);
            isVisited.push_back(0);
            distance.push_back(0);

            _nbNodes += 1;

            return true;
        }

        return false;
    }

    //void addEdge(u_int32_t from, u_int32_t to, float weight){
        //_nodes[from] = getAdjListNode(to, weight, _nodes[from]);
        //_nodes[to] = getAdjListNode(from, weight, _nodes[to]);
    //}

    vector<u_int32_t> _neighbors;

    bool addEdge_checkDuplicate(u_int32_t from, u_int32_t to, float weight, bool isSuccessor){

        //if(from == 1651 || to == 1651){
        //    cout << "TOFL: " << from << " " << to << endl;
        //}
        //cout << "add edge: " << from << " " << to << endl;
        //Check if edge exists
        
        //adjNode* node = _nodes[from];
        //while (node != nullptr) {
        //    if(node->val == to) return false;
        //    node = node->next;
        //}



        //cout << _neighbors.size() << endl;

        _nodes[from] = getAdjListNode(to, weight, _nodes[from], isSuccessor);
        _nbEdges += 1;
        //_nodes[to] = getAdjListNode(from, weight, _nodes[to]);

        return true;
    }


    void create(vector<GraphEdge>* edges, u_int64_t nbNodes, GraphInfo* graphInfo)  {


        _nbNodes = nbNodes;
        _nbEdges = edges->size();
        _graphInfo = graphInfo;

        //prev.resize(_nbNodes, 0);
        //isVisited.resize(_nbNodes, 0);
        //distance.resize(_nbNodes, 0);
        // allocate new node
        //_nodes = new adjNode*[_nbNodes]();
        
        // initialize head pointer for all vertices
        //for (int i = 0; i < _nbNodes; ++i)
        //    _nodes[i] = nullptr;

        // construct directed graph by adding edges to it


    }
    
    ~AdjGraph() {
        for (int i = 0; i < _nbNodes; i++){
            delete[] _nodes[i];
        }

        //for (int i = 0; i < _nbNodes; i++){
        //    delete[] _nodes[i];
        //}
        //delete[] _nodes;
    }


    void collectNeighbors(u_int32_t s, u_int32_t maxDistance, vector<u_int32_t>& neighbors, u_int32_t maxNeighbors){

        neighbors.clear();
        //if(_nodes[s] == nullptr) return;

        //neighbors.push_back(s);
        unordered_set<u_int32_t> isVisitedSet;

        for(size_t i=0; i<_nbNodes; i++){
            //isVisited[i] = false;
            distance[i] = 0;
            //prev[i] = -1;
        }

    
        queue<u_int32_t> queue;
    
        //isVisited[s] = true;
        distance[s] = 0;
        queue.push(s);
    
        while(!queue.empty()){
            
            u_int32_t n = queue.front();
            queue.pop();

            adjNode* node = _nodes[n];

            while (node != nullptr) {

                u_int32_t nn = node->val;

                //if (isVisited[nn]){
                if(isVisitedSet.find(nn) != isVisitedSet.end()){
                    node = node->next;
                    continue;
                }

                if(maxNeighbors > 0 && neighbors.size() >= maxNeighbors){
                    break;
                }

                isVisitedSet.insert(nn);
                //isVisited[nn] = true;
                distance[nn] = distance[n] + 1;
                neighbors.push_back(nn);


                if(distance[nn] >= maxDistance){
                    node = node->next;
                    continue;
                }

                //cout << "Push: " << nn << " " << distance[nn] << endl;
                queue.push(nn);

                node = node->next;
            }

            if(maxNeighbors > 0 && neighbors.size() >= maxNeighbors){
                break;
            }

        }



    }

    void collectNeighbors(u_int32_t s, u_int32_t maxDistance, vector<u_int32_t>& neighbors, u_int32_t maxNeighbors, unordered_set<u_int32_t>& visitedNodes){

        neighbors.clear();
        //if(_nodes[s] == nullptr) return;

        //neighbors.push_back(s);
        unordered_set<u_int32_t> isVisitedSet;

        for(size_t i=0; i<_nbNodes; i++){
            //isVisited[i] = false;
            distance[i] = 0;
            //prev[i] = -1;
        }

    
        queue<u_int32_t> queue;
    
        //isVisited[s] = true;
        distance[s] = 0;
        queue.push(s);
    
        while(!queue.empty()){
            
            u_int32_t n = queue.front();
            queue.pop();

            adjNode* node = _nodes[n];

            while (node != nullptr) {

                u_int32_t nn = node->val;

                //if (isVisited[nn]){
                if(isVisitedSet.find(nn) != isVisitedSet.end()){
                    node = node->next;
                    continue;
                }

                if(visitedNodes.find(nn) != visitedNodes.end()){
                    node = node->next;
                    continue;
                }

                if(maxNeighbors > 0 && neighbors.size() >= maxNeighbors){
                    break;
                }

                isVisitedSet.insert(nn);
                //isVisited[nn] = true;
                distance[nn] = distance[n] + 1;
                neighbors.push_back(nn);


                if(distance[nn] >= maxDistance){
                    node = node->next;
                    continue;
                }

                //cout << "Push: " << nn << " " << distance[nn] << endl;
                queue.push(nn);

                node = node->next;
            }

            if(maxNeighbors > 0 && neighbors.size() >= maxNeighbors){
                break;
            }

        }



    }

    
    void display_AdjList(u_int64_t nodeID, GraphInfo* graphInfo){
        adjNode* node = _nodes[nodeID];
        while (node != nullptr) {
            cout << "(" << graphInfo->id_to_readIndex(nodeID) << ", " << graphInfo->id_to_readIndex(node->val)
                << ", " << node->weight << ") ";
            node = node->next;
        }
        cout << endl;
    }

    u_int32_t shortest_path(u_int32_t src, u_int32_t dest, vector<u_int32_t>& path){

        for(size_t i=0; i<_nbNodes; i++){
            isVisited[i] = false;
            distance[i] = 0;
            prev[i] = -1;
        }
        // Keep track of visited nodes

        // Initialize initial distances as 0 for all nodes

        queue <int> queue1;
        distance[src] = 0;
        queue1.push(src);
        isVisited[src] = true;
        bool found = false;

        while (!queue1.empty() && !found){

            u_int64_t node_current = queue1.front();
            //cout << id_to_name(node_current) << endl;

            queue1.pop();

            adjNode* node = _nodes[node_current];
            while (node != nullptr) {

                u_int64_t node_neighbor = node->val;
                //cout << id_to_name(node_neighbor) << endl;

                if (isVisited[node_neighbor]){
                    node = node->next;
                    continue;
                }

                distance[node_neighbor] = distance[node_current] + 1;
                queue1.push(node_neighbor);
                isVisited[node_neighbor] = true;
                prev[node_neighbor] = node_current;

                if(node_neighbor == dest){
                    found = true;
                    break;
                }

                node = node->next;
            }



        }

        if(found){

            path.clear();
            u_int64_t n = dest;
            while(n != src){
                path.push_back(n);
                //cout << id_to_name(n) << endl;
                n = prev[n];

            }
            path.push_back(src);

            return distance[dest];
        }

        return -1;
    }





    void computeConnectedComponents(vector<vector<u_int32_t>>& components){


        for(size_t i=0; i<_nbNodes; i++){
            isVisited[i] = false;
        }

        for(size_t n=0; n<_nbNodes; n++){
            if(isVisited[n]) continue;


            components.push_back(vector<u_int32_t>());
            vector<u_int32_t>& component = components[components.size()-1];

            queue <int> queue;

            queue.push(n);
            isVisited[n] = true;
            component.push_back(n);

            while (!queue.empty()){

                u_int64_t node_current = queue.front();

                queue.pop();

                adjNode* node = _nodes[node_current];
                while (node != nullptr) {

                    u_int64_t node_neighbor = node->val;

                    if (isVisited[node_neighbor]){
                        node = node->next;
                        continue;
                    }

                    queue.push(node_neighbor);

                    isVisited[node_neighbor] = true;
                    component.push_back(node_neighbor);

                    node = node->next;
                }



            }

        }



    }

};
*/


/*
// graph implementation
int main()
{
    // graph edges array.
    graphEdge edges[] = {
        // (x, y, w) -> edge from x to y with weight w
        {0,1,2},{0,2,4},{1,4,3},{2,3,2},{3,1,4},{4,3,3}
    };
    int N = 6;      // Number of vertices in the graph
    // calculate number of edges
    int n = sizeof(edges)/sizeof(edges[0]);
    // construct graph
    DiaGraph diagraph(edges, n, N);
    // print adjacency list representation of graph
    cout<<"Graph adjacency list "<<endl<<"(start_vertex, end_vertex, weight):"<<endl;
    for (int i = 0; i < N; i++)
    {
        // display adjacent vertices of vertex i
        display_AdjList(diagraph.head[i], i);
    }
    return 0;
}
*/

struct NodeData{
    u_int32_t _abundance;
    u_int32_t _length;
    //u_int32_t _quality;
};

class BiGraph{



    /*
    adjNode* getAdjListNode(u_int32_t to, float weight, adjNode* head)   {
        adjNode* newNode = new adjNode;
        newNode->val = to;
        newNode->weight = weight;

        newNode->next = head;
        return newNode;
    }
    */

public:

    u_int64_t _nbNodes;
    u_int64_t _nbEdges;
    //u_int32_t _nodeIndex;

    vector<vector<AdjNode>> _nodes;
    vector<bool> isVisited;
    vector<u_int32_t> distance;
    vector<u_int32_t> prev;

    vector<NodeData> _nodeDatas;
    //vector<u_int32_t> _nodeLengths;
    
    BiGraph(u_int32_t nbNodes){
        _nbNodes = nbNodes * 2;
        _nodes.resize(_nbNodes);

        prev.resize(_nbNodes, 0);
        isVisited.resize(_nbNodes, 0);
        distance.resize(_nbNodes, 0);

        _nbEdges = 0;
       // _nodeIndex = 0;
    }

    ~BiGraph() {
        //for (int i = 0; i < _nbNodes; i++){
        //    delete[] _nodes[i];
        //}
    }

    /*
    u_int16_t getOverlap(u_int32_t nodeIndex_from, u_int32_t nodeIndex_to){
        if(nodeIndex_from == nodeIndex_to) return 0;
        for(AdjNode& node : _nodes[nodeIndex_from]){
            if(node._index == nodeIndex_to){
                return node._overlap;
            }
        }
        
        for(AdjNode& node : _nodes[nodeIndex_to]){
            if(node._index == nodeIndex_from){
                return node._overlap;
            }
        }

        //cout << "BiGraph::getOverlap : Overlap dosn't exist " << BiGraph::nodeIndex_to_nodeName(nodeIndex_from) << " " << BiGraph::nodeIndex_to_nodeName(nodeIndex_to) << endl;
    
        //getchar();



        return 0;
    }
    */
    /*
    u_int16_t getOverlap_prev(u_int32_t nodeIndex_from, u_int32_t nodeIndex_to){
        for(AdjNode& node : _nodes[nodeIndex_from]){
            if(node._index == nodeIndex_to){
                return node._overlap;
            }
        }
        
        cout << "BiGraph::getOverlap : Overlap dosn't exist " << BiGraph::nodeIndex_to_nodeName(nodeIndex_from) << " " << BiGraph::nodeIndex_to_nodeName(nodeIndex_to) << endl;
        getchar();
        return 0;
    }
    */
    static u_int32_t nodeName_to_nodeIndex(u_int32_t nodeName, bool isOrientationForward){
        if(isOrientationForward){
            return nodeName * 2;
        }
        else{
            return nodeName * 2 + 1;
        }
    }

    static u_int32_t nodeIndex_to_nodeName(u_int32_t nodeIndex, bool& isOrientationForward){
        isOrientationForward = (nodeIndex % 2) == 0;
        return nodeIndex / 2;
    }

    static u_int32_t nodeIndex_to_nodeName(u_int32_t nodeIndex){
        //isOrientationForward = (nodeIndex % 2) == 0;
        return nodeIndex / 2;
    }

    bool addEdge(u_int32_t from, bool fromOrient, u_int32_t to, bool toOrient, u_int16_t overlap){

        u_int32_t nodeIndex_from = nodeName_to_nodeIndex(from, fromOrient);
        u_int32_t nodeIndex_to = nodeName_to_nodeIndex(to, toOrient);

        _nodes[nodeIndex_from].push_back({nodeIndex_to, false});//getAdjListNode(nodeIndex_to, weight, _nodes[nodeIndex_from]);
        _nbEdges += 1;

        return true;
    }

    bool addEdge(u_int32_t nodeIndex_from, u_int32_t nodeIndex_to){

        //u_int32_t nodeIndex_from = nodeName_to_nodeIndex(from, fromOrient);
        //u_int32_t nodeIndex_to = nodeName_to_nodeIndex(to, toOrient);
        u_int16_t overlap = 200;

        _nodes[nodeIndex_from].push_back({nodeIndex_to, false});//getAdjListNode(nodeIndex_to, weight, _nodes[nodeIndex_from]);
        _nbEdges += 1;

        return true;
    }

    bool addEdge_debug(u_int32_t from, u_int32_t to){

        u_int32_t nodeIndex_from = nodeName_to_nodeIndex(from, true);
        u_int32_t nodeIndex_to = nodeName_to_nodeIndex(to, true);
        _nodes[nodeIndex_from].push_back({nodeIndex_to, 0}); //= getAdjListNode(nodeIndex_to, 0, _nodes[nodeIndex_from]);
        _nbEdges += 1;

        nodeIndex_from = nodeName_to_nodeIndex(to, false);
        nodeIndex_to = nodeName_to_nodeIndex(from, false);
        _nodes[nodeIndex_from].push_back({nodeIndex_to, 0});  // = getAdjListNode(nodeIndex_to, 0, _nodes[nodeIndex_from]);
        _nbEdges += 1;

        return true;
    }

    void clearEdgeRemoved(){
        for(size_t i=0; i<_nodes.size(); i++){
            vector<AdjNode>& successors = _nodes[i];
            for(AdjNode& node : successors){
                node._isRemoved = false;
            }
        }
    }

    /*
    bool removeEdge(u_int32_t from, bool fromOrient, u_int32_t to, bool toOrient){

        u_int32_t nodeIndex_from = nodeName_to_nodeIndex(from, fromOrient);
        u_int32_t nodeIndex_to = nodeName_to_nodeIndex(to, toOrient);

        vector<AdjNode>& successors = _nodes[nodeIndex_from];
        // Traversing through the first vector list
        // and removing the second element from it
        for (size_t i=0; i < successors.size(); i++) {
            if (successors[i]._index == nodeIndex_to) {
                successors.erase(successors.begin() + i);
                return true;
            }
        }
    
        return false;
    }
    */

    bool setEdgeRemoved(u_int32_t nodeIndex_from, u_int32_t nodeIndex_to, bool isRemoved){

        vector<AdjNode>& successors = _nodes[nodeIndex_from];

        bool found = false;
        for (size_t i=0; i < successors.size(); i++) {
            if (successors[i]._index == nodeIndex_to) {
                //if(successors[i]._isRemoved) continue;
                successors[i]._isRemoved = isRemoved;
                //cout << "Removed: " <<  i << " " << successors[i]._isRemoved << " " << BiGraph::nodeIndex_to_nodeName(nodeIndex_from) << " " << BiGraph::nodeIndex_to_nodeName(nodeIndex_to) << endl;
                found = true;
                
                //successors.erase(successors.begin() + i);
                _nbEdges -= 1;
                //return true; //On ne breakp as car les aplindrome genere chhacun 2 edge avec la meme source -> dest
            }
        }
    
        return found;
    }

    bool isEdgeRemoved(u_int32_t nodeIndex_from, u_int32_t nodeIndex_to){

        vector<AdjNode>& successors = _nodes[nodeIndex_from];

        for (size_t i=0; i < successors.size(); i++) {
            if (successors[i]._index == nodeIndex_to) {
                if(successors[i]._isRemoved) return true;
            }
        }
    
        return false;
    }

    bool edgeExists(u_int32_t nodeIndex_from, u_int32_t nodeIndex_to){

        vector<AdjNode>& successors = _nodes[nodeIndex_from];
        // Traversing through the first vector list
        // and removing the second element from it
        for (size_t i=0; i < successors.size(); i++) {
            if (successors[i]._index == nodeIndex_to) {
                return true;
            }
        }
    
        return false;
    }

    static string nodeToString(u_int32_t nodeIndex){
        bool orient;
        u_int32_t nodeName = nodeIndex_to_nodeName(nodeIndex, orient);
        string orientStr = "";
        if(orient){
            orientStr = "+";
        }
        else{
            orientStr = "-";
        }
        //string unitigName = unitigIndex_to_unitigName(nodeName);
        string unitigName = to_string(nodeName); //unitigIndex_to_unitigName(nodeName);
        return unitigName + orientStr;
    }



    void collectNeighbors(u_int32_t s, u_int32_t maxDistance, vector<u_int32_t>& neighbors, u_int32_t maxNeighbors, unordered_set<u_int32_t>& visitedNodes){

        neighbors.clear();
        //if(_nodes[s] == nullptr) return;

        //neighbors.push_back(s);
        unordered_set<u_int32_t> isVisitedSet;

        for(size_t i=0; i<_nbNodes; i++){
            //isVisited[i] = false;
            distance[i] = 0;
            //prev[i] = -1;
        }

    
        queue<u_int32_t> queue;
    
        //isVisited[s] = true;
        distance[s] = 0;
        queue.push(s);
    
        while(!queue.empty()){
            
            u_int32_t n = queue.front();
            queue.pop();

            vector<AdjNode>& nodes = _nodes[n];
            //& node = _nodes[n];

            //while (node != nullptr) {
            for(AdjNode& node : nodes){

                u_int32_t nn = node._index;

                //if (isVisited[nn]){
                if(isVisitedSet.find(nn) != isVisitedSet.end()){
                    //node = node->next;
                    continue;
                }

                if(visitedNodes.find(nn) != visitedNodes.end()){
                    //node = node->next;
                    continue;
                }

                if(maxNeighbors > 0 && neighbors.size() >= maxNeighbors){
                    break;
                }

                isVisitedSet.insert(nn);
                //isVisited[nn] = true;
                distance[nn] = distance[n] + 1;
                neighbors.push_back(nn);


                if(distance[nn] >= maxDistance){
                    //node = node->next;
                    continue;
                }

                //cout << "Push: " << nn << " " << distance[nn] << endl;
                queue.push(nn);

                //node = node->next;
            }

            if(maxNeighbors > 0 && neighbors.size() >= maxNeighbors){
                break;
            }

        }



    }




};

class UnitigGraph{

public:

    class Node{
        public:

        u_int32_t _unitigIndex;
        bool _isPalindrome;
        vector<u_int32_t> _nodes;
        //u_int32_t _startNode;
        //u_int32_t _endNode;
        float _abundance;
        u_int32_t _length;
        //u_int32_t _nbNodes;
        //vector<u_int32_t> _nodes;
        vector<Node*> _successors;
        vector<Node*> _predecessors;
        vector<float> _abundances;
        //u_int32_t _sortingIndex;

        u_int32_t _nodeIndexStart_subStart;
        u_int32_t _nodeIndexStart_subEnd;
        u_int32_t _nodeIndexEnd_subStart;
        u_int32_t _nodeIndexEnd_subEnd;

        Node(u_int32_t unitigIndex, const vector<u_int32_t>& nodes, float abundance, u_int32_t length, const vector<NodeData>& nodeDatas){
            _unitigIndex = unitigIndex;
            _nodes = nodes;
            _abundance = abundance;
            _length = length;
            _isPalindrome = false;

            for(u_int32_t nodeIndex : nodes){
                _abundances.push_back(nodeDatas[BiGraph::nodeIndex_to_nodeName(nodeIndex)]._abundance);
            }

            _abundance = computeMedianAbundance(_abundances);

            //_sortingIndex = startNode();
        }

        u_int64_t nodeSum(){
            u_int64_t sum = 0;
            for(u_int32_t nodeIndex : _nodes){
                sum += nodeIndex;
            }
            return sum;
        }

        u_int64_t successorSum(){
            
            u_int64_t sum = 0;
            for(Node* nn : _successors){
                sum += nn->startNode() + nn->endNode();
            }
            return sum;
            /*
            u_int64_t sum = 0;
            for(Node* nn : _successors){
                for(u_int32_t nodeIndex :  nn->_nodes){
                    sum += nodeIndex;
                }
            }
            return sum;
            */
        }

        u_int64_t predecessorSum(){

            u_int64_t sum = 0;
            for(Node* nn : _predecessors){
                sum += nn->startNode() + nn->endNode();
            }
            return sum;
            /*
            u_int64_t sum = 0;
            for(Node* nn : _predecessors){
                for(u_int32_t nodeIndex :  nn->_nodes){
                    sum += nodeIndex;
                }
            }
            return sum;
            */
        }

        u_int32_t startNode(){
            return _nodes[0];
        }

        u_int32_t endNode(){
            return _nodes[_nodes.size()-1];
        }

        bool isCircular(){
            return _nodes.size() > 1 && (_successors.size() == 1) && (_predecessors.size() == 1) && (_successors[0] == this) && (_predecessors[0] == this);  //_nodes.size() > 1 && (startNode() == endNode());
        }
        void mergeWith(Node* node2, u_int32_t kminmerLength, u_int32_t kminmerLengthNonOverlap, const vector<NodeData>& nodeDatas){


            //_length += (node2->_nodes.size()*kminmerLengthNonOverlap);
            for(u_int32_t nodeIndex : node2->_nodes){
                _abundances.push_back(nodeDatas[BiGraph::nodeIndex_to_nodeName(nodeIndex)]._abundance);
            }

            _abundance = computeMedianAbundance(_abundances);
            
            _nodes.insert( _nodes.end(), node2->_nodes.begin(), node2->_nodes.end());

            _length = kminmerLength + ((_nodes.size()-1)*kminmerLengthNonOverlap);

            /*
            for(u_int32_t nodeIndex : _nodes){
                if(nodeIndex == 585233){
                    cout << "haha " << _length << endl;
                    getchar();
                }
                if(nodeIndex == 11822266){
                    cout << "houhou " << _length << endl;
                    getchar();
                }
            }*/

        }

        double compute_median_float(vector<float>& scores){
            size_t size = scores.size();

            if (size == 0){
                return 0;  // Undefined, really.
            }
            else{
                sort(scores.begin(), scores.end());
                if (size % 2 == 0){
                    return (scores[size / 2 - 1] + scores[size / 2]) / 2.0; //scores[size / 2 - 1];
                }
                else {
                    return scores[size / 2];
                }
            }
        }

        double computeMedianAbundance(vector<float> abundances){

            return compute_median_float(abundances);

            /*
            size_t n = abundances.size();

            if (n % 2 == 0) {

                nth_element(abundances.begin(), abundances.begin() + n / 2, abundances.end());
                nth_element(abundances.begin(), abundances.begin() + (n - 1) / 2, abundances.end());

                return (double)(abundances[(n - 1) / 2]+ abundances[n / 2])/ 2.0;
            }

            else {

                nth_element(abundances.begin(), abundances.begin() + n / 2, abundances.end());

                return (double)abundances[n / 2];
            }
            */
        }


    };

	static bool NodeComparatorByLength(UnitigGraph::Node*a, UnitigGraph::Node*b){

        if(a->_nodes.size() == b->_nodes.size()){
            return a->startNode() > b->startNode();
        }
        return a->_nodes.size() > b->_nodes.size();
        //if(a._startNodeIndex == b._startNodeIndex){
        //    return a._length < b._length;
        //}
        //return a._startNodeIndex < b._startNodeIndex;
        //if(a._length == b._length){
        //    return a._startNodeIndex < b._startNodeIndex;
        //}
		//return a._length < b._length;
	}

	static bool NodeComparator(UnitigGraph::Node*a, UnitigGraph::Node*b){
        return a->startNode() < b->startNode();
	}

    //class Edge{
    //    u_int32_t _indexFrom;
    //    u_int32_t _indexTo;
    //};


    u_int32_t _currentUnitigIndex;
    vector<Node*> _nodes;
    vector<u_int32_t> _nodeIndex_to_unitigIndex;
	float _kminmerLength;
	float _kminmerOverlapMean;
	float _kminmerLengthNonOverlap;
    const vector<NodeData> _nodeDatas;

    class GraphChangedIndex{
        public:
        unordered_set<u_int32_t> _superBubbleIndex2_notSuperbubble;
        unordered_map<u_int32_t, vector<u_int32_t>> _superBubbleIndex2_notSuperbubbleNodes;
    };

    GraphChangedIndex _graphChangedIndex;

    UnitigGraph(u_int32_t nbNodes, float kminmerLength, float kminmerOverlapMean, float kminmerLengthNonOverlap, const vector<NodeData>& nodeDatas) : _nodeDatas(nodeDatas){
        _currentUnitigIndex = 0;
        _kminmerLength = kminmerLength;
        _kminmerOverlapMean = kminmerOverlapMean;
        _kminmerLengthNonOverlap = kminmerLengthNonOverlap;
        _nodeIndex_to_unitigIndex.resize(nbNodes, -1);
    }

    ~UnitigGraph(){
        for(size_t i=0; i<_nodes.size(); i++){
            delete _nodes[i];
        }
    }

    void addNode(const vector<u_int32_t>& nodes, float abundance, u_int32_t length){

        u_int32_t unitigIndex = _currentUnitigIndex;

        Node* node = new Node(unitigIndex, nodes, abundance, length, _nodeDatas);

        for(u_int32_t nodeIndex : nodes){
            _nodeIndex_to_unitigIndex[nodeIndex] = unitigIndex;
        }

        _nodes.push_back(node);
        _currentUnitigIndex += 1;
    }

    void addSuccessor(u_int32_t fromUnitigIndex, u_int32_t toUnitigIndex){
        _nodes[fromUnitigIndex]->_successors.push_back(_nodes[toUnitigIndex]);
        //Edge edge = new Edge(fromUnitigIndex, toUnitigIndex);
    }

    void addPredecessors(u_int32_t fromUnitigIndex, u_int32_t toUnitigIndex){
        _nodes[fromUnitigIndex]->_predecessors.push_back(_nodes[toUnitigIndex]);
        //Edge edge = new Edge(fromUnitigIndex, toUnitigIndex);
    }

    void removePredecessor(Node* fromNode, Node* toNode){
        


        /*
        for(UnitigGraph::Node* node : _progressiveAbundanceFilter->_validNodes){
            if(node->startNode() == 2713422){
                cout << (_unitigGraph->_graphChangedIndex._superBubbleIndex2_notSuperbubble.find(node->_unitigIndex) != _unitigGraph->_graphChangedIndex._superBubbleIndex2_notSuperbubble.end()) << endl;
            }
            if(node->startNode() == 2713423){
                cout << (_unitigGraph->_graphChangedIndex._superBubbleIndex2_notSuperbubble.find(node->_unitigIndex) != _unitigGraph->_graphChangedIndex._superBubbleIndex2_notSuperbubble.end()) << endl;
            }
            if(node->endNode() == 2713422){
                cout << (_unitigGraph->_graphChangedIndex._superBubbleIndex2_notSuperbubble.find(node->_unitigIndex) != _unitigGraph->_graphChangedIndex._superBubbleIndex2_notSuperbubble.end()) << endl;
            }
            if(node->endNode() == 2713423){
                cout << (_unitigGraph->_graphChangedIndex._superBubbleIndex2_notSuperbubble.find(node->_unitigIndex) != _unitigGraph->_graphChangedIndex._superBubbleIndex2_notSuperbubble.end()) << endl;
            }
        }
        */

        fromNode->_predecessors.erase(std::remove(fromNode->_predecessors.begin(), fromNode->_predecessors.end(), toNode), fromNode->_predecessors.end());
        toNode->_successors.erase(std::remove(toNode->_successors.begin(), toNode->_successors.end(), fromNode), toNode->_successors.end());

        Node* fromNode_rc = unitigIndex_toReverseDirection(fromNode);
        Node* toNode_rc = unitigIndex_toReverseDirection(toNode);

        fromNode_rc->_successors.erase(std::remove(fromNode_rc->_successors.begin(), fromNode_rc->_successors.end(), toNode_rc), fromNode_rc->_successors.end());
        toNode_rc->_predecessors.erase(std::remove(toNode_rc->_predecessors.begin(), toNode_rc->_predecessors.end(), fromNode_rc), toNode_rc->_predecessors.end());
    
        /*
        if(fromNode->startNode() == 2713422 || fromNode->startNode() == 2713423 || fromNode->endNode() == 2713422 || fromNode->endNode() == 2713423){

            cout << "lala" << endl;
        }
        
        if(toNode->startNode() == 2713422 || toNode->startNode() == 2713423 || toNode->endNode() == 2713422 || toNode->endNode() == 2713423){

            cout << "lala" << endl;
        }
        if(fromNode_rc->startNode() == 2713422 || fromNode_rc->startNode() == 2713423 || fromNode_rc->endNode() == 2713422 || fromNode_rc->endNode() == 2713423){

            cout << "lala" << endl;
        }
        if(toNode_rc->startNode() == 2713422 || toNode_rc->startNode() == 2713423 || toNode_rc->endNode() == 2713422 || toNode_rc->endNode() == 2713423){

            cout << "lala" << endl;
        }
        */
        nodeChanged(fromNode);
        nodeChanged(toNode);
        //nodeChanged(fromNode_rc);
        //nodeChanged(toNode_rc);
    }

    /*
    void removePredecessor(Node* fromNode, Node* toNode, ofstream& debugFile){

        //debugFile << "RP" << fromNode->startNode() << " " << fromNode->endNode() << " " << toNode->startNode() << " " << toNode->endNode() << endl;
        fromNode->_predecessors.erase(std::remove(fromNode->_predecessors.begin(), fromNode->_predecessors.end(), toNode), fromNode->_predecessors.end());
        toNode->_successors.erase(std::remove(toNode->_successors.begin(), toNode->_successors.end(), fromNode), toNode->_successors.end());

        Node* fromNode_rc = unitigIndex_toReverseDirection(fromNode);
        Node* toNode_rc = unitigIndex_toReverseDirection(toNode);

        //debugFile << "RP" << fromNode_rc->startNode() << " " << fromNode_rc->endNode() << " " << toNode_rc->startNode() << " " << toNode_rc->endNode() << endl;

        fromNode_rc->_successors.erase(std::remove(fromNode_rc->_successors.begin(), fromNode_rc->_successors.end(), toNode_rc), fromNode_rc->_successors.end());
        toNode_rc->_predecessors.erase(std::remove(toNode_rc->_predecessors.begin(), toNode_rc->_predecessors.end(), fromNode_rc), toNode_rc->_predecessors.end());
    }
    */
    
    void nodeChanged(Node* node){

        Node* node_rc = unitigIndex_toReverseDirection(node);

        for(u_int32_t unitigIndex : _graphChangedIndex._superBubbleIndex2_notSuperbubbleNodes[node->_unitigIndex]){
            _graphChangedIndex._superBubbleIndex2_notSuperbubble.erase(unitigIndex);
        }
        _graphChangedIndex._superBubbleIndex2_notSuperbubbleNodes[node->_unitigIndex].clear();


        for(u_int32_t unitigIndex : _graphChangedIndex._superBubbleIndex2_notSuperbubbleNodes[node_rc->_unitigIndex]){
            _graphChangedIndex._superBubbleIndex2_notSuperbubble.erase(unitigIndex);
        }
        _graphChangedIndex._superBubbleIndex2_notSuperbubbleNodes[node_rc->_unitigIndex].clear();
    }

    void removeNode(Node* node){


        //if(node->_unitigIndex > 10000000){
        //    cout << "rlala" << endl;
        //    getchar();
        //}
        
        for(Node* nn : node->_successors){
            nn->_predecessors.erase(std::remove(nn->_predecessors.begin(), nn->_predecessors.end(), node), nn->_predecessors.end());
            nodeChanged(nn);
        }
        for(Node* nn : node->_predecessors){
            nn->_successors.erase(std::remove(nn->_successors.begin(), nn->_successors.end(), node), nn->_successors.end());
            nodeChanged(nn);
        }

        Node* node_rc = unitigIndex_toReverseDirection(node);

        for(Node* nn : node_rc->_successors){
            nn->_predecessors.erase(std::remove(nn->_predecessors.begin(), nn->_predecessors.end(), node_rc), nn->_predecessors.end());
            //nodeChanged(nn);
        }
        for(Node* nn : node_rc->_predecessors){
            nn->_successors.erase(std::remove(nn->_successors.begin(), nn->_successors.end(), node_rc), nn->_successors.end());
            //nodeChanged(nn);
        }


        nodeChanged(node);
        //nodeChanged(node_rc);

        node->_unitigIndex = -1;
        node_rc->_unitigIndex = -1;


        //_nodes[node->_unitigIndex] = nullptr;
        //_nodes[node_rc->_unitigIndex] = nullptr;
        //delete node;
        //delete node_rc;
        //node = nullptr;
        //node_rc = nullptr;
    }

    void recompact(Node* node, ofstream& debugFile){

        //debugFile << "recompact: " << node->startNode() << endl;

        while(true){

            bool isChanged = false;

            if(node->_successors.size() == 1){

                Node* nodeSuccessor = node->_successors[0];

                if(nodeSuccessor->_predecessors.size() == 1){

                    if(node->_successors[0] != nodeSuccessor->_predecessors[0]){
                        isChanged = true;
                        mergeNode(node, nodeSuccessor);
                        //debugFile << "recompact: " << node->startNode() << " " << nodeSuccessor->startNode() << endl;
                    }

                }
            }

            /*
            if(node->_predecessors.size() == 1){

                Node* nodePredecessor = node->_predecessors[0];

                isChanged = true;
                mergeNode(node, nodeSuccessor);
                node1->_successors.erase(std::remove(node1->_successors.begin(), node1->_successors.end(), nodeSuccessor), node1->_successors.end());
                node1->_successors = nodeSuccessor->_successors;

                delete successor;
            }
            */

            if(!isChanged) break;
        }

        //debugFile << "recompact: " << node->startNode() << " " << node->_nodes.size() << " " << node->_successors.size() << " " << node->_predecessors.size() << endl;

        /*
        cout << "a" << endl;
        for(Node* nn : node->_successors){
            //cout << nn->_unitigIndex << endl;
            for(Node* nnn : nn->_predecessors){
                //cout << (nnn == nullptr) << endl;
                if(nnn->_unitigIndex > 10000000){
                    cout << "derp" << endl;
                    getchar();
                }
            }
        }
        for(Node* nn : node->_predecessors){
            //cout << nn->_unitigIndex << endl;
            for(Node* nnn : nn->_successors){
                //cout << (nnn == nullptr) << endl;
                if(nnn->_unitigIndex > 10000000){
                    cout << "derp" << endl;
                    getchar();
                }
            }
        }
        cout << "b" << endl;
        */

    }

    void recompact(Node* node){

        while(true){

            bool isChanged = false;

            //if(node->_predecessors.size() == 1 && node->_predecessors[0]->_successors.size() == 1) break;

            if(node->_successors.size() == 1){

                Node* nodeSuccessor = node->_successors[0];

                if(nodeSuccessor->_predecessors.size() == 1){

                    if(node->_successors[0] != nodeSuccessor->_predecessors[0]){
                        isChanged = true;
                        mergeNode(node, nodeSuccessor);
                    }

                }
            }

            /*
            if(node->_predecessors.size() == 1){

                Node* nodePredecessor = node->_predecessors[0];

                isChanged = true;
                mergeNode(node, nodeSuccessor);
                node1->_successors.erase(std::remove(node1->_successors.begin(), node1->_successors.end(), nodeSuccessor), node1->_successors.end());
                node1->_successors = nodeSuccessor->_successors;

                delete successor;
            }
            */

            if(!isChanged) break;
        }

        /*
        cout << "a" << endl;
        for(Node* nn : node->_successors){
            //cout << nn->_unitigIndex << endl;
            for(Node* nnn : nn->_predecessors){
                //cout << (nnn == nullptr) << endl;
                if(nnn->_unitigIndex > 10000000){
                    cout << "derp" << endl;
                    getchar();
                }
            }
        }
        for(Node* nn : node->_predecessors){
            //cout << nn->_unitigIndex << endl;
            for(Node* nnn : nn->_successors){
                //cout << (nnn == nullptr) << endl;
                if(nnn->_unitigIndex > 10000000){
                    cout << "derp" << endl;
                    getchar();
                }
            }
        }
        cout << "b" << endl;
        */

    }

    
    void mergeNode(Node* node1, Node* node2){

        /*
        //if(node2->_unitigIndex == 236531 || node2->_unitigIndex == 236532){
        if(node1->_abundance < 100 && node2->_abundance > 1000){
            cout << "merged" << endl;
            cout << node1->_unitigIndex << " " << node1->_abundance << endl;
            cout << node2->_unitigIndex << " " << node2->_abundance << endl;
            cout << node1->_successors.size() << " " << node1->_predecessors.size() << endl;
            for(Node* nn : node1->_successors){
                cout << nn->_unitigIndex << endl;
            }
            for(Node* nn : node1->_predecessors){
                cout << nn->_unitigIndex << endl;
            }
            getchar();
        }
        */
       
        Node* node1_rc = unitigIndex_toReverseDirection(node1);
        Node* node2_rc = unitigIndex_toReverseDirection(node2);

        nodeChanged(node1);
        //nodeChanged(node1_rc);
        nodeChanged(node2);
        //nodeChanged(node2_rc);

        /*
        for(Node* nn : node1->_successors){
            nodeChanged(nn);
            nodeChanged(unitigIndex_toReverseDirection(nn));
        }
        for(Node* nn : node1->_predecessors){
            nodeChanged(nn);
            nodeChanged(unitigIndex_toReverseDirection(nn));
        }
        for(Node* nn : node2->_successors){
            nodeChanged(nn);
            nodeChanged(unitigIndex_toReverseDirection(nn));
        }
        for(Node* nn : node2->_predecessors){
            nodeChanged(nn);
            nodeChanged(unitigIndex_toReverseDirection(nn));
        }
        for(Node* nn : node1_rc->_successors){
            nodeChanged(nn);
            nodeChanged(unitigIndex_toReverseDirection(nn));
        }
        for(Node* nn : node1_rc->_predecessors){
            nodeChanged(nn);
            nodeChanged(unitigIndex_toReverseDirection(nn));
        }
        for(Node* nn : node2_rc->_successors){
            nodeChanged(nn);
            nodeChanged(unitigIndex_toReverseDirection(nn));
        }
        for(Node* nn : node2_rc->_predecessors){
            nodeChanged(nn);
            nodeChanged(unitigIndex_toReverseDirection(nn));
        }
        */

        //cout << "\tMerging " << node1->_unitigIndex << " with " << node2->_unitigIndex << endl;
        node1->mergeWith(node2, _kminmerLength, _kminmerLengthNonOverlap, _nodeDatas);

        /*
        cout << "----" << endl;
        for(Node* successor : node1->_successors){
            cout << successor->_unitigIndex << endl;
        }
        cout << "----" << endl;
        for(Node* successor : node1_rc->_predecessors){
            cout << successor->_unitigIndex << endl;
        }
        */
        
        //node1->_successors.erase(std::remove(node1->_successors.begin(), node1->_successors.end(), node2), node1->_successors.end());
        node1->_successors = node2->_successors;

        for(Node* nn : node2->_successors){


            //for(Node* nnP : nn->_predecessors){
            //    cout << "haha: " << nnP->_unitigIndex << endl;
            //}

            for(size_t i=0; i<nn->_predecessors.size(); i++){
                if(nn->_predecessors[i] == node2){
                    nn->_predecessors[i] = node1;
                    break;
                }
            }

            //nn->_predecessors.erase(std::remove(nn->_predecessors.begin(), nn->_predecessors.end(), node2), nn->_predecessors.end());
            //nn->_predecessors.push_back(node1);
            /*
            for(Node* nnP : nn->_predecessors){
                if(nnP == node2){
                    cout << "lala1" << endl;
                    nnP = node1;
                    break;
                }
            }
            */

            //for(Node* nnP : nn->_predecessors){
            //    cout << "houhou: " << nnP->_unitigIndex << endl;
            //}

        }

        node1_rc->_length = node1->_length;
        node1_rc->_abundance = node1->_abundance;
        node1_rc->_abundances = node1->_abundances;
        node1_rc->_nodes = node1->_nodes;
        std::reverse(node1_rc->_nodes.begin(), node1_rc->_nodes.end());
        for(size_t i=0; i<node1_rc->_nodes.size(); i++){
            node1_rc->_nodes[i] = nodeIndex_toReverseDirection(node1_rc->_nodes[i]);
        }


        //node1_rc->_predecessors.erase(std::remove(node1_rc->_predecessors.begin(), node1_rc->_predecessors.end(), node2_rc), node1_rc->_predecessors.end());
        node1_rc->_predecessors = node2_rc->_predecessors;

        for(Node* nn : node2_rc->_predecessors){
            for(size_t i=0; i<nn->_successors.size(); i++){
                if(nn->_successors[i] == node2_rc){
                    nn->_successors[i] = node1_rc;
                    break;
                }
            }
        }


        //cout << node1->_unitigIndex << " " << node1_rc->_unitigIndex << " " << node2->_unitigIndex << " " << node2_rc->_unitigIndex << endl;

        //cout << "\tResult: " << node1->_unitigIndex << "  " << node2->_unitigIndex << endl;
        //cout << "\tResult: " << node1_rc->_unitigIndex << "  " << node2_rc->_unitigIndex << endl;

        //_nodes[node2->_unitigIndex] = nullptr;
        //_nodes[node2_rc->_unitigIndex] = nullptr;




        node2_rc->_unitigIndex = -1;
        node2->_unitigIndex = -1;



        //cout << "\tResult: " << node1->_unitigIndex << "  " << node2->_unitigIndex << endl;
        //cout << "\tResult: " << node1_rc->_unitigIndex << "  " << node2_rc->_unitigIndex << endl;
        //delete node2;
        //delete node2_rc;

        /*
        node1->_nodes.insert( node1->_nodes.end(), node2->_nodes.begin(), node2->_nodes.end());
        for(u_int32_t nodeIndex : node1._nodes){

        }
        */
        /*
        cout << "----" << endl;
        for(Node* successor : node1->_successors){
            cout << "S: " << successor->_unitigIndex << endl;
        }
        for(Node* successor : node1->_predecessors){
            cout << "P: " << successor->_unitigIndex << endl;
        }
        cout << "----" << endl;
        for(Node* successor : node1_rc->_successors){
            cout << "S: " << successor->_unitigIndex << endl;
        }
        for(Node* successor : node1_rc->_predecessors){
            cout << "P: " << successor->_unitigIndex << endl;
        }
        */

    }
    //void removeEdgeTip(Node* fromNode, Node* toNode){
        
    //}

    static u_int32_t nodeIndex_toReverseDirection(u_int32_t nodeIndex){
        bool orient;
        u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex, orient);

        u_int32_t nodeIndex_prev;
        if(orient){
            nodeIndex_prev = BiGraph::nodeName_to_nodeIndex(nodeName, false);
        }
        else{
            nodeIndex_prev = BiGraph::nodeName_to_nodeIndex(nodeName, true);
        }

        return nodeIndex_prev;
    }

    u_int32_t unitigIndex_toReverseDirection(u_int32_t unitigIndex){
        if(unitigIndex % 2 == 0){
            return unitigIndex+1;
        }
        else{
            return unitigIndex-1;
        }
    }

    Node* unitigIndex_toReverseDirection(Node* node){
        if(node->_unitigIndex % 2 == 0){
            return _nodes[node->_unitigIndex+1];
        }
        else{
            return _nodes[node->_unitigIndex-1];
        }
    }

    u_int64_t computeChecksum_successor(){

        u_int64_t checksumTotal = 0;

        for(Node* node : _nodes){

            if(node->_unitigIndex == -1) continue;

            //u_int64_t checksum = 0;
            //for(u_int32_t nodeIndex : node->_nodes){
            //    checksum += nodeIndex;
            //}
            //checksum *= node->_nodes.size();

            //checksumTotal += checksum;

            checksumTotal += node->startNode() + node->endNode();

            for(Node* nn : node->_successors){
                checksumTotal += (nn->startNode() + nn->endNode()) * node->_successors.size();
            }
            for(Node* nn : node->_predecessors){
                checksumTotal += (nn->startNode() + nn->endNode()) * node->_predecessors.size();
            }
        }

        return checksumTotal;
    }

    u_int64_t computeChecksum_nodeIndex(){

        u_int64_t checksumTotal = 0;

        for(Node* node : _nodes){

            if(node->_unitigIndex == -1) continue;

            u_int64_t checksum = 0;
            for(u_int32_t nodeIndex : node->_nodes){
                checksum += nodeIndex;
            }
            checksum *= node->_nodes.size();

            checksumTotal += checksum;
        }

        return checksumTotal;
    }

    u_int64_t computeChecksum_abundance(){

        u_int64_t checksumTotal = 0;

        for(Node* node : _nodes){

            if(node->_unitigIndex == -1) continue;

            u_int64_t checksum = 0;
            for(u_int32_t ab : node->_abundances){
                checksum += ab;
            }
            checksum *= node->_nodes.size();

            checksumTotal += checksum;
        }

        return checksumTotal;
    }

    u_int64_t nbUnitigs(){
        u_int64_t nbUnitigs = 0;
        for(size_t i=0; i<_nodes.size(); i++){

            UnitigGraph::Node* node = _nodes[i];
            if(node->_unitigIndex == -1) continue;
            nbUnitigs += 1;
        }

        return nbUnitigs;
    }

    void save(const string& outputFilename, float referenceAbundance){

        //_logFile << "Saving unitig graph: " << outputFilename << endl;
        
		unordered_set<u_int32_t> writtenUnitigs;

        unordered_set<u_int32_t> selectedUnitigIndex;

        ofstream colorFile(outputFilename + "_color.csv");
        colorFile << "Name,Color" << endl;

        ofstream file_nodeNameToUnitigIndex(outputFilename + ".nodeToUnitig");

        for(Node* node : _nodes){
            //if(node == nullptr) continue;
            if(node->_unitigIndex == -1) continue;
            //cout << node->_unitigIndex << endl;

			if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(node->startNode())) != writtenUnitigs.end()) continue;
			if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(node->endNode())) != writtenUnitigs.end()) continue;

			writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(node->startNode()));
			writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(node->endNode()));

            selectedUnitigIndex.insert(node->_unitigIndex);

            for(u_int32_t nodeIndex : node->_nodes){
                u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
		        file_nodeNameToUnitigIndex.write((const char*)&nodeName, sizeof(nodeName));
		        file_nodeNameToUnitigIndex.write((const char*)&node->_unitigIndex, sizeof(node->_unitigIndex));
            }

            if(referenceAbundance != 0){
                if(node->_abundance / referenceAbundance < 1.3){
                    colorFile << node->_unitigIndex << "," << "blue" << endl;
                }
                else{
                    colorFile << node->_unitigIndex << "," << "red" << endl;
                }
            }
        }

        file_nodeNameToUnitigIndex.close();

        
        ofstream outputFile(outputFilename);


        unordered_set<u_int32_t> isVisited;



        for(Node* node : _nodes){
            //if(node == nullptr) continue;
            if(node->_unitigIndex == -1) continue;


            if(selectedUnitigIndex.find(node->_unitigIndex) == selectedUnitigIndex.end()) continue;
            if(isVisited.find(node->_unitigIndex) != isVisited.end()) continue;

            queue<u_int32_t> queue;
            queue.push(node->_unitigIndex);

            while(queue.size() > 0){
                u_int32_t unitigIndex = queue.front();
                queue.pop();

                //cout << unitigIndex << endl;
                //cout << unitigIndex << endl;
                //cout << (_nodes[unitigIndex] == nullptr) << endl;
                if(isVisited.find(unitigIndex) != isVisited.end()) continue;
                isVisited.insert(unitigIndex);

                string ori = "+";
                if(selectedUnitigIndex.find(unitigIndex) == selectedUnitigIndex.end()){
                    unitigIndex = unitigIndex_toReverseDirection(unitigIndex);
                    ori = "-";
                }

                //isVisited.insert(unitigIndex_toReverseDirection(unitigIndex));
                //linkedUnitigIndex.insert(unitigIndex);

                outputFile << "S" << "\t" << unitigIndex << "\t" << "*" << "\t" << "LN:i:" << _nodes[unitigIndex]->_length << "\t" << "dp:i:" << _nodes[unitigIndex]->_abundance << endl;


                //vector<u_int32_t> successors;
                //getSuccessors_unitig(unitigIndex, 0, successors);

                for(Node* nodeSuccessor : _nodes[unitigIndex]->_successors){
                    

                    u_int32_t unitigIndexN = nodeSuccessor->_unitigIndex;

                    string ori2 = "+";
                    if(selectedUnitigIndex.find(unitigIndexN) == selectedUnitigIndex.end()){
                        unitigIndexN = unitigIndex_toReverseDirection(unitigIndexN);
                        ori2 = "-";
                    }
                    
                    u_int32_t overlap = 600;
                    outputFile << "L" << "\t" << unitigIndex << "\t" << ori << "\t" << unitigIndexN << "\t" << ori2 << "\t" << overlap << "M" << endl;
                    queue.push(unitigIndexN);

                    //if(unitigIndexN > 100000){
                        //cout << "wtf succ: " << unitigIndex << " -> " << unitigIndexN << endl;
                        //getchar();
                    //}
                }

                
                //vector<u_int32_t> predecessors;
                //getPredecessors_unitig(unitigIndex, 0, predecessors);
                //for(u_int32_t unitigIndexN : predecessors){

                for(Node* nodePredecessor : _nodes[unitigIndex]->_predecessors){


                    u_int32_t unitigIndexN = nodePredecessor->_unitigIndex;

                    string ori2 = "+";
                    if(selectedUnitigIndex.find(unitigIndexN) == selectedUnitigIndex.end()){
                        unitigIndexN = unitigIndex_toReverseDirection(unitigIndexN);
                        ori2 = "-";
                    }
                    
                    u_int32_t overlap = 600;
                    outputFile << "L" << "\t" << unitigIndexN << "\t" << ori2 << "\t" << unitigIndex << "\t" << ori << "\t" << overlap << "M" << endl;
                    queue.push(unitigIndexN);


                    //if(unitigIndexN > 100000){
                    //    cout << "wtf pred: " << unitigIndex << " -> " << unitigIndexN << endl;
                    //    getchar();
                    //}
                }
                

            }
        }



        outputFile.close();
        colorFile.close();

        //_logFile << "\tdone" << endl;
    }
};

/*
class Graph2{

public:

    u_int64_t _nbNodes;
    u_int64_t _nbEdges;

    vector<vector<AdjNode>> _successors;
    vector<vector<AdjNode>> _predecessors;
    
    Graph2(){
        //_nbNodes = nbNodes;
        //_successors.resize(_nbNodes);
        //_predecessors.resize(_nbNodes);

        _successors.clear();
        _predecessors.clear();
        _nbNodes = 0;
        _nbEdges = 0;
    }

    unordered_map<u_int32_t, u_int32_t> _nodeName_to_nodeIndex;
    vector<u_int32_t> _nodeIndex_to_nodeName;


    //vector<u_int64_t> _unitigs_length;

    void addNode(u_int32_t nodeName){
        if (_nodeName_to_nodeIndex.find(nodeName) != _nodeName_to_nodeIndex.end()) return;

        _nodeName_to_nodeIndex[nodeName] = _nodeIndex_to_nodeName.size();
        _nodeIndex_to_nodeName.push_back(nodeName);
    }

    u_int32_t nodeName_to_nodeIndex(u_int32_t nodeName){
        return _nodeName_to_nodeIndex[nodeName];
    }

    u_int32_t nodeIndex_to_nodeName(u_int32_t nodeIndex){
        return _nodeIndex_to_nodeName[nodeIndex];
    }

    void addEdge(u_int32_t from, u_int32_t to){

        u_int32_t nodeIndex_from = nodeName_to_nodeIndex(from);
        u_int32_t nodeIndex_to = nodeName_to_nodeIndex(to);
        _successors[nodeIndex_from].push_back({nodeIndex_to, 0}); //= getAdjListNode(nodeIndex_to, 0, _nodes[nodeIndex_from]);
        _predecessors[nodeIndex_to].push_back({nodeIndex_from, 0}); //= getAdjListNode(nodeIndex_to, 0, _nodes[nodeIndex_from]);
        _nbEdges += 1;
    }

    void removeEdge(u_int32_t from, u_int32_t to){

        u_int32_t nodeIndex_from = nodeName_to_nodeIndex(from);
        u_int32_t nodeIndex_to = nodeName_to_nodeIndex(to);

        vector<AdjNode>& successors = _successors[nodeIndex_from];
        
        for (size_t i=0; i < successors.size(); i++) {
            if (successors[i]._index == nodeIndex_to) {
                successors.erase(successors.begin() + i);
                break;
            }
        }
    
        vector<AdjNode>& predecessors = _predecessors[nodeIndex_to];
        
        for (size_t i=0; i < predecessors.size(); i++) {
            if (predecessors[i]._index == nodeIndex_from) {
                predecessors.erase(predecessors.begin() + i);
                break;
            }
        }

    }
    
    u_int32_t in_degree(u_int32_t v){
        return _predecessors[nodeName_to_nodeIndex(v)].size();
    }
    
    u_int32_t out_degree(u_int32_t v){
        return _successors[nodeName_to_nodeIndex(v)].size();
    }

    void getSuccessors(u_int32_t v, vector<u_int32_t>& successors){
        successors.clear();
        for(AdjNode& node : _successors[nodeName_to_nodeIndex(v)]){
            successors.push_back(nodeIndex_to_nodeName(node._index));
        }
    }
    
    void getPredecessors(u_int32_t v, vector<u_int32_t>& predecessors){
        predecessors.clear();
        for(AdjNode& node : _predecessors[nodeName_to_nodeIndex(v)]){
            predecessors.push_back(nodeIndex_to_nodeName(node._index));
        }
    }

};

*/
/*
class UnitigGraph{
    
    adjNode* getAdjListNode(u_int32_t to, float weight, adjNode* head)   {
        adjNode* newNode = new adjNode;
        newNode->val = to;
        newNode->weight = weight;

        newNode->next = head;
        return newNode;
    }

public:

    u_int64_t _nbNodes;
    u_int64_t _nbEdges;
    //u_int32_t _nodeIndex;

    vector<adjNode*> _nodes;
    //vector<bool> isVisited;
    //vector<u_int32_t> distance;
    //vector<u_int32_t> prev;

    UnitigGraph(u_int32_t nbNodes){
        _nbNodes = nbNodes;
        _nodes.resize(_nbNodes, nullptr);

        //prev.resize(_nbNodes, 0);
        //isVisited.resize(_nbNodes, 0);
        //distance.resize(_nbNodes, 0);

        _nbEdges = 0;
    }

    ~UnitigGraph() {
        for (int i = 0; i < _nbNodes; i++){
            delete[] _nodes[i];
        }
    }

    bool addEdge(u_int32_t from, bool fromOrient, u_int32_t to, bool toOrient, float weight){

        //if(from == 798 || to == 798){
        //  cout << unitigIndex_to_unitigName(from) << " " << fromOrient << " " << unitigIndex_to_unitigName(to) << " " << toOrient << endl;
        //}
        if((fromOrient && toOrient)){
            //if(from == 798 || to == 798){
            //    cout << "1:     " << unitigIndex_to_unitigName(from) << " -> " << unitigIndex_to_unitigName(to) << endl;
            //}
            _nodes[from] = getAdjListNode(to, weight, _nodes[from]);
            _nbEdges += 1;
            return true; 
        }
        else if((fromOrient && !toOrient)){
            //if(from == 798 || to == 798){
            //    cout << "2:    " << unitigIndex_to_unitigName(to) << " -> " << unitigIndex_to_unitigName(from) << endl;
            //}
            _nodes[to] = getAdjListNode(from, weight, _nodes[to]);
            _nbEdges += 1;
            return true; 
        }


        return false;
    }


};
*/



#endif