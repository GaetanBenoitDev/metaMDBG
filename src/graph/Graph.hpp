


#ifndef MDBG_METAG_GRAPH
#define MDBG_METAG_GRAPH


#include <iostream>
#include <unordered_map>
#include <vector>
#include <queue>

using namespace std;
typedef pair<int, int> iPair;



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

/*
    void create(vector<GraphNode>* nodesInfos){
        _unitigs_length.resize(nodesInfos->size(), 0);
        for(size_t i=0; i<nodesInfos->size(); i++){
            //cout << (*nodesInfos)[i]._nodeID << " " << (*nodesInfos)[i]._length << endl;
            _unitigs_length[(*nodesInfos)[i]._nodeID] = (*nodesInfos)[i]._length;
        }
    }
*/
/*
    u_int64_t _current_node_id;
    unordered_map<string, u_int64_t> _name_to_id;
    vector<string> _id_to_name;


    vector<u_int64_t> _unitigs_length;

    void addNode(const string& name){
        if (_name_to_id.find(name) != _name_to_id.end()) return;

        _name_to_id[name] = _current_node_id;
        _id_to_name.push_back(name);
        _current_node_id += 1;
    }

    u_int64_t name_to_id(const string& name){
        return _name_to_id[name];
    }

    string& id_to_name(u_int64_t id){
        return _id_to_name[id];
    }

    void create(vector<GraphNode>* nodesInfos){
        _unitigs_length.resize(nodesInfos->size(), 0);
        for(size_t i=0; i<nodesInfos->size(); i++){
            //cout << (*nodesInfos)[i]._nodeID << " " << (*nodesInfos)[i]._length << endl;
            _unitigs_length[(*nodesInfos)[i]._nodeID] = (*nodesInfos)[i]._length;
        }
    }
    */
};











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

        /*
        collectNeighbors(from, 1, _neighbors);
        cout << _neighbors.size() << endl;
        for(u_int32_t n : _neighbors){
            if(n == to){
                return false;
            }
        }*/

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
        /*
        for (unsigned i = 0; i < _nbEdges; i++)  {
            int start_ver = (*edges)[i].start_ver;
            int end_ver = (*edges)[i].end_ver;
            int weight = (*edges)[i].weight;
            // insert in the beginning
            adjNode* newNode = getAdjListNode(end_ver, weight, directionFrom, (*edges)[i].directionTo, (*edges)[i].isBidirection, _nodes[start_ver]);
             
            // point head pointer to new node
            _nodes[start_ver] = newNode;
        }

        for (unsigned i = 0; i < _nbEdges; i++)  {
            int start_ver = (*edges)[i].end_ver;
            int end_ver = (*edges)[i].start_ver;
            int weight = (*edges)[i].weight;
            // insert in the beginning

            adjNode* newNode = getAdjListNode(end_ver, weight, (*edges)[i].directionTo, (*edges)[i].directionFrom, true, _nodes[start_ver]);
             
            // point head pointer to new node
            _nodes[start_ver] = newNode;
        }*/

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
    /*
    // print all adjacent vertices of given vertex
    //void display_AdjList(adjNode* ptr, int i)
    void display_AdjList(u_int64_t nodeID){
        adjNode* node = _nodes[nodeID];
        while (node != nullptr) {
            cout << "(" << id_to_name(nodeID) << ", " << id_to_name(node->val)
                << ", " << node->overlap << ") ";
            node = node->next;
        }
        cout << endl;
    }
    */

   /*
	static bool unitigLength_sorter(GraphNode& x, GraphNode& y) { return x._length > y._length; }

    void getUnitigSortedByLength(vector<GraphNode>& result){

        result.clear();

        //vector<GraphNode> nodes;
        for(size_t i=0; i<_unitigs_length.size(); i++){
            result.push_back({i, _unitigs_length[i]});
        }

	    //std::sort(result.begin(), result.end(), unitigLength_sorter);

        //for(size_t i=0; i<result.size(); i++){

        //    cout << result[i]._nodeID << " " << result[i]._length << endl;

	    //}


    }*/

    /*
    void collectNeighbors(u_int64_t nodeID, int maxDepth, vector<u_int64_t>& neighbors){

        for(size_t i=0; i<_nbNodes; i++){
            isVisited[i] = false;
        }

        neighbors.clear();
        neighbors.push_back(nodeID);
        isVisited[nodeID] = true;
        collectNeighbors_aux(nodeID, maxDepth, 0, neighbors, isVisited);
        //neighbors.erase(neighbors.begin());
    }

    void collectNeighbors_aux(u_int64_t nodeID, int maxDepth, int currentDepth, vector<u_int64_t>& neighbors, vector<bool>& isVisited){
        if(currentDepth >= maxDepth) return;

        adjNode* node = _nodes[nodeID];
        while (node != nullptr) {

            u_int64_t node_neighbor = node->val;
            //cout << "Neigh: " << node_neighbor << endl;

            if(!isVisited[node->val]){
                isVisited[node->val] = true;
                neighbors.push_back(node_neighbor);
                collectNeighbors_aux(node_neighbor, maxDepth, currentDepth+1, neighbors, isVisited);
            }
            //if(find(neighbors.begin(), neighbors.end(), node_neighbor) == neighbors.end()) {
            //    neighbors.push_back(node_neighbor);
            //    collectNeighbors_aux(node_neighbor, maxDepth, currentDepth+1, neighbors);
            //}

            node = node->next;
        }
    }

    



    */
};



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

class BiGraph{
    
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
    vector<bool> isVisited;
    vector<u_int32_t> distance;
    vector<u_int32_t> prev;

    BiGraph(u_int32_t nbNodes){
        _nbNodes = nbNodes * 2;
        _nodes.resize(_nbNodes, nullptr);

        prev.resize(_nbNodes, 0);
        isVisited.resize(_nbNodes, 0);
        distance.resize(_nbNodes, 0);

        _nbEdges = 0;
       // _nodeIndex = 0;
    }

    ~BiGraph() {
        for (int i = 0; i < _nbNodes; i++){
            delete[] _nodes[i];
        }
    }

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

    bool addEdge(u_int32_t from, bool fromOrient, u_int32_t to, bool toOrient, float weight){

        u_int32_t nodeIndex_from = nodeName_to_nodeIndex(from, fromOrient);
        u_int32_t nodeIndex_to = nodeName_to_nodeIndex(to, toOrient);

        _nodes[nodeIndex_from] = getAdjListNode(nodeIndex_to, weight, _nodes[nodeIndex_from]);
        _nbEdges += 1;

        return true;
    }

    string nodeToString(u_int32_t nodeIndex){
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


};

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
        /*
        if((fromOrient && toOrient) || (fromOrient && !toOrient) || (fromOrient && !toOrient)){
            _nodes[from] = getAdjListNode(to, weight, _nodes[from]);
            _nbEdges += 1;
            return true;
        }
        else if(!fromOrient && !toOrient){
            _nodes[to] = getAdjListNode(from, weight, _nodes[to]);
            _nbEdges += 1;
            return true;
        }
        */

        return false;
    }


};


#endif