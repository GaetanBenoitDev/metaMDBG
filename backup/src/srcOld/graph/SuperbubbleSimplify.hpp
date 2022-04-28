


#ifndef MDBG_METAG_SUPERBUBBLE
#define MDBG_METAG_SUPERBUBBLE

class GraphSimplify;

#include "../Commons.hpp"


//#include "GraphSimplify.hpp"








/*

    //https://networkx.org/documentation/networkx-1.9/_modules/networkx/algorithms/components/strongly_connected.html#strongly_connected_components
    void getStronglyConnectedComponents(vector<vector<u_int32_t>>& sccs){

        sccs.clear();

        compact();
        cout << "Nb unitigs: " << _unitigs.size() << endl;

        unordered_map<u_int32_t, u_int32_t> preorder;
        unordered_map<u_int32_t, u_int32_t> lowlink;
        unordered_set<u_int32_t> scc_found;
        vector<u_int32_t> scc_queue;
        vector<u_int32_t> queue;

        u_int64_t i = 0;

        
        for(Unitig& unitig_source : _unitigs){
            //for(u_int32_t nodeIndex : _isNodeValid2){
            if(scc_found.find(unitig_source._index) == scc_found.end()){
                

                queue.clear();
                queue.push_back(unitig_source._index);

                while(queue.size() > 0){
                    u_int32_t v = queue[queue.size()-1];

                    if(preorder.find(v) == preorder.end()){
                        i += 1;
                        preorder[v] = i;
                    }

                    bool done = true;

                    vector<u_int32_t> v_nbrs;
                    getSuccessors_unitig(v, v_nbrs); //_unitigs[v]._endNode
                    
                    for(u_int32_t w : v_nbrs){
                        if(preorder.find(w) == preorder.end()){
                            queue.push_back(w);
                            done = false;
                            break;
                        }
                    }

                    if(done){
                        lowlink[v] = preorder[v];
                        for(u_int32_t w : v_nbrs){
                            if(scc_found.find(w) == scc_found.end()){
                                if(preorder[w] > preorder[v]){
                                    lowlink[v] = min(lowlink[v], lowlink[w]);
                                }
                                else{
                                    lowlink[v] = min(lowlink[v], preorder[w]);
                                }
                            }
                        }

                        queue.pop_back();

                        if(lowlink[v] == preorder[v]){
                            scc_found.insert(v);
                            vector<u_int32_t> scc = {v};
                            while(scc_queue.size() > 0 && preorder[scc_queue[scc_queue.size()-1]] > preorder[v]){
                                u_int32_t k = scc_queue[scc_queue.size()-1];
                                scc_queue.pop_back();
                                scc_found.insert(k);
                                scc.push_back(k);
                            }
                            sccs.push_back(scc);
                        }
                        else{
                            scc_queue.push_back(v);
                        }
                    }
                }
            }
        }

        
        cout << "Nb sccs: " << sccs.size() << endl;
		for(vector<u_int32_t>& scc : sccs){
			cout << "----" << endl;
			for(u_int32_t unitigIndex : scc){
				vector<u_int32_t> nodes;
				getUnitigNodes(_unitigs[unitigIndex], nodes);
				for(u_int32_t nodeIndex : nodes){
					cout << BiGraph::nodeIndex_to_nodeName(nodeIndex) << endl;
				}
			}
		}
        
    }

    void subgraph(Graph2* graph, const unordered_set<u_int32_t>& nodes){
        graph = new Graph2();

        for(u_int32_t v : nodes){
            graph->addNode(v);

            vector<u_int32_t> successors;
            getSuccessors_unitig(v, successors);
            for(u_int32_t unitigIndex : successors){
                graph->addEdge(v, unitigIndex);
            }
        }
    }






    class AuxiliaryGraph{
        
    public:

        unordered_set<u_int32_t> vertices;
        u_int32_t number;
        u_int32_t artificial_sink;
        u_int32_t artificial_source;
        //string artificial_sink;
        //string artificial_source;
        Graph2* _g;

        AuxiliaryGraph(u_int32_t n){
            vertices.clear();
            number = (-n*2) - 1;
            artificial_sink =  number;//"artificial_sink_" + to_string(number);
            artificial_source = number-1; //"artificial_source_" + to_string(number);
            _g = nullptr;
        }

        void add(u_int32_t v){
            vertices.insert(v);
        }
        
        void add_edge(u_int32_t v, u_int32_t v2){
            _g->addEdge(v, v2);
        }
        
        void copy_graph(GraphSimplify* g){
            g->subgraph(_g, vertices);
            //_g = new Graph2();
            //for(u_int32_t v : vertices){
            //    _g->addNode(v);
            //}
            _g->addNode(artificial_sink);
            _g->addNode(artificial_source);
        }
        
        void connect2sink(u_int32_t v){
            _g->addEdge(v, artificial_sink);
        }
        
        void connect2source(u_int32_t v){
            _g->addEdge(artificial_source, v);
        }
        
        u_int32_t in_degree(u_int32_t v){
            return _g->in_degree(v);
        }
        
        u_int32_t out_degree(u_int32_t v){
            return _g->out_degree(v);
        }
        
        bool source_connected(){
            return _g->out_degree(artificial_source) != 0;
        }
        
        bool sink_connected(){
            return _g->in_degree(artificial_sink) != 0;
        }
        
        void remove_edge(u_int32_t f, u_int32_t t){
            _g->removeEdge(f, t);
        }
        
        void getSuccessors(u_int32_t v, vector<u_int32_t>& successors){
            _g->getSuccessors(v, successors);
        }
        
        void getPredecessors(u_int32_t v, vector<u_int32_t>& predecessors){
            _g->getPredecessors(v, predecessors);
        }
        
        //void nodes(){
            //return self.g.nodes
        //}

        u_int32_t a(){
            return artificial_source;
        }

        u_int32_t b(){
            return artificial_sink;
        }

        u_int32_t nr(){
            return number;
        }
        





    };















    class SuperbubbleSimplify{

    public:

        GraphSimplify* _g;

        SuperbubbleSimplify(GraphSimplify* g){
            _g = g;
        }

        void execute(){

            cout << "Start superbubble detection" << endl;

            AuxiliaryGraph* dag;
            vector<AuxiliaryGraph*> scc;

            get_strongly_connected_component(_g, dag, scc);

            for(AuxiliaryGraph* c : scc){
                create_auxiliary_graph(c, _g);
            }











            cout << "Todo: delete all memory allocation" << endl;

        }

        void get_strongly_connected_component(GraphSimplify* g, AuxiliaryGraph* dag, vector<AuxiliaryGraph*>& non_singelton){

            vector<vector<u_int32_t>> sccs;
            g->getStronglyConnectedComponents(sccs);

            
            dag = new AuxiliaryGraph(0);

            u_int32_t n = 1;
            
            for(vector<u_int32_t>& scc : sccs){
                bool isNonSingleton = false;
                if(scc.size() > 1){
                    isNonSingleton = true;
                }
                else if(scc.size() == 1){
                    u_int32_t unitigIndex = scc[0];
                    vector<u_int32_t> successors;
                    g->getSuccessors_unitig(unitigIndex, successors);
                    if(std::find(successors.begin(), successors.end(), unitigIndex) != successors.end()){
                        isNonSingleton = true;
                    }
                }

                if(isNonSingleton){
                    AuxiliaryGraph* nsscc = new AuxiliaryGraph(n);
                    n += 1;
                    for(u_int32_t v : scc){
                        nsscc->add(v);
                    }
                    non_singelton.push_back(nsscc);
                }
                else{
                    dag->add(scc[0]);
                }
            }

        }

        void create_auxiliary_graph(AuxiliaryGraph* c, GraphSimplify* g){
            c->copy_graph(g);
            
            for(u_int32_t v : c->vertices){
                if(g->in_degree(v) > c->in_degree(v) || c->in_degree(v) == 0){
                    c->connect2source(v);
                }
                if(g->out_degree(v) > c->out_degree(v) || c->out_degree(v) == 0){
                    c->connect2sink(v);
                }
            }

        }


    };



*/


#endif