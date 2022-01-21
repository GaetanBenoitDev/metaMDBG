
//	"reprise: retester IsImproved pendant pathExploration, si on revient au noeud source sans amelioration on empeche ce path
//	hifiasm ground truth: utiliser contig plutot que unitig pour avoir un meilleur ordering
//	"

//PB: quand zone complexe tres fragmenté en unitig (succession d'unitig tres court), le bestPrevUnitigRank n'est pas très faible si ça se joue à un unitig prêt, ajouter une notion 
//collectPossibleSuccessors: ici on ne veut prune que les chemins dont on est sûr à 100% qu'il ne peuvent pas etre successeur, donc plutot une methode permissive, les loop font qu'on peut rater le bon successeur si on est trop stringent
//Idée: Utiliser strongly connected component -> si solved, on parcours le chemin résolu et on vérifie qu'il n'y a pas d'alternative dû a un génome incomplet (grande tips probablement), si chemin alternatif, on cancel le chemin et on l'assemble avec un connected component classique

/*
		CollectPossibleSuccessor: Mieux gérer la récursivité et le currentDepth
		"reprise: isReachable2: ajouter le maxVisitable?
		Opti: CollectPossibleSUccessors: reset du nbVisitedTimes quand node non visité: on peut faire ça uniquement si on a changé d'unitigs par rapport au node précédent 
		Detection de superbubble:
			- Ajouter un quick patch :( boucle infini)
			- ajouter le isReachable pour prendre en compte les cycle interne, peut etre tricky
		Detecter la source et sink des patterns complexes en les delimitant par les long unitigs, puis developper une methode pour trouer un chemin dedans
		Low abundant genome / gap: parvenir à lier les gaps, il y a deux type de gap:
			- gap qui n'existe pas dans les données (low abundant genome avec une partie non séquencé), Il faut bridge deux partie du graphe avec le nbSharedReads (potentiellement des tips)
			- gap dû a un cutoff min trop élevé: ici le successeur est juste masqué par le isNodeValid2
		"
		Graph::getOverlap(): pas d'overlap parfois ou dans le mauvais sens, a check
		- ToBasespace: mdbg_nodes: on en a besoin pour obtenir le nodeName des kminmer, mais dans ce MDBG, on peut enlever tous les nodes qui ne font pas parti d'un contig pour sauver de la mémoire
		- SUperbubble et COmpelx area detection: on veut enlever les branches qui n'appartiennent pas au contexte actuelle (actuellement on verifie juste que les noeuds de l'area partage des reads avec le dernier noeud visité mais ça fait qu'on est limité à une certaine taille d'area)
		- Take bubble: si on prend un chemin arbitraire dans une bubble, il faut marqué le/les autre chemins comme visité pour ne qu'ils soient détecté comme successeurs ensuite
		- Complex area: quand on skip une complex area, il faut ajouté ses nodes dans _binnedNodes pour ne pa que l'assemblage puissent partir de'un de ses unitigs
		- cutoff plus intelligent: utiliser l'abondance de tous les long unitigs pour filtrer les erreurs, plutot qu'un seul unitig, vu que l'abondance varie sur l'ensemble du genome
			- ou utiliser un filtre de base assez bas, puis filtrer on the fly, mais ça va reduire la taille des unitigs de base
		- getOverlap: on devrait bannir cette methode et toujour utiliser getSuccessors_overlap pour obtenir les infos des edge en plus des successeur
		- Complex area detection issue: quand on lance un ecoli avec 40 de coverage, il y a des pb d'assemblage a cause des erreurs,
Gros changement non testé:
	- attention on a ajouter les visitedNodes dans nodeExplored collectpossiblesuccessor (ça va etre a test) qu'on explore bien tout malgré ça"

	- cutoff onthe fly = assemblyAbundance / 2, remettre superbubble detection, isTip on he fly egalement ?"
	- dans detection superbubble classique, on pourrait jouter la gestion des cycle des la phase de cleaning
	- minimizer contig: pas besoin de recuperer/ecrire les supporting reads
	- ToBasespace: recuperer les sequenceModel et variant en une seule passes sur les data, peut etre plus encessaire avec l'approche multi k?
	- Kminmer / kmer parser: verifier qu'on a la taille de seuqnece suffisante (pour les kmer) et nb minimizers suffisant (pour kminmer) pour les générer (k+1)
	- Gfa format: créer notre propre format binaire, mais attention gfa est pratique pour debug
	- Palindrome: augmenter la densité au endroit des positions bannis, pour conserver la meme taille d'overlaps (idée seb)
	- essayer open syncmer a la place de minimzier ?
	- Superbubble qui n'a qu'une branche variante (l'autre coté va direct de la source vers sink), actuellement on supprime la branche variante mais on pourrait voir laquel garder en fonction du score checkm ou quast
	- Chimeric reads: a mettre dans readselection direct
	- graphSImplify: indexer des nodeNames plutot que des nodeIndex (_isNodeValid), car on ne conserve jamais qu'un seul des deux coté d'un nodename
	- metaflye: plein de bonne idée a check (identification superbubble, roundabout)
	- ToBasespace: correction: on peut appliquer la correction d'un kminmer dés qu'on a atteint les 20 variants, comme ça on peut free la mémoire de ces variants dés que possible 
	- Unsupported edges: on peut distinguer les fake edges (qui n'existe dans aucun reads), des bons (provenant des especes rares), en chackant l'abondance de l'edge. Les fakes edge devrait avoir une abondance extrement elevé

	REPRISE:
		- on veut detecter le plus de bubble possible
		- on veut verifier qu'un branchement est une bubble avant quoi que ce soit pour accelerer le processus d'assemblage
*/
//MetaBinner: https://www.biorxiv.org/content/10.1101/2021.07.25.453671v1


#ifndef MDBG_METAG_ASSEMBLY2
#define MDBG_METAG_ASSEMBLY2

//#define PRINT_ADDED_NODES
//#define PRINT_PREV_RANK

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
#include "toBasespace/ToBasespaceOnTheFly.hpp";


const u_int32_t LONG_UNITIG_LENGTH = 10000;



class PathExplorer{

public: 

	vector<u_int32_t> _nodePath;
	vector<u_int32_t> _prevNodes;
	GraphSimplify* _graph;
	unordered_set<u_int32_t> _visitedNodes;
	vector<UnitigData>& _unitigDatas;
	unordered_map<u_int32_t, u_int16_t> _nbVisitedTimes;
	u_int32_t _sourceAbundance;


	PathExplorer(const vector<u_int32_t>& prevNodes, unordered_set<u_int32_t>& visitedNodes, vector<UnitigData>& unitigDatas, GraphSimplify* graph, unordered_map<u_int32_t, u_int16_t>& nbVisitedTimes, u_int32_t sourceAbundance) : _unitigDatas(unitigDatas){
		_prevNodes = prevNodes;
		_visitedNodes = visitedNodes;
		_graph = graph;
		_nbVisitedTimes = nbVisitedTimes;
		_sourceAbundance = sourceAbundance;
	}

	void getNextNode(u_int32_t current_nodeIndex, bool forward, vector<u_int32_t>& nextNodes){


		nextNodes.clear();
	

		if(_graph->isEdgeNode(current_nodeIndex)){

			Unitig& unitig = _graph->nodeIndex_to_unitig(current_nodeIndex);

			u_int16_t maxVisitables = ceil(unitig._abundance / _sourceAbundance) + 1;
			if(maxVisitables < 1) maxVisitables = 1;
			
			//cout << "YO: " << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << " " << _nbVisitedTimes[current_nodeIndex] << " " << maxVisitables << endl;
			if(_nbVisitedTimes[current_nodeIndex] >= maxVisitables){
				//cout << "STAP" << endl;
				return;
			}

			
			_nbVisitedTimes[current_nodeIndex] += 1;
		}

		//if(isInInfiniteCycle(current_nodeIndex, forward)){
		//	return;
		//}

		vector<u_int32_t> successors;
		if(forward){
			_graph->getSuccessors(current_nodeIndex, 0, successors);
		}
		else{
			_graph->getPredecessors(current_nodeIndex, 0, successors);
		}



		/*
		bool needScc = false;
		for(u_int32_t nodeIndex : successors){
			vector<u_int32_t> predecessors;

			if(forward){
				_graph->getPredecessors(nodeIndex, 0, predecessors);
			}
			else{
				_graph->getSuccessors(nodeIndex, 0, predecessors);
			}

			for(u_int32_t p : predecessors){
				if(p == current_nodeIndex) continue;
				needScc = true;
			}
		}

		if(successors.size() >= 2){
			needScc = true;
		}

		if(needScc){

			Unitig& unitig = _graph->nodeIndex_to_unitig(current_nodeIndex);

        	_graph->disconnectSubGraph(unitig._index, 5000, forward);

			vector<vector<u_int32_t>> components;
			_graph->getStronglyConnectedComponents(unitig._index, forward, components);


			ofstream file_scc("/home/gats/workspace/run/overlap_test_201_multik/scc.csv");
			file_scc << "Name,Colour" << endl;

			int lala = 0;
			unordered_map<u_int32_t, vector<u_int32_t>> unitigToSccs;
			for(vector<u_int32_t>& component : components){
				for(u_int32_t unitigIndex : component){
					unitigToSccs[unitigIndex] = component;

					for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
						file_scc << BiGraph::nodeIndex_to_nodeName(nodeIndex) << "," << lala << endl;
					}
				}
				lala += 1;
			}

			vector<u_int32_t> successors2;

			for(u_int32_t nodeIndex : successors){
				unordered_set<u_int32_t> sccSuccessors;
				_graph->getSccSuccessors(unitigToSccs[_graph->nodeIndex_to_unitigIndex(nodeIndex)], forward, sccSuccessors);

				for(u_int32_t sccSuccessor : sccSuccessors){
					cout << BiGraph::nodeIndex_to_nodeName(_graph->_unitigs[sccSuccessor]._startNode) << endl;
					successors2.push_back(_graph->_unitigs[sccSuccessor]._startNode);
				}
			}

			successors.clear();
			successors = successors2;

			_graph->_isNodeInvalid_tmp.clear();
			//getchar();
			//return;
		}
		*/
		/*
		bool isBubble = true;
		for(const SuccessorData& successor : data_successors){
			if(!graph->_isBubble[successor._nodeIndexSuccessor]){
				isBubble = false;
			}
		}
		if(isBubble){
			needScc = false;
		}

		if(!needScc) return -1;
		*/





		//cout << "Succ: " << successors.size() << endl;
		if(successors.size() == 0){
			
		}
		else if(successors.size() == 1){

			nextNodes.push_back(successors[0]);
		}
		else{

			computeBestSuccessors_byUnitigRank(successors, nextNodes);
		}


		//cout << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << " " << nextNodes.size() << endl;


	}

	
	bool isInInfiniteCycle(u_int32_t current_nodeIndex, bool forward){
		/*
		//cout << "Check infinite cycle" << endl;
		Unitig& unitig = _graph->nodeIndex_to_unitig(current_nodeIndex);
		if(unitig._startNode != current_nodeIndex) return false;

		unordered_map<u_int32_t, DataSuccessorPath> successorPaths;
		collectPossibleSuccessors(current_nodeIndex, _graph, _unitigDatas, forward, successorPaths, true, -1, _visitedNodes);

		//cout << "Check infinite cycle: " << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << " " <<  successorPaths.size() << endl;
		if(successorPaths.size() == 0){
			cout << "Infinite cycle" << endl;
			//getchar();
		} 

		return successorPaths.size() == 0;
		*/
		return false;
	}


	void computeBestSuccessors_byUnitigRank(vector<u_int32_t>& successors, vector<u_int32_t>& possibleSuccessors){

		#ifdef PRINT_PREV_RANK
			cout << endl;
		#endif
		//u_int64_t longestPath = 0;
		//u_int64_t longestPathNodeIndex = -1;


		possibleSuccessors.clear();
		//for(u_int32_t nodeIndex : successors){
		//	possibleSuccessors.push_back(nodeIndex);
		//}

		u_int32_t prevRank = 0;
		size_t i = 0;

		u_int32_t lastNode = -1;
		u_int32_t pathLength = 0;


		vector<u_int32_t> removedNodeIndex;

		while(true){
			


			int prevIndex = _prevNodes.size() - prevRank - 1;
			if(prevIndex < 0 ) break;

			u_int32_t prev_nodeIndex = _prevNodes[prevIndex];
			u_int32_t prev_nodeName = BiGraph::nodeIndex_to_nodeName(prev_nodeIndex);

			if(i == 0){
				//cout << prev_nodeName << " " << _graph->_nodeLengths[prev_nodeName] << endl;
				pathLength = _graph->_nodeLengths[prev_nodeName];
			}
			else{
				u_int16_t overlapLength = _graph->_graphSuccessors->getOverlap(prev_nodeIndex, lastNode);
				pathLength += (_graph->_nodeLengths[prev_nodeName] - overlapLength);
			}

			
			#ifdef PRINT_PREV_RANK
				string str_debug = "";
				str_debug += "" + _graph->_graphSuccessors->nodeToString(prev_nodeIndex);

				for(u_int32_t& successor : successors){
					u_int32_t successor_nodeName = BiGraph::nodeIndex_to_nodeName(successor);
					u_int32_t nbSharedReads = Utils::computeSharedReads(_unitigDatas[prev_nodeName], _unitigDatas[successor_nodeName]);
					str_debug += "    " + BiGraph::nodeToString(successor) + " " + to_string(nbSharedReads);
				}
				str_debug += "    " + to_string(pathLength);
				cout << "\t" << str_debug << endl;
			#endif

			/*
			if(pathLength > 5000){
				for(u_int32_t& successor : successors){
					u_int32_t successor_nodeName = BiGraph::nodeIndex_to_nodeName(successor);
					u_int32_t nbSharedReads = Utils::computeSharedReads(_unitigDatas[prev_nodeName], _unitigDatas[successor_nodeName]);
					if(nbSharedReads >= 2){
						possibleSuccessors.push_back(successor);
					}

				}
				break;
			}
			*/

			removedNodeIndex.clear();


			if(pathLength > 9000){
				for(u_int32_t successor : successors){
					u_int32_t successor_nodeName = BiGraph::nodeIndex_to_nodeName(successor);
					u_int32_t nbSharedReads = Utils::computeSharedReads(_unitigDatas[prev_nodeName], _unitigDatas[successor_nodeName]);
					if(nbSharedReads >= 1){
						possibleSuccessors.push_back(successor);
						//longestPathNodeIndex = successor;
					}
					else{
						//removedNodeIndex.push_back(successor);
					}

				}
				break;
			}


			//for(u_int32_t successor : removedNodeIndex){
			//	possibleSuccessors.erase(std::remove(possibleSuccessors.begin(), possibleSuccessors.end(), successor), possibleSuccessors.end());
			//}

			//if(pathLength > 5000) break;

			prevRank += 1;
			i += 1;
			lastNode = prev_nodeIndex;


		}

		#ifdef PRINT_PREV_RANK
			cout << "\tNb winners: " << possibleSuccessors.size() << endl;
		#endif
		//cout << successors.size() << endl;
		//possibleSuccessors.push_back(longestPathNodeIndex);
		//exit(1);

	}



	void nodeExplored(u_int32_t nodeIndex){
		
		_nodePath.push_back(nodeIndex);
		_prevNodes.push_back(nodeIndex);
		//_exploredNodes.push_back(nodeIndex);

		PathExplorer::clampPrevNodes(_prevNodes, _unitigDatas);
		u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);


		//_currentPathLength += graph->_nodeLengths[nodeName];

		//Unitig& unitig = _graph->nodeIndex_to_unitig(nodeIndex);

		//if(nodeIndex == unitig._startNode || nodeIndex == unitig._endNode){

			//u_int16_t maxVisitables = ceil(unitig._abundance / _sourceAbundance);
			//if(maxVisitables < 1) maxVisitables = 1;
			
			//if(_nbVisitedTimes[current_nodeIndex] >= maxVisitables){
			//	return;
			//}

			
		//	_nbVisitedTimes[nodeIndex] += 1;
		//}

		//_visitedNodes.insert(nodeName);
		//_binnedNodes.insert(current_unitigIndex);
		//cout << "Node explored: " << graph->_graphSuccessors->nodeToString(nodeIndex) << " " << graph->_nodeAbundances[nodeName]  << endl;

		//_nbVisitedTimes[current_unitigIndex] += 1;
		//cout << _nbVisitedTimes[current_unitigIndex] << endl;
		//return utg_nodeIndex;
	}



	static void clampPrevNodes(vector<u_int32_t>& prevNodes, vector<UnitigData>& _unitigDatas){

		//cout << "clamp" << endl;
		if(prevNodes.size() > 10){
			u_int32_t nodeName_current = BiGraph::nodeIndex_to_nodeName(prevNodes[prevNodes.size()-1]);
			UnitigData& unitigData_current = _unitigDatas[nodeName_current];

			while(true){
				if(prevNodes.size() == 0) break;
				u_int32_t nodeIndex = prevNodes[0];
				u_int32_t nodeName_first = BiGraph::nodeIndex_to_nodeName(prevNodes[0]);

				const UnitigData& unitigData_first = _unitigDatas[nodeName_first];

				//cout << nodeName_current << " " << nodeName_first << " " << Utils::shareAnyRead(unitigData_current, unitigData_first) << " " << unitigData_current._readIndexes.size() << " " << unitigData_first._readIndexes.size() << endl;
				if(Utils::shareAnyRead(unitigData_current, unitigData_first)){
					break;
				}
				else{
					prevNodes.erase(prevNodes.begin());
				}

			}
		}

	}
};








class PathTree{

public:

	GraphSimplify* _graph;
	u_int32_t _source_nodeIndex;
	u_int32_t _source_abundance;
	vector<UnitigData>& _unitigDatas;
	string _inputDir;

	PathTree(const string& inputDir, u_int32_t source_nodeIndex, u_int32_t source_abundance, vector<UnitigData>& unitigDatas, GraphSimplify* graph) : _unitigDatas(unitigDatas){
		_inputDir = inputDir;
		_source_abundance = source_abundance;
		_source_nodeIndex = source_nodeIndex;
		_graph = graph;
	}


	void findSuccessors_unitig(u_int32_t source_unitigIndex, unordered_set<u_int32_t>& isVisited, unordered_set<u_int32_t>& longSuccessors){

		longSuccessors.clear();
		isVisited.clear();

		//u_int32_t source_unitigIndex = _graph->nodeIndex_to_unitigIndex(source_nodeIndex);

        queue<u_int32_t> queue;
        queue.push(source_unitigIndex);

		while(!queue.empty()){
			
            u_int64_t unitigIndex = queue.front();
            queue.pop();
			
			if(isVisited.find(unitigIndex) != isVisited.end()) continue;

			isVisited.insert(unitigIndex);

			const Unitig& unitig = _graph->_unitigs[unitigIndex];
			if(unitigIndex != source_unitigIndex && unitig._length > LONG_UNITIG_LENGTH){
				longSuccessors.insert(unitigIndex);
				continue;
			}

			vector<u_int32_t> successors;
			_graph->getSuccessors_unitig(unitigIndex, 0, successors);

			for(u_int32_t successor : successors){
				queue.push(successor);
			}

		}

	}

	void execute2(unordered_set<u_int32_t>& visitedNodesAll){

		ofstream file_visitedNodes(_inputDir + "/binning_results.csv");
		file_visitedNodes << "Name,Colour" << endl;

		visitedNodesAll.clear();

		unordered_set<u_int32_t> isVisited;
		u_int32_t source_unitigIndex = _graph->nodeIndex_to_unitigIndex(_source_nodeIndex);

		queue<u_int32_t> queue;
		queue.push(source_unitigIndex);

		while(!queue.empty()){

            u_int64_t unitigIndex = queue.front();
            queue.pop();


			for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
				
				if(visitedNodesAll.find(nodeName) == visitedNodesAll.end()){
					visitedNodesAll.insert(nodeName);
					file_visitedNodes << nodeName << "," << "red" << endl;
				}
			}	


			cout << "Start unitig: " << BiGraph::nodeIndex_to_nodeName(_graph->_unitigs[unitigIndex]._startNode) << endl;

			if(isVisited.find(unitigIndex) != isVisited.end()) continue;
			isVisited.insert(unitigIndex);

			unordered_set<u_int32_t> successors_unitig;
			unordered_set<u_int32_t> visitedNodes;
			findSuccessors_unitig(unitigIndex, visitedNodes, successors_unitig);

			cout << "Nb successors unitigs: " << successors_unitig.size() << endl;
			cout << "\t";
			for(u_int32_t unitigIndex : successors_unitig){
				cout << BiGraph::nodeIndex_to_nodeName(_graph->_unitigs[unitigIndex]._startNode) << " (" << (isVisited.find(unitigIndex) != isVisited.end()) << ")  ";
			}
			cout << endl;


			if(successors_unitig.size() == 0){
				cout << "No successors :(" << endl;
				continue;
			}
			else if(successors_unitig.size() == 1){
				for(u_int32_t unitigIndex : visitedNodes){
					for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
						u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
						
						if(visitedNodesAll.find(nodeName) == visitedNodesAll.end()){
							file_visitedNodes << nodeName << "," << "red" << endl;
						}
						visitedNodesAll.insert(nodeName);
					}	
				}

				for(u_int32_t successor_unitigIndex : successors_unitig){
					queue.push(successor_unitigIndex);
				}

			}
			else{
				unordered_set<u_int32_t> reachableUnitigs;
				isReachable(unitigIndex, successors_unitig, reachableUnitigs, visitedNodesAll, file_visitedNodes);

				cout << "Nb reachable unitigs: " << reachableUnitigs.size() << endl;
				cout << "\t";
				for(u_int32_t unitigIndex : reachableUnitigs){
					cout << BiGraph::nodeIndex_to_nodeName(_graph->_unitigs[unitigIndex]._startNode) << " ";
				}
				cout << endl;

				for(u_int32_t successor_unitigIndex : reachableUnitigs){
					//u_int32_t successor_nodeIndex = _graph->_unitigs[successor_unitigIndex]._startNode;

					queue.push(successor_unitigIndex);
				}
			}



		}

		//for(u_int32_t nodeName : visitedNodesAll){
		//	file_visitedNodes << nodeName << "," << "red" << endl;
		//}

		file_visitedNodes.close();

		cout << "Visited nodes: " << visitedNodesAll.size() << endl;
	}


	void isReachable(u_int32_t source_unitigIndex, unordered_set<u_int32_t>& targetUnitigIndex, unordered_set<u_int32_t>& reachableUnitigIndex, unordered_set<u_int32_t>& visitedNodesAll, ofstream& file_visitedNodes){

		reachableUnitigIndex.clear();
		unordered_set<u_int32_t> visitedNodesLocal;


		u_int32_t source_nodeIndex = _graph->_unitigs[source_unitigIndex]._startNode;

		//u_int32_t current_nodeIndex = _source_nodeIndex;

		bool forward = true;
		vector<PathExplorer> queue;

		vector<u_int32_t> prevNodes;
		unordered_set<u_int32_t> visitedNodes;
		unordered_map<u_int32_t, u_int16_t> nbVisitedTimes;
		unordered_set<u_int32_t> blockedNodes;

		PathExplorer pathExplorer_source (prevNodes, visitedNodes, _unitigDatas, _graph, nbVisitedTimes, _source_abundance);
		pathExplorer_source.nodeExplored(source_nodeIndex);

		queue.push_back(pathExplorer_source);

		u_int64_t nbIters = 0;
		u_int64_t prevVisited = 0;

		while(!queue.empty()){

			if(reachableUnitigIndex.size() >= targetUnitigIndex.size()) break;

			PathExplorer pathExplorer = queue[queue.size()-1];
			queue.pop_back();
			//cout << queue.size() << endl;
			u_int32_t current_nodeIndex = pathExplorer._prevNodes[pathExplorer._prevNodes.size()-1];
			
			u_int32_t current_unitigIndex = _graph->nodeIndex_to_unitigIndex(current_nodeIndex);
			if(targetUnitigIndex.find(current_unitigIndex) != targetUnitigIndex.end()){

				if(reachableUnitigIndex.find(current_unitigIndex) == reachableUnitigIndex.end()){
					reachableUnitigIndex.insert(current_unitigIndex);
					cout << "Reached: " << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << endl;
					
					
					for(u_int32_t nodeIndex : pathExplorer._nodePath){
						u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
						
						if(visitedNodesAll.find(nodeName) == visitedNodesAll.end()){
							file_visitedNodes << nodeName << "," << "red" << endl;
						}
						visitedNodesAll.insert(nodeName);
					}

				}

				continue;
			}


			if(current_unitigIndex != source_unitigIndex && _graph->_unitigs[current_unitigIndex]._length > LONG_UNITIG_LENGTH) continue;

			//cout << BiGraph::nodeToString(current_nodeIndex) << endl;

			/*
			if(_graph->isEdgeNode(current_nodeIndex)){
				
				Unitig& unitig = _graph->nodeIndex_to_unitig(current_nodeIndex);
				if(unitig._length > 20000){
					if(visitedNodes.find(current_nodeIndex) != visitedNodes.end()) continue;
					visitedNodes.insert(current_nodeIndex);
				}
			}
			*/
			
			cout << blockedNodes.size() << " " << visitedNodesAll.size() << endl;
			if(blockedNodes.find(BiGraph::nodeIndex_to_nodeName(current_nodeIndex)) != blockedNodes.end()) continue;

			if(prevVisited == visitedNodesLocal.size()){
				nbIters += 1;
				if(nbIters > 1000){
					blockedNodes.insert(BiGraph::nodeIndex_to_nodeName(current_nodeIndex));
					nbIters = 0;
				} 
				//cout << nbIters << endl;
			}
			else{
				prevVisited = visitedNodesLocal.size();
				nbIters = 0;
			}


			while(true){

				vector<u_int32_t> nextNodes;
				pathExplorer.getNextNode(current_nodeIndex, forward, nextNodes);

				//cout << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << " " << nextNodes.size() << endl;

				if(nextNodes.size() == 0){
					break;
				}
				else if(nextNodes.size() == 1){



					visitedNodesLocal.insert(BiGraph::nodeIndex_to_nodeName(current_nodeIndex));

					u_int32_t current_unitigIndex = _graph->nodeIndex_to_unitigIndex(current_nodeIndex);
					if(targetUnitigIndex.find(current_unitigIndex) != targetUnitigIndex.end()){

						if(reachableUnitigIndex.find(current_unitigIndex) == reachableUnitigIndex.end()){
							reachableUnitigIndex.insert(current_unitigIndex);
							cout << "Reached: " << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << endl;
							
							
							for(u_int32_t nodeIndex : pathExplorer._nodePath){
								u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
								
								if(visitedNodesAll.find(nodeName) == visitedNodesAll.end()){
									file_visitedNodes << nodeName << "," << "red" << endl;
								}
								visitedNodesAll.insert(nodeName);
							}
							
						}

						break;
					}

					if(current_unitigIndex != source_unitigIndex && _graph->_unitigs[current_unitigIndex]._length > LONG_UNITIG_LENGTH) break;

					current_nodeIndex = nextNodes[0];
					pathExplorer.nodeExplored(current_nodeIndex);
					
					
					//if(_graph->isEdgeNode(current_nodeIndex)){
					//
					//	Unitig& unitig = _graph->nodeIndex_to_unitig(current_nodeIndex);
					//	if(unitig._length > 10000){
					//		if(visitedNodes.find(current_nodeIndex) != visitedNodes.end()) break;
					//		visitedNodes.insert(current_nodeIndex);
					//	}
					//}
					

					//if(visitedNodes.find(current_nodeIndex) != visitedNodes.end()) break;

					//if(visitedNodesAll.find(BiGraph::nodeIndex_to_nodeName(current_nodeIndex)) == visitedNodesAll.end()){
					//	file_visitedNodes << BiGraph::nodeIndex_to_nodeName(current_nodeIndex) << "," << "red" << endl;
					//}

					visitedNodesLocal.insert(BiGraph::nodeIndex_to_nodeName(current_nodeIndex));
					
					#ifdef PRINT_ADDED_NODES
						cout << "\tAdd node: " << BiGraph::nodeToString(current_nodeIndex) << endl;
					#endif


					//if(current_nodeIndex == _source_nodeIndex) break;
				}
				else{
					for(u_int32_t nodeIndex : nextNodes){

						
						#ifdef PRINT_ADDED_NODES
							cout << "\tAdd node (branch): " << BiGraph::nodeToString(nodeIndex) << endl;
						#endif

						//if(_graph->isEdgeNode(nodeIndex)){
					
						//	Unitig& unitig = _graph->nodeIndex_to_unitig(nodeIndex);
						//	if(unitig._length > 10000){
						//		if(visitedNodes.find(nodeIndex) != visitedNodes.end()) continue;
						//		visitedNodes.insert(nodeIndex);
						//	}
						//}
						


						PathExplorer pathExplorerNext = pathExplorer;
						pathExplorerNext.nodeExplored(nodeIndex);
						
						//if(visitedNodesAll.find(BiGraph::nodeIndex_to_nodeName(nodeIndex)) == visitedNodesAll.end()){
						//	file_visitedNodes << BiGraph::nodeIndex_to_nodeName(nodeIndex) << "," << "red" << endl;
						//}
						visitedNodesLocal.insert(BiGraph::nodeIndex_to_nodeName(nodeIndex));

						//if(visitedNodes.find(nodeIndex) != visitedNodes.end()) continue;
						//if(nodeIndex == _source_nodeIndex) continue;

						//visitedNodes.insert(nodeIndex);

						queue.push_back(pathExplorerNext);


					}

					break;
				}

			}
		}

		//cout << "Visited nodes: " << visitedNodesAll.size() << endl;

		

		//for(u_int32_t nodeName : visitedNodesAll){
		//	file_visitedNodes << nodeName << "," << "red" << endl;
		//}

		//file_visitedNodes.close();
	}

};




































































class Assembly2 : public Tool{

public:

	string _inputDir;
	string _truthInputFilename;
	string _outputFilename;
	string _outputFilename_complete;
	bool _debug;
	string _inputFilename_unitigNt;
	string _inputFilename_unitigCluster;

	float _minimizerDensity;
    size_t _minimizerSize;
    size_t _kminmerSize;
	vector<UnitigData> _unitigDatas;

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
	ToBasespaceOnTheFly _toBasespace;
	string _gfaFilename;

	Assembly2(): Tool (){

	}

	~Assembly2(){

	}

	void execute (){

		_outputContigFile = gzopen(_outputFilename.c_str(),"wb");
		_outputContigFile_complete = gzopen(_outputFilename_complete.c_str(),"wb");

		loadGraph();
		execute_detectSpecies_byCutoff();
		
		gzclose(_outputContigFile);
		gzclose(_outputContigFile_complete);
	}

	void parseArgs(int argc, char* argv[]){



		cxxopts::Options options("Assembly", "");
		options.add_options()
		(ARG_OUTPUT_DIR, "", cxxopts::value<string>())
		(ARG_INPUT_FILENAME_CONTIG, "", cxxopts::value<string>()->default_value(""))
		(ARG_INPUT_FILENAME_TRUTH, "", cxxopts::value<string>()->default_value(""))
		(ARG_DEBUG, "", cxxopts::value<bool>()->default_value("false"))
		(ARG_INPUT_FILENAME_UNITIG_NT, "", cxxopts::value<string>()->default_value(""))
		(ARG_INPUT_FILENAME_UNITIG_CLUSTER, "", cxxopts::value<string>()->default_value(""));



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
		_filename_readMinimizers = _inputDir + "/read_data.gz";
		_filename_hifiasmGroundtruth = _inputDir + "/hifiasmGroundtruth.gz";
		_filename_outputContigs = _inputDir + "/contigs.min.gz";
		_filename_solidKminmers = _inputDir + "/solid.min.gz";

	}

	void loadGraph(){

		_gfaFilename = _inputDir + "/minimizer_graph.gfa";
		string gfa_filename_noUnsupportedEdges = _inputDir + "/minimizer_graph_noUnsupportedEdges.gfa";
		//string gfa_filename_unitigs = _inputDir + "/minimizer_graph_unitigs.gfa";
		string mdbg_filename = _inputDir + "/mdbg_nodes.gz";




		
		cout << _gfaFilename << endl;
		_mdbg = new MDBG(_kminmerSize);
		_mdbg->load(mdbg_filename);
		cout << "Nb nodes: " <<  _mdbg->_dbg_nodes.size() << endl;


		if(_truthInputFilename != ""){
			extract_truth_kminmers();
		}


		file_groundTruth = ofstream(_inputDir + "/binning_results.csv");
		file_groundTruth << "Name,Colour" << endl;

		file_groundTruth_hifiasmContigs = ofstream(_inputDir + "/binning_results_hifiasm_" + to_string(_kminmerSize) + ".csv");
		file_groundTruth_hifiasmContigs << "Name,Colour" << endl;

		//if(_debug){
            //gfa_filename = _inputDir + "/minimizer_graph_debug.gfa";
		//}
		
		GraphSimplify* graphSimplify = new GraphSimplify(_gfaFilename, _inputDir, 0);
		_graph = graphSimplify;
		
		
		//Generate unitigs
		cout << "Indexing reads" << endl;
		_unitigDatas.resize(_mdbg->_dbg_nodes.size());
		_graph->clear(0);
		_graph->compact(false, _unitigDatas);
		removeUnsupportedEdges(_gfaFilename, gfa_filename_noUnsupportedEdges, _graph);
		delete _mdbg;
		cout << "done" << endl;
	


		//cout << graphSimplify->_graphSuccessors->_nbEdges << endl;
		//graphSimplify->execute(5, _kminmerSize);
		//graphSimplify->debug_writeGfaErrorfree(1000, PathExplorer::computeAbundanceCutoff(1000, 0, CutoffType::ERROR), -1, _kminmerSize, false, true, false, _unitigDatas);
		_graph->debug_writeGfaErrorfree(500, 500, -1, _kminmerSize, false, true, false, _unitigDatas);




		_toBasespace.create(_inputDir);


	}

	void execute_detectSpecies_byCutoff(){



		//size_t globalAbundanceCutoff_min = 3;
		
		

		//detectIsolatedSpecies(gfa_filename, toBasespace);
		//exit(1);
		u_int64_t nbUnitigs = 0;
		unordered_set<u_int32_t> writtenNodeNames;
		unordered_map<u_int32_t, u_int32_t> nodeName_to_unitigName;
		

		string clusterDir = _inputDir + "/" + "binByCutoff";
		fs::path path(clusterDir);
	    if(!fs::exists (path)){
            fs::create_directory(path);
        } 

		ofstream fileHifiasmAll(clusterDir + "/component_hifiasm_all.csv");
		fileHifiasmAll << "Name,Colour" << endl;
		ofstream fileComponentNodeAll(clusterDir + "/component_node_all.csv");
		fileComponentNodeAll << "Name,Colour" << endl;

		u_int32_t clusterIndex = 0;
		u_int32_t contaminatedIndex = 0;
		unordered_set<string> written_unitigName;
		float prevCutoff = -1;

		vector<vector<u_int32_t>> clusters;
		vector<vector<u_int32_t>> clustersValid;
		bool reloadState = true;


		for(Unitig& unitig : _graph->_startingUnitigstest){

			
			u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(unitig._startNode);

			/*
			bool isValid = false;
			for(u_int32_t nodeIndex : unitig._nodes){
				if(_hifiasm_startingNodenames.find(BiGraph::nodeIndex_to_nodeName(nodeIndex)) != _hifiasm_startingNodenames.end()){
					isValid = true;
				}
			}
			if(!isValid) continue;
			*/

			float cutoff = (unitig._abundance*0.9)*0.5;


			float realCutoff = 0;
			for(const SaveState2& saveState : _graph->_cachedGraphStates){
				if(saveState._abundanceCutoff_min > cutoff) break;
				realCutoff = saveState._abundanceCutoff_min;
			}


			if(realCutoff != prevCutoff || reloadState){
				_graph->loadState2(cutoff, unitig._startNode, _unitigDatas);
				prevCutoff = realCutoff;
				reloadState = false;
			}
			

			//float cutoff = unitig._abundance*0.1;
			//_graph->loadState2(cutoff, unitig._startNode, _unitigDatas);

        	//unordered_set<u_int32_t> component;
			//for(u_int32_t nodeIndex : _graph->_isNodeValid2){
			//	component.insert(_graph->nodeIndex_to_unitigIndex(nodeIndex));
			//}

        	unordered_set<u_int32_t> component;
        	_graph->getConnectedComponent(unitig._startNode, component);

			if(component.size() <= 1) continue;

			vector<u_int32_t> cluster;
			//unordered_set<u_int32_t> validNodes;

			for(u_int32_t unitigIndex : component){
				for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
					u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
					cluster.push_back(nodeName);
					//validNodes.insert(nodeName);
				}
			}

			bool isNewCluster = true;
			std::sort(cluster.begin(), cluster.end());

			cout << endl << "Unitig: " << unitig._index << endl;

			for(const vector<u_int32_t>& existingCluster : clusters){
				float jaccardDistance = Utils::computeJaccardDistance(existingCluster, cluster);
				cout << "Existing: " << jaccardDistance << endl;
				//float sharedRate = ((float) nbSharedNodes) / ((float)cluster.size());
				//cout << sharedRate << endl;
				if(jaccardDistance < 0.33){
					isNewCluster = false;
					break;
				}
			}
			
			for(const vector<u_int32_t>& existingCluster : clustersValid){
				float jaccardDistance = Utils::computeJaccardDistance(existingCluster, cluster);
				cout << "Valid: " << jaccardDistance << endl;
				//float sharedRate = ((float) nbSharedNodes) / ((float)cluster.size());
				//cout << sharedRate << endl;
				if(jaccardDistance < 0.6){
					isNewCluster = false;
					break;
				}
			}

			if(!isNewCluster) continue;
		
			//cout << clusterIndex << " " << component.size() << endl;
			
			/*
			if(validNodes.size() > 0){
				
				reloadState = true;

				cout << endl << "Start new cluster: " << endl;
				
				const string& filenameComponentNodes = clusterDir + "/component_" + to_string(clusterIndex) + "_nodes.csv";
				ofstream fileComponenetNodes(filenameComponentNodes);
				fileComponenetNodes << "Name,Colour" << endl;
				for(u_int32_t nodeName : visitedNodes){
					fileComponenetNodes << nodeName << "," << clusterIndex << endl;
				}
				fileComponenetNodes.close();
				
				//exit(1);
				clusters.push_back(cluster);
			*/

			clusters.push_back(cluster);
			reloadState = true;

			bool hasValidContamination = false;
			float componentCompletneness = 0;
			vector<u_int32_t> componentNodesValid;


			//for(float cutoffRate : {0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1}){			
			for(float cutoffRate : {0.5, 0.4, 0.3, 0.2, 0.1}){
					
				cout << "Cutoff rate: " << cutoffRate << endl;

				float cutoff = (unitig._abundance*0.9)*cutoffRate;
				_graph->loadState2(cutoff, unitig._startNode, _unitigDatas);


				unordered_set<u_int32_t> component;
				_graph->getConnectedComponent(unitig._startNode, component);

				//unordered_set<u_int32_t> component;
				//for(u_int32_t nodeIndex : _graph->_isNodeValid2){
				//	component.insert(_graph->nodeIndex_to_unitigIndex(nodeIndex));
				//}

				if(component.size() <= 1) continue;





				vector<u_int32_t> cluster;
				//unordered_set<u_int32_t> validNodes;

				for(u_int32_t unitigIndex : component){
					for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
						u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
						cluster.push_back(nodeName);
						//validNodes.insert(nodeName);
					}
				}

				std::sort(cluster.begin(), cluster.end());

				bool isNewCluster = true;
				for(const vector<u_int32_t>& existingCluster : clustersValid){
					float jaccardDistance = Utils::computeJaccardDistance(existingCluster, cluster);
					cout << jaccardDistance << endl;
					//float sharedRate = ((float) nbSharedNodes) / ((float)cluster.size());
					//cout << sharedRate << endl;
					if(jaccardDistance < 0.6){
						isNewCluster = false;
						break;
					}
				}
				if(!isNewCluster) break;


				/*
				cout << "subgraph extract" << endl;
				//_graph->compact(false, _unitigDatas);
				_graph->extractReadpathSubgraph(_graph->nodeIndex_to_unitigIndex(BiGraph::nodeName_to_nodeIndex(963565, true)));
				exit(1);
				*/

				const string& filenameUnitigColor= clusterDir + "/component_" + to_string(clusterIndex) + "_unitigColor.csv";
				ofstream fUnitigColor(filenameUnitigColor);
				fUnitigColor << "Name,Colour" << endl;

				cout << "Generate contigs" << endl;
				const string& basespaceUnitigFilename = clusterDir + "/component_" + to_string(clusterIndex) + "_unitigs.fasta";
				//gzFile basespaceUnitigFile = gzopen(basespaceUnitigFilename.c_str(),"wb");
				ofstream basespaceUnitigFile = ofstream(basespaceUnitigFilename);

				vector<u_int32_t> componentNodes;

				//u_int32_t contigIndex = 0;
				unordered_set<u_int32_t> writtenUnitigs;
				for(u_int32_t unitigIndex : component){
					
					const Unitig& u = _graph->_unitigs[unitigIndex];
					u_int32_t contigIndex = u._index;
					//if(u._length < 10000) continue;
					
					if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._startNode)) != writtenUnitigs.end()) continue;
					if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._endNode)) != writtenUnitigs.end()) continue;

					writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._startNode));
					writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._endNode));

					string unitigSequence;
					_toBasespace.createSequence(u._nodes, unitigSequence);

					string header = ">ctg" + to_string(contigIndex) + '\n';
					basespaceUnitigFile << header;
					//gzwrite(basespaceUnitigFile, (const char*)&header[0], header.size());
					unitigSequence +=  '\n';
					basespaceUnitigFile << unitigSequence;
					//gzwrite(basespaceUnitigFile, (const char*)&unitigSequence[0], unitigSequence.size());
					
					for(u_int32_t nodeIndex : u._nodes){
						u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
						componentNodes.push_back(nodeName);
						fUnitigColor << nodeName << "," << contigIndex << endl;
					}

					//cout << unitigSequence.size() << endl;

					contigIndex += 1;
				}

				fUnitigColor.close();
				//gzclose(basespaceUnitigFile);
				basespaceUnitigFile.close();

				
				string command_annotation = "python3 /home/gats/workspace/tools/binner/thirdparty/scg/metacoag_main.py " + basespaceUnitigFilename + " " + filenameUnitigColor + " --isSingleSpecies";
				cout << command_annotation << endl;
				int ret = system(command_annotation.c_str());
				if(ret != 0){
					cerr << "Command failed: " << ret << endl;
					exit(ret);
				}
				
				//char isSingleSpeciesResult;
				//ifstream fileIsSingleSpecies(basespaceUnitigFilename + "_isSingleSpecies.txt");
				//fileIsSingleSpecies.get(isSingleSpeciesResult);
				//fileIsSingleSpecies.close();

				//cout << isSingleSpeciesResult << endl;


				indexScgClusters(basespaceUnitigFilename);
				float contamination;
				float completeness;
				computeContamination(component, completeness, contamination);
				cout << "Completeness: " << completeness << endl;
				cout << "Contamination: " << contamination << endl;

				//computeScgCentrality(clusterDir + "/contaminated_" + to_string(contaminatedIndex) + ".gfa.centrality.csv", component);
				//exit(1);

				//if(contamination <= 0.01){



				if(contamination < 0.05){

					/*
					vector<vector<u_int32_t>> components;
					_graph->getStronglyConnectedComponents(unitig._index, true, components);
					ofstream file_scc1(clusterDir + "/component_" + to_string(clusterIndex) + ".gfa" + ".scc_1.csv");
					file_scc1 << "Name,Colour" << endl;
					int sccCluster = 0;
					for(vector<u_int32_t>& component : components){
						for(u_int32_t unitigIndex : component){
							for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
								file_scc1 << BiGraph::nodeIndex_to_nodeName(nodeIndex) << "," << sccCluster << endl;
							}
						}
						sccCluster += 1;
					}
					file_scc1.close();

					_graph->getStronglyConnectedComponents(_graph->unitigIndex_toReverseDirection(unitig._index), true, components);
					ofstream file_scc2(clusterDir + "/component_" + to_string(clusterIndex) + ".gfa" + ".scc_2.csv");
					file_scc2 << "Name,Colour" << endl;
					sccCluster = 0;
					for(vector<u_int32_t>& component : components){
						for(u_int32_t unitigIndex : component){
							for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
								file_scc2 << BiGraph::nodeIndex_to_nodeName(nodeIndex) << "," << sccCluster << endl;
							}
						}
						sccCluster += 1;
					}
					file_scc2.close();
					*/

					componentNodesValid = componentNodes;
					componentCompletneness = completeness;
					hasValidContamination = true;
				}
				else{

					/*
					//if(contamination > 0.7){

						unordered_set<u_int32_t> validNodes;
						for(u_int32_t unitigIndex : component){
							for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
								validNodes.insert(BiGraph::nodeIndex_to_nodeName(nodeIndex));
							}
						}

						string outputFilename = clusterDir + "/contaminated_" + to_string(contaminatedIndex) + ".gfa";
						GfaParser::rewriteGfa_withoutNodes(_gfaFilename, outputFilename, validNodes, _graph->_isEdgeRemoved, _graph->_graphSuccessors);
			
						string command_annotation = "python3 /home/gats/workspace/tools/binner/thirdparty/scg/metacoag_main.py " + basespaceUnitigFilename + " " + filenameUnitigColor;
						cout << command_annotation << endl;
						int ret = system(command_annotation.c_str());
						if(ret != 0){
							cerr << "Command failed: " << ret << endl;
							exit(ret);
						}

						indexScgClusters(basespaceUnitigFilename);
						computeScgCentrality(outputFilename, component);

						contaminatedIndex += 1;
					//}

					break;
					*/
				}

			}

			//if(!hasValidContamination){
			//	clusters.push_back(cluster);
			//}

			if(hasValidContamination && componentCompletneness > 0.5){

				

				std::sort(componentNodesValid.begin(), componentNodesValid.end());
				clustersValid.push_back(componentNodesValid);



				cout << "Found valid component: " << componentCompletneness << endl;

				const string& filenameHifiasm = clusterDir + "/component_" + to_string(clusterIndex) + "_hifiasm.csv";
				unordered_set<string> hifiasmWrittenUnitigs;
				ofstream fileHifiasm(filenameHifiasm);
				fileHifiasm << "Name,Colour" << endl;

				unordered_set<u_int32_t> componentNodes;
				for(u_int32_t nodeName : componentNodesValid){
					componentNodes.insert(nodeName);

					fileComponentNodeAll << nodeName << "," << clusterIndex << endl;

					if(_truthInputFilename != ""){
						

						if(_evaluation_hifiasmGroundTruth_nodeName_to_unitigName.find(nodeName) != _evaluation_hifiasmGroundTruth_nodeName_to_unitigName.end()){
							for(string& unitigName : _evaluation_hifiasmGroundTruth_nodeName_to_unitigName[nodeName]){
								if(hifiasmWrittenUnitigs.find(unitigName) != hifiasmWrittenUnitigs.end()) continue;
								fileHifiasm << unitigName << "," << clusterIndex << endl;
								fileHifiasmAll << unitigName << "," << clusterIndex << endl;
								hifiasmWrittenUnitigs.insert(unitigName);
							}
						}

					}

				}

				fileHifiasm.close();

				string outputFilename = clusterDir + "/component_" + to_string(clusterIndex) + ".gfa";
				GfaParser::rewriteGfa_withoutNodes(_gfaFilename, outputFilename, componentNodes, _graph->_isEdgeRemoved, _graph->_graphSuccessors);
		




				clusterIndex += 1;
				//getchar();
			}
			else{
				

			}






				
			//}

			//cout << "cluster done" << endl;
			
		}





		//getchar();

		file_groundTruth.close();
		file_groundTruth_hifiasmContigs.close();
		file_kminmersContigs.close();

		fileHifiasmAll.close();
		fileComponentNodeAll.close();

	}

	struct UnitigScgCluster{
		u_int32_t _unitigIndex;
		u_int32_t _scgIndex;
		u_int32_t _scgCluster;
	};

	struct ScgCluster{
		u_int32_t _scgIndex;
		u_int32_t _scgCluster;
	};

	unordered_map<u_int32_t, vector<ScgCluster>> _scgClusters;

	void indexScgClusters(const string& contigFilename){

		_scgClusters.clear();

        vector<string>* fields = new vector<string>();

    	const string& scgClusterFilename = contigFilename + ".scgCluster.txt";

		string line;
		ifstream infile(scgClusterFilename);

		while (std::getline(infile, line)){
			
			GfaParser::tokenize(line, fields, '\t');

			u_int32_t unitigIndex = stoull((*fields)[0]);
			u_int32_t scgIndex = stoull((*fields)[1]);
			u_int32_t clusterIndex = stoull((*fields)[2]);

			_scgClusters[unitigIndex].push_back({scgIndex, clusterIndex});

			//cout << unitigIndex << " " << scgIndex << " " << clusterIndex << endl;
		}

		delete fields;
	}

	void computeContamination(unordered_set<u_int32_t>& unitigs, float& completeness, float& contamination){

		completeness = 0;
		contamination = 1;

		unordered_map<u_int32_t, vector<u_int32_t>> scgCounts;

		for(u_int32_t unitigIndex : unitigs){

			if(_scgClusters.find(unitigIndex) == _scgClusters.end()) continue;

			const vector<ScgCluster>& scgClusters = _scgClusters[unitigIndex];

			for(const ScgCluster& scgCluster : scgClusters){

				if(scgCounts.find(scgCluster._scgIndex) == scgCounts.end()){
					scgCounts[scgCluster._scgIndex].push_back(scgCluster._scgCluster);
				}
				else{
					vector<u_int32_t>& clusters = scgCounts[scgCluster._scgIndex];
					if(std::find(clusters.begin(), clusters.end(), scgCluster._scgCluster) == clusters.end()){
						clusters.push_back(scgCluster._scgCluster);
					}
				}
			}
		}

		u_int32_t nbContaminatedScg = 0;
		for(auto& it : scgCounts){

			u_int32_t scgIndex = it.first;
			const vector<u_int32_t>& scgCluters = it.second;

			if(scgCluters.size() > 1){
				nbContaminatedScg += 1;
			}

			cout << scgIndex << " " << scgCluters.size() << endl;
		}

		completeness = ((float) scgCounts.size()) / ((float) 107);
		contamination = ((float) nbContaminatedScg) / ((float) 107); //((float) scgCounts.size());
	}

	void execute_assembly(){


		u_int64_t nbUnitigs = 0;
		unordered_set<u_int32_t> writtenNodeNames;
		unordered_map<u_int32_t, u_int32_t> nodeName_to_unitigName;



		//detectIsolatedSpecies(gfa_filename, toBasespace);
		//exit(1);


		string clusterDir = _inputDir + "/" + "components";
		fs::path path(clusterDir);
	    if(!fs::exists (path)){
            fs::create_directory(path);
        } 

		u_int32_t clusterIndex = 0;
		unordered_set<string> written_unitigName;
		float prevCutoff = -1;

		vector<vector<u_int32_t>> clusters;

		for(Unitig& unitig : _graph->_startingUnitigstest){

			u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(unitig._startNode);
			//if(nodeName_to_cluster.find(nodeName) != nodeName_to_cluster.end()) continue;

			float cutoff = unitig._abundance*0.1;

			//float minUnitigAbundance = unitig._abundance*0.8;

			float realCutoff = 0;
			for(const SaveState2& saveState : _graph->_cachedGraphStates){
				//cout << saveState._abundanceCutoff_min << endl;
				if(saveState._abundanceCutoff_min > cutoff) break;
				realCutoff = saveState._abundanceCutoff_min;
			}


			if(realCutoff != prevCutoff){
				_graph->loadState2(cutoff, unitig._startNode, _unitigDatas);
				prevCutoff = realCutoff;
			}

        	//unordered_set<u_int32_t> component;
			//for(u_int32_t nodeIndex : _graph->_isNodeValid2){
			//	component.insert(_graph->nodeIndex_to_unitigIndex(nodeIndex));
			//}

        	unordered_set<u_int32_t> component;
        	_graph->getConnectedComponent(unitig._startNode, component);

			if(component.size() <= 1) continue;

			vector<u_int32_t> cluster;
			unordered_set<u_int32_t> validNodes;
			//unordered_set<u_int32_t> clusterUnitigNames;

			for(u_int32_t unitigIndex : component){
				
				//if(graphSimplify->_unitigs[unitigIndex]._length < 10000) continue;
				//if(graphSimplify->_unitigs[unitigIndex]._abundance < minUnitigAbundance) continue;

				for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
					u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
					//if(nodeName_to_cluster.find(nodeName) != nodeName_to_cluster.end()) continue;
					cluster.push_back(nodeName);
					//nodeName_to_cluster[nodeName] = clusterIndex;
					validNodes.insert(nodeName);


				}
			}

			bool isNewCluster = true;
			std::sort(cluster.begin(), cluster.end());

			//if(clusters.size() == 0){
				//isNewCluster = true;
			//}
			//else{
			for(const vector<u_int32_t>& existingCluster : clusters){
				float jaccardDistance = Utils::computeJaccardDistance(existingCluster, cluster);
				cout << jaccardDistance << endl;
				//float sharedRate = ((float) nbSharedNodes) / ((float)cluster.size());
				//cout << sharedRate << endl;
				if(jaccardDistance < 0.1){
					isNewCluster = false;
					break;
				}
			}
			//}
			if(!isNewCluster) continue;

			cout << clusterIndex << " " << component.size() << endl;
			
			if(component.size() > 50 && validNodes.size() > 0){

				cout << "Start tree exploration" << endl;
				unordered_set<u_int32_t> visitedNodes;
				PathTree tree(_inputDir, unitig._startNode, unitig._abundance, _unitigDatas, _graph);
				tree.execute2(visitedNodes);
				//exit(1);

				const string& filenameComponentNodes = clusterDir + "/component_" + to_string(clusterIndex) + "_nodes.csv";
				ofstream fileComponenetNodes(filenameComponentNodes);
				fileComponenetNodes << "Name,Colour" << endl;
				for(u_int32_t nodeName : visitedNodes){
					fileComponenetNodes << nodeName << "," << clusterIndex << endl;
				}
				fileComponenetNodes.close();
				
				//exit(1);
				clusters.push_back(cluster);
				
				string outputFilename = clusterDir + "/component_" + to_string(clusterIndex) + ".gfa";
				GfaParser::rewriteGfa_withoutNodes(_gfaFilename, outputFilename, validNodes, _graph->_isEdgeRemoved, _graph->_graphSuccessors);
			



				const string& filenameUnitigColor= clusterDir + "/component_" + to_string(clusterIndex) + "_unitigColor.csv";
				ofstream fUnitigColor(filenameUnitigColor);
				fUnitigColor << "Name,Colour" << endl;

				cout << "Generate contigs" << endl;
				const string& basespaceUnitigFilename = clusterDir + "/component_" + to_string(clusterIndex) + "_unitigs.fasta";
				//gzFile basespaceUnitigFile = gzopen(basespaceUnitigFilename.c_str(),"wb");
				ofstream basespaceUnitigFile = ofstream(basespaceUnitigFilename);

				u_int32_t contigIndex = 0;
				unordered_set<u_int32_t> writtenUnitigs;
				for(u_int32_t unitigIndex :component ){
					
					const Unitig& u = _graph->_unitigs[unitigIndex];
					
					//if(u._length < 10000) continue;
					
					if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._startNode)) != writtenUnitigs.end()) continue;
					if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._endNode)) != writtenUnitigs.end()) continue;

					writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._startNode));
					writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._endNode));

					string unitigSequence;
					_toBasespace.createSequence(u._nodes, unitigSequence);

					string header = ">ctg" + to_string(contigIndex) + '\n';
					basespaceUnitigFile << header;
					//gzwrite(basespaceUnitigFile, (const char*)&header[0], header.size());
					unitigSequence +=  '\n';
					basespaceUnitigFile << unitigSequence;
					//gzwrite(basespaceUnitigFile, (const char*)&unitigSequence[0], unitigSequence.size());
					
					for(u_int32_t nodeIndex : u._nodes){
						u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
						fUnitigColor << nodeName << "," << contigIndex << endl;
					}



					contigIndex += 1;
				}

				fUnitigColor.close();
				//gzclose(basespaceUnitigFile);
				basespaceUnitigFile.close();

				
				string command_annotation = "python3 /home/gats/workspace/tools/binner/src/scg/metacoag_main.py " + basespaceUnitigFilename + " " + filenameUnitigColor;
				cout << command_annotation << endl;
				int ret = system(command_annotation.c_str());
				if(ret != 0){
					cerr << "Command failed: " << ret << endl;
					exit(ret);
				}

			}

			clusterIndex += 1;
			cout << "cluster done" << endl;
			getchar();
		}





		//getchar();

		file_groundTruth.close();
		file_groundTruth_hifiasmContigs.close();
		file_kminmersContigs.close();



	}

	void detectIsolatedSpecies(){

		u_int64_t nbUnitigs = 0;
		unordered_set<u_int32_t> writtenNodeNames;

		string clusterDir = _inputDir + "/" + "components";
		fs::path pathComponents(clusterDir);
	    if(!fs::exists (pathComponents)){
            fs::create_directory(pathComponents);
        } 

		string singleSpeciesDir = _inputDir + "/" + "singleSpecies";
		fs::path pathSingleSpecies(singleSpeciesDir);
	    if(!fs::exists (pathSingleSpecies)){
            fs::create_directory(pathSingleSpecies);
        } 

		u_int32_t clusterIndex = 0;
		unordered_set<string> written_unitigName;
		float prevCutoff = -1;

		vector<vector<u_int32_t>> clusters;

		for(Unitig& unitig : _graph->_startingUnitigstest){

			u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(unitig._startNode);
			//if(nodeName_to_cluster.find(nodeName) != nodeName_to_cluster.end()) continue;

			float cutoff = unitig._abundance*0.1;

			//float minUnitigAbundance = unitig._abundance*0.8;

			float realCutoff = 0;
			for(const SaveState2& saveState : _graph->_cachedGraphStates){
				//cout << saveState._abundanceCutoff_min << endl;
				if(saveState._abundanceCutoff_min > cutoff) break;
				realCutoff = saveState._abundanceCutoff_min;
			}


			if(realCutoff != prevCutoff){
				_graph->loadState2(cutoff, unitig._startNode, _unitigDatas);
				prevCutoff = realCutoff;
			}

			if(_graph->_isNodeValid2.find(unitig._startNode) == _graph->_isNodeValid2.end()) continue;

        	unordered_set<u_int32_t> component;
        	_graph->getConnectedComponent(unitig._startNode, component);

			if(component.size() <= 1) continue;

			vector<u_int32_t> cluster;
			unordered_set<u_int32_t> validNodes;
			//unordered_set<u_int32_t> clusterUnitigNames;

			for(u_int32_t unitigIndex : component){
				
				//if(graphSimplify->_unitigs[unitigIndex]._length < 10000) continue;
				//if(graphSimplify->_unitigs[unitigIndex]._abundance < minUnitigAbundance) continue;

				for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
					u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
					//if(nodeName_to_cluster.find(nodeName) != nodeName_to_cluster.end()) continue;
					cluster.push_back(nodeName);
					//nodeName_to_cluster[nodeName] = clusterIndex;
					validNodes.insert(nodeName);


				}
			}

			bool isNewCluster = true;
			std::sort(cluster.begin(), cluster.end());

			//if(clusters.size() == 0){
				//isNewCluster = true;
			//}
			//else{
			for(const vector<u_int32_t>& existingCluster : clusters){
				float jaccardDistance = Utils::computeJaccardDistance(existingCluster, cluster);
				cout << jaccardDistance << endl;
				//float sharedRate = ((float) nbSharedNodes) / ((float)cluster.size());
				//cout << sharedRate << endl;
				if(jaccardDistance < 0.1){
					isNewCluster = false;
					break;
				}
			}
			//}
			if(!isNewCluster) continue;

			
			cout << clusterIndex << " " << component.size() << endl;
			
			if(validNodes.size() > 0){

				const string& filenameUnitigColor = clusterDir + "/component_" + to_string(clusterIndex) + "_unitigColor.csv";
				ofstream fUnitigColor(filenameUnitigColor);
				fUnitigColor << "Name,Colour" << endl;

				cout << "Generate contigs" << endl;
				const string& basespaceUnitigFilename = clusterDir + "/component_" + to_string(clusterIndex) + "_unitigs.fasta";
				//gzFile basespaceUnitigFile = gzopen(basespaceUnitigFilename.c_str(),"wb");
				ofstream basespaceUnitigFile = ofstream(basespaceUnitigFilename);

				u_int32_t contigIndex = 0;
				unordered_set<u_int32_t> writtenUnitigs;
				for(u_int32_t unitigIndex :component ){
					
					const Unitig& u = _graph->_unitigs[unitigIndex];
					
					//if(u._length < 10000) continue;
					
					if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._startNode)) != writtenUnitigs.end()) continue;
					if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._endNode)) != writtenUnitigs.end()) continue;

					writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._startNode));
					writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._endNode));

					string unitigSequence;
					_toBasespace.createSequence(u._nodes, unitigSequence);

					string header = ">ctg" + to_string(contigIndex) + '\n';
					basespaceUnitigFile << header << endl;
					//gzwrite(basespaceUnitigFile, (const char*)&header[0], header.size());
					unitigSequence +=  '\n';
					basespaceUnitigFile << unitigSequence << endl;
					//gzwrite(basespaceUnitigFile, (const char*)&unitigSequence[0], unitigSequence.size());
					
					for(u_int32_t nodeIndex : u._nodes){
						u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
						fUnitigColor << nodeName << "," << contigIndex << endl;
					}



					contigIndex += 1;
				}

				fUnitigColor.close();
				//gzclose(basespaceUnitigFile);
				basespaceUnitigFile.close();

				
				string command_annotation = "python3 /home/gats/workspace/tools/binner/thirdparty/scg/metacoag_main.py " + basespaceUnitigFilename + " " + filenameUnitigColor + " --isSingleSpecies";
				cout << command_annotation << endl;
				int ret = system(command_annotation.c_str());
				if(ret != 0){
					cerr << "Command failed: " << ret << endl;
					exit(ret);
				}

				char isSingleSpeciesResult;
				ifstream fileIsSingleSpecies(basespaceUnitigFilename + "_isSingleSpecies.txt");
				fileIsSingleSpecies.get(isSingleSpeciesResult);
				fileIsSingleSpecies.close();

				cout << isSingleSpeciesResult << endl;

				if(isSingleSpeciesResult == '1'){
									
					string outputFilename = clusterDir + "/component_" + to_string(clusterIndex) + ".gfa";
					GfaParser::rewriteGfa_withoutNodes(_gfaFilename, outputFilename, validNodes, _graph->_isEdgeRemoved, _graph->_graphSuccessors);
			
					string toFilename = singleSpeciesDir + "/component_" + to_string(clusterIndex) + ".gfa";
					fs::path pathFrom(outputFilename);
					fs::path pathTo(toFilename);
					fs::copy(pathFrom, pathTo);
				}



				/*
				//cout << "Start tree exploration" << endl;
				//unordered_set<u_int32_t> visitedNodes;
				//PathTree tree(_inputDir, unitig._startNode, unitig._abundance, _unitigDatas, graphSimplify);
				//tree.execute2(visitedNodes);
				//exit(1);

				const string& filenameComponentNodes = clusterDir + "/component_" + to_string(clusterIndex) + "_nodes.csv";
				ofstream fileComponenetNodes(filenameComponentNodes);
				fileComponenetNodes << "Name,Colour" << endl;
				for(u_int32_t nodeName : visitedNodes){
					fileComponenetNodes << nodeName << "," << clusterIndex << endl;
				}
				fileComponenetNodes.close();
				*/

				//exit(1);
				clusters.push_back(cluster);






			}

			clusterIndex += 1;
			cout << "cluster done" << endl;
			
		}

	}

	void indexReads_read(const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex){//}, const vector<KmerVec>& kminmers_k3, const vector<ReadKminmer>& kminmersInfos_k3){

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

		
		KminmerParser parser(_filename_readMinimizers, _minimizerSize, _kminmerSize);
		//auto fp = std::bind(&Assembly::indexReads_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
		auto fp = std::bind(&Assembly2::indexReads_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
		parser.parse(fp);
		
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
	ofstream file_groundTruth;
	ofstream file_groundTruth_hifiasmContigs;

	ofstream file_kminmersContigs;
	

	float computeAbundanceCutoff_min(u_int32_t abundance){
		return abundance / 4.0;
	}

	MinimizerParser* _minimizerParser;
	u_int32_t _extract_truth_kminmers_read_position;
	unordered_map<u_int32_t, vector<string>> _evaluation_hifiasmGroundTruth_nodeName_to_unitigName;
	vector<u_int32_t> _evaluation_hifiasmGroundTruth_path;
	unordered_set<u_int32_t> _hifiasm_startingNodenames;

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

		_extract_truth_kminmers_read_position = 0;
		_minimizerParser = new MinimizerParser(_minimizerSize, _minimizerDensity);
		
		auto fp = std::bind(&Assembly2::extract_truth_kminmers_read, this, std::placeholders::_1, std::placeholders::_2);
		ReadParser readParser(_truthInputFilename, true);
		readParser.parse(fp);

		delete _minimizerParser;
	}

	void computeScgCentrality(const string& outputFilenameBase, unordered_set<u_int32_t>& unitigs){
		//"shortest path weighted: est ce quon peut erreter l'algo des qu'on a trouver toutes les destination ? ou est-ce que les destination peuvent etre ameliorer, a check via des cout"

		u_int32_t nbPaths = 0;
		unordered_map<u_int32_t, u_int32_t> nbPathsTrought;
		
    	unordered_set<DbgEdge, hash_pair> sourceDests;

		for(auto& it1 : _scgClusters){
			u_int32_t unitigIndex1 = it1.first;
			const vector<ScgCluster>& scgClusters1 = it1.second;
			for(const ScgCluster& scgCluster1: scgClusters1){
				
				
				for(auto& it2 : _scgClusters){
					u_int32_t unitigIndex2 = it2.first;
					const vector<ScgCluster>& scgClusters2 = it2.second;
					for(const ScgCluster& scgCluster2: scgClusters2){
						
						
						if(scgCluster1._scgIndex == scgCluster2._scgIndex && scgCluster1._scgCluster == scgCluster2._scgCluster) continue;
						if(unitigIndex1 == unitigIndex2) continue;

						
						DbgEdge edge = {unitigIndex1, unitigIndex2};
						edge = edge.normalize();

						sourceDests.insert(edge);
					}

				}

				
			}

		}

		cout << sourceDests.size() << endl;

		unordered_set<u_int32_t> processedUnitigs;

		while(true){
		
			vector<DbgEdge> toProcess;
			u_int32_t sourceUnitigIndex = -1;

			for(const DbgEdge& sourceDest : sourceDests){

				u_int32_t unitigIndex1 = sourceDest._from;
				u_int32_t unitigIndex2 = sourceDest._to;


				if(sourceUnitigIndex == -1){

					if(processedUnitigs.find(unitigIndex1) == processedUnitigs.end()){
						sourceUnitigIndex = unitigIndex1;
					}
					else if(processedUnitigs.find(unitigIndex2) == processedUnitigs.end()){
						sourceUnitigIndex = unitigIndex2;
					}
				}

				if(sourceUnitigIndex == -1) continue;

				if(sourceUnitigIndex == unitigIndex1 || sourceUnitigIndex == unitigIndex2){
					toProcess.push_back(sourceDest);
				}

			}



			if(toProcess.size() == 0) break;
			processedUnitigs.insert(sourceUnitigIndex);

			//cout << "Source unitig index: " << sourceUnitigIndex << " " << toProcess.size() << endl;

			unordered_map<u_int32_t, vector<u_int32_t>> destPaths1;
			unordered_map<u_int32_t, vector<u_int32_t>> destPaths2;
			vector<u_int32_t> unitigIndexDests;

			for(const DbgEdge& sourceDest : toProcess){

				if(sourceDest._from != sourceUnitigIndex){
					destPaths1[sourceDest._from] = {};
					destPaths2[sourceDest._from] = {};
					unitigIndexDests.push_back(sourceDest._from);

					//cout << "Path: " << sourceUnitigIndex << " " << sourceDest._from << endl;
				}
				else if(sourceDest._to != sourceUnitigIndex){
					destPaths1[sourceDest._to] = {};
					destPaths2[sourceDest._to] = {};
					unitigIndexDests.push_back(sourceDest._to);
					//cout << "Path: " << sourceUnitigIndex << " " << sourceDest._to << endl;
				}
			}

			//cout << destPaths1.size() << endl;
			//cout << destPaths2.size() << endl;
			_graph->shortestPath_unitig_weighted(sourceUnitigIndex, destPaths1, false, false);
			_graph->shortestPath_unitig_weighted(_graph->unitigIndex_toReverseDirection(sourceUnitigIndex), destPaths2, false, false);

			for(u_int32_t unitigIndex : unitigIndexDests){

				//cout << destPaths1[unitigIndex].size() << " " << destPaths2[unitigIndex].size() << endl;
				/*
				vector<u_int32_t> path1;
				const vector<u_int32_t>& path1_forward = destPaths1[unitigIndex];
				const vector<u_int32_t>& path1_reverse = destPaths1[_graph->unitigIndex_toReverseDirection(unitigIndex)];
				if(path1_forward.size() == 0){
					path1 = path1_reverse;
				}
				else{
					path1 = path1_forward;
				}

				vector<u_int32_t> path2;
				const vector<u_int32_t>& path2_forward = destPaths2[unitigIndex];
				const vector<u_int32_t>& path2_reverse = destPaths2[_graph->unitigIndex_toReverseDirection(unitigIndex)];
				if(path2_forward.size() == 0){
					path2 = path2_reverse;
				}
				else{
					path2 = path2_forward;
				}
				*/
			
				u_int32_t path1_nbNodes = 0;
				u_int32_t path2_nbNodes = 0;
				for(u_int32_t unitigIndex : destPaths1[unitigIndex]){
					path1_nbNodes += _graph->_unitigs[unitigIndex]._nbNodes;
				}
				for(u_int32_t unitigIndex : destPaths2[unitigIndex]){
					path2_nbNodes += _graph->_unitigs[unitigIndex]._nbNodes;
				}


				if(path1_nbNodes < path2_nbNodes){
					
					for(u_int32_t unitigIndex : destPaths1[unitigIndex]){
						for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
							nbPathsTrought[BiGraph::nodeIndex_to_nodeName(nodeIndex)] += 1;
						}
					}
				}
				else{

					for(u_int32_t unitigIndex : destPaths2[unitigIndex]){
						for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
							nbPathsTrought[BiGraph::nodeIndex_to_nodeName(nodeIndex)] += 1;
						}
					}
				}

				nbPaths += 1;

			}
		}

		//cout << nbPathsTrought.size() << endl;

		/*
		vector<UnitigScgCluster> allScgClusters;

		for(auto& it : _scgClusters){
			
			u_int32_t unitigIndex = it.first;
			const vector<ScgCluster>& scgClusters = it.second;

			for(const ScgCluster& scgCluster: scgClusters){
				allScgClusters.push_back({unitigIndex, scgCluster._scgIndex, scgCluster._scgCluster});
			}

		}

		u_int32_t nbPaths = 0;
		unordered_map<u_int32_t, u_int32_t> nbPathsTrought;


		for(size_t i=0; i<allScgClusters.size(); i++){

			cout << i << " " << allScgClusters.size() << endl;

			const UnitigScgCluster& scgCluster1 = allScgClusters[i];
			unordered_map<u_int32_t, vector<u_int32_t>> destPaths1;
			unordered_map<u_int32_t, vector<u_int32_t>> destPaths2;

			vector<u_int32_t> unitigIndexDests;

			for(size_t j=0; j<allScgClusters.size(); j++){
			
				const UnitigScgCluster& scgCluster2 = allScgClusters[j];

				if(scgCluster1._scgIndex == scgCluster2._scgIndex && scgCluster1._scgCluster == scgCluster2._scgCluster) continue;
				if(scgCluster1._unitigIndex == scgCluster2._unitigIndex) continue;

				destPaths1[scgCluster2._unitigIndex] = {};
				destPaths1[_graph->unitigIndex_toReverseDirection(scgCluster2._unitigIndex)] = {};
				destPaths2[scgCluster2._unitigIndex] = {};
				destPaths2[_graph->unitigIndex_toReverseDirection(scgCluster2._unitigIndex)] = {};

				unitigIndexDests.push_back(scgCluster2._unitigIndex);
			}

			_graph->shortestPath_unitig_weighted(scgCluster1._unitigIndex, destPaths1, false, false);
			_graph->shortestPath_unitig_weighted(_graph->unitigIndex_toReverseDirection(scgCluster1._unitigIndex), destPaths2, false, false);

			for(u_int32_t unitigIndex : unitigIndexDests){

				vector<u_int32_t> path1;
				const vector<u_int32_t>& path1_forward = destPaths1[unitigIndex];
				const vector<u_int32_t>& path1_reverse = destPaths1[_graph->unitigIndex_toReverseDirection(unitigIndex)];
				if(path1_forward.size() == 0){
					path1 = path1_reverse;
				}
				else{
					path1 = path1_forward;
				}

				vector<u_int32_t> path2;
				const vector<u_int32_t>& path2_forward = destPaths2[unitigIndex];
				const vector<u_int32_t>& path2_reverse = destPaths2[_graph->unitigIndex_toReverseDirection(unitigIndex)];
				if(path2_forward.size() == 0){
					path2 = path2_reverse;
				}
				else{
					path2 = path2_forward;
				}

				u_int32_t path1_nbNodes = 0;
				u_int32_t path2_nbNodes = 0;
				for(u_int32_t unitigIndex : path1){
					path1_nbNodes += _graph->_unitigs[unitigIndex]._nbNodes;
				}
				for(u_int32_t unitigIndex : path2){
					path2_nbNodes += _graph->_unitigs[unitigIndex]._nbNodes;
				}


				if(path1_nbNodes < path2_nbNodes){
					
					for(u_int32_t unitigIndex : path1){
						for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
							nbPathsTrought[BiGraph::nodeIndex_to_nodeName(nodeIndex)] += 1;
						}
					}
				}
				else{

					for(u_int32_t unitigIndex : path2){
						for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
							nbPathsTrought[BiGraph::nodeIndex_to_nodeName(nodeIndex)] += 1;
						}
					}
				}

				nbPaths += 1;

			}

		}

		*/


		
		float min_value = 2;
		float max_value = -1;

		for(auto& it : nbPathsTrought){
			//u_int32_t componentIndex = unitigComponentIndex[i];
			//u_int32_t nbPaths = componentNbPaths[componentIndex];
			
			//float min_value = 2;
			//float max_value = -1;

			float centrality = it.second / ((float)nbPaths);
			if(centrality < min_value) min_value = centrality;
			if(centrality > max_value) max_value = centrality;
		}

		//cout << min_value << " " << max_value << endl;

		ofstream output_file(outputFilenameBase + ".centrality.csv");
		output_file << "Name,Colour" << endl;
		
		for(auto& it : nbPathsTrought){
			
			u_int32_t nodeName = it.first;
			u_int32_t nbPathNode = it.second;
			
			
			float centrality = nbPathNode / ((float)nbPaths);
			float centrality_normalized = (centrality-min_value) / (max_value-min_value);
			output_file << nodeName << "," << centrality_normalized << endl;
		}

		output_file.close();
		
	}
};



#endif
