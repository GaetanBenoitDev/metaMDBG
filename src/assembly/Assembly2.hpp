
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



#ifndef MDBG_METAG_ASSEMBLY2
#define MDBG_METAG_ASSEMBLY2

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






class PathExplorer{

public: 

	vector<u_int32_t> _prevNodes;
	GraphSimplify* _graph;
	unordered_set<u_int32_t> _visitedNodes;
	vector<UnitigData>& _unitigDatas;


	PathExplorer(const vector<u_int32_t>& prevNodes, unordered_set<u_int32_t>& visitedNodes, vector<UnitigData>& unitigDatas, GraphSimplify* graph) : _unitigDatas(unitigDatas){
		_prevNodes = prevNodes;
		_visitedNodes = visitedNodes;
		_graph = graph;
	}

	void getNextNode(u_int32_t current_nodeIndex, bool forward, vector<u_int32_t>& nextNodes){


		nextNodes.clear();
	
		if(isInInfiniteCycle(current_nodeIndex, forward)){
			return;
		}

		vector<u_int32_t> successors;
		if(forward){
			_graph->getSuccessors(current_nodeIndex, 0, successors);
		}
		else{
			_graph->getPredecessors(current_nodeIndex, 0, successors);
		}



		if(successors.size() == 0){
			
		}
		else if(successors.size() == 1){

			nextNodes.push_back(successors[0]);
		}
		else{

			computeBestSuccessors_byUnitigRank(successors, nextNodes);
		}

			


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

		possibleSuccessors.clear();


		u_int32_t prevRank = 0;
		size_t i = 0;

		u_int32_t lastNode = -1;
		u_int32_t pathLength = 0;


		while(true){
			


			int prevIndex = _prevNodes.size() - prevRank - 1;
			if(prevIndex < 0 ) break;

			u_int32_t prev_nodeIndex = _prevNodes[prevIndex];
			u_int32_t prev_nodeName = BiGraph::nodeIndex_to_nodeName(prev_nodeIndex);

			if(i == 0){
				pathLength = _graph->_nodeLengths[prev_nodeName];
			}
			else{
				u_int16_t overlapLength = _graph->_graphSuccessors->getOverlap(prev_nodeIndex, lastNode);
				pathLength += (_graph->_nodeLengths[prev_nodeName] - overlapLength);
			}

			
			#ifdef PRINT_PREV_RANK
				string str_debug = "";
				str_debug += "" + pathLengthgraph->_graphSuccessors->nodeToString(prev_nodeIndex);

				for(u_int32_t& successor : successors){
					u_int32_t successor_nodeName = BiGraph::nodeIndex_to_nodeName(successor);
					u_int32_t nbSharedReads = Utils::computeSharedReads(_unitigDatas[prev_nodeName], _unitigDatas[successor_nodeName]);
					str_debug += "    " + graph->_graphSuccessors->nodeToString(successor._nodeIndex) + " " + to_string(nbSharedReads);
				}
				str_debug += "    " + to_string(pathLength);
				cout << "\t" << str_debug << endl;
			#endif


			if(pathLength > 9000){
				for(u_int32_t& successor : successors){
					u_int32_t successor_nodeName = BiGraph::nodeIndex_to_nodeName(successor);
					u_int32_t nbSharedReads = Utils::computeSharedReads(_unitigDatas[prev_nodeName], _unitigDatas[successor_nodeName]);
					if(nbSharedReads >= 2){
						possibleSuccessors.push_back(successor);
					}

				}
				break;
			}


			prevRank += 1;


		}


	}



	void nodeExplored(u_int32_t nodeIndex){
		
		_prevNodes.push_back(nodeIndex);
		//_exploredNodes.push_back(nodeIndex);

		PathExplorer::clampPrevNodes(_prevNodes, _unitigDatas);
		u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
		//_currentPathLength += graph->_nodeLengths[nodeName];

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

	PathTree(u_int32_t source_nodeIndex, u_int32_t source_abundance, vector<UnitigData>& unitigDatas, GraphSimplify* graph) : _unitigDatas(unitigDatas){
		_source_abundance = source_abundance;
		_source_nodeIndex = source_nodeIndex;
		_graph = graph;
	}


	void execute(){

		//u_int32_t current_nodeIndex = _source_nodeIndex;

		bool forward = true;
		vector<PathExplorer> queue;

		vector<u_int32_t> prevNodes;
		unordered_set<u_int32_t> visitedNodes;

		PathExplorer pathExplorer_source (prevNodes, visitedNodes, _unitigDatas, _graph);
		pathExplorer_source.nodeExplored(_source_nodeIndex);

		queue.push_back(pathExplorer_source);
		cout << queue.size() << endl;

		while(!queue.empty()){

			cout << queue.size() << endl;

			PathExplorer pathExplorer = queue[queue.size()-1];
			queue.pop_back();

			u_int32_t current_nodeIndex = pathExplorer._prevNodes[pathExplorer._prevNodes.size()-1];

			while(true){

				vector<u_int32_t> nextNodes;
				pathExplorer.getNextNode(current_nodeIndex, forward, nextNodes);

				if(current_nodeIndex == _source_nodeIndex) break;

				if(nextNodes.size() == 0){
					break;
				}
				else if(nextNodes.size() == 1){
					pathExplorer.nodeExplored(nextNodes[0]);
				}
				else{
					for(u_int32_t nodeIndex : nextNodes){
						pathExplorer.nodeExplored(nodeIndex);
						queue.push_back(pathExplorer);
					}
					break;
				}

			}
		}

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


	Assembly2(): Tool (){

	}

	~Assembly2(){

	}

	void execute (){

		_outputContigFile = gzopen(_outputFilename.c_str(),"wb");
		_outputContigFile_complete = gzopen(_outputFilename_complete.c_str(),"wb");

		execute_assembly();
		
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


	void execute_assembly(){

		//size_t globalAbundanceCutoff_min = 3;
		
		string gfa_filename = _inputDir + "/minimizer_graph.gfa";
		string gfa_filename_noUnsupportedEdges = _inputDir + "/minimizer_graph_noUnsupportedEdges.gfa";
		//string gfa_filename_unitigs = _inputDir + "/minimizer_graph_unitigs.gfa";
		string mdbg_filename = _inputDir + "/mdbg_nodes.gz";




		
		cout << gfa_filename << endl;
		_mdbg = new MDBG(_kminmerSize);
		_mdbg->load(mdbg_filename);
		cout << "Nb nodes: " <<  _mdbg->_dbg_nodes.size() << endl;





		file_groundTruth = ofstream(_inputDir + "/binning_results.csv");
		file_groundTruth << "Name,Colour" << endl;

		file_groundTruth_hifiasmContigs = ofstream(_inputDir + "/binning_results_hifiasm_" + to_string(_kminmerSize) + ".csv");
		file_groundTruth_hifiasmContigs << "Name,Colour" << endl;

		if(_debug){
            gfa_filename = _inputDir + "/minimizer_graph_debug.gfa";
		}
		
		GraphSimplify* graphSimplify = new GraphSimplify(gfa_filename, _inputDir, 0);
		
		
		
		//Generate unitigs
		cout << "Indexing reads" << endl;
		_unitigDatas.resize(_mdbg->_dbg_nodes.size());
		graphSimplify->clear(0);
		graphSimplify->compact(false, _unitigDatas);
		removeUnsupportedEdges(gfa_filename, gfa_filename_noUnsupportedEdges, graphSimplify);
		delete _mdbg;
		cout << "done" << endl;
	
		//cout << graphSimplify->_graphSuccessors->_nbEdges << endl;
		//graphSimplify->execute(5, _kminmerSize);
		//graphSimplify->debug_writeGfaErrorfree(1000, PathExplorer::computeAbundanceCutoff(1000, 0, CutoffType::ERROR), -1, _kminmerSize, false, true, false, _unitigDatas);
		graphSimplify->debug_writeGfaErrorfree(100, 100, -1, _kminmerSize, false, true, false, _unitigDatas);















		u_int64_t nbUnitigs = 0;
		unordered_set<u_int32_t> writtenNodeNames;
		unordered_map<u_int32_t, u_int32_t> nodeName_to_unitigName;


		ToBasespaceOnTheFly toBasespace(_inputDir);



		string clusterDir = _inputDir + "/" + "components";
		fs::path path(clusterDir);
	    if(!fs::exists (path)){
            fs::create_directory(path);
        } 

		u_int32_t clusterIndex = 0;
		unordered_set<string> written_unitigName;
		float prevCutoff = -1;

		vector<vector<u_int32_t>> clusters;

		for(Unitig& unitig : graphSimplify->_startingUnitigstest){

			u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(unitig._startNode);
			//if(nodeName_to_cluster.find(nodeName) != nodeName_to_cluster.end()) continue;

			float cutoff = unitig._abundance*0.2;

			//float minUnitigAbundance = unitig._abundance*0.8;

			float realCutoff = 0;
			for(const SaveState2& saveState : graphSimplify->_cachedGraphStates){
				//cout << saveState._abundanceCutoff_min << endl;
				if(saveState._abundanceCutoff_min > cutoff) break;
				realCutoff = saveState._abundanceCutoff_min;
			}


			if(realCutoff != prevCutoff){
				graphSimplify->loadState2(cutoff, unitig._startNode, _unitigDatas);
				prevCutoff = realCutoff;
			}

        	unordered_set<u_int32_t> component;
        	graphSimplify->getConnectedComponent(unitig._startNode, component);

			if(component.size() <= 1) continue;

			vector<u_int32_t> cluster;
			unordered_set<u_int32_t> validNodes;
			//unordered_set<u_int32_t> clusterUnitigNames;

			for(u_int32_t unitigIndex : component){
				
				//if(graphSimplify->_unitigs[unitigIndex]._length < 10000) continue;
				//if(graphSimplify->_unitigs[unitigIndex]._abundance < minUnitigAbundance) continue;

				for(u_int32_t nodeIndex : graphSimplify->_unitigs[unitigIndex]._nodes){
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
				u_int64_t nbSharedNodes = Utils::computeSharedElements(existingCluster, cluster);
				float sharedRate = ((float) nbSharedNodes) / ((float)cluster.size());
				cout << sharedRate << endl;
				if(sharedRate > 0.9){
					isNewCluster = false;
					break;
				}
			}
			//}
			if(!isNewCluster) continue;

			cout << clusterIndex << " " << component.size() << endl;
			
			if(component.size() > 50 && validNodes.size() > 0){

				cout << "Start tree exploration" << endl;
				PathTree tree(unitig._startNode, 0, _unitigDatas, graphSimplify);
				tree.execute();
				exit(1);

				clusters.push_back(cluster);
				
				string outputFilename = clusterDir + "/component_" + to_string(clusterIndex) + ".gfa";
				GfaParser::rewriteGfa_withoutNodes(gfa_filename, outputFilename, validNodes, graphSimplify->_isEdgeRemoved, graphSimplify->_graphSuccessors);
			



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
					
					const Unitig& u = graphSimplify->_unitigs[unitigIndex];
					
					//if(u._length < 10000) continue;
					
					if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._startNode)) != writtenUnitigs.end()) continue;
					if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._endNode)) != writtenUnitigs.end()) continue;

					writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._startNode));
					writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._endNode));

					string unitigSequence;
					toBasespace.createSequence(u._nodes, unitigSequence);

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

				
				string command_annotation = "python3 /home/gats/workspace/tools/binner/thirdparty/scg/metacoag_main.py " + basespaceUnitigFilename + " " + filenameUnitigColor;
				cout << command_annotation << endl;
				int ret = system(command_annotation.c_str());
				if(ret != 0){
					cerr << "Command failed: " << ret << endl;
					exit(ret);
				}

			}

			clusterIndex += 1;
		}





		//getchar();

		file_groundTruth.close();
		file_groundTruth_hifiasmContigs.close();
		file_kminmersContigs.close();



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

};



#endif
