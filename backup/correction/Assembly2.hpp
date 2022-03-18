
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

/*
Solve repeat: attention startNode et endNode d'un meme unitig pourrait etre un meme successor/predecessor si l'unitig est couvert par plusieurs long reads (il faut choisir le bon dans ce cas)
ToBasespaceOnTheFly: en théorie, on peut reconstruire la partie manquante des contigs basespace (les k-1 premiers nodes), en utilisant les predecesseurs de ce noeud de départ

- singleton erroneous unitigs dans le graphe (si une erreur est isolé, ça créé qq noeuds attaché erroneous, on pourrait ptet supprimé les mini component peut couverte à voir)
*/

//:/mnt/chris-native/chris/Projects/AD_HiFi
//mnt/gpfs/Hackathon/meren_hifi

#ifndef MDBG_METAG_ASSEMBLY3
#define MDBG_METAG_ASSEMBLY3

#define PRINT_DEBUG_PREVRANK

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










class Assembly3 : public Tool{

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
	MDBG* _mdbgNoFilter;
	GraphSimplify* _graph;
	//ToBasespaceOnTheFly _toBasespace;
	//ContigFeature _contigFeature;
	string _gfaFilename;
	u_int64_t _nbReads;

	Assembly3(): Tool (){

	}

	~Assembly3(){

	}

	void execute (){

		loadGraph();
		//generateContigs2(_inputDir + "/contigs.nodepath.gz", _inputDir + "/contigs.fasta.gz");

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

		ifstream file_data(_inputDir + "/data.txt");
		file_data.read((char*)&_nbReads, sizeof(_nbReads));
		file_data.close();

		cout << endl;
		cout << "Input dir: " << _inputDir << endl;
		//cout << "Output filename: " << _outputFilename << endl;
		cout << "Minimizer length: " << _minimizerSize << endl;
		cout << "Kminmer length: " << _kminmerSize << endl;
		cout << "Density: " << _minimizerDensity << endl;
		cout << "Nb reads: " << _nbReads << endl;
		cout << endl;

		//getchar();

		_outputFilename = _inputDir + "/minimizer_contigs.gz";
		_outputFilename_complete = _inputDir + "/minimizer_contigs_complete.gz";
		_filename_readMinimizers = _inputDir + "/read_data.txt";
		_filename_hifiasmGroundtruth = _inputDir + "/hifiasmGroundtruth.gz";
		_filename_outputContigs = _inputDir + "/contigs.min.gz";
		_filename_solidKminmers = _inputDir + "/solid.min.gz";


		

	}

	unordered_map<u_int64_t, vector<u_int32_t>> _debug_readPath;
    vector<bool> _isBubble;
	unordered_map<u_int32_t, vector<u_int64_t>> _nodeName_to_kminmerSequence;
	

	void loadGraph(){

		_gfaFilename = _inputDir + "/minimizer_graph.gfa";
		string gfa_filename_noUnsupportedEdges = _inputDir + "/minimizer_graph_noUnsupportedEdges.gfa";
		//string gfa_filename_unitigs = _inputDir + "/minimizer_graph_unitigs.gfa";
		string mdbg_filename = _inputDir + "/mdbg_nodes.gz";




		
		_mdbgNoFilter = new MDBG(_kminmerSize);
		_mdbgNoFilter->load(_inputDir + "/mdbg_nodes_noFilter.gz");

		cout << _gfaFilename << endl;
		_mdbg = new MDBG(_kminmerSize);
		_mdbg->load(mdbg_filename);
		cout << "Nb nodes: " <<  _mdbg->_dbg_nodes.size() << endl;

		//extractContigKminmers2(_inputDir + "/contigs.fasta.gz");

		if(_truthInputFilename != ""){
			extract_truth_kminmers();
		}


		file_groundTruth = ofstream(_inputDir + "/binning_results.csv");
		file_groundTruth << "Name,Colour" << endl;

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
		//generateContigs();
		removeUnsupportedEdges(_gfaFilename, gfa_filename_noUnsupportedEdges, _graph);
		

		cout << "done" << endl;
	

		//if(_kminmerSize == 4){
			//cout << graphSimplify->_graphSuccessors->_nbEdges << endl;
			//graphSimplify->execute(5, _kminmerSize);
			//graphSimplify->debug_writeGfaErrorfree(1000, PathExplorer::computeAbundanceCutoff(1000, 0, CutoffType::ERROR), -1, _kminmerSize, false, true, false, _unitigDatas);


			cout << "Cleanning graph 1" << endl;
			_graph->debug_writeGfaErrorfree(0, 0, -1, _kminmerSize, false, true, false, _unitigDatas, true, false, true, false, true);
			_isBubble = _graph->_isBubble;
			
			//cout << "Cleanning graph 2" << endl;
			//_graph->debug_writeGfaErrorfree(0, 0, -1, _kminmerSize, false, true, false, _unitigDatas, false, false, false, false);

			//!

			cout << "done" << endl;
			_graph->loadState2(0, -1, _unitigDatas);
		//}

		
		
		cout << _graph->_isNodenameRoundabout.size() << endl;
		
		
		ofstream file_correction(_inputDir + "/roundabout.csv");
		file_correction << "Name,Colour" << endl;
		for(u_int32_t nodeName : _graph->_isNodenameRoundabout){
			if(!_isBubble[BiGraph::nodeName_to_nodeIndex(nodeName, true)]){
				file_correction << nodeName << ",red" << endl;
			}
		}
		for(u_int32_t nodeIndex : _graph->_isNodeValid2){
			if(_graph->_isBubble[nodeIndex]){
				//file_correction << BiGraph::nodeIndex_to_nodeName(nodeIndex) << ",red" << endl;
			}
		}
		file_correction.close();
		
		//cout << _graph->nodeIndex_to_unitig(BiGraph::nodeName_to_nodeIndex(12871, true))._nbNodes << " " << _graph->nodeIndex_to_unitig(BiGraph::nodeName_to_nodeIndex(12871, true))._length << endl;
		//getchar();
		
		
		//_graph->loadState2(100, -1, _unitigDatas);
		//_graph->saveGraph(_inputDir + "/minimizer_graph_contigs.gfa");
		_graph->loadState2(0, -1, _unitigDatas);

		_partitionDir = _inputDir + "/" + "partitions";
		fs::path path(_partitionDir);
		if(fs::exists (path)){
			fs::remove_all(path);
		} 
		fs::create_directory(path);

		//_graph->loadState2(0, -1, _unitigDatas);
		_graph->saveGraph(_inputDir + "/minimizer_graph_sub.gfa");
		//_graph->loadState2(7, -1, _unitigDatas);
		//_graph->saveGraph(_inputDir + "/minimizer_graph_sub_2.gfa");
		/*
		if(_kminmerSize == 4){
			//cout << "A ENLEVER" << endl;
			_graph->loadState2(10, -1, _unitigDatas);
			_graph->saveGraph(_inputDir + "/minimizer_graph_sub_2.gfa");
		}
		else{
			//cout << "A ENLEVER" << endl;
			_graph->loadState2(2.5, -1, _unitigDatas);
			_graph->saveGraph(_inputDir + "/minimizer_graph_sub_2.gfa");
		}
		*/

		//cout << Utils::computeSharedReads(_unitigDatas[8228], _unitigDatas[8229]) << endl;
		//cout << Utils::computeSharedReads(_unitigDatas[6323], _unitigDatas[2778]) << endl;
		//cout << Utils::computeSharedReads(_unitigDatas[8228], _unitigDatas[2778]) << endl;
		//cout << Utils::computeSharedReads(_unitigDatas[6323], _unitigDatas[8229]) << endl;
		
		//cout << Utils::computeSharedReads(_unitigDatas[3076], _unitigDatas[4346]) << endl;
		//cout << Utils::computeSharedReads(_unitigDatas[4506], _unitigDatas[4886]) << endl;
		//cout << Utils::computeSharedReads(_unitigDatas[3076], _unitigDatas[4886]) << endl;
		//cout << Utils::computeSharedReads(_unitigDatas[3031], _unitigDatas[8043]) << endl;
		//cout << Utils::computeSharedReads(_unitigDatas[62], _unitigDatas[1106]) << endl;
		//cout << Utils::computeSharedReads(_unitigDatas[9895], _unitigDatas[933]) << endl;
		//cout << Utils::computeSharedReads(_unitigDatas[12321], _unitigDatas[933]) << endl;
		
		//for(auto& it : _evaluation_hifiasmGroundTruth_nodeNamePosition){
		//	if(it.second == 4790){
		//		cout << it.first << " " << _unitigDatas[it.first]._readIndexes.size() << endl;
		//	}
		//}


		//getchar();
		//debug_checkReads();
		//cout << endl << endl;
		//partitionReads();

		extractKminmerSequences();
		correctReads();
		//getchar();

		delete _mdbg;

		//_toBasespace.create(_inputDir);


		//_graph->loadState2(0, -1, _unitigDatas);
		//generateFasta(_inputDir + "/allContigs.fasta.gz");
		//generateFastaExpanded(_inputDir + "/allContigsExpanded.fasta.gz");

	}

	void extractKminmerSequences (){

		cout << "Extracting kminmer sequences" << endl;
		
		
		ifstream kminmerFile(_inputDir + "/kminmerData_min.txt");

		while (true) {

			u_int16_t size;
			kminmerFile.read((char*)&size, sizeof(size));

			if(kminmerFile.eof())break;

			vector<u_int64_t> minimizerSeq;
			minimizerSeq.resize(size);
			kminmerFile.read((char*)&minimizerSeq[0], size*sizeof(u_int64_t));

			u_int32_t nodeName;
			u_int32_t length;
			u_int32_t lengthStart;
			u_int32_t lengthEnd;
			//bool isReversed = false;

			kminmerFile.read((char*)&nodeName, sizeof(nodeName));
			kminmerFile.read((char*)&length, sizeof(length));
			kminmerFile.read((char*)&lengthStart, sizeof(lengthStart));
			kminmerFile.read((char*)&lengthEnd, sizeof(lengthEnd));

			_nodeName_to_kminmerSequence[nodeName] = minimizerSeq;

		}

		kminmerFile.close();
	}

	void indexReads_read(const vector<u_int64_t>& minimizers, const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex){//}, const vector<KmerVec>& kminmers_k3, const vector<ReadKminmer>& kminmersInfos_k3){

		if(_indexingContigs) readIndex += 2000000000ull;
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
	
	bool _indexingContigs;

	void removeUnsupportedEdges(const string& gfaFilename, const string& gfa_filename_noUnsupportedEdges, GraphSimplify* graph){

		_indexingContigs = false;
		
		KminmerParser parser(_filename_readMinimizers, _minimizerSize, _kminmerSize, true);
		//auto fp = std::bind(&Assembly::indexReads_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
		auto fp = std::bind(&Assembly3::indexReads_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
		parser.parse(fp);
		
		//_indexingContigs = true;
		//KminmerParser parser2(_inputDir + "/contig_data.txt", _minimizerSize, _kminmerSize, false);
		//auto fp = std::bind(&Assembly::indexReads_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
		//auto fp2 = std::bind(&Assembly3::indexReads_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
		//parser2.parse(fp2);
		
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
	ofstream file_groundTruth;
	//ofstream file_groundTruth_hifiasmContigs;

	ofstream file_kminmersContigs;
	


	MinimizerParser* _minimizerParser;
	u_int32_t _extract_truth_kminmers_read_position;
	unordered_map<u_int32_t, vector<string>> _evaluation_hifiasmGroundTruth_nodeName_to_unitigName;
	vector<u_int32_t> _evaluation_hifiasmGroundTruth_path;
	unordered_set<u_int32_t> _hifiasm_startingNodenames;
	ofstream _file_groundTruth_hifiasm_position;
	EncoderRLE _encoderRLE;

	void extract_truth_kminmers_read(const Read& read){

		u_int64_t readIndex = read._index;
		//ottalSize += strlen(read->seq.s);


		string rleSequence;
		vector<u_int64_t> rlePositions;
		_encoderRLE.execute(read._seq.c_str(), read._seq.size(), rleSequence, rlePositions);

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
				//if("utg009434l" == string(read->name.s)){
				//	cout << nodeName << endl;
				//	_hifiasm_startingNodenames.insert(nodeName);
				//}

				_evaluation_hifiasmGroundTruth_nodeName_to_unitigName[_mdbg->_dbg_nodes[vec]._index].push_back(read._header);
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
		
		auto fp = std::bind(&Assembly3::extract_truth_kminmers_read, this, std::placeholders::_1);
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

    struct NodeNameAbundance{
        u_int32_t _nodeName;
        u_int32_t _abundance;
    };
    
    struct NodeNameAbundance_Comparator {
        bool operator()(NodeNameAbundance const& p1, NodeNameAbundance const& p2){
            return p1._abundance > p2._abundance;
        }
    };

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


					cout << endl << "Starting unitig: " << BiGraph::nodeIndex_to_nodeName(unitig._startNode) << " " << unitig._length << endl;


					for(u_int32_t nodeIndex : unitig._nodes){
						processedNodeNames.insert(BiGraph::nodeIndex_to_nodeName(nodeIndex));
					}


					vector<u_int32_t> nodePath = unitig._nodes;
					vector<u_int64_t> nodePath_supportingReads;
					//assembly.solveBin2(unitig._startNode, unitig._abundance, _graph, 0, 0, false, nodePath, nodePath_supportingReads, 0);




					string unitigSequenceExpanded;
					//_toBasespace.createSequence(nodePath, unitigSequenceExpanded);
					cout << "\tExpanded contigs: " << nodePath.size() << endl; //<< " " << unitigSequenceExpanded.size() << endl;

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
					cout << "\tNb valid contigs: " << validContigs << endl;
					
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

					cout << "----" << endl;
					
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


	unordered_map<u_int32_t, ofstream> _readPartitions;

	
	void partitionReads(){

		//cout << "Partitionning reads" << endl;

		u_int32_t cutoffLevel = 0;
		for(const SaveState2& saveState : _graph->_cachedGraphStates){
			const string filename = _partitionDir + "/part_" + to_string(cutoffLevel) + ".gz";
			_readPartitions[cutoffLevel] = ofstream(filename);
			cutoffLevel += 1;
		}

		_graph->loadState2(0, -1, _unitigDatas);

		
		KminmerParser parser(_filename_readMinimizers, _minimizerSize, _kminmerSize, true);
		auto fp = std::bind(&Assembly3::partitionReads_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
		parser.parse(fp);

		for(auto& it : _readPartitions){
			it.second.close();
		}
	}

	void partitionReads_read(const vector<u_int64_t>& minimizers, const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex){//}, const vector<KmerVec>& kminmers_k3, const vector<ReadKminmer>& kminmersInfos_k3){
	
		vector<u_int32_t> nodePath;

		for(const KmerVec& vec : kminmers){
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()){
				//cout << "XXXXX ";
				continue;
			}
			
			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
			nodePath.push_back(nodeName);
			//cout << nodeName << " ";
		}
		//cout << endl;

		double n = 0;
		double sum = 0;
		unordered_set<u_int32_t> writtenUnitigs;

		vector<float> readpathAbudance_values;
		for(u_int32_t nodeName : nodePath){

			u_int32_t nodeIndex = BiGraph::nodeName_to_nodeIndex(nodeName, false);
			u_int32_t unitigIndex = _graph->nodeIndex_to_unitigIndex(nodeIndex);
			const Unitig& u = _graph->_unitigs[unitigIndex];

			if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._startNode)) != writtenUnitigs.end()) continue;
			if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._endNode)) != writtenUnitigs.end()) continue;

			writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._startNode));
			writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._endNode));

			for(u_int32_t nodeIndex : u._nodes){
				if(_graph->_isNodenameRoundabout.find(BiGraph::nodeIndex_to_nodeName(nodeIndex)) != _graph->_isNodenameRoundabout.end()) continue;
				if(_graph->_isBubble[nodeIndex]) continue;

				readpathAbudance_values.push_back(u._abundance);
				//cout << u._abundance << " ";
				n += 1;
				sum += u._abundance;
			}
		}
		//cout << endl;
		
		float cutoff = 0;

		if(n > 0){
			float readPathAbundance = sum / n; ////Utils::compute_median_float(readpathAbudance_values);
			//cout << "Read path abundance: " << readPathAbundance << " " << (sum / n) << endl;
			//if(readPathAbundance < 10) getchar();

			//if(readPathAbundance > 50) getchar();
			cutoff = readPathAbundance * 0.2;
			//cout << cutoff << endl;
		}

		//cutoff = 4;
		//if(_kminmerSize != 4) 
		//cutoff = 0;

		u_int32_t cutoffLevel = 0;
		float realCutoff = 0;
		for(const SaveState2& saveState : _graph->_cachedGraphStates){
			if(saveState._abundanceCutoff_min >= cutoff) break;
			realCutoff = saveState._abundanceCutoff_min;
			cutoffLevel += 1;
		}

		ofstream& file = _readPartitions[cutoffLevel];
		
		u_int32_t size = minimizers.size();
		file.write((const char*)&size, sizeof(size));
		file.write((const char*)&minimizers[0], size * sizeof(u_int64_t));
	}


	//ofstream file_correction;
	//gzFile _outputFile_correctedReads;
	ofstream _file_uncorrectedReads;

	void correctReads(){

		_nbreadsLala = 0;
		_nbUncorrectedReads = 0;

		_graph->loadState2(0, -1, _unitigDatas);
		_graph->_isBubble = _isBubble;

		//const string& outputFilename_correctedReads = _inputDir + "/correctedReads_" + to_string(_kminmerSize) + ".min.gz";
		//_outputFile_correctedReads = gzopen(outputFilename_correctedReads.c_str(),"wb");

		_file_uncorrectedReads = ofstream(_inputDir + "/read_uncorrected.txt");
		
		/*
		u_int32_t cutoffLevel = 0;
		for(const SaveState2& saveState : _graph->_cachedGraphStates){

			cout << cutoffLevel << " " << saveState._abundanceCutoff_min << endl;
			//getchar();

			_graph->loadState2(saveState._abundanceCutoff_min, -1, _unitigDatas);
			//_graph->saveGraph(_inputDir + "minimizer_graph_sub.gfa");
			const string readPartitionFilename = _partitionDir + "/part_" + to_string(cutoffLevel) + ".gz";
			
			KminmerParser parser(readPartitionFilename, _minimizerSize, _kminmerSize, false);
			auto fp = std::bind(&Assembly3::correctReads_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
			parser.parse(fp);





			//_readPartitions[cutoffLevel] = gzopen(filename.c_str(),"wb");
			cutoffLevel += 1;
		}
		*/



		
		KminmerParser parser(_filename_readMinimizers, _minimizerSize, _kminmerSize, true);
		auto fp = std::bind(&Assembly3::correctReads_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
		parser.parse(fp);
		

		//file_correction.close();
		//gzclose(_outputFile_correctedReads);
		_file_uncorrectedReads.close();

		cout << "Nb reads: " << _nbreadsLala << endl;
		cout << "Nb uncorrected reads: " << _nbUncorrectedReads << endl;
		//getchar();
		
	}

	u_int64_t _nbreadsLala;
	u_int64_t _nbUncorrectedReads = 0;
	unordered_map<u_int32_t, vector<u_int64_t>> _nodeIndex_to_kminmerSequence;

	void correctReads_read(const vector<u_int64_t>& minimizers, const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex){//}, const vector<KmerVec>& kminmers_k3, const vector<ReadKminmer>& kminmersInfos_k3){

		//u_int32_t readSize = minimizers.size();
		//_file_uncorrectedReads.write((const char*)&readSize, sizeof(readSize));
		//_file_uncorrectedReads.write((const char*)&minimizers[0], readSize*sizeof(u_int64_t));
		//return;

		/*
		if(readIndex >= _nbReads){
			u_int32_t readSize = minimizers.size();
			_file_uncorrectedReads.write((const char*)&readSize, sizeof(readSize));
			_file_uncorrectedReads.write((const char*)&minimizers[0], readSize*sizeof(u_int64_t));
			return;
		}
		*/

		bool print_read = false;
		if(readIndex % 10000 == 0) print_read = true;

		//file_correction = ofstream(_inputDir + "/correction.csv");
		//file_correction << "Name,Color" << endl;

		if(print_read) cout << "-------------" << endl;
		if(print_read) cout << readIndex << endl;
		//if(print_read) cout << minimizers.size() << endl;


		/*
		bool isHere = false;
		vector<u_int32_t> nodePath;
		vector<u_int32_t> nodePath_withMissing;
	
		if(print_read) cout << "\t";
		for(const KmerVec& vec : kminmers){
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()){
				nodePath_withMissing.push_back(-1);
				if(print_read) cout << "XXXXX ";
				continue;
			}
			
			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;

			if(_graph->_isBubble[BiGraph::nodeName_to_nodeIndex(nodeName, true)]){
				if(print_read) cout << nodeName << "B ";
				//continue;
			}
			else{
				if(print_read) cout << nodeName << " ";
			}

			
			if(nodeName == 2630) isHere = true;
			//if(readIndex == 1510){
			//file_correction << nodeName << "," << "red" << endl;
			//}
			nodePath.push_back(nodeName);
			nodePath_withMissing.push_back(nodeName);
		}
		if(print_read) cout << endl;

		if(print_read){
			cout << "\tRead minimizers: " << endl;
			cout << "\t";
			for(u_int64_t minimizer : minimizers){
				cout << minimizer << " ";
			}
			cout << endl;
		}

		//if(_kminmerSize == 7 && isHere) getchar();
		*/

		if(print_read){
			cout << "\tRead original: " << endl;
			cout << "\t";
			for(u_int64_t m : minimizers){
				cout << m << " ";
			}
			cout << endl;

			cout << "\tRead nodes: " << endl;
			cout << "\t";
		}


		_nodeIndex_to_kminmerSequence.clear();
		vector<u_int32_t> nodepath;

		for(size_t i=0; i<kminmers.size(); i++){

			const KmerVec& vec = kminmers[i];
			const ReadKminmer& info = kminmersInfos[i];

			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()){
				if(print_read) cout << "XXXXX ";
				continue;
			}
			
			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
			u_int32_t nodeIndex = BiGraph::nodeName_to_nodeIndex(nodeName, !info._isReversed);

			vector<u_int64_t> kminmerSeq = _nodeName_to_kminmerSequence[nodeName];

			if(_graph->_isBubble[BiGraph::nodeName_to_nodeIndex(nodeName, true)]){
				if(print_read) cout << nodeName << "B ";
				//continue;
			}
			else{
				if(print_read) cout << nodeName << " ";
			}

			if(info._isReversed){
				std::reverse(kminmerSeq.begin(), kminmerSeq.end());
			}

			nodepath.push_back(nodeIndex);
			//_nodeIndex_to_kminmerSequence[nodeIndex] = kminmerSeq;

			//cout << BiGraph::nodeToString(nodeIndex) << ": ";
			//for(u_int64_t m : kminmerSeq){
			//	cout << m << " ";
			//}
			//cout << endl;
		}
		if(print_read) cout << endl;


		/*
		double n = 0;
		double sum = 0;
		unordered_set<u_int32_t> writtenUnitigs;

		vector<float> readpathAbudance_values;
		for(u_int32_t nodeName : nodePath){

			u_int32_t nodeIndex = BiGraph::nodeName_to_nodeIndex(nodeName, false);
			u_int32_t unitigIndex = _graph->nodeIndex_to_unitigIndex(nodeIndex);
			const Unitig& u = _graph->_unitigs[unitigIndex];

			if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._startNode)) != writtenUnitigs.end()) continue;
			if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._endNode)) != writtenUnitigs.end()) continue;

			writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._startNode));
			writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._endNode));

			for(u_int32_t nodeIndex : u._nodes){
				readpathAbudance_values.push_back(u._abundance);
				//cout << u._abundance << " ";
				n += 1;
				sum += u._abundance;
			}
		}
		//cout << endl;
		
		float readPathAbundance = Utils::compute_median_float(readpathAbudance_values); //sum / n;
		cout << "Read path abundance: " << readPathAbundance << " " << (sum / n) << endl;
		*/


		vector<u_int32_t> nodePathSolid;
		vector<u_int64_t> readpath;
		applyReadCorrection(nodepath, readpath, print_read, minimizers, kminmers, kminmersInfos, nodePathSolid);


		if(print_read){
			cout << "\tRead original: " << endl;
			cout << "\t";
			for(u_int64_t m : minimizers){
				cout << m << " ";
			}
			cout << endl;
		}

		if(print_read){
			cout << "\tRead minimizers corrected: " << endl;
			cout << "\t";
			for(u_int64_t minimizer : readpath){
				cout << minimizer << " ";
			}
			cout << endl;
		}

		if(nodePathSolid.size() != 2){
			if(print_read) cout << "\tcorrection failed" << endl;
			readpath = minimizers;
			_nbUncorrectedReads += 1;
		}
		else{
			extendReadpath(minimizers, kminmers, nodePathSolid, readpath, print_read);
		}
		
		//readpath = minimizers;
		//vector<u_int32_t> contigpath;

		if(print_read){
			cout << "\tRead minimizers corrected extended: " << endl;
			cout << "\t";
			for(u_int64_t minimizer : readpath){
				cout << minimizer << " ";
			}
			cout << endl;
		}


		_nbreadsLala += 1;
		/*

		bool isPathErroneous = false;
		//cout << "lala" << endl;
		//cout << readpath.size() << endl;
		for(u_int32_t nodeIndex : readpath){
			//cout << BiGraph::nodeIndex_to_nodeName(nodeIndex) << " ";

			if(!_truthInputFilename.empty()){
				if(_evaluation_hifiasmGroundTruth_nodeNamePosition.find(BiGraph::nodeIndex_to_nodeName(nodeIndex)) == _evaluation_hifiasmGroundTruth_nodeNamePosition.end()){
					isPathErroneous = true;
				}
			}
		}
		*/
		//cout << endl;

		/*
		if(isPathErroneous){
			for(u_int32_t nodeIndex : readpath){

				if(!_truthInputFilename.empty()){
					if(_evaluation_hifiasmGroundTruth_nodeNamePosition.find(BiGraph::nodeIndex_to_nodeName(nodeIndex)) == _evaluation_hifiasmGroundTruth_nodeNamePosition.end()){

						//cout << "XXXXX ";
					}
					else{
						//cout << BiGraph::nodeIndex_to_nodeName(nodeIndex) << " ";
					}
				}
			}
			//cout << endl;

			//cout << "error" << endl;
			//getchar();
		}
		*/

		//file_correction.flush();

		//346, 1734, 1770, 2136, 2323, 2364
		if(readIndex == 346){

			/*
			bool isHere = false;
			for(u_int32_t nodeIndex : _graph->_isNodeValid2){
				if(BiGraph::nodeIndex_to_nodeName(nodeIndex) == 0){
					isHere = true;
					break;
				}
			}
			cout << isHere << endl;
			*/
			//getchar();

			//file_correction.flush();
		}
		//getchar();

		//if(contigpath.size() == 0) return;
		//cout << "decommenter ici dans v final, contig size 0 a discard" << endl;

		//for(u_int32_t nodeIndex : contigpath){
		//	cout << BiGraph::nodeIndex_to_nodeName(nodeIndex) << " ";
		//}
		//cout << endl;
		//getchar();

		/*
		bool printRead = false;
		for(u_int32_t nodeIndex : contigpath){
			u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);

			if(_evaluation_hifiasmGroundTruth_nodeNamePosition.find(nodeName) == _evaluation_hifiasmGroundTruth_nodeNamePosition.end()) continue;

			if(_evaluation_hifiasmGroundTruth_nodeNamePosition[nodeName] == 4790){
				printRead = true;
			}
		}
		
		if(readIndex == 3427){
			cout << readIndex << ": " << endl;
			for(u_int32_t nodeIndex : _debug_readPath[readIndex]){
				cout << nodeIndex << " ";
			}
			cout << endl;
			
			for(u_int32_t nodeIndex : contigpath){
				cout << BiGraph::nodeIndex_to_nodeName(nodeIndex) << " ";
			}
			cout << endl;
			//getchar();
		}
		*/
		

		checkCorrectedSequence(readpath, print_read);

		//if(readIndex == 6) getchar();
		//getchar();
		
		u_int32_t readSize = readpath.size();
		_file_uncorrectedReads.write((const char*)&readSize, sizeof(readSize));
		_file_uncorrectedReads.write((const char*)&readpath[0], readSize*sizeof(u_int64_t));
		
		/*
		if(print_read){
			cout << readpath.size() << endl;
			cout << "\t";
			for(u_int32_t nodeIndex : readpath){
				cout << BiGraph::nodeIndex_to_nodeName(nodeIndex) << " ";
			}
			cout << endl;
		}

		_nbreadsLala += 1;

		if(readpath.size() < nodePath_withMissing.size() || nodePath.size() == 0 || readpath.size() == 0){
			u_int32_t readSize = minimizers.size();
			_file_uncorrectedReads.write((const char*)&readSize, sizeof(readSize));
			_file_uncorrectedReads.write((const char*)&minimizers[0], readSize*sizeof(u_int64_t));

			_nbUncorrectedReads += 1;
		}
		else{
			u_int64_t size = readpath.size();
			gzwrite(_outputFile_correctedReads, (const char*)&size, sizeof(size));
			gzwrite(_outputFile_correctedReads, (const char*)&readpath[0], size * sizeof(u_int32_t));
		}
		*/
		//vector<ReadIndexType> unitigIndexex;

		/*
		for(const KmerVec& vec : kminmers){
			//if(mdbg->_dbg_nodes[vec]._index == 55479) cout << "AAAAA" << endl;

			//cout << mdbg->_dbg_nodes[vec]._index << endl;
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()) continue;

			//if(vec.isPalindrome()) cout << _mdbg->_dbg_nodes[vec]._index << endl;
			//getchar();


			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
			//u_int32_t nodeIndex = BiGraph::nodeName_to_nodeIndex(nodeName, true);
			//if(_graph->_isNodeValid2.find(nodeIndex) == _graph->_isNodeValid2.end()) continue;
			//if(_nodeData.find(nodeName) == _nodeData.end()) continue;

			//UnitigData& unitigData = _nodeData[nodeName];
			UnitigData& unitigData = _unitigDatas[nodeName];
			//if(std::find(unitigData._readIndexes.begin(), unitigData._readIndexes.end(), readIndex) != unitigData._readIndexes.end()) continue;
			unitigData._readIndexes.push_back(readIndex);
			//cout << "indexing : " << readIndex << endl;
		}
		*/



	}

	void fillUncorrectedArea(u_int32_t position_source, u_int32_t nodeIndex_dest, const vector<u_int64_t>& readMinimizers, vector<u_int64_t>& readpath, const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, bool isFirstIteration, bool print_read){
		
		bool adding = false;


		for(size_t i=position_source; i<kminmers.size(); i++){

			const KmerVec& vec = kminmers[i];
			const ReadKminmer& info = kminmersInfos[i];

			

			if(!adding){
				
				u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
				u_int32_t nodeIndex = BiGraph::nodeName_to_nodeIndex(nodeName, !info._isReversed);

				if(readpath.size() == 0){
					if(isFirstIteration){
						for(size_t j=i; j<i+_kminmerSize; j++){
							u_int64_t minimizer = readMinimizers[j];
							if(print_read)  cout << "\tFill uncorrected area (start): " << j << " " << minimizer << endl;
							readpath.push_back(minimizer);
						}
					}
				}

				if(print_read)  cout << "Source position: " << nodeName << " " << i << " " << readMinimizers[i] << endl;
				adding = true;
			}
			else if(adding){
				u_int64_t minimizer = readMinimizers[i+_kminmerSize-1];
				readpath.push_back(minimizer);
				if(print_read)  cout << "\tFill uncorrected area: " << i << " " << minimizer << endl;
			}

			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()) continue;
			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
			u_int32_t nodeIndex = BiGraph::nodeName_to_nodeIndex(nodeName, !info._isReversed);

			if(nodeIndex == nodeIndex_dest) return;

		}

	}


	void applyReadCorrection(const vector<u_int32_t>& nodePath, vector<u_int64_t>& readpath, bool print_read, const vector<u_int64_t>& readMinimizers, const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, vector<u_int32_t>& nodePathSolid){
		


		float minSupportingReads = 2;
		u_int32_t maxDepth = nodePath.size()*2;

		readpath.clear();
		


		vector<u_int32_t> nodePath_errorFree;
		for(u_int32_t nodeIndex : nodePath){

			u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);

			if(_graph->_isNodeValid2.find(nodeIndex) == _graph->_isNodeValid2.end()) continue;
			
			//if(!_isBubble[nodeIndex]){
				//if(_graph->_isNodenameRoundabout.find(nodeName) != _graph->_isNodenameRoundabout.end()) continue;
			//}

			nodePath_errorFree.push_back(nodeIndex);
		}

		if(print_read){
			cout << "\tRead error free:" << endl;
			cout << "\t";
			for(u_int32_t nodeIndex : nodePath_errorFree){
				cout << BiGraph::nodeIndex_to_nodeName(nodeIndex) << " ";
			}
			cout << endl;
		}

		if(nodePath_errorFree.size() == 0) return;






		vector<u_int32_t> nodePath_errorFree_fixed;

		for(long i=0; i<((long)nodePath_errorFree.size())-1; i++){
			
			u_int32_t nodeIndex_source = nodePath_errorFree[i];
			u_int32_t nodeIndex_dest = nodePath_errorFree[i+1];

			nodePath_errorFree_fixed.push_back(nodeIndex_source);

			vector<u_int32_t> path;
			//_graph->shortestPath_nodeName(nodeName_source, nodeName_dest, path, true, true, maxDepth, _unitigDatas, unallowedNodeNames, nodeIndex_source, nodePath_anchor);
    		_graph->findUniquePath(nodeIndex_source, nodeIndex_dest, path, false, false, maxDepth);

			if(path.size() == 0){
				continue;
			}
			else{
				//isFix = true;
			}

			std::reverse(path.begin(), path.end());

			for(u_int32_t nodeIndex : path){
				nodePath_errorFree_fixed.push_back(nodeIndex);
			}
			
		}
		nodePath_errorFree_fixed.push_back(nodePath_errorFree[nodePath_errorFree.size()-1]);

		if(print_read){
			cout << "\tRead error free (fixed):" << endl;
			cout << "\t";
			for(u_int32_t nodeIndex : nodePath_errorFree_fixed){
				cout << BiGraph::nodeIndex_to_nodeName(nodeIndex) << " ";
			}
			cout << endl;
			//if(isFix) getchar();
		}


		//vector<u_int32_t> nodePath_anchor = nodePath_errorFree_fixed;
		//vector<u_int32_t> nodePath_connected = nodePath_errorFree_fixed;
		
		//u_int32_t prevNodename = -1;
		//u_int32_t nodeIndex_source = -1;
		/*
		//for(u_int32_t nodeName : nodePath_errorFree_fixed)



		
		//u_int32_t nodeIndex_dest = -1;

		for(long i=0; i<((long)nodePath_connected.size())-1; i++){
			
			u_int32_t nodeName_source = nodePath_connected[i];
			u_int32_t nodeName_dest = nodePath_connected[i+1];


			if(nodeIndex_source != -1){
				if(tryFindPathDirect(nodeIndex_source, nodeName_dest, readpath)){
					nodeIndex_source = readpath[readpath.size()-1];
					continue;
				}
			}




			bool foundPath = false;


			priority_queue<NodeNameAbundance, vector<NodeNameAbundance> , NodeNameAbundance_Comparator> queue;
			for(u_int32_t nodeName : nodePath_errorFree_fixed){
				queue.push({nodeName, _graph->getNodeUnitigAbundance(BiGraph::nodeName_to_nodeIndex(nodeName, true))});
			}

			while(!queue.empty()){

				NodeNameAbundance nodeNameAbundance = queue.top();
				queue.pop();

				u_int32_t anchorNodeName = nodeNameAbundance._nodeName;

				unordered_set<u_int32_t> unallowedNodeNames;
				for(long j=0; j<nodePath_connected.size(); j++){
					if(nodePath_connected[j] == nodeName_source || nodePath_connected[j] == nodeName_dest) continue;
					unallowedNodeNames.insert(nodePath_connected[j]);
				}

				//cout << "\tSearching path: " << nodeName_source << " " << nodeName_dest << endl;
				//cout << Utils::computeSharedReads(_unitigDatas[nodeName_source], _unitigDatas[nodeName_dest]) << endl;

				//cout << nodePath_anchor.size() << endl;
				vector<u_int32_t> path;
				_graph->shortestPath_nodeName(nodeName_source, nodeName_dest, path, true, true, maxDepth, _unitigDatas, unallowedNodeNames, nodeIndex_source, anchorNodeName);

				if(path.size() == 0){
					foundPath = false;
					break;
				}

				foundPath = true;
				//cout << "done" << endl;

				std::reverse(path.begin(), path.end());

				//cout << "\t";
				if(readpath.size() == 0){
					for(u_int32_t nodeIndex : path){
						readpath.push_back(nodeIndex);
					}
				}
				else{ //Discard source
					for(size_t i=1; i<path.size(); i++){ 
						//cout << BiGraph::nodeIndex_to_nodeName(nodeIndex) << " ";
						readpath.push_back(path[i]);
					}
				}
				prevNodename = nodeName_source;
				
				if(readpath.size() > 0){
					nodeIndex_source = readpath[readpath.size()-1];
				}
				//cout << endl;

				if(foundPath){
					//nodePath_anchor = nodePath_errorFree_fixed;
					break;
				}

			}


		}

		*/

		u_int32_t nodeIndex_source;
		u_int32_t nodeIndexFinal = -1;

		unordered_set<u_int32_t> isNodenameSolid;
		for(u_int32_t nodeIndex : nodePath_errorFree_fixed){
			isNodenameSolid.insert(BiGraph::nodeIndex_to_nodeName(nodeIndex));
		}

		vector<u_int32_t> readpath_nodes;
		bool isLastSuccess = false;
		size_t i=0;
		bool isFirstIteration = true;

		for(size_t i=0; i<kminmers.size()-1; i++){

			const KmerVec& vec = kminmers[i];
			const ReadKminmer& info = kminmersInfos[i];

			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()) continue;
			
			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
			if(isNodenameSolid.find(nodeName) == isNodenameSolid.end()) continue;

			u_int32_t nodeIndex_source = BiGraph::nodeName_to_nodeIndex(nodeName, !info._isReversed);

			//for(long i=0; i<((long)nodePath_connected.size())-1; i++){
			
			//u_int32_t nodeName_source = nodeName; //nodePath_connected[i];
			u_int32_t nodeIndex_dest = -1;
			for(size_t j=i+1; j<kminmers.size(); j++){
				const KmerVec& vec = kminmers[j];
				const ReadKminmer& info = kminmersInfos[j];
				if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()) continue;
				u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
				if(isNodenameSolid.find(nodeName) == isNodenameSolid.end()) continue;
				nodeIndex_dest = BiGraph::nodeName_to_nodeIndex(nodeName, !info._isReversed);
				break;
			}

			if(print_read)  cout << BiGraph::nodeIndex_to_nodeName(nodeIndex_source) << " -> " << BiGraph::nodeIndex_to_nodeName(nodeIndex_dest) << endl;

			if(tryFindPathDirect(nodeIndex_source, nodeIndex_dest, readpath_nodes)){
				isLastSuccess = true;
				


				if(readpath.size() == 0){

					bool orientation;
					u_int32_t nodeName_start = BiGraph::nodeIndex_to_nodeName(nodeIndex_source, orientation);
					vector<u_int64_t> minimizerSeq = _nodeName_to_kminmerSequence[nodeName_start];

					if(orientation){
						if(isFirstIteration){
							for(u_int64_t m : minimizerSeq){
								readpath.push_back(m);
								if(print_read)  cout << "Add: " << m << endl;
							}
						}
						else{
							if(print_read)  cout << "Add: " << minimizerSeq[minimizerSeq.size()-1] << endl;
							readpath.push_back(minimizerSeq[minimizerSeq.size()-1]);
						}
					}
					else{
						if(isFirstIteration){
							std::reverse(minimizerSeq.begin(), minimizerSeq.end());
							for(u_int64_t m : minimizerSeq){
								readpath.push_back(m);
								if(print_read)  cout << "Add: " << m << endl;
							}
						}
						else{
							if(print_read)  cout << "Add: " << minimizerSeq[0] << endl;
							readpath.push_back(minimizerSeq[0]);
						}

					}
				}


				bool orientation;
				u_int32_t nodeName_dest = BiGraph::nodeIndex_to_nodeName(nodeIndex_dest, orientation);
				vector<u_int64_t> minimizerSeq = _nodeName_to_kminmerSequence[nodeName_dest];

				if(orientation){
					if(print_read)  cout << "Add: " << minimizerSeq[minimizerSeq.size()-1] << endl;
					readpath.push_back(minimizerSeq[minimizerSeq.size()-1]);
				}
				else{
					if(print_read)  cout << "Add: " << minimizerSeq[0] << endl;
					readpath.push_back(minimizerSeq[0]);
				}

				//lastSolidNodeIndex = nodeIndex_source;

				nodeIndex_source = nodeIndex_dest; //readpath_nodes[readpath_nodes.size()-1];
				nodeIndexFinal = nodeIndex_dest;


				isFirstIteration = false;
				//continue;
			}
			else{
				if(print_read) cout << "failed" << endl;
				//lastSolidNodeIndex = -1;
				isLastSuccess = false;
				fillUncorrectedArea(i, nodeIndex_dest, readMinimizers, readpath, kminmers, kminmersInfos, isFirstIteration, print_read);
				nodeIndexFinal = nodeIndex_dest;
				isFirstIteration = false;
				//getchar();
			}


			/*
			if(print_read){
				cout << "\t\tIter " << i << " " << nodeName_source << " " << nodeName_dest << endl; 
				cout << "\t\t";
				for(u_int64_t minimizer : readpath){
					cout << minimizer << " ";
				}
				cout << endl;
			}
			*/

			/*
			unordered_set<u_int32_t> unallowedNodeNames;
			for(long j=0; j<nodePath_connected.size(); j++){
				if(nodePath_connected[j] == nodeName_source || nodePath_connected[j] == nodeName_dest) continue;
				unallowedNodeNames.insert(nodePath_connected[j]);
			}

			bool foundPath = false;

			//while(true){



				//cout << "\tSearching path: " << nodeName_source << " " << nodeName_dest << endl;
				//cout << Utils::computeSharedReads(_unitigDatas[nodeName_source], _unitigDatas[nodeName_dest]) << endl;

				//cout << nodePath_anchor.size() << endl;
				vector<u_int32_t> path;
				_graph->shortestPath_nodeName(nodeName_source, nodeName_dest, path, true, true, maxDepth, _unitigDatas, unallowedNodeNames, nodeIndex_source, nodePath_anchor);

				if(path.size() == 0){
					//foundPath = false;
					//break;
				}

				foundPath = true;
				//cout << "done" << endl;

				std::reverse(path.begin(), path.end());

				//cout << "\t";
				if(readpath.size() == 0){
					for(u_int32_t nodeIndex : path){
						readpath.push_back(nodeIndex);
					}
				}
				else{ //Discard source
					for(size_t i=1; i<path.size(); i++){ 
						//cout << BiGraph::nodeIndex_to_nodeName(nodeIndex) << " ";
						readpath.push_back(path[i]);
					}
				}
				prevNodename = nodeName_source;
				
				if(readpath.size() > 0){
					nodeIndex_source = readpath[readpath.size()-1];
				}
				//cout << endl;

				

			//}
			*/

		}
		
		
		nodePathSolid.clear();
		nodePathSolid.push_back(nodePath_errorFree_fixed[0]);
		if(nodeIndexFinal != -1){
			nodePathSolid.push_back(nodeIndexFinal);
		}

		/*
		if(isLastSuccess && nodeIndex_source != -1){
			bool orientation;
			u_int32_t nodeName_start = BiGraph::nodeIndex_to_nodeName(nodeIndex_source, orientation);
			vector<u_int64_t> minimizerSeq = _nodeName_to_kminmerSequence[nodeName_start];

			if(orientation){
				readpath.push_back(minimizerSeq[minimizerSeq.size()-1]);
				cout << "Add: " << minimizerSeq[minimizerSeq.size()-1] << endl;
			}
			else{
				readpath.push_back(minimizerSeq[0]);
				cout << "Add: " << minimizerSeq[0] << endl;
			}
		}
		*/

		//nodePathSolid = readpath_nodes;

		/*
		if(print_read){
			cout << "\tRead solid:" << endl;
			cout << "\t";
			for(u_int32_t nodeIndex : readpath_nodes){
				cout << BiGraph::nodeIndex_to_nodeName(nodeIndex) << " ";
			}
			cout << endl;
			//if(isFix) getchar();
		}*/

		if(print_read){
			cout << "\tRead corrected (minimizers):" << endl;
			cout << "\t";
			for(u_int64_t minimizer : readpath){
				cout << minimizer << " ";
			}
			cout << endl;



		}
/*
				if(foundPath){
					nodePath_anchor = nodePath_errorFree_fixed;
					break;
				}
				else{
					size_t distStart = 0;
					size_t distEnd = 0;
					for(u_int32_t nodeName : nodePath_connected){
						if(nodeName == nodeName_source || nodeName == nodeName_dest) break;
						distStart += 1;
					}
					for(long i=nodePath_connected.size()-1; i>=0; i--){
						u_int32_t nodeName = nodePath_connected[i];
						if(nodeName == nodeName_source || nodeName == nodeName_dest) break;
						distEnd += 1;
					}

					bool isSolution = true;

					if(distStart < distEnd){
						u_int32_t nodeName = nodePath_anchor[0];
						if(nodeName == nodeName_source || nodeName == nodeName_dest){
							nodeName = nodePath_anchor[nodePath_anchor.size()-1];
							if(nodeName == nodeName_source || nodeName == nodeName_dest){
								isSolution = false;
							}
							else{
								nodePath_anchor.pop_back();
							}
						}
						else{
							nodePath_anchor.erase(nodePath_anchor.begin());
						}
					}
					else{
						u_int32_t nodeName = nodePath_anchor[nodePath_anchor.size()-1];
						if(nodeName == nodeName_source || nodeName == nodeName_dest){
							nodeName = nodePath_anchor[0];
							if(nodeName == nodeName_source || nodeName == nodeName_dest){
								isSolution = false;
							}
							else{
								nodePath_anchor.pop_back();
							}
						}
						else{
							nodePath_anchor.erase(nodePath_anchor.begin());
						}
					}

					if(!isSolution){
						nodePath_anchor = nodePath_errorFree_fixed;
						break;
					}
				}*/
		



		//if(readpath.size() < nodePath.size())
		
		/*
		long supportSize = 6;
		if(nodePath_existing.size() < supportSize) return;

		u_int64_t anchorPos = -1;
		vector<u_int32_t> nodePath_anchor;
		for(long i=0; i< ((long)nodePath_existing.size())-supportSize; i++){
			
			vector<u_int32_t> nodeNames;
			for(size_t j=i; j<i+supportSize; j++){
				nodeNames.push_back(nodePath_existing[j]);
			}
			//u_int32_t nodeName_source = nodePath_existing[i];
			//u_int32_t nodeName_dest = nodePath_existing[i+supportSize];

			u_int64_t nbSupportingReads = Utils::computeSharedReads(nodeNames, _unitigDatas);
			//cout << i << " " << nbSupportingReads << endl;
			//cout << nodeName_source << " " << nodeName_dest << " " << nbSupportingReads << endl;
			
			if(nbSupportingReads < minSupportingReads) continue;

			for(size_t j=i; j<i+supportSize; j++){
				nodePath_anchor.push_back(nodePath_existing[j]);
			}

			anchorPos = i;
			break;
		}

		
		if(anchorPos == -1) return;

		//vector<u_int32_t> nodePath_supported_left;

		//LEFT
		long ii = anchorPos - 1;
		while(true){

			if(ii < 0) break;

			vector<u_int32_t> nodeNames;
			for(size_t j=ii; j<ii+supportSize; j++){
				nodeNames.push_back(nodePath_existing[j]);
			}
			
			u_int64_t nbSupportingReads = Utils::computeSharedReads(nodeNames, _unitigDatas);

			//cout << "left " << ii << " " << nbSupportingReads << endl;

			if(nbSupportingReads < minSupportingReads){
				nodePath_existing.erase(nodePath_existing.begin()+ii);
			}

			ii -= 1;
		}

		//RIGHT
		ii = anchorPos + 1;
		while(true){

			if(ii > ((long)nodePath_existing.size())-supportSize) break;

			vector<u_int32_t> nodeNames;
			for(size_t j=ii; j<ii+supportSize; j++){
				nodeNames.push_back(nodePath_existing[j]);
			}
			
			u_int64_t nbSupportingReads = Utils::computeSharedReads(nodeNames, _unitigDatas);

			//cout << "right " << ii << " " << nbSupportingReads << endl;

			if(nbSupportingReads < minSupportingReads){
				nodePath_existing.erase(nodePath_existing.begin()+ii+supportSize-1);
				continue; //! the ii incrmeent is handled by erase
			}

			ii += 1;
		}
		*/
		
		/*
		cout << "\tRead supported:" << endl;
		cout << "\t";
		for(u_int32_t nodeName : nodePath_existing){
			cout << nodeName << " ";
		}
		cout << endl;
		*/
		/*
		bool hasFoundAnchor = false;
		while(true){

			cout << "----" << endl;
			size_t j=0;
			bool isUnsupported = false;
			nodePath_supported.clear();

			for(long i=0; i< ((long)nodePath_existing.size())-supportSize; i++){
				
				j = i;
				u_int32_t nodeName_source = nodePath_existing[i];
				u_int32_t nodeName_dest = nodePath_existing[i+supportSize];

				u_int64_t nbSupportingReads = Utils::computeSharedReads(_unitigDatas[nodeName_source], _unitigDatas[nodeName_dest]);
				cout << nodeName_source << " " << nodeName_dest << " " << nbSupportingReads << endl;
				
				if(nbSupportingReads < 3){ 
					isUnsupported = true;
					break;
				}

				hasFoundAnchor = true;

				if(nodePath_supported.size() == 0){
					for(size_t j=i; j<i+supportSize; j++){
						nodePath_supported.push_back(nodePath_existing[j]);
					}
				}
				
				nodePath_supported.push_back(nodeName_dest);
				//cout << "\tSupporting reads: " << nodeName_source << " " << nodeName_dest << endl;
				//cout << "\t" << Utils::computeSharedReads(_unitigDatas[nodeName_source], _unitigDatas[nodeName_dest]) << endl;
			}

			if(isUnsupported){
				if(hasFoundAnchor){
					cout << "Removed: " << nodePath_existing[j+supportSize] << endl;
					nodePath_existing.erase(nodePath_existing.begin() + j + supportSize);
				}
				else{
					cout << "Removed: " << nodePath_existing[j] << endl;
					nodePath_existing.erase(nodePath_existing.begin() + j);
				}
			}
			else{
				break;
			}
		}
		*§
		/*
		for(long i=0; i< ((long)nodePath_existing.size())-supportSize; i++){
			
			u_int32_t nodeName_source = nodePath_existing[i];
			u_int32_t nodeName_dest = nodePath_existing[i+supportSize];

			u_int64_t nbSupportingReads = Utils::computeSharedReads(_unitigDatas[nodeName_source], _unitigDatas[nodeName_dest]);
			if(nbSupportingReads < 3) continue;

			if(nodePath_supported.size() == 0){
				for(size_t j=i; j<i+supportSize; j++){
					nodePath_supported.push_back(nodePath_existing[j]);
				}
			}
			
			nodePath_supported.push_back(nodeName_dest);
			//cout << "\tSupporting reads: " << nodeName_source << " " << nodeName_dest << endl;
			//cout << "\t" << Utils::computeSharedReads(_unitigDatas[nodeName_source], _unitigDatas[nodeName_dest]) << endl;
		}
		*/

		/*
		cout << "\tRead supported:" << endl;
		cout << "\t";
		for(u_int32_t nodeName : nodePath_supported){
			cout << nodeName << " ";
		}
		cout << endl;




		//cout << "Computing components" << endl;
		u_int32_t componentIndex = 0;
		//u_int32_t currentComponentIndex = 0;
		unordered_map<u_int32_t, u_int32_t> nodeName_to_component;
		//unordered_set<u_int32_t> isConnected;
		//vector<vector<u_int32_t>> components;
		
		for(size_t i=0; i<nodePath_supported.size(); i++){

			u_int32_t nodeName_source = nodePath_supported[i];

			if(nodeName_to_component.find(nodeName_source) == nodeName_to_component.end()){
				nodeName_to_component[nodeName_source] = componentIndex;
				componentIndex += 1;
			}

			u_int32_t currentComponentIndex = nodeName_to_component[nodeName_source];
			//currentComponentIndex = componentIndex;

			//if(isConnected.find(nodeName_source) != isConnected.end()) continue;

			//vector<u_int32_t> component;
			//component.push_back(nodeName_source);
			//isConnected.insert(nodeName_source);

			for(size_t j=i+1; j<nodePath_supported.size(); j++){
			
				u_int32_t nodeName_dest = nodePath_supported[j];

				//if(nodeName_to_component.find(nodeName_dest) != nodeName_to_component.end()){
				//	if(nodeName_to_component[nodeName_dest] == componentIndex) continue;
				//}

				//if(isConnected.find(nodeName_dest) != isConnected.end()) continue;

				//cout << "\tSearching path: " << nodeName_source << " " << nodeName_dest << endl;
				vector<u_int32_t> path;
				_graph->shortestPath_nodeName(nodeName_source, nodeName_dest, path, true, true, 20);
				//cout << "done" << endl;

				if(path.size() == 0) continue;
				
				nodeName_to_component[nodeName_dest] = currentComponentIndex;

				//component.push_back(nodeName_dest);
				//isConnected.insert(nodeName_dest);
				
			}

			//components.push_back(component);
		}

		unordered_map<u_int32_t, vector<u_int32_t>> components;
		for(auto& it : nodeName_to_component){
			components[it.second].push_back(it.first);
		}

		vector<u_int32_t> maxComponent;
		u_int32_t maxComponentSize = 0;
		if(components.size() != 1){
			cout << "\tComponents: " << components.size() << endl;
			for(auto& it : components){
				const vector<u_int32_t>& c = it.second;
				cout << "\t\t";
				for(u_int32_t nodeName : c) cout << nodeName << " ";
				cout << endl;
			}
			//getchar();
		} 
		
		for(auto& it : components){
			//cout << "\t\t" << c.size() << endl;
			const vector<u_int32_t>& c = it.second;

			if(c.size() > maxComponentSize){
				maxComponentSize = c.size();
				maxComponent = c;
			}
		}

		unordered_set<u_int32_t> componentIndexed;
		for(u_int32_t nodeName : maxComponent) componentIndexed.insert(nodeName);
		
		vector<u_int32_t> nodePath_connected;
		for(u_int32_t nodeName : nodePath_supported){
			if(componentIndexed.find(nodeName) == componentIndexed.end()) continue;
			nodePath_connected.push_back(nodeName);
		}
		*/


		/*
		vector<u_int32_t> nodePath_connected = nodePath_existing;
		u_int32_t prevNodename = -1;
		u_int32_t nodeIndex_source = -1;
		//u_int32_t nodeIndex_dest = -1;

		for(long i=0; i<((long)nodePath_connected.size())-1; i++){
			
			u_int32_t nodeName_source = nodePath_connected[i];
			u_int32_t nodeName_dest = nodePath_connected[i+1];

			unordered_set<u_int32_t> unallowedNodeNames;
			for(long j=0; j<nodePath_connected.size(); j++){
				if(nodePath_connected[j] == nodeName_source || nodePath_connected[j] == nodeName_dest) continue;
				unallowedNodeNames.insert(nodePath_connected[j]);
			}

			//cout << "\tSearching path: " << nodeName_source << " " << nodeName_dest << endl;
			//cout << Utils::computeSharedReads(_unitigDatas[nodeName_source], _unitigDatas[nodeName_dest]) << endl;

			vector<u_int32_t> path;
			_graph->shortestPath_nodeName(nodeName_source, nodeName_dest, path, true, true, maxDepth, _unitigDatas, unallowedNodeNames, nodeIndex_source, nodePath_anchor);


			//cout << "done" << endl;

			std::reverse(path.begin(), path.end());

			//cout << "\t";
			if(readpath.size() == 0){
				for(u_int32_t nodeIndex : path){
					readpath.push_back(nodeIndex);
				}
			}
			else{ //Discard source
				for(size_t i=1; i<path.size(); i++){ 
					//cout << BiGraph::nodeIndex_to_nodeName(nodeIndex) << " ";
					readpath.push_back(path[i]);
				}
			}
			prevNodename = nodeName_source;
			
			if(readpath.size() > 0){
				nodeIndex_source = readpath[readpath.size()-1];
			}
			//cout << endl;
		}
		//cout << endl;
		*/
		/*
		cout << nodePath_connected.size() << " " << nodePath_connected.size()-6 << endl;
		for(long i=0; i< ((long)nodePath_connected.size())-6; i++){
			
			u_int32_t nodeName_source = nodePath_connected[i];
			u_int32_t nodeName_dest = nodePath_connected[i+6];
			cout << "\tSupporting reads: " << nodeName_source << " " << nodeName_dest << endl;
			cout << "\t" << Utils::computeSharedReads(_unitigDatas[nodeName_source], _unitigDatas[nodeName_dest]) << endl;
		}
		*/

	}



	bool tryFindPathDirect(u_int32_t nodeIndex_source, u_int32_t nodeIndex_dest, vector<u_int32_t>& readpath){

		u_int32_t unitigIndex_source = _graph->nodeIndex_to_unitigIndex(nodeIndex_source);

		vector<u_int32_t> successors;
		_graph->getSuccessors(nodeIndex_source, 0, successors);

		for(u_int32_t nodeIndexSuccessor: successors){

			u_int32_t unitigIndex = _graph->nodeIndex_to_unitigIndex(nodeIndexSuccessor);
			if(unitigIndex != unitigIndex_source) continue;

			if(nodeIndexSuccessor == nodeIndex_dest) return true;
			/*
			if(BiGraph::nodeIndex_to_nodeName(nodeIndexSuccessor) == nodeName_dest){
				if(readpath.size() == 0){
					readpath.push_back(nodeIndex_source);
					readpath.push_back(nodeIndexSuccessor);
					return true;
				}
				else{
					readpath.push_back(nodeIndexSuccessor);
					return true;
				}
			}
			*/

		}

		return false;
	}

	void extendReadpath(const vector<u_int64_t>& readMinimizers, const vector<KmerVec>& kminmers, const vector<u_int32_t>& nodePathSolid, vector<u_int64_t>& readpath, bool print_read){


		if(nodePathSolid.size() == 0) return;

		unordered_set<u_int32_t> unallowedNodeNames;
		for(u_int32_t nodeIndex : nodePathSolid){
			unallowedNodeNames.insert(BiGraph::nodeIndex_to_nodeName(nodeIndex));
		}


		if(print_read){
			cout << "\tExtending read path" << endl;
			cout << "\t";
			for(u_int32_t nodeIndex : nodePathSolid){
				cout << BiGraph::nodeIndex_to_nodeName(nodeIndex) << " ";
			}
			cout << endl;
		}

		//if(readpath.size() == 0) return;

		//cout << "\tNode path original:" << endl;
		//cout << "\t";
		//for(u_int32_t nodeName : originalNodePath){
		//	cout << nodeName << " ";
		//}
		//cout << endl;
		/*
		unordered_set<u_int32_t> readPathIndexed;
		for(size_t i=0; i<readpath.size(); i++){
			//cout << BiGraph::nodeIndex_to_nodeName(readpath[i]) << endl;
			readPathIndexed.insert(BiGraph::nodeIndex_to_nodeName(readpath[i]));
		}

		//cout << "lala" << endl;

		long nbExtendLeft = 0;
		for(size_t i=0; i<originalNodePath.size(); i++){
			//cout << originalNodePath[i] << " " << (readPathIndexed.find(originalNodePath[i]) != readPathIndexed.end()) << endl;
			if(readPathIndexed.find(originalNodePath[i]) == readPathIndexed.end()){
				nbExtendLeft += 1;
			}
			else{
				break;
			}
		}

		long nbExtendRight = 0;
		for(long i=originalNodePath.size()-1; i >= 0; i--){
			if(readPathIndexed.find(originalNodePath[i]) == readPathIndexed.end()){
				nbExtendRight += 1;
			}
			else{
				break;
			}
		}

		nbExtendLeft += 2;
		nbExtendRight += 2;
		*/

		u_int32_t nodeNameSolidLeft = BiGraph::nodeIndex_to_nodeName(nodePathSolid[0]);

		u_int64_t solidPositionLeft = 0;
		for(long i=0; i < kminmers.size(); i++){
			const KmerVec& vec = kminmers[i];

			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()){
				continue;
			}
			
			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;

			if(nodeName == nodeNameSolidLeft){
				
				solidPositionLeft = i;
				break;
			}
		}

		u_int32_t nodeNameSolidRight = BiGraph::nodeIndex_to_nodeName(nodePathSolid[nodePathSolid.size()-1]);

		u_int64_t solidPositionRight = 0;
		for(long i=kminmers.size()-1; i>=0; i--){
			const KmerVec& vec = kminmers[i];

			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()){
				continue;
			}
			
			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;

			if(nodeName == nodeNameSolidRight){
				
				solidPositionRight = i;
				break;
			}
		}


		u_int64_t nbExtendLeft = solidPositionLeft;
		u_int64_t nbExtendRight = readMinimizers.size() - (solidPositionRight + _kminmerSize - 1);
		if(print_read){
			cout << "\tLeft solid position: " << nodeNameSolidLeft << " " << solidPositionLeft << endl;
			cout << "\tRight solid position: " << nodeNameSolidRight << " " << solidPositionRight << endl;
		}

		nbExtendLeft += 2;
		nbExtendRight += 2;
		
		if(print_read) cout << "\tExtend: " << nbExtendLeft << " " << nbExtendRight << endl;

		solidPositionLeft -= 1;
		solidPositionRight += _kminmerSize;

		extendReadpath2(readpath, nodePathSolid, true, nbExtendRight, unallowedNodeNames, solidPositionRight, readMinimizers, print_read);
		extendReadpath2(readpath, nodePathSolid, false, nbExtendLeft, unallowedNodeNames, solidPositionLeft, readMinimizers, print_read);
		
		/*
		vector<u_int32_t> prevNodes;

		contigpath.clear();
		if(readpath.size() == 0) return;

		u_int32_t nodeIndexStart = readpath[0];
		const Unitig& unitigStart = _graph->nodeIndex_to_unitig(nodeIndexStart);

		u_int32_t prevNode = -1;
		for(u_int32_t nodeIndex : unitigStart._nodes){
			if(nodeIndex == nodeIndexStart) break;
			prevNode = nodeIndex;
			prevNodes.push_back(nodeIndex);
		}
		//if(prevNode != -1) contigpath.push_back(prevNode);
		
		
		long size = prevNodes.size();
		long start_i = max(0l, size-n);
		for(size_t i=start_i; i<prevNodes.size(); i++){
			contigpath.push_back(prevNodes[i]);
		}
		
		for(u_int32_t nodeIndex : readpath){
			contigpath.push_back(nodeIndex);
		}

		u_int32_t nodeIndexEnd = readpath[readpath.size()-1];
		const Unitig& unitigEnd = _graph->nodeIndex_to_unitig(nodeIndexEnd);

		vector<u_int32_t> nextNodes;
		u_int32_t nextNode = -1;
		bool canAdd = false;
		for(u_int32_t nodeIndex : unitigEnd._nodes){
			if(nodeIndex == nodeIndexEnd){
				canAdd = true;
				continue;
			}
			if(!canAdd) continue;
			nextNodes.push_back(nodeIndex);
			//nextNode = nodeIndex;
			//break;
		}
		
		//if(nextNode != -1) contigpath.push_back(nextNode);
		for(size_t i=0; i<n && i<nextNodes.size(); i++){
			contigpath.push_back(nextNodes[i]);
		}
		*/

	}

	struct SuccessorData3{
		u_int32_t _nodeIndex;
		bool _isSupported;
	};

	void extendReadpath2(vector<u_int64_t>& readpath, const vector<u_int32_t>& nodePath, bool forward, size_t nbSuccessors, unordered_set<u_int32_t>& unallowedNodeNames, u_int64_t extendPosition, const vector<u_int64_t>& readMinimizers, bool print_read){

		if(print_read){
			if(forward){
				cout << "\tExtending right" << endl;
			}
			else{
				cout << "\tExtending left" << endl;
			}
		}

		if(nodePath.size() == 0) return;

		u_int32_t nodeIndex = -1;
		if(forward){
			nodeIndex = nodePath[nodePath.size()-1];
		}
		else{
			nodeIndex = nodePath[0];
		}

		//bool orientation;
		//u_int32_t nodeName_start = BiGraph::nodeIndex_to_nodeName(nodeIndex, orientation);
		//cout << nodeName_start << endl;
		//vector<u_int64_t> minimizerSeq = _nodeName_to_kminmerSequence[nodeName_start];
		//cout << orientation << endl;
		//for(u_int64_t minimizer : minimizerSeq){
		//	cout << minimizer << " ";
		//}
		//cout << endl;

		bool failed = false;

		for(size_t i=0; i<nbSuccessors; i++){

			vector<u_int32_t> successors;
			if(forward){
				_graph->getSuccessors(nodeIndex, 0, successors);
			}
			else{
				_graph->getPredecessors(nodeIndex, 0, successors);
			}


			if(successors.size() == 0){
				failed = true;
				break;
			}
			else if(successors.size() == 1){
				u_int32_t nodeIndexSuccessor = successors[0];
				if(unallowedNodeNames.find(BiGraph::nodeIndex_to_nodeName(nodeIndexSuccessor)) != unallowedNodeNames.end()){
					failed = true;
					break;
				}
				if(_graph->nodeIndex_to_unitigIndex(nodeIndex) != _graph->nodeIndex_to_unitigIndex(nodeIndexSuccessor)){
					failed = true;
					break;
				}
				
				
				bool orientation;
				u_int32_t nodeName_start = BiGraph::nodeIndex_to_nodeName(nodeIndexSuccessor, orientation);
				vector<u_int64_t> minimizerSeq = _nodeName_to_kminmerSequence[nodeName_start];
				if(!forward){
					std::reverse(minimizerSeq.begin(), minimizerSeq.end());
				}

				u_int64_t minimizer = -1;
				if(orientation){
					minimizer = minimizerSeq[minimizerSeq.size()-1];
				}
				else{
					minimizer = minimizerSeq[0];
				}


				//cout << "Add: " << nodeName_start << " " <<  minimizer << endl;
				//for(u_int64_t m : minimizerSeq){
				//	cout << m << " ";
				//}
				//cout << endl;

				if(forward){
					readpath.push_back(minimizer);
					extendPosition += 1;
					//nodePath.push_back(nodeIndexSuccessor);
				}
				else{
					//for(u_int64_t m : minimizerSeq){
					//	cout << m <<  " ";
					//}
					//cout << endl;
					readpath.insert(readpath.begin(), minimizer);
					extendPosition -= 1;
					//nodePath.insert(nodePath.begin(), nodeIndexSuccessor);
				}

				nodeIndex = nodeIndexSuccessor;

			}
			else{
				failed = true;
				break;
				/*
				nodeIndex = determineBestSupportedSuccessors(nodePath, forward, unallowedNodeNames, successors);
				if(nodeIndex == -1) return;

				if(unallowedNodeNames.find(BiGraph::nodeIndex_to_nodeName(nodeIndex)) != unallowedNodeNames.end()) return;

				if(forward){
					nodePath.push_back(nodeIndex);
				}
				else{
					nodePath.insert(nodePath.begin(), nodeIndex);
				}
				*/
			}
			
		}
		//vector<u_int32_t> prevNodes;
		//if()

		/*
		if(failed){
			

			cout << "failed" << endl;
			if(forward){
				cout << extendPosition << " " << readMinimizers.size() << endl;
				//cout << nbRemainingExtension << " " << extendPosition << endl;
				for(size_t i=extendPosition; i<readMinimizers.size(); i++){
					readpath.push_back(readMinimizers[i]);
					cout << "Add: " << i << " " << readMinimizers[i] << endl;
				}
			}
			else{
				//cout << nbRemainingExtension << " " << extendPosition << endl;
				for(long i=extendPosition; i>=0; i--){
					readpath.insert(readpath.begin(), readMinimizers[i]);
					cout << "Add: " << i << " " << readMinimizers[i] << endl;
				}
			}
		}
		*/
	}

	u_int32_t determineBestSupportedSuccessors(vector<u_int32_t>& nodePath, bool forward, unordered_set<u_int32_t>& unallowedNodeNames, const vector<u_int32_t>& successors){

		#ifdef PRINT_DEBUG_PREVRANK
			cout << "\t-------------------------" << endl;
		#endif
		
		vector<SuccessorData3> possibleSuccessors;
		for(u_int32_t nodeIndex_successor : successors){
			possibleSuccessors.push_back({nodeIndex_successor, true});
		}

		size_t prevRank = 0;

		while(true){


			long prevRankIndex = -1; 

			if(forward){
				prevRankIndex = nodePath.size() - prevRank - 1;
				if(prevRankIndex < 0) return -1;
			}
			else{
				prevRankIndex = prevRank;
				if(prevRankIndex >= nodePath.size()) return -1;
			}

			//cout << "\t" << prevRankIndex << endl;
			u_int32_t nodeIndex_prev = nodePath[prevRankIndex];
			u_int32_t nodeName_prev = BiGraph::nodeIndex_to_nodeName(nodeIndex_prev);

			#ifdef PRINT_DEBUG_PREVRANK
				cout << "\t" << prevRank << " " << nodeName_prev << "    ";
			#endif

			bool isFinished = false;
			for(SuccessorData3& successorData : possibleSuccessors){
				if(!successorData._isSupported) continue;

				u_int32_t nodeName_successor = BiGraph::nodeIndex_to_nodeName(successorData._nodeIndex);

				u_int64_t nbShared = Utils::computeSharedReads(_unitigDatas[nodeName_successor], _unitigDatas[nodeName_prev]);
				
				#ifdef PRINT_DEBUG_PREVRANK
					cout << nodeName_successor << " " << nbShared << "    ";
				#endif

				if(nbShared <= 2){
					successorData._isSupported = false;
				}

			}

			#ifdef PRINT_DEBUG_PREVRANK
				cout << endl;
			#endif

			u_int32_t supporteNodeIndex = -1;
			size_t nbSupportedSuccessors = 0;
			for(SuccessorData3& successorData : possibleSuccessors){
				if(successorData._isSupported){
					nbSupportedSuccessors += 1;
					supporteNodeIndex = successorData._nodeIndex;
				}
			}

			if(nbSupportedSuccessors == 1){
				return supporteNodeIndex;
			}

			prevRank += 1;

			if(prevRank > 10) return -1;
		}

		return -1;
	}

	/*
	void readpathToContig(const vector<u_int32_t>& readpath, vector<u_int32_t>& contigpath){

		long n = 2;
		vector<u_int32_t> prevNodes;

		contigpath.clear();
		if(readpath.size() == 0) return;

		u_int32_t nodeIndexStart = readpath[0];
		const Unitig& unitigStart = _graph->nodeIndex_to_unitig(nodeIndexStart);

		u_int32_t prevNode = -1;
		for(u_int32_t nodeIndex : unitigStart._nodes){
			if(nodeIndex == nodeIndexStart) break;
			prevNode = nodeIndex;
			prevNodes.push_back(nodeIndex);
		}
		//if(prevNode != -1) contigpath.push_back(prevNode);
		
		
		long size = prevNodes.size();
		long start_i = max(0l, size-n);
		for(size_t i=start_i; i<prevNodes.size(); i++){
			contigpath.push_back(prevNodes[i]);
		}
		
		for(u_int32_t nodeIndex : readpath){
			contigpath.push_back(nodeIndex);
		}

		u_int32_t nodeIndexEnd = readpath[readpath.size()-1];
		const Unitig& unitigEnd = _graph->nodeIndex_to_unitig(nodeIndexEnd);

		vector<u_int32_t> nextNodes;
		u_int32_t nextNode = -1;
		bool canAdd = false;
		for(u_int32_t nodeIndex : unitigEnd._nodes){
			if(nodeIndex == nodeIndexEnd){
				canAdd = true;
				continue;
			}
			if(!canAdd) continue;
			nextNodes.push_back(nodeIndex);
			//nextNode = nodeIndex;
			//break;
		}
		
		//if(nextNode != -1) contigpath.push_back(nextNode);
		for(size_t i=0; i<n && i<nextNodes.size(); i++){
			contigpath.push_back(nextNodes[i]);
		}

	}
	*/

	void debug_checkReads(){

		cout << "Checking reads" << endl;

		KminmerParser parser(_filename_readMinimizers, _minimizerSize, _kminmerSize, true);
		auto fp = std::bind(&Assembly3::debug_checkReads_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
		parser.parse(fp);
	}

	void debug_checkReads_read(const vector<u_int64_t>& minimizers, const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex){//}, const vector<KmerVec>& kminmers_k3, const vector<ReadKminmer>& kminmersInfos_k3){
	
		vector<u_int32_t> nodePath;
		bool printRead = false;

		for(const KmerVec& vec : kminmers){
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()){
				nodePath.push_back(0);
				continue;
			}
			
			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
			nodePath.push_back(nodeName);

			if(_evaluation_hifiasmGroundTruth_nodeNamePosition.find(nodeName) == _evaluation_hifiasmGroundTruth_nodeNamePosition.end()) continue;

			if(_evaluation_hifiasmGroundTruth_nodeNamePosition[nodeName] == 4790){
				printRead = true;
			}
		}

		if(printRead || readIndex == 3427){ //|| readIndex == 2010
			cout << readIndex << endl;
			_debug_readPath[readIndex] = nodePath;
			//cout << readIndex << ": " << endl;
			//for(u_int32_t nodeName : nodePath){
			//	cout << nodeName << " ";
			//}
			//cout << endl;
		}

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

	void checkCorrectedSequence(const vector<u_int64_t> minimizers, bool print_read){
		
		//MDBG* mdbg = new MDBG(_kminmerSize);
		//mdbg->load(_inputDir + "/mdbg_nodes_init.gz");

		bool notGood = false;

		vector<KmerVec> kminmers; 
		vector<ReadKminmer> kminmersInfo;
		vector<u_int64_t> minimizersPos; 
		vector<u_int64_t> rlePositions;
		MDBG::getKminmers(_minimizerSize, _kminmerSize, minimizers, minimizersPos, kminmers, kminmersInfo, rlePositions, -1, false);

		u_int64_t nbFoundMinimizers = 0;
		for(size_t i=0; i<kminmers.size(); i++){
			KmerVec& vec = kminmers[i];
			

			if(_mdbgNoFilter->_dbg_nodes.find(vec) == _mdbgNoFilter->_dbg_nodes.end()){
				//if(i==2){ nbFailed += 1; }
				cout << "Not good: " << i << endl;
				cout << vec._kmers[0] << " " << vec._kmers[1] << " " << vec._kmers[2] << " " << vec._kmers[3] << endl;

				notGood = true;
				continue;	
			}

			//u_int32_t nodeName = mdbg->_dbg_nodes[vec]._index;

			nbFoundMinimizers += 1;

		}

		if(print_read){
			cout << "\tNb minimizers: " << kminmers.size() << endl;
			cout << "\tFound minimizers: " << nbFoundMinimizers << endl;
			//if(kminmers.size() != nbFoundMinimizers) getchar();
			if(notGood) getchar();
		}
	}




};





#endif

