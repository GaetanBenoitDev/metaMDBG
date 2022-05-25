
//contig=81 graph=4 abundance=31: 33 28 30    1 ()
//contig=81 graph=4 abundance=31: 41 35 21    3 (no composition)
//contig=81 graph=4 abundance=51: 47 32 22    3 (no composition)
//contig=81 graph=4 abundance=61: 49 29 23    2 (no composition) 
//contig=81 graph=4 abundance=31: 26 14 18    3 (with abundance correction, ab < 1 => ab = 0)
//contig=41 graph=41 abundance=31: 29 19 16    4
//contig=41 graph=41 abundance=31: 29 37 25    5 (interrupted)
//contig=10, graph=4 abundance=31: 31 52 25    9
//

///mnt/gpfs/seb/downloads/orkun_hifi_25_03_22 



//./bin/mdbgAsmMeta binPass  -o ~/workspace/run/overlap_test_multik_AD/pass_k81 -c ~/workspace/run/overlap_test_multik_AD/contigs_81.fasta.gz --bi lala.bin --bo lala2.bin --firstpass --eval -a ~/workspace/run/overlap_test_multik_AD/contigCoverages_k81.tsv --itruth ~/workspace/run/hifiasm_meta/AD2W20.asm.r_utg.fasta

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
*/

//:/mnt/chris-native/chris/Projects/AD_HiFi
//mnt/gpfs/Hackathon/meren_hifi

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
#include "toBasespace/ToBasespaceOnTheFly.hpp"
//#include "assembly/Assembly.hpp"

const u_int32_t LONG_UNITIG_LENGTH = 10000;
const u_int32_t MIN_SUPPORTED_READS = 2;
const u_int32_t MIN_REPEAT_SIDE_LENGTH = 10000;









class Assembly2 : public Tool{

public:


	string _inputDir;
	string _truthInputFilename;
	string _outputFilename;
	string _outputFilename_complete;
	bool _debug;
	string _inputFilename_unitigNt;
	string _inputFilename_unitigCluster;
	string _filename_abundance;
	bool _isFirstPass;
	string _filename_inputBinning;
	string _filename_outputBinning;
	bool _computeBinStats;

	float _minimizerDensity;
    size_t _minimizerSize;
    size_t _kminmerSize;
	size_t _kminmerSizeFirst;
    size_t _kminmerSizePrev;
    size_t _kminmerSizeLast;
	float _minimizerSpacingMean;
	float _kminmerLengthMean;
	float _kminmerOverlapMean;
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
	MDBG* _mdbgInit;
	GraphSimplify* _graph;
	ToBasespaceOnTheFly _toBasespace;
	ContigFeature _contigFeature;
	string _gfaFilename;
	ofstream _fileOutput_contigBin;
	int _nbCores;

	u_int64_t _nbHighQualityBins;
	u_int64_t _nbMedQualityBins;
	u_int64_t _nbLowQualityBins;
	u_int64_t _nbContaminatedBins;

	unordered_map<u_int32_t, u_int32_t> _contigIndex_to_componentIndex;
	unordered_map<u_int32_t, vector<u_int32_t>> _componentIndex_to_contigIndexes;
	bool _isFinalPass;
	string _binOutputDir;

	Assembly2(): Tool (){

	}

	~Assembly2(){

	}

	void execute (){

		_nbHighQualityBins = 0;
		_nbMedQualityBins = 0;
		_nbLowQualityBins = 0;
		_nbContaminatedBins = 0;

		//vector<float> lala1 = {0.0704131, 0, 0}; 
		//vector<float> lala2 = {4.79126, 6.7738, 2.99172};
	
		//ContigFeatures f1 = {0, {}, lala1, lala1};
		//ContigFeatures f2 = {0, {}, lala2, lala2};

		//cout << _contigFeature.isIntra(f1, f2) << endl;
		//exit(1);
		/*
		unordered_map<string, vector<u_int32_t>> contigToNodes;
		ifstream file_contigToNode(_inputDir + "/contigs.fasta.gz" + ".nodes.csv");
        vector<string>* fields = new vector<string>();
		string line;
		while (std::getline(file_contigToNode, line)){
			
			GfaParser::tokenize(line, fields, ';');

			u_int32_t nodeName = stoull((*fields)[0]);
			u_int32_t contigIndex = stoull((*fields)[1]);
			//cout << nodeName << " " << contigIndex << endl;
			contigToNodes["ctg" + to_string(contigIndex)].push_back(nodeName);
		}
		delete fields;
		file_contigToNode.close();

		_contigFeature.loadAbundanceFile_metabat(_filename_abundance, contigToNodes);
		*/
	

		//_contigFeature.computeFastaComposition("/home/gats/workspace/run/overlap_test_AD_k7/binByCutoff/component_16_unitigs.fasta");
		
		//exit(1);

		//unordered_map<string, vector<u_int32_t>> contigToNodesLala;
		//_contigFeature.loadAbundanceFile_metabat(_filename_abundance, contigToNodesLala);

		_outputContigFile = gzopen(_outputFilename.c_str(),"wb");
		_outputContigFile_complete = gzopen(_outputFilename_complete.c_str(),"wb");

		loadGraph();

		//generateContigs2(_inputDir + "/contigs.min.gz", _inputDir + "/contigs.fasta.gz");
		//return;




		
		//_contigFeature.loadAbundanceFile_nodename(_filename_abundance);
		
		/*
		_graph->loadState2(0, -1, _unitigDatas);
		unordered_map<string, vector<u_int32_t>> contigToNodes;
		for(const Unitig& u : _graph->_unitigs){
			if(u._index%2 == 0){
				for(u_int32_t nodeIndex : u._nodes){
					contigToNodes["ctg" + to_string(u._index)].push_back(BiGraph::nodeIndex_to_nodeName(nodeIndex));
				}
			}
		}
		_contigFeature.loadAbundanceFile_metabat(_filename_abundance, contigToNodes);
		*/


		/*
		unordered_set<u_int32_t> component;
		_graph->getConnectedComponent_unitig(141512, component);
		unordered_set<u_int32_t> contigIndexes;
		for(u_int32_t unitigIndex : component){
			
			u_int32_t contigIndex = _contigFeature.nodepathToContigIndex(_graph->_unitigs[unitigIndex]._nodes);
			contigIndexes.insert(contigIndex);
		}

		_fileOutput_contigBin = ofstream(_filename_outputBinning);
			
		u_int32_t binIndex = 0;
		cout << "Nb contigs: " << contigIndexes.size();
		for(u_int32_t contigIndex : contigIndexes){
			cout << contigIndex << " " << endl;
			if(contigIndex == -1) continue;
			_fileOutput_contigBin.write((const char*)&contigIndex, sizeof(contigIndex));
			_fileOutput_contigBin.write((const char*)&binIndex, sizeof(binIndex));
		}
		

		_fileOutput_contigBin.close();
		*/

		execute_binning2();
		mergeBins();
		//execute_detectSpecies_byCutoff();
		
		gzclose(_outputContigFile);
		gzclose(_outputContigFile_complete);

		//cout << _nbHighQualityBins << " " << _nbMedQualityBins << " " << _nbLowQualityBins << "    " << _nbContaminatedBins << endl;

	}

	void parseArgs(int argc, char* argv[]){



		cxxopts::Options options("Assembly", "");
		options.add_options()
		(ARG_OUTPUT_DIR, "", cxxopts::value<string>())
		(ARG_INPUT_FILENAME_CONTIG, "", cxxopts::value<string>()->default_value(""))
		(ARG_INPUT_FILENAME_BINNING, "", cxxopts::value<string>()->default_value(""))
		(ARG_OUTPUT_FILENAME_BINNING, "", cxxopts::value<string>()->default_value(""))
		(ARG_INPUT_FILENAME_TRUTH, "", cxxopts::value<string>()->default_value(""))
		(ARG_DEBUG, "", cxxopts::value<bool>()->default_value("false"))
		(ARG_INPUT_FILENAME_UNITIG_NT, "", cxxopts::value<string>()->default_value(""))
		(ARG_INPUT_FILENAME_UNITIG_CLUSTER, "", cxxopts::value<string>()->default_value(""))
		(ARG_FIRST_PASS, "", cxxopts::value<bool>()->default_value("false"))
		(ARG_INPUT_FILENAME_ABUNDANCE, "", cxxopts::value<string>()->default_value(""))
		(ARG_EVAL, "", cxxopts::value<bool>()->default_value("false"))
		(ARG_NB_CORES, "", cxxopts::value<int>()->default_value(NB_CORES_DEFAULT));



		if(argc <= 1){
			cout << options.help() << endl;
			exit(0);
		}

		cxxopts::ParseResult result;

		try{
			result = options.parse(argc, argv);

			_inputDir = result[ARG_OUTPUT_DIR].as<string>();
			_filename_inputContigs = result[ARG_INPUT_FILENAME_CONTIG].as<string>();
			_filename_inputBinning = result[ARG_INPUT_FILENAME_BINNING].as<string>();
			_filename_outputBinning = result[ARG_OUTPUT_FILENAME_BINNING].as<string>();
			_truthInputFilename = result[ARG_INPUT_FILENAME_TRUTH].as<string>();
			_inputFilename_unitigNt = result[ARG_INPUT_FILENAME_UNITIG_NT].as<string>();
			_inputFilename_unitigCluster = result[ARG_INPUT_FILENAME_UNITIG_CLUSTER].as<string>();
			_debug = result[ARG_DEBUG].as<bool>();
			_computeBinStats = result[ARG_EVAL].as<bool>();
			_filename_abundance = result[ARG_INPUT_FILENAME_ABUNDANCE].as<string>();
			_isFirstPass = result[ARG_FIRST_PASS].as<bool>();
			_nbCores = result[ARG_NB_CORES].as<int>();
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
		gzread(file_parameters, (char*)&_kminmerSizeFirst, sizeof(_kminmerSizeFirst));
		gzread(file_parameters, (char*)&_minimizerSpacingMean, sizeof(_minimizerSpacingMean));
		gzread(file_parameters, (char*)&_kminmerLengthMean, sizeof(_kminmerLengthMean));
		gzread(file_parameters, (char*)&_kminmerOverlapMean, sizeof(_kminmerOverlapMean));
		gzread(file_parameters, (char*)&_kminmerSizePrev, sizeof(_kminmerSizePrev));
		gzread(file_parameters, (char*)&_kminmerSizeLast, sizeof(_kminmerSizeLast));
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
		_filename_readMinimizers = _inputDir + "/read_data.txt";
		_filename_hifiasmGroundtruth = _inputDir + "/hifiasmGroundtruth.gz";
		_filename_outputContigs = _inputDir + "/contigs.min.gz";
		_filename_solidKminmers = _inputDir + "/solid.min.gz";
	}

	void loadGraph(){

		_gfaFilename = _inputDir + "/minimizer_graph.gfa";
		string gfa_filename_noUnsupportedEdges = _inputDir + "/minimizer_graph_noUnsupportedEdges.gfa";
		//string gfa_filename_unitigs = _inputDir + "/minimizer_graph_unitigs.gfa";
		string mdbg_filename = _inputDir + "/mdbg_nodes.gz";


		if(_truthInputFilename != ""){
			_mdbgInit = new MDBG(4);

			string path = _inputDir;
			while(path[path.size()-1] == '/') path.pop_back();
			fs::path p(path);
			cout << path << endl;
			cout << p.parent_path().string() + "/mdbg_nodes_init.gz" << endl;
			_mdbgInit->load(p.parent_path().string() + "/mdbg_nodes_init.gz");
			extract_truth_kminmers();
		}

		
		cout << _gfaFilename << endl;
		_mdbg = new MDBG(_kminmerSize);
		_mdbg->load(mdbg_filename);
		cout << "Nb nodes: " <<  _mdbg->_dbg_nodes.size() << endl;

		extractContigKminmers2();
		if(!_isFirstPass){
			_contigFeature.loadContigBins(_filename_inputBinning);
		}
		else{
			fs::remove(_filename_inputBinning);
		}



		//cout << "lala " << _evaluation_hifiasmGroundTruth_position.size() << endl;

		file_groundTruth = ofstream(_inputDir + "/binning_results.csv");
		file_groundTruth << "Name,Colour" << endl;

		//file_groundTruth_hifiasmContigs = ofstream(_inputDir + "/binning_results_hifiasm_" + to_string(_kminmerSize) + ".csv");
		//file_groundTruth_hifiasmContigs << "Name,Colour" << endl;

		//if(_debug){
            //gfa_filename = _inputDir + "/minimizer_graph_debug.gfa";
		//}
		
		GraphSimplify* graphSimplify = new GraphSimplify(_gfaFilename, _inputDir, 0, _kminmerSize, _nbCores, _kminmerLengthMean, _kminmerOverlapMean);
		_graph = graphSimplify;
		

		_graph->clear(0);
		_graph->compact(false, _unitigDatas);
		//_graph->removeErrors_4(_kminmerSize, _unitigDatas);

		//_graph->_contigFeature = &_contigFeature;
		//_graph->debug_writeGfaErrorfree(0, 0, -1, _kminmerSize, false, true, false, _unitigDatas, true, false, false, false, false, true, _mdbg, _minimizerSize, _nbCores, false, false);

		delete _mdbg;
		//Generate unitigs
		//cout << "Indexing reads" << endl;
		//_unitigDatas.resize(_mdbg->_dbg_nodes.size());
		//_graph->clear(0);
		//_graph->compact(false, _unitigDatas);
		//removeUnsupportedEdges(_gfaFilename, gfa_filename_noUnsupportedEdges, _graph);



		//cout << "done" << endl;
	


		//cout << graphSimplify->_graphSuccessors->_nbEdges << endl;
		//graphSimplify->execute(5, _kminmerSize);
		//graphSimplify->debug_writeGfaErrorfree(1000, PathExplorer::computeAbundanceCutoff(1000, 0, CutoffType::ERROR), -1, _kminmerSize, false, true, false, _unitigDatas);




		//_toBasespace.create(_inputDir);


		//_graph->loadState2(0, -1, _unitigDatas);
		//generateFasta(_inputDir + "/allContigs.fasta.gz");
		//generateFastaExpanded(_inputDir + "/allContigsExpanded.fasta.gz");

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
		infile.close();
	}

	void writeScgDebug(const string& contigFilename, unordered_set<u_int32_t>& component){

		unordered_set<u_int32_t> scgs;
		scgs.insert(0);
		/*
		for(const auto& it : _scgClusters){
			const vector<ScgCluster>& scgClusters = it.second;
			for(const ScgCluster& scgCluster : scgClusters){
				scgs.insert(scgCluster._scgIndex);
			}
		}
		*/
				
		ofstream file(contigFilename + ".scg.csv");
		string header = "Name";

		for(size_t i=0; i<scgs.size(); i++){
			header += ",Color";//",SCG_" + to_string(i);
		}
		file << header << endl;

		unordered_set<u_int32_t> writtenUnitigs;
				
		for(u_int32_t unitigIndex : component){
					
			const Unitig& u = _graph->_unitigs[unitigIndex];
			u_int32_t contigIndex = u._index;
			//if(u._length < 10000) continue;
			
			if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._startNode)) != writtenUnitigs.end()) continue;
			if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._endNode)) != writtenUnitigs.end()) continue;

			writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._startNode));
			writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._endNode));

			if(_scgClusters.find(unitigIndex) == _scgClusters.end()){
				for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
					u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
					string line = "";
					line += to_string(nodeName);
					for(size_t i=0; i<scgs.size(); i++){
						line += ",-1";
					}
					file << line << endl;
				}

			}

			unordered_map<u_int32_t, u_int32_t> unitigScgs;
			for(const ScgCluster& scgCluster : _scgClusters[unitigIndex]){
				unitigScgs[scgCluster._scgIndex] = scgCluster._scgCluster;
			}

			for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
				string line = "";
				line += to_string(nodeName);
				for(size_t i=0; i<scgs.size(); i++){
					if(unitigScgs.find(i) == unitigScgs.end()){
						line += ",-1";
					}
					else{
						line += "," + to_string(unitigScgs[i]);
					}
				}
				line += "\n";
				file << line << endl;
			}

		}

		file.close();

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



	void indexReads_read(const vector<u_int64_t>& minimizers, const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex){//}, const vector<KmerVec>& kminmers_k3, const vector<ReadKminmer>& kminmersInfos_k3){

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

		
		KminmerParser parser(_filename_readMinimizers, _minimizerSize, _kminmerSize, true);
		//auto fp = std::bind(&Assembly::indexReads_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
		auto fp = std::bind(&Assembly2::indexReads_read, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
		parser.parse(fp);
		
		if(_kminmerSize > 3) return;

		vector<DbgEdge> removedEdges;

        for(const Unitig& unitig: _graph->_unitigs){

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

	void extract_truth_kminmers(){

		_file_groundTruth_hifiasm_position.open(_inputDir + "/groundtruth_hifiasm_position.csv");
		_file_groundTruth_hifiasm_position << "Name,Position" << endl;

		_extract_truth_kminmers_read_position = 0;
		_minimizerParser = new MinimizerParser(_minimizerSize, _minimizerDensity);
		
		auto fp = std::bind(&Assembly2::extract_truth_kminmers_read, this, std::placeholders::_1);
		ReadParser readParser(_truthInputFilename, true, false);
		readParser.parse(fp);

		_file_groundTruth_hifiasm_position.close();

		//delete _minimizerParser;
	}

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
		MDBG::getKminmers(_minimizerSize, 4, minimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, readIndex, false);

		for(size_t i=0; i<kminmers.size(); i++){

			KmerVec& vec = kminmers[i];
			//if(_mdbg->_dbg_nodes.find(kminmers[i]) == _mdbg->_dbg_nodes.end()) continue;

			//u_int32_t nodeName = _mdbg->_dbg_nodes[kminmers[i]]._index;
			if(_mdbgInit->_dbg_nodes.find(vec) != _mdbgInit->_dbg_nodes.end()){

				u_int32_t nodeName = _mdbgInit->_dbg_nodes[vec]._index;
				//2262408
				//if("utg009434l" == string(read->name.s)){
				//	cout << nodeName << endl;
				//	_hifiasm_startingNodenames.insert(nodeName);
				//}

				_evaluation_hifiasmGroundTruth_nodeName_to_unitigName[_mdbgInit->_dbg_nodes[vec]._index].push_back(read._header);
				_evaluation_hifiasmGroundTruth_path.push_back(_mdbgInit->_dbg_nodes[vec]._index);

				if(_evaluation_hifiasmGroundTruth_nodeNamePosition.find(_mdbgInit->_dbg_nodes[vec]._index) == _evaluation_hifiasmGroundTruth_nodeNamePosition.end()){
					_evaluation_hifiasmGroundTruth_nodeNamePosition[_mdbgInit->_dbg_nodes[vec]._index] = _extract_truth_kminmers_read_position;
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








	int computeBinStats(const string& outputDir, const vector<string>& bin, const string& filename_binStats){

		int lala = 0;
		const string& basespaceUnitigFilename = outputDir + "/contigsTmp.fasta";
		ofstream basespaceUnitigFile = ofstream(basespaceUnitigFilename);

		for(const string& sequence : bin){



			string header = ">ctg" + to_string(lala);
			basespaceUnitigFile << header << endl;
			basespaceUnitigFile << sequence << endl;
			
			lala += 1;
		}

		basespaceUnitigFile.close();

		string command_annotation = "python3 /home/gats/workspace/scripts/annotation/annotation2/computeBinStats.py -t 4 " + basespaceUnitigFilename + " " + filename_binStats;
		cout << command_annotation << endl;
		int ret = system(command_annotation.c_str());
		if(ret != 0){
			cerr << "Command failed: " << ret << endl;
			//exit(ret);
		}

		return ret;
		
	}

	void binByReadpath(u_int32_t source_nodeIndex, unordered_set<u_int32_t>& processedNodeNames, unordered_set<u_int32_t>& processedContigIndex, const string& clusterDir, const string& filename_binStats, ofstream& fileHifiasmAll, ofstream& fileComponentNodeAll, u_int64_t& clusterIndex, u_int32_t& binIndex, u_int64_t lengthThreshold, u_int32_t& binIndexValid){


		u_int32_t source_unitigIndex = _graph->nodeIndex_to_unitigIndex(source_nodeIndex);
		const Unitig& unitig_model = _graph->_unitigs[source_unitigIndex];

		u_int32_t contigIndex_model = _contigFeature.nodepathToContigIndex(unitig_model._nodes);
		if(contigIndex_model == -1){
			//cout << "\tNo contig index" << endl;
			return;
		}

		//string unitigSequence_model;
		//_toBasespace.createSequence(unitig_model._nodes, unitigSequence_model);
		//_contigFeature.nodepathToContigSequence(unitig_model._nodes, unitigSequence_model, contigIndex_model);
		//if(unitigSequence_model.empty()) return;	


		if(processedContigIndex.find(contigIndex_model) != processedContigIndex.end()){
			//cout << "\tAlready binned" << endl;
			return;
		}



		//unordered_map<u_int32_t, u_int32_t> lala;
		unordered_set<u_int32_t> componentContigIndex;

		//vector<string> bin;

		unordered_set<u_int32_t> binnedContigIndex;

      	unordered_set<u_int32_t> isVisited;
        queue<u_int32_t> queue;

        queue.push(source_unitigIndex);

		//unordered_set<u_int32_t> binContigIndexes_set;
		vector<u_int32_t> binContigIndexes;
		u_int32_t currentBinIndex = _contigFeature.contigIndexToBinIndex(contigIndex_model);

		//16300

		//cout << "Contig index: " << contigIndex_model << " " << currentBinIndex << endl;
		//cout << _contigFeature._contigIndex_to_binIndex.size() << endl;
		//cout << _contigFeature._contigIndex_to_binIndex[contigIndex_model] << endl;

		if(currentBinIndex == -1){
			//cout << "No bin" << endl;
			if(!_isFirstPass) return;
			binContigIndexes.push_back(contigIndex_model);
		}
		else{
			//cout << "Has bin: " << _contigFeature._binIndex_to_contigIndex[currentBinIndex].size() << endl;
			for(u_int32_t contigIndex : _contigFeature._binIndex_to_contigIndex[currentBinIndex]){
				binContigIndexes.push_back(contigIndex);
			}
		}

		//cout << "Current: " << contigIndex_model << endl;
		for(u_int32_t contigIndex : binContigIndexes){
			assignContigToBin(contigIndex, binIndex);
			processedContigIndex.insert(contigIndex);
			componentContigIndex.insert(contigIndex);
			//binContigIndexes_set.insert(contigIndex);
		}

		//cout << "lala" << " " << binContigIndexes.size() << endl;
		//getchar();
		//bin.push_back(unitigSequence_model);
		//cout << "\tContig index model: " << contigIndex_model << endl;

		//vector<float> compositionModel;
		//_contigFeature.nodepathToComposition(unitig_model._nodes, compositionModel);

		//vector<float> abundancesModel;
		//vector<float> abundancesModel_var;
		//bool isAbValid = _contigFeature.sequenceToAbundance(unitig_model._nodes, abundancesModel, abundancesModel_var);
		//if(!isAbValid) return;
			
		//ContigFeatures contigFeatureModel = {contigIndex_model, compositionModel, abundancesModel, abundancesModel_var};

		unordered_set<u_int32_t> nonIntraUnitigs;

		unordered_set<u_int32_t> writtenUnitigs;
		//u_int32_t unitigIndex_model = _graph->nodeIndex_to_unitigIndex(unitig._startNode);
		//const Unitig& unitig_model = _graph->_unitigs[unitigIndex_model];
		writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(unitig_model._startNode));
		writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(unitig_model._endNode));

		//unordered_set<u_int32_t> allComponentNodenames;




		//for(u_int32_t nodeIndex : unitig_model._nodes){
		//	u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
		//	allComponentNodenames.insert(nodeName);
		//}
		/*
		unordered_set<u_int32_t> componentDebug;

		while(!queue.empty()){

            u_int64_t unitigIndex = queue.front();
            queue.pop();

            if (isVisited.find(unitigIndex) != isVisited.end()) continue;

            isVisited.insert(unitigIndex);
            isVisited.insert(_graph->unitigIndex_toReverseDirection(unitigIndex));





			size_t componentNbNodes = 9999999;
			size_t componentLength = 100000;

			unordered_set<u_int32_t> component;
			unordered_set<u_int32_t> componentNbNodes2;
			//_graph->getConnectedComponent_readpath(unitigIndex, _unitigDatas, 1, component);
			//_graph->getConnectedComponent_unitig_nt(unitigIndex, componentLength, component);
			_graph->getConnectedComponent_unitig(unitigIndex, componentNbNodes, component);
			//cout << "\tComponent size: " << componentNbNodes2.size() << " " << component.size() << endl;
			componentDebug = component;

			//component = componentNbNodes2;





			
			//------------------------------------------------------------------- short reads
			ofstream file_coverages(_inputDir + "/coverages.csv");
        	file_coverages << "Name,Abundance" << endl;

			ofstream file_coverages_component(_inputDir + "/coverages_component.csv");
        	file_coverages_component << "Name,Color" << endl;

			_graph->_isNodeInvalid_tmp.clear();

			
			unordered_set<u_int32_t> unitigIndexCovLala;
			for(u_int32_t unitigIndex : component){
				
				const Unitig& u = _graph->_unitigs[unitigIndex];
				
				if(unitigIndexCovLala.find(BiGraph::nodeIndex_to_nodeName(u._startNode)) != unitigIndexCovLala.end()) continue;
				if(unitigIndexCovLala.find(BiGraph::nodeIndex_to_nodeName(u._endNode)) != unitigIndexCovLala.end()) continue;

				unitigIndexCovLala.insert(BiGraph::nodeIndex_to_nodeName(u._startNode));
				unitigIndexCovLala.insert(BiGraph::nodeIndex_to_nodeName(u._endNode));




				float abRatio = 1;
				vector<float> abundances;
				vector<float> abundances_var;
				bool isAbValid = _contigFeature.sequenceToAbundance(u._nodes, abundances, abundances_var);
				if(isAbValid){


					vector<float> means_f1 = abundancesModel;
					vector<float> means_f2 = abundances;

					float mean_ratio = 0;
					for(size_t i=0; i<means_f1.size(); i++){
						if(means_f2[i] > 0){
							float ratio = means_f1[i] / means_f2[i];
							mean_ratio += ratio;
							//cout << ratio << " ";
						}
						else{
							//cout << "0" << " ";
						}
					}
					//cout << endl;

					mean_ratio /= means_f1.size();

					for(size_t i=0; i<means_f1.size(); i++){
						means_f2[i] *= mean_ratio;
					}


					bool hasValidAb = true;

					for(size_t i=0; i<means_f2.size(); i++){
						if(abundancesModel[i] <= 2) continue;
						float abOffset = abundancesModel[i] * abRatio;
						if(means_f2[i] < abundancesModel[i]-abOffset || means_f2[i] > abundancesModel[i]+abOffset){
							hasValidAb = false;
							break;
						}
					}

					if(hasValidAb){
						for(u_int32_t nodeIndex : u._nodes){
							file_coverages_component << BiGraph::nodeIndex_to_nodeName(nodeIndex) << ",green" << endl;
						}
					}
					else{
						for(u_int32_t nodeIndex : u._nodes){
							file_coverages_component << BiGraph::nodeIndex_to_nodeName(nodeIndex) << ",red" << endl;
							_graph->_isNodeInvalid_tmp.insert(_graph->nodeIndex_to_unitigIndex(nodeIndex));
							_graph->_isNodeInvalid_tmp.insert(_graph->nodeIndex_to_unitigIndex(_graph->nodeIndex_toReverseDirection(nodeIndex)));
						}

					}

				}


				unordered_set<u_int32_t> contigCovs;
				for(u_int32_t nodeIndex : u._nodes){
					u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);

					if(_contigFeature._nodeName_to_contigIndex.find(nodeName) == _contigFeature._nodeName_to_contigIndex.end()) continue;
					
					u_int32_t contigIndex = _contigFeature._nodeName_to_contigIndex[nodeName];

					if(contigCovs.find(contigIndex) != contigCovs.end()) continue;

					const vector<float>& contigAbundances = _contigFeature._contigCoverages[contigIndex];

					contigCovs.insert(contigIndex);

					string str = "";
					for(u_int32_t ab : contigAbundances){
						str += to_string(ab) + ";";
					}

					file_coverages << nodeName << "," << str << endl;
				}

				//vector<float> abundances;
				//vector<float> abundancesVar;
				//_contigFeature.sequenceToAbundance(u._nodes, abundances, abundancesVar);



			}

			file_coverages.close();
			file_coverages_component.close();
			
			component.clear();
        	_graph->getConnectedComponent_unitig_nt(unitigIndex, componentLength, component);
			//---------------------------- Short reads End
			
		
			//cout << "Component size filtered: " << component.size() << endl;

			
			//cout << "\t";
			//for(u_int32_t count : abundancesModel) cout << count << " ";
			//cout << endl;
			//getchar();
			

			//for(u_int32_t utgIndex : component) queue.push(utgIndex);




			//"Reprise: better shortest path from set of nodes, si on entre dans un noeud qui a une distance plus faible que celle actuelle, pas besoin de le visiter"
			//"Extraire les bins (contigs) et les evaluer avec checkm"
			//"evaluer metaflye sur Human data"
			//if(unitigIndex != source_unitigIndex && unitigIndex != _graph->unitigIndex_toReverseDirection(source_unitigIndex)){
			vector<u_int32_t> validUnitigs;
			unordered_set<u_int32_t> existingContigIndexes;
			
			for(u_int32_t unitigIndex : component){
				const Unitig& u = _graph->_unitigs[unitigIndex];

				//if(nonIntraUnitigs.find(unitigIndex) != nonIntraUnitigs.end()) continue;



				if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._startNode)) != writtenUnitigs.end()) continue;
				if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._endNode)) != writtenUnitigs.end()) continue;

				writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._startNode));
				writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._endNode));

				//if(isContigBinned(u._nodes, processedNodeNames)) continue;
				//if(u._length < 2500) continue;
				//if(u._nbNodes < _kminmerSize*2) continue;
				//if(u._abundance < _minUnitigAbundance) continue;
				//if(u._nbNodes <= _kminmerSize*2) continue;
				
				

				


				//string unitigSequence;
				u_int32_t contigIndex = _contigFeature.nodepathToContigIndex(u._nodes);
				if(contigIndex == -1) continue;

				//_toBasespace.createSequence(u._nodes, unitigSequence);
				//_contigFeature.nodepathToContigSequence(u._nodes, unitigSequence, contigIndex);
				//if(unitigSequence.empty()) continue;

				if(componentContigIndex.find(contigIndex) != componentContigIndex.end()){
					continue;
				}

				if(processedContigIndex.find(contigIndex) != processedContigIndex.end()){
					//validUnitigs.push_back(unitigIndex);
					continue;
				}

				existingContigIndexes.insert(contigIndex);
			}
			*/





			//cout << endl << "\tUnitig: " << BiGraph::nodeIndex_to_nodeName(u._startNode) << " " << u._length << " " << u._nodes.size() << " " << unitigSequence.size() << endl;
			
			/*
			vector<float> composition;
			bool hasComposition = _contigFeature.nodepathToComposition(u._nodes, composition);
			//_contigFeature.sequenceToComposition(unitigSequence, composition);
			
			vector<float> abundances;
			vector<float> abundancesVar;
			bool hasAbundances = _contigFeature.sequenceToAbundance(u._nodes, abundances, abundancesVar);

			ContigFeatures contigFeature = {contigIndex, composition, abundances, abundancesVar};
			*/

			u_int32_t componentIndex = _contigIndex_to_componentIndex[contigIndex_model];
			const vector<u_int32_t> existingContigIndexes = _componentIndex_to_contigIndexes[componentIndex];

			for(u_int32_t contigIndex : existingContigIndexes){

				vector<u_int32_t> contigIndexes;
				u_int32_t newBinIndex = _contigFeature.contigIndexToBinIndex(contigIndex);
				if(newBinIndex == -1){
					if(_contigFeature._contigLengths[contigIndex] >= lengthThreshold){
						if(componentContigIndex.find(contigIndex) == componentContigIndex.end()){
							contigIndexes.push_back(contigIndex);
						}
					}
				}
				else{
					for(u_int32_t contigIndex : _contigFeature._binIndex_to_contigIndex[newBinIndex]){
						if(_contigFeature._contigLengths[contigIndex] < lengthThreshold) continue;
						//if(binContigIndexes_set.find(contigIndex) != binContigIndexes_set.end()) continue; //already in current bin
						if(componentContigIndex.find(contigIndex) != componentContigIndex.end()) continue;
						contigIndexes.push_back(contigIndex);
					}
				}

				if(contigIndexes.size() == 0) continue;

				//if(_contigFeature.isIntra(contigFeatureModel, contigFeature, hasComposition, hasAbundances)){
				if(_contigFeature.isIntra(binContigIndexes, contigIndexes)){

					/*
					cout << "\tComposition dist: " << -log10(_contigFeature.computeCompositionProbability(compositionModel, composition)) << " " << _contigFeature.cal_tnf_dist(compositionModel, composition, contigIndex_model, contigIndex) << endl; // << " " << _contigFeature.computeEuclideanDistance(abundancesModel, abundances_init) << " " << _contigFeature.computeProbability(contigFeatureModel, contigFeature_init)  << endl; 
					//cout << "\tComposition dist (extended): " << _contigFeature.computeEuclideanDistance(compositionModel_init, composition) << " " << _contigFeature.computeEuclideanDistance(abundancesModel, abundances) << " " << _contigFeature.computeProbability(contigFeatureModel, contigFeature)  << endl; 
					int nnz = 0;
					cout << "\tMetabat Abudance: " << _contigFeature.cal_abd_dist(contigFeatureModel, contigFeature, nnz) << endl;// << " " << _contigFeature.cal_abd_dist(contigFeatureModel, contigFeature, nnz) << endl;
					cout << "\tMetabat Abudance new: " << _contigFeature.cal_abd_dist_new(contigFeatureModel, contigFeature, nnz) << endl; // << " " << _contigFeature.cal_abd_dist_new(contigFeatureModel, contigFeature, nnz) << endl;
					cout << "\tMetacoag abundance: " << -log10(_contigFeature.computeAbundanceProbability_new(abundancesModel, abundances)) << endl; // << " " << -log10(_contigFeature.computeAbundanceProbability_new(abundancesModel, abundances)) << endl;
					cout << "\tCorrelation: " << _contigFeature.computeAbundanceCorrelation(abundancesModel, abundances) << endl; // << " " << _contigFeature.computeAbundanceCorrelation(abundancesModel, abundances) << endl;
					//cout << "\tCorrelation: " << _contigFeature.computeAbundanceCorrelation_new(abundancesModel, abundances_init) << " " << _contigFeature.computeAbundanceCorrelation_new(abundancesModel, abundances) << endl;
					cout << "\t";
					for(u_int32_t count : abundancesModel) cout << count << " ";
					cout << endl;
					cout << "\t";
					for(u_int32_t count : abundances) cout << count << " ";
					cout << endl;
					cout << "\t>" << u._length << ": " << contigIndex << endl;
					*/
				
					//lala[BiGraph::nodeIndex_to_nodeName(u._startNode)] += 1;
					//lala[BiGraph::nodeIndex_to_nodeName(u._endNode)] += 1;

					//if(lala[BiGraph::nodeIndex_to_nodeName(u._startNode)] > 1 || lala[BiGraph::nodeIndex_to_nodeName(u._endNode)] > 1){
					//	cout << "nani??" << endl;
					//	getchar();
					//}


					//binnedUnitigs.insert(u._index);

					for(u_int32_t contigIndex : contigIndexes){
						assignContigToBin(contigIndex, binIndex);
						binContigIndexes.push_back(contigIndex);
						processedContigIndex.insert(contigIndex);
						componentContigIndex.insert(contigIndex);

						const vector<u_int32_t>& contigNodeNames = _contigFeature._contigIndex_to_nodeName[contigIndex];
						
						//add first contig node
						for(long i=0; i<contigNodeNames.size(); i++){
							u_int32_t nodeName = contigNodeNames[i];
							u_int32_t nodeIndex = BiGraph::nodeName_to_nodeIndex(nodeName, true);
							if(_graph->_isNodeValid2.find(nodeIndex) == _graph->_isNodeValid2.end()) continue;

							//validUnitigs.push_back(_graph->nodeIndex_to_unitigIndex(nodeIndex));
							break;
						}

						//add last contig node
						for(long i=contigNodeNames.size()-1; i>=0; i--){
							u_int32_t nodeName = contigNodeNames[i];
							u_int32_t nodeIndex = BiGraph::nodeName_to_nodeIndex(nodeName, true);
							if(_graph->_isNodeValid2.find(nodeIndex) == _graph->_isNodeValid2.end()) continue;

							//validUnitigs.push_back(_graph->nodeIndex_to_unitigIndex(nodeIndex));
							break;
						}

					}
					

					//for(u_int32_t nodeIndex : u._nodes){
					//	u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
					//	processedNodeNames.insert(nodeName);
					//	allComponentNodenames.insert(nodeName);
					//}


					
					//string unitigSequence_init;
					//_toBasespace.createSequence(u._nodes, unitigSequence_init);
					//bin.push_back(unitigSequence);

					//queue.push(unitigIndex);




				}
				else{
					//nonIntraUnitigs.insert(u._index);
					//nonIntraUnitigs.insert(_graph->unitigIndex_toReverseDirection(u._index));
				}
				
			}

			//for(u_int32_t unitigIndex : validUnitigs){
				//queue.push(unitigIndex);
			//}


		//}

		/*
		vector<string> bin;
		for(u_int32_t contigIndex : binContigIndexes){
			const string& sequence = _contigFeature._contigSequences[contigIndex];
			bin.push_back(sequence);
		}

		u_int64_t lengthTotal = 0;
		for(const string& contig : bin){
			lengthTotal += contig.size();
		}

		//cout << _contigFeature._binningThreshold << " " << (_contigFeature._binningThreshold == 0.65f) << endl;

		if(_computeBinStats ){
			if(lengthTotal > 20000 && _isFinalPass){//_contigFeature._binningThreshold == 0.65f && lengthThreshold == 10000){

				//cout << binIndexValid << " " << binContigIndexes.size() << " " << lengthTotal << endl;
				dumpBin(binIndexValid, binContigIndexes);
				
				binIndexValid += 1;
				//cout << "bin done" << endl;

				
				cout << _nbHighQualityBins << " " << _nbMedQualityBins << " " << _nbLowQualityBins << "    " << _nbContaminatedBins << endl;

				int ret = computeBinStats(clusterDir, bin, filename_binStats);

				if(ret == 0){

					string lastLine = "";
					string line = "";
					ifstream file_binStats(filename_binStats);
					while (std::getline(file_binStats, line)){
						if(line.empty()) continue;
						lastLine = line;
					}
					file_binStats.close();
					
					vector<string>* fields = new vector<string>();
					GfaParser::tokenize(lastLine, fields, ',');
					float completeness = stof((*fields)[0]);
					float contamination = stof((*fields)[1]);
					delete fields;
					cout << "Result: " << completeness << " " << contamination << endl;

					if(completeness > 0.8 && contamination < 0.05){


						if(!_truthInputFilename.empty()){
						
							unordered_set<string> hifiasmUnitigNames;
							
							for(u_int32_t contigIndex : binContigIndexes){


								const string& sequence = _contigFeature._contigSequences[contigIndex];

								string rleSequence;
								vector<u_int64_t> rlePositions;
								_encoderRLE.execute(sequence.c_str(), sequence.size(), rleSequence, rlePositions);

								vector<u_int64_t> minimizers;
								vector<u_int64_t> minimizers_pos;
								_minimizerParser->parse(rleSequence, minimizers, minimizers_pos);

								vector<KmerVec> kminmers; 
								vector<ReadKminmer> kminmersInfo;
								MDBG::getKminmers(_minimizerSize, 4, minimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, 0, false);

								for(size_t i=0; i<kminmers.size(); i++){

									KmerVec& vec = kminmers[i];

									//cout << (_mdbgInit->_dbg_nodes.find(vec) != _mdbgInit->_dbg_nodes.end()) << endl;
									//if(_mdbg->_dbg_nodes.find(kminmers[i]) == _mdbg->_dbg_nodes.end()) continue;

									//u_int32_t nodeName = _mdbg->_dbg_nodes[kminmers[i]]._index;
									if(_mdbgInit->_dbg_nodes.find(vec) != _mdbgInit->_dbg_nodes.end()){

										u_int32_t nodeName = _mdbgInit->_dbg_nodes[vec]._index;
										if(_evaluation_hifiasmGroundTruth_nodeName_to_unitigName.find(nodeName) != _evaluation_hifiasmGroundTruth_nodeName_to_unitigName.end()){
											for(string& unitigName : _evaluation_hifiasmGroundTruth_nodeName_to_unitigName[nodeName]){
												if(hifiasmUnitigNames.find(unitigName) != hifiasmUnitigNames.end()) continue;
												hifiasmUnitigNames.insert(unitigName);
												fileHifiasmAll << unitigName << "," << clusterIndex << endl;
											}
										}
									}
								}
							}

								
						}

						//for(u_int32_t nodeName : allComponentNodenames){
						//	fileComponentNodeAll << nodeName << "," << clusterIndex << endl;
						//}
						//fileComponentNodeAll.flush();

						fileHifiasmAll.flush();

						//if(clusterIndex >= 4) getchar();
						clusterIndex += 1;
					}

					if(contamination > 0.05){

						
						//unordered_set<u_int32_t> validNodes;
						//for (u_int32_t unitigIndex : componentDebug){
						//	for(u_int32_t nodeIndex : _graph->_unitigs[unitigIndex]._nodes){
						//		u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
						//		//if(_debug_groundTruthNodeNames.find(nodeName) == _debug_groundTruthNodeNames.end()) continue;
						//		validNodes.insert(nodeName);
						//	}
						//}
						//string outputFilename = _inputDir + "/minimizer_graph_binning.gfa";
						//GfaParser::rewriteGfa_withoutNodes(_gfaFilename, outputFilename, validNodes, _graph->_isEdgeRemoved, _graph->_graphSuccessors);
			

						

						//for(u_int32_t contigIndex : binContigIndexes){
						//	const vector<float>& coverages = _contigFeature._contigCoverages[contigIndex];
						//	for(float ab : coverages){
						//		cout << ab << " ";
						//	}
						//	cout << endl;
						//}

						for(size_t i=0; i<binContigIndexes.size(); i++){
							for(size_t j=i+1; j<binContigIndexes.size(); j++){
								for(float ab : _contigFeature._contigCoverages[binContigIndexes[i]]){
									cout << ab << " ";
								}
								cout << endl;
								for(float ab : _contigFeature._contigCoverages[binContigIndexes[j]]){
									cout << ab << " ";
								}
								cout << endl;

								int nnz = 0;

								vector<float> means_f1 = _contigFeature._contigCoverages[binContigIndexes[i]];
								vector<float> means_f2 = _contigFeature._contigCoverages[binContigIndexes[j]];
								//cout << endl;
								float mean_ratio = 0;
								for(size_t i=0; i<means_f1.size(); i++){
									if (means_f1[i] > 0 || means_f2[i] > 0) {
										if(means_f2[i] > 0){
											float ratio = means_f1[i] / means_f2[i];
											mean_ratio += ratio;
											//cout << ratio << " ";
										}
										nnz += 1;
									}
									else{
										//cout << "0" << " ";
									}
								}
								//cout << endl;
								
								mean_ratio /= nnz;
								//cout << mean_ratio << endl;
								
								for(size_t i=0; i<means_f1.size(); i++){
									means_f2[i] *= mean_ratio;
									cout << means_f2[i] << " ";
								}
								cout << endl;


								cout << _contigFeature.computeDistance(binContigIndexes[i], binContigIndexes[j]) << endl;
							}
						}

						//getchar();
					}

					float qualityScore = completeness - 5*contamination;
				
					if(contamination > 0.05){
						_nbContaminatedBins += 1;
					}
					else if (completeness >= 0.9 and contamination <= 0.05)
						_nbHighQualityBins += 1;
					else if (completeness >= 0.7 and contamination <= 0.1)
						_nbMedQualityBins += 1;
					else if (qualityScore >= 0.5)
						_nbLowQualityBins += 1;



					cout << _nbHighQualityBins << " " << _nbMedQualityBins << " " << _nbLowQualityBins << "    " << _nbContaminatedBins << endl;

				}
				else{
					exit(1);
				}
				
			}
		}
		*/
		
		binIndex += 1;
		//getchar();
	}

	/*
	bool dumpBin(const u_int64_t binIndex, const vector<u_int32_t>& binContigIndexes){


		u_int64_t lengthTotal = 0;
		for(u_int32_t contigIndex : binContigIndexes){
			const string& sequence = _contigFeature._contigSequences[contigIndex];
			lengthTotal += sequence.size();
		}

		if(lengthTotal < 20000) return false;

		const string& filename = _binOutputDir + "/bin_" + to_string(binIndex) + ".fasta";
		ofstream file = ofstream(filename);

		for(u_int32_t contigIndex : binContigIndexes){
			const string& sequence = _contigFeature._contigSequences[contigIndex];

			string header = ">ctg" + to_string(contigIndex);
			file << header << endl;
			file << sequence << endl;
			
		}

		file.close();

		return true;
	}
	*/


	
	unordered_map<u_int32_t, u_int32_t> _contigIndex_to_binIndex;
	void assignContigToBin(u_int32_t contigIndex, u_int32_t binIndex){

		_contigIndex_to_binIndex[contigIndex] = binIndex;
		//if(_contigFeature._contigLengths[contigIndex] > 500000)
		//cout << "Assign: " << contigIndex << " -> " << binIndex << endl;
		//_fileOutput_contigBin.write((const char*)&contigIndex, sizeof(contigIndex));
		//_fileOutput_contigBin.write((const char*)&binIndex, sizeof(binIndex));
		//_fileOutput_contigBin.flush();

		//if(contigIndex == 209) getchar();
	}
	

	void dumpBins(){
		
		_fileOutput_contigBin = ofstream(_filename_outputBinning);

		for(const auto& it : _contigIndex_to_binIndex){
			u_int32_t contigIndex = it.first;
			u_int32_t binIndex = it.second;
			_fileOutput_contigBin.write((const char*)&contigIndex, sizeof(contigIndex));
			_fileOutput_contigBin.write((const char*)&binIndex, sizeof(binIndex));
		}

		_fileOutput_contigBin.close();
	}


	ofstream _fileTestLala;


	ofstream _file_contigToNode;

	void extractContigKminmers2 (){

		_file_contigToNode = ofstream(_inputDir + "/nodeToContig.csv");
		_file_contigToNode << "Name,Color" << endl;
		_contigFeature.loadAbundanceFile_metabat(_filename_abundance);

		ReadParser parser(_filename_inputContigs, true, _minimizerSize, _kminmerSize, _minimizerDensity);
		auto fp = std::bind(&Assembly2::extractContigKminmers_read2, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6);
		parser.parseKminmers(fp);

		_file_contigToNode.close();

		//getchar();
		//vector<u_int32_t> nodePath = {933376, 1651014, 1762772, 732878};
		//vector<float> lala1;
		//vector<float> lala2;
		//_contigFeature.sequenceToAbundance(nodePath, lala1, lala2);

		//exit(1);

		/*
		if(!_truthInputFilename.empty()){
		
			ofstream file_contigToHifiasmUnitig(_inputDir + "/hifiasm_contig.csv");
			file_contigToHifiasmUnitig << "Name,Color" << endl;

			unordered_set<string> hifiasmUnitigNames;
			
			for(auto& it : _contigFeature._contigSequences){

				u_int32_t contigIndex = it.first;
				const string& sequence = it.second;

				string rleSequence;
				vector<u_int64_t> rlePositions;
				_encoderRLE.execute(sequence.c_str(), sequence.size(), rleSequence, rlePositions);

				vector<u_int64_t> minimizers;
				vector<u_int64_t> minimizers_pos;
				_minimizerParser->parse(rleSequence, minimizers, minimizers_pos);

				vector<KmerVec> kminmers; 
				vector<ReadKminmer> kminmersInfo;
				MDBG::getKminmers(_minimizerSize, 4, minimizers, minimizers_pos, kminmers, kminmersInfo, rlePositions, 0, false);

				for(size_t i=0; i<kminmers.size(); i++){

					KmerVec& vec = kminmers[i];

					//cout << (_mdbgInit->_dbg_nodes.find(vec) != _mdbgInit->_dbg_nodes.end()) << endl;
					//if(_mdbg->_dbg_nodes.find(kminmers[i]) == _mdbg->_dbg_nodes.end()) continue;

					//u_int32_t nodeName = _mdbg->_dbg_nodes[kminmers[i]]._index;
					if(_mdbgInit->_dbg_nodes.find(vec) != _mdbgInit->_dbg_nodes.end()){

						u_int32_t nodeName = _mdbgInit->_dbg_nodes[vec]._index;
						if(_evaluation_hifiasmGroundTruth_nodeName_to_unitigName.find(nodeName) != _evaluation_hifiasmGroundTruth_nodeName_to_unitigName.end()){
							for(string& unitigName : _evaluation_hifiasmGroundTruth_nodeName_to_unitigName[nodeName]){
								if(hifiasmUnitigNames.find(unitigName) != hifiasmUnitigNames.end()) continue;
								hifiasmUnitigNames.insert(unitigName);
								file_contigToHifiasmUnitig << unitigName << "," << "red" << endl;
							}
						}
					}
				}
			}

			file_contigToHifiasmUnitig.close();
				
		}
		*/

		/*
		cout << _contigFeature._contigSequences[31].size() << endl;
		cout << _contigFeature._contigSequences[48].size() << endl;
		cout << _contigFeature._contigSequences[46].size() << endl;
		cout << _contigFeature._contigSequences[39].size() << endl;
		cout << _contigFeature._contigSequences[30].size() << endl;
		cout << _contigFeature._contigSequences[33].size() << endl;

		cout << _contigFeature.computeDistance(31, 48) << endl;
		cout << _contigFeature.computeDistance(31, 46) << endl;
		cout << _contigFeature.computeDistance(31, 39) << endl;
		cout << _contigFeature.computeDistance(31, 30) << endl;
		cout << _contigFeature.computeDistance(31, 33) << endl;
		getchar();
		*/
	}


	void extractContigKminmers_read2(const vector<KmerVec>& kminmers, const vector<ReadKminmer>& kminmersInfos, u_int64_t readIndex, u_int64_t datasetIndex, const string& header, const string& seq){

		//cout << readIndex << " " << kminmers.size() << endl;
		if(seq.size() < 2500) return;
		//cout << seq.size() << endl;
		//if(readIndex == 20010) getchar(); 

		//cout << readIndex << " " << seq.size() << endl;

		vector<float> composition;
		_contigFeature.sequenceToComposition(seq, composition);
		_contigFeature._contigCompositions[readIndex] = composition;

		_contigFeature._contigLengths[readIndex] = seq.size();
		//unordered_set<u_int32_t> nodeNames;

		for(size_t i=0; i<kminmersInfos.size(); i++){
			

			KmerVec vec = kminmers[i];
			//for(KmerVec& vec : kminmers){
			if(_mdbg->_dbg_nodes.find(vec) == _mdbg->_dbg_nodes.end()) continue;
			
			u_int32_t nodeName = _mdbg->_dbg_nodes[vec]._index;
			//_file_contigToNode << nodeName << "," << readIndex << endl;
			//if(nodeName == 933376) cout << "HAHAHAHA" << endl;

			//if(_contigFeature._nodeName_to_contigIndex.find(nodeName) == _contigFeature._nodeName_to_contigIndex.end()) continue;

			//_contigFeature._nodeNameDuplicate[nodeName] += 1;
			_contigFeature._nodeName_to_contigIndex[nodeName].push_back(readIndex);
			_contigFeature._contigIndex_to_nodeName[readIndex].push_back(nodeName);
			//_fileTestLala << nodeName << "," << "0" << endl;
			//_kminmerCounts[nodeName] = _countsInit;
			
			//nodeNames.insert(nodeName);
		}
		
		//std::sort(nodeNames.begin(), nodeNames.end());
		//_contigFeature._contigNodes[readIndex] = nodeNames;

	}





	static bool UnitigComparator_ByLength2(const UnitigLength &a, const UnitigLength &b){
		return a._length > b._length;
	}

	float _minUnitigAbundance;
	
	void execute_binning2(){






		_isFinalPass = false;

		vector<vector<string>> bins;
		
		_binOutputDir = _inputDir + "/binning";
	    if(fs::exists (_binOutputDir)){
			fs::remove_all(_binOutputDir);
        } 
		fs::create_directory(_binOutputDir);

		string clusterDir = _inputDir + "/" + "binGreedy";
		fs::path path(clusterDir);
	    if(!fs::exists (path)){
            fs::create_directory(path);
        } 

		string filename_binStats = clusterDir + "/binStats.txt";
		ofstream file_binStats(filename_binStats);

		ofstream fileHifiasmAll(clusterDir + "/component_hifiasm_all.csv");
		fileHifiasmAll << "Name,Colour" << endl;
		ofstream fileComponentNodeAll(clusterDir + "/component_node_all.csv");
		fileComponentNodeAll << "Name,Colour" << endl;


		u_int64_t clusterIndex = 0;
		u_int32_t contaminatedIndex = 0;
		unordered_set<string> written_unitigName;
		float prevCutoff = -1;

		vector<vector<u_int32_t>> clusters;
		vector<vector<u_int32_t>> clustersValid;
		bool reloadState = true;

		unordered_set<u_int32_t> processedNodeNames;
		unordered_set<u_int32_t> processedContigIndex;
		u_int64_t processedUnitigs = 0;

		u_int32_t binIndex = 0;

		ofstream file_bin_all(_inputDir + "/bin_all.csv");
		file_bin_all << "Name,Color" << endl;


		_fileOutput_contigBin = ofstream(_filename_outputBinning);
		_fileOutput_contigBin.close();


		vector<float> allCutoffs;

		for(const SaveState2& saveState : _graph->_cachedGraphStates){
			allCutoffs.push_back(saveState._abundanceCutoff_min);
			cout << saveState._abundanceCutoff_min << endl;
		}
		
		std::reverse(allCutoffs.begin(), allCutoffs.end());

		allCutoffs = {0};

		
		for(float cutoff : allCutoffs){

			_graph->loadState2(cutoff, -1, _unitigDatas);


			_contigIndex_to_componentIndex.clear();
			_componentIndex_to_contigIndexes.clear();
			unordered_set<u_int32_t> writtenUnitigs;
			unordered_set<u_int32_t> writtenContigs;
			u_int32_t componentIndex = 0;

			/*
			for(const Unitig& u : _graph->_unitigs){
				
				if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._startNode)) != writtenUnitigs.end()) continue;
				if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._endNode)) != writtenUnitigs.end()) continue;



				unordered_set<u_int32_t> component;
				_graph->getConnectedComponent_unitig(u._index, component);

				for(u_int32_t unitigIndex : component){
					
					const Unitig& u = _graph->_unitigs[unitigIndex];
					
					writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._startNode));
					writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._endNode));
					
					u_int32_t contigIndex = _contigFeature.nodepathToContigIndex(u._nodes);
					if(contigIndex == -1) continue;

					if(writtenContigs.find(contigIndex) != writtenContigs.end()) continue;

					_contigIndex_to_componentIndex[contigIndex] = componentIndex;
					_componentIndex_to_contigIndexes[componentIndex].push_back(contigIndex);

					writtenContigs.insert(contigIndex);
				}

				componentIndex += 1;
			}
			*/


			for(const Unitig& u: _graph->_unitigs){
				
				if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._startNode)) != writtenUnitigs.end()) continue;
				if(writtenUnitigs.find(BiGraph::nodeIndex_to_nodeName(u._endNode)) != writtenUnitigs.end()) continue;

				writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._startNode));
				writtenUnitigs.insert(BiGraph::nodeIndex_to_nodeName(u._endNode));

				u_int32_t contigIndex = _contigFeature.nodepathToContigIndex(u._nodes);
				if(contigIndex == -1) continue;
				
				if(writtenContigs.find(contigIndex) != writtenContigs.end()) continue;
				writtenContigs.insert(contigIndex);

				_contigIndex_to_componentIndex[contigIndex] = componentIndex;
				_componentIndex_to_contigIndexes[componentIndex].push_back(contigIndex);
			}

			cout << "Nb components: " << componentIndex << endl;
			//getchar();

			vector<UnitigLength> startingUnitigs;

			_minUnitigAbundance = cutoff / 0.2;

			vector<float> binningThresholds;
			if(_contigFeature._useAbundance){
				binningThresholds = {0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65};
			}
			else{
				binningThresholds = {0.99};
			}

			vector<u_int64_t> lengthThresholds = {100000, 50000, 10000}; 

			//0.01, 0.05, 0.5, 1.0, 2.0
			for(float binningThreshold : binningThresholds){
				for(u_int64_t lengthThreshold : lengthThresholds){ //, 10000, 2500
					
					u_int32_t binIndexValid = 0;

					if(binningThreshold == binningThresholds[binningThresholds.size()-1] && lengthThreshold == lengthThresholds[lengthThresholds.size()-1]){
						_isFinalPass = true;
					}
					_contigFeature._binningThreshold = binningThreshold;
					processedNodeNames.clear();
					processedContigIndex.clear();

					//cout << "STARTING PASS" << endl;
					//getchar();

					processedUnitigs = 0;
					startingUnitigs.clear();

			for(const Unitig& unitig: _graph->_unitigs){
				//if(unitig._nbNodes <= _kminmerSize*2) continue;
				//if(unitig._length < 100000) continue;

				if(unitig._abundance < _minUnitigAbundance) continue;

				u_int32_t contigIndex;
				u_int32_t contigLength;

				//string unitigSequence;
				//_toBasespace.createSequence(unitig_model._nodes, unitigSequence_model);
				bool isValid = _contigFeature.nodepathToContigLength(unitig._nodes, contigIndex, contigLength);
				if(!isValid) continue;
				if(contigLength < lengthThreshold) continue;

				if(processedContigIndex.find(contigIndex) != processedContigIndex.end()) continue;

				startingUnitigs.push_back({contigLength, unitig._abundance, unitig._startNode});

				//if(contigIndex == 209){
				//	cout << "loulou" << endl;
				//	getchar();
				//}
			}

			std::sort(startingUnitigs.begin(), startingUnitigs.end(), UnitigComparator_ByLength2);


				//if(binningThreshold != 0.01){
				dumpBins();
				_contigIndex_to_binIndex.clear();
				//_fileOutput_contigBin.close();
				_contigFeature.loadContigBins(_filename_outputBinning);
				//_fileOutput_contigBin = ofstream(_filename_outputBinning);
				cout << binningThreshold << endl;
					//getchar();
				//}


				
				for(const UnitigLength& unitigLength : startingUnitigs){

					u_int32_t source_nodeIndex = unitigLength._startNodeIndex;
					const Unitig& unitig = _graph->nodeIndex_to_unitig(source_nodeIndex);

					//cout << cutoff << " " << processedUnitigs << " " << startingUnitigs.size() << "     " << unitig._length << endl;
					processedUnitigs += 1;

					//if(isContigBinned(unitig._nodes, processedNodeNames)) continue;
					//if(isContigBinned(unitig._nodes, processedNodeNames)) continue;


					u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(unitig._startNode);

					//cout << processedUnitigs << "/" << startingUnitigs.size() << "    " << BiGraph::nodeIndex_to_nodeName(unitig._startNode) << " " << unitigLength._length << endl;




					//u_int32_t unitigIndex_model = _graph->nodeIndex_to_unitigIndex(unitig._startNode);
					//const Unitig& unitig_model = _graph->_unitigs[unitigIndex_model];


					//for(u_int32_t nodeIndex : unitig._nodes){
					//	processedNodeNames.insert(BiGraph::nodeIndex_to_nodeName(nodeIndex));
					//}

					binByReadpath(unitig._startNode, processedNodeNames, processedContigIndex, clusterDir, filename_binStats, fileHifiasmAll, fileComponentNodeAll, clusterIndex, binIndex, lengthThreshold, binIndexValid);


				}

				//cout << "END PASS" << endl;
				//getchar();
				//cout << "length done" << endl;
				//getchar();
				//_isFirstPass = false;
			}
			}

		}



		dumpBins();

		file_groundTruth.close();
		//file_groundTruth_hifiasmContigs.close();
		file_kminmersContigs.close();

		fileHifiasmAll.close();
		fileComponentNodeAll.close();
		file_binStats.close();
		file_bin_all.close();
		_fileOutput_contigBin.close();

	}

	void mergeBins(){

	}


};





#endif