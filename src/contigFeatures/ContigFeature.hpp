

#ifndef MDBG_METAG_CONTIGFEATURE
#define MDBG_METAG_CONTIGFEATURE

#include "../Commons.hpp"
#include <boost/math/distributions.hpp>

typedef double Distance;
typedef boost::math::normal_distribution<Distance> Normal;

static Distance minCV = 1;
static Distance minCVSum = 2;
const float MU_INTRA = 0;
const float SIGMA_INTRA = 0.01037897 / 2.0;
const float MU_INTER = 0.0676654;
const float SIGMA_INTER = 0.03419337;
const float VERY_SMALL_DOUBLE = 0.0000000001;
const float MAX_WEIGHT = std::numeric_limits<float>::max();

struct ContigFeatures{
	u_int32_t _unitigIndex;
	vector<float> _composition;
	vector<float> _abundance;
	vector<float> _abundanceVar;
};


struct BinData{
	vector<float> _sequenceModel;
};

class ContigFeature{

public:

	size_t _kmerSize;
	KmerModel* _kmerModel;
	u_int64_t _compositionVectorSize;
	vector<size_t> _kmerToCompositionIndex;
	string _filename_abundance;

	//vector<vector<float>> _comps;

	size_t _nbContigs;
	unordered_map<string, BinData> _bins;
	//unordered_map<string, vector<float>> _contigCompositions;

	u_int64_t _nbDatasets;
	unordered_map<string, vector<u_int32_t>> _contigToNodenames;

	unordered_map<u_int32_t, vector<u_int32_t>> _nodenameCounts;
	
	unordered_map<u_int32_t, vector<float>> _nodenameAbundanceMean;
	unordered_map<u_int32_t, vector<float>> _nodenameAbundanceVar;
	unordered_map<u_int32_t, u_int32_t> _nodeNameDuplicate;
	unordered_map<u_int32_t, string> _contigSequences;
	unordered_map<u_int32_t, u_int32_t> _contigIndex_to_binIndex;
	unordered_map<u_int32_t, vector<u_int32_t>> _binIndex_to_contigIndex;

	ContigFeature(){

		_nbDatasets = 0;

		_kmerSize = 4;
		_kmerModel = new KmerModel(_kmerSize);

		//int kmerSize = 4;
		//ModelCanonical model (kmerSize);
		//ModelCanonical::Iterator itKmer (model);

		unordered_map<u_int64_t, u_int64_t> setlala;

		int compositionIndex = 0;

		for(u_int64_t i=0 ;i<pow(4, _kmerSize); i++){
			KmerCanonical kmer;
			u_int64_t rev = KmerModel::revcomp(i, _kmerSize);
			//kmer += (u_int64_t) i;
			//coutu_int64_t v = i;
			kmer.set(i, rev);
			u_int64_t kmerVal = kmer.value();
			//KmerCanonical kmer_min = min(revcomp(kmer, kmerSize), kmer);

			if(setlala.find(kmerVal) == setlala.end()){
				setlala[kmerVal] = compositionIndex;
				_kmerToCompositionIndex.push_back(compositionIndex);
				compositionIndex += 1;
			}
			else{
				_kmerToCompositionIndex.push_back(setlala[kmerVal]);
			}
		}

		_compositionVectorSize = setlala.size();
		cout << "Kmer composition size: " << _compositionVectorSize << endl;

		//exit(1);
		setup();
	}

	//unordered_map<u_int32_t, unordered_set<u_int32_t>> _contigNodes;
	unordered_map<u_int32_t, vector<float>> _contigCoverages;
	unordered_map<u_int32_t, vector<float>> _contigCompositions;
	unordered_map<u_int32_t, vector<float>> _contigCoveragesVar;
	unordered_map<u_int32_t, vector<u_int32_t>> _nodeName_to_contigIndex;
	unordered_map<u_int32_t, vector<u_int32_t>> _contigIndex_to_nodeName;

	float _w_intra;
	float _w_inter;

	void setup(){
		float p_intra = 0.1;
		float p_inter = 0.01;
		float bin_threshold = -log10(p_intra);
		float break_threshold = -log10(p_inter);

		float w_intra = bin_threshold * (_nbDatasets + 1);
		float w_inter = break_threshold * (_nbDatasets + 1);

		//cout << bin_threshold << " " << break_threshold << endl;
		//cout << w_intra << " " << w_inter << endl;

		_w_intra = w_intra;
		_w_inter = w_inter;

		//_w_intra = 12;
		//_w_intra = 1;
		cout << "W_intra: " << bin_threshold << " " << break_threshold << endl;
	}

	void loadContigBins(const string& filename){
		
	    if(!fs::exists (filename)){
			cout << "No contig bin file" << endl;
			return;
		}

		ifstream file_contigBin(filename);

		while(true){

			u_int32_t contigIndex;
			u_int32_t binIndex;

			file_contigBin.read((char*)&contigIndex, sizeof(contigIndex));

			if(file_contigBin.eof()) break;
			
			file_contigBin.read((char*)&binIndex, sizeof(binIndex));

			_contigIndex_to_binIndex[contigIndex] = binIndex;
			_binIndex_to_contigIndex[binIndex].push_back(contigIndex);

			cout << contigIndex << " -> " << binIndex << endl;
		}


		file_contigBin.close();
	}

	void loadAbundanceFile(const string& filename){

		ifstream infile(filename);

		vector<string>* fields = new vector<string>();
		string line;
		std::getline(infile, line); //skip header

		while (std::getline(infile, line)){
			
			GfaParser::tokenize(line, fields, '\t');

			string contigName = (*fields)[0];

			size_t pos = contigName.find("ctg");
			contigName.erase(pos, 3);
			u_int32_t contigIndex = stoull(contigName);


			//cout << contigIndex << ": ";
			vector<float> coverages;
			for(size_t i=1; i<fields->size(); i++){
				const string& field = (*fields)[i];
				if(field.empty()) continue;
				//cout << i << " " << (*fields)[i] << endl;
				float ab = stof((*fields)[i]);
				//if(ab == 0) ab = 1;
				coverages.push_back(ab);
				//cout << ab << " ";
			}
			//cout << endl;

			_contigCoverages[contigIndex] = coverages;
			_nbDatasets = coverages.size();
			//u_int32_t nodeName = stoull((*fields)[0]);
			//u_int32_t contigIndex = stoull((*fields)[1]);

			//_contigToNodenames["ctg" + to_string(contigIndex)].push_back(nodeName);

			//cout << unitigIndex << " " << scgIndex << " " << clusterIndex << endl;


		}

		delete fields;
		infile.close();

	}
	
	void loadAbundanceFile_nodename(const string& filename){

		ifstream infile(filename);

		vector<string>* fields = new vector<string>();
		string line;
		std::getline(infile, line); //skip header

		while (std::getline(infile, line)){
			
			GfaParser::tokenize(line, fields, '\t');

			u_int32_t nodeName = stoull((*fields)[0]);

			vector<u_int32_t> coverages;
			for(size_t i=1; i<fields->size(); i++){
				//const string& field = (*fields)[i];
				//if(field.empty()) continue;
				//cout << i << " " << (*fields)[i] << endl;
				u_int32_t ab = stoull((*fields)[i]);
				//if(ab == 0) ab = 1;
				coverages.push_back(ab);
			}

			_nodenameCounts[nodeName] = coverages;
			_nbDatasets = coverages.size();
			//u_int32_t nodeName = stoull((*fields)[0]);
			//u_int32_t contigIndex = stoull((*fields)[1]);

			//_contigToNodenames["ctg" + to_string(contigIndex)].push_back(nodeName);

			//cout << unitigIndex << " " << scgIndex << " " << clusterIndex << endl;


		}

		delete fields;
		infile.close();

		setup();
	}

	void loadAbundanceFile_metabat(const string& filename){

		if(filename.empty()){
			cout << "No abundance file provided" << endl;
			return;
		}

		ifstream infile(filename);

		vector<string>* fields = new vector<string>();
		string line;
		std::getline(infile, line); //skip header

		while (std::getline(infile, line)){
			
			GfaParser::tokenize(line, fields, '\t');

			string contigName = (*fields)[0];

			size_t pos = contigName.find("ctg");
			contigName.erase(pos, 3);
			u_int32_t contigIndex = stoull(contigName);

			//vector<u_int32_t> nodes;
			//nodes = contigToNodenames[contigName];
			//if(contigToNodeIndex.find(contigName) != contigToNodeIndex.end()){
			//}

			vector<float> abundanceMean;
			vector<float> abundanceVar;

			for(size_t i=1; i<fields->size(); i++){ //metabat i=3
				const string& field = (*fields)[i];
				if(field.empty()) continue;
				//cout << i << " " << (*fields)[i] << endl;
				float ab = stof((*fields)[i]);
				//if(ab == 0) ab = 1;
				if(i%2 == 0){
					abundanceVar.push_back(ab);
					//abundanceVar.push_back(10);
					//cout << "Var: " << ab << endl;
				}
				else{
					//abundanceVar.push_back(ab);
					abundanceMean.push_back(ab);
					//cout << "Mean: " << ab << endl;
				}
			}

			_contigCoverages[contigIndex] = abundanceMean;
			_contigCoveragesVar[contigIndex] = abundanceVar;

			_nbDatasets = abundanceMean.size();

			/*
			for(u_int32_t nodeName : nodes){
				//u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
				_nodenameAbundanceMean[nodeName] = abundanceMean;
				_nodenameAbundanceVar[nodeName] = abundanceVar;
			}
			*/

			//u_int32_t nodeName = stoull((*fields)[0]);
			//u_int32_t contigIndex = stoull((*fields)[1]);

			//_contigToNodenames["ctg" + to_string(contigIndex)].push_back(nodeName);

			//cout << unitigIndex << " " << scgIndex << " " << clusterIndex << endl;
			//exit(1);

		}

		delete fields;
		infile.close();

		setup();
	}
	
	u_int32_t contigIndexToBinIndex(u_int32_t contigIndex){
		
		if(_contigIndex_to_binIndex.find(contigIndex) == _contigIndex_to_binIndex.end()) return -1;

		return _contigIndex_to_binIndex[contigIndex];
	}
	void sequenceToComposition(const string& sequence, vector<float>& composition){
		
		composition.resize(_compositionVectorSize);
		u_int64_t unitigLength = sequence.size();

		if(sequence.size() < _kmerSize) return;

		//cout << sequence.size() << endl;
		//size_t i=0;
		//itKmer.setData (sequence.getData());

		vector<u_int64_t> kmers;
		_kmerModel->iterate(sequence.c_str(), sequence.size(), kmers);

		//size_t overlap_size = 0;

		for (u_int64_t kmer : kmers){
			composition[_kmerToCompositionIndex[kmer]] += 1;
		}

		/*
		float rsum = 0;
		for(size_t i = 0; i < composition.size(); ++i) {
			rsum += composition[i] * composition[i];
		}
		rsum = sqrt(rsum);
		for(size_t i = 0; i < composition.size(); ++i) {
			composition[i] /= rsum;
		}
		*/
		
		float sum = 0;
		for(size_t i = 0; i < composition.size(); ++i) {
			sum += composition[i];
		}

		if(sum == 0) return;

		for(size_t i = 0; i < composition.size(); ++i) {
			composition[i] /= sum;
		}
		
	}

	bool nodepathToComposition(const vector<u_int32_t>& nodepath, vector<float>& composition){
		
		u_int32_t contigIndex = nodepathToContigIndex(nodepath);
		if(contigIndex == -1) return false;

		composition = _contigCompositions[contigIndex];
		
		return true;
		/*
		vector<vector<float>> values_mean;
		values_mean.resize(_compositionVectorSize);
		//vector<vector<float>> values_var;
		//values_var.resize(_compositionVectorSize);

		composition.clear();
		composition.resize(_compositionVectorSize, 0);

		size_t n = 0;
		for(u_int32_t nodeIndex : sequence){
			u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
			
			if(_nodeName_to_contigIndex.find(nodeName) == _nodeName_to_contigIndex.end()) continue;
			
			u_int32_t contigIndex = _nodeName_to_contigIndex[nodeName];

			const vector<float>& contigComposition = _contigCompositions[contigIndex];

			for(size_t i=0; i<contigComposition.size(); i++){
				//values_mean[i].push_back(contigComposition[i]);
				composition[i] += contigComposition[i];
			}
			
			n += 1;
		}

		if(n == 0) return false;

		for(size_t i=0; i<composition.size(); i++){
			//composition[i] = Utils::compute_median_float(values_mean[i]);
			composition[i] /= n;;
		}

		return true;
		*/
	}

	bool nodepathToContigSequence(const vector<u_int32_t>& nodepath, string& sequence, u_int32_t& contigIndexResult){
		
		u_int32_t contigIndex = nodepathToContigIndex(nodepath);
		if(contigIndex == -1) return false;

		contigIndexResult = contigIndex;
		sequence = _contigSequences[contigIndex];
		
		return true;

		/*
		unordered_map<u_int32_t, u_int32_t> contigCounts;

		u_int32_t existingContigIndex = -1;

		contigIndexResult = -1;
		sequence.clear();

		for(u_int32_t nodeIndex : nodepath){
			u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
			
			if(_nodeName_to_contigIndex.find(nodeName) == _nodeName_to_contigIndex.end()) continue;
			
			u_int32_t contigIndex = _nodeName_to_contigIndex[nodeName];

			
			if(existingContigIndex == -1){
				existingContigIndex = contigIndex;
			}
			else{
				if(existingContigIndex != contigIndex){ //repeated node
					sequence.clear();
					return false;
				}
			}

			contigIndexResult = contigIndex;
			sequence = _contigSequences[contigIndex];
			
			contigCounts[contigIndex]  += 1;
		}

		u_int32_t maxCount = 0;
		for(const auto& it : contigCounts){
			if(it.second > maxCount){
				maxCount = it.second;
			}
		}

		u_int32_t maxContigIndex = -1;
		u_int32_t nbMaxCount = 0;
		for(const auto& it : contigCounts){
			if(it.second == maxCount){
				nbMaxCount += 1;
				maxContigIndex = it.first;
			}
		}

		//cout << nbMaxCount << endl;

		if(nbMaxCount == 1){
			contigIndexResult = maxContigIndex;
			sequence = _contigSequences[maxContigIndex];
			return true;
		}

		return false;
		*/
	}

	u_int32_t nodepathToContigIndex(const vector<u_int32_t>& nodepath){

		unordered_map<u_int32_t, u_int32_t> contigCounts;



		for(u_int32_t nodeIndex : nodepath){
			u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
			
			if(_nodeName_to_contigIndex.find(nodeName) == _nodeName_to_contigIndex.end()) continue;
			
			const vector<u_int32_t>& contigIndexes = _nodeName_to_contigIndex[nodeName];

			for(u_int32_t contigIndex : contigIndexes){
				contigCounts[contigIndex]  += 1;
			}
		}

		u_int32_t maxCount = 0;
		for(const auto& it : contigCounts){
			if(it.second > maxCount){
				maxCount = it.second;
			}
		}

		u_int32_t maxContigIndex = -1;
		u_int32_t nbMaxCount = 0;
		for(const auto& it : contigCounts){
			if(it.second == maxCount){
				nbMaxCount += 1;
				maxContigIndex = it.first;
			}
		}

		//cout << nbMaxCount << endl;

		if(nbMaxCount == 1){
			return maxContigIndex;
		}

		return -1;

	}

	bool sequenceToAbundance(const vector<u_int32_t>& nodepath, vector<float>& abundances, vector<float>& abundancesVar){

		abundances.clear();
		abundancesVar.clear();

		u_int32_t contigIndex = nodepathToContigIndex(nodepath);
		if(contigIndex == -1) return false;

		abundances = _contigCoverages[contigIndex];
		abundancesVar = _contigCoveragesVar[contigIndex];
		
		float sum = 0.0;
		for(float val : abundances){
			sum += val;
		}
		if(sum == 0) return false;


		return true;


		/*
		vector<vector<float>> values_mean;
		values_mean.resize(_nbDatasets);
		vector<vector<float>> values_var;
		values_var.resize(_nbDatasets);

		abundances.clear();
		abundances.resize(_nbDatasets, 0);
		abundancesVar.clear();
		abundancesVar.resize(_nbDatasets, 0);

		size_t n = 0;
		for(u_int32_t nodeIndex : sequence){
			u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
			
			
			//cout << nodeName << " " << (_nodeName_to_contigIndex.find(nodeName) != _nodeName_to_contigIndex.end()) << endl;
			if(_nodeName_to_contigIndex.find(nodeName) == _nodeName_to_contigIndex.end()) continue;
			//if(_nodeNameDuplicate.find(nodeName) != _nodeNameDuplicate.end() && _nodeNameDuplicate[nodeName] != 1) continue;
			u_int32_t contigIndex = _nodeName_to_contigIndex[nodeName];

			const vector<float>& contigAbundances = _contigCoverages[contigIndex];
			const vector<float>& contigAbundancesVar = _contigCoveragesVar[contigIndex];

			for(size_t i=0; i<contigAbundances.size(); i++){
				abundances[i] += contigAbundances[i];
				abundancesVar[i] += contigAbundancesVar[i];
				//values_mean[i].push_back(contigAbundances[i]);
			}
			
			n += 1;
		}

		if(n == 0) return false;

		for(size_t i=0; i<abundances.size(); i++){
			abundances[i] /= n;
			if(abundances[i] < 1) abundances[i] = 0;
			//abundances[i] = Utils::compute_median_float(values_mean[i]);
			//abundancesVar[i] = abundances[i];
			abundancesVar[i] /= n;
			if(abundancesVar[i] < 1) abundancesVar[i] = 0;
		}


		float sum = 0.0;
		for(float val : abundances){
			sum += val;
		}
		if(sum == 0) return false;

		return true;
		*/
		/*
		if(n <= 1){
			for(size_t i=0; i<abundancesVar.size(); i++){
				abundancesVar[i] = abundances[i];
			}
		}
		else{
			for(u_int32_t nodeIndex : sequence){
				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
				
				if(_nodeName_to_contigIndex.find(nodeName) == _nodeName_to_contigIndex.end()) continue;
				
				u_int32_t contigIndex = _nodeName_to_contigIndex[nodeName];

				const vector<float>& contigAbundances = _contigCoverages[contigIndex];

				//const vector<u_int32_t>& counts = _nodenameCounts[nodeName];

				for(size_t i=0; i<contigAbundances.size(); i++){
					abundancesVar[i] += (contigAbundances[i] - abundances[i]) * (contigAbundances[i] - abundances[i]);
				}
			}

			for(u_int32_t i=0; i<abundances.size(); i++){
				abundancesVar[i] /= ((float) n-1);
			}
		}
		*/


		/* METABAT 
		for(u_int32_t nodeIndex : sequence){
			u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
			
			if(_nodenameAbundanceMean.find(nodeName) == _nodenameAbundanceMean.end()) continue;
			
			const vector<float>& means = _nodenameAbundanceMean[nodeName];

			for(size_t i=0; i<means.size(); i++){
				//forabundances[i] += counts[i];
				values_mean[i].push_back(means[i]);
			}

			const vector<float>& vars = _nodenameAbundanceVar[nodeName];

			for(size_t i=0; i<vars.size(); i++){
				//forabundances[i] += counts[i];
				values_var[i].push_back(vars[i]);
			}

		}

		for(u_int32_t i=0; i<abundances.size(); i++){
			abundances[i] = Utils::compute_median_float(values_mean[i]);
		}
		for(u_int32_t i=0; i<abundancesVar.size(); i++){
			abundancesVar[i] = Utils::compute_median_float(values_var[i]);
		}
		*/ //METABAT

		/*
		vector<vector<u_int32_t>> values;
		values.resize(_nbDatasets);

		//cout << "----" << endl;
		abundances.clear();
		abundances.resize(_nbDatasets, 0);
		abundancesVar.clear();
		abundancesVar.resize(_nbDatasets, 0);
		//vector<float> abundances(_nbDatasets, 0);

		for(u_int32_t nodeIndex : sequence){
			u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
			
			const vector<u_int32_t>& counts = _nodenameCounts[nodeName];
			for(size_t i=0; i<counts.size(); i++){
				//abundances[i] += counts[i];
				values[i].push_back(counts[i]);
			}
		}

		
		for(u_int32_t i=0; i<abundances.size(); i++){
			//for(size_t j=0; j<sequence.size(); j++){
			//	cout << values[i][j] << " ";
			//}
			//cout << endl;
			//cout << abundances[i] << " " << sequence.size() << " " << (abundances[i] / sequence.size()) << endl;
			//abundances[i] /= ((float) sequence.size());
			abundances[i] = Utils::compute_median(values[i]);
			if(abundances[i] == 0) abundances[i] = VERY_SMALL_DOUBLE;
		}
		*/
	
		/*
		if(sequence.size() > 1){
			for(u_int32_t nodeIndex : sequence){
				u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
				
				const vector<u_int32_t>& counts = _nodenameCounts[nodeName];

				for(size_t i=0; i<counts.size(); i++){
					abundancesVar[i] += (counts[i] - abundances[i]) * (counts[i] - abundances[i]);
				}
			}

			for(u_int32_t i=0; i<abundances.size(); i++){
				abundancesVar[i] /= ((float) sequence.size()-1);
			}
		}
		*/

		//cout << sequence.size() << endl;
		
	}


	/*
	void computeFastaComposition(const string& sequenceFilename){


		float p_intra = 0.01;
		float p_inter = 0.01;
		float bin_threshold = -log10(p_intra);
		float break_threshold = -log10(p_inter);

		float w_intra = bin_threshold * (_nbDatasets + 1);
		float w_inter = break_threshold * (_nbDatasets + 1);

		cout << bin_threshold << " " << break_threshold << endl;
		cout << w_intra << " " << w_inter << endl;

		ifstream infile("/home/gats/workspace/run/overlap_test_AD_k7/binByCutoff/component_16_unitigColor.csv");
		ofstream outfile("/home/gats/workspace/run/overlap_test_AD_k7/binByCutoff/component_16_unitigClusterNew.csv");
		outfile << "Name,Color" << endl;

		vector<string>* fields = new vector<string>();
		string line;
		std::getline(infile, line); //skip header

		while (std::getline(infile, line)){
			
			GfaParser::tokenize(line, fields, ',');

			u_int32_t nodeName = stoull((*fields)[0]);
			u_int32_t contigIndex = stoull((*fields)[1]);

			_contigToNodenames["ctg" + to_string(contigIndex)].push_back(nodeName);

			//cout << unitigIndex << " " << scgIndex << " " << clusterIndex << endl;
		}

		delete fields;
		infile.close();


		auto fp = std::bind(&ContigFeature::computeFastaComposition_read, this, std::placeholders::_1, std::placeholders::_2);
		ReadParser readParser(sequenceFilename, true);
		readParser.parse(fp);

		cout << _contigCompositions.size() << " " << _bins.size() << endl;
		for(auto& it : _contigCompositions){
			const string& name = it.first;
			const vector<float>& composition = it.second;

			const vector<float>& contigCoverages = _contigCoverages[name];

			if(_bins.find(name) != _bins.end()) continue;

			float maxProb = 9999999;
			string maxBinName = "";

			cout << "Contig: " << name << endl;

			for(auto& it : _bins){

				const string& binName = it.first;
				const BinData& binData = it.second;

				const vector<float>& binCoverages = _contigCoverages[binName];

				float prob_comp = computeCompositionProbability(binData._sequenceModel, composition);
				//float log_prob = -log10(prob_comp);
				//float prob_comp = computeEuclideanDistance(binData._sequenceModel, composition);
				float prob_cov = computeAbundanceProbability(binCoverages, contigCoverages);
				//prob_cov = abs(prob_comp);
				//float log_prob = prob_comp;


				float prob_product = prob_comp * prob_cov;

				float log_prob = 0;

				if (prob_product > 0.0){
					log_prob = - (log10(prob_comp) + log10(prob_cov));
				}

				cout << "\t" << binName << ": " << log_prob << "    " << (-log10(prob_comp)) << " " << (-log10(prob_cov)) << "    " << computeEuclideanDistance(binCoverages, contigCoverages) << endl;
				
				cout << "\t";
				for(float ab : binCoverages) cout << ab << " ";
				cout << endl;
				cout << "\t";
				for(float ab : contigCoverages) cout << ab << " ";
				cout << endl;

				if(log_prob <= w_intra){
					if(log_prob != 0 && log_prob < maxProb){
						maxProb = log_prob;
						maxBinName = binName;
					}
				}

			}

			if(!maxBinName.empty()){
				cout << name << " -> " << maxBinName << endl;
				for(u_int32_t nodeName : _contigToNodenames[name]){
					outfile << nodeName << "," << maxBinName << endl;
				}
				//getchar();

			}



			//if(name == "ctg4893") getchar();


		}


		for(auto& it : _bins){
			const string& binName = it.first;
			for(u_int32_t nodeName : _contigToNodenames[binName]){
				outfile << nodeName << "," << binName << endl;
			}
		}
		
		outfile.close();
		//cout << _nbContigs << endl;
		exit(1);

	}


	void computeFastaComposition_read(kseq_t* read, u_int64_t readIndex){

		string sequence = string(read->seq.s, strlen(read->seq.s));
		
		vector<float> composition;
		sequenceToComposition(sequence, composition);

		string header = string(read->name.s, strlen(read->name.s));

		if(header == "ctg1539" || header == "ctg482" || header == "ctg1629" || header == "ctg166" || header == "ctg286" || header == "ctg648" || header == "ctg1165"){
		//if(header == "ctg314" || header == "ctg1" || header == "ctg625" || header == "ctg1464" || header == "ctg671" || header == "ctg888" || header == "ctg1165"){
			if(_bins.find(header) == _bins.end()){
				_bins[header] = {composition};
				cout << header << endl;
			}
		}

		_contigCompositions[header] = composition;

		//if( )
		//for(float v : composition){
		//	cout << v << " ";
		//}
		//cout << endl;

		//exit(1);
		//_comps.push_back(composition);
		_nbContigs += 1;
	}
	*/

	inline float computeEuclideanDistance(const vector<float>& v1, const vector<float>& v2){
		float sum = 0;
		for(int i=0; i<v1.size(); i++) {
			sum += pow(v2[i] - v1[i], 2);
		}
		return sqrt(sum);
	}

	float computeNormpdf(float x, float mean, float sd){
		float var = pow(sd, 2);
    	float denom = sd*pow(2*M_PI, 0.5);
		float num = exp(-pow((x-mean), 2)/(2*var));
		return num/denom;
	}


	float computeCompositionProbability(const vector<float>& v1, const vector<float>& v2){
		float dist = computeEuclideanDistance(v1, v2);
		float gaus_intra = computeNormpdf(dist, MU_INTRA, SIGMA_INTRA);
		float gaus_inter = computeNormpdf(dist, MU_INTER, SIGMA_INTER);
		return gaus_intra/(gaus_intra+gaus_inter);
	}

	double computeAbundanceCorrelation(const vector<float>& cov1, const vector<float>& cov2) {

		size_t i, ii;
		double sum_xsq = 0.0;
		double sum_ysq = 0.0;
		double sum_cross = 0.0;
		double ratio;
		double delta_x, delta_y;
		double mean_x = 0.0, mean_y = 0.0;
		double r = 0.0;

		size_t s = 0; //skipped

		for (i = 0; i < cov1.size(); ++i) {
			double m1 = cov1[i]; //ABD(r1,i);
			double m2 = cov2[i]; //is_small ? small_ABD(r2, i) : ABD(r2, i);
			//m2 = rand()%500;

			//cout << m1 << " " << m2 << endl;
			ii = i - s;

			if(ii == 0) {
				mean_x = m1;
				mean_y = m2;
				continue;
			}

			ratio = ii / (ii + 1.0);
			delta_x = m1 - mean_x;
			delta_y = m2 - mean_y;
			sum_xsq += delta_x * delta_x * ratio;
			sum_ysq += delta_y * delta_y * ratio;
			sum_cross += delta_x * delta_y * ratio;
			mean_x += delta_x / (ii + 1.0);
			mean_y += delta_y / (ii + 1.0);
		}

		r = sum_cross / (sqrt(sum_xsq) * sqrt(sum_ysq));

		return r;
	}



	float computeAbundanceProbability(const vector<float>& cov1, const vector<float>& cov2){

		float poisson_prod_1 = 1;
		float poisson_prod_2 = 1;

		for(size_t i=0; i<cov1.size(); i++){
			
			float poisson_pmf_1 = exp((cov1[i] * log(cov2[i])) - lgamma(cov1[i] + 1.0) - cov2[i]);

			float poisson_pmf_2 = exp((cov2[i] * log(cov1[i])) - lgamma(cov2[i] + 1.0) - cov1[i]);

			if (poisson_pmf_1 < VERY_SMALL_DOUBLE){
				poisson_pmf_1 = VERY_SMALL_DOUBLE;
			}

			if (poisson_pmf_2 < VERY_SMALL_DOUBLE){
				poisson_pmf_2 = VERY_SMALL_DOUBLE;
			}

			poisson_prod_1 = poisson_prod_1 * poisson_pmf_1;

			poisson_prod_2 = poisson_prod_2 * poisson_pmf_2;
		}

		return min(poisson_prod_1, poisson_prod_2);
	}
	
	float computeAbundanceProbability_new(const vector<float>& cov1, const vector<float>& cov2){

		vector<float> means_f1 = cov1;
		vector<float> means_f2 = cov2;

		//for(size_t i=0; i<f1._abundance.size(); i++){
		//	if(means_f1[i] < minCV) means_f1[i] = 0;
		//	if(means_f2[i] < minCV) means_f2[i] = 0;
		//}

		int nnz = 0;

		float mean_ratio = 0;
		for(size_t i=0; i<means_f1.size(); i++){
			if (means_f1[i] > 0 || means_f2[i] > 0) {
				if(means_f2[i] > 0){
					float ratio = means_f1[i] / means_f2[i];
					mean_ratio += ratio;
				}
				nnz += 1;
			}
		}

		if(nnz == 0) return 1;

		mean_ratio /= nnz;
		
		for(size_t i=0; i<means_f1.size(); i++){
			means_f2[i] *= mean_ratio;
		}


		float poisson_prod_1 = 1;
		float poisson_prod_2 = 1;

		for(size_t i=0; i<means_f1.size(); i++){

			float m1 = std::max(means_f1[i], (float)0.000001);
			float m2 = std::max(means_f2[i], (float)0.000001);

			float poisson_pmf_1 = exp((m1 * log(m2)) - lgamma(m1 + 1.0) - m2);

			float poisson_pmf_2 = exp((m2 * log(m1)) - lgamma(m2 + 1.0) - m1);

			if (poisson_pmf_1 < VERY_SMALL_DOUBLE){
				poisson_pmf_1 = VERY_SMALL_DOUBLE;
			}

			if (poisson_pmf_2 < VERY_SMALL_DOUBLE){
				poisson_pmf_2 = VERY_SMALL_DOUBLE;
			}

			poisson_prod_1 = poisson_prod_1 * poisson_pmf_1;

			poisson_prod_2 = poisson_prod_2 * poisson_pmf_2;
		}

		return min(poisson_prod_1, poisson_prod_2);
	}


	float computeProbability(const ContigFeatures& f1, const ContigFeatures& f2){

		float prob_comp = computeCompositionProbability(f1._composition, f2._composition);
		//return - log10(prob_comp);

		float prob_cov = computeAbundanceProbability(f1._abundance, f2._abundance);

		float prob_product = prob_comp * prob_cov;

		float log_prob = 0;

		if (prob_product > 0.0){
			log_prob = - (log10(prob_comp) + log10(prob_cov));
		}
		else{
			log_prob = MAX_WEIGHT;
		}

		return log_prob;

	}

	/*
	bool isIntra(const ContigFeatures& f1, const ContigFeatures& f2, bool hasComposition, bool hasAbundances){


		float compositionProb = -log10(computeCompositionProbability(f1._composition, f2._composition));
		if(isinf(compositionProb)) return false;
		return compositionProb < 0.05;


		if(!hasAbundances) return false;
		
		if(hasComposition){
			float compositionProb = -log10(computeCompositionProbability(f1._composition, f2._composition));
			
			if(!isinf(compositionProb)){
				//if(compositionProb > 10) return false;
			}
			//if(isnan(compositionProb)) return false;

		}

		
		int nnz = 0;

		//cout << cal_abd_dist_new(f1, f2 ,nnz) << endl;
		//cout << isinf(cal_abd_dist_new(f1, f2 ,nnz)) << endl;
		
		float dist = cal_abd_dist_new(f1, f2 ,nnz);
		//cout << isinf(dist) << endl;
		//cout << fpclassify(dist) << endl;
		//cout << (fpclassify(dist) == FP_INFINITE) << endl;


		if(isinf(dist)) return false;
		if(isnan(dist)) return false;

		float cor = computeAbundanceCorrelation(f1._abundance, f2._abundance);
		if(isinf(cor)) return false;
		if(isnan(cor)) return false;
		
		float tnf_dist = cal_tnf_dist(f1._composition, f2._composition, f1._unitigIndex, f2._unitigIndex);

		//(1-tnf_dist)
		return  cor * (1-dist) > 0.65;


		//return  computeAbundanceCorrelation(f1._abundance, f2._abundance) > 0.95 && -log10(computeAbundanceProbability(f1._abundance, f2._abundance)) < 20;
		//int nnz = 0;
		//return  computeAbundanceCorrelation(f1._abundance, f2._abundance) > 0.95 && cal_abd_dist(f1, f2 ,nnz) < 0.3;// && computeCompositionProbability(f1._composition, f2._composition) < 2;

		//float prob_cov = computeAbundanceProbability(f1._abundance, f2._abundance);
		//return computeProbability(f1, f2) < 30;
		

		//int nnz = 0;
		//return computeAbundanceCorrelation(f1._abundance, f2._abundance) > 0.98 && cal_abd_dist(f1, f2 ,nnz) < 0.3;
		//float prob = computeProbability(f1 , f2);
		//cout << "\tProb: " << prob << " " << _w_intra << endl;
		//return prob < _w_intra;
	}
	*/

	bool isIntra(const vector<u_int32_t>& bin1, const vector<u_int32_t>& bin2){

		cout << "\tIs intra: " << bin1.size() << " " << bin2.size() << endl;

		float distance = computeDistance(bin1, bin2);

		if(isinf(distance)) return false;

		//return distance < 0.015;
		//return distance < 0.015;
		return distance < 2;

		//(1-tnf_dist)
		//return  cor * (1-dist) > 0.65;
	}

	float computeDistance(const vector<u_int32_t>& bin1, const vector<u_int32_t>& bin2){

		double distance_max = 0;
		double distance_sum = 0;
		double distance_n = 0;

		for(u_int32_t contigIndex1 : bin1){
			for(u_int32_t contigIndex2 : bin2){
				float distance = computeDistance(contigIndex1, contigIndex2);
				distance_sum += distance;
				distance_n += 1;

				if(distance > distance_max){
					distance_max = distance;
				}
			}
		}

		cout << "\tComposition distance mean: " << (distance_sum / distance_n) << endl;

		//return distance_sum / distance_n;
		return distance_max;
	}

	float computeDistance(u_int32_t contigIndex1, u_int32_t contigIndex2){

		const vector<float>& composition1 = _contigCompositions[contigIndex1];
		const vector<float>& composition2 = _contigCompositions[contigIndex2];

		float compositionProb = -log10(computeCompositionProbability(composition1, composition2));

		cout << "\tComposition distance: " << compositionProb << endl;
		//const vector<float>& abundance1 = _contigCoverages[contigIndex1];
		//const vector<float>& abundance1_var = _contigCoveragesVar[contigIndex1];
		//const vector<float>& abundance2 = _contigCoverages[contigIndex2];
		//const vector<float>& abundance2_var = _contigCoveragesVar[contigIndex2];

		//int nnz = 0;
		//float dist_abundance = cal_abd_dist_new(abundance1, abundance1_var,  abundance2, abundance2_var, nnz);
		//if(compositionProb < 0.05){
			//cout << dist_abundance << endl;
		//}
		
		return compositionProb;
	}

	// for normal distributions
	Distance cal_abd_dist2(Normal& p1, Normal& p2) {
		Distance k1, k2, tmp, d = 0;

		Distance m1 = p1.mean();
		Distance m2 = p2.mean();
		Distance v1 = p1.standard_deviation();
		v1 = v1 * v1;
		Distance v2 = p2.standard_deviation();
		v2 = v2 * v2;

		//normal_distribution
		if (fabs(v2 - v1) < 1e-4) {
			k1 = k2 = (m1 + m2) / 2;
		} else {
			tmp = sqrt(v1 * v2 * ((m1 - m2) * (m1 - m2) - 2 * (v1 - v2) * log(sqrt(v2 / v1))));
			k1 = (tmp - m1 * v2 + m2 * v1) / (v1 - v2);
			k2 = (tmp + m1 * v2 - m2 * v1) / (v2 - v1);
		}

		if (k1 > k2) {
			tmp = k1;
			k1 = k2;
			k2 = tmp;
		}
		if (v1 > v2) {
			std::swap(p1, p2);
		}

		if (k1 == k2)
			d = (fabs(boost::math::cdf(p1, k1) - boost::math::cdf(p2, k1)));
		else
			d = (fabs(boost::math::cdf(p1, k2) - boost::math::cdf(p1, k1) + boost::math::cdf(p2, k1) - boost::math::cdf(p2, k2)));

		return d;

	}

	Distance cal_abd_dist(const ContigFeatures& f1, const ContigFeatures& f2, int& nnz) {

		float distSum = 0.0f;
		nnz = 0;

		for (size_t i=0; i < f1._abundance.size(); ++i) {

			Distance d = 0;
			Distance m1 = f1._abundance[i];
			Distance m2 = f2._abundance[i];
			if (m1 > minCV || m2 > minCV) {
				//nz = true;
				m1 = std::max(m1, (double)0.000001);
				m2 = std::max(m2, (double)0.000001);
				if (m1 != m2) {
					Distance v1 = f1._abundanceVar[i] < 1 ? 1 : f1._abundanceVar[i];
					Distance v2 = f2._abundanceVar[i] < 1 ? 1 : f2._abundanceVar[i];

					//cout << m1 << " " << m2 << " " << v1 << " " << v2 << endl;
					Normal p1(m1, sqrt(v1)), p2(m2, sqrt(v2));
					d = cal_abd_dist2(p1, p2);
				}
				nnz += 1;
			}
			distSum += std::min(std::max(d, 1e-6), 1. - 1e-6);

		}

		//return POW(EXP(distSum), 1.0 / nnz);
		return distSum / nnz;
	}

	Distance cal_abd_dist_new(const vector<float>& abundance1, const vector<float>& abundance1_var, const vector<float>& abundance2, const vector<float>& abundance2_var, int& nnz) {

		vector<float> means_f1 = abundance1;
		vector<float> means_f2 = abundance2;

		//for(size_t i=0; i<f1._abundance.size(); i++){
		//	if(means_f1[i] < minCV) means_f1[i] = 0;
		//	if(means_f2[i] < minCV) means_f2[i] = 0;
		//}

		nnz = 0;

		//cout << endl;
		float mean_ratio = 0;
		for(size_t i=0; i<abundance1.size(); i++){
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

		if(nnz == 0) return 1;
		
		mean_ratio /= nnz;
		//cout << mean_ratio << endl;
		
		for(size_t i=0; i<abundance1.size(); i++){
			means_f2[i] *= mean_ratio;
		}


		float distSum = 0.0f;
		nnz = 0;

		for (size_t i=0; i < abundance1.size(); ++i) {

			Distance d = 0;
			Distance m1 = means_f1[i];
			Distance m2 = means_f2[i];
			if (m1 > minCV || m2 > minCV) {
				//nz = true;
				m1 = std::max(m1, (double)0.000001);
				m2 = std::max(m2, (double)0.000001);
				if (m1 != m2) {
					Distance v1 = abundance1_var[i] < 1 ? 1 : abundance1_var[i];
					Distance v2 = abundance2_var[i] < 1 ? 1 : abundance2_var[i];

					//cout << m1 << " " << m2 << " " << v1 << " " << v2 << endl;
					Normal p1(m1, sqrt(v1)), p2(m2, sqrt(v2));
					d = cal_abd_dist2(p1, p2);
				}
				nnz += 1;
			}
			distSum += std::min(std::max(d, 1e-6), 1. - 1e-6);

		}

		//return POW(EXP(distSum), 1.0 / nnz);
		return distSum / nnz;
	}

	Distance cal_tnf_dist(const vector<float>& c1, const vector<float>& c2, u_int32_t contigIndex1, u_int32_t contigIndex2) {

		size_t length1 = _contigSequences[contigIndex1].size();
		size_t length2 = _contigSequences[contigIndex2].size();
		//if(d1._length < 2500) return 1;
		//if(d2._length < 2500) return 1;
		/*
		Distance d = 0;

		for (size_t i = 0; i < nTNF; ++i) {
			d += (TNF(r1,i) - TNF(r2,i)) * (TNF(r1,i) - TNF(r2,i)); //euclidean distance
		}

		d = SQRT(d);
		*/

		Distance d =  computeEuclideanDistance(c1, c2);
		//Distance d = getCompositionDistance(r1, r2);
		//cout << d << endl;
		Distance b,c; //parameters

		size_t ctg1 = std::min(length1, (size_t)500000);
		size_t ctg2 = std::min(length2, (size_t)500000);
		//cout << ctg1 << " " << ctg2 << endl;
		Distance lw11 = log10(std::min(ctg1, ctg2));
		Distance lw21 = log10(std::max(ctg1, ctg2));
		Distance lw12 = lw11 * lw11;
		Distance lw13 = lw12 * lw11;
		Distance lw14 = lw13 * lw11;
		Distance lw15 = lw14 * lw11;
		Distance lw16 = lw15 * lw11;
		Distance lw17 = lw16 * lw11;
		Distance lw22 = lw21 * lw21;
		Distance lw23 = lw22 * lw21;
		Distance lw24 = lw23 * lw21;
		Distance lw25 = lw24 * lw21;
		Distance lw26 = lw25 * lw21;

		Distance prob;

		b = 46349.1624324381 + -76092.3748553155*lw11 + -639.918334183*lw21 + 53873.3933743949*lw12 + -156.6547554844*lw22 + -21263.6010657275*lw13 + 64.7719132839*lw23 +
				5003.2646455284*lw14 + -8.5014386744*lw24 + -700.5825500292*lw15 + 0.3968284526*lw25 + 54.037542743*lw16 + -1.7713972342*lw17 + 474.0850141891*lw11*lw21 +
				-23.966597785*lw12*lw22 + 0.7800219061*lw13*lw23 + -0.0138723693*lw14*lw24 + 0.0001027543*lw15*lw25;
		c = -443565.465710869 + 718862.10804858*lw11 + 5114.1630934534*lw21 + -501588.206183097*lw12 + 784.4442123743*lw22 + 194712.394138513*lw13 + -377.9645994741*lw23 +
				-45088.7863182741*lw14 + 50.5960513287*lw24 + 6220.3310639927*lw15 + -2.3670776453*lw25 + -473.269785487*lw16 + 15.3213264134*lw17 + -3282.8510348085*lw11*lw21 +
				164.0438603974*lw12*lw22 + -5.2778800755*lw13*lw23 + 0.0929379305*lw14*lw24 + -0.0006826817*lw15*lw25;

		//logistic model
		prob = 1.0 / ( 1 + exp(-(b + c * d)) );

		if(prob >= .1) { //second logistic model
			b = 6770.9351457442 + -5933.7589419767*lw11 + -2976.2879986855*lw21 + 3279.7524685865*lw12 + 1602.7544794819*lw22 + -967.2906583423*lw13 + -462.0149190219*lw23 +
					159.8317289682*lw14 + 74.4884405822*lw24 + -14.0267151808*lw15 + -6.3644917671*lw25 + 0.5108811613*lw16 + 0.2252455343*lw26 + 0.965040193*lw12*lw22 +
					-0.0546309127*lw13*lw23 + 0.0012917084*lw14*lw24 + -1.14383e-05*lw15*lw25;
			c = 39406.5712626297 + -77863.1741143294*lw11 + 9586.8761567725*lw21 + 55360.1701572325*lw12 + -5825.2491611377*lw22 + -21887.8400068324*lw13 + 1751.6803621934*lw23 +
					5158.3764225203*lw14 + -290.1765894829*lw24 + -724.0348081819*lw15 + 25.364646181*lw25 + 56.0522105105*lw16 + -0.9172073892*lw26 + -1.8470088417*lw17 +
					449.4660736502*lw11*lw21 + -24.4141920625*lw12*lw22 + 0.8465834103*lw13*lw23 + -0.0158943762*lw14*lw24 + 0.0001235384*lw15*lw25;
			prob = 1.0 / ( 1 + exp(-(b + c * d)) );
			prob = prob < .1 ? .1 : prob;
		}

		//cout << prob << " " << d << " " << b << " " << c << endl;
		return prob;
	}

	/*
	double computeAbundanceCorrelation_new(const vector<float>& cov1, const vector<float>& cov2) {

		vector<float> means_f1 = cov1;
		vector<float> means_f2 = cov2;

		float mean_ratio = 0;
		for(size_t i=0; i<means_f1.size(); i++){
			if(means_f2[i] > 0){
				float ratio = means_f1[i] / means_f2[i];
				mean_ratio += ratio;
				cout << ratio << " ";
			}
			else{
				cout << "0" << " ";
			}
		}
		cout << endl;

		mean_ratio /= means_f1.size();

		for(size_t i=0; i<means_f1.size(); i++){
			means_f2[i] *= mean_ratio;
		}

		size_t i, ii;
		double sum_xsq = 0.0;
		double sum_ysq = 0.0;
		double sum_cross = 0.0;
		double ratio;
		double delta_x, delta_y;
		double mean_x = 0.0, mean_y = 0.0;
		double r = 0.0;

		size_t s = 0; //skipped

		for (i = 0; i < means_f1.size(); ++i) {
			double m1 = means_f1[i]; //ABD(r1,i);
			double m2 = means_f2[i]; //is_small ? small_ABD(r2, i) : ABD(r2, i);
			//m2 = rand()%500;

			//cout << m1 << " " << m2 << endl;
			ii = i - s;

			if(ii == 0) {
				mean_x = m1;
				mean_y = m2;
				continue;
			}

			ratio = ii / (ii + 1.0);
			delta_x = m1 - mean_x;
			delta_y = m2 - mean_y;
			sum_xsq += delta_x * delta_x * ratio;
			sum_ysq += delta_y * delta_y * ratio;
			sum_cross += delta_x * delta_y * ratio;
			mean_x += delta_x / (ii + 1.0);
			mean_y += delta_y / (ii + 1.0);
		}

		r = sum_cross / (sqrt(sum_xsq) * sqrt(sum_ysq));

		return r;
	}
	*/

	/*
	void computeUnitigComposition(IBank* inbank, size_t k, ModelCanonical::Iterator& itKmer, const vector<size_t>& kmer_to_compositionIndex, u_int64_t compositionVectorSize, vector<UnitigData>& unitigCompositions){


		//_datasetSequences = new vector<DatasetSequence*>();

		Iterator<Sequence>* itSeq = createIterator<Sequence> (
                                                          inbank->iterator(),
                                                          inbank->estimateNbItems(),
                                                          "Sorting sequences"
                                                          );

    	LOCAL (itSeq);
		{
			
			std::vector<Iterator<Sequence>*> itBanks =  itSeq->getComposition();

			for (size_t i=0; i<itBanks.size(); i++)
			{
				itSeq = createIterator<Sequence> (itBanks[i], inbank->estimateNbItemsBanki(i), "lala");

				for (itSeq->first(); !itSeq->isDone(); itSeq->next()){

					Sequence& sequence = itSeq->item();
					//if(sequence.getDataSize() < 100000) continue; //!!!!!!!!!!!!!!!!!!!!!!!!!

					string header = sequence.getComment();
					size_t pos = header.find("ctg");
					header.erase(pos, 3);
     				//str.replace(index, 3, "def");
					
					u_int32_t unitigName = stoull(header);
     				//if (index == std::string::npos) break;



					vector<float> composition;
					unitigToComposition(sequence, itKmer, kmer_to_compositionIndex, compositionVectorSize, composition);
					//cout << i << endl;
					unitigCompositions.push_back({i, unitigName, composition, {}, {}, sequence.getDataSize()});
					//DatasetSequence* datasetSequence = new DatasetSequence(sequence, i, 0, itKmer, kmer_to_compositionIndex, compositionVectorSize); //calc_sequence_score(sequence, k)
					//unitigIndex += 1;

				}
			}
			
		}

		cout << "Nb contigs: " << unitigCompositions.size() << endl;

	}
	*/

};


#endif 




