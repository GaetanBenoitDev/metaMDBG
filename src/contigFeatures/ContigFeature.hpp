

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
	unordered_map<string, vector<float>> _contigCompositions;

	u_int64_t _nbDatasets;
	unordered_map<string, vector<u_int32_t>> _contigToNodenames;

	unordered_map<u_int32_t, vector<u_int32_t>> _nodenameCounts;

	ContigFeature(){


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
	}

	unordered_map<string, vector<float>> _contigCoverages;
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
		cout << "W_intra: " << _w_intra << endl;
	}

	void loadAbundanceFile(const string& filename){

		ifstream infile(filename);

		vector<string>* fields = new vector<string>();
		string line;
		std::getline(infile, line); //skip header

		while (std::getline(infile, line)){
			
			GfaParser::tokenize(line, fields, '\t');

			string contigName = (*fields)[0];

			vector<float> coverages;
			for(size_t i=1; i<fields->size(); i++){
				const string& field = (*fields)[i];
				if(field.empty()) continue;
				//cout << i << " " << (*fields)[i] << endl;
				float ab = stof((*fields)[i]);
				if(ab == 0) ab = 1;
				coverages.push_back(ab);
			}

			_contigCoverages[contigName] = coverages;
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


	void sequenceToAbundance(const vector<u_int32_t>& sequence, vector<float>& abundances, vector<float>& abundancesVar){

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
				abundances[i] += counts[i];
				//values[i].push_back(counts[i]);
			}
		}

		
		for(u_int32_t i=0; i<abundances.size(); i++){
			//for(size_t j=0; j<sequence.size(); j++){
			//	cout << values[i][j] << " ";
			//}
			//cout << endl;
			//cout << abundances[i] << " " << sequence.size() << " " << (abundances[i] / sequence.size()) << endl;
			abundances[i] /= ((float) sequence.size());
			//abundances[i] = Utils::compute_median(values[i]);
		}

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


		//cout << sequence.size() << endl;

	}



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

			/*
			float prob_comp = computeCompositionProbability(_sequenceModel, composition);
			float log_prob = 0;

			float prob_product = prob_comp;

			if(prob_product > 0){
				log_prob = -log10(prob_comp); //- (math.log(prob_comp, 10) + math.log(prob_cov, 10));
			}
			else{
				//log_prob = MAX_WEIGHT
			}

			//log_prob_sum += log_prob

			if(prob_comp > 0.1){
				cout << header << " " << computeEuclideanDistance(_sequenceModel, composition) << " " <<  prob_comp << " " << log_prob << endl;
			}
			*/

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
		/*
		for(size_t i=0; i<_comps.size(); i++){			
			for(size_t j=i+1; j<_comps.size(); j++){
				cout << computeCompositionProbability(_comps[i], _comps[j]) << endl;
				exit(1);
			}
		}
		*/
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
		/*
		if(_sequenceModel.size() != 0){

			float prob_comp = computeCompositionProbability(_sequenceModel, composition);
			float log_prob = 0;

			float prob_product = prob_comp;

			if(prob_product > 0){
				log_prob = -log10(prob_comp); //- (math.log(prob_comp, 10) + math.log(prob_cov, 10));
			}
			else{
				//log_prob = MAX_WEIGHT
			}

			//log_prob_sum += log_prob

			if(prob_comp > 0.1){
				cout << header << " " << computeEuclideanDistance(_sequenceModel, composition) << " " <<  prob_comp << " " << log_prob << endl;
			}
		}
		*/
		//if( )
		//for(float v : composition){
		//	cout << v << " ";
		//}
		//cout << endl;

		//exit(1);
		//_comps.push_back(composition);
		_nbContigs += 1;
	}

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

	bool isIntra(const ContigFeatures& f1, const ContigFeatures& f2){
		float prob = computeProbability(f1 , f2);
		//cout << "\tProb: " << prob << " " << _w_intra << endl;
		return prob < _w_intra;
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
				m1 = std::max(m1, (Distance) 1e-6);
				m2 = std::max(m2, (Distance) 1e-6);
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




