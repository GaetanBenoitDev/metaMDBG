

#ifndef MDBG_METAG_CONTIGFEATURE
#define MDBG_METAG_CONTIGFEATURE

#include "../Commons.hpp"

const float MU_INTRA = 0;
const float SIGMA_INTRA = 0.01037897 / 2.0;
const float MU_INTER = 0.0676654;
const float SIGMA_INTER = 0.03419337;

struct BinData{
	vector<float> _sequenceModel;
};

class ContigFeature{

public:

	size_t _kmerSize;
	KmerModel* _kmerModel;
	u_int64_t _compositionVectorSize;
	vector<size_t> _kmerToCompositionIndex;

	//vector<vector<float>> _comps;

	size_t _nbContigs;
	unordered_map<string, BinData> _bins;
	unordered_map<string, vector<float>> _contigCompositions;

	unordered_map<string, vector<u_int32_t>> _contigToNodenames;


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





	void computeFastaComposition(const string& sequenceFilename){


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

			if(_bins.find(name) != _bins.end()) continue;

			float maxProb = 1.0;
			string maxBinName = "";

			for(auto& it : _bins){

				const string& binName = it.first;
				const BinData& binData = it.second;

				
				//float prob_comp = computeCompositionProbability(binData._sequenceModel, composition);
				//float log_prob = -log10(prob_comp);
				float prob_comp = computeEuclideanDistance(binData._sequenceModel, composition);
				float log_prob = prob_comp;

				cout << "\t" << binName << ": " << log_prob << endl;
				if(log_prob < maxProb){
					maxProb = log_prob;
					maxBinName = binName;
				}
			}

			cout << name << " -> " << maxBinName << endl;

			for(u_int32_t nodeName : _contigToNodenames[name]){
				outfile << nodeName << "," << maxBinName << endl;
			}
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




