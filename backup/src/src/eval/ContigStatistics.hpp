

#ifndef MDBG_METAG_CONTIGSTATISTICS
#define MDBG_METAG_CONTIGSTATISTICS

#include "Commons.hpp"

class ContigStatistics{

public:

	string _outputDir;
	const vector<u_int32_t>& _nodeAbundances;
	const vector<u_int32_t>& _contigNodes;
	u_int32_t _contigIndex;

	ContigStatistics(const string& outputDir, const vector<u_int32_t>& nodeAbundances, const vector<u_int32_t>& contigNodes, u_int32_t contigIndex): _nodeAbundances(nodeAbundances), _contigNodes(contigNodes)  {
       
		_outputDir = outputDir + "/eval";

	    fs::path path(_outputDir);
	    if(!fs::exists (path)){
            fs::create_directory(path);
        }

		_contigIndex = contigIndex;
	}

	void execute(){
		
		ofstream outputFile_1(_outputDir + "/" + to_string(_contigIndex) + "_abundances_1" + ".txt");

		for(u_int32_t nodeIndex : _contigNodes){
			u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex);
			u_int32_t nodeAbundance = _nodeAbundances[nodeName];
			outputFile_1 << nodeAbundance << " ";

		}

		outputFile_1.close();

		for(int t : {10, 100, 1000}){

			ofstream outputFile(_outputDir + "/" + to_string(_contigIndex) + "_abundances_" + to_string(t) + ".txt");

			int contigPos = 0;

			for(u_int32_t nodeIndex : _contigNodes){

				vector<u_int32_t> abundances;


				for(int i=-t/2; i<t/2; i++){
					int pos = contigPos + i;
					if(pos < 0){
						pos = _contigNodes.size() + pos;
					}
					else if(pos >= _contigNodes.size()){
						pos = pos - _contigNodes.size();
					}

					//cout << i << " " << pos << " " << _contigNodes.size() << endl;

					u_int32_t nodeIndex2 = _contigNodes[pos];
					u_int32_t nodeName2 = BiGraph::nodeIndex_to_nodeName(nodeIndex2);
					u_int32_t nodeAbundance2 = _nodeAbundances[nodeName2];
					abundances.push_back(nodeAbundance2);
				}

				u_int32_t median = Utils::compute_median(abundances);
				outputFile << median << " ";

				contigPos += 1;
			}
			
			outputFile.close();
		}


	}

};

#endif