

#ifndef MDBG_METAG_TOBASESPACEONTHEFLY
#define MDBG_METAG_TOBASESPACEONTHEFLY

#include "../Commons.hpp"
#include "../utils/DnaBitset.hpp"

class ToBasespaceOnTheFly{

public:

	string _inputDir;
	unordered_map<u_int32_t, DnaBitset*> _nodeNameToSequence_left;
	unordered_map<u_int32_t, DnaBitset*> _nodeNameToSequence_right;

	ToBasespaceOnTheFly(){


	}

	void create(const string& inputDir){
		_inputDir = inputDir;

		cout << "Loading kminmer sequences" << endl;
		
		loadKminmerSequences();

		cout << "Nb kminmers sequences: " << _nodeNameToSequence_left.size() << " " << _nodeNameToSequence_right.size() << endl;
	}

	void loadKminmerSequences(){
		
		string filename_kminmerSequences = _inputDir + "/kminmerSequences";
		
		loadKminmerSequences_aux(filename_kminmerSequences + "_left.gz", _nodeNameToSequence_left);
		loadKminmerSequences_aux(filename_kminmerSequences + "_right.gz", _nodeNameToSequence_right);
	}

	void loadKminmerSequences_aux(const string& filename, unordered_map<u_int32_t, DnaBitset*>& nodeNameToSequence){
		
		gzFile file = gzopen(filename.c_str(), "rb");

		while(true){

			u_int32_t nodeName;
			//u_int64_t readIndex;
			u_int16_t sequenceLength;
			string sequence;
			gzread(file, (char*)&nodeName, sizeof(nodeName));
			
			if(gzeof(file)) break;

			//gzread(file, (char*)&readIndex, sizeof(readIndex));
			gzread(file, (char*)&sequenceLength, sizeof(sequenceLength));
			sequence.resize(sequenceLength);
			gzread(file, (char*)&sequence[0], sequenceLength);

			//cout << nodeName << " " << sequenceLength << " " << sequence << endl;

			nodeNameToSequence[nodeName] = new DnaBitset(sequence);
		}

		gzclose(file);
	}

	void createSequence(const vector<u_int32_t>& nodePath, string& contigSequence){

		//cout << nodePath.size() << endl;

		contigSequence = "";

		for(size_t i=0; i<nodePath.size(); i++){
			
			//cout << i << endl;

			u_int32_t nodeIndex = nodePath[i];
			bool orientation;
			u_int32_t nodeName = BiGraph::nodeIndex_to_nodeName(nodeIndex, orientation);
			//u_int64_t readIndex = supportingReads[i];
			//ContigNode contigNode = {nodeName, readIndex};

			if(orientation){
				
				//if(_nodeNameToSequence_right.find(nodeName) == _nodeNameToSequence_right.end()){
				//	cout << "pas normal" << endl;
				//}

				DnaBitset* dna = _nodeNameToSequence_right[nodeName];
				char* dnaStr = dna->to_string();
				string seq = string(dnaStr, dna->m_len);
				free(dnaStr);


				contigSequence += seq;
				
			}
			else{
				
				//if(_nodeNameToSequence_left.find(nodeName) == _nodeNameToSequence_left.end()){
				//	cout << "pas normal" << endl;
				//}

				DnaBitset* dna = _nodeNameToSequence_left[nodeName];
				char* dnaStr = dna->to_string();
				string seq = string(dnaStr, dna->m_len);
				free(dnaStr);

				Utils::revcomp(seq);
				contigSequence += seq;

			}

			/*
			if(i == 0){
				if(orientation){ //+
					cout << "Entire" << endl;

					contigSequence += correctedSequence;
				}
				else{
					cout << "Entire RC" << endl;

					Utils::revcomp(correctedSequence);
					contigSequence += correctedSequence;
				}
			}
			else {
				if(orientation){
					
					contigSequence += correctedSequence;
					

				}
				else{
					
					Utils::revcomp(correctedSequence);
					contigSequence += correctedSequence;

				}
			}
			*/
			
			
		}
	}



};


#endif 



