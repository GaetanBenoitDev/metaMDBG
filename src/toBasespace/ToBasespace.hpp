
//Real data
//./bin/Bloocoo  -i ~/workspace/data/overlap_test/input.txt -l 21 -k 3 -d 0.005 -o ~/workspace/run/overlap_test/ -ihifiasm ~/workspace/run/hifiasm_meta/AD_components/big/component_3.fasta -idir ~/workspace/run/overlap_test/

//Simulation
///bin/Bloocoo  -i ~/workspace/data/overlap_test/input.txt -l 16 -k 3 -d 0.005 -o ~/workspace/run/overlap_test_3/ -ihifiasm ~/workspace/data/overlap_test/genome_2371_20x/truth_input.txt -idir ~/workspace/run/overlap_test_3/

#ifndef MDBG_METAG_TOBASESPACE
#define MDBG_METAG_TOBASESPACE

#include "../Commons.hpp"




class ToBasespace : public Tool{
    
public:

	string _inputContigFilename;
	string _outputContigFilename;

    ToBasespace ();
    void execute ();


};	


#endif 



