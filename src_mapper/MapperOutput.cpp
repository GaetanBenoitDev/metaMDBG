
#include "../src/Commons.hpp"
#include "../src/graph/GfaParser.hpp"

using namespace std;

void displayHelp(){
	cout << "Usage: ./simkaMin [option]" << endl;

	cout << endl << "[Distance computation options]" << endl;
	cout << "\tsketch      : transform datasets in small sketches of k-mers and their abundance" << endl;
	cout << "\tdistance    : compute Jaccard and Bray-Curtis distances between sketches" << endl;

	cout << endl << "[Distance matrix manipulation options]" << endl;
	cout << "\texport      : export distance matrices stored in binary format" << endl;
	//cout << "\tmatrix-update       : update existing distance matrices" << endl;

	cout << endl << "[Sketch options]" << endl;
	cout << "\tappend      : merge multiple sketch files into a single one" << endl;
	cout << "\tinfo        : list datasets contained in a sketch file" << endl;

	cout << endl;
}


class MappingIndex{

public:

	string _inputFilename_contigs;
	string _inputFilename_reads;
	unordered_map<string, u_int32_t> _contigName_to_contigIndex;
	unordered_map<string, u_int64_t> _readName_to_readIndex;
	//unordered_map<u_int64_t, vector<Alignment>> _alignments;
	//vector<string> _contigSequences;
	//vector<vector<vector<Window>>> _contigWindowSequences;
	//unordered_map<ContigRead, u_int32_t, ContigRead_hash> _alignmentCounts;

	MappingIndex(const string& contigFilename, const string& readFilename){
		_inputFilename_contigs = contigFilename;
		_inputFilename_reads = readFilename;
		indexContigName();
		indexReadName();
		cout << "done" << endl;
	}

	void indexContigName(){
		
		cout << "Indexing contig names" << endl;

		auto fp = std::bind(&MappingIndex::indexContigName_read, this, std::placeholders::_1);
		ReadParser readParser(_inputFilename_contigs, true, false);
		readParser.parse(fp);
	}

	void indexContigName_read(const Read& read){

		string minimap_name;

		auto find = read._header.find(' ');
		if(find == std::string::npos){
			minimap_name = read._header;
		}
		else{
			minimap_name = read._header.substr(0, find);
		}

		//cout << minimap_name << endl;
		_contigName_to_contigIndex[minimap_name] = read._index;
	}

	void indexReadName(){

		cout << "Indexing read names" << endl;
		cout << _inputFilename_reads << endl;
		auto fp = std::bind(&MappingIndex::indexReadName_read, this, std::placeholders::_1);
		ReadParser readParser(_inputFilename_reads, false, false);
		readParser.parse(fp);

	}

	void indexReadName_read(const Read& read){
		if(read._index % 100000 == 0) cout << read._index << endl;

		string minimap_name;

		auto find = read._header.find(' ');
		if(find == std::string::npos){
			minimap_name = read._header;
		}
		else{
			minimap_name = read._header.substr(0, find);
		}

		//string minimap_name = read._header.substr(0, read._header.find(' '));
		//cout << minimap_name << endl;
		_readName_to_readIndex[minimap_name] = read._index;
	}

};


int main (int argc, char* argv[])
{

	double errorThreshold = 0.3;

	//"/home/gats/workspace/run/overlap_test_AD/contigs_43_corrected.fasta.gz", "/home/gats/workspace/data/AD/HiFi/input.txt"
	string contigFilename = string(argv[1]);
	string readFilename = string(argv[2]);
	string outputFilename = string(argv[3]);
	cout << contigFilename << endl;
	cout << readFilename << endl;
	cout << outputFilename << endl;

	ofstream outputFile(outputFilename);
	vector<string>* fields = new vector<string>();
	vector<string>* fields_optional = new vector<string>();
	MappingIndex mappingIndex(contigFilename, readFilename);

    string lineInput;
	while (getline(cin,lineInput)) {
		//cout << lineInput << endl;
		//getchar();


		GfaParser::tokenize(lineInput, fields, '\t');

		//cout << line << endl;

		const string& readName = (*fields)[0];
		const string& contigName = (*fields)[5];

		u_int32_t readStart = stoull((*fields)[2]);
		u_int32_t readEnd = stoull((*fields)[3]);
		u_int64_t contigStart = stoull((*fields)[7]);
		u_int64_t contigEnd = stoull((*fields)[8]);

		u_int64_t nbMatches = stoull((*fields)[9]);
		u_int64_t alignLength = stoull((*fields)[10]);
		u_int64_t queryLength = stoull((*fields)[1]);

		bool strand = (*fields)[4] == "-";
		//float score = (double) nbMatches / (double) queryLength;
		//float score2 = (double) nbMatches / (double) alignLength;

		u_int32_t contigIndex = mappingIndex._contigName_to_contigIndex[contigName];
		u_int64_t readIndex = mappingIndex._readName_to_readIndex[readName];
		//u_int64_t length = std::max(readEnd - readStart, contigEnd - contigStart);


		//cout << contigIndex << " " << readIndex << endl;

		double length = max((double)(contigEnd - contigStart), (double)(readEnd - readStart));
		double error = 1 - min((double)(contigEnd - contigStart), (double)(readEnd - readStart)) / length;
		//cout << error << " " << errorThreshold << endl;
		
		if(error > errorThreshold) continue;

		float divergence = 0;
		//cout << error << endl;
		//getchar();
		//cout << contigIndex << " " << readIndex << endl;

		for(size_t i=12; i<fields->size(); i++){

			//cout << (*fields)[i] << endl;

			GfaParser::tokenize((*fields)[i], fields_optional, ':');

			if((*fields_optional)[0] == "dv"){
				divergence = std::stof((*fields_optional)[2]);

				/*
				//cout << (*fields_optional)[2] << endl;
				//cout << contigName << " " << readName << " " << (alignLength/queryLength*100) << " " << (divergence*100) << "     " << queryLength << " " << targetLength << endl;
				if(divergence < 0.05){
					string name = readName;
					size_t pos = name.find("ctg");
					name.erase(pos, 3);
					u_int32_t contigIndex = stoull(name);
					//cout << "Duplicate: " << contigIndex << endl;

					_duplicatedContigIndex.insert(contigIndex);
				}
				*/
			}

		}


		float score = divergence; //nbMatches / 

		outputFile.write((const char*)&contigIndex, sizeof(contigIndex));
		outputFile.write((const char*)&contigStart, sizeof(contigStart));
		outputFile.write((const char*)&contigEnd, sizeof(contigEnd));
		outputFile.write((const char*)&readIndex, sizeof(readIndex));
		outputFile.write((const char*)&readStart, sizeof(readStart));
		outputFile.write((const char*)&readEnd, sizeof(readEnd));
		outputFile.write((const char*)&strand, sizeof(strand));
		outputFile.write((const char*)&score, sizeof(score));
		//Alignment align = {contigIndex, readIndex, strand, readStart, readEnd, contigStart, contigEnd, score, length};

	}

	outputFile.close();

    return EXIT_SUCCESS;
}
