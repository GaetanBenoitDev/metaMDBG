
#ifndef MDBG_METAG_MINIMIZERCHAINER
#define MDBG_METAG_MINIMIZERCHAINER

#include "../Commons.hpp"

struct AlignmentResult2{

	public:

	
	ReadType _referenceReadIndex;
	ReadType _queryReadIndex;
	//u_int32_t _cigarReferenceStart;
	//u_int32_t _cigarQueryStart;
	//string _cigar;
	bool _isQueryReversed;
	float _chainingScore;
	int32_t _nbMatches;
	int32_t _nbMissmatches;
	int32_t _nbDeletions;
	int32_t _nbInsertions;
	float _identity;
	u_int32_t _overHangStart;
	u_int32_t _overHangEnd;
	u_int32_t _alignLength;
	vector<pair<u_int32_t, u_int32_t>> _alignments;
	//float _divergence;
	u_int32_t _referenceStart;
	u_int32_t _referenceEnd;
	u_int32_t _queryStart;
	u_int32_t _queryEnd;
	int64_t _queryLength;
	int64_t _referenceLength;

	AlignmentResult2(){

		_referenceReadIndex = 0;
		_queryReadIndex = 0;
		_chainingScore = 0;
		_nbMatches = 0;
		_nbMissmatches = 0;
		_nbDeletions = 0;
		_nbInsertions = 0;
		_identity = 0;
		_overHangStart = 0;
		_overHangEnd = 0;
		_alignLength = 0;
		//_divergence = 1;
		_referenceStart = 0;
		_referenceEnd = 0;
		_queryStart = 0;
		_queryEnd = 0;
		_queryLength = 0;
		_referenceLength = 0;
	}

	bool isMaximalMapping(const int64_t maxOverhang){

		const int64_t queryStart = _queryStart;
		const int64_t queryEnd = _queryEnd;
		const int64_t referenceStart = _referenceStart;
		const int64_t referenceEnd = _referenceEnd;

		if (((queryStart < maxOverhang) || (referenceStart < maxOverhang)) && ((queryEnd+maxOverhang > _queryLength) || (referenceEnd+maxOverhang > _referenceLength))) return true;
		return false;
		
	}

};

struct Anchor{

	int32_t _referencePosition;
	int32_t _queryPosition;
	bool _isReversed;
	int32_t _referencePositionIndex;
	int32_t _queryPositionIndex;

	Anchor(const int32_t referencePosition, const int32_t queryPosition, bool isReversed, int32_t referencePositionIndex, int32_t queryPositionIndex){
		_referencePosition = referencePosition;
		_queryPosition = queryPosition;
		_isReversed = isReversed;
		_referencePositionIndex = referencePositionIndex;
		_queryPositionIndex = queryPositionIndex;
	}
};

struct Chain{
	float _chainingScore;
	vector<size_t> _anchorInterval;
};



class MinimizerChainer {
    
public:





	bool _print_debug;
	size_t _minimizerSize;

	MinimizerChainer(size_t minimizerSize){
		_print_debug = false;
		_minimizerSize = minimizerSize;
	}

	AlignmentResult2 computeChainingAlignment(vector<Anchor>& anchors, const MinimizerRead& referenceRead, const MinimizerRead& queryRead, const int64_t& maxChainingBand){
		
		//size_t nbAnchors = _anchorIndex;
		//cout << "nb anchors: " << nbAnchors << endl;
		//return computeAlignment(anchors, referenceRead, queryRead, minimizerAligner);
		//_print_debug = queryRead._readIndex == 231;


		AlignmentResult2 alignmentResult;
		/*
		vector<Anchor> anchors;

		for(size_t i=0; i<queryRead._minimizers.size(); i++){
			
			MinimizerType queryMinimizer = queryRead._minimizers[i];

			if(referenceReadMinimizerPositionMap.find(queryMinimizer) == referenceReadMinimizerPositionMap.end()) continue;

			u_int32_t queryPosition = queryRead._minimizersPos[i];
			bool queryIsReversed = queryRead._readMinimizerDirections[i];

			const vector<MinimizerPosition>& referenceMinimizerPositions = referenceReadMinimizerPositionMap[queryMinimizer];
			
			for(const MinimizerPosition& referenceMinimizerPosition : referenceMinimizerPositions){
				anchors.push_back({referenceMinimizerPosition._position, queryPosition, referenceMinimizerPosition._isReversed != queryIsReversed, referenceMinimizerPosition._positionIndex, i});
			}
			//MinimizerPosition queryMinmizerPosition = {position, i, isReversed};
			//referenceReadMinimizerPosition[minimizer].push_back(minmizerPosition);
			
		}
		*/
		
		if(anchors.size() < 3) return alignmentResult;
		//if(anchors.size() > 150) continue;


		//const MinimizerRead& referenceReadHighDensity = read;//_parent._readWorkers[referenceReadIndex]._referenceRead;
		//const MinimizerRead& referenceReadLowDensity = //Utils::getLowDensityMinimizerRead(referenceReadHighDensity, _parent._minimizerDensityAssembly);


		std::sort(anchors.begin(), anchors.end(), [](const Anchor& a, const Anchor& b){
			if(a._referencePosition == b._referencePosition){
				return a._queryPosition < b._queryPosition;
			}
			return a._referencePosition < b._referencePosition;
		});
		
		if(_print_debug){
			cout << endl << endl << endl << endl << endl;
			cout << "Reference read: " << referenceRead._readIndex << endl;
			for(size_t i=0; i<referenceRead._minimizers.size(); i++){
				cout << "\t" << i << "\t" << referenceRead._minimizers[i] << endl;
			}
			cout << "Query read: " << queryRead._readIndex << endl;
			for(size_t i=0; i<queryRead._minimizers.size(); i++){
				cout << "\t" << i << "\t" << queryRead._minimizers[i] << endl;
			}
			cout << endl;

			cout << "Anchors: " << endl;
			for(size_t i=0; i<anchors.size(); i++){
				const Anchor& anchor = anchors[i];
				cout << "\t" << i << "\t" << anchor._referencePosition << "\t" << anchor._queryPosition << "\t" << anchor._referencePositionIndex << "\t" << anchor._queryPositionIndex << "\t" << anchor._isReversed << "\t" << referenceRead._minimizers[anchor._referencePositionIndex] << endl;
			}

		}
		




		bool printChaining = false; //_print_debug;

		//if(anchors.size() < 2) return alignmentResult;



		
		/*
		unordered_map<MinimizerType, vector<MinimizerPosition>> minimizer_to_referencePositions;

		for(u_int32_t i=0; i<referenceRead._minimizers.size(); i++){

			MinimizerType minimizer = referenceRead._minimizers[i];
			u_int32_t position = referenceRead._minimizersPos[i];
			bool isReversed = referenceRead._readMinimizerDirections[i];

			MinimizerPosition minmizerPosition = {referenceRead._readIndex, position, i, isReversed};
			minimizer_to_referencePositions[minimizer].push_back(minmizerPosition);
		}

		//u_int64_t queryReadLength = queryRead._minimizersPos[queryRead._minimizersPos.size()-1];

		vector<Anchor> anchors;

		for(u_int32_t i=0; i<queryRead._minimizers.size(); i++){

			MinimizerType minimizer = queryRead._minimizers[i];
			u_int32_t queryMinimizerPosition = queryRead._minimizersPos[i];
			bool queryMinimizerIsReversed = queryRead._readMinimizerDirections[i];
			//u_int32_t queryPosition = position;
			
			//if(isQueryReversed){
			//	queryPosition = queryReadLength - position;
			//}

			if(minimizer_to_referencePositions.find(minimizer) == minimizer_to_referencePositions.end()) continue;

			for(const MinimizerPosition& referenceMinimizerPosition : minimizer_to_referencePositions[minimizer]){
				anchors.push_back({referenceMinimizerPosition._position, queryMinimizerPosition, referenceMinimizerPosition._isReversed != queryMinimizerIsReversed, referenceMinimizerPosition._positionIndex, i});
			}
			
		}

		std::sort(anchors.begin(), anchors.end(), [](const Anchor& a, const Anchor& b){
			if(a._referencePosition == b._referencePosition){
				return a._queryPosition < b._queryPosition;
			}
			return a._referencePosition < b._referencePosition;
		});
		*/

		//vector<Point2> F;
		vector<Chain> chainIntervals = chainAnchors(anchors, 0.01, 0.5, referenceRead, queryRead, maxChainingBand);
		/*
		//std::reverse(F.begin(), F.end());
		//std::reverse(maxNode.begin(), maxNode.end());

		if(printChaining){
			cout << "Chain score: " << chainScore << endl;

			for(size_t i=0; i<F.size(); i++){
				cout << "\t" << F[i]._fromIndex << "\t" << F[i]._score << "\t" << F[i]._toIndex << endl;
			}
			cout << endl;
			for(float y : maxNode){
				cout << "\t" << y << endl;
			}


			cout << endl;
			for(float y : maxNode){
				cout << "\t" << anchors[y]._referencePositionIndex << " " << anchors[y]._queryPositionIndex << " " << anchors[y]._referencePosition << " " << anchors[y]._queryPosition << " " << referenceRead._minimizers[anchors[y]._referencePositionIndex] << " " << queryRead._minimizers[anchors[y]._queryPositionIndex] << endl;
			}

		}

		if(maxNode.size() <= 2) return alignmentResult;
		*/

		if(chainIntervals.size() == 0) return alignmentResult;

		const Chain& chain = chainIntervals[0];
		const vector<size_t>& chainInterval = chain._anchorInterval;
		if(chainInterval.size() <= 3) return alignmentResult;

		const Anchor& firstAnchor = anchors[chainInterval[0]];
		const Anchor& secondAnchor = anchors[chainInterval[1]];
		const Anchor& lastAnchor = anchors[chainInterval[chainInterval.size()-1]];

		bool isQueryReversed = firstAnchor._queryPositionIndex > lastAnchor._queryPositionIndex;
		//if(!isQueryReversed) return alignmentResult;


		if(printChaining){
			cout << "\tIs reversed: " << isQueryReversed << endl;
		}


		int64_t referenceLength = referenceRead._readLength;// ._minimizersPos[referenceRead._minimizersPos.size()-1];
		int64_t queryLength = queryRead._readLength;// ._minimizersPos[queryRead._minimizersPos.size()-1];

		
		int32_t nbStartingMissmatches = 0;
		int64_t overHangStart = 0;
		if(isQueryReversed){
			overHangStart = min((int64_t)referenceRead._minimizersPos[firstAnchor._referencePositionIndex], (int64_t)queryLength - queryRead._minimizersPos[firstAnchor._queryPositionIndex-1]);
			nbStartingMissmatches = min(firstAnchor._referencePositionIndex, (int32_t) queryRead._minimizers.size()-firstAnchor._queryPositionIndex-1);
		}
		else{
			overHangStart = min((int64_t)referenceRead._minimizersPos[firstAnchor._referencePositionIndex], (int64_t) queryRead._minimizersPos[firstAnchor._queryPositionIndex]);
			nbStartingMissmatches = min(firstAnchor._referencePositionIndex, firstAnchor._queryPositionIndex);
		}

		int32_t nbEndingMissmatches = 0;
		int64_t overHangEnd = 0;
		if(isQueryReversed){
			overHangEnd = min((int64_t) referenceLength - referenceRead._minimizersPos[lastAnchor._referencePositionIndex-1], (int64_t)queryRead._minimizersPos[lastAnchor._queryPositionIndex]);
			nbEndingMissmatches = min((int32_t) referenceRead._minimizers.size()-lastAnchor._referencePositionIndex-1, lastAnchor._queryPositionIndex);
		}
		else{
			overHangEnd = min((int64_t) referenceLength - referenceRead._minimizersPos[lastAnchor._referencePositionIndex-1], (int64_t) queryLength - queryRead._minimizersPos[lastAnchor._queryPositionIndex-1]);
			nbEndingMissmatches = min(referenceRead._minimizers.size()-lastAnchor._referencePositionIndex-1, queryRead._minimizers.size()-lastAnchor._queryPositionIndex-1);
		}


		//nbStartingMissmatches = 0;
		//nbEndingMissmatches = 0;
		
		if(printChaining) cout << "Nb starting missmatches: " << nbStartingMissmatches << endl;
		if(printChaining) cout << "Nb ending missmatches: " << nbEndingMissmatches << endl;
		//cout << "Nb starting missmatches: " << nbStartingMissmatches << endl;
		//cout << "Nb ending missmatches: " << nbEndingMissmatches << endl;

		int64_t nbMissmatchesMiddle = 0;
		//int64_t totalNbMatches = 0;
		//int64_t totalNbMissmatches = 0;
		//int64_t totalNbInsertions = 0;

		//string cigar = "";




		//cout << nbStartingMissmatches << " " << nbEndingMissmatches << " " << overHangStart << " " << overHangEnd << endl;
		//getchar();
		u_int32_t alignStart = -1;
		u_int32_t alignEnd = -1;

		u_int32_t referencePosition = anchors[chainInterval[0]]._referencePositionIndex - nbStartingMissmatches;
		u_int32_t queryPosition = anchors[chainInterval[0]]._queryPositionIndex;

		if(isQueryReversed){
			queryPosition += nbStartingMissmatches;
		}
		else{
			queryPosition -= nbStartingMissmatches;
		}

		for(int i=0; i<nbStartingMissmatches; i++){
			//cigar += "M";
			alignmentResult._alignments.push_back({referencePosition, queryPosition});
			alignmentResult._nbMissmatches += 1;
			referencePosition += 1;
			if(isQueryReversed){
				queryPosition -= 1;
			}
			else{
				queryPosition += 1;
			}
		}

		for(size_t i=0; i<chainInterval.size()-1; i++){

			//size_t currentAnchorIndex = maxNode[i];
			//size_t nextAnchorIndexTo = maxNode[i+1];

			//if(printChaining) cout << currentAnchorIndex << " " << nextAnchorIndexTo << endl;
			const Anchor& currentAnchor = anchors[chainInterval[i]];
			const Anchor& nextAnchor = anchors[chainInterval[i+1]];
			//cout << "LOL: " << chainInterval[i] << " " << currentAnchor._referencePositionIndex << " " << nextAnchor._referencePositionIndex << endl;
			//cout << "\t" << i << endl;

			if(printChaining){
				cout << "\t" << i << "\t" << chainInterval[i] << "\t" << currentAnchor._referencePositionIndex << "\t" << currentAnchor._queryPositionIndex << "\t" << currentAnchor._referencePosition << "\t" << currentAnchor._queryPosition << "\t" << referenceRead._minimizers[currentAnchor._referencePositionIndex] << "\t" << queryRead._minimizers[currentAnchor._queryPositionIndex] << "\t" << currentAnchor._isReversed << "\t" << nextAnchor._isReversed << endl;
				cout << "\t" << i+1 << "\t" << chainInterval[i+1] << "\t" << nextAnchor._referencePositionIndex << "\t" << nextAnchor._queryPositionIndex << "\t" << nextAnchor._referencePosition << "\t" << nextAnchor._queryPosition << "\t" << referenceRead._minimizers[currentAnchor._referencePositionIndex] << "\t" << queryRead._minimizers[currentAnchor._queryPositionIndex] << "\t" << currentAnchor._isReversed << "\t" << nextAnchor._isReversed << endl;
			}


			int referenceGapIndexSize = nextAnchor._referencePositionIndex - currentAnchor._referencePositionIndex - 1;
			int queryGapIndexSize = 0;//
			
			if(isQueryReversed){
				queryGapIndexSize = currentAnchor._queryPositionIndex - nextAnchor._queryPositionIndex - 1;
			}
			else{
				queryGapIndexSize = nextAnchor._queryPositionIndex - currentAnchor._queryPositionIndex - 1;
			}

			
			//if(printChaining) cout << "\t" << referenceGapIndexSize << " " << queryGapIndexSize << endl;
			
			int nbMissmatches = min(referenceGapIndexSize, queryGapIndexSize);
			int nbInsertions = 0;
			int nbDeletions = 0;
			
			if(referenceGapIndexSize > queryGapIndexSize){
				nbDeletions = referenceGapIndexSize - nbMissmatches;
			}
			else{
				nbInsertions = queryGapIndexSize - nbMissmatches;
			}
			

			//if(printChaining) cout << "\t" << nbMissmatches << " " << nbInsertions << " " << nbDeletions << endl;


			//if(nbMissmatches < 0 || nbInsertions < 0 || nbDeletions < 0){
			//	cout << "issue" << endl;
			//	continue;
			//}

			if(printChaining) cout << "\tMatch        : " << referenceRead._minimizers[currentAnchor._referencePositionIndex] << endl;
			//cigar += "M";
			alignmentResult._alignments.push_back({referencePosition, queryPosition});
			referencePosition += 1;
			if(isQueryReversed){
				queryPosition -= 1;
			}
			else{
				queryPosition += 1;
			}
			alignmentResult._nbMatches += 1;

			if(alignStart == -1) alignStart = currentAnchor._referencePosition;
			alignEnd = nextAnchor._referencePosition;

			//int64_t nbGapsReference = nbMissmatches + nbDeletions;
			//int64_t nbGapsQuery = nbMissmatches + nbInsertions;

			alignmentResult._nbMissmatches += nbMissmatches;
			nbMissmatchesMiddle += nbMissmatches;
			alignmentResult._nbDeletions += nbDeletions;
			alignmentResult._nbInsertions += nbInsertions;

			for(size_t i=0; i<nbMissmatches; i++){
				if(printChaining)cout << "\tMissmatch reference: " << referenceRead._minimizers[currentAnchor._referencePositionIndex+i+1] << endl;// << " " << queryRead._minimizers[currentAnchor._queryPositionIndex+i+1] << endl;
				//cigar += "M";
				alignmentResult._alignments.push_back({referencePosition, (u_int32_t)-1});
				referencePosition += 1;
			}

			for(size_t i=0; i<nbDeletions; i++){
				if(printChaining) cout << "\tDeletion : " << referenceRead._minimizers[currentAnchor._referencePositionIndex+i+1] << endl;
				//cigar += "D";
				alignmentResult._alignments.push_back({referencePosition, (u_int32_t)-1});
				referencePosition += 1;
			}

			/*
			for(size_t i=0; i<nbGapsReference; i++){
				if(printChaining) cout << "\tDeletion: " << referenceRead._minimizers[currentAnchor._referencePositionIndex+i+1] << endl;
				//cigar += "D";
				alignmentResult._alignments.push_back({referencePosition, (u_int32_t)-1});
				//referencePosition += 1;
			}
			*/

			/*
			for(size_t i=0; i<nbGapsQuery; i++){
				if(printChaining){
					if(isQueryReversed){
						cout << "\tInsertion: " << endl;//<< queryRead._minimizers[queryRead._minimizers.size()-currentAnchor._queryPositionIndex-i-1] << endl;
					}
					else{
						cout << "\tInsertion: " << endl; //<< queryRead._minimizers[currentAnchor._queryPositionIndex+i+1] << endl;
					}
				}
				alignmentResult._alignments.push_back({(u_int32_t)-1, queryPosition});
				if(isQueryReversed){
					queryPosition -= 1;
				}
				else{
					queryPosition += 1;
				}
				//cigar += "I";
			}
			*/

			for(size_t i=0; i<nbMissmatches; i++){
				if(printChaining)cout << "\tMissmatch query: " << endl;// << " " << queryRead._minimizers[currentAnchor._queryPositionIndex+i+1] << endl;
				//cigar += "M";
				alignmentResult._alignments.push_back({(u_int32_t)-1, queryPosition});
				if(isQueryReversed){
					queryPosition -= 1;
				}
				else{
					queryPosition += 1;
				}
			}

			for(size_t i=0; i<nbInsertions; i++){
				if(printChaining){
					if(isQueryReversed){
						cout << "\tInsertion: " << endl;//<< queryRead._minimizers[queryRead._minimizers.size()-currentAnchor._queryPositionIndex-i-1] << endl;
					}
					else{
						cout << "\tInsertion: " << endl; //<< queryRead._minimizers[currentAnchor._queryPositionIndex+i+1] << endl;
					}
				}
				alignmentResult._alignments.push_back({(u_int32_t)-1, queryPosition});
				if(isQueryReversed){
					queryPosition -= 1;
				}
				else{
					queryPosition += 1;
				}
			}

			/*
			for(size_t i=0; i<nbDeletions; i++){
				if(printChaining) cout << "\tDeletion : " << referenceRead._minimizers[currentAnchor._referencePositionIndex+i+1] << endl;
				//cigar += "D";
				alignmentResult._alignments.push_back({referencePosition, (u_int32_t)-1});
				referencePosition += 1;
				alignmentResult._nbDeletions += 1;
			}
			for(size_t i=0; i<nbInsertions; i++){
				if(printChaining){
					if(isQueryReversed){
						cout << "\tInsertion: " << endl;//<< queryRead._minimizers[queryRead._minimizers.size()-currentAnchor._queryPositionIndex-i-1] << endl;
					}
					else{
						cout << "\tInsertion: " << endl; //<< queryRead._minimizers[currentAnchor._queryPositionIndex+i+1] << endl;
					}
				}
				alignmentResult._alignments.push_back({(u_int32_t)-1, queryPosition});
				if(isQueryReversed){
					queryPosition -= 1;
				}
				else{
					queryPosition += 1;
				}
				//cigar += "I";
				alignmentResult._nbInsertions += 1;
			}
			for(size_t i=0; i<nbMissmatches; i++){
				if(printChaining)cout << "\tMissmatch: " << referenceRead._minimizers[currentAnchor._referencePositionIndex+i+1] << endl;// << " " << queryRead._minimizers[currentAnchor._queryPositionIndex+i+1] << endl;
				//cigar += "M";
				alignmentResult._alignments.push_back({referencePosition, queryPosition});
				referencePosition += 1;
				if(isQueryReversed){
					queryPosition -= 1;
				}
				else{
					queryPosition += 1;
				}
				alignmentResult._nbMissmatches += 1;
				nbMissmatchesMiddle += 1;
			}
			*/
		}

		alignmentResult._alignments.push_back({referencePosition, queryPosition});
		referencePosition += 1;
		if(isQueryReversed){
			queryPosition -= 1;
		}
		else{
			queryPosition += 1;
		}
		//cigar += "M";
		alignmentResult._nbMatches += 1;

		for(int i=0; i<nbEndingMissmatches; i++){
			alignmentResult._alignments.push_back({referencePosition, queryPosition});
			referencePosition += 1;
			if(isQueryReversed){
				queryPosition -= 1;
			}
			else{
				queryPosition += 1;
			}
			//cigar += "M";
			alignmentResult._nbMissmatches += 1;
		}

		//cout << totalNbMatches << " " << totalNbMissmatches << endl;
		//if(queryRead._readIndex == 231) getchar();

		//if(isMainAlignment){

			//cout << "score check disabled" << endl;
			//if(_parent._print_debug){

			//}
			//else{
		u_int64_t referenceSize = alignmentResult._nbMatches + alignmentResult._nbMissmatches + alignmentResult._nbDeletions;
		u_int64_t querySize = alignmentResult._nbMatches + alignmentResult._nbMissmatches + alignmentResult._nbInsertions;

		//u_int64_t querySize = 0;
		//if(isQueryReversed){
		//	querySize = firstAnchor._queryPositionIndex - lastAnchor._queryPositionIndex + 1;
		//}
		//else{
		//	querySize = lastAnchor._queryPositionIndex - firstAnchor._queryPositionIndex + 1;
		//}

		//cout << querySize << " " << (alignmentResult._nbMatches + nbMissmatchesMiddle + alignmentResult._nbInsertions) << endl;
		//if(referenceSize > querySize){
		//	if(alignmentResult._nbMatches - alignmentResult._nbMissmatches - alignmentResult._nbDeletions < 5) return alignmentResult;
		//}
		//else{
		//	if(alignmentResult._nbMatches - alignmentResult._nbMissmatches - alignmentResult._nbInsertions < 5) return alignmentResult;
		//}
		
		//if(alignmentResult._nbMatches < 5) return alignmentResult;
		//if(overHangStart > 2000) return alignmentResult;
		//if(overHangEnd > 2000) return alignmentResult;
		//cout << alignmentResult._nbMatches << " " << alignmentResult._nbMissmatches << endl;
		//if(alignmentResult._nbMatches - alignmentResult._nbMissmatches < 5) return alignmentResult;

		//float divergence = computeDivergence();
		//if(divergence > _parent._maxDivergence) return alignmentResult;
		//if(chain._chainingScore < 40) return alignmentResult;
		//}
		//}

		
	
		long double nbSeeds = min(referenceSize, querySize);
		float divergence = 0;

		//cout << nbMatches << " " << nbSeeds << endl;
		//double nbSeeds = alignmentResult._nbMatches + nbMissmatchesMiddle + alignmentResult._nbInsertions + alignmentResult._nbDeletions;

		if(_print_debug) cout << "Divergence: " << alignmentResult._nbMatches << " " << querySize << endl;
		if(_print_debug) cout << alignmentResult._nbMatches << "\t" << nbMissmatchesMiddle << "\t" << alignmentResult._nbDeletions << "\t" << alignmentResult._nbInsertions << endl;
		if(alignmentResult._nbMatches == nbSeeds){
			divergence = 0;
		}
		else if(alignmentResult._nbMatches == 0){
			divergence = 1;
		}
		else{
			divergence = 1.0 - pow((alignmentResult._nbMatches / nbSeeds), 1.0/_minimizerSize);
		}
		
		//if(overHangStart > 2000 || overHangEnd > 2000){
		//	divergence = 1;
		//}
		//if(divergence > 0.04) return alignmentResult;
		//if(alignmentResult._nbMatches < 20) return alignmentResult;

		//float divergence = computeDivergence(anchors, chain, referenceReadLowDensity, referenceReadHighDensity, queryReadLowDensity, queryReadHighDensity, isQueryReversed);


		//cout << cigar << endl;

		//string cigarCompressed = compressCigar(cigar);

		alignmentResult._queryLength = queryLength;
		alignmentResult._referenceLength = referenceLength;
		
		alignmentResult._referenceReadIndex = referenceRead._readIndex;
		alignmentResult._queryReadIndex = queryRead._readIndex;
		//alignmentResult._cigarReferenceStart = cigarReferenceStart;
		//alignmentResult._cigarQueryStart = cigarQueryStart;
		//alignmentResult._cigar = cigarCompressed;
		alignmentResult._isQueryReversed = isQueryReversed;
		alignmentResult._chainingScore = chain._chainingScore;
		alignmentResult._identity = 1.0-divergence;
		//alignmentResult._divergence = divergence;
		alignmentResult._overHangStart = overHangStart;
		alignmentResult._overHangEnd = overHangEnd;
		alignmentResult._alignLength = alignEnd - alignStart;
		alignmentResult._referenceStart = firstAnchor._referencePosition;
		alignmentResult._referenceEnd = lastAnchor._referencePosition;
		if(isQueryReversed){
			//alignmentResult._queryStart = queryRead._readLength - anchors[0]._queryPosition;
			//alignmentResult._queryEnd = queryRead._readLength - anchors[anchors.size()-1]._queryPosition;
			alignmentResult._queryStart = lastAnchor._queryPosition;
			alignmentResult._queryEnd = firstAnchor._queryPosition;
		}
		else{
			alignmentResult._queryStart = firstAnchor._queryPosition;
			alignmentResult._queryEnd = lastAnchor._queryPosition;
		}
		//cout << referenceRead._readIndex << " " << queryRead._readIndex << " " << (1.0-divergence) << " " << overHangStart << " " << overHangEnd << " " << (alignEnd - alignStart) << endl;
		//getchar();
		
		//#pragma omp critical
		//{
		//	cout << alignmentResult._alignLength << endl;
		//}
		/*
		//if(isMainAlignment){
		unordered_map<u_int32_t, u_int32_t> highDensityReferenceIndex;
		unordered_map<u_int32_t, u_int32_t> highDensityQueryIndex;

		for(size_t i=0; i<referenceReadHighDensity._minimizers.size(); i++){
			highDensityReferenceIndex[referenceReadHighDensity._minimizersPos[i]] = i;
		}
		for(size_t i=0; i<queryReadHighDensity._minimizers.size(); i++){
			highDensityQueryIndex[queryReadHighDensity._minimizersPos[i]] = i;
		}

		
		float divergence = computeDivergence(alignmentResult, highDensityReferenceIndex, highDensityQueryIndex, referenceReadLowDensity, referenceReadHighDensity, queryReadLowDensity, queryReadHighDensity);
		if(_parent._print_debug) cout << "\tDivergence: " << divergence << endl;

		if(divergence > _parent._maxDivergence){
			//alignmentResult._cigar.clear(); //return null alignment
		}
		*/

		normalizeAlignment(alignmentResult._alignments, referenceRead._minimizers, queryRead._minimizers);

		return alignmentResult;
	}

	/*
	u_int64_t _anchorIndex;
	vector<float> _scores;
	vector<float> _parents;
	vector<Anchor> _anchors;

	void addAnchor(const int32_t& referencePosition, const int32_t& queryPosition, const bool& isReversed, const int32_t& referencePositionIndex, const int32_t& queryPositionIndex){
		
		//cout << "add anchor " << _anchorIndex << " " << _anchors.size() << endl;
		if(_anchorIndex >= _anchors.size()){
			//cout << "\tpush" << endl;
			_anchors.push_back({referencePosition, queryPosition, isReversed, referencePositionIndex, queryPositionIndex});
		}
		else{
			//cout << "\tset" << endl;
			Anchor& anchor = _anchors[_anchorIndex];
			anchor._referencePosition = referencePosition;
			anchor._queryPosition = queryPosition;
			anchor._isReversed = isReversed;
			anchor._referencePositionIndex = referencePositionIndex;
			anchor._queryPositionIndex = queryPositionIndex;
		}
		_anchorIndex += 1;

		//cout << "add anchor done" << endl;
	}
	*/

	vector<Chain> chainAnchors(const vector<Anchor>& anchors, float a, float b, const MinimizerRead& referenceRead, const MinimizerRead& queryRead, const int64_t& maxChainingBand) {

		//size_t nbAnchors = _anchorIndex;
		vector<Chain> chainIntervals;

		//cout << "---- " << points.size() << endl;
		float w = 20;//_minimizerSize;// 20;// _minimizerSize; //_minimizerSize*5; //_parent._minimizerSize;
		vector<float> scores(anchors.size(), 0);
		vector<size_t> parents(anchors.size(), 0);

		//if(nbAnchors > _scores.size()){
		//	_scores.resize(nbAnchors);
		//	_parents.resize(nbAnchors);
		//}

		//for(size_t i=0; i<nbAnchors; i++){
		//	_scores[i] = 0;
		//	_parents[i] = 0;
		//}

		//for (var i in global_list_of_point) {
		for (size_t i=0; i<anchors.size(); i++) {

			//if(_parent._print_debug) cout << "\t--- " << i << endl;
			
			//cout << i << " " << points.size() << endl;
			//size_t j = i;
			//var j = parseInt(i)

			float bestScore = 0;
			size_t bestPrevIndex = i;

			//cout << "a" << endl;
			argmaxPosition(anchors, scores, i, w, a, b, bestScore, bestPrevIndex, referenceRead, queryRead, maxChainingBand);

			
			if(bestPrevIndex != i) {
				scores[i] = bestScore;
				parents[i] = bestPrevIndex;
				//cout << "Add chain part: " << i << " " << bestPrevIndex << endl;
				//chain_part.union(i, best_prev_index);
			}
			else{
				scores[i] = w;
				parents[i] = -1;
			}
			/*
			if(bestScore < w){
				scores[i] = w;
				parents[i] = -1;
			}
			else{
				scores[i] = bestScore;
				parents[i] = bestPrevIndex;
			}
			*/

		}


		//for(size_t i=0; i<scores.size(); i++) {
		//	cout << i << "\t" << scores[i] << "\t" << parents[i] << endl;
		//}

		//if(_print_debug){
		//	cout << endl;
		//	for(size_t i=0; i<scores.size(); i++) {
		//		cout << "\tlala: " << i << " " << scores[i] << " " << getRoot(i, parents) << endl;
		//		F.push_back({i, scores[i], getRoot(i, parents)});
		//	}
		//}
		
		float maxScore = 0;
		size_t bestIndex = -1;

		for(size_t i=0; i<anchors.size(); i++) {
			if(scores[i] > maxScore){
				maxScore = scores[i];
				bestIndex = i;
			}
			//cout << "\tlala: " << i << " " << scores[i] << " " << getRoot(i, parents) << endl;
			//F.push_back({i, scores[i], getRoot(i, parents)});
		}

		if(_print_debug){
			cout << "\tParents:" << endl;
			for(size_t i=0; i<anchors.size(); i++){
				cout << "\t\t" << i << "\t" << anchors[i]._isReversed << "\t" << parents[i] << "\t" << scores[i] << endl;
			}
		}

		//cout << bestIndex << endl;

		vector<size_t> interval;

		size_t nullVal = -1;
		size_t index = bestIndex;
		while (index != nullVal) {
			interval.push_back(index);
			index = parents[index];
			//num_anchors += 1;
		}

		//cout << bestIndex << " " << index << endl;
		if(_print_debug){
			cout << "\tBest interval:" << endl;
			for(size_t i=0; i<interval.size(); i++){
				cout << "\t\t" << interval[i] << endl;
			}
		}

		std::reverse(interval.begin(), interval.end());
		chainIntervals.push_back({maxScore, interval});
		//std::sort(F.begin(), F.end(), [](const Point2& a, const Point2& b){
		//	if(a._score == b._score){
		//		return a._fromIndex > b._fromIndex;
		//	}
		//	return a._score > b._score;
		//});



		return chainIntervals; //F[0]._score;
	}

	//size_t _maxChainingBand;

	void argmaxPosition(const vector<Anchor>& anchors, vector<float>& scores, size_t i, float w, float a, float b, float& bestScore, size_t& bestPrevIndex, const MinimizerRead& referenceRead, const MinimizerRead& queryRead, const int64_t& maxChainingBand) {
		
		//cout << "---" << i << endl;
		bestScore = 0;
		bestPrevIndex = i;

		//cout << "\tglurp" << endl;
		//size_t bestJ = 0;
		//float v = 0;
		const Anchor& xi = anchors[i];
		//console.log("\t", xi);
		for (int64_t j = i-1; j >= 0; j--) {
		//for (int64_t k =0; k < i; k++) {

			//cout << "\t\t" << j << endl;

			//cout << "\t\t" << i << " " << k << " " << i-k << " " << maxBand << endl;
			if(i-j > maxChainingBand) break;
			const Anchor& xj = anchors[j];


			//if(xi._queryPosition - xj._queryPosition > 2500) continue;
			if(xi._isReversed != xj._isReversed) continue;
			if(xi._referencePosition == xj._referencePosition || xi._queryPosition == xj._queryPosition) continue;
			//cout << "\t" << i << " " << k << " " << xi._isReversed << " " << xj._isReversed << " " << xi._referencePosition << "-" << xi._queryPosition << " " << xj._referencePosition << "-" << xj._queryPosition << endl;




			int32_t d_q;// = abs(xi._queryPosition - xj._queryPosition);
			if(xi._isReversed){
				d_q = xj._queryPosition - xi._queryPosition;
			}
			else{
				d_q = xi._queryPosition - xj._queryPosition;
			}
			int32_t d_r = xi._referencePosition - xj._referencePosition;

			//cout << "\t" << xi._referencePosition << " " << xj._referencePosition << " " << d_r << endl;
			//if(d_r > 40) continue;
			//if(d_q > 40) continue;

			//if(xi._isReversed) {
			//	d_r = aprpf64 - acrpf64;
			//} else {
			//	d_r = ;
			//}
			//if(_parent._print_debug) cout << "\t\t" << j << "\t" << d_q << "\t" << d_r << "\t" << abs(d_r - d_q) << "\t" << referenceRead._minimizers[xi._referencePositionIndex] << "\t" << referenceRead._minimizers[xj._referencePositionIndex] << endl;
			//cout << "\t" << d_r << " " << d_q << endl;

			//if d_q > D_MAX_LIN_LENGTH || d_r > D_MAX_LIN_LENGTH {
			if(d_q > 5000 || d_r > 5000) {
				continue;
			}

			if(d_r <= 0) {
				continue;
			}

			int32_t gap = abs(d_r - d_q);
			//cout << "\t" << xi._referencePosition << "\t" << xj._referencePosition << "\t" << gap << endl;

			if(gap > 100) {
				continue;
			}
			
			if(xi._isReversed){
				if(xi._queryPosition > xj._queryPosition) continue;
			}
			else{
				if(xi._queryPosition < xj._queryPosition) continue;
			}
			
			//cout << w << " " << min(smallerDistance(xi,xj,xi._isReversed),w) << " " << gamma(xi,xj,w,a,b) << endl;
			//float anchor_score = min(smallerDistance(xi,xj,xi._isReversed),w) - gamma(xi,xj,w,a,b,xi._isReversed); //w - gap;
			float anchor_score = w - gap;

			float new_score = scores[j] + anchor_score;
			if (new_score > bestScore) {
				bestScore = new_score;
				bestPrevIndex = j;
				//if(_parent._print_debug) cout << "\t\tBest score:\t" << bestScore << "\t" << anchor_score << "\t" << j << "\t" << referenceRead._minimizers[xj._referencePositionIndex] << endl;
			}

			//map_params.anchor_score - gap


			//if (xi._referencePosition-xj._referencePosition < maxGapLength && xi._queryPosition-xj._queryPosition < maxGapLength && xi._queryPosition > xj._queryPosition) {
			//float newv = F[j] + 20 - gap;//+ 20 - gapLenth(xi,xj);//F[k] + min(smallerDistance(xi,xj),w) - gamma(xi,xj,w,a,b);
			//if (v < newv) {
			//	v = newv;
			//	bestJ = j;
			//}
			//}
		}
		//cout << "\tBest: " << bestPrevIndex << " " << bestScore << endl;
		//cout << "\tglurpppp" << endl;
		//jReturn = bestJ;
		//vReturn = v;
	}



	float smallerDistance(const Anchor& xi, const Anchor& xj, bool isQueryReversed) {
		if(isQueryReversed){
			return min((xj._queryPosition-xi._queryPosition), (xi._referencePosition-xj._referencePosition));
		}
		else{
			return min((xi._queryPosition-xj._queryPosition), (xi._referencePosition-xj._referencePosition));
		}
	}

	float gamma(const Anchor& xi, const Anchor& xj, float w, float a, float b, bool isQueryReversed) {
		float g = gapLength(xi,xj,isQueryReversed);
		return a*w*g + b*max(log2(g), (float)0);

	}

	//double computeDistance(const Point& p1, const Point& p2){

	//	double referenceDistance = p2.x - p1.x;
	//	double queryDistance = p2.y - p1.y;

	//	return abs(referenceDistance - queryDistance);
	//}

	float gapLength(const Anchor& xi, const Anchor& xj, bool isQueryReversed) {
		if(isQueryReversed){
			return abs((xj._queryPosition-xi._queryPosition)-(xi._referencePosition-xj._referencePosition));
		}
		else{
			return abs((xi._queryPosition-xj._queryPosition)-(xi._referencePosition-xj._referencePosition));
		}

	}

	size_t getRoot(size_t i, const vector<size_t>& parents) {
		size_t y = i;
		size_t z = parents[y];
		//cout << z << endl;
		while (z != -1) {
			y = z;
			z = parents[y];

			//cout << z << endl;
			//getchar();
		}

		return y;
	}



	void normalizeAlignment(vector<pair<u_int32_t, u_int32_t>>& alignment, const vector<MinimizerType>& referenceMinimizers, const vector<MinimizerType>& queryMinimizers){


		//while(true){

		//	bool isModification = false;

		for(size_t i=0; i<alignment.size(); i++){

			//bool isMod = false;

			const auto& it = alignment[i];
			u_int32_t referencePosition = it.first;
			u_int32_t queryPosition = it.second;

			//cout << endl;
			//cout << referencePosition << " " << referenceMinimizers.size() << endl;
			//cout << queryPosition << " " << queryMinimizers.size() << endl;
			//cout << endl;

			if(referencePosition == -1){ //insertion
				//cout << "I: " << "-1" << "\t" << queryPosition << "\t" << queryMinimizers[queryPosition] << endl;

				int64_t j = getNextMatchInReference(i, alignment);
				if(j == -1) continue;

				int64_t refPos = alignment[j].first;

				if(referenceMinimizers[refPos] == queryMinimizers[queryPosition]){
					//isModification = true;
					alignment[i].first = refPos;
					alignment[j].first = -1;
				}

				if(alignment[j].first == -1 && alignment[j].second == -1){
					alignment.erase(alignment.begin()+j);
				}

				//printAlignment(alignment, referenceMinimizers, queryMinimizers);
			}
			else if(queryPosition == -1){ //deletion

				int64_t j = getNextMatchInQuery(i, alignment);
				if(j == -1) continue;

				int64_t queryPos = alignment[j].second;

				if(referenceMinimizers[referencePosition] == queryMinimizers[queryPos]){
					//isModification = true;
					alignment[i].second = queryPos;
					alignment[j].second = -1;
				}

				if(alignment[j].first == -1 && alignment[j].second == -1){
					alignment.erase(alignment.begin()+j);
				}

				//printAlignment(alignment, referenceMinimizers, queryMinimizers);

				//cout << "D: " << referencePosition << "\t" << "-1" << "\t" << referenceMinimizers[referencePosition] << endl;
			}
			else if(referenceMinimizers[referencePosition] == queryMinimizers[queryPosition]){ //Match
				//cout << "M: " << referencePosition << "\t" << queryPosition << "\t" << referenceMinimizers[referencePosition] << "\t" << queryMinimizers[queryPosition] << endl;
			}
			else { //missmatch
				//cout << "X: " << referencePosition << "\t" << queryPosition << "\t" << referenceMinimizers[referencePosition] << "\t" << queryMinimizers[queryPosition] << endl;
			}


			//cout << "hum" << endl;
			//getchar();

			//if(isModification) break;
		}	

		//	if(!isModification) break;
		//}



	}

	int64_t getNextMatchInReference(size_t i, const vector<pair<u_int32_t, u_int32_t>>& alignment){
		for(size_t j=i; j<alignment.size(); j++){
			if(alignment[j].first != -1) return j;
		}

		return -1;
	}

	int64_t getNextMatchInQuery(size_t i, const vector<pair<u_int32_t, u_int32_t>>& alignment){
		for(size_t j=i; j<alignment.size(); j++){
			if(alignment[j].second != -1) return j;
		}

		return -1;
	}


};

#endif