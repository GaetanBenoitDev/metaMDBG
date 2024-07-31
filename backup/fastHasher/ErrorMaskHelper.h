#ifndef ErrorMaskHelper_h
#define ErrorMaskHelper_h

#include <vector>
#include <string>
#include <tuple>
#include "MBGCommon.h"

namespace MBG
{

std::pair<SequenceCharType, SequenceLengthType> multiRLECompress(const SequenceCharType& str, const SequenceLengthType& poses, const size_t maxMaskLength);
size_t maxCode();
SequenceCharType revCompRLE(const SequenceCharType& str);
CharType complement(const CharType original);
size_t codeMotifLength(const uint16_t code);

}

#endif
