

#ifndef MDBG_METAG_KMER
#define MDBG_METAG_KMER


#include <vector>
#include <string>
#include <cmath>

//#include "hasher.hpp"
//#include "enumerator.hpp"
//#include "fastmod.h"
//#include "../src/utils/fastHasher/FastHasher.h"



using namespace std;

typedef unsigned __int128 u_int128_t;
typedef u_int64_t MinimizerType;



//complement of one NT
const unsigned char comp_NT[4] = {  2,3,0,1  };
const char bin2NT[4] = {'A','C','T','G'};
const char binrev[4] = {2,3,0,1};

//reverse complement of 4NT,  ie one byte
const unsigned char revcomp_4NT[256] = {
 0xaa,
 0xea,
 0x2a,
 0x6a,
 0xba,
 0xfa,
 0x3a,
 0x7a,
 0x8a,
 0xca,
 0xa,
 0x4a,
 0x9a,
 0xda,
 0x1a,
 0x5a,
 0xae,
 0xee,
 0x2e,
 0x6e,
 0xbe,
 0xfe,
 0x3e,
 0x7e,
 0x8e,
 0xce,
 0xe,
 0x4e,
 0x9e,
 0xde,
 0x1e,
 0x5e,
 0xa2,
 0xe2,
 0x22,
 0x62,
 0xb2,
 0xf2,
 0x32,
 0x72,
 0x82,
 0xc2,
 0x2,
 0x42,
 0x92,
 0xd2,
 0x12,
 0x52,
 0xa6,
 0xe6,
 0x26,
 0x66,
 0xb6,
 0xf6,
 0x36,
 0x76,
 0x86,
 0xc6,
 0x6,
 0x46,
 0x96,
 0xd6,
 0x16,
 0x56,
 0xab,
 0xeb,
 0x2b,
 0x6b,
 0xbb,
 0xfb,
 0x3b,
 0x7b,
 0x8b,
 0xcb,
 0xb,
 0x4b,
 0x9b,
 0xdb,
 0x1b,
 0x5b,
 0xaf,
 0xef,
 0x2f,
 0x6f,
 0xbf,
 0xff,
 0x3f,
 0x7f,
 0x8f,
 0xcf,
 0xf,
 0x4f,
 0x9f,
 0xdf,
 0x1f,
 0x5f,
 0xa3,
 0xe3,
 0x23,
 0x63,
 0xb3,
 0xf3,
 0x33,
 0x73,
 0x83,
 0xc3,
 0x3,
 0x43,
 0x93,
 0xd3,
 0x13,
 0x53,
 0xa7,
 0xe7,
 0x27,
 0x67,
 0xb7,
 0xf7,
 0x37,
 0x77,
 0x87,
 0xc7,
 0x7,
 0x47,
 0x97,
 0xd7,
 0x17,
 0x57,
 0xa8,
 0xe8,
 0x28,
 0x68,
 0xb8,
 0xf8,
 0x38,
 0x78,
 0x88,
 0xc8,
 0x8,
 0x48,
 0x98,
 0xd8,
 0x18,
 0x58,
 0xac,
 0xec,
 0x2c,
 0x6c,
 0xbc,
 0xfc,
 0x3c,
 0x7c,
 0x8c,
 0xcc,
 0xc,
 0x4c,
 0x9c,
 0xdc,
 0x1c,
 0x5c,
 0xa0,
 0xe0,
 0x20,
 0x60,
 0xb0,
 0xf0,
 0x30,
 0x70,
 0x80,
 0xc0,
 0x0,
 0x40,
 0x90,
 0xd0,
 0x10,
 0x50,
 0xa4,
 0xe4,
 0x24,
 0x64,
 0xb4,
 0xf4,
 0x34,
 0x74,
 0x84,
 0xc4,
 0x4,
 0x44,
 0x94,
 0xd4,
 0x14,
 0x54,
 0xa9,
 0xe9,
 0x29,
 0x69,
 0xb9,
 0xf9,
 0x39,
 0x79,
 0x89,
 0xc9,
 0x9,
 0x49,
 0x99,
 0xd9,
 0x19,
 0x59,
 0xad,
 0xed,
 0x2d,
 0x6d,
 0xbd,
 0xfd,
 0x3d,
 0x7d,
 0x8d,
 0xcd,
 0xd,
 0x4d,
 0x9d,
 0xdd,
 0x1d,
 0x5d,
 0xa1,
 0xe1,
 0x21,
 0x61,
 0xb1,
 0xf1,
 0x31,
 0x71,
 0x81,
 0xc1,
 0x1,
 0x41,
 0x91,
 0xd1,
 0x11,
 0x51,
 0xa5,
 0xe5,
 0x25,
 0x65,
 0xb5,
 0xf5,
 0x35,
 0x75,
 0x85,
 0xc5,
 0x5,
 0x45,
 0x95,
 0xd5,
 0x15,
 0x55
};



class KmerDirect
{
public:
	/** Returns the value of the kmer.
	 * \return the kmer value as a Type object. */
	const u_int64_t& value  () const { return _value;   }

	/** This is a dummy function that always returns the value of the kmer in the forward direction, even when "which" is 1.
	 * it's provided for API compatibility with KmerCanonical
	 * We could compute revcomp(_value) when which=1, but KmerDirect doesn't know about its k-mer size.
	 * \param[in] which: dummy parameter
	 * \return the kmer value as a Type object. */
	const u_int64_t& value  (int which) const { if (which==1){ std::cout << "unsupported call to value(which) for KmerDirect" << std::endl; exit(1); } 
											return _value;   }

	/*_ Comparison operator between two instances.
		* \param[in] t : object to be compared to
		* \return true if the values are the same, false otherwise. */
	bool operator< (const KmerDirect& t) const  { return this->_value < t._value; };

	/** Set the value of the kmer
	 * \param[in] val : value to be set. */
	void set (const u_int64_t& val) { _value=val; }

	/** Tells whether the kmer is valid or not. It may be invalid if some unwanted
	 * nucleotides characters (like N) have been used to build it.
	 * \return true if valid, false otherwise. */
	bool isValid () const { return _isValid; }

	/* compatibility with KmerCanonical API */
	bool which () const { return true; }

	/* compatibility with KmerCanonical API */
	//Strand strand() const { return STRAND_FORWARD;  }

	/* compatibility with KmerCanonical API */
	const u_int64_t& forward() const { return value(); }

	/** Returns the reverse complement value of this canonical kmer.
	 * \return the reverse complement value */
	//const Type& revcomp() const { return  }


	u_int64_t _value;
	bool _isValid;
	//friend class ModelDirect;

	/** Extract a mmer from a kmer. This is done by using a mask on the kmer.
	 * \param[in] mask : mask to be applied to the current kmer
	 * \param[in] size : shift size (needed for some kmer classes but not all)
	 * \param[in] mmer_lut : lookup table of minimizers
	 * \return the extracted kmer.
	 */

	KmerDirect extract      (const u_int64_t& mask, size_t size, u_int64_t * mmer_lut)  {  KmerDirect output;  output.set (mmer_lut[(this->value() & mask)]);  return output;  }
	KmerDirect extractShift (const u_int64_t& mask, size_t size, u_int64_t * mmer_lut)  {  KmerDirect output = extract(mask,size,mmer_lut);  _value = _value >> 2;  return output;  }
};



class KmerCanonical
{
public:


	const u_int64_t& value  () const { return table[(int)choice];   }

	/** Returns the value of the kmer.
	 * \return the kmer value as a Type object. */
	const u_int64_t& value  (int& direction) const { direction = choice; return table[(int)choice];   }

	/** Returns the value of the kmer.
	 * \param[in] which: forward or reverse strand
	 * \return the kmer value as a Type object. */
	//const u_int64_t& value  (int which) const { return table[which];   }

	/** Comparison operator between two instances.
	 * \param[in] t : object to be compared to
	 * \return true if the values are the same, false otherwise. */
	bool operator< (const KmerDirect& t) const  { return this->value() < t.value(); };

	/** Set the value of the kmer. IMPORTANT: Not really a forward/revcomp couple,
	 * but may be useful for the minimizer default value.
	 * \param[in] val : value to be set (set to both forward and reverse complement). */
	void set (const u_int64_t& val)
	{
		table[0]=val;
		table[1]=val;
		choice = 0;
	}

	/** Set the forward/revcomp attributes. The canonical form is computed here.
	 * \param[in] forward : forward value
	 * \param[in] revcomp : reverse complement value.
	 */
	void set (const u_int64_t& forward, const u_int64_t& revcomp)
	{
		table[0]=forward;
		table[1]=revcomp;
		updateChoice ();
	}
	
	/** Tells whether the kmer is valid or not. It may be invalid if some unwanted
	 * nucleotides characters (like N) have been used to build it.
	 * \return true if valid, false otherwise. */
	bool isValid () const { return _isValid; }

	/** Returns the forward value of this canonical kmer.
	 * \return the forward value */
	const u_int64_t& forward() const { return table[0]; }

	/** Returns the reverse complement value of this canonical kmer.
	 * \return the reverse complement value */
	const u_int64_t& revcomp() const { return table[1]; }

	/** Tells which strand is used for the kmer.
	 * \return true if the kmer value is the forward value, false if it is the reverse complement value
	 */
	bool which () const { return choice==0 ? true : false; }

	/** Tells which strand is used.
	 * \return the used strand. */
	//Strand strand() const { return which() ? STRAND_FORWARD : STRAND_REVCOMP; }

	/* tells whether a kmer and its revcomp are identical */
	bool isPalindrome () const { return table[0] == table[1]; }

	u_int64_t table[2];  char choice;
	
	bool _isValid;
	void updateChoice () { choice = (table[0] < table[1]) ? 0 : 1; }
	//friend class ModelCanonical;

	/** Extract a mmer from a kmer. This is done by using a mask on the kmer.
	 * \param[in] mask : mask to be applied to the current kmer
	 * \param[in] size : shift size (needed for some kmer classes but not all)
	 * \param[in] mmer_lut : lookup table of minimizers
	 * \return the extracted kmer.
	 */
	KmerCanonical extract (const u_int64_t& mask, size_t size, u_int64_t* mmer_lut)
	{

		KmerCanonical output;
		
		output.set(mmer_lut[(this->table[0] & mask)]); //no need to recomp updateChoice with this
		//mmer_lut takes care of revcomp and forbidden mmers
		//output.set (this->table[0] & mask, (this->table[1] >> size) & mask);
		//output.updateChoice();
		return output;
	}


	KmerCanonical extractShift (const u_int64_t& mask, size_t size, u_int64_t * mmer_lut)
	{
		KmerCanonical output = extract (mask, size,mmer_lut);
		table[0] = table[0] >> 2;   table[1] = table[1] << 2;  updateChoice();
		return output;
	}
	
};

class KmerModel{
public:

    typedef std::pair<char,char> ConvertChar;
    struct ConvertASCII    { static ConvertChar get (const char* buffer, size_t idx)  { return ConvertChar((buffer[idx]>>1) & 3, (buffer[idx]>>3) & 1); }};

	
	size_t  _kmerSize;
	u_int64_t  _kmerMask;
	u_int64_t _revcompTable[4];

	typedef KmerCanonical Kmer;
	//typedef Kmer<span>::KmerCanonical Kmer;

	KmerModel(size_t kmerSize){

		_kmerSize = kmerSize;

		u_int64_t un;
		un = 1;
		_kmerMask = (un << (_kmerSize*2)) - un;

		size_t shift = 2*(_kmerSize-1);

		/** The _revcompTable is a shortcut used while computing revcomp recursively. */
		/** Important: don't forget the Type cast, otherwise the result in only on 32 bits. */
		for (size_t i=0; i<4; i++)   {  u_int64_t tmp; tmp = comp_NT[i];  _revcompTable[i] = tmp << shift;  }
	}

	//template<class Convert>
	int polynom (const char* seq, u_int64_t& kmer, size_t startIndex)  const
	{
		ConvertChar c;
		int badIndex = -1;

		kmer = 0;
		for (size_t i=0; i<_kmerSize; ++i)
		{
			//cout << seq[i] << endl;
			c = ConvertASCII::get(seq,i+startIndex);

			kmer = (kmer<<2) + c.first;

			if (c.second)  { badIndex = i; }
		}

		return badIndex;
	}

	inline static u_int64_t revcomp(const u_int64_t& x, size_t sizeKmer)
    {
        u_int64_t res = x;

        unsigned char* kmerrev  = (unsigned char *) (&(res));
        unsigned char* kmer     = (unsigned char *) (&(x));

        for (size_t i=0; i<8; ++i)  {  kmerrev[8-1-i] = revcomp_4NT [kmer[i]];  }

        return (res >> (2*( 32 - sizeKmer))) ;
    }

	u_int64_t reverse (const u_int64_t& kmer)  const  { return revcomp (kmer, this->_kmerSize); }
	/*
	std::string toString (u_int64_t val, size_t sizeKmer) const
    {
        char seq[sizeKmer+1];
        char bin2NT[4] = {'A','C','T','G'};

        for (size_t i=0; i<sizeKmer; i++)  {  seq[sizeKmer-i-1] = bin2NT [(*this)[i]];  }
        seq[sizeKmer]='\0';
        return seq;
    }*/

	bool iterate (const char* seq, size_t length, vector<u_int64_t>& tmpKmers, vector<u_int8_t>& tmpKmersDirection) const{

		int32_t nbKmers = length - _kmerSize + 1;
		if (nbKmers <= 0)  { return false; }

		tmpKmers.resize(nbKmers);
		tmpKmersDirection.resize(nbKmers);

		Kmer result;
		int indexBadChar = first(seq, result, 0);
		//int indexBadChar = static_cast<const ModelCanonical*>(this)->template first<Convert> (seq, result, 0);

		size_t idxComputed = 0;

		/*
		cout << result.value() << endl;

		u_int64_t kmerVal = result.value();
		Type un;
		un.setVal(kmerVal);

		cout << un.toString(_kmerSize) << endl;

		Type deux;
		deux.setVal(revcomp(kmerVal, _kmerSize));

		cout << deux.toString(_kmerSize) << endl;

		exit(1);
		//this->notification<Callback> (result, idxComputed, callback);
		*/

		int direction;
		tmpKmers[idxComputed] = result.value(direction);
		tmpKmersDirection[idxComputed] = direction;

		if(!result.isValid()) tmpKmers[idxComputed] = -1; //Kmer with max value will be skipped as minimizers
		idxComputed += 1;
		
		for (size_t idx=_kmerSize; idx<length; idx++)
		{
			ConvertChar c = ConvertASCII::get (seq, idx);

			if (c.second)  { indexBadChar = _kmerSize-1; }
			else           { indexBadChar--;     }

			next(c.first, result, indexBadChar<0);
			tmpKmers[idxComputed] = result.value(direction);
			tmpKmersDirection[idxComputed] = direction;
			if(!result.isValid()) tmpKmers[idxComputed] = -1; //Kmer with max value will be skipped as minimizers

			//if(!result.isValid()) cout << "img" << endl;
			//tmpKmers.push_back(result.value());
			//cout << result.value() << endl;
			//static_cast<const ModelCanonical*>(this)->template next<Convert> (c.first, result, indexBadChar<0);

			//this->notification<Callback> (result, ++idxComputed, callback);
			idxComputed += 1;
		}

		return true;
	}

	int first (const char* seq, Kmer& value, size_t startIndex)   const
	{

		int result = polynom(seq, value.table[0], startIndex);
		value._isValid = result < 0;
		value.table[1] = this->reverse (value.table[0]);
		value.updateChoice();
		return result;
	}

	void  next (char c, Kmer& value, bool isValid)   const
	{
		value.table[0] = ( (value.table[0] << 2) +  c                          ) & this->_kmerMask;
		value.table[1] = ( (value.table[1] >> 2) +  this->_revcompTable[(int)c]) & this->_kmerMask;
		value._isValid = isValid;

		value.updateChoice();
	}

	/*
	uint64_t getHash(const u_int64_t &k) const
	{
		return hash1(k, 0);
	}

	void getHash2(const u_int64_t &k) const
	{
		hash2(k, 1LL);
	}*/

};
















class KmerDirect128
{
public:
	/** Returns the value of the kmer.
	 * \return the kmer value as a Type object. */
	const u_int128_t& value  () const { return _value;   }

	/** This is a dummy function that always returns the value of the kmer in the forward direction, even when "which" is 1.
	 * it's provided for API compatibility with KmerCanonical
	 * We could compute revcomp(_value) when which=1, but KmerDirect doesn't know about its k-mer size.
	 * \param[in] which: dummy parameter
	 * \return the kmer value as a Type object. */
	const u_int128_t& value  (int which) const { if (which==1){ std::cout << "unsupported call to value(which) for KmerDirect" << std::endl; exit(1); } 
											return _value;   }

	/*_ Comparison operator between two instances.
		* \param[in] t : object to be compared to
		* \return true if the values are the same, false otherwise. */
	bool operator< (const KmerDirect128& t) const  { return this->_value < t._value; };

	/** Set the value of the kmer
	 * \param[in] val : value to be set. */
	void set (const u_int128_t& val) { _value=val; }

	/** Tells whether the kmer is valid or not. It may be invalid if some unwanted
	 * nucleotides characters (like N) have been used to build it.
	 * \return true if valid, false otherwise. */
	bool isValid () const { return _isValid; }

	/* compatibility with KmerCanonical API */
	bool which () const { return true; }

	/* compatibility with KmerCanonical API */
	//Strand strand() const { return STRAND_FORWARD;  }

	/* compatibility with KmerCanonical API */
	const u_int128_t& forward() const { return value(); }

	/** Returns the reverse complement value of this canonical kmer.
	 * \return the reverse complement value */
	//const Type& revcomp() const { return  }


	u_int128_t _value;
	bool _isValid;
	//friend class ModelDirect;

	/** Extract a mmer from a kmer. This is done by using a mask on the kmer.
	 * \param[in] mask : mask to be applied to the current kmer
	 * \param[in] size : shift size (needed for some kmer classes but not all)
	 * \param[in] mmer_lut : lookup table of minimizers
	 * \return the extracted kmer.
	 */

	KmerDirect128 extract      (const u_int128_t& mask, size_t size, u_int128_t * mmer_lut)  {  KmerDirect128 output;  output.set (mmer_lut[(this->value() & mask)]);  return output;  }
	KmerDirect128 extractShift (const u_int128_t& mask, size_t size, u_int128_t * mmer_lut)  {  KmerDirect128 output = extract(mask,size,mmer_lut);  _value = _value >> 2;  return output;  }
};



class KmerCanonical128
{
public:

	/** Returns the value of the kmer.
	 * \return the kmer value as a Type object. */
	const u_int128_t& value  () const { return table[(int)choice];   }

	/** Returns the value of the kmer.
	 * \param[in] which: forward or reverse strand
	 * \return the kmer value as a Type object. */
	const u_int128_t& value  (int which) const { return table[which];   }

	/** Comparison operator between two instances.
	 * \param[in] t : object to be compared to
	 * \return true if the values are the same, false otherwise. */
	bool operator< (const KmerDirect128& t) const  { return this->value() < t.value(); };

	/** Set the value of the kmer. IMPORTANT: Not really a forward/revcomp couple,
	 * but may be useful for the minimizer default value.
	 * \param[in] val : value to be set (set to both forward and reverse complement). */
	void set (const u_int128_t& val)
	{
		table[0]=val;
		table[1]=val;
		choice = 0;
	}

	/** Set the forward/revcomp attributes. The canonical form is computed here.
	 * \param[in] forward : forward value
	 * \param[in] revcomp : reverse complement value.
	 */
	void set (const u_int128_t& forward, const u_int128_t& revcomp)
	{
		table[0]=forward;
		table[1]=revcomp;
		updateChoice ();
	}
	
	/** Tells whether the kmer is valid or not. It may be invalid if some unwanted
	 * nucleotides characters (like N) have been used to build it.
	 * \return true if valid, false otherwise. */
	bool isValid () const { return _isValid; }

	/** Returns the forward value of this canonical kmer.
	 * \return the forward value */
	const u_int128_t& forward() const { return table[0]; }

	/** Returns the reverse complement value of this canonical kmer.
	 * \return the reverse complement value */
	const u_int128_t& revcomp() const { return table[1]; }

	/** Tells which strand is used for the kmer.
	 * \return true if the kmer value is the forward value, false if it is the reverse complement value
	 */
	bool which () const { return choice==0 ? true : false; }

	/** Tells which strand is used.
	 * \return the used strand. */
	//Strand strand() const { return which() ? STRAND_FORWARD : STRAND_REVCOMP; }

	/* tells whether a kmer and its revcomp are identical */
	bool isPalindrome () const { return table[0] == table[1]; }

	u_int128_t table[2];  char choice;
	
	bool _isValid;
	void updateChoice () { choice = (table[0] < table[1]) ? 0 : 1; }
	//friend class ModelCanonical;

	/** Extract a mmer from a kmer. This is done by using a mask on the kmer.
	 * \param[in] mask : mask to be applied to the current kmer
	 * \param[in] size : shift size (needed for some kmer classes but not all)
	 * \param[in] mmer_lut : lookup table of minimizers
	 * \return the extracted kmer.
	 */
	KmerCanonical128 extract (const u_int128_t& mask, size_t size, u_int128_t* mmer_lut)
	{

		KmerCanonical128 output;
		
		output.set(mmer_lut[(this->table[0] & mask)]); //no need to recomp updateChoice with this
		//mmer_lut takes care of revcomp and forbidden mmers
		//output.set (this->table[0] & mask, (this->table[1] >> size) & mask);
		//output.updateChoice();
		return output;
	}


	KmerCanonical128 extractShift (const u_int128_t& mask, size_t size, u_int128_t * mmer_lut)
	{
		KmerCanonical128 output = extract (mask, size,mmer_lut);
		table[0] = table[0] >> 2;   table[1] = table[1] << 2;  updateChoice();
		return output;
	}
	
	/*
	u_int8_t  operator[]  (size_t idx) const    {  
		u_int128_t v = value();
        return (v[idx/32] >> (2*(idx % 32))) & 3; }

	std::string toString (size_t kmerSize) const
    {
        char seq[kmerSize+1];
        char bin2NT[4] = {'A','C','T','G'};

        for (size_t i=0; i<kmerSize; i++)  {  seq[kmerSize-i-1] = bin2NT [(*this)[i]];  }
        seq[kmerSize]='\0';
        return seq;
    }
	*/

};

class KmerModel128{
public:

    typedef std::pair<char,char> ConvertChar;
    struct ConvertASCII    { static ConvertChar get (const char* buffer, size_t idx)  { return ConvertChar((buffer[idx]>>1) & 3, (buffer[idx]>>3) & 1); }};

	
	size_t  _kmerSize;
	u_int128_t  _kmerMask;
	u_int128_t _revcompTable[4];

	typedef KmerCanonical128 Kmer;
	//typedef Kmer<span>::KmerCanonical128 Kmer;

	KmerModel128(size_t kmerSize){

		_kmerSize = kmerSize;

		u_int128_t un;
		un = 1;
		_kmerMask = (un << (_kmerSize*2)) - un;

		size_t shift = 2*(_kmerSize-1);

		/** The _revcompTable is a shortcut used while computing revcomp recursively. */
		/** Important: don't forget the Type cast, otherwise the result in only on 32 bits. */
		for (size_t i=0; i<4; i++)   {  u_int128_t tmp; tmp = comp_NT[i];  _revcompTable[i] = tmp << shift;  }
	}

	//template<class Convert>
	int polynom (const char* seq, u_int128_t& kmer, size_t startIndex)  const
	{
		ConvertChar c;
		int badIndex = -1;

		kmer = 0;
		for (size_t i=0; i<_kmerSize; ++i)
		{
			//cout << seq[i] << endl;
			c = ConvertASCII::get(seq,i+startIndex);

			kmer = (kmer<<2) + c.first;

			if (c.second)  { badIndex = i; }
		}

		return badIndex;
	}

	inline static u_int128_t revcomp(const u_int128_t& x, size_t sizeKmer)
    {
        u_int128_t res = x;

        unsigned char* kmerrev  = (unsigned char *) (&(res));
        unsigned char* kmer     = (unsigned char *) (&(x));

        for (size_t i=0; i<8; ++i)  {  kmerrev[8-1-i] = revcomp_4NT [kmer[i]];  }

        return (res >> (2*( 32 - sizeKmer))) ;
    }

	u_int128_t reverse (const u_int128_t& kmer)  const  { return revcomp (kmer, this->_kmerSize); }
	


	bool iterate (const char* seq, size_t length, vector<u_int128_t>& tmpKmers) const{

		int32_t nbKmers = length - _kmerSize + 1;
		if (nbKmers <= 0)  { return false; }

		tmpKmers.resize(nbKmers);

		Kmer result;
		int indexBadChar = first(seq, result, 0);
		//int indexBadChar = static_cast<const ModelCanonical*>(this)->template first<Convert> (seq, result, 0);

		size_t idxComputed = 0;

		/*
		cout << result.value() << endl;

		u_int128_t kmerVal = result.value();
		Type un;
		un.setVal(kmerVal);

		cout << un.toString(_kmerSize) << endl;

		Type deux;
		deux.setVal(revcomp(kmerVal, _kmerSize));

		cout << deux.toString(_kmerSize) << endl;

		exit(1);
		//this->notification<Callback> (result, idxComputed, callback);
		*/
		tmpKmers[idxComputed] = result.value();
		if(!result.isValid()) tmpKmers[idxComputed] = -1; //Kmer with max value will be skipped as minimizers
		idxComputed += 1;
		
		for (size_t idx=_kmerSize; idx<length; idx++)
		{
			ConvertChar c = ConvertASCII::get (seq, idx);

			if (c.second)  { indexBadChar = _kmerSize-1; }
			else           { indexBadChar--;     }

			next(c.first, result, indexBadChar<0);
			tmpKmers[idxComputed] = result.value();
			if(!result.isValid()) tmpKmers[idxComputed] = -1; //Kmer with max value will be skipped as minimizers

			//if(!result.isValid()) cout << "img" << endl;
			//tmpKmers.push_back(result.value());
			//cout << result.value() << endl;
			//static_cast<const ModelCanonical*>(this)->template next<Convert> (c.first, result, indexBadChar<0);

			//this->notification<Callback> (result, ++idxComputed, callback);
			idxComputed += 1;
		}

		return true;
	}

	int first (const char* seq, Kmer& value, size_t startIndex)   const
	{

		int result = polynom(seq, value.table[0], startIndex);
		value._isValid = result < 0;
		value.table[1] = this->reverse (value.table[0]);
		value.updateChoice();
		return result;
	}

	void  next (char c, Kmer& value, bool isValid)   const
	{
		value.table[0] = ( (value.table[0] << 2) +  c                          ) & this->_kmerMask;
		value.table[1] = ( (value.table[1] >> 2) +  this->_revcompTable[(int)c]) & this->_kmerMask;
		value._isValid = isValid;

		value.updateChoice();
	}


	/*
	uint64_t getHash(const u_int64_t &k) const
	{
		return hash1(k, 0);
	}

	void getHash2(const u_int64_t &k) const
	{
		hash2(k, 1LL);
	}*/

};









/*
template <typename Hasher>
struct mod_sampling {
    static std::string name() { return "mod_sampling"; }

    mod_sampling(uint64_t w, uint64_t k, uint64_t t, uint64_t seed)
        : m_w(w), m_k(k), m_t(t), m_seed(seed), m_enum_tmers(w + k - t, t, seed) {
        m_M_w = fastmod::computeM_u32(m_w);
    }

    /// Sample from a single window.
    uint64_t sample(char const* window) {

		//cout << "\t---" << endl;
		//cout << string(window) << endl;
        const uint64_t num_tmers = (m_w + m_k - 1) - m_t + 1;
        uint64_t p = -1;
        typename Hasher::hash_type min_hash(-1);
        // Find the leftmost tmer with minimal hash.
        for (uint64_t i = 0; i != num_tmers; ++i) {

            //char const* tmer = window + i;

			string kmerForward(window + i, m_t);
			string kmerReverse(window + i, m_t);
			toReverseComplement(kmerReverse);

			auto hash = Hasher::hash(kmerForward.c_str(), m_w, m_t, m_seed);
			auto hashRev = Hasher::hash(kmerReverse.c_str(), m_w, m_t, m_seed);

			if(hashRev < hash){
				hash = hashRev;
			}
			//string s1 = string(window + i, m_t);
			//string s2 = string(window + i, m_t);
			//toReverseComplement(s2);

			//cout << "\t" << i << " " << kmerForward << " " << kmerReverse << " " << hash << endl;
			//auto hashRev = Hasher::hash(kmerReverse.c_str(), m_w, m_t, m_seed);

			//u_int8_t direction = 0;
			//if(hashRev < hash){
			//	hash = hashRev;
			//	direction = 1;
			//}

			//cout << kmerForward << " " << hash << endl;
			//cout << kmerReverse << " " << hashRev << endl;

			//getchar();

			//if(kmerForward == "TTGCACTCATGAA"){
			//	cout << hash << " " << (int)direction << endl;
			//}
			//if(kmerForward == "TTCATGAGTGCAA"){
			//	cout << hash << " " << (int)direction << endl;
			//}
			//cout << i << " " << m_w << " " << m_t << endl;
			//getchar();
			//string kmer_forward(window + i, );

			//cout << window + i
            
            if (hash < min_hash) {
                min_hash = hash;
                p = i;
            }
        }
        assert(p < num_tmers);
        uint64_t pos = fastmod::fastmod_u32(p, m_M_w, m_w);  // p % m_w

		//cout << "\t" << pos << endl;
        return pos;
    }

    /// Sample from a stream.
    /// If `clear`, this is the first call.
    uint64_t sample(char const* window, bool clear) {
        m_enum_tmers.eat(window, clear);
        uint64_t p = m_enum_tmers.next();
        return fastmod::fastmod_u32(p, m_M_w, m_w);  // p % m_w
    }

	void toReverseComplement(string& seq) {

		//string seq_rc;

		//seq_rc.clear();
		//seq_rc.reserve(seq.size());

		size_t size = seq.size();

		std::reverse(seq.begin(), seq.end());
		for (std::size_t i = 0; i < size; ++i){
			
			switch (seq[i]){
			case 'A':
				seq[i] = 'T';
				break;    
			case 'C':
				seq[i] = 'G';
				break;
			case 'G':
				seq[i] = 'C';
				break;
			case 'T':
				seq[i] = 'A';
				break;
			//default:
			//	seq[i] = seq[i];
			//	break;
			}
		}


	}

private:
    uint64_t m_w, m_k, m_t, m_seed;
    uint64_t m_M_w;
    minimizers::enumerator<Hasher> m_enum_tmers;
};
*/






class MinimizerParser{

public:

	u_int16_t _minimizerSize;
	KmerModel* _kmerModel;
	u_int32_t _seed;
	u_int64_t* _hash_otpt;
	double _minimizerBound;
	//double _minimizerBound_int32;
	size_t _trimBps;

	MinimizerParser(u_int16_t minimizerSize, double minimizerDensity){
		_minimizerSize = minimizerSize;
		_kmerModel = new KmerModel(_minimizerSize);
		_seed = 42;
		_hash_otpt = new u_int64_t[2];
		MinimizerType maxHashValue = -1;
		_minimizerBound = minimizerDensity * maxHashValue;

		//u_int32_t maxHashValue_int32 = -1;
		//_minimizerBound_int32 = minimizerDensity * maxHashValue_int32;
		_trimBps = 1; 
	}

	~MinimizerParser(){
		delete[] _hash_otpt;
	}

	void parse(const string& seq, vector<MinimizerType>& minimizers, vector<u_int32_t>& minimizersPos, vector<u_int8_t>& minimizersDirection){

		minimizers.clear();
		minimizersPos.clear();
		minimizersDirection.clear();

		//parseMod(seq, 380, 13, minimizers, minimizersPos, minimizersDirection);
		//return;
		//parse_windowed(seq, minimizers, minimizersPos, 200);
		//return;
		//0 0
		//100000 3278692
		//200000 6550207
		//300000 9827077


		vector<u_int64_t> kmers;
		vector<u_int8_t> kmerDirections;
		_kmerModel->iterate(seq.c_str(), seq.size(), kmers, kmerDirections);

		if(kmers.size() == 0) return;

		for(u_int64_t pos=_trimBps; pos<kmers.size()-_trimBps; pos++){

			//cout << itKmer->value().getVal() << endl;
			//cout << itKmer->value() << endl;
			//cout << pos << " " << (rleSequence.size()-_minimizerSize) << endl;

			/*
			if(pos == 0){
				//cout << "lala1" << endl;
				//pos += 1;
				continue;
			}
			else if(pos == kmers.size()-1; //seq.size()-_minimizerSize){
				//cout << "lala2" << endl;
				continue;
			}*/

			//if(!itKmer->value().isValid()) continue;
			//kmer_type kmerMin = min(itKmer->value(), revcomp(itKmer->value(), _kmerSize));
			//if(lala < 100 ) cout << model.toString(itKmer->value()) << endl;
			//lala += 1;
			u_int64_t kmerValue = kmers[pos];
			//MinimizerType minimizer = revhash(kmerValue);
			MurmurHash3_x64_128 ((const char*)&kmerValue, sizeof(kmerValue), _seed, _hash_otpt);
			MinimizerType minimizer = _hash_otpt[0];

			//u_int32_t minimizer_32bit = minimizer;;

			//if(minimizerCounts[minimizer] > 1000) cout << minimizer << endl;
			//double kmerHashed_norm = ((double) minimizer) / maxHashValue;
			if(minimizer < _minimizerBound){//_minimizerDensity){
			//if(minimizer_32bit < _minimizerBound_int32){//_minimizerDensity){


				minimizers.push_back(minimizer);
				minimizersPos.push_back(pos);
				minimizersDirection.push_back(kmerDirections[pos]);

				
				//minimizerCounts[minimizer] += 1;
				//nbMinimizers += 1;
				//"reprise: systeme de functor pour les reads et kmers?"
				

			}
			
			//pos += 1;
		}

	}

	uint64_t revhash(uint64_t x) {

		x = ((x >> 32) ^ x) * 0xD6E8FEB86659FD93;
		x = ((x >> 32) ^ x) * 0xD6E8FEB86659FD93;
		x = ((x >> 32) ^ x);
		return x;
	}



	uint64_t unrevhash(uint64_t x) {

		x = ((x >> 32) ^ x) * 0xCFEE444D8B59A89B;
		x = ((x >> 32) ^ x) * 0xCFEE444D8B59A89B;
		x = ((x >> 32) ^ x);
		return x;
	}

	/*
	void parseMod(const string& seq, uint64_t w, uint64_t k, vector<MinimizerType>& minimizers, vector<u_int32_t>& minimizersPos, vector<u_int8_t>& minimizersDirection){                                                       

   		//mod_sampling(uint64_t w, uint64_t k, uint64_t t, uint64_t seed)
        //: m_w(w), m_k(k), m_t(t), m_seed(seed), m_enum_tmers(w + k - t, t, seed) {
        //m_M_w = 
		//unordered_map<u_int64_t, u_int64_t> minimizer_to_abundance;

		vector<u_int64_t> hashes;

		vector<u_int64_t> kmers;
		vector<u_int8_t> kmerDirections;
		_kmerModel->iterate(seq.c_str(), seq.size(), kmers, kmerDirections);

		if(kmers.size() == 0) return;
		if(kmers.size() < w) return;

		for(u_int64_t pos=0; pos<kmers.size(); pos++){

			//cout << itKmer->value().getVal() << endl;
			//cout << itKmer->value() << endl;
			//cout << pos << " " << (rleSequence.size()-_minimizerSize) << endl;


			//if(!itKmer->value().isValid()) continue;
			//kmer_type kmerMin = min(itKmer->value(), revcomp(itKmer->value(), _kmerSize));
			//if(lala < 100 ) cout << model.toString(itKmer->value()) << endl;
			//lala += 1;
			u_int64_t kmerValue = kmers[pos];
			MurmurHash3_x64_128 ((const char*)&kmerValue, sizeof(kmerValue), _seed, _hash_otpt);
			u_int64_t minimizer = _hash_otpt[0];
			hashes.push_back(minimizer);

			//cout << "\t" << pos << " " << hashes[pos] << endl;
			//minimizer_to_abundance[hashes[pos]] += 1;
		}



		//vector<u_int64_t> hashesFiltered;
		//vector<u_int8_t> kmerDirectionsFilered;

		//for(size_t i=0; i<hashes.size(); i++){
		//	if(minimizer_to_abundance[hashes[i]] > 10) continue;

		//	hashesFiltered.push_back(hashes[i]);

		//}

		

		size_t m_M_w = fastmod::computeM_u32(w);

    	const uint64_t r = 4;
        const uint64_t t = r + ((k - r) % w);
		//mod_sampling<minimizers::hasher64_type> sampler(w, k, t, 42);
		//const uint64_t l = w + k - 1;  // num. symbols in window

		unordered_set<u_int32_t> sampledPositions;

    	const uint64_t l = w + k - 1;  // num. symbols in window
		const uint64_t num_windows = seq.size() - l;
		//cout << "Nb hashes: " << seq.size() << " " << hashes.size() << endl;
		//cout << "Nb window: " << num_windows << endl;
		//const uint64_t num_windows = hashes.size() - w + 1;
		//size_t tot_num_windows += num_windows;


		for (size_t i = 0; i < num_windows; i++) {

			size_t p = modSampleWindow(hashes, i, w, k, t, m_M_w);
			//cout << "---- " << i << " " << p << endl;
			sampledPositions.insert(p+i);
		}
		

		for(u_int32_t pos : sampledPositions){
			//if(pos == 0) continue;
			//if(pos + k >= seq.size()) continue;
			minimizersPos.push_back(pos);
		}

		std::sort(minimizersPos.begin(), minimizersPos.end());

		for(u_int64_t pos : minimizersPos){
			//if(minimizer_to_abundance[hashes[pos]] > 10) continue;

			minimizers.push_back(hashes[pos]);
			minimizersDirection.push_back(kmerDirections[pos]);

			//cout << "\t" << pos << " " << hashes[pos] << endl;
			//cout << m << " " << hashes[m] << endl;
		}

		//cout << kmers.size() << " " << minimizers.size() << " " << (long double) minimizers.size() /kmers.size() << endl;

		//modSample(seq, hashes, w, k);

	}
	*/
	/*
	void modSample(const string& seq, const vector<u_int64_t>& hashes, size_t w, size_t k, vector<MinimizerType>& minimizers, vector<u_int32_t>& minimizersPos, vector<u_int8_t>& minimizersDirection){

		//for(u_int64_t m : hashes){
			//cout << "\t" << m << endl;
		//}
		//size_t m_w = w;
		size_t m_M_w = fastmod::computeM_u32(w);

    	const uint64_t r = 4;
        const uint64_t t = r + ((k - r) % w);
		//mod_sampling<minimizers::hasher64_type> sampler(w, k, t, 42);
		//const uint64_t l = w + k - 1;  // num. symbols in window

		unordered_set<u_int32_t> sampledPositions;

    	const uint64_t l = w + k - 1;  // num. symbols in window
		const uint64_t num_windows = seq.size() - l + 1;
		//cout << "Nb hashes: " << seq.size() << " " << hashes.size() << endl;
		//cout << "Nb window: " << num_windows << endl;
		//const uint64_t num_windows = hashes.size() - w + 1;
		//size_t tot_num_windows += num_windows;
		for (size_t i = 0; i < num_windows; i++) {

			size_t p = modSampleWindow(hashes, i, w, k, t, m_M_w);
			//cout << "---- " << i << " " << p << endl;
			sampledPositions.insert(p+i);
		}
		

		for(u_int32_t pos : sampledPositions){
			//if(pos == 0) continue;
			//if(pos + k >= seq.size()) continue;
			minimizersPos.push_back(pos);
		}

		std::sort(minimizersPos.begin(), minimizersPos.end());

		for(u_int64_t pos : minimizersPos){
			minimizers.push_back(hashes[pos]);

			//cout << m << " " << hashes[m] << endl;
		}


	}
	*/
	/*
	u_int64_t modSampleWindow(const vector<u_int64_t>& hashes, size_t startI, size_t m_w, size_t m_k, size_t m_t, size_t m_M_w){

        const uint64_t num_tmers = (m_w + m_k - 1) - m_t + 1;
		//cout << "Num tmers: " << num_tmers << endl;

        u_int64_t p = -1;
		u_int64_t minHash = -1;

		for(size_t i=startI; i<startI+num_tmers; i++){

			//cout << i << " " << hashes.size() << " " << hashes[i] << endl;
            if (hashes[i] < minHash) {
                minHash = hashes[i];
                p = i-startI;
            }

		}

		//cout << "Min: " << minHash << " " << p << endl;
		//getchar();
        return fastmod::fastmod_u32(p, m_M_w, m_w);  // p % m_w
    }

	*/
	/*
	void parse_windowed(const string& seq, vector<u_int64_t>& minimizers, vector<u_int64_t>& minimizersPos, size_t w){

		minimizers.clear();
		minimizersPos.clear();

		//0 0
		//100000 3278692
		//200000 6550207
		//300000 9827077

		vector<u_int64_t> kmers;
		_kmerModel->iterate(seq.c_str(), seq.size(), kmers);

		if(kmers.size() == 0) return;

		//unordered_set<u_int64_t> selectedMinimizers;

		int lastSelectedW = 0;
		vector<u_int64_t> kmersCurrent;
		for(u_int64_t pos=1; pos<kmers.size()-1; pos++){

			u_int64_t kmerValue = kmers[pos];
			MurmurHash3_x64_128 ((const char*)&kmerValue, sizeof(kmerValue), _seed, _hash_otpt);
			u_int64_t kemrHashed = _hash_otpt[0];

			kmersCurrent.push_back(kemrHashed);

			if(pos >= w){
				u_int64_t minimizer = -1; 
				for(u_int64_t kmer : kmersCurrent){
					if(kmer < minimizer){
						minimizer = kmer;
					}
				}
				//cout << kmers[0] << endl;
				//auto it = minimizers_0_set.find(kmers[0]);
				//minimizers_0_set.erase(it);

				kmersCurrent.erase(kmersCurrent.begin());

				lastSelectedW += 1;
				//if
				//selectedMinimizers.insert(minimizer);
				//minimizersPos.push_back(pos);
				if(minimizers.size() == 0 || minimizer != minimizers[minimizers.size()-1] || lastSelectedW >= w){
					minimizers.push_back(minimizer);
					minimizersPos.push_back(pos);
					lastSelectedW = 0;
				}
				//	minimizers_0.push_back(minimizer);
				//}
				//cout << "Min: " << (*minimizers_0_set.begin()) << endl;
			}

		}


	}
	*/

	/*

	void parseSyncmers(const string& readSeq, vector<u_int64_t>& minimizers, vector<u_int64_t>& minimizersPos){ //size_t kmerSize, size_t windowSize, size_t minimizerSize


		size_t kmerSize = 13;
		size_t windowSize = 1000;
		size_t smerSize = 5;

		minimizers.clear();
		minimizersPos.clear();

		vector<u_int16_t> seq;

		for(char nt : readSeq){
			switch(nt)
			{
				case 'a':
				case 'A':
					seq.push_back(0);
					break;
				case 'c':
				case 'C':
					seq.push_back(1);
					break;
				case 'g':
				case 'G':
					seq.push_back(2);
					break;
				case 't':
				case 'T':
					seq.push_back(3);
					break;
			}
		}



		if (seq.size() < kmerSize) return;
		// keep the same smerOrder to reduce mallocs which destroy multithreading performance
		thread_local std::vector<std::tuple<size_t, uint64_t>> smerOrder;
		smerOrder.resize(0);
		std::vector<size_t> positions;
		//if (endSmers.size() > 0)
		//{
		//	findSyncmerPositions(seq, kmerSize, kmerSize - windowSize + 1, smerOrder, [this](uint64_t hash) { return endSmers[hash % endSmers.size()]; }, [this, &positions](size_t pos)
		//	{
		//		assert(positions.size() == 0 || pos > positions.back());
		//		assert(positions.size() == 0 || pos - positions.back() <= windowSize);
		//		positions.push_back(pos);
		//	});
		//}
		//else
		//{
			findSyncmerPositions(seq, kmerSize, windowSize, smerSize, smerOrder, minimizers, minimizersPos); //kmerSize, kmerSize - windowSize + 1

			//{
				// assert(positions.size() == 0 || pos > positions.back());
				// assert(positions.size() == 0 || pos - positions.back() <= windowSize);
			//	positions.push_back(pos);
			//});
		//}
		//callback(read, seq, poses, rawSeq, positions);

		//cout << endl << endl << endl;
		for(u_int64_t pos : minimizersPos){
			//cout << pos << " " << (readSeq.size() - (pos + windowSize + smerSize - 1))<< endl;
		}
		
		string readSeqRev = readSeq;
		toReverseComplement(readSeqRev);

		for(u_int64_t pos : minimizersPos){

			string minimizerSequence = readSeq.substr(pos, kmerSize);
			size_t revPos = (readSeq.size() - (pos + windowSize + smerSize - 1));
			string minimizerSequenceRev = readSeqRev.substr(revPos, kmerSize);

			//cout << minimizerSequence << " " << minimizerSequenceRev << endl;

			u_int64_t h1 = hashLala(minimizerSequence);
			u_int64_t h2 = hashLala(minimizerSequenceRev);

			if(h1 < h2){
				minimizers.push_back(h1);
				//cout << h1 << endl;
			}
			else{
				minimizers.push_back(h2);
				//cout << h2 << endl;
			}

			//	VectorView<CharType> minimizerSequence { seq, pos, pos + kmerSize };
			//	size_t revPos = seq.size() - (pos + kmerSize);
			//	VectorView<CharType> revMinimizerSequence { revSeq, revPos, revPos + kmerSize };
			//	HashType fwHash = hash(minimizerSequence, revMinimizerSequence);



		}
		

	}


	static void toReverseComplement(string& seq) {


		size_t size = seq.size();

		std::reverse(seq.begin(), seq.end());
		for (std::size_t i = 0; i < size; ++i){
			
			switch (seq[i]){
			case 'A':
				seq[i] = 'T';
				break;    
			case 'C':
				seq[i] = 'G';
				break;
			case 'G':
				seq[i] = 'C';
				break;
			case 'T':
				seq[i] = 'A';
				break;
			//default:
			//	seq[i] = seq[i];
			//	break;
			}
		}


	}

	u_int64_t hashLala(const string& seq){
		MurmurHash3_x64_128 (seq.c_str(), seq.size(), _seed, _hash_otpt);
		u_int64_t kemrHashed = _hash_otpt[0];

		return kemrHashed;
	}

	//template <typename F, typename EdgeCheckFunction>
	void findSyncmerPositions(const vector<u_int16_t>& sequence, size_t kmerSize, size_t windowSize, size_t smerSize, std::vector<std::tuple<size_t, uint64_t>>& smerOrder, vector<u_int64_t>& minimizers, vector<u_int64_t>& minimizersPos)
	{

		//vector<u_int64_t> hashes(sequence.size(), -1);

		//if (sequence.size() < kmerSize) return;
		//assert(smerSize <= kmerSize);
		//size_t windowSize = kmerSize - smerSize + 1;
		//assert(windowSize >= 1);


		MBG::FastHasher fwkmerHasher { smerSize };
		for (size_t i = 0; i < smerSize; i++)
		{
			fwkmerHasher.addChar(sequence[i]);
		}
		auto thisHash = fwkmerHasher.hash();
		//if (endSmer(thisHash)) thisHash = 0;
		smerOrder.emplace_back(0, thisHash);
		for (size_t i = 1; i < windowSize; i++)
		{
			size_t seqPos = smerSize+i-1;
			fwkmerHasher.addChar(sequence[seqPos]);
			fwkmerHasher.removeChar(sequence[seqPos-smerSize]);
			uint64_t hash = fwkmerHasher.hash();
			//hashes[i] = hash;
			//if (endSmer(hash)) hash = 0;
			while (smerOrder.size() > 0 && std::get<1>(smerOrder.back()) > hash) smerOrder.pop_back();
			smerOrder.emplace_back(i, hash);
		}
		if ((std::get<0>(smerOrder.front()) == 0) || (std::get<1>(smerOrder.back()) == std::get<1>(smerOrder.front()) && std::get<0>(smerOrder.back()) == windowSize-1))
		{
			minimizersPos.push_back(0);
			//minimizers.push_back(0);
			//cout << 0 << endl;
			//callback(0);
		}
		for (size_t i = windowSize; smerSize+i-1 < sequence.size(); i++)
		{
			size_t seqPos = smerSize+i-1;
			fwkmerHasher.addChar(sequence[seqPos]);
			fwkmerHasher.removeChar(sequence[seqPos-smerSize]);
			uint64_t hash = fwkmerHasher.hash();
			//hashes[i] = hash;
			//if (endSmer(hash)) hash = 0;
			// even though pop_front is used it turns out std::vector is faster than std::deque ?!
			// because pop_front is O(w), but it is only called in O(1/w) fraction of loops
			// so the performace penalty of pop_front does not scale with w!
			// and std::vector's speed in normal, non-popfront operation outweighs the slow pop_front
			while (smerOrder.size() > 0 && std::get<0>(smerOrder.front()) <= i - windowSize) smerOrder.erase(smerOrder.begin());
			while (smerOrder.size() > 0 && std::get<1>(smerOrder.back()) > hash) smerOrder.pop_back();
			smerOrder.emplace_back(i, hash);
			if ((std::get<0>(smerOrder.front()) == i-windowSize+1) || (std::get<1>(smerOrder.back()) == std::get<1>(smerOrder.front()) && std::get<0>(smerOrder.back()) == i))
			{
				minimizersPos.push_back(i-windowSize+1);
				//minimizers.push_back(hashes[i]);
				//cout << i-windowSize+1 << " " << hashes[i] << endl;
				//callback(i-windowSize+1);
			}
		}

	}
	*/

};






class MinimizerParser128{

public:

	u_int16_t _minimizerSize;
	KmerModel128* _kmerModel;
	u_int32_t _seed;
	u_int128_t _hash_otpt;
	double _minimizerBound;

	MinimizerParser128(u_int16_t minimizerSize, double minimizerDensity){
		_minimizerSize = minimizerSize;
		_kmerModel = new KmerModel128(_minimizerSize);
		_seed = 42;
		//_hash_otpt = new u_int64_t[2];
		u_int128_t maxHashValue = -1;
		_minimizerBound = minimizerDensity * maxHashValue;
	}

	~MinimizerParser128(){
		//delete[] _hash_otpt;
	}

	void parse(const string& seq, vector<u_int128_t>& minimizers, vector<u_int64_t>& minimizersPos){

		minimizers.clear();
		minimizersPos.clear();

		//0 0
		//100000 3278692
		//200000 6550207
		//300000 9827077

		vector<u_int128_t> kmers;
		_kmerModel->iterate(seq.c_str(), seq.size(), kmers);

		if(kmers.size() == 0) return;

		for(u_int64_t pos=1; pos<kmers.size()-1; pos++){

			//cout << _kmerModel->toString(kmers[pos]) << endl;
			//cout << itKmer->value().getVal() << endl;
			//cout << itKmer->value() << endl;
			//cout << pos << " " << (rleSequence.size()-_minimizerSize) << endl;

			/*
			if(pos == 0){
				//cout << "lala1" << endl;
				//pos += 1;
				continue;
			}
			else if(pos == kmers.size()-1; //seq.size()-_minimizerSize){
				//cout << "lala2" << endl;
				continue;
			}*/

			//if(!itKmer->value().isValid()) continue;
			//kmer_type kmerMin = min(itKmer->value(), revcomp(itKmer->value(), _kmerSize));
			//if(lala < 100 ) cout << model.toString(itKmer->value()) << endl;
			//lala += 1;
			u_int128_t kmerValue = kmers[pos];
			MurmurHash3_x64_128 ((const char*)&kmerValue, sizeof(kmerValue), _seed, &_hash_otpt);
			u_int128_t minimizer = _hash_otpt;



			//if(minimizerCounts[minimizer] > 1000) cout << minimizer << endl;
			//double kmerHashed_norm = ((double) minimizer) / maxHashValue;
			if(minimizer < _minimizerBound){//_minimizerDensity){


				minimizers.push_back(minimizer);
				minimizersPos.push_back(pos);

				
				//minimizerCounts[minimizer] += 1;
				//nbMinimizers += 1;
				//"reprise: systeme de functor pour les reads et kmers?"
				

			}
			
			//pos += 1;
		}

	}

};


#endif