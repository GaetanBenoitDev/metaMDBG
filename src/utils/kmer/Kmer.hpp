

#ifndef MDBG_METAG_KMER
#define MDBG_METAG_KMER


#include <vector>
#include <string>
#include <cmath>

using namespace std;

typedef unsigned __int128 u_int128_t;



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

	/** Returns the value of the kmer.
	 * \return the kmer value as a Type object. */
	const u_int64_t& value  () const { return table[(int)choice];   }

	/** Returns the value of the kmer.
	 * \param[in] which: forward or reverse strand
	 * \return the kmer value as a Type object. */
	const u_int64_t& value  (int which) const { return table[which];   }

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

	bool iterate (const char* seq, size_t length, vector<u_int64_t>& tmpKmers) const{

		int32_t nbKmers = length - _kmerSize + 1;
		if (nbKmers <= 0)  { return false; }

		tmpKmers.resize(nbKmers);

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
















class MinimizerParser{

public:

	u_int16_t _minimizerSize;
	KmerModel* _kmerModel;
	u_int32_t _seed;
	u_int64_t* _hash_otpt;
	double _minimizerBound;
	size_t _trimBps;

	MinimizerParser(u_int16_t minimizerSize, double minimizerDensity){
		_minimizerSize = minimizerSize;
		_kmerModel = new KmerModel(_minimizerSize);
		_seed = 42;
		_hash_otpt = new u_int64_t[2];
		u_int64_t maxHashValue = -1;
		_minimizerBound = minimizerDensity * maxHashValue;
		_trimBps = 1; 
	}

	~MinimizerParser(){
		delete[] _hash_otpt;
	}

	void parse(const string& seq, vector<u_int64_t>& minimizers, vector<u_int64_t>& minimizersPos){

		minimizers.clear();
		minimizersPos.clear();

		//parse_windowed(seq, minimizers, minimizersPos, 200);
		//return;
		//0 0
		//100000 3278692
		//200000 6550207
		//300000 9827077

		vector<u_int64_t> kmers;
		_kmerModel->iterate(seq.c_str(), seq.size(), kmers);

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
			MurmurHash3_x64_128 ((const char*)&kmerValue, sizeof(kmerValue), _seed, _hash_otpt);
			u_int64_t minimizer = _hash_otpt[0];



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