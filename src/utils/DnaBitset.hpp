

#ifndef MDBG_METAG_DNABITSET
#define MDBG_METAG_DNABITSET

#include <string.h>
#include <unistd.h>
#include <string>
#define BASE_MASK 0x3 /* binary: 11 */
using namespace std;

/* useful constants */
enum
{
    BASE_A = 0x0, /* binary: 00 */
    BASE_C = 0x1, /* binary: 01 */
    BASE_G = 0x2, /* binary: 10 */
    BASE_T = 0x3, /* binary: 11 */
};
/*
class DnaBitset
{
public:



    DnaBitset(const string& sequence){

        const char* dna_str = sequence.c_str();
        const size_t dna_len = sequence.size();
        m_len = dna_len;
 
        size_t dna_bytes = (dna_len / 4) + (dna_len % 4 != 0);
        _bitsetSize = dna_bytes;

        m_data = new uint8_t[dna_bytes];
 
        //std::memset(m_data, 0, dna_bytes);
        memset(m_data, 0, dna_bytes); 

        for (size_t i = 0; i < dna_len; ++i)
        {
            uint8_t shift = 6 - 2 * (i % 4);
 
            switch (dna_str[i])
            {
                case 'A':
                    m_data[i / 4] |= BASE_A << shift;
                    break;
                case 'C':
                    m_data[i / 4] |= BASE_C << shift;
                    break;
                case 'G':
                    m_data[i / 4] |= BASE_G << shift;
                    break;
                case 'T':
                    m_data[i / 4] |= BASE_T << shift;
                    break;
                //default:
                  //  throw std::invalid_argument("invalid DNA base");
            }
 
            shift = (shift == 0) ? 6 : shift - 2;
        }
    }
 
     DnaBitset(uint8_t* dna_byte, u_int32_t bitsetSize, u_int32_t sequenceSize){
        _bitsetSize = bitsetSize;
        m_data = dna_byte;
        m_len = sequenceSize;
    }

    ~DnaBitset()
    {
        delete[] m_data;
    }
 
    char* to_string() const
    {
        char* dna_str = new char[m_len + 1];
 
        for (size_t i = 0; i < m_len; ++i)
        {
            uint8_t shift = 6 - 2 * (i % 4);
            uint8_t mask = BASE_MASK << shift;
 
            uint8_t base = (m_data[i / 4] & mask) >> shift;
 
            switch (base)
            {
                case BASE_A:
                    dna_str[i] = 'A';
                    break;
                case BASE_C:
                    dna_str[i] = 'C';
                    break;
                case BASE_G:
                    dna_str[i] = 'G';
                    break;
                case BASE_T:
                    dna_str[i] = 'T';
                    break;
                //default:
                  //  throw std::runtime_error("invalid DNA base");
            }
        }
 
        dna_str[m_len] = '\0';
        return dna_str;
    }
 
    u_int32_t m_len;
    uint8_t* m_data;
    u_int32_t _bitsetSize;
    
};
*/
class DnaBitset2
{
public:

    DnaBitset2(){
        _m_len = 0;
    }
    
    DnaBitset2(const string& sequence){

        _m_len = 0;
        
        const char* dna_str = sequence.c_str();
        const size_t dna_len = sequence.size();
        _m_len = dna_len;
 
        /* number of bytes necessary to store dna_str as a bitset */
        size_t dna_bytes = getBinarySize(dna_len); //(dna_len / 4) + (dna_len % 4 != 0);

        //m_data = new uint8_t[dna_bytes];
        _m_data.resize(dna_bytes, 0);

        //std::memset(m_data, 0, dna_bytes);
        //memset(m_data, 0, dna_bytes); 

        /* for each base of the DNA sequence */
        for (size_t i = 0; i < dna_len; ++i)
        {
            uint8_t shift = 6 - 2 * (i % 4);
 
            switch (dna_str[i])
            {
                case 'A':
                    _m_data[i / 4] |= BASE_A << shift;
                    break;
                case 'C':
                    _m_data[i / 4] |= BASE_C << shift;
                    break;
                case 'G':
                    _m_data[i / 4] |= BASE_G << shift;
                    break;
                case 'T':
                    _m_data[i / 4] |= BASE_T << shift;
                    break;
                //default:
                  //  throw std::invalid_argument("invalid DNA base");
            }
 
            shift = (shift == 0) ? 6 : shift - 2;
        }

        //if(sequence != to_string()){
        //    cout << sequence << endl;
        //    cout << to_string() << endl;
        //    cout << sequence.size() << endl;
        //    cout << to_string().size() << endl;
        //    cout << "derp" << endl;
        //    getchar();
        //}

        //#pragma omp critical
        //{
        //    cout << sequence << endl;
        //    cout << to_string() << endl;
        //}
    }
 
    //DnaBitset2(uint8_t* dna_byte, u_int32_t sequenceSize){
    //    m_data = dna_byte;
    //    m_len = sequenceSize;
    //}

    ~DnaBitset2(){
        _m_data.clear();
        //delete[] m_data;
    }
 
    string to_string() const{
        
        if(_m_len == 0) return "";

        string dna_str(_m_len, 'A');
        //char* dna_str = new char[m_len + 1];
 
        /* for each base of the DNA sequence */
        for (size_t i = 0; i < _m_len; ++i)
        {
            uint8_t shift = 6 - 2 * (i % 4);
            uint8_t mask = BASE_MASK << shift;
 
            /* get the i-th DNA base */
            uint8_t base = (_m_data[i / 4] & mask) >> shift;
 
            switch (base)
            {
                case BASE_A:
                    dna_str[i] = 'A';
                    break;
                case BASE_C:
                    dna_str[i] = 'C';
                    break;
                case BASE_G:
                    dna_str[i] = 'G';
                    break;
                case BASE_T:
                    dna_str[i] = 'T';
                    break;
                //default:
                  //  throw std::runtime_error("invalid DNA base");
            }
        }
 
        //dna_str[_m_len] = '\0';
        return dna_str;
    }
 
    size_t getBinarySize(const size_t dna_len){
        size_t dna_bytes = (dna_len / 4) + (dna_len % 4 != 0);
        return dna_bytes;
    }

    u_int32_t _m_len;
    vector<uint8_t> _m_data;
    
};

#endif