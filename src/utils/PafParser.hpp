
#ifndef MDBG_METAG_PAFPARSER
#define MDBG_METAG_PAFPARSER

#include <string>
#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include "zlib.h"
#include <cstdlib>
#include <functional>
using namespace std;



class PafParser {  

public:

	gzFile _file;

	vector<char> buffer_;
	uint32_t buffer_ptr_;
	uint32_t buffer_bytes_;


	PafParser(const string& filename){

		_file = gzopen(filename.c_str(), "rb");
		buffer_ptr_ = 0;
		buffer_bytes_ = 0;
		buffer_ = vector<char>(65536, 0);  // 64 kB

	}

	~PafParser(){
		//cout << "close gz file)" << endl;
		gzclose(_file);
	}


	const std::vector<char>& buffer() const {
		return buffer_;
	}

	std::uint32_t buffer_ptr() const {
		return buffer_ptr_;
	}

	std::uint32_t buffer_bytes() const {
		return buffer_bytes_;
	}

	bool Read() {
		buffer_ptr_ = 0;
		buffer_bytes_ = gzread(_file, buffer_.data(), buffer_.size());
		return buffer_bytes_ < buffer_.size();
	}

	void parse(const std::function<void(string)>& fun){

	
		bool is_eof = false;

		string line = "";

		while (true) {
			
			auto buffer_ptr = this->buffer_ptr();
			for (; buffer_ptr < this->buffer_bytes(); ++buffer_ptr) {

				//cout << buffer_ptr << endl;
				auto c = this->buffer()[buffer_ptr];
				line += c;
				//cout << c << endl;
				if (c == '\n') {
					line.pop_back();
					//cout << line << endl;
					
					fun(line);
					line = "";

				}
			}

			if (is_eof) {
				break;
			}
			is_eof = this->Read();
		}

	}

};





#endif 