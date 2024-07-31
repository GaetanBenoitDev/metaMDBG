#ifndef VectorView_h
#define VectorView_h

#include <string_view>
#include <vector>

namespace MBG
{

template <typename T>
class VectorView
{
public:
	VectorView(const std::vector<T>& data, size_t startpos, size_t endpos) :
	data(data),
	startpos(startpos),
	endpos(endpos)
	{
	}
	size_t size() const
	{
		return endpos-startpos;
	}
	const T& operator[](size_t pos) const
	{
		return data[startpos+pos];
	}
	const T* begin()
	{
		return data.data() + startpos;
	}
	const T* end()
	{
		return data.data() + endpos;
	}
	std::string_view getView()
	{
		return std::string_view { (char*)(data.data() + startpos), size() * sizeof(T) };
	}
private:
	const std::vector<T>& data;
	size_t startpos;
	size_t endpos;
	friend unsigned __int128 hash(VectorView<uint16_t> sequence);
	friend unsigned __int128 hash(VectorView<uint16_t> sequence, VectorView<uint16_t> reverseSequence);
};

}

#endif
