#pragma once
#ifndef ParallelStream_h__
#define ParallelStream_h__

#include <sstream>

class ParallelStream {
	std::ostringstream stdStream;

public:
	ParallelStream() {}
	template <class T>
	ParallelStream& operator<<(const T& inData) {
		stdStream << inData;
		return *this;
	}
	std::string toString() const {
		return stdStream.str();
	}
};
#endif // ParallelStream_h__