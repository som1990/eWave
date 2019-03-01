#pragma once
#ifndef Wavelet_h__
#define Wavelet_h__

#include "SurfaceWavelet_PCH.h"

class Wavelet {
public:
	Wavelet()
		:m_data(0), m_dimensions({ 0,0,0,0 }) {}
	~Wavelet();

	void resize(int xSize, int ySize, int thetaSize, int kSize);

	float &operator ()(int i0, int i1, int i2, int i3);
	float const &operator ()(int i0, int i1, float i2, int i3) const;
	//1 - xSize 2 - ySize 3 - thetaSize 4 - kSize
	int dimSize(const int option) const;

private:
	std::vector<float> m_data;
	std::array<int, 4> m_dimensions; 


};


#endif // Wavelet_h__