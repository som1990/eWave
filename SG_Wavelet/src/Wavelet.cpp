#include "Wavelet.h"
#include <cassert>

Wavelet::~Wavelet()
{

}

void Wavelet::resize(int xSize, int ySize, int thetaSize, int kSize)
{
	m_dimensions = std::array<int, 4> {xSize, ySize, thetaSize, kSize};
	m_data.resize(xSize*ySize*thetaSize*kSize);
	std::cout << "m_data size (xSize * ySize * thetaSize * kSize ) :  " << xSize << " * " << ySize << " * " << thetaSize << " * " << kSize << " = " << m_data.size() << std::endl;
	std::cout << "RESIZE COMPLETE!" << std::endl;
}

float & Wavelet::operator()(int i0, int i1, int i2, int i3)
{
	assert(i0 >= 0 && i0 < dimSize(dimOption::X) && i1 >= 0 && i1 < dimSize(dimOption::Y) && i2 >= 0 &&
		i2 < dimSize(dimOption::THETA) && i3 >= 0 && i3 < dimSize(dimOption::K));
	return m_data[i0 + i1 * dimSize(0) + i2 * dimSize(0)*dimSize(1) + i3 * dimSize(0)*dimSize(1)*dimSize(2)];
}

float const & Wavelet::operator()(int i0, int i1, float i2, int i3) const
{
	return m_data[i0 + i1 * dimSize(0) + i2 * dimSize(0)*dimSize(1) + i3 * dimSize(0)*dimSize(1)*dimSize(2)];
}

int Wavelet::dimSize(const int option) const
{
	return m_dimensions[option];
}
