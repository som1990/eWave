#include "ProfileBuffer.h"


std::array<float, 4> ProfileBuffer::operator() (float p) const {

	const int N = m_data.size();
	//Makes sure the values outside the buffer wrap around
	auto extended_buffer = [=](int i)->std::array<float, 4> {
		return m_data[pos_modulo(i, N)];
	};

	//Linear interpolation performed
	auto interpolated_buffer = LinearInterpolation(extended_buffer);

	//rescale p to interval [0,1)
	return interpolated_buffer(N*p / m_period); 
}

float ProfileBuffer::dispersionRelation(float k) const
{
	const float GRAVITY = 9.81;
	return sqrt(k*GRAVITY);
}

std::array<float, 4> ProfileBuffer::trochoidalWave(float phase, float knum) const
{
	float s = sin(phase);
	float c = cos(phase);
	return std::array<float, 4>{-s, c, -knum * c, -knum * s};
}

//bi-cubic interpolation.
float ProfileBuffer::cubic_bump(float x) const
{
	if (abs(x) >= 1)
		return 0.0f;
	else
		return x * x *(2 * abs(x) - 3) + 1;
}

