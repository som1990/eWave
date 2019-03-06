#pragma once
#ifndef ProfileBuffer_h__
#define ProfileBuffer_h__

#include <array>
#include <vector>

#include "math/Math.h"

class ProfileBuffer {
public:
	ProfileBuffer() {}


	template <typename Spectrum>
	void precompute(Spectrum &spectrum, float time, float zeta_min, float zeta_max,
		int resolution = 2048, int periodicity = 2, int integration_nodes = 100)
	{
		m_period = periodicity * pow(2, zeta_max);

#pragma omp parallel for
		for (int i = 0; i < resolution; i++)
		{
			constexpr float tau = 6.28318530718;
			float p = (i * m_period) / resolution;

			m_data[i] = integrate(integration_nodes, zeta_min, zeta_max, [&](float zeta) {

				float waveLength = pow(2, zeta); //WaveLength from 1 to 1024 in our case
				float waveNumber = tau / waveLength; // lambda = 2*PI / WaveLength 

				float phase1 = waveNumber * p - dispersionRelation(waveNumber) * time;
				float phase2 = waveNumber * (p - m_period) -
					           dispersionRelation(waveNumber) * time;

				float weight1 = p / m_period;
				float weight2 = 1 - weight1;

				//periodic
				return waveLength * spectrum(zeta) * 
					( weight1 * trochoidalWave(phase1, waveNumber) +
				      weight2 * trochoidalWave(phase2, waveNumber));

			});
		}
	}


	std::array<float, 4> operator() (float p) const;
private:

	float dispersionRelation(float k) const;

	std::array<float, 4> trochoidalWave(float phase, float knum) const;

	float cubic_bump(float x) const;

public:
	float m_period;

	std::vector<std::array<float, 4>> m_data;
};
#endif // ProfileBuffer_h__