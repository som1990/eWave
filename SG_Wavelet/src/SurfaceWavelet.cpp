#include "SurfaceWavelet.h"

#include <iostream>
#include <algorithm>
using namespace std;

#define PI 3.14159265359

SurfaceWavelet::SurfaceWavelet() {

}

SurfaceWavelet::SurfaceWavelet(int gridX, int gridY)
	: simGridX(gridX/m_dx), simGridY(gridY/m_dy) {
}

SurfaceWavelet::~SurfaceWavelet() {
	free(m_height);
	free(res_height);
}

const float SurfaceWavelet::get_height(const int &index) const
{
	float res = res_height[index];
	return res;
}

const float SurfaceWavelet::get_height(const int &i, const int &j) const
{
	int index = i + simGridX * j;
	float res = res_height[index];
	//float res = m_height[index];
	return res;
}

void SurfaceWavelet::addAmplitude(float *&source_height)
{
	//std::cout << "Adding Amplitude from paint" << std::endl;
	 
	for (int i = 0; i < simGridX; i++)
	{
		#pragma omp parallel for  
		for (int j = 0; j < simGridY; j++)
		{
			int index = i + j * simGridX;
			//Adding a equal force around all angles
		
			for (int theta = 0; theta < m_Amplitude.dimSize(dimOption::THETA); theta++)
			{
				m_Amplitude(i, j, theta, 0) += source_height[index]*m_dt;
			}
			
		}
	}
	//std::cout << "Finished adding Amplitude" << std::endl;
}

void SurfaceWavelet::genHeightMap(float *&out_height)
{
	
	for (int i = 0; i < simGridX; i++)
	{
		#pragma omp parallel for 
		for (int j = 0; j < simGridY; j++)
		{
			int index = i + j * simGridX;
			float h = 0;
			Vec2 pos = Vec2{ float(i),float(j) };
			int numDir = m_Amplitude.dimSize(dimOption::THETA);
			int numAngleSamples = 4 * numDir; 
			float dx = numDir * 2 * PI / numAngleSamples;
			
			for (int c = 0; c < numAngleSamples; c += 1)
			{
				float angle = 2.0 * PI * c / (numAngleSamples);
				Vec2 k_dir = Vec2{cos(angle), sin(angle)};
				float p = k_dir * pos ;
				
				for (int k = 0; k < m_Amplitude.dimSize(dimOption::K); k++)
				{
					h += dx*m_Amplitude( int(pos[dimOption::X]), int(pos[dimOption::Y]), int(c/4.0) , k ); // Needs a profile buffer and Interpolation
				}
			}
			out_height[index] += h;
		}
	}
}

void SurfaceWavelet::initFields(int nGridX, int nGridY, int thetaSamples, int kSamples)
{
	this->simGridX = nGridX / m_dx;
	this->simGridY = nGridY / m_dy;

	int size = simGridX * simGridY;
	
	cout << "simGridX: " << simGridX << " SimGrid Y: " << simGridY << endl;
	m_height = (float*)calloc(size, sizeof(float)) ;
//	m_height = new float[size];
	res_height = (float*)calloc(size, sizeof(float));
//	this->res_height = new float[size];
//	res_height = (float*)malloc(sizeof(float)*size);
	cout << "Array Size: " << size << endl;
	m_Amplitude.resize(simGridX, simGridY, thetaSamples, kSamples);
	cout << "4D Grid Resize Complete!" << endl;

}

void SurfaceWavelet::timeStep()
{
	
	
}

void SurfaceWavelet::propogate(float *&source_height, float dt)
{
	this->m_dt = dt;
	addAmplitude(source_height);
	timeStep();


	genHeightMap (m_height);

	res_height = m_height;


}
