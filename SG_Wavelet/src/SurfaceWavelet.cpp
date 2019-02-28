#include "SurfaceWavelet.h"

#include <iostream>
#include <algorithm>
using namespace std;

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
	//float res = res_height[index];
	float res = m_height[index];
	return res;
}

void SurfaceWavelet::initFields(int nGridX, int nGridY)
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
}

void SurfaceWavelet::propogate(float *&source_height, float dt)
{
	this->m_dt = dt;
	
	res_height = m_height;


}


