
#ifndef SurfaceWavelet_h__
#define SurfaceWavelet_h__


#include "SurfaceWavelet_PCH.h"
#include "Wavelet.h"

class SurfaceWavelet 
{
public:
	SurfaceWavelet();
	SurfaceWavelet(int gridX, int gridY);
	~SurfaceWavelet();

	//Returns height for main.cpp to access
	const float get_height(const int &i, const int &j) const;
	const float get_height(const int &index) const;
	void initFields(int nGridX, int nGridY, int thetaSamples, int kSamples);		
	void propogate(float *&source_height, float dt);

private:
	int simGridX, simGridY;
	float m_dx = 1.0f, m_dy = 1.0f;
	
	float m_dt;
	float *m_height;
	float *res_height;

	Wavelet m_Amplitude;

	void timeStep();
	void addAmplitude(float *&source_height);
	void genHeightMap(float *&out_height);
	void precomputeProfileBuffer();

};


#endif // SurfaceWavelet_h__