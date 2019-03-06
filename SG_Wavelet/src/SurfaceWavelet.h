
#ifndef SurfaceWavelet_h__
#define SurfaceWavelet_h__


#include "SurfaceWavelet_PCH.h"
#include "Wavelet.h"
#include "ProfileBuffer.h"

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
	void timeStep();
	void addAmplitude(float *&source_height);
	void genHeightMap(float *&out_height);
	void precomputeProfileBuffer();

	float idxToPos(const int idx, const int dim) const;
	Vec4 idxToPos(std::array<int, 4> idx) const;

private:
	int simGridX, simGridY;
	float m_dx = 1.0f, m_dy = 1.0f;
	
	float m_dt;
	float *m_height;
	float *res_height;

	Wavelet m_Amplitude;

	std::array<float, 4> m_xmin;
	std::array<float, 4> m_xmax;
	std::array<float, 4> m_delta;
	std::array<float, 4> m_idx;

	std::vector<float> m_groupSpeeds;
	std::vector<ProfileBuffer> m_profileBuffer;

	
	

};


#endif // SurfaceWavelet_h__