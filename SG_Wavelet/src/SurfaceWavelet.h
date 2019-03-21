
#ifndef SurfaceWavelet_h__
#define SurfaceWavelet_h__


#include "SurfaceWavelet_PCH.h"
#include "Wavelet.h"
#include "ProfileBuffer.h"
#include "Spectrum.h"

class SurfaceWavelet 
{
public:
	typedef std::array<int, 4> Idx;
	
	struct Settings {
		float xSize = 512;
		float ySize = 512;

		float max_zeta = 0.01;

		float min_zeta = 10;

		float n_x = 100;
		float n_y = 100;


		int n_theta = 16;

		int n_zeta = 1;

		float initial_time = 100;

		enum SpectrumType {
			PiersonMoskowitz = 0,
			Phillips

		} spectrumType = PiersonMoskowitz;
		float wind_direction = 0;
		float wind_speed = 1;
	};


	SurfaceWavelet(Settings s);
	SurfaceWavelet(int gridX, int gridY);
	~SurfaceWavelet();

	//Returns height for main.cpp to access
	const float get_height(const int &i, const int &j) const;
	const float get_height(const int &index) const;
	void initFields(int nGridX, int nGridY, int thetaSamples, int kSamples);		
	void addPointAmplitude(float i, float j, float val);
	void propogate(float *&source_height, float dt);
	

private:
	void timeStep();

	void AdvectionStep(const float dt);
	void WaveDiffusionStep(const float dt);

	Vec4 boundaryConditions(const Vec4& pos4) const;

	void addAmplitude(float *&source_height);
	void genHeightMap(float *&out_height);
	void precomputeProfileBuffer();
	void precomputeGroupSpeeds();

	float idxToPos(const int idx, const int dim) const;
	Vec4 idxToPos(std::array<int, 4> idx) const;

	float posToGrid(float pos, int dim) const;
	Vec4 posToGrid(Vec4 pos4) const;

	int posToIdx(float pos, int dim) const;
	Idx posToIdx(Vec4 pos4) const;
	
	auto correctFringeValues() const;

	auto interpolateAmplitude() const;
	float getAmplitude(Vec4 pos4) const;

	int gridDim(const int dim) const;
	float defaultAmplitude(const int itheta, const int izeta) const;

private:
	int simGridX, simGridY;
	float m_dx = 1.0f, m_dy = 1.0f;
	
	float m_dt, m_time = 0;
	float *m_height;
	float *res_height;

	Wavelet m_Amplitude, new_Amplitude;

	std::array<float, 4> m_xmin;
	std::array<float, 4> m_xmax;
	std::array<float, 4> m_delta;
	std::array<float, 4> m_idx;

	std::vector<float> m_groupSpeeds;
	std::vector<ProfileBuffer> m_profileBuffer;
	Spectrum m_spectrum;
	
};


#endif // SurfaceWavelet_h__