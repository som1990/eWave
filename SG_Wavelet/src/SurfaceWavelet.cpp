#include "SurfaceWavelet.h"

#include "ParallelStream.h"
#include <iostream>
#include <algorithm>
using namespace std;

#define PI 3.14159265359

static void printVector(std::vector<float> &arr, string name )
{
	int size = arr.size();
	cout << name << ": " ;
	for (int i = 0; i < size; i++)
	{
		cout << "     " << arr[i] << ", ";
	}
	cout << endl;
}

SurfaceWavelet::SurfaceWavelet(Settings s) : m_spectrum(10) {
	m_Amplitude.resize(s.n_x, s.n_y, s.n_theta, s.n_zeta);
	new_Amplitude.resize(s.n_x, s.n_y, s.n_theta, s.n_zeta);
	cout << "4D Grid Resize Complete!" << endl;
	float zeta_min = m_spectrum.minZeta();
	float zeta_max = m_spectrum.maxZeta();
	simGridX = s.n_x;
	simGridY = s.n_y;

	m_xmin = { 0.0,0.0,0.0, zeta_min };
	m_xmax = { s.xSize, s.ySize, 2 * PI, zeta_max };

	for (int i = 0; i < 4; i++)
	{
		m_delta[i] = (m_xmax[i] - m_xmin[i]) / m_Amplitude.dimSize(i);
		m_idx[i] = 1.0 / m_delta[i];
	}
	m_spectrum.setWindDir(s.wind_direction);
	m_spectrum.setWindSpeed(s.wind_speed);
	m_spectrum.SpectrumType = s.spectrumType;
	m_time = s.initial_time;
	m_profileBuffer.resize(s.n_zeta);
	precomputeGroupSpeeds();
	printVector(m_groupSpeeds, "GroupSpeed");

	int size = simGridX * simGridY;

	cout << "simGridX: " << simGridX << " SimGrid Y: " << simGridY << endl;
	m_height = (float*)calloc(size, sizeof(float));
	res_height = (float*)calloc(size, sizeof(float));
	
	cout << "Array Size: " << size << endl;
	

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

static float mag(Vec2 a) {
	return (sqrt(a[0] * a[0] + a[1] * a[1]));
}

float sdBox(Vec2 &p, Vec2 &b)
{
	Vec2 d = Vec2{ fabs(p[0]),fabs(p[1]) } - b;
	return mag(Vec2{ max(p[0],0.0f),max(p[1],0.0f) }) + min(max(d[0], d[1]), 0.0f);
}

Vec4 SurfaceWavelet::boundaryConditions(const Vec4 &pos4) const
{
	return pos4;
	Vec4 boundaryCenter = (m_xmax - m_xmin) * 0.5;
	Vec2 cen = Vec2{ boundaryCenter[0], boundaryCenter[1] };
	Vec2 pos = Vec2{ pos4[dimOption::X], pos4[dimOption::Y] } - cen ;

	

	Vec2 b = Vec2{ m_xmax[0], m_xmax[1] };
	//sim boundary level set
	float ls = sdBox(pos, b);
	
	if (ls >= 0.0)// no reflection if point is inside the sim grid.
	{
		return pos4;
	}
	//Generating Boundary conditions based on levelsets
	cout << "pos: " << pos4[0] << ", " << pos4[1] << endl;
	cout << "ls: " << ls << endl;
	return Vec4 {0, 0, 0, 0};

	/*
	//In our case the boundary are the edges of window
	Vec2 p_min = Vec2{ m_xmin[dimOption::X],m_xmin[dimOption::Y] };
	Vec2 p_max = Vec2{ m_xmax[dimOption::X],m_xmax[dimOption::Y] };

	Vec2 lowCheck = pos - p_min;
	Vec2 highCheck = p_max - pos;

	//Within the boundaries of window.
	if (lowCheck[0] >= 0 && highCheck[0] > 0 && lowCheck[1] >= 0 && highCheck[1] > 0)
	{
		return pos4;
	}

	Vec2 n{ 0,0 };
	
	if (lowCheck[0] < 0) n[0] = 1.0f;
	else if (highCheck[0] <= 0) n[0] = -1.0f;
	if (lowCheck[1] < 0) n[1] = 1.0f;
	else if (highCheck[1] <= 0) n[1] = -1.0f;
	*/
	


}

void SurfaceWavelet::addPointAmplitude(float x, float y, float val)
{
	int ix = posToIdx(x, dimOption::X);
	int iy = posToIdx(y, dimOption::Y);

	if (ix >= 0 && ix < gridDim(dimOption::X) && iy >= 0 && iy < gridDim(dimOption::Y)) {
		for (int theta = 0; theta < gridDim(dimOption::K); theta++)
		{
			m_Amplitude(ix, iy, theta, 0) += val;
		}

	}
}

void SurfaceWavelet::addAmplitude(float *&source_height)
{
	//std::cout << "Adding Amplitude from paint" << std::endl;
	 
	for (int i = 1; i < simGridX-1; i++)
	{
		
		#pragma omp parallel for  
		for (int j = 1; j < simGridY-1; j++)
		{
			int index = i + j * simGridX;
			float gradY = (source_height[i + (j + 1)*simGridX] - source_height[i + (j - 1)*simGridX]) * 0.5;
			float gradX = (source_height[(i + 1) + j*simGridX] - source_height[(i - 1) + j*simGridX]) * 0.5;
			//Adding a equal force around all angles
			if (source_height[index] > 0)
			{
				float dir = 0;
				float ngradX = 0;
				float ngradY = 0;
				float magnitude = sqrtf(gradX*gradX + gradY * gradY);
				//int theta = 0;
				if (magnitude > 0)
				{
					ngradX = gradX / magnitude;
					ngradY = gradY / magnitude;
					dir = acos(-ngradX);
					//theta = int(dir / m_delta[dimOption::THETA]);
					
				}
				else
				{
					magnitude = 0.0f;
				}
				
				

				//cout << (ParallelStream() << "\nnGradX, nGradY: " << ngradX << ", " << ngradY).toString() << endl;
				//cout << (ParallelStream() << "\ngradX, gradY: " << gradX << ", " << gradY).toString() << endl;
				//cout << (ParallelStream() << "\nDir, theta, magnitude: " << dir << ", " << theta << ", " << magnitude).toString() << endl;
				
				
				
				for (int theta = 0; theta < m_Amplitude.dimSize(dimOption::THETA); theta++)
				{
					float angle = (2.0*PI*theta) / m_Amplitude.dimSize(dimOption::THETA);
					float dotprod = cos(angle) * -ngradX + sin(angle) * -ngradY;
					//float dotprod = 1.0f;
					//clamp(dotprod, 0.0f, 1.0f);
					//cout << (ParallelStream() << "dotProd: " << dotprod).toString() << endl;
					m_Amplitude(i, j, theta, 0) += source_height[index] * 1000; //* dotprod * 1000;
				}
				
			}
		}
	}
	//std::cout << "Finished adding Amplitude" << std::endl;
}

void SurfaceWavelet::genHeightMap(float *&out_height)
{
	float windAngle = m_spectrum.getWindDir();
	Vec2 windVec = Vec2{ cos(windAngle), sin(windAngle) };
	for (int i = 0; i < simGridX; i++)
	{
		#pragma omp parallel for 
		for (int j = 0; j < simGridY; j++)
		{
			int index = i + j * simGridX;
			float h = 0;
			Vec2 pos = Vec2{ idxToPos(i,dimOption::X),idxToPos(j,dimOption::Y) };
			//cout << (ParallelStream() << "\npos: " << pos[0] << ", " << pos[1]).toString() << endl;
			int numDir = gridDim(dimOption::THETA);
			int numAngleSamples = 4 * numDir; 
			float dx = numDir * 2 * PI / numAngleSamples;
			for (int k = 0; k < gridDim(dimOption::K); k++)
			{
				float zeta = idxToPos(k, dimOption::K);
				auto &profile = m_profileBuffer[k];
				for (int c = 0; c < numAngleSamples; c += 1)
				{
					float angle = 2.0 * PI * c / (numAngleSamples);
					Vec2 k_dir = Vec2{cos(angle), sin(angle)};
					float p = k_dir * pos ;
					
					Vec4 precompProfileSpectrum = profile(p);
					Vec4 wave_data =
						dx * getAmplitude({ pos[dimOption::X],pos[dimOption::Y], angle, zeta }) * precompProfileSpectrum;
					
					if (m_spectrum.SpectrumType == 0)
						h += wave_data[1]; //cos(phase) gives you the vertical displacement
					else if (m_spectrum.SpectrumType == 1) {
						
						float k_dot_w = k_dir[0] * windVec[0] + k_dir[1] * windVec[1];
						k_dot_w = max(k_dot_w, 0.0f);
						h += k_dot_w * wave_data[1];
					}
				}
			}
			out_height[index] = h;
		}
	}
}

void SurfaceWavelet::precomputeProfileBuffer()
{
	for (int izeta = 0; izeta < gridDim(dimOption::K); izeta++)
	{
		float zeta_min = idxToPos(izeta, dimOption::K) - 0.5 * m_delta[dimOption::K];
		float zeta_max = idxToPos(izeta, dimOption::K) + 0.5 * m_delta[dimOption::K];
		m_profileBuffer[izeta].precompute(m_spectrum,m_time, zeta_min, zeta_max );
	}
}

void SurfaceWavelet::precomputeGroupSpeeds()
{
#define GRAV 9.81
	m_groupSpeeds.resize(gridDim(dimOption::K));
	for(int izeta = 0; izeta < gridDim(dimOption::K); izeta++) {
		float zeta_min = idxToPos(izeta, dimOption::K) - 0.5 * m_delta[dimOption::K];
		float zeta_max = idxToPos(izeta, dimOption::K) + 0.5 * m_delta[dimOption::K];

		auto result = integrate(100, zeta_min, zeta_max, [&](float zeta)-> float {
			float wavelength = pow(2, zeta);
			float waveNumber = 2 * PI / wavelength;
			/* Integrate the dW/dK to get the group speed.
			//cg = dW/dK
			// W = sqrt(Gravity * K)
			// dW = 0.5 * sqrt(Gravity/K) * dK
			*/
			float cg = 0.5 * sqrt(GRAV / waveNumber);
			return cg;
		});

		m_groupSpeeds[izeta] = result;
	}

}

float SurfaceWavelet::idxToPos(const int idx, const int dim) const
{
	return m_xmin[dim] + (idx + 0.5) * m_delta[dim];
}


Vec4 SurfaceWavelet::idxToPos(std::array<int, 4> idx) const
{
	return Vec4{ idxToPos(idx[dimOption::X],dimOption::X), idxToPos(idx[dimOption::Y],dimOption::Y),
				 idxToPos(idx[dimOption::THETA], dimOption::THETA), idxToPos(idx[dimOption::K],dimOption::K) };
}

float SurfaceWavelet::posToGrid(float pos, int dim) const
{
	return (pos - m_xmin[dim])*m_idx[dim] - 0.5;
}


Vec4 SurfaceWavelet::posToGrid(Vec4 pos4) const
{
	return Vec4{ posToGrid(pos4[dimOption::X], dimOption::X), posToGrid(pos4[dimOption::Y],dimOption::Y), 
		         posToGrid(pos4[dimOption::THETA], dimOption::THETA), posToGrid(pos4[dimOption::K], dimOption::K) };
}

int SurfaceWavelet::posToIdx(float pos, int dim) const
{
	return round(posToGrid(pos, dim));
}

SurfaceWavelet::Idx SurfaceWavelet::posToIdx(Vec4 pos4) const
{
	return Idx{ posToIdx(pos4[dimOption::X], dimOption::X), posToIdx(pos4[dimOption::Y],dimOption::Y),
			    posToIdx(pos4[dimOption::THETA], dimOption::THETA), posToIdx(pos4[dimOption::K],dimOption::K) };
}



auto SurfaceWavelet::correctFringeValues() const
{
	return [this](int ix, int iy, int itheta, int izeta) {
		itheta = pos_modulo(itheta, gridDim(dimOption::THETA));

		if (izeta < 0 || izeta >= gridDim(dimOption::K)) { return 0.0f; }

		if (ix < 0 || ix >= gridDim(dimOption::X) || iy < 0 || iy >= gridDim(dimOption::Y))
		{
			return defaultAmplitude(itheta, izeta);
		}

		return m_Amplitude(ix, iy, itheta, izeta);
	};
}

auto SurfaceWavelet::interpolateAmplitude() const
{
	auto extended_grid = correctFringeValues();
	
	auto domain = [this](int ix, int iy, int itheta, int izeta) -> bool {
		return true;
	};

	auto interpolation = InterpolationDimWise(
		LinearInterpolation, LinearInterpolation, LinearInterpolation,
		ConstantInterpolation);
	
	auto interpolated_grid =
		DomainInterpolation(interpolation, domain)(extended_grid);

	return [interpolated_grid, this](Vec4 pos4) mutable {
		Vec4 ipos4 = posToGrid(pos4);
		return interpolated_grid(ipos4[dimOption::X], ipos4[dimOption::Y], ipos4[dimOption::THETA], ipos4[dimOption::K]);
	};
}

float SurfaceWavelet::getAmplitude(Vec4 pos4) const
{
	return interpolateAmplitude()(pos4);
}

int SurfaceWavelet::gridDim(const int dim) const
{
	return m_Amplitude.dimSize(dim);
}

float SurfaceWavelet::defaultAmplitude(const int itheta, const int izeta) const
{
	if (itheta == 5 * gridDim(dimOption::THETA) / 16)
		return 0.1;
	return 0.0;
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
	AdvectionStep(m_dt);
	WaveDiffusionStep(m_dt);

	precomputeProfileBuffer();
	m_time += m_dt;
}

void SurfaceWavelet::AdvectionStep(const float dt)
{
	auto extended_grid = correctFringeValues();

	for (int i = 0; i < gridDim(dimOption::X); i++)
	{
		#pragma omp parallel for
		for (int j = 0; j < gridDim(dimOption::Y); j++)
		{
			Vec2 pos = Vec2{ idxToPos(i,dimOption::X), idxToPos(j, dimOption::Y) };

			for (int theta = 0; theta < gridDim(dimOption::THETA); theta++)
			{
				for (int zeta = 0; zeta < gridDim(dimOption::K); zeta++)
				{
					Vec4 pos4 = idxToPos({ i,j,theta,zeta });
					//cout << (ParallelStream() << "position: " << pos4[0] << ", " << pos4[1]).toString()<< endl;
					
					float waveAngle = pos4[dimOption::THETA];
					Vec2 vel = m_groupSpeeds[zeta] * Vec2{ cosf(waveAngle),sinf(waveAngle) };
					
					//SL Advection
					Vec4 prev_pos = pos4;
					prev_pos[dimOption::X] -= dt * vel[dimOption::X];
					prev_pos[dimOption::Y] -= dt * vel[dimOption::Y];

					prev_pos = boundaryConditions(prev_pos);
					new_Amplitude(i, j, theta, zeta) = getAmplitude(prev_pos);
					/*if (fabs(new_Amplitude(i, j, theta, zeta)) > 1e-10 )
						cout << (ParallelStream() << "Amplitude: " << new_Amplitude(i, j, theta, zeta)).toString() << endl;
					*/
				}
			}

		}
	}
	std::swap(new_Amplitude, m_Amplitude);
}

void SurfaceWavelet::WaveDiffusionStep(const float dt)
{
	auto grid = correctFringeValues();
	
	for (int i = 0; i < gridDim(dimOption::X); i++)
	{
	#pragma omp parallel for
		for (int j = 0; j < gridDim(dimOption::Y); j++)
		{
			for (int theta = 0; theta < gridDim(dimOption::THETA); theta++)
			{
				for (int zeta = 0; zeta < gridDim(dimOption::K); zeta++)
				{
					Vec4 pos4 = idxToPos({ i,j,theta,zeta });
					float gamma = 2 * 0.025 * m_groupSpeeds[zeta] * dt * m_idx[dimOption::X];

					//within 2 grid nodes away from boundary
					if (pos4[0] >= 2 && (m_xmax[0] - pos4[0]) >= 2 && pos4[1] >= 2 && (m_xmax[1] - pos4[1]) >= 2)
					{
						new_Amplitude(i, j, theta, zeta) = (1 - gamma) * grid(i, j, theta, zeta) + 
							gamma * 0.5 * 
								(grid(i, j, theta + 1, zeta)+
								 grid(i,j, theta - 1, zeta));
					}
					else
					{
						m_Amplitude(i, j, theta, zeta) = grid(i, j, theta, zeta);
					}
				
				}
			}
		}
	}
	std::swap(new_Amplitude, m_Amplitude);
}

void SurfaceWavelet::propogate(float *&source_height, float dt)
{
	this->m_dt = dt;
	addAmplitude(source_height);
	timeStep();


	genHeightMap (m_height);

	res_height = m_height;


}
