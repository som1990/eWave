//******************************************************************
// FileName: eWaveSim.cpp
// Author: Soumitra Goswami
// Date: 5/08/2018
// Last Update: 5/23/2018
// 
// Description: Implementation of all the methods defined in "eWaveSim.h"
//****************************************************************** 
#include "eWaveSim.h"
#define REAL 0
#define IMAG 1
#define PI 3.14159265359
#include <iostream>

using namespace std;

//************************************
// Method:    initmaps2
// FullName:  eWaveSim::initmaps2
// Access:    private 
// Returns:   void
// Parameter: float(* & map2) - array to populate
// Parameter: int size - size of the array
// Parameter: float value - value used to populate the array.
// 
// Description: Quick function to populate the array with a given value.
//************************************

void eWaveSim::initmaps(float (*&map2), int size, float value)
{

#pragma omp parallel for
	for (int i = 0; i < size; i++)
	{
		(map2)[i] = value;
	}
}

//Constructor function
eWaveSim::eWaveSim()
{
}


eWaveSim::eWaveSim(int gridX, int gridY)
	: simGridX(gridX/m_dx),simGridY(gridY/m_dy)
{
}

//CleanUp of data
eWaveSim::~eWaveSim()
{
	delete [] m_height;
	delete [] res_height;
	delete [] prev_height;
	delete [] vel_potential;
	delete [] prev_velPotential;
	delete [] cPrev_height;
	delete [] cPrev_vel;
	delete [] cH;
	delete [] cVel;
	delete [] cHfft;
	delete [] cVelfft;
	delete [] driftVel;
	delete [] m_ambientWaves;
	delete [] m_ambWaveSource;
	
}


//************************************
// Method:    initFields
// FullName:  eWaveSim::initFields
// Access:    public 
// Returns:   void

// Parameter: int nGridX - number of simulation grids in X
// Parameter: int nGridY - number of simulation grids in Y
// 
// Description: This function initializes all the fields to be ready to use based on the 
//              simulation gridsize provided
//************************************

void eWaveSim::initFields(int nGridX, int nGridY)
{
	this->simGridX = nGridX/m_dx;
	this->simGridY = nGridY/m_dy;
	
	int size = simGridX*simGridY;

	m_height = new float[size];
	prev_height = (float*) malloc(sizeof(float)*size);
	vel_potential = new float[size];
	prev_velPotential = (float*) malloc(sizeof(float)*size);
	res_height = new float[size];
	driftVel = new float[size*2];
	m_ambientWaves = new float[size];
	m_ambWaveSource = new float[size];


	cH = (fftwf_complex*)malloc(sizeof(fftwf_complex)*size);
	cVel = (fftwf_complex*)malloc(sizeof(fftwf_complex)*size);
	cHfft = (fftwf_complex*)malloc(sizeof(fftwf_complex)*size);
	cVelfft = (fftwf_complex*)malloc(sizeof(fftwf_complex)*size);
	cPrev_height = (fftwf_complex*)malloc(sizeof(fftwf_complex)*size);
	cPrev_vel = (fftwf_complex*)malloc(sizeof(fftwf_complex)*size);

	initmaps(m_ambWaveSource, size, 0.0);
	initmaps(m_ambientWaves, size, 0.0);
	initmaps(m_height, size, 0.0);
	//initmaps2(res_height, size, 0.0);
	initmaps(vel_potential, size, 0.0);


}

// Return function returning the heightfield at i and j
float eWaveSim::height(int &i, int &j)
{
	int index = i + simGridX*j;
	float res = res_height[index];
	return res;
}

// Return function returning the heightfield at index position
float eWaveSim::height(int &index)
{
	float res = res_height[index];
	return res;
}


//************************************
// Method:    addingSources
// FullName:  eWaveSim::addingSources
// Access:    private 
// Returns:   void

// Parameter: float *& source_height - source field for height
// 
// Description: Source field gets added to the height field in this function.
//				We use a subtle and gradual addition by using delta_time.
//				height = height + source*time*amplifier.
//************************************

void eWaveSim::addingSources(float *&source_height)
{
	int size = simGridX*simGridY;
#pragma omp parallel for
	for (int i = 0; i < size; i++)
	{
		//I was getting a subtle effect and hence multiplied it with a constant
		float temp = source_height[i]*m_dt*2;
		float temp2 = m_height[i];
		m_height[i] = temp2 + temp;
	}
}

//************************************
// Method:    applyObstruction
// FullName:  eWaveSim::applyObstruction
// Access:    private 
// Returns:   void

// Parameter: float * & sourceObstruction - obstruction map to use
// 
// Description: We multiply the obstruction map with height and vel_potential
//              so that they don't get considered in the simulation.
//              ->height = height * obstructionMap
//              ->velPotential = velPotential * obstructionMap
//************************************

void eWaveSim::applyObstruction(float *&sourceObstruction)
{
	for (int j = 0; j < simGridY; j++)
	{
#pragma omp parallel for
		for (int i = 0; i < simGridX; i++)
		{
			int index = i + j*simGridX;
			m_height[index] *= sourceObstruction[index];
			vel_potential[index] *= sourceObstruction[index];
		}
	}
}
//FOURIER TRANSFORM FUNCTIONS

/*
void eWaveSim::fft(float*&in, float *&out)
{
	fftwf_plan fft_plan;
	int size = simGridX*simGridY;
	out = (float *) malloc(sizeof(float)*size);
	fft_plan = fftwf_plan_r2r_1d(size, in, out, FFTW_R2HC, FFTW_ESTIMATE);
	fftwf_execute(fft_plan);
	fftwf_destroy_plan(fft_plan);
	fftwf_cleanup();
}
*/

//Forward Fourier Transform Function
void eWaveSim::fft(fftwf_complex *&in, fftwf_complex *&out)
{
	fftwf_plan fft_plan;
	int size = simGridX*simGridY;
	fft_plan = fftwf_plan_dft_1d(size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftwf_execute(fft_plan);
	fftwf_destroy_plan(fft_plan);
	fftwf_cleanup();

}
/*
void eWaveSim::ifft(float *&in, float *&out)
{
	fftwf_plan ifft_plan;
	int size = simGridX*simGridY;
	out = (float*)malloc(size*sizeof(float));
	ifft_plan = fftwf_plan_r2r_1d(size, in, out,FFTW_HC2R, FFTW_ESTIMATE);
	fftwf_execute(ifft_plan);
	fftwf_destroy_plan(ifft_plan);
	fftwf_cleanup();

//	for (int i = 0; i < size; i++)
//	{
//		out[i] = out[i] / size;
//	}
	
}
*/

// Inverse Fourier transform functions
void eWaveSim::ifft(fftwf_complex *&in, fftwf_complex *&out)
{
	fftwf_plan ifft_plan;
	int size = simGridX*simGridY;
	ifft_plan = fftwf_plan_dft_1d(size, in, out,FFTW_BACKWARD,FFTW_ESTIMATE);
	fftwf_execute(ifft_plan);
	fftwf_destroy_plan(ifft_plan);
	fftwf_cleanup();

	//Normalizing the inverse fft.
	for (int j = 0; j < simGridY; j++)
	{
#pragma omp parallel for
		for (int i = 0; i < simGridX; i++)
		{
			int index = i + j*simGridX;
			out[index][REAL] = out[index][REAL] / size;
			out[index][IMAG] = out[index][IMAG] / size;
		}
	}

}






//************************************
// Method:    calc_eWave
// FullName:  eWaveSim::calc_eWave
// Access:    private 
// Returns:   void
// Qualifier:
// Parameter: fftwf_complex * & h - newHeightField
// Parameter: fftwf_complex * & v - newVelocityPotentialField
// Parameter: fftwf_complex * & oldh - oldHeightField
// Parameter: fftwf_complex * & oldv - oldVelocityField
// 
// Description: We run the eWave equations on the FFTransformed height(h') and velocityPotential(v')
//              ->h'(k,t+dT) = cos(w(k)*dT)*h'(k,t) + (k/w(k))*sin(w(k)*dt)*v'(k,t)  eq.19
//              ->v'(k,t+dT) = cos(w(k)*dT)*v'(k,t) + (g/w(k))*sin(w(k)*dt)*h'(k,t)  eq.20
//				
//              here, 
//              k is the absolute magnitude of the Fourier Vector
//              w(k) = sqrt(g*k)
//              g = gravity 
//************************************
void eWaveSim::calc_eWave(fftwf_complex *&h, fftwf_complex *&v, fftwf_complex *&oldh, fftwf_complex *&oldv)
{
	int size = simGridX*simGridY;


	for (int j = 0; j < simGridY; j++)
	{
		for (int i = 0; i < simGridX; i++)
		{

			//Direction 1 from Siggraph Course 2004
			// 		float n = i - int((simGridX) / 2.0);
			// 		float m = i - int((simGridY) / 2.0);
			float dkx = 2 * PI / m_Lx;
			float dky = 2 * PI / m_Ly;

			float kx, ky;
			if (i <= simGridX / 2.0)
				kx = i*dkx;
			else
				kx = (simGridX - i )*dkx;

			if (j <= simGridY / 2.0)
				ky = j*dky;
			else
				ky = ( simGridY - j)*dky;
			//  	float m = simGridY;
			// 		float m_kx = 2 * PI*n / simGridX;
			// 		float m_ky = 2 * PI*m / simGridY;
			//  	m_k = sqrtf(m_kx*m_kx + m_ky*m_ky);

			//cout << "w(k): " << w << ", k: " << m_k << endl;
			//cout << "h(r): " << h[i][REAL] << ", h(imag): " << h[i][IMAG] << endl;
			//cout << "vel(r): " << v[i][REAL] << ", vel(imag): " << v[i][IMAG] << endl;

			int index = i + j*simGridX;

			//Direction 2 with (absolute magnitude of Fourier vector)

			float m_k = sqrtf(kx*kx + ky*ky);

			
			//calculating w(k) = sqrt(g*k)

			float w = sqrtf(gravity*m_k);

			float coswk = cos(w*m_dt);
			float sinwk = sin(w*m_dt);
			//Equations 19 and 20 of the eWave Note/
			// h - FFTransformed height. Init = 0
			// v - FFTranformed Vel potential. Init = 0
			if (m_k != 0)
			{
				h[index][REAL] = coswk*oldh[index][REAL] + (m_k / w)*sinwk*oldv[index][REAL]; // eq 19
				h[index][IMAG] = coswk*oldh[index][IMAG] + (m_k / w)*sinwk*oldv[index][IMAG]; // eq 19

				v[index][REAL] = coswk*oldv[index][REAL] - (gravity / w)*sinwk*oldh[index][REAL]; // eq 20	
				v[index][IMAG] = coswk*oldv[index][IMAG] - (gravity / w)*sinwk*oldh[index][IMAG]; // eq 20
			}
		}
	}
}

//Real to complex conversion

void eWaveSim::convert_r2c(float* &r, fftwf_complex *&c)
{
	int size = simGridX*simGridY;
	
	for (int j = 0; j < simGridY; j++)
	{
		for (int i = 0; i < simGridX; i++)
		{
			int index = i + j*simGridX;
			c[index][REAL] = r[index];
			c[index][IMAG] = 0.0f;

		}
	}
}

//Complex to real conversion
void eWaveSim::convert_c2r(fftwf_complex *&c, float *&r)
{
	int size = simGridX*simGridY;
	for (int j = 0; j < simGridY; j++)
	{
		for (int i = 0; i < simGridX; i++)
		{
			int index = i + j*simGridX;
			r[index] = c[index][REAL];
		}
	}
}
//---------------------------------------------------------------------------------
//DEBUG PRINT FUNCTIONS
void printlist(float *&m_list, int size, char* c)
{
	for (int i = 0; i < size; i++)
	{
		cout << c << ": " << m_list[i] << endl;
	}
}

void printlist(fftwf_complex *&m_list, int size, char* c)
{
	float f_max = -INT_MAX;
	for (int i = 0; i < size; i++)
	{
		if (f_max < m_list[i][IMAG])
			f_max = m_list[i][IMAG];
		//cout << c << "[r]: " << m_list[i][REAL] << endl;
		//if(fabs(m_list[i][IMAG]) >= 1e-10)
			//cout << c << "[i]: " << m_list[i][IMAG] << endl;
	}
	cout << c << "[IMAG][max]: " << f_max << endl;
}

//---------------------------------------------------------------------------------
//UPDATE FUNCTION

//************************************
// Method:    propogate
// FullName:  eWaveSim::propogate
// Access:    public 
// Returns:   void
// Qualifier:
// Parameter: float * & source_height - painted heightMap
// Parameter: float * & sourceObstruction - painted obstructionMap
// Parameter: float dt - Simulation time.
// 
// Description: This is the update function and goes through all the steps necessary for the program
//             - Initializes and updates all the necessary variables
//             - Adds the sourceHeightMap to the heightField
//             - Converts the real arrays from floats to fftwf_complex data types
//             - Creates a duplicate of these arrays to be utilized in the eWaveEquations
//             - Runs the eWave Equations
//             - Runs an inverse FFT to the updated heightfield 
//             - Converts the complex arrays used by the FFTs back to real arrays.
//             - Cleanup all the variables
//             - Updates the resultantHeight to be used by the display Functions  
//************************************

void eWaveSim::propogate(float *&source_height, float *&sourceObstruction, float dt)
{
	this->m_dt = dt;
	gravity = 9.8f;

	int size = simGridX*simGridY;

	addingSources(source_height);
	applyObstruction(sourceObstruction);
	convert_r2c(m_height, cH);
	convert_r2c(vel_potential, cVel);

	//fft(m_height, h_fft);
	//fft(vel_potential, vel_fft);
	fft(cH, cHfft);
	fft(cVel, cVelfft);
	
	//memcpy(previous_height, h_fft, size*sizeof(float));
	//memcpy(prev_vel, vel_fft, size*sizeof(float));
	
	memcpy(cPrev_height, cHfft, size*sizeof(fftwf_complex));
	memcpy(cPrev_vel, cVelfft, size * sizeof(fftwf_complex));


	calc_eWave(cHfft, cVelfft, cPrev_height, cPrev_vel);
	//printlist(cHfft, size, "chfft");
	//ifft(h_fft, m_height);
	//ifft(vel_fft, vel_potential);
	//ifft(cHfft, m_height);
	//ifft(cVelfft, vel_potential);
	ifft(cHfft, cH);
	ifft(cVelfft, cVel);

 //	printlist(cH, size, "ch");
	convert_c2r(cH, m_height);
	convert_c2r(cVel, vel_potential);
	//printlist(cH, size,"cheight");
	//printlist(m_height,size, "height");
	

	res_height = m_height;
}


