//******************************************************************
// FileName: eWaveSim.cpp
// Author: Soumitra Goswami
// Date: 5/08/2018
// Last Update: 7/23/2018
// 
// Description: Implementation of all the methods defined in "eWaveSim.h"
//****************************************************************** 
#include "eWaveSim.h"
#define REAL 0
#define IMAG 1
#define PI 3.14159265359
#include <iostream>
#include <utility>
using namespace std;


void bufferSwap(fftwf_complex* &a, fftwf_complex* &b)
{
	fftwf_complex* temp = a;
	a = b;
	b = temp;
}

void bufferSwap(float* &a, float* &b)
{
	float* temp;
	temp = a;
	a = b;
	b = temp;
}

static float fclamp(float x,float minVal, float maxVal)
{
	x = (x<minVal) ? minVal : ((x>maxVal) ? maxVal : x);
	return x;
}

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

void eWaveSim::initmaps(float (*&map2), int size, float value,int dim)
{

	if (dim == 1)
	{
		#pragma omp parallel for
		for (int i = 0; i < size; i++)
		{
			(map2)[i] = value;
		}
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
	delete [] obs_height;
	delete [] vel_potential;
	delete [] prev_velPotential;
	delete [] obs_velPot;
	delete [] driftVel;
	delete [] m_ambientWaves;
	delete [] m_ambWaveSource;


	free(cH);
	free(cVel);
	free(cHfft);
	free(cVelfft);
    free(cPrev_height);
	free(cPrev_vel);
	
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
	res_height = new float[size];
	prev_height = new float[size];
	obs_height = new float[size];
	vel_potential = new float[size];
	prev_velPotential = new float[size];
	obs_velPot = new float[size];
	driftVel = new float[size*2];
	m_ambientWaves = new float[size];
	m_ambWaveSource = new float[size];


	cH = (fftwf_complex*)malloc(sizeof(fftwf_complex)*size);
	cVel = (fftwf_complex*)malloc(sizeof(fftwf_complex)*size);
	cHfft = (fftwf_complex*)malloc(sizeof(fftwf_complex)*size);
	cVelfft = (fftwf_complex*)malloc(sizeof(fftwf_complex)*size);
	cPrev_height = (fftwf_complex*)calloc(size,sizeof(fftwf_complex));
	cPrev_vel = (fftwf_complex*)calloc(size,sizeof(fftwf_complex));


	
	initmaps(m_ambWaveSource, size, 0.0,1);
	initmaps(m_ambientWaves, size, 0.0,1);
	initmaps(m_height, size, 0.0,1);
	initmaps(obs_height, size, 0.0,1);
	initmaps(obs_velPot, size, 0.0,1);
	initmaps(vel_potential, size, 0.0,1);
	initDriftVel(driftVel);
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
// Method:    initDriftVel
// FullName:  eWaveSim::initDriftVel
// Access:    private 
// Returns:   void
// Qualifier:
// Parameter: float * & dVel - drift velocity to generate
// Description:	Generates a flow map moving horizontally. Can be 
//				Substituted for an actual hand generated flow map.
//************************************
void eWaveSim::initDriftVel(float *&dVel)
{
	int size = simGridX*simGridY*2;
#pragma omp parallel for
	for (int i = 0; i < size; i = i + 2)
	{
		dVel[i + 0] = -1.0f * driftVelScale;
		dVel[i + 1] = 0;
	}
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
		float temp2 = m_height[i] + obs_height[i];
		m_height[i] = temp2 + temp + m_ambientWaves[i] + m_ambWaveSource[i];
		vel_potential[i] += obs_velPot[i];
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
//					->height = height * obstructionMap
//					->velPotential = velPotential * obstructionMap
//				We also generate the ambient waves source function to cancel out ambient 
//				waves inside the obstruction volume
//************************************

void eWaveSim::applyObstruction(float *&sourceObstruction)
{
	for (int j = 0; j < simGridY; j++)
	{
#pragma omp parallel for
		for (int i = 0; i < simGridX; i++)
		{
			int index = i + j*simGridX;
//			m_height[index] *= sourceObstruction[index];
//			vel_potential[index] *= sourceObstruction[index];
			//Generating ambient Waves
			m_ambWaveSource[index] = fclamp(1.0f - sourceObstruction[index],0.0f,1.0f) * -m_ambientWaves[index];
			obs_height[index] = fclamp(1.0f - sourceObstruction[index], 0.0f, 1.0f) * -m_height[index];
			obs_velPot[index] = fclamp(1.0f - sourceObstruction[index], 0.0f, 1.0f) * -vel_potential[index];
		}
	}
}

//************************************
// Method:    boundaryConditions
// FullName:  eWaveSim::boundaryConditions
// Access:    private 
// Returns:   float -	returns 1.0 if inside the boundary, or a value of
//						[0,1] within the boundary determined by a function.
// Qualifier:
// Parameter: float x - horizontal position on the Sim Grid 
// Parameter: float y - vertical position on the Sim Grid
// 
// Description:	Dampens the periodicity of the waves by smoothing the height to zero 
//				at the edges of the grid based on user controlled function.
//				function of padding f(d) = pow(d/d0,alpha) where f(d) falls in [0,1] range
//					d - the distance from edge	
//					d0 - the padding distance
//					alpha - exponential decay of the function
//				What works well in this case is a padding of 10% of the grid dimensions and
//				alpha to be 0.05
//************************************
float eWaveSim::boundaryConditions(float x, float y)
{
	float x0 = x; float x1 = simGridX - x;
	float y0 = y; float y1 = simGridY - y;

	float leftBound, rightBound, topBound, bottomBound;
	leftBound = rightBound = topBound = bottomBound = 1.0;

	//10% of GridSize
	float hori_Padding = .1f * simGridX;
	float vert_Padding = .1f * simGridY;

	if (x0 < 0.0f)
		leftBound = 0;
	else
		leftBound = fclamp(pow(x0 / hori_Padding, bound_alpha), 0.0, 1.0);
	if (x1 < 0.0f)
		rightBound = 0;
	else
		rightBound = fclamp(pow(x1 / hori_Padding, bound_alpha), 0.0, 1.0);
	if (y0 < 0.0f)
		bottomBound = 0;
	else
		bottomBound = fclamp(pow(y0 / vert_Padding, bound_alpha), 0.0, 1.0);
	if (y1 < 0.0f)
		topBound = 0;
	else
		topBound = fclamp(pow(y1 / vert_Padding, bound_alpha), 0.0, 1.0);

	return topBound*bottomBound*rightBound*leftBound;
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
				kx = (simGridX - i)*dkx;

			if (j <= simGridY / 2.0)
				ky = j*dky;
			else
				ky = (simGridY - j)*dky;
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



//****************************************************************************
// Method:    SLadvection
// FullName:  eWaveSim::SLadvection
// Access:    private 
// Returns:   void
// Qualifier:
// Parameter: float * & vel - velocity to use to advect the field
// Parameter: float * & newField - final field after advection 
// Parameter: float * & oldField - field before advection
// Parameter: int dim - dimensions of the variable.
// 
// Description: Semi Lagrangian advection is a stable method to propogate elements on
//				the grid using the stored grid velocity.
//				The method involves:
//					- back tracing each grid point to it's previous state based on the velocity at the grid.
//					  in our case its the grid points (i,j) we are tracing back
//					  newPos(x,y) = pos(i,j)*gridCellSize - velocity(i,j)*delta_T
//					- Once we have a new position. We bilinear interpolate based on the surrounding 4 grid points.
//					- We then use the resulting value and store it in our new field.
//					
//				This proves to be pretty stable and reliable. One shortcoming of the method is that we lose information
//				because of interpolating each time. 
//***************************************************************************
void eWaveSim::SLadvection(float *&vel, float *&newField, float* &oldField, int dim)
{
	for (int j = 0; j < simGridY; j++)
	{
#pragma omp parallel for
		for (int i = 0; i < simGridX; i++)
		{
			int index = i + ((simGridX)*j);

			//Semi-Lagrangian trace back
			float x = i*m_dx - (vel[index * 2 + 0] * m_dt*1.0/m_Lx);
			float y = j*m_dy - (vel[index * 2 + 1] * m_dt*1.0/m_Ly);

			//Checking limits 
			if (x < 0.5) { x = 0.5; } if (x >(simGridX - 0.5)) { x = simGridX - 0.5; }
			int i_prev = (int)x; int i_next = i_prev + 1;
			i_next = (i_next >= simGridX) ? (simGridX - 1) : i_next;

			if (y < 0.5) { y = 0.5; } if (y >(simGridY - 0.5)) { y = simGridY - 0.5; }
			int j_prev = (int)y; int j_next = j_prev + 1;
			j_next = (j_next >= simGridY) ? (simGridY - 1) : j_next;
			

			//Bilinear interpolation
			//Finding position inside the grid rescaled from 0-1
			float s1 = x - i_prev; float s0 = 1.0 - s1;
			float t1 = y - j_prev; float t0 = 1.0 - t1;

			int ij = index;
			int i0j0 = i_prev + (j_prev*(simGridX));
			int i0j1 = i_prev + (j_next*(simGridX));
			int i1j0 = i_next + (j_prev*(simGridX));
			int i1j1 = i_next + (j_next*(simGridX));

			//Generating dampening in the boundaries of the grid
			float bound = 1.0f;
			bound = boundaryConditions(x, y);

			if (dim == 1)
			{
				newField[ij] = s0*(t0*oldField[i0j0] + t1*oldField[i0j1]) +
					s1*(t0*oldField[i1j0] + t1*oldField[i1j1]);
				newField[ij] *= bound;
			}

			if (dim == 2)
			{
				newField[(ij * 2) + 0] = s0*(t0*oldField[i0j0 * 2 + 0] + t1*oldField[i0j1 * 2 + 0]) +
					s1*(t0*oldField[i1j0 * 2 + 0] + t1*oldField[i1j1 * 2 + 0]);
				newField[(ij * 2) + 0] *= bound;
				
				newField[(ij * 2) + 1] = s0*(t0*oldField[i0j0 * 2 + 1] + t1*oldField[i0j1 * 2 + 1]) +
					s1*(t0*oldField[i1j0 * 2 + 1] + t1*oldField[i1j1 * 2 + 1]);
				newField[(ij * 2) + 1] *= bound;
			}

			if (dim == 3)
			{
				newField[(ij * 3) + 0] = s0*(t0*oldField[i0j0 * 3 + 0] + t1*oldField[i0j1 * 3 + 0]) +
					s1*(t0*oldField[i1j0 * 3 + 0] + t1*oldField[i1j1 * 3 + 0]);
				newField[(ij * 3) + 0] *= bound;

				newField[(ij * 3) + 1] = s0*(t0*oldField[i0j0 * 3 + 1] + t1*oldField[i0j1 * 3 + 1]) +
					s1*(t0*oldField[i1j0 * 3 + 1] + t1*oldField[i1j1 * 3 + 1]);
				newField[(ij * 3) + 1] *= bound;

				newField[(ij * 3) + 2] = s0*(t0*oldField[i0j0 * 3 + 2] + t1*oldField[i0j1 * 3 + 2]) +
					s1*(t0*oldField[i1j0 * 3 + 2] + t1*oldField[i1j1 * 3 + 2]);
				newField[(ij * 3) + 2] *= bound;
			}
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
	//test
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
	applyObstruction(sourceObstruction);
	addingSources(source_height);
	
	convert_r2c(m_height, cH);
	convert_r2c(vel_potential, cVel);

	//fft(m_height, h_fft);
	//fft(vel_potential, vel_fft);
	fft(cH, cHfft);
	fft(cVel, cVelfft);
	//printlist(cHfft, size, "cHfft");

	bufferSwap(cPrev_height,cHfft);
	bufferSwap(cPrev_vel,cVelfft);
	//memcpy(cPrev_height, cHfft, size*sizeof(fftwf_complex));
	//memcpy(cPrev_vel, cVelfft, size * sizeof(fftwf_complex));

	//printlist(cPrev_height, size, "prevHeight");
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
	
	bufferSwap(prev_height, m_height);
	bufferSwap(prev_velPotential, vel_potential);
	//memcpy(prev_height, m_height, sizeof(float)*size);
	//memcpy(prev_velPotential, vel_potential, sizeof(float)*size);

	SLadvection(driftVel, m_height, prev_height, 1);
	SLadvection(driftVel, vel_potential, prev_velPotential, 1);

	res_height = m_height;
}


