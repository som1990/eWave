//*******************************************
// File: eWaveSim.h
// Author: Soumitra Goswami
// Description: Class file contains all the methods to simulate the interactive heightmap 
//              change based on the eWave note by Dr. Jerry Tessendorf.
// Date: 5/08/2018
// Last Edited: 6/22/2018 
// 
//******************************************* 


#ifndef eWaveSim_h__
#define eWaveSim_h__

#include <vector>
#include "fftw3.h"

class eWaveSim
{
public:
	//Constructors
	eWaveSim();
	eWaveSim(int gridX, int gridY);
	~eWaveSim();

	//Return functions for the main.cpp to access
	float height(int &i, int &j);
	float height(int &index);
	
	float ambientWaves(int i, int j);
	float ambientWaves(int index);

	//Initialize all the fields to have them ready to use.
	void initFields(int nGridX, int nGridY);

	//This is the update function and goes through all the steps necessary for the program
	void propogate(float *&source_height, float *&sourceObstruction, float dt);
	
	

private:
	//Quick function to assign a value to the entire array
	static void initmaps(float(*&map2), int size, float value);
	//Adds source_height user added to the heightfield
	void addingSources(float *&source_height);
	//Multiplies the height and velocityPotential with the obstructionfield
	void applyObstruction(float *&sourceObstruction);
	

	//Fast Fourier Transform functions using fftw libraries
	//void fft(float*&in, float *&out);
	void fft(fftwf_complex *&in, fftwf_complex *&out);
	//void ifft(float *&in, float *&out);
	void ifft(fftwf_complex *&in, fftwf_complex *&out);
	
	//Simulation of the heightmap based on the eWave note.

	void calc_eWave(fftwf_complex *&h, fftwf_complex *&v, fftwf_complex *&oldh, fftwf_complex *&oldv);
	
	//Conversions between real and complex variables
	void convert_r2c(float* &r, fftwf_complex *&c);
	void convert_c2r(fftwf_complex *&c, float *&r);
	
	//advection of the components based on velocity provided
	void SLadvection(float *&vel, float *&newField, float* &oldField, int dim);

	//Generating Fake Drift Velocity
	void initDriftVel(float *&dVel);

	//Introducing Dampening of waves on the boundaries
	float boundaryConditions(float x, float y);

	//Variables used
private:
	float driftVelScale = 100.0f;
	float gravity = 9.8f;
	
	//fft variables
	float m_Lx=1.0f, m_Ly=1.0f;
	
	//boundary variables
	
	float bound_alpha = 0.05;
	
	//grid variables
	float m_dx = 1.0f, m_dy = 1.0f;
	int simGridX, simGridY;
	//required variables
	float m_dt;
	float *m_height;
	float *prev_height; 
	float *vel_potential;
	float *prev_velPotential;
	float *driftVel;
	float *m_ambientWaves;
	float *m_ambWaveSource;

	float *obs_height;
	float *obs_velPot;
	fftwf_complex *cPrev_height;
	fftwf_complex *cPrev_vel;
	fftwf_complex *cHfft, *cVelfft;
	fftwf_complex *cH, *cVel;

	float *res_height;
};

#endif // eWaveSim_h__
