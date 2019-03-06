#pragma once
#ifndef Spectrum_h__
#define Spectrum_h__

#include <cmath>
#include "SurfaceWavelet_PCH.h"
class Spectrum {
public:
	Spectrum(double windSpeed);
	Spectrum(double windSpeed, float windDirection);

	//Max reasonable zeta value to consider for wavelength
	double maxZeta() const;
	
	//Min reasonable zeta value to consider for wavelength
	double minZeta() const;

	//Returns spectrum wave density based on the wavelength (= pow(2,zeta))
	double operator() (double zeta) const;//Pierson Moskowitz 

	double operator() (double zeta, float theta);

public:
	double m_windSpeed = 1;
	float m_windDirection = 0;
	int SpectrumType = 0; //0 - Pierson Moskowitz  1 - Phillips Spectrum without windDirectoinCalculation
	
};


#endif // Spectrum_h__