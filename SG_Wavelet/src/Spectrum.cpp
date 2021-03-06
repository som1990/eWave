#include "Spectrum.h"

Spectrum::Spectrum(double windSpeed) : m_windSpeed(windSpeed) {
	SpectrumType = 0;
}

Spectrum::Spectrum(double windSpeed, float windDirection) 
: m_windSpeed(windSpeed), m_windDirection(windDirection){
	SpectrumType = 1;
}

double Spectrum::maxZeta() const
{
	return log2(10);
}

double Spectrum::minZeta() const
{
	return log2(0.03);
}

//Pierson-Moskowitz spectrum
double Spectrum::operator()(double zeta) const
{
#define PI 3.14159265359
#define GRAV 9.81
	float waveLength = pow(2, zeta);
	float k_length = 2.0 * PI / waveLength;
	double w = sqrt(GRAV*k_length);

	if (SpectrumType == 0)
	{
	#define ALPHA 0.0081
		double w_p = 0.879 * GRAV / m_windSpeed;
		
		
		double temp = w / w_p;
		double A = ALPHA * pow(GRAV,2) * pow(w,-5);
		// double A = pow(1.1, 1.5*zeta);//pow(2,1.5*zeta)original// According to the provided program but doesn't match the spectrum function
		//std::cout << A << std::endl;
		//double A = 0.00002613369 * pow(2, 2.5*zeta); // This is my calculation
		//double B = exp(-1.8038897788078411 * pow(4, zeta) / pow(m_windSpeed, 4));//This is right according to my calculation
		double B = exp(-1.25 * pow(temp,4));
		//std::cout << 0.139098 * sqrt(A * B) << std::endl;
		return 0.139098 * sqrt(A * B); //Random Gausian*random number (0,1) * Sqrt(Spectrum)

	}
	else if (SpectrumType == 1) // Phillips Spectrum without direction of wind consideration. Will get implemented during gathering of profileBuffer of the spectrum
	{
		float wind_angle = PI * m_windDirection / 180.0f;
		
		float damping = 1 / 1000; //Filtering out Really small waves 

		Vec2 wind_dir = Vec2{ cos(wind_angle),sin(wind_angle) };

		float L = m_windSpeed * m_windSpeed / GRAV;
		float phil = expf(-1 / (k_length*k_length*L*L)) / pow(k_length, 4);

		return 0.139098 * sqrt(phil*expf(-k_length * damping*damping));
	}
}



