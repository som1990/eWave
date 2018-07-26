#pragma once
class DubBuffer
{
public:
	DubBuffer(float* &a, float* &b);
	~DubBuffer();

	void switchBuffer();

	float* first() const ;
	float* second() const ;
private:
	float* m_first;
	float* m_second;
};

