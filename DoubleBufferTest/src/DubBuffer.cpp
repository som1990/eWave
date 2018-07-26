#include "DubBuffer.h"



DubBuffer::DubBuffer(float* &a, float* &b)
	:m_first(a),m_second(b)
{
}


DubBuffer::~DubBuffer()
{
}

void DubBuffer::switchBuffer()
{
	float* temp;
	temp = m_first;
	m_first = m_second;
	m_second = temp;

}

float* DubBuffer::first() const
{
	return m_first;
}

float* DubBuffer::second() const
{
	return m_second;
}

