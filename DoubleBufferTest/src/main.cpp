#include <iostream>
#include <utility>
#include <string>
#include "DubBuffer.h"

void main(int argc, char* argv)
{
	std::cout << "Hello World" << std::endl;

	float * a, *b;
	int size = 3;
	a = (float*)malloc(size * sizeof(float));
	b = (float*)malloc(size * sizeof(float));

	for (int i = 0; i < size; i++)
	{
		a[i] = i;
		b[i] = size + i;
	}

	

	DubBuffer test(a , b);
	
	float* tA = test.first();
	float* tB = test.second();
	std::cout << "first: ";
	for (int i = 0; i < size; i++)
	{ 
		std::cout << a[i] << ",";
	}
	std::cout << std::endl;
	std::cout << "second: ";
	for (int i = 0; i < size; i++)
	{
		std::cout << b[i] << ",";
	}
	std::cout << std::endl;
	
	std::swap(a, b);
	test.switchBuffer();
	tA = test.first();
	tB = test.second();

	std::cout << "first: ";
	for (int i = 0; i < size; i++)
	{
		std::cout << a[i] << ",";
	}
	std::cout << std::endl;
	std::cout << "second: ";
	for (int i = 0; i < size; i++)
	{
		std::cout << b[i] << ",";
	}
	std::cout << std::endl;
	
	
	std::string temp;
	std::cin >> temp;
}