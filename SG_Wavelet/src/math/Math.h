#pragma once


#include "ArrayAlgebra.h"
#include "Integration.h"

#include "interpolation/DomainInterpolation.h"
#include "interpolation/Interpolation.h"

constexpr int pos_modulo(int n, int d) { return (n % d + d) % d; }

template<class T>
constexpr const T& clamp(const T& v, const T& lo, const T& hi)
{
	return max(lo, min(v, hi));
}