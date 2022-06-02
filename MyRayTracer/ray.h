#ifndef RAY_H
#define RAY_H

#include "vector.h"

class Ray
{
public:
	Ray(const Vector& o, const Vector& dir, double tm = 0.0) : origin(o), direction(dir), time(tm) {};

	Vector origin;
	Vector direction;
	double time;
};
#endif