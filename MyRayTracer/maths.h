#ifndef __MATHS__
#define __MATHS__

#include <stdlib.h>
#include "vector.h"

#define PI				3.141592653589793238462f

// prototypes

unsigned int float_to_int(double x);
double min(double x0, double x1);
double max(double x0, double x1);
double clamp(const double x, const double min, const double max);
int	rand_int(void);
float rand_float(void);
double rand_double(void);
double rand_double(double min, double max);
Vector rnd_unit_disk(void);
Vector rnd_unit_sphere(void);
void set_rand_seed(const int seed);
uint8_t u8fromfloat(float x);
float u8tofloat(uint8_t x);

// inlined functions

//--------------------------------------------------------float to min
inline unsigned int
float_to_int(double x) {
	return ((x) >= 0 ? (unsigned int)((x)+0.5) : (unsigned int)((x)-0.5));
}

// ----------------------------------------------------------------- min

inline double
min(double x0, double x1) {
	return ((x0 < x1) ? x0 : x1);
}


// ----------------------------------------------------------------- max

inline double
max(double x0, double x1) {
	return ((x0 > x1) ? x0 : x1);
}

// ---------------------------------------------------- clamp

inline double
clamp(const double x, const double min, const double max) {
	return (x < min ? min : (x > max ? max : x));
}


// ---------------------------------------------------- rand_int
// a wrapper for rand()

inline int
rand_int(void) {
	return(rand());
}


// ---------------------------------------------------- rand_float

inline float
rand_float(void) {
	return((float)rand() / ((float)RAND_MAX+1.0));
}


// ---------------------------------------------------- rand_double

inline double
rand_double(void) {
	return((double)rand() / ((double)RAND_MAX + 1.0));
}

// ---------------------------------------------------- rand_double(min, max)

inline double
rand_double(double min, double max) {
	return min + (max - min)*rand_double();
}

// ---------------------------------------------------- rnd_unit_disk

inline Vector rnd_unit_disk(void) {
	Vector p;
	do {
		p = Vector(rand_float(), rand_float(), 0.0) * 2 - Vector(1.0, 1.0, 0.0);
	} while (p * p >= 1.0);
	return p;
}

// ---------------------------------------------------- rnd_unit_sphere
inline Vector rnd_unit_sphere(void) {
	Vector p;
	do {
		p = Vector(rand_float(), rand_float(), rand_float()) * 2 - Vector(1.0, 1.0, 0.0);
	} while (p * p >= 1.0);
	return p;
}

// ---------------------------------------------------- set_rand_seed
inline void
set_rand_seed(const int seed) {
	srand(seed);
}

// ---------------------------------------------------- float to byte (unsigned char)
inline uint8_t u8fromfloat(float x)
{
	//return (uint8_t)(x * 255.99f);
	return ((x * 255.99f) >= 255.0f ? 255 : (uint8_t)(x * 255.99f));
}

// ---------------------------------------------------- byte (unsigned char) to float
inline  float u8tofloat(uint8_t x)
{
	return (float)(x / 255.99f);
}

#endif
