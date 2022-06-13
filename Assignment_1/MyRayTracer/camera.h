#ifndef CAMERA_H
#define CAMERA_H

#include <cmath>
#include <stdio.h>
using namespace std;

#include "vector.h"
#include "ray.h"
#include "maths.h"

class Camera
{

private:
	Vector eye, at, up;
	float fovy, vnear, vfar, plane_dist, focal_ratio, aperture;
	float w, h;
	int res_x, res_y;
	Vector u, v, n;
	float time0, time1;

public:
	Vector GetEye() { return eye; }
	int GetResX() { return res_x; }
	int GetResY() { return res_y; }
	float GetFov() { return fovy; }
	float GetPlaneDist() { return plane_dist; }
	float GetFar() { return vfar; }
	float GetAperture() { return aperture; }

	Camera(Vector from, Vector At, Vector Up, float angle, float hither, float yon, int ResX, int ResY, float Aperture_ratio, float Focal_ratio, float t0 = 0.0, float t1 = 0.0) {
		time0 = t0;
		time1 = t1;
		eye = from;
		at = At;
		up = Up;
		fovy = angle;
		vnear = hither;
		vfar = yon;
		res_x = ResX;
		res_y = ResY;
		focal_ratio = Focal_ratio;

		// set the camera frame uvn
		n = (eye - at);
		plane_dist = n.length();
		n = n / plane_dist;

		u = up % n;
		u = u / u.length();

		v = n % u;

		//Dimensions of the vis window
		h = 2 * plane_dist * tan((PI * angle / 180) / 2.0f);
		w = ((float)res_x / res_y) * h;

		aperture = Aperture_ratio * (w / res_x); //Lens aperture = aperture_ratio * pixel_size

		printf("\nwidth=%f height=%f fov=%f, viewplane distance=%f, pixel size=%.3f\n", w, h, fovy, plane_dist, w / res_x);
		if (Aperture_ratio != 0) printf("\nDepth-Of-Field effect enabled with a lens aperture = %.1f\n", Aperture_ratio);
	}

	void SetEye(Vector from) {
		eye = from;
		// set the camera frame uvn
		n = (eye - at);
		plane_dist = n.length();
		n = n / plane_dist;
		u = up % n;
		u = u / u.length();
		v = n % u;
	}

	//motion blur
	void SetTime(float t0, float t1) {
		time0 = t0;
		time1 = t1;
	}

	Ray PrimaryRay(const Vector& pixel_sample, bool MOTION_BLUR) //  Rays cast from the Eye to a pixel sample which is in Viewport coordinates
	{
		// Powerpoint Whitted Ray-Tracing: Practice, slide 31
		Vector vector_x = u * w * (pixel_sample.x / res_x - 0.5f);
		Vector vector_y = v * h * (pixel_sample.y / res_y - 0.5f);
		Vector vector_z = n * -plane_dist;

		Vector ray_dir;

		ray_dir = (vector_x + vector_y + vector_z).normalize();

		float time = 0.0;
		if (MOTION_BLUR) {
			time = time0 + rand_float() * (time1 - time0);
		}

		return Ray(eye, ray_dir, time);
	}

	Ray PrimaryRay(const Vector& lens_sample, const Vector& pixel_sample, bool MOTION_BLUR) // DOF: Rays cast from  a thin lens sample to a pixel sample
	{
		// Powerpoint Distribution Ray-Tracing, Slide 39
		Vector ray_dir;
		Vector eye_offset;

		Vector ps;

		ps.x = w * (pixel_sample.x / res_x - 0.5f);
		ps.y = h * (pixel_sample.y / res_y - 0.5f);

		Vector p;

		p.x = ps.x * focal_ratio;
		p.y = ps.y * focal_ratio;

		//d = normalize((px - lsx)u + (py - lsy)v - fz
		ray_dir = (u * (p.x - lens_sample.x) + v * (p.y - lens_sample.y) + n * focal_ratio * -plane_dist);
		ray_dir = ray_dir.normalize();

		//eye_offset = eye + lsx * u + lsy * v
		eye_offset = eye + (u * lens_sample.x) + (v * lens_sample.y);


		float time = 0.0;
		if (MOTION_BLUR) {
			time = time0 + rand_float() * (time1 - time0);
		}

		return Ray(eye_offset, ray_dir, time);
	}
};

#endif