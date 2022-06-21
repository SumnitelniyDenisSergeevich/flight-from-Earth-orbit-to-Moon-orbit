#pragma once 

#include "vector.h"
#include "matrix.h"
#include "ephemeris.h"


#include "log_duration.h"

struct Right {
	Vector r = { 0.,0.,0. };
	Vector v = { 0.,0.,0. };
	Vector f = { 0.,0.,0. };// ускорение
	double time = 0.;
};

void CheckConnection(const Vector& r, const Vector& planet_r, const double t);

void RungeCutte8(Vector& radius, Vector& speed, double t0, double t_end, double step, const dph::EphemerisRelease& de405);
Right RungeCutte4(const Vector& radius, const Vector& speed, const double t0, const double t_end, const double step, const dph::EphemerisRelease& de405);

Vector CalcSumY(const std::vector<Right>& fk);
Vector CalcSumYk(const std::vector<Right>& fk);
Vector CalcSumVy(const std::vector<Right>& fk);
Vector CalcSumVyk(const std::vector<Right>& fk);
void CalcVisiblePlanetAngle(const Vector& spececraft_r, const Vector& planet_r, const double planet_radius, const double time);
Vector CalcTimeDerivativeVelocity(const Vector& spacecraft_radius_vector, const Vector& moon_radius_vector, const Vector& sun_radius_vector, const double time);
Vector CalcGravitationalAcceleration(const double planet_gravitational_parameter, const Vector& planet_radius_vector, const Vector& spacecraft_radius_vector);
void CreateUV(Matrix& U, Matrix& V, const Vector& spacecraft_r, size_t size);
void FillUdVd(std::vector<Matrix>& U_d, std::vector<Matrix>& V_d, const Matrix& U, const Matrix& V, size_t d_size);

Vector Fk_rk4(const Vector& radius, const double time, const dph::EphemerisRelease& de405);

std::vector<Right> Fk(const Vector& radius, const Vector& speed, const double step, const double time, const dph::EphemerisRelease& de405);

Vector Abs(const Vector& v);
double MaxElement(const Vector& v1, const Vector& v2);

void PrintAngles(const double time, const double visible_angle, const double latitude, const double longitude);
void PrintAngles(std::ofstream& out, const double visible_angle, const double latitude, const double longitude);

template <std::size_t SizeT>
Vector CalcHarmonicAcceleratin(const double grav_param, const double equat_r, const Vector& spacecraft_r, size_t d_size,
				const std::array<std::array<double, SizeT>, SizeT>& Cnn, const std::array<std::array<double, SizeT>, SizeT>& Snn) {
	Vector harmonic_accel;
	size_t size = d_size + 1;
	std::vector<Matrix> U_d;
	std::vector<Matrix> V_d;
	for (size_t i = 0; i < 3; ++i) {
		U_d.emplace_back(d_size);
		V_d.emplace_back(d_size);
	}
	{
		Matrix U(size);
		Matrix V(size);
		CreateUV(U, V, spacecraft_r, size);
		FillUdVd(U_d, V_d, U, V, d_size);
	}
	for (size_t k = 0; k < 3; ++k) {
		double sum2 = 0.;
		double equat_r_pow = equat_r;
		for (int row = 1; row < d_size; ++row) {
			double sum1 = 0.;
			for (int col = 0; col <= row; ++col) {
				sum1 += Cnn.at(row).at(col) * U_d.at(k).at(row).at(col) + Snn.at(row).at(col) * V_d.at(k).at(row).at(col);
			}
			sum2 += sum1 * equat_r_pow;
			equat_r_pow *= equat_r;
		}
		harmonic_accel.PushBack(grav_param * sum2);
	}
	return harmonic_accel;
}



