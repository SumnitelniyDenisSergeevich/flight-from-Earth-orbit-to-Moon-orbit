#pragma once 

#include "vector.h"
#include "matrix.h"

#include "ephemeris.h"

#include <iostream> // delete
#include <fstream>

struct Right {
	Vector r = { 0.,0.,0. };
	Vector v = { 0.,0.,0. };
	Vector f = { 0.,0.,0. };// ускорение
	double time = 0.;
};

void RungeCutte8(Vector& radius, Vector& speed, double t, double step, const dph::EphemerisRelease& de405);
Vector CalcSumY(const std::vector<Right>& fk);
Vector CalcSumYk(const std::vector<Right>& fk);
Vector CalcSumVy(const std::vector<Right>& fk);
Vector CalcSumVyk(const std::vector<Right>& fk);
Vector CalcTimeDerivativeVelocity(const Vector& spacecraft_radius_vector, const Vector& moon_radius_vector, const Vector& sun_radius_vector, const double time);
Vector CalcGravitationalAcceleration(const double planet_gravitational_parameter, const Vector& planet_radius_vector, const Vector& spacecraft_radius_vector);
std::pair<Matrix,Matrix> CreateUV(const Vector& spacecraft_r, size_t size);
void FillUdVd(std::vector<Matrix>& U_d, std::vector<Matrix>& V_d, const Matrix& U, const Matrix& V, size_t d_size);
std::vector<Right> Fk(const Vector& radius, const Vector& speed, const double step, const double time, const dph::EphemerisRelease& de405);

Vector Abs(const Vector& v);
double MaxElement(const Vector& v1, const Vector& v2);

template <std::size_t SizeT>
Vector CalcHarmonicAcceleratin(const double grav_param, const double equat_r, const Vector& spacecraft_r, size_t d_size,
				const std::array<std::array<double, SizeT>, SizeT>& Cnn, const std::array<std::array<double, SizeT>, SizeT>& Snn) {
	Vector harmonic_accel(3);
	size_t size = d_size + 1;
	std::vector<Matrix> U_d = { Matrix(d_size), Matrix(d_size), Matrix(d_size) };
	std::vector<Matrix> V_d = { Matrix(d_size), Matrix(d_size), Matrix(d_size) };
	{
		std::ofstream out("output.txt"); // delete
		auto [U, V] = CreateUV(spacecraft_r, size);
		FillUdVd(U_d, V_d, U, V, d_size);
	}
	double sum1;
	double sum2;
	for (size_t k = 0; k < 3; ++k) {
		sum2 = 0.;
		for (int row = 1; row < d_size; ++row) {
			sum1 = 0.;
			for (int col = 0; col <= row; ++col) {// или <
				sum1 += Cnn.at(row).at(col) * U_d.at(k).at(row).at(col) + Snn.at(row).at(col) * V_d.at(k).at(row).at(col);
			}
			sum2 += equat_r * sum1;
		}
		harmonic_accel[k] = grav_param * sum2;
	}
	return harmonic_accel;
}



