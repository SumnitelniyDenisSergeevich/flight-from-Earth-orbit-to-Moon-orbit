#include "integrate.h"
#include "constants.h"
#include "efemeris_enum.h"
#include "change_coord_system.h"


#include <iostream> // delete
#include <vector>
#include <iomanip>
#include <algorithm>
using namespace std;


void RungeCutte8(Vector& radius, Vector& speed, double t, double step, const dph::EphemerisRelease& de405) {
	const double TOL = 1.0e-12;
	vector<Right> fk;
	double t0 =  2459935.375;
	while (t0 + step / 86400.  <= t) {
		Vector TE_r(3);
		Vector TE_v(3);
		fk = std::move(Fk(radius, speed, step, t0, de405));
		Vector r_plus_dr = radius + CalcSumY(fk) * step;
		Vector r_plus_drk = radius + CalcSumYk(fk) * step;
		Vector v_plus_dv = speed + CalcSumVy(fk) * step;
		Vector v_plus_dvk = speed + CalcSumVyk(fk) * step;
		TE_r =Abs(r_plus_drk - r_plus_dr);
		TE_v =Abs(v_plus_dvk - v_plus_dv);
		double max_TE = MaxElement(TE_r, TE_v);
		if(max_TE > TOL) {
			step *= 0.7;
			continue;
		}
		else if (max_TE < 0.01679616 * TOL) {
			step *= 1.4;
			continue;
		}
		t0 += step / 86400.;// перевод в  JED
		radius = r_plus_dr;
		speed = v_plus_dv;
	}
	step = (t - t0)*86400.;
	fk = std::move(Fk(radius, speed, step, t0, de405));
	Vector r_plus_dr = radius + CalcSumY(fk) * step;
	Vector v_plus_dv = speed + CalcSumVy(fk) * step;
	t0 += step / 86400.;// перевод в  JED
	radius = r_plus_dr;
	speed = v_plus_dv;

	cout.precision(15);
	radius.Print();
	speed.Print();
	cout << t0 << endl;
}


Vector CalcSumY(const vector<Right>& fk) {
	Vector r_plus_dr(3);
	for (size_t k = 0; k < 3; ++k) {
		double sumy = 0.;
		for (size_t i = 0; i < 11; ++i) {
			sumy += fk.at(i).v.at(k) * rkf8_ce[i];
		}
		r_plus_dr[k] = sumy;
	}
	return r_plus_dr;
}

Vector CalcSumYk(const vector<Right>& fk) {
	Vector r_plus_drk(3);
	for (size_t k = 0; k < 3; ++k) {
		double sumyk = 0.;
		for (size_t i = 0; i < 13; ++i) {
			sumyk += fk.at(i).v.at(k) * rkf8_ce_cap[i];
		}
		r_plus_drk[k] = sumyk;
	}
	return r_plus_drk;
}

Vector CalcSumVy(const vector<Right>& fk) {
	Vector v_plus_dv(3);
	for (size_t k = 0; k < 3; ++k) {
		double sumv = 0.;
		for (size_t i = 0; i < 11; ++i) {
			sumv += fk.at(i).f.at(k) * rkf8_ce[i];
		}
		v_plus_dv[k] = sumv;
	}
	return v_plus_dv;
}

Vector CalcSumVyk(const vector<Right>& fk) {
	Vector v_plus_dv(3);
	for (size_t k = 0; k < 3; ++k) {
		double sumvk = 0.;
		for (size_t i = 0; i < 13; ++i) {
			sumvk += fk.at(i).f.at(k) * rkf8_ce_cap[i];
		}
		v_plus_dv[k] = sumvk;
	}
	return v_plus_dv;
}

void CalcVisiblePlanetAngle(const Vector& spececraft_r, const Vector& planet_r, const double planet_radius) {
	Vector dist = planet_r - spececraft_r;
	double visible_planet_angle = 2. * (atan(planet_radius / dist.Module())) /*/ pi * 180.*/;
	double latitude = atan(dist.at(2) / sqrt(dist.at(0) * dist.at(0) + dist.at(1) * dist.at(1)));
	double longitude = atan2(dist.at(1), dist.at(0));
	PrintAngles(visible_planet_angle, latitude, longitude);
}

void PrintAngles(const double visible_angle, const double latitude, const double longitude) {
	cout << visible_angle << "  " << latitude << "  " << longitude << endl;
}

void PrintAngles(std::ofstream& out, const double visible_angle, const double latitude, const double longitude) {
	out << visible_angle << "  " << latitude << "  " << longitude << endl;
}

vector<Right> Fk(const Vector& radius, const Vector& speed, const double step, const double time, const dph::EphemerisRelease& de405) {
	vector<Right> fk;
	fk.reserve(13);
	Right right;
	for (size_t i = 0; i < 13; ++i) {
		right.time = time + rkf8_alpha[i] * step /*/ 86400.*/; // из c в JED
		for (size_t j = 0; j < 3; ++j) {
			double sum_y = 0., sum_v = 0.;
			for (size_t k = 0; k < i; ++k) {
				sum_y += rkf8_beta[i][k] * fk.at(k).v[j];
				sum_v += rkf8_beta[i][k] * fk.at(k).f[j];
			}
			right.r[j] = radius.at(j) + sum_y * step;
			right.v[j] = speed.at(j) + sum_v * step;
		}
		//------------------------------------------------------------------------------------
		double resultArray[3]{};
		de405.calculateBody(dph::Calculate::POSITION, dph::Body::MOON, dph::Body::EARTH, right.time, resultArray);
		Vector moon_radius_vector = { resultArray[0], resultArray[1], resultArray[2] };
		de405.calculateBody(dph::Calculate::POSITION, dph::Body::SUN, dph::Body::EARTH, right.time, resultArray);
		Vector sun_radius_vector = { resultArray[0], resultArray[1], resultArray[2] };
		//------------------------------------------------------------------------------------
		right.f = CalcTimeDerivativeVelocity(right.r, moon_radius_vector, sun_radius_vector,right.time);
		fk.push_back(right);
	}
	return fk;
}

Vector CalcTimeDerivativeVelocity(const Vector& spacecraft_r, const Vector& moon_radius_vector, const Vector& sun_radius_vector, const double time) {
	Vector result(3);
	Vector grav_accel_moon = CalcGravitationalAcceleration(mum, spacecraft_r, moon_radius_vector);
	Vector grav_accel_sun = CalcGravitationalAcceleration(mus, spacecraft_r, sun_radius_vector); 

	Vector r_grinvich = J200toGCS(spacecraft_r, time);
	Vector harmonic_accel_earth_grinvich = CalcHarmonicAcceleratin(mu, Rz, r_grinvich, 33, Cnn, Snn);
	Vector harmonic_accel_earth = GCStoJ2000(harmonic_accel_earth_grinvich, time);
	//Vector harmonic_accel_moon = CalcHarmonicAcceleratin(mum, R_Moon, r_grinvich, 76, Cnn_moon, Snn_moon);

	double spacecraft_r_module3 = spacecraft_r.Module() * spacecraft_r.Module() * spacecraft_r.Module();
	result = spacecraft_r * (-mu) / (spacecraft_r_module3) +grav_accel_moon + grav_accel_sun +harmonic_accel_earth;
	return result;
}

Vector CalcGravitationalAcceleration(const double planet_gravitational_parameter, const Vector& spacecraft_radius_vector, const Vector& planet_radius_vector) {
	Vector temp = planet_radius_vector - spacecraft_radius_vector;
	double temp_module3 = temp.Module() * temp.Module() * temp.Module();
	double planet_module3 = planet_radius_vector.Module() * planet_radius_vector.Module() * planet_radius_vector.Module();
	return (planet_gravitational_parameter * ((temp / temp_module3)  - (planet_radius_vector / planet_module3)));
}

std::pair<Matrix,Matrix> CreateUV(const Vector& spacecraft_r, size_t size) {
	Matrix U(size);
	Matrix V(size);
	double spacecraft_mod_2 = spacecraft_r.Module() * spacecraft_r.Module();
	U[0][0] = 1. / spacecraft_r.Module();
	for (int k = 0; k + 1 < size; ++k) {
		U[k + 1][k + 1] = (2 * k + 1) / spacecraft_mod_2 * (spacecraft_r.at(0) * U.at(k).at(k) - spacecraft_r.at(1) * V.at(k).at(k));
		V[k + 1][k + 1] = (2 * k + 1) / spacecraft_mod_2 * (spacecraft_r.at(0) * V.at(k).at(k) + spacecraft_r.at(1) * U.at(k).at(k));
	}
	U[1][0] = spacecraft_r.at(2) / spacecraft_mod_2 * U.at(0).at(0);
	V[1][0] = spacecraft_r.at(2) / spacecraft_mod_2 * V.at(0).at(0);
	for (int row = 1; row + 1 < size; ++row) {
		for (int col = 0; col <= row; ++col) {
			U[row + 1][col] = (2 * row + 1) * spacecraft_r.at(2) * U.at(row).at(col) / (row - col + 1.)  / spacecraft_mod_2
				- (row + col) * U.at(row - 1).at(col) / (row - col + 1.) / spacecraft_mod_2;

			V[row + 1][col] = (2 * row + 1) * spacecraft_r.at(2) * V.at(row).at(col) /(row - col + 1.)  / spacecraft_mod_2
				- (row + col) * V.at(row - 1).at(col) / (row - col + 1.) / spacecraft_mod_2 ;
		}
	}
	return make_pair(U, V);
}

void FillUdVd(vector<Matrix>& U_d, vector<Matrix>& V_d, const Matrix& U, const Matrix& V, size_t d_size){
	for (int row = 0; row < d_size; ++row) {
		U_d[0][row][0] = -0.5 * U.at(row + 1).at(1) - 0.5 * U.at(row + 1).at(1);
		V_d[0][row][0] = -0.5 * V.at(row + 1).at(1) + 0.5 * V.at(row + 1).at(1);

		U_d[1][row][0] = -0.5 * V.at(row + 1).at(1) - 0.5 * V.at(row + 1).at(1);
		V_d[1][row][0] = 0.5 * U.at(row + 1).at(1) - 0.5 * U.at(row + 1).at(1);

		U_d[2][row][0] = -(row + 1) * U.at(row + 1).at(0);
		V_d[2][row][0] = -(row + 1) * V.at(row + 1).at(0);
	}
	for (int row = 0; row < d_size; ++row) {
		for (int col = 1; col < d_size; ++col) {
			U_d[0][row][col] = -0.5 * U.at(row + 1).at(col + 1) + 0.5 * (row - col + 2) * (row - col + 1) * U.at(row + 1).at(col - 1);
			V_d[0][row][col] = -0.5 * V.at(row + 1).at(col + 1) + 0.5 * (row - col + 2) * (row - col + 1) * V.at(row + 1).at(col - 1);

			U_d[1][row][col] = -0.5 * V.at(row + 1).at(col + 1) - 0.5 * (row - col + 2) * (row - col + 1) * V.at(row + 1).at(col - 1);
			V_d[1][row][col] = 0.5 * U.at(row + 1).at(col + 1) + 0.5 * (row - col + 2) * (row - col + 1) * U.at(row + 1).at(col - 1);

			U_d[2][row][col] = -(row - col + 1) * U.at(row + 1).at(col);
			V_d[2][row][col] = -(row - col + 1) * V.at(row + 1).at(col);
		}
	}
}

Vector Abs(const Vector& v) {
	Vector result(v.Size());
	for (size_t i = 0; i < v.Size(); ++i) {
		result[i] = std::fabs(v.at(i));
	}
	return result;
}

double MaxElement(const Vector& v1, const Vector& v2) {
	return std::max(*std::max_element(v1.begin(), v1.end()), *std::max_element(v2.begin(), v2.end()));
}