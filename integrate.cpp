#include "integrate.h"
#include "constants.h"
#include "efemeris_enum.h"
#include "change_coord_system.h"

#include <iostream> // delete
#include <fstream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <random>
#include <cmath>

using namespace std;

ofstream runge_cutte8("rughe_cutte8_EARTH_TO_MOON_ORBIT_END.csv"s);
ofstream check_connection("check_connection_EARTH_TO_MOON_ORBIT_END.csv"s);
ofstream optical_measurements("optical_measurements_EARTH_TO_MOON_ORBIT_END.csv"s);

void CheckConnection(const Vector& r, const Vector& radiotelescope_r, const double t) {
	int i = 0;
	static vector<bool> start_flag(4, false);
	static vector<double> start_communication(4, 0.);
	Vector r_gcs = J200toGCS(r, t);
	const vector<string> KIP_id = { "Уссурийск"s, "Байконур"s,"Медвежьи озера","Евпатория" };
	//for (const auto& [key, vec] : CONST_NIP_COOR) {
	for(const string id : KIP_id){
		double cos_alpha = СosineAngleBetweenVectors(r_gcs - CONST_NIP_COOR.at(id), CONST_NIP_COOR.at(id));
		double alpha = acos(cos_alpha);
		double beta = (pi / 2 - alpha) * rtd;
		if (beta >= 7.) {
			if (!start_flag[i]) {
				start_flag[i] = true;
				start_communication[i] = t;
			}
		}
		else {
			if (start_flag[i]) {
				start_flag[i] = false;
				check_connection << id <<';' << ConvertUlianeToDate(start_communication[i]) << ';' << to_string(start_communication[i])
						<< ';' << ConvertUlianeToDate(t) << ';' << to_string(t) << endl;
			}
		}
		++i;
	}
}

//void CheckConnection(const Vector& r, const Vector& radiotelescope_r, const double t) {
//	double cos_alpha = СosineAngleBetweenVectors(r - radiotelescope_r, radiotelescope_r);
//	double alpha = acos(cos_alpha);
//	double beta = (pi/2 - alpha) *rtd;
//	static bool start_flag = false;
//	static double start_communication = 0.;
//	if (beta >= 7.) {
//		if (!start_flag) {
//			start_flag = true;
//			start_communication = t;
//		}
//	}
//	else {
//		if (start_flag) {
//			start_flag = false;
//			check_connection << start_communication*86400 << ';' << t * 86400 << endl;
//		}
//	}
//}

void RungeCutte8(Vector& radius, Vector& speed, double t0, double t_end, double step, const dph::EphemerisRelease& de405) {
	optical_measurements << "Угловой диаметр" << ';' << "широта" << ';' << "долгота"  << ';' << "Время" << endl;
	runge_cutte8 << "x KA" << ';' << "y KA" << ';' << "z KA" << ';' << "Время" << ';' << "x MOON" << ';' << "y MOON" << ';' << "z MOON" << endl;
	check_connection << "Место связи"s << ';' << "Начало коммуникации" << ';' << "Конец коммуникации" << endl;

	const double TOL = 1.0e-8;
	vector<Right> fk;
	bool flag_angle = false;
	int z = 0;
	int k = 0;
	double step_rk4;
	double communication_step = 1.;
	//{
	//	double resultArray[3]{};
	//	de405.calculateBody(dph::Calculate::POSITION, dph::Body::MOON, dph::Body::EARTH, t0, resultArray);
	//	Vector moon_radius_vector = { resultArray[0], resultArray[1], resultArray[2] };

	//	runge_cutte8 << radius.at(0) << ';' << radius.at(1) << ';' << radius.at(2) << ';' << ConvertUlianeToDate(t0) << ';' <<
	//		moon_radius_vector.at(0) << ';' << moon_radius_vector.at(1) << ';' << moon_radius_vector.at(2) << endl;
	//}
	Vector given_vector{ CONST_NIP_COOR.at("Химки")};//задаётся
	CheckConnection(radius, given_vector, t0);
	while (t0 + step / 86400.  <= t_end) {
		fk = Fk(radius, speed, step, t0, de405);
		Vector r_plus_dr = radius + CalcSumY(fk) * step;
		Vector r_plus_drk = radius + CalcSumYk(fk) * step;
		Vector v_plus_dv = speed + CalcSumVy(fk) * step;
		Vector v_plus_dvk = speed + CalcSumVyk(fk) * step;
		Vector TE_r = Abs(r_plus_drk - r_plus_dr);
		Vector TE_v = Abs(v_plus_dvk - v_plus_dv);
		double max_TE = MaxElement(TE_r, TE_v);
		if (z == 3) {

		}
		else if (max_TE > TOL) {
			++z;
			step *= 0.7;
			continue;
		}
		else if (max_TE < 0.01679616 * TOL) {
			++z;
			step *= 1.4;
			continue;
		}
		//---------------------------------------------------------angles-------------------------------------------
		double d_sec = (t0 * 24. - floor(t0 * 24.)) * 3600;
		if (d_sec + step > 3600 || flag_angle) { // выходит за рамки часа, тогда надо 3 минуты поинтегрировать
			if (!flag_angle) {
				step_rk4 = 3600. - d_sec;
			}
			flag_angle = true;
			while (k != 180 && step_rk4 < step) {
				Right rp = RungeCutte4(radius, speed, t0, t0 + step_rk4 / 86400., step_rk4, de405);

				double resultArray[3]{};
				de405.calculateBody(dph::Calculate::POSITION, dph::Body::MOON, dph::Body::EARTH, rp.time, resultArray);//change end_time -> rp.time
				Vector moon_radius_vector = { resultArray[0], resultArray[1], resultArray[2] };
				if (rp.r.Module() > 40000 && (moon_radius_vector - rp.r).Module() > 10000) {
					CalcVisiblePlanetAngle(rp.r, moon_radius_vector, R_Moon, rp.time);
				}

				step_rk4 += 1.;
				++k;
			}
			if (k == 180) {
				flag_angle = false;
				k = 0;
			}
			else {
				step_rk4 -= step;
			}
		}
		//---------------------------------------------------------angles-------------------------------------------
		
		//---------------------------------------------------------radio-------------------------------------------
		while (communication_step < step) {
			Right rp = RungeCutte4(radius, speed, t0, t0 + communication_step / 86400., communication_step, de405);
			CheckConnection(rp.r, given_vector, t0);
			++communication_step;
		}
		communication_step -= step;
		//---------------------------------------------------------radio-------------------------------------------
		z = 0;
		t0 += step / 86400.;// перевод в  JD
		cout << t0 << '\t' << step << endl;
		radius = std::move(r_plus_dr);
		speed = std::move(v_plus_dv);
		/*{
			double resultArray[3]{};
			de405.calculateBody(dph::Calculate::POSITION, dph::Body::MOON, dph::Body::EARTH, t0, resultArray);
			Vector moon_radius_vector = { resultArray[0], resultArray[1], resultArray[2] };
			runge_cutte8 << radius.at(0) << ';' << radius.at(1) << ';' << radius.at(2) << ';' << ConvertUlianeToDate(t0) << ';' <<
				moon_radius_vector.at(0) << ';' << moon_radius_vector.at(1) << ';' << moon_radius_vector.at(2) << endl;
		}*/
	}
	step = (t_end - t0)*86400.;
	fk = Fk(radius, speed, step, t0, de405);
	Vector r_plus_dr = radius + CalcSumY(fk) * step;
	Vector v_plus_dv = speed + CalcSumVy(fk) * step;
	t0 += step / 86400.;// перевод в  JED
	radius = std::move(r_plus_dr);
	speed = std::move(v_plus_dv);
	CheckConnection(radius, { 0,0,0 }, t0);

	cout.precision(15);
	radius.Print();
	speed.Print();
	cout << t0 << endl;
}

Right RungeCutte4Step(const Right& rp, const double step, const dph::EphemerisRelease& de405) {
	Right ans;
	Vector k1 = rp.f * step;
	Vector v1 = rp.v * step;
	ans.time = rp.time + step /86400. / 2.;
	ans.r = rp.r + v1 / 2.;
	ans.v = rp.v + k1 / 2.;
	ans.f = Fk_rk4(ans.r, ans.time, de405);

	Vector k2 = ans.f * step;
	Vector v2 = ans.v * step;

	ans.r = rp.r + v2 / 2.;
	ans.v = rp.v + k2 / 2.;
	ans.f = Fk_rk4(ans.r, ans.time, de405);

	Vector k3 = ans.f * step;
	Vector v3 = ans.v * step;

	ans.time = rp.time + step / 86400.;
	ans.r = rp.r + v3;
	ans.v = rp.v + k3;

	Vector k4 = ans.f * step;
	Vector v4 = ans.v * step;

	ans.r = rp.r + (v1 + 2. * v2 + 2. * v3 + v4) / 6.;
	ans.v = rp.v + (k1 + 2. * k2 + 2. * k3 + k4) / 6.;
	ans.f = Fk_rk4(ans.r, ans.time, de405);
	return ans;
}

Right RungeCutte4(const Vector& radius, const Vector& speed, const double t_start, const double t_end, const double step, const dph::EphemerisRelease& de405) {
	double t0 = t_start;
	Right rp = { radius, speed, Fk_rk4(radius, t0, de405), t0 };
	/*while (t0 + step/86400. < t_end) { // просто шаг, без цикла, так как нет начального шага.
		rp = RungeCutte4Step(rp, step, de405);
		t0 += step / 86400.;
	}
	double end_step = (t_end - t0) * 86400.;*/
	return RungeCutte4Step(rp, step, de405);
}


Vector CalcSumY(const vector<Right>& fk) {
	Vector r_plus_dr;
	for (size_t k = 0; k < 3; ++k) {
		double sumy = 0.;
		for (size_t i = 0; i < 11; ++i) {
			sumy += fk.at(i).v.at(k) * rkf8_ce[i];
		}
		r_plus_dr.PushBack(sumy);
	}
	return r_plus_dr;
}

Vector CalcSumYk(const vector<Right>& fk) {
	Vector r_plus_drk;
	for (size_t k = 0; k < 3; ++k) {
		double sumyk = 0.;
		for (size_t i = 0; i < 13; ++i) {
			sumyk += fk.at(i).v.at(k) * rkf8_ce_cap[i];
		}
		r_plus_drk.PushBack(sumyk);
	}
	return r_plus_drk;
}

Vector CalcSumVy(const vector<Right>& fk) {
	Vector v_plus_dv;
	for (size_t k = 0; k < 3; ++k) {
		double sumv = 0.;
		for (size_t i = 0; i < 11; ++i) {
			sumv += fk.at(i).f.at(k) * rkf8_ce[i];
		}
		v_plus_dv.PushBack(sumv);
	}
	return v_plus_dv;
}

Vector CalcSumVyk(const vector<Right>& fk) {
	Vector v_plus_dv;
	for (size_t k = 0; k < 3; ++k) {
		double sumvk = 0.;
		for (size_t i = 0; i < 13; ++i) {
			sumvk += fk.at(i).f.at(k) * rkf8_ce_cap[i];
		}
		v_plus_dv.PushBack(sumvk);
	}
	return v_plus_dv;
}

void CalcVisiblePlanetAngle(const Vector& spececraft_r, const Vector& planet_r, const double planet_radius, const double time) {
	Vector dist = planet_r - spececraft_r;
	
	//-----------------------------------noise----------------------------------------------
	//std::random_device rd{};
	//std::mt19937 gen{ rd() };

	//std::normal_distribution<> alpha{ 1.1,2.5 }; // в угловых секундах
	//std::normal_distribution<> betta{ 4.5,7.4 };

	//-----------------------------------noise----------------------------------------------
	double visible_planet_angle = 2. * (atan(planet_radius / dist.Module())) * rtd /* + betta(gen) / 3600. */ ; // переводим в градусы
	double latitude = atan(dist.at(2) / sqrt(dist.at(0) * dist.at(0) + dist.at(1) * dist.at(1))) * rtd /*+ alpha(gen) / 3600.*/;//широта - прямое склонение
	double longitude = atan2(dist.at(1), dist.at(0)) * rtd /*+ betta(gen) / 3600.*/;// долгота - прямое восхождение
	optical_measurements << to_string(visible_planet_angle) << ';' << to_string(latitude) << ';' << to_string(longitude)
		<< ';' << ConvertUlianeToDate(time) << ';' << to_string(time) << endl;
}

void PrintAngles(const double time, const double visible_angle, const double latitude, const double longitude) {
	double koef = 1 / pi * 180.;
	cout << time << "  " << visible_angle * koef << "  " << latitude * koef << "  " << longitude * koef << endl;
}

void PrintAngles(std::ofstream& out, const double visible_angle, const double latitude, const double longitude) {
	out << visible_angle << "  " << latitude << "  " << longitude << endl;
}

Vector Fk_rk4(const Vector& radius, const double time, const dph::EphemerisRelease& de405) {
	Vector fk;
	double resultArray[3]{};
	de405.calculateBody(dph::Calculate::POSITION, dph::Body::MOON, dph::Body::EARTH, time, resultArray);
	Vector moon_radius_vector = { resultArray[0], resultArray[1], resultArray[2] };
	de405.calculateBody(dph::Calculate::POSITION, dph::Body::SUN, dph::Body::EARTH, time, resultArray);
	Vector sun_radius_vector = { resultArray[0], resultArray[1], resultArray[2] };
	fk = CalcTimeDerivativeVelocity(radius, moon_radius_vector, sun_radius_vector, time);
	return fk;
}

vector<Right> Fk(const Vector& radius, const Vector& speed, const double step, const double time, const dph::EphemerisRelease& de405) {
	vector<Right> fk;
	fk.reserve(13);
	for (size_t i = 0; i < 13; ++i) {
		Right right;
		right.time = time + rkf8_alpha[i] * step;
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
		fk.push_back(std::move(right));
	}
	return fk;
}

Vector CalcTimeDerivativeVelocity(const Vector& spacecraft_r, const Vector& moon_radius_vector, const Vector& sun_radius_vector, const double time) {
	Vector grav_accel_moon = CalcGravitationalAcceleration(mum, spacecraft_r, moon_radius_vector);
	Vector grav_accel_sun = CalcGravitationalAcceleration(mus, spacecraft_r, sun_radius_vector);
	Vector r_grinvich = J200toGCS(spacecraft_r, time);
	Vector harmonic_accel_earth_grinvich = CalcHarmonicAcceleratin(mu, Rz, r_grinvich, 33, Cnn, Snn);
	Vector harmonic_accel_earth = GCStoJ2000(harmonic_accel_earth_grinvich, time);
	//Vector harmonic_accel_moon = CalcHarmonicAcceleratin(mum, R_Moon, r_grinvich, 76, Cnn_moon, Snn_moon);

	return spacecraft_r * (-mu) / (spacecraft_r.Module() * spacecraft_r.Module() * spacecraft_r.Module())+grav_accel_moon + grav_accel_sun + harmonic_accel_earth;
}

Vector CalcGravitationalAcceleration(const double planet_gravitational_parameter, const Vector& spacecraft_radius_vector, const Vector& planet_radius_vector) {
	Vector temp = planet_radius_vector - spacecraft_radius_vector;
	double temp_module3 = temp.Module() * temp.Module() * temp.Module();
	double planet_module3 = planet_radius_vector.Module() * planet_radius_vector.Module() * planet_radius_vector.Module();
	return (planet_gravitational_parameter * ((temp / temp_module3)  - (planet_radius_vector / planet_module3)));
}

void CreateUV(Matrix& U, Matrix& V, const Vector& spacecraft_r, size_t size) {
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