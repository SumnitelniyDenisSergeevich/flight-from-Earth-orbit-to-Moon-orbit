#include "change_coord_system.h"
#include "constants.h"
#include <cmath>

double UTCtoTDB(const double utc_time) {
	return utc_time + (32.184 + 37) / 86400.;
}

double Smean(const double utc_time) {
	double JD0 = std::round(utc_time - 0.5) + 0.5;
	double T_ut1 = (JD0 - 2451545.0) / 36525;
	double T2_ut1 = T_ut1 * T_ut1;
	double T3_ut1 = T2_ut1 * T_ut1; 
	return 1.7533685592 + 628.3319706889 * T_ut1 + 6.7707139e-6 * T2_ut1 - 4.50876e-10 * T3_ut1 + 1.002737909350795 * 2 * pi * (utc_time - JD0);
}

std::pair<Vector,Vector> J200toGCS(const Vector& r, const Vector& v, const double t) {
	double T = (UTCtoTDB(t) - 2451545.0) / 36525.0;  //Ttdb число юлианских столетий
	double T2 = T * T;
	double T3 = T2 * T;
	double eps0 = 0.4090928042 - 0.2269655E-3 * T - 2.86E-9 * T2 + 8.80E-9 * T3;
	auto [psi, eps] = std::move(CalcPsiEps(eps0, T, T2, T3));
	Matrix P = std::move(Create_Matrix_P(T, T2, T3));
	Matrix N = std::move(Create_Matrix_N(psi, eps0, eps));
	double xp, yp; // обновляются ежедневно
	xp = 0;
	yp = 0;
	Matrix M = std::move(Create_Matrix_M(xp, yp));

	auto [S, S1] = std::move(Create_Matrix_S(t, psi, cos(eps)));

	Matrix MSNP = M * S * N * P;
	Matrix MS1NP = M * S1 * N * P;

	Vector r_result = MSNP * r;
	Vector v_result = MSNP * v + MS1NP * r;
	return { r_result, v_result };
}

Vector J200toGCS(const Vector& r, const double t) {
	double T = (UTCtoTDB(t) - 2451545.0) / 36525.0;  //Ttdb число юлианских столетий
	double T2 = T * T;
	double T3 = T2 * T;
	double eps0 = 0.4090928042 - 0.2269655E-3 * T - 2.86E-9 * T2 + 8.80E-9 * T3;
	auto [psi, eps] = std::move(CalcPsiEps(eps0, T, T2, T3));
	Matrix P = std::move(Create_Matrix_P(T, T2, T3));
	Matrix N = std::move(Create_Matrix_N(psi, eps0, eps));
	double xp, yp; // обновляются ежедневно
	xp = 0;
	yp = 0;
	Matrix M = std::move(Create_Matrix_M(xp, yp));

	auto [S, S1] = std::move(Create_Matrix_S(t, psi, cos(eps)));

	Matrix MSNP = M * S * N * P;
	Matrix MS1NP = M * S1 * N * P;

	Vector r_result = MSNP * r;
	return r_result;
}

std::pair<double, double> CalcPsiEps(const double eps0, const double T, const double T2, const double T3) {
	double Ml = 2.355548394 + (1325 * pi_2 + 3.470890873) * T + 1.517952E-4 * T2 + 3.103E-7 * T3;
	double Ms = 6.240035940 + (99 * pi_2 + 6.266610600) * T - 2.797400E-6 * T2 - 5.820E-8 * T3;
	double ul = 1.627901934 + (1342 * pi_2 + 1.431476084) * T - 6.427170E-5 * T2 + 5.340E-8 * T3;
	double Ds = 5.198469514 + (1236 * pi_2 + 5.360106500) * T - 3.340860E-5 * T2 + 9.220E-8 * T3;
	double Uzl = 2.182438624 - (5 * pi_2 + 2.341119397) * T + 3.614290E-5 * T2 + 3.880E-8 * T3;

	double eps = 0, psi = 0;
	for (int i = 0; i <= 105; i++) {
		psi += (coef_A[i] + coef_B[i] * T) * sin(coef_a1[i] * Ml + coef_a2[i] * Ms + coef_a3[i] * ul + coef_a4[i] * Ds + coef_a5[i] * Uzl);
		eps += (coef_C[i] + coef_D[i] * T) * cos(coef_a1[i] * Ml + coef_a2[i] * Ms + coef_a3[i] * ul + coef_a4[i] * Ds + coef_a5[i] * Uzl);
	}
	psi = psi * 1E-4 / 3600. / 180. * pi;// перевод в радианты
	eps = eps * 1E-4 / 3600. / 180. * pi;// перевод в радианты
	eps += eps0;
	return { psi,eps };
}

Matrix Create_Matrix_P(const double T, const double T2, const double T3) {
	double ksi = 0.011180860 * T + 1.464E-6 * T2 + 8.70E-8 * T3;
	double zeta = 0.011180860 * T + 5.308E-6 * T2 + 8.90E-8 * T3;
	double teta = 0.009717173 * T - 2.068E-6 * T2 - 2.02E-7 * T3;

	double c_ksi = cos(ksi);
	double s_ksi = sin(ksi);
	double c_zeta = cos(zeta);
	double s_zeta = sin(zeta);
	double c_teta = cos(teta);
	double s_teta = sin(teta);

	Matrix P(3);

	P[0][0] = c_ksi * c_zeta * c_teta - s_ksi * s_zeta;
	P[0][1] = -s_ksi * c_zeta * c_teta - c_ksi * s_zeta;
	P[0][2] = -c_zeta * s_teta;

	P[1][0] = c_ksi * s_zeta * c_teta + s_ksi * c_zeta;
	P[1][1] = -s_ksi * s_zeta * c_teta + c_ksi * c_zeta;
	P[1][2] = -s_zeta * s_teta;

	P[2][0] = c_ksi * s_teta;
	P[2][1] = -s_ksi * s_teta;
	P[2][2] = c_teta;
	return P;
}


Matrix Create_Matrix_N( const double psi, const double eps0, const double eps) {
	double c_psi = cos(psi);
	double s_psi = sin(psi);
	double c_eps = cos(eps);
	double s_eps = sin(eps);
	double c_eps0 = cos(eps0);
	double s_eps0 = sin(eps0);

	Matrix N(3);

	N[0][0] = c_psi;
	N[0][1] = -s_psi * c_eps0;
	N[0][2] = -s_psi * s_eps0;

	N[1][0] = s_psi * c_eps;
	N[1][1] = c_psi * c_eps * c_eps0 + s_eps * s_eps0;
	N[1][2] = c_psi * c_eps * s_eps0 - s_eps * c_eps0;

	N[2][0] = s_psi * s_eps;
	N[2][1] = c_psi * s_eps * c_eps0 - c_eps * s_eps0;
	N[2][2] = c_psi * s_eps * s_eps0 + c_eps * c_eps0;
	return N;
}

Matrix Create_Matrix_M(const double xp, const double yp) {
	Matrix M(3);

	M[0][0] = cos(xp);
	M[0][1] = 0;
	M[0][2] = sin(xp);

	M[1][0] = sin(xp) * sin(yp);
	M[1][1] = cos(yp);
	M[1][2] = -cos(xp) * sin(yp);

	M[2][0] = -sin(xp) * cos(yp);
	M[2][1] = sin(yp);
	M[2][2] = cos(xp) * cos(yp);
	return M;
}


std::pair<Matrix,Matrix> Create_Matrix_S(const double t, const double psi, const double c_eps) {
	double w = -7.292115e-5;
	double s = Smean(t) + psi * c_eps;

	double c_l = cos(s);
	double s_l = sin(s);
	Matrix S(3);
	S[0][0] = c_l;
	S[0][1] = s_l;
	S[0][2] = 0;	

	S[1][0] = -s_l;
	S[1][1] = c_l;
	S[1][2] = 0;

	S[2][0] = 0;  
	S[2][1] = 0;   
	S[2][2] = 1;
	Matrix S1(3);

	S1[0][0] = w * s_l;
	S1[0][1] = -w * c_l; 
	S1[0][2] = 0;

	S1[1][0] = w * c_l; 
	S1[1][1] = w * s_l; 
	S1[1][2] = 0;

	S1[2][0] = 0;    
	S1[2][1] = 0;    
	S1[2][2] = 0;

	return { S, S1 };
}