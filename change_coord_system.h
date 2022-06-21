#pragma once

#include "vector.h"
#include "matrix.h"

std::pair<Vector, Vector>  J200toGCS(const Vector& r, const Vector& v, const double t);
Vector J200toGCS(const Vector& r, const double t);
Vector GCStoJ2000(const Vector& v, const double t);
Matrix Create_Matrix_P(const double T, const double T2, const double T3);
Matrix Create_Matrix_N(const double psi, const double eps0, const double eps);
Matrix Create_Matrix_M(const double xp, const double yp);
Matrix Create_Matrix_S(const double t, const double psi, const double c_eps);
std::pair<Matrix, Matrix> Create_Matrix_S_S1(const double t, const double psi, const double c_eps);
std::pair<double, double> CalcPsiEps(const double eps0, const double T, const double T2, const double T3);

double ConvertDateToUliane(const size_t d, const size_t mon, const size_t y, const size_t h = 12, const size_t m = 0, const size_t s = 0, const size_t ms = 0);
std::string ConvertUlianeToDate(const double uliane_date);