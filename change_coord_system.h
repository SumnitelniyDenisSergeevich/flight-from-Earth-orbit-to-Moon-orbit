#pragma once

#include "vector.h"
#include "matrix.h"

std::pair<Vector, Vector>  J200toGCS(const Vector& r, const Vector& v, const double t);
Vector J200toGCS(const Vector& r, const double t);
Matrix Create_Matrix_P(const double T, const double T2, const double T3);
Matrix Create_Matrix_N(const double psi, const double eps0, const double eps);
Matrix Create_Matrix_M(const double xp, const double yp);
std::pair<Matrix, Matrix> Create_Matrix_S(const double t, const double psi, const double c_eps);
std::pair<double, double> CalcPsiEps(const double eps0, const double T, const double T2, const double T3);