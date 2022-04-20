#pragma once
#include <vector>
#include <initializer_list>
#include <fstream> //delete

#include "vector.h"

class Matrix {
public:
	Matrix(size_t size, double val = 0.); // создание квадратной матрицы
	Matrix(size_t rows, size_t cols, double val = 0.); // создание прямоугольной матрицы матрицы
	Matrix(std::initializer_list<std::vector<double>>& inputs); // создание матрицы по list прямоугольной матрицы матрицы

	std::vector<double>& operator[](size_t i);
	const std::vector<double>& at(size_t i) const;

	Matrix operator*(const Matrix& right_matrix) const;
	Vector operator*(const Vector& v) const;

	void Print(std::ofstream& out) const;
private:
	std::vector<std::vector<double>> matr_;
	size_t rows_, cols_;
};