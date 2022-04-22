#include "matrix.h"
#include <iostream> //delete

using namespace std;

Matrix::Matrix(size_t size, double val): rows_(size), cols_(size) {
	matr_.resize(size);
	for (vector<double>& vec : matr_) {
		vec.resize(size, val);
	}
}

Matrix::Matrix(size_t rows, size_t cols, double val) : rows_(rows), cols_(cols) {
	matr_.resize(rows);
	for (vector<double>& vec : matr_) {
		vec.resize(cols, val);
	}
}

Matrix::Matrix(std::initializer_list<vector<double>>& inputs) : matr_(inputs), rows_(inputs.size()), cols_(matr_.at(0).size()) {
	
}

std::vector<double>& Matrix::operator[](size_t i) {
	if (i > rows_) {
		throw invalid_argument("out of range");
	}
	return matr_[i];
}

const std::vector<double>& Matrix::at(size_t i) const {
	if (i > rows_) {
		throw invalid_argument("out of range");
	}
	return matr_.at(i);
}

Matrix Matrix::operator*(const Matrix& right_matrix) const {
	if (cols_ != right_matrix.rows_) {
		throw invalid_argument("left_cols != right_rows");
	}
	Matrix result(this->rows_, right_matrix.cols_);

	for (size_t j = 0; j < result.rows_; ++j) {
		for (size_t i = 0; i < result.cols_; ++i) {
			for (size_t k = 0; k < cols_; ++k)
				result[j][i] += matr_[j][k] * right_matrix.matr_[k][i];
		}
	}
	return result;
}

Vector Matrix::operator*(const Vector& v) const {
	if (cols_ != v.Size()) {
		throw invalid_argument("left_cols != right_rows");
	}
	Vector result(rows_);

	for (size_t j = 0; j < rows_; ++j) {
		for (size_t k = 0; k < cols_; ++k)
			result[j] += matr_[j][k] * v.at(k);
	}
	return result;
}

Matrix Matrix::Transponir() const {
	Matrix result = *this;
	for (size_t row = 1; row < rows_; ++row) {
		for (size_t col = 0; col < row; ++col) {
			std::swap(result[row][col], result[col][row]);
		}
	}
	return result;
}

void Matrix::Print() const {
	for (const vector<double>& vec : matr_) {
		for (double val : vec) {
			cout << val << "  ";
		}
		cout << endl;
	}
}


void Matrix::Print(std::ofstream& out) const {
	for (const vector<double>& vec : matr_) {
		for (double val : vec) {
			out << val << "  ";
		}
		out << endl;
	}
}