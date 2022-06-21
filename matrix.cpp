#include "matrix.h"
#include <iostream> //delete

using namespace std;

Matrix::Matrix() : rows_(0), cols_(0) {
}

Matrix::Matrix(size_t size) : rows_(size), cols_(size) {
	matr_.resize(size);
	for (vector<double>& vec : matr_) {
		vec.resize(size);
	}
}

Matrix::Matrix(size_t rows, size_t cols) : rows_(rows), cols_(cols) {
	matr_.resize(rows);
	for (vector<double>& vec : matr_) {
		vec.resize(cols);
	}
}

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

Matrix::Matrix(Matrix&& matr) noexcept {
	cols_ = matr.cols_;
	rows_ = matr.rows_;
	matr_ = std::move(matr.matr_);
}

Matrix::Matrix(const Matrix& matr) {
	cols_ = matr.cols_;
	rows_ = matr.rows_;
	matr_ = matr.matr_;
}

std::vector<double>& Matrix::operator[](size_t i) {
	return matr_[i];
}

const std::vector<double>& Matrix::at(size_t i) const {
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

void Matrix::operator=(Matrix&& tmp) noexcept {
	cols_ = tmp.cols_;
	rows_ = tmp.rows_;
	matr_ = std::move(tmp.matr_);
}

void Matrix::operator=(const Matrix& tmp) {
	cols_ = tmp.cols_;
	rows_ = tmp.rows_;
	matr_ = tmp.matr_;
}

const Matrix& Matrix::Transponir() {
	for (size_t row = 1; row < rows_; ++row) {
		for (size_t col = 0; col < row; ++col) {
			std::swap(matr_[row][col], matr_[col][row]);
		}
	}
	return *this;
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