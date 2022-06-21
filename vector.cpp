#include "vector.h"
#include <iostream> //delete

using namespace std;

Vector::Vector(std::initializer_list<double> inputs) : vec_(inputs) {
}
Vector::Vector(size_t size, double val):vec_(size,val) {
}
Vector::Vector(size_t size) : vec_(size) {
}
Vector::Vector() {
}
Vector::Vector(Vector&& vec) noexcept {
	mod_.reset();
	vec_ = std::move(vec.vec_);
}
Vector::Vector(const Vector& vec) {
	mod_.reset();
	vec_ = vec.vec_;
}

const double Vector::at(size_t i) const {
	return vec_.at(i);
}

double& Vector::operator[](size_t i) {
	mod_.reset();
	return vec_[i];
}

Vector Vector::operator*(const Vector& right_vec) const {  // векторное произведение для vec.size() = 3;
	if (vec_.size() != right_vec.Size() || vec_.size() != 3) {
		throw invalid_argument("vectors size != 3");
	}
	Vector result(vec_.size());
	result[0] = vec_.at(1) * right_vec.vec_.at(2) - vec_.at(2) * right_vec.vec_.at(1);
	result[1] = -1 * (vec_.at(0) * right_vec.vec_.at(2) - vec_.at(2) * right_vec.vec_.at(0));
	result[2] = vec_.at(0) * right_vec.vec_.at(1) - vec_.at(1) * right_vec.vec_.at(0);
	return result;
}

Vector Vector::operator-(const Vector& right_vec) const {
	if (vec_.size() != right_vec.Size()) {
		throw invalid_argument("vectors has diffrent sizes");
	}
	Vector result(vec_.size());
	for (size_t i = 0; i < vec_.size(); ++i) {
		result[i] = vec_.at(i) - right_vec.vec_.at(i);
	}
	return result;
}

Vector Vector::operator+(const Vector& right_vec) const {
	if (vec_.size() != right_vec.Size()) {
		throw invalid_argument("vectors has diffrent sizes");
	}
	Vector result(vec_.size());
	for (size_t i = 0; i < vec_.size(); ++i) {
		result[i] = vec_.at(i) + right_vec.vec_.at(i);
	}
	return result;
}

Vector Vector::operator/(const double val) const {
	Vector result(vec_.size());
	for (size_t i = 0; i < vec_.size(); ++i) {
		result[i] = vec_.at(i) / val ;
	}
	return result;
}

Vector Vector::operator*(const double val) const {
	Vector result(vec_.size());
	for (size_t i = 0; i < vec_.size(); ++i) {
		result[i] = vec_.at(i) * val;
	}
	return result;
}

Vector operator*(const double val, const Vector& vec) {
	return vec * val;
}

double Vector::Module() const {
	if (mod_) {
		return *mod_;
	}
	double result = 0.;
	for (const double val : vec_) {
		result += val * val;
	}
	mod_ = sqrt(result);
	return *mod_;
}

void Vector::operator=(Vector&& vec) noexcept {
	mod_.reset();
	vec_ = std::move(vec.vec_);
}
void Vector::operator=(const Vector& vec) {
	mod_.reset();
	vec_ = vec.vec_;
}

double ScalarMultiply(const Vector& left_vec, const Vector& right_vec)
{
	if (left_vec.vec_.size() != right_vec.vec_.size()) {
		throw invalid_argument("diffrent ranges");
	}
	double result = 0.;
	for (size_t i = 0; i < left_vec.vec_.size(); ++i) {
		result += left_vec.vec_.at(i) * right_vec.vec_.at(i);
	}
	return result;
}

double СosineAngleBetweenVectors(const Vector& left_vec, const Vector& right_vec){
	return ScalarMultiply(left_vec, right_vec) / left_vec.Module() / right_vec.Module();
}

void Vector::PushBack(const double val) {
	vec_.push_back(val);
}

size_t Vector::Size() const {
	return vec_.size();
}

void Vector::Print() const {
	for (const double val : vec_) {
		cout << val << "  ";
	}
	cout << endl;
}
