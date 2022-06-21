#pragma once
#include <vector>
#include <initializer_list>
#include <optional>

class Vector {
public:
	Vector(std::initializer_list<double> inputs);
	Vector(size_t size, double val);
	Vector(size_t size);
	Vector();
	Vector(Vector&& vec) noexcept;
	Vector(const Vector& vec);

	const double at(size_t i) const;
	double& operator[](size_t i);
	Vector operator*(const Vector& right_vec) const; // векторное произведение
	Vector operator-(const Vector& right_vec) const;
	Vector operator+(const Vector& right_vec) const;
	Vector operator/(const double val) const;
	Vector operator*(const double val) const;
	friend Vector operator*(const double val, const Vector& vec);
	void operator=(Vector&& tmp) noexcept;
	void operator=(const Vector& tmp);

	double Module() const;
	friend double ScalarMultiply(const Vector& left_vec, const Vector& right_vec);				//скалярное произведение
	friend double СosineAngleBetweenVectors(const Vector& left_vec, const Vector& right_vec);
	void PushBack(const double val);

	size_t Size() const;
	void Print() const;
	
	inline std::vector<double>::iterator begin() {return vec_.begin();}
	inline std::vector<double>::const_iterator begin() const {return vec_.begin();}
	inline const std::vector<double>::const_iterator cbegin() const {return vec_.cbegin();}
	inline std::vector<double>::iterator end() {return vec_.end();}
	inline std::vector<double>::const_iterator end() const {return vec_.end();}
	inline std::vector<double>::const_iterator cend() const {return vec_.cend();}
private:
	std::vector<double> vec_;
	mutable std::optional<double> mod_ = std::nullopt;
};

double ScalarMultiply(const Vector& left_vec, const Vector& right_vec);

double СosineAngleBetweenVectors(const Vector& left_vec, const Vector& right_vec);
