#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>
#include <stdexcept>
#include <random>
#include <complex>

template<typename T>
class Matrix {
private:
    std::size_t rows;
    std::size_t cols;
    std::vector<std::vector<T>> data;

public:
    // Constructor with value initialization
    Matrix(std::size_t rows, std::size_t cols, T value = T()) 
        : rows(rows)
        , cols(cols)
        , data(rows, std::vector<T>(cols, value)) {}
    
    // Constructor with random values
    Matrix(std::size_t rows_, std::size_t cols_, T min, T max) 
        : rows(rows_)
        , cols(cols_)
        , data(rows_, std::vector<T>(cols_)) {
        std::random_device rd;
        std::mt19937 gen(rd());
        
        if constexpr (std::is_same_v<T, std::complex<double>>) {
            std::uniform_real_distribution<double> dis_real(min.real(), max.real());
            std::uniform_real_distribution<double> dis_imag(min.imag(), max.imag());
            
            for (std::size_t i = 0; i < rows; ++i) {
                for (std::size_t j = 0; j < cols; ++j) {
                    data[i][j] = std::complex<double>(dis_real(gen), dis_imag(gen));
                }
            }
        } else {
            std::uniform_real_distribution<double> dis(std::real(min), std::real(max));
            
            for (std::size_t i = 0; i < rows; ++i) {
                for (std::size_t j = 0; j < cols; ++j) {
                    data[i][j] = T(dis(gen));
                }
            }
        }
    }
    
    // Element access operator
    T& operator()(std::size_t i, std::size_t j) {
        return data[i][j];
    }
    
    const T& operator()(std::size_t i, std::size_t j) const {
        return data[i][j];
    }
    
    // Addition operator
    Matrix operator+(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols)
            throw std::invalid_argument("Matrix dimensions must match for addition.");
        Matrix result(rows, cols);
        for (std::size_t i = 0; i < rows; ++i)
            for (std::size_t j = 0; j < cols; ++j)
                result(i, j) = data[i][j] + other(i, j);
        return result;
    }
    
    // Subtraction operator
    Matrix operator-(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols)
            throw std::invalid_argument("Matrix dimensions must match for subtraction.");
        Matrix result(rows, cols);
        for (std::size_t i = 0; i < rows; ++i)
            for (std::size_t j = 0; j < cols; ++j)
                result(i, j) = data[i][j] - other(i, j);
        return result;
    }
    
    // Multiplication operator
    Matrix operator*(const Matrix& other) const {
        if (cols != other.rows)
            throw std::invalid_argument("Matrix dimensions must match for multiplication.");
        Matrix result(rows, other.cols, T());
        for (std::size_t i = 0; i < rows; ++i)
            for (std::size_t j = 0; j < other.cols; ++j)
                for (std::size_t k = 0; k < cols; ++k)
                    result(i, j) += data[i][k] * other(k, j);
        return result;
    }
    
    // Scalar multiplication
    Matrix operator*(const T& scalar) const {
        Matrix result(rows, cols);
        for (std::size_t i = 0; i < rows; ++i)
            for (std::size_t j = 0; j < cols; ++j)
                result(i, j) = data[i][j] * scalar;
        return result;
    }
    
    // Scalar division
    Matrix operator/(const T& scalar) const {
        Matrix result(rows, cols);
        for (std::size_t i = 0; i < rows; ++i)
            for (std::size_t j = 0; j < cols; ++j)
                result(i, j) = data[i][j] / scalar;
        return result;
    }
    
    // Trace calculation
    T Trace() const {
        if (rows != cols)
            throw std::invalid_argument("Trace is defined only for square matrices.");
        T trace = T();
        for (std::size_t i = 0; i < rows; ++i)
            trace += data[i][i];
        return trace;
    }
    
    // Calculate determinant for 3x3 matrix
    T determinant() const {
        if (rows != 3 || cols != 3) {
            throw std::invalid_argument("Determinant calculation supported only for 3x3 matrices");
        }
        
        return data[0][0] * (data[1][1] * data[2][2] - data[1][2] * data[2][1]) -
               data[0][1] * (data[1][0] * data[2][2] - data[1][2] * data[2][0]) +
               data[0][2] * (data[1][0] * data[2][1] - data[1][1] * data[2][0]);
    }
    
    // Friend function for commutative scalar multiplication
    template<typename U>
    friend Matrix<U> operator*(const U& scalar, const Matrix<U>& matrix);
};

// Commutative scalar multiplication
template<typename T>
Matrix<T> operator*(const T& scalar, const Matrix<T>& matrix) {
    return matrix * scalar;
}

#endif // MATRIX_HPP