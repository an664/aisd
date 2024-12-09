#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <stdexcept>
#include <random>
#include <complex>
#include <iostream>

template<typename T>
class Matrix {
private:
    std::size_t rows;
    std::size_t cols;
    T** data;

    // Static constant for floating-point comparison
    static constexpr double EPSILON = 1e-10;
    
    // Helper method for comparing values with precision
    static bool AreEqual(const T& a, const T& b) {
        if constexpr (std::is_same_v<T, std::complex<double>>) {
            return std::abs(a.real() - b.real()) < EPSILON && 
                   std::abs(a.imag() - b.imag()) < EPSILON;
        } else {
            return std::abs(a - b) < EPSILON;
        }
    }

    void AllocateMemory() {
        data = new T*[rows];
        for (std::size_t i = 0; i < rows; ++i) {
            data[i] = new T[cols];
        }
    }

    void FreeMemory() {
        if (data) {
            for (std::size_t i = 0; i < rows; ++i) {
                delete[] data[i];
            }
            delete[] data;
        }
    }

public:
    Matrix(std::size_t rows_, std::size_t cols_, T value = T()) 
        : rows(rows_), cols(cols_) {
        if (rows == 0 || cols == 0) {
            throw std::invalid_argument("Matrix dimensions cannot be zero");
        }
        
        AllocateMemory();
        for (std::size_t i = 0; i < rows; ++i) {
            for (std::size_t j = 0; j < cols; ++j) {
                data[i][j] = value;
            }
        }
    }

    Matrix(const Matrix& other) : rows(other.rows), cols(other.cols) {
        AllocateMemory();
        for (std::size_t i = 0; i < rows; ++i) {
            for (std::size_t j = 0; j < cols; ++j) {
                data[i][j] = other.data[i][j];
            }
        }
    }

    Matrix& operator=(const Matrix& other) {
        if (this != &other) {
            Matrix temp(other);
            std::swap(rows, temp.rows);
            std::swap(cols, temp.cols);
            std::swap(data, temp.data);
        }
        return *this;
    }

    ~Matrix() {
        FreeMemory();
    }

    Matrix(std::size_t rows_, std::size_t cols_, T min, T max) 
        : rows(rows_), cols(cols_) {
        if (rows == 0 || cols == 0) {
            throw std::invalid_argument("Matrix dimensions cannot be zero");
        }
        
        AllocateMemory();
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
    
    T& operator()(std::size_t i, std::size_t j) {
        return data[i][j];
    }
    
    const T& operator()(std::size_t i, std::size_t j) const {
        return data[i][j];
    }
    
    Matrix operator+(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols)
            throw std::invalid_argument("Matrix dimensions must match for addition.");
        Matrix result(rows, cols);
        for (std::size_t i = 0; i < rows; ++i)
            for (std::size_t j = 0; j < cols; ++j)
                result(i, j) = data[i][j] + other(i, j);
        return result;
    }
    
    Matrix operator-(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols)
            throw std::invalid_argument("Matrix dimensions must match for subtraction.");
        Matrix result(rows, cols);
        for (std::size_t i = 0; i < rows; ++i)
            for (std::size_t j = 0; j < cols; ++j)
                result(i, j) = data[i][j] - other(i, j);
        return result;
    }
    
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
    
    Matrix operator*(const T& scalar) const {
        Matrix result(rows, cols);
        for (std::size_t i = 0; i < rows; ++i)
            for (std::size_t j = 0; j < cols; ++j)
                result(i, j) = data[i][j] * scalar;
        return result;
    }
    
    Matrix operator/(const T& scalar) const {
        Matrix result(rows, cols);
        for (std::size_t i = 0; i < rows; ++i)
            for (std::size_t j = 0; j < cols; ++j)
                result(i, j) = data[i][j] / scalar;
        return result;
    }
    
    T Trace() const {
        if (rows != cols)
            throw std::invalid_argument("Trace is defined only for square matrices.");
        T trace = T();
        for (std::size_t i = 0; i < rows; ++i)
            trace += data[i][i];
        return trace;
    }
    
    T determinant() const {
        if (rows != 3 || cols != 3) {
            throw std::invalid_argument("Determinant calculation supported only for 3x3 matrices");
        }
        
        return data[0][0] * (data[1][1] * data[2][2] - data[1][2] * data[2][1]) -
               data[0][1] * (data[1][0] * data[2][2] - data[1][2] * data[2][0]) +
               data[0][2] * (data[1][0] * data[2][1] - data[1][1] * data[2][0]);
    }
    
    template<typename U>
    friend Matrix<U> operator*(const U& scalar, const Matrix<U>& matrix);
    
    std::size_t GetRows() const { return rows; }
    std::size_t GetCols() const { return cols; }

    // Equality operator
    bool operator==(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) {
            return false;
        }
        
        for (std::size_t i = 0; i < rows; ++i) {
            for (std::size_t j = 0; j < cols; ++j) {
                if (!AreEqual(data[i][j], other.data[i][j])) {
                    return false;
                }
            }
        }
        return true;
    }

    // Inequality operator
    bool operator!=(const Matrix& other) const {
        return !(*this == other);
    }

    // Stream output operator
    friend std::ostream& operator<<(std::ostream& os, const Matrix& matrix) {
        for (std::size_t i = 0; i < matrix.rows; ++i) {
            for (std::size_t j = 0; j < matrix.cols; ++j) {
                os << matrix.data[i][j];
                if (j < matrix.cols - 1) {
                    os << " ";
                }
            }
            if (i < matrix.rows - 1) {
                os << "\n";
            }
        }
        return os;
    }
};

template<typename T>
Matrix<T> operator*(const T& scalar, const Matrix<T>& matrix) {
    return matrix * scalar;
}

#endif // MATRIX_HPP