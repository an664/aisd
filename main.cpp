#include <iostream>
#include <cmath>
#include <limits>
#include <complex>
#include "matrix.hpp"

// Clear input buffer
void ClearInputBuffer() {
    std::cin.clear();
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

// Safe number input
template<typename T>
T GetNumber(const std::string& prompt) {
    T value;
    while (true) {
        std::cout << prompt;
        if (std::cin >> value) {
            ClearInputBuffer();
            return value;
        }
        std::cout << "Input error. Please enter a number.\n";
        ClearInputBuffer();
    }
}

// Vector input
template<typename T>
std::vector<T> InputVector(const std::string& vectorName) {
    std::vector<T> vec(3);
    std::cout << "\nEnter vector " << vectorName << ":\n";
    vec[0] = GetNumber<T>("x = ");
    vec[1] = GetNumber<T>("y = ");
    vec[2] = GetNumber<T>("z = ");
    return vec;
}

template<typename T>
void PrintVector(const std::vector<T>& vec, const std::string& vectorName) {
    std::cout << vectorName << " = (" << vec[0] << ", " << vec[1] << ", " << vec[2] << ")\n";
}

template<typename T>
bool AreVectorsCoplanar(const std::vector<T>& a, const std::vector<T>& b, const std::vector<T>& c) {
    Matrix<T> m(3, 3);
    
    for (std::size_t i = 0; i < 3; ++i) {
        m(i, 0) = a[i];
        m(i, 1) = b[i];
        m(i, 2) = c[i];
    }
    
    const double epsilon = 1e-10;
    return std::abs(std::abs(m.determinant())) < epsilon;
}

// Специализация GetNumber для комплексных чисел
template<>
std::complex<double> GetNumber(const std::string& prompt) {
    double real, imag;
    while (true) {
        std::cout << prompt << " (real imag): ";
        if (std::cin >> real >> imag) {
            ClearInputBuffer();
            return std::complex<double>(real, imag);
        }
        std::cout << "Input error. Please enter two numbers for real and imaginary parts.\n";
        ClearInputBuffer();
    }
}

template<typename T>
void PrintMatrix(const Matrix<T>& m, const std::string& name) {
    std::cout << "\nMatrix " << name << ":\n";
    for (std::size_t i = 0; i < m.GetRows(); ++i) {
        for (std::size_t j = 0; j < m.GetCols(); ++j) {
            std::cout << m(i, j) << " ";
        }
        std::cout << "\n";
    }
}

void ShowInputChoiceMenu() {
    std::cout << "\nChoose input method:\n";
    std::cout << "1. Manual input\n";
    std::cout << "2. Random values\n";
    std::cout << "Your choice: ";
}

template<typename T>
Matrix<T> GetRandomMatrix(std::size_t rows, std::size_t cols) {
    if constexpr (std::is_same_v<T, std::complex<double>>) {
        return Matrix<T>(rows, cols, std::complex<double>(-10.0, -10.0), std::complex<double>(10.0, 10.0));
    } else {
        return Matrix<T>(rows, cols, -10.0, 10.0);
    }
}

template<typename T>
std::vector<T> GetRandomVector() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(-10.0, 10.0);
    
    std::vector<T> vec(3);
    if constexpr (std::is_same_v<T, std::complex<double>>) {
        for (std::size_t i = 0; i < 3; ++i) {
            vec[i] = std::complex<double>(dis(gen), dis(gen));
        }
    } else {
        for (std::size_t i = 0; i < 3; ++i) {
            vec[i] = dis(gen);
        }
    }
    return vec;
}

void RealMatrixOperations() {
    ShowInputChoiceMenu();
    char inputChoice;
    std::cin >> inputChoice;
    ClearInputBuffer();

    Matrix<double> m1(2, 2), m2(2, 2);
    
    if (inputChoice == '1') {
        // Manual input
        std::cout << "\nEnter first matrix elements:\n";
        for (std::size_t i = 0; i < 2; ++i) {
            for (std::size_t j = 0; j < 2; ++j) {
                m1(i, j) = GetNumber<double>("Enter element [" + std::to_string(i) + "][" + std::to_string(j) + "]: ");
            }
        }
        std::cout << "\nEnter second matrix elements:\n";
        for (std::size_t i = 0; i < 2; ++i) {
            for (std::size_t j = 0; j < 2; ++j) {
                m2(i, j) = GetNumber<double>("Enter element [" + std::to_string(i) + "][" + std::to_string(j) + "]: ");
            }
        }
    } else {
        // Random values
        m1 = GetRandomMatrix<double>(2, 2);
        m2 = GetRandomMatrix<double>(2, 2);
    }

    // Display matrices
    PrintMatrix(m1, "M1");
    PrintMatrix(m2, "M2");

    // Addition
    auto sum = m1 + m2;
    PrintMatrix(sum, "M1 + M2");

    // Subtraction
    auto diff = m1 - m2;
    PrintMatrix(diff, "M1 - M2");

    // Matrix multiplication
    auto prod = m1 * m2;
    PrintMatrix(prod, "M1 * M2");

    // Scalar multiplication
    double scalar = 2.0;
    auto scaledM1 = m1 * scalar;
    auto scaledM2 = scalar * m2; // Commutability check
    PrintMatrix(scaledM1, "M1 * 2.0");
    PrintMatrix(scaledM2, "2.0 * M2");

    // Division by scalar
    auto dividedM1 = m1 / scalar;
    PrintMatrix(dividedM1, "M1 / 2.0");

    // Trace
    std::cout << "Trace of M1: " << m1.Trace() << std::endl;
    std::cout << "Trace of M2: " << m2.Trace() << std::endl;
}

template<typename T>
void CheckVectorCoplanarity() {
    ShowInputChoiceMenu();
    char inputChoice;
    std::cin >> inputChoice;
    ClearInputBuffer();

    std::vector<T> a, b, c;
    
    if (inputChoice == '1') {
        // Manual input
        a = InputVector<T>("a");
        b = InputVector<T>("b");
        c = InputVector<T>("c");
    } else {
        // Random values
        a = GetRandomVector<T>();
        b = GetRandomVector<T>();
        c = GetRandomVector<T>();
    }

    std::cout << "\nVectors:\n";
    PrintVector(a, "a");
    PrintVector(b, "b");
    PrintVector(c, "c");

    if (AreVectorsCoplanar(a, b, c)) {
        std::cout << "\nVectors are COPLANAR\n";
    } else {
        std::cout << "\nVectors are NOT coplanar\n";
    }
}

void ComplexMatrixOperations() {
    ShowInputChoiceMenu();
    char inputChoice;
    std::cin >> inputChoice;
    ClearInputBuffer();

    Matrix<std::complex<double>> m1(2, 2), m2(2, 2);
    
    if (inputChoice == '1') {
        // Manual input
        std::cout << "\nEnter first matrix elements:\n";
        for (std::size_t i = 0; i < 2; ++i) {
            for (std::size_t j = 0; j < 2; ++j) {
                m1(i, j) = GetNumber<std::complex<double>>("Enter element [" + std::to_string(i) + "][" + std::to_string(j) + "]");
            }
        }
        std::cout << "\nEnter second matrix elements:\n";
        for (std::size_t i = 0; i < 2; ++i) {
            for (std::size_t j = 0; j < 2; ++j) {
                m2(i, j) = GetNumber<std::complex<double>>("Enter element [" + std::to_string(i) + "][" + std::to_string(j) + "]");
            }
        }
    } else {
        // Random values
        m1 = GetRandomMatrix<std::complex<double>>(2, 2);
        m2 = GetRandomMatrix<std::complex<double>>(2, 2);
    }

    // Display matrices
    PrintMatrix(m1, "M1");
    PrintMatrix(m2, "M2");

    // Addition
    auto sum = m1 + m2;
    PrintMatrix(sum, "M1 + M2");

    // Subtraction
    auto diff = m1 - m2;
    PrintMatrix(diff, "M1 - M2");

    // Matrix multiplication
    auto prod = m1 * m2;
    PrintMatrix(prod, "M1 * M2");

    // Scalar multiplication
    std::complex<double> scalar(2.0, 1.0);
    auto scaledM1 = m1 * scalar;
    auto scaledM2 = scalar * m2; // Commutability check
    PrintMatrix(scaledM1, "M1 * (2+i)");
    PrintMatrix(scaledM2, "(2+i) * M2");

    // Division by scalar
    auto dividedM1 = m1 / scalar;
    PrintMatrix(dividedM1, "M1 / (2+i)");

    // Trace
    std::cout << "Trace of M1: " << m1.Trace() << std::endl;
    std::cout << "Trace of M2: " << m2.Trace() << std::endl;
}

// Specialization PrintVector for complex numbers
template<>
void PrintVector(const std::vector<std::complex<double>>& vec, const std::string& vectorName) {
    std::cout << vectorName << " = (";
    for (std::size_t i = 0; i < vec.size(); ++i) {
        std::cout << vec[i];
        if (i < vec.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << ")\n";
}

int main() {
    std::cout << "Matrix Operations and Coplanarity Check\n";
    std::cout << "=======================================\n";

    char choice;
    do {
        std::cout << "\nChoose operation type:\n";
        std::cout << "1. Real number operations\n";
        std::cout << "2. Complex number operations\n";
        std::cout << "3. Vector coplanarity check (real numbers)\n";
        std::cout << "4. Vector coplanarity check (complex numbers)\n";
        std::cout << "5. Exit\n";
        std::cout << "Your choice: ";
        
        std::cin >> choice;
        ClearInputBuffer();

        switch (choice) {
            case '1': {
                RealMatrixOperations();
                break;
            }
            case '2': {
                ComplexMatrixOperations();
                break;
            }
            case '3': {
                CheckVectorCoplanarity<double>();
                break;
            }
            case '4': {
                CheckVectorCoplanarity<std::complex<double>>();
                break;
            }
            case '5':
                std::cout << "\nThank you for using the program!\n";
                break;
            default:
                std::cout << "\nInvalid choice. Please try again.\n";
        }
    } while (choice != '5');

    return 0;
}