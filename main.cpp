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

template<typename T>
bool AreVectorsCoplanar(const Matrix<T>& a, const Matrix<T>& b, const Matrix<T>& c) {
    // Check if vectors are 3x1
    if (a.GetRows() != 3 || a.GetCols() != 1 ||
        b.GetRows() != 3 || b.GetCols() != 1 ||
        c.GetRows() != 3 || c.GetCols() != 1) {
        throw std::invalid_argument("Vectors must be 3x1 matrices");
    }
    
    // Create matrix from vectors
    Matrix<T> m(3, 3);
    
    // Fill matrix with vectors as columns
    for (std::size_t i = 0; i < 3; ++i) {
        m(i, 0) = a(i, 0);
        m(i, 1) = b(i, 0);
        m(i, 2) = c(i, 0);
    }
    
    // Vectors are coplanar if determinant is zero
    const double epsilon = 1e-10;
    return std::abs(std::abs(m.determinant())) < epsilon;
}

// GetNumber for complex numbers
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
        m1 = Matrix<double>(2, 2, -10.0, 10.0);
        m2 = Matrix<double>(2, 2, -10.0, 10.0);
    }

    // Display matrices
    std::cout << "\nMatrix M1:\n" << m1 << "\n";
    std::cout << "\nMatrix M2:\n" << m2 << "\n";

    // Operations
    std::cout << "\nM1 + M2:\n" << (m1 + m2) << "\n";
    std::cout << "\nM1 - M2:\n" << (m1 - m2) << "\n";
    std::cout << "\nM1 * M2:\n" << (m1 * m2) << "\n";
    
    double scalar = 2.0;
    std::cout << "\nM1 * " << scalar << ":\n" << (m1 * scalar) << "\n";
    std::cout << scalar << " * M2:\n" << (scalar * m2) << "\n";
    std::cout << "\nM1 / " << scalar << ":\n" << (m1 / scalar) << "\n";
    
    std::cout << "\nTrace of M1: " << m1.Trace() << "\n";
    std::cout << "Trace of M2: " << m2.Trace() << "\n";
}

void CheckVectorCoplanarity() {
    ShowInputChoiceMenu();
    char inputChoice;
    std::cin >> inputChoice;
    ClearInputBuffer();

    Matrix<double> a(3, 1), b(3, 1), c(3, 1);
    
    if (inputChoice == '1') {
        // Manual input
        std::cout << "\nEnter vector a:\n";
        for (std::size_t i = 0; i < 3; ++i) {
            a(i, 0) = GetNumber<double>("a" + std::to_string(i + 1) + ": ");
        }
        
        std::cout << "\nEnter vector b:\n";
        for (std::size_t i = 0; i < 3; ++i) {
            b(i, 0) = GetNumber<double>("b" + std::to_string(i + 1) + ": ");
        }
        
        std::cout << "\nEnter vector c:\n";
        for (std::size_t i = 0; i < 3; ++i) {
            c(i, 0) = GetNumber<double>("c" + std::to_string(i + 1) + ": ");
        }
    } else {
        // Random values
        a = Matrix<double>(3, 1, -10.0, 10.0);
        b = Matrix<double>(3, 1, -10.0, 10.0);
        c = Matrix<double>(3, 1, -10.0, 10.0);
    }

    // Display vectors
    std::cout << "\nVectors:\n";
    std::cout << "a:\n" << a << "\n";
    std::cout << "b:\n" << b << "\n";
    std::cout << "c:\n" << c << "\n";

    // Check coplanarity
    if (AreVectorsCoplanar(a, b, c)) {
        std::cout << "\nVectors are coplanar\n";
    } else {
        std::cout << "\nVectors are not coplanar\n";
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
                m1(i, j) = GetNumber<std::complex<double>>("Enter element [" + std::to_string(i) + "][" + std::to_string(j) + "] (real imag): ");
            }
        }
        std::cout << "\nEnter second matrix elements:\n";
        for (std::size_t i = 0; i < 2; ++i) {
            for (std::size_t j = 0; j < 2; ++j) {
                m2(i, j) = GetNumber<std::complex<double>>("Enter element [" + std::to_string(i) + "][" + std::to_string(j) + "] (real imag): ");
            }
        }
    } else {
        // Random values
        m1 = Matrix<std::complex<double>>(2, 2, std::complex<double>(-10.0, -10.0), std::complex<double>(10.0, 10.0));
        m2 = Matrix<std::complex<double>>(2, 2, std::complex<double>(-10.0, -10.0), std::complex<double>(10.0, 10.0));
    }

    // Display matrices
    std::cout << "\nMatrix M1:\n" << m1 << "\n";
    std::cout << "\nMatrix M2:\n" << m2 << "\n";

    // Operations
    std::cout << "\nM1 + M2:\n" << (m1 + m2) << "\n";
    std::cout << "\nM1 - M2:\n" << (m1 - m2) << "\n";
    std::cout << "\nM1 * M2:\n" << (m1 * m2) << "\n";
    
    std::complex<double> scalar(2.0, 1.0);
    std::cout << "\nM1 * " << scalar << ":\n" << (m1 * scalar) << "\n";
    std::cout << scalar << " * M2:\n" << (scalar * m2) << "\n";
    std::cout << "\nM1 / " << scalar << ":\n" << (m1 / scalar) << "\n";
    
    std::cout << "\nTrace of M1: " << m1.Trace() << "\n";
    std::cout << "Trace of M2: " << m2.Trace() << "\n";
}

void CheckComplexVectorCoplanarity() {
    ShowInputChoiceMenu();
    char inputChoice;
    std::cin >> inputChoice;
    ClearInputBuffer();

    Matrix<std::complex<double>> a(3, 1), b(3, 1), c(3, 1);
    
    if (inputChoice == '1') {
        // Manual input
        std::cout << "\nEnter vector a:\n";
        for (std::size_t i = 0; i < 3; ++i) {
            a(i, 0) = GetNumber<std::complex<double>>("a" + std::to_string(i + 1) + " (real imag): ");
        }
        
        std::cout << "\nEnter vector b:\n";
        for (std::size_t i = 0; i < 3; ++i) {
            b(i, 0) = GetNumber<std::complex<double>>("b" + std::to_string(i + 1) + " (real imag): ");
        }
        
        std::cout << "\nEnter vector c:\n";
        for (std::size_t i = 0; i < 3; ++i) {
            c(i, 0) = GetNumber<std::complex<double>>("c" + std::to_string(i + 1) + " (real imag): ");
        }
    } else {
        // Random values
        a = Matrix<std::complex<double>>(3, 1, std::complex<double>(-10.0, -10.0), std::complex<double>(10.0, 10.0));
        b = Matrix<std::complex<double>>(3, 1, std::complex<double>(-10.0, -10.0), std::complex<double>(10.0, 10.0));
        c = Matrix<std::complex<double>>(3, 1, std::complex<double>(-10.0, -10.0), std::complex<double>(10.0, 10.0));
    }

    // Display vectors
    std::cout << "\nVectors:\n";
    std::cout << "a:\n" << a << "\n";
    std::cout << "b:\n" << b << "\n";
    std::cout << "c:\n" << c << "\n";

    // Check coplanarity
    if (AreVectorsCoplanar(a, b, c)) {
        std::cout << "\nVectors are coplanar\n";
    } else {
        std::cout << "\nVectors are not coplanar\n";
    }
}

void TestMatrixOperations() {
    std::cout << "\n=== Testing Matrix Operations ===\n";

    // Test matrix creation and equality
    Matrix<double> m1(2, 2, 1.0);
    Matrix<double> m2(2, 2, 1.0);
    Matrix<double> m3(2, 2, 2.0);

    std::cout << "\nTesting equality operator:\n";
    std::cout << "m1 == m2 (should be true): " << (m1 == m2) << "\n";
    std::cout << "m1 == m3 (should be false): " << (m1 == m3) << "\n";

    // Test matrix arithmetic
    std::cout << "\nTesting arithmetic operations:\n";
    std::cout << "Matrix m1:\n" << m1 << "\n";
    std::cout << "Matrix m3:\n" << m3 << "\n";
    
    std::cout << "m1 + m3:\n" << (m1 + m3) << "\n";
    std::cout << "m3 - m1:\n" << (m3 - m1) << "\n";
    std::cout << "m1 * 2:\n" << (m1 * 2.0) << "\n";

    // Test matrix multiplication
    Matrix<double> m4(2, 3, 1.0);
    Matrix<double> m5(3, 2, 2.0);
    std::cout << "\nTesting matrix multiplication:\n";
    std::cout << "Matrix m4 (2x3):\n" << m4 << "\n";
    std::cout << "Matrix m5 (3x2):\n" << m5 << "\n";
    std::cout << "m4 * m5:\n" << (m4 * m5) << "\n";

    // Test trace
    Matrix<double> square(3, 3, 2.0);
    std::cout << "\nTesting trace:\n";
    std::cout << "Matrix:\n" << square << "\n";
    std::cout << "Trace: " << square.Trace() << "\n";

    // Test determinant
    std::cout << "\nTesting determinant (3x3):\n";
    std::cout << "Determinant: " << square.determinant() << "\n";

    // Test complex numbers
    std::cout << "\nTesting complex matrices:\n";
    Matrix<std::complex<double>> c1(2, 2, std::complex<double>(1.0, 1.0));
    Matrix<std::complex<double>> c2(2, 2, std::complex<double>(2.0, -1.0));
    std::cout << "Complex matrix c1:\n" << c1 << "\n";
    std::cout << "Complex matrix c2:\n" << c2 << "\n";
    std::cout << "c1 + c2:\n" << (c1 + c2) << "\n";
    std::cout << "c1 * c2:\n" << (c1 * c2) << "\n";
}

int main() {
    std::cout << "Matrix Operations and Coplanarity Check\n";
    std::cout << "=======================================\n";

    char choice;
    do {
        std::cout << "\nChoose operation type:\n";
        std::cout << "1. Real matrix operations\n";
        std::cout << "2. Complex matrix operations\n";
        std::cout << "3. Real vector coplanarity check\n";
        std::cout << "4. Complex vector coplanarity check\n";
        std::cout << "5. Run tests\n";
        std::cout << "6. Exit\n";
        std::cout << "Your choice: ";
        
        std::cin >> choice;
        ClearInputBuffer();

        switch (choice) {
            case '1':
                RealMatrixOperations();
                break;
            case '2':
                ComplexMatrixOperations();
                break;
            case '3':
                CheckVectorCoplanarity();
                break;
            case '4':
                CheckComplexVectorCoplanarity();
                break;
            case '5':
                TestMatrixOperations();
                break;
            case '6':
                std::cout << "\nThank you for using the program!\n";
                break;
            default:
                std::cout << "\nInvalid choice. Please try again.\n";
        }
    } while (choice != '6');

    return 0;
}