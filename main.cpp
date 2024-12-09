#include <iostream>
#include <cmath>
#include "matrix.hpp"

// Clear input buffer function
void ClearInputBuffer() {
    std::cin.clear();
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

// Get number from user with validation
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

// Print matrix function
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

// Check if three vectors are coplanar
bool AreVectorsCoplanar(const Matrix<double>& a, const Matrix<double>& b, const Matrix<double>& c) {
    // Create matrix 3x3 from vectors
    Matrix<double> vectors(3, 3);
    
    // Fill matrix with vectors as columns
    for (std::size_t i = 0; i < 3; ++i) {
        vectors(i, 0) = a(i, 0);
        vectors(i, 1) = b(i, 0);
        vectors(i, 2) = c(i, 0);
    }
    
    // Vectors are coplanar if determinant is zero
    return std::abs(vectors.determinant()) < 1e-10;
}

// Test matrix operations
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
}

int main() {
    try {
        // Run tests first
        TestMatrixOperations();
        
        std::cout << "\n=== Coplanar Vectors Check ===\n";
        // Create three 3x1 vectors
        Matrix<double> a(3, 1);
        Matrix<double> b(3, 1);
        Matrix<double> c(3, 1);

        // Input vectors
        std::cout << "Enter vector a (3 components):\n";
        for (std::size_t i = 0; i < 3; ++i) {
            a(i, 0) = GetNumber<double>("a" + std::to_string(i + 1) + ": ");
        }

        std::cout << "Enter vector b (3 components):\n";
        for (std::size_t i = 0; i < 3; ++i) {
            b(i, 0) = GetNumber<double>("b" + std::to_string(i + 1) + ": ");
        }

        std::cout << "Enter vector c (3 components):\n";
        for (std::size_t i = 0; i < 3; ++i) {
            c(i, 0) = GetNumber<double>("c" + std::to_string(i + 1) + ": ");
        }

        // Print vectors
        PrintMatrix(a, "a");
        PrintMatrix(b, "b");
        PrintMatrix(c, "c");

        // Check if vectors are coplanar
        if (AreVectorsCoplanar(a, b, c)) {
            std::cout << "\nVectors are coplanar.\n";
        } else {
            std::cout << "\nVectors are not coplanar.\n";
        }

    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}