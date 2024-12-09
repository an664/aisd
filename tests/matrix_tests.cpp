#include <gtest/gtest.h>
#include <cmath>
#include <complex>
#include "../matrix.hpp"

// Basic constructor tests
TEST(MatrixTest, ConstructorTest) {
    Matrix<double> m(3, 3, 1.0);
    EXPECT_EQ(m(0, 0), 1.0);
    EXPECT_EQ(m(1, 1), 1.0);
    EXPECT_EQ(m(2, 2), 1.0);
}

// Default constructor test (with value)
TEST(MatrixTest, DefaultConstructorTest) {
    Matrix<double> m1(2, 3, 1.5);
    Matrix<double> m2(2, 3, 1.5);
    Matrix<double> m3(2, 3, 2.0);
    
    EXPECT_EQ(m1.GetRows(), 2);
    EXPECT_EQ(m1.GetCols(), 3);
    EXPECT_EQ(m1, m2);      
    EXPECT_NE(m1, m3);      
}

// Random constructor test
TEST(MatrixTest, RandomConstructorTest) {
    Matrix<double> m1(3, 2, -1.0, 1.0);
    Matrix<double> m2(m1);  
    
    EXPECT_EQ(m1.GetRows(), 3);
    EXPECT_EQ(m1.GetCols(), 2);
    EXPECT_EQ(m1, m2);      
    
    m2(0, 0) += 1.0;       
    EXPECT_NE(m1, m2);     
}

// Copy constructor test
TEST(MatrixTest, CopyConstructorTest) {
    Matrix<double> m1(2, 2, 3.0);
    Matrix<double> m2(m1);
    Matrix<double> m3(2, 2, 3.0);
    
    EXPECT_EQ(m1, m2);      
    EXPECT_EQ(m1, m3);      
    
    m2(0, 0) = 5.0;
    EXPECT_NE(m1, m2);     
}

// Complex number constructor test
TEST(MatrixTest, ComplexConstructorTest) {
    using Complex = std::complex<double>;
    Matrix<Complex> m1(2, 2, Complex(1.0, 2.0));
    Matrix<Complex> m2(2, 2, Complex(1.0, 2.0));
    Matrix<Complex> m3(2, 2, Complex(1.0, 2.1));
    
    EXPECT_EQ(m1, m2);      
    EXPECT_NE(m1, m3);      
}

// Random complex constructor test
TEST(MatrixTest, RandomComplexConstructorTest) {
    using Complex = std::complex<double>;
    Complex min(-1.0, -1.0), max(1.0, 1.0);
    Matrix<Complex> m1(2, 2, min, max);
    Matrix<Complex> m2(m1);
    
    EXPECT_EQ(m1, m2);      
    
    m2(0, 0) = Complex(2.0, 2.0);
    EXPECT_NE(m1, m2);      
}

// Invalid size constructor test
TEST(MatrixTest, InvalidConstructorTest) {
    EXPECT_THROW(Matrix<double>(0, 1), std::invalid_argument);
    EXPECT_THROW(Matrix<double>(1, 0), std::invalid_argument);
}

// Addition test
TEST(MatrixTest, AdditionTest) {
    Matrix<double> m1(2, 2, 1.0);
    Matrix<double> m2(2, 2, 2.0);
    Matrix<double> result = m1 + m2;
    EXPECT_EQ(result(0, 0), 3.0);
    EXPECT_EQ(result(0, 1), 3.0);
    EXPECT_EQ(result(1, 0), 3.0);
    EXPECT_EQ(result(1, 1), 3.0);
}

// Subtraction test
TEST(MatrixTest, SubtractionTest) {
    Matrix<double> m1(2, 2, 3.0);
    Matrix<double> m2(2, 2, 1.0);
    Matrix<double> result = m1 - m2;
    EXPECT_EQ(result(0, 0), 2.0);
    EXPECT_EQ(result(0, 1), 2.0);
    EXPECT_EQ(result(1, 0), 2.0);
    EXPECT_EQ(result(1, 1), 2.0);
}

// Multiplication test
TEST(MatrixTest, MultiplicationTest) {
    Matrix<double> m1(2, 2, 2.0);
    Matrix<double> m2(2, 2, 3.0);
    Matrix<double> result = m1 * m2;
    EXPECT_EQ(result(0, 0), 12.0);
    EXPECT_EQ(result(0, 1), 12.0);
    EXPECT_EQ(result(1, 0), 12.0);
    EXPECT_EQ(result(1, 1), 12.0);
}

// Scalar multiplication test
TEST(MatrixTest, ScalarMultiplicationTest) {
    Matrix<double> m(2, 2, 2.0);
    Matrix<double> result1 = m * 3.0;
    Matrix<double> result2 = 3.0 * m;
    
    EXPECT_EQ(result1(0, 0), 6.0);
    EXPECT_EQ(result1(1, 1), 6.0);
    EXPECT_EQ(result2(0, 0), 6.0);
    EXPECT_EQ(result2(1, 1), 6.0);
}

// Scalar division test
TEST(MatrixTest, ScalarDivisionTest) {
    Matrix<double> m(2, 2, 6.0);
    Matrix<double> result = m / 2.0;
    EXPECT_EQ(result(0, 0), 3.0);
    EXPECT_EQ(result(0, 1), 3.0);
    EXPECT_EQ(result(1, 0), 3.0);
    EXPECT_EQ(result(1, 1), 3.0);
}

// Trace test
TEST(MatrixTest, TraceTest) {
    Matrix<double> m(3, 3, 2.0);
    EXPECT_EQ(m.Trace(), 6.0);
}

// Complex number tests
TEST(MatrixTest, ComplexTest) {
    using Complex = std::complex<double>;
    Matrix<Complex> m(2, 2, Complex(1.0, 1.0));
    
    // Test addition
    Matrix<Complex> result = m + m;
    EXPECT_EQ(result(0, 0), Complex(2.0, 2.0));
    
    // Test multiplication
    result = m * m;
    EXPECT_EQ(result(0, 0), Complex(0.0, 4.0));
    
    // Test trace
    EXPECT_EQ(m.Trace(), Complex(2.0, 2.0));
}

// Coplanarity tests
class CoplanarityTest : public ::testing::Test {
protected:
    bool TestCoplanarity(const std::vector<double>& a, 
                        const std::vector<double>& b, 
                        const std::vector<double>& c) {
        Matrix<double> m(3, 3);
        for (std::size_t i = 0; i < 3; ++i) {
            m(i, 0) = a[i];
            m(i, 1) = b[i];
            m(i, 2) = c[i];
        }
        const double epsilon = 1e-10;
        return std::abs(m.determinant()) < epsilon;
    }
};

TEST_F(CoplanarityTest, XYPlaneVectors) {
    std::vector<double> a = {1.0, 0.0, 0.0};
    std::vector<double> b = {0.0, 1.0, 0.0};
    std::vector<double> c = {1.0, 1.0, 0.0};
    
    EXPECT_TRUE(TestCoplanarity(a, b, c));
}

TEST_F(CoplanarityTest, NonCoplanarVectors) {
    std::vector<double> a = {1.0, 0.0, 0.0};
    std::vector<double> b = {0.0, 1.0, 0.0};
    std::vector<double> c = {0.0, 0.0, 1.0};
    
    EXPECT_FALSE(TestCoplanarity(a, b, c));
}

TEST_F(CoplanarityTest, ZeroVectors) {
    std::vector<double> zero = {0.0, 0.0, 0.0};
    
    EXPECT_TRUE(TestCoplanarity(zero, zero, zero));
}

TEST_F(CoplanarityTest, ParallelVectors) {
    std::vector<double> a = {1.0, 2.0, 3.0};
    std::vector<double> b = {2.0, 4.0, 6.0};
    std::vector<double> c = {3.0, 6.0, 9.0};
    
    EXPECT_TRUE(TestCoplanarity(a, b, c));
}

TEST_F(CoplanarityTest, LargeValues) {
    std::vector<double> a = {1e8, 0.0, 0.0};
    std::vector<double> b = {0.0, 1e8, 0.0};
    std::vector<double> c = {1e8, 1e8, 0.0};
    
    EXPECT_TRUE(TestCoplanarity(a, b, c));
}

// Determinant test
TEST(MatrixTest, DeterminantTest) {
    Matrix<double> m(3, 3, 0.0);
    
    // Identity matrix
    m(0, 0) = 1.0; m(0, 1) = 0.0; m(0, 2) = 0.0;
    m(1, 0) = 0.0; m(1, 1) = 1.0; m(1, 2) = 0.0;
    m(2, 0) = 0.0; m(2, 1) = 0.0; m(2, 2) = 1.0;
    
    EXPECT_DOUBLE_EQ(m.determinant(), 1.0);
}

// Equality operator test
TEST(MatrixTest, EqualityTest) {
    Matrix<double> m1(2, 2, 1.0);
    Matrix<double> m2(2, 2, 1.0);
    Matrix<double> m3(2, 2, 1.1);
    Matrix<double> m4(3, 3, 1.0);

    EXPECT_TRUE(m1 == m2);
    EXPECT_FALSE(m1 == m3);
    EXPECT_FALSE(m1 == m4);
}

// Inequality operator test
TEST(MatrixTest, InequalityTest) {
    Matrix<double> m1(2, 2, 1.0);
    Matrix<double> m2(2, 2, 1.0);
    Matrix<double> m3(2, 2, 1.1);

    EXPECT_FALSE(m1 != m2);
    EXPECT_TRUE(m1 != m3);
}

// Floating point comparison precision test
TEST(MatrixTest, FloatingPointComparisonTest) {
    Matrix<double> m1(2, 2, 1.0);
    Matrix<double> m2(2, 2, 1.0 + 1e-11);  
    Matrix<double> m3(2, 2, 1.0 + 1e-9);   
    
    EXPECT_EQ(m1, m2);      
    EXPECT_NE(m1, m3);      
}

// Complex number equality test
TEST(MatrixTest, ComplexEqualityTest) {
    using Complex = std::complex<double>;
    Matrix<Complex> m1(2, 2, Complex(1.0, 1.0));
    Matrix<Complex> m2(2, 2, Complex(1.0, 1.0));
    Matrix<Complex> m3(2, 2, Complex(1.0, 1.1));

    EXPECT_TRUE(m1 == m2);
    EXPECT_FALSE(m1 == m3);
}

// Stream output operator test
TEST(MatrixTest, StreamOutputTest) {
    Matrix<double> m(2, 2, 1.5);
    std::stringstream ss;
    ss << m;
    std::string expected = "1.5 1.5\n1.5 1.5";
    EXPECT_EQ(ss.str(), expected);
}

// Complex stream output test
TEST(MatrixTest, ComplexStreamOutputTest) {
    using Complex = std::complex<double>;
    Matrix<Complex> m(2, 2, Complex(1.0, 2.0));
    std::stringstream ss;
    ss << m;
    EXPECT_FALSE(ss.str().empty());
}

// Different size matrices equality test
TEST(MatrixTest, DifferentSizeEqualityTest) {
    Matrix<double> m1(2, 2, 1.0);
    Matrix<double> m2(2, 3, 1.0);
    Matrix<double> m3(3, 2, 1.0);

    EXPECT_FALSE(m1 == m2);
    EXPECT_FALSE(m1 == m3);
    EXPECT_FALSE(m2 == m3);
}

// Edge cases for equality test
TEST(MatrixTest, EdgeCasesEqualityTest) {
    Matrix<double> m1(1, 1, 0.0);
    Matrix<double> m2(1, 1, -0.0);
    Matrix<double> m3(1, 1, 1e-11);

    EXPECT_TRUE(m1 == m2);
    EXPECT_TRUE(m1 == m3);
}

// Test scalar multiplication commutativity
TEST(MatrixTest, ScalarMultiplicationCommutativityTest) {
    Matrix<double> m1(2, 2, 1.0);
    double scalar = 2.5;
    
    Matrix<double> result1 = m1 * scalar;
    Matrix<double> result2 = scalar * m1;
    
    EXPECT_EQ(result1, result2);
    
    // Проверяем для комплексных чисел
    using Complex = std::complex<double>;
    Matrix<Complex> m2(2, 2, Complex(1.0, 1.0));
    Complex complex_scalar(2.0, 1.0);
    
    Matrix<Complex> result3 = m2 * complex_scalar;
    Matrix<Complex> result4 = complex_scalar * m2;
    
    EXPECT_EQ(result3, result4);
} 