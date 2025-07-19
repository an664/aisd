# Matrix Operations Project

This repository contains a small C++ project showcasing a templated `Matrix` class together with an interactive demo application and unit tests based on Google Test.

## Features
- Generic `Matrix<T>` implementation supporting real and `std::complex<double>` values
- Construction with constants or random values
- Basic arithmetic: addition, subtraction, multiplication, division
- Matrix multiplication and scalar operations
- Trace and 3×3 determinant functions
- Coplanarity check for three 3D vectors
- Unit tests covering the above functionality

## Building
Ensure you have a C++17 compiler. To build the demo application run:

```bash
g++ -std=c++17 main.cpp -o matrix_demo
```

To compile and run the unit tests you need Google Test installed (Debian/Ubuntu package `libgtest-dev`):

```bash
g++ -std=c++17 tests/matrix_tests.cpp -lgtest -lgtest_main -lpthread -o matrix_tests
./matrix_tests
```

## Usage
Run `./matrix_demo` to start the program. A menu allows you to perform real or complex matrix operations, check coplanarity of vectors or execute a short demo of the API.

## Repository layout
- `matrix.hpp` &ndash; templated matrix implementation
- `main.cpp` &ndash; interactive console program and helper routines
- `tests/` &ndash; Google Test suite

This project is intended as a simple introduction to template programming and basic linear algebra with C++.
