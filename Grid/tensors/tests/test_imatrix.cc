#define BOOST_TEST_MODULE test_tensors
#include <Grid/GridCore.h>
#include <boost/test/included/unit_test.hpp>

// For a terser syntax
template <int N>
using SquareMatrix = std::array<std::array<std::complex<double>, N>, N>;

template <int N>
Grid::iMatrix<Grid::ComplexD, N> create_iMatrix(SquareMatrix<N> data) {
  Grid::iMatrix<Grid::ComplexD, N> M;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      M._internal[i][j] = data[i][j];
    }
  }
  return M;
}

template <int N>
Grid::iVector<Grid::ComplexD, N> create_iVector(std::array<std::complex<double>, N> data) {
  Grid::iVector<Grid::ComplexD, N> V;
  for (int i = 0; i < N; i++) {
    V._internal[i] = data[i];
  }
  return V;
}


struct TestMatrices {
  Grid::iMatrix<Grid::ComplexD, 3> m3a, m3b;
  Grid::iMatrix<Grid::ComplexD, 4> m4a, m4b;
  TestMatrices()
      : m3a(create_iMatrix<3>({{{1, 2, 3}, //
                                {4, 5, 6}, //
                                {7, 8, 9}}})),
        m3b(create_iMatrix<3>({{{1, 2, 3}, //
                                {1, 2, 1}, //
                                {3, 2, 1}}})),
        m4a(create_iMatrix<4>({{{1, 2, 3, 4},    //
                                {5, 6, 7, 8},    //
                                {9, 10, 11, 12}, //
                                {13, 14, 15, 16}}})),
        m4b(create_iMatrix<4>({{{1, 2, 3, 4}, //
                                {1, 2, 1, 2}, //
                                {0, 1, 1, 0}, //
                                {4, 3, 2, 1}}})){};
};

struct TestVectors {
  Grid::iVector<Grid::ComplexD, 3> v3a, v3b;
  Grid::iVector<Grid::ComplexD, 4> v4a, v4b;
  TestVectors()
    : v3a(create_iVector<3>({{1, 2, 3}})),
      v3b(create_iVector<3>({{3, 2, 4}})),
      v4a(create_iVector<4>({{1, 2, 3, 4}})),
      v4b(create_iVector<4>({{3, 0, 1, 4}})){};
};


BOOST_AUTO_TEST_SUITE(test_identity)

BOOST_AUTO_TEST_CASE(test_iMatrix_id_3) {
  Grid::iMatrix<Grid::ComplexD, 3> id = static_cast<Grid::ComplexD>(1.);
  Grid::iMatrix<Grid::ComplexD, 3> id_verify;
  for (int i = 0; i < 3; i++) {
    id_verify._internal[i][i] = 1.;
  }
  BOOST_TEST(id == id_verify);
}

BOOST_AUTO_TEST_CASE(test_iMatrix_id_4) {
  Grid::iMatrix<Grid::ComplexD, 4> id = static_cast<Grid::ComplexD>(1.);
  Grid::iMatrix<Grid::ComplexD, 4> id_verify;
  for (int i = 0; i < 4; i++) {
    id_verify._internal[i][i] = 1.;
  }
  BOOST_TEST(id == id_verify);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(test_arithmetic)

BOOST_FIXTURE_TEST_CASE(test_iMatrix_add_3, TestMatrices) {
  auto result = m3a + m3b;

  std::complex<double> expect_array[3][3] = {{2, 4, 6}, //
                                             {5, 7, 7}, //
                                             {10, 10, 10}};
  BOOST_TEST(result._internal == expect_array);
}

BOOST_FIXTURE_TEST_CASE(test_iMatrix_add_4, TestMatrices) {
  auto result = m4a + m4b;

  std::complex<double> expect_array[4][4] = {{2, 4, 6, 8},    //
                                             {6, 8, 8, 10},   //
                                             {9, 11, 12, 12}, //
                                             {17, 17, 17, 17}};
  BOOST_TEST(result._internal == expect_array);
}

BOOST_FIXTURE_TEST_CASE(test_iMatrix_sub_3, TestMatrices) {
  auto result = m3a - m3b;

  std::complex<double> expect_array[3][3] = {{0, 0, 0}, //
                                             {3, 3, 5},
                                             {4, 6, 8}};
  BOOST_TEST(result._internal == expect_array);
}

BOOST_FIXTURE_TEST_CASE(test_iMatrix_sub_4, TestMatrices) {
  auto result = m4a - m4b;

  std::complex<double> expect_array[4][4] = {{0, 0, 0, 0},   //
                                             {4, 4, 6, 6},   //
                                             {9, 9, 10, 12}, //
                                             {9, 11, 13, 15}};
  BOOST_TEST(result._internal == expect_array);
}

BOOST_FIXTURE_TEST_CASE(test_iMatrix_mul_3, TestMatrices) {
  auto result = m3a * m3b;

  std::complex<double> expect_array[3][3] = {{12, 12, 8},  //
                                             {27, 30, 23}, //
                                             {42, 48, 38}};
  BOOST_TEST(result._internal == expect_array);
}

BOOST_FIXTURE_TEST_CASE(test_iMatrix_mul_4, TestMatrices) {
  auto result = m4a * m4b;

  std::complex<double> expect_array[4][4] = {{19, 21, 16, 12}, //
                                             {43, 53, 44, 40}, //
                                             {67, 85, 72, 68}, //
                                             {91, 117, 100, 96}};
  BOOST_TEST(result._internal == expect_array);
}

BOOST_FIXTURE_TEST_CASE(test_iMatrix_minus_3, TestMatrices) {
  auto result = -m3a;

  std::complex<double> expect_array[3][3] = {{-1,-2,-3},  //
                                             {-4,-5,-6}, //
                                             {-7,-8,-9}};
  BOOST_TEST(result._internal == expect_array);
}

BOOST_FIXTURE_TEST_CASE(test_iMatrix_minus_4, TestMatrices) {
  auto result = -m4a;

  std::complex<double> expect_array[4][4] = {{- 1,- 2,- 3,- 4}, //
                                             {- 5,- 6,- 7,- 8}, //
                                             {- 9,-10,-11,-12}, //
                                             {-13,-14,-15,-16}};
  BOOST_TEST(result._internal == expect_array);
}

BOOST_FIXTURE_TEST_CASE(test_iMatrix_scalar_mul_3, TestMatrices) {
  auto result = m3a * 2.0;
  std::complex<double> expect_array[3][3] = {{2, 4, 6},
                                             {8, 10, 12},
                                             {14, 16, 18}};
  BOOST_TEST(result._internal == expect_array);
}

BOOST_FIXTURE_TEST_CASE(test_iMatrix_scalar_mul_4, TestMatrices) {
  auto result = m4a * 2.0;
  std::complex<double> expect_array[4][4] = {{2, 4, 6, 8},
                                             {10, 12, 14, 16},
                                             {18, 20, 22, 24},
                                             {26, 28, 30, 32}};
  BOOST_TEST(result._internal == expect_array);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(test_innerproduct)

BOOST_FIXTURE_TEST_CASE(test_iMatrix_innerproduct_3, TestMatrices) {
  std::complex<double> result = innerProductD(m3a, m3b);
  std::complex<double> expect = 80;
  BOOST_TEST(result == expect);
}

BOOST_FIXTURE_TEST_CASE(test_iMatrix_innerproduct_4, TestMatrices) {
  std::complex<double> result = innerProductD(m4a, m4b);
  std::complex<double> expect = 231;
  BOOST_TEST(result == expect);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(test_outerproduct)

BOOST_FIXTURE_TEST_CASE(test_iVector_outerproduct_3, TestVectors) {
  auto result = outerProduct(v3a, v3b);
  std::complex<double> expect_array[3][3] = {{3, 2, 4},
                                             {6, 4, 8},
                                             {9, 6, 12}};
  BOOST_TEST(result._internal == expect_array);
}

BOOST_FIXTURE_TEST_CASE(test_iVector_outerproduct_4, TestVectors) {
  auto result = outerProduct(v4a, v4b);
  std::complex<double> expect_array[4][4] = {{3, 0, 1, 4},
                                             {6, 0, 2, 8},
                                             {9, 0, 3, 12},
                                             {12, 0, 4, 16}};
  BOOST_TEST(result._internal == expect_array);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(test_transpose)

BOOST_FIXTURE_TEST_CASE(test_iMatrix_transpose_3, TestMatrices) {
  auto result = transpose(m3a);
  std::complex<double> expect_array[3][3] = {{1, 4, 7},
                                             {2, 5, 8},
                                             {3, 6, 9}};
  BOOST_TEST(result._internal == expect_array);
}

BOOST_FIXTURE_TEST_CASE(test_iMatrix_transpose_4, TestMatrices) {
  auto result = transpose(m4a);
  std::complex<double> expect_array[4][4] = {{1, 5, 9, 13},
                                             {2, 6, 10, 14},
                                             {3, 7, 11, 15},
                                             {4, 8, 12, 16}};
  BOOST_TEST(result._internal == expect_array);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(test_trace)

BOOST_FIXTURE_TEST_CASE(test_iMatrix_trace_3, TestMatrices) {
  std::complex<double> result = trace(m3a);
  std::complex<double> expect = 15;
  BOOST_TEST(result == expect);
}

BOOST_FIXTURE_TEST_CASE(test_iMatrix_trace_4, TestMatrices) {
  std::complex<double> result = trace(m4a);
  std::complex<double> expect = 34;
  BOOST_TEST(result == expect);
}

// TODO: Test Tensor_index.h functions based on some understanding of what they do

// TODO: Test Ta function based on some understanding of what it does.
// Appears to be 0.5(M - M*) - 1/2Nc Tr (M - M*), but result is not pure imaginary off diagonal.

BOOST_AUTO_TEST_SUITE_END()
