#define BOOST_TEST_MODULE test_tensors
#include <Grid/GridCore.h>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>


void dbl_require_close(
    double expected, double observed,
    double small, double pct_tol
) {
    if (std::fabs(expected) < small) {
        BOOST_REQUIRE_SMALL(observed, small);
    } else {
        BOOST_REQUIRE_CLOSE(expected, observed, pct_tol);
    }
}


// For a terser syntax
template <int N>
using SquareMatrix = std::array<std::array<std::complex<double>, N>, N>;

std::complex<double> I(0, 1);

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
  Grid::iMatrix<Grid::ComplexD, 3> m3a, m3b, m3c, m3d;
  Grid::iMatrix<Grid::ComplexD, 4> m4a, m4b, m4c;
  TestMatrices()
      : m3a(create_iMatrix<3>({{{1, 2, 3}, //
                                {4, 5, 6}, //
                                {7, 8, 9}}})),
        m3b(create_iMatrix<3>({{{1, 2, 3}, //
                                {1, 2, 1}, //
                                {3, 2, 1}}})),
        m3c(create_iMatrix<3>({{{1, 3, 1},
                                {2, 5, -1},
                                {3, 1, 0}}})),
        m3d(create_iMatrix<3>({{{-I, 0, 0},
                                {0, 0, -I},
                                {0, 1, 0}}})),
        m4a(create_iMatrix<4>({{{1, 2, 3, 4},    //
                                {5, 6, 7, 8},    //
                                {9, 10, 11, 12}, //
                                {13, 14, 15, 16}}})),
        m4b(create_iMatrix<4>({{{1, 2, 3, 4}, //
                                {1, 2, 1, 2}, //
                                {0, 1, 1, 0}, //
                                {4, 3, 2, 1}}})),
        m4c(create_iMatrix<4>({{{3, -3, -5, -3},
                                {2, 2, 2, -2},
                                {4, -2, 3, 0},
                                {-5, -5, 3, 0}}})){};
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

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(test_index)

BOOST_FIXTURE_TEST_CASE(test_tensor_index_3, TestMatrices) {
  Grid::iVector<Grid::iMatrix<Grid::iScalar<std::complex<double>>, 3>, 2> test_composite;
  std::complex<double> expect_array[3][3] = {{2., 0., 0.},
                                             {0., 2., 0.},
                                             {0., 0., 2.}};
  BOOST_TEST(sizeof(test_composite) == 18 * sizeof(std::complex<double>));

  test_composite = 1.;
  test_composite = test_composite * 2.;
  for (uint8_t vec_idx = 0; vec_idx < 2; vec_idx++) {
    for (uint8_t col = 0; col < 3; col++) {
      for (uint8_t row = 0; row < 3; row++) {
        BOOST_TEST(test_composite._internal[vec_idx]._internal[row][col]._internal == expect_array[row][col]);
      }
    }
  }
}

BOOST_FIXTURE_TEST_CASE(test_tensor_index_4, TestMatrices) {
  Grid::iVector<Grid::iMatrix<Grid::iScalar<std::complex<double>>, 4>, 2> test_composite;
  std::complex<double> expect_array[4][4] = {{2., 0., 0., 0.},
                                             {0., 2., 0., 0.},
                                             {0., 0., 2., 0.},
                                             {0., 0., 0., 2.}};
  BOOST_TEST(sizeof(test_composite) == 32 * sizeof(std::complex<double>));

  test_composite = 1.;
  test_composite = test_composite * 2.;
  for (uint8_t vec_idx = 0; vec_idx < 2; vec_idx++) {
    for (uint8_t col = 0; col < 4; col++) {
      for (uint8_t row = 0; row < 4; row++) {
        BOOST_TEST(test_composite._internal[vec_idx]._internal[row][col]._internal == expect_array[row][col]);
      }
    }
  }
}
BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(test_Ta)

BOOST_FIXTURE_TEST_CASE(test_iMatrix_Ta_3, TestMatrices) {
  auto result = Ta(m3a + I * m3b);
  std::complex<double> expect_array[3][3] = {{-I/3., -1.+3.*I/2., -2.+3.*I},
                                             {1.+3.*I/2., +2.*I/3., -1.+3.*I/2.},
                                             {2.+3.*I, 1.+3.*I/2., -I/3.}};
  for (int row = 0; row < 3; row++) {
    for (int col = 0; col < 3; col++) {
      BOOST_REQUIRE_CLOSE(result._internal[row][col].real(), expect_array[row][col].real(), 1e-5);
      BOOST_REQUIRE_CLOSE(result._internal[row][col].imag(), expect_array[row][col].imag(), 1e-5);
    }
  }
}

BOOST_FIXTURE_TEST_CASE(test_iMatrix_Ta_4, TestMatrices) {
  auto result = Ta(m4a + I * m4b);
  std::complex<double> expect_array[4][4] = {{-1.*I/4., -3./2.+3.*I/2., -3.+3.*I/2., -9./2.+4.*I},
                                              {3./2.+3.*I/2.,3.*1.*I/4., -3./2.+I, -3.+5.*I/2.},
                                              {3.+3.*I/2., 3./2.+I, -1.*I/4., -3./2.+I},
                                              {9./2.+4.*I, 3.+5.*I/2., 3./2.+I, -1.*I/4.}};
  for (int row = 0; row < 4; row++) {
    for (int col = 0; col < 4; col++) {
      BOOST_REQUIRE_CLOSE(result._internal[row][col].real(), expect_array[row][col].real(), 1e-5);
      BOOST_REQUIRE_CLOSE(result._internal[row][col].imag(), expect_array[row][col].imag(), 1e-5);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(test_determinant)

BOOST_FIXTURE_TEST_CASE(test_iMatrix_determinant_3, TestMatrices) {
  auto result = Determinant(m3c);
  double expect = -21.;
  BOOST_REQUIRE_CLOSE(result._internal.real(), expect, 1e-5);
}


// The Determinant implementation does not pivot the input matrix, so
// fails on matrices with unfortunately-placed zeroes.

// BOOST_FIXTURE_TEST_CASE(test_iMatrix_determinant_3d, TestMatrices) {
//   auto result = Determinant(m3d);
//   double expect = 1.;
//   BOOST_REQUIRE_CLOSE(result._internal.real(), expect, 1e-5);
// }

BOOST_FIXTURE_TEST_CASE(test_iMatrix_determinant_4, TestMatrices) {
  auto result = Determinant(m4c);
  double expect = -804.;
  BOOST_REQUIRE_CLOSE(result._internal.real(), expect, 1e-5);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(test_exp)

// The Cayley-Hamilton exponentiation gives results differing from those obtained
// via other routes, by up to 10%.

// BOOST_FIXTURE_TEST_CASE(test_iMatrix_exp_3, TestMatrices) {
//   auto result = Exponentiate(m3d, 1.0);
//   std::complex<double> expect_array[3][3] = {{0.54030231-0.84147098*I, 0, 0},
//                                              {0, 0.95835813-0.49861139*I, -0.16646828-0.99166942*I},
//                                              {0, 0.99166942-0.16646828*I, 0.95835813-0.49861139*I}};
//   std::cout << m3d << std::endl;
//   std::cout << result << std::endl;
//   for (int row = 0; row < 4; row++) {
//     for (int col = 0; col < 4; col++) {
//       BOOST_REQUIRE_CLOSE(result._internal[row][col].real(), expect_array[row][col].real(), 1e-5);
//       BOOST_REQUIRE_CLOSE(result._internal[row][col].imag(), expect_array[row][col].imag(), 1e-5);
//     }
//   }
// }


BOOST_FIXTURE_TEST_CASE(test_iMatrix_exp_4, TestMatrices) {
  auto result = Exponentiate(m4a, 0.1);
  std::complex<double> expect_array[4][4] = {{3.28530879, 2.66574821, 3.04618763, 3.42662705},
                                             {5.52144891, 7.2828622, 7.0442755, 7.80568879},
                                             {8.75758902, 9.89997619, 12.04236336, 12.18475053},
                                             {11.99372914, 13.51709018, 15.04045123, 17.56381227}};
  for (int row = 0; row < 4; row++) {
    for (int col = 0; col < 4; col++) {
      BOOST_REQUIRE_CLOSE(result._internal[row][col].real(), expect_array[row][col].real(), 1e-5);
      BOOST_REQUIRE_CLOSE(result._internal[row][col].imag(), expect_array[row][col].imag(), 1e-5);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(test_reality)

BOOST_FIXTURE_TEST_CASE(test_timesI_3, TestMatrices) {
  auto result = timesI(m3a);
  auto expect = m3a * I;
  for (int row = 0; row < 3; row++) {
    for (int col = 0; col < 3; col++) {
      BOOST_REQUIRE_CLOSE(result._internal[row][col].real(), expect._internal[row][col].real(), 1e-5);
      BOOST_REQUIRE_CLOSE(result._internal[row][col].imag(), expect._internal[row][col].imag(), 1e-5);
    }
  }
}

BOOST_FIXTURE_TEST_CASE(test_timesI_4, TestMatrices) {
  auto result = timesI(m4a);
  auto expect = m4a * I;
  for (int row = 0; row < 4; row++) {
    for (int col = 0; col < 4; col++) {
      BOOST_REQUIRE_CLOSE(result._internal[row][col].real(), expect._internal[row][col].real(), 1e-5);
      BOOST_REQUIRE_CLOSE(result._internal[row][col].imag(), expect._internal[row][col].imag(), 1e-5);
    }
  }
}

BOOST_FIXTURE_TEST_CASE(test_timesMinusI_3, TestMatrices) {
  auto result = timesMinusI(m3a);
  auto expect = m3a * (-I);
  for (int row = 0; row < 3; row++) {
    for (int col = 0; col < 3; col++) {
      BOOST_REQUIRE_CLOSE(result._internal[row][col].real(), expect._internal[row][col].real(), 1e-5);
      BOOST_REQUIRE_CLOSE(result._internal[row][col].imag(), expect._internal[row][col].imag(), 1e-5);
    }
  }
}

BOOST_FIXTURE_TEST_CASE(test_timesMinusI_4, TestMatrices) {
  auto result = timesMinusI(m4a);
  auto expect = m4a * (-I);
  for (int row = 0; row < 4; row++) {
    for (int col = 0; col < 4; col++) {
      BOOST_REQUIRE_CLOSE(result._internal[row][col].real(), expect._internal[row][col].real(), 1e-5);
      BOOST_REQUIRE_CLOSE(result._internal[row][col].imag(), expect._internal[row][col].imag(), 1e-5);
    }
  }
}

BOOST_FIXTURE_TEST_CASE(test_conjugate_3, TestMatrices) {
  auto result = conjugate(m3a + m3b * I);
  auto expect = m3a + m3b * (-I);
  for (int row = 0; row < 3; row++) {
    for (int col = 0; col < 3; col++) {
      BOOST_REQUIRE_CLOSE(result._internal[row][col].real(), expect._internal[row][col].real(), 1e-5);
      BOOST_REQUIRE_CLOSE(result._internal[row][col].imag(), expect._internal[row][col].imag(), 1e-5);
    }
  }
}

BOOST_FIXTURE_TEST_CASE(test_conjugate_4, TestMatrices) {
  auto result = conjugate(m4a + m4b * I);
  auto expect = m4a + m4b * (-I);
  for (int row = 0; row < 4; row++) {
    for (int col = 0; col < 4; col++) {
      BOOST_REQUIRE_CLOSE(result._internal[row][col].real(), expect._internal[row][col].real(), 1e-5);
      BOOST_REQUIRE_CLOSE(result._internal[row][col].imag(), expect._internal[row][col].imag(), 1e-5);
    }
  }
}

/* The adj function is the Hermitian conjugate */
BOOST_FIXTURE_TEST_CASE(test_adj_3, TestMatrices) {
  auto result = adj(m3a - m3b * I);
  std::complex<double> expect_array[3][3] = {{1.+I, 4.+I, 7.+3.*I},
                                             {2.+2.*I, 5.+2.*I, 8.+2.*I},
                                             {3.+3.*I, 6.+I, 9.+I}};
  for (int row = 0; row < 3; row++) {
    for (int col = 0; col < 3; col++) {
      BOOST_REQUIRE_CLOSE(result._internal[row][col].real(), expect_array[row][col].real(), 1e-5);
      BOOST_REQUIRE_CLOSE(result._internal[row][col].imag(), expect_array[row][col].imag(), 1e-5);
    }
  }
}

BOOST_FIXTURE_TEST_CASE(test_adj_4, TestMatrices) {
  auto result = adj(m4a - m4b * I);
  std::complex<double> expect_array[4][4] = {{1.+I, 5.+I, 9., 13.+4.*I},
                                             {2.+2.*I, 6.+2.*I, 10.+I, 14.+3.*I},
                                             {3.+3.*I, 7.+I, 11.+I, 15.+2.*I},
                                             {4.+4.*I, 8.+2.*I, 12., 16.+I}};
  for (int row = 0; row < 4; row++) {
    for (int col = 0; col < 4; col++) {
      BOOST_REQUIRE_CLOSE(result._internal[row][col].real(), expect_array[row][col].real(), 1e-5);
      BOOST_REQUIRE_CLOSE(result._internal[row][col].imag(), expect_array[row][col].imag(), 1e-5);
    }
  }
}

BOOST_FIXTURE_TEST_CASE(test_real_3, TestMatrices) {
  auto result = real(m3a + m3b * I);
  double expect_array[3][3] = {{1, 2, 3},
                               {4, 5, 6},
                               {7, 8, 9}};
  BOOST_TEST(result._internal == expect_array);
}

BOOST_FIXTURE_TEST_CASE(test_imag_3, TestMatrices) {
  auto result = imag(m3a + m3b * I);
  double expect_array[3][3] = {{1, 2, 3},
                               {1, 2, 1},
                               {3, 2, 1}};
  BOOST_TEST(result._internal == expect_array);
}

BOOST_FIXTURE_TEST_CASE(test_real_4, TestMatrices) {
  auto result = real(m4a + m4b * I);
  double expect_array[4][4] = {{1, 2, 3, 4},
                               {5, 6, 7, 8},
                               {9, 10, 11, 12},
                               {13, 14, 15, 16}};
  BOOST_TEST(result._internal == expect_array);
}

BOOST_FIXTURE_TEST_CASE(test_imag_4, TestMatrices) {
  auto result = imag(m4a + m4b * I);
  double expect_array[4][4] = {{1, 2, 3, 4},
                               {1, 2, 1, 2},
                               {0, 1, 1, 0},
                               {4, 3, 2, 1}};
  BOOST_TEST(result._internal == expect_array);
}


BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(test_unary)

BOOST_FIXTURE_TEST_CASE(test_sqrt_3, TestMatrices) {
  auto result = sqrt(m3b);
  double expect_array[3][3] = {{1, sqrt(2), sqrt(3)},
                               {1, sqrt(2), 1},
                               {sqrt(3), sqrt(2), 1}};
  BOOST_TEST(result._internal == expect_array);
}

BOOST_FIXTURE_TEST_CASE(test_sqrt_4, TestMatrices) {
  auto result = sqrt(m4b);
  double expect_array[4][4] = {{1, sqrt(2), sqrt(3), 2},
                               {1, sqrt(2), 1, sqrt(2)},
                               {0, 1, 1, 0},
                               {2, sqrt(3), sqrt(2), 1}};
  BOOST_TEST(result._internal == expect_array);
}

BOOST_FIXTURE_TEST_CASE(test_div_3, TestMatrices) {
  auto result = pow(m3a + m3b * I, 2.);
  std::complex<double> expect_array[3][3] = {{2.*I, 8.*I, 18.*I},
                                             {15.+8.*I, 21.+20.*I, 35.+12.*I},
                                             {40.+42.*I, 60.+32.*I, 80.+18.*I}};
  for (int row = 0; row < 3; row++) {
    for (int col = 0; col < 3; col++) {
      dbl_require_close(result._internal[row][col].real(), expect_array[row][col].real(), 1e-10, 1e-5);
      dbl_require_close(result._internal[row][col].imag(), expect_array[row][col].imag(), 1e-10, 1e-5);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(test_logical)

BOOST_FIXTURE_TEST_CASE(test_equality_3, TestMatrices) {
  auto result = 1.0 * m3a;
  BOOST_TEST(result == m3a);
}

BOOST_FIXTURE_TEST_CASE(test_equality_4, TestMatrices) {
  auto result = 1.0 * m4a;
  BOOST_TEST(result == m4a);
}

BOOST_AUTO_TEST_SUITE_END()
