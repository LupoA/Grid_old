#define BOOST_TEST_MODULE test_tensors
#include <boost/test/included/unit_test.hpp>
#include <Grid/GridCore.h>


template <int N>
Grid::iMatrix<Grid::ComplexD, N> create_iMatrix(std::complex<double> data[N][N]) {
  Grid::iMatrix<Grid::ComplexD, N> M;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      M._internal[i][j] = data[i][j];
    }
  }
  return M;
}


template <int N>
struct arb_mat {
  arb_mat() {
    M = new Grid::iMatrix<Grid::ComplexD, N>;
  }
  ~arb_mat() {
    delete M;
  }
  Grid::iMatrix<Grid::ComplexD, N>* M;
};

struct arb_mat_3a : arb_mat<3> {
  arb_mat_3a() : arb_mat<3>() {
    std::complex<double> data[3][3] = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    *M = create_iMatrix<3>(data);
  }
};

struct arb_mat_3b : arb_mat<3> {
  arb_mat_3b() : arb_mat<3>() {
    std::complex<double> data[3][3] = {{1, 2, 3}, {1, 2, 1}, {3, 2, 1}};
    *M = create_iMatrix<3>(data);
  }
};

struct arb_mat_4a : arb_mat<4> {
  arb_mat_4a() : arb_mat<4>() {
    std::complex<double> data[4][4] = {{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}};
    *M = create_iMatrix<4>(data);
  }
};

struct arb_mat_4b : arb_mat<4> {
  arb_mat_4b() : arb_mat<4>() {
    std::complex<double> data[4][4] = {{1, 2, 3, 4}, {1, 2, 1, 2}, {0, 1, 1, 0}, {4, 3, 2, 1}};
    *M = create_iMatrix<4>(data);
  }
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

BOOST_AUTO_TEST_CASE(test_iMatrix_add_3) {
  arb_mat_3a M1;
  arb_mat_3b M2;
  Grid::iMatrix<Grid::ComplexD, 3> result = *M1.M + *M2.M;

  std::complex<double> expect_array[3][3] = {{2, 4, 6}, {5, 7, 7}, {10, 10, 10}};
  BOOST_TEST(result._internal == expect_array);
}

BOOST_AUTO_TEST_CASE(test_iMatrix_add_4) {
  arb_mat_4a M1;
  arb_mat_4b M2;
  Grid::iMatrix<Grid::ComplexD, 4> result = *M1.M + *M2.M;

  std::complex<double> expect_array[4][4] = {{2, 4, 6, 8}, {6, 8, 8, 10}, {9, 11, 12, 12}, {17, 17, 17, 17}};
  BOOST_TEST(result._internal == expect_array);
}

BOOST_AUTO_TEST_CASE(test_iMatrix_sub_3) {
  arb_mat_3a M1;
  arb_mat_3b M2;
  Grid::iMatrix<Grid::ComplexD, 3> result = *M1.M - *M2.M;

  std::complex<double> expect_array[3][3] = {{0, 0, 0}, {3, 3, 5}, {4, 6, 8}};
  BOOST_TEST(result._internal == expect_array);
}

BOOST_AUTO_TEST_CASE(test_iMatrix_sub_4) {
  arb_mat_4a M1;
  arb_mat_4b M2;
  Grid::iMatrix<Grid::ComplexD, 4> result = *M1.M - *M2.M;

  std::complex<double> expect_array[4][4] = {{0, 0, 0, 0}, {4, 4, 6, 6}, {9, 9, 10, 12}, {9, 11, 13, 15}};
  BOOST_TEST(result._internal == expect_array);
}

BOOST_AUTO_TEST_CASE(test_iMatrix_mul_3) {
  arb_mat_3a M1;
  arb_mat_3b M2;
  Grid::iMatrix<Grid::ComplexD, 3> result = *M1.M * *M2.M;

  std::complex<double> expect_array[3][3] = {{12, 12, 8}, {27, 30, 23}, {42, 48, 38}};
  BOOST_TEST(result._internal == expect_array);
}

BOOST_AUTO_TEST_CASE(test_iMatrix_mul_4) {
  arb_mat_4a M1;
  arb_mat_4b M2;
  Grid::iMatrix<Grid::ComplexD, 4> result = *M1.M * *M2.M;

  std::complex<double> expect_array[4][4] = {{19, 21, 16, 12}, {43, 53, 44, 40}, {67, 85, 72, 68}, {91, 117, 100, 96}};
  BOOST_TEST(result._internal == expect_array);
}




BOOST_AUTO_TEST_SUITE_END()
