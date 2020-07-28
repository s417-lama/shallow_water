#include <iostream>
#include <fstream>
#include <memory>

extern "C" {
  void dgtsv_(const int& N, const int& NRHS,
              double* DL, double* D, double* DU,
              double* B, const int& LDB, int* INFO);
};

static void dgtsv(const int& N, double* DL, double* D, double* DU, double* B) {
  int info;
  dgtsv_(N, 1, DL, D, DU, B, N, &info);
}

constexpr double g = 9.81;

void solve_2d(double* h, const double* h_1, const double* h_2, const double* b,
              const double dx, const double dt, const int N) {
  std::unique_ptr<double[]> d(new double[N]);
  std::unique_ptr<double[]> D(new double[N]);
  std::unique_ptr<double[]> DL(new double[N]);
  std::unique_ptr<double[]> DU(new double[N]);

  for (int i = 0; i < N; i++) {
    d[i] = h_1[i] - b[i];
  }

  auto e = [&](const int& i) {
    if (i == 0) {
      return 1.0 + g * dt * dt * (d[0] + d[1]) / (2.0 * dx * dx);
    } else if (i == N - 1) {
      return 1.0 + g * dt * dt * (d[N-2] + d[N-1]) / (2.0 * dx * dx);
    } else {
      return 1.0 + g * dt * dt * (d[i-1] + 2.0 * d[i] + d[i+1]) / (2.0 * dx * dx);
    }
  };
  auto f = [&](const int& i) {
    if (i == N - 1) {
      return 0.0;
    } else {
      return -g * dt * dt * (d[i] + d[i+1]) / (2.0 * dx * dx);
    }
  };

  for (int i = 0; i < N; i++) {
    D[i] = e(i);
    DL[i] = f(i);
    DU[i] = f(i);
    h[i] = 2.0 * h_1[i] - h_2[i];
  }

  // Ax = b
  //   - D, DL, DU: tridiagonal matrix A
  //   - h: input as b and output as x
  dgtsv(N, DL.get(), D.get(), DU.get(), h);
}

void shallow_water_2d(int N) {
  constexpr double dt = 0.1;
  constexpr double dx = 0.1;
  constexpr int ITER = 1000;
  constexpr int output_intvl = 10;

  std::unique_ptr<double[]> b(new double[N]);
  std::unique_ptr<double[]> h(new double[N]);
  std::unique_ptr<double[]> h_1(new double[N]);
  std::unique_ptr<double[]> h_2(new double[N]);

  // init
  for (int i = 0; i < N; i++) {
    b[i] = 0.0;
    /* h[i] = 1.0; */
    h[i] = (i == N / 2) ? 0.5 : 1.0;
    h_1[i] = h[i];
    h_2[i] = h[i];
  }

  // main loop
  for (int it = 0; it < ITER; it++) {
    solve_2d(h.get(), h_1.get(), h_2.get(), b.get(), dx, dt, N);

    // output
    if (it % output_intvl == 0) {
      char filename[20];
      std::snprintf(filename, 20, "out2d/h_%04d.txt", it / output_intvl);
      std::ofstream ofs(filename);
      for(int i = 0; i < N; i++) {
        ofs << dx * i << "," << h[i] << std::endl;
      }
      std::cout << "Output to: " << filename << std::endl;
    }

    // h <- h_2
    // h_1 <- h
    // h_2 <- h_1
    h_1.swap(h_2);
    h.swap(h_1);
  }
}

void shallow_water_3d(int NX, int NY) {
  constexpr double dt = 0.01;
  constexpr double dx = 0.1;
  constexpr double dy = 0.1;
  /* constexpr double dx = 1.0; */
  /* constexpr double dy = 1.0; */
  constexpr int ITER = 1000;
  constexpr int output_intvl = 10;

  std::unique_ptr<double[]> b(new double[NX*NY]);
  std::unique_ptr<double[]> h(new double[NX*NY]);
  std::unique_ptr<double[]> h_1(new double[NX*NY]);
  std::unique_ptr<double[]> h_2(new double[NX*NY]);

  // init
  for (int y = 0; y < NY; y++) {
    for (int x = 0; x < NX; x++) {
      int i = y * NX + x;
      b[i] = 0.0;
      /* h[i] = 1.0; */
      /* h[i] = (NX / 2 - 5 < x && x < NX / 2 + 5 && NY / 2 - 5 < y && y < NY / 2 + 5) ? 0.5 : 1.0; */
      h[i] = (NX / 3 < x && x < NX / 3 + 20 && NY / 3 < y && y < NY / 3 + 20) ? 0.5 : 3.0;
      /* h[i] = (x < 10 && y < 10) ? 0.5 : 1.0; */
      h_1[i] = h[i];
      h_2[i] = h[i];
    }
  }

  // main loop
  for (int it = 0; it < ITER; it++) {
    // solve x
    for (int y = 0; y < NY; y++) {
      double* h_ = &h.get()[y * NX];
      double* h_1_ = &h_1.get()[y * NX];
      double* h_2_ = &h_2.get()[y * NX];
      double* b_ = &b.get()[y * NX];
      solve_2d(h_, h_1_, h_2_, b_, dx, dt, NX);
    }

    // solve y
    std::unique_ptr<double[]> b_(new double[NY]);
    std::unique_ptr<double[]> h_(new double[NY]);
    std::unique_ptr<double[]> h_1_(new double[NY]);
    std::unique_ptr<double[]> h_2_(new double[NY]);
    for (int x = 0; x < NX; x++) {
      for (int y = 0; y < NY; y++) {
        int i = y * NX + x;
        b_[y] = b[i];
        h_1_[y] = h_1[i];
        h_2_[y] = h_2[i];
      }
      solve_2d(h_.get(), h_1_.get(), h_2_.get(), b_.get(), dy, dt, NY);
      for (int y = 0; y < NY; y++) {
        int i = y * NX + x;
        h[i] = (h[i] + h_[y]) * 0.5;
      }
    }

    // output
    if (it % output_intvl == 0) {
      char filename[20];
      std::snprintf(filename, 20, "out3d/h_%04d.txt", it / output_intvl);
      std::ofstream ofs(filename);
      for(int y = 0; y < NY; y++) {
        for(int x = 0; x < NX; x++) {
          int i = y * NX + x;
          ofs << dx * x << "," << dy * y << "," << h[i] << std::endl;
        }
      }
      std::cout << "Output to: " << filename << std::endl;
    }

    // h <- h_2
    // h_1 <- h
    // h_2 <- h_1
    h_1.swap(h_2);
    h.swap(h_1);
  }
}

int main(int argc, char* argv[]) {
  /* shallow_water_2d(1000); */
  shallow_water_3d(100, 100);
}
