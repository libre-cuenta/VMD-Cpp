#define _USE_MATH_DEFINES

#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
#include <functional>
#include <algorithm>
#include <array>
#include <string_view>
#include <numeric>

using namespace std;

std::vector<std::complex<double>> fftshift(std::vector<std::complex<double>> f_hat) {
    int T = f_hat.size();
    std::vector<std::complex<double>> f_hat_shifted(T);
    int ltemp = T / 2;
    for (int i = 0; i < ltemp; i++) {
        f_hat_shifted[i] = f_hat[i + ltemp];
        f_hat_shifted[i + ltemp] = f_hat[i];
    }
    return f_hat_shifted;
}

template<typename T> std::vector<std::complex<double>> fft(T f) {
    int N = f.size();
    std::vector<std::complex<double>> f_hat(N);
    for (int k = 0; k < N; k++) {
        std::complex<double> sum = 0;
        for (int n = 0; n < N; n++) {
            sum += f[n] * std::exp(std::complex<double>(0, -2 * M_PI * k * n / N));
        }
        f_hat[k] = sum;
    }
    return f_hat;
}

template<typename T> std::vector<std::complex<double>> ifft(T f) {
    int N = f.size();
    std::vector<std::complex<double>> f_hat(N);
    for (int k = 0; k < N; k++) {
        std::complex<double> sum = 0;
        for (int n = 0; n < N; n++) {
            sum += f[n] * std::exp(std::complex<double>(0, 2 * M_PI * k * n / N));
        }
        f_hat[k] = sum;
    }
    return f_hat;
}

std::vector<double> linspace(double start, double end, int num) {
    std::vector<double> result(num);
    double step = (end - start) / (num - 1);
    for (int i = 0; i < num; i++) {
        result[i] = start + i * step;
    }
    return result;
}

std::vector<double> arange(double start, double end, double step) {
    std::vector<double> result;
    for (double i = start; i < end; i += step) {
        result.push_back(i);
    }
    return result;
}

std::vector<double> dot_product(std::vector<double> a, std::vector<double> b) {
    std::vector<double> result(a.size());
    for (int i = 0; i < a.size(); i++) {
        result[i] = a[i] * b[i];
    }
    return result;
}

std::vector<double> abs(std::vector<double> a) {
    std::vector<double> result(a.size());
    for (int i = 0; i < a.size(); i++) {
        result[i] = std::abs(a[i]);
    }
    return result;
}

double sum(std::vector<double> a) {
    double result = 0;
    for (int i = 0; i < a.size(); i++) {
        result += a[i];
    }
    return result;
}

std::vector<double> sort(std::vector<double> a) {
    std::sort(a.begin(), a.end());
    return a;
}

std::vector<double> exp(std::vector<double> a) {
    std::vector<double> result(a.size());
    for (int i = 0; i < a.size(); i++) {
        result[i] = std::exp(a[i]);
    }
    return result;
}

std::vector<double> random(int num) {
    std::vector<double> result(num);
    for (int i = 0; i < num; i++) {
        result[i] = (double)rand() / RAND_MAX;
    }
    return result;
}


std::vector<std::vector<std::vector<complex<double>>>> VMD(std::vector<double> f, double alpha, double tau, int K, bool DC, int init, double tol) {
    if (f.size() % 2 != 0) {
        f.pop_back();
    }

    double fs = 1.0 / f.size();
    int ltemp = f.size() / 2;
    std::vector<double> fMirr(ltemp);
    for (int i = 0; i < ltemp; i++) {
        fMirr[i] = f[ltemp - i - 1];
    }
    fMirr.insert(fMirr.end(), f.begin(), f.end());
    for (int i = 0; i < ltemp; i++) {
        fMirr.push_back(f[f.size() - i - 1]);
    }

    int T = fMirr.size();
    std::vector<double> t = linspace(1, T, T);

    std::vector<double> freqs(T);
    for (int i = 0; i < T; i++) {
        freqs[i] = t[i] - 0.5 - (1.0 / T);
    }

    int Niter = 500;

    std::vector<double> Alpha(K, alpha);

    std::vector<std::complex<double>> f_hat = fftshift(fft(fMirr));
    std::vector<std::complex<double>> f_hat_plus = f_hat;
    for (int i = 0; i < T / 2; i++) {
        f_hat_plus[i] = 0;
    }

    std::vector<std::vector<double>> omega_plus(Niter, std::vector<double>(K));
    if (init == 1) {
        for (int i = 0; i < K; i++) {
            omega_plus[0][i] = (0.5 / K) * i;
        }
    }
    else if (init == 2) {
        std::vector<double> random_nums = random(K);
        std::vector<double> exp_vals;
        for (int i = 0; i < K; i++) 
        {
            exp_vals[i] = exp((log(fs) + (log(0.5) - log(fs))) * random_nums[i]);
        }
        std::vector<double> sorted_exp_vals = sort(exp_vals);
        for (int i = 0; i < K; i++) {
            omega_plus[0][i] = sorted_exp_vals[i];
        }
    }
    else {
        std::vector<double> zeros(K);
        omega_plus[0] = zeros;
    }

    if (DC) {
        omega_plus[0][0] = 0;
    }

    std::vector<std::vector<std::complex<double>>> lambda_hat(Niter, std::vector<std::complex<double>>(freqs.size()));

    double uDiff = tol + std::numeric_limits<double>::epsilon();
    int n = 0;
    double sum_uk = 0;

    std::vector<vector<std::vector<std::complex<double>>>> u_hat_plus(Niter, std::vector<std::vector<std::complex<double>>>(freqs.size(), std::vector<std::complex<double>>(K)));

    while (uDiff > tol && n < Niter - 1) {
        int k = 0;
        std::vector<double> sum_uk(T, 0.0);
        for (int i = 0; i < T; i++) {
            sum_uk[i] = u_hat_plus[n][i][K - 1].real() + sum_uk[i] - u_hat_plus[n][i][0].real();
        }

        for (int i = 0; i < T; i++) {
            u_hat_plus[n + 1][i][k] = (f_hat_plus[i] - sum_uk[i] - (lambda_hat[n][i] / 2.0)) / (1.0 + Alpha[k] * (freqs[i] - omega_plus[n][k]) * (freqs[i] - omega_plus[n][k]));
        }

        if (!DC) {
            double sum = 0.0;
            for (int i = T / 2; i < T; i++) {
                sum += std::abs(u_hat_plus[n + 1][i][k] * u_hat_plus[n + 1][i][k]);
            }
            omega_plus[n + 1][k] = std::inner_product(freqs.begin() + T / 2, freqs.end(), u_hat_plus[n + 1].begin() + T / 2, 0.0) / sum;
        }

        for (int k = 1; k < K; k++) {
            std::vector<double> sum_uk(T, 0.0);
            for (int i = 0; i < T; i++) {
                sum_uk[i] = u_hat_plus[n + 1][i][k - 1].real() + sum_uk[i] - u_hat_plus[n][i][k].real();
            }

            for (int i = 0; i < T; i++) {
                u_hat_plus[n + 1][i][k] = (f_hat_plus[i] - sum_uk[i] - (lambda_hat[n][i] / 2.0)) / (1.0 + Alpha[k] * (freqs[i] - omega_plus[n][k]) * (freqs[i] - omega_plus[n][k]));
            }

            double sum = 0.0;
            for (int i = T / 2; i < T; i++) {
                sum += std::abs(u_hat_plus[n + 1][i][k] * u_hat_plus[n + 1][i][k]);
            }
            omega_plus[n + 1][k] = std::inner_product(freqs.begin() + T / 2, freqs.end(), u_hat_plus[n + 1].begin() + T / 2, 0.0) / sum;
        }

        for (int i = 0; i < T; i++) {
            lambda_hat[n + 1][i] = lambda_hat[n][i] + tau * (std::accumulate(u_hat_plus[n + 1][i].begin(), u_hat_plus[n + 1][i].end(), 0.0) - f_hat_plus[i]);
        }

        n = n + 1;

        uDiff = std::numeric_limits<double>::epsilon();
        for (int i = 0; i < K; i++) {
            uDiff = uDiff + (1.0 / T) * std::inner_product(u_hat_plus[n][i].begin(), u_hat_plus[n][i].end(), u_hat_plus[n - 1][i].begin(), std::complex<double>(0.0, 0.0), std::plus<std::complex<double>>(), [](const std::complex<double>& a, const std::complex<double>& b) { return std::conj(a - b) * (a - b); }).real();
        }
        uDiff = std::abs(uDiff);
    }

    Niter = std::min(Niter, n);
    std::vector<std::vector<std::complex<double>>> omega(Niter, std::vector<complex<double>>(K, 0.0));
    for (int i = 0; i < Niter; i++) {
        for (int k = 0; k < K; k++) {
            omega[i][k] = omega_plus[i][k];
        }
    }

    std::vector<int> idxs(T / 2);
    std::iota(idxs.begin(), idxs.end(), 1);
    std::reverse(idxs.begin(), idxs.end());

    std::vector<std::vector<std::complex<double>>> u_hat(T, std::vector<std::complex<double>>(K, std::complex<double>(0.0, 0.0)));
    for (int i = T / 2; i < T; i++) {
        u_hat[i] = u_hat_plus[Niter - 1][i];
    }
    for (int i = 0; i < T / 2; i++) {
        for (int j = 0; j < K; j++) {
            u_hat[idxs[i]][j] = std::conj(u_hat_plus[Niter - 1][i][j]);
        }
    }
    for (int j = 0; j < K; j++) {
        u_hat[0][j] = std::conj(u_hat[T - 1][j]);
    }

    std::vector<std::vector<complex<double>>> u(K, std::vector<complex<double>>(t.size(), 0.0));
    for (int i = 0; i < T; i++) {
        for (int k = 0; k < K; k++) {
            u[k][i] = 0;
        }
    }
    for (int k = 0; k < K; k++) {
        u[k] = ifft(vector<complex<double>>(fftshift(u_hat[k])));
    }

    std::vector<std::vector<complex<double>>> u_final(K, vector<complex<double>>(T / 2));
    for (int k = 0; k < K; k++) {
        for (int i = T / 4; i < 3 * T / 4; i++) {
            u_final[k][i - T / 4] = u[k][i];
        }
    }

    std::vector<std::vector<complex<double>>> u_hat_final(K, vector<complex<double>>(T / 2));
    for (int k = 0; k < K; k++) {
        u_hat_final[k] = fftshift(fft(u[k]));
    }

    vector<vector<vector<complex<double>>>> output = { u_final , u_hat_final , omega };

    return output;
}


