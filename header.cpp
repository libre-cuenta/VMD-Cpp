
#include "header.hpp"

template<typename SignalDT>
vector<SignalDT> fftshift(vector<SignalDT> f_hat) {
    int T = f_hat.size();
    std::vector f_hat_shifted(T);
    int ltemp = T / 2;
    for (int i = 0; i < ltemp; i++) {
        f_hat_shifted[i] = f_hat[i + ltemp];
        f_hat_shifted[i + ltemp] = f_hat[i];
    }
    return f_hat_shifted;
}

template<typename SignalDT> 
std::vector<std::complex<double>> fft(vector<SignalDT> f) {
    int N = f.size();
    std::vector<std::complex<double>> f_hat(N);
    std::complex<double> sum = 0;
    for (int k = 0; k < N; k++) {
        sum = 0
        for (int n = 0; n < N; n++) {
            sum += f[n] * std::exp(std::complex<double>(0, -2 * M_PI * k * n / N));
        }
        f_hat[k] = sum;
    }
    return f_hat;
}

template<typename SignalDT> 
std::vector<std::complex<double>> ifft(vector<SignalDT> f) {
    int N = f.size();
    std::vector<std::complex<double>> f_hat(N);
    std::complex<double> sum = 0;
    for (int k = 0; k < N; k++) {
        sum = 0;
        for (int n = 0; n < N; n++) {
            sum += f[n] * std::exp(std::complex<double>(0, 2 * M_PI * k * n / N));
        }
        f_hat[k] = sum;
    }
    return f_hat;
}

inline std::vector<double> linspace(int start, int end, int num) {
    std::vector<double> result(num);
    double step = (end - start) / (num - 1);
    for (int i = 0; i < num; i++) {
        result[i] = start + i * step;
    }
    return result;
}

template<typename T>
inline std::vector<double> arange(T start, T end, T step) {
    std::vector<double> result;
    for (double i = start; i < end; i += step) {
        result.push_back(i);
    }
    return result;
}

template<typename T>
inline std::vector<T> dot_product(vector<T> a, vector<T> b) {
    std::vector<double> result(a.size());
    for (int i = 0; i < a.size(); i++) {
        result[i] = a[i] * b[i];
    }
    return result;
}

template<typename T>
inline std::vector<double> abs(T a) {
    std::vector<double> result(a.size());
    for (int i = 0; i < a.size(); i++) {
        result[i] = std::abs(a[i]);
    }
    return result;
}

template<typename T>
inline double sum(T a) {
    double result = 0;
    for (int i = 0; i < a.size(); i++) {
        result += a[i];
    }
    return result;
}

template<typename T>
inline std::vector<double> sort(T a) {
    std::sort(a.begin(), a.end());
    return a;
}

template<typename T>
inline std::vector<double> exp(T a) {
    std::vector<double> result(a.size());
    for (int i = 0; i < a.size(); i++) {
        result[i] = std::exp(a[i]);
    }
    return result;
}

inline std::vector<double> random(int num) {
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937_64 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);

    std::vector<double> result(num);
    for (int i = 0; i < num; i++) {
        result[i] = dis(gen);
    }
    return result;
}

