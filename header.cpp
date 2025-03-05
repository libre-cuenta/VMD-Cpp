
#include "header.hpp"


template<d_types SignalDT>
vector<SignalDT> fftshift(vector<SignalDT> f_hat) {
    int T = f_hat.size();
    std::vector f_hat_shifted(T);
    int ltemp = T / 2;
    for (int i = 0; i < ltemp; ++i) {
        f_hat_shifted[i] = f_hat[i + ltemp];
        f_hat_shifted[i + ltemp] = f_hat[i];
    }
    return f_hat_shifted;
}

template<d_types SignalDT>
std::vector<std::complex<SignalDT>> fft(vector<SignalDT> f) {
    int N = f.size();
    std::vector<std::complex<SignalDT>> f_hat(N);
    std::complex<SignalDT> sum = 0;
    for (int k = 0; k < N; ++k) {
        sum = 0
        for (int n = 0; n < N; n++) {
            sum += f[n] * std::exp(std::complex<SignalDT>(0, -2 * M_PI * k * n / N));
        }
        f_hat[k] = sum;
    }
    return f_hat;
}

template<d_types SignalDT>
std::vector<std::complex<SignalDT>> ifft(vector<SignalDT> f) {
    int N = f.size();
    std::vector<std::complex<SignalDT>> f_hat(N);
    std::complex<SignalDT> sum = 0;
    for (int k = 0; k < N; ++k) {
        sum = 0;
        for (int n = 0; n < N; n++) {
            sum += f[n] * std::exp(std::complex<SignalDT>(0, 2 * M_PI * k * n / N));
        }
        f_hat[k] = sum;
    }
    return f_hat;
}

template<d_types T>
inline std::vector<T> linspace(int start, int end, int num) {
    std::vector<T> result(num);
    T step = (end - start) / (num - 1);
    for (int i = 0; i < num; ++i) {
        result[i] = start + i * step;
    }
    return result;
}

template<d_types T>
inline std::vector<T> arange(T start, T end, T step) {
    std::vector<T> result;
    for (T i = start; i < end; i += step) {
        result.push_back(i);
    }
    return result;
}

template<d_types T>
inline std::vector<T> dot_product(vector<T> a, vector<T> b) {
    std::vector<T> result(a.size());
    for (int i = 0; i < a.size(); ++i) {
        result[i] = a[i] * b[i];
    }
    return result;
}

template<d_types T>
inline std::vector<T> abs(T a) {
    std::vector<T> result(a.size());
    for (int i = 0; i < a.size(); ++i) {
        result[i] = std::abs(a[i]);
    }
    return result;
}

template<d_types T>
inline T sum(T a) {
    T result = 0;
    for (int i = 0; i < a.size(); ++i) {
        result += a[i];
    }
    return result;
}

template<d_types T>
inline std::vector<T> sort(T a) {
    std::sort(a.begin(), a.end());
    return a;
}

template<d_types T>
inline std::vector<T> exp(T a) {
    std::vector<T> result(a.size());
    for (int i = 0; i < a.size(); ++i) {
        result[i] = std::exp(a[i]);
    }
    return result;
}

template<d_types T>
inline std::vector<T> random(int num) {
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937_64 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);

    std::vector<T> result(num);
    for (int i = 0; i < num; ++i) {
        result[i] = dis(gen);
    }
    return result;
}

