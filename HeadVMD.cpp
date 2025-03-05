
#include "HeadVMD.hpp"

template<d_types SignalDT> F<SignalDT>::F(vector<SignalDT> f)
{
    //if (f.size() % 2 != 0) {
    //    f.pop_back();
    //}
    //double fs = 1.0 / f.size();
    //int ltemp = f.size() / 2;
    //std::vector<double> fMirr(ltemp);
    //for (int i = 0; i < ltemp; ++i) {
    //    fMirr[i] = f[ltemp - i - 1];
    //}???
    //fMirr.insert(fMirr.end(), f.begin(), f.end());
    //for (int i = 0; i < ltemp; ++i) {
    //    fMirr.push_back(f[f.size() - i - 1]);
    //}
    //int T = fMirr.size();

    if (f.size() % 2 != 0) {
        f.pop_back();
    }

    fs = 1.0 / f.size();
    int ltemp = f.size() / 2;
    std::vector fMirr(ltemp);
    for (int i = 0; i < ltemp; ++i) {
        fMirr[i] = f[ltemp - i - 1];
    }
    fMirr.insert(fMirr.end(), f.begin(), f.end());
    for (int i = 0; i < ltemp; ++i) {
        fMirr.push_back(f[f.size() - i - 1]);
    }

    this->T = fMirr.size();

    //std::vector<std::complex<double>> f_hat = fftshift(fft(fMirr));
    //std::vector<std::complex<double>> f_hat_plus = f_hat;
    //for (int i = 0; i < T / 2; ++i) {
    //    f_hat_plus[i] = 0;
    //}
    std::vector<std::complex<SignalDT>> f_hat = fftshift<complex<SignalDT>>(fft<SignalDT>(fMirr));
    f_hat_plus = f_hat;
    for (int i = 0; i < T / 2; ++i) {
        f_hat_plus[i] = 0;
    }

};
template<d_types SignalDT> vector<complex<SignalDT>> F<SignalDT>:: operator()() const
{
    return this->f_hat_plus;
};
template<d_types SignalDT> complex<SignalDT>& F<SignalDT>:: operator[](int i) const
{
    return f_hat_plus[i];
};
template<d_types SignalDT> size_t F<SignalDT>::size() const
{
    return f_hat_plus.size();
};
template<d_types SignalDT> vector<complex<SignalDT>> F<SignalDT>::getF() const
{
    return this->f_hat_plus;
}
template<d_types SignalDT> int F<SignalDT>::getT() const
{
    return T;
};
template<d_types SignalDT> SignalDT F<SignalDT>::getFs() const
{
    return fs;
};


Fregs::Fregs(int T)
{
    //std::vector<double> t = linspace(1, T, T);
    //std::vector<double> freqs(T);
    //for (int i = 0; i < T; ++i) {
    //    freqs[i] = t[i] - 0.5 - (1.0 / T);
    //}

    std::vector<double> t = linspace<double>(1, T, T);

    freqs.assign(T, 0.0);
    for (int i = 0; i < T; ++i) {
        freqs[i] = t[i] - 0.5 - (1.0 / T);
    }
};
vector<double> Fregs::operator()() const
{
    return this->freqs;
};
double Fregs:: operator[](int i) const
{
    return freqs[i];
};
size_t Fregs::size() const
{
    return freqs.size();
};


alpha::alpha(int K, double a)
{
    //std::vector<double> Alpha(K, alpha);
    Alpha.assign(K, a);
};
vector<double> alpha::operator()() const
{
    return this->Alpha;
};
double alpha::operator[](int i) const
{
    return Alpha[i];
};
size_t alpha::size() const
{
    return Alpha.size();
};


NumbIter::NumbIter() {};
int NumbIter::operator()() const
{
    return N;
};
NumbIter& NumbIter::operator += (const int& counter)
{
    this->N += counter;
    return *this;   // возвращаем ссылку на текущий объект
};
void NumbIter::SetMinToN(int n)
{
    this->N = min(n, N);
};


//double uDiff = tol + std::numeric_limits<double>::epsilon();
 Epsilon::Epsilon(double tol) : udiff(numeric_limits<double>::epsilon() + tol) {};
 double& Epsilon::operator()()
{
    return this->udiff;
};
 void Epsilon::update_udiff(int T, int K, int n, Ud& u_hat_plus)
{
    for (int i = 0; i < K; ++i) {
        this->udiff = this->udiff + (1.0 / T) * std::inner_product(u_hat_plus()[n][i].begin(), u_hat_plus()[n][i].end(), u_hat_plus()[n - 1][i].begin(), std::complex<double>(0.0, 0.0), std::plus<std::complex<double>>(), [](const std::complex<double>& a, const std::complex<double>& b) { return std::conj(a - b) * (a - b); }).real();
    }
    this->udiff = std::abs(this->udiff);
};


 Omega::Omega(int Niter, int K, int init, double fs, bool DC)
{
    this->omega_plus.assign(Niter, std::vector<double>(K));
    if (init == 1) {
        for (int i = 0; i < K; ++i) {
            this->omega_plus[0][i] = (0.5 / K) * i;
        }
    }
    else if (init == 2) {
        std::vector<double> random_nums = random(K);
        std::vector<double> exp_vals;
        for (int i = 0; i < K; ++i)
        {
            exp_vals[i] = exp((log(fs) + (log(0.5) - log(fs))) * random_nums[i]);
        }
        std::vector<double> sorted_exp_vals = sort(exp_vals);
        for (int i = 0; i < K; ++i) {
            this->omega_plus[0][i] = sorted_exp_vals[i];
        }
    }
    else {
        std::vector<double> zeros(K);
        this->omega_plus[0] = zeros;
    }

    if (DC) {
        this->omega_plus[0][0] = 0;
    }
};
 vector<vector<double>> Omega::operator()() const
{
    return this->omega_plus;
};
 vector<double>& Omega::operator[](int i)
{
    return omega_plus[i];
};
 bool Omega::isFinite()
{
    return check_finite;
};
 void Omega::update_omega_plus(int T, int k, int n, Ud& u_hat_plus, Fregs& fregs)
{
    double sum = 0.0;
    for (int i = T / 2; i < T; ++i) {
        sum += std::abs(u_hat_plus[n + 1][i][k] * u_hat_plus[n + 1][i][k]);
    }
    this->omega_plus[n + 1][k] = std::inner_product(fregs().begin() + T / 2, fregs().end(), u_hat_plus[n + 1].begin() + T / 2, 0.0) / sum;
    if (!isfinite(this->omega_plus[n + 1][k]))
    {
        this->check_finite = false;
    }
};
 void Omega::count_omega(int Niter, int K)
{
    this->omega.assign(Niter, std::vector<complex<double>>(K, 0.0));
    for (int i = 0; i < Niter; ++i) {
        for (int k = 0; k < K; ++k) {
            this->omega[i][k] = this->omega_plus[i][k];
        }
    }
};


 Lamb::Lamb(int Niter, vector<double> freqs)
{
    //std::vector<std::vector<std::complex<double>>> lambda_hat(Niter, std::vector<std::complex<double>>(freqs.size()));
    lambda.assign(Niter, std::vector<std::complex<double>>(freqs.size()));
};
 vector<vector<complex<double>>> Lamb::operator()() const
{
    return this->lambda;
};
 vector<complex<double>>& Lamb::operator[](int i)
{
    return lambda[i];
};
 void Lamb::update_lamb(int T, int n, double tau, Ud& u_hat_plus, vector<complex<double>>& f_hat_plus)
{
    for (int i = 0; i < T; ++i) {
        lambda[n + 1][i] = lambda[n][i] + tau * (std::accumulate(u_hat_plus[n + 1][i].begin(), u_hat_plus[n + 1][i].end(), 0.0) - f_hat_plus[i]);
    }
};


 template<d_types SignalDT> Sum_UK<SignalDT>::Sum_UK(int T)
{
    sum_uk.assign(T, complex<SignalDT>(0.0, 0.0));
};
 template<d_types SignalDT>
 vector<complex<SignalDT>> Sum_UK<SignalDT>::operator()() const
{
    return this->sum_uk;
};
 template<d_types SignalDT>
 complex<SignalDT>& Sum_UK<SignalDT>::operator[](int i)
{
    return sum_uk[i];
};
 template<d_types SignalDT>
 void Sum_UK<SignalDT>::update_sum_first(int T, int K, int n, Ud<SignalDT>& u_hat_plus)
{
    for (int i = 0; i < T; ++i) {
        sum_uk[i] = u_hat_plus[n][i][K - 1] + sum_uk[i] - u_hat_plus[n][i][0];
    }
};
 template<d_types SignalDT>
 void Sum_UK<SignalDT>::update_sum_uk(int T, int k, int n, Ud<SignalDT>& u_hat_plus)
{
    for (int i = 0; i < T; ++i) {
        sum_uk[i] = u_hat_plus[n + 1][i][k - 1] + sum_uk[i] - u_hat_plus[n][i][k];
    }
};

 template<d_types SignalDT>
 Ud<SignalDT>::Ud(int Niter, int K, vector<SignalDT> freqs)
 {
     //std::vector<vector<std::vector<std::complex<double>>>> u_hat_plus(Niter, std::vector<std::vector<std::complex<double>>>(freqs.size(), std::vector<std::complex<double>>(K)));
     u_hat_plus.assign(Niter, vector<vector<complex<SignalDT>>>(freqs.size(), vector<complex<SignalDT>>(K)));
 };
 template<d_types SignalDT>
 vector<vector<vector<complex<SignalDT>>>> Ud<SignalDT>::operator()() const
{
    return this->u_hat_plus;
};
 template<d_types SignalDT>
 vector<vector<complex<SignalDT>>>& Ud<SignalDT>::operator[](int i)
{
    return u_hat_plus[i];
};
 template<d_types SignalDT>
 bool Ud<SignalDT>::isFinite()
{
    return check_finite;
};
 template<d_types SignalDT>
 void Ud<SignalDT>::update_uhat(int T, int k, int n, vector<complex<SignalDT>>& f_hat_plus, Sum_UK<SignalDT>& sum_uk, Lamb& lambda_hat, alpha& Alpha, Fregs& freqs, Omega& omega_plus)
{
    for (int i = 0; i < T; ++i) {
        u_hat_plus[n + 1][i][k] = (f_hat_plus[i] - sum_uk[i] - (lambda_hat[n][i] / 2.0)) / (1.0 + Alpha[k] * (freqs[i] - omega_plus[n][k]) * (freqs[i] - omega_plus[n][k]));
        if (!isfinite(u_hat_plus[n + 1][i][k]))
        {
            check_finite = false;
            break;
        }
    }
};
 template<d_types SignalDT>
 void Ud<SignalDT>::u(int N, int T, int K)
{
    std::vector<int> idxs(T / 2);
    std::iota(idxs.begin(), idxs.end(), 1);
    std::reverse(idxs.begin(), idxs.end());

    std::vector<std::vector<std::complex<SignalDT>>> u_hat(T, std::vector<std::complex<SignalDT>>(K, std::complex<SignalDT>(0.0, 0.0)));
    for (int i = T / 2; i < T; ++i) {
        u_hat[i] = u_hat_plus[N - 1][i];
    }
    for (int i = 0; i < T / 2; ++i) {
        for (int j = 0; j < K; ++j) {
            u_hat[idxs[i]][j] = std::conj(u_hat_plus[N - 1][i][j]);
        }
    }
    for (int j = 0; j < K; ++j) {
        u_hat[0][j] = std::conj(u_hat[T - 1][j]);
    }

    std::vector<std::vector<complex<SignalDT>>> u_buf(K, std::vector<complex<SignalDT>>(T, 0.0));
    for (int i = 0; i < T; ++i) {
        for (int k = 0; k < K; ++k) {
            u_buf[i][k] = 0;
        }
    }
    for (int k = 0; k < K; ++k) {
        u_buf[k] = ifft<complex<SignalDT>>(vector<complex<SignalDT>>(fftshift<complex<SignalDT>>(u_hat[k])));
    }

    u_final.assign(K, vector<complex<SignalDT>>(T / 2));
    for (int k = 0; k < K; ++k) {
        for (int i = T / 4; i < 3 * T / 4; ++i) {
            u_final[k][i - T / 4] = u_buf[k][i];
        }
    }

    u_hat_final.assign(K, vector<complex<SignalDT>>(T / 2));
    for (int k = 0; k < K; ++k) {
        u_hat_final[k] = fftshift<complex<SignalDT>>(fft<complex<SignalDT>>(u_buf[k]));
    }
};
 template<d_types SignalDT>
 vector<vector<complex<SignalDT>>> Ud<SignalDT>::getU() const
{
    return u_final;
};
 template<d_types SignalDT>
 vector<vector<complex<SignalDT>>> Ud<SignalDT>::getU_hat() const
{
    return u_hat_final;
};

 template<d_types SignalDT>
VMDDecomposition<SignalDT>::VMDDecomposition(vector<vector<complex<SignalDT>>> u_hat, vector<vector<complex<SignalDT>>> u, vector<vector<SignalDT>> omega) : u_hat(u_hat), u(u), omega(omega) {};
template<d_types SignalDT>
vector<vector<complex<SignalDT>>> VMDDecomposition<SignalDT>::GetUhat() const
{
    return u_hat;
};
template<typename SignalDT>
vector<vector<complex<SignalDT>>> VMDDecomposition<SignalDT>::GetU() const
{
    return u;
};
template<typename SignalDT>
vector<vector<SignalDT>> VMDDecomposition<SignalDT>::GetOmega() const
{
    return omega;
};
