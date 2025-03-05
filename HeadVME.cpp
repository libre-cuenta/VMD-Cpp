
//F<complex<double>>& f_hat_onesided, Ud& u_hat_d, OmegaAxis& omega_axis, Omega& omega_d, Lamb& lamb
#include "HeadVME.hpp"


template<d_types SignalDT> F<SignalDT>::F(vector<SignalDT> signal)
{
    int save_T = signal.size();
    int T = save_T;

    vector<SignalDT> f_mir(2 * T);
    for (int i = 0; i < T / 2; i++) {
        f_mir[i] = signal[i];
    }
    for (int i = T / 2; i < 3 * T / 2; i++) {
        f_mir[i] = signal[T - i - 1];
    }
    for (int i = 3 * T / 2; i < 2 * T; i++) {
        f_mir[i] = signal[i - T / 2];
    }


    vector<SignalDT> f(T);
    for (int i = 0; i < T; i++) {
        f[i] = f_mir[T - i - 1];
    }
    T = f.size();

    vector<complex<double>> f_hat(T);
    for (int i = 0; i < T; i++) {
        f_hat[i] = complex<double>(0, 0);
    }
    f_hat = fftshift<complex<double>>(fft<SignalDT>(f));
    f_hat_onesided = f_hat;
    for (int i = 0; i < T / 2; i++) {
        f_hat_onesided[i] = complex<double>(0, 0);
    }

};
template<d_types SignalDT> vector<complex<double>> F<SignalDT>::operator()() const
{
    return this->f_hat_onesided;
};
template<d_types SignalDT> complex<double> F<SignalDT>::operator[](int i) const
{
    return f_hat_onesided[i];
};
template<d_types SignalDT> vector<complex<double>> F<SignalDT>::getF() const
{
    return this->f_hat_onesided;
};
template<d_types SignalDT> size_t F<SignalDT>::size() const
{
    return f_hat_onesided.size();
};
template<d_types SignalDT> int F<SignalDT>::getT() const
{
    return T;
};


complex<double> Epsilon::count_udiff(complex<double> udiff, int T, int n, Ud& u_hat_d)
{
    //auto count_udiff = [&](double udiff) {
    //    complex<double> result = udiff;
    //    for (int i = 0; i < T; i++) {
    //        result += ((u_hat_d[n][i] - u_hat_d[n-1][i]) * conj(u_hat_d[n][i] - u_hat_d[n-1][i]) *= (1 / T));
    //    }
    //    return result;
    //};
    complex<double> result = udiff;
    for (int i = 0; i < T; i++) {
        result += ((u_hat_d[n][i] - u_hat_d[n - 1][i]) * conj(u_hat_d[n][i] - u_hat_d[n - 1][i]) *= (1.0 / T));
    }
    return result;
};
Epsilon::Epsilon(double tol) : udiff(tol + numeric_limits<double>::epsilon()) {};
complex<double> Epsilon::getUdiff() const
{
    return udiff;
};
void Epsilon::update_udiff(complex<double> udiff, int T, int n, Ud& u_hat_d)
{
    complex<double> udiff_buf;
    udiff_buf = count_udiff(udiff, T, n, u_hat_d);
    this->udiff = abs(udiff_buf);
};


OmegaAxis::OmegaAxis(int T)
{
    f.assign(T, 0.0);
    //vector<double> t(T);
    //for (int i = 0; i < T; i++) {
    //    t[i] = i / (double)T;
    //}
    std::vector<double> t = linspace(1, T, T);

    for (int i = 0; i < T; i++) {
        f[i] = t[i] - 0.5 - (1.0 / T);
    }
};
vector<double> OmegaAxis::operator()() const
{
    return this->f;
};
double OmegaAxis::operator[](int i) const
{
    return f[i];
};
size_t OmegaAxis::size() const
{
    return f.size();
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
    //N = min(N, n) - 1;
    this->N = min(n, N) - 1;
};


vector<complex<double>> Lamb::count_lamb(int T, int n, double tau, double Alpha, vector<complex<double>>& f_hat_onesided, Ud& u_hat_d, OmegaAxis& omega_axis, Omega& omega_d)
{
    //auto count_lamb = [&]() {
    //    vector<complex<double>> result(T);
    //    for (int i = 0; i < T; i++) {
    //        result[i] = lamb[n][i] + (tau * (f_hat_onesided[i] - (u_hat_d[n+1][i] + (((pow(Alpha, 2) * pow(omega_axis[i] - omega_d[n+1], 4)) * (f_hat_onesided[i] - u_hat_d[n+1][i])) / (2 * pow(Alpha, 2) + 1) * pow(omega_axis[i] - omega_d[n+1], 4)))));
    //    }
    //    return result;
    //};
    vector<complex<double>> result(T);
    for (int i = 0; i < T; i++) {
        result[i] = lamb[n][i] + (tau * (f_hat_onesided[i] - (u_hat_d[n + 1][i] + (((pow(Alpha, 2) * pow(omega_axis[i] - omega_d[n + 1], 4)) * (f_hat_onesided[i] - u_hat_d[n + 1][i])) / (2 * pow(Alpha, 2) + 1) * pow(omega_axis[i] - omega_d[n + 1], 4)))));
    }
    return result;
};
Lamb::Lamb(int N, int T)
{
    //vector<vector<complex<double>>> lamb(N, vector<complex<double>>(T));
    vector<complex<double>> zeros;
    zeros.assign(T, complex<double>(0, 0));
    lamb.assign(N, zeros);
};
void Lamb::update_lamb(int T, int n, double tau, double Alpha, vector<complex<double>>& f_hat_onesided, Ud& u_hat_d, OmegaAxis& omega_axis, Omega& omega_d)
{
    lamb[n + 1] = count_lamb(T, n, tau, Alpha, f_hat_onesided, u_hat_d, omega_axis, omega_d);
};
vector<vector<complex<double>>> Lamb::operator()() const
{
    return this->lamb;
};
vector<complex<double>>& Lamb::operator[](int i)
{
    return lamb[i];
};
size_t Lamb::size_1() const
{
    return lamb.size();
};
size_t Lamb::size_2() const
{
    return lamb[0].size();
};


complex<double> Omega::count_omegad(int T, int n, Ud& u_hat_d, OmegaAxis& omega_axis)
{
    //auto count_omegad = [&]() {
    //    complex<double> A = 0;
    //    complex<double> B = 0;
    //    for (int i = T/2; i < T; i++) {
    //        A += omega_axis[i] * pow(abs(u_hat_d[n+1][i]), 2);
    //        B += pow(abs(u_hat_d[n+1][i]), 2);
    //    }
    //    if (!isfinite(A / B))
    //    {
    //        check_finite = false;
    //    }
    //    return A / B;
    //};
    complex<double> A = 0;
    complex<double> B = 0;
    for (int i = T / 2; i < T; i++) {
        A += omega_axis[i] * pow(abs(u_hat_d[n + 1][i]), 2);
        B += pow(abs(u_hat_d[n + 1][i]), 2);
    }
    if (!isfinite(A / B))
    {
        check_finite = false;
    }
    return A / B;
};
Omega::Omega(int N, double omega_int, double fs)
{
    //vector<complex<double>> omega_d(N);
    //for (int i = 0; i < N; i++) {
    //    omega_d[i] = complex<double>(0, 0);
    //}
    //omega_d[0] = complex<double>(omega_int / fs, 0);
    omega_d.assign(N, complex<double>(0, 0));
    omega_d[0] = complex<double>(omega_int / fs, 0);

};
bool Omega::isFinite()
{
    return check_finite;
};
vector<complex<double>> Omega::operator()() const
{
    return this->omega_d;
};
complex<double>& Omega::operator[](int i)
{
    return omega_d[i];
};
void Omega::update_omega(int T, int n, Ud& u_hat_d, OmegaAxis& omega_axis)
{
    omega_d[n + 1] = count_omegad(T, n, u_hat_d, omega_axis);
};
void Omega::minOmega(int n)
{
    //vector<complex<double>> omega(N);
    //for (int i = 0; i < N; i++) {
    //    omega[i] = omega_d[i];
    //}
    omega_d.resize(n);
};

template<d_types SignalDT>
vector<complex<SignalDT>> Ud<SignalDT>::count_uhat(int T, int n, double Alpha, vector<complex<SignalDT>>& f_hat_onesided, OmegaAxis& omega_axis, Omega& omega_d, Lamb& lamb)
{
    //auto count_uhat = [&]() {
    //    vector<complex<double>> result(T);
    //    for (int i = 0; i < T; i++) {
    //        result[i] = (f_hat_onesided[i] + (u_hat_d[n][i] * pow(Alpha, 2) * pow(omega_axis[i] - omega_d[n], 4)) + (lamb[n][i] / 2.0)) / (((pow(Alpha, 2) + 1) * pow(omega_axis[i] - omega_d[n], 4)) * ((2 * pow(Alpha, 2) + 1) * pow(omega_axis[i] - omega_d[n], 4)));
    //        if (!isfinite(result[i]))
    //        {
    //            check_finite = false;
    //            break;
    //        }
    //    }
    //    return result;
    //};
    vector<complex<double>> result(T);
    for (int i = 0; i < T; i++) {
        result[i] = (f_hat_onesided[i] + (u_hat_d[n][i] * pow(Alpha, 2) * pow(omega_axis[i] - omega_d[n], 4)) + (lamb[n][i] / 2.0)) / (((pow(Alpha, 2) + 1) * pow(omega_axis[i] - omega_d[n], 4)) * ((2 * pow(Alpha, 2) + 1) * pow(omega_axis[i] - omega_d[n], 4)));
        if (!isfinite(result[i]))
        {
            check_finite = false;
            break;
        }
    }
    return result;
};
template<d_types SignalDT>
Ud<SignalDT>::Ud(int N, int T)
{
    //vector<vector<complex<double>>> u_hat_d(N, vector<complex<double>>(T));
    vector<complex<double>> zeros;
    zeros.assign(T, complex<double>(0, 0));
    u_hat_d.assign(N, zeros);
};
template<d_types SignalDT>
void Ud<SignalDT>::update_u_hat(int T, int n, double Alpha, vector<complex<SignalDT>>& f_hat_onesided, OmegaAxis& omega_axis, Omega& omega_d, Lamb& lamb)
{
    u_hat_d[n + 1] = count_uhat(T, n, Alpha, f_hat_onesided, omega_axis, omega_d, lamb);
};
template<d_types SignalDT>
bool Ud<SignalDT>::isFinite()
{
    return check_finite;
};
template<d_types SignalDT>
vector<vector<complex<SignalDT>>> Ud<SignalDT>::operator()() const
{
    return this->u_hat_d;
};
template<d_types SignalDT>
vector<complex<SignalDT>>& Ud<SignalDT>::operator[](int i)
{
    return u_hat_d[i];
};
template<d_types SignalDT>
size_t Ud<SignalDT>::size_1() const
{
    return u_hat_d.size();
};
template<d_types SignalDT>
size_t Ud<SignalDT>::size_2() const
{
    return u_hat_d[0].size();
};
template<d_types SignalDT>
void Ud<SignalDT>::ud(int T, int N)
{
    //vector<complex<double>> u_hatd(T);
    //for (int i = 0; i < T; i++) {
    //    u_hatd[i] = complex<double>(0, 0);
    //}
    //for (int i = T/2; i < T; i++) {
    //    u_hatd[i] = u_hat_d[N][i];
    //}
    //for (int i = T/2; i > 0; i--) {
    //    u_hatd[T - i] = conj(u_hat_d[N][i]);
    //}
    //u_hatd[0] = conj(u_hatd[T-1]);
    //vector<complex<double>> u_d(T);
    //for (int i = 0; i < T; i++) {
    //    u_d[i] = 0;
    //}
    //u_d = fftshift(ifft<vector<complex<double>>>(u_hatd));
    //vector<complex<double>> u_d_final(T/2);
    //for (int i = T/4; i < 3*T/4; i++) {
    //    u_d_final[i - T/4] = u_d[i];
    //}
    //vector<complex<double>> u_hatd_final = fftshift(fft<vector<complex<double>>>(u_d));
    vector<complex<SignalDT>> u_hatd(T);
    for (int i = 0; i < T; i++) {
        u_hatd[i] = complex<SignalDT>(0, 0);
    }
    for (int i = T / 2; i < T; i++) {
        u_hatd[i] = u_hat_d[N][i];
    }
    for (int i = T / 2; i > 0; i--) {
        u_hatd[T - i] = conj(u_hat_d[N][i]);
    }
    u_hatd[0] = conj(u_hatd[T - 1]);

    vector<complex<SignalDT>> u_d(T);
    for (int i = 0; i < T; i++) {
        u_d[i] = 0;
    }
    u_d = fftshift<complex<SignalDT>>(ifft<complex<SignalDT>>(u_hatd));
    u_d_final.assign(T / 2, complex<SignalDT>(0, 0));
    for (int i = T / 4; i < 3 * T / 4; i++) {
        u_d_final[i - T / 4] = u_d[i];
    }
    u_hatd_final = fftshift<complex<SignalDT>>(fft<complex<SignalDT>>(u_d));

};
template<d_types SignalDT>
vector<complex<SignalDT>> Ud<SignalDT>::getU_D() const
{
    return u_d_final;
};
template<d_types SignalDT>
vector<complex<SignalDT>> Ud<SignalDT>::getU_hatD() const
{
    return u_hatd_final;
};

template<d_types SignalDT>
VMEDecomposition<SignalDT>::VMEDecomposition(vector<complex<SignalDT>> u_hatd, vector<complex<SignalDT>> u_d, vector<complex<SignalDT>> omega_d) : u_hatd(u_hatd), u_d(u_d), omega_d(omega_d) {};

template<d_types SignalDT>
vector<complex<SignalDT>> VMEDecomposition<SignalDT>::GetUhatD() const
{
    return u_hatd;
};
template<d_types SignalDT>
vector<complex<SignalDT>> VMEDecomposition<SignalDT>::GetUD() const
{
    return u_d;
};
template<d_types SignalDT>
vector<complex<SignalDT>> VMEDecomposition<SignalDT>::GetOmegaD() const
{
    return omega_d;
};
