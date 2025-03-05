#include "header.hpp"
//����� python-������: Vinicius Rezende Carvalho
template<d_types SignalDT>
class F;        //���������������� ��������� ���

template<d_types SignalDT>
class Ud;       //������������ (������������� �������, ����)
class Fregs;    //���������� ��������� �������
class alpha;    //�����-�������� ��� ������ ������������
class NumbIter; //����������� ���������� ��������
class Epsilon;  //������� ��������� �������� (������������ �������� ������������� �������� ������(����� ��������) � ������� ��������)
class Omega;    //�����-�������� ��� ������ ������������
class Lamb;     //������-��������

template<d_types SignalDT>
class Sum_UK;   //����� ������������

template<d_types SignalDT>
class VMDDecomposition; //����� ������������


template<d_types SignalDT>
class F
{
private:
    vector<complex<SignalDT>> f_hat_plus;
    int T;
    double fs;
public:
    //������������ � ����������� �������, ��������� �������� ������� ������� (fs) � ����� ������� �
    F(vector<SignalDT> f);

    vector<complex<SignalDT>> operator()() const;
    complex<SignalDT>& operator[](int i) const;
    F<SignalDT>& operator=(vector<complex<SignalDT>> value) = delete;
    size_t size() const;
    vector<complex<SignalDT>> getF() const;
    int getT() const;
    SignalDT getFs() const;
    ~F() = default;
};


class Fregs
{
private:
    vector<double> freqs;
public:
    //��������� ���������� ������� ������
    Fregs(int T);

    vector<double> operator()() const;
    double operator[](int i) const;
    Fregs& operator=(vector<double> value) = delete;
    size_t size() const;
    ~Fregs() = default;
};


class alpha
{
private:
    vector<double> Alpha;
public:
    //�������� ������� �����-�������� ��� ���������� ��� �
    alpha(int K, double a);

    vector<double> operator()() const;
    double operator[](int i) const;
    alpha& operator=(vector<double> value) = delete;
    size_t size() const;
    ~alpha() = default;
};


class NumbIter
{
private:
    //int Niter = 500;
    int N = 500;
public:
    NumbIter();

    int operator()() const;
    NumbIter& operator += (const int& counter);
    void SetMinToN(int n);
    ~NumbIter() = default;
};


class Epsilon
{
private:
    double udiff;
public:
    //double uDiff = tol + std::numeric_limits<double>::epsilon();
    //������������� ������� ���������
    Epsilon(double tol);

    double& operator()();
    //���������� ������� ��������� ����� ��������� ��������
    void update_udiff(int T, int k, int n, Ud& u_hat_plus);
    ~Epsilon() = default;
};


class Omega
{
private:
    vector<vector<double>> omega_plus;//����������� ������� ����� ������ ��������
    vector<vector<complex<double>>> omega; //����������� �������
    bool check_finite = true;
public:
    //������������� �������� ����� ��� ������ ����
    Omega(int Niter, int K, int init, double fs, bool DC);

    vector<vector<double>> operator()() const;
    vector<double>& operator[](int i);
    //�������� ������������
    bool isFinite();
    //���������� �������� �����
    void update_omega_plus(int T, int k, int n, Ud& u_hat_plus, Fregs& fregs);
    void count_omega(int Niter, int K);
    ~Omega() = default;
};


class Lamb
{
private:
    vector<vector<complex<double>>> lambda;
public:
    //������������� ������
    Lamb(int Niter, vector<double> freqs);

    vector<vector<complex<double>>> operator()() const;
    vector<complex<double>>& operator[](int i);
    //���������� ������
    void update_lamb(int T, int n, double tau, Ud& u_hat_plus, vector<complex<double>>& f_hat_plus);
    ~Lamb() = default;
};


template<d_types SignalDT>
class Sum_UK
{
private:
    vector<complex<SignalDT>> sum_uk;
public:
    //������� ����� �����
    Sum_UK<SignalDT>(int T);

    vector<complex<SignalDT>> operator()() const;
    complex<SignalDT>& operator[](int i);
    //���������� ����� ����� (��� ������ � �����������)
    void update_sum_first(int T, int K, int n, Ud& u_hat_plus);
    void update_sum_uk(int T, int k, int n, Ud& u_hat_plus);
    ~Sum_UK() = default;
};

template<d_types SignalDT>
class Ud
{
private:
    vector<vector<vector<complex<SignalDT>>>> u_hat_plus; //�������� ����� ����� ������ ��������
    bool check_finite = true;
    vector<vector<complex<SignalDT>>> u_hat_final; //���� �� ��������� �������
    vector<vector<complex<SignalDT>>> u_final; //���� �� ��������� �������
public:
    Ud(int Niter, int K, vector<SignalDT> freqs);

    vector<vector<vector<complex<SignalDT>>>> operator()() const;
    vector<vector<complex<SignalDT>>>& operator[](int i);
    //�������� ������������
    bool isFinite();
    //���������� ����
    void update_uhat(int T, int k, int n, vector<complex<SignalDT>>& f_hat_plus, Sum_UK<SignalDT>& sum_uk, Lamb& lambda_hat, alpha& Alpha, Fregs& freqs, Omega& omega_plus);
    //��������� �������� ������������
    void u(int N, int T, int k);

    vector<vector<complex<SignalDT>>> getU() const;
    vector<vector<complex<SignalDT>>> getU_hat() const;
    ~Ud() = default;
};

template<d_types SignalDT>
class VMDDecomposition
{
private:
    vector<vector<complex<SignalDT>>> u_hat;//���� �� ��������� �������
    vector<vector<complex<SignalDT>>> u;//���� �� ��������� �������
    vector<vector<SignalDT>> omega;//����������� �������
public:
    VMDDecomposition(vector<vector<complex<SignalDT>>> u_hat, vector<vector<complex<SignalDT>>> u, vector<vector<SignalDT>> omega);
    vector<vector<complex<SignalDT>>> GetUhat() const;
    vector<vector<complex<SignalDT>>> GetU() const;
    vector<vector<SignalDT>> GetOmega() const;
    ~VMDDecomposition() = default;
};

