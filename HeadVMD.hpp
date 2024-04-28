#include "header.hpp"
//����� python-������: Vinicius Rezende Carvalho
template<typename SignalDT>
class F;        //���������������� ��������� ���

class Ud;       //������������ (������������� �������, ����)
class Fregs;    //���������� ��������� �������
class alpha;    //�����-�������� ��� ������ ������������
class NumbIter; //����������� ���������� ��������
class Epsilon;  //������� ��������� �������� (������������ �������� ������������� �������� ������(����� ��������) � ������� ��������)
class Omega;    //�����-�������� ��� ������ ������������
class Lamb;     //������-��������
class Sum_UK;   //����� ������������
class VMDDecomposition; //����� ������������


template<typename SignalDT>
class F
{
private:
    vector<complex<double>> f_hat_plus;
    int T;
    double fs;
public:
    //������������ � ����������� �������, ��������� �������� ������� ������� (fs) � ����� ������� �
    F(vector<SignalDT> f);

    vector<complex<double>> operator()() const;
    complex<double>& operator[](int i);
    F<SignalDT>& operator=(vector<complex<double>> value);
    size_t size() const;
    vector<complex<double>> getF() const;
    int getT() const;
    double getFs() const;
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
    Fregs& operator=(vector<double> value);
    size_t size() const;
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
    alpha& operator=(vector<double> value);
    size_t size() const;
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
};
class Sum_UK
{
private:
    vector<complex<double>> sum_uk;
public:
    //������� ����� �����
    Sum_UK(int T);

    vector<complex<double>> operator()() const;
    complex<double>& operator[](int i);
    //���������� ����� ����� (��� ������ � �����������)
    void update_sum_first(int T, int K, int n, Ud& u_hat_plus);
    void update_sum_uk(int T, int k, int n, Ud& u_hat_plus);
};
class Ud
{
private:
    vector<vector<vector<complex<double>>>> u_hat_plus; //�������� ����� ����� ������ ��������
    bool check_finite = true;
    vector<vector<complex<double>>> u_hat_final; //���� �� ��������� �������
    vector<vector<complex<double>>> u_final; //���� �� ��������� �������
public:
    Ud(int Niter, int K, vector<double> freqs);

    vector<vector<vector<complex<double>>>> operator()() const;
    vector<vector<complex<double>>>& operator[](int i);
    //�������� ������������
    bool isFinite();
    //���������� ����
    void update_uhat(int T, int k, int n, vector<complex<double>>& f_hat_plus, Sum_UK& sum_uk, Lamb& lambda_hat, alpha& Alpha, Fregs& freqs, Omega& omega_plus);
    //��������� �������� ������������
    void u(int N, int T, int k);

    vector<vector<complex<double>>> getU() const;
    vector<vector<complex<double>>> getU_hat() const;
};
class VMDDecomposition
{
private:
    vector<vector<complex<double>>> u_hat;//���� �� ��������� �������
    vector<vector<complex<double>>> u;//���� �� ��������� �������
    vector<vector<double>> omega;//����������� �������
public:
    VMDDecomposition(vector<vector<complex<double>>> u_hat, vector<vector<complex<double>>> u, vector<vector<double>> omega);
    vector<vector<complex<double>>> GetUhat() const;
    vector<vector<complex<double>>> GetU() const;
    vector<vector<double>> GetOmega() const;
};

template<typename SignalDT> 
VMDDecomposition VMD(std::vector<SignalDT> f, double alpha_var, double tau, int k, bool DC, int init, double tol);
