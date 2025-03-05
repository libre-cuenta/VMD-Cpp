#include "header.hpp"
//автор python-версии: Vinicius Rezende Carvalho
template<d_types SignalDT>
class F;        //предобработанный временной ряд

template<d_types SignalDT>
class Ud;       //декомпозиции (внутримодовые функции, моды)
class Fregs;    //дискретная частотная область
class alpha;    //альфа-значение для каждой декомпозиции
class NumbIter; //ограничение количества итераций
class Epsilon;  //условие остановки итераций (максимальное значение относительной разности нового(после итерации) и старого значения)
class Omega;    //омега-значение для каждой декомпозиции
class Lamb;     //лямбда-значение

template<d_types SignalDT>
class Sum_UK;   //сумма декомпозиций

template<d_types SignalDT>
class VMDDecomposition; //вывод декомпозиций


template<d_types SignalDT>
class F
{
private:
    vector<complex<SignalDT>> f_hat_plus;
    int T;
    double fs;
public:
    //иницализация и предобратка сигнала, получение значение средней частоты (fs) и длины сигнала Т
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
    //генерация дискретной области частот
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
    //создания массива альфа-значений для количества мод К
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
    //инициализация условия остановки
    Epsilon(double tol);

    double& operator()();
    //обновление условия остановки после очередной итерации
    void update_udiff(int T, int k, int n, Ud& u_hat_plus);
    ~Epsilon() = default;
};


class Omega
{
private:
    vector<vector<double>> omega_plus;//центральные частоты после каждой итерации
    vector<vector<complex<double>>> omega; //центральные частоты
    bool check_finite = true;
public:
    //инициализация значения омега для каждой моды
    Omega(int Niter, int K, int init, double fs, bool DC);

    vector<vector<double>> operator()() const;
    vector<double>& operator[](int i);
    //проверка переполнения
    bool isFinite();
    //обновления значений омега
    void update_omega_plus(int T, int k, int n, Ud& u_hat_plus, Fregs& fregs);
    void count_omega(int Niter, int K);
    ~Omega() = default;
};


class Lamb
{
private:
    vector<vector<complex<double>>> lambda;
public:
    //инициализация лямбда
    Lamb(int Niter, vector<double> freqs);

    vector<vector<complex<double>>> operator()() const;
    vector<complex<double>>& operator[](int i);
    //обновление лямбда
    void update_lamb(int T, int n, double tau, Ud& u_hat_plus, vector<complex<double>>& f_hat_plus);
    ~Lamb() = default;
};


template<d_types SignalDT>
class Sum_UK
{
private:
    vector<complex<SignalDT>> sum_uk;
public:
    //подсчёт суммы модов
    Sum_UK<SignalDT>(int T);

    vector<complex<SignalDT>> operator()() const;
    complex<SignalDT>& operator[](int i);
    //обновление суммы модов (для первой и последующих)
    void update_sum_first(int T, int K, int n, Ud& u_hat_plus);
    void update_sum_uk(int T, int k, int n, Ud& u_hat_plus);
    ~Sum_UK() = default;
};

template<d_types SignalDT>
class Ud
{
private:
    vector<vector<vector<complex<SignalDT>>>> u_hat_plus; //значения модов после каждой итерации
    bool check_finite = true;
    vector<vector<complex<SignalDT>>> u_hat_final; //моды во частотной области
    vector<vector<complex<SignalDT>>> u_final; //моды во временной области
public:
    Ud(int Niter, int K, vector<SignalDT> freqs);

    vector<vector<vector<complex<SignalDT>>>> operator()() const;
    vector<vector<complex<SignalDT>>>& operator[](int i);
    //проверка переполнения
    bool isFinite();
    //обновление моды
    void update_uhat(int T, int k, int n, vector<complex<SignalDT>>& f_hat_plus, Sum_UK<SignalDT>& sum_uk, Lamb& lambda_hat, alpha& Alpha, Fregs& freqs, Omega& omega_plus);
    //получение итоговой декомпозиции
    void u(int N, int T, int k);

    vector<vector<complex<SignalDT>>> getU() const;
    vector<vector<complex<SignalDT>>> getU_hat() const;
    ~Ud() = default;
};

template<d_types SignalDT>
class VMDDecomposition
{
private:
    vector<vector<complex<SignalDT>>> u_hat;//моды во частотной области
    vector<vector<complex<SignalDT>>> u;//моды во временной области
    vector<vector<SignalDT>> omega;//центральные частоты
public:
    VMDDecomposition(vector<vector<complex<SignalDT>>> u_hat, vector<vector<complex<SignalDT>>> u, vector<vector<SignalDT>> omega);
    vector<vector<complex<SignalDT>>> GetUhat() const;
    vector<vector<complex<SignalDT>>> GetU() const;
    vector<vector<SignalDT>> GetOmega() const;
    ~VMDDecomposition() = default;
};

