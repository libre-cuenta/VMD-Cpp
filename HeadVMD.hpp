#include "header.hpp"
//автор python-версии: Vinicius Rezende Carvalho
template<typename SignalDT>
class F;        //предобработанный временной ряд

class Ud;       //декомпозиции (внутримодовые функции, моды)
class Fregs;    //дискретная частотная область
class alpha;    //альфа-значение для каждой декомпозиции
class NumbIter; //ограничение количества итераций
class Epsilon;  //условие остановки итераций (максимальное значение относительной разности нового(после итерации) и старого значения)
class Omega;    //омега-значение для каждой декомпозиции
class Lamb;     //лямбда-значение
class Sum_UK;   //сумма декомпозиций
class VMDDecomposition; //вывод декомпозиций


template<typename SignalDT>
class F
{
private:
    vector<complex<double>> f_hat_plus;
    int T;
    double fs;
public:
    //иницализация и предобратка сигнала, получение значение средней частоты (fs) и длины сигнала Т
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
    //генерация дискретной области частот
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
    //создания массива альфа-значений для количества мод К
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
    //инициализация условия остановки
    Epsilon(double tol);

    double& operator()();
    //обновление условия остановки после очередной итерации
    void update_udiff(int T, int k, int n, Ud& u_hat_plus);
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
};
class Sum_UK
{
private:
    vector<complex<double>> sum_uk;
public:
    //подсчёт суммы модов
    Sum_UK(int T);

    vector<complex<double>> operator()() const;
    complex<double>& operator[](int i);
    //обновление суммы модов (для первой и последующих)
    void update_sum_first(int T, int K, int n, Ud& u_hat_plus);
    void update_sum_uk(int T, int k, int n, Ud& u_hat_plus);
};
class Ud
{
private:
    vector<vector<vector<complex<double>>>> u_hat_plus; //значения модов после каждой итерации
    bool check_finite = true;
    vector<vector<complex<double>>> u_hat_final; //моды во частотной области
    vector<vector<complex<double>>> u_final; //моды во временной области
public:
    Ud(int Niter, int K, vector<double> freqs);

    vector<vector<vector<complex<double>>>> operator()() const;
    vector<vector<complex<double>>>& operator[](int i);
    //проверка переполнения
    bool isFinite();
    //обновление моды
    void update_uhat(int T, int k, int n, vector<complex<double>>& f_hat_plus, Sum_UK& sum_uk, Lamb& lambda_hat, alpha& Alpha, Fregs& freqs, Omega& omega_plus);
    //получение итоговой декомпозиции
    void u(int N, int T, int k);

    vector<vector<complex<double>>> getU() const;
    vector<vector<complex<double>>> getU_hat() const;
};
class VMDDecomposition
{
private:
    vector<vector<complex<double>>> u_hat;//моды во частотной области
    vector<vector<complex<double>>> u;//моды во временной области
    vector<vector<double>> omega;//центральные частоты
public:
    VMDDecomposition(vector<vector<complex<double>>> u_hat, vector<vector<complex<double>>> u, vector<vector<double>> omega);
    vector<vector<complex<double>>> GetUhat() const;
    vector<vector<complex<double>>> GetU() const;
    vector<vector<double>> GetOmega() const;
};

template<typename SignalDT> 
VMDDecomposition VMD(std::vector<SignalDT> f, double alpha_var, double tau, int k, bool DC, int init, double tol);
