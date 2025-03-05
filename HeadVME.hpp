#include "header.hpp"
//F<complex<double>>& f_hat_onesided, Ud& u_hat_d, OmegaAxis& omega_axis, Omega& omega_d, Lamb& lamb

template<d_types SignalDT>
class F;           //предобработанный временной ряд

class Ud;          //декомпозиция (внутримодовая функция)
class OmegaAxis;   //дискретная частотная область
class NumbIter;    //ограничение количества итераций
class Epsilon;     //условие остановки итераций (максимальное значение относительной разности нового(после итерации) и старого значения)
class Omega;       //омега-значение для декомпозиции
class Lamb;        //лямбда-значение
class VMEDecomposition; //вывод декомпозиции


template<d_types SignalDT>
class F
{
private:
    vector<complex<double>> f_hat_onesided;
    int T;
public:
    //иницализация и предобратка сигнала, получение значение средней частоты (fs) и длины сигнала Т
    F(vector<SignalDT> signal);

    vector<complex<double>> operator()() const;
    complex<double> operator[](int i) const;
    F<SignalDT>& operator=(vector<complex<double>> value) = delete;
    vector<complex<double>> getF() const;
    size_t size() const;
    int getT() const;
    ~F() = default;
};

class Epsilon
{
private:
    complex<double> udiff;
    complex<double> count_udiff(complex<double> udiff, int T, int n, Ud& u_hat_d);
public:
    //double eps = numeric_limits<double>::epsilon();
    //double udiff = tol + eps;
    //инициализация условия остановки
    Epsilon(double tol);

    complex<double> getUdiff() const;
    //обновление условия остановки после очередной итерации
    void update_udiff(complex<double> udiff, int T, int n, Ud& u_hat_d);
    ~Epsilon() = default;
};
class OmegaAxis
{
private:
    vector<double> f;
public:
    //генерация дискретной области частот
    OmegaAxis(int T);

    vector<double> operator()() const;
    double operator[](int i) const;
    OmegaAxis& operator=(vector<double> value) = delete;
    size_t size() const;
    ~OmegaAxis() = default;
};
class NumbIter
{
private:
    //int N = 500;
    int N = 500;
public:
    NumbIter();

    int operator()() const;
    NumbIter& operator += (const int& counter);
    void SetMinToN(int n);
    ~NumbIter() = default;
};
class Lamb
{
private:
    vector<vector<complex<double>>> lamb;
    vector<complex<double>> count_lamb(int T, int n, double tau, double Alpha, vector<complex<double>>& f_hat_onesided, Ud& u_hat_d, OmegaAxis& omega_axis, Omega& omega_d);
public:
    //инициализация лямбда
    Lamb(int N = 500, int T);
    //обновление лямбда
    void update_lamb(int T, int n, double tau, double Alpha, vector<complex<double>>& f_hat_onesided, Ud& u_hat_d, OmegaAxis& omega_axis, Omega& omega_d);
    vector<vector<complex<double>>> operator()() const;
    vector<complex<double>>& operator[](int i);
    size_t size_1() const;
    size_t size_2() const;
    ~Lamb() = default;
};
class Omega
{
private:
    vector<complex<double>> omega_d;//центральные частоты после каждой итерации
    bool check_finite = true;
    complex<double> count_omegad(int T, int n, Ud& u_hat_d, OmegaAxis& omega_axis);
public:
    //инициализация значения омега для каждой моды
    Omega(int N = 500, double omega_int, double fs);
    //проверка переполнения
    bool isFinite();
    vector<complex<double>> operator()() const;
    complex<double>& operator[](int i);
    //обновления значений омега
    void update_omega(int T, int n, Ud& u_hat_d, OmegaAxis& omega_axis);
    //финальное обновление значение омега
    void minOmega(int n);
    ~Omega() = default;
};

template<d_types SignalDT>
class Ud
{
private:
    vector<vector<complex<SignalDT>>> u_hat_d; //значения модов после каждой итерации
    bool check_finite = true;
    vector<complex<SignalDT>> u_hatd_final; //моды во частотной области
    vector<complex<SignalDT>> u_d_final; //моды во временной области
    vector<complex<SignalDT>> count_uhat(int T, int n, double Alpha, vector<complex<SignalDT>>& f_hat_onesided, OmegaAxis& omega_axis, Omega& omega_d, Lamb& lamb);
public:
    Ud(int N = 500, int T);
    //обновление моды
    void update_u_hat(int T, int n, double Alpha, vector<complex<SignalDT>>& f_hat_onesided, OmegaAxis& omega_axis, Omega& omega_d, Lamb& lamb);
    //проверка переполнения
    bool isFinite();
    vector<vector<complex<SignalDT>>> operator()() const;
    vector<complex<SignalDT>>& operator[](int i);
    size_t size_1() const;
    size_t size_2() const;
    //получение итоговой декомпозиции
    void ud(int T, int N);

    vector<complex<SignalDT>> getU_D() const;
    vector<complex<SignalDT>> getU_hatD() const;
    ~Ud() = default;
};

template<d_types SignalDT>
class VMEDecomposition
{
private:
    vector<complex<SignalDT>> u_hatd; //моды во частотной области
    vector<complex<SignalDT>> u_d; //моды во временной области
    vector<complex<SignalDT>> omega_d; //центральные частоты
public:
    VMEDecomposition(vector<complex<SignalDT>> u_hatd, vector<complex<SignalDT>> u_d, vector<complex<SignalDT>> omega_d);
    vector<complex<SignalDT>> GetUhatD() const;
    vector<complex<SignalDT>> GetUD() const;
    vector<complex<SignalDT>> GetOmegaD() const;
    ~VMEDecomposition() = default;
};
