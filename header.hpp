#pragma once
#define _USE_MATH_DEFINES
#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
#include <functional>
#include <algorithm>
#include <array>
#include <string_view>
#include <random>
#include <numeric>

using namespace std;

//»спользовалс€ следующий источник: Expert C++
// First published: April 2020
//Second edition : August 2023
//Production reference : 1280723
//Published by Packt Publishing Ltd.
//Grosvenor House
//11 St PaulТs Square
//Birmingham
//B3 1RB, UK.
//ISBN 978 - 1 - 80461 - 783 - 0
// јвторы:
//Marcelo Guerra Hahn
//Araks Tigranyan
//John Asatryan
//Vardan Grigoryan
//Shunguang Wu


//сдвиг элементов на половину длины массива
//SignalDT (тип исходного сигнала) - int, unsigned, double, complex
template<typename SignalDT>
vector<SignalDT> fftshift(vector<SignalDT> f_hat);

//быстрое преобразование фурье
template<typename SignalDT> 
std::vector<std::complex<double>> fft(vector<SignalDT> f);

//обратное преобразование фурье
template<typename SignalDT> 
std::vector<std::complex<double>> ifft(vector<SignalDT> f);

//генераци€ временных интервалов (от start до end с заданным числом шагов num)
inline std::vector<double> linspace(int start, int end, int num);

//генераци€ временных интервалов (от start до end с шагом step)
template<typename T>
inline std::vector<double> arange(T start, T end, T step);

template<typename T>
inline std::vector<T> dot_product(vector<T> a, vector<T> b);

template<typename T>
inline std::vector<double> abs(T a);

template<typename T>
inline double sum(T a);

template<typename T>
inline std::vector<double> sort(T a);

template<typename T>
inline std::vector<double> exp(T a);

inline std::vector<double> random(int num);

