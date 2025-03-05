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
#include <concepts>
#include <algorithm>

using namespace std;

//������������� ��������� ��������: Expert C++
// First published: April 2020
//Second edition : August 2023
//Production reference : 1280723
//Published by Packt Publishing Ltd.
//Grosvenor House
//11 St Paul�s Square
//Birmingham
//B3 1RB, UK.
//ISBN 978 - 1 - 80461 - 783 - 0
// ������:
//Marcelo Guerra Hahn
//Araks Tigranyan
//John Asatryan
//Vardan Grigoryan
//Shunguang Wu

template <typename SignalDT>
concept d_types = requires(SignalDT a, SignalDT b) {
	{ a == b } -> std::convertible_to<double>;
	{ a != b } -> std::convertible_to<double>;
};

//����� ��������� �� �������� ����� �������
//SignalDT (��� ��������� �������) - int, unsigned, double, complex
template <d_types SignalDT>
vector<SignalDT> fftshift(vector<SignalDT> f_hat);

//������� �������������� �����
template <d_types SignalDT>
std::vector<std::complex<SignalDT>> fft(vector<SignalDT> f);

//�������� �������������� �����
template <d_types SignalDT> 
std::vector<std::complex<SignalDT>> ifft(vector<SignalDT> f);

//��������� ��������� ���������� (�� start �� end � �������� ������ ����� num)
template <d_types T>
inline std::vector<T> linspace(int start, int end, int num);

//��������� ��������� ���������� (�� start �� end � ����� step)
template <d_types T>
inline std::vector<T> arange(T start, T end, T step);

template <d_types T>
inline std::vector<T> dot_product(vector<T> a, vector<T> b);

template <d_types T>
inline std::vector<T> abs(T a);

template <d_types T>
inline T sum(T a);

template <d_types T>
inline std::vector<T> sort(T a);

templat e<d_types T>
inline std::vector<T> exp(T a);

template <d_types T>
inline std::vector<T> random(int num);

