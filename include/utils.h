#ifndef UTILS_H
#define UTILS_H


#include<iostream>
#include <string>
#include <tuple>
#include <vector>


int calc_n_elem(const std::string &path);

std::string define_filepath(std::string &base);

template <class T>
T** alloc2d(const int ni, const int nj);

template <class T>
T*** alloc3d(const int ni, const int nj, const int nk);

template <class T>
T**** alloc4d(const int ni, const int nj, const int nk, const int n);

template <class T>
void free2d(T ***arr, const int ni, const int nj);

template <class T>
void free3d(T ****arr, const int ni, const int nj, const int nk);

template <class T>
void free4d(T *****arr, const int ni, const int nj, const int nk, const int n);

template <class T>
void write_array(T *arr, const int n,  const std::string &path);

template <class T>
std::tuple< T*,double> linspace(T xi, T xf, int n);
#endif
