#ifndef UTILS_H
#define UTILS_H


#include<iostream>
#include <string>
#include <tuple>
#include <vector>


int calc_n_elem(const std::string &path);

std::string define_filepath(std::string &base);

template <class T>
void write_array(T *arr, const int n,  const std::string &path);

template <class T>
std::tuple< T*,double> linspace(T xi, T xf, int n);
#endif
