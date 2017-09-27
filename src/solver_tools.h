#ifndef SOLVER_TOOLS_H
#define SOLVER_TOOLS_H

#include "global_var.h"
#include"general_tools.h"

void construct_tridiag_example(const int n, double *&a, double *&b, double *&c
                               , double *&d);

void test_cyclic_tridiag(const int& n);

void tridiag_solve_mh(const int n, const double *a, const double *b
                     , const double *c, const double *d, double *x);

void cyclic_tridiag_solve_mh(const int n, const double *a, const double *b
                     , const double *c, const double *d, double *x);

void test_tridaig(const int& n);

#endif
