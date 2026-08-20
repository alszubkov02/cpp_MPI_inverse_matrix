#include "mpi.h"
#include <stdio.h>
#include <string.h>
#include <sys/resource.h>
#include <cmath>

#define EPS 1e-14

void get_block(double* a, double* b, int n, int m, int f, int l, int row, int col);
void put_block(double* a, double* b, int n, int m, int f, int l, int row, int col);

double get_cpu_time();
int local_to_global(int m, int p, int k, int i_loc);

int get_rows(int n, int m, int p, int k);
int get_max_rows(int n, int m, int p);
void init_a(double* a, int n, int m, int p, int k, int s);
void init_x(double* x, int n, int m, int p, int k);

int read_a(double* a, int n, int m, int p, int k, const char* name, double* buf, MPI_Comm com);
void print_a(double* a, int n, int m, int p, int k, double* buf, int max_print, MPI_Comm com);

double matrix_norm(double* a,int n, int m, int p, int k, MPI_Comm com);
void multiply(double* A, double* B,double* C, int na, int ma,int nb, int mb);
void solution(int argc, char* argv[], MPI_Comm com, int p, int k);
double residual(double* a, double* x, int n, int m, int p, int k, MPI_Comm com, double* block1, double* block2, double* block3, double* buffer, double* ans);