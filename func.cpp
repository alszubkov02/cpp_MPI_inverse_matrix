#include "func.h"
#include <immintrin.h>


double get_cpu_time(){
    struct rusage buf;
    getrusage(RUSAGE_THREAD, &buf);
    return buf.ru_utime.tv_sec + buf.ru_utime.tv_usec * 1e-6;
}

int local_to_global(int m, int p, int k, int i_loc){
    int i_loc_m = i_loc / m;
    int i_glob_m = i_loc_m * p + k;
    return i_glob_m * m + i_loc % m;
}

int global_to_local(int m, int p, int i_glob){
    int i_glob_m = i_glob / m;
    int i_loc_m = i_glob_m / p ;
    return i_loc_m * m + i_glob % m;
}

int my_max(int i, int j){
    return i > j ? i : j;
}

double func(int s, int n, int i, int j){
    if(s == 1) return n - my_max(i, j);
    else if(s == 2) return my_max(i, j) + 1;
    else if(s == 3) return std::fabs(i-j);
    else return 1. / (i + j + 1);
}

int get_rows(int n, int m, int p, int k){
    int b = (n + m - 1) / m;
    return (k >= b%p ? b/p : (b + p - 1)/p);
}

int get_max_rows(int n, int m, int p){
    int b = (n + m - 1)/m;
    return (b + p - 1)/p;
}


void init_a(double* a, int n, int m, int p, int k, int s){
    int rows = get_rows(n, m, p, k);
    for(int i_loc_m = 0; i_loc_m < rows; ++i_loc_m){
        for(int i_loc = i_loc_m * m; i_loc < (i_loc_m + 1) * m; ++i_loc){
            int i_glob = local_to_global(m, p, k, i_loc);
            if(i_glob >= n) break;
            for(int j = 0; j < n; ++j){
                a[n*i_loc + j] = func(s, n, i_glob, j);
            }
        }
    }
}

void init_x(double* x, int n, int m, int p, int k){
    int rows = get_rows(n, m, p, k);
    for(int i_loc_m = 0; i_loc_m < rows; ++i_loc_m){
        for(int i_loc = i_loc_m * m; i_loc < (i_loc_m + 1) * m; ++i_loc){
            int i_glob = local_to_global(m, p, k, i_loc);
            if(i_glob >= n) break;
            for(int j = 0; j < n; ++j){
                x[n*i_loc + j] = (i_glob == j) ? 1.0 : 0.0;
            }
        }
    }
}



int read_array(FILE* fp, double* buf, int len){
    for(int i = 0; i < len; ++i){
        if(fscanf(fp, "%lf", buf + i) != 1){
            return -1;
        }
    }
    return 0;
}

int read_a(double* a, int n, int m, int p, int k, const char* name, double* buf, MPI_Comm com){

    int main_k = 0;
    FILE* fp = nullptr;
    int err = 0;
    if(k == main_k){
        fp = fopen(name, "r");
        if(fp == nullptr) err = 1;
    }
    MPI_Bcast(&err, 1, MPI_INT, main_k, com);
    if(err) return err;
    memset(buf, 0, n*m*sizeof(double));

    int max_b = (n + m - 1) / m;
    for(int b = 0; b < max_b; ++b){
        int owner = b % p;
        int num_rows = ((b + 1) * m  <= n ? m : n - b * m);
        int b_loc = b / p;
        if(k == main_k){
            err += read_array(fp, buf, n*num_rows);
            if(owner == k) memcpy(a + b_loc*m*n, buf, n*num_rows*sizeof(double)); 
            else MPI_Send(buf, n*num_rows, MPI_DOUBLE, owner, 0 , com);
        } 
        else if(k == owner){
            MPI_Status st;
            MPI_Recv(a + b_loc*m*n, n*num_rows, MPI_DOUBLE, main_k, 0 , com, &st);
        }
    }

    if(k == main_k) fclose(fp);

    MPI_Bcast(&err, 1, MPI_INT, main_k, com);
    if(err) return err;
    return 0;
}




int print_array(double* a, int n, int m, int printed_rows, int max_print){
    if(printed_rows >= max_print) return 0;

    int p_n = (n > max_print ? max_print : n);
    int p_m = printed_rows + m < max_print ? m : max_print - printed_rows;
    for(int i = 0; i < p_m; ++i){
        for(int j = 0; j < p_n; ++j){
            printf(" %10.3e", a[i*n + j]);
        }
        printf("\n");
    }
    return p_m;
}

void print_a(double* a, int n, int m, int p, int k, double* buf, int max_print, MPI_Comm com) {
    int main_k = 0;
    int max_b = (n + m - 1) / m;
    int printed_rows = 0;
    for(int b = 0; b < max_b; ++b){
        int owner = b % p;
        int num_rows = (b + 1) * m <= n ? m : n - b * m; 
        int b_loc = b / p;
        if(k == main_k){
            if(owner == main_k) printed_rows += print_array(a + b_loc*m*n, n, num_rows, printed_rows, max_print);
            else{
                MPI_Status st;
                MPI_Recv(buf, n*num_rows, MPI_DOUBLE, owner, 0/*tag*/, com, &st);
                printed_rows += print_array(buf, n, num_rows, printed_rows, max_print);
            }
        } 
        else if(owner == k){
            MPI_Send(a + b_loc*m*n, n*num_rows, MPI_DOUBLE, main_k, 0/*tag*/, com);
        }
    }
}

double matrix_norm(double* a, int n, int m, int p, int k, MPI_Comm com) {
    int rows = get_rows(n, m, p, k);
    
    double norm = -1;
    double global_sum;
    double sum;
    for (int j = 0; j < n; ++j) {
        sum = 0;
        for (int i_loc_m = 0; i_loc_m < rows; ++i_loc_m) {
            for (int i_loc = i_loc_m * m; i_loc < (i_loc_m + 1) * m; ++i_loc) {
                int i_glob = local_to_global(m, p, k, i_loc);
                if (i_glob >= n) {
                    break;
                }

                sum += std::fabs(a[n*i_loc + j]);
            }
        }

        MPI_Allreduce(&sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, com);
        norm = (norm > global_sum) ? norm : global_sum;
    }

    return norm;
}


void add_matrix(int n, int m, double* a, double* b){
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < m; ++j){
            a[m*i + j] += b[m*i + j];           
        }
    }
}


double residual(double* a, double* x, int n, int m, int p, int k, MPI_Comm com, double* block1, double* block2, double* block3, double* buffer, double* ans){
    int f = n / m;
    int l = n - f * m;
    int h = l ? f + 1 : f;
    int rows = get_rows(n, m, p, k);
    double answer = 0.0;
    int step = 0;
    for(int i_glob_m = 0; i_glob_m < f; ++i_glob_m){

        int i_loc_m = i_glob_m / p;
        int main_k = i_glob_m % p;

        if(k == main_k){
            memcpy(buffer, a + m*n*i_loc_m, n*m*sizeof(double));
        }
        memset(ans, 0.0, n*m*sizeof(double));
        MPI_Bcast(buffer, n*m, MPI_DOUBLE, main_k, com);
        for(int j = 0; j < h; ++j){
            for(int i = 0; i < rows; ++i){
                int multiplier_rows = (m + local_to_global(m, p, k, i*m) <= n ? m : l);
                int block_cols = j < f ? m : l;
                
                get_block(buffer, block1, n, m, f, l, 0, i*p+k);
                get_block(x,      block2, n, m, f, l, i, j);
                multiply(block1, block2, block3, m, multiplier_rows, multiplier_rows, block_cols);
                get_block(ans, block2, n, m, f, l, 0, j);
                add_matrix(m, block_cols, block2, block3);
                put_block(ans, block2, n, m, f, l, 0, j);
            }
        } 
        MPI_Allreduce(ans, buffer, m*n, MPI_DOUBLE, MPI_SUM, com);
        for(int i = 0; i < m; ++i){
            double sum = 0.0;
            for(int j = 0; j < n; ++j){
                if(j == i + step) sum += fabs(buffer[i*n + j]-1.0);
                else sum += fabs(buffer[i*n + j]);
            }
            answer = std::max(answer, sum);
        }
        step += m;
    }
    return answer;
}