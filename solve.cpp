#include "func.h"

void get_block(double* a, double* b, int n, int m, int f, int l, int row, int col){
    int h = (row < f ? m : l);
    int w = (col < f ? m : l);
    double* block = a + m * (row * n + col);
    for(int i = 0; i < h; i ++){
        for(int j = 0; j < w; j ++){
            b[i*w+j] = block[i*n+j];
        }
    }
}

void put_block(double* a, double* b, int n, int m, int f, int l, int row, int col){
    int h = (row < f ? m : l);
    int w = (col < f ? m : l);
    double* block = a + m * (row * n + col);
    for(int i = 0; i < h; i ++){
        for(int j = 0; j < w; j ++){
            block[i * n + j] = b[i * w + j];
        }
    }
}


void put_vector(int i, int m, int k, int l, double* b_i, double* b) {
    int length = i < k ? m : l;
    for (int p = 0; p < length; ++p) {
        b[m*i + p] = b_i[p];
    }
}

void subtract_matrix(int n, int m, double* a, double* b) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            a[m*i + j] -= b[m*i + j];           
        }
    }
}

void my_swap(double& a, double& b){
    double tmp = a;
    a = b;
    b = tmp;
}
void swap_rows(double* a, int n, int i, int j){
    int k;
    for(k = 0; k < n; k++){
        my_swap(a[i*n+k],a[j*n+k]);
    }
}

int reverse(double* C, double* D, int n,double norm){
    double leader;
    int ind_lead;
    double tmp;
    int i,j,k;
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            D[i*n+j] = (i == j) ? 1.0 : 0.0;
        }
    }

    for (k = 0; k < n; k++) {
        leader = C[k*n+k];
        ind_lead = k;
        for(j = k; j < n; j++){
            if(fabs(C[j*n+k]) > fabs(leader)){
                leader = C[j*n+k];
                ind_lead = j;
            }
        }
        if (fabs(leader) < EPS*norm)  return -1;
        swap_rows(C,n,k,ind_lead);
        swap_rows(D,n,k,ind_lead);

        tmp = C[k*n+k];

        for (j = k; j < n; j++) C[k*n+j] /= tmp;
        for (j = 0; j < n; j++) D[k*n+j] /= tmp;

        for (i = 0;i < n; i++){
            if(i==k) continue;
            tmp = C[i*n+k];
            for(j = k; j < n; j++) C[i*n+j] -= tmp * C[k*n+j];
            for(j = 0; j < n; j++) D[i*n+j] -= tmp * D[k*n+j];
        }
    }

    return 0;
}

void multiply(double* A, double* B,double* C, int na, int ma, int nb, int mb){
    int v3 = na%3;
    int h3 = mb%3;
    int i = 0, j = 0, k = 0;
    for(i = 0; i < na; i++){
        for(j = 0; j < mb; j++){
            C[i*mb+j] = 0.0;
        }
    }
    double s, tmp;
    double s00, s01, s02,
           s10, s11, s12,
           s20, s21, s22;
    for(i = 0; i < v3; i++){
        for(j = 0; j < h3; j++){
            s = 0;
            for(k = 0; k < nb; k++){
                s += A[i*ma+k]*B[k*mb+j];
            }
            C[i*mb+j] = s;
        }
        for(;j < mb;j+=3){
            s00 = 0;
            s01 = 0;
            s02 = 0;
            for(k = 0; k < ma; k++){
                tmp = A[i*ma+k];
                s00 += tmp*B[k*mb+j];
                s01 += tmp*B[k*mb+j+1];
                s02 += tmp*B[k*mb+j+2];
            }
            C[i*mb+j] += s00;
            C[i*mb+j+1] += s01;
            C[i*mb+j+2] += s02;
        }
    }
    for(;i < na; i+=3){
        for(j = 0; j < h3; j++){
            s00 = 0;
            s10 = 0;
            s20 = 0;
            for(k = 0; k < ma; k++){
                tmp = B[k*mb+j];
                s00 += A[i*ma+k]*tmp;
                s10 += A[(i+1)*ma+k]*tmp;
                s20 += A[(i+2)*ma+k]*tmp;
            }
            C[i*mb+j] += s00;
            C[(i+1)*mb+j] += s10;
            C[(i+2)*mb+j] += s20;
        }
        for(;j < mb; j+=3){
            s00 = 0;
            s01 = 0;
            s02 = 0;
            s10 = 0;
            s11 = 0;
            s12 = 0;
            s20 = 0;
            s21 = 0;
            s22 = 0;
            for(k = 0; k < ma; k++){
                s00 += A[i*ma+k]*B[k*mb+j];
                s01 += A[i*ma+k]*B[k*mb+j+1];
                s02 += A[i*ma+k]*B[k*mb+j+2];
                s10 += A[(i+1)*ma+k]*B[k*mb+j];
                s11 += A[(i+1)*ma+k]*B[k*mb+j+1];
                s12 += A[(i+1)*ma+k]*B[k*mb+j+2];
                s20 += A[(i+2)*ma+k]*B[k*mb+j];
                s21 += A[(i+2)*ma+k]*B[k*mb+j+1];
                s22 += A[(i+2)*ma+k]*B[k*mb+j+2];
            }
            C[i*mb+j] += s00;
            C[i*mb+j+1] += s01;
            C[i*mb+j+2] += s02;
            C[(i+1)*mb+j] += s10;
            C[(i+1)*mb+j+1] += s11;
            C[(i+1)*mb+j+2] += s12;
            C[(i+2)*mb+j] += s20;
            C[(i+2)*mb+j+1] += s21;
            C[(i+2)*mb+j+2] += s22;
        }
    }
}



void free_all(double* a, double* x, double* buf_a, double* buf_x, double* block1, double* block2, double* block3){
    delete[] a;
    delete[] x;
    delete[] buf_a;
    delete[] buf_x;
    delete[] block1;
    delete[] block2;
    delete[] block3;
}

int get_max_cols(int n, int m, int p){
    int b = (n + m - 1) / m;
    return (b + p - 1) / p * m;
}


void solution(int argc, char* argv[], MPI_Comm com, int p, int k){
    const int task = 20;
    const char* file_name = nullptr;
    int n, m, r, s;
    double r1 = 0.0, r2 = 0.0, t1 = 0.0, t2 = 0.0;

    if(!(argc == 5 || argc == 6) || sscanf(argv[1], "%d", &n) != 1 || sscanf(argv[2], "%d", &m) != 1
     	|| sscanf(argv[3], "%d", &r) != 1 || sscanf(argv[4], "%d", &s) != 1){
        printf("Usage: %s n m r s (file)\n", argv[0]);
    	return;
	}

    int rows = get_rows(n, m, p, k);
    double* a = new double[rows*m*n];
    double* buf_a = new double[m*n];

    if(s) init_a(a, n, m, p, k, s);
    else{
        file_name = argv[5];
        if(read_a(a, n, m, p, k, file_name, buf_a, com)){    
            if(k == 0){
                printf("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n",
                argv[0], task, -1., -1., 0., 0., s, n, m, p);
            }
            delete[] a;
            delete[] buf_a;
            return;
        }
    }

    if(k == 0){
        printf("Matrix :\n");
    }
    print_a(a, n, m, p, k, buf_a, r, com);

    double* x = new double[rows*m*n];
    double* buf_x = new double[m*n];
    double* block1 = new double[m*m];
    double* block2 = new double[m*m];
    double* block3 = new double[m*m];

    init_x(x, n, m, p, k);

    double a_norm = matrix_norm(a, n, m, p, k, com);
    t1 = get_cpu_time();

    int f = n / m;
    int l = n - f * m;
    int h = l ? f + 1 : f;

    int err = 0;
    for(int i_glob_m = 0; i_glob_m < f; ++i_glob_m){

        int i_loc_m = i_glob_m / p;
        int main_k = i_glob_m % p;

        if(k == main_k){
            get_block(a, block1, n, m, f, l, i_loc_m, i_glob_m);
            if(reverse(block1, block2, m, a_norm) != 0){
                err = -1;
            }           
            for(int j_glob_m = i_glob_m + 1; j_glob_m < f; ++j_glob_m){
                get_block(a, block1, n, m, f, l, i_loc_m, j_glob_m);
                multiply(block2, block1, block3, m, m, m, m);
                put_block(a, block3, n, m, f, l, i_loc_m, j_glob_m);
            }

            for(int j_glob_m = 0; j_glob_m < f; ++j_glob_m){
                get_block(x, block1, n, m, f, l, i_loc_m, j_glob_m);
                multiply(block2, block1, block3, m, m, m, m);
                put_block(x, block3, n, m, f, l, i_loc_m, j_glob_m);
            }
            
            if(l){
                get_block(a, block1, n, m, f, l, i_loc_m, f);
                multiply(block2, block1, block3, m, m, m, l);
                put_block(a, block3, n, m, f, l, i_loc_m, f);

                get_block(x, block1, n, m, f, l, i_loc_m, f);
                multiply(block2, block1, block3, m, m, m, l);
                put_block(x, block3, n, m, f, l, i_loc_m, f);
            }

            memcpy(buf_a, a + m*n*i_loc_m, n*m*sizeof(double));
            memcpy(buf_x, x + m*n*i_loc_m, n*m*sizeof(double)); 
        } 

        MPI_Bcast(buf_a, n*m, MPI_DOUBLE, main_k, com);
        MPI_Bcast(buf_x, n*m, MPI_DOUBLE, main_k, com); 
        MPI_Bcast(&err, 1, MPI_INT, main_k, com);

        if(err){
            if(k == 0){
                printf("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n",
                argv[0], task, -1., -1., t1, 0., s, n, m, p);
            }
            free_all(a, x, buf_a, buf_x, block1, block2, block3);
            return;
        }
        
        int offset = (main_k <= k) ? k - main_k : p + k - main_k; 
        int loc_start_m = (i_glob_m + offset) / p;
        if(k == main_k) loc_start_m++;

        for(int i = 0; i < rows; ++i){
            if(k == main_k && i == i_loc_m) continue;
            get_block(a, block1, n, m, f, l, i, i_glob_m);
            int multiplier_rows = (m + local_to_global(m, p, k, i * m) <= n ? m : l);

            for(int j = i_glob_m + 1; j < h; ++j){
                int block_cols = j < f ? m : l;
                get_block(buf_a, block2, n, m, f, l, 0, j);
                multiply(block1, block2, block3, multiplier_rows, m, m, block_cols);
                get_block(a, block2, n, m, f, l, i, j);
                subtract_matrix(multiplier_rows, block_cols, block2, block3);
                put_block(a, block2, n, m, f, l, i, j);
            }

            for(int j = 0; j < h; ++j){
                int block_cols = j < f ? m : l;
                get_block(buf_x, block2, n, m, f, l, 0, j);
                multiply(block1, block2, block3, multiplier_rows, m, m, block_cols);
                get_block(x, block2, n, m, f, l, i, j);
                subtract_matrix(multiplier_rows, block_cols, block2, block3);
                put_block(x, block2, n, m, f, l, i, j);
            }
        }
    }
    if(l){
        int var = f/p; 
        if(f == k + var*p){
            get_block(a, block1, n, m, f, l, var, f);  
            if(reverse(block1, block2, l, a_norm) != 0){
                err = -1;
            }
            for(int j = 0; j < f; ++j){
                get_block(x, block1, n, m, f, l, var, j);
                multiply(block2, block1, block3, l, l, l, m);
                put_block(x, block3, n, m, f, l, var, j);
            }
            get_block(x, block1, n, m, f, l, var, f);
            multiply(block2, block1, block3, l, l, l, l);
            put_block(x, block3, n, m, f, l, var, f);

            memcpy(buf_x, x + m*n*var, n*m*sizeof(double)); 
        }

        MPI_Bcast(&err, 1, MPI_INT, f - var*p, com);

        if(err){
            if(k == 0){
                printf("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n",
                argv[0], task, -1., -1., t1, 0., s, n, m, p);
            }

            free_all(a, x, buf_a, buf_x, block1, block2, block3);
            return;
        }

        MPI_Bcast(buf_x, n*m, MPI_DOUBLE, f - var*p, com);

        for(int i = 0; i < rows; ++i){
            if(k == f - var*p && i == var) continue;
            get_block(a, block1, n, m, f, l, i, f);
            //int multiplier_rows = (m + local_to_global(m, p, k, i * m) <= n ? m : l);


            for(int j = 0; j < h; ++j){
                int block_cols = j < f ? m : l;
                get_block(buf_x, block2, n, m, f, l, 0, j);
                multiply(block1, block2, block3, m, l, l, block_cols);
                get_block(x, block2, n, m, f, l, i, j);
                subtract_matrix(m, block_cols, block2, block3);
                put_block(x, block2, n, m, f, l, i, j);
            }
        }
    }

    t1 = get_cpu_time() - t1;
    MPI_Barrier(com);

    if(n <= 11000){
        if(s) init_a(a, n, m, p, k, s);
        else{
            if(read_a(a, n, m, p, k, file_name, buf_a, com)){       
                if(k == 0){
                    printf("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n",
                    argv[0], task, -1., -1., t1, 0., s, n, m, p);
                }
                free_all(a, x, buf_a, buf_x, block1, block2, block3);
                return;
            }
        }
        MPI_Barrier(com);
        t2 = get_cpu_time();
        r1 = residual(a, x, n, m, p, k, com, block1, block2, block3, buf_a, buf_x);
        r2 = residual(x, a, n, m, p, k, com, block1, block2, block3, buf_a, buf_x);
        t2 = get_cpu_time() - t2;
    }

    if(k == 0){
        printf("Inverse matrix :\n");
    }
    print_a(x, n, m, p, k, buf_x, r, com);


    if(k == 0){
        printf("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n", argv[0], task, r1, r2, t1, t2, s, n, m, p);
    }

    free_all(a, x, buf_a, buf_x, block1, block2, block3);
}
