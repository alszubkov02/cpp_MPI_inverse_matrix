#include "func.h"

int main(int argc, char* argv[]){
    int p, k;
    MPI_Comm com = MPI_COMM_WORLD;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(com, &p);
    MPI_Comm_rank(com, &k);
    solution(argc, argv, com, p, k); 
    MPI_Finalize();
    return 0;
}