#ifndef NUM_ITERATIONS
#define NUM_ITERATIONS 2
#endif

#define LOG(rank, str) fprintf(stderr, "#R%d: %s\n", rank, str)

#define SPEC_RESTRICT __restrict__
// #define SPEC_RESTRICT restrict

#include <mpi.h>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <string>
#include <cmath>
#include <unistd.h>
#include <sys/syscall.h>
#include "math.h"
#include "extrae_include.h"

static int my_rank_id   = 0;
static int num_procs    = -1;
int numberOfTasks       = 0;
int matrixSize          = 100;

void initialize_matrix(double *mat, int matrixSize, double val) {
	for(int i=0; i<matrixSize*matrixSize; i++) {
		mat[i]= val;
	}
}

void compute_matrix_matrix(double * SPEC_RESTRICT a, double * SPEC_RESTRICT b, double * SPEC_RESTRICT c, int matrixSize) {
    EXTRAE_ENTER(EVENT_COMPUTE);
    for(int i=0;i<matrixSize;i++) {
        for(int j=0;j<matrixSize;j++) {
            c[i*matrixSize+j]=0;
            for(int k=0;k<matrixSize;k++) {
                c[i*matrixSize+j] += a[i*matrixSize+k] * b[k*matrixSize+j];
            }
        }
    }
    EXTRAE_EXIT(EVENT_COMPUTE);
}

bool check_test_matrix(double *c, double val, int matrixSize) {
	for(int i=0;i<matrixSize;i++) {
		for(int j=0;j<matrixSize;j++) {
			if(fabs(c[i*matrixSize+j] - val) > 1e-3) {
				printf("#R%d (OS_TID:%ld): Error in matrix entry (%d,%d) expected:%f but value is %f\n", my_rank_id, syscall(SYS_gettid),i,j,val,c[i*matrixSize+j]);
				return false;
			}
		}
	}
	return true;
}

void printHelpMessage() {
    if(my_rank_id == 0) {
        std::cout << "Usage: mpiexec -n np ./matrixExample matrixSize nt_(0) ... nt_(np-1) " << std::endl;
        std::cout << "    Arguments: " << std::endl;
        std::cout << "        matrixSize:   Number of elements of the matrixSize x matrixSize matrices" << std::endl;
        std::cout << "        nt_(i):       Number of tasks for process i " << std::endl;
    }
}

int parse_command_line_args(int argc, char **argv) {
    if(argc==2) {
        matrixSize      = atoi( argv[1] );
        numberOfTasks   = 200;
    } else if(argc==num_procs+2) {
        matrixSize      = atoi( argv[1] ); 
        numberOfTasks   = atoi( argv[my_rank_id+2] ); 
    } else { 
        printHelpMessage();
        return 1;
    }
    return 0;
}

int main(int argc, char **argv)
{
	int iMyRank, iNumProcs;
	int provided;
	int requested = MPI_THREAD_MULTIPLE;
	MPI_Init_thread(&argc, &argv, requested, &provided);
	MPI_Comm_size(MPI_COMM_WORLD, &iNumProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &iMyRank);
    my_rank_id = iMyRank;
    num_procs = iNumProcs;

	double fTimeStart, fTimeEnd;
	bool pass = true;

#ifdef TRACE_EXTRAE
    Extrae_init();
    REGISTER_EXTRAE();
#endif

    int ret_code = parse_command_line_args(argc, argv);
    if (ret_code != 0) {
        return ret_code;
    }
	
    std::string msg = "will create "+std::to_string(numberOfTasks)+" tasks";
    LOG(iMyRank, msg.c_str());

	double **matrices_a, **matrices_b, **matrices_c;
	matrices_a = new double*[numberOfTasks];
	matrices_b = new double*[numberOfTasks];
	matrices_c = new double*[numberOfTasks];

    #pragma omp parallel for
	for(int i=0; i<numberOfTasks; i++) {
        int cur_size = matrixSize;

        EXTRAE_ENTER(EVENT_ALLOCATE);
 		matrices_a[i] = new double[(long)cur_size*cur_size];
    	matrices_b[i] = new double[(long)cur_size*cur_size];
    	matrices_c[i] = new double[(long)cur_size*cur_size];
        EXTRAE_EXIT(EVENT_ALLOCATE);
        
        EXTRAE_ENTER(EVENT_INIT);
        initialize_matrix(matrices_a[i], cur_size, 1);
        initialize_matrix(matrices_b[i], cur_size, 1);
        initialize_matrix(matrices_c[i], cur_size, 0);
        EXTRAE_EXIT(EVENT_INIT);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    fTimeStart=MPI_Wtime();
    #pragma omp parallel
    {
        EXTRAE_ENTER(EVENT_SIMULATION);
        for(int iter = 0; iter < NUM_ITERATIONS; iter++) {
		
            #pragma omp for
            for(int i=0; i<numberOfTasks; i++) {
                int cur_size = matrixSize;

                EXTRAE_ENTER(EVENT_CREATE);
                #pragma omp task default(shared) firstprivate(i)
                {
                    compute_matrix_matrix(matrices_a[i], matrices_b[i], matrices_c[i], cur_size);
                }
                EXTRAE_EXIT(EVENT_CREATE);
            }

            #pragma omp single
            MPI_Barrier(MPI_COMM_WORLD);
        }
        EXTRAE_EXIT(EVENT_SIMULATION);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    fTimeEnd=MPI_Wtime();
    double wTimeHost = fTimeEnd-fTimeStart;

    if( iMyRank==0 ) {
        printf("#R%d: Computations with normal tasking took %.5f\n", iMyRank, wTimeHost);
    }

    EXTRAE_ENTER(EVENT_VALIDATION);
    LOG(iMyRank, "Validation:");
    pass = true;
    if(numberOfTasks>0) {
        for(int t=0; t<numberOfTasks; t++) {
            pass &= check_test_matrix(matrices_c[t], matrixSize, matrixSize);
        }
        if(pass)
            LOG(iMyRank, "TEST SUCCESS");
        else
            LOG(iMyRank, "TEST FAILED");
    }
    EXTRAE_EXIT(EVENT_VALIDATION);

    //deallocate matrices
    for(int i=0; i<numberOfTasks; i++) {
    	delete[] matrices_a[i];
    	delete[] matrices_b[i];
    	delete[] matrices_c[i];
    }

    delete[] matrices_a;
    delete[] matrices_b;
    delete[] matrices_c;

    MPI_Finalize();

#ifdef TRACE_EXTRAE
    Extrae_fini();
#endif
    return 0;
}
