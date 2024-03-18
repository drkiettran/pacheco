#include <iostream>
#include <string>
#include <mpi.h>

void Get_input(int my_rank, int comm_sz, double* a_p, double* b_p, int* n_p) {
    int dest;
    if (my_rank == 0) {
        std::cout << "Enter, b, and n" << std::endl;
        scanf("%lf %lf %d", a_p, b_p, n_p);
        for (dest = 1; dest < comm_sz; dest++) {
            MPI_Send(a_p, 1, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
            MPI_Send(b_p, 1, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
            MPI_Send(n_p, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
        }
    } else {
        MPI_Recv(a_p, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(b_p, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(n_p, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

double f(double x) {
    return x*x + 1.0;
}


double Trap(double left_endpt, double right_endpt, int trap_count, double base_len) {
    double estimate, x;
    int i;

    estimate = (f(left_endpt) + f(right_endpt))/2.0;
    for (i = 1; i < trap_count; i++) {
        x = left_endpt + i * base_len;
        estimate += f(x);
    }
    return estimate * base_len;
}

/**
* Version 1 from Quinn textbook.
*
* - Each process is assigned a local_n that is a subset of n evenly, 
*   including process 0.
* - Ideally, the number of processes should be divisble by n.
*/
int version_1(int argc, char* argv[])
{
    MPI_Comm comm;
    int num_processes, proc_num;
    int n = 40, local_n;
    double a = 0.0, b = 3.0, h, local_a, local_b;
    double local_int, total_int;
    int source;

    std::cout << "Pacheco's version 1 ...\n\n";

    MPI_Init(&argc, &argv);

    comm = MPI_COMM_WORLD;

    MPI_Comm_size(comm, &num_processes);
    MPI_Comm_rank(comm, &proc_num);
    
    h = (b-a)/n;
    local_n = n / num_processes;
    
    local_a = a + proc_num * local_n * h;
    local_b = local_a + local_n * h;
    local_int = Trap(local_a, local_b, local_n, h);
    std::cout << ">>> proc_num: " << proc_num 
            << "; local_a: " << local_a
            << "; local_b: " << local_b
            << "; local_n: " << local_n
            << "; local_int: " << local_int << std::endl;


    if (proc_num != 0) {
        MPI_Send(&local_int, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        std::cout << ">>\tproc_num: " << proc_num << " sent ..." << std::endl;
    } else {
        total_int = local_int;
        for (source = 1; source < num_processes; source++) {
            MPI_Recv(&local_int, 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            total_int += local_int;
            std::cout << "\t<< proc_num " << source << " received ..."
                      << "total_int: " << total_int << std::endl;
        }
    }

    if (proc_num == 0) {
        printf("with n = %d trapezoids, our estimate\n", n);
        printf("of the integral from %f to %f - %.15e\n", a, b, total_int);
    }
    MPI_Finalize();
    return 0;
}

/**
 * n: number of segments of [a,b]
 * size: number of cores/processes in the world.
 * 
 * each core process n/size segments. if n is not divisible by size, 
 * left over segments are processed by rank 0 (i.e., size 6; n = 40; remainder (40,n)=4)
*/
int version_1_revised(int argc, char* argv[])
{
    MPI_Comm comm;
    int size, rank;
    int n = 40, local_n;
    double a = 0.0, b = 3.0, h, local_a, local_b;
    double local_int, total_int;
    int source;

    std::cout << "\n\nPacheco's version 1 revised ...\n\n";
    MPI_Init(&argc, &argv);

    comm = MPI_COMM_WORLD;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);
    
    h = (b-a)/n;
    local_n = n / size;
    
    local_a = a + rank * local_n * h;
    local_b = local_a + local_n * h;
    local_int = Trap(local_a, local_b, local_n, h);
    std::cout << ">>> size: " << size
            << "; local_a: " << local_a
            << "; local_b: " << local_b
            << "; local_n: " << local_n
            << "; local_int: " << local_int << std::endl;


    if (rank != 0) {
        MPI_Send(&local_int, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        std::cout << "\t>>proc_num: " << size << " sent ..." << std::endl;
    } else {
        total_int = local_int;
        // if remainder exists, process the rest by process 0.
        int remainder_tasks = n % size;
        if (remainder_tasks > 0) {
            std::cout << "*** remainder tasks ... " << remainder_tasks << std::endl;
            local_a = a + size * local_n * h;
            local_b = b;
            total_int += Trap(local_a, local_b, remainder_tasks, h);
        }

        for (source = 1; source < size; source++) {
            MPI_Recv(&local_int, 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            total_int += local_int;
            std::cout << "\t<< rank " << source << " received ..."
                      << "total_int: " << total_int << std::endl;
        }
    }

    if (rank == 0) {
        printf("with n = %d trapezoids, our estimate\n", n);
        printf("of the integral from %f to %f - %.15e\n", a, b, total_int);
    }
    MPI_Finalize();
    return 0;
}

int main (int argc, char* argv[]) {

    // version_1(argc, argv);
    return version_1_revised(argc, argv);
}