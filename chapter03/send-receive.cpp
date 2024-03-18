#include <iostream>
#include <string>
#include <mpi.h>

/*
    - n processes are sprung up.
    - process rank 0 sends out a greeting message to 
      the rest of the processes.
    - each of the processes rank 1 to n-1, receive a greeting
      message from process 0; then, print the message.
*/
int send_receive(int argc, char* argv[])
{
    MPI_Comm comm;
    int num_processes, proc_num;

    MPI_Init(&argc, &argv);

    comm = MPI_COMM_WORLD;

    MPI_Comm_size(comm, &num_processes);
    MPI_Comm_rank(comm, &proc_num);
    
    
    if (proc_num == 0) {
        std::cout << "==> There are " << num_processes << " processes and I am process number " << proc_num << std::endl;
        std::cout << "<< I am sending out messages to others ..." << std::endl;
        for (int dest_proc = 1; dest_proc < num_processes; dest_proc++) {
            std::string msg = "Hello, " + std::to_string(dest_proc) + "!";
            
            MPI_Send(msg.c_str(), msg.length(), MPI_CHAR, dest_proc, dest_proc, comm);
        }
    } else {
        std::cout << ">> I am process: " << proc_num << " and receiving a message from process 0" << std::endl;
        char buf[100];
        MPI_Status status;

        memset(buf, 0, sizeof(buf));
        MPI_Recv(&buf, sizeof(buf), MPI_CHAR, 0, proc_num, comm, &status);
        std::cout << ">> I received " << buf << " from process 0 tag " << status.MPI_TAG << std::endl;
    }
    MPI_Finalize();
    return 0;
}

int main(int argc, char* argv[]) {
    return send_receive(argc, argv);
}
