#ifndef func_mpi_h
#define func_mpi_h
////----------------------------------------------
#include <mpi.h>
#include "func_c.h"
bool msgTest(MPI_Request *req){
    int flag; MPI_Status stat;
    MPI_Test(req, &flag, &stat);
    if(flag && stat.MPI_SOURCE!=MPI_ANY_SOURCE) return true;
    return false;
}
//end of func_mpi.h
#endif
