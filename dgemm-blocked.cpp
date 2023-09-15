// const char *dgemm_desc = "Blocked dgemm.";
// /* This routine performs a dgemm operation
//  *  C := C + A * B
//  * where A, B, and C are n-by-n matrices stored in column-major format.
//  * On exit, A and B maintain their input values. */
// void square_dgemm_blocked(int n, int block_size, double *A, double *B, double *C)
// {
//    // insert your code here
//    for (int i0 = 0; i0 < n; i0 += block_size)
//    {

//       for (int j0 = 0; j0 < n; j0 += block_size)
//       {

//          for (int k0 = 0; k0 < n; k0 += block_size)
//          {

//             for (int j1 = j0; j1 < j0 + block_size; ++j1)
//             {

//                for (int i1 = i0; i1 < i0 + block_size; ++i1)
//                {

//                   for (int k1 = k0; k1 < k0 + block_size; ++k1)
//                   {
//                      C[n * i1 + j1] += A[n * k1 + j1] * B[n * i1 + k1];
//                   }
//                }
//             }
//          }
//       }
//    }
// }
const char *dgemm_desc = "Blocked dgemm.";
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>
#include <string.h>
/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are n-by-n matrices stored in column-major format.
 * On exit, A and B maintain their input values. */
void square_dgemm_blocked(int n, int block_size, double *A, double *B, double *C)
{
   std::vector<double> buf(3 * block_size * block_size);
   double *Ccopy = buf.data() + 0;
   double *Bcopy = Ccopy + block_size * block_size;
   double *Acopy = Bcopy + block_size * block_size;

   for (int i0 = 0; i0 < n; i0 += block_size)
   {
      for (int j0 = 0; j0 < n; j0 += block_size)
      {
         for (int copyIndex = 0; copyIndex < block_size; copyIndex++)
         {
            memcpy(&Ccopy[copyIndex * block_size], &C[n * copyIndex + i0 * n + j0], sizeof(double) * block_size);
         }
         for (int k0 = 0; k0 < n; k0 += block_size)
         {
            for (int copyIndex = 0; copyIndex < block_size; copyIndex++)
            {
               memcpy(&Acopy[copyIndex * block_size], &A[n * copyIndex + k0 * n + j0], sizeof(double) * block_size);
            }
            for (int copyIndex = 0; copyIndex < block_size; copyIndex++)
            {
               memcpy(&Bcopy[copyIndex * block_size], &B[n * copyIndex + k0 + i0 * n], sizeof(double) * block_size);
            }

            for (int i1 = 0; i1 < block_size; i1++)
            {
               for (int j1 = 0; j1 < block_size; j1++)
               {
                  for (int k1 = 0; k1 < block_size; k1++)
                  {
                     Ccopy[i1 * block_size + j1] += Acopy[k1 * block_size + j1] * Bcopy[i1 * block_size + k1];
                  }
               }
            }
            for (int copyIndex = 0; copyIndex < block_size; copyIndex++)
            {
               memcpy(&C[n * copyIndex + i0 * n + j0], &Ccopy[copyIndex * block_size], sizeof(double) * block_size);
            }
         }
      }
   }
}