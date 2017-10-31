/*

-Course-
DAT300: Data-driven support for cyber-physical systems

-Project-
Intrusion Detection for Industrial Control Networks

-Group 8-
Malama Kasanda - malama@student.chalmers.se
Vaios Taxiarchis - vaios@student.chalmers.se

*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_blas.h>

int main(){

   /* Start timer */
   int N=0,L=0,r=0,l=0;
   double t=0;

   printf("[ PASAD Algorithm ]\n");
   printf("Enter length N of subseries: ");
   scanf("%d", &N);
   printf("Enter lag L of subseries: ");
   scanf("%d", &L);
   printf("Enter statistical dimension r: ");
   scanf("%d", &r);
   printf("Enter threshold theta: ");
   scanf("%lf", &t);

   /* Arrays and variables */
   int sL=4801;
   double s[sL];
   double *X,*U,*product,*product2,*product3,*product4;
   double product_Xt_P_X;
   double product_Xt_X;
   int i=0,j=0;
   double dist=0.0;

   X = (double *)malloc( L*sizeof( double ));
   U = (double *)malloc( L*r*sizeof( double ));
   product = (double *)malloc( sizeof( double ));
   product2 = (double *)malloc( r*sizeof( double ));
   product3 = (double *)malloc( r*sizeof( double ));
   product4 = (double *)malloc( sizeof( double ));

   /* Basic program to read all contents from a file */
   /* File Descriptor to read from .txt */
   FILE *file;
   char ch='a';
   int flag=0;
   printf(" [Initialization started...]\n");
   /*************************************************************/
   /* Read series matrix */
   /*************************************************************/
   file=fopen("s.txt", "r");
   printf("  >Reading s[%dx1] matrix from .txt file...",sL);
   /* Read all values to an array */
  for(i = 0; i < sL; i++){
      /* Use %lf format specifier, to read doubles with scanf */
      if (!fscanf(file, "%lf", &s[i]))
          break;
      //printf("(%d,1)=%lf\n",i,s[i]);
      /* Break the inner loop and set flag=true */
      ch = getc(file);
      if(ch == EOF){
        flag=1;
        break;
      }
      /* Break the inner loop when find \n */
      else if(ch - '0'== -38){break;}
      l++;
    }

   /* Close File Descriptor */
   fclose(file);
   printf("completed!\n");

   /*************************************************************/
   /* Read U matrix  */
   /*************************************************************/
   ch='a';
   flag=0;
   file=fopen("U.txt", "r");
   printf("  >Reading U[%dx%d] matrix from .txt file...",L,r);
   /* Read all values to an array */
   for(i = 0; i < L; i++){
      /* Break the loop if the flag is true */
      if(flag==1){break;}
      for(j = 0; j < r; j++){
        /* Use %lf format specifier, to read doubles with scanf */
        if (!fscanf(file, "%lf", &U[i+j]))
          break;
        //printf("(%d,%d)=%lf\t",i,j,U[i][j]);

        /* Break the inner loop and set flag=true */
        ch = getc(file);
        if(ch == EOF){
          flag=1;
          break;
        }
        /* Break the inner loop when find \n */
        else if(ch - '0'== -38){break;}
      }
      //printf("\n");
   }
   /* Close File Descriptor */
   fclose(file);
   printf("completed!\n");
   printf(" [Matrix multiplication started...]\n");
   /* Starting point of for loop */
    int count=0;
    count++;

    for(i=N-L+2;i<N+1;i++)
    {
      X[count] = s[i];
      count++;
    }
    count=0;

    /* Generate gsl_matrix_views */
    gsl_matrix_view A = gsl_matrix_view_array(X, 1, L);
    gsl_matrix_view B = gsl_matrix_view_array(X, L, 1);
    gsl_matrix_view C = gsl_matrix_view_array(product, 1, 1);

    /* Compute Xt * X using library */
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                 1.0, &A.matrix, &B.matrix,
                 0.0, &C.matrix);

    product_Xt_X = product[0];

    int start=N+1;
    clock_t begin = clock();

    /* Testing phase */
    for(i=start;i<sL;i++)
    {
        /* Xt[1xL] * X[Lx1] */
        /* Remove the square of the fist value in X and add the square of the new s */
        product_Xt_X = product_Xt_X - X[0]*X[0] + s[i]*s[i];

        /* Generate test vector Xtest */
        /* Shifting the array X by one to the left then adding the new s */
        for(j=0;j<L-1;j++)
        {
          X[j] = X[j+1];
        }
        X[L-1] = s[i];

        count++;
        /* Compute all products */
        /* array [ (columns) x (rows) ] */

        product_Xt_P_X = 0;
        /* P_1 = Xt[1xL] * U[Lxr] * U[rxL] * X[Lx1] */

        /* Generate gsl_matrix_views */
        gsl_matrix_view D = gsl_matrix_view_array(X, 1, L);
        gsl_matrix_view E = gsl_matrix_view_array(U, L, r);
        gsl_matrix_view F = gsl_matrix_view_array(product2, 1, r);

        /* Compute Xt * U using library */
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                        1.0, &D.matrix, &E.matrix,
                        0.0, &F.matrix);

        /* Generate gsl_matrix_views */
        gsl_matrix_view G = gsl_matrix_view_array(U, r, L);
        gsl_matrix_view H = gsl_matrix_view_array(X, L, 1);
        gsl_matrix_view I = gsl_matrix_view_array(product3, r, 1);

        /* Compute X * Ut using library */
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                        1.0, &G.matrix, &H.matrix,
                        0.0, &I.matrix);

        /* Generate gsl_matrix_views */
        gsl_matrix_view J = gsl_matrix_view_array(product2, 1, r);
        gsl_matrix_view K = gsl_matrix_view_array(product3, r, 1);
        gsl_matrix_view L = gsl_matrix_view_array(product4, 1, 1);

        /* Compute using library */
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                        1.0, &J.matrix, &K.matrix,
                        0.0, &L.matrix);

        product_Xt_P_X = product4[0];


    //dist = (product_Xt_X - product_Xt_P_X)/L;

    if(dist>t){
    //  printf("  >>>d(%d)=%lf (GENERATE ALARM!!!)\n", i, dist);
    }
    else{
    // printf("  >>>d(%d)=%lf\n", i, dist);
    }

   }
   printf("  >Number of total iterations: %d\n", i);
   printf(" [Matrix multiplication completed!]\n");

   printf(" [Program exits]\n");
   /* Stop timer */
   clock_t end = clock();

   /* Print the time elapsed */
   printf("Time elapsed: %f milliseconds\n", 1000*((double)(end - begin) / CLOCKS_PER_SEC));

   return 0;
}
