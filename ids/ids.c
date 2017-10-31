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
   double U[L][r];
   double X[L];
   double product_Xt_P_X;
   double product_Xt_X;
   int i=0,j=0,c=0,d=0,k=0;
   double sum=0.0,dist=0.0;

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
        if (!fscanf(file, "%lf", &U[i][j]))
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
   printf(" [Initialization completed!]\n");

    printf(" [Matrix multiplication started...]\n");
   /* Starting point of for loop */
    int count=0;

    X[count] = 0;
    count++;

    for(i=N-L+2;i<N+1;i++)
    {
        X[count] = s[i];
        count++;
    }
    count=0;

    for (c = 0; c < L; c++)
    {
        sum = sum + X[c]*X[c];
    }
    product_Xt_X = sum;
    sum = 0;

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
        /* P_1 = Xt[1xL] * U[Lxr] */
        for (d = 0; d < r; d++)
        {
            for (k = 0; k < L; k++)
            {
                sum = sum + X[k]*U[k][d];
            }
            product_Xt_P_X = product_Xt_P_X + sum*sum;
            sum = 0;
        }

    dist = (product_Xt_X - product_Xt_P_X)/L;

    if(dist>t){
    //  printf("  >>>d(%d)=%lf (GENERATE ALARM!!!)\n", i, dist);
    }
    else{
    //printf("  >>>d(%d)=%lf\n", i, dist);
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
