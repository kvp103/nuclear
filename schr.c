/* Purpose: To compute a wavefunction moving in 1-dimensional space with a chosen potential
 * Input: Changing the potential function. Grid size can be left unchanged/hardcoded
 * Output: a text file with the density function of the moving wavepacket
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#define hbarc 197.4      // in MeV.fm
#define hbar 6.59        // in MeV.zs^{-1}
#define nmass 938        // nucleon mass in MeV.c^{-2}
#define size 450
#define maxX 2
#define minX -2
#define maxT 0.2
#define tsteps 15

int main(int argc, char*argv[]){

/*############ Starter things ################################################*/
    double x[size];
    double t[tsteps];
    double mass = 12*nmass;
    double E = 80;
    double k0 = sqrt(2*mass*E)/hbarc;
    double sigma = 0.1;
    double x0;
    double dx = (double)(maxX-minX)/(size-1);  
    
	int i;
    for(i=0; i<size; i++) { 
        x[i] = minX+i*dx;
    }
    
    x0 = x[(int) (double) (size/ (double) 4)]; // starting position for wavepacket
    
	for(i=0; i<tsteps; i++) {
        t[i] = i*(maxT/tsteps);
    } 

/*############ Potential  ####################################################*/
    double v[size];

    for(i=0; i<size; i++) {
		 // Free particle
        v[i] = 0;                           
/*      
		// Exp wall  
		v[i] = 10000*exp(-25*pow(x[i],2));      

		// Step potential
		if(x[i]<0) {
            v[i] = 0;
        } else {
            v[i] = 6000;                        
        } 
*/
    }

/*############ Hamiltonian eigen decomp ######################################*/
    gsl_matrix *mat = gsl_matrix_alloc (size,size);
    gsl_vector *eval = gsl_vector_alloc (size);
    gsl_matrix *evec = gsl_matrix_alloc (size,size);
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc (size);
    double evals[size];
    double evecs[size][size];

	int j;
    for(i=0;i<size;i++) {
        for(j=0;j<size;j++) {
            gsl_matrix_set(mat, i, j, 0);
            if(i==j) {
                gsl_matrix_set(mat, i, j, pow(hbarc,2)/(pow(dx,2)*mass)+v[i]);
            } else if (j==i-1 && i!=0) {
                gsl_matrix_set(mat, i, j, -pow(hbarc,2)/(pow(dx,2)*2*mass));
            } else if (j==i+1 && i!= size-1) {
                gsl_matrix_set(mat, i, j, -pow(hbarc,2)/(pow(dx,2)*2*mass));
            }              
        }
    }

    gsl_eigen_symmv (mat, eval, evec, w);

    for(i=0;i<size;i++) {
        evals[i] = gsl_vector_get(eval, i);
        for(j=0;j<size;j++) {
            evecs[i][j] = gsl_matrix_get(evec, i, j);
        }
    }

    gsl_eigen_symmv_free (w);
    gsl_vector_free (eval);
    gsl_matrix_free (evec);
    gsl_matrix_free (mat);

/*############ Summing #######################################################*/
    long double complex psi0[size];
    long double complex ip[size];
    long double complex PSI[size][tsteps];
    long double complex sum=0;

    // Set initial wavepacket Psi0
    for(i=0; i<size; i++) {
        psi0[i] = cexpl(1*I*k0*x[i]-(pow((x[i]-x0),2)/pow(sigma,2)));
    }

    // Computes inner product < \psi_\mu | \psi_0 >
    for(i=0; i<size; i++) {
        for(j=0; j<size; j++) {
            sum = conj(evecs[j][i])*psi0[j] + sum;
        }
        ip[i] = sum;
        sum=0;
    }
    sum=0;

    // makes big wavefunction 
	int k;
    for(i=0; i<tsteps; i++) {
        for(j=0; j<size; j++) {
            for(k=0; k<size; k++) {     
                sum = cexpl(-1*I*evals[k]*t[i]/hbar)*ip[k]*evecs[j][k] + sum;
            }
            PSI[j][i] = sum;
            sum=0;
        }
    }

/*############ Density #######################################################*/
    double complex dens[size][tsteps];
    
    for(i=0; i<size; i++) {
        for(j=0; j<tsteps; j++) {
            dens[i][j] = creal(conj(PSI[i][j])*PSI[i][j]);
        }
    } 

/*############ Output wavefunction density values ############################*/
    FILE *output;
    output = fopen("dens.dat", "w");

    for(i=0; i<size; i++) {
        fprintf(output, "%10lf    ", x[i]);
        for(j=0; j<tsteps; j++) {
            fprintf(output, "%e   ", creal(dens[i][j]));
        }
        fprintf(output, "\n"); 
    }

    fclose(output);
    return 0;
}
