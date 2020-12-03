/*
Purpose: To compute an internuclear Akyuz--Winther potential (and barrier) 
Input: Two nuclei specified by Z and A numbers (as input when running the code)
Output: Text file "AW.dat" with several potential values at each separation distance.
		Potential values include Coulomb (C), Akyuz--Winther (AW), Woods--Saxon (WS) 
		and Total potential (AW+C)
*/


// Libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Physics constants
#define a 0.63              // fm
#define k 8.988e9           // N m^2 C^{-2}
#define ep 1.602e-19        // C
#define ev 6.242e18         // convert J --> eV

// Grid constants
#define MAX 21              // fm
#define MIN 2               // fm
#define dx 0.02             // fm

// Global variables
float size; // to perform division to calculate size of array
int size2; // convert size to integer to create array

// Function list
void potentials(int A1, int Z1, int A2, int Z2, double *Vt);
void delete(double Vt[]);
void barrier(double Vt[]);

//#########################################//
//               MAIN                      //
//#########################################//

int main(int argc, char*argv[]){
    int A1 = atoi(argv[1]);
    int Z1 = atoi(argv[2]);
    int A2 = atoi(argv[3]);
    int Z2 = atoi(argv[4]);

    size = (MAX-MIN)/dx+1;
    size2 = (int) size;
    double Vt[size2]; 		//make total potential initial array

	// Compute potential function and save to file
    potentials(A1, Z1, A2, Z2, Vt);
	
	// Remove inner potential
    delete(Vt);
	
	// Print barrier to screen
    barrier(Vt);
    return 0;
}

//#########################################//
//             FUNCTIONS                   //
//#########################################//

//###########################################################################//
void potentials(int A1, int Z1, int A2, int Z2, double *Vt) {

    double r0, r1, r2, rr;
    double S0, ainv;
    double r[size2], Vaw[size2], Vc[size2], Vws[size2];

	// Create radial separation vector
    int i;
    for(i=0; i<size2; i++) {
        r[i] = MIN+i*dx;
    }

	// Compute parameters
    r1 = (1.2*pow(A1,pow(3,-1)))-0.35; 	//radius nucleus 1
    r2 = (1.2*pow(A2,pow(3,-1)))-0.35;  //radius nucleus 2
    rr = (r1*r2)/(r1+r2); 				//reduced radius
    r0 = r1+r2;
    ainv = 1.16*(1+0.48*(pow(A1,pow(-3,-1))+pow(A2,pow(-3,-1))));
    S0 = 65.4;

	// Calculate AW pot and coulomb pot. Total pot is their sum
    for(i=0; i<size2; i++) {
        Vc[i] = k*ep*ep*Z1*Z2/(r[i]*1e-15)*ev*1e-6;			//Coulomb
        Vaw[i] = -S0*rr*exp(ainv*(r0-r[i]));				//Akyuz--Winther
        Vws[i] = -100/(1+exp((r[i]-1.2*(pow(A1,pow(3,-1))+pow(A2,pow(3,-1))))/0.6));
															//Woods-Saxon
        Vt[i] = Vaw[i]+Vc[i];								//Total
    }


    // Write potential function to text file
	FILE *output;
	output = fopen("AW.dat", "w");

	for(i=0; i<size2; i++) {
		fprintf(output, "%6.2lf    %10.4lf    %10.4lf    %10.4lf    %10.4lf \n",
				r[i], Vt[i], Vc[i], Vaw[i], Vws[i]);
	}

	fclose(output);

}

//###########################################################################//
void delete(double Vt[]) {
// delete crap entries
    int i;
    for(i=0; i<size2-1; i++) {
        if(Vt[i] > Vt[i+1]) {
            Vt[i] =0;
        } else if(Vt[i] < Vt[i+1]) {
            break;
        }
    }
}

//###########################################################################//
void barrier(double Vt[]) {
// find potential barrier
    int i;
    double max=0;
    for(i=0; i<size2-1; i++) {
        if(Vt[i+1] > Vt[i]) {
            max = Vt[i+1];
        }
    }
    printf("--------------------------- \n");
    printf("    Vb is %.3lf MeV\n", max);
    printf("--------------------------- \n");
}
