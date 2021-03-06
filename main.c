#include <iostream>
#include <cmath>

#define pi 3.14159265359

// gauss points and weights
static double gp[] = {0.0, -0.2695431559523450, 0.2695431559523450, -0.5190961292068118, 0.5190961292068118,
	                  -0.7301520055740494, 0.7301520055740494, -0.8870625997680953, 0.8870625997680953,
	                  -0.9782286581460570, 0.9782286581460570};

static double gw[] = {0.2729250867779006, 0.2628045445102467, 0.2628045445102467, 0.2331937645919905, 
                      0.2331937645919905, 0.1862902109277343, 0.1862902109277343, 0.1255803694649046,
                      0.1255803694649046, 0.0556685671161737, 0.0556685671161737};

static const int np = 10;


inline double fr(double theta, double a, double z){

	double temp = cos(theta);
	double r = sqrt( 1 + a * a - 2 * a * temp + z * z);
	double result;
	return result = r * temp + temp * temp * log(r + a - temp);		

}


inline double fz(double theta, double a, double z){
	
	double temp = cos(theta);
	double r = sqrt( 1 + a * a - 2 * a * temp + z * z);
	double temp1 = sin(theta);
	double result = z * log(r + a - temp);
	result += (z/fabs(z) / 2.0) * temp * log((r - fabs(z) )/ (r + fabs(z)));
	result -= (z/fabs(z)) * temp1 * atan((fabs(z) * (a - temp)) / (r * temp1));	

	return result;	
}

inline double F_rho(double a, double z){

	double pi2 = pi/2;
	double theta;

	double integral = 0.0;	
	for(int j = 0; j < np; j++){
		theta = gp[j] * pi2 + pi2;
		integral += fz(theta, a, z) * gw[j];
	}

}

inline double F_z(double a, double z){
	
	double pi2 = pi/2;
	double theta;

	double integral = 0.0;	
	for(int j = 0; j < np; j++){
		theta = gp[j] * pi2 + pi2;
		integral += fr(theta, a, z) * gw[j];
	}	
}

double *B_rho(double a1, double a2, double z1, double z2, double *rho, double *z, size_t num)
{
	double A1, A2, Z1, Z2, K; // how to set K

	double *result = new double[num]; 

	for(int i = 0; i < num ; i++){

		A1 = a1 / rho[i]; A2 = a2 / rho[i]; Z1 = (z1 - z[i])/ rho[i]; Z2 = (z2  - z[i] )/ rho[i];
		result[i] = F_rho(A2, Z2) - F_rho(A2, Z1) - F_rho(A1, Z2) - F_rho(A1, Z1);

	}
	return result;
}


double *B_z(double a1, double a2, double z1, double z2, double *rho, double *z, size_t num)
{
	double A1, A2, Z1, Z2, K; // how to set K

	double *result = new double[num];

	for(int i = 0; i < num ; i++){
		A1 = a1 / rho[i]; A2 = a2 / rho[i]; Z1 = (z1 - z[i])/ rho[i], Z2 = (z2  - z[i] )/ rho[i];	
		result[i] = F_z(A2, Z2) - F_z(A2, Z1) - F_z(A1, Z2) - F_z(A1, Z1);
	}

	return result;
}

int main()
{

	int ncoil = 8; 
	size_t num = 200;

	double *rho = new double[num];
	double *z = new double[num]; 

	double *a1s = new double[ncoil];
	double *a2s = new double[ncoil];
	double *z1s = new double[ncoil];
	double *z2s = new double[ncoil];

	double **allresults_r = new double*[ncoil];
	double **allresults_z = new double*[ncoil];	

	// Some Initialization here, arbitrary assignments now
	for(int i = 0; i < num; i++){
		rho[i] = i * 0.1;
		z[i] = i * 0.3;
	}
    // no idea how to set a1,a2, z1, z2, all the same here 
	for (int i = 0; i < ncoil; ++i)
	{
		a1s[i] = -1.0;
		a2s[i] = 1.0;
		z1s[i] = -1.0;
		z2s[i] = 1.0; 
	}

	for (int i = 0; i < ncoil; ++i)
	{
		allresults_r[i] = B_rho(a1s[i], a2s[i], z1s[i], z2s[i], rho, z, num);
		allresults_z[i] = B_z(a1s[i], a2s[i], z1s[i], z2s[i], rho, z, num);
	}


	// Free memory

	delete []rho;
	delete []z;
	delete []a1s;
	delete []a2s;
	delete []z1s;
	delete []z2s;


	for (int i = 0; i < ncoil; ++i){
		delete [] allresults_r[i];
		delete [] allresults_z[i];
	}

	delete [] allresults_z;
	delete [] allresults_r;
}
