#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<complex.h>

#include "funciones_Matlab.h"
#include "funciones.h"

#define MAT(vec,i,j,num_col) (vec[(j) + num_col*(i)])    //macro para manejo de matrices arroja un solo valor
//=============================== Funciones MATLAB ==============================================
//sumIvn: Suma todos los valores de un array
double sumIvn(double *matriz, int tam){
	double resultado=0;
	int i;
	for (i=0;i<tam;i++){
		resultado=resultado+matriz[i];	
	}
	return resultado;	
}

// zerosIvn: Genera un vector de zeros del tamaño de filas,columnas
void zerosIvn(double *vector, int filas, int columnas){
	int i,j;
	for (i=0;i<filas;i++){
		for (j=0;j<columnas;j++){ 
			MAT(vector,i,j,columnas) = 0; 
		}
	}
}

void linspaceIvn(double a,double b, int N,double *salida){
	double delta;
	int i;
	
	delta = (b-a)/(N-1);
	salida[0]=a;
	salida[N-1]=b;
	for (i=1;i<(N-1);i++){
		salida[i]= salida[i-1] + delta;
	}
}

double minIvn(double *vector, int tam){
	int i;
	double minimo;
	minimo=vector[0];
	
	for (i=1;i<tam;i=i+1){
		if (vector[i]<minimo){
			minimo=vector[i];			
		}
	}
	return minimo;
}

double maxIvn(double *vector, int tam){
	int i;
	double maximo;
	maximo=vector[0];
	
	for (i=1;i<tam;i=i+1){
		if (vector[i]>maximo){
			maximo=vector[i];			
		}
	}
	return maximo;
}
/*  URL: https://gist.github.com/svdamani/1015c5c4b673c3297309
	Shreevardhan Damani */
void splineCoefIvn(double *x, double *a, int tam, double *coef){
	/** Step 0 */
    int i, j, n;
	n = tam-1;
	
	double h[n], A[n], l[n + 1], u[n + 1], z[n + 1], c[n + 1], b[n], d[n];
    /** Step 1 */
    for (i = 0; i <n; ++i){
		h[i] = x[i + 1] - x[i];
	}	
	
			/** Step 2 */
    for (i = 1; i <= n - 1; ++i){
        A[i] = 3 * (a[i + 1] - a[i]) / h[i] - 3 * (a[i] - a[i - 1]) / h[i - 1];
	}
	
    /** Step 3 */
    l[0] = 1;
    u[0] = 0;
    z[0] = 0;

    /** Step 4 */
    for (i = 1; i <= n - 1; ++i) {
        l[i] = 2 * (x[i + 1] - x[i-1]) - h[i - 1] * u[i - 1];
        u[i] = h[i] / l[i];
        z[i] = (A[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    /** Step 5 */
    l[n] = 1;
    z[n] = 0;
    c[n] = 0;

    /** Step 6 */
    for (j = n - 1; j >= 0; --j) {
        c[j] = z[j] - u[j] * c[j + 1];
        b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3;
        d[j] = (c[j + 1] - c[j]) / (3 * h[j]);
    }
	
	/** Step 7 */
    for (i = 0; i < n; ++i){
        //printf("%2d %8.2f %8.2f %8.2f %8.2f\n", i, a[i], b[i], c[i], d[i]);
		MAT(coef,i,0,4) = d[i]; MAT(coef,i,1,4) = c[i]; MAT(coef,i,2,4) = b[i]; MAT(coef,i,3,4) = a[i];
	}
}
	
/* Author: Ivan Obregon*/
/* Returns a vector with the evaluate points using a natural cubic spline coefficients matrix */
void evalCoefSpline(double *vec, int n_vec, double *coef, double *x, int n_x, /*Out: */ double *salida){
	int i,j;
	for (i=0;i<n_vec ;i++){
		for (j=0;j<n_x-1;j++){
			
			if 	(vec[i]<x[0]){
				salida[i] = MAT(coef,0,0,4)*pow(vec[i] - x[0],3) + MAT(coef,0,1,4)*pow(vec[i] - x[0],2) + MAT(coef,0,2,4)*pow(vec[i] - x[0],1) + MAT(coef,0,3,4);
				//printf("\n menor");
				break;
			}else if (vec[i]>x[n_x-1]){
				salida[i] = MAT(coef,n_x-2,0,4)*pow(vec[i] - x[n_x-1],3) + MAT(coef,n_x-2,1,4)*pow(vec[i] - x[n_x-1],2) + MAT(coef,n_x-2,2,4)*pow(vec[i] - x[n_x-1],1) + MAT(coef,n_x-2,3,4);
				//printf("\n mayor");
				break;
			}else if ( (vec[i]>=x[j]) && (vec[i]<=x[j+1]) ){
				salida[i] = MAT(coef,j,0,4)*pow(vec[i] - x[j],3) + MAT(coef,j,1,4)*pow(vec[i] - x[j],2) + MAT(coef,j,2,4)*pow(vec[i] - x[j],1) + MAT(coef,j,3,4);
				//printf("\n mitad");
				break;
			}
		}
	}
}

//Normalizar=1 Aplica la normalizacion de los MU, MU0=0 y MU1=1 Para que no afecte
int polyfitIvn(double* XX_temp, double* YY,int countOfElements,  int order, double* coefficients, int Normalizar, double MU0, double MU1)
{
	double *XX = (double *) malloc(countOfElements*sizeof(double));
	if (Normalizar == 1){
		int ij;
		for (ij=0;ij<countOfElements;ij++){
			XX[ij] = (XX_temp[ij] - MU0)/MU1;
		}
	}else{
		memcpy(XX,XX_temp,countOfElements*sizeof(double));
	}
	
	double *dependentValues=XX;
	double *independentValues=YY;
    // Declarations...
    // ----------------------------------
    enum {maxOrder = 20};
    double B[maxOrder+1] = {0.0f};
    double P[((maxOrder+1) * 2)+1] = {0.0f};
    double A[(maxOrder + 1)*2*(maxOrder + 1)] = {0.0f};
    double x, y, powx;
    unsigned int ii, jj, kk;
    // Verify initial conditions....
    // ----------------------------------
    // This method requires that the countOfElements > 
    // (order+1) 
    if (countOfElements <= order)
        return -1;
    // This method has imposed an arbitrary bound of
    // order <= maxOrder.  Increase maxOrder if necessary.
    if (order > maxOrder)
        return -1;
    // Begin Code...
    // ----------------------------------
    // Identify the column vector
    for (ii = 0; ii < countOfElements; ii++)
    {
        x    = dependentValues[ii];
        y    = independentValues[ii];
        powx = 1;

        for (jj = 0; jj < (order + 1); jj++)
        {
            B[jj] = B[jj] + (y * powx);
            powx  = powx * x;
        }
    }
    // Initialize the PowX array
    P[0] = countOfElements;
    // Compute the sum of the Powers of X
    for (ii = 0; ii < countOfElements; ii++)
    {
        x    = dependentValues[ii];
        powx = dependentValues[ii];

        for (jj = 1; jj < ((2 * (order + 1)) + 1); jj++)
        {
            P[jj] = P[jj] + powx;
            powx  = powx * x;
        }
    }
    // Initialize the reduction matrix
    //
    for (ii = 0; ii < (order + 1); ii++)
    {
        for (jj = 0; jj < (order + 1); jj++)
        {
            A[(ii * (2 * (order + 1))) + jj] = P[ii+jj];
        }

        A[(ii*(2 * (order + 1))) + (ii + (order + 1))] = 1;
    }
    // Move the Identity matrix portion of the redux matrix
    // to the left side (find the inverse of the left side
    // of the redux matrix
    for (ii = 0; ii < (order + 1); ii++)
    {
        x = A[(ii * (2 * (order + 1))) + ii];
        if (x != 0)
        {
            for (kk = 0; kk < (2 * (order + 1)); kk++)
            {
                A[(ii * (2 * (order + 1))) + kk] = 
                    A[(ii * (2 * (order + 1))) + kk] / x;
            }

            for (jj = 0; jj < (order + 1); jj++)
            {
                if ((jj - ii) != 0)
                {
                    y = A[(jj * (2 * (order + 1))) + ii];
                    for (kk = 0; kk < (2 * (order + 1)); kk++)
                    {
                        A[(jj * (2 * (order + 1))) + kk] = 
                            A[(jj * (2 * (order + 1))) + kk] -
                            y * A[(ii * (2 * (order + 1))) + kk];
                    }
                }
            }
        }
        else
        {
            // Cannot work with singular matrices
            return -1;
        }
    }
    // Calculate and Identify the coefficients
    for (ii = 0; ii < (order + 1); ii++)
    {
        for (jj = 0; jj < (order + 1); jj++)
        {
            x = 0;
            for (kk = 0; kk < (order + 1); kk++)
            {
                x = x + (A[(ii * (2 * (order + 1))) + (kk + (order + 1))] *
                    B[kk]);
            }
            coefficients[order - ii] = x;
        }
    }
	free(XX);
    return 0;
}

double polyvalIvn(double *vecPoly, int tamPoly, double X, double MU0, double MU1){  //MU0=0 y MU1=1 Para que no afecte
	double nuevoX;
	double acum = 0;
	int i;
	nuevoX = (X - MU0)/(MU1);
	
	for (i=0;i<tamPoly;i=i+1){
		acum = acum + vecPoly[i]*pow(nuevoX,(tamPoly-i-1));
	}
	
	return acum;	
}

double meanIvn(double *vector, int tamVector){
	int i;
	double mean;
	mean=0;
	for (i=0;i<tamVector;i=i+1){
		mean = mean + vector[i];
	}
	mean = mean /(tamVector);
	
	return mean;
}

double stdIvn(double *vector, int tamVector){
	int i;
	double mean;
	double std;
	mean=0;
	for (i=0;i<tamVector;i=i+1){
		mean = mean + vector[i];
	}
	mean = mean /(tamVector);
	
	for (i=0;i<tamVector;i=i+1){
		std = std + pow(vector[i] - mean,2);
	}
	std = sqrt((1.0/(tamVector-1))*std);
	
	return std;	
}

//=========================================================================================
// Solo sirven para binarios guardados en double
void getMatrixSizeIvn(int *filas, int *columnas, char *nombreBin){
  double *tam = (double *) malloc(2*sizeof(double));
  FILE *archivoBin;
  archivoBin = fopen(nombreBin,"rb");
  fread(tam,sizeof(double),2,archivoBin);
  fclose(archivoBin);
  *filas=(int)tam[0];
  *columnas=(int)tam[1];
  free(tam);
}

void getMatrixIvn(double *vector, char *nombreBin,int tam){
  FILE *archivoBin;
  archivoBin = fopen(nombreBin,"rb");
  fseek ( archivoBin , 2*sizeof(double) , SEEK_SET );
  fread(vector,sizeof(double),tam,archivoBin);
  fclose(archivoBin);
} 

void roots3Ivn(double *polinomio, double *raices_real, double *raices_imag){
  double a,b,c,d;
  double complex raiz1;
  double complex raiz2;
  double complex raiz3;
  a=polinomio[0];
  b=polinomio[1];
  c=polinomio[2];
  d=polinomio[3];
  //printf("a=%f b=%f  c=%f  d=%f\n",a,b,c,d);
  raiz1 = ((-b)/(3*a)) - (((cpow(2,1.0/3.0))*(-(b*b) + 3*a*c))/(3*a*cpow(-2*(b*b*b) + 9*a*b*c - 27*a*a*d + csqrt(4*cpow(-(b*b) + 3*a*c,3) + cpow(-2*(b*b*b) + 9*a*b*c - 27*a*a*d,2)),1.0/3.0)))+  (cpow(-2*(b*b*b) + 9*a*b*c - 27*a*a*d + csqrt(4*cpow(-(b*b) + 3*a*c,3) + cpow(-2*(b*b*b) + 9*a*b*c - 27*a*a*d,2)),1.0/3.0)/((3*pow(2,(1.0/3.0))*a)));
  raiz2 = ((-b)/(3*a)) + ((  ((1 + csqrt(3)*I))*(-(b*b) + 3*a*c))/(3*a*(cpow(2,(2.0/3.0)))*cpow(-2*(b*b*b) + 9*a*b*c - 27*a*a*d + csqrt(4*cpow(-(b*b) + 3*a*c,3) + cpow(-2*(b*b*b) + 9*a*b*c - 27*a*a*d,2)),1.0/3.0)))  -  ((1-csqrt(3)*I)*cpow(-2*(b*b*b) + 9*a*b*c - 27*a*a*d + csqrt(4*cpow(-(b*b) + 3*a*c,3) + cpow(-2*(b*b*b) + 9*a*b*c - 27*a*a*d,2)),1.0/3.0)/((6*pow(2,(1.0/3.0))*a)));
  raiz3 = ((-b)/(3*a)) + ((  ((1 - csqrt(3)*I))*(-(b*b) + 3*a*c))/(3*a*(cpow(2,(2.0/3.0)))*cpow(-2*(b*b*b) + 9*a*b*c - 27*a*a*d + csqrt(4*cpow(-(b*b) + 3*a*c,3) + cpow(-2*(b*b*b) + 9*a*b*c - 27*a*a*d,2)),1.0/3.0)))  -  ((1+csqrt(3)*I)*cpow(-2*(b*b*b) + 9*a*b*c - 27*a*a*d + csqrt(4*cpow(-(b*b) + 3*a*c,3) + cpow(-2*(b*b*b) + 9*a*b*c - 27*a*a*d,2)),1.0/3.0)/((6*pow(2,(1.0/3.0))*a)));
  raices_real[0]=creal(raiz1);
  raices_real[1]=creal(raiz2);
  raices_real[2]=creal(raiz3);
  raices_imag[0]=cimag(raiz1);
  raices_imag[1]=cimag(raiz2);
  raices_imag[2]=cimag(raiz3);
}
//=========rootsIvnAll===================================================================================
void Print_results(double x, double y, int *cont, double *salida_real, double *salida_imag){
//Label: l10
  double a1; int i;
  //printf("%13.8f  %13.8f \n", x,y);
  
  salida_real[cont[0]] = x;
  salida_imag[cont[0]] = y;
  cont[0] = cont[0] + 1;
//calculate error estimation (optional)
}

//Solve x^2+p1x+q=0
void Solve2(double *x, double *y, double p1, double q, int *cont, double *salida_real, double *salida_imag){
//Label: l10
  double d;
  d=p1*p1-4*q;
  if (d<0) goto l10;
  d=sqrt(d); y[0]=0.0;
  x[0]=(-p1+d)/2.0; Print_results(x[0],y[0],cont,salida_real,salida_imag);
  x[0]=(-p1-d)/2.0; Print_results(x[0],y[0],cont,salida_real,salida_imag);
  return;
l10: d=sqrt(-d)/2.0; x[0]=-p1/2.0;
  y[0]=d; Print_results(x[0],y[0],cont,salida_real,salida_imag);
  y[0]=-d; Print_results(x[0],y[0],cont,salida_real,salida_imag);
}

void rootsIvnAll(double *poli, int n, double *salida_real, double *salida_imag){
  double A[20],B[20],C[20],P[20];
  double a1,b1,d,e,f,p1,q,r,t,u,v,x,y,z;	
  int    i,j,k,l;
  int 	 cont = 0;
	
//Enter polynomial A(i) and copy in P(i)
  for (i=0; i<=n; i++) {
    A[i] = poli[i];
	P[i]=A[i];
  }
  //printf("\n             ROOTS                ERROR\n");
//Init section
  p1=0.0; q=0.0; k=100; e=0.00001;
//Factoring main loop
l100:if (n<=2) goto l500;
  j=0;
l200: if (j>k) goto l1000;
  j++;
//calculate B(i) and C(i)
  B[0]=0.0; B[1]=0.0; C[0]=0.0; C[1]=0.0;
  for (i=2; i<=n+2; i++) {
    B[i]=A[i-2]-p1*B[i-1]-q*B[i-2];
    C[i]=-B[i-1]-p1*C[i-1]-q*C[i-2];
  }	
//calculate dp=a1 and dq=b1
  x=B[n+1]; y=B[n+2]; z=C[n];
  t=C[n+1]; u=C[n+2];
  d=t*t-z*(u+x);
  if (d==0) goto l1000;
  a1=(z*y-x*t)/d;
  b1=(-x*(q*z+p1*t)-y*t)/d;
//New p1 and q
  p1=p1+a1; q=q+b1;
  f=(fabs(a1)+fabs(b1))/(fabs(p1)+fabs(q));
  if (f>e) goto l200;
//A factor has been found
  Solve2(&x,&y,p1,q,&cont,salida_real,salida_imag);
//Update polynomial
  n=n-2;
  for (i=0; i<=n; i++) A[i]=B[i+2];
  goto l100;
l500: //Last factor, first or second degree
  if (n==2) goto l550;
  x=-A[1]/A[0]; y=0.0; //first degree
  Print_results(x,y,&cont,&salida_real[0],&salida_imag[0]);
  printf("\n");
  return;
l550: p1=A[1]/A[0]; q=A[2]/A[0];
  Solve2(&x,&y,p1,q,&cont,salida_real,salida_imag);
  printf("\n");
  return;
l1000: printf("\n Process not convergent !\n\n");
  for (l=0;l<n;l++){
	salida_imag[l] = -1.0;
	salida_real[l] = -1.0;
  }
}

//==================================Fsolve1===============================================================================
double f1(double x, double Miu_co, double Miu_do, double VdispO, double Vm)
{
	double ViscO = Miu_co;
	double salida;
	salida = ((x/Miu_co)*( pow((((2*x/Miu_co) + (5*Miu_do/Miu_co))/(2 + (5*Miu_do/Miu_co))),(3.0/2.0)))) - ( pow((1 - (VdispO/Vm)),(-2.5)));
	return salida;
}
double df1(double x, double Miu_co, double Miu_do, double VdispO, double Vm)
{
	double ViscO = Miu_co;
	double salida;
	double h=0.001;
	salida = (f1(x+h,Miu_co,Miu_do,VdispO,Vm) - f1(x,Miu_co,Miu_do,VdispO,Vm))/h;
	return salida;
}
double fsolveIvn1(double x0, double Miu_co, double Miu_do, double VdispO, double Vm){
    double x1;
    x0 = Miu_co;
    double allerr = 0.001;
    double maxmitr = 500;
    int itr;
    double h;
    //printf("derivada =%5.15f \n",df(Miu_co,Miu_co,Miu_do,VdispO,Vm));
    for (itr=1; itr<=maxmitr; itr++)
    {
        h=f1(x0,Miu_co,Miu_do,VdispO,Vm)/df1(x0,Miu_co,Miu_do,VdispO,Vm);
        x1=x0-h;
        //printf(" At Iteration no. %3d, x = %9.6f\n", itr, x1);
        if (fabs(h) < allerr)
        {
            //printf("After %3d iterations, root = %8.6f\n", itr, x1);
        }
        x0=x1;
    }
    return x1;
}
//==========================================================================================================================
//==================================Fsolve2===============================================================================
double f2(double x, double Miu_co, double Miu_do, double VdispW, double Vm)
{
	double ViscO = Miu_co;
	double salida;
	salida = ((x/Miu_co)*( pow((((2*x/Miu_co) + (5*Miu_do/Miu_co))/(2 + (5*Miu_do/Miu_co))),(3.0/2.0)))) - ( pow((1 - (VdispW/(1 - Vm))),(-2.5)));
	return salida;
}
double df2(double x, double Miu_co, double Miu_do, double VdispW, double Vm)
{
	double ViscO = Miu_co;
	double salida;
	double h=0.001;
	salida = (f2(x+h,Miu_co,Miu_do,VdispW,Vm) - f2(x,Miu_co,Miu_do,VdispW,Vm))/h;
	return salida;
}
double fsolveIvn2(double x0, double Miu_co, double Miu_do, double VdispW, double Vm){
    double x1;
    x0 = Miu_co;
    double allerr = 0.001;
    double maxmitr = 500;
    int itr;
    double h;
    //printf("derivada =%5.15f \n",df(Miu_co,Miu_co,Miu_do,VdispW,Vm));
    for (itr=1; itr<=maxmitr; itr++)
    {
        h=f2(x0,Miu_co,Miu_do,VdispW,Vm)/df2(x0,Miu_co,Miu_do,VdispW,Vm);
        x1=x0-h;
        //printf(" At Iteration no. %3d, x = %9.6f\n", itr, x1);
        if (fabs(h) < allerr)
        {
            //printf("After %3d iterations, root = %8.6f\n", itr, x1);
        }
        x0=x1;
    }
    return x1;
}
//==========================================================================================================================
//==================================Fsolve3===============================================================================
double f3(double x, double Miu_co, double Miu_do, double VdispO, double Vm)
{
	double ViscO = Miu_co;
	double salida;
	salida = ((x/Miu_co)*( pow((((2*x/Miu_co) + (5*Miu_do/Miu_co))/(2 + (5*Miu_do/Miu_co))),(3.0/2.0)))) - ((9.0/8.0)*(( pow((VdispO/Vm),(1.0/3.0)))/(1 - ( pow((VdispO/Vm),(1.0/3.0))))));
	return salida;
}
double df3(double x, double Miu_co, double Miu_do, double VdispO, double Vm)
{
	double ViscO = Miu_co;
	double salida;
	double h=0.001;
	salida = (f3(x+h,Miu_co,Miu_do,VdispO,Vm) - f3(x,Miu_co,Miu_do,VdispO,Vm))/h;
	return salida;
}
double fsolveIvn3(double x0, double Miu_co, double Miu_do, double VdispO, double Vm){
    double x1;
    x0 = Miu_co;
    double allerr = 0.001;
    double maxmitr = 500;
    int itr;
    double h;
    for (itr=1; itr<=maxmitr; itr++)
    {
        h=f3(x0,Miu_co,Miu_do,VdispO,Vm)/df3(x0,Miu_co,Miu_do,VdispO,Vm);
        x1=x0-h;
        //printf(" At Iteration no. %3d, x = %9.6f\n", itr, x1);
        if (fabs(h) < allerr)
        {
            //printf("After %3d iterations, root = %8.6f\n", itr, x1);
        }
        x0=x1;
    }
    return x1;
}
//==========================================================================================================================
//==================================Fsolve4===============================================================================
double f4(double x, double Miu_co, double Miu_do, double VdispW, double Vm)
{
	double ViscO = Miu_co;
	double salida;
	salida = ((x/Miu_co)*( pow((((2*x/Miu_co) + (5*Miu_do/Miu_co))/(2 + (5*Miu_do/Miu_co))),(3.0/2.0)))) - ((9.0/8.0)*(( pow((VdispW/(1 - Vm)),(1.0/3.0)))/(1 - ( pow((VdispW/(1 - Vm)),(1.0/3.0))))));
	return salida;
}
double df4(double x, double Miu_co, double Miu_do, double VdispW, double Vm)
{
	double salida;
	double h=0.001;
	salida = (f4(x+h,Miu_co,Miu_do,VdispW,Vm) - f4(x,Miu_co,Miu_do,VdispW,Vm))/h;
	return salida;
}
double fsolveIvn4(double x0, double Miu_co, double Miu_do, double VdispW, double Vm){
    double x1;
    x0 = Miu_co;
    double allerr = 0.001;
    double maxmitr = 500;
    int itr;
    double h;
    for (itr=1; itr<=maxmitr; itr++)
    {
        h=f4(x0,Miu_co,Miu_do,VdispW,Vm)/df4(x0,Miu_co,Miu_do,VdispW,Vm);
        x1=x0-h;
        //printf(" At Iteration no. %3d, x = %9.6f\n", itr, x1);
        if (fabs(h) < allerr)
        {
            //printf("After %3d iterations, root = %8.6f\n", itr, x1);
        }
        x0=x1;
    }
    return x1;
}
//==========================================================================================================================
//==================================Fsolve5===============================================================================
double f5(double x, double Miu_co, double Miu_do, double VdispO, double Vm)
{
	double ViscO = Miu_co;
	double salida;
	salida = ((x/Miu_co)*( pow((((2*x/Miu_co) + (5*Miu_do/Miu_co))/(2 + (5*Miu_do/Miu_co))),(3.0/2.0)))) - ((9.0/8.0)*(( pow((VdispO/Vm),(1.0/3.0)))/(1 - ( pow((VdispO/Vm),(1.0/3.0))))));
	return salida;
}
double df5(double x, double Miu_co, double Miu_do, double VdispO, double Vm)
{
	double ViscO = Miu_co;
	double salida;
	double h=0.001;
	salida = (f5(x+h,Miu_co,Miu_do,VdispO,Vm) - f5(x,Miu_co,Miu_do,VdispO,Vm))/h;
	return salida;
}
double fsolveIvn5(double x0, double Miu_co, double Miu_do, double VdispO, double Vm){
    double x1;
    x0 = Miu_co;
    double allerr = 0.001;
    double maxmitr = 500;
    int itr;
    double h;
    for (itr=1; itr<=maxmitr; itr++)
    {
        h=f5(x0,Miu_co,Miu_do,VdispO,Vm)/df5(x0,Miu_co,Miu_do,VdispO,Vm);
        x1=x0-h;
        //printf(" At Iteration no. %3d, x = %9.6f\n", itr, x1);
        if (fabs(h) < allerr)
        {
            //printf("After %3d iterations, root = %8.6f\n", itr, x1);
        }
        x0=x1;
    }
    return x1;
}

//========================================================================================================================
//==================================Fsolve6===============================================================================
double f6(double x, double Miu_co, double Miu_do, double VdispW, double Vm)
{
	double ViscO = Miu_co;
	double salida;
	salida = ((x/Miu_co)*(pow((((2*x/Miu_co) + (5*Miu_do/Miu_co))/(2 + (5*Miu_do/Miu_co))),(3.0/2.0)))) - ((9.0/8.0)*(( pow((VdispW/(1 - Vm)),(1.0/3.0)))/(1 - ( pow((VdispW/(1 - Vm)),(1.0/3.0))))));
	return salida;
}
double df6(double x, double Miu_co, double Miu_do, double VdispW, double Vm)
{
	double ViscO = Miu_co;
	double salida;
	double h=0.001;
	salida = (f6(x+h,Miu_co,Miu_do,VdispW,Vm) - f6(x,Miu_co,Miu_do,VdispW,Vm))/h;
	return salida;
}
double fsolveIvn6(double x0, double Miu_co, double Miu_do, double VdispW, double Vm){
    double x1;
    x0 = Miu_co;
    double allerr = 0.001;
    double maxmitr = 500;
    int itr;
    double h;
    for (itr=1; itr<=maxmitr; itr++)
    {
        h=f6(x0,Miu_co,Miu_do,VdispW,Vm)/df6(x0,Miu_co,Miu_do,VdispW,Vm);
        x1=x0-h;
        //printf(" At Iteration no. %3d, x = %9.6f\n", itr, x1);
        if (fabs(h) < allerr)
        {
            //printf("After %3d iterations, root = %8.6f\n", itr, x1);
        }
        x0=x1;
    }
    return x1;
}

//================================Aparece en inductor.m===================================================================
//==================================Fzero1===============================================================================
double f7(double x, double Re, double RugRel)
{
	double salida;
	salida = ((1/(pow(x,0.5))) + (2*log10((RugRel/3.7) + (2.51/(Re*(pow(x,0.5)))))));
	return salida;
}
double df7(double x, double Re, double RugRel)
{
	double salida;
	double h=0.001;
	salida = (f7(x+h,Re,RugRel) - f7(x,Re,RugRel))/h;
	return salida;
}
double fsolveIvn7(double x0, double Re, double RugRel){
    double x1;
    double allerr = 0.001;
    double maxmitr = 500;
    int itr;
    double h;
    for (itr=1; itr<=maxmitr; itr++)
    {
        h=f7(x0,Re,RugRel)/df7(x0,Re,RugRel);
        x1=x0-h;
        //printf(" At Iteration no. %3d, x = %9.6f\n", itr, x1);
        if (fabs(h) < allerr)
        {
            //printf("After %3d iterations, root = %8.6f\n", itr, x1);
        }
        x0=x1;
    }
    return x1;
}

//================================Aparece en inductor.m===================================================================
//==================================Fzero2===============================================================================
double f8(double x, double CoefCabeza, double Phi, double Cf)
{
	double salida;
	salida = (CoefCabeza*( pow((sin(x)),2))/Phi) + (Cf*Phi) - (Cf*cos(x));
	return salida;
}
double df8(double x, double CoefCabeza, double Phi, double Cf)
{
	double salida;
	double h=0.001;
	salida = (f8(x+h,CoefCabeza,Phi,Cf) - f8(x,CoefCabeza,Phi,Cf))/h;
	return salida;
}
double fsolveIvn8(double x0, double CoefCabeza, double Phi, double Cf){
    double x1;
    double allerr = 0.001;
    double maxmitr = 500;
    int itr;
    double h;
    for (itr=1; itr<=maxmitr; itr++)
    {
        h = f8(x0,CoefCabeza,Phi,Cf)/df8(x0,CoefCabeza,Phi,Cf);
        x1 = x0-h;
        //printf(" At Iteration no. %3d, x = %9.6f\n", itr, x1);
        if (fabs(h) < allerr)
        {
            //printf("After %3d iterations, root = %8.6f\n", itr, x1);
        }
        x0 = x1;
    }
    return x1;
}

void JacobianMezclaFunciones1(double F_Oil, double Do_w, double Dw_o, double C_1, double C_2, double DiametroTubing, double m, double n, double ViscCrudo, double ViscAgua, double We1, double We2, double p, /*Out:*/ double *jaco, double *funEval){
	double a00,a01,a02;
	double a10,a11,a12;
	double a20,a21,a22;
	//df1/dF_Oil
	a00 = (C_1*C_2*DiametroTubing*pow(F_Oil,(m - 1))*m*pow((ViscCrudo/ViscAgua),p))/pow(We1,n);
	//df1/dDo_w
	a01 = -1;
	//df1/dDw_o
	a02 = 0;
	//df2/dF_oil
	a10 = 0;
	//df2/dDo_w
	a11 = 0;
	//df2/dDw_o
	a12 = -1;
	//df3/dF_oil
	a20 = -1;
	//df3/dDo_w
	a21 = 1/(Do_w + Dw_o) - Do_w/pow((Do_w + Dw_o),2);
	//df3/dDw_o
	a22 = -Do_w/pow((Do_w + Dw_o),2);
	
	jaco[0]=a00; jaco[1]=a01; jaco[2]=a02;
	jaco[3]=a10; jaco[4]=a11; jaco[5]=a12;
	jaco[6]=a20; jaco[7]=a21; jaco[8]=a22;
	
	funEval[0] = (C_1*(1 + (C_2*pow(F_Oil,m)))*(pow(We1,-n))*DiametroTubing*pow((ViscCrudo/ViscAgua),p)) - Do_w;
	funEval[1] = (C_1*(pow(We2,-n))*DiametroTubing*pow((ViscAgua/ViscCrudo),p)) - Dw_o;
	funEval[2] = (Do_w/(Dw_o + Do_w)) - F_Oil;
}

void JacobianMezclaFunciones2(double F_Agua, double Do_w, double Dw_o, double C_1, double C_2, double DiametroTubing, double m, double n, double ViscCrudo, double ViscAgua, double We1, double We2, double p, /*Out:*/ double *jaco, double *funEval){
	double a00,a01,a02;
	double a10,a11,a12;
	double a20,a21,a22;
	//df1/dF_Agua
	a00 = 0;
	//df1/dDo_w
	a01 = -1;
	//df1/dDw_o
	a02 = 0;
	
	//df2/dF_Agua
	a10 = (C_1*C_2*DiametroTubing*pow(F_Agua,(m - 1))*m*pow((ViscCrudo/ViscAgua),p))/pow(We1,n);
	//df2/dDo_w
	a11 = 0;
	//df2/dDw_o
	a12 = -1;
	
	//df3/dF_Agua
	a20 = -1;
	//df3/dDo_w
	a21 = 1/(Do_w + Dw_o) - Do_w/pow((Do_w + Dw_o),2);
	//df3/dDw_o
	a22 = -Do_w/pow((Do_w + Dw_o),2);
	
	jaco[0]=a00; jaco[1]=a01; jaco[2]=a02;
	jaco[3]=a10; jaco[4]=a11; jaco[5]=a12;
	jaco[6]=a20; jaco[7]=a21; jaco[8]=a22;
	
	funEval[0] = (C_1*DiametroTubing*pow((ViscAgua/ViscCrudo),p))/pow(We2,n) - Do_w;
	funEval[1] = (C_1*DiametroTubing*(C_2*pow(F_Agua,m) + 1)*pow((ViscCrudo/ViscAgua),p))/pow(We1,n) - Dw_o;
	funEval[2] = Do_w/(Do_w + Dw_o) - F_Agua;
}

//Invierte matrix y la guarda en la misma variable
void inversaJacobian3x3(double *matrix){
	double temp[9];
	double det;
	memcpy(&temp[0],matrix,3*3*sizeof(double));
	
	//Determinante
	det = temp[0]*temp[4]*temp[8] + temp[1]*temp[5]*temp[6] + temp[2]*temp[3]*temp[7] - temp[2]*temp[4]*temp[6] - temp[1]*temp[3]*temp[8] - temp[0]*temp[5]*temp[7];
	
	matrix[0] = temp[4]*temp[8]-temp[7]*temp[5];   matrix[1] = -1*(temp[1]*temp[8]-temp[7]*temp[2]);  matrix[2] = temp[1]*temp[5]-temp[4]*temp[2];
	matrix[3] = -1*(temp[3]*temp[8]-temp[6]*temp[5]);   matrix[4] = temp[0]*temp[8]-temp[6]*temp[2];  matrix[5] = -1*(temp[0]*temp[5]-temp[3]*temp[2]);
	matrix[6] = temp[3]*temp[7]-temp[6]*temp[4];   matrix[7] = -1*(temp[0]*temp[7]-temp[6]*temp[1]);  matrix[8] = temp[0]*temp[4]-temp[3]*temp[1];
	
	matrix[0] = matrix[0]/(det);  matrix[1] = matrix[1]/(det);  matrix[2] = matrix[2]/(det);
	matrix[3] = matrix[3]/(det);  matrix[4] = matrix[4]/(det);  matrix[5] = matrix[5]/(det);
	matrix[6] = matrix[6]/(det);  matrix[7] = matrix[7]/(det);  matrix[8] = matrix[8]/(det);
	
}


void solveIvn1(int itmax, double tol, double F_Oil, double Do_w, double Dw_o, double C_1, double C_2, double DiametroTubing, double m, double n, double ViscCrudo, double ViscAgua, double We1, double We2, double p, /*Out: */double *solucion){
	double F_Oil1,Do_w1,Dw_o1;
	double F_Oil2,Do_w2,Dw_o2;
	int i;
	double jacobian[9];
	double funEval[3];
	
	F_Oil1 = F_Oil;
	Do_w1 = Do_w;
	Dw_o1 =Dw_o;
	
	for (i=0;i<itmax;i=i+1){
		JacobianMezclaFunciones1(F_Oil1,Do_w1,Dw_o1,C_1,C_2,DiametroTubing,m,n,ViscCrudo,ViscAgua,We1,We2,p, /*Out:*/ &jacobian[0], &funEval[0]);
		inversaJacobian3x3(&jacobian[0]);
		F_Oil2 = F_Oil1-1*( jacobian[0]*funEval[0] + jacobian[1]*funEval[1] + jacobian[2]*funEval[2] );
		Do_w2 = Do_w1-1*( jacobian[3]*funEval[0] + jacobian[4]*funEval[1] + jacobian[5]*funEval[2] );
		Dw_o2 = Dw_o1-1*( jacobian[6]*funEval[0] + jacobian[7]*funEval[1] + jacobian[8]*funEval[2] );
		
		if (fabs((F_Oil2 - F_Oil1)/F_Oil2) < tol){
			F_Oil1 = F_Oil2;
			Do_w1 = Do_w2;
			Dw_o1 = Dw_o2;
			break; 
		}
		F_Oil1 = F_Oil2;
		Do_w1 = Do_w2;
		Dw_o1 = Dw_o2;
	}
	
	solucion[0] = F_Oil1;
	solucion[1] = Do_w1;
	solucion[2] = Dw_o1;	
}

void solveIvn2(int itmax, double tol, double F_Agua, double Do_w, double Dw_o, double C_1, double C_2, double DiametroTubing, double m, double n, double ViscCrudo, double ViscAgua, double We1, double We2, double p, /*Out: */double *solucion){
	double F_Agua1,Do_w1,Dw_o1;
	double F_Agua2,Do_w2,Dw_o2;
	int i;
	double jacobian[9];
	double funEval[3];
	
	F_Agua1 = F_Agua;
	Do_w1 = Do_w;
	Dw_o1 =Dw_o;
	
	for (i=0;i<itmax;i=i+1){
		JacobianMezclaFunciones2(F_Agua1,Do_w1,Dw_o1,C_1,C_2,DiametroTubing,m,n,ViscCrudo,ViscAgua,We1,We2,p, /*Out:*/ &jacobian[0], &funEval[0]);
		inversaJacobian3x3(&jacobian[0]);
		F_Agua2 = F_Agua1-1*( jacobian[0]*funEval[0] + jacobian[1]*funEval[1] + jacobian[2]*funEval[2] );
		Do_w2 = Do_w1-1*( jacobian[3]*funEval[0] + jacobian[4]*funEval[1] + jacobian[5]*funEval[2] );
		Dw_o2 = Dw_o1-1*( jacobian[6]*funEval[0] + jacobian[7]*funEval[1] + jacobian[8]*funEval[2] );
		
		if (fabs((F_Agua2 - F_Agua1)/F_Agua2) < tol){
			F_Agua1 = F_Agua2;
			Do_w1 = Do_w2;
			Dw_o1 = Dw_o2;
			break; 
		}
		F_Agua1 = F_Agua2;
		Do_w1 = Do_w2;
		Dw_o1 = Dw_o2;
	}
	
	solucion[0] = F_Agua1;
	solucion[1] = Do_w1;
	solucion[2] = Dw_o1;	
}