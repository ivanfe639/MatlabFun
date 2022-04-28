#ifndef funciones_Matlab
#define funciones_Matlab

	double sumIvn(double *matriz, int tam);

	void zerosIvn(double *vector, int filas, int columnas);

	void linspaceIvn(double a,double b, int N,double *salida);

	double minIvn(double *vector, int tam);

	double maxIvn(double *vector, int tam);

	void splineCoefIvn(double *x, double *a, int tam, double *coef);

	void evalCoefSpline(double *vec, int n_vec, double *coef, double *x, int n_x, /*Out: */ double *salida);

	int polyfitIvn(double* XX_temp, double* YY,int countOfElements,  int order, double* coefficients, int Normalizar, double MU0, double MU1);

	double polyvalIvn(double *vecPoly, int tamPoly, double X, double MU0, double MU1);

	double meanIvn(double *vector, int tamVector);

	double stdIvn(double *vector, int tamVector);

	void getMatrixSizeIvn(int *filas, int *columnas, char *nombreBin);

	void getMatrixIvn(double *vector, char *nombreBin,int tam);

	void roots3Ivn(double *polinomio, double *raices_real, double *raices_imag);

	void Print_results(double x, double y, int *cont, double *salida_real, double *salida_imag);

	void Solve2(double *x, double *y, double p1, double q, int *cont, double *salida_real, double *salida_imag);

	void rootsIvnAll(double *poli, int n, double *salida_real, double *salida_imag);

	double f1(double x, double Miu_co, double Miu_do, double VdispO, double Vm);

	double df1(double x, double Miu_co, double Miu_do, double VdispO, double Vm);

	double fsolveIvn1(double x0, double Miu_co, double Miu_do, double VdispO, double Vm);

	double f2(double x, double Miu_co, double Miu_do, double VdispW, double Vm);

	double df2(double x, double Miu_co, double Miu_do, double VdispW, double Vm);

	double fsolveIvn2(double x0, double Miu_co, double Miu_do, double VdispW, double Vm);

	double f3(double x, double Miu_co, double Miu_do, double VdispO, double Vm);

	double df3(double x, double Miu_co, double Miu_do, double VdispO, double Vm);

	double fsolveIvn3(double x0, double Miu_co, double Miu_do, double VdispO, double Vm);

	double f4(double x, double Miu_co, double Miu_do, double VdispW, double Vm);

	double df4(double x, double Miu_co, double Miu_do, double VdispW, double Vm);

	double fsolveIvn4(double x0, double Miu_co, double Miu_do, double VdispW, double Vm);

	double f5(double x, double Miu_co, double Miu_do, double VdispO, double Vm);

	double df5(double x, double Miu_co, double Miu_do, double VdispO, double Vm);

	double fsolveIvn5(double x0, double Miu_co, double Miu_do, double VdispO, double Vm);

	double f6(double x, double Miu_co, double Miu_do, double VdispW, double Vm);

	double df6(double x, double Miu_co, double Miu_do, double VdispW, double Vm);

	double fsolveIvn6(double x0, double Miu_co, double Miu_do, double VdispW, double Vm);

	double f7(double x, double Re, double RugRel);

	double df7(double x, double Re, double RugRel);

	double fsolveIvn7(double x0, double Re, double RugRel);

	double f8(double x, double CoefCabeza, double Phi, double Cf);

	double df8(double x, double CoefCabeza, double Phi, double Cf);

	double fsolveIvn8(double x0, double CoefCabeza, double Phi, double Cf);

	void JacobianMezclaFunciones1(double F_Oil, double Do_w, double Dw_o, double C_1, double C_2, double DiametroTubing, double m, double n, double ViscCrudo, double ViscAgua, double We1, double We2, double p, /*Out:*/ double *jaco, double *funEval);

	void JacobianMezclaFunciones2(double F_Agua, double Do_w, double Dw_o, double C_1, double C_2, double DiametroTubing, double m, double n, double ViscCrudo, double ViscAgua, double We1, double We2, double p, /*Out:*/ double *jaco, double *funEval);

	void inversaJacobian3x3(double *matrix);

	void solveIvn1(int itmax, double tol, double F_Oil, double Do_w, double Dw_o, double C_1, double C_2, double DiametroTubing, double m, double n, double ViscCrudo, double ViscAgua, double We1, double We2, double p, /*Out: */double *solucion);

	void solveIvn2(int itmax, double tol, double F_Agua, double Do_w, double Dw_o, double C_1, double C_2, double DiametroTubing, double m, double n, double ViscCrudo, double ViscAgua, double We1, double We2, double p, /*Out: */double *solucion);
	
#endif