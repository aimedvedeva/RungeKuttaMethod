#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;
#define e 2.718281828459045235
ofstream file;


double fun2(double x, double y) {
	return powl(y, 2)*powl(e, x) - 2 * y;
}
double fun(double x, double y) {
	return (1 / x) * (-powl(x, 5) * powl(y, 3) * powl(e, x) - 2 * y);
}

double testFun(double x, double y) {
	return x + y;
}

double *GetEvenGrid(int numOfDots, double a, double b, double *delta) {
	*delta = (b - a) / (double)(numOfDots);
	double *grid = new double[numOfDots + 1];
	grid[0] = a;
	for (int i = 1; i < numOfDots + 1; i++) {
		grid[i] = grid[i - 1] + (*delta);
	}
	return grid;
}

// Runge-Kutta method
double GetR(double yH, double yh, double r) {
	return abs(yH - yh) / (powl(2, r - 1) - 1);
}

double *RungeKuttaMethod(double a, double b, double x0, double y0, int n, double(*fun) (double, double), double **grid) {
	double *funGrid = new double[n + 1];
	funGrid[0] = y0;

	double delta;
	*grid = GetEvenGrid(n, a, b, &delta);

	for (int i = 1; i < 4; i++) {
		double k0, k1, k2, k3;
		k0 = fun((*grid)[i - 1], funGrid[i - 1]);
		k1 = fun((*grid)[i - 1] + delta / (double)2, funGrid[i - 1] + k0 * delta / (double)2);
		k2 = fun((*grid)[i - 1] + delta / (double)2, funGrid[i - 1] + k1 * delta / (double)2);
		k3 = fun((*grid)[i - 1] + delta, funGrid[i - 1] + k2 * delta);

		funGrid[i] = funGrid[i - 1] + (delta / (double)6) * (k0 + 2 * k1 + 2 * k2 + k3);
	}

	return funGrid;
}

double *RungeKuttaMethodE(double a, double b, double x0, double y0, int *n, double(*fun) (double, double), double **grid, double E, double r, double *delta) {
	int hn = *n * 2;
	double *yh = NULL;
	double *xh = NULL; 
	double R;
	bool isMudified = false;

	for (int i = 1; i < *n + 1; i++) {
		yh = RungeKuttaMethod(a, b, x0, y0, hn, fun, &xh);

        R = GetR(yh[i * 2], yh[i], r);

		if (R > E){
			isMudified = true;
			hn *= 2;
		}

		if (i == *n) {
			if (isMudified) {
				i = 1;
				isMudified = false;
			}
			else {
				break;
			}
		}
		delete[] xh;
		delete[] yh;
	}
	*n = hn;
	*grid = xh;
	 *delta = (b - a) / *n;
	return yh;
}

double *AdamsMethod(double a, double b, int *n, double(*fun) (double, double), double *E /* out */ , double x0, double y0, double **grid) {
	double num = (double)19 / (double)270;

	double delta = 0.1;

	double *funGrid = RungeKuttaMethod(a, b, x0, y0, *n, fun, grid);

	double corrector, predictor;

	for (int i = 4; i < *n + 1; i++) {

	   predictor = funGrid[i - 1] + (delta / 24) * (55 * fun((*grid)[i - 1], funGrid[i - 1]) 
			                                              - 59 * fun((*grid)[i - 2], funGrid[i - 2]) 
			                                              + 37 * fun((*grid)[i - 3], funGrid[i - 3]) 
														  - 9  * fun((*grid)[i - 4], funGrid[i - 4]));
	
	   //carry out several iterations to improve accuracy
		double tmp = predictor;
		
		for (int j = 0; j < 2; j++) {
			corrector = funGrid[i - 1] + (delta / 24) * ( 9 * fun((*grid)[i], tmp)
                                        			   + 19 * fun((*grid)[i - 1], funGrid[i - 1])
													   - 5  * fun((*grid)[i - 2], funGrid[i - 2])
													   +      fun((*grid)[i - 3], funGrid[i - 3]));
			tmp = corrector;
		}
		
		funGrid[i] = corrector;
     }

	*E = 19.0 / 270.0 * abs(predictor - corrector);
	return funGrid;
}

void PrintResult(double *funGrid, double *grid, int n) {
	file << n + 1 << endl;

	for (int i = 0; i < n + 1; i++) {
		file << funGrid[i] << endl;
	}

	file << endl;

	for (int i = 0; i < n + 1; i++) {
		file << grid[i] << endl;
	}

}

double EuclideanNorm(double *vector, int lenght) {
	double sum = 0;
	for (int i = 0; i < lenght; i++) {
		sum += powl(vector[i], 2);
	}
	sum = sqrt(sum);
	return sum;
}

int main() {
	/*double a = 0;
	double b = 1;
	double x0 = 0;
	double y0 = 0;

	
	//number of partition gaps --> H = 2 && h = 4
	int n = 9;
	double E = 0;
	double *grid;
	double delta;

	double *funGrid = AdamsMethod(a,b, &n, testFun, &E, x0, y0, &grid);*/
	double a = 0;
	double b = 1;
	double x0 = 0;
	double y0 = 1;
	int n = 10;
	double *grid = NULL;
	double delta = 0;
	double E = 0;
	double *funGrid = AdamsMethod(a,b, &n, testFun, &E, x0, y0, &grid);

	file.open("output.txt");
	PrintResult(funGrid, grid, n);

	/*cout << E << endl;
	x0 = 0;
	y0 = 1 + 0.01;
	double *funGrid1 = AdamsMethod(a, b, &n, testFun, &E, x0, y0, &grid);

	for (int i = 0; i < n + 1; i++) {
		funGrid[i] -= funGrid1[i];
	}

	cout << EuclideanNorm(funGrid, n + 1) << endl;*/
	/*
	cout << E << endl;
	file.open("output.txt");
	PrintResult(funGrid, grid, n);
	*/
	delete[] grid;
	delete[] funGrid;
	file.close();
}

