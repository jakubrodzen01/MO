#include <iostream>
#include <Windows.h>
#include <math.h>
#include <iomanip>

#define warunek_stopu 1.0e-12
#define max_iter 100
#define N 4

using namespace std;

double estymator(double* x, double* x2);
double residuum(double A[][N], double* b, double* x2);
void wypiszWektor(double*);
void wypiszMacierz(double A[][N]);
void MetodaGaussaSeidel(double A[][N], double*, double*);
void MetodaSOR(double A[][N], double*, double*, double omega);
void MetodaJacobiego(double A[][N], double*, double*);

int main()
{
	double omega = 0.5;
	double A[N][N]{
	{100, -1, 2, -3},
	{1, 200, -4, 5},
	{-2, 4, 300, -6},
	{3, -5, 6, 400} };
	double b[N]{ 116, -226, 912, -1174 };
	double x[N]{ 2,2,2,2 };
	wypiszMacierz(A);

	cout << endl;
	cout << "____________________________________________________________________________________________________________________________________________"<<endl;
	cout << " Wektor b:" << endl;
	wypiszWektor(b);
	cout << endl;
	cout << "____________________________________________________________________________________________________________________________________________"<<endl;
	cout << " Wektor x (przyblizenie poczatkowe):" << endl;
	wypiszWektor(x);
	cout << endl;
	cout << "____________________________________________________________________________________________________________________________________________"<<endl;
	//MetodaJacobiego(A, b, x);
	//MetodaGaussaSeidel(A, b, x);
	MetodaSOR(A, b, x, omega);
	return 0;
}
void wypiszMacierz(double A[][N])
{
	cout << "Macierz A:" << endl;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			cout << left << setw(5) << A[i][j];
		}
		cout << endl;
	}
}
void wypiszWektor(double* c)
{
	for (int i = 0; i < N; i++)
	{
		cout << c[i] << endl;
	}
}
double estymator(double* x, double* x2)
{
	double ESTYMATOR = 0.0;
	double maks_EST = 0.0;
	for (int i = 0; i < N; i++)
	{
		ESTYMATOR = fabs(x2[i] - x[i]);
		if (ESTYMATOR > maks_EST)
		{
			maks_EST = ESTYMATOR;
		}
		ESTYMATOR = 0.0;
	}
	return maks_EST;
}
// R(0) = b - A * x(0), dla kazdego wiersza 
double residuum(double A[][N], double* b, double* x)
{
	double RESIDUUM = 0.0;
	double maks_RES = 0.0;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			RESIDUUM += A[i][j] * x[j];
		}
		RESIDUUM = fabs(b[i] - RESIDUUM);
		if (RESIDUUM > maks_RES)
		{
			maks_RES = RESIDUUM;
		}

		RESIDUUM = 0.0;
	}
	return maks_RES;
}
void MetodaJacobiego(double A[][N], double* b, double* x)
{
	double suma = 0.0;
	double EST = 0.0;
	double RESIDUUM = 0.0;
	double* NowyX = new double[N];

	cout << endl << endl << "Metoda MetodaJacobiego" << endl;
	cout << left << setw(15) << "Iteracja" << setw(20) << "x0" << setw(20) << "x1" << setw(20) << "x2"<<setw(20)<<setw(20)<<"x3"<<setw(20)<<"ESTYMATOR"<<setw(20)<<"RESIDUUM"<< endl << endl;
		for (int ITERACJA = 0; ITERACJA < max_iter; ITERACJA++)
		{
			for (int i = 0; i < N; i++)
			{
				suma = 0.0;
				for (int j = 0; j < N; j++)
				{
					if (j != i)
						suma += A[i][j] * x[j];		//(L+U) * x
				}
				NowyX[i] = (b[i] - suma) / A[i][i];
			}
			EST = estymator(x, NowyX);
			RESIDUUM = residuum(A, b, NowyX);
			for (int i = 0; i < N; i++)
			{
				x[i] = NowyX[i];
			}

			cout << left << setw(15) << ITERACJA << setw(20) << NowyX[0] << setw(20) << NowyX[1] << setw(20) << NowyX[2] << setw(20) << NowyX[3] << setw(20) << EST << setw(20) << RESIDUUM << endl;
			if (EST < warunek_stopu && RESIDUUM < warunek_stopu)
			{
				break;
			}
		}
}
//metoda gaussa seidela nowe wartosci x przypisuje natychmiast po ich znalezieniu zastepuja stare i sia uzywane do obliczenia xi+1, inaczej jest z metoda jacobiego
//gdzie przyblizenia skladowych rozwiazania sa wykorzystywane dopiero w nastepnej iteracji i mozna je obliczac jednoczesnie
void MetodaGaussaSeidel(double A[][N], double* b, double* x)
{
	double* PoprzedniX = new double[N];
	double suma = 0.0;
	double EST = 0.0, RESIDUUM = 0.0;
	cout << endl << endl << "Metoda Gaussa-Seidela" << endl;
	cout << left << setw(15) << "Iteracja" << setw(20) << "x0" << setw(20) << "x1" << setw(20) << "x2"<<setw(20)<<setw(20)<<"x3"<<setw(20)<<"ESTYMATOR"<<setw(20)<<"RESIDUUM"<< endl << endl;
		for (int ITERACJA = 0; ITERACJA < max_iter; ITERACJA++)
		{
			for (int i = 0; i < N; i++)
			{
				suma = 0.0;
				for (int j = 0; j < N; j++)
				{
					if (j != i)
						suma += A[i][j] * x[j];			
				}
				PoprzedniX[i] = x[i];
				x[i] = (b[i] - suma) / A[i][i];
			}
			EST = estymator(PoprzedniX, x);
			RESIDUUM = residuum(A, b, x);
			cout << left << setw(15) << ITERACJA << setw(20) << x[0] << setw(20) << x[1] << setw(20) << x
				[2] << setw(20) << x[3] << setw(20) << EST << setw(20) << RESIDUUM << endl;
			if (EST < warunek_stopu && RESIDUUM < warunek_stopu)
			{
				break;
			}

		}
}
void MetodaSOR(double A[][N], double* b, double* x, double omega)
{
	double EST = 0.0;
	double RESIDUUM = 0.0;
	double* NowyX = new double[N];
	double* PoprzedniX = new double[N];
	double suma = 0.0;

	cout << endl << endl << "Metoda MetodaSOR" << endl;
	cout << left << setw(15) << "Iteracja" << setw(20) << "x0" << setw(20) << "x1" << setw(20) << "x2"<<setw(20)<<setw(20)<<"x3"<<setw(20)<<"ESTYMATOR"<<setw(20)<<"RESIDUUM"<< endl << endl;
		for (int ITERACJA = 0; ITERACJA < max_iter; ITERACJA++)
		{
			for (int i = 0; i < N; i++)
			{
				suma = 0.0;
				for (int j = 0; j < N; j++)
					if (j != i)
						suma += A[i][j] * x[j];
				PoprzedniX[i] = x[i];
				NowyX[i] = (1.0 - omega) * x[i] + (omega / A[i][i]) * (b[i] - suma);
				x[i] = NowyX[i];
			}
			EST = estymator(PoprzedniX, NowyX);
			RESIDUUM = residuum(A, b, NowyX);
			cout << left << setw(15) << ITERACJA << setw(20) << NowyX[0] << setw(20) << NowyX[1] << setw(20) << NowyX[2] << setw(20) << NowyX[3] << setw(20) << EST << setw(20) << RESIDUUM << endl;
			if (EST < warunek_stopu && RESIDUUM < warunek_stopu)
			{
				break;
			}
		}
}
