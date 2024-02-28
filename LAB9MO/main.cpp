#include <iostream>
#include <cmath>
#include <fstream>
#include <Windows.h>
#include <iomanip>

using namespace std;

const double p = 1.0, q = 0.0, r = -4.0, fi = 0.0, psi = 1.0, teta = 0.0, przedzial_poczatek = 0.0, przedzial_koniec = 1.0, alfa = 0.0, beta = 1.0, gamma = -1.0;

double wzorAnalityczny(double x);
double s(double x);
void algorytmThomasa(double* l, double* d, double* u, double* b, double* x, int N);
double maksymalnyBlad(double* blad, int N);
double dyskretyzacjaNumerowa(double h, int N);
double dyskretyzacjaTrzypunktowa(double h, int N);

int main()
{
	fstream bledy, dane;
	double h, dkt, dn;
	bledy.open("bledy.txt", fstream::out);
	for (int N = 10; N < 100000; N += 10) {
		h = (przedzial_koniec - przedzial_poczatek) / (N - 1);
		dkt = log10(dyskretyzacjaTrzypunktowa(h, N));
		dn = log10(dyskretyzacjaNumerowa(h, N));
		bledy << log10(h) << " " << dkt << " " << dn << endl;
	}
	bledy.close();
	system("pause");
	return 0;
}

void algorytmThomasa(double* l, double* d, double* u, double* b, double* x, int N)
{
	double* tb = new double[N];
	tb[0] = b[0];
	double* td = new double[N];
	td[0] = d[0];


	for (int i = 1; i < N; i++) {
		td[i] = d[i] - l[i - 1] * (u[i - 1] / td[i - 1]);
	}
	for (int i = 1; i < N; i++) {
		tb[i] = b[i] - l[i - 1] * tb[i - 1] / td[i - 1];
	}

	x[N - 1] = tb[N - 1] / td[N - 1];
	for (int i = N - 2; i >= 0; i--) {
		x[i] = (tb[i] - u[i] * x[i + 1]) / td[i];
	}
	delete[] tb;
	delete[] td;
}

double wzorAnalityczny(double x)
{
	return (exp(2 - 2 * x) - 4 * exp(4 - 2 * x) + 4 * exp(2 * x) - exp(2 + 2 *x) - x + x * exp(4)) / (4 - 4 * exp(4));
}

double s(double x)
{
	return -x;
}

double maksymalnyBlad(double* blad, int N)
{
	double maks_blad = fabs(blad[0]);
	for (int i = 0; i < N; i++) {
		if (fabs(blad[i]) > maks_blad) {
			maks_blad = fabs(blad[i]);
		}
	}
	return maks_blad;
}

double dyskretyzacjaNumerowa(double h, int N)
{
	double* blad;
	double* l, * d, * u, * b, * x;
	double xpb = przedzial_poczatek, xp = przedzial_poczatek;
	l = new double[N];
	d = new double[N];
	u = new double[N];
	b = new double[N];
	x = new double[N];
	blad = new double[N];
	//postać wektorowo-macierzowa
	u[0] = alfa / h;
	d[0] = beta - alfa / h;
	b[0] = -gamma;
	for (int i = 1; i < N - 1; i++) {
		l[i - 1] = p / (h * h) + r / 12.0;
		d[i] = (-2.0 * p) / (h * h) + r * (10.0 / 12.0);
		u[i] = p / (h * h) + r / 12.0;
		b[i] = -s(xpb + i * h - h) / 12.0 - (10.0 / 12.0) * s(xpb + i * h) - s(xpb + i * h + h) / 12.0;
	}
	l[N - 2] = -fi / h;
	d[N - 1] = -fi / h + psi;
	b[N - 1] = -teta;
	algorytmThomasa(l, d, u, b, x, N);
	//zapis do pliku i wypisanie w konsoli
	for (int i = 0; i < N; i++)
	{
		blad[i] = fabs(x[i] - wzorAnalityczny(xpb));
		xpb += h;
	}
	if (N == 20) {
		fstream plik;
		plik.open("wynikiNumerow.txt", fstream::out);
		cout << "Wyniki metody Numerowa:" << endl;
		cout << left << setw(15) << "i" << setw(20) << "x" << setw(20) << "Un" << setw(20) << "U(x)" << endl;
		for (int i = 0; i < N; i++)
		{
			plik << xp << " " << x[i] << " " << wzorAnalityczny(xp) << " " << endl;
			cout << left << setw(15) << i << setw(20) << xp << setw(20) << x[i] << setw(20) << wzorAnalityczny(xp) << endl;
			xp += h;
		}
		plik.close();
	}
	double maxblad = maksymalnyBlad(blad, N);
	delete[] l;
	delete[] d;
	delete[] u;
	delete[] x;
	delete[] b;
	delete[] blad;
	return maxblad;
}
double dyskretyzacjaTrzypunktowa(double h, int N)
{
	double* blad;
	double* l, * d, * u, * b, * x;
	double xpb = przedzial_poczatek, xp = przedzial_poczatek;
	l = new double[N];
	d = new double[N];
	u = new double[N];
	b = new double[N];
	x = new double[N];
	blad = new double[N];
	//postać wektorowo-macierzowa
	u[0] = alfa / h;
	d[0] = beta - alfa / h;
	b[0] = -gamma;
	for (int i = 1; i < N - 1; i++) {
		l[i - 1] = p / (h * h) - q / (2.0 * h);
		d[i] = (-2.0 * p) / (h * h) + r;
		u[i] = p / (h * h) + q / (2.0 * h);
		b[i] = -s(xpb + i * h);
	}
	l[N - 2] = -fi / h;
	d[N - 1] = -fi / h + psi;
	b[N - 1] = -teta;
	algorytmThomasa(l, d, u, b, x, N);
	for (int i = 0; i < N; i++) {
		blad[i] = fabs(x[i] - wzorAnalityczny(xpb));
		xpb += h;
	}
	//zapis do pliku i wypisanie wartosci w konsoli
	if (N == 20) {
		fstream plik;
		plik.open("wyniki_konwencjonalna.txt", fstream::out);
		cout << "Wyniki metody konwencjonalnej trzypunktowej:" << endl;
		cout << left << setw(15) << "i" << setw(20) << "x" << setw(20) << "Un" << setw(20) << "U(x)" << endl;
		for (int i = 0; i < N; i++)
		{
			blad[i] = fabs(x[i] - wzorAnalityczny(xp));
			plik << xp << " " << x[i] << " " << wzorAnalityczny(xp) << " " << endl;
			cout << left << setw(15) << i << setw(20) << xp << setw(20) << x[i] << setw(20) << wzorAnalityczny(xp) << endl;
			xp += h;
		}
		cout << endl;
		plik.close();
	}
	double maxblad = maksymalnyBlad(blad, N);
	delete[] l;
	delete[] d;
	delete[] u;
	delete[] x;
	delete[] b;
	delete[] blad;
	return maxblad;
}
