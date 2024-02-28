#include <iostream>
#include <cmath>

#define TOLX 10e-8
#define TOLF 10e-8
#define ITERATIONS 30

int liczba_iteracji = 1;

using namespace std;

double fun_a(double x) { return sin(x / 4.0) * sin(x / 4.0) - x; }
double fun_b(double x) { return tan(2.0 * x) - 1.0 - x; }
double fun_a_fi(double x) { return sin(x / 4.0) * sin(x / 4.0); }
double fun_b_fi(double x) { return tan(2.0 * x) - 1.0; }	//fi(x) jest rozbiezna
double fun_b_fi_pic(double x) { return atan(x+1)/2; }
double fun_a_der(double x) { return 1.0 / 4.0 * sin(x / 2.0) - 1.0; }
double fun_b_der(double x) { return -1.0 + (2.0 / (cos(2.0 * x) * cos(2.0 * x))); }
//double fun_a_fi_der(double x) { return sin(x / 2.0) / 4.0; }
//double fun_b_fi_der(double x) { return 2.0 / (cos(2.0 * x) * cos(2.0 * x)); }


double picard(double (*fun)(double), double (*fun_fi)(double), double x);
double bisekcja(double (*fun)(double), double a, double b);
double newton(double (*fun)(double), double (*fun_der)(double), double x);
double sieczne(double (*fun)(double), double a, double b);


int main()
{
	cout << "Metoda picarda dla fun_a:" << endl;
	picard(fun_a, fun_a_fi, 3.0);
	cout << endl;
	liczba_iteracji = 1;

	cout << "Metoda picarda dla fun_b:" << endl;
	picard(fun_b, fun_b_fi, 3.0);
	cout << endl;
	liczba_iteracji = 1;

	cout << "Metoda bisekcji dla fun_a:" << endl;
	bisekcja(fun_a, -1.0, 5.0);
	cout << endl;
	liczba_iteracji = 1;

	cout << "Metoda bisekcji dla fun_b:" << endl;
	bisekcja(fun_b, 0.0, 1.0);
	cout << endl;
	liczba_iteracji = 1;

	cout << "Metoda Newtona dla fun_a:" << endl;
	newton(fun_a, fun_a_der, 3.0);
	cout << endl;
	liczba_iteracji = 1;

	cout << "Metoda Newtona dla fun_b:" << endl;
	newton(fun_b, fun_b_der, 3.0);
	cout << endl;
	liczba_iteracji = 1;

	cout << "Metoda siecznych dla fun_a:" << endl;
	sieczne(fun_a, 3.0, 2.5);
	cout << endl;
	liczba_iteracji = 1;

	cout << "Metoda siecznych dla fun_b:" << endl;
	sieczne(fun_b, 30.0, 20.5);
	cout << endl;
	liczba_iteracji = 1;

	return 0;
}

double picard(double (*fun)(double), double (*fun_fi)(double), double x)
{
	double xn = fun_fi(x);

	if (fabs(xn) >= 1.0)
	{
		printf("Nie jest zbiezna\n");
		return -1;
	}

	double estymator = xn - x;
	printf("%d.\txn=%-12.6g\te.bledu=%-12.6g\tresiduum=%-12.6g\n", liczba_iteracji, xn, fabs(estymator), fabs(fun(xn)));
	
	if (fabs(estymator) <= TOLX && fabs(fun(xn)) <= TOLF)
	{
		printf("Osiagnieto tolerancje bledu %-12.6g\n", fabs(estymator));
		printf("Osiagnieto tolerancje residuum %-12.6g\n", fabs(fun(xn)));
		return xn;
	}
	/*
	if (fabs(fun(xn)) <= TOLF)
	{
		printf("Osiagnieto tolerancje residuum %-12.6g\n", fabs(fun(xn)));
		return xn;
	}
	*/
	if (liczba_iteracji < ITERATIONS)
	{
		liczba_iteracji++;
		return picard(fun, fun_fi, xn);
	}

	printf("Osiagnieto limit iteracji %d\n", liczba_iteracji);
	return xn;
}

double bisekcja(double (*fun)(double), double a, double b)
{
	double xn = (a + b) / 2.0;
	double x = fun(xn);
	double estymator = (b - a) / 2.0;			//polowa dlugosci podprzedzialu
	printf("%d.\txn=%-12.6g\ta=%-12.6g\tb=%-12.6g\te.bledu=%-12.6g\tresiduum=%-12.6g\n", liczba_iteracji, xn, a, b, fabs(estymator), fabs(x));
	
	if (fabs(estymator) <= TOLX && fabs(x) <= TOLF)
	{
		printf("Osiagnieto tolerancje bledu %-12.6g\n", fabs(estymator));
		printf("Osiagnieto tolerancje residuum %-12.6g\n", fabs(x));
		return xn;
	}
	/*
	if (fabs(x) <= TOLF)
	{
		printf("Osiagnieto tolerancje residuum %-12.6g\n", fabs(fun_a(x)));
		return xn;
	}
	*/
	if (liczba_iteracji < ITERATIONS)
	{
		if (fun(a) * x > 0.0)			//maja ten sam znak to xn nowy poczatek przedzialu 
			a = xn;
		else 
			b = xn;
		liczba_iteracji++;
		return bisekcja(fun, a, b);
	}

	printf("Osiagnieto limit iteracji %d\n", liczba_iteracji);
	return xn;
}

double newton(double (*fun)(double), double (*fun_der)(double), double x)
{
	if (fun_der(x) == 0.0)
	{
		printf("Nie jest zbiezna");
		return -1;
	}

	double xn = x - ((fun(x)) / (fun_der(x)));
	double estymator = xn - x;
	printf("%d.\txn=%-12.6g\te.bledu=%-12.6g\tresiduum=%-12.6g\n", liczba_iteracji, xn, fabs(estymator), fabs(fun(xn)));

	if (fabs(estymator) <= TOLX && fabs(fun(xn)) <= TOLF)
	{
		printf("Osiagnieto tolerancje bledu %-12.6g\n", fabs(estymator));
		printf("Osiagnieto tolerancje residuum %-12.6g\n", fabs(fun(xn)));
		return xn;
	}
	/*
	if (fabs(xn) <= TOLF)
	{
		printf("Osiagnieto tolerancje residuum %-12.6g\n", fabs(xn));
		return xn;
	}
	*/
	if (liczba_iteracji < ITERATIONS)
	{
		liczba_iteracji++;
		return newton(fun, fun_der, xn);
	}

	printf("Osiagnieto limit iteracji %d\n", liczba_iteracji);
	return xn;
}

double sieczne(double (*fun)(double), double a, double b)
{
	double xn = b - (fun(b) / (((fun(b)) - (fun(a))) / (b - a)));
	double estymator = xn - b;
	printf("%d.\txn+2=%-12.6g\txn+1=%-12.6g\txn=%-12.6g\te.bledu=%-12.6g\tresiduum = % -12.6g\n", liczba_iteracji, xn, b, a, fabs(estymator), fabs(fun(xn)));
	
	if (fabs(estymator) <= TOLX && fabs(fun(xn)) <= TOLF)
	{
		printf("Osiagnieto tolerancje bledu %-12.6g\n", fabs(estymator));
		printf("Osiagnieto tolerancje residuum %-12.6g\n", fabs(fun(xn)));
		return xn;
	}
	/*
	if (fabs(xn) <= TOLF)
	{
		printf("Osiagnieto tolerancje residuum %-12.6g\n", fabs(xn));
		return xn;
	}
	*/
	if (liczba_iteracji < ITERATIONS)
	{
		liczba_iteracji++;
		return sieczne(fun, b, xn);
	}

	printf("Osiagnieto limit iteracji %d\n", liczba_iteracji);
	return xn;
}
