#include<iostream>
#include<cmath>

using namespace std;

#define TOLX 1e-12
#define TOLF 1e-12
#define ITERATIONS 50

double f1(double x, double y, double z)		{return x * x + y * y + z * z - 2.0;}
double f2(double x, double y)				{return x * x + y * y - 1.0;}
double f3(double x, double y)				{return x * x - y;}

double d1(double x, double y)				{return (2.0 * x * x * y + x * x - y * y - 1.0) / (4.0 * x * y + 2.0 * x);}
double d2(double y)							{return (y * y + y - 1.0) / (2.0 * y + 1.0);}
double d3(double z)							{return (z * z - 1.0) / (2.0 * z);}

int main() 
{
	double est1, est2, est3, res1, res2, res3;
	//punkty poczatkowe
	double ax = 5.0, ay = 7.0, az = 14.0;
	double ax1, ay1, az1;
	int iter = 1;

	printf("i: 0 | x: %f | y: %f | z: %f\n", ax, ay, az);

	while (iter <= ITERATIONS) 
	{
		ax1 = ax - d1(ax, ay);
		ay1 = ay - d2(ay);
		az1 = az - d3(az);

		est1 = fabs(ax1 - ax);
		est2 = fabs(ay1 - ay);
		est3 = fabs(az1 - az);
		res1 = fabs(f1(ax1, ay1, az1));
		res2 = fabs(f2(ax1, ay1));
		res3 = fabs(f3(ax1, ay1));
		
		printf("i: %1d | x: %f | y: %f | z: %f | est(x): %.4e | est(y): %.4e | est(z): %.4e | f1(xn): %.4e | f2(xn): %.4e | f3(xn): %.4e \n", iter, ax1, ay1, az1, est1, est2, est3, res1, res2, res3); 
		
		if ((est1 <= TOLX) && (est2 <= TOLX) && (est3 <= TOLX) && (res1 <= TOLF) && (res2 <= TOLF) && (res3 <= TOLF))
		{
			cout << "Osiagnieto kryterium dokladnosci TOLX!" << endl;
			cout << "Osiagnieto kryterium wiarygodnosci TOLF!" << endl;
			return 0;
		}

		/*
		if ((res1 <= TOLF) && (res2 <= TOLF) && (res3 <= TOLF))
		{
			cout << "Osiagnieto kryterium wiarygodnosci TOLF!" << endl;
			return 0;
		}
		*/

		iter++;
		ax = ax1;
		ay = ay1;
		az = az1;
	}
	cout << "Osiagnieto limit iteracji" << endl;
}