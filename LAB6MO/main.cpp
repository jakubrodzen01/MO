#include <iostream>
#include <cmath>
#include <iomanip>
#define SIZE 6
using namespace std;

void print_vector(double* vector) 
{
	for (int i = 0; i < SIZE; i++)
		cout << setw(12) << vector[i];
	cout << endl;
}

void print_matrix(double* lower, double* diagonal, double* upper) 
{
	int l_counter = 0, u_counter = 0;
	for (int i = 0; i < SIZE; i++) 
	{
		for (int j = 0; j < SIZE; j++) 
		{
			if (i == j)
				cout << setw(12) << diagonal[i];
			else if (i - j == 1)
				cout << setw(12) << lower[l_counter++];
			else if (i - j == -1)
				cout << setw(12) << upper[u_counter++];
			else
				cout << setw(12) << 0.0;
		}
		cout << endl;
	}
}
void thomas_in_matrix(double* lower, double* upper, double* eta)	//L U D
{
	int l_counter = 0;
	int u_counter = 1;
	for (int i = 1; i < SIZE; i++)
		eta[i] -= lower[l_counter++] * upper[u_counter++ - 1] / eta[i - 1];		//eta(i) = d(i) - l(i)/eta(i-1)u(i-1)
}
void thomas_in_vector(double* lower, double* eta, double* r)		//L D b
{
	int l_counter = 0;
	for (int i = 1; i < SIZE; i++)
		r[i] -= lower[l_counter++] * r[i - 1] / eta[i - 1];						//r(i) = b(i) - l(i)/eta(i-1)r(i-1)
}
void solve_x(double* eta, double* upper, double* r, double* x)		//D, U, b, x
{
	for (int i = SIZE - 1; i >= 0; i--) 
	{
		if (i == SIZE - 1)
			x[i] = r[i] / eta[i];
		else
			x[i] = (r[i] - (upper[i] * x[i + 1])) / eta[i];
	}
}
int main() 
{
	double AL[SIZE - 1]{ 1.0 / 3.0, 1.0 / 5.0, 1.0 / 7.0, 1.0 / 9.0, 1.0 / 11.0 };
	double AD[SIZE]{ 10.0, 20.0, 30.0, 30.0, 20.0, 10.0 };
	double AU[SIZE - 1]{ 1.0 / 2.0, 1.0 / 4.0, 1.0 / 6.0, 1.0 / 8.0, 1.0 / 10.0 };
	double b[SIZE]{ 31.0, 165.0 / 4.0, 917.0 / 30.0, 851.0 / 28.0, 3637.0 / 90.0, 332.0 / 11.0 };
	double x[SIZE]{ 0.0, };

	cout << "Macierz A: \n";
	print_matrix(AL, AD, AU); 
	cout << endl << "\nWektor B: \n";
	print_vector(b);

	thomas_in_matrix(AL, AU, AD);
	thomas_in_vector(AL, AD, b);

	solve_x(AD, AU, b, x);

	cout << "\n========= Po obliczeniach =========" << endl << "Wektor eta: \n";
	print_vector(AD);
	cout << "\nWektor U: " << endl;
	print_vector(AU);
	cout << "\nWektor B: " << endl;
	print_vector(b);
	cout << "\n\nWynik obliczen - Wektor X: \n";
	print_vector(x);
	return 0;
}
