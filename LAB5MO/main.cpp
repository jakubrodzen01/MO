#include <iostream>
#include <cmath>
#define SIZE 4

using namespace std;

void print_matrix(double matrix[][SIZE]) 
{
	for (int i = 0; i < SIZE; i++) 
	{
		for (int j = 0; j < SIZE; j++)
			cout << matrix[i][j] << "\t";
		cout << endl;
	}
}

void print_vector(double vector[]) 
{
	for (int i = 0; i < SIZE; i++)
		cout << vector[i] << "\t";
	cout << endl;
}

void swapRows(double table[][SIZE], int index, int i) 
{
	for (int j = 0; j < SIZE; j++)
		swap(table[index][j], table[i][j]);
}

void LU(double matrix[][SIZE], double matrixL[][SIZE], int index[]) 
{
	for (int i = 0; i < SIZE; i++)			//wykonuje sie tyle ile wierszy 
	{
		if (matrix[i][i] == 0.0)				//sprawdza czy wartość do zerowania nie jest już równa zero 
		{
			double max = 0.0;
			int newIndex = 0;
			for (int j = i; j < SIZE; j++) 
			{
				if (fabs(matrix[j][i]) > max && matrix[j][i] != 0)	//szuka w kolejnych wierszach najwiekszej wartosci bezwzglednej
				{
					newIndex = j;
					max = matrix[j][i];
				}
			}
			swapRows(matrixL, index[i], newIndex);		//zamienia ten wiersz z elementem zerowym i znaleziony wiersz w tablicy L
			swapRows(matrix, index[i], newIndex);		//zamienia ten wiersz z elementem zerowym i znaleziony wiersz w tablicy A
			swap(index[newIndex], index[i]);			
		}

		matrixL[i][i] = 1.0;

		for (int j = i + 1; j < SIZE; j++)				//rozklad gaussa
		{
			matrixL[j][i] = matrix[j][i] / matrix[i][i];		//wspolczynnik
			for (int k = 0; k <= SIZE; k++)
				matrix[j][k] -= (matrixL[j][i] * matrix[i][k]);	
		}
	}
}

void solveY(double matrix[][SIZE], double result[], double vector[])		//Ly=b
{
	for (int i = 0; i < SIZE; i++) 
	{
		result[i] = vector[i];
		double sum = 0.0;
		for (int j = 0; j < i; j++)
			sum += result[j] * matrix[i][j];
		result[i] -= sum;
	}
}
void solveX(double matrix[][SIZE], double result[], double vector[])		//Ux=y
{
	for (int i = SIZE - 1; i >= 0; i--) 
	{
		result[i] = vector[i];
		double sum = 0.0;
		for (int j = SIZE - 1; j > i; j--)
			result[i] -= result[j] * matrix[i][j];
		result[i] /= matrix[i][i];
	}
}
int main() 
{
	double A[SIZE][SIZE]{
	{1.0, -20.0, 30.0, -4.0},
	{2.0, -40.0, -6.0, 50.0},
	{9.0, -180.0, 11.0, -12.0},
	{-16.0, 15.0, -140.0, 13.0} };

	double B[SIZE]{ 35.0, 104.0, -366.0, -354.0 };

	int index_in_vector[SIZE]{ 0, 1, 2, 3 };

	double L[SIZE][SIZE] = { 0.0 };
	double Y[SIZE]{ 0.0 };
	double X[SIZE]{ 0.0 };

	cout << "Macierz A: " << endl;
	print_matrix(A);
	cout << endl << "Wektor B: " << endl;
	print_vector(B);

	LU(A, L, index_in_vector);

	for (int i = 0; i < SIZE - 1; i++)		//dopasowuje kolejnosc elementow w B do nowych indeksow
		swap(B[i], B[index_in_vector[i]]);

	cout << endl << "Macierz L po dekompozycji LU: " << endl;
	print_matrix(L);
	cout << endl << "Macierz U (A) po dekompozycji LU: " << endl;
	print_matrix(A);
	cout << endl << "Wektor B po zamianie wierszy: " << endl;
	print_vector(B);

	solveY(L, Y, B);
	solveX(A, X, Y);

	cout << endl << "Wektor Y: " << endl;
	print_vector(Y);
	cout << endl << "Wektor X: " << endl;
	print_vector(X);
}
