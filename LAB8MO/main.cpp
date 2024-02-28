#include <iostream>
#include <math.h>
#include <fstream>

#define PI 3.141592653589793238462643383279
#define ILOSC 100
#define WYNIKI 10

using namespace std;

template <typename T> void zapisz(T** tab, string nazwaPliku);
template<typename T> inline T fx(T x) {
	// sin(x)
	return sin(x);
}
template<typename T> inline T dfx(T x) {
	//pochodna sin(x)
	return cos(x);
}
template <typename T> inline T progresywna2(T x, T h) {
	return (fx(x + h) - fx(x)) / h;
}
template <typename T> inline T centralna2(T x, T h) {
	return (fx(x + h) - fx(x - h)) / (static_cast <T>(2.0) * h);
}
template <typename T> inline T wsteczna2(T x, T h) {
	return (fx(x) - fx(x - h)) / h;;
}
template <typename T> inline T progresywna3(T x, T h) {
	return (static_cast <T>(-3.0) / static_cast <T>(2.0) * fx(x) + static_cast <T>(2.0)
		* fx(x + h) - static_cast <T>(1.0) / static_cast <T>(2.0) * fx(x + h + h)) / (h);
}
template <typename T> inline T wsteczna3(T x, T h) {
	return (static_cast <T>(3.0) / static_cast <T>(2.0) * fx(x) - static_cast <T>(2.0)
		* fx(x - h) + static_cast <T>(1.0) / static_cast <T>(2.0) * fx(x - h - h)) / (h);
}
template <class T> void licz(string nazwaPliku) {
	T** blad = new T * [ILOSC];
	T h = 0.2; //pierwszy krok
	T poczatek = 0;
	T koniec = PI / static_cast <T>(2.0);
	T srodek = (poczatek + koniec) / static_cast <T>(2.0);
	for (int i = 0; i < ILOSC; i++) {
		blad[i] = new T[WYNIKI];
		blad[i][0] = log10(h);
		//lewa
		blad[i][1] = log10(fabs(progresywna2(poczatek, h) - dfx(poczatek)));
		blad[i][2] = log10(fabs(progresywna3(poczatek, h) - dfx(poczatek)));
		//srodek
		blad[i][3] = log10(fabs(progresywna2(srodek, h) - dfx(srodek)));
		blad[i][4] = log10(fabs(progresywna3(srodek, h) - dfx(srodek)));
		blad[i][5] = log10(fabs(centralna2(srodek, h) - dfx(srodek)));
		blad[i][6] = log10(fabs(wsteczna2(srodek, h) - dfx(srodek)));
		blad[i][7] = log10(fabs(wsteczna3(srodek, h) - dfx(srodek)));
		//prawa
		blad[i][8] = log10(fabs(wsteczna2(koniec, h) - dfx(koniec)));
		blad[i][9] = log10(fabs(wsteczna3(koniec, h) - dfx(koniec)));
		h = h / static_cast <T>(1.2); //zmniejszamy krok
	}
	zapisz(blad, nazwaPliku);
}
template <typename T> void zapisz(T** tab, string nazwaPliku) {
	fstream plik;
	plik << scientific;
	plik.precision(16);
	plik.open(nazwaPliku, ios::out);
	if (!plik.good()) {
		cout << "Blad otwarcia pliku " << endl;
		return;
	}
	for (int i = 0; i < ILOSC; i++) {
		plik << tab[i][0] << " ";
		for (int j = 1; j < WYNIKI; j++) {
			if (!isinf(tab[i][j])) {
				plik << tab[i][j] << " ";
			}
			else {
				plik << 0.0 << " ";
			}
		}
		plik << endl;
	}
	plik.close();
	plik.open("rzedy_dokladnosci_" + nazwaPliku, ios::out);
	plik << "Rzedy dokladnosci: " << endl;
	for (int i = 1; i < 10; i++) {
		plik << (tab[10][i] - tab[1][i]) / (tab[10][0] - tab[1][0]) << endl;
	}
	plik.close();
}
int main()
{
	cout << scientific;
	cout.precision(10);
	licz<float>("float.txt");
	licz<double>("double.txt");
}
