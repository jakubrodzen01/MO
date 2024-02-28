#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

double BME(double, double, double, double);
double PME(double, double, double, double);
double PMT(double, double, double, double);

void get_actual(double, double, double);
void errors(double, double, double, double, double);

double actual(double t) {
	return 1.0 - exp(-10.0 * (t + atan(t)));
}

double f_pom(double t) {
	return (10.0 * t * t + 20.0) / (t * t + 1.0);
}

double f(double t, double yt) {
	return (1.0 - yt) * f_pom(t);
}

int main() {
	double delta_t = 0.001;
	double start = 0.0, finish = 1.0;
	double cond = 0.0;

	get_actual(delta_t, start, finish);
	BME(delta_t, start, finish, cond);
	PME(delta_t, start, finish, cond);
	PMT(delta_t, start, finish, cond);

	double lim = 1e-10;

	errors(delta_t, start, finish, cond, lim);
	return 0;
}

void get_actual(double delta_t, double start, double finish) {
	ofstream actual_values;
	actual_values.open("wartosci_dokladne.dat");
	while (start < finish) {
		actual_values << start << " " << actual(start) << endl;
		start += delta_t;
	}
}

double BME(double delta_t, double start, double finish, double yk) {
	ofstream BME_val;
	BME_val.open("BME.dat");
	double yk1;
	double maks_err = 0.0, local_err;
	while (start <= finish) {
		BME_val << start << " " << yk << endl;
		yk1 = yk + f(start, yk) * delta_t;
		start += delta_t;
		yk = yk1;
		local_err = abs(yk - actual(start));
		if (maks_err < local_err) {
			maks_err = local_err;
		}
	}
	return maks_err;
}

double PME(double delta_t, double start, double finish, double yk) {
	ofstream PME_val;
	PME_val.open("PME.dat");
	double yk1;
	double maks_err = 0.0, local_err;
	while (start <= finish) {
		PME_val << start << " " << yk << endl;
		yk1 = (yk + f_pom(start + delta_t) * delta_t) / (1.0 + f_pom(start + delta_t) * delta_t);
		start += delta_t;
		yk = yk1;
		local_err = abs(yk - actual(start));
		if (maks_err < local_err) {
			maks_err = local_err;
		}
	}
	return maks_err;
}

double PMT(double delta_t, double start, double finish, double yk) {
	ofstream PMT_val;
	PMT_val.open("PMT.dat");
	double yk1;
	double maks_err = 0.0, local_err;
	while (start <= finish) {
		PMT_val << start << " " << yk << endl;
		yk1 = (yk + (delta_t / 2.0 * f_pom(start)) + (delta_t / 2.0 * f_pom(start + delta_t)) - yk * delta_t / 2.0 * f_pom(start)) / ((delta_t / 2.0) * f_pom(start + delta_t) + 1.0);
		start += delta_t;
		yk = yk1;
		local_err = abs(yk - actual(start));
		if (maks_err < local_err) {
			maks_err = local_err;
		}
	}
	return maks_err;
}

void errors(double delta_t, double start, double finish, double cond, double lim) {
	ofstream BME_err, PME_err, PMT_err;
	BME_err.open("BME_err.dat");
	PME_err.open("PME_err.dat");
	PMT_err.open("PMT_err.dat");
	double BME_err_val, PME_err_val, PMT_err_val;
	while (delta_t >= lim) {
		BME_err_val = BME(delta_t, start, finish, cond);
		PME_err_val = PME(delta_t, start, finish, cond);
		PMT_err_val = PMT(delta_t, start, finish, cond);
		BME_err << log10(delta_t) << " " << log10(BME_err_val) << endl;
		PME_err << log10(delta_t) << " " << log10(PME_err_val) << endl;
		PMT_err << log10(delta_t) << " " << log10(PMT_err_val) << endl;
		delta_t /= 2.0;
	}
}
