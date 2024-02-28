#include <iostream>
#include <cfloat>

using namespace std;

int main() {
  int mantysa = 0;
  double d_epsilon = 1.0;
  float f_epsilon = 1.0f;

  while(1.f + (f_epsilon / 2.0f) > 1.0f) {
    f_epsilon /= 2.0f;
    mantysa++;
  }
  cout << "Dla float: " << endl;
  cout << "Mantysa = " << mantysa << endl;
  cout << "Epsilon = " << f_epsilon << endl;

  mantysa = 0;

  while(1 + (d_epsilon / 2.0) > 1.0) {
    d_epsilon /= 2.0;
    mantysa++;
  }
  cout << "Dla double: " << endl;
  cout << "Mantysa = " << mantysa << endl;
  cout << "Epsilon = " << d_epsilon << endl;

  cout << endl;
  cout << FLT_EPSILON << endl;
  cout << DBL_EPSILON << endl;

  return 0;
}
