#include <iostream>
#include <fstream>
#include <vector>


using namespace std;

class condition {
public: float xi, yi, omega;

condition(float y, float x, float omega): yi(y), xi(x), omega(omega) {}


condition evolution (float dt) {
    condition cj (yi - dt * xi * omega, xi + dt * yi, omega);
    return cj;
}

};


void data(int n, float x0, float y0, float omega, float dt) {

    condition c0 (x0, y0, omega);

    vector<condition> d;

    d.push_back(c0);
    for (int i = 1; i < n; ++i) {
        condition c1 = c0.evolution(dt);
        d.push_back(c1);
        c0 = c1;
        // condition::evolution(d(i-1), omega, deltat);
    }

    ofstream fout("03_10.txt");
    for (int i = 0; i < n; ++i) {
        fout << d[i].xi << " " << d[i].yi << endl;
    }
    fout.close();

}

int main() {
    long n;
    cin >> n;
    float x0, y0, omega, dt;
    cin >> x0 >> y0 >> omega >> dt;

    data (n, x0, y0, omega, dt);

}



