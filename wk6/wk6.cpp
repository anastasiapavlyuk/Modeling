#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <cmath>


using namespace std;

class condition {
public: double xi, yi, omega;

condition(double y, double x, double omega): yi(y), xi(x), omega(omega) {}

condition operator+ (condition other) {
    return condition (this->yi+other.yi, this->xi+other.xi, this->omega);
}

condition operator* (double a) {
    return condition (a * this->yi, a * this->xi, this->omega);
}
};

condition f(condition ci, double x) {
    condition cj( - ci.xi * ci.omega, ci.yi, ci.omega);
    return cj;
}

class solver {

public: 
    long n;
    double y0, x0, omega, dt;
    vector <condition> data;

    solver(long n, double y0, double x0, double omega, double dt):
    n(n), y0(y0), x0(x0), omega(omega), dt(dt) {

        data.push_back(condition(y0, x0, omega));

        for (long i = 1; i < n; ++i) {
            condition c1 = data[i-1] + (f(data[i-1], dt) * dt);
            condition c2 = data[i-1] + ((f(data[i-1], dt) + f(c1, dt))*(dt/2));
            data.push_back(c2);
        }
    }

    /*
    void write() {
        ofstream fout(this->name);
        fout << "v x max dt n \n";
        for (long i = 0; i < n; ++i) {
            if (data[i].xi >= data[i-1].xi && data[i].xi >= data[i+1].xi) {
                fout << data[i].yi << " " << data[i].xi << " " << data[i].xi <<" " << dt << " " << n << endl;
            }
            else {
                fout << data[i].yi << " " << data[i].xi << " " << 0 <<" " << dt << " " << n << endl;
            }
        }
        fout.close();
    }



    void clear(){
        this->data.clear();
    }
    */

};

int main() {
    int n = 10;

    long long Time[n];

    for (int i = 0; i < n; ++i){
        auto startTime = chrono::high_resolution_clock::now();
        solver c1(pow(10, i), 0, 6400000, 1, 0.1);
        auto endTime = chrono::high_resolution_clock::now();
        Time[i] = chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count();
        cout << " " << Time[i] << endl;
    }

    return 0;
}

