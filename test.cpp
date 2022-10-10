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

condition heunev (float h) {
    condition inter (yi - h * xi * omega, xi + h * yi, omega);
    condition cj (yi - 0.5 * h * (xi*omega + inter.xi*omega), xi + 0.5 * h * (yi + inter.yi), omega);
    return cj;
}

};


void data(int n, float x0, float y0, float omega, float dt, string name) {

    condition c0 (x0, y0, omega);

    vector<condition> d;

    d.push_back(c0);
    for (int i = 1; i < n; ++i) {
        condition c1 = c0.evolution(dt);
        d.push_back(c1);
        c0 = c1;
    }

    vector<float> max;
    
    for (int i=1; i < n; ++i) {
        if (d[i].yi >= d[i-1].yi & d[i].yi >= d[i+1].yi) {
            max.push_back(d[i].yi/x0);
        }
        else {
            max.push_back(0);
        }
    }

    ofstream fout(name);
    fout << "x v max dx n \n";
    for (int i = 0; i < n; ++i) {
        fout << d[i].xi << " " << d[i].yi << " " << max[i] <<" " << dt << " " << n << endl;
    }

    fout.close();
}

void dataheun(int n, float x0, float y0, float omega, float h, string name) {

    condition c0 (x0, y0, omega);

    vector<condition> d;

    d.push_back(c0);
    for (int i = 1; i < n; ++i) {
        condition c1 = c0.heunev(h);
        d.push_back(c1);
        c0 = c1;
    }

    vector<float> max;
    max.push_back(1);
    for (int i=1; i < n; ++i) {
        if (d[i].yi >= d[i-1].yi && d[i].yi >= d[i+1].yi) {
            max.push_back(d[i].yi/x0);
        }
        else {
            max.push_back(0);
        }
    }

    ofstream fout(name);
    fout << "v x max h n \n";
    for (int i = 0; i < n; ++i) {
        fout << d[i].xi << " " << d[i].yi << " " << max[i] << " " << h << " " << n << endl;
    }

    fout.close();
}




int main() {

    float x0 = 6400000, y0 = 0, omega = 1, h = 0.1, dt = 0.1;

    data(20000, x0, y0, omega, 0.001, "e6400000_0_0001_20000.txt");
    data(20000, x0, y0, omega, 0.005, "e6400000_0_0005_20000.txt");
    data(20000, x0, y0, omega, 0.01,  "e6400000_0_001_20000.txt"); 
    data(20000, x0, y0, omega, 0.02,  "e6400000_0_002_20000.txt"); 
    data(20000, x0, y0, omega, 0.03,  "e6400000_0_003_20000.txt"); 
    data(20000, x0, y0, omega, 0.04,  "e6400000_0_004_20000.txt"); 
    data(20000, x0, y0, omega, 0.05,  "e6400000_0_005_20000.txt"); 
    data(20000, x0, y0, omega, 0.06,  "e6400000_0_006_20000.txt"); 
    data(20000, x0, y0, omega, 0.07,  "e6400000_0_007_20000.txt"); 
    data(20000, x0, y0, omega, 0.08,  "e6400000_0_008_20000.txt"); 
    data(20000, x0, y0, omega, 0.09,  "e6400000_0_009_20000.txt"); 
    data(20000, x0, y0, omega, 0.1,   "e6400000_0_01_20000.txt"); 
    data(20000, x0, y0, omega, 0.2,   "e6400000_0_02_20000.txt"); 
    data(20000, x0, y0, omega, 0.3,   "e6400000_0_03_20000.txt"); 
    data(20000, x0, y0, omega, 0.4,   "e6400000_0_04_20000.txt"); 
    data(20000, x0, y0, omega, 0.5,   "e6400000_0_05_20000.txt");
    data(20000, x0, y0, omega, 0.6,   "e6400000_0_06_20000.txt");
    data(20000, x0, y0, omega, 0.7,   "e6400000_0_07_20000.txt");
    data(20000, x0, y0, omega, 0.8,   "e6400000_0_08_20000.txt");
    data(20000, x0, y0, omega, 0.9,   "e6400000_0_09_20000.txt");
    data(20000, x0, y0, omega, 1,     "e6400000_0_1_20000.txt");
    

    dataheun(20000, x0, y0, omega, 0.001, "h6400000_0_0001_20000.txt");
    dataheun(20000, x0, y0, omega, 0.005, "h6400000_0_0005_20000.txt");
    dataheun(20000, x0, y0, omega, 0.01,  "h6400000_0_001_20000.txt"); 
    dataheun(20000, x0, y0, omega, 0.02,  "h6400000_0_002_20000.txt"); 
    dataheun(20000, x0, y0, omega, 0.03,  "h6400000_0_003_20000.txt"); 
    dataheun(20000, x0, y0, omega, 0.04,  "h6400000_0_004_20000.txt"); 
    dataheun(20000, x0, y0, omega, 0.05,  "h6400000_0_005_20000.txt"); 
    dataheun(20000, x0, y0, omega, 0.06,  "h6400000_0_006_20000.txt"); 
    dataheun(20000, x0, y0, omega, 0.07,  "h6400000_0_007_20000.txt"); 
    dataheun(20000, x0, y0, omega, 0.08,  "h6400000_0_008_20000.txt"); 
    dataheun(20000, x0, y0, omega, 0.09,  "h6400000_0_009_20000.txt"); 
    dataheun(20000, x0, y0, omega, 0.1,   "h6400000_0_01_20000.txt"); 
    dataheun(20000, x0, y0, omega, 0.2,   "h6400000_0_02_20000.txt"); 
    dataheun(20000, x0, y0, omega, 0.3,   "h6400000_0_03_20000.txt"); 
    dataheun(20000, x0, y0, omega, 0.4,   "h6400000_0_04_20000.txt"); 
    dataheun(20000, x0, y0, omega, 0.5,   "h6400000_0_05_20000.txt");
    dataheun(20000, x0, y0, omega, 0.6,   "h6400000_0_06_20000.txt");
    dataheun(20000, x0, y0, omega, 0.7,   "h6400000_0_07_20000.txt");
    dataheun(20000, x0, y0, omega, 0.8,   "h6400000_0_08_20000.txt");
    dataheun(20000, x0, y0, omega, 0.9,   "h6400000_0_09_20000.txt");
    dataheun(20000, x0, y0, omega, 1,     "h6400000_0_1_20000.txt");
   
    return 0;
}



