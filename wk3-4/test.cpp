#include <iostream>
#include <fstream>
#include <vector>
#include <string>
//#include "/Users/anastasiapavluk/json/single_include/nlohmann/json.hpp"

//using json = nlohmann::json;

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
    string type, name;
    vector <condition> data;

    solver(long n, double y0, double x0, double omega, double dt, string type, string name):
    n(n), y0(y0), x0(x0), omega(omega), dt(dt), type(type), name(name) {

        data.push_back(condition(y0, x0, omega));

        if (type == "euler") {
            for (long i = 1; i < n; ++i) {
                data.push_back(data[i-1]+f(data[i-1], dt)*dt);
            }
        }

        if (type == "heun") {
            for (long i = 1; i < n; ++i) {
                condition c1 = data[i-1] + (f(data[i-1], dt) * dt);
                condition c2 = data[i-1] + ((f(data[i-1], dt) + f(c1, dt))*(dt/2));
                data.push_back(c2);
            }
        }

        if (type == "rg45") {
            for (long i = 1; i < n; ++i) {
                condition k1 = f(data[i-1], dt);
                condition k2 = f(data[i-1] + (k1 * (dt/2)), dt);
                condition k3 = f(data[i-1] + (k2 * (dt/2)), dt);
                condition k4 = f(data[i-1] + (k3 * dt), dt);
                condition c = data[i-1] + (k1 + (k2*2) + (k3*2) + k4) * (dt/6);
                data.push_back(c);
            }
        }
    }

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
};


int main() {
    double x0 = 6400000, y0 = 0, omega = 1, h = 0.1, dt = 0.1;

    vector <solver> Data;

    Data.push_back(solver(2000, 0, 6400000, 1, 0.1, "euler", "e6400000_0_01_2000.txt"));
    Data.push_back(solver(500, 0, 6400000, 1, 0.1, "euler", "e6400000_0_01_500.txt"));
    Data.push_back(solver(2000, 0, 6400000, 1, 0.5, "euler", "e6400000_0_05_2000.txt"));
    Data.push_back(solver(2000, 0, 6400000, 1, 0.01, "euler", "e6400000_0_001_2000.txt"));

    Data.push_back(solver(2000, 0, 6400000, 1, 0.1, "heun", "h6400000_0_01_2000.txt"));
    Data.push_back(solver(100000, 0, 6400000, 1, 0.1, "heun", "h6400000_0_01_100000.txt"));
    Data.push_back(solver(2000, 0, 6400000, 1, 0.5, "heun", "h6400000_0_05_2000.txt"));
    Data.push_back(solver(2000, 0, 6400000, 1, 0.01, "heun", "h6400000_0_001_2000.txt"));

    Data.push_back(solver(20000, 0, 6400000, 1, 0.001, "euler", "e6400000_0_0001_20000.txt"));
    Data.push_back(solver(20000, 0, 6400000, 1, 0.005, "euler", "e6400000_0_0005_20000.txt"));

    for (double i = 0.01; i <= 0.08;  i += 0.01) {
        string s = to_string(i).substr(0, 1+to_string(i).find_last_not_of('0')), s1="";
        for (unsigned j = 0; j < s.size(); ++j) if (s[j] != '.') s1 += s[j];
        Data.push_back(solver(20000, 0, 6400000, 1, i, "euler", "e6400000_0_"+s1+"_20000.txt"));
        cout << "e6400000_0_"+s1+"_20000.txt" << endl;
    }

    for (double i = 0.01; i <= 0.09;  i += 0.01) {
        string s = to_string(i).substr(0, 1+to_string(i).find_last_not_of('0')), s1="";
        for (unsigned j = 0; j < s.size(); ++j) if (s[j] != '.') s1 += s[j];
        Data.push_back(solver(20000, 0, 6400000, 1, i, "heun", "h6400000_0_"+s1+"_20000.txt"));
    }

    for (double i = 0.1; i <= 0.4;  i += 0.1) {
        string s = to_string(i).substr(0, 1+to_string(i).find_last_not_of('0')), s1="";
        for (unsigned j = 0; j < s.size(); ++j) if (s[j] != '.') s1 += s[j];
        Data.push_back(solver(20000, 0, 6400000, 1, i, "heun", "h6400000_0_"+s1+"_20000.txt"));
    }

    for (double i = 0.01; i <= 0.09;  i += 0.01) {
        string s = to_string(i).substr(0, 1+to_string(i).find_last_not_of('0')), s1="";
        for (unsigned j = 0; j < s.size(); ++j) if (s[j] != '.') s1 += s[j];
        Data.push_back(solver(20000, 0, 6400000, 1, i, "rg45", "rg6400000_0_"+s1+"_20000.txt"));
    }
    

    for (int i = 0; i < Data.size(); ++i) {
        Data[i].write();
    }
/*
    double d = 0.01;
    string s = to_string(d);
    cout << s << endl;

    cout << s.substr(0, s.find_last_not_of('0')+1) << endl;
*/
    //cout << config.dump(4);
    return 0;
}



