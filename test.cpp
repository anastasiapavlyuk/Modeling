#include <iostream>
#include <fstream>
#include <vector>


using namespace std;

class State{
public:
    float x, v;
};

class Solver {
public: 
    float xi, yi, omega;

    State state;


    Solver(float y, float x, float omega): yi(y), xi(x), omega(omega) {
        state.v = y;
        state.x = x;
    }


    Solver evolution (float dt) {
        // Solver cj (yi - dt * xi * omega, xi + dt * yi, omega);
        // return cj;
        return Solver {yi - dt * xi * omega, xi + dt * yi, omega};
    }

    void euler_step (float dt){
        float v = state.v - dt * state.x * omega;
        float x = state.x + dt * state.v;
        this->state.x = x;
        this->state.v = v;
    }

};


void solve(int n, float x0, float y0, float omega, float dt) {

    Solver c0 (x0, y0, omega);

    vector<State> d;

    d.push_back(c0.state);
    for (int i = 1; i < n; ++i) {
        // condition c1 = c0.evolution(dt);
        c0.euler_step(dt);
        d.push_back(c0.state);
        // c0 = c1;
        // condition::evolution(d(i-1), omega, deltat);
    }

    {
        ofstream fout("03_10.txt");
        for (auto &el: d) {
        // int i = 0; i < n; ++i) {
            fout << el.x << " " << el.v << endl;
        }
        fout.close();
    }
}

#include <array>
int main(int argc, char *argv[]) {
    long n;
    n = 10000;
    // cin >> n;
    float x0, y0, omega, dt;
    x0 = 64e5;
    y0 = 0.0;
    omega = 1;
    dt = 1.0e-3;

    if (argc == 9){
        std::cout << "using cmd arguments:" << std::endl;
        for(int el=0; el<9; el++){
            std::cout << argv[el]<< ", ";
        }
        for(int el=1; el<9; el+=2){
            if(std::string(argv[el])=="-x"){
                x0 = atof(argv[el+1]);
            }
            if(std::string(argv[el])=="-v"){
                y0 = atof(argv[el+1]);
            }
            if(std::string(argv[el])=="-n"){
                n = atol(argv[el+1]);
            }
            if(std::string(argv[el])=="-dt"){
                dt = atof(argv[el+1]);
            }
        }
    }
    else {
        std::cout << "using default initial state as  number of arguments != 9" << std::endl;

    }
    // cin >> x0 >> y0 >> omega >> dt;


    solve (n, x0, y0, omega, dt);

}



