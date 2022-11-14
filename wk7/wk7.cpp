#include <iostream>
#include <array>
#include <cmath>
#include <fstream>
#include <vector>

class State {
public:

    std::array<double, 2> state;
    double w;

    State(): w(0) {
        state[0] = 0;
        state[1] = 0;
    }
    State (double x, double v, double w): w(w) {
        state[0] = x;
        state[1] = v;
    }
    
    State (const std::array<double, 2> s, double w): w(w), state(s) {}

    double & operator [](int i) {
        return state[i];
    }

    State operator + (State other) {
        return State ((*this)[0]+other[0], (*this)[1]+other[1], (*this).w);
    }

     State operator * (double a) {
        return State ((*this)[0]*a, (*this)[1]*a, (*this).w);
     }
};

template <class S>
class Evolution {
public:
    Evolution() {}
    S pp(S s) {
        return S (s[1], -s.w*sin(s[0]), s.w);
    }
};

template <class S, class Ev>
class Euler {

public:
    Ev f;
    Euler() {}
    S next_step (S& state, double dt) {
        S new_state = state + f.pp(state) * dt;
        return new_state;
    }

    S n_step (S & state, double dt, size_t n) {
        S new_state = next_step(state, dt);
        for (int i = 1; i < 0; ++i) {
            new_state = next_state(new_state, dt);
        }
        return new_state;
    }
};

template <class S, class Ev>
class Heun {

public:
    Ev f;
    Heun () {} 
    S next_step (S& state, double dt) {
        S new_state1 = state + f.pp(state) * dt;
        S new_state2 = state + (f.pp(state) + f.pp(new_state1)) * 0.5*dt;
        return new_state2;
    }

    S n_step (S & state, double dt, size_t n) {
        S new_state = next_step(state, dt);
        for (int i = 1; i < 0; ++i) {
            new_state = next_state(new_state, dt);
        }
        return new_state;
    }
};

template <class S, class Ev>
class Rg45 {

public:
    Ev f;
    Rg45 () {} 
    S next_step (S state, double dt) {
        S new_state1 = f.pp(state);
        S new_state2 = f.pp(state + new_state1 * 0.5*dt);
        S new_state3 = f.pp(state + new_state2 * 0.5*dt);
        S new_state4 = f.pp(state + new_state1 * dt);
        S new_state = state + (new_state1 + new_state2*2 + new_state3*2 + new_state4) * (dt/6);
        return new_state;
    }

    S n_step (S & state, double dt, size_t n) {
        S new_state = next_step(state, dt);
        for (int i = 1; i < 0; ++i) {
            new_state = next_state(new_state, dt);
        }
        return new_state;
    }
};


template <class S, class M>
class Solver {
public:

std::vector<S> data;

Solver(State s0, double dt, double time) {
    data.push_back(s0);
    M f;
    for (double i = 0; i < time; i += dt) {
        s0 = f.next_step(s0, dt);
        data.push_back(s0);
    }
}

Solver(State s0, double dt, double time, size_t n) {
    data.push_back(s0);
    M f;
    for (double i = 0; i < time; i += dt) {
        s0 = f.n_step(s0, dt, n);
        data.push_back(s0);
    }
}

    void write (std::string file){

        std::ofstream out(file, std::ios::binary);
        assert(out.good());

        for (long i = 0; i < data.size(); ++i) {
            out.write((char*)&data[i], sizeof(S));
        }
        out.close();
        std::cout << "data was written to " << file << std::endl;
    }
};



int main() {
    
    double x = M_PI_4, v = 0, w = 1, dt = 0.01, time = 10;
    std::string type = "Euler";

    std::string s;
    std::ifstream in("config.txt");
    in >> s; x    = std::stod(s);
    in >> s; v    = std::stod(s);
    in >> s; w    = std::stod(s);
    in >> s; dt   = std::stod(s);
    in >> s; time = std::stod(s);
    in >> s; type = s;
    in.close();
    
    bool OK = 0;

    State s1(x, v, w);
    if (type == "Euler") {
        Solver<State, Euler<State, Evolution<State> > > sol1 (s1, dt, time);
        sol1.write("data.binary");
        OK = 1;
    }
    if (type == "Heun") {
        Solver<State, Heun<State, Evolution<State> > > sol1 (s1, dt, time);
        sol1.write("data.binary");
        OK = 1;
    }
    if (type == "Rg45") {
        Solver<State, Rg45<State, Evolution<State> > > sol1 (s1, dt, time);
        sol1.write("data.binary");
        OK = 1;
    }
    if(OK == 0) {
        std::cout << "I don`t know this method";
    }

    return 0;
}