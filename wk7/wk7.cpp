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

class State_loss {
public:
    std::array<double, 2> state;
    double w;
    double g;

    State_loss(): w(0), g(0) {
        state[0] = 0;
        state[1] = 0;
    }
    
    State_loss(double x, double v, double w, double g): w(w), g(g) {
        state[0] = x;
        state[1] = v;
    }
    
    State_loss(const std::array<double, 2> s, double w, double g): state(s), w(w), g(g) {}

    double & operator [](int i) {
        return state[i];
    }

    State_loss operator + (State_loss other) {
        return State_loss ((*this)[0]+other[0], (*this)[1]+other[1], (*this).w, (*this).g);
    }

     State_loss operator * (double a) {
        return State_loss ((*this)[0]*a, (*this)[1]*a, (*this).w, (*this).g);
     }
};

class State_F {
public:
    std::array<double, 2> state;
    double w;
    double g;
    double time_s;
    double dt_s;
    double A;
    double w0;

    State_F(): w(0), g(0), time_s(0), dt_s(0.01), A(1), w0(1) {
        state[0] = 0;
        state[1] = 0;
    }

    State_F(double x, double v, double w, double g, double time, double dt, double a, double w0): w(w), g(g), time_s(time), dt_s(dt), A(a), w0(w0) {
        state[0] = x;
        state[1] = v;
    }

    State_F(double x, double v, double w, double g, double w0): w(w), g(g), time_s(0), dt_s(0.01), A(1), w0(w0) {
        state[0] = x;
        state[1] = v;
    }

    State_F(const std::array<double, 2> s, double w, double g, double time, double dt, double a, double w0): state(s), w(w), g(g), time_s(time), dt_s(dt), A(a), w0(w0) {}

    double & operator [](int i) {
        return state[i];
    }

    State_F operator + (State_F other) {
        return State_F ((*this)[0]+other[0], (*this)[1]+other[1], (*this).w, (*this).g, (*this).time_s, (*this).dt_s, (*this).A, (*this).w0);
    }

     State_F operator * (double a) {
        return State_F ((*this)[0]*a, (*this)[1]*a, (*this).w, (*this).g, (*this).time_s, (*this).dt_s, (*this).A, (*this).w0);
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

template <class S>
class Evolution_loss {
public:
    Evolution_loss() {}
    S pp(S s) {
        return S (s[1], - s.g*s[1] - s.w*s[0], s.w, s.g);
    }
};

template <class S>
class Evolution_F {
public:
    Evolution_F() {}

    S pp(S s) {
        return S (s[1], - s.g*s[1] - s.w*s[0] + s.A * cos(s.w0 * s.time_s * s.dt_s), s.w, s.g, s.time_s+=1, s.dt_s, s.A, s.w0);
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

Solver(S s0, double dt, double time) {
    data.push_back(s0);
    M f;
    for (double i = 0; i < time; i += dt) {
        s0 = f.next_step(s0, dt);
        data.push_back(s0);
    }
}

Solver(S s0, double dt, double time, size_t n) {
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
    
    double x = M_PI_4, v = 0, w = 1, g = 2, dt = 0.01, time = 10, dt_s = 0.01, time_s = 0, w0 = 6;
    std::string type = "Euler";

    std::string s1;
    std::ifstream in("config1.txt");
    in >> s1; x      = std::stod(s1);
    in >> s1; v      = std::stod(s1);
    in >> s1; w      = std::stod(s1);
    in >> s1; dt     = std::stod(s1);
    in >> s1; time   = std::stod(s1);
    in >> s1; type   = s1;
    in.close();
    
    bool OK = 0;

    State S1(x, v, w);
    if (type == "Euler") {
        Solver<State, Euler<State, Evolution<State> > > sol1 (S1, dt, time);
        sol1.write("data1.binary");
        OK = 1;
    }
    if (type == "Heun") {
        Solver<State, Heun<State, Evolution<State> > > sol1 (S1, dt, time);
        sol1.write("data1.binary");
        OK = 1;
    }
    if (type == "Rg45") {
        Solver<State, Rg45<State, Evolution<State> > > sol1 (S1, dt, time);
        sol1.write("data1.binary");
        OK = 1;
    }
    if(OK == 0) {
        std::cout << "I don`t know this method";
    }

    std::string s2;
    in.open("config2.txt");  
    in >> s2; x      = std::stod(s2);
    in >> s2; v      = std::stod(s2);
    in >> s2; w      = std::stod(s2);
    in >> s2; dt     = std::stod(s2);
    in >> s2; time   = std::stod(s2);
    in >> s2; type   = s2;
    in >> s2; g      = std::stod(s2);
    in.close();

    State_loss S2(x, v, w, g);
    if (type == "Euler") {
        Solver<State_loss, Euler<State_loss, Evolution_loss<State_loss> > > sol2 (S2, dt, time);
        sol2.write("data2.binary");
        OK = 1;
    }
    if (type == "Heun") {
        Solver<State_loss, Heun<State_loss, Evolution_loss<State_loss> > > sol2 (S2, dt, time);
        sol2.write("data2.binary");
        OK = 1;
    }
    if (type == "Rg45") {
        Solver<State_loss, Rg45<State_loss, Evolution_loss<State_loss> > > sol2 (S2, dt, time);
        sol2.write("data2.binary");
        OK = 1;
    }
    if(OK == 0) {
        std::cout << "I don`t know this method";
    }

    std::string s3;
    in.open("config3.txt");  
    in >> s3; x      = std::stod(s3);
    in >> s3; v      = std::stod(s3);
    in >> s3; w      = std::stod(s3);
    in >> s3; dt     = std::stod(s3);
    in >> s3; time   = std::stod(s3);
    in >> s3; type   = s3;
    in >> s3; g      = std::stod(s3);
    in >> s3; time_s = std::stod(s3);
    in >> s3; dt_s   = std::stod(s3);
    in >> s3; w0     = std::stod(s3);
    in.close();

    State_F S3(x, v, w, g, w0);
    if (type == "Euler") {
        Solver<State_F, Euler<State_F, Evolution_F<State_F> > > sol3 (S3, dt, time);
        sol3.write("data3.binary");
        OK = 1;
    }
    if (type == "Heun") {
        Solver<State_F, Heun<State_F, Evolution_F<State_F> > > sol3 (S3, dt, time);
        sol3.write("data3.binary");
        OK = 1;
    }
    if (type == "Rg45") {
        Solver<State_F, Rg45<State_F, Evolution_F<State_F> > > sol3 (S3, dt, time);
        sol3.write("data3.binary");
        OK = 1;
    }
    if(OK == 0) {
        std::cout << "I don`t know this method";
    }

    return 0;
}