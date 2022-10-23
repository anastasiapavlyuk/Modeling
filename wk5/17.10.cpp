#include <iostream>
#include <cmath>



float sum(long n) { // подсчет суммы гармонического ряда
    float y = 0;
    for (float i = 1; i< std::pow(10,n); ++i) { 
        y += 1/i; 
        }
    std::cout << y << std::endl;
    return y;
}

float fact(long x) {  // формула факториала
    float y = 1;
    for (float i = 1; i <= x; ++i) {
        y *= i;
    }
    return y;
}

float macloren (float x, long n) {  // формула маклорена для синуса
    float y = 0;
    for (long i = 1; i < std::pow(10, n); ++i) {
        y += std::pow(-1, i-1) * std::pow(x, 2*i-1) / fact(2*i - 1);
    }
    std::cout << "teor" << " " << std::sin(x) << " " << "calculation" << " " << y << std::endl;
    return y;
}

float integr (float min, max, float step) {
    float y = 0;
    for (float i = min; i < std::pow(10, max); i += step) {
        y += (std::exp(-std::pow(i+step, 2)) + std::exp(-std::pow(i, 2))) * step / 2;
    }
    std::cout << "Teor" << " " << std::sqrt(M_PI)/2 << "calculation" << " " << y << std::endl;
    return y;
}



int main() {
    unsigned n = 25;
    float x = 0;

    for (unsigned i = 1; i < n; ++i) {
        x += 1/(std::pow(2, i));
        std::cout << i << " " << 1/(std::pow(2, i)) << " " << x << std::endl;
    }

    // становится тождественной единицей после 21 шага

    for (int i = 1; i < 8; ++i) {
        std::cout << i << " ";
        sum(i);
    } 

    // сходится на 10^7 к 15,4037, дальше не хватет вычислительной мощности

    macloren(M_PI_2, 2); 

    // сходится за ~ 100 шагов

    integr(0, 1, 0.01);
    integr(0, 1, 0.1);
    integr(0, 1, 0.5);
    integr(0, 2, 1);
    integr(0, 3, 1);
    integr(0, 7, 1);

    // сходится за менее чем 10 шагов при точности не более, чем 0.5; при низкой точности (step = 1)
    // сходится до третьей значащей цифры вне зависимости от количества шагов

    return 0;
}