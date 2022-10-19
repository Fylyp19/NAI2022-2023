#include <iostream>
#include <vector>
#include <functional>
#include <random>
#include <string>
#include <cmath>
#include <map>

std::mt19937 mt_generator((std::random_device())());
using domain_t = std::vector<double>;
using namespace std;

double beale(double x, double y){
    return pow((1.5 - x + x * y), 2) + pow((2.25 - x + x * pow(y, 2)), 2) + pow((2.625 - x + x * pow(y, 3)), 2);
};

double booth(double x, double y){
    return pow(x + 2*y - 7, 2) + pow(2*x + y - 5, 2);
}

double matyas(double x, double y){
    return 0.26*(pow(x,2)+pow(y,2))-0.48*x*y;
}

double func(string name, double x, double y){
    if(name == "beale"){
        return beale(x,y);
    } else if(name == "booth"){
        return booth(x,y);
    } else if(name == "matyas"){
        return matyas(x,y);
    }
}

domain_t hill_climbing(const std::function<double(domain_t)> &f, domain_t start_point, std::function<std::vector<domain_t>(domain_t)> get_close_points, int max_iterations) {
    domain_t best_p = start_point;
    for (int iteration = 0; iteration < max_iterations; iteration++) {
        auto close_points = get_close_points(best_p);
        auto best_neighbour = *std::min_element(close_points.begin(), close_points.end(), [f](auto a, auto b){return f(a) > f(b);});
        if (f(best_neighbour) < f(best_p)) best_p = best_neighbour;
    }
    return best_p;
}

int main(int argc, char **argv) {
    cout << "Funkcja: " << argv[1] << "\nKrance dziedziny: " << argv[2] << ' ' << argv[3] << endl;
    try{
        string name = argv[1];
        double min = atof(argv[2]);
        double max = atof(argv[3]);
        auto numbers_f_v = [&name](domain_t x) { return func(name, x[0], x[1]); };
        auto get_random_point = [&min, &max]() -> domain_t {
            std::uniform_real_distribution<double> distr(min, max);
            return {distr(mt_generator), distr(mt_generator)};
        };
        auto get_close_points_random = [&min, &max](domain_t p0) -> std::vector<domain_t> {
            std::uniform_real_distribution<double> distr(min, max);
            return {{distr(mt_generator), distr(mt_generator)}};
        };
        auto best1 = hill_climbing(numbers_f_v, get_random_point(), get_close_points_random, 1000000);
        std::cout << "Liczby minimalne = " << best1[0] << " " << best1[1] << std::endl;

    } catch (std::out_of_range aor){
        cout << "Uzupełnij";
    }
    //Beale działa
    //Booth działa
    //Matyas działa

    /*Beale Zakres [-4.5,4.5]
    double beale_x = 3;
    double beale_y = 0.5;
    cout << beale(beale_x, beale_y) << endl;

    Booth Zakres [-10,10]
    double booth_x = 1;
    double booth_y = 3;
    cout << booth(booth_x, booth_y) << endl;

    Matyas Zakres [-10, 10]
    double matyas_x = 0;
    double matyas_y = 0;
    cout << matyas(matyas_x, matyas_y) << endl;
    */
    return 0;
}
