#include <functional>
#include <iostream>
#include <list>
#include <optional>
#include <random>
#include <vector>
#include <string>
#include <math.h>

# define M_PI           3.14159265358979323846  /* pi */

/**
 * @brief random number generator
 *
 * @param std::random_device
 * @return std::mt19937
 */
std::mt19937 mt_generator((std::random_device()) ());
/**
 * @brief goal function domain.
 */
using domain_t = std::vector<double>;

std::ostream &operator<<(std::ostream &o, domain_t &d) {
    o << d[0] << " " << d[1];
    return o;
}

/**
 * @brief calculate minimum point using hill climbing algorithm
 *
 * @param f goal function
 * @param start_point the start point for calculations
 * @param get_close_points function generating neighbours
 * @param max_iterations number of iterations
 * @return domain_t the domain ponint where the function f has minimum
 */
domain_t hill_climbing(const std::function<double(domain_t)> &f, domain_t start_point, std::function<std::vector<domain_t>(domain_t)> get_close_points, int max_iterations) {
    domain_t best_p = start_point;
    for (int iteration = 0; iteration < max_iterations; iteration++) {
        auto close_points = get_close_points(best_p);
        auto best_neighbour = *std::min_element(close_points.begin(), close_points.end(), [f](auto a, auto b){return f(a) > f(b);});
        if (f(best_neighbour) < f(best_p)) best_p = best_neighbour;
    }
    return best_p;
}


double get_random_cud_variable(domain_t point){
    std::uniform_real_distribution dist(0.0,1.0);
    double random_u = dist(mt_generator);
    double f_u;
    if (random_u<=point[0]){
        f_u = 0;
    }else if(point[0]< random_u && random_u <= point[1]){
        f_u = (point[0] - random_u)/(point[1] - random_u);
    }else{
        f_u = 1;
    }
    return f_u;
};

double calc(double best, double neighbour, int iter) {
    double upper = -1 * abs(neighbour - best);
    double lower =  1/iter;
    double result = exp(upper / lower);
    return result;
}

domain_t simulated_annealing_method(
        const std::function<double(domain_t)> &f, domain_t start_point,
        std::function<std::vector<domain_t>(domain_t)> get_close_points,
        int max_iterations) {
    using namespace std;
    domain_t best_p = start_point;
    double uk = get_random_cud_variable(best_p);
    for (int iteration = 1; iteration < max_iterations; iteration++) {
        domain_t best_neighbour = std::vector<double>(uk);
        std::transform (best_p.begin(), best_p.end(), best_neighbour.begin(), best_p.begin(), std::plus<double>());
        if (f(best_neighbour) <= f(best_p)) {
            best_p = best_neighbour;
        } else {
            if (uk < calc(f(best_p), f(best_neighbour), iteration)) {
                best_p = best_neighbour;
            }
        }
    }
    return best_p;

}

/**
 * @brief full review method that will check every domain point
 *
 * @param f goal function
 * @param domain_generator the function that will generate consecutive points
 * from the domain and it will return empty when there are no more points to
 * check
 * @return domain_t the point where f has its minimum
 */
domain_t brute_force_method(
        const std::function<double(domain_t)> &f,
        const std::function<std::optional<domain_t>()> &domain_generator) {
    auto best_p = domain_generator();
    for (auto current_p = best_p; current_p.has_value();
         current_p = domain_generator()) {
        if (f(current_p.value()) < f(best_p.value())) {
            best_p = current_p;
        }
    }
    return best_p.value_or(domain_t());
}

double beale(double x, double y) {
    return pow((1.5 - x + x * y), 2) + pow((2.25 - x + x * pow(y, 2)), 2) + pow((2.625 - x + x * pow(y, 3)), 2);
};

double booth(double x, double y) {
    return pow(x + 2 * y - 7, 2) + pow(2 * x + y - 5, 2);
}

double matyas(double x, double y) {
    return 0.26 * (pow(x, 2) + pow(y, 2)) - 0.48 * x * y;
}

double himmelblau(double x, double y) {
    double part_one = pow(pow(x,2)+y-11,2);
    double part_two = pow(x+pow(y,2)-7,2);
    return part_one+part_two;
}

double cross_in_tray(double x, double y){
    double part_one = sin(x)*sin(y);
    double part_two = abs(100-((sqrt(pow(x,2)+ pow(y,2)))/M_PI));
    return -0.0001*pow((abs(part_one*exp(part_two))+1),0.1);
}

double func(std::string name, double x, double y) {
    if (name == "beale") {
        return beale(x, y);
    } else if (name == "booth") {
        return booth(x, y);
    } else if (name == "matyas") {
        return matyas(x, y);
    } else if (name == "himmelblau"){
        return himmelblau(x,y);
    } else if (name == "himmelblau"){
        return himmelblau(x,y);
    } else if (name == "cross"){
        return cross_in_tray(x,y);
    }
}

int main(int argc, char **argv) {
    using namespace std;
    cout << "Funkcja: " << argv[1] << "\nKrance dziedziny: " << argv[2] << ' ' << argv[3] << "\nLiczba iteracji: " <<argv[4] << endl;
    try {
        string name = argv[1];
        double min = atof(argv[2]);
        double max = atof(argv[3]);
        int iterations = atof(argv[4]);
        auto numbers_f_v = [&name](domain_t x) { return func(name, x[0], x[1]); };
        auto get_random_point = [&min, &max]() -> domain_t {
            std::uniform_real_distribution<double> distr(min, max);
            return {distr(mt_generator), distr(mt_generator)};
        };
        auto get_close_points_random = [&min, &max](domain_t p0) -> std::vector<domain_t> {
            std::uniform_real_distribution<double> distr(min, max);
            return {{distr(mt_generator), distr(mt_generator)}};
        };
        const double precision = 1.0 / 128;
        auto numbers_f_generator = [precision, &min, &max]() -> std::optional<domain_t> {
            static domain_t p = {min, min};
            int i = 0;
            for (i; i < p.size(); i++) {
                p[i] = p[i] + precision;
                if (p[i] < max) return std::optional(p);
                p[i] = min;
            }
            return {};
        };
        auto best0 =
                hill_climbing(numbers_f_v, get_random_point(), get_close_points_random, iterations);
        std::cout << "# hill climbing x = " << best0[0] << " " << best0[1] << std::endl;

        auto best1 =
                simulated_annealing_method(numbers_f_v, get_random_point(), get_close_points_random, iterations);
        std::cout << "# simulated annealing x = " << best1[0] << " " << best1[1] << std::endl;


        auto best3 = brute_force_method(numbers_f_v, numbers_f_generator);
        std::cout << "# brute x = " << best3[0] << " " << best3[1] << std::endl;

        return 0;
    } catch (std::out_of_range aor) {
        cout << "Uzupełnij";
    }

    //Beale działa
    //Booth działa
    //Matyas działa

    /*Beale Zakres [-4.5,4.5]
    double beale_x = 3;
    double beale_y = 0.5;
    cout << beale(beale_x, beale_y) << endl;

    Himmelblau Zakres [-5,5]
    double himmelblau_x = {3, -2.805118, -3.779310, 3.584428};
    double himmelblau_y = {2, 3.131312, -3.283186, -1.848126};

    Cross-in-tray [-10,10]
    double cross-in-tray_x = {1.34941, 1.34941, -1.34941, -1.34941};
    double cross-in-tray_y = {-1.34941, 1.34941, 1.34941,- 1.34941};

    Booth Zakres [-10,10]
    double booth_x = 1;
    double booth_y = 3;
    cout << booth(booth_x, booth_y) << endl;

    Matyas Zakres [-10, 10]
    double matyas_x = 0;
    double matyas_y = 0;
    cout << matyas(matyas_x, matyas_y) << endl;
     */

}

//Brudnopis
/*
 * auto sphere_f_v = [](domain_t x) { return x[0] * x[0] + x[1] * x[1]; };
    auto sphere_f_generator = [precision]() -> std::optional<domain_t> {
        static domain_t p = {-10, -10};
        int i = 0;
        for (i; i < p.size(); i++) {
            p[i] = p[i] + precision;
            if (p[i] < 10) return std::optional(p);
            p[i] = -10;
        }
        return {};
    };

 auto get_close_points_random = [](domain_t p0) -> std::vector<domain_t> {
        std::uniform_real_distribution<double> distr(-10, 10);
        return {{distr(mt_generator), distr(mt_generator)}};
    };
    const double precision = 1.0 / 16;

     auto best0 =
            tabu_method(rastrigin_f_v, get_random_point(), get_close_points, 1000);
    std::cout << "# tabu x = " << best0[0] << " " << best0[1] << std::endl;
    auto best1 = hill_climbing(sphere_f_v, get_random_point(),
     get_close_points_random, 1000000);
     std::cout << "# hill_climbing x = " << best1[0] << " " << best1[1] <<
     std::endl; auto best2 = brute_force_method(sphere_f_v,
     sphere_f_generator); std::cout << "# hill_climbing x = " << best2[0] << "
     " << best2[1] << std::endl;
    */