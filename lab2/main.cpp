#include <iostream>
#include <vector>
#include <functional>

double goal_function(std::vector<double> x) {
    return 0;
}

using domain_t = std::vector<double>;

std::random_device rd;
std::mt19937 mt_generator(rd());
domain_t hill_climbing(std::function<double(domain_t)>f, domain_t minimal_d, domain_t maximum_d){
    std::uniform_int_distribution<int> dist(0,9);
    for(int i = 0; i < minimal_d.size())
};

//std::vector<double> brute_force(std::function<double(std::vector<double)> f){
//
//}
/**
 * domain - generate domain points. Throws exception when all the points were returned
 */
auto brute_force = [](auto f, auto domain) {
    auto current_p = domain();
    auto best_point = current_p;
    try {
        while (true) {
            if (f(current_p) < f(best_point)) {
                best_point = current_p;
            }
            current_p = domain();
        }
    } catch (std::exception &e) {
    }
    return best_point;
};

int main() {
    auto sphere_f = [](double x) { return x * x; };
    double current_sphere_x = -10;
    auto sphere_generator = [&]() {
        current_sphere_x += 1.0/128.0;
        if (current_sphere_x >= 10) throw std::invalid_argument("finished");
        return current_sphere_x;
    };
    auto best_point = brute_force(sphere_f, sphere_generator);
    std::cout << "best x = " << best_point << std::endl;
    return 0;
}
