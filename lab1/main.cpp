#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

int main(int argc, char **argv) {
    using namespace std;
    map<string, function<double(vector<double>)> > func;
    func["add"] = [](auto arg) { return arg.at(0) + arg.at(1); };
    func["sin"] = [](auto arg) { return sin(arg.at(0)); };
    func["mod"] = [](auto arg) { return (int) arg.at(0) % (int) arg.at(1); };
    try {
        vector<string> func_args(argv, argv + 2);
        vector<string> num_args(argv + 2, argv + argc);
        auto selected_f = func_args.at(1);
        std::vector<double> doubleVector(num_args.size());
        std::transform(num_args.begin(), num_args.end(), doubleVector.begin(), [](const std::string &val) {
            return std::stod(val);
        });
        auto result = func.at(selected_f);
        cout << "Rozwiązanie = " << result(doubleVector);
    } catch (std::out_of_range aor) {
        cout << "Podaj brakujące argumenty.";
    }

    return 0;
}