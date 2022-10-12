#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <cmath>


int main(int argc, char **argv) {
    using namespace std;
    map<string, function<double(vector<double>)> > func;
    func["sin"] = [](auto arg) { return sin(arg.at(0)); };
    func["add"] = [](auto arg) { return arg.at(0) + arg.at(1); };
    func["mod"] = [](auto arg) { return (int) arg.at(0) % (int) arg.at(1); };
    try {
        vector<string> func_args(argv, argv + 2);
        vector<string> val_args(argv + 2, argv + argc);
        auto selected_f = func_args.at(1);
        std::vector<double> doubleVector(val_args.size());
        std::transform(val_args.begin(), val_args.end(), doubleVector.begin(), [](const std::string &val) {
            return std::stod(val);
        });
        auto result = func.at(selected_f);
        cout << "Rozwiazanie = " << result(doubleVector);
    } catch (std::out_of_range aor) {
        cout << "Uzupełnij polecenie.";
    }

    return 0;
}