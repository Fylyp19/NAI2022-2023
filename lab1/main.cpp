#include <any>
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <cmath>

using mojafunkcja_t = std::function<std::string(std::string)>;

void calculate(std::string x, mojafunkcja_t fun) {
    using namespace std;
    cout << sin(stod(x)) << endl;
}

int main(int argc, char **argv) {
    using namespace std;

    map<string, function<double(vector<double>)>> calculator;
    calculator["add"] = [](auto arg) {return arg.at(0) + arg.at(1);};
    calculator["mod"] = [](auto arg) { (int)arg.at(0) % (int)arg.at(1);};
    calculator["sin"] = [](auto x) { return x.at(0);};


    try {
        vector<string> task_args(argv, argv + argc);
        vector<string> numbers_args(argv+2, argv+argc);
        auto selected_f = task_args.at(1);
        auto x = argumenty.at(2);
        //calculate(x,calculator.at(selected_f));
        cout << calculate(x, selected_f) << endl;
    } catch (std::out_of_range aor) {
        cout << "Podaj argument. Dostepne to: ";
        for (auto [k, v] : calculator) cout << " " << k;
        cout << endl;
    }
    return 0;
}
