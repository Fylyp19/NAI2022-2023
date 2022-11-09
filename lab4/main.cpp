#include <iostream>
#include <vector>
#include <functional>
#include <random>
#include <iterator>
#include <sstream>
#include <bitset>

using namespace std;
std::random_device rd;
std::mt19937 mt_generator{rd()};
auto genetic_algorithm = [](
        auto initial_population, auto fitness, auto term_condition,
        auto selection, double p_crossover,
        auto crossover, double p_mutation, auto mutation) {

    uniform_real_distribution<double> uniform(0.0, 1.0);
    auto population = initial_population;
    vector<double> population_fit = fitness(population);
    while (!term_condition(population, population_fit)) {
        auto parents_indexes = selection(population_fit);
        decltype(population) new_population;
        for (int i = 0; i < parents_indexes.size(); i += 2) {
            decltype(initial_population) offspring = {population[i], population[i + 1]};
            if (uniform(mt_generator) < p_crossover) {
                offspring = crossover(offspring);
            }
            for (auto chromosome: offspring) new_population.push_back(chromosome);
        }
        for (auto &chromosome: new_population) {
            chromosome = mutation(chromosome, p_mutation);
        }
        population = new_population;
        population_fit = fitness(population);
    }
    //cout << geno_to_feno(population)  << endl;
    return population;
};
using chromosome_t = std::vector<int>;
using population_t = std::vector<chromosome_t>;

string v_to_str(vector<int> v){
    stringstream ss;
    copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, ""));
    string s = ss.str();
    return s.substr(0,s.length());
}

string bi_str_to_dec_str(string s){
    unsigned long long dec_numb = std::stoull(s,0,2);
    return to_string(dec_numb);
}

std::vector<double> geno_to_feno(std::vector<int> geno){
    using namespace std;

    vector<double> d_v;

    vector<int> f_w_v;
    vector<int> f_f_v;
    vector<int> s_w_v;
    vector<int> s_f_v;

    f_w_v.insert(f_w_v.end(), make_move_iterator(geno.begin()),make_move_iterator(geno.begin()+9));
    geno.erase(geno.begin(), geno.begin()+9);

    f_f_v.insert(f_f_v.end(), make_move_iterator(geno.begin()),make_move_iterator(geno.begin()+42));
    geno.erase(geno.begin(), geno.begin()+42);

    s_w_v.insert(s_w_v.end(), make_move_iterator(geno.begin()),make_move_iterator(geno.begin()+9));
    geno.erase(geno.begin(), geno.begin()+9);

    s_f_v.insert(s_f_v.end(), make_move_iterator(geno.begin()),make_move_iterator(geno.begin()+42));
    geno.erase(geno.begin(), geno.begin()+42);

    cout << v_to_str(f_w_v) << endl;
    cout << v_to_str(f_f_v)<< endl;
    cout << v_to_str(s_w_v) << endl;
    cout << v_to_str(s_f_v)<< endl;


    string f_w_int = bi_str_to_dec_str(v_to_str(f_w_v));
    cout << f_w_int << endl;
    string f_f_int = bi_str_to_dec_str(v_to_str(f_f_v));
    cout << f_f_int << endl;
    string s_w_int = bi_str_to_dec_str(v_to_str(s_w_v));
    cout << s_w_int << endl;
    string s_f_int = bi_str_to_dec_str(v_to_str(s_f_v));
    cout << s_f_int << endl;

    double f = stod(f_w_int+'.'+f_f_int);
    double s = stod(s_w_int+'.'+s_f_int);

    d_v.push_back(f);
    d_v.push_back(s);

    return d_v;
}

double eggholder(double x1, double x2) {

    double term1 = -(x2+47) * sin(sqrt(abs(x2+x1/2+47)));
    double term2 = x1 * sin(sqrt(abs(x1-(x2+47))));
    return term1 - term2;
}

std::vector<double> fitness_function(population_t pop) {
    std::vector<double> result;
    for(vector<int> part: pop) {
        vector<int> geno = part;
        vector<double> feno = geno_to_feno(geno);
        double egg = eggholder(feno[0], feno[1]);
        result.push_back((512 - abs(egg)) * 0.001);
        cout << result[0] << endl;
    }
    return result;
}

std::vector<int> selection_empty(std::vector<double> fitnesses) {
    return {};
}

std::vector<chromosome_t> crossover_empty(std::vector<chromosome_t> parents) {
    return parents;
}

chromosome_t mutation_empty(chromosome_t parents, double p_mutation) {
    return parents;
}




std::vector<int> binary_generator(int length){
    using namespace std;
    mt19937 mersenne_engine {rd()};
    uniform_int_distribution<int> dist {0, 1};
    auto gen = [&dist, &mersenne_engine](){
        return dist(mersenne_engine);
    };
    vector<int> vec(length);
    generate(begin(vec), end(vec),gen);
    return vec;
}

int main() {
    using namespace std;
    int length = 100+(23211%10)*2;
    cout << length << endl;

    vector<int> genotype = binary_generator(length);
    population_t population;
    population.push_back(genotype);

    for (int i: genotype) {
        std::cout << i;
    }

    cout << endl;
    cout << endl;

    //vector<double> feno = geno_to_feno(genotype);
    //for(double i : feno){
    //    cout << i << endl;
    //}

    //cout << eggholder(feno[0],feno[1]) << endl;

    auto result = genetic_algorithm(population,
                                    fitness_function,
                                    [](auto a, auto b) { return true; },
                                    selection_empty, 1.0,
                                    crossover_empty,
                                    0.01, mutation_empty);
    /*for (chromosome_t chromosome: result) {
        cout << "[";
        for (int p: chromosome) {
            cout << p;
        }
        cout << "] ";
    }*/

    return 0;
}


/*
if(placeholder_part[i] != 0){
mini.push_back(placeholder_part[i]);

} else{
placeholder.insert(placeholder.end(),mini);
mini.clear();
}
 */