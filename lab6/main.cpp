// zaczynamy od tego
#include <algorithm>
#include <functional>
#include <iostream>
#include <random>
#include <vector>
#include <cassert>
#include <sstream>
#include <string>
#include <iterator>

std::random_device rd;
std::mt19937 mt_generator(rd());
using switcher_t = bool;
using switcher_iter_v_t = int;
using chromosome_t = std::vector<int>;
using population_t = std::vector<chromosome_t>;
using fitness_f = std::function<double(const chromosome_t &)>;
using term_condition_f = std::function<bool(const population_t &, const std::vector<double> &)>;
using selection_f = std::function<int(const std::vector<double> &)>;
using crossover_f = std::function<std::vector<chromosome_t>(const std::vector<chromosome_t> &, double)>;
using mutation_f = std::function<chromosome_t(const chromosome_t, double)>;

void datas(std::vector<double> fitness){
    using namespace std;
    cout << "Srednia: " << accumulate(fitness.begin(),fitness.end(),0.0)/ fitness.size() << endl;
    cout << "Max: " << *max_element(fitness.begin(),fitness.end()) << endl;
    cout << "Min: " << *min_element(fitness.begin(),fitness.end()) << endl;

}

population_t genetic_algorithm(population_t initial_population,
                               switcher_t switcher,
                               switcher_iter_v_t switcher_iter_v,
                               fitness_f fitness,
                               term_condition_f term_condition,
                               selection_f selection, double p_crossover,
                               crossover_f crossover, double p_mutation,
                               mutation_f mutation) {
    using namespace std;
    uniform_real_distribution<double> uniform(0.0, 1.0);
    auto population = initial_population;
    vector<double> population_fit(population.size());
    transform(population.begin(), population.end(), population_fit.begin(),fitness);
    int switcher_iter = 0;
    while (!term_condition(population, population_fit)) {
        vector<int> parents_indexes(population.size());
        population_t new_population(population.size());
        // calculate fitness
        transform(population_fit.begin(), population_fit.end(),
                  parents_indexes.begin(),
                  [&](auto e) { return selection(population_fit); });
        if(switcher == true) {
            if(switcher_iter % switcher_iter_v == 0) {
                datas(population_fit);
            }
        }
        switcher_iter++;
        // perform crossover operations
        for (int i = 0; i < parents_indexes.size() - 1; i += 2) {
            vector<chromosome_t> offspring = {population[parents_indexes[i]], population[parents_indexes[i + 1]]};
            if (uniform(mt_generator) < p_crossover) {
                offspring = crossover(offspring, p_crossover);
            }
            new_population[i] = offspring[0];
            new_population[i + 1] = offspring[1];
        }
        for (auto &chromosome : new_population) {
            chromosome = mutation(chromosome, p_mutation);
        }
        population = new_population;
        std::transform(population.begin(), population.end(), population_fit.begin(),
                       fitness);
    }
    if(switcher == true) {
        datas(population_fit);
    }
    return population;
};

int selection_empty(std::vector<double> fitnesses) {
    return 0;
}

//void

int selection_tournament_2(std::vector<double> fitnesses) {
    std::uniform_int_distribution<int> uniform(0, fitnesses.size()-1);
    int a = uniform(mt_generator);
    int b = uniform(mt_generator);
    return (fitnesses[a]>fitnesses[b])?a:b;
}

int selection_roulette(std::vector<double> fitnesses) {
    using namespace std;
    double sum_of_f;
    for(int i = 0; i < fitnesses.size(); i++){
        sum_of_f += fitnesses[i];
    }
    //cout << sum_of_f << endl;
    std::vector<int> p_chromosom;
    for(int i = 0; i < fitnesses.size(); i++){
        p_chromosom.push_back((fitnesses[i])/sum_of_f);
    }
    std::uniform_int_distribution<int> uniform(0, sum_of_f);
    double a = uniform(mt_generator);
    int offset_down, offset_up;
    int pick = 0;
    for(int i = 0; i < fitnesses.size(); i++){
        offset_up +=  p_chromosom[i];
        if(offset_down < a && a < offset_up){
            pick = i;
            break;
        }
        offset_down += p_chromosom[i];
    }
    return fitnesses[pick];
}

std::vector<chromosome_t> crossover_empty(std::vector<chromosome_t> parents) {
    return parents;
}
std::vector<chromosome_t> crossover_two_point(std::vector<chromosome_t> parents) {
    using namespace std;
    uniform_int_distribution<int> locus(0,parents.at(0).size()-1);
    int a = locus(mt_generator);
    int b = locus(mt_generator);
    if (a > b) swap(a,b);
    auto children = parents;
    for (int i = a; i < b; i++) {
        swap(children[0][i],children[1][i]);
    }
    return children;
}

std::vector<chromosome_t> crossover_one_point(std::vector<chromosome_t> parents, double pc) {
    using namespace std;
    uniform_int_distribution<int> locus(0,parents.at(0).size()-1);
    int a = locus(mt_generator);
    auto children = parents;
    for (int i = a; i < parents.at(0).size(); i++) {
        swap(children[0][i],children[1][i]);
    }
    return children;
}
chromosome_t mutation_empty(const chromosome_t parent, double p_mutation) {
    return parent;
}
chromosome_t mutation_one_point(const chromosome_t parent, double p_mutation) {
    using namespace std;
    uniform_real_distribution<double> uni(0.0,1.0);
    // mutation??
    if (uni(mt_generator) < p_mutation) {
        uniform_int_distribution<int> locus(0,parent.size()-1);
        chromosome_t child = parent;
        auto l = locus(mt_generator);
        child[l] = 1 - child[l];
        return child;
    } else
        return parent;
}

chromosome_t mutation_many_points(const chromosome_t parent, double p_mutation) {
    using namespace std;
    uniform_real_distribution<double> n(0.0, parent.size());
    uniform_real_distribution<double> uni(0.0, 1.0);
    if (uni(mt_generator) < p_mutation) {
        uniform_int_distribution<int> locus(0, parent.size() - 1);
        chromosome_t child = parent;
        int a = n(mt_generator);
        for (int i = 0; i < a; i++) {
            auto l = locus(mt_generator);
            child[l] = 1 - child[l];
        }
        return child;
    } else{
        return parent;
    }
}

chromosome_t probabilistic_mutation(const chromosome_t parent, double p_mutation) {
    using namespace std;
    uniform_real_distribution<double> uni(0.0, 1.0);
    chromosome_t child = parent;
    for(int i = 0; i < parent.size(); i++){
        if (uni(mt_generator) < p_mutation) {
            if(parent[i] == 0){
                child[i] = 1;
            } else if(parent[i] == 1){
                child[i] = 0;
            }
        }
    }
    return child;
}
/************************************************************************
 *
 * MINIMAL SET OF FUNCTIONS DEFINING PROBLEM
 *
 *************************************************************************/
/**
 * @brief calculates fitness for every population element
 *
 * @param pop the population
 * @return std::vector<double> fitness function values
 */

std::ostream &operator<<(std::ostream &o, const chromosome_t chromosome) {
    for (const int p : chromosome) {
        o << p;
    }
    return o;
}
std::ostream &operator<<(std::ostream &o,
                         std::pair<population_t, fitness_f> pop) {
    for (const auto p : pop.first) {
        o << "{" << p << " " << (pop.second(p)) << "} ";
    }
    return o;
}

std::string v_to_str(std::vector<int> v) {
    using namespace std;
    stringstream ss;
    copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, ""));
    string s = ss.str();
    return s.substr(0, s.length());
}

std::string bi_str_to_dec_str(std::string s){
    using namespace std;
    unsigned long long dec_numb = stoull(s,0,2);
    return to_string(dec_numb);
}


std::vector<double> geno_to_feno(const chromosome_t &chromosome){
    using namespace std;
    vector<int> geno = chromosome;

    //for (int i = 0; i < geno.size(); i++){
    //    cout << geno[i] << endl;
    //}

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

    //cout << v_to_str(f_w_v) << endl;
    //cout << v_to_str(f_f_v)<< endl;
    //cout << v_to_str(s_w_v) << endl;
    //cout << v_to_str(s_f_v)<< endl;


    string f_w_int = bi_str_to_dec_str(v_to_str(f_w_v));
    //cout << f_w_int << endl;
    string f_f_int = bi_str_to_dec_str(v_to_str(f_f_v));
    //cout << f_f_int << endl;
    string s_w_int = bi_str_to_dec_str(v_to_str(s_w_v));
    //cout << s_w_int << endl;
    string s_f_int = bi_str_to_dec_str(v_to_str(s_f_v));
    //cout << s_f_int << endl;

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

double fitness_function(const chromosome_t &chromosome) {
    return std::accumulate(chromosome.begin(), chromosome.end(), 0);
}

double my_fitness(const chromosome_t &chromosome) {
    using namespace std;
    double result;
    auto geno = chromosome;
    vector<double> feno = geno_to_feno(geno);
    double egg = eggholder(feno[0], feno[1]);
    cout << egg << endl;
    result = abs((512 - abs(egg))) * 0.01;
    return result;
}

std::vector<chromosome_t> generate_initial_population(int n) {
    std::vector<chromosome_t> ret(n);
    std::uniform_int_distribution<int> uniform(0, 1);
    std::transform(ret.begin(), ret.end(), ret.begin(), [&](auto e) {
        chromosome_t c(102);
        for (int i = 0; i < c.size(); i++) c[i] = uniform(mt_generator);
        return c;
    });
    return ret;
}

int main(int argc, char **argv) {

    using namespace std;

    int size_of_population = atoi(argv[1]);
    int iterations = atoi(argv[2]);
    double p_crossover = atof(argv[3]);
    double p_mutation = atof(argv[4]);

    bool switch_p = false;
    string switch_point = argv[5];
    int switch_point_v = atoi(argv[6]);

    if(switch_point == "true"){
        switch_p = true;
    }

    population_t population = generate_initial_population(size_of_population);
    auto result = genetic_algorithm(
            population, switch_p, switch_point_v, my_fitness,
            [&iterations](auto a, auto b) {
                static int i = 0;
                i++;
                //cout << i << ": " << make_pair(a, fitness_function) << endl;
                return i >= iterations;
            },
            selection_tournament_2, p_crossover, crossover_one_point, p_mutation, probabilistic_mutation);
    cout << make_pair(result, fitness_function);
    cout << endl;
    return 0;
}