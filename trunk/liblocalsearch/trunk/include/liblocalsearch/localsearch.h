#include <vector>
#include <cstdlib>
#include <ctime>

unsigned int random_number(unsigned int range=RAND_MAX);

class Neighbourhood {
  public:
    virtual ~Neighbourhood() {}
    virtual void set_solution(std::vector<unsigned int> solution) = 0;
    virtual std::vector<unsigned int> get_solution() = 0;
    virtual bool has_next_neighbour_solution() = 0;
    virtual const std::vector<unsigned int> &next_neighbour_solution() = 0;
};

class TwoOptNeighbourhood : public Neighbourhood {
  private:
    std::vector<unsigned int> solution_;
    std::vector<unsigned int> neighbour_;
    unsigned int prev_i_, prev_j_, i_,j_;
  public:
    void set_solution(std::vector<unsigned int> solution);
    std::vector<unsigned int> get_solution();
    bool has_next_neighbour_solution();
    const std::vector<unsigned int> &next_neighbour_solution();
    void swap(std::vector<unsigned int> &solution, unsigned int i, unsigned int j);
};

class EvaluationFunction {
  public:
    virtual ~EvaluationFunction() {}
    virtual double eval_solution(const std::vector<unsigned int> &solution) = 0;
};

class PerturbationFunction {
  public:
    virtual ~PerturbationFunction() {}
    virtual std::vector<unsigned int> perturbate(const std::vector<unsigned int> &solution) = 0;
};

class LocalSearch {
  private:
    virtual void search_neighbourhood() = 0;
  protected:
    std::vector<unsigned int> best_so_far_solution_;
    double best_so_far_quality_;
    EvaluationFunction *eval_func_;
    Neighbourhood *neighbourhood_;
  public:
    LocalSearch(std::vector<unsigned int> initial_solution, EvaluationFunction &eval_func, Neighbourhood &neighbourhood);
    virtual ~LocalSearch() {}
    std::vector<unsigned int> get_best_so_far_solution();
    double get_best_so_far_quality();
    void set_solution(std::vector<unsigned int> solution);
    double eval_solution(const std::vector<unsigned int> &solution);
    void search_iterations(int iterations);
    void search_iterations_without_improve(int iterations_without_improve);
};

class HillClimbing : public LocalSearch {
  private:
    void search_neighbourhood();
  public:
    HillClimbing(std::vector<unsigned int> initial_solution, EvaluationFunction &eval_func, Neighbourhood &neighbourhood);
};

class IterativeLocalSearch {
  private:
    LocalSearch *local_search_;
    PerturbationFunction *perturbation_func_;
  public:
    IterativeLocalSearch(LocalSearch *local_search, PerturbationFunction *perturbation_func);
    void run(int iterations=5);
};
