#include <string>
#include <vector>
#include <list>
#include <map>
#include <cmath>
#include "graph.h"

class PheromoneMatrix : protected Matrix<double> {
  protected:
    double evaporation_rate_;
    double initial_pheromone_;
  public:
    PheromoneMatrix(int vertices, double evaporation_rate, double initial_pheromone);
    double get(unsigned int v, unsigned int w);
    virtual void add(unsigned int v, unsigned int w, double amount);
    virtual void evaporate(unsigned int v, unsigned int w);
    void evaporate_all();
    double get_evaporation_rate();
    unsigned int size();
};

class MaxMinPheromoneMatrix : public PheromoneMatrix {
  private:
    double max_;
    double min_;
  public:
    MaxMinPheromoneMatrix(int vertices, double evaporation_rate, double initial_pheromone);
    void set_min(double min);
    void set_max(double max);
    void add(unsigned int v, unsigned int w, double amount);
    void evaporate(unsigned int v, unsigned int w);
};

class ACSPheromoneMatrix : public PheromoneMatrix {
  private:
    double epsilon_;
  public:
    ACSPheromoneMatrix(int vertices, double evaporation_rate, double initial_pheromone);
    void set_epsilon(double epsilon);
    void local_pheromone_update(unsigned int v, unsigned int w);
};

class Tour {
  private:
    std::vector<unsigned int> *vertices_;
    unsigned int capacity_;
    double length_;
  public:
    Tour(unsigned int vertices);
    ~Tour();
    double get_length();
    void set_length(double length);
    const std::vector<unsigned int> &get_vertices();
    void set_vertices(const std::vector<unsigned int> &vertices);
    Tour &add_vertex(unsigned int vertex);
    unsigned int &operator[](const unsigned int vertex) const;
    unsigned int size() const;
    unsigned int capacity();
    void clear();
    Tour &operator=(const Tour &t);
    bool operator<(const Tour &t);
};

class OptimizationProblem {
  public:
    virtual unsigned int get_max_tour_size() = 0;
    virtual unsigned int number_of_vertices() = 0;
    virtual std::map<unsigned int,double> get_feasible_start_vertices() = 0;
    virtual std::map<unsigned int,double> get_feasible_neighbours(unsigned int vertex) = 0;
    virtual double eval_tour(const std::vector<unsigned int> &tour) = 0;
    virtual double pheromone_update(unsigned int v, double tour_length) = 0;
    virtual void added_vertex_to_tour(unsigned int vertex) = 0;
    virtual bool is_tour_complete(const std::vector<unsigned int> &tour) = 0;
    virtual std::vector<unsigned int> apply_local_search(const std::vector<unsigned int> &tour) { return tour; }
    virtual void cleanup() = 0;
};

class Ant {
  protected:
    struct MultiMapComp {
      bool operator()(const double &p1, const double &p2) const;
    };
    Tour *tour;
    void update_tour_length(OptimizationProblem &op);
    void add_vertex_to_tour(OptimizationProblem &op, unsigned int vertex);
    std::multimap<double,unsigned int,MultiMapComp> get_feasible_vertices(OptimizationProblem &op, PheromoneMatrix &pheromones, double alpha, double beta);
    unsigned int current_vertex();
    unsigned int choose_next_vertex_with_likelihood(std::multimap<double,unsigned int,MultiMapComp> probabilities);
  public:
    Ant(unsigned int vertices);
    Ant(const Ant &ant);
    Ant &operator=(const Ant &ant);
    bool operator<(const Ant &ant);
    ~Ant();
    double get_tour_length();
    std::vector<unsigned int> get_vertices();
    void reset();
    void apply_local_search(OptimizationProblem &op);
    virtual void construct_solution(OptimizationProblem &op, PheromoneMatrix &pheromones, double alpha, double beta) {}
    void construct_rational_solution(OptimizationProblem &op, PheromoneMatrix &pheromones, double alpha, double beta);
    void construct_random_proportional_solution(OptimizationProblem &op, PheromoneMatrix &pheromones, double alpha, double beta);
    virtual void offline_pheromone_update(OptimizationProblem &op, PheromoneMatrix &pheromones, double weight=1.0) {}
};

class SimpleAnt : public Ant {
  public:
    SimpleAnt(unsigned int vertices);
    SimpleAnt(const SimpleAnt &ant);
    SimpleAnt &operator=(const SimpleAnt &ant);
    void construct_solution(OptimizationProblem &op, PheromoneMatrix &pheromones, double alpha, double beta);
    void offline_pheromone_update(OptimizationProblem &op, PheromoneMatrix &pheromones, double weight=1.0);
};

class ACSAnt : public Ant {
  private:
    double q0_;
  public:
    ACSAnt(unsigned int vertices);
    ACSAnt(const ACSAnt &ant);
    ACSAnt &operator=(const ACSAnt &ant);
    void construct_solution(OptimizationProblem &op, PheromoneMatrix &pheromones, double alpha, double beta);
    void construct_pseudorandom_proportional_solution(OptimizationProblem &op, PheromoneMatrix &pheromones, double alpha, double beta);
    void offline_pheromone_update(OptimizationProblem &op, PheromoneMatrix &pheromones, double weight=1.0);
    void set_q0(double q0);
};

enum LocalSearchType { NONE, ITERATION_BEST, ALL };

class AntColonyConfiguration {
  public:
    unsigned int number_of_ants;
    double alpha;
    double beta;
    bool stagnation_measure;
    double evaporation_rate;
    double initial_pheromone;
    LocalSearchType local_search;

    AntColonyConfiguration();
};

class ElitistAntColonyConfiguration : public AntColonyConfiguration {
  public:
    double elitist_weight;
    ElitistAntColonyConfiguration();
};

class RankBasedAntColonyConfiguration : public AntColonyConfiguration {
  public:
    unsigned int elitist_ants;
    RankBasedAntColonyConfiguration();
};

class MaxMinAntColonyConfiguration : public AntColonyConfiguration {
  public:
    unsigned int best_so_far_frequency;
    double a;
    MaxMinAntColonyConfiguration();
};

class ACSAntColonyConfiguration : public AntColonyConfiguration {
  public:
    double q0;
    double epsilon;
    ACSAntColonyConfiguration();
};

template<class T=Ant, class P=PheromoneMatrix> class AntColony {
  private:
    void construct_ants_solutions() {
      for(typename std::list<T>::iterator it=ants_->begin();it!=ants_->end();it++) {
        T &ant = (*it);
        ant.construct_solution(*problem_, *pheromones_, alpha_, beta_);
      }
    }

    void apply_local_search() {
      ants_->sort();
      if(local_search_ == ITERATION_BEST) {
        typename std::list<T>::iterator it_best = ants_->begin();
        (*it_best).apply_local_search(*problem_);
      } else if(local_search_ == ALL) {
        for(typename std::list<T>::iterator it=ants_->begin();it!=ants_->end();it++) {
          T &ant = (*it);
          ant.apply_local_search(*problem_);
        }
      }
    }

    void reset_ants() {
      for(typename std::list<T>::iterator it=ants_->begin();it!=ants_->end();it++) {
        it->reset();
      }
    }

    void update_best_tours_no_ls() {
      ants_->sort();
      (*best_iteration_no_ls_) = ants_->front();
      if(best_so_far_no_ls_->get_tour_length() > best_iteration_no_ls_->get_tour_length()) {
        (*best_so_far_no_ls_) = (*best_iteration_no_ls_);
      }
    }

    void update_best_tours() {
      ants_->sort();
      (*best_iteration_) = ants_->front();
      if(best_so_far_->get_tour_length() > best_iteration_->get_tour_length()) {
        (*best_so_far_) = (*best_iteration_);
      }
    }

    void update_stagnation_measure() {
      double average = 0.0;
      double standard_deviation = 0.0;
      for(typename std::list<T>::iterator it=ants_->begin();it!=ants_->end();it++) {
        average += it->get_tour_length();
      }
      average = average / ants_->size();

      for(typename std::list<T>::iterator it=ants_->begin();it!=ants_->end();it++) {
        standard_deviation += pow(it->get_tour_length() - average, 2);
      }
      standard_deviation *= 1.0 / ants_->size();
      standard_deviation = sqrt(standard_deviation);
      stagnation_measure_ = standard_deviation / average;
    }

    virtual void update_pheromones() = 0;
  protected:
    P *pheromones_;
    double alpha_;
    double beta_;
    bool compute_stagnation_measure_;
    LocalSearchType local_search_;
    double stagnation_measure_;
    std::list<T> *ants_;
    OptimizationProblem *problem_;
    T *best_so_far_;
    T *best_iteration_;
    T *best_so_far_no_ls_;
    T *best_iteration_no_ls_;

  public:
    AntColony(OptimizationProblem *problem, const AntColonyConfiguration &config) {
      problem_ = problem;
      ants_ = new std::list<T>(config.number_of_ants, T(problem->get_max_tour_size()));
      pheromones_ = new P(problem->number_of_vertices()+1, config.evaporation_rate, config.initial_pheromone);
      alpha_ = config.alpha;
      beta_ = config.beta;
      compute_stagnation_measure_ = config.stagnation_measure;
      local_search_ = config.local_search;
      stagnation_measure_ = 0.0;
      best_so_far_ = new T(problem->get_max_tour_size());
      best_iteration_ = new T(problem->get_max_tour_size());
      best_so_far_no_ls_ = new T(problem->get_max_tour_size());
      best_iteration_no_ls_ = new T(problem->get_max_tour_size());
    }

    ~AntColony() {
      delete problem_;
      delete ants_;
      delete pheromones_;
      delete best_so_far_;
      delete best_iteration_;
      delete best_so_far_no_ls_;
      delete best_iteration_no_ls_;
    };

    void run() {
      construct_ants_solutions();
      update_best_tours_no_ls();
      apply_local_search();
      update_best_tours();
      if(compute_stagnation_measure_) {
        update_stagnation_measure();
      }
      update_pheromones();
      reset_ants();
    }

    std::vector<unsigned int> get_best_tour() {
      return best_so_far_->get_vertices();
    }

    std::vector<unsigned int> get_best_tour_in_iteration() {
      return best_iteration_->get_vertices();
    }

    double get_best_tour_length() {
      return best_so_far_->get_tour_length();
    }

    double get_best_tour_length_in_iteration() {
      return best_iteration_->get_tour_length();
    }

    std::vector<unsigned int> get_best_tour_no_ls() {
      return best_so_far_no_ls_->get_vertices();
    }

    std::vector<unsigned int> get_best_tour_in_iteration_no_ls() {
      return best_iteration_no_ls_->get_vertices();
    }

    double get_best_tour_length_no_ls() {
      return best_so_far_no_ls_->get_tour_length();
    }

    double get_best_tour_length_in_iteration_no_ls() {
      return best_iteration_no_ls_->get_tour_length();
    }

    double get_stagnation_measure() {
      return stagnation_measure_;
    }
};

class SimpleAntColony : public AntColony<SimpleAnt> {
  public:
    SimpleAntColony(OptimizationProblem *problem, const AntColonyConfiguration &config);
  protected:
    void update_pheromones();
};

class ElitistAntColony : public AntColony<SimpleAnt> {
  private:
    double elitist_weight_;
  public:
    ElitistAntColony(OptimizationProblem *problem, const ElitistAntColonyConfiguration &config);
  protected:
    void update_pheromones();
};

class RankBasedAntColony : public AntColony<SimpleAnt> {
  private:
    unsigned int elitist_ants_;
  public:
    RankBasedAntColony(OptimizationProblem *problem, const RankBasedAntColonyConfiguration &config);
  protected:
    void update_pheromones();
};

class MaxMinAntColony : public AntColony<SimpleAnt, MaxMinPheromoneMatrix> {
  private:
    unsigned int best_so_far_frequency_;
    double a_;
  public:
    MaxMinAntColony(OptimizationProblem *problem, const MaxMinAntColonyConfiguration &config);
  protected:
    void update_pheromones();
};

class ACSAntColony : public AntColony<ACSAnt, ACSPheromoneMatrix> {
  public:
    ACSAntColony(OptimizationProblem *problem, const ACSAntColonyConfiguration &config);
  protected:
    void update_pheromones();
};
