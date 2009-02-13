#include <iostream>
#include <fstream>
#include <acotreewidth/decomp.h>

EliminationGraph::EliminationGraph(const Graph &graph) {
  size__ = graph.number_of_vertices();
  a__ = new unsigned int*[size__];
  initial_degrees__ = new unsigned int[size__];
  degrees__ = new unsigned int[size__];
  t__ = new bool*[size__];
  eliminated__ = new bool[size__];
  memset(eliminated__,false,sizeof(bool)*size__);
  nr_eliminations__ = 0;

  for(unsigned int k=0;k<size__;k++) {
    a__[k] = new unsigned int[size__];
    t__[k] = new bool[size__];
    memset(t__[k],false,sizeof(bool)*size__);
  }


  for(unsigned int i=0;i<size__;i++) {
    std::vector<unsigned int> neighbours = graph.get_neighbours(i);
    for(unsigned int j=0;j<neighbours.size();j++) {
      a__[i][j] = neighbours[j];
      t__[i][neighbours[j]] = true;
      t__[neighbours[j]][i] = true;
    }
    initial_degrees__[i] = neighbours.size();
  }
  memcpy(degrees__, initial_degrees__, sizeof(unsigned int)*size__);
}

EliminationGraph::~EliminationGraph() {
  for(unsigned int i=0;i<size__;i++) {
    delete[] a__[i];
    delete[] t__[i];
  }
  delete[] a__;
  delete[] initial_degrees__;
  delete[] degrees__;
  delete[] t__;
  delete[] eliminated__;
}

void EliminationGraph::eliminate(unsigned int vertex) {
  for(unsigned int i=0;i<degrees__[vertex];i++) {
    int n1 = a__[vertex][i];
    if(!eliminated__[n1]) {
      for(unsigned int j=i+1;j<degrees__[vertex];j++) {
        int n2 = a__[vertex][j];
        if(!eliminated__[n2] && !t__[n1][n2]) {
          t__[n1][n2] = true;
          t__[n2][n1] = true;
          a__[n1][degrees__[n1]] = n2;
          a__[n2][degrees__[n2]] = n1;
          degrees__[n1]++;
          degrees__[n2]++;
        }
      }
    }
    t__[vertex][n1] = false;
    t__[n1][vertex] = false;
  }

  nr_eliminations__++;
  eliminated__[vertex] = true;
}

unsigned int EliminationGraph::get_degree(unsigned int vertex) const {
  unsigned int degree = 0;
  for(unsigned int i=0;i<degrees__[vertex];i++) {
    if(!eliminated__[a__[vertex][i]]) {
      degree++;
    }
  }
  return degree;
}

unsigned int EliminationGraph::min_fill(unsigned int vertex) const {
  unsigned int min_fill = 0;
  for(unsigned int i=0;i<degrees__[vertex];i++) {
    int n1 = a__[vertex][i];
    if(!eliminated__[n1]) {
      for(unsigned int j=i+1;j<degrees__[vertex];j++) {
        int n2 = a__[vertex][j];
        if(!eliminated__[n1] && !eliminated__[n2] && !t__[n1][n2]) {
          min_fill++;
        }
      }
    }
  }
  return min_fill;
}

std::vector<unsigned int> EliminationGraph::get_neighbours(unsigned int vertex) const {
  std::vector<unsigned int> neighbours;
  for(unsigned int i=0;i<degrees__[vertex];i++) {
    if(!eliminated__[a__[vertex][i]]) {
      neighbours.push_back(a__[vertex][i]);
    }
  }
  return neighbours;
}

unsigned int EliminationGraph::number_of_vertices() {
  return size__ - nr_eliminations__;
}

void EliminationGraph::rollback() {
  nr_eliminations__ = 0;
  memset(eliminated__,false,sizeof(bool)*size__);
  for(unsigned int k=0;k<size__;k++) {
    memset(t__[k],false,sizeof(bool)*size__);
  }

  for(unsigned int i=0;i<size__;i++) {
    for(unsigned int j=0;j<initial_degrees__[i];j++) {
      t__[i][a__[i][j]] = true;
      t__[a__[i][j]][i] = true;
    }
  }
  memcpy(degrees__, initial_degrees__, sizeof(unsigned int)*size__);
}

std::vector<unsigned int> get_max_clique_positions(EliminationGraph &elim_graph, const std::vector<unsigned int> &solution) {
  unsigned int max_clique = 0;
  std::vector<unsigned int> max_clique_positions;
  for(unsigned int i=0;i<solution.size();i++) {
    unsigned int w = elim_graph.get_degree(solution[i]);

    if(w == max_clique) {
      max_clique_positions.push_back(i);
    } else if(max_clique < w) {
      max_clique_positions.clear();
      max_clique_positions.push_back(i);
      max_clique = w;
    } 

    elim_graph.eliminate(solution[i]);
    if(solution.size()-i <= max_clique) {
      break;
    }
  }
  elim_graph.rollback();
  return max_clique_positions;
}

MaxCliqueRandomNeighbour::MaxCliqueRandomNeighbour(EliminationGraph &graph, std::vector<unsigned int> solution) {
  elim_graph_ = &graph;
  solution_ = solution;
  neighbour_ = solution;
  has_next_neighbour_ = true;
}

MaxCliqueRandomNeighbour::~MaxCliqueRandomNeighbour() {
}

void MaxCliqueRandomNeighbour::set_solution(std::vector<unsigned int> solution) {
  solution_ = solution;
  neighbour_ = solution;
  has_next_neighbour_ = true;
}

std::vector<unsigned int> MaxCliqueRandomNeighbour::get_solution() {
  return solution_;
}

bool MaxCliqueRandomNeighbour::has_next_neighbour_solution() {
  return has_next_neighbour_;
}

std::vector<unsigned int> MaxCliqueRandomNeighbour::next_neighbour_solution() {
  std::vector<unsigned int> max_cliques = get_max_clique_positions(*elim_graph_, solution_);
  unsigned int max_clique_pos = 0;
  unsigned int other_pos = 0;
  while(max_clique_pos == other_pos) {
    max_clique_pos = max_cliques[random_number(max_cliques.size())];
    other_pos = random_number(solution_.size());
  }

  unsigned int tmp = neighbour_[max_clique_pos];
  neighbour_[max_clique_pos] = neighbour_[other_pos];
  neighbour_[other_pos] = tmp;
  has_next_neighbour_ = false;
  return neighbour_;
}

MaxCliqueNeighbourhood::MaxCliqueNeighbourhood(EliminationGraph &graph, std::vector<unsigned int> solution) {
  elim_graph_ = &graph;
  solution_ = solution;
  std::vector<unsigned int> max_cliques = get_max_clique_positions(*elim_graph_, solution_);
  max_clique_pos_ = max_cliques[random_number(max_cliques.size())];
  swap_pos_ = 0;
}

MaxCliqueNeighbourhood::~MaxCliqueNeighbourhood() {
}

void MaxCliqueNeighbourhood::set_solution(std::vector<unsigned int> solution) {
  solution_ = solution;
  std::vector<unsigned int> max_cliques = get_max_clique_positions(*elim_graph_, solution_);
  max_clique_pos_ = max_cliques[random_number(max_cliques.size())];
  swap_pos_ = 0;
}

std::vector<unsigned int> MaxCliqueNeighbourhood::get_solution() {
  return solution_;
}

bool MaxCliqueNeighbourhood::has_next_neighbour_solution() {
  return swap_pos_ < solution_.size();
}

std::vector<unsigned int> MaxCliqueNeighbourhood::next_neighbour_solution() {
  std::vector<unsigned int> neighbour = solution_;
  unsigned int tmp = neighbour[max_clique_pos_];
  neighbour[max_clique_pos_] = neighbour[swap_pos_];
  neighbour[swap_pos_] = tmp;
  swap_pos_++;
  return neighbour;
}

DecompLocalSearch::DecompLocalSearch(std::vector<unsigned int> initial_solution, EvaluationFunction &eval_func, Neighbourhood &neighbourhood) : LocalSearch(initial_solution, eval_func, neighbourhood) {
}

void DecompLocalSearch::search_neighbourhood() {
  const std::vector<unsigned int> &solution = neighbourhood_->next_neighbour_solution();
  double quality = eval_func_->eval_solution(solution);
  if(quality > best_so_far_quality_) {
    best_so_far_solution_ = solution;
    best_so_far_quality_ = quality;
  }
  neighbourhood_->set_solution(best_so_far_solution_);
}

double Heuristic::min_degree(const EliminationGraph &graph, unsigned int vertex) {
  return 1.0 / (graph.get_degree(vertex) + 1);
}

double Heuristic::min_fill(const EliminationGraph &graph, unsigned int vertex) {
  return 1.0 / (graph.min_fill(vertex) + 1);
}

FileNotFoundException::FileNotFoundException(const char *filepath) {
  filepath_ = filepath;
}

const char *FileNotFoundException::what() const throw() {
  return filepath_;
}

HyperGraph &Parser::parse_hypertreelib(const char *filepath) throw(FileNotFoundException) {
  HyperGraph *graph = new HyperGraph();
  std::ifstream file(filepath);
  char buf[1024*8];

  if(!file) {
    throw FileNotFoundException(filepath);
  }

  std::map<std::string, unsigned int> vertex_ids;
  int edge_nr = 0;
  int vertex_nr = 0;
  while(file.good()) {
    std::string line;
    std::string atom;
    std::size_t found;
    file.getline(buf, 1024*8);
    line = std::string(buf);
    //strip whitespaces
    found = line.find(" ");
    while(found != std::string::npos) {
      line = line.replace(int(found), 1, "");
      found = line.find(" ");
    }

    if (line.length() > 0 && line[0] != '%' && line[0] != '\r' && line[0] != '<') {
      //strip ')' character
      found = line.find(")");
      line.replace(int(found), 1, "");
      //strip '.' character
      found = line.find(".");
      if(found != std::string::npos) {
        line.replace(int(found), 1, "");
      }

      found = line.find("(");
      atom = line.substr(0, int(found)); 
      line.replace(0, int(found)+1, "");
      found = line.find(",");
      std::vector<unsigned int> vars;
      while(found != std::string::npos) {
        std::string var = line.substr(0, int(found));
        if(vertex_ids.find(var) == vertex_ids.end()) {
          vars.push_back(vertex_nr);
          graph->set_vertex_label(vertex_nr, var);
          vertex_ids[var] = vertex_nr;
          vertex_nr++;
        } else {
          vars.push_back(vertex_ids[var]);
        }
        line.replace(0, int(found) + 1, "");
        found = line.find(",");
      }
      graph->add_hyperedge(edge_nr, vars);
      graph->set_edge_label(edge_nr, atom);
      edge_nr++;
    }
  }

  file.close();
  return *graph;
}
