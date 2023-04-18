#ifndef QAPT600_H
#define QAPT600_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdlib.h>     /* malloc, calloc, realloc, free */
#include <math.h>       /* round, floor, ceil, trunc */
#include <time.h>

using namespace std;

class QAPT
{
public:
   QAPT(string instance);
   ~QAPT();

   int run_ts(int nb_iterations);
   bool isVerbose;
   int n,max_time,lower_taboo_list_size, higher_taboo_list_size;
   int nrep,nr_iteration_before_aspiration,nr_iter_resizeTL;

protected:
private:
   const int n_max    = 240;       // maximal size of the problem  
   const int infinite = 999999999;

   vector<vector<int>> a,b;
   vector<vector<int>> taboo_list;   
   vector<int> solution,best_solution; 
   clock_t cpu1,cpu2;
   int cost, cpu, init_sol_seed, optimum;    
   int best_cost,taboo_list_seed,current_taboo_list_size;
   int aspirating_iteration, current_iteration;

   void read_problem(string instance, bool isVerbose);
   int  unif(int &seed, int low, int high);
   void swap(int &a, int &b);
   void delta_full_computation(int i, int j, vector<vector<int>>& delta, vector<int> &p);
   void delta_short_computation(int r,int s,int i,int j, vector<vector<int>> &delta, vector<int> &p);
   void initialize(vector<vector<int>> &delta, vector<int> &p);
   void find_best_move(int &u,int &v, vector<vector<int>> &delta, vector<int> &p);
   void perform_one_move(vector<vector<int>> &delta, vector<int> &p);
   void improve_qap_solution(
      int n,  // problem size  
      vector<int> &p,   // initial solution -> best solution  
      int &cost,  // -> cost of the best solution  
      int optimum,
      int nr_iterations);
   void print_vec(const vector<int> &v);
   void print_vec2d(const vector<vector<int>> &v);
   void check_sol(const vector<int> &sol, const int &cost);
};

#endif // QAPT600_H
