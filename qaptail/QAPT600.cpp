// Taboo search procedure for asymmetrical quadratic assignment problems  
// Author : Eric Taillard  
// Date : 94/11/23  
// Implementation used in : E. Taillard, "Comparison of iterative searches  
// for the quadratic assignment problem", Publication 989, Center for  
// research on transportation, University of Montreal, Montreal, 1994  
#include "QAPT600.h"

// ctor
QAPT::QAPT(string instance)
{
   read_problem(instance,0);
}

//dtor
QAPT::~QAPT()
{
}

int QAPT::run_ts(int nb_iterations)
{
   int i,irep;


   init_sol_seed = 123456789;
   best_cost= infinite;
   cpu1 = clock();
   for(irep=0;irep<nrep;irep++)
   {
      cost = infinite;
      for(int i=0; i<n; i++)
      {  taboo_list.push_back(vector<int>());
         for(int j=0; j<n; j++)
            taboo_list[i].push_back(0);
      }

      // initial solution
      for (i=0;i<n;i++) solution.push_back(i);
      // random solution
      for (i=0;i<n-1;i++)
         swap(solution[i], solution[unif(init_sol_seed, i, n-1)]);

      // the search
      improve_qap_solution(n, solution, cost, optimum, nb_iterations);
   }

   if(isVerbose)
   {  cpu=(double)(clock()-cpu1)/CLOCKS_PER_SEC;   
      cout << "cpu=" << cpu << " cost:" << best_cost << endl;
   }

   return best_cost;
}

void QAPT::initialize(vector<vector<int>> &delta, vector<int> &p)
{
   int i,j;

   cost = 0;
   
   for (i = 0; i<n; i++)
      for (j = 0;j<n;j++) 
      {  cost = cost + a[i][j] * b[p[i]][p[j]];
         if (i < j)
            delta_full_computation(i,j,delta,p);
      }
  
   best_cost = cost;
   taboo_list_seed = 123456789;
   current_taboo_list_size = 
      unif(taboo_list_seed, lower_taboo_list_size, higher_taboo_list_size);
}

int QAPT::unif(int &seed, int low, int high)
{
   const int m = 2147483647, a = 16807, b = 127773, c = 2836;
   int  kl,res;
   double value_0_1; // real double precision (coded on 64 bits)  

   res = rand() % (high-low+1) + low;
   goto lend;

   kl = seed / b ;
   seed = a * (seed % b) - kl * c ;
   if (seed < 0) seed = seed + m ;
   value_0_1 = seed/m;
   res = low + trunc((high - low + 1) * value_0_1);
lend:
   return res;
};

void QAPT::read_problem(string instance, bool isVerbose)
{
   int i,j;
   string path = instance;
   ifstream fin;
   fin.exceptions ( std::ifstream::failbit | std::ifstream::badbit ); 
   string line;
   if(isVerbose)
      cout << "INSTANCE: "<<path<<endl;
   try 
   {  fin.open(path);
                  
      getline(fin,line);
      optimum = atoi(line.c_str() );
      getline(fin,line);
      getline(fin,line);
      n = atoi(line.c_str() );
      if(isVerbose)
        cout << "n:" << n << " opt:" << optimum << endl;
      getline(fin,line);

      // initialize distance and flow matrices
      vector<int> myvector;
      for(j = 0; j<n; j++)
          myvector.push_back(0);
      for(i = 0; i<n; i++)
      {
          a.push_back(myvector);
          b.push_back(myvector);
      }

      // read distance and flow matrices
      for(i = 0; i<n; i++)
      {  for(j = 0;j<n;j++)
          {  fin >> a[i][j];
            if(isVerbose) cout << a[i][j] << ' ';  
          };
          if(isVerbose) cout << endl;
      }
      getline(fin,line);
      if(isVerbose) cout << endl;

      for (i = 0;i<n;i++) 
      {  for(j = 0;j<n;j++) 
          {  fin >> b[i][j];
            if(isVerbose) cout << b[i][j] << ' ';  
          };
          if(isVerbose) cout << endl;
      }

      fin.close();  
   }
   catch (const ifstream::failure& e)
   {  cout << "ERROR opening/reading file";
   }
}

void QAPT::swap(int &a, int &b)
{
   int temp;
   temp = a; a = b; b = temp;
}

// computes the difference of solution values if units u and v are permuted  
void QAPT::delta_full_computation(int i, int j, vector<vector<int>> &delta, vector<int> &p)
{
   int k, sum;

   sum = 0;
   for (k=0;k<n;k++)
      if ( (k!=i) && (k!=j) )
         sum = sum + a[k][i]*(b[p[k]][p[j]]-b[p[k]][p[i]]) +
                     a[k][j]*(b[p[k]][p[i]]-b[p[k]][p[j]]) + 
                     a[i][k]*(b[p[j]][p[k]]-b[p[i]][p[k]]) +
                     a[j][k]*(b[p[i]][p[k]]-b[p[j]][p[k]]);

   sum = sum + a[i][i]*(b[p[j]][p[j]]-b[p[i]][p[i]]) +
               a[i][j]*(b[p[j]][p[i]]-b[p[i]][p[j]]) +
               a[j][i]*(b[p[i]][p[j]]-b[p[j]][p[i]]) +
               a[j][j]*(b[p[i]][p[i]]-b[p[j]][p[j]]);
   delta[i][j] = sum;
}

// idem but needs the value of delta[u,v] before the exchange of i and j   
void QAPT::delta_short_computation(int r,int s,int i,int j, vector<vector<int>> &delta, vector<int> &p)
{
  delta[i][j] = delta[i][j] + 
        (a[r][i]-a[r][j]+a[s][j]-a[s][i])*(b[p[s]][p[i]]-b[p[s]][p[j]]+
        b[p[r]][p[j]]-b[p[r]][p[i]]) +
        (a[i][r]-a[j][r]+a[j][s]-a[i][s])*(b[p[i]][p[s]]-b[p[j]][p[s]]+
        b[p[j]][p[r]]-b[p[i]][p[r]]);
}

void QAPT::find_best_move(int &u,int &v, vector<vector<int>> &delta, vector<int> &p)
{
  int i, j, delta_min;
  bool aspired, taboo;

  delta_min = infinite;
  u = infinite; // in case no allowed moves exist  
  aspirating_iteration = current_iteration - nr_iteration_before_aspiration;
  
  for (i=0;i<n-1;i++) 
    for (j=i+1;j<n;j++) 
    {
      taboo = (taboo_list[i][p[j]] >= current_iteration) &&
              (taboo_list[j][p[i]] >= current_iteration);

      aspired = ( (taboo_list[i][p[j]] < aspirating_iteration) &&
                  (taboo_list[j][p[i]] < aspirating_iteration) ) ||
                ((cost + delta[i][j] < best_cost) && taboo);

      if ( ((delta[i][j] < delta_min) && !taboo) || aspired ) 
      {
        u = i; v = j;
        if (aspired) delta_min = -infinite;
        else 
           delta_min = delta[i][j];
      }
    }
};

void QAPT::perform_one_move(vector<vector<int>> &delta, vector<int> &p)
{
   int i, j, u, v;

   find_best_move(u,v,delta,p);

   if (u != infinite) 
   {  cost = cost + delta[u][v];
      taboo_list[u][p[u]] = current_iteration + current_taboo_list_size;
      taboo_list[v][p[v]] = current_iteration + current_taboo_list_size;
      swap(p[u], p[v]);
   }
   else return; //writeln('no allowed move !');  
   check_sol(p,cost);

   if (cost < best_cost)
   {
      best_cost = cost;
      best_solution = p;
      cpu2 = clock();
      if(isVerbose) 
         cout << ">> iter "<<current_iteration<<" new cost "<<best_cost<<" t.cpu "<<cpu2 - cpu1<<endl;
   }

   for (i = 0;i<n-1;i++) 
      for (j = i+1;j< n;j++)
         if ( (i!=u) && (i!=v) && (j!=u) && (j!=v) )
            delta_short_computation(u,v, i,j,delta,p);
         else 
            delta_full_computation(i,j,delta,p);

   if (current_iteration % (nr_iter_resizeTL) == 0)
      current_taboo_list_size = unif(taboo_list_seed, lower_taboo_list_size, higher_taboo_list_size);
}

void QAPT::improve_qap_solution(
     int n,             // problem size  
     vector<int> &p,    // initial solution -> best solution  
     int &cost,         // -> cost of the best solution  
     int optimum,
     int nb_iterations)
{
   double tcpu;

   vector<vector<int> > delta(n, vector<int>(n, 0));  // value of the difference between two solutions
   // improve qap solution  
   initialize(delta,p);
   
   current_iteration = 0;
   do
   {  current_iteration++;
      perform_one_move(delta,p);
      tcpu= (double) (clock()-cpu1)/CLOCKS_PER_SEC;
      //cout << "iter "<<current_iteration<<" best cost:"<<best_cost<<" tcpu:"<<tcpu<<endl;
   } while ( (best_cost > optimum) && (tcpu<max_time) && (current_iteration < nb_iterations));

   if(isVerbose) 
      cout << "Final iter "<<current_iteration<<" best cost:"<<best_cost<<" tcpu:"<<tcpu<<endl;
   cost = best_cost;
   p = best_solution;
}

void QAPT::print_vec(const vector<int> &v)
{
  for(unsigned int i=0;i<v.size();i++)
    cout << v[i] << " ";
  cout << endl;
}

void QAPT::print_vec2d(const vector<vector<int>> &v)
{
  for(unsigned int i=0;i<v.size();i++)
  { for(unsigned int j=0;j<v[i].size();j++)
      cout << v[i][j] << " ";
    cout << endl;
  }
}

void QAPT::check_sol(const vector<int> &sol, const int &cost)
{
  int i,j,sum = 0;
  for(i=0;i<n;i++)
    for(j=i;j<n;j++)
        sum += a[i][i]*b[sol[j]][sol[j]] +
               a[i][j]*b[sol[j]][sol[i]] +
               a[j][i]*b[sol[i]][sol[j]] +
               a[j][j]*b[sol[i]][sol[i]];
  if(sum!=cost)
    cout << "********** ALERT!! cost mismatch: cost="<<cost<<" sum="<<sum<<endl;
}