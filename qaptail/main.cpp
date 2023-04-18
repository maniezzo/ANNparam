#include "QAPT600.h"

int main(int argc, char** argv) 
{  int cost,n,maxiter;
   clock_t start, end;
   string instance;
   double minTLcoef,maxTLcoef;

   // initialize random seed: 
   srand(550); // (time(NULL));

   if(argc > 1)
      /* IRACE PARAMETERS
      :: %%1 is the candidate configuration number
      :: %%2 is the instance ID
      :: %%3 is the seed
      */
      instance = argv[4];           // --inst instance name
   else
      instance="instances//TAI15A.DAT";

   QAPT* QAP = new QAPT(instance);
   n = QAP->n;

   if(argc > 1)
   {
      maxiter  = n*n*atoi(argv[8]);     // --maxiter
      minTLcoef = atof(argv[10]);        // --mintl  min tabu list length
      QAP->lower_taboo_list_size = trunc(minTLcoef*n);
      maxTLcoef = atof(argv[12]);       // --maxtl  max tabu list length
      QAP->higher_taboo_list_size = trunc(maxTLcoef*n) + 4;
      QAP->nr_iteration_before_aspiration = atoi(argv[14])*n*n;   // --aspir
      QAP->nr_iter_resizeTL = atoi(argv[16])*QAP->higher_taboo_list_size; // --resize iterations before tabu list resizing
      QAP->max_time = atoi(argv[18]);   // --maxt max cpu time, in sec
      QAP->isVerbose = false;            // --verbose
      QAP->nrep = atoi(argv[20]);       // --nrep
      if(QAP->isVerbose)
         for(int i=0;i<argc;i++)
            cout << i <<") "<<argv[i]<<endl;
   }
   else
   {  QAP->isVerbose = true;
      QAP->lower_taboo_list_size = trunc(0.78*n);
      QAP->higher_taboo_list_size = trunc(2.12*n) + 4;
      QAP->nr_iteration_before_aspiration = 4*n*n;
      QAP->nr_iter_resizeTL = 5*QAP->higher_taboo_list_size; // iterations before tabu list resizing
      QAP->max_time = 10; // max cpu time, in sec
      QAP->nrep = 2;
      maxiter = 440;
   }

   start = clock();
   cost = QAP->run_ts(maxiter); 
   end  = clock();
   double ttot = (double)(end-start)/CLOCKS_PER_SEC;
   cout << cost;

   if(QAP->isVerbose)
   {  cout << "\n<ENTER> to exit"; 
      getchar();
   }

   delete QAP;
}