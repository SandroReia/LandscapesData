#include<cmath>
#include<iostream>
#include<string.h>
#include<cstring>
#include<stdio.h>
#include<stdlib.h>
#include<cstdio>
#include<iomanip> 				
#include<fstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort_vector.h>
#include <ios>
#include <errno.h>
#include <stdbool.h>

#include <unistd.h>

using namespace std;

#define PI 3.1415927

#define NIL (0)    

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

long idum;
gsl_rng *gerador;


struct populacao{

  int N;
  double mut;
  int *step;
  int nstep;
  int *step_dom;
  int nstep_dom;
  double p_imit;
  int *fixpos;
  int *pilhapos;
  int nstep_real;
  int *step_real;
};

struct caminhada{

  int nstep;
  double fitness_final;
  int sequence_final;
};



struct sequencia{
  unsigned int ind;
  int L;
  int ind_min;
  int ind_max;
  int dham;
  int Nmax;
  int *max;
  double roughness;
  double rough_local;
  int *path;
  int path_size;
  double *w;
  double *epsilon;
  double fmax;
  double roughness_path;
  double *fitness;
  int *z;
  int **C;
  double fmin;
  char **seq;
  int *indice_order;
  int *rank;
};

void fitnessinit(struct sequencia *sequence, struct populacao *popul);
unsigned int pot(int a1, int a2);
void landscape(struct sequencia *sequence, struct populacao *popul);
void landscape_stat(struct sequencia *sequence, struct populacao *popul);
void dynamics(struct sequencia *sequence, struct populacao *popul, struct caminhada *walk);
double Err(double aa,int bb);
double fitness(int ind, struct sequencia *sequence);
unsigned int hamming_distance(unsigned int sequence1, unsigned  int sequence2);
double path_divergence(int* path1, int npath1, int* path2, int npath2, struct sequencia *sequence);
double average_path_divergence(int** paths, int* path_size, int npaths, int npath_pairs, double *W, struct sequencia *sequence);

int i1=0;

int main(int ac, char **av)
{
//   ofstream arqteste;
//   arqteste.open("test.dat");
  FILE *ptt1, *ptt2, *ptt3;

  sequencia sequence;

  populacao popul;

  caminhada walk;

  int nwalks, max1, cont_path_total, cont_path[50], cont_path_pred[50], *cont_newpath, i, j, k, cont_walk, verif, **nstep_path, **predictability, sum, ***step_path, sum_step, indice_walk, **step_path_div, *nstep_path_div,
    sum_path, contador_caminho[50], sum_length[200], ind_maior, minimum_path, cont_minimum, lim_inf, nwalks_min;

  double soma1, *P2, *path_divergence, med_step, *W, *accessibility, med_length[200], maior;

  long seed;  //define a seed

  char arq1[1000], arq2[1000], arq3[1000];

  if (ac!=2)
    {
      cout  <<  "start the program like this:\n" << av[0]
            << " <N> <U> <s> <sd> <alpha> <Nlevel> <tmax> <Tterm> <config> <Nbins> <Nmax> <semente> \n"
            << endl;
      exit (-1);
    }

  j=0;
  
  nwalks_min = atoi (av[++j]);

  cout  << "#invocation: ";
  for (int i=0; i<ac; i++){
    cout << av[i] << " ";
  }
  cout << endl;

  seed = time (NULL) * getpid();  //set the seed to system time
  gerador=gsl_rng_alloc(gsl_rng_mt19937);  
  gsl_rng_set(gerador, seed); //give the seed to random generator
  
  srand(time(NULL));

  max1 = 160000;

  sequence.fitness = new double[160000];

  sequence.z = new int[160000];

  sequence.rank = new int[160000];

  sequence.indice_order = new int[160000];

  sequence.C = new int*[160000];
  for( i=0; i<160000; i++ )
    sequence.C[i] = new int[100];

  sequence.seq = new char*[160000];
  for( i=0; i<160000; i++ )
    sequence.seq[i] = new char[4];

  sequence.max = new int[1000];

  landscape(&sequence,&popul);

  landscape_stat(&sequence,&popul);
  cout << sequence.Nmax << endl;
  
  
  nstep_path = new int*[sequence.Nmax];
  for( i=0; i<sequence.Nmax; i++ )
    nstep_path[i] = new int[100001];
  
  step_path = new int**[sequence.Nmax];
  for( i=0; i<sequence.Nmax; i++ )
    step_path[i] = new int*[100001];
  for( i=0; i<sequence.Nmax; i++ )
    for( j=0; j<100001; j++ )
      step_path[i][j] = new int[100];

  nstep_path_div = new int[100001];

  step_path_div = new int*[100001];
  for( j=0; j<100001; j++ )
    step_path_div[j] = new int[100];
 
  predictability = new int*[sequence.Nmax];
  for( i=0; i<sequence.Nmax; i++ )
    predictability[i] = new int[100001];

  W = new double[100001];

  sequence.path = new int[100];

  popul.step_real = new int[100];

  cont_newpath = new int[sequence.Nmax];

  // cont_path = new int[sequence.Nmax];

  P2 = new double[sequence.Nmax];

  accessibility = new double[sequence.Nmax];

  path_divergence = new double[sequence.Nmax];

  cont_path_total = 0;
  for( k=0; k<sequence.Nmax; k++ )
    {
      cont_newpath[k] = 0;
      cont_path[k] = 0;
      cont_path_pred[k] = 0;
      sum_length[k] = 0;
    }

  for( k=0; k<sequence.Nmax; k++ )
    contador_caminho[k] = 0;
  
  cont_walk = 0;
  sum_step = 0;

  cout << sequence.Nmax << endl;

  minimum_path = 0;
  nwalks = 0;
  while( (minimum_path<nwalks_min) )
    {
      nwalks++;
      if( (nwalks%1000000)==0 )
      	cout << nwalks << endl;
      
      fitnessinit(&sequence,&popul);
      
      dynamics(&sequence,&popul,&walk);

      cont_path_total++;
      
      sum_step += popul.nstep_real;

      for( k=0; k<sequence.Nmax; k++ )
	if( walk.sequence_final==sequence.max[k] )
	  {
	    indice_walk = k;
	    break;
	  }
      
      cont_path[indice_walk]++;

      sum_length[indice_walk] += popul.nstep_real;

      cont_minimum = 0;
      for( k=0; k<sequence.Nmax; k++ )
	cont_minimum += (cont_path[k]>0);

      if( cont_minimum==12 )
	{
	  minimum_path = 10000000;      
	  for( k=0; k<sequence.Nmax; k++ )
	if( (cont_path[k]<minimum_path) && (cont_path[k]>0))
	  minimum_path = cont_path[k];
	  cout << minimum_path << endl;
	}
      /*   sprintf(arq1,"GB1-DISTLENGTH-nwalks%d.dat",nwalks);
      ptt1 = fopen(arq1,"a");
      fprintf(ptt1,"%d,%d\n",indice_walk,popul.nstep_real);
      fclose(ptt1);*/

      if( cont_newpath[indice_walk]<=100000 )
	{
	  cont_path_pred[indice_walk] = cont_path[indice_walk];
	  
	  verif = 1;
	  for( i=0; i<cont_newpath[indice_walk]; i++ )
	    {
	      if( nstep_path[indice_walk][i]==popul.nstep_real )
		{
		  sum = 0;
		  for( j=0; j<=popul.nstep_real; j++ )
		    {
		      if( popul.step_real[j]==step_path[indice_walk][i][j] )
			sum++;
		    }
		  if( sum==(popul.nstep_real+1) )
		    {
		      predictability[indice_walk][i]++;
		      verif = 0;
		      break;
		    }
		}
	    }
	  if( verif==1 )
	    { 
	      nstep_path[indice_walk][cont_newpath[indice_walk]] = popul.nstep_real;
	      for( k=0; k<=popul.nstep_real; k++ )
		{
		  step_path[indice_walk][cont_newpath[indice_walk]][k] = popul.step_real[k];
		  sequence.path[k] = step_path[indice_walk][cont_newpath[indice_walk]][k];
		}
	      sequence.path_size = popul.nstep_real;
	      predictability[indice_walk][cont_newpath[indice_walk]] = 1;
	      cont_newpath[indice_walk]++;     
	    }
	}
      
    }


 
  for( k=0; k<sequence.Nmax; k++ )
    {
      soma1 = 0;
      for( i=0; i<cont_newpath[k]; i++ )
	{
	  W[i] = (double)predictability[k][i]/cont_path_pred[k];
	  soma1 += W[i];

	  if( W[i]>maior )
	    {
	      maior = W[i];
	      ind_maior = i;
	    }
	}

      P2[k] = 0;
      for( i=0; i<cont_newpath[k]; i++ )
	P2[k] += W[i]*W[i];
      
      accessibility[k] = (double)cont_path[k]/nwalks;
      
      for( i=0; i<cont_newpath[k]; i++ )
	{
	  for( j=0; j<=nstep_path[k][i]; j++ )
	    step_path_div[i][j] = step_path[k][i][j];
	  nstep_path_div[i] = nstep_path[k][i];
	}
      
      path_divergence[k] = average_path_divergence(step_path_div, nstep_path_div, cont_newpath[k], cont_newpath[k], W, &sequence);	

      med_length[k] = (double)sum_length[k]/cont_path[k];

      if( cont_newpath[k]>0 )
	{

	  int k1;
	  int n_size = cont_newpath[k];
	  size_t i;
	  size_t n = cont_newpath[k];
	  gsl_vector * v = gsl_vector_alloc(n);
	  
	  gsl_permutation * perm = gsl_permutation_alloc(n);
	  gsl_permutation * rank = gsl_permutation_alloc(n);
	  
	  for( k1=0; k1<n_size; k1++ )
	    gsl_vector_set(v,k1,W[k1]);
	  
	  gsl_sort_vector_index (perm, v);
	  gsl_permutation_inverse (rank, perm);
	  
	  for (i = 0; i < n; i++)
	    sequence.rank[i] = rank->data[i];
	  
	  for (i = 0; i < n; i++)
	    sequence.indice_order[sequence.rank[i]] = i;
	   
	  gsl_permutation_free (perm);
	  gsl_permutation_free (rank);
	  gsl_vector_free (v);
	  
	  sprintf(arq1,"FREQ-GB1-FREQUENCY-nwalks%d.dat",nwalks);
	  ptt1 = fopen(arq1,"a");
	  lim_inf = n - 500;
	  if( lim_inf<0 )
	    lim_inf = 0;
	  for( i=lim_inf; i<n; i++ )
	    {
	      fprintf(ptt1,"%d    %d    %g  \t",k,sequence.indice_order[i],W[sequence.indice_order[i]]); 
	      for( j=0; j<=nstep_path[k][sequence.indice_order[i]]; j++ )
		fprintf(ptt1,"%c%c%c%c  \t",sequence.seq[step_path[k][sequence.indice_order[i]][j]][0],sequence.seq[step_path[k][sequence.indice_order[i]][j]][1],sequence.seq[step_path[k][sequence.indice_order[i]][j]][2],sequence.seq[step_path[k][sequence.indice_order[i]][j]][3]);
	      fprintf(ptt1,"\n");
	    }
	  fclose(ptt1);
	}

    }
  
  med_step = (double)sum_step/nwalks;
  
 sprintf(arq1,"FREQ-GB1-LENGTH-nwalks%d.dat",nwalks);
  ptt1 = fopen(arq1,"a");
  fprintf(ptt1,"%g   \t",med_step);
  for( k=0; k<sequence.Nmax; k++ )
    fprintf(ptt1,"%g   \t",med_length[k]);
  fclose(ptt1);
  
  sprintf(arq1,"FREQ-GB1-Predictability-nwalks%d.dat",nwalks);
  ptt1 = fopen(arq1,"a");
  for( k=0; k<sequence.Nmax; k++ )
    fprintf(ptt1," %g \t",P2[k]);
  fclose(ptt1);
  
  
  sprintf(arq1,"FREQ-GB1-Accessibility-nwalks%d.dat",nwalks);
  ptt1 = fopen(arq1,"a");
  for( k=0; k<sequence.Nmax; k++ )
    fprintf(ptt1,"%g \t",accessibility[k]);
  fclose(ptt1);
  

  sprintf(arq1,"FREQ-GB1-Divergence-nwalks%d.dat",nwalks);
  ptt1 = fopen(arq1,"a");
  for( k=0; k<sequence.Nmax; k++ )
    fprintf(ptt1,"%g \t",path_divergence[k]);
  fclose(ptt1);

  sprintf(arq1,"FREQ-GB1-LocalOptimum-nwalks.dat");
  ptt1 = fopen(arq1,"a");
  for( k=0; k<sequence.Nmax; k++ )
    fprintf(ptt1,"%g \t",sequence.fitness[sequence.max[k]]);
  fclose(ptt1);


  sprintf(arq1,"FREQ-Relacao-GB1.dat");
  ptt1 = fopen(arq1,"a");
  for( k=0; k<sequence.Nmax; k++ )
    fprintf(ptt1,"%d->%c%c%c%c \n",k,sequence.seq[sequence.max[k]][0],sequence.seq[sequence.max[k]][1],sequence.seq[sequence.max[k]][2],sequence.seq[sequence.max[k]][3]);
    fclose(ptt1);

  
}



void landscape(sequencia *sequence, populacao *popul)
{

  FILE *fp;

  char seq[10];

  int seq_numb, dist, indice, indice1, degree, neighbor;

  int A[5], j;

  double fitness;

  fp = fopen("elife_seq_number.txt", "r");
  j = 0;

  while (!feof(fp))
    {
      fscanf(fp,"%d\t%c%c%c%c\n",&seq_numb,&seq[0],&seq[1],&seq[2],&seq[3]);

      for( j=0; j<4; j++ )
	sequence->seq[seq_numb][j] = seq[j];
      //      printf("%d\n", seq_numb);
      // printf("%c,%c,%c,%c,%c,%c,%c,%c,%c\n",
      //seq[j][0],seq[j][1],seq[j][2],seq[j][3],seq[j][4],seq[j][5],seq[j][6],seq[j][7],seq[j][8]);
      // j++;
    }

  /*  dist = 0;
  for( j=0; j<9; j++ )
    dist += (seq[637][j]!=seq[638][j]);
    cout << "dist=" << dist << endl;*/


  ifstream arq_leit("elife_sequence_fitness.txt");
  while( !arq_leit.eof() )
    {      
      arq_leit >> indice >> fitness;
      sequence->fitness[indice] = fitness;
    }

  ifstream arq_leit1("elife_sequence_degree.txt");
  while( !arq_leit1.eof() )
    {      
      arq_leit1 >> indice >> degree;
      sequence->z[indice] = degree;
    }

  ifstream arq_leit2("elife_sequence_neighbors_correta.txt");
  while( !arq_leit2.eof() )
    {      
      arq_leit2 >> indice >> indice1 >> neighbor;
      sequence->C[indice][indice1] = neighbor;
    }

  
}


void fitnessinit(sequencia *sequence, populacao *popul)
{

  sequence->ind = 0;

}  

 void dynamics(sequencia *sequence, populacao *popul, caminhada *walk)
{

  int cont, i, j, sum, indice, ind, nmaior, stop, nstep, k;
  
  int *pilha;
  
  double fmaior;

  pilha = new int[100];

  ind = sequence->ind;

  nmaior = 0;
  fmaior = 0;

  stop = 0;
  nstep = 0;

  popul->step_real[0] = ind;
  while(stop==0)
    { 
      nmaior = 0;
      for( j=0; j<sequence->z[ind]; j++ )
	{
	  if( sequence->fitness[sequence->C[ind][j]]>sequence->fitness[ind] )
	    {
	      pilha[nmaior] = sequence->C[ind][j];
	      nmaior++;
	    }
	}
   
      if( nmaior>0 )
	{
	  nstep++;
	  ind = (int)(gsl_ran_flat(gerador, 0., 1.)*nmaior);
	  ind = pilha[ind];
	  popul->step_real[nstep] = ind;
	}
      else
	{
	  stop = 1;
	  walk->nstep = nstep;
	  walk->fitness_final = sequence->fitness[ind];
	  walk->sequence_final = ind;
	  popul->nstep_real = nstep;
	}
    }

  delete[] pilha;
  
}

double fitness(int ind, sequencia *sequence)
{
  int k, config, indice, segmento1;
  double fitness;

  fitness = sequence->fitness[ind];

  return(fitness);
}


void landscape_stat(sequencia *sequence, populacao *popul)
{

  int max1, ind, ind_min, ind_max, Nmax, verif, k;

  double fit_min, fit_max, fitness;
  
  max1 = 160000;

  fit_min = 100000;
  fit_max = 0;

  for( ind=0; ind<max1; ind++ )
    {
      fitness = sequence->fitness[ind];
      
      if( fitness<fit_min )
	{
	  fit_min = fitness;
	  ind_min = ind;
	}
      if( fitness>fit_max )
	{
	  fit_max = fitness;
	  ind_max = ind;
	}
    }

  sequence->ind_min = ind_min;
  sequence->ind_max = ind_max;

  sequence->fmax = fit_max;
  sequence->fmin = fit_min;
  
  Nmax = 0;
  for( ind=0; ind<max1; ind++ )
    {
      verif = 0;
      for( k=0; k<sequence->z[ind]; k++ )
	if( sequence->fitness[sequence->C[ind][k]]>sequence->fitness[ind] )
	  {
	    verif = 1;
	    break;
	  }

      if( verif==0 )
	{
	  sequence->max[Nmax] = ind;
	  Nmax++;
	}
    }
  sequence->Nmax = Nmax;
  
}

double Err(double aa,int bb)
{
  double erro;
  
  erro = pow((double)aa/bb,0.5);
  
  return(erro);
}




unsigned int hamming_distance(unsigned int sequence1, unsigned  int sequence2)
{
	unsigned int dist = 0;
	unsigned int val  = sequence1^sequence2;
	
	while(val)
	{
		++dist;
		val &= (val-1);
	}
	
	return dist;
}


double path_divergence(int* path1, int npath1, int* path2, int npath2, sequencia *sequence) 
{
  int calc_dist, k;
    int distance = 0;
    for(int i=0; i<npath1; i++)
    {
      calc_dist = 0;
      for( k=0; k<4; k++ )
	calc_dist += (sequence->seq[path1[i]][k]!=sequence->seq[path2[0]][k]);
      int min_dist = calc_dist;
      for(int j=0; j<npath2; j++)
        {
	  if(min_dist == 0) break;

	  calc_dist = 0;
	  for( k=0; k<4; k++ )
	    calc_dist += (sequence->seq[path1[i]][k]!=sequence->seq[path2[j]][k]);
	  
	  int dist = calc_dist;
	  min_dist = (min_dist < dist)?min_dist:dist;
        }
      distance += min_dist;
    }
    return double(distance)/npath1;
}

#include <algorithm>
double average_path_divergence(int** paths, int* path_size, int npaths, int npaths_maximum, double *W, sequencia *sequence)
{
    //Choose the smallest between the number of paths and the maximum
    npaths = (npaths<npaths_maximum)?npaths:npaths_maximum;
    
    //Take a random order for the paths (we will then consider the first npaths)
    int arr[npaths];
    for(int i = 0; i < npaths; ++i) arr[i] = i;
    random_shuffle(arr, arr+npaths);
    
    double divergence = 0.0;
    for(int i=0; i<npaths; i++) for(int j=0; j<npaths; j++)
    {
        if(i != j)
	  divergence += W[arr[i]]*W[arr[j]]*path_divergence(paths[arr[i]], path_size[arr[i]], paths[arr[j]], path_size[arr[j]], sequence);
    }
    
    return divergence;
}



unsigned int pot(int a1, int a2)
{
   int i;

   unsigned int produto;

   produto = 1;
   for( i=1; i<=a2; i++ )
      produto = produto*a1;
   return(produto);
}

