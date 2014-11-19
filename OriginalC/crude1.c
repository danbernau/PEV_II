
 //          - Source Treminal Network Unrealiability -
 //                     Crude Monte Carlo
 //  Leslie Murray - FCEIA, Universidad Nacional de Rosario, Argentina
 //                      - June, 2008 -

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"mt19937ar.c"

#define U genrand_real1()  // a random number (mersenne twister)

// structure link holds: the adjacency value (adj) in a scope of 1 or 0
// and the individual reliability (rlb) of the corresponding link. It
// also holds a value (smp) to be used whenever a network configuration
// is randomly sampled, being 1 if the link is operative and 0 otherwise.
// While (adj) is to build the adjacency matrix of the network, (smp)
// is to build the adjacency matrix of the network after random fails.
typedef struct link{
  int adj;
  double rlb;
  int smp;
}lnk, *pt_lnk;

// structure network holds: the source node (s) and the target node (t),
// the number of nodes (numnodes) and links (numlinks) and a matrix of
// links (I), allocated as I[numnodes+1][numnodes+1] just to be able to
// refer to the nodes by their real "name", i.e. avoiding the use of 0
// wich is reserved for some algoritms tricks
typedef struct network{
  int s;
  int t;
  int numnodes;
  int numlinks;
  struct link **I;
} net, *pt_net;

// adj is a matrix to hold the adjacency lists. Even when a waste of
// memory, every lists size is "numnodes"
int **list;

// An array to keep track of the visited nodes trough the DFS search
int *visited;

// A variable set to "1" if there is a path connecting "s" and "t"
// and 0 otherwise
int connected;

// "Initialize" takes the text file with the network info, as input,
// and returns a pointer to the network data structure (net)
pt_net Initialize(char* filename);

// "Phi" is the structure function of network nt. It returns 1 if
// nodes s and t are connected and 0 otherwise. It makes the
// determination after ramndom failure, i.e., over the (smp)
// adjacency matrix
int Phi(pt_net nt);

// "DFS" performs Depth First Search from node "s" and stops as soon
// as node "t" is visited
void DFS(int node, pt_net nt);

// "X" returns 1 if the link between node1 and node2 of network nt
// is operative, and 0 otherwise
int X(int node1, int node2, pt_net nt);

// "Fail" sets a random value of 0 or 1 in (smp) for all the links
// of network nt, i.e. it builds up a "failed" network and, upon
// these values, creates:
// - an instance of the "adjacency matrix"
// - an instance of the "adjacency lists"
void Fail(pt_net nt);

int seed;           // seed for the random numbers generator

/********************************************************************
* The executable "e" runs as:  ./e <File> <Size>, where:
* File: is the text file with the network (graph) info
* Size: is the number of Monte Carlo trials.
* Seed: is the seed for the random numbers generator
********************************************************************/
int main(int argc, char *argv[]){
  int i,j,k,S,size;
  double X,V,Q,t;
  pt_net n;
  clock_t  ti;


  if(argc<4){
    printf("\n Some input data is missing! Run as:\n");
    printf("\n ./executable <File> <Size> <Seed>\n\n");
    exit(1);
  }
  if(argc>4){
    printf("\n You've entered more data than necessary! Run as:\n\n");
    printf("\n ./executable <File> <Size> <Seed>\n\n");
    exit(1);
  }
  if((size=atoi(argv[2]))<1){
    printf("\n The number of trials can not be less than 1! Run as:\n");
    printf("\n ./executable <File> <Size> <Seed>\n\n");
    exit(1);
  }

  if((seed=atoi(argv[3]))<0){
    printf("\n The seed can not be negative! Run as:\n");
    printf("\n ./executable <File> <Size> <Seed>\n\n");
    exit(1);
  }

  n=Initialize(argv[1]);

  //---------------- THE CRUDE MONTE CARLO ALGORITHM ----------------
  ti = clock();

  S=0;
  X=0.0;
  V=0.0;
  for(k=0;k<size;k++){
    Fail(n);
    X=(1-Phi(n));
    S+=X;
    V+=X*X;
  }
  Q=(double)S/size;
  V=(V/size-Q*Q)/(size-1);

  t=(double)(clock()-ti)/CLOCKS_PER_SEC;
 //------------------ PRINT OF OUTPUT AND RESULTS ------------------
 printf("\n  CRUDE MONTE CARLO:");
 printf("\n  Network: %s   Replications: %d   ExecTime=%f",argv[1],size,t);
 printf("\n\n  Q =%1.16f = %1.2e  (Unreliability)",Q,Q);
 printf("\n  V =%1.16f = %1.2e  (Variance)",V,V);
 printf("\n  SD=%1.16f = %1.2e  (Standard Deviation)",sqrt(V),sqrt(V));
 printf("\n  RE=%1.16f = %1.2f%%     (Relative Error)",sqrt(V)/Q,100*sqrt(V)/Q);
 printf("\n  ------------------------------------------------------------\n\n");
 //-----------------------------------------------------------------

  return 0;

}




/********************************************************************
* Initialize reads from the file that holds the network data all in-
* formation necessary to build up and allocate the network data
* structure net. It also initializes the random numbers generator.
*********************************************************************/
pt_net Initialize(char* filename){

  pt_net pt_n;
  FILE* fp;
  int i,j,node1,node2,num;
  double reliability;

  // Allocate one unit of the structure net to hold the network info
  if((pt_n=(pt_net)calloc(1,sizeof(net)))==NULL){
    printf("\n Fail attempting to allocate memory...\n");
    exit(1);
  }

  // Open the file with the network info and scan for the source and
  // target node, the number of nodes and the number of links at the
  // top of it
  if((fp=fopen(filename,"r"))==NULL){
   printf("\n Fail attempting to open a disk file...\n");
    exit(1);
  }
  fscanf(fp,"%d",&pt_n->s);
  fscanf(fp,"%d",&pt_n->t);
  fscanf(fp,"%d",&pt_n->numnodes);
  fscanf(fp,"%d",&pt_n->numlinks);

  // Allocate matrix I with a size of (numnodes+1)*(numnodes+1), and link it
  // to the corresponding pointer of network net

  if((pt_n->I=(pt_lnk*)calloc(pt_n->numnodes+1,sizeof(struct link*)))==NULL){
    printf("\n Fail attempting to allocate memory...\n");
    exit(1);
  }
  for(i=0;i<=pt_n->numnodes;i++){
    if((pt_n->I[i]=(pt_lnk)calloc(pt_n->numnodes+1,sizeof(struct link)))==NULL){
      printf("\n Fail attempting to allocate memory...\n");
      exit(1);
    }
  }

  // Initialize matrix I with 0s for the adjacencies and 0.0s for the
  // reliabilities of every link
  for(i=1;i<=pt_n->numnodes;i++)
    for(j=1;j<=pt_n->numnodes;j++){
      pt_n->I[i][j].adj=0;
      pt_n->I[i][j].rlb=0.0;
      pt_n->I[i][j].smp=0;
    }

  // load matrix I values from the file with the network info and set
  // the corresponding values of adjacency and reliability
  for(i=1;i<=pt_n->numnodes;i++){
    fscanf(fp, "%d%d", &node1, &num);
    for(j=1;j<=num;j++){
      fscanf(fp, "%d%lf", &node2, &reliability);
        pt_n->I[node1][node2].adj=1;
        pt_n->I[node1][node2].rlb=reliability;
    }
  }

  fclose(fp);

  // Select a seed value and initialize the pseudorandom numbers
  // generator Mersenne Twister
  init_genrand(seed);


  // Allocate matrix lists of size (numnodes+1)*(numnodes+1), same
  // as matrix I
  if((list=(int**)calloc(pt_n->numnodes+1,sizeof(int*)))==NULL){
    printf("\n Fail attempting to allocate memory...\n");
    exit(1);
  }
  for(i=0;i<=pt_n->numnodes;i++){
    if((list[i]=(int*)calloc(pt_n->numnodes,sizeof(int)))==NULL){
      printf("\n Fail attempting to allocate memory...\n");
      exit(1);
    }
  }

  // Initialize matrix lists, filling it with "0"
  for(i=1;i<=pt_n->numnodes;i++)
    for(j=1;j<=pt_n->numnodes;j++){
      list[i][j]=0;
    }

  // Allocate the visited array and set 0 for all nodes (not visited)
  if((visited=(int*)calloc(pt_n->numnodes,sizeof(int)))==NULL){
    printf("\n Fail attempting to allocate memory...\n");
    exit(1);
  }

  for(i=1;i<=pt_n->numnodes;i++)
    visited[i]=0;

  connected=0;

  return pt_n;
}


/********************************************************************
* Phi is the structure function. Input D is the diameter value for
* the evaluation of Diameter Constrained Reliability and nt is a
* pointer to the network under consideration. Phi returns 1 if any
* pair of nodes of the network is linked by a path of lenght  D or
* less, and 0 otherwise (Phi uses the adjacency matrix after network
* fail: smp)
*********************************************************************/
int Phi(pt_net nt){

  int k;

  connected=0;

  for(k=1;k<=nt->numnodes;k++)
    visited[k]=0;

  // DFS() makes a Depth First Search from node "s" and:
  // - returns 1 if node "t" is reached
  // - returns 0 if node "t" is not reached

  DFS(nt->s, nt);

  return connected;
}


void DFS(int node, pt_net nt){

  int k;

  if(connected)
    return;

  if(node == nt->t){
    connected=1;
    return;
  }

  visited[node]=1;

  k=1;
  while(list[node][k]>0){
    if(!visited[list[node][k]])
      DFS(list[node][k],nt);
    k++;
  }
  
  return;
}

/********************************************************************
* Depending on the probability distribution of the link between node1
* and node 2, "X" returns 1 if such link is operative and 0 otherwise
* It is accepted that function "X" is only required for links that
* actuallly exists, i.e. links for which I[node1][node2].adj = 1
*********************************************************************/
int X(int node1, int node2, pt_net nt){

  if(U < nt->I[node1][node2].rlb)
    return 1;
  else
    return 0;
}

/********************************************************************
* Depending on the probability distribution of every link, "Fail"
* sets a random value of 0 or 1 in (smp) for every existing link of
* network nt. This way a randomly failed network is build up and saved
* into matrix I[][].smp and also in lists[][]
*********************************************************************/
void Fail(pt_net nt){

  int i,j,k;

  for(i=1;i<=nt->numnodes;i++)
    for(j=i+1;j<=nt->numnodes;j++)
      if(nt->I[i][j].adj){
        nt->I[i][j].smp=(nt->I[i][j].adj && X(i,j,nt));
        nt->I[j][i].smp=nt->I[i][j].smp;
      }
      // It is accepeted that the network graph is "undirected" and
      // that the reliability of every link has the same value in
      // both directions. It is therefore unnecesary to sample all
      // the matrix elements, it suffices to do it for all the
      // elements above the diagonal and to set the same value to
      // the simetric element (the diagonal elements are all 0)


// IMPRIMO LA MATRIZ DE ADYACENCIA
//  for(i=1;i<=nt->numnodes;i++){
//    for(j=1;j<=nt->numnodes;j++){
//      printf("%d ",nt->I[i][j].smp);
//    }
//    printf("\n");
//  }
// printf("\n");

  // LIMPIO LAS LISTAS DE ADYACENCIA ANTERIORES
  for(i=1;i<=nt->numnodes;i++){
    k=1;
    while(list[i][k]>0){
      list[i][k++]=0;
    }
  }

  // CREO LAS LISTAS DE ADYACENCIA
  for(i=1;i<=nt->numnodes;i++){
    k=1;
    for(j=1;j<=nt->numnodes;j++){
      if(nt->I[i][j].smp){
        list[i][k++]=j;
      }
    }
  }

// IMPRIMO LAS LISTAS DE ADYACENCIA
//  for(i=1;i<=nt->numnodes;i++){
//    for(j=1;j<=nt->numnodes;j++){
//      printf("%d ",list[i][j]);
//    }
//    printf("\n");
//  }

  return;
}



