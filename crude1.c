#define U genrand_real1()  \
\

/*1:*/
#line 120 "./crude1.w"


/*2:*/
#line 140 "./crude1.w"


#include<stdio.h>  
#include<stdlib.h>  
#include<math.h>  
#include<time.h>  
#include"mt19937ar.c" 

/*:2*/
#line 122 "./crude1.w"
;
/*4:*/
#line 162 "./crude1.w"


typedef struct link{
	int adj;

	double rlb;

	int smp;

}lnk,*pt_lnk;

typedef struct network{
	int s;

	int t;

	int numnodes;

	int numlinks;

	struct link**I;

}net,*pt_net;

int**list;


int*visited;


int connected;


int seed;

double breadth = 2.0;
double l;
double beta;
int indicator;


/*:4*/
#line 123 "./crude1.w"
;
/*5:*/
#line 211 "./crude1.w"


pt_net Initialize(char*filename);


int X(int node1,int node2,pt_net nt);


void Fail(pt_net nt);


void DFS(int node,pt_net nt);


int Phi(pt_net nt);


/*:5*/
#line 124 "./crude1.w"
;
int main(int argc,char*argv[])
{
/*6:*/
#line 235 "./crude1.w"

/*Start of actual implementation*/
	int i,j,k,size,X;
	double V,Q,t,S;
	pt_net n;
	clock_t ti;

/*:6*/
#line 127 "./crude1.w"
	;
/*7:*/
#line 254 "./crude1.w"


	if(argc<4){
		printf("\n Some input data is missing! Run as:\n");
		printf("\n ./executable <File> <Size> <Seed>\n\n");
		exit(1);
	}
	if(argc> 4){
		printf("\n You've entered more data than necessary! Run as:\n\n");
		printf("\n ./executable <File> <Size> <Seed>\n\n");
		exit(1);
	}
	if((size= atoi(argv[2]))<1){
		printf("\n The number of trials can not be less than 1! Run as:\n");
		printf("\n ./executable <File> <Size> <Seed>\n\n");
		exit(1);
	}

	if((seed= atoi(argv[3]))<0){
		printf("\n The seed can not be negative! Run as:\n");
		printf("\n ./executable <File> <Size> <Seed>\n\n");
		exit(1);
	}

/*:7*/
#line 128 "./crude1.w"
	;
/*8:*/
#line 289 "./crude1.w"


	n= Initialize(argv[1]);
/*minimum cut set in relation to nodes in network -> how critical is the min cut set in relation to overall network size?*/
	ti= clock();

	S= 0.0;
	X= 0.0;
	V= 0.0;

/* start of monte carlo logic*/
	for(k= 0;k<size;k++){
/* reset l for current monte carlo trial*/
			l = 1.0;

/* set link random fails and builds the adjacency lists after that fail */
		Fail(n);

/* Phi: if node t is reached Phi ( ) returns 1, otherwise it returns 0 */
/* Deviation metric*/
		X= (1-Phi(n));
		if (X==1){
			S+=l;
			indicator++;
		}
		V+= l*l;
		//		V+= X*X;
	}
/* compare successes to size as number of monte carlo trials*/
/* sum/N*/
	Q= (double)S/size;
/*estimator of the variance*/
	V= (V/size-Q*Q)/(size-1);

	t= (double)(clock()-ti)/CLOCKS_PER_SEC;


/*:8*/
#line 129 "./crude1.w"
	;
/*9:*/
#line 319 "./crude1.w"


	printf("\n  Network: %s   Replications: %d   ExecTime=%f",argv[1],size,t);
	printf("\n  ******************* CRUDE MONTE CARLO (1) ******************");
	printf("\n       Unreliability    Q = %1.16f = %1.2e",Q,Q);
	printf("\n       Success    	S = %d",indicator);
	printf("\n       Variance         V = %1.16f = %1.2e",V,V);
	printf("\n       Std. Dev.       SD = %1.16f = %1.2e",sqrt(V),sqrt(V));
	printf("\n       Relative Error  RE = %1.16f = %1.2f%%",sqrt(V)/Q,
		100*sqrt(V)/Q);
	printf("\n  ************************************************************\n\n");



/*:9*/
#line 130 "./crude1.w"
	;
	return 0;
}
/*10:*/
#line 336 "./crude1.w"

/*11:*/
#line 345 "./crude1.w"

pt_net Initialize(char*filename){

	pt_net pt_n;
	FILE*fp;
	int i,j,node1,node2,num;
	double reliability;


	if((pt_n= (pt_net)calloc(1,sizeof(net)))==NULL){
		printf("\n Fail attempting to allocate memory...\n");
		exit(1);
	}



	if((fp= fopen(filename,"r"))==NULL){
		printf("\n Fail attempting to open a disk file...\n");
		exit(1);
	}
	fscanf(fp,"%d",&pt_n->s);
	fscanf(fp,"%d",&pt_n->t);
	fscanf(fp,"%d",&pt_n->numnodes);
	fscanf(fp,"%d",&pt_n->numlinks);


	if((pt_n->I= (pt_lnk*)calloc(pt_n->numnodes+1,sizeof(struct link*)))==NULL){
		printf("\n Fail attempting to allocate memory...\n");
		exit(1);
	}
	for(i= 0;i<=pt_n->numnodes;i++){
		if((pt_n->I[i]= (pt_lnk)calloc(pt_n->numnodes+1,sizeof(struct link)))==NULL){
			printf("\n Fail attempting to allocate memory...\n");
			exit(1);
		}
	}

	/* relation of min cut set to network size */
		//beta = breadth/pt_n->numnodes;
		beta = 0.5;

	for(i= 1;i<=pt_n->numnodes;i++)
		for(j= 1;j<=pt_n->numnodes;j++){
			pt_n->I[i][j].adj= 0;
			pt_n->I[i][j].rlb= 0.0;
			pt_n->I[i][j].smp= 0;
		}

		for(i= 1;i<=pt_n->numnodes;i++){
			fscanf(fp,"%d%d",&node1,&num);
			for(j= 1;j<=num;j++){
				fscanf(fp,"%d%lf",&node2,&reliability);
				pt_n->I[node1][node2].adj= 1;
				pt_n->I[node1][node2].rlb= reliability;
			}
		}

		fclose(fp);

		init_genrand(seed);

		if((list= (int**)calloc(pt_n->numnodes+1,sizeof(int*)))==NULL){
			printf("\n Fail attempting to allocate memory...\n");
			exit(1);
		}
		for(i= 0;i<=pt_n->numnodes;i++){
			if((list[i]= (int*)calloc(pt_n->numnodes,sizeof(int)))==NULL){
				printf("\n Fail attempting to allocate memory...\n");
				exit(1);
			}
		}


		for(i= 1;i<=pt_n->numnodes;i++)
			for(j= 1;j<=pt_n->numnodes;j++){
				list[i][j]= 0;
			}


			if((visited= (int*)calloc(pt_n->numnodes,sizeof(int)))==NULL){
				printf("\n Fail attempting to allocate memory...\n");
				exit(1);
			}

			for(i= 1;i<=pt_n->numnodes;i++)
				visited[i]= 0;

			connected= 0;




			return pt_n;
		}

/*:11*/
#line 337 "./crude1.w"
		;
/*12:*/
#line 456 "./crude1.w"

		int X(int node1,int node2,pt_net nt){
		/* returns 1 if link node1 –node2 is operative, 0 otherwise */
			double tempL;
			double epsilon;
			epsilon = 1-nt->I[node1][node2].rlb;
			if(U<beta){
								tempL = (1-epsilon)/(1-beta);

				l *= tempL;
				return 1;
			}else{
								tempL = epsilon/beta;
				l *= tempL;
				return 0;
			}
		}

				void Fail(pt_net nt){

					int i,j,k;

					for(i= 1;i<=nt->numnodes;i++)
						for(j= i+1;j<=nt->numnodes;j++)
							if(nt->I[i][j].adj){
								nt->I[i][j].smp= (nt->I[i][j].adj&&X(i,j,nt));
								nt->I[j][i].smp= nt->I[i][j].smp;
							}

							for(i= 1;i<=nt->numnodes;i++){
								k= 1;
								while(list[i][k]> 0){
									list[i][k++]= 0;
								}
							}

							for(i= 1;i<=nt->numnodes;i++){
								k= 1;
								for(j= 1;j<=nt->numnodes;j++){
									if(nt->I[i][j].smp){
										list[i][k++]= j;
									}
								}
							}
							return;
						}

/*:12*/
#line 338 "./crude1.w"
						;
/*13:*/
#line 502 "./crude1.w"

						int Phi(pt_net nt){

							int k;

							connected= 0;

							for(k= 1;k<=nt->numnodes;k++)
								visited[k]= 0;

							DFS(nt->s,nt);

							return connected;
						}

						void DFS(int node,pt_net nt){

							int k;

							if(connected)
								return;

							if(node==nt->t){
								connected= 1;
								return;
							}

							visited[node]= 1;

							k= 1;
							while(list[node][k]> 0){
								if(!visited[list[node][k]])
									DFS(list[node][k],nt);
								k++;
							}

							return;
						}

/*:13*/
#line 339 "./crude1.w"
						;

/*:10*/
#line 133 "./crude1.w"
						;

/*:1*/
