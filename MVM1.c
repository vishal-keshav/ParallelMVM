#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

//Policy
//Each process have right to read input file
//Each process reads its own part and do computation
//MPI Allgather is used to recollect computed value
//Each process checks for break condition
//Finally each have the final result

/******************************************Global variables**************************************************/
int totalnodes,myid;
int mpi_err;
int n_Rows,n_Cols,c_Rows;//n_Rows = number of rows, c_Rows = no. of rows in which computation is done by a process
double *mat,*vec,*ans_0,*ans_n;//1D structure of matrix
double **b;//pointer to different rows of the matrix
double *ans_part;//broadcast part

/*******************************Non-algorithmic function declaration****************************************/
void MPI_initialize(int *, char ***);
void file_seeker(FILE *fp,int n);
void show_data(double *data,int n);
double norm(double *v1,double *v2,int n);

/*******************************************Main function***************************************************/
int main(int argc,char *argv[]){
//Variables Declaration
	int i,k,j;

	int f_err;//File error msg
	//int otherid;//Id of the process for which root have to wait for broadcast data
	//MPI_Status status;
	//int tag = 007;//a random tag
	//double tol=0.001;

	//File pointers(Single file is pointed)
	FILE *file_reader;
	//FILE *file_writer;

	MPI_initialize(&argc,&argv);

	file_reader = fopen("mvm_data.txt","r");
	f_err = fscanf(file_reader, "%d %d\n", &n_Rows,&n_Cols);

	//Division of work - Work should be integer divisible
	c_Rows = n_Rows/totalnodes;
	b=(double **)malloc(sizeof(double *)*c_Rows);
	mat=(double *)malloc(sizeof(double)*c_Rows*n_Cols);
	vec=(double *)malloc(sizeof(double)*n_Cols);
	ans_0=(double *)malloc(sizeof(double)*n_Rows);
	ans_n=(double *)malloc(sizeof(double)*n_Rows);
	ans_part=(double *)malloc(sizeof(double)*c_Rows);
	
	for(j=0;j<c_Rows;j++){
    	b[j] = &mat[j*n_Cols];
    }
	//Scan vector
    for(i=0;i<n_Cols;i++){
    	f_err=fscanf(file_reader,"%lf",&vec[i]);
    }
    //Scan its part of matrix
    file_seeker(file_reader,c_Rows*n_Cols*myid);
    for(j=0;j<c_Rows;j++){
		 for(k=0;k<n_Cols;k++){
		 	f_err = fscanf(file_reader,"%lf",&b[j][k]);
		}
	}
	//Initialize guess ans
	for(i=0;i<n_Rows;i++){
		ans_n[i]=1;
	}
	//Computation - iteratively
	//do{
		for(i=0;i<n_Rows;i++){
			ans_0[i] = ans_n[i];
		}
		for(i=0;i<c_Rows;i++){
			ans_part[i]=0;
			for(j=0;j<n_Cols;j++){
				ans_part[i] = ans_part[i]+b[i][j]*vec[j];
			}
		}
		MPI_Allgather(ans_part,c_Rows,MPI_DOUBLE,ans_n,c_Rows,MPI_DOUBLE,MPI_COMM_WORLD);

	//}while(norm(ans_n,ans_0,n_Rows)>tol);

	//Print ans
	if(myid==0){
		for(i=0;i<n_Rows;i++){
			printf("%lf\t",ans_n[i]);
		}
	}
}

/*******************************************Main function ends***************************************************/

/********************Function defenition************************/
void MPI_initialize(int *argc,char ***argv){
	mpi_err = MPI_Init(argc,argv);
  mpi_err = MPI_Comm_size( MPI_COMM_WORLD,&totalnodes);
	mpi_err = MPI_Comm_rank(MPI_COMM_WORLD,&myid);
}

void file_seeker(FILE *fp,int n){
	int i,f_err;
	double temp;
	for(i=0;i<n;i++){
		f_err=fscanf(fp,"%lf",&temp);
	}
}

void show_data(double *data,int n){
	int i;
	for(i=0;i<n;i++){
		printf("%lf\t",data[i]);
	}
	printf("\n");
}

double norm(double *v1,double *v2,int n){
	double out=0;
	double temp;
	int i;
	for(i=0;i<n;i++){
		temp = abs(v1[i]-v2[i]);
		out = out + temp*temp;
	}
	out = sqrt(out);
	return out;
}
/*********************Function defenition ends*******************/