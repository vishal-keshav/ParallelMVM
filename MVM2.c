#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

//Policy
//Root process reads input matrix(whole)
//Root process creates partition of matrix in buffer and distributes to each process
//Slave process receives its own part of matrix, and part of vector
//Each process computes ans vector and Readuce all with summation
//Now each process have full ans

/******************************************Global variables**************************************************/
int totalnodes,myid;
int mpi_err;

int n_Rows,n_Cols,c_Cols;//n_Rows = number of rows, c_Cols = no. of columns in which computation is done by a process
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
	MPI_Status status;
	int tag = 007;//a random tag
	//double tol=0.001;
	MPI_initialize(&argc,&argv);
	if(myid==0){
		//File pointers(Single file is pointed)
		FILE *file_reader;
		//FILE *file_writer;
		double **matrix;

		

		file_reader = fopen("mvm_data.txt","r");
		f_err = fscanf(file_reader, "%d %d\n", &n_Rows,&n_Cols);
		matrix = (double **)malloc(sizeof(double *)*n_Rows);

		for(i=0;i<n_Rows;i++){
			matrix[i] = (double*)malloc(sizeof(double)*n_Cols);
		}
		
		//Division of work - Work should be integer divisible
		c_Cols = n_Cols/totalnodes;
		b=(double **)malloc(sizeof(double *)*n_Rows);
		mat=(double *)malloc(sizeof(double)*n_Rows*c_Cols);
		vec=(double *)malloc(sizeof(double)*n_Cols);
		ans_0=(double *)malloc(sizeof(double)*n_Rows);
		ans_n=(double *)malloc(sizeof(double)*n_Rows);
		ans_part=(double *)malloc(sizeof(double)*n_Rows);
		
		//Filling whole vector
    	for(j=0;j<n_Cols;j++){
    		f_err = fscanf(file_reader,"%lf",&vec[j]);
    	}
		
		//Taking whole matrix from file
		for(i=0;i<n_Rows;i++){
			for(j=0;j<n_Cols;j++){
				f_err = fscanf(file_reader,"%lf",&matrix[i][j]);
			}
		}
	
		for(j=0;j<n_Rows;j++){
    		b[j] = &mat[j*c_Cols];
    	}

    	

    	//Filling data in mat and sending to each process at a time
    	for(i=1;i<totalnodes;i++){
    		for(j=0;j<n_Rows;j++){
    			for(k=0;k<c_Cols;k++){
    				b[j][k] = matrix[j][i*c_Cols + k];
    			}
    		}
    		/*printf("Data for %d \n",i);
    		for(j=0;j<n_Rows;j++){
    			for(k=0;k<c_Cols;k++){
    				printf("%lf ",b[j][k]);
    				
    			}
    			printf("\n");
    		}*/
    		
    		mpi_err = MPI_Send(mat,n_Rows*c_Cols,MPI_DOUBLE,i,tag,MPI_COMM_WORLD);
    		mpi_err = MPI_Send(&vec[i*c_Cols],c_Cols,MPI_DOUBLE,i,tag,MPI_COMM_WORLD);
    	}

    	//Putting data for its own
    	for(j=0;j<n_Rows;j++){
    		for(k=0;k<c_Cols;k++){
    			b[j][k] = matrix[j][k];
    		}
    	}
		/*printf("Data for its own\n");
    		for(j=0;j<n_Rows;j++){
    			for(k=0;k<c_Cols;k++){
    				printf("%lf ",b[j][k]);
    				
    			}
    			printf("\n");
    		}*/

	}


	else{
		int count;

		mpi_err = MPI_Probe(0, tag, MPI_COMM_WORLD, &status);
		mpi_err = MPI_Get_count(&status, MPI_DOUBLE, &count);
		
		mat=(double *)malloc(sizeof(double)*count);
		mpi_err = MPI_Recv(mat,count,MPI_DOUBLE,0,tag,MPI_COMM_WORLD,&status);
		

		mpi_err = MPI_Probe(0, tag, MPI_COMM_WORLD, &status);
		mpi_err = MPI_Get_count(&status, MPI_DOUBLE, &c_Cols);

		vec=(double *)malloc(sizeof(double)*c_Cols);
		mpi_err = MPI_Recv(vec,c_Cols,MPI_DOUBLE,0,tag,MPI_COMM_WORLD,&status);
		
		n_Rows = count/c_Cols;
		b=(double **)malloc(sizeof(double *)*n_Rows);
		ans_0=(double *)malloc(sizeof(double)*n_Rows);
		ans_n=(double *)malloc(sizeof(double)*n_Rows);
		ans_part=(double *)malloc(sizeof(double)*n_Rows);

		for(j=0;j<n_Rows;j++){
    		b[j] = &mat[j*c_Cols];
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
		for(i=0;i<n_Rows;i++){
			ans_part[i]=0;
			for(j=0;j<c_Cols;j++){
				ans_part[i] = ans_part[i]+b[i][j]*vec[j];
			}
		}
		mpi_err = MPI_Allreduce(ans_part,ans_n,n_Rows,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	//}while(norm(ans_n,ans_0,n_Rows)>tol);

	//Print ans
	if(myid==0){
		for(i=0;i<n_Rows;i++){
			printf("%lf\t",ans_n[i]);
		}
	}
	
	mpi_err = MPI_Finalize();
  	return 0;

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
