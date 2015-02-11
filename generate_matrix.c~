//Generate Input file with matrix of size specified in input

#include<stdio.h>
#include<stdlib.h>

int main(int argc,char *argv[]){
	FILE *file_pointer = fopen("mvm_data.txt","w");
	int m,n,i,j;

	printf("Give matrix size, Rows first, then column\n");
	scanf("%d %d",&m,&n);

	if(file_pointer==NULL){
		printf("Cannot open file");
	}
	fprintf(file_pointer,"%d %d\n",m,n);
	for(i=0;i<n-1;i++){
		fprintf(file_pointer,"%d ",rand()%100);
	}
	fprintf(file_pointer,"%d\n",rand()%100);
	for(i=0;i<m;i++){
		for(j=0;j<n-1;j++){
			fprintf(file_pointer,"%d ",rand()%100);
		}
		fprintf(file_pointer,"%d",rand()%100);
		fprintf(file_pointer,"\n");
	}
	printf("File created");
	return 0;
}
