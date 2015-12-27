#include<stdio.h>
#include<math.h>
#include<cuda.h>
#define threads 8


/* 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~README~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
1. To run the code : 
    
   nvcc -deviceemu cuda1.cu  (in emulation mode)
   a.exe 8                   (exefile sizeofstring)


2. Optimizations not performed 
3. Configuration for each kernel function will change according to size
4. Code not tested for various other sizes
5. Parallel max algo may be slower than expected

*/



/* 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~CODE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
kernel functions
1 : creation of substitution matrix (H)    ( Time complexity : O( Nlog(N) )  )
2 : Finding maximum in H matrix            ( Time complexity : O( Nlog(N) )  )
3 : creation of bactrace matrix (pr_dest)  (stores destination location from each cell for performing bactracking )   ( Time complexity : O( N )  )

sequential code : Finds aligned sequences (backtracking)    ( Time complexity : O( N )  )
*/



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  1 : creation of H matrix  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~// 

__global__ void matrix_num(char *X,char *Y,int* FF, int* HH1, int* HH ,int* EE,int N , int i, int Ge, int Gs, int S, int D)
{    
    //find Row and Column corresponding to a data element for each thread
    int loc = blockIdx.y * blockDim.y + threadIdx.y;

	int inloc = i * (N+1) + loc + 1;		//value to be updated
	

	///// done with FF ////
	loc = (i-1)*(N+1)+loc+1;	
	
	FF[inloc] = FF[loc] - Ge ;
	if(HH[loc] - Gs - Ge  > FF[inloc])
		FF[inloc] = HH[loc] - Gs - Ge ;
		
	///// done with HH1 ///
	if( X[i-1]==Y[blockIdx.y * blockDim.y + threadIdx.y ] )
		HH1[ inloc ] = HH[ loc-1 ] + S ; 	// loc-1 is diagonal above
	else 
		HH1[ inloc ] = HH[ loc-1 ] + D ;  
	
	if( HH1[ inloc ] < 0 ) HH1[ inloc ] = 0 ;
	if( HH1[ inloc ] < FF[inloc] )  HH1[ inloc ] = FF[inloc] ;	
	EE[inloc]  = HH1[inloc] ;
}

__global__ void matrix_num3(int* HH1, int* EE, int* HH, int N , int i, int Gs)
{
	int inloc = i*(N+1) + blockIdx.y * blockDim.y + threadIdx.y + 1;
	HH[inloc] = HH1[inloc]  ;
	if( HH[inloc] < EE[inloc] - Gs )
		HH[inloc] = EE[inloc] - Gs ; 	
}

__global__ void prefixsum(int* HH1, int* EE , int index ,int N, int Ge){
	
	/// loc+1 is the memory

	int k,loc = index*(N+1)+blockIdx.y * blockDim.y + threadIdx.y, j ;
	int temp ;

	j=N;

	/// upsweep
	for(j = 2; j < N; j*=2){
		if( (loc % j) == 0 && (loc+1+j-1)/(N+1)==(loc)/(N+1) )
			EE[loc+1+j-1] = (EE[loc+1+j/2-1] > EE[loc+1+j-1] + Ge * j/2) ? EE[loc+1+j/2-1] : EE[loc+1+j-1] + Ge * j/2 ;
	}

	if( (loc+1)%(N+1) == N ){
    	EE[loc+1] = 0 ;	
    }

}

__global__ void prefixsum2(int* HH1, int* EE , int index ,int N, int Ge){

	int k,loc = index*(N+1)+blockIdx.y * blockDim.y + threadIdx.y, j = N;
	int temp ;    
	

    /// downsweep
	for(j = j ; j > 1 ; j/=2 ) {
		if((loc % j) == 0)  {   //// intially the first thread works 

			temp = EE[loc + 1 + j/2 -1];
			EE[loc + 1 + j/2 -1] = EE[loc + 1 + j - 1];
			EE[loc + 1 + j - 1] = ((temp > EE[loc + 1 + j - 1]) ? temp  : EE[loc + 1 + j - 1]) -	  Ge * j/2;
		} 
	}

}



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 2 : parallel max ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

__global__ void findStore(int* HH , int* store , int N , int i)
{
	int locx = blockIdx.x * blockDim.x + threadIdx.x+1 ;
	int locy = blockIdx.y * blockDim.y + threadIdx.y+1 ;

	store[ locy*(N+1) + locx ] = 0 ;

	if(  HH[ i*(N+1)+locy ] >=  HH[ i*(N+1)+locx ]  ) store[locy*(N+1)+locx] = 1 ;
		
}

__global__ void initAnd(int* array , int N)
{
	array[ blockIdx.y * blockDim.y + threadIdx.y+1] = 1 ;
}

__global__ void parAnding(int* store , int* array , int N)
{
	int locx = blockIdx.x * blockDim.x + threadIdx.x+1 ;
        int locy = blockIdx.y * blockDim.y + threadIdx.y+1 ;

	array[locy] = array[locy] & store[ locy*(N+1)+locx ] ; 
}

__global__ void giveMax(int* HH,int* array,int N, int i,int *maxX, int *maxY)
{
	int locy = blockIdx.y * blockDim.y + threadIdx.y+1 ;
	if( array[locy]==1 && HH[i*(N+1)+locy]  > HH[*maxY * (N+1) + *maxX] ) { *maxX = locy ; *maxY = i ; }  
}



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 3 : back trace ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//


__global__ void traceback(char*XX , char* YY , int* HH , int* dest , int N , int S, int D, int Gs, int Ge)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
	int j = blockIdx.y * blockDim.y + threadIdx.y + 1;

	int k ;
	
	int loc=i*(N+1)+j;
	dest[loc]=0 ;

	if(i>0 && j>0){

		if( XX[i-1]==YY[j-1] && HH[(i-1)*(N+1)+j-1]+S == HH[loc] ) dest[loc] = (i-1)*(N+1)+j-1 ;
	
		else if( XX[i-1]!=YY[j-1] && HH[(i-1)*(N+1)+j-1]+D == HH[loc] ) dest[loc] = (i-1)*(N+1)+j-1 ;

		else {

			for( k=i-1; k>=0 ; k-- ) {
				if ( HH[ k*(N+1)+j ] - Gs - (i-k)*Ge   == HH[loc] ) { dest[loc] = k*(N+1)+j ; break ; }
			}

			for( k=j-1; k>=0 ; k--) {
				if ( HH[ i*(N+1)+k ] - Gs - (j-k)*Ge == HH[loc] ) { dest[loc] = i*(N+1)+k ; break ; }
			} 
		}

 	}

}





//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  MAIN FUNCTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//



int main(int argc, char *argv[])
{
	int N, i, j;				//N == size of square matrix
	int Ge=1, Gs=8, S=5, D=-3;
	N = atoi(argv[1]);	
	
	int *F,*FF,*H1,*HH1,*E,*EE,*H,*HH,*dest,*store,*array;
	char *X,*Y,*a,*b,*seq1,*seq2;

	int *pr_dest ;

	size_t size=sizeof(int)* (N+1) * (N+1);
	size_t size_str= sizeof(char)* N ;
	size_t sizeN=sizeof(int)* (N+1) ;

	freopen("oupt.txt","w",stdout);
	
    //allocate host side memory
	a=(char*)malloc(size_str);
	b=(char*)malloc(size_str);
	seq1=(char*)malloc(size_str);
	seq2=(char*)malloc(size_str);
	F=(int*)malloc(size);
	H1=(int*)malloc(size);
	E=(int*)malloc(size);
	H=(int*)malloc(size);

	//dest=(int*)malloc(size);
	pr_dest=(int*)malloc(size);
	//store=(int*)malloc(size);
	//array=(int*)malloc(size);


	int *maxXX = 0, *maxYY=0, *maxX = 0, *maxY = 0 ;
	maxX = (int*) calloc(1, sizeof(int));
	maxY = (int*) calloc(1, sizeof(int));
	cudaMalloc(&maxXX, sizeof(int));
	cudaMalloc(&maxYY, sizeof(int));
	cudaMemcpy(maxXX, maxX, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(maxYY, maxY, sizeof(int), cudaMemcpyHostToDevice);

	FILE *fp = fopen("seq.txt","r") ;
	
	fgets(a,20000,fp) ;
	fgets(b,20000,fp) ;


    //allocate device memory
	cudaMalloc(&X,size_str);  //	printf("\nAfter cudaMalloc for X\t%s\n",cudaGetErrorString(cudaGetLastError()));
	cudaMalloc(&Y,size_str);  //	printf("\nAfter cudaMalloc for Y\t%s\n",cudaGetErrorString(cudaGetLastError()));
	cudaMalloc(&FF,size);     //	printf("\nAfter cudaMalloc for FF\t%s\n",cudaGetErrorString(cudaGetLastError()));
	cudaMalloc(&HH1,size);    //	printf("\nAfter cudaMalloc for HH1\t%s\n",cudaGetErrorString(cudaGetLastError()));
	cudaMalloc(&EE,size);     //	printf("\nAfter cudaMalloc for EE\t%s\n",cudaGetErrorString(cudaGetLastError()));
	cudaMalloc(&HH,size);     //	printf("\nAfter cudaMalloc for HH\t%s\n",cudaGetErrorString(cudaGetLastError()));
    cudaMalloc(&dest,size);   //    printf("\nAfter cudaMalloc for HH\t%s\n",cudaGetErrorString(cudaGetLastError()));
	cudaMalloc(&store,size);  //    printf("\nAfter cudaMalloc for HH\t%s\n",cudaGetErrorString(cudaGetLastError()));
	cudaMalloc(&array,sizeN); //    printf("\nAfter cudaMalloc for HH\t%s\n",cudaGetErrorString(cudaGetLastError()));


	
	for(i=0;i<=N;i++){   for(j=0;j<=N;j++){  F[i*(N+1)+j]=0;  H1[i*(N+1)+j]=0;  E[i*(N+1)+j]=0;  H[i*(N+1)+j]=0;  }  }
		
    cudaMemcpy(X,a,size_str,cudaMemcpyHostToDevice);
	cudaMemcpy(Y,b,size_str,cudaMemcpyHostToDevice);
 	cudaMemcpy(FF,F,size,cudaMemcpyHostToDevice);
 	cudaMemcpy(HH1,H1,size,cudaMemcpyHostToDevice);
 	cudaMemcpy(EE,E,size,cudaMemcpyHostToDevice);
 	cudaMemcpy(HH,H,size,cudaMemcpyHostToDevice);

//	printf("\nAfter HostToDevice Memcpy\n%s\n",cudaGetErrorString(cudaGetLastError()));








    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ calculate execution configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	    	
	dim3 blocksize(1,threads);		        //each block contains 16 * 16 (=256) threads 
	int k = (N/threads) + (N%threads != 0);
	dim3 gridsize(1,k);			//creating just sufficient no of blocks
    

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ kernel call 1 : H-matrix calculation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~// 
    //GPU timer code
   	
   	float time1;
   	cudaEvent_t start1,stop1;			
	cudaEventCreate(&start1);		
	cudaEventCreate(&stop1);
	cudaEventRecord(start1,0);

	/////for(i=0;i<=N;i++){for(j=0;j<=N;j++){ H1[i*(N+1)+j]=5; }} cudaMemcpy(HH1,H1,size,cudaMemcpyHostToDevice);
	
	for(i=1;i<=N;i++)
	{
		matrix_num <<< gridsize, blocksize >>> (X, Y, FF, HH1,HH,EE, N, i, Ge, Gs, S, D);
		prefixsum <<< gridsize, blocksize >>> (HH1, EE ,i, N, Ge); 
		prefixsum2 <<< gridsize, blocksize >>> (HH1, EE ,i, N, Ge); 	
		matrix_num3 <<< gridsize, blocksize >>>  (HH1 , EE , HH ,N ,i ,Gs );	
	}

	cudaEventRecord(stop1,0);
	cudaEventSynchronize(stop1);
	cudaEventElapsedTime(&time1,start1,stop1);			//time taken in kernel call calculated
	cudaEventDestroy(start1);
	cudaEventDestroy(stop1);




    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ calculate execution configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
   
    dim3 block(threads,threads);                      //each block contains 16 * 16 (=256) threads
    k = (N/threads) + (N%threads != 0);
    dim3 grid(k,k);        


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ kernel call 2 : Parallel Max calculation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~// 
    //GPU timer code
   	
   	float time2;
   	cudaEvent_t start2,stop2;			
	cudaEventCreate(&start2);		
	cudaEventCreate(&stop2);
	cudaEventRecord(start2,0);

	for(i=1;i<=N;i++)
	{
		findStore  <<< grid, block >>> (HH, store, N, i);
		initAnd <<< gridsize, blocksize >>> (array , N );
		parAnding <<< grid, block >>> (store , array , N) ;
		giveMax   <<< gridsize, blocksize >>> (HH , array , N ,i , maxXX , maxYY) ; 	

	}

	cudaEventRecord(stop2,0);
	cudaEventSynchronize(stop2);
	cudaEventElapsedTime(&time2,start2,stop2);			//time taken in kernel call calculated
	cudaEventDestroy(start2);
	cudaEventDestroy(stop2);

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ kernel call 3 : Traceback Matrix calculation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~// 
    //GPU timer code

    float time3;
   	cudaEvent_t start3,stop3;			
	cudaEventCreate(&start3);		
	cudaEventCreate(&stop3);
	cudaEventRecord(start3,0);

    traceback <<< grid , block >>> (X,Y,HH,dest,N,S,D,Gs,Ge);

    cudaEventRecord(stop3,0);
	cudaEventSynchronize(stop3);
	cudaEventElapsedTime(&time3,start3,stop3);			//time taken in kernel call calculated
	cudaEventDestroy(start3);
	cudaEventDestroy(stop3);

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//



	cudaMemcpy(E,EE,size,cudaMemcpyDeviceToHost);
	cudaMemcpy(H,HH,size,cudaMemcpyDeviceToHost);		
	cudaMemcpy(pr_dest,dest,size,cudaMemcpyDeviceToHost);	
	cudaMemcpy(maxX, maxXX, sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(maxY, maxYY, sizeof(int), cudaMemcpyDeviceToHost);






    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Backtrace ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

	int i1,j1,t,l,count;
    i=*maxY; j=*maxX; t=0; l=0;

    while(i>0 && j>0) {
        i1=pr_dest[i*(N+1)+j]/(N+1);
        j1=pr_dest[i*(N+1)+j]%(N+1);
 		
        if((i-i1) > (j-j1) ) {
            count = j1-i1;
            while(count-- > 0){
                seq1[l++] = a[i1];
                seq2[t++] = '-';
            }
        }
	
		else if(i-i1 < j-j1){
            count = i1-j1;
            while(count-- > 0){
                seq2[t++] = b[j1];
                seq1[l++] = '-';
            }
        }
	
	    else{			
            seq1[l++]=a[i1];
            seq2[t++]=b[j1];
        }

        i= i1;
        j= j1;
    }
	

//	while(j>=0){   seq2[t++]=b[j--];   seq1[l++]='-';   }
//  while(i>=0){   seq1[l++]=a[i--];   seq2[t++]='-';   }





//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//


	printf("Substitution Matrix \n\n");

	for(i=0;i<=N;i++)
	{
		for(j=0;j<=N;j++){
			printf("%d\t",H[i*(N+1)+j]);
		}
		printf("\n");
	}


    printf("\n\nBacktrack Path Matrix \n\n");
	for(i=0;i<=N;i++){
        for(j=0;j<=N;j++){
            printf("%d\t",H[i*(N+1)+j]);
            printf("(%d,%d)\t",     pr_dest[i*(N+1)+j]/(N+1) , pr_dest[i*(N+1)+j]%(N+1) );
        }
        printf("\n");
    }



    printf("\n\nOriginal Sequences\n\n");
	for(j=0;j<N;j++){
        printf("%c",a[j]);
    }
    printf("\n");

    for(j=0;j<N;j++){
        printf("%c",b[j]);
    }
    printf("\n");
  
    printf("\n\nMaximally locally Alligned Sequence\n\n");
	for(j=l-1;j>=0;j--){
        printf("%c",seq1[j]);
    }
    printf("\n");

    for(j=t-1;j>=0;j--){
        printf("%c",seq2[j]);
    }
    printf("\n");

    printf("\n\nMaximum\n");
    printf("\n%d %d H[x][y]=%d\n",*maxY,*maxX,H[*maxY*(N+1)+*maxX]);


    printf("\nTime taken :\n");
    printf("Time to find substitution matrix = %f (ms)\n",time1);
    printf("Time to find maximum = %f (ms)\n",time2);
    printf("Time to find traceback = %f (ms)\n",time3);


    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DEBUGGING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PREFIX SUM DEBUGGING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//


	
	//printf("\n");
	//H1[0]=0 ;  H1[1]=2 ; H1[2]=3 ; H1[3]=4 ; H1[4]=5 ; H1[5]=6 ; H1[6]=7 ; H1[7]=8 ; H1[8]=9 ; cudaMemcpy(HH1,H1,size,cudaMemcpyHostToDevice);
	//E[0]=0 ;  E[1]=2 ; E[2]=3 ; E[3]=4 ; E[4]=5 ; E[5]=6 ; E[6]=7 ; E[7]=8 ; E[8]=9 ; cudaMemcpy(EE,E,size,cudaMemcpyHostToDevice);
	//E[0]=0 ;  E[1]=2 ; E[2]=4 ; E[3]=4 ; E[4]=8 ; E[5]=6 ; E[6]=8 ; E[7]=8 ; E[8]=0 ; cudaMemcpy(EE,E,size,cudaMemcpyHostToDevice);
	//for(i=0;i<=N;i++) { printf( "%d\t",E[i] ) ; } printf("\n");

	//prefixsum <<< gridsize, blocksize >>> (HH1, EE , 0 , N, Ge); 
	//prefixsum2 <<< gridsize, blocksize >>> (HH1, EE , 0, N, Ge); 	
   	//cudaMemcpy(E,EE,size,cudaMemcpyDeviceToHost);
   	//for(i=0;i<=N;i++) { printf( "%d\t",E[i] ) ; }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//	

    


    cudaFree(FF); cudaFree(HH1); cudaFree(EE) ; cudaFree(HH); cudaFree(X); cudaFree(Y); cudaFree(array); cudaFree(store); cudaFree(dest);	



    return 0;
}


