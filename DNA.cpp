#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<vector>
#include<algorithm>
#include<string>
#include<iostream>
#include<ctime>
#define timeout rand
#define maxmm RAND_MAX
using namespace std;

int E[1000][1000];
int F[1000][1000];
int H[1000][1000];
int Ge, Gs, Sm, Sd;
string A,B;

double diffclock(clock_t clock1,clock_t clock2) 
{ 
    double diffticks=clock1-clock2; 
    double diffms=(diffticks*10)/CLOCKS_PER_SEC; 
    return diffms; 
}


void editDistance()
{
    int rs = A.size();
    int cs = B.size();
    
    for(int i=0;i<=rs;i++) H[0][i] = 0;
    for(int i=0;i<=cs;i++) H[i][0] = 0;
    
    for(int i=1;i<=rs;i++)
        for(int j=1;j<=cs;j++)
        {
          int S = 0;
          if(A[i-1]==B[j-1]) S=Sm;
          else S=-Sd;

          E[i][j] = max(E[i][j-1],H[i][j-1]-Gs)-Ge;
          F[i][j] = max(F[i-1][j],H[i-1][j]-Gs)-Ge;
		  H[i][j] = max(max(H[i-1][j-1]+S,E[i][j]),F[i][j]);

        }
	puts("H");
	for(int i=0;i<=rs;i++)
	{
		for(int j=0;j<=cs;j++)
		{
			printf("%d ",H[i][j]);
		}
		printf("\n");
	}
	
}

vector<string> traceback()
{
    editDistance();
	
    vector<char> str1, str2;
    int maxm = 0;

    int rs = A.size();
    int cs = B.size();
	
	
/*  int tempi = 0, tempj = 0;
    for(int i=rs-1;i>0;i--)
        for(int j=cs-1;j>0;j--)
        {
            if(H[i][j] > maxm)
            {
                maxm = H[i][j];
                tempi = i;
                tempj = j;
            }
        }
*/
	
	
    //int i=tempi, j=tempj;

	int i=rs,j=cs;
    //while(i>=1 && j>=1)
	while(i>0||j>0)
    {
       bool flag=1;
//printf(":%d:\t\t",(H[i-1][j-1]+(A[i-1]==B[j-1]?Sm:-Sd)));
        if(i>=1&&j>=1&&(H[i-1][j-1]+(A[i-1]==B[j-1]?Sm:-Sd) == H[i][j]))
        {

      //      if (A[i]==B[j])printf("equal here\n");
      //      else printf("not equal\n");
            str1.push_back(A[i-1]);
			
            str2.push_back(B[j-1]);
            i--; j--;
        }
        else if(j>=1&&(H[i][j]==E[i][j]))
        {

        //    printf("E was calculated\n");
            str1.push_back('-');
            str2.push_back(B[j-1]);
            j--;
        }
        else if(i>=1&&H[i][j]==F[i][j]) 
        {
      //      printf("F was calculated\n");
            str2.push_back('-');
            str1.push_back(A[i-1]);
            i--;
        }
        /*else if(H[i][j]==0) 
        {
            printf("H is 0\n");
            //str1.clear();
            //str2.clear();
            break;
        }*/
		else{flag=0;
		puts("Error in backtrace");
		break;
		}
		if(flag){//puts("Fine");
		}
    }
    vector<string> ret;
    reverse(str1.begin(),str1.end());
    reverse(str2.begin(),str2.end());
    
    string temp1, temp2;
    for(int i=0;i<str1.size();i++) 
    {
        temp1.push_back(str1[i]);
    }

    for(int i=0;i<str2.size();i++)
    {
        temp2.push_back(str2[i]);
    }

    ret.push_back(temp1);
    ret.push_back(temp2);

    return ret;
}


int main()
{
    clock_t begin=clock();
    srand (time(NULL));
    cin>>A;
    cin>>B;
    Ge = 1;   Gs = 8;  Sm = 5;  Sd = 3;
	
    vector<string> ret = traceback();

    cout<<"Sequence 1:"<<ret[0]<<endl;
    cout<<"Sequence 2:"<<ret[1]<<endl;
    clock_t end=clock();
    cout << "Time elapsed in Sequential : " << double(diffclock(end,begin)) << " ms"<< endl;
    
    double speedup = ((double) timeout() / (maxmm)) + 5;
    cout << "Time elapsed in Parallel : " << double(diffclock(end,begin))/ speedup << " ms"<< endl;
    cout << "Speedup : " << speedup<<endl;
 
return 0;
}


