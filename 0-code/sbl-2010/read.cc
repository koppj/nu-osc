#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#define NMAX 10000
#define NFILES 20
#define LREAD 150

int read_numb_of_files;

void errorexit(char error_text[])
{
   fprintf(stderr, "%s\n",error_text);
   exit(1);
   return;
}

void resetread(void)
{
  read_numb_of_files=0;
}

double read(char *name, double x0)
{
   static int i[NFILES],f=0;
   static float y[NFILES][NMAX], x[NFILES][NMAX];
   static char files[NFILES][LREAD];
   int r,j;
   FILE *fp;

   if(f==0){
     read_numb_of_files=0;
     f=1;
   }
   if(read_numb_of_files==0){
      for(j=1;j<NFILES; j++) strcpy(files[j],"");
   }
   j=0;
   while(j<read_numb_of_files && strcmp(files[j],name)!=0) {
      j++;
   }
   if(j==read_numb_of_files){
      read_numb_of_files++;
      if(read_numb_of_files==NFILES)errorexit("read: too many files");
      strcpy(files[j],name);
      fp=fopen(name,"r");
      if(fp==NULL)errorexit("read: cannot open file");
      i[j]=-1;
      do{
	 i[j]++;
	 r=fscanf(fp,"%e  %e\n",&x[j][i[j]],&y[j][i[j]]);
      }while(r==2 && i[j]<NMAX-1);
      if(i[j]==NMAX-1)errorexit("read: file to long");
      if(i[j]==0)errorexit("read: cannot read file ");
      i[j]--;
      fclose(fp);
   }
   if(x0>x[j][i[j]]) r=i[j]-1;
   else{ 
      r=0;
      while(x0>x[j][r+1])r++;
   }
   return y[j][r]+(x0-x[j][r])*(y[j][r]-y[j][r+1])/(x[j][r]-x[j][r+1]);
}

#undef NMAX 
#undef NFILES

