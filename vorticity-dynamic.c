#include <stdio.h>
#include  <stdlib.h> // important for calloc function
#include <float.h>
#include <math.h>
/* *************************************************************************************************************************************************************** *\
                                                      Dynamic memory creation Function
1. calloc function initialises the matrix value to 0.
2. malloc put random value in the array elements.

\* *************************************************************************************************************************************************************** */
//DynaMic MemoRy Allocation 
void freematrix(double **a,int rows)
{
    for(int i=0;i<rows;i++)
        free(a[i]);//removing columns
    free(a);//removing row
}
double ** allocate2d( int rows, int columns)
{
  double **create_matrix;
  create_matrix=(double **)calloc(rows,sizeof(double *));
  for(int i=0;i<rows;i++)
  {
    create_matrix[i]=(double *)calloc(columns,sizeof(double));
    if(!create_matrix[i])
    {
      freematrix(create_matrix,rows);
    }
  }

  return create_matrix;
}
/* *************************************************************************************************************************************************************** *\


\* *************************************************************************************************************************************************************** */
//FuNcTiOnS
double diff (double **arr, double **arr1, int m, int n)
{
	double max=0;
	int i, j;
	double temp=0;
	for(i=1;i<m-1;i++)
	{
		for(j=1;j<n-1;j++)
		{
			temp= (arr[i][j]-arr1[i][j])>0?(arr[i][j]-arr1[i][j]):(-1*(arr[i][j]-arr1[i][j]));
			if(temp>max)
			{
				max=temp;
			}
		}
	}
	return max;
}
/*
double diff (double *arr, double *arr1, int m, int n)
{
	double max=0;
	int i, j;
	double temp=0;
	for(i=1;i<m-1;i++)
	{
		for(j=1;j<n-1;j++)
		{
			temp= (*((arr+i*n) + j)-*((arr1+i*n)+j))>0?(*((arr+i*n) + j)-*((arr1+i*n)+j)):(-1*(*((arr+i*n) + j)-*((arr1+i*n)+j)));
			if(temp>max)
			{
				max=temp;
			}
		}
	}
	return max;
}
*/
float deltaT(float nu_m,float deltax,float deltay)
{
 float result;
 result = 0.5/(nu_m/(deltax*deltax) + nu_m/(deltay*deltay));
 return result;
}
int main()
{
//VaRiaBlEs
 int i,j;
 // geometry length
 float length=5;
 float height=1;
 //number of grid points
 const int nx=50;
 const int ny=20;
 //spacing
 float dx=length/(nx-1);
 float dy=height/(ny-1);
 float nu=0.01; //kinematic viscosity
 float dt=0.001;//deltaT(nu,dx,dy);//time step
 float U=1;//left wall velocity inlet
 float q=0.5;//discretisation scheme
 float tol=0.0001;
 //initializing matrix
 double **psi;
 double  **omega, **u, **v;
 double **psi_prev,**omega_prev,**u_prev,**v_prev;
 double **psi_iter;
 psi=allocate2d(nx,ny);
 psi_prev=allocate2d(nx,ny);
 psi_iter=allocate2d(nx,ny);
 omega_prev=allocate2d(nx,ny);
 omega=allocate2d(nx,ny);
 u=allocate2d(nx,ny);
 u_prev=allocate2d(nx,ny);
 v=allocate2d(nx,ny);
 v_prev=allocate2d(nx,ny);

 /*
 %
 double psi[nx][ny],omega[nx][ny],u[nx][ny],v[nx][ny];
 double psi_prev[nx][ny],omega_prev[nx][ny],u_prev[nx][ny],v_prev[nx][ny];
 double psi_iter[nx][ny];
 %
 */
 //inittialisation of velocity with Boundary conditions initial
 for(i=0;i<nx;i++)
   {
     for(j=0;j<ny;j++)
       {
	 if(j==0)
	   {
	     u[i][j]=0;
	     v[i][j]=0;
	     psi[i][j]=0;
	     omega[i][j]=0;//initialise
	   }
	 else  if (j==ny-1)
	   {
	     u[i][j]=0;
	     v[i][j]=0;
	     psi[i][j]=0+U*height;
	     omega[i][j]=0;
	   }
	 else if (i==0)
	   {
	     u[i][j]=U;
	     v[i][j]=0;
	     psi[i][j]=0+U*dy*j;
	     omega[i][j]=0;
	   }
	 else if (i==nx-1)
	   {
	     u[i][j]=U;
	     v[i][j]=0;
		 psi[i][j]=0;
	     omega[i][j]=0;
	   }
	 else
	   {
	     u[i][j]=0;
	     v[i][j]=0;
	     psi[i][j]=0;//check?
	     omega[i][j]=0;
	   }
       }
   }
      for(i=0;i<nx;i++)
	{
		omega[i][ny-1]=(-7*psi[i][ny-1]+8*psi[i][ny-2]-psi[i][ny-3])/(2*dy*dy);//top wall
		omega[i][0]=(7*psi[i][0]-8*psi[i][1]+psi[i][2])/(2*dy*dy);//bottom wall
	}

     for(i=0;i<nx;i++)
   {
     for (j=0;j<ny;j++)
       {
	 psi_prev[i][j]=psi[i][j];
	 // omega_prev[i][j]=omega[i][j];
	 u_prev[i][j]=u[i][j];
	 v_prev[i][j]=v[i][j];
	 omega_prev[i][j]=omega[i][j];
     psi_iter[i][j]=psi[i][j];
       }
   }
float del=1;
float psi_error=1;
float t=1;
while(del>tol)
{
        for(j=1;j<ny-1;j++)
    {
        omega[0][j]=2*(psi[0][j]-psi[1][j])/(dx*dx);//inlet BC
    }
    for(i=1;i<nx;i++)
    {
        for(j=0;j<ny;j++)
        {
            if (i==(nx-1))
            {
                omega[i][j]=omega[i-1][j];
            }
            else if (j==0)
            {
                omega[i][j]=(7*psi[i][0]-8*psi[i][1]+psi[i][2])/(2*dy*dy);
            }
            else if (j==(ny-1))
            {
                omega[i][j]=(-1)*(-7*psi[i][ny-1]+8*psi[i][ny-2]-psi[i][ny-3])/(2*dy*dy);
            }
            else{
                omega[i][j]=omega_prev[i][j]+dt*nu*((omega_prev[i+1][j]-2*omega_prev[i][j]+omega_prev[i-1][j])/(dx*dx) + (omega_prev[i][j+1]-2*omega_prev[i][j]+omega_prev[i][j-1])/(dy*dy))-dt*u_prev[i][j]*(omega_prev[i+1][j]-omega_prev[i-1][j])/(2*dx) -dt*v_prev[i][j]*(omega_prev[i][j+1]-omega_prev[i][j-1])/(2*dy);
            }
        }
    }
   //added
    for (i=1;i<nx-1;i++)
   {
	for(j=1;j<ny-1;j++)
	{
	   if (i==1 || i==nx-2)
	   {
		   omega[i][j]=omega[i][j];
	   }
	   else{
		   omega[i][j]=omega[i][j] - dt*q*((u_prev[i][j]>0?0:u_prev[i][j])*(omega_prev[i-1][j]-3*omega_prev[i][j]+3*omega_prev[i+1][j]-omega_prev[i+2][j])/(3*dx))-dt*q*((u_prev[i][j]>0?u_prev[i][j]:0)*(omega_prev[i-2][j]-3*omega_prev[i-1][j]+3*omega_prev[i][j]-omega_prev[i+1][j])/(3*dx));
	   }
    }
   }
   //term 2
   for (i=1;i<nx-1;i++)
   {
	for(j=1;j<ny-1;j++)
	{
	   if (j==1 || j==ny-2)
	   {
		   omega[i][j]=omega[i][j];
	   }
	   else{
		   omega[i][j]=omega[i][j] - dt*q*((v_prev[i][j]>0?0:v_prev[i][j])*(omega_prev[i][j-1]-3*omega_prev[i][j]+3*omega_prev[i][j+1]-omega_prev[i][j+2])/(3*dy))-dt*q*((v_prev[i][j]>0?v_prev[i][j]:0)*(omega_prev[i][j-2]-3*omega_prev[i][j-1]+3*omega_prev[i][j]-omega_prev[i][j+1])/(3*dy));
	   }
    }
   }
   //last


    psi_error=1;
    while(psi_error>tol)
    {
        for(i=1;i<nx;i++)
        {
            for(j=0;j<ny;j++)
            {
                if (i==(nx-1))
                {
                    psi[i][j]=psi[i-1][j];
                }
                else if (j==0)
                {
                    psi[i][j]=0;
                }
                else if (j==(ny-1))
                {
                    psi[i][j]=U*height;
                }
                else{
                    psi[i][j]=0.5*(dx*dx*dy*dy/(dx*dx+dy*dy))*(omega[i][j] + (psi[i+1][j]+psi[i-1][j])/(dx*dx) + (psi[i][j+1]+psi[i][j-1])/(dy*dy));
                }
            }
        }
        //added
        for(j=1;j<ny-1;j++)
        {
          psi[0][j]=(4*psi[1][j]-psi[2][j])/3;
        }
        psi_error=diff(psi,psi_iter,nx,ny);
        for(i=0;i<nx;i++)
		{
			for (j=0;j<ny;j++)
			{
			psi_iter[i][j]=psi[i][j];
			}
		}
    }
    for(i=1;i<nx;i++)
    {
        for(j=0;j<ny;j++)
        {
            if(i==(nx-1))
            {
                u[i][j]=u[i-1][j];
            }
            else if (j==(ny-1))
            {
                u[i][j]=0;
            }
            else if (j==0)
            {
                u[i][j]=0;
            }
            else{
                u[i][j]=(psi[i][j+1]-psi[i][j-1])/(2*dy);
            }
        }
    }
    for(i=1;i<nx;i++)
    {
        for(j=0;j<ny;j++)
        {
            if(i==(nx-1))
            {
                v[i][j]=v[i-1][j];
            }
            else if (j==(ny-1))
            {
                v[i][j]=0;
            }
            else if (j==0)
            {
                v[i][j]=0;
            }
            else{
                v[i][j]=(-1)*(psi[i+1][j]-psi[i-1][j])/(2*dx);
            }
        }
    }
    del=diff(omega,omega_prev,nx,ny);
     for(i=0;i<nx;i++)
   {
     for (j=0;j<ny;j++)
       {
	 psi_prev[i][j]=psi[i][j];
	 // omega_prev[i][j]=omega[i][j];
	 u_prev[i][j]=u[i][j];
	 v_prev[i][j]=v[i][j];
	 omega_prev[i][j]=omega[i][j];
     psi_iter[i][j]=psi[i][j];
       }
   }
t=t+1.0;
//printf("%f",t);
}

//check

 FILE *fp;
 fp=fopen("test.txt","w");

  for(i=0;i<nx;i++)
  {
      for(j=0;j<ny;j++)
      {
          fprintf(fp,"%f \t",u[i][j]);
          //printf("%f \t",omega[i][j]);
      }
      //printf("\n");
      fprintf(fp,"\n");
  }
  fclose(fp);

  FILE *fp1;
  fp1=fopen("Plot.dat","w");
  for(j=0;j<ny;j++)
    {
        fprintf(fp1,"%f \t %f",u[nx-1][j],j*dy);
        fprintf(fp1,"\n");
        //printf("%f \t",omega[i][j]);
    }
  fclose(fp1);

  FILE *fp2;
  fp2=fopen("Str.txt","w");
  for(i=0;i<nx;i++)
  {
      for(j=0;j<ny;j++)
      {
          fprintf(fp2,"%f \t %f \t %f \t",i*dx,dy*j,psi[i][j]);
          //printf("%f \t",omega[i][j]);
          fprintf(fp2,"\n");
      }
      //printf("\n");
      
  }
  fclose(fp2);

return 0;
}
