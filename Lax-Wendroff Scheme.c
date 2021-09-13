 
 #include<stdio.h>
 #include<stdlib.h>
 #include<conio.h>
 #include<math.h>
 #define nx 200
 #define m  200
 #define L 1
 
 int main()
 {
 	int count = 0,i;
    double dx,x[m],u_old[m],u_new[m],lamda = 0.8;
 	FILE *g,*p;
 	g = fopen("final(n=40).dat","w");
 	p= fopen("initial.dat(n=40)","w");
 	//****Grid generation****//
 	dx = L/(nx-1);
 	for(i=1;i<=nx;i++)
 	{
 		x[i] = (i-1)*dx;	
	}
	//**initialize**//
	for(i=1;i<=(nx+1)/2;i++)
	{
		u_old[i] = 0.0;
	}
	for(i=(nx+1)/2+1;i<=nx;i++)
	{
		u_old[i] = 1.0;
	}
	//***writing the initial solution to corresponding output file***//
	for(i=1;i<=nx;i++)
	{
		fprintf(g,"\n%lf\t%lf",x[i],u_old[i]);
	}
	//***time marching***//
	while(count<=40)
	{
		for(i=2;i<=nx-1;i++)
		{
			u_new[i] = u_old[i] - 0.5*lamda*(u_old[i+1] - u_old[i-1]) + 0.5*lamda*lamda*(u_old[i+1] - 2*u_old[i] + u_old[i-1]);
		}
		for(i=2;i<=nx-1;i++)
		{
			u_old[i] = u_new[i];
		}
		printf("\n%d",count);
	count++;
	}
	//***writng the final solution to the***//
	for(i=1;i<=nx;i++)
	{
		fprintf(p,"\n%lf\t%lf",x[i],u_old[i]);
	}
	fclose(g);
	fclose(p);
	return (0);
 }
