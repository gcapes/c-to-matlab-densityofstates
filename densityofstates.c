#include <stdio.h>
#include <math.h>
#include <malloc.h>
int n, m; /* Nanotube parameters */
int dr;
int Nhex; /* Number of hexagons in the 1D unit cell */
double PI=3.141592653;
double circum; /* Circumference of the nanotube in unit of a0 */
double period; /* Period of the nanotube along its axis */
double Emax=9.0; /* Energy range to take into account */
double E0=2.9; /* Hopping parameter without strain */
double a0=2.49; /* Length of the lattice vector without strain */
double E1, E2, E3; /* New hopping parameters with strain */
double r1, r2, r3; /* Strained bond lengths */
float straint; /* Strain along the nanotube axis */
float strainc; /* Strain along the nanotube circumference */
int gcdn(int n, int m); /* Program to calculate the gcdn of n and m */
void density(); /* Subprogram to calculate the electronic DOS */

void main(){
	printf("Enter index n: ");
	scanf("%d", &n);
	printf("Enter index m: ");
	scanf("%d", &m);
	printf("Enter strain along axis in percent: ");
	scanf("%f", &straint);
	printf("Enter strain along circumference in percent: ");
	scanf("%f", &strainc);
	printf("\n");
	
	circum=sqrt(n*n+m*m+n*m);
	dr=gcdn(2*n+m,2*m+n);
	Nhex=2*(n*n+m*m+n*m)/dr;
	period=(1+straint/100)*sqrt(3)*circum*a0/dr;
	r1=sqrt((n+m)*(n+m)*a0*a0*(1+strainc/100)*(1+strainc/100)/4/
			circum/circum+(n-m)*(n-m)*a0*a0*(1+straint/100)*
			(1+straint/100)/12/circum/circum);
	r2=sqrt(m*m*a0*a0*(1+strainc/100)*(1+strainc/100)/4/circum/circum+
			(2*n+m)*(2*n+m)*a0*a0*(1+straint/100)*(1+straint/100)
			/12/circum/circum);
	r3=sqrt(n*n*a0*a0*(1+strainc/100)*(1+strainc/100)/4/circum/circum+(n+2
			*m)*(n+2*m)*a0*a0*(1+straint/100)*(1+straint/100)/12/
			circum/circum);
	E1=E0*a0*a0/3/r1/r1;
	E2=E0*a0*a0/3/r2/r2;
	E3=E0*a0*a0/3/r3/r3;
	density();
}

void density(){
	int i, j, p, iter, channel;
	double interp, kt, energy0, energy1, result, slope;
	double med1, med2, med3, med4;
	double *data;
	FILE *fp2;
	fp2=fopen("density.txt","w");
	
	printf("Enter number of points in the graph N: ");
	scanf("%d", &iter);
	printf("Enter number of points analyzed p: ");
	scanf("%d", &p);
	
	data=(double *)(malloc(sizeof(double)*iter));
	interp=PI/period/(double)p;
	med1=PI/circum/circum;
	med2=sqrt(3)/2/circum*a0;
	
	for(i=0; i<iter; i++){
		data[i]=0;
	}
	
	for(i=1; i<Nhex+1; i++){
		energy0=sqrt(fabs(E1*E1+E2*E2+E3*E3+2*E1*E2*cos(med1*i*(n+2*m)-
			n*med2*(1+straint/100)*(1-p)*interp)+2*E1*E3*cos(med1*i*
			(2*n+m)+m*med2*(1+straint/100)*(1-p)*interp)+2*E2*E3*
			cos(med1*i*(n-m)+(n+m)*med2*(1+straint/100)*(1-p)*interp)));
		
		for(j=1-p; j<p; j++){
			kt=(j+1)*interp;
			med3=n*kt*med2*(1+straint/100);
			med4=m*kt*med2*(1+straint/100);
			
			energy1=sqrt(fabs(E1*E1+E2*E2+E3*E3+2*E1*E2*cos(med1*i*(n+2*m)
					-med3)+2*E1*E3*cos(med1*i*(2*n+m)+med4)+2*E2*E3
					*cos(med1*i*(n-m)+med3+med4)));
			
			if((energy0>=0) && (energy0<Emax)){
				channel=(int)(energy0/Emax*iter);
				slope=fabs(energy1-energy0)/interp;
				
				if(slope==0){
					result=100;
				}
				else{
					result=1/slope;
				}
				data[channel]+=result;
			}
			energy0=energy1;
		}
	}
	
	fprintf(fp2, "Energy\t DOS\t 2*Energy\n");
	
	for(i=0; i<iter; i++){
		fprintf(fp2, "%f\t %f\t %f\n", i*Emax/iter, data[i]*2/Nhex/(2*p-1),
				2*i*Emax/iter);
	}
	
	free(data);
	fclose(fp2);

	printf("\n");
	main();
}

int gcdn(int n, int m){
	int result, i;
	if(m==0){
		result=n;
	}
	if(m>0){
		for(i=1; i<m+1; i++){
			if( ((float)n/(float)i==(int)((float)n/(float)i)) &&
					((float)m/(float)i==(int)((float)m/(float)i)) ){
				result=i;
			}
		}
	}
	return(result);
}
