% Porting of the original C code into MATLAB
% Gerard Capes <gerard.capes@manchester.ac.uk> 14 December 2017
% Run the main function like this:
    % DOS(5,5,0,0,1000,1000);

function DOS(n,m,straint,strainc,iter,p)
    % n, m are nanotube parameters
    % straint: Strain along the nanotube axis (%)
    % strainc: Strain along the nanotube circumference (%)
    % iter is number of points in the graph
    % p is number of points analyzed

    % Declaration of variables
    % Nhex = Number of hexagons in the 1D unit cell
    % circum = Circumference of the nanotube in unit of a0
    % period = Period of the nanotube along its axis

    E0=2.9; % Hopping parameter without strain
    a0=2.49; % Length of the lattice vector without strain
    % E1, E2, E3 = New hopping parameters with strain
    % r1, r2, r3 = Strained bond lengths
    

    circum=sqrt(n*n+m*m+n*m);
    dr=gcdn(2*n+m,2*m+n);
    Nhex=2*(n*n+m*m+n*m)/dr;
    period=(1+straint/100)*sqrt(3)*circum*a0/dr;
    r1=sqrt((n+m)*(n+m)*a0*a0*(1+strainc/100)*(1+strainc/100)/4/ ...
        circum/circum+(n-m)*(n-m)*a0*a0*(1+straint/100)* ...
        (1+straint/100)/12/circum/circum);
    r2=sqrt(m*m*a0*a0*(1+strainc/100)*(1+strainc/100)/4/circum/circum+ ...
        (2*n+m)*(2*n+m)*a0*a0*(1+straint/100)*(1+straint/100) ...
        /12/circum/circum);
    r3=sqrt(n*n*a0*a0*(1+strainc/100)*(1+strainc/100)/4/circum/circum+ ...
        (n+2*m)*(n+2*m)*a0*a0*(1+straint/100)*(1+straint/100)/12/ ...
        circum/circum);
    E1=E0*a0*a0/3/r1/r1;
    E2=E0*a0*a0/3/r2/r2;
    E3=E0*a0*a0/3/r3/r3;

    density(iter,p);
    
    function density(iter,p)
        % Nested function to calculate the electronic density of states
        % Nested so it has access to variables in DOS
        
        % iter is number of points in the graph
        % p is number of points analyzed

        Emax=9.0; % Energy range to take into account
        PI=3.141592653;

        fp2=fopen('density.txt','w');

        data=zeros(iter,1);
        interp=PI/period/p;
        med1=PI/circum/circum;
        med2=sqrt(3)/2/circum*a0;

        for i=1:1:Nhex
            energy0=sqrt(abs(E1*E1+E2*E2+E3*E3+2*E1*E2*cos(med1*i*(n+2*m)- ...
                n*med2*(1+straint/100)*(1-p)*interp)+2*E1*E3*cos(med1*i* ...
                (2*n+m)+m*med2*(1+straint/100)*(1-p)*interp)+2*E2*E3* ...
                cos(med1*i*(n-m)+(n+m)*med2*(1+straint/100)*(1-p)*interp)));

            for j=1-p:1:p-1
                kt=(j+1)*interp;
                med3=n*kt*med2*(1+straint/100);
                med4=m*kt*med2*(1+straint/100);

                energy1=sqrt(abs(E1*E1+E2*E2+E3*E3+2*E1*E2*cos(med1*i*(n+2*m) ...
                    -med3)+2*E1*E3*cos(med1*i*(2*n+m)+med4)+2*E2*E3 ...
                    *cos(med1*i*(n-m)+med3+med4)));

                if((energy0>=0) && (energy0<Emax))
                    % +1 because MATLAB indexes from 1
                    channel=floor(energy0/Emax*iter)+1; 
                    slope=abs(energy1-energy0)/interp;

                    if(slope==0)
                        result=100;
                    else
                        result=1/slope;
                    end
                    data(channel)=data(channel)+result;
                end
                energy0=energy1;
            end
        end

        fprintf(fp2, 'Energy\t DOS\t 2*Energy\n');

        for(i=1:1:iter) % MATLAB indexes from 1
            fprintf(fp2, '%f\t %f\t %f\n', (i-1)*Emax/iter, data(i)*2/Nhex ...
                /(2*p-1),2*(i-1)*Emax/iter);
        end

        fclose(fp2);
    end
end

function [result]=gcdn(n,m)
    % Program to calculate the gcdn of n and m
    if(m==0)
        result=n;
    end
    if(m>0)
        for(i=1:1:m)
            if(mod(n,i)==0 && mod(m,i)==0)
                result=i;
            end
        end
    end
end

