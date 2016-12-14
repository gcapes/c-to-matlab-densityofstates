% %%%%%% Calculate the energy and density of states %%%%%
clc
clear all 

E0=2.9; % Hopping parameter without strain */
a0=0.242; %Ao Length of the lattice vector without strain */
n=5; 
m=5; % Nanotube parameters

straint=0.0; %Strain along the nanotube axis
strainc=0; %Strain along the nanotube circumference
circum=sqrt(n*n+m*m+n*m);
dr=gcd(2*n+m,2*m+n);
Nhex=2*(n*n+m*m+n*m)/dr
period=(1+straint/100)*sqrt(3)*circum*a0/dr;
r1=sqrt(((n+m)^2*a0^2*(1+strainc/100)^2)/(4*circum^2)+((n-m)^2*a0^2*(1+straint/100)^2)/(12*circum^2));
r2=sqrt((m^2*a0^2*(1+strainc/100)^2)/(4*circum^2)+((2*n+m)^2*a0^2*(1+straint/100)^2)/(12*circum^2));
r3=sqrt((n^2*a0^2*(1+strainc/100)^2)/(4*circum^2)+((n+2*m)^2*a0^2*(1+straint/100)^2)/(12*circum^2));
E1=E0*(a0/(sqrt(3)*r1))^2;
E2=E0*(a0/(sqrt(3)*r2))^2;
E3=E0*(a0/(sqrt(3)*r3))^2;

T=period;
h=0.01;%
kt=-pi/T:h:pi/T;
% Calculate the energy eV
for k=1:length(kt)
	for j=1:Nhex
		e1(k,j)=sqrt(abs(E1^2+E2^2+E3^2+2*E1*E2*cos(pi*j*((n+2*m)/circum^2)+((sqrt(3)*n*a0/...
		(2*circum))*(1+straint)*kt(k)))+2*E1*E3*cos((pi*j*(n-m)/circum^2)-(sqrt(3)*a0*(n+m)/...
		(2*circum))*(1+straint)*kt(k))+2*E1*E3*cos((pi*j*(2*n+m)/circum^2)-((sqrt(3)*m*a0/...
		(2*circum))*(1+straint)*kt(k)))));

		e2(k,j)=-sqrt(abs(E1^2+E2^2+E3^2+2*E1*E2*cos(pi*j*((n+2*m)/circum^2)+((sqrt(3)*n*a0/...
		(2*circum))*(1+straint)*kt(k)))+2*E1*E3*cos((pi*j*(n-m)/circum^2)-(sqrt(3)*a0*(n+m)/...
		(2*circum))*(1+straint)*kt(k))+2*E1*E3*cos((pi*j*(2*n+m)/circum^2)-((sqrt(3)*m*a0/...
		(2*circum))*(1+straint)*kt(k)))));
	end 
end;
figure(1)
plot(kt,e1,kt,e2)
