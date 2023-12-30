function error = funerr5_6_hpc76(para)

global xx12 yy12 eps ntotvect

solver = 'bvp4c'; 
bvpsolver = fcnchk(solver);   
    
% Problem parameter, shared with nested functions

TD=0.01;

L = 600; va = 0.5; vr = 0.5; 

%gamma10 = 0.093; gamma01 = 0.041;
gamma10 = para(1); gamma01 = para(2);

%gammaar = 3.1*0.00001; gammara = 6.9*0.00001; 
gammaar = para(3); gammara = para(4);
 
%gammaona = 2.75*0.0001; gammaonr = 2.75*0.0001;
gammaona=para(5); gammaonr=para(6);

%gammaoffa = 0.001*4.45*0.001;gammaoffr = 0.001*4.45*0.001;
gammaoffa = para(7); gammaoffr = para(8);

%Dfree = 3; T12 =2.16*100000;
Dfree=0.5; T12 = 5.01e5;

%jtot0=0.1; nfree0=0.5;
jtot0=para(9); nfree0=para(10);  A5=para(11);
 
% Initial mesh - duplicate the interface point x = 1.
xinit1 = [0, 0.25*L, 0.5*L, 0.75*L, L]; 

% Constant initial guess for the solution
yinit1 = [1;1;1;1];

% The initial profile
sol1 = bvpinit(xinit1,yinit1);
  sol1 = bvpsolver(@f1,@(YL,YR)bc1(YL,YR,jtot0,nfree0,A5),sol1);  
%sol1 = bvpsolver(@f1,@bc1,sol1);  
  
na=sol1.y(1,:);
nr=sol1.y(2,:);
nfree=sol1.y(3,:);

na0=(gamma10*(gamma01 + gammaoffr + gammara)*sol1.y(1,:)...
    + gamma10*gammara*sol1.y(2,:)... 
  +((gamma01 + gammaoffr)*gammaona + (gammaona + gammaonr)*gammara)...
*sol1.y(3,:))/((gamma01 + gammaar + gammaoffa)...
*(gamma01 + gammaoffr) + (gamma01 + gammaoffa)*gammara);

nr0=(gamma10*gammaar*sol1.y(1,:) + gamma10*(gamma01 + gammaar...
    + gammaoffa)*sol1.y(2,:)...
  +(gammaar*gammaona + (gamma01 + gammaar + gammaoffa)*gammaonr)...
  *sol1.y(3,:))/((gamma01 + gammaar + gammaoffa)*(gamma01 + gammaoffr)...
  + (gamma01 + gammaoffa)*gammara);

ntot=na+nr+na0+nr0+nfree;
jtot=-Dfree*sol1.y(4,:)+va*sol1.y(1,:)-vr*sol1.y(2,:);
nmt=na+nr+na0+nr0;
nmtproc=100*nmt./ntot;

xx11=xx12;

for i=1:76
ntotp(i)=interp1(sol1.x,ntot,xx11(i));
end

for i=1:76
vvp(i)=interp1(sol1.x,jtot./ntot,xx11(i));
end

nmtprocp=interp1(sol1.x,nmtproc,xx11(38));

dev1=0;
for i=1:76
    dev1=dev1+(ntotp(i)-yy12(i))^2;
end

for i=1:76
eps(i)=yy12(i)-ntotp(i);
ntotvect(i)=ntotp(i);
end

dev2=0;
for i=1:76
dev2=dev2+(vvp(i)-0.0926)^2;
end

dev3=(nmtprocp-100)^2;

error=dev1+10*dev2+10*dev3;
 
 %ylabel('$$\rm{j}$$','Interpreter','latex','FontSize',12);
    
  % -----------------------------------------------------------------------
  % Nested functions -- n and eta are provided by the outer function.
  %
%-----------------------------
function dydx = f1(x,y,region)
  % Derivative function -- share n and eta with the outer function.
    dydx = zeros(4,1);
    % The definition depends on the region.
    
    % y(1)=na
    % y(2)=nr
    % y(3)=nfree
    % y(4)=nfree'
    
 dydx(1)=(-(gamma10*((gammaar + gammaoffa)*(gamma01 + gammaoffr)...
     + gammaoffa*gammara)*y(1))... 
    +gamma01*gamma10*gammara*y(2) + gamma01...
*((gamma01 + gammaoffr)*gammaona + (gammaona + gammaonr)*gammara)*y(3))...
   /((gamma01 + gammaar + gammaoffa)*(gamma01 + gammaoffr)*va...
   + (gamma01 + gammaoffa)*gammara*va);
 
 dydx(2)= (gamma10*y(2) - (gamma01*(gamma10*gammaar*y(1)... 
    +gamma10*(gamma01 + gammaar + gammaoffa)*y(2)... 
 +(gammaar*gammaona + (gamma01 + gammaar + gammaoffa)*gammaonr)*y(3)))...
   /((gamma01 + gammaar + gammaoffa)*(gamma01 + gammaoffr)...
   + (gamma01 + gammaoffa)*gammara))/vr;
    
  dydx(3)=y(4);
  
  dydx(4)= (-(gamma10*(gammaar*gammaoffr + gammaoffa*(gamma01...
      + gammaoffr + gammara))*T12*y(1))... 
  -gamma10*((gamma01 + gammaar + gammaoffa)*gammaoffr...
  + gammaoffa*gammara)*T12*y(2)... 
   +(gamma01*(gammaoffr*gammaona + gammaoffa*gammaonr...
   + gamma01*(gammaona + gammaonr)... 
  +gammaar*(gammaona + gammaonr) + (gammaona + gammaonr)*gammara)*T12... 
   +((gamma01 + gammaar + gammaoffa)*(gamma01 + gammaoffr)...
   + (gamma01 + gammaoffa)*gammara)...
    *log(2))*y(3))...
   /(Dfree*((gamma01 + gammaar + gammaoffa)...
   *(gamma01 + gammaoffr) + (gamma01 + gammaoffa)*gammara)*T12);
    
  end
  % -----------------------------------------------------------------------
    
     function res = bc1(YL,YR,jtot0,nfree0,A5)
  % Boundary (and internal) conditions
     res = [-Dfree*YL(4)+va*YL(1)-vr*YL(2)-jtot0 % x = 0
      YL(3)-nfree0           % x = 0
           
      YR(4)  % x = L
-Dfree*YR(4)+va*YR(1)-vr*YR(2)-A5*(1-exp(-log(2)/(T12*gammaar)))*va*YR(1)];       
  end
end
    
   