function error = funerr5_6_figure(para)

global Dfree L T12 va vr gamma10 gamma01 gammaona gammaoffr jtot0 nfree0 ...
    A5 xx12 yy12 eps ntotvect

A = importdata('fig3Da.txt');
xx12=A(:,1); yy12=A(:,2);

x_vr_0_05=[];
na_vr_0_05=[];
nr_vr_0_05=[];
na0_vr_0_05=[];
nr0_vr_0_05=[];
nfree_vr_0_05=[];
ntot_vr_0_05=[];
nmtproc_vr_0_05=[];
jtot_vr_0_05=[];
v_vr_0_05=[];
jtot_vr_0_05=[];
jretr_vr_0_05=[];
janter_vr_0_05=[];
jdif_vr_0_05=[];

x_Dfree_0_05=[];
na_Dfree_0_05=[];
nr_Dfree_0_05=[];
na0_Dfree_0_05=[];
nr0_Dfree_0_05=[];
nfree_Dfree_0_05=[];
ntot_Dfree_0_05=[];
nmtproc_Dfree_0_05=[];
jtot_Dfree_0_05=[];
v_Dfree_0_05=[];
jtot_Dfree_0_05=[];
jretr_Dfree_0_05=[];
janter_Dfree_0_05=[];
jdif_Dfree_0_05=[];

load('full_Dfree_0_05.mat');
load('full_vr_0_05.mat');

solver = 'bvp4c'; 
bvpsolver = fcnchk(solver);   
    
% Problem parameter, shared with nested functions

TD=0.01;

L = 600; va = 0.5; 
vr = 0.5; 
% vr = 0.05; 

%gamma10 = 0.093; gamma01 = 0.041;
gamma10 = para(1); gamma01 = para(2);

%gammaar = 3.1*0.00001; gammara = 6.9*0.00001; 
gammaar = para(3); gammara = para(4);
 
%gammaona = 2.75*0.0001; gammaonr = 2.75*0.0001;
gammaona=para(5); gammaonr=para(6);

%gammaoffa = 0.001*4.45*0.001;gammaoffr = 0.001*4.45*0.001;
gammaoffa = para(7); gammaoffr = para(8);

%Dfree = 3; T12 =2.16*100000;
Dfree=0.5; 
% Dfree=0.05; 

T12 = 5.01e5;

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
janter=va*sol1.y(1,:);
jretr=-vr*sol1.y(2,:);
jdif=-Dfree*sol1.y(4,:);

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

close all

figure(1)
 hold on
    h1=plot(sol1.x,na);
    h2=plot(sol1.x,nr);
    h3=plot(sol1.x,na0);
    h4=plot(sol1.x,nr0);
    h5=plot(sol1.x,nfree);
    h=legend('n^*_a, D_{free}=0.5 \mum^2/s, v_r=0.5 \mum/s',...
        'n^*_r, D_{free}=0.5 \mum^2/s, v_r=0.5 \mum/s',...
        'n^*_{a0}, D_{free}=0.5 \mum^2/s, v_r=0.5 \mum/s',...
        'n^*_{r0}, D_{free}=0.5 \mum^2/s, v_r=0.5 \mum/s',...
        'n^*_{free}, D_{free}=0.5 \mum^2/s, v_r=0.5 \mum/s',...
        'Location','northwest');
    set(h1,'linewidth',1.5,'color','r','LineStyle','-');
    set(h2,'linewidth',1.5,'color','b','LineStyle','--');
    set(h3,'linewidth',1.5,'color','g','LineStyle','-.');
    set(h4,'linewidth',1.5,'color','c','LineStyle','-');
    set(h5,'linewidth',1.5,'color','m','LineStyle','--');
    ylabel('n^*, full model','FontSize',12)
    h20=xlabel('x, \mum','FontSize',12);
ylim([0 2.5])
 hold off
 
 figure(2)
 hold on
    h1=plot(x_Dfree_0_05,na_Dfree_0_05);
    h2=plot(x_Dfree_0_05,nr_Dfree_0_05);
    h3=plot(x_Dfree_0_05,na0_Dfree_0_05);
    h4=plot(x_Dfree_0_05,nr0_Dfree_0_05);
    h5=plot(x_Dfree_0_05,nfree_Dfree_0_05);
    h=legend('n^*_a, D_{free}=0.05 \mum^2/s, v_r=0.5 \mum/s',...
        'n^*_r, D_{free}=0.05 \mum^2/s, v_r=0.5 \mum/s',...
        'n^*_{a0}, D_{free}=0.05 \mum^2/s, v_r=0.5 \mum/s',...
        'n^*_{r0}, D_{free}=0.05 \mum^2/s, v_r=0.5 \mum/s',...
        'n^*_{free}, D_{free}=0.05 \mum^2/s, v_r=0.5 \mum/s',...
        'Location','northwest');
    set(h1,'linewidth',1.5,'color','r','LineStyle','-');
    set(h2,'linewidth',1.5,'color','b','LineStyle','--');
    set(h3,'linewidth',1.5,'color','g','LineStyle','-.');
    set(h4,'linewidth',1.5,'color','c','LineStyle','-');
    set(h5,'linewidth',1.5,'color','m','LineStyle','--');
    ylabel('n^*, full model','FontSize',12)
    h20=xlabel('x, \mum','FontSize',12);
ylim([0 2.5])
 hold off
 
 figure(3)
 hold on
    h1=plot(x_vr_0_05,na_vr_0_05);
    h2=plot(x_vr_0_05,nr_vr_0_05);
    h3=plot(x_vr_0_05,na0_vr_0_05);
    h4=plot(x_vr_0_05,nr0_vr_0_05);
    h5=plot(x_vr_0_05,nfree_vr_0_05);
    h=legend('n^*_a, D_{free}=0.5 \mum^2/s, v_r=0.05 \mum/s',...
        'n^*_r, D_{free}=0.5 \mum^2/s, v_r=0.05 \mum/s',...
        'n^*_{a0}, D_{free}=0.5 \mum^2/s, v_r=0.05 \mum/s',...
        'n^*_{r0}, D_{free}=0.5 \mum^2/s, v_r=0.05 \mum/s',...
        'n^*_{free}, D_{free}=0.5 \mum^2/s, v_r=0.05 \mum/s',...
        'Location','northwest');
    set(h1,'linewidth',1.5,'color','r','LineStyle','-');
    set(h2,'linewidth',1.5,'color','b','LineStyle','--');
    set(h3,'linewidth',1.5,'color','g','LineStyle','-.');
    set(h4,'linewidth',1.5,'color','c','LineStyle','-');
    set(h5,'linewidth',1.5,'color','m','LineStyle','--');
    ylabel('n^*, full model','FontSize',12)
    h20=xlabel('x, \mum','FontSize',12);
ylim([0 2.5])
 hold off
 
 figure(4)
 hold on
 h2=plot(sol1.x,janter);
 h3=plot(sol1.x,jretr);
 h4=plot(sol1.x,jdif);
 h5=plot(sol1.x,jtot);
 
 hold off
 set(h2,'linewidth',1.5,'color','r','LineStyle','-');
    set(h3,'linewidth',1.5,'color','b','LineStyle','--');
    set(h4,'linewidth',1.5,'color','g','LineStyle','-.');
    set(h5,'linewidth',1.5,'color','m','LineStyle',':');
    
 h=legend('j^*_{anter}, D_{free}=0.5 \mum^2/s, v_r=0.5 \mum/s',...
     'j^*_{retr}, D_{free}=0.5 \mum^2/s, v_r=0.5 \mum/s',...
     'j^*_{dif}, D_{free}=0.5 \mum^2/s, v_r=0.5 \mum/s',...
     'j^*_{tot}, D_{free}=0.5 \mum^2/s, v_r=0.5 \mum/s',...
 'Location','south');
 ylabel('j^*, full model','FontSize',12)
    h20=xlabel('x, \mum','FontSize',12);
    
    figure(5)
 hold on
 
 h6=plot(x_Dfree_0_05,janter_Dfree_0_05);
 h7=plot(x_Dfree_0_05,jretr_Dfree_0_05);
 h8=plot(x_Dfree_0_05,jdif_Dfree_0_05);
 h9=plot(x_Dfree_0_05,jtot_Dfree_0_05);
 
 hold off
    
    set(h6,'linewidth',1.5,'color','r','LineStyle','-');
    set(h7,'linewidth',1.5,'color','b','LineStyle','--');
     set(h8,'linewidth',1.5,'color','g','LineStyle','-.');
     set(h9,'linewidth',1.5,'color','m','LineStyle',':');
     
 h=legend('j^*_{anter}, D_{free}=0.05 \mum^2/s, v_r=0.5 \mum/s',...
     'j^*_{retr}, D_{free}=0.05 \mum^2/s, v_r=0.5 \mum/s',...
     'j^*_{dif}, D_{free}=0.05 \mum^2/s, v_r=0.5 \mum/s',...
     'j^*_{tot}, D_{free}=0.05 \mum^2/s, v_r=0.5 \mum/s',...
 'Location','south');
 ylabel('j^*, full model','FontSize',12)
    h20=xlabel('x, \mum','FontSize',12);
    
    figure(6)
 hold on
 
 h10=plot(x_vr_0_05,janter_vr_0_05);
 h11=plot(x_vr_0_05,jretr_vr_0_05);
  h12=plot(x_vr_0_05,jdif_vr_0_05);
 h13=plot(x_vr_0_05,jtot_vr_0_05);
 
 hold off
 
     set(h10,'linewidth',1.5,'color','r','LineStyle','-');
    set(h11,'linewidth',1.5,'color','b','LineStyle','--');
     set(h12,'linewidth',1.5,'color','g','LineStyle','-.');
     set(h13,'linewidth',1.5,'color','m','LineStyle',':');
     
 h=legend('j^*_{anter}, D_{free}=0.5 \mum^2/s, v_r=0.05 \mum/s',...
     'j^*_{retr}, D_{free}=0.5 \mum^2/s, v_r=0.05 \mum/s',...
     'j^*_{dif}, D_{free}=0.5 \mum^2/s, v_r=0.05 \mum/s',...
     'j^*_{tot}, D_{free}=0.5 \mum^2/s, v_r=0.05 \mum/s',...
 'Location','south');
 ylabel('j^*, full model','FontSize',12)
    h20=xlabel('x, \mum','FontSize',12);
    
    figure(7)
 hold on
h7=plot(sol1.x,ntot);
h8=plot(x_Dfree_0_05,ntot_Dfree_0_05);
h9=plot(x_vr_0_05,ntot_vr_0_05);

 set(h7,'linewidth',1.5,'color','r','LineStyle','-');
  set(h8,'linewidth',1.5,'color','g','LineStyle','--');
  set(h9,'linewidth',1.5,'color','b','LineStyle','-.');
x8={0.001 50 100 150 200 250 300 350 400 450 500 550 600};
y8={TD*18.5*300/55 TD*13*300/55 TD*7*300/55 TD*7*300/55 ...
    TD*7*300/55 TD*7*300/55 ...
    TD*7*300/55 TD*7*300/55 TD*7*300/55 TD*7*300/55 ...
    TD*7*300/55 TD*17*300/55 TD*25*300/55};

plot(xx12, yy12,'o','color','k')

   h=legend( 'D_{free}=0.5 \mum^2/s, v_r=0.5 \mum/s',...
   'D_{free}=0.05 \mum^2/s, v_r=0.5 \mum/s',...
       'D_{free}=0.5 \mum^2/s, v_r=0.05 \mum/s',...
       'data from Fig. 3D of Black et al.','Location','northeast');
    
    ylabel('n^*_{tot}, full model','FontSize',12)
    h20=xlabel('x, \mum','FontSize',12);
  %   axis([0 L 0 150]);
 hold off
 
% x_vr_0_05=sol1.x; na_vr_0_05=na; nr_vr_0_05=nr; na0_vr_0_05=na0; 
% nr0_vr_0_05=nr0;
% nfree_vr_0_05=nfree; ntot_vr_0_05=ntot; nmtproc_vr_0_05=nmtproc;jtot_vr_0_05=jtot;
% v_vr_0_05=jtot./ntot; janter_vr_0_05=janter; jretr_vr_0_05=jretr; 
% jdif_vr_0_05=jdif;
% save('full_vr_0_05.mat','x_vr_0_05','na_vr_0_05','nr_vr_0_05','na0_vr_0_05','nr0_vr_0_05',...
%     'nfree_vr_0_05','ntot_vr_0_05','nmtproc_vr_0_05',...
%     'v_vr_0_05',...
% 'janter_vr_0_05','jretr_vr_0_05','jdif_vr_0_05','jtot_vr_0_05');

% x_Dfree_0_05=sol1.x; na_Dfree_0_05=na; nr_Dfree_0_05=nr; na0_Dfree_0_05=na0; 
% nr0_Dfree_0_05=nr0;
% nfree_Dfree_0_05=nfree; ntot_Dfree_0_05=ntot; nmtproc_Dfree_0_05=nmtproc;jtot_Dfree_0_05=jtot;
% v_Dfree_0_05=jtot./ntot;
% janter_Dfree_0_05=janter; jretr_Dfree_0_05=jretr; jdif_Dfree_0_05=jdif;
% save('full_Dfree_0_05.mat','x_Dfree_0_05','na_Dfree_0_05','nr_Dfree_0_05','na0_Dfree_0_05','nr0_Dfree_0_05',...
%     'nfree_Dfree_0_05','ntot_Dfree_0_05','nmtproc_Dfree_0_05',...
%     'v_Dfree_0_05',...
% 'janter_Dfree_0_05','jretr_Dfree_0_05','jdif_Dfree_0_05','jtot_Dfree_0_05');
    

function dydx = f1(x,y,region)
 
    dydx = zeros(4,1);
    
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