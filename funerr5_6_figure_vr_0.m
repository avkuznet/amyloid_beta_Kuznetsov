function error = funerr5_6_figure_vr_0(para)

global Dfree L T12 va vr gamma10 gamma01 gammaona gammaoffr jtot0 nfree0 ...
    A5 xx12 yy12 eps ntotvect

A = importdata('fig3Da.txt');
xx12=A(:,1); yy12=A(:,2);

x0_05=[];
na0_05=[];
nr0_05=[];
na00_05=[];
nr00_05=[];
nfree0_05=[];
ntot0_05=[];
nmtproc0_05=[];
jtot0_05=[];
v0_05=[];
janter0_05=[];
jdif0_05=[];
jtot0_05=[];

load('Dfree_0_05.mat');

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
Dfree=0.5; 
% Dfree=0.05; 
T12 = 5.01e5;

%jtot0=0.1; nfree0=0.5;
jtot0=para(9); nfree0=para(10);  A5=para(11);
 
% Initial mesh - duplicate the interface point x = 1.
xinit1 = [0, 0.25*L, 0.5*L, 0.75*L, L]; 

% Constant initial guess for the solution
yinit1 = [1;1;1];

% The initial profile
sol1 = bvpinit(xinit1,yinit1);
  sol1 = bvpsolver(@f1,@(YL,YR)bc1(YL,YR,jtot0,nfree0,A5),sol1);  
%sol1 = bvpsolver(@f1,@bc1,sol1);  
  
na=sol1.y(1,:);
nfree=sol1.y(2,:);

% jtot0117=jtot0
% jtotL117=A5*(1-exp(-log(2)/(T12*gammaar)))*na(end)*va
% pause

nr=(gamma10.*(gamma01+gammaar+gammaoffa).*gammaoffr+gamma10.*( ...
  gamma01+gammaoffa).*gammara).^(-1).*(gamma01.*gamma10.* ...
  gammaar.*sol1.y(1,:)+gamma01.*(gammaar.*gammaona+(gamma01+gammaar+ ...
  gammaoffa).*gammaonr).*sol1.y(2,:));

na0=((gamma01+gammaar+gammaoffa).*gammaoffr+(gamma01+gammaoffa) ...
  .*gammara).^(-1).*(gamma10.*(gammaoffr+gammara).*sol1.y(1,:)+ ...
  gammaoffr.*gammaona.*sol1.y(2,:)+(gammaona+gammaonr).*gammara.* ...
  sol1.y(2,:));

nr0=((gamma01+gammaar+gammaoffa).*gammaoffr+(gamma01+gammaoffa) ...
  .*gammara).^(-1).*(gamma10.*gammaar.*sol1.y(1,:)+gammaar.*gammaona.* ...
  sol1.y(2,:)+(gamma01+gammaar+gammaoffa).*gammaonr.*sol1.y(2,:));

ntot=na+nr+na0+nr0+nfree;

jtot=(-Dfree*sol1.y(3,:)+va*sol1.y(1,:));
jdif=-Dfree*sol1.y(3,:);
janter=va*sol1.y(1,:);

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
    h=legend('n^{*(0)}_a, D_{free}=0.5 \mum^2/s','n^{*(0)}_r, D_{free}=0.5 \mum^2/s',...
        'n^{*(0)}_{a0}, D_{free}=0.5 \mum^2/s','n^{*(0)}_{r0}, D_{free}=0.5 \mum^2/s',...
        'n^{*(0)}_{free}, D_{free}=0.5 \mum^2/s','Location','northwest');
    set(h1,'linewidth',1.5,'color','r','LineStyle','-');
    set(h2,'linewidth',1.5,'color','b','LineStyle','--');
    set(h3,'linewidth',1.5,'color','g','LineStyle','-.');
    set(h4,'linewidth',1.5,'color','c','LineStyle','-');
    set(h5,'linewidth',1.5,'color','m','LineStyle','--');
    ylabel('n^{*(0)}, v^*_r=0','FontSize',12)
    h20=xlabel('x, \mum','FontSize',12);
ylim([0 8])
 hold off
 
 figure(2)
 hold on
    h1=plot(x0_05,na0_05);
    h2=plot(x0_05,nr0_05);
    h3=plot(x0_05,na00_05);
    h4=plot(x0_05,nr00_05);
    h5=plot(x0_05,nfree0_05);
    h=legend('n^{*(0)}_a, D_{free}=0.05 \mum^2/s','n^{*(0)}_r, D_{free}=0.05 \mum^2/s',...
        'n^{*(0)}_{a0}, D_{free}=0.05 \mum^2/s','n^{*(0)}_{r0}, D_{free}=0.05 \mum^2/s',...
        'n^{*(0)}_{free}, D_{free}=0.05 \mum^2/s','Location','northwest');
    set(h1,'linewidth',1.5,'color','r','LineStyle','-');
    set(h2,'linewidth',1.5,'color','b','LineStyle','--');
    set(h3,'linewidth',1.5,'color','g','LineStyle','-.');
    set(h4,'linewidth',1.5,'color','c','LineStyle','-');
    set(h5,'linewidth',1.5,'color','m','LineStyle','--');
    ylabel('n^{*(0)}, v^*_r=0','FontSize',12)
    h20=xlabel('x, \mum','FontSize',12);
ylim([0 8])
 hold off
 
 figure(3)
 hold on
 h2=plot(sol1.x,janter);
 h3=plot(sol1.x,jdif);
 h4=plot(sol1.x,jtot);
 
 h5=plot(x0_05,janter0_05);
 h6=plot(x0_05,jdif0_05);
 h7=plot(x0_05,jtot0_05);
 
 hold off
 set(h2,'linewidth',1.5,'color','r','LineStyle','-');
    set(h3,'linewidth',1.5,'color','b','LineStyle','-');
    set(h4,'linewidth',1.5,'color','g','LineStyle','-');
    set(h5,'linewidth',1.5,'color','r','LineStyle','--');
    set(h6,'linewidth',1.5,'color','b','LineStyle','--');
    set(h7,'linewidth',1.5,'color','g','LineStyle','--');
 h=legend('j^*_{anter}, D_{free}=0.5 \mum^2/s',...
     'j^*_{dif}, D_{free}=0.5 \mum^2/s',...
     'j^*_{tot}, D_{free}=0.5 \mum^2/s',...
     'j^*_{anter}, D_{free}=0.05 \mum^2/s',...
     'j^*_{dif}, D_{free}=0.05 \mum^2/s',...
     'j^*_{tot}, D_{free}=0.05 \mum^2/s',...
 'Location','south');
 ylabel('j^*, v^*_r=0','FontSize',12)
    h20=xlabel('x, \mum','FontSize',12);
 
 figure(4)
 hold on
h7=plot(sol1.x,ntot);
h8=plot(x0_05,ntot0_05);
 set(h7,'linewidth',1.5,'color','r','LineStyle','-');
  set(h8,'linewidth',1.5,'color','g','LineStyle','--');
x8={0.001 50 100 150 200 250 300 350 400 450 500 550 600};
y8={TD*18.5*300/55 TD*13*300/55 TD*7*300/55 TD*7*300/55 ...
    TD*7*300/55 TD*7*300/55 ...
    TD*7*300/55 TD*7*300/55 TD*7*300/55 TD*7*300/55 ...
    TD*7*300/55 TD*17*300/55 TD*25*300/55};

plot(xx12, yy12,'o','color','k')

   h=legend('D_{free}=0.5 \mum^2/s','D_{free}=0.05 \mum^2/s'...
       ,'data from Fig. 3D of Black et al.','Location','northeast');
    ylabel('n^{*(0)}_{tot}, v^*_r=0','FontSize',12)
    h20=xlabel('x, \mum','FontSize',12);
  ylim([0 6])
 hold off
 
% x0_05=sol1.x; na0_05=na; nr0_05=nr; na00_05=na0; nr00_05=nr0;
% nfree0_05=nfree; ntot0_05=ntot; nmtproc0_05=nmtproc;jtot0_05=jtot;
% v0_05=jtot./ntot;
% jtot0_05=jtot;
% jdif0_05=jdif;
% janter0_05=janter;
% save('Dfree_0_05.mat','x0_05','na0_05','nr0_05','na00_05','nr00_05',...
%     'nfree0_05','ntot0_05','nmtproc0_05',...
%     'v0_05','jtot0_05','jdif0_05','janter0_05');
 
    
function dydx = f1(x,y,region)
 
    dydx = zeros(3,1);
    
    % y(1)=na
    % y(2)=nfree
    % y(3)=nfree'
    
 dydx(1)=((-1).*gamma10.*(gammaar.*gammaoffr+gammaoffa.*(gammaoffr+ ...
  gammara)).*y(1)+gamma01.*(gammaoffr.*gammaona+(gammaona+ ...
  gammaonr).*gammara).*y(2)).*((gamma01+gammaar+gammaoffa).* ...
  gammaoffr.*va+(gamma01+gammaoffa).*gammara.*va).^(-1);

dydx(2)=y(3);

dydx(3)=Dfree.^(-1).*((gamma01+gammaar+gammaoffa).*gammaoffr+( ...
  gamma01+gammaoffa).*gammara).^(-1).*T12.^(-1).*((-1).* ...
  gamma10.*(gammaar.*gammaoffr+gammaoffa.*(gammaoffr+gammara)) ...
  .*y(1).*T12+gamma01.*(gammaoffr.*gammaona+(gammaona+gammaonr) ...
  .*gammara).*y(2).*T12+((gamma01+gammaar+gammaoffa).* ...
  gammaoffr+(gamma01+gammaoffa).*gammara).*y(2).*log(2));
    
  end
  % -----------------------------------------------------------------------
    
      function res = bc1(YL,YR,jtot0,nfree0,A5)
  % Boundary (and internal) conditions
  % y(1)=na
    % y(2)=nfree
    % y(3)=nfree'
     res = [-Dfree*YL(3)+va*YL(1)-jtot0 % x = 0
         YL(2)-nfree0           % x = 0
  -Dfree*YR(3)+va*YR(1)-A5*(1-exp(-log(2)/(T12*gammaar)))*va*YR(1)];          
%       YR(3)  % x = L
   
       
  end

end