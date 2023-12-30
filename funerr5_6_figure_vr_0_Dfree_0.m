function error = funerr5_6_figure_vr_0_Dfree_0(para)

global Dfree L T12 va vr gamma10 gamma01 gammaona gammaoffr jtot0 nfree0 ...
    A5 xx12 yy12 eps ntotvect

A = importdata('fig3Da.txt');
xx12=A(:,1); yy12=A(:,2);

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
 
xanal=linspace(0,L,2000);

na=2.^((-1).*gamma10.*(gammaar.*gammaoffr+gammaoffa.*( ...
  gammaoffr+gammara)).*xanal.*(gamma01.*(gammaoffr.*gammaona+( ...
  gammaona+gammaonr).*gammara).*T12.*va+((gamma01+gammaar+ ...
  gammaoffa).*gammaoffr+(gamma01+gammaoffa).*gammara).*va.* ...
  log(2)).^(-1)).*jtot0.*va.^(-1);

nr=gamma01.*na.*(gammaar.*gammaona.*T12+(gammaar+gammaoffa).* ...
  gammaonr.*T12+gammaar.*log(2)).*(gamma01.*(gammaoffr.* ...
  gammaona+(gammaona+gammaonr).*gammara).*T12+((gamma01+ ...
  gammaar+gammaoffa).*gammaoffr+(gamma01+gammaoffa).*gammara) ...
  .*log(2)).^(-1);

na0=gamma10.*na.*(gammaoffr.*gammaona.*T12+(gammaona+gammaonr).* ...
  gammara.*T12+(gammaoffr+gammara).*log(2)).*(gamma01.*( ...
  gammaoffr.*gammaona+(gammaona+gammaonr).*gammara).*T12+(( ...
  gamma01+gammaar+gammaoffa).*gammaoffr+(gamma01+gammaoffa).* ...
  gammara).*log(2)).^(-1);

nr0=gamma10.*na.*(gammaar.*gammaona.*T12+(gammaar+gammaoffa).* ...
  gammaonr.*T12+gammaar.*log(2)).*(gamma01.*(gammaoffr.* ...
  gammaona+(gammaona+gammaonr).*gammara).*T12+((gamma01+ ...
  gammaar+gammaoffa).*gammaoffr+(gamma01+gammaoffa).*gammara) ...
  .*log(2)).^(-1);

nfree=gamma10.*(gammaar.*gammaoffr+gammaoffa.*(gammaoffr+gammara)) ...
  .*na.*T12.*(gamma01.*(gammaoffr.*gammaona+(gammaona+ ...
  gammaonr).*gammara).*T12+((gamma01+gammaar+gammaoffa).* ...
  gammaoffr+(gamma01+gammaoffa).*gammara).*log(2)).^(-1);

ntot=na+nr+na0+nr0+nfree;
jtot=va*na;
janter=va*na;
nmt=na+nr+na0+nr0;
nmtproc=100*nmt./ntot;

close all

figure(1)
 hold on
    h1=plot(xanal,na);
    h2=plot(xanal,nr);
    h3=plot(xanal,na0);
    h4=plot(xanal,nr0);
    h5=plot(xanal,nfree);
   h=legend('n^{*\{0\}}_a','n^{*\{0\}}_r','n^{*\{0\}}_{a0}',...
       'n^{*\{0\}}_{r0}','n^{*\{0\}}_{free}',...
       'Location','northwest');
    set(h1,'linewidth',1.5,'color','r','LineStyle','-');
    set(h2,'linewidth',1.5,'color','b','LineStyle','--');
    set(h3,'linewidth',1.5,'color','g','LineStyle','-.');
    set(h4,'linewidth',1.5,'color','c','LineStyle','-');
    set(h5,'linewidth',1.5,'color','m','LineStyle','--');
    ylabel('n^{*\{0\}}, v^*_r=0, D^*_{free}=0','FontSize',12)
    h20=xlabel('x, \mum','FontSize',12);
 %   axis([0 L 0 1]);
 hold off
 
  figure(2)
 hold on
 h2=plot(xanal,janter);
 
 hold off
 set(h2,'linewidth',1.5,'color','r','LineStyle','-');
 
 h=legend('j^*_{anter}', 'Location','south');
 ylabel('j^*_{anter}','FontSize',12)
    h20=xlabel('x, \mum','FontSize',12);
    ylim([0.05 0.06])
 
 figure(3)
 hold on
h7=plot(xanal,ntot);
x8={0.001 50 100 150 200 250 300 350 400 450 500 550 600};
y8={TD*18.5*300/55 TD*13*300/55 TD*7*300/55 TD*7*300/55 ...
    TD*7*300/55 TD*7*300/55 ...
    TD*7*300/55 TD*7*300/55 TD*7*300/55 TD*7*300/55 ...
    TD*7*300/55 TD*17*300/55 TD*25*300/55};
xfull=xanal; nafull=na; nrfull=nr; na0full=na0; nr0full=nr0;
nfreefull=nfree; ntotfull=ntot; nmtprocfull=nmtproc;jtotfull=jtot;
vfull=jtot./ntot;
save('fullmoodel.mat','xfull','nafull','nrfull','na0full','nr0full',...
    'nfreefull','ntotfull','nmtprocfull',...
    'jtotfull','vfull');

% plot([x8{:}], [y8{:}],'o','color','k')
plot(xx12, yy12,'o','color','k')


   h=legend('model with kinesin-only transport'...
       ,'data from Fig. 3D of Black et al.','Location','northeast');
    set(h7,'linewidth',1.5,'color','r','LineStyle','-');
    ylabel('n^{*\{0\}}_{tot}, v^*_r=0, D^*_{free}=0','FontSize',12)
    h20=xlabel('x, \mum','FontSize',12);
  %   axis([0 L 0 150]);
 hold off
 
function dydx = f1(x,y,region)
  
    dydx = zeros(2,1);
    
    % y(1)=na
    % y(2)=nr
    
 dydx(1)=(gamma01.*(gammaoffr.*gammaona+gammaoffa.*gammaonr+gamma01.* ...
  (gammaona+gammaonr)+gammaar.*(gammaona+gammaonr)+(gammaona+ ...
  gammaonr).*gammara).*T12.*va+((gamma01+gammaar+gammaoffa).*( ...
  gamma01+gammaoffr)+(gamma01+gammaoffa).*gammara).*va.*log(2) ...
  ).^(-1).*(gamma01.*gamma10.*y(2).*(gammaoffr.*gammaona.*T12+( ...
  gammaona+gammaonr).*gammara.*T12+gammara.*log(2))+gamma10.* ...
  y(1).*((-1).*gamma01.*(gammaar.*gammaona+(gammaar+gammaoffa).* ...
  gammaonr).*T12+(-1).*((gammaar+gammaoffa).*(gamma01+ ...
  gammaoffr)+gammaoffa.*gammara).*log(2)));

dydx(2)=(gamma01.*(gammaoffr.*gammaona+gammaoffa.*gammaonr+gamma01.* ...
  (gammaona+gammaonr)+gammaar.*(gammaona+gammaonr)+(gammaona+ ...
  gammaonr).*gammara).*T12.*vr+((gamma01+gammaar+gammaoffa).*( ...
  gamma01+gammaoffr)+(gamma01+gammaoffa).*gammara).*vr.*log(2) ...
  ).^(-1).*((-1).*gamma01.*gamma10.*y(1).*(gammaar.*gammaona.* ...
  T12+(gammaar+gammaoffa).*gammaonr.*T12+gammaar.*log(2))+ ...
  gamma10.*y(2).*(gamma01.*(gammaoffr.*gammaona+(gammaona+ ...
  gammaonr).*gammara).*T12+((gamma01+gammaar+gammaoffa).* ...
  gammaoffr+(gamma01+gammaoffa).*gammara).*log(2)));
    
  end
  % -----------------------------------------------------------------------
    
      function res = bc1(YL,YR,jtot0,nfree0,A5)
  % Boundary (and internal) conditions
  % y(1)=na
    % y(2)=nr
     res = [va*YL(1)-vr*YL(2)-jtot0 % x = 0
 va*YR(1)-vr*YR(2)-A5*(1-exp(-log(2)/(T12*gammaar)))*va*YR(1)];
%       YL(3)-nfree0           % x = 0           
%       YR(4)  % x = L
 
   
       
  end

end