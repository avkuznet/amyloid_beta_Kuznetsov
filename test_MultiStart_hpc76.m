clear

global xx12 yy12 eps ntotvect

tStart = tic;

A = importdata('fig3Da.txt');
xx12=A(:,1); yy12=A(:,2);

close all

z(1)=1.602439e-02; z(2)=4.082819e-03; z(3)=6.673352e-17;
z(4)=1.479517e-10; z(5)=6.325496e-04; z(6)=5.055525e-12;
z(7)=6.673503e-17; z(8)=2.000000e+00; z(9)=5.613625e-02;
z(10)=7.161740e-01; z(11)=5.811685e-01;

para0=[z(1) z(2) z(3) z(4) z(5) z(6) z(7) z(8) z(9) z(10) ...
    z(11)];

nrandpoints=5000;
rpts = RandomStartPointSet('NumStartPoints',nrandpoints);
tpoints = CustomStartPointSet(para0);
allpts = {tpoints,rpts};

problem = createOptimProblem('fmincon',...
    'objective',@funerr5_6_hpc76,...
'x0',para0,'lb',[0 0 0 0 0 0 0 5 0 0 0],...
    'ub',[0.1 0.1 0.1 0.1 0.1 0.1 0.1 20 0.1 1 1],'options',...
    optimoptions(@fmincon,'Display','iter',...
    'TolX',1.0000e-10,'TolFun',1.0000e-100));

z(1)=1.602439e-02; z(2)=4.082819e-03; z(3)=6.673352e-17;
z(4)=1.479517e-10; z(5)=6.325496e-04; z(6)=5.055525e-12;
z(7)=6.673503e-17; z(8)=2.000000e+00; z(9)=5.613625e-02;
z(10)=7.161740e-01; z(11)=5.811685e-01;


ms = MultiStart;
[x,f,flag,outpt,allmins] = run(ms,problem,allpts)

save('full76_5000.mat','x')
format long

fileID = fopen('full76_5000.txt','w');
fprintf(fileID,'%d random points\n',nrandpoints);
fprintf(fileID,'%6s %6s %6s %6s %6s %6s %6s\n','par1',...
    'par2','par3','par4','par5','par6');
fprintf(fileID,'%8d %8d %8d %8d %8d %8d\n',x(1),x(2),x(3),x(4),...
    x(5),x(6));
fprintf(fileID,'%6s %6s %6s %6s %6s %6s %6s\n','par7',...
    'par8','par9','par10','par11');
fprintf(fileID,'%8d %8d %8d %8d %8d\n',x(7),x(8),x(9),x(10),x(11));
fprintf(fileID,'f_min= %8d\n',f);

tEnd = toc(tStart);
 fprintf(fileID,'%d minutes and %f seconds\n',floor(tEnd/60),rem(tEnd,60));
 tEndminutes=tEnd/60
fclose(fileID);
