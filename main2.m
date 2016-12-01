function main2(mesh_num)
%% Parameters
clc
clearvars -except mesh_num
close all

% mesh_num = 1;
CFL = 100;
limit = 10e-7;
load(['1mesh' num2str(mesh_num) '.mat'])
U1 = U;

% U(:,1) = Uinf(1);
% U(:,2) = Uinf(2);
% U(:,3) = Uinf(3);
% U(:,4) = Uinf(4);

%% Import path 
addpath('givens/flux')
str = ['givens/meshes/' num2str(mesh_num)];
addpath(str)

%% Area and Time and Update U
[row col] = size(E);
for i = 1:row
    v1 = E(i,1);
    v2 = E(i,2);
    v3 = E(i,3);
    x1 = V(v1,1);
    y1 = V(v1,2);
    x2 = V(v2,1);
    y2 = V(v2,2);
    x3 = V(v3,1);
    y3 = V(v3,2);
    X = [x1 x2 x3];
    Y = [y1 y2 y3];
    A(i,:) = polyarea(X,Y);
end

%% Solver
Linf = [];
clmainarr = [];
clslatarr = [];
clflaparr = [];
cltotarr = [];
R(:) = 1;
j = 1;
while max(max(abs(R))) > limit;
    % for k = 1:3
    R(:) = 0;
    dT(:) = 0;
    dT_den(:) = 0;
    
    %% R update from IE and BE & Time and Update U
    [R dT_den clmain clslat clflap] = residual2(U,R,IE,BE,E,V,dT_den,Uinf,gam,A);
    cltot = clmain + clslat + clflap;
    dT = 2*A*CFL./dT_den;
    UFE = U - dT./A.*R;
    RFE = residual2(UFE,R,IE,BE,E,V,dT_den,Uinf,gam,A);
    U = 1/2*(U + UFE - dT./A.*RFE);
    
    %% Plot
    Linf = vertcat(Linf,max(abs(R)));
    Linf1 = Linf(:,1);
    Linf2 = Linf(:,2);
    Linf3 = Linf(:,3);
    Linf4 = Linf(:,4);
    clmainarr = [clmainarr clmain];
    clslatarr = [clslatarr clslat];
    clflaparr = [clflaparr clflap];
    cltotarr = [cltotarr cltot];
    if j == 1
        figure(2)
        subplot(2,1,1)
        h1 = semilogy(Linf1,'r');
        hold on
        h2 = semilogy(Linf2,'b');
        h3 = semilogy(Linf3,'g');
        h4 = semilogy(Linf4,'y');
        xlabel('Iteration')
        ylabel('Residual')
        legend('\rho','\rhou','\rhov','\rhoE')
        j = j + 1;
        subplot(2,1,2)
        h5 = plot(clmainarr,'r');
        hold on
        h6 = plot(clslatarr,'b');
        h7 = plot(clflaparr,'g');
        h8 = plot(cltotarr,'y');
        xlabel('Iteration')
        ylabel('c_l')
        legend('c_lmain','c_lslat','c_lflap','c_ltotal')
    elseif mod(j,10) == 0;
        h1.YData = Linf1;
        h2.YData = Linf2;
        h3.YData = Linf3;
        h4.YData = Linf4;
        h5.YData = clmainarr;
        h6.YData = clslatarr;
        h7.YData = clflaparr;
        h8.YData = cltotarr;
        drawnow
    end
    j = j+1;
end
subplot(2,1,1)
title(['2nd Order, Convergence Iteration = ' num2str(j)])
save(['2mesh' num2str(mesh_num) '.mat'])





