%% Parameter
clear all
clc
mesharr = 0:3;
order = 2; % change this for different orders (either 1 or 2)
if order == 1
    sty = 'flat';
elseif order == 2
    sty = 'interp';
end

%% Import path
addpath('givens/flux')

%%
for k = 1:numel(mesharr);
    mesh = mesharr(k);
    %     mesh = 3;
    
    clear var -except [order mesh];
    
    file = [num2str(order) 'mesh' num2str(mesh) '.mat'];
    load(file);
    
    %% Inf
    rhoEinf = Uinf(4);
    rhouinf = Uinf(2);
    rhovinf = Uinf(3);
    qinf = 1/2*sqrt(rhouinf^2+rhovinf^2);
    pinf = (gam-1)*(rhoEinf-qinf);
    
    %% cp
    figure(1)
    subplot(2,2,k)
    [row col] = size(BE);
    cp1 = 0;
    cp2 = 0;
    cp3 = 0;
    for i = 1:row
        n1 = BE(i,1);
        n2 = BE(i,2);
        eL = BE(i,3);
        ind = BE(i,4);
        fL = BE(i,5);
        xA = V(n1,1);
        yA = V(n1,2);
        xB = V(n2,1);
        yB = V(n2,2);
        dl = sqrt((xA-xB)^2+(yA-yB)^2);
        n = [(yB-yA) (xA-xB)]/dl;
        if ind == 2
            rhoE = U(eL,4);
            rhou = U(eL,2);
            rhov = U(eL,3);
            q = 1/2*sqrt(rhou^2+rhov^2);
            p = (gam-1)*(rhoE-q);
            cp = (p-pinf)/q;
            xavg = 1/2*(xA+xB);
            plot(xavg,cp,'xr')
            hold on
            if cp > cp1
                cp1 = cp;
            end
        elseif ind == 3
            rhoE = U(eL,4);
            rhou = U(eL,2);
            rhov = U(eL,3);
            q = 1/2*sqrt(rhou^2+rhov^2);
            p = (gam-1)*(rhoE-q);
            cp = (p-pinf)/q;
            xavg = 1/2*(xA+xB);
            plot(xavg,cp,'xb')
            hold on
            if cp > cp2
                cp2 = cp;
            end
        elseif ind == 4
            rhoE = U(eL,4);
            rhou = U(eL,2);
            rhov = U(eL,3);
            q = 1/2*sqrt(rhou^2+rhov^2);
            p = (gam-1)*(rhoE-q);
            cp = (p-pinf)/q;
            xavg = 1/2*(xA+xB);
            plot(xavg,cp,'xg')
            hold on
            if cp > cp3
                cp3 = cp;
            end
        end
    end
    maxcp = max([cp1 cp2 cp3]);
    set(gca,'Ydir','reverse')
    axis([-0.3, 0.9, -1, 4])
    xlabel('x')
    ylabel('c_p')
    title(['Order ' num2str(order) ', Mesh ' num2str(mesh) ', c_p plot, max c_p(main,slat,flap) = '...
        num2str(cp1) ',' num2str(cp2) ',' num2str(cp3)])
    
    %% Convergence
    figure(2)
    subplot(2,numel(mesharr),k)
    h1 = semilogy(1:10:numel(Linf1),Linf1(1:10:end),'r');
    hold on
    h2 = semilogy(1:10:numel(Linf2),Linf2(1:10:end),'b');
    h3 = semilogy(1:10:numel(Linf3),Linf3(1:10:end),'g');
    h4 = semilogy(1:10:numel(Linf4),Linf4(1:10:end),'y');
    xlabel('Iteration')
    ylabel('Residual')
    legend('\rho','\rhou','\rhov','\rhoE')
    j = j + 1;
    
    subplot(2,numel(mesharr),k+numel(mesharr))
    h5 = plot(1:10:numel(clmainarr),clmainarr(1:10:end),'r');
    hold on
    h6 = plot(1:10:numel(clslatarr),clslatarr(1:10:end),'b');
    h7 = plot(1:10:numel(clflaparr),clflaparr(1:10:end),'g');
    h8 = plot(1:10:numel(cltotarr),cltotarr(1:10:end),'y');
    xlabel('Iteration')
    ylabel('c_l')
    legend('c_lmain','c_lslat','c_lflap','c_ltotal')
    
    subplot(2,numel(mesharr),k)
    title(['Order ' num2str(order) ', Mesh ' num2str(mesh) ', Converged Iter. = ' num2str(j)])
    
    subplot(2,numel(mesharr),k+numel(mesharr))
    title(['Order ' num2str(order) ', Mesh ' num2str(mesh) ', Converged c_ptot = ' num2str(cltotarr(end))])
    
    %% Report
    [row col] = size(BE);
    Fmain = [0 0];
    Fslat = [0 0];
    Fflap = [0 0];
    for i = 1:row
        n1 = BE(i,1);
        n2 = BE(i,2);
        eL = BE(i,3);
        ind = BE(i,4);
        fL = BE(i,5);
        xA = V(n1,1);
        yA = V(n1,2);
        xB = V(n2,1);
        yB = V(n2,2);
        dl = sqrt((xA-xB)^2+(yA-yB)^2);
        n = [(yB-yA) (xA-xB)]/dl;
        if ind == 2
            rhoE = U(eL,4);
            rhou = U(eL,2);
            rhov = U(eL,3);
            p = (gam-1)*(rhoE-1/2*sqrt(rhou^2+rhov^2));
            Fmain = Fmain + p*n*dl;
        elseif ind == 3
            rhoE = U(eL,4);
            rhou = U(eL,2);
            rhov = U(eL,3);
            p = (gam-1)*(rhoE-1/2*sqrt(rhou^2+rhov^2));
            Fslat = Fslat + p*n*dl;
        elseif ind == 4
            rhoE = U(eL,4);
            rhou = U(eL,2);
            rhov = U(eL,3);
            p = (gam-1)*(rhoE-1/2*sqrt(rhou^2+rhov^2));
            Fflap = Fflap + p*n*dl;
        end
    end
    chord = 0.5588;
    den = 1/2*sqrt(Uinf(2)^2+Uinf(3)^2)*chord;
    clmain = Fmain(2)/den;
    clslat = Fslat(2)/den;
    clflap = Fflap(2)/den;
    Ftot = Fmain+Fslat+Fflap;
    cltot = Ftot(2)/den;
    cd = Ftot(1)/den;
    
    str = {};
    str{1} = ['Order ' num2str(order) ', Mesh ' num2str(mesh)];
    str{2} = ['cl_main = ' num2str(clmain)];
    str{3} = ['cl_slat = ' num2str(clslat)];
    str{4} = ['cl_flap = ' num2str(clflap)];
    str{5} = ['cl_total = ' num2str(cltot)];
    str{6} = ['cd = ' num2str(cd)];
    for j = 1:numel(str)
        disp(str{j})
    end
    
    %% Boundary Flow
    [FrhoIE FrhoBE dlIE dlBE] = massflow(order,IE,BE,V,U,Uinf,A,E);
    
    %% Streamline
    [row col] = size(V);
    psi = zeros(row,1)+NaN;
    bc = psi;
    % Set BE nodes(psi) to zero
    [row col] = size(BE);
    for i = 1:row
        n1 = BE(i,1);
        n2 = BE(i,2);
        ind = BE(i,4);
        if ind ==2
            psi(n1) = 0;
            psi(n2) = 0;
        end
        bc(n1) = ind;
        bc(n2) = ind;
    end
    % set IE psi
    while norm(double(isnan(psi))) > 0
        [row col] = size(IE);
        for i = 1:row
            n1 = IE(i,1);
            n2 = IE(i,2);
            if xor(isnan(psi(n1)), isnan(psi(n2)))
                if isnan(psi(n1))
                    psi(n1) = psi(n2) - FrhoIE(i)*dlIE(i);
                elseif isnan(psi(n2))
                    psi(n2) = psi(n1) + FrhoIE(i)*dlIE(i);
                end
            end
        end
    end
    %% Flow quantities
    rho = U(:,1);
    u = U(:,2)./rho;
    v = U(:,3)./rho;
    vmag = sqrt(u.^2 + v.^2);
    En = U(:,4)./rho;
    pres = (gam-1)*(rho.*En-1/2*rho.*vmag.^2);
    c = sqrt(gam*pres./rho);
    M = vmag./c;
    
    p = V;
    temp = edge(:,1)*0;
    e = [edge,temp,temp,temp,temp,temp];
    temp = E(:,1)*0;
    t = [E temp];
    
    %% Plot Streamline
    psi2 = psi;
    psi2(psi~=0) = NaN;
    
    figure(3);
    lvl = 50;
    lvl2 = 500;
    subplot(2,2,k);
    
    x = V(:,1);
    y = V(:,2);
    v = psi;
    xarr = linspace(-0.3, 0.9,lvl2);
    yarr = linspace(-0.35, 0.25,lvl2);
    [xq,yq] = meshgrid(xarr, yarr);
    vq = griddata(x,y,v,xq,yq);
    contour(xq,yq,vq,lvl)
    hold on
    
    edge = vertcat(BE(:,1:2));
    [row col] = size(edge);
    for i = 1:row
        n1 = edge(i,1);
        n2 = edge(i,2);
        X = [V(n1,1) V(n2,1)];
        Y = [V(n1,2) V(n2,2)];
        line(X,Y, 'Color', [0 0 0]);
    end
    
    xlabel('x')
    ylabel('y')
    title(['Order ' num2str(order) ', Mesh ' num2str(mesh) ', Streamline'])
    axis('equal')
    axis([-0.3, 0.9, -0.35, 0.25]);
    colorbar
    colormap('jet')
    caxis([-.1 .1])
    
    %% Streamline 2
    
    figure(7);
    lvl = 50;
    s = subplot(2,2,k);
    
    h = pdeplot(p',e',t', ...
        'XYData', psi, ...
        'ZData', psi, ...
        'Mesh', 'on', ...
        'ColorMap', 'jet', ...
        'XYStyle', 'interp');
    
    xlabel('x')
    ylabel('y')
    title(['Order ' num2str(order) ', Mesh ' num2str(mesh) ', Streamline 2'])
    axis('equal')
    axis([-0.3, 0.9, -0.35, 0.25]);
    colorbar
    colormap('jet')
    caxis([-.1 .1])
    
    
    %% Plot Pressure
    figure(4)
    subplot(2,2,k)
    
    h = pdeplot(p',e',t', ...
        'XYData', pres, ...
        'ZData', pres, ...
        'Mesh', 'on', ...
        'ColorMap', 'jet', ...
        'XYStyle', sty);
    
    xlabel('x')
    ylabel('y')
    title(['Order ' num2str(order) ', Mesh ' num2str(mesh) ', Pressure'])
    axis('equal')
    axis([-0.3, 0.9, -0.35, 0.25]);
    colorbar
    colormap('jet')
    
    %% Plot Mach Magnitude
    figure(5)
    subplot(2,2,k)
    
    h = pdeplot(p',e',t', ...
        'XYData', M, ...
        'ZData', M, ...
        'Mesh', 'on', ...
        'ColorMap', 'jet', ...
        'XYStyle', sty);
    
    xlabel('x')
    ylabel('y')
    title(['Order ' num2str(order) ', Mesh ' num2str(mesh) ', Mach'])
    axis('equal')
    axis([-0.3, 0.9, -0.35, 0.25]);
    colorbar
    colormap('jet')
    
    
    %% Plot Mesh
    figure(6);
    subplot(2,2,k)
    
    %     edge = vertcat(IE(:,1:2),BE(:,1:2));
    %     [row col] = size(edge);
    %     for i = 1:row
    %         n1 = edge(i,1);
    %         n2 = edge(i,2);
    %         X = [V(n1,1) V(n2,1)];
    %         Y = [V(n1,2) V(n2,2)];
    %         line(X,Y);
    %     end
    
    h = pdeplot(p',e',t');
    set(h,'Color',[0 0 1])
    
    axis('equal')
    axis([-0.3, 0.9, -0.35, 0.25]);
    xlabel('x')
    ylabel('y')
    title(['Mesh ' num2str(mesh) ', Mesh'])
    
    %% Flow rate
    psi2 = psi(find(bc==2,1)); % main
    psi3 = psi(find(bc==3,1)); % slat
    psi4 = psi(find(bc==4,1)); % flap
    mdot23 = abs(psi2-psi3);
    mdot34 = abs(psi3-psi4);
    str = ['mdot(slat,main) = ' num2str(mdot23)];
    disp(str)
    str = ['mdot(flap,main) = ' num2str(mdot34)];
    disp(str)       
    
    disp(' ')
end

%% Flow rate
% order 1
od = 1;
x =[0.01051 0.02817 0.4883 0.4908 0.01051 0.02873 0.4883 0.4929 0.01051 0.02853 0.4883 0.4909 0.01051 0.02678 0.4883 0.4924];
y = [0.00286 0.0005632 0.01712 0.01042 0.00286 0.001265 0.01712 0.01139 0.00286 0.001022 0.01712 0.01046 0.00286 -0.001466 0.01712 0.01121];
psi = [0.02082 0 0 -0.001617 0.002645 0 0 -0.001789 0.003379 0 0 -0.001896 0.003989 0 0 -0.001996];
for i = 1:4:numel(x)
    v1 = [x(i) y(i)];
    v2 = [x(i+1) y(i+1)];
    v3 = [x(i+2) y(i+2)];
    v4 = [x(i+3) y(i+3)];
    psi1 = psi(i);
    psi2 = psi(i+1);
    psi3 = psi(i+2);
    psi3 = psi(i+3);
    d12 = norm(v1-v2);
    d34 = norm(v3-v4);
end



%%
f1 = figure(1);
f2 = figure(2);
f3 = figure(3);
f4 = figure(4);
f5 = figure(5);
f6 = figure(6);
f7 = figure(7);
% 
% close(f1)
% close(f2)
% close(f3)
% close(f4)
% close(f5)
% close(f6)
% close(f7)



