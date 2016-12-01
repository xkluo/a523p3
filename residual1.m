function [R dT_den clmain clslat clflap] = residual1(U,R,IE,BE,V,dT_den,Uinf,gam)
% Note: R must be a zero matrix

%% R update from IE
[row col] = size(IE);
for i = 1:row
    n1 = IE(i,1);
    n2 = IE(i,2);
    eL = IE(i,3);
    eR = IE(i,4);
    fL = IE(i,5);
    fR = IE(i,6);
    xA = V(n1,1);
    yA = V(n1,2);
    xB = V(n2,1);
    yB = V(n2,2);
    dl = sqrt((xA-xB)^2+(yA-yB)^2);
    n = [(yB-yA) (xA-xB)]/dl;
    UL = U(eL,:);
    UR = U(eR,:);
    [F smag] = flux(UL, UR, n);
    R(eL,:) = R(eL,:) + F*dl;
    R(eR,:) = R(eR,:) - F*dl;
    dT_den(eL,:) = dT_den(eL,:) + smag*dl;
    dT_den(eR,:) = dT_den(eR,:) + smag*dl;
end

%% R update from BE
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
    if ind == 1
        UL = U(eL,:);
        UR = Uinf;
        [F smag] = flux(UL, UR, n);
        R(eL,:) = R(eL,:) + F*dl;
        dT_den(eL,:) = dT_den(eL,:) + smag*dl;
    else
        UL = U(eL,:);
        rho = UL(1);
        u = UL(2)/rho;
        v = UL(3)/rho;
        En = UL(4)/rho;
        vb = norm([u v]);
        pb = (gam-1)*(rho*En-1/2*rho*vb^2);
        F = [0 pb*n(1) pb*n(2) 0];
        R(eL,:) = R(eL,:) + F*dl;
    end
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
