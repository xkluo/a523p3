function [FrhoIE FrhoBE dlIE dlBE] = massflow(order,IE,BE,V,U,Uinf,A,E)

%% Order 1
if order == 1
    % IE
    [row col] = size(IE);
    FrhoIE = zeros(row,1);
    dlIE = zeros(row,1);
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
        FrhoIE(i) = F(1);
        dlIE(i) = dl;
    end
    % BE
    [row col] = size(BE);
    FrhoBE = zeros(row,1);
    dlBE = zeros(row,1);
    %%
%     for i = 1:row
%         n1 = BE(i,1);
%         n2 = BE(i,2);
%         eL = BE(i,3);
%         ind = BE(i,4);
%         fL = BE(i,5);
%         xA = V(n1,1);
%         yA = V(n1,2);
%         xB = V(n2,1);
%         yB = V(n2,2);
%         dl = sqrt((xA-xB)^2+(yA-yB)^2);
%         n = [(yB-yA) (xA-xB)]/dl;
%         if ind == 1
%             UL = U(eL,:);
%             UR = Uinf;
%             [F smag] = flux(UL, UR, n);
%             FrhoBE(i) = F(1);
%             dlBE(i) = dl;
%         else
%             FrhoBE(i) = 0;
%             dlBE(i) = dl;
%         end
%     end
end

%% Order 2

if order == 2
    %% Gradient
    GU = cell(size(U));
    for i = 1:numel(GU);
        GU{i} = [0 0];
    end
    % update from IE
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
        Uhat = 1/2*(UL+UR);
        for j = 1:numel(Uhat)
            GU{eL,j} = GU{eL,j} + Uhat(j)*n*dl;
            GU{eR,j} = GU{eR,j} - Uhat(j)*n*dl;
        end
    end
    % update from BE
    [row col] = size(BE);
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
            Uhat = 1/2*(UL+UR);
            for j = 1:numel(Uhat)
                GU{eL,j} = GU{eL,j} + Uhat(j)*n*dl;
            end
        else
            UL = U(eL,:);
            UR = UL;
            Uhat = 1/2*(UL+UR);
            for j = 1:numel(Uhat)
                GU{eL,j} = GU{eL,j} + Uhat(j)*n*dl;
            end
        end
    end
    % Area normalization
    for i = 1:numel(GU);
        GU{i} = GU{i}/A(i);
    end
    %% IE
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
        midpoint = 1/2*[xA+xB yA+yB];
        % eL
        nodes = E(eL,:);
        centroid = 0;
        for j = 1:numel(nodes)
            centroid = centroid + V(nodes(j),:);
        end
        centroid = centroid/numel(nodes);
        dXL = midpoint -  centroid;
        % eR
        nodes = E(eR,:);
        centroid = 0;
        for j = 1:numel(nodes)
            centroid = centroid + V(nodes(j),:);
        end
        centroid = centroid/numel(nodes);
        dXR = midpoint -  centroid;
        % Obtain midpoint value
        UL = U(eL,:);
        UR = U(eR,:);
        for j = 1:numel(UL)
            UL(j) = UL(j) + sum(GU{eL,j}.*dXL);
            UR(j) = UR(j) + sum(GU{eR,j}.*dXR);
        end
        [F smag] = flux(UL, UR, n);
        FrhoIE(i) = F(1);
        dlIE(i) = dl;
    end
    % BE
    [row col] = size(BE);
    FrhoBE = zeros(row,1);
    dlBE = zeros(row,1);
end

