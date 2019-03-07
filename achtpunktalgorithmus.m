function EF = achtpunktalgorithmus(Korrespondenzen,k)
%% Vorbereitung
    x1 = [Korrespondenzen(1:2,:);ones(1,size(Korrespondenzen,2))];
    x2 = [Korrespondenzen(3:4,:);ones(1,size(Korrespondenzen,2))];
    if nargin>1
        x1 = K\x1;
        x2 = K\x2;
    end
    
%% Achtpunktalgorithmus    
    A=[ [ x1(1,:) ; x1(1,:) ; x1(1,:) ] .* x2(1:3,:) ; 
        [ x1(2,:) ; x1(2,:) ; x1(2,:) ] .* x2(1:3,:) ; 
        [ x1(3,:) ; x1(3,:) ; x1(3,:) ] .* x2(1:3,:) ]';
    [~,~,V] = svd(A);
    G = reshape(V(:,9),[3,3]);
    [Ug,Sg,Vg] = svd(G);
    % output Essentielen Matrix
    if nargin > 1
        sigma = diag([1,1,0]);
        EF = Ug*sigma*Vg';
    % output Fundamentalmatrix
    else
        Sg(3,3) = 0;
        EF = Ug*Sg*Vg';
    end
end