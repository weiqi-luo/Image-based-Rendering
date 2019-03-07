function [Korrespondenzen_robust,largest_set_EF] = ransac(I1,I2,Korrespondenzen,varargin)
%% Input parser
    P = inputParser;
    % geschaetzte Wahrscheinlichkeit, dass ein zufaellig gewaehltes Korrespondenzpunktpaar ein Ausreisser ist
    P.addOptional('epsilon', 0.5, @(x) isnumeric(x) && (x>0) && (x<1));
    % gewuenschte Wahrscheinlichkeit, dass der Algorithmus einen Satz Korrespondenzpunktpaare liefert, 
    % in dem sich kein Ausreisser befindet
    P.addOptional('p', 0.5, @(x) isnumeric(x) && (x>0) && (x<1));
    % Toleranz, dass innerhalb dessen ein Korrespondenzpaar als zum Modell passend bewertet wird
    P.addOptional('tolerance', 0.01, @isnumeric);
    % Anzahl der benoetigten Puntke
    P.addOptional('k', 8, @isnumeric);
     % Plot oder nicht
    P.addOptional('do_plot', true, @islogical);
    % den Input lesen
    P.parse(varargin{:});
    epsilon = P.Results.epsilon;
    p = P.Results.p;
    tolerance = P.Results.tolerance;
    k = P.Results.k;
    do_plot = P.Results.do_plot;
    x1_pixel = [Korrespondenzen(1:2,:);ones(1,size(Korrespondenzen,2))];
    x2_pixel = [Korrespondenzen(3:4,:);ones(1,size(Korrespondenzen,2))];
    
%% Initialisation
    % die Anzahl der Korrespondenzen im bisher groessten Consensus-Set
    largest_set_size = 0;
    % die Sampson-Distanz des bisher groessten Consensus-Set 
    largest_set_dist = inf;
    % die Iterationszahl
    S = ceil( log(1-p) / log(1-(1-epsilon)^k) );
    e3_ = [0,-1,0;1,0,0;0,0,0];
     
%% ransac algorithm
    for s=1:S
        rand_index = randperm(size(Korrespondenzen,2),k);
        current_set_EF = achtpunktalgorithmus( Korrespondenzen(:,rand_index) );
        current_sd = sum(x2_pixel.*(current_set_EF*x1_pixel)).^2 ./ (sum((e3_*current_set_EF*x1_pixel).^2) + sum((e3_*current_set_EF'*x2_pixel).^2));
        curent_set_ind = current_sd<tolerance;
        current_set_size = sum( curent_set_ind );
        current_set_dist = sum( current_sd(curent_set_ind) );
        if current_set_size>largest_set_size || (current_set_size==largest_set_size && current_set_dist<largest_set_dist)
            largest_set_dist = current_set_dist;
            largest_set_size = current_set_size;
            largest_set_EF = current_set_EF;
            largest_set_ind = curent_set_ind;
        end        
    end
    Korrespondenzen_robust = Korrespondenzen(:,largest_set_ind);
    
%% plot the result
    if do_plot
        disp('die Iterationszahl S von RanSaC:')
        disp(S)
        figure('Name','ransac');
        I=[I1;I2];
        colormap('gray');
        imagesc(I);
        hold on;
        Korrespondenzen_ = Korrespondenzen_robust(4,:)+size(I1,1);
        x1 = Korrespondenzen_robust(1,:);
        y1 = Korrespondenzen_robust(2,:);
        x2 = Korrespondenzen_robust(3,:);
        y2 = Korrespondenzen_;
        plot(x1,y1,'o','Color','yellow');
        plot(x1,y1,'.','Color','yellow');
        plot(x2,y2,'o','Color','green');
        plot(x2,y2,'.','Color','green');       
        x = [Korrespondenzen_robust(1,:);Korrespondenzen_robust(3,:)];
        y = [Korrespondenzen_robust(2,:);Korrespondenzen_];
        plot(x,y,'LineWidth',1);
        
        t = cell(1,size(Korrespondenzen_robust,2));
        for n = 1:size(Korrespondenzen_robust,2)
            t{1,n}=n;
        end
        text(x1+20,y1+20,t,'Color','white');
        hold off;     
    end
end