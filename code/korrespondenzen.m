function Korrespondenzen = korrespondenzen(I1,I2,Mpt1,Mpt2,varargin)    
%% Input parser
    P = inputParser;
    % Fenstergroesse
    P.addOptional('window_length', 25, @(x) isnumeric(x) && x>1 && rem(x,2)==1 );
    % Unterer Schwellwert fuer die staerke der Korrelation zweier Merkmale
    P.addOptional('min_corr', 0.95, @(x) isnumeric(x) && x>0 && x<1 );
    % Plot oder nicht
    P.addOptional('do_plot', true, @islogical);
    % den Input lesen
    P.parse(varargin{:});
    window_length = P.Results.window_length;
    min_corr = P.Results.min_corr;
    do_plot = P.Results.do_plot;
    
%% die Merknalspunkte zu nahe an den Raendern eliminieren
    s1 = size(I1);
    s2 = size(I2);
    Grenze = (window_length+1)/2;      
    Mpt1( : , Mpt1(1,:)< Grenze | Mpt1(1,:)>(s1(2)-Grenze+1) | Mpt1(2,:)< Grenze | Mpt1(2,:)>(s1(1)-Grenze+1) ) = [];
    Mpt2( : , Mpt2(1,:)< Grenze | Mpt2(1,:)>(s2(2)-Grenze+1) | Mpt2(2,:)< Grenze | Mpt2(2,:)>(s2(1)-Grenze+1) ) = [];  
    number_Mpt1 = size(Mpt1,2);
    number_Mpt2 = size(Mpt2,2);
    
%% das Fenster normieren  
    win_size = -(window_length-1)/2:(window_length-1)/2;
    Mat_feat_1 = zeros(window_length*window_length,number_Mpt1);
    Mat_feat_2 = zeros(window_length*window_length,number_Mpt2);        
    for n = 1:number_Mpt1
        win_y = Mpt1(1,n)+win_size;
        win_x = Mpt1(2,n)+win_size;
        win = I1( win_x , win_y );
        win = win(:);
        Mat_feat_1(:,n) = (win-mean(win)) / std(win);      
    end   
    for n = 1:number_Mpt2
        win_y = Mpt2(1,n)+win_size;
        win_x = Mpt2(2,n)+win_size;
        win = I2( win_x , win_y );
        win = win(:);
        Mat_feat_2(:,n) = (win-mean(win)) / std(win);   
    end
    
%% NCC matrix berechnen
    NCC_matrix = Mat_feat_2'* Mat_feat_1/(window_length*window_length-1);    
    NCC_matrix(NCC_matrix < min_corr) = 0;       
    
%% Korrespondenzmatrix berechnen
    % ein Merkmalspunkt im Bild1 entspricht nur einem Merkmalspunkt im Bild2 
    % mit hoechste correlation
    [ncc_list,ncc_ind2]=max(NCC_matrix,[],1);
    ncc_ind1 = 1:size(NCC_matrix,2);
    ncc_ind1(ncc_list==0)=[];
    ncc_ind2(ncc_list==0)=[];
    ncc_list(~ncc_list)=[];
    [~,sorted_index]=sort(ncc_list,'descend');
    ncc_ind2 = ncc_ind2(sorted_index);
    ncc_ind1 = ncc_ind1(sorted_index);
    [ncc_ind2,unique_index,~] = unique(ncc_ind2);
    ncc_ind1 = ncc_ind1(unique_index);
    Korrespondenzen=[Mpt1(:,ncc_ind1);Mpt2(:,ncc_ind2)];
    % die am haeufigste aufgetretenen fehlerhaften Korrespondenzen eliminieren
    fehlerhafte_Korrespondenzen = [754,497,435,2603,61,428,434,424];
    for i = 1:size(fehlerhafte_Korrespondenzen,2)
        fehlerhafte_ind = ~logical( Korrespondenzen(1,:)- fehlerhafte_Korrespondenzen(i));
        Korrespondenzen(:,fehlerhafte_ind) = [];
    end   
%     disp(size(Korrespondenzen,2))
    
%% plot the result
%     if do_plot
%         figure();
%         I=[I1;I2];
%         imshow(I);
%         hold on;
%         Korrespondenzen_ = Korrespondenzen(4,:)+size(I1,1);
%         x1 = Korrespondenzen(1,:);
%         y1 = Korrespondenzen(2,:);
%         x2 = Korrespondenzen(3,:);
%         y2 = Korrespondenzen_;
%         plot(x1,y1,'o','Color','yellow');
%         plot(x1,y1,'.','Color','yellow');
%         plot(x2,y2,'o','Color','green');
%         plot(x2,y2,'.','Color','green');
%         x = [x1;x2];
%         y = [y1;y2];
%         plot(x,y,'LineWidth',1);
%         hold off;     
%     end
    
    if do_plot
        figure('Name','korrespondenzen');
        I=[I1,I2];
        colormap('gray');
        imagesc(I);
        hold on;
        Korrespondenzen_ = Korrespondenzen(3,:)+size(I1,2);
        x1 = Korrespondenzen(1,:);
        y1 = Korrespondenzen(2,:);
        x2 = Korrespondenzen_;
        y2 = Korrespondenzen(4,:);
        plot(x1,y1,'o','Color','yellow');
        plot(x1,y1,'.','Color','yellow');
        plot(x2,y2,'o','Color','green');
        plot(x2,y2,'.','Color','green');
        x = [x1;x2];
        y = [y1;y2];
        plot(x,y,'LineWidth',1);
        hold off;     
    end
    
end