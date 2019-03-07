function Merkmale = harris_detektor(input_image, varargin)
%% Input parser
    P = inputParser;
    % die Groesse des Bildsegments
    P.addOptional('segment_length', 15, @isnumeric);
    % gewichtet zwischen Ecken- und Kantenprioritaet
    P.addOptional('k', 0.05, @isnumeric);
    % der Schwellenwert zur Detektion einer Ecke
    P.addOptional('tau', 1000000, @isnumeric);
    % Kachelgroesse
    P.addOptional('tile_size', [200,200], @isnumeric);
    % maximale Anzahl an Merkmale innerhalb einer Kachel
    P.addOptional('N', 10, @isnumeric);
    % minimaler Pixelabstand zweier Merkmale
    P.addOptional('min_dist', 20, @isnumeric);
    % Plot oder nicht
    P.addOptional('do_plot', true, @islogical);
    % den Input lesen
    P.parse(varargin{:});
    segment_length = P.Results.segment_length;
    k = P.Results.k;
    tau = P.Results.tau;
    N = P.Results.N;
    min_dist = P.Results.min_dist;
    do_plot = P.Results.do_plot;
    tile_size = P.Results.tile_size;
    if numel(tile_size) == 1
        tile_size=[tile_size,tile_size];
    end
%     input_image = double(input_image);

%% sobel filter    
    wx = [1,0,-1;2,0,-2;1,0,-1];
    wy = [1,2,1;0,0,0;-1,-2,-1];
    Fx = conv2(input_image,wx,'same');
    Fy = conv2(input_image,wy,'same');
 
%% gewichted harris matrix
    w=fspecial('gaussian',[segment_length,1],segment_length/5);  
    G11 = conv2(w,w,Fx.*Fx,'same');
    G22 = conv2(w,w,Fy.*Fy,'same');
    G12 = conv2(w,w,Fx.*Fy,'same');

%% harris detector
    H = G11.*G22 - G12.^2 - k*((G11+G22).^2);   
    % die fehlerhafte Ausschlaege am Rand zu 0 setzen  
    size_image = size(input_image);    
    c = ceil(segment_length/2);
    mesh = zeros(size_image);   
    mesh(c:size_image(1)-c,c:size_image(2)-c) = 1;
    corners = mesh .* H;
    % die Messwerte unterhalb des Schwellwerts tau zu 0 setzen
    corners(corners<tau) = 0;
    % einen Nullrand der Breite min_dist um die Matrix corner fuegen
    m1 = zeros(size_image(1),min_dist);
    m2 = zeros(min_dist,size_image(2)+2*min_dist);
    corners = [m1,corners,m1];
    corners = [m2;corners;m2];
    % alle Merkmale der staerke nach absteigend sortieren
    [sorted_list,sorted_index] = sort(corners(:),'descend');
    sorted_index(sorted_list==0) = [];
    % eine kreisfoermige matrix aufstellen
    [x,y] = meshgrid(-min_dist:min_dist,-min_dist:min_dist);
    Cake = sqrt(x.^2 + y.^2)> min_dist;
    % Vorbereitung
    sorted_number = numel(sorted_index);
    AKKA = zeros(ceil(size_image(1)/tile_size(1)),ceil(size_image(2)/tile_size(2)));
    Merkmale = zeros(2,min(numel(AKKA)*N,sorted_number));
    feature_count = 1;
    % die Distanz zwischen Merkmale kontrollieren
    for  current_point = 1:sorted_number
        current_index = sorted_index(current_point);
        if(corners(current_index)==0)
            continue;
        else
            [row,col]=ind2sub(size(corners),current_index);
        end
        AKKA_row = floor((row-min_dist-1) / tile_size(1))+1;
        AKKA_col = floor((col-min_dist-1) / tile_size(2))+1;
        AKKA(AKKA_row,AKKA_col) = AKKA(AKKA_row,AKKA_col)+1;
        corners(row-min_dist:row+min_dist,col-min_dist:col+min_dist) = corners(row-min_dist:row+min_dist,col-min_dist:col+min_dist).*Cake;
        % die Zahl der Merkmale innerhalb einer Kachel kontrollieren 
        if AKKA(AKKA_row,AKKA_col)==N
            corners( (AKKA_row-1)*tile_size(1)+min_dist+1 : min( size(corners,1) , AKKA_row*tile_size(1)+min_dist ), ...
                    (AKKA_col-1)*tile_size(2)+min_dist+1 : min( size(corners,2) , AKKA_col*tile_size(2)+min_dist ) ) = 0;   
        end
        % den Merkmalspunkt speichern
        Merkmale(:,feature_count)=[col-min_dist;row-min_dist];
        feature_count = feature_count+1;
    end
    Merkmale = Merkmale( : , 1:feature_count-1 );
     
%% plot the result
    if do_plot
        figure('Name','harris_detector');
        colormap('gray');
        imagesc(input_image);
        hold on;
        scatter(Merkmale(1,:),Merkmale(2,:),'r','+');
    end
end