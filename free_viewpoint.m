function [output_image]  = free_viewpoint(image1, image2, varargin)
% This function generates an image from a virtual viewpoint between two
% real images. The output image has the same size aiss the input images.
%% Input parser
P = inputParser;
% choose the image set
P.addOptional('choose_img', false, @(x) islogical(x) );
% direct load the generated disparity map to save time
P.addOptional('load_disparityMap', false, @(x) islogical(x) );
% optimize the disparity map per segmentation and plane fitting
P.addOptional('do_optimization', false, @(x) islogical(x) );
% Fenstergroesse
P.addOptional('p', 0.5, @(x) isnumeric(x) && x>=0 && x<=1 );
% Unterer Schwellwert fuer die staerke der Korrelation zweier Merkmale
P.addOptional('down_ratio', 0.3, @(x) isnumeric(x) && x>0 && x<=1 );
% Plot oder nicht
P.addOptional('disparity_range', [-426,450], @(x) isnumeric(x) && x(1)<x(2));
P.addOptional('Np', 1400, @(x) isnumeric(x));
% den Input lesen
P.parse(varargin{:});
choose_img = P.Results.choose_img;
load_disparityMap = P.Results.load_disparityMap;
do_optimization = P.Results.do_optimization;    
p = P.Results.p;
down_ratio = P.Results.down_ratio;
disparity_range = P.Results.disparity_range; 
Np= P.Results.Np;     

    %% ein RGB-Bild in ein Graustufenbild umwandeln
    I1_gray = double(rgb_to_gray(image1));
    I2_gray = double(rgb_to_gray(image2));

    %% Harris-Detektor implementieren, um die Merkmalspunkte zu extrahieren
    Merkmale1 = harris_detektor(I1_gray, ...
        'segment_length', 15, 'k', 0.05, 'tau', 100000, ...
        'tile_size', [200,200], 'N', 10, 'min_dist', 20, 'do_plot', false);
    Merkmale2 = harris_detektor(I2_gray,  ...
        'segment_length', 15, 'k', 0.05, 'tau', 100000, ...
        'tile_size', [200,200], 'N', 10, 'min_dist', 20, 'do_plot', false);

    %% die Korrespondenzen nach NCC finden
    Korrespondenzen = korrespondenzen(I1_gray, I2_gray, Merkmale1, Merkmale2, ...
        'window_length', 25, 'min_corr', 0.95, 'do_plot', false);

    %% RanSaC Algorithmus implementieren, um die robuste Korrespondenzen zu erhalten
    [Korrespondenzen_robust,EF_robust] = ransac(I1_gray, I2_gray, Korrespondenzen, ...
        'epsilon', 0.7, 'p', 0.8, 'tolerance', 0.01, 'k', 8, 'do_plot', false);

    %% Achtpunktalgorithmus implementieren, um die Essentielen Matrix / Fundamentalmatrix zu schaetzen
    %load('K.mat');
    %achtpunktalgorithmus(Korrespondenzen_robust,K);

    %% find two linear transformation of the projective coordinates to map the epipole to infinity in the x-axis direction
    % [t1, t2] = estimateUncalibratedRectification(EF_robust, [Korrespondenzen_robust(2,:)',Korrespondenzen_robust(1,:)'], [Korrespondenzen_robust(4,:)',Korrespondenzen_robust(3,:)'],size(image2_gray));
    [t1, t2] = epipolar_rectification(EF_robust, [Korrespondenzen_robust(2,:)',Korrespondenzen_robust(1,:)'], [Korrespondenzen_robust(4,:)',Korrespondenzen_robust(3,:)'],size(I2_gray));

    %% apply the transformation to the entire image, so that all epipolar lines correspond to the horizontal scan lines
    % [I1_Rect_,I2_Rect_] = rectifyStereoImages(image1,image2,t1,t2,'OutputView','valid');
    [I1_gray_rect,I2_gray_rect] = image_rectification(I1_gray,I2_gray,t1,t2,'OutputView','valid','do_plot', false);
    [I1_rgb_rect,I2_rgb_rect]   = image_rectification(image1,image2,t1,t2,'OutputView','valid','do_plot', false);

    %% downsample the image    
    I1_gray_rect_down = imresize(I1_gray_rect,down_ratio);
    I2_gray_rect_down = imresize(I2_gray_rect,down_ratio);
    
if ~load_disparityMap 
    %% compute the disparity map using semi-global algorithm
    if p<=0.5
        disparity_map = disparity_computation(I1_gray_rect_down,I2_gray_rect_down, ...
            round([disparity_range(1),disparity_range(2)]*down_ratio), 'do_plot', false);
    else
        disparity_map = disparity_computation(I2_gray_rect_down,I1_gray_rect_down, ...
            round([-disparity_range(2),-disparity_range(1)]*down_ratio), 'do_plot', false);
        disparity_map = -disparity_map;
    end

%% or to save time, direct load the generated disparity map and rectified RGB image
else
    down_ratio = 0.5;
    if choose_img
        I1_rgb_rect = imread('img/L1_rect.png');
        I2_rgb_rect = imread('img/R1_rect.png');
        if p<=0.5
            load 'img/L1_250_310.mat';
        else
            load 'img/R1_310_250.mat';
            disparity_map = -disparity_map;
        end
    else
        I1_rgb_rect = imread('img/L2_rect.png');
        I2_rgb_rect = imread('img/R2_rect.png');
        if p<=0.5
            load 'img/L2_213_225.mat';
        else
            load 'img/R2_225_213.mat';
            disparity_map = -disparity_map;
        end 
    end
end  


%% upsample the disparity map
disparity_map = imresize(disparity_map,1/down_ratio);

%% imshow the disparity map
figure;
imshow(disparity_map,down_ratio*disparity_range);
title('disparity map');
colorbar;

%% optimise the disparity map per segmentation and plane fitting
if do_optimization
    
    %% segment the lines using Hough argorithm    
    %% segement the clusters using regiongrowing algorithm
    [lines1, lines2] = line_segmentation(image1, image2, choose_img, 'do_plot', false);
    
    if p<=0.5  
        line  = lines2;
        [wall,book,box1,box2,table] = cluster_segmentation(image2, 'start_point', [1500,3000], 'reg_maxdist', 0.2, 'down_ratio', 0.2, ...
            'hsize', 15, 'sigma', 2.5, 'do_plot', true);    
        mesh = logical(line+wall);
        [~, mesh] = image_rectification(lines1,mesh,t1,t2,'OutputView','valid','do_plot', false);
    else
        line  =lines1;
        [wall,book,box1,box2,table] = cluster_segmentation(image1, 'start_point', [980,738], 'reg_maxdist', 0.2, 'down_ratio', 0.2, ...
            'hsize', 15, 'sigma', 2.5, 'do_plot', true);
        mesh = logical(line+wall);
        [mesh, ~] = image_rectification(mesh,lines2,t1,t2,'OutputView','valid','do_plot', false);
    end
    mesh = logical(mesh);
    figure('Name','mesh'),imshow(mesh);
    
    
%     load 'img/wall.mat';
%     load 'img/book.mat';

    %% optimize the desparity of the wall
    if choose_img
     wall_disparity = 180*down_ratio;
    else
     wall_disparity = 280*down_ratio;
    end
    
    
    [x,y] = meshgrid( 1:size(disparity_map,2) , 1:size(disparity_map,1) );
    [xq,yq] = meshgrid( 1: size(disparity_map,2)/(size(mesh,2)): size(disparity_map,2), ...
                        1: size(disparity_map,1)/(size(mesh,1)): size(disparity_map,1) );
    disparity_map = interp2(x,y,disparity_map,xq,yq);


    disparity_map(mesh)=-wall_disparity;
    
    %% optimise the disparity map by means of plane fitting using ransac algorithm
%     disparity_map = ransac_plane(disparity_map,wall,3,'do_plot', true);
%     disparity_map = ransac_plane(disparity_map,book,3,'do_plot', true);
%     disparity_map = ransac_plane(disparity_map,box1,3,'do_plot', false);
%     disparity_map = ransac_plane(disparity_map,box2,3,'do_plot', false);
end

%% imshow the disparity map
figure;
imshow(disparity_map,down_ratio*disparity_range);
title('disparity map');
colorbar;

%% 3D reconstruction
reconstruction_3D(I1_rgb_rect, I2_rgb_rect, disparity_map, 'do_plot', false);

%% generate the image in the new viewpoint using depth-image-based rendering algorithm
if p<=0.5
    output_image_ = DIBR(I1_rgb_rect, disparity_map, p, Np,down_ratio,'do_plot', false);
else
    output_image_ = DIBR(I2_rgb_rect, disparity_map, p, Np,down_ratio,'do_plot', false);
end
    
%% upsample the image to the original size
output_image_ = double(output_image_);
[x,y] = meshgrid( 1:size(output_image_,2) , 1:size(output_image_,1) );
[xq,yq] = meshgrid( 1: size(output_image_,2)/(size(image1,2)+1): size(output_image_,2), ...
                    1: size(output_image_,1)/(size(image1,1)+1): size(output_image_,1) );
output_image(:,:,1) = interp2(x,y,output_image_(:,:,1),xq,yq);
output_image(:,:,2) = interp2(x,y,output_image_(:,:,2),xq,yq);
output_image(:,:,3) = interp2(x,y,output_image_(:,:,3),xq,yq);
output_image = uint8(output_image);
figure,imshow(output_image);title(['Virtual View 2000*3000 with p equal to ',num2str(p)]);
end

