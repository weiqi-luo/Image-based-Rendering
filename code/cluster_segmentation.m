function [wall,book,box1,box2,table] = cluster_segmentation(I,varargin)
%% Input parser
P = inputParser;
P.addOptional('start_point', [10,10], @isnumeric);
P.addOptional('reg_maxdist', 0.2, @isnumeric);
P.addOptional('down_ratio', 0.2, @isnumeric);
P.addOptional('hsize', 15, @isnumeric);
P.addOptional('sigma', 2.5, @isnumeric);
P.addOptional('do_plot', true, @islogical);
P.parse(varargin{:});
reg_maxdist = P.Results.reg_maxdist;
down_ratio = P.Results.down_ratio;
table = round(P.Results.start_point*down_ratio);
hsize = P.Results.hsize;
sigma = P.Results.sigma;
do_plot = P.Results.do_plot;
I = rgb_to_gray(I);

%% detect the edge
G=fspecial('gaussian',hsize,sigma);
Img_smooth=conv2(I,G,'same');  % smooth image by Gaussiin convolution
[Ix,Iy]=gradient(Img_smooth);
f=Ix.^2+Iy.^2;
g = 1./(1+f);  % the traditional edge stop function
g = imresize(g,down_ratio);

%% segment the wall using regiongrowing algorithm
wall = false(size(g));
for ii = round([size(g,2)/50, size(g,2)/2, size(g,2)*95/100])
    j = regiongrowing(g, 5, ii, reg_maxdist);
    wall = j+wall;
end
wall = imresize(logical(wall),1/down_ratio); 

% imshow(wall)
% wall = 200*double(wall);
% wall = regiongrowing(wall, 10, 10,500);
% imshow(wall);

book = [122,172];
box1 = [223,378];
box2 = [216,337];
table = [233,87];
book = regiongrowing(g, book(1), book(2), reg_maxdist);
box1 = regiongrowing(g, box1(1), box1(2), reg_maxdist);
box2 = regiongrowing(g, box2(1), box2(2), reg_maxdist);
table = regiongrowing(g, table(1), table(2), reg_maxdist);
book = imresize(logical(book),1/down_ratio); 
box1 = imresize(logical(box1),1/down_ratio); 
box2 = imresize(logical(box2),1/down_ratio); 
table = imresize(logical(table),1/down_ratio); 

% hold on;
%     plot(table(2),table(1),'o');
% hold off;

%% plot
if do_plot
    figure(),imshow(g);
    figure('Name','wall'),imshow(wall);
%     figure('Name','book'),imshow(book);
%     figure('Name','box1'),imshow(box1);
%     figure('Name','box2'),imshow(box2);
%     figure('Name','table'),imshow(table);
end

% % G=fspecial('gaussian',1,2.5);
% % g=conv2(wall,G,'same');  % smooth image by Gaussiin convolution
% % g= imresize(g,0.1)
% % figure,imshow(g);
% % 
% % % g = imresize(,down_ratio);
% % [Ix,Iy]=gradient(g);
% % f=Ix.^2+Iy.^2;
% % g = 1./(1+f);  % the traditional edge stop function
% % 
% % % g = imresize(g,1/(down_ratio*0.1));
% % figure,imshow(g);


end

function J=regiongrowing(I,x,y,reg_maxdist)
% This function performs "region growing" in an image from a specified
% seedpoint (x,y)
%
% J = regiongrowing(I,x,y,t) 
% 
% I : input image 
% J : logical output image of region
% x,y : the position of the seedpoint (if not given uses function getpts)
% t : maximum intensity distance (defaults to 0.2)
%
% The region is iteratively grown by comparing all unallocated neighbouring pixels to the region. 
% The difference between a pixel's intensity value and the region's mean, 
% is used as a measure of similarity. The pixel with the smallest difference 
% measured this way is allocated to the respective region. 
% This process stops when the intensity difference between region mean and
% new pixel become larger than a certain treshold (t)
%
% Example:
%
% I = im2double(imread('medtest.png'));
% x=198; y=359;
% J = regiongrowing(I,x,y,0.2); 
% figure, imshow(I+J);
%
% Author: D. Kroon, University of Twente

if(exist('reg_maxdist','var')==0), reg_maxdist=0.2; end
if(exist('y','var')==0), figure, imshow(I,[]); [y,x]=getpts; y=round(y(1)); x=round(x(1)); end

J = zeros(size(I)); % Output 
Isizes = size(I); % Dimensions of input image

reg_mean = I(x,y); % The mean of the segmented region
reg_size = 1; % Number of pixels in region

% Free memory to store neighbours of the (segmented) region
neg_free = 10000; neg_pos=0;
neg_list = zeros(neg_free,3); 

pixdist=0; % Distance of the region newest pixel to the regio mean

% Neighbor locations (footprint)
neigb=[-1 0; 1 0; 0 -1;0 1];

% Start regiogrowing until distance between regio and posible new pixels become
% higher than a certain treshold
while(pixdist<reg_maxdist&&reg_size<numel(I))

    % Add new neighbors pixels
    for j=1:4
        % Calculate the neighbour coordinate
        xn = x +neigb(j,1); yn = y +neigb(j,2);
        
        % Check if neighbour is inside or outside the image
        ins=(xn>=1)&&(yn>=1)&&(xn<=Isizes(1))&&(yn<=Isizes(2));
        
        % Add neighbor if inside and not already part of the segmented area
        if(ins&&(J(xn,yn)==0)) 
                neg_pos = neg_pos+1;
                neg_list(neg_pos,:) = [xn yn I(xn,yn)]; J(xn,yn)=1;
        end
    end

    % Add a new block of free memory
    if(neg_pos+10>neg_free), neg_free=neg_free+10000; neg_list((neg_pos+1):neg_free,:)=0; end
    
    % Add pixel with intensity nearest to the mean of the region, to the region
    dist = abs(neg_list(1:neg_pos,3)-reg_mean);
    [pixdist, index] = min(dist);
    J(x,y)=2; reg_size=reg_size+1;
    
    % Calculate the new mean of the region
    reg_mean= (reg_mean*reg_size + neg_list(index,3))/(reg_size+1);
    
    % Save the x and y coordinates of the pixel (for the neighbour add proccess)
    x = neg_list(index,1); y = neg_list(index,2);
    
    % Remove the pixel from the neighbour (check) list
    neg_list(index,:)=neg_list(neg_pos,:); neg_pos=neg_pos-1;
end

% Return the segmented area as logical matrix
J=J>1;
end