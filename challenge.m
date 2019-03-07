%% Computer Vision Challenge
close all; clear;
% Groupnumber:
Group_number = 15;

% Groupmembers:
% members = {'Max Mustermann', 'Johannes Daten'};
members = {};

% Email-Adress (from Moodle!):
% mail = {'ga99abc@tum.de', 'daten.hannes@tum.de'};
mail = {};

%% Load images
choose_L1R1 = false; % choose image set L1,R1 or L2,R2
if choose_L1R1
    I1 = imread('img/L1.jpg');
    I2 = imread('img/R1.jpg');
    dispaiy_range = [-500,620];
    Np = 700;
else
    I1 = imread('img/L2.jpg');
    I2 = imread('img/R2.jpg');
    dispaiy_range = [-426,450];
    Np = 1400 ;
end

%% Free Viewpoint Rendering
% start execution timer -> tic;
tic
% Free Viewpoint Rendering
output_img = free_viewpoint(I1, I2, 'choose_img', choose_L1R1, 'load_disparityMap',false, ...
    'do_optimization', true, 'p', 0.5, 'down_ratio', 0.3, 'disparity_range', dispaiy_range,'Np',Np);
% stop execution timer -> toc;
elapsed_time = toc

%% Display Output
% Display Virtual View


