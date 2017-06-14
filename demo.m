clear;clc;
dbstop if error

addpath(genpath('guidedfilter'));
addpath(genpath('preprocess'));


%% Parameters
nD = 64;         % the number of depth labels
tau = 25;        % tau in Eq.1
sigma = 100;     % sigma in Eq.3

input_path  = 'input/Buddha/';
output_path = 'output/'; 
mkdir(output_path);

%% image load
img = double(imread(strcat(input_path,'lf.bmp')));             % the 4D light field image
run (strcat(input_path,'depth_opt.m'));                        % parameter of the light field image

%% micro-lens matching cost and depth computation
depth_label = CostVolume(img, output_path, nD, tau, sigma, opts); 

%% view synthesis
add_angular_number = 1;                                        % the angular resolution of light field image increase 2*add_novel

tic;
ViewSynthesis(img, output_path, depth_label, add_angular_number, nD, opts);
toc;

%% super-resolution
scale = 2;                                                      % the scale of the super-resolved image

sgn_downsampling = 0;                                           % 1: the spatial resolution is first downsampling
if sgn_downsampling == 1
    [img, opts] = SpatialDownsampling(img, scale, opts);        
    depth_label = CostVolume(img, output_path, nD, tau, sigma, opts); 
end

tic;
SuperResolution(img, output_path, depth_label, scale, sigma, nD, opts);
toc;
