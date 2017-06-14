function [new_lf, opts] = SpatialDownsampling(lf, scale, opts)
 
% CONTACT:
% Shuo Zhang (shuo.zhang@buaa.edu.cn)

% TERMS OF USE : 
% the spacial resolution of light field image is resized to scale X 


[h, w, nB] = size(lf);
height = h/opts.NumView;
width = w/opts.NumView;
angular_resolution = opts.NumView;


new_h = floor(height/scale)*angular_resolution;
new_w = floor(width/scale)*angular_resolution;
new_lf = zeros(new_h, new_w, nB);


for i = 1:angular_resolution
    for j = 1:angular_resolution
        view = lf(i:angular_resolution:new_h*scale, j:angular_resolution:new_w*scale, :);
        new_view = view(1:scale:end, 1:scale:end, :);
%            new_view = imresize(view, 1/scale);
        new_lf(i:angular_resolution:end, j:angular_resolution:end, :) = new_view;
    end
end

opts.Dmin = opts.Dmin/scale;
opts.Dmax = opts.Dmax/scale;
end