function [img_out, opts] =  AngularResolutionChange(img_in, opts, newView)

% CONTACT:
% Shuo Zhang (shuo.zhang@buaa.edu.cn)

% TERMS OF USE : 
% the angular resolution of light field image is changed to newView 


[height, width, nB] = size(img_in);

shift = floor((opts.NumView - newView)/2);

img_new = zeros(height/(opts.NumView)*newView, width/(opts.NumView)*newView, nB);
for v = 1:newView
    for u = 1:newView
        img_new(u: newView:end, v: newView :end, :) = img_in(shift + u : opts.NumView:end, shift + v : opts.NumView:end, :);
    end
end

new_mid = floor(newView/2);
old_mid = floor(opts.NumView/2);
opts.Dmin = opts.Dmin/old_mid*new_mid;
opts.Dmax = opts.Dmax/old_mid*new_mid;
opts.NumView = newView;

img_out = img_new;


end