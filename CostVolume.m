function labels_min = CostVolume(img, output_path, nD, tau,  sigma, opts)

% CONTACT:
% Shuo Zhang (shuo.zhang@buaa.edu.cn)

%% image reading and parameter setting

Dmin = opts.Dmin;                                              % the minimum disparity between the border view and the central view
Dmax = opts.Dmax;                                              % the maximum disparity between the border view and the central view
NumView = opts.NumView;                                        % the angular resolution

[h, w, nB] = size(img);
height = h / NumView;                                          % the height of the sub-aperture image
width = w / NumView;                                           % the height of the sub-aperture image
midView = round(NumView/2);                                    % the location of the central sub-aperture image
view_RGB = img(midView:NumView:end, midView:NumView:end,:);    % the central sub-aperture image

%% the preparation for sub-pixel location 
u = NumView:-1:1;
Dstep = (Dmax - Dmin) * 1.0 / (nD - 1);                        % the disparity difference of adjacent depth label 

shift = zeros(nD,NumView);
down = zeros(nD,NumView);
up = zeros(nD,NumView);
weight = zeros(nD,NumView);
weight1 = zeros(NumView,NumView,nB,nD);
weight2 = zeros(NumView,NumView,nB,nD);
weight3 = zeros(NumView,NumView,nB,nD);
weight4 = zeros(NumView,NumView,nB,nD);  

for depth = 1:nD
    shift(depth,:) = (Dmin+(depth-1)*Dstep)-(Dmin+(depth-1)*Dstep)*2/((NumView-1)*1.0)*(u-1);
%     down(depth,:) = sign(shift(depth,:)).*floor(abs(shift(depth,:))); 
%     up(depth,:) = sign(shift(depth,:)).*ceil(abs(shift(depth,:))); 
%     weight(depth,:) = abs(shift(depth,:) - down(depth,:));
%     
    down(depth,:) = floor(shift(depth,:) ); 
    up(depth,:) = ceil(shift(depth,:)  ); 
    weight(depth,:) = shift(depth,:)  - floor(shift(depth,:)) ;
   
    weight1(:,:,:,depth) =  repmat(mtimes(weight(depth,:)',weight(depth,:)),[1,1,nB]);
    weight2(:,:,:,depth) =  repmat(mtimes((1-weight(depth,:))',weight(depth,:)),[1,1,nB]);
    weight3(:,:,:,depth) =  repmat(mtimes(weight(depth,:)',(1-weight(depth,:))),[1,1,nB]);
    weight4(:,:,:,depth) =  repmat(mtimes((1-weight(depth,:))',(1-weight(depth,:))),[1,1,nB]);
end

%% local cost volume calculate
cost = zeros(height, width, nD);
reverseStr = '';
tic;
for x = 1:height
    for y = 1:width
       %% corresponding micro-lens image
        img_mic = img((x-1)*NumView+1:(x-1)*NumView+NumView,(y-1)*NumView+1:(y-1)*NumView+NumView,:);  
        
       %% mcro-lens-based consistency metric (MCM) in Eq.3
        central_point = view_RGB(x,y,:);
        diff = (img_mic - repmat(central_point,[NumView, NumView])).^2;
        reliable = exp(-sum(diff,3)/(sigma^2));
        
        for depth = 1:nD 
           %% different regions extracted in central sub-aperture image
            location_x_down = x+down(depth,:);
            location_x_up   = x+up(depth,:);
            location_y_down = y+down(depth,:);
            location_y_up   = y+up(depth,:);
            x_l = (location_x_down>0)&(location_x_up<=(height));
            y_l = (location_y_down>0)&(location_y_up<=(width));
            match11 = view_RGB(uint16(location_x_down(x_l)), uint16(location_y_down(y_l)), :);
            match12 = view_RGB(uint16(location_x_up(x_l)),   uint16(location_y_down(y_l)), :);
            match13 = view_RGB(uint16(location_x_down(x_l)), uint16(location_y_up(y_l)),   :);
            match14 = view_RGB(uint16(location_x_up(x_l)),   uint16(location_y_up(y_l)),   :);
            match = times(weight1(x_l,y_l,:,depth),match14)+...
                times(weight2(x_l,y_l,:,depth),match13)+...
                times(weight3(x_l,y_l,:,depth),match12)+...
                times(weight4(x_l,y_l,:,depth),match11);
         
           %% micro-lens-based matching cost calculation in Eq.4
            diff_abs = sum(abs(match - (img_mic(x_l,y_l,:))),3);
            value_xsum = min(diff_abs, tau).*reliable(x_l,y_l);
            num = sum(sum(reliable(x_l,y_l)));
            cost(x,y,depth) = sum(value_xsum(:))./sum(num(:));
            
        end
    end
    msg = sprintf('Processing: %d/%d done!\n', x, height )  ;
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
end
toc;
%% local depth map calculation using winner-takes-all strategy 
[~,labels_min] = min(cost,[],3);
dis = uint8((256/(nD))*labels_min);
imwrite(dis,strcat(output_path,'initial_depth.bmp'));

%% guided filter for the initial cost volume
r = 9;
eps = 0.0001; 
cost_filter = zeros(height, width, nD);
for d=1:nD 
    p = cost(:,:,d);
    q = guidedfilter_color(double(view_RGB), double(p), r, eps);
    cost_filter(:,:,d) = q;
end

[~,labels_min] = min(cost_filter,[],3);
save_img = uint8((256/(nD))*labels_min);
imwrite(save_img,strcat(output_path,'depth_map.bmp'));

end


