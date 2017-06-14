
function ViewSynthesis(img, output_path, depth_label, add_angular_number, nD, opts)

% CONTACT:
% Shuo Zhang (shuo.zhang@buaa.edu.cn)

% TERMS OF USE : 
% the angular resolution of light field image is first decrease 2*add_angular_number and
% then reconstructed to the original size.

%% image reading and parameter setting
T1 = -1;                % T1 in Eq.7;
                        % -1: if the depth is small for close objects;
                        % 1: if the depth is large for close objects

T2 = 1000;            % T2 in Eq.8;


% the angular resolution of light field image is first decrease 2*add_angular_number
[img, opts ] = AngularResolutionChange (img, opts, (opts.NumView-2*add_angular_number));                     
imwrite(uint8(img), strcat(output_path,'X_original_downsampling_lf.bmp'));

numView = opts.NumView;                                        % the angular resolution

[h, w, nB] = size(img);
height = h / numView;                                          % the height of the sub-aperture image
width = w / numView;                                           % the height of the sub-aperture image
midView = round(numView/2);                                    % the location of the central sub-aperture image
view_RGB = img(midView:numView:end, midView:numView:end,:);    % the central sub-aperture image


%% the preparation for sub-pixel location 
novel_numView = numView + add_angular_number*2;                         % the new angular resolution
Dmin = opts.Dmin/(midView-1)*(midView-1+add_angular_number); 
Dmax = opts.Dmax/(midView-1)*(midView-1+add_angular_number);

u = novel_numView:-1:1;
midView = round(novel_numView/2);
Dstep=(Dmax - Dmin) * 1.0 / (nD - 1);

shift = zeros(nD,novel_numView);
down = zeros(nD,novel_numView);
up = zeros(nD,novel_numView);
weight = zeros(nD,novel_numView);
weight1 = zeros(novel_numView,novel_numView,nD);
weight2 = zeros(novel_numView,novel_numView,nD);
weight3 = zeros(novel_numView,novel_numView,nD);
weight4 = zeros(novel_numView,novel_numView,nD);

for depth = 1:nD
    shift(depth,:) = ((Dmin+(depth-1)*Dstep)-(Dmin+(depth-1)*Dstep)*2/((novel_numView-1)*1.0)*(u-1));
    down(depth,:) = sign(shift(depth,:)).*floor(abs(shift(depth,:)));
    up(depth,:) = sign(shift(depth,:)).*ceil(abs(shift(depth,:))); 
    weight(depth,:) = abs(shift(depth,:) - down(depth,:));
    weight1(:,:,depth) =  mtimes(weight(depth,:)',weight(depth,:));
    weight2(:,:,depth) =  mtimes((1-weight(depth,:))',weight(depth,:));
    weight3(:,:,depth) =  mtimes(weight(depth,:)',(1-weight(depth,:)));
    weight4(:,:,depth) =  mtimes((1-weight(depth,:))',(1-weight(depth,:)));
end

%% depth information calculate
reverseStr = ''  ;
img_novel = zeros(height*novel_numView, width*novel_numView, nB);
conf_novel = zeros(height*novel_numView, width*novel_numView);
depth_label = double(depth_label);
for x = 1:height; 
    for y = 1:width ;
        
        depth = depth_label(x,y); % depth of the reference point
        
        location_x_down = x+down(depth,:);
        location_x_up   = x+up(depth,:);
        location_y_down = y+down(depth,:);
        location_y_up   = y+up(depth,:);
        x_l = (location_x_up>0)&(location_x_up<=(height));
        y_l = (location_y_up>0)&(location_y_up<=(width));
        
       %% resized micro-lens image, M^{h+k} in Eq.9
        img_mic = img((x-1)*numView+1:(x-1)*numView+numView,(y-1)*numView+1:(y-1)*numView+numView,1:nB);
        img_mic_resize = imresize(img_mic, novel_numView/numView, 'nearest');   
        
       %% region extracted from central sub-aperture image, V^{h+k} in Eq.9
        match_new = zeros(novel_numView, novel_numView, nB);
        
        match11 = view_RGB(uint16(location_x_down(x_l)), uint16(location_y_down(y_l)), :);
        match12 = view_RGB(uint16(location_x_up(x_l)),   uint16(location_y_down(y_l)), :);
        match13 = view_RGB(uint16(location_x_down(x_l)), uint16(location_y_up(y_l)),   :);
        match14 = view_RGB(uint16(location_x_up(x_l)),   uint16(location_y_up(y_l)),   :);
        match = times(repmat(weight1(x_l,y_l,depth),[1,1,nB]),match14)+...
            times(repmat(weight2(x_l,y_l,depth),[1,1,nB]),match13)+...
            times(repmat(weight3(x_l,y_l,depth),[1,1,nB]),match12)+...
            times(repmat(weight4(x_l,y_l,depth),[1,1,nB]),match11);
        match_new(x_l,y_l,:) = match;
        
       %% depth of the region extracted from central sub-aperture image, D^{h+k} in Eq.9
        depth_new = zeros(novel_numView, novel_numView);
       
        match11 = depth_label(uint16(location_x_down(x_l)), uint16(location_y_down(y_l)), :);
        match12 = depth_label(uint16(location_x_up(x_l)),   uint16(location_y_down(y_l)), :);
        match13 = depth_label(uint16(location_x_down(x_l)), uint16(location_y_up(y_l)),   :);
        match14 = depth_label(uint16(location_x_up(x_l)),   uint16(location_y_up(y_l)),   :);
        match = times(weight1(x_l,y_l,depth),match14)+...
            times(weight2(x_l,y_l,depth),match13)+...
            times(weight3(x_l,y_l,depth),match12)+...
            times(weight4(x_l,y_l,depth),match11);
        depth_new(x_l,y_l,:) = match;
        
       %% Newly-appeared region, O_1^{h+k}       
        weight_d = (depth_new - depth)< T1; % Eq.7
                
       %% Newly-occluded region, O_2^{h}
        diff = (match_new(add_angular_number+1:end-add_angular_number,add_angular_number+1:end-add_angular_number,:) - img_mic).^2; 
        weight_s = sum(diff,3) < T2;        % Eq.8
        weight = WeightExpand(weight_s);    
    
        %% the expanded micro-lens image, M^{h+k} in Eq.9 
        weight = weight&(1-weight_d); % final occlusion estimation
               match_new = (1-repmat(weight,[1,1,nB])).*img_mic_resize + repmat(weight,[1,1,nB]).*match_new;

        conf_novel((x-1)*novel_numView+1:(x-1)*novel_numView+novel_numView,...
            (y-1)*novel_numView+1:(y-1)*novel_numView+novel_numView)=weight;

        img_novel((x-1)*novel_numView+1:(x-1)*novel_numView+novel_numView,...
            (y-1)*novel_numView+1:(y-1)*novel_numView+novel_numView,1:nB) = match_new;
        
    end
    msg = sprintf('Processing: %d/%d done!\n', x, height);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
end

img_novel = uint8(img_novel);
imwrite(img_novel, strcat(output_path, num2str(add_angular_number),'add_lf.bmp'));

view_left = img_novel(midView:novel_numView:end,1:novel_numView:end,:);        
imwrite(view_left, strcat(output_path, num2str(add_angular_number),'add_synthesized_left_view.bmp'));

view_up = img_novel(1:novel_numView:end,midView:novel_numView:end,:);        
imwrite(view_up, strcat(output_path, num2str(add_angular_number),'add_synthesized_up_view.bmp'));

view_left_upper = img_novel(1:novel_numView:end,1:novel_numView:end,:);        
imwrite(view_left_upper, strcat(output_path, num2str(add_angular_number),'add_synthesized_left_upper_view.bmp'));

end
