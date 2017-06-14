
function SuperResolution(img, output_path, depth_label, scale, sigma, nD, opts)

% CONTACT:
% Shuo Zhang (shuo.zhang@buaa.edu.cn)

% TERMS OF USE : 
% the spatial resolution of light field image is first downsampling and
% then reconstructed to the original size.

%% image reading and parameter setting

% img = double(imread(strcat(input_path,'lf.png')));             % the 4D light field image
% run (strcat(input_path,'depth_opt.m'));                        % parameter of the light field image

Dmin = opts.Dmin;                                              % the minimum disparity between the border view and the central view
Dmax = opts.Dmax;                                              % the maximum disparity between the border view and the central view
NumView = opts.NumView;                                        % the angular resolution

[h, w, nB] = size(img);
height = h / NumView;                                          % the height of the sub-aperture image
width = w / NumView;                                           % the height of the sub-aperture image
midView = round(NumView/2);                                    % the location of the central sub-aperture image
view_RGB = img(midView:NumView:end, midView:NumView:end,:);    % the central sub-aperture image

imwrite(uint8(view_RGB), strcat(output_path, num2str(scale),'X_original_downsampling_image.bmp'));

%% the preparation for sub-pixel location 
u = NumView:-1:1;
Dstep = (Dmax - Dmin) * 1.0 / (nD - 1);                        % the disparity difference of adjacent depth label 

shift = zeros(nD,NumView);
down = zeros(nD,NumView);
up = zeros(nD,NumView);
weight = zeros(nD,NumView);
weight1 = zeros(NumView,NumView,nD);
weight2 = zeros(NumView,NumView,nD);
weight3 = zeros(NumView,NumView,nD);
weight4 = zeros(NumView,NumView,nD);  

for depth = 1:nD
    shift(depth,:) =( (Dmin+(depth-1)*Dstep)-(Dmin+(depth-1)*Dstep)*2/((NumView-1)*1.0)*(u-1) )* scale;
%     down(depth,:) = sign(shift(depth,:)).*floor(abs(shift(depth,:))); 
%     up(depth,:) = sign(shift(depth,:)).*ceil(abs(shift(depth,:))); 
%     weight(depth,:) = abs(shift(depth,:) - down(depth,:));
    
    down(depth,:) = floor(shift(depth,:) ); 
    up(depth,:) = ceil(shift(depth,:)  ); 
    weight(depth,:) = shift(depth,:)  - floor(shift(depth,:)) ;
    
    weight1(:,:,depth) =  mtimes(weight(depth,:)',weight(depth,:));
    weight2(:,:,depth) =  mtimes((1-weight(depth,:))',weight(depth,:));
    weight3(:,:,depth) =  mtimes(weight(depth,:)',(1-weight(depth,:)));
    weight4(:,:,depth) =  mtimes((1-weight(depth,:))',(1-weight(depth,:)));
end

%% super-resolved sub-aperture image calculation

confidence = ones(height, width);
count = zeros(height*scale, width*scale, nB);
super_resolved = zeros(height*scale, width*scale, nB);

reverseStr = ''  ;
for x = 1:height;  
    for y = 1:width   
  
        %% corresponding micro-lens image
        img_mic = img((x-1)*NumView+1:(x-1)*NumView+NumView,(y-1)*NumView+1:(y-1)*NumView+NumView,:);  
        
        %% mcro-lens-based consistency metric (MCM) in Eq.3
        central_point = view_RGB(x,y,:);
        diff = (img_mic - repmat(central_point,[NumView, NumView])).^2;
        reliable = exp(-sum(diff,3)/(sigma^2));
        
        %% the calculated depth information
        depth = depth_label(x,y);
        
        %% the pixel location computation in view images.
        location_x_down = x*scale+down(depth,:);
        location_x_up   = x*scale+up(depth,:);
        location_y_down = y*scale+down(depth,:);
        location_y_up   = y*scale+up(depth,:);
        region = double((location_x_down>0)&(location_x_up<=(height*scale)))'*double((location_y_down>0)&(location_y_up<=(width*scale)));
        
        value1 = times(repmat(confidence(x,y)*reliable.*weight4(:,:,depth),[1,1,nB]), img_mic(:,:,:));
        value2 = times(repmat(confidence(x,y)*reliable.*weight3(:,:,depth),[1,1,nB]), img_mic(:,:,:));
        value3 = times(repmat(confidence(x,y)*reliable.*weight2(:,:,depth),[1,1,nB]), img_mic(:,:,:));
        value4 = times(repmat(confidence(x,y)*reliable.*weight1(:,:,depth),[1,1,nB]), img_mic(:,:,:));
        
        for i = 1:NumView
            for j = 1:NumView
                if region(i,j)==0 
                    continue;
                end
                super_resolved(uint16(location_x_down(i)),uint16(location_y_down(j)), :) = ...
                    super_resolved(uint16(location_x_down(i)),uint16(location_y_down(j)), :) + value1(i,j,:);
                super_resolved(uint16(location_x_up(i)),uint16(location_y_down(j)), :)   = ....
                    super_resolved(uint16(location_x_up(i)),uint16(location_y_down(j)), :) + value2(i,j,:);
                super_resolved(uint16(location_x_down(i)),uint16(location_y_up(j)), :)   = ...
                    super_resolved(uint16(location_x_down(i)),uint16(location_y_up(j)), :) + value3(i,j,:);
                super_resolved(uint16(location_x_up(i)),uint16(location_y_up(j)), :)     = ...
                    super_resolved(uint16(location_x_up(i)),uint16(location_y_up(j)), :)+ value4(i,j,:);
                count(uint16(location_x_down(i)),uint16(location_y_down(j)), :)    = ...
                    count(uint16(location_x_down(i)),uint16(location_y_down(j)), :) + confidence(x,y)*reliable(i,j)*weight4(i,j,depth);
                count(uint16(location_x_up(i)),uint16(location_y_down(j)), :)      = ...
                    count(uint16(location_x_up(i)),uint16(location_y_down(j)), :) + confidence(x,y)*reliable(i,j)*weight3(i,j,depth);
                count(uint16(location_x_down(i)),uint16(location_y_up(j)), :)      = ...
                    count(uint16(location_x_down(i)),uint16(location_y_up(j)), :) + confidence(x,y)*reliable(i,j)*weight2(i,j,depth);
                count(uint16(location_x_up(i)),uint16(location_y_up(j)), :)        = ...
                    count(uint16(location_x_up(i)),uint16(location_y_up(j)), :)+ confidence(x,y)*reliable(i,j)*weight1(i,j,depth);
            end
        end
    end
    msg = sprintf('Processing: %d/%d done!\n', x, height )  ;
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
end

super_resolved = uint8(super_resolved./count);
imwrite(super_resolved, strcat(output_path, num2str(scale),'X_super_resolved_image.bmp'));

view_RGB_bicubic = imresize(view_RGB,scale/1,'bicubic');
imwrite(uint8(view_RGB_bicubic), strcat(output_path, num2str(scale),'X_bicubic_image.bmp'));

end
