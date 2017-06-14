function new_weight = WeightExpand(weight)
	[h, w] = size(weight);
	temp = ones(h+2,w+2);
	temp(2:end-1,2:end-1) = double(weight);

	temp = filter2(fspecial('average',3),temp);
	new_weight = temp>0.5;
	new_weight(2:end-1,2:end-1) = weight;
end