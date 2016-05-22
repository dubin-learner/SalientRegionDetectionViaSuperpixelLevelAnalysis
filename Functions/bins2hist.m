function bins_hist = bins2hist(superpixel, bins_num)
%----------------------------------------------------------------------
% Calculate the Histogram of a superpixel
% 	bins_hist = bins2hist(superpixel, bins_num)
%
% Input Parameters List:
% superpixel:   A superpixel Area
% bins_num:     The number of quantification level
% 
% Output Parameters:
% bins_hist:    The Histogram of a superpixel
%----------------------------------------------------------------------
    bins_hist = zeros(1, bins_num);
    for i1 = 1:bins_num
        bins_hist(i1) = length( find(superpixel == i1) );
    end
    total = length(superpixel);
    bins_hist = bins_hist/total;
end