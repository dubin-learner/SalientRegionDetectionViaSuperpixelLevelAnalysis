function salMap = sp2pixel(img, spnum, superpixels, adjloop, salMap_sp)
%----------------------------------------------------------------------
% Calculate Saliency Map from Superpixel level to Pixel level
% 	salMap = sp2pixel(img, spnum, superpixels, adjloop, salMap_sp)
%
% Input Parameters List:
% img:          Original Grayscale Image
% spnum:        The number of superpixel
% superpixels:  The index of superpixels in the Size of Original Image
% adjloop:      Adjacent Matrix of Superpixels
% salMap_sp:    Saliency Map in Superpixel Level
% 
% Output Parameters:
% salMap:       Saliency Map in Pixel Level
%----------------------------------------------------------------------

%% Superpixel level to Pixel level
% Generate a color quantization table
    bins_num = 8;
    img_bin = gray2bins(img, bins_num);
    
    sp_bin_hist = zeros(spnum, bins_num);
    for sp_i = 1:spnum
        superpixel_gray = img_bin(superpixels == sp_i);
        sp_bin_hist(sp_i, :) = bins2hist(superpixel_gray, bins_num);
    end
    
    salMap = zeros(size(img));
    for sp_i = 1:spnum
        p_i = superpixels == sp_i;
        bin_p_i = img_bin(p_i);

        neigbors = find( adjloop(sp_i, :) == 1 );
        N = length(neigbors);
    
        % sp_sal should be the Saliency Map of A superpixel and its neigbors
        % (1:N) are the Saliency of neigbors, (N + 1) for itself.
        sp_sal = salMap_sp([neigbors, sp_i]);
    
        tmp1 = sp_bin_hist(neigbors, bin_p_i);
        tmp2 = sp_bin_hist(sp_i, bin_p_i);
    
        dividend = sp_sal(1:N)*tmp1 + tmp2*sp_sal(N + 1);
        divisor = sum(tmp1) + tmp2;
    
        salMap(p_i == 1) = dividend./divisor;
    end
end