function bins_level = gray2bins(grayscale, bins_num)
%----------------------------------------------------------------------
% Quantify grayscale Image to several quantification level
% 	bins_level = gray2bins(grayscale, bins_num)
%
% Input Parameters List:
% grayscale:    Original Grayscale Image
% bins_num:     The number of quantification level (default = 16)
% 
% Output Parameters:
% bins_level:   The Image after quantification
%----------------------------------------------------------------------
    if bins_num == 0
        bins_num = 16;
    end
    bins_num = 256/bins_num;
    bins_level = (double(grayscale) + 1)/bins_num;
    bins_level = ceil(bins_level);
end