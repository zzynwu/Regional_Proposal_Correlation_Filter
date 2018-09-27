function [row, col] = subPixelPeak(row_prop, col_prop, response)
    [M,N] = size(response);
    ind_l = row_prop + (col_prop - 2) .* M;
    ind_c = row_prop + (col_prop - 1) .* M;
    ind_r = row_prop + col_prop .* M;  
    left = response(ind_l);
    center = response(ind_c);
    right = response(ind_r);
    divisor = 2 * center - left - right;
    if(divisor == 0)
        delta_r = 0;
    else
        delta_r = 0.5 * (right - left) ./ divisor;
    end
    row = row_prop + delta_r;

    ind_l = (row_prop - 1) + (col_prop - 1) .* M;
    ind_c = row_prop + (col_prop - 1) .* M;
    ind_r = (row_prop + 1) + (col_prop - 1) .* M;   
    left = response(ind_l);
    center = response(ind_c);
    right = response(ind_r);
    divisor = 2 * center - left - right;
    if(divisor == 0)
        delta_c = 0;
    else
        delta_c = 0.5 * (right - left) ./ divisor;
    end      
    col = col_prop + delta_c;
end
