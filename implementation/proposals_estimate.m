function [row,col] = proposals_estimate(row_props, col_props, val_props, response_pwp, currentScaleFactor, area_resize_factor)
    [m,n] = size(response_pwp);
    K = numel(row_props);
    response = zeros(K,1); 
    val = zeros(K,1);
    row_resize = round(row_props * area_resize_factor * currentScaleFactor);
    col_resize = round(col_props * area_resize_factor * currentScaleFactor);
    row_final = row_resize + 1 + floor((m-1)/2);
    col_final = col_resize + 1 + floor((n-1)/2); 
    for i =1:K
        if( row_final(i) < m  &&  row_final(i) > 0  && col_final(i) < n  && col_final(i) > 0)
            val(i) = response_pwp(row_final(i), col_final(i));
            response(i) = val_props(i) * val(i);
        else
            val(i) = median(response_pwp(:));
            response(i) = val_props(i) * val(i);
        end
    end
    [~,maxID] = max(response);
    row = row_props(maxID);
    col = col_props(maxID);
end
    
    
    
    