function [row, col, val, n] = get_proposals(response, param)
    Fmax = max(response(:));
    response = fftshift(response);
    % find local peak in large response area 
    ind = response < Fmax * param.thresh;
    response(ind) = 0;
    localPeak = imregionalmax(response);    
    proposals = response .* localPeak;
    % find  row/col of local peak and do clustering
    [row,col,val] = find(proposals(:,:));
    n = numel(val);
    if(n > 1)
        [~,I] = sort(val,'descend');     
        if( n > 10 ) 
            I = I(1:10);
            I = sort(I);
        end
        row = row(I);
        col = col(I);
        val = val(I);
        row_1 = circshift(row,-1);
        col_1 = circshift(col,-1);
        distant = (row - row_1).^2 + (col - col_1).^2;
        change = find(distant <= param.dist);
        change = change(1:end-1);
%         row(change) = round(row(change)/2 + row_1(change)/2);
%         col(change) = round(col(change)/2 + col_1(change)/2);
        row(mod(change,n)+1) = [];
        col(mod(change,n)+1) = [];
    end
end

    
    
    

