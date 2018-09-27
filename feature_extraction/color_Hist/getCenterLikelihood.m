function [center_likelihood, possible] = getCenterLikelihood(object_likelihood, m)
%GETCENTERLIKELIHOOD computes the sum over rectangles of size M.
% CENTER_LIKELIHOOD is the 'colour response'
    [h,w] = size(object_likelihood);
    n1 = h - m(1) + 1;
    n2 = w - m(2) + 1;

%% intergralImage
    SAT = integralImage(object_likelihood);
    i = 1:n1;
    j = 1:n2;
    center_likelihood = (SAT(i,j) + SAT(i+m(1), j+m(2)) - SAT(i+m(1), j) - SAT(i, j+m(2))) / prod(m);
    
%% Possible Object area
    Peak = max(center_likelihood(:));
    index = center_likelihood >= Peak * 0.7;
    p_area = sum(index(:));
    index = center_likelihood < Peak * 0.4;
    ip_area = sum(index(:));
    whole_area = n1 * n2;
    if( p_area > 0.6 * whole_area || ip_area > 0.6 * whole_area) %ip_area  0.8
        possible = 0;
    else
        possible = 1;
    end
end
