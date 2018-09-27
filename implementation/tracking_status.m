function [Conf, param] = tracking_status(response, num_props , param)
    Peak = max(response(:));
    mu = mean2(response);
    sigma = std2(response);
    PSR = (Peak - mu) / sigma;
%     Fmin = min(response(:)); 
%     PSR = (Fmax - Fmin)^2/mean2(power(response - Fmin,2));
%     area = find(response(:) > (Fmax/2));
%     area = numel(area);
%     total = numel(response(:));   
    if(param.frames == 0)
       param.avg_PSR = PSR;
       param.avg_Peak = Peak;
    end        
    ratio_PSR = PSR/param.avg_PSR;
    ratio_Peak = Peak/param.avg_Peak;
    ratio_Props = 1/num_props;
    Conf = ratio_PSR + ratio_Peak + ratio_Props;
    
    param.total_PSR = param.total_PSR + PSR;
    param.total_Peak = param.total_Peak + Peak;
    param.frames = param.frames + 1;
    param.avg_PSR = param.total_PSR / param.frames;
    param.avg_Peak = param.total_Peak / param.frames; 
    
    param.criteria = [ratio_PSR, ratio_Peak, ratio_Props, param.frames];
end