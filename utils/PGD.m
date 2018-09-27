function hf_adj = PGD(model_xf, S_xx, yf, use_sz, small_filter_sz, params)
    S_xx = sum(conj(model_xf) .* model_xf, 3);
    S_xy = sum(bsxfun(@times, yf, conj(model_xf)),3);
    hf = S_xy ./ (S_xx + 1e-4);
    hf_prox = padding_fs(hf, use_sz, small_filter_sz); 
    Wf = hf_prox;
    Vf = hf_prox;
    i = 1;
    while(i <= params.admm_iterations)
        deltaW = S_xx .* conj(Wf) - S_xy;
        Vf_pre = Vf;
        Vf =  Wf - params.admm_lambda .* deltaW;  
        Vf_prox = padding_fs(Vf, use_sz, small_filter_sz);
        Wf = Vf_prox - (i-1)/(i+2).*(Vf_prox - Vf_pre);
        i = i+1;
    end
    hf_adj = hf; %Vf_prox;
end

function hf_prox = padding_fs(hf , use_sz, small_filter_sz) 
    h = ifft2(hf);
    [sx,sy,h] = get_subwindow_no_window(h, floor(use_sz/2) , small_filter_sz);
    t = single(zeros(use_sz(1), use_sz(2), size(h,3)));
    t(sx,sy,:) = h; 
    hf_prox = fft2(t); 
end
