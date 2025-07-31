function out = make_H_v2(my_diag_elements, my_sigma_squared, other)

    
    H = diag(my_diag_elements);
    
    for i=1:(numel(my_diag_elements)-1)
        H(i, i+1) = my_sigma_squared;
        H(i+1, i) = my_sigma_squared;
    end
    
    out = H;
    
    
    