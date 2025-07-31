function out = make_K_v2(diffusion, drift, NSTATES)


    if NSTATES < 3
        error('nstates need to be >3')
    end    
    

    m = eye(NSTATES) * diffusion * -1;
    % this encodes for the drift to HIGHER states; if we wont want that,
    % we would set the drift to a NEGATIVE value.
    off_diag_upper = diffusion/2 - drift/2;
    off_diag_lower = diffusion/2 + drift/2;
    
    
    m(logical(triu(ones(size(m)),1) - triu(ones(size(m)),2))) = off_diag_upper;
    m(logical(triu(ones(size(m)),-1) - triu(ones(size(m))))) = off_diag_lower;
    m(2,1) = diffusion;
    m(end-1,end) = diffusion;
    
    
    
    out = m;
    
    
    