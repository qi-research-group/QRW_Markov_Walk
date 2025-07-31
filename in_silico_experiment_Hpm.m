% figure out some basics:

% some code regarding gaussians etc.
% I am guessing we just use that...

do_these = {...
    [7.5, 0.75], 1, 'high', '+', 'pos';...
    [7.5, 0.75], 0, 'high', 'c', 'const';...
    [7.5, 0.75], -1, 'high', '-', 'neg';...
    [2.5, 0.75], 1, 'low', '+', 'pos';...
    [2.5, 0.75], 0, 'low', 'c', 'const';...
    [2.5, 0.75], -1, 'low', '-', 'neg';...
    };

for i=1:size(do_these, 1)
    
    trust_lvl_S0 = do_these{i, 1}(1);
    trust_lvl_S0_sd = do_these{i, 1}(2);
    
    slope_H = do_these{i, 2};
    highlow_text = do_these{i, 3};
    slope_text = do_these{i, 4};
    slope_text2 = do_these{i, 5};
    
    
    
    % let's make an S0 of:
    % S0=[0 0 0 0 0 0 0 0.5 1 0.5]';
    % S0=[0.5 1 0.5 0 0 0 0 0 0 0 ]';
    % S0 = S0/sqrt(S0'*S0);
    
    rating_states = [0.5:9.5];
    % DARK INCANTATION:
    f = @(x, m, s) 1/sqrt(2*pi*s^2) * exp (-1*(x-m)^2 / (2*s^2));
    S0_vec = @(m, s) arrayfun(@(x) f(x, m, s), [0.5:9.5])';
    S0_vec_n = @(m, s) S0_vec(m, s) / sqrt(S0_vec(m, s)' * S0_vec(m, s));
    
    
    S0 = S0_vec_n(trust_lvl_S0, trust_lvl_S0_sd);
    
    
    d = [-4.5:4.5];
    d = slope_H*d;
    if slope_H == 0
        d = [1 1 1 1 1 1 1 1 1 1];
    end
    
    a=1;
    b=1;
    
    H = [
        d(1)     b     0     0     0     0     0     0   0   0
        a    d(2)     b     0     0     0     0     0   0   0
        0     a    d(3)     b     0     0     0     0   0   0
        0     0     a     d(4)     b     0     0    0   0   0
        0     0     0     a     d(5)     b     0    0   0   0
        0     0     0     0     a     d(6)     b    0   0   0
        0     0     0     0     0     a     d(7)    b   0   0
        0     0     0     0     0     0     a     d(8)    b   0
        0     0     0     0     0     0     0      a     d(9) b
        0    0      0     0     0     0     0     0     a     d(10)
        ];
    
    time_scaling = 0.1;
    dt = 1 * time_scaling;
    
    U = expm(-1i*dt*H); % this is a constant U
    
    St = S0;
    state_evolution = S0;
    reliability_as_per_qmodel = [];
    NT=200;
    for i=(1:NT)
        
        U = expm(-1i*dt*H); % this is a constant U
        St = U*St;
        
        reliability_as_per_qmodel(end+1) = rating_states * (abs(St).^2);
        
        state_evolution = [state_evolution St];
        
        
    end
    
    
    
    ul = 0.05; % upp limit
    dl = 0.20; % down limit
    ll = 0.15; % left limit
    rl = 0.03; % right limit
    sh = 0.01; % spacing horz
    sv = 0.01; % spacing vert
    
    NY = 10; % how manycols
    NX = 1;  % how many rows
    
    
    p = [];
    
    for i=1:(NX*NY)
        
        [ix,iy] = ind2sub([NX, NY],i);
        
        xpos = ll + (ix-1) * (1-ll-rl-(NX-1)*sh)/NX + (ix-1) * sh;
        xsize = (1-ll-rl-(NX-1)*sh)/NX;
        
        ypos = 1 - ul - (NY-iy+1) * (1-ul-dl-(NY-1)*sv)/NY - (NY-iy) * sv;
        ysize = (1-ul-dl-(NY-1)*sv)/NY;
        
        p(i,:) = [xpos ypos xsize ysize];
    end
    
    
    
    rating_labels = {'[0-10]', '[10-20]', '[20-30]', '[30-40]', '[40-50]', '[50-60]', ...
        '[60-70]', '[70-80]', '[80-90]', '[90-100]'};
    
    % m=state_evolution.*conj(state_evolution);
    m = abs(state_evolution);
    fh=figure('color','w');
    ahs = [];
    for i=1:10
        ahs(i) = axes('parent',fh, 'position', p(i,:));
        set(ahs(i),'nextplot','add');
        plot(((1:(NT+1))-1)*time_scaling, m(i,:),'k');
        ylim([0 1]);
        xlim([0 (NT+1)]*time_scaling);
        if i~=1
            set(ahs(i),'xticklabel',[]);
            set(ahs(i),'ytick',[0 1]);
            set(ahs(i),'yticklabel',{'', '1'});
        else
            set(ahs(i),'ytick',[0 1]);
            set(ahs(i),'yticklabel',[0  1]);
            
            % xlabel('time');
            set(ahs(i),'xticklabel',[]);
        end
        set(ahs(i),'ylim', [0 1]);
        % set(ahs(i),'yticklabel',[0 1],'ytick',[0 1]);
        % ylabel(rating_labels{i},'rotation', 60);
        
    end
    new_ah = axes('parent', fh, 'position', [ll - 0.13, dl, 0.09, 1-ul-dl]);
    bar(S0,'parent', new_ah,'k');
    set(new_ah,'xlim',[0.5 10.5]);
    set(new_ah,'xtick', [1:10], 'xticklabel', rating_labels,'xticklabelrotation', 90);
    th=text(new_ah, 5.5, 1, 'Initial State Vector (S0)');
    set(th,'horizontalalignment','center','verticalalignment','top');
    set(th,'rotation',90);
    box off
    camroll(90);
    
    
    new_ah2 = axes('parent', fh, 'position', [ll, dl-0.13, 1-rl-ll, 0.1]);
    ph=plot(((1:NT)-1)*time_scaling, reliability_as_per_qmodel, 'k');
    set(ph,'parent', new_ah2);
    set(new_ah2,'xlim',[0 (NT)]*time_scaling);
    set(new_ah2,'ylim',[0 10]);
    set(new_ah2,'ytick',[0 5 10]);
    set(new_ah2,'yticklabel',[0 50 100]);
    % set(new_ah,'xtick', [1:10], 'xticklabel', rating_labels,'xticklabelrotation', 90);
    th=text(new_ah, 5.5, 1, 'Initial State Vector (S0)');
    set(th,'horizontalalignment','center','verticalalignment','top');
    set(th,'rotation',90);
    ylabel('rating');
    xlabel('time');
    
    
    
    anno_ax = axes('parent', fh,'position', [0 0 1 1],'visible','off');
    th=text(anno_ax, 0.5, 0.99, sprintf(' H%s Quantum state vector (St) evolution (amplitude)', slope_text));
    set(th,'fontsize', 9, 'horizontalalignment','center','verticalalignment','top');
    
    set(fh,'paperunits','centimeters');
    set(fh,'papersize',[20 15]);
    set(fh,'paperposition',[0 0 20 15]);
    
    to_save_file = sprintf('figures_examples/new_example_%sS0_H%s', highlow_text, slope_text);
    print('-djpeg','-r300', to_save_file);
    print('-dsvg', to_save_file);
    
    
    
end

