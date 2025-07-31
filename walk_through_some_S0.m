% figure out some basics:

fh=figure;



WIDEN_IT = 1;
time_scaling = 0.1;
catch_m_i1 = 2;
catch_m_i2 = 2;
states = [-3:3];
slope = 1;
offset = [-3 -2 -1 0 1 2 3] * 0;

for i=1:7
    subplot(7,1,i);
    hold on;
    
    S0 = zeros(1,7);
    S0(i) = 1;
    
    if WIDEN_IT
        if i>1
            S0(i-1) = 0.5;
        end
        if i<7
            S0(i+1) = 0.5;
        end
    end

    % set S0 constant, but now try different slopes.
    % S0 = [1 1 1 1 1 1 1];

    S0 = S0/sqrt(S0*conj(S0'));
    S0 = reshape(S0, numel(S0), 1);

    dt = 1 * time_scaling;
    d = states * slope + offset(i);
    a=1;
    b=1;
    
    H = [    
       d(1)     b     0     0     0     0     0
        a    d(2)     b     0     0     0     0
        0     a    d(3)     b     0     0     0
        0     0     a     d(4)     b     0     0
        0     0     0     a     d(5)     b     0
        0     0     0     0     a     d(6)     b
        0     0     0     0     0     a     d(7)
     ];
 
    U = expm(-1i*dt*H);
    U1 = U;
    mv = [-3:3];
    
    
    meas=[]; m=S0; for j=1:1000; m(:, end+1) = U*m(:, end); meas(end+1)=mv * (m(:, end).*conj(m(:,end))); end;

    disp([i catch_m_i1]);
    if i==catch_m_i1
            disp('hallo');
            m1=m;
    end
    
    
    plot(meas);

    a=1;
    b=1;
    
    d = -1*d;


    H = [    
       d(1)     b     0     0     0     0     0
        a    d(2)     b     0     0     0     0
        0     a    d(3)     b     0     0     0
        0     0     a     d(4)     b     0     0
        0     0     0     a     d(5)     b     0
        0     0     0     0     a     d(6)     b
        0     0     0     0     0     a     d(7)
     ];
 
    U = expm(-1i*dt*H);
    U2 = U;
    meas=[]; m=S0; for j=1:1000; m(:, end+1) = U*m(:, end); meas(end+1)=mv * (m(:, end).*conj(m(:,end))); end;
    
    if i==catch_m_i2
        m2=m;
    end
    plot(meas, 'r');
    ylim([-3 3]);
    
    title(sprintf('S0 = %2.2f   ', S0'));
    if i==catch_m_i1
        use_this_title = sprintf('S0 = %2.2f   ', S0');
    end
    
    
    
    
end


nfh = figure; 
tmp_mat1 = m1.*conj(m1);
tmp_mat2 = m2.*conj(m2);

for i=1:7
    subplot(7,1,i); 
    yd = tmp_mat1(i,:);
    plot(yd); xlim([0 size(m1, 2)]); 
    text_a = mean(yd);

    hold on;
    yd = tmp_mat2(i,:);
    plot(yd);xlim([0 size(m2, 2)]);
    text_b = mean(yd);
    
    
    title(sprintf('LVL: %d   H+, mean= %2.2f   H-, mean= %2.2f', states(i), text_a, text_b));
end


anno_ax = axes('parent',nfh,'visible','off');
set(anno_ax, 'position', [0 0 1 1]);
th = text(anno_ax, 0.5, 0.999, use_this_title);
set(th,'horizontalalignment','center','verticalalignment','top');
set(th,'fontweight','bold','fontsize',10);

