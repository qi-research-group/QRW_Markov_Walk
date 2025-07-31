% this scripts loads saved_model_parameters.mat from x1_ script
% and then we will explore stuff and make some graphs too.

load saved_model_parameters_Markov.mat
p=saved_model_parameters_Markov;





cost_f=[];
for i=1:numel(p)
    cost_f = [cost_f p{i}{2}];
end
markov_costf = cost_f;


% slight helper:
get_s_num = @(i) sprintf('%02d', str2num(p{i}{5}(9:11)));

% Create figure
fh=figure('color','w');
ah=axes('parent',fh); 
set(ah,'nextplot','add');


% Plot individual data points
x_vals = rand(size(cost_f)) * 0.3 + 0.85;
plot(x_vals, cost_f, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
for i=1:numel(x_vals)
    th=text(x_vals(i), cost_f(i), [' ' get_s_num(i)]);
    set(th,'horizontalalignment','left','verticalalignment','middle');
end



% Compute statistics
data_min = min(cost_f);
data_max = max(cost_f);
data_median = median(cost_f);
data_q1 = quantile(cost_f, 0.25);
data_q3 = quantile(cost_f, 0.75);

my_left_pos = 0.7;
my_right_pos = 1.3;

% Draw box
rectangle('Position', [my_left_pos, data_q1, my_right_pos-my_left_pos, data_q3 - data_q1], ...
    'EdgeColor', 'k', 'LineWidth', 1.5, 'FaceColor', 'none');

% Draw whiskers
line([1 1], [data_min data_q1], 'Color', 'k', 'LineWidth', 1.5); % Lower whisker
line([1 1], [data_q3 data_max], 'Color', 'k', 'LineWidth', 1.5); % Upper whisker

% Draw median line
line([my_left_pos my_right_pos], [data_median data_median], 'Color', 'k', 'LineWidth', 1.5);



% Compute and plot spread (boxplot equivalent)
% boxplot(cost_f, 'Colors', 'k', 'Whisker', 1.5, 'Symbol', 'ko', ...
%     'Widths', 0.9, 'Labels', {''}, 'OutlierSize', 1);

set(ah,'xcolor','none');
xlim([0.5 1.5]);
ylim([0 1.5]);

box off;

ylabel('MAE modeled - predicted reliability rating');
title('Final cost model fit all 34 participants');




% now - let's check model params:
fh=figure('color','w');

    
ul = 0.05; % upp limit
dl = 0.20; % down limit
ll = 0.15; % left limit
rl = 0.03; % right limit
sh = 0.05; % spacing horz
sv = 0.15; % spacing vert

NY = 2; % how manycols
NX = 2;  % how many rows


p = [];

for i=1:(NX*NY)

    [ix,iy] = ind2sub([NX, NY],i);

    xpos = ll + (ix-1) * (1-ll-rl-(NX-1)*sh)/NX + (ix-1) * sh;
    xsize = (1-ll-rl-(NX-1)*sh)/NX;

    ypos = 1 - ul - (NY-iy+1) * (1-ul-dl-(NY-1)*sv)/NY - (NY-iy) * sv;
    ysize = (1-ul-dl-(NY-1)*sv)/NY;

    p(i,:) = [xpos ypos xsize ysize];
end



m=[];
for i=1:numel(saved_model_parameters_Markov)
    m=[m; saved_model_parameters_Markov{i}{1}];
end

my_ylabels = {'\gamma','\mu_{+}','\mu_{-}','\sigma'};

ahs = [];
for i=1:4
    ahs(i) = axes('parent',fh,'position', p(5-i,:));
end
ahs = ahs([2 1 4 3]);
    
for i=1:4
    % ahs(i) = axes('parent',fh, 'position', p(i,:));
    set(fh,'currentaxes', ahs(i));
    x_labels = arrayfun(get_s_num, 1:34,'uniformoutput',false);
    ph=plot(m(:, i),'k');
    set(ahs,'xtick', 1:34, 'xticklabel',x_labels,'xticklabelrotation', 90);
    xlim([0 35]);
    % ylabel(my_ylabels{i});
    xl = get(gca,'xlim');
    yl = get(gca,'ylim');
    
    th = text(xl(2) - diff(xl)*0.15, yl(2), my_ylabels{i});
    set(th,'fontsize', 11);
    
    
    box off;
end




% Customize appearance

% ylabel('Values');
% title('Data Spread with Outliers and Individual Points');