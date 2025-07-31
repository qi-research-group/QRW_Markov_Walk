

% define our subjects
subs = arrayfun(@(x) sprintf('%02d',x), [1 3:20 22:36], 'uniformoutput', false);


% get all of our behav data; history and ratings:
d_h = dir('../../history_files/history_*txt');
d_r = dir('../../history_files/ratings_*txt');



h=[];
for i=1:numel(d_h)
    this_h = load([d_h(i).folder filesep d_h(i).name]);
    h=[h;this_h];
end

r=[];
for i=1:numel(d_r)
    this_r = load([d_r(i).folder filesep d_r(i).name]);
    r=[r;this_r];
end



% now we need to define some things. % our delta-t, for example
% for now we just set it as 3.5.
% we CAN/probably should incorporate delta-t from the data itself.

dt = 3.5;
nstates = 10;

% define some random walk(s)
h_match = 
h_mismatch = 



%    m = np.diag(diag).astype(np.double)
%     
%     # set the off diagonals, too:
%     row_idx = np.arange(1, m.shape[0]).tolist() + np.arange(m.shape[0] - 1).tolist()
%     col_idx = np.arange(m.shape[0] - 1).tolist() + np.arange(1, m.shape[0]).tolist()
%     m[row_idx, col_idx] = drift
%     
%     # set the dissonance (z)
%     row_idx = [0, m.shape[0]-1] # zero-row, last column
%     col_idx = [m.shape[0]-1, 0] # last row, zero column
%     m[row_idx, col_idx] = dissonance
%     
%     return m