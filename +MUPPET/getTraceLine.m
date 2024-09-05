function [t_trace, f_trace] = getTraceLine(t_stft, f_stft, logpsdm, varargin)
% Find an ideal trace line through a call in a spectrogram by using a
% "shortest path" searching algorithm (Dijkstra's algorithm).
%
% NOTES
% This function uses the following prefix convention for variable names:
%   t = time value(s)
%   f = frequency value(s)
%   i = index of time value(s)
%   j = index of frequency value(s)
%   ij = linear index through time-frequency matrix
% 
%
% Written by Wilfried Beslin
% Last updated 2024-09-05
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEV NOTES:
% - consider adding trace smoothing options as optional input arguments.
% In this case though, might also need to return the interpolated 
% magnitudes. It will also no longer be possible to match the trace
% frequencies with the discrete frequenciy vector of the PSD matrix.
% - 2024-08-23: will try to use multiple Dijkstra paths to determine the
% final trace line insted of just one - will see if this improves results.
% - 2024-09-03: I had attempted to allow negative weights to represent
% nodes that are below the average noise level. The Dijkstra algorithm does
% not support negative weights, but according to the MATLAB documentation,
% the "shortestpathtree" function falls back to the Bellman-Ford algorithm
% instead when negative weights are present. However, this algorithm is not
% behaving as I expect it to and results in erroneous trace lines.
% Therefore I've decided to force all weights to be positive in this
% function and continue using Dijkstra only.

    p = inputParser;
    p.addParameter('PenaltyCoefficient', 0.01, @(val)validateattributes(val,{'numeric'},{'scalar','nonnegative'}));
    p.addParameter('PenaltyExponent', 3, @(val)validateattributes(val,{'numeric'},{'scalar','nonnegative'}));
    p.addParameter('ClippingThreshold', Inf, @(val)validateattributes(val,{'numeric'},{'scalar'}));
    p.addParameter('MaxTimeGap', Inf, @(val)validateattributes(val,{'numeric'},{'scalar','nonnegative'}));
    p.addParameter('MaxFreqGap', Inf, @(val)validateattributes(val,{'numeric'},{'scalar','nonnegative'}));
    p.addParameter('NumAveragingPaths', 1, @(val)validateattributes(val,{'numeric'},{'scalar','positive'}));

    p.parse(varargin{:})
    penalty_coeff = p.Results.PenaltyCoefficient;
    penalty_exp = p.Results.PenaltyExponent;
    clipping_th = p.Results.ClippingThreshold;
    max_t_gap = p.Results.MaxTimeGap;
    max_f_gap = p.Results.MaxFreqGap;
    num_ave = p.Results.NumAveragingPaths;
    
    % define the penalty function for jumps across fequency
    penalty_fcn = @(df) penalty_coeff.*(df.^penalty_exp) + 1;
    
    % get counts from spectrogram
    [nf, nt] = size(logpsdm);

    % find the peak of the spectrogram
    [logpsd_peak, ij_peak] = max(logpsdm(:));
    [j_peak, i_peak] = ind2sub([nf,nt], ij_peak);
    
    % identify all nodes that are above threshold
    is_good_node = logpsdm >= clipping_th;
    
    % create a matrix of node weights based on PSD intensity. 
    % For this, the spectrogram must be scaled such that the peak is 0 dB,
    % to ensure that all weights can be positive (required for Dijkstra)
    logpsdm_anal = logpsdm - logpsd_peak;
    all_node_weights = -logpsdm_anal;
    
    % create directional networks representing all possible paths from the
    % peak to the end (forward) and beginning (backward) of the spectrogram
    
    %%% forward
    fw_node_weights = all_node_weights(:,(i_peak+1):end);
    fw_net = makePathNetwork(j_peak, fw_node_weights, f_stft, penalty_fcn);
    
    %%% backward
    bw_node_weights = all_node_weights(:,(i_peak-1):-1:1);
    bw_net = makePathNetwork(j_peak, bw_node_weights, f_stft, penalty_fcn);
    
    % find the best paths through the forward and backward networks
    j_trace_fw = getBestPath(fw_net, fw_node_weights, num_ave);
    j_trace_bw = getBestPath(bw_net, bw_node_weights, num_ave);
    
    % combine the forward and backward paths into a single trace line
    j_trace_full = [fliplr(j_trace_bw), j_peak, j_trace_fw];
    i_trace_full = 1:nt;
    
    t_trace_full = t_stft(i_trace_full);
    f_trace_full = f_stft(j_trace_full);
    
    % create a logical matrix identifying the full trace nodes
    is_full_trace_node = false(nf,nt);
    is_full_trace_node(sub2ind([nf,nt], j_trace_full, i_trace_full)) = true;
    
    % identify significant trace nodes
    is_good_trace_node = is_full_trace_node & is_good_node;
    l_t_good_trace = any(is_good_trace_node, 1);
    
    % get the start and end points of the final trace based on which
    % points cross the threshold and tolerance for gaps
    [i_trace_start_t, i_trace_stop_t] = getSequencEndPoints(t_trace_full, l_t_good_trace, i_peak, max_t_gap);
    [i_trace_start_f, i_trace_stop_f] = getSequencEndPoints(f_trace_full', l_t_good_trace, i_peak, max_f_gap);
    i_trace_start = max([i_trace_start_t, i_trace_start_f]);
    i_trace_stop = min([i_trace_stop_t, i_trace_stop_f]);
    
    t_trace = t_trace_full(i_trace_start:i_trace_stop);
    f_trace = f_trace_full(i_trace_start:i_trace_stop);
end


%% makePathNetwork --------------------------------------------------------
function path_net = makePathNetwork(r_start, jumpnode_weights, penalty_var, penalty_fcn)
% Create a weighted directional graph representing every possible path from
% a starting point through to the end of a rectangular network of nodes.
% Intended to represent (partial) trace lines through a spectrogram.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Here, "jump node" means any node except the starting node.
    % Jump nodes are arranged in a rectangular matrix, and paths flow
    % through every node from left to right.

    % get counts from the jump node matrix
    [n_rows, n_jumpnode_cols] = size(jumpnode_weights);
    n_jumpnodes = prod([n_rows,n_jumpnode_cols]);
    n_jumpnodes_connected = n_jumpnodes - n_rows;
    
    % create matrix of indices for all jump nodes
    jumpnode_idxmat = reshape(1:n_jumpnodes, n_rows, n_jumpnode_cols) + 1;
    
    % create matrices representing jumps from source nodes to target nodes.
    % Rows = target nodes, Columns = source nodes.
    
    %%% source
    jumpnode_src_mat = repmat((1:n_jumpnodes_connected)+1, n_rows, 1);
    full_src_mat = [ones(n_rows,1), jumpnode_src_mat];
    
    %%% target
    jumpnode_tgt_mat = repmat((0:(n_rows-1))', 1, n_jumpnodes_connected) + n_rows.*repelem(1:(n_jumpnode_cols-1), 1, n_rows) + 2;
    full_tgt_mat = [jumpnode_idxmat(:,1), jumpnode_tgt_mat];
    
    % create matrix of original weights for each jump
    % (this consists simply of the weight of the corresponding target node)
    initial_edge_weight_mat = repelem(jumpnode_weights(:,2:end), 1, n_rows);
    initial_edge_weight_mat = [jumpnode_weights(:,1), initial_edge_weight_mat];
    
    % create matrix of penalty factors for each jump
    penalty_factor_refmat = abs(repmat(penalty_var, 1, n_rows) - repmat(penalty_var', n_rows, 1));
    edge_penalty_factor_mat = repmat(penalty_factor_refmat, 1, n_jumpnode_cols-1);
    edge_penalty_factor_mat = [abs(penalty_var - penalty_var(r_start)), edge_penalty_factor_mat];
    edge_penalty_mat = penalty_fcn(edge_penalty_factor_mat);
    
    % combine the original weight matrix with the penalty matrix
    final_edge_weight_mat = initial_edge_weight_mat .* edge_penalty_mat;
    
    % create the source node, target node, and edge weight vectors that
    % will be used to construct the directional graph object
    src_vec = reshape(full_src_mat, 1, numel(full_src_mat));
    tgt_vec = reshape(full_tgt_mat, 1, numel(full_tgt_mat));
    weight_vec = reshape(final_edge_weight_mat, 1, numel(final_edge_weight_mat));
    
    % create digraph object
    path_net = digraph(src_vec, tgt_vec, weight_vec);
end


%% getBestPath ------------------------------------------------------------
function [j_best_path, j_all_paths, is_candidate_path] = getBestPath(path_net, node_weights, num_ave)
% Returns an estimated most likely path through a rectangular flow network
% representing a time-frequency spectrogram by averaging the top few 
% "shortest" paths from a source node to the target nodes at the end of the
% network (or just using the top-ranking path if "averaging" one line).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    n_nodes = height(path_net.Nodes);
    [nf_net, nt_net] = size(node_weights);

    % get target nodes
    tgt = ((n_nodes-nf_net)+1):n_nodes;
    
    % get shortest paths to all target nodes
    [all_paths, path_dists] = shortestpathtree(path_net, 1, tgt, 'OutputForm','cell');
    all_path_mat = cell2mat(all_paths);
    
    % get the row position of nodes in every path
    [j_all_paths,~] = ind2sub([nf_net,nt_net], all_path_mat(:,2:end)-1);
    
    % get rank of paths based on total weight, where smaller is better
    [~, path_rank] = sort(path_dists, 'ascend');
    
    % determine the best path from the top candidate paths
    is_candidate_path = false(size(path_rank));
    is_candidate_path(path_rank(1:num_ave)) = true;
    if num_ave == 1
        % just use the top path
        j_best_path = j_all_paths(is_candidate_path,:);
    else
        % step through each column and create the final path by getting a
        % weighted average row position of the top nodes
        j_best_path = zeros(1,nt_net);
        for ii = 1:nt_net
            j_ii = j_all_paths(is_candidate_path, ii);
            weights_ii = node_weights(j_ii, ii);

            j_ave_ii = sum(j_ii.*weights_ii)/sum(weights_ii);
            j_best_path(ii) = round(j_ave_ii);
        end

        %** DEBUG
        %%% Here, I'm planning to plot:
        %%% - all paths
        %%% - the top paths (5th percentile weight)
        %%% - the "best" averaged path
        %{
        fig = figure(2);
        clf(fig);
        ax = axes();
        ax.NextPlot = 'add';

        surf(-node_weights, 'EdgeColor','none');
        axis(ax, 'xy');
        view(ax, 0, 90);

        % plot lines
        idx_top = find(is_candidate_path);
        idx_not_top = find(~is_candidate_path);

        %%% plot non-top lines
        for ii = 1:numel(idx_not_top)
            idx_line = idx_not_top(ii);
            plot3(1:nt_net, j_all_paths(idx_line,:), ones(1,nt_net), 'wo-')
        end

        %%% plot top lines
        for ii = 1:numel(idx_top)
            idx_line = idx_top(ii);
            plot3(1:nt_net, j_all_paths(idx_line,:), ones(1,nt_net), 'co-', 'LineWidth',1.5)
        end

        % plot "best" line
        plot3(1:nt_net, j_best_path, ones(1,nt_net), 'rx-', 'LineWidth',2)

        %title(plot_name)

        keyboard
        %}
        %** END DEBUG
    end
end


%% getSequenceEndpoints ---------------------------------------------------
function [i_start, i_stop] = getSequencEndPoints(full_seq, is_good, i_origin, max_gap)
% This function determines where gaps occur in a sequence, and returns a
% start and end point about an origin point point that depends on the siZe
% of the gaps (large gaps cut the sequence short).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % backward analysis
    i_bw = i_origin:-1:1;
    j_stop_bw = findEndPoint(full_seq(i_bw), is_good(i_bw), max_gap);
    i_start = i_origin - j_stop_bw + 1;
    
    % forward analysis
    i_fw = i_origin:numel(full_seq);
    j_stop_fw = findEndPoint(full_seq(i_fw), is_good(i_fw), max_gap);
    i_stop = i_origin + j_stop_fw - 1;
    
    
    % NESTED FUNCTION: findEndPoint .......................................
    function j_stop = findEndPoint(seq_sub, is_good_sub, max_gap)
        
        % get differences between points
        seq_sub_diff = abs(diff(seq_sub(is_good_sub)));
        
        % find the first instance where the difference is greater than the
        % max allowed. Here I am assuming that the sequence never starts in
        % a gap.
        j_stop = find(seq_sub_diff > max_gap, 1, 'first');
        if isempty(j_stop)
            j_stop = find(is_good_sub, 1, 'last');
        end
    end
end