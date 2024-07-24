function [t_trace, f_trace] = getTraceLine(t_stft, f_stft, psdm, varargin)
% Find an ideal trace line through a call in a spectrogram by using a
% "shortest path" searching algorithm (Dijkstra's algorithm)
%
% Written by Wilfried Beslin
% Last updated 2024-07-24
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEV NOTES:
% - consider adding trace smoothing options as optional input arguments.
% In this case though, might also need to return the interpolated 
% magnitudes. It will also no longer be possible to match the trace
% frequencies with the discrete frequenciy vector of the PSD matrix.

    p = inputParser;
    p.addParameter('PenaltyCoefficient', 0.01, @(val)validateattributes(val,{'numeric'},{'scalar','nonnegative'}));
    p.addParameter('PenaltyExponent', 3, @(val)validateattributes(val,{'numeric'},{'scalar','nonnegative'}));
    p.addParameter('ClippingThreshold', Inf, @(val)validateattributes(val,{'numeric'},{'scalar'}));
    p.addParameter('MaxTimeGap', Inf, @(val)validateattributes(val,{'numeric'},{'scalar','nonnegative'}));
    p.addParameter('MaxFreqGap', Inf, @(val)validateattributes(val,{'numeric'},{'scalar','nonnegative'}));

    p.parse(varargin{:})
    penalty_coeff = p.Results.PenaltyCoefficient;
    penalty_exp = p.Results.PenaltyExponent;
    clipping_weight_th = abs(p.Results.ClippingThreshold);
    max_t_gap = p.Results.MaxTimeGap;
    max_f_gap = p.Results.MaxFreqGap;
    
    % define the penalty function for jumps across fequency
    penalty_fcn = @(df) penalty_coeff.*(df.^penalty_exp) + 1;
    
    % get counts from spectrogram
    [nf, nt] = size(psdm);

    % find the peak of the spectrogram
    [~, i_peak] = max(psdm(:));
    [i_f_peak, i_t_peak] = ind2sub(size(psdm), i_peak);
    
    % create a matrix of node weights based on PSD intensity
    all_node_weights = 10*log10(psdm);
    all_node_weights = abs(all_node_weights - all_node_weights(i_peak));
    
    % identify nodes that are above the weight threshold
    is_good_node = all_node_weights <= clipping_weight_th;
    
    % create directional networks representing all possible paths from the
    % peak to the end (forward) and beginning (backward) of the spectrogram
    
    %%% forward
    fw_node_weights = all_node_weights(:,(i_t_peak+1):end);
    fw_net = makePathNetwork(i_f_peak, fw_node_weights, f_stft, penalty_fcn);
    
    %%% backward
    bw_node_weights = all_node_weights(:,(i_t_peak-1):-1:1);
    bw_net = makePathNetwork(i_f_peak, bw_node_weights, f_stft, penalty_fcn);
    
    % find the path in the forward graph with smallest total weight
    fw_trace_nodes = getShortestPath(fw_net, nf);
    
    % do the same for the backwards network
    bw_trace_nodes = getShortestPath(bw_net, nf);
    
    % translate the trace node IDs into time-frequency points
    [i_f_trace_fw, ~] = ind2sub([nf,nt-i_t_peak], fw_trace_nodes(2:end) - 1);
    [i_f_trace_bw, ~] = ind2sub([nf,i_t_peak-1], bw_trace_nodes(2:end) - 1);
    
    i_f_trace_full = [fliplr(i_f_trace_bw), i_f_peak, i_f_trace_fw];
    i_t_trace_full = 1:nt;
    
    t_trace_full = t_stft(i_t_trace_full);
    f_trace_full = f_stft(i_f_trace_full);
    
    % create a logical matrix identifying the full trace nodes
    is_full_trace_node = false(nf,nt);
    is_full_trace_node(sub2ind([nf,nt], i_f_trace_full, i_t_trace_full)) = true;
    
    % identify significant trace nodes
    is_good_trace_node = is_full_trace_node & is_good_node;
    l_t_good_trace = any(is_good_trace_node, 1);
    
    % get the start and end points of the final trace based on which
    % points cross the threshold and tolerance for gaps
    [i_t_trace_start_t, i_t_trace_stop_t] = getSequencEndPoints(t_trace_full, l_t_good_trace, i_t_peak, max_t_gap);
    [i_t_trace_start_f, i_t_trace_stop_f] = getSequencEndPoints(f_trace_full', l_t_good_trace, i_t_peak, max_f_gap);
    i_t_trace_start = max([i_t_trace_start_t, i_t_trace_start_f]);
    i_t_trace_stop = min([i_t_trace_stop_t, i_t_trace_stop_f]);
    
    t_trace = t_trace_full(i_t_trace_start:i_t_trace_stop);
    f_trace = f_trace_full(i_t_trace_start:i_t_trace_stop);
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


%% getShortestPath --------------------------------------------------------
function shortest_path = getShortestPath(path_net, n_cols)
% Return the "shortest" path from a source node to the edge of a
% rectangular flow network.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    n_nodes = height(path_net.Nodes);

    % get target nodes
    tgt = ((n_nodes-n_cols)+1):n_nodes;
    
    % get shortest paths to all target nodes
    [paths, dists] = shortestpathtree(path_net, 1, tgt, 'OutputForm','cell');
    
    % identify the shortest of all paths
    [~, i_shortest] = min(dists);
    shortest_path = paths{i_shortest};
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