function [inh_comm_ids, comm_resp_frac, comm_qvals] = identifyInhibitoryCommunities( ...
        community_ids, resp_mask, fdr_level)
% IDENTIFYINHIBITORYCOMMUNITIES  Flag Louvain communities enriched for responsive units.
%
% One-sided hypergeometric enrichment test per community — is its responsive
% count higher than expected by chance under the population-wide responsive
% rate, given the community's size — with Benjamini-Hochberg FDR correction
% across communities. Scales correctly with community size (a few responsive
% units in a 10-unit community is far less surprising than the same count in
% a 200-unit community), unlike a fixed relative/absolute fraction threshold
% which lets small, noisy communities clear an arbitrary bar by chance.
%
% INPUTS:
%   community_ids - (N x 1) or (1 x N) Louvain community assignment per unit
%   resp_mask     - (N x 1) or (1 x N) logical, true = responsive
%   fdr_level     - Benjamini-Hochberg FDR level (e.g. 0.05)
%
% OUTPUTS:
%   inh_comm_ids   - (1 x K) community IDs passing the FDR-corrected test
%   comm_resp_frac - (1 x n_comm) responsive fraction per community
%   comm_qvals     - (1 x n_comm) FDR-corrected p-values (q-values) per community

community_ids = community_ids(:);
resp_mask     = logical(resp_mask(:));
n_comm        = max(community_ids);

comm_size     = zeros(1, n_comm);
comm_resp_cnt = zeros(1, n_comm);
for c = 1:n_comm
    c_mask           = community_ids == c;
    comm_size(c)     = sum(c_mask);
    comm_resp_cnt(c) = sum(resp_mask(c_mask));
end
comm_resp_frac = comm_resp_cnt ./ max(comm_size, 1);

% One-sided hypergeometric enrichment: P(X >= k) for X ~ Hyge(N_all, n_resp, comm_size)
N_all     = numel(community_ids);
n_resp    = sum(resp_mask);
raw_pvals = ones(1, n_comm);
for c = 1:n_comm
    if comm_resp_cnt(c) == 0; continue; end
    raw_pvals(c) = 1 - hygecdf(comm_resp_cnt(c) - 1, N_all, n_resp, comm_size(c));
end

comm_qvals   = StatisticalTest.fdrCorrection(raw_pvals(:))';
inh_comm_ids = find(comm_qvals < fdr_level);
end
