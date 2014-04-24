function pkstats = dpeak_heuristic(x, varargin)
% DPEAK_HEURISTIC Detect peaks for a single analyte in lxb data.
%   PKSTATS = DPEAK_HEURISTIC(X) detect peaks for intensities X. X
%   is a 1d vector of fluorescent intensities for the given analyte.
%   PKSTATS is a structure array with the detected peaks and support
%   information.
%
%   PKSTATS = DPEAK_HEURISTIC(X, 'Param1, 'Value1',...)
%   Specify optional parameter/value pairs:
%
%   'logxform' : boolean, log2 transform the data before detecting peaks.
%       Default is true.
%   'lowthresh' : scalar, log2 lower threshold for expression. Values below
%       lowthresh are ignored. Default is 4.
%   'highthresh' : scalar, log2 upper threshold for expression. Values
%       above highthresh are ignored. Default is 15.
%       false.
%   'minbead : integer, Minimum good beads to make a call. Default is 10.
%   'pkmethod' : string, Peak calling method. Valid options are
%       {'kmeans_viable', 'median'}. Default is kmeans_viable.
%       'kmeans_viable': Applies k-means with k=2,3,4. Picks the K based on
%           the clustering that yields the support gradient closest to
%           expect_support_pct. In addition it applies a heuristic to
%           remove non-viable peaks that lack support.
%       'median': Reports the median value of the distribution. This is the
%           default for 'notduo' analytes and for analytes which lack
%           sufficient support to make a call.
%   'debug' : Print debugging info. Logical [true] 
%   'expect_npeak' : integer, Expected number of peaks. Default is 2.
%   'expect_support_pct': vector, Expected bead support gradient. Default
%       is [65 35]s


toolName = mfilename;
    
pnames = {'debug', ...
    'highthresh', ...
    'ksmethod', ...
    'logxform', ...
    'lowthresh', ...
    'max_k', ...
    'merge_close_peaks', ...
    'minbead', ...
    'min_peak_support', ...
    'min_peak_support_pct', ...
    'ncomp', ...
    'opt_support_pct', ...
    'out', ...
    'overwrite', ...
    'pkmethod', ...
    'rpt', ...
    'sortbyht', ...
    'subtractbg', ...
    'title'};

dflts = { false, ...
    14, ...
    'kmeans_viable',...
    true, ...
    4, ...
    4, ...
    true, ...
    20, ...
    10, ...
    10, ...
    2, ...
    [65 35], ...
    pwd, ...
    false, ...
    'kspeak', ...    
    'hist', ...      
    false, ...
    false, ...
    ''};    

arg = parse_args(pnames, dflts, varargin{:});

nbead = length(x);

%convert x to log scale
if (arg.logxform)
    x = safe_log2(x - arg.subtractbg);
end

%censor beads with high / low expression
badbeads = (x < arg.lowthresh | x > arg.highthresh);
xcensored=x;
xcensored(badbeads)=[];
ngoodbead = length(xcensored);

%default value is median of good beads
medexp = median(xcensored);
if (arg.logxform)
   medexp = round(pow2(medexp));
end
% bad or missing data defaults to one (zero on log scale)
if isempty(medexp) || isnan(medexp) || isequal(nbead,0)
    medexp = 1;
end

%default stats
pkstats =  struct('pkexp', medexp,...
    'pksupport', ngoodbead, 'pksupport_pct',100,'pkheight', 1,...
    'totbead', nbead, 'ngoodbead', ngoodbead,...
    'medexp', medexp, ...
    'method','median');

if (ngoodbead < arg.minbead || isequal(arg.ksmethod, 'median'))
    % dont detect peaks just return default
else
    dbg(arg.debug, arg.pkmethod);
    switch (arg.pkmethod)
            
        case 'kmeans_viable'
            % turn off emptycluster warning
            warning('off', 'stats:kmeans:EmptyCluster')
            % optimal K
            best_k = 1;
            % centroids
            best_c = log2(medexp);
            % cluster identity for each K
            best_idx = ones(ngoodbead, 1);            
            % diff between measured and ideal support pct
            best_sup_dist = sum(abs([100 0] - arg.opt_support_pct));
            % Pick an optimal k
            for k=2:arg.max_k
                [idx, c, sumd, d] = kmeans(xcensored, k, ...
                    'emptyaction', 'drop', 'replicates', 5);                
                nc = nnz(~isnan(c));
                if nc < 2
                    % NaN centroids found, skip
                    dbg(arg.debug, 'k=%d NaN centroids found, skipping', k)
                    continue                    
                else
                    sup = accumarray(idx, ones(size(idx)))';
                    srtsup = sort(sup, 'descend');
                    sup_dist = sum(abs(100 * srtsup(1:2) / ngoodbead - arg.opt_support_pct));
                    isviable = sup >= arg.min_peak_support & 100*sup/ngoodbead >= arg.min_peak_support_pct;
                    nviable = sum(isviable);
                    if nviable > 1 && (sup_dist < best_sup_dist)
                        best_sup_dist = sup_dist;                        
                        best_idx = idx;
                        best_c = c;
                        best_c(~isviable) = nan;
                        best_k = nviable;
                    end
                end
            end
            
            % Keep viable clusters
            if nnz(isnan(best_c))
                keep = find(~isnan(best_c));
                best_c = best_c(keep);                
                tmp_idx = nan(ngoodbead, 1);
                for ii=1:best_k                    
                    tmp_idx(best_idx==keep(ii)) = ii;
                end
                % censor non-viable clusters
                valid_beads = ~isnan(tmp_idx);
                xcensored = xcensored(valid_beads);
                best_idx = tmp_idx(valid_beads);
                ngoodbead = nnz(valid_beads);
            end
                        
            dbg(arg.debug, 'optk=%d, sup_dist=%2.2f\n', ...
                best_k, best_sup_dist);

            pksup = accumarray(best_idx, ones(size(best_idx)))';
            pksup_pct = 100*pksup/ngoodbead';
            % median of clusters instead of the mean
            pkexp = zeros(1, best_k);
            for ii=1:best_k
                pkexp(1, ii) = median(xcensored(best_idx==ii));
            end
                        
            % compute smoothed kernel, lookup pk heights
            [f, xi] = ksdensity(xcensored);
            pkheight = interp1(xi, f, pkexp);
            [srt_pksup_pct, srtidx] = sort(pksup_pct, 'descend');
            srt_pkexp = pkexp(srtidx);
            srt_pkheight = pkheight(srtidx);
            srt_pksup = pksup(srtidx);
            
            pkstats = struct(...
                'pkexp', srt_pkexp,...
                'pksupport', srt_pksup,...
                'pksupport_pct', srt_pksup_pct,...
                'pkheight', srt_pkheight,...
                'totbead', nbead,...
                'ngoodbead', ngoodbead,...
                'medexp', medexp,...
                'method', arg.pkmethod);
        otherwise
            error('Unknown pkmethod:%s',arg.pkmethod)            
            
    end
        
    % convert expression to linear scale
    if (arg.logxform)
        pkstats.pkexp = round(pow2(pkstats.pkexp));        
    end
                
end

end

