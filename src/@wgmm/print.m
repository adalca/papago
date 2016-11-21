function print(wg)
    s = '--- WGMM object:\n';

    s = sprintf('%soperating methods:', s);
    s = sprintf('%s%30s:%s\n', s, 'initialization', wg.initmethod);
    s = sprintf('%s%30s:%s\n', s, 'covariance update', wg.covarUpdateMethod);
    s = sprintf('%s%30s:%s\n', s, 'mu update', wg.muUpdateMethod);
    s = sprintf('%s%30s:%s\n', s, 'logp update', wg.logpUpdateMethod);
    s = sprintf('%s%30s:%s\n', s, 'covariance reconstruction', wg.covarReconMethod);
    s = sprintf('%s%30s:%s\n', s, 'covariance merge', wg.covarMergeMethod);
    fprintf(s);
    sys.structmemory(wg)
end