# Update the crossing block on phenotype (burnin and PS control).
# The nRecycled oldest parents are replaced by the best S4 clones.

old <- parents
recycled <- selectInd(S4, nInd = nRecycled, use = 'pheno')
parents <- c(old[-(1:nRecycled)], recycled)
nKept <- sum(old@id %in% parents@id)
