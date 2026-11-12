library(OptSampling)
check_support <- function(pi, sigma) {
  q <- systematic_support(pi, sigma)
  stopifnot(abs(sum(q$probabilities) - 1) < 1e-10,
            all(lengths(q$samples) == round(sum(pi))))
  inc <- numeric(length(pi))
  for (j in seq_along(q$samples)) inc[q$samples[[j]]] <- inc[q$samples[[j]]] + q$probabilities[j]
  stopifnot(max(abs(inc - pi)) < 1e-9)
  q
}
set.seed(72)
for (pi in list(rep(.25, 12), c(.2,.35,.6,.85,1), rep(1,4))) {
  X <- matrix(runif(length(pi)*2, 1, 4), ncol=2)
  fit <- vns_systematic(pi, X, verbose=FALSE)
  q <- check_support(pi, fit$sigma)
  exact <- sum(vapply(q$samples, function(s)
    sum(((colSums(X[s,,drop=FALSE]/pi[s])-colSums(X))/colSums(X))^2), numeric(1))*q$probabilities)
  stopifnot(abs(exact-fit$final_objective)<1e-8,
            fit$final_objective <= fit$initial_objective+1e-9)
  stopifnot(length(draw_systematic_sample(fit, seed=2))==round(sum(pi)))
}
pi <- rep(.25,12)
X <- matrix(runif(24,1,4),ncol=2)
starts <- replicate(3,sample.int(12),simplify=FALSE)
a <- vns_multistart(pi,X,starts=starts,time_budgets=c(1,2),verbose=FALSE)
b <- vns_multistart(pi,X,starts=starts,time_budgets=c(1,2),verbose=FALSE)
stopifnot(identical(a$sigma,b$sigma),
 a$final_objective == min(vapply(a$multistart$fits, function(z) z$final_objective,numeric(1))))
manual <- lapply(starts, function(s) vns_systematic(pi,X,start=s,verbose=FALSE))
stopifnot(abs(a$final_objective-min(vapply(manual,function(z) z$final_objective,numeric(1))))<1e-9)
# Force continuation by exhausting accepted-move budgets.
c <- vns_multistart(pi,X,starts=starts[1],max_iter=0,time_budgets=c(1,2),verbose=FALSE)
stopifnot(nrow(c$multistart$attempts)==2)
for (bad in list(0,NA_real_,-Inf,numeric()))
 stopifnot(inherits(try(vns_multistart(pi,X,time_budgets=bad,verbose=FALSE),silent=TRUE),'try-error'))
# Check local certificate against every transposition with unequal probabilities.
pi <- c(.2,.35,.6,.85,1)
X <- matrix(seq_len(10),ncol=2)
f <- vns_systematic(pi,X,verbose=FALSE)
stopifnot(f$diagnostics$type2_certified)
for(i in 1:4) for(j in (i+1):5) {
 s <- f$sigma; s[c(i,j)] <- s[c(j,i)]
 stopifnot(vns_objective(s,pi,X) >= f$final_objective-1e-9)
}
# Both spatial criteria and rank-deficient balance.
coords <- cbind(seq_len(12),sin(seq_len(12)))
for(metric in c('sb','ib')) {
 control <- list(metric=metric)
 if(metric=='ib') control$W <- 1/(1+as.matrix(dist(coords)))
 f <- vns_systematic(rep(.25,12),coords,objective='spread',objective_control=control,verbose=FALSE)
 stopifnot(is.finite(f$final_objective), f$final_objective<=f$initial_objective+1e-9)
}
f <- vns_systematic(rep(.25,12),cbind(coords,coords[,1]),
 objective_control=list(metric='mahalanobis'),verbose=FALSE)
stopifnot(is.finite(f$final_objective))
pi <- c(.18,.31,.43,.52,.67,.89)
y <- c(2,5,3,8,7,11)
f <- cos_design(y,pi,pricing='exact',seed=12,verbose=FALSE)
stopifnot(f$certified, f$max_inclusion_error<1e-7,
 abs(sum(f$probabilities)-1)<1e-8,
 abs(evaluate_cos_design(f,y)$mse-f$objective)<1e-8,
 length(draw_cos_sample(f))==3)
# Independent global oracle: enumerate all feasible basic designs (20 choose 6).
samples <- t(combn(6,3))
A <- sapply(seq_len(nrow(samples)),function(j) as.numeric(1:6 %in% samples[j,]))
cost <- apply(samples,1,function(s) (sum(y[s]/pi[s])-sum(y))^2)
best <- Inf
for(cols in combn(ncol(A),6,simplify=FALSE)) {
 weights <- tryCatch(solve(A[,cols],pi),error=function(e) NULL)
 if(!is.null(weights) && min(weights)>-1e-9) best <- min(best,sum(weights*cost[cols]))
}
stopifnot(abs(best-f$objective)<1e-6)
