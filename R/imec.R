## imec.R -- exact interventional Markov equivalence class machinery.
##
## Replaces the covered-arc random walk of rcar() by exact operations on the
## I-equivalence class of the current DAG:
##
##   .iessgraph()   the I-essential graph (Hauser & Buehlmann, 2012)
##   .uccgs()       its undirected (chain) components
##   .amos()        the acyclic moral orientations of one component
##   imec.size()    |[D]_I|, exactly
##   imec.sample()  one draw, uniform on [D]_I
##   imec.list()    every member of [D]_I
##
## The last three rest on Proposition 4 of Wienoebst, Bannach & Liskiewicz
## (JMLR 2023), after Hauser & Buehlmann (2012): the undirected components of an
## I-essential graph are chordal, and a DAG lies in the I-MEC iff it is obtained
## by acyclic moral orientations of those components independently of each
## other.  So the class factorises over components, |[D]_I| is the product of
## their AMO counts, and independent uniform draws per component are uniform on
## the class.  analysis/35_prop4_check.R verifies this against brute force.
##
## Graphs are p x p logical matrices with A[i, j] = TRUE meaning i -> j; an
## undirected edge is A[i, j] && A[j, i].

## an arrow u -> v is protected by the targets when some target tells the two
## apart, i.e. exactly one of u and v belongs to some I in the family of targets.
.protects <- function(targets, u, v)
    any(vapply(targets, function(I) ((u %in% I) + (v %in% I)) == 1L, TRUE))

## I-essential graph of a DAG: Chickering's protected-arrow labelling, with
## target protection added.  Mirrors EssentialGraph::replaceUnprotected() of
## pcalg (src/greedy.cpp:981); validated against pcalg::dag2essgraph in
## tests/testthat/test-imec.R.
.iessgraph <- function(A, targets = list(integer(0))) {
    p <- nrow(A)
    G <- A                                   # working graph, mutated in place
    PROTECTED <- 1L; UNDECIDABLE <- 2L; NOT_PROTECTED <- 3L

    arrows <- which(G & !t(G), arr.ind = TRUE)
    if (!nrow(arrows)) return(G)
    flag <- integer(nrow(arrows))
    for (e in seq_len(nrow(arrows)))
        flag[e] <- if (.protects(targets, arrows[e, 1], arrows[e, 2])) PROTECTED
                   else UNDECIDABLE

    adj <- function(u, v) G[u, v] || G[v, u]
    ## an arrow in a v-structure is protected
    for (v in seq_len(p)) {
        inn <- which(arrows[, 2] == v)
        if (length(inn) > 1)
            for (a in seq_along(inn)) for (b in seq_len(a - 1)) {
                s1 <- arrows[inn[a], 1]; s2 <- arrows[inn[b], 1]
                if (!adj(s1, s2)) { flag[inn[a]] <- PROTECTED; flag[inn[b]] <- PROTECTED }
            }
    }

    ## `active` mirrors pcalg's arrowFlags map: an arrow leaves it when it is
    ## turned into an undirected edge, and the configuration scans below must
    ## not see those any more.
    active <- rep(TRUE, nrow(arrows))
    repeat {
        und <- which(flag == UNDECIDABLE & active)
        if (!length(und)) break
        adj <- function(u, v) G[u, v] || G[v, u]
        isPar <- function(u, v) G[u, v] && !G[v, u]
        isNbr <- function(u, v) G[u, v] && G[v, u]
        for (e in und) {
            a <- arrows[e, 1]; b <- arrows[e, 2]; f <- NOT_PROTECTED
            ## (a)  c -> a  with c not adjacent to b
            for (g in which(arrows[, 2] == a & active)) {
                if (f == PROTECTED) break
                cc <- arrows[g, 1]
                if (!adj(cc, b)) f <- if (flag[g] == PROTECTED) PROTECTED else UNDECIDABLE
            }
            ## (c)  a -> c -> b
            for (g in which(arrows[, 2] == b & active)) {
                if (f == PROTECTED) break
                cc <- arrows[g, 1]
                if (isPar(a, cc)) {
                    h <- which(arrows[, 1] == a & arrows[, 2] == cc & active)
                    f <- if (flag[g] == PROTECTED && length(h) && flag[h] == PROTECTED)
                             PROTECTED else UNDECIDABLE
                }
            }
            ## (d)  c1 -> b, c2 -> b, both neighbours of a, c1 and c2 non-adjacent
            inn <- which(arrows[, 2] == b & active)
            for (i1 in seq_along(inn)) {
                if (f == PROTECTED) break
                for (i2 in seq_len(i1 - 1)) {
                    if (f == PROTECTED) break
                    c1 <- arrows[inn[i1], 1]; c2 <- arrows[inn[i2], 1]
                    if (isNbr(a, c1) && isNbr(a, c2) && !adj(c1, c2))
                        f <- if (flag[inn[i1]] == PROTECTED && flag[inn[i2]] == PROTECTED)
                                 PROTECTED else UNDECIDABLE
                }
            }
            flag[e] <- f      # in place, so later arrows in this pass see it
        }
        ## replace the unprotected arrows by lines; count what got settled
        settled <- 0L
        drop <- which(flag == NOT_PROTECTED & active)
        for (e in drop) {
            G[arrows[e, 2], arrows[e, 1]] <- TRUE          # becomes undirected
            active[e] <- FALSE
        }
        settled <- length(drop) + sum(flag[und] == PROTECTED)
        if (settled == 0L) stop("invalid graph passed to .iessgraph()")
    }
    G
}

## undirected chain components of an I-essential graph, as vertex index vectors;
## singletons are dropped, since they contribute one orientation each
.uccgs <- function(E) {
    U <- E & t(E)
    p <- nrow(U); seen <- rep(FALSE, p); out <- list()
    for (v in seq_len(p)) {
        if (seen[v] || !any(U[v, ])) next
        comp <- v; seen[v] <- TRUE; frontier <- v
        while (length(frontier)) {
            w <- frontier[1]; frontier <- frontier[-1]
            nb <- which(U[w, ] & !seen)
            seen[nb] <- TRUE; comp <- c(comp, nb); frontier <- c(frontier, nb)
        }
        if (length(comp) > 1) out[[length(out) + 1L]] <- sort(comp)
    }
    out
}

## every acyclic moral orientation of the component `vs` of the undirected graph
## U, as a list of p x p logical matrices.  Enumerated over vertex orderings and
## de-duplicated: every AMO of a chordal graph is induced by at least one
## ordering, and orienting along an ordering can only fail morality, never
## acyclicity, so no acyclicity test is needed.
.amos <- function(U, vs, max.size = 8L) {
    m <- length(vs)
    if (m > max.size) return(NULL)                       # caller falls back
    p <- nrow(U)
    sub <- U[vs, vs, drop = FALSE]
    ei <- which(sub & upper.tri(sub), arr.ind = TRUE)
    seen <- new.env(parent = emptyenv()); out <- list()
    perms <- .perms(m)
    for (r in seq_len(nrow(perms))) {
        pos <- integer(m); pos[perms[r, ]] <- seq_len(m)
        ## orient every edge from the earlier to the later vertex
        A <- matrix(FALSE, m, m)
        A[cbind(ifelse(pos[ei[,1]] < pos[ei[,2]], ei[,1], ei[,2]),
                ifelse(pos[ei[,1]] < pos[ei[,2]], ei[,2], ei[,1]))] <- TRUE
        if (.has.immorality(A)) next
        k <- paste(which(A), collapse = ",")
        if (!is.null(seen[[k]])) next
        assign(k, TRUE, envir = seen)
        full <- matrix(FALSE, p, p); full[vs, vs] <- A
        out[[length(out) + 1L]] <- full
    }
    out
}

.perms <- function(m) {
    if (m == 1L) return(matrix(1L, 1, 1))
    sub <- .perms(m - 1L)
    do.call(rbind, lapply(seq_len(m), function(i)
        cbind(i, matrix(c(seq_len(m)[-i])[sub], nrow = nrow(sub)))))
}

.has.immorality <- function(A) {
    m <- nrow(A)
    for (v in seq_len(m)) {
        pa <- which(A[, v])
        if (length(pa) > 1)
            for (i in seq_along(pa)) for (j in seq_len(i - 1))
                if (!A[pa[i], pa[j]] && !A[pa[j], pa[i]]) return(TRUE)
    }
    FALSE
}

## ---------------------------------------------------------------------------
## public interface

## One draw, uniform on [D]_I, by Clique-Picking (src/cliquepick.c): the class
## factorises over the chain components of the I-essential graph, and the C
## routine returns a topological order whose induced orientation of each
## component is uniform among its acyclic moral orientations.  Both engines
## call this same code, which is what keeps their random streams aligned.
##
## The one limit is an implementation limit, and it is per chain component:
## vertex sets within a component are 64-bit masks, so a component of more
## than 64 vertices makes the draw decline and the caller fall back. Nothing
## bounds the size of the class, and nothing bounds p -- C_cp_amo_sample()
## returns a full order of all p vertices, the ones outside any component
## appended in index order -- so a DAG of any size whose components are each
## within 64 vertices is sampled exactly.
imec.sample <- function(A, targets = list(integer(0))) {
    p <- nrow(A)
    E <- .iessgraph(A, targets)
    U <- E & t(E)
    ord <- .Call(C_cp_amo_sample, U)
    if (is.null(ord)) return(NULL)
    pos <- integer(p); pos[ord] <- seq_len(p)
    out <- E & !t(E)                               # the directed part is fixed
    ij <- which(U & upper.tri(U), arr.ind = TRUE)
    if (nrow(ij))
        out[cbind(ifelse(pos[ij[, 1]] < pos[ij[, 2]], ij[, 1], ij[, 2]),
                  ifelse(pos[ij[, 1]] < pos[ij[, 2]], ij[, 2], ij[, 1]))] <- TRUE
    out
}

## |[D]_I| by Clique-Picking, in polynomial time. NA when a chain component
## exceeds the 64-vertex implementation limit described above; nothing else
## is bounded, the class size and p included.
imec.size <- function(A, targets = list(integer(0))) {
    E <- .iessgraph(A, targets)
    .Call(C_cp_amo_count, E & t(E))
}

## Every member of [D]_I, or NULL when there are more than `max.members` of
## them. The count is taken first, so nothing is built for a class that will
## be refused; each member is then decoded directly from its index in the
## class ("unranking"), one step per member, never by generating permutations.
imec.list <- function(A, targets = list(integer(0)), max.members = Inf) {
    p <- nrow(A)
    E <- .iessgraph(A, targets)
    U <- E & t(E)
    orders <- .Call(C_cp_amo_list, U, as.double(max.members))
    if (is.null(orders)) return(NULL)
    dirpart <- E & !t(E)
    ij <- which(U & upper.tri(U), arr.ind = TRUE)
    lapply(orders, function(ord) {
        pos <- integer(p); pos[ord] <- seq_len(p)
        out <- dirpart
        if (nrow(ij))
            out[cbind(ifelse(pos[ij[, 1]] < pos[ij[, 2]], ij[, 1], ij[, 2]),
                      ifelse(pos[ij[, 1]] < pos[ij[, 2]], ij[, 2], ij[, 1]))] <- TRUE
        out
    })
}

## The enumerate-and-deduplicate listing, retained ONLY as a test reference for
## the Clique-Picking enumerator: it walks the m! vertex orderings of each
## component, so it is unusable in the search but is an independent check on
## small components.
imec.list.ref <- function(A, targets = list(integer(0)), max.size = 8L) {
    p <- nrow(A)
    E <- .iessgraph(A, targets)
    U <- E & t(E)
    comps <- .uccgs(E)
    amos <- lapply(comps, function(vs) .amos(U, vs, max.size))
    if (any(vapply(amos, is.null, TRUE))) return(NULL)
    out <- list(E & !t(E))
    for (j in seq_along(amos)) {
        nxt <- vector("list", length(out) * length(amos[[j]])); i <- 0L
        for (base in out) for (o in amos[[j]]) { i <- i + 1L; nxt[[i]] <- base | o }
        out <- nxt
    }
    out
}

## ---------------------------------------------------------------------------
## glue for the search: conversions and the bookkeeping hcmc() maintains

.adj.from.dag <- function(dag, p) {
    A <- matrix(FALSE, p, p)
    em <- graph::edgeMatrix(dag)
    if (ncol(em)) A[cbind(em[1, ], em[2, ])] <- TRUE
    A
}
.dag.from.adj <- function(A, vnames) {
    g <- graph::graphNEL(vnames, edgemode = "directed")
    E <- which(A, arr.ind = TRUE)
    if (nrow(E)) g <- graph::addEdge(vnames[E[, 1]], vnames[E[, 2]], g)
    g
}
.anc.from.adj <- function(A, vnames) {
    a <- A
    repeat { b <- a | ((a %*% a) > 0); if (all(b == a)) break; a <- b }
    dimnames(a) <- list(vnames, vnames)
    a
}
.pasets.from.adj <- function(A, vidx.nodes) {
    out <- replicate(nrow(A), integer(0), simplify = FALSE)
    for (j in seq_len(nrow(A)))
        out[[vidx.nodes[j]]] <- as.integer(vidx.nodes[which(A[, j])])
    out
}

## A uniform draw from the I-equivalence class of `dag`, returned in the same
## shape as rcar(): the DAG plus the ancestor matrix and parent sets, which are
## rebuilt because an arbitrary jump invalidates the incremental updates.
## Falls back to rcar() only when a chain component exceeds the 64 vertices the
## sampler represents as a bitmask.
isample.move <- function(dag, targets, utargets, anc, pasets, vidx, vnames,
                         vidx.nodes, r = 20L) {
    p <- length(vnames)
    A <- .adj.from.dag(dag, p)
    B <- imec.sample(A, targets)
    if (is.null(B))
        return(c(rcar(dag, r, utargets, anc, pasets, vidx), list(fallback = TRUE)))
    list(dag = .dag.from.adj(B, vnames), anc = .anc.from.adj(B, vnames),
         pasets = .pasets.from.adj(B, vidx.nodes), fallback = FALSE)
}

## Every member of the I-equivalence class of `dag`, each with its bookkeeping.
## NULL when the class has more than `max.members` of them; the count is taken
## before anything is built.
imec.members <- function(dag, targets, vnames, vidx.nodes, max.members = Inf) {
    p <- length(vnames)
    A <- .adj.from.dag(dag, p)
    L <- imec.list(A, targets, max.members)
    if (is.null(L)) return(NULL)
    lapply(L, function(B) list(dag = .dag.from.adj(B, vnames),
                               anc = .anc.from.adj(B, vnames),
                               pasets = .pasets.from.adj(B, vidx.nodes)))
}
