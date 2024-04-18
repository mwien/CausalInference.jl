function propose_move(g)
    n = nv(g)
    while true
        u = rand(1:n)
        v = rand(1:n)
        u == v && continue 
        if !isadjacent(g, u, v) 
            rand() < 0.5 && !has_path(g, u, v) && return 1, v, u # insert v -> u 
            !has_path(g, v, u) && return 1, u, v # insert u -> v
        else 
            if has_edge(g, v, u) 
                tmp = u
                u = v 
                v = tmp
            end
            rand() < 0.5 && return 0, u, v # delete u -> v
            !has_a_path(g, setdiff(outneighbors(g, u), v), v) && return 2, u, v # reverse u -> v
        end
    end
end

function deltascore(g, score, mv, u, v) 
    mv == 0 && return -Δscore(score, setdiff(inneighbors(g, v), u), u, v)
    mv == 1 && return Δscore(score, inneighbors(g, v), u, v)
    mv == 2 && return -Δscore(score, setdiff(inneighbors(g, v), u), u, v) + Δscore(score, inneighbors(g, u), v, u)
end

function makemove!(g, mv, u, v) 
    mv == 0 && rem_edge!(g, u, v)
    mv == 1 && add_edge!(g, u, v)
    if mv == 2 
        rem_edge!(g, u, v) 
        add_edge!(g, v, u)
    end
end

function basicannealer(N, G, score, iterations)
    g = G 
    s = 0.0
    optg = G
    opts = s
    for iter in 1:iterations 
        temperature = (1.0 - iter/iterations) * 10
        mv, u, v = propose_move(g)
        delta = deltascore(g, score, mv, u, v)
        if exp(delta / temperature) >= rand() 
            makemove!(g, mv, u, v)
            s += delta 
            if s > opts
                optg = copy(g)
                opts = s
            end
        end
    end
    return alt_cpdag(optg)
end

function basicannealer(X::AbstractMatrix; method=:gaussian_bic, penalty=0.5, iterations=1_000_000)
    (N, n) = size(X)
    n ≥ N && @warn "High dimensional data (n ≤ p), ges might not terminate."
    if method == :gaussian_bic
        C = Symmetric(cov(X, dims = 1, corrected = false))
        S = GaussianScore(C, N, penalty)
        return basicannealer(N, DiGraph(n), S, iterations)
    elseif method == :gaussian_bic_raw
        S = GaussianScoreQR(X, penalty)
        return basicannealer(N, DiGraph(n), S, iterations)
    else 
        throw(ArgumentError("method=$method"))
    end
end

basicannealer(X; method=:gaussian_bic, penalty=0.5, iterations = 1_000_000) = basicannealer(Tables.matrix(X); method, penalty, iterations)

function preorientannealer(N, G, score, iterations, preorientations)
    g = G 
    s = 0.0
    optg = G
    opts = s
    iter = 1
    while true
        iter > iterations && break
        temperature = (1.0 - iter/iterations) * 10
        mv, u, v = propose_move(g)
        mv == 1 && preorientations[u, v] == 0 && continue 
        mv == 2 && preorientations[v, u] == 0 && continue
        iter += 1
        delta = deltascore(g, score, mv, u, v)
        if exp(delta / temperature) >= rand() 
            makemove!(g, mv, u, v)
            s += delta 
            if s > opts
                optg = copy(g)
                opts = s
            end
        end
    end
    return alt_cpdag(optg)
end

function preorientannealer(X::AbstractMatrix, preorientations; method=:gaussian_bic, penalty=0.5, iterations=1_000_000)
    (N, n) = size(X)
    n ≥ N && @warn "High dimensional data (n ≤ p), ges might not terminate."
    if method == :gaussian_bic
        C = Symmetric(cov(X, dims = 1, corrected = false))
        S = GaussianScore(C, N, penalty)
        return preorientannealer(N, DiGraph(n), S, iterations, preorientations)
    elseif method == :gaussian_bic_raw
        S = GaussianScoreQR(X, penalty)
        return preorientannealer(N, DiGraph(n), S, iterations, preorientations)
    else 
        throw(ArgumentError("method=$method"))
    end
end

preorientannealer(X, preorientations; method=:gaussian_bic, penalty=0.5, iterations = 1_000_000) = preorientannealer(Tables.matrix(X), preorientations; method, penalty, iterations)


