function modify_graph!(g, u, v, mod)
    if mod == :insert 
        add_edge!(g, u, v)
    end
    if mod == :delete 
        rem_edge!(g, u, v)
    end
    if mod == :turn 
        rem_edge!(g, u, v)
        add_edge!(g, v, u)
    end
end
 
function move(g, u, v, mod, local_score) 
    Δscore = -local_score(parents(g, v), v)
    mod == :turn && (Δscore -= local_score(parents(g, u), u))
    cg = copy(g)
    modify_graph!(cg, u, v, mod)
    Δscore += local_score(parents(cg, v), v)
    mod == :turn && (Δscore += local_score(parents(cg, u), u))
    return Δscore, cg 
end

function next_graph(g, local_score) 
    moves = []
    for u in vertices(g), v in vertices(g) 
        u == v && continue
        if has_edge(g, u, v)
            push!(moves, move(g, u, v, :delete, local_score))
            push!(moves, move(g, u, v, :turn, local_score))
        elseif !has_edge(g, v, u)
            push!(moves, move(g, u, v, :insert, local_score))
        end
    end
    filter!(mv -> !is_cyclic(mv[2]), moves)
    return argmax(first, moves)
end


function hill_climber(n::Integer, local_score, g=DiGraph(n))
    while true
        Δscore, nextg = next_graph(g, local_score)
        if Δscore > 10^-9
            g = nextg
        else 
            break
        end
    end
    return alt_cpdag(g)
end

function hill_climber(X::AbstractMatrix; method=:gaussian_bic, penalty=0.5, init_graph=DiGraph())
    (N, n) = size(X)
    length(vertices(init_graph)) < n && (init_graph = DiGraph(n))
    if method == :gaussian_bic
        C = Symmetric(cov(X, dims = 1, corrected = false))
        S = GaussianScore(C, N, penalty)
        return hill_climber(n, (p, v) -> local_score(S, p, v), init_graph)
    elseif method == :gaussian_bic_raw
        S = GaussianScoreQR(X, penalty)
        return hill_climber(n, (p, v) -> local_score(S, p, v), init_graph)
    else 
        throw(ArgumentError("method=$method"))
    end
end
hill_climber(X; method=:gaussian_bic, penalty=0.5, init_graph=DiGraph()) = hill_climber(Tables.matrix(X); method, penalty, init_graph)


