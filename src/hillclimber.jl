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
    #println(moves)
    filter!(mv -> !is_cyclic(mv[2]), moves)
    #println(moves)
    return argmax(first, moves)
end


# also call with complete DAG and random edge orientations
function hill_climber(n, local_score, g=DiGraph(n)) 
    while true
        Δscore, nextg = next_graph(g, local_score)
        if Δscore > 10^-9
            g = nextg
        else 
            break
        end
    end
    # println(length(edges(g)))
    return g
end

function hill_climber(X::AbstractMatrix; method=:gaussian_bic, penalty=0.5)
    (_, n) = size(X)
    if method == :gaussian_bic
        C = Symmetric(cov(X, dims = 1, corrected = false))
        S = GaussianScore(C, n, penalty)
        return hill_climber(n, (p, v) -> local_score(S, p, v))
    elseif method == :gaussian_bic_raw
        S = GaussianScoreQR(X, penalty)
        return hill_climber(n, (p, v) -> local_score(S, p, v))
    else 
        throw(ArgumentError("method=$method"))
    end
end
hill_climber(X; method=:gaussian_bic, penalty=0.5) = hill_climber(Tables.matrix(X); method, penalty)
