function propose_move(g)
    n = nv(g)
    while true
        u = rand(1:n)
        v = rand(1:n)
        u == v && continue 
        if !isadjacent(g, u, v) 
            rand() < 0.5 && continue
            !has_path(g, v, u) && return 1, u, v
        else 
            rand() < 0.5 && return 0, u, v
            !has_path(g, v, u) && return 2, u, v
        end
    end
end

function deltascore(g, score, mv, x, y) 
    mv == 0 && return Δscore(score, inneighbors(g, y), x, y)
    mv == 1 && return -Δscore(score, setdiff(inneighbors(g, y), x), x, y)
    mv == 2 && return -Δscore(score, setdiff(inneighbors(g, y), x), x, y) + Δscore(score, inneighbors(g, x), y, x)
end

function makemove!(g, mv, u, v) 
    mv == 0 && rem_edge!(g, u, v)
    mv == 1 && add_edge!(g, u, v)
    mv == 2 && (rem_edge!(g, u, v) || add_edge!(g, v, u))
end

function basicannealer(n, G = DiGraph(n), score)
    iterations = 1_000_000
    g = G 
    s = 0
    optg = G
    opts = s
    for iter in 1:iterations 
        temperature = 1.0 - iter/iterations 
        mv, u, v = propose_move(g)
        delta = deltascore(g, score, mv, u, v)
        if exp(delta / temperature) >= rand() 
            makemove!(g, mv, u, v)
            s += delta 
            if s > opts 
                optg = g
                opts = s
            end
        end
    end
    return g
end
