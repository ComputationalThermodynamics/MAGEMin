# NR 12-09-24 Script to test nullspace minimization with generic solution models 

"""
    eye(i::Int64, j::Int64)
"""
function eye(i::Int64, j::Int64)
    return i == j ? 1.0 : 0.0
end


"""
    G_excess( W, v, p, symmetry )
"""
function get_G_ex(  W           :: Vector{Float64},
                    v           :: Vector{Float64}, 
                    p           :: Vector{Float64}, 
                    symmetry    :: Int64                )
    
    if symmetry == 1
        G_ex  = 0.0
        it      = 1
        for j in 1:n_em-1
            for k in j+1:n_em
                G_ex += (eye(i,j) - p[j]) * (eye(i,k) - p[k]) * W[it]
                it += 1
            end
        end
    else
        sum_v = 0.0
        for i in 1:n_em
            sum_v += p[i] * v[i]
        end

        mat_phi = zeros(n_em)
        for i in 1:n_em
            mat_phi[i] = (p[i] * v[i]) / sum_v
        end

        G_ex  = 0.0
        it      = 1
        for j in 1:n_em-1
            for k in j+1:n_em
                excess += (eye(i,j) - mat_phi[j]) * (eye(i,k) - mat_phi[k]) * (W[it] * 2.0 * v[i] / (v[j] + v[k]))
                it += 1
            end
        end        
    end

    return G_ex
end
