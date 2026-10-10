# Coordinate helpers shared by the matched-comparison bridge and its tests.
# Existing shuffling and historical observables are unchanged.

"""Physical face (x,y), in units of one primal edge, from one-based table indices."""
function comparison_face_xy(order::Integer, row::Integer, column::Integer)
    1 <= row <= 2order + 1 || throw(ArgumentError("invalid face row"))
    1 <= column <= order || throw(ArgumentError("invalid face column"))
    steps = row - 1
    return (x=column - 1 - fld(steps, 2), y=order + 1 - column - cld(steps, 2))
end

"""Edge weights in the existing Boolean domino matrix: red=a, yellow=b.

Figure 1 of Duits--Van Peski assigns random a to red and b to yellow;
blue/green have weight one. Family array indices match the creation blocks.
"""
function comparison_edge_weights(a::AbstractMatrix, b::AbstractMatrix)
    size(a) == size(b) || throw(ArgumentError("family sizes differ"))
    n, m = size(a)
    n == m || throw(ArgumentError("families must be square"))
    weights = ones(Float64, 2n, 2n)
    for i in 1:n, j in 1:n
        weights[2i, 2j-1] = a[i,j]
        weights[2i, 2j] = b[i,j]
    end
    return weights
end

"""Gamma reduction adapted to the tiling matrix's row orientation.

The paper recurrence and this repository's shuffle use opposite first-axis
orientations. Reflect input family rows into paper coordinates, apply the
recurrence, then reflect every creation table back. Without both reflections,
exact order-two Gibbs probabilities disagree. Historical routines remain
unchanged for reproducibility; new comparison code calls this adapter.
"""
function comparison_gamma_probabilities(a::AbstractMatrix, b::AbstractMatrix)
    paper = gamma_disordered_probabilities(reverse(a,dims=1),reverse(b,dims=1))
    return [reverse(level,dims=1) for level in paper]
end

"""Packed independent replica coins with the same orientation correction."""
function comparison_gamma_choice_pair(rng, a::AbstractMatrix, b::AbstractMatrix)
    paper = gamma_disordered_creation_choice_pair(rng,reverse(a,dims=1),reverse(b,dims=1))
    return (first=[reverse(level,dims=1) for level in paper.first],
            second=[reverse(level,dims=1) for level in paper.second])
end
