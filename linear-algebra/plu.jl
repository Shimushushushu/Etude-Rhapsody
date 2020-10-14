function identical_matrix(n::Integer)
    ret = zeros(n, n)
    for itr1 in 1 : n
        ret[itr1, itr1] = 1
    end
    return ret
end

function exchange_row!(mat::Matrix, a::Integer, b::Integer)
    # if $a$ or $b$ are out of bound, do nothing
    if a > size(mat, 1) || b > size(mat, 1)
        return mat
    end
    # exchange row $a$ and row $b$
    for itr1 in 1 : size(mat, 2)
        mat[a, itr1], mat[b, itr1] = mat[b, itr1], mat[a, itr1]
    end
    return mat
end

function similarity_transformation!(
        mat::Matrix, a::Integer, b::Integer)
    # if $a$ or $b$ are out of bound, do nothing
    if a > size(mat, 1) || b > size(mat, 1)
        return mat
    end
    # taking advantage of the matrix, the transformation
    # is just exchanging part of the row $a$ and row $b$
    for itr1 in 1 : (min(a, b) - 1)
        mat[a, itr1], mat[b, itr1] = mat[b, itr1], mat[a, itr1]
    end
    return mat
end

function plu(mat::Matrix)
    # make a copy of the original matrix, which will finally become $U$
    U = float(copy(mat))
    # # U = copy(mat)
    # it's assumed that the matrix has $m$ rows and $n$ columns
    m, n = size(U)
    # initialize two identical matrices as $P$ and $L$
    P = identical_matrix(m)
    L = identical_matrix(m)
    # $i$ and $j$ are used to represent iterators of row and column,
    # respectively
    i = j = 1
    while i <= m && j <= n
        # find the entry with largest absolute value
        imax, vmax = 0.0, 0.0
        for itr1 in i : m
            if vmax < abs(U[itr1, j])
                imax, vmax = itr1, abs(U[itr1, j])
            end
        end
        (vmax > 0) ? (k = imax) : (j += 1; continue)
        # left multiply $U$ and $P$ with a permutation matrix
        if k != i
            exchange_row!(U, i, k)
            exchange_row!(P, i, k)
        end
        # left and right multiply $L$ with the same matrix
        similarity_transformation!(L, i, k)
        # `forward elimination'
        for itr1 in (i + 1) : m
            scaler = U[itr1, j] / U[i, j]
            for itr2 in j : n
                U[itr1, itr2] -= U[i, itr2] * scaler
            end
            # right multiply $L$ with an elimination matrix
            L[itr1, i] = scaler
        end
        i += 1
        j += 1
    end
    return P, L, U
end