function identical_matrix(n::Integer)
    ret = zeros(n, n)
    for itr1 in 1 : n
        ret[itr1, itr1] = 1
    end
    return ret
end

function inner_product(mat::Matrix, a::Integer, b::Integer)
    # if $a$ or $b$ are out of bound, return $0$
    if a > size(mat, 2) || b > size(mat, 2)
        return 0
    end
    # initialize the result
    ret = 0
    # calculate the inner (dot) product of column $a$ and $b$
    for itr1 in 1 : size(mat, 1)
        ret += mat[itr1, a] * mat[itr1, b]
    end
    return ret
end

function column_elimination!(mat::Matrix, a::Integer, b::Integer, k::Number, row::Integer = 1)
    for itr1 in row : size(mat, 1)
        mat[itr1, b] += mat[itr1, a] * k
    end
    return mat
end

function mgs(mat::Matrix)
    # make a copy of the matrix as $Q$
    Q = float(copy(mat))
    # it's assumed that the matrix has $m$ rows and $n$ columns
    m, n = size(mat)
    # make a identical matrix as $R$
    R = identical_matrix(n)
    for itr1 in 1 : n
        # calculate the norm of the column vector
        norm = sqrt(inner_product(Q, itr1, itr1))
        # if the column is $0$, just skip it
        norm == 0 && continue
        # normalize the column vector
        for itr2 in 1 : m
            Q[itr2, itr1] /= norm
        end
        # multiply the $itr1$th row with the norm
        for itr2 in itr1 : n
            R[itr1, itr2] *= norm
        end
        # substract the projections of the current vector
        # from the remaining vectors
        for itr2 in (itr1 + 1) : n
            scaler = inner_product(Q, itr1, itr2)
            column_elimination!(Q, itr1, itr2, -scaler)
            # add the corresponding vectors to some rows
            for itr3 in itr2 : n
                R[itr1, itr3] += R[itr2, itr3] * scaler
            end
        end
    end
    return Q, R
end