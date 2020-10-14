function cholesky(mat::Matrix)
    # if the matrix is no a square matrix, return one with zeros
    if size(mat, 1) != size(mat, 2)
        return zeros(size(mat, 1), size(mat, 1))
    end
    # since verifing that the matrix is positive definite is so
    # time-consuming, the step is skipped
    # it's assumed that the matrix has $n$ rows and columns
    n = size(mat, 1)
    # $L$ is the matrix to be returned
    L = zeros(n, n)
    # the formula is proved by induction, basically
    for itr1 in 1 : n, itr2 in itr1 : n
        if itr2 == itr1
            L[itr2, itr1] = mat[itr2, itr1]
            for itr3 in 1 : (itr1 - 1)
                L[itr2, itr1] -= L[itr2, itr3] ^ 2
            end
            L[itr2, itr1] = sqrt(L[itr2, itr1])
        else
            L[itr2, itr1] = mat[itr2, itr1]
            for itr3 in 1 : (itr1 - 1)
                L[itr2, itr1] -= L[itr2, itr3] * L[itr1, itr3]
            end
            L[itr2, itr1] /= L[itr1, itr1]
        end
    end
    return L
end