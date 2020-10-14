function rref(mat::Matrix)
    # make a copy of the original matrix
    ret = float(copy(mat))
    # # ret = copy(mat)
    # initialize an array to memorize the pivots
    pivot = []
    # it's assumed that the matrix has $m$ rows and $n$ columns
    m, n = size(ret)
    # $i$ and $j$ are used to represent iterators of row and column,
    # respectively
    i = j = 1
    # part of `forward elimination'
    while i <= m && j <= n
        # find the entry with largest absolute value
        imax, vmax = 0.0, 0.0
        for itr1 in i : m
            if vmax < abs(ret[itr1, j])
                imax, vmax = itr1, abs(ret[itr1, j])
            end
        end
        # but this one seems to work pretty fine?
        (vmax > 0) ? (k = imax) : (j += 1; continue)
        # # if vmax > 0
        # #     k = imax
        # # else 
        # #     j += 1
        # #     continue
        # # end
        # add the pivot position to the array
        push!(pivot, (i, j))
        # exchange the rows
        for itr1 in i : n
            ret[i, itr1], ret[k, itr1] = ret[k, itr1], ret[i, itr1]
        end
        # `forward elimination'
        for itr1 in (i + 1) : m
            scaler = ret[itr1, j] / ret[i, j]
            for itr2 in (j + 1) : n
                ret[itr1, itr2] -= ret[i, itr2] * scaler
            end
            # manually set this entry to 0, in case of using
            # LU factorization
            ret[itr1, j] = 0
        end
        i += 1
        j += 1
    end
    # part of `back subsitution'
    for itr1 in length(pivot) : -1 : 1
        # get the pivot
        i, j = pivot[itr1]
        # `back subsitution`
        for itr2 in 1 : (i - 1)
            scaler = ret[itr2, j] / ret[i, j]
            for itr3 in j : n
                ret[itr2, itr3] -= ret[i, itr3] * scaler
            end
        end
        # make entry on the pivot be 1
        # and a little trick in order not to create a new variable
        for itr2 in n : -1 : j
            ret[i, itr2] /= ret[i, j]
        end
    end
    return ret
end