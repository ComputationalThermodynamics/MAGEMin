#=~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#   Project      : MAGEMin_C
#   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
#   Developers   : Nicolas Riel, Boris Kaus
#   Contributors : Dominguez, H., Assunção J., Green E., Berlie N., and Rummel L.
#   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
#   Contact      : nriel[at]uni-mainz.de
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ =#
# Set of routines taken from MAGEMinApp to perform adaptive mesh refinement (AMR) on a 2D grid. 
# The AMR is performed by splitting the cells when a reaction line is detected.


"""
    AMR_data

    Mutable structure holding the state of an Adaptive Mesh Refinement (AMR) grid.

    Fields
    ------
    cells : Vector{Vector{Int64}}
        Current quadrilateral cells, each defined by 4 point indices (counter-clockwise).
    ncells : Vector{Vector{Int64}}
        Newly generated cells from the last refinement step.
    points : Vector{Vector{Float64}}
        All grid point coordinates [X, Y].
    npoints : Vector{Vector{Float64}}
        Newly generated point coordinates from the last refinement step.
    npoints_ig : Vector{Tuple}
        Parent point index tuples for each new point (used for initial guess interpolation).
    hash_map : Dict{Vector{Float64}, Int}
        Map from point coordinates to their index, used to avoid duplicate points.
    bnd_cells : Vector{Tuple}
        Boundary information for kept cells: tuples of (cell_index, midpoint_indices...).
    split_cell_list : Vector{Int64}
        Indices of cells to be split in the next refinement pass.
    keep_cell_list : Vector{Int64}
        Indices of cells to be kept unchanged in the next refinement pass.
    Xrange : Vector{Float64}
        X-axis bounds of the domain [Xmin, Xmax].
    Yrange : Vector{Float64}
        Y-axis bounds of the domain [Ymin, Ymax].
"""
mutable struct AMR_data
    cells           :: Vector{Vector{Int64}}
    ncells          :: Vector{Vector{Int64}}
    points          :: Vector{Vector{Float64}}
    npoints         :: Vector{Vector{Float64}}
    npoints_ig      :: Vector{Tuple}
    hash_map        :: Dict{Vector{Float64}, Int}
    bnd_cells       :: Vector{Tuple}
    split_cell_list :: Vector{Int64}
    keep_cell_list  :: Vector{Int64}
    Xrange          :: Vector{Float64}
    Yrange          :: Vector{Float64}
end

"""
    compute_index(value, min_value, delta)

    Compute the 1-based grid index of a coordinate value along one axis.

    Parameters
    ----------
    value : Float64
        Coordinate value to index.
    min_value : Float64
        Minimum value of the axis.
    delta : Float64
        Grid spacing along the axis.

    Returns
    -------
    index : Int64
        1-based index of `value` in the uniform grid.
"""
function compute_index(value, min_value, delta)
    return Int64(round((value - min_value + delta) / delta))
end


"""
    initialize_AMR(Xrange, Yrange, igs)

    Create an initial uniform quadrilateral grid for Adaptive Mesh Refinement.

    The grid has `2^igs` cells per dimension, covering the rectangular domain defined by `Xrange` and `Yrange`.

    Parameters
    ----------
    Xrange : Union{Tuple{Float64,Float64}, Vector{Float64}}
        X-axis bounds of the domain [Xmin, Xmax].
    Yrange : Union{Tuple{Float64,Float64}, Vector{Float64}}
        Y-axis bounds of the domain [Ymin, Ymax].
    igs : Int64
        Initial subdivision level; produces `2^igs × 2^igs` cells.

    Returns
    -------
    data : AMR_data
        Initialized AMR data structure with the uniform grid.
"""
function initialize_AMR(Xrange,Yrange,igs)

    nc_per_dim      = 2^igs              #number of cells per dimension
    np_per_dim      = nc_per_dim+1       #number of points per dimension

    points          = Vector{Vector{Float64}}(undef, 0)
    cells           = Vector{Vector{Int64}}(undef, 0)
    npoints         = Vector{Vector{Float64}}(undef, 0)
    npoints_ig      = Vector{Tuple}(undef, 0)   
    ncells          = Vector{Vector{Int64}}(undef, 0)
    hash_map        = Dict{Vector{Float64}, Int}()
    bnd_cells       = Vector{Tuple}(undef, 0)
    split_cell_list = Vector{Int64}(undef, 0)
    keep_cell_list  = Vector{Int64}(undef, 0)

    # initialize rectinilinear grid
    for j in 1:np_per_dim
        for i in 1:np_per_dim
            push!(points, [(j-1)*(Xrange[2]-Xrange[1])/(np_per_dim-1) + Xrange[1], (i-1)*(Yrange[2]-Yrange[1])/(np_per_dim-1) + Yrange[1]])
        end
    end
    for j in 1:nc_per_dim
        for i in 1:nc_per_dim
            push!(cells, [i+(j-1)*np_per_dim, i+1+(j-1)*np_per_dim, i+1+j*np_per_dim, i+j*np_per_dim])
        end
    end

    data = AMR_data(    cells,
                        ncells,
                        points,
                        npoints,
                        npoints_ig,
                        hash_map,
                        bnd_cells,
                        split_cell_list,
                        keep_cell_list,
                        [Xrange[1],Xrange[2]],[Yrange[1],Yrange[2]])
    return data
end

"""
    all_identical(arr::Vector{UInt64})

    Check whether all elements in a `UInt64` vector are identical.

    Parameters
    ----------
    arr : Vector{UInt64}
        Vector of hash values to compare.

    Returns
    -------
    result : Bool
        `true` if all elements equal `arr[1]`, `false` otherwise.
"""
function all_identical(arr::Vector{UInt64})
    return all(x -> x == arr[1], arr)
end


"""
    split_and_keep(data, Hash_XY)

    Categorise cells into the split or keep list based on phase assembly uniformity.

    For each cell in `data.keep_cell_list` (from the previous refinement step) and
    any remaining cells beyond that list, the function collects the phase-assembly
    hash values at all corner points (and any boundary midpoints already stored in
    `data.bnd_cells`). If all hashes are identical the cell is stable and added to
    `keep_cell_list`; otherwise it straddles a reaction boundary and is added to
    `split_cell_list`.

    Parameters
    ----------
    data : AMR_data
        AMR state structure. `data.cells`, `data.keep_cell_list`, and
        `data.bnd_cells` are read; `data.split_cell_list` and
        `data.keep_cell_list` are overwritten.
    Hash_XY : Vector{UInt64}
        Per-point phase-assembly hash values, indexed by point index.

    Returns
    -------
    data : AMR_data
        Updated AMR state with refreshed `split_cell_list` and `keep_cell_list`.
"""
function split_and_keep(data,Hash_XY)
    kp                      = length(data.keep_cell_list)
    data.split_cell_list    = []
    data.keep_cell_list     = []
    tot                     = length(data.cells)

    for i=1:kp
        tmp0 = Vector{UInt64}(undef, 0)
        for j=1:4
            tmp0 = push!(tmp0,Hash_XY[data.cells[i][j]])
        end
        nb = length(data.bnd_cells[i])
        if nb > 1
            for j=2:nb
                tmp0 = push!(tmp0,Hash_XY[data.bnd_cells[i][j]])
            end
        end
        if all_identical(tmp0)
            push!(data.keep_cell_list, i)
        else
            push!(data.split_cell_list, i)
        end

    end

    tmp = Vector{UInt64}(undef, 4)
    for i=kp+1:tot
        for j=1:4
            tmp[j] = Hash_XY[data.cells[i][j]]
        end

        if all_identical(tmp)
            push!(data.keep_cell_list, i)
        else
            push!(data.split_cell_list, i)
        end
    end

    return data
end

"""
    AMR(data)

    Perform one adaptive mesh refinement step on the quadrilateral grid.

    Each cell in `data.split_cell_list` is subdivided into four child cells by
    inserting five new points: one at the cell centroid and one at each edge
    midpoint. Duplicate midpoints shared by adjacent cells are detected via
    `data.hash_map` and reused rather than duplicated. Cells in
    `data.keep_cell_list` are retained unchanged; their edge midpoints that were
    created by neighbouring splits are recorded in `data.bnd_cells` for use in
    the next `split_and_keep` call.

    Parameters
    ----------
    data : AMR_data
        AMR state structure. `data.cells`, `data.points`, `data.hash_map`,
        `data.split_cell_list`, and `data.keep_cell_list` are read and updated.

    Returns
    -------
    data : AMR_data
        Updated AMR state with:
        - `data.cells`     — combined list of kept cells followed by new child cells.
        - `data.points`    — extended point list including all new midpoints.
        - `data.ncells`    — child cells generated in this step.
        - `data.npoints`   — new point coordinates generated in this step.
        - `data.npoints_ig`— parent-point index tuples for each new point, used for
                             initial-guess interpolation.
        - `data.bnd_cells` — boundary midpoint records for each kept cell.
        - `data.hash_map`  — updated coordinate-to-index map.
"""
function AMR(data)
    npoints         = Vector{Vector{Float64}}(undef, 0)
    npoints_ig      = Vector{Tuple}(undef, 0)
    ncells          = Vector{Vector{Int64}}(undef, 0)

    tp              = length(data.points)
    ns              = length(data.split_cell_list)
    for i=1:ns

        tmp     = data.points[data.cells[data.split_cell_list[i]][1]]/2.0 + data.points[data.cells[data.split_cell_list[i]][3]]/2.0
        if haskey(data.hash_map, tmp)
            p = data.hash_map[tmp]
        else
            push!(npoints, tmp)
            push!(npoints_ig, (data.cells[data.split_cell_list[i]][1], data.cells[data.split_cell_list[i]][2], data.cells[data.split_cell_list[i]][3], data.cells[data.split_cell_list[i]][4]))
            
            tp += 1
            data.hash_map[tmp] = tp
            p  = tp
        end
        c       = p;

        tmp = data.points[data.cells[data.split_cell_list[i]][1]]/2.0 + data.points[data.cells[data.split_cell_list[i]][2]]/2.0
        if haskey(data.hash_map, tmp)
            p = data.hash_map[tmp]
        else
            push!(npoints, tmp)
            push!(npoints_ig, (data.cells[data.split_cell_list[i]][1], data.cells[data.split_cell_list[i]][2]))
            
            tp += 1
            data.hash_map[tmp] = tp
            p  = tp
        end
        w       = p;

        tmp = data.points[data.cells[data.split_cell_list[i]][2]]/2.0 + data.points[data.cells[data.split_cell_list[i]][3]]/2.0
        if haskey(data.hash_map, tmp)
            p = data.hash_map[tmp]
        else
            push!(npoints, tmp)
            push!(npoints_ig, (data.cells[data.split_cell_list[i]][2], data.cells[data.split_cell_list[i]][3]))
           
            tp += 1
            data.hash_map[tmp] = tp
            p  = tp
        end
        n       = p;

        tmp = data.points[data.cells[data.split_cell_list[i]][3]]/2.0 + data.points[data.cells[data.split_cell_list[i]][4]]/2.0
        if haskey(data.hash_map, tmp)
            p = data.hash_map[tmp]
        else
            push!(npoints, tmp)
            push!(npoints_ig, (data.cells[data.split_cell_list[i]][3], data.cells[data.split_cell_list[i]][4]))
            
            tp += 1
            data.hash_map[tmp] = tp
            p  = tp
        end
        e       = p;

        tmp = data.points[data.cells[data.split_cell_list[i]][4]]/2.0 + data.points[data.cells[data.split_cell_list[i]][1]]/2.0
        if haskey(data.hash_map, tmp)
            p = data.hash_map[tmp]
        else
            push!(npoints, tmp)
            push!(npoints_ig, (data.cells[data.split_cell_list[i]][4], data.cells[data.split_cell_list[i]][1]))
            
            tp += 1
            data.hash_map[tmp] = tp
            p  = tp
        end
        s       = p;

        push!(ncells, [data.cells[data.split_cell_list[i]][1], w, c, s])
        push!(ncells, [w, data.cells[data.split_cell_list[i]][2], n, c])
        push!(ncells, [c, n, data.cells[data.split_cell_list[i]][3], e])
        push!(ncells, [s, c, e, data.cells[data.split_cell_list[i]][4]])
    end

    data.points     = vcat(data.points, npoints)

    data.bnd_cells  = Vector{Tuple}(undef, 0)
    nk              = length(data.keep_cell_list) 
    ix = [1 2; 2 3; 3 4; 4 1]
    for i=1:nk
        tmp_bnd = (i,)
        for k = 1:4
            tmp = data.points[data.cells[data.keep_cell_list[i]][ix[k,1]]]/2.0 + data.points[data.cells[data.keep_cell_list[i]][ix[k,2]]]/2.0
            if haskey(data.hash_map, tmp)
                p        = data.hash_map[tmp]
                tmp_bnd = (tmp_bnd...,p)
            end
        end

        push!(data.bnd_cells,tmp_bnd)

    end

    data.cells = vcat(data.cells[data.keep_cell_list], ncells)

    data.ncells     = ncells
    data.npoints    = npoints
    data.npoints_ig = npoints_ig

    return data
end