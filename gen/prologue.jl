#
# START OF PROLOGUE
#
using MAGEMin_C, MAGEMin_jll
const HASH_JEN = 0;


function __init__()
    if isfile("libMAGEMin.dylib")
        global libMAGEMin = joinpath(pwd(),"libMAGEMin.dylib")
        println("Using locally compiled version of libMAGEMin.dylib")
    else
        global libMAGEMin = MAGEMin_jll.libMAGEMin
        println("Using libMAGEMin.dylib from MAGEMin_jll")
    end
end

#
# END OF PROLOGUE
#