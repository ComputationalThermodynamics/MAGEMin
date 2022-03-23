#
# START OF PROLOGUE
#


const HASH_JEN = 0;

if isfile("libMAGEMin.dylib")
    libMAGEMin = joinpath(pwd(),"libMAGEMin.dylib")
    println("Loading local libMAGEMin dynamic library")
else
    using MAGEMin_jll
    export MAGEMin_jll
    println("Loading libMAGEMin dynamic library from MAGEMin_jll")
end




#
# END OF PROLOGUE
#