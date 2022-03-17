# This installs the precompiled MAGEMin_jll binary to the current working directory
import Pkg
Pkg.add("MAGEMin_jll")

# Use MAGEMin_jll
using MAGEMin_jll

function write_environmental_variables_file()

    open("environmental_variables.m", "w") do io
        pth = MAGEMin_jll.__init__();           # location of required dynamic libraries
        write(io, "path_dylib = '$(pth)'; \n")

        pth = MAGEMin_jll.PATH[]                # location of binaries
        write(io, "path_bin = '$(pth)'; \n")
        
        pth = Sys.BINDIR                # location of binaries
        write(io, "path_julia = '$(pth)'; \n")
        
    end
    return nothing
end

write_environmental_variables_file()

println("Succesfully downloaded the MAGEMin_jll package")
println("And installed the environmental_variables.m file in: $(pwd())")

# Show the version number
run(`$(MAGEMin()) --version`)
