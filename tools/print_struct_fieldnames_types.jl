#=~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#   Project      : MAGEMin_C
#   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
#   Developers   : Nicolas Riel, Boris Kaus
#   Contributors : Moccetti, N. B., Dominguez, H., Assunção J., Green E., Dolejš, D., Berlie N., and Rummel L.
#   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
#   Contact      : nriel[at]uni-mainz.de
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ =#
function print_struct_fieldnames_types(obj; indent = "")
    T = typeof(obj)
    if isstructtype(T)
        println(indent, T)
        for field in fieldnames(T)
            field_type = fieldtype(T, field)
            println(indent, "├─ ", field, " :: ", field_type)
            value = getfield(obj, field)
            # Only recurse if the value is a struct (not primitive types, arrays, etc.)
            if isstructtype(field_type) && field_type != T && typeof(value) <: AbstractDict == false
                print_struct_fieldnames_types(value; indent = indent * "   ")
            end
        end
    end
end

