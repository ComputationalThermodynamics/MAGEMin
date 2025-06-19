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

