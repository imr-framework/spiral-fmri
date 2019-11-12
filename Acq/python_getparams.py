def struct2python(params_MATLAB_struct):
    struct = params_MATLAB_struct[0,0]
    fieldNames = struct.dtype.fields.keys()
    output = {}
    for key in fieldNames:
        struct_el = struct[key]
        is_nested = len(struct_el.dtype)
        if is_nested == 0:
            while struct_el.dtype == object:
                struct_el = struct_el[0]
            output.update({key: struct_el[0]})
        else:
            nested_struct_fieldNames = struct_el.dtype.fields.keys()
            nested_output = {}
            for key2 in nested_struct_fieldNames:
                nested_struct_el = struct_el[key2]
                while nested_struct_el.dtype == object:
                    nested_struct_el = nested_struct_el[0]
                nested_output.update({key2: nested_struct_el[0]})
            output.update({key: nested_output})
    params_python_struct = output

    return params_python_struct