def float_formatting(float_value):
    decimal_part = abs(float_value - int(float_value))
    if 0 == decimal_part:
        digits_after_point = 1
    else:
        #print "decimal_part=", decimal_part
        #print "log(decimal_part)/log(10.0)=", log(decimal_part)/log(10.0)
        digits_after_point = 1 - floor(log(decimal_part)/log(10.0))
        #print "digits_after_point=", digits_after_point
    format_string_float_value = "\"str_float_value = '%1." + str(digits_after_point) + "f' % (" + str(float_value) + ")\""
    #format_string_float_value = "str_float_value = {:1." + str(digits_after_point) + "f}.format(" + str(float_value) + ")"
    exec(format_string_float_value)
    try:
        return str_float_value
    except:
        return str(float_value)

def suffix(t, r0, a0):
    return "_t=" + float_formatting(t) + "_r0=" + float_formatting(r0) + "_a0=" + float_formatting(a0)

