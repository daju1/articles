def float_formatting(float_value):
    print ("float_value =", float_value)
    decimal_part = abs(float_value - int(float_value))
    print ("decimal_part =", decimal_part)
    if 0 == decimal_part:
        digits_after_point = 0
    else:
        #print "decimal_part=", decimal_part
        print ("log(decimal_part)=", log(decimal_part))
        print ("log(10.0)=", log(10.0))
        print ("log(decimal_part)/log(10.0)=", log(decimal_part)/log(10.0))
        print ("floor(log(decimal_part)/log(10.0)=", floor(log(decimal_part)/log(10.0)))
        digits_after_point = 2 - floor(log(decimal_part)/log(10.0))
    print ("digits_after_point=", digits_after_point)
    format_string_float_value = "str_float_value = '%1." + str(digits_after_point) + "f' % (" + str(float_value) + ")"
    format_string = "{:1." + str(digits_after_point) + "f}"
    format_string = "%1." + str(digits_after_point) + "f"
    print ("format_string =", format_string)
    res = format_string.format(float_value)
    res = format_string % (float_value)
    print ("res =", res)
    return res

def suffix(t, r0, a0):
    return "_r0=" + float_formatting(r0) + "_a0=" + float_formatting(a0) + "_t=" + float_formatting(t)

def suffix_trvaR(t, r0, v0, a0, R0):
    return "_R0=" + float_formatting(R0) + "_r0=" + float_formatting(r0) + "_v0=" + float_formatting(v0) + "_a0=" + float_formatting(a0) + "_t=" + float_formatting(t)

