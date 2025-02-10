def print_stack(stack):
    print ("print_stack")
    for frame in stack:
        file = frame[1]
        line = frame[2]
        func_name = frame[3]
        code_context = frame[4]
        print ("file = ", file)
        print ("line = ", line)
        print ("func_name = ", func_name)
        print ("code_context = ", code_context)

def get_integrand_view(f):
    return f(x)

class my_dummy_integral:
    f = None
    a = None
    b = None
    def __init__(self, f, a, b):
        print ("my_dummy_integral ", f, a, b)
        self.f = f
        self.a = a
        self.b = b

def num_int(f, a, b):
    from scipy import integrate

    log_fn = 'num_int.txt'

    to_call_integration = True

    if type(f) is my_dummy_integral:
        to_call_integration = False

    import inspect
    stack = inspect.stack()
    for frame in stack:
        func_name = frame[3]
        # print ("func_name = ", func_name)
        if ('get_integrand_view' == func_name):
            to_call_integration = False
            break;

    try:
        print ("integrand = ", get_integrand_view(f))
    except Exception as ex:
        print ("Exception while print get_integrand_view ex = ", ex)
    print ("a = ", a)
    print ("b = ", b)
    print ("to_call_integration = ", to_call_integration)

    # file = open('field_of_deyna_cylinder_num_int.txt', 'a')
    file = open(log_fn, 'a')
    file.write('\n')

    file.write("f = " + str(f))
    file.write('\n')
    try:
        file.write("integrand = " + str(get_integrand_view(f)))
    except Exception as ex:
        file.write("Exception while write get_integrand_view ex = " + str(ex))
    file.write('\n')
    file.write("a = " + str(a) + ", b = " + str(b))
    file.write('\n')
    file.write("to_call_integration = " + str(to_call_integration))
    file.write('\n\n')
    file.close()

    if not to_call_integration:
        print ("")
        return my_dummy_integral(f,a,b)

    try:
        integral = integrate.quad(f, a, b)
        # integral = numerical_integral(f, a, b)

        print ("integral = ", integral)
        print ("")

        # file = open(log_fn, 'a')
        file = open(log_fn, 'a')
        file.write('\n')
        file.write("integral = " + str(integral))
        file.write('\n\n')
        file.close()

        result = integral[0]
        return result

    except Exception as ex:

        print ("Exception ex = ", str(ex))
        # print ("f = ", f)
        if to_call_integration:
            try:
                print ("integrand = ", get_integrand_view(f))
            except Exception as ex2:
                print ("Exception while print get_integrand_view ex2 = ", ex2)
            print ("a = ", a)
            print ("b = ", b)
            stack = inspect.stack()
            print_stack(stack)

        file = open(log_fn, 'a')
        file.write('\n')
        file.write("Exception ex = " + str(ex))
        file.write('\n')
        file.write("f = " + str(f))
        file.write('\n')
        try:
            file.write("integrand = " + str(get_integrand_view(f)))
        except Exception as ex2:
            file.write("Exception while write get_integrand_view ex2 = " + str(ex2))
        file.write('\n')
        file.write("a = " + str(a) + ", b = " + str(b))
        file.write('\n\n')
        file.write('\n\n')
        file.close()

        if 'unable to simplify to float approximation' == str(ex):
            raise ex

        raise ex

        integral = numerical_integral(f, a, b)

        print ("integral = ", integral)

        file = open(log_fn, 'a')
        file.write('\n')
        file.write("integral = " + str(integral))
        file.write('\n\n')
        file.close()

        result = integral[0]
        print ("result = ", result)
        return result

