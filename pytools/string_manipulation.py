import numpy as np


def simplify_string(v):

    #print("simplify", v)

    # cast ints to ints
    if float(v).is_integer():
        s = str(int(v))
    else:
        s = str(v)

    #print("a", s)

    # drop leading zero 
    if v > 1.0: 
        if s[0] == "0" and len(s) > 1:
            s = s[1:]

    #print("b", s)

    # drop decimal separator
    s = s.replace(".", "")

    #print("c", s)

    return s

def simplify_large_num(v):

    # TODO could also fallback to regular if
    #1e-5
    #0.01
    #1e5
    #100

    s = np.format_float_scientific(v, 
                                   unique=True,
                                   trim="-",
                                   sign=False,
                                   exp_digits=False,
                                   )
    # remove + after e is there
    s = s.replace("+", "")

    return s

    # TODO manual version

    s = "{:.1e}".format(v)
    #print('simp:', s)

    # check if there are only zeros after the decimal point
    is_int = True
    for i in range(2,len(s)):
        #print('   ', i, s[i])

        if s[i] == "e":
            break

        if not( s[i] == "0"):
            is_int = False

    # strip out zeros
    if is_int:
        s = s[0] + s[i:]

    #print('simp:', s)

    # remove + after e is there
    s = s.replace("+", "")

    #found_e = False
    #for i in range(2,len(s)):
    #    if found_e:
    #    if s[i] == "e": found_e = True


    return s


