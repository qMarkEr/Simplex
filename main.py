import simplex as simp
from scipy import optimize

def lossFunction(eq, start_point):
        point = [0, 0]
        point[0] = eq[0] * 2 * start_point[0] + eq[2] * start_point[1] + eq[3]
        point[1] = eq[1] * 2 * start_point[1] + eq[2] * start_point[0] + eq[4]
        return point

def sub(vec1, vec2):
    res = []
    for (i, j) in zip(vec1, vec2):
        res.append(i - j)
    return res

def mul_points(pnt1, pnt2):
    return [pnt1[0] * pnt2[0], pnt1[0] * pnt2[1] + pnt1[1] * pnt2[0], pnt1[1] * pnt2[1]]

def expansion(f, eq):
    f1_t = f[0].copy()
    f1_t.insert(0, 0)
    f2_t = f[1].copy()
    f2_t.insert(0, 0)
    return [[eq[0] * i for i in mul_points(f[0], f[0])],
            [eq[1] * i for i in mul_points(f[1], f[1])],
            [eq[2] * i for i in mul_points(f[0], f[1])],
            [eq[3] * i for i in f1_t],
            [eq[4] * i for i in f2_t]]

def add(vec1, vec2):
    res = []
    for (i, j) in zip(vec1, vec2):
        res.append(i + j)
    return res

def derivative(func):
    res = func[0]
    for i in func[1::]:
        res = add(res, i)
    res[0] *= 2
    del res[-1]
    return res

def calcFunc(eq, x):
    return eq[0] * x[0]**2 + eq[1] * x[1]**2 + eq[2] * x[0] * x[1] + eq[3] * x[0] + eq[4] * x[1]

def convert_to_float(frac_str):
    try:
        return float(frac_str)
    except ValueError:
        num, denom = frac_str.split('/')
        try:
            leading, num = num.split(' ')
            whole = float(leading)
        except ValueError:
            whole = 0
        frac = float(num) / float(denom)
        return whole - frac if whole < 0 else whole + frac

def readConstrains(file):
    eq = []
    A = []
    b = []
    f = open(file)
    for i in f:
        i = i.replace(' ', '')
        if i[0] == 'F':
            eq = list(map(convert_to_float, i.split(':')[1].split(',')))
        if i[0] == 'L':
            if len(i.split(':')[1].split(',')) == 3:
                A.append(list(map(convert_to_float, i.split(':')[1].split(',')[:-1:])))
                b.append(convert_to_float(i.split(':')[1].split(',')[-1]))
            else:
                return None, None, None
    return eq, A, b

equation, A_start, b_start = readConstrains("C:/Users/Markk/Downloads/Telegram Desktop/Simplex/constrains.txt")
# x1sq x2sq x1x2 x1 x2
if equation != None and A_start != None and b_start != None:
    if len(equation) != 5 or len(A_start) != len(b_start) or len(b_start) != 2:
        print("weirdo constrains")
    else:   
        start_point = [1, 1]
        eps = 0.1
        prev_aprx = [0, 0]
        iter = 0

        print(f"x{iter} = {start_point}")
        print(f"f(x{iter}) = {calcFunc(equation, start_point)}")

        while (sub(prev_aprx, start_point)[0]**2 + sub(prev_aprx, start_point)[1]**2)**0.5 > eps:
            print("-"*20)
            iter += 1
            print("Iteration: ", iter)

            A = [row[:] for row in A_start]
            b = b_start.copy()
            c = lossFunction(equation, start_point)
            for (i, j) in zip(A, b):
                i.insert(0, j)
            c.insert(0, 0)
            simp.equalizeLessOrEqual(A, c, 2)
            basis = simp.determineBasis(A)

            res_point = simp.SolveDefault(A, c, basis, 2)

            res_point = sub(res_point, start_point)
            f_diff = [[res_point[0], start_point[0]], [res_point[1], start_point[1]]]
            exp = expansion(f_diff, equation)

            der = derivative(exp)
            l = der[1] / -der[0]

            current_aprx = add(start_point, [i * l for i in res_point])

            print(f"x{iter} = {current_aprx}")
            prev_aprx = start_point

            print(f"f(x{iter}) = {calcFunc(equation, current_aprx)}")
            start_point = current_aprx
else:
    print("wierdooo")