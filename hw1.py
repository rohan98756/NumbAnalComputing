import sympy as sp

x = sp.symbols('x')

equations_dict = {}

def load_equations(filename):
    with open(filename, 'r') as file:   
        for line in file:
            parts = line.strip().split(', ')
            equation = parts[0]  
            x0_part = parts[1]   
            a_part = parts[2]
            b_part = parts[3]     
            
            x0 = float(eval(x0_part.split('=')[1].strip()))
            a = float(eval(a_part.split('=')[1].strip()))
            b = float(eval(b_part.split('=')[1].strip()))
            
            equations_dict[equation.strip()] = {'x0': x0, 'a': a , 'b': b}

    return equations_dict

def f_of_x(equation, x_value):
    expr = sp.sympify(equation)
    evaluated_expr = expr.evalf(subs={x: x_value})
    return evaluated_expr

def fprime_of_x(equation, x_value):
    expr = sp.sympify(equation)
    derivative = sp.diff(expr, x)

    return derivative.evalf(subs={x: x_value})

def newtons_method(equation, intial_guess):
    newton_count = 0
    x0=intial_guess

    eq = sp.sympify(equation)
    derv = fprime_of_x(eq, intial_guess)

    if(abs(derv)<(1e-6)): 
        print("Derivative is zero at the initial guess. Root finding may fail.")
        return None, None
    
    x1 =  x0 - (f_of_x(eq, x0)/derv)
  
    while(abs(x1-x0) >= (1e-6)):
        x0 = x1 
        derv = fprime_of_x(eq, x1)

        if(abs(derv)<(1e-6)): 
            print("Newton's Method failed to converge from the given initial guess.")
            return None, None
        
        x1 = x1 - (f_of_x(eq, x1)/derv)
        newton_count+=1

    return x1,newton_count 

def secant_method(eq, initial_guess):
    secant_count = 0
    x1 = initial_guess
    x2 = initial_guess + 0.1

    while True:

        denom = (f_of_x(eq, x2) - f_of_x(eq, x1))

        if(abs(denom) < 1e-6):
            print("Division by zero. Secant method failed to converge.")
            return None, None
        
        xnext = x2 - (f_of_x(eq, x2)) * ((x2 - x1) / denom)
        secant_count+=1

        if(abs(xnext-x2) < (1e-6)): 
            break

        x1 = x2 
        x2 = xnext
           
    return xnext, secant_count

def bisection_method(eq, a, b):
    bisection_count=0
    
    # For my own reference: a < b and f(a)* f(b) < 0

    if f_of_x(eq, a) * f_of_x(eq, b) >= 0:
        print("Values of a and b are not valid")
        return None, None

    while abs(b-a) >= 1e-6:

        t = (a+b)/2

        if(f_of_x(eq, t) == 0.0): break
        
        if(f_of_x(eq, t) * f_of_x(eq,a) < 0):
            b=t
        else:
            a=t
        bisection_count+=1

    return (a+b)/2, bisection_count
       
load_equations("equations.txt")

for equation, values in equations_dict.items():
    x0 = values['x0']  
    a = values['a']
    b = values['b']    
    eq_sym = sp.sympify(equation)

    print("=" * 60)
    print(f"Equation: {sp.latex(eq_sym)}")
    print(f"Value of x0: {x0}", "\n")

    # I hardcoded values of a and b by graphing each function on Desmos
    print(f"Value of a: {a}")
    print(f"Value of b: {b}" , "\n")

    newton_root, newton_count = newtons_method(equation, x0)
    sec_root, sec_count = secant_method(equation, x0)
    bisection_root,  bisection_count = bisection_method(equation, a, b)

   # print(f"Newton's method root is: {newton_root}")
    #print(f"Number of iterations: {newton_count}")
    print(f"Secant's method root is: {sec_root}")
    print(f"Number of iterations: {sec_count}")
    print(f"Bisection's method root is: {bisection_root}")
    print(f"Number of iterations: {bisection_count}")

    print("=" * 60, "\n")


