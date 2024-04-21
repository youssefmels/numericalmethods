import sympy as sp
import numpy as np

def is_continuous_numeric(func, interval, num_points=1000):
    x_values = np.linspace(interval[0], interval[1], num_points)
    y_values = func(x_values)

    smooth_transitions = all(np.isfinite(y_values))

    division_by_zero = not any(np.isclose(y_values, 0, atol=1e-10))
    
    print(f"Continuous: {smooth_transitions}")
    # print(f"Division by zero: {division_by_zero}")
    if isinstance(y_values, complex):
        raise ValueError("Function has negative root")

    if not smooth_transitions:
        print(f"Non-finite values at x: {x_values[np.isinf(y_values) | np.isnan(y_values)]}")


    return smooth_transitions and division_by_zero

def has_single_variable(func):
    x = sp.symbols('x')
    expression = sp.sympify(func)
    variables = sp.symbols(sp.preorder_traversal(expression))

    return (len(variables) == 1) and (x in variables)

def is_workable(func, a, b):
    if not is_continuous_numeric(func, [a, b]):
        return False

    fa, fb = func(a), func(b)

    if fa * fb > 0:
        return False

    roots_in_interval = 1 if fa * fb < 0 else 0

    print(f"Function values at endpoints: f({a}) = {fa}, f({b}) = {fb}")
    print(f"Roots in interval: {roots_in_interval}")

    return roots_in_interval == 1


def bisection_with_workability(func, a, b, d, max_iter=100):
    x = sp.symbols('x')
    if not is_workable(func, a, b):
        raise ValueError("Function is not workable in the given interval.")
    elif sp.degree(sp.sympify(func), x) > 2:
        raise ValueError("Function has more than a single root")

    approx_iterations = int(np.ceil((np.log(b - a) - np.log(0.5*10**-d) / np.log(2))-1))
    print(f"Number of Iterations= {approx_iterations}")

    iterations = 0

    prev_c_prev = None

    while iterations < max_iter and iterations < approx_iterations:
        c = (a + b) / 2
        error = 0
        if(func(c) * func(a) >= 0):
            error = abs((c-a)/c)
        else:
            error = abs((c-b)/c)
        fc = func(c)

        if iterations > 0 and prev_c_prev is not None or error <= 0.5*10**-d:
            str_c = f"{c:.15f}"
            str_prev_c_prev = f"{prev_c_prev:.15f}"
            if str_c[:d] == str_prev_c_prev[:d]:
                print(f"Converged to solution {c} in {iterations} iterations with {d} similar decimal digits.")
                return c

        if  fc < 0:
            a = c
        elif fc == 0:
            print(f"Found exact solution at {c}.")
            return c
        elif fc > 0:
            b = c

        prev_c_prev = c
        iterations += 1
        print(f"relative error is {error}")
        print(f"Iteration {iterations}: [a, b] = [{a}, {b}], c = {c}, f(c) = {fc}")

    print("Maximum number of iterations reached or achieved approximate iterations.")
    return None

equation_str = input("Enter the equation in terms of 'x': ")

x = sp.symbols('x')
user_equation = sp.sympify(equation_str)
user_function = sp.lambdify(x, user_equation, 'numpy')

interval_a = float(input("Enter the lower bound of the interval (a): "))
interval_b = float(input("Enter the upper bound of the interval (b): "))
d = int(input("Enter the number of similar digits for iterations (d): "))

try:
    result = bisection_with_workability(user_function, interval_a, interval_b, d)
    print(f"The root is: {result}")
except ValueError as e:
    print(f"Error: {e}")

