import math
import matplotlib.pyplot as plt

# f(x) = e^x
def f_of_x(x_value):
    return math.exp(x_value)

def newton_interpolation(n):
    num_of_points = n + 1
    x_values = [-1 + (2 * i) / n for i in range(num_of_points)] 
    y_values = [f_of_x(x) for x in x_values]


    divided_difference_table = [[0 for _ in range(num_of_points)] for _ in range(num_of_points)]

 
    for i in range(num_of_points):
        divided_difference_table[i][0] = y_values[i]

    for j in range(1, num_of_points):
        for i in range(num_of_points - j):
            divided_difference_table[i][j] = (
                (divided_difference_table[i + 1][j - 1] - divided_difference_table[i][j - 1]) /
                (x_values[i + j] - x_values[i])
            )

    return [divided_difference_table[0][i] for i in range(num_of_points)]  

def polynomial(n, x, coeffs, x_nodes):
    result = coeffs[0]
    for i in range(1, len(coeffs)):
        term = coeffs[i]
        for j in range(i):
            term *= (x - x_nodes[j])
        result += term
    return result

def maximum_error(coeffs, n):
    max_error = 0
    x_nodes = [-1 + (2 * k) / n for k in range(n + 1)]

    for k in range(501):
        tk = -1 + ((2 * k) / 500)
        error = abs(f_of_x(tk) - polynomial(n, tk, coeffs, x_nodes))
        max_error = max(max_error, error)

    return max_error

n_values = [2, 4, 8, 16, 32]
points = []

for n in n_values:
    coeffs = newton_interpolation(n)
    error = maximum_error(coeffs, n)
    points.append((n, error))

for point in points:
    x_value, y_value = point
    print(f"The maximum error for the n={x_value} degree Newton interpolant is: {y_value:.15e}")

x_coords, y_coords = zip(*points)

# Log scale
plt.scatter(x_coords, y_coords, color='blue', label='Error')
plt.plot(x_coords, y_coords, color='orange')

plt.title('log scale - Maximum Error versus n')
plt.xlabel('n')
plt.ylabel('Maximum Error')
plt.yscale('log')  
plt.grid()
plt.legend()
plt.show()

# Regular scale
plt.scatter(x_coords, y_coords, color='blue', label='Error')
plt.plot(x_coords, y_coords, color='orange')

plt.title('regular scale - Maximum Error versus n')
plt.xlabel('n')
plt.ylabel('Maximum Error')
plt.grid()
plt.legend()
plt.show()





