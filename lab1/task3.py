import numpy as np

x_leq = lambda x, y, a, b: x < 1 - ((b - y) ** 2) / (4 * (1 - a))

y_leq = lambda x, y, a, b: y < 1 - ((a - x) ** 2) / (4 * (1 - b))

x_geq = lambda x, y, a, b: x > ((b - y) ** 2) / (4 * a)

y_geq = lambda x, y, a, b: y > ((a - x) ** 2) / (4 * b)


hit = lambda x, y, a, b: x_leq(x, y, a, b) and y_leq(x, y, a, b) and x_geq(x, y, a, b) and y_geq(x, y, a, b)


def big_boss():
    print("enter n:")
    n = int(input())
    a = np.linspace(0, 1, n)
    b = np.linspace(0, 1, n)
    x = np.linspace(0, 1, n)
    y = np.linspace(0, 1, n)

    sum = 0

    for i_a in range(1, n - 1):
        for j_b in range(1, n - 1):
            cur_sum = 0
            for i in range(n):
                for j in range(n):
                    if hit(x[i], y[j], a[i_a], b[j_b]):
                        cur_sum += 1
            sum += cur_sum

    print("sum of area = ", sum, "\n")

    print("total sum = ", n ** 2, "\n")

    print("percent = ", sum / ((n ** 2) * ((n - 2) ** 2)) * 100, "%")



big_boss()