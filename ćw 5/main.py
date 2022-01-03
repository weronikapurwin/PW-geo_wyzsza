import numpy as np

def hirvonen(x, y, z, a, e2):
    r = (x ** 2 + y ** 2) ** 0.5
    fi = np.arctan((z/r) * (1-e2) ** -1)

    n = a/(1-e2 * np.sin(fi) ** 2)
    h = r/np.cos(fi) - n
    fi2 = np.arctan((z / r) * (1 - e2 * (n / (n + h))) ** -1)
    sek = np.deg2rad(0.00005/3600)
    while abs(fi2-fi) >= sek:
        fi = fi2
        n = a / (1 - e2 * np.sin(fi) ** 2)
        h = r / np.cos(fi) - n
        fi2 = np.arctan((z / r) * (1 - e2 * (n / (n + h))) ** -1)
    n = a / (1 - e2 * np.sin(fi2) ** 2)
    h = r / np.cos(fi2) - n
    lam = np.arctan(y/x)

    # x = (n+h)*np.cos(fi2)*np.cos(lam)
    # y = (n+h)*np.cos(fi2)*np.sin(lam)
    # z = (n*(1-e2)+h)*np.sin(fi2)
    # print(x,y,z)
    return fi2, lam, h

def trans(x, y, z):
    kappa = 0.8407728 * 10 ** -6
    alfa = np.deg2rad(-0.35867/3600)
    beta = np.deg2rad(-0.05283/3600)
    gamma = np.deg2rad(0.84354/3600)
    x0 = -33.4297
    y0 = 146.5746
    z0 = 76.2865
    matrix = np.array([[kappa, gamma, -beta],
              [-gamma, kappa, alfa],
              [beta, -alfa, kappa]])
    matrix_p = np.array([[x],
                [y],
                [z]])
    matrix_0 = np.array([[x0],
                [y0],
                [z0]])
    matrix_w = matrix_p + matrix @ matrix_p +matrix_0
    print(matrix_w)
    return matrix_w


if __name__ == "__main__":
    a = 6378137  # metry
    e2 = 0.0066943800290
    fia = np.deg2rad(50.25)
    lama = np.deg2rad(20.75)
    fib = np.deg2rad(50)
    lamb = np.deg2rad(20.75)

    hirvonen(1,2,3,a,e2)
    trans(1,2,3)