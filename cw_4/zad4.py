import numpy as np
from shapely.geometry import Polygon

def to_GK(fi, lam, lam0):
    fi = np.deg2rad(fi)
    lam = np.deg2rad(lam)
    e2 = 0.0066943800290
    a = 6378137

    b2 = a ** 2 * (1 - e2)
    e2_prim = (a ** 2 - b2) / b2
    lam0 = np.deg2rad(lam0)
    t = np.tan(fi)
    N = a / np.sqrt(1 - e2 * np.sin(fi) ** 2)
    eta2 = e2_prim * (np.cos(fi) ** 2)
    delta_lam = lam - lam0

    A0 = 1 - e2 / 4 - 3 * e2 ** 2 / 64 - 5 * e2 ** 3 / 256
    A2 = 3 / 8 * (e2 + e2 ** 2 / 4 + 15 * e2 ** 3 / 128)
    A4 = 15 / 256 * (e2 ** 2 + 3 * e2 ** 3 / 4)
    A6 = 35 * e2 ** 3 / 3072

    sigma = a * (A0 * fi - A2 * np.sin(2 * fi) + A4 * np.sin(4 * fi) - A6 * np.sin(6 * fi))

    x = sigma + delta_lam ** 2 / 2 * N * np.sin(fi) * np.cos(fi) * (1 + delta_lam ** 2 / 12 * (np.cos(fi) ** 2) * (5 -
        t ** 2 + 9 * eta2 + 4 * eta2 ** 2) + delta_lam ** 4 / 360 * (np.cos(fi) ** 4) * (61 - 58 * t ** 2 + t ** 4 +
        270 * eta2 - 330 * eta2 * t ** 2))

    y = delta_lam * N * np.cos(fi) * (
                1 + delta_lam ** 2 / 6 * np.cos(fi) ** 2 * (1 - t ** 2 + eta2) + delta_lam ** 4 / 120
                * np.cos(fi) ** 4 * (5 - 18 * t ** 2 + t ** 4 + 14 * eta2 - 58 * eta2 * t ** 2))
    return x, y


def from_GK(x, y, lam0):
    e2 = 0.0066943800290
    a = 6378137
    lam0 = np.deg2rad(lam0)
    a2 = a ** 2
    b2 = a2 * (1 - e2)
    e2_prim = (a2 - b2) / b2

    A0 = 1 - e2 / 4 - 3 * e2 ** 2 / 64 - 5 * e2 ** 3 / 256
    A2 = 3 / 8 * (e2 + e2 ** 2 / 4 + 15 * e2 ** 3 / 128)
    A4 = 15 / 256 * (e2 ** 2 + 3 * e2 ** 3 / 4)
    A6 = 35 * e2 ** 3 / 3072

    fi = x / (a * A0)
    sigma = a * (A0 * fi - A2 * np.sin(2 * fi) + A4 * np.sin(4 * fi) - A6 * np.sin(6 * fi))

    while True:
        fi2 = fi + (x - sigma) / (a * A0)

        sigma = a * (A0 * fi2 - A2 * np.sin(2 * fi2) + A4 * np.sin(4 * fi2) - A6 * np.sin(6 * fi2))
        N = a / (1 - e2 * np.sin(fi2) ** 2) ** 0.5
        M = (a * (1 - e2)) / (((1 - e2 * np.sin(fi2)) ** 2) ** 3) ** 0.5
        t = np.tan(fi2)
        eta2 = e2_prim * np.cos(fi2) ** 2

        if abs(fi2 - fi) < np.deg2rad(0.000001 / 3600):
            break

        fi = fi2

    fi = fi2 - y ** 2 * t / (2 * M * N) * (1 - y ** 2 / (12 * N ** 2) * (5 + 3 * t ** 2 + eta2 - 9 * eta2 *
                                                                         t ** 2 - 4 * eta2 ** 2) + y ** 4 / (
                                                       360 * N ** 4) * (61 + 90 * t ** 2 + 45 * t ** 4))

    lam = lam0 + y / (N * np.cos(fi2) * (1 - y ** 2 / (6 * N ** 2) * (1 + 2 * t ** 2 + eta2) + y ** 4 / (120 * N ** 4)
                                         * (5 + 28 * t ** 2 + 24 * t ** 4 + 6 * eta2 + 8 * eta2 * t ** 2)))
    fi = np.rad2deg(fi)
    lam = np.rad2deg(lam)
    return fi, lam


def to_1992(x, y):
    m = 0.9993
    x = m * x - 5300000
    y = m * y + 500000
    return x, y

def from_1992(x, y):
    m = 0.9993
    x = (x + 5300000)/m
    y = (y - 500000)/m
    return x, y


def to_2000(x, y, lam):
    lam = np.deg2rad(lam)
    m = 0.999923
    nr = 1
    if np.deg2rad(13.5) <= lam <= np.deg2rad(16.5):
        nr = 5
    elif np.deg2rad(16.5) <= lam <= np.deg2rad(19.5):
        nr = 6
    elif np.deg2rad(19.5) <= lam <= np.deg2rad(22.5):
        nr = 7
    elif np.deg2rad(22.5) <= lam <= np.deg2rad(25.5):
        nr = 8

    x = m * x
    y = m * y + nr * 1000000 + 500000

    return x, y

def from_2000(x, y, nr):
    m = 0.999923
    x = x/m
    y = (y - 500000 - nr * 1000000)/m
    return x, y

def skala_1992(xgk, ygk):
    a = 6378137
    e2 = 0.0066943800290
    x, y = from_1992(xgk, ygk)
    fi, lam = from_GK(x, y, 19)
    M = a * (1 - e2)/(1 - e2 * np.sin(fi) ** 2) ** (3/2)
    N = a / (1 - e2 * np.sin(fi) ** 2) ** 0.5

    Q = np.sqrt(M * N)

    mgk = 1 + y ** 2/(2 * Q ** 2) + y ** 2/(24 * Q ** 4)
    m92 = 0.9993 * mgk
    kappa = (1 - m92)*1000

    return m92, kappa

def skala_2000(xgk, ygk):
    m0 = 0.999923
    x, y = from_2000(xgk, ygk, 7)
    fi, lam = from_GK(x, y, 21)
    M = a * (1 - e2) / (1 - e2 * np.sin(fi) ** 2) ** (3 / 2)
    N = a / (1 - e2 * np.sin(fi) ** 2) ** 0.5
    Q = np.sqrt(M * N)

    mgk = 1 + y ** 2 / (2 * Q ** 2) + y ** 2 / (24 * Q ** 4)
    m2000 = m0 * mgk
    kappa = (1 - m2000)*1000

    return m2000, kappa

def skala_GK(xgk, ygk):
    fi, lam = from_GK(xgk, ygk, 19)
    M = a * (1 - e2) / (1 - e2 * np.sin(fi) ** 2) ** (3 / 2)
    N = a / (1 - e2 * np.sin(fi) ** 2) ** 0.5
    Q = np.sqrt(M * N)
    mgk = 1 + ygk ** 2 / (2 * Q ** 2) + ygk ** 2 / (24 * Q ** 4)
    kappa = (1 - mgk)*1000
    return mgk, kappa

def skala_pola(m, kappa):
    m = m ** 2
    kappa = (1-m) * 10000
    return m, kappa


if __name__ == "__main__":
    A = [50.25, 20.75]
    B = [50, 20.75]
    C = [50.25, 21.25]
    D = [50, 21.25]
    srodkowy = [50.12527054195198, 21.00065090208314]
    sr_szero = [50.125, 21]

    a = 6378137  # metry
    e2 = 0.0066943800290  # bez jednostek

    x_kg, y_kg = to_GK(A[0], A[1], 21)
    fi, lam = from_GK(x_kg, y_kg, 19)
    x_1992, y_1992 = to_1992(x_kg, y_kg)
    x_2000, y_2000 = to_2000(x_kg, y_kg, A[1])
    print(x_kg, y_kg, '\n', x_1992, y_1992, '\n', fi, lam, '\n', x_2000, y_2000)

    x1_kg, y1_kg = to_GK(D[0], D[1], 19)
    x1_1992, y1_1992 = to_1992(x1_kg, y1_kg)
    x1_2000, y1_2000 = to_2000(x1_kg, y1_kg, D[1])
    print('\n', x1_kg, y1_kg, '\n', x1_1992, y1_1992, '\n', x1_2000, y1_2000)

    xgk = [5570120.59683, 5542315.02536,5543273.89185, 5571077.96006, 5570120.59683]
    ygk = [124812.22774, 125464.20084, 161308.28340, 160469.90666, 124812.22774]

    x1992 = [266221.51242, 238435.40484, 239393.60013, 267178.20549, 266221.51242]
    y1992 = [624724.85918, 625376.37591, 661195.36761, 660357.57773, 624724.85918]

    x2000 = [5568256.02985, 5540450.34986, 5540450.34986, 5568256.02985, 5568256.02985]
    y2000 = [7482170.56245, 7482077.45154,  7517922.54845,7517829.43754, 7482170.56245]

    pgongk = Polygon(zip(xgk, ygk))
    pgon1992 = Polygon(zip(x1992, y1992))
    pgon2000 = Polygon(zip(x2000, y2000))
    print("GK", pgongk.area, '\n', "1992", pgon1992.area, '\n', "2000", pgon2000.area)

    m1992_A, kappa1992_A = skala_1992(266221.51242, 624724.85918)
    m1992_B, kappa1992_B = skala_1992(238435.40484, 625376.37591)
    m1992_C, kappa1992_C = skala_1992(267178.20549, 660357.57773)
    m1992_D, kappa1992_D = skala_1992(239393.60013, 661195.36761)
    m1992_srodkowy, kappa1992_srodkowy = skala_1992(252808.42569, 642959.83145)
    m1992_srszer, kappa1992_srszer = skala_1992(252777.11061, 642914.12934)

    m2000_A, kappa2000_A = skala_2000(5568256.02985, 7482170.56245)
    m2000_B, kappa2000_B = skala_2000(5540450.34986, 7482077.45154)
    m2000_C, kappa2000_C = skala_2000(5568256.02985, 7517829.43754)
    m2000_D, kappa2000_D = skala_2000(5540450.34986, 7517922.54845)
    m2000_srszer, kappa2000_srszer = skala_2000(5554323.10967,7500000.0)
    m2000_srodkowy, kappa2000_srodkowy = skala_2000(5554353.20033,7500046.54196)

    mGK_A, kappaGK_A = skala_GK(5570120.59683, 124812.22774)
    mGK_B, kappaGK_B = skala_GK(5542315.02536,125464.20084)
    mGK_C, kappaGK_C = skala_GK(5571077.96006, 160469.90666)
    mGK_D, kappaGK_D = skala_GK(5543273.89185, 161308.28340)
    mGK_srszer, kappaGK_srszer = skala_GK(5556666.77736, 143014.23931)
    mGK_srodkowy, kappaGK_srodkowy = skala_GK(5556698.11437,143059.97343)

    GK_m_A, GK_kappa_A = skala_pola(mGK_A, kappaGK_A)
    GK_m_B, GK_kappa_B= skala_pola(mGK_B, kappaGK_B)
    GK_m_C, GK_kappa_C = skala_pola(mGK_C, kappaGK_C)
    GK_m_D, GK_kappa_D= skala_pola(mGK_D, kappaGK_D)
    GK_m_srszer, GK_kappa_srszer = skala_pola(mGK_srszer, kappaGK_srszer)
    GK_m_srodkowy, GK_kappa_srodkowy = skala_pola(mGK_srodkowy, kappaGK_srodkowy)

    m_A_2000, kappa_A_2000 = skala_pola(m2000_A, kappa2000_A)
    m_B_2000, kappa_B_2000 = skala_pola(m2000_B, kappa2000_B)
    m_C_2000, kappa_C_2000 = skala_pola(m2000_C, kappa2000_C)
    m_D_2000, kappa_D_2000 = skala_pola(m2000_D, kappa2000_D)
    m_srszer_2000, kappa_srszer_2000 = skala_pola(m2000_srszer, kappa2000_srszer)
    m_srodkowy_2000, kappa_srodkowy_2000 = skala_pola(m2000_srodkowy, kappa2000_srodkowy)

    m_A_1992, kappa_A_1992 = skala_pola(m1992_A, kappa1992_A)
    m_B_1992, kappa_B_1992 = skala_pola(m1992_B, kappa1992_B)
    m_C_1992, kappa_C_1992 = skala_pola(m1992_C, kappa1992_C)
    m_D_1992, kappa_D_1992 = skala_pola(m1992_D, kappa1992_D)
    m_srszer_1992, kappa_srszer_1992 = skala_pola(m1992_srszer, kappa1992_srszer)
    m_srodkowy_1992, kappa_srodkowy_1992 = skala_pola(m1992_srodkowy, kappa1992_srodkowy)
