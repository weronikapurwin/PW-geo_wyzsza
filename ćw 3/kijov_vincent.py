import numpy as np

def punkt_sr_szer(a, b):
    return (a + b) / 2

def dzies_na_stop(x):
    h = np.floor(x)
    m = np.floor((x - h) * 60)
    s = round((x - h - m/60) * 3600, 5)
    return h, m, s


def vincent(a, e2, fi1, lam1, fi2, lam2):
    b = a * (1 - e2) ** 0.5
    f = 1 - b / a

    fi1 = np.deg2rad(fi1)
    lam1 = np.deg2rad(lam1)
    fi2 = np.deg2rad(fi2)
    lam2 = np.deg2rad(lam2)

    delta_lam = lam2 - lam1

    U1 = np.arctan((1 - f) * np.tan(fi1))
    U2 = np.arctan((1 - f) * np.tan(fi2))

    L = delta_lam
    sekunda = np.deg2rad(0.000001 / 3600)

    while True:
        sin_sigma = ((np.cos(U2) * np.sin(L)) ** 2 +
                      (np.cos(U1) * np.sin(U2) - np.sin(U1) * np.cos(U2) * np.cos(L)) ** 2) ** 0.5

        cos_sigma = np.sin(U1) * np.sin(U2) + np.cos(U1) * np.cos(U2) * np.cos(L)

        sigma = np.arctan(sin_sigma / cos_sigma)

        sin_alfa = np.cos(U1) * np.cos(U2) * np.sin(L) / np.sin(sigma)

        cos_alfa_kwadrat = 1 - sin_alfa ** 2

        cos_2sigma = cos_alfa_kwadrat ** 0.5 - (2 * np.sin(U1) * np.sin(U2) / cos_alfa_kwadrat)

        C = f / 16 * cos_alfa_kwadrat * (4 + f * (4 - 3 * cos_alfa_kwadrat))

        L2 = delta_lam + (1 - C) * f * sin_alfa * (sigma + C * sin_sigma * (cos_2sigma + C * cos_sigma *
                                                                    (-1 + 2 * cos_2sigma ** 2)))
        if abs(L2 - L) < sekunda:
            break

        L = L2

    u2 = (a ** 2 - b ** 2) / b ** 2 * cos_alfa_kwadrat
    A = 1 + u2 / 16384 * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
    B = u2 / 1024 * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))
    delta_sigma = B * sin_sigma * (
            cos_2sigma + 1 / 4 * B * (cos_sigma * (-1 + 2 * cos_2sigma ** 2) - 1 / 6 * B * cos_2sigma
                                      * (-3 + 4 * sin_sigma ** 2) * (-3 + 4 * cos_2sigma ** 2)))
    s_AB = b * A * (sigma - delta_sigma)

    num_ab = np.cos(U2) * np.sin(L2)
    denom_ab = np.cos(U1) * np.sin(U2) - np.sin(U1) * np.cos(U2) * np.cos(L2)
    az_AB = np.rad2deg(np.arctan(num_ab / denom_ab))
    if (denom_ab < 0):
        az_AB += 180
    elif denom_ab > 0 and num_ab < 0:
        az_AB += 360
        return az_AB

    num_ba = np.cos(U2) * np.sin(L2)
    denom_ba = -np.sin(U1) * np.cos(U2) + np.cos(U1) * np.sin(U2) * np.cos(L2)
    az_BA = np.rad2deg(np.arctan(num_ba / denom_ba) + np.pi)
    if (denom_ba < 0):
        az_BA += 180
    elif denom_ba > 0 and num_ba < 0:
        az_BA += 360
        return az_BA

    return s_AB, az_AB, az_BA


def kijovi(a, e2, fi, lam, AzAB, s):
    n = int(s / 1000)
    ds = 1000
    fi = np.deg2rad(fi)
    lam = np.deg2rad(lam)
    AzAB = np.deg2rad(AzAB)

    for i in range(n):
        M = (a * (1 - e2)) / ((1 - e2 * np.sin(fi)) ** 3) ** 0.5
        N = a / (1 - (e2 * np.sin(fi) ** 2)) ** 0.5

        dfi = ds * np.cos(AzAB) / M
        dAzAB = np.sin(AzAB) * np.tan(fi) * ds / N

        fi_m = fi + 0.5 * dfi
        Az_m = AzAB + 0.5 * dAzAB

        Mm = (a * (1 - e2)) / ((1 - e2 * np.sin(fi_m)) ** 3) ** 0.5
        Nm = a / (1 - (e2) * np.sin(fi_m) ** 2) ** 0.5

        dfim = np.cos(Az_m) * ds / Mm
        dlam = np.sin(Az_m) * ds / (Nm * np.cos(fi_m))
        dAz_m = np.sin(Az_m) * np.tan(fi_m) * ds / Nm

        fi = fi + dfim
        lam = lam + dlam
        AzAB =AzAB + dAz_m

    ds = s % 1000
    M = (a * (1 - e2)) / ((1 - e2 * np.sin(fi)) ** 3) ** 0.5
    N = a / (1 - (e2 * np.sin(fi) ** 2)) ** 0.5

    dfi = ds * np.cos(AzAB) / M
    dAzAB = np.sin(AzAB) * np.tan(fi) * ds / N

    fi_m = fi + 0.5 * dfi
    Az_m = AzAB + 0.5 * dAzAB

    Mm = (a * (1 - e2)) / ((1 - e2 * np.sin(fi_m)) ** 3) ** 0.5
    Nm = a / (1 - (e2) * np.sin(fi_m) ** 2) ** 0.5

    dfim = np.cos(Az_m) * ds / Mm
    dlam = np.sin(Az_m) * ds / (Nm * np.cos(fi_m))
    dAz_m = np.sin(Az_m) * np.tan(fi_m) * ds / Nm

    fi = fi + dfim
    lam = lam + dlam
    AzAB = AzAB + dAz_m

    fi = np.rad2deg(fi)
    lam = np.rad2deg(lam)
    #AzBA = np.rad2deg(AzAB)
    AzAB = np.rad2deg(AzAB)
    return fi, lam, AzAB

def pole(fi1, lam1, fi2, lam2):
    e2 = 0.0066943800290
    a = 6378137
    e = e2 ** 0.5
    b = a * (1 - e2) ** 0.5

    fi1 = np.deg2rad(fi1)
    fi2 = np.deg2rad(fi2)
    lam1 = np.deg2rad(lam1)
    lam2 = np.deg2rad(lam2)

    integral1 = np.sin(fi1)/(1 - e2 * np.sin(fi1) ** 2) + 1/(2 * e) * np.log((1 + e * np.sin(fi1))/
                                                                              (1 - e * np.sin(fi1)))
    integral2 = np.sin(fi2)/(1 - e2 * np.sin(fi2) ** 2) + 1/(2 * e) * np.log((1 + e * np.sin(fi2))/
                                                                              (1 - e * np.sin(fi2)))
    integral = integral1 - integral2
    p = b ** 2 * (lam2 - lam1) * 0.5 * integral
    return p


if __name__ == "__main__":
    A = [50.25, 20.75]
    B = [50, 20.75]
    C = [50.25, 21.25]
    D = [50, 21.25]

    a = 6378137  # metry
    e2 = 0.0066943800290  # bez jednostek

    print("współrzędne punktu średniej szerokości:", "fi:", punkt_sr_szer(B[0], C[0]), "lam:", punkt_sr_szer(B[1], C[1]))

    print("algorytm VINCENTEGO")
    s, a_ab, a_ba = vincent(a,e2,A[0], A[1], D[0], D[1])
    print("długość linii geodezyjnej:", round(s, 5), '\n', "Azymut AD:", round(a_ab, 5), '\n', "Azymut DA:", round(a_ba,5))

    print("algorytm KIJOVI")
    fi, lam, az_12 = kijovi(a, e2, A[0], A[1], a_ab, s)
    print("współrzędne punktu środkowego:", "fi:", round(fi, 5), "lam:", round(lam, 5),
          '\n', "Azymut AD:", az_12)

    print("algorytm VINCENTEGO znowu")
    x = punkt_sr_szer(B[0], C[0])
    y = punkt_sr_szer(B[1], C[1])
    s1, a_xy, a_yx = vincent(a,e2, x, y, fi, lam)
    print("Różnica odległości:", round(s1, 5),"Azymut wprost:", a_xy,'\n' "Azymut odwrotny:", a_yx)

    p = pole(A[0], A[1], D[0], D[1])
    print("pole:", p)



