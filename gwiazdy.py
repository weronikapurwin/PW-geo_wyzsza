import numpy as np


# Gwiazda z gwiazdozbioru Byka - Elnath
# Wawa / Nairobi (stolica Kenii) / Puerto Montt (miasto w Chile)

def gregorian_to_julian(y, m, d, h):
    jd = np.floor(365.25 * (y + 4716)) + np.floor(30.6001 * (m + 1)) + d + h / 24 - 1537.5
    if m <= 2:
        y = y - 1
        m = m + 12
        return np.floor(365.25 * (y + 4716)) + np.floor(30.6001 * (m + 1)) + d + h / 24 - 1537.5
    else:
        return jd


def GMST(y, m, d, h):
    T = (gregorian_to_julian(y, m, d, h) - 2451545) / 36525
    g = 280.46061837 + 360.98564736629 * (gregorian_to_julian(y, m, d, h) - 2451545.0) + 0.000387933 * T ** 2 - T ** 3 \
        / 38710000
    g = g % 360
    return g


def kat_godzinny(y, m, d, h, lam, alfa):
    g = GMST(y, m, d, 0)
    UT1 = h * 1.002737909350795
    S = UT1 * 15 + lam + g
    t = S - alfa * 15
    return t

def odleg_zenit(fi, dek, t):
    z = np.sin(fi) * np.sin(dek) + np.cos(fi) * np.cos(dek) * np.cos(t)
    return z

def azymut(fi, dek, t): #tuuu dorobić
    tag_A = (-np.cos(dek) * np.sin(t))/(np.cos(fi) * np.sin(dek) - np.sin(fi) * np.sin(dek) * np.cos(t))
    return tag_A

if __name__ == "__main__":
    print(gregorian_to_julian(1600, 1, 1, 0))

    #t = kat_godzinny(y, m, d, h, lam, alfa)
