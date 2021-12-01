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
    g = 280.46061837 + 360.98564736629 * (gregorian_to_julian(y, m, d, h) - 2451545.0) + 0.000387933 * T ** 3 / 38710000
    g = g % 360
    return g

def kat_godzinny(y, m, d, h,lam):
    jd = gregorian_to_julian(y, m, d, 0)
    g = 1 #tu zaimplementuj
    UT1 = h * 1.002737909350795
    S = UT1 * 15 + lam + g


if __name__ == "__main__":
    print(gregorian_to_julian(1600, 1, 1, 0))
