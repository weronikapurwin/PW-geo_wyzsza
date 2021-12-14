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

def odleg_zenit(rek, dek, t):
    z = np.sin(rek) * np.sin(dek) + np.cos(rek) * np.cos(dek) * np.cos(t)
    return z

def azymut(rek, dek, t): #tuuu dorobić
    tag_A = (-np.cos(dek) * np.sin(t))/(np.cos(rek) * np.sin(dek) - np.sin(rek) * np.sin(dek) * np.cos(t))
    return tag_A

def godzina_na_dzies(h, m, s):
    return h+ m/60 +s/3600

if __name__ == "__main__":
    # współrzedne miejsc
    Wawa = [52.24920, 21.00030]
    Nairobi = [-1.28333, 36.81667]
    Puerto_Montt = [-41.4693, -72.94237]

    godz = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]

    # Rektascensja przeliczona już z godzin na stopnie
    rek = [81.92083, 81.92083, 81.92083, 81.92083, 81.92083, 81.92087, 81.92087, 81.92087, 81.92087, 81.92087,
           81.92087, 81.92087, 81.92087, 81.92087, 81.92087, 81.92087, 81.92087, 81.92087, 81.92087, 81.92087,
           81.92087, 81.92087, 81.92087, 81.92087, 81.92092]

    dek = [28.62483, 28.62483, 28.62483, 28.62483, 28.62483, 28.62483, 28.62483, 28.62483, 28.62483, 28.62483,
           28.62483, 28.62483, 28.62483, 28.62483, 28.62483, 28.62483, 28.62483, 28.62483, 28.62483, 28.62483,
           28.62483, 28.62483, 28.62483, 28.62483, 28.62483]

    print(15 * godzina_na_dzies(5, 27, 41.02), godzina_na_dzies(28, 37, 29.4))
    print(len(rek))
    print(len(dek))
    t_wawa = kat_godzinny(2021, 11, 30, 0, Wawa[1], rek) #?
    print(t_wawa)
