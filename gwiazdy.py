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
    for value in alfa:
        value *= 15
    t = S - alfa
    return t

def odleg_zenit(fi, dek, t):
    fi = np.deg2rad(fi)
    dek = np.deg2rad(dek)
    t = np.deg2rad(t)
    z = np.sin(fi) * np.sin(dek) + np.cos(fi) * np.cos(dek) * np.cos(t)
    return z

def azymut(fi, dek, t): #tuuu dorobić
    fi = np.deg2rad(fi)
    dek = np.deg2rad(dek)
    t = np.deg2rad(t)
    tag_A = (-np.cos(dek) * np.sin(t))/(np.cos(fi) * np.sin(dek) - np.sin(fi) * np.cos(dek) * np.cos(t))
    return tag_A

def godzina_na_dzies(h, m, s):
    return h + m/60 + s / 3600

if __name__ == "__main__":
    # współrzedne miejsc
    Wawa = [52.24920, 21.00030]
    Nairobi = [-1.28333, 36.81667]
    Puerto_Montt = [-41.4693, -72.94237]

    godz = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
    godz = np.asarray(godz)
    # Rektascensja przeliczona już z godzin na stopnie
    rek = [81.92083, 81.92083, 81.92083, 81.92083, 81.92083, 81.92087, 81.92087, 81.92087, 81.92087, 81.92087,
           81.92087, 81.92087, 81.92087, 81.92087, 81.92087, 81.92087, 81.92087, 81.92087, 81.92087, 81.92087,
           81.92087, 81.92087, 81.92087, 81.92087, 81.92092]

    dek = [28.62483, 28.62483, 28.62483, 28.62483, 28.62483, 28.62483, 28.62483, 28.62483, 28.62483, 28.62483,
           28.62483, 28.62483, 28.62483, 28.62483, 28.62483, 28.62483, 28.62483, 28.62483, 28.62483, 28.62483,
           28.62483, 28.62483, 28.62483, 28.62483, 28.62483]

    print(15 * godzina_na_dzies(0, 32, 0))
    t_wawa = kat_godzinny(2021, 11, 30, godz, Wawa[1], rek)
    t_nairobi = kat_godzinny(2021, 11, 30, godz, Nairobi[1], rek)
    t_puerto_montt = kat_godzinny(2021, 11, 30, godz, Puerto_Montt[1], rek)

    i = 0
    while i < 25:
        print("kąt godzinny dla Warszawy:", t_wawa[i],"godzina:", godz[i],)
        i = i+1
