import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm


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
    z = np.rad2deg(np.arccos(z))
    return z


def azymut(fi, dek, t):
    fi = np.deg2rad(fi)
    dek = np.deg2rad(dek)
    t = np.deg2rad(t)
    licz = -np.cos(dek) * np.sin(t)
    mian = np.cos(fi) * np.sin(dek) - np.sin(fi) * np.cos(dek) * np.cos(t)
    tag_A = np.arctan(licz / mian)
    tag_A = np.rad2deg(tag_A)
    if mian < 0:
        tag_A += 180
    elif mian > 0 and licz < 0:
        tag_A += 360
    return tag_A

def wysokosc(fi, dek, t):
    h = 90 - odleg_zenit(fi, dek, t)
    return h

def transf_wspol(fi, dek, t):
    x = np.sin(np.deg2rad(odleg_zenit(fi, dek, t))) * np.cos(np.deg2rad(azymut(fi, dek, t)))
    y = np.sin(np.deg2rad(odleg_zenit(fi, dek, t))) * np.sin(np.deg2rad(azymut(fi, dek, t)))
    z = np.cos(np.deg2rad(odleg_zenit(fi, dek, t)))
    return x, y, z

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

    t_wawa = kat_godzinny(2021, 11, 30, godz, Wawa[1], rek)
    t_nairobi = kat_godzinny(2021, 11, 30, godz, Nairobi[1], rek)
    t_puerto_montt = kat_godzinny(2021, 11, 30, godz, Puerto_Montt[1], rek)

    i = 0
    while i < 25:
        # print("kąt godzinny dla Warszawy:", t_wawa[i],"godzina:", godz[i],)
        i = i + 1

    b = np.vectorize(azymut)
    print("azymut dla wawy:", b(Wawa[0], dek, t_wawa))
    c = np.vectorize(wysokosc)
    print("wysokosc dla wawy:", c(Wawa[0], dek, t_wawa))
    d = np.vectorize(transf_wspol)
    x_wa, y_wa, z_wa = d(Wawa[0], dek, t_wawa)
    x_pm, y_pm, z_pm = d(Puerto_Montt[0], dek, t_puerto_montt)
    x_n, y_n, z_n = d(Nairobi[0], dek, t_nairobi)
    print(len(x_wa))

    u, v = np.mgrid[0:2 * np.pi:30j, 0:np.pi:20j]
    x = np.cos(u) * np.sin(v)
    y = np.sin(u) * np.sin(v)
    z = np.cos(v)

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.set_title("Ruch gwiazdy na niebie")
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.plot_surface(x, y, z, alpha=0.05, facecolors=cm.jet(z/np.amax(z)))
    ax.plot3D(x_wa, y_wa, z_wa, 'grey')
    ax.scatter3D(x_wa, y_wa, z_wa, c='black', cmap='cividis')

    plt.show()
