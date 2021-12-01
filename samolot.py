import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt


# przeliczanie wspolrzednych fi, lam, h na wspolrzedne x, y, z
def geo_to_xyz(fi, lam, h, a, e2):
    fi = np.deg2rad(fi)
    lam = np.deg2rad(lam)
    N = a / (1 - (e2) * np.sin(fi) ** 2) ** 0.5
    x = (N + h) * np.cos(fi) * np.cos(lam)
    y = (N + h) * np.cos(fi) * np.sin(lam)
    z = ((N * (1 - e2) + h) * np.sin(fi))
    return np.array((x, y, z))


# przeliczanie wspolrzednych fi, lam, h na wspolrzedne n, e, u
def geo_to_neu(F1, L1, H1, F2, L2, H2):
    pkt_1 = geo_to_xyz(F1, L1, H1, a, e2)
    pkt_2 = geo_to_xyz(F2, L2, H2, a, e2)

    R = np.array(
        [[-np.sin(F1) * np.cos(L1), -np.sin(L1), np.cos(F1) * np.cos(L1)],
         [-np.sin(F1) * np.sin(L1), np.cos(L1), np.cos(F1) * np.sin(L1)],
         [np.cos(F1), 0, np.sin(F1)]
         ])
    R = R.transpose()

    D = np.array([[pkt_2[0] - pkt_1[0]],
                  [pkt_2[1] - pkt_1[1]],
                  [pkt_2[2] - pkt_1[2]]]).squeeze()
    neu = R @ D
    return neu


# funckja liczaca odleglosc skosna (w metrach)
def odleg_skos(N, E, U):
    return (N ** 2 + E ** 2 + U ** 2) ** 0.5


# funckja liczaca azymut (w stopniach)
def azymut(E, N):
    if N == 0:
        return 0
    else:
        az = np.rad2deg(np.arctan(E / N))
        if (N < 0 and E > 0) or (N < 0 and E < 0):
            return az + 180
        elif N > 0 and E < 0:
            return az + 360
        else:
            return az


def odleg_zenit(N, E, U):
    if odleg_skos(N, E, U) != 0:
        return U / odleg_skos(N, E, U)
    else:
        pass


if __name__ == "__main__":
    a = 6378137  # metry
    e2 = 0.0066943800290  # bez jednostek

    # wspolrzedne lotniska wylotu
    F = 48.9959
    L = 2.5517
    H = 404.00

    # wczytanie danych lotu
    data = np.genfromtxt('dane.txt', delimiter=';', dtype="U75")
    dane_fi = data[:, 0:1].astype(float).squeeze()
    dane_lam = data[:, 1:2].astype(float).squeeze()
    dane_h = data[:, 2:3].astype(float).squeeze()

    # obliczenia
    n, e, u = geo_to_neu(F, L, H, dane_fi, dane_lam, dane_h)
    print("odległość skośna:", odleg_skos(n, e, u))
    b = np.vectorize(azymut)
    print("azymut:", b(e, n))
    c = np.vectorize(odleg_zenit)
    print("odległość zenitalna:", c(n, e, u))

    # wyswietlenie trasy lotu 2d
    airport_data = pd.read_csv('dane2.txt', delimiter=';')
    airport_gdf = gpd.GeoDataFrame(airport_data, geometry=gpd.points_from_xy(airport_data['lam'],
                                                                             airport_data['fi']))
    airport_gdf.plot(markersize=1.5, figsize=(10, 10))
    plt.show()

    # wyswietlenie trasy lotu 3d
    fig = plt.figure()
    ax = plt.axes(projection="3d")
    ax.plot3D(*geo_to_neu(F, L, H, dane_fi, dane_lam, dane_h), 'red')
    ax.scatter3D(*geo_to_neu(F, L, H, dane_fi, dane_lam, dane_h), cmap='cividis')
    
    # punkt, gdy samolot znika za horyzontem
    for a, b, c in zip(n, e, u):
        if c > 0:
            ax.plot3D(a, b, c, "red")
            ax.text(a, b, c, '%s' % ("samolot znika za horyzontem"), size=10, zorder=1, color='k')
            print("samolot zniknie za horyzontem we współrzędnych NEU: ", round(a, 3), round(b, 3), round(c, 3))
            break
    ax.set_xlabel('$N$')
    ax.set_ylabel('$E$')
    plt.show()
