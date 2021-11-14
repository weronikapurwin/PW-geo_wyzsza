import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt

def geo_to_xyz(fi, lam, h, a, e2):
    fi = np.deg2rad(fi)
    lam = np.deg2rad(lam)
    N = a / (1 - (e2) * np.sin(fi) ** 2) ** 0.5
    x = (N + h) * np.cos(fi) * np.cos(lam)
    y = (N + h) * np.cos(fi) * np.sin(lam)
    z = ((N * (1 - e2) + h) * np.sin(fi))
    return np.array((x, y, z))


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

# def odleg_skos(N,E,U):

if __name__ == "__main__":
    a = 6378137  # metry
    e2 = 0.0066943800290  # bez jednostek

    # wspolrzedne lotniska wylotu
    F = 48.9959
    L = 2.5517
    H = 404.00

    # wspolrzedne lotniska w neu
    data = np.genfromtxt('dane.txt',delimiter=';', dtype="U75")
    dane_fi = data[:, 0:1].astype(float).squeeze()
    dane_lam = data[:, 1:2].astype(float).squeeze()
    dane_h = data[:, 2:3].astype(float).squeeze()
    geo_to_neu(F, L, H, dane_fi, dane_lam, dane_h)
    
    airport_data = pd.read_csv('dane2.txt', delimiter=';')
    airport_gdf = gpd.GeoDataFrame(airport_data, geometry=gpd.points_from_xy(airport_data['lam'],
                                                                             airport_data['fi']))
    airport_gdf.plot(markersize=1.5, figsize=(10, 10))
    plt.show()

    fig = plt.figure()
    ax = plt.axes(projection="3d")
    ax.plot3D(*geo_to_neu(F, L, H, dane_fi, dane_lam, dane_h), 'red')
    ax.scatter3D(*geo_to_neu(F, L, H, dane_fi, dane_lam, dane_h), cmap='cividis')
    ax.set_xlabel('$N$')
    ax.set_ylabel('$E$')
    plt.show()
