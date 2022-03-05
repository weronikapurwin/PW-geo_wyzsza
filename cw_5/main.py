import numpy as np

def dzies_na_stop(x):
    h = int(x)
    m = int((x - h) * 60)
    s = round((x - h - m/60) * 3600, 5)
    h = str(h) + "°"
    m =str(m) + "'"
    s = str(s) + "''"
    hms = str(h+m+s)
    return hms

def to_xyz(fi, lam, h):
    fi = np.deg2rad(fi)
    lam = np.deg2rad(lam)
    e2 = 0.0066943800290
    a = 6378137
    n = a / np.sqrt(1 - e2 * np.sin(fi) ** 2)
    x = (n + h) * np.cos(fi) * np.cos(lam)
    y = (n + h) * np.cos(fi) * np.sin(lam)
    z = (n * (1 - e2) + h) * np.sin(fi)
    return x, y, z

def hirvonen(x, y, z, a, e2):
    r = (x ** 2 + y ** 2) ** 0.5
    fi = np.arctan((z/r) * (1-e2) ** -1)

    n = a/np.sqrt(1-e2 * np.sin(fi) ** 2)
    h = r/np.cos(fi) - n
    fi2 = np.arctan((z / r) * (1 - e2 * (n / (n + h))) ** -1)
    sek = np.deg2rad(0.00005/3600)
    while abs(fi2-fi) >= sek:
        fi = fi2
        n = a / np.sqrt(1 - e2 * np.sin(fi) ** 2)
        h = r / np.cos(fi) - n
        fi2 = np.arctan((z / r) * (1 - e2 * (n / (n + h))) ** -1)
    n = a / np.sqrt(1 - e2 * np.sin(fi2) ** 2)
    h = r / np.cos(fi2) - n
    lam = np.arctan(y/x)

    fi2 = np.rad2deg(fi2)
    lam = np.rad2deg(lam)

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
    matrix_p = np.array([x,
                y,
                z])
    matrix_0 = np.array([x0,
                y0,
                z0])
    matrix_w = matrix_p + matrix @ matrix_p + matrix_0

    return matrix_w


if __name__ == "__main__":
    a = 6378137  # metry
    e2 = 0.0066943800290

    a1 = 6378245
    e21 = 0.0066934215520

    A = [50.25, 20.75]
    B = [50, 20.75]
    C = [50.25, 21.25]
    D = [50, 21.25]
    srodkowy = [50.12527054195198, 21.00065090208314]
    sr_szero = [50.125, 21]
    h = 0

    xa, ya, za = to_xyz(A[0], A[1], h)
    xb, yb, zb = to_xyz(B[0], B[1], h)
    xc, yc, zc = to_xyz(C[0], C[1], h)
    xd, yd, zd = to_xyz(D[0], D[1], h)
    xsrodkowy, ysrodkowy, zsrodkowy = to_xyz(srodkowy[0], srodkowy[1], h)
    xsr_szero, ysr_szero, zsr_szero = to_xyz(sr_szero[0], sr_szero[1], h)

    fi_A, lam_A, h_A = trans(xa, ya, za)
    fi_B, lam_B, h_B = trans(xb, yb, zb)
    fi_C, lam_C, h_C = trans(xc, yc, zc)
    fi_D, lam_D, h_D = trans(xd, yd, zd)
    fi_srodkowy, lam_srodkowy, h_srodkowy = trans(xsrodkowy, ysrodkowy, zsrodkowy)
    fi_sr_szero, lam_sr_szero, h_sr_szero = trans(xsr_szero, ysr_szero, zsr_szero)

    print("A:", fi_A, lam_A, h_A, '\n',
          "B:", fi_B, lam_B, h_B, '\n',
          "C:", fi_C, lam_C, h_C, '\n',
          "D:", fi_D, lam_D, h_D, '\n',
          "srodkowy:", fi_srodkowy, lam_srodkowy, h_srodkowy, '\n',
          "średniej szerokości:", fi_sr_szero, lam_sr_szero, h_sr_szero, '\n',)

    fia, lama, ha = hirvonen(fi_A, lam_A, h_A, a1, e21)
    fib, lamb, hb = hirvonen(fi_B, lam_B, h_B, a1, e21)
    fic, lamc, hc = hirvonen(fi_C, lam_C, h_C, a1, e21)
    fid, lamd, hd = hirvonen(fi_D, lam_D, h_D, a1, e21)
    fisrod, lamsrod, hsrod = hirvonen(fi_srodkowy, lam_srodkowy, h_srodkowy, a1, e21)
    fisrszer, lamsrszer, hsrszer = hirvonen(fi_sr_szero, lam_sr_szero, h_sr_szero, a1, e21)

    #zamiana dzies na stopnie
    fia = dzies_na_stop(float(fia))
    lama = dzies_na_stop(float(lama))

    fib = dzies_na_stop(float(fib))
    lamb = dzies_na_stop(float(lamb))

    fic = dzies_na_stop(float(fic))
    lamc = dzies_na_stop(float(lamc))

    fid = dzies_na_stop(float(fid))
    lamd = dzies_na_stop(float(lamd))

    fisrod = dzies_na_stop(float(fisrod))
    lamsrod = dzies_na_stop(float(lamsrod))

    fisrszer = dzies_na_stop(float(fisrszer))
    lamsrszer = dzies_na_stop(float(lamsrszer))

    print("Punkt A:", fia, lama, round(ha, 3), '\n',
          "Punkt B:", fib, lamb, round(hb, 3), '\n',
          "Punkt C:", fic, lamc, round(hc, 3), '\n',
          "Punkt D:", fid, lamd, round(hd, 3), '\n',
          "Punkt środkowy:", fisrod, lamsrod, round(hsrod, 3), '\n',
          "Punkt średniej szerokości:", fisrszer, lamsrszer, round(hsrszer, 3), '\n')
    

