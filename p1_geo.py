import sys
from math import sin, cos, sqrt, atan, atan2, degrees, radians, tan
from numpy import array

o = object()


class Transformacje:
    def __init__(self, model: str = "wgs84"):
        """
        Parametry elipsoid:
            a - duża półoś elipsoidy - promień równikowy
            b - mała półoś elipsoidy - promień południkowy
            flat - spłaszczenie
            ecc2 - mimośród^2
        + WGS84: https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        + Inne powierzchnie odniesienia: https://en.wikibooks.org/wiki/PROJ.4#Spheroid
        + Parametry planet: https://nssdc.gsfc.nasa.gov/planetary/factsheet/index.html
        """
        if model == "wgs84":
            self.a = 6378137.0  # semimajor_axis
            self.b = 6356752.31424518  # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "mars":
            self.a = 3396900.0
            self.b = 3376097.80585952
        # elif model == "krasowski":
        #     self.a = 6378245.0
        #     self.b = 6356863.01877307
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flat = (self.a - self.b) / self.a
        # eccentricity  WGS84:0.0818191910428
        self.ecc = sqrt(2 * self.flat - self.flat ** 2)
        self.ecc2 = (2 * self.flat - self.flat ** 2)  # eccentricity**2

    def xyz2plh(self, X, Y, Z, output='dec_degree'):
        """
        Algorytm Hirvonena - algorytm transformacji współrzędnych ortokartezjańskich (x, y, z)
        na współrzędne geodezyjne szerokość, długość i wysokośc elipsoidalna (phi, lam, h). 
        W wyniku 3-4-krotnej iteracji wyznaczenia wsp. phi można przeliczyć współrzędne z dokładnoscią ok 1 cm.     
        Parameters
        ----------
        X, Y, Z : FLOAT
             współrzędne w układzie orto-kartezjańskim, 

        Returns
        -------
        lat
            [stopnie dziesiętne] - szerokość geodezyjna
        lon
            [stopnie dziesiętne] - długość geodezyjna.
        h : 
            [metry] - wysokość elipsoidalna
        output [STR] - optional, defoult 
            dec_degree - decimal degree
            dms - degree, minutes, sec
        """
        r = sqrt(X**2 + Y**2)           # promień
        # pierwsze przybliilizenie
        lat_prev = atan(Z / (r * (1 - self.ecc2)))
        lat = 0
        while abs(lat_prev - lat) > 0.000001/206265:
            lat_prev = lat
            N = self.a / sqrt(1 - self.ecc2 * sin(lat_prev)**2)
            h = r / cos(lat_prev) - N
            lat = atan((Z/r) * (((1 - self.ecc2 * N/(N + h))**(-1))))
        lon = atan(Y/X)
        N = self.a / sqrt(1 - self.ecc2 * (sin(lat))**2)
        h = r / cos(lat) - N
        if output == "dec_degree":
            return degrees(lat), degrees(lon), h
        elif output == "dms":
            lat = self.deg2dms(degrees(lat))
            lon = self.deg2dms(degrees(lon))
            return f"{lat[0]:02d}:{lat[1]:02d}:{lat[2]:.2f}", f"{lon[0]:02d}:{lon[1]:02d}:{lon[2]:.2f}", f"{h:.3f}"
        else:
            raise NotImplementedError(f"{output} - output format not defined")

    def plh2xyz(self, phi, lam, h):
        """
        Zadanie odwrotne do Algorytmu Hirvonena; 
        transformacja współrzędnych geodezyjnych i wysokosci elipsoidalnej (phi, lam, h),
        na współrzędne ortokartezjańskie (x, y, z).
        
        Rn - długosc odcinka normalnej, mierzona od punktu P0 do punktu przecięcia z osią obrotu elipsoidy;
        promień krzywizny przekroju poprzecznego (pierwszego wertykału) elipsoidy w punkcie P0.
        
        q - pionowe przesunięcie srodka krzywizny przekroju poprzecznego wzgledem srodka elipsoidy

        Parameters
        ----------
        Phi
            [stopnie dziesiętne] - szerokość geodezyjna
        Lam
            [stopnie dziesiętne] - długośc geodezyjna.
        h : 
            [metry] - wysokość elipsoidalna

        Returns
        -------
        X, Y, Z : FLOAT
            Współrzędne w układzie orto-kartezjańskim

        """
        phi = radians(phi)
        lam = radians(lam)
        Rn = self.a / sqrt(1-self.ecc2 * sin(phi)**2)
        q = Rn * self.ecc2 * sin(phi)

        X = (Rn + h) * cos(phi) * cos(lam)
        Y = (Rn + h) * cos(phi) * sin(lam)
        Z = (Rn + h) * sin(phi) - q
        return X, Y, Z

    def xyz2neu(self, x, y, z, x_0, y_0, z_0):
        """
        Transformacja współrzędnych do układu topocentrycznego
        Współrzędne w układzie topocentrycznym otrzymujemy 
        przez przesunięcie początku układu współrzędnych do wybranego punktu (x_0, y_0, z_0) - translacja, a następnie rotację.
        Paramtery rotacji zależne są od szerokosci i długosci geodezyjnej punktu (phi, lam)
        
        R - macierz obrotu (rotacji)

        Parameters
        ----------
        x, y, z : FLOAT
            Współrzędne geocentryczne punktów, których pozycje względem wybranego punktu obliczamy (koniec wektora )
        x_0, y_0, z_0 : FLOAT
            Współrzędne geocentryczne wybranego punktu (początek wektora)

        Returns
        -------
        N, E, U : FLOAT
            Topocentryczne współrzędne satelitów

        """
        phi, lam, _ = [radians(coord) for coord in self.xyz2plh(x, y, z)]

        R = array([[-sin(lam), -sin(phi)*cos(lam), cos(phi)*cos(lam)],
                   [cos(lam),  -sin(phi)*sin(lam), cos(phi)*sin(lam)],
                   [0,           cos(phi),       sin(phi)]])

        xyz_t = array([[x - x_0],
                       [y - y_0],
                       [z - z_0]])

        [[E], [N], [U]] = R.T @  xyz_t

        return N, E, U

    def sigma_p(self, phi):

        A0 = 1 - (self.ecc2/4) - ((3*(self.ecc2**2))/64) - \
            ((5*(self.ecc2**3))/256)
        A2 = (3/8) * (self.ecc2 + (self.ecc2**2)/4 + (15*(self.ecc2**3))/128)
        A4 = (15/256) * (self.ecc2**2 + (3*(self.ecc2**3))/4)
        A6 = (35*(self.ecc2**3))/3072

        sigma = self.a * (A0 * phi - A2*sin(2*phi) +
                          A4*sin(4*phi) - A6*sin(6*phi))

        return sigma

    def plh2gk(self, phi, lam, lam0):
        """
            Transformacja współrzędnych geodezyjnych (phi, lam), na współrzędne w odwzorowaniu Gaussa-Krugera (x_gk, y_gk).
    
        Parameters
        ----------
        Phi: FLOAT
            Szerokość geodezyjna
        Lam: FLOAT
            Długość geodezyjna.
        Lam0: FLOAT
            Długosć geodezyjna południka osiowego

        Returns
        -------
        x_gk, y_gk : FLOAT
            Współrzędne w odwzorowaniu Gaussa-Krugera

        """
        a2 = self.a**2
        b2 = (a2 * (1 - self.ecc2))
        e_prim_2 = (a2 - b2)/b2
        dl = lam - lam0
        t = tan(phi)
        eta2 = e_prim_2 * (cos(phi)**2)
        N = self.a / sqrt(1 - self.ecc2 * sin(phi)**2)

        sigma = self.sigma_p(phi)

        x_gk = sigma + (dl ** 2 / 2) * N * sin(phi) * cos(phi) \
            * (1 + (dl ** 2 / 12) * cos(phi) ** 2
               * (5 - t ** 2 + 9 * eta2 + 4 * eta2 ** 2)
               + (dl ** 4 / 360) * cos(phi) ** 4
               * (61 - 58 * t ** 2 + t ** 4 + 270 * eta2 - 330 * eta2 * t ** 2))

        y_gk = dl * N * cos(phi) * (1 + (dl ** 2 / 6) * cos(phi) ** 2
                                    * (1 - t ** 2 + eta2)
                                    + (dl ** 4 / 120) * cos(phi) ** 4
                                    * (5 - 18 * t ** 2 + t ** 4 + 14 * eta2 - 58 * eta2 * t ** 2))
        return x_gk, y_gk

    def gk1992(self, x_gk, y_gk):
        """
            Transformacja współrzędnych w odwzorowaniu Gaussa-Krugera 
            do współrzędnych w układzie PL-1992.

        Parameters
        ----------
        x_gk, y_gk : FLOAT
            Współrzędne w odwzorowaniu Gaussa-Krugera

        Returns
        -------
        x_1992, y_1992 : FLOAT
            Współrzędne w układzie PL-1992

        """
        m0 = 0.9993
        x_1992 = x_gk * m0 - 5300000
        y_1992 = y_gk * m0 + 500000
        return x_1992, y_1992

    def plh1992(self, phi, lam):
        """
            Transformacja współrzędnych geodezyjnych (phi, lam), na współrzędne w układzie PL-1992 (x_1992, y_1992).
            
        Parameters
        ----------
        Phi: FLOAT
            Szerokość geodezyjna
        Lam: FLOAT
            Długość geodezyjna.

        Returns
        -------
        x_1992, y_1992 : FLOAT
            Współrzędne w układzie PL-1992

        """

        l0 = radians(19)
        x_gk, y_gk = self.plh2gk(phi, lam, l0)
        x_1992, y_1992 = self.gk1992(x_gk, y_gk)

        return x_1992, y_1992

    def gk2000(self,  x_gk, y_gk, zone):
        """
        Transformacja współrzędnych w odwzorowaniu Gaussa-Krugera 
        do współrzędnych w układzie PL-2000.

        Parameters
        ----------
        x_gk, y_gk : FLOAT
            Współrzędne w odwzorowaniu Gaussa-Krugera
        zone : STR
            Numer strefy odzworowawczej

        Returns
        -------
        x_2000, y_2000 : FLOAT
            Współrzędne w układzie PL-2000

        """
        m0 = 0.999923
        x_2000 = x_gk * m0
        y_2000 = y_gk * m0 + zone * 1000000 + 500000
        return x_2000, y_2000

    def plh2000(self, phi, lam):
        """
            Transformacja współrzędnych geodezyjnych (phi, lam), na współrzędne w układzie PL-2000 (x_2000, y_2000).

        Parameters
        ----------
        Phi: FLOAT
            Szerokość geodezyjna
        Lam: FLOAT
            Długość geodezyjna.

        Returns
        -------
        x_2000, y_2000 : FLOAT
            Współrzędne w układzie PL-2000

        """

        strefa = 0
        lam0 = 0
        center_zone_lam = (degrees(lam) // 3)

        if degrees(lam) % 3 == 0.0:
            left_zone_boundary = center_zone_lam * 3 - 1.5
            right_zone_boundary = center_zone_lam * 3 + 1.5

            if radians(left_zone_boundary) <= lam < radians(right_zone_boundary):
                lam0 = radians(center_zone_lam * 3)
                strefa = degrees(lam0) / 3

            x_gk, y_gk = self.plh2gk(phi, lam, lam0)
        else:
            if degrees(lam) - center_zone_lam * 3 < 1.5:
                left_zone_boundary = ((center_zone_lam - 1) * 3 + 1.5)
                right_zone_boundary = center_zone_lam * 3 + 1.5
                if radians(left_zone_boundary) <= lam < radians(right_zone_boundary):
                    lam0 = radians(center_zone_lam * 3)
                    strefa = degrees(lam0) / 3
            else:
                left_zone_boundary = center_zone_lam * 3 + 1.5
                right_zone_boundary = ((center_zone_lam + 1) * 3 + 1.5)
                if radians(left_zone_boundary) <= lam < radians(right_zone_boundary):
                    lam0 = radians((center_zone_lam + 1) * 3)
                    strefa = degrees(lam0) / 3

            x_gk, y_gk = self.plh2gk(phi, lam, lam0)

        x_2000, y_2000 = self.gk2000(x_gk, y_gk, strefa)
        return x_2000, y_2000

    
if __name__ == "__main__":

    # podanie długości nagłówka oraz modelu elipsoidy
    for i in range(len(sys.argv)):
        if sys.argv[i] == '--naglowek':
            if i + 1 < len(sys.argv):
                naglowek = int(sys.argv[i + 1])

        if sys.argv[i] == '--model':
            if i + 1 < len(sys.argv):
                model = sys.argv[i + 1]
                geo = Transformacje(model=model)

    # Wartości domyślne nagłówka i modelu
    if '--naglowek' not in sys.argv:
        naglowek = 1

    if '--model' not in sys.argv:
        model = 'wgs84'
        geo = Transformacje(model='wgs84')

    # podanie ścieżki pliku z danymi
    input_file_path = sys.argv[-1]

    # wydruk argumentów
    # print(sys.argv)
    print('Długość nagłówka:', naglowek)
    print('Wybrany model:', model)
    print('Wybrana operacja:' )

    # licznik, który sprawdza, czy nie została podana więcej niż jedna funkcja
    flags = ['--xyz2plh', '--plh2xyz',
             '--xyz2neu',
             '--plh1992', '--plh2000',
             '--xyz1992', '--xyz2000']  # lista flag-funkcji
    count_flags = 0
    chosen_flag = None  # zmienna przechowująca wybraną flagę
    
    for flag in flags:
        if flag in sys.argv:
            count_flags += 1
            chosen_flag = flag
    
    if count_flags > 1:
        raise ValueError('Możesz podać tylko jedną flagę.')
    elif count_flags == 1:
        print(chosen_flag)
    
    elif count_flags == 0:
        raise ValueError('Nie wybrano żadnej operacji')
  
    # --xyz2plh

    if '--xyz2plh' in sys.argv:

        with open(input_file_path, 'r') as f:
            lines = f.readlines()
            coords_lines = lines[naglowek:]

            coords_plh = []

            for coord_line in coords_lines:
                coord_line = coord_line.strip('\n')
                x_str, y_str, z_str = coord_line.split(',')
                x, y, z = (float(x_str), float(y_str), float(z_str))
                phi, lam, h = geo.xyz2plh(x, y, z)
                coords_plh.append([phi, lam, h])

        with open('result_xyz2plh.txt', 'w') as f:
            f.write('phi [deg], lam [deg], h [m]\n')

            for coords_list in coords_plh:
                line = ','.join([str(coord) for coord in coords_list])
                f.writelines(line + '\n')

    # --plh2xyz

    elif '--plh2xyz' in sys.argv:

        with open(input_file_path, 'r') as f:
            lines = f.readlines()
            coords_lines = lines[naglowek:]

            coords_xyz = []

            for coord_line in coords_lines:
                coord_line = coord_line.strip('\n')
                phi_str, lam_str, h_str = coord_line.split(',')
                phi, lam, h = (float(phi_str,), float(lam_str), float(h_str))
                x, y, z = geo.plh2xyz(phi, lam, h)
                coords_xyz.append([x, y, z])

        with open('result_plh2xyz.txt', 'w') as f:
            f.write('X [m], Y [m], Z [m]\n')

            for coords_list in coords_xyz:
                line = ','.join([str(coord) for coord in coords_list])
                f.writelines(line + '\n')

    # --xyz2neu

    elif '--xyz2neu' in sys.argv:

        with open(input_file_path, 'r') as f:
            lines = f.readlines()
            coords_lines = lines[naglowek:]

            coords_plh = []

            for coord_line in coords_lines:
                coord_line = coord_line.strip('\n')
                x_str, y_str, z_str = coord_line.split(',')
                x, y, z = (float(x_str), float(y_str), float(z_str))
                x_0, y_0, z_0 = [float(coord) for coord in sys.argv[-4:-1]]
                n, e, u = geo.xyz2neu(x, y, z, x_0, y_0, z_0)
                coords_plh.append([n, e, u])

        with open('result_xyz2neu.txt', 'w') as f:
            f.write('n [m], e [m], u [m]\n')

            for coords_list in coords_plh:
                line = ','.join([str(coord) for coord in coords_list])
                f.writelines(line + '\n')

    # --plh1992

    elif '--plh1992' in sys.argv:

        with open(input_file_path, 'r') as f:
            lines = f.readlines()
            coords_lines = lines[naglowek:]

            coords_xy = []

            for coord_line in coords_lines:
                coord_line = coord_line.strip('\n')
                phi_str, lam_str, h_str = coord_line.split(',')
                phi, lam, h = (radians(float(phi_str)), radians(float(lam_str)), float(h_str))
                x_1992, y_1992 = geo.plh1992(phi, lam)
                coords_xy.append([x_1992, y_1992])

        with open('result_plh1992.txt', 'w') as f:
            f.write('x [m], y [m]\n')

            for coords_list in coords_xy:
                line = ','.join([str(coord) for coord in coords_list])
                f.writelines(line + '\n')

    # --plh2000

    elif '--plh2000' in sys.argv:

        with open(input_file_path, 'r') as f:
            lines = f.readlines()
            coords_lines = lines[naglowek:]

            coords_xy = []

            for coord_line in coords_lines:
                coord_line = coord_line.strip('\n')
                phi_str, lam_str, h_str = coord_line.split(',')
                phi, lam, h = (radians(float(phi_str,)), radians(float(lam_str)), float(h_str))
                x_2000, y_2000 = geo.plh2000(phi, lam)
                coords_xy.append([x_2000, y_2000])

        with open('result_plh2000.txt', 'w') as f:
            f.write('x [m], y [m]\n')

            for coords_list in coords_xy:
                line = ','.join([str(coord) for coord in coords_list])
                f.writelines(line + '\n')

    # --xyz1992

    elif '--xyz1992' in sys.argv:

        with open(input_file_path, 'r') as f:
            lines = f.readlines()
            coords_lines = lines[naglowek:]

            coords_xy = []

            for coord_line in coords_lines:
                coord_line = coord_line.strip('\n')
                x_str, y_str, z_str = coord_line.split(',')
                x, y, z = (float(x_str), float(y_str), float(z_str))
                phi, lam, h = [radians(coord) for coord in geo.xyz2plh(x, y, z)]
                x_1992, y_1992 = geo.plh1992(phi, lam)
                coords_xy.append([x_1992, y_1992])

        with open('result_xyz1992.txt', 'w') as f:
            f.write('x [m], y [m]\n')

            for coords_list in coords_xy:
                line = ','.join([str(coord) for coord in coords_list])
                f.writelines(line + '\n')

    # --xyz2000

    elif '--xyz2000' in sys.argv:

        with open(input_file_path, 'r') as f:
            lines = f.readlines()
            coords_lines = lines[naglowek:]

            coords_xy = []

            for coord_line in coords_lines:
                coord_line = coord_line.strip('\n')
                x_str, y_str, z_str = coord_line.split(',')
                x, y, z = (float(x_str), float(y_str), float(z_str))
                phi, lam, h = [radians(coord) for coord in geo.xyz2plh(x, y, z)]
                x_2000, y_2000 = geo.plh2000(phi, lam)
                coords_xy.append([x_2000, y_2000])

        with open('result_xyz2000.txt', 'w') as f:
            f.write('x [m], y [m]\n')

            for coords_list in coords_xy:
                line = ','.join([str(coord) for coord in coords_list])
                f.writelines(line + '\n')