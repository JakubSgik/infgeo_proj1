import sys
from math import sin, cos, sqrt, atan, atan2, degrees, radians, tan
from numpy import array
import numpy as np

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
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "mars":
            self.a = 3396900.0
            self.b = 3376097.80585952
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flat = (self.a - self.b) / self.a
        self.ecc = sqrt(2 * self.flat - self.flat ** 2) # eccentricity  WGS84:0.0818191910428 
        self.ecc2 = (2 * self.flat - self.flat ** 2) # eccentricity**2


    
    def xyz2plh(self, X, Y, Z, output = 'dec_degree'):
        """
        Algorytm Hirvonena - algorytm transformacji współrzędnych ortokartezjańskich (x, y, z)
        na współrzędne geodezyjne długość szerokość i wysokośc elipsoidalna (phi, lam, h). Jest to proces iteracyjny. 
        W wyniku 3-4-krotneej iteracji wyznaczenia wsp. phi można przeliczyć współrzędne z dokładnoscią ok 1 cm.     
        Parameters
        ----------
        X, Y, Z : FLOAT
             współrzędne w układzie orto-kartezjańskim, 

        Returns
        -------
        lat
            [stopnie dziesiętne] - szerokość geodezyjna
        lon
            [stopnie dziesiętne] - długośc geodezyjna.
        h : TYPE
            [metry] - wysokość elipsoidalna
        output [STR] - optional, defoulf 
            dec_degree - decimal degree
            dms - degree, minutes, sec
        """
        r   = sqrt(X**2 + Y**2)           # promień
        lat_prev = atan(Z / (r * (1 - self.ecc2)))    # pierwsze przybliilizenie
        lat = 0
        while abs(lat_prev - lat) > 0.000001/206265:    
            lat_prev = lat
            N = self.a / sqrt(1 - self.ecc2 * sin(lat_prev)**2)
            h = r / cos(lat_prev) - N
            lat = atan((Z/r) * (((1 - self.ecc2 * N/(N + h))**(-1))))
        lon = atan(Y/X)
        N = self.a / sqrt(1 - self.ecc2 * (sin(lat))**2);
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
        h : TYPE
            [metry] - wysokość elipsoidalna

        Returns
        -------
        X, Y, Z : FLOAT
             współrzędne w układzie orto-kartezjańskim

        """
        phi = radians(phi)
        lam = radians(lam)
        Rn = self.a / sqrt(1-self.ecc2 * sin(phi)**2)
        q = Rn * self.ecc2 * sin(phi)
        
        X = (Rn + h) * cos(phi) * cos(lam)
        Y = (Rn + h) * cos(phi) * sin(lam)
        Z = (Rn + h) * sin(phi) - q
        return X, Y, Z
    
    def xyz2neu(self, x, y, z):
        x0 = 3664940.500
        y0 = 1409153.590
        z0 = 5009571.170
        
        xyz = array((x, y, z))
        xyz0 = array((x0, y0, z0))
        
        dX = xyz - xyz0
        
        phi = radians(52.09727221841272)
        lam = radians(21.03153333279777)

        R = array([[-sin(phi) * cos(phi), -sin(lam), cos(phi) * cos(lam)],
                      [-sin(phi) * sin(lam), cos(lam), cos(phi) * sin(lam)],
                      [cos(phi), 0, sin(phi)]])
        dx = R.T @ dX
        n, e, u = dx
        return n, e, u

    def xyz2neu2(self, x, y, z, x_0, y_0, z_0):
        phi, lam, _ = [radians(coord) for coord in self.xyz2plh(x, y, z)]
        
        R = np.array([[-sin(lam), -sin(phi)*cos(lam), cos(phi)*cos(lam)],
                      [cos(lam), -sin(phi)*sin(lam), cos(phi)*sin(lam)],
                      [        0,           cos(phi),       sin(phi)]])
        
        xyz_t = np.array([[x - x_0],
                          [y - y_0],
                          [z - z_0]])
        
        [[E], [N], [U]] = R.T @  xyz_t
        
        return N, E, U
        
        
#informacje
# zrobilem dokumentacje do plh2xyz zgodnie ze strona 26 z tego zrodla:
# http://www.geonet.net.pl/images/2002_12_uklady_wspolrz.pdf
#wydaje mi sie natomiast ze x0, y0, z0 to powinien podawac uzytkownik

#ponizej moja proba zrobienia transformacji, oczywiscie nieudolna
    # #BL(GRS80, WGS84, ew. Krasowski) -> 2000
    # def bl22000(self, model, b, l, h, numer_pasa_odwzorowawczego):
    #     #kroki - bl krasowski -> xyz krasowski -> xyz gk -> bl gk -> xy 2000
        
    #     #przeliczanie fi lambda ha -> xyz
    #     #pierwszy wertykał
    #     N = self.a/sqrt(1 - self.ecc2 * sin(b)**2)
    #     X = (N + h) * cos(b) * cos(l)
    #     Y = (N + h) * cos(b) * sin(l)
    #     Z = (N * (1 - self.ecc2) + h) * sin(b)
    #     #tylko teraz czy to ma sens skoro mamy juz funkcje plh2xyz
    #     #http://www.geonet.net.pl/images/2002_12_uklady_wspolrz.pdf
    #     #na tej stronie na str 11 mozna przeczytac ze mozna odrazu przejsc z bl krasowskiego na bl gk zaniedbujac wplyw wysokosci
    #     #tylko jak???
        
    #     #tutaj mamy xyz krasowskiego, musimy przejsc do gk
    #     #to chyba trzeba transformacjami z macierzami, tak jak w cwiczeniu 8 z gw sem 3
        
    #     m0 = 0.999923
    #     x2000 = x_gk * m0
    #     y2000 = y_gk * m0 + numer_pasa_odwzorowawczego * 1000000 + 500000
        
    #     return (x2000, y2000)
        
        
    #     #jeszcze jedna wazna rzecz jest taka ze nie musimy robic tych parametrow elipsoidy w argumentach funkcji, 
    #     #bo modele sa juz zbudowane w funkcji __init__ zatem wystarczy zeby uzytkownik wybieral model
        
    def sigma_p(self, model, f):
        
        A0 = 1 - (self.ecc2/4) - ((3*(self.ecc2**2))/64) - ((5*(self.ecc2**3))/256)
        A2 = (3/8) * (self.ecc2 + (self.ecc2**2)/4 + (15*(self.ecc2**3))/128)
        A4 = (15/256) * (self.ecc2**2 + (3*(self.ecc2**3))/4)
        A6 = (35*(self.ecc2**3))/3072

        sigma = self.a * (A0 * f - A2*sin(2*f) + A4*sin(4*f) - A6*sin(6*f))

        return sigma
    
    def plh2gk(self, model, f, l, l0):
        a2 = self.a**2
        b2 = (a2 * (1 - self.ecc2))
        e_prim_2 = (a2 - b2)/b2
        dl = l - l0
        t = tan(f)
        eta2 = e_prim_2 * (cos(f)**2)
        N = self.a / sqrt(1 - self.ecc2 * sin(f)**2)

        sigma = Transformacje.sigma_p(self, model, f)

        x_gk = sigma + (dl ** 2 / 2) * N * sin(f) * cos(f) * (
            1 + (dl ** 2 / 12) * cos(f) ** 2 *
            (5 - t ** 2 + 9 * eta2 + 4 * eta2 ** 2)
            + (dl ** 4 / 360) * cos(f) ** 4 * (61 - 58 * t **
                                                  2 + t ** 4 + 270 * eta2 - 330 * eta2 * t ** 2)
        )

        y_gk = dl * N * cos(f) * (
            1 + (dl ** 2 / 6) * cos(f) ** 2 * (1 - t ** 2 + eta2)
            + (dl ** 4 / 120) * cos(f) ** 4 * (5 - 18 * t **
                                                  2 + t ** 4 + 14 * eta2 - 58 * eta2 * t ** 2)
        )

        return x_gk, y_gk        

    def gk1992(self, x_gk, y_gk):
        m0 = 0.9993
        x_1992 = x_gk * m0 - 5300000
        y_1992 = y_gk * m0 + 500000
        return x_1992, y_1992    

    def plh1992(self, model, phi, lam):
        phi = radians(phi)
        lam = radians(lam)
        l0 = radians(19)
        x_gk, y_gk = Transformacje.plh2gk(self, model, phi, lam, l0)
        x_1992, y_1992 = Transformacje.gk1992(self, x_gk, y_gk)
        
        return x_1992, y_1992   
    
    def konwersja_gk_2000(x_gk, y_gk, nr_strefy):
        m0 = 0.999923
        x_2000 = x_gk * m0
        y_2000 = y_gk * m0 + nr_strefy * 1000000 + 500000
        return x_2000, y_2000
    
    def konwersja_deg_2000(f, l):
        strefa = 0
        l0 = 0
        if l % np.radians(3) == 0.0:
            podzielone_l = (l // np.radians(3))
            
            lewa_granica_strefy = podzielone_l*3 - 1.5
            prawa_granica_strefy = podzielone_l*3 + 1.5
            
            if np.radians(lewa_granica_strefy) <= l < np.radians(prawa_granica_strefy):
                l0 = np.radians(podzielone_l*3)
                strefa = np.degrees(l0)/3
            
            x_gk, y_gk = konwersja_deg_gk(f, l, l0)
        else:
            podzielone_l = (l // np.radians(3))
        
            lewa_granica_strefy = podzielone_l*3 + 1.5
            prawa_granica_strefy = ((podzielone_l+1)*3 + 1.5)
        
            if np.radians(lewa_granica_strefy) <= l < np.radians(prawa_granica_strefy):
                l0 = np.radians(podzielone_l*3 + 3)
                strefa = np.degrees(l0)/3
            x_gk, y_gk = konwersja_deg_gk(f, l, l0)

        x_2000, y_2000 = konwersja_gk_2000(x_gk, y_gk, strefa)
        return x_2000, y_2000, ((int(np.degrees(l0)), int(strefa)))
if __name__ == "__main__":
    # utworzenie obiektu
    geo = Transformacje(model = "wgs84")
    # dane XYZ geocentryczne
    # X = 3664940.500; Y = 1409153.590; Z = 5009571.170
    # phi, lam, h = geo.xyz2plh(X, Y, Z)
    # print(phi, lam, h)
    # phi, lam, h = geo.xyz2plh2(X, Y, Z)
    # print(phi, lam, h)
    
    print(sys.argv)
    input_file_path = sys.argv[-1]
    
    if '--xyz2plh' in sys.argv and '--phl2xyz' and '--xyz2neu' in sys.argv:
        print('możesz podać tylko jedną falgę')
        
    # --xyz2neu
    
    elif '--xyz2plh' in sys.argv:

        with open(input_file_path, 'r') as f:
            lines = f.readlines()
            coords_lines = lines[1:]
            #print(coords_lines)
            
            coords_plh = []
            
            for coord_line in coords_lines:
                coord_line = coord_line.strip('\n')
                x_str, y_str, z_str = coord_line.split(',')
                x, y, z = (float(x_str), float(y_str), float(z_str))
                phi, lam, h = geo.xyz2plh(x, y, z)
                coords_plh.append([phi, lam, h])
            
            
        with open('result_xyz2plh.txt', 'w') as f:
            f.write('phi[deg], lam[deg], h[m]\n')
            
            for coords_list in coords_plh:
                line = ','.join([str(coord) for coord in coords_list])
                f.writelines(line + '\n')
     
    # --plh2xyz    
     
    elif '--plh2xyz' in sys.argv:

        with open(input_file_path, 'r') as f:
            lines = f.readlines()
            coords_lines = lines[1:]
            #print(coords_lines)
            
            coords_xyz = []
            
            for coord_line in coords_lines:
                coord_line = coord_line.strip('\n')
                phi_str, lam_str, h_str = coord_line.split(',')
                phi, lam, h = (float(phi_str,), float(lam_str), float(h_str))
                x, y, z = geo.plh2xyz(phi, lam, h)
                coords_xyz.append([x, y, z])
            
            
        with open('result_plh2xyz.txt', 'w') as f:
            f.write('x[m], y[m], z[m]\n')
            
            for coords_list in coords_xyz:
                line = ','.join([str(coord) for coord in coords_list])
                f.writelines(line + '\n')
                
                
    elif '--xyz2neu' in sys.argv:

        with open(input_file_path, 'r') as f:
            lines = f.readlines()
            coords_lines = lines[1:]
            #print(coords_lines)
            
            coords_plh = []
            
            for coord_line in coords_lines:
                coord_line = coord_line.strip('\n')
                x_str, y_str, z_str = coord_line.split(',')
                x, y, z = (float(x_str), float(y_str), float(z_str))
                n,e,u = geo.xyz2neu(x, y, z)
                coords_plh.append([n,e,u])
            
            
        with open('result_xyz2neu.txt', 'w') as f:
            f.write('n [m], e [m], u [m]\n')
            
            for coords_list in coords_plh:
                line = ','.join([str(coord) for coord in coords_list])
                f.writelines(line + '\n')   
                
    # --plh1992

    elif '--plh1992' in sys.argv:
        model = sys.argv[-2]
        with open(input_file_path, 'r') as f:
            lines = f.readlines()
            coords_lines = lines[1:]
            #print(coords_lines)
            
            coords_xy = []
            
            for coord_line in coords_lines:
                coord_line = coord_line.strip('\n')
                phi_str, lam_str, h_str = coord_line.split(',')
                phi, lam, h = (float(phi_str,), float(lam_str), float(h_str))
                # ------------------------------ DO ZMIANY NA PLH->1992
                x, y = geo.plh1992(model, phi, lam)
                coords_xy.append([x, y])
            
            
        with open('result_plh1992.txt', 'w') as f:
            f.write('x[m], y[m]\n')
            
            for coords_list in coords_xy:
                line = ','.join([str(coord) for coord in coords_list])
                f.writelines(line + '\n')
                
    elif '--xyz2neu2' in sys.argv:

        with open(input_file_path, 'r') as f:
            lines = f.readlines()
            coords_lines = lines[1:]
            #print(coords_lines)
            
            coords_plh = []
            
            for coord_line in coords_lines:
                coord_line = coord_line.strip('\n')
                x_str, y_str, z_str = coord_line.split(',')
                x, y, z = (float(x_str), float(y_str), float(z_str))
                x_0, y_0, z_0 = [float(coord) for coord in sys.argv[-4:-1]]
                n,e,u = geo.xyz2neu2(x, y, z, x_0, y_0, z_0)
                coords_plh.append([n,e,u])
            
            
        with open('result_xyz2neu2.txt', 'w') as f:
            f.write('n [m], e [m], u [m]\n')
            
            for coords_list in coords_plh:
                line = ','.join([str(coord) for coord in coords_list])
                f.writelines(line + '\n')   
        
    
