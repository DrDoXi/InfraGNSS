# -*- coding: utf-8 -*-
"""
Crée on Tue Jan  1 09:48:38 2021

@author: Driss L'hamdouchi , lhamdouchidrisss2@gmail.com

Calcul de la position des satellites GPS à partir du fichier de navigation Rinex
 et traçage de leurs positions dans une période de temps spécifiée par l'utilisateur

Test on Rinex v2.10, v3.02 with GPS navigation file

Requirements:

"""
#--------------------------------------------------Parse_rinex_navigation_file-------------------------------------

""" 
Cette partie du code s'intéresse à la fonction read_rinex
qui prend comme élément d'entrée le fichier de navigation rinex 
et génère comme élément de sortie un objet dictionnaire contenant toutes les données, 
ce qui va permettra ensuite d'accéder à ces données pour le traitement 
et les calculs des positions qui suivent
"""

from io import TextIOWrapper

class EndOfFile(Exception):
    pass

class ErrorOBSRecord(Exception):
    pass

def _split_neg_num(number, start_index=0):
    index_minus = number.find('-', start_index)
    fixed = []

    if index_minus > 0 and not number[index_minus-1].isalpha():
        num1 = number[:index_minus]
        num2 = number[index_minus:]
        fixn1 = _split_neg_num(num1)
        fixn2 = _split_neg_num(num2)
        
        if fixn1 is None:
            fixed.append(num1)
        else:
            for i in fixn1:
                fixed.append(i)

        if fixn2 is None:
            fixed.append(num2)
        else:
            for i in fixn2:
                fixed.append(i)

        return fixed

    else:    
        if index_minus != -1:
            return _split_neg_num(number, index_minus+1)
        else:
            return None

def _fix_negative_num(nums: list) -> list:
    fixed_nums = []
    for num in nums:
        fixed_num = _split_neg_num(num)
        if fixed_num is not None:
            for fn in fixed_num:
                fixed_nums.append(fn)
        else:
            fixed_nums.append(num)

    return fixed_nums

def skip_header(rinex_file: TextIOWrapper) -> int:
    noLine = 0
    while True:
        noLine += 1
        if 'END' in rinex_file.readline():
            break

    return noLine

def read_PRN_EPOCH_SV_CLK(nums: list) -> dict:
    if len(nums) != 10:
        raise ErrorOBSRecord(f'PRN_EPOCH_SV_CLK read error str: {str(nums)}')
    
    return {
        'PRN': nums[0],
        'EPOCH': {
            'YEAR': nums[1],
            'MONTH': nums[2],
            'DAY': nums[3],
            'HOUR': nums[4],
            'MINUTE': nums[5],
            'SECOND': nums[6],
        },
        'SV_clock_bias':float(nums[7]),
        'SV_clock_drift':float(nums[8]),
        'SV_clock_drift_rate':float(nums[9])
    }
    
def read_BROADCAST_ORBIT_1(nums: list) -> dict:
    if len(nums) != 4:
        raise ErrorOBSRecord(f'BROADCAST_ORBIT_1 read error str: {str(nums)}')

    return {
        'IODE': float(nums[0]),
        'Crs': float(nums[1]),
        'Delta_n': float(nums[2]),
        'M0': float(nums[3])
    }

def read_BROADCAST_ORBIT_2(nums: list) -> dict:
    if len(nums) != 4:
        raise ErrorOBSRecord(f'BROADCAST_ORBIT_2 read error str: {str(nums)}')

    return {
        'Cuc': float(nums[0]),
        'e_Eccentricity': float(nums[1]),
        'Cus': float(nums[2]),
        'sqrt_A': float(nums[3])
    }

def read_BROADCAST_ORBIT_3(nums: list) -> dict:
    if len(nums) != 4:
        raise ErrorOBSRecord(f'BROADCAST_ORBIT_3 read error str: {str(nums)}')

    return {
        'Toe': float(nums[0]),
        'Cic': float(nums[1]),
        'OMEGA': float(nums[2]),
        'Cis': float(nums[3])
    }

def read_BROADCAST_ORBIT_4(nums: list) -> dict:
    if len(nums) != 4:
        raise ErrorOBSRecord(f'BROADCAST_ORBIT_4 read error str: {str(nums)}')

    return {
        'i0': float(nums[0]),
        'Crc': float(nums[1]),
        'omega': float(nums[2]),
        'OMEGA_DOT': float(nums[3])
    }

def read_BROADCAST_ORBIT_5(nums: list) -> dict:
    if len(nums) != 4:
        raise ErrorOBSRecord(f'BROADCAST_ORBIT_5 read error str: {str(nums)}')

    return {
        'IDOT': float(nums[0]),
        'Codes_L2_channel': float(nums[1]),
        'GPS_week': float(nums[2]),
        'L2_P': float(nums[3])
    }

def read_BROADCAST_ORBIT_6(nums: list) -> dict:
    if len(nums) != 4:
        raise ErrorOBSRecord(f'BROADCAST_ORBIT_6 read error str: {str(nums)}')

    return {
        'SV_accuracy': float(nums[0]),
        'SV_health': float(nums[1]),
        'TGD': float(nums[2]),
        'IODC': float(nums[3])
    }

def read_BROADCAST_ORBIT_7(nums: list) -> dict:
    if len(nums) != 2:
        raise ErrorOBSRecord(f'BROADCAST_ORBIT_7 read error str: {str(nums)}')

    return {
        'TTM': float(nums[0]),
        'Fit_interval': float(nums[1])
    }

def _next_line(rinex_file: TextIOWrapper) -> list:
    line = rinex_file.readline()
    if not line or line.isspace():
        raise EndOfFile

    nums = [num for num in line.strip().replace('D', 'e').split(' ') if num != '']
    fixed_nums = _fix_negative_num(nums)

    return fixed_nums

def _extract_data(rinex_file: TextIOWrapper) -> dict:    
    ext_data = {}
    nr_sat = 0
    nr_line = 0
    while True:
        try:
            str_data = []

            for _ in range(8):
                nr_line += 1
                data_from_string = _next_line(rinex_file)
                str_data.append(data_from_string)
                
            ex_data_l1 = read_PRN_EPOCH_SV_CLK(str_data[0])
            key = f"{nr_sat}_{str(ex_data_l1['PRN'])}"
            nr_sat += 1

            ext_data[key] = ex_data_l1
            ext_data[key].update(read_BROADCAST_ORBIT_1(str_data[1]))
            ext_data[key].update(read_BROADCAST_ORBIT_2(str_data[2]))
            ext_data[key].update(read_BROADCAST_ORBIT_3(str_data[3]))
            ext_data[key].update(read_BROADCAST_ORBIT_4(str_data[4]))
            ext_data[key].update(read_BROADCAST_ORBIT_5(str_data[5]))
            ext_data[key].update(read_BROADCAST_ORBIT_6(str_data[6]))
            ext_data[key].update(read_BROADCAST_ORBIT_7(str_data[7]))

        except EndOfFile:
            break
        
        except ErrorOBSRecord as eobsr:
            print(f'Error: OBS Record {nr_sat}, Data: {str_data}, NoLine: {nr_line}', eobsr)
            break
        
    return ext_data

def read_rinex(filename: str) -> dict:
    ext_data = None
    with open(filename, 'r') as rinex_file:
        skipped_lines = skip_header(rinex_file)
        ext_data = _extract_data(rinex_file)
    return ext_data

#--------------------------------------------------Parse_rinex_observation_file--------------------------------------------------
"""
Lecture de l'en-tête du fichier rinex d'observation 
et extraction de la position approximative de la station permanente en XYZ 
dans le référentiel Earth Centered Earth Fixed ECEF
 afin d'utiliser cette information lors du traçage des orbites des satellites

"""
def RINEXREADER(file):
    head = True
    header = {}
    with open(file) as f:
        
        for line in f:
            if not ((line[0:36].strip()== 'other post-header comments skipped') or (line[28:34].strip()== '4  1')):
                if head:
                    lines=(line[0:60],line[60:])
                    HT= HTYPER((lines[1]))
                    header = ASSIGNDIC(header,lines[0],HT)
                    if HT == 20:
                        head = False
                else:
                    continue

    f.close()
    fielddict_file = open("obsHead_data.txt","w")
    pprint(header, fielddict_file)
    fielddict_file.close()
    return header

def HTYPER(ss): 
    ss = ss.strip().split()   
    if ss[len(ss)-1] == 'TYPE':
        if ss[len(ss)-3] == 'VERSION':
            return 1
        else:
            return 8
    elif ss[len(ss)-1] == 'COMMENT':
        return 2
    elif ss[len(ss)-1] == 'DATE':
        return 3
    elif ss[len(ss)-1] == 'NAME':
        return 4
    elif ss[len(ss)-1] == 'NUMBER':
        return 5
    elif ss[len(ss)-1] == 'AGENCY':
        return 6
    elif ss[len(ss)-1] == 'VERS':
        return 7
    elif ss[len(ss)-1] == 'XYZ':
        return 9
    elif ss[len(ss)-1] == 'H/E/N':
        return 10
    elif ss[len(ss)-1] == 'L1/2':
        return 11
    elif ss[len(ss)-1] == 'OBSERV':
        return 12
    elif ss[len(ss)-1] == 'INTERVAL':
        return 13
    elif ss[len(ss)-1] == 'OBS':
        if ss[len(ss)-2] == 'FIRST':
            return 14
        elif ss[len(ss)-2] == 'LAST':
            return 15
        else:
            return 19
    elif ss[len(ss)-1] == 'APPL':
        return 16
    elif ss[len(ss)-1] == 'SECONDS':
        return 17
    elif ss[len(ss)-1] == 'SATELLITES':
        return 18
    elif ss[len(ss)-1] == 'HEADER':
        return 20
    else:
        print ("ERROR, DOES NOT EXIST")

def ASSIGNDIC(HE, S, HT):
    
    """
This function ASSIGNDIC() takes in 3 values the main dictionary HE, the current string line being,
and the Header type HT being currently read.
"""  

    if HT == 1:
        RVDT = {}
        RVDT['VER'] = S[0:20].strip()
        RVDT['OBSTYP'] = S[20:40].strip()
        RVDT['SATSYS'] = S[40:60].strip()
        #print RVDT
        HE['RVDT']=RVDT
        return HE
    elif HT == 2:
        if not 'COMMENT' in HE.keys():
            HE['COMMENT']= S.strip()
        else: 
            t = HE['COMMENT']
            
            HE['COMMENT']=t+'\n'+S
        return HE
    elif HT == 3:
        PRBD = {}
        PRBD['PGEN']=S[0:20].strip()
        PRBD['RUNBY']=S[20:40].strip()
        PRBD['DATE']=S[40:60].strip()
        HE['PRBD'] = PRBD
        #print PRBD #get rid later
        return HE
    elif HT == 4:
        HE['MRKR'] = S.strip()
        return HE
    elif HT == 5:
        HE['MKNUM'] = S.strip()
        return HE
    elif HT == 6:
        DONEBY = {}
        DONEBY['OBSV'] = S[0:20].strip()
        DONEBY['AGEN'] = S[20:60].strip()
        HE['DONEBY'] = DONEBY
        #print DONEBY
        return HE    
    elif HT == 7:
        RECV = {}
        RECV['NUM'] = S[0:20].strip()
        RECV['TYP'] = S[20:40].strip()
        RECV['VERS'] = S[40:60].strip()
        #print RECV
        HE['RECV'] = RECV
        return HE
    elif HT == 8:
        ANT = {} 
        ANT['NUM'] = S[0:20].strip()
        ANT['TYP'] = S[20:40].strip()
        #print ANT
        HE['ANT'] = ANT   
        return HE
    elif HT == 9:
        POS = {}
        POS['X'] = float(S[0:15].strip())
        POS['Y'] = float(S[15:30].strip())
        POS['Z'] = float(S[30:45].strip())
        #print POS
        HE['POS'] = POS
        return HE
    elif HT == 10:
        ANTDEL = {}
        ANTDEL['HT'] = S[0:15].strip()
        ANTDEL['EAEC'] = S[15:30].strip()
        ANTDEL['NOEC'] = S[30:45].strip()
        #print ANTDEL
        HE['ANTDEL'] = ANTDEL
        return HE
    elif HT == 12:
        if not 'OBSTYP' in HE.keys():
            OBSTYP = {}
            OBSTYP['NUM']=int(S[0:6])            
            lt=[]
            
            if(OBSTYP['NUM']>=8):
                for x in range(0, 8):
                    lt.append(S[6+x*6:12+x*6].strip())
            else:
                for x in range(0, OBSTYP['NUM']):
                    lt.append(S[6+x*6:12+x*6].strip())
            OBSTYP['OBS']=lt
            HE['OBSTYP']=OBSTYP
        else: 
            t = HE['OBSTYP']
            for x in range(0, t['NUM']-8):
                    t['OBS'].append(S[6+x*6:12+x*6].strip())        
            HE['OBSTYP']=t
        return HE
    elif HT == 14:
        TFIRST = {}
        TFIRST['YEAR'] = S[0:6].strip()
        TFIRST['MON'] = S[6:12].strip()
        TFIRST['DAY'] = S[12:18].strip()
        TFIRST['HR'] = S[18:24].strip()
        TFIRST['MIN'] = S[24:30].strip()
        TFIRST['SEC'] = S[30:43].strip()
        TFIRST['TS'] = S[43:51].strip()
        # print TFIRST
        HE['TFIRST'] = TFIRST
        return HE   
    else:
        return HE

#--------------------------------------------------Satpos--------------------------------------------------
"""

fonctions mathématiques pour calculer la position des satellites à partir du fichier rinex
navigation , la fonction calculate_satpos prend comme argument l'objet dictionnaire généré
par la fonction read_rinex définie ci-dessus

"""
import math

GM = 3.986005 * 10**14
OMEGA_e = 7.292115 * 10**(-5)

def _calculate_tk(t, toe) -> float:
    """
    Calculer le temps tk à partir de l'époque de référence des 
    éphémérides toe (t et toe sont exprimés en secondes dans la semaine GPS) :

    """
    tk = t - toe
    if tk > 302400.0:
        tk = tk - 604800.0
    elif tk < -302400.0:
        tk = tk + 604800.0
    return tk

def _calculate_Ek(Mk, e):
    """
    Résoudre (itérativement) l'équation de Kepler 
    pour l'anomalie d'excentricité Ek :

    """
    Ek = Mk
    temp = Ek
    while math.fabs(Ek-temp) >= 1e-10:
        temp = Ek
        Ek = Mk + e*math.sin(Ek)
    return Ek

def calculate_satpos(sat_rinex: dict,t) -> tuple:

    #-----déclaration des variables----- 
    A = sat_rinex['sqrt_A']**2
    n_0 = math.sqrt(GM/ (A**3))
    n = n_0 + sat_rinex['Delta_n']
    e = sat_rinex['e_Eccentricity']
    tk = _calculate_tk(t,sat_rinex['Toe'])
    Mk = sat_rinex['M0'] + n * tk
    Ek = _calculate_Ek(Mk, e)
    
    #----Calculez l’anomalie vraie vk :----
    vk = math.atan2(math.sqrt(1-e*e) * math.sin(Ek), math.cos(Ek) - e)

    #----Calculer l'argument de latitude uk à partir de l'argument du périgée ω,
    #----de l'anomalie vraie vk et des corrections cuc et cus :
    phi_k = vk + sat_rinex['omega']
    d_uk = sat_rinex['Cuc'] * math.cos(2*phi_k) + sat_rinex['Cus'] * math.sin(2*phi_k)
    uk = phi_k + d_uk

    #----Calculer la distance radiale rk en tenant compte des corrections crc et crs :
    d_rk = sat_rinex['Crc'] * math.cos(2*phi_k) + sat_rinex['Crs'] * math.sin(2*phi_k)
    rk = A * (1 - e * math.cos(Ek)) + d_rk

    #----Calculer l'inclinaison ik du plan orbital à partir de l'inclinaison io au temps de référence toe,
    #----et des corrections cic et cis :
    d_ik = sat_rinex['Cic'] * math.cos(2*phi_k) + sat_rinex['Cis'] * math.sin(2*phi_k)
    ik = sat_rinex['i0'] + d_ik + sat_rinex['IDOT'] * tk

    #----Calculer la longitude du nœud ascendant λk (par rapport à Greenwich).
    #----Ce calcul utilise l'ascension droite du début de la semaine en cours (Ωo),
    #----la correction de la variation apparente du temps sidéral à Greenwich 
    #----entre le début de la semaine et l'heure de référence tk=t−toe,
    #----et le changement de longitude de l'ascendant nœud à partir de la pointe de temps de référence :
    omega_k = sat_rinex['OMEGA'] + (sat_rinex['OMEGA_DOT'] - OMEGA_e) * tk - OMEGA_e * sat_rinex['Toe'] 

    #----Calculer les coordonnées dans le référentiel CTS 
    #----en appliquant trois rotations (autour de uk, ik et λk) :

    #----rotations autour de uk :
    xk_prim = rk * math.cos(uk)
    yk_prim = rk * math.sin(uk)

    #----rotations autour de ik et λk :
    xk = xk_prim * math.cos(omega_k) - yk_prim * math.cos(ik) * math.sin(omega_k)
    yk = xk_prim * math.sin(omega_k) + yk_prim * math.cos(ik) * math.cos(omega_k)
    zk = yk_prim * math.sin(ik)

    return (xk, yk, zk)

def calculate_positions(sat_data:dict,Tstart,Tfinish,intervale):
    sat_pos = {}

    for key, sat_rinex in sat_data.items(): 
        for i in range(int(((getSecs(Tfinish)-getSecs(Tstart))//intervale))):
            xk, yk, zk = calculate_satpos(sat_rinex,getSecs(Tstart)+i*intervale)
            sat_pos[key+'_At epoch  '+str(weeksecondstoutc(sat_rinex['GPS_week'],getSecs(Tstart)+i*intervale,0))] = {
            'x': xk,
            'y': yk,
            'z': zk
            }
    return sat_pos

def select_best_ephemeride(sat_data:dict,Tmoy):

    sat_data_selected={}
    
    sat_in_sat_data=[]

    Sat_index_near_time=[]

    for each in sat_data.keys() :
        if each.split('_')[1] not in sat_in_sat_data :
            sat_in_sat_data.append(each.split('_')[1])

    # print(sat_in_sat_data)


    for satellite_num in sat_in_sat_data:
        L=[]
        for key, sat_rinex in sat_data.items():
            
            if key.split('_')[1]== satellite_num :

                {key:abs(sat_rinex['Toe']-Tmoy)}
 
                L.append([key,abs(sat_rinex['Toe']-Tmoy)])
        
        min= 0

        for i in range(len(L)):
            if L[i][1] < L[min][1]:
                min = i

        Sat_index_near_time.append(L[min][0])
           
    # print(Sat_index_near_time)

    for ele in Sat_index_near_time :
        try :
          sat_data_selected[ele] = sat_data[ele]
        except :
            continue
    
    return sat_data_selected   

#--------------------------------------------------getSecs----------------------------------------------------

from math import floor,fmod
import datetime

def getSecs(epoch):
    """
    renvoie les secondes GPS de la semaine avec la date entrée

    """
    t = epoch.split(":")
    yy = 2000+float(t[0])
    mm = float(t[1])
    dd = float(t[2])
    hh = float(t[3])
    mins = float(t[4])
    ss = float(t[5])
    hh = hh+(mins/60)+(ss/3600)
    if mm <= 2:
        yy = yy-1
        mm = mm+12
    jd = floor(365.25*(yy+4716))+floor(30.6001*(mm+1))+dd+(hh/24)-1537.5
    a = floor(jd+0.5)
    b = a+1537
    c = floor((b-122.1)/365.25)
    e = floor(365.25*c)
    f = floor((b-e)/30.6001)
    d = b-e-floor(30.6001*f)+(fmod((jd+0.5),1))
    day_of_week = fmod(floor(jd+0.5),7)
    secs = (fmod(d,1)+day_of_week+1)*86400
    return secs

def weeksecondstoutc(gpsweek,gpsseconds,leapseconds):

    datetimeformat = "%Y-%m-%d %H:%M:%S"
    epoch = datetime.datetime.strptime("1980-01-06 00:00:00",datetimeformat)
    elapsed = datetime.timedelta(days=(gpsweek*7),seconds=(gpsseconds-leapseconds))
    return datetime.datetime.strftime(epoch + elapsed,datetimeformat)

#--------------------------------------------------Main----------------------------------------------------
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from pymap3d import ecef2geodetic
import numpy as np
from pprint import pprint

import csv



def show(sat_pos:dict,sat_data,n,obs_data):
    
    sat=sat_pos.keys()
    for each in sat :
        if each.split('_')[1]==str(n) : 
                print('satellite : ','G'+each.split('_')[1]+' '+each.split('_')[2],'x :',sat_pos[each]['x'],'y :',sat_pos[each]['y'],'z :',sat_pos[each]['z'])

def plotsat(sat_pos,n,fileObs):
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    xp =[]
    yp =[]
    zp =[]


    phi = np.linspace(0,2*np.pi, 256).reshape(256, 1) # the angle of the projection in the xy-plane
    theta = np.linspace(0, np.pi, 256).reshape(-1, 256) # the angle from the polar axis, ie the polar angle
    a = 6378137.1
    b = 6356752.314140

    # Transformation formulae for a spherical coordinate system.
    x = a*np.sin(theta)*np.cos(phi)
    y = a*np.sin(theta)*np.sin(phi)
    z = b*np.cos(theta)

    # fig = plt.figure(figsize=plt.figaspect(1))  # Square figure
    ax.plot_surface(x, y, z, color='b')


    sat=sat_pos.keys()

    for each in sat :

        if int(each.split('_')[1]) == n:
              xp.append(sat_pos[each]['x'])
              yp.append(sat_pos[each]['y'])
              zp.append(sat_pos[each]['z'])
            
    xApr=[RINEXREADER(fileObs)['POS']['X']]
    yApr=[RINEXREADER(fileObs)['POS']['Y']]
    zApr=[RINEXREADER(fileObs)['POS']['Z']]
    
    ax.scatter(xp, yp, zp, c='r', marker='o')
    ax.scatter(xApr, yApr, zApr, c='r', marker='o')
    
    xp = np.array(xp)
    yp = np.array(yp)
    zp = np.array(zp)


    max_range = np.array([xp.max()-xp.min(), yp.max()-yp.min(), zp.max()-zp.min()]).max() / 2.0

    mid_x = (xp.max()+xp.min()) * 0.5
    mid_y = (yp.max()+yp.min()) * 0.5
    mid_z = (zp.max()+zp.min()) * 0.5
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range) 



    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    text='la position du satellite  :G'+str(n)
    plt.title(text)

    plt.show()

print('\n------------ Logiciel de calcul de position satellite à partit du message de navigation ------------ ')
print("------------ @author: Driss L'hamdouchi , lhamdouchidrisss2@gmail.com ------------------------------")
menu_options = {
    1: 'Choisir/changer le fichier rinex de navigation  : ',
    2: "Choisir/changer le fichier rinex d'observation  : ",
    3: "Régler Tstart ,Tfinish et l'intervalle",
    4: 'Afficher les satellite',
    5: 'Afficher et tracer la postion des satellites',
    6: 'Exit',
}

def print_menu():

    print('\n')

    for key in menu_options.keys():
        print (key, '--', menu_options[key] )
    
    print('\n')

    try :
        print("Fichier de navigation :",fileNav)
    except :
        print("Fichier de navigation : Pas encore choisi")
    try :
        print("Fichier d'observation :",fileObs)
    except :
        print("Fichier d'observation : Pas encore choisi")
    try:
        print("Tstart :",Tstart)
        print("Tfinish :",Tfinish)
        print("intervale :",intervale)      
    except :
        print("Tstart : Pas encore choisi")
        print("Tfinish : Pas encore choisi")
        print("intervale : Pas encore choisi")

def option1():
     global sat_data,fileNav
     print('\n\nVeuillez écrire le nom du fichier rinex de navigation')
     fileNav=input("Fichier de navigation :  ")
     sat_data = read_rinex(fileNav)
     fielddict_file = open("sat_data_brut.txt","w")
     pprint(sat_data, fielddict_file)
     fielddict_file.close()

def option2():
     global obs_data,fileObs
     print("\n\nVeuillez écrire le nom du fichier rinex d'observation")
     fileObs=input("Fichier d'observation :  ")
     obs_data=RINEXREADER(fileObs) 
     
def option3():
     global sat_pos ,sat_data,Tstart,Tfinish,intervale,Tmoy
     print('Veuillez écrire les époches sous la format suivante 21:12:14:00:00:00')
     Tstart=input(("time start  : "))
     Tfinish=input(("time finish  : "))
     intervale=int(input("intervalle en secondes : "))
     Tmoy=int(((getSecs(Tfinish)-getSecs(Tstart))/2)+getSecs(Tstart))


def option5():
    global sat_pos,sat_data,obs_data,fileObs,Tmoy
    print('\n\nVeuillez écrire le Numéro du satellite :')
    sat_pos = calculate_positions(select_best_ephemeride(sat_data,Tmoy),Tstart,Tfinish,intervale)
    fielddict_file = open("sat_pos.txt","w")
    pprint(sat_pos, fielddict_file)
    fielddict_file.close() 
    with open('output.csv', 'w',newline='') as output:
          writer = csv.writer(output)
          writer.writerow(['Satellite',"époque","X","Y","Z"])
          for key, value in sat_pos.items():
                    
                    writer.writerow(["G"+key.split('_')[1],key.split('_')[2].split('  ')[1], value['x'], value['y'], value['z']])

    n=int(input("Entrez votre choix : ")[1:])
    show(sat_pos,sat_data,n,obs_data)


    plotsat(sat_pos,n,fileObs)

def option4():
    global sat_data
    sat_in_sat_data=[]
    for each in sat_data.keys() :
      if each.split('_')[1] not in sat_in_sat_data :
        sat_in_sat_data.append(each.split('_')[1])

    for ele in  sat_in_sat_data :
        if len(ele)==1 :

          print('G'+'0'+ele,end='//')
        
        else :
            print('G'+ele,end='//')



if __name__=='__main__':
    while(True):
        print_menu()
        option = ''
        print('\n')
        try:
            option = int(input('Entrez votre choix  : '))
        except:
            print('Mauvaise entrée. Veuillez entrer un nombre...')
        #Check what choice was entered and act accordingly
        if option == 1:
           option1()
        elif option == 2:
            option2()
        elif option == 3:
            option3()
        elif option == 4:
            option4()
        elif option == 5:
            option5()
        elif option == 6:
            print('Message de remerciement avant de quitter')
            exit()
        else:
            print('Option invalide. Veuillez entrer un nombre entre 1 et 6.')












