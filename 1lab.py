Здесь удалил и добавил строчку
import matplotlib.pyplot as plt
import numpy as np
import math as m
import matplotlib as mpl
import seaborn as sns
import pandas as pd
sns.set()

Добавил зжесь строчку что бы можно было увидеть изменпения 



def main():
    R = 10 #Дано
    Ptx = 0.01 #Дано излучаемая мощность
    dF = 20 * 10**6 #Дано
    Kn = 3 #Дано
    K = 29 #Дано
    f0 = 2400 #Дано
    D_summ_Eq_bl, D_min_Eq_bl, D_mean_Eq_bl = [],[],[]
    D_summ_Max_Thr, D_min_Max_Thr, D_mean_Max_Thr  = [],[],[]
    D_summ_Prop_Fair, D_min_Prop_Fair, D_mean_Prop_Fair= [],[],[]
    Narray = []
    for i in range(7):
        N = 2**i
        dmin1, dsred1, dsum1  = 0,0,0
        dmin2, dsred2, dsum2  = 0,0,0
        dmin3, dsred3, dsum3  = 0,0,0
        Key = 0
        while Key < 100:
            x, y, rarr = gen(R, N)
            # Считаем С
            Ci = []
            for i in range(len(rarr)):
                L = 10 ** ((20 * m.log(f0, 10) + K * m.log(rarr[i], 10) - 28) / 10)
                SNR = (Ptx / L) / (dF * 340 * (1.38 / 10 ** 23) * Kn)
                C = dF * m.log(1 + SNR, 2)
                Ci.append(C)
            #1 алгоритм
            sum = 0  # сумма С
            for i in range(len(Ci)):
                sum = sum + 1 / Ci[i]
            Di_Eq_Bl = 1 / sum  # суммарная скорость
            #2 алгоритм
            Di_Max_Thr = 0
            for i in range(len(Ci)):
                if Ci[i] > Di_Max_Thr:
                    Di_Max_Thr = Ci[i]
            #3 алгоритм
            Di_Prop_Fair = []
            for i in range(len(Ci)):
                Di_Prop_Fair.append(Ci[i] / len(Ci))
            # Расчитываем величины
            Di_Eq_Bl_min = Di_Eq_Bl
            Di_Eq_Bl_sred = Di_Eq_Bl
            Di_Eq_Bl_summ = Di_Eq_Bl * N
            if N == 1:
                Di_Max_Thr_min = Di_Max_Thr
            else:
                Di_Max_Thr_min = 0
            Di_Max_Thr_sred = Di_Max_Thr / N
            Di_Max_Thr_summ = Di_Max_Thr
            Di_Prop_Fair_min = np.min(Di_Prop_Fair)
            Di_Prop_Fair_sred = np.mean(Di_Prop_Fair)
            Di_Prop_Fair_summ = np.sum(Di_Prop_Fair)
            # Записываем величины в соответствующие массивы
            dmin1 += Di_Eq_Bl_min
            dsred1 += Di_Eq_Bl_sred
            dsum1 += Di_Eq_Bl_summ
            dmin2 += Di_Max_Thr_min
            dsred2 += Di_Max_Thr_sred
            dsum2 += Di_Max_Thr_summ
            dmin3 += Di_Prop_Fair_min
            dsred3 += Di_Prop_Fair_sred
            dsum3 += Di_Prop_Fair_summ
            Key += 1
        D_summ_Eq_bl.append(dsum1/100)
        D_min_Eq_bl.append(dmin1/100)
        D_mean_Eq_bl.append(dsred1/100)
        D_summ_Max_Thr.append(dsum2/100)
        D_min_Max_Thr.append(dmin2/100)
        D_mean_Max_Thr.append(dsred2/100)
        D_summ_Prop_Fair.append(dsum3/100)
        D_min_Prop_Fair.append(dmin3/100)
        D_mean_Prop_Fair.append(dsred3/100)
        Narray.append(N)
    #Графики
    plt.figure(1)
    plt.semilogx(Narray, D_summ_Eq_bl, label='EB')
    plt.semilogx(Narray, D_summ_Max_Thr, label='MT')
    plt.semilogx(Narray, D_summ_Prop_Fair, label='PF')
    plt.legend()
    plt.title('sum')

    plt.figure(2)
    plt.semilogx(Narray, D_min_Eq_bl, label='EB')
    plt.semilogx(Narray, D_min_Max_Thr, label='MT')
    plt.semilogx(Narray, D_min_Prop_Fair, label='PF')
    plt.legend()
    plt.title('min')

    plt.figure(3)
    plt.semilogx(Narray, D_mean_Eq_bl, label='EB')
    plt.semilogx(Narray, D_mean_Max_Thr, label='MT')
    plt.semilogx(Narray, D_mean_Prop_Fair, label='PF')
    plt.title('mean')
    plt.legend()
    plt.show()


def gen(R, N):
    r = []
    xarr = []
    yarr = []
    for i in range(N):
        score = 1
        while score:
            x = rand.uniform(-1 * R, R)
            y = rand.uniform(-1 * R, R)
            if (x ** 2 + y ** 2) < R ** 2:
                score = 0
                r.append((x ** 2 + y ** 2) ** (1 / 2))
                xarr.append(x)
                yarr.append(y)
    return xarr, yarr, r


def paint():
    R = 10
    N = 1500
    x, y, rarr = gen(R, N)
    theta = np.linspace(0, 2 * np.pi, 100)
    xr = R * np.cos(theta)
    yr = R * np.sin(theta)
    plt.figure(1)
    plt.plot(xr, yr)
    sns.scatterplot(x=x, y=y)
    plt.axis('scaled')
    plt.show()

def hist():
    N = 1000000
    Rad = 30
    x, y, R = gen(Rad, N)
    plt.figure(1)
    sns.distplot(R, kde=False, bins=30)
    plt.show()

def emp_func():
    N = 1000
    R = 30
    a1 = 0
    a = []
    b = []
    for i in range(N):
        a1 += R/N
        a.append(a1)
        b.append(a1**2/R**2)
    x, y, r = gen(R, N)
    r = np.sort(r)
    plt.plot(np.sort(r), np.linspace(0,1,len(r), endpoint=False))
    plt.plot(a,b)
    plt.show()
    
#emp_func()
#hist()
#main()
paint()
