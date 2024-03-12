import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.constants import g
matplotlib.rc('font', size=12)
plt.figure(figsize=(10,8))

L = 30      # x меняется на отрезке [-L, L]
T = 5       # t меняется на [0, T]
N = 100     # число узлов сетки по x
M = 100     # число узлов сетки по t
h = 2*L/N   # длина шага по x
tau = T/M   # длина шага по t
x = np.linspace(-L, L, N + 1)
t = np.linspace(0, T, M + 1)

# Входные данные: u - скорость, H - высота
# L - слева от нуля, R - справа от нуля
u_L = 0
u_R = 0
H_L = 3
H_R = 1

# Начальные профили u и H, они имеют вид ступеньки
u_0 = (x <= 0)*u_L + (x > 0)*u_R
H_0 = (x <= 0)*H_L + (x > 0)*H_R

#plt.plot(x, H_0, label='t = 0.0', linewidth=2)

# Вычисляем инварианты Римана, они тоже имеют вид ступеньки
R_L = u_L + 2*np.sqrt(g*H_L)
R_R = u_R + 2*np.sqrt(g*H_R)
Q_L = u_L - 2*np.sqrt(g*H_L)
Q_R = u_R - 2*np.sqrt(g*H_R)
R = u_0 + 2*np.sqrt(g*H_0)
Q = u_0 - 2*np.sqrt(g*H_0)
Q[N//2] = Q[N//2 + 1]

# Подготовим массивы

# R,Q в серединах ячеек сетки на данном временном слое (n)
# То есть считается, что R_half[i] - это как бы R[i + 1/2]
R_half = np.concatenate((R_L*np.ones(N//2), R_R*np.ones(N//2)))
Q_half = np.concatenate((Q_L*np.ones(N//2), Q_R*np.ones(N//2)))

# R,Q в центрах расчётных ячеек, (n+1/2)-й временной слой
R_center = np.zeros(N)
Q_center = np.zeros(N)

# R,Q на верхнем временном слое (n+1) в узлах сетки
R_next = np.zeros(N + 1)
Q_next = np.zeros(N + 1)

# R,Q на верхнем временном слое (n+1) в центрах ячеек
R_half_next = np.zeros(N)
Q_half_next = np.zeros(N)


for n in range(M):
    # это что-то типа граничных условий слева и справа
    R_next[ 0] = R_L
    Q_next[-1] = Q_R
    
    # вычисляем лямбды, т.е. скорости переноса инвариантов Римана R,Q
    lambda1 = (3*R_half + Q_half)/4
    lambda2 = (R_half + 3*Q_half)/4
    
    # первая и вторая фазы схемы Кабаре для R
    for i in range(N):
        # первая фаза схемы Кабаре, вычисление значения в центре ячейки
        R_center[i] = R_half[i] - tau/2*lambda1[i]*(R[i+1] - R[i])/h
        # вторая фаза, экстраполяция, 
        # вычисление узлового значения на следующем временном слое
        R_next[i + 1] = 2*R_center[i] - R[i]
        
        # нелинейная коррекция потоков
        R_max = np.max([R[i], R_half[i], R[i+1]])
        R_min = np.min([R[i], R_half[i], R[i+1]])
        if R_next[i + 1] > R_max:
            R_next[i + 1] = R_max
        if R_next[i + 1] < R_min:
            R_next[i + 1] = R_min
            
    # то же самое для Q, но в другом направлении, т.к. lambda2 < 0
    for i in range(N - 1, -1, -1):
        # первая фаза
        Q_center[i] = Q_half[i] - tau/2*lambda2[i]*(Q[i+1] - Q[i])/h
        # вторая фаза
        Q_next[i] = 2*Q_center[i] - Q[i + 1]
        
        # нелинейная коррекция потоков
        Q_max = np.max([Q[i], Q_half[i], Q[i+1]])
        Q_min = np.min([Q[i], Q_half[i], Q[i+1]])
        if Q_next[i] > Q_max:
            Q_next[i] = Q_max
        if Q_next[i] < Q_min:
            Q_next[i] = Q_min
    
    for i in range(N):
        # вычисляем новые лямбды для данной ячейки
        lambda1_center = (3*R_center[i] + Q_center[i])/4
        lambda2_center = (R_center[i] + 3*Q_center[i])/4
        # третья фаза схемы Кабаре для R и Q
        R_half_next[i] = R_center[i] - tau/2*lambda1_center*(R_next[i + 1] - R_next[i])/h
        Q_half_next[i] = Q_center[i] - tau/2*lambda2_center*(Q_next[i + 1] - Q_next[i])/h
    
    # копируем данные для следующей итерации
    R = R_next.copy()
    Q = Q_next.copy()
    R_half = R_half_next.copy()
    Q_half = Q_half_next.copy()
    
    # выводим некоторые промежуточные графики
    #if n > 0 and n%100 == 0:
    if n == M-1:
        # обратно вычисляем u,H по R,Q
        u_next = 0.5*(R_next + Q_next)
        H_next = (R_next - Q_next)**2/(16*g)
        # выводим H
        plt.plot(x, H_next, label=f't = {(n+1)*T/M}', linewidth=2)

plt.plot(x, (2*np.sqrt(g*H_L) - x/T)**2/(9*g))
plt.scatter([-T*np.sqrt(g*H_L)], [H_L])

plt.grid(True)
plt.xlabel('x')
plt.ylabel('H')
plt.legend()