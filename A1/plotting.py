import matplotlib.pyplot as plt

strong4 = [1.103426, 0.738230, 0.619121, 0.522384, 0.334046, 0.374576, 0.220918, 0.259324, 0.220100, 0.170415]
threads = [1, 2, 4, 5, 8, 10, 16, 20, 25, 32]

weak1 = [0.267334, 0.185316, 0.109104, 0.053138, 0.159252, 0.125832, 0.078559, 0.063434, 0.051859, 0.058583]
weak2 = [0.532881, 0.369374, 0.286632, 0.262433, 0.228227, 0.261248, 0.160429, 0.130075, 0.104282, 0.125486] 
weak4 = [1.060042, 0.747509, 0.532513, 1.030408, 0.466463, 0.529466, 0.323659, 0.259794, 0.217751, 0.182035] 
weak8 = [2.242319, 1.499842, 1.133611, 1.024943, 0.914709, 0.927197, 0.637910, 0.536501, 0.485546, 0.359222] 

s4 = [e[0]/e for e in strong4]

plt.plot(threads, s4)
plt.show()

