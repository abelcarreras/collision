import numpy as np

m1 = 1
m2 = 2

v1 = 2
v2 = -2

#coef = [-m1 - m2**2/m1, m2**2/m1*2*v2 - 2*m1*v1, m2*v2**2-m2**2*v2**2/m1 - 2*m1*v1*v2]

mm = m1/m2

coef = [mm**2 + mm, -2*v1*mm**2 - 2*mm*v2, mm**2*v1**2+ 2*mm*v1*v2 - mm*v1**2]

res = np.roots(coef)

print(res)

v1f = res[0]

v2f = mm * (v1 - v1f) + v2
#v2f = np.sqrt(mm * (v1**2 - v1f**2) + v2**2)

print('solution')
print('v1f={:10.5f} v2f={:10.5f}'.format(v1f, v2f))
print('')


pi = m1*v1 + m2*v2
pf = m1*v1f + m2*v2f

print(pi, pf)

ki = m1*v1**2 + m2*v2**2
kf = m1*v1f**2 + m2*v2f**2

print(ki, kf)