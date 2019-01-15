import numpy as np

_test = False


def collision_two_masses_1d(mass_1, mass_2, vel_1, vel_2):

    mm = mass_1 / mass_2
    coefficients = [mm ** 2 + mm,
                    -2 * vel_1 * mm ** 2 - 2 * mm * vel_2,
                    mm ** 2 * vel_1 ** 2 + 2 * mm * vel_1 * vel_2 - mm * vel_1 ** 2]
    res = np.roots(coefficients)

    if np.abs(res[0] - vel_1) < np.abs(res[1] - vel_1):
        vel_1f = res[1]
    else:
        vel_1f = res[0]

    vel_2f = mm * (vel_1 - vel_1f) + vel_2

    print('v1={:10.5f} v2={:10.5f}'.format(vel_1, vel_2))
    print('v1f={:10.5f} v2f={:10.5f}'.format(vel_1f, vel_2f))

    return vel_1f, vel_2f


def collision_two_masses_nd(mass_1, mass_2, vel_1, vel_2):

    vel_1f = []
    vel_2f = []
    for v1, v2 in zip(vel_1, vel_2):
        v1f, v2f = collision_two_masses_1d(mass_1, mass_2, v1, v2)
        vel_1f.append(v1f)
        vel_2f.append(v2f)

    return vel_1f, vel_2f

m1 = 1
m2 = 2

v1 = 3
v2 = -2

v1f, v2f = collision_two_masses_1d(m1, m2, v1, v2)

print('solution')
print('v1f={:10.5f} v2f={:10.5f}'.format(v1f, v2f))
print('')

if _test:
    print('test')
    pi = m1*v1 + m2*v2
    pf = m1*v1f + m2*v2f

    print('momentum', pi, pf)

    ki = m1*v1**2 + m2*v2**2
    kf = m1*v1f**2 + m2*v2f**2

    print('kinetic', ki, kf)
