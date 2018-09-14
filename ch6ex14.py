# Example 6.14
# One-Dimensional boundary value problem
from pylab import*
from tri_diag import*

clf()
x1 = 0.
x2 = 1.
xe = linspace(x1, x2, 1000)
ue = -6.*xe + np.power(xe, 3) + 5.*(1./sin(1.))*sin(xe)

ffunc = lambda x: np.power(x, 3)

A, b, x4 = makeAb(4, ffunc)
u4 = tri_diag(A, b)
x4b = hstack((x1, x4, x2))
u4b = hstack((0., u4, 0.))

A, b, x8 = makeAb(8, ffunc)
u8 = tri_diag(A, b)
x8b = hstack((x1, x8, x2))
u8b = hstack((0., u8, 0.))

A, b, x16 = makeAb(16, ffunc)
u16 = tri_diag(A, b)
x16b = hstack((x1, x16, x2))
u16b = hstack((0., u16, 0.))

plot(x4b, u4b, "o-", mfc = "none", label = "N=4")
plot(x8b, u8b, "d-", mfc = "none", label = "N=8")
plot(x16b, u16b, "k+", mfc = "none", label = "N=16")
plot(xe, ue, "k", label = "exact")
legend(loc = "lower left")
xlabel("x")
title("u(x)")
grid("on")
show()
