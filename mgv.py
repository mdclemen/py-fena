# 1d multigrid v-cycle for ch5ex9
from pylab import*

def mgv1n(u, b):
    n = len(u)
    N = n-1
    h = 1./N

    if (n==3):
        # solve exactly
        z = asmatrix(zeros([3,1]))
        z[1,0] = -h**2/2.*b[1,0]
    else:
        # one step of GS
        for j in range(1,n-1):
            u[j,0] = (u[j+1,0]+u[j-1,0])/2. - h**2/2.*b[j,0]

        # compute residual
        r = asmatrix(zeros([n,1]))
        r[1:n-1,0] = b[1:n-1,0] + (2.*u[1:n-1,0]-u[0:n-2,0]-u[2:n,0]) / h**2

        # restrict residual
        m = (n+1)//2
        rhat = asmatrix(zeros([m,1]))
        rhat[1:m-1,0] = (r[1:n-3:2,0] + 2*r[2:n-2:2,0] + r[3:n-1:2,0])/4.

        # recur
        uhat = mgv1n(asmatrix(zeros([size(rhat),1])),rhat)

        # interpolate
        ucor = asmatrix(zeros([n,1]))
        ucor[2:n-2:2,0] = uhat[1:m-1,0]
        ucor[1:n-1:2,0] = .5*(uhat[0:m-1,0]+uhat[1:m,0])
        z = u + ucor

        # one more step of GS
        for j in range(1,n-1):
            z[j,0] = (z[j+1,0]+z[j-1,0])/2. - h**2/2.*b[j,0]

    return z

def mgv1r(u, b):
# simple restriction
    n = len(u)
    N = n-1
    h = 1./N

    if (n==3):
        # solve exactly
        z = asmatrix(zeros([3,1]))
        z[1,0] = -h**2/2.*b[1,0]
    else:

        # one step of GS
        for j in range(1,n-1):
            u[j,0] = (u[j+1,0]+u[j-1,0])/2. - h**2/2*b[j,0]

        # compute residual
        r = asmatrix(zeros([n,1]))
        r[1:n-1,0] = b[1:n-1,0] + (2*u[1:n-1,0]-u[0:n-2,0]-u[2:n,0]) / h**2

        # restrict residual
        m = (n+1)//2
        rhat = asmatrix(zeros([m,1]))
        rhat[1:m-1,0] = r[2:n-2:2,0]

        # recur
        uhat = mgv1n(asmatrix(zeros([size(rhat),1])),rhat)

        # interpolate
        ucor = asmatrix(zeros([n,1]))
        ucor[2:n-2:2,0] = uhat[1:m-1,0]
        ucor[1:n-1:2,0] = .5*(uhat[0:m-1,0]+uhat[1:m,0])
        z = u + ucor

        # one more step of GS
        for j in range(1,n-1):
            z[j,0] = (z[j+1,0]+z[j-1,0])/2. - h**2/2.*b[j,0]

    return z

def mgv2d(U, B, gs_iter):
# 2d multigrid v-cycle for ch5ex10
    n, n1 = shape(B)

    assert n == n1, "mgv2d: input not square"

    N = n-1
    h = 2/N

    if (n == 3):
        # final solution
        Z = asmatrix(zeros([n,n]))
        Z[1,1] = -h**2/4.*B[1,1]
    else:
        # k GS iterations
        for niter in range(gs_iter):
            for j in range(1,n-1):
                for i in range(1,n-1):
                    U[i,j] = (U[i+1,j]+U[i-1,j]+U[i,j+1]+U[i,j-1])/4. - h**2/4.*B[i,j]

        # compute residual
        cen = slice(1, n-1, 1)
        dec = slice(0, n-2, 1)
        inc = slice(2, n, 1)
        R = asmatrix(zeros([n,n]))
        R[cen,cen] = B[cen,cen] - 1./h**2*( \
            -4.*U[cen,cen] \
            +U[dec,cen] \
            +U[inc,cen] \
            +U[cen,dec] \
            +U[cen,inc] )

        # restrict residual
        cen_s = slice(2, n-2, 2)
        dec_s = slice(1, n-3, 2)
        inc_s = slice(3, n-1, 2)

        m = (n+1)//2
        Rhat = asmatrix(zeros([m,m]))
        Rhat[1:m-1,1:m-1] = 1./16.*( \
            R[dec_s,dec_s] + \
            R[dec_s,inc_s] + \
            R[inc_s,dec_s] + \
            R[inc_s,inc_s] + \
            2*R[inc_s,cen_s] + \
            2*R[dec_s,cen_s] + \
            2*R[cen_s,inc_s] + \
            2*R[cen_s,dec_s] + \
            4*R[cen_s,cen_s])

        # recur
        Uhat = mgv2d(asmatrix(zeros([m,m])),Rhat,gs_iter)

        # interpolate
        Ucor = asmatrix(zeros([n,n]))
        Ucor[cen_s,cen_s] = Uhat[1:m-1,1:m-1]
        Ucor[1:n-1:2,cen_s] = .5*(Uhat[0:m-1,1:m-1]+Uhat[1:m,1:m-1])
        Ucor[cen_s,1:n-1:2] = .5*(Uhat[1:m-1,0:m-1]+Uhat[1:m-1,1:m])
        Ucor[1:n-1:2,1:n-1:2] = .25*( \
            Uhat[0:m-1,0:m-1] + \
            Uhat[1:m,0:m-1] + \
            Uhat[0:m-1,1:m] + \
            Uhat[1:m,1:m])

        Z = U + Ucor

        # k steps of GS
        for niter in range(gs_iter):
            for j in range(1,n-1):
                for i in range(1,n-1):
                    Z[i,j] = (Z[i+1,j]+Z[i-1,j]+Z[i,j+1]+Z[i,j-1])/4. - h**2/4.*B[i,j]

    return Z
