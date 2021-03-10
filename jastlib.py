import numpy as np
import matplotlib.pyplot as plt

def get_sphere_distribution(n, dmin, Ls, maxiter=1e4, allow_wall=True):
    """Get random points in a box with given dimensions and minimum separation.

    Parameters:

    - n: number of points
    - dmin: minimum distance
    - Ls: dimensions of box, shape (3,) array
    - maxiter: maximum number of iterations.
    - allow_wall: whether to allow points on wall;
       (if False: points need to keep distance dmin/2 from the walls.)

    Return:

    - ps: array (n, 3) of point positions,
      with 0 <= ps[:, i] < Ls[i]
    - n_iter: number of iterations
    - dratio: average nearest-neighbor distance, divided by dmin.

    Note: with a fill density (sphere volume divided by box volume) above about
    0.53, it takes very long. (Random close-packed spheres have a fill density
    of 0.64).

    Author: Han-Kwang Nienhuys (2020)
    Copying: BSD, GPL, LGPL, CC-BY, CC-BY-SA
    See Stackoverflow: https://stackoverflow.com/a/62895898/6228891
    """
    Ls = np.array(Ls).reshape(3)
    if not allow_wall:
        Ls -= dmin

    # filling factor; 0.64 is for random close-packed spheres
    # This is an estimate because close packing is complicated near the walls.
    # It doesn't work well for small L/dmin ratios.
    sphere_vol = np.pi/6*dmin**3
    box_vol = np.prod(Ls + 0.5*dmin)
    fill_dens = n*sphere_vol/box_vol
    if fill_dens > 0.64:
        msg = f'Too many to fit in the volume, density {fill_dens:.3g}>0.64'
        raise ValueError(msg)

    # initial try
    ps = np.random.uniform(size=(n, 3)) * Ls

    # distance-squared matrix (diagonal is self-distance, don't count)
    dsq = ((ps - ps.reshape(n, 1, 3))**2).sum(axis=2)
    dsq[np.arange(n), np.arange(n)] = np.infty

    for iter_no in range(int(maxiter)):
        # find points that have too close neighbors
        close_counts = np.sum(dsq < dmin**2, axis=1)  # shape (n,)
        n_close = np.count_nonzero(close_counts)
        if n_close == 0:
            break

        # Move the one with the largest number of too-close neighbors
        imv = np.argmax(close_counts)

        # new positions
        newp = np.random.uniform(size=3)*Ls
        ps[imv]= newp

        # update distance matrix
        new_dsq_row = ((ps - newp.reshape(1, 3))**2).sum(axis=-1)
        dsq[imv, :] = dsq[:, imv] = new_dsq_row
        dsq[imv, imv] = np.inf
    else:
        raise RuntimeError(f'Failed after {iter_no+1} iterations.')

    if not allow_wall:
        ps += dmin/2

    dratio = (np.sqrt(dsq.min(axis=1))/dmin).mean()
    return ps, iter_no+1, dratio

def generateElsAroundPoints(n,LS,dmin):
    """
    Parameters:

    - n: number of points
    - LS: list of position of all atoms
    - dmin: minimum intra block distance
    - shift: inter block distance

    Return:

    - r: array (n, 3) of point positions,

    """

    xs = None
    for Ls in LS:

        # Get list of random points around Ls
        distrib,a,b = get_sphere_distribution(n,dmin,Ls)

        if xs is None:
            xs = distrib[:,0]
            ys = distrib[:,1]
            zs = distrib[:,2]
        else:
            xs = np.concatenate((xs,distrib[:,0]))
            ys = np.concatenate((ys,distrib[:,1]))
            zs = np.concatenate((zs,distrib[:,2]))

    return((np.array((xs,ys,zs))).T)

def getCoefList(Nord,Natom):
    assert(Nord < 11)
    dict = {
        0 : lambda x,y:x-y-2,
        1 : lambda x,y:x-y,
        2 : lambda x,y:x-y,
        3 : lambda x,y:x-y,
        4 : lambda x,y:x-y,
        5 : lambda x,y:x-y,
        6 : lambda x,y:x-y,
        7 : lambda x,y:x-y,
        8 : lambda x,y:x-y,
        9 : lambda x,y:x-y,
        10 : lambda x,y:x-y,
        11 : lambda x,y:x-y,
    }
    count = 0
    for p in range(2,Nord+1):
        for k in range(p-1,-1,-1):
            lmax = dict[k](p,k)
            for l in range(lmax,-1,-1):
                if (p-k-l) & 1 is 0:
                    count += 1
    coeflista = np.random.rand(Nord+1,Natom)
    coeflistb = np.random.rand(Nord+1)
    coeflistc = np.random.rand(count,Natom)
    return (coeflista.reshape((Nord+1)*Natom),coeflistb,coeflistc.reshape(count*Natom))
    #return (coeflista,coeflistb,coeflistc)

def get_sphere_distribution(n, dmin, Ls, maxiter=1e4, allow_wall=True):
    """Get random points in a box with given dimensions and minimum separation.

    Parameters:

    - n: number of points
    - dmin: minimum distance
    - Ls: dimensions of box, shape (3,) array
    - maxiter: maximum number of iterations.
    - allow_wall: whether to allow points on wall;
       (if False: points need to keep distance dmin/2 from the walls.)

    Return:

    - ps: array (n, 3) of point positions,
      with 0 <= ps[:, i] < Ls[i]
    - n_iter: number of iterations
    - dratio: average nearest-neighbor distance, divided by dmin.

    Note: with a fill density (sphere volume divided by box volume) above about
    0.53, it takes very long. (Random close-packed spheres have a fill density
    of 0.64).

    Author: Han-Kwang Nienhuys (2020)
    Copying: BSD, GPL, LGPL, CC-BY, CC-BY-SA
    See Stackoverflow: https://stackoverflow.com/a/62895898/6228891
    """
    Ls = np.array(Ls).reshape(3)
    if not allow_wall:
        Ls -= dmin

    # filling factor; 0.64 is for random close-packed spheres
    # This is an estimate because close packing is complicated near the walls.
    # It doesn't work well for small L/dmin ratios.
    sphere_vol = np.pi/6*dmin**3
    box_vol = np.prod(Ls + 0.5*dmin)
    fill_dens = n*sphere_vol/box_vol
    if fill_dens > 0.64:
        msg = f'Too many to fit in the volume, density {fill_dens:.3g}>0.64'
        raise ValueError(msg)

    # initial try
    ps = np.random.uniform(size=(n, 3)) * Ls

    # distance-squared matrix (diagonal is self-distance, don't count)
    dsq = ((ps - ps.reshape(n, 1, 3))**2).sum(axis=2)
    dsq[np.arange(n), np.arange(n)] = np.infty

    for iter_no in range(int(maxiter)):
        # find points that have too close neighbors
        close_counts = np.sum(dsq < dmin**2, axis=1)  # shape (n,)
        n_close = np.count_nonzero(close_counts)
        if n_close == 0:
            break

        # Move the one with the largest number of too-close neighbors
        imv = np.argmax(close_counts)

        # new positions
        newp = np.random.uniform(size=3)*Ls
        ps[imv]= newp

        # update distance matrix
        new_dsq_row = ((ps - newp.reshape(1, 3))**2).sum(axis=-1)
        dsq[imv, :] = dsq[:, imv] = new_dsq_row
        dsq[imv, imv] = np.inf
    else:
        raise RuntimeError(f'Failed after {iter_no+1} iterations.')

    if not allow_wall:
        ps += dmin/2

    dratio = (np.sqrt(dsq.min(axis=1))/dmin).mean()
    return ps, iter_no+1, dratio


def scalingee(r,kappa=1.0):
    return (numpy.ones_like(r) - numpy.exp(-kappa*r))/kappa

def scalingen(r,kappa=1.0):
    return numpy.exp(-kappa*r)

if False:
  Nord = 5
  L1 = 2.0
  n = 2 # number of points
  dmin = 0.1 # min dist between points
  Ls = np.array([L1,L1,L1]) # lengths of the box
  shift = -10.0
  kappa = 2.0
  filename_atom = str(n) + "_geometry.txt"
  filename_elec = str(n)
  filename_coeffs = str(n) + "_jast_coeffs.txt"
  (coeffsa, coeffsb, coeffsc) = getCoefList(Nord,n)
  coeffsall = np.concatenate((coeffsa,coeffsb,coeffsc))
  print(coeffsa.shape,coeffsb.shape,coeffsc.shape)

  atomList,_,_ = get_sphere_distribution(n, dmin, Ls, maxiter=1e4, allow_wall=True)
  #print(atomList)

  L1 = 5.0
  n = 5 # number of points
  dmin = 0.1 # min dist between points
  Ls = np.array([L1,L1,L1]) # lengths of the box
  shift = -10.0
  kappa = 2.0
  filename_elec = filename_elec + "_" + str(n) + "_elec_coord.txt"

  #rlist = generateBlockRandomPointsAtShftApart(n,L1,dmin,shift)
  rlist = generateElsAroundPoints(n,atomList,dmin)


  # Save file
  np.savetxt(filename_elec,rlist)
  np.savetxt(filename_atom,atomList)
  np.savetxt(filename_coeffs,coeffsall)

  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')

  xs = rlist.T[0]
  ys = rlist.T[1]
  zs = rlist.T[2]
  ax.scatter(xs, ys, zs, marker='o')

  plt.show()

  rijScaled = np.array([[(lambda xval, yval: np.linalg.norm(xval-yval))(xval,yval) for yval in rlist] for xval in rlist])

  plt.imshow(rijScaled)
  plt.colorbar()
  plt.show()
