import strawberryfields.backends.gaussianbackend.gaussiancircuit as gc
from strawberryfields.decompositions import takagi
import hafnian as haf
import numpy as np


def tmsq(state, i, j, r):
    """ Given a gaussiancircuit object it applies a two mode squeezing operator
    by amount r between modes i and j using the decomposition of this operation
    in terms of beamsplitters and (single mode) squeezers.

    Args: 

    state (gaussiancircuit): A gaussiancircuit object
    i,j (integers): The two modes in which to apply the squeezing operation
    r (real): Squeezing parameter
    """
    state.beamsplitter(np.pi/4, 0, i, j)
    state.squeeze(-r, 0, i)
    state.squeeze(r, 0, j)
    state.beamsplitter(-np.pi/4, 0, i, j)


def matelem(l, m, n, U, Up, ls, alpha):
    """ Calculates a Fock matrix element <m|W(alpha,U,ls,Up|n> of the Gaussian
    unitary W specified by alpha, U, ls, Up.

    Args:

    l (integer): Number of modes
    m,n : Lists of integers of length l
    U, Up: Unitary matrices (square numpy arrays) of length l
    ls: Array of real floats specifying the squeezing parameters
    alpha: Array of complex floats specifying the displacements
    """
    assert l == len(m)
    assert l == len(n)
    assert U.shape == (l, l)
    assert Up.shape == (l, l)
    assert len(ls) == l
    assert len(alpha) == l

    idl = np.identity(l)
    # Define extended unitaries that are identities in the second half of the mode
    Ue = np.block([[U, 0*idl], [0*idl, idl]])
    Uep = np.block([[Up, 0*idl], [0*idl, idl]])

    # Define the ts of the squeezing parameters
    ts = np.arcsinh(np.sqrt(1.0*np.array(n)))

    # Now we generate the circuit in Fig 4.(b)
    nmodes = 2*l
    state = gc.GaussianModes(nmodes, hbar=2)
    for i, t in enumerate(ts):
        tmsq(state, i, i+l, -t)

    state.apply_u(Uep.conj())

    for i, l in enumerate(ls):
        state.squeeze(-l, 0, i)

    state.apply_u(Ue.conj())

    # Shortcircuited Bloch-Messiah using Takagi
    Mt = state.mmat
    Nt = state.mmat
    lt, ut = takagi(Mt, 15)

    # Define the lambda tilde and the u tilde
    lt = -0.5*np.arcsinh(2*lt)
    ut = ut.conj()

    alphat = np.array(list(alpha)+list(np.zeros_like(alpha)))
    B = ut @ np.diag(np.tanh(lt)) @ ut.T
    zeta = alphat - B @ alphat.conj()

    pref = -0.5*alphat.conj()@zeta

    p = m+n

    # Calculating prefactors
    R = 1.0/np.prod((np.tanh(ts)**n)/np.cosh(ts))
    prefns = np.sqrt(np.prod(np.array([np.math.factorial(i) for i in p])))
    T = np.exp(pref)/(prefns*np.sqrt(np.prod(np.cosh(lt))))

    # Calculating the multiser S_p
    sp = []
    for k, pval in enumerate(p):
        for i in range(pval):
            sp.append(k)

    # Generate Bp with possibly repeated rows and columns
    Bp = B[:, sp][sp, :]

    # Generate zetap with possibly repeated entries
    zetap = zeta[sp]

    # Calculate Bt
    Bt = Bp-np.diag(np.diag(Bp))+np.diag(zetap)

    if sum(p) % 2 != 0:
        Bt = np.pad(Bt, pad_width=((0, 1), (0, 1)), mode='constant')
        Bt[-1, -1] = 1.0

    if Bt.shape == (0, 0):
        amp = 1.0
    else:
        amp = haf.haf_complex(Bt, loop=True)

    mu = R*T*amp

    return mu
