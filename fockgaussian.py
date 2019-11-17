"""
fockgaussian
============
Provives a simple function to calculate the Fock matrix elements of Gaussian
unitary using loop hafnians.
"""


import numpy as np
from thewalrus import hafnian
from strawberryfields.decompositions import takagi
import strawberryfields.backends.gaussianbackend.gaussiancircuit as gc


# pylint: disable=invalid-name


def tmsq(state, i, j, r):
    """ Given a gaussiancircuit object it applies a two mode squeezing operator
    by amount r between modes i and j using the decomposition of this operation
    in terms of beamsplitters and (single mode) squeezers.

    Args:
    	state (gaussiancircuit): A gaussiancircuit object
    	i,j (integers): The two modes in which to apply the squeezing operation
    	r (float): Squeezing parameter
    """
    state.beamsplitter(np.pi / 4, 0, i, j)
    state.squeeze(-r, 0, i)
    state.squeeze(r, 0, j)
    state.beamsplitter(-np.pi / 4, 0, i, j)


# pylint: disable=too-many-arguments, too-many-locals
def matelem(l, m, n, U, Up, ls, alpha):
    """ Calculates a Fock matrix element <m|W(alpha,U,ls,Up)|n> of the Gaussian
    unitary W specified by alpha, U, ls, Up.

    Args:
    	l (integer): Number of modes
    	m (list): List of integers specifying the input Fock states
    	n (list): List of integers specifying the output Fock states
    	U (array): Unitary matrix of size l
    	Up (array): Unitary matrix of size l
    	ls (array): Squeezing parameters
    	alpha (array): Complex displacements
    Returns:
    	(complex): Value of the required matrix element
    """
    assert l == len(m)
    assert l == len(n)
    assert U.shape == (l, l)
    assert Up.shape == (l, l)
    assert len(ls) == l
    assert len(alpha) == l

    idl = np.identity(l)
    # Define extended unitaries that are identities in the second half of the mode
    Ue = np.block([[U, 0 * idl], [0 * idl, idl]])
    Uep = np.block([[Up, 0 * idl], [0 * idl, idl]])

    # Define the ts of the squeezing parameters
    # pylint: disable=assignment-from-no-return
    ts = np.arcsinh(np.sqrt(1.0 * np.array(n)))

    # Now we generate the circuit in Fig 4.(b)
    nmodes = 2 * l
    state = gc.GaussianModes(nmodes)
    for i, t in enumerate(ts):
        tmsq(state, i, i + l, -t)

    state.apply_u(Uep)

    for i, lval in enumerate(ls):
        state.squeeze(-lval, 0, i)

    state.apply_u(Ue)

    # Shortcircuited Bloch-Messiah using Takagi
    Mt = state.mmat
    lt, ut = takagi(Mt, 15)

    # Define the lambda tilde and the u tilde
    lt = -0.5 * np.arcsinh(2 * lt)
    ut = ut.conj()

    alphat = np.array(list(alpha) + list(np.zeros_like(alpha)))
    B = ut @ np.diag(np.tanh(lt)) @ ut.T
    zeta = alphat - B @ alphat.conj()

    pref = -0.5 * alphat.conj() @ zeta

    p = m + n

    # Calculating prefactors
    R = 1.0 / np.prod((np.tanh(ts) ** n) / np.cosh(ts))
    prefns = np.sqrt(np.prod(np.array([np.math.factorial(i) for i in p])))
    T = np.exp(pref) / (prefns * np.sqrt(np.prod(np.cosh(lt))))

    # Calculating the multiset S_p
    sp = []
    for k, pval in enumerate(p):
        for i in range(pval):
            sp.append(k)

    # Generate Bp with possibly repeated rows and columns
    Bp = B[:, sp][sp, :]

    # Generate zetap with possibly repeated entries
    zetap = zeta[sp]

    # Calculate Bt
    np.fill_diagonal(Bp, zetap)

    if Bp.shape == (0, 0):
        amp = 1.0
    else:
        amp = hafnian(Bp, loop=True)

    mu = R * T * amp

    return mu
