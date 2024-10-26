import cupy as cp

def norm_KR(A, tol=1e-4, x0=None, delta=0.1, Delta=3, fl=0):
    """
    norm_KR A balancing algorithm for symmetric matrices using CuPy for GPU acceleration.
    
    Args:
        A (cupy.ndarray): Symmetric and nonnegative matrix.
        tol (float): Error tolerance.
        x0 (cupy.ndarray, optional): Initial guess vector.
        delta (float): Lower bound for balancing vectors.
        Delta (float): Upper bound for balancing vectors.
        fl (int): If 1, prints intermediate convergence statistics.
        
    Returns:
        x (cupy.ndarray): The vector x such that diag(x)*A*diag(x) is close to doubly stochastic.
        res (list): Residual error at each outer iteration.
        
    Performs KR normalization. The function is a translation of the 'Matlab' code provided in the 2012 manuscript. 
    Knight PA, Ruiz D. A fast algorithm for matrix balancing. IMA Journal of Numerical Analysis. Oxford University Press; 2012;
    """

    # Ensure A is a CuPy array
    A = cp.asarray(A)  # Convert A to CuPy array if it’s not already
    
    # Initialize parameters and variables
    n = A.shape[0]
    e = cp.ones((n, 1), dtype=A.dtype)  # CuPy array of ones
    res = []
    if x0 is None:
        x0 = e
    else:
        x0 = cp.asarray(x0)  # Convert x0 to CuPy array if it’s not already
    
    g = 0.9
    etamax = 0.1
    eta = etamax
    stop_tol = tol * 0.5
    x = x0
    rt = tol**2
    v = x * (A @ x)
    rk = 1 - v
    rho_km1 = rk.T @ rk
    rout = rho_km1
    rold = rout
    MVP = 0  # Count matrix-vector products
    i = 0  # Outer iteration count

    if fl == 1:
        print("it   in. it   res")

    # Outer iteration
    while rout > rt:
        i += 1
        k = 0
        y = e
        innertol = max([eta**2 * rout, rt])

        # Inner iteration by Conjugate Gradient (CG)
        while rho_km1 > innertol:
            k += 1
            if k == 1:
                Z = rk / v
                p = Z
                rho_km1 = rk.T @ Z
            else:
                beta = rho_km1 / rho_km2
                p = Z + beta * p

            # Update search direction
            w = x * (A @ (x * p)) + v * p
            alpha = rho_km1 / (p.T @ w)
            ap = alpha * p

            # Test distance to boundary of cone
            ynew = y + ap
            if cp.min(ynew) <= delta:
                if delta == 0:
                    break
                ind = cp.where(ap < 0)
                gamma = cp.min((delta - y[ind]) / ap[ind])
                y = y + gamma * ap
                break

            if cp.max(ynew) >= Delta:
                ind = cp.where(ynew > Delta)
                gamma = cp.min((Delta - y[ind]) / ap[ind])
                y = y + gamma * ap
                break

            y = ynew
            rk = rk - alpha * w
            rho_km2 = rho_km1
            Z = rk / v
            rho_km1 = rk.T @ Z

        # Update variables for next iteration
        x = x * y
        v = x * (A @ x)
        rk = 1 - v
        rho_km1 = rk.T @ rk
        rout = rho_km1
        MVP += k + 1

        # Update inner iteration stopping criterion
        rat = rout / rold
        rold = rout
        res_norm = cp.sqrt(rout)
        eta_o = eta
        eta = g * rat
        if g * eta_o**2 > 0.1:
            eta = max([eta, g * eta_o**2])
        eta = max([min([eta, etamax]), stop_tol / res_norm])

        if fl == 1:
            print(f"{i:3d} {k:6d} {res_norm:.3e}")
            res.append(res_norm.item())

    print(f"Matrix-vector products = {MVP}")
    return x
