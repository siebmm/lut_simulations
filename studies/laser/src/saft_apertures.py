"""Independent dipole launch and traction-free reception in an isotropic halfspace.

Local y points into the solid; convention exp[i omega(t-p*x-q*y)].
"""
import numpy as np
ANGLES=np.linspace(0,80,1601)


def polarization(p, q, mode, cl, cs):
    c = cl if mode == "L" else cs
    return np.stack((c*p, c*q), axis=-1) if mode == "L" else np.stack((-c*q, c*p), axis=-1)


def traction(p, q, mode, cl, cs):
    e = polarization(p, q, mode, cl, cs)
    mu, lam = cs**2, cl**2 - 2*cs**2  # common density cancels
    return np.stack((mu*(q*e[..., 0]+p*e[..., 1]),
                     lam*(p*e[..., 0]+q*e[..., 1])+2*mu*q*e[..., 1]), axis=-1)


def boundary(p, cl, cs):
    p = np.asarray(p)
    ql = np.sqrt(1/cl**2-p*p+0j).conj()
    qs = np.sqrt(1/cs**2-p*p+0j).conj()
    matrix = np.stack((traction(p, ql, "L", cl, cs), traction(p, qs, "S", cl, cs)), axis=-1)
    return matrix, ql, qs


def reflection(theta, mode, cl, cs):
    """Modal reflected amplitudes and total surface displacement/velocity vector."""
    p = np.sin(theta)/(cl if mode == "L" else cs)
    matrix, ql, qs = boundary(p, cl, cs)
    qi = ql if mode == "L" else qs
    incident = traction(p, -qi, mode, cl, cs)
    reflected = np.linalg.solve(matrix, -incident[..., None])[..., 0]
    surface = polarization(p, -qi, mode, cl, cs)
    surface += reflected[..., 0, None]*polarization(p, ql, "L", cl, cs)
    surface += reflected[..., 1, None]*polarization(p, qs, "S", cl, cs)
    return reflected, surface


def directivity(band, cl, cs, width_m=.0005, pulse_fwhm_s=6.5e-9):
    """Band-RMS source amplitude, including actual FIR, dipole and integrated pulse.

    Angular stationarity contributes cos(theta)/sqrt(c) to 2D far-field velocity.
    Spatial dipole and integrated-Gaussian factors leave omega*p^2 in power.
    Constants common to an individual mode/band cancel on normalization.
    """
    from saft_imaging import kernel, DT
    theta = np.deg2rad(ANGLES)
    h = kernel(*band[1:])
    nfft = 16384
    f = np.fft.fftfreq(nfft, DT)
    positive = f > 0
    omega = 2*np.pi*f[positive]
    response = abs(np.fft.fft(h, nfft)[positive])**2
    sigma = pulse_fwhm_s/(2*np.sqrt(2*np.log(2)))
    result = {}
    for mode, c, index in [("L", cl, 0), ("S", cs, 1)]:
        p = np.sin(theta)/c
        matrix, _, _ = boundary(p, cl, cs)
        rhs = np.zeros((len(p), 2), complex); rhs[:, 0] = 1
        amplitude = np.linalg.solve(matrix, rhs[..., None])[..., 0][:, index]
        integral = np.trapezoid(omega[None, :]*response[None, :]*np.exp(
            -omega[None, :]**2*(p[:, None]**2*width_m**2/4+sigma**2)), omega, axis=-1)
        launch = np.cos(theta)*abs(amplitude*p)*np.sqrt(integral)
        _, surface = reflection(theta, mode, cl, cs)
        result[mode] = {"source": (launch/launch.max()).tolist()}
        for component, column in [("vx", 0), ("vy", 1)]:
            receive = abs(surface[:, column])
            result[mode][component] = (receive/receive.max()).tolist()
    return result


def angular_gain(angle_radians, curve):
    return np.interp(np.abs(np.rad2deg(angle_radians)), ANGLES, curve, right=0)
