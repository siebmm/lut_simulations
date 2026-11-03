"""Independent finite-support analytic filtering for the weld physics check."""
import numpy as np
DT=1e-8
FS=1/DT
TAPS=657
HALF=TAPS//2
BETA=.1102*(55-8.7)


def kernel(low_MHz, high_MHz):
    """Finite complex FIR: real part even, imaginary part odd; positive-band gain 2.

    Nominal edges are approximately half-amplitude. The real part is a
    symmetric bandpass; its quadrature partner shares its exact support.
    """
    if not 0 < low_MHz < high_MHz < FS / 2e6:
        raise ValueError("Invalid band edges")
    n = np.arange(-HALF, HALF + 1)
    width, centre = (high_MHz - low_MHz) * 1e6, (high_MHz + low_MHz) * .5e6
    h = 2 * width / FS * np.sinc(width * n / FS) * np.exp(2j * np.pi * centre * n / FS)
    h *= np.kaiser(TAPS, BETA)
    h *= 2 / np.sum(h * np.exp(-2j * np.pi * centre * n / FS)).real
    return h


def filter_record(values, h):
    """Linear convolution, centred without delay; padded edges are later excluded."""
    values = np.asarray(values, dtype=np.float64)
    h = np.asarray(h, dtype=np.complex128)
    if values.ndim != 2 or h.ndim != 1 or len(h) % 2 != 1 or not np.isfinite(values).all():
        raise ValueError("Expected finite receiver-by-time data and an odd FIR")
    length = 1 << (values.shape[-1] + len(h) - 2).bit_length()
    result = np.fft.ifft(np.fft.fft(values, n=length, axis=-1) * np.fft.fft(h, n=length), axis=-1)
    half = len(h) // 2
    return result[:, half:half + values.shape[-1]]
