"""
Paquete src - Módulos principales del análisis Ising
"""

from .config import *
from .reader import IsingDataReader
from .plotter import IsingPlotter
from .utils import (
    format_q_str,
    try_multiple_filenames,
    read_exact,
    ensure_dir,
    tanh_model,
    fit_tanh,
    estimate_critical_temperature,
    calculate_critical_exponent,
    smooth_array
)

__all__ = [
    'IsingDataReader',
    'IsingPlotter',
    'format_q_str',
    'try_multiple_filenames',
    'read_exact',
    'ensure_dir',
    'tanh_model',
    'fit_tanh',
    'estimate_critical_temperature',
    'calculate_critical_exponent',
    'smooth_array',
]