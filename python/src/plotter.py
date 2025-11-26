"""
Sistema base para graficado científico de datos Ising.
"""

import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from typing import Tuple, Optional

from config import (
    FIGURE_DPI, FIGURE_DPI_HIGH, DEFAULT_FIGSIZE, LABELS, 
    MARKER_SIZE, LINE_WIDTH, MARKER_STYLE, CMAP_SPIN,
    OUTPUT_DIR
)

import os


class IsingPlotter:
    """Base para sistemas de graficado con formato científico"""
    
    def __init__(self):
        """Inicializa directorios de salida"""
        self.output_dir = OUTPUT_DIR
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Crear subdirectorios
        self.subdirs = {
            'para': self.output_dir / 'paramagnetismo',
            'hyst': self.output_dir / 'ferromagnetismo',
            'trans': self.output_dir / 'transicion',
            'snap': self.output_dir / 'snapshots',
            'comp': self.output_dir / 'comparacion_z',
            'relax': self.output_dir / 'relajacion',
            'analysis': self.output_dir / 'analisis',
            'tc': self.output_dir / 'temperatura_critica',
            'm_vs_t': self.output_dir / 'm_vs_T',
        }
        
        for subdir in self.subdirs.values():
            subdir.mkdir(parents=True, exist_ok=True)
        
        # Configurar matplotlib para estilo científico
        self._configure_matplotlib()
    
    @staticmethod
    def _configure_matplotlib() -> None:
        """Configura matplotlib para salida científica de calidad"""
        try:
            plt.style.use('seaborn-v0_8-darkgrid')
        except:
            # Si seaborn no está disponible, usar default
            pass
        
        plt.rcParams.update({
            'figure.dpi': FIGURE_DPI,
            'savefig.dpi': FIGURE_DPI,
            'font.size': 10,
            'font.family': 'sans-serif',
            'axes.labelsize': 11,
            'axes.titlesize': 12,
            'xtick.labelsize': 9,
            'ytick.labelsize': 9,
            'legend.fontsize': 10,
            'lines.linewidth': LINE_WIDTH,
            'lines.markersize': MARKER_SIZE,
            'grid.alpha': 0.3,
        })
    
    def save_figure(self, output_dir: str, filename: str, 
                   dpi: Optional[int] = None) -> Path:
        """Guarda figura con metadata científica"""
        dpi = dpi or FIGURE_DPI
        
        if output_dir not in self.subdirs:
            raise ValueError(f"Directorio '{output_dir}' no reconocido")
        
        output_path = self.subdirs[output_dir] / filename
        
        plt.savefig(output_path, dpi=dpi, bbox_inches='tight', 
                   facecolor='white', edgecolor='none')
        print(f"✓ {output_path}")
        
        return output_path
    
    @staticmethod
    def create_figure(figsize: Tuple[int, int] = DEFAULT_FIGSIZE,
                     title: Optional[str] = None) -> Tuple:
        """Crea figura con formato estándar"""
        fig, ax = plt.subplots(figsize=figsize)
        if title:
            fig.suptitle(title, fontsize=13, fontweight='bold')
        return fig, ax
    
    @staticmethod
    def create_grid_figure(rows: int, cols: int,
                          figsize: Tuple[int, int] = (12, 10),
                          title: Optional[str] = None) -> Tuple:
        """Crea grilla de subplots"""
        fig, axes = plt.subplots(rows, cols, figsize=figsize)
        if title:
            fig.suptitle(title, fontsize=13, fontweight='bold')
        
        # Asegurar que siempre es un array
        if rows * cols == 1:
            axes = np.array([axes])
        
        return fig, axes.flatten()
    
    @staticmethod
    def format_scientific_label(label_key: str, **kwargs) -> str:
        """Obtiene etiqueta formateada con unidades física (k_B = 1)"""
        base_label = LABELS.get(label_key, label_key)
        
        # Substituciones adicionales si es necesario
        for key, value in kwargs.items():
            base_label = base_label.replace(f'{{{key}}}', str(value))
        
        return base_label
    
    @staticmethod
    def add_scientific_text(ax, x: float, y: float, text: str,
                           **kwargs) -> None:
        """Añade texto con formato científico"""
        ax.text(x, y, text, transform=ax.transAxes,
               fontsize=10, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
               **kwargs)
    
    @staticmethod
    def plot_with_uncertainty(ax, x: np.ndarray, y: np.ndarray,
                             dy: Optional[np.ndarray] = None,
                             label: Optional[str] = None,
                             **plot_kwargs) -> None:
        """Grafica con barras de error"""
        ax.plot(x, y, MARKER_STYLE, label=label, **plot_kwargs)
        if dy is not None:
            ax.fill_between(x, y - dy, y + dy, alpha=0.2)
    
    @staticmethod
    def add_critical_temperature_line(ax, Tc: float, 
                                     label: Optional[str] = None,
                                     color: str = 'red',
                                     linestyle: str = '--',
                                     alpha: float = 0.7) -> None:
        """Añade línea vertical de temperatura crítica"""
        if label is None:
            label = f'$T_c$ = {Tc:.3f}'
        ax.axvline(Tc, color=color, linestyle=linestyle, 
                  alpha=alpha, linewidth=2, label=label)