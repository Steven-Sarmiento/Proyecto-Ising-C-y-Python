#!/usr/bin/env python3
"""
Script principal para generar gráficas de análisis Ising.
Versión mejorada con formato científico, unidades explícitas y mejor organización.

Estructura: Una figura por RED, con subplots para cada Q

Ejecutar desde: PROYECTO_ISING_HIBRIDO/python/
Comando: python main.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / 'src'))

from src.config import (
    PARA_DATA_DIR, HYST_DATA_DIR, TRANS_DATA_DIR,
    SNAP_DATA_DIR, OUTPUT_DIR, Q_VALUES, LATTICES, TC_DATA_DIR, RELAX_DATA_DIR
)
from src.reader import IsingDataReader
from src.utils import format_q_str, try_multiple_filenames, calculate_blocking_error, calculate_hysteresis_area

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit



def print_header(title: str) -> None:
    """Imprime encabezado formateado"""
    print(f"\n{'='*70}")
    print(f"📊 {title}")
    print(f"{'='*70}\n")


def verify_data_exists() -> bool:
    """Verifica que los datos existen antes de continuar"""
    print("\n" + "="*70)
    print("🔍 VERIFICANDO DATOS")
    print("="*70 + "\n")
    
    dirs_to_check = {
        "Paramagnetismo": PARA_DATA_DIR,
        "Ferromagnetismo": HYST_DATA_DIR,
        "Transiciones": TRANS_DATA_DIR,
        "Snapshots": SNAP_DATA_DIR,
    }
    
    all_exist = True
    for name, path in dirs_to_check.items():
        if path.exists():
            bin_files = list(path.glob('*.bin'))
            status = "✓" if bin_files else "⚠ (vacío)"
            print(f"{status} {name:25s} → {len(bin_files)} archivos")
        else:
            print(f"✗ {name:25s} → NO EXISTE")
            all_exist = False
    
    print("\n" + "="*70)
    if not all_exist:
        print("❌ Faltan carpetas de datos\n")
        return False
    
    print("✅ Todas las carpetas de datos existen\n")
    return True


def plot_paramagnetism() -> None:
    """Genera gráficas de magnetización vs campo con BARRAS DE ERROR"""
    print_header("1. PARAMAGNETISMO - m(H)")
    
    reader = IsingDataReader()
    
    for lattice in LATTICES:
        print(f"  Generando para {lattice.name_str.upper()}...")
        
        fig, axes = plt.subplots(1, 3, figsize=(16, 5))
        fig.suptitle(
            f'Magnetización vs Campo - {lattice.name_str.upper()} ($z = {lattice.coordination_number}$)\n'
            r'Sistema natural: $k_B = 1$, $J = 1$',
            fontsize=14, fontweight='bold', y=1.02
        )
        
        # Intentar leer datos RAW primero
        raw_filename = PARA_DATA_DIR / f'{lattice.name_str}_q50.bin' # Nombre ejemplo para buscar
        # Nota: El nombre del archivo depende de q, lo buscamos dentro del loop
        
        for idx, q in enumerate(Q_VALUES):
            q_str = format_q_str(q)
            # Buscamos el archivo específico para este q
            filename = PARA_DATA_DIR / f'{lattice.name_str}_{q_str}.bin'
            
            # A) Intentar leer RAW
            data_raw = reader.read_paramagnetism_raw(str(filename))
            
            # Variables para plotear
            H_axis, m_means, m_errs, T_labels = [], [], [], []
            has_error = False
            
            ax = axes[idx]

            if data_raw is not None:
                # PROCESAMIENTO RAW (Con error)
                print(f"    [q={q}] Datos RAW detectados. Calculando errores...")
                T_arr = data_raw['T']
                H_arr = data_raw['H']
                m_matrix = data_raw['m_raw'] # [T][H][Samples]
                
                H_axis = H_arr
                T_labels = T_arr
                
                # Procesar cada temperatura
                for t_i in range(len(T_arr)):
                    m_per_H = []
                    err_per_H = []
                    for h_i in range(len(H_arr)):
                        series = m_matrix[t_i, h_i, :] # Serie de tiempo
                        m_mean = np.mean(series)
                        # Calcular error con método de bloques
                        err = calculate_blocking_error(series)
                        
                        m_per_H.append(m_mean)
                        err_per_H.append(err)
                    
                    m_means.append(m_per_H)
                    m_errs.append(err_per_H)
                
                has_error = True
                
            else:
                # PROCESAMIENTO ANTIGUO (Sin error)
                data = reader.read_paramagnetism(str(filename))
                if data is None:
                    ax.text(0.5, 0.5, 'No disponible', ha='center', transform=ax.transAxes)
                    continue
                
                H_axis = data['H']
                T_labels = data['T']
                m_means = data['m'] # Ya viene promediado
                m_errs = [None]*len(T_labels)
                has_error = False

            # --- GRAFICAR ---
            cmap = plt.get_cmap('viridis')
            colors = cmap(np.linspace(0, 1, len(T_labels)))
            
            for i, T in enumerate(T_labels):
                label_txt = f'$T = {T:.2f}$'
                
                if has_error:
                    ax.errorbar(H_axis, m_means[i], yerr=m_errs[i], fmt='o-', 
                               label=label_txt, color=colors[i], 
                               ecolor='black', capsize=2, elinewidth=1, alpha=0.8, markersize=4)
                else:
                    ax.plot(H_axis, m_means[i], 'o-', label=label_txt, 
                           color=colors[i], alpha=0.8)
            
            ax.set_xlabel(r'Campo Magnético Efectivo $H$', fontsize=11)
            ax.set_ylabel(r'Magnetización $m$', fontsize=11)
            ax.set_title(f'$q = {q}$', fontsize=12)
            ax.legend(fontsize=9)
            ax.grid(True, alpha=0.3)
            
        plt.tight_layout()
        out_path = OUTPUT_DIR / 'paramagnetismo' / f'01a_m_vs_H_{lattice.name_str}.png'
        out_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(out_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"    ✓ Guardado: {out_path.name}")


def plot_susceptibility() -> None:
    """Genera gráficas de susceptibilidad con BARRAS DE ERROR (Estimación Estadística)"""
    print_header("2. SUSCEPTIBILIDAD MAGNÉTICA - χ(H)")
    
    reader = IsingDataReader()
    
    for lattice in LATTICES:
        print(f"  Generando para {lattice.name_str.upper()}...")
        
        fig, axes = plt.subplots(1, 3, figsize=(16, 5))
        fig.suptitle(
            f'Susceptibilidad Magnética - {lattice.name_str.upper()} ($z = {lattice.coordination_number}$)\n'
            r'$Sistema natural: k_B = 1$',
            fontsize=14, fontweight='bold', y=1.02
        )
        
        for idx, q in enumerate(Q_VALUES):
            q_str = format_q_str(q)
            filename = PARA_DATA_DIR / f'{lattice.name_str}_{q_str}.bin'
            
            # A) Intentar leer RAW
            data_raw = reader.read_paramagnetism_raw(str(filename))
            
            ax = axes[idx]
            H_axis, chi_vals, chi_errs, T_labels = [], [], [], []
            has_error = False
            
            if data_raw is not None:
                print(f"    [q={q}] Datos RAW para Susceptibilidad...")
                T_arr = data_raw['T']
                H_arr = data_raw['H']
                m_matrix = data_raw['m_raw']
                n_samples = data_raw['n_samples']
                
                H_axis = H_arr
                T_labels = T_arr
                
                for t_i, T in enumerate(T_arr):
                    chi_per_H = []
                    chi_err_per_H = []
                    
                    for h_i in range(len(H_arr)):
                        series = m_matrix[t_i, h_i, :]
                        
                        # 1. Calcular Susceptibilidad: Var(m) / T
                        # Nota: Como no tenemos N exacto aquí, calculamos la "Susceptibilidad por spin"
                        # que es proporcional a la varianza de la serie intensiva.
                        variance = np.var(series, ddof=1)
                        chi = variance / T if T > 0 else 0
                        
                        # 2. Estimar Error de Chi
                        # Usamos el truco: error_relativo(chi) approx sqrt(2 / N_eff)
                        # N_eff lo sacamos comparando error blocking vs error simple
                        
                        std_naive = np.std(series, ddof=1) / np.sqrt(n_samples)
                        std_blocking = calculate_blocking_error(series)
                        
                        if std_naive > 0 and std_blocking > 0:
                            # Factor de correlación tau
                            # sigma_block^2 = sigma_naive^2 * (2*tau)
                            tau_factor = (std_blocking / std_naive)**2
                            N_eff = n_samples / tau_factor
                            
                            # Error de la varianza para datos gaussianos
                            relative_error = np.sqrt(2.0 / N_eff)
                            chi_err = chi * relative_error
                        else:
                            chi_err = 0.0
                            
                        chi_per_H.append(chi)
                        chi_err_per_H.append(chi_err)
                        
                    chi_vals.append(chi_per_H)
                    chi_errs.append(chi_err_per_H)
                has_error = True
                
            else:
                # PROCESAMIENTO ANTIGUO
                data = reader.read_paramagnetism(str(filename))
                if data is None:
                    ax.text(0.5, 0.5, 'No disponible', ha='center', transform=ax.transAxes)
                    continue
                H_axis = data['H']
                T_labels = data['T']
                chi_vals = data['chi']
                chi_errs = [None]*len(T_labels)
                has_error = False

            # --- GRAFICAR ---
            colors = plt.get_cmap('plasma')(np.linspace(0, 1, len(T_labels)))
            for i, T in enumerate(T_labels):
                if has_error:
                    ax.errorbar(H_axis, chi_vals[i], yerr=chi_errs[i], fmt='s-', 
                               label=f'$T = {T:.2f}$', color=colors[i], 
                               ecolor='black', capsize=2, alpha=0.8, markersize=4)
                else:
                    ax.plot(H_axis, chi_vals[i], 's-', label=f'$T = {T:.2f}$', 
                           color=colors[i], alpha=0.8)
            
            ax.set_xlabel(r'Campo Magnético Efectivo $H$', fontsize=11)
            ax.set_ylabel(r'Susceptibilidad $\chi$', fontsize=11)
            ax.set_title(f'$q = {q}$', fontsize=12)
            ax.legend(fontsize=9)
            ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        out_path = OUTPUT_DIR / 'paramagnetismo' / f'02_chi_vs_H_{lattice.name_str}.png'
        plt.savefig(out_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"    ✓ Guardado: {out_path.name}")

def plot_hysteresis() -> None:
    """Genera gráficas de histéresis magnética con múltiples temperaturas"""
    print_header("3. FERROMAGNETISMO - HISTÉRESIS m(H)")
    
    reader = IsingDataReader()
    
    # Importar las nuevas constantes
    from src.config import T_VALUES_HYSTERESIS
    from src.utils import build_filename_candidates
    
    # Diccionario para almacenar áreas (para tabla resumen)
    areas_table = {}
    
    for lattice in LATTICES:
        print(f"  Generando para {lattice.name_str.upper()}...")
        
        # Crear figura: 3 filas (q) × 2 columnas (T)
        fig, axes = plt.subplots(3, 2, figsize=(14, 15))
        fig.suptitle(
            f'Histéresis Magnética - {lattice.name_str.upper()} ($z = {lattice.coordination_number}$)\n'
            r'Sistema natural: $k_B = 1$, $J = 1$',
            fontsize=14, fontweight='bold', y=0.995
        )
        
        areas_table[lattice.name_str] = {}
        any_data = False
        
        for q_idx, q in enumerate(Q_VALUES):
            areas_table[lattice.name_str][q] = {}
            
            for t_idx, T in enumerate(T_VALUES_HYSTERESIS):
                ax = axes[q_idx, t_idx]
                
                # Construir candidatos de nombres de archivo
                q_str = q_str = str(int(q * 100))
                T_str = str(int(T * 10))
                
                candidates = build_filename_candidates(
                    HYST_DATA_DIR,
                    'hysteresis_{lattice}_q{q}_T{T}.bin',
                    lattice=lattice.name_str,
                    q=q_str,
                    T=T_str
                )
                
                # Buscar archivo
                filepath = try_multiple_filenames(candidates)
                
                if filepath is None:
                    ax.text(0.5, 0.5, f'No disponible\nq={q}, T={T}',
                           ha='center', va='center', transform=ax.transAxes, fontsize=10)
                    ax.set_xticks([])
                    ax.set_yticks([])
                    areas_table[lattice.name_str][q][T] = None
                    continue
                
                # Leer datos
                data = reader.read_hysteresis(str(filepath))
                if data is None:
                    ax.text(0.5, 0.5, f'Error lectura\nq={q}, T={T}',
                           ha='center', va='center', transform=ax.transAxes, fontsize=10)
                    ax.set_xticks([])
                    ax.set_yticks([])
                    areas_table[lattice.name_str][q][T] = None
                    continue
                
                any_data = True
                
                # Plotear histéresis
                ax.plot(data['H_down'], data['m_down'], 'o-', 
                       label='Descendente (H↓)', linewidth=2, markersize=5, 
                       color='#1f77b4', alpha=0.8)
                ax.plot(data['H_up'], data['m_up'], 's-', 
                       label='Ascendente (H↑)', linewidth=2, markersize=5, 
                       color='#ff7f0e', alpha=0.8)
                
                # Rellenar áreas
                ax.fill_between(data['H_down'], data['m_down'], 0, 
                              alpha=0.2, color='#1f77b4')
                
                H_up_sorted = np.sort(data['H_up'])
                m_up_sorted = data['m_up'][np.argsort(data['H_up'])]
                ax.fill_between(H_up_sorted, m_up_sorted, 0, 
                              alpha=0.2, color='#ff7f0e')
                
                # Calcular área
                area = calculate_hysteresis_area(data['H_down'], data['m_down'], 
                                                  data['H_up'], data['m_up'])
                areas_table[lattice.name_str][q][T] = area
                
                ax.text(0.05, 0.95, f'Área = {area:.3f}', transform=ax.transAxes,
                        fontsize=9, verticalalignment='top',
                        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
                
                # Etiquetas
                ax.set_xlabel(r'Campo Magnético $H$ ($J$)', fontsize=10)
                ax.set_ylabel(r'Magnetización $m$', fontsize=10)
                ax.set_title(f'$q = {q}$, $T = {T:.1f}$', fontsize=11, pad=8)
                ax.legend(fontsize=9, loc='lower right')
                ax.grid(True, alpha=0.3, linestyle='--')
                ax.axhline(y=0, color='k', linewidth=0.5, alpha=0.3)
                ax.axvline(x=0, color='k', linewidth=0.5, alpha=0.3)
        
        if not any_data:
            print(f"    ⚠ Sin datos de histéresis para {lattice.name_str}")
            plt.close()
            continue
        
        plt.tight_layout()
        out_path = OUTPUT_DIR / 'ferromagnetismo' / f'03_hysteresis_{lattice.name_str}.png'
        out_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='white')
        plt.close()
        print(f"    ✓ 03_hysteresis_{lattice.name_str}.png")
    
    # Imprimir tabla resumen
    print("\n  " + "="*60)
    print("  TABLA RESUMEN: ÁREAS DE HISTÉRESIS")
    print("  " + "="*60)
    for lattice in LATTICES:
        print(f"\n  {lattice.name_str.upper()}:")
        print(f"  {'q':<8} | {'T=2.0':<12} | {'T=5.0':<12}")
        print("  " + "-"*40)
        for q in Q_VALUES:
            row = f"  {q:<8} |"
            for T in T_VALUES_HYSTERESIS:
                area = areas_table.get(lattice.name_str, {}).get(q, {}).get(T)
                if area is not None:
                    row += f" {area:<12.4f}|"
                else:
                    row += f" {'N/A':<12}|"
            print(row)
    print("  " + "="*60 + "\n")

def plot_phase_transitions() -> None:
    """
    Genera gráficas de transiciones de fase con BARRAS DE ERROR,
    Línea de Tc y REGIÓN DE INCERTIDUMBRE.
    Versión DEBUG para forzar la aparición de líneas.
    """
    print_header("4. TRANSICIONES DE FASE (Con Tc y Región Crítica)")
    
    reader = IsingDataReader()
    
    # ============================================================
    # 1. CARGAR DATOS DE Tc
    # ============================================================
    tc_file = TC_DATA_DIR / 'tc_analysis.bin'
    tc_map = {}
    
    tc_data = reader.read_tc_analysis(str(tc_file))
    if tc_data:
        print("    ✓ Datos de Tc cargados.")
        lattices_in = tc_data['lattices']
        qs_in = tc_data['q_values']
        matrix = tc_data['tc_matrix']
        
        # Crear mapa simple
        for i, lat in enumerate(lattices_in):
            tc_map[lat] = {}
            for j, q_val in enumerate(qs_in):
                tc_map[lat][float(q_val)] = float(matrix[i][j])
        
        # DEBUG: Mostrar qué tenemos en memoria
        print(f"    DEBUG: Redes disponibles en Tc: {list(tc_map.keys())}")
        #print(f"    EBUG: Qs disponibles para 'square': {list(tc_map.get('square', {}).keys())}")
    else:
        print("    ⚠ No se encontró archivo de análisis Tc.")

    # ============================================================
    # 2. CONFIGURAR FIGURAS
    # ============================================================
    fig_m, axes_m = plt.subplots(2, 2, figsize=(13, 10))
    fig_m.suptitle(r'Magnetización vs Temperatura con q=0.8 \n (con $T_c$ y Región Crítica)', fontsize=14, fontweight='bold', y=0.995)
    axes_m = axes_m.flatten()
    
    fig_e, axes_e = plt.subplots(2, 2, figsize=(13, 10))
    fig_e.suptitle(r'Energía vs Temperatur con q = 0.8 \n (con $T_c$ y Región Crítica)', fontsize=14, fontweight='bold', y=0.995)
    axes_e = axes_e.flatten()
    
    # INTENTAMOS ENCONTRAR EL Q = 0.8 (o el más cercano)
    target_q = 0.8 
    
    for idx, lattice in enumerate(LATTICES):
        print(f"  Procesando {lattice.name_str}...")
        filename = TRANS_DATA_DIR / f'transicion_{lattice.name_str}.bin'
        
        # --- Lectura de Datos (RAW o Normal) ---
        data_raw = reader.read_phase_transition_raw(str(filename))
        
        T, m_mean, m_err, E_mean, E_err = [], [], [], [], []
        has_error_bars = False
        
        if data_raw is not None:
            T = data_raw['T']
            m_raw_matrix = data_raw['m_raw']
            e_raw_matrix = data_raw['E_raw']
            for i in range(len(T)):
                m_series = np.abs(m_raw_matrix[i])
                m_mean.append(np.mean(m_series))
                m_err.append(calculate_blocking_error(m_series))
                e_series = e_raw_matrix[i]
                E_mean.append(np.mean(e_series))
                E_err.append(calculate_blocking_error(e_series))
            has_error_bars = True
        else:
            data = reader.read_phase_transition(str(filename))
            if data is None:
                axes_m[idx].text(0.5, 0.5, 'N/A'); axes_e[idx].text(0.5, 0.5, 'N/A'); continue
            T = data['T']
            m_mean = np.abs(data['m'])
            E_mean = data['E']
            m_err = E_err = None
            has_error_bars = False

        # --- BUSCAR TC DE FORMA FLEXIBLE ---
        tc_val = None
        if lattice.name_str in tc_map:
            # Buscar el q más cercano a 0.8 disponible en el mapa
            available_qs = list(tc_map[lattice.name_str].keys())
            if available_qs:
                # Encuentra el q que minimiza la diferencia con 0.8
                closest_q = min(available_qs, key=lambda x: abs(x - target_q))
                
                if abs(closest_q - target_q) < 0.1: # Tolerancia de 0.1
                    tc_candidate = tc_map[lattice.name_str][closest_q]
                    if tc_candidate > 0.1: # Solo si Tc es válido (>0)
                        tc_val = tc_candidate
                        print(f"    ➜ Encontrada Tc={tc_val:.2f} para {lattice.name_str} (usando q={closest_q:.1f})")
                    else:
                        print(f"    ⚠ Tc es 0.0 para {lattice.name_str}, posible error en C++")
                else:
                    print(f"    ⚠ No hay datos cercanos a q={target_q} (Disponibles: {available_qs})")
            else:
                print(f"    ⚠ Lista de Q vacía para {lattice.name_str}")
        else:
            print(f"    ⚠ {lattice.name_str} no está en tc_map")

        # --- GRAFICAR MAGNETIZACIÓN ---
        ax = axes_m[idx]
        if has_error_bars:
            ax.errorbar(T, m_mean, yerr=m_err, fmt='o-', linewidth=1.5, markersize=4,
                       color='#2ca02c', label='|m| ± error', ecolor='black', capsize=3, alpha=0.8)
        else:
            ax.plot(T, m_mean, 'o-', linewidth=2.5, markersize=6,
                   color='#2ca02c', label='|m|', alpha=0.8)

        # DIBUJAR LINEAS Y SOMBRA (Si encontramos Tc)
        if tc_val is not None:
            # Línea
            ax.axvline(tc_val, color='blue', linestyle='--', alpha=0.9, linewidth=1.5, 
                      label=f'$T_c \\approx {tc_val:.2f}$')
            # Sombra
            uncertainty = 0.05  # Ajustar según el análisis de Tc
            ax.axvspan(tc_val - uncertainty, tc_val + uncertainty, 
                      color='blue', alpha=0.15, label='Región Crítica')
        
        ax.set_xlabel(r'$T$ ($J/k_B$)'); ax.set_ylabel(r'$|m|$')
        ax.set_title(f'{lattice.name_str.upper()} ($z={lattice.coordination_number}$)')
        ax.grid(True, alpha=0.3, linestyle='--'); ax.legend(fontsize=9)

        # --- GRAFICAR ENERGÍA ---
        ax = axes_e[idx]
        if has_error_bars:
            ax.errorbar(T, E_mean, yerr=E_err, fmt='o-', linewidth=1.5, markersize=4,
                       color='#d62728', label='E ± error', ecolor='black', capsize=3, alpha=0.8)
        else:
            ax.plot(T, E_mean, 'o-', linewidth=2.5, markersize=6,
                   color='#d62728', label='E', alpha=0.8)

        if tc_val is not None:
            ax.axvline(tc_val, color='blue', linestyle='--', alpha=0.9, linewidth=1.5,
                       label=f'$T_c \\approx {tc_val:.2f}$')
            ax.axvspan(tc_val - uncertainty, tc_val + uncertainty, color='blue', alpha=0.1, 
                       label='Región Crítica')

        ax.set_xlabel(r'$T$ ($J/k_B$)'); ax.set_ylabel(r'$E$ ($J$)')
        ax.set_title(f'{lattice.name_str.upper()}')
        ax.grid(True, alpha=0.3, linestyle='--'); ax.legend(fontsize=9)

    # GUARDAR
    plt.figure(fig_m.number); plt.tight_layout()
    out_path_m = OUTPUT_DIR / 'transicion' / '04_magnetizacion_transicion.png'
    out_path_m.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path_m, dpi=150, bbox_inches='tight')
    plt.close(fig_m)
    print(f"  ✓ Guardado: {out_path_m.name}")
    
    plt.figure(fig_e.number); plt.tight_layout()
    out_path_e = OUTPUT_DIR / 'transicion' / '04_energia_transicion.png'
    out_path_e.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path_e, dpi=150, bbox_inches='tight')
    plt.close(fig_e)
    print(f"  ✓ Guardado: {out_path_e.name}")

def plot_snapshots() -> None:
    """Genera visualizaciones de snapshots agrupadas por lattice (3 q × 3 T)"""
    print_header("5. VISUALIZACIÓN DE CONFIGURACIONES")
    
    reader = IsingDataReader()
    snap_count = 0
    
    for lattice in LATTICES:
        print(f"  Generando para {lattice.name_str.upper()}...")
        
        T_values = [1.0, 3.0, 5.0]
        Q_list = Q_VALUES
        
        # Crear figura agrupada por lattice
        if lattice.name_str == 'bcc':
            # BCC: 9×8 (3 q × 3 T × 8 capas)
            fig, axes = plt.subplots(9, 8, figsize=(20, 18))
            fig.suptitle(
                f'Configuración de Spins - BCC (z=8, 3 q × 3 T × 8 capas)\n'
                r'Sistema natural: $k_B = 1$, $J = 1$',
                fontsize=16, fontweight='bold', y=0.995
            )
        elif lattice.name_str == 'chain':
            # Chain: 3 filas (q) x 3 columnas (T)
            # Ajustamos figsize para que las tiras se vean alargadas y bonitas
            fig, axes = plt.subplots(3, 3, figsize=(15, 6)) # Altura reducida para enfatizar formato "cinta"
            fig.suptitle(
                f'Configuración de Spins - Chain (z=2, 3 q × 3 T)\n'
                r'Sistema natural: $k_B = 1$, $J = 1$',
                fontsize=16, fontweight='bold', y=0.995
            )
        elif lattice.name_str == 'honeycomb':
            fig, axes = plt.subplots(3, 3, figsize=(15, 12))
            fig.suptitle(
                f'Configuración de Spins - Honeycomb (z=3, 3 q × 3 T)\n'
                r'Sistema natural: $k_B = 1$, $J = 1$',
                fontsize=16, fontweight='bold', y=0.995
            )
        elif lattice.name_str == 'square':
            fig, axes = plt.subplots(3, 3, figsize=(15, 12))
            fig.suptitle(
                f'Configuración de Spins - Square (z=4, 3 q × 3 T)\n'
                r'Sistema natural: $k_B = 1$, $J = 1$',
                fontsize=16, fontweight='bold', y=0.995
            )
        
        any_data = False
        
        # =====================================================
        # LÓGICA ESPECÍFICA POR RED
        # =====================================================
        
        if lattice.name_str == 'bcc':
            for q_idx, q in enumerate(Q_list):
                q_str = format_q_str(q)
                for t_idx, T in enumerate(T_values):
                    row = q_idx * 3 + t_idx
                    
                    for col, k in enumerate(range(8)):
                        ax = axes[row, col]
                        candidates = [
                            SNAP_DATA_DIR / f'snapshot_{lattice.name_str}_T{int(T)}.bin',
                            SNAP_DATA_DIR / f'snapshot_{lattice.name_str}_q{q_str}_T{int(T)}.bin',
                            SNAP_DATA_DIR / f'snapshot_{lattice.name_str}_q{int(q*100)}_T{int(T)}.bin',
                        ]
                        filepath = try_multiple_filenames(candidates)
                        
                        if not filepath or not (data := reader.read_snapshot(str(filepath))):
                            ax.set_axis_off()
                            continue
                        
                        any_data = True
                        nx, ny, nz = data['shape']
                        grid = np.zeros((nx, ny), dtype=np.int8)
                        for spin_val, coord in zip(data['spins'], data['coordinates']):
                            i, j, coord_k = coord
                            if coord_k == k:
                                grid[i, j] = int(spin_val)
                        
                        ax.imshow(grid, cmap='RdBu', vmin=-1, vmax=1, origin='lower')
                        ax.set_xticks([]); ax.set_yticks([])
                        if col == 0: ax.set_ylabel(f'q={q}\nT={T}', fontsize=7)
                        if row == 8: ax.set_xlabel(f'k={k}', fontsize=7)
                        if row == 0: ax.set_title(f'k={k}', fontsize=8)

        elif lattice.name_str == 'chain':
            # =====================================================
            # CHAIN: MODO "CÓDIGO DE BARRAS" (HEATMAP STRIP)
            # =====================================================
            for q_idx, q in enumerate(Q_list):
                q_str = format_q_str(q)
                for t_idx, T in enumerate(T_values):
                    row, col = q_idx, t_idx
                    ax = axes[row, col]
                    
                    candidates = [
                        SNAP_DATA_DIR / f'snapshot_{lattice.name_str}_T{int(T)}.bin',
                        SNAP_DATA_DIR / f'snapshot_{lattice.name_str}_q{q_str}_T{int(T)}.bin',
                        SNAP_DATA_DIR / f'snapshot_{lattice.name_str}_q{int(q*100)}_T{int(T)}.bin',
                    ]
                    
                    filepath = try_multiple_filenames(candidates)
                    if not filepath:
                        ax.text(0.5, 0.5, 'No Data', ha='center', va='center')
                        ax.axis('off')
                        continue
                        
                    data = reader.read_snapshot(str(filepath))
                    if not data:
                        ax.text(0.5, 0.5, 'Error', ha='center', va='center')
                        ax.axis('off')
                        continue
                    
                    any_data = True
                    temp = data.get('T', T)
                    shape = data['shape']
                    
                    # 1. Reconstruir array 1D
                    chain_values = np.zeros(shape[0], dtype=np.int8)
                    for spin_val, coord in zip(data['spins'], data['coordinates']):
                        chain_values[coord[0]] = int(spin_val)
                    
                    # 2. Convertir a matriz (1, L) para imshow
                    # Esto crea una "tira" de 1 pixel de alto por L de ancho
                    chain_grid = chain_values.reshape(1, -1)
                    
                    # 3. Graficar como Mapa de Calor
                    # aspect='auto' es CRUCIAL: fuerza a que la tira llene todo el ancho del plot
                    im = ax.imshow(chain_grid, cmap='RdBu', vmin=-1, vmax=1, aspect='auto', interpolation='nearest')
                    
                    # 4. Limpieza estética
                    ax.set_yticks([]) # Quitamos eje Y porque no aporta nada en 1D
                    
                    if row == 0: ax.set_title(f'T = {temp:.1f}', fontsize=12)
                    if col == 0: ax.set_ylabel(f'q = {q}', fontsize=12, fontweight='bold')
                    if row == 2: ax.set_xlabel('Sitios del sistema (i)', fontsize=10)
                    else: ax.set_xticks([]) # Solo mostrar eje X abajo del todo

        elif lattice.name_str in ['square', 'honeycomb']:
            # =====================================================
            # SQUARE / HONEYCOMB: Lógica estándar 2D
            # =====================================================
            for q_idx, q in enumerate(Q_list):
                q_str = format_q_str(q)
                for t_idx, T in enumerate(T_values):
                    row, col = q_idx, t_idx
                    ax = axes[row, col]
                    
                    candidates = [
                        SNAP_DATA_DIR / f'snapshot_{lattice.name_str}_T{int(T)}.bin',
                        SNAP_DATA_DIR / f'snapshot_{lattice.name_str}_q{q_str}_T{int(T)}.bin',
                        SNAP_DATA_DIR / f'snapshot_{lattice.name_str}_q{int(q*100)}_T{int(T)}.bin',
                    ]
                    
                    filepath = try_multiple_filenames(candidates)
                    if not filepath or not (data := reader.read_snapshot(str(filepath))):
                        ax.text(0.5, 0.5, 'N/A', ha='center')
                        ax.axis('off')
                        continue

                    any_data = True
                    grid = np.zeros(data['shape'], dtype=np.int8)
                    for spin_val, coord in zip(data['spins'], data['coordinates']):
                        grid[coord[0], coord[1]] = int(spin_val)
                    
                    ax.imshow(grid, cmap='RdBu', vmin=-1, vmax=1, origin='lower')
                    
                    if row == 0: ax.set_title(f'T = {data.get("T", T):.1f}', fontsize=12)
                    if col == 0: ax.set_ylabel(f'q = {q}', fontsize=12, fontweight='bold')
                    ax.set_xticks([]); ax.set_yticks([])

        if not any_data:
            print(f"    ⚠ Sin datos de snapshots para {lattice.name_str}")
            plt.close()
            continue
        
        plt.tight_layout()
        out_path = OUTPUT_DIR / 'snapshots' / f'05_snapshot_{lattice.name_str}_q.png'
        out_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='white')
        plt.close()
        print(f"    ✓ 05_snapshot_{lattice.name_str}_q.png")


def plot_reduced_magnetization() -> None:
    """
    Genera la gráfica de m vs H/T con BARRAS DE ERROR y realiza un AJUSTE
    cuantitativo calculando Chi-cuadrado para validar la ley.
    """
    print_header("1.5. LEY DE ESTADOS CORRESPONDIENTES (con Ajuste Chi^2)")
    
    reader = IsingDataReader()
    
    # Modelo teórico para el ajuste
    def tanh_func(x, a, b):
        return a * np.tanh(b * x)
    
    for lattice in LATTICES:
        print(f"  Generando y ajustando para {lattice.name_str.upper()}...")
        
        fig, axes = plt.subplots(1, 3, figsize=(16, 5))
        fig.suptitle(
            f'Ley de Estados Correspondientes $m$ vs $H/T$ - {lattice.name_str.upper()} ($z = {lattice.coordination_number}$)\n'
            r'Modelo: $m = a \cdot \tanh(b \cdot H/T)$',
            fontsize=14, fontweight='bold', y=1.02
        )
        
        for idx, q in enumerate(Q_VALUES):
            q_str = format_q_str(q)
            filename = PARA_DATA_DIR / f'{lattice.name_str}_{q_str}.bin'
            data_raw = reader.read_paramagnetism_raw(str(filename))
            
            ax = axes[idx]
            
            # Listas para acumular TODOS los puntos (colapso) para el ajuste
            all_x = []
            all_y = []
            all_err = []
            
            if data_raw is not None:
                T_arr = data_raw['T']
                H_arr = data_raw['H']
                m_matrix = data_raw['m_raw']
                
                cmap = plt.get_cmap('viridis')
                colors = cmap(np.linspace(0, 1, len(T_arr)))
                
                # 1. RECOLECTAR Y GRAFICAR DATOS
                for t_i, T in enumerate(T_arr):
                    if T == 0: continue
                    
                    # Eje X reducido y Eje Y
                    x_reduced = H_arr / T
                    m_means = []
                    m_errs = []
                    
                    for h_i in range(len(H_arr)):
                        series = m_matrix[t_i, h_i, :]
                        val = np.mean(series)
                        err = calculate_blocking_error(series)
                        
                        m_means.append(val)
                        m_errs.append(err)
                        
                        # Acumular para el ajuste global (solo si el error es válido > 0)
                        if err > 1e-9: 
                            all_x.append(H_arr[h_i] / T)
                            all_y.append(val)
                            all_err.append(err)
                    
                    ax.errorbar(x_reduced, m_means, yerr=m_errs, fmt='o', 
                               color=colors[t_i], ecolor='black', capsize=2, 
                               markersize=3, alpha=0.5, label=f'T={T:.1f}' if idx==0 else "")

                # 2. REALIZAR EL AJUSTE Y CALCULAR CHI-CUADRADO
                if len(all_x) > 10:
                    x_fit_data = np.array(all_x)
                    y_fit_data = np.array(all_y)
                    sigma_data = np.array(all_err)
                    
                    try:
                        # Ajuste no lineal ponderado por el error (sigma)
                        popt, pcov = curve_fit(tanh_func, x_fit_data, y_fit_data, 
                                             p0=[1.0, 1.0], sigma=sigma_data, absolute_sigma=True)
                        
                        a_opt, b_opt = popt
                        
                        # Calcular Chi-Cuadrado
                        y_model = tanh_func(x_fit_data, *popt)
                        residuals = y_fit_data - y_model
                        chi_sq = np.sum((residuals / sigma_data) ** 2)
                        
                        # Chi-Cuadrado Reducido (dof = N - numero_parametros)
                        dof = len(x_fit_data) - 2
                        chi_sq_red = chi_sq / dof
                        
                        # Graficar la línea ajustada
                        x_smooth = np.linspace(min(all_x), max(all_x), 200)
                        ax.plot(x_smooth, tanh_func(x_smooth, *popt), 'r--', linewidth=2, 
                               label=f'Ajuste')
                        
                        # Mostrar estadísticas en un cuadro de texto
                        textstr = '\n'.join((
                            r'$a=%.3f$' % (a_opt, ),
                            r'$b=%.3f$' % (b_opt, ),
                            r'$\chi^2_\nu=%.2f$' % (chi_sq_red, )))
                        
                        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
                        ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=9,
                               verticalalignment='top', bbox=props)
                        
                    except Exception as e:
                        print(f"    ⚠ Falló el ajuste en {lattice.name_str} q={q}: {e}")
            else:
                ax.text(0.5, 0.5, 'Sin datos RAW\nNo hay ajuste', ha='center', transform=ax.transAxes)

            ax.set_xlabel(r'$H/T$', fontsize=11)
            ax.set_ylabel(r'$m$', fontsize=11)
            ax.set_title(f'$q = {q}$', fontsize=12)
            ax.grid(True, alpha=0.3)
            #ax.set_xlim(-2, 2)
            if idx == 2: ax.legend(loc='lower right', fontsize=8) # Leyenda solo en el último
            
        plt.tight_layout()
        out_path = OUTPUT_DIR / 'paramagnetismo' / f'01b_ley_de_estados_correspondientes_fit_{lattice.name_str}.png'
        plt.savefig(out_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"    ✓ Guardado: {out_path.name}")


def plot_relaxation() -> None:
    """
    Genera gráficas de Relajación Dinámica (Energía vs Tiempo).
    - Una imagen por Topología (Lattice).
    - Subplots (columnas) para cada valor de q.
    - Curvas de distintos colores para cada Temperatura T.
    - Solo considera J=1 (Ferromagnético).
    """
    print_header("6. RELAJACIÓN DINÁMICA - AGRUPADO POR TOPOLOGÍA")
    
    reader = IsingDataReader()
    
    # Importamos constantes
    from src.config import T_VALUES_RELAXATION
    from src.utils import build_filename_candidates
    
    # Directorio de salida
    output_folder = OUTPUT_DIR / 'relajacion'
    output_folder.mkdir(parents=True, exist_ok=True)
    
    # Parámetros fijos
    J = 1  # Solo nos interesa el caso ferromagnético
    
    # Iterar por cada red (Una imagen por red)
    for lattice in LATTICES:
        print(f"  Generando panel para {lattice.name_str.upper()}...")
        
        # Configuramos la figura: 1 fila, N columnas (una por cada q)
        n_cols = len(Q_VALUES)
        fig, axes = plt.subplots(nrows=1, ncols=n_cols, figsize=(5 * n_cols, 6), 
                                 sharey=True, constrained_layout=True)
        
        # Si solo hay un q, axes no es una lista, lo convertimos para que el bucle funcione
        if n_cols == 1:
            axes = [axes]
            
        fig.suptitle(
            f'Relajación Dinámica - {lattice.name_str.upper()} ($z={lattice.coordination_number}$)\n'
            r'Régimen Ferromagnético ($J=1$)',
            fontsize=16, fontweight='bold'
        )
        
        found_data_for_lattice = False
        
        # 1. Obtener el objeto del mapa de color de forma segura
        cmap = plt.get_cmap('plasma')
            
        # 2. Generar los colores
        colors = cmap(np.linspace(0, 0.85, len(T_VALUES_RELAXATION)))

        # Iterar sobre cada valor de q (Cada q es un subplot)
        for i, q in enumerate(Q_VALUES):
            ax = axes[i]
            q_str = str(int(q * 100)) # Corrección de formato (0.5 -> "50")
            
            ax.set_title(f'Probabilidad de recableado $q = {q}$', fontsize=12)
            
            found_data_in_subplot = False
            
            # Iterar sobre Temperaturas (Curvas en el subplot)
            for t_idx, T in enumerate(T_VALUES_RELAXATION):
                T_str = str(int(T * 10))
                
                candidates = build_filename_candidates(
                    RELAX_DATA_DIR,
                    'relax_{lattice}_q{q}_T{T}_J{J}.bin',
                    lattice=lattice.name_str,
                    q=q_str,
                    T=T_str,
                    J=str(J)
                )
                
                filepath = try_multiple_filenames(candidates)
                
                if filepath is None:
                    continue
                
                data = reader.read_relaxation(str(filepath))
                
                if data is not None:
                    found_data_in_subplot = True
                    found_data_for_lattice = True
                    
                    E = data['E']
                    steps = np.arange(len(E))
                    color = colors[t_idx]
                    
                    # 1. Traza cruda (muy transparente, fondo)
                    ax.plot(steps, E, '-', linewidth=0.5, alpha=0.1, color=color)
                    
                    # 2. Media móvil (suavizado)
                    window = 100 # Ventana de suavizado (ajusta si quieres más/menos suave)
                    if len(E) > window:
                        E_smooth = np.convolve(E, np.ones(window)/window, mode='valid')
                        # Ajustar eje x para centrar
                        steps_smooth = steps[window-1:]
                        ax.plot(steps_smooth, E_smooth, '-', linewidth=2, 
                                label=f'T={T:.1f}', color=color)
            
            # Decoración del subplot
            ax.set_xlabel('Pasos de Monte Carlo (MCS)', fontsize=11)
            if i == 0:
                ax.set_ylabel('Energía Total $E$ ($J$)', fontsize=11)
            
            ax.grid(True, alpha=0.3, linestyle='--')
            
            if found_data_in_subplot:
                # Solo ponemos leyenda en el primer subplot para no saturar
                if i == 0: 
                    ax.legend(title="Temperatura ($T$)", fontsize=9, loc='lower right')
            else:
                ax.text(0.5, 0.5, "Sin datos", ha='center', va='center', transform=ax.transAxes)

        # Guardar resultado
        if found_data_for_lattice:
            out_path = output_folder / f'06_relax_combined_{lattice.name_str}.png'
            plt.savefig(out_path, dpi=150, facecolor='white')
            print(f"    ✓ Guardado: {out_path.name}")
        else:
            print(f"    ⚠ No se encontraron datos para {lattice.name_str}")
        
        plt.close(fig)



def main() -> int:
    """Ejecuta el análisis completo"""
    
    print("\n" + "="*70)
    print("🔬 ANÁLISIS DEL MODELO ISING HÍBRIDO")
    print(r"Sistema de unidades natural: $k_B = 1$, $J = 1$")
    print("="*70)
    
    if not verify_data_exists():
        return 1
    
    try:
        plot_paramagnetism()
        plot_reduced_magnetization()
        plot_susceptibility()
        plot_hysteresis()
        plot_phase_transitions()
        plot_snapshots()
        plot_relaxation()
        
        print("\n" + "="*70)
        print("✅ ANÁLISIS COMPLETADO")
        print(f"📁 Gráficas guardadas en: {OUTPUT_DIR}")
        print("="*70 + "\n")
        
        return 0
        
    except Exception as e:
        print(f"\n❌ ERROR: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())