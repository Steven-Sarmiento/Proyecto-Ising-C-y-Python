# Simulación Monte Carlo del Modelo de Ising en Redes Cristalinas

Este repositorio contiene la implementación computacional y el análisis de resultados del "Miniproyecto de Física Estadística", enfocado en la simulación del modelo de Ising con dinámica de Metropolis en diversas topologías y grados de dilución magnética.

## 👤 Autor

**Jorge Steven Sarmiento Arboleda** Instituto de Física, Universidad de Antioquia, Medellín, Colombia  
Curso: Física Estadística (0302470)

---

## 1. Descripción del Proyecto

El proyecto estudia el comportamiento termodinámico de espines clásicos interactuantes bajo el Hamiltoniano de Ising:

$$\mathcal{H} = -J \sum_{\langle i, j \rangle} s_i s_j - H \sum_i s_i$$

Se investigan sistemáticamente los efectos de la **dimensionalidad** y la **dilución magnética** ($q$) en las transiciones de fase y la histéresis.

### Parámetros del Estudio

- **Topologías ($z$):** Cadena 1D ($z=2$), Honeycomb 2D ($z=3$), Cuadrada 2D ($z=4$) y BCC 3D ($z=8$).
- **Tamaños del sistema ($L$):** $L=500$ (1D), $L=25$ (2D), $L=8$ (3D).
- **Dilución ($q$):** Probabilidad de ocupación de espines $q \in \{0.5, 0.8, 1.0\}$.
- **Regímenes:**
  - Paramagnético ($J=0$)
  - Ferromagnético ($J=1$)

---

## 📁 Estructura del Proyecto

```text
proyecto_ising_hibrido/
├── cpp/                       # Núcleo de simulación (C++)
│   ├── ising_model.hpp        # Clase IsingModel
│   ├── ising_model.cpp        # Lógica de Metropolis y observables
│   ├── main.cpp               # Orquestador de simulaciones
│   ├── ising_sim              # Ejecutable compilado
│   └── data/                  # Salida de datos binarios
├── python/                    # Análisis y Graficado (Python)
│   ├── main.py                # Script principal de generación de gráficas
│   ├── src/
│   │   ├── config.py          # Rutas y constantes físicas
│   │   ├── reader.py          # Lectura de binarios C++
│   │   ├── plotter.py         # Motores de graficado científico
│   │   └── utils.py           # Análisis estadístico (Blocking Method)
└── graficas/                  # Resultados visuales (Output)
    ├── paramagnetismo/
    ├── ferromagnetismo/       # Ciclos de histéresis
    ├── transicion/            # Curvas M vs T y Energía
    ├── snapshots/             # Configuraciones de espines
    └── relajacion/            # Evolución temporal
## 🚀 Instalación y Uso

### 1. Compilar C++

```bash
cd cpp
g++ -o ising_sim main.cpp ising_model.cpp -std=c++17 -O2
```

### 2. Ejecutar Simulaciones

```bash
./ising_sim
```

Esto genera los datos binarios en `cpp/data/`.

### 3. Configurar Python

```bash
cd ../python
python3 -m venv venv
source venv/bin/activate  # Linux/Mac
pip install -r requirements.txt
```

### 4. Generar Gráficas

```bash
python main.py
```

Las gráficas se guardan en `python/graficas/`.

## 📦 Dependencias

**C++:**

- Compilador con soporte C++17 (g++ 7+ o clang++ 5+)
- Librería estándar (`<filesystem>`, `<random>`, `<vector>`)

**Python (`requirements.txt`):**

```text
numpy
matplotlib
scipy
```

## 📐 Sistema de Unidades

El código usa unidades naturales:

- $k_B = 1$
- $J = 1$ (constante de acoplamiento)
- Temperatura en unidades de $J/k_B$
- Campo magnético en unidades de $J$

## 🧮 Métodos Implementados

### Simulaciones (C++)

| Simulación | Descripción | Termalización |
|------------|-------------|---------------|
| Relajación | Evolución temporal de E desde estado inicial | No |
| Snapshots | Configuraciones de espines en equilibrio | Sí (escala con z) |
| Paramagnetismo | Curvas m vs H para J=0 | Mínima |
| Histéresis | Ciclos m vs H para J=1 | Solo al inicio |
| Transición de fase | Curvas m vs T, estimación de Tc | Sí (escala con z) |

### Análisis (Python)

- **Ley de Estados Correspondientes**: Colapso de datos m vs H/T
- **Ajuste Tanh**: $m = a \cdot \tanh(b \cdot H/T)$ con bondad de ajuste $\chi^2$
- **Temperatura Crítica**: Estimada por máxima pendiente de m(T)

## 📊 Ejemplos de Resultados

### Paramagnetismo

- Curvas m vs H para diferentes T
- Verificación de $m = \tanh(H/T)$

### Ferromagnetismo

- Ciclos de histéresis con área no nula
- Transición de fase con identificación de Tc

### Efecto de la Topología

- Mayor z → Mayor Tc
- Mayor z → Histéresis más pronunciada

## 📝 Notas Técnicas

- Datos almacenados en formato binario little-endian (float32)
- Condiciones de frontera periódicas en todas las redes
- Algoritmo de Metropolis con selección aleatoria de sitios
- Errores estadísticos calculados por método de blocking

## 📄 Licencia

Este código es parte de un trabajo académico. Se permite su uso y modificación con la debida atribución.

Este proyecto utilizó herramientas de IA (Claude, DeepSeek, Gemini, ChatGPT y NotebookLM) como soporte auxiliar en la depuración de código y refinamiento de redacción, bajo supervisión humana crítica.

---

**Versión**: 1.0  
**Fecha**: 2025
