#include "ising_model.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <filesystem>
#include <cstring>
#include <algorithm>
#include <map>

namespace fs = std::filesystem;

// ====================================================================
// GUARDADO DE DATOS BINARIOS
// ====================================================================

void save_paramagnetism_data(const std::string& filename,
                            const IsingModel& model,
                            const std::vector<float>& temperatures,
                            const std::vector<float>& fields,
                            const std::vector<std::vector<float>>& m_data,
                            const std::vector<std::vector<float>>& chi_data) {
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    uint32_t L = model.get_L();
    float q = model.get_q();
    float J = model.get_J();
    uint32_t n_T = temperatures.size();
    uint32_t n_H = fields.size();
    
    std::string lattice_str = model.get_lattice_type_str();
    char lattice_arr[32] = {0};
    std::strncpy(lattice_arr, lattice_str.c_str(), 31);
    
    file.write(reinterpret_cast<char*>(&L), sizeof(uint32_t));
    file.write(reinterpret_cast<char*>(&q), sizeof(float));
    file.write(reinterpret_cast<char*>(&J), sizeof(float));
    file.write(reinterpret_cast<char*>(&n_T), sizeof(uint32_t));
    file.write(reinterpret_cast<char*>(&n_H), sizeof(uint32_t));
    file.write(lattice_arr, 32);
    
    file.write(reinterpret_cast<const char*>(temperatures.data()),
               temperatures.size() * sizeof(float));
    file.write(reinterpret_cast<const char*>(fields.data()),
               fields.size() * sizeof(float));
    
    for (const auto& row : m_data) {
        file.write(reinterpret_cast<const char*>(row.data()),
                   row.size() * sizeof(float));
    }
    
    for (const auto& row : chi_data) {
        file.write(reinterpret_cast<const char*>(row.data()),
                   row.size() * sizeof(float));
    }
    
    file.close();
    std::cout << "✓ " << filename << " (" << fs::file_size(filename) / 1024.0 << " KB)" << std::endl;
}

void save_hysteresis_data(const std::string& filename,
                         const std::vector<float>& H_down, const std::vector<float>& m_down,
                         const std::vector<float>& H_up, const std::vector<float>& m_up) {
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    uint32_t n = H_down.size();
    file.write(reinterpret_cast<char*>(&n), sizeof(uint32_t));
    file.write(reinterpret_cast<const char*>(H_down.data()), n * sizeof(float));
    file.write(reinterpret_cast<const char*>(m_down.data()), n * sizeof(float));
    file.write(reinterpret_cast<const char*>(H_up.data()), n * sizeof(float));
    file.write(reinterpret_cast<const char*>(m_up.data()), n * sizeof(float));
    
    file.close();
    std::cout << "  ✓ " << filename << std::endl;
}

// ====================================================================
// NUEVA FUNCIÓN: Guardado de DATOS CRUDOS (RAW) para Barras de Error
// ====================================================================
void save_phase_transition_raw(const std::string& filename,
                               const std::vector<float>& temperatures,
                               const std::vector<std::vector<float>>& m_raw,
                               const std::vector<std::vector<float>>& e_raw) {
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    uint32_t n_T = (uint32_t)temperatures.size();
    // Asumimos que todas las T tienen el mismo número de muestras
    uint32_t n_samples = m_raw.empty() ? 0 : (uint32_t)m_raw[0].size(); 
    
    // 1. Escribir dimensiones
    file.write(reinterpret_cast<const char*>(&n_T), sizeof(uint32_t));
    file.write(reinterpret_cast<const char*>(&n_samples), sizeof(uint32_t));
    
    // 2. Escribir vector de Temperaturas
    file.write(reinterpret_cast<const char*>(temperatures.data()), n_T * sizeof(float));
    
    // 3. Escribir matriz de Magnetización (aplanada fila por fila)
    for (const auto& row : m_raw) {
        file.write(reinterpret_cast<const char*>(row.data()), row.size() * sizeof(float));
    }
    
    // 4. Escribir matriz de Energía (aplanada fila por fila)
    for (const auto& row : e_raw) {
        file.write(reinterpret_cast<const char*>(row.data()), row.size() * sizeof(float));
    }
    
    file.close();
    std::cout << "  ✓ (RAW) " << filename << " [" << n_samples << " muestras/T]" << std::endl;
}


void save_snapshot(const std::string& filename, const IsingModel& model, float temperature) {
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    auto config = model.get_configuration();
    auto shape = model.get_shape();
    
    uint32_t n_dims = shape.size();
    file.write(reinterpret_cast<char*>(&n_dims), sizeof(uint32_t));
    
    for (int dim : shape) {
        uint32_t d = dim;
        file.write(reinterpret_cast<char*>(&d), sizeof(uint32_t));
    }
    
    float T = temperature;
    file.write(reinterpret_cast<char*>(&T), sizeof(float));
    
    uint32_t n_occupied = config.size();
    file.write(reinterpret_cast<char*>(&n_occupied), sizeof(uint32_t));
    
    // Escribir valores de spins
    file.write(reinterpret_cast<const char*>(config.data()), n_occupied * sizeof(int8_t));
    
    // =====================================================
    // NUEVO: Escribir coordenadas de sitios ocupados
    // =====================================================
    if (model.get_dimension() <= 2) {
        // Para 1D y 2D: guardar pares (i, j)
        const auto& occupied_sites = model.get_occupied_2d();
        for (const auto& site : occupied_sites) {
            uint32_t i = site.first;
            uint32_t j = site.second;
            file.write(reinterpret_cast<char*>(&i), sizeof(uint32_t));
            file.write(reinterpret_cast<char*>(&j), sizeof(uint32_t));
        }
    } else {
        // Para 3D: guardar triplas (i, j, k)
        const auto& occupied_sites = model.get_occupied_3d();
        for (const auto& site : occupied_sites) {
            uint32_t i = site[0];
            uint32_t j = site[1];
            uint32_t k = site[2];
            file.write(reinterpret_cast<char*>(&i), sizeof(uint32_t));
            file.write(reinterpret_cast<char*>(&j), sizeof(uint32_t));
            file.write(reinterpret_cast<char*>(&k), sizeof(uint32_t));
        }
    }
    
    file.close();
}

// ====================================================================
// NUEVAS FUNCIONES PARA ANÁLISIS COMPLETO
// ====================================================================

// guarda vector<float> energies en binario: uint32_t n + n floats
void save_relaxation_data(const std::string& filename, const std::vector<float>& energies) {
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) throw std::runtime_error("Cannot open file: " + filename);
    uint32_t n = (uint32_t)energies.size();
    file.write(reinterpret_cast<const char*>(&n), sizeof(uint32_t));
    file.write(reinterpret_cast<const char*>(energies.data()), n * sizeof(float));
    file.close();
}

// Guarda datos de ajuste tanh para análisis cuantitativo
void save_tanh_fit_data(const std::string& filename, 
                       const std::vector<float>& H_over_T,
                       const std::vector<float>& m_experimental,
                       const std::vector<float>& m_fitted,
                       float a_param, float b_param, float r_squared) {
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) throw std::runtime_error("Cannot open file: " + filename);
    
    uint32_t n = H_over_T.size();
    file.write(reinterpret_cast<char*>(&n), sizeof(uint32_t));
    file.write(reinterpret_cast<const char*>(H_over_T.data()), n * sizeof(float));
    file.write(reinterpret_cast<const char*>(m_experimental.data()), n * sizeof(float));
    file.write(reinterpret_cast<const char*>(m_fitted.data()), n * sizeof(float));
    file.write(reinterpret_cast<const char*>(&a_param), sizeof(float));
    file.write(reinterpret_cast<const char*>(&b_param), sizeof(float));
    file.write(reinterpret_cast<const char*>(&r_squared), sizeof(float));
    
    file.close();
}

// Guarda datos de Tc para análisis de transición de fase
void save_tc_analysis_data(const std::string& filename,
                          const std::vector<std::string>& lattices,
                          const std::vector<float>& q_values,
                          const std::vector<std::vector<float>>& tc_matrix) {
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) throw std::runtime_error("Cannot open file: " + filename);
    
    uint32_t n_lattices = lattices.size();
    uint32_t n_q = q_values.size();
    
    file.write(reinterpret_cast<char*>(&n_lattices), sizeof(uint32_t));
    file.write(reinterpret_cast<char*>(&n_q), sizeof(uint32_t));
    
    // Guardar nombres de redes (cada uno 32 chars)
    for (const auto& lattice : lattices) {
        char lattice_arr[32] = {0};
        std::strncpy(lattice_arr, lattice.c_str(), 31);
        file.write(lattice_arr, 32);
    }
    
    // Guardar valores de q
    file.write(reinterpret_cast<const char*>(q_values.data()), n_q * sizeof(float));
    
    // Guardar matriz Tc[lattice][q]
    for (const auto& row : tc_matrix) {
        file.write(reinterpret_cast<const char*>(row.data()), n_q * sizeof(float));
    }
    
    file.close();
}

// ====================================================================
// NUEVA FUNCIÓN: Guardado RAW para Paramagnetismo (3D: T, H, Samples)
// ====================================================================
void save_paramagnetism_raw(const std::string& filename,
                            const std::vector<float>& temperatures,
                            const std::vector<float>& fields,
                            const std::vector<std::vector<std::vector<float>>>& raw_data) {
    // raw_data tiene dimensiones: [Indice_T][Indice_H][Indice_Muestra]
    
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    uint32_t n_T = (uint32_t)temperatures.size();
    uint32_t n_H = (uint32_t)fields.size();
    
    // Asumimos que todos tienen el mismo n_samples, tomamos el del primero
    uint32_t n_samples = 0;
    if (n_T > 0 && n_H > 0 && !raw_data[0].empty()) {
        n_samples = (uint32_t)raw_data[0][0].size();
    }
    
    // 1. Escribir dimensiones
    file.write(reinterpret_cast<const char*>(&n_T), sizeof(uint32_t));
    file.write(reinterpret_cast<const char*>(&n_H), sizeof(uint32_t));
    file.write(reinterpret_cast<const char*>(&n_samples), sizeof(uint32_t));
    
    // 2. Escribir ejes (Temperaturas y Campos)
    file.write(reinterpret_cast<const char*>(temperatures.data()), n_T * sizeof(float));
    file.write(reinterpret_cast<const char*>(fields.data()), n_H * sizeof(float));
    
    // 3. Escribir la matriz gigante 3D aplanada
    // Orden: Para cada T -> Para cada H -> Escribir vector de muestras
    for (size_t i = 0; i < n_T; ++i) {
        for (size_t j = 0; j < n_H; ++j) {
            const auto& samples = raw_data[i][j];
            if (samples.size() != n_samples) {
                std::cerr << "⚠️ Advertencia: Inconsistencia en tamaño de muestras en T=" << i << ", H=" << j << std::endl;
            }
            file.write(reinterpret_cast<const char*>(samples.data()), samples.size() * sizeof(float));
        }
    }
    
    file.close();
    std::cout << "  ✓ (RAW) " << filename << " [" << n_samples << " muestras/punto]" << std::endl;
}


void simulate_paramagnetism(const std::string& output_dir) {
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "PARAMAGNETISMO (J=0) - MODO RAW DATA" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    // --- PARÁMETROS ---
    int L = 8;              
    int n_samples = 10000;     // Aumentar a 1000+ para producción
    int thermal_steps = 0;  // No se necesita termalización para J=0
    // ------------------
    
    std::vector<float> q_values = {0.5f, 0.8f, 1.0f};
    std::vector<float> temperatures = {1.0f, 3.0f, 5.0f};
    std::vector<float> fields;
    
    for (int i = 0; i <= 40; ++i) {
        fields.push_back(-2.0f + i * 0.1f);
    }
    
    std::vector<std::string> lattices = {"chain", "honeycomb", "square", "bcc"};
    
    fs::create_directories(output_dir + "/paramagnetismo/");
    fs::create_directories(output_dir + "/relax/");
    
    for (const auto& lattice : lattices) {
        std::cout << "\n" << lattice << ":" << std::endl;
        
        // Nota: Para J=0, z no afecta la física (espines independientes)
        // pero lo incluimos por consistencia con otras simulaciones
        
        for (float q : q_values) {
            std::cout << "  q=" << std::fixed << std::setprecision(1) << q << " ... ";
            std::cout.flush();
            
            auto t_start = std::chrono::high_resolution_clock::now();
            
            IsingModel model(L, q, lattice, 0.0f); // J=0 para paramagnetismo
            
            // Estructura 3D: [Temp_idx][Field_idx][Sample_vector]
            std::vector<std::vector<std::vector<float>>> all_raw_data;
            
            for (float T : temperatures) {
                model.set_T(T);
                
                std::vector<std::vector<float>> raw_data_for_this_T;
                
                for (float H : fields) {
                    model.set_H(H);
                    model.initialize_random();
                    model.thermalize(thermal_steps);
                    
                    std::vector<float> m_history;
                    m_history.reserve(n_samples);
                    
                    for (int measure = 0; measure < n_samples; ++measure) {
                        model.metropolis_step();
                        m_history.push_back(model.magnetization());
                    }
                    
                    raw_data_for_this_T.push_back(m_history);
                }
                
                all_raw_data.push_back(raw_data_for_this_T);
            }
            
            std::string filename = output_dir + "/paramagnetismo/" + lattice + 
                                  "_q" + std::to_string((int)(q * 10)) + "0.bin";
            
            save_paramagnetism_raw(filename, temperatures, fields, all_raw_data);
            
            auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::high_resolution_clock::now() - t_start).count();
            std::cout << "(" << elapsed << " ms)" << std::endl;
        }
    }
}


// void simulate_ferromagnetism_hysteresis(const std::string& output_dir) {
//     std::cout << "\n" << std::string(70, '=') << std::endl;
//     std::cout << "FERROMAGNETISMO - HISTÉRESIS (J=1)" << std::endl;
//     std::cout << std::string(70, '=') << std::endl;
    
//     std::vector<float> q_values = {0.5f, 0.8f, 1.0f};
//     std::vector<float> temperatures = {2.0f, 5.0f};  // Dos valores de temperatura
//     float H_max = 3.0f;
//     int n_H_points = 20;
//     int steps_per_H = 100;
    
//     std::vector<std::string> lattices = {"chain", "honeycomb", "square", "bcc"};
    
//     fs::create_directories(output_dir + "/ferromagnetismo/");
//     fs::create_directories(output_dir + "/relax/");
    
//     for (const auto& lattice : lattices) {
//         std::cout << "\n" << lattice << ":" << std::endl;
        
//         // L variable según topología para mejor estadística
//         int L;
//         if (lattice == "chain") L = 500;
//         else if (lattice == "bcc") L = 8;
//         else L = 25;  // honeycomb, square
        
//         // Número de coordinación según topología
//         int z;
//         if (lattice == "chain") z = 2;
//         else if (lattice == "honeycomb") z = 3;
//         else if (lattice == "square") z = 4;
//         else if (lattice == "bcc") z = 8;
//         else z = 4;
        
//         std::cout << "  (L=" << L << ", z=" << z << ")" << std::endl;
        
//         // Escalar steps_per_H con z
//         int steps_per_H_scaled = steps_per_H * z / 2;
        
//         for (float q : q_values) {
//             for (float T : temperatures) {  // Nuevo loop para temperaturas
//                 std::cout << "  q=" << std::fixed << std::setprecision(1) << q 
//                          << " T=" << std::fixed << std::setprecision(1) << T << " ... ";
//                 std::cout.flush();
                
//                 auto t_start = std::chrono::high_resolution_clock::now();
                
//                 // ============================================
//                 // RELAJACIÓN INICIAL PARA HISTÉRESIS
//                 // ============================================
//                 IsingModel model_relax(L, q, lattice, 1.0f);
//                 model_relax.set_T(T);
//                 model_relax.set_H(H_max);
//                 model_relax.initialize_ordered(1);
                
//                 std::vector<float> E_t;
//                 for (int step = 0; step < 100; step++) {
//                     model_relax.metropolis_step();
//                     E_t.push_back(model_relax.energy());
//                 }
                
//                 std::string relax_fname = output_dir + "/relax/relax_" + lattice + 
//                                          "_q" + std::to_string((int)(q*100)) + 
//                                          "_T" + std::to_string((int)(T*10)) + "_J1_hyst.bin";
//                 save_relaxation_data(relax_fname, E_t);
                
//                 // ============================================
//                 // HISTÉRESIS PRINCIPAL
//                 // ============================================
//                 IsingModel model(L, q, lattice, 1.0f);
//                 model.set_T(T);
                
//                 // Preparación: saturar en H_max
//                 model.set_H(H_max);
//                 model.initialize_ordered(1);
//                 model.thermalize(500);
                
//                 // ══════════════════════════════════════════
//                 // RAMA DESCENDENTE: H_max → -H_max
//                 // ══════════════════════════════════════════
//                 std::vector<float> H_down, m_down;
                
//                 for (int i = 0; i < n_H_points; ++i) {
//                     float H = H_max - i * (2.0f * H_max) / (n_H_points - 1);
//                     model.set_H(H);
                    
//                     for (int step = 0; step < steps_per_H_scaled; step++) {
//                         model.metropolis_step();
//                     }
                    
//                     H_down.push_back(H);
//                     m_down.push_back(model.magnetization());
//                 }
                
//                 // ══════════════════════════════════════════
//                 // RAMA ASCENDENTE: -H_max → +H_max
//                 // ══════════════════════════════════════════
//                 std::vector<float> H_up, m_up;
                
//                 for (int i = 0; i < n_H_points; ++i) {
//                     float H = -H_max + i * (2.0f * H_max) / (n_H_points - 1);
//                     model.set_H(H);
                    
//                     for (int step = 0; step < steps_per_H_scaled; step++) {
//                         model.metropolis_step();
//                     }
                    
//                     H_up.push_back(H);
//                     m_up.push_back(model.magnetization());
//                 }
                
//                 // Modificar nombre del archivo para incluir temperatura
//                 std::string filename = output_dir + "/ferromagnetismo/hysteresis_" + lattice + 
//                                       "_q" + std::to_string((int)(q * 100)) + 
//                                       "_T" + std::to_string((int)(T * 10)) + ".bin";
                
//                 save_hysteresis_data(filename, H_down, m_down, H_up, m_up);
                
//                 auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
//                     std::chrono::high_resolution_clock::now() - t_start).count();
//                 std::cout << "(" << elapsed << " ms)" << std::endl;
//             }
//         }
//     }
// }

void simulate_ferromagnetism_hysteresis(const std::string& output_dir) {
    using clock = std::chrono::high_resolution_clock;
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "FERROMAGNETISMO - HISTÉRESIS (OPTIMIZADO, J=1)" << std::endl;
    std::cout << std::string(70, '=') << std::endl;

    // ---------- parámetros (ajustables) ----------
    const float H_max = 3.0f;
    const int n_H_points = 80;        // resolución razonable y rápida
    const int equil_sweeps = 30;      // sweeps para equilibrar por punto H
    const int meas_sweeps = 10;       // sweeps para medir y promediar por punto H
    const int n_realizations = 3;     // promedios con semillas independientes
    const int saturate_sweeps = 5;    // sweeps para saturación inicial (±Hmax)
    // ---------------------------------------------

    std::vector<float> q_values = {0.5f, 0.8f, 1.0f};
    std::vector<float> temperatures = {2.0f, 5.0f};
    std::vector<std::string> lattices = {"chain", "honeycomb", "square", "bcc"};

    fs::create_directories(output_dir + "/ferromagnetismo/");

    for (const auto& lattice : lattices) {
        int L;
        if (lattice == "chain") L = 500;
        else if (lattice == "bcc") L = 8;
        else L = 25; // honeycomb, square

        std::cout << "\n" << lattice << " (L=" << L << ") :" << std::endl;

        for (float q : q_values) {
            for (float T : temperatures) {
                std::cout << "  q=" << q << "  T=" << T << " ... " << std::flush;
                auto t_case_start = clock::now();

                // construir malla H (descendente: +Hmax -> -Hmax)
                std::vector<float> H_grid(n_H_points);
                for (int i = 0; i < n_H_points; ++i)
                    H_grid[i] = H_max - i * (2.0f * H_max) / (n_H_points - 1);

                std::vector<float> m_down(n_H_points, 0.0f);
                std::vector<float> m_up(n_H_points, 0.0f);

                // realizaciones independientes
                for (int run = 0; run < n_realizations; ++run) {
                    uint32_t seed = static_cast<uint32_t>(std::random_device{}());

                    // Creamos modelo (uso mismo API que antes)
                    IsingModel model(L, q, lattice, 1.0f, seed);
                    model.set_T(T);

                    // número de intentos por sweep (aprox. N ocupado)
                    // int n_attempts = std::max(1, (int)model.get_n_occupied()); // ← COMENTADO (no se usa)

                    // --- saturar en +Hmax ---
                    model.set_H(H_max);
                    model.initialize_ordered(1);
                    for (int s = 0; s < saturate_sweeps; ++s) model.metropolis_step();

                    // --- rama descendente: +Hmax -> -Hmax ---
                    for (int ih = 0; ih < n_H_points; ++ih) {
                        model.set_H(H_grid[ih]);

                        // equilibrar: equil_sweeps sweeps
                        for (int s = 0; s < equil_sweeps; ++s) model.metropolis_step();

                        // medir: meas_sweeps sweeps (promedio)
                        double m_acc = 0.0;
                        for (int s = 0; s < meas_sweeps; ++s) {
                            model.metropolis_step();
                            m_acc += model.magnetization();
                        }
                        m_down[ih] += static_cast<float>(m_acc / meas_sweeps);
                    }

                    // --- saturar en -Hmax (iniciamos una configuración saturada negativa) ---
                    model.set_H(-H_max);
                    model.initialize_ordered(-1);
                    for (int s = 0; s < saturate_sweeps; ++s) model.metropolis_step();

                    // --- rama ascendente: -Hmax -> +Hmax ---
                    for (int ih = 0; ih < n_H_points; ++ih) {
                        // recorrer H_grid en orden inverso (desde -Hmax hacia +Hmax)
                        float Hval = H_grid[n_H_points - 1 - ih];
                        model.set_H(Hval);

                        for (int s = 0; s < equil_sweeps; ++s) model.metropolis_step();

                        double m_acc = 0.0;
                        for (int s = 0; s < meas_sweeps; ++s) {
                            model.metropolis_step();
                            m_acc += model.magnetization();
                        }
                        // almacenar en la posición correspondiente (coincide con H_grid index)
                        m_up[n_H_points - 1 - ih] += static_cast<float>(m_acc / meas_sweeps);
                    }
                } // end runs

                // promedio sobre realizaciones
                for (int i = 0; i < n_H_points; ++i) {
                    m_down[i] /= static_cast<float>(n_realizations);
                    m_up[i]   /= static_cast<float>(n_realizations);
                }

                // Guardar (compatible con reader)
                std::string file = output_dir +
                    "/ferromagnetismo/hysteresis_" + lattice +
                    "_q" + std::to_string(int(q*100)) +
                    "_T" + std::to_string(int(T*10)) +
                    "_corrected.bin";

                save_hysteresis_data(file, H_grid, m_down, H_grid, m_up);

                auto t_case_end = clock::now();
                auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_case_end - t_case_start).count();
                std::cout << "  ✓  (≈ " << ms << " ms)" << std::endl;
            }
        }
    }
}

void simulate_phase_transition(const std::string& output_dir) {
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "TRANSICIÓN DE FASE (J=1) - MODO RAW DATA" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    int n_samples = 10000; 
    std::vector<float> q_values = {0.5f, 0.8f, 1.0f};
    std::vector<float> temperatures;
    for (float T = 0.5f; T <= 9.0f; T += 0.1f) {
        temperatures.push_back(T);
    }
    
    std::vector<std::string> lattices = {"chain", "honeycomb", "square", "bcc"};
    
    fs::create_directories(output_dir + "/transicion/");
    fs::create_directories(output_dir + "/relax/");
    fs::create_directories(output_dir + "/temperatura_critica/");
    
    std::vector<std::vector<float>> tc_matrix(4, std::vector<float>(3, 0.0f));
    
    for (int lat_idx = 0; lat_idx < lattices.size(); lat_idx++) {
        const auto& lattice = lattices[lat_idx];
        
        // L variable según topología
        int L;
        if (lattice == "chain") L = 500;
        else if (lattice == "bcc") L = 8;
        else L = 25;  // honeycomb, square
        
        // Número de coordinación según topología
        int z;
        if (lattice == "chain") z = 2;
        else if (lattice == "honeycomb") z = 3;
        else if (lattice == "square") z = 4;
        else if (lattice == "bcc") z = 8;
        else z = 4;
        
        // Termalización aumentada para BCC
        int therm_steps;
        if (lattice == "bcc") {
            therm_steps = 250 * z * 16;  // 32000 pasos
        } else {
            therm_steps = 250 * z;
        }

        std::cout << "\n" << lattice << " (L=" << L << ", z=" << z << ", therm=" << therm_steps << "):" << std::endl;
        
        for (int q_idx = 0; q_idx < q_values.size(); q_idx++) {
            float q = q_values[q_idx];
            std::cout << "  q=" << q << " ... ";
            std::cout.flush();
            
            auto t_start = std::chrono::high_resolution_clock::now();
            
            IsingModel model(L, q, lattice, 1.0f);
            model.set_H(0.0f);
            
            std::vector<std::vector<float>> m_raw_data;
            std::vector<std::vector<float>> e_raw_data;
            std::vector<float> m_avg_temp; 
            
            for (float T : temperatures) {
                model.set_T(T);
                model.initialize_random();
                model.thermalize(therm_steps);
                
                std::vector<float> m_history;
                std::vector<float> e_history;
                m_history.reserve(n_samples);
                e_history.reserve(n_samples);
                
                double m_sum_local = 0.0;
                
                for (int i = 0; i < n_samples; i++) {
                    model.metropolis_step();
                    
                    float m_val = std::abs(model.magnetization());
                    float e_val = model.energy();
                    
                    m_history.push_back(m_val);
                    e_history.push_back(e_val);
                    m_sum_local += m_val;
                }
                
                m_raw_data.push_back(m_history);
                e_raw_data.push_back(e_history);
                m_avg_temp.push_back(m_sum_local / n_samples);
            }
            
            // Estimar Tc
            float Tc_estimated = 0.0f;
            float max_slope = 0.0f;
            
            for (size_t i = 1; i < m_avg_temp.size() - 1; i++) {
                float slope = std::abs(m_avg_temp[i+1] - m_avg_temp[i-1]) / 
                             (temperatures[i+1] - temperatures[i-1]);
                if (slope > max_slope) {
                    max_slope = slope;
                    Tc_estimated = temperatures[i];
                }
            }
            
            tc_matrix[lat_idx][q_idx] = Tc_estimated;
            std::cout << "Tc≈" << std::fixed << std::setprecision(2) << Tc_estimated << " ";
            
            std::string filename = output_dir + "/transicion/transicion_" + lattice + ".bin";
            save_phase_transition_raw(filename, temperatures, m_raw_data, e_raw_data);
            
            auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::high_resolution_clock::now() - t_start).count();
            std::cout << "(" << elapsed << " ms)" << std::endl;
        }
    }
    
    std::string tc_filename = output_dir + "/temperatura_critica/tc_analysis.bin";
    save_tc_analysis_data(tc_filename, lattices, q_values, tc_matrix);
    std::cout << "\n✓ Análisis de Tc guardado." << std::endl;
}
// ====================================================================
// Relajación para múltiples T y J (compatible con Python)
// ====================================================================
void simulate_relaxation(const std::string& output_dir) {
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "RELAJACIÓN DINÁMICA (Solo Ferromagnético J=1)" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    // NOTA: Eliminamos 'int L = 8' de aquí porque ahora es variable
    
    std::vector<float> q_values = {0.5f, 0.8f, 1.0f};  
    std::vector<float> temperatures = {1.0f, 3.0f, 5.0f};
    
    // CAMBIO 1: Solo simulamos J=1 (ahorramos 50% de tiempo)
    std::vector<int> J_values = {1}; 
    
    int n_steps = 10000;  
    
    std::vector<std::string> lattices = {"chain", "honeycomb", "square", "bcc"};
    
    fs::create_directories(output_dir + "/relax/");
    
    for (const auto& lattice : lattices) {
        
        // CAMBIO 2: Asignación dinámica de L para mantener N ~ 500-600
        int L;
        if (lattice == "bcc") {
            L = 8;          // 3D: 8^3 = 512 espines
        } else if (lattice == "chain") {
            L = 500;        // 1D: 500 espines (Necesario para ver algo en 1D)
        } else {
            L = 25;         // 2D (Square, Honeycomb): 25^2 = 625 espines
        }
        
        std::cout << lattice << " (L=" << L << "): ";
        
        for (float q : q_values) {
            for (int J : J_values) {
                for (float T : temperatures) {
                    
                    // Aquí se crea el modelo con el L específico para esta red
                    IsingModel model(L, q, lattice, (float)J);
                    model.set_T(T);
                    model.set_H(0.0f);
                    model.initialize_random(); // Arranque en caliente (T infinita)
                    
                    // Vector para guardar energía
                    std::vector<float> E_t;
                    E_t.reserve(n_steps);
                    
                    for (int step = 0; step < n_steps; step++) {
                        model.metropolis_step();
                        

                        E_t.push_back(model.energy());
                    }
                
                    // Formato de nombre coincidente con tu Python
                    std::string q_str = std::to_string((int)(q * 100)); 
                    std::string filename = output_dir + "/relax/relax_" + lattice + 
                                          "_q" + q_str +
                                          "_T" + std::to_string((int)(T * 10)) +
                                          "_J" + std::to_string(J) + ".bin";
                    
                    save_relaxation_data(filename, E_t);
                    
                    // Feedback visual corto en consola
                    std::cout << "[q" << q_str << " T" << (int)(T*10) << "] ";
                }
            }
        }
        std::cout << std::endl;
    }
}

void simulate_snapshots(const std::string& output_dir) {
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "SNAPSHOTS" << std::endl;
    std::cout << std::string(70, '=') << std::endl;

    std::vector<float> q_values = {0.5f, 0.8f, 1.0f};
    std::vector<float> T_values = {1.0f, 3.0f, 5.0f};
    std::vector<std::string> lattices = {"chain", "honeycomb", "square", "bcc"};
    
    fs::create_directories(output_dir + "/snapshots/");
    
    for (const auto& lattice : lattices) {
        std::cout << lattice << ": ";
        
        // L variable según topología
        int L;
        if (lattice == "chain") L = 500;
        else if (lattice == "bcc") L = 8;
        else L = 25;  // honeycomb, square
        
        // Número de coordinación según topología
        int z;
        if (lattice == "chain") z = 2;
        else if (lattice == "honeycomb") z = 3;
        else if (lattice == "square") z = 4;
        else if (lattice == "bcc") z = 8;
        else z = 4;
        
        // Termalización escalada con aumento para BCC
        int therm_steps;
        if (lattice == "bcc") {
            therm_steps = L * L * L * z * 4;  // 16384 pasos
        } else if (lattice == "chain") {
            therm_steps = L * z;
        } else {
            therm_steps = L * L * z;
        }
        
        for (float q : q_values) {
            for (float T : T_values) {
                IsingModel model(L, q, lattice, 1.0f);
                model.set_T(T);
                model.initialize_random();
                model.thermalize(therm_steps);
                
                // Formato q: q*10 con 0 al final → q50, q80, q100
                std::string q_str = std::to_string((int)(q * 10)) + "0";
                
                std::string filename = output_dir + "/snapshots/snapshot_" + lattice + 
                                      "_q" + q_str +
                                      "_T" + std::to_string((int)T) + ".bin";
                
                save_snapshot(filename, model, T);
                std::cout << "q" << q << "T" << (int)T << " ";
            }
        }
        std::cout << std::endl;
    }
}


// ====================================================================
// MAIN
// ====================================================================
int main() {
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "SIMULACIÓN MONTE CARLO - MODELO DE ISING (C++)" << std::endl;
    std::cout << "Universidad de Antioquia - Física Estadística" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    std::string output_dir = "data";
    fs::create_directories(output_dir);
    
    // Información de configuración
    std::cout << "\nConfiguración:" << std::endl;
    std::cout << "  Directorio de salida: " << output_dir << "/" << std::endl;
    std::cout << "  Redes: chain, honeycomb, square, bcc" << std::endl;
    std::cout << "  Valores de q: 0.5, 0.8, 1.0" << std::endl;
    
    auto t_total_start = std::chrono::high_resolution_clock::now();
    
    try {
        // Orden: de más rápido a más lento
        //simulate_snapshots(output_dir);
        //simulate_relaxation(output_dir);
        //simulate_paramagnetism(output_dir);
        simulate_ferromagnetism_hysteresis(output_dir);
        //simulate_phase_transition(output_dir);  // 
        
        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
            std::chrono::high_resolution_clock::now() - t_total_start).count();
        
        std::cout << "\n" << std::string(70, '=') << std::endl;
        std::cout << "✅ SIMULACIONES COMPLETADAS" << std::endl;
        std::cout << "Tiempo total: " << elapsed << " segundos" << std::endl;
        std::cout << "Datos guardados en: " << output_dir << "/" << std::endl;
        std::cout << std::string(70, '=') << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "❌ ERROR: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}