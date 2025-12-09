import numpy as np

def generate_h_file(samples_txt, positions_txt, output_h, 
                    samples_array_name="float_mic_data_buff",
                    positions_array_name="mics_positions"):
    # Wczytaj dane próbek
    samples = np.loadtxt(samples_txt, delimiter=',')
    M, L = samples.shape
    print(f"Wczytano próbki o wymiarach: {M}x{L}")
    
    # Wczytaj pozycje mikrofonów
    positions = np.loadtxt(positions_txt, delimiter=',')
    coords, num_mics = positions.shape
    print(f"Wczytano pozycje mikrofonów o wymiarach: {coords}x{num_mics}")
    
    # Sprawdź zgodność liczby mikrofonów
    if M != num_mics:
        print(f"UWAGA: Liczba mikrofonów w próbkach ({M}) różni się od liczby w pozycjach ({num_mics})")
    
    # Generuj plik .h
    with open(output_h, 'w') as f:
        # Nagłówek
        f.write("#ifndef TEST_DATA_H\n")
        f.write("#define TEST_DATA_H\n\n")
        f.write("#include \"arm_math.h\"\n\n")
        
        # Deklaracja tablicy próbek
        f.write(f"const float32_t {samples_array_name}[{M}][{L}] = {{\n")
        
        # Zapisz próbki
        for i in range(M):
            f.write("    {")
            for j in range(L):
                f.write(f"{samples[i, j]:.8e}f")
                if j < L - 1:
                    f.write(", ")
                # Nowa linia co 4 elementy dla czytelności
                if (j + 1) % 4 == 0 and j < L - 1:
                    f.write("\n     ")
            f.write("}")
            if i < M - 1:
                f.write(",\n")
            else:
                f.write("\n")
        
        f.write("};\n\n")
        
        # Definicja struktury coords_3D_t
        f.write("typedef struct coords_3D_t {\n")
        f.write("    float32_t x;\n")
        f.write("    float32_t y;\n")
        f.write("    float32_t z;\n")
        f.write("} coords_3D_t;\n\n")
        
        # Deklaracja tablicy pozycji mikrofonów
        f.write(f"const coords_3D_t {positions_array_name}[{num_mics}] = {{\n")
        
        # Zapisz pozycje - transpozycja, żeby iterować po mikrofonach
        for j in range(num_mics):
            f.write("    {")
            f.write(f"{positions[0, j]:.8e}f, {positions[1, j]:.8e}f, {positions[2, j]:.8e}f")
            f.write("}")
            if j < num_mics - 1:
                f.write(",\n")
            else:
                f.write("\n")
        
        f.write("};\n\n")
        f.write("#endif // TEST_DATA_H\n")
    
    print(f"Plik {output_h} został wygenerowany pomyślnie!")

# Przykład użycia
if __name__ == "__main__":
    samples_file = "./matlab_sim/ura_rx.txt"  # Plik z próbkami 16x512
    positions_file = "./matlab_sim/mics.txt"  # Plik z pozycjami 3x16
    output_file = "./stm_code/asl_demo/Core/Inc/test_data.h"
    
    generate_h_file(samples_file, positions_file, output_file)