# generate_header.py

import numpy as np

# Nazwa pliku wejściowego i wyjściowego
input_file = "./matlab_sim/ura_rx.txt"
output_file = "./stm_code/asl_demo/Core/Inc/test_data.h"

# Wczytanie danych z pliku txt
# Zakładamy, że dane są oddzielone spacjami lub tabulatorami
data = np.loadtxt(input_file, delimiter=',')

# Sprawdzenie wymiarów
# if data.shape != (1024, 16):
#     raise ValueError(f"Oczekiwano wymiarów 1024x16, a otrzymano {data.shape}")

# Otwieramy plik .h do zapisu
with open(output_file, "w") as f:
    f.write("float mic_data_buff[16][1024] = {\n")
    
    # Transponujemy dane, żeby odpowiadały formatowi [16][1024]
    data = data.T

    for i, row in enumerate(data):
        f.write("    {")
        # zapisujemy wartości jako float z literą 'f'
        f.write(", ".join(f"{val:.6f}f" for val in row))
        f.write("}")
        if i != len(data) - 1:
            f.write(",\n")
        else:
            f.write("\n")
    f.write("};\n")

#print(f"Plik {output_file} został wygenerowany pomyślnie.")
