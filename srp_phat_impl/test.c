/* === prolly some meta info here - TODO === */


/* === includes === */

#include <stdio.h>
// #include <windows.h>
// #include <psapi.h>
#include <stdint.h>
#include <math.h>
#include <limits.h>

#include "kiss_fftr.h"


/* === typedefs === */

typedef struct coords_3D_t {
    float x;
    float y;
    float z;
} coords_3D_t;

typedef struct angle_pair_t {
    float az;
    float el;
} angle_pair_t;

/**
 * @brief this struct holds TDoA between reference mic (0th idx) and assigned mic (mix_idx);
 *        such implementation helps reducing memory cost of program
 * 
 * @param tau calculated relative TDoA
 * @param mic_idx microphone index
 * @param angle_pair_idx indicates for which DoA angles pair {az, el} indeex in *all_possible_doas_lut*
 *                       this tau is relevant
 */
typedef struct tdoa_per_mic_per_doa_t {
    float tau;
    uint8_t mic_idx;
    uint16_t angle_pair_idx;
} tdoa_per_mic_per_doa_t;

/**
 * @brief enumerator used to indicate whether data acquisition from SPI was successful or not
 */
typedef enum spi_mic_acq_enum_t {
    MIC_ACQ_OK = 0,
    MIC_ACQ_ERR = 1
} spi_mic_acq_enum_t;


/* === defines === */

#define M 16                    /* num of mics */
#define Ny 4                    /* num of mics in y column */
#define Nz 4                    /* num of mics in z row */
#define C 343                   /* speed of sound in m/s */
#define d 0.042                 /* spacing between mics */
#define F_MIN 300               /* minimum freq considered in Hz */
#define F_MAX 3000              /* minimum freq considered in Hz */
#define F_S 16000               /* incoming signal sampling rate */
#define L 1024                  /* x_m microphone's signal frame length - TO BE DETERMINED */
#define R 10                    /* angualar grid resolution in degrees */
#define AZ_LOW -90              /* min azimuth angle conidered in degrees */
#define AZ_HIGH 90              /* max azimuth angle conidered in degrees */
#define EL_LOW -90              /* min elevation angle considered in degrees */
#define EL_HIGH 90              /* max elevation angle considered in degrees */
#define PI 3.141592653589793238 /* pi approximation */
#define EPS 0.0000000001        /* epsilon value for numerical stability */

#define F_NYQ (F_S / 2)                                     /* Nyquist frequency */
#define REAL_FFT_LEN ((L/2) + 1)                            /* length of the output of real-valued signal's FFT, which is of use */
#define ALL_AZIMUTHS ((AZ_HIGH - AZ_LOW)/R + 1)             /* size of azimuth grid */
#define ALL_ELEVATIONS ((EL_HIGH - EL_LOW)/R + 1)           /* size of elevation grid */
#define ALL_POSSIBLE_DOAS (ALL_AZIMUTHS * ALL_ELEVATIONS)   /* size of all DoAs grid (i.e. every possible pair of {az,el}) */


/* === constant globals === */

// const coords_3D_t mics_positions[M] = {
//     { 0,   0, 0 }, { 0,   0, d }, { 0,   0, 2*d }, { 0,   0, 3*d },
//     { 0,   d, 0 }, { 0,   d, d }, { 0,   d, 2*d }, { 0,   d, 3*d },
//     { 0, 2*d, 0 }, { 0, 2*d, d }, { 0, 2*d, 2*d }, { 0, 2*d, 3*d },
//     { 0, 3*d, 0 }, { 0, 3*d, d }, { 0, 3*d, 2*d }, { 0, 3*d, 3*d }
// }; //<< 192 B


/* === non-constant globals === */
coords_3D_t mics_positions[M];
/**
 * @brief a look-up table storing every possible DoA as pair {az, el}
 */
static angle_pair_t all_possible_doas_lut[ALL_POSSIBLE_DOAS]; //<< 2,888 B
/**
 * @brief a look-up table storing every mic's relative TDoA for every possible DoA
 *        in relation to 0th microphone - this way we save memory in trade for slightly 
 *        more calculations in the run-time
 */
static tdoa_per_mic_per_doa_t all_relative_tdoas_lut[ALL_POSSIBLE_DOAS * M]; //<< 46,208 B with padding, 40,432 B optimized

/**
 * @brief a look-up table storing every frequency bin of interest
 */
static float fft_bins_lut[REAL_FFT_LEN];//<< 2,052 B
/**
 * @brief indicators of minimum and maximum bins of FFT for freqs of interest
 */
static uint16_t fft_bin_low = 0, fft_bin_high = 0;//<< 2 * 2 B = 4 B

/**
 * @brief kiss_fft library config field, used to run FFT - this is a real-input FFT, outputs cpx
 */
kiss_fftr_cfg kissfftrcfg = NULL; //<< this is a ptr, idk how much it allocs

/**
 * @brief preallocated buffer for incoming real signal data from mics (thru SPI)
 */
float mic_data_buff[M][L] = { 0 };//<<  65,536 B

/**
 * @brief preallocated buffer for incoming real signal data FFT 
 */
static kiss_fft_cpx mic_data_fft_buff[M][REAL_FFT_LEN]; //<< 65,664 B

/**
 * @brief preallocated buffer for SRP-PHAT map (as function of DoA angles pair)
 * 
 */
static float srp_phat_map_per_doa[ALL_POSSIBLE_DOAS];//<< 1,444 B


/* === functions === */

/**
 * @brief basic impl of conversion from degrees to radians
 * 
 * @param arg_deg argument in degrees 
 * @return float argument converted to radians
 */
static inline float deg2rad(float arg_deg) {
    return (arg_deg * (float)(PI / 180.0f));
}

/**
 * @brief basic impl of conversion from radians to degrees
 * 
 * @param arg_rad argument in radians
 * @return float argument converted to degrees
 */
static inline float rad2deg(float arg_rad) {
    return (arg_rad * (float)(180.0f / PI));
}

/**
 * @brief this function prepares look-up tables: *all_possible_doas_lut* and *all_relative_tdoas_lut*
 */
static void prepare_luts(void) {
    {// prep mics
        FILE *f = fopen("../mics.txt", "r");
        if (!f) {
            perror("plik z majkami nie dzialaaaa\n");
            return;
        }

        char line[2048];
        int coord = 0;

        while (coord < 3 && fgets(line, sizeof(line), f)) {
            char *ptr = line;

            for (int mic = 0; mic < 16; mic++) {
                if (coord == 0) {
                    mics_positions[mic].x = strtof(ptr, &ptr);
                }
                else if (coord == 1) {
                    mics_positions[mic].y = strtof(ptr, &ptr);
                }
                else {
                    mics_positions[mic].z = strtof(ptr, &ptr);
                }

                if (*ptr == ',') ptr++;
            }

            coord++;
        }

        fclose(f);
    }

    coords_3D_t sphere_vec;
    uint16_t doa_idx = 0;
    uint32_t tdoa_idx = 0;
    uint16_t el_idx = 0;
    uint8_t mic_idx = 0;

    /* for every azimuth degree step */
    for (uint16_t az_idx = 0; az_idx < ALL_AZIMUTHS; az_idx++) {
        /* and for every elevation degree step */
        for (el_idx = 0; el_idx < ALL_ELEVATIONS; el_idx++) {
            /* define an entry in degree-grid */
            all_possible_doas_lut[doa_idx].az = (float) (AZ_LOW + (az_idx * R));
            all_possible_doas_lut[doa_idx].el = (float) (EL_LOW + (el_idx * R));
            
            /* Convert to radians ONCE for this DoA */
            float az_rad = deg2rad(all_possible_doas_lut[doa_idx].az);
            float el_rad = deg2rad(all_possible_doas_lut[doa_idx].el);
            
            /* 0th element is our ref, so it's zero all the time */
            all_relative_tdoas_lut[tdoa_idx].angle_pair_idx = doa_idx;
            all_relative_tdoas_lut[tdoa_idx].mic_idx = (uint8_t) 0;
            all_relative_tdoas_lut[tdoa_idx].tau = (float) 0;
            tdoa_idx++;

            /* Calculate spherical direction vector according to MATLAB convention:
             * d_hat(1) = cos(theta) * cos(phi)   where theta=elevation, phi=azimuth
             * d_hat(2) = cos(theta) * sin(phi)
             * d_hat(3) = sin(theta)
             */
            sphere_vec.x = cosf(el_rad) * cosf(az_rad);
            sphere_vec.y = cosf(el_rad) * sinf(az_rad);
            sphere_vec.z = sinf(el_rad);

            /* every other mic's relative TDoA is calculated relatively to 0th mic */
            for (mic_idx = 1; mic_idx < M; mic_idx++) {
                all_relative_tdoas_lut[tdoa_idx].angle_pair_idx = doa_idx;
                all_relative_tdoas_lut[tdoa_idx].mic_idx = mic_idx;
                all_relative_tdoas_lut[tdoa_idx].tau = (float) (
                    (sphere_vec.x * (mics_positions[mic_idx].x - mics_positions[0].x)) 
                  + (sphere_vec.y * (mics_positions[mic_idx].y - mics_positions[0].y))
                  + (sphere_vec.z * (mics_positions[mic_idx].z - mics_positions[0].z))
                ) / ((float) C);
                tdoa_idx++;
            }
            
            doa_idx++;
        }
    }
    
    /* prepare every FFT bin and indicators */
    fft_bin_low = 0;
    fft_bin_high = 0;
    
    for (uint16_t k = 0; k < REAL_FFT_LEN; k++) {
        fft_bins_lut[k] = k * ((float)F_S / (float)L);

        // Count bins below F_MIN
        if (fft_bins_lut[k] < F_MIN) {
            fft_bin_low++;
        }

        // Count bins up to AND including F_MAX
        if (fft_bins_lut[k] <= F_MAX) {
            fft_bin_high = k + 1;  // Set to NEXT index (for < comparison in loop)
        }
    }
    
    // Debug output
    // printf("FFT bins: low=%d (%.1f Hz), high=%d (%.1f Hz)\n", 
    //        fft_bin_low, fft_bins_lut[fft_bin_low],
    //        fft_bin_high-1, fft_bins_lut[fft_bin_high-1]);
}

/**
 * @brief this functions acquires data from microphones using SPI and converts collected signals 
 *        to kissfft_cpx array format
 * 
 * @return spi_mic_acq_enum_t enum indication of errors or success
 */
spi_mic_acq_enum_t acquire_mic_data(void) {
    /* static array used to write data from microphones */
    spi_mic_acq_enum_t ret_code = MIC_ACQ_ERR;

    /* 
        fill with data from SPI
    */
    // for now, read simulated data
    FILE *f = fopen("../ura_rx.txt", "r");
    
    if (!f) {
        perror("error przy otwarciu pliku ;(");
    }
    else {
        char line[2048];   // linia z 16 floatami oddzielonymi przecinkami
        int frame = 0;

        while (frame < L && fgets(line, sizeof(line), f)) {

            char *ptr = line;
            for (int mic = 0; mic < M; mic++) {
                mic_data_buff[mic][frame] = strtof(ptr, &ptr);  // konwersja float
                if (*ptr == ',') ptr++;              // pomiń przecinek
            }

            frame++;
        }

        fclose(f);
        
        ret_code = MIC_ACQ_OK; // <-- this shall be asigned accordingly to situation
    }

    return ret_code;
}

/**
 * @brief helper function to run FFT for real-valued signals acquired from mics
 * 
 */
static void run_rfft(void) {
    for (uint8_t m = 0; m < M; m++) {
        kiss_fftr(kissfftrcfg, &mic_data_buff[m][0], &mic_data_fft_buff[m][0]);
    }
}

/**
 * @brief absolute value of complex number in kiss_fft_cpx convention
 * 
 * @param c complex to compute abs from
 * @return kiss_fft_scalar real-valued abs(c)
 */
static kiss_fft_scalar cpx_abs(const kiss_fft_cpx c) {
    return sqrtf((c.r * c.r) + (c.i * c.i));
} 

/**
 * @brief complex conjugation in kiss_fft_cpx convention
 * 
 * @param c complex to compute conj from
 * @return kiss_fft_cpx complex-valued conj(c)
 */
static kiss_fft_cpx cpx_conj(const kiss_fft_cpx c) {
    kiss_fft_cpx ret = {
        .r = (kiss_fft_scalar) c.r,
        .i = (kiss_fft_scalar)((-1) * c.i)
    };
    return ret;
} 

/**
 * @brief complex multiplication of two complex numbers in kiss_fft_cpx convention
 * 
 * @param a first complex
 * @param b second complex
 * @return kiss_fft_cpx complex-valued result of a*b = (a.r + a.i) * (b.r + b.i) = ...
 */
static kiss_fft_cpx cpx_mult(const kiss_fft_cpx a, const kiss_fft_cpx b) {
    kiss_fft_cpx ret = {
        .r = (kiss_fft_scalar) ((a.r * b.r) - (a.i * b.i)),
        .i = (kiss_fft_scalar) ((a.r * b.i) + (a.i * b.r))
    };
    return ret;
}

/**
 * @brief implements the proper SRP-PHAT algorithm search and fills srp_phat_map_per_doa array 
 *        with computed values - after this call, the array maximum's index is the same
 *        as estimated DoA's in all_possible_doas_lut array
 */
void srp_phat(void) {
    kiss_fft_cpx gcc_phat_sum;
    kiss_fft_cpx tmp_cross_spec;
    float tmp_cross_spec_abs = 0;
    float tmp_theta = 0;
    float cos_theta = 0;
    float sin_theta = 0;
    float tau_ml = 0;
    uint16_t doa_idx = 0;
    
    // Initialize map
    for (uint16_t i = 0; i < ALL_POSSIBLE_DOAS; i++) {
        srp_phat_map_per_doa[i] = 0.0f;
    }
    
    /* For each DoA... */
    for (uint16_t doa = 0; doa < ALL_POSSIBLE_DOAS; doa++) {
        // Calculate base index for this DoA in the TDoA LUT
        uint32_t doa_base_idx = (uint32_t)doa * (uint32_t)M;
        
        /* For each unique microphone pair (m, l) where l > m */
        for (uint8_t mic_m = 0; mic_m < M; mic_m++) {
            for (uint8_t mic_l = mic_m + 1; mic_l < M; mic_l++) {
                
                gcc_phat_sum.r = 0.0f;
                gcc_phat_sum.i = 0.0f;
                
                /* TDoA for this particular pair and DoA */
                tau_ml = all_relative_tdoas_lut[doa_base_idx + mic_l].tau 
                       - all_relative_tdoas_lut[doa_base_idx + mic_m].tau;
                
                /* For each frequency bin in [F_MIN, F_MAX] */
                for (uint16_t k = fft_bin_low; k < fft_bin_high; k++) {
                    // Cross-spectrum
                    tmp_cross_spec = cpx_mult(mic_data_fft_buff[mic_l][k], 
                                              cpx_conj(mic_data_fft_buff[mic_m][k]));
                    tmp_cross_spec_abs = cpx_abs(tmp_cross_spec);
                    
                    /* Enter if magnitude is non-zero */
                    if (tmp_cross_spec_abs > EPS) {
                        /* GCC-PHAT: normalize by magnitude */
                        tmp_cross_spec.r /= tmp_cross_spec_abs;
                        tmp_cross_spec.i /= tmp_cross_spec_abs;
                        
                        /* Steering vector: e^(j * 2 * pi * f * tau) */
                        tmp_theta = 2.0f * (float)PI * fft_bins_lut[k] * tau_ml;
                        cos_theta = cosf(tmp_theta);
                        sin_theta = sinf(tmp_theta);
                        
                        /* Multiply GCC-PHAT with steering vector */
                        gcc_phat_sum.r += tmp_cross_spec.r * cos_theta - tmp_cross_spec.i * sin_theta;
                        gcc_phat_sum.i += tmp_cross_spec.r * sin_theta + tmp_cross_spec.i * cos_theta;
                    }
                }
                
                /* Add power to SRP map for this DoA */
                srp_phat_map_per_doa[doa] += (gcc_phat_sum.r * gcc_phat_sum.r) 
                                            + (gcc_phat_sum.i * gcc_phat_sum.i);
            }
        }
    }
}

/**
 * @brief this func returns index of the found maximum in already-computed SRP-PHAT map
 * 
 * @return uint32_t index of the found maximum
 */
uint32_t find_max_doa_idx() {
    uint32_t ret = 0;
    for (uint32_t i = 0; i < ALL_POSSIBLE_DOAS; i++) {
        // if (srp_phat_map_per_doa[i] > 2000000) {
        //     printf("az = %+.2f,\tel = %+.2f,\tval = %10.0f\n", all_possible_doas_lut[i].az, all_possible_doas_lut[i].el, srp_phat_map_per_doa[i]);
        // }
        if (srp_phat_map_per_doa[i] > srp_phat_map_per_doa[ret]) {
            ret = i;
        }
    }

    return ret;
}


// void print_peak_memory_usage() {
//     PROCESS_MEMORY_COUNTERS pmc;
//     if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) {
//         printf("PeakWorkingSetSize: %llu bytes\n", (unsigned long long)pmc.PeakWorkingSetSize);
//         printf("Current WorkingSetSize: %llu bytes\n", (unsigned long long)pmc.WorkingSetSize);
//     }
// }


/* === main === */
int main() {
    kissfftrcfg = kiss_fftr_alloc(L, 0, NULL, NULL);
    uint32_t max = 0;

    if (kissfftrcfg == NULL) {
        printf("error fft\n");
        return 1;
    }

    prepare_luts();
    // for(uint8_t i = 0; i < M; i++) {
    //     printf("%f\t%f\t%f\n", mics_positions[i].x, mics_positions[i].y, mics_positions[i].z);
    // }

    if (acquire_mic_data() == MIC_ACQ_OK) {
        run_rfft();
        // for (int i = 0; i < REAL_FFT_LEN; i++) {
        //         printf("fft_mic[%d][%d] = %f\n", 1, i, cpx_abs(mic_data_fft_buff[1][i]));
        // }
        srp_phat();
        max = find_max_doa_idx();
        printf("znalezione maksimum w: { %f; %f }\n", all_possible_doas_lut[max].az, all_possible_doas_lut[max].el);
    }
    else {
        printf("error przyczytaniu\n");
    }

    KISS_FFT_FREE(kissfftrcfg);

    //print_peak_memory_usage();  // wypisze RAM po zakończeniu programu

    return 0;
}

/*
================ trash, don't mind this pls ======================
/ *

        for (uint32_t mic_l_offset = 1; mic_l_offset < M; mic_l_offset++) {

            for (uint32_t l = 0; l < ((L/2) + 1); l++) {
                tmp_cross_spec = cpx_mult(mic_data_fft_buff[mic_m + mic_l_offset][l], cpx_conj(mic_data_fft_buff[mic_m][l]));
                tmp_cross_spec_abs = cpx_abs(tmp_cross_spec);
                gcc_phat_per_doa[l].r = tmp_cross_spec.r / tmp_cross_spec_abs;
                gcc_phat_per_doa[l].i = tmp_cross_spec.i / tmp_cross_spec_abs;

                for (uint16_t f_idx = 0; f_idx < REAL_FFT_LEN; f_idx++) {
                    
                }
            }
        }

* /

    test prints for 'prepare_luts()':
    
        for (uint8_t i = 0; i < 25; i++) {
            printf("all_possible_doas_lut[%d].az = %f, all_possible_doas_lut[%d].el = %f\n", i, all_possible_doas_lut[i].az, i, all_possible_doas_lut[i].el);

            printf("all_relative_tdoas_lut[%d].angle_pair_idx = %d, all_relative_tdoas_lut[%d].mic_idx = %d, all_relative_tdoas_lut[%d].tau = %f\n", 
                    i, all_relative_tdoas_lut[i].angle_pair_idx, i, all_relative_tdoas_lut[i].mic_idx, i, all_relative_tdoas_lut[i].tau);
        }



        8*4 = 32

        micK[tN] = 4 B  ----> PCM : int32!!!!! 

        ramka1: [ mic1[t1] mic2[t1] mic3[t1] ... ]
        ramka2: [ mic1[t2] mic2[t2] mic3[t2] ... ]
        ...


*/




