/* USER CODE BEGIN Header */
/**
  ******************************************************************************
  * @file           : main.c
  * @brief          : Main program body
  ******************************************************************************
  * @attention
  *
  * Copyright (c) 2025 STMicroelectronics.
  * All rights reserved.
  *
  * This software is licensed under terms that can be found in the LICENSE file
  * in the root directory of this software component.
  * If no LICENSE file comes with this software, it is provided AS-IS.
  *
  ******************************************************************************
  */
/* USER CODE END Header */
/* Includes ------------------------------------------------------------------*/
#include "main.h"
#include "test_data.h"
#include "kiss_fftr.h"
/* Private includes ----------------------------------------------------------*/
/* USER CODE BEGIN Includes */
/* USER CODE END Includes */

/* Private typedef -----------------------------------------------------------*/
/* USER CODE BEGIN PTD */
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

/* USER CODE END PTD */

/* Private define ------------------------------------------------------------*/
/* USER CODE BEGIN PD */
#define M 16                    /* num of mics */
#define Ny 4                    /* num of mics in y column */
#define Nz 4                    /* num of mics in z row */
#define C 343                   /* speed of sound in m/s */
#define d 0.042                 /* spacing between mics */
#define F_MIN 300               /* minimum freq considered in Hz */
#define F_MAX 3000              /* minimum freq considered in Hz */
#define F_S 16000               /* incoming signal sampling rate */
#define L 1024                   /* x_m microphone's signal frame length - TO BE DETERMINED */
#define R 10                    /* angualar grid resolution in degrees */
#define AZ_LOW -90              /* min azimuth angle conidered in degrees */
#define AZ_HIGH 90              /* max azimuth angle conidered in degrees */
#define EL_LOW -45              /* min elevation angle considered in degrees */
#define EL_HIGH 45              /* max elevation angle considered in degrees */
#define PI 3.141592653589793238 /* pi approximation */
#define EPS 0.000000001           /* epsilon value for numerical stability */

#define F_NYQ (F_S / 2)                                     /* Nyquist frequency */
#define REAL_FFT_LEN ((L/2) + 1)                            /* length of the output of real-valued signal's FFT, which is of use */
#define ALL_AZIMUTHS ((AZ_HIGH - AZ_LOW)/R + 1)             /* size of azimuth grid */
#define ALL_ELEVATIONS ((EL_HIGH - EL_LOW)/R + 1)           /* size of elevation grid */
#define ALL_POSSIBLE_DOAS (ALL_AZIMUTHS * ALL_ELEVATIONS)   /* size of all DoAs grid (i.e. every possible pair of {az,el}) */

#define CIC_GAIN_F ((float) (1ULL << 30))  /* CIC has its internal gain like (RM)^N, we have (64*1)^5 = 2^30 */
/*
  CIC params: stages 5, diff delay 1, dec ratio 64
  TODO:
  - add windowing b4 fft
  - remove DC??? depends if it's zero-mean
  - add scaling down -- quantization mode: full precision
  - 
*/
/* USER CODE END PD */

/* Private macro -------------------------------------------------------------*/
/* USER CODE BEGIN PM */

/* USER CODE END PM */

/* Private variables ---------------------------------------------------------*/

SPI_HandleTypeDef hspi1;
DMA_HandleTypeDef handle_GPDMA1_Channel0;

/* USER CODE BEGIN PV */
/**
 * @brief simple enum for SPI status indications
 */
enum { SPI_CPLT, SPI_BUSY } spi_status;

coords_3D_t mics_positions[M] = {
  { 0, (-2.0f)*d, 2*d }, { 0, (-2.0f)*d, d }, {  0, (-2.0f)*d, (-1.0f)*d }, {  0, (-2.0f)*d, (-2.0f)*d },
  { 0, (-1.0f)*d, 2*d }, { 0, (-1.0f)*d, d }, {  0, (-1.0f)*d, (-1.0f)*d }, {  0, (-1.0f)*d, (-2.0f)*d },
  { 0,         d, 2*d }, { 0,         d, d }, {  0,         d, (-1.0f)*d }, {  0,         d, (-2.0f)*d },
  { 0,  (2.0f)*d, 2*d }, { 0,  (2.0f)*d, d }, {  0,  (2.0f)*d, (-1.0f)*d }, {  0,  (2.0f)*d, (-2.0f)*d }
};
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
 * @brief buffer for incoming mic data from SPI
 */
uint32_t spi_data_buff[M * L] = { 0U };
#ifndef TEST_DATA_H
/**
 * @brief preallocated buffer for incoming real signal data from mics (thru SPI)
 */
float mic_data_buff[M][L] = { 0 };//<<  65,536 B
#endif

/**
 * @brief preallocated buffer for incoming real signal data FFT 
 */
static kiss_fft_cpx mic_data_fft_buff[M][REAL_FFT_LEN]; //<< 65,664 B

/**
 * @brief preallocated buffer for SRP-PHAT map (as function of DoA angles pair)
 * 
 */
static float srp_phat_map_per_doa[ALL_POSSIBLE_DOAS];//<< 1,444 B

/* USER CODE END PV */

/* Private function prototypes -----------------------------------------------*/
void SystemClock_Config(void);
static void MX_GPIO_Init(void);
static void MX_GPDMA1_Init(void);
static void MX_SPI1_Init(void);
/* USER CODE BEGIN PFP */

/* USER CODE END PFP */

/* Private user code ---------------------------------------------------------*/
/* USER CODE BEGIN 0 */
/**
 * @brief helper function for data buffer conversion process (SPI(CIC): uint32_t -> kissfft: float32_t)
 * 
 */
void convert_buffer(void) {
    static int32_t tmp_24bit = 0;
    
    for (uint32_t l = 0U; l < L; l++) {
        for (uint32_t m = 0U; m < M; m++) {
            // spi_data_buff[m + M * l] <--> mic_data_buff[m][l]
            /* get rid of the last 8-bit marker */
            tmp_24bit = (int32_t) (spi_data_buff[m + M * l] >> 8);
            
            /* extend the sign of CIC sample */
            if (tmp_24bit & 0x00800000) {
                tmp_24bit |= 0xFF000000;
            }
            
            /* compensate the CIC gain - from vivado it seems to be not done internally */
            mic_data_buff[m][l] = ((float) tmp_24bit) / CIC_GAIN_F;
        }
    }
}

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
    coords_3D_t sphere_vec;
    uint16_t doa_idx = 0;
    uint32_t tdoa_idx = 0;
    uint16_t el_idx = 0;
    uint8_t mic_idx = 0;
    float az_rad = 0.0f;
    float el_rad = 0.0f;

    /* for every azimuth degree step */
    for (uint16_t az_idx = 0; az_idx < ALL_AZIMUTHS; az_idx++) {
        /* and for every elevation degree step */
        for (el_idx = 0; el_idx < ALL_ELEVATIONS; el_idx++) {
            /* define an entry in degree-grid */
            all_possible_doas_lut[doa_idx].az = (float) (AZ_LOW + (az_idx * R));
            all_possible_doas_lut[doa_idx].el = (float) (EL_LOW + (el_idx * R));
            
            /* Convert to radians ONCE for this DoA */
            az_rad = deg2rad(all_possible_doas_lut[doa_idx].az);
            el_rad = deg2rad(all_possible_doas_lut[doa_idx].el);
            
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
    uint32_t doa_base_idx = 0;
    
    // Initialize map
    for (uint16_t i = 0; i < ALL_POSSIBLE_DOAS; i++) {
        srp_phat_map_per_doa[i] = 0.0f;
    }
    
    /* For each DoA... */
    for (uint16_t doa = 0; doa < ALL_POSSIBLE_DOAS; doa++) {
        // Calculate base index for this DoA in the TDoA LUT
        doa_base_idx = (uint32_t)doa * (uint32_t)M;
        
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
uint32_t find_max_doa_idx(void) {
    uint32_t ret = 0;
    for (uint32_t i = 0; i < ALL_POSSIBLE_DOAS; i++) {
        if (srp_phat_map_per_doa[i] > srp_phat_map_per_doa[ret]) {
            ret = i;
        }
    }

    return ret;
}

/**
 * @brief custom handler for full buffer DMA event
 * 
 * @param hspi SPI handler
 */
void HAL_SPI_RxCpltCallback(SPI_HandleTypeDef *hspi) {
	  spi_status = SPI_CPLT;
}
/* USER CODE END 0 */

/**
  * @brief  The application entry point.
  * @retval int
  */
int main(void)
{

  /* USER CODE BEGIN 1 */

  /* USER CODE END 1 */

  /* MCU Configuration--------------------------------------------------------*/

  /* Reset of all peripherals, Initializes the Flash interface and the Systick. */
  HAL_Init();

  /* USER CODE BEGIN Init */

  /* USER CODE END Init */

  /* Configure the system clock */
  SystemClock_Config();

  /* USER CODE BEGIN SysInit */
  prepare_luts();
  #ifndef TEST_DATA_H
  	  spi_status = SPI_BUSY;
  #else
  	  spi_status = SPI_CPLT;
  #endif

  /* USER CODE END SysInit */

  /* Initialize all configured peripherals */
  MX_GPIO_Init();
  MX_GPDMA1_Init();
  MX_SPI1_Init();
  /* USER CODE BEGIN 2 */
  kissfftrcfg = kiss_fftr_alloc(L, 0, NULL, NULL);

  if (kissfftrcfg == NULL) {
	  HardFault_Handler();
  }

  #ifndef TEST_DATA_H
    HAL_SPI_Receive_DMA(&hspi1, (uint8_t*) spi_data_buff, (uint16_t) (M * L * sizeof(uint32_t)));
  #endif

  uint32_t doa_idx = 0U;
  angle_pair_t doa_est;
  /* USER CODE END 2 */

  /* Infinite loop */
  /* USER CODE BEGIN WHILE */
  while (1)
  {
      if (spi_status == SPI_CPLT) {
		  #ifndef TEST_DATA_H
             convert_buffer();
          #endif
          run_rfft();
          srp_phat();
          doa_idx = find_max_doa_idx();
          doa_est = all_possible_doas_lut[doa_idx];

		  #ifndef TEST_DATA_H
              spi_status = SPI_BUSY;
              HAL_SPI_Receive_DMA(&hspi1, (uint8_t*) spi_data_buff, (uint16_t) (M * L * sizeof(uint32_t)));
		  #endif
      }
    /* USER CODE END WHILE */

    /* USER CODE BEGIN 3 */
  }
  /* USER CODE END 3 */
}

/**
  * @brief System Clock Configuration
  * @retval None
  */
void SystemClock_Config(void)
{
  RCC_OscInitTypeDef RCC_OscInitStruct = {0};
  RCC_ClkInitTypeDef RCC_ClkInitStruct = {0};

  /** Configure the main internal regulator output voltage
  */
  __HAL_PWR_VOLTAGESCALING_CONFIG(PWR_REGULATOR_VOLTAGE_SCALE3);

  while(!__HAL_PWR_GET_FLAG(PWR_FLAG_VOSRDY)) {}

  /** Initializes the RCC Oscillators according to the specified parameters
  * in the RCC_OscInitTypeDef structure.
  */
  RCC_OscInitStruct.OscillatorType = RCC_OSCILLATORTYPE_HSE;
  RCC_OscInitStruct.HSEState = RCC_HSE_BYPASS;
  RCC_OscInitStruct.PLL.PLLState = RCC_PLL_ON;
  RCC_OscInitStruct.PLL.PLLSource = RCC_PLL1_SOURCE_HSE;
  RCC_OscInitStruct.PLL.PLLM = 1;
  RCC_OscInitStruct.PLL.PLLN = 16;
  RCC_OscInitStruct.PLL.PLLP = 2;
  RCC_OscInitStruct.PLL.PLLQ = 3;
  RCC_OscInitStruct.PLL.PLLR = 2;
  RCC_OscInitStruct.PLL.PLLRGE = RCC_PLL1_VCIRANGE_3;
  RCC_OscInitStruct.PLL.PLLVCOSEL = RCC_PLL1_VCORANGE_WIDE;
  RCC_OscInitStruct.PLL.PLLFRACN = 0;
  if (HAL_RCC_OscConfig(&RCC_OscInitStruct) != HAL_OK)
  {
    Error_Handler();
  }

  /** Initializes the CPU, AHB and APB buses clocks
  */
  RCC_ClkInitStruct.ClockType = RCC_CLOCKTYPE_HCLK|RCC_CLOCKTYPE_SYSCLK
                              |RCC_CLOCKTYPE_PCLK1|RCC_CLOCKTYPE_PCLK2
                              |RCC_CLOCKTYPE_PCLK3;
  RCC_ClkInitStruct.SYSCLKSource = RCC_SYSCLKSOURCE_PLLCLK;
  RCC_ClkInitStruct.AHBCLKDivider = RCC_SYSCLK_DIV1;
  RCC_ClkInitStruct.APB1CLKDivider = RCC_HCLK_DIV1;
  RCC_ClkInitStruct.APB2CLKDivider = RCC_HCLK_DIV1;
  RCC_ClkInitStruct.APB3CLKDivider = RCC_HCLK_DIV1;

  if (HAL_RCC_ClockConfig(&RCC_ClkInitStruct, FLASH_LATENCY_3) != HAL_OK)
  {
    Error_Handler();
  }

  /** Configure the programming delay
  */
  __HAL_FLASH_SET_PROGRAM_DELAY(FLASH_PROGRAMMING_DELAY_1);
}

/**
  * @brief GPDMA1 Initialization Function
  * @param None
  * @retval None
  */
static void MX_GPDMA1_Init(void)
{

  /* USER CODE BEGIN GPDMA1_Init 0 */

  /* USER CODE END GPDMA1_Init 0 */

  /* Peripheral clock enable */
  __HAL_RCC_GPDMA1_CLK_ENABLE();

  /* GPDMA1 interrupt Init */
    HAL_NVIC_SetPriority(GPDMA1_Channel0_IRQn, 0, 0);
    HAL_NVIC_EnableIRQ(GPDMA1_Channel0_IRQn);

  /* USER CODE BEGIN GPDMA1_Init 1 */

  /* USER CODE END GPDMA1_Init 1 */
  /* USER CODE BEGIN GPDMA1_Init 2 */

  /* USER CODE END GPDMA1_Init 2 */

}

/**
  * @brief SPI1 Initialization Function
  * @param None
  * @retval None
  */
static void MX_SPI1_Init(void)
{

  /* USER CODE BEGIN SPI1_Init 0 */

  /* USER CODE END SPI1_Init 0 */

  /* USER CODE BEGIN SPI1_Init 1 */

  /* USER CODE END SPI1_Init 1 */
  /* SPI1 parameter configuration*/
  hspi1.Instance = SPI1;
  hspi1.Init.Mode = SPI_MODE_SLAVE;
  hspi1.Init.Direction = SPI_DIRECTION_2LINES_RXONLY;
  hspi1.Init.DataSize = SPI_DATASIZE_8BIT;
  hspi1.Init.CLKPolarity = SPI_POLARITY_HIGH;
  hspi1.Init.CLKPhase = SPI_PHASE_1EDGE;
  hspi1.Init.NSS = SPI_NSS_HARD_INPUT;
  hspi1.Init.FirstBit = SPI_FIRSTBIT_MSB;
  hspi1.Init.TIMode = SPI_TIMODE_DISABLE;
  hspi1.Init.CRCCalculation = SPI_CRCCALCULATION_DISABLE;
  hspi1.Init.CRCPolynomial = 0x7;
  hspi1.Init.NSSPMode = SPI_NSS_PULSE_DISABLE;
  hspi1.Init.NSSPolarity = SPI_NSS_POLARITY_LOW;
  hspi1.Init.FifoThreshold = SPI_FIFO_THRESHOLD_01DATA;
  hspi1.Init.MasterSSIdleness = SPI_MASTER_SS_IDLENESS_00CYCLE;
  hspi1.Init.MasterInterDataIdleness = SPI_MASTER_INTERDATA_IDLENESS_00CYCLE;
  hspi1.Init.MasterReceiverAutoSusp = SPI_MASTER_RX_AUTOSUSP_DISABLE;
  hspi1.Init.MasterKeepIOState = SPI_MASTER_KEEP_IO_STATE_DISABLE;
  hspi1.Init.IOSwap = SPI_IO_SWAP_DISABLE;
  hspi1.Init.ReadyMasterManagement = SPI_RDY_MASTER_MANAGEMENT_INTERNALLY;
  hspi1.Init.ReadyPolarity = SPI_RDY_POLARITY_HIGH;
  if (HAL_SPI_Init(&hspi1) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN SPI1_Init 2 */

  /* USER CODE END SPI1_Init 2 */

}

/**
  * @brief GPIO Initialization Function
  * @param None
  * @retval None
  */
static void MX_GPIO_Init(void)
{
  GPIO_InitTypeDef GPIO_InitStruct = {0};
  /* USER CODE BEGIN MX_GPIO_Init_1 */

  /* USER CODE END MX_GPIO_Init_1 */

  /* GPIO Ports Clock Enable */
  __HAL_RCC_GPIOE_CLK_ENABLE();
  __HAL_RCC_GPIOC_CLK_ENABLE();
  __HAL_RCC_GPIOF_CLK_ENABLE();
  __HAL_RCC_GPIOH_CLK_ENABLE();
  __HAL_RCC_GPIOA_CLK_ENABLE();
  __HAL_RCC_GPIOB_CLK_ENABLE();
  __HAL_RCC_GPIOD_CLK_ENABLE();
  __HAL_RCC_GPIOG_CLK_ENABLE();

  /*Configure GPIO pin Output Level */
  HAL_GPIO_WritePin(GPIOF, GPIO_PIN_4, GPIO_PIN_RESET);

  /*Configure GPIO pin Output Level */
  HAL_GPIO_WritePin(GPIOB, GPIO_PIN_0, GPIO_PIN_RESET);

  /*Configure GPIO pin Output Level */
  HAL_GPIO_WritePin(GPIOG, GPIO_PIN_4, GPIO_PIN_RESET);

  /*Configure GPIO pin : PC13 */
  GPIO_InitStruct.Pin = GPIO_PIN_13;
  GPIO_InitStruct.Mode = GPIO_MODE_IT_RISING;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  HAL_GPIO_Init(GPIOC, &GPIO_InitStruct);

  /*Configure GPIO pin : PF4 */
  GPIO_InitStruct.Pin = GPIO_PIN_4;
  GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_LOW;
  HAL_GPIO_Init(GPIOF, &GPIO_InitStruct);

  /*Configure GPIO pins : RMII_MDC_Pin RMII_RXD0_Pin RMII_RXD1_Pin */
  GPIO_InitStruct.Pin = RMII_MDC_Pin|RMII_RXD0_Pin|RMII_RXD1_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_AF_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_HIGH;
  GPIO_InitStruct.Alternate = GPIO_AF11_ETH;
  HAL_GPIO_Init(GPIOC, &GPIO_InitStruct);

  /*Configure GPIO pins : RMII_REF_CLK_Pin RMII_MDIO_Pin RMII_CRS_DV_Pin */
  GPIO_InitStruct.Pin = RMII_REF_CLK_Pin|RMII_MDIO_Pin|RMII_CRS_DV_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_AF_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_HIGH;
  GPIO_InitStruct.Alternate = GPIO_AF11_ETH;
  HAL_GPIO_Init(GPIOA, &GPIO_InitStruct);

  /*Configure GPIO pin : VBUS_SENSE_Pin */
  GPIO_InitStruct.Pin = VBUS_SENSE_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_ANALOG;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  HAL_GPIO_Init(VBUS_SENSE_GPIO_Port, &GPIO_InitStruct);

  /*Configure GPIO pin : PB0 */
  GPIO_InitStruct.Pin = GPIO_PIN_0;
  GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_LOW;
  HAL_GPIO_Init(GPIOB, &GPIO_InitStruct);

  /*Configure GPIO pins : UCPD_CC1_Pin UCPD_CC2_Pin */
  GPIO_InitStruct.Pin = UCPD_CC1_Pin|UCPD_CC2_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_ANALOG;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  HAL_GPIO_Init(GPIOB, &GPIO_InitStruct);

  /*Configure GPIO pin : RMII_TXD1_Pin */
  GPIO_InitStruct.Pin = RMII_TXD1_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_AF_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_HIGH;
  GPIO_InitStruct.Alternate = GPIO_AF11_ETH;
  HAL_GPIO_Init(RMII_TXD1_GPIO_Port, &GPIO_InitStruct);

  /*Configure GPIO pins : PD8 PD9 */
  GPIO_InitStruct.Pin = GPIO_PIN_8|GPIO_PIN_9;
  GPIO_InitStruct.Mode = GPIO_MODE_AF_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_LOW;
  GPIO_InitStruct.Alternate = GPIO_AF7_USART3;
  HAL_GPIO_Init(GPIOD, &GPIO_InitStruct);

  /*Configure GPIO pin : PG4 */
  GPIO_InitStruct.Pin = GPIO_PIN_4;
  GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_LOW;
  HAL_GPIO_Init(GPIOG, &GPIO_InitStruct);

  /*Configure GPIO pin : UCPD_FLT_Pin */
  GPIO_InitStruct.Pin = UCPD_FLT_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_IT_RISING;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  HAL_GPIO_Init(UCPD_FLT_GPIO_Port, &GPIO_InitStruct);

  /*Configure GPIO pins : USB_FS_N_Pin USB_FS_P_Pin */
  GPIO_InitStruct.Pin = USB_FS_N_Pin|USB_FS_P_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_AF_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_LOW;
  GPIO_InitStruct.Alternate = GPIO_AF10_USB;
  HAL_GPIO_Init(GPIOA, &GPIO_InitStruct);

  /*Configure GPIO pins : RMII_TXT_EN_Pin RMI_TXD0_Pin */
  GPIO_InitStruct.Pin = RMII_TXT_EN_Pin|RMI_TXD0_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_AF_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_HIGH;
  GPIO_InitStruct.Alternate = GPIO_AF11_ETH;
  HAL_GPIO_Init(GPIOG, &GPIO_InitStruct);

  /*Configure GPIO pins : ARD_D1_TX_Pin ARD_D0_RX_Pin */
  GPIO_InitStruct.Pin = ARD_D1_TX_Pin|ARD_D0_RX_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_AF_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_LOW;
  GPIO_InitStruct.Alternate = GPIO_AF8_LPUART1;
  HAL_GPIO_Init(GPIOB, &GPIO_InitStruct);

  /* USER CODE BEGIN MX_GPIO_Init_2 */

  /* USER CODE END MX_GPIO_Init_2 */
}

/* USER CODE BEGIN 4 */

/* USER CODE END 4 */

/**
  * @brief  This function is executed in case of error occurrence.
  * @retval None
  */
void Error_Handler(void)
{
  /* USER CODE BEGIN Error_Handler_Debug */
  /* User can add his own implementation to report the HAL error return state */
  __disable_irq();
  while (1)
  {
  }
  /* USER CODE END Error_Handler_Debug */
}
#ifdef USE_FULL_ASSERT
/**
  * @brief  Reports the name of the source file and the source line number
  *         where the assert_param error has occurred.
  * @param  file: pointer to the source file name
  * @param  line: assert_param error line source number
  * @retval None
  */
void assert_failed(uint8_t *file, uint32_t line)
{
  /* USER CODE BEGIN 6 */
  /* User can add his own implementation to report the file name and line number,
     ex: printf("Wrong parameters value: file %s on line %d\r\n", file, line) */
  /* USER CODE END 6 */
}
#endif /* USE_FULL_ASSERT */
