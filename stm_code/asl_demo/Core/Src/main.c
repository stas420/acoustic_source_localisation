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

/* Private includes ----------------------------------------------------------*/
/* USER CODE BEGIN Includes */
#include "arm_math.h"
/* USER CODE END Includes */

/* Private typedef -----------------------------------------------------------*/
/* USER CODE BEGIN PTD */
typedef struct coords_3D_t {
    float32_t x;
    float32_t y;
    float32_t z;
} coords_3D_t;

typedef struct angle_pair_t {
    float32_t az;
    float32_t el;
} angle_pair_t;

typedef struct tdoa_per_mic_per_doa_t {
    float32_t tau;              // Changed to float for consistency
    uint8_t mic_idx;
    uint16_t angle_pair_idx;
} tdoa_per_mic_per_doa_t;
/* USER CODE END PTD */

/* Private define ------------------------------------------------------------*/
/* USER CODE BEGIN PD */
#define M 16                    /* num of mics */
#define Ny 4                    /* num of mics in y column */
#define Nz 4                    /* num of mics in z row */
#define C 343.0f                /* speed of sound in m/s */
#define d 0.042f                /* spacing between mics */
#define F_MIN 300.0f            /* minimum freq considered in Hz */
#define F_MAX 3000.0f           /* maximum freq considered in Hz */
#define F_S 16000.0f            /* incoming signal sampling rate */
#define L 512                   /* x_m microphone's signal frame length */
#define R 10                    /* angular grid resolution in degrees */
#define AZ_LOW -90              /* min azimuth angle considered in degrees */
#define AZ_HIGH 90              /* max azimuth angle considered in degrees */
#define EL_LOW -90              /* min elevation angle considered in degrees */
#define EL_HIGH 90              /* max elevation angle considered in degrees */

#ifndef PI
    #define PI 3.141592653589793f /* pi approximation */
#endif

#define EPS 1e-10f              /* epsilon value for numerical stability */

#define F_NYQ (F_S / 2.0f)                                  /* Nyquist frequency */
#define REAL_FFT_LEN ((L/2) + 1)
#define ALL_AZIMUTHS ((AZ_HIGH - AZ_LOW)/R + 1)             /* size of azimuth grid */
#define ALL_ELEVATIONS ((EL_HIGH - EL_LOW)/R + 1)           /* size of elevation grid */
#define ALL_POSSIBLE_DOAS (ALL_AZIMUTHS * ALL_ELEVATIONS)   /* size of all DoAs grid */

#define Q31_SIZE_BYTES (sizeof(int32_t))
#define SPI_DATA_BUFF_LEN (M * L * Q31_SIZE_BYTES)

/* Scaling factor for CIC filter output (PDM to PCM conversion)
 * CIC filters typically have gain = (R*M)^N where R=decimation, M=differential delay, N=stages
 * Adjust this based on your actual CIC configuration to prevent clipping
 * Typical range: 1.0f (no scaling) to 1e-6f (heavy decimation)
 */
#define CIC_SCALE_FACTOR (1.0f / 2147483648.0f)  /* TODO - possible change here, needs research */
/* USER CODE END PD */

/* Private macro -------------------------------------------------------------*/
/* USER CODE BEGIN PM */

/* USER CODE END PM */

/* Private variables ---------------------------------------------------------*/

SPI_HandleTypeDef hspi1;

/* USER CODE BEGIN PV */
const coords_3D_t mics_positions[M] = {
    {0.0f, -1.5f * d, -1.5f * d},
    {0.0f, -0.5f * d, -1.5f * d},
    {0.0f,  0.5f * d, -1.5f * d},
    {0.0f,  1.5f * d, -1.5f * d},

    {0.0f, -1.5f * d, -0.5f * d},
    {0.0f, -0.5f * d, -0.5f * d},
    {0.0f,  0.5f * d, -0.5f * d},
    {0.0f,  1.5f * d, -0.5f * d},

    {0.0f, -1.5f * d,  0.5f * d},
    {0.0f, -0.5f * d,  0.5f * d},
    {0.0f,  0.5f * d,  0.5f * d},
    {0.0f,  1.5f * d,  0.5f * d},

    {0.0f, -1.5f * d,  1.5f * d},
    {0.0f, -0.5f * d,  1.5f * d},
    {0.0f,  0.5f * d,  1.5f * d},
    {0.0f,  1.5f * d,  1.5f * d}
};

angle_pair_t all_possible_doas_lut[ALL_POSSIBLE_DOAS];
tdoa_per_mic_per_doa_t all_relative_tdoas_lut[ALL_POSSIBLE_DOAS * M];

uint16_t low_fft_bin_idx = 0;
uint16_t high_fft_bin_idx = 0;
float32_t fft_bins_lut[REAL_FFT_LEN];

uint8_t spi_mic_data_buff[SPI_DATA_BUFF_LEN];
float32_t float_mic_data_buff[M][L];
float32_t float_mic_data_fft_buff[M][L];
float32_t hann_window[L];

float32_t srp_phat_map[ALL_POSSIBLE_DOAS];

arm_rfft_fast_instance_f32 rfft_f32;
/* USER CODE END PV */

/* Private function prototypes -----------------------------------------------*/
void SystemClock_Config(void);
static void MX_GPIO_Init(void);
static void MX_SPI1_Init(void);
/* USER CODE BEGIN PFP */

/* USER CODE END PFP */

/* Private user code ---------------------------------------------------------*/
/* USER CODE BEGIN 0 */
/**
 * @brief Convert degrees to radians
 */
static inline float32_t deg2rad(float32_t arg_deg) {
    return (float32_t)(arg_deg * (PI / 180.0f));
}

/**
 * @brief Convert radians to degrees
 */
static inline float32_t rad2deg(float32_t arg_rad) {
    return (float32_t) (arg_rad * (180.0f / PI));
}

/**
 * @brief Prepare lookup tables for DoA grid and TDoA values
 * @note Called once during initialization
 */
static void prepare_luts(void) {
    coords_3D_t sphere_vec;
    uint16_t doa_idx = 0U;
    uint32_t tdoa_idx = 0U;
    uint8_t mic_idx = 0U;
    float32_t az_rad, el_rad;
    float32_t tmp_tau;

    /* prep Hann window */
    for (uint16_t i = 0; i < L; i++) {
        hann_window[i] = 0.5f * (1.0f - cosf(2.0f * PI * (float32_t)i / (float32_t)(L - 1)));
    }

    /* prep DoA grid and TDoA for each mic pair */
    for (uint16_t az_idx = 0; az_idx < ALL_AZIMUTHS; az_idx++) {
        for (uint16_t el_idx = 0; el_idx < ALL_ELEVATIONS; el_idx++) {
            /* DoA in degrees */
            all_possible_doas_lut[doa_idx].az = (float32_t)(AZ_LOW + (az_idx * R));
            all_possible_doas_lut[doa_idx].el = (float32_t)(EL_LOW + (el_idx * R));

            az_rad = deg2rad(all_possible_doas_lut[doa_idx].az);
            el_rad = deg2rad(all_possible_doas_lut[doa_idx].el);

            /* the reference mic is of 0th index, it has TDoA = 0 */
            all_relative_tdoas_lut[tdoa_idx].angle_pair_idx = doa_idx;
            all_relative_tdoas_lut[tdoa_idx].mic_idx = 0U;
            all_relative_tdoas_lut[tdoa_idx].tau = 0.0f;
            tdoa_idx++;

            /* prep spherical direction vector as per MATLAB convention:
             * sphere_vec.x = cos(el) * cos(az)
             * sphere_vec.y = cos(el) * sin(az)
             * sphere_vec.z = sin(el)
             *
             * they say cosf/sinf are slower than arm_cos_f32/sin
             * but this is our LUT, we need accuracy for once, right?
             */
            sphere_vec.x = cosf(el_rad) * cosf(az_rad);
            sphere_vec.y = cosf(el_rad) * sinf(az_rad);
            sphere_vec.z = sinf(el_rad);

            /* Compute relative TDoA for each mic relative to mic 0 */
            for (mic_idx = 1; mic_idx < M; mic_idx++) {
                all_relative_tdoas_lut[tdoa_idx].angle_pair_idx = doa_idx;
                all_relative_tdoas_lut[tdoa_idx].mic_idx = mic_idx;

                /* TDoA = (direction_vector · position_difference) / speed_of_sound */
                tmp_tau = (sphere_vec.x * (mics_positions[mic_idx].x - mics_positions[0].x) +
                          sphere_vec.y * (mics_positions[mic_idx].y - mics_positions[0].y) +
                          sphere_vec.z * (mics_positions[mic_idx].z - mics_positions[0].z)) / C;

                all_relative_tdoas_lut[tdoa_idx].tau = tmp_tau;
                tdoa_idx++;
            }

            doa_idx++;
        }
    }

    /*
     * TODO - check if omega pre calc is good
     */
    /* Prepare FFT bin frequencies and determine frequency range of interest */
    low_fft_bin_idx = (uint16_t)ceilf(F_MIN * (float32_t)L / F_S);
    high_fft_bin_idx = (uint16_t)floorf(F_MAX * (float32_t)L / F_S) + 1U;

    for (uint16_t k = 0U; k < (uint16_t) REAL_FFT_LEN; k++) {
        fft_bins_lut[k] = k * (F_S / (float32_t)L);
    }
}

/**
 * @brief Convert SPI buffer (int32_t CIC output) to float32_t format
 * @note CIC filter output is typically in Q31 format (32-bit signed integer)
 *       We convert to float [-1.0, 1.0] for FPU-optimized processing
 */
static void convert_buffer(void) {
    union {
        uint8_t u8[4];
        int32_t i32;
    } conv;

    uint32_t mic_frame_i = 0U;
    uint32_t spi_i = 0U;

    /* incoming SPI data format is as follows:
     * [mic0_sample0, mic1_sample0, ..., mic15_sample0, mic0_sample1, ...]
     * each sample is 4 bytes - a 32-bit signed integer from CIC filter
     */
    while (spi_i < SPI_DATA_BUFF_LEN) {
        for (uint32_t mic_i = 0U; mic_i < M; mic_i++) {
        	/* MSB data format? */
            conv.u8[3] = spi_mic_data_buff[spi_i];
            spi_i++;
            conv.u8[2] = spi_mic_data_buff[spi_i];
            spi_i++;
            conv.u8[1] = spi_mic_data_buff[spi_i];
            spi_i++;
            conv.u8[0] = spi_mic_data_buff[spi_i];
            spi_i++;

            float_mic_data_buff[mic_i][mic_frame_i] = ((float)conv.i32) * CIC_SCALE_FACTOR;
        }
        mic_frame_i++;
    }
}

/**
 * @brief get complex FFT bin (Re, Im) from arm_rfft_fast_f32 output
 * @param fft_buff Pointer to FFT output buffer
 * @param k Bin index (0 to L/2)
 * @param out_re_im Output array [Re, Im]
 *
 * @note arm_rfft_fast_f32 output format for L=512:
 *       { Re[0], Im[0], Re[1], Im[1], ... Re[L], Im[L] }
 *       DC (k=0) and Nyquist (k=256) have no imaginary part
 */
static inline void get_fft_bin(const float32_t* fft_buff, uint16_t k, float32_t* out_re_im) {
    if (k == 0) {
        out_re_im[0] = fft_buff[0];
        out_re_im[1] = fft_buff[1];  // Im[0] is typically 0 but stored at index 1
    } else if (k == (L / 2)) {
        out_re_im[0] = fft_buff[L];  // Re[N/2] is at index L (not L/2!)
        out_re_im[1] = 0.0f;
    } else {
        out_re_im[0] = fft_buff[2 * k];      // Re[k]
        out_re_im[1] = fft_buff[2 * k + 1];  // Im[k]
    }
}

/**
 * @brief the core SRP-PHAT algorithm for DoA estimation
 * @return estimated DoA as pair of angles { azimuth, elevation }
 */
static angle_pair_t srp_phat(void) {
    angle_pair_t ret_doa;
    uint32_t doa_idx = 0U;
    uint32_t max_doa_idx = 0U;

    /* -- 1. RFFT for all microphones' data -- */
    for (uint8_t m = 0U; m < M; m++) {
    	/* but also window the signal, for better data */
    	arm_mult_f32(&float_mic_data_buff[m][0], hann_window, &float_mic_data_buff[m][0], (uint32_t) L);
        arm_rfft_fast_f32(&rfft_f32, float_mic_data_buff[m], float_mic_data_fft_buff[m], 0);
    }

    /* -- 2. SRP-PHAT map init -- */
    arm_fill_f32(0.0f, srp_phat_map, (uint32_t) ALL_POSSIBLE_DOAS);

    /* -- 3. calc  SRP-PHAT for each DoA -- */
    float32_t bin_m[2], bin_l[2];
    float32_t cross_spec[2];
    float32_t cross_spec_mag;
    float32_t gcc_phat_sum[2];
    float32_t tau_ml, phase, omega;
    float32_t cos_phase, sin_phase;
    float32_t steered_re, steered_im;

    /* for each doa considered... */
    for (uint32_t doa = 0U; doa < ALL_POSSIBLE_DOAS; doa++) {
        doa_idx = doa * M;

        /* ...and for each mic pair... */
        for (uint8_t mic_m = 0U; mic_m < M; mic_m++) {
            for (uint8_t mic_l = mic_m + 1U; mic_l < M; mic_l++) {

                /* prep GCC-PHAT func */
                gcc_phat_sum[0] = 0.0f;
                gcc_phat_sum[1] = 0.0f;

                /* ...and for each FFT bin of interest */
                for (uint16_t k = low_fft_bin_idx; k < high_fft_bin_idx; k++) {

                    /* get correct FFT bins for both mics */
                    get_fft_bin(float_mic_data_fft_buff[mic_m], k, bin_m);
                    get_fft_bin(float_mic_data_fft_buff[mic_l], k, bin_l);

                    /* computing GCC-PHAT cross-spectrum (nominator) */
                    cross_spec[0] = bin_l[0] * bin_m[0] + bin_l[1] * bin_m[1];
                    cross_spec[1] = bin_l[1] * bin_m[0] - bin_l[0] * bin_m[1];

                    /* magnitude for PHAT normalization - this or that? */
                    //(void) arm_sqrt_f32(cross_spec[0] * cross_spec[0] + cross_spec[1] * cross_spec[1],
                    //					&cross_spec_mag);
                    arm_cmplx_mag_f32(cross_spec, &cross_spec_mag, (uint32_t) 1U);

                    /* skip if cross spectrum's magnitude is close to 0, leave it that way */
                    if (cross_spec_mag > EPS) {

                        /* apply PHAT normalization */
                        cross_spec[0] /= cross_spec_mag;
                        cross_spec[1] /= cross_spec_mag;

                        /* calc TDoA for this mic pair and DoA */
                        tau_ml =   all_relative_tdoas_lut[doa_idx + mic_l].tau
                        		 - all_relative_tdoas_lut[doa_idx + mic_m].tau;

                        /* prep phase shift */
                        omega = 2.0f * PI * fft_bins_lut[k];
                        phase = fmodf(omega * tau_ml, 2.0f * PI);

                        if (phase > PI) {
                        	phase -= 2.0f * PI;
                        }

                        if (phase < -PI) {
                        	phase += 2.0f * PI;
                        }

                        /* prep steering vec: e^(j*phase) = cos(phase) + j*sin(phase) */
                        cos_phase = arm_cos_f32(phase);
                        sin_phase = arm_sin_f32(phase);

                        /* apply steering vec */
                        steered_re = cross_spec[0] * cos_phase - cross_spec[1] * sin_phase;
                        steered_im = cross_spec[0] * sin_phase + cross_spec[1] * cos_phase;

                        /* add to this GCC-PHAT function */
                        gcc_phat_sum[0] += steered_re;
                        gcc_phat_sum[1] += steered_im;
                    }
                }

                /* add its power to SRP-PHAT map */
                srp_phat_map[doa] +=   gcc_phat_sum[0] * gcc_phat_sum[0]
									 + gcc_phat_sum[1] * gcc_phat_sum[1];
            }
        }
    }

    /* -- 4. find maximum in SRP-PHAT map and return it -- */
    float32_t max_power;
    arm_max_f32(srp_phat_map, ALL_POSSIBLE_DOAS, &max_power, &max_doa_idx);
    ret_doa = all_possible_doas_lut[max_doa_idx];

    return ret_doa;
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
  arm_status arm_stat;

  #if L == 512
    arm_stat = arm_rfft_fast_init_512_f32(&rfft_f32);
  #else
	arm_stat = arm_rfft_fast_init_f32(&rfft_f32, (uint16_t) L);
  #endif

  /* well, if RFFT doesn't work, we don't have anything to play with */
  if (arm_stat != ARM_MATH_SUCCESS) {
	  Error_Handler();
  }
  /* USER CODE END SysInit */

  /* Initialize all configured peripherals */
  MX_GPIO_Init();
  MX_SPI1_Init();
  /* USER CODE BEGIN 2 */
  const uint32_t spi_timeout = 100U;
  HAL_StatusTypeDef status = HAL_OK;
  angle_pair_t doa;
  /* USER CODE END 2 */

  /* Infinite loop */
  /* USER CODE BEGIN WHILE */
  while (1)
  {
	status = HAL_SPI_Receive(&hspi1, spi_mic_data_buff, (uint32_t) SPI_DATA_BUFF_LEN, spi_timeout);

	if (status == HAL_OK) {
		convert_buffer();
		doa = srp_phat();
		/* TODO: send DoA thru UART! */
	}
	else {

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
  hspi1.Init.DataSize = SPI_DATASIZE_32BIT;
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
