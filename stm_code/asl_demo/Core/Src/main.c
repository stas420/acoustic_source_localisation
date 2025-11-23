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
    float x;
    float y;
    float z;
} coords_3D_t;

typedef struct angle_pair_t {
    float az;
    float el;
} angle_pair_t;

typedef struct tdoa_per_mic_per_doa_t {
    q31_t tau;
    uint8_t mic_idx;
    uint16_t angle_pair_idx;
} tdoa_per_mic_per_doa_t;

typedef struct q31_cpx_t {
	q31_t r;
	q31_t i;
} q31_cpx_t;
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
#define L 512                   /* x_m microphone's signal frame length - TO BE DETERMINED */
#define R 10                    /* angualar grid resolution in degrees */
#define AZ_LOW -90              /* min azimuth angle conidered in degrees */
#define AZ_HIGH 90              /* max azimuth angle conidered in degrees */
#define EL_LOW -90              /* min elevation angle considered in degrees */
#define EL_HIGH 90              /* max elevation angle considered in degrees */

#ifndef PI
	#define PI 3.141592653589793238 /* pi approximation */
#endif

#define EPS ((q31_t) 1)        /* epsilon value for numerical stability */

#define F_NYQ (F_S / 2)                                     /* Nyquist frequency */
#define REAL_FFT_LEN ((L/2) + 1)                            /* length of the output of real-valued signal's FFT, which is of use */
#define ALL_AZIMUTHS ((AZ_HIGH - AZ_LOW)/R + 1)             /* size of azimuth grid */
#define ALL_ELEVATIONS ((EL_HIGH - EL_LOW)/R + 1)           /* size of elevation grid */
#define ALL_POSSIBLE_DOAS (ALL_AZIMUTHS * ALL_ELEVATIONS)   /* size of all DoAs grid (i.e. every possible pair of {az,el}) */

#define Q31_SIZE_BYTES (sizeof(q31_t))
#define SPI_DATA_BUFF_LEN (M * L * Q31_SIZE_BYTES)
/* USER CODE END PD */

/* Private macro -------------------------------------------------------------*/
/* USER CODE BEGIN PM */

/* USER CODE END PM */

/* Private variables ---------------------------------------------------------*/

SPI_HandleTypeDef hspi1;

/* USER CODE BEGIN PV */
const coords_3D_t mics_positions[M] = {
	{ 0,   0, 3*d }, { 0,   0, 2*d }, { 0,   0, d }, { 0,   0, 0 },
	{ 0,   d, 3*d }, { 0,   d, 2*d }, { 0,   d, d }, { 0,   d, 0 },
	{ 0, 2*d, 3*d }, { 0, 2*d, 2*d }, { 0, 2*d, d }, { 0, 2*d, 0 },
	{ 0, 3*d, 3*d }, { 0, 3*d, 2*d }, { 0, 3*d, d }, { 0, 3*d, 0 },
};

angle_pair_t all_possible_doas_lut[ALL_POSSIBLE_DOAS];
tdoa_per_mic_per_doa_t all_relative_tdoas_lut[(ALL_POSSIBLE_DOAS * M)];

uint16_t low_fft_bin_idx = 0;
uint16_t high_fft_bin_idx = 0;
float fft_bins_lut[REAL_FFT_LEN];

uint8_t spi_mic_data_buff[SPI_DATA_BUFF_LEN];
q31_t q31_mic_data_buff[M][L];
q31_t q31_mic_data_fft_buff[M][L];

float srp_phat_mic_per_doa[ALL_POSSIBLE_DOAS];

arm_rfft_instance_q31 rfft;
/* USER CODE END PV */

/* Private function prototypes -----------------------------------------------*/
void SystemClock_Config(void);
static void MX_GPIO_Init(void);
static void MX_SPI1_Init(void);
/* USER CODE BEGIN PFP */

/* USER CODE END PFP */

/* Private user code ---------------------------------------------------------*/
/* USER CODE BEGIN 0 */
static inline float deg2rad(float arg_deg) {
    return (arg_deg * (float)(PI / 180.0f));
}

static inline float rad2deg(float arg_rad) {
    return (arg_rad * (float)(180.0f / PI));
}

static void prepare_luts(void) {
    coords_3D_t sphere_vec;
    uint16_t doa_idx = 0U;
    uint32_t tdoa_idx = 0U;
    uint16_t el_idx = 0U;
    uint8_t mic_idx = 0U;
    float az_rad = 0.0f;
    float el_rad = 0.0f;
    float tmp_tau = 0.0f;

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
            all_relative_tdoas_lut[tdoa_idx].mic_idx = (uint8_t) 0U;
            all_relative_tdoas_lut[tdoa_idx].tau = (q31_t) 0U;
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
                tmp_tau = (float) (
                    (sphere_vec.x * (mics_positions[mic_idx].x - mics_positions[0].x))
                  + (sphere_vec.y * (mics_positions[mic_idx].y - mics_positions[0].y))
                  + (sphere_vec.z * (mics_positions[mic_idx].z - mics_positions[0].z))
                ) / ((float) C);
                arm_float_to_q31(&tmp_tau, &(all_relative_tdoas_lut[tdoa_idx].tau), (uint32_t) 1U);
                tdoa_idx++;
            }

            doa_idx++;
        }
    }

    /* prepare every FFT bin and indicators */
    low_fft_bin_idx = 0U;
    high_fft_bin_idx = 0U;

    for (uint16_t k = 0U; k < REAL_FFT_LEN; k++) {
        fft_bins_lut[k] = k * ((float)F_S / (float)L);

        // Count bins below F_MIN
        if (fft_bins_lut[k] < F_MIN) {
        	low_fft_bin_idx++;
        }

        // Count bins up to AND including F_MAX
        if (fft_bins_lut[k] <= F_MAX) {
        	high_fft_bin_idx = k + 1U;  // Set to NEXT index (for < comparison in loop)
        }
    }
}

static void convert_buffer(void) {
	static union {
		uint8_t u8[4];
		q31_t q31;
	} conv;

	uint32_t mic_frame_i = 0U;
	uint32_t spi_i = 0U;

	while (spi_i < (uint32_t) SPI_DATA_BUFF_LEN) {
		for (uint32_t mic_i = 0U; mic_i < (uint32_t) M; mic_i++) {
			conv.u8[0] = spi_mic_data_buff[spi_i];
			spi_i++;
			conv.u8[1] = spi_mic_data_buff[spi_i];
			spi_i++;
			conv.u8[2] = spi_mic_data_buff[spi_i];
			spi_i++;
			conv.u8[3] = spi_mic_data_buff[spi_i];
			spi_i++;

			q31_mic_data_buff[mic_i][mic_frame_i] = conv.q31;
		}
		mic_frame_i++;
	}
}

static angle_pair_t srp_phat(void) {
	angle_pair_t ret_doa;

	/* -- 1. rfft -- */
	for (uint8_t m = 0U; m < M; m++) {
		arm_rfft_q31(&rfft, &q31_mic_data_buff[m][0], &q31_mic_data_fft_buff[m][0]);
	}

	/* -- 2. SRP-PHAT map init -- */
	for (uint32_t i = 0U; i < ALL_POSSIBLE_DOAS; i++) {
		srp_phat_mic_per_doa[i] = 0.0f;
	}

	/* -- 3. compute the core map of SRP-PHAT -- */
	uint32_t doa_idx = 0U;
	q31_t tau_ml = (q31_t) 0;
	q31_t gcc_phat_sum[2] = { 0 };  // Akumulator dla pary mic [Re, Im]
	q31_t tmp_cpx_conj[2] = { 0 };
	q31_t tmp_cross_spec[2] = { 0 };
	q31_t tmp_cross_spec_abs = (q31_t) 0;
	q31_t exp_phase[2] = { 0 };
	q31_t result[2] = { 0 };
	q31_t phase_q31 = 0;
	float tau_float = 0.0f;
	float phase = 0.0f;
	float omega = 0.0f;
	q63_t tmp_gcc_pow = 0;

	/* for each DoA... */
	for (uint32_t doa = 0U; doa < ALL_POSSIBLE_DOAS; doa++) {
		doa_idx = doa * ((uint32_t) M);

		/* ...and for each mic pair ml, l > m... */
		for (uint8_t mic_m = 0U; mic_m < M; mic_m++) {
			for (uint8_t mic_l = mic_m + 1U; mic_l < M; mic_l++) {

				/* Reset akumulatora GCC-PHAT dla tej pary mikrofonów */
				gcc_phat_sum[0] = 0;
				gcc_phat_sum[1] = 0;

				/* ...and for each freq bin of interest... */
				for (uint16_t k = low_fft_bin_idx; k < high_fft_bin_idx; k++) {
					/* ...we do a GCC-PHAT: */
					/* pre-calc conj */
					arm_cmplx_conj_q31(&q31_mic_data_fft_buff[mic_m][2*k], tmp_cpx_conj, 1U);
					/* calc cross-spectrum */
					arm_cmplx_mult_cmplx_q31(&q31_mic_data_fft_buff[mic_l][2*k],
											 tmp_cpx_conj,
											 tmp_cross_spec, 1U);
					/* calc this cross-spectrum magnitude */
					arm_cmplx_mag_q31(tmp_cross_spec, &tmp_cross_spec_abs, 1U);

					/* for numerical stability - if it's zero or near-zero, skip */
					if (tmp_cross_spec_abs > EPS) {
						/* perform PHAT part - normalizacja */
						arm_divide_q31(tmp_cross_spec[0], tmp_cross_spec_abs, &tmp_cross_spec[0], NULL);
						arm_divide_q31(tmp_cross_spec[1], tmp_cross_spec_abs, &tmp_cross_spec[1], NULL);

						/* oblicz tau_ml dla tej pary mikrofonów i tego DoA */
						tau_ml = __QSUB(all_relative_tdoas_lut[doa_idx + mic_l].tau,
										all_relative_tdoas_lut[doa_idx + mic_m].tau);

						/* Oblicz fazę: omega * tau */
						omega = 2.0f * PI * fft_bins_lut[k];
						arm_q31_to_float(&tau_ml, &tau_float, 1U);
						phase = omega * tau_float;

						/* Oblicz e^(j*phase) w Q31 */
						arm_float_to_q31(&phase, &phase_q31, 1U);
						exp_phase[0] = arm_cos_q31(phase_q31);  // Re
						exp_phase[1] = arm_sin_q31(phase_q31);  // Im

						/* Mnożenie: tmp_cross_spec * exp_phase */
						arm_cmplx_mult_cmplx_q31(tmp_cross_spec, exp_phase, result, 1U);

						/* Akumuluj wynik (sumuj po częstotliwościach) */
						gcc_phat_sum[0] = __QADD(gcc_phat_sum[0], result[0]);
						gcc_phat_sum[1] = __QADD(gcc_phat_sum[1], result[1]);
					}
				}

				/* Po zsumowaniu po częstotliwościach, oblicz moc (|gcc_phat_sum|^2) */
				arm_power_q31(gcc_phat_sum, 2U, &tmp_gcc_pow);

				/* Dodaj do mapy SRP dla tego DoA (konwersja Q63 -> float) */
				srp_phat_mic_per_doa[doa] += (float)((double)tmp_gcc_pow / ((double)(1LL << 63)));
			}
		}
	}

	/* -- 4. return maximum from computed map -- */
	float tmp_dummy = 0.0f;
	arm_max_f32(srp_phat_mic_per_doa, (uint32_t) ALL_POSSIBLE_DOAS, &tmp_dummy, &doa_idx);
	ret_doa = all_possible_doas_lut[doa_idx];

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

  #if L == 512
    (void) arm_rfft_init_512_q31(&rfft, 0U, 0U);
  #else
	(void) arm_rfft_init_q31(&rfft, (uint32_t) L, 0U, 0U);
  #endif

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
		// send doa thru UART
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
