/***************************************************************************//**
 * @file
 * @brief A very simple demonstration of different power modes.
 * @version 5.6.1
 *******************************************************************************
 * # License
 * <b>Copyright 2015 Silicon Labs, Inc. http://www.silabs.com</b>
 *******************************************************************************
 *
 * This file is licensed under the Silabs License Agreement. See the file
 * "Silabs_License_Agreement.txt" for details. Before using this software for
 * any purpose, you must agree to the terms of that agreement.
 *
 ******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include "em_device.h"
#include "em_chip.h"
#include "em_cmu.h"
#include "em_emu.h"
#include "em_wdog.h"
#include "rtcdriver.h"
#include "bsp_trace.h"
#include "matrix_ops.h"
#include "qp.h"
#include "qp_solvers.h"

#define TMP_FILENAME "tmp_test_file"
#define SHELL_REF_TEST_CMD "python qp_ref.py " TMP_FILENAME
#define TEST_PRECISION 1e-6
#define NUM_INVERSION_TEST_RUNS 8

#define P_RAND_ENTRY_MIN -1e4
#define P_RAND_ENTRY_MAX 1e4
#define Q_RAND_ENTRY_MIN -1e3
#define Q_RAND_ENTRY_MAX 1e3
#define X_RAND_ENTRY_MIN -1e3
#define X_RAND_ENTRY_MAX 1e3

#define NUM_OPT_TEST_RUNS 8

#define ITERATIONS 1e4

#define ABS(x) (x < 0 ? -x : x)
/** Counts 1ms timeTicks */
volatile uint32_t msTicks;

/** Timer used for bringing the system back to EM0. */
RTCDRV_TimerID_t xTimerForWakeUp;

/***************************************************************************//**
 * @brief SysTick_Handler
 * Interrupt Service Routine for system tick counter
 ******************************************************************************/
void SysTick_Handler(void)
{
  msTicks++;       /* increment counter necessary in Delay()*/
}

/***************************************************************************//**
 * @brief SysTick_Disable
 * Disable systick interrupts
 ******************************************************************************/
void SysTick_Disable(void)
{
  SysTick->CTRL = 0x0000000;
}

/***************************************************************************//**
 * @brief Delays number of msTick Systicks (typically 1 ms)
 * @param dlyTicks Number of ticks to delay
 ******************************************************************************/
void Delay(uint32_t dlyTicks)
{
  uint32_t curTicks;

  curTicks = msTicks;
  while ((msTicks - curTicks) < dlyTicks) ;
}

/***************************************************************************//**
 * @brief  Main function
 ******************************************************************************/
int main(void)
{
  //WDOG_Init_TypeDef wInit = WDOG_INIT_DEFAULT;
  //int i;

  /* Chip revision alignment and errata fixes */
  CHIP_Init();

  /* If first word of user data page is non-zero, enable Energy Profiler trace */
  BSP_TraceProfilerSetup();

  /* Initialize RTC timer. */
  RTCDRV_Init();
  RTCDRV_AllocateTimer(&xTimerForWakeUp);

  /* Watchdog setup - Use defaults, excepts for these :*/
  //wInit.em2Run = true;
  //wInit.em3Run = true;
  //wInit.perSel = wdogPeriod_4k; /* 4k 1kHz periods should give ~4 seconds in EM3 */

  /* Do the demo forever. */

while(1)
{
	RTCDRV_StartTimer(xTimerForWakeUp, rtcdrvTimerTypeOneshot, 10000, NULL, NULL);

  /* EM0 - 1 sec HFRCO */
  CMU_ClockSelectSet(cmuClock_HF, cmuSelect_HFRCO);
  /* Setup SysTick Timer for 1 msec interrupts  */
  if (SysTick_Config(CMU_ClockFreqGet(cmuClock_CORE) / 1000)) {
    while (1) ;
  }
  Delay(5000);

  /* Turn off systick */
  SysTick_Disable();

  /* EM1 - 1 sec */
  //RTCDRV_StartTimer(xTimerForWakeUp, rtcdrvTimerTypeOneshot, 1000, NULL, NULL);
  //EMU_EnterEM1();

  /* EM2 - 1 sec */

  //EMU_EnterEM2(true);
  EMU_EnterEM3(true);

}

  /* Start watchdog */
  //WDOG_Init(&wInit);

  /* Enter EM3 - Watchdog will reset chip (and confuse debuggers) */
  //EMU_EnterEM3(true);

  /* We will never reach this point */
  return 0;
}
