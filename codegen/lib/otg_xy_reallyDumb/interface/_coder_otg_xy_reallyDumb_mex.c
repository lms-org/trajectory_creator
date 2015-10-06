/*
 * _coder_otg_xy_reallyDumb_mex.c
 *
 * Code generation for function 'otg_xy_reallyDumb'
 *
 */

/* Include files */
#include "mex.h"
#include "_coder_otg_xy_reallyDumb_api.h"

/* Function Declarations */
static void otg_xy_reallyDumb_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

/* Variable Definitions */
emlrtContext emlrtContextGlobal = { true, false, EMLRT_VERSION_INFO, NULL, "otg_xy_reallyDumb", NULL, false, {2045744189U,2170104910U,2743257031U,4284093946U}, NULL };
void *emlrtRootTLSGlobal = NULL;

/* Function Definitions */
static void otg_xy_reallyDumb_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  const mxArray *outputs[4];
  const mxArray *inputs[18];
  int n = 0;
  int nOutputs = (nlhs < 1 ? 1 : nlhs);
  int nInputs = nrhs;
  emlrtStack st = { NULL, NULL, NULL };
  /* Module initialization. */
  otg_xy_reallyDumb_initialize(&emlrtContextGlobal);
  st.tls = emlrtRootTLSGlobal;
  /* Check for proper number of arguments. */
  if (nrhs != 18) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, mxINT32_CLASS, 18, mxCHAR_CLASS, 17, "otg_xy_reallyDumb");
  } else if (nlhs > 4) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, mxCHAR_CLASS, 17, "otg_xy_reallyDumb");
  }
  /* Temporary copy for mex inputs. */
  for (n = 0; n < nInputs; ++n) {
    inputs[n] = prhs[n];
  }
  /* Call the function. */
  otg_xy_reallyDumb_api(inputs, outputs);
  /* Copy over outputs to the caller. */
  for (n = 0; n < nOutputs; ++n) {
    plhs[n] = emlrtReturnArrayR2009a(outputs[n]);
  }
  /* Module finalization. */
  otg_xy_reallyDumb_terminate();
}

void otg_xy_reallyDumb_atexit_wrapper(void)
{
   otg_xy_reallyDumb_atexit();
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* Initialize the memory manager. */
  mexAtExit(otg_xy_reallyDumb_atexit_wrapper);
  /* Dispatch the entry-point. */
  otg_xy_reallyDumb_mexFunction(nlhs, plhs, nrhs, prhs);
}
/* End of code generation (_coder_otg_xy_reallyDumb_mex.c) */
