/*
 * File: _coder_otg_xy_reallyDumb_api.c
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 06-Oct-2015 16:22:34
 */

/* Include Files */
#include "_coder_otg_xy_reallyDumb_api.h"

/* Function Declarations */
static real_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *v1, const
  char_T *identifier);
static real_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *dataVeh,
  const char_T *identifier, emxArray_real_T *y);
static void d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y);
static const mxArray *emlrt_marshallOut(const real_T u);
static const mxArray *b_emlrt_marshallOut(const emxArray_real_T *u);
static real_T e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);
static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret);
static void emxInit_real_T(const emlrtStack *sp, emxArray_real_T **pEmxArray,
  int32_T numDimensions, boolean_T doPush);
static void emxFree_real_T(emxArray_real_T **pEmxArray);

/* Function Definitions */

/*
 * Arguments    : emlrtContext *aContext
 * Return Type  : void
 */
void otg_xy_reallyDumb_initialize(emlrtContext *aContext)
{
  emlrtStack st = { NULL, NULL, NULL };

  emlrtCreateRootTLS(&emlrtRootTLSGlobal, aContext, NULL, 1);
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void otg_xy_reallyDumb_terminate(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void otg_xy_reallyDumb_atexit(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  otg_xy_reallyDumb_xil_terminate();
}

/*
 * Arguments    : const mxArray *prhs[20]
 *                const mxArray *plhs[4]
 * Return Type  : void
 */
void otg_xy_reallyDumb_api(const mxArray *prhs[20], const mxArray *plhs[4])
{
  emxArray_real_T *dataVeh;
  emxArray_real_T *x;
  emxArray_real_T *y;
  real_T v1;
  real_T d1;
  real_T kj;
  real_T kT;
  real_T ks;
  real_T kd;
  real_T dT;
  real_T Tmin;
  real_T Tmax;
  real_T safetyS;
  real_T safetyD;
  real_T dt;
  real_T m;
  real_T kappa;
  real_T b_y0;
  real_T phi;
  real_T vx0;
  real_T ax0;
  real_T w;
  real_T T;
  real_T flag;
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  emlrtHeapReferenceStackEnterFcnR2012b(&st);
  emxInit_real_T(&st, &dataVeh, 2, true);
  emxInit_real_T(&st, &x, 2, true);
  emxInit_real_T(&st, &y, 2, true);
  prhs[9] = emlrtProtectR2012b(prhs[9], 9, false, -1);

  /* Marshall function inputs */
  v1 = emlrt_marshallIn(&st, emlrtAliasP(prhs[0]), "v1");
  d1 = emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "d1");
  kj = emlrt_marshallIn(&st, emlrtAliasP(prhs[2]), "kj");
  kT = emlrt_marshallIn(&st, emlrtAliasP(prhs[3]), "kT");
  ks = emlrt_marshallIn(&st, emlrtAliasP(prhs[4]), "ks");
  kd = emlrt_marshallIn(&st, emlrtAliasP(prhs[5]), "kd");
  dT = emlrt_marshallIn(&st, emlrtAliasP(prhs[6]), "dT");
  Tmin = emlrt_marshallIn(&st, emlrtAliasP(prhs[7]), "Tmin");
  Tmax = emlrt_marshallIn(&st, emlrtAliasP(prhs[8]), "Tmax");
  c_emlrt_marshallIn(&st, emlrtAlias(prhs[9]), "dataVeh", dataVeh);
  safetyS = emlrt_marshallIn(&st, emlrtAliasP(prhs[10]), "safetyS");
  safetyD = emlrt_marshallIn(&st, emlrtAliasP(prhs[11]), "safetyD");
  dt = emlrt_marshallIn(&st, emlrtAliasP(prhs[12]), "dt");
  m = emlrt_marshallIn(&st, emlrtAliasP(prhs[13]), "m");
  kappa = emlrt_marshallIn(&st, emlrtAliasP(prhs[14]), "kappa");
  b_y0 = emlrt_marshallIn(&st, emlrtAliasP(prhs[15]), "y0");
  phi = emlrt_marshallIn(&st, emlrtAliasP(prhs[16]), "phi");
  vx0 = emlrt_marshallIn(&st, emlrtAliasP(prhs[17]), "vx0");
  ax0 = emlrt_marshallIn(&st, emlrtAliasP(prhs[18]), "ax0");
  w = emlrt_marshallIn(&st, emlrtAliasP(prhs[19]), "w");

  /* Invoke the target function */
  otg_xy_reallyDumb(v1, d1, kj, kT, ks, kd, dT, Tmin, Tmax, dataVeh, safetyS,
                    safetyD, dt, m, kappa, b_y0, phi, vx0, ax0, w, &flag, x, y,
                    &T);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(flag);
  plhs[1] = b_emlrt_marshallOut(x);
  plhs[2] = b_emlrt_marshallOut(y);
  plhs[3] = emlrt_marshallOut(T);
  y->canFreeData = false;
  emxFree_real_T(&y);
  x->canFreeData = false;
  emxFree_real_T(&x);
  dataVeh->canFreeData = false;
  emxFree_real_T(&dataVeh);
  emlrtHeapReferenceStackLeaveFcnR2012b(&st);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *v1
 *                const char_T *identifier
 * Return Type  : real_T
 */
static real_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *v1, const
  char_T *identifier)
{
  real_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  y = b_emlrt_marshallIn(sp, emlrtAlias(v1), &thisId);
  emlrtDestroyArray(&v1);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real_T
 */
static real_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = e_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *dataVeh
 *                const char_T *identifier
 *                emxArray_real_T *y
 * Return Type  : void
 */
static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *dataVeh,
  const char_T *identifier, emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  d_emlrt_marshallIn(sp, emlrtAlias(dataVeh), &thisId, y);
  emlrtDestroyArray(&dataVeh);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                emxArray_real_T *y
 * Return Type  : void
 */
static void d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y)
{
  f_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const real_T u
 * Return Type  : const mxArray *
 */
static const mxArray *emlrt_marshallOut(const real_T u)
{
  const mxArray *y;
  const mxArray *m0;
  y = NULL;
  m0 = emlrtCreateDoubleScalar(u);
  emlrtAssign(&y, m0);
  return y;
}

/*
 * Arguments    : const emxArray_real_T *u
 * Return Type  : const mxArray *
 */
static const mxArray *b_emlrt_marshallOut(const emxArray_real_T *u)
{
  const mxArray *y;
  static const int32_T iv0[2] = { 0, 0 };

  const mxArray *m1;
  y = NULL;
  m1 = emlrtCreateNumericArray(2, iv0, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m1, (void *)u->data);
  emlrtSetDimensions((mxArray *)m1, u->size, 2);
  emlrtAssign(&y, m1);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real_T
 */
static real_T e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId)
{
  real_T ret;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 0U, 0);
  ret = *(real_T *)mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 *                emxArray_real_T *ret
 * Return Type  : void
 */
static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret)
{
  int32_T iv1[2];
  boolean_T bv0[2];
  int32_T i0;
  static const boolean_T bv1[2] = { false, true };

  int32_T iv2[2];
  for (i0 = 0; i0 < 2; i0++) {
    iv1[i0] = 3 + -4 * i0;
    bv0[i0] = bv1[i0];
  }

  emlrtCheckVsBuiltInR2012b(sp, msgId, src, "double", false, 2U, iv1, bv0, iv2);
  ret->size[0] = iv2[0];
  ret->size[1] = iv2[1];
  ret->allocatedSize = ret->size[0] * ret->size[1];
  ret->data = (real_T *)mxGetData(src);
  ret->canFreeData = false;
  emlrtDestroyArray(&src);
}

/*
 * Arguments    : const emlrtStack *sp
 *                emxArray_real_T **pEmxArray
 *                int32_T numDimensions
 *                boolean_T doPush
 * Return Type  : void
 */
static void emxInit_real_T(const emlrtStack *sp, emxArray_real_T **pEmxArray,
  int32_T numDimensions, boolean_T doPush)
{
  emxArray_real_T *emxArray;
  int32_T i;
  *pEmxArray = (emxArray_real_T *)emlrtMallocMex(sizeof(emxArray_real_T));
  if (doPush) {
    emlrtPushHeapReferenceStackR2012b(sp, (void *)pEmxArray, (void (*)(void *))
      emxFree_real_T);
  }

  emxArray = *pEmxArray;
  emxArray->data = (real_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int32_T *)emlrtMallocMex((uint32_T)(sizeof(int32_T)
    * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

/*
 * Arguments    : emxArray_real_T **pEmxArray
 * Return Type  : void
 */
static void emxFree_real_T(emxArray_real_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_real_T *)NULL) {
    if (((*pEmxArray)->data != (real_T *)NULL) && (*pEmxArray)->canFreeData) {
      emlrtFreeMex((void *)(*pEmxArray)->data);
    }

    emlrtFreeMex((void *)(*pEmxArray)->size);
    emlrtFreeMex((void *)*pEmxArray);
    *pEmxArray = (emxArray_real_T *)NULL;
  }
}

/*
 * File trailer for _coder_otg_xy_reallyDumb_api.c
 *
 * [EOF]
 */
