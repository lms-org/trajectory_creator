/*
 * File: _coder_otg_smart_xy_api.c
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 23-Oct-2015 14:10:50
 */

/* Include Files */
#include "_coder_otg_smart_xy_api.h"

/* Function Declarations */
static real_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *absTOL,
  const char_T *identifier);
static real_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static int32_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *maxIter,
  const char_T *identifier);
static int32_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *dataVeh,
  const char_T *identifier))[3];
static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[3];
static const mxArray *emlrt_marshallOut(const real_T u);
static const mxArray *b_emlrt_marshallOut(const emxArray_real_T *u);
static const mxArray *c_emlrt_marshallOut(const real_T u_data[], const int32_T
  u_size[2]);
static real_T g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);
static int32_T h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId);
static real_T (*i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[3];
static void emxInit_real_T(const emlrtStack *sp, emxArray_real_T **pEmxArray,
  int32_T numDimensions, boolean_T doPush);
static void emxFree_real_T(emxArray_real_T **pEmxArray);

/* Function Definitions */

/*
 * Arguments    : emlrtContext *aContext
 * Return Type  : void
 */
void otg_smart_xy_initialize(emlrtContext *aContext)
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
void otg_smart_xy_terminate(void)
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
void otg_smart_xy_atexit(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  otg_smart_xy_xil_terminate();
}

/*
 * Arguments    : const mxArray *prhs[20]
 *                const mxArray *plhs[8]
 * Return Type  : void
 */
void otg_smart_xy_api(const mxArray *prhs[20], const mxArray *plhs[8])
{
  real_T (*T_data)[3];
  emxArray_real_T *x;
  emxArray_real_T *y;
  real_T absTOL;
  int32_T maxIter;
  real_T v1;
  real_T d1;
  real_T kj;
  real_T kT;
  real_T ks;
  real_T kd;
  real_T (*dataVeh)[3];
  real_T safetyS;
  real_T safetyD;
  real_T kappaMax;
  real_T aOrthMax;
  int32_T m;
  real_T kappa;
  real_T b_y0;
  real_T phi;
  real_T vx0;
  real_T ax0;
  real_T w;
  real_T TOL;
  int32_T T_size[2];
  real_T flagAll;
  real_T flag3;
  real_T flag2;
  real_T flag1;
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  T_data = (real_T (*)[3])mxMalloc(sizeof(real_T [3]));
  emlrtHeapReferenceStackEnterFcnR2012b(&st);
  emxInit_real_T(&st, &x, 2, true);
  emxInit_real_T(&st, &y, 2, true);
  prhs[8] = emlrtProtectR2012b(prhs[8], 8, false, -1);

  /* Marshall function inputs */
  absTOL = emlrt_marshallIn(&st, emlrtAliasP(prhs[0]), "absTOL");
  maxIter = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "maxIter");
  v1 = emlrt_marshallIn(&st, emlrtAliasP(prhs[2]), "v1");
  d1 = emlrt_marshallIn(&st, emlrtAliasP(prhs[3]), "d1");
  kj = emlrt_marshallIn(&st, emlrtAliasP(prhs[4]), "kj");
  kT = emlrt_marshallIn(&st, emlrtAliasP(prhs[5]), "kT");
  ks = emlrt_marshallIn(&st, emlrtAliasP(prhs[6]), "ks");
  kd = emlrt_marshallIn(&st, emlrtAliasP(prhs[7]), "kd");
  dataVeh = e_emlrt_marshallIn(&st, emlrtAlias(prhs[8]), "dataVeh");
  safetyS = emlrt_marshallIn(&st, emlrtAliasP(prhs[9]), "safetyS");
  safetyD = emlrt_marshallIn(&st, emlrtAliasP(prhs[10]), "safetyD");
  kappaMax = emlrt_marshallIn(&st, emlrtAliasP(prhs[11]), "kappaMax");
  aOrthMax = emlrt_marshallIn(&st, emlrtAliasP(prhs[12]), "aOrthMax");
  m = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[13]), "m");
  kappa = emlrt_marshallIn(&st, emlrtAliasP(prhs[14]), "kappa");
  b_y0 = emlrt_marshallIn(&st, emlrtAliasP(prhs[15]), "y0");
  phi = emlrt_marshallIn(&st, emlrtAliasP(prhs[16]), "phi");
  vx0 = emlrt_marshallIn(&st, emlrtAliasP(prhs[17]), "vx0");
  ax0 = emlrt_marshallIn(&st, emlrtAliasP(prhs[18]), "ax0");
  w = emlrt_marshallIn(&st, emlrtAliasP(prhs[19]), "w");

  /* Invoke the target function */
  otg_smart_xy(absTOL, maxIter, v1, d1, kj, kT, ks, kd, *dataVeh, safetyS,
               safetyD, kappaMax, aOrthMax, m, kappa, b_y0, phi, vx0, ax0, w,
               &flag1, &flag2, &flag3, &flagAll, x, y, *T_data, T_size, &TOL);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(flag1);
  plhs[1] = emlrt_marshallOut(flag2);
  plhs[2] = emlrt_marshallOut(flag3);
  plhs[3] = emlrt_marshallOut(flagAll);
  plhs[4] = b_emlrt_marshallOut(x);
  plhs[5] = b_emlrt_marshallOut(y);
  plhs[6] = c_emlrt_marshallOut(*T_data, T_size);
  plhs[7] = emlrt_marshallOut(TOL);
  y->canFreeData = false;
  emxFree_real_T(&y);
  x->canFreeData = false;
  emxFree_real_T(&x);
  emlrtHeapReferenceStackLeaveFcnR2012b(&st);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *absTOL
 *                const char_T *identifier
 * Return Type  : real_T
 */
static real_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *absTOL,
  const char_T *identifier)
{
  real_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  y = b_emlrt_marshallIn(sp, emlrtAlias(absTOL), &thisId);
  emlrtDestroyArray(&absTOL);
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
  y = g_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *maxIter
 *                const char_T *identifier
 * Return Type  : int32_T
 */
static int32_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *maxIter,
  const char_T *identifier)
{
  int32_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  y = d_emlrt_marshallIn(sp, emlrtAlias(maxIter), &thisId);
  emlrtDestroyArray(&maxIter);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : int32_T
 */
static int32_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  int32_T y;
  y = h_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *dataVeh
 *                const char_T *identifier
 * Return Type  : real_T (*)[3]
 */
static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *dataVeh,
  const char_T *identifier))[3]
{
  real_T (*y)[3];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  y = f_emlrt_marshallIn(sp, emlrtAlias(dataVeh), &thisId);
  emlrtDestroyArray(&dataVeh);
  return y;
}
/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real_T (*)[3]
 */
  static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
  const emlrtMsgIdentifier *parentId))[3]
{
  real_T (*y)[3];
  y = i_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
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
 * Arguments    : const real_T u_data[]
 *                const int32_T u_size[2]
 * Return Type  : const mxArray *
 */
static const mxArray *c_emlrt_marshallOut(const real_T u_data[], const int32_T
  u_size[2])
{
  const mxArray *y;
  static const int32_T iv1[2] = { 0, 0 };

  const mxArray *m2;
  y = NULL;
  m2 = emlrtCreateNumericArray(2, iv1, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m2, (void *)u_data);
  emlrtSetDimensions((mxArray *)m2, u_size, 2);
  emlrtAssign(&y, m2);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real_T
 */
static real_T g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
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
 * Return Type  : int32_T
 */
static int32_T h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId)
{
  int32_T ret;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "int32", false, 0U, 0);
  ret = *(int32_T *)mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real_T (*)[3]
 */
static real_T (*i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[3]
{
  real_T (*ret)[3];
  int32_T iv2[1];
  iv2[0] = 3;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 1U, iv2);
  ret = (real_T (*)[3])mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
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
 * File trailer for _coder_otg_smart_xy_api.c
 *
 * [EOF]
 */
