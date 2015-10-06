/* 
 * File: _coder_otg_xy_reallyDumb_api.h 
 *  
 * MATLAB Coder version            : 2.7 
 * C/C++ source code generated on  : 06-Oct-2015 11:15:45 
 */

#ifndef ___CODER_OTG_XY_REALLYDUMB_API_H__
#define ___CODER_OTG_XY_REALLYDUMB_API_H__
/* Include Files */ 
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"

/* Type Definitions */ 
#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T
struct emxArray_real_T
{
    real_T *data;
    int32_T *size;
    int32_T allocatedSize;
    int32_T numDimensions;
    boolean_T canFreeData;
};
#endif /*struct_emxArray_real_T*/
#ifndef typedef_emxArray_real_T
#define typedef_emxArray_real_T
typedef struct emxArray_real_T emxArray_real_T;
#endif /*typedef_emxArray_real_T*/

/* Function Declarations */ 
extern void otg_xy_reallyDumb_initialize(emlrtContext *aContext);
extern void otg_xy_reallyDumb_terminate(void);
extern void otg_xy_reallyDumb_atexit(void);
extern void otg_xy_reallyDumb_api(const mxArray *prhs[18], const mxArray *plhs[4]);
extern void otg_xy_reallyDumb(real_T S[3], real_T D[4], real_T kj, real_T kT, real_T ks, real_T kd, real_T dT, real_T Tmin, real_T Tmax, emxArray_real_T *dataVeh, real_T safetyS, real_T safetyD, real_T dt, real_T ds, real_T m, real_T kappa, real_T b_y0, real_T phi, real_T *flag1, real_T *flag2, emxArray_real_T *x, emxArray_real_T *y);
extern void otg_xy_reallyDumb_xil_terminate(void);

#endif
/* 
 * File trailer for _coder_otg_xy_reallyDumb_api.h 
 *  
 * [EOF] 
 */
