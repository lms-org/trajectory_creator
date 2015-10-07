/* 
 * File: _coder_otg_smart_xy_api.h 
 *  
 * MATLAB Coder version            : 2.7 
 * C/C++ source code generated on  : 07-Oct-2015 17:18:09 
 */

#ifndef ___CODER_OTG_SMART_XY_API_H__
#define ___CODER_OTG_SMART_XY_API_H__
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
extern void otg_smart_xy_initialize(emlrtContext *aContext);
extern void otg_smart_xy_terminate(void);
extern void otg_smart_xy_atexit(void);
extern void otg_smart_xy_api(const mxArray *prhs[20], const mxArray *plhs[8]);
extern void otg_smart_xy(real_T absTOL, int16_T maxIter, real_T v1, real_T d1, real_T kj, real_T kT, real_T ks, real_T kd, real_T dataVeh[3], real_T safetyS, real_T safetyD, real_T kappaMax, real_T aOrthMax, int16_T m, real_T kappa, real_T b_y0, real_T phi, real_T vx0, real_T ax0, real_T w, real_T *flag1, real_T *flag2, real_T *flag3, real_T *flagAll, emxArray_real_T *x, emxArray_real_T *y, real_T T_data[], int32_T T_size[2], real_T *TOL);
extern void otg_smart_xy_xil_terminate(void);

#endif
/* 
 * File trailer for _coder_otg_smart_xy_api.h 
 *  
 * [EOF] 
 */
