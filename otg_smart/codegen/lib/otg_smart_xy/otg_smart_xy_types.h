/* 
 * File: otg_smart_xy_types.h 
 *  
 * MATLAB Coder version            : 2.7 
 * C/C++ source code generated on  : 07-Oct-2015 17:18:09 
 */

#ifndef __OTG_SMART_XY_TYPES_H__
#define __OTG_SMART_XY_TYPES_H__

/* Include Files */ 
#include "rtwtypes.h"

/* Type Definitions */ 
#ifndef struct_emxArray__common
#define struct_emxArray__common
struct emxArray__common
{
    void *data;
    int *size;
    int allocatedSize;
    int numDimensions;
    boolean_T canFreeData;
};
#endif /*struct_emxArray__common*/
#ifndef typedef_emxArray__common
#define typedef_emxArray__common
typedef struct emxArray__common emxArray__common;
#endif /*typedef_emxArray__common*/
#ifndef struct_emxArray_boolean_T_1x2
#define struct_emxArray_boolean_T_1x2
struct emxArray_boolean_T_1x2
{
    boolean_T data[2];
    int size[2];
};
#endif /*struct_emxArray_boolean_T_1x2*/
#ifndef typedef_emxArray_boolean_T_1x2
#define typedef_emxArray_boolean_T_1x2
typedef struct emxArray_boolean_T_1x2 emxArray_boolean_T_1x2;
#endif /*typedef_emxArray_boolean_T_1x2*/
#ifndef struct_emxArray_boolean_T_1x3
#define struct_emxArray_boolean_T_1x3
struct emxArray_boolean_T_1x3
{
    boolean_T data[3];
    int size[2];
};
#endif /*struct_emxArray_boolean_T_1x3*/
#ifndef typedef_emxArray_boolean_T_1x3
#define typedef_emxArray_boolean_T_1x3
typedef struct emxArray_boolean_T_1x3 emxArray_boolean_T_1x3;
#endif /*typedef_emxArray_boolean_T_1x3*/
#ifndef struct_emxArray_boolean_T_4
#define struct_emxArray_boolean_T_4
struct emxArray_boolean_T_4
{
    boolean_T data[4];
    int size[1];
};
#endif /*struct_emxArray_boolean_T_4*/
#ifndef typedef_emxArray_boolean_T_4
#define typedef_emxArray_boolean_T_4
typedef struct emxArray_boolean_T_4 emxArray_boolean_T_4;
#endif /*typedef_emxArray_boolean_T_4*/
#ifndef struct_emxArray_creal_T_1x3
#define struct_emxArray_creal_T_1x3
struct emxArray_creal_T_1x3
{
    creal_T data[3];
    int size[2];
};
#endif /*struct_emxArray_creal_T_1x3*/
#ifndef typedef_emxArray_creal_T_1x3
#define typedef_emxArray_creal_T_1x3
typedef struct emxArray_creal_T_1x3 emxArray_creal_T_1x3;
#endif /*typedef_emxArray_creal_T_1x3*/
#ifndef struct_emxArray_creal_T_1x4
#define struct_emxArray_creal_T_1x4
struct emxArray_creal_T_1x4
{
    creal_T data[4];
    int size[2];
};
#endif /*struct_emxArray_creal_T_1x4*/
#ifndef typedef_emxArray_creal_T_1x4
#define typedef_emxArray_creal_T_1x4
typedef struct emxArray_creal_T_1x4 emxArray_creal_T_1x4;
#endif /*typedef_emxArray_creal_T_1x4*/
#ifndef struct_emxArray_creal_T_1x5
#define struct_emxArray_creal_T_1x5
struct emxArray_creal_T_1x5
{
    creal_T data[5];
    int size[2];
};
#endif /*struct_emxArray_creal_T_1x5*/
#ifndef typedef_emxArray_creal_T_1x5
#define typedef_emxArray_creal_T_1x5
typedef struct emxArray_creal_T_1x5 emxArray_creal_T_1x5;
#endif /*typedef_emxArray_creal_T_1x5*/
#ifndef struct_emxArray_creal_T_4
#define struct_emxArray_creal_T_4
struct emxArray_creal_T_4
{
    creal_T data[4];
    int size[1];
};
#endif /*struct_emxArray_creal_T_4*/
#ifndef typedef_emxArray_creal_T_4
#define typedef_emxArray_creal_T_4
typedef struct emxArray_creal_T_4 emxArray_creal_T_4;
#endif /*typedef_emxArray_creal_T_4*/
#ifndef struct_emxArray_int32_T_1x3
#define struct_emxArray_int32_T_1x3
struct emxArray_int32_T_1x3
{
    int data[3];
    int size[2];
};
#endif /*struct_emxArray_int32_T_1x3*/
#ifndef typedef_emxArray_int32_T_1x3
#define typedef_emxArray_int32_T_1x3
typedef struct emxArray_int32_T_1x3 emxArray_int32_T_1x3;
#endif /*typedef_emxArray_int32_T_1x3*/
#ifndef struct_emxArray_int32_T_4
#define struct_emxArray_int32_T_4
struct emxArray_int32_T_4
{
    int data[4];
    int size[1];
};
#endif /*struct_emxArray_int32_T_4*/
#ifndef typedef_emxArray_int32_T_4
#define typedef_emxArray_int32_T_4
typedef struct emxArray_int32_T_4 emxArray_int32_T_4;
#endif /*typedef_emxArray_int32_T_4*/
#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T
struct emxArray_real_T
{
    double *data;
    int *size;
    int allocatedSize;
    int numDimensions;
    boolean_T canFreeData;
};
#endif /*struct_emxArray_real_T*/
#ifndef typedef_emxArray_real_T
#define typedef_emxArray_real_T
typedef struct emxArray_real_T emxArray_real_T;
#endif /*typedef_emxArray_real_T*/
#ifndef struct_emxArray_real_T_1x3
#define struct_emxArray_real_T_1x3
struct emxArray_real_T_1x3
{
    double data[3];
    int size[2];
};
#endif /*struct_emxArray_real_T_1x3*/
#ifndef typedef_emxArray_real_T_1x3
#define typedef_emxArray_real_T_1x3
typedef struct emxArray_real_T_1x3 emxArray_real_T_1x3;
#endif /*typedef_emxArray_real_T_1x3*/
#ifndef struct_emxArray_real_T_1x4
#define struct_emxArray_real_T_1x4
struct emxArray_real_T_1x4
{
    double data[4];
    int size[2];
};
#endif /*struct_emxArray_real_T_1x4*/
#ifndef typedef_emxArray_real_T_1x4
#define typedef_emxArray_real_T_1x4
typedef struct emxArray_real_T_1x4 emxArray_real_T_1x4;
#endif /*typedef_emxArray_real_T_1x4*/
#ifndef struct_emxArray_real_T_1x5
#define struct_emxArray_real_T_1x5
struct emxArray_real_T_1x5
{
    double data[5];
    int size[2];
};
#endif /*struct_emxArray_real_T_1x5*/
#ifndef typedef_emxArray_real_T_1x5
#define typedef_emxArray_real_T_1x5
typedef struct emxArray_real_T_1x5 emxArray_real_T_1x5;
#endif /*typedef_emxArray_real_T_1x5*/
#ifndef struct_emxArray_real_T_3
#define struct_emxArray_real_T_3
struct emxArray_real_T_3
{
    double data[3];
    int size[1];
};
#endif /*struct_emxArray_real_T_3*/
#ifndef typedef_emxArray_real_T_3
#define typedef_emxArray_real_T_3
typedef struct emxArray_real_T_3 emxArray_real_T_3;
#endif /*typedef_emxArray_real_T_3*/

#endif
/* 
 * File trailer for otg_smart_xy_types.h 
 *  
 * [EOF] 
 */
