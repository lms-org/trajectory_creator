/*
 * File: repmat.c
 *
 * MATLAB Coder version            : 2.7
 * C/C++ source code generated on  : 23-Oct-2015 13:28:58
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "otg_smart_xy.h"
#include "repmat.h"
#include "otg_smart_xy_emxutil.h"

/* Function Definitions */

/*
 * Arguments    : const emxArray_real_T *a
 *                emxArray_real_T *b
 * Return Type  : void
 */
void repmat(const emxArray_real_T *a, emxArray_real_T *b)
{
  int unnamed_idx_1;
  int ibmat;
  int itilerow;
  unnamed_idx_1 = b->size[0] * b->size[1];
  b->size[0] = 2;
  b->size[1] = a->size[1];
  emxEnsureCapacity((emxArray__common *)b, unnamed_idx_1, (int)sizeof(double));
  unnamed_idx_1 = a->size[1];
  if (unnamed_idx_1 == 0) {
  } else {
    for (unnamed_idx_1 = 0; unnamed_idx_1 + 1 <= a->size[1]; unnamed_idx_1++) {
      ibmat = unnamed_idx_1 << 1;
      for (itilerow = 0; itilerow < 2; itilerow++) {
        b->data[ibmat + itilerow] = a->data[unnamed_idx_1];
      }
    }
  }
}

/*
 * File trailer for repmat.c
 *
 * [EOF]
 */
