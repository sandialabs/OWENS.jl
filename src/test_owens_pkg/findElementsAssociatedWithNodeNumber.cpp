//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: findElementsAssociatedWithNodeNumber.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 13-Apr-2020 09:25:21
//

// Include Files
#include "findElementsAssociatedWithNodeNumber.h"
#include "find.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include <cmath>
#include <math.h>
#include <string.h>

// Function Definitions

//
// findElementsAssociatedWithNodeNumber calculates reaction force at a node
//  **********************************************************************
//  *                   Part of the SNL OWENS Toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [elList,localNode] = findElementsAssociatedWithNodeNumber(nodeNum,...
//                         conn,jointData)
//
//    This function finds elements associated with a node number through
//    mesh connectivity or joint constraints
//
//    input:
//    nodeNum    = node number joint constraints are desired at
//    conn       = object containing mesh connectivity
//    jointData  = object containing joint informatino
// Arguments    : double nodeNum
//                const emxArray_real_T *conn
//                const emxArray_real_T *jointData
//                emxArray_real_T *elList
//                emxArray_real_T *localNode
// Return Type  : void
//
void c_findElementsAssociatedWithNod(double nodeNum, const emxArray_real_T *conn,
  const emxArray_real_T *jointData, emxArray_real_T *elList, emxArray_real_T
  *localNode)
{
  unsigned int b_index;
  int i;
  emxArray_boolean_T *tf;
  int exponent;
  int i1;
  boolean_T res[2];
  int idx;
  double absx;
  boolean_T b;
  boolean_T b1;
  int j;
  double b_r;
  emxArray_int32_T *res1;
  int b_exponent;
  int c_exponent;
  double d;
  boolean_T exitg1;
  signed char ii_data[2];
  double s;
  int d_exponent;

  //
  //    output:
  //    elList     = array containing a list of element numbers associated with
  //                 nodeNum
  //    localNode  = array containing the local node number that correspond to
  //                 nodeNum in the list of associated elements
  // search joint constraints
  b_index = 1U;

  // get number of elements in mesh
  i = elList->size[0] * elList->size[1];
  elList->size[0] = 1;
  elList->size[1] = 1;
  emxEnsureCapacity_real_T(elList, i);
  elList->data[0] = 0.0;
  i = localNode->size[0] * localNode->size[1];
  localNode->size[0] = 1;
  localNode->size[1] = 1;
  emxEnsureCapacity_real_T(localNode, i);
  localNode->data[0] = 0.0;
  if (jointData->size[0] != 0) {
    emxInit_boolean_T(&tf, 1);

    // first see if specified node is a slave node in a joint constraint
    i = jointData->size[0] - 1;
    i1 = tf->size[0];
    tf->size[0] = jointData->size[0];
    emxEnsureCapacity_boolean_T(tf, i1);
    idx = jointData->size[0];
    for (i1 = 0; i1 < idx; i1++) {
      tf->data[i1] = false;
    }

    for (j = 0; j <= i; j++) {
      frexp(0.5, &b_exponent);
      if (std::abs(1.0 - jointData->data[j + jointData->size[0] * 2]) <
          1.1102230246251565E-16) {
        tf->data[j] = true;
      }
    }

    emxInit_int32_T(&res1, 1);

    // search joint data slave nodes for node number
    // if it is, change it to the corresponding master node
    e_eml_find(tf, res1);
    if (res1->size[0] != 0) {
      nodeNum = jointData->data[(jointData->size[0] + jointData->size[0]) - 1];
    }

    i = jointData->size[0] - 1;
    i1 = tf->size[0];
    tf->size[0] = jointData->size[0];
    emxEnsureCapacity_boolean_T(tf, i1);
    idx = jointData->size[0];
    for (i1 = 0; i1 < idx; i1++) {
      tf->data[i1] = false;
    }

    for (j = 0; j <= i; j++) {
      absx = std::abs(nodeNum / 2.0);
      if ((!rtIsInf(absx)) && (!rtIsNaN(absx))) {
        if (absx <= 2.2250738585072014E-308) {
          b_r = 4.94065645841247E-324;
        } else {
          frexp(absx, &exponent);
          b_r = std::ldexp(1.0, exponent - 53);
        }
      } else {
        b_r = rtNaN;
      }

      d = jointData->data[j + jointData->size[0]];
      if ((std::abs(nodeNum - d) < b_r) || (rtIsInf(d) && rtIsInf(nodeNum) &&
           ((d > 0.0) == (nodeNum > 0.0)))) {
        tf->data[j] = true;
      }
    }

    e_eml_find(tf, res1);

    // search joint data master nodes for node number
    emxFree_boolean_T(&tf);
    if (res1->size[0] != 0) {
      i = elList->size[0] * elList->size[1];
      elList->size[0] = 1;
      elList->size[1] = static_cast<int>((static_cast<double>(res1->size[0]) *
        static_cast<double>(conn->size[0])));
      emxEnsureCapacity_real_T(elList, i);
      idx = static_cast<int>((static_cast<double>(res1->size[0]) * static_cast<
        double>(conn->size[0])));
      for (i = 0; i < idx; i++) {
        elList->data[i] = 0.0;
      }

      d = static_cast<double>(res1->size[0]) * static_cast<double>(conn->size[0]);
      i = static_cast<int>(d);
      i1 = localNode->size[0] * localNode->size[1];
      localNode->size[0] = i;
      localNode->size[1] = i;
      emxEnsureCapacity_real_T(localNode, i1);
      idx = i * i;
      for (i = 0; i < idx; i++) {
        localNode->data[i] = 0.0;
      }

      i = res1->size[0];
      for (j = 0; j < i; j++) {
        // loop over joints
        i1 = conn->size[0];
        for (exponent = 0; exponent < i1; exponent++) {
          s = jointData->data[(res1->data[j] + jointData->size[0] * 2) - 1];
          res[0] = false;
          res[1] = false;
          absx = std::abs(s / 2.0);
          b = !rtIsInf(absx);
          b1 = !rtIsNaN(absx);
          if (b && b1) {
            if (absx <= 2.2250738585072014E-308) {
              b_r = 4.94065645841247E-324;
            } else {
              frexp(absx, &d_exponent);
              b_r = std::ldexp(1.0, d_exponent - 53);
            }
          } else {
            b_r = rtNaN;
          }

          d = conn->data[exponent];
          if ((std::abs(s - d) < b_r) || (rtIsInf(d) && rtIsInf(s) && ((d > 0.0)
                == (s > 0.0)))) {
            res[0] = true;
          }

          if (b && b1) {
            if (absx <= 2.2250738585072014E-308) {
              b_r = 4.94065645841247E-324;
            } else {
              frexp(absx, &d_exponent);
              b_r = std::ldexp(1.0, d_exponent - 53);
            }
          } else {
            b_r = rtNaN;
          }

          d = conn->data[exponent + conn->size[0]];
          if ((std::abs(s - d) < b_r) || (rtIsInf(d) && rtIsInf(s) && ((d > 0.0)
                == (s > 0.0)))) {
            res[1] = true;
          }

          // finds indices of nodeNum in connectivity of element i
          idx = 0;
          b_exponent = 0;
          exitg1 = false;
          while ((!exitg1) && (b_exponent < 2)) {
            if (res[b_exponent]) {
              idx++;
              ii_data[idx - 1] = static_cast<signed char>((b_exponent + 1));
              if (idx >= 2) {
                exitg1 = true;
              } else {
                b_exponent++;
              }
            } else {
              b_exponent++;
            }
          }

          // finds the local node number of element i that corresponds to nodeNum 
          if (1 > idx) {
            idx = 0;
          }

          if (idx != 0) {
            // assigns to an elementList and localNode list
            idx = static_cast<int>(b_index) - 1;
            elList->data[idx] = static_cast<double>(exponent) + 1.0;
            localNode->data[idx] = ii_data[0];
            b_index++;
          }
        }
      }
    }

    emxFree_int32_T(&res1);
  }

  i = conn->size[0];
  for (exponent = 0; exponent < i; exponent++) {
    // loop over elements
    res[0] = false;
    res[1] = false;
    absx = std::abs(nodeNum / 2.0);
    b = !rtIsInf(absx);
    b1 = !rtIsNaN(absx);
    if (b && b1) {
      if (absx <= 2.2250738585072014E-308) {
        b_r = 4.94065645841247E-324;
      } else {
        frexp(absx, &c_exponent);
        b_r = std::ldexp(1.0, c_exponent - 53);
      }
    } else {
      b_r = rtNaN;
    }

    d = conn->data[exponent];
    if ((std::abs(nodeNum - d) < b_r) || (rtIsInf(d) && rtIsInf(nodeNum) && ((d >
           0.0) == (nodeNum > 0.0)))) {
      res[0] = true;
    }

    if (b && b1) {
      if (absx <= 2.2250738585072014E-308) {
        b_r = 4.94065645841247E-324;
      } else {
        frexp(absx, &c_exponent);
        b_r = std::ldexp(1.0, c_exponent - 53);
      }
    } else {
      b_r = rtNaN;
    }

    d = conn->data[exponent + conn->size[0]];
    if ((std::abs(nodeNum - d) < b_r) || (rtIsInf(d) && rtIsInf(nodeNum) && ((d >
           0.0) == (nodeNum > 0.0)))) {
      res[1] = true;
    }

    // finds indices of nodeNum in connectivity of element i
    idx = 0;
    b_exponent = 0;
    exitg1 = false;
    while ((!exitg1) && (b_exponent < 2)) {
      if (res[b_exponent]) {
        idx++;
        ii_data[idx - 1] = static_cast<signed char>((b_exponent + 1));
        if (idx >= 2) {
          exitg1 = true;
        } else {
          b_exponent++;
        }
      } else {
        b_exponent++;
      }
    }

    // finds the local node number of element i that corresponds to nodeNum
    if (1 > idx) {
      idx = 0;
    }

    if (idx != 0) {
      // assigns to an elementList and localNode list
      idx = static_cast<int>(b_index) - 1;
      elList->data[idx] = static_cast<double>(exponent) + 1.0;
      localNode->data[idx] = ii_data[0];
      b_index++;
    }
  }
}

//
// File trailer for findElementsAssociatedWithNodeNumber.cpp
//
// [EOF]
//
