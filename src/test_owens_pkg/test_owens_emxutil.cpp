//
// File: test_owens_emxutil.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 16-Apr-2020 09:21:06
//

// Include Files
#include "test_owens_emxutil.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include <cstdlib>
#include <cstring>
#include <string.h>

// Function Declarations
static void emxCopy_real_T(emxArray_real_T **dst, emxArray_real_T * const *src);
static void emxExpand_struct_T(d_emxArray_struct_T *emxArray, int fromIndex, int
  toIndex);
static void emxExpand_struct_T1(f_emxArray_struct_T *emxArray, int fromIndex,
  int toIndex);
static void emxExpand_struct_T2(e_emxArray_struct_T *emxArray, int fromIndex,
  int toIndex);
static void emxFreeStruct_struct_T8(s_struct_T *pStruct);
static void emxInitStruct_struct_T8(r_struct_T *pStruct);
static void emxInitStruct_struct_T9(s_struct_T *pStruct);
static void emxTrim_struct_T(d_emxArray_struct_T *emxArray, int fromIndex, int
  toIndex);
static void emxTrim_struct_T1(f_emxArray_struct_T *emxArray, int fromIndex, int
  toIndex);
static void emxTrim_struct_T2(e_emxArray_struct_T *emxArray, int fromIndex, int
  toIndex);

// Function Definitions

//
// Arguments    : emxArray_real_T **dst
//                emxArray_real_T * const *src
// Return Type  : void
//
static void emxCopy_real_T(emxArray_real_T **dst, emxArray_real_T * const *src)
{
  int numElDst;
  int numElSrc;
  int i;
  numElDst = 1;
  numElSrc = 1;
  for (i = 0; i < (*dst)->numDimensions; i++) {
    numElDst *= (*dst)->size[i];
    numElSrc *= (*src)->size[i];
  }

  for (i = 0; i < (*dst)->numDimensions; i++) {
    (*dst)->size[i] = (*src)->size[i];
  }

  emxEnsureCapacity_real_T(*dst, numElDst);
  for (i = 0; i < numElSrc; i++) {
    (*dst)->data[i] = (*src)->data[i];
  }
}

//
// Arguments    : d_emxArray_struct_T *emxArray
//                int fromIndex
//                int toIndex
// Return Type  : void
//
static void emxExpand_struct_T(d_emxArray_struct_T *emxArray, int fromIndex, int
  toIndex)
{
  int i;
  for (i = fromIndex; i < toIndex; i++) {
    emxInitStruct_struct_T5(&emxArray->data[i]);
  }
}

//
// Arguments    : f_emxArray_struct_T *emxArray
//                int fromIndex
//                int toIndex
// Return Type  : void
//
static void emxExpand_struct_T1(f_emxArray_struct_T *emxArray, int fromIndex,
  int toIndex)
{
  int i;
  for (i = fromIndex; i < toIndex; i++) {
    emxInitStruct_struct_T6(&emxArray->data[i]);
  }
}

//
// Arguments    : e_emxArray_struct_T *emxArray
//                int fromIndex
//                int toIndex
// Return Type  : void
//
static void emxExpand_struct_T2(e_emxArray_struct_T *emxArray, int fromIndex,
  int toIndex)
{
  int i;
  for (i = fromIndex; i < toIndex; i++) {
    emxInitStruct_struct_T4(&emxArray->data[i]);
  }
}

//
// Arguments    : s_struct_T *pStruct
// Return Type  : void
//
static void emxFreeStruct_struct_T8(s_struct_T *pStruct)
{
  emxFree_real_T(&pStruct->pBC);
  emxFree_real_T(&pStruct->isConstrained);
  emxFree_real_T(&pStruct->map);
  emxFree_real_T(&pStruct->redVectorMap);
}

//
// Arguments    : r_struct_T *pStruct
// Return Type  : void
//
static void emxInitStruct_struct_T8(r_struct_T *pStruct)
{
  pStruct->bladeNum.size[0] = 0;
  pStruct->h.size[0] = 0;
  pStruct->nodeNum.size[0] = 0;
  pStruct->elementNum.size[0] = 0;
  pStruct->remaining.size[0] = 0;
  pStruct->remaining.size[1] = 0;
}

//
// Arguments    : s_struct_T *pStruct
// Return Type  : void
//
static void emxInitStruct_struct_T9(s_struct_T *pStruct)
{
  emxInit_real_T(&pStruct->pBC, 2);
  emxInit_real_T(&pStruct->isConstrained, 1);
  emxInit_real_T(&pStruct->map, 1);
  emxInit_real_T(&pStruct->redVectorMap, 1);
}

//
// Arguments    : d_emxArray_struct_T *emxArray
//                int fromIndex
//                int toIndex
// Return Type  : void
//
static void emxTrim_struct_T(d_emxArray_struct_T *emxArray, int fromIndex, int
  toIndex)
{
  int i;
  for (i = fromIndex; i < toIndex; i++) {
    emxFreeStruct_struct_T4(&emxArray->data[i]);
  }
}

//
// Arguments    : f_emxArray_struct_T *emxArray
//                int fromIndex
//                int toIndex
// Return Type  : void
//
static void emxTrim_struct_T1(f_emxArray_struct_T *emxArray, int fromIndex, int
  toIndex)
{
  int i;
  for (i = fromIndex; i < toIndex; i++) {
    emxFreeStruct_struct_T6(&emxArray->data[i]);
  }
}

//
// Arguments    : e_emxArray_struct_T *emxArray
//                int fromIndex
//                int toIndex
// Return Type  : void
//
static void emxTrim_struct_T2(e_emxArray_struct_T *emxArray, int fromIndex, int
  toIndex)
{
  int i;
  for (i = fromIndex; i < toIndex; i++) {
    emxFreeStruct_struct_T5(&emxArray->data[i]);
  }
}

//
// Arguments    : coder_internal_sparse *pStruct
// Return Type  : void
//
void c_emxFreeStruct_coder_internal_(coder_internal_sparse *pStruct)
{
  emxFree_real_T(&pStruct->d);
  emxFree_int32_T(&pStruct->colidx);
  emxFree_int32_T(&pStruct->rowidx);
}

//
// Arguments    : coder_internal_sparse *pStruct
// Return Type  : void
//
void c_emxInitStruct_coder_internal_(coder_internal_sparse *pStruct)
{
  emxInit_real_T(&pStruct->d, 1);
  emxInit_int32_T(&pStruct->colidx, 1);
  emxInit_int32_T(&pStruct->rowidx, 1);
}

//
// Arguments    : l_struct_T *dst
//                const l_struct_T *src
// Return Type  : void
//
void emxCopyStruct_struct_T(l_struct_T *dst, const l_struct_T *src)
{
  dst->NElem = src->NElem;
  dst->FlipN = src->FlipN;
  emxCopy_real_T(&dst->QCx, &src->QCx);
  emxCopy_real_T(&dst->QCy, &src->QCy);
  emxCopy_real_T(&dst->QCz, &src->QCz);
  emxCopy_real_T(&dst->tx, &src->tx);
  emxCopy_real_T(&dst->ty, &src->ty);
  emxCopy_real_T(&dst->tz, &src->tz);
  emxCopy_real_T(&dst->CtoR, &src->CtoR);
  emxCopy_real_T(&dst->PEx, &src->PEx);
  emxCopy_real_T(&dst->PEy, &src->PEy);
  emxCopy_real_T(&dst->PEz, &src->PEz);
  emxCopy_real_T(&dst->tEx, &src->tEx);
  emxCopy_real_T(&dst->tEy, &src->tEy);
  emxCopy_real_T(&dst->tEz, &src->tEz);
  emxCopy_real_T(&dst->nEx, &src->nEx);
  emxCopy_real_T(&dst->nEy, &src->nEy);
  emxCopy_real_T(&dst->nEz, &src->nEz);
  emxCopy_real_T(&dst->sEx, &src->sEx);
  emxCopy_real_T(&dst->sEy, &src->sEy);
  emxCopy_real_T(&dst->sEz, &src->sEz);
  emxCopy_real_T(&dst->ECtoR, &src->ECtoR);
  emxCopy_real_T(&dst->EAreaR, &src->EAreaR);
  emxCopy_real_T(&dst->iSect, &src->iSect);
}

//
// Arguments    : n_struct_T *dst
//                const n_struct_T *src
// Return Type  : void
//
void emxCopyStruct_struct_T1(n_struct_T *dst, const n_struct_T *src)
{
  dst->NElem = src->NElem;
  dst->TtoC = src->TtoC;
  emxCopy_real_T(&dst->MCx, &src->MCx);
  emxCopy_real_T(&dst->MCy, &src->MCy);
  emxCopy_real_T(&dst->MCz, &src->MCz);
  emxCopy_real_T(&dst->CtoR, &src->CtoR);
  emxCopy_real_T(&dst->PEx, &src->PEx);
  emxCopy_real_T(&dst->PEy, &src->PEy);
  emxCopy_real_T(&dst->PEz, &src->PEz);
  emxCopy_real_T(&dst->sEx, &src->sEx);
  emxCopy_real_T(&dst->sEy, &src->sEy);
  emxCopy_real_T(&dst->sEz, &src->sEz);
  emxCopy_real_T(&dst->ECtoR, &src->ECtoR);
  emxCopy_real_T(&dst->EAreaR, &src->EAreaR);
  dst->BIndS = src->BIndS;
  dst->EIndS = src->EIndS;
  dst->BIndE = src->BIndE;
  dst->EIndE = src->EIndE;
}

//
// Arguments    : m_struct_T *dst
//                const m_struct_T *src
// Return Type  : void
//
void emxCopyStruct_struct_T2(m_struct_T *dst, const m_struct_T *src)
{
  emxCopy_real_T(&dst->N, &src->N);
  emxCopy_real_T(&dst->T, &src->T);
  emxCopy_real_T(&dst->M25, &src->M25);
}

//
// Arguments    : emxArray_boolean_T *emxArray
//                int oldNumel
// Return Type  : void
//
void emxEnsureCapacity_boolean_T(emxArray_boolean_T *emxArray, int oldNumel)
{
  int newNumel;
  int i;
  void *newData;
  if (oldNumel < 0) {
    oldNumel = 0;
  }

  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }

  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }

    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i *= 2;
      }
    }

    newData = std::calloc(static_cast<unsigned int>(i), sizeof(boolean_T));
    if (emxArray->data != NULL) {
      std::memcpy(newData, emxArray->data, sizeof(boolean_T) * oldNumel);
      if (emxArray->canFreeData) {
        std::free(emxArray->data);
      }
    }

    emxArray->data = (boolean_T *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
}

//
// Arguments    : emxArray_char_T *emxArray
//                int oldNumel
// Return Type  : void
//
void emxEnsureCapacity_char_T(emxArray_char_T *emxArray, int oldNumel)
{
  int newNumel;
  int i;
  void *newData;
  if (oldNumel < 0) {
    oldNumel = 0;
  }

  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }

  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }

    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i *= 2;
      }
    }

    newData = std::calloc(static_cast<unsigned int>(i), sizeof(char));
    if (emxArray->data != NULL) {
      std::memcpy(newData, emxArray->data, sizeof(char) * oldNumel);
      if (emxArray->canFreeData) {
        std::free(emxArray->data);
      }
    }

    emxArray->data = (char *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
}

//
// Arguments    : emxArray_creal_T *emxArray
//                int oldNumel
// Return Type  : void
//
void emxEnsureCapacity_creal_T(emxArray_creal_T *emxArray, int oldNumel)
{
  int newNumel;
  int i;
  void *newData;
  if (oldNumel < 0) {
    oldNumel = 0;
  }

  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }

  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }

    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i *= 2;
      }
    }

    newData = std::calloc(static_cast<unsigned int>(i), sizeof(creal_T));
    if (emxArray->data != NULL) {
      std::memcpy(newData, emxArray->data, sizeof(creal_T) * oldNumel);
      if (emxArray->canFreeData) {
        std::free(emxArray->data);
      }
    }

    emxArray->data = (creal_T *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
}

//
// Arguments    : emxArray_int32_T *emxArray
//                int oldNumel
// Return Type  : void
//
void emxEnsureCapacity_int32_T(emxArray_int32_T *emxArray, int oldNumel)
{
  int newNumel;
  int i;
  void *newData;
  if (oldNumel < 0) {
    oldNumel = 0;
  }

  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }

  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }

    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i *= 2;
      }
    }

    newData = std::calloc(static_cast<unsigned int>(i), sizeof(int));
    if (emxArray->data != NULL) {
      std::memcpy(newData, emxArray->data, sizeof(int) * oldNumel);
      if (emxArray->canFreeData) {
        std::free(emxArray->data);
      }
    }

    emxArray->data = (int *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
}

//
// Arguments    : emxArray_int8_T *emxArray
//                int oldNumel
// Return Type  : void
//
void emxEnsureCapacity_int8_T(emxArray_int8_T *emxArray, int oldNumel)
{
  int newNumel;
  int i;
  void *newData;
  if (oldNumel < 0) {
    oldNumel = 0;
  }

  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }

  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }

    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i *= 2;
      }
    }

    newData = std::calloc(static_cast<unsigned int>(i), sizeof(signed char));
    if (emxArray->data != NULL) {
      std::memcpy(newData, emxArray->data, sizeof(signed char) * oldNumel);
      if (emxArray->canFreeData) {
        std::free(emxArray->data);
      }
    }

    emxArray->data = (signed char *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
}

//
// Arguments    : emxArray_real_T *emxArray
//                int oldNumel
// Return Type  : void
//
void emxEnsureCapacity_real_T(emxArray_real_T *emxArray, int oldNumel)
{
  int newNumel;
  int i;
  void *newData;
  if (oldNumel < 0) {
    oldNumel = 0;
  }

  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }

  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }

    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i *= 2;
      }
    }

    newData = std::calloc(static_cast<unsigned int>(i), sizeof(double));
    if (emxArray->data != NULL) {
      std::memcpy(newData, emxArray->data, sizeof(double) * oldNumel);
      if (emxArray->canFreeData) {
        std::free(emxArray->data);
      }
    }

    emxArray->data = (double *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
}

//
// Arguments    : emxArray_struct_T *emxArray
//                int oldNumel
// Return Type  : void
//
void emxEnsureCapacity_struct_T(emxArray_struct_T *emxArray, int oldNumel)
{
  int newNumel;
  int i;
  void *newData;
  if (oldNumel < 0) {
    oldNumel = 0;
  }

  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }

  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }

    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i *= 2;
      }
    }

    newData = std::calloc(static_cast<unsigned int>(i), sizeof(e_struct_T));
    if (emxArray->data != NULL) {
      std::memcpy(newData, emxArray->data, sizeof(e_struct_T) * oldNumel);
      if (emxArray->canFreeData) {
        std::free(emxArray->data);
      }
    }

    emxArray->data = (e_struct_T *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
}

//
// Arguments    : b_emxArray_struct_T *emxArray
//                int oldNumel
// Return Type  : void
//
void emxEnsureCapacity_struct_T1(b_emxArray_struct_T *emxArray, int oldNumel)
{
  int newNumel;
  int i;
  void *newData;
  if (oldNumel < 0) {
    oldNumel = 0;
  }

  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }

  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }

    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i *= 2;
      }
    }

    newData = std::calloc(static_cast<unsigned int>(i), sizeof(f_struct_T));
    if (emxArray->data != NULL) {
      std::memcpy(newData, emxArray->data, sizeof(f_struct_T) * oldNumel);
      if (emxArray->canFreeData) {
        std::free(emxArray->data);
      }
    }

    emxArray->data = (f_struct_T *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
}

//
// Arguments    : d_emxArray_struct_T *emxArray
//                int oldNumel
// Return Type  : void
//
void emxEnsureCapacity_struct_T2(d_emxArray_struct_T *emxArray, int oldNumel)
{
  int newNumel;
  int i;
  void *newData;
  if (oldNumel < 0) {
    oldNumel = 0;
  }

  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }

  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }

    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i *= 2;
      }
    }

    newData = std::calloc(static_cast<unsigned int>(i), sizeof(l_struct_T));
    if (emxArray->data != NULL) {
      std::memcpy(newData, emxArray->data, sizeof(l_struct_T) * oldNumel);
      if (emxArray->canFreeData) {
        std::free(emxArray->data);
      }
    }

    emxArray->data = (l_struct_T *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }

  if (oldNumel > newNumel) {
    emxTrim_struct_T(emxArray, newNumel, oldNumel);
  } else {
    if (oldNumel < newNumel) {
      emxExpand_struct_T(emxArray, oldNumel, newNumel);
    }
  }
}

//
// Arguments    : f_emxArray_struct_T *emxArray
//                int oldNumel
// Return Type  : void
//
void emxEnsureCapacity_struct_T3(f_emxArray_struct_T *emxArray, int oldNumel)
{
  int newNumel;
  int i;
  void *newData;
  if (oldNumel < 0) {
    oldNumel = 0;
  }

  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }

  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }

    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i *= 2;
      }
    }

    newData = std::calloc(static_cast<unsigned int>(i), sizeof(n_struct_T));
    if (emxArray->data != NULL) {
      std::memcpy(newData, emxArray->data, sizeof(n_struct_T) * oldNumel);
      if (emxArray->canFreeData) {
        std::free(emxArray->data);
      }
    }

    emxArray->data = (n_struct_T *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }

  if (oldNumel > newNumel) {
    emxTrim_struct_T1(emxArray, newNumel, oldNumel);
  } else {
    if (oldNumel < newNumel) {
      emxExpand_struct_T1(emxArray, oldNumel, newNumel);
    }
  }
}

//
// Arguments    : e_emxArray_struct_T *emxArray
//                int oldNumel
// Return Type  : void
//
void emxEnsureCapacity_struct_T4(e_emxArray_struct_T *emxArray, int oldNumel)
{
  int newNumel;
  int i;
  void *newData;
  if (oldNumel < 0) {
    oldNumel = 0;
  }

  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }

  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }

    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i *= 2;
      }
    }

    newData = std::calloc(static_cast<unsigned int>(i), sizeof(m_struct_T));
    if (emxArray->data != NULL) {
      std::memcpy(newData, emxArray->data, sizeof(m_struct_T) * oldNumel);
      if (emxArray->canFreeData) {
        std::free(emxArray->data);
      }
    }

    emxArray->data = (m_struct_T *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }

  if (oldNumel > newNumel) {
    emxTrim_struct_T2(emxArray, newNumel, oldNumel);
  } else {
    if (oldNumel < newNumel) {
      emxExpand_struct_T2(emxArray, oldNumel, newNumel);
    }
  }
}

//
// Arguments    : c_emxArray_struct_T *emxArray
//                int oldNumel
// Return Type  : void
//
void emxEnsureCapacity_struct_T5(c_emxArray_struct_T *emxArray, int oldNumel)
{
  int newNumel;
  int i;
  void *newData;
  if (oldNumel < 0) {
    oldNumel = 0;
  }

  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }

  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }

    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i *= 2;
      }
    }

    newData = std::calloc(static_cast<unsigned int>(i), sizeof(g_struct_T));
    if (emxArray->data != NULL) {
      std::memcpy(newData, emxArray->data, sizeof(g_struct_T) * oldNumel);
      if (emxArray->canFreeData) {
        std::free(emxArray->data);
      }
    }

    emxArray->data = (g_struct_T *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
}

//
// Arguments    : emxArray_uint32_T *emxArray
//                int oldNumel
// Return Type  : void
//
void emxEnsureCapacity_uint32_T(emxArray_uint32_T *emxArray, int oldNumel)
{
  int newNumel;
  int i;
  void *newData;
  if (oldNumel < 0) {
    oldNumel = 0;
  }

  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }

  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }

    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i *= 2;
      }
    }

    newData = std::calloc(static_cast<unsigned int>(i), sizeof(unsigned int));
    if (emxArray->data != NULL) {
      std::memcpy(newData, emxArray->data, sizeof(unsigned int) * oldNumel);
      if (emxArray->canFreeData) {
        std::free(emxArray->data);
      }
    }

    emxArray->data = (unsigned int *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
}

//
// Arguments    : h_struct_T *pStruct
// Return Type  : void
//
void emxFreeStruct_struct_T(h_struct_T *pStruct)
{
  emxFree_real_T(&pStruct->nodeNum);
  emxFree_real_T(&pStruct->x);
  emxFree_real_T(&pStruct->y);
  emxFree_real_T(&pStruct->z);
  emxFree_real_T(&pStruct->elNum);
  emxFree_real_T(&pStruct->conn);
}

//
// Arguments    : i_struct_T *pStruct
// Return Type  : void
//
void emxFreeStruct_struct_T1(i_struct_T *pStruct)
{
  emxFree_struct_T(&pStruct->props);
  emxFree_real_T(&pStruct->elLen);
  emxFree_real_T(&pStruct->psi);
  emxFree_real_T(&pStruct->theta);
  emxFree_real_T(&pStruct->roll);
  emxFree_boolean_T(&pStruct->rotationalEffects);
}

//
// Arguments    : j_struct_T *pStruct
// Return Type  : void
//
void emxFreeStruct_struct_T2(j_struct_T *pStruct)
{
  emxFree_real_T(&pStruct->displ_s);
  emxFree_real_T(&pStruct->displdot_s);
  emxFree_real_T(&pStruct->displddot_s);
}

//
// Arguments    : k_struct_T *pStruct
// Return Type  : void
//
void emxFreeStruct_struct_T3(k_struct_T *pStruct)
{
  emxFree_struct_T1(&pStruct->elStrain);
  emxFree_real_T(&pStruct->displ_sp1);
  emxFree_real_T(&pStruct->displddot_sp1);
  emxFree_real_T(&pStruct->displdot_sp1);
}

//
// Arguments    : l_struct_T *pStruct
// Return Type  : void
//
void emxFreeStruct_struct_T4(l_struct_T *pStruct)
{
  emxFree_real_T(&pStruct->QCx);
  emxFree_real_T(&pStruct->QCy);
  emxFree_real_T(&pStruct->QCz);
  emxFree_real_T(&pStruct->tx);
  emxFree_real_T(&pStruct->ty);
  emxFree_real_T(&pStruct->tz);
  emxFree_real_T(&pStruct->CtoR);
  emxFree_real_T(&pStruct->PEx);
  emxFree_real_T(&pStruct->PEy);
  emxFree_real_T(&pStruct->PEz);
  emxFree_real_T(&pStruct->tEx);
  emxFree_real_T(&pStruct->tEy);
  emxFree_real_T(&pStruct->tEz);
  emxFree_real_T(&pStruct->nEx);
  emxFree_real_T(&pStruct->nEy);
  emxFree_real_T(&pStruct->nEz);
  emxFree_real_T(&pStruct->sEx);
  emxFree_real_T(&pStruct->sEy);
  emxFree_real_T(&pStruct->sEz);
  emxFree_real_T(&pStruct->ECtoR);
  emxFree_real_T(&pStruct->EAreaR);
  emxFree_real_T(&pStruct->iSect);
}

//
// Arguments    : m_struct_T *pStruct
// Return Type  : void
//
void emxFreeStruct_struct_T5(m_struct_T *pStruct)
{
  emxFree_real_T(&pStruct->N);
  emxFree_real_T(&pStruct->T);
  emxFree_real_T(&pStruct->M25);
}

//
// Arguments    : n_struct_T *pStruct
// Return Type  : void
//
void emxFreeStruct_struct_T6(n_struct_T *pStruct)
{
  emxFree_real_T(&pStruct->MCx);
  emxFree_real_T(&pStruct->MCy);
  emxFree_real_T(&pStruct->MCz);
  emxFree_real_T(&pStruct->CtoR);
  emxFree_real_T(&pStruct->PEx);
  emxFree_real_T(&pStruct->PEy);
  emxFree_real_T(&pStruct->PEz);
  emxFree_real_T(&pStruct->sEx);
  emxFree_real_T(&pStruct->sEy);
  emxFree_real_T(&pStruct->sEz);
  emxFree_real_T(&pStruct->ECtoR);
  emxFree_real_T(&pStruct->EAreaR);
}

//
// Arguments    : t_struct_T *pStruct
// Return Type  : void
//
void emxFreeStruct_struct_T7(t_struct_T *pStruct)
{
  emxFree_char_T(&pStruct->aeroloadfile);
  emxFree_char_T(&pStruct->owensfile);
  emxFreeStruct_struct_T8(&pStruct->BC);
  emxFree_real_T(&pStruct->joint);
  emxFree_real_T(&pStruct->jointTransform);
  emxFree_real_T(&pStruct->reducedDOFList);
}

//
// Arguments    : emxArray_boolean_T **pEmxArray
// Return Type  : void
//
void emxFree_boolean_T(emxArray_boolean_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_boolean_T *)NULL) {
    if (((*pEmxArray)->data != (boolean_T *)NULL) && (*pEmxArray)->canFreeData)
    {
      std::free((*pEmxArray)->data);
    }

    std::free((*pEmxArray)->size);
    std::free(*pEmxArray);
    *pEmxArray = (emxArray_boolean_T *)NULL;
  }
}

//
// Arguments    : emxArray_char_T **pEmxArray
// Return Type  : void
//
void emxFree_char_T(emxArray_char_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_char_T *)NULL) {
    if (((*pEmxArray)->data != (char *)NULL) && (*pEmxArray)->canFreeData) {
      std::free((*pEmxArray)->data);
    }

    std::free((*pEmxArray)->size);
    std::free(*pEmxArray);
    *pEmxArray = (emxArray_char_T *)NULL;
  }
}

//
// Arguments    : emxArray_creal_T **pEmxArray
// Return Type  : void
//
void emxFree_creal_T(emxArray_creal_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_creal_T *)NULL) {
    if (((*pEmxArray)->data != (creal_T *)NULL) && (*pEmxArray)->canFreeData) {
      std::free((*pEmxArray)->data);
    }

    std::free((*pEmxArray)->size);
    std::free(*pEmxArray);
    *pEmxArray = (emxArray_creal_T *)NULL;
  }
}

//
// Arguments    : emxArray_int32_T **pEmxArray
// Return Type  : void
//
void emxFree_int32_T(emxArray_int32_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_int32_T *)NULL) {
    if (((*pEmxArray)->data != (int *)NULL) && (*pEmxArray)->canFreeData) {
      std::free((*pEmxArray)->data);
    }

    std::free((*pEmxArray)->size);
    std::free(*pEmxArray);
    *pEmxArray = (emxArray_int32_T *)NULL;
  }
}

//
// Arguments    : emxArray_int8_T **pEmxArray
// Return Type  : void
//
void emxFree_int8_T(emxArray_int8_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_int8_T *)NULL) {
    if (((*pEmxArray)->data != (signed char *)NULL) && (*pEmxArray)->canFreeData)
    {
      std::free((*pEmxArray)->data);
    }

    std::free((*pEmxArray)->size);
    std::free(*pEmxArray);
    *pEmxArray = (emxArray_int8_T *)NULL;
  }
}

//
// Arguments    : emxArray_real_T **pEmxArray
// Return Type  : void
//
void emxFree_real_T(emxArray_real_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_real_T *)NULL) {
    if (((*pEmxArray)->data != (double *)NULL) && (*pEmxArray)->canFreeData) {
      std::free((*pEmxArray)->data);
    }

    std::free((*pEmxArray)->size);
    std::free(*pEmxArray);
    *pEmxArray = (emxArray_real_T *)NULL;
  }
}

//
// Arguments    : emxArray_struct_T **pEmxArray
// Return Type  : void
//
void emxFree_struct_T(emxArray_struct_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_struct_T *)NULL) {
    if (((*pEmxArray)->data != (e_struct_T *)NULL) && (*pEmxArray)->canFreeData)
    {
      std::free((*pEmxArray)->data);
    }

    std::free((*pEmxArray)->size);
    std::free(*pEmxArray);
    *pEmxArray = (emxArray_struct_T *)NULL;
  }
}

//
// Arguments    : b_emxArray_struct_T **pEmxArray
// Return Type  : void
//
void emxFree_struct_T1(b_emxArray_struct_T **pEmxArray)
{
  if (*pEmxArray != (b_emxArray_struct_T *)NULL) {
    if (((*pEmxArray)->data != (f_struct_T *)NULL) && (*pEmxArray)->canFreeData)
    {
      std::free((*pEmxArray)->data);
    }

    std::free((*pEmxArray)->size);
    std::free(*pEmxArray);
    *pEmxArray = (b_emxArray_struct_T *)NULL;
  }
}

//
// Arguments    : c_emxArray_struct_T **pEmxArray
// Return Type  : void
//
void emxFree_struct_T2(c_emxArray_struct_T **pEmxArray)
{
  if (*pEmxArray != (c_emxArray_struct_T *)NULL) {
    if (((*pEmxArray)->data != (g_struct_T *)NULL) && (*pEmxArray)->canFreeData)
    {
      std::free((*pEmxArray)->data);
    }

    std::free((*pEmxArray)->size);
    std::free(*pEmxArray);
    *pEmxArray = (c_emxArray_struct_T *)NULL;
  }
}

//
// Arguments    : d_emxArray_struct_T **pEmxArray
// Return Type  : void
//
void emxFree_struct_T3(d_emxArray_struct_T **pEmxArray)
{
  int numEl;
  int i;
  if (*pEmxArray != (d_emxArray_struct_T *)NULL) {
    if ((*pEmxArray)->data != (l_struct_T *)NULL) {
      numEl = 1;
      for (i = 0; i < (*pEmxArray)->numDimensions; i++) {
        numEl *= (*pEmxArray)->size[i];
      }

      for (i = 0; i < numEl; i++) {
        emxFreeStruct_struct_T4(&(*pEmxArray)->data[i]);
      }

      if ((*pEmxArray)->canFreeData) {
        std::free((*pEmxArray)->data);
      }
    }

    std::free((*pEmxArray)->size);
    std::free(*pEmxArray);
    *pEmxArray = (d_emxArray_struct_T *)NULL;
  }
}

//
// Arguments    : e_emxArray_struct_T **pEmxArray
// Return Type  : void
//
void emxFree_struct_T4(e_emxArray_struct_T **pEmxArray)
{
  int numEl;
  int i;
  if (*pEmxArray != (e_emxArray_struct_T *)NULL) {
    if ((*pEmxArray)->data != (m_struct_T *)NULL) {
      numEl = 1;
      for (i = 0; i < (*pEmxArray)->numDimensions; i++) {
        numEl *= (*pEmxArray)->size[i];
      }

      for (i = 0; i < numEl; i++) {
        emxFreeStruct_struct_T5(&(*pEmxArray)->data[i]);
      }

      if ((*pEmxArray)->canFreeData) {
        std::free((*pEmxArray)->data);
      }
    }

    std::free((*pEmxArray)->size);
    std::free(*pEmxArray);
    *pEmxArray = (e_emxArray_struct_T *)NULL;
  }
}

//
// Arguments    : f_emxArray_struct_T **pEmxArray
// Return Type  : void
//
void emxFree_struct_T5(f_emxArray_struct_T **pEmxArray)
{
  int numEl;
  int i;
  if (*pEmxArray != (f_emxArray_struct_T *)NULL) {
    if ((*pEmxArray)->data != (n_struct_T *)NULL) {
      numEl = 1;
      for (i = 0; i < (*pEmxArray)->numDimensions; i++) {
        numEl *= (*pEmxArray)->size[i];
      }

      for (i = 0; i < numEl; i++) {
        emxFreeStruct_struct_T6(&(*pEmxArray)->data[i]);
      }

      if ((*pEmxArray)->canFreeData) {
        std::free((*pEmxArray)->data);
      }
    }

    std::free((*pEmxArray)->size);
    std::free(*pEmxArray);
    *pEmxArray = (f_emxArray_struct_T *)NULL;
  }
}

//
// Arguments    : emxArray_uint32_T **pEmxArray
// Return Type  : void
//
void emxFree_uint32_T(emxArray_uint32_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_uint32_T *)NULL) {
    if (((*pEmxArray)->data != (unsigned int *)NULL) && (*pEmxArray)
        ->canFreeData) {
      std::free((*pEmxArray)->data);
    }

    std::free((*pEmxArray)->size);
    std::free(*pEmxArray);
    *pEmxArray = (emxArray_uint32_T *)NULL;
  }
}

//
// Arguments    : h_struct_T *pStruct
// Return Type  : void
//
void emxInitStruct_struct_T(h_struct_T *pStruct)
{
  emxInit_real_T(&pStruct->nodeNum, 1);
  emxInit_real_T(&pStruct->x, 1);
  emxInit_real_T(&pStruct->y, 1);
  emxInit_real_T(&pStruct->z, 1);
  emxInit_real_T(&pStruct->elNum, 1);
  emxInit_real_T(&pStruct->conn, 2);
}

//
// Arguments    : i_struct_T *pStruct
// Return Type  : void
//
void emxInitStruct_struct_T1(i_struct_T *pStruct)
{
  emxInit_struct_T(&pStruct->props, 2);
  emxInit_real_T(&pStruct->elLen, 1);
  emxInit_real_T(&pStruct->psi, 1);
  emxInit_real_T(&pStruct->theta, 1);
  emxInit_real_T(&pStruct->roll, 1);
  emxInit_boolean_T(&pStruct->rotationalEffects, 2);
}

//
// Arguments    : j_struct_T *pStruct
// Return Type  : void
//
void emxInitStruct_struct_T2(j_struct_T *pStruct)
{
  emxInit_real_T(&pStruct->displ_s, 1);
  emxInit_real_T(&pStruct->displdot_s, 1);
  emxInit_real_T(&pStruct->displddot_s, 1);
}

//
// Arguments    : k_struct_T *pStruct
// Return Type  : void
//
void emxInitStruct_struct_T3(k_struct_T *pStruct)
{
  emxInit_struct_T1(&pStruct->elStrain, 2);
  emxInit_real_T(&pStruct->displ_sp1, 1);
  emxInit_real_T(&pStruct->displddot_sp1, 1);
  emxInit_real_T(&pStruct->displdot_sp1, 1);
}

//
// Arguments    : m_struct_T *pStruct
// Return Type  : void
//
void emxInitStruct_struct_T4(m_struct_T *pStruct)
{
  emxInit_real_T(&pStruct->N, 2);
  emxInit_real_T(&pStruct->T, 2);
  emxInit_real_T(&pStruct->M25, 2);
}

//
// Arguments    : l_struct_T *pStruct
// Return Type  : void
//
void emxInitStruct_struct_T5(l_struct_T *pStruct)
{
  emxInit_real_T(&pStruct->QCx, 1);
  emxInit_real_T(&pStruct->QCy, 1);
  emxInit_real_T(&pStruct->QCz, 1);
  emxInit_real_T(&pStruct->tx, 1);
  emxInit_real_T(&pStruct->ty, 1);
  emxInit_real_T(&pStruct->tz, 1);
  emxInit_real_T(&pStruct->CtoR, 1);
  emxInit_real_T(&pStruct->PEx, 1);
  emxInit_real_T(&pStruct->PEy, 1);
  emxInit_real_T(&pStruct->PEz, 1);
  emxInit_real_T(&pStruct->tEx, 1);
  emxInit_real_T(&pStruct->tEy, 1);
  emxInit_real_T(&pStruct->tEz, 1);
  emxInit_real_T(&pStruct->nEx, 1);
  emxInit_real_T(&pStruct->nEy, 1);
  emxInit_real_T(&pStruct->nEz, 1);
  emxInit_real_T(&pStruct->sEx, 1);
  emxInit_real_T(&pStruct->sEy, 1);
  emxInit_real_T(&pStruct->sEz, 1);
  emxInit_real_T(&pStruct->ECtoR, 1);
  emxInit_real_T(&pStruct->EAreaR, 1);
  emxInit_real_T(&pStruct->iSect, 1);
}

//
// Arguments    : n_struct_T *pStruct
// Return Type  : void
//
void emxInitStruct_struct_T6(n_struct_T *pStruct)
{
  emxInit_real_T(&pStruct->MCx, 1);
  emxInit_real_T(&pStruct->MCy, 1);
  emxInit_real_T(&pStruct->MCz, 1);
  emxInit_real_T(&pStruct->CtoR, 1);
  emxInit_real_T(&pStruct->PEx, 1);
  emxInit_real_T(&pStruct->PEy, 1);
  emxInit_real_T(&pStruct->PEz, 1);
  emxInit_real_T(&pStruct->sEx, 1);
  emxInit_real_T(&pStruct->sEy, 1);
  emxInit_real_T(&pStruct->sEz, 1);
  emxInit_real_T(&pStruct->ECtoR, 1);
  emxInit_real_T(&pStruct->EAreaR, 1);
}

//
// Arguments    : t_struct_T *pStruct
// Return Type  : void
//
void emxInitStruct_struct_T7(t_struct_T *pStruct)
{
  emxInit_char_T(&pStruct->aeroloadfile, 2);
  emxInit_char_T(&pStruct->owensfile, 2);
  emxInitStruct_struct_T8(&pStruct->bladeData);
  emxInitStruct_struct_T9(&pStruct->BC);
  emxInit_real_T(&pStruct->joint, 2);
  pStruct->outFilename.size[0] = 0;
  pStruct->outFilename.size[1] = 0;
  emxInit_real_T(&pStruct->jointTransform, 2);
  emxInit_real_T(&pStruct->reducedDOFList, 2);
}

//
// Arguments    : emxArray_boolean_T **pEmxArray
//                int numDimensions
// Return Type  : void
//
void emxInit_boolean_T(emxArray_boolean_T **pEmxArray, int numDimensions)
{
  emxArray_boolean_T *emxArray;
  int i;
  *pEmxArray = (emxArray_boolean_T *)std::malloc(sizeof(emxArray_boolean_T));
  emxArray = *pEmxArray;
  emxArray->data = (boolean_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)std::malloc(sizeof(int) * numDimensions);
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

//
// Arguments    : emxArray_char_T **pEmxArray
//                int numDimensions
// Return Type  : void
//
void emxInit_char_T(emxArray_char_T **pEmxArray, int numDimensions)
{
  emxArray_char_T *emxArray;
  int i;
  *pEmxArray = (emxArray_char_T *)std::malloc(sizeof(emxArray_char_T));
  emxArray = *pEmxArray;
  emxArray->data = (char *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)std::malloc(sizeof(int) * numDimensions);
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

//
// Arguments    : emxArray_creal_T **pEmxArray
//                int numDimensions
// Return Type  : void
//
void emxInit_creal_T(emxArray_creal_T **pEmxArray, int numDimensions)
{
  emxArray_creal_T *emxArray;
  int i;
  *pEmxArray = (emxArray_creal_T *)std::malloc(sizeof(emxArray_creal_T));
  emxArray = *pEmxArray;
  emxArray->data = (creal_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)std::malloc(sizeof(int) * numDimensions);
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

//
// Arguments    : emxArray_int32_T **pEmxArray
//                int numDimensions
// Return Type  : void
//
void emxInit_int32_T(emxArray_int32_T **pEmxArray, int numDimensions)
{
  emxArray_int32_T *emxArray;
  int i;
  *pEmxArray = (emxArray_int32_T *)std::malloc(sizeof(emxArray_int32_T));
  emxArray = *pEmxArray;
  emxArray->data = (int *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)std::malloc(sizeof(int) * numDimensions);
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

//
// Arguments    : emxArray_int8_T **pEmxArray
//                int numDimensions
// Return Type  : void
//
void emxInit_int8_T(emxArray_int8_T **pEmxArray, int numDimensions)
{
  emxArray_int8_T *emxArray;
  int i;
  *pEmxArray = (emxArray_int8_T *)std::malloc(sizeof(emxArray_int8_T));
  emxArray = *pEmxArray;
  emxArray->data = (signed char *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)std::malloc(sizeof(int) * numDimensions);
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

//
// Arguments    : emxArray_real_T **pEmxArray
//                int numDimensions
// Return Type  : void
//
void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions)
{
  emxArray_real_T *emxArray;
  int i;
  *pEmxArray = (emxArray_real_T *)std::malloc(sizeof(emxArray_real_T));
  emxArray = *pEmxArray;
  emxArray->data = (double *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)std::malloc(sizeof(int) * numDimensions);
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

//
// Arguments    : emxArray_struct_T **pEmxArray
//                int numDimensions
// Return Type  : void
//
void emxInit_struct_T(emxArray_struct_T **pEmxArray, int numDimensions)
{
  emxArray_struct_T *emxArray;
  int i;
  *pEmxArray = (emxArray_struct_T *)std::malloc(sizeof(emxArray_struct_T));
  emxArray = *pEmxArray;
  emxArray->data = (e_struct_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)std::malloc(sizeof(int) * numDimensions);
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

//
// Arguments    : b_emxArray_struct_T **pEmxArray
//                int numDimensions
// Return Type  : void
//
void emxInit_struct_T1(b_emxArray_struct_T **pEmxArray, int numDimensions)
{
  b_emxArray_struct_T *emxArray;
  int i;
  *pEmxArray = (b_emxArray_struct_T *)std::malloc(sizeof(b_emxArray_struct_T));
  emxArray = *pEmxArray;
  emxArray->data = (f_struct_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)std::malloc(sizeof(int) * numDimensions);
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

//
// Arguments    : c_emxArray_struct_T **pEmxArray
//                int numDimensions
// Return Type  : void
//
void emxInit_struct_T2(c_emxArray_struct_T **pEmxArray, int numDimensions)
{
  c_emxArray_struct_T *emxArray;
  int i;
  *pEmxArray = (c_emxArray_struct_T *)std::malloc(sizeof(c_emxArray_struct_T));
  emxArray = *pEmxArray;
  emxArray->data = (g_struct_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)std::malloc(sizeof(int) * numDimensions);
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

//
// Arguments    : d_emxArray_struct_T **pEmxArray
//                int numDimensions
// Return Type  : void
//
void emxInit_struct_T3(d_emxArray_struct_T **pEmxArray, int numDimensions)
{
  d_emxArray_struct_T *emxArray;
  int i;
  *pEmxArray = (d_emxArray_struct_T *)std::malloc(sizeof(d_emxArray_struct_T));
  emxArray = *pEmxArray;
  emxArray->data = (l_struct_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)std::malloc(sizeof(int) * numDimensions);
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

//
// Arguments    : e_emxArray_struct_T **pEmxArray
//                int numDimensions
// Return Type  : void
//
void emxInit_struct_T4(e_emxArray_struct_T **pEmxArray, int numDimensions)
{
  e_emxArray_struct_T *emxArray;
  int i;
  *pEmxArray = (e_emxArray_struct_T *)std::malloc(sizeof(e_emxArray_struct_T));
  emxArray = *pEmxArray;
  emxArray->data = (m_struct_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)std::malloc(sizeof(int) * numDimensions);
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

//
// Arguments    : f_emxArray_struct_T **pEmxArray
//                int numDimensions
// Return Type  : void
//
void emxInit_struct_T5(f_emxArray_struct_T **pEmxArray, int numDimensions)
{
  f_emxArray_struct_T *emxArray;
  int i;
  *pEmxArray = (f_emxArray_struct_T *)std::malloc(sizeof(f_emxArray_struct_T));
  emxArray = *pEmxArray;
  emxArray->data = (n_struct_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)std::malloc(sizeof(int) * numDimensions);
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

//
// Arguments    : emxArray_uint32_T **pEmxArray
//                int numDimensions
// Return Type  : void
//
void emxInit_uint32_T(emxArray_uint32_T **pEmxArray, int numDimensions)
{
  emxArray_uint32_T *emxArray;
  int i;
  *pEmxArray = (emxArray_uint32_T *)std::malloc(sizeof(emxArray_uint32_T));
  emxArray = *pEmxArray;
  emxArray->data = (unsigned int *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)std::malloc(sizeof(int) * numDimensions);
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

//
// File trailer for test_owens_emxutil.cpp
//
// [EOF]
//
