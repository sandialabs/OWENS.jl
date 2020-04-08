//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: readMesh.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 08-Apr-2020 17:30:34
//

// Include Files
#include "readMesh.h"
#include "fileManager.h"
#include "getSplitLine.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include "test_owens_emxutil.h"
#include <string.h>

// Function Definitions

//
// readMesh  reads mesh file and stores data in mesh object
//  **********************************************************************
//  *                   Part of the SNL OWENS Toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [mesh] = readMesh(filename)
//
//    This function reads the mesh file and stores data in the mesh object.
//
//       input:
//       filename      = string containing mesh filename
// Arguments    : const emxArray_char_T *filename
//                emxArray_real_T *mesh_nodeNum
//                double *mesh_numEl
//                double *mesh_numNodes
//                emxArray_real_T *mesh_x
//                emxArray_real_T *mesh_y
//                emxArray_real_T *mesh_z
//                emxArray_real_T *mesh_elNum
//                emxArray_real_T *mesh_conn
// Return Type  : void
//
void b_readMesh(const emxArray_char_T *filename, emxArray_real_T *mesh_nodeNum,
                double *mesh_numEl, double *mesh_numNodes, emxArray_real_T
                *mesh_x, emxArray_real_T *mesh_y, emxArray_real_T *mesh_z,
                emxArray_real_T *mesh_elNum, emxArray_real_T *mesh_conn)
{
  emxArray_real_T *temp;
  emxArray_real_T *nodeNum;
  emxArray_real_T *conn;
  signed char fileid;
  int i;
  int i1;
  int i2;
  int loop_ub;
  emxArray_real_T *x;
  emxArray_real_T *y;
  emxArray_real_T *z;
  emxArray_real_T *elNum;
  emxArray_real_T *b_temp;
  int b_i;
  emxInit_real_T(&temp, 1);
  emxInit_real_T(&nodeNum, 1);
  emxInit_real_T(&conn, 2);

  //       output:
  //       mesh          = object containing mesh data
  fileid = c_cfopen(filename, "rb");

  // open mesh file
  //  temp = fscanf(fid,'%i',2);   %read in number of nodes and number of elements 
  getSplitLine(static_cast<double>(fileid), temp);
  i = static_cast<int>(temp->data[0]);
  i1 = nodeNum->size[0];
  nodeNum->size[0] = i;
  emxEnsureCapacity_real_T(nodeNum, i1);
  i1 = static_cast<int>(temp->data[1]);
  i2 = conn->size[0] * conn->size[1];
  conn->size[0] = i1;
  conn->size[1] = 2;
  emxEnsureCapacity_real_T(conn, i2);
  loop_ub = i1 << 1;
  for (i2 = 0; i2 < loop_ub; i2++) {
    conn->data[i2] = 0.0;
  }

  emxInit_real_T(&x, 1);
  emxInit_real_T(&y, 1);
  emxInit_real_T(&z, 1);
  emxInit_real_T(&elNum, 1);
  i2 = elNum->size[0];
  elNum->size[0] = i1;
  emxEnsureCapacity_real_T(elNum, i2);
  i2 = x->size[0];
  x->size[0] = i;
  emxEnsureCapacity_real_T(x, i2);
  i2 = y->size[0];
  y->size[0] = i;
  emxEnsureCapacity_real_T(y, i2);
  i2 = z->size[0];
  z->size[0] = i;
  emxEnsureCapacity_real_T(z, i2);
  emxInit_real_T(&b_temp, 1);
  for (b_i = 0; b_i < i; b_i++) {
    //  read in node number and node coordinates
    getSplitLine(static_cast<double>(fileid), b_temp);
    nodeNum->data[b_i] = b_temp->data[0];
    x->data[b_i] = b_temp->data[1];
    y->data[b_i] = b_temp->data[2];
    z->data[b_i] = b_temp->data[3];
  }

  for (b_i = 0; b_i < i1; b_i++) {
    //  read in element number and connectivity list
    getSplitLine(static_cast<double>(fileid), b_temp);
    elNum->data[b_i] = b_temp->data[0];
    conn->data[b_i] = b_temp->data[2];
    conn->data[b_i + conn->size[0]] = b_temp->data[3];
  }

  emxFree_real_T(&b_temp);
  cfclose(static_cast<double>(fileid));

  // close mesh file
  i = mesh_nodeNum->size[0];
  mesh_nodeNum->size[0] = nodeNum->size[0];
  emxEnsureCapacity_real_T(mesh_nodeNum, i);
  loop_ub = nodeNum->size[0];
  for (i = 0; i < loop_ub; i++) {
    mesh_nodeNum->data[i] = nodeNum->data[i];
  }

  emxFree_real_T(&nodeNum);

  // store data in mesh object
  *mesh_numEl = temp->data[1];
  *mesh_numNodes = temp->data[0];
  i = mesh_x->size[0];
  mesh_x->size[0] = x->size[0];
  emxEnsureCapacity_real_T(mesh_x, i);
  loop_ub = x->size[0];
  emxFree_real_T(&temp);
  for (i = 0; i < loop_ub; i++) {
    mesh_x->data[i] = x->data[i];
  }

  emxFree_real_T(&x);
  i = mesh_y->size[0];
  mesh_y->size[0] = y->size[0];
  emxEnsureCapacity_real_T(mesh_y, i);
  loop_ub = y->size[0];
  for (i = 0; i < loop_ub; i++) {
    mesh_y->data[i] = y->data[i];
  }

  emxFree_real_T(&y);
  i = mesh_z->size[0];
  mesh_z->size[0] = z->size[0];
  emxEnsureCapacity_real_T(mesh_z, i);
  loop_ub = z->size[0];
  for (i = 0; i < loop_ub; i++) {
    mesh_z->data[i] = z->data[i];
  }

  emxFree_real_T(&z);
  i = mesh_elNum->size[0];
  mesh_elNum->size[0] = elNum->size[0];
  emxEnsureCapacity_real_T(mesh_elNum, i);
  loop_ub = elNum->size[0];
  for (i = 0; i < loop_ub; i++) {
    mesh_elNum->data[i] = elNum->data[i];
  }

  emxFree_real_T(&elNum);
  i = mesh_conn->size[0] * mesh_conn->size[1];
  mesh_conn->size[0] = conn->size[0];
  mesh_conn->size[1] = 2;
  emxEnsureCapacity_real_T(mesh_conn, i);
  loop_ub = conn->size[0] * conn->size[1];
  for (i = 0; i < loop_ub; i++) {
    mesh_conn->data[i] = conn->data[i];
  }

  emxFree_real_T(&conn);
}

//
// readMesh  reads mesh file and stores data in mesh object
//  **********************************************************************
//  *                   Part of the SNL OWENS Toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [mesh] = readMesh(filename)
//
//    This function reads the mesh file and stores data in the mesh object.
//
//       input:
//       filename      = string containing mesh filename
// Arguments    : const emxArray_char_T *filename
//                emxArray_real_T *mesh_nodeNum
//                double *mesh_numEl
//                double *mesh_numNodes
//                emxArray_real_T *mesh_x
//                emxArray_real_T *mesh_y
//                emxArray_real_T *mesh_z
//                emxArray_real_T *mesh_elNum
//                emxArray_real_T *mesh_conn
// Return Type  : void
//
void readMesh(const emxArray_char_T *filename, emxArray_real_T *mesh_nodeNum,
              double *mesh_numEl, double *mesh_numNodes, emxArray_real_T *mesh_x,
              emxArray_real_T *mesh_y, emxArray_real_T *mesh_z, emxArray_real_T *
              mesh_elNum, emxArray_real_T *mesh_conn)
{
  emxArray_real_T *temp;
  emxArray_real_T *nodeNum;
  emxArray_real_T *conn;
  signed char fileid;
  int i;
  int i1;
  int i2;
  int loop_ub;
  emxArray_real_T *x;
  emxArray_real_T *y;
  emxArray_real_T *z;
  emxArray_real_T *elNum;
  emxArray_real_T *b_temp;
  int b_i;
  emxInit_real_T(&temp, 1);
  emxInit_real_T(&nodeNum, 1);
  emxInit_real_T(&conn, 2);

  //       output:
  //       mesh          = object containing mesh data
  fileid = b_cfopen(filename, "rb");

  // open mesh file
  //  temp = fscanf(fid,'%i',2);   %read in number of nodes and number of elements 
  getSplitLine(static_cast<double>(fileid), temp);
  i = static_cast<int>(temp->data[0]);
  i1 = nodeNum->size[0];
  nodeNum->size[0] = i;
  emxEnsureCapacity_real_T(nodeNum, i1);
  i1 = static_cast<int>(temp->data[1]);
  i2 = conn->size[0] * conn->size[1];
  conn->size[0] = i1;
  conn->size[1] = 2;
  emxEnsureCapacity_real_T(conn, i2);
  loop_ub = i1 << 1;
  for (i2 = 0; i2 < loop_ub; i2++) {
    conn->data[i2] = 0.0;
  }

  emxInit_real_T(&x, 1);
  emxInit_real_T(&y, 1);
  emxInit_real_T(&z, 1);
  emxInit_real_T(&elNum, 1);
  i2 = elNum->size[0];
  elNum->size[0] = i1;
  emxEnsureCapacity_real_T(elNum, i2);
  i2 = x->size[0];
  x->size[0] = i;
  emxEnsureCapacity_real_T(x, i2);
  i2 = y->size[0];
  y->size[0] = i;
  emxEnsureCapacity_real_T(y, i2);
  i2 = z->size[0];
  z->size[0] = i;
  emxEnsureCapacity_real_T(z, i2);
  emxInit_real_T(&b_temp, 1);
  for (b_i = 0; b_i < i; b_i++) {
    //  read in node number and node coordinates
    getSplitLine(static_cast<double>(fileid), b_temp);
    nodeNum->data[b_i] = b_temp->data[0];
    x->data[b_i] = b_temp->data[1];
    y->data[b_i] = b_temp->data[2];
    z->data[b_i] = b_temp->data[3];
  }

  for (b_i = 0; b_i < i1; b_i++) {
    //  read in element number and connectivity list
    getSplitLine(static_cast<double>(fileid), b_temp);
    elNum->data[b_i] = b_temp->data[0];
    conn->data[b_i] = b_temp->data[2];
    conn->data[b_i + conn->size[0]] = b_temp->data[3];
  }

  emxFree_real_T(&b_temp);
  cfclose(static_cast<double>(fileid));

  // close mesh file
  i = mesh_nodeNum->size[0];
  mesh_nodeNum->size[0] = nodeNum->size[0];
  emxEnsureCapacity_real_T(mesh_nodeNum, i);
  loop_ub = nodeNum->size[0];
  for (i = 0; i < loop_ub; i++) {
    mesh_nodeNum->data[i] = nodeNum->data[i];
  }

  emxFree_real_T(&nodeNum);

  // store data in mesh object
  *mesh_numEl = temp->data[1];
  *mesh_numNodes = temp->data[0];
  i = mesh_x->size[0];
  mesh_x->size[0] = x->size[0];
  emxEnsureCapacity_real_T(mesh_x, i);
  loop_ub = x->size[0];
  emxFree_real_T(&temp);
  for (i = 0; i < loop_ub; i++) {
    mesh_x->data[i] = x->data[i];
  }

  emxFree_real_T(&x);
  i = mesh_y->size[0];
  mesh_y->size[0] = y->size[0];
  emxEnsureCapacity_real_T(mesh_y, i);
  loop_ub = y->size[0];
  for (i = 0; i < loop_ub; i++) {
    mesh_y->data[i] = y->data[i];
  }

  emxFree_real_T(&y);
  i = mesh_z->size[0];
  mesh_z->size[0] = z->size[0];
  emxEnsureCapacity_real_T(mesh_z, i);
  loop_ub = z->size[0];
  for (i = 0; i < loop_ub; i++) {
    mesh_z->data[i] = z->data[i];
  }

  emxFree_real_T(&z);
  i = mesh_elNum->size[0];
  mesh_elNum->size[0] = elNum->size[0];
  emxEnsureCapacity_real_T(mesh_elNum, i);
  loop_ub = elNum->size[0];
  for (i = 0; i < loop_ub; i++) {
    mesh_elNum->data[i] = elNum->data[i];
  }

  emxFree_real_T(&elNum);
  i = mesh_conn->size[0] * mesh_conn->size[1];
  mesh_conn->size[0] = conn->size[0];
  mesh_conn->size[1] = 2;
  emxEnsureCapacity_real_T(mesh_conn, i);
  loop_ub = conn->size[0] * conn->size[1];
  for (i = 0; i < loop_ub; i++) {
    mesh_conn->data[i] = conn->data[i];
  }

  emxFree_real_T(&conn);
}

//
// File trailer for readMesh.cpp
//
// [EOF]
//
