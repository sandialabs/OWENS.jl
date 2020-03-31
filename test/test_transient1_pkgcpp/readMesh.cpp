//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: readMesh.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 31-Mar-2020 09:49:33
//

// Include Files
#include "readMesh.h"
#include "fileManager.h"
#include "getSplitLine.h"
#include "rt_nonfinite.h"
#include "test_transient1.h"
#include "test_transient1_emxutil.h"
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
  signed char fileid;
  int i;
  int i1;
  int b_i;
  int loop_ub;
  emxArray_real_T *b_temp;
  emxInit_real_T(&temp, 1);

  //       output:
  //       mesh          = object containing mesh data
  fileid = c_cfopen(filename, "rb");

  // open mesh file
  //  temp = fscanf(fid,'%i',2);   %read in number of nodes and number of elements 
  getSplitLine(static_cast<double>(fileid), temp);
  i = static_cast<int>(temp->data[0]);
  i1 = mesh_nodeNum->size[0];
  mesh_nodeNum->size[0] = i;
  emxEnsureCapacity_real_T(mesh_nodeNum, i1);
  i1 = static_cast<int>(temp->data[1]);
  b_i = mesh_conn->size[0] * mesh_conn->size[1];
  mesh_conn->size[0] = i1;
  mesh_conn->size[1] = 2;
  emxEnsureCapacity_real_T(mesh_conn, b_i);
  loop_ub = i1 << 1;
  for (b_i = 0; b_i < loop_ub; b_i++) {
    mesh_conn->data[b_i] = 0.0;
  }

  b_i = mesh_elNum->size[0];
  mesh_elNum->size[0] = i1;
  emxEnsureCapacity_real_T(mesh_elNum, b_i);
  b_i = mesh_x->size[0];
  mesh_x->size[0] = i;
  emxEnsureCapacity_real_T(mesh_x, b_i);
  b_i = mesh_y->size[0];
  mesh_y->size[0] = i;
  emxEnsureCapacity_real_T(mesh_y, b_i);
  b_i = mesh_z->size[0];
  mesh_z->size[0] = i;
  emxEnsureCapacity_real_T(mesh_z, b_i);
  emxInit_real_T(&b_temp, 1);
  for (b_i = 0; b_i < i; b_i++) {
    //  read in node number and node coordinates
    getSplitLine(static_cast<double>(fileid), b_temp);
    mesh_nodeNum->data[b_i] = b_temp->data[0];
    mesh_x->data[b_i] = b_temp->data[1];
    mesh_y->data[b_i] = b_temp->data[2];
    mesh_z->data[b_i] = b_temp->data[3];
  }

  for (b_i = 0; b_i < i1; b_i++) {
    //  read in element number and connectivity list
    getSplitLine(static_cast<double>(fileid), b_temp);
    mesh_elNum->data[b_i] = b_temp->data[0];
    mesh_conn->data[b_i] = b_temp->data[2];
    mesh_conn->data[b_i + mesh_conn->size[0]] = b_temp->data[3];
  }

  emxFree_real_T(&b_temp);
  cfclose(static_cast<double>(fileid));

  // close mesh file
  // store data in mesh object
  *mesh_numEl = temp->data[1];
  *mesh_numNodes = temp->data[0];
  emxFree_real_T(&temp);
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
  signed char fileid;
  int i;
  int i1;
  int b_i;
  int loop_ub;
  emxArray_real_T *b_temp;
  emxInit_real_T(&temp, 1);

  //       output:
  //       mesh          = object containing mesh data
  fileid = b_cfopen(filename, "rb");

  // open mesh file
  //  temp = fscanf(fid,'%i',2);   %read in number of nodes and number of elements 
  getSplitLine(static_cast<double>(fileid), temp);
  i = static_cast<int>(temp->data[0]);
  i1 = mesh_nodeNum->size[0];
  mesh_nodeNum->size[0] = i;
  emxEnsureCapacity_real_T(mesh_nodeNum, i1);
  i1 = static_cast<int>(temp->data[1]);
  b_i = mesh_conn->size[0] * mesh_conn->size[1];
  mesh_conn->size[0] = i1;
  mesh_conn->size[1] = 2;
  emxEnsureCapacity_real_T(mesh_conn, b_i);
  loop_ub = i1 << 1;
  for (b_i = 0; b_i < loop_ub; b_i++) {
    mesh_conn->data[b_i] = 0.0;
  }

  b_i = mesh_elNum->size[0];
  mesh_elNum->size[0] = i1;
  emxEnsureCapacity_real_T(mesh_elNum, b_i);
  b_i = mesh_x->size[0];
  mesh_x->size[0] = i;
  emxEnsureCapacity_real_T(mesh_x, b_i);
  b_i = mesh_y->size[0];
  mesh_y->size[0] = i;
  emxEnsureCapacity_real_T(mesh_y, b_i);
  b_i = mesh_z->size[0];
  mesh_z->size[0] = i;
  emxEnsureCapacity_real_T(mesh_z, b_i);
  emxInit_real_T(&b_temp, 1);
  for (b_i = 0; b_i < i; b_i++) {
    //  read in node number and node coordinates
    getSplitLine(static_cast<double>(fileid), b_temp);
    mesh_nodeNum->data[b_i] = b_temp->data[0];
    mesh_x->data[b_i] = b_temp->data[1];
    mesh_y->data[b_i] = b_temp->data[2];
    mesh_z->data[b_i] = b_temp->data[3];
  }

  for (b_i = 0; b_i < i1; b_i++) {
    //  read in element number and connectivity list
    getSplitLine(static_cast<double>(fileid), b_temp);
    mesh_elNum->data[b_i] = b_temp->data[0];
    mesh_conn->data[b_i] = b_temp->data[2];
    mesh_conn->data[b_i + mesh_conn->size[0]] = b_temp->data[3];
  }

  emxFree_real_T(&b_temp);
  cfclose(static_cast<double>(fileid));

  // close mesh file
  // store data in mesh object
  *mesh_numEl = temp->data[1];
  *mesh_numNodes = temp->data[0];
  emxFree_real_T(&temp);
}

//
// File trailer for readMesh.cpp
//
// [EOF]
//
