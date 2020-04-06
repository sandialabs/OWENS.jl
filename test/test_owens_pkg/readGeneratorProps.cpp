//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: readGeneratorProps.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 06-Apr-2020 16:48:15
//

// Include Files
#include "readGeneratorProps.h"
#include "fileManager.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include <string.h>

// Function Definitions

//
// readGeneratorProps reads generator properties from file
//  **********************************************************************
//  *                   Part of the SNL OWENS Toolkit                    *
//  * Developed by Sandia National Laboratories Wind Energy Technologies *
//  *             See license.txt for disclaimer information             *
//  **********************************************************************
//    [genprops] = readGeneratorProps(generatorfilename)
//
//    This function reads generator properties from file.
//
//    input:
//    generatorfilenanme  = string containing generator property file name
//
//    output:
//    genprops          = model object containing generator properties
// Arguments    : const emxArray_char_T *generatorfilename
// Return Type  : double
//
double readGeneratorProps(const emxArray_char_T *generatorfilename)
{
  double genprops;
  boolean_T b_bool;
  int fid;
  signed char fileid;
  int exitg1;
  static const char b_cv[3] = { 'a', 'l', 'l' };

  b_bool = false;
  if ((generatorfilename->size[0] == 1) && (generatorfilename->size[1] == 3)) {
    fid = 0;
    do {
      exitg1 = 0;
      if (fid < 3) {
        if (generatorfilename->data[fid] != b_cv[fid]) {
          exitg1 = 1;
        } else {
          fid++;
        }
      } else {
        b_bool = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }

  if (b_bool) {
    fid = 0;
  } else {
    fileid = b_cfopen(generatorfilename, "rb");
    fid = fileid;
  }

  // open generator property file
  if (fid == -1) {
    genprops = 0.0;

    // if generator property file does not exist, set object to null
  } else {
    // if file can be opened
    //          genprops.ratedTorque = fscanf(fid,'%f',1); %store rated torque
    //          dum = fgetl(fid);
    //          genprops.zeroTorqueGenSpeed = fscanf(fid,'%f',1); %store zero torque generator zpeed 
    //          dum = fgetl(fid);
    //          genprops.pulloutRatio = fscanf(fid,'%f',1); %store pullout ratio 
    //          dum = fgetl(fid);
    //          genprops.ratedGenSlipPerc= fscanf(fid,'%f',1); %store rated generator slip percentage 
    //          dum = fgetl(fid);
    //
    //          fclose(fid); %close generator propery file
  }

  return genprops;
}

//
// File trailer for readGeneratorProps.cpp
//
// [EOF]
//
