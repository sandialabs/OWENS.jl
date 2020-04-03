//
// Trial License - for use to evaluate programs for possible purchase as
// an end-user only.
// File: main.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 02-Apr-2020 17:04:27
//

//***********************************************************************
// This automatically generated example C++ main file shows how to call
// entry-point functions that MATLAB Coder generated. You must customize
// this file for your application. Do not modify this file directly.
// Instead, make a copy of this file, modify it, and integrate it into
// your development environment.
//
// This file initializes entry-point function arguments to a default
// size and value before calling the entry-point functions. It does
// not store or use any values returned from the entry-point functions.
// If necessary, it does pre-allocate memory for returned values.
// You can use this file as a starting point for a main function that
// you can deploy in your application.
//
// After you copy the file, and before you deploy it, you must make the
// following changes:
// * For variable-size function arguments, change the example sizes to
// the sizes that your application requires.
// * Change the example values of function arguments to the values that
// your application requires.
// * If the entry-point functions return values, store these values or
// otherwise use them as required by your application.
//
//***********************************************************************

// Include Files
#include "main.h"
#include "rt_nonfinite.h"
#include "test_owens.h"
#include "test_owens_terminate.h"
#include <string.h>

// Function Declarations
static boolean_T argInit_boolean_T();
static void main_test_owens();

// Function Definitions

//
// Arguments    : void
// Return Type  : boolean_T
//
static boolean_T argInit_boolean_T()
{
  return true;
}

//
// Arguments    : void
// Return Type  : void
//
static void main_test_owens()
{
  boolean_T test_transient_tmp;

  // Initialize function 'test_owens' input arguments.
  test_transient_tmp = argInit_boolean_T();

  // Call the entry-point 'test_owens'.
  test_owens(test_transient_tmp, test_transient_tmp);
}

//
// Arguments    : int argc
//                const char * const argv[]
// Return Type  : int
//
int main(int, const char * const [])
{
  // The initialize function is being called automatically from your entry-point function. So, a call to initialize is not included here.
  // Invoke the entry-point functions.
  // You can call entry-point functions multiple times.
  main_test_owens();

  // Terminate the application.
  // You do not need to do this more than one time.
  test_owens_terminate();
  return 0;
}

//
// File trailer for main.cpp
//
// [EOF]
//
