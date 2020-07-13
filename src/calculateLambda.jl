function [lambda] = calculateLambda(theta1,theta2,theta3)
#calculateLambda Calculates transformation matrix from element to hub frame
# **********************************************************************
# *                   Part of the SNL OWENS Toolkit                    *
# * Developed by Sandia National Laboratories Wind Energy Technologies *
# *             See license.txt for disclaimer information             *
# **********************************************************************
#   [lambda] = calculateLambda(theta1,theta2,theta3 )
#
#   This function calculates a transformation matrix to transform the
#   element degree of freedom vector (12 DOFs) from the hub frame to
#   the element frame. The transformation matrix is constructed via the
#   direction cosine matrices of a 3-2-1 Euler rotation sequence.
#
#      input:
#      theta1        = angle (rad) of rotation for 1st rotation
#                      of 3-2-1 sequence
#      theta2        = angle (rad) of rotation for 2nd rotation
#                      of 3-2-1 sequence
#      theta3        = angle (rad) of rotation for 3rd rotation
#                      of 3-2-1 sequence

#      output:
#      lambda        = 12 x 12 transformation matrix

# dcm that is created is [dcm] = [M1(theta3)][M2(theta2)][M3(theta1)]

[dcm] = calculateLambdaSlim(theta1,theta2,theta3);


lambda = zeros(12);
lambda(1:3,1:3) = dcm;
lambda(4:6,4:6) = dcm;
lambda(7:9,7:9) = dcm;
lambda(10:12,10:12) = dcm;

end

