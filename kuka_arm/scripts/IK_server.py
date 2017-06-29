#!/usr/bin/env python

# Copyright (C) 2017 Electric Movement Inc.
#
# This file is part of Robotic Arm: Pick and Place project for Udacity
# Robotics nano-degree program
#
# All Rights Reserved.

# Author: Harsh Pandya

# import modules
import rospy
import tf
from kuka_arm.srv import *
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint
from geometry_msgs.msg import Pose
from mpmath import *
from sympy import *


def handle_calculate_IK(req):
    #rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))

    if handle_calculate_IK.globals is None:
        handle_calculate_IK.globals = {}
        q = symbols('q0:7')
        d = symbols('d0:8')
        a = symbols('a0:6')

        w_x, w_y, w_z = symbols('w_x w_y w_z')
        #w_tx, w_ty, w_tz = symbols('w_tx, w_ty, w_tz')

        T1_0 = Matrix([[ cos(q[1]), -sin(q[1]),        0,        0],
                       [ sin(q[1]), -cos(q[1]),        0,        0],
                       [       0,        0,        1,     d[1]],
                       [       0,        0,        0,       1]])

        T2_1 = Matrix([[ cos(q[2]), -sin(q[2]),        0,     a[1]],
                       [       0,        0,        1,        0],
                       [-sin(q[2]), -cos(q[2]),        0,        0],
                       [       0,        0,        0,        1]])

        T3_2 = Matrix([[ cos(q[3] + pi/2), -sin(q[3] + pi/2),        0,        a[2]],
                       [ sin(q[3] + pi/2),  cos(q[3] + pi/2),        0,           0],
                       [              0,               0,        1,           0],
                       [              0,               0,        0,           1]])

        T4_3 = Matrix([[ cos(q[4]), -sin(q[4]),        0,    a[3]],
                       [       0,        0,        1,        d[4]],
                       [-sin(q[4]), -cos(q[4]),        0,       0],
                       [       0,        0,        0,           1]])

        T5_4 = Matrix([[ cos(q[5]), -sin(q[5]),        0,       0],
                       [       0,        0,       -1,       0],
                       [ sin(q[5]),  cos(q[5]),        0,       0],
                       [       0,        0,        0,       1]])

        T6_5 = Matrix([[ cos(q[6]), -sin(q[6]),        0,       0],
                       [       0,        0,        1,       0],
                       [-sin(q[6]), -cos(q[6]),        0,       0],
                       [       0,        0,        0,       1]])

        T6_G = Matrix([[     1,     0,     0,    0],
                       [     0,     1,     0,    0],
                       [     0,     0,     1, d[7]],
                       [     0,     0,     0,    1]])

        T4_1 = simplify(T2_1 * T3_2 * T4_3)

        w_tx = cos(q[1])*w_x + sin(q[1])*w_y
        w_ty = 0
        w_tz = w_z - d[1]

        print(T4_1, w_tx, w_ty, w_tz)

        return


    if len(req.poses) < 1:
        print "No valid poses received"
        return -1
    else:
        # Initialize service response
        joint_trajectory_list = []
        for x in xrange(0, len(req.poses)):
            # IK code starts here
            joint_trajectory_point = JointTrajectoryPoint()

            # Define DH param symbols



            # Joint angle symbols



            # Modified DH params



            # Define Modified DH Transformation matrix



            # Create individual transformation matrices



            # Extract end-effector position and orientation from request
	    # px,py,pz = end-effector position
	    # roll, pitch, yaw = end-effector orientation
            px = req.poses[x].position.x
            py = req.poses[x].position.y
            pz = req.poses[x].position.z

            (roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
                [req.poses[x].orientation.x, req.poses[x].orientation.y,
                    req.poses[x].orientation.z, req.poses[x].orientation.w])

            # Calculate joint angles using Geometric IK method




            # Populate response for the IK request
            # In the next line replace theta1,theta2...,theta6 by your joint angle variables
	    joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5, theta6]
	    joint_trajectory_list.append(joint_trajectory_point)

        rospy.loginfo("length of Joint Trajectory List: %s" % len(joint_trajectory_list))
        return CalculateIKResponse(joint_trajectory_list)


handle_calculate_IK.globals = None


def IK_server():
    # initialize node and declare calculate_ik service
    rospy.init_node('IK_server')
    s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
    print "Ready to receive an IK request"
    rospy.spin()

if __name__ == "__main__":
    handle_calculate_IK(None)
    #IK_server()
