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
import numpy as np
import math


# Class for calculating the RMSE of the predicted gripper position
class RMSE:
    def __init__(self):
        self.acc = None
        self.n = 0

    def add(self, x, x_gt):
        x = np.array((x[0], x[1], x[2]), np.float64)
        x_gt = np.array((x_gt[0], x_gt[1], x_gt[2]), np.float64)

        if self.acc is None:
            self.acc = np.zeros_like(x_gt)

        self.acc += np.square(x - x_gt)
        self.n += 1.0

    def get(self):
        return np.sqrt(self.acc / self.n)


def vec(m):
    return np.array((m[0], m[1], m[2]))

def make_sympy_vector(px,py,pz):
    return Matrix([[px],[py],[pz],[1]])


# Some helper functions
def safe_asin(expr):
    return asin(Max(-1.0, Min(1.0, expr)))

def normalize(expr):
    return atan2(sin(expr), cos(expr))

def normalize_vector(v):
    s = sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    return v / s


# Define the symbols for the DH parameters:
q0,q1,q2,q3,q4,q5,q6 = symbols('q0:7')
d0,d1,d2,d3,d4,d5,d6,d7 = symbols('d0:8')
a0,a1,a2,a3,a4,a5 = symbols('a0:6')


# Define the static DH parameters:
parameters = {a0:0,a1:0.35,a2:1.25,a3:-0.054,a4:0,a5:0,
     d0:0,d1:0.75,d2:0,d3:0,d4:1.5,d5:0,d6:0,d7:0.303,
     q0:0,q1:0,q2:0,q3:0,q4:0,q5:0,q6:0
    }


# Symbols for the gripper position:
g_x, g_y, g_z = symbols('g_x g_y g_z')


# Symbols for the wrist position:
w_x, w_y, w_z = symbols('w_x w_y w_z')


# Symbols for rotation matrices:
phi_x, phi_y, phi_z = symbols('phi_x phi_y phi_z')


# Atomic rotation matrices:
R_x = Matrix(
    [[ 1,          0,          0,  0 ],
     [ 0,  cos(phi_x), -sin(phi_x),  0 ],
     [ 0,  sin(phi_x),  cos(phi_x),  0 ],
     [ 0,          0,          0,  1 ]])

R_y = Matrix(
    [[ cos(phi_y), 0, sin(phi_y), 0 ],
     [          0, 1,          0, 0 ],
     [-sin(phi_y), 0, cos(phi_y), 0 ],
     [          0, 0,          0, 1 ]])

R_z = Matrix(
    [[ cos(phi_z), -sin(phi_z), 0, 0 ],
     [ sin(phi_z),  cos(phi_z), 0, 0 ],
     [        0,         0, 1, 0 ],
     [        0,         0, 0, 1 ]])



# Components of the correction matrix from DH to URDF frame:
R_corr_z = Matrix(
    [[cos(pi), -sin(pi), 0, 0 ],
     [sin(pi),  cos(pi), 0, 0 ],
     [      0,        0, 1, 0 ],
     [      0,        0, 0, 1 ]])

R_corr_y = Matrix(
    [[ cos(-pi/2), 0, sin(-pi/2), 0 ],
     [          0, 1,          0, 0 ],
     [-sin(-pi/2), 0, cos(-pi/2), 0 ],
     [          0, 0,          0, 1 ]])

R_corr = R_corr_z * R_corr_y
R_zyx = simplify(R_z * R_y * R_x)


# Individiual link frame transformation matrices:
T1_0 = Matrix(
  [[ cos(q1), -sin(q1),   0,   0],
   [ sin(q1),  cos(q1),   0,   0],
   [       0,        0,   1,  d1],
   [       0,        0,   0,   1]])

T2_1 = Matrix(
  [[ cos(q2-pi/2), -sin(q2-pi/2),   0,  a1],
   [            0,             0,   1,   0],
   [-sin(q2-pi/2), -cos(q2-pi/2),   0,   0],
   [            0,             0,   0,  1]])

T3_2 = Matrix(
  [[ cos(q3), -sin(q3),  0,  a2],
   [ sin(q3),  cos(q3),  0,   0],
   [       0,        0,  1,   0],
   [       0,        0,  0,   1]])

T4_3 = Matrix(
  [[ cos(q4), -sin(q4),  0, a3],
   [       0,        0,  1, d4],
   [-sin(q4), -cos(q4),  0,  0],
   [       0,        0,  0,  1]])

T5_4 = Matrix(
  [[ cos(q5), -sin(q5),  0,  0],
   [       0,        0, -1,  0],
   [ sin(q5),  cos(q5),  0,  0],
   [       0,        0,  0,  1]])

T6_5 = Matrix(
  [[ cos(q6), -sin(q6),  0,  0],
   [       0,        0,  1,  0],
   [-sin(q6), -cos(q6),  0,  0],
   [       0,        0,  0,  1]])

TG_6 = Matrix(
  [[ 1,  0,  0,  0],
   [ 0,  1,  0,  0],
   [ 0,  0,  1, d7],
   [ 0,  0,  0,  1]])

T3_0 = simplify(T1_0 * T2_1 * T3_2)
T4_1 = simplify(T2_1 * T3_2 * T4_3)
T4_0 = simplify(T1_0 * T4_1)
T6_0 = T1_0 * T2_1 * T3_2 * T4_3 * T5_4 * T6_5
T6_3 = simplify(T4_3 * T5_4 * T6_5)
T_total = T6_0 * TG_6 * R_corr

# Axis along the gripper tool:
N_gripper = R_zyx.col(2)
gripper_pos = Matrix([g_x, g_y, g_z, 1])
wrist_pos = gripper_pos - N_gripper * d7

# Transformed wrist center position:
w_tx = cos(q1)*w_x + sin(q1)*w_y
w_ty = 0
w_tz = w_z - d1


# Zero vector for forward kinematics:
p_0 = Matrix([[0],[0],[0],[1]])


# Expressions for calculating q3:
rho_sq = pow(w_tx - a1,2) + pow(w_tz,2)
sigma_sq = a2**2 + a3**2 + d4**2
theta3_sym = pi - safe_asin((rho_sq - sigma_sq)/(2*a2*sqrt(a3**2+d4**2))) - atan2(a3, -d4)


# Expressions for calculating q2:
theta2_1_sym = atan2(w_tz,w_tx - a1)
theta2_2_sym = -atan2(-a3*sin(q3)-d4*cos(q3), a2 + a3*cos(q3)-d4*sin(q3))
theta2_sym = pi/2 - theta2_1_sym - theta2_2_sym


# Objects for accumulating RMSE values:
rmse_wrist = RMSE()
rmse_gripper = RMSE()


def inverse_kinematics_wrist(wrist_pos):
    wx, wy, wz, ww = wrist_pos
    theta1 = atan2(wy,wx)
    print("theta1:", theta1)

    parameters[w_x] = wx
    parameters[w_y] = wy
    parameters[w_z] = wz
    parameters[q1] = theta1
    parameters[q4] = 0.0

    theta3 = normalize(theta3_sym.evalf(subs=parameters))
    print("theta3:", theta3)

    parameters[q3] = theta3
    theta2 = theta2_sym.evalf(subs=parameters)

    print("theta2:", theta2)
    parameters[q2] = theta2

    return theta1,theta2,theta3


def inverse_kinematics_gripper(R_gripper):
    R_gripper = T3_0.inv().evalf(subs=parameters) * R_gripper

    theta5 = math.acos(R_gripper[1,2])
    print "theta5:", theta5
    parameters[q5] = theta5

    if theta5 < 0.01:
        theta4 = 0.0
        theta6 = math.acos(R_gripper[0,0])
        parameters[q4] = theta4
        parameters[q6] = theta6
    else:
        theta6 = math.acos(R_gripper[1,0]/sin(theta5))
        parameters[q6] = theta6
        theta4 = math.atan2(R_gripper[2,2], -R_gripper[0,2])
        parameters[q4] = theta4

    return theta4,theta5,theta6


def forward_kinematics(angles):
    parameters[q1] = angles[0]
    parameters[q2] = angles[1]
    parameters[q3] = angles[2]
    parameters[q4] = 0.0

    if len(angles) == 6:
        parameters[q4] = angles[3]
        parameters[q5] = angles[4]
        parameters[q6] = angles[5]
        return (T_total*p_0).evalf(subs=parameters)
    else:
        return (T1_0 * T2_1 * T3_2 * T4_3 * p_0).evalf(subs=parameters)


def calc_gripper_rotation_matrix(quaternion):
    if False:
        return Matrix(tf.transformations.quaternion_matrix(quaternion)) * R_corr
    else:
        yaw_val, pitch_val, roll_val = tf.transformations.euler_from_quaternion(quaternion, 'rzyx')
        return (R_zyx * R_corr).evalf(subs={phi_x:roll_val, phi_y:pitch_val, phi_z:yaw_val})


def calc_wrist_position(grapper_pos, R):
    return (grapper_pos - R.col(2) * d7).evalf(subs=parameters)


def handle_calculate_IK(req):
    if len(req.poses) < 1:
        print "No valid poses received"
        return -1
    else:
        joint_trajectory_list = []

        for x in xrange(0, len(req.poses)):
            gripper_pos = Matrix([req.poses[x].position.x, req.poses[x].position.y, req.poses[x].position.z, 1.0])
            quaternion = Matrix([[req.poses[x].orientation.x, req.poses[x].orientation.y,
                req.poses[x].orientation.z, req.poses[x].orientation.w]])

            R_gripper = calc_gripper_rotation_matrix(quaternion)
            wrist_pos = calc_wrist_position(gripper_pos, R_gripper)
            theta1,theta2,theta3 = inverse_kinematics_wrist(wrist_pos)
            theta4,theta5,theta6 = inverse_kinematics_gripper(R_gripper)

            wrist_pos_pred = forward_kinematics([theta1, theta2, theta3])
            gripper_pos_pred = forward_kinematics([theta1, theta2, theta3, theta4, theta5, theta6])

            print "Tool position :", gripper_pos
            print "Wrist position:", wrist_pos
            print "Resulting wrist position: ", wrist_pos_pred
            print "Resulting Tool position:", gripper_pos_pred

            rmse_wrist.add(wrist_pos_pred, wrist_pos)
            rmse_gripper.add(gripper_pos_pred, gripper_pos)

            print "RMSE wrist:", rmse_wrist.get()
            print "RMSE gripper:", rmse_gripper.get()

            joint_trajectory_point = JointTrajectoryPoint()
            joint_trajectory_point.positions = [theta1,theta2,theta3,theta4,theta5,theta6]
            joint_trajectory_list.append(joint_trajectory_point)

        rospy.loginfo("length of Joint Trajectory List: %s" % len(joint_trajectory_list))
        return CalculateIKResponse(joint_trajectory_list)


def IK_server():
    # initialize node and declare calculate_ik service
    rospy.init_node('IK_server')
    s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
    print "Ready to receive an IK request"
    rospy.spin()

if __name__ == "__main__":
    TESTING = False
    if TESTING:
        test_cases = [
            (
                [2.16135,-1.42635,1.55109],
                [0.708611,0.186356,-0.157931,0.661967],
                [1.89451,-1.44302,1.69366],
                [-0.65,0.45,-0.36,0.95,0.79,0.49]
            ),
            (
                [-0.56754,0.93663,3.0038],
                [0.62073, 0.48318,0.38759,0.480629],
                [-0.638,0.64198,2.9988],
                [-0.79,-0.11,-2.33,1.94,1.14,-3.68]
            ),
            (
                [-1.3863,0.02074,0.90986],
                [0.01735,-0.2179,0.9025,0.371016],
                [-1.1669,-0.17989,0.85137],
                [-2.99,-0.12,0.94,4.06,1.29,-4.12]
            )
        ]

        for t in test_cases:
            print "="*60
            (px,py,pz),(qx,qy,qz,qw),(wx,wy,wz),angles = t

            gripper_pos = Matrix([px, py, pz, 1.0])
            quaternion = Matrix([qx, qy, qz, qw])

            print "Gripper axis:", normalize_vector(Matrix([px, py, pz]) - Matrix([wx, wy, wz]))
            R_gripper = calc_gripper_rotation_matrix(quaternion)
            print "R_gripper:", R_gripper
            wrist_pos = calc_wrist_position(gripper_pos, R_gripper)
            angles_pred = inverse_kinematics_wrist(wrist_pos) + inverse_kinematics_gripper(R_gripper)
            parameters[q1] = angles[0]
            parameters[q2] = angles[1]
            parameters[q3] = angles[2]
            parameters[q4] = angles[3]
            parameters[q5] = angles[4]
            parameters[q6] = angles[5]
            print "T6_0  :", T6_0.evalf(subs=parameters)
            print "p     :", (px,py,pz)
            print "w     :", (wx,wy,wz)
            print "wrist_pos:", wrist_pos
            print "p_forward:", forward_kinematics(angles)
            print "w_pred", forward_kinematics(angles_pred[0:3])
            print "p_pred", forward_kinematics(angles_pred)
            print "Angles:", angles
            print "Angles pred:", angles_pred
            print
            print

        rmse = RMSE()
        phi = 0
        r = 3.0
        z = 1.5
        for z in np.arange(1.0, 1.51, 0.25):
            for phi in np.arange(0.0, 2.0 * math.pi, 0.25 * math.pi):
                wx,wy,wz = r * math.cos(phi), r * math.sin(phi), z
                angles = inverse_kinematics_wrist(wx,wy,wz)
                w_pred = vec(forward_kinematics(angles))
                rmse.add(w_pred, (wx,wy,wz))
                print "="*70
                print "z = %.2f, phi = %.2f" % (z,phi)
                print "="*70
                print "T4_0: ", T4_0
                print "w: ", wx,wy,wz
                print "Angles: ", angles
                print "w_pred:", w_pred
                print "RMSE: ", rmse.get()
                print
    else:
        IK_server()
