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


def gripper_pos(px,py,pz):
    return Matrix([[px],[py],[pz],[1]])

    #rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))

q0,q1,q2,q3,q4,q5,q6 = symbols('q0:7')
d0,d1,d2,d3,d4,d5,d6,d7 = symbols('d0:8')
a0,a1,a2,a3,a4,a5 = symbols('a0:6')
p_g = symbols('p_gx p_gy p_gz')

parameters = {a0:0,a1:0.35,a2:1.25,a3:-0.054,a4:0,a5:0,
     d0:0,d1:0.75,d2:0,d3:0,d4:1.5,d5:0,d6:0,d7:0.303,
     q0:0,q1:0,q2:0,q3:0,q4:0,q5:0,q6:0
    }

w_x, w_y, w_z = symbols('w_x w_y w_z')

alpha, beta, gamma = symbols('alpha beta gamma')
#w_tx, w_ty, w_tz = symbols('w_tx, w_ty, w_tz')

T1_0 = Matrix(
  [[ cos(q1), -sin(q1),   0,   0],
   [ sin(q1), -cos(q1),   0,   0],
   [       0,        0,   1,  d1],
   [       0,        0,   0,   1]])

T2_1 = Matrix(
  [[ cos(q2-pi/2), -sin(q2-pi/2),   0,  a1],
   [            0,             0,   1,   0],
   [-sin(q2-pi/2), -cos(q2-pi/2),   0,   0],
   [            0,             0,   0,   1]])

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

T4_1 = simplify(T2_1 * T3_2 * T4_3)
T6_0 = simplify(T1_0 * T2_1 * T3_2 * T4_3 * T5_4 * T6_5)

T6_4_euler = Matrix(
    [[cos(alpha)*cos(beta), cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma), cos(alpha)*sin(beta)*cos(gamma)+sin(alpha)*sin(gamma), 0],
     [sin(alpha)*cos(beta), sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma), sin(alpha)*sin(beta)*cos(gamma)-cos(alpha)*sin(gamma), 0],
     [          -sin(beta),                                    cos(beta)*sin(gamma),                                  cos(beta)*cos(gamma), 0],
     [                   0,                                                       0,                                                     0, 1]])

TG_4_euler = simplify(T6_4_euler * TG_6);
T4_G_euler = TG_4_euler.inv();


TG_3 = simplify(T4_3 * T5_4 * T6_5 * TG_6)
theta5_sym = solve(TG_3[1,2] - T6_4_euler[1,2], q5)
print(TG_3)
print(theta5_sym)


w_tx = cos(q1)*w_x + sin(q1)*w_y
w_ty = 0
w_tz = w_z - d1

p_gripper = Matrix([[p_g[0]], [p_g[1]], [p_g[2]], [1]])

T4_0 = simplify(T1_0 * T4_1)
p_0 = Matrix([[0],[0],[0],[1]])

#print(simplify(T2_1))

#print((T1_0*p_0).evalf(subs=parameters))
#print((T1_0*T2_1*p_0).evalf(subs=parameters))
#print((T1_0*T2_1*T3_2*p_0).evalf(subs=parameters))
#print((T1_0*T2_1*T3_2*T4_3*p_0).evalf(subs=parameters))

#print(T4_1)
rho1_sq = simplify((T4_1[0,3]-a1)**2+T4_1[2,3]**2)
#rho2_sq = simplify((T4_1[0,3]-a1)**2-T4_1[2,3]**2)
print(rho1_sq)
#print(rho2_sq)



#rho_sq = simplify(pow(T4_1[0,3]-a1,2) + pow(T4_1[2,3],2))
rho_sq = pow(w_tx-a1,2) + pow(w_tz,2)
sigma = (a2**2 + a3**2 + d4**2 - rho_sq) / (2*a2)

theta3_sym = asin(Min(1.0, Max(-1.0, sigma / sqrt(a3**2 + d4**2)))) - atan2(-a3,d4)
#theta3_sym = -2 * atan2(sqrt(d4**+a3**2-sigma**2) + d4, a3 - sigma)
theta2_2_sym = asin(d4/sqrt(rho_sq)*sin(q3+pi/2))
theta2_1_sym = atan2(w_tz,w_tx - a1)
theta2_sym = pi/2 - (theta2_1_sym + theta2_2_sym)

alpha_res = alpha - q1
beta_res = beta - q2 - q3


theta5_sym = acos(cos(alpha_res)*cos(beta_res))
theta6_sym = asin((-sin(alpha_res)*cos(gamma) + cos(alpha_res)*sin(beta_res)*sin(gamma)) / sin(q5))
theta4_sym = asin(sin(alpha_res)*cos(beta_res) / sin(q5))

w_t = Matrix([[w_tx],[w_ty],[w_tz],[1]])

w_tt = simplify(T4_3.inv()*w_t)
print(w_tt)
    #print(rho_sq)
    #print(rho2_sq)
    #print(solVve(rho_sq - rho2_sq, q3, manual=True))



    #print(T4_G_euler * p_gripper)

    #print(solve(-sin(q1)*w_x + cos(q1)*w_y, q1))

    #print(T4_1)
    #print(w_tx)
    #print(w_ty)
    #print(w_tz)
    #print(T4_G_euler)

    # handle_calculate_IK.globals = {
    #   "q":q, "d":d, "a":a,
    #   "alpha":alpha, "beta":beta, "gamma":gamma,
    #   "w_x":w_x, "w_y":w_y, "w_z":w_z,
    #   "w_tx":w_tx, "w_ty":w_ty, "w_tz":w_tz,
    #   "T4_1":T4_1, "T4_G_euler":T4_G_euler,
    #   "rho_sq":rho_sq, "sigma":sigma,
    #   "theta2_sym":theta2_sym,
    #   "theta3_sym":theta3_sym,
    #   "parameters":parameters
    #   }


def inverse_kinematics_wrist(wx,wy,wz):
    theta1 = atan2(wy,wx)
    print("theta1:", theta1)

    print("rho_sq:", rho_sq.evalf(subs=parameters))
    parameters[w_x] = wx
    parameters[w_y] = wy
    parameters[w_z] = wz
    print("rho_sq:", rho_sq.evalf(subs=parameters))
    parameters[q1] = theta1
    print("sigma:", sigma.evalf(subs=parameters))
    print("sigma/sqrt(a3**2+d4**2):", (sigma/sqrt(a3**2+d4**2)).evalf(subs=parameters))

    print(theta3_sym)
    theta3 = theta3_sym.evalf(subs=parameters)
    print("theta3:", theta3)

    parameters[q3] = theta3
    theta2 = theta2_sym.evalf(subs=parameters)

    print("theta2:", theta2)
    print("w_t:", w_t.evalf(subs=parameters))

    return theta1,theta2,theta3


def inverse_kinematics_gripper():
    theta5 = theta5_sym.evalf(subs=parameters)
    parameters[q5] = theta5

    if theta5 < 0.01:
        theta4 = 0.0
        theta6 = parameters[gamma]
        parameters[q4] = theta4
        parameters[q6] = theta6
    else:
        theta6 = theta6_sym.evalf(subs=parameters)
        parameters[q6] = theta6
        theta4 = theta4_sym.evalf(subs=parameters)
        parameters[q4] = theta4

    return theta4,theta5,theta6



def forward_kinematics(angles):
    parameters[q1] = angles[0]
    parameters[q2] = angles[1]
    parameters[q3] = angles[2]

    if len(angles) == 6:
        parameters[q4] = angles[3]
        parameters[q5] = angles[4]
        parameters[q6] = angles[5]
        return (T6_0*p_0).evalf(subs=parameters)
    else:
        return (T4_0*p_0).evalf(subs=parameters)


def handle_calculate_IK(req):
    if len(req.poses) < 1:
        print "No valid poses received"
        return -1
    else:
        # Initialize service response
        joint_trajectory_list = []
        # T4_G_euler = g["T4_G_euler"]
        # alpha = g["alpha"]
        # beta = g["beta"]
        # gamma = g["gamma"]
        # d0,d1,d2,d3,d4,d5,d6,d7 = g["d"]
        # w_x, w_y, w_z = g["w_x"], g["w_y"], g["w_z"]
        # q0,q1,q2,q3,q4,q5,q6 = g["q"]
        # parameters = g["parameters"]
        # theta2_sym = g["theta2_sym"]
        # theta3_sym = g["theta3_sym"]
        # w_tx, w_ty, w_tz = g["w_tx"], g["w_ty"], g["w_tz"]
        # w_t = Matrix([[w_tx], [w_ty], [w_tz], [1]])

        for x in xrange(0, len(req.poses)):
            px = req.poses[x].position.x
            py = req.poses[x].position.y
            pz = req.poses[x].position.z

            (roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
                [req.poses[x].orientation.x, req.poses[x].orientation.y,
                    req.poses[x].orientation.z, req.poses[x].orientation.w])

            parameters[gamma] = roll
            parameters[beta] = pitch
            parameters[alpha] = yaw

            p_gripper = gripper_pos(px,py,pz)
            wx_val,wy_val,wz_val,ww_val = simplify((T4_G_euler * p_gripper).evalf(subs=parameters))
            theta1,theta2,theta3 = inverse_kinematics_wrist(wx_val, wy_val, wz_val)
            theta4,theta5,theta6 = inverse_kinematics_gripper()

            print "Tool position :", px, py, pz
            print "Wrist position:", wx_val, wy_val, wz_val
            print "Resulting wrist position: ", forward_kinematics([theta1, theta2, theta3])
            print "Resulting Tool position:", forward_kinematics([theta1, theta2, theta3, theta4, theta5, theta6])

            joint_trajectory_point = JointTrajectoryPoint()
            joint_trajectory_point.positions = [theta1,theta2,theta3,theta4,theta5,theta6]
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
    # print(forward_kinematics([0,0,0,0,0,0]))
    # wx = 1.85
    # wy = 0
    # wz = 1.946
    #
    # wx,wy,wz = 2.20378066044923, -0.217688238728822, 0.263026020380731
    if False:
        import math
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
                print "w: ", wx,wy,wz
                print "Angles: ", angles
                print "w_pred:", w_pred
                print "RMSE: ", rmse.get()
                print
    else:
        IK_server()
