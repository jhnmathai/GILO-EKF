// MIT License
//
// Copyright (c) 2024 Yibin Wu, Tiziano Guadagnino
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
#include "imuPropagation.hpp"
#include "rotation.hpp"
#include <sophus/so3.hpp>

namespace lio_ekf {

void insMechanization(const BodyState &pvapre, BodyState &pvacur,
                      const IMU &imupre, const IMU &imucur) {

  velUpdate(pvapre, pvacur, imupre, imucur);
  posUpdate(pvapre, pvacur, imupre, imucur);
  attUpdate(pvapre, pvacur, imupre, imucur);
}

void velUpdate(const BodyState &pvapre, BodyState &pvacur, const IMU &imupre, const IMU &imucur) {

    Eigen::Vector3d d_vfb, d_vfn, d_vgn, gl, midvel, midpos;
    Eigen::Vector3d temp1, temp2, temp3;
    Eigen::Matrix3d cnn, I33 = Eigen::Matrix3d::Identity();
    Eigen::Quaterniond qne, qee, qnn, qbb, q1, q2;

    // 计算地理参数，子午圈半径和卯酉圈半径，地球自转角速度投影到n系, n系相对于e系转动角速度投影到n系，重力值
    // calculate geographic parameters, Meridian and Mao unitary radii,
    // earth rotational angular velocity projected to n-frame,
    // rotational angular velocity of n-frame to e-frame projected to n-frame, and gravity
    Eigen::Vector2d rmrn = Earth::meridianPrimeVerticalRadius(pvapre.pose.translation()(0));
    Eigen::Vector3d wie_n, wen_n;
    wie_n << WGS84_WIE * cos(pvapre.pose.translation()[0]), 0, -WGS84_WIE * sin(pvapre.pose.translation()[0]);
    wen_n << pvapre.vel[1] / (rmrn[1] + pvapre.pose.translation()[2]), -pvapre.vel[0] / (rmrn[0] + pvapre.pose.translation()[2]),
        -pvapre.vel[1] * tan(pvapre.pose.translation()[0]) / (rmrn[1] + pvapre.pose.translation()[2]);
    double gravity = Earth::gravity(pvapre.pose.translation());

    // 旋转效应和双子样划桨效应
    // rotational and sculling motion
    temp1 = imucur.dtheta.cross(imucur.dvel) / 2;
    temp2 = imupre.dtheta.cross(imucur.dvel) / 12;
    temp3 = imupre.dvel.cross(imucur.dtheta) / 12;

    // b系比力积分项
    // velocity increment due to the specific force
    d_vfb = imucur.dvel + temp1 + temp2 + temp3;

    // 比力积分项投影到n系
    // velocity increment dut to the specfic force projected to the n-frame
    temp1 = (wie_n + wen_n) * imucur.dt / 2;
    cnn   = I33 - Rotation::skewSymmetric(temp1);
    d_vfn = cnn * pvapre.pose.rotationMatrix() * d_vfb;

    // 计算重力/哥式积分项
    // velocity increment due to the gravity and Coriolis force
    gl << 0, 0, gravity;
    d_vgn = (gl - (2 * wie_n + wen_n).cross(pvapre.vel)) * imucur.dt;

    // 得到中间时刻速度
    // velocity at k-1/2
    midvel = pvapre.vel + (d_vfn + d_vgn) / 2;

    // 外推得到中间时刻位置
    // position extrapolation to k-1/2
    qnn = Rotation::rotvec2quaternion(temp1);
    temp2 << 0, 0, -WGS84_WIE * imucur.dt / 2;
    qee = Rotation::rotvec2quaternion(temp2);
    qne = Earth::qne(pvapre.pose.translation());
    qne = qee * qne * qnn;
    midpos[2] = pvapre.pose.translation()[2] - midvel[2] * imucur.dt / 2;
    midpos    = Earth::blh(qne, midpos[2]);

    // 重新计算中间时刻的rmrn, wie_e, wen_n
    // recompute rmrn, wie_n, and wen_n at k-1/2
    rmrn = Earth::meridianPrimeVerticalRadius(midpos[0]);
    wie_n << WGS84_WIE * cos(midpos[0]), 0, -WGS84_WIE * sin(midpos[0]);
    wen_n << midvel[1] / (rmrn[1] + midpos[2]), -midvel[0] / (rmrn[0] + midpos[2]),
        -midvel[1] * tan(midpos[0]) / (rmrn[1] + midpos[2]);

    // 重新计算n系下平均比力积分项
    // recompute d_vfn
    temp3 = (wie_n + wen_n) * imucur.dt / 2;
    cnn   = I33 - Rotation::skewSymmetric(temp3);
    d_vfn = cnn * pvapre.pose.rotationMatrix() * d_vfb;

    // 重新计算重力、哥式积分项
    // recompute d_vgn
    gl << 0, 0, Earth::gravity(midpos);
    d_vgn = (gl - (2 * wie_n + wen_n).cross(midvel)) * imucur.dt;

    // 速度更新完成
    // velocity update finish
    pvacur.vel = pvapre.vel + d_vfn + d_vgn;
}

void posUpdate(const BodyState &pvapre, BodyState &pvacur, const IMU &imupre, const IMU &imucur) {

    Eigen::Vector3d temp1, temp2, midvel, midpos;
    Eigen::Quaterniond qne, qee, qnn;

    // 重新计算中间时刻的速度和位置
    // recompute velocity and position at k-1/2
    midvel = (pvacur.vel + pvapre.vel) / 2;
    midpos = pvapre.pose.translation() + Earth::DRi(pvapre.pose.translation()) * midvel * imucur.dt / 2;

    // 重新计算中间时刻地理参数
    // recompute rmrn, wie_n, wen_n at k-1/2
    Eigen::Vector2d rmrn;
    Eigen::Vector3d wie_n, wen_n;
    rmrn = Earth::meridianPrimeVerticalRadius(midpos[0]);
    wie_n << WGS84_WIE * cos(midpos[0]), 0, -WGS84_WIE * sin(midpos[0]);
    wen_n << midvel[1] / (rmrn[1] + midpos[2]), -midvel[0] / (rmrn[0] + midpos[2]),
        -midvel[1] * tan(midpos[0]) / (rmrn[1] + midpos[2]);

    // 重新计算 k时刻到k-1时刻 n系旋转矢量
    // recompute n-frame rotation vector (n(k) with respect to n(k-1)-frame)
    temp1 = (wie_n + wen_n) * imucur.dt;
    qnn   = Rotation::rotvec2quaternion(temp1);
    // e系转动等效旋转矢量 (k-1时刻k时刻，所以取负号)
    // e-frame rotation vector (e(k-1) with respect to e(k)-frame)
    temp2 << 0, 0, -WGS84_WIE * imucur.dt;
    qee = Rotation::rotvec2quaternion(temp2);

    // 位置更新完成
    // position update finish
    qne           = Earth::qne(pvapre.pose.translation());
    qne           = qee * qne * qnn;
    pvacur.pose.translation()[2] = pvapre.pose.translation()[2] - midvel[2] * imucur.dt;
    pvacur.pose.translation()    = Earth::blh(qne, pvacur.pose.translation()[2]);
}

void attUpdate(const BodyState &pvapre, BodyState &pvacur, const IMU &imupre, const IMU &imucur) {

    Eigen::Quaterniond qne_pre, qne_cur, qne_mid, qnn, qbb;
    Eigen::Vector3d temp1, midpos, midvel;

    // 重新计算中间时刻的速度和位置
    // recompute velocity and position at k-1/2
    midvel = (pvapre.vel + pvacur.vel) / 2;
    qne_pre   = Earth::qne(pvapre.pose.translation());
    qne_cur   = Earth::qne(pvacur.pose.translation());
    temp1     = Rotation::quaternion2vector(qne_cur.inverse() * qne_pre);
    qne_mid   = qne_pre * Rotation::rotvec2quaternion(temp1 / 2).inverse();
    midpos[2] = (pvacur.pose.translation()[2] + pvapre.pose.translation()[2]) / 2;
    midpos    = Earth::blh(qne_mid, midpos[2]);

    // 重新计算中间时刻地理参数
    // recompute rmrn, wie_n, wen_n at k-1/2
    Eigen::Vector2d rmrn;
    Eigen::Vector3d wie_n, wen_n;
    rmrn = Earth::meridianPrimeVerticalRadius(midpos[0]);
    wie_n << WGS84_WIE * cos(midpos[0]), 0, -WGS84_WIE * sin(midpos[0]);
    wen_n << midvel[1] / (rmrn[1] + midpos[2]), -midvel[0] / (rmrn[0] + midpos[2]),
        -midvel[1] * tan(midpos[0]) / (rmrn[1] + midpos[2]);

    // 计算n系的旋转四元数 k-1时刻到k时刻变换
    // n-frame rotation vector (n(k-1) with respect to n(k)-frame)
    temp1 = -(wie_n + wen_n) * imucur.dt;
    qnn   = Rotation::rotvec2quaternion(temp1);

    // 计算b系旋转四元数 补偿二阶圆锥误差
    // b-frame rotation vector (b(k) with respect to b(k-1)-frame)
    // compensate the second-order coning correction term.
    temp1 = imucur.dtheta + imupre.dtheta.cross(imucur.dtheta) / 12;
    qbb   = Rotation::rotvec2quaternion(temp1);

    // 姿态更新完成
    // attitude update finish
    // pvacur.pose.unit_quaternion()   = qnn * pvapre.pose.unit_quaternion() * qbb;
    pvacur.pose.rotationMatrix()   = Rotation::quaternion2matrix(pvacur.pose.unit_quaternion());
    // pvacur.pose.eulerAngles(2, 1, 0) = Rotation::matrix2euler(pvacur.pose.rotationMatrix());
}

void imuInterpolate(const IMU &imu1, IMU &imu2, const double timestamp, IMU &midimu) {

        if (imu1.time > timestamp || imu2.time < timestamp) {
            return;
        }

        double lamda = (timestamp - imu1.time) / (imu2.time - imu1.time);

        midimu.time   = timestamp;
        midimu.dtheta = imu2.dtheta * lamda;
        midimu.dvel   = imu2.dvel * lamda;
        midimu.dt     = timestamp - imu1.time;

        imu2.dtheta = imu2.dtheta - midimu.dtheta;
        imu2.dvel   = imu2.dvel - midimu.dvel;
        imu2.dt     = imu2.dt - midimu.dt;
    }

void imuCompensate(IMU &imu, ImuError &imuerror) {

  // compensate the imu bias
  imu.dtheta -= imuerror.gyrbias * imu.dt;
  imu.dvel -= imuerror.accbias * imu.dt;

  
  // imu.dtheta = imu.dtheta.cwiseProduct(gyrscale.cwiseInverse());
  // imu.dvel   = imu.dvel.cwiseProduct(accscale.cwiseInverse());
}




} // namespace lio_ekf