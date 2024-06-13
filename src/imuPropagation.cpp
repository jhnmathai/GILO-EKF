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

/**
 * @brief INS mechanization.
 *
 * @param pvapre 
 * @param pvacur 
 * @param imupre 
 * @param imucur 
 */
void insMechanization(const BodyState &pvapre, BodyState &pvacur,
                      const IMU &imupre, const IMU &imucur) {

  // perform velocity update, position updata and attitude update in sequence,
  // irreversible order

  // Variables:
  // t -> position
  // v -> velocity
  // a -> acceleration
  // \omega -> angular velocity
  // R -> orientation in the global frame
  // s -> integration interval (time step)
  // g -> gravity vector in global frame
  // Exp -> exponential mapping to the SO(3) group
  // b_g -> gyroscope bias
  // b_a -> accelerometer bias

  // Position update eq.
  // t_k = t_{k-1} + 1/2 (v_{k-1} + v_k) s_k

  // Velocity update eq.
  // v_k = v_{k-1} + g s_k + R_{k-1} [a_k s_k + 1/2 (w_k \times a_k) s_k^2 + ...
  // ... + 1/12 (w_{k-1} \times a_k - w_k \times a_{k-1}) s_k^2]

  // Orientation update eq.
  // R_k = R_{k-1} Exp (w_k s_k + 1/12 (w_{k-1} \times w_k) s_k^2)

  // Bias update eqs.  (Constant)
  // b_{g, k} = b_{g, k-1}
  // b_{a, k} = b_{a, k-1}


  Eigen::Vector3d d_vfb, d_vfn, d_vgn, gl;
  Eigen::Vector3d temp1, temp2, temp3;
  Eigen::Vector3d imucur_dvel, imucur_dtheta, imupre_dvel, imupre_dtheta;

  // This is a way of simplifying the velocity update equation.
  imucur_dvel = imucur.linear_acceleration * imucur.dt;
  imucur_dtheta = imucur.angular_velocity * imucur.dt;
  imupre_dvel = imupre.linear_acceleration * imupre.dt;
  imupre_dtheta = imupre.angular_velocity * imupre.dt;

  // rotational and sculling motion
  temp1 = imucur_dtheta.cross(imucur_dvel) / 2;   // 1/2 (\theta_k \times v_k)
  temp2 = imupre_dtheta.cross(imucur_dvel) / 12;  // 1/12 (\theta_k \times v_k)
  temp3 = imupre_dvel.cross(imucur_dtheta) / 12;  // 1/12 (v_k \times \theta_k)

  // velocity increment due to the specific force
  // v_k + 1/2 (\theta_k \times v_k) + 1/12 (\theta_k \times v_k) + 1/12 (v_k \times \theta_k)
  d_vfb = imucur_dvel + temp1 + temp2 + temp3;

  // velocity increment dut to the specfic force projected to the n-frame
  d_vfn = pvapre.pose.rotationMatrix() * d_vfb;

  // velocity increment due to the gravity and Coriolis force
  gl << 0, 0, NormG;
  d_vgn = gl * imucur.dt;

  // velocity update finish -> v_{k-1} + temp Term + gravity vector in global frame
  pvacur.vel = pvapre.vel + d_vfn + d_vgn;

  Eigen::Vector3d midvel;

  // recompute velocity and position at k-1/2
  midvel = (pvacur.vel + pvapre.vel) / 2;  // 1/2 (v_k + v_{k-1})
  pvacur.pose.translation() += midvel * imucur.dt;  // Updates the position in the BodyState object
  Eigen::Vector3d rot_bframe;

  // b-frame rotation vector (b(k) with respect to b(k-1)-frame)
  // compensate the second-order coning correction term.
  // rot_brame = \theta_k + 1/12 (\theta_{k-1} \times \theta_k)
  rot_bframe = imucur_dtheta + imupre_dtheta.cross(imucur_dtheta) / 12;

  // R_k = R_{k-1} * Exp (rot_bframe)
  pvacur.pose.so3() = pvapre.pose.so3() * Sophus::SO3d::exp(rot_bframe);
}


/**
 * @brief Simple IMU measurement interpolation.
 * 
 * @param imu1 
 * @param imu2 
 * @param timestamp 
 * @param midimu 
 */
void imuInterpolate(const IMU &imu1, IMU &imu2, const double timestamp,
                    IMU &midimu) {

  if (imu1.timestamp > timestamp || imu2.timestamp < timestamp) {
    return;
  }

  // \lambda = (t - t_{imu_1}) / (t_{imu_2} - t_{imu_1})
  double lambda =
      (timestamp - imu1.timestamp) / (imu2.timestamp - imu1.timestamp);

  // t_{imu_{1/2}} = t
  midimu.timestamp = timestamp;

  // \Delta t_{imu_{1/2/}} = t - t_{imu_1}
  midimu.dt = timestamp - imu1.timestamp;
  // \omega_{imu_{1/2}} = \lambda * \omega_{imu_2} + (1 - \lambda) * \omega_{imu_1}
  midimu.angular_velocity =
      lambda * imu2.angular_velocity + (1 - lambda) * imu1.angular_velocity;

  // a_{imu_{1/2}} = \lambda * a_{imu_2} + (1 - \lambda) * a_{imu_1}
  midimu.linear_acceleration = lambda * imu2.linear_acceleration +
                               (1 - lambda) * imu1.linear_acceleration;

  // \Delta t_{imu_2} = \Delta t_{imu_2} - \Delta t_{imu_{1/2}}
  imu2.dt = imu2.dt - midimu.dt;
}

/**
 * @brief Adds the Bias correction to the IMU measurements.
 * 
 * @param imu 
 * @param imuerror 
 */
void imuCompensate(IMU &imu, ImuError &imuerror) {

  // compensate the imu bias error
  // \omega = \omega_{true} - b_g
  imu.angular_velocity -= imuerror.gyrbias;
  // f = f_{true} - b_a
  imu.linear_acceleration -= imuerror.accbias;
}

} // namespace lio_ekf