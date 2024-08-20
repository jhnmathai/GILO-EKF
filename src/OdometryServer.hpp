// This file was heavily inspired by KISS-ICP, for this reason we report here
// the original LICENSE
// MIT License
//
// Copyright (c) 2022 Ignacio Vizzo, Tiziano Guadagnino, Benedikt Mersch, Cyrill
// Stachniss.
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
#pragma once

#include <deque>
#include <fstream>
#include <mutex>
#include <nav_msgs/Path.h>
#include <ros/ros.h>
#include <sensor_msgs/Imu.h>
#include <sensor_msgs/PointCloud2.h>
#include <tf2_ros/transform_broadcaster.h>

#include "kiss_icp/pipeline/KissICP.hpp"
#include "lio_ekf.hpp"
#include "lio_types.hpp"

namespace lio_ekf {

class OdometryServer {
public:
  /// OdometryServer constructor
  OdometryServer(const ros::NodeHandle &nh, const ros::NodeHandle &pnh);

  // buffer imu data
  void imu_cbk(const sensor_msgs::Imu::ConstPtr &msg_in);
  void lidar_cbk(const sensor_msgs::PointCloud2ConstPtr &msg);

  void writeResults(std::ofstream &odo);

  void publishMsgs();

  std::deque<std_msgs::Header> lidar_header_buffer_;

private:
  lio_ekf::LIOPara lio_para_; // parameters for lidar-inertial fusion
  lio_ekf::LIOEKF lio_ekf_;   // lidar-inertial processor

  bool data_synced_ = false;

  // output dir
  std::string outputdir;
  std::ofstream odomRes_, odomRes_tum_,  odomRes_blh_;

  /// Ros node stuff
  ros::NodeHandle nh_;
  ros::NodeHandle pnh_;
  int queue_size_{1};
  std::string lid_topic, imu_topic;

  std::mutex mtx_buffer_;

  double last_timestamp_imu_, last_timestamp_lidar_;

  std::deque<std::vector<Eigen::Vector3d>> lidar_buffer_;
  std::deque<lio_ekf::IMU> imu_buffer_;
  std::deque<double> lidar_time_buffer_;
  std::deque<std::vector<double>> points_per_scan_time_buffer_;

  /// Tools for broadcasting TFs.
  tf2_ros::TransformBroadcaster tf_broadcaster_;

  /// Data subscribers.
  ros::Subscriber pointcloud_sub_;
  ros::Subscriber imu_sub_;

  /// Data publishers.
  ros::Publisher odom_publisher_;
  ros::Publisher traj_publisher_;
  nav_msgs::Path path_msg_;
  ros::Publisher frame_publisher_;
  ros::Publisher kpoints_publisher_;
  ros::Publisher map_publisher_;

  /// Global/map coordinate frame.
  std::string odom_frame_{"odom"};
  std::string pointcloud_frame_{"pointcloud_frame"};
};

} // namespace lio_ekf
