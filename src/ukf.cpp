#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

#define ALMOST_ZERO 0.0001

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

namespace {
  double NormalizeAngle(double angle) {
    if (angle < -M_PI || angle > M_PI) {
      angle = fmod(angle, 2. * M_PI);
    }
    if (angle < -M_PI) angle += 2. * M_PI;
    if (angle > M_PI) angle -= 2. * M_PI;
    return angle;
  }
}
/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.0;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_ = false;
  time_us_ = 0;

  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;

  weights_ = VectorXd(2 * n_aug_ + 1);
  double denom = lambda_ + n_aug_;
  weights_(0) = lambda_ / denom;
  weights_.tail(2 * n_aug_) = VectorXd::Constant(2 * n_aug_, 1 / (2 * denom));

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_) {
    cout << "UKF: " << endl;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      auto rho = meas_package.raw_measurements_[0];
      auto phi = meas_package.raw_measurements_[1];
      auto px = cos(phi) * rho;
      auto py = sin(phi) * rho;
      x_ << px, py, 0, 0, 0;
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }

    P_.fill(0.0);
    P_(0, 0) = 0.5;
    P_(1, 1) = 0.5;
    P_(2, 2) = 10.5;
    P_(3, 3) = 10.5;
    P_(4, 4) = 10.2;

    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  if (dt < ALMOST_ZERO) {
    cout << "this measurement took place at the same time as the last one, skipping update..." << endl;
    cout << meas_package.timestamp_ << endl;
    cout <<  time_us_ << endl;
    return;
  }
  time_us_ = meas_package.timestamp_;

  // step 2: predict the ground_truth state at the next step
  Prediction(dt);
  // step 3: mapping the prediction to the measurement space happens here
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    cout << "in radar update..." << endl;
    UpdateRadar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    cout << "in lidar update..." << endl;
    UpdateLidar(meas_package);
  }

  // print the output
  cout << "x_: " << x_ << endl;
  cout << "P_: " << P_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  // step 1: generate sigma points with augmented x
  MatrixXd Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  VectorXd x_aug = VectorXd(n_aug_);
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  P_aug.fill(0.0);
  P_aug.topLeftCorner(5, 5) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;

  MatrixXd L = P_aug.llt().matrixL();

  Xsig_aug_.col(0) = x_aug;
  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug_.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug_.col(i + n_aug_ + 1) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }
  // step 2: predict sigma points for time step k + 1
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    double p_x = Xsig_aug_(0, i);
    double p_y = Xsig_aug_(1, i);
    double v = Xsig_aug_(2, i);
    double yaw = Xsig_aug_(3, i);
    double yawd = Xsig_aug_(4, i);
    double nu_a = Xsig_aug_(5, i);
    double nu_yawdd = Xsig_aug_(6, i);

    double px_pred, py_pred;

    if (fabs(yawd) > ALMOST_ZERO) {
      px_pred = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
      py_pred = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
    } else {
      px_pred = p_x + v * delta_t * cos(yaw);
      py_pred = p_y + v * delta_t * sin(yaw);
    }

    double v_pred = v;
    double yaw_pred = yaw + yawd * delta_t;
    double yawd_pred = yawd;

    px_pred += 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_pred += 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_pred += nu_a * delta_t;

    yaw_pred += 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_pred += nu_yawdd * delta_t;

    Xsig_pred_(0, i) = px_pred;
    Xsig_pred_(1, i) = py_pred;
    Xsig_pred_(2, i) = v_pred;
    Xsig_pred_(3, i) = yaw_pred;
    Xsig_pred_(4, i) = yawd_pred;
  }
  // step 3: compute mean and covariance matrix from sigma points
  // store them in x_ and P_
  x_ = Xsig_pred_ * weights_;
  P_.fill(0.0);
  MatrixXd Xsig_diff = Xsig_pred_.colwise() - x_;
  for (int i = 0; i < Xsig_pred_.cols(); i++) {
    Xsig_diff(3, i) = NormalizeAngle(Xsig_diff(3, i));
    P_ += weights_(i) * Xsig_diff.col(i) * Xsig_diff.col(i).transpose();
  }
  cout << "predicted state and covariance matrix: " << endl;
  cout << x_ << endl;
  cout << P_ << endl;
  cout << "prediction ends..." << endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  int n_z = 2;
  auto z = meas_package.raw_measurements_;
  MatrixXd H = MatrixXd(n_z, n_x_);
  H << 1, 0, 0, 0, 0,
      0, 1, 0, 0, 0;
  VectorXd y = z - H * x_;
  MatrixXd Ht = H.transpose();
//  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  VectorXd z_pred = VectorXd(n_z);
//  MatrixXd S = MatrixXd(n_z, n_z);
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_laspx_ * std_laspx_, 0,
          0, std_laspy_ * std_laspy_;
  MatrixXd S = H * P_ * Ht + R;
  MatrixXd K = P_ * Ht * S.inverse();
  MatrixXd I = MatrixXd::Identity(n_x_, n_x_);

  x_ += K * y;
  P_ = (I - K * H) * P_;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  // transform sigma points to radar measurement space
  int n_z = 3;
  auto z = meas_package.raw_measurements_;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  // z_pred is the predicted radar reading
  // S is its covariance matrix
  VectorXd z_pred = VectorXd(n_z);
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);

  for (int i = 0; i < Xsig_pred_.cols(); i++) {
    VectorXd z_sig_i = VectorXd(n_z);
    VectorXd x_sig_i = Xsig_pred_.col(i);
    double px = x_sig_i(0);
    double py = x_sig_i(1);
    double psi = x_sig_i(3);
    double v = x_sig_i(2);
    z_sig_i(0) = sqrt(px * px + py * py);
    if (px == 0.0) {
      cout << "px is zero, skipping this update..." << endl;
      return;
    }
    z_sig_i(1) = atan2(py, px);
    if (z_sig_i(0) < ALMOST_ZERO) {
      z_sig_i(2) = 0.0;
    } else {
      z_sig_i(2) = (px * cos(psi) * v + py * sin(psi) * v) / z_sig_i(0);
    }
    Zsig.col(i) = z_sig_i;
  }
  // calculate mean predicted measurement
  z_pred = Zsig * weights_;
  // calculate innovation covariance matrix S
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_radr_ * std_radr_, 0, 0,
      0, std_radphi_ * std_radphi_, 0,
      0, 0, std_radrd_ * std_radrd_;
  MatrixXd Zsig_diff = Zsig.colwise() - z_pred;
  for (int i = 0; i < Zsig.cols(); i++) {
    // angles in Zsig_diff are normalized
    Zsig_diff(1, i) = NormalizeAngle(Zsig_diff(1, i));
    S += weights_(i) * Zsig_diff.col(i) * Zsig_diff.col(i).transpose();
  }
  S += R;

  // update step
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  MatrixXd Xsig_diff = Xsig_pred_.colwise() - x_;
  for (int i = 0; i < Zsig.cols(); i++) {
    Zsig_diff(1, i) = NormalizeAngle(Zsig_diff(1, i));
    Xsig_diff(3, i) = NormalizeAngle(Xsig_diff(3, i));
    Tc += weights_(i) * Xsig_diff.col(i) * Zsig_diff.col(i).transpose();
  }

  VectorXd z_residual = z - z_pred;
  z_residual(1) = NormalizeAngle(z_residual(1));
  MatrixXd K = Tc * S.inverse();
  x_ += K * z_residual;
  P_ -= K * S * K.transpose();
}
