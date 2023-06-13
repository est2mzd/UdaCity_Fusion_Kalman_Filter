#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
using Eigen::MatrixXd;
using Eigen::VectorXd;

#include <time.h>

/**
 * Initializes Unscented Kalman filter
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

  /**  -------------------------------------------
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   * -------------------------------------------
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 2;
  

  //time_us_        = 0.0;
  is_initialized_ = false;
  n_x_            = 5; // State Dimension
  n_aug_          = 7; // Augumented State Dimension  
  lambda_         = 3 - n_aug_;
  
  weights_    = VectorXd(2 * n_aug_ + 1);
  weights_.fill(0.5 / (lambda_ + n_aug_));
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  
  x_.fill(0.0);
  P_.fill(0.0);

  Xsig_pred_ = MatrixXd(n_x_, 1 + n_aug_*2);
  Xsig_pred_.fill(0.0)  ;
}

//===========================================================//

UKF::~UKF() {}

//===========================================================//

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   * [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
   */

  //-------- Step-0: Initialize State Variable X and P --------//
  if(!is_initialized_)
  {
    InitializeState(meas_package);
    return;
  }

  //-------- Step-1: Predict with a CTRV (=Constant Turn Rate Velocity) Model --------//
  double delta_t = (meas_package.timestamp_ - time_us_) / 1.0e+6; // sec
  Prediction(delta_t);

  //-------- Step-2: Measurement Update --------//
  UpdateLidar(meas_package);
  UpdateRadar(meas_package);
  //std::cout << "End : ProcessMeasurement()" << std::endl;
}


//===========================================================//

void UKF::InitializeState(MeasurementPackage meas_package) 
{
  //std::cout << "---------- Start: InitializeState" << std::endl;
  //---- Lidar
  if(meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER)
  {
    // raw_measurements_ = marker.x, marker.y;
    // x_ = [pos1 pos2 vel_abs yaw_angle yaw_rate]
    x_.fill(0.0);
    x_(0) = meas_package.raw_measurements_(0); // px
    x_(1) = meas_package.raw_measurements_(1); // py

    P_.fill(0.0);
    P_(0,0) = 1.0;//std_laspx_*std_laspx_;
    P_(1,1) = 1.0;//std_laspy_*std_laspy_;
    P_(2,2) = P_(3,3) = P_(4,4) = 0.5;//1.0;
  }
  
  //---- Radar
  if(meas_package.sensor_type_ == MeasurementPackage::SensorType::RADAR)
  {
    // raw_measurements_ = marker.rho, marker.phi, marker.rho_dot;
    // x_ = [pos1 pos2 vel_abs yaw_angle yaw_rate]
    
    double rho  = meas_package.raw_measurements_(0);
    double phi  = meas_package.raw_measurements_(1);
    double rhod = meas_package.raw_measurements_(2);

    x_.fill(0.0);
    x_(0)  = rho * cos(phi);
    x_(1)  = rho * sin(phi);
    x_(2) = fabs(meas_package.raw_measurements_(2)); // vel_abs
    x_(3) = 0.0; // yaw
    x_(4) = 0.0; // yawd

    P_.fill(0.0);
    P_(0,0) = 1.0;//std_radr_ *std_radr_;
    P_(1,1) = 1.0;//std_radrd_*std_radrd_;
    P_(2,2) = 0.5;//std_radrd_*std_radrd_;
    P_(3,3) = 0.5;//std_radphi_;
    P_(4,4) = 0.5;//std_radphi_;
  }
  
  std::cout << "----------- Initilization ------------" << std::endl;
  std::cout << "x_" << std::endl << x_ << std::endl;
  std::cout << "P_" << std::endl << P_ << std::endl;

  is_initialized_ = true;
  time_us_ = meas_package.timestamp_;        // us

  //std::cout << "---------- End: InitializeState" << std::endl;
}


//===========================================================//

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

  // Step-1-1. Generate Sigma Points 
  MatrixXd x_sig_aug = MatrixXd(n_aug_, 1 + 2*n_aug_);
  MatrixXd p_aug     = MatrixXd(n_aug_, n_aug_);
  GenerateAugmentedSigmaPoints(x_sig_aug, p_aug);

  // Step-1-2. Predict Sigma Points >>> Update Xsig_pred
  PredictSigmaPoints(x_sig_aug, delta_t);

  // Step-1-3. Calculate Predicted Mean and Covariance           
  CalculatePredictedMeanAndCovariance();

}


//===========================================================//

void UKF::GenerateAugmentedSigmaPoints(MatrixXd& x_sig_aug, MatrixXd& p_aug)
{
  //std::cout << "---------- Start: Step-1-1: GenerateAugmentedSigmaPoints" << std::endl;
  // Create augmented mean vector
  VectorXd x_aug_mean = VectorXd(n_aug_);
  x_aug_mean.head(n_x_) = x_; // [pos1, pos2, vel_abs, yaw_angle, yaw_rate]
  x_aug_mean(n_x_)      =  0; // longitudinal acceleration noise
  x_aug_mean(n_x_+1)    =  0; // yaw acceleration noise

  // Create augmented covariance matrix
  p_aug.fill(0.0);
  p_aug.topLeftCorner(n_x_, n_x_) = P_;
  p_aug(n_x_,   n_x_)   = std_a_*std_a_;
  p_aug(n_x_+1, n_x_+1) = std_yawdd_*std_yawdd_;

  // Create square root matrix
  MatrixXd L = p_aug.llt().matrixL();
  double   coef_of_L = sqrt(lambda_ + n_aug_);

  // Create augmented sigma points
  x_sig_aug.col(0) = x_aug_mean;
  for(int col=0; col < n_aug_; ++col)
  {
    x_sig_aug.col(col+1)        = x_aug_mean + coef_of_L * L.col(col);
    x_sig_aug.col(col+1+n_aug_) = x_aug_mean - coef_of_L * L.col(col);
  }
  ////std::cout << "---------- End: Step-1-1: GenerateAugmentedSigmaPoints" << std::endl;
}

//===========================================================//

void UKF::PredictSigmaPoints(const MatrixXd& x_sig_aug, double delta_t)
{
  ////std::cout << "---------- Start: Step-1-2: PredictSigmaPoints" << std::endl;
  for(int col=0; col < 1 + 2*n_aug_; ++col)
  {
    // Extract variables
    double px       = x_sig_aug(0, col);
    double py       = x_sig_aug(1, col);
    double v        = x_sig_aug(2, col);
    double yaw      = x_sig_aug(3, col);
    double yawd     = x_sig_aug(4, col);
    double nu_a     = x_sig_aug(5, col);
    double nu_yawdd = x_sig_aug(6, col);

    // Predict sigma points : x & y, avoid division by zero
    if(fabs(yawd) > 1.0e-3)
    {
      double theta      = yaw + yawd * delta_t;
      Xsig_pred_(0,col) = px + v / yawd * ( sin(theta) - sin(yaw)); // x
      Xsig_pred_(1,col) = py + v / yawd * (-cos(theta) + cos(yaw)); // y
    }
    else
    {
      Xsig_pred_(0,col) = px + v * cos(yaw) * delta_t; // x
      Xsig_pred_(1,col) = py + v * sin(yaw) * delta_t; // y
    }

    // Predict sigma points : others
    Xsig_pred_(2,col) = v;    // v
    Xsig_pred_(3,col) = yaw + yawd * delta_t; // yaw
    Xsig_pred_(4,col) = yawd; // yawd 

    // Add noise
    Xsig_pred_(0,col) += 0.5 * delta_t * delta_t * cos(yaw) * nu_a; // x
    Xsig_pred_(1,col) += 0.5 * delta_t * delta_t * sin(yaw) * nu_a; // y
    Xsig_pred_(2,col) += delta_t * nu_a; // v
    Xsig_pred_(3,col) += 0.5 * delta_t * delta_t * nu_yawdd; // yaw
    Xsig_pred_(4,col) += delta_t * nu_yawdd; // yawd
  }
  //std::cout << "---------- End: Step-1-2: PredictSigmaPoints" << std::endl;
}

//===========================================================//

void UKF::CalculatePredictedMeanAndCovariance()
{
  //std::cout << "---------- Start: Step-1-3: CalculatePredictedMeanAndCovariance" << std::endl;
  // Initialization
  x_.fill(0.0);
  P_.fill(0.0);
  //
  for(int col = 0; col < Xsig_pred_.cols(); ++col)
  {
    x_ += weights_(col) * Xsig_pred_.col(col);
  }
  //
  for(int col = 0; col < Xsig_pred_.cols(); ++col)
  {
    VectorXd diff_x = Xsig_pred_.col(col) - x_;
    diff_x(3)       = ModifyAngle(diff_x(3));
    P_ += weights_(col) * diff_x * diff_x.transpose();
  }
  //std::cout << "---------- End: Step-1-3: CalculatePredictedMeanAndCovariance" << std::endl;
}

//===========================================================//

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  if(meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER && use_laser_)
  {
    //std::cout << "=================== LASER ==================== " << std::endl;
    UpdateCommon(meas_package);
  }
}


//===========================================================//

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  if(meas_package.sensor_type_ == MeasurementPackage::SensorType::RADAR && use_radar_)
  {
    //std::cout << "=================== RADAR ==================== " << std::endl;
    UpdateCommon(meas_package);
  }
}



//===========================================================//

void UKF::UpdateCommon(MeasurementPackage meas_package)
{
  // Step-2-1 : Transform sigma points into measurement space
  //std::cout << " ---------- Step-2-1 ---------- SensorType= " << meas_package.sensor_type_ << std::endl;
  MatrixXd z_sig_meas = TransformSigmaPointsFromPredictionToMeasurementSpace(meas_package.sensor_type_);
  //std::cout << "z_sig_meas = " << std::endl << z_sig_meas << std::endl;

  // Step-2-2: Calculate mean predicted measurement
  //std::cout << " ---------- Step-2-2 ---------- SensorType= " << meas_package.sensor_type_ << std::endl;
  VectorXd z_sig_meas_mean = CalculateMeanPredictedMeasurement(z_sig_meas);
  //std::cout << "z_sig_meas_mean = " << std::endl << z_sig_meas_mean << std::endl;

  // Step-2-3 : Calculate innovation covariance matrix S
  //std::cout << " ---------- Step-2-3 ---------- SensorType= " << meas_package.sensor_type_ << std::endl;
  MatrixXd S = CalculateInnovationCovarianceMatrix(z_sig_meas, z_sig_meas_mean, meas_package.sensor_type_);
  //std::cout << "S = " << std::endl << S << std::endl;

  // Step-2-4: Calculate cross correlation matrix
  //std::cout << " ---------- Step-2-4 ---------- SensorType= " << meas_package.sensor_type_ << std::endl;
  MatrixXd Tc = CalculateCrossCorrelationMatrix(z_sig_meas, z_sig_meas_mean, meas_package.sensor_type_);
  //std::cout << "Tc = " << std::endl << Tc << std::endl;

  // Step-2-5: Calculate Kalman gain K
  //std::cout << " ---------- Step-2-5 ---------- SensorType= " << meas_package.sensor_type_ << std::endl;
  MatrixXd K = CalculateKalmanGain(Tc,S);
  //std::cout << "K = " << std::endl << K << std::endl;

  // Step-2-6: Update state mean and covariance matrix
  //std::cout << " ---------- Step-2-6 ---------- SensorType= " << meas_package.sensor_type_ << std::endl;
  UpdateStateMeanAndCovarianceMatrix(meas_package, z_sig_meas_mean, K, S);

  //std::cout << "x_ = " << std::endl << x_ << std::endl;
  //std::cout << "P_ = " << std::endl << P_ << std::endl;
}

//===========================================================//

MatrixXd UKF::TransformSigmaPointsFromPredictionToMeasurementSpace(MeasurementPackage::SensorType sensor_type)
{
  int n_z = GetStateDegreeOfFreedom(sensor_type); // Get DOF
  MatrixXd z_sig_meas(n_z, 1+2*n_aug_);
  z_sig_meas.fill(0.0);

  if(sensor_type == MeasurementPackage::SensorType::LASER)
  {
    for(int col=0; col < z_sig_meas.cols(); ++col)
    {
      z_sig_meas(0, col) = Xsig_pred_(0, col); // x
      z_sig_meas(1, col) = Xsig_pred_(1, col); // y
    }
  }
  else if(sensor_type == MeasurementPackage::SensorType::RADAR)
  {
    for(int col=0; col < z_sig_meas.cols(); ++col)
    {
      double px   = Xsig_pred_(0, col);
      double py   = Xsig_pred_(1, col);
      double v    = Xsig_pred_(2, col);
      double yaw  = Xsig_pred_(3, col);
      double rho  = sqrt(px*px + py*py); // range

      z_sig_meas(0, col) = rho;
      z_sig_meas(1, col) = atan2(py, px); // yaw

      if(rho > 0.001)
      {
        z_sig_meas(2, col) = (px*cos(yaw)*v + py*sin(yaw)*v) / rho; // v
      }
      else
      {
        z_sig_meas(2, col) = (px*cos(yaw)*v + py*sin(yaw)*v) / 0.001; // v
      }
    }
  }

  return z_sig_meas;
}

//===========================================================//

VectorXd UKF::CalculateMeanPredictedMeasurement(MatrixXd z_sig_meas)
{
  int n_z = z_sig_meas.rows();
  VectorXd z_sig_meas_mean(n_z);
  z_sig_meas_mean.fill(0);

  for(int col=0; col < z_sig_meas.cols(); ++col)
  {
    z_sig_meas_mean += weights_(col) * z_sig_meas.col(col);
  }

  return z_sig_meas_mean;
}

//===========================================================//

MatrixXd UKF::CalculateInnovationCovarianceMatrix(MatrixXd z_sig_meas, VectorXd z_sig_meas_mean, MeasurementPackage::SensorType sensor_type)
{
  //std::cout << "   A " << std::endl;
  int n_z = GetStateDegreeOfFreedom(sensor_type);
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);

  //std::cout << "   B " << std::endl;
  for(int col=0; col < z_sig_meas.cols(); ++col)
  {
    VectorXd diff_z = z_sig_meas.col(col) - z_sig_meas_mean;

    if(sensor_type == MeasurementPackage::SensorType::RADAR) // RADAR
    {
      diff_z(1) = ModifyAngle(diff_z(1));
    }

    S += weights_(col) * diff_z * diff_z.transpose();
  }

  //std::cout << "   C " << std::endl;
  MatrixXd R = MatrixXd(n_z, n_z);
  R.fill(0.0);

  if(sensor_type == MeasurementPackage::SensorType::LASER) // Lidar
  {
    R(0,0) = std_laspx_*std_laspx_;
    R(1,1) = std_laspy_*std_laspy_;
  }
  else // Radar
  {
    R(0,0) = std_radr_*std_radr_;
    R(1,1) = std_radphi_*std_radphi_;
    R(2,2) = std_radrd_*std_radrd_;
  }

  //std::cout << "   D " << std::endl;
  S += R;
  return S;
}


//===========================================================//

MatrixXd UKF::CalculateCrossCorrelationMatrix(MatrixXd z_sig_meas, VectorXd z_sig_meas_mean, 
                                                MeasurementPackage::SensorType sensor_type)
{
  int n_z = GetStateDegreeOfFreedom(sensor_type);
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  for(int col = 0; col < z_sig_meas.cols(); ++col)
  {
    VectorXd diff_x = Xsig_pred_.col(col) - x_;
    VectorXd diff_z = z_sig_meas.col(col) - z_sig_meas_mean;
    diff_x(3) = ModifyAngle(diff_x(3));

    if(sensor_type == MeasurementPackage::SensorType::RADAR)
    {
      diff_z(1) = ModifyAngle(diff_z(1));
    }

    // sum
    Tc += weights_(col) * diff_x * diff_z.transpose();
  }

  return Tc;
}

//===========================================================//

MatrixXd UKF::CalculateKalmanGain(MatrixXd Tc, MatrixXd S)
{
  return Tc * S.inverse();
}

//===========================================================//

void UKF::UpdateStateMeanAndCovarianceMatrix(MeasurementPackage meas_package, VectorXd z_sig_meas_mean,
                                                  MatrixXd K, MatrixXd S)
{
  int n_z = GetStateDegreeOfFreedom(meas_package.sensor_type_);
  VectorXd z_meas = VectorXd(n_z);

  for(int i=0; i < n_z; ++i)
  {
    z_meas(i) = meas_package.raw_measurements_(i);
  }

  //std::cout << "   C " << std::endl;
  VectorXd diff_z = z_meas - z_sig_meas_mean;

  //std::cout << "   D " << std::endl;
  if(meas_package.sensor_type_ == MeasurementPackage::SensorType::RADAR)
  {
    diff_z(1) = ModifyAngle(diff_z(1));
  }

  //std::cout << "   E " << std::endl;
  x_ += K * diff_z;
  //std::cout << "   F " << std::endl;
  P_ -= K * S * K.transpose();
  //std::cout << "   G " << std::endl;
}


//===========================================================//

double UKF::ModifyAngle(double angle_rad)
{
  if(angle_rad > M_PI)
  {
    angle_rad -= 2.0*M_PI;
  }
  //
  if(angle_rad < -M_PI)
  {
    angle_rad += 2.0*M_PI;
  }  
  //
  return angle_rad;
}

//===========================================================//

int UKF::GetStateDegreeOfFreedom(MeasurementPackage::SensorType sensor_type)
{
  int n_z = 2; // for Lidar
  
  if(sensor_type == MeasurementPackage::SensorType::RADAR)
  {
    n_z = 3;
  }

  return n_z;
}