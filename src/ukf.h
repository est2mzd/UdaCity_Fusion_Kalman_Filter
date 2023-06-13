#ifndef UKF_H
#define UKF_H

#include "Eigen/Dense"
#include "measurement_package.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF
{
public:
  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * InitializeState
   * @param meas_package The latest measurement data of either radar or laser
   */
  void InitializeState(MeasurementPackage meas_package);

  /**
   * Step-1: Prediction Predicts sigma points, the state, and the state covariance matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * @brief Step-1-1: GenerateAugmentedSigmaPoints
   *       (Ref = 18. Augmentation Assignment 2)
   * @param[out] x_sig_aug Augmented Sigma Points 
   *                      [pos1, pos2, vel_abs, yaw_angle, yaw_rate, noise_a, noise_yawdd]
   * @param[out] p_aug     Augmanetd Covariance Matrix
   */
  void GenerateAugmentedSigmaPoints(MatrixXd& x_sig_aug, MatrixXd& p_aug);

  /**
   * @brief Step-1-2: PredictSigmaPoints
   *        (Ref = 21. Sigma Point Prediction Assignment 2)
   * @param[in]  x_sig_aug Augmented Sigma Points
   *                      [pos1, pos2, vel_abs, yaw_angle, yaw_rate, noise_a, noise_yawdd]
   * @param[out] void Xsig_pred will be modified
   */
  void PredictSigmaPoints(const MatrixXd& x_sig_aug, double delta_t);


  /**
   * @brief Step-1-3: CalculatePredictedMeanAndCovariance matrix
   *       (Ref = 24. Predicted Mean and Covariance Assignment 2)
   * @param[in]  void Xsig_pred_ and weights_ will be used
   * @param[out] void x_ and P_ will be modified
   */
  void CalculatePredictedMeanAndCovariance();

  /**
   * Step-2: Updates the state and the state covariance matrix
   * @param meas_package The measurement at k+1
   */
  void UpdateCommon(MeasurementPackage meas_package);

  /**
   * Step-2L: Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Step-2R: Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);


  /**
   * @brief Step-2-1: TransformSigmaPointsFromPredictionToMeasurementSpace for Lidar and Radar
   *        Ref = 27. Predict Radar Measurement Assignment 2
   * @param[in]  meas_package
   * @param[out] z_sig_meas Sigma points in measurement space
   */
  MatrixXd TransformSigmaPointsFromPredictionToMeasurementSpace(MeasurementPackage::SensorType sensor_type);


  /**
   * @brief Step-2-2: CalculateMeanPredictedMeasurement for Lidar and Radar
   *        Ref = 27. Predict Radar Measurement Assignment 2
   * @param[in]  z_sig_meas      Sigma points in measurement space
   * @param[out] z_sig_meas_mean Mean Sigma points in measurement space
   */
  VectorXd CalculateMeanPredictedMeasurement(MatrixXd z_sig_meas);


  /**
   * @brief Step-2-3: CalculateInnovationCovarianceMatrix for Lidar and Radar
   *        Ref = 27. Predict Radar Measurement Assignment 2
   * @param[in]  z_sig_meas      Sigma points in measurement space
   * @param[in]  z_sig_meas_mean Mean Sigma points in measurement space
   * @param[out] S Innovation Covariance Matrix
   */
  MatrixXd CalculateInnovationCovarianceMatrix(MatrixXd z_sig_meas, VectorXd z_sig_meas_mean, MeasurementPackage::SensorType sensor_type);


  /**
   * @brief Step-2-4: CalculateCrossCorrelationMatrix for Lidar and Radar
   *        Ref = 30. UKF Update Assignment 2
   * @param[in] meas_package
   * @param[out] Tc Cross Correlation Matrix
   */
  MatrixXd CalculateCrossCorrelationMatrix(MatrixXd z_sig_meas, VectorXd z_sig_meas_mean, MeasurementPackage::SensorType sensor_type);

  /**
   * @brief Step-2-5: CalculateKalmanGain for Lidar and Radar
   *        Ref = 30. UKF Update Assignment 2
   * @param[in]  Tc Cross Correlation Matrix
   * @param[in]  S Predicted Measurement Covariance
   * @param[out] K Kalman Gain
   */
  MatrixXd CalculateKalmanGain(MatrixXd Tc, MatrixXd S);

  /**
   * @brief Step-2-6: UpdateStateMeanAndCovarianceMatrix for Lidar and Radar
   * @param[in ] meas_package sensing data
   * @param[in ] z_sig_meas_mean mean predicted measurement
   * @param[out] void  x_ and P_ will be updated
   */
  void UpdateStateMeanAndCovarianceMatrix(MeasurementPackage meas_package, VectorXd z_sig_meas_mean,
                                                  MatrixXd K, MatrixXd S);


  /**
   * ModifyAngle : -pi to +pi
   * @param[in] angle_rad angle
   */
  double ModifyAngle(double angle_rad);


  /**
   * GetStateDegreeOfFreedom
   * @param[in] n_z degree of freedom : Lidar or Radar
   */
  int GetStateDegreeOfFreedom(MeasurementPackage::SensorType sensor_type);


  // initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  // if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  // if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  // state vector: [pos1, pos2, vel_abs, yaw_angle, yaw_rate] in SI units and rad
  Eigen::VectorXd x_;

  // state covariance matrix
  Eigen::MatrixXd P_;

  // predicted sigma points matrix
  Eigen::MatrixXd Xsig_pred_;

  // time when the state is true, in us
  long long time_us_;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  // Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  // Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  // Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  // Radar measurement noise standard deviation radius in m
  double std_radr_;

  // Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  // Radar measurement noise standard deviation radius change in m/s
  double std_radrd_;

  // Weights of sigma points
  Eigen::VectorXd weights_;

  // State dimension
  int n_x_;

  // Augmented state dimension
  int n_aug_;

  // Sigma point spreading parameter
  double lambda_;
};

#endif // UKF_H