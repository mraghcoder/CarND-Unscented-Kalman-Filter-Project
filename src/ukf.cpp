#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

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
  x_[2] = 0.0;
  x_[3] = 0.0;
  x_[4] = 0.0;
  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_<<0.0225, 0, 0, 0, 0,
      0, 0.0225, 0, 0, 0,
      0, 0, 0.1, 0, 0,
      0, 0, 0, 0.1, 0,
      0, 0, 0, 0, 0.1;
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.07; //30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.97; //M_PI/3.0; //M_PI/4; // 30;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15; //0.3;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15; //0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;//0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03; //0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3; //0.3;

  ///* State dimension
  n_x_ = 5;

  ///* Augmented state dimension
  n_aug_ = n_x_ + 2;

  ///* Sigma point spreading parameter
  lambda_ = 3 - n_aug_ ;

  ///* predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);
  Xsig_pred_.fill(0.0);

  ///* Weights of sigma points
  weights_ = VectorXd(2*n_aug_+1);
  float w = 1.0/(2*(lambda_+n_aug_));
  weights_.fill(w);
  weights_(0) = lambda_/(lambda_+n_aug_);

  // cout<<"Weights"<<weights_<<endl;
  // cout<<"W sum"<<weights_.sum()<<endl;

  NIS_radar_ = 0.0;
  NIS_laser_ = 0.0;

  is_initialized_ = false;


  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
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
  if (!is_initialized_)
  {
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
      float rho = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];
      float rho_d = meas_package.raw_measurements_[2];
      x_[0] = rho*cos(phi);
      x_[1] = rho*sin(phi);
      //x_[2] = 2.0; //sqrt(rho_d*cos(phi)*rho_d*cos(phi)+rho_d*sin(phi)*rho_d*sin(phi));
      //x_[3] = 0.01;
      //x_[4] = 0.001;
      // P_<<1, 0, 0, 0, 0,
      //     0, 1, 0, 0, 0,
      //     0, 0, 0.5, 0, 0,
      //     0, 0, 0, 0.5, 0,
      //     0, 0, 0, 0, 1;
      //P_ = P_*0.8;
      cout<<"Radar Init"<<endl;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER){
      x_[0] = meas_package.raw_measurements_[0];
      x_[1] = meas_package.raw_measurements_[1];
      // x_[2] = 0.5;
      // x_[3] = 0.001;
      // x_[4] = 0.001;
      // P_<<1, 0, 0, 0, 0,
      //     0, 1, 0, 0, 0,
      //     0, 0, 0.5, 0, 0,
      //     0, 0, 0, 0.5, 0,
      //     0, 0, 0, 0, 1;
      //P_ = P_*1.2;
      cout<<"Laser Init"<<endl;
    }
    previous_timestamp_ = meas_package.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  } // !is_initialized_

  double delta_t = (meas_package.timestamp_ - previous_timestamp_)/1000000.0;

  // if (delta_t > 0.001){
  //   if ((meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) || (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_))
  //   {
  //     Prediction(delta_t);
  //   }
  //   else
  //   {
  //     printf("No meas\n" );
  //     return;
  //   }
  // }else{
  //   printf("Small dt\n" );
  //   return ;
  // }
  // previous_timestamp_ = meas_package.timestamp_;

  if (delta_t > 0.001 && meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
  {
    // printf("Radar UPdate\n" );
    Prediction(delta_t);
    UpdateRadar(meas_package);
    //return;
  }
  else if (delta_t > 0.001 && meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
  {
    // printf("Laser UPdate\n" );
    Prediction(delta_t);
    UpdateLidarKF(meas_package);
    //return;
  }
  else
  {
    printf("No meas\n" );
    return;

  }
  previous_timestamp_ = meas_package.timestamp_;
  return;
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
  // Augmentation
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // Set Process noise cov matrix
  MatrixXd Q = MatrixXd(2,2);
  Q << std_a_*std_a_, 0,
        0, std_yawdd_*std_yawdd_;

  //create augmented mean state
  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_;
  x_aug(5) = 0.0;
  x_aug(6) = 0.0;
  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug.bottomRightCorner(2,2) = Q;
  //create square root matrix
  MatrixXd P_aug_sqrt = P_aug.llt().matrixL();
  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  double scale_fac = sqrt(lambda_+n_aug_);
  for (int i=0; i<n_aug_; i++)
  {
      Xsig_aug.col(i+1) = x_aug + scale_fac*P_aug_sqrt.col(i);
      Xsig_aug.col(n_aug_+i+1) =  x_aug - scale_fac*P_aug_sqrt.col(i);
  }

  // Predict Sigma points

  //float v, yaw, yaw_rate, nu_a, nu_yawa;
  //VectorXd U = VectorXd(n_x_);
  //VectorXd N = VectorXd(n_x_);
  for (int i=0; i<2*n_aug_; i++)
  {
    float v = Xsig_aug(2,i); //fabs(Xsig_aug(2,i));
    //printf("V:%f\n", v);
    float yaw = Xsig_aug(3,i);
    float yaw_rate = Xsig_aug(4,i);
    float nu_a = Xsig_aug(5,i);
    float nu_yawa = Xsig_aug(6,i);

    VectorXd U = VectorXd(n_x_);
    VectorXd N = VectorXd(n_x_);

    U.fill(0.0);
    if (fabs(yaw_rate) > 0.001)
    {
      //cout <<"Using yaw rate" << endl;
      U(0) = v*(sin(yaw+yaw_rate*delta_t) - sin(yaw))/yaw_rate;
      U(1) = v*(cos(yaw)-cos(yaw+yaw_rate*delta_t))/yaw_rate;
      U(3) = yaw_rate*delta_t;
    }else
    {
      //cout <<"Yaw rate 0"<<endl;
      U(0) = v*cos(yaw)*delta_t;
      U(1) = v*sin(yaw)*delta_t;
    }

    N(0) = 0.5*delta_t*delta_t*cos(yaw)*nu_a;
    N(1) = 0.5*delta_t*delta_t*sin(yaw)*nu_a;
    N(2) = delta_t*nu_a;
    N(3) = 0.5*delta_t*delta_t*nu_yawa;
    N(4) = delta_t*nu_yawa;

    //write predicted sigma points into right column
    Xsig_pred_.col(i) = Xsig_aug.col(i).head(5) + U + N;
    //angle wrap
    // while (Xsig_pred_(3,i)> M_PI) Xsig_pred_(3,i)-=2.*M_PI;
    // while (Xsig_pred_(3,i)<-M_PI) Xsig_pred_(3,i)+=2.*M_PI;
    // float th = Xsig_pred_(3,i);
    // Xsig_pred_(3,i) = atan2(sin(th), cos(th));
  }

  VectorXd x = VectorXd(n_x_);
  x.fill(0.0);
  MatrixXd P = MatrixXd(n_x_, n_x_);
  P.fill(0.0);

  for (int i=0; i<2*n_aug_+1; i++){
      x = x + weights_(i)*Xsig_pred_.col(i);
  }
  //angle wrap
  float th = x(3);
  x(3) = atan2(sin(th), cos(th));
  // while (x(3)> M_PI) x(3)-=2.*M_PI;
  // while (x(3)<-M_PI) x(3)+=2.*M_PI;

  for (int i=0; i<2*n_aug_+1; i++){
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    //angle wrap
    float th = x_diff(3);
    x_diff(3) = atan2(sin(th), cos(th));
    // while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    // while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P = P + weights_(i)*x_diff*x_diff.transpose();
  }
  // x[2] = fabs(x[2]);
  x_ = x;
  P_ = P;
  // cout << "Delta_t"<< delta_t << endl;
  cout << "X (pred)" << x << endl;
  cout << "P (pred)" << P << endl;

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
  //Predict Lidar Measurements
  int n_z = 2;
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  Zsig = Xsig_pred_.block(0,0,n_z, 2*n_aug_+1);
  z_pred[0] = x_[0];
  z_pred[1] = x_[1];

  MatrixXd R = MatrixXd(n_z, n_z);
  R <<std_laspx_*std_laspx_, 0,
      0, std_laspy_*std_laspy_;

  S = CalculateS(z_pred, Zsig, R);

  // Cross Correlation Matrix T
  MatrixXd T = MatrixXd(n_x_, n_z);

  T = CalculateT(z_pred, Zsig);

  // Update measurements
  VectorXd z = VectorXd(n_z);
  z[0] = meas_package.raw_measurements_[0];
  z[1] = meas_package.raw_measurements_[1];

  //calculate Kalman gain K;
  MatrixXd K = MatrixXd(n_x_, n_z);
  MatrixXd Si = S.inverse();
  K = T*Si;

  //update state mean and covariance matrix
  VectorXd zdiff = z - z_pred;
  x_ = x_ + K*zdiff;
  //angle wrap
  float th = x_(3);
  x_(3) = atan2(sin(th), cos(th));

  // while (x_(3)> M_PI) x_(3)-=2.*M_PI;
  // while (x_(3)<-M_PI) x_(3)+=2.*M_PI;
  MatrixXd SKt;
  SKt = S*K.transpose();
  P_ = P_ - K*SKt;

  // NIS
  VectorXd SiD = Si*zdiff;
  NIS_laser_ = zdiff.adjoint()*SiD;

  // cout << "NIS L:" << NIS_laser_ << endl;

  cout << "X (L update)" << x_ << endl;
  cout << "P (L update)" << P_ << endl;

}
void UKF::UpdateLidarKF(MeasurementPackage meas_package) {
  int n_z = 2;
  MatrixXd H  = MatrixXd(2,5);
  H<<1,0,0,0,0,
     0,1,0,0,0;

  // Update measurements
  VectorXd z = VectorXd(n_z);
  z[0] = meas_package.raw_measurements_[0];
  z[1] = meas_package.raw_measurements_[1];

  VectorXd y = z - H*x_;

  MatrixXd R = MatrixXd(n_z, n_z);
  R <<std_laspx_*std_laspx_, 0,
      0, std_laspy_*std_laspy_;

  MatrixXd PHt = P_*H.transpose();
  MatrixXd S = H*PHt + R;
  MatrixXd Si = S.inverse();
  MatrixXd K = PHt*Si ;

  x_ = x_ + (K * y);
  //angle wrap
  float th = x_(3);
  x_(3) = atan2(sin(th), cos(th));

	MatrixXd I = MatrixXd::Identity(n_x_, n_x_);
	P_ = (I - K * H) * P_;

  // NIS
  VectorXd SiD = Si*y;
  NIS_laser_ = y.adjoint()*SiD;

  // cout << "NIS L:" << NIS_laser_ << endl;

  cout << "X (L update)" << x_ << endl;
  cout << "P (L update)" << P_ << endl;

  return;

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
  //Predict Radar Measurements
  int n_z = 3;
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);


  //transform sigma points into measurement space
  for (int i=0; i<2*n_aug_+1; i++){
      float px = Xsig_pred_(0,i);
      float py = Xsig_pred_(1,i);
      float v = Xsig_pred_(2,i);
      float yaw = Xsig_pred_(3,i);
      float c1 = sqrt(px*px + py*py);
      /*
      if (fabs(px) < 0.0001){
          px = (px>0 ? 0.0001 : -0.0001);
      }*/
      float c2 = px*cos(yaw)*v + py*sin(yaw)*v;
      if (fabs(c1) < 0.0001)
      {
        //return;
        //Zsig(2,i) = 0.0;
        printf("Small c1\n" );
        //return;
        c1 = 0.0001;
      }
      if (fabs(px) < 0.0001)
      {
        px = px>0? 0.0001 : -0.0001;
      }
      Zsig(0,i) = c1 ;
      Zsig(1,i) = atan2(py, px);
      Zsig(2,i) = c2/c1;
      // else
      // {
      //   Zsig(2,i) = c2/c1;
      // }
  }
  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for (int i=0; i<2*n_aug_+1; i++){
      z_pred = z_pred + weights_(i)*Zsig.col(i);
  }
  float th = z_pred(1);
  z_pred(1) = atan2(sin(th), cos(th));

  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_radr_*std_radr_,0,0,
      0,std_radphi_*std_radphi_,0,
      0,0,std_radrd_*std_radrd_;

  S = CalculateS(z_pred, Zsig, R);

  // Cross Correlation Matrix T
  MatrixXd T = MatrixXd(n_x_, n_z);

  T = CalculateT(z_pred, Zsig);

  // Update measurements
  VectorXd z = VectorXd(n_z);
  z[0] = meas_package.raw_measurements_[0];
  z[1] = meas_package.raw_measurements_[1];
  z[2] = meas_package.raw_measurements_[2];

  //calculate Kalman gain K;
  MatrixXd K = MatrixXd(n_x_, n_z);
  MatrixXd Si = S.inverse();
  K = T*Si;

  //update state mean and covariance matrix
  VectorXd zdiff = z-z_pred;
  th = zdiff(1);
  zdiff(1) = atan2(sin(th), cos(th));

  x_ = x_ + K*zdiff;
  //angle wrap
  th = x_(3);
  x_(3) = atan2(sin(th), cos(th));
  // while (x_(3)> M_PI) x_(3)-=2.*M_PI;
  // while (x_(3)<-M_PI) x_(3)+=2.*M_PI;
  MatrixXd SKt;
  SKt = S*K.transpose();
  P_ = P_ - K*SKt;

  // NIS
  VectorXd SiD = Si*zdiff;
  NIS_radar_ = zdiff.adjoint()*SiD;

  // cout << "NIS R:" << NIS_radar_ << endl;
  cout << "X (R update)" << x_ << endl;
  cout << "P (R update)" << P_ << endl;

}

MatrixXd UKF::CalculateS(const VectorXd& z_pred, const MatrixXd& Zsig, const MatrixXd& R) {
  int n_z = z_pred.rows();
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    if (n_z == 3)
    {
      //angle wrap for Radar meas
      float th = z_diff(1);
      z_diff(1) = atan2(sin(th), cos(th));
      // while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
      // while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    }


    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  S = S + R;
  return S;
}

MatrixXd UKF::CalculateT(const VectorXd& z_pred, const MatrixXd& Zsig)
{
  int n_z = z_pred.rows();
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i=0; i<2*n_aug_+1; i++)
  {
      VectorXd xdiff = Xsig_pred_.col(i)-x_;
      float th = xdiff(3);
      xdiff(3) = atan2(sin(th), cos(th));
      VectorXd zdiff = Zsig.col(i) - z_pred;
      // Angle wrap for Radar measurements
      if (n_z == 3)
      {
        th = zdiff(1);
        zdiff(1) = atan2(sin(th), cos(th));
        // while (zdiff(1) > M_PI) zdiff(1) -= 2.*M_PI;
        // while (zdiff(1) < -M_PI) zdiff(1) += 2.*M_PI;
      }


      Tc = Tc + weights_(i)*xdiff*zdiff.transpose();


  }
  return Tc;
}
