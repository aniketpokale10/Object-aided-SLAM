#define numEigV_c 10

#include <bits/stdc++.h>
#include <ceres/loss_function.h>
#include <ceres/iteration_callback.h>
#include <ceres/rotation.h>

// Contains definitions of various problem structs
#include "problemStructs.hpp"
// Contains various cost function struct specifications
#include "costFunctions.hpp"

using namespace std;
#define num_frames 10

int numEigV = 10;
// double Actual_angle_axis[num_frames][3] = {{0,0,0}, {0.3,0,0.2}, {0.5,0.2,0.0}};
double Actual_angle_axis[num_frames][3] ={{0,0,0},
                                          {-0.174485, 0, 0}, 
                                          {-0.349044, 0, 0}, 
                                          {-0.523606, 0, 0}, 
                                          {-0.69816, 0, 0},
                                          {-0.872614, 0, 0},
                                          {-1.04717, 0, 0}, 
                                          {-1.22175, 0, 0},
                                          {-1.3963, 0, 0},
                                          {-1.5708, 0, 0}};
                                           
double Actual_translation[num_frames][3] = {{2,0,5.2}, {2,0,5.4}, {2,0,5.6}, {2,0,5.8},{2,0,6}, {2,0,6.2}, {2,0,6.4},{2,0,6.6},{2,0,6.8},{2,0,7}};
double Actual_K[] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
double *X_bar;
double Actual_lambda[] = {-0.0346848,0.000143055,0.0428355 ,0.0533419 ,0.000372525,-0.0395549,0.00494395,-0.0947407,0.0189058 ,-0.0387325};
double *V; // Eigenvectors

// new Point cloud by JG
double *lambdas_jg;
double *k_jg;
double *rot_jg;
double *trans_jg;
double *observation_jg;
double *point_cloud_pts;
double *point_cloud;
int len_point_cloud_pts, num_point_cloud_pts;
double Actual_error;
double *point_cloud_without_error;
double *point_cloud_world;

MultiViewMultiObjectAdjustmentProblem totalProblem;

void pose_shapeAdjustment(double *observation, double *K, double *trans, double *rot, double *lambdas, bool point_cloud, double *rot_actual, double *trans_actual);

// To add error to the point cloud
double errorGen() {
	// return 0;
	return (double)((rand()%20-10.0)/50.0);
}

// To add errors to the image projection of the point cloud
double errorGenImage() {
  // return 0;
  return (double)((rand()%10-5.0)/1000.0);
}

// Computation of the error of the original point cloud(with errors) and the new optimized point cloud
double total_error (double *optimized_pt_cloud){
  double error = 0.0;
  double *pt;
  pt = new double[num_point_cloud_pts*3];
  copy ( point_cloud_without_error, point_cloud_without_error+num_point_cloud_pts*3, pt );
  for(int i = 0; i < num_point_cloud_pts; ++i) {
    float x=abs(pt[3*i+0]-optimized_pt_cloud[3*i+0]); float y=abs(pt[3*i+1]-optimized_pt_cloud[3*i+1]); float z=abs(pt[3*i+2]-optimized_pt_cloud[3*i+2]);
    error+= sqrt(x*x +y*y + z*z);
  }
  error = error/double(num_point_cloud_pts);
  return error;
}

// Calculation of the errors in the ground truth of the Rotation and the observed rotation -- experiment purpose
double R_error(double *rot, double * trans, double *rot_actual, double *trans_actual){
  double error = 0.0;
  for(int i=0;i<3;i++){
    error += abs(rot_actual[i]-rot[i])*abs(rot_actual[i]-rot[i]); 
  }
  error = sqrt(error);
  return error;
}

// Calculation of the errors in the ground truth of the Translation and the observed rotation -- experiment purpose
double T_error(double *rot, double * trans, double *rot_actual, double *trans_actual){
  double error = 0.0;
  for(int i=0;i<3;i++){
    error += abs(trans_actual[i]-trans[i])*abs(trans_actual[i]-trans[i]); 
  }
  error = sqrt(error);
  return error;
}

// Observations of the chair - generation
void generateObservation(double *observations, double *rot_actual, double *trans_actual) {
  double *points = new double[3*10]; // 10 = #points
  for(int i = 0; i < 30; i++) points[i] = X_bar[i];
  for(int vec = 0; vec < 10; vec++)
  for(int point = 0; point < 30; point++) {
    points[point] += Actual_lambda[vec]* V[30*vec + point];
  }
  for(int i = 0; i < 10; i++) {
    ceres::AngleAxisRotatePoint(rot_actual, &points[3*i], &points[3*i]);
    for(int j = 0; j < 3; j++) points[3*i+j]+=trans_actual[j];
  }
  for(int i = 0; i < 10; i++) {
    observations[2*i] = points[3*i]/points[3*i+2]+errorGenImage();
    observations[2*i+1] = points[3*i+1]/points[3*i+2]+errorGenImage();
  }
}

// Generate 2D observations of point cloud
void generatePointCloudPoints(double *rot, double *trans) {
  point_cloud_pts = new double[num_point_cloud_pts*2];
  for(int index=0; index<num_point_cloud_pts;index++) {
    double pt[] = {point_cloud_without_error[index*3+0]*1.0, point_cloud_without_error[index*3+1]*1.0, point_cloud_without_error[index*3+2]*1.0};
    ceres::AngleAxisRotatePoint(rot, &pt[0], &pt[0]);
    for(int i=0;i<3;i++) pt[i] += trans[i];
    point_cloud_pts[index*2] = pt[0]/pt[2]+errorGenImage();
    point_cloud_pts[index*2+1] = pt[1]/pt[2]+errorGenImage();
  }
}

//Generate point cloud in world coordinates
void generateWorldPoints(){
  point_cloud_world = new double[num_point_cloud_pts*3];
  for(int i=0;i<4;i++)
  for(int j=0;j<4;j++)
  for(int k=0;k<4;k++) {
    int index = k+4*j+16*i;
    point_cloud_world[index*3+0]=double(i)/5+0.01 + 2;
    point_cloud_world[index*3+1]=double(j)/5+0.01;
    point_cloud_world[index*3+2]=double(k)/5+0.01 + 5.2;
  }
}

// Convert point cloud from world to chair coordinate system
void pointCloudToChairFrame (double *rot, double *trans){
  point_cloud = new double[num_point_cloud_pts*3];
  point_cloud_without_error = new double[num_point_cloud_pts*3];

  for(int i=0;i<3;i++) rot[i] = -rot[i];

  for(int i=0;i<num_point_cloud_pts;i++){
    for(int j=0;j<3;j++) point_cloud_without_error[3*i+j] = point_cloud_world[3*i+j] - trans[j]; 

    ceres::AngleAxisRotatePoint(rot,point_cloud_without_error + 3*i,point_cloud_without_error + 3*i);
    
    double x=errorGen(), y=errorGen(), z=errorGen();
    point_cloud[3*i+0] = point_cloud_without_error[3*i+0] + x; point_cloud[3*i+1] = point_cloud_without_error[3*i+1] + y; point_cloud[3*i+2] = point_cloud_without_error[3*i+2] + z; 
    Actual_error+=sqrt(x*x+y*y+z*z);
  }  
  Actual_error = Actual_error/double(num_point_cloud_pts);  
}

// Setting up the env and initializing the important variables
void initialize() {

  Actual_error = 0.0;
  char inputFileName[] = "ceres_input_chair_multiViewmultiObjectAdjuster.txt";
  FILE *fptr=totalProblem.loadFile(inputFileName);
	if(fptr==NULL){
		std::cerr << "ERROR: Unable to open file " << inputFileName << std::endl;
    return;
	}
  V = totalProblem.getV();
  // Actual_lambda = totalProblem.getLambdas();
  cout << "Actual lambda" << endl;
  for(int i = 0; i < 10; i++) cout << Actual_lambda[i] << " ";cout<<endl;
  X_bar = totalProblem.getX_bar();
  //k_jg = Actual_K;
  k_jg = new double[9];
  for(int i = 0; i < 9; i++) k_jg[i] = Actual_K[i];
  rot_jg = new double[3];
  trans_jg = new double[3];
  lambdas_jg = new double[10];
  observation_jg = new double[20];
  // INITIALIZE
  for(int i = 0; i < 3; i++) {
    rot_jg[i] = 0;
    trans_jg[i] = 0;
  }

  rot_jg[0] = 0.1;
  rot_jg[2] = 0.1;
  trans_jg[2] = 5;
  trans_jg[0] = 1.6;

  double *obs = totalProblem.observations();
  copy ( obs, obs+20, observation_jg );
  generateObservation(observation_jg,Actual_angle_axis[0], Actual_translation[0]);
  double *lambdas = totalProblem.getLambdas();
  copy ( lambdas, lambdas+10, lambdas_jg );
}

// this is carried out on each frame
void reinitialize(double *rot_actual, double *trans_actual) {

  double *obs = totalProblem.observations();
  copy ( obs, obs+20, observation_jg );

  generateObservation(observation_jg, rot_actual, trans_actual); 
}

int main() {

  len_point_cloud_pts = 4;
  num_point_cloud_pts = len_point_cloud_pts*len_point_cloud_pts*len_point_cloud_pts;
  double *all_rots = new double[num_point_cloud_pts*3];
  double *all_trans = new double[num_point_cloud_pts*3];
  double *all_point_cloud_pts = new double[num_point_cloud_pts*3*num_frames];
  double *BArots = new double[num_frames-1];
  double *BAtrans = new double[num_frames-1];

  srand(time(0));
  initialize();  
	generateWorldPoints();


  ceres::Problem MultiViewProblem; //declare multiview problem to optimize for all frames at the end after 10 frames


  for(int i = 0; i < num_frames; i++) {
    if(i == 0){
      pose_shapeAdjustment(observation_jg, k_jg, trans_jg, rot_jg, lambdas_jg,false, Actual_angle_axis[i], Actual_translation[i]);
      // shape_adjustment(observation_jg, k_jg, trans_jg, trans_jg, lambdas_jg);
      pointCloudToChairFrame(rot_jg,trans_jg);

      cout << endl;
      cout << "Initial Point CLoud with error\n";
      for(int j = 0; j < num_point_cloud_pts; ++j) {
         cout << point_cloud[j*3+0] << " " << point_cloud[j*3+1] << " " << point_cloud[j*3+2] << endl;
       }

      cout<<"Initial point cloud ERROR = "<<Actual_error<<endl;

    }

    reinitialize(Actual_angle_axis[i], Actual_translation[i]);
    generatePointCloudPoints(Actual_angle_axis[i], Actual_translation[i]); //Generate 2D images of the point cloud points
    cout<<"\nCombined\n"<<endl;
    
    pose_shapeAdjustment(observation_jg, k_jg, trans_jg, rot_jg, lambdas_jg,false, Actual_angle_axis[i], Actual_translation[i]);
 
    copy(rot_jg, rot_jg+3, all_rots + 3*i);
    copy(trans_jg, trans_jg+3, all_trans + 3*i);
    copy(point_cloud_pts, point_cloud_pts + num_point_cloud_pts*2, all_point_cloud_pts + num_point_cloud_pts*2*i);

  }

  //Adding multiview error
  for(int i=0;i<num_frames;i++){
    for(int j = 0; j < num_point_cloud_pts; ++j) {

    ceres::CostFunction *MultiView = new ceres::AutoDiffCostFunction<BAError, 2,3,3,3>(
        new BAError(all_point_cloud_pts + num_point_cloud_pts*2*i + 2*j, k_jg));
      MultiViewProblem.AddResidualBlock(MultiView, new ceres::HuberLoss(0.1) , all_rots + 3*i, all_trans + 3*i, point_cloud + 3*j);
      // MultiViewProblem.SetParameterBlockConstant(rot_jg);
      // MultiViewProblem.SetParameterBlockConstant(trans_jg);
    }
  }

  ceres::Solver::Options MultiViewoptions;
  ceres::Solver::Summary MultiViewsummary;
  // MultiViewoptions.sparse_linear_algebra_library_type = ceres::SparseLinearAlgebraLibraryType::SUITE_SPARSE;
  MultiViewoptions.linear_solver_type = ceres::DENSE_SCHUR;
  MultiViewoptions.preconditioner_type = ceres::JACOBI;
  MultiViewoptions.use_inner_iterations = true;
  MultiViewoptions.minimizer_progress_to_stdout = false;
  ceres::Solve(MultiViewoptions, &MultiViewProblem, &MultiViewsummary);


  cout << endl;
  cout << "Point CLoud after MultiView optimization\n";
  for(int i = 0; i < num_point_cloud_pts; ++i) {
     cout << point_cloud[i*3+0] << " " << point_cloud[i*3+1] << " " << point_cloud[i*3+2] << endl;
   }
  cout << endl;

  cout << "Point cloud Error after MultiView= " << total_error(point_cloud) << endl;

  cout << "Initial error: " << Actual_error << endl;
}


void pose_shapeAdjustment(double *observations, double *K, double *trans, double *rot, double *lambdas, bool PointCloud, double *rot_actual, double *trans_actual)
{
    ceres::Problem pose_shapeProblem;
    const int numPointCloud = 64;
    double *pt;

    int numObs = 10;
    for(int i = 0; i < numObs; ++i)
    {
      double *curEigVec = new double[3*10];
      for(int j = 0; j < numEigV; ++j){
        curEigVec[3*j+0] = V[3*numObs*j + 3*i + 0];
        curEigVec[3*j+1] = V[3*numObs*j + 3*i + 1];
        curEigVec[3*j+2] = V[3*numObs*j + 3*i + 2];
      }

        ceres::CostFunction *PnP = new ceres::AutoDiffCostFunction<PnPError, 2,3,3>(
        new PnPError(X_bar+3*i, observations+2*i, curEigVec, K, 1, lambdas, numEigV));
      pose_shapeProblem.AddResidualBlock(PnP, new ceres::HuberLoss(0.1) , rot, trans);

      ceres::CostFunction *lambdaError = new ceres::AutoDiffCostFunction<LambdaReprojectionError, 2,3,3,numEigV_c>(
        new LambdaReprojectionError(X_bar+3*i, observations+2*i, curEigVec, K, 1, numEigV));
      pose_shapeProblem.AddResidualBlock(lambdaError, new ceres::HuberLoss(0.6) , rot, trans, lambdas);

    }

    ceres::Solver::Options pose_shapeoptions;
    ceres::Solver::Summary pose_shapesummary;
    pose_shapeoptions.linear_solver_type = ceres::DENSE_SCHUR;
    pose_shapeoptions.preconditioner_type = ceres::JACOBI;
    pose_shapeoptions.use_inner_iterations = true;
    pose_shapeoptions.minimizer_progress_to_stdout = false;
    ceres::Solve(pose_shapeoptions, &pose_shapeProblem, &pose_shapesummary);

    // cout << "Rotational Matrix" << endl;
    // for(int i = 0; i < 3; i++) cout << rot[i] << " "; cout << endl;

    // cout << "Translation" << endl;
    // cout << trans[0] << " " << trans[1] << " " << trans[2] << std::endl;

    cout << "lambdas:\n";
    for(int i =0; i <  numEigV; ++i)
    cout << lambdas[i] << " "; cout<<endl;
    copy(lambdas,lambdas+10,lambdas_jg); // Current values of lambdas used as initialization for next frame 

    if(PointCloud){

      pt = new double[numPointCloud*3];
      copy ( point_cloud, point_cloud+num_point_cloud_pts*3, pt );
      for(int i = 0; i < num_point_cloud_pts; ++i) {
        ceres::CostFunction *PointCloudfunct = new ceres::AutoDiffCostFunction<PointCloudError, 2, 3>(
            new PointCloudError(pt + 3 * i, point_cloud_pts + 2 * i, K, rot, trans));
        // pose_shapeProblem.AddResidualBlock(PointCloud, new ceres::HuberLoss(0.9), rot, trans, pt + 3 * i);
        pose_shapeProblem.AddResidualBlock(PointCloudfunct, new ceres::HuberLoss(0.9), pt + 3 * i);
       }

       // pose_shapeoptions.minimizer_progress_to_stdout = true;
       ceres::Solve(pose_shapeoptions, &pose_shapeProblem, &pose_shapesummary);


      cout << endl;
      cout << "Point CLoud after optimization\n";
      for(int i = 0; i < numPointCloud; ++i) {
         cout << pt[i*3+0] << " " << pt[i*3+1] << " " << pt[i*3+2] << endl;
       }
      cout << endl;

      cout << "Point cloud Error = " << total_error(pt) << endl;
      copy ( pt, pt+num_point_cloud_pts*3, point_cloud );
    }
    cout << "\n" << "R error = " << R_error(rot,trans, rot_actual, trans_actual) << endl << endl;
    cout << "\n" << "T error = " << T_error(rot,trans, rot_actual, trans_actual) << endl << endl;

    copy ( rot, rot+3, rot_jg );
    copy ( trans, trans+3, trans_jg );
}