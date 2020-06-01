#include <ceres/ceres.h>
#include <ceres/rotation.h>
// #include <fstream>
// #include <iostream>



class MultiViewMultiObjectAdjustmentProblem{

public:

	// Return the number of views (each view can have multiple observation)
	int getNumViews() const { return numViews_; }
	// Return the number of observations
	int getNumObs() const { return numObs_; }
	// Return the number of keypoints observed
	int getNumPts() const { return numPts_; }
	// Return number of eigen vectors
	int getNumEigV() { return numEigV_; }
	// Return the total number of observations of all frames
	int getTotalObs() const { return totalObs_; }
	// Return the total number of observations of all frames
	int getTotalObject() { return totalObject_; }
	// Return the center of the car
	double* getCarCenter() { return carCenter_; }
	// Return a pointer to the observation vector
	double* observations() const { return observations_; }
	// Return a pointer to the observation weights vector
	double* observationWeights() const { return observationWeights_; }
	// Return a pointer to the camera intrinsics
	double* getK() const { return K_; }
	// Return a pointer to the mean 3D locations
	double* getX_bar() const { return X_bar_; }
	// Return a pointer to the top numEigV_ eigenvectors
	double* getV() { return V_; }
	// Return a pointer to the weights (lambdas)
	double* getLambdas() { return lambdas_; }
	// Return Initial Transalation
	double* getTransIni() { return trans_ini_; }
	// Return Initial Rotation
	double* getRotIni() { return rot_ini_; }
	// Return the pointCloudDataPoints
	double* getPointCloudPoints() { return point_cloud_pts; }

	// Read bare bone data from input file
	FILE* loadFile(const char *fileName){
		FILE *fptr = fopen(fileName, "r");
		if(fptr == NULL){
			return fptr;
		}

		fscanfOrDie(fptr, "%d", &numViews_);
		fscanfOrDie(fptr, "%d", &numPts_);
		fscanfOrDie(fptr, "%d", &numObs_);
		fscanfOrDie(fptr, "%d", &numEigV_);
		fscanfOrDie(fptr, "%d", &totalObs_);
		fscanfOrDie(fptr, "%d", &totalObject_);

		K_ = new double[9];
		for(int i = 0; i < 9; ++i){
			fscanfOrDie(fptr, "%lf", K_ + i);
		}


		X_bar_ = new double[3*numObs_];
		for(int i = 0; i < numObs_; ++i){
			for(int j = 0; j < 3; ++j){
				fscanfOrDie(fptr, "%lf", X_bar_ + i*3 + j);
			}
		}

		V_ = new double[numEigV_*3*numPts_];
		for(int i = 0; i < numEigV_; ++i){
			for(int j = 0; j < numPts_; ++j){
				fscanfOrDie(fptr, "%lf", V_ + i*3*numPts_ + 3*j + 0);
				fscanfOrDie(fptr, "%lf", V_ + i*3*numPts_ + 3*j + 1);
				fscanfOrDie(fptr, "%lf", V_ + i*3*numPts_ + 3*j + 2);
			}
		}
		carCenter_ = new double[3*totalObject_];
		trans_ini_ = new double[3*totalObject_];
		rot_ini_ = new double[9*totalObject_];
		lambdas_ = new double[numEigV_*totalObject_];
		observations_ = new double[2*numObs_*totalObs_];
		observationWeights_ = new double[numObs_*totalObs_];

		observation_pointer=0;
		observation_weight_pointer=0;

		return fptr;

	}

	FILE* readNextObs(FILE *fptr,int *ans)
	{

		int frame_id, object_id, initial_flag, lambdas_flag;

		fscanfOrDie(fptr, "%d", &frame_id); //Frame Number
		fscanfOrDie(fptr, "%d", &object_id ); //Object Id
		fscanfOrDie(fptr, "%d", &initial_flag ); //Does Rin and Tin requires Re-initialisation
		fscanfOrDie(fptr, "%d", &lambdas_flag ); //Does lambdas need initialization

		ans[0]=frame_id;
		ans[1]=object_id;
		ans[2]=initial_flag;
		ans[3]=lambdas_flag;

		// Observations 2D

		for(int i = 0; i < numObs_; ++i){
			for(int j = 0; j < 2; ++j){
				fscanfOrDie(fptr, "%lf", observations_ + observation_pointer + i*2 + j);
			}
		}


		observation_pointer=observation_pointer+2*numObs_;

		// Observation weights

		for(int i = 0; i < numObs_; ++i){
			fscanfOrDie(fptr, "%lf", observationWeights_ + observation_weight_pointer + i);
		}

		observation_weight_pointer=observation_weight_pointer+numObs_;




		if (initial_flag==1)
		{

			fscanfOrDie(fptr, "%lf", trans_ini_+3*object_id+0);
			fscanfOrDie(fptr, "%lf", trans_ini_+3*object_id+1);
			fscanfOrDie(fptr, "%lf", trans_ini_+3*object_id+2);



			for(int i = 0; i < 9; ++i){
				fscanfOrDie(fptr, "%lf", rot_ini_+9*object_id+ i);
			}

			fscanfOrDie(fptr, "%lf", carCenter_ + 3*object_id + 0);
			fscanfOrDie(fptr, "%lf", carCenter_ + 3*object_id + 1);
			fscanfOrDie(fptr, "%lf", carCenter_ + 3*object_id + 2);


		}


		if (lambdas_flag==1)
		{
			for(int i = 0; i < numEigV_; ++i){
				fscanfOrDie(fptr, "%lf", lambdas_+ (numEigV_*object_id) + i);
			}
		}

		return fptr;

	}

	FILE* loadPointCloud(const char *fileName){
		FILE *fptr = fopen(fileName, "r");
		if(fptr == NULL){
			return fptr;
		}
		for(int i = 0; i < 192256; i++) {
			fscanfOrDie(fptr, "%lf", point_cloud_pts + i);
		}
	}



private:

	// Helper function to read in one value to a text file
	template <typename T>
	void fscanfOrDie(FILE *fptr, const char *format, T *value){
		int numScanned = fscanf(fptr, format, value);
		if(numScanned != 1){
			LOG(FATAL) << "Invalid data file";
		}
	}

	// Private variables

	// Number of views
	int numViews_;
	// Number of keypoints
	int numPts_;
	// Number of observations
	int numObs_;
	// Number of eigen vectors
	int numEigV_;
	// Number of Total Observation in all frames
	int totalObs_;
	// Maximum Number of Objects in environment
	int totalObject_;

	// Car Center
	double *carCenter_;
	// Camera intrinsics
	double *K_;
	// Observation vector
	double *observations_;
	// Observation weight vector
	double *observationWeights_;
	// 3D point
	double *X_bar_;

	// Top numEigV_ eigenvectors for the shape, i.e., the deformation basis vectors
	double *V_;

	// Weights for the eigenvectors
	double *lambdas_;

	// Rotation estimate (would be initialized with some random value)
	double *rot_;
	// Translation estimate (would be initialized with some random value)
	double *trans_;
	// Initial Translation Estimate (Useless)
	double *trans_ini_;
	// Initial Rotation Estimate (Useless)
	double *rot_ini_;
	// Point Cloud points
	double *point_cloud_pts;


	//pointer
	int observation_pointer;
	int observation_weight_pointer;
};



class MultiViewAdjustmentProblem{

public:

	// Return the number of views
	int getNumViews() const { return numViews_; }
	// Return the number of observations
	int getNumObs() const { return numObs_; }
	// Return the number of keypoints observed
	int getNumPts() const { return numPts_; }
	// Return number of eigen vectors
	int getNumEigV() { return numEigV_; }
	// Return the center of the car
	double* getCarCenter() { return carCenter_; }
	// Return a pointer to the observation vector
	double* observations() const { return observations_; }
	// Return a pointer to the observation weights vector
	double* observationWeights() const { return observationWeights_; }
	// Return a pointer to the camera intrinsics
	double* getK() const { return K_; }
	// Return a pointer to the mean 3D locations
	double* getX_bar() const { return X_bar_; }
	// Return a pointer to the top numEigV_ eigenvectors
	double* getV() { return V_; }
	// Return a pointer to the weights (lambdas)
	double* getLambdas() { return lambdas_; }
	// Return Initial Transalation
	double* getTransIni() { return trans_ini_; }
	// Return Initial Rotation
	double* getRotIni() { return rot_ini_; }

	// Read data from input file
	bool loadFile(const char *fileName){
		FILE *fptr = fopen(fileName, "r");
		if(fptr == NULL){
			return false;
		}

		// numViews, numPts, numObs, numFaces
		fscanfOrDie(fptr, "%d", &numViews_);
		fscanfOrDie(fptr, "%d", &numPts_);
		fscanfOrDie(fptr, "%d", &numObs_);
		fscanfOrDie(fptr, "%d", &numEigV_);

		// Center of the car
		carCenter_ = new double[3];
		fscanfOrDie(fptr, "%lf", carCenter_+0);
		fscanfOrDie(fptr, "%lf", carCenter_+1);
		fscanfOrDie(fptr, "%lf", carCenter_+2);


		// translation Initialization
		trans_ini_ = new double[3];

		fscanfOrDie(fptr, "%lf", trans_ini_+0);
		fscanfOrDie(fptr, "%lf", trans_ini_+1);
		fscanfOrDie(fptr, "%lf", trans_ini_+2);

		// Rotation Initialization
		rot_ini_ = new double[9];

		for(int i = 0; i < 9; ++i){
			fscanfOrDie(fptr, "%lf", rot_ini_ + i);
		}

		// Camera intrinsics
		K_ = new double[9];
		for(int i = 0; i < 9; ++i){
			fscanfOrDie(fptr, "%lf", K_ + i);
		}

		// Observations
		observations_ = new double[2*numObs_*numViews_];
		for(int i = 0; i < numObs_*numViews_; ++i){
			for(int j = 0; j < 2; ++j){
				fscanfOrDie(fptr, "%lf", observations_ + i*2 + j);
			}
		}

		// Observation weights
		observationWeights_ = new double[numObs_*numViews_];
		for(int i = 0; i < numObs_*numViews_; ++i){
			fscanfOrDie(fptr, "%lf", observationWeights_ + i);
		}

		// Mean locations
		X_bar_ = new double[3*numObs_];
		for(int i = 0; i < numObs_; ++i){
			for(int j = 0; j < 3; ++j){
				fscanfOrDie(fptr, "%lf", X_bar_ + i*3 + j);
			}
		}

		// Read in the top numEigV_ eigenvectors for the shape
		// Size allocation: numEigV_ vecs * 3 coordinates per vex * 14 keypoints (numPts_)
		V_ = new double[numEigV_*3*numPts_];
		for(int i = 0; i < numEigV_; ++i){
			for(int j = 0; j < numPts_; ++j){
				fscanfOrDie(fptr, "%lf", V_ + i*3*numPts_ + 3*j + 0);
				fscanfOrDie(fptr, "%lf", V_ + i*3*numPts_ + 3*j + 1);
				fscanfOrDie(fptr, "%lf", V_ + i*3*numPts_ + 3*j + 2);
			}
		}

		// Read in the initial values for lambdas
		lambdas_ = new double[numEigV_];
		for(int i = 0; i < numEigV_; ++i){
			fscanfOrDie(fptr, "%lf", lambdas_ + i);
		}
/*
		// Read in the rotation estimate (from PnP) (column-major ordered rotation matrix)
		rot_ = new double[9*numViews_];
		for(int i = 0; i < 9*numViews_; ++i){
			fscanfOrDie(fptr, "%lf", rot_ + i);
		}

		// Read in the translation estimate (from PnP)
		trans_ = new double[3*numViews_];
		for(int i = 0; i < 3*numViews_; ++i){
			fscanfOrDie(fptr, "%lf", trans_ + i);
		}

*/

		// // Printing out data (for verification)
		// std::cout << "numViews: " << numViews_ << std::endl;
		// std::cout << "numPoints: " << numPts_ << std::endl;
		// std::cout << "numObs: " << numObs_ << std::endl;
		// std::cout << "K: " << K_[0] << " " << K_[1] << " " << K_[2] << " " << K_[3] << " " \
		// 	<< K_[4] << " " << K_[5] << " " << K_[6] << " " << K_[7] << " " << K_[8] << std::endl;
		// for(int i = 0; i < numObs_; ++i){
		// 	std::cout << "Obs: " << observations_[0+2*i] << " " << observations_[1+2*i] << std::endl;
		// }
		// for(int i = 0; i < numObs_; ++i){
		// 	std::cout << "ObsWeight: " << observationWeights_[i] << std::endl;
		// }
		// for(int i = 0; i < numObs_; ++i){
		// 	std::cout << "3D Point: " << X_bar_[0+3*i] << " " << X_bar_[1+3*i] << " " \
		// 	<< X_bar_[2+3*i] << std::endl;
		// }

		return true;

	}

private:

	// Helper function to read in one value to a text file
	template <typename T>
	void fscanfOrDie(FILE *fptr, const char *format, T *value){
		int numScanned = fscanf(fptr, format, value);
		if(numScanned != 1){
			LOG(FATAL) << "Invalid data file";
		}
	}

	// Private variables

	// Number of views
	int numViews_;
	// Number of keypoints
	int numPts_;
	// Number of observations
	int numObs_;
	// Number of eigen vectors
	int numEigV_;

	// Car Center
	double *carCenter_;
	// Camera intrinsics
	double *K_;
	// Observation vector
	double *observations_;
	// Observation weight vector
	double *observationWeights_;
	// 3D point
	double *X_bar_;

	// Top numEigV_ eigenvectors for the shape, i.e., the deformation basis vectors
	double *V_;

	// Weights for the eigenvectors
	double *lambdas_;

	// Rotation estimate (would be initialized with some random value)
	double *rot_;
	// Translation estimate (would be initialized with some random value)
	double *trans_;
	// Initial Translation Estimate (Useless)
	double *trans_ini_;
	// Initial Rotation Estimate (Useless)
	double *rot_ini_;

};



class SingleViewPoseAdjustmentProblem{

public:

	// Return the number of views
	int getNumViews() const { return numViews_; }
	// Return the number of observations
	int getNumObs() const { return numObs_; }
	// Return the number of keypoints observed
	int getNumPts() const { return numPts_; }
	// Return number of eigen vectors
	int getNumEigV() { return numEigV_; }
	// Return the center of the car
	double* getCarCenter() { return carCenter_; }
	// Return the height of the car
/*	double getCarHeight() const { return h_; }
	// Return the width of the car
	double getCarWidth() const { return w_; }
	// Return the length of the car
	double getCarLength() const { return l_; }
*/	// Return a pointer to the observation vector
	double* observations() const { return observations_; }
	// Return a pointer to the observation weights vector
	double* observationWeights() const { return observationWeights_; }
	// Return a pointer to the camera intrinsics
	double* getK() const { return K_; }
	// Return a pointer to the mean 3D locations
	double* getX_bar() const { return X_bar_; }
	// Return a pointer to the top numEigV_ eigenvectors
	double* getV() { return V_; }
	// Return a pointer to the weights (lambdas)
	double* getLambdas() { return lambdas_; }
	// Return Translation Initialization
	double* getTransIni() { return trans_ini_; }
	// Return Rotation Initialization
	double* getRotIni() { return rot_ini_; }

	// Read data from input file
	bool loadFile(const char *fileName){
		FILE *fptr = fopen(fileName, "r");
		if(fptr == NULL){
			return false;
		}

		// numViews, numPts, numObs
		fscanfOrDie(fptr, "%d", &numViews_);
		fscanfOrDie(fptr, "%d", &numPts_);
		fscanfOrDie(fptr, "%d", &numObs_);
		fscanfOrDie(fptr, "%d", &numEigV_);

		// Center of the car
		carCenter_ = new double[3];
		fscanfOrDie(fptr, "%lf", carCenter_+0);
		fscanfOrDie(fptr, "%lf", carCenter_+1);
		fscanfOrDie(fptr, "%lf", carCenter_+2);
/*
		// Height, Width, and Length of the Car
		fscanfOrDie(fptr, "%lf", &h_);
		fscanfOrDie(fptr, "%lf", &w_);
		fscanfOrDie(fptr, "%lf", &l_);
*/
		// Translation Initialization
		trans_ini_ = new double[3];

		fscanfOrDie(fptr, "%lf", trans_ini_+0);
		fscanfOrDie(fptr, "%lf", trans_ini_+1);
		fscanfOrDie(fptr, "%lf", trans_ini_+2);


		// Rotation Initialization
		rot_ini_ = new double[9];

		for(int i = 0; i < 9; ++i){
			fscanfOrDie(fptr, "%lf", rot_ini_ + i);
		}


		// Camera intrinsics
/*		K_ = new double[9*numViews_];
		for(int i = 0; i < 9*numViews_; ++i){
			fscanfOrDie(fptr, "%lf", K_ + i);
		}
*/
		K_ = new double[9];
		for(int i = 0; i < 9; ++i){
			fscanfOrDie(fptr, "%lf", K_ + i);
		}

		// Observations
		observations_ = new double[2*numObs_];
		for(int i = 0; i < numObs_; ++i){
			for(int j = 0; j < 2; ++j){
				fscanfOrDie(fptr, "%lf", observations_ + i*2 + j);
			}
		}

		// Observation weights
		observationWeights_ = new double[numObs_];
		for(int i = 0; i < numObs_; ++i){
			fscanfOrDie(fptr, "%lf", observationWeights_ + i);
		}

		// Mean locations
		X_bar_ = new double[3*numObs_];
		for(int i = 0; i < numObs_; ++i){
			for(int j = 0; j < 3; ++j){
				fscanfOrDie(fptr, "%lf", X_bar_ + i*3 + j);
			}
		}

		// Read in the top numEigV_ eigenvectors for the shape
		// Size allocation: numEigV_ vecs * 3 coordinates per vex * 14 keypoints (numPts_)
		V_ = new double[numEigV_*3*numPts_];
		for(int i = 0; i < numEigV_; ++i){
			for(int j = 0; j < numPts_; ++j){
				fscanfOrDie(fptr, "%lf", V_ + i*3*numPts_ + 3*j + 0);
				fscanfOrDie(fptr, "%lf", V_ + i*3*numPts_ + 3*j + 1);
				fscanfOrDie(fptr, "%lf", V_ + i*3*numPts_ + 3*j + 2);
			}
		}

		// Read in the initial values for lambdas
		lambdas_ = new double[numEigV_];
		for(int i = 0; i < numEigV_; ++i){
			fscanfOrDie(fptr, "%lf", lambdas_ + i);
		}



		// // Printing out data (for verification)
		// std::cout << "numViews: " << numViews_ << std::endl;
		// std::cout << "numPoints: " << numPts_ << std::endl;
		// std::cout << "numObs: " << numObs_ << std::endl;
		// std::cout << "K: " << K_[0] << " " << K_[1] << " " << K_[2] << " " << K_[3] << " " \
		// 	<< K_[4] << " " << K_[5] << " " << K_[6] << " " << K_[7] << " " << K_[8] << std::endl;
		// for(int i = 0; i < numObs_; ++i){
		// 	std::cout << "Obs: " << observations_[0+2*i] << " " << observations_[1+2*i] << std::endl;
		// }
		// for(int i = 0; i < numObs_; ++i){
		// 	std::cout << "ObsWeight: " << observationWeights_[i] << std::endl;
		// }
		// for(int i = 0; i < numObs_; ++i){
		// 	std::cout << "3D Point: " << X_bar_[0+3*i] << " " << X_bar_[1+3*i] << " " \
		// 	<< X_bar_[2+3*i] << std::endl;
		// }

		return true;

	}

private:

	// Helper function to read in one value to a text file
	template <typename T>
	void fscanfOrDie(FILE *fptr, const char *format, T *value){
		int numScanned = fscanf(fptr, format, value);
		if(numScanned != 1){
			LOG(FATAL) << "Invalid data file";
		}
	}

	// Private variables

	// Number of views
	int numViews_;
	// Number of keypoints
	int numPts_;
	// Number of observations
	int numObs_;
	// Number of eigen vectors
	int numEigV_;

	// Center of the car
	double *carCenter_;
	// Dimensions of the car
	//double h_, w_, l_;

	// Camera intrinsics
	double *K_;
	// Observation vector
	double *observations_;
	// Observation weight vector
	double *observationWeights_;
	// 3D point
	double *X_bar_;

	// Top numEigV_ eigenvectors for the shape, i.e., the deformation basis vectors
	double *V_;

	// Weights for the eigenvectors
	double *lambdas_;

	// Initial Transaltion
	double *trans_ini_;

	// Initial Rotation
	double *rot_ini_;

};




// Read a shape adjustment problem (a single view one)
class SingleViewShapeAdjustmentProblem{

public:

	// Return the number of views
	int getNumViews() const { return numViews_; }
	// Return the number of observations
	int getNumObs() const { return numObs_; }
	// Return the number of keypoints observed
	int getNumPts() const { return numPts_; }
	// Return number of eigen vectors
	int getNumEigV() { return numEigV_; }
	// Return the center of the car
	double* getCarCenter() { return carCenter_; }
	// Return the height of the car
/*	double getCarHeight() const { return h_; }
	// Return the width of the car
	double getCarWidth() const { return w_; }
	// Return the length of the car
	double getCarLength() const { return l_; }
*/	// Return a pointer to the observation vector
	double* observations() const { return observations_; }
	// Return a pointer to the observation weights vector
	double* observationWeights() const { return observationWeights_; }
	// Return a pointer to the camera intrinsics
	double* getK() const { return K_; }
	// Return a pointer to the mean 3D locations
	double* getX_bar() const { return X_bar_; }
	// Return a pointer to the top numEigV_ eigenvectors
	double* getV() { return V_; }
	// Return a pointer to the weights (lambdas)
	double* getLambdas() { return lambdas_; }
	// Return a pointer to the rotation estimated (from PnP)
	double* getRot() { return rot_; }
	// Return a pointer to the translation estimated (from PnP)
	double* getTrans() { return trans_; }
	// Return Initial Transalation
	double* getTransIni() { return trans_ini_; }
	// Return Initial Rotation
	double* getRotIni() { return rot_ini_; }

	// Read data from input file
	bool loadFile(const char *fileName){
		FILE *fptr = fopen(fileName, "r");
		if(fptr == NULL){
			return false;
		}

		// numViews, numPts, numObs, numFaces
		fscanfOrDie(fptr, "%d", &numViews_);
		fscanfOrDie(fptr, "%d", &numPts_);
		fscanfOrDie(fptr, "%d", &numObs_);
		fscanfOrDie(fptr, "%d", &numEigV_);

		// Center of the car
		carCenter_ = new double[3];
		fscanfOrDie(fptr, "%lf", carCenter_+0);
		fscanfOrDie(fptr, "%lf", carCenter_+1);
		fscanfOrDie(fptr, "%lf", carCenter_+2);

		// Height, Width, and Length of the Car
/*		fscanfOrDie(fptr, "%lf", &h_);
		fscanfOrDie(fptr, "%lf", &w_);
		fscanfOrDie(fptr, "%lf", &l_);
*/

		// translation Initialization
		trans_ini_ = new double[3];

		fscanfOrDie(fptr, "%lf", trans_ini_+0);
		fscanfOrDie(fptr, "%lf", trans_ini_+1);
		fscanfOrDie(fptr, "%lf", trans_ini_+2);

		// Rotation Initialization
		rot_ini_ = new double[9];

		for(int i = 0; i < 9; ++i){
			fscanfOrDie(fptr, "%lf", rot_ini_ + i);
		}

		// Camera intrinsics
		K_ = new double[9];
		for(int i = 0; i < 9; ++i){
			fscanfOrDie(fptr, "%lf", K_ + i);
		}

		// Observations
		observations_ = new double[2*numObs_*numViews_];
		for(int i = 0; i < numObs_*numViews_; ++i){
			for(int j = 0; j < 2; ++j){
				fscanfOrDie(fptr, "%lf", observations_ + i*2 + j);
			}
		}

		// Observation weights
		observationWeights_ = new double[numObs_*numViews_];
		for(int i = 0; i < numObs_*numViews_; ++i){
			fscanfOrDie(fptr, "%lf", observationWeights_ + i);
		}

		// Mean locations
		X_bar_ = new double[3*numObs_];
		for(int i = 0; i < numObs_; ++i){
			for(int j = 0; j < 3; ++j){
				fscanfOrDie(fptr, "%lf", X_bar_ + i*3 + j);
			}
		}

		// Read in the top numEigV_ eigenvectors for the shape
		// Size allocation: numEigV_ vecs * 3 coordinates per vex * 14 keypoints (numPts_)
		V_ = new double[numEigV_*3*numPts_];
		for(int i = 0; i < numEigV_; ++i){
			for(int j = 0; j < numPts_; ++j){
				fscanfOrDie(fptr, "%lf", V_ + i*3*numPts_ + 3*j + 0);
				fscanfOrDie(fptr, "%lf", V_ + i*3*numPts_ + 3*j + 1);
				fscanfOrDie(fptr, "%lf", V_ + i*3*numPts_ + 3*j + 2);
			}
		}

		// Read in the initial values for lambdas
		lambdas_ = new double[numEigV_];
		for(int i = 0; i < numEigV_; ++i){
			fscanfOrDie(fptr, "%lf", lambdas_ + i);
		}

		// Read in the rotation estimate (from PnP) (column-major ordered rotation matrix)
		rot_ = new double[9*numViews_];
		for(int i = 0; i < 9*numViews_; ++i){
			fscanfOrDie(fptr, "%lf", rot_ + i);
		}

		// Read in the translation estimate (from PnP)
		trans_ = new double[3*numViews_];
		for(int i = 0; i < 3*numViews_; ++i){
			fscanfOrDie(fptr, "%lf", trans_ + i);
		}



		// // Printing out data (for verification)
		// std::cout << "numViews: " << numViews_ << std::endl;
		// std::cout << "numPoints: " << numPts_ << std::endl;
		// std::cout << "numObs: " << numObs_ << std::endl;
		// std::cout << "K: " << K_[0] << " " << K_[1] << " " << K_[2] << " " << K_[3] << " " \
		// 	<< K_[4] << " " << K_[5] << " " << K_[6] << " " << K_[7] << " " << K_[8] << std::endl;
		// for(int i = 0; i < numObs_; ++i){
		// 	std::cout << "Obs: " << observations_[0+2*i] << " " << observations_[1+2*i] << std::endl;
		// }
		// for(int i = 0; i < numObs_; ++i){
		// 	std::cout << "ObsWeight: " << observationWeights_[i] << std::endl;
		// }
		// for(int i = 0; i < numObs_; ++i){
		// 	std::cout << "3D Point: " << X_bar_[0+3*i] << " " << X_bar_[1+3*i] << " " \
		// 	<< X_bar_[2+3*i] << std::endl;
		// }

		return true;

	}

private:

	// Helper function to read in one value to a text file
	template <typename T>
	void fscanfOrDie(FILE *fptr, const char *format, T *value){
		int numScanned = fscanf(fptr, format, value);
		if(numScanned != 1){
			LOG(FATAL) << "Invalid data file";
		}
	}

	// Private variables

	// Number of views
	int numViews_;
	// Number of keypoints
	int numPts_;
	// Number of observations
	int numObs_;
	// Number of eigen vectors
	int numEigV_;

	// Center of the car
	double *carCenter_;
	// Dimensions of the car
//	double h_, w_, l_;

	// Camera intrinsics
	double *K_;
	// Observation vector
	double *observations_;
	// Observation weight vector
	double *observationWeights_;
	// 3D point
	double *X_bar_;

	// Top numEigV_ eigenvectors for the shape, i.e., the deformation basis vectors
	double *V_;

	// Weights for the eigenvectors
	double *lambdas_;

	// Rotation estimate (from PnP)
	double *rot_;
	// Translation estimate (from PnP)
	double *trans_;
	// Initial Translation Estimate (Useless)
	double *trans_ini_;
	// Initial Rotation Estimate (Useless)
	double *rot_ini_;

};
