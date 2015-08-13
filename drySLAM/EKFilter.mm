//
//  EKFilter.m
//  EKFMonoSLAM
//
//  Created by Shreenivaas Devarajan on 3/12/15.
//  Copyright (c) 2015 dryslam. All rights reserved.
//

#import "EKFilter.hpp"
#import <opencv2/opencv.hpp>



//@implementation EKFilter
//Setting initial Velocity Values
//double* x_k_k = (double*)malloc(13*sizeof(double));
//double* p_k_k = (double*)malloc(13*13*sizeof(double));
//double delta_t = 1.0;
//double linear_acceleration_noise_covariance;
//double angular_acceleration_noise_covariance;

//Global Function Declarations

double* v2q(double omega1_mod, double omega2_mod, double omega3_mod, double *qwt);

 double * dq3_by_dq1(double *qOld, double* dq3_by_dq1RES);

//double * dqomegadt_by_domega(double *dqomegadt_by_domegaRES, double omega1, double omega2, double omega3 ,double delta_t);

 double* qprod(double* qWR, double* v2qRES, double* qprodRES);


/***** END OF GLOBAL FUNCTION DECLARATIONS *****/

double* generate_random_6D_sphere(double* randSphere6D, int nPointsRand){
    
    // Actual Code for the function begins here.
    
    if(randSphere6D == NULL)
        return FALSE;
    
    if(nPointsRand < 0)
        return FALSE;
    
    nPointsRand = 1000;
    
    return randSphere6D;
    
    //return true;
}

 void initialize_x_and_p() {
     
     //global variable dependency needs to be fixed
     /*int i, j, k;
     for (i = 0; i < 13; i++)
     {
         x_k_k[i] = 0;
     }
     
     for (j = 0; j < 13; j ++) {
         p_k_k[(13*j)+k] = 0;
     }
    double* initialize_x(double* x_k_k);
     x_k_k = initialize_x(x_k_k);
    double * initialize_p(double* p_k_k);
     p_k_k = initialize_p(p_k_k);*/
}

int ekf_filter (NSArray * varargin, ... ){
    int n = sizeof(varargin);
    if (n == 0) {
        //double ftype[][] = {};
        
    }
    return n;
    // Code for the measurement vector and other vector initialization goes here
}



double * EKF_PredictionState (NSArray * xk_kcap, NSArray * pk_kcap) {
    
    //Using the initial Velocity values to set the paramaters for State and Covariance estimate
    
    //double x_k_k[1][13] = {0, 0, 0, 1, 0, 0, 0, v_0, v_0, v_0, w_0, w_0, w_0};
    //double p_k_k[13][13] = {};
    //p_k_k[0][0] = 2^-52;
    //p_k_k[1][1] = 2^-52;
    //p_k_k[2][2] = 2^-52;
    //p_k_k[3][3] = 2^-52;
    //p_k_k[4][4] = 2^-52;
    //p_k_k[5][5] = 2^-52;
    //p_k_k[6][6] = 2^-52;
    //p_k_k[7][7] = std_v_o*std_v_o;
    //p_k_k[8][8] = std_v_o*std_v_o;
    //p_k_k[9][9] = std_v_o*std_v_o;
    //p_k_k[10][10] = std_w_o*std_w_o;
    //p_k_k[11][11] = std_w_o*std_w_o;
    //p_k_k[12][12] = std_w_o*std_w_o;
    
    
    
    // Initialize the values for noise parameters
    double sigma_a = 0.007;
    double sigma_alpha = 0.007;
    double sigma_image_noise = 1;
    
    // For State Transition Equation Derivatives
    
    //void dfv_by_dxv();
    
    
    double xk;
    //double pk;
    //return x_k_k; ERROR
    return NULL;
    // Code to predict the State vector of the EKF
};

 /*double  EKF_Measurement() {
    //Code to return the vector of measurements of camera and other parameters
    // Implementation of 'fv' function call in MATLAB code
    double rW [1][3] = {};
    double qWR [1][4] = {};
    double vW[1][3] = {};
    double wW[1][3] = {};
    rW [1][0] = x_k_k[1][0];
    rW [1][1] = x_k_k[1][1];
    rW [1][2] = x_k_k[1][2];
    qWR[1][0] = x_k_k[1][3];
    qWR[1][1] = x_k_k[1][4];
    qWR[1][2] = x_k_k[1][5];
    qWR[1][3] = x_k_k[1][6];
    vW [1][0] = x_k_k[1][7];
    vW [1][1] = x_k_k[1][8];
    vW [1][2] = x_k_k[1][9];
    wW [1][0] = x_k_k[1][10];
    wW [1][1] = x_k_k[1][11];
    wW [1][2] = x_k_k[1][12];
    
    if ((type = @"constant_orientation"))
    {
        wW[1][0] = 0;
        wW[1][1] = 0;
        wW[1][2] = 0;
        x_k_km1[1][0] = rW[1][2]+vW[1][2];
        x_k_km1[2][0] = qWR[1][3];
        x_k_km1[3][0] = vW[1][2];
        x_k_km1[4][0] = wW[1][2];
        
        
    }
    else if ((type = @"constant_position"))
    {
        vW[1][0] = 0;
        vW[1][1] = 0;
        vW[1][2] = 0;
        
    }
    
    else if ((type = @"constant_position_and_orientation"))
    {
        wW[1][0] = 0;
        wW[1][1] = 0;
        wW[1][2] = 0;
        vW[1][0] = 0;
        vW[1][1] = 0;
        vW[1][2] = 0;
        
        x_k_km1[1][0] = rW[1][2];
        x_k_km1[2][0] = qWR[1][3];
        x_k_km1[3][0] = vW[1][2];
        x_k_km1[4][0] = wW[1][2];
        
        
    }
    
    else if ((type = @"constant_position_and_orientation_location_noise"))
    {
        wW[1][0] = 0;
        wW[1][1] = 0;
        wW[1][2] = 0;
        vW[1][0] = 0;
        vW[1][1] = 0;
        vW[1][2] = 0;
        
        x_k_km1[1][0] = rW[1][2];
        x_k_km1[2][0] = qWR[1][3];
        x_k_km1[3][0] = vW[1][2];
        x_k_km1[4][0] = wW[1][2];
        
        
    }
    
    else
    {
        // CODE FOR COSTANT VELOCITY TYPE GOES HERE
        x_k_km1[3][0] = vW[1][2];
        x_k_km1[4][0] = wW[1][2];
        
    }
    
    double z;*/
     
    // double* func_Q(double * x_k_k[1][13],double * zeroessixbyone[6][1], double  pn[6][6],double delta_t, //char* type, double* Q[13][13]);
    //return x_k_km1;
     double p;
     /**** DUPLICATE CODE OF func_Q FUNCTION. NEEDS TO BE REMOVED AT THE END AFTER CHECKS
     /*double* omegaOld[1][3] = {};
     double* qOld[1][4] = {};

     if((type=@"constant_position_and_orientation_location_noise")) {
         *omegaOld[0][0] = x_k_k[0][10];
         *omegaOld[0][1] = x_k_k[0][11];
         *omegaOld[0][2] = x_k_k[0][12];
         *qOld[0][0] = x_k_k[0][3];
         *qOld[0][1] = x_k_k[0][4];
         *qOld[0][2] = x_k_k[0][5];
         *qOld[0][3] = x_k_k[0][6];
         
         //Equivalent of declaration of sparse(zeros(13,6)) in MATLAB code
         double *G[13][6] = {};
         *G[0][0] = 1*delta_t;
         *G[0][1] = 0*delta_t;
         *G[0][2] = 0*delta_t;
         *G[1][0] = 0*delta_t;
         *G[1][1] = 1*delta_t;
         *G[1][2] = 0*delta_t;
         *G[2][0] = 0*delta_t;
         *G[2][1] = 0*delta_t;
         *G[2][2] = 1*delta_t;
         
         //Calling the q2tr method passing qOld matrix as the input and obtaining t which is a 4 x 4 matrix as the output
         
         double *t[4][4] = {};
         double *rpy[1][4];
         double *d[4][3];
         double *q2tr(double *qOld[1][4], double *t[4][4]);
         double *tr2rpy(double* t[4][4], double*rpy[1][3]);
         double *dq_by_deuler(double *rpy[1][3], double *d[4][3]);
         
         *G[3][3] = *d[0][0];
         *G[3][4] = *d[0][1];
         *G[3][5] = *d[0][2];
         *G[4][3] = *d[1][0];
         *G[4][4] = *d[1][1];
         *G[4][5] = *d[1][2];
         *G[5][3] = *d[2][0];
         *G[5][4] = *d[2][1];
         *G[5][5] = *d[2][2];
         *G[6][3] = *d[3][0];
         *G[6][4] = *d[3][1];
         *G[6][5] = *d[3][2];
         
         
     }
     
     else {
         *omegaOld[0][0] = x_k_k[0][10];
         *omegaOld[0][1] = x_k_k[0][11];
         *omegaOld[0][2] = x_k_k[0][12];
         *qOld[0][0] = x_k_k[0][3];
         *qOld[0][1] = x_k_k[0][4];
         *qOld[0][2] = x_k_k[0][5];
         *qOld[0][3] = x_k_k[0][6];
         
         double *qwt[1][4] = {};
         double* v2q(double *omegaOld[1][3], double *qwt[1][4]);
         
          double *G[13][6] = {};
          *G[7][0] = 1;
         *G[7][1] = 0;
         *G[7][2] = 0;
         *G[8][0] = 0;
         *G[8][1] = 1;
         *G[8][2] = 0;
         *G[9][0] = 0;
         *G[9][1] = 0;
         *G[9][2] = 1;
         
         *G[10][3] = 1;
         *G[10][4] = 0;
         *G[10][5] = 0;
         *G[11][3] = 0;
         *G[11][4] = 1;
         *G[11][5] = 0;
         *G[12][3] = 0;
         *G[12][4] = 0;
         *G[12][5] = 1;
         
         *G[0][0] = 1*delta_t;
         *G[0][1] = 0*delta_t;
         *G[0][2] = 0*delta_t;
         *G[1][0] = 0*delta_t;
         *G[1][1] = 1*delta_t;
         *G[1][2] = 0*delta_t;
         *G[2][0] = 0*delta_t;
         *G[2][1] = 0*delta_t;
         *G[2][2] = 1*delta_t;
         
         double* dq3_by_dq1RES[4][4] = {};
         double *dqomegadt_by_domegaRES[4][3] = {};
         
         //Need to write the functions dq3_by_dq1(qOld) and qomegadt_by_domega(omegaOld,dt)
         double * dq3_by_dq1(double *qOld[1][4], double* dq3_by_dq1RES[4][4]);
         double * dqomegadt_by_domega(double *omegaOld[1][3],double delta_t, double *dqomegadt_by_domegaRES[4][3]);
         //The Matrix multiplication is given below. The multiplied result is stored in rows 4 through 7 and columns 10 through 12. It is the same as in the MATLAB file.
         
         //NEED TO EDIT THE MATRIX INDICES TO PRODUCE THE CORRECT VALUES UPON MULTIPLICATION
         
         *G[3][3] = *dq3_by_dq1RES[0][0]**dqomegadt_by_domegaRES[0][0] + *dq3_by_dq1RES[1][1]**dqomegadt_by_domegaRES[2][0] + *dq3_by_dq1RES[1][2]**dqomegadt_by_domegaRES[3][0] + *dq3_by_dq1RES[1][3]**dqomegadt_by_domegaRES[4][0];
         *G[3][4] = *dq3_by_dq1RES[1][0]**dqomegadt_by_domegaRES[1][1] + *dq3_by_dq1RES[1][1]**dqomegadt_by_domegaRES[2][1] + *dq3_by_dq1RES[1][2]**dqomegadt_by_domegaRES[3][1] + *dq3_by_dq1RES[1][3]**dqomegadt_by_domegaRES[4][1];
         *G[3][5] = *dq3_by_dq1RES[1][0]**dqomegadt_by_domegaRES[1][2] + *dq3_by_dq1RES[1][1]**dqomegadt_by_domegaRES[2][2] + *dq3_by_dq1RES[1][2]**dqomegadt_by_domegaRES[3][2] + *dq3_by_dq1RES[1][3]**dqomegadt_by_domegaRES[4][2];
         *G[4][3] = *dq3_by_dq1RES[2][0]**dqomegadt_by_domegaRES[1][0] + *dq3_by_dq1RES[2][1]**dqomegadt_by_domegaRES[2][0] + *dq3_by_dq1RES[2][2]**dqomegadt_by_domegaRES[3][0] + *dq3_by_dq1RES[2][3]**dqomegadt_by_domegaRES[4][0];
         *G[4][4] = *dq3_by_dq1RES[2][0]**dqomegadt_by_domegaRES[1][1] + *dq3_by_dq1RES[2][1]**dqomegadt_by_domegaRES[2][1] + *dq3_by_dq1RES[2][2]**dqomegadt_by_domegaRES[3][1] + *dq3_by_dq1RES[2][3]**dqomegadt_by_domegaRES[4][1];
         *G[4][5] = *dq3_by_dq1RES[2][0]**dqomegadt_by_domegaRES[1][2] + *dq3_by_dq1RES[2][1]**dqomegadt_by_domegaRES[2][2] + *dq3_by_dq1RES[2][2]**dqomegadt_by_domegaRES[3][2] + *dq3_by_dq1RES[2][3]**dqomegadt_by_domegaRES[4][2];
         *G[5][3] = *dq3_by_dq1RES[3][0]**dqomegadt_by_domegaRES[1][0] + *dq3_by_dq1RES[3][1]**dqomegadt_by_domegaRES[2][0] + *dq3_by_dq1RES[3][2]**dqomegadt_by_domegaRES[3][0] + *dq3_by_dq1RES[3][3]**dqomegadt_by_domegaRES[4][0];
         *G[5][4] = *dq3_by_dq1RES[3][0]**dqomegadt_by_domegaRES[1][1] + *dq3_by_dq1RES[3][1]**dqomegadt_by_domegaRES[2][1] + *dq3_by_dq1RES[3][2]**dqomegadt_by_domegaRES[3][1] + *dq3_by_dq1RES[3][3]**dqomegadt_by_domegaRES[4][1];
         *G[5][5] = *dq3_by_dq1RES[3][0]**dqomegadt_by_domegaRES[1][2] + *dq3_by_dq1RES[3][1]**dqomegadt_by_domegaRES[2][2] + *dq3_by_dq1RES[3][2]**dqomegadt_by_domegaRES[3][2] + *dq3_by_dq1RES[3][3]**dqomegadt_by_domegaRES[4][2];
         *G[6][3] = *dq3_by_dq1RES[4][0]**dqomegadt_by_domegaRES[1][0] + *dq3_by_dq1RES[4][1]**dqomegadt_by_domegaRES[2][0] + *dq3_by_dq1RES[4][2]**dqomegadt_by_domegaRES[3][0] + *dq3_by_dq1RES[4][3]**dqomegadt_by_domegaRES[4][0];
         *G[6][4] = *dq3_by_dq1RES[4][0]**dqomegadt_by_domegaRES[1][1] + *dq3_by_dq1RES[4][1]**dqomegadt_by_domegaRES[2][1] + *dq3_by_dq1RES[4][2]**dqomegadt_by_domegaRES[3][1] + *dq3_by_dq1RES[4][3]**dqomegadt_by_domegaRES[4][1];
         *G[6][5] = *dq3_by_dq1RES[4][0]**dqomegadt_by_domegaRES[1][2] + *dq3_by_dq1RES[4][1]**dqomegadt_by_domegaRES[2][2] + *dq3_by_dq1RES[4][2]**dqomegadt_by_domegaRES[3][2] + *dq3_by_dq1RES[4][3]**dqomegadt_by_domegaRES[4][2];
        //Declaring and Initializing G Transpose Matrix as the transpose of G Matrix
         double *G_Transpose[6][13];
         
         int m,n;
         for(m = 0; m <=5; m++) {
             for(n = 0; n <=12; n++) {
                  G_Transpose[m][n] = G[n][m];
             }
         
        }
         
         //G
         // The Matrix Multiplication for Q goes here
         
        // return Q;
         
         //Declaration for Q_int1 which results in the multiplication of the matrices G[13][6] and pn[6][6]
         
        double* Q_int[13][6] = {};
        
         //Matrix Multiplication code goes here for Q_int1[13][6]
         //Q_int[0][0] = *G[0][0]**pn[0][0] + *G[0][1]**pn[1][0] + *G[0][2]**pn[2][0];
         
         
         
         
     }*/
     
     
    //return p;
 //};

 double EKF_getmeasurementvector()  {
    // Method to return the Measurement Vector
    double z_ic;
    return z_ic;
};

 double EKF_StateUpdateMeasurement(double * z_i,double *xk_kcap) {
    // Method to update the State Value of the EKF
    double x_i;
    return x_i;
};

 double EKF_Update() {
    // Updates the EKF with the Low and High Innovation Inliers at the end
    double x;
    //double p;
    
    //Code from the MATLAB function ekf_update_li_inliers
    void update();
    return x;
    
    //Code from the MATLAB function ekf_update_hi_inliers
    void update();
};







void dq0_by_domegaA() {
    
}













//@end