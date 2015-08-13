//
//  MapManagement.m
//  drySLAM
//
//  Created by Anthony Rodriguez on 4/17/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//

#include "EKFPrediction.hpp"

//Main gateway:
bool ekf_prediction(Filter* filter, NSMutableArray* features_info){
    double* x_temp = filter->x_k_km1;
    double* p_temp = filter->p_k_km1;
    double* x_new = (double*)malloc(filter->state_size*sizeof(double));
    double* p_new = (double*)malloc(filter->state_size*filter->state_size*sizeof(double));
    bool success = predict_state_and_covariance(filter->x_k_k, filter->p_k_k, filter->state_size, x_new, p_new, filter->type, filter->std_a, filter->std_alpha);
    filter->x_k_km1 = x_new;
    filter->p_k_km1 = p_new;
    free(x_temp);
    free(p_temp);
    return success;
}

/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
//INTERNAL FUNCTIONS:

// Function dq3_by_dq2 used for forming the matrix q1.R which is a 1 x 4 matrix. Have used q1_R instead of q1.R etc as semantics do not allow such a declaration
double *dq3_by_dq2(double* qwt, double*dfv_by_dxvRES) {
    // Method to calculate the Jacobians
    double q1_R, q1_X, q1_Y, q1_Z, minusq1_X, minusq1_Y, minusq1_Z; //Equivalent of q1.R, q1.X, q1.Y, q1.Z in MATLAB code
    q1_R = qwt[0];
    q1_X = qwt[1];
    q1_Y = qwt[2];
    q1_Z = qwt[3];
    minusq1_X = q1_X*(-1);
    minusq1_Y = q1_Y*(-1);
    minusq1_Z = q1_Z*(-1);
    
    //Updating the values to the dfv_by_dxvRES matrix for the rows and columns 4, 5, 6 and 7
    dfv_by_dxvRES[42] = q1_R;
    dfv_by_dxvRES[43] = minusq1_X;
    dfv_by_dxvRES[44] = minusq1_Y;
    dfv_by_dxvRES[45] = minusq1_Z;
    dfv_by_dxvRES[55] = q1_X;
    dfv_by_dxvRES[56] = q1_R;
    dfv_by_dxvRES[57] = q1_Z;
    dfv_by_dxvRES[58] = minusq1_Y;
    dfv_by_dxvRES[68] = q1_Y;
    dfv_by_dxvRES[69] = minusq1_Z;
    dfv_by_dxvRES[70] = q1_R;
    dfv_by_dxvRES[71] = q1_X;
    dfv_by_dxvRES[81] = q1_Z;
    dfv_by_dxvRES[82] = q1_Y;
    dfv_by_dxvRES[83] = minusq1_X;
    dfv_by_dxvRES[84] = q1_R;
    
    
    //return fv_by_dxvRES; build ERROR
    return dfv_by_dxvRES;
}

double * dq3_by_dq1(double *q2_in, double* dq3_by_dq1RES) {
    double q2_R, q2_X, q2_Y, q2_Z, minusq2_X, minusq2_Y, minusq2_Z;
    q2_R = q2_in[0];
    q2_X = q2_in[1];
    q2_Y = q2_in[2];
    q2_Z = q2_in[3];
    
    minusq2_X = q2_X*(-1);
    minusq2_Y = q2_Y*(-1);
    minusq2_Z = q2_Z*(-1);
    
    dq3_by_dq1RES[0] = q2_R;
    dq3_by_dq1RES[1] = minusq2_X;
    dq3_by_dq1RES[2] = minusq2_Y;
    dq3_by_dq1RES[3] = minusq2_Z;
    dq3_by_dq1RES[4] = q2_X;
    dq3_by_dq1RES[5] = q2_R;
    dq3_by_dq1RES[6] = minusq2_Z;
    dq3_by_dq1RES[7] = q2_Y;
    dq3_by_dq1RES[8] = q2_Y;
    dq3_by_dq1RES[9] = q2_Z;
    dq3_by_dq1RES[10] = q2_R;
    dq3_by_dq1RES[11] = minusq2_X;
    dq3_by_dq1RES[12] = q2_Z;
    dq3_by_dq1RES[13] = minusq2_Y;
    dq3_by_dq1RES[14] = q2_X;
    dq3_by_dq1RES[15] = q2_R;
    
    return dq3_by_dq1RES;//ERROR
    // return NULL;
}

double * initialize_x(double* x_k_k) {
    //double x_k_k[0] = {0, 0, 0, 1, 0, 0, 0, v_0, v_0, v_0, w_0, w_0, w_0};
    //double x_k_k[1][13];
    //*x_k_k[0] = 0;
    
    x_k_k[0] = x_k_k[1] = x_k_k[2] = x_k_k[4] = x_k_k[5] = x_k_k[6] = 0;
    x_k_k[3] = 1;
    x_k_k[7] = x_k_k[8] = x_k_k[9] = 0; // Since v_0 is zero
    x_k_k[10] = x_k_k[11] = x_k_k[12] = (10)^-15; //Since 1e-15 is 10^-15
    return x_k_k;
}

double * initialize_p(double* p_k_k) {
    double std_v_o = 0.025;
    double std_w_o = 0.025;
    
    // Modified code in order to better use the Matrices in C
    
    p_k_k[0]=-2^-52;
    p_k_k[14] = 2^-52;
    p_k_k[28] = 2^-52;
    p_k_k[42] = 2^-52;
    p_k_k[56] = 2^-52;
    p_k_k[70] = 2^-52;
    p_k_k[84] = 2^-52;
    p_k_k[98] = std_v_o*std_v_o;
    p_k_k[112] = std_v_o*std_v_o;
    p_k_k[126] = std_v_o*std_v_o;
    p_k_k[140] = std_w_o*std_w_o;
    p_k_k[154] = std_w_o*std_w_o;
    p_k_k[168] = std_w_o*std_w_o;
    //End of the modified code for the function
    
    return p_k_k;
}


double * dq_by_deuler (double* tr2rpyRES, double* dq_by_deulerRES) {
    //double *G;
    // Code for the function begins here
    double phi = tr2rpyRES[0];
    double theta = tr2rpyRES[1];
    double psi = tr2rpyRES[2];
    
    dq_by_deulerRES[0] = (0.5)*((-1)*sin(phi/2) +cos(phi/2));
    dq_by_deulerRES[1] = (0.5)*((-1)*sin(theta/2) +cos(theta/2));
    dq_by_deulerRES[2] = (0.5)*((-1)*sin(psi/2) +cos(psi/2));
    dq_by_deulerRES[3] = (0.5)*(cos(phi/2) + sin(phi/2));
    dq_by_deulerRES[4] = (0.5)*((-1)*sin(theta/2) -cos(theta/2));
    dq_by_deulerRES[5] = (0.5)*((-1)*sin(psi/2) -cos(psi/2));
    dq_by_deulerRES[6] = (0.5)*((-1)*sin(phi/2) +cos(phi/2));
    dq_by_deulerRES[7] = (0.5)*((-1)*sin(theta/2) +cos(theta/2));
    dq_by_deulerRES[8] = (0.5)*((-1)*sin(psi/2) +cos(psi/2));
    dq_by_deulerRES[9] = (0.5)*((-1)*sin(phi/2) -cos(phi/2));
    dq_by_deulerRES[10] = (0.5)*((-1)*sin(theta/2) -cos(theta/2));
    dq_by_deulerRES[11] = (0.5)*(cos(psi/2) + sin(psi/2));
    
    
    
    return dq_by_deulerRES;
}
// Implementation of the v2q function as in MATLAB
double * v2q(double omega1, double omega2, double omega3, double* q) {
    //double value1 = omega;
    double squared_value1 = omega1*omega1;
    // double value2 = v[1];
    double squared_value2 = omega2*omega2;
    //double value3 = v[2];
    double squared_value3 = omega3*omega3;
    //double *q = (double *)malloc(1*4*sizeof(double));
    
    // Normalized function in MATLAB code in v2q.m file is implemented here. Theta stores the normalized value
    double theta = sqrt(squared_value1 + squared_value2 + squared_value3);
    // eps function in MATLAB evaluates to 2^-52
    if (theta < ((2)^-52)) {
        //double *q[1][4] = {1, 0, 0, 0};
        
        q[0] = 1;
        q[1] = q[2] = q[3] = 0;
    }
    else {
        //double *v_n[1][3] =  {};
        double *v_n = (double *)malloc(1*3*sizeof(double));
        v_n[0] = omega1/theta;
        v_n[1] = omega2/theta;
        v_n[2] = omega3/theta;
        double *quaternion(double *v_n, double theta, double *q);
        q = quaternion(v_n, theta, q);
        
    }
    return q; //ERROR
    //return NULL;
}

double *quaternion(double *v_n, double theta, double *q) {
    //return *q; // Returns the 1 x 4 q array which takes into account the camera movements (ERROR)
    double value1 = v_n[0];
    double squared_value1 = value1*value1;
    double value2 = v_n[1];
    double squared_value2 = value2*value2;
    double value3 = v_n[2];
    double squared_value3 = value3*value3;
    
    double norm_v_n = sqrt(squared_value1 + squared_value2 + squared_value3);
    // Assigning values to the q matrix to form the quaternion
    q[0] = cos(theta/2);
    q[1] = sin(theta/2)*value1/norm_v_n;
    q[2] = sin(theta/2)*value2/norm_v_n;
    q[3] = sin(theta/2)*value3/norm_v_n;
    return q;
}

double *q2tr(double *qOld, double *t) {
    double s = qOld[0];
    double x = qOld[1];
    double y = qOld[2];
    double z = qOld[3];
    
    //double *r[3][3] = {};
    
    // Declaring the matrix r which is a 3 x 3 matrix and setting the values as specified in the MATLAB code
    
    double* r = (double *)malloc(3*3*sizeof(double *));
    
    r[0] = 1 -2*((y*y) + (z*z));
    r[1] = 2*((x*y) - (s*z));
    r[2] = 2*((x*z) + (s*y));
    r[3] = 2*((x*y) + (s*z));
    r[4] = 1 -2*((x*x) + (z*z));
    r[5] = 2*((y*z) - (s*x));
    r[6] = 2*((x*z) - (s*y));
    r[7] = 2*((y*z) + (s*x));
    r[8] = 1 -2*((x*x) + (y*y));
    
    // Setting the diagonal values of the t matrix to be 1, making it an identity matrix
    
    
    
    t[0] = 1;
    t[5] = 1;
    t[10] = 1;
    t[15] = 1;
    
    // Adding values to the t matrix from the r matrix and setting the last value in the matrix as 1.
    
    
    /* *t[0][0] = *r[0][0];
     *t[0][1] = *r[0][1];
     *t[0][2] = *r[0][2];
     *t[1][0] = *r[1][0];
     *t[1][1] = *r[1][1];
     *t[1][2] = *r[1][2];
     *t[2][0] = *r[2][0];
     *t[2][1] = *r[2][1];
     *t[2][2] = *r[2][2];
     
     *t[3][3] = 1;*/
    
    t[0] = r[0];
    t[1] = r[1];
    t[2] = r[2];
    t[4] = r[3];
    t[5] = r[4];
    t[6] = r[5];
    t[8] = r[6];
    t[9] = r[7];
    t[10] = r[8];
    
    t[15] = 1;
    
    
    
    return t;
    
    
    
    
    
}

double *tr2rpy(double* q2trRES, double* tr2pryRES) {
    
    // Taking the negative values first in order to use them in the if condition below
    
    double neg_r2c0 = (-1)*q2trRES[8];
    double neg_r1c2 = (-1)*q2trRES[6];
    
    if((abs((q2trRES[0] < (2^-52) ))) && (abs((q2trRES[4] < (2^-52) ))) ) {
        tr2pryRES[0] = 0;
        
        //Declaration of the negative MATRIX elements
        
        tr2pryRES[1] = atan2(neg_r2c0, tr2pryRES[0]);
        tr2pryRES[2] = atan2(neg_r1c2, tr2pryRES[5]);
    }
    
    else {
        tr2pryRES[0] = atan2(q2trRES[4], q2trRES[0]);
        double sp = sin(tr2pryRES[0]);
        double cp = cos(tr2pryRES[0]);
        tr2pryRES[1] = atan2(neg_r2c0, ((cp*tr2pryRES[0]) + (sp*tr2pryRES[4])));
        tr2pryRES[2] = atan2(((sp*tr2pryRES[2]) - (cp*tr2pryRES[6])), ((cp*tr2pryRES[5]) - sp*tr2pryRES[2]));
        
        
    }
    return tr2pryRES;
    
}

/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/

bool predict_state_and_covariance(double* x_k, double* p_k, int state_size, double* x_km1_k, double* p_km1_k, NSString* type, double SD_A_component_filter, double SD_alpha_component_filter) {
    //This is not a good place to declare a function...
    // double * fv (double *x_k_k[1][13], double delta_t, char* type, double SD_A_component_filter, double SD_alpha_component_filter);
    // double * x_km1_k[1][12];
    if (x_k==NULL || p_k==NULL || x_km1_k == NULL || p_km1_k==NULL || type==NULL)
        return FALSE;
    double delta_t = DELTA_T;
    double *xv = (double *) malloc(1*13*sizeof(double ));
    int i;
    for (i = 0; i <13; i++) {
        xv[i] = x_k[i];
    }
    
    /* Function declaration and other initialization for fv function */
    
    // For Camera Motion Prediction, we call the fv method passing the x_k matrix from the first till 13th value as an input.
    
    
    // Matrix Declaration for Xv_km1_k which is a 1 x 13 array
    
    double *Xv_km1_k = (double *)malloc(13*sizeof(double ));
    // Matrix Initialization for Xv_km1_k
    int e;
    for (e = 0; e < 13; e++) {
        Xv_km1_k[e] = 0;
    }
    
    double* fv (double *xv_km1_k, double* xv, double delta_t, NSString* type, double std_a, double std_alphaa);
    
    Xv_km1_k = fv(Xv_km1_k, xv, delta_t, type, SD_A_component_filter, SD_alpha_component_filter);
    
    // Setting the values to the predicted features matrix which is x_km1_K. This code performs features prediction and is translated from the MATLAB code.
    
    // THE BELOW ITERATION RETURNS THE PREDICTED STATE MATRIX VALUES
    
    int g;
    for (g = 0; g < state_size; g++ )
    {
        if (g < 13) {
            x_km1_k[g] = Xv_km1_k[g];
        }
        else {
            x_km1_k[g] = x_k[g];
        }
    }
    
    
    
    
    
    
    
    
    // Declaration  of 6 x 1 all zero matrix for the dfv_by_dxv method function input argument
    
    double *dfv_by_dxv_zero_matrix = (double *)malloc(6*1*sizeof(double));
    int j;
    for (j = 0; j  < 6; j++) {
        dfv_by_dxv_zero_matrix[(1*j)+0] = 0;
    }
    
    // Declaration of the dfv_by_dxvRES matrix which is a 13 x 13 matrix to store the result of the function, dfv_by_dxv
    double *dfv_by_dxvRES = (double *) malloc(13*13*sizeof(double));
    int x,y;
    for (x = 0; x <13; x++) {
        for (y = 0; y <13; y++) {
            dfv_by_dxvRES[(13*x)+y] = 0;
        }
        
    }
    
    double *dfv_by_dxv(double *dfv_by_dxvRES, double *xv, double *dfv_by_dxv_zero_matrix, double delta_t, NSString* type);
    
    //dfv_by_dxvRES = dfv_by_dxv(dfv_by_dxvRES, xv, dfv_by_dxv_zero_matrix, delta_t, type);
    dfv_by_dxvRES = dfv_by_dxv(dfv_by_dxvRES, xv, dfv_by_dxv_zero_matrix, delta_t, type);
    
    
    // State Noise Calculation as calculated in the MATLAB file
    double linear_acceleration_noise_covariance = (SD_A_component_filter*delta_t)*(SD_A_component_filter*delta_t);
    double angular_acceleration_noise_covariance = (SD_alpha_component_filter*delta_t)*(SD_alpha_component_filter*delta_t);
    
    //return x_km1_k;
    //Declaration and initialization of the Matrix Pn. The matrix is initialized with all zeroes and the diagonal elements carry the values as shown below.
    
    //******This syntax causes a crash******** -AR
    //To declare a 2d array, either do double* Pn = (double*)malloc(6*6*sizeof(double));
    //                              OR double Pn[6][6];
    double* Pn = (double*)malloc(6*6*sizeof(double));
    for (int i=0; i<36; ++i){
        Pn[i] = 0.0;
    }
    Pn[0] = linear_acceleration_noise_covariance;
    Pn[7] = linear_acceleration_noise_covariance;
    Pn[14] = linear_acceleration_noise_covariance;
    Pn[21] = angular_acceleration_noise_covariance;
    Pn[28] = angular_acceleration_noise_covariance;
    Pn[35] = angular_acceleration_noise_covariance;
    
    // Declaring and initializing the Q matrix as a matrix of all zeros which is passed to the func_Q function to update
    double* Q = (double*)malloc(13*13*sizeof(double));
    
    int v,s;
    for (v = 0; v < 13; v++)
    {
        for (s = 0; s < 13; s++) {
            Q[(13*v)+s] = 0;
        }
    }
    //double *func_Q_zero_matrix[6][1] = {0};
    double *func_Q_zero_matrix = (double*)malloc(6*1*sizeof(double));
    int u;
    for (u = 0; u <=5; u++) {
        func_Q_zero_matrix[1*u] = 0;
    }
    double *func_Q (double *xv, double* func_Q_zero_matrix, double *Pn, double delta_t, NSString* type, double *Q);
    
    Q = func_Q(xv, func_Q_zero_matrix, Pn, delta_t, type, Q);
    
    //Matrix Declaration for F. The dfv_by_dxvRES is the F matrix declared in MATLAB code. The transpose is taken for the Covariance calculation.
    
    cv::Mat F_mat(13, 13, CV_64FC1, dfv_by_dxvRES, cv::Mat::AUTO_STEP);
    
    cv::Mat F_transpose;
    cv::transpose(F_mat, F_transpose);
    
    // The Q matrix is converted to an openCV matrix
    
    cv::Mat Q_mat(13, 13, CV_64FC1, Q, cv::Mat::AUTO_STEP);
    
    // The conversion of the p_k matrix to an openCV matrix
    
    cv::Mat p_k_mat(state_size, state_size, CV_64FC1, p_k, cv::Mat::AUTO_STEP);
    
    int remainder = state_size-13;
    
    double* p_k_13By13 = (double*)malloc(13*13*sizeof(double));
    for (int i=0; i<13; ++i){
        double* ptr1 = (p_k + state_size*i);
        double* ptr2 = (p_k_13By13 + 13*i);
        for (int j=0; j<13; ++j){
            ptr2[j] = ptr1[j];
        }
    }
    cv::Mat p_k_13By13_mat(13, 13, CV_64FC1, p_k_13By13, cv::Mat::AUTO_STEP);
    
    double* p_k_13ByRemainder = (double*)malloc(13*remainder*sizeof(double));
    for (int i=0; i<13; ++i){
        double* ptr1 = (p_k+state_size*i + 13);
        double* ptr2 = (p_k_13ByRemainder + remainder*i);
        for (int j=0; j<remainder; ++j){
            ptr2[j] = ptr1[j];
        }
    }
    cv::Mat p_k_13ByRemainder_mat(13, remainder, CV_64FC1, p_k_13ByRemainder, cv::Mat::AUTO_STEP);
    
    double* p_k_RemainderBy13 = (double*)malloc(13*remainder*sizeof(double));
    for (int i=0; i<remainder; ++i){
        double* ptr1 = (p_k + state_size*(13+i));
        double* ptr2 = (p_k_RemainderBy13 + 13*i);
        for (int j=0; j<13; ++j){
            ptr2[j] = ptr1[j];
        }
    }
    cv::Mat p_k_RemainderBy13_mat(remainder, 13, CV_64FC1, p_k_RemainderBy13, cv::Mat::AUTO_STEP);
    
    
    double* p_k_RemainderByRemainder = (double*)malloc(remainder*remainder*sizeof(double));
    for (int i=0; i<remainder; ++i){
        double* ptr1 = (p_k + state_size*(13+i) + 13);
        double* ptr2 = (p_k_RemainderByRemainder + remainder*i);
        for (int j=0; j<remainder; ++j){
            ptr2[j] = ptr1[j];
        }
    }
    cv::Mat p_k_RemainderByRemainder_mat(remainder, remainder, CV_64FC1, p_k_RemainderByRemainder, cv::Mat::AUTO_STEP);
    
    cv::Mat product1 = F_mat * p_k_13By13_mat * F_transpose + Q_mat;
    cv::Mat product2 = F_mat * p_k_13ByRemainder_mat;
    cv::Mat product3 = p_k_RemainderBy13_mat * F_transpose;
    cv::Mat product4 = p_k_RemainderByRemainder_mat;
    
    //reassemble into p_km1_k
    for (int i=0; i<13; ++i){
        double* ptrL = product1.ptr<double>(i);
        double* ptrR = product2.ptr<double>(i);
        double* p_k_rowPtr = p_km1_k + state_size*i;
        for (int j=0; j<13; ++j){
            p_k_rowPtr[j] = ptrL[j];
        }
        for (int j=13; j<state_size; ++j){
            p_k_rowPtr[j] = ptrR[j-13];
        }
    }
    for (int i=13; i<state_size; ++i){
        double* ptrL = product3.ptr<double>(i-13);
        double* ptrR = product4.ptr<double>(i-13);
        double* p_k_rowPtr = p_km1_k + state_size*i;
        for (int j=0; j<13; ++j){
            p_k_rowPtr[j] = ptrL[j];
        }
        for (int j=13; j<state_size; ++j){
            p_k_rowPtr[j] = ptrR[j-13];
        }
    }
    /* End of Matrix Conversion code */
    // This is implemented in MATLAB code, but we do not need the size of the array as state size is one of the inputs to the function.
    //Note: I'm pretty sure the above doesn't work, but the value "state_size" should be what you need here^
    
    //int size_P_k = sizeof(p_k[0]);
    
    // The code for the formation of P_km1_k matrix goes here
    
    free(p_k_13By13);
    free(p_k_13ByRemainder);
    free(p_k_RemainderByRemainder);
    free(p_k_RemainderBy13);
    product1.release();
    product2.release();
    product3.release();
    product4.release();
    
    free(Xv_km1_k);
    free(dfv_by_dxvRES);
    free(xv);
    free(func_Q_zero_matrix);
    free(Q);
    
    return true;
}


/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/


// Method called when the camera motion is to be estimated
double* fv (double *xv_km1_k, double* xv, double delta_t, NSString* type, double std_a, double std_alphaa) {
    //double *xkkm1[13][1] = {};
    //double* x_k_km1[4][4] = {};

    // Actual code as observed in MATLAB
    double * rW = (double *)malloc(3*sizeof(double));
    double * qWR = (double *)malloc(4*sizeof(double));
    double * vW = (double *)malloc(3*sizeof(double));
    double * wW = (double *)malloc(3*sizeof(double));
    
    rW[0] = xv[0];
    rW[1] = xv[1];
    rW[2] = xv[2];
    
    qWR[0] = xv[3];
    qWR[1] = xv[4];
    qWR[2] = xv[5];
    qWR[3] = xv[6];
    
    vW[0] = xv[7];
    vW[1] = xv[8];
    vW[2] = xv[9];
    
    wW[0] = xv[10];
    wW[1] = xv[11];
    wW[2] = xv[12];
    /*double* rW[1][3] = {x_k_k[1][0], x_k_k[1][1], x_k_k[1][2]};
     double* qWR[1][4] = {x_k_k[1][3], x_k_k[1][4], x_k_k[1][5], x_k_k[1][6]};
     double* vW[1][3] = {x_k_k[1][7], x_k_k[1][8], x_k_k[1][9]};
     double* wW[1][3] = {x_k_k[1][10], x_k_k[1][11], x_k_k[1][12]};*/
    
    if([type isEqualToString: @"constant_orientation"]) {
        wW[0] = 0;
        wW[1] = 0;
        wW[2] = 0;
        /**x_k_km1[1][0] = *rW[1][0]*delta_t + *vW[1][0]*delta_t;
         *x_k_km1[1][1] = *rW[1][1]*delta_t + *vW[1][1]*delta_t;
         *x_k_km1[1][2] = *rW[1][2]*delta_t + *vW[1][2]*delta_t;
         *x_k_km1[2][0] = *qWR[1][0];
         *x_k_km1[2][1] = *qWR[1][1];
         *x_k_km1[2][2] = *qWR[1][2];
         *x_k_km1[2][3] = *qWR[1][3];
         *x_k_km1[3][0] = *vW[1][0];
         *x_k_km1[3][1] = *vW[1][1];
         *x_k_km1[3][2] = *vW[1][2];
         *x_k_km1[4][0] = *wW[1][0];
         *x_k_km1[4][1] = *wW[1][1];
         *x_k_km1[4][2] = *wW[1][2];*/
        
        xv_km1_k[0] = rW[0] + (vW[0]*delta_t);
        xv_km1_k[1] = rW[1] + (vW[1]*delta_t);
        xv_km1_k[2] = rW[2] + (vW[2]*delta_t);
        xv_km1_k[3] = qWR[0];
        xv_km1_k[4] = qWR[1];
        xv_km1_k[5] = qWR[2];
        xv_km1_k[6] = qWR[3];
        xv_km1_k[7] = vW[0];
        xv_km1_k[8] = vW[1];
        xv_km1_k[9] = vW[2];
        xv_km1_k[10] = wW[0];
        xv_km1_k[11] = wW[1];
        xv_km1_k[12] = wW[2];
        
        
        
        
    }
    
    if([type isEqualToString: @"constant_position"]) {
        vW[0] = 0;
        vW[1] = 0;
        vW[2] = 0;
        xv_km1_k[0] = rW[0];
        xv_km1_k[1] = rW[1];
        xv_km1_k[2] = rW[2];
        //double *input_arg_for_v2q[1][3] = *wW[1][3]*delta_t;
        
        //Method call for v2q function as given in MATLAB code. Passing the value of wW Matrix as wW*delta_1 = wW*1 = wW. Hence, passing the array directly
        //double *q[1][4] = {};
        //double *v2q(double *wW[1][3], double *q[1][4]);
        
        double * v2qRES = (double *)malloc(4*sizeof(double ));
        int i;
        for (i = 0; i < 4; i++) {
            v2qRES[i] = 0;
        }
        
        double wW_mod1 = wW[0]*delta_t;
        double wW_mod2 = wW[1]*delta_t;
        double wW_mod3 = wW[2]*delta_t;
        
        v2qRES = v2q(wW_mod1, wW_mod2, wW_mod3, v2qRES);
        
        
        
        double* qprodRES = (double *)malloc(4*sizeof(double));
        qprodRES[0] = qprodRES[1] = qprodRES[2] = qprodRES[3] = 0;
        
        
        qprodRES = qprod(qWR, v2qRES, qprodRES); // Returns a 1 x 4 array qprodRES as the output
        
        xv_km1_k[3] = qprodRES[0];
        xv_km1_k[4] = qprodRES[1];
        xv_km1_k[5] = qprodRES[2];
        xv_km1_k[6] = qprodRES[3];
        
        
        // The code for assigning to the xv_km1_k matrix goes here
        
        xv_km1_k[7] = vW[0];
        xv_km1_k[8] = vW[1];
        xv_km1_k[9] = vW[2];
        xv_km1_k[10] = wW[0];
        xv_km1_k[11] = wW[1];
        xv_km1_k[12] = wW[2];
        
    }
    
    if([type isEqualToString: @"constant_position_and_orientation"] || [type isEqualToString: @"constant_position_and_orientation_location_noise"]) {
        vW[0] = 0;
        vW[1] = 0;
        vW[2] = 0;
        wW[0] = 0;
        wW[1] = 0;
        wW[2] = 0;
        xv_km1_k[0] = rW[0];
        xv_km1_k[1] = rW[1];
        xv_km1_k[2] = rW[2];
        xv_km1_k[3] = qWR[0];
        xv_km1_k[4] = qWR[1];
        xv_km1_k[5] = qWR[2];
        xv_km1_k[6] = qWR[3];
        xv_km1_k[7] = vW[0];
        xv_km1_k[8] = vW[1];
        xv_km1_k[9] = vW[2];
        xv_km1_k[10] = wW[0];
        xv_km1_k[11] = wW[1];
        xv_km1_k[12] = wW[2];
    }
    
    if([type isEqualToString: @"constant_velocity"]) {
        xv_km1_k[0] = (rW[0] + (vW[0]*delta_t));
        xv_km1_k[1] = (rW[1] + (vW[1]*delta_t));
        xv_km1_k[2] = (rW[2] + (vW[2]*delta_t));
        
        double * v2qRES = (double *)malloc(4*sizeof(double ));
        int i;
        for (i = 0; i < 4; i++) {
            v2qRES[i] = 0;
        }
        
        double wW_mod1 = wW[0]*delta_t;
        double wW_mod2 = wW[1]*delta_t;
        double wW_mod3 = wW[2]*delta_t;
        
        v2qRES = v2q(wW_mod1, wW_mod2, wW_mod3, v2qRES);
        
        
        
        double* qprodRES = (double *)malloc(4*sizeof(double));
        qprodRES[0] = qprodRES[1] = qprodRES[2] = qprodRES[3] = 0;
        
        
        qprodRES = qprod(qWR, v2qRES, qprodRES); // Returns a 1 x 4 array qprodRES as the output
        
        xv_km1_k[3] = qprodRES[0];
        xv_km1_k[4] = qprodRES[1];
        xv_km1_k[5] = qprodRES[2];
        xv_km1_k[6] = qprodRES[3];
        
        
        
        xv_km1_k[7] = vW[0];
        xv_km1_k[8] = vW[1];
        xv_km1_k[9] = vW[2];
        xv_km1_k[10] = wW[0];
        xv_km1_k[11] = wW[1];
        xv_km1_k[12] = wW[2];
        
        
        /**x_k_km1[3][0] = *vW[1][0];
         *x_k_km1[3][1] = *vW[1][1];
         *x_k_km1[3][2] = *vW[1][2];
         *x_k_km1[4][0] = *wW[1][0];
         *x_k_km1[4][1] = *wW[1][1];
         *x_k_km1[4][2] = *wW[1][2];*/
        
        return xv_km1_k;
    }
    return xv_km1_k;
    
}

/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/

double* func_Q(double * xv,double *func_Q_zero_matrix, double*  pn,double delta_t, NSString* type, double* Q) {
    double* res;
    
    /* Actual Code for func_Q */
    
    // Declaring and Initializing the arrays for omegaOld and qOld based on the values from xv
    
    double* omegaOld = (double *)malloc(3*sizeof(double *));
    double* qOld = (double *)malloc(4*sizeof(double *));
    
    omegaOld[0] = xv[10];
    omegaOld[1] = xv[11];
    omegaOld[2] = xv[12];
    
    qOld[0] = xv[3];
    qOld[1] = xv[4];
    qOld[2] = xv[5];
    qOld[3] = xv[6];
    
    double* G = (double *)malloc(13*6*sizeof(double));
    int u,v;
    for (u = 0; u < 13; u++)
    {
        for (v = 0; v < 6; v++)
        {
            G[(6*u)+v] = 0;
        }
    }
    
    
    
    if([type isEqualToString: @"constant_position_and_orientation_location_noise"]) {
        
        
        //Equivalent of declaration of sparse(zeros(13,6)) in MATLAB code
        /*double *G[13][6] = {};
         *G[0][0] = 1*delta_t;
         *G[0][1] = 0*delta_t;
         *G[0][2] = 0*delta_t;
         *G[1][0] = 0*delta_t;
         *G[1][1] = 1*delta_t;
         *G[1][2] = 0*delta_t;
         *G[2][0] = 0*delta_t;
         *G[2][1] = 0*delta_t;
         *G[2][2] = 1*delta_t;*/
        
        
        
        // Making the first 3 rows and columns similar to an identity matrix in G
        
        G[0] = G[7] = G[14] = 1*delta_t;
        
        // Now initializing the matrix q2trRES which is a 4 x 4 matrix obtained from the result of q2tr function. qold (1 x 4 matrix) is the input to the function.
        
        double* q2trRES = (double *)malloc(4*4*sizeof(double *));
        
        int a,b;
        for (a = 0; a < 4; a++)
        {
            for (b = 0; b < 4; b++)
            {
                q2trRES[(4*a)+b] = 0;
            }
        }
        
        double* q2tr(double* qOld, double* q2trRES);
        
        q2trRES = q2tr(qOld, q2trRES); // Returns a 4 x 4 matrix as result
        
        // Declaring and initializing the matrix tr2rpyRES which is a 1 x 3 matrix of all zeros. This is passed along with q2trRES matrix (4 x 4) to update the tr2rpy matrix
        
        double* tr2rpyRES = (double *)malloc(3*sizeof(double));
        
        tr2rpyRES[0] = tr2rpyRES[1] = tr2rpyRES[2] = 0;
        
        // Function declaration for tr2pry which returns tr2rpyRES (1 x 3 matrix)
        
        double* tr2pry(double* q2trRES, double* tr2rpyRES);
        
        tr2rpyRES = tr2rpy(q2trRES, tr2rpyRES); // Returns a 1 x 3 matrix as result
        
        // Initialization of the dq_by_deulerRES matrix (4 x 3 matrix( which is the result of dq_by_deuler function which takes tr2rpyRES (1 x 3 matrix) as input
        
        double* dq_by_deulerRES = (double *)malloc(4*3*sizeof(double));
        
        int y,k;
        
        for (y = 0; y < 4; y++)
        {
            for (k = 0; k < 3; k++)
            {
                dq_by_deulerRES[(3*y)+k] = 0;
            }
        }
        
        double* dq_by_deuler(double* tr2rpyRES, double* dq_by_deulerRES);
        
        dq_by_deulerRES = dq_by_deuler(tr2rpyRES, dq_by_deulerRES);  // Returns a 4 x 3 matrix as output
        
        /* Assigning the values of dq_by_deuler matrix to G from rows 4 through 7 and columns 4 through 6 */
        
        // For Row 4
        G[21] = dq_by_deulerRES[0];
        G[22] = dq_by_deulerRES[1];
        G[23] = dq_by_deulerRES[2];
        
        // For Row 5
        G[27] = dq_by_deulerRES[3];
        G[28] = dq_by_deulerRES[4];
        G[29] = dq_by_deulerRES[5];
        
        //For Row 6
        G[33] = dq_by_deulerRES[6];
        G[34] = dq_by_deulerRES[7];
        G[35] = dq_by_deulerRES[8];
        
        // For Row 7
        G[39] = dq_by_deulerRES[9];
        G[40] = dq_by_deulerRES[10];
        G[41] = dq_by_deulerRES[11];
        
        
        
        
        //Calling the q2tr method passing qOld matrix as the input and obtaining t which is a 4 x 4 matrix as the output
        
        /* double *t[4][4] = {};
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
         *G[6][5] = *d[3][2];*/
        
        
    }
    
    else {
        
        
        //double *qwt[1][4] = {};
        //double* v2q(double *omegaOld[1][3], double *qwt[1][4]);
        
        //Declaration of v2qRES array which has dimensions 1 x 4. The array is updated by the v2q function.
        
        double* v2qRES = (double *)malloc(4*sizeof(double *));
        v2qRES[0] = v2qRES[1] = v2qRES[2] = v2qRES[3] = 0;
        
        double omegaOld_mod1 = omegaOld[0]*delta_t;
        double omegaOld_mod2 = omegaOld[1]*delta_t;
        double omegaOld_mod3 = omegaOld[2]*delta_t;
        
        v2qRES = v2q(omegaOld_mod1, omegaOld_mod2, omegaOld_mod3, v2qRES);
        
        // Making the diagonal elements of rows 8 through 10 and columns 1 through 3 as an identity matrix in G matrix
        
        G[42] = G[49] = G[56] = 1;
        
        // Making the diagonal elements of rows 11 through 13 and columns 4 through 6 as an identity matrix in G matrix
        
        G[63] = G[70] = G[77] = 1;
        
        // Making the diagonal elements of rows 1 through 3 and columns 1 through 3 as an identity matrix in G matrix. Each value is multiplied by delta_t value
        
        G[0] = G[7] = G[14] = 1*delta_t;
        
        
        //double* dq3_by_dq1RES[4][4] = {};
        //double *dqomegadt_by_domegaRES[4][3] = {};
        
        // Declaring and initializing the dq3_by_dq1RES and dqomegadt_by_domegaRES matrices which are 4 x 4 and 4 x 3 matrices respectively
        
        double* dq3_by_dq1RES = (double *)malloc(4*4*sizeof(double *));
        int y,u;
        for (y = 0; y < 4; y++)
        {
            for (u = 0; u < 4; u++)
            {
                dq3_by_dq1RES[(4*y)+u] = 0;
            }
        }
        
        double* dqomegadt_by_domegaRES = (double *)malloc(4*3*sizeof(double *));
        int w, e;
        for (w = 0; w < 4; w++)
        {
            for (e = 0; e < 3; e++)
            {
                dqomegadt_by_domegaRES[(3*w)+e] = 0;
            }
        }
        
        // Calling the functions to returns the updated values for matrices
        
        dq3_by_dq1RES = dq3_by_dq1(qOld, dq3_by_dq1RES);
        
        dqomegadt_by_domegaRES = dqomegadt_by_domega(dqomegadt_by_domegaRES, omegaOld_mod1, omegaOld_mod2, omegaOld_mod3, delta_t);
        
        //The Matrix multiplication is given below. The multiplied result is stored in rows 4 through 7 and columns 10 through 12. It is the same as in the MATLAB file.
        
        
        //NEED TO EDIT THE MATRIX INDICES TO PRODUCE THE CORRECT VALUES UPON MULTIPLICATION
        
        G[21] = dq3_by_dq1RES[0]*dqomegadt_by_domegaRES[0] + dq3_by_dq1RES[1]*dqomegadt_by_domegaRES[3] + dq3_by_dq1RES[2]*dqomegadt_by_domegaRES[6] + dq3_by_dq1RES[3]*dqomegadt_by_domegaRES[9];
        
        G[22] = dq3_by_dq1RES[0]*dqomegadt_by_domegaRES[1] + dq3_by_dq1RES[1]*dqomegadt_by_domegaRES[4] + dq3_by_dq1RES[2]*dqomegadt_by_domegaRES[7] + dq3_by_dq1RES[3]*dqomegadt_by_domegaRES[10];
        
        G[23] = dq3_by_dq1RES[0]*dqomegadt_by_domegaRES[2] + dq3_by_dq1RES[1]*dqomegadt_by_domegaRES[5] + dq3_by_dq1RES[2]*dqomegadt_by_domegaRES[8] + dq3_by_dq1RES[3]*dqomegadt_by_domegaRES[11];
        
        G[27] = dq3_by_dq1RES[4]*dqomegadt_by_domegaRES[0] + dq3_by_dq1RES[5]*dqomegadt_by_domegaRES[3] + dq3_by_dq1RES[6]*dqomegadt_by_domegaRES[6] + dq3_by_dq1RES[7]*dqomegadt_by_domegaRES[9];
        
        G[28] = dq3_by_dq1RES[4]*dqomegadt_by_domegaRES[1] + dq3_by_dq1RES[5]*dqomegadt_by_domegaRES[4] + dq3_by_dq1RES[6]*dqomegadt_by_domegaRES[7] + dq3_by_dq1RES[7]*dqomegadt_by_domegaRES[10];
        
        G[29] = dq3_by_dq1RES[4]*dqomegadt_by_domegaRES[2] + dq3_by_dq1RES[5]*dqomegadt_by_domegaRES[5] + dq3_by_dq1RES[6]*dqomegadt_by_domegaRES[8] + dq3_by_dq1RES[7]*dqomegadt_by_domegaRES[11];
        
        G[33] = dq3_by_dq1RES[8]*dqomegadt_by_domegaRES[0] + dq3_by_dq1RES[9]*dqomegadt_by_domegaRES[3] + dq3_by_dq1RES[10]*dqomegadt_by_domegaRES[6] + dq3_by_dq1RES[11]*dqomegadt_by_domegaRES[9];
        
        G[34] = dq3_by_dq1RES[8]*dqomegadt_by_domegaRES[1] + dq3_by_dq1RES[9]*dqomegadt_by_domegaRES[4] + dq3_by_dq1RES[10]*dqomegadt_by_domegaRES[7] + dq3_by_dq1RES[11]*dqomegadt_by_domegaRES[10];
        
        G[35] = dq3_by_dq1RES[8]*dqomegadt_by_domegaRES[2] + dq3_by_dq1RES[9]*dqomegadt_by_domegaRES[5] + dq3_by_dq1RES[10]*dqomegadt_by_domegaRES[8] + dq3_by_dq1RES[11]*dqomegadt_by_domegaRES[11];
        
        G[39] = dq3_by_dq1RES[12]*dqomegadt_by_domegaRES[0] + dq3_by_dq1RES[13]*dqomegadt_by_domegaRES[3] + dq3_by_dq1RES[14]*dqomegadt_by_domegaRES[6] + dq3_by_dq1RES[15]*dqomegadt_by_domegaRES[9];
        
        G[40] = dq3_by_dq1RES[12]*dqomegadt_by_domegaRES[1] + dq3_by_dq1RES[13]*dqomegadt_by_domegaRES[4] + dq3_by_dq1RES[14]*dqomegadt_by_domegaRES[7] + dq3_by_dq1RES[15]*dqomegadt_by_domegaRES[10];
        
        G[41] = dq3_by_dq1RES[12]*dqomegadt_by_domegaRES[2] + dq3_by_dq1RES[13]*dqomegadt_by_domegaRES[5] + dq3_by_dq1RES[14]*dqomegadt_by_domegaRES[8] + dq3_by_dq1RES[15]*dqomegadt_by_domegaRES[11];
        
    }
    /*//Declaring and Initializing G Transpose Matrix as the transpose of G Matrix
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
     
     double* Q_int1[13][6] = {};
     
     //Matrix Multiplication code goes here for Q_int1[13][6]
     *Q_int1[0][0] = *G[0][0]*pn[0][0] + *G[0][1]*pn[1][0] + *G[0][2]*pn[2][0] + *G[0][3]*pn[3][0] + *G[0][4]*pn[4][0] + *G[0][5]*pn[5][0];
     *Q_int1[0][1] = *G[0][0]*pn[1][0] + *G[0][1]*pn[1][1] + *G[0][2]*pn[2][1] + *G[0][3]*pn[3][1] + *G[0][4]*pn[4][1] + *G[0][5]*pn[5][1];
     
     // Rest of the Multiplication code goes here
     
     
     
     
     }*/
    
    //Converting the G (13 x 6) and Pn (6 x 6) matrices to OpenCV matrices in order to simply matrix multiplication process.
    
    // Also created a transpose of G matrix using OpenCV.
    //Q_mat matrix stores the multiplication result.
    
    cv::Mat G_mat(13, 6, CV_64FC1, G, cv::Mat::AUTO_STEP);
    cv::Mat Pn_mat(6, 6, CV_64FC1, pn, cv::Mat::AUTO_STEP);
    
    cv::Mat G_transpose;
    cv::transpose(G_mat, G_transpose);
    
    cv::Mat Q_mat = G_mat*Pn_mat*G_transpose;
    
    //Code for storing the values in Q_mat matrix to Q matrix and return the Q matrix.
    //Generic data structure:
    int sizeOfQ = Q_mat.rows*Q_mat.cols;
    
    double* Q_ans = (double*)malloc(sizeOfQ*sizeof(double));
    
    //Create a pointer to data in Q_mat:
    double* mxPtr = Q_mat.ptr<double>(0);
    //Iterate through and assign mat values to our data:
    for (int i=0; i<(sizeOfQ); ++i)
        Q_ans[i] = mxPtr[i];
    //Free the data from Q (also make sure to free() any data that you malloc'ed in here and won't use again)
    Q_mat.release();
    
    
    return Q_ans;
}

/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/

double* qprod(double* qWR, double* v2qRES, double* qprodRES) {
    
    // Variable a stores the result of qWR[0]
    double a = qWR[0];
    
    double* v = (double*)malloc(3*1*sizeof(double));
    double* u = (double*)malloc(3*1*sizeof(double));
    
    v[0] = qWR[1];
    v[1] = qWR[2];
    v[2] = qWR[3];
    
    double x = v2qRES[0];
    
    u[0] = v2qRES[1];
    u[1] = v2qRES[2];
    u[2] = v2qRES[3];
    
    double prod_a_x = a*x;
    
    double* v_transpose = (double*)malloc(1*3*sizeof(double));
    v_transpose[0] = v[0];
    v_transpose[1] = v[1];
    v_transpose[2] = v[2];
    
    double vtrans_prod_u = v_transpose[0]*u[0] + v_transpose[1]*u[1] + v_transpose[2]*u[2];
    
    qprodRES[0] = (prod_a_x - vtrans_prod_u);
    
    double* a_prod_umat = (double *)malloc(3*1*sizeof(double));
    
    a_prod_umat[0] = a*u[0];
    a_prod_umat[1] = a*u[1];
    a_prod_umat[2] = a*u[2];
    
    double* x_prod_vmat = (double*)malloc(3*1*sizeof(double));
    
    x_prod_vmat[0] = x*v[0];
    x_prod_vmat[1] = x*v[1];
    x_prod_vmat[2] = x*v[2];
    
    double* au_plus_xv = (double*)malloc(3*1*sizeof(double));
    
    au_plus_xv[0] = a_prod_umat[0] + x_prod_vmat[0];
    au_plus_xv[1] = a_prod_umat[1] + x_prod_vmat[1];
    au_plus_xv[2] = a_prod_umat[2] + x_prod_vmat[2];
    
    double* au_plus_xv_transpose = (double*)malloc(1*3*sizeof(double));
    
    au_plus_xv_transpose[0] = au_plus_xv[0];
    au_plus_xv_transpose[1] = au_plus_xv[1];
    au_plus_xv_transpose[2] = au_plus_xv[2];
    
    // Code for the Cross Product Calculations
    
    double* cross_v_and_u = (double*)malloc(3*1*sizeof(double));
    
    cross_v_and_u[0] = ((v[1]*u[2]) - (v[2]*u[1]));
    cross_v_and_u[1] = ((v[2]*u[0]) - (v[0]*u[2]));
    cross_v_and_u[2] = ((v[0]*u[1]) - (v[1]*u[0]));
    
    double* cross_v_and_u_transpose = (double*)malloc(1*3*sizeof(double));
    
    cross_v_and_u_transpose[0] = cross_v_and_u[0];
    cross_v_and_u_transpose[1] = cross_v_and_u[1];
    cross_v_and_u_transpose[2] = cross_v_and_u[2];
    
    double* qprodsecond_term = (double*)malloc(1*3*sizeof(double));
    
    qprodsecond_term[0] = au_plus_xv_transpose[0] + cross_v_and_u_transpose[0];
    qprodsecond_term[1] = au_plus_xv_transpose[1] + cross_v_and_u_transpose[1];
    qprodsecond_term[2] = au_plus_xv_transpose[2] + cross_v_and_u_transpose[2];
    
    qprodRES[1] = qprodsecond_term[0];
    qprodRES[2] = qprodsecond_term[1];
    qprodRES[3] = qprodsecond_term[2];
    
    
    return qprodRES;
}


/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/


// Method to calculate the Jacobian Derivatives. Implemented as given in the MATLAB code
//double* dfv_by_dxv(double *x_k_k[1][13], double *dfv_by_dxv_zero_matrix[6][1], double delta_t, char* type) {
double* dfv_by_dxv(double* dfv_by_dxvRES, double* xv, double* u, double dt, NSString* type){
    
    if (dfv_by_dxvRES==NULL || xv==NULL || u==NULL || type==NULL) {
        return FALSE;
    }
    double delta_t = DELTA_T;
    
    /* if(type != @"constant_velocity" || type != @"constant_orientation" || type != @"constant_position" || type != @"constant_position_and_orientation")
     return FALSE;*/
    
    if(!([type isEqualToString:@"constant_velocity"]|| [type isEqualToString:@"constant_orientation"] || [type isEqualToString:@"constant_position"] || [type isEqualToString:@"constant_position_and_orientation"]))
        return FALSE;
    // Code for the State Transition Equation Derivative and other Calculations
    // double *rettype;
    // List of functions to call
    /*void dq3_by_dq2();
     void dq3_by_dq1();
     void dqomegadt_by_domega();
     void dq0_by_domegaA();*/
    
    // End of the list
    
    // Code translated from the dfv_by_dxv.m MATLAB file
    //double *omegaOld[1][3] = {x_k_k[0], x_k_k[1], x_k_k[2]};
    //double *omegaOld = (double*)malloc(1*3*sizeof(double));
    /*omegaOld[0] = xv[0];
     omegaOld[1] = xv[1];
     omegaOld[2] = xv[2];*/
    double omega1 = xv[10];
    double omega2 = xv[11];
    double omega3 = xv[12];
    //double qOld[4] = {xv[3], xv[4], xv[5], xv[6]};
    double *qOld = (double*)malloc(1*4*sizeof(double));
    qOld[0] = xv[3];
    qOld[1] = xv[4];
    qOld[2] = xv[5];
    qOld[3] = xv[6];
    // double* dq3_by_dq1RES[4][4] = {};
    double *dq3_by_dq1RES = (double*)malloc(4*4*sizeof(double));
    
    int c, d;
    for (c = 0; c <=3; c++) {
        for (d = 0; d <=3; d++) {
            dq3_by_dq1RES[(4*c)+d] = 0;
        }
    }
    //double *dqomegadt_by_domegaRES[4][3] = {};
    double *dqomegadt_by_domegaRES = (double*)malloc(4*3*sizeof(double));
    int i, j;
    for(i = 0; i <4; i++) {
        for (j = 0; j <3; j++) {
            dqomegadt_by_domegaRES[(3*i)+j] = 0;
        }
    }
    
    
    //This step creates a 13x13 Identity Matrix known as fv_by_dxvRES
    // double *dfv_by_dxvRES[13][13] = {};
    /* double *dfv_by_dxvRES = (double*)malloc(13*13*sizeof(double));
     int m,n;
     //Made changes to the loop to ensure that there is no error. This loop has to be tested
     // BAD ACCESS crash - probably because the memory was never allocated
     for(m = 0; m <=12; m++) {
     for (n = 0; n <=12; n++) {
     dfv_by_dxvRES[(12*m)+n] = 0;
     }
     
     }*/
    
    
    //END OF THE LOOP
    //double *qwt[1][4] = {};
    double *qwt = (double*)malloc(1*4*sizeof(double));
    qwt[0] = qwt[1] = qwt[2] = qwt[3] = 0.0;
    
    //dfv_by_dxvRES[1] = 0;
    
    // Setting identity matrix values to the dfv_by_dxv matrix as given in the MATLAB code
    int l,n;
    for (l = 0; l < 13; l++)
    {
        for (n = 0; n < 13; n++)
        {
            if(n ==l) {
                dfv_by_dxvRES[(13*l)+l] = 1;
            }
            else {
                dfv_by_dxvRES[(13*l)+n] = 0;
            }
        }
        
    }
    
    double omega1_mod = omega1*delta_t;
    double omega2_mod = omega2*delta_t;
    double omega3_mod = omega3*delta_t;
    
    // double* v2q(double omega1_mod, double omega2_mod, double omega3_mod, double *qwt);
    
    qwt = v2q(omega1_mod, omega2_mod, omega3_mod, qwt);
    double *dq3_by_dq2(double* qwt, double*dfv_by_dxvRES);
    
    dfv_by_dxvRES = dq3_by_dq2(qwt, dfv_by_dxvRES);
    
    //Assigning the identity matrix values to the dfv_by_dxvRES matrix as given in the MATLAB code
    if([type  isEqualToString: @"constant_velocity"]) {
        dfv_by_dxvRES[7] = 1*delta_t;
        dfv_by_dxvRES[21] = 1*delta_t;
        dfv_by_dxvRES[35] = 1*delta_t;
        //dfv_by_dxvRES(4:7,11:13) = dq3_by_dq1(qOld)*dqomegadt_by_domega(omegaOld,dt); implementation goes here
        //Need to write the functions dq3_by_dq1(qOld) and qomegadt_by_domega(omegaOld,dt)
        
        
        dq3_by_dq1RES = dq3_by_dq1(qOld, dq3_by_dq1RES);
        
        
        dqomegadt_by_domegaRES = dqomegadt_by_domega(dqomegadt_by_domegaRES, omega1, omega2, omega3, delta_t);
        //The Matrix multiplication is given below. The multiplied result is stored in rows 4 through 7 and columns 10 through 12. It is the same as in the MATLAB file
        
        dfv_by_dxvRES[49] = dq3_by_dq1RES[0]*dqomegadt_by_domegaRES[0] + dq3_by_dq1RES[1]*dqomegadt_by_domegaRES[3] + dq3_by_dq1RES[2]*dqomegadt_by_domegaRES[6] + dq3_by_dq1RES[3]*dqomegadt_by_domegaRES[9];
        dfv_by_dxvRES[50] = dq3_by_dq1RES[0]*dqomegadt_by_domegaRES[1] + dq3_by_dq1RES[1]*dqomegadt_by_domegaRES[4] + dq3_by_dq1RES[2]*dqomegadt_by_domegaRES[7] + dq3_by_dq1RES[3]*dqomegadt_by_domegaRES[10];
        dfv_by_dxvRES[51] = dq3_by_dq1RES[0]*dqomegadt_by_domegaRES[2] + dq3_by_dq1RES[1]*dqomegadt_by_domegaRES[5] + dq3_by_dq1RES[2]*dqomegadt_by_domegaRES[8] + dq3_by_dq1RES[3]*dqomegadt_by_domegaRES[11];
        dfv_by_dxvRES[62] = dq3_by_dq1RES[4]*dqomegadt_by_domegaRES[0] + dq3_by_dq1RES[5]*dqomegadt_by_domegaRES[3] + dq3_by_dq1RES[6]*dqomegadt_by_domegaRES[6] + dq3_by_dq1RES[7]*dqomegadt_by_domegaRES[9];
        dfv_by_dxvRES[63] = dq3_by_dq1RES[4]*dqomegadt_by_domegaRES[1] + dq3_by_dq1RES[5]*dqomegadt_by_domegaRES[4] + dq3_by_dq1RES[6]*dqomegadt_by_domegaRES[7] + dq3_by_dq1RES[7]*dqomegadt_by_domegaRES[10];
        dfv_by_dxvRES[64] = dq3_by_dq1RES[4]*dqomegadt_by_domegaRES[2] + dq3_by_dq1RES[5]*dqomegadt_by_domegaRES[5] + dq3_by_dq1RES[6]*dqomegadt_by_domegaRES[8] + dq3_by_dq1RES[7]*dqomegadt_by_domegaRES[11];
        dfv_by_dxvRES[75] = dq3_by_dq1RES[8]*dqomegadt_by_domegaRES[0] + dq3_by_dq1RES[9]*dqomegadt_by_domegaRES[3] + dq3_by_dq1RES[10]*dqomegadt_by_domegaRES[6] + dq3_by_dq1RES[11]*dqomegadt_by_domegaRES[9];
        dfv_by_dxvRES[76] = dq3_by_dq1RES[8]*dqomegadt_by_domegaRES[1] + dq3_by_dq1RES[9]*dqomegadt_by_domegaRES[4] + dq3_by_dq1RES[10]*dqomegadt_by_domegaRES[7] + dq3_by_dq1RES[11]*dqomegadt_by_domegaRES[10];
        dfv_by_dxvRES[77] = dq3_by_dq1RES[8]*dqomegadt_by_domegaRES[2] + dq3_by_dq1RES[9]*dqomegadt_by_domegaRES[5] + dq3_by_dq1RES[10]*dqomegadt_by_domegaRES[8] + dq3_by_dq1RES[11]*dqomegadt_by_domegaRES[11];
        dfv_by_dxvRES[88] = dq3_by_dq1RES[12]*dqomegadt_by_domegaRES[0] + dq3_by_dq1RES[13]*dqomegadt_by_domegaRES[3] + dq3_by_dq1RES[14]*dqomegadt_by_domegaRES[6] + dq3_by_dq1RES[15]*dqomegadt_by_domegaRES[9];
        dfv_by_dxvRES[89] = dq3_by_dq1RES[12]*dqomegadt_by_domegaRES[1] + dq3_by_dq1RES[13]*dqomegadt_by_domegaRES[4] + dq3_by_dq1RES[14]*dqomegadt_by_domegaRES[7] + dq3_by_dq1RES[15]*dqomegadt_by_domegaRES[10];
        dfv_by_dxvRES[90] = dq3_by_dq1RES[12]*dqomegadt_by_domegaRES[2] + dq3_by_dq1RES[13]*dqomegadt_by_domegaRES[5] + dq3_by_dq1RES[14]*dqomegadt_by_domegaRES[8] + dq3_by_dq1RES[15]*dqomegadt_by_domegaRES[11];
        
        
    }
    
    else if([type isEqualToString: @"constant_orientation"]) {
        dfv_by_dxvRES[49] = 0;
        dfv_by_dxvRES[50] = 0;
        dfv_by_dxvRES[51] = 0;
        dfv_by_dxvRES[62] = 0;
        dfv_by_dxvRES[63] = 0;
        dfv_by_dxvRES[64] = 0;
        dfv_by_dxvRES[75] = 0;
        dfv_by_dxvRES[76] = 0;
        dfv_by_dxvRES[77] = 0;
        dfv_by_dxvRES[88] = 0;
        dfv_by_dxvRES[89] = 0;
        dfv_by_dxvRES[90] = 0;
        
        dfv_by_dxvRES[140] = 0;
        dfv_by_dxvRES[141] = 0;
        dfv_by_dxvRES[142] = 0;
        dfv_by_dxvRES[153] = 0;
        dfv_by_dxvRES[154] = 0;
        dfv_by_dxvRES[155] = 0;
        dfv_by_dxvRES[166] = 0;
        dfv_by_dxvRES[167] = 0;
        dfv_by_dxvRES[168] = 0;
    }
    
    else if([type isEqualToString: @"constant_position"]) {
        dfv_by_dxvRES[7] = 0;
        dfv_by_dxvRES[8] = 0;
        dfv_by_dxvRES[9] = 0;
        dfv_by_dxvRES[20] = 0;
        dfv_by_dxvRES[21] = 0;
        dfv_by_dxvRES[22] = 0;
        dfv_by_dxvRES[33] = 0;
        dfv_by_dxvRES[34] = 0;
        dfv_by_dxvRES[35] = 0;
        dfv_by_dxvRES[98] = 0;
        dfv_by_dxvRES[99] = 0;
        dfv_by_dxvRES[100] = 0;
        dfv_by_dxvRES[111] = 0;
        dfv_by_dxvRES[112] = 0;
        dfv_by_dxvRES[113] = 0;
        dfv_by_dxvRES[124] = 0;
        dfv_by_dxvRES[125] = 0;
        dfv_by_dxvRES[126] = 0;
    }
    
    else if([type isEqualToString: @"constant_position_and_orientation"]) {
        dfv_by_dxvRES[49] = 0;
        dfv_by_dxvRES[50] = 0;
        dfv_by_dxvRES[51] = 0;
        dfv_by_dxvRES[62] = 0;
        dfv_by_dxvRES[63] = 0;
        dfv_by_dxvRES[64] = 0;
        dfv_by_dxvRES[75] = 0;
        dfv_by_dxvRES[76] = 0;
        dfv_by_dxvRES[77] = 0;
        dfv_by_dxvRES[88] = 0;
        dfv_by_dxvRES[89] = 0;
        dfv_by_dxvRES[90] = 0;
        dfv_by_dxvRES[7] = 0;
        dfv_by_dxvRES[8] = 0;
        dfv_by_dxvRES[9] = 0;
        dfv_by_dxvRES[20] = 0;
        dfv_by_dxvRES[21] = 0;
        dfv_by_dxvRES[22] = 0;
        dfv_by_dxvRES[33] = 0;
        dfv_by_dxvRES[34] = 0;
        dfv_by_dxvRES[35] = 0;
        dfv_by_dxvRES[140] = 0;
        dfv_by_dxvRES[141] = 0;
        dfv_by_dxvRES[142] = 0;
        dfv_by_dxvRES[153] = 0;
        dfv_by_dxvRES[154] = 0;
        dfv_by_dxvRES[155] = 0;
        dfv_by_dxvRES[166] = 0;
        dfv_by_dxvRES[167] = 0;
        dfv_by_dxvRES[168] = 0;
        dfv_by_dxvRES[98] = 0;
        dfv_by_dxvRES[99] = 0;
        dfv_by_dxvRES[100] = 0;
        dfv_by_dxvRES[111] = 0;
        dfv_by_dxvRES[112] = 0;
        dfv_by_dxvRES[113] = 0;
        dfv_by_dxvRES[124] = 0;
        dfv_by_dxvRES[125] = 0;
        dfv_by_dxvRES[126] = 0;
        return dfv_by_dxvRES;
        
        
        // Calculate commonly used Jacobian part dq(omega * delta_t) by domega code goes here
        
    }
    
    else
    {
        //return FALSE;
    }
    
    
    // return rettype;*/
    return dfv_by_dxvRES;
    //return true;
    //return FALSE;
}

/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/

//double * dqomegadt_by_domega(double *omega[1][3],double delta_t, double *dqomegadt_by_domegaRES[4][3]) {
double* dqomegadt_by_domega(double *dqomegadt_by_domegaRES, double omega1, double omega2, double omega3, double delta_t){
    
    if(dqomegadt_by_domegaRES == NULL)
        return FALSE;
    
    double value1 = omega1;
    double squared_value1 = value1*value1;
    double value2 = omega2;
    double squared_value2 = value2*value2;
    double value3 = omega3;
    double squared_value3 = value3*value3;
    
    
    // The modulus of omega is obtained as the square root of the sum of the aquared values
    double omegamod = sqrt(squared_value1 + squared_value2 + squared_value3);
    
    // Use generic ancillary functions to calculate components of Jacobian
    //dq0_by_domegaARES=(-delta_t / 2.0) * (omegaA / omega) * sin(omega * delta_t / 2.0);
    dqomegadt_by_domegaRES[0] = (((-1)*delta_t)/2)*(omega1/omegamod)*sin((omegamod*delta_t)/2);
    dqomegadt_by_domegaRES[1] = (((-1)*delta_t)/2)*(omega2/omegamod)*sin((omegamod*delta_t)/2);
    dqomegadt_by_domegaRES[2] = (((-1)*delta_t)/2)*(omega3/omegamod)*sin((omegamod*delta_t)/2);
    
    dqomegadt_by_domegaRES[3] = (delta_t/2)*omega1*omega1/(omegamod*omegamod)*cos((omegamod*delta_t)/2)+(1/omegamod)*(1-omega1*omega1/(omegamod*omegamod))*sin((omegamod*delta_t)/2);
    dqomegadt_by_domegaRES[4] = ((omega1*omega2)/(omegamod*omegamod))*(((delta_t)/2)*cos((omegamod*delta_t)/2) -((1/omegamod)*sin((omegamod*delta_t)/2)));
    dqomegadt_by_domegaRES[5] = ((omega1*omega3)/(omegamod*omegamod))*(((delta_t)/2)*cos((omegamod*delta_t)/2) -((1/omegamod)*sin((omegamod*delta_t)/2)));
    dqomegadt_by_domegaRES[6] = ((omega2*omega1)/(omegamod*omegamod))*(((delta_t)/2)*cos((omegamod*delta_t)/2) -((1/omegamod)*sin((omegamod*delta_t)/2)));
    dqomegadt_by_domegaRES[7] = (delta_t/2)*omega2*omega2/(omegamod*omegamod)*cos((omegamod*delta_t)/2)+(1/omegamod)*(1-omega2*omega2/(omegamod*omegamod))*sin((omegamod*delta_t)/2);
    
    dqomegadt_by_domegaRES[8] = ((omega2*omega3)/(omegamod*omegamod))*(((delta_t)/2)*cos((omegamod*delta_t)/2) -((1/omegamod)*sin((omegamod*delta_t)/2)));
    dqomegadt_by_domegaRES[9] = ((omega3*omega1)/(omegamod*omegamod))*(((delta_t)/2)*cos((omegamod*delta_t)/2) -((1/omegamod)*sin((omegamod*delta_t)/2)));
    dqomegadt_by_domegaRES[10] = ((omega3*omega2)/(omegamod*omegamod))*(((delta_t)/2)*cos((omegamod*delta_t)/2) -((1/omegamod)*sin((omegamod*delta_t)/2)));
    dqomegadt_by_domegaRES[11] = (delta_t/2)*omega3*omega3/(omegamod*omegamod)*cos((omegamod*delta_t)/2)+(1/omegamod)*(1-omega3*omega3/(omegamod*omegamod))*sin((omegamod*delta_t)/2);
    
    //return dqomegadt_by_domegaRES;
    return dqomegadt_by_domegaRES;
}