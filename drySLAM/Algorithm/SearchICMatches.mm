//
//  MapManagement.m
//  drySLAM
//
//  Created by Anthony Rodriguez on 4/17/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//

#include "SearchICMatches.hpp"
#include "math.h"


bool search_IC_matches(Filter* filter, NSMutableArray* features_info, Camera* cam, cv::Mat* im)
{
    int i;
    
    predict_camera_measurements(features_info, filter, cam);
    calculateDerivatives(features_info, filter, cam);
    
    for (i=0; i<(int)[features_info count]; i++)
    {
        Feature *current = (Feature *)features_info[i];
        
        if (current->h != NULL)
        {
            int size = filter->state_km1_size;
            
            cv::Mat *H = new cv::Mat(size,size,CV_64FC1,current->H);
            cv::Mat *H_transpose =new cv::Mat(size,size,CV_64FC1,0.0);
            cv::transpose(*H, *H_transpose);
            
            cv::Mat *p_k_km1 =  new cv::Mat(size,size,CV_64FC1,filter->p_k_km1);
            cv::Mat *R = new cv::Mat(size,size,CV_64FC1,current->R);
            cv::Mat *S = new cv::Mat(size,size,CV_64FC1,0.0);
            *S = (*H) * (*p_k_km1) * (*H_transpose) + (*R);
        }
        
        predict_features_appearance(features_info, filter, cam);
        
        matching(im, features_info, cam);
    }
    
    return true;
}


/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/

void calculateDerivatives(NSMutableArray* features_info, Filter* filter,Camera* cam)
{
    double *x_v;
    double *x_features0;
    
    double *y;
    double *x_features;
    
    int i,j;
    
    for (i=0;i<13;i++)
    {
        x_v=&(filter->x_k_km1[i]);
        x_v++;
    }
    
    for (i=13;filter->x_k_km1+i;i++)
    {
        x_features0=&(filter->x_k_km1[i]);
        x_features0++;
    }
    
    for (i=0;i<[features_info count];i++)
    {
        Feature *current=features_info[i];
        
        if (sizeof(current->h)) //if ~isempty(features_info(i).h)
        {
            if ([current->type isEqualToString:@"cartesian"])
            {
                for (j=0;j<3;j++)
                {
                    y=&(x_features0[j]); //syntax
                    y++;
                }
                
                for (j=3;x_features0+j;j++)
                {
                    x_features=&(x_features0[j]);
                    x_features++;
                }
                
                //current->H=sparse();
            }
            
            else
            {
                for (j=0;j<6;j++)
                {
                    y=&(x_features0[j]);
                    y++;
                }
                
                for (j=6;x_features0+j;j++)
                {
                    x_features=&(x_features0[j]);
                    x_features++;
                }
                
                //current->H=sparse();
            }
        }
        
        else
        {
            if ([current->type isEqualToString:@"cartesian"])
            {
                for (j=3;x_features0+j;j++)
                {
                    x_features=&(x_features0[j]);
                    x_features++;
                }
            }
            
            if ([current->type isEqualToString:@"inversedepth"])
            {
                for (j=6;x_features0+j;j++)
                {
                    x_features=&(x_features0[j]);
                    x_features++;
                }
            }
        }
    }
}


/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/

bool calculate_hi_inverse_depth(double* Hi, double* inverse_depth_features_index, double* cartesian_features_index,
                                double* Xv_km1_k, double* yi, Camera* cam,
                                int i, NSMutableArray* features_info){
    if (features_info==nil || [features_info count]<i || Hi==NULL || cartesian_features_index==NULL || Xv_km1_k==NULL || yi==NULL || cam==nil)
        return false;
    double* zi = (double*)malloc(2*sizeof(double));
    zi[0] = ((Feature*)(features_info[i]))->h[0];
    zi[1] = ((Feature*)(features_info[i]))->h[1];
    if (zi==NULL)
        return false;
    int numberOfFeatures = (int)[features_info count];
    int* id_featuresIdx = (int*)malloc(numberOfFeatures*sizeof(int));
    int* cartesian_featuresIdx = (int*)malloc(numberOfFeatures*sizeof(int));
    int id_sum = 0;
    int id_sumBeforeIndex_i = 0;
    int cart_sum = 0;
    int cart_sumBeforeIndex_i = 0;
    
    for (int j=0; j<numberOfFeatures; ++j){
        Feature* feat = (Feature*)features_info[j];
        if (feat->type==nil)
            return false;
        else if ([feat->type isEqualToString:@"inversedepth"]){
            if (j<i)
                ++id_sumBeforeIndex_i;
            id_sum = id_sum+1;
            id_featuresIdx[j] = 1;
            cartesian_featuresIdx[j] = 0;
        }
        //TODO: Fix this!
        else if ([feat->type isEqualToString:@"cartesian"]){
            if (j<i)
                ++cart_sumBeforeIndex_i;
            cart_sum = cart_sum+1;
            id_featuresIdx[j] = 0;
            cartesian_featuresIdx[j] = 1;
            NSLog(@"Error: Cartesian operations not supported (calculate_hi_inverse_depth)");
            return false;
        }
        else{
            NSLog(@"Error: calculate_hi_inverse_depth encountered undefined type");
            return false;
        }
    }
    
    //This should be how Hi is created in calling method:
    //double* Hi = (double*)malloc(2*(13+3*id_sum+6*cart_sum)*sizeof(double));
    int hiNumCols = (13+3*cart_sum+6*id_sum);
    int hiSize = 2*hiNumCols;
    for (int i=0; i<hiSize; ++i){
        Hi[i] = 0.0;
    }
    
    double* Hi_1to13 = (double*)malloc(26*sizeof(double));
    for (int i=0; i<26; ++i){
        Hi_1to13[i] = 0.0;
    }
    double* Hi_y = (double*)malloc(12*sizeof(double));
    for (int i=0; i<12; ++i){
        Hi_y[i] = 0.0;
    }
    
    bool success1 = dh_dxv(Hi_1to13, cam, Xv_km1_k, yi, zi);
    bool success2 = dh_dy(Hi_y, cam, Xv_km1_k, yi, zi);
    
    int insertionIdx = 12+3*cart_sumBeforeIndex_i + 6*id_sumBeforeIndex_i +1;
    
    //Fill first 13 cols of both rows:
    Hi[0] = Hi_1to13[0];
    Hi[1] = Hi_1to13[1];
    Hi[2] = Hi_1to13[2];
    Hi[3] = Hi_1to13[3];
    Hi[4] = Hi_1to13[4];
    Hi[5] = Hi_1to13[5];
    Hi[6] = Hi_1to13[6];
    Hi[7] = Hi_1to13[7];
    Hi[8] = Hi_1to13[8];
    Hi[9] = Hi_1to13[9];
    Hi[10] = Hi_1to13[10];
    Hi[11] = Hi_1to13[11];
    Hi[12] = Hi_1to13[12];
    Hi[hiNumCols] = Hi_1to13[13];
    Hi[hiNumCols + 1] = Hi_1to13[14];
    Hi[hiNumCols + 2] = Hi_1to13[15];
    Hi[hiNumCols + 3] = Hi_1to13[16];
    Hi[hiNumCols + 4] = Hi_1to13[17];
    Hi[hiNumCols + 5] = Hi_1to13[18];
    Hi[hiNumCols + 6] = Hi_1to13[19];
    Hi[hiNumCols + 7] = Hi_1to13[20];
    Hi[hiNumCols + 8] = Hi_1to13[21];
    Hi[hiNumCols + 9] = Hi_1to13[22];
    Hi[hiNumCols + 10] = Hi_1to13[23];
    Hi[hiNumCols + 11] = Hi_1to13[24];
    Hi[hiNumCols + 12] = Hi_1to13[25];
    //Fill i col and 5 rows immediately after:
    Hi[insertionIdx] = Hi_y[0];
    Hi[insertionIdx+1] = Hi_y[1];
    Hi[insertionIdx+2] = Hi_y[2];
    Hi[insertionIdx+3] = Hi_y[3];
    Hi[insertionIdx+4] = Hi_y[4];
    Hi[insertionIdx+5] = Hi_y[5];
    int r2i = hiNumCols+insertionIdx;
    Hi[r2i] = Hi_y[6];
    Hi[r2i+1] = Hi_y[7];
    Hi[r2i+2] = Hi_y[8];
    Hi[r2i+3] = Hi_y[9];
    Hi[r2i+4] = Hi_y[10];
    Hi[r2i+5] = Hi_y[11];
    
    free(zi);
    free(id_featuresIdx);
    free(cartesian_featuresIdx);
    free(Hi_1to13);
    free(Hi_y);
    
    
    return success1&&success2;
}


/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/

bool calculate_hi_cartesian(double* Hi, double* inverse_depth_features_index, double* cartesian_features_index,
                            double* Xv_km1_k, double* yi, Camera* cam,
                            int i, NSMutableArray* features_info){
    
    
    
    return false;
}

/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
//Miscellaneous private functions for dh_dxv:

//called by dh_dhrl
bool dhd_dhu(double* a, Camera* cam, double* zi_d){
    //a is a 2x2 matrix
    //zi_d is a 2-element vector
    cv::Mat* inv_a = new cv::Mat(2,2,CV_64FC1,0.0);
    JacobUndistorFM(cam, zi_d, inv_a);
    cv::invert(*inv_a, *inv_a);
    double* mxPtr = inv_a->ptr<double>(0);
    a[0] = mxPtr[0];
    a[1] = mxPtr[1];
    a[2] = mxPtr[2];
    a[3] = mxPtr[3];
    inv_a->release();
    return true;
}

//called by dh_dhrl
bool dhu_dhrl(double* a, Camera* cam, double* Xv_km1_k, double* yi){
    //a is a 2x3 matrix
    double f = cam->f;
    double ku = 1/cam->dy;
    double kv = 1/cam->dy;
    double* rw = (double*)malloc(3*sizeof(double));
    rw[0] = Xv_km1_k[0];
    rw[1] = Xv_km1_k[1];
    rw[2] = Xv_km1_k[2];
    double* Xv_km1_k_4to7 = (double*)malloc(4*sizeof(double));
    Xv_km1_k_4to7[0] = Xv_km1_k[3];
    Xv_km1_k_4to7[1] = Xv_km1_k[4];
    Xv_km1_k_4to7[2] = Xv_km1_k[5];
    Xv_km1_k_4to7[3] = Xv_km1_k[6];
    double* R = (double*)malloc(9*sizeof(double));
    bool success = q2r(Xv_km1_k_4to7, R);
    cv::Mat Rrw(3,3,CV_64FC1, R, cv::Mat::AUTO_STEP);
    cv::invert(Rrw, Rrw);
    double theta = yi[3];
    double phi = yi[4];
    double rho = yi[5];
    double* mi_data = (double*)malloc(3*sizeof(double));
    mi_data[0] = cos(phi)*sin(theta);
    mi_data[1] = -sin(phi);
    mi_data[2] = cos(phi)*cos(theta);
    cv::Mat mi(1,3,CV_64FC1,mi_data,cv::Mat::AUTO_STEP);
    cv::transpose(mi,mi);
    
    double* yi_1to3_minus_rw = (double*)malloc(3*sizeof(double));
    yi_1to3_minus_rw[0] = (yi[0]-rw[0])*rho;
    yi_1to3_minus_rw[1] = (yi[1]-rw[1])*rho;
    yi_1to3_minus_rw[2] = (yi[2]-rw[2])*rho;
    
    double* mxPtr = mi.ptr<double>(0);
    mxPtr[0] = yi_1to3_minus_rw[0]+mxPtr[0];
    mxPtr[1] = yi_1to3_minus_rw[1]+mxPtr[1];
    mxPtr[2] = yi_1to3_minus_rw[2]+mxPtr[2];
    
    cv::Mat hc = Rrw*mi;
    mxPtr = hc.ptr<double>(0);
    double hcx = mxPtr[0];
    double hcy = mxPtr[1];
    double hcz = mxPtr[2];
    
    a[0] = (f*ku)/hcz;
    a[1] = 0;
    a[2] = -1*hcx*f*ku/(hcz*hcz);
    a[3] = 0;
    a[4] = f*kv/hcz;
    a[5] = -1*hcy*f*kv/(hcz*hcz);
    
    free(rw);
    free(Xv_km1_k_4to7);
    free(R);
    Rrw.release();
    free(mi_data);
    mi.release();
    free(yi_1to3_minus_rw);
    hc.release();
    
    return success;
}

//called by dh_drw and dh_dy
bool dh_dhrl(double* a, Camera* cam, double* Xv_km1_k, double* yi, double* zi){
    //a is a 2x3 matrix
    //dhd_dhu is a 2x2 matrix
    //dhu_dhrl is a 2x3 matrix
    double* dhd_dhu_data = (double*)malloc(4*sizeof(double));
    double* dhu_dhrl_data = (double*)malloc(6*sizeof(double));
    
    bool success1 = dhd_dhu(dhd_dhu_data, cam, zi);
    bool success2 = dhu_dhrl(dhu_dhrl_data, cam, Xv_km1_k, yi);
    
    cv::Mat dhd_dhu_mat(2,2,CV_64FC1,dhd_dhu_data,cv::Mat::AUTO_STEP);
    cv::Mat dhu_dhrl_mat(2,3,CV_64FC1,dhu_dhrl_data,cv::Mat::AUTO_STEP);
    cv::Mat ans = dhd_dhu_mat*dhu_dhrl_mat;
    
    double* mxPtr = ans.ptr<double>(0);
    
    a[0] = mxPtr[0];
    a[1] = mxPtr[1];
    a[2] = mxPtr[2];
    a[3] = mxPtr[3];
    a[4] = mxPtr[4];
    a[5] = mxPtr[5];
    
    free(dhd_dhu_data);
    free(dhu_dhrl_data);
    //dhd_dhu_mat.release();
    //dhu_dhrl_mat.release();
    ans.release();
    
    return success1&&success2;
}

//called by dh_drw
bool dhrl_drw(double* a, double* Xv_km1_k, double* yi){
    //a is a 3x3 matrix
    double* Xv_4to7 = (double*)malloc(4*sizeof(double));
    Xv_4to7[0] = Xv_km1_k[3];
    Xv_4to7[1] = Xv_km1_k[4];
    Xv_4to7[2] = Xv_km1_k[5];
    Xv_4to7[3] = Xv_km1_k[6];
    double* R = (double*)malloc(9*sizeof(double));
    bool success = q2r(Xv_4to7, R);
    cv::Mat R_mat(3, 3, CV_64FC1, R, cv::Mat::AUTO_STEP);
    cv::invert(R_mat,R_mat);
    
    double scalarCoeff = -1*yi[5];
    double* mxPtr = R_mat.ptr<double>(0);
    
    a[0] = mxPtr[0]*scalarCoeff;
    a[1] = mxPtr[1]*scalarCoeff;
    a[2] = mxPtr[2]*scalarCoeff;
    a[3] = mxPtr[3]*scalarCoeff;
    a[4] = mxPtr[4]*scalarCoeff;
    a[5] = mxPtr[5]*scalarCoeff;
    a[6] = mxPtr[6]*scalarCoeff;
    a[7] = mxPtr[7]*scalarCoeff;
    a[8] = mxPtr[8]*scalarCoeff;
    
    free(Xv_4to7);
    free(R);
    R_mat.release();
    
    return success;
}

//called by dh_dxv
bool dh_drw(double* Hi1, Camera* cam, double* Xv_km1_k, double* yi, double* zi){
    //Hi1 is a 2x3 matrix
    //dh_dhrl is a 2x3 matrix
    //dhrl_drw is a 3x3 matrix
    //Multiplication: (2x3)x(3x3) = (2x3)
    double* dh_dhrl_data = (double*)malloc(6*sizeof(double));
    double* dhrl_drw_data = (double*)malloc(9*sizeof(double));
    
    bool success1 = dh_dhrl(dh_dhrl_data, cam, Xv_km1_k, yi, zi);
    bool success2 = dhrl_drw(dhrl_drw_data, Xv_km1_k, yi);
    
    cv::Mat dh_dhrl_mat(2, 3, CV_64FC1, dh_dhrl_data, cv::Mat::AUTO_STEP);
    cv::Mat dhrl_drw_mat(3, 3, CV_64FC1, dhrl_drw_data, cv::Mat::AUTO_STEP);
    cv::Mat ans = dh_dhrl_mat*dhrl_drw_mat;
    
    double* mxPtr = ans.ptr<double>(0);
    Hi1[0] = mxPtr[0];
    Hi1[1] = mxPtr[1];
    Hi1[2] = mxPtr[2];
    Hi1[3] = mxPtr[3];
    Hi1[4] = mxPtr[4];
    Hi1[5] = mxPtr[5];
    
    free(dh_dhrl_data);
    free(dhrl_drw_data);
    dh_dhrl_mat.release();
    dhrl_drw_mat.release();
    ans.release();
    
    return success1&&success2;
}

//called by dh_dqwr
//dependency: dRq_times_a_by_dq
bool dhrl_dqwr(double* a, double* Xv_km1_k, double* yi){
    //a is a 3x4 matrix
    double rw[3] = {Xv_km1_k[0], Xv_km1_k[1], Xv_km1_k[2]};
    //qconj is taken care of in this declaration:
    double qwr[4] = {Xv_km1_k[3], -1*Xv_km1_k[4], -1*Xv_km1_k[5], -1*Xv_km1_k[6]};
    double lambda = yi[5];
    double phi = yi[4];
    double theta = yi[3];
    double mi[3] = {cos(phi)*sin(theta), -sin(phi), cos(phi)*cos(theta)};
    double* qconj = (double*)malloc(4*sizeof(double));
    qconj[0] = qwr[0];
    qconj[1] = qwr[1];
    qconj[2] = qwr[2];
    qconj[3] = qwr[3];
    double* arg2 = (double*)malloc(3*sizeof(double));
    arg2[0] = (yi[0] - rw[0])*lambda + mi[0];
    arg2[1] = (yi[1] - rw[1])*lambda + mi[1];
    arg2[2] = (yi[2] - rw[2])*lambda + mi[2];
    cv::Mat* dRq_times_a_by_dqRES_mat = new cv::Mat(3, 4, CV_64FC1, 0.0);
    
    cv::Mat* aMat = new cv::Mat(3, 1, CV_64FC1, arg2, cv::Mat::AUTO_STEP);
    
    bool success = dRqTimesABydq(qconj, aMat, dRq_times_a_by_dqRES_mat);
    
    double* dqbar_by_dq = (double*)malloc(16*sizeof(double));
    for (int i=0; i<4; ++i){
        for (int j=0; j<4; ++j){
            if (j==i){
                if (i==0)
                    dqbar_by_dq[i*4+j] = 1;
                else
                    dqbar_by_dq[i*4+j] = -1;
            }
            else
                dqbar_by_dq[i*4+j] = 0;
        }
    }
    

    cv::Mat dqbar_mat(4, 4, CV_64FC1, dqbar_by_dq, cv::Mat::AUTO_STEP);
    cv::Mat a_mat = *dRq_times_a_by_dqRES_mat*dqbar_mat;
    double* mxPtr = a_mat.ptr<double>(0);
    
    a[0] = mxPtr[0]; a[1] = mxPtr[1]; a[2] = mxPtr[2]; a[3] = mxPtr[3];
    a[4] = mxPtr[4]; a[5] = mxPtr[5]; a[6] = mxPtr[6]; a[7] = mxPtr[7];
    a[8] = mxPtr[8]; a[9] = mxPtr[9]; a[10] = mxPtr[10]; a[11] = mxPtr[11];
    
    free(arg2);
    dRq_times_a_by_dqRES_mat->release();
    //aMat->release();
    free(dqbar_by_dq);

    dqbar_mat.release();
    a_mat.release();
    
    return success;
}

//called by dh_dxv
bool dh_dqwr(double* Hi2, Camera* cam, double* Xv_km1_k, double* yi, double* zi){
    //Hi2 is a 2x4 matrix
    //dh_dhrl is a 2x3 matrix
    //dhrl_dqwr is a 3x4 matrix
    double* dh_dhrl_data = (double*)malloc(6*sizeof(double));
    double* dhrl_dqwr_data = (double*)malloc(12*sizeof(double));
    
    bool success1 = dh_dhrl(dh_dhrl_data, cam, Xv_km1_k, yi, zi);
    bool success2 = dhrl_dqwr(dhrl_dqwr_data, Xv_km1_k, yi);
    
    cv::Mat dh_dhrl_mat(2, 3, CV_64FC1, dh_dhrl_data, cv::Mat::AUTO_STEP);
    cv::Mat dhrl_dqwr_mat(3, 4, CV_64FC1, dhrl_dqwr_data, cv::Mat::AUTO_STEP);
    cv::Mat ans = dh_dhrl_mat*dhrl_dqwr_mat;
    
    double* mxPtr = ans.ptr<double>(0);
    Hi2[0] = mxPtr[0]; Hi2[1] = mxPtr[1]; Hi2[2] = mxPtr[2]; Hi2[3] = mxPtr[3];
    Hi2[4] = mxPtr[4]; Hi2[5] = mxPtr[5]; Hi2[6] = mxPtr[6]; Hi2[7] = mxPtr[7];
    //Hi2[8] = mxPtr[8]; Hi2[9] = mxPtr[9]; Hi2[10] = mxPtr[10]; Hi2[11] = mxPtr[11];
    
    
    return success1&&success2;
}

/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
//Miscellaneous private functions for dh_dy:

//called by dh_dy:
bool dhrl_dy(double* a, double* Xv_km1_k, double* yi){
    double rw[3] = {Xv_km1_k[0], Xv_km1_k[1], Xv_km1_k[2]};
    double* Xv_4to7 = (double*)malloc(4*sizeof(double));
    Xv_4to7[0] = Xv_km1_k[3];
    Xv_4to7[1] = Xv_km1_k[4];
    Xv_4to7[2] = Xv_km1_k[5];
    Xv_4to7[3] = Xv_km1_k[6];
    double* R = (double*)malloc(9*sizeof(double));
    bool success1 = q2r(Xv_4to7, R);
    cv::Mat Rrw(3, 3, CV_64FC1, R, cv::Mat::AUTO_STEP);
    cv::invert(Rrw, Rrw);
    double lambda = yi[5];
    double phi = yi[4];
    double theta = yi[3];
    
    //redundancies:
    double cosPhi = cos(phi);
    double sinPhi = sin(phi);
    double cosTheta = cos(theta);
    double sinTheta = sin(theta);
    
    double* dmi_dthetai_data = (double*)malloc(3*sizeof(double));
    dmi_dthetai_data[0] = cosPhi*cosTheta;
    dmi_dthetai_data[1] = 0;
    dmi_dthetai_data[2] = -1*cosPhi*sinTheta;
    cv::Mat dmi_dthetai_mat(1, 3, CV_64FC1, dmi_dthetai_data, cv::Mat::AUTO_STEP);
    cv::transpose(dmi_dthetai_mat, dmi_dthetai_mat);
    dmi_dthetai_mat = Rrw*dmi_dthetai_mat;
    
    double* dmi_dphii_data = (double*)malloc(3*sizeof(double));
    dmi_dphii_data[0] = -1*sinPhi*sinTheta;
    dmi_dphii_data[1] = -1*cosPhi;
    dmi_dphii_data[2] = -1*sinPhi*cosTheta;
    cv::Mat dmi_dphii_mat(1, 3, CV_64FC1, dmi_dphii_data, cv::Mat::AUTO_STEP);
    cv::transpose(dmi_dphii_mat, dmi_dphii_mat);
    dmi_dphii_mat = Rrw*dmi_dphii_mat;
    
    double* mxPtr_Rrw = Rrw.ptr<double>(0);
    double* mxPtr_dmi_dthetai = dmi_dthetai_mat.ptr<double>(0);
    double* mxPtr_dmi_dphii = dmi_dphii_mat.ptr<double>(0);
    
    double* yi_1to3_minus_rw = (double*)malloc(3*sizeof(double));
    yi_1to3_minus_rw[0] = yi[0]-rw[0];
    yi_1to3_minus_rw[1] = yi[1]-rw[1];
    yi_1to3_minus_rw[2] = yi[2]-rw[2];
    cv::Mat rrw_times_yi_1to3_minus_rw(3, 1, CV_64FC1, yi_1to3_minus_rw, cv::Mat::AUTO_STEP);
    rrw_times_yi_1to3_minus_rw = Rrw*rrw_times_yi_1to3_minus_rw;
    double* mxPtr_rrw_times_yi = rrw_times_yi_1to3_minus_rw.ptr<double>(0);
    
    /*************first row************/
    a[0] = lambda*mxPtr_Rrw[0];
    a[1] = lambda*mxPtr_Rrw[1];
    a[2] = lambda*mxPtr_Rrw[2];
    
    a[3] = mxPtr_dmi_dthetai[0];
    
    a[4] = mxPtr_dmi_dphii[0];
    
    a[5] = mxPtr_rrw_times_yi[0];
     /*************second row************/
    a[6] = lambda*mxPtr_Rrw[3];
    a[7] = lambda*mxPtr_Rrw[4];
    a[8] = lambda*mxPtr_Rrw[5];
    
    a[9] = mxPtr_dmi_dthetai[1];
    
    a[10] = mxPtr_dmi_dphii[1];
    
    a[11] = mxPtr_rrw_times_yi[1];
    /*************third row************/
    a[12] = lambda*mxPtr_Rrw[6];
    a[13] = lambda*mxPtr_Rrw[7];
    a[14] = lambda*mxPtr_Rrw[8];
    
    a[15] = mxPtr_dmi_dthetai[2];
    
    a[16] = mxPtr_dmi_dphii[2];
    
    a[17] = mxPtr_rrw_times_yi[2];
    
    free(Xv_4to7);
    free(R);
    Rrw.release();
    free(dmi_dthetai_data);
    dmi_dthetai_mat.release();
    free(dmi_dphii_data);
    dmi_dphii_mat.release();
    free(yi_1to3_minus_rw);
    rrw_times_yi_1to3_minus_rw.release();
    
    return success1;
}

/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/

bool dh_dxv(double* Hi_features, Camera* cam, double* Xv_km1_k, double* yi, double* zi){
    //Concatenate two matrices by column:
    //dh_drw returns a 2x3 matrix
    //dh_dqwr returns a 3x4 matrix
    //The last 2x6 entries are just zeros
    double* dh_drw_dat = (double*)malloc(2*3*sizeof(double));
    double* dh_dqwr_dat = (double*)malloc(2*4*sizeof(double));
    
    bool success1 = dh_drw(dh_drw_dat, cam, Xv_km1_k, yi, zi);
    bool success2 = dh_dqwr(dh_dqwr_dat, cam, Xv_km1_k, yi, zi);
    
    //first row
    Hi_features[0] = dh_drw_dat[0];
    Hi_features[1] = dh_drw_dat[1];
    Hi_features[2] = dh_drw_dat[2];
    
    Hi_features[3] = dh_dqwr_dat[0];
    Hi_features[4] = dh_dqwr_dat[1];
    Hi_features[5] = dh_dqwr_dat[2];
    Hi_features[6] = dh_dqwr_dat[3];

    Hi_features[7] = 0.0;
    Hi_features[8] = 0.0;
    Hi_features[9] = 0.0;
    Hi_features[10] = 0.0;
    Hi_features[11] = 0.0;
    Hi_features[12] = 0.0;
    
    //second row
    Hi_features[13] = dh_drw_dat[3];
    Hi_features[14] = dh_drw_dat[4];
    Hi_features[15] = dh_drw_dat[5];
    
    Hi_features[16] = dh_dqwr_dat[4];
    Hi_features[17] = dh_dqwr_dat[5];
    Hi_features[18] = dh_dqwr_dat[6];
    Hi_features[19] = dh_dqwr_dat[7];
    
    Hi_features[20] = 0.0;
    Hi_features[21] = 0.0;
    Hi_features[22] = 0.0;
    Hi_features[23] = 0.0;
    Hi_features[34] = 0.0;
    Hi_features[25] = 0.0;
    
    free(dh_drw_dat);
    free(dh_dqwr_dat);
    return success1&&success2;
}

bool dh_dy(double* Hi_features, Camera* cam, double* Xv_km1_k, double* yi, double* zi){
    //dh_dhrl produces a 2x3 matrix
    //dhrl_dy produces a 3x6 matrix
    //Hi_features is a 2x6 matrix
    
    double* dh_dhrl_data = (double*)malloc(6*sizeof(double));
    double* dhrl_dy_data = (double*)malloc(18*sizeof(double));
    bool success1 = dh_dhrl(dh_dhrl_data, cam, Xv_km1_k, yi, zi);
    bool success2 = dhrl_dy(dhrl_dy_data, Xv_km1_k, yi);
    
    cv::Mat dh_dhrl_mat(2, 3, CV_64FC1, dh_dhrl_data, cv::Mat::AUTO_STEP);
    cv::Mat dhrl_dy_mat(3, 6, CV_64FC1, dhrl_dy_data, cv::Mat::AUTO_STEP);
    cv::Mat Hii = dh_dhrl_mat*dhrl_dy_mat;
    
    double* mxPtr = Hii.ptr<double>(0);
    
    Hi_features[0] = mxPtr[0];
    Hi_features[1] = mxPtr[1];
    Hi_features[2] = mxPtr[2];
    Hi_features[3] = mxPtr[3];
    Hi_features[4] = mxPtr[4];
    Hi_features[5] = mxPtr[5];
    Hi_features[6] = mxPtr[6];
    Hi_features[7] = mxPtr[7];
    Hi_features[8] = mxPtr[8];
    Hi_features[9] = mxPtr[9];
    Hi_features[10] = mxPtr[10];
    Hi_features[11] = mxPtr[11];
    
    free(dh_dhrl_data);
    free(dhrl_dy_data);
    dh_dhrl_mat.release();
    dhrl_dy_mat.release();
    Hii.release();
    
    return success1&&success2;
}


/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/


bool inverseDepth2Cartesian(double* inverse_depth, double* cartesian)
{
    double *rw = (double *)malloc(1*3*sizeof(double));
    rw[0] = *(inverse_depth+0);
    rw[1] = *(inverse_depth+1);
    rw[2] = *(inverse_depth+2);
    
    double theta = *(inverse_depth+3);
    double phi = *(inverse_depth+4);
    double rho = *(inverse_depth+5);
        
    double *m = (double *)malloc(1*3*sizeof(double));
    m[0] = cos(phi)*sin(theta);
    m[1] = (-1)*sin(phi);
    m[2] = cos(phi)*cos(theta);
        
    cartesian[0] = rw[0] + 1/rho*m[0];
    cartesian[1] = rw[1] + 1/rho*m[1];
    cartesian[2] = rw[2] + 1/rho*m[2];
        
    return true;
}

/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/


bool predict_features_appearance(NSMutableArray* features_info, Filter* filter, Camera* cam){
    
    if(features_info == NULL || filter == NULL || cam == NULL )
        return FALSE;
    
    
    double* r_wc = (double*)malloc(3*sizeof(double));
    r_wc[0] = filter->x_k_km1[0];
    r_wc[1] = filter->x_k_km1[1];
    r_wc[2] = filter->x_k_km1[2];
    
    double* R_wc = (double*)malloc(3*3*sizeof(double));
    
    /* Instead of writing the function q2r, the following code performs the same as the function, hence written here
     */
    
    double x = filter->x_k_km1[4];
    double y = filter->x_k_km1[5];
    double z = filter->x_k_km1[6];
    double r = filter->x_k_km1[3];
    
    R_wc[0] = ((r*r)+(x*x)-(y*y) - (z*z));
    R_wc[1] = 2*((x*y) - (r*z));
    R_wc[2] = 2*((z*x)+(r*y));
    R_wc[3] = 2*((x*y)+(r*z));
    R_wc[4] = ((r*r)-(x*x)+(y*y)-(z*z));
    R_wc[5] = 2*((y*z)-(r*x));
    R_wc[6] = 2*((z*x)-(r*y));
    R_wc[7] = 2*((y*z)+(r*x));
    R_wc[8] = ((r*r)-(x*x)-(y*y)+(z*z));
    
    /* End of the code which is originally written for q2r function */
    
    //double* x_k_k_rest_of_features_0 = (double*)malloc((filter->state_size - 13)*sizeof(double));
    
    int yStartIdx = 13;
    
    /*int j;
    for ( j = 13; j < filter->state_size; j++) {
        x_k_k_rest_of_features_0[j-13] = filter->x_k_k[j]; // Size is (state_size - 13)
        xResidual[j-13] = filter->x_k_k[j];
    }*/
    
    double* yDat = (double*)malloc(6*sizeof(double));
    
    NSLog(@"r_wc %lf %lf %lf", r_wc[0], r_wc[1], r_wc[2]);
    NSLog(@"R_wc %lf %lf %lf %lf", R_wc[0], R_wc[1], R_wc[2]);
    
    for (int i=0; i<(int)[features_info count]; i++)
    {
        Feature *current=features_info[i];
        double* XYZ_w = (double*)malloc(3*sizeof(double));
        
        if ([current->type isEqualToString:@"cartesian"])
        {
            //TODO: Implement this for cartesian
            /*XYZ_w[0] = x_k_k_rest_of_features_0[0];
            XYZ_w[1] = x_k_k_rest_of_features_0[1];
            XYZ_w[2] = x_k_k_rest_of_features_0[2];*/
            //double* x_k_k_rest_of_features_1 = (double*)malloc(1*(filter->state_size-16)*sizeof(double));
            //for (int y = 3; y < (filter->state_size-13); y++) {
               // x_k_k_rest_of_features_1[y-3] = x_k_k_rest_of_features_0[y];
            //}
            
        }
        
        if ([current->type isEqualToString:@"inversedepth"]) {
            //double* yDat = (double*)malloc(6*sizeof(double));
            /*yDat[0] = x_k_k_rest_of_features_0[0];
            yDat[1] = x_k_k_rest_of_features_0[1];
            yDat[2] = x_k_k_rest_of_features_0[2];
            yDat[3] = x_k_k_rest_of_features_0[3];
            yDat[4] = x_k_k_rest_of_features_0[4];
            yDat[5] = x_k_k_rest_of_features_0[5];*/
           
            yDat[0] = filter->x_k_km1[yStartIdx];
            yDat[1] = filter->x_k_km1[yStartIdx+1];
            yDat[2] = filter->x_k_km1[yStartIdx+2];
            yDat[3] = filter->x_k_km1[yStartIdx+3];
            yDat[4] = filter->x_k_km1[yStartIdx+4];
            yDat[5] = filter->x_k_km1[yStartIdx+5];
            
            //double* x_k_k_rest_of_features = (double*)malloc(1*(filter->state_size-19)*sizeof(double));
            yStartIdx+=6;
            
            //for (int h = 6; h < (filter->state_size-13); h++) {
               // x_k_k_rest_of_features[h-6] = x_k_k_rest_of_features_0[h];
            //}
            inverseDepth2Cartesian(yDat, XYZ_w);
        }
        
        if(current->h != NULL) {
            pred_patch_fc(cam, current, R_wc, r_wc, XYZ_w, current->patch_when_matching);
        }
        
        free(XYZ_w);
    }
    
    
    NSLog(@"r_wc %lf %lf %lf", r_wc[0], r_wc[1], r_wc[2]);
    NSLog(@"R_wc %lf %lf %lf %lf", R_wc[0], R_wc[1], R_wc[2]);
      
    
    
    free(r_wc);
    free(R_wc);
    //free(x_k_k_rest_of_features_0);
    
    
    
    return true;
    
}

bool pred_patch_fc(Camera* cam, Feature* patch_p_f_ini, double* R_Wk, double* r_Wk, double* XYZ_w, cv::Mat* patch_pred)
{
    if(cam == NULL || patch_p_f_ini == NULL || R_Wk == NULL || r_Wk == NULL || XYZ_w == NULL || patch_pred == NULL)
        return FALSE;
    
    // Actual code for the function begins here
    double* uv_p_pred = (double*)malloc(2*sizeof(double));
    uv_p_pred[0] = patch_p_f_ini->h[0];
    uv_p_pred[1] = patch_p_f_ini->h[1];
    int halfW_pred = patch_p_f_ini->half_patch_size_when_matching;
    
    if(((uv_p_pred[0] > halfW_pred) && (uv_p_pred[0] < (cam->nCols - halfW_pred)) && (uv_p_pred[1] > halfW_pred) && (uv_p_pred[1] < (cam->nRows - halfW_pred)))) {
        
        double* uv_p_f = (double*)malloc(2*sizeof(double));
        
        uv_p_f[0] = patch_p_f_ini->uv_when_initialized[0]; // 1 x 2 Matrix
        uv_p_f[1] = patch_p_f_ini->uv_when_initialized[1];
        double* R_Wk_p_f = (double*)malloc(3*3*sizeof(double));
        for (int u = 0; u < 9; u++)
            R_Wk_p_f[u] = patch_p_f_ini->R_wc_when_initialized[u];
        double* r_Wk_p_f = (double*)malloc(3*sizeof(double));
        r_Wk_p_f[0] = patch_p_f_ini->r_wc_when_initialized[0];
        r_Wk_p_f[1] = patch_p_f_ini->r_wc_when_initialized[1];
        r_Wk_p_f[2] = patch_p_f_ini->r_wc_when_initialized[2];
        
        
        /* THE FOLLOWING CODE HAS TO BE WRITTEN HERE. BELOW ARE THE LINES */
         
        
        
        double* H_Wk_p_f = (double*)malloc(4*4*sizeof(double));
        /* Writing the code in order to execute the matrix multiplication as shown in MATLAB code for this matrix*/
        
        H_Wk_p_f[0] = R_Wk_p_f[0];
        H_Wk_p_f[1] = R_Wk_p_f[1];
        H_Wk_p_f[2] = R_Wk_p_f[2];
        H_Wk_p_f[3] = R_Wk_p_f[0]*r_Wk_p_f[0] + R_Wk_p_f[1]*r_Wk_p_f[1] + R_Wk_p_f[2]*r_Wk_p_f[2];
        H_Wk_p_f[4] = R_Wk_p_f[3];
        H_Wk_p_f[5] = R_Wk_p_f[4];
        H_Wk_p_f[6] = R_Wk_p_f[5];
        H_Wk_p_f[7] = R_Wk_p_f[3]*r_Wk_p_f[0] + R_Wk_p_f[4]*r_Wk_p_f[1] + R_Wk_p_f[5]*r_Wk_p_f[2];
        H_Wk_p_f[8] = R_Wk_p_f[6];
        H_Wk_p_f[9] = R_Wk_p_f[7];
        H_Wk_p_f[10] = R_Wk_p_f[8];
        H_Wk_p_f[11] = R_Wk_p_f[6]*r_Wk_p_f[0] + R_Wk_p_f[7]*r_Wk_p_f[1] + R_Wk_p_f[8]*r_Wk_p_f[2];
        H_Wk_p_f[12] = H_Wk_p_f[13] = H_Wk_p_f[14] = 0;
        H_Wk_p_f[15] = 1;
        
        /* End of the Matrix Multiplication Result
         
         /* Matrix Multiplication code for H_Wk */
        
        double* H_Wk = (double*)malloc(4*4*sizeof(double));
        
        H_Wk[0] = R_Wk[0];
        H_Wk[1] = R_Wk[1];
        H_Wk[2] = R_Wk[2];
        H_Wk[3] = R_Wk[0]*r_Wk[0] + R_Wk[1]*r_Wk[1] + R_Wk[2]*r_Wk[2];
        H_Wk[4] = R_Wk[3];
        H_Wk[5] = R_Wk[4];
        H_Wk[6] = R_Wk[5];
        H_Wk[7] = R_Wk[3]*r_Wk[0] + R_Wk[4]*r_Wk[1] + R_Wk[5]*r_Wk[2];
        H_Wk[8] = R_Wk[6];
        H_Wk[9] = R_Wk[7];
        H_Wk[10] = R_Wk[8];
        H_Wk[11] = R_Wk[6]*r_Wk[0] + R_Wk[7]*r_Wk[1] + R_Wk[8]*r_Wk[2];
        H_Wk[12] = H_Wk[13] = H_Wk[14] = 0;
        H_Wk[15] = 1;
        
        /* End of Matrix Multiplication equivalent code */
        
        // Converting the matrices to OpenCV format so that Inverse and Matrix MUltiplication can be carried out
        
         cv::Mat H_Wk_p_f_mat(4, 4, CV_64FC1, H_Wk_p_f, cv::Mat::AUTO_STEP);
         cv::Mat H_Wk_mat(4, 4, CV_64FC1, H_Wk, cv::Mat::AUTO_STEP);
        
         cv::Mat H_Wk_p_f_inverse;
         cv::invert(H_Wk_p_f_mat, H_Wk_p_f_inverse);
        
        // Matrix Multiplication is carried out in OpenCV and rewritten back
        
        cv::Mat H_kpf_k_mat = H_Wk_p_f_inverse*H_Wk_mat;
        
        int H_kpf_k_size = H_kpf_k_mat.rows*H_kpf_k_mat.cols;
        
        // Declaration of the original matrix H_kpf_k in double* format to store the data back to C array
        
        double* H_kpf_k = (double*)malloc(H_kpf_k_size*sizeof(double));
        
        //Create a pointer to data in Q_mat:
        double* mxPtr = H_kpf_k_mat.ptr<double>(0);
        //Iterate through and assign mat values to our data:
        for (int i=0; i<(H_kpf_k_size); ++i)
            H_kpf_k[i] = mxPtr[i];
        
         cv::Mat patch_p_f = patch_p_f_ini->patch_when_initialized->clone();
        
        int halfW_fea = patch_p_f_ini->half_patch_size_when_initialized;
        double d = cam->dx;
        double f = cam->f;
        
        double* n1 = (double*)malloc(3*sizeof(double));
        n1[0] = -((-1)*(uv_p_f[0] - cam->Cx));
        n1[1] = -((-1)*(uv_p_f[1] - cam->Cy));
        n1[2] = (-1)*(f/d);
        
        double* n2_init = (double*)malloc(3*sizeof(double));
        n2_init[0] = -((-1)*(uv_p_pred[0] - cam->Cx));
        n2_init[1] = -((-1)*(uv_p_pred[1] - cam->Cy));
        n2_init[2] = (-1)*(f/d);
        
        // Creating a matrix n2 modified so that it can have a fourth value as 1
        
         double* n2_init_mod = (double*)malloc(4*sizeof(double));
        n2_init_mod[0] = n2_init[0];
        n2_init_mod[1] = n2_init[1];
        n2_init_mod[2] = n2_init[2];
        n2_init_mod[3] = 1;
        // Declaring n2_mod as an openCV matrix so that matrix multiplication is carried out easily
        
        cv::Mat n2_mod_mat(4, 1, CV_64FC1, n2_init_mod, cv::Mat::AUTO_STEP);
        //cv::Mat n2_mod_transpose;
        //cv::transpose(n2_mod_mat, n2_mod_transpose);
        
        cv::Mat n2_final = H_kpf_k_mat*n2_mod_mat;
        
        double* n2 = (double*)malloc(4*sizeof(double));
        int n2_final_size = n2_final.rows*n2_final.cols;
        
        double* pointer_n2_cvmat = n2_final.ptr<double>(0);
        for (int j = 0; j <(n2_final_size); ++j)
            n2[j] = pointer_n2_cvmat[j];
        
        double* n2_intermediate_0 = (double*)malloc(4*sizeof(double));
        n2_intermediate_0[0] = n2[0]/n2[3];
        n2_intermediate_0[1] = n2[1]/n2[3];
        n2_intermediate_0[2] = n2[2]/n2[3];
        n2_intermediate_0[3] = 1;
        
        double* n2_intermediate_1 = (double*)malloc(3*sizeof(double));
        n2_intermediate_1[0] = n2_intermediate_0[0];
        n2_intermediate_1[1] = n2_intermediate_0[1];
        n2_intermediate_1[2] = n2_intermediate_0[2];
        
        double* n1_intermediate = (double*)malloc(3*sizeof(double));
        double n1_norm = sqrt((n1[0]*n1[0])+(n1[1]*n1[1])+(n1[2]*n1[2]));
        double* n2_final_matrix = (double*)malloc(3*sizeof(double));
        double n2_norm = sqrt((n2_intermediate_1[0]*n2_intermediate_1[0])+(n2_intermediate_1[1]*n2_intermediate_1[1])+ (n2_intermediate_1[2]*n2_intermediate_1[2]));
        
        n1_intermediate[0] = n1[0]/n1_norm;
        n1_intermediate[1] = n1[1]/n1_norm;
        n1_intermediate[2] = n1[2]/n1_norm;
        
        n2_final_matrix[0] = n2_intermediate_1[0]/n2_norm;
        n2_final_matrix[1] = n2_intermediate_1[1]/n2_norm;
        n2_final_matrix[2] = n2_intermediate_1[2]/n2_norm;
        
        double* n_intermediate = (double*)malloc(3*sizeof(double));
        n_intermediate[0] = n1_intermediate[0] + n2_final_matrix[0];
        n_intermediate[1] = n1_intermediate[1] + n2_final_matrix[1];
        n_intermediate[2] = n1_intermediate[2] + n2_final_matrix[2];
        
        double n_intermediate_norm = sqrt((n_intermediate[0]*n_intermediate[0])+(n_intermediate[1]*n_intermediate[1])+(n_intermediate[2]*n_intermediate[2]));
        
       double* n = (double*)malloc(3*sizeof(double));
        n[0] = n_intermediate[0]/n_intermediate_norm;
        n[1] = n_intermediate[1]/n_intermediate_norm;
        n[2] = n_intermediate[2]/n_intermediate_norm;
        
        // Declaration for the modified XYZ_w matrix with an extra value of 1 making the size as 4
        
        double* XYZ_w_mod = (double*)malloc(4*sizeof(double));
        XYZ_w_mod[0] = XYZ_w[0];
        XYZ_w_mod[1] = XYZ_w[1];
        XYZ_w_mod[2] = XYZ_w[2];
        XYZ_w_mod[3] = 1;
        
         cv::Mat XYZ_w_mod_mat(4, 1, CV_64FC1, XYZ_w_mod, cv::Mat::AUTO_STEP);
         cv::Mat XYZ_kpf_mat = H_Wk_p_f_inverse*XYZ_w_mod_mat;
        
        double* XYZ_kpf = (double*)malloc(4*sizeof(double));
        int XYZ_kpf_mat_final_size = XYZ_kpf_mat.rows*XYZ_kpf_mat.cols;
        
        double* pointer_XYZ_kpf_cvmat = XYZ_kpf_mat.ptr<double>(0);
        for (int k = 0; k <(XYZ_kpf_mat_final_size); ++k)
            XYZ_kpf[k] = pointer_XYZ_kpf_cvmat[k];
        
        XYZ_kpf[0] = XYZ_kpf[0]/XYZ_kpf[3];
        XYZ_kpf[1] = XYZ_kpf[1]/XYZ_kpf[3];
        XYZ_kpf[2] = XYZ_kpf[2]/XYZ_kpf[3];
        XYZ_kpf[3] = 1;
        
        d = (-1)*((n[0]*XYZ_kpf[0])+(n[1]*XYZ_kpf[1])+(n[2]*XYZ_kpf[2]));
        
        double* uv_p_pred_patch = (double*)malloc(2*sizeof(double));
        
        uv_p_pred_patch[0] = uv_p_pred_patch[1] = 0; // Initializing the matrix to all zero
        
        double* H_kpf_k_mod0 = (double*)malloc(3*3*sizeof(double));
        H_kpf_k_mod0[0] = H_kpf_k[0];
        H_kpf_k_mod0[1] = H_kpf_k[1];
        H_kpf_k_mod0[2] = H_kpf_k[2];
        H_kpf_k_mod0[3] = H_kpf_k[4];
        H_kpf_k_mod0[4] = H_kpf_k[5];
        H_kpf_k_mod0[5] = H_kpf_k[6];
        H_kpf_k_mod0[6] = H_kpf_k[8];
        H_kpf_k_mod0[7] = H_kpf_k[9];
        H_kpf_k_mod0[8] = H_kpf_k[10];
        
        double* H_kpf_k_mod1 = (double*)malloc(3*sizeof(double));
        H_kpf_k_mod1[0] = H_kpf_k[3];
        H_kpf_k_mod1[1] = H_kpf_k[7];
        H_kpf_k_mod1[2] = H_kpf_k[11];
        
        rotate_with_dist_fc_c2c1(cam, uv_p_f, H_kpf_k_mod0, H_kpf_k_mod1, n, d, uv_p_pred_patch); // Returns uv_p_pred_patch as result which is a 1 x 2 matrix

        /* Code to create a meshgrid. This code replicates the matrix formation for u_pred and v_pred matrices */
        
        int numrows = (uv_p_pred_patch[1] + halfW_pred) - (uv_p_pred_patch[1]- halfW_pred) + 1;
        int numcols = (uv_p_pred_patch[0] + halfW_pred) - (uv_p_pred_patch[0] - halfW_pred) + 1;
        
        //Declaration for the matrices u_pred and v_pred with the obtained number of columns and rows. Each of these matrices has equal number of rows and columns respectively
        
        double* u_pred = (double*)malloc(numrows*numcols*sizeof(double));
        double* v_pred = (double*)malloc(numrows*numcols*sizeof(double));
        
        int a,b;
        double startvalue_u_matrix = uv_p_pred_patch[0] - halfW_pred;
        double startvalue_v_matrix = uv_p_pred_patch[1] - halfW_pred;
        for (a = 0; a < numrows; a++) {
            for (b =0; b < numcols; b++) {
                int index =  ((numcols*a)+b);
                if(b == 0) {
                     startvalue_u_matrix = uv_p_pred_patch[0] - halfW_pred;
                      u_pred[index] = startvalue_u_matrix;
                }
                else {
                    u_pred[index] = startvalue_u_matrix;
                }
                startvalue_u_matrix = startvalue_u_matrix+1.0;;
            }
        }
        for (a = 0; a<numrows; ++a){
            startvalue_v_matrix = uv_p_pred_patch[1] - halfW_pred + 1.0*a;
        
            for (b=0; b<numcols; ++b){
                int idx = ((numcols*a)+b);
                v_pred[idx] = startvalue_v_matrix;
            }
        }
        int nr_uv_pred = ((halfW_pred*2))+1;
        int nrows_uv_pred = (nr_uv_pred)*(nr_uv_pred);
        
        double* uv_pred = (double*)malloc(nrows_uv_pred*2*sizeof(double));
        int y;
        //NSLog(@"******************************UVPRED*****************************");
        for (y = 0; y < nrows_uv_pred; y++)
        {
            uv_pred[(2*y)+0] = u_pred[y];
            uv_pred[(2*y)+1] = v_pred[y];
            //NSLog(@"%lf %lf", uv_pred[(2*y)+0], uv_pred[(2*y)+1]);
        }
        
        double* uv_pred_imak_dist = (double*)malloc(nrows_uv_pred*2*sizeof(double));
        
        
        rotate_with_dist_fc_c1c2(cam, uv_pred, nrows_uv_pred, H_kpf_k_mod0, H_kpf_k_mod1, n, d, uv_pred_imak_dist);
        
        
        
        //Selecting the first and second columns of the uv_pred_imak_dist matrix and doing algebraic manipulations
        
        int r;
        //NSLog(@"*****************UVPREDIMAK***********************");
        for (r = 0; r < nrows_uv_pred; r++) {
            uv_pred_imak_dist[(2*r)+0] = uv_pred_imak_dist[(2*r)+0] - (uv_p_f[0] - halfW_fea  - 1);
            uv_pred_imak_dist[(2*r)+1] = uv_pred_imak_dist[(2*r)+1] - (uv_p_f[1] - halfW_fea  - 1);
            //NSLog(@"%lf %lf", uv_pred_imak_dist[(2*r)+0], uv_pred_imak_dist[(2*r)+1] );
        }
    

        
        
        // Declaration of the matrices u_pred_imak_dist and v_pred_imak_dist which store the column 1 and 2 values of uv_pred_imak_dist matrix
        
        double* u_pred_imak_dist = (double*)malloc(nr_uv_pred*nr_uv_pred*sizeof(double));
        double* v_pred_imak_dist = (double*)malloc(nr_uv_pred*nr_uv_pred*sizeof(double));
        
        int h,g;
        int uv_counter = 0;
        for(h = 0; h <13; h++) {
            for(g = 0; g <13; g++) {
                u_pred_imak_dist[(13*h)+g] = uv_pred_imak_dist[uv_counter];
                ++uv_counter;
                v_pred_imak_dist[(13*h)+g] = uv_pred_imak_dist[uv_counter];
                ++uv_counter;
            }
        }
        

       
        double* patch_pred_data = (double*)malloc(13*13*sizeof(double));
        //MATLAB CODE FOR INTERPOLATOR
        //patch_pred=interp2(u_fea,v_fea,double(patch_p_f),u_pred_imak_dist,v_pred_imak_dist);
        interpPatchPred(patch_p_f.ptr<uchar>(0), u_pred_imak_dist, v_pred_imak_dist, patch_pred_data);
        
        double* mxPtrPP = patch_pred->ptr<double>(0);
        for (int i=0; i<169; ++i){
            mxPtrPP[i] = patch_pred_data[i];
        }
        //patch_pred = new cv::Mat(13, 13, CV_8UC1, patch_pred_data, cv::Mat::AUTO_STEP);
        
        free(n1);
        free(n1_intermediate);
        free(n2_init);
        free(n2_final_matrix);
        free(n2_init_mod);
        free(n2_intermediate_0);
        free(n2_intermediate_1);
        free(n);
        free(XYZ_kpf);
        //free(XYZ_w);
        free(XYZ_w_mod);
        free(H_kpf_k);
        free(H_kpf_k_mod0);
        free(H_kpf_k_mod1);
        //free(H_Wk);
        free(H_Wk_p_f);
    }
    else {
        // If the condition fails, then make patch_pred as all zero matrix
        patch_pred->release();
        patch_pred = new cv::Mat(13, 13, CV_8UC1, 0);
    }
    
    
    
    return true;
}


/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/

bool rotate_with_dist_fc_c2c1(Camera* cam, double* uv_c1, double* R_c2c1, double* t_c2c1, double* n, double d, double* uv_c2){
    cv::Mat *K_mat = new cv::Mat(3,3,CV_64FC1,cam->K, cv::Mat::AUTO_STEP);
    
    cv::Mat *invK = new cv::Mat(3,3,CV_64FC1,0.0);
    
    cv::invert(*K_mat, *invK);
    
    
    
    double *uv_c1_und = (double *)malloc(2*sizeof(double));

    undistortFM(uv_c1, 1, cam, uv_c1_und);
    
    
    
    cv::Mat *R_c2c1_mat = new cv::Mat(3,3,CV_64FC1,R_c2c1, cv::Mat::AUTO_STEP);
    
    cv::Mat *t_c2c1_mat = new cv::Mat(3,1,CV_64FC1,t_c2c1, cv::Mat::AUTO_STEP);
    
    cv::Mat *n_mat = new cv::Mat(1,3,CV_64FC1,n, cv::Mat::AUTO_STEP);
    
    
    
    cv::Mat *A = new cv::Mat(3,1,CV_64FC1,0.0);
    
    A->at<double>(0,0) = *(uv_c1_und+0);
    
    A->at<double>(1,0) = *(uv_c1_und+1);
    
    A->at<double>(2,0) = 1;
    
    
    
    cv::Mat *B = new cv::Mat(3,3,CV_64FC1,0.0);
    
    *B = (*K_mat) * ((*R_c2c1_mat)-(*t_c2c1_mat)*(*n_mat)/d) * (*invK);
    
    cv::Mat *invB = new cv::Mat(3,3,CV_64FC1,0.0);
    
    cv::invert(*B, *invB);
    
    
    
    cv::Mat *uv_c2_und = new cv::Mat(3,1,CV_64FC1,0.0);
    
    *uv_c2_und = (*invB) * (*A);
    
    
    
    double *uv_c2_0 = (double *)malloc(1*2*sizeof(double));
    
    uv_c2_0[0] = uv_c2_und->at<double>(0,0)/uv_c2_und->at<double>(2,0);
    
    uv_c2_0[1] = uv_c2_und->at<double>(1,0)/uv_c2_und->at<double>(2,0);
    
    distortFM(uv_c2, uv_c2_0, 1, cam);
    
    K_mat->release();
    invK->release();
    free(uv_c1_und);
    R_c2c1_mat->release();
    n_mat->release();
    A->release();
    B->release();
    invB->release();
    free(uv_c2_0);
    
    return true;
}


/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/

bool rotate_with_dist_fc_c1c2(Camera* cam, double* uv_c2, int uv_c2_length, double* R_c1c2, double* t_c1c2, double* n, double d, double* uv_c1)
{
    cv::Mat K_mat(3,3,CV_64FC1,cam->K, cv::Mat::AUTO_STEP);
    cv::Mat invK(3,3,CV_64FC1,0.0);
    cv::invert(K_mat, invK);
    double *uv_c2_und = (double *)malloc(2*uv_c2_length*sizeof(double));
    cv::Mat uv_c2_transpose(uv_c2_length, 2, CV_64FC1, uv_c2, cv::Mat::AUTO_STEP);
    cv::transpose(uv_c2_transpose, uv_c2_transpose);
    undistortFM(uv_c2_transpose.ptr<double>(0), uv_c2_length, cam, uv_c2_und);
    //uv_c2_und is 2x169
    //cv::Mat uv_c2_und_actual(2, 169, CV_64FC1, uv_c2_und);
    //cv::transpose(uv_c2_und_actual, uv_c2_und_actual);
    
    //double* undPtr = uv_c2_und_actual.ptr<double>(0);
    //for (int i=0; i<2*uv_c2_length; ++i){
    //    uv_c2_und[i] = undPtr[0];
   // }
    
    cv::Mat R_c1c2_mat(3,3,CV_64FC1,R_c1c2, cv::Mat::AUTO_STEP);
    cv::Mat t_c1c2_mat(3,1,CV_64FC1,t_c1c2, cv::Mat::AUTO_STEP);
    cv::Mat n_transpose_mat(1,3,CV_64FC1,n, cv::Mat::AUTO_STEP);
    
    //This forms a lengthx3 matrix
    double* uv_c2_3_data = (double*)malloc(3*uv_c2_length*sizeof(double));
    for (int i=0; i<uv_c2_length; ++i){
        uv_c2_3_data[i] = uv_c2_und[i];
        uv_c2_3_data[uv_c2_length+i] = uv_c2_und[uv_c2_length+i];
        uv_c2_3_data[2*uv_c2_length+i] = 1.0;
    }
    cv::Mat uv_c1_und_lastTerm(3, 169,CV_64FC1,uv_c2_3_data, cv::Mat::AUTO_STEP);
    //turn it into a 3xlength matrix
    //cv::transpose(uv_c1_und_lastTerm, uv_c1_und_lastTerm);
    
    cv::Mat uv_c1_und(3,uv_c2_length,CV_64FC1,0.0);
    uv_c1_und = K_mat * ((R_c1c2_mat)-(t_c1c2_mat)*(n_transpose_mat)/d) * (invK) * uv_c1_und_lastTerm;
    
    double* uv_c1_und_transpose = (double*)malloc(uv_c2_length*2*sizeof(double));
    double* uv_c1_und_ptr = uv_c1_und.ptr<double>(0);
    for (int i=0; i<uv_c2_length; ++i){
        uv_c1_und_transpose[i] = uv_c1_und_ptr[i]/uv_c1_und_ptr[2*uv_c2_length+i];
        uv_c1_und_transpose[uv_c2_length+i] = uv_c1_und_ptr[uv_c2_length+i]/uv_c1_und_ptr[2*uv_c2_length+i];
       // NSLog(@"%lf %lf %lf", uv_c1_und_ptr[i], uv_c1_und_ptr[uv_c2_length+i], uv_c1_und_ptr[2*uv_c2_length+i]);
    }
    //
    uv_c1_und = cv::Mat(1, uv_c2_length, CV_64FC1, uv_c1_und_transpose, cv::Mat::AUTO_STEP);
    cv::transpose(uv_c1_und, uv_c1_und);
    uv_c1_und_ptr = uv_c1_und.ptr<double>(0);
    
    distortFM(uv_c1, uv_c1_und_transpose, uv_c2_length, cam);
   
    cv::Mat uv_c1_mat(2, uv_c2_length, CV_64FC1, uv_c1, cv::Mat::AUTO_STEP);
    cv::transpose(uv_c1_mat, uv_c1_mat);
    uv_c1_und_ptr = uv_c1_mat.ptr<double>(0);
    //NSLog(@"*****************UV_C1***********************");
    for (int i=0; i<2*uv_c2_length; ++i){
        uv_c1[i] = uv_c1_und_ptr[i];
       // NSLog(@"%lf ", uv_c1[i]);
    }
    
    K_mat.release();
    invK.release();
    free(uv_c2_und);
    R_c1c2_mat.release();
    t_c1c2_mat.release();
    n_transpose_mat.release();
    uv_c1_und.release();
    uv_c1_und_lastTerm.release();
    
    return true;
}

/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/

bool interpPatchPred(uchar* patch_p_f, double* Xq, double* Yq, double* patch_pred){
    double xSamplePt = 0.0;
    double ySamplePt = 0.0;
    int xSamplePtFloor = 0;
    int ySamplePtFloor = 0;
    double yMaxVal = 0.0;
    double xMaxVal = 0.0;
    double maxVal = 0.0;
    double minVal = 0.0;
    double deltaX = 0.0;
    double deltaY = 0.0;
    double deltaXY = 0.0;
    double distX = 0.0;
    double distY = 0.0;
    double distXY = 0.0;
    double interpDeltaX = 0.0;
    double interpDeltaY = 0.0;
    double interpDeltaXY = 0.0;
    for (int r=0; r<13; ++r){
        for (int c=0; c<13; ++c){
            xSamplePt = Xq[r*13+c]-1;
            xSamplePtFloor = (int)xSamplePt;
            ySamplePt = Yq[r*13+c]-1;
            ySamplePtFloor = (int)ySamplePt;
            minVal = (double)patch_p_f[ySamplePtFloor*41+xSamplePtFloor];
            xMaxVal = (double)patch_p_f[ySamplePtFloor*41+(xSamplePtFloor+1)];
            yMaxVal = (double)patch_p_f[(ySamplePtFloor+1)*41+xSamplePtFloor];
            maxVal = (double)patch_p_f[(ySamplePtFloor+1)*41+(xSamplePtFloor+1)];
            //By default, these deltas are distance/unit...in the original grid, unit is 1, so no need to divide:
            deltaX = xMaxVal-minVal;
            deltaY = yMaxVal-minVal;
            deltaXY = maxVal-minVal;
            distX = xSamplePt-(double)xSamplePtFloor;
            distY = ySamplePt-(double)ySamplePtFloor;
            distXY = sqrt(distX*distX+distY*distY);
            interpDeltaX = deltaX*distX;
            interpDeltaY = deltaY*distY;
            interpDeltaXY = deltaXY*distXY;
            
            //The following is an accumulation of strategies for calculating the predicted patch (the one which generally works best is used)
            
            //patch_pred[r*13+c] = minVal+ sqrt(interpDeltaX*interpDeltaX + interpDeltaY*interpDeltaY);
            //patch_pred[r*13+c] = minVal+interpDeltaXY;
            patch_pred[r*13+c] = minVal+interpDeltaX+interpDeltaY;
            
            //Current best strategy
            /*if (deltaX>deltaY)
                patch_pred[r*13+c] = minVal+(interpDeltaXY+interpDeltaX)/2;
            else
                patch_pred[r*13+c] = minVal+(interpDeltaXY+interpDeltaY)/2;
            */
            //Take average of X/Ycomponent and XYabsolute values
            //patch_pred[r*13+c] = minVal+ (interpDeltaX+interpDeltaY + interpDeltaXY)/2;
            
            /*********Least-squares regression impl**************/
            /*
            double bDat[6] = {1,0,0,1,1,1};
            cv::Mat B(3,2, CV_64FC1, bDat, cv::Mat::AUTO_STEP);
            cv::Mat Btranspose;
            cv::transpose(B, Btranspose);
            cv::Mat B0Inv;
            cv::invert(Btranspose*B, B0Inv);
            //B is now Moore-Penrose pseudoinverse:
            B = B0Inv*Btranspose;
            double cDat[3] = {deltaX, deltaY, deltaXY};
            cv::Mat cMat(3,1, CV_64FC1, cDat, cv::Mat::AUTO_STEP);
            cv::Mat zMat = B*cMat;
            double* zPtr = zMat.ptr<double>(0);
            patch_pred[r*13+c] = minVal + deltaX*zPtr[0] + deltaY*zPtr[1];*/
        }
    }
    return true;
}



/***************************************************************************************************************************************/
/***************************************************************************************************************************************/
/***************************************************************************************************************************************/


bool matching(cv::Mat* im, NSMutableArray* features_info, Camera* cam){
    int i,j,k,ii,jj;
    
    float correlation_threshold = 0.80;
    float chi_095_2 = 5.9915;
    
    for (i=0; i<(int)[features_info count]; i++){
        Feature *current = (Feature *)features_info[i];
        
        if (current->h != NULL ){
            NSLog(@"Evaluating h for index %d",i);
            double *h = current->h;
            cv::Mat *S = new cv::Mat(2,2,CV_64FC1,current->S, cv::Mat::AUTO_STEP);
            
            int half_patch_size_when_matching = current->half_patch_size_when_matching;
            int pixels_in_the_matching_patch = pow(2*half_patch_size_when_matching+1,2);
            
            cv::Mat *eigS = new cv::Mat(2,1,CV_64FC1,0.0);
            eigen(*S,*eigS);

            //Make sure this works (Matlab just treats as a scalar)
            if (eigS->at<double>(0,0)<100 && eigS->at<double>(1,0)<100){
                NSLog(@"Eig for index %d within limits", i);
                cv::Mat *invS = new cv::Mat(2,2,CV_64FC1,0.0);
                invert(*S,*invS);
                
                int half_search_region_size_x = ceil(2*sqrt(S->at<double>(0,0)));
                int half_search_region_size_y = ceil(2*sqrt(S->at<double>(1,1)));
                
                //cv::Mat *patches_for_correlation = new cv::Mat(pixels_in_the_matching_patch,(2*half_search_region_size_x+1)*(2*half_search_region_size_y+1)+1,CV_64FC1,0.0);
                
                /*
                for (j=0; j<pixels_in_the_matching_patch; j++)
                {
                    patches_for_correlation->at<double>(j,0) =
                }
                */
                //number of rows and columns in each patch (pre-defining)
                //Each patch is (half_patch_size_when_matching)x(half_patch_size_when_matching)
                int searchXSize = (2*half_search_region_size_x+1);
                int searchYSize = (2*half_search_region_size_y+1);
                
                //Using a vector here so we can treat each matrix individually
                //vector of candidate patches to test against predicted patch
                std::vector<cv::Mat*> patches_for_correlation;
                //vector of x/y coordinates corresponding to elements of patches_for_correlation
                std::vector<int> patchRowCoordinates;
                std::vector<int> patchColCoordinates;
                
                //int *match_candidates = (int*)malloc(2*(searchYSize*searchXSize+1)*sizeof(int));
                //int index_patches_for_correlation = 0;
                
                //iterate through each possible image in the search region of im
                //More precisely, iterate through each coordinate (j,k) of subimage in im
                for (j = round(current->h[0])-half_search_region_size_x; j < round(current->h[0])+half_search_region_size_x; j++){
                    for (k = round(current->h[1])-half_search_region_size_y; k < round(current->h[1])+half_search_region_size_y; k++){
                        cv::Mat *nu = new cv::Mat(2,1,CV_64FC1,0.0);
                        nu->at<double>(0,0) = j-current->h[0];
                        nu->at<double>(1,0) = k-current->h[1];
                        
                        cv::Mat *nu_transpose = new cv::Mat(1,2,CV_64FC1,0.0);
                        cv::transpose(*nu, *nu_transpose);
                        cv::Mat *A = new cv::Mat(1,1,CV_64FC1,0.0);
                        *A = (*nu_transpose)*(*invS)*(*nu);
                        nu->release();
                        nu_transpose->release();
                        if (A->at<double>(0,0) < chi_095_2){
                            //NSLog(@"img at coordinate %d,%d satisfies chi2 constraint", j,k);
                            if ((j>half_patch_size_when_matching) && (j<(cam->nCols-half_patch_size_when_matching-1)) &&
                                (k>half_patch_size_when_matching) && (k<(cam->nRows-half_patch_size_when_matching-1))){
                                int rowsInPatch = (2*half_patch_size_when_matching+1);
                                int colsInPatch = rowsInPatch; //Or, alternatively: (2*half_patch_size_when_matching+1)
                                uchar *image_patch = (uchar*)malloc(rowsInPatch*colsInPatch*sizeof(uchar));
                                uchar* imPtr;
                                uchar* imPatchPtr = image_patch;
                                //copy subimage to image_patch, push that patch into vector
                                for (int r = k-half_patch_size_when_matching; r<=k+half_patch_size_when_matching; ++r){
                                    imPtr = im->ptr<uchar>(r);
                                    for (int c = j-half_patch_size_when_matching; c<=j+half_patch_size_when_matching; ++c){
                                        *imPatchPtr = imPtr[c];
                                        ++imPatchPtr;
                                    }
                                }
                                patches_for_correlation.push_back(new cv::Mat(rowsInPatch, colsInPatch, CV_8UC1, image_patch, cv::Mat::AUTO_STEP));
                                free(image_patch);
                                //match_candidates[index_patches_for_correlation] = k;
                                //match_candidates[index_patches_for_correlation+(patchYSize*patchXSize+1)] = i;
                                patchRowCoordinates.push_back(k);
                                patchColCoordinates.push_back(j);
                                
                               // match_candidates->at<int>(0,index_patches_for_correlation-2) = j;
                               // match_candidates->at<int>(1,index_patches_for_correlation-2) = k;
                            }
                        }
                        //else{
                        //    NSLog(@"img at coordinate %d,%d does not satisfy chi2 constraint", j,k);
                        //}
                        A->release();
                    }
                }
                //Correlation coefficients:
                float bestCorrelation = -2.f;
                float initCorrelation = -2.f;
                int bestIdx = 0;
                cv::Mat templ = (*current->patch_when_matching).clone();
                templ.convertTo(templ, CV_8UC1);
                {
                    cv::Mat result;
                    cv::matchTemplate(templ, templ, result, CV_TM_CCOEFF);
                    float* correlationPtr = result.ptr<float>(0);
                    initCorrelation = *correlationPtr;
                    NSLog(@"Unity correlation %f", *correlationPtr);
                }
                for (j=0; j<patches_for_correlation.size(); ++j){
                    cv::Mat result;
                    cv::Mat compareImg = *patches_for_correlation.at(j);
                    cv::matchTemplate(compareImg, templ, result, CV_TM_CCOEFF);
                    float* correlationPtr = result.ptr<float>(0);
                    *correlationPtr = *correlationPtr/initCorrelation;
                    NSLog(@"result %d correlation %f", j, *correlationPtr);
                    if (*correlationPtr>bestCorrelation){
                        bestCorrelation = *correlationPtr;
                        bestIdx = j;
                    }
                    result.release();
                    patches_for_correlation.at(j)->release();
                }
                //Normalize
                NSLog(@"DEBUG: Best match: index %d correlation %f coordinates %d %d", bestIdx, bestCorrelation, patchRowCoordinates.at(bestIdx), patchColCoordinates.at(bestIdx));
                if (bestCorrelation>1.0){
                    NSLog(@"ERROR: Found a correlation better than unity");
                    return false;
                }
                
                //update z for this feature
                if (bestCorrelation>correlation_threshold){
                    current->individually_compatible = 1;
                    if (current->z==NULL){
                        current->z = (int*)malloc(2*sizeof(int));
                    }
                    current->z[0] = patchColCoordinates.at(bestIdx);
                    current->z[1] = patchRowCoordinates.at(bestIdx);
                }
                
                templ.release();
                invS->release();
            }
            else{
                NSLog(@"Eig for index %d not within limits...Skipping", i);
            }
            S->release();
            eigS->release();
        }
        else{
            NSLog(@"skipping NULL h at index %d", i);
        }
    }
    
    return true;
}