
//
//  RANSAC.m
//  RANSAC
//
//  Created by guest on 3/30/15.
//  Copyright (c) 2015 dryslam. All rights reserved.
//

#import "RANSAC.hpp"
#import <opencv2/opencv.hpp>
#import <opencv2/highgui/highgui_c.h>
#import <opencv2/highgui/cap_ios.h>
#import <math.h>
#include <stdio.h>
#include <stdlib.h>







//Helper functions:



/***************************************************************************************************************************************/
/********************************************************FUNCTIONAL METHODS*************************************************************/
/***************************************************************************************************************************************/









/***************************************************************************************************************************************/
/**********************************************METHODS THAT MUST BE IMPLEMENTED*********************************************************/
/***************************************************************************************************************************************/




/***************************************************************************************************************************************/
/********************************************************AUXILIARY METHODS**************************************************************/
/***************************************************************************************************************************************/


















//bool rescue_hi_inliers(Filter* filter, NSMutableArray* features_info, Camera* cam){
//    return true;
//}


////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////
/////I.mapManagement
//I.1)deleteFeatures


/*
 //I.2)updateFeatureInfo
 bool updateFeaturesInfo(NSMutableArray* features_info)
 {
 int i;
 
 for (i=0;i<(int)[features_info count];i++)
 {
 Feature *current=features_info[i];
 if (sizeof(current->h))     //~isempty(features_info(i).h)
 current->times_predicted=current->times_measured+1;
 if (current->low_innovation_inlier||current->high_innovation_inlier)
 current->times_measured++;
 
 current->individually_compatible=0;
 current->low_innovation_inlier=0;
 current->high_innovation_inlier=0;
 
 current->h[0]=0;
 current->h[1]=0;
 
 //current->z=NULL;
 //current->H=cv::Mat;
 
 current->S[0][0]=0;
 current->S[0][1]=0;
 current->S[1][0]=0;
 current->S[1][1]=0;
 }
 
 return true;
 }
*/
//I.3)inverseDepth2Cartesian(corresponding to inversedepth_2_cartesian)
/*
void m(double *a,double *m)
{
    
}
*/
/*
void m(double a,double b,cv::Mat m)
{
    m.at<double>(0,0)=cos(b)*sin(a);
    m.at<double>(1,0)=-sin(b);
    m.at<double>(2,0)=cos(b)*cos(a);
}
/*
void inverse2cartesian(cv::Mat inverseDepth,cv::Mat cartesian)
{
    cv::Mat rw=cv::Mat(3,inverseDepth.cols,CV_64FC1);
    cv::Mat theta=cv::Mat(1,inverseDepth.cols,CV_64FC1);
    cv::Mat phi=cv::Mat(1,inverseDepth.cols,CV_64FC1);
    cv::Mat rho=cv::Mat(1,inverseDepth.cols,CV_64FC1);
    
    cv::Mat m=cv::Mat(3,inverseDepth.cols,CV_64FC1);
    
    inverseDepth.rowRange(0,2).copyTo(rw);
    inverseDepth.row(3).copyTo(theta);
    inverseDepth.row(4).copyTo(phi);
    inverseDepth.row(5).copyTo(rho);
    
    int i;
    
    for (i=0;i<inverseDepth.cols;i++)
    {
        m.at<double>(0,i)=cosf(phi.at<double>(0,i))*sinf(theta.at<double>(0,i));
        m.at<double>(1,i)=-sinf(phi.at<double>(0,i));
        m.at<double>(2,i)=cosf(phi.at<double>(0,i))*cosf(theta.at<double>(0,i));
        
        cartesian.at<double>(0,i)=rw.at<double>(0,i)+(1/rho.at<double>(0,i)*m.at<double>(0,i));
        cartesian.at<double>(1,i)=rw.at<double>(1,i)+(1/rho.at<double>(0,i)*m.at<double>(1,i));
        cartesian.at<double>(2,i)=rw.at<double>(2,i)+(1/rho.at<double>(0,i)*m.at<double>(2,i));
    }
}
*/
/*
void inverseDepth_2_Cartesian(NSMutableArray* features_info,Filter* filter)
{
    double linearityIndexThreshold=0.1;
    
    cv::Mat X=cv::Mat(filter->state_size,1,CV_64FC1,&filter->x_k_k); //Please check later.
    cv::Mat P=cv::Mat(filter->state_size,filter->state_size,CV_64FC1,&filter->p_k_k);
    
    double std_rho,rho,std_d;
    double theta,phi;
    
    cv::Mat mi=cv::Mat(3,1,CV_64FC1);
    cv::Mat x_c1=cv::Mat(3,1,CV_64FC1);
    cv::Mat x_c2=cv::Mat(3,1,CV_64FC1);
    cv::Mat p=cv::Mat(3,X.cols,CV_64FC1);
    cv::Mat p_0=cv::Mat(1,3,CV_64FC1);
    
    double d_c2p,cos_alpha;
    
    int i,j;
    
    for (i=0;i<(int)[features_info count];i++)
    {
        Feature *current=features_info[i];
        
        if ([current->type isEqualToString:@"inversedepth"])
        {
            int initialPositionOfFeature=14;
            
            for (j=0;j<i-2;j++)
            {
                if ([current->type isEqualToString:@"cartesian"])
                    initialPositionOfFeature+=3;
                if ([current->type isEqualToString:@"inversedepth"])
                    initialPositionOfFeature+=6;
            }
            
            std_rho=sqrtf(P.at<double>(initialPositionOfFeature+4,initialPositionOfFeature+4));
            rho=X.at<double>(initialPositionOfFeature+4,0);
            std_d=std_rho/powf(rho,2);
            
            theta=X.at<double>(initialPositionOfFeature+2,0);
            phi=X.at<double>(initialPositionOfFeature+3,0);
            m(theta,phi,mi);
            
            X.rowRange(initialPositionOfFeature-1,initialPositionOfFeature+1).copyTo(x_c1);
            X.rowRange(0,3).copyTo(x_c2);
            
            inverse2cartesian(X.rowRange(initialPositionOfFeature-1,initialPositionOfFeature+4),p);
            
            d_c2p=cv::norm(p-x_c2);
            
            //Build error: variables aren't defined
            if (linnearityIndex<linearityIndexThreshold)
             {
             
             current->type=@"cartesian";
             }
        }
    }
    
}
*/
/*
 bool inversedepth_2_cartesian(NSMutableArray* features_info_size, double Cartesian,
 double *x_k_k, double* p_k_k, int state_size)
 {
 int J=sizeof(*(inverseDepth+0))/sizeof(*((inverseDepth+0)+0)); //
 
 double rw[3][J];
 double theta[J],stheta[J],ctheta[J];
 double phi[J],sphi[J],cphi[J];
 double rho[J],srho[J],crho[J];
 
 double m[3][J];
 
 int i,j;
 
 for (j=0;j<J;j++)
 {
 theta[j]=*((inverseDepth+3)+j);
 stheta[j]=sin(theta[j]);
 ctheta[j]=cos(theta[j]);
 
 phi[j]=*((inverseDepth+4)+j);
 sphi[j]=sin(phi[j]);
 cphi[j]=cos(phi[j]);
 
 rho[j]=*((inverseDepth+5)+j);
 srho[j]=sin(rho[j]);
 crho[j]=cos(rho[j]);
 
 for (i=0;i<3;i++)
 {
 rw[i][j]=*((inverseDepth+i)+j);
 }
 
 m[0][j]=cphi[j]*stheta[j];
 m[1][j]=-sphi[j];
 m[2][j]=cphi[j]*ctheta[j];
 
 *((Cartesian+i)+j)=rw[i][j]+1/rho[j]*m[0][j];
 }
 
 return true;
 }
 */

//I.4)initializeFeatures
/*structs replaced with classes
 void initializeFeatures(NSMutableArray* featureInfo,Filter* filter,int step,Camera* cam,double *im,int numberOfFeaturesToInitialized)
 {
 
 }
 */

/*see header for correct definition
 void initializeOneFeature(NSMutableArray* featureInfo,Filter* filter,int step,Camera* cam,double *im_k)
 {
 
 }*/
/*
void q2r(double *q, cv::Mat* R)
//completed.
//Plaese initialize R as a 3*3 array.
{
    double x=*(q+1);
    double y=*(q+2);
    double z=*(q+3);
    double r=*q;
    
    R->at<double>(0,0)=r*r+x*x-y*y-z*z;
    R->at<double>(1,2)=2*(x*y-r*z);
    R->at<double>(1,3)=2*(z*x+r*y);
    R->at<double>(2,1)=2*(x*y+r*z);
    R->at<double>(2,2)=r*r-x*x+y*y-z*z;
    R->at<double>(2,3)=2*(y*z-r*x);
    R->at<double>(3,1)=2*(z*x-r*y);
    R->at<double>(3,2)=2*(y*z+r*x);
    R->at<double>(3,3)=r*r-x*x-y*y+z*z;
}
*/
//void hu(cv::Mat yi,Camera *cam,cv::Mat uv_u)
/*
 bool hu(double *uv_u,double *yi,int yi_ncols,Camera *cam)
 {
 double u0=cam->Cx;
 double v0=cam->Cy;
 double f=cam->f;
 double ku=1/cam->dx;
 double kv=1/cam->dy;
 
 int i;
 
 for (i=0;i<yi_ncols;i++)
 {
 uv_u[i+0*yi_ncols]=u0+f*ku*yi[i+0*yi_ncols]/yi[i+2*yi_ncols];
 uv_u[i+1*yi_ncols]=v0+f*kv*yi[i+1*yi_ncols]/yi[i+2*yi_ncols];
 }
 
 return true;
 }
 */
/*
 bool distortFM(double *uvd,double *uv,int uv_ncols,Camera* cam)
 {
 double k1=cam->k1;
 double k2=cam->k2;
 double Cx=cam->Cx;
 double Cy=cam->Cy;
 double dx=cam->dx;
 double dy=cam->dy;
 */
 //cv::Mat xu=cv::Mat(1,uv_mat.cols,CV_64FC1);
 //cv::Mat yu=cv::Mat(1,uv_mat.cols,CV_64FC1);
 //cv::Mat ru=cv::Mat(1,uv_mat.cols,CV_64FC1);
 //cv::Mat rd=cv::Mat(1,uv_mat.cols,CV_64FC1);
 //cv::Mat D=cv::Mat(1,uv_mat.cols,CV_64FC1);
 //cv::Mat xd=cv::Mat(1,uv_mat.cols,CV_64FC1);
 //cv::Mat yd=cv::Mat(1,uv_mat.cols,CV_64FC1);
 //double xu,yu,ru,rd,D,xd,yd;
 /*
 double f;
 double f_p;
 
 int i,j;
 
 for (i=0;i<uv_ncols;i++)
 {
 xu=(uv[i+0*uv_ncols]-Cx)*dx;
 yu=(uv[i+1*uv_ncols]-Cy)*dy;
 
 ru=sqrt(pow(xu,2)+pow(yu,2));
 rd=ru/(1+k1*pow(ru,2)+k2*pow(ru,4));
 
 for (j=0;j<10;j++)
 {
 f=rd+k1*pow(rd,3)+k2*pow(rd,5)-ru;
 f_p=1+3*k1*pow(rd,2)+5*k2*pow(rd,4);
 rd=rd-f/f_p;
 }
 
 D=1+k1*pow(rd,2)+k2*pow(rd,4);
 xd=xu/D;
 yd=yu/D;
 
 uvd[i+0*uv_ncols]=xd/dx+Cx;
 uvd[i+1*uv_ncols]=yd/dy+Cy;
 }
 
 return true;
 }
*/
/*
void distortOnePoint(Camera* cam,double *uvu,double *uvd)
//completed. Please initialize uvd as a 1*2 array.
{
    double k1=cam->k1;
    double k2=cam->k2;
    double Cx=cam->Cx;
    double Cy=cam->Cy;
    double dx=cam->dx;
    double dy=cam->dy;
    
    double xu=(*uvu-Cx)*dx;
    double yu=(*(uvu+1)-Cy)*dy;
    
    double ru=sqrt(xu*xu+yu*yu);
    double rd=ru/(1+k1*pow(ru,2)+k2*pow(ru,4));
    
    int numberOfIterations=20;
    int i;
    double f,f_p;
    
    for (i=0;i<numberOfIterations;i++)
    {
        f=rd+k1*pow(rd,3)+k2*pow(rd,5)-ru;
        f_p=1+3*k1*pow(rd,2)+5*k2*pow(rd,4);
        rd-=f/f_p;
    }
    
    double D=1+k1*pow(rd,2)+k2*pow(rd,4);
    double xd=xu/D;
    double yd=yu/D;
    
    *uvd=xd/dx+Cx;
    *((uvd+0)+1)=yd/dy+Cy;
    
}
*/
/*
int hiCartesian(NSMutableArray* features_info,Camera* cam,cv::Mat yi,cv::Mat t_wc,cv::Mat r_wc,cv::Mat zi)
{
    const double pi=3.1415926;
    
    cv::Mat r_cw=cv::Mat(r_wc.rows,r_wc.cols,CV_64FC1);
    cv::Mat hrl=cv::Mat(r_wc.rows,yi.cols,CV_64FC1);
    cv::Mat uv_u=cv::Mat(2,yi.cols,CV_64FC1);
    cv::Mat uv_d=cv::Mat(2,yi.cols,CV_64FC1);
    
    bool a0,a1,a2,a3;
    bool b[hrl.cols];
    
    int i;
    
    cv::invert(r_wc,r_cw);
    hrl=r_cw*(yi-t_wc);
    
    for (i=0;i<hrl.cols;i++)
    {
        a0=atan2f(hrl.at<double>(0,i), hrl.at<double>(2,i))*180/pi<-60;
        a1=atan2f(hrl.at<double>(0,i), hrl.at<double>(2,i))*180/pi>+60;
        a2=atan2f(hrl.at<double>(0,i), hrl.at<double>(2,i))*180/pi<-60;
        a3=atan2f(hrl.at<double>(0,i), hrl.at<double>(2,i))*180/pi>+60;
        
        b[i]=a0||a1||a2||a3;
        
        if (b[i]==0)
            return 0;
    }
    
    //    hu(hrl,cam,uv_u);
    
    // distortFM(uv_u,cam,uv_d);
    
    if ((uv_d.at<double>(0)>0)&&(uv_d.at<double>(0)<cam->nCols)&&  //Please check later.
        (uv_d.at<double>(1)>0)&&(uv_d.at<double>(1)<cam->nRows))
    {
        uv_d.copyTo(zi);
        return 0;
    }
    
    return 0;
}
*/
/*
int hiInverseDepth(NSMutableArray* features_info,Camera* cam,cv::Mat yinit,cv::Mat t_wc,cv::Mat r_wc,cv::Mat zi)
{
    const double pi=3.1415926;
    
    cv::Mat r_cw=cv::Mat(r_wc.cols,r_wc.rows,CV_64FC1);
    cv::Mat yi=cv::Mat(1,3,CV_64FC1);
    cv::Mat mi=cv::Mat(3,1,CV_64FC1);
    cv::Mat hrl=cv::Mat(r_wc.cols,yi.cols,CV_64FC1);
    cv::Mat uv_u=cv::Mat(2,yi.cols,CV_64FC1);
    cv::Mat uv_d=cv::Mat(2,yi.cols,CV_64FC1);
    
    double theta,phi,rho;
    
    bool a0,a1,a2,a3;
    bool b[hrl.cols];
    
    int i;
    
    cv::transpose(r_wc,r_cw);
    
    for (i=0;i<3;i++)
    {
        yi.at<double>(0,i)=yinit.at<double>(0,i);
    }
    
    theta=yinit.at<double>(0,3);
    phi=yinit.at<double>(0,4);
    rho=yinit.at<double>(0,5);
    
    m(theta,phi,mi);
    
    hrl=r_cw*((yi-t_wc)*rho+mi);
    
    for (i=0;i<hrl.cols;i++)
    {
        a0=atan2f(hrl.at<double>(0,i), hrl.at<double>(2,i))*180/pi<-60;
        a1=atan2f(hrl.at<double>(0,i), hrl.at<double>(2,i))*180/pi>+60;
        a2=atan2f(hrl.at<double>(0,i), hrl.at<double>(2,i))*180/pi<-60;
        a3=atan2f(hrl.at<double>(0,i), hrl.at<double>(2,i))*180/pi>+60;
        
        b[i]=a0||a1||a2||a3;
        
        if (b[i]==0)
            return 0;
    }
    
    //hu(hrl,cam,uv_u);
    
    //distortFM(uv_u,cam,uv_d);
    
    if ((uv_d.at<double>(0)>0)&&(uv_d.at<double>(0)<cam->nCols)&&  //Please check later.
        (uv_d.at<double>(1)>0)&&(uv_d.at<double>(1)<cam->nRows))
    {
        uv_d.copyTo(zi);
        return 0;
    }
    
    return 0;
}
*/
/*
bool predict_camera_measurements(NSMutableArray *features_info,Filter *filter,Camera *cam)
{
    cv::Mat t_wc=cv::Mat(1,3,CV_64FC1);
    cv::Mat x_k_k_0=cv::Mat(1,4,CV_64FC1);
    cv::Mat r_wc=cv::Mat(3,3,CV_64FC1);
    cv::Mat yi=cv::Mat(1,3,CV_64FC1);
    cv::Mat hi=cv::Mat(2,1,CV_64FC1);  //Please check later.
    
    int i,j;
    int J=0;
    
    for (j=0;filter->x_k_k+j;j++)
        J++;
    
    cv::Mat features0=cv::Mat(1,J-13,CV_64FC1);
    cv::Mat features=cv::Mat(1,features0.cols-3,CV_64FC1);
    
    for (i=0;i<(int)[features_info count];i++)
    {
        Feature *current=features_info[i];
        
        if ([current->type isEqualToString:@"cartesian"])
        {
            for (j=0;j<3;j++)
                yi.at<double>(0,i)=features0.at<double>(0,i);
            for (j=3;j<features0.cols;j++)
                features0.at<double>(0,i-3)=features0.at<double>(0,i);
            
            hiCartesian(features_info,cam,yi,t_wc,r_wc,hi);
            
            if ((current->h[0])&&(current->h[1]))
            {
                current->h[0]=hi.at<double>(0,0);
                current->h[1]=hi.at<double>(1,0);
            }
        }
        
        if ([current->type isEqualToString:@"inversedepth"])
        {
            for (j=0;j<6;j++)
                yi.at<double>(0,i)=features0.at<double>(0,i);
            for (j=6;j<features0.cols;j++)
                features0.at<double>(0,i-6)=features0.at<double>(0,i);
            
            hiCartesian(features_info,cam,yi,t_wc,r_wc,hi);
            
            if ((current->h[0])&&(current->h[1]))
            {
                current->h[0]=hi.at<double>(0,0);
                current->h[1]=hi.at<double>(1,0);
            }
        }
    }
    return true;
}
*/
/*
 structs replaced with classes
 void fastCornerDetect9(double *im, double threshold,double *coords)
 {
 
 }
 
 void fastNonMax(double *im,double barrier,double c,double ret)
 {
 
 }
 
 void addFeaturesInveerseDepth(double *uvd,double *X,double *P,Camera* cam,double *std_px1,double *initialRho,double *stdRho,double *XRES,double *PRES,double *newFeature)
 {
 int nNewFeat=sizeof(*(uvd+0))/sizeof(*((uvd+0)+0)); //nNewFeat=size(uvd,2)
 
 if (nNewFeat==0)
 {
 
 }
 }
 */

/*
void undistortFM(cv::Mat uvd,Camera *cam,cv::Mat uvu)
{
    double Cx=cam->Cx;
    double Cy=cam->Cy;
    double k1=cam->k1;
    double k2=cam->k2;
    double dx=cam->dx;
    double dy=cam->dy;
    
    int j;
    
    cv::Mat xu=cv::Mat(1,uvd.cols,CV_64FC1);
    cv::Mat yu=cv::Mat(1,uvd.cols,CV_64FC1);
    cv::Mat ru=cv::Mat(1,uvd.cols,CV_64FC1);
    cv::Mat rd=cv::Mat(1,uvd.cols,CV_64FC1);
    cv::Mat D=cv::Mat(1,uvd.cols,CV_64FC1);
    cv::Mat xd=cv::Mat(1,uvd.cols,CV_64FC1);
    cv::Mat yd=cv::Mat(1,uvd.cols,CV_64FC1);
    
    for (j=0;j<uvd.cols;j++)
    {
        xd.at<double>(0,j)=(uvd.at<double>(0,j)-Cx)*dx;
        yd.at<double>(0,j)=(uvd.at<double>(0,j)-Cy)*dy;
        rd.at<double>(0,j)=sqrtf(powf(xd.at<double>(0,j),2)+powf(yd.at<double>(0,j),2));
        D.at<double>(0,j)=1+k1*powf(rd.at<double>(0,j),2)+k2*powf(rd.at<double>(0,j),4);
        xu.at<double>(0,j)=xd.at<double>(0,j)*D.at<double>(0,j);
        yu.at<double>(0,j)=yd.at<double>(0,j)*D.at<double>(0,j);
        
        uvu.at<double>(0,j)=xu.at<double>(0,j)/dx+Cx;
        uvu.at<double>(1,j)=yu.at<double>(0,j)/dy+Cy;
    }
}
*/
/*
void undistortOnePoint(double *uvd,Camera* cam,double *uvu)
//complete. Please initialize uvu as a 1*2 array.
{
    double Cx=cam->Cx;
    double Cy=cam->Cy;
    double k1=cam->k1;
    double k2=cam->k2;
    double dx=cam->dx;
    double dy=cam->dy;
    
    double ud=*uvd;
    double vd=*(uvd+1);
    double rd=sqrt(pow(dx*(ud-Cx),2)+pow(dy*(vd-Cy),2));
    
    *((uvu+0)+0)=Cx+(ud-Cx)*(1+k1*pow(rd,2)+k2*pow(rd,4));
    *((uvu+0)+1)=Cy+(vd-Cy)*(1+k1*pow(rd,2)+k2*pow(rd,4));
    
}
*/
/*
bool compute_hypothesis_support(cv::Mat xi, int state_size, Camera* cam, cv::Mat state_vector_pattern,
                                cv::Mat z_id, cv::Mat z_euc, int threshold, int* hypothesis_support,
                                NSMutableArray* positions_li_inliers_id, NSMutableArray* positions_li_inliers_euc)
{
    int hypothesisSupport=0;
    int n_id=0;
    
    double u0=cam->Cx;
    double v0=cam->Cy;
    double f=cam->f;
    double ku=1/cam->dx;
    double kv=1/cam->dy;
    
    int i;
    int sum=0;
    
    if (z_id.cols) //if ~isempty(z_id)
    {
        n_id=z_id.cols;
        
        cv::Mat ri=cv::Mat(3*n_id,1,CV_64FC1);
        cv::Mat anglesi=cv::Mat(2*n_id,1,CV_64FC1);
        cv::Mat rhoi=cv::Mat();
        
        
    }
    
    hypothesis_support+=sum;
    
    return true;
}*/
/*
void hinv(cv::Mat uvd,cv::Mat Xv,Camera* cam,double initialRho,cv::Mat newFeature)
{
    double fku=cam->k->at<double>(0,0);
    double fkv=cam->k->at<double>(1,1);
    double U0=cam->k->at<double>(0,2);
    double V0=cam->k->at<double>(1,2);
    
    double r_W[3];
    double q_WR[4];
    
    cv::Mat h_LR=cv::Mat(3,1,CV_64FC1, 0.0);
    cv::Mat* q_WR_0= new cv::Mat;
    *q_WR_0 = cv::Mat(3,3,CV_64FC1, 0.0);
    cv::Mat uv=cv::Mat(1,2,CV_64FC1, 0.0); //Please check later.
    cv::Mat n=cv::Mat(3,1,CV_64FC1, 0.0);
    
    double u,v;
    
    int i;
    
    undistortFM(uvd,cam,uv);
    u=uv.at<double>(0,0); //Please check later.
    v=uv.at<double>(0,1); //Please check later.
    
    for (i=0;i<3;i++)
        r_W[i]=Xv.at<double>(0,i);
    for (i=3;i<7;i++)
        q_WR[i-3]=Xv.at<double>(0,i);
    
    h_LR.at<double>(0,0)=-(U0-u)/fku;
    h_LR.at<double>(1,0)=-(V0-v)/fkv;
    h_LR.at<double>(2,0)=1;
    
    q2r(q_WR,q_WR_0);
    
    n=*q_WR_0*h_LR;
    
    for (i=0;i<3;i++)
    {
        newFeature.at<double>(0,i)=*(r_W+i);
        newFeature.at<double>(1,i)=atan2f(n.at<double>(0,0), n.at<double>(2,0));
        newFeature.at<double>(2,i)=atan2f(-n.at<double>(1,0),sqrtf(powf(n.at<double>(0,0),2)+powf(n.at<double>(2,0),2)));
        newFeature.at<double>(3,i)=initialRho;
    }
}
*/
/*
void dRBydq0(double *q,cv::Mat dRBydq0RES)
//completed.
//Plaese initialize dRBydq0RES as a3*3 array.
{
    double q0=*q;
    double qx=*(q+1);
    double qy=*(q+2);
    double qz=*(q+3);
    
    dRBydq0RES.at<double>(0,0)=2*q0;
    dRBydq0RES.at<double>(0,1)=-2*qz;
    dRBydq0RES.at<double>(0,2)=2*qy;
    dRBydq0RES.at<double>(1,0)=2*qz;
    dRBydq0RES.at<double>(1,1)=2*q0;
    dRBydq0RES.at<double>(1,2)=-2*qx;
    dRBydq0RES.at<double>(2,0)=-2*qy;
    dRBydq0RES.at<double>(2,1)=2*qx;
    dRBydq0RES.at<double>(2,2)=-2*qx;
}
*/
/*
void dRBydqx(double *q,cv::Mat dRBydqxRES)
//completed.
//Plaese initialize dRBydqxRES as a3*3 array.
{
    double q0=*q;
    double qx=*(q+1);
    double qy=*(q+2);
    double qz=*(q+3);
    
    dRBydqxRES.at<double>(0,0)=2*qx;
    dRBydqxRES.at<double>(0,1)=2*qy;
    dRBydqxRES.at<double>(0,2)=2*qz;
    dRBydqxRES.at<double>(1,0)=2*qy;
    dRBydqxRES.at<double>(1,1)=-2*qx;
    dRBydqxRES.at<double>(1,2)=-2*q0;
    dRBydqxRES.at<double>(2,0)=2*qz;
    dRBydqxRES.at<double>(2,1)=2*q0;
    dRBydqxRES.at<double>(2,2)=-2*qx;
}
*/
/*
void dRBydqy(double *q,cv::Mat dRBydqyRES)
//completed.
//Plaese initialize dRBydy0RES as a3*3 array.
{
    double q0=*q;
    double qx=*(q+1);
    double qy=*(q+2);
    double qz=*(q+3);
    
    dRBydqyRES.at<double>(0,0)=-2*qy;
    dRBydqyRES.at<double>(0,1)=2*qx;
    dRBydqyRES.at<double>(0,2)=2*q0;
    dRBydqyRES.at<double>(1,0)=2*qx;
    dRBydqyRES.at<double>(1,1)=2*qy;
    dRBydqyRES.at<double>(1,2)=2*qz;
    dRBydqyRES.at<double>(2,0)=-2*q0;
    dRBydqyRES.at<double>(2,1)=2*qz;
    dRBydqyRES.at<double>(2,2)=-2*qy;
}
*/
/*
void dRBydqz(double *q,cv::Mat dRBydqzRES)
//completed.
//Plaese initialize dRBydy0RES as a3*3 array.
{
    double q0=*q;
    double qx=*(q+1);
    double qy=*(q+2);
    double qz=*(q+3);
    
    dRBydqzRES.at<double>(0,0)=-2*qz;
    dRBydqzRES.at<double>(0,1)=-2*q0;
    dRBydqzRES.at<double>(0,2)=2*qx;
    dRBydqzRES.at<double>(1,0)=2*q0;
    dRBydqzRES.at<double>(1,1)=-2*qz;
    dRBydqzRES.at<double>(1,2)=2*qy;
    dRBydqzRES.at<double>(2,0)=2*qx;
    dRBydqzRES.at<double>(2,1)=2*qy;
    dRBydqzRES.at<double>(2,2)=2*qz;
}
*/
/*

*/
/*
void JacobUndistorFM(Camera* cam,cv::Mat* uvd,cv::Mat* J_unidistor)
//complete. Please initialize initial J_unidistor as a 2*2 array.
{
    double k1=cam->k1;
    double k2=cam->k2;
    double Cx=cam->Cx;
    double Cy=cam->Cy;
    double dx=cam->dx;
    double dy=cam->dy;
    
    double ud=uvd->at<double>(0,0);
    double vd=uvd->at<double>(0,1);
    
    double xd=(ud-Cx)*dx;
    double yd=(ud-Cy)*dy;
    
    double rd2=xd*xd+yd*yd;
    double rd4=rd2*rd2;
    
    double uu_ud=(1+k1*rd2+k2*rd4)+(ud-Cx)*(k1+2*k2*rd2)*(2*(ud-Cx)*dx*dx);
    double vu_vd=(1+k1*rd2+k2*rd4)+(vd-Cy)*(k1+2*k2*rd2)*(2*(vd-Cy)*dy*dy);
    double uu_vd=(ud-Cx)*(k1+2*k2*rd2)*(2*(vd-Cy)*dy*dy);
    double vu_ud=(vd-Cy)*(k1+2*k2*rd2)*(2*(ud-Cx)*dx*dx);
    
    J_unidistor->at<double>(0,0)=uu_ud;
    J_unidistor->at<double>(0,1)=uu_vd;
    J_unidistor->at<double>(1,0)=vu_ud;
    J_unidistor->at<double>(1,1)=vu_vd;
}
*/
/*

*/
/*
void addFeatureToInfoVector(NSMutableArray* features_info,double *uv,double *im_k,double *X_RES,int step,double *newFeature)
{
    int new_position=(int)[features_info count]+1;
    double X_RES0[4];
    
    int a=*(uv+0);
    int b=*(uv+1);
    
    Feature *next;
    //This is initialized when you create a feature
    //Make sure to call when creating = [[Feature alloc]init];
    //next->R_wc_when_initialized=cv::Mat::zeros(3,3,CV_64FC1);
    
    int i,j;
    
    for (i=0;i<13;i++)
    {
        for (j=0;j<13;j++)
        {
            //next->patch_when_initialized.at<double>(i,j)=
            //*(im_k+b-20)=0;
            //*(im_k+a-20)=0;
        }
    }
    
    //This is done upon allocation
    //next->patch_when_matching=cv::Mat::zeros(41,41,CV_64FC1);
    
    for (i=0;i<3;i++)
        *(next->r_wc_when_initialized+i)=X_RES[i];
    
    for (i=3;i<7;i++)
        X_RES0[i-3]=*(X_RES+i);
    q2r((double *)X_RES0,next->R_wc_when_initialized);
    
    next->uv_when_initialized[0]=*(uv+0); //Please check later.
    next->uv_when_initialized[1]=*(uv+1); //Please check later.
    
    next->half_patch_size_when_initialized=20;
    next->half_patch_size_when_matching=6;
    next->times_predicted=0;
    next->times_measured=0;
    next->init_frame=step;
    
    next->init_measurement[0][0]=*(uv+0); //Please check later.
    next->init_measurement[1][0]=*(uv+1); //Please check later.
    
    next->type=@"inversedepth";
    
    for (i=0;i<6;i++)
        next->yi[i]=*(newFeature+i);
    
    next->individually_compatible=0;
    next->low_innovation_inlier=0;
    next->high_innovation_inlier=0;
    
    next->z[0] = 0;
    next->z[1] = 0;
    
    next->h[0]=0;
    next->h[1]=0;
    
    //features_info(new_position).H=[];
    
    next->S[0][0]=0;
    next->S[0][1]=0;
    next->S[1][0]=0;
    next->S[1][1]=0;
    
    next->state_size=6;
    
    next->measurement_size=2;
    
    next->R[0][0]=1;
    next->R[0][1]=0;
    next->R[1][0]=0;
    next->R[1][1]=1;
    
    [features_info insertObject:(id)next
                        atIndex:(NSInteger)new_position];
}
*/
////////////////////////////////////////////////////////////////////////////////////////////////////
/////II.searchICMatches
//II.1) PredictCameraMeasurements
/*
 void predictCameraMeasurements(NSMutableArray* features_info, Filter* filter, Camera* cam,int stepnumber)
 {
 /*
 int end=sizeof(filter.x_k_k)/sizeof(double);
 int i;
 
 double t_wc[3]={0};
 double r_wc0[3]={0};
 double r_wc[3][3]={0};
 double features0[end-13]={0};
 
 double yi[3]={0};
 double features[end-3]={0};
 
 for (i=0;i<3;i++)
 {
 t_wc[i]=filter.x_k_k[i];
 }
 
 for (i=3;i<7;i++)
 {
 r_wc0[i-3]=filter.x_k_k[i];
 }
 q2r((double *)r_wc0,(double *)r_wc);
 
 for (i=13;i<end;i++)
 {
 features[i-13]=filter.x_k_k[i];
 }
 
 for (i=0;i<stepnumber;i++)
 {
 int j;
 
 if (featureInfo.type[i][0]=='c')
 {
 for (j=0;j<3;j++)
 {
 yi[j]=features0[j];
 }
 
 for (j=3;j<end;j++)
 {
 features[j-3]=features0[j];
 }
 
 hiCartesian(featureInfo,cam,(double *)yi,(double *)t_wc,(double *)r_wc,(double *)hi);
 
 if (hi) //if (~isempty(hi))
 featuresInfo.h[i]=
 }
 
 if (featureInfo.type[i][0]=='i')
 {
 
 }
 }
 
 }
 */
//II.2)calculateDerivatives
/*
void dhrl_dy(double *Xv_km1_k,cv::Mat* a)
{
    double *Xv_km1_k_0;
    cv::Mat* a0= new cv::Mat;
    *a0 = cv::Mat::zeros(3,3,CV_64FC1);
    
    int i;
    
    for (i=3;i<7;i++)
    {
        Xv_km1_k_0=&(Xv_km1_k[i]);
        Xv_km1_k_0++;
    }
    
    q2r(Xv_km1_k_0,a0);
    
    cv::invert(*a0, *a);
    a0->release();
}
*/
/*
void dhrl_drw(double *Xv_km1_k,cv::Mat* a)
{
    double *Xv_km1_k_0;
    cv::Mat* a0= new cv::Mat;
    *a0 = cv::Mat::zeros(3,3,CV_64FC1);
    
    int i;
    
    for (i=3;i<7;i++)
    {
        Xv_km1_k_0=&(Xv_km1_k[i]);
        Xv_km1_k_0++;
    }
    
    q2r(Xv_km1_k_0,a0);
    
    cv::invert(*a0, *a);
    *a=-1*(*a);
}
*/
/*
void dhd_dhu(Camera *cam,cv::Mat* zi_d,cv::Mat* a)
{
    cv::Mat* a0= new cv::Mat;
    *a0 = cv::Mat::zeros(2,2,CV_64FC1);
    JacobUndistorFM(cam,zi_d, a0);
    cv::invert(*a0,*a);
    
    a0->release();
}
*/
/*
void dhu_dhrl(Camera *cam,double *Xv_km1_k,cv::Mat* yi,cv::Mat* a)
{
    double f=cam->f;
    double ku=1/cam->dx;
    double kv=1/cam->dy;
    
    double *Xv_km1_k_0;
    
    cv::Mat rw=cv::Mat(3,1,CV_64FC1);
    cv::Mat* Rrw_0= new cv::Mat;
    *Rrw_0 = cv::Mat::zeros(1,4,CV_64FC1);
    cv::Mat Rrw=cv::Mat::zeros(3,3,CV_64FC1);
    cv::Mat hrl=cv::Mat(3,1,CV_64FC1);
    
    int i;
    
    for (i=0;i<3;i++)
        rw.at<double>(i,0)=*(Xv_km1_k+i);
    
    for (i=3;i<7;i++)
    {
        Xv_km1_k_0=&(Xv_km1_k[i]);
        Xv_km1_k_0++;
    }
    
    q2r(Xv_km1_k,Rrw_0);
    cv::invert(*Rrw_0,Rrw);
    
    hrl=Rrw*(*yi-rw);
    
    a->at<double>(0,0)=f*ku/hrl.at<double>(2,0);
    a->at<double>(0,1)=0;
    a->at<double>(0,2)=-hrl.at<double>(0,0)*f*ku/(powf(hrl.at<double>(2,0),2));
    a->at<double>(1,0)=0;
    a->at<double>(1,1)=f*kv/hrl.at<double>(2,0);
    a->at<double>(1,2)=-hrl.at<double>(1,0)*f*kv/(powf(hrl.at<double>(2,0),2));
    
    rw.release();
    Rrw_0->release();
    Rrw.release();
    hrl.release();
}
*/
/*
void dhrl_dqwr(double *Xv_km1_k,cv::Mat* yi,cv::Mat* a)
{
    cv::Mat* Xv_km1_k_0= new cv::Mat;
    *Xv_km1_k_0 = cv::Mat(1,3,CV_64FC1);
    
    double Xv_km1_k_1[4];
    
    int i;
    
    for (i=0;i<3;i++)
        Xv_km1_k_0->at<double>(0,i)=*(Xv_km1_k+i);
    
    Xv_km1_k_1[0]=Xv_km1_k[0];
    Xv_km1_k_1[1]=-Xv_km1_k[1];
    Xv_km1_k_1[2]=-Xv_km1_k[2];
    Xv_km1_k_1[3]=-Xv_km1_k[3];
    
    cv::Mat* arg = new cv::Mat;
    *arg = *yi-*Xv_km1_k_0;
    
    dRqTimesABydq((double *)Xv_km1_k_1,arg,a);
    
    //a*=dqbar_by_dq;
    arg->release();
    Xv_km1_k_0->release();
    
}
*/
/*
void dh_dhrl(Camera *cam,double *Xv_km1_k,cv::Mat* yi,cv::Mat* zi,cv::Mat* a)
{
    cv::Mat* a1= new cv::Mat;
    *a1 = cv::Mat::zeros(2,2,CV_64FC1);
    cv::Mat* a2= new cv::Mat;
    *a2 = cv::Mat::zeros(2,3,CV_64FC1);
    
    dhd_dhu(cam,zi,a1);
    dhu_dhrl(cam,Xv_km1_k,yi,a2);
    
    *a=*a1**a2;
    a1->release();
    a2->release();
}
*/
/*
void dh_drw(Camera *cam,double *Xv_km1_k,cv::Mat* yi,cv::Mat* zi,cv::Mat* a)
{
    cv::Mat* a1= new cv::Mat;
    *a1 = cv::Mat::zeros(2,3,CV_64FC1);
    cv::Mat* a2= new cv::Mat;
    *a2 = cv::Mat::zeros(3,3,CV_64FC1);
    
    dh_dhrl(cam,Xv_km1_k,yi,zi,a1);
    dhrl_drw(Xv_km1_k,a2);
    
    *a=*a1**a2;
    a1->release();
    a2->release();
}
*/
/*
void dh_dqwr(Camera *cam,double *Xv_km1_k,cv::Mat* yi,cv::Mat* zi,cv::Mat* a)
{
    cv::Mat* a1= new cv::Mat;
    *a1 = cv::Mat::zeros(2,3,CV_64FC1);
    cv::Mat* a2= new cv::Mat;
    *a2 = cv::Mat::zeros(3,3,CV_64FC1);
    
    dhu_dhrl(cam,Xv_km1_k,yi,a1);
    dhrl_dqwr(Xv_km1_k,yi,a2);
    
    *a=*a1**a2;
    a1->release();
    a2->release();
}
*/
//void dh_dxv(Camera *cam,double *Xv_km1_k,cv::Mat* yi,cv::Mat* zi,cv::Mat* Hi_features)
/*
void dh_dxv(Camera *cam, double *Xv_km1_k, double* yi, double* zi, double* Hi_features)
{
    cv::Mat* a1= new cv::Mat;
     *a1 = cv::Mat::zeros(2,3,CV_64FC1);
     cv::Mat* a2= new cv::Mat;
     *a2 = cv::Mat::zeros(2,3,CV_64FC1);
     cv::Mat* a3= new cv::Mat;
     cv::Mat::ones(2,6,CV_64FC1);
     
     dh_drw(cam,Xv_km1_k,yi,zi,a1);
     dh_dqwr(cam,Xv_km1_k,yi,zi,a2);
     
     a1->copyTo(Hi_features->colRange(0,2));
     a2->copyTo(Hi_features->colRange(3,5));
     a3->copyTo(Hi_features->colRange(6,11));
     
     a1->release();
     a2->release();
     a3->release();
}
*/
/*
void dh_dy(Camera *cam,double *Xv_km1_k,cv::Mat* yi,cv::Mat* zi,cv::Mat* a)
{
    cv::Mat* a1= new cv::Mat;
    *a1 = cv::Mat::zeros(2,3,CV_64FC1);
    cv::Mat* a2= new cv::Mat;
    *a2 = cv::Mat::zeros(3,3,CV_64FC1);
    
    dh_dhrl(cam,Xv_km1_k,yi,zi,a1);
    dhrl_dy(Xv_km1_k,a2);
    
    *a=*a1**a2;
    a1->release();
    a2->release();
}
*/
/*
void calculateHiCartesian(NSMutableArray *features_info,Camera *cam,double *Xv_km1_k,cv::Mat* yi, int step,cv::Mat* Hi)
{
    Feature *current=features_info[step-1];
    cv::Mat* zi= new cv::Mat;
    *zi = cv::Mat(1,2,CV_64FC1); //Please check later.
    
    int numberOfFeatures=(int)[features_info count];
    int *inverseDepthFeaturesIndex;
    int *cartesianFeaturesIndex;
    int indexOfInsertion;
    
    cv::Mat* a1= new cv::Mat;
    *a1 = cv::Mat(2,12,CV_64FC1);
    cv::Mat* a2= new cv::Mat;
    *a2 = cv::Mat(2,3,CV_64FC1);
    
    int i;
    int sum1=0,sum2=0;
    
    zi->at<double>(0,0)=current->h[0];
    zi->at<double>(0,1)=current->h[1];
    
    for (i=0;i<numberOfFeatures;i++)
    {
        Feature *current=features_info[i];
        
        if ([current->type isEqualToString:@"inversedepth"]) //if strncmp(features_info(j).type, 'inversedepth', 1)
        {
            *(inverseDepthFeaturesIndex+i)=1;
            sum1++;
        }
        else
            *(inverseDepthFeaturesIndex+i)=0;
        
        if ([current->type isEqualToString:@"cartesian"])
        {
            *(cartesianFeaturesIndex+i)=1;
            sum2++;
        }
        else
            *(cartesianFeaturesIndex+i)=0;
    }
    
    cv::Mat Hi_0=cv::Mat::zeros(2,13+3*sum1+6*sum2,CV_64FC1);
    
    //dh_dxv(cam,Xv_km1_k,yi,zi,a1);
    a1->copyTo(Hi_0.colRange(0,12));
    
    if (inverseDepthFeaturesIndex[step-2]==1)
        sum1--;
    if (cartesianFeaturesIndex[step-2]==1)
        sum2--;
    indexOfInsertion=13+3*sum2+6*sum1+1;
    
    dh_dy(cam,Xv_km1_k,yi,zi,a2);
    a2->copyTo(Hi_0.colRange(indexOfInsertion-1,indexOfInsertion+1));
    
    Hi_0.copyTo(*Hi);
}
*/
/*
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
                    y=&(x_features0[j]);
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
*/
//II.3)predictFeaturesAppearance
/*
void pred_pathch_fc(NSMutableArray* featureInfo,double *R_Wk,double *r_Wk,double *XYZ_w,Camera* cam,double *patchPred)
{
    
}
 */
/*
 void rotateWithDistFCC2C1(Camera* cam,double *uv_c2,double *R_c1c2,double *t_c1c2,double *n,int d,double *uv_c1)
 {
 
 }
 void rotateWithDistFCC1C2(Camera *cam,double *uv_c1,double *R_c2c1,double *t_c2c1,double *n,int d,double *uv_c2)
 {
 double *uv_c1_und;
 double *uv_c2_und;
 
 undistortFM(transpose());
 }
 */
/*
void predictFeaturesAppearance(NSMutableArray* features_info,Filter* filter,Camera* cam)
{
    double r_wc[3];
    double R_wc0[4];
    double R_wc[3][3];
    double *x_k_k_rest_of_features0;
    double *x_k_k_rest_of_features;
    double XYZ_w[3];
    double y[6];
    
    int i,j;
    int J=0;
    
    for (i=0;i<3;i++)
        r_wc[i]=filter->x_k_k[i];
    
    for (i=3;i<7;i++)
        R_wc0[i-3]=filter->x_k_k[i];
    //q2r((double *)R_wc0,(double *)R_wc);
    
    for (j=0;filter->x_k_k+j;j++)
        J++;
    
    for (j=13;j<J;j++)
    {
        x_k_k_rest_of_features0=&(filter->x_k_k[j]);
        x_k_k_rest_of_features0++;
    }
    
    for (i=0;i<(int)[features_info count];i++)
    {
        Feature *current=features_info[i];
        
        if ([current->type isEqualToString:@"cartesian"])
        {
            for (j=0;j<3;j++)
                XYZ_w[j]=*(x_k_k_rest_of_features0+i);
            
            for (j=3;j<J;j++)
            {
                x_k_k_rest_of_features=x_k_k_rest_of_features0+j;
                x_k_k_rest_of_features++;
            }
        }
        
        if ([current->type isEqualToString:@"inversedepth"])
        {
            for (j=0;j<6;j++)
                y[j]=*(x_k_k_rest_of_features0+i);
            
            for (j=6;j<J;j++)
            {
                x_k_k_rest_of_features=x_k_k_rest_of_features0+j;
                x_k_k_rest_of_features++;
            }
            
            //            inversedepth2cartesian((double *)y,(double *)XYZ_w);
        }
        
        if (sizeof(current->h)) //if ~isempty(features_info(i).h)
        {
            //            pred_pathch_fc(features_info,(double *)R_wc,(double *)R_wc,(double *)XYZ_w,cam,(double *)patchPred);
        }
    }
}
*/
//II.4)matching
/*
void matching(NSMutableArray* features_info,Camera* cam,double *im)
{
    double correlation_threshold=0.80;
    double chi_095_2=5.9915;
    double chi_099_2=9.2103;
    
    int i;
    
    for (i=0;i<(int)[features_info count];i++)
    {
        //features_info *current=features_info[i];
        
        //if (current->h) //if ~isempty(features_info(i_feature).h)
    }
    
}
 */
////////////////////////////////////////////////////////////////////////////////////////////////////
/////III.RANSACHypotheses
//III.1)generateStateVectorPattern
//Note: THIS IS THE CORRECT VERSION



//III.2)selectRandomMatch


//III.3)computeHypothesisSupportFast
/*

*/
 
//III.4)setAsMostSupportedHypothesis


////////////////////////////////////////////////////////////////////////////////////////////////////
/////MainLoop
void mapManagement(NSMutableArray* features_info,Filter* filter,Camera* cam,double *im,int minNUmberOfFeaturesInImage,int step)      //I
{
    
}
/*
void searchICMatches(NSMutableArray* featureInfo,Filter* filter,Camera* cam,double *im,int stepnumber)
//Correct later.
{
  int i;
  
  predictCameraMeasurements(featureInfo,filter,cam,stepnumber);
  
  calculateDerivatives(featureInfo,filter,cam);
  
  for (i=0;i<steps;i++) //for i=1:length(features_info)
  {
  if (featureInfo.h[i]) //~isempty(features_info(i).h)
  featureInfo.S[i]=featureInfo.H[i]*filter.p_k_km1*featureInfo.H[i]+featureInfo.R[i];
  //features_info(i).S = features_info(i).H*get_p_k_km1(filter)*features_info(i).H' + features_info(i).R
  }
  
  predictFeaturesAppearance(featureInfo,filter,cam);
  
  matching(featureInfo,cam,im);
  
}
*/
/*
void RANSACHypothesis(NSMutableArray* features_info,Filter* filter,Camera* cam)     //III
{
    double pAtLeastOneSpuriousFree=0.99;
    double threshold=filter->std_z;
    int numberOfHypotheses=1000;
    int maxHypothesisSupport=0;
    
    //Runtime error: Using a bunch of undeclared variables here
    //generate_state_vector_pattern(features_info,filter,z_id,z_euc)
    
    //   for (i=0;i<numberOfHypotheses;i++)
    //   {
    //       selectRandomMatch(featureInfo,double *position)
    //   }
}
*/
/*
 void normJac(double *q,cv::Mat J)
 {
 double r,x,y,z,j;
 
 r=*(q+0);
 x=*(q+1);
 y=*(q+2);
 z=*(q+3);
 
 j=powf(r*r+x*+y*y+z*z,-3/2);
 
 J.at<double>(0,0)=j*(x*x+y*y+z*z);
 J.at<double>(0,1)=j*(-r*x);
 J.at<double>(0,2)=j*(-r*y);
 J.at<double>(0,3)=j*(-r*z);
 J.at<double>(1,0)=j*(-x*r);
 J.at<double>(1,1)=j*(r*r+y*y+z*z);
 J.at<double>(1,2)=j*(-x*y);
 J.at<double>(1,3)=j*(-x*z);
 J.at<double>(2,0)=j*(-y*r);
 J.at<double>(2,1)=j*(-y*x);
 J.at<double>(2,2)=j*(r*r+x*x+z*z);
 J.at<double>(2,3)=j*(-y*z);
 J.at<double>(3,0)=j*(-z*r);
 J.at<double>(3,1)=j*(-z*x);
 J.at<double>(3,2)=j*(-z*y);
 J.at<double>(3,3)=r*r+x*x+y*y;
 }*/

void ComputeHypothesisSupport(double *xi, Camera *cam, double *state_vector_pattern, int nRows, double *z_id, int n_id, double *z_euc, double throshold, int hypothesis_support, int positions_li_inliers_id, int positions_li_inliers_euc)
{
    double u0 = cam->Cx;
    double v0 = cam->Cy;
    double f = cam->f;
    double ku = 1/cam->dx;
    double kv = 1/cam->dy;
    
    int i;
    
    if (z_id != NULL)
    {
        int ri_0_size = 0;
        int anglesi_0_size = 0;
        int rhoi_0_size = 0;
        
        for (i=0; i<nRows; i++)
        {
            if (state_vector_pattern[i*3+0] != 0)
                ri_0_size++;
            if (state_vector_pattern[i*3+1] != 0)
                anglesi_0_size++;
            if (state_vector_pattern[i*3+2] != 0)
                rhoi_0_size++;
        }
        
        double* ri_0 = (double *)malloc(ri_0_size*1*sizeof(double));
        double* anglesi_0 = (double *)malloc(anglesi_0_size*1*sizeof(double));
        double* rhoi_0 = (double *)malloc(rhoi_0_size*1*sizeof(double));
        
        for (i=0; i<ri_0_size; i++)
        {
            
        }
        
    }
}


void RescueHIInliers()
{
    
}

void EKFUpdateHIIlier()
{
    
}







