//
//  RecordViewController.m
//  drySLAM
//
//  Created by Anthony Rodriguez on 3/8/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//

#import "RecordViewController.h"
//#import <opencv2/highgui/cap_ios.h>
#import <Accelerate/Accelerate.h>
#import <AssetsLibrary/AssetsLibrary.h>
#import <CoreGraphics/CoreGraphics.h>
#import <CoreImage/CoreImage.h>
#import <CoreMedia/CoreMedia.h>
#import <CoreVideo/CoreVideo.h>
#import <QuartzCore/QuartzCore.h>
#import <UIKit/UIKit.h>
#import <Foundation/Foundation.h>
#import <AVFoundation/AVFoundation.h>

#import "SLAMModelAdapter.h"


//Calibrate for which service we are using
#define OPENCV 0
#define AVFRAMEWORK 1

#define FPS 30

//using namespace cv;

extern NSMutableArray* fileNamesArray;
extern NSString* filePath;

@interface RecordViewController ()

@end

@implementation RecordViewController

#if (AVFRAMEWORK)
@synthesize previewLayer;
#endif

- (void)viewDidLoad {
    [super viewDidLoad];
#if (OPENCV)
    // Do any additional setup after loading the view from its nib.
    self.videoCamera = [[CvVideoCamera alloc]initWithParentView:self.imageView];
    self.videoCamera.defaultAVCaptureDevicePosition = AVCaptureDevicePositionBack;
    //TODO Figure out what resolution,fps,colorspace we need:
    self.videoCamera.defaultAVCaptureSessionPreset = AVCaptureSessionPreset640x480;
    self.videoCamera.defaultAVCaptureVideoOrientation = AVCaptureVideoOrientationLandscapeRight;
    self.videoCamera.defaultFPS = 30;
    self.videoCamera.grayscaleMode = YES;
    self->cameraIsCapturing = NO;
    self.videoCamera.delegate = self;
    [self.imageView setHidden:NO];
    [self.AVView setHidden:YES];
    [self.imageView setNeedsDisplay];
#endif
    
    
#if (AVFRAMEWORK)
    [self.imageView setHidden:YES];
    [self.AVView setHidden:NO];
    self->captureSession = [[AVCaptureSession alloc]init];
    AVCaptureDevice* videoDevice = [AVCaptureDevice defaultDeviceWithMediaType:AVMediaTypeVideo];
    if (videoDevice){
        NSError* error;
        videoInputDevice = [AVCaptureDeviceInput deviceInputWithDevice:videoDevice error:&error];
        if (!error){
            if ([captureSession canAddInput:videoInputDevice])
                [captureSession addInput:videoInputDevice];
            else
                NSLog(@"Couldnt add video input");
        }
        else{
            NSLog(@"Couldnt create video input");
        }
    }
    else{
        NSLog(@"Couldnt create video capture device");
    }
    [videoDevice setActiveVideoMaxFrameDuration:CMTimeMake(1,FPS)];
    [videoDevice setActiveVideoMinFrameDuration:CMTimeMake(1,FPS)];
    [videoDevice unlockForConfiguration];
    AVCaptureInput* cameraInput = [[AVCaptureDeviceInput alloc] initWithDevice:videoDevice error:nil];
    AVCaptureVideoDataOutput* videoOutput = [[AVCaptureVideoDataOutput alloc]init];
    dispatch_queue_t captureQueue = dispatch_queue_create("captureQueue",DISPATCH_QUEUE_SERIAL);
    
    [videoOutput setSampleBufferDelegate:self queue:captureQueue];
    
    videoOutput.videoSettings = [NSDictionary dictionaryWithObjectsAndKeys:[NSNumber numberWithUnsignedInt:kCVPixelFormatType_420YpCbCr8BiPlanarFullRange],(id)kCVPixelBufferPixelFormatTypeKey,nil];
    [captureSession setSessionPreset:AVCaptureSessionPresetMedium];
    if ([captureSession canAddInput:cameraInput])
        [captureSession addInput:cameraInput];
    if ([captureSession canAddOutput:videoOutput])
        [captureSession addOutput:videoOutput];
    
    [self setPreviewLayer:[[AVCaptureVideoPreviewLayer alloc] initWithSession:captureSession]];
    [self.previewLayer setVideoGravity:AVLayerVideoGravityResizeAspectFill];
    movieFileOutput = [[AVCaptureMovieFileOutput alloc]init];
    Float64 totalSeconds = 60;
    int32_t preferredTimeScale = FPS;
    CMTime maxDuration = CMTimeMakeWithSeconds(totalSeconds, preferredTimeScale);
    movieFileOutput.maxRecordedDuration = maxDuration;
    movieFileOutput.minFreeDiskSpaceLimit = 1024*1024;
    if ([captureSession canAddOutput:movieFileOutput])
        [captureSession addOutput:movieFileOutput];
    [self CameraSetOutputProperties];
    if ([captureSession canSetSessionPreset:AVCaptureSessionPreset640x480]){
        [captureSession setSessionPreset:AVCaptureSessionPreset640x480];
        NSLog(@"Low-res preset set");
    }
    else
        [captureSession setSessionPreset:AVCaptureSessionPresetMedium];
    self.previewLayer.frame = self.AVView.frame;
    [self.AVView.layer addSublayer:previewLayer];
    CGRect bounds = self.AVView.bounds;
    previewLayer.bounds = bounds;
    previewLayer.position = CGPointMake(CGRectGetMidX(bounds), CGRectGetMidY(bounds));
    previewLayer.orientation = AVCaptureVideoOrientationLandscapeRight;
#endif
    
    [self.startStopButton setTitle:@"Start" forState:UIControlStateNormal];
    [self.feedbackLabel setText:@"Waiting..."];
}

-(void)viewWillDisappear:(BOOL)animated{
#if (OPENCV)
    [self.videoCamera stop];
#endif
#if (AVFRAMEWORK)
    [captureSession stopRunning];
    captureSession = nil;
    movieFileOutput = nil;
    videoInputDevice = nil;
#endif
    [self.feedbackLabel setText:@"Waiting..."];
}

-(void)CameraSetOutputProperties{
    AVCaptureConnection* captureConnection = [movieFileOutput connectionWithMediaType:AVMediaTypeVideo];
    if ([captureConnection isVideoOrientationSupported]){
        AVCaptureVideoOrientation orientation = AVCaptureVideoOrientationLandscapeRight;
        [captureConnection setVideoOrientation:orientation];
        NSLog(@"Orientation set");
    }
}


-(IBAction)startStopButtonPressed:(id)sender{
    if (self->cameraIsCapturing==YES){
#if (OPENCV)
        [self.videoCamera stop];
       // videoWriter.release();
#endif
#if (AVFRAMEWORK)
        [captureSession stopRunning];
        [movieFileOutput stopRecording];
        [SLAMModelAdapter updateFileNamesArray];
#endif
        self->cameraIsCapturing = NO;
        [self.startStopButton setTitle:@"Start" forState:UIControlStateNormal];
            [self.feedbackLabel setText:@"Waiting..."];
    }
    else{
#if (OPENCV)
        [self.videoCamera start];
        //videoWriter = cv::VideoWriter("captureSession.avi", CV_FOURCC('M','J','P','G'), FPS, cvSize(640,480),true);
#endif
#if (AVFRAMEWORK)
        [captureSession startRunning];
        //Create temporary URL to record to
        NSString *outputPath = [[NSString alloc] initWithFormat:@"%@%@", filePath, @"/recordedVideo.mov"];
        NSURL *outputURL = [[NSURL alloc] initFileURLWithPath:outputPath];
        NSLog(@"Record video to %@", outputPath);
        NSFileManager *fileManager = [NSFileManager defaultManager];
        if ([fileManager fileExistsAtPath:outputPath])
        {
            NSError *error;
            if ([fileManager removeItemAtPath:outputPath error:&error] == NO)
            {
                //Error - handle if requried
            }
        }
        outputPath = nil;
        //Start recording
        [movieFileOutput startRecordingToOutputFileURL:outputURL recordingDelegate:self];
        outputURL = nil;
#endif
        self->cameraIsCapturing = YES;
        [self.startStopButton setTitle:@"Stop" forState:UIControlStateNormal];
        NSString* status = [NSString stringWithFormat:@"Recording to %@", (NSString*)fileNamesArray[0]];
        
        [self.feedbackLabel setText:status];
    }
}

-(void)captureOutput:(AVCaptureFileOutput *)captureOutput didFinishRecordingToOutputFileAtURL:(NSURL *)outputFileURL fromConnections:(NSArray *)connections error:(NSError *)error{
    NSLog(@"didFinishRecording");
    if ([error code] != noErr){
        NSLog(@"Error occured in recording");
    }
}

//Unused unless using opencv
-(void)processImage:(cv::Mat &)image{
    
}

-(IBAction)backButtonPressed:(id)sender{
    [self dismissViewControllerAnimated:YES completion:NULL];
    [self.presentingViewController dismissViewControllerAnimated:YES completion:NULL];
}

- (void)didReceiveMemoryWarning {
    [super didReceiveMemoryWarning];
    // Dispose of any resources that can be recreated.
}

/*
#pragma mark - Navigation

// In a storyboard-based application, you will often want to do a little preparation before navigation
- (void)prepareForSegue:(UIStoryboardSegue *)segue sender:(id)sender {
    // Get the new view controller using [segue destinationViewController].
    // Pass the selected object to the new view controller.
}
*/

@end
