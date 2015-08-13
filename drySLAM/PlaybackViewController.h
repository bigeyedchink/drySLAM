//
//  PlaybackViewController.h
//  drySLAM
//
//  Created by Anthony Rodriguez on 3/8/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//


#import <UIKit/UIKit.h>
@import MediaPlayer;
@import GLKit;

@interface graphicsView : GLKView{
    CAEAGLLayer* _eaglLayer;
    EAGLContext* _context;
    GLuint _positionSlot;
    GLuint _colorSlot;
    GLuint _colorRenderBuffer;
}

@property(retain) EAGLContext* glContext;
-(id)initWithFrame:(CGRect)frame;
+(Class)layerClass;
-(void)setupLayer;
-(void)setupContext;

@end

/*********************************************************/

@interface PlaybackViewController : UIViewController <UITableViewDelegate, UITableViewDataSource>{
    //graphicsView* _glView;
}

@property(weak,nonatomic) IBOutlet UIButton* backButton;
@property(weak,nonatomic) IBOutlet UITableView* tableView;

//Left view
@property (strong, nonatomic) MPMoviePlayerController *moviePlayer;
@property(weak,nonatomic) IBOutlet UIView* mpvcView;

@property(weak,nonatomic) IBOutlet UIView* graphicsWindow;
@property(retain,nonatomic) IBOutlet graphicsView* glView;


-(IBAction)backButtonPressed:(id)sender;
-(IBAction)videoContextPressed:(id)sender;
-(IBAction)graphicsContextPressed:(id)sender;

@end

