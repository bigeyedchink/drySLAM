//
//  PlaybackViewController.m
//  drySLAM
//
//  Created by Anthony Rodriguez on 3/8/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//

#import "PlaybackViewController.h"
#import "SLAMModelAdapter.h"
#import "VisualGraphicsViewController.h"



typedef struct{
    float Position[3];
    float Color[4];
} Vertex;

const Vertex Vertices[] = {
    {{1,-1,0}, {1,0,0,1}},
    {{1,1,0}, {0,1,0,1}},
    {{-1,1,0}, {0,0,1,1}},
    {{-1,-1,0}, {0,0,0,1}}
};

const GLubyte Indices[] = {
    0, 1, 2,
    2, 3, 0
};

@interface graphicsView (){
    //uses cocos3d
    GLuint _projectionUniform;
}
@end

@implementation graphicsView

-(void)setupVBOs {
    GLuint vertexBuffer;
    glGenBuffers(1, &vertexBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(Vertices), Vertices, GL_STATIC_DRAW);
    GLuint indexBuffer;
    glGenBuffers(1, &indexBuffer);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBuffer);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(Indices), Indices, GL_STATIC_DRAW);
}

-(GLuint)compileShader:(NSString*)shaderName withType:(GLenum)shaderType{
    NSString* shaderPath = [[NSBundle mainBundle] pathForResource:shaderName ofType:@"glsl"];
    NSError* error;
    NSString* shaderString = [NSString stringWithContentsOfFile:shaderPath encoding:NSUTF8StringEncoding error:&error];
    if (!shaderString){
        NSLog(@"Error loading shader: %@", error.localizedDescription);
        exit(1);
    }
    GLuint shaderHandle = glCreateShader(shaderType);
    
    const char* shaderStringUTF8 = [shaderString UTF8String];
    int shaderStringLength = (int)[shaderString length];
    glShaderSource(shaderHandle, 1, &shaderStringUTF8, &shaderStringLength);
    glCompileShader(shaderHandle);
    GLint compileSuccess;
    glGetShaderiv(shaderHandle, GL_COMPILE_STATUS, &compileSuccess);
    if (compileSuccess==GL_FALSE){
        GLchar messages[256];
        glGetShaderInfoLog(shaderHandle, sizeof(messages), 0, &messages[0]);
        NSString* messageString = [NSString stringWithUTF8String:messages];
        NSLog(@"%@", messageString);
        exit(1);
    }
    
    
    return shaderHandle;
}

-(void)compileShaders{
    GLuint vertexShader = [self compileShader:@"SimpleVertex" withType:GL_VERTEX_SHADER];
    GLuint fragmentShader = [self compileShader:@"SimpleFragment" withType:GL_FRAGMENT_SHADER];
    GLuint programHandle = glCreateProgram();
    glAttachShader(programHandle, vertexShader);
    glAttachShader(programHandle, fragmentShader);
    glLinkProgram(programHandle);
    GLint linkSuccess;
    glGetProgramiv(programHandle, GL_LINK_STATUS, &linkSuccess);
    if (linkSuccess==GL_FALSE){
        GLchar messages[256];
        glGetProgramInfoLog(programHandle, sizeof(messages), 0, &messages);
        NSString* messageString = [NSString stringWithUTF8String:messages];
        NSLog(@"%@", messageString);
        exit(1);
    }
    glUseProgram(programHandle);
    _positionSlot = glGetAttribLocation(programHandle, "Position");
    _colorSlot = glGetAttribLocation(programHandle, "SourceColor");
    glEnableVertexAttribArray(_positionSlot);
    glEnableVertexAttribArray(_colorSlot);
    
    _projectionUniform = glGetUniformLocation(programHandle, "Projection");
}

-(id)initWithFrame:(CGRect)frame{
    //self = [super init];
    NSLog(@"glView initializing...");
    self = [super initWithFrame:frame];
    if (self){
        [self setupLayer];
        [self setupContext];
        [self setupRenderBuffer];
        [self setupFrameBuffer];
        [self compileShaders];
        [self setupVBOs];
        [self render];
    }
    NSLog(@"Initialization complete");
    return self;
}

+(Class)layerClass{
    return [CAEAGLLayer class];
}

-(void)setupLayer{
    _eaglLayer = (CAEAGLLayer*) self.layer;
    _eaglLayer.opaque = YES;
    NSLog(@"Layer set up");
}

-(void)setupContext{
    EAGLRenderingAPI api = kEAGLRenderingAPIOpenGLES2;
    _context = [[EAGLContext alloc]initWithAPI:api];
    if (!_context){
        NSLog(@"Failed to initialize OpenGL ES 2.0 context");
        exit(1);
    }
    if (![EAGLContext setCurrentContext:_context]){
        NSLog(@"Failed to set currentOpenGL context");
        exit(1);
    }
    NSLog(@"Context set up");
}

-(void)setupRenderBuffer{
    glGenRenderbuffers(1, &_colorRenderBuffer);
    glBindRenderbuffer(GL_RENDERBUFFER, _colorRenderBuffer);
    [_context renderbufferStorage:GL_RENDERBUFFER fromDrawable:_eaglLayer];
    NSLog(@"Renderbuffer set up");
}

-(void)drawRect:(CGRect)rect{
    glClearColor(0.f,0.f,0.1f,1.0f);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
}

-(void)setupFrameBuffer{
    GLuint framebuffer;
    glGenFramebuffers(1, &framebuffer);
    glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, _colorRenderBuffer);
    NSLog(@"Framebuffer set up");
}

-(void)render{
    glClearColor(0, 104.0/255.0, 55.0/255.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);
    //Projection information (move scope in 3-space)
   // CC3GLMatrix* projection = [CC3GLMatrix matrix];
    float h = 4.0f * self.frame.size.height/self.frame.size.width;
    NSLog(@"H:%lf\tW:%lf", self.frame.size.height, self.frame.size.width);
    //[projection populateFromFrustumLeft:-1 andRight:1 andBottom:-1 andTop:1 andNear:5 andFar:10];
    //pass data to vertex shader (converts matrix into an openGL array)
   // glUniformMatrix4fv(_projectionUniform, 1, 0, projection.glMatrix);
    
    glViewport(0, 0, self.frame.size.width, self.frame.size.height);
    glVertexAttribPointer(_positionSlot, 3, GL_FLOAT,GL_FALSE, sizeof(Vertex), 0);
    glVertexAttribPointer(_colorSlot, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (GLvoid*)(sizeof(float)*3));
    glDrawElements(GL_TRIANGLES, sizeof(Indices)/sizeof(Indices[0]), GL_UNSIGNED_BYTE, 0);
    
    [_context presentRenderbuffer:GL_RENDERBUFFER];
    NSLog(@"Matrix rendered");
}

-(void) dealloc{
    _context = nil;
}

@end

/*******************************************************************************************************/

@interface PlaybackViewController ()

@end

@implementation PlaybackViewController

@synthesize glView=_glView;
VisualGraphicsViewController* vgvc;

- (void)viewDidLoad {
    [super viewDidLoad];
    NSLog(@"ViewDIdLoad");
    // Do any additional setup after loading the view from its nib.
    [SLAMModelAdapter updateFileNamesArray];
    [self.tableView setDelegate:self];
    CGRect bounds = [self.graphicsWindow bounds];
    NSLog(@"Calling glView init...");
    self.glView = [[graphicsView alloc] initWithFrame:bounds];
    [self.glView setCenter:CGPointMake(CGRectGetMidX(bounds), CGRectGetMidY(bounds))];
    [self.graphicsWindow addSubview:_glView];
    [self.tableView reloadData];
    [self.view setNeedsDisplay];
}

-(void) drawRect:(CGRect)rect{

}

- (NSInteger)tableView:(UITableView *)tableView numberOfRowsInSection:(NSInteger)section{
    [SLAMModelAdapter updateFileNamesArray];
    return [fileNamesArray count];
}

//Allow scrolling through files
- (UITableViewCell *)tableView:(UITableView *)tableView cellForRowAtIndexPath:(NSIndexPath *)indexPath
{
    static NSString *simpleTableIdentifier = @"SimpleTableCell";
    
    UITableViewCell *cell = [tableView dequeueReusableCellWithIdentifier:simpleTableIdentifier];
    
    if (cell == nil) {
        cell = [[UITableViewCell alloc] initWithStyle:UITableViewCellStyleDefault reuseIdentifier:simpleTableIdentifier];
    }
    [SLAMModelAdapter updateFileNamesArray];
    //Save this for later:
    /*fileNames = nil;
     fileNames = [[NSMutableArray alloc]init];
     NSArray* directoryContent = [[NSFileManager defaultManager]contentsOfDirectoryAtPath:filePath error:nil];
     numFiles = (int)[directoryContent count];
     for(int i=0; i<numFiles; ++i){
     NSString* addStr = [(NSString*)[directoryContent objectAtIndex:i] copy];
     [fileNames addObject:addStr];
     }*/
    cell.textLabel.text = (NSString*)[fileNamesArray objectAtIndex:indexPath.row];
    return cell;
}

- (void) tableView: (UITableView*) tableView didSelectRowAtIndexPath: (NSIndexPath *) indexPath
{
    NSIndexPath* idxPath = [self.tableView indexPathForSelectedRow];
    int idx = (int)idxPath.row;
    NSString* fileName = fileNamesArray[idx];
    NSString* key = [NSString stringWithFormat:@".mov"];
    NSRange range = [fileName  rangeOfString: key options: NSCaseInsensitiveSearch];
    if (range.location == NSNotFound){
        NSLog(@"Selected file isn't video! :(");
        return;
    }
    else
        NSLog(@"Video %@ selected", fileName);
    NSString* fullFilePath = [NSString stringWithFormat:@"%@/%@", filePath, fileName];
    NSURL *videoURL = [[NSURL alloc] initFileURLWithPath:fullFilePath];
    _moviePlayer =  [[MPMoviePlayerController alloc] initWithContentURL:videoURL];
    [[NSNotificationCenter defaultCenter] addObserver:self
                                             selector:@selector(moviePlayBackDidFinish:)
                                                 name:MPMoviePlayerPlaybackDidFinishNotification
                                               object:_moviePlayer];
    _moviePlayer.controlStyle = MPMovieControlStyleDefault;
    _moviePlayer.shouldAutoplay = YES;
    [_moviePlayer.view setFrame:self.mpvcView.frame];
    CGRect bounds = self.mpvcView.bounds;
    _moviePlayer.view.bounds = bounds;
    [_moviePlayer.view setCenter:CGPointMake(CGRectGetMidX(bounds), CGRectGetMidY(bounds))];
    [self.mpvcView addSubview:_moviePlayer.view];
    [_moviePlayer setFullscreen:NO animated:YES];
}

- (void) moviePlayBackDidFinish:(NSNotification*)notification {
    MPMoviePlayerController *player = [notification object];
    [[NSNotificationCenter defaultCenter]
     removeObserver:self
     name:MPMoviePlayerPlaybackDidFinishNotification
     object:player];
    
    if ([player
         respondsToSelector:@selector(setFullscreen:animated:)])
    {
        [player.view removeFromSuperview];
    }
}

-(IBAction)backButtonPressed:(id)sender{
    [self dismissViewControllerAnimated:YES completion:NULL];
    [self.presentingViewController dismissViewControllerAnimated:YES completion:NULL];
}

-(IBAction)graphicsContextPressed:(id)sender{
    NSLog(@"Graphics context pressed");
    /*CGRect bounds = [self.view bounds];
    NSLog(@"Calling glView init...");
    self.glView = [[graphicsView alloc] initWithFrame:bounds];
    [self.glView setCenter:CGPointMake(CGRectGetMidX(bounds), CGRectGetMidY(bounds))];
    self.view = self.graphicsWindow;*/
    vgvc = [[VisualGraphicsViewController alloc]init];
    [self presentViewController:vgvc animated:YES completion:NULL];
    
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
