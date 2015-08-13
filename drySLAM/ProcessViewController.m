//
//  ProcessViewController.m
//  drySLAM
//
//  Created by Anthony Rodriguez on 3/8/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//

#import "ProcessViewController.h"
#import "SLAMModelAdapter.h"

extern NSMutableArray* fileNamesArray;
extern NSString* filePath;

@interface ProcessViewController ()

@end

@implementation ProcessViewController

- (void)viewDidLoad {
    [super viewDidLoad];
    // Do any additional setup after loading the view from its nib.
    [SLAMModelAdapter updateFileNamesArray];
    [self.tableView setDelegate:self];
    [self.tableView reloadData];
    [self.view setNeedsDisplay];
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

-(IBAction)processButtonPressed:(id)sender{
    [self.processButton setTitle:@"Abort" forState:UIControlStateNormal];
    NSString* status = [NSString stringWithFormat:@"Processing %@...", fileNamesArray[0]];
    [self.feedbackLabel setText:status];
    status = nil;
    //call Model through SLAMModelAdapter here
    //TODO: Set to a busy state and launch processing thread
    [SLAMModelAdapter processDataWithRANSAC];
    [self processingFinished];
}

-(void)processingFinished{
    //called when processing thread is complete
    [self.processButton setTitle:@"Process" forState:UIControlStateNormal];
    [self.feedbackLabel setText:@"Processing complete!"];
    [self.tableView reloadData];
    [self.view setNeedsDisplay];
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
