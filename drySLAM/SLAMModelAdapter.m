//
//  SLAMModelAdapter.m
//  drySLAM
//
//  Created by Anthony Rodriguez on 3/13/15.
//  Copyright (c) 2015 ECE473 573 Common. All rights reserved.
//

#import <Foundation/Foundation.h>
#include "SLAMModelAdapter.h"


@interface SLAMModelAdapter ()

@end

@implementation SLAMModelAdapter


+(void)updateFileNamesArray{
    fileNamesArray = nil;
    fileNamesArray = [[NSMutableArray alloc]init];
    //TODO: Read files from directory
    [SLAMModelAdapter setFilePath];
    NSArray* directoryContent = [[NSFileManager defaultManager] contentsOfDirectoryAtPath:filePath error:nil];
    int nFiles = (int)[directoryContent count];
    for (int i=0; i<nFiles; ++i){
        NSString* addStr = [(NSString*)[directoryContent objectAtIndex:i] copy];
        [fileNamesArray addObject:addStr];
    }
    [fileNamesArray addObject:@"dummyFile.SLAM"];
}

+(void)processDataWithRANSAC{
    //TODO: Atually call model here
    for (int i=0; i<10000; ++i);
    
    //Write a data file that can be visibly shown:
    NSString* fileName = [NSString stringWithFormat:@"testFile.SLAM"];
    
    NSString* fullPath = [filePath stringByAppendingPathComponent:fileName];
    if (![[NSFileManager defaultManager]fileExistsAtPath:fullPath]){
        [[NSFileManager defaultManager] createFileAtPath:fullPath contents:nil attributes:nil];
        NSString* header = @"%There's nothing here!";
        [header writeToFile:fullPath atomically:NO encoding:NSUTF8StringEncoding error:nil];
    }
    else{
        NSLog(@"Attempt to write duplicate file");
    }
    [self updateFileNamesArray];
}

+(void)setFilePath{
    NSArray* pathList = NSSearchPathForDirectoriesInDomains(NSDocumentDirectory, NSUserDomainMask, YES);
    filePath = (NSString*)[pathList objectAtIndex:0];
}

@end