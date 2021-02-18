#import pandas as pd
import os
import sys
import math

CHKPT_NUM = 10000000
#10000000



def main(argv):

    if len(argv) - 1 != 2:
        print("eva.py <bowtie.sam> <myfile.txt>")

    print("Import from the file: ", argv[1], argv[2])

    
    #myFile = open(argv[2], 'r')

    bowtieLineCnt = 0
    startData = False


    column_names = ["chr", "loc"]
    bowtieDict = {}
    #bowtiedf = pd.DataFrame()

    bowiteFile = open(argv[1], 'r')

    for line in bowiteFile:
        #if bowtieLineCnt == 200:
            #break

        if startData == True:
            bowtieLineCnt += 1

        if bowtieLineCnt % CHKPT_NUM == 0 and bowtieLineCnt != 0:
            print("check pt", int(bowtieLineCnt/CHKPT_NUM), "Line:", bowtieLineCnt)

        

        line = line.split("\n")
        lineContent = line[0].split("\t")


        if lineContent[0] == "@PG":
            startData = True
            continue
        elif startData == False:
            continue

        # Start processing
        bowtieDict[lineContent[0]] = [lineContent[2], lineContent[3]]


    bowiteFile.close()



    print("\nGetting to process my file")
    
    tmpcnt = 0
    curName = ''

    myLineCnt = 0

    myCorAlign = 0
    myIncorAlign = 0
    totalProc = 0

    CantFindInBowtie = 0 

    for line in open(argv[2], 'r'):

        #if totalProc == 1:
            #break

        line = line.split("\n")

        myLineCnt = myLineCnt + 1
        if myLineCnt % 1000000 == 0:
            print("Process check pt", int(myLineCnt/1000000))

        if line[0][0] == '@':
            curName = line[0][1:]
            totalProc = totalProc + 1

        else:
            lineContent = line[0].split(" ")

            bowtieEntry = None
            furtherProcess =  False
            bowtieEntry = []
            #print(curName)
            if curName in bowtieDict:
                bowtieEntry = bowtieDict[curName]

                if bowtieEntry[0] == '*':
                    furtherProcess = False
                    CantFindInBowtie = CantFindInBowtie + 1
                else:
                    furtherProcess = True
                
            else:
                CantFindInBowtie = CantFindInBowtie + 1



            canAlign = False
            if furtherProcess == True:
                #print("Mine: ", lineContent)
                for i in range(0, len(lineContent), 2):
                    if lineContent[i] == bowtieEntry[0] and lineContent[i+1][:-2] == bowtieEntry[1][:-2]:
                        myCorAlign = myCorAlign + 1
                        canAlign = True
                        #print("TTTTTTTTUre")
                        break
                
                if canAlign == False:
                    myIncorAlign = myIncorAlign + 1
                    #print("FFFFFFFFFFFFFalse")


            #print()

    

    print("Correct Align: ", myCorAlign)
    print("Incorrect Align:", myIncorAlign)
    print("Cant Find In Bowtie: ", CantFindInBowtie)
    print("Total: ", totalProc)

        

        



if __name__ == '__main__':
    main(sys.argv)



"""
def main(argv):

    if len(argv) - 1 != 2:
        print("eva.py <bowtie.sam> <myfile.txt>")

    print("Import from the file: ", argv[1], argv[2])

    
    #myFile = open(argv[2], 'r')

    bowtieLineCnt = 0
    startData = False


    column_names = ["chr", "loc"]
    bowtieDict = {}
    #bowtiedf = pd.DataFrame()

    bowiteFile = open(argv[1], 'r')

    for line in bowiteFile:
        #if bowtieLineCnt == 200:
            #break

        if startData == True:
            bowtieLineCnt += 1

        if bowtieLineCnt % CHKPT_NUM == 0 and bowtieLineCnt != 0:
            #if bowtieLineCnt == CHKPT_NUM:
                #bowtiedf = pd.DataFrame.from_dict(bowtieDict, orient='index', columns=column_names)

                #print(bowtiedf)

            #else:
                #print(bowtieDict)
                #newdf = pd.DataFrame.from_dict(bowtieDict, orient='index',  columns=column_names)
                #bowtiedf = pd.concat([bowtiedf, newdf])

                #print(bowtiedf)
                
            
            #bowtieDict = {}

            print("check pt", int(bowtieLineCnt/CHKPT_NUM), "Line:", bowtieLineCnt)

        

        line = line.split("\n")
        lineContent = line[0].split("\t")


        if lineContent[0] == "@PG":
            startData = True
            continue
        elif startData == False:
            continue

        # Start processing
        bowtieDict[lineContent[0]] = [lineContent[2], lineContent[3]]
        
        if startData == True:
            print(bowtieDict)
            break 
        

    bowiteFile.close()

    #print("\nGenerate final df")
    # the final one
    #if bowtieDict:
        #newdf = pd.DataFrame.from_dict(bowtieDict, orient='index', columns=column_names)
        #bowtiedf = pd.concat([bowtiedf, newdf])

    print("\nGetting to process my file")
    
    tmpcnt = 0
    curName = ''

    myLineCnt = 0

    myCorAlign = 0
    myIncorAlign = 0
    totalProc = 0

    CantFindInBowtie = 0 

    for line in open(argv[2], 'r'):

        #if totalProc == 1:
            #break

        line = line.split("\n")

        myLineCnt = myLineCnt + 1
        if myLineCnt % 1000000 == 0:
            print("Process check pt", int(myLineCnt/1000000))

        if line[0][0] == '@':
            curName = line[0][1:]
            totalProc = totalProc + 1

        else:
            lineContent = line[0].split(" ")

            bowtieEntry = None
            furtherProcess =  False
            bowtieEntry = []
            print(curName)
            if curName in bowtieDict:
                furtherProcess = True
                bowtieEntry = bowtieDict[curName]
                
            else:
                CantFindInBowtie = CantFindInBowtie + 1
            
            try:
                bowtieEntry = bowtiedf.loc[curName]
                furtherProcess = True

                #print(bowtieEntry)
            except KeyError:
                CantFindInBowtie = CantFindInBowtie + 1
            

            canAlign = False
            if furtherProcess == True:
                #print("Mine: ", lineContent)
                for i in range(0, len(lineContent), 2):
                    if lineContent[i] == bowtieEntry[0] and lineContent[i+1][:-2] == bowtieEntry[1][:-2]:
                        myCorAlign = myCorAlign + 1
                        canAlign = True
                        #print("TTTTTTTTUre")
                        break
                
                if canAlign == False:
                    myIncorAlign = myIncorAlign + 1
                    #print("FFFFFFFFFFFFFalse")


            #print()

    

    print("Correct Align: ", myCorAlign)
    print("Incorrect Align:", myIncorAlign)
    print("Cant Find In Bowtie: ", CantFindInBowtie)
    print("Total: ", totalProc)
"""