import sys
import os
import shutil

def main():
    datasets = ['1000M1', '1000M4']
    reps = 10
    
    # first 16S.M
    path = '16S.M/'
    files = os.listdir(path)
    for f in files:
        if f.find('pastajob1') != -1:
            ind = f.find('pastajob1') + len('pastajob1')
            #print(path+'pastajob'+f[ind:])
            os.rename(path+f, path+'pastajob'+f[ind:])

    #for i in range(0, reps):
    #    path = 

if __name__ == "__main__":
    main()
