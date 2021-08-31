import os
import sys

from shutil import copy
from tabulate import tabulate
class ShellPath :
    '''
    Designate the position of shell scripts and 
    check the shell scripts registered with gmd.
    '''
    def __init__(self) :
        self.shellpath = "%s/shell"%(os.path.split(os.path.dirname(__file__))[0])
        self.files = [f for f in os.listdir(self.shellpath)]

    def register_shell(self,shellpath):
        '''
        Save shell file
        '''
        if os.path.isfile(shellpath):
            filename = os.path.basename(shellpath)
            filepath = os.path.abspath("%s/%s"%(self.shellpath,filename))
            if os.path.isfile(filepath):
                print("%s is already enrolled\n"%(filename))
            else :
                copy(shellpath,filepath)
                print("%s Successfully Save "%(filepath))
        else :
            print("%s isn't file"%(shellpath))

    def check(self):
        '''
        Check shell file registered 
        '''
        shell = [] ; numberlist =[] ;number = 1
        for f in self.files :
            if f.split(".")[-1] == 'sh' :
                shell.append([number, f])
                numberlist.append(number)
                number += 1
        print("\n\n",tabulate(shell, tablefmt='grid'),end='\n\n') 
        return shell, numberlist

    def remove(self) :
        '''
        Remove the shell file registered
        '''
        shell, number = self.check()

        while True :
            r = int(input("Please enter the script name for removal >> "))
            if r in number :
                break
        os.remove("{}/{}".format(self.shellpath, shell[r-1][-1]))
        print("{} scirpt is deleted".format(shell[r-1][-1]))