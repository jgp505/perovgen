import os
import sys
import yaml

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

    def generateshell(self, shell) :

        while True :
            regist = input(str("Please enter Y to register shell script, otherwise enter N >> "))
            if regist == 'Y' :
                path = self.shellpath + "/{}.sh".format(shell['shell_name'])
                break
            elif regist == 'N' :
                path = "{}.sh".format(shell['shell_name'])
                break
            else :
                print("Please Enter Y or N")

        with open(path,'w') as sh :
            sh.write("#!/bin/sh\n")
            sh.write("# control options #\n")
            sh.write("#PBS -N {} \n".format(shell['shell_name']))
            sh.write("#PBS -l nodes={}:ppn={}:{}\n".format(shell['node'],shell['ppn'],shell['node_name']))
            sh.write("########\n")
            sh.write("#PBS -q {}\n".format(shell['node_name']))
            sh.write("#PBS -o out.log\n")
            sh.write("#PBS -j oe\n")
    
            sh.write("\n# PATH & EXE\n")
            sh.write("EXE='{}'\n".format(shell['vasp_std']))
            sh.write("#EXE='{}'\n".format(shell['vasp_ncl']))
            sh.write("#EXE='{}'\n".format(shell['vasp_gam']))
    
            sh.write("\n#\n")
            sh.write("NUMBER=`cat $PBS_NODEFILE | wc -l`\n")
            sh.write("cd $PBS_O_WORKDIR\n")
    
            sh.write("\n# run \n")
            sh.write("echo job started at `date` >> time\n")
            sh.write("{} -np $NUMBER -machinefile $PBS_NODEFILE $EXE  > $PBS_JOBNAME.out\n".format(shell['mpi_command']))
            sh.write("echo job ended at `date` >> time\n")
    
        sh.close()


                