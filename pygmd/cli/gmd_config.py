import os
import sys

from perovgen.pygmd.shell import ShellPath
import yaml

def argsregistershell(args):
    sp = ShellPath()
    shellname = os.path.abspath(args[0])
    sp.register_shell(shellname)

def argsshellcheck(args):
    ShellPath().check()

def config(args):
    if args.shell :
        argsregistershell(args.shell)
    elif args.check :
        argsshellcheck(args.check)
    elif args.remove :
        ShellPath().remove()
    elif args.generate :
        if args.generate[0] == 'shell.yaml' :
            with open(args.generate[0],'r') as f:
                shell = yaml.load(f, Loader=yaml.FullLoader)
            print(shell)
            ShellPath().generateshell(shell=shell)
        else :
            print("Please load shell.yaml file")
            
    elif args.first : 
        path = os.path.abspath(__file__).split(os.sep)[:-3]
        path.append("pdos")
        path = os.sep.join(path)
        
        for f in os.listdir(path) :
            os.chmod("%s/%s"%(path,f),0o0777)
        #os.chmod()