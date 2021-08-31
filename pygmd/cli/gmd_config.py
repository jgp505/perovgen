import os
import sys

from perovgen.pygmd.shell import ShellPath

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

