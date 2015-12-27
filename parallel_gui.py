#!/usr/bin/python
from Tkinter import *
import os
import subprocess
import tkFileDialog

seq1=''
seq2=''
ans=''

root=top = Tk()

from tkFileDialog import askopenfilename
from tkFileDialog import *


def openfile1():
   global seq1,seq2;
   filename = askopenfilename(filetypes=[("Text files","*.txt")])
   f = open(filename)

   seq1=f.read()
   print seq1

def openfile2():
   global seq1,seq2;
   filename = askopenfilename(filetypes=[("Text files","*.txt")])
   f = open(filename)
   seq2=f.read()
   print seq2

def export():
    global ans
    """name=asksaveasfilename(mode='w',defaultextension=".txt")
    text2save=str(text.get(0.0,END))
    name.write(ans)
    name.close"""


    #def file_save():
    f = tkFileDialog.asksaveasfile(mode='w', defaultextension=".txt")
    if f is None: # asksaveasfile return `None` if dialog closed with "cancel".
        return
    text2save = str(text.get(1.0, END)) # starts from `1.0`, not `0.0`
    f.write(text2save)
    f.close() # `()` was missing.

menubar = Menu(top)
filemenu = Menu(menubar, tearoff=0)
#filemenu.add_command(label="Import from file", command=openfile)

submenu = Menu(filemenu, tearoff=0)
submenu.add_command(label="Sequence 1", command=openfile1)
submenu.add_command(label="Sequence 2", command=openfile2)


filemenu.add_cascade(label='Import', menu=submenu, underline=0)
filemenu.add_separator()


filemenu.add_command(label="Export", command=export)
filemenu.add_command(label="Exit", command=root.quit)

menubar.add_cascade(label="File", menu=filemenu)

root.config(menu=menubar)


top.title('Parallel DNA Sequence Alignment')

frame = Frame(top)
frame.pack(side=TOP)

bottomframe = Frame(top)
bottomframe.pack( side = BOTTOM )


def set1():
    global seq1
    seq1=E1.get()
def set2():
    global seq2
    seq2=E2.get()

L1 = Label(top, text="Sequence 1")
L1.pack( side = LEFT)
E1 = Entry(top, bd = 5)

E1.pack(side = LEFT)
bset1 = Button(text="Set 1", width=10, command=set1)
bset1.pack(side=LEFT)


L2 = Label(top, text="Sequence 2")
L2.pack( side = LEFT)
E2 = Entry(top, bd = 5)

E2.pack(side = LEFT)
bset2 = Button(text="Set 2", width=10, command=set2)
bset2.pack(side=LEFT)

L3 = Label(top, text="Output")
L3.pack();
text = Text(bottomframe,height=25,width=130)
text.insert(INSERT, "The output space")

text.pack()


def callback():
    global seq1,seq2
    #seq1=E1.get()
    #seq2=E2.get()
    ans=seq1+seq2+'\n'
    #print seq1+ seq2
    text.delete("1.0",END)
    #z=align(seq1,seq2)
    command="g++ -o dna_run DNA.cpp"


    f = open('input.txt','w')
    f.write(seq1+'\n'+seq2+'\n')
    f.close()



    runner="./dna_run < input.txt"

    p = subprocess.Popen(command, shell=True, bufsize=0, stdout=subprocess.PIPE, universal_newlines=True)
    p.wait()
    p = subprocess.Popen(runner, shell=True, bufsize=0, stdout=subprocess.PIPE, universal_newlines=True)
    for line in iter(p  .stdout.readline, ''):
            ans+=line+'\n'
    p.wait()
    text.insert(END,ans)
    #print ans
    #ans=z[1]+'\n'+z[2]+'\n'
    #text.insert(END,ans)


b = Button(bottomframe, text="Align", width=10, command=callback)
b.pack()


top.mainloop()

