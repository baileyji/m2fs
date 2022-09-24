#!/usr/bin/env python
from Tkinter import *
from m2fs.focus import RFOCUS, BFOCUS

def rfocus(temp):
    try:
        temp=float(temp)
        return '%.1f'%RFOCUS(temp)
    except Exception:
        return 'Enter Temp'

def bfocus(temp):
    try:
        temp=float(temp)
        return '%.1f'%BFOCUS(temp)
    except Exception:
        return 'Enter Temp'

master = Tk()

Label(master, text="R").grid(row=1, column=0)
Label(master, text="B").grid(row=2, column=0)

Label(master, text="Cradle Temp").grid(row=0, column=1)
Label(master, text="Focus Value").grid(row=0, column=2)

rfocusv = StringVar()
bfocusv = StringVar()
Label(master, textvariable=rfocusv).grid(row=1,column=2)
Label(master, textvariable=bfocusv).grid(row=2,column=2)

rtemp = StringVar()
rtemp.trace('w',
            lambda name, index, mode, sv=rtemp: rfocusv.set(rfocus(sv.get())))
Entry(master, textvariable=rtemp, width=4).grid(row=1, column=1)

btemp = StringVar()
btemp.trace('w',
            lambda name, index, mode, sv=btemp: bfocusv.set(bfocus(sv.get())))
Entry(master, textvariable=btemp, width=4).grid(row=2, column=1)


rfocusv.set('Enter Temp')
bfocusv.set('Enter Temp')

master.wm_title("M2FS Focus")

mainloop()
