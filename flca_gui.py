import os
import tkinter as tk
from tkinter import filedialog as fd
import flca
import time

indir0='/home/usr/data/reference'
indir='/home/usr/data/origin'
outdir='/home/usr/data/aligned'
firstfits='/home/usr/data/first_fits_file.fits'
x0='300'; x1='700'; y0='100'; y1='400'
time0a='-23'; time0b='-8'; time1a=time0a; time1b=time0b
every='1'; skip='10'; xoffset='5'; yoffset='5'
sigma='20'; threshold='0'; kr='0.4'; biascorrect=1; interpolate=0;
if os.path.exists('./input.txt'):
    f=open('./input.txt', 'r')
    lines=f.readlines()
    f.close()
    tmp = []
    for i in range(len(lines)):
        if lines[i][0] == '#': continue
        tmp1 = lines[i].split()
        tmp.extend(tmp1)
    if tmp.count('*indir0'):
        ind = tmp.index('*indir0')
        if tmp[ind+1][0] != '*': indir0=tmp[ind+1]
    if tmp.count('*indir'):
        ind = tmp.index('*indir')
        if tmp[ind+1][0] != '*': indir=tmp[ind+1]
    if tmp.count('*outdir'):
        ind = tmp.index('*outdir')
        if tmp[ind+1][0] != '*': outdir=tmp[ind+1]
    if tmp.count('*firstfits'):
        ind = tmp.index('*firstfits')
        if tmp[ind+1][0] != '*': firstfits=tmp[ind+1]
    if tmp.count('*x0'):
        ind = tmp.index('*x0')
        if tmp[ind+1][0] != '*': x0=int(tmp[ind+1])
    if tmp.count('*x1'):
        ind = tmp.index('*x1')
        if tmp[ind+1][0] != '*': x1=int(tmp[ind+1])
    if tmp.count('*y0'):
        ind = tmp.index('*y0')
        if tmp[ind+1][0] != '*': y0=int(tmp[ind+1])
    if tmp.count('*y1'):
        ind = tmp.index('*y1')
        if tmp[ind+1][0] != '*': y1=int(tmp[ind+1])
    if tmp.count('*time0a'):
        ind = tmp.index('*time0a')
        if tmp[ind+1][0] != '*': time0a=int(tmp[ind+1])
    if tmp.count('*time0b'):
        ind = tmp.index('*time0b')
        if tmp[ind+1][0] != '*': time0b=int(tmp[ind+1])
    if tmp.count('*time1a'):
        ind = tmp.index('*time1a')
        if tmp[ind+1][0] != '*': time1a=int(tmp[ind+1])
    if tmp.count('*time1b'):
        ind = tmp.index('*time1b')
        if tmp[ind+1][0] != '*': time1b=int(tmp[ind+1])
    if tmp.count('*every'):
        ind = tmp.index('*every')
        if tmp[ind+1][0] != '*': every=int(tmp[ind+1])
    if tmp.count('*skip'):
        ind = tmp.index('*skip')
        if tmp[ind+1][0] != '*': skip=int(tmp[ind+1])
    if tmp.count('*xoffset'):
        ind = tmp.index('*xoffset')
        if tmp[ind+1][0] != '*': xoffset=int(tmp[ind+1])
    if tmp.count('*yoffset'):
        ind = tmp.index('*yoffset')
        if tmp[ind+1][0] != '*': yoffset=int(tmp[ind+1])
    if tmp.count('*sigma'):
        ind = tmp.index('*sigma')
        if tmp[ind+1][0] != '*': sigma=int(tmp[ind+1])
    if tmp.count('*threshold'):
        ind = tmp.index('*threshold')
        if tmp[ind+1][0] != '*': threshold=tmp[ind+1]
    if tmp.count('*kr'):
        ind = tmp.index('*kr')
        if tmp[ind+1][0] != '*':
            kr=float(tmp[ind+1])
    if tmp.count('*biascorrect'): biascorrect=1
    if tmp.count('*interpolate'): interpolate=1
    
def browse_indir0():
    tmp = fd.askdirectory(title='Select reference folder')
    if not tmp:
        return
    entin0.delete(0, tk.END)
    entin0.insert(0, tmp)
def browse_indir():
    tmp = fd.askdirectory(title='Select input folder')
    if not tmp:
        return
    entin.delete(0, tk.END)
    entin.insert(0, tmp)
def browse_outdir():
    tmp = fd.askdirectory(title='Select output folder')
    if not tmp:
        return
    entout.delete(0, tk.END)
    entout.insert(0, tmp)
def browse_firstf():
    tmp = fd.askopenfilename(title='Select the first fits',
                             filetypes=[("Fits Files", "*.f?ts"), ("All Files", "*.*")])
    if not tmp:
        return
    entfit.delete(0, tk.END)
    entfit.insert(0, tmp)
def run_onec():
    starttime = time.time()
    indir = entin.get()
    outdir = entout.get()
    if indir == '' or outdir == '':
        print("indir and outdir must be defined.")
        return
    if os.path.exists(outdir) == False:
        os.mkdir(outdir)
        print(f"{outdir} is created")
    firstfits = entfit.get()
    x0 = int(entx0.get())
    x1 = int(entx1.get())
    y0 = int(enty0.get())
    y1 = int(enty1.get())
    every = int(entevery.get())
    skip = int(entskip.get())
    xoffset = int(entxoffset.get())
    yoffset = int(entyoffset.get())
    sigma = int(entsigma.get())
    threshold = entthreshold.get()
    kr = float(entkr.get())
    biascorrect = biastk.get()
    interpolate = inptk.get()
    error = flca.onechannel(indir, outdir, firstfits, x0, x1, y0, y1, every, skip, xoffset, 
               yoffset, sigma, threshold, kr, biascorrect, interpolate)
    if error == 0:
        endtime = time.time()
        print('Finish')
        print(f"The time cost is {round(endtime-starttime)} s.")
    else:
        print("Error!")
def run_twoc():
    starttime = time.time()
    indir0 = entin0.get()
    indir = entin.get()
    outdir = entout.get()
    if indir0 == '' or indir == '' or outdir == '':
        print('indir0, indir and outdir must be defined if twochannel is turned on.')
        return
    if os.path.exists(outdir) == False:
        os.mkdir(outdir)
        print(f"{outdir} is created")
    x0 = int(entx0.get())
    x1 = int(entx1.get())
    y0 = int(enty0.get())
    y1 = int(enty1.get())
    time0a = entt0a.get()
    time0b = entt0b.get()
    if time0a == '' or time0b == '':
        print('time0a and time0b must be defined.')
        return
    else:
        time0a = int(time0a)
        time0b = int(time0b)
    time1a = entt1a.get()
    time1b = entt1b.get()
    if time1a == '':
        if time1b == '':
            time1a = time0a
            time1b = time0b
        else:
            time1b = int(time1b)
            time1a = time1b-(time0b-time0a)
    else:
        time1a = int(time1a)
        if time1b == '':
            time1b = time1a+time0b-time0a
        else:
            time1b = int(time1b)
    skip = int(entskip.get())
    xoffset = int(entxoffset.get())
    yoffset = int(entyoffset.get())
    sigma = int(entsigma.get())
    threshold = entthreshold.get()
    kr = float(entkr.get())
    biascorrect = biastk.get()
    interpolate = inptk.get()
    error = flca.doublechannel(indir0, indir, outdir, x0, x1, y0, y1, time0a, time0b, time1a, time1b, skip, xoffset, 
               yoffset, sigma, threshold, kr, biascorrect, interpolate)
    if error == 0:
        endtime = time.time()
        print('Finish')
        print(f"The time cost is {round(endtime-starttime)} s.")
    else:
        print("Error!")

win = tk.Tk(className="flca_gui")
win.geometry("500x500")
frm0=tk.Frame(master=win)
frm0.grid(row=0, column=0, padx=20)

cvsid = tk.Canvas(frm0, height=60, width=400)
imgid = tk.PhotoImage(file='id.png')
img = cvsid.create_image(220,30,anchor='center', image=imgid)
cvsid.grid(row=0, column=0, sticky='sew')

frm00 = tk.Frame(master=frm0, relief=tk.SUNKEN, borderwidth=3)
frm00.grid(row=1, column=0, pady=6)
lblin = tk.Label(master=frm00, text='Input folder')
btnin = tk.Button(master=frm00, text="Browse...", command=browse_indir)
entin = tk.Entry(master=frm00, width=40)
entin.insert(0, indir)
lblin.grid(row=0, column=0, sticky='e', pady=3)
entin.grid(row=0, column=1, pady=3)
btnin.grid(row=0, column=2, pady=3, padx=3)
lblout = tk.Label(master=frm00, text='Output folder')
btnout = tk.Button(master=frm00, text="Browse...", command=browse_outdir)
entout = tk.Entry(master=frm00, width=40)
entout.insert(0, outdir)
lblout.grid(row=1, column=0, sticky='e', pady=3)
entout.grid(row=1, column=1, pady=3)
btnout.grid(row=1, column=2, pady=3, padx=3)
lblfit = tk.Label(master=frm00, text='First fits')
btnfit = tk.Button(master=frm00, text="Browse...", command=browse_firstf)
entfit = tk.Entry(master=frm00, width=40)
entfit.insert(0, firstfits)
lblfit.grid(row=2, column=0, sticky='e', pady=3)
entfit.grid(row=2, column=1, pady=3)
btnfit.grid(row=2, column=2, pady=3, padx=3)

frm01 = tk.Frame(master=frm0, borderwidth=3)
frm01.grid(row=2, column=0, pady=6)
frm010 = tk.Frame(frm01)
frm010.grid(row=0, column=0, padx=7)
lblx0 = tk.Label(frm010, text='x0')
entx0 = tk.Entry(frm010, width=4)
lblx0.grid(row=0,column=0)
entx0.grid(row=0,column=1)
entx0.insert(0, x0)
frm011 = tk.Frame(frm01)
frm011.grid(row=0, column=1, padx=7)
lblx1 = tk.Label(frm011, text='x1')
entx1 = tk.Entry(frm011, width=4)
entx1.insert(0, x1)
lblx1.grid(row=0,column=0)
entx1.grid(row=0,column=1)
frm012 = tk.Frame(frm01)
frm012.grid(row=0, column=2, padx=7)
lbly0 = tk.Label(frm012, text='y0')
enty0 = tk.Entry(frm012, width=4)
enty0.insert(0, y0)
lbly0.grid(row=0,column=0)
enty0.grid(row=0,column=1)
frm013 = tk.Frame(frm01)
frm013.grid(row=0, column=3, padx=7)
lbly1 = tk.Label(frm013, text='y1')
enty1 = tk.Entry(frm013, width=4)
enty1.insert(0, y1)
lbly1.grid(row=0,column=0)
enty1.grid(row=0,column=1)
#btncheck = tk.Button(master=frm01, text="Check")
#btncheck.grid(row=0, column=4, sticky='e', padx=7)

frm02 = tk.Frame(master=frm0, borderwidth=3)
frm02.grid(row=3, column=0, pady=6)
frm020 = tk.Frame(frm02)
frm020.grid(row=0, column=0, padx=7)
lblevery = tk.Label(frm020, text='every')
entevery = tk.Entry(frm020, width=4)
lblevery.grid(row=0,column=0)
entevery.grid(row=0,column=1)
entevery.insert(0, every)
frm021 = tk.Frame(frm02)
frm021.grid(row=0, column=1, padx=7)
lblskip = tk.Label(frm021, text='skip')
entskip = tk.Entry(frm021, width=4)
entskip.insert(0, skip)
lblskip.grid(row=0,column=0)
entskip.grid(row=0,column=1)
frm022 = tk.Frame(frm02)
frm022.grid(row=0, column=2, padx=7)
lblxoffset = tk.Label(frm022, text='xoffset')
entxoffset = tk.Entry(frm022, width=4)
entxoffset.insert(0, xoffset)
lblxoffset.grid(row=0,column=0)
entxoffset.grid(row=0,column=1)
frm023 = tk.Frame(frm02)
frm023.grid(row=0, column=3, padx=7)
lblyoffset = tk.Label(frm023, text='yoffset')
entyoffset = tk.Entry(frm023, width=4)
entyoffset.insert(0, yoffset)
lblyoffset.grid(row=0,column=0)
entyoffset.grid(row=0,column=1)

#quiettk = tk.IntVar(value=quiet); 
biastk = tk.IntVar(value=biascorrect)
inptk = tk.IntVar(value=interpolate)
frm03 = tk.Frame(master=frm0, borderwidth=3)
frm03.grid(row=4, column=0, pady=6)
frm030 = tk.Frame(frm03)
frm030.grid(row=0, column=0, padx=7)
lblsigma = tk.Label(frm030, text='sigma')
entsigma = tk.Entry(frm030, width=4)
lblsigma.grid(row=0,column=0)
entsigma.grid(row=0,column=1)
entsigma.insert(0, sigma)
frm031 = tk.Frame(frm03)
frm031.grid(row=0, column=1, padx=7)
lblthreshold = tk.Label(frm031, text='threshold')
entthreshold = tk.Entry(frm031, width=7)
entthreshold.insert(0, threshold)
lblthreshold.grid(row=0,column=0)
entthreshold.grid(row=0,column=1)
frm032 = tk.Frame(frm03)
frm032.grid(row=0, column=2, padx=7)
lblkr = tk.Label(frm032, text='kr')
entkr = tk.Entry(frm032, width=4)
entkr.insert(0, kr)
lblkr.grid(row=0,column=0)
entkr.grid(row=0,column=1)
cbtinp = tk.Checkbutton(frm03, text='interpolate', variable=inptk, onvalue=1, offvalue=0)
cbtinp.grid(row=0,column=4, padx=7)
cbtbias = tk.Checkbutton(frm03, text='biascorrect', variable=biastk, onvalue=1, offvalue=0)
cbtbias.grid(row=0,column=3, padx=7)

'''frm04 = tk.Frame(frm0)
frm04.grid(row=5, column=0, ipadx=5, ipady=5, padx=10, pady=6)
btncreate = tk.Button(frm04, text='Create input.txt')
btncreate.grid(row=0, column=0, ipadx=10)
btnrun = tk.Button(frm04, text='Create and Run')
btnrun.grid(row=0, column=1, padx=10, ipadx=10)'''
btnrun = tk.Button(frm0, text='Single-channel coalignment', command=run_onec)
btnrun.grid(row=5, column=0, ipadx=10)

frm05 = tk.Frame(master=frm0, relief=tk.SUNKEN, borderwidth=3)
frm05.grid(row=6, column=0, pady=6)
lblin0 = tk.Label(master=frm05, text='Reference folder')
btnin0 = tk.Button(master=frm05, text="Browse...", command=browse_indir0)
entin0 = tk.Entry(master=frm05, width=40)
entin0.insert(0, indir0)
lblin0.grid(row=0, column=0, sticky='e', pady=3)
entin0.grid(row=0, column=1, pady=3)
btnin0.grid(row=0, column=2, pady=3, padx=3)

frm06 = tk.Frame(master=frm0, borderwidth=3)
frm06.grid(row=7, column=0, pady=6)
frm060 = tk.Frame(frm06)
frm060.grid(row=0, column=0, padx=7)
lblt0a = tk.Label(frm060, text='time0a')
entt0a = tk.Entry(frm060, width=4)
lblt0a.grid(row=0,column=0)
entt0a.grid(row=0,column=1)
entt0a.insert(0, time0a)
frm061 = tk.Frame(frm06)
frm061.grid(row=0, column=1, padx=7)
lblt0b = tk.Label(frm061, text='time0b')
entt0b = tk.Entry(frm061, width=4)
entt0b.insert(0, time0b)
lblt0b.grid(row=0,column=0)
entt0b.grid(row=0,column=1)
frm062 = tk.Frame(frm06)
frm062.grid(row=0, column=2, padx=7)
lblt1a = tk.Label(frm062, text='time1a')
entt1a = tk.Entry(frm062, width=4)
entt1a.insert(0, time1a)
lblt1a.grid(row=0,column=0)
entt1a.grid(row=0,column=1)
frm063 = tk.Frame(frm06)
frm063.grid(row=0, column=3, padx=7)
lblt1b = tk.Label(frm063, text='time1b')
entt1b = tk.Entry(frm063, width=4)
entt1b.insert(0, time1b)
lblt1b.grid(row=0,column=0)
entt1b.grid(row=0,column=1)

btnrun2 = tk.Button(frm0, text='Double-channel coalignment', command=run_twoc)
btnrun2.grid(row=8, column=0, ipadx=10)

#frm1=tk.Frame(master=win, width=400, height=400, relief=tk.SUNKEN, borderwidth=3)
#frm1.grid(row=0, column=1)


win.mainloop()
