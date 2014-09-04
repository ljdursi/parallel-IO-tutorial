#!/bin/env python
"""
./plots.py [--help] file
Plots the supplied hdf5 file, assuming it consists of a 2d density array
and a 2d array of 2d velocities.
"""

import tables
import netCDF4
import numpy
import matplotlib
import pylab
import sys
import getopt


def readHDF5file(filename):
    densData = None
    velData = None
    ydims = None
    xdims = None

    h5file=tables.openFile(filename,mode="r")
#    densData = h5file.root.ArrayData.dens.read()
#    velData = h5file.root.ArrayData.vel.read()

    for node in h5file.walkNodes('/',classname='Array'):
      if (node.name == "dens"):
        densData = node.read()
      if (node.name == "vel"):
        velData = node.read()

    print 'veldata shape = ', velData.shape
    if ((velData.shape)[2] <= 3):
        velData = numpy.transpose(velData,(2,1,0))
        densData = numpy.transpose(densData,(1,0))
        print 'Transposing...'

    return ("HDF5", densData,velData,ydims,xdims)

def readNetCDF4file(filename):
    densData = None
    velData = None
    ydims = None
    xdims = None
    file=netCDF4.Dataset(filename,"r")
    ycoordname = 'Y coordinate'
    xcoordname = 'X coordinate'
    densname   = 'Density'
    velname    = 'Velocity'

    if ycoordname in file.variables:
        ydims=file.variables['Y coordinate'][:]
    if xcoordname in file.variables:
        xdims=file.variables['X coordinate'][:]
    if densname in file.variables:
        densData = file.variables['Density'][:,:]
        print 'densData shape = ', densData.shape
    if velname in file.variables:
        velData = file.variables['Velocity'][:,:,:]

    if not velData is None:
        print 'veldata shape = ', velData.shape
        if ((velData.shape)[2] <= 3):
            velData = numpy.transpose(velData,(2,1,0))
            densData = numpy.transpose(densData,(1,0))
            print 'Transposing...'

    file.close()
    return ("NetCDF4", densData,velData,xdims,ydims)

def getData(filename="data.h5"):
    densData = None
    velData  = None
    ydims = None
    xdims = None
    ext = ((filename.split("."))[-1]).lower()
    if ext == "h5" or ext == "hdf5" or ext == "hdf" :
        (filetype, densData, velData, xdims, ydims) = readHDF5file(filename)
    elif ext == "nc" or ext == "nc4" or ext == "netcdf"  or ext=="ncdf":
        (filetype, densData, velData, xdims, ydims) = readNetCDF4file(filename)
    return (filetype,densData, velData, xdims, ydims)


def plot2darray(filename="data.h5"):
    (filetype,densData, velData, xdims, ydims) = getData(filename)

    if densData is None and velData is None:
        print "No Data in file "+filename
        return

    print "Plotting ",filename
    if not densData is None:
        densDataT = numpy.transpose(densData)
    if not velData is None:
        vx = numpy.transpose(velData[0,:,:])
        vy = numpy.transpose(velData[1,:,:])

    nx = None
    ny = None
    if not densData is None:
        size = densData.shape
        nx = size[0]
        ny = size[1]
    elif not velData is None:
        size = vx.shape
        nx = size[1]
        ny = size[0]
    
    if xdims is None:
        hdims = numpy.arange(0,nx)
    else:
        hdims = xdims

    if ydims is None:
        vdims = numpy.arange(0,ny)
    else:
        vdims = ydims 

    X,Y = pylab.meshgrid(hdims,vdims)

    # want about 20 arrows per dim
    narrows = 20 
    everyX = nx/narrows
    everyY = nx/narrows

    fig = pylab.figure(1)
    ax = fig.add_subplot(111,aspect='equal')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title(filetype+' 2D Arrays Output: '+filename)
    ax.set_xlim(min(hdims), max(hdims))
    ax.set_ylim(min(vdims), max(vdims))

    if not densData is None:
        ax.contourf(hdims,vdims,densDataT)

    if not velData is None:
        maxv = numpy.max(numpy.sqrt(vx*vx + vy*vy))
        ax.quiver(X[::everyX,::everyY], Y[::everyX, ::everyY], vx[::everyX,::everyY],vy[::everyX,::everyY], scale=4.*maxv)

    pylab.show()


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "h", ["help"])
        except getopt.error, msg:
             raise Usage(__doc__)

        for o,a in opts:
            if o in ("-h", "--help"):
                print __doc__
                sys.exit(0)
        
        if args == []:
            plot2darray()
        for arg in args:
            plot2darray(arg)            
        sys.exit(0)

    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, "for help use --help"
        return 2


if __name__ == "__main__":
    sys.exit(main())

