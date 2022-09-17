#! / usr / bin / env python

#------------------------------------------------------------------------------

#$Id: TestLSQs.py 9880 2013-06-04 16:15:16Z NicholasCross $

"""

   TestLSQ: test the least squares fitting.


   @author: N.J.G. Cross

   @org:    WFAU, IfA, University of Edinburgh



"""

#------------------------------------------------------------------------------

from __future__      import division, print_function

import matplotlib.pyplot    as plt

import numpy

import math

from scipy import optimize




#------------------------------------------------------------------------------


class BlackBodyTester(object):

    """ Tests current and future corrections

    """



    def run(self):

        """

        Take values of Vega - 0 mag in all bands



        """

        VegaToAB = numpy.array([0.1,0.2,0.93,1.37,1.90])



        lambArray = numpy.array([4.5e-7, 8.0e-7, 1.2e-6,1.6e-6,2.2e-6])

        magArray=numpy.array([0.,0.,0.,0.,0.])

        magErrArray=numpy.array([0.2,0.1,0.1,0.15,0.2])

        fluxArray =  numpy.power(10.,-0.4*(magArray+VegaToAB-8.90))

        fluxErrArray = fluxArray*(2.5/math.log(10.))*magErrArray

        data = numpy.zeros([lambArray.size,3])

        data[:,0]=lambArray

        data[:,1]=fluxArray

        data[:,2]=fluxErrArray



        print("DATA: ",data)



        parameters = [1.e-27,10000.]



        fit,chiSquared = leastSquaresFit(_blackbodyModel, parameters, data)

        print("FIT: ",fit,chiSquared)

        model=[]

        for ii in range(lambArray.size):

            expVal=_blackbodyModel(lambArray[ii],fit[0],fit[1])

            print("MODEL: ",lambArray[ii],expVal)

            model.append(expVal)



        fig = plt.figure()



        plt.errorbar(1.e6*lambArray, fluxArray, yerr=fluxErrArray, fmt='o',label='both limits (default)')

        plt.plot(1.e6*lambArray,numpy.array(model),'-')

        plt.xlabel("Wavelength / micrometres")

        plt.show()



    #------------------------------------------------------------------------------



def leastSquaresFit(function, parameters, data):

    """ Non-linear least-squares fitting function using Levenburg-Marquardt

    """

    xdata = numpy.array([dp[0] for dp in data])

    ydata = numpy.array([dp[1] for dp in data])

    if len(data) > 0 and len(data[0]) > 2:

        sigma = numpy.array([dp[2] for dp in data])

    else:

        sigma = None

    # @NOTE optimize.curve_fit() has parameter absolute_sigma, but this is fed into

    #       optimize.leastsq() which does not

    fit, _covMat, infoDict, _comment, _rank = optimize.curve_fit(function,

        xdata, ydata, p0=parameters, sigma=sigma, full_output=True) # absolute_sigma=True)

    # Calc chiSq

    chiSquared = numpy.sum(infoDict['fvec'] ** 2)


    return (fit, chiSquared)

#------------------------------------------------------------------------------



def _blackbodyModel(wavelength, norm_coeff, temp_coeff):

    """ Returns a blackbody model

    """

    h=6.626e-34

    c=2.998e8

    k=1.380e-23

    expConst=h*c/(wavelength*k)

    return norm_coeff*numpy.power(wavelength,-5)*1./(numpy.exp((expConst/temp_coeff)-1.))




#------------------------------------------------------------------------------


if __name__ == '__main__':

        # Define specific command-line interface settings required by CU13



    cu = BlackBodyTester()

    cu.run()
