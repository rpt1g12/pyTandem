from scipy import interpolate

class interf(object):
    """Returns an interpolating function"""
    def __init__(self, x,y):
        super(interf, self).__init__()
        xspln=interpolate.splrep(x,y,s=0)
        self.xspln = xspln
    def f(self,x):
        """Returns interpolated function value"""
        y=interpolate.splev(x,self.xspln,der=0)
        return y
        
