import numpy as np

class vdf:
    ''' 
    vdf object class
    -----------------
    this is a base class for looking up velocity 
    distribution function values. 

    METHOD vdf.evaluate(arg)
    the vdf.evaluate(arg) function returns evaluations 
    of the vdf, where are is assumed to be a 3-vector or 
    array of 3-vectors 

    METHOD vdf.set_func(arg)
    assigns the evaluation method to the function of 
    choice, where are is the [string] name of the function

    METHOD vdf.set_udata(arg)
    stores any support data that is necessary for 
    calling the evaluation function

    One sample implimentation is included: 
    'maxwellian_vdf' -- evaluates a 3D maxwellian, assuming 
                        that the udata supplies a 
                        structure with tags {vx, vy, vz, w, n}
    '''
    def __init__(self, func_name, kwargs):
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])
        
        if(func_name == 'maxwellian_vdf'):
            self.set_function = self.maxwellian_vdf

    def maxwellian_vdf(self, vxyz):
        vecdim = vxyz.shape

        if(vecdim[0] == 3):
            vx, vy, vz = vxyz
        else:
            vx, vy, vz = vxyz.T

        g = (self.n / (np.sqrt(np.pi)*self.w)**3) *\
            np.exp(-((vx-self.ux)**2 + (vy-self.uy)**2 + (vz-self.uz)**2 )/(self.w**2))
        return g
    
    def evaluate(self, vxyz):
        return self.set_function(vxyz)


if __name__=='__main__':
    # would be more orderly to extend vdf with maxwellian_vdf and then
    # have these parameters as object vars instead of in the pointer
    # because maxwellian_vdf will break now if you don't give it the
    # right parameter format. Regardless...
    g = vdf('maxwellian_vdf', {'ux':100., 'uy':0., 'uz':0., 'w':50., 'n':100.})
    print(g.evaluate(np.array([[80, 0, 0], [90, 0, 0], [100, 0, 0], [110, 0, 0]])))
