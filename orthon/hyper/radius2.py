from rsf.proj import *

def radius2(high, low,               # initial high-resolution and legacy images
            niter,                   # number of corrections
            c,                       # step length for radius corrections. Can be type int or float for constant c 
                                     #      or type array for changing c.
            bias=-15, clip=30,       # bias and clip for display
            rect1=40, rect2=80,      # radius for local frequency calculation
            maxrad=1000,             # maximum allowed radius
            theor=True,              # use theoretical smoothing radius
            scale=9,                 # scale for theoretical smoothing radius
            initial=10,              # initial value for contant smoothing radius 
            minval=0, maxval=25,     # minval and maxval for freqdif-filt display
            titlehigh="Hires",       # high title
            titlelow="Legacy",       # low title
            it=0):                   # correction number (from newarp.py)

    if type(c) is float or type(c) is int:
        c = [c]*niter

    def seisplot(name):
        return 'grey title="%s" screenratio=1.2'%name
    
    locfreq = '''iphase order=10 rect1=%d rect2=%d hertz=y complex=y | 
                 put label="Frequency" unit=Hz'''%(rect1,rect2)

    def locfreqplot(name):
        return 'grey mean=y color=j scalebar=y title="%s" screenratio=1.2'%name

    freqdif = 'add scale=-1,1 ${SOURCES[1]} | put label=Frequency'

    def freqdifplot(num):
        #return '''grey allpos=y color=j scalebar=y mean=y title="Difference in Local Frequencies %d" ''' %num
        #return '''grey allpos=y color=j scalebar=y mean=y title="Difference in Local Frequencies %d" 
        #          clip=%d bias=%d minval=%d maxval=%d''' %(num,clip,bias,minval,maxval) 
        return '''grey allpos=y color=j scalebar=y mean=y title="Difference in Local Frequencies %d" 
                  clip=%d bias=%d minval=%d maxval=%d screenratio=1.2''' %(num,22,-15,-15,8) 

    specplot = '''cat axis=2 ${SOURCES[1]} | 
                  scale axis=1 | window max1=125 |
                  graph title="Normalized Spectra" label2="Amplitude" unit2=""'''

    def rectplot(name):
        return 'grey color=j mean=y title="%s" scalebar=y barlabel=Radius barunit=samples screenratio=1.2'%name

    smooth = 'nsmooth1 rect=${SOURCES[1]}'

    #Result(high, seisplot(titlehigh))
    #Result(low, seisplot(titlelow))

    Flow('high-freq%i'%it,high,locfreq)
    Result('high-freq%i'%it,locfreqplot('%s Local Frequency'%titlehigh))

    Flow('low-freq%i'%it,low,locfreq)
    Result('low-freq%i'%it,locfreqplot('%s Local Frequency'%titlelow))

    Flow('freqdif%i'%it,'low-freq%i high-freq%i'%(it,it),freqdif)

    Result('freqdif%i'%it,freqdifplot(-1))

    # initial smoothing radius
    if (theor):
        from math import pi
        Flow('rect0%i'%it,'low-freq%i high-freq%i'%(it,it),'''math f1=${SOURCES[1]} 
              output="sqrt(%g*(1/(input*input)-1/(f1*f1)))/%g" '''%(scale,2*pi*0.001))
    else:
        Flow('rect0%i'%it,'low-freq%i'%it,'math output=%f'%initial)

    Result('rect0%i'%it,rectplot("Smoothing Radius 0"))

    Flow('high-smooth0%i'%it,'%s rect0%i'% (high,it),smooth)
    Result('high-smooth0%i'%it, seisplot("%s Smooth 0"%titlehigh))

    Flow('high-spec%i'%it,high,'spectra all=y')
    Flow('low-spec%i'%it,low,'spectra all=y')
    Flow('high-smooth-spec0%i'%it,'high-smooth0%i'%it,'spectra all=y')
    Result('nspectra%i'%it,'high-spec%i low-spec%i'%(it,it),specplot)
    Result('high-smooth-spec0%i'%it,'high-smooth-spec0%i low-spec%i'%(it,it),specplot)
    
    Flow('high-smooth-freq0%i'%it,'high-smooth0%i'%it,locfreq)
    Result('high-smooth-freq0%i'%it,locfreqplot("%s Local Frequency Smoothed %d" %(titlehigh,0)))
           
    Flow('freqdif-filt0%i'%it,'low-freq%i high-smooth-freq0%i'%(it,it),freqdif) 
    Result('freqdif-filt0%i'%it,freqdifplot(0)) 

    #prog=Program('radius2.c') 
    for i in range(1, niter+1): 
        j = i-1
        Flow('rect%d%i rect-low%d%i rect-high%d%i'%(i,it,i,it,i,it),'rect%d%i freqdif-filt%d%i'%(j,it,j,it),
             'radius2 freq=${SOURCES[1]} low=${TARGETS[1]} high=${TARGETS[2]} c=%f'%c[j])
        Result('rect%d%i'%(i,it),rectplot("Smoothing Radius %d")%i)
        Result('rect-low%d%i'%(i,it),rectplot("Smoothing Radius %d low")%i)
        Result('rect-high%d%i'%(i,it),rectplot("Smoothing Radius %d high")%i)
        
        Flow('high-smooth%d%i'%(i,it),'%s rect-high%d%i'%(high,i,it),smooth)
        Result('high-smooth%d%i'%(i,it), seisplot('%s Smooth %d'%(titlehigh,i)))

        Flow('low-smooth%d%i'%(i,it),'%s rect-low%d%i'%(low,i,it),smooth)
        Result('low-smooth%d%i'%(i,it), seisplot('%s Smooth %d'%(titlelow,i)))

        Flow('high-smooth-spec%d%i'%(i,it),'high-smooth%d%i'%(i,it),'spectra all=y')
        Flow('low-smooth-spec%d%i'%(i,it),'low-smooth%d%i'%(i,it),'spectra all=y')
        Result('high-smooth-spec%d%i'%(i,it),'high-smooth-spec%d%i low-smooth-spec%d%i'%(i,it,i,it),specplot)
        
        Flow('high-smooth-freq%d%i'%(i,it),'high-smooth%d%i'%(i,it),locfreq)
        Result('high-smooth-freq%d%i'%(i,it),locfreqplot('%s Local Frequency Smoothed %d'%(titlehigh,i)))

        Flow('low-smooth-freq%d%i'%(i,it),'low-smooth%d%i'%(i,it),locfreq)
        Result('low-smooth-freq%d%i'%(i,it),locfreqplot('%s Local Frequency Smoothed %d'%(titlelow,i)))

        Flow('freqdif-filt%d%i'%(i,it),'low-smooth-freq%d%i high-smooth-freq%d%i'%(i,it,i,it),freqdif)
        Result('freqdif-filt%d%i'%(i,it),freqdifplot(i))
       
