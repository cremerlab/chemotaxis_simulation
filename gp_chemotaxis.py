###################################################################################
#Chemotaxis simulation
###################################################################################
#Code to simulate chemotaxis driven expansion and population growth,
#following the Growth-Expansion Model.
#Code can also used to study competition and selection within the 
#expanding documentions.
#The Growth-Expansion model and its biological context is introduced
#our manuscript:
#J.Cremer, T.Honda, Y.Tang, J. Wong-Ng, M.Vergassola, T.Hwa 
#"Chemotaxis as a navigation strategy to thrive in nutrient-replete environments"
#
#Code can also be used to study competition and selection dynamics, as
#in our manuscript:
#W.Liu, J.Cremer, D.Li, T.Hwa, C.Liu
#"An evolutionary stable strategy to colonize  spatially extended habitats"
#
#August 2019, Jonas Cremer and all coauthors.
#
###################################################################################
#Simulations using Python 2.7 and the partial differential equation solver FiPy
#Code available via GitHub at. See provided README file for additional information.
###################################################################################

import datetime
import time
import os
import numpy as np
from fipy import *
import numpy as np
import json


#general scalingfacotrs
rhofactor=0.5 #how to convert rho to  OD, rho=0.5 for our strain.
nfactor=1.
mfactor=1.         
tfactor=1./3600.

parsimsweep1={}
parsimsweep2={}
#load JSON parameter file
def load_parameters_simulation(pathin):
    print "look for parameter files:"
    print os.path.join(os.getcwd(),pathin)
    print os.path.join(os.getcwd(),os.pardir)
    pardir=os.path.abspath(os.path.join(os.getcwd(),os.pardir))
    print os.path.join(os.getcwd(),pathin)
    try:        
        json_data=open(os.path.join(os.getcwd(),pathin))
    except:
        json_data=open(os.path.join(pardir,pathin))
    try:
        dictpar = json.load(json_data)
    except:
        dictpar = json.load(json_data)
        print "problem with json format; check if file exist: "+str(os.path.join(pardir,pathin))
        dictpar={}
    json_data.close() 
    return dictpar
    
#run simulations, for one set of parameters
def run_simulation(p,loaddata=False,datastart=[]):
    print "Start Simulation with parameters: "
    print json.dumps(p)

    #####
    #setup simultion according to simpulation parameters (argument p)
    ####    
    
    #setup mesh and numerical parameters. Depending on dimension, mesh is 1D, 2D, or 1Dcylindrical.
    if p['nx']==-2:
        nxc=1000 #use for sizes up to 40mm
        
    else:
        nxc=p['nx']
        
    dx=(p['size'])/float(nxc)
    if p['nx']==-3:
        dx=0.08/2000.
        nxc=int(p['size']/dx)
    
    if p['dimensionmode']=="two":   
        mesh = Grid2D(dx= dx, dy=dx, nx=nxc, ny=nxc)
        xsimused=np.linspace(0,p['size'],nxc)
    elif p['dimensionmode']=="radial":  
        mesh = CylindricalGrid1D(dr=dx, nr=nxc, Lr=p['size'], origin=(p['radius_cut'], ))# add communicator if needed, ,http://www.ctcms.nist.gov/fipy/fipy/generated/fipy.meshes.html
        xsimused=np.linspace(0,p['size'],nxc)
    else:
        mesh = Grid1D(dx= dx, nx=nxc)
        xsimused=np.linspace(0,p['size'],nxc)
    

    
    runmode=p['runmode'] #weber, noweber,...
    dynamicmode=p['dynamicsmode'] #sets the dynamics of consumption etc....
   
    try:
        widthprofile_nutrients=p["widthprofile_nutrients"]
        initialnutrientprofile=p["initialnutrientprofile"]
    except:
        initialnutrientprofile=0
   
    #give names to important variables
    nx=nxc
    dt = p['dt']
    D_nutrients=p['D_nutrients'] #assume same diffusion constant for all nutrients....
    timetotal=p['time'] #runtime [s]
    vsteps=p['viewsteps'] #number of outputs
    
    #get parameters for nutrients
    try:
        startnprofile=p["startnprofile"]
    except:
        startnprofile="constant"
    
    #get parameters for different strains
    lambdac=[]
    lambdadeltac=[]
    rhostart=[]
    Km1=[]
    Y1=[]
    Km2=[]
    Kml2=[]
    Y2=[]
    Km3=[]
    Y3=[]
    xhi=[]
    xhiconstant=[]
    D_chemotaxis=[]
    weber_offset=[]
    weber_offsetmax=[]
    Y2k1=[]
    Y2k2=[]
    
    CClog=[]
    CClog1=[]
    deltal1=[]
    deltal2=[]
    hilllog=[]
    hilllogCT=[]
    xhiTWO=[]
    Y3k1=[]
    Y3k2=[]
        
    try:
        weber_offsetlogalpha=p["weber_offset_logalpha"]
    except:
        weber_offsetlogalpha=1.

    
    usefacegrad=True #always use facegrad for radial coordinates.
    try:
        chemgrad_mode=p["chemgrad_mode"]
    except:
        chemgrad_mode="gradient" #gradient is standard
    
    try:
        lambda_frontd=p["lambda_front_delta"]
    except:
        lambda_frontd=0.
    try:
        lambda_back=p["lambda_back"]
    except:
        try:
            lambda_back=p["lambda_front"]-lambda_frontd
        except:
            lambda_back=0.
    try:
        lambdaodchange=p["lambda_odchange"]
    except:
        lambdaodchange=0.
          
    try: 
        competitionmode=p["competition_mode"]
    except:
        competitionmode=False
        
    try: 
        epsiloncompetitionmodestr=p["competition_mode_epsilon"].split(",")
        epsiloncompetitionmode=[]
        for istri in range(0,len(epsiloncompetitionmodestr)):
            epsiloncompetitionmode.append(float(epsiloncompetitionmodestr[istri]))
    except:
        epsiloncompetitionmode=[]
    
    try:
        evolutionmode=p["evolution_mode"]
    except:
        evolutionmode=False
    if evolutionmode:
        evolutionstd=p["evolution_std"]
        evolutionvar=p["evolution_var"]
        evolutionvalue=p["evolution_value"]
        evolutionselectiondimension=p["evolution_selectiondimension"]
        evolutionposition=p["evolution_position"]
        
        if evolutionposition==-1: #position is size of simulated plate
            evolutionposition=0.9999*p['size']
        evolutiondensity=p["evolution_density"]
        try:
            evolutionstep=p["evolution_step"]
        except:
            evolutionstep=0
        try:
            evolutionsinglemutation=p["evolution_singlemutation"]
        except:
            evolutionsinglemutation=-1
    if evolutionmode==False:
        n_strains=int(p['number_strains'])
    else:
        n_strains=int(p['evolution_strains'])
        
    try:
        tworingmode=False
        tworingmodec=p['tworingmode']
        if tworingmodec==1:
            tworingmode=True
    except:
        tworingmode=False
    if tworingmode:
        n_nutrients=3
    else:
        n_nutrients=int(p['number_nutrients'])
    try:
        scalexhiDwithgrowth=p['scalewithgrowthxhiD']
    except:
        scalexhiDwithgrowth=-1
    try:
        scalemodexhiD=p['scalemodexhiD']*1
    except:
        scalemodexhiD=-1
    try:
        scalemodetwo=p['scalemodetwo']*1 #if 2nd strain is a bit slower than first #only used if othersclemodes not used
    except:
        scalemodetwo=-1
    try:
        scalemodeDxhi=p['scalemodeDxhi']*1
    except:
        scalemodeDxhi=-1
    print scalemodexhiD
    try:
        fixed_totalrho=p['rhototfixed']
    except:
        fixed_totalrho=-1.
    try:
        same_rho=p['same_rho']
    except:
        same_rho=False
        
    try:
        odcuttofflogistic=p['od_cutoff_logistic']
    except:
        odcuttofflogistic=1.
        
        
        
    try:
        change_xhiB=p['change_xhiB'] #if larger than zero, that is the time in seconds where dynamics is flippted to 2nd strain having different xhi (only works for n_strains=2)
    except:
        change_xhiB=-1
    try:
        change_mode=p['change_mode']
    except:
        change_mode=0

   
    #go through different strains and set parameters
    for isc in range(1,n_strains+1):
        
        if evolutionmode==True or competitionmode>0:
            iscuse=1
        else:
            iscuse=isc
        #growthmode: same. sets starting density and standard growth condition to same value as first strain
        if p['growthmode']=="same":
            iscuse=1
        
        if p['growthmode']=="samegrowth":
        
            lambdac.append(p['lambda_'+str(1)])
        else:
            
                lambdac.append(p['lambda_'+str(iscuse)])
           
        
        
        try:
            lambdadeltac.append(p['lambdadelta_'+str(iscuse)])
        except:
            lambdadeltac.append(1.2)
        if same_rho==True:
            rhostart.append(p['rhostart_1'])
        else:            
            rhostart.append(p['rhostart_'+str(iscuse)])
        Km1.append(p['Km1_'+str(iscuse)])
        Y1.append(p['Y1_'+str(iscuse)])
        Km2.append(p['Km2_'+str(iscuse)])
        try:
            Kml2.append(p['Kml2_'+str(iscuse)])
        except:
            Kml2.append(p['Km2_'+str(iscuse)])
        
        Y2.append(p['Y2_'+str(iscuse)])
        Km3.append(p['Km3_'+str(iscuse)])
        Y3.append(p['Y3_'+str(iscuse)])
        xhiconstant.append(p['xhiconstant_'+str(iscuse)]) #used for non-weber scaling
        
        if competitionmode>0: #set values of D and chi here
             
            if scalemodexhiD>=0:
                if len(epsiloncompetitionmode)==0:
                    print "define epsiloncompetitonmode"                    
                    
                    if n_strains==3:
                        Dlist=np.array([117.211721172,135.413541354,154.455445545])*np.power(10.,-12.)
                    elif n_strains==5:
                        Dlist=np.array([99.9499949995,117.211721172,135.413541354,154.455445545,174.297429743])*np.power(10.,-12.)
                    else:
                        errordefineDlistforthatstrainnumber
                else:
                    Dlist=np.array(epsiloncompetitionmode)*np.power(10.,-12.)
                    
                    
                D_chemotaxis.append(Dlist[isc-1])
                xhi.append(D_chemotaxis[-1]*scalemodexhiD)
                
            else:
                error_competitionmodeparametersonlydefinedforscalemodexhiD
            
        elif change_xhiB>0 and n_strains==2:
            xhi.append(p['xhi_'+str(1)]) #same for both strains....
            D_chemotaxis.append(p['D_chemotaxis_'+str(iscuse)])
        else:
            if scalexhiDwithgrowth>0:
                xhioverD=p['xhi_'+str(iscuse)]/p['D_chemotaxis_'+str(iscuse)]
                conversionfc=np.power(10.0,-12.)
                if scalexhiDwithgrowth==1:
                    D_chemotaxis.append((p['D_chemotaxis_offset_'+str(iscuse)]+p['D_chemotaxis_slope_'+str(iscuse)]*lambdac[isc-1]*3600.)*conversionfc)
                    xhi.append(xhioverD*D_chemotaxis[isc-1]) #used for weber scaling
                elif scalexhiDwithgrowth==2:
                    curslopec=p['D_chemotaxis_slope_'+str(iscuse)]
                    #print curslopec
                    #print p['D_chemotaxis_'+str(iscuse)]
                    #print p['D_chemotaxis_offset_'+str(iscuse)]
                    curoffsetc=p['D_chemotaxis_'+str(iscuse)]/conversionfc-curslopec*p['D_chemotaxis_offset_'+str(iscuse)] #offset parameter is taken as growth rate in 1/h
                    D_chemotaxis.append((curoffsetc+curslopec*lambdac[isc-1]*3600.)*conversionfc)
                    xhi.append(xhioverD*D_chemotaxis[isc-1]) #used for weber scaling
                    print "******* USED OFFSET, SLOP, XHI AND D VALUES*****"
                    print curoffsetc                    
                    print curslopec
                    print xhi
                    print D_chemotaxis
                    print "**********************************"
                else:
                    error 
            
            elif scalemodexhiD>=0:
                D_chemotaxis.append(p['D_chemotaxis_'+str(iscuse)])
                xhi.append(D_chemotaxis[-1]*scalemodexhiD)
            elif scalemodeDxhi>0:
                xhi.append(p['xhi_'+str(iscuse)]) #used for weber scaling
            
                D_chemotaxis.append(xhi[-1]*scalemodeDxhi)
            elif scalemodetwo>0:
                if isc==1:
                    D_chemotaxis.append(p['D_chemotaxis_'+str(iscuse)])
                    xhi.append(p['xhi_'+str(iscuse)]) #used for weber scaling
                else:
                    D_chemotaxis.append(p['D_chemotaxis_1']*scalemodetwo)
                    xhi.append(p['xhi_1']*scalemodetwo)
            else:
                D_chemotaxis.append(p['D_chemotaxis_'+str(iscuse)])
                xhi.append(p['xhi_'+str(iscuse)]) #used for weber scaling
            
        weber_offset.append(p['weber_offset_'+str(iscuse)])
        try:
            weber_offsetmax.append(p['weber_offsetmax_'+str(iscuse)])
        except:
            weber_offsetmax.append(1)
            
        if dynamicmode in ["AACarbonGRC","AACarbonGRC2","AACarbonGRC3","AACarbonGRC3starv1","AACarbonGRC3starv2","AACarbonGRC3starv3","AACarbonGRC3starv4","AACarbonGRC3logcut","AACarbonGRC3twoatr","AACarbonGRC3neqa","AACarbonGRC3neqalogcut","AACarbonGRC3NL","AACarbonGRC4","AACarbonGRClog","AACarbonGRClog2"]:
            Y2k1.append(p["Y2k1_"+str(iscuse)])
            Y2k2.append(p["Y2k2_"+str(iscuse)])
            
            if p['lambda_'+str(iscuse)]==-2:
                lambdac[isc-1]=Y2k2[-1]#add fit here....
            else:
                pass
                
            try:
                Y3k1.append(p["Y3k1_"+str(iscuse)])
                Y3k2.append(p["Y3k2_"+str(iscuse)])
            except:
                pass
        
        if tworingmode or dynamicmode in ["AACarbonGRC3twoatr"] :
            
                if scalemodexhiD>=0:
                    #D_chemotaxisTWO.append(p['D_chemotaxisTWO_'+str(iscuse)])
                    xhiTWO.append(D_chemotaxis[-1]*scalemodexhiD)
                elif scalemodeDxhi>0:
                    xhiTWO.append(p['xhiTWO_'+str(iscuse)]) #used for weber scaling
                
                    #D_chemotaxisTWO.append(xhiTWO[-1]*scalemodeDxhi)
                else:
                    #D_chemotaxisTWO.append(p['D_chemotaxisTWO_'+str(iscuse)])
                    xhiTWO.append(p['xhiTWO_'+str(iscuse)]) #used for weber scaling
            
                
        
        if dynamicmode in ["AACarbonGRClog","AACarbonGRClog2"]:
            try:
                CClog.append(p["CClog_"+str(iscuse)])
            except:
                CClog.append(1)
                
        elif dynamicmode in ["AACarbonGRClog2"]:
            CClog1.append(p["CClog1_"+str(iscuse)])
            deltal1.append(p['deltal1_'+str(iscuse)])
            deltal2.append(p['deltal2_'+str(iscuse)])
            hilllog.append(p['hilllog_'+str(iscuse)])
            hilllogCT.append(p['hilllog_'+str(iscuse)])

    print "testoutput used D chemotax sand xhi"
    print D_chemotaxis
    print xhi
    
    if fixed_totalrho>0 and n_strains>1:
        rhostart[0]=fixed_totalrho-rhostart[1]
        print "used rho1: "+str(rhostart[0])        
    if dynamicmode=="AACarbonLOGISTIC":
            rhomaxlogistic=p["rhomax_logistic"]
            
    #for evolution mode, overwrite some parameter settings to set wanted properties of competing strains
    if evolutionmode==True: 
        
        for isc in range(0,n_strains):

          if  evolutionvalue<0:
              
              if n_strains==5:
                  xcurlist=np.array([34.074, 78.58, 135.63, 202.2, 276.376])*np.power(10.,-12.0)
                  D_chemotaxis[isc]=D_chemotaxis[isc]#10.00*np.power(10.,-12.0)
                  xhi[isc]=xcurlist[isc]*scalemodexhiD
              elif n_strains==11:
                  xcurlist=[]
                  xcurlist=np.array([40.,56., 72., 88., 104., 120., 136., 152., 168., 184., 200.])*np.power(10.0,-12.0)
                  D_chemotaxis[isc]=xcurlist[isc]
                  xhi[isc]=xcurlist[isc]*scalemodexhiD
              else:
                  errorwrongstrainnumber
          
          else:
            
            if evolutionsinglemutation<0:   
                if evolutionvar=="xhi":
                    varc=evolutionvalue-evolutionstd+2*evolutionstd*isc/float(n_strains-1)
                    xhi[isc]=varc
                else:
                    print evolutionvarnotfound
            else:#set xhi values according to set values
                if isc<2:
                    if evolutionvar=="xhi":
                        varc=p["xhi_"+str(isc+1)]
                        
                        xhi[isc]=varc
                    else:
                        print evolutionvarnotfound
                else:
                    dynamicsforevolutionsinglemutationnotcompletlydefinedformorethantwostrains
       

    #read in initial conditions according to parameter setings
    nstart=[]
    for inn in range(1,n_nutrients+1):
        nstart.append(p['n'+str(inn)+'_start'])
    
    
    fisheralpha=p['fisheralpha']
    #setup needed variables
    n=[] #used for nutrients
    rho=[] #used for cell density
    #rhogf=[]
    m=[]
    phi=[]
    for inn in range(0,n_nutrients):
        phi.append([])
    
    if runmode in ['noweber','weber','weber2','weber-precise','weber-log','Keller']:
        
        
        #####
        #set initial and boundary conditions nutrient and attractant field
        ####
        for inn in range(0,n_nutrients):
            n.append(CellVariable(name="nutrient "+str(inn+1), mesh=mesh,hasOld=True,value=nstart[inn]))
            n.append(CellVariable(name="gradient "+str(inn+1), mesh=mesh,hasOld=True,value=nstart[inn])) #use to handle gradient
            
        if runmode =='weber-log':
            for inn in range(0,n_nutrients):
                for inns in range(0,n_strains):
                    phi[inn].append(CellVariable(name="phi "+str(inn+1), mesh=mesh,hasOld=True,value=0.))
                    
                    
        if initialnutrientprofile==1:
                
                widthprofile_nutrients
                for inn in range(0,n_nutrients):
                        valuecmax = nstart[inn]
                        for ix in range(0,nxc):
                            xc=xsimused[ix]
                            nccv=valuecmax*np.exp(-xc*xc/(2*widthprofile_nutrients*widthprofile_nutrients))
                            n[2*inn][ix]=nccv
        else:
            for inn in range(0,n_nutrients):
                n[2*inn].value = nstart[inn]
        
        
        #Set boundary conditions for different strains
        for isc in range(0,n_strains):
            rho.append(CellVariable(name="strain "+str(isc+1), mesh=mesh,hasOld=True,value=0.))        
            rho[isc].faceGrad.constrain([0.], mesh.facesLeft)      
            rho[isc].faceGrad.constrain([0.], mesh.facesRight)
            rho[isc].value = 0.
            
        if evolutionmode==False:
            
            print "run simulation not in evolution mode"
            nr=int(p['startring_radius']/dx)
            for isc in range(0,n_strains):
                if p['dimensionmode'] in ["one","one_explicit","1dradial","radial"]:
                    for nrx in range(1,1+nr,1):
                        rho[isc][nrx]=rhostart[isc]
                else:
                    errordimmodenotfound
                        #rho[1][nrx]=rho2start
                        #set linear gradient to test dynamics, remove later
            #generate nutrient fields so that plotting script is working
        else: #in evolution model, simulations are repeatedly run
                print "run simulation in evolution mode"
                #if there is already a prevous run, then initial composition is determined from previous run
                if loaddata==True:
                    rhoinitial=np.zeros([n_strains])
                    xctocheck=datastart[0]
                    if evolutionposition>=0: 
                     nrxcpos=xctocheck.shape[0]-2
                     for nrxc in range(0,xctocheck.shape[0]):
                        
                        if xctocheck[nrxc]>=evolutionposition:
                            nrxcpos=nrxc
                            break
                    elif evolutionposition==-3:
                        nrxcpos=xctocheck.shape[0]-2
                    else:
                        errorevolpostionvaluewrong
                        
                    print "last x"
                    print xctocheck[-1]
                    print "....loading from state before..."
                    print "evolution position"
                    print evolutionposition
                    print "index x"
                    print nrxcpos
                    print "used density"
                    print datastart[1].shape
                    print "density of different strains at that position..."
                    print datastart[1][0,0,:,-1,nrxcpos]
                    
                    if evolutionselectiondimension==0:

    
                        sumxrhoin=0
                        for nrxc in range(0,n_strains):
                            curvc=datastart[1][0,0,nrxc,-1,nrxcpos] #check, get value at right position.....
                            rhoinitial[nrxc]=curvc
                            sumxrhoin=sumxrhoin+curvc
                    elif evolutionselectiondimension==1:

                        #take along line until that position....
                        sumxrhoin=0
                        for nrxc in range(0,n_strains):
                            rhoinitial[nrxc]=0
                            for ixxc in range(1,nrxcpos):
                                curvc=datastart[1][0,0,nrxc,-1,ixxc] #check, get value at right position.....
                                rhoinitial[nrxc]=rhoinitial[nrxc]+curvc
                                sumxrhoin=sumxrhoin+curvc
                    elif evolutionselectiondimension==2: #in runs before this was ==1 and ==1 no (notweighted) was not included

                        #take along line until that position....
                        sumxrhoin=0
                        for nrxc in range(0,n_strains):
                            rhoinitial[nrxc]=0
                            for ixxc in range(1,nrxcpos):
                                curvc=datastart[1][0,0,nrxc,-1,ixxc]*xctocheck[ixxc] #check, get value at right position.....
                                rhoinitial[nrxc]=rhoinitial[nrxc]+curvc
                                sumxrhoin=sumxrhoin+curvc
                    elif evolutionselectiondimension==3:

                        #take along line until that position....
                        sumxrhoin=0
                        for nrxc in range(0,n_strains):
                            rhoinitial[nrxc]=0
                            for ixxc in range(1,nrxcpos):
                                curvc=datastart[1][0,0,nrxc,-1,ixxc]*xctocheck[ixxc]*xctocheck[ixxc] #check, get value at right position.....
                                rhoinitial[nrxc]=rhoinitial[nrxc]+curvc
                                sumxrhoin=sumxrhoin+curvc
                        
                        
                    #normalize
                    rhoinitial=evolutiondensity*rhoinitial/sumxrhoin
                        
                    
                    #set initial abundance
                    for nrxc in range(0,n_strains,1):
                        nr=int(p['startring_radius']/dx)
                        if p['dimensionmode'] in ["one","one_explicit","1dradial","radial"]:
                            for nrx in range(1,1+nr,1):
                                rho[nrxc][nrx]=rhoinitial[nrxc]
                        else:
                            errordimensionmodenotfound
                               
                    p['evolution_step']=evolutionstep+1#last point is end time of previous simulation
                    
                else:
                    print "start evolution without using previous simulation"
                    for nrxc in range(0,n_strains,1):
                        nr=int(p['startring_radius']/dx)
                        if p['dimensionmode'] in ["one","one_explicit","1dradial","radial"]:
                            if evolutionsinglemutation <0:
                                for nrx in range(1,1+nr,1):
                                    rho[nrxc][nrx]=evolutiondensity/float(n_strains)
                            else: #set last strain to single point mutation
                                if nrxc==n_strains-1:
                                    for nrx in range(1,1+nr,1):
                                        rho[nrxc][nrx]=evolutionsinglemutation
                                else:
                                    for nrx in range(1,1+nr,1):
                                        rho[nrxc][nrx]=evolutiondensity/float(n_strains-1)
                                    
                                    
                        else:
                            errordimensionmodenotfound
                        p['evolution_step']=0
                    
            
        
        #boundary condditions nutrients
        for inn in range(0,n_nutrients):
            n[inn*2].faceGrad.constrain([0.], mesh.facesLeft)    
            n[inn*2].faceGrad.constrain([0.], mesh.facesRight)
        
        #####################
        #update gradient
        #####################
        if runmode =='weber-log':
            for inn in range(0,n_nutrients):
                for inns in range(0,n_strains):
                    if dynamicmode in ["AACarbonGRC3neqa","AACarbonGRC3neqalogcut"]:
                        #phi[inn][inns][:]=weber_offsetlogalpha*np.log(weber_offset[inns]+n[0*inn]).faceGrad.value[0,:-1]
                        phi[inn][inns][:]=weber_offsetlogalpha*np.log((1.+n[0]/weber_offset[inns])/(1.+n[0]/weber_offsetmax[inns])).faceGrad.value[0,:-1]   
                        #n[inn+1]=n[inn].faceGrad.value[0,:-1]
                        phi[inn][inns][0]=0
                    
                    else:
                    #   phi[inn][inns][:]=weber_offsetlogalpha*np.log(weber_offset[inns]+n[2*inn]).faceGrad.value[0,:-1]
                        phi[inn][inns][:]=weber_offsetlogalpha*np.log((1.+n[2*inn]/weber_offset[inns])/(1.+n[2*inn]/weber_offsetmax[inns])).faceGrad.value[0,:-1]   
                        #n[inn+1]=n[inn].faceGrad.value[0,:-1]
                        phi[inn][inns][0]=0
                
        if usefacegrad:
            for inn in range(0,n_nutrients):
                n[2*inn+1][:]=n[2*inn].faceGrad.value[0,:-1]   
                #n[inn+1]=n[inn].faceGrad.value[0,:-1]
                n[2*inn+1][0]=0
            #n[1][:]=n[0].faceGrad.value[0,:-1]
        else:
            n[1][:]=n[0].grad.value[0,:] #several nutrients only defined for facegrad...
            n[1][0]=0    #correct initial value       
        
        
    if runmode in ['fisherwave']: 
        rho.append(CellVariable(name="strain "+str(0+1), mesh=mesh,hasOld=True,value=0.))                  
        rho[0].faceGrad.constrain([0.], mesh.facesLeft)      
        rho[0].faceGrad.constrain([0.], mesh.facesRight)
        rho[0].value = 0.
     
        nr=int(p['startring_radius']/dx)
        
        if p['dimensionmode'] in ["one","one_explicit","1dradial","radial"]:
                for nrx in range(1,1+nr,1):
                    rho[0][nrx]=rhostart[0]
                   
            
        n.append(CellVariable(name="nutrient 1", mesh=mesh,hasOld=True,value=0))
        n.append(np.zeros(nx))  #gradient field, not needed here

        
    length_n=len(n)
    length_rho=len(rho)
    length_m=len(m)
    

    ##############
    #setup final equations depending on mode
    ##############


    if runmode in ['fisherwave']:
        
        eqrho=TransientTerm(var=rho[0]) == DiffusionTerm(coeff=D_chemotaxis[0], var=rho[0])+ImplicitSourceTerm(lambdac[0]*(1.-fisheralpha*rho[0]),var=rho[0])
    if runmode in ['Keller']:
        
        if p['dimensionmode'] in ["one","two","radial"]:
            eqn1 = (TransientTerm(var=n[0]) == ImplicitSourceTerm((-1)*lambdac[0]*1./(n[0]+Km1[0])/Y1[0],var=n[0])+DiffusionTerm(coeff=D_nutrients,var=n[0]))
            eqrho=[]
            eqrho.append((TransientTerm(var=rho[0]) == DiffusionTerm(coeff=D_chemotaxis[0], var=rho[0])+ExponentialConvectionTerm(coeff=(-1)*n[3]*xhi[0]*(1./(n[2]+weber_offset[0]))*[[1.]],var=rho[0])))
          
        
    if runmode in ['noweber','weber','weber2','weber-precise','weber-log']:
        if p['dimensionmode'] in ["one","two","radial"]:
                 
            if dynamicmode=="standard":
                if n_strains==1: 
                    eqn1 = (TransientTerm(var=n[0]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km[0])/Y[0],var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[0]))
                elif n_strains==2:
                    eqn1 = (TransientTerm(var=n[0]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km[0])/Y[0],var=rho[0])+ImplicitSourceTerm((-1)*lambdac[1]*n[0]/(n[0]+Km[1])/Y[1],var=rho[1])+DiffusionTerm(coeff=D_nutrients,var=n[0]))
                if n_nutrients>1:
                    if n_strains==1: 
                        eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km[0])/Y[0],var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
                    elif n_strains==2:
                        eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km[0])/Y[0],var=rho[0])+ImplicitSourceTerm((-1)*lambdac[1]*n[0]/(n[0]+Km[1])/Y[1],var=rho[1])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
            elif dynamicmode=="AACarbon": #n_nutrients should be 2 anyway
                if n_strains==1: 
                    eqn1 = (TransientTerm(var=n[0]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])/Y1[0],var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[0]))
                    eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])/Y2[0]*n[2]/(n[2]+Km2[0]),var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
                elif n_strains==2: 
                    eqn1 = (TransientTerm(var=n[0]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])/Y1[0],var=rho[0])+ImplicitSourceTerm((-1)*lambdac[1]*n[0]/(n[0]+Km1[1])/Y1[1],var=rho[1])+DiffusionTerm(coeff=D_nutrients,var=n[0]))
                    eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])/Y2[0]*n[2]/(n[2]+Km2[0]),var=rho[0])+ImplicitSourceTerm((-1)*lambdac[1]*n[0]/(n[0]+Km1[1])/Y2[1]*n[2]/(n[2]+Km2[1]),var=rho[1])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
                elif n_strains==3:
                    eqn1 = (TransientTerm(var=n[0]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])/Y1[0],var=rho[0])+ImplicitSourceTerm((-1)*lambdac[1]*n[0]/(n[0]+Km1[1])/Y1[1],var=rho[1])+ImplicitSourceTerm((-1)*lambdac[2]*n[0]/(n[0]+Km1[2])/Y1[2],var=rho[2])+DiffusionTerm(coeff=D_nutrients,var=n[0]))
                    eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])/Y2[0]*n[2]/(n[2]+Km2[0]),var=rho[0])+ImplicitSourceTerm((-1)*lambdac[1]*n[0]/(n[0]+Km1[1])/Y2[1]*n[2]/(n[2]+Km2[1]),var=rho[1])+ImplicitSourceTerm((-1)*lambdac[2]*n[0]/(n[0]+Km1[2])/Y2[2]*n[2]/(n[2]+Km2[2]),var=rho[2])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
                elif n_strains==5:
                    eqn1 = (TransientTerm(var=n[0]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])/Y1[0],var=rho[0])+ImplicitSourceTerm((-1)*lambdac[1]*n[0]/(n[0]+Km1[1])/Y1[1],var=rho[1])+ImplicitSourceTerm((-1)*lambdac[2]*n[0]/(n[0]+Km1[2])/Y1[2],var=rho[2])+ImplicitSourceTerm((-1)*lambdac[3]*n[0]/(n[0]+Km1[3])/Y1[3],var=rho[3])+ImplicitSourceTerm((-1)*lambdac[4]*n[0]/(n[0]+Km1[4])/Y1[4],var=rho[4])+DiffusionTerm(coeff=D_nutrients,var=n[0]))
                    eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])/Y2[0]*n[2]/(n[2]+Km2[0]),var=rho[0])+ImplicitSourceTerm((-1)*lambdac[1]*n[0]/(n[0]+Km1[1])/Y2[1]*n[2]/(n[2]+Km2[1]),var=rho[1])+ImplicitSourceTerm((-1)*lambdac[2]*n[0]/(n[0]+Km1[2])/Y2[2]*n[2]/(n[2]+Km2[2]),var=rho[2])+ImplicitSourceTerm((-1)*lambdac[3]*n[0]/(n[0]+Km1[3])/Y2[3]*n[2]/(n[2]+Km2[3]),var=rho[3])+ImplicitSourceTerm((-1)*lambdac[4]*n[0]/(n[0]+Km1[4])/Y2[4]*n[2]/(n[2]+Km2[4]),var=rho[4])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
                elif n_strains==7:
                    eqn1 = (TransientTerm(var=n[0]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])/Y1[0],var=rho[0])+ImplicitSourceTerm((-1)*lambdac[1]*n[0]/(n[0]+Km1[1])/Y1[1],var=rho[1])+ImplicitSourceTerm((-1)*lambdac[2]*n[0]/(n[0]+Km1[2])/Y1[2],var=rho[2])+ImplicitSourceTerm((-1)*lambdac[3]*n[0]/(n[0]+Km1[3])/Y1[3],var=rho[3])+ImplicitSourceTerm((-1)*lambdac[4]*n[0]/(n[0]+Km1[4])/Y1[4],var=rho[4])+ImplicitSourceTerm((-1)*lambdac[5]*n[0]/(n[0]+Km1[5])/Y1[5],var=rho[5])+ImplicitSourceTerm((-1)*lambdac[6]*n[0]/(n[0]+Km1[6])/Y1[6],var=rho[6])+DiffusionTerm(coeff=D_nutrients,var=n[0]))
                    eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])/Y2[0]*n[2]/(n[2]+Km2[0]),var=rho[0])+ImplicitSourceTerm((-1)*lambdac[1]*n[0]/(n[0]+Km1[1])/Y2[1]*n[2]/(n[2]+Km2[1]),var=rho[1])+ImplicitSourceTerm((-1)*lambdac[2]*n[0]/(n[0]+Km1[2])/Y2[2]*n[2]/(n[2]+Km2[2]),var=rho[2])+ImplicitSourceTerm((-1)*lambdac[3]*n[0]/(n[0]+Km1[3])/Y2[3]*n[2]/(n[2]+Km2[3]),var=rho[3])+ImplicitSourceTerm((-1)*lambdac[4]*n[0]/(n[0]+Km1[4])/Y2[4]*n[2]/(n[2]+Km2[4]),var=rho[4])+ImplicitSourceTerm((-1)*lambdac[5]*n[0]/(n[0]+Km1[5])/Y2[5]*n[2]/(n[2]+Km2[5]),var=rho[5])+ImplicitSourceTerm((-1)*lambdac[6]*n[0]/(n[0]+Km1[6])/Y2[6]*n[2]/(n[2]+Km2[6]),var=rho[6])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
                
                elif n_strains==11:
                    eqn1 = (TransientTerm(var=n[0]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])/Y1[0],var=rho[0])+ImplicitSourceTerm((-1)*lambdac[1]*n[0]/(n[0]+Km1[1])/Y1[1],var=rho[1])+ImplicitSourceTerm((-1)*lambdac[2]*n[0]/(n[0]+Km1[2])/Y1[2],var=rho[2])+ImplicitSourceTerm((-1)*lambdac[3]*n[0]/(n[0]+Km1[3])/Y1[3],var=rho[3])+ImplicitSourceTerm((-1)*lambdac[4]*n[0]/(n[0]+Km1[4])/Y1[4],var=rho[4])+ImplicitSourceTerm((-1)*lambdac[5]*n[0]/(n[0]+Km1[5])/Y1[5],var=rho[5])+ImplicitSourceTerm((-1)*lambdac[6]*n[0]/(n[0]+Km1[6])/Y1[6],var=rho[6])+ImplicitSourceTerm((-1)*lambdac[7]*n[0]/(n[0]+Km1[7])/Y1[7],var=rho[7])+ImplicitSourceTerm((-1)*lambdac[8]*n[0]/(n[0]+Km1[8])/Y1[8],var=rho[8])+ImplicitSourceTerm((-1)*lambdac[9]*n[0]/(n[0]+Km1[9])/Y1[9],var=rho[9])+ImplicitSourceTerm((-1)*lambdac[10]*n[0]/(n[0]+Km1[10])/Y1[10],var=rho[10])+DiffusionTerm(coeff=D_nutrients,var=n[0]))
                    eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])/Y2[0]*n[2]/(n[2]+Km2[0]),var=rho[0])+ImplicitSourceTerm((-1)*lambdac[1]*n[0]/(n[0]+Km1[1])/Y2[1]*n[2]/(n[2]+Km2[1]),var=rho[1])+ImplicitSourceTerm((-1)*lambdac[2]*n[0]/(n[0]+Km1[2])/Y2[2]*n[2]/(n[2]+Km2[2]),var=rho[2])+ImplicitSourceTerm((-1)*lambdac[3]*n[0]/(n[0]+Km1[3])/Y2[3]*n[2]/(n[2]+Km2[3]),var=rho[3])+ImplicitSourceTerm((-1)*lambdac[4]*n[0]/(n[0]+Km1[4])/Y2[4]*n[2]/(n[2]+Km2[4]),var=rho[4])+ImplicitSourceTerm((-1)*lambdac[5]*n[0]/(n[0]+Km1[5])/Y2[5]*n[2]/(n[2]+Km2[5]),var=rho[5])+ImplicitSourceTerm((-1)*lambdac[6]*n[0]/(n[0]+Km1[6])/Y2[6]*n[2]/(n[2]+Km2[6]),var=rho[6])+ImplicitSourceTerm((-1)*lambdac[7]*n[0]/(n[0]+Km1[7])/Y2[7]*n[2]/(n[2]+Km2[7]),var=rho[7])+ImplicitSourceTerm((-1)*lambdac[8]*n[0]/(n[0]+Km1[8])/Y2[8]*n[2]/(n[2]+Km2[8]),var=rho[8])+ImplicitSourceTerm((-1)*lambdac[9]*n[0]/(n[0]+Km1[9])/Y2[9]*n[2]/(n[2]+Km2[9]),var=rho[9])+ImplicitSourceTerm((-1)*lambdac[10]*n[0]/(n[0]+Km1[10])/Y2[10]*n[2]/(n[2]+Km2[10]),var=rho[10])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
                else:
                    numberofstrainsunkown
            elif dynamicmode=="AACarbonROMAN": #n_nutrients should be 2 anyway
                if n_strains==1: 
                    eqn1 = (TransientTerm(var=n[0]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])/Y1[0],var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[0]))
                    eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])/Y2[0]*n[2]/(n[2]+Km2[0]),var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
                elif n_strains==2: 
                    eqn1 = (TransientTerm(var=n[0]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])/Y1[0],var=rho[0])+ImplicitSourceTerm((-1)*lambdac[1]*n[0]/(n[0]+Km1[1])/Y1[1],var=rho[1])+DiffusionTerm(coeff=D_nutrients,var=n[0]))
                    eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])/Y2[0]*n[2]/(n[2]+Km2[0]),var=rho[0])+ImplicitSourceTerm((-1)*lambdac[1]*n[0]/(n[0]+Km1[1])/Y2[1]*n[2]/(n[2]+Km2[1]),var=rho[1])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
                else:
                    numberofstrainsunkown
            elif dynamicmode=="AACarbonOD": #n_nutrients should be 2 anyway
                if n_strains==1: 
                    eqn1 = (TransientTerm(var=n[0]) == ImplicitSourceTerm((-1)*(lambda_back+lambda_frontd*(1-rho[0]*rho[0]/(lambdaodchange*lambdaodchange+rho[0]*rho[0])))*n[0]/(n[0]+Km1[0])/Y1[0],var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[0]))
                    eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*(lambda_back+lambda_frontd*(1-rho[0]*rho[0]/(lambdaodchange*lambdaodchange+rho[0]*rho[0])))*n[0]/(n[0]+Km1[0])/Y2[0]*n[2]/(n[2]+Km2[0]),var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
                elif n_strains==2: 
                    eqn1 = (TransientTerm(var=n[0]) == ImplicitSourceTerm((-1)*(lambda_back+lambda_frontd*(1-(rho[0]+rho[1])*(rho[0]+rho[1])/(lambdaodchange*lambdaodchange+(rho[0]+rho[1])*(rho[0]+rho[1]))))*n[0]/(n[0]+Km1[0])/Y1[0],var=rho[0])+ImplicitSourceTerm((-1)*(lambda_back+lambda_frontd*(1-(rho[0]+rho[1])*(rho[0]+rho[1])/(lambdaodchange*lambdaodchange+(rho[0]+rho[1])*(rho[0]+rho[1]))))*n[0]/(n[0]+Km1[1])/Y1[1],var=rho[1])+DiffusionTerm(coeff=D_nutrients,var=n[0]))
                    eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*(lambda_back+lambda_frontd*(1-(rho[0]+rho[1])*(rho[0]+rho[1])/(lambdaodchange*lambdaodchange+(rho[0]+rho[1])*(rho[0]+rho[1]))))*n[0]/(n[0]+Km1[0])/Y2[0]*n[2]/(n[2]+Km2[0]),var=rho[0])+ImplicitSourceTerm((-1)*(lambda_back+lambda_frontd*(1-(rho[0]+rho[1])*(rho[0]+rho[1])/(lambdaodchange*lambdaodchange+(rho[0]+rho[1])*(rho[0]+rho[1]))))*n[0]/(n[0]+Km1[1])/Y2[1]*n[2]/(n[2]+Km2[1]),var=rho[1])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
                else:
                    numberofstrainsunkown
            elif dynamicmode=="AACarbonOD4": #n_nutrients should be 2 anyway
                if n_strains==1: 
                    eqn1 = (TransientTerm(var=n[0]) == ImplicitSourceTerm((-1)*(lambda_back+lambda_frontd*(1-rho[0]*rho[0]*rho[0]*rho[0]/(lambdaodchange*lambdaodchange*lambdaodchange*lambdaodchange+rho[0]*rho[0]*rho[0]*rho[0])))*n[0]/(n[0]+Km1[0])/Y1[0],var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[0]))
                    eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*(lambda_back+lambda_frontd*(1-rho[0]*rho[0]*rho[0]*rho[0]/(lambdaodchange*lambdaodchange*lambdaodchange*lambdaodchange+rho[0]*rho[0]*rho[0]*rho[0])))*n[0]/(n[0]+Km1[0])/Y2[0]*n[2]/(n[2]+Km2[0]),var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
                elif n_strains==2: 
                    eqn1 = (TransientTerm(var=n[0]) == ImplicitSourceTerm((-1)*(lambda_back+lambda_frontd*(1-(rho[0]+rho[1])*(rho[0]+rho[1])*(rho[0]+rho[1])*(rho[0]+rho[1])/(lambdaodchange*lambdaodchange*lambdaodchange*lambdaodchange+(rho[0]+rho[1])*(rho[0]+rho[1])*(rho[0]+rho[1])*(rho[0]+rho[1]))))*n[0]/(n[0]+Km1[0])/Y1[0],var=rho[0])+ImplicitSourceTerm((-1)*(lambda_back+lambda_frontd*(1-(rho[0]+rho[1])*(rho[0]+rho[1])*(rho[0]+rho[1])*(rho[0]+rho[1])/(lambdaodchange*lambdaodchange*lambdaodchange*lambdaodchange+(rho[0]+rho[1])*(rho[0]+rho[1])*(rho[0]+rho[1])*(rho[0]+rho[1]))))*n[0]/(n[0]+Km1[1])/Y1[1],var=rho[1])+DiffusionTerm(coeff=D_nutrients,var=n[0]))
                    eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*(lambda_back+lambda_frontd*(1-(rho[0]+rho[1])*(rho[0]+rho[1])*(rho[0]+rho[1])*(rho[0]+rho[1])/(lambdaodchange*lambdaodchange*lambdaodchange*lambdaodchange+(rho[0]+rho[1])*(rho[0]+rho[1])*(rho[0]+rho[1])*(rho[0]+rho[1]))))*n[0]/(n[0]+Km1[0])/Y2[0]*n[2]/(n[2]+Km2[0]),var=rho[0])+ImplicitSourceTerm((-1)*(lambda_back+lambda_frontd*(1-(rho[0]+rho[1])*(rho[0]+rho[1])*(rho[0]+rho[1])*(rho[0]+rho[1])/(lambdaodchange*lambdaodchange*lambdaodchange*lambdaodchange+(rho[0]+rho[1])*(rho[0]+rho[1])*(rho[0]+rho[1])*(rho[0]+rho[1]))))*n[0]/(n[0]+Km1[1])/Y2[1]*n[2]/(n[2]+Km2[1]),var=rho[1])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
                else:
                    numberofstrainsunkown
                    
                    
         
            elif dynamicmode=="AACarbonLOGISTIC": #n_nutrients should be 2 anyway
                if n_strains==1: 
                    eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*lambdac[0]*(1.-rho[0]/rhomaxlogistic)/Y2[0]*n[2]/(n[2]+Km2[0]),var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
                elif n_strains==2: 
                    eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*lambdac[0]*(1.-(rho[0]+rho[1])/rhomaxlogistic)/Y2[0]*n[2]/(n[2]+Km2[0]),var=rho[0])+ImplicitSourceTerm((-1)*lambdac[1]*(1.-(rho[0]+rho[1])/rhomaxlogistic)/Y2[1]*n[2]/(n[2]+Km2[1]),var=rho[1])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
                elif n_strains==5:
                    eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*lambdac[0]*(1.-(rho[0]+rho[1]+rho[2]+rho[3]+rho[4])/rhomaxlogistic)/Y2[0]*n[2]/(n[2]+Km2[0]),var=rho[0])+ImplicitSourceTerm((-1)*lambdac[1]*(1.-(rho[0]+rho[1]+rho[2]+rho[3]+rho[4])/rhomaxlogistic)/Y2[1]*n[2]/(n[2]+Km2[1]),var=rho[1])+ImplicitSourceTerm((-1)*lambdac[2]*(1.-(rho[0]+rho[1]+rho[2]+rho[3]+rho[4])/rhomaxlogistic)/Y2[2]*n[2]/(n[2]+Km2[2]),var=rho[2])+ImplicitSourceTerm((-1)*lambdac[3]*(1.-(rho[0]+rho[1]+rho[2]+rho[3]+rho[4])/rhomaxlogistic)/Y2[3]*n[2]/(n[2]+Km2[3]),var=rho[3])+ImplicitSourceTerm((-1)*lambdac[4]*(1.-(rho[0]+rho[1]+rho[2]+rho[3]+rho[4])/rhomaxlogistic)/Y2[4]*n[2]/(n[2]+Km2[4]),var=rho[4])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
                elif n_strains==11:
                    eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*lambdac[0]*(1.-(rho[0]+rho[1]+rho[2]+rho[3]+rho[4]+rho[5]+rho[6]+rho[7]+rho[8]+rho[9]+rho[10])/rhomaxlogistic)/Y2[0]*n[2]/(n[2]+Km2[0]),var=rho[0])+ImplicitSourceTerm((-1)*lambdac[1]*(1.-(rho[0]+rho[1]+rho[2]+rho[3]+rho[4]+rho[5]+rho[6]+rho[7]+rho[8]+rho[9]+rho[10]))/Y2[1]*n[2]/(n[2]+Km2[1]),var=rho[1])+ImplicitSourceTerm((-1)*lambdac[2]*(1.-(rho[0]+rho[1]+rho[2]+rho[3]+rho[4]+rho[5]+rho[6]+rho[7]+rho[8]+rho[9]+rho[10]))/Y2[2]*n[2]/(n[2]+Km2[2]),var=rho[2])+ImplicitSourceTerm((-1)*lambdac[3]*(1.-(rho[0]+rho[1]+rho[2]+rho[3]+rho[4]+rho[5]+rho[6]+rho[7]+rho[8]+rho[9]+rho[10]))/Y2[3]*n[2]/(n[2]+Km2[3]),var=rho[3])+ImplicitSourceTerm((-1)*lambdac[4]*(1.-(rho[0]+rho[1]+rho[2]+rho[3]+rho[4]+rho[5]+rho[6]+rho[7]+rho[8]+rho[9]+rho[10]))/Y2[4]*n[2]/(n[2]+Km2[4]),var=rho[4])+ImplicitSourceTerm((-1)*lambdac[5]*(1.-(rho[0]+rho[1]+rho[2]+rho[3]+rho[4]+rho[5]+rho[6]+rho[7]+rho[8]+rho[9]+rho[10]))/Y2[5]*n[2]/(n[2]+Km2[5]),var=rho[5])+ImplicitSourceTerm((-1)*lambdac[6]*(1.-(rho[0]+rho[1]+rho[2]+rho[3]+rho[4]+rho[5]+rho[6]+rho[7]+rho[8]+rho[9]+rho[10]))/Y2[6]*n[2]/(n[2]+Km2[6]),var=rho[6])+ImplicitSourceTerm((-1)*lambdac[7]*(1.-(rho[0]+rho[1]+rho[2]+rho[3]+rho[4]+rho[5]+rho[6]+rho[7]+rho[8]+rho[9]+rho[10]))/Y2[7]*n[2]/(n[2]+Km2[7]),var=rho[7])+ImplicitSourceTerm((-1)*lambdac[8]*(1.-(rho[0]+rho[1]+rho[2]+rho[3]+rho[4]+rho[5]+rho[6]+rho[7]+rho[8]+rho[9]+rho[10]))/Y2[8]*n[2]/(n[2]+Km2[8]),var=rho[8])+ImplicitSourceTerm((-1)*lambdac[9]*(1.-(rho[0]+rho[1]+rho[2]+rho[3]+rho[4]+rho[5]+rho[6]+rho[7]+rho[8]+rho[9]+rho[10]))/Y2[9]*n[2]/(n[2]+Km2[9]),var=rho[9])+ImplicitSourceTerm((-1)*lambdac[10]*(1.-(rho[0]+rho[1]+rho[2]+rho[3]+rho[4]+rho[5]+rho[6]+rho[7]+rho[8]+rho[9]+rho[10]))/Y2[10]*n[2]/(n[2]+Km2[10]),var=rho[10])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
                else:
                    numberofstrainsunkown
                                    
                    
            elif dynamicmode=="AACarbonGRC": #n_nutrients should be 2 anyway
                if n_strains==1: 
                    eqn1 = (TransientTerm(var=n[0]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])/Y1[0],var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[0]))
                    eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])*(Y2k1[0]+Y2k2[0]*lambdac[0]*n[0]/(n[0]+Km1[0]))*n[2]/(n[2]+Km2[0]),var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
                elif n_strains==2: 
                    eqn1 = (TransientTerm(var=n[0]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])/Y1[0],var=rho[0])+ImplicitSourceTerm((-1)*lambdac[1]*n[0]/(n[0]+Km1[1])/Y1[1],var=rho[1])+DiffusionTerm(coeff=D_nutrients,var=n[0]))
                    eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])*(Y2k1[0]+Y2k2[0]*lambdac[0]*n[0]/(n[0]+Km1[0]))*n[2]/(n[2]+Km2[0]),var=rho[0])+ImplicitSourceTerm((-1)*lambdac[1]*n[0]/(n[0]+Km1[1])*(Y2k1[1]+Y2k2[1]*lambdac[1]*n[0]/(n[0]+Km1[1]))*n[2]/(n[2]+Km2[1]),var=rho[1])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
                else:
                    errorc
            elif dynamicmode=="AACarbonGRC2": #n_nutrients should be 2 anyway
                if n_strains==1: 
                    eqn1 = (TransientTerm(var=n[0]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])/Y1[0],var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[0]))
                    eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*(lambdac[0]+lambdadeltac[0]*n[2]/(n[2]+Kml2[0]))*n[0]/(n[0]+Km1[0])*(Y2k1[0]+Y2k2[0]*(lambdac[0]+lambdadeltac[0]*n[2]/(n[2]+Kml2[0]))*n[0]/(n[0]+Km1[0]))*n[2]/(n[2]+Km2[0]),var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
                #elif n_strains==2: 
                #    eqn1 = (TransientTerm(var=n[0]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])/Y1[0],var=rho[0])+ImplicitSourceTerm((-1)*lambdac[1]*n[0]/(n[0]+Km1[1])/Y1[1],var=rho[1])+DiffusionTerm(coeff=D_nutrients,var=n[0]))
                #    eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])*(Y2k1[0]+Y2k2[0]*lambdac[0]*n[0]/(n[0]+Km1[0]))*n[2]/(n[2]+Km2[0]),var=rho[0])+ImplicitSourceTerm((-1)*lambdac[1]*n[0]/(n[0]+Km1[1])*(Y2k1[1]+Y2k2[1]*lambdac[1]*n[0]/(n[0]+Km1[1]))*n[2]/(n[2]+Km2[1]),var=rho[1])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
                else:
                    errorc
            elif dynamicmode=="AACarbonGRC3": #n_nutrients should be 2 anyway
                if tworingmode and n_strains==1:
                    eqn1 = (TransientTerm(var=n[0]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])/Y1[0],var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[0]))
                    eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*(Y2k1[0]*lambdac[0]*n[0]/(n[0]+Km1[0])+Y2k2[0])*n[2]/(n[2]+Km2[0]),var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
                    eqn3 = (TransientTerm(var=n[4]) == ImplicitSourceTerm((-1)*(Y3k1[0]*lambdac[0]*n[0]/(n[0]+Km1[0])+Y3k2[0])*n[4]/(n[4]+Km3[0]),var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[4]))
               
                elif n_strains==1: 
                    eqn1 = (TransientTerm(var=n[0]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])/Y1[0],var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[0]))
                    eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*(Y2k1[0]*lambdac[0]*n[0]/(n[0]+Km1[0])+Y2k2[0])*n[2]/(n[2]+Km2[0]),var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
                else:
                    errorc
            elif dynamicmode=="AACarbonGRC3starv1": #n_nutrients should be 2 anyway
                if n_strains==1: 
                    eqn1 = (TransientTerm(var=n[0]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])/Y1[0],var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[0]))
                    eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*(Y2k1[0]*lambdac[0]*n[0]/(n[0]+Km1[0])+Y2k2[0])*n[2]/(n[2]+Km2[0]),var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
                else:
                    errorc
            elif dynamicmode=="AACarbonGRC3starv3": #n_nutrients should be 2 anyway
                if n_strains==1: 
                    eqn1 = (TransientTerm(var=n[0]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])/Y1[0],var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[0]))
                    eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*(Y2k1[0]*lambdac[0]*n[0]/(n[0]+Km1[0])+Y2k2[0])*n[2]/(n[2]+Km2[0]),var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
                else:
                    errorc
            elif dynamicmode=="AACarbonGRC3starv2": #n_nutrients should be 2 anyway
                if n_strains==1: 
                    eqn1 = (TransientTerm(var=n[0]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])/Y1[0],var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[0]))
                    eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*(Y2k1[0]*lambdac[0]*n[0]/(n[0]+Km1[0])+Y2k2[0])*n[2]/(n[2]+Km2[0]),var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
                else:
                    errorc
            elif dynamicmode=="AACarbonGRC3starv4": #n_nutrients should be 2 anyway
                if n_strains==1: 
                    eqn1 = (TransientTerm(var=n[0]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])/Y1[0],var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[0]))
                    eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*(Y2k1[0]*lambdac[0]*n[0]/(n[0]+Km1[0])+Y2k2[0])*n[2]/(n[2]+Km2[0]),var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
                else:
                    errorc
                        
            elif dynamicmode=="AACarbonGRC3logcut": #n_nutrients should be 2 anyway
                if n_strains==1: 
                    eqn1 = (TransientTerm(var=n[0]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])/Y1[0],var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[0]))
                    eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*(Y2k1[0]*lambdac[0]*n[0]/(n[0]+Km1[0])+Y2k2[0])*n[2]/(n[2]+Km2[0]),var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
                else:
                    errorc
           
            elif dynamicmode=="AACarbonGRC3twoatr": #n_nutrients should be 2 anyway
                if n_strains==1:
                    eqn1 = (TransientTerm(var=n[0]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])/Y1[0],var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[0]))
                    eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*(Y2k1[0]*lambdac[0]*n[0]/(n[0]+Km1[0])+Y2k2[0])*n[2]/(n[2]+Km2[0]),var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
                    
                else:
                    errorc
           
            elif dynamicmode=="AACarbonGRC3neqa": #n_nutrients should be 2 anyway
                
                    eqn1 = (TransientTerm(var=n[0]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])/Y1[0],var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[0]))
                                       
                   
            elif dynamicmode=="AACarbonGRC3neqalogcut": #n_nutrients should be 2 anyway
                
                    eqn1 = (TransientTerm(var=n[0]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])/Y1[0],var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[0]))
                
            
            elif dynamicmode=="AACarbonGRC3NL": #n_nutrients should be 2 anyway
                if n_strains==1: 
                    eqn1 = (TransientTerm(var=n[0]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])/Y1[0],var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[0]))
                    eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*(Y2k1[0]*lambdac[0]*n[0]/(n[0]+Km1[0])+Y2k2[0])*n[2]/(n[2]+Km2[0]),var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
                else:
                    errorc
            elif dynamicmode=="AACarbonGRC4": #onlty attractant...
                if n_strains==1: 
                    eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*(Y2k1[0])*n[2]/(n[2]+Km2[0]),var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
                else:
                    errorc
            
                        
            elif dynamicmode=="AACarbonGRClog": #n_nutrients should be 2 anyway
                if n_strains==1: 
                    eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*lambdac[0]*(1.-rho[0]/CClog[0])*(Y2k1[0]+Y2k2[0]*lambdac[0]*(1.-rho[0]/CClog[0]))*n[2]/(n[2]+Km2[0]),var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
                else:
                    errorc
            elif dynamicmode=="AACarbonGRClog2": #n_nutrients should be 2 anyway
                if n_strains==1: 
                    eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*(lambdac[0]*np.power(CClog[0],hilllog[0])/(np.power(CClog[0],hilllog[0])+np.power(rho[0],hilllog[0]))+deltal1[0]*np.power(CClog1[0],hilllog[0])/(np.power(CClog1[0],hilllog[0])+np.power(rho[0],hilllog[0]))+deltal2[0]*np.power(Km2[0],hilllogCT[0])/(np.power(Km2[0],hilllogCT[0])+np.power(n[2],hilllogCT[0]))
                    )*(Y2k1[0]+Y2k2[0]*(lambdac[0]*np.power(CClog[0],hilllog[0])/(np.power(CClog[0],hilllog[0])+np.power(rho[0],hilllog[0]))+deltal1[0]*np.power(CClog1[0],hilllog[0])/(np.power(CClog1[0],hilllog[0])+np.power(rho[0],hilllog[0]))+deltal2[0]*np.power(Km2[0],hilllogCT[0])/(np.power(Km2[0],hilllogCT[0])+np.power(n[2],hilllogCT[0]))
                    ))*n[2]/(n[2]+Km2[0]),var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
                else:
                    errorc
                    
                    
               
                
            else:
                errordynamicmodenotknown
                
            eqrho=[]            
            for isc in range(0,n_strains):   
                
                if runmode=='weber':
                    if dynamicmode=="standard":
                        
                        eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*n[1]*xhi[isc]*(1./(n[0]+weber_offset))*[[1.]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*n[0]/(n[0]+Km[isc]),var=rho[isc])))
                    elif dynamicmode=="AACarbon":
                        
                        eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*n[3]*xhi[isc]*(1./(n[2]+weber_offset[isc]))*[[1.]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                    elif dynamicmode=="AACarbonROMAN":
                        
                        eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*(xhi[isc]*np.tanh(n[3]/weber_offset[isc])*(n[3]/np.abs(n[3])))*[[1.]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                    
                    elif dynamicmode=="AACarbonOD":
                        if n_strains==1:
                            #growth rate (lambda_back+lambda_frontd*(1-rho[0]*rho[0]/(lambdaodchange*lambdaodchange+rho[0]*rho[0])))
                            eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*n[3]*xhi[isc]*(1./(n[2]+weber_offset[isc]))*[[1.]],var=rho[isc])+ImplicitSourceTerm((lambda_back+lambda_frontd*(1-rho[0]*rho[0]/(lambdaodchange*lambdaodchange+rho[0]*rho[0])))*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                        elif n_strains==2:
                            #growth rate (lambda_back+lambda_frontd*(1-rho[0]*rho[0]/(lambdaodchange*lambdaodchange+rho[0]*rho[0])))
                            eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*n[3]*xhi[isc]*(1./(n[2]+weber_offset[isc]))*[[1.]],var=rho[isc])+ImplicitSourceTerm((lambda_back+lambda_frontd*(1-(rho[0]+rho[1])*(rho[0]+rho[1])/(lambdaodchange*lambdaodchange+(rho[0]+rho[1])*(rho[0]+rho[1]))))*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                        else:
                            errornstrains
                    elif dynamicmode=="AACarbonOD4":
                        if n_strains==1:
                            #growth rate (lambda_back+lambda_frontd*(1-rho[0]*rho[0]/(lambdaodchange*lambdaodchange+rho[0]*rho[0])))
                            eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*n[3]*xhi[isc]*(1./(n[2]+weber_offset[isc]))*[[1.]],var=rho[isc])+ImplicitSourceTerm((lambda_back+lambda_frontd*(1-rho[0]*rho[0]*rho[0]*rho[0]/(lambdaodchange*lambdaodchange*lambdaodchange*lambdaodchange+rho[0]*rho[0]*rho[0]*rho[0])))*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                        elif n_strains==2:
                            #growth rate (lambda_back+lambda_frontd*(1-rho[0]*rho[0]/(lambdaodchange*lambdaodchange+rho[0]*rho[0])))
                            eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*n[3]*xhi[isc]*(1./(n[2]+weber_offset[isc]))*[[1.]],var=rho[isc])+ImplicitSourceTerm((lambda_back+lambda_frontd*(1-(rho[0]+rho[1])*(rho[0]+rho[1])*(rho[0]+rho[1])*(rho[0]+rho[1])/(lambdaodchange*lambdaodchange*lambdaodchange*lambdaodchange+(rho[0]+rho[1])*(rho[0]+rho[1])*(rho[0]+rho[1])*(rho[0]+rho[1]))))*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                        else:
                            errornstrains
                     
                      
                    elif dynamicmode=="AACarbonGRC":                        
                        eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*n[3]*xhi[isc]*(1./(n[2]+weber_offset[isc]))*[[1.]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*n[0]/(n[0]+Km1[isc]),var=rho[isc])))                                        
                        
                    
                    else:
                        errordynamicmodenotknown
                elif runmode=='weber2':
                    if dynamicmode=="standard":
                        
                        eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*n[1]*xhi[isc]*(1./(n[0]*n[0]+weber_offset))*[[1.]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*n[0]/(n[0]+Km[isc]),var=rho[isc])))
                    elif dynamicmode=="AACarbon":
                        
                        eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*n[3]*xhi[isc]*(1./(n[2]*n[2]+weber_offset[isc]))*[[1.]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                    elif dynamicmode=="AACarbonGRC":                        
                        eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*n[3]*xhi[isc]*(1./(n[2]*n[2]+weber_offset[isc]))*[[1.]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*n[0]/(n[0]+Km1[isc]),var=rho[isc])))                                        
                    else:
                        errordynamicmodenotknown
                        
                elif runmode=='noweber':
                    if dynamicmode=="standard":
                        eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*xhiconstant[isc]*n[1]*[[1]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*n[0]/(n[0]+Km[isc]),var=rho[isc])))
                    elif dynamicmode=="AACarbon":
                        eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*xhiconstant[isc]*n[3]*[[1]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                    else:
                        errordynamicmodenotknown
                elif runmode=='weber-log':
                    if dynamicmode=="AACarbon":
                        eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*xhi[isc]*phi[1][isc]*[[1]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                    elif dynamicmode=="AACarbonROMAN":
                        eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*(xhi[isc]*np.tanh(n[3]/weber_offset[isc]))*[[1.]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                    
                    elif dynamicmode=="AACarbonOD":
                        if n_strains==1:
                            #growth rate (lambda_back+lambda_frontd*(1-rho[0]*rho[0]/(lambdaodchange*lambdaodchange+rho[0]*rho[0])))
                            eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*xhi[isc]*phi[1][isc]*[[1]],var=rho[isc])+ImplicitSourceTerm((lambda_back+lambda_frontd*(1-rho[0]*rho[0]/(lambdaodchange*lambdaodchange+rho[0]*rho[0])))*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                     
                        elif n_strains==2:
                            #growth rate (lambda_back+lambda_frontd*(1-rho[0]*rho[0]/(lambdaodchange*lambdaodchange+rho[0]*rho[0])))
                            eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*xhi[isc]*phi[1][isc]*[[1]],var=rho[isc])+ImplicitSourceTerm((lambda_back+lambda_frontd*(1-(rho[0]+rho[1])*(rho[0]+rho[1])/(lambdaodchange*lambdaodchange+(rho[0]+rho[1])*(rho[0]+rho[1]))))*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                    
                        else:
                            errornstrains
                    elif dynamicmode=="AACarbonOD4":
                        if n_strains==1:
                            #growth rate (lambda_back+lambda_frontd*(1-rho[0]*rho[0]/(lambdaodchange*lambdaodchange+rho[0]*rho[0])))
                            eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*xhi[isc]*phi[1][isc]*[[1]],var=rho[isc])+ImplicitSourceTerm((lambda_back+lambda_frontd*(1-rho[0]*rho[0]*rho[0]*rho[0]/(lambdaodchange*lambdaodchange*lambdaodchange*lambdaodchange+rho[0]*rho[0]*rho[0]*rho[0])))*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                     
                        elif n_strains==2:
                            #growth rate (lambda_back+lambda_frontd*(1-rho[0]*rho[0]/(lambdaodchange*lambdaodchange+rho[0]*rho[0])))
                            eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*xhi[isc]*phi[1][isc]*[[1]],var=rho[isc])+ImplicitSourceTerm((lambda_back+lambda_frontd*(1-(rho[0]+rho[1])*(rho[0]+rho[1])*(rho[0]+rho[1])*(rho[0]+rho[1])/(lambdaodchange*lambdaodchange*lambdaodchange*lambdaodchange+(rho[0]+rho[1])*(rho[0]+rho[1])*(rho[0]+rho[1])*(rho[0]+rho[1]))))*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                    
                        else:
                            errornstrains
                     
                     
                     
                        
                    elif dynamicmode=="AACarbonGRC":
                        if tworingmode and n_strains==1:
                            eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*(xhi[isc]*phi[1][isc]*[[1]]+xhiTWO[isc]*phi[2][isc]*[[1]]),var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                            
                        elif n_strains==1:
                            eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*xhi[isc]*phi[1][isc]*[[1]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                        else:
                            error
                            
                    
                    
                    
                    elif dynamicmode=="AACarbonGRC3":
                        if tworingmode and n_strains==1:
                            eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*(xhi[isc]*phi[1][isc]*[[1]]+xhiTWO[isc]*phi[2][isc]*[[1]]),var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                        elif n_strains==1:
                            eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*xhi[isc]*phi[1][isc]*[[1]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                        else:
                            error
                    elif dynamicmode=="AACarbonGRC3starv1":
                        
                        if n_strains==1:
                            eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*np.clip((-41.4219114219*rho[isc]*rho[isc]+49.3244755245*rho[isc]+3.48666666667)*(-41.4219114219*rho[isc]*rho[isc]+49.3244755245*rho[isc]+3.48666666667)/18.1702564531/18.1702564531,0,10000)*xhi[isc]*phi[1][isc]*[[1]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                            #eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*xhi[isc]*phi[1][isc]*[[1]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                        
                        else:
                            error
                    elif dynamicmode=="AACarbonGRC3starv3":
                        if n_strains==1:
                            eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=(np.clip((-41.4219114219*rho[isc]*rho[isc]+49.3244755245*rho[isc]+3.48666666667)*(-41.4219114219*rho[isc]*rho[isc]+49.3244755245*rho[isc]+3.48666666667)/(18.1702564531*18.1702564531),0,10000))*D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*np.clip((-41.4219114219*rho[isc]*rho[isc]+49.3244755245*rho[isc]+3.48666666667)*(-41.4219114219*rho[isc]*rho[isc]+49.3244755245*rho[isc]+3.48666666667)/(18.1702564531*18.1702564531),0,10000)*xhi[isc]*phi[1][isc]*[[1]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                        else:
                            error
                    elif dynamicmode=="AACarbonGRC3starv2":
                        if n_strains==1:
                            lmaxcstarv=0.69/3600.
                            nhccstarv1=5.0
                            eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*(1.-np.power(lambdac[isc]*n[0]/(n[0]+Km1[isc]),nhccstarv1)/(np.power(0.5*lmaxcstarv,nhccstarv1)+np.power(lambdac[isc]*n[0]/(n[0]+Km1[isc]),nhccstarv1)))*xhi[isc]*phi[1][isc]*[[1]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                        else:
                            error
                    elif dynamicmode=="AACarbonGRC3starv4":
                        if n_strains==1:
                            lmaxcstarv=0.69/3600.
                            nhccstarv1=5.0
                            eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=(1.-np.power(lambdac[isc]*n[0]/(n[0]+Km1[isc]),nhccstarv1)/(np.power(0.5*lmaxcstarv,nhccstarv1)+np.power(lambdac[isc]*n[0]/(n[0]+Km1[isc]),nhccstarv1)))*D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*(1.-np.power(lambdac[isc]*n[0]/(n[0]+Km1[isc]),nhccstarv1)/(np.power(0.5*lmaxcstarv,nhccstarv1)+np.power(lambdac[isc]*n[0]/(n[0]+Km1[isc]),nhccstarv1)))*xhi[isc]*phi[1][isc]*[[1]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                            #eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=*D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*np.clip((-41.4219114219*rho[isc]*rho[isc]+49.3244755245*rho[isc]+3.48666666667)*(-41.4219114219*rho[isc]*rho[isc]+49.3244755245*rho[isc]+3.48666666667)/(18.1702564531*18.1702564531),0,10000)*xhi[isc]*phi[1][isc]*[[1]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                       
                        else:
                            error
                   
                    elif dynamicmode=="AACarbonGRC3logcut":
                        if n_strains==1:
                            rhomaxcc=odcuttofflogistic
                            eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*xhi[isc]*phi[1][isc]*[[1]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*(1-rho[isc]/rhomaxcc)*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                        else:
                            error
                    elif dynamicmode=="AACarbonGRC3twoatr":
                        if n_strains==1:
                            eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*(xhi[isc]*phi[1][isc]*[[1]]+xhiTWO[isc]*phi[0][isc]*[[1]]),var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                        #elif n_strains==1:
                        #    eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*xhi[isc]*phi[1][isc]*[[1]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                        else:
                            error
                    elif dynamicmode=="AACarbonGRC3neqa":
                        if n_strains==1:
                            eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*xhi[isc]*phi[0][isc]*[[1]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                        else:
                            error
                                            
                            
                    elif dynamicmode=="AACarbonGRC3neqalogcut":
                        if n_strains==1:
                            rhomaxcc=odcuttofflogistic
                            eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*xhi[isc]*phi[0][isc]*[[1]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*(1-rho[isc]/rhomaxcc)*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                            
                        else:
                            error
                        
                            
                            
                    
                    elif dynamicmode=="AACarbonGRC3NL":
                        if tworingmode and n_strains==1:
                            eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*(xhi[isc]*phi[1][isc]*[[1]]+xhiTWO[isc]*phi[2][isc]*[[1]]),var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                        elif n_strains==1:
                            eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*xhi[isc]*phi[1][isc]*[[1]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                        else:
                            error
               #eli                   
                    elif dynamicmode=="AACarbonGRC4":
                        eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*xhi[isc]*phi[1][isc]*[[1]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*n[2]/(n[2]+Km2[isc]),var=rho[isc])))
                                        
                    elif dynamicmode=="AACarbonGRC2":
                        eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*xhi[isc]*phi[1][isc]*[[1]],var=rho[isc])+ImplicitSourceTerm((lambdac[isc]+lambdadeltac[isc]*n[2]/(n[2]+Kml2[isc]))*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                                        
                    elif dynamicmode=="AACarbonGRClog":
                        eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*xhi[isc]*phi[1][isc]*[[1]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*(1.-rho[isc]/CClog[isc]),var=rho[isc])))
                    elif dynamicmode=="AACarbonGRClog2":
                        eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*xhi[isc]*phi[1][isc]*[[1]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*np.power(CClog[isc],hilllog[isc])/(np.power(CClog[isc],hilllog[isc])+np.power(rho[isc],hilllog[isc]))+deltal1[isc]*np.power(CClog1[isc],hilllog[isc])/(np.power(CClog1[isc],hilllog[isc])+np.power(rho[isc],hilllog[isc]))+deltal2[isc]*np.power(Km2[isc],hilllogCT[isc])/(np.power(Km2[isc],hilllogCT[isc])+np.power(n[2],hilllogCT[isc])),var=rho[isc])))
                                        
                    #eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*xhi[isc]*phi[1][isc]*[[1]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                                                             
                                        
                    elif dynamicmode=="AACarbonLOGISTIC":
                        if n_strains==1: 
                            eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*xhi[isc]*phi[1][isc]*[[1]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*(1.-rho[0]/rhomaxlogistic),var=rho[isc])))
                        elif n_strains==2: 
                            eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*xhi[isc]*phi[1][isc]*[[1]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*(1.-(rho[0]+rho[1])/rhomaxlogistic),var=rho[isc])))
                        elif n_strains==5: 
                            eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*xhi[isc]*phi[1][isc]*[[1]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*(1.-(rho[0]+rho[1]+rho[2]+rho[3]+rho[4])/rhomaxlogistic),var=rho[isc])))
                        elif n_strains==11: 
                            eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*xhi[isc]*phi[1][isc]*[[1]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*(1.-(rho[0]+rho[1]+rho[2]+rho[3]+rho[4]+rho[5]+rho[6]+rho[7]+rho[8]+rho[9]+rho[10])/rhomaxlogistic),var=rho[isc])))
                        else:
                            numberstrainsnotknown
                    else:
                        errordynamicmodenotknown
                   
            
        
        
     
    
      
    
    #define equation array by combination of different equations
    if runmode in ['noweber','weber','weber2','weber-precise','weber-log']:
        if n_strains==1:
            if n_nutrients==1:
                eq = eqrho[0] & eqn1 #eqn1 & eqrho1
            elif n_nutrients==2:
                if dynamicmode in ["AACarbonLOGISTIC","AACarbonGRClog","AACarbonGRClog2","AACarbonGRC4"]:
                    eq = eqrho[0] & eqn2
                elif dynamicmode in ["AACarbonGRC3neqa","AACarbonGRC3neqalogcut"]:
                    eq = eqrho[0] & eqn1
                else:
                    eq = eqrho[0] & eqn1 & eqn2
            elif n_nutrients==3:
                
                    eq = eqrho[0] & eqn1 & eqn2 & eqn3
                
        elif n_strains==2:
            if n_nutrients==1:
                eq = eqrho[0] & eqrho[1] & eqn1
            elif n_nutrients==2:
                if dynamicmode=="AACarbonLOGISTIC":
                    eq = eqrho[0] & eqrho[1] & eqn2
                else:
                    eq = eqrho[0] & eqrho[1] & eqn1 & eqn2
            else:
                eror
        elif n_strains==11:
            if n_nutrients==2:
                if dynamicmode=="AACarbonLOGISTIC":
                    eq = eqrho[0] & eqrho[1] & eqrho[2] & eqrho[3] & eqrho[4] & eqrho[5] & eqrho[6] & eqrho[7] & eqrho[8] & eqrho[9] & eqrho[10] & eqn2 
                else:
                    eq = eqrho[0] & eqrho[1] & eqrho[2] & eqrho[3] & eqrho[4] & eqrho[5] & eqrho[6] & eqrho[7] & eqrho[8] & eqrho[9] & eqrho[10] & eqn1 & eqn2 
            else:
                numnutrientsnotdefined
        elif n_strains==3:
            if n_nutrients==2:
                if dynamicmode=="AACarbonLOGISTIC":
                    eq = eqrho[0] & eqrho[1] & eqrho[2] & eqn2 
                else:
                    eq = eqrho[0] & eqrho[1] & eqrho[2] & eqn1 & eqn2 
            else:
                numnutrientsnotdefined
        elif n_strains==7:
            if n_nutrients==2:
                if dynamicmode=="AACarbonLOGISTIC":
                    eq = eqrho[0] & eqrho[1] & eqrho[2] & eqrho[3] & eqrho[4] & eqrho[5] & eqrho[6] & eqn2 
                else:
                    eq = eqrho[0] & eqrho[1] & eqrho[2] & eqrho[3] & eqrho[4] & eqrho[5] & eqrho[6] & eqn1 & eqn2 
            else:
                numnutrientsnotdefined
        
        elif n_strains==5:
            if n_nutrients==2:
                if dynamicmode=="AACarbonLOGISTIC":
                    eq = eqrho[0] & eqrho[1] & eqrho[2] & eqrho[3] & eqrho[4] & eqn2 
                else:
                    eq = eqrho[0] & eqrho[1] & eqrho[2] & eqrho[3] & eqrho[4] & eqn1 & eqn2 
            else:
                numnutrientsnotdefined
        else:
            numstraisnnotdefined
            
    elif runmode in ['fisherwave']:
        eq=eqrho
    elif runmode in ['Keller']:
        eq = eqrho[0] & eqn1
    else:
        "simulation mode not known"
    
    
    
    tsteps=int(timetotal/dt)
    #store over time
    length_n=len(n)
    length_m=len(m)
    length_rho=len(rho)
    
    try:
        every_x=int(p['store_everyx'])
    except:
        every_x=1
    
    #prepare numpy arrays for output
    if p['dimensionmode']=="two":
        nxv=nx*nx
        nxvgrad=nx*nx
    elif p['runmode'] in ["weber","weber2","weber-precise",'weber-log','Keller']:
        nxv=nx
        nxvgrad=nx-1
    else:
        nxv=nx
        nxvgrad=nx
    
    nxvad=int(nx/float(every_x))
    n_output=np.zeros((length_n,vsteps+1,nxvad))
    rho_output=np.zeros((length_rho,vsteps+1,nxvad))
    
    try:
        m_output=np.zeros((length_m,vsteps+1,nxvad))
    except:
        m_output=np.zeros((1,vsteps+1,nxvad))
    for j in range(0,length_n):
        n_output[j,0,:]=n[j][::every_x,]
    
    for j in range(0,length_rho):
        rho_output[j,0,:]=rho[j][::every_x,]
    for j in range(0,length_m):
        m_output[j,0,:]=m[j][::every_x,]
    
        
    print "View "+str(0)+" of "+str(vsteps) +": time "+str(0.)+" of "+str(timetotal)+" at cpu time: "+str(datetime.datetime.now())

        
    #prepare everything for t/space output
    t_output=np.zeros((vsteps+1))
    x_output=np.zeros(nxvad)
    for ils in range(0,nxvad):
        x_output[ils]=p['radius_cut']+dx*(ils+0.5)*every_x
    
    stepssingleview=int(tsteps/vsteps)
    timec=0
    t_output[0]=timec
    
    if  3>2:
        dynamicsnotyetupdated=True
       
        
        
        for vstep in range(1,vsteps+1):
            
            
            for step in range(1,stepssingleview+1):
                
               
                
                #make sure values are up to date            
                for j in range(0,length_n/2):
                    n[j*2].updateOld()
                for j in range(0,length_m):
                    m[j].updateOld()
                #update p variables which only occur in source term and which are not modelred implicitly
                for j in range(0,length_rho):
                    rho[j].updateOld()
                
                for j in range(0,length_n/2):
                    if usefacegrad==False:
                        n[2*j+1][:]=n[2*j].grad.value[0,:]
                    
                    
                   
                    else:
                        
                        n[2*j+1][:]=n[2*j].faceGrad.value[0,1:]
                        n[2*j+1][:]=n[2*j+1][:]*0.5*(np.sign(n[2*j][:]-0.00001)+1.) #make cutoff when concentration falls below very low levels, levels much lower than weber offset (that id handleld via the dif equations)
                                                
                         
                    if runmode =='weber-log':
                        for inn in range(0,n_nutrients):
                            for inns in range(0,n_strains):
                                
                                if dynamicmode in ["AACarbonGRC3neqa","AACarbonGRC3neqalogcut"]:
                                    #phi[inn][inns][:]=weber_offsetlogalpha*np.log(weber_offset[inns]+n[0*inn]).faceGrad.value[0,:-1]
                                    phi[inn][inns][:]=weber_offsetlogalpha*np.log((1.+n[0]/weber_offset[inns])/(1.+n[0]/weber_offsetmax[inns])).faceGrad.value[0,:-1]   
                                    #n[inn+1]=n[inn].faceGrad.value[0,:-1]
                                    phi[inn][inns][0]=0
                    
                                elif dynamicmode=="AACarbonGRC3NL":
                                    phi[inn][inns][:]=weber_offsetlogalpha*np.log((1.+n[2*inn]/weber_offset[inns])/(1.+n[2*inn]/weber_offsetmax[inns])).faceGrad.value[0,:-1]                                                       

                                else:
                                    phi[inn][inns][:]=weber_offsetlogalpha*np.log((1.+n[2*inn]/weber_offset[inns])/(1.+n[2*inn]/weber_offsetmax[inns])).faceGrad.value[0,:-1]                                                       
                                #phi[inn][inns][:]=weber_offsetlogalpha*np.log(weber_offset[inns]+n[2*inn]).faceGrad.value[0,:-1]  
                                                                                   
                                                                
                                #n[inn+1]=n[inn].faceGrad.value[0,:-1]
                                phi[inn][inns][0]=0
                                #phi[inn][inns][0]=0  
                                phi[inn][inns][-1]=np.power(10.,-15) 
             
                                
                    
                    ##here
                    n[2*j+1][0]=0  
                    #n[2*j+1][-2]=0   
                    n[2*j+1][-1]=np.power(10.,-15) 
                    #n[2*j+1][-2]=np.power(10.,-15)  
                 
                 
                eq.solve(dt=dt)
                #update time
                timec=timec+dt
                
            
             
            print "\r View "+str(vstep)+" of "+str(vsteps) +": time "+str(timec)+" of "+str(timetotal)+" at cpu time: "+str(datetime.datetime.now())
            
            if np.nanmin(n[0])<0:
                print "min n"            
                print np.nanmin(n[0])
                
           
            
            #put back again
            for j in range(0,int(length_n/2)):
                n_output[2*j,vstep,:]=n[2*j][::every_x,]
            for j in range(0,int(length_n/2)):
                #n_output[2*j+1,vstep,:]=n[2*j+1][0:-1,]
                n_output[2*j+1,vstep,:]=n[2*j+1][::every_x,]
            for j in range(0,length_rho):
                rho_output[j,vstep,:]=rho[j][::every_x,]
            #print rho_output[0,vstep,0:150]
            for j in range(0,length_m):
                m_output[j,vstep,:]=m[j][::every_x,]
           
            rhos=0
            for irl in range(0,rho_output.shape[2]):
                rhos=rhos+rho_output[0,vstep,irl]

            t_output[vstep]=timec
        
            #switching dynamics 
            if change_xhiB>0 and change_xhiB<=timec and dynamicsnotyetupdated==True:
                print "Swtich dynamics..."
                print timec
                dynamicsnotyetupdated=False
                eqrho=[]
                for isc in range(0,n_strains):   
                    if dynamicmode=="AACarbon":
                        eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*xhi[isc]*phi[1][isc]*[[1]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                                        
                    elif dynamicmode=="AACarbonOD":
                        errornotdefinedforchange
                    elif dynamicmode=="AACarbonOD4":
                        errornotdefinedforchange
                        
                    elif dynamicmode=="AACarbonGRC":
                        eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*xhi[isc]*phi[1][isc]*[[1]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                    elif dynamicmode=="AACarbonGRC3":
                        if change_mode==0:      #all
                               eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*xhi[isc]*phi[1][isc]*[[1]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                               eqn1 = (TransientTerm(var=n[0]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])/Y1[0],var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[0]))
                               eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*(Y2k1[0]*lambdac[0]*n[0]/(n[0]+Km1[0])+Y2k2[0])*n[2]/(n[2]+Km2[0]),var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
                        elif change_mode==1:     #no diffuiosn
                               eqrho.append((TransientTerm(var=rho[isc]) == ExponentialConvectionTerm(coeff=(-1)*xhi[isc]*phi[1][isc]*[[1]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                               eqn1 = (TransientTerm(var=n[0]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])/Y1[0],var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[0]))
                               eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*(Y2k1[0]*lambdac[0]*n[0]/(n[0]+Km1[0])+Y2k2[0])*n[2]/(n[2]+Km2[0]),var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
                        elif change_mode==2:     #no growth
                               eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*xhi[isc]*phi[1][isc]*[[1]],var=rho[isc])))
                               eqn1 = (TransientTerm(var=n[0]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])/Y1[0],var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[0]))
                               eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*(Y2k1[0]*lambdac[0]*n[0]/(n[0]+Km1[0])+Y2k2[0])*n[2]/(n[2]+Km2[0]),var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
                        elif change_mode==3:    #no growth and diffusion
                               eqrho.append((TransientTerm(var=rho[isc]) == ExponentialConvectionTerm(coeff=(-1)*xhi[isc]*phi[1][isc]*[[1]],var=rho[isc])))
                               eqn1 = (TransientTerm(var=n[0]) == ImplicitSourceTerm((-1)*lambdac[0]*n[0]/(n[0]+Km1[0])/Y1[0],var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[0]))
                               eqn2 = (TransientTerm(var=n[2]) == ImplicitSourceTerm((-1)*(Y2k1[0]*lambdac[0]*n[0]/(n[0]+Km1[0])+Y2k2[0])*n[2]/(n[2]+Km2[0]),var=rho[0])+DiffusionTerm(coeff=D_nutrients,var=n[2]))
                    

                        #eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*xhi[isc]*phi[1][isc]*[[1]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                    elif dynamicmode=="AACarbonGRC4":
                        error
                        eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*xhi[isc]*phi[1][isc]*[[1]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                                        
                    elif dynamicmode=="AACarbonGRC2":
                        eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*xhi[isc]*phi[1][isc]*[[1]],var=rho[isc])+ImplicitSourceTerm((lambdac[isc]+lambdadeltac[isc]*n[2]/(n[2]+Kml2[isc]))*n[0]/(n[0]+Km1[isc]),var=rho[isc])))
                                        
                    elif dynamicmode=="AACarbonGRClog":
                        #logistic growth: 
                        eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*xhi[isc]*phi[1][isc]*[[1]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*(1.-rho[isc]/CClog[isc]),var=rho[isc])))
                    elif dynamicsmode=="AACarbonGRClog2":
                        error_notyet_defined
                    elif dynamicmode=="AACarbonLOGISTIC":
                        eqrho.append((TransientTerm(var=rho[isc]) == DiffusionTerm(coeff=D_chemotaxis[isc], var=rho[isc])+ExponentialConvectionTerm(coeff=(-1)*xhi[isc]*phi[1][isc]*[[1]],var=rho[isc])+ImplicitSourceTerm(lambdac[isc]*(1.-(rho[0]+rho[1])/rhomaxlogistic),var=rho[isc])))
                #redifine dynamics....
                if dynamicmode in ["AACarbonLOGISTIC","AACarbonGRC3","AACarbonGRC3logcut"]:
                    eq = eqrho[0] & eqn1 & eqn2
                else:
                    eq = eqrho[0] & eqrho[1] & eqn1 & eqn2
#
    return [x_output,t_output,n_output,rho_output,m_output]
  
#function to merge different parameter scans into one output
def merge_simulation(name_sim, clusterrun=False,tritoncluster=False,peregrinecluster=False,diagonalonly=-1,changetodifference=False,differentstrainnumbers=-1):
    name_sim=name_sim.split(".par")[0]
    name_sim=name_sim.split(".run")[0]
    
    #deine were merged files should be stores
    simpathc=os.path.join('/USERS/jonascremer/Desktop/humangutsim/simulations/testsimulations',name_sim) 
    simpathcout=os.path.join('/USERS/jonascremer/Desktop/humangutsim/simulations/testsimulations_merged',name_sim) 
        
    if not os.path.exists(simpathcout):
            os.makedirs(simpathcout)
    else:
            pass
   
    filetagarray=[""]
    present_repeats=[]
    listfiles=os.listdir(simpathc)
    for il in listfiles:
        if len(il.split("repeat_"))>1:
            if il[-3:]=='npz':
                present_repeats.append(il.split("repeat_")[1].split(".")[0])
    
    present_repeats= set(present_repeats)
    present_repeats= list(present_repeats)
    if len(present_repeats)>0:
        filetagarray=[]
        for il in present_repeats:
            filetagarray.append("_repeat_"+il)
                
    #check for filetags 
    for curfiletag in filetagarray:
        print "current filetage "+curfiletag
        #open parameter one
        firstcounter=-1
        
        for simcounter in range(1,1000):
            
            datafile=os.path.join(simpathc,'data_'+str(simcounter)+curfiletag+'.npz')
                                        
            parfile=os.path.join(simpathc,'par_'+str(simcounter)+curfiletag+'.txt')
            parsweepAfile=os.path.join(simpathc,'parsweep_A_'+str(simcounter)+curfiletag+'.txt')
            parsweepBfile=os.path.join(simpathc,'parsweep_B_'+str(simcounter)+curfiletag+'.txt')
            gooncheck=False
            try:        
                npzfile=np.load(datafile)
                print "Including file: "+datafile
                gooncheck=True
                if firstcounter<0:
                    firstcounter=simcounter
            except:
                if simcounter<10:
                    print "File "+datafile+" not found."
            if gooncheck==True:
                if diagonalonly>0:
                    x=npzfile['x'][::4]
                    n=npzfile['n'][:,:,:,:,::4]
                    rho=npzfile['rho'][:,:,:,:,::4]
                    m=npzfile['m'][:,:,:,:,::4]
                    
                else:
                    n=npzfile['n']
                    rho=npzfile['rho']
                    m=npzfile['m']
                    x=npzfile['x']
                    
                if differentstrainnumbers>0: #makre sure that rho is adjusted when different strain numbers are included....
                    rhocc=np.zeros([rho.shape[0],rho.shape[1],differentstrainnumbers,rho.shape[3],rho.shape[4]])
                    rhocc[:,:,:rho.shape[2],:,:]=rho
                    rho=rhocc
                    
                t=npzfile['t']
                #for first sim, generate new arrays
                if simcounter==firstcounter:
                    t_old=t
                    x_old=x
                    n_old=n
                    rho_old=rho
                    m_old=m
                    #p_old=pHc
                  
                #read in parameter files
                try:
                    json_data=open(parfile)
                    dictpar = json.load(json_data)
                    json_data.close()
                except:
                    print "Parameter file "+parfile+" not found"
                    #implement here:
                    
                #for first sim, store parfile as new one, for following, check if parfile is the same
                if simcounter==firstcounter:
                    dictpar_old=dictpar
                    #store parameter file as new one
                else:
                    shared_items = set(dictpar_old.items()) & set(dictpar.items())
                    
                    if len(shared_items)<len(shared_items):
                        print "shared items"
                        print shared_items
                        print "not same parameter file"
                        print "Break at simulation"+str(isimc)
                        break
                    
                #open parsweep files
                try:
                    json_data=open(parsweepAfile)
                    dictparA = json.load(json_data)
                    json_data.close()
                except:
                    print "Parameter file "+parsweepAfile+" not found"
                    break
                #try:
                if 3>2:
                    print parsweepBfile                                        
                    json_data=open(parsweepBfile)
                    dictparB = json.load(json_data)
                    if changetodifference:
                        dictparB[dictparB.keys()[0]]=float(round((dictparB[dictparB.keys()[0]]-dictparA[dictparA.keys()[0]])*1000000000000.,2))
                        print dictparB[dictparB.keys()[0]]
                        
                        
                    json_data.close()
                  
                if simcounter==firstcounter:
                    dictparB_old=dictparB
                    dictparA_old=dictparA
                    print dictparA_old
                    print dictparA_old.keys()[0]
                    keyA_new=dictparA_old.keys()[0].split(':')[0]
                    keyB_new=dictparB_old.keys()[0].split(':')[0]
                if simcounter>firstcounter:
                    #t stays same t_new=t
                    #x stays same x_new=x
                    keyA=dictparA.keys()[0].split(':')[0]
                    keyB=dictparB.keys()[0].split(':')[0]
                    
                    if ((keyA != keyA_new) or (keyB != keyB_new) ):
                        print "Error: varied parameteres are different"
                        return False
                    
                    #find unique parameter values
                    valueA_unique=list(set(dictparA.values()+dictparA_old.values()))
                    valueB_unique=list(set(dictparB.values()+dictparB_old.values()))
                    
                    #generate new numpy arrays with right dimension
                    lx=len(valueA_unique)
                    ly=len(valueB_unique)
                    
                    try:
                        valueA_uniquenp=np.array([float(xlc) for xlc in valueA_unique])
                    except:
                        valueA_uniquenp=np.array([xlc for xlc in valueA_unique])
                    try:
                        valueB_uniquenp=np.array([float(xlc) for xlc in valueB_unique])
                    except:
                        valueB_uniquenp=np.array([xlc for xlc in valueB_unique])
                    valueA_uniquenp.sort()
                    valueB_uniquenp.sort()
                    valueA_unique=valueA_uniquenp.tolist()
                    valueB_unique=valueB_uniquenp.tolist()
                    
                    maxxdim=4000
                    n_new=np.zeros([lx,ly,n_old.shape[2],n_old.shape[3],maxxdim])
                    m_new=np.zeros([lx,ly,m_old.shape[2],m_old.shape[3],maxxdim])
                    rho_new=np.zeros([lx,ly,rho_old.shape[2],rho_old.shape[3],maxxdim])
                    
                    n_new[:]=np.nan
                    rho_new[:]=np.nan
                    m_new[:]=np.nan
                    
                    dictparA_new={}
                    dictparB_new={}
                    
                    counterA=-1
                    #go through different values and take old values to new
                    
                    for valueA in valueA_unique:
                        counterA=counterA+1
                        dictparA_new[keyA+":"+str(counterA)]=valueA
                        
                        counterB=-1
                        for valueB in valueB_unique:
                            
                            
                            counterB=counterB+1
                            dictparB_new[keyB+":"+str(counterB)]=valueB
                            #check if values are present in simulation to add      
                            valueexist=False
                            valueexistA=False
                            for varc in dictparA.values():
                                try:
                                    if (varc/valueA< 1.000000001 and varc/valueA>.999999999) or varc==valueA:
                                        valueexistA=True
                                        varcA=varc
                                except:
                                    
                                    
                                    
                                    try:
                                     if (varc-valueA<.00000000001 and valueA-varc<.0000000001):
                                            valueexistA=True
                                            varcA=varc   
                                    except:
                                        if (varc==valueA):
                                            valueexistA=True
                                            varcA=varc  
                            if valueexistA:
                                rankA=dictparA.keys()[dictparA.values().index(varcA)].split(":")[1]
                                valueexistB=False
                                
                                for varc in dictparB.values():
                                    
                                        
                                        try:
                                            if (varc/valueB< 1.0000001 and varc/valueB>.9999999999) or varc==valueB:
                                                valueexistB=True
                                                varcB=varc
                                        except:
                                            try:
                                                if varc-valueB<.000000001 and valueB-varc<.000000001:
                                                    valueexistB=True
                                                    varcB=varc
                                            except:
                                                if (varc==valueB):
                                                    valueexistB=True
                                                    varcB=varc
                          
                                if valueexistB:
                                   
                                    valueexist=True
                                    #print varc
                                    #print dictparB.keys()[dictparB.values().index(varcB)]
                                    rankB=dictparB.keys()[dictparB.values().index(varcB)].split(":")[1]
                                    rankB=int(rankB)
                                    rankA=int(rankA)
                                    n_dimc=n.shape[4]
                                    
                                    
                                    #adjusted to allow different sizes...
                                    n_new[counterA,counterB,:,:,:n_dimc]=n[rankA,rankB,:,:,:]
                                    m_new[counterA,counterB,:,:,:n_dimc]=m[rankA,rankB,:,:,:]
                                    rho_new[counterA,counterB,:,:,:n_dimc]=rho[rankA,rankB,:,:,:]
                                    
                                  
                            if valueexist==True:
                                pass
                            else:
                                valueexistA=False
                                for varc in dictparA_old.values():
                                    
                                    try:
                                        if (varc==valueA) or (varc/valueA< 1.0001 and varc/valueA>.9999):
                                            valueexistA=True
                                            varcA=varc
                                    except:
                                        try:
                                            if varc-valueA<.00001 and valueA-varc<.0001:
                                                valueexistA=True
                                                varcA=varc
                                        except: #to handle stringa
                                            if varc==valueA:
                                                valueexistA=True
                                                varcA=varc
                                                
                                if  valueexistA:
                                    rankA=dictparA_old.keys()[dictparA_old.values().index(varcA)].split(":")[1]
                                    #check if list exists in old directory
                                    valueexistB=False
                                    for varc in dictparB_old.values():
                                        
                                             try:
                                                if (varc==valueB) or varc/valueB< 1.0001 and varc/valueB>.9999:
                                                    valueexistB=True
                                                    varcB=varc
                                             except:
                                                try:
                                                    if varc-valueB<.00001 and valueB-varc<.0001:
                                                        valueexistB=True
                                                        varcB=varc
                                                except:
                                                    pass #to make sure things are working if parameter value is a string
                                                    
                                            
                                    if valueexistB:
                                        rankB=dictparB_old.keys()[dictparB_old.values().index(varcB)].split(":")[1]
                                        rankA=int(rankA)
                                        rankB=int(rankB)
                                        n_dimcc=n_old.shape[4]
                                        n_new[counterA,counterB,:,:,:n_dimcc]=n_old[rankA,rankB,:,:,:]
                                        m_new[counterA,counterB,:,:,:n_dimcc]=m_old[rankA,rankB,:,:,:]
                                        rho_new[counterA,counterB,:,:,:n_dimcc]=rho_old[rankA,rankB,:,:,:]
                                        
                                        valuefound=True
                                    else:
                                        pass
                                    
                            
                    dictparA_old=dictparA_new
                    dictparB_old=dictparB_new
                    n_old=n_new
                  
                    m_old=m_new
                    rho_old=rho_new
                    
        
        if diagonalonly<10:
            print "Check for missing parameters in grid..."
            #after last simulation added, show which values are missing in grid
            for ixl in range(0,n_old.shape[0]):
                for iyl in range(0,n_old.shape[1]):
                    if np.isnan(n_old[ixl,iyl,0,0,0]):
                        keysoutA=dictparA_old.keys()
                        keysoutB=dictparB_old.keys()
                        for ikeylA in range(0,len(keysoutA)):
                            for ikeylB in range(0,len(keysoutB)):
                                if (int(keysoutA[ikeylA].split(':')[1])==ixl) and (int(keysoutB[ikeylB].split(':')[1])==iyl):
                                   print "Simulation still missing for parameters A: "+str(dictparA_old[keysoutA[ikeylA]])+" B: "+str(dictparB_old[keysoutB[ikeylB]])
        
        #save files 
        np.savez(os.path.join(simpathcout,'data'+curfiletag+'.npz'),n=n_old, rho=rho_old,m=m_old, x=x_old,t=t_old)
                 
        file_ = open(os.path.join(simpathcout,'par'+curfiletag+'.txt'), 'w')
        file_.write(json.dumps(dictpar_old))
        file_.close()  
        
        
        file_ = open(os.path.join(simpathcout,'parsweep_A'+curfiletag+'.txt'), 'w')
        file_.write(json.dumps(dictparA_old))
        file_.close()    
        
        if changetodifference:
            for ik in dictparB_old.keys():
                dictparB_old[ik]=dictparB_old[ik]/1000000000000.
        file_ = open(os.path.join(simpathcout,'parsweep_B'+curfiletag+'.txt'), 'w')
        file_.write(json.dumps(dictparB_old))
        file_.close()
     
    
    return True

def merge_repeatsimulation(name_sim, clusterrun=False,tritoncluster=False,peregrinecluster=False):
    name_sim=name_sim.split(".par")[0]
    name_sim=name_sim.split(".run")[0]
    if full_prog=="full":
        simpathcout=os.path.join('/Users/jonascremer/Documents/Science/chemotaxis',name_sim) 
        simpathcouttime=os.path.join('/Users/jonascremer/Documents/Science/chemotaxis',name_sim+"_TM")
        
        if not os.path.exists(simpathcout):
            os.makedirs(simpathcout)
        else:
            pass
        
        if not os.path.exists(simpathcouttime):
            os.makedirs(simpathcouttime)
        else:
            pass

    
    filetagarray=[""]
    present_repeats=[]
    listfiles=os.listdir(simpathcout)
    print "found files in folder"
    print simpathcout
    print listfiles
    
    for il in listfiles:
        if len(il.split("repeat_"))>1:
            if il[-3:]=='npz':
                present_repeats.append(il.split("repeat_")[1].split(".")[0])
    
    present_repeats= set(present_repeats)
    present_repeats= list(present_repeats)
    print "Repeats presented:"
    print present_repeats
    present_repeats.sort()
    print present_repeats
    
    #unique repeats
    if len(present_repeats)>0:
        filetagarray=[]
        for il in present_repeats:
            filetagarray.append("_repeat_"+il)
    repc=[]
    for ilv in present_repeats:
        repc.append(int(ilv.split('_')[-1])+1)
        #old int(present_repeats[-1].split('_')[-1])+1
    if len(filetagarray) != max(repc):
        print "error, numbers of repats wrong, something is missing"
        return False
    
    #go through repeats and determine number of time steps for total arrays
    print "go through files to determine total time"
    timestepsc=1
    for curfiletag in filetagarray:
        datafile=os.path.join(simpathcout,'data'+curfiletag+'.npz')
        print datafile
        npzfile=np.load(datafile)
        timestepsc=timestepsc+npzfile["t"].shape[0]-1
        #determine other dimensions
    
    curfiletag =  filetagarray[0]  #use first tile here.... 
    datafile=os.path.join(simpathcout,'data'+curfiletag+'.npz')  
    npzfile_shape=np.load(datafile)['rho'].shape
    numparA= npzfile_shape[0]  
    numparB= npzfile_shape[1]
    numX=npzfile_shape[-1]
    #that should be correct dimension...
    
    #generate empty arrays to fill in timesteps
    
    
    firstfilec=True
    timecountcT=0
    timebefore=0
    
    firsttag=True
    for curfiletag in filetagarray: 
        datafile=os.path.join(simpathcout,'data'+curfiletag+'.npz')
        
        npzfile=np.load(datafile)
        
        #load parameter file
        parfile=os.path.join(simpathcout,'par'+curfiletag+'.txt')
        parfileAsweep=os.path.join(simpathcout,'parsweep_A'+curfiletag+'.txt')
        parfileBsweep=os.path.join(simpathcout,'parsweep_B'+curfiletag+'.txt')
        json_data=open(parfile)
        dictparcc = json.load(json_data)
        json_data.close()
        json_data=open(parfileAsweep)
        dictparccA = json.load(json_data)
        json_data.close()
        json_data=open(parfileBsweep)
        dictparccB = json.load(json_data)
        json_data.close()
        avaluec=[]
        for aic in dictparccA:
            avaluec.append(dictparccA[aic])
        bvaluec=[]
        for bic in dictparccB:
            bvaluec.append(dictparccB[bic])
        bvaluec.sort()
        avaluec.sort()
        if firsttag:
            #get all a and b values present
            bvaluecfirst=bvaluec
            avaluecfirst=avaluec
            
            if dictparcc['evolution_mode']==True:
                rhocT=np.zeros([numparA,numparB,dictparcc['evolution_strains'],timestepsc,numX])  
                ncT=np.zeros([numparA,numparB,4,timestepsc,numX])  
                mcT=np.zeros([numparA,numparB,0,timestepsc,numX])
            else:
                rhocT=np.zeros([numparA,numparB,2,timestepsc,numX])  
                ncT=np.zeros([numparA,numparB,2,timestepsc,numX])  
                if dictparcc['runmode']=='onespecies':
                    mcT=np.zeros([numparA,numparB,0,timestepsc,numX])  
                else:
                    mcT=np.zeros([numparA,numparB,4,timestepsc,numX])  
            #if dictparcc['runmode']=='onespecies':
            #    pcT=np.zeros([numparA,numparB,0,timestepsc,numX])  
            #else:
            #    pcT=np.zeros([numparA,numparB,2,timestepsc,numX])  
            xcT=np.zeros([numX])    
            timecT=np.zeros([timestepsc])
            rhocT[:]=np.nan
            ncT[:]=np.nan
            mcT[:]=np.nan
            #pcT[:]=np.nan
            xcT[:]=np.nan
            timecT[:]=np.nan
            
            
    
            firsttag=False
        endtime=dictparcc["time"]
        
        timestepsc=npzfile["t"].shape[0]
        if firstfilec==True:
            
            firstfilec=False
            xcT=npzfile['x']
            timecT[timecountcT:timecountcT+timestepsc]=npzfile['t'][:]+timebefore #adjust add time
            rhocT[:,:,:,timecountcT:timecountcT+timestepsc,:]=npzfile['rho'][:,:,:,:,:]
            mcT[:,:,:,timecountcT:timecountcT+timestepsc,:]=npzfile['m'][:,:,:,:,:]
            ncT[:,:,:,timecountcT:timecountcT+timestepsc,:]=npzfile['n'][:,:,:,:,:]
            #pcT[:,:,:,timecountcT:timecountcT+timestepsc,:]=npzfile['p'][:,:,:,:,:]
            timecountcT = timecountcT+timestepsc
            timebefore=npzfile['t'][-1]
            #print "rho 0, 0 par comb, first timestep"
            #print npzfile['rho'][1,0,0,-1,:]
            #print "other initial states"
            #print npzfile['rho'][0,0,0,-1,:]
            #print npzfile['rho'][0,1,-1,:]
            #print npzfile['rho'][1,1,0,-1,:]
            #for ir in range(0,npzfile['t'].shape[0]):
            #    print str(ir)+": "+str(npzfile['t'][ir])
        else:
            #print "notfirst"
            #print curfiletag
            #length to add: 97
            timestepsc=timestepsc-1
            #print timecountcT+timestepsc
            #print timecountcT
            #print timecT.shape
            try:
                timecT[timecountcT:timecountcT+timestepsc]=npzfile['t'][1:]+timebefore #adjust add time
                rhocT[:,:,:,timecountcT:timecountcT+timestepsc,:]=npzfile['rho'][:,:,:,1:,:]
                mcT[:,:,:,timecountcT:timecountcT+timestepsc,:]=npzfile['m'][:,:,:,1:,:]
                ncT[:,:,:,timecountcT:timecountcT+timestepsc,:]=npzfile['n'][:,:,:,1:,:]
                #pcT[:,:,:,timecountcT:timecountcT+timestepsc,:]=npzfile['p'][:,:,:,1:,:]
            except:
                print "Some parameter combinations are missing "+datafile
                bindxc=[]
                for ilcb in range(0,len(bvaluecfirst)):
                    
                    if bvaluecfirst[ilcb] in bvaluec:
                        bindxc.append(ilcb)
                aindxc=[]
                for ilca in range(0,len(avaluecfirst)):
                    if avaluecfirst[ilca] in avaluec:
                        aindxc.append(ilca)
                if len(avaluec) != len(avaluecfirst):
                    
                    print str(len(avaluecfirst)-len(avaluec)) +" parameter scans for A still missing"
                    print "Have: "+str(avaluec)
                    print "Should: "+str(avaluecfirst)
                if len(bvaluec) != len(bvaluecfirst):
                    print str(len(bvaluecfirst)-len(bvaluec)) +" parameter scans for B still missing"
                    print "Have: "+str(bvaluec)
                    print "Should: "+str(bvaluecfirst)
                
                if len(bindxc)==1:
                    timecT[timecountcT:timecountcT+timestepsc]=npzfile['t'][1:]+timebefore #adjust add time
                    rhocT[aindxc,bindxc,:,timecountcT:timecountcT+timestepsc,:]=npzfile['rho'][:,0,:,1:,:]
                    mcT[aindxc,bindxc,:,timecountcT:timecountcT+timestepsc,:]=npzfile['m'][:,0,:,1:,:]
                    ncT[aindxc,bindxc,:,timecountcT:timecountcT+timestepsc,:]=npzfile['n'][:,0,:,1:,:]
                    
                else:
                    timecT[timecountcT:timecountcT+timestepsc]=npzfile['t'][1:]+timebefore #adjust add time
                    rhocT[aindxc,bindxc,:,timecountcT:timecountcT+timestepsc,:]=npzfile['rho'][:,:,:,1:,:]
                    mcT[aindxc,bindxc,:,timecountcT:timecountcT+timestepsc,:]=npzfile['m'][:,:,:,1:,:]
                    ncT[aindxc,bindxc,:,timecountcT:timecountcT+timestepsc,:]=npzfile['n'][:,:,:,1:,:]
                    
            
            #adjust timestep
            timecountcT = timecountcT+timestepsc
            timebefore=timebefore+npzfile['t'][-1]
            #compare first with timestep before
            
            try:
                diff_repeatrho=npzfile['rho'][:,:,:,0,:]-laststeprho
                diff_repeatm=npzfile['m'][:,:,:,0,:]-laststepm
                diff_repeatn=npzfile['n'][:,:,:,0,:]-laststepn
                print "Differences found for the following parameters (# of differences in rho, m, and n)..."
                
                for ila in range(0,diff_repeatrho.shape[0]):
                    for ilb in range(0,diff_repeatrho.shape[1]):
                    
                        nerr=np.count_nonzero(diff_repeatrho[ila,ilb,:,:])
                        if nerr>0:
                            print str(nerr)+" for parameter ranks "+str(ila)+" "+str(ilb)
                print "End checking."
            except:
                print "Difference checking not possible because of differences in dimensions (parameter scans still incomplete)"
            
        laststeprho=npzfile['rho'][:,:,:,-1,:]
        laststepm=npzfile['m'][:,:,:,-1,:]
        laststepn=npzfile['n'][:,:,:,-1,:]
        
       
        
    curfiletag =  filetagarray[0] 
    #load parameter files and save as general one  
    parfile=os.path.join(simpathcout,'par'+curfiletag+'.txt')
    json_data=open(parfile)
    dictparcc = json.load(json_data)
    json_data.close()
    #adjust endtime and timesteps
    dictparcc["time"]=endtime
    #save it...
    file_ = open(os.path.join(simpathcouttime,'par'+'.txt'), 'w')
    file_.write(json.dumps(dictparcc))
    file_.close()  
    
    parsweepAfile=os.path.join(simpathcout,'parsweep_A'+curfiletag+'.txt')
    json_data=open(parsweepAfile)
    dictparcc = json.load(json_data)
    json_data.close()
    #save it
    file_ = open(os.path.join(simpathcouttime,'parsweep_A'+'.txt'), 'w')
    file_.write(json.dumps(dictparcc))
    file_.close()    
        
    
    parsweepBfile=os.path.join(simpathcout,'parsweep_B'+curfiletag+'.txt')
    json_data=open(parsweepBfile)
    dictparcc = json.load(json_data)
    json_data.close()
    #save it
    file_ = open(os.path.join(simpathcouttime,'parsweep_B'+'.txt'), 'w')
    file_.write(json.dumps(dictparcc))
    file_.close()

    #save main data    
    print "mereged main data saved to :"+simpathcouttime
    np.savez(os.path.join(simpathcouttime,'data'+'.npz'),n=ncT, rho=rhocT,m=mcT, x=xcT,t=timecT)
     
        
    
    return True    

#run parameter scans
def run_sweep(name_sim,parameters,x_values,xpar,y_values,ypar, scanparameter=False,nametag=-1):
    timestart=time.time()
    #option scanparameter if yes, data is not overwritten, but new files are generated
    #start simulation with parameter sweep
    simpathc=os.path.join(os.getcwd(),name_sim) 
    
    #generate folder to store results
    if not os.path.exists(simpathc):
        try:
            os.makedirs(simpathc)
        except:
            pass
    else:
        print "ERROR: FOLDER FOR SIMULATION EXISTS ALREADY"
        
    
    simcounter=0
    parametersweep_1={}
    parametersweep_2={}

    if (x_values.shape[0]>=1) and (xpar in parameters) and (np.isnan(x_values.shape[0])==False):
        dox=True
        print "go through x"
    else:
        if (xpar in parameters)==False:
            print "Error: parameter to scan not found"
        #    error
        dox=False
    
        
    if (y_values.shape[0]>=1) and (ypar in parameters) and (np.isnan(y_values.shape[0])==False) and dox==True:
        doy=True
        print "go through y"
    else:
        if (ypar in parameters)==False:
            print "Error: parameter to scan not found"
        #    error
        doy=False
              

    if dox==False and doy==False:              
        [x,t,n,rho,m]=run_simulation(parameters) #,p
        
    
    for isimx in range(0,len(x_values)):
        #sweep through parameter
        if dox:        
            parameters[xpar]=x_values[isimx]
            parametersweep_1[xpar+":"+str(isimx)]=x_values[isimx]
        if doy==False:        
            [x,t,n,rho,m]=run_simulation(parameters)
            if simcounter==0:
                nout=np.zeros([len(x_values),1,n.shape[0],n.shape[1],n.shape[2]])                    
                rhoout=np.zeros([len(x_values),1,rho.shape[0],rho.shape[1],rho.shape[2]])                    
                mout=np.zeros([len(x_values),1,m.shape[0],m.shape[1],m.shape[2]]) 
                #pout=np.zeros([len(x_values),1,p.shape[0],p.shape[1],p.shape[2]])                   
                
            nout[isimx,0,:,:,:]=n
            rhoout[isimx,0,:,:,:]=rho
            mout[isimx,0,:,:,:]=m
            simcounter=simcounter+1
                
        for isimy in range(0,len(y_values)):
            if doy:            
                
                parameters[ypar]=y_values[isimy]
                parametersweep_2[ypar+":"+str(isimy)]=y_values[isimy]
                
                
                
                if ('evolution_mode' in parameters) and parameters['evolution_mode']==True:
                    print "load simulation data"
                    
                    if nametag==None:
                        countcc=isimx*len(x_values)+isimy+1
                    else:
                        countcc=nametag

                    #to set here, where is data of previous simulations stored
                    dd=simpathc=os.path.join(os.getcwd(),name_sim)
                    dirdata=os.path.join(dd,name_sim)
                    print "dirdata to load..."
                    print dirdata
                    timerepeatcc=0
                    iltimerepend=-1
                    iltimerep=100
                    while iltimerep <200 and iltimerep>=0:
                        filec="data_"+str(countcc)+'_repeat_'+str(iltimerep)+".npz"
                        
                        skippfile=True
                        try:
                            #print os.path.join(dirdata,filec)
                            curfile=np.load(os.path.join(dirdata,filec))
                            iltimerepend=iltimerep
                            fileusec=os.path.join(dirdata,filec)
                            skippfile=False                            
                            print "found previous simulation"
                            print os.path.join(dirdata,filec)
                            iltimerep=1000
                            
                        except:
                            iltimerepend=iltimerep
                        iltimerep=iltimerep-1
                     
                    if skippfile==False:
                            curfile=np.load(fileusec)
                            print "*******"
                            print "file found and used"
                            print curfile.files
                            timerepeatcc=iltimerepend+1
                            #one of the following arrays ['m', 'n', 'p', 't', 'rho', 'x']
                            nc=curfile['n']
                            mc=curfile['m']
                            rhoc=curfile['rho']
                            xc=curfile['x']
                            #parameter from run before
                            json_dataoldpar=open(os.path.join(dirdata,"par_"+str(countcc)+'_repeat_'+str(iltimerepend)+".txt"))
                            paroldrun = json.load(json_dataoldpar)
                            #evolutionstepc=paroldrun["evolution_step"]
                            try:
                                evolutionstepc=paroldrun["evolution_step"]
                            except:
                                evolutionstepc=0
                            datastart=[xc,rhoc,nc,mc,iltimerepend]
                            
                            #here89
                            [x,t,n,rho,m]=run_simulation(parameters,loaddata=True,datastart=datastart)
                            
                    elif skippfile==True:
                            [x,t,n,rho,m]=run_simulation(parameters)
                else:
                    
                    [x,t,n,rho,m]=run_simulation(parameters)
                    
                
             
                #start simulation
                simcounter=simcounter+1
        
            
                #store to output file
                #if file is new, generate array            
                if isimy==0 and isimx==0:
                     nout=np.zeros([len(x_values),len(y_values),n.shape[0],n.shape[1],n.shape[2]])                    
                     rhoout=np.zeros([len(x_values),len(y_values),rho.shape[0],rho.shape[1],rho.shape[2]])                    
                     mout=np.zeros([len(x_values),len(y_values),m.shape[0],m.shape[1],m.shape[2]]) 
                nout[isimx,isimy,:,:,:]=n
                rhoout[isimx,isimy,:,:,:]=rho
                mout[isimx,isimy,:,:,:]=m
                
                    
    #store files if no x and no y parameter
    
    if simcounter==0:
        nout=np.zeros([1,1,n.shape[0],n.shape[1],n.shape[2]])                    
        rhoout=np.zeros([1,1,rho.shape[0],rho.shape[1],rho.shape[2]])                    
        mout=np.zeros([1,1,m.shape[0],m.shape[1],m.shape[2]])    
        nout[0,0,:,:,:]=n
        rhoout[0,0,:,:,:]=rho
        mout[0,0,:,:,:]=m
        
        nout=np.expand_dims(n, axis=0)
        rhoout=np.expand_dims(rho,axis=0)
        mout=np.expand_dims(n, axis=0)
        
    if scanparameter==False:
        #save data to file
    
        np.savez(os.path.join(simpathc,'data.npz'),n=nout, rho=rhoout,m=mout, x=x,t=t)#,p=pout
        nametag=""        
    else:
        
        if nametag == 'None':
            #determine name for new simulation
            counter_sim=0
            for icsim in range(0,100):
                counter_sim=counter_sim+1
                
                nametag=str(counter_sim)
                name_dataout=os.path.join(simpathc,'data_'+nametag+'.npz')
                if os.path.exists(name_dataout):
                    pass
                else:
                    break

        else:
            try:
                if parameters["evolution_mode"]==True:
                    nametag=nametag+'_repeat_'+str(timerepeatcc)
                    print "nametag not none"
                    print nametag
            except:
                pass
 
        
        
        if nametag == 'None':
            #determine name for new simulation
            counter_sim=0
            for icsim in range(0,100):
                counter_sim=counter_sim+1
                
                nametag=str(counter_sim)
                if parameters["loaddata"]==True and parameters["runmode"] in ['onespecies','twospecies_crossfeeding','onespecies_buffer','twospecies_buffer','twospecies_buffer_simple','twospecies_buffer_twonutrients']:
                    nametag=nametag+'_repeat_'+str(timerepeatcc)
                    print "nametag none"
                    print nametag
                    
               
             
                name_dataout=os.path.join(simpathc,'data_'+nametag+'.npz')
                if os.path.exists(name_dataout):
                    pass
                else:
                    break
   
        name_dataout=os.path.join(simpathc,'data_'+nametag+'.npz')
        print "merged saving to:"
        print name_dataout
        np.savez(name_dataout,n=nout, rho=rhoout,m=mout, x=x,t=t) #p=pout,
            



    #save parameters
    if scanparameter==False:    
        file_ = open(os.path.join(simpathc,'par.txt'), 'w')
    else:
        file_ = open(os.path.join(simpathc,'par_'+nametag+'.txt'), 'w')
        
    file_.write(json.dumps(parameters))
    file_.close()
    #save information about parameter sweep
    if dox:
        if scanparameter==False: 
            file_ = open(os.path.join(simpathc,'parsweep_A.txt'), 'w')
        else:
            file_ = open(os.path.join(simpathc,'parsweep_A_'+nametag+'.txt'), 'w')
        file_.write(json.dumps(parametersweep_1))
        file_.close()
    if doy:
        if scanparameter==False: 
            file_ = open(os.path.join(simpathc,'parsweep_B.txt'), 'w')
        else:
            file_ = open(os.path.join(simpathc,'parsweep_B_'+nametag+'.txt'), 'w')
        file_.write(json.dumps(parametersweep_2))
        file_.close()
        
        
    timeend=time.time()
    print("Simulation time: "+str(int((timeend-timestart)/60.))+"min")
    
    
    return(x,t,nout,rhoout,mout)
    
def load_simulation(name_sim):
    try:
        namefile=os.path.join(name_sim,'data.npz')
        
        npzfile=np.load(namefile)
    except:
        namefile=os.path.join(name_sim,'data.npz')
        
        npzfile=np.load(namefile)
        
        print "File "+namefile+" not found argh"
        return [False,0,0,0,0,0,0,{},{},{}]
    
    t=npzfile['t']
    x=npzfile['x']
    
    #try:  
    if 3>2:
        n=npzfile['n']
        rho=npzfile['rho']
        m=npzfile['m']
        #p=npzfile['p']
       
        
    #read in parameter files
    try:
        namefile=os.path.join(name_sim,'par.txt')
        json_data=open(namefile)
        dictpar = json.load(json_data)
        json_data.close()
    except:
        print "Parameter file "+namefile+" not found"
        dictpar={}
    #try reading files for parameter-sweeps
    try:
        namefile=os.path.join(name_sim,'parsweep_A.txt')
        json_data=open(namefile)
        pardictx = json.load(json_data)
        json_data.close()
    except:
        try:#check also old file naming.
            namefile=os.path.join(name_sim,'parsweep_1.txt')
            json_data=open(namefile)
            pardictx = json.load(json_data)
            json_data.close()
        except:
            
            pardictx={}
    try:
        namefile=os.path.join(name_sim,'parsweep_B.txt')
        json_data=open(namefile)
        pardicty = json.load(json_data)
        json_data.close()
    except:
        try:#check alo old naming of parameter file
            namefile=os.path.join(name_sim,'parsweep_2.txt')
            json_data=open(namefile)
            pardicty = json.load(json_data)
            json_data.close()
        except:
            pardicty={}
           

    return [True,x,t,n,rho,m,dictpar,pardictx,pardicty]#,p



