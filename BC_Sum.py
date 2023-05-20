import numpy as np
import sys
import math

def BC_Sum(fun,BC,**List):
    
    if BC == 'PBC' or BC == 'aPBC':
        
        #########    Check parameter integrity and self consistency
           ######    Check L or Lm exist
        if (List.get('L') is None) == False:
            if List.get('L') < 0:
                print("ERROR: System size L < 0 is impossible. If you have imput Lm thus this L is unnecessary. ")
                sys.exit(1)
            if (List.get('Lm') is None) == False:
                print("WARNING: Too much parameter are given! It is 'Co-exist of L and Lm' that may cause conflictions. ")
                sys.exit(1)
        else:
            if List.get('Lm') is None:
                print("ERROR: Parameter are NOT given! It is 'Losing L or Lm' which is necessary for PBC/aPBC. ")
                sys.exit(1)
            if List.get('Lm') < 0:
                print("ERROR: System size Lm < 0 is impossible.")
                sys.exit(1)
           ######    Check p2_Limit or p2_error exist
        if (List.get('p2_Limit') is None):
            print("ERROR: Parameter are NOT given! It is 'Losing p2_Limit' which is necessary. ")
            sys.exit(1)
        else:
            if List.get('p2_Limit') < 0:
                print("ERROR: Upper limit of summed momentum p2_Limit < 0 is impossible.")
                sys.exit(1)
                
        
        #########    Pre-Caculation
           ######    Caculate Lm
        if (List.get('L') is None) == False:
            Lm = List.get('L')/197.32
        else:# if (List.get('Lm') is None) == False
            Lm = List.get('Lm')
        
           ######    define SFun
        if BC == 'PBC':
            def SFun(nx,ny,nz):
                return fun((((2*nx)*np.pi/Lm)**2 + ((2*ny)*np.pi/Lm)**2 + ((2*nz)*np.pi/Lm)**2)**0.5)
        else:
            def SFun(nx,ny,nz):
                return fun((((2*nx+1)*np.pi/Lm)**2 + ((2*ny+1)*np.pi/Lm)**2 + ((2*nz+1)*np.pi/Lm)**2)**0.5)
            
        #########    Sum
        result = 0 #Sum pool
        if BC == 'PBC':
            evs = List.get('p2_Limit')
            nl = int((evs*Lm/np.pi)/2)
            ri = 0
            result += SFun(ri,ri,ri)
            nkl = evs**2 - 2*((2*ri)*np.pi/Lm)**2
            if nkl >= 0:
                nkl = int(((nkl)**0.5*Lm/np.pi)/2) 
            
                for rk in range(ri+1,nkl+1,1):
                    njl = (evs**2 - ((2*ri)*np.pi/Lm)**2 - ((2*rk)*np.pi/Lm)**2)
                    if njl >= 0:
                        njl = int(((njl)**0.5*Lm/np.pi)/2)
                        result += 6*SFun(ri,ri,rk)
                    
                        for rj in range(ri+1,njl+1,1):
                            
                            result += 12*SFun(ri,rj,rk)
            if nl>=1:
                for ri in range(1,nl+1):
                    if 3*((2*ri)*np.pi/Lm)**2 <= evs**2:
                        result += 8*SFun(ri,ri,ri)
                        nkl = evs**2 - 2*((2*ri)*np.pi/Lm)**2
                        if nkl >= 0:
                            nkl = int(((nkl)**0.5*Lm/np.pi)/2) 
                        
                            for rk in range(ri+1,nkl+1,1):
                                njl = (evs**2 - ((2*ri)*np.pi/Lm)**2 - ((2*rk)*np.pi/Lm)**2)
                                if njl >= 0:
                                    njl = int(((njl)**0.5*Lm/np.pi)/2)
                                
                                    for rj in range(ri,njl+1,1):
                                        result += 24*SFun(ri,rj,rk)
            result = result/Lm**3
        else:#if BC == 'aPBC':
            evs = List.get('p2_Limit')
            nl = int((evs*Lm/np.pi-1)/2)
            for ri in range(nl+1):
                if 3*((2*ri+1)*np.pi/Lm)**2 <= evs**2:
                    result += SFun(ri,ri,ri)
                    nkl = evs**2 - 2*((2*ri+1)*np.pi/Lm)**2
                    if nkl >= 0:
                        nkl = int(((nkl)**0.5*Lm/np.pi-1)/2) 
                    
                        for rk in range(ri+1,nkl+1,1):
                            njl = (evs**2 - ((2*ri+1)*np.pi/Lm)**2 - ((2*rk+1)*np.pi/Lm)**2)
                            if njl >= 0:
                                njl = int(((njl)**0.5*Lm/np.pi-1)/2)
                            
                                for rj in range(ri,njl+1,1):
                                    result += 3*SFun(ri,rj,rk)
            result = 8*result/Lm**3
        
    return result
