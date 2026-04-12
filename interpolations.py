from ngsolve import VOL
def getVolumeOneNode2d(psit, psif1, psif2):
    a=psit/(psit-psif1)
    b=psit/(psit-psif2)
    V=a*b
    if psit>0:
        return 1-V
    else:
        return V    

def getVolumeFraction2d(psi0, psi1, psi2):
    if psi0>0:
        if psi1>0:
            if psi2>0:
                return 0
            else:
                return getVolumeOneNode2d(psi2, psi0, psi1)
        else:
            if psi2>0:
                return getVolumeOneNode2d(psi1, psi0, psi2)
            else:
                return getVolumeOneNode2d(psi0, psi1, psi2)
    else:
        if psi1>0:
            if psi2>0:
                return getVolumeOneNode2d(psi0, psi1, psi2)
            else:
                return getVolumeOneNode2d(psi1, psi0, psi2)
        else:
            if psi2>0:
                return getVolumeOneNode2d(psi2, psi0, psi1)
            else:
                return 1

def interpolateLevelSetToElems(levelset_p1, func_p0, mesh, material):
    for el in mesh.Elements(VOL):
        if el.mat == material:
            V=levelset_p1.space
            L=func_p0.space
            dofs=V.GetDofNrs(el)
            psi0 = levelset_p1.vec[dofs[0]]
            psi1 = levelset_p1.vec[dofs[1]]
            psi2 = levelset_p1.vec[dofs[2]]
            func_p0.vec[L.GetDofNrs(el)[0]] = getVolumeFraction2d(psi0, psi1, psi2)  #s... volume fraction of negative part of triangle 