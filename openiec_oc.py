"""According to given conditions and input parameters, calculate corresponding interfacial energies.
"""
import os
import pyOC
import numpy as np
import matplotlib.pyplot as plt
from pyOC import opencalphad as oc
from pyOC import PhaseStatus as phStat
from pyOC import GridMinimizerStatus as gmStat


##Calculate the phase equilibrium, Gibbs energy and Chemical potentials of the associated components in the bulk phases
class CoherentGibbsEnergy(object):
    """
    Equilibrium calculation for sing phase or two  or more phases.

    Parameters 
    ----------
    T: float- Given temperature.
    db : Database- Database containing the relevant parameters.
    comp: str- Name of pure component.
    phasename: list- Name of phase model to build.
    composition: float- Composition of the element in the phases
    -----------
    Returns
    G: Molar Gibbs energy of the bulk phases
    Mu: Chemical potential of the component in bulk phases
    Equilibrium Phases 
    Equilibrium phases composition
    Mole fraction of componets in phase equilibrium
    Sitefractions of compoent in equilibrium
    -----------
    """
    ##Calculate equilibrium of the bulk phases

def Gm_bulkphase(temperature,elementMolarAmounts):
	# set temperature
	oc.setTemperature(temperature)
	oc.setElementMolarAmounts(elementMolarAmounts)
	Gm_bulkphase={}	
	Gm_bulkphase = oc.calculateEquilibrium(gmStat.On)
	G=oc.getScalarResult('G')
	phase_description = oc.getPhasesAtEquilibrium().getPhaseConstituentComposition()
	chemicalpotential=oc.getComponentAssociatedResult('MU')

######
# Functions for evaluating the equilibrium and chemical potential of the components (global and local) 
######

class MolarVolume(object):
    """
    Evaluate partial molar volumes by an approximation of the first order volume derivative by a second-order finite difference formula
    Parameters (additonal, not listed above)
    -----------
    mass density laws (from Barrachin2004)
    """
######
# Functions for partial molar volumes evaluation
######
## convert constituent molar fractions to mass fractions
## function to calculate molar volume from constituent
## evaluate partial molar volumes by an approximation of the first order volume derivative by a second-order finite difference formula
##Evaluate molar interfacial area of the component

class SigmaCoherentInterface(object):
     """
    Contruct the interfacial energy of the coherent interface and their partial quantaties.

    Parameters
    -----------
    alphafuncs: list of functions - The phase with a single sublattice.
    betafuncs: tuple of functions - The phase with two sublattices.
    mueq: list - The chemical potentials of the equilibrium state at the given compositions.
    vmis: list - The partial molar volumes of components.
    """
def infenergy(self, x):
    """
        Compute coherent partial interfacial energies of various components.
        Parameters
        ----------
        x: list - Interfacial composition.
    """
def objective(self, x):
        """
        Compute the absolute value of the differences between the partial interfacial energy.
        
        Parameters
        ----------
        x: list - Interfacial composition.
    """
def SigmaCoherent(T, x0, db, comps, phasenames, purevms, intervms=[], limit=[0, 1.0], dx=0.01):
    """
    Contruct the interfacial energy of the coherent interface and their partial quantaties.
    Inputs
    ----------
    Calculate the coherent interfacial energy in alloys.

    Parameters
    -----------
    T: float - Given temperature.
    x0: list - Initial alloy composition.
    db: Database - Database containing the relevant parameters.
    comps : list - Names of components to consider in the calculation.
    phasenames : list - Names of phase model to build.    
    limit: list - The limit of composition for searching interfacial composition in equilibrium.
    purevms: list - The molar volumes of pure components.
    dx: float - The step of composition for searching interfacial composition in equilibrium.

    Returns:   
    -----------
    Components：list of str - Given components.
    Temperature: float - Given temperature.
    Initial_Alloy_Composition: list - Given initial alloy composition.
    Interfacial_Composition: list - Interfacial composition of the grid minimization.
    Partial_Interfacial_Energies: list - Partial interfacial energies of components.
    Interfacial_Energy: float - Requested interfacial energies.
    Return type: xarray Dataset
    """   
    
# oc setup
## setting verbosity (True or False - default), if set to yes, in particular, when getters are called the returned values are displayed in a comprehensive way
oc.setVerbosity(True)
## tdb filepath
tdbFile=os.environ.get('OCDATA')+'/feouzr.tdb'
#reading tdb
elems=('O', 'U', 'ZR')
oc.readtdb(tdbFile,elems)
## suspend all phases except the liquid one
##oc.setPhasesStatus(('C1_FCC',),phStat.Suspended)
oc.setPhasesStatus(('LIQUID', 'C1_FCC'),phStat.Entered)
## set pressure
oc.setPressure(1E5)

# mass density laws (from Barrachin2004)
coriumMassDensityLaws = {
	'U1'   : lambda T: 17270.0-1.358*(T-1408),
	'ZR1'  : lambda T: 6844.51-0.609898*T+2.05008E-4*T**2-4.47829E-8*T**3+3.26469E-12*T**4,
	'O2U1' : lambda T: 8860.0-9.285E-1*(T-3120),
	'O2ZR1': lambda T: 5150-0.445*(T-2983),
	'O1'   : lambda T: 1.141 # set to meaningless value but ok as, no 'free' oxygen in the considered mixtures
}

# temperature and composition for which partial molar volumes are to be evaluated
temperature=3000
elementMolarAmounts = {
	'U' : 0.343298,
	'O' : 0.414924,
	'ZR': 0.241778
}


