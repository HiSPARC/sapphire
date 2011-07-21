"""
Classes corresponding to corsika blocks and sub-blocks

The classes in this module correspond one-to-one with the sub-blocks (and blocks) as specified in corsika's user manual.

Authors
-------
Javier Gonzalez <jgonzalez@ik.fzk.de>
"""

import struct
import CorsikaUnits as units
import numpy

# All sizes are in bytes

# in 32 bit, one field is a float
gFieldSize = struct.calcsize('f')

# one block contains 6552 fields plus one header and one trailer field
gBlockFormat = '6554f'
gBlockSize = struct.calcsize(gBlockFormat)
gBlockPaddingSize = struct.calcsize('f')

# Each block contains 21 sub-blocks
# Each sub-block consists of 312 fields
# the first of which _might_ be a string id.
gSubblockFormat = '4s311f'
gSubblockSize = struct.calcsize(gSubblockFormat)
gSubBlocksPerBlock = 21

# Each particle record sub-block contains a fixed
# number of particle records
# With the thinned option, each of these is 8 fields long
# for a total of 39 records per sub block
gParticleFormat = '8f'
gParticlesPerSubblock = 39
gParticleRecordSize = struct.calcsize(gParticleFormat)


# a couple of sanity checks for the formats
if (gBlockSize-2*gBlockPaddingSize)/gSubblockSize != gSubBlocksPerBlock:
    raise Exception('The block format ({block}) and sub-block format ({sub_block}) do not agree! block size is {block_size} and sub-block size is {sub_block_size}. Block size should be {subblocks_per_block} times the sub-block size plus padding (usually 8 bytes).'.format(block=gBlockFormat, sub_block=gSubblockFormat, block_size=gBlockSize, sub_block_size=gSubblockSize, subblocks_per_block=gSubBlocksPerBlock))


if gSubblockSize/gParticleRecordSize != gParticlesPerSubblock:
    raise Exception('The sub_block format ({sub_block}) and particle format ({particle}) do not agree! sub-block size is {sub_block_size} and particle record size is {particle_size}. Sub-block size should be {particles_per_subblock} times the particle record size.'.format(sub_block=gSubblockFormat, particle=gParticleFormat, sub_block_size=gSubblockSize, particle_size=gParticleRecordSize, particles_per_subblock=gParticlesPerSubblock))


# From here on, things should not depend on the field size as everything is

class RunHeader:
    """
    Class representing the run header sub-block
    as specified in Corsika user manual, section 10.2
    """
    def __init__(self, subblock):
        self.fId = subblock[0]
        self.fRunNumber = subblock[1]
        self.fDateStart = subblock[2]
        self.fVersion = subblock[3]

        self.fObservationLevels = subblock[4]
        self.fObservationHeight = numpy.array(subblock[5:15])*units.cm

        self.fSpectralSlope = subblock[15]
        self.fEMin = subblock[16]*units.GeV
        self.fEMax = subblock[17]*units.GeV

        self.fFlagEGS4 = subblock[18]
        self.fFlagNKG = subblock[19]

        self.fCutoffHadrons = subblock[20]*units.GeV
        self.fCutoffMuons = subblock[21]*units.GeV
        self.fCutoffElectrons = subblock[22]*units.GeV
        self.fCutoffPhotons = subblock[23]*units.GeV

        self.fConstAATM = numpy.array(subblock[254:259])
        self.fConstBATM = numpy.array(subblock[259:264])
        self.fConstCATM = numpy.array(subblock[264:269])
        self.fConstNFLAIN = subblock[269]
        self.fConstNFLDIF = subblock[270]
        self.fConstNFLPI = subblock[271]
        self.fConstNFLCHE = subblock[272]
    def __str__(self):
        return '''Run header:
  id: {run_n}
  date: {date}
  version: {version}
  '''.format(run_n=self.fRunNumber,
             date=self.fDateStart,
             version=self.fVersion)

class EventHeader:
    """
    Class representing the event header sub-block
    as specified in Corsika user manual, section 10.2
    """
    def __init__(self, subblock):
        self.fId = subblock[0]
        self.fEventNumber = subblock[1]
        self.fParticleId = subblock[2]
        self.fEnergy = subblock[3]*units.GeV
        self.fStartingAltitude = subblock[4]*units.g/units.cm2
        self.fFirstTarget = subblock[5]
        self.fZFirst = subblock[6]*units.cm2
        self.fPx = subblock[7]*units.GeV
        self.fPy = subblock[8]*units.GeV
        self.fPz = subblock[9]*units.GeV
        self.fTheta = subblock[10]*units.rad
        self.fPhi = subblock[11]*units.rad

        self.fNofRandomSequences = subblock[12]
        self.fRandomSequences = numpy.array(subblock[13:43])
        self.fRunNumber = subblock[43]
        self.fDateStart = subblock[44]
        self.fVersion = subblock[45]

        self.fObservationLevels = subblock[46]
        self.fObservationHeight = numpy.array(subblock[47:57])*units.cm

        self.fSpectralSlope = subblock[57]
        self.fEMin = subblock[58]*units.GeV
        self.fEMax = subblock[59]*units.GeV

        self.fCutoffHadrons = subblock[60]*units.GeV
        self.fCutoffMuons = subblock[61]*units.GeV
        self.fCutoffElectrons = subblock[62]*units.GeV
        self.fCutoffPhotons = subblock[63]*units.GeV

        self.fNFLAIN = subblock[64]
        self.fNFLDIF = subblock[65]
        self.fNFLPI0 = subblock[66]
        self.fNFLPIF = subblock[67]
        self.fNFLCHE = subblock[68]
        self.fNFRAGM = subblock[69]

        self.fBx = subblock[70]*units.micro*units.tesla
        self.fBz = subblock[71]*units.micro*units.tesla
        self.fFlagEGS4 = subblock[72]
        self.fFlagNKG = subblock[73]

        self.fFlagGeisha = subblock[74]
        self.fFlagVenus = subblock[75]
        self.fFlagCerenkov = subblock[76]
        self.fFlagNeutrino = subblock[77]
        self.fFlagCurved = subblock[78]
        self.fFlagComputer = subblock[79]
        self.fThetaMin = subblock[80]*units.rad
        self.fThetaMax = subblock[81]*units.rad
        self.fPhiMin = subblock[82]*units.rad
        self.fPhiMax = subblock[83]*units.rad

        self.fCerenkovBunch = subblock[84]
        self.fCerenkovNumberX = subblock[85]
        self.fCerenkovNumberY = subblock[86]
        self.fCerenkovGridX = subblock[87]*units.cm
        self.fCerenkovGridY = subblock[88]*units.cm
        self.fCerenkovDetectorX = subblock[89]*units.cm
        self.fCerenkovDetectorY = subblock[90]*units.cm
        self.fCerenkovOutputFlag = subblock[91]

        self.fArrayRotation = subblock[92]*units.rad
        self.fFlagExtraMuonInformation = subblock[93]

        self.fMultipleScatteringStep = subblock[94]
        self.fCerenkovBandwidthMin = subblock[95]*units.nanometer
        self.fCerenkovBandwidthMax = subblock[96]*units.nanometer
        self.fUsersOfEvent = subblock[97]
        self.fCoreX = numpy.array(subblock[98:118])*units.cm
        self.fCoreY = numpy.array(subblock[118:138])*units.cm

        self.fFlagSIBYLL = subblock[138]
        self.fFlagSIBYLLCross = subblock[139]
        self.fFlagQGSJET = subblock[140]
        self.fFlagQGSJETCross = subblock[141]
        self.fFlagDPMJET = subblock[142]
        self.fFlagDPMJETCross = subblock[143]
        self.fFlagVENUSCross = subblock[144]
        self.fFlagMuonMultiple = subblock[145]
        self.fNKGRadialRange = subblock[146]*units.cm
        self.fEFractionThinningH = subblock[147]
        self.fEFractionThinningEM = subblock[148]
        self.fWMaxHadronic = subblock[149]
        self.fWMaxEM = subblock[150]
        self.fRMaxThinning = subblock[151]*units.cm
        self.fInnerAngle = subblock[152]*units.rad
        self.fOuterAngle = subblock[153]*units.rad

        self.fHighEnergyLowEnergyTransition = subblock[154]*units.GeV
        self.fFlagSkimmingIncidence = subblock[155]
        self.fAltitudeHorizontalShowerAxis = subblock[156]*units.cm
        self.fStartingHeight = subblock[157]*units.cm
        self.fFlagCharm = subblock[158]
        self.fFlagHadronOrigin = subblock[159]
        self.fFlagObseLevelCurvature = subblock[167]

    def __str__(self):
        return '''Event header:
  id: {event_n}
  primary: {primary}
  energy: {energy} EeV
  direction: ({zenith}, {azimuth})
  '''.format(event_n=self.fEventNumber,
             primary=self.fParticleId,
             energy=self.fEnergy/units.EeV,
             zenith=self.fTheta/units.degree,
             azimuth=self.fPhi/units.degree)

class RunTrailer:
    """
    Class representing the run end sub-block
    as specified in Corsika user manual, section 10.2.
    """
    def __init__(self, subblock):
        self.fID = subblock[0]
        self.fRunNumber = subblock[1]
        self.fEventsProcessed = subblock[2]
    def __str__(self):
        return '''Run trailer:
  id: {run_n}
  events: {events}
'''.format(run_n=self.fRunNumber,
           events=self.fEventsProcessed)

class EventTrailer:
    """
    Class representing the event end sub-block
    as specified in Corsika user manual, section 10.2.
    """
    def __init__(self, subblock):
        self.fID = subblock[0]
        self.fEventNumber = subblock[1]

        self.fPhotons = subblock[2]
        self.fElectrons = subblock[3]
        self.fHadrons = subblock[4]
        self.fMuons = subblock[5]
        self.fParticles = subblock[6]

        # NKG output
        self.fLateral1X = numpy.array(subblock[7:28])/units.cm2
        self.fLateral1Y = numpy.array(subblock[28:49])/units.cm2
        self.fLateral1XY = numpy.array(subblock[49:70])/units.cm2
        self.fLateral1YX = numpy.array(subblock[70:91])/units.cm2

        self.fLateral2X = numpy.array(subblock[91:112])/units.cm2
        self.fLateral2Y = numpy.array(subblock[112:133])/units.cm2
        self.fLateral2XY = numpy.array(subblock[133:154])/units.cm2
        self.fLateral2YX = numpy.array(subblock[154:175])/units.cm2

        self.fElectronNumber = numpy.array(subblock[175:185])
        self.fAge = numpy.array(subblock[185:195])
        self.fDistances = numpy.array(subblock[195:205])*units.cm2
        self.fLocalAge1 =numpy.array( subblock[205:215])

        self.fLevelHeightMass = numpy.array(subblock[215:225])
        self.fLevelHeightDistance = numpy.array(subblock[225:235])
        self.fDistanceBinsAge = numpy.array(subblock[235:245])*units.cm2
        self.fLocalAge2 = numpy.array(subblock[245:255])

        # Longitudinal distribution
        self.fLongitudinalPar = numpy.array(subblock[255:261])
        self.fChi2 = subblock[261]

        # Added according to the CORSIKA manual
        self.fWeightedPhotons = subblock[262]
        self.fWeightedElectrons = subblock[263]
        self.fWeightedHadrons = subblock[264]
        self.fWeightedMuons = subblock[265]
        self.fPreshowerEMParticles = subblock[266]

    def __str__(self):
        return '''Event trailer:
  id: {event_n}
  particles: {particles}
  '''.format(event_n=self.fEventNumber,
             particles=self.fParticles)


class ParticleData:
    """
    Class representing the particle data sub-block
    as specified in Corsika user manual, section 10.2.

    The number of CherenkovData records in a sub-block depends on
    compilation options.
    """
    def __init__(self, subblock):
        self.fDescription = subblock[0]
        self.fPx = subblock[1]*units.GeV
        self.fPy = subblock[2]*units.GeV
        self.fPz = subblock[3]*units.GeV
        self.fX = subblock[4]*units.cm2
        self.fY = subblock[5]*units.cm2
        self.fTorZ = subblock[6]*units.ns
        self.fWeight = subblock[7]

    def IsParticle():
        return 0 < self.fDescription and self.fDescription < 100000
    def IsNucleus():
        return 100000 <= self.fDescription and self.fDescription < 9900000
    def IsCherenkov():
        return 9900000 <= self.fDescription

    def __str__(self):
       return '''Particle:
  id: {description}
  momentum: {momentum}
  position: {position}
  time: {time}
  weight {weight}
  '''.format(description=self.fDescription,
             momentum=(self.fPx*units.GeV, self.fPy*units.GeV, self.fPz*units.GeV),
             position=(self.fX*units.m, self.fY*units.m),
             time=self.fTorZ,
             weight=self.fWeight)

class CherenkovData:
    """
    Class representing the cherenkov photon sub-block
    as specified in Corsika user manual, section 10.2.

    The number of CherenkovData records in a sub-block depends on
    compilation options.
    """
    def __init__(self, subblock):
        self.fPhotonsInBunch = subblock[0]
        self.fX = subblock[1]
        self.fY = subblock[2]
        self.fU = subblock[3]
        self.fV = subblock[4]
        self.fT = subblock[5]
        self.fProductionHeight = subblock[6]
        self.fWeight = subblock[7]

    def __str__(self):
       return '''Cherenkov:
  n: {n}
  position: {position} m
  (u,v): {direction}
  time: {time} ns
  prod. height: {height} m
  '''.format(n=self.fPhotonsInBunch,
             direction=(self.fU, self.fV),
             position=(self.fX*units.m, self.fY*units.m),
             time=self.fT/units.nanosecond,
             height=self.fWeight/units.m)
