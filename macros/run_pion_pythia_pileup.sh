#!/bin/csh 
echo "Executing run_pion_pythia_pileup.sh"
set workdir = /sphenix/user/shlee/coresoftware/patreco
if ($?_CONDOR_SCRATCH_DIR) then
set workdir = ${_CONDOR_SCRATCH_DIR}
endif
cd $workdir
source /opt/sphenix/core/bin/sphenix_setup.csh -n

# nevent, number of tracks in an event, zsigma, zsigmapileup
root -l -b -q Fun4All_G4_sPHENIX.C\(3,$4,$5,$6,$7,\"/sphenix/user/shlee/hepmc_pythia/minbias/phpythia8_200gev_mb_00$1$2$3.dat.gz\",\"/sphenix/user/shlee/hepmc_pythia/minbias_pileup/phpythia8_200gev_mb_00$1$2$3.dat.gz\",\"tracktree_pythia_pileup_$4_$5_$6_$7_$1$2$3.root\",\"multi_zvtx_pileup_$4_$5_$6_$7_$1$2$3.root\" \)

# use ittf pileup file
root -l -b -q Fun4All_G4_sPHENIX.C\(3,$4,$5,$6,$7,\"/sphenix/user/shlee/hepmc_pythia/minbias/phpythia8_200gev_mb_00$1$2$3.dat.gz\",\"/sphenix/user/shlee/hepmc_pythia/minbias_pileup/phpythia8_200gev_mb_00$1$2$3.dat.gz\",\"tracktree_pythia_pileup_$4_$5_$6_$7_$1$2$3.root\",\"multi_zvtx_pileup_$4_$5_$6_$7_$1$2$3.root\" \)

echo "The run_pion_pythia_pileup.sh script has finished executing, exiting".
exit
