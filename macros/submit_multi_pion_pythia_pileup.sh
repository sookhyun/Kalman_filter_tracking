#!/bin/csh 

set num_of_tracks = 1

set runs_in_group = 10

foreach ggroup (0 1 2 3 4 5 6 7 8 9)
foreach group (0 1 2 3 4 5 6 7 8 9)

set zsigma = 10.0
set zsigmapileup = 40.0
set collisionrate = 12000000;

    mkdir -p /sphenix/user/shlee/coresoftware/patreco18/LOGS
    mkdir -p /sphenix/user/shlee/coresoftware/patreco18/output

    set run = 0

    while ( $run < $runs_in_group )
    
      ./submit_pion_pythia_pileup.sh $ggroup $group $run $num_of_tracks $zsigma $zsigmapileup $collisionrate

      @ run++

    end 

end
end
