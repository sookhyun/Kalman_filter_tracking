# Kalman_filter_tracking
This repository contains updates to the initial implementation of the **Kalman Filter** and **Cellular Automaton** methods for the determination of *multiple collision event vertices* with a resolving power of a few ${\micro}m$ at the sPHENIX experiment. 

The sPHENIX charged particle tracking system consists of 3 layers of MVTX (Monolithic Active Pixel Sensor based VerTeX identification detector), 2 layers of INTT (INTermediate strip pixel silicon sensor Tracker) and 68 layers of TPC (gas-electron-multiplier based Time Projection Chamber). 

Initial implementation utilized triplets in the innermost MVTX to reconstruct all collision vertices present at a given crossing. 

$\bullet$ inital commit - https://github.com/sPHENIX-Collaboration/coresoftware/pull/190 

$\bullet$ compatibility update - https://github.com/sPHENIX-Collaboration/coresoftware/pull/687

Note that the code here is outdated and incompatible with current sPHENIX software, which has evolved continuously since the development of this code. The main algorithms for vertex finding using 3 hits on MVTX layers are implemented in PHG4InitZVertexing and the ones for extending triplets to quadruplets by adding an INTT hit are found in PHG4PatternReco in source directory. Both heavily utilize CellularAutomaton and HelixKalmanFilter. 

