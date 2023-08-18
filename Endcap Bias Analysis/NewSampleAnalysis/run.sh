#!/bin/bash

source /cvmfs/ilc.desy.de/key4hep/luxe_setup.sh

#tar -xzf bibData.tar.gz

# Run the geolocation code

wget -r -np -Q100000000k -P srv/data/ -R "index.html*" http://stash.osgconnect.net/collab/project/snowmass21/data/muonc/fmeloni/LegacyProductions/before29Jul23/DataMuC_MuColl_v1/muonGun/recoBIB/ --progress=bar:force 2>&1 | tail -f -n +6


PYTHONPATH=src:$PYTHONPATH python DataGathering.py 
