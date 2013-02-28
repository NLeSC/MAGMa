typew={"AROMATIC":3.0,\
       "DOUBLE":2.0,\
       "TRIPLE":3.0,\
       "SINGLE":1.0}
global missingfragmentpenalty
heterow={False:2,True:1}
missingfragmentpenalty=10

mims={'H':1.0078250321,\
      'C':12.0000000,\
      'N':14.0030740052,\
      'O':15.9949146221,\
      'F':18.99840320,\
      'P':30.97376151,\
      'S':31.97207069,\
      'Cl':34.96885271,\
      'Br':78.9183376,\
      'I':126.904468}

Hmass=mims['H']     # Mass of hydrogen atom
elmass=0.0005486

fragmentation_engine='rdkit' # choose rdkit or cdk