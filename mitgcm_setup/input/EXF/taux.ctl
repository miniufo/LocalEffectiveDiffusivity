dset ^taux.dat
undef 0
options big_endian
title sin wind stress
xdef 560 linear 140 0.1
ydef 400 linear  15 0.1
zdef   1 linear   0 1
tdef   1 linear 00z1Jan2000 1dy
vars 1
taux 0 99 taux
endvars

