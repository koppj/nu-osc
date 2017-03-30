#!/bin/bash 

CL='0.95'

cntxmgr  -c$CL -f"XYZ:%lf %lf %lf %*f  %*f " $1 > cnt1.dat
cntxmgr  -c$CL -f"XYZ:%lf %lf %*f %lf  %*f " $1 > cnt2.dat
cntxmgr  -c$CL -f"XYZ:%lf %lf %*f %*f  %lf " $1 > cnt3.dat

GRAPH='0'
LW='2'
LS='1'

#add2xmgr  -g$GRAPH=[md=3,lw=$LW,lc=11,ls=$LS]cnt1.dat $2 > tt.agr 

add2xmgr  -g$GRAPH=[lw=$LW,lc=2,ls=$LS]cnt3.dat \
          -g$GRAPH=[md=3,lw=$LW,lc=4,ls=$LS]cnt1.dat \
          -g$GRAPH=[md=3,lw=$LW,lc=1,ls=$LS]cnt2.dat \
$2 > tt.agr 

rm -f cnt*

xmgrace tt.agr
