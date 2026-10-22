#!/bin/csh -f

set RR = "-R-124.1/-114/31.5/41"
#gmt makecpt -Cjet -Qo -T1/1000/1 -Z -D > vel.cpt
cp strain.cpt vel.cpt
gmt makecpt -Cgray -T-1.8/.5/.1 -Z > grad.cpt
set uu = GNSS_um.grd
set vv = GNSS_vm.grd
set output = GNSS_strain.ps
set uu = InSAR_GNSS_um.grd
set vv = InSAR_GNSS_vm.grd
set output = InSAR_GNSS_strain.ps
gmt set PS_MEDIA A2
gmt set FONT_ANNOT_PRIMARY 30
set scale = -Si0.00001i

# compute strain
gmt grdsample dem.grd -R$uu -Gdem_tmp.grd
gmt grdmath dem_tmp.grd 0 GT 0. NAN = mask.grd
gmt grdgradient $uu -A270 -fg -Gexx.grd
gmt grdgradient $uu -A180 -fg -Gexy.grd
gmt grdgradient $vv -A270 -fg -Geyx.grd
gmt grdgradient $vv -A180 -fg -Geyy.grd

gmt grdmath exy.grd eyx.grd ADD 2 DIV = e12.grd
gmt grdmath exx.grd exx.grd MUL eyy.grd eyy.grd MUL ADD e12.grd e12.grd MUL 2 MUL ADD SQRT 1E9 MUL = e_2nd.grd

gmt grdmath exx.grd 1E9 MUL mask.grd MUL = exx.grd
gmt grdmath exy.grd 1E9 MUL mask.grd MUL = exy.grd
gmt grdmath eyx.grd 1E9 MUL mask.grd MUL = eyx.grd
gmt grdmath eyy.grd 1E9 MUL mask.grd MUL = eyy.grd
gmt grdmath exy.grd eyx.grd ADD 2 DIV = e12.grd

# compute principal strain axis
gmt grdmath 2 e12.grd MUL exx.grd eyy.grd SUB ATAN2 2 DIV = a.grd
gmt grdmath exx.grd eyy.grd SUB 2 POW 0.25 MUL e12.grd 2 POW ADD SQRT exx.grd eyy.grd ADD 2 DIV ADD = e1.grd
gmt grdmath exx.grd eyy.grd SUB 2 POW 0.25 MUL e12.grd 2 POW ADD SQRT NEG exx.grd eyy.grd ADD 2 DIV ADD = e2.grd

gmt grdmath e1.grd a.grd COS MUL = e1_x.grd
gmt grdmath e1.grd a.grd SIN MUL = e1_y.grd
gmt grdmath e2.grd a.grd SIN NEG MUL = e2_x.grd
gmt grdmath e2.grd a.grd COS MUL = e2_y.grd

gmt grdgradient dem_tmp.grd -A315 -Nt.7 -Gdem_grad.grd -Vq
gmt grdimage e_2nd.grd -Cvel.cpt -Idem_grad.grd -JX10i $RR -K -Bxa4f1 -Bya2f1 -Bwesn -Q > $output
gmt pscoast -J -R -N1/0.5p,black -N2/0.3p,black -S153/204/255 -W1/0.5p,black -N3/0.5p,black -Dh -K -O >> $output
gmt psxy Historic.gmt -W0.1p,100/100/100 -J -R -K -O -P >> $output
gmt psxy Holocene.gmt -W0.1p,100/100/100 -J -R -K -O -P >> $output
gmt psxy data_coverage.gmt -W1p,black,-- -J -R -K -O -P >> $output

gmt grdclip e1.grd  -Gmask.grd -Sb0/0 -Sa100/100
gmt grdmath mask.grd e1_x.grd MUL = tx.grd
gmt grdmath mask.grd e1_y.grd MUL = ty.grd
gmt grdvector tx.grd ty.grd -J -R -Q+jc -W2.5,black -I0.2 -K -O $scale >> $output
gmt grdclip e1.grd  -Gmask.grd -Sa0/0 -Sb-100/-100
gmt grdmath mask.grd e1_x.grd MUL = tx.grd
gmt grdmath mask.grd e1_y.grd MUL = ty.grd
gmt grdvector tx.grd ty.grd -J -R -Q+jc -W2,blue -I0.2 -K -O $scale >> $output
gmt grdclip e2.grd  -Gmask.grd -Sb0/0 -Sa100/100
gmt grdmath mask.grd e2_x.grd MUL = tx.grd
gmt grdmath mask.grd e2_y.grd MUL = ty.grd
gmt grdvector tx.grd ty.grd -J -R -Q+jc -W2.5,black -I0.2 -K -O $scale >> $output
gmt grdclip e2.grd  -Gmask.grd -Sa0/0 -Sb-100/-100
gmt grdmath mask.grd e2_x.grd MUL = tx.grd
gmt grdmath mask.grd e2_y.grd MUL = ty.grd
gmt grdvector tx.grd ty.grd -J -R -Q+jc -W2,blue -I0.2 -K -O $scale >> $output

gmt psscale -Dx2c/3c+w8c/0.5c+h -Cvel.cpt -B+l"nanostrain/yr" --FONT_LABEL=30 -O -Q >> $output

gmt psconvert $output -Tf -A

#open $output






