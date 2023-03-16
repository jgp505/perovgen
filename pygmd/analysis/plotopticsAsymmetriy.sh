awk 'BEGIN{i=0} /<dielectricfunction>/,\
                /<\/dielectricfunction>/ \
                 {if ($1=="<r>") {a[i]=$2 ; b[i]=$3 ; c[i]=$4 ; d[i]=$5 ; i=i+1}} \
     END{for (j=0;j<i/2;j++) print a[j],b[j],b[j+i/2],c[j],c[j+i/2],d[j],d[j+i/2]}' vasprun.xml > optics.dat

cat >plotfile<<!

# set term postscript enhanced eps colour lw 2 "Helvetica" 20
# set output "optics.eps"

set xrange [0:15]

plot "optics.dat" using (\$1):(\$3)  w l  lt -1 lw 1      lc  3 title "IPA real", \
     "optics.dat" using (\$1):(-\$2) w l  lt  0 lw 2      lc  3 title "IAP imag"
!

#gnuplot -persist plotfile
