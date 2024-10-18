#!/bin/bash
rm -r fig
mkdir fig
cat <<EOF> plot.gs
'reinit'
'set display color white'
'c'
"open shallow_water.ctl"
'set gxout shaded'
t=0
pull n
while(t<=15000)
'set grid off'
'set grads off'
'set t 't
'color -.4 .4 .1 -kind navy->white->darkred'
'd n'

*'!sleep 0.01'
'set gxout vector'
'set ccolor 1'


'set arrlab on '
'set arrscl 1 0.75'
'set arrowhead 0.05'
'd skip(u,4);skip(v,4)'
*'printim f't'.png'

*'!sleep 0.025'

say t
t=t+10
endwhile
pull n
'quit'
EOF

grads -l -c plot.gs
mv *png fig
