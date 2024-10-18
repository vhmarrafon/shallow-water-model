#!/bin/bash

cat <<EOF> plot.gs
'reinit'
'set display color white'
'c'
"open shallow_water.ctl"
'set gxout shaded'
t=1

while(t<=4000)
'set grid off'
'set grads off'
'set t 't
'color -.4 .4 .1 -kind navy->white->darkred'
'd n'

'!sleep 0.01'
'set gxout vector'
'set ccolor 1'


'set arrlab on '
'set arrscl 0.5 1.5'
'set arrowhead 0.05'
'd skip(u,4);skip(v,4)'

'!sleep 0.025'

say t
t=t+5
endwhile
pull n
'quit'
EOF

grads -l -c plot.gs
