#!/bin/bash

cat <<EOF> plot.gs
'reinit'
'set display color white'
'c'
"open ADV${1}.ctl"
t=1
while(t<=88)
'c'
'set t 't
'd c'
say t
t=t+1
endwhile
pull n
'quit'
EOF

grads -l -c plot.gs
