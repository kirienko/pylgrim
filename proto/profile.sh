#!/bin/bash
python -m cProfile -o rover_pos.prof rover_position.py
PNG_name='profile_results.png'
if test -f ${PNG_name}  ; then
    i=1
    while test -f ${PNG_name}.${i} ; do
        let i++
    done
    PNG_name=${PNG_name}.${i}
fi
gprof2dot -f pstats rover_pos.prof | dot -Tpng -o ${PNG_name}
echo ${PNG_name} created
xdg-open ${PNG_name}
