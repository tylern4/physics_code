ls | grep "\.hbook" $1 | cut -c1-10 | while read x; do h2root $x".hbook" $x".root"; done
