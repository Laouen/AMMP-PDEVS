# Tuned version of matemarote/utils/svg2png.sh to crop svgs from excesive
# transparency
if [ -z "$1" ]; then
    files="*.svg"
else
    files=$1
fi
{ for x in *.svg; do inkscape -D --export-png="$x.png" "$x"; done }
rename -f 's/(.*).svg.png/$1.png/' *.svg.png
