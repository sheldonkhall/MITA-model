#/bin/bash

NOW=$(date)
FNAME="test results $NOW.txt"
touch "$FNAME"
for f in tests/*.py; do
    NOW=$(date +"%T")
    echo "$NOW Start of: $f" >> "$FNAME"
    ipython "$f" >> "$FNAME"
    NOW=$(date +"%T")
    echo "$NOW End of: $f" >> "$FNAME"
done
