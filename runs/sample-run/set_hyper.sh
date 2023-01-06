#!/bin/bash
START=$(date)
set -e
time python3 integrate_hyper.py -hfd cv_outputs.txt -ind cv_integrated.txt -hl enet/hyperlist_enet.p
time python3 set_hyper.py -ind cv_integrated.txt -r hyper/hyper_df.txt -o hyper/best_hyper.p -hl enet/hyperlist_enet.p -tn enet 
END=$(date)
echo set_hyper.sh,$START,$END,$SECONDS >> timing/result_time.csv
