time python integrate_hyper.py -hfd cv_outputs.txt -ind cv_integrated.txt -hl lasso/hyperlist_lasso.p
time python set_hyper.py -ind cv_integrated.txt -r hyper/hyper_df.txt -o hyper/best_hyper.p -hl lasso/hyperlist_lasso.p -tn lasso 
