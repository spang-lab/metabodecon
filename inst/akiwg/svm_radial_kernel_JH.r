# Note for SVM package e1071 is required
# Note for feature selection the package multtest is required can only be installed from bioconductordouble
# Note for balanced sampling the functions included in the extra file fun_sampling.r are required
# Note for calculation of ROC curves the functions prediction and performance from package ROCR are required
install.packages("ROCR") # for performance measures
library(ROCR)
# requires package multtest; for that we first we need to install Bioconductor
    if (!requireNamespace("BiocManager", quietly = TRUE))

    install.packages("BiocManager")

    BiocManager::install()
# now install multtest for multiple testing adjustments
BiocManager::install("multtest")
library("multtest")
install.packages("e1071")
library("e1071")

fsel <- function (x,y,n) # x data, y labels, n how many returned
{
	dat<-t(x) # transpose back to normal order only for SVM
	labels <-as.integer(y)-1 #labels
	tscores <- mt.teststat(dat,labels,test="t") #tscores
	sel<-order(-abs(tscores))[1:n] #select best
	return(sel)  #return selection
}


double_nested_cv <- function (x_out, y_out, leave_out, balanced_sampling, max_feat)
{
	# Erklärung der Parameter:

	# x_out: 		Datenmatrix mit Samples in Zeilen und Features in Spalten
	# y_out: 		Vektor mit Klassenlabels as factors
	# leave_out: 		Umfang der Subsets, die herausgelassen werden
	# balanced_sampling: 	Boolsche Variable zum an-/ausschalten der Stratified Crossvalidation
	# max_feat:		Höchstzahl Features bis zu der optimiert wird

	# noch offen: Höchstzahl Features
	
	rp <- 1
	avg_acc <- as.numeric(rp)
	
	sel_fin_all <- as.numeric(nrow(x_out)/leave_out*rp) 
	sel_cost_all <- as.numeric(nrow(x_out)/leave_out*rp) 
	sel_gamma_all <- as.numeric(nrow(x_out)/leave_out*rp) 
	
	r <- 1
	rc <- 1
	pprob <- numeric(nrow(x_out)*rp)
	lab <- numeric(nrow(x_out)*rp)
	predicted_labels <- numeric(nrow(x_out)*rp) 
	cat("Initialisierung...\n")
	for (j_out in 1:rp)
	{
		cat("Wiederholungsschleife: ",j_out,"\n")
		
		if (balanced_sampling)
		{
			heaps_out <- makeheaps(y_out, leave_out)
		}
		else
		{
			perm_out <- sample(1:nrow(x_out))
			heaps_out <- matrix(perm_out, nrow = leave_out)
		}
		
		k_out <- nrow(x_out)/leave_out
		
		acc_out <- numeric(k_out) 

		for (j in 1:k_out)
		{
			cat("Split vom ganzen Datensatz: ",j,"\n") 
			test_out <- heaps_out[,j]

			# Loop for tuning of variable selection

			# alt: tune_no <- 60
			tune_no <- max_feat
			
			tot_acc_in <- as.numeric(tune_no)
			tot_costs <- as.numeric(tune_no)
			tot_gammas <- as.numeric(tune_no)
			tot_acc_in_start <- 0
			for (m in 1:tune_no)
			{
				cat("Featurezahl: ",m,"\n")
				no_start <- 0
				no_sel <- no_start+m 
				x <- x_out[-test_out,]
				y <- y_out[-test_out]
				
				if (balanced_sampling)
				{
					heaps <- makeheaps(y, leave_out)
				}
				else
				{
					perm <- sample(1:nrow(x))
					heaps <- matrix(perm,nrow=leave_out)
				}

				k <- nrow(x)/leave_out

				acc <- numeric(k)
				cost_output <- numeric(k)
				gamma_output <- numeric(k)
				for (n in 1:k)
				{
					cat("Split fuer Featureselection: ",n,"\n")
					test <- heaps[,n]
		
					# Loop for cost-parameter tuning

					all_gamma <- rep(2^(-3:-8), 6)
					all_cost <- rep(2^(0:5), each = 6)
					parameters <- rbind(all_gamma, all_cost)
					tune_no_in <- dim(parameters)[2]
					tot_acc_in_in <- as.numeric(tune_no_in)
					tot_acc_in_in_start <- 0
					for (i in 1:tune_no_in)
					{
						gamma <- parameters[1,i]
						cost <- parameters[2,i]
						x_in <- x[-test,]
						y_in <- y[-test]
						
						if (balanced_sampling)
						{
							heaps_in <- makeheaps(y_in, leave_out)
						}
						else
						{
							perm_in <- sample(1:nrow(x_in))
							heaps_in <- matrix(perm_in,nrow = leave_out)
						}
						
						k_in <- nrow(x_in)/leave_out
						
						acc_in <- numeric(k_in)
						for (p in 1:k_in)
						{
							test_in <- heaps_in[,p]
							sel_in <- fsel(x_in[-test_in,], y_in[-test_in], no_sel)
							SVM <- svm(x_in[-test_in,sel_in], y_in[-test_in], kernel="radial", cost=cost, gamma=gamma)
							p_in <- predict(SVM, x_in[test_in, sel_in])
							acc_in[p] <- 100*sum(p_in==y_in[test_in])/length(test_in)
						}
						tot_acc_in_in[i] <- mean(acc_in)
						# cat("cost optimization loop", i, "selected features", no_sel, "accuracy", mean(acc_in),"\n")
						if (tot_acc_in_in[i] > tot_acc_in_in_start)
						{
							tot_acc_in_in_start <- tot_acc_in_in[i]
							cost_fin <- cost
							gamma_fin <- gamma
						} 
					}
					sel <- fsel(x[-test,], y[-test], no_sel)
					cat("selected features: ",sel,"\n","optimal costparameter:",cost_fin,"optimal gamma:",gamma_fin,"\n")
					SVM <- svm(x[-test,sel], y[-test], kernel="radial", cost=cost_fin, gamma=gamma_fin)
					p <- predict(SVM, x[test, sel])
					acc[n] <- 100*sum(p==y[test])/length(test)
					cost_output[n] <- cost_fin
					gamma_output[n] <- gamma_fin
				}
	
				tot_acc_in[m] <- mean(acc)
				tot_costs[m] <- median(cost_output)
				tot_gammas[m] <- median(gamma_output)
				cat("feature optimization loop", m, "optimal cost", tot_costs[m], "optimal gamma", tot_gammas[m], "accuracy", tot_acc_in[m],"\n")
				if (tot_acc_in[m] > tot_acc_in_start)
				{
					tot_acc_in_start <- tot_acc_in[m]
					sel_fin <- no_sel
					sel_cost <- tot_costs[m]
					sel_gamma <- tot_gammas[m]
				}
			}
			sel_fin_all[r] <- sel_fin
			sel_cost_all[r] <- sel_cost
			sel_gamma_all[r] <- sel_gamma
			r = r + 1 
			cat("outer loop", j, "repeat loop", j_out, "sel_fin", sel_fin, "sel_cost", sel_cost, "sel_gamma", sel_gamma, "tot_acc_in_start", tot_acc_in_start,"\n")
			sel <- fsel(x_out[-test_out,], y_out[-test_out],sel_fin)
			SVM <- svm(x_out[-test_out,sel], y_out[-test_out], kernel = "radial", cost = sel_cost, gamma = sel_gamma, probability = TRUE)
			p <- predict(SVM, x_out[test_out, sel],probability=TRUE)
			
			for (zaehler in 0:(leave_out-1))
			
			{
				pprob[rc+zaehler] <- attributes(p)$probabilities[zaehler+1,2]
				lab[rc+zaehler] <- as.numeric(y_out[test_out[zaehler+1]])-1
				predicted_labels[rc+zaehler] <- as.numeric(p[zaehler+1])-1
			}
		        acc_out[j] <- 100*sum(p==y_out[test_out])/length(test_out)
		        print(round(acc_out[j]))
			
			rc = rc + leave_out
			
		}
		avg_acc[j_out] <- mean(acc_out)
		print(round(acc_out))
   
   
	}
	sel_fin_end <-median(sel_fin_all) #averaged tuning parameter
	sel_cost_end <- median(sel_cost_all)
	sel_gamma_end <- median(sel_gamma_all) 
	pred<-prediction(pprob,lab) # ROC curve
   	perf<-performance(pred,measure="tpr",x.measure="fpr")
   	plot(perf,col=rainbow(10))
   	perf2<-performance(pred,measure="auc")
   	auc<-as.numeric(attributes(perf2)$y.values)
   	cat("Area under the ROC curve",auc,"\n")
   	return(list(avg.tot.acc=avg_acc, medi.sel.fin=sel_fin_end, medi.sel.cost=sel_cost_end, medi.sel.gamma=sel_gamma_end, true.labels=lab, probabilites=pprob, predicted.labels=predicted_labels, samples = heaps_out))
}

# Main function
# read urine bucket table, data already normalized relative to creatinine we do PQN in addition.
bt<-read.table("C:/Users/grw28475/Metabolomics/AKI/Test_wolfram_urine/buckettable/bucket_table",as.is=T, header=T,row.names=1)
bt1<-as.matrix(t(bt))
kick.out<-apply(bt1,1,function(z){all(z==0)})
bt2<-bt1[!kick.out,]
# exclude patients 47,54, 83 as diagnosis is unclear, Note Patient 42 already excluded as no sample available
bt3<-bt2[,-c(46,53,82)] # 106 samples with 701 buckets
#PQN normalization
reference <- apply(bt3,1,median)
quotient <- bt3/reference
quotient.median <- apply(quotient,2,median)
bt3.pqn <- t(t(bt3)/quotient.median)
bt4.pqn<-t(bt3.pqn) # transpose back for SVM
bt5.pqn<-bt4.pqn[1:105,] # remove last sample so that we can do leave 5 out cross-validation
#Labels
l1<-c(0,1,0,1,1,1,0,1,0,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,1,0,1,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,1,1,0,0,1,1,1,0,0,0,1,0,1,0,0,0,0,0,1,1,0,1,0,1,0,0,0,1,1,0,0,0,1,0,0,1,0,0,0,0) #labels for samples
l3<-c(0,1,0,1,1,1,0,1,0,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,1,0,1,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,1,1,0,0,1,1,1,0,0,0,1,0,1,0,0,0,0,0,1,1,0,1,0,1,0,0,0,1,1,0,0,0,1,0,0,1,0,0,0) #labels for 105 samples last sample removed for leave 5 out cross validation
#l2<-as.factor(l1) # for classification
l4<-as.factor(l3)
res<-double_nested_cv(bt5.pqn,l4,5,balanced_sampling=FALSE,5) #leave 5 out cross val with parameter optimization without balanced sampling max 5 features

