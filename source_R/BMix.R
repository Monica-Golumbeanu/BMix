##########################
# Wrapper for the learn_model function, applies the model on the input data
#
# Input: Signals_folder - path of the directory where the mismatch
#                         frequency profiles are located
#        Output_folder - path of the folder where the results are saved
#        CovTHR - minimum coverage to consider
#        PosteriorTHR - posterior lower bound used for classification
#        Separates_strands - 1 or 0, depending on whether the parameters should be learned for both strands separately or not
#
# Written by Monica Golumbeanu
# monica.golumbeanu@bsse.ethz.ch
##########################

library(plyr)
library(nloptr)

get_compressed_data_frame=function(DataArray, cov_thr) {
  mism_sum_vec = read.table(file = DataArray[1], row.names = NULL, col.names = c("chr","pos","X","Y"), skip = 2, sep="\t")
  mism_sum_vec = mism_sum_vec[mism_sum_vec$Y>=cov_thr,]
  compressed_vec = count(mism_sum_vec, c("X", "Y"))
  
  if (length(DataArray)>1) {
    for (i in 2:length(DataArray)) {
      mism_sum_vec = read.table(file = DataArray[i], row.names = NULL, col.names = c("chr","pos","X","Y"), skip = 2, sep="\t")
      mism_sum_vec = mism_sum_vec[mism_sum_vec$Y>=cov_thr,]
      compressed_vec_temp = count(mism_sum_vec, c("X", "Y"))
      compressed_vec_temp = rbind(compressed_vec_temp,compressed_vec)
      compressed_vec = aggregate(freq ~ X + Y, data = compressed_vec_temp, sum)
    }
  }
  return(compressed_vec)
}

c_binopdf = function(x, xc, p) {
# ~ binomial pdf
  p = (p^x) * ((1-p)^xc)
  return(p)
}

c_lbinopdf = function(x, xc, p) {
  # ~ binomial pdf
  p = log(p)*x + log(1-p)*xc
  return(p)
}

Calc_nLL = function(compressed_vec, epsilon, gamma, nu, psi){
# calculates the likelihood of the positive data
  epsl = log(2.2251E-308);
  theta = (1-gamma)*epsilon + (1-3*epsilon)*gamma;
  compressed_diff = compressed_vec$Y-compressed_vec$X;
  L = nu     * c_binopdf(compressed_vec$X, compressed_diff, 1-3*epsilon) + 
  (1-nu) * (1-psi) * c_binopdf(compressed_vec$X, compressed_diff, epsilon) + 
  (1-nu) *    psi  * c_binopdf(compressed_vec$X, compressed_diff, theta);
  nLL = -sum( compressed_vec$freq * sapply(log(L), function(x) max(epsl, x)));
  return(nLL)
}

NCalc_nLL = function(compressed_vec, epsilon, nu){
# calculates the likelihood of the negative data
  epsl = log(2.2251E-308);
  compressed_diff = compressed_vec$Y-compressed_vec$X;
  L = nu     * c_binopdf(compressed_vec$X, compressed_diff, 1-3*epsilon) + 
  (1-nu) * c_binopdf(compressed_vec$X, compressed_diff, epsilon)
  nLL = -sum( compressed_vec$freq * sapply(log(L), function(x) max(epsl, x)));
  return(nLL)
}

sGet_like_grad=function(compressed_vec, epsilon, gamma, nu, psi) {
# gradient of the positive data likelihood
  nLL =  Calc_nLL(compressed_vec, epsilon, gamma, nu, psi);
  EPS=1E-5;
  GnLL = c( Calc_nLL(compressed_vec, (epsilon+EPS), gamma, nu, psi), 
            Calc_nLL(compressed_vec, epsilon, (gamma+EPS), nu, psi), 
            Calc_nLL(compressed_vec, epsilon, gamma, (nu+EPS), psi), 
            Calc_nLL(compressed_vec, epsilon, gamma, nu, (psi+EPS)))
  GnLL = (GnLL-nLL)/EPS;
  return(GnLL)
}

sGet_likeN_grad=function(compressed_vec, epsilon, nu) {
# gradient of the negative data likelihood
  nLL =  NCalc_nLL(compressed_vec, epsilon, nu);
  EPS=1E-5;
  GnLL = c( NCalc_nLL(compressed_vec, (epsilon+EPS), nu),  
            nLL,
            NCalc_nLL(compressed_vec, epsilon, (nu+EPS)),
            nLL)
  GnLL = (GnLL-nLL)/EPS;
  return(GnLL)
}

Get_total_like = function(z, compressed_vec, compressed_neg_vec) {
# total likelihood
  LL  =  Calc_nLL(compressed_vec, z[1], z[2], z[3], z[4]);
  if (length(compressed_neg_vec) > 0) {
    nLL =  NCalc_nLL(compressed_neg_vec, z[1], z[3]);
  } else {
    nLL = 0;
  }
  return(LL+nLL)
}

Get_total_like_grad = function(z, compressed_vec, compressed_neg_vec) {
# gradient of the total likelihood
  grad1 =  sGet_like_grad(compressed_vec, z[1], z[2], z[3], z[4]);
  if (length(compressed_neg_vec) > 0) {
    grad2 =  sGet_likeN_grad(compressed_neg_vec, z[1], z[3]);
  } else {
    grad2 = 0;
  }
  GnLL =  grad1 + grad2;
  return(GnLL)
}

sClassify = function(X,Y, epsilon, gamma, nu, psi) {
  theta = (1-gamma)*epsilon + (1-3*epsilon)*gamma;
  Xc = Y-X; 
  # Calculate posteriors
  lpC = cbind(log(nu)+c_lbinopdf(X, Xc, 1-3*epsilon), 
              log(1-nu) + log(1-psi) + c_lbinopdf(X, Xc, epsilon),
              log(1-nu)+log(psi)+c_lbinopdf(X, Xc, theta))
  
  max_lpC = apply(lpC,MARGIN = 1,function(x) max(x))
  lpC = lpC - max_lpC
  pC = exp(lpC)
  pC = pC / apply(pC,MARGIN = 1, sum)
  MAP = apply(pC,MARGIN = 1, function(x) max(x))
  MAP_I = apply(pC,MARGIN = 1,function(x) which(x == max(x)))
  return(cbind(MAP_I,MAP))
}

plot_classification = function (figure_file_name, classified_data, PosteriorTHR){ 
  # plot the classified points
  pdf(width=5, height=5, figure_file_name, useDingbats=F)
  F1 = (classified_data$MLE_ClassI==1)
  F2 = (classified_data$MLE_ClassI==2)
  F3 = (classified_data$MLE_ClassI==3)
  F4 = (classified_data$MLE_pC<PosteriorTHR)
  xl = range(c(classified_data$Y[F1],classified_data$Y[F2],classified_data$Y[F3]))
  yl = c(0,1)
  plot(classified_data$Y[F1],classified_data$X[F1]/classified_data$Y[F1], 
       log="x", xlim=xl, ylim=yl, pch=20, col = "darkorange", xaxt = 'n', yaxt='n',
       xlab="",ylab="")
  par(new=TRUE)
  plot(classified_data$Y[F2],classified_data$X[F2]/classified_data$Y[F2], 
       log="x", xlim=xl, ylim=yl, pch=20, col = "darkred", xaxt = 'n', yaxt='n',
       xlab="",ylab="")
  par(new=TRUE)
  plot(classified_data$Y[F3],classified_data$X[F3]/classified_data$Y[F3], 
       log="x", xlim=xl, ylim=yl, pch=20, col = "navy", xaxt = 'n', yaxt='n',
       xlab="",ylab="")
  par(new=TRUE)
  plot(classified_data$Y[F4],classified_data$X[F4]/classified_data$Y[F4], 
       log="x", xlim=xl, ylim=yl, pch=20, col = "darkgrey",
       xlab="coverage",ylab="mismatch frequency", cex.axis=1.5, cex.lab=1.5)
  legend("topright", cex = 0.5, legend = c("variant", "background", "cross-link loci", "not classified"), bty = "n", lty = 1, lwd=4, col = c("darkorange","darkred","navy","darkgrey"))
  par(new=FALSE)
  dev.off()
}

learn_model = function(DataPathArray, NegDataPathArray, OutDir, CovTHR, PosteriorTHR) { 
  cPosData = get_compressed_data_frame(DataPathArray, CovTHR)
  if (length(NegDataPathArray) > 0) {
    cNegData = get_compressed_data_frame(NegDataPathArray, CovTHR)
  } else {
    cNegData = NULL;
  }
  
  # create output and figure folders
  dir.create(OutDir, showWarnings = FALSE)
  dir.create(paste(OutDir,"Figures/", sep=""), showWarnings = FALSE)
  
  # learn parameters
  params0 = c(0.1,0.1,0.1,0.1)
  opts <- list("algorithm"="NLOPT_LD_MMA","xtol_rel"=1.0e-5, "xtol_abs"=1.0e-5)
  res <- nloptr( x0=params0, eval_f=Get_total_like, eval_grad_f=Get_total_like_grad, lb = c(1E-3, 1E-3, 1E-6, 1E-3), ub = c(0.25, 1, 1, 1), opts=opts, compressed_vec = cPosData, compressed_neg_vec = cNegData)

  print("Estimated parameters:")
  print(paste("epsilon", res$solution[1], sep=" "))
  print(paste("gamma", res$solution[2], sep=" "))
  print(paste("nu", res$solution[3], sep=" "))
  print(paste("psi", res$solution[4], sep=" "))
  parameters_frame = data.frame(res$solution,row.names = c("epsilon","gamma","nu","psi"))
  write.table(parameters_frame, file = paste(OutDir, "parameters.txt"), row.names=TRUE, col.names=FALSE,append = TRUE, sep="\t", quote = FALSE)
  
  # classify
  classification = data.frame(cbind(cPosData$X, cPosData$Y, sClassify(cPosData$X, cPosData$Y, res$solution[1], res$solution[2], res$solution[3], res$solution[4])))
  colnames(classification) = c("X", "Y", "MLE_ClassI","MLE_pC")
  # write the results to file
  for (i in 1:length(DataPathArray)) {
    mism_sum_vec = read.table(file = DataPathArray[i], row.names = NULL, col.names = c("Chr","Position","X","Y"), skip = 2, sep="\t")
    mism_sum_vec = mism_sum_vec[mism_sum_vec$Y>=CovTHR,]
    cDATA = data.frame(mism_sum_vec,id... = seq_len(nrow(mism_sum_vec)))
    total = merge(classification,cDATA, by = c("X","Y"))
    ordered_data = order(total$id...)
    good_columns = colnames(total) != "id..."
    classified_data = total[ordered_data, c("Chr","Position","X","Y","MLE_ClassI","MLE_pC")]
    write.table(classified_data,file = paste(OutDir,strsplit(basename(DataPathArray[i]),"txt"),"results",sep="") ,row.names = FALSE,col.names = TRUE,sep="\t", quote = FALSE)
    figure_file_name = paste(OutDir,"Figures/Classification_", strsplit(basename(DataPathArray[i]),"txt"),"pdf", sep="")
    plot_classification(figure_file_name, classified_data, PosteriorTHR)
  }
}

args = commandArgs(TRUE)
print(args)
Signals_folder = args[1];
Output_folder = args[2];
CovTHR = as.numeric(args[3]);
PosteriorTHR = as.numeric(args[4]);
Separate_strands = as.numeric(args[5]);

if (Separate_strands == 1) {
  # learn parameters for the forward strand
  DataPathArray = paste(Signals_folder, "TC_f.txt",sep="")
  NegDataPathArray = c(paste(Signals_folder, "AC_f.txt",sep=""),paste(Signals_folder, "GC_f.txt",sep=""))
  learn_model(DataPathArray, NegDataPathArray, Output_folder, CovTHR, PosteriorTHR);
  
  # learn parameters for the reverse strand
  DataPathArray = paste(Signals_folder, "AG_r.txt",sep="")
  NegDataPathArray = c(paste(Signals_folder, "TG_r.txt", sep=""),paste(Signals_folder, "CG_r.txt",sep=""))
  learn_model(DataPathArray, NegDataPathArray, Output_folder, CovTHR, PosteriorTHR);
} else { 
  DataPathArray = c(paste(Signals_folder, "TC_f.txt",sep=""), paste(Signals_folder, "AG_r.txt",sep=""))
  NegDataPathArray = c(paste(Signals_folder, "AC_f.txt",sep=""), paste(Signals_folder, "GC_f.txt",sep=""), 
                       paste(Signals_folder, "TG_r.txt",sep=""), paste(Signals_folder, "CG_r.txt",sep=""))
  learn_model(DataPathArray, NegDataPathArray, Output_folder, CovTHR, PosteriorTHR);
}


