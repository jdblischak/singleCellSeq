# Function from Ding et al. 2015.
# Downloaded GRM-0.2.1.tgz on 2015-06-29.
# Download site: http://wanglab.ucsd.edu/star/GRM/
# publication: http://bioinformatics.oxfordjournals.org/content/early/2015/03/22/bioinformatics.btv122.full
library("MASS") # for function gamma.dispersion

### 1.gamma regression and make prediction
# ercc_response file: the fpkm reads for ERCC across all the samples
#             format: each row is each ERCC, each column is each sample
#                     especially sorted by the ERCC names
# ercc file: the add-in standard ERCC concentration, gotten from the experiments directly
#    format: each row is each ERCC, each column is each sample
#            the first column is ERCC name, the second column is ERCC standard concentration.
#            especially sorted by the ERCC names
# gene file: the fpkm reads for genes across all the sampels
#    format: each row is each gene, each column is each sample
#            especially the column names should be in the same order with ercc_response file
# filename: the vector for all the samples' names
#    format: it should be a 1*n vector in character variable type

gammareg=function(ercc_response, ercc, gene, filename,SE=FALSE){
  ## check files
  # check filename in character type
  filename=as.matrix(filename)
  filename=as.character(filename)
  # check ercc_response and ercc share the same ercc
  spike=ercc[,1]
  spike=as.character(spike)
  controlset=ercc[spike %in% rownames(ercc_response),]
  control=controlset[,2]
  logctr=log(control)
  gene_name=rownames(gene)
  name=append(spike,gene_name)
  total=dim(ercc_response)[1]+dim(gene)[1]
  all_predictresult=matrix(0, total, dim(ercc_response)[2])
  all_error=matrix(0, total, dim(ercc_response)[2])
  all_predictresult=data.frame(all_predictresult)
  all_error=data.frame(all_error)
  sename=c()
  pvalue=c()
  i=1
  # do regression sample by sample
  for (i in 1: dim(ercc_response)[2]){

    ##matrix preparation to input
    response=ercc_response[,i]
    genelist=gene[,i]
    set=cbind(control, logctr, response)
    set=data.frame(set)
    subset=set[which(response>0),]
    logresponse=log(subset$response)
    subset=cbind(subset, logresponse)
    subset1=subset[which(subset$logctr>0),]

    ##gamma regression

    # give a reasonable initial points
    nullmodel=glm(control~1, family=Gamma(link="log"), data=subset1)
    # regression
    fit1=glm(control~I(logresponse)+I(logresponse^2), family=Gamma(link="log"), data=subset1, start=c(coef(nullmodel)[1],0,0))

    #Log likelihood ratio test for the parameters
    LRT=anova(nullmodel, fit1, test="Chisq")
    p=LRT$Pr[2]
    pvalue=c(pvalue, p)

    ## make prediction for each variable in one sample

    # make predict dataset with index name
    predictlist=append(response, genelist)
    predictlist=as.vector(predictlist)
    pred1=c()
    nonzeroindex=c()
    # looking for nonzero value
    j=1
    for ( j in 1: length(predictlist)) {
      if(predictlist[j]!=0){
        pred1=c(pred1,predictlist[j])
        nonzeroindex=c(nonzeroindex, j)
      }
    }

    # fpkm expression expression equals to non-zero
    pred1=as.vector(pred1)
    logset_pred=log(pred1)

    # make prediction
    disper=gamma.dispersion(fit1)
    preddata=data.frame(logresponse=logset_pred)
    pred=predict.glm(fit1,newdata=preddata,se.fit=TRUE, dispersion=disper)
    p=pred$fit
    pred_value1=exp(p)

    # if standard error(log scale) needed
    if (SE){
      seresult=rep(10000,total);
      seresult=data.frame(seresult);

      # standard error calculated from predict.glm function
      # if the original fpkm is zero, then we set the standard error
      # for these genes or ERCCs are 10000
      s1=pred$se.fit;

      # standard error output
      k=1;
      for (k in 1: length(nonzeroindex)){
        seresult[nonzeroindex[k],1]=s1[[k]]
      }
      seresult=data.frame(seresult);
      colname=filename[i];
      secolname=paste(colname,"se",sep="_");
      sename=c(sename, secolname);
      all_error[,i]=seresult;
    }
    # predict result as original expression scale for each item in one sample
    nullpredict=rep(0,total)
    nullpredict=data.frame(nullpredict)
    q=1
    for (q in 1: length(nonzeroindex)){
      nullpredict[nonzeroindex[q],1]=pred_value1[[q]]
    }

    predict_result=nullpredict
    predict_result=as.matrix(predict_result)


    #  predict result for all the samples
    all_predictresult[, i]=predict_result

    print(i)
  }
  ## result output

  # predict value
  all_predictresult=data.frame(all_predictresult)
  all_error=data.frame(all_error)
  rownames(all_predictresult)=name
  colnames(all_predictresult)=filename

  # pvalue for each model fit in each sample
  pvalue=data.frame(pvalue)
  rownames(pvalue)=filename

  if (SE){
    all_error=data.frame(all_error)
    rownames(all_error)=name
    colnames(all_error)=sename
    return(list(errorset=all_error,predictset=all_predictresult, pvalue_set=pvalue))
  }else{
    return(list(predictset=all_predictresult,pvalue_set=pvalue))
  }
}
