
### Brownian Bridge transformation
BBT <-
function (x){

(cumsum(x) - seq_along(x)/length(x) * sum(x))/sqrt(length(x) * 
    var(x))
}
### get break points form PvarBreakTest
BreakPoints <-
function (object) 
object$BreakPoints

### get the mean of p-variation acording of BBT(x), then H0 is TRUE (according to n).
getMean <-
function (n, bMean = MeanCoef) 
{
    unname(bMean[1] + bMean[2] * n^bMean[3])
}

### get the standart deviation of p-variation acording of BBT(x), then H0 is TRUE (according to n).
getSd <-
function (n, bSd = SdCoef) 
{
    unname(bSd[1] + bSd[2] * n^bSd[3])
}

### Normalyse p-variation acording n.
NormalisePvar <-
function (x, n=length(x), bMean = MeanCoef, bSd = SdCoef) 
{
    (x - getMean(n, bMean))/getSd(n, bSd)
}

### the test of structial breaks by p-variation.
PvarBreakTest <-
function (x, TimeLabel = as.vector(time(x)), alpha = 0.05, FullInfo = TRUE) 
{
    StartTime <- proc.time()
    dname <- deparse(substitute(x))
    NAInd = is.na(x)
    if (any(NAInd)) {
        warning("NA values was removed.")
        x = x[!NAInd]
        TimeLabel = TimeLabel[!NAInd]
    }
    n = length(x)
    if (n < 20) 
        stop("The size must be greater than 20.")
    if (n < 100) 
        warning("The test migth by biased when n<100.")
    CriticalValue = PvarQuantile(n, prob = 1 - alpha)
	
	if( sd(x)==0 ){
		y = x
	}else{
		y = BBT(x)
    }
	
	PvarStat = pvar(y, p = 4, TimeLabel = TimeLabel, info = FullInfo)
    if (FullInfo) {
        Stat = unname(PvarStat$value)
        p.value = PvarPvalue(n, Stat)
        reject = Stat >= CriticalValue[1]
        if (reject) {
            BreakPos = abs(PvarStat$Partition/n - 0.5)
            accept = abs(PvarStat$Partition/n - 0.5) < 0.4
            if (!any(accept)) {
                accept[which.min[BreakPos]] = TRUE
            }
            BreakPoints = PvarStat$Partition[accept]
        }
        else {
            BreakPoints = NULL
        }
        ans = list(
		    Stat = c(Statistics = Stat)
			, CriticalValue = c(`Critical Value` = CriticalValue)
			, alpha = c(alpha = alpha)
			, p.value = c(`~p.value` = p.value)
			, reject = reject
			, dname = dname
			, p = unname(PvarStat$p)
			, x = x
			, y = y
			, TimeLabel = TimeLabel
			, BreakPoints = BreakPoints
			, Partition = PvarStat$Partition
			, SplitPoints = PvarStat$SplitPoints
			, info = PvarStat$info
		)
        FinishTime <- proc.time()
        TakesTime = unname((FinishTime - StartTime)[3])
        ans$info$TakesTime = TakesTime
        class(ans) <- "PvarBreakTest"
        return(ans)
    }
    else {
        Stat = PvarStat
        p.value = PvarPvalue(n, Stat)
        ans = c(Statistics = Stat, `Critical Value` = CriticalValue, 
            alpha = alpha, `~p.value` = p.value)
        return(ans)
    }
}

### get critical values accorging to n and prob.
PvarQuantile <-
function (n, prob = c(0.9, 0.95, 0.99), DF = PvarQuantileDF) 
{
    intervals = cut(prob, breaks = DF$prob, include.lowest = TRUE, 
        right = TRUE)
    Quant = DF$Quant[as.numeric(intervals) + 1]
    if (any(is.na(Quant))) 
        warning("`prob` must be between 0 and 1")
    ans = Quant * getSd(n) + getMean(n)
    unname(ans)
}

### gets approcimate p-value, of the test.
PvarPvalue <-
function (n, stat, DF = PvarQuantileDF) 
{
    NormStat = NormalisePvar(stat, n)
    intervals = cut(NormStat, breaks = c(-Inf, DF$Quant, Inf), include.lowest = TRUE, 
        right = FALSE)
    ind = as.numeric(intervals) - 1
    ind[ind == 0] = 1
    1 - DF$prob[ind]
}


#######################################################################################
############################# PvarBreakTest S3 methods ################################
#######################################################################################


plot.PvarBreakTest <-
function (x, main1 = "Data", main2 = "Brownian Bridge transformation", 
    ylab1 = x$dname, ylab2 = "BBT(" %.% x$dname %.% ")", sub2=NULL, 
	col.PP = 3, cex.PP = 0.5, col.BP = 2, cex.BP = 1, cex.DP = 0.5, ...) 
{
    Time = x$TimeLabel
    op <- par(mfrow = c(2, 1), mar = c(5.1, 4.1, 2.1, 2.1))
    
	plot(Time, x$x, type = "p", pch = 19, cex = cex.DP, ylab = ylab1, 
        main = main1, ...)
    BP = c(1, x$BreakPoints, length(x$x))
    for (i in 2:length(BP)) {
        Eseg = mean(x$x[BP[i - 1]:BP[i]])
        colseg = i%%2 + 2
        segments(x0 = Time[BP[i - 1]], y0 = Eseg, x1 = Time[BP[i]], y1 = Eseg, 
            lwd = 3, col = colseg)
    }
	
	if(is.null(sub2)){	
		if(x$reject){
			if(length(x$BreakPoints)==1){
				sub2 = "Program suggests " %.% length(x$BreakPoints) %.% " break point."
			}else{
				sub2 = "Program suggests " %.% length(x$BreakPoints) %.% " break points."
			}	
		}else{
			sub2 = "Program didn't find structural breaks at the confidens level of " %.% formatC(head(x$alpha, 1)) %.% "."
		}
	}
	
    plot(Time, x$y, type = "l", ylab = ylab2, main = main2, sub=sub2, ...)
    points(x$TimeLabel[x$Partition], x$y[x$Partition], cex = cex.PP, 
        pch = 19, col = col.PP, bg = col.PP)
    points(x$TimeLabel[x$BreakPoints], x$y[x$BreakPoints], cex = cex.BP, 
        pch = 19, col = col.BP, bg = col.BP)
    par(op)
}

print.PvarBreakTest <-
function (x, ...) 
{
    cat("       PvarBreakTest \n\n")
    cat("H0: there is no structural break. \n")
    cat("Results: ")
    if (x$reject) {
        cat("H0 is rejected at the confidens level of " %.% formatC(head(x$alpha, 
            1)) %.% ".\n")
    }
    else {
        cat("H0 is accepted at the confidens level of " %.% formatC(head(x$alpha, 
            1)) %.% ".\n")
    }
    cat("Data: " %.% x$dname %.% ", n=" %.% length(x$x) %.% ".\n")
    cat("Test's output: \n")
    print(c(x$Stat, x$CriticalValue, x$alpha, x$p.value))
}

summary.PvarBreakTest <-
function (object, ...) 
{
    class(object) <- c("summary.PvarBreakTest", "PvarBreakTest")
    object
}

print.summary.PvarBreakTest <-
function (x, ...) 
{
    cat("The summary of PvarBreakTest:\n")
    cat("H0: there is no structural break. \n")
    cat("Results: ")
    if (x$reject) {
        cat("H0 is rejected at the confidens level of " %.% formatC(head(x$alpha, 
            1)) %.% ".\n")
    }
    else {
        cat("H0 is accepted at the confidens level of " %.% formatC(head(x$alpha, 
            1)) %.% ".\n")
    }
    if (x$reject) {
        if (length(x$BreakPoints) > 6) {
            cat("Suggesting " %.% length(x$BreakPoints) %.% " break points: " %.% 
                paste(formatC(head(x$BreakPoints, 6)), collapse = ", ") %.% 
                ", ...\n")
        }
        else {
            if (length(x$BreakPoints) == 1) {
                cat("Suggesting " %.% length(x$BreakPoints) %.% 
                  " break point: " %.% paste(formatC(head(x$BreakPoints, 
                  6)), collapse = ", ") %.% ".\n")
            }
            else {
                cat("Suggesting " %.% length(x$BreakPoints) %.% 
                  " break points: " %.% paste(formatC(head(x$BreakPoints, 
                  6)), collapse = ", ") %.% ".\n")
            }
        }
    }
    cat("Data: " %.% x$dname %.% ", n=" %.% length(x$x) %.% ".\n")
    cat("Test's output: \n")
    print(c(x$Stat, x$CriticalValue, x$alpha, x$p.value))
    cat("Takes time: " %.% formatC(x$info$TakesTime) %.% "\n")
    cat("p-avriation calculation info:\n")
    print(unlist(x$info)[-length(x$info)])
    if (length(x$x) > 6) {
        cat("\nData vector (n=" %.% length(x$x) %.% "): " %.% 
            paste(formatC(head(x$x, 6)), collapse = ", ") %.% 
            ", ...\n")
    }
    else {
        cat("\nData vector (n=" %.% length(x$x) %.% "): " %.% 
            paste(formatC(head(x$x, 6)), collapse = ", ") %.% 
            ".\n")
    }
}




