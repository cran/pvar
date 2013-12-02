
#######################################################################################
############################### pvar public functions #################################
#######################################################################################

### finds Sp sum.
Sum_p <-
function (x, p, lag = 1) 
{
    sum((abs(diff(x, lag)))^p)
}


### The main p-variation calculation function
pvar <-
function (x, p, TimeLabel=as.vector(time(x)), info = TRUE, sizeN = 7) 
{

	StartTime <- proc.time()
	dname <- deparse(substitute(x))
	
	NAInd = is.na(x)
    if (any(NAInd)) {
        warning("NA values was removed.")
        x = x[!NAInd]
        TimeLabel = TimeLabel[!NAInd]
    }
	
	sizeN = as.integer(sizeN)
	if(sizeN<1){
		warning("sizeN must be positive odd number. sizeN changed to 7")
		sizeN = 7
	}

	if(sizeN%%2==0){
		sizeN = sizeN-1
		warning("sizeN must be positive odd number. sizeN changed to "%.%sizeN)
	}

	### Check posible errors
	
	if (length(p) != 1) {
		if(length(p)<1) stop("The length of 'p' is zero.")
		warning("The 'p' must be a scalar. Only first element is used.")
		p = p[1]
    }	
	
    if (p <= 0) {
        stop("'p' must be positive.")
    }
	
	if (length(x) < 2) {
		if(length(x) < 1)
			stop("The length of 'x' is zero.")
			
        pvar.value = 0
		FinishTime <- proc.time()
		TakesTime = unname((FinishTime-StartTime)[3])
		if (info) {
			SplitPoints = c(1,length(x))
			if(p==1){
				MP = ChangePointsId(x)
			}else{
				MP = seq_along(x)
			}
			Partition = MP
			sizeN = length(x)		
			ans = list(
				value = c("p-variation" = pvar.value)
				,x = x
				,p = p
				,dname = dname
				,TimeLabel = TimeLabel
				,Partition = Partition
				,SplitPoints = SplitPoints
				,info = list(					
					PartitionN = length(Partition)
					,AnalysedPointsN = length(MP)
					,SegmentsN = length(SplitPoints)-1
					,sizeN = sizeN
					,TakesTime = TakesTime
				)
			)		
			class(ans) <- "pvar"
			return(ans)
		}
		else {
			return(pvar.value)
		}
    }
	
	### check taime label
	if(length(TimeLabel)!=length(x)){
		TimeLabel = seq_along(x)
		warning("TimeLabel must have the same length as 'x'. TimeLabel changed to `seq_along(x)`")
	}
	
	
	### if p<=1, then maximums is spliting x into as small as posible.
    if (p <= 1) {
		
		pvar.value = Sum_p(x, p)
		FinishTime <- proc.time()
		TakesTime = unname((FinishTime-StartTime)[3])
		if (info) {
			#SplitPoints = c(1,length(x))
			MP = ChangePointsId(x)
			SplitPoints = MP[SplitByExtremum(x[MP])]
			if(p==1){
				MP = ChangePointsId(x)
			}else{
				MP = seq_along(x)
			}
			Partition = MP
			sizeN = length(x)		
			ans = list(
				value = c("p-variation" = pvar.value)
				,x = x
				,p = p
				,dname = dname
				,TimeLabel = TimeLabel
				,Partition = Partition
				,SplitPoints = SplitPoints
				,info = list(					
					PartitionN = length(Partition)
					,AnalysedPointsN = length(MP)
					,SegmentsN = length(SplitPoints)-1
					,sizeN = sizeN
					,TakesTime = TakesTime
				)
			)		
			class(ans) <- "pvar"
			return(ans)
		}
		else {
			return(pvar.value)
		}
    }
	
	
	### drom meaningless points: 1) drop monotonic points; 2) Drop midle unimportent points.
    MP = ChangePointsId(x)
    if (info) {
        MP = MP[DropMidAllB(x[MP], p, sizeN)]
        xx = x[MP]
    }
    else {
        xx = DropMidAll(x[MP], p, sizeN)		# the number "7" mus be the same as in function 'pvarPseudoMonotonic'
    }
	
	### Split vector into Pseudo-Monotonic segments
    Len = length(xx)
    # s = Split_MinMax_req(xx, add = 1)
    s = SplitByExtremum(xx)
    Len_s = length(s)
    ChB = list(TRUE)
    pvarVec = list(0)
	
	### calculate p-var for each segment
    for (i in 2:(Len_s)) {
        pvarQ = pvarPseudoMonotonic(x = xx, p = p, a = s[i - 1], b = s[i], 
            sizeN = sizeN)
        ChB[[i]] = pvarQ$ChB
        pvarVec[[i]] = pvarQ$pvarVec
    }
    ChB = unlist(ChB)
    pvarVec = unlist(pvarVec)
    
	### the main anser:
	pvar.value = Sum_p(xx[ChB], p)
	
	FinishTime <- proc.time()
	TakesTime = unname((FinishTime-StartTime)[3])
	
	### Check possible error:
    if (abs(sum(pvarVec[ChB]) - Sum_p(xx[ChB], p)) > 1/10^13) {
        warning("Posible error(001),'pvar': The Sum_p is not the same as result form pvar.")
    }
    if (info) {
        SplitPoints = MP[s]
        Partition = MP[ChB]
        if (abs(Sum_p(xx[ChB], p) - Sum_p(x[Partition], p)) != 0) {
            Diference = (Sum_p(xx[ChB], p) - Sum_p(x[Partition], p))
            warning("Posible error(002),'pvar': The Sum_p is not the same in pvar points and global points, the diff: " %.% 
                Diference)
        }
        ans = list(
			value = c("p-variation" = pvar.value)
			,x = x
			,p = p
			,dname = dname
			,TimeLabel = TimeLabel
			,Partition = Partition
			,SplitPoints = SplitPoints
			,info = list(
				PartitionN = length(Partition)
				,AnalysedPointsN = length(MP)
				,SegmentsN = length(SplitPoints)-1
				,sizeN = sizeN
				,TakesTime = TakesTime
			)
		)		
        class(ans) <- "pvar"
        return(ans)
    }
    else {
        return(pvar.value)
    }
}



### spit vector according to min and max.
Split_MinMax <-
function (x) 
{
    sort(unique(c(1, which.max(x), which.min(x), length(x))))
}

### spit vector recursively according to min and max.
Split_MinMax_req <-
function (x, add = 1) 
{
    s = Split_MinMax(x)
    Len = length(s)
    resList = list()
    if (Len > 2) {
        for (i in 1:(Len - 1)) {
            resList[[i]] = Recall(x[(s[i]:s[i + 1])], s[i])
        }
        res = sort(unique(unlist(resList))) + add - 1
        return(res)
    }
    else {
        return(s + add - 1)
    }
}



SplitByExtremumFromStart <- function(x, a=1, b=length(x), aMin = TRUE){
	# a = which.min(x)
	# b = length(x)
	# aMin = TRUE

#	ExtremVec = rep(NA, b-a+1)
	ExtremVec = rep(NA, 10)		
	#ExtremVec[1] = a #pradinio tasko neitraukinejam, kad nereiketu ismetineti
	
	i=1
	while(a<b & i<=length(x) ){
		i = i + 1
		
		if(length(ExtremVec)<i)
			ExtremVec = append(ExtremVec, rep(NA, length(ExtremVec)))
		
		if(aMin){
			aNew = which.max(x[a:b]) + a - 1
			aMin = FALSE
		}else{
			aNew = which.min(x[a:b]) + a - 1
			aMin = TRUE
		}
		
		if(aNew==a)
			aNew = b
			
		a = aNew
		ExtremVec[i] = a 	
	}
	
	ExtremVec[!is.na(ExtremVec)]
	
}


SplitByExtremum <- function(x){	
	if(length(x)==0)
		return(vector("numeric"))
	MinID = which(x==min(x))
	MaxID = which(x==max(x))
	SplitPoints = c(1, MinID, MaxID, length(x))
	SplitPointsDup = duplicated(SplitPoints)
	
	SplitPoints = SplitPoints[!SplitPointsDup]
	SplitPointsOrder = order(SplitPoints)
	SplitPoints = SplitPoints[SplitPointsOrder]
	
	FunctInd = c(NA, rep(TRUE,length(MinID)), rep(FALSE,length(MaxID)), NA)	# jei TRUE, tai minimumo taskas
	FunctInd = FunctInd[!SplitPointsDup][SplitPointsOrder]
	
	if(length(SplitPoints)==1){
		res = list(1)
	}else{
		res = list(rev(1 + SplitPoints[2] - SplitByExtremumFromStart(rev(x[SplitPoints[1]:SplitPoints[2]]), a=1, b=SplitPoints[2], aMin = FunctInd[2])), SplitPoints[2])
		if (length(SplitPoints)>2)
			res = c(res,sapply(2:(length(SplitPoints)-1), function(i) SplitByExtremumFromStart(x, a=SplitPoints[i], b=SplitPoints[i+1], aMin = FunctInd[i])))
	}
	
	res = unlist(res)	
	res
	
}	




### Finds the ids of change points in vector
ChangePointsId <-
function (x) 
{
    Len = length(x)
    id = seq_along(x)
    x = sign(diff(x))
    id = id[x != 0]
    x = diff(x[id])
    c(1, id[which(x != 0) + 1], Len)
}
# ChangePointsId(c(0,0,0))

#######################################################################################
################################## pvar S3 methods ####################################
#######################################################################################
print.pvar <- 
function(x, ...)
{
	print(x$value)
}


print.pvar <- 
function(x, ...)
{
	print(x$value)
}

summary.pvar <- 
function(object, ...){
	class(object) <- c("summary.pvar", "pvar")
	object
}

print.summary.pvar <- 
function(x, ...){
	cat("The summary of p-variation:\n")
	cat("Value: "%.%formatC(x$value)%.%", p = "%.%x$p%.%"\n")
	cat("Data: " %.% x$dname %.% ", n=" %.% length(x$x) %.% "\n")
	cat("Takes time: "%.%formatC(x$info$TakesTime)%.%"\n")
	cat("Info:\n")
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
	
	if (length(x$Partition) > 6) {
        cat("With "%.%length(x$Partition)%.%" meaningful points: "%.% 
            paste(formatC(head(x$Partition, 6)), collapse = ", ") %.% 
            ", ...\n")
    }
    else {
        cat("With "%.%length(x$Partition)%.%" meaningful points: "%.% 
            paste(formatC(head(x$Partition, 6)), collapse = ", ") %.% 
            ".\n")
    }
	
} 


plot.pvar <- 
function(x, 
	main = "p-variation",
	ylab = x$dname, 
	sub="p="%.%round(x$p,5)%.%", p-variation: "%.% formatC(x$value, 5, format = "f"), 
	col.PP = 3, cex.PP = 0.5, col.SP = 2, cex.SP = 1,	...)
{
	Time = x$TimeLabel
	plot(Time,  x$x, type="l", sub=sub, ylab=ylab, main=main,...)
	
	points(x$TimeLabel[x$Partition],x$x[x$Partition], cex=cex.PP, pch=19, col=col.PP, bg=col.PP)
	points(x$TimeLabel[x$SplitPoints], x$x[x$SplitPoints], cex=cex.SP, pch=19, col=col.SP, bg=col.SP)

}

# `==.pvar` <-
# function (x, y){ 
	# x$value==y$value
# }

Ops.pvar <-
function(e1, e2) 
{
    if (nargs() == 1) 
        stop("unary ", .Generic, " not defined for pvar objects")
		
    boolean <- switch(.Generic, `<` = , `>` = , `==` = , `!=` = , `<=` = , `>=` = TRUE, FALSE)
    if (boolean) 
		return(eval(call(.Generic, unname(e1$value), unname(e2$value))))
		
	# if(.Generic=="+")
		# return(AddPvar(e1, e2))

    stop(.Generic, " not defined for pvar objects")
}


# sum.pvar <- 
# function(..., na.rm = FALSE, dname = "data"){
	# ArgL = list(...)
	
	# if (na.rm){
		# for(i in length(ArgL):1){
			# if(is.na(ArgL[[i]][1]))
				# ArgL[[i]] = NULL			
		# }
	# }

	# do.call("AddPvar", c(ArgL, JoinIfPossible = TRUE, dname = "data"))
# }


######################################################################################
######################################################################################
############################ The end of public functions #############################
######################################################################################
######################################################################################




#######################################################################################
########################### General helping Functions #################################
#######################################################################################

### add two strings
`%.%` <-
function (x, y) 
paste(x, y, sep = "")

### The sum of moving window
MovingSum <-
function (x, n = 2) 
filter(x, rep(1, n), sides = 1)[-seq_len(n - 1)]

#######################################################################################
############################## Generating functions ###################################
#######################################################################################

### Brownian bridege
rbridge <-
function (end = 1, frequency = 1000) 
{
    z <- rwiener(end = end, frequency = frequency)
    ts(z - time(z) * as.vector(z)[frequency], start = 0, frequency = frequency)
}

### Wiener proces
rwiener <-
function (end = 1, frequency = 1000) 
{
    z <- c(0, cumsum(rnorm(end * frequency)/sqrt(frequency)))
     ts(z, start = 0, end=1, frequency = frequency )
}

#######################################################################################
############################## Iner pvar functions ###################################
#######################################################################################


### Similar functions to Sum_p:
Sum_pt <-
function (x, p, a, b) 
{
    (abs(x[a] - x[b]))^p
}
vec_p <-
function (x, p, lag = 1) 
{
    (abs(diff(x, lag)))^p
}
CumSum_p <-
function (x, p, lag = 1) 
{
    cumsum((abs(diff(x, lag)))^p)
}

### the function of p-variation, but only to pseudo-monotonic quasy monotonic functions.  
pvarPseudoMonotonic <-
function (x, p, a = 1, b = length(x), sizeN = 1) 
{
	### if the is little points, we now that they all menaingfull:
    if (abs(b - a) <= sizeN) {	
        pvarVec = vec_p(x[a:b], p)
        ans = list(ChB = rep(TRUE, b - a), pvarVec = pvarVec)
        return(ans)
    }
	
	### find analysing points
    APr = FindAnalysingPointsOP(x[a:b])
    AP = APr + a - 1
    
	### for each segment splint into 2 more and calcuate recursevly as ir become small enough.
	ChB = list()
    pvarVec = list()
    for (i in 1:(length(AP) - 2)) {
        SAP = Split_MinMax(x[AP[i]:AP[i + 1]]) + AP[i] - 1
        pvarSeg1 = Recall(x, p, SAP[1], SAP[2])
        pvarSeg2 = Recall(x, p, SAP[2], SAP[3])
        ChB[[i]] = c(pvarSeg1$ChB, pvarSeg2$ChB)
        pvarVec[[i]] = c(pvarSeg1$pvarVec, pvarSeg2$pvarVec)
    }
    
	### add the last point:
	ChB[[i + 1]] = TRUE
    pvarVec[[i + 1]] = Sum_pt(x, p, AP[i + 1], AP[i + 2])
 
	ChB = unlist(ChB)
    pvarVec = unlist(pvarVec)
   
	### join intervals, by examining points in APr
    for (i in (length(APr) - 2):1) {
        nr = which(ChB[APr[i]:length(ChB)]) + APr[i] - 1
        lost = cumsum(pvarVec[nr])
        gain = abs(x[AP[i]] - x[nr + a])^p
        balance = gain - lost
        RSeg = nr[which.max(balance)]
        if (RSeg - APr[i] > 1) {
            ChB[(APr[i]):(RSeg - 1)] = FALSE
            pvarVec[RSeg] = Sum_pt(x, p, AP[i], a + RSeg)
        }
    }
    list(ChB = ChB, pvarVec = pvarVec)
}



### The functions of droping mieningless points
DropAllMid3 <-
function (x, p) 
{
    dn = 3
    if (length(x) > dn) {
        DropInt = which(MovingSum((abs(diff(x, lag = 1)))^p, 
            dn) < (abs(diff(x, lag = dn)))^p)
        while (length(DropInt) > 0) {
            Drop = as.vector(sapply(1:(dn - 1), FUN = function(ind) DropInt + 
                ind))
            x = x[-Drop]
            if (length(x) > dn) {
                DropInt = which(MovingSum((abs(diff(x, lag = 1)))^p, 
                  dn) < (abs(diff(x, lag = dn)))^p)
            }
            else {
                DropInt = vector("numeric", 0)
            }
        }
    }
    return(x)
}
DropAllMid3B <-
function (x, p) 
{
    dn = 3
    Len = length(x)
    DropB = rep(TRUE, Len)
    if (Len > dn) {
        ind = 1:(dn - 1)
        DropInt = which(MovingSum((abs(diff(x, lag = 1)))^p, 
            dn) < (abs(diff(x, lag = dn)))^p)
        while (length(DropInt) > 0) {
            Drop = rep(DropInt, each = length(ind)) + ind
            x = x[-Drop]
            DropB[DropB][Drop] = FALSE
            if (length(x) > dn) {
                DropInt = which(MovingSum((abs(diff(x, lag = 1)))^p, 
                  dn) < (abs(diff(x, lag = dn)))^p)
            }
            else {
                return(DropB)
            }
        }
    }
    return(DropB)
}
DropMid <-
function (x, p, dn) 
{
    dn = min(dn, length(x) - 1)
    DropInt = which(MovingSum((abs(diff(x, lag = 1)))^p, dn) < 
        (abs(diff(x, lag = dn)))^p)
    if (length(DropInt) > 0) {
        Drop = as.vector(sapply(1:(dn - 1), FUN = function(ind) DropInt + 
            ind))
        x = x[-Drop]
        new = TRUE
    }
    else new = FALSE
    return(list(x = x, new = new))
}
DropMidAll <-
function (x, p, dn) 
{
    DM = list(x = x, new = TRUE)
    if (dn >= 5) {
        while (DM$new) {
            x = DropMidAll(DM$x, p, dn - 2)
            DM = DropMid(x, p, dn = dn)
        }
    }
    else {
        if (dn > 2) {
            DM$x = DropAllMid3(x, p)
        }
    }
    return(DM$x)
}
DropMidAllB <-
function (x, p, dn) 
{
    DM = list(x = x, new = TRUE)
    DropB = rep(TRUE, length(x))
    if (dn >= 5) {
        while (DM$new) {
            DropB[DropB][!DropMidAllB(x[DropB], p, dn - 2)] = FALSE
            DM = DropMidB(x[DropB], p, dn = dn)
            DropB[DropB][!DM$DropB] = FALSE
        }
    }
    else {
        if (dn > 2) {
            DropB = DropAllMid3B(x, p)
        }
    }
    return(DropB)
}
DropMidB <-
function (x, p, dn) 
{
    dn = min(dn, length(x) - 1)
    ind = 1:(dn - 1)
    DropB = rep(TRUE, length(x))
    DropInt = which(MovingSum((abs(diff(x, lag = 1)))^p, dn) < 
        (abs(diff(x, lag = dn)))^p)
    if (length(DropInt) > 0) {
        Drop = rep(DropInt, each = length(ind)) + ind
        DropB[Drop] = FALSE
        new = TRUE
    }
    else new = FALSE
    return(list(DropB = DropB, new = new))
}

### finds analysing poits that maigth be meaningul
FindAnalysingPoints <-
function (x) 
{
    Len = length(x)
    if (!all(range(x)==sort(x[c(1,Len)]))) {
        stop("Error(003),FindAnalysingPoints: The vector is not Pseudo-Monotonic.")
    }
    if (x[1] < x[Len]) {
        minmaxf = cummax(x)
    }
    else {
        minmaxf = cummin(x)
    }
    CumPMPB = (x == minmaxf)
    CumPMPB[1] = FALSE
    return(which(CumPMPB))
}
### similar to FindAnalysingPoints, but from oposite direction, and includes ands.
FindAnalysingPointsOP <-
function (x) 
{
    Len = length(x)
	if (!all(range(x)==sort(x[c(1,Len)]))) {
        stop("Error(003),FindAnalysingPointsOP: The vector is not Pseudo-Monotonic.")
    }
    if (x[1] < x[Len]) {
        minmaxf = rev(cummin(rev(x)))
    }
    else {
        minmaxf = rev(cummax(rev(x)))
    }
    CumPMPB = (x == minmaxf)
    return(which(CumPMPB))
}

#######################################################################################
################# NEW functions: pvar optrations (under development) ##################
#######################################################################################

	Take1Variable <- function(Li, variable) Li[[variable]] 
	TakeVariables <- function(variable, L ){
		sapply(L, Take1Variable, variable=variable, simplify = FALSE)
	}
	

	
	### Vidine sunkcija sudedanti dvi p-variacijs (brud force), t.y. tiesmukai apjungiamos reiksmes ir info (jokio perskaiciavimo)
	### pertikrinama, ar sungti intervalai yra tame paciame taske		
	MergePvar <- function(..., JoinIfPossible = TRUE, dname = NULL){
		
		ArgL = list(...)
	
		if(length(ArgL)<1)
			stop("Function requires at least 1 argument.")
	
	
		### patikrinma ar abu reikiamos klases ir vienodu p
		AllClass = sapply(ArgL, class)
		if(!all(AllClass=="pvar"))
			stop("The `class` of arguments mus be `pvar`")
		
		AllP = unlist(TakeVariables("p", L=ArgL))
		if(!all(AllP==AllP[1]))
			stop("Adding two `pvars` is meaningfull only with same `p`")
		
		p = AllP[1]
		if(is.null(dname))
			dname = paste(TakeVariables("dname", L=ArgL), collapse=" + ")
		
		
		
		xL = TakeVariables("x", L=ArgL)		
		PartitionL = TakeVariables("Partition", L=ArgL)		
		
		#PartitionLNew = PartitionL
		if(class(xL)=="list"){
			LastPartitionPoint = tail(PartitionL[[1]],1)
			LastX = tail(xL[[1]],1)
			for(i in 2:length(xL)){
				if(JoinIfPossible & (head(xL[[i]],1)==LastX)){ # numetam x[i], jei x_i yra sujungimo taskas
					LastX = tail(xL[[i]],1)
					xL[[i]] = xL[[i]][-1]
					
					PartitionL[[i]] = PartitionL[[i]] + LastPartitionPoint - 1
					LastPartitionPoint = tail(PartitionL[[i]],1)
					PartitionL[[i]] = PartitionL[[i]][-1]
				}else{
					PartitionL[[i]] = PartitionL[[i]] + LastPartitionPoint
					LastPartitionPoint = tail(PartitionL[[i]],1)
				}
			}
			Partition = unlist(PartitionL)
			x = unlist(xL)
		}else{
			Partition = as.vector(PartitionL)
			x = as.vector(xL)
		}

		
				# xLNew = xL
		# #PartitionLNew = PartitionL
		# if(class(xL)=="list"){
			# for(i in 2:length(xL)){
				# LastPartition = tail(PartitionL[[i-1]],1)
				# if(
				# if(JoinIfPossible & (head(xL[[i]],1)==tail(xL[[i-1]],1))){ # numetam x[i], jei x_i yra sujungimo taskas
					# xLNew[[i]] = xL[[i]][-1]
					# PartitionL[[i]] = PartitionL[[i]][-1] + tail(PartitionL[[i-1]],1) - 1
				# }else{
					# PartitionL[[i]] = PartitionL[[i]] + tail(PartitionL[[i-1]],1) 
				# }
			# }
			# Partition = unlist(PartitionLNew)
			# x = unlist(xLNew)
		# }else{
			# Partition = as.vector(PartitionLNew)
			# x = as.vector(xLNew)
		# }

		
		
		
		
				
		 SplitPoints = Partition[SplitByExtremum(x[Partition])]		
		
		# MP = ChangePointsId(x[Partition])
		# SplitPoints = Partition[MP[SplitByExtremum(x[Partition][MP])]]
		
		# MP = ChangePointsId(x)
		# SplitPoints = MP[SplitByExtremum(x[MP])]
		
		
		TimeLabel = seq_along(x)
			
		
		pvar.value = unname(Sum_p(x[Partition],p))
		
		
		ans = list(
				value = c("p-variation" = pvar.value)
				,x = x
				,p = p
				,dname = dname
				,TimeLabel = TimeLabel
				,Partition = Partition
				,SplitPoints = SplitPoints
				,info = list(					
					PartitionN = length(Partition)
					,AnalysedPointsN = NA
					,SegmentsN = length(SplitPoints)-1
					,sizeN = NA
					,TakesTime = -1
				)
			)	
		class(ans) <- "pvar"
        return(ans)

	}
	

	
	### subset of p-var	
	PvarSubset <- function(pv, a=1, b=length(pv$x), dname=NULL){

		if(a>b)
			stop("The length of 'x' is zero")
		
		### patikrinma ar abu reikiamos klases ir vienodu p
		if(!all(class(pv)=="pvar"))
			stop("The `class` of arguments mus be `pvar`")
	
		p = pv$p
		x = pv$x[a:b]
		if(is.null(dname))
			dname = pv$dname %.% ", subset(a="%.%a%.%", b="%.%b%.%")"
		
		PartitionOK = pv$Partition[a <= pv$Partition &  pv$Partition <= b]

		if(length(PartitionOK)<2){
			if(length(PartitionOK)==1){
				ats = MergePvar(pvar(pv$x[a:PartitionOK], p), pvar(pv$x[PartitionOK:b], p) )
				ats$TimeLabel = pv$TimeLabel[a:b]
				return(ats)
			}else{
				ats = pvar(pv$x[a:b], p, TimeLabel=pv$TimeLabel[a:b])
				return(ats)
			}		
		}
		
		pv1 = pvar(pv$x[a:head(PartitionOK,1)], p)
		pv2 = pvar(pv$x[tail(PartitionOK,1):b], p)


		PartitionL = list()
		PartitionL[[1]] = pv1$Partition 
		PartitionL[[2]] = PartitionOK[-1] - PartitionOK[1] + tail(PartitionL[[1]],1)
		PartitionL[[3]] = pv2$Partition[-1] + tail(PartitionL[[2]],1) - 1
		Partition = unlist(PartitionL)
		
		
		SplitPoints = Partition[SplitByExtremum(x[Partition])]
		
		# MP = ChangePointsId(x[Partition])
		# SplitPoints = Partition[MP[SplitByExtremum(x[Partition][MP])]]
		
		# MP = ChangePointsId(x)
		# SplitPoints = MP[SplitByExtremum(x[MP])]		
		
		pvar.value = unname(Sum_p(x[Partition],p))
		
		ans = list(
				value = c("p-variation" = pvar.value)
				,x = x
				,p = p
				,dname = dname
				,TimeLabel = pv$TimeLabel[a:b]
				,Partition = Partition
				,SplitPoints = SplitPoints
				,info = list(					
					PartitionN = length(Partition)
					,AnalysedPointsN = NA
					,SegmentsN = length(SplitPoints)-1
					,sizeN = NA
					,TakesTime = -1
				)
			)

		class(ans) <- "pvar"
        return(ans)			
	}

	
	
	Add2Pvar <- function(pv1, pv2, JoinIfPossible = TRUE, dname = NULL){
	
	
		### patikrinma ar abu reikiamos klases ir vienodu p
		if(!(class(pv1)==class(pv2)& class(pv2)=="pvar"))
			stop("The `class` of arguments mus be `pvar`")
		
		if(pv1$p!=pv2$p)
			stop("Adding two `pvars` is meaningfull only with same `p`")
		p = pv1$p
	
		pv1range = range(pv1$x)
		pv2range = range(pv2$x)

		XSplitPoints1 = pv1$x[pv1$SplitPoints]
		XSplitPoints2 = pv2$x[pv2$SplitPoints]
		
		### mums reikia pirmo (paskutinio) splitPointo, kuris uzeina uz kito - p-variacijos rezio
		#	tai yra pirmas splitPoint taskas kuris yra uz jungiamosios dalies range srities. Ir toks grynas (t.y. be trukiu).
		# Jie tokio nera, tai tiks ir siap pirmas taskas
		
		JoinP1 = tail(c(pv1$SplitPoints[1], pv1$SplitPoints[XSplitPoints1<pv2range[1] | XSplitPoints1>pv2range[2]]),1)
		JoinP2 = head(c(pv2$SplitPoints[XSplitPoints2<pv1range[1] | XSplitPoints2>pv1range[2]], tail(pv2$SplitPoints,1)),1) 
	
	
	
		### patikrinam ar sujungiam x'us tiesmukai ir paruosiame vidurine dali
		### efekyvumo delei galetume tik partition points nagrineti
		if(JoinIfPossible & (tail(pv1$x,1)==head(pv2$x,1))){
			xMid = c(pv1$x[JoinP1:length(pv1$x)], pv2$x[2:JoinP2])
		}else{
			xMid = c(pv1$x[JoinP1:length(pv1$x)], pv2$x[1:JoinP2])
		}
		
		

		pvMid = pvar(xMid, p)
	
		ans = MergePvar(PvarSubset(pv1,b=JoinP1),pvMid, PvarSubset(pv2,a=JoinP2))
	
		if(is.null(dname))
			ans$dname <- deparse(substitute(pv1)) %.%" + "%.% deparse(substitute(pv2))
		
		ans	
	}
	
	### yra rimtu bugu
	AddPvar <- function(..., JoinIfPossible = TRUE, dname = "data"){
			
		ArgL = list(...)
	
		if(length(ArgL)<1)
			stop("Function requires at least 1 argument.")
	
	
		### patikrinma ar abu reikiamos klases ir vienodu p
		AllClass = sapply(ArgL, class)
		if(!all(AllClass=="pvar"))
			stop("The `class` of arguments mus be `pvar`")
		
		AllP = unlist(TakeVariables("p", L=ArgL))
		if(!all(AllP==AllP[1]))
			stop("Adding two `pvars` is meaningfull only with same `p`")
			
		ans = ArgL[[1]]	
		if(length(ArgL)>1){
			for(i in 2:length(ArgL)){
				ans = Add2Pvar(ans,ArgL[[i]], JoinIfPossible=JoinIfPossible)
			}
		}		
		
		ans$dname = dname
		
		return(ans)				
	
	}
	
	
	
	
	
	
