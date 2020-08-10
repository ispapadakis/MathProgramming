#' ---
#' title: "MINIMIZE TRANSPORTATION COST USING R"
#' author: "Yanni Papadakis"
#' date: "August 2020"
#' ---

#' <ul> <li>Solves a Linear Program Using lpSolveAPI</li>
#' <li>Solves a Quadratic Program Using quadprog</li>
#' <li>Solves a Non-Convex (Global Optimization) Problem Using DEoptim</li> </ul>

#' Note that low-cost for-profit solvers are much faster and necessary for large problems.<br>
#' Inspired by Example 2.1.8 in C. Floudas Textbook

library(lpSolveAPI, quietly = TRUE)
library(knitr, quietly = TRUE)

#' # CASE A: Linear Cost Structure

#' ### Load Problem Data
df = read.csv("transprob.csv", row.names=1)

#' Label Demand and Supply Points
demand.points = paste("D",1:4,sep="")
supply.points = paste("S",1:6,sep="")

#' Ensure total supply equals total demand. <br>
#' Not necessary for linear formulation, but helpful later.
stopifnot(sum(df[supply.points,'LS']) == sum(df['LD',demand.points]))

kable(df,caption="Cost Structure for Linear Problem and Supply/Demand Limits")

#' Ingest raw data to R objects
cost.data = data.matrix(df[supply.points,demand.points])
n = length(supply.points)
m = length(demand.points)

#' ### Build Linear Program
lps.model = make.lp(0,n*m)
name.lp(lps.model, "Transportation Problem")

r = 0
#' Supply Constraints
tmpl = seq(0,n*m-1,n)
for(i in 1:n){
  r = r + 1
  add.constraint(lps.model, rep(1,m), "<=", df[i,"LS"], tmpl+i)
}
#' Demand Constraints
tmpl = seq(1,n)
for(j in 1:m){
  r = r + 1
  add.constraint(lps.model, rep(1,n), ">=", df["LD",j], tmpl+n*(j-1))
}
dimnames(lps.model)[[1]] = c(supply.points, demand.points)
dimnames(lps.model)[[2]] = c(outer(supply.points,demand.points,paste,sep='_'))

#' Objective
set.objfn(lps.model, c(cost.data))

#' Call solve function from lpSolveAPI
print(solve(lps.model))

#' If status code is 0, optimal solution was found.<br>

#' ### Optimal Solution

#' Optimal Objective Value
print(get.objective(lps.model))
#' Optimal Flow from S to D Points
v = get.variables(lps.model)
constr = get.constraints(lps.model)

soln = rbind(
  cbind(
    matrix(round(v,3),n,m, dimnames=list(supply.points,demand.points)),
    TS=constr[1:n],
    LS=df[1:n,"LS"]
  ),
  TD=c(constr[n + (1:m)],0,0),
  LD=c(df["LD",1:m],0,0)
)
kable(soln,caption="Flow from S point to D Point, T: Total, L: Limit")

#' Save LP
write.lp(lps.model,"transprob.lp")

#' Get Dual Solution
x  = get.dual.solution(lps.model)
#' Dual Prices for Supply
print(data.frame(dp=x[1+(1:n)],row.names=supply.points))
#' Dual Prices for Demand
print(data.frame(dp=x[1+n+(1:m)],row.names=demand.points))

#' Cost Differential to Change Optimal Flow
print(matrix(x[-(1:(1+n+m))],n,m,dimnames=list(supply.points,demand.points)))


#' # CASE B: Quadratic Costs - Negative Externality (e.g. congestion costs)

#' Introducing a Positive Definite Term in the Cost Function

Q = data.frame(
  D1 = c(7,12,13,7,4,17),
  D2 = c(4,9,12,9,10,9),
  D3 = c(6,14,8,16,21,8),
  D4 = c(8,7,4,8,13,4),
  row.names = supply.points
)
kable(Q)

#' ### Ingest data in objects agreeing with quadprog library

#' D-Matrix
Dmat = 2 * diag(c(data.matrix(Q)))
#' d-vector
dvec = c(cost.data)
#' #### A-Matrix
Amat = matrix(0,n+m,n*m)
r = 0
#' Supply Constraints
tmpl = seq(0,n*m-1,n)
for(i in 1:n){
  r = r + 1
  Amat[r,tmpl+i] = 1
}
#' Demand Constraints
tmpl = seq(1,n)
for(j in 1:m){
  r = r + 1
  Amat[r,tmpl+n*(j-1)] = 1
}
Amat = rbind(Amat,diag(n*m))

#' b-vector
bvec = c(unlist(df[1:n,"LS"]), unlist(df["LD",1:m]), rep(0, n*m))

#' ### Minimize Cost
qp = quadprog:::solve.QP(Dmat,-dvec,t(Amat),bvec)

kable(matrix(round(qp$solution,3),n,m, dimnames=list(supply.points,demand.points)),
      caption="Solution for Quadratic Cost Structure (Case B)")

#' Compare to Solution Arising When Cost Structure is Linear (Case A)
kable(matrix(v,n,m, dimnames=list(supply.points,demand.points)),
      caption="Linear Cost Solution (Case A)")

#' Now, more S/D combinations are utilized and optimal flow is not integral.

#' # CASE C: Quadratic Costs - Cost Discount (e.g. economies of scale)

#' Just as in Floudas's Example<br>
#' Quadratic Term is not Positive Definite. But, problem is separable.<br>

#' <p style="font-size:20px">Quick Solution Strategy: Use Global Optimization.</p>

#' This is not computationally efficient.
#' There are a better ways to do this, 
#' but using a generic solver is tempting when problem size is small.
#' 

#' ### Set Up Functions to help the DEoptim formulation


#' Couldn't find an element-wise min function, so here is a home-made one
elementwise.min = function(x,y) ifelse(x>y,y,x)

#' Assign upper flow limits depending on the corresponding demand and supply constraints
ul = outer(df[supply.points,'LS'],
           as.numeric(df['LD',demand.points]),
           elementwise.min)
#' Make sure flow is not large enough to make any S/D transporation cost negative
ul = mapply(elementwise.min, as.data.frame(cost.data / Q), as.data.frame(ul))


#' Cost Evaluation Function
cost.eval = function(x){
  sum(x * (cost.data - x*Q) )
}

#' ### Decision Variables

#' Initial Solution: Chosen arbitrarily and leaving room for improvement
x.init = matrix(0,n,m)
x.init[1,1] = 2
x.init[1,2] = 6
x.init[2,1] = 3
x.init[2,4] = 21
x.init[3,1] = 19
x.init[3,3] = 1
x.init[4,1] = 2
x.init[4,2] = 22
x.init[5,1] = 3
x.init[5,2] = 1
x.init[5,3] = 12
x.init[6,2] = 12


#' Consider a perturbation flow,
#' which if added to the initial solution flow can give any desired feasible flow.<br>
#' The dimension of the perturbation is (n-1) X (m-1)<br>

#' Here is a function to give a new solution, given a perturbation and the initial solution
perb2full = function(u,x0) {
  try(if(!is.matrix(x0)) stop("pmat is of class matrix"))
  d = dim(x0)
  out = matrix(0,d[1],d[2])
  out[-d[1],-d[2]] = matrix(u,d[1]-1,d[2]-1)
  out[-d[1], d[2]]= -apply(out[-d[1],-d[2]],1,sum)
  out[d[1],] = -apply(out[-d[1],],2,sum)
  x0+out
}


#' Problem Cost (adds penalty for negative flows)
probcost = function(u, penalty=cost.data+1000) {
  x = perb2full(u,x.init)
  cost.eval(x) - sum( penalty[x<0] * x[x<0] )
}

#' ### Minimize Cost

#' Cost of Initial Solution
probcost(rep(0,(n-1)*(m-1)))

#' Compute Optimum (Be Patient!)
set.seed(2020)
print(Sys.time())
start.time = Sys.time()
soln = DEoptim:::DEoptim(
  probcost,
  -c(x.init[-n,-m]),
  c(ul[-n,-m]-x.init[-n,-m]),
  control=DEoptim:::DEoptim.control(NP=150,itermax=9000,trace=FALSE,strategy=5,c=0.4),
  fnMap = round # requires integrality
)
print(Sys.time()-start.time)

#' ### Min Cost Flow Obtained Using the Global Optimization Heuristic

#' Near-Optimal Cost
print(soln$optim$bestval)

#' Near-Optimal Flow
round(perb2full(soln$optim$bestmem,x.init),2)


