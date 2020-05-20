tenonradial <- function(formula, data, subset,
                        rts = c("C", "NI", "V"), 
                        base = c("output", "input"), 
                        ref = NULL, data.ref = NULL, subset.ref = NULL,
                        full.solution = TRUE,
                        print.level = 1)
{
  lmdConstr <- TRUE
 if( !is.null(ref) & is.null(data.ref) ){
  warning("If you use variable names in 'ref', 'data.ref' is required", call. = FALSE)
 }
 
 if( is.null(ref) & !is.null(data.ref) ){
  stop("If you use 'data.ref', 'ref' is required", call. = FALSE)
 }
 
 winw <- getOption("width")
 if (winw > 100+5){
  winw <- 100
 }
 else if (winw > 90+5) {
  winw <- 90
 }
 else if (winw > 80+5) {
  winw <- 80
 }
 else if (winw > 70+5) {
  winw <- 70
 }
 else if (winw > 60+5) {
  winw <- 60
 }
 else if (winw > 50+5) {
  winw <- 50
 }
 else {
  winw <- 0
 }
 
 t1 <- .prepareYX(formula = formula, data = data, subset = subset, rts = rts,
                  base = base, ref = ref,	data.ref = data.ref, subset.ref = subset.ref,
                  print.level = print.level, type = "RM", winw = winw, sysnframe = sys.nframe())

 te <- .teNonrad2(t(t1$y), t(t1$x), ncol(t1$y), ncol(t1$x), nrow(t1$y),
              t(t1$y.ref), t(t1$x.ref), nrow(t1$y.ref), t1$myrts, t1$mybase,
              ifqh = FALSE, print.level = print.level,
              lmdConstr = lmdConstr, full.solution = full.solution)
 # cat("\n Printing te \n")
 # print(te)
 # cat.print(cbind(te$te, colMeans(te$te.all)))
 # return(te)
 # te <- ifelse(abs(te - 1) < 1e-12, 1, te)
 # te <- ifelse(te == -998, NA, te)
 if(print.level >= 3){
  cat("\n")
 }
 if(print.level >= 2){
  cat(paste("",rep("_", (winw-10)/1),"", sep = ""), "\n\n", sep = "")
  cat("Summary of efficiencies:\n\n", sep = "")
  .su(te, print = FALSE)
 }
 if(full.solution){
   tymch <- list(call = match.call(), model = "tenonradial", K = nrow(t1$y), M = ncol(t1$y), N = ncol(t1$x), Kref = nrow(t1$y.ref), rts = t1$rts.string, base = t1$base.string, te = te$te, te.detail = te$te.detail, intensity = te$intensity, esample = t1$esample, esample.ref = t1$esample.ref)
 }
 else {
   tymch <- list(call = match.call(), model = "tenonradial", K = nrow(t1$y), M = ncol(t1$y), N = ncol(t1$x), Kref = nrow(t1$y.ref), rts = t1$rts.string, base = t1$base.string, te = te, esample = t1$esample, esample.ref = t1$esample.ref)
 }
 
 class(tymch) <- "npsf"
 return(tymch)
}



#