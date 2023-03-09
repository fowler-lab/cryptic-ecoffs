* ECOFFS.DO

* Estimate Ecoffs for CRyPTIC dataset

cap log close
log using ecoffs.log, replace
noi di "Run by ASW on $S_DATE $S_TIME"
noi version

clear
insheet using CRyPTIC-MIC-distributions.csv
rename v1 drug 
rename v2 type 
rename v3 plate
rename v4 well 
rename v5 micstr 
rename v6 n
gen byte signmic=0
replace signmic=1 if strpos(micstr,">")>0
* csv has replaced <= with ? 
replace micstr=subinstr(micstr,"?","<=",.)
replace signmic=-1 if strpos(micstr,"<=")>0
label define signlbl -1 "<=" 0 "=" 1 ">"
label values signmic signlbl
gen str20 temp=micstr
replace temp=subinstr(temp,">","",.)
replace temp=subinstr(temp,"<=","",.)
gen mic=real(temp)
assert mic<.
drop temp
* generate interval censored log2 mic treating actual mic as observed
gen left = log10(mic)/log10(2) if signmic>=0
gen right = log10(mic)/log10(2) if signmic<=0
* PROBLEM: this goes down to -infinity and up to +infinity for left and right censored observations, which then won't converge
* constrain at +/- 5 DD below/above minimum respectively - across plates
for any min max: egen X=X(mic),by(drug)
replace left=log10(min)/log10(2)-5 if left>=.
replace right=log10(max)/log10(2)+5 if right>=.
* generate interval censored log2 mic treating observed mic as interval censored too
* NB not exactly 1DD between - argh - do by plate so get correct
* doesn't work anyway see below but keep to remember not to try again! 
sort drug type plate mic signmic
by drug type plate mic signmic: assert _N==1
gen LEFT=left
gen RIGHT=right
by drug type plate: replace LEFT=RIGHT[_n-1]+0.01 if signmic==0 & _n>1
compress
save rawdata, replace

* set up a file to store all the means and SDs from 2 component mixture models
cap postclose ecoffs
postfile ecoffs str3 drug str8 type str8 method comp1_b1 comp1_v1 str3 fail comp2_cons comp2_b1 comp2_v1 comp2_b2 comp2_v2 using ecoffs_output, replace

* loop over drugs, over types (raw vs SSSS only)
foreach drug in RIF INH EMB RFB BDQ CFZ DLM ETH LZD LEV MXF AMI KAN PAS {
* starred out line is test for one drug - swap star to just run one at a time
*foreach drug in RIF {
	noi di _n(2) _dup(80) "*" _n _dup(80) "*" _n "`drug' - ALL raw (UKMYC5 & UKMYC6)" _n _dup(80) "*" _n _dup(80) "*" 
        noi di _n "Check MIC thresholds match what expected from plate"
        noi tab micstr plate [fweight=n] if drug=="`drug'" & type=="raw"
        noi di _n "Observations"
        noi tab mic signmic [fweight=n] if drug=="`drug'" & type=="raw"

	noi di _n _dup(80) "=" _n "`drug' - ALL raw (UKMYC5 & UKMYC6) - one component model" _n _dup(80) "=" 
        qui intreg left right [fweight=n] if drug=="`drug'" & type=="raw"
        noi intreg
        est store component1
        local comp1_b1=_b[_cons]
        local comp1_v1=exp(_b[lnsigma:_cons])

        * two component model: this won't always work - but even when convergence is not achieved 
        * it seems to give a return code of 0 - which it shouldn't! check for failure of convergence via missing SE instead
        cap fmm 2, difficult iterate(10000): intreg left right [fweight=n] if drug=="`drug'" & type=="raw"
        est store component2
	* get estimates to store in spreadsheet
        * wildtype will be class 1: need to get mean and SD on log2 scale and convert to percentile
        * mean is simple _b[1.Class], variance is a bugger - have to get from matrix
        mat b=e(b)
        local comp2_cons=b[1,2]
        local comp2_b1=b[1,3]
        local comp2_b2=b[1,4]
        local comp2_v1=b[1,5]
        local comp2_v2=b[1,6]
        * non-covergence indicated by missing VARIANCE of the mean OR variance estimates of the distribution parameters
        mat v=e(V)
	mat v=vecdiag(v)
        local V1=v[1,3]
        local V2=v[1,4]
        local V3=v[1,5]
        local V4=v[1,6]
        if `V1'<=0|`V2'<=0|`V3'<=0|`V4'<=0 {
	   local fail "yes"
	   noi di _n _dup(80) "=" _n "`drug' - ALL raw (UKMYC5 & UKMYC6) - TWO component model - DID NOT CONVERGE - DON'T USE RESULTS" _n _dup(80) "=" 
	}
        else {
	   local fail "no" 
	   noi di _n _dup(80) "=" _n "`drug' - ALL raw (UKMYC5 & UKMYC6) - TWO component model CONVERGED - use lower AIC/BIC" _n _dup(80) "=" 
	}
        noi est stats component*
        noi est replay component2
        if "`fail'"=="no" {
           noi di _n "WT comprises " %4.1f 100-100*exp(`comp2_cons')/(exp(`comp2_cons')+1) "% of the population"
           noi di    "Mean (var) log2 MIC in WT is " `comp2_b1' " (" `comp2_v1' ")"
           noi di    "Mean (95th PCT) MIC in WT is " 2^`comp2_b1' " (" 2^(`comp2_b1'+1.96*sqrt(`comp2_v1')) ")"
        }
        else {
           noi di _n "Mean (var) log2 MIC in single population is " `comp1_b1' " (" `comp1_v1' ")"
           noi di    "Mean (95th PCT) MIC in single population is " 2^`comp1_b1' " (" 2^(`comp1_b1'+1.96*sqrt(`comp1_v1')) ")"
        }
        post ecoffs ("`drug'") ("raw") ("actual") (`comp1_b1') (`comp1_v1') ("`fail'") (`comp2_cons') (`comp2_b1') (`comp2_v1') (`comp2_b2') (`comp2_v2')

/* won't converge using intervals for everything - not enough data - keep as placeholder to remember not to try this again!
	noi di _n _dup(80) "=" _n "`drug' - ALL raw (UKMYC5 & UKMYC6) BASED ON ALL MICS TREATED AS INTERVAL CENSORED" _n _dup(80) "=" 
        qui fmm 2: intreg LEFT RIGHT [fweight=n] if drug=="`drug'" & type=="raw"
*/

	noi di _n _dup(80) "=" _n "`drug' - SSSS and SSUS (UKMYC5 & UKMYC6) BASED ON ACTUAL MIC - ONE COMPONENT ONLY" _n _dup(80) "=" 
        noi di _n "Observations"
        noi tab mic signmic [fweight=n] if drug=="`drug'" & inlist(type,"SSSS","SSUS")
        * one component model ONLY relevant here
        qui intreg left right [fweight=n] if drug=="`drug'" & inlist(type,"SSSS","SSUS")
        noi intreg
        local b1=_b[_cons]
        local v1=exp(_b[lnsigma:_cons])
        noi di _n "Mean (var) log2 MIC in SS*S is " `b1' " (" `v1' ")"
        noi di    "Mean (95th PCT) MIC in SS*S is " 2^`b1' " (" 2^(`b1'+1.96*sqrt(`v1')) ")"
        post ecoffs ("`drug'") ("SS*S") ("actual") (`b1') (`v1') ("") (.) (.) (.) (.) (.)
}

postclose ecoffs

use ecoffs_output, clear
* calculate on MIC scale
gen comp1_mic_mean=2^comp1_b1
gen comp1_mic_p95=2^(comp1_b1+1.96*sqrt(comp1_v1))
gen comp2_mic1_mean=2^comp2_b1
gen comp2_mic1_p95=2^(comp2_b1+1.96*sqrt(comp2_v1))
gen comp2_mic2_mean=2^comp2_b2
gen comp2_mic2_p95=2^(comp2_b2+1.96*sqrt(comp2_v2))
gen comp2_perc2=100*exp(comp2_cons)/(exp(comp2_cons)+1) 
gen comp2_perc1=100-comp2_perc2
for any mean p95: replace comp2_mic1_X=comp1_mic_X if type=="SS*S"
format comp* %8.3f
order drug type comp1_mic* comp1_b comp1_v fail comp2_perc1 comp2_mic1* comp2_perc2 comp2_mic2* comp2_b* comp2_v*
drop method
compress
export excel using ecoffs_output.xlsx, replace firstrow(var)
save ecoffs_output, replace

