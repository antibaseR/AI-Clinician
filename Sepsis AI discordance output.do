clear 

local path "/Users/chenxinlei/Library/Mobile Documents/com~apple~CloudDocs/Projects/AI Clinician"

log using "`path'/Sepsis AI discordance output Jul-3-2024.log", replace
use "`path'/Data/Sepsis AI - ICU Cohort - Discordance Scores.dta", clear

merge 1:1 interval empi_id hosp_id using "`path'/Data/Sepsis AI - ICU Cohort - Treated Variable.dta"


*--- 1. Data management and preliminary analysis
/* check number of obervations per patient */
codebook empi_id
sort empi_id interval
bys empi_id: gen run=_n
bys empi_id: egen maxrun=max(run)


tab maxrun if run==1
summ maxrun, d
hist maxrun if run==1

/* generate age groups */
gen agegp=.
replace agegp=1 if age<40
replace agegp=2 if age>=40 & age<50
replace agegp=3 if age>=50 & age<60
replace agegp=4 if age>=60 & age<70
replace agegp=5 if age>=70 & age<80
replace agegp=6 if age>=80 & age<90
replace agegp=7 if age>=90 & age<.

label define agelb 1 "<40" 2 "40s" 3 "50s" 4 "60s" 5 "70s" 6 "80s" 7 "90+"
label values agegp agelb

/* generate Elixhauser groups */
gen elixgp=.
replace elixgp=1 if elix==0
replace elixgp=2 if elix>=1 & elix<=2
replace elixgp=3 if elix>=3 & elix<=4
replace elixgp=4 if elix>=5 & elix<.

label define elixlb 1 "0" 2 "1-2" 3 "3-4" 4 "5+"
label values elixgp elixlb

gen sqrtsofa=sqrt(sofa_total)

/* generate sofa groups by tertiles */
xtile sofagp = sofa_total, nq(3)
tabulate sofagp

/* check whether y=dead_90 varying within patient
   -- Among N=25,175 patients,
      n=16,713 (66.39%) patients have all dead_90 values =0
      n= 7,470 (29.67%) patients have all dead_90 values =1 */
bys empi_id: egen mdead_90=mean(dead_90)
tab mdead_90 if run==1

/* assign integer numerical scores to the state-level discordance scores */
gen score=.
replace score=1 if state_discordance==0
replace score=2 if state_discordance==1
replace score=3 if state_discordance>1 & state_discordance <1.5
replace score=4 if state_discordance==2
replace score=5 if state_discordance>2 & state_discordance <2.5
replace score=6 if state_discordance>2.5 & state_discordance <3
replace score=7 if state_discordance==3
replace score=8 if state_discordance>3 & state_discordance <3.5
replace score=9 if state_discordance>3.5 & state_discordance <4
replace score=10 if state_discordance==4
replace score=11 if state_discordance>4 & state_discordance <4.2
replace score=12 if state_discordance>4.2 & state_discordance <4.4
replace score=13 if state_discordance>4.4 & state_discordance <5
replace score=14 if state_discordance==5
replace score=15 if state_discordance>5 & state_discordance <6

label var score "ordinal discordance score"
label define scorelb 1 "0" 2 "1" 3 "1.414214" 4 "2" 5 "2.236068" 6 "2.828427" 7 "3" 8 "3.162278" 9 "3.605551" 10 "4" 11 "4.123106" 12 "4.24264" 13 "4.472136" 14 "5" 15 "5.656854"
label values score scorelb 

/* create binary diconcordance scores */
forvalues i=1/14 {
	gen score_bin`i'=0
	replace score_bin`i'=1 if score >`i' & score !=.
	replace score_bin`i'=. if score ==.
}


/* check whether score varying within patient
   Among N=25,175 patients, n=512 (2.03%) had no varying state-level discordance scores */
bys empi_id: egen minscore=min(score)
bys empi_id: egen maxscore=max(score)
gen diffscore = maxscore - minscore

tab diffscore if run==1


*--- 2. Analyis

/*
* Models stratified by SOFA scores
forvalues j=1/3 {
	forvalues i=1/14 {
		melogit score_bin`i' dead_90 if sofagp==`j' || empi_id:, nolog
		estimates store crude_primary_`i'
		estat icc
		margins, at(dead_90=(0 1))
	}
	forvalues i=1/14 {
		melogit score_bin`i' dead_90 interval lactate age gender elix if sofagp==`j' || empi_id:, nolog
		estimates store all_primary_`i'
		lrtest crude_primary_`i' all_primary_`i'
		estat icc
		margins, at(dead_90=(0 1)) atmeans 
	}
}


* Models among vasopressor-only
forvalues i=1/12 {
		melogit score_bin`i' dead_90 if norepi_equiv!=0 & sofagp==1 || empi_id:, nolog
		*estimates store crude_primary_`i'
		*estat icc
		margins, at(dead_90=(0 1))
		}

forvalues i=1/12 {
		melogit score_bin`i' dead_90 interval lactate age gender elix if norepi_equiv!=0 & sofagp==1 || empi_id:, nolog
		*estimates store all_primary_`i'
		*lrtest crude_primary_`i' all_primary_`i'
		estat icc
		margins, at(dead_90=(0 1)) atmeans 
	}




forvalues j=2/3 {
	
	forvalues i=1/14 {
		melogit score_bin`i' dead_90 if norepi_equiv!=0 & sofagp==`j' || empi_id:, nolog
		*estimates store crude_primary_`i'
		*estat icc
		margins, at(dead_90=(0 1))
		}
		
	
	forvalues i=1/14 {
		melogit score_bin`i' dead_90 interval lactate age gender elix if norepi_equiv!=0 & sofagp==`j' || empi_id:, nolog
		*estimates store all_primary_`i'
		*lrtest crude_primary_`i' all_primary_`i'
		estat icc
		margins, at(dead_90=(0 1)) atmeans 
	}

}


* Models among intravenous fluid only
forvalues j=1/3{
	
	forvalues i=1/14 {
		melogit score_bin`i' dead_90 if fluid!=0 & sofagp==`j' || empi_id:, nolog
		*estimates store crude_primary_`i'
		*estat icc
		margins, at(dead_90=(0 1))
	}
	
	forvalues i=1/14 {
		melogit score_bin`i' dead_90 interval lactate age gender elix if fluid!=0 & sofagp==`j' || empi_id:, nolog
		estimates store all_primary_`i'
		*lrtest crude_primary_`i' all_primary_`i'
		*estat icc
		margins, at(dead_90=(0 1)) atmeans 
	}
	
}
*/


* Models only among treated
keep if treated==1
forvalues j=1/3 {
	forvalues i=1/14 {
		melogit score_bin`i' dead_90 if sofagp==`j' || empi_id:, nolog
		estimates store crude_primary_`i'
		estat icc
		margins, at(dead_90=(0 1))
	}
	forvalues i=1/14 {
		melogit score_bin`i' dead_90 interval lactate age gender elix if sofagp==`j' || empi_id:, nolog
		estimates store all_primary_`i'
		lrtest crude_primary_`i' all_primary_`i'
		estat icc
		margins, at(dead_90=(0 1)) atmeans 
	}
}




log close

