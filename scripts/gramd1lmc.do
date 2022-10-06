
* Melanoma specific survival (MSS) analysis of GRAMD1B in Leeds Melanoma Cohort (LMC) - STATA v16
* Figure 4C, Suresh et al 2022, Cancer Discovery
* Data accession number: EGAS00001002922
* Jeremie Nsengimana, Newcastle University, October 2022
* jeremie.nsengimana@newcastle.ac.uk

clear 
clear mata
clear matrix
set maxvar 30000
use "daslgenelist.dta", clear
list if gene=="GRAMD1B" /* to find out the probe name present in data: it is ilmn_3237376*/

use "LMCdaslgx.dta",clear
hist ilmn_3237376
hist breslow
gen logbres=log(breslow)
hist logbres

stset ostime if mss=="yes", fail(died) /* ignoring deaths from non-melanoma causes */
stcox ilmn_3237376   
xtile gramd1b=ilmn_3237376
stcox gramd1b 
xi:stcox i.gramd1b*logbres

gen bres2split3mm=1 if breslow<=3
replace bres2split3mm=2 if breslow>3&breslow<.

xi:stcox i.gramd1b*i.bres2split3mm 
sts gr, by(gramd1b bres2split3)
sts gr, by(gramd1b bres2split3) tmax(13) title("MSS by GRAMD1B and tumour thickness",color(b)) ///
 graphregion(lcolor(white)) graphregion(fcolor(white)) ylabel(,nogrid) xtitle(" Follow up time(years)") ///
 text(0.1 3 "Interaction P=0.005" 0.3 9 "low GRAMD1B/Breslow>3mm"0.5 10 "high GRAMD1B/breslow>3mm", size(small)) ///
 text(0.68 9.5 "high GRAMD1B/Breslow<=3mm"0.85 11 "low GRAMD1B/breslow<3mm", size(small)) ///
 ytitle("Survival fraction") legend(off)
 gr export KMplotgramd1b_breslow.pdf 

 
