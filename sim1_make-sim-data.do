/// Simulation 1 - Data-generating script

// Set working directory
cd "Z:\jointmodels\projects\Informative observation process"
clear

// Simulate data
local B = 1000
local N = 200
local expN = 2500
local a0 = 0
local a1 = 1
local a2 = 0.2
local sigmau2 = 1
local sigmav2 = 0.5
local sigmage2 = 0.1
local sigmae2 = 1
local p = 1.05
local beta = 1
local dgm = 1
set seed 286826724

// Scenarios generated from a joint model
clear
foreach gamma in 0 1.50 {
	foreach lambda in 0.10 0.30 1.00 {
		nois _dots 0, title(Lambda: `lambda'; Gamma: `gamma') reps(`B')
		forvalues i = 1/`B' {
			quietly {
				clear
				set obs `N'
				generate id = _n
				generate trt = runiform() > 0.5
				generate u = rnormal(0, `=sqrt(`sigmau2')')
				generate v = rnormal(0, `=sqrt(`sigmav2')')
				// admin censoring time
				generate adcens = runiform(5, 10)
				// recurrent event/observation times
				expand `expN'
				sort id
				generate t0 = 1
				bysort id: replace t0 = 0 if _n > 1
				survsim obtime, lambdas(`lambda') gammas(`p') cov(trt `beta' u 1) distribution(weibull)
				generate t = 0
				bysort id: generate t2 = sum(obtime)
				bysort id: replace t = t2[_n - 1] if _n > 1
				drop if t > adcens
				bysort id: replace t2 = adcens if _n == _N
				generate ind = 1
				bysort id: replace ind = 0 if t2 == adcens
				drop t0 t2
				// generate longitudinal outcome
				generate xb = `a0' + trt * `a1' + t * `a2' + `gamma' * u + v
				generate y = rnormal(xb, `=sqrt(`sigmae2')')
				// generate number of measurements per individual
				bysort id: generate nn = _N
				// generate cumulative number of measurements per individual
				bysort id: generate cumnn = _n
				// generate centered number of measurements per individual
				summarize nn if cumnn == 1
				generate cnn = nn - r(mean)
				// make Andersen-Gill counting process style data
				sort id t
				bysort id: generate t1 = t[_n+1]
				sort id t
				replace t1 = adcens if t1 == .
				// estimate normalised inverse intensity of visit weights
				stset t1, fail(ind) exit(time .) id(id) enter(t)
				stcox trt, nohr robust
				stset, clear
				predict lp, xb
				generate w = 1 / exp(lp)
				summarize w
				generate wn = w - r(mean) + 1
				bysort id: generate iivw = wn[_n - 1]
				replace iivw = 1 if iivw == .
				drop lp w wn t1
				// add gamma, lambda, DGM
				generate dgm = `dgm'
				generate i = `i'
				// compress
				compress
				// save dataset
				save "Data/df_`dgm'_`i'", replace
			}
			nois _dots `i' 0
		}
		local dgm = `dgm' + 1
	}
}

// Additional scenarios where the observation process is generated from a Gamma distribution
clear
foreach shape in 2.00 {
	foreach theta in 0.00 2.00 {
		nois _dots 0, title(Shape: `shape'; Theta: `theta') reps(`B')
		forvalues i = 1/`B' {
			quietly {
				clear
				set obs `N'
				generate id = _n
				generate trt = runiform() > 0.5
				generate v = rnormal(0, `=sqrt(`sigmav2')')
				generate ge = rnormal(0, `=sqrt(`sigmage2')')
				// admin censoring time
				generate adcens = runiform(5, 10)
				// recurrent event/observation times from a Gamma distribution depending on trt if theta != 0
				expand `expN'
				sort id
				generate t0 = 1
				bysort id: replace t0 = 0 if _n > 1
				local gamma = 0
				generate obtime = rgamma(`shape', exp(`theta' * (-1 * `beta') * trt + ge))
				generate t = 0
				bysort id: generate t2 = sum(obtime)
				bysort id: replace t = t2[_n - 1] if _n > 1
				drop if t > adcens
				bysort id: replace t2 = adcens if _n == _N
				generate ind = 1
				bysort id: replace ind = 0 if t2 == adcens
				drop t0 t2
				// generate longitudinal outcome
				generate xb = `a0' + trt * `a1' + t * `a2' + v
				generate y = rnormal(xb, `=sqrt(`sigmae2')')
				// generate number of measurements per individual
				bysort id: generate nn = _N
				// generate cumulative number of measurements per individual
				bysort id: generate cumnn = _n
				// generate centered number of measurements per individual
				summarize nn if cumnn == 1
				generate cnn = nn - r(mean)
				// make Andersen-Gill counting process style data
				sort id t
				bysort id: generate t1 = t[_n+1]
				sort id t
				replace t1 = adcens if t1 == .
				// estimate normalised inverse intensity of visit weights
				stset t1, fail(ind) exit(time .) id(id) enter(t)
				stcox trt, nohr robust
				stset, clear
				predict lp, xb
				generate w = 1 / exp(lp)
				summarize w
				generate wn = w - r(mean) + 1
				bysort id: generate iivw = wn[_n - 1]
				replace iivw = 1 if iivw == .
				drop lp w wn t1
				// add gamma, lambda, DGM
				generate dgm = `dgm'
				generate i = `i'
				// compress
				compress
				// save dataset
				save "Data/df_`dgm'_`i'", replace
			}
			nois _dots `i' 0
		}
		local dgm = `dgm' + 1
	}
}

// Additional scenarios where the Gamma distribution depends on the previous value of Y AND treatment
clear
local shape = 2.0
local theta = 2.0
local prevY = 0.2
local nwide = 250
nois _dots 0, title(Gamma depending on previous values of Y) reps(`B')
forvalues i = 1/`B' {
	quietly {
		clear
		set obs `N'
		generate id = _n
		generate trt = runiform() > 0.5
		generate v = rnormal(0, `=sqrt(`sigmav2')')
		generate ge = rnormal(0, `=sqrt(`sigmage2')')
		generate adcens = runiform(5, 10)
		generate t0 = 0
		generate xb0 = `a0' + trt * `a1' + t0 * `a2' + v
		generate y0 = rnormal(xb0, `=sqrt(`sigmae2')')
		generate obtime0 = rgamma(`shape', exp(`theta' * (-1 * `beta') * trt + `prevY' * y0 + ge))
		forvalues j = 1(1)`nwide' {
			egen t`j' = rowtotal(obtime*) if `j' > 0
			generate xb`j' = `a0' + trt * `a1' + t`j' * `a2' + v
			generate y`j' = rnormal(xb`j', `=sqrt(`sigmae2')')
			generate obtime`j' = rgamma(`shape', exp(`theta' * (-1 * `beta') * trt + `prevY' * y`j' + ge))
		}
		reshape long xb y obtime t, i(id) j(obs)
		drop if t > adcens
		bysort id: generate t2 = t[_n + 1]
		generate ind = 1
		replace ind = 0 if t2 == .
		drop t2 obs
		// generate number of measurements per individual
		bysort id: generate nn = _N
		// generate cumulative number of measurements per individual
		bysort id: generate cumnn = _n
		// generate centered number of measurements per individual
		summarize nn if cumnn == 1
		generate cnn = nn - r(mean)
		// make Andersen-Gill counting process style data
		sort id t
		bysort id: generate t1 = t[_n+1]
		sort id t
		replace t1 = adcens if t1 == .
		// estimate normalised inverse intensity of visit weights
		stset t1, fail(ind) exit(time .) id(id) enter(t)
		stcox trt, nohr robust
		stset, clear
		predict lp, xb
		generate w = 1 / exp(lp)
		summarize w
		generate wn = w - r(mean) + 1
		bysort id: generate iivw = wn[_n - 1]
		replace iivw = 1 if iivw == .
		drop lp w wn t1
		// add gamma, lambda, DGM
		generate dgm = `dgm'
		generate i = `i'
		// compress
		compress
		// save dataset
		save "Data/df_`dgm'_`i'", replace
	}
	nois _dots `i' 0
}
local dgm = `dgm' + 1

// Final scenario: JM with sparse observation process + pre-defined visits each year
// df with predefined visit times
clear
set obs `N'
generate id = _n
expand 10
sort id
bysort id: generate t2 = _n
save _tmp, replace
// rest of the df
clear
foreach gamma in 3.00 {
	foreach lambda in 0.05 {
		nois _dots 0, title(Lambda: `lambda'; Gamma: `gamma') reps(`B')
		forvalues i = 1/`B' {
			quietly {
				clear
				set obs `N'
				generate id = _n
				generate trt = runiform() > 0.5
				generate u = rnormal(0, `=sqrt(`sigmau2')')
				generate v = rnormal(0, `=sqrt(`sigmav2')')
				// admin censoring time
				generate adcens = runiform(5, 10)
				// recurrent event/observation times
				expand `expN'
				sort id
				survsim obtime, lambdas(`lambda') gammas(`p') cov(trt `beta' u 1) distribution(weibull)
				bysort id: generate t2 = sum(obtime)
				// add predefined visits
				append using _tmp
				sort id t2
				bysort id: egen u_tmp = mean(u)
				replace u = u_tmp if u == .
				bysort id: egen v_tmp = mean(v)
				replace v = v_tmp if v == .
				bysort id: egen trt_tmp = mean(trt)
				replace trt = trt_tmp if trt == .
				bysort id: egen adcens_tmp = mean(adcens)
				replace adcens = adcens_tmp if adcens == .
				// generate time variables
				generate t0 = 1
				bysort id: replace t0 = 0 if _n > 1
				generate t = 0
				bysort id: replace t = t2[_n - 1] if _n > 1
				drop if t > adcens
				replace obtime = t2 - t if obtime == .
				bysort id: replace t2 = adcens if _n == _N
				generate ind = 1
				bysort id: replace ind = 0 if t2 == adcens
				drop t0 t2 u_tmp v_tmp trt_tmp adcens_tmp
				// generate longitudinal outcome
				generate xb = `a0' + trt * `a1' + t * `a2' + `gamma' * u + v
				generate y = rnormal(xb, `=sqrt(`sigmae2')')
				// generate number of measurements per individual
				bysort id: generate nn = _N
				// generate cumulative number of measurements per individual
				bysort id: generate cumnn = _n
				// generate centered number of measurements per individual
				summarize nn if cumnn == 1
				generate cnn = nn - r(mean)
				// make Andersen-Gill counting process style data
				sort id t
				bysort id: generate t1 = t[_n+1]
				sort id t
				replace t1 = adcens if t1 == .
				// estimate normalised inverse intensity of visit weights
				stset t1, fail(ind) exit(time .) id(id) enter(t)
				stcox trt, nohr robust
				stset, clear
				predict lp, xb
				generate w = 1 / exp(lp)
				summarize w
				generate wn = w - r(mean) + 1
				bysort id: generate iivw = wn[_n - 1]
				replace iivw = 1 if iivw == .
				drop lp w wn t1
				// add gamma, lambda, DGM
				generate dgm = `dgm'
				generate i = `i'
				// compress
				compress
				// save dataset
				save "Data/df_`dgm'_`i'", replace
			}
			nois _dots `i' 0
		}
		local dgm = `dgm' + 1
	}
}
