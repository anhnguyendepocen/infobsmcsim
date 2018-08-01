clear

// Simulate data
local B = 1000
local N = 200
local a0 = 0
local a1 = 1
local a2 = 0.2
local sigmau2 = 1
local sigmav2 = 0.5
local sigmae2 = 1
local p = 1.05
local beta = 1
local dgm = 1
set seed 705862548
foreach gamma in 0 1.50 {
	foreach lambda in 0.10 0.30 1.00 {
		nois _dots 0, title(Lambda: `lambda'; Gamma: `gamma') reps(`B')
		forvalues i = 1/`B' {
			quietly {
				clear
				set obs `N'
				gen id = _n
				gen trt = runiform() > 0.5
				gen u = rnormal(0, `=sqrt(`sigmau2')')
				gen v = rnormal(0, `=sqrt(`sigmav2')')
				// admin censoring time
				gen adcens = runiform(5, 10)
				// recurrent event/observation times
				expand 2000
				sort id
				generate t0 = 1
				bysort id: replace t0 = 0 if _n > 1
				survsim obtime, lambdas(`lambda') gammas(`p') cov(trt `beta' u 1) distribution(weibull)
				generate t = 0
				bysort id: gen double t2 = sum(obtime)
				bysort id: replace t = t2[_n - 1] if _n > 1
				drop if t > adcens
				bysort id: replace t2 = adcens if _n == _N
				generate ind = 1
				bysort id: replace ind = 0 if t2 == adcens
				drop t0 t2
				// generate longitudinal outcome
				gen double xb = `a0' + trt*`a1' + t * `a2' + `gamma' * u + v
				gen double y = rnormal(xb, `=sqrt(`sigmae2')')
				// gen number of measurements per individual
				bysort id: gen nn = _N
				// gen cumulative number of measurements per individual
				bysort id: gen cumnn = _n
				// gen centered number of measurements per individual
				summ nn if cumnn == 1
				gen cnn = nn - r(mean)
				// make Andersen-Gill counting process style data
				sort id t
				bysort id: gen t1 = t[_n+1]
				sort id t
				replace t1 = adcens if t1 == .
				// estimate normalised inverse intensity of visit weights
				stset t1, fail(ind) exit(time .) id(id) enter(t)
				stcox trt, nohr vce(cluster id)
				stset, clear
				predict lp, xb
				gen w = 1 / exp(lp)
				summ w
				gen wn = w - r(mean) + 1
				bysort id: gen iivw = wn[_n - 1]
				replace iivw = 1 if iivw == .
				drop lp w wn t1
				// add gamma, lambda, DGM
				gen gamma = `gamma'
				gen lambda = `lambda'
				gen dgm = `dgm'
				gen i = `i'
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
