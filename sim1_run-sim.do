// Set working directory
cd "Z:\jointmodels\projects\Informative observation process"
clear

// Local values
local totdgm = 10
local B = 1000

// Create postfile
postfile sim1 i model dgm est b se error using sim1, replace

// Do not put results of estimation commands in r(.)
set coeftabresults off

forvalues dgm = 1/`totdgm' {
	nois _dots 0, title(DGM: `dgm') reps(`B')
	forvalues i = 1/`B' {
		quietly {
			use "Data/df_`dgm'_`i'", replace
			capture gsem	(obtime <- trt M1[id]@1, family(weibull, failure(ind)))	///
				(y <- trt t M1[id] M2[id]@1, family(gaussian)), 					///
				covstruct(M1[id] M2[id], diag)
			if (_rc > 0 ) {
				// gsem failed
				post sim1 (`i') (1) (`dgm') (1) (.) (.) (99)
				post sim1 (`i') (1) (`dgm') (2) (.) (.) (99)
				post sim1 (`i') (1) (`dgm') (3) (.) (.) (99)
				post sim1 (`i') (1) (`dgm') (4) (.) (.) (99)
				post sim1 (`i') (1) (`dgm') (5) (.) (.) (99)
				post sim1 (`i') (1) (`dgm') (6) (.) (.) (99)
				post sim1 (`i') (1) (`dgm') (7) (.) (.) (99)
			}
			else {
				// gsem ok
				post sim1 (`i') (1) (`dgm') (1) ([y]_b[_cons]) ([y]_se[_cons]) (0)
				post sim1 (`i') (1) (`dgm') (2) ([y]_b[trt]) ([y]_se[trt]) (0)
				post sim1 (`i') (1) (`dgm') (3) ([y]_b[t]) ([y]_se[t]) (0)
				post sim1 (`i') (1) (`dgm') (4) ([/]_b[var(M1[id])]) ( [/]_se[var(M1[id])]) (0)
				post sim1 (`i') (1) (`dgm') (5) ([/]_b[var(M2[id])]) ( [/]_se[var(M2[id])]) (0)
				post sim1 (`i') (1) (`dgm') (6) ([/]_b[var(e.y)]) ([/]_se[var(e.y)]) (0)
				post sim1 (`i') (1) (`dgm') (7) ([y]_b[M1[id]]) ([y]_se[M1[id]]) (0)
			}
			// fit mixed model adjusting for centered number of measurements
			capture mixed y trt t cnn || id:
			if (_rc > 0) {
				// mixed failed
				post sim1 (`i') (2) (`dgm') (1) (.) (.) (99)
				post sim1 (`i') (2) (`dgm') (2) (.) (.) (99)
				post sim1 (`i') (2) (`dgm') (3) (.) (.) (99)
				post sim1 (`i') (2) (`dgm') (4) (.) (.) (99)
				post sim1 (`i') (2) (`dgm') (5) (.) (.) (99)
				post sim1 (`i') (2) (`dgm') (6) (.) (.) (99)
				post sim1 (`i') (2) (`dgm') (7) (.) (.) (99)
			}
			else {
				// mixed ok
				post sim1 (`i') (2) (`dgm') (1) ([y]_b[_cons]) ([y]_se[_cons]) (0)
				post sim1 (`i') (2) (`dgm') (2) ([y]_b[trt]) ([y]_se[trt]) (0)
				post sim1 (`i') (2) (`dgm') (3) ([y]_b[t]) ([y]_se[t]) (0)
				post sim1 (`i') (2) (`dgm') (4) (.) (.) (0)
				post sim1 (`i') (2) (`dgm') (5) (exp([lns1_1_1]_b[_cons]) ^ 2) (sqrt([lns1_1_1]_se[_cons] ^ 2 * (2 * exp([lns1_1_1]_b[_cons]) ^ 2) ^ 2)) (0)
				post sim1 (`i') (2) (`dgm') (6) (exp([lnsig_e]_b[_cons]) ^ 2) (sqrt([lnsig_e]_se[_cons] ^ 2 * (2 * exp([lnsig_e]_b[_cons]) ^ 2) ^ 2)) (0)
				post sim1 (`i') (2) (`dgm') (7) (.) (.) (0)

			}
			// fit mixed model adjusting for cumulative number of measurements
			capture mixed y trt t cumnn || id:
			if (_rc > 0) {
				// mixed failed
				post sim1 (`i') (3) (`dgm') (1) (.) (.) (99)
				post sim1 (`i') (3) (`dgm') (2) (.) (.) (99)
				post sim1 (`i') (3) (`dgm') (3) (.) (.) (99)
				post sim1 (`i') (3) (`dgm') (4) (.) (.) (99)
				post sim1 (`i') (3) (`dgm') (5) (.) (.) (99)
				post sim1 (`i') (3) (`dgm') (6) (.) (.) (99)
				post sim1 (`i') (3) (`dgm') (7) (.) (.) (99)
			}
			else {
				// mixed ok
				post sim1 (`i') (3) (`dgm') (1) ([y]_b[_cons]) ([y]_se[_cons]) (0)
				post sim1 (`i') (3) (`dgm') (2) ([y]_b[trt]) ([y]_se[trt]) (0)
				post sim1 (`i') (3) (`dgm') (3) ([y]_b[t]) ([y]_se[t]) (0)
				post sim1 (`i') (3) (`dgm') (4) (.) (.) (0)
				post sim1 (`i') (3) (`dgm') (5) (exp([lns1_1_1]_b[_cons]) ^ 2) (sqrt([lns1_1_1]_se[_cons] ^ 2 * (2 * exp([lns1_1_1]_b[_cons]) ^ 2) ^ 2)) (0)
				post sim1 (`i') (3) (`dgm') (6) (exp([lnsig_e]_b[_cons]) ^ 2) (sqrt([lnsig_e]_se[_cons] ^ 2 * (2 * exp([lnsig_e]_b[_cons]) ^ 2) ^ 2)) (0)
				post sim1 (`i') (3) (`dgm') (7) (.) (.) (0)

			}
			// fit standard mixed model with iiv weights
			capture mixed y trt t || id: [pw = iivw]
			if (_rc > 0) {
				// mixed failed
				post sim1 (`i') (4) (`dgm') (1) (.) (.) (99)
				post sim1 (`i') (4) (`dgm') (2) (.) (.) (99)
				post sim1 (`i') (4) (`dgm') (3) (.) (.) (99)
				post sim1 (`i') (4) (`dgm') (4) (.) (.) (99)
				post sim1 (`i') (4) (`dgm') (5) (.) (.) (99)
				post sim1 (`i') (4) (`dgm') (6) (.) (.) (99)
				post sim1 (`i') (4) (`dgm') (7) (.) (.) (99)
			}
			else {
				// mixed ok
				post sim1 (`i') (4) (`dgm') (1) ([y]_b[_cons]) ([y]_se[_cons]) (0)
				post sim1 (`i') (4) (`dgm') (2) ([y]_b[trt]) ([y]_se[trt]) (0)
				post sim1 (`i') (4) (`dgm') (3) ([y]_b[t]) ([y]_se[t]) (0)
				post sim1 (`i') (4) (`dgm') (4) (.) (.) (0)
				post sim1 (`i') (4) (`dgm') (5) (exp([lns1_1_1]_b[_cons]) ^ 2) (sqrt([lns1_1_1]_se[_cons] ^ 2 * (2 * exp([lns1_1_1]_b[_cons]) ^ 2) ^ 2)) (0)
				post sim1 (`i') (4) (`dgm') (6) (exp([lnsig_e]_b[_cons]) ^ 2) (sqrt([lnsig_e]_se[_cons] ^ 2 * (2 * exp([lnsig_e]_b[_cons]) ^ 2) ^ 2)) (0)
				post sim1 (`i') (4) (`dgm') (7) (.) (.) (0)
			}
			// fit standard mixed model
			capture mixed y trt t || id:
			if (_rc > 0) {
				// mixed failed
				post sim1 (`i') (5) (`dgm') (1) (.) (.) (99)
				post sim1 (`i') (5) (`dgm') (2) (.) (.) (99)
				post sim1 (`i') (5) (`dgm') (3) (.) (.) (99)
				post sim1 (`i') (5) (`dgm') (4) (.) (.) (99)
				post sim1 (`i') (5) (`dgm') (5) (.) (.) (99)
				post sim1 (`i') (5) (`dgm') (6) (.) (.) (99)
				post sim1 (`i') (5) (`dgm') (7) (.) (.) (99)
			}
			else {
				// mixed ok
				post sim1 (`i') (5) (`dgm') (1) ([y]_b[_cons]) ([y]_se[_cons]) (0)
				post sim1 (`i') (5) (`dgm') (2) ([y]_b[trt]) ([y]_se[trt]) (0)
				post sim1 (`i') (5) (`dgm') (3) ([y]_b[t]) ([y]_se[t]) (0)
				post sim1 (`i') (5) (`dgm') (4) (.) (.) (0)
				post sim1 (`i') (5) (`dgm') (5) (exp([lns1_1_1]_b[_cons]) ^ 2) (sqrt([lns1_1_1]_se[_cons] ^ 2 * (2 * exp([lns1_1_1]_b[_cons]) ^ 2) ^ 2)) (0)
				post sim1 (`i') (5) (`dgm') (6) (exp([lnsig_e]_b[_cons]) ^ 2) (sqrt([lnsig_e]_se[_cons] ^ 2 * (2 * exp([lnsig_e]_b[_cons]) ^ 2) ^ 2)) (0)
				post sim1 (`i') (5) (`dgm') (7) (.) (.) (0)
			}
			// fit gee with iiv weights
			capture glm y trt t [pw = iivw], family(gaussian) link(identity) vce(cluster id)
			if (_rc > 0) {
				// glm failed
				post sim1 (`i') (6) (`dgm') (1) (.) (.) (99)
				post sim1 (`i') (6) (`dgm') (2) (.) (.) (99)
				post sim1 (`i') (6) (`dgm') (3) (.) (.) (99)
				post sim1 (`i') (6) (`dgm') (4) (.) (.) (99)
				post sim1 (`i') (6) (`dgm') (5) (.) (.) (99)
				post sim1 (`i') (6) (`dgm') (6) (.) (.) (99)
				post sim1 (`i') (6) (`dgm') (7) (.) (.) (99)
			}
			else {
				// glm ok
				post sim1 (`i') (6) (`dgm') (1) ([y]_b[_cons]) ([y]_se[_cons]) (0)
				post sim1 (`i') (6) (`dgm') (2) ([y]_b[trt]) ([y]_se[trt]) (0)
				post sim1 (`i') (6) (`dgm') (3) ([y]_b[t]) ([y]_se[t]) (0)
				post sim1 (`i') (6) (`dgm') (4) (.) (.) (0)
				post sim1 (`i') (6) (`dgm') (5) (.) (.) (0)
				post sim1 (`i') (6) (`dgm') (6) (.) (.) (0)
				post sim1 (`i') (6) (`dgm') (7) (.) (.) (0)
			}
			// fit gee without iiv weights
			capture glm y trt t, family(gaussian) link(identity) vce(cluster id)
			if (_rc > 0) {
				// glm failed
				post sim1 (`i') (7) (`dgm') (1) (.) (.) (99)
				post sim1 (`i') (7) (`dgm') (2) (.) (.) (99)
				post sim1 (`i') (7) (`dgm') (3) (.) (.) (99)
				post sim1 (`i') (7) (`dgm') (4) (.) (.) (99)
				post sim1 (`i') (7) (`dgm') (5) (.) (.) (99)
				post sim1 (`i') (7) (`dgm') (6) (.) (.) (99)
				post sim1 (`i') (7) (`dgm') (7) (.) (.) (99)
			}
			else {
				// glm ok
				post sim1 (`i') (7) (`dgm') (1) ([y]_b[_cons]) ([y]_se[_cons]) (0)
				post sim1 (`i') (7) (`dgm') (2) ([y]_b[trt]) ([y]_se[trt]) (0)
				post sim1 (`i') (7) (`dgm') (3) ([y]_b[t]) ([y]_se[t]) (0)
				post sim1 (`i') (7) (`dgm') (4) (.) (.) (0)
				post sim1 (`i') (7) (`dgm') (5) (.) (.) (0)
				post sim1 (`i') (7) (`dgm') (6) (.) (.) (0)
				post sim1 (`i') (7) (`dgm') (7) (.) (.) (0)
			}
		}
		nois _dots `i' 0
	}
}
postclose sim1
use sim1, clear
compress
label define model 1 "true" 2 "me(cnn)" 3 "me(cumnn)" 4 "me(iivw)" 5 "me" 6 "gee(iivw)" 7 "gee"
label values model model
label define est 1 "a0" 2 "a1" 3 "a2" 4 "v(u)" 5 "v(v)" 6 "v(e)" 7 "gamma"
label values est est
save, replace

// summarise simulations
simsum b if (est == 1), id(i) se(se) true(0) methodvar(model) mcse by(dgm) bias cover
simsum b if (est == 2), id(i) se(se) true(1) methodvar(model) mcse by(dgm) bias cover
simsum b if (est == 3), id(i) se(se) true(0.2) methodvar(model) mcse by(dgm) bias cover
simsum b if (est == 4), id(i) se(se) true(1) methodvar(model) mcse by(dgm) bias cover
simsum b if (est == 5), id(i) se(se) true(0.5) methodvar(model) mcse by(dgm) bias cover
simsum b if (est == 6), id(i) se(se) true(1) methodvar(model) mcse by(dgm) bias cover
simsum b if (est == 7 & (dgm == 1 | dgm == 2 | dgm == 3)), id(i) se(se) true(0) methodvar(model) mcse by(dgm) bias cover
simsum b if (est == 7 & (dgm == 4 | dgm == 5 | dgm == 6)), id(i) se(se) true(1.5) methodvar(model) mcse by(dgm) bias cover
