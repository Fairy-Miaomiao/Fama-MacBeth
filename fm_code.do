  * 读取25个投资组合收益率数据
  insheet using "Rp.csv", noname clear
  drop if _n==1 // 删除原表头
  list in 1/5
  rename v1 t
  forvalues i=2(1)26{
	local j = `i'-1
	destring v`i', replace
    rename v`i' rp`j'
	replace  rp`j' =  rp`j'/100 //原始数据是百分号前面的部分
  } // rpi对应25个投资组合收益率，调编号
  save Rp_temp, replace
  // use Rp_temp, replace
  * 读取市场组合风险溢价数据
  insheet using "Mkt-RF.csv", name clear
  rename v1 t
  list in 1/5
  tostring t, replace // 转为字符串格式，为了merge时数据类型统一
  replace mktrf = mktrf/100 //原始数据是百分号前面的部分
  replace rf = rf/100 //原始数据是百分号前面的部分
  drop smb hml
  save RMkt_temp,replace
  
  * 拼接数据
  // use RMkt_temp, clear
  merge 1:1 t using Rp_temp.dta
  keep if _merge==3
  drop _merge
  forvalues i=1(1)25{
	gen rpe`i' = rp`i'-rf // 投资组合超额收益
  }
  list in 1/5
  save merged, replace
  
  // use merged, replace
  gen i = _n
  reshape long rpe, i(i) j(port_num) //reshape成long型数据
  rename i time_index 
  keep port_num t rpe mktrf 
  sort port_num t rpe mktrf 
  
  gen t_ = substr(t,1,4)+"-"+substr(t,5,6)
  destring t, replace
  replace t = monthly(t_,"YM")
  format t %tm //日期格式
  drop t_
  
  xtset port_num t 
  label var port_num "portfolio number"
  label var mktrf "excess market return"
  label var rpe "excess portfolio return"
  save long_, replace


  *******************************************************************
  * pass1 1930.1-1938.11：25*48次时序回归 （1930.1-1934.12->1933.12-1938.11）
  * 估计beta_it i=1,2...25，窗口为五年，每次向后移动一个月
  use long_, clear
  format mktrf rpe %4.3f  
  list in 1131/1135
  //winsor2 mktrf rpe, replace cuts(1 99) //缩尾处理
  tabstat mktrf rpe if (t>=ym(1930,1) & t<=ym(1938,12)) ,stats(mean sd min max p1 p25 p50 p75 p99) column(s) long  format(%4.3f)

  bys port_num: asreg rp mktrf if (t>=ym(1930,1) & t<=ym(1938,12)) , wind(t 60) rmse se newey(4) 
  //分组时序回归，每个投资组合为一组
  // 窗口60月；保存rmse se；回报 Newey&West 调整标准误
  
  keep if _Nobs==60 // 保留样本数为60的
  format mktrf rpe _rmse _R2 _adjR2 _b_mktrf _b_cons _se_mktrf    _se_cons %4.2f  
  format _Nobs %4.0f
  list in 46/50
  
   tabstat  _b_mktrf _se_mktrf _R2 _rmse, by(port_num) stat(mean) format(%4.2f) 
   // pass1 回归结果在时序上的均值

  rename _b_mktrf beta
  keep port_num t rpe beta 
  // save pass1, replace
  bys port_num: gen Lbeta = L.beta //为了截面回归更方便，直接将自变量取滞后项
  label var Lbeta "beta_i,t-1" 
  sort port_num t rpe Lbeta 
  save pass1,replace
  
  *******************************************************************
  * pass2 1935.1-1938.12：48次截面回归
  use pass1, clear
  drop if (t==ym(1934,12))
  tsset port_num t  
  // xtset port_num t  
  
	* before regression
  local c1 "if t==ym(1935,1)"
  local c2 "if t==ym(1938,1)"
  twoway (scatter rpe Lbeta `c1')(lfit rpe Lbeta `c1')   ///
  (scatter rpe Lbeta `c2')(lfit rpe Lbeta `c2'),   ///
	xtitle("Estimated Beta 1 Month Before (Lbeta)") ytitle("Portfolio Excess Return (rpe)") /// 
	legend(label(1 "1935m1") label(2 "1935m1 fitted") label(3 "1938m1") label(4 "1938m1 fitted")) ///
	note("Data Source: Ken French Data Library")
  
  format rpe beta Lbeta %4.2f
  list in 47/51

  global regvar "rpe Lbeta"
	*xtfmb
	xtfmb $regvar
		est store m1
		
	xtfmb $regvar, verbose lag(4) 
	// verbose显示每期回归的lambda_t
	//考虑四阶自相关的newey&west correction
		est store m2
	
	*asreg
	asreg $regvar, fmb
		estadd scalar r2_a = e(adjr2) 
		// asreg fmb汇报的adjr2需要手动添加
		est store m3
			
	asreg $regvar, fmb first newey(4) save(FirstStage) 
	// 考虑四阶自相关的newey&west correction
	// 保存lambda_t至FirstStage.dta  		
		estadd scalar r2_a = e(adjr2)
		est store m4
		
	local mt "xtfmb xtfmb_newey asreg asreg_newey"
	local m "m1 m2 m3 m4"
	local t "Fama-Macbeth second pass regression"
	local n "Data Source: Ken French Data Library"

	esttab `m', title(`t') mtitle(`mt') star(* 0.1 ** 0.05 *** 0.01) ///
	r2 ar2 label nogap addnote(`n') b(%4.3f) t(%4.3f)

		
	*通过asreg FirstStage数据获得和表格一样的结果
	use FirstStage.dta, clear 
	format _b_Lbeta _Cons _R2 _adjR2 %4.2f 
	list in 1/5
	tabstat _b _C _R2 _ad , stats(mean sd) format(%4.3f )
	//描述性统计,_b均值和估计结果一致
	
	qui sum _R2 //_R2 95%置信区间
	disp %4.3f r(mean)-1.96*r(sd)/sqrt(r(N)) //.20344001
	disp %4.3f r(mean)+1.96*r(sd)/sqrt(r(N)) //.31000475

	local vars "_b_Lbeta"
	foreach v of varlist `vars'{
		qui sum `v'
		scalar lambda_`v' = r(mean) 
		// 每期斜率估计量均值
		scalar se_`v' = r(sd)/sqrt(r(N))
		scalar t_`v' = lambda_`v'/se_`v'
		disp in green "Second Pass result of " in yellow "`v'"
		disp %4.3f in green "coefficient "  in yellow lambda_`v'
		disp %4.3f in green "standard error "  in yellow  se_`v'
		disp %4.3f in green "t statistics "   in yellow t_`v'
	}
	
	* 分组回归：以每个时期为一组
	use pass1, clear
	statsby _b _se, by(t): reg $regvar 
	// 每期做一次截面回归
	local vars "_b_Lbeta"
	foreach v of varlist `vars'{
		qui sum `v'
		scalar lambda_`v' = r(mean)
		scalar se_`v' = r(sd)/sqrt(r(N))
		scalar t_`v' = lambda_`v'/se_`v'
		disp in green "Second Pass result of " in yellow "`v'"
		disp in green "coefficient "  in yellow %4.3f lambda_`v'
		disp in green "standard error "  in yellow %4.3f  se_`v'
		disp in green "t statistics "   in yellow %4.3f t_`v'
	}
	
	* 关于fm
	  * pass1：和statsby对各投资组合分组回归等价	
	use long_, replace
	tsset port_num t 
	fm rp mktrf if t >= ym(1934,12) & t <= ym(1938,12), byfm(port_num)
	ereturn list
	save temp_fm , replace
	statsby _b _se, by(port_num): reg rp mktrf if t >= ym(1934,12) & t <= ym(1938,12)
	save statsby_pass1, replace
	sum _b_mktrf
	scalar _b_mean = r(mean)
	scalar _b_se_mean = r(se)
	
	  * pass2
	use pass1, clear
    tsset port_num t 
	fm $regvar, byfm(t) // R2计算有问题
	

