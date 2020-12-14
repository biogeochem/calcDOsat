#' Solubility Equations for Oxygen in Water
#'
#' This function allows you to calculate dissolved oxygen from USGS Technical Memorandum 2011.03
#' @param dat dataframe with your values
#' @param method choose from either Weiss ('weiss'), Benson-Krause ('benson_krause'), or Garcia-Gordon ('garcia_gordon')
#' @param water_temp_C column name with sample water temperature in C
#' @param baro_P_kPa column name with atmospheric pressure at time of sampling
#' @param spec_cond_uS.cm column name with sample specific conductance in uS/cm
#' @param meas_do_mg.L column name with sample water dissolved oxygen value in mg/L
#' @keywords dissolved oxygen, USGS, percent saturation
#' @export
#' @examples
#' dat <- data.frame(Temp_C = 20, Baro_kPa = 101.325, Spec_Cond_uS.cm = 0, DO_mg.L = 8)
#' calcDOsat(dat, 'weiss', 'Temp_C', 'Baro_kPa', 'Spec_Cond_uS.cm', 'DO_mg.L')


calcDOsat <- function(dat, method, water_temp_C, baro_P_kPa, spec_cond_uS.cm, meas_do_mg.L){
  #Conversion Calculations
  temp_K <- dat[[water_temp_C]] + 273.15;
  press_mmHg <- dat[[baro_P_kPa]] * 7.50061575;
  press_atm <- dat[[baro_P_kPa]] * 0.00986923;
  S_ppt <- ((5.572*10^-4) * dat[[spec_cond_uS.cm]]) + ((2.02*10^-9) * dat[[spec_cond_uS.cm]]^2);
  u_weiss <- 10^(8.10765 - (1750.286/(235 + dat[[water_temp_C]])) );
  theta <- 0.000975 - ((1.426*10^-5) * dat[[water_temp_C]]) + ((6.436*10^-8) * (dat[[water_temp_C]]^2));
  u_bk <- exp(11.8571 - (3840.7/temp_K) - (216961 / (temp_K^2)) );
  Ts <- log((298.15 - dat[[water_temp_C]]) / (273.15 + dat[[water_temp_C]]));

  if(method == 'weiss') {
    term_DO <- 1.42905 * exp(-173.4292 + (249.6339 * (100/temp_K)) + (143.3483 * (log(temp_K/100))) - (21.8492 * (temp_K/100)) );
    term_Fs <- exp(S_ppt * (-0.033096 + (0.014259 * (temp_K/100)) - (0.0017000 * ((temp_K/100)^2)) ) );
    term_Fp <-  ((press_mmHg - u_weiss) / (760 - u_weiss)) ;
    DO_sat <- term_DO * term_Fs * term_Fp  ;
    dat$calc_DO_sat_perc <- (dat[[meas_do_mg.L]] / DO_sat) * 100 ;

    return(dat)

  } else {

    if(method == 'benson_krause') {
      term_DO <- exp(-139.34411 + ((1.575701*10^5)/temp_K) - ((6.642308*10^7)/(temp_K^2)) + ((1.243800*10^10)/(temp_K^3)) - ((8.621949*10^11)/(temp_K^4)) );
      term_Fs <- exp(-S_ppt * (0.017674 - (10.754/temp_K) + (2140.7/(temp_K^2)) ) );
      term_Fp <- ((press_atm - u_bk) * (1 - (theta * press_atm))) / ((1 - u_bk) * (1 - theta));
      DO_sat <- term_DO * term_Fs * term_Fp  ;
      dat$calc_DO_sat_perc <- (dat[[meas_do_mg.L]] / DO_sat) * 100 ;

      return(dat)

    } else {
      if(method == 'garcia_gordon'){
        term_DO <- 1.42905 * exp(2.00907 + (3.22014 * Ts) + (4.05010 * (Ts^2)) + (4.94457 * (Ts^3)) - (0.256847 * (Ts^4)) + (3.88767 * (Ts^5)) );
        term_Fs <- exp( ( (-0.00624523 - (0.00737614 * Ts) - (0.0103410 * (Ts^2)) - (0.00817083 * (Ts^3)) ) * S_ppt) - ((4.88682*10^-7)*(S_ppt^2)) );
        term_Fp <- ((press_atm - u_bk) * (1 - (theta * press_atm))) / ((1 - u_bk) * (1 - theta)) ;
        DO_sat <- term_DO * term_Fs * term_Fp  ;
        dat$calc_DO_sat_perc <- (dat[[meas_do_mg.L]] / DO_sat) * 100 ;

        return(dat)
      }
    }
  }
}


