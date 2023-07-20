prepareInitials <- function(init_beta,
                            init_scale_horizontal, init_scale_vertical,
                            init_a1, init_b1, init_c1_coef, init_d1,
                            init_a2, init_b2, init_c2_coef, init_d2,
                            beta_fix, scale_horizontal_fix,
                            scale_vertical_fix,
                            a1_fix, b1_fix, c1_fix, d1_fix,
                            a2_fix, b2_fix, c2_fix, d2_fix){

  if(c1_fix & c2_fix){
    if(d1_fix & d2_fix){

      theta0 <- c(init_beta, init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_a2, init_b2)

    }
  }else{
    if(!beta_fix){
      if(a1_fix & b1_fix & a2_fix & b2_fix){
        if(scale_horizontal_fix & scale_vertical_fix){
          if(d1_fix & d2_fix){

            theta0 <- c(init_beta, init_c1_coef, init_c2_coef)

          }else{

            theta0 <- c(init_beta, init_c1_coef, init_d1, init_c2_coef, init_d2)

          }
        }else{
          if(d1_fix & d2_fix){

            theta0 <- c(init_beta, init_scale_horizontal, init_scale_vertical, init_c1_coef, init_c2_coef)

          }else{

            theta0 <- c(init_beta, init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_c1_coef, init_d1, init_a2, init_b2, init_c2_coef, init_d2)

          }
        }
      }else{
        if(scale_horizontal_fix & scale_vertical_fix){
          if(d1_fix & d2_fix){

            theta0 <- c(init_beta, init_a1, init_b1, init_c1_coef, init_a2, init_b2, init_c2_coef)

          }else{

            theta0 <- c(init_beta, init_a1, init_b1, init_c1_coef, init_d1, init_a2, init_b2, init_c2_coef, init_d2)

          }
        }else if(!scale_horizontal_fix & !scale_vertical_fix){
          if(d1_fix & d2_fix){

            theta0 <- c(init_beta, init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_c1_coef, init_a2, init_b2, init_c2_coef)

          }else{

            theta0 <- c(init_beta, init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_c1_coef, init_d1, init_a2, init_b2, init_c2_coef, init_d2)

          }
        }
      }
    }else{
      if(a1_fix & b1_fix & a2_fix & b2_fix){
        if(scale_horizontal_fix & scale_vertical_fix){
          if(d1_fix & d2_fix){

            theta0 <- c(init_c1_coef, init_c2_coef)

          }else{

            theta0 <- c(init_c1_coef, init_d1, init_c2_coef, init_d2)

          }
        }else{
          if(d1_fix & d2_fix){

            theta0 <- c(init_scale_horizontal, init_scale_vertical, init_c1_coef, init_c2_coef)

          }else{

            theta0 <- c(init_scale_horizontal, init_scale_vertical, init_c1_coef, init_d1, init_c2_coef, init_d2)

          }
        }
      }else{
        if(scale_horizontal_fix & scale_vertical_fix){
          if(d1_fix & d2_fix){

            theta0 <- c(init_a1, init_b1, init_c1_coef, init_a2, init_b2, init_c2_coef)

          }else{

            theta0 <- c(init_a1, init_b1, init_c1_coef, init_d1, init_a2, init_b2, init_c2_coef, init_d2)

          }
        }else{
          if(d1_fix & d2_fix){

            theta0 <- c(init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_c1_coef, init_a2, init_b2, init_c2_coef)

          }else{

            theta0 <- c(init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_c1_coef, init_d1, init_a2, init_b2, init_c2_coef, init_d2)

          }
        }
      }
    }
  }
  return(theta0)
}

