prepareResults_mle <- function(fittedModel,
                      BETA, SCALE_HORIZONTAL, SCALE_VERTICAL,
                      A1, B1, C1_coef, D1, A2, B2, C2_coef, D2,
                      est_sd,
                      beta_fix, scale_horizontal_fix,
                      scale_vertical_fix,
                      a1_fix, b1_fix, c1_fix, d1_fix,
                      a2_fix, b2_fix, c2_fix, d2_fix,
                      splines_degree,
                      inner_knots1, inner_knots2, nb1, nb2){

  if(c1_fix & c2_fix){
    if(d1_fix & d2_fix){

      results <- list(convergence_code = fittedModel$code, iterations = fittedModel$iterations, loglikelihood_value = -fittedModel$minimum, theta = fittedModel$estimate, gradient = fittedModel$gradient,
                      est_beta = BETA, est_beta_sd = est_sd[1], est_scale_horizontal = SCALE_HORIZONTAL, est_scale_horizontal_sd = est_sd[2],
                      est_scale_vertical = SCALE_VERTICAL, est_scale_vertical_sd = est_sd[3], est_a1 = A1, est_a1_sd = est_sd[4], est_b1 = B1, est_b1_sd = est_sd[5],
                      est_c1_coef = C1_coef, est_d1 = D1, est_a2 = A2, est_a2_sd = est_sd[6], est_b2 = B2, est_b2_sd = est_sd[7], est_c2_coef = C2_coef, est_d2 = D2,
                      splines_degree = splines_degree, inner_knots1 = inner_knots1, inner_knots2 = inner_knots2)
    }
  }else{
    if(!beta_fix){
      if(a1_fix & b1_fix & a2_fix & b2_fix){
        if(scale_horizontal_fix & scale_vertical_fix){
          if(d1_fix & d2_fix){

          }else{

          }
        }else{
          if(d1_fix & d2_fix){

            results <- list(convergence_code = fittedModel$code, iterations = fittedModel$iterations, loglikelihood_value = -fittedModel$minimum, theta = fittedModel$estimate, gradient = fittedModel$gradient,
                            est_beta = BETA, est_beta_sd = est_sd[1], est_scale_horizontal = SCALE_HORIZONTAL, est_scale_horizontal_sd = est_sd[2],
                            est_scale_vertical = SCALE_VERTICAL, est_scale_vertical_sd = est_sd[3], est_a1 = A1, est_b1 = B1,
                            est_c1_coef = C1_coef, est_c1_coef_sd = est_sd[3 + 1:nb1], est_d1 = D1, est_a2 = A2, est_b2 = B2,
                            est_c2_coef = C2_coef, est_c2_coef_sd = est_sd[3 + nb1 + 1:nb2], est_d2 = D2,
                            splines_degree = splines_degree, inner_knots1 = inner_knots1, inner_knots2 = inner_knots2)

          }else{

          }
        }
      }else{
        if(scale_horizontal_fix & scale_vertical_fix){
          if(d1_fix & d2_fix){

            results <- list(convergence_code = fittedModel$code, iterations = fittedModel$iterations, loglikelihood_value = -fittedModel$minimum, theta = fittedModel$estimate, gradient = fittedModel$gradient,
                            est_beta = BETA, est_beta_sd = est_sd[1], est_scale_horizontal = SCALE_HORIZONTAL, est_scale_vertical = SCALE_VERTICAL,
                            est_a1 = A1, est_a1_sd = est_sd[2], est_b1 = B1, est_b1_sd = est_sd[3],
                            est_c1_coef = C1_coef, est_c1_coef_sd = est_sd[3 + 1:nb1], est_d1 = D1, est_a2 = A2, est_a2_sd = est_sd[3 + nb1 + 1],
                            est_b2 = B2, est_b2_sd = est_sd[3 + nb1 + 2], est_c2_coef = C2_coef, est_c2_coef_sd = est_sd[3 + nb1 + 2 + 1:nb2], est_d2 = D2,
                            splines_degree = splines_degree, inner_knots1 = inner_knots1, inner_knots2 = inner_knots2)

          }else{

          }
        }else if(!scale_horizontal_fix & !scale_vertical_fix){
          if(d1_fix & d2_fix){

            results <- list(convergence_code = fittedModel$code, iterations = fittedModel$iterations, loglikelihood_value = -fittedModel$minimum, theta = fittedModel$estimate, gradient = fittedModel$gradient,
                            est_beta = BETA, est_beta_sd = est_sd[1], est_scale_horizontal = SCALE_HORIZONTAL, est_scale_horizontal_sd = est_sd[2],
                            est_scale_vertical = SCALE_VERTICAL, est_scale_vertical_sd = est_sd[3],
                            est_a1 = A1, est_a1_sd = est_sd[4], est_b1 = B1, est_b1_sd = est_sd[5],
                            est_c1_coef = C1_coef, est_c1_coef_sd = est_sd[5 + 1:nb1], est_d1 = D1, est_a2 = A2, est_a2_sd = est_sd[5 + nb1 + 1],
                            est_b2 = B2, est_b2_sd = est_sd[5 + nb1 + 2], est_c2_coef = C2_coef, est_c2_coef_sd = est_sd[5 + nb1 + 2 + 1:nb2], est_d2 = D2,
                            splines_degree = splines_degree, inner_knots1 = inner_knots1, inner_knots2 = inner_knots2)

          }else{

          }
        }
      }
    }else{
      if(a1_fix & b1_fix & a2_fix & b2_fix){
        if(scale_horizontal_fix & scale_vertical_fix){
          if(d1_fix & d2_fix){

            results <- list(convergence_code = fittedModel$code, iterations = fittedModel$iterations, loglikelihood_value = -fittedModel$minimum, theta = fittedModel$estimate, gradient = fittedModel$gradient,
                            est_beta = BETA, est_scale_horizontal = SCALE_HORIZONTAL, est_scale_vertical = SCALE_VERTICAL, est_a1 = A1, est_b1 = B1,
                            est_c1_coef = C1_coef, est_c1_coef_sd = est_sd[1:nb1], est_d1 = D1, est_a2 = A2, est_b2 = B2,
                            est_c2_coef = C2_coef, est_c2_coef_sd = est_sd[nb1 + 1:nb2], est_d2 = D2,
                            splines_degree = splines_degree, inner_knots1 = inner_knots1, inner_knots2 = inner_knots2)

          }else{

          }
        }else{
          if(d1_fix & d2_fix){

          }else{

          }
        }
      }else{
        if(scale_horizontal_fix & scale_vertical_fix){
          if(d1_fix & d2_fix){

          }else{

          }
        }else{
          if(d1_fix & d2_fix){

          }else{

          }
        }
      }
    }
  }
  return(results)
}

prepareResults_wls <- function(fittedModel,
                               BETA, SCALE_HORIZONTAL, SCALE_VERTICAL,
                               A1, B1, C1_coef, D1, A2, B2, C2_coef, D2,
                               est_sd,
                               beta_fix, scale_horizontal_fix,
                               scale_vertical_fix,
                               a1_fix, b1_fix, c1_fix, d1_fix,
                               a2_fix, b2_fix, c2_fix, d2_fix,
                               splines_degree,
                               inner_knots1, inner_knots2, nb1, nb2){

  if(c1_fix & c2_fix){
    if(d1_fix & d2_fix){

      results <- list(convergence_code = fittedModel$code, iterations = fittedModel$iterations, minimum = fittedModel$minimum, theta = fittedModel$estimate, gradient = fittedModel$gradient,
                      est_beta = BETA, est_beta_sd = est_sd[1], est_scale_horizontal = SCALE_HORIZONTAL, est_scale_horizontal_sd = est_sd[2],
                      est_scale_vertical = SCALE_VERTICAL, est_scale_vertical_sd = est_sd[3], est_a1 = A1, est_a1_sd = est_sd[4], est_b1 = B1, est_b1_sd = est_sd[5],
                      est_c1_coef = C1_coef, est_d1 = D1, est_a2 = A2, est_a2_sd = est_sd[6], est_b2 = B2, est_b2_sd = est_sd[7], est_c2_coef = C2_coef, est_d2 = D2,
                      splines_degree = splines_degree, inner_knots1 = inner_knots1, inner_knots2 = inner_knots2)
    }
  }else{
    if(!beta_fix){
      if(a1_fix & b1_fix & a2_fix & b2_fix){
        if(scale_horizontal_fix & scale_vertical_fix){
          if(d1_fix & d2_fix){

          }else{

          }
        }else{
          if(d1_fix & d2_fix){

            results <- list(convergence_code = fittedModel$code, iterations = fittedModel$iterations, minimum = fittedModel$minimum, theta = fittedModel$estimate, gradient = fittedModel$gradient,
                            est_beta = BETA, est_beta_sd = est_sd[1], est_scale_horizontal = SCALE_HORIZONTAL, est_scale_horizontal_sd = est_sd[2],
                            est_scale_vertical = SCALE_VERTICAL, est_scale_vertical_sd = est_sd[3], est_a1 = A1, est_b1 = B1,
                            est_c1_coef = C1_coef, est_c1_coef_sd = est_sd[3 + 1:nb1], est_d1 = D1, est_a2 = A2, est_b2 = B2,
                            est_c2_coef = C2_coef, est_c2_coef_sd = est_sd[3 + nb1 + 1:nb2], est_d2 = D2,
                            splines_degree = splines_degree, inner_knots1 = inner_knots1, inner_knots2 = inner_knots2)

          }else{

          }
        }
      }else{
        if(scale_horizontal_fix & scale_vertical_fix){
          if(d1_fix & d2_fix){

            results <- list(convergence_code = fittedModel$code, iterations = fittedModel$iterations, minimum = fittedModel$minimum, theta = fittedModel$estimate, gradient = fittedModel$gradient,
                            est_beta = BETA, est_beta_sd = est_sd[1], est_scale_horizontal = SCALE_HORIZONTAL, est_scale_vertical = SCALE_VERTICAL,
                            est_a1 = A1, est_a1_sd = est_sd[2], est_b1 = B1, est_b1_sd = est_sd[3],
                            est_c1_coef = C1_coef, est_c1_coef_sd = est_sd[3 + 1:nb1], est_d1 = D1, est_a2 = A2, est_a2_sd = est_sd[3 + nb1 + 1],
                            est_b2 = B2, est_b2_sd = est_sd[3 + nb1 + 2], est_c2_coef = C2_coef, est_c2_coef_sd = est_sd[3 + nb1 + 2 + 1:nb2], est_d2 = D2,
                            splines_degree = splines_degree, inner_knots1 = inner_knots1, inner_knots2 = inner_knots2)

          }else{

          }
        }else if(!scale_horizontal_fix & !scale_vertical_fix){
          if(d1_fix & d2_fix){

            results <- list(convergence_code = fittedModel$code, iterations = fittedModel$iterations, minimum = fittedModel$minimum, theta = fittedModel$estimate, gradient = fittedModel$gradient,
                            est_beta = BETA, est_beta_sd = est_sd[1], est_scale_horizontal = SCALE_HORIZONTAL, est_scale_horizontal_sd = est_sd[2],
                            est_scale_vertical = SCALE_VERTICAL, est_scale_vertical_sd = est_sd[3],
                            est_a1 = A1, est_a1_sd = est_sd[4], est_b1 = B1, est_b1_sd = est_sd[5],
                            est_c1_coef = C1_coef, est_c1_coef_sd = est_sd[5 + 1:nb1], est_d1 = D1, est_a2 = A2, est_a2_sd = est_sd[5 + nb1 + 1],
                            est_b2 = B2, est_b2_sd = est_sd[5 + nb1 + 2], est_c2_coef = C2_coef, est_c2_coef_sd = est_sd[5 + nb1 + 2 + 1:nb2], est_d2 = D2,
                            splines_degree = splines_degree, inner_knots1 = inner_knots1, inner_knots2 = inner_knots2)

          }else{

          }
        }
      }
    }else{
      if(a1_fix & b1_fix & a2_fix & b2_fix){
        if(scale_horizontal_fix & scale_vertical_fix){
          if(d1_fix & d2_fix){

            results <- list(convergence_code = fittedModel$code, iterations = fittedModel$iterations, minimum = fittedModel$minimum, theta = fittedModel$estimate, gradient = fittedModel$gradient,
                            est_beta = BETA, est_scale_horizontal = SCALE_HORIZONTAL, est_scale_vertical = SCALE_VERTICAL, est_a1 = A1, est_b1 = B1,
                            est_c1_coef = C1_coef, est_c1_coef_sd = est_sd[1:nb1], est_d1 = D1, est_a2 = A2, est_b2 = B2,
                            est_c2_coef = C2_coef, est_c2_coef_sd = est_sd[nb1 + 1:nb2], est_d2 = D2,
                            splines_degree = splines_degree, inner_knots1 = inner_knots1, inner_knots2 = inner_knots2)

          }else{

          }
        }else{
          if(d1_fix & d2_fix){

          }else{

          }
        }
      }else{
        if(scale_horizontal_fix & scale_vertical_fix){
          if(d1_fix & d2_fix){

          }else{

          }
        }else{
          if(d1_fix & d2_fix){

          }else{

          }
        }
      }
    }
  }
  return(results)
}

