computeSD <- function(fittedModel,
                      beta_fix, scale_horizontal_fix,
                      scale_vertical_fix,
                      a1_fix, b1_fix, c1_fix, d1_fix,
                      a2_fix, b2_fix, c2_fix, d2_fix,
                      beta_scaling, horizontal_scale_scaling,
                      vertical_scale_scaling, a1_scaling, b1_scaling,
                      c1_coef_scaling, d1_scaling,
                      a2_scaling, b2_scaling,
                      c2_coef_scaling, d2_scaling, nb1, nb2){

  est_param <- fittedModel$estimate
  est_gradient <- fittedModel$gradient
  est_hessian <- fittedModel$hessian
  fisher_info<-solve(est_hessian)

  if(!beta_fix){

    if(scale_horizontal_fix & scale_vertical_fix){

      if(a1_fix & b1_fix & a2_fix & b2_fix){

        if(c1_fix & c2_fix){

          if(d1_fix & d2_fix){

            j <- diag(c(beta_scaling * exp(-est_param[1] * beta_scaling) / (1 + exp(-est_param[1] * beta_scaling))^2),
                      nrow = length(est_param), ncol = length(est_param))

          }else if(!d1_fix & !d2_fix){

            j <- diag(c(beta_scaling * exp(-est_param[1] * beta_scaling) / (1 + exp(-est_param[1] * beta_scaling))^2,
                        d1_scaling,
                        d2_scaling),
                      nrow = length(est_param), ncol = length(est_param))

          }

        }else if(!c1_fix & !c2_fix){

          if(d1_fix & d2_fix){

            j <- diag(c(beta_scaling * exp(-est_param[1] * beta_scaling) / (1 + exp(-est_param[1] * beta_scaling))^2,
                        rep(c1_coef_scaling, nb1),
                        rep(c2_coef_scaling, nb2)),
                      nrow = length(est_param), ncol = length(est_param))

          }else if(!d1_fix & !d2_fix){

            j <- diag(c(beta_scaling * exp(-est_param[1] * beta_scaling) / (1 + exp(-est_param[1] * beta_scaling))^2,
                        rep(c1_coef_scaling, nb1),
                        d1_scaling,
                        rep(c2_coef_scaling, nb2),
                        d2_scaling),
                      nrow = length(est_param), ncol = length(est_param))

          }

        }

      }else if(!a1_fix & !b1_fix & !a2_fix & !b2_fix){

        if(c1_fix & c2_fix){

          if(d1_fix & d2_fix){

            j <- diag(c(beta_scaling * exp(-est_param[1] * beta_scaling) / (1 + exp(-est_param[1] * beta_scaling))^2,
                        exp(est_param[2]) * a1_scaling,
                        b1_scaling,
                        a2_scaling,
                        b2_scaling),
                      nrow = length(est_param), ncol = length(est_param))

          }else if(!d1_fix & !d2_fix){

            j <- diag(c(beta_scaling * exp(-est_param[1] * beta_scaling) / (1 + exp(-est_param[1] * beta_scaling))^2,
                        exp(est_param[2]) * a1_scaling,
                        b1_scaling,
                        d1_scaling,
                        a2_scaling,
                        b2_scaling,
                        d2_scaling),
                      nrow = length(est_param), ncol = length(est_param))

          }

        }else if(!c1_fix & !c2_fix){

          if(d1_fix & d2_fix){

            j <- diag(c(beta_scaling * exp(-est_param[1] * beta_scaling) / (1 + exp(-est_param[1] * beta_scaling))^2,
                        exp(est_param[2]) * a1_scaling,
                        b1_scaling,
                        rep(c1_coef_scaling, nb1),
                        a2_scaling,
                        b2_scaling,
                        rep(c2_coef_scaling, nb2)),
                      nrow = length(est_param), ncol = length(est_param))

          }else if(!d1_fix & !d2_fix){

            j <- diag(c(beta_scaling * exp(-est_param[1] * beta_scaling) / (1 + exp(-est_param[1] * beta_scaling))^2,
                        exp(est_param[2]) * a1_scaling,
                        b1_scaling,
                        rep(c1_coef_scaling, nb1),
                        d1_scaling,
                        a2_scaling,
                        b2_scaling,
                        rep(c2_coef_scaling, nb2),
                        d2_scaling),
                      nrow = length(est_param), ncol = length(est_param))

          }

        }

      }

    }else if(!scale_horizontal_fix & !scale_vertical_fix){

      if(a1_fix & b1_fix & a2_fix & b2_fix){

        if(c1_fix & c2_fix){

          if(d1_fix & d2_fix){

            j <- diag(c(beta_scaling * exp(-est_param[1] * beta_scaling) / (1 + exp(-est_param[1] * beta_scaling))^2,
                        exp(est_param[2]) * horizontal_scale_scaling,
                        exp(est_param[3]) * vertical_scale_scaling),
                      nrow = length(est_param), ncol = length(est_param))

          }else if(!d1_fix & !d2_fix){

            j <- diag(c(beta_scaling * exp(-est_param[1] * beta_scaling) / (1 + exp(-est_param[1] * beta_scaling))^2,
                        exp(est_param[2]) * horizontal_scale_scaling,
                        exp(est_param[3]) * vertical_scale_scaling,
                        d1_scaling,
                        d2_scaling),
                      nrow = length(est_param), ncol = length(est_param))

          }

        }else if(!c1_fix & !c2_fix){

          if(d1_fix & d2_fix){

            j <- diag(c(beta_scaling * exp(-est_param[1] * beta_scaling) / (1 + exp(-est_param[1] * beta_scaling))^2,
                        exp(est_param[2]) * horizontal_scale_scaling,
                        exp(est_param[3]) * vertical_scale_scaling,
                        rep(c1_coef_scaling, nb1),
                        rep(c2_coef_scaling, nb2)),
                      nrow = length(est_param), ncol = length(est_param))

          }else if(!d1_fix & !d2_fix){

            j <- diag(c(beta_scaling * exp(-est_param[1] * beta_scaling) / (1 + exp(-est_param[1] * beta_scaling))^2,
                        exp(est_param[2]) * horizontal_scale_scaling,
                        exp(est_param[3]) * vertical_scale_scaling,
                        rep(c1_coef_scaling, nb1),
                        d1_scaling,
                        rep(c2_coef_scaling, nb2),
                        d2_scaling),
                      nrow = length(est_param), ncol = length(est_param))

          }

        }

      }else if(!a1_fix & !b1_fix & !a2_fix & !b2_fix){
        if(c1_fix & c2_fix){
          if(d1_fix & d2_fix){

            j <- diag(c(beta_scaling * exp(-est_param[1] * beta_scaling) / (1 + exp(-est_param[1] * beta_scaling))^2,
                        exp(est_param[2]) * horizontal_scale_scaling,
                        exp(est_param[3]) * vertical_scale_scaling,
                        exp(est_param[4]) * a1_scaling,
                        b1_scaling,
                        a2_scaling,
                        b2_scaling),
                      nrow = length(est_param), ncol = length(est_param))

          }else if(!d1_fix & !d2_fix){

            j <- diag(c(beta_scaling * exp(-est_param[1] * beta_scaling) / (1 + exp(-est_param[1] * beta_scaling))^2,
                        exp(est_param[2]) * horizontal_scale_scaling,
                        exp(est_param[3]) * vertical_scale_scaling,
                        exp(est_param[4]) * a1_scaling,
                        b1_scaling,
                        d1_scaling,
                        a2_scaling,
                        b2_scaling,
                        d2_scaling),
                      nrow = length(est_param), ncol = length(est_param))

          }

        }else if(!c1_fix & !c2_fix){

          if(d1_fix & d2_fix){

            j <- diag(c(beta_scaling * exp(-est_param[1] * beta_scaling) / (1 + exp(-est_param[1] * beta_scaling))^2,
                        exp(est_param[2]) * horizontal_scale_scaling,
                        exp(est_param[3]) * vertical_scale_scaling,
                        exp(est_param[4]) * a1_scaling,
                        b1_scaling,
                        rep(c1_coef_scaling, nb1),
                        a2_scaling,
                        b2_scaling,
                        rep(c2_coef_scaling, nb2)),
                      nrow = length(est_param), ncol = length(est_param))

          }else if(!d1_fix & !d2_fix){

            j <- diag(c(beta_scaling * exp(-est_param[1] * beta_scaling) / (1 + exp(-est_param[1] * beta_scaling))^2,
                        exp(est_param[2]) * horizontal_scale_scaling,
                        exp(est_param[3]) * vertical_scale_scaling,
                        exp(est_param[4]) * a1_scaling,
                        b1_scaling,
                        rep(c1_coef_scaling, nb1),
                        d1_scaling,
                        a2_scaling,
                        b2_scaling,
                        rep(c2_coef_scaling, nb2),
                        d2_scaling),
                      nrow = length(est_param), ncol = length(est_param))

          }

        }

      }

    }

  }else if(beta_fix){

    if(scale_horizontal_fix & scale_vertical_fix){

      if(a1_fix & b1_fix & a2_fix & b2_fix){

        if(c1_fix & c2_fix){

          if(!d1_fix & !d2_fix){

            j <- diag(c(d1_scaling,
                        d2_scaling),
                      nrow = length(est_param), ncol = length(est_param))

          }

        }else if(!c1_fix & !c2_fix){

          if(d1_fix & d2_fix){

            j <- diag(c(rep(c1_coef_scaling, nb1),
                        rep(c2_coef_scaling, nb2)),
                      nrow = length(est_param), ncol = length(est_param))

          }else if(!d1_fix & !d2_fix){

            j <- diag(c(rep(c1_coef_scaling, nb1),
                        d1_scaling,
                        rep(c2_coef_scaling, nb2),
                        d2_scaling),
                      nrow = length(est_param), ncol = length(est_param))

          }

        }

      }else if(!a1_fix & !b1_fix & !a2_fix & !b2_fix){

        if(c1_fix & c2_fix){

          if(d1_fix & d2_fix){

            j <- diag(c(exp(est_param[1]) * a1_scaling,
                        b1_scaling,
                        a2_scaling,
                        b2_scaling),
                      nrow = length(est_param), ncol = length(est_param))

          }else if(!d1_fix & !d2_fix){

            j <- diag(c(exp(est_param[1]) * a1_scaling,
                        b1_scaling,
                        d1_scaling,
                        a2_scaling,
                        b2_scaling,
                        d2_scaling),
                      nrow = length(est_param), ncol = length(est_param))

          }

        }else if(!c1_fix & !c2_fix){

          if(d1_fix & d2_fix){

            j <- diag(c(exp(est_param[1]) * a1_scaling,
                        b1_scaling,
                        rep(c1_coef_scaling, nb1),
                        a2_scaling,
                        b2_scaling,
                        rep(c2_coef_scaling, nb2)),
                      nrow = length(est_param), ncol = length(est_param))

          }else if(!d1_fix & !d2_fix){

            j <- diag(c(exp(est_param[1]) * a1_scaling,
                        b1_scaling,
                        rep(c1_coef_scaling, nb1),
                        d1_scaling,
                        a2_scaling,
                        b2_scaling,
                        rep(c2_coef_scaling, nb2),
                        d2_scaling),
                      nrow = length(est_param), ncol = length(est_param))

          }

        }

      }

    }else if(!scale_horizontal_fix & !scale_vertical_fix){

      if(a1_fix & b1_fix & a2_fix & b2_fix){

        if(c1_fix & c2_fix){

          if(d1_fix & d2_fix){

            j <- diag(c(exp(est_param[1]) * horizontal_scale_scaling,
                        exp(est_param[2]) * vertical_scale_scaling),
                      nrow = length(est_param), ncol = length(est_param))

          }else if(!d1_fix & !d2_fix){

            j <- diag(c(exp(est_param[1]) * horizontal_scale_scaling,
                        exp(est_param[2]) * vertical_scale_scaling,
                        d1_scaling,
                        d2_scaling),
                      nrow = length(est_param), ncol = length(est_param))

          }

        }else if(!c1_fix & !c2_fix){

          if(d1_fix & d2_fix){

            j <- diag(c(exp(est_param[1]) * horizontal_scale_scaling,
                        exp(est_param[2]) * vertical_scale_scaling,
                        rep(c1_coef_scaling, nb1),
                        rep(c2_coef_scaling, nb2)),
                      nrow = length(est_param), ncol = length(est_param))

          }else if(!d1_fix & !d2_fix){

            j <- diag(c(exp(est_param[1]) * horizontal_scale_scaling,
                        exp(est_param[2]) * vertical_scale_scaling,
                        rep(c1_coef_scaling, nb1),
                        d1_scaling,
                        rep(c2_coef_scaling, nb2),
                        d2_scaling),
                      nrow = length(est_param), ncol = length(est_param))

          }

        }

      }else if(!a1_fix & !b1_fix & !a2_fix & !b2_fix){
        if(c1_fix & c2_fix){
          if(d1_fix & d2_fix){

            j <- diag(c(exp(est_param[1]) * horizontal_scale_scaling,
                        exp(est_param[2]) * vertical_scale_scaling,
                        exp(est_param[3]) * a1_scaling,
                        b1_scaling,
                        a2_scaling,
                        b2_scaling),
                      nrow = length(est_param), ncol = length(est_param))

          }else if(!d1_fix & !d2_fix){

            j <- diag(c(exp(est_param[1]) * horizontal_scale_scaling,
                        exp(est_param[2]) * vertical_scale_scaling,
                        exp(est_param[3]) * a1_scaling,
                        b1_scaling,
                        d1_scaling,
                        a2_scaling,
                        b2_scaling,
                        d2_scaling),
                      nrow = length(est_param), ncol = length(est_param))

          }

        }else if(!c1_fix & !c2_fix){

          if(d1_fix & d2_fix){

            j <- diag(c(exp(est_param[1]) * horizontal_scale_scaling,
                        exp(est_param[2]) * vertical_scale_scaling,
                        exp(est_param[3]) * a1_scaling,
                        b1_scaling,
                        rep(c1_coef_scaling, nb1),
                        a2_scaling,
                        b2_scaling,
                        rep(c2_coef_scaling, nb2)),
                      nrow = length(est_param), ncol = length(est_param))

          }else if(!d1_fix & !d2_fix){

            j <- diag(c(exp(est_param[1]) * horizontal_scale_scaling,
                        exp(est_param[2]) * vertical_scale_scaling,
                        exp(est_param[3]) * a1_scaling,
                        b1_scaling,
                        rep(c1_coef_scaling, nb1),
                        d1_scaling,
                        a2_scaling,
                        b2_scaling,
                        rep(c2_coef_scaling, nb2),
                        d2_scaling),
                      nrow = length(est_param), ncol = length(est_param))

          }

        }

      }

    }

  }

  variance <- j %*% fisher_info %*% t(j)
  est_sd <- sqrt(diag(variance))

  return(est_sd)

}
