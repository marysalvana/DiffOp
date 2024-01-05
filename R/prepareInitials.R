prepareInitials <- function(init_beta,
                            init_scale_horizontal, init_scale_vertical,
                            init_a1, init_b1, init_c1, init_d1,
                            init_a2, init_b2, init_c2, init_d2,
                            beta_fix, scale_horizontal_fix,
                            scale_vertical_fix,
                            a1_fix, b1_fix, c1_fix, d1_fix,
                            a2_fix, b2_fix, c2_fix, d2_fix){

  if(!beta_fix){
    if(scale_horizontal_fix & scale_vertical_fix){
      if(a1_fix & b1_fix & a2_fix & b2_fix){
        if(c1_fix & c2_fix){
          if(d1_fix & d2_fix){
            theta0 <- init_beta
          }else if(!d1_fix & !d2_fix){
            theta0 <- c(init_beta, init_d1, init_d2)
          }else{
            stop("d1 and d2 must either both be FIXED or ESTIMATED. Try again.")
          }
        }else if(!c1_fix & !c2_fix){
          if(d1_fix & d2_fix){
            theta0 <- c(init_beta, init_c1, init_c2)
          }else if(!d1_fix & !d2_fix){
            theta0 <- c(init_beta, init_c1, init_d1, init_c2, init_d2)
          }else{
            stop("d1 and d2 must either both be FIXED or ESTIMATED. Try again.")
          }
        }else{
          stop("c1 and c2 must either both be FIXED or ESTIMATED. Try again.")
        }
      }else if(!a1_fix & !b1_fix & !a2_fix & !b2_fix){
        if(c1_fix & c2_fix){
          if(d1_fix & d2_fix){
            theta0 <- c(init_beta, init_a1, init_b1, init_a2, init_b2)
          }else if(!d1_fix & !d2_fix){
            theta0 <- c(init_beta, init_a1, init_b1, init_d1, init_a2, init_b2, init_d2)
          }else{
            stop("d1 and d2 must either both be FIXED or ESTIMATED. Try again.")
          }
        }else if(!c1_fix & !c2_fix){
          if(d1_fix & d2_fix){
            theta0 <- c(init_beta, init_a1, init_b1, init_c1, init_a2, init_b2, init_c2)
          }else if(!d1_fix & !d2_fix){
            theta0 <- c(init_beta, init_a1, init_b1, init_c1, init_d1, init_a2, init_b2, init_c2, init_d2)
          }else{
            stop("d1 and d2 must either both be FIXED or ESTIMATED. Try again.")
          }
        }else{
          stop("c1 and c2 must either both be FIXED or ESTIMATED. Try again.")
        }
      }else{
        stop("a1, b1, a2, and b2 must either all be FIXED or ESTIMATED. Try again.")
      }
    }else if(!scale_horizontal_fix & !scale_vertical_fix){
      if(a1_fix & b1_fix & a2_fix & b2_fix){
        if(c1_fix & c2_fix){
          if(d1_fix & d2_fix){
            theta0 <- c(init_beta, init_scale_horizontal, init_scale_vertical)
          }else if(!d1_fix & !d2_fix){
            theta0 <- c(init_beta, init_scale_horizontal, init_scale_vertical, init_d1, init_d2)
          }else{
            stop("d1 and d2 must either both be FIXED or ESTIMATED. Try again.")
          }
        }else if(!c1_fix & !c2_fix){
          if(d1_fix & d2_fix){
            theta0 <- c(init_beta, init_scale_horizontal, init_scale_vertical, init_c1, init_c2)
          }else if(!d1_fix & !d2_fix){
            theta0 <- c(init_beta, init_scale_horizontal, init_scale_vertical, init_c1, init_d1, init_c2, init_d2)
          }else{
            stop("d1 and d2 must either both be FIXED or ESTIMATED. Try again.")
          }
        }else{
          stop("c1 and c2 must either both be FIXED or ESTIMATED. Try again.")
        }
      }else if(!a1_fix & !b1_fix & !a2_fix & !b2_fix){
        if(c1_fix & c2_fix){
          if(d1_fix & d2_fix){
            theta0 <- c(init_beta, init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_a2, init_b2)
          }else if(!d1_fix & !d2_fix){
            theta0 <- c(init_beta, init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_d1, init_a2, init_b2, init_d2)
          }else{
            stop("d1 and d2 must either both be FIXED or ESTIMATED. Try again.")
          }
        }else if(!c1_fix & !c2_fix){
          if(d1_fix & d2_fix){
            theta0 <- c(init_beta, init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_c1, init_a2, init_b2, init_c2)
          }else if(!d1_fix & !d2_fix){
            theta0 <- c(init_beta, init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_c1, init_d1, init_a2, init_b2, init_c2, init_d2)
          }else{
            stop("d1 and d2 must either both be FIXED or ESTIMATED. Try again.")
          }
        }else{
          stop("c1 and c2 must either both be FIXED or ESTIMATED. Try again.")
        }
      }else{
        stop("a1, b1, a2, and b2 must either all be FIXED or ESTIMATED. Try again.")
      }
    }else if(scale_horizontal_fix & !scale_vertical_fix){
      if(a1_fix & b1_fix & a2_fix & b2_fix){
        if(c1_fix & c2_fix){
          if(d1_fix & d2_fix){
            theta0 <- c(init_beta, init_scale_vertical)
          }else if(!d1_fix & !d2_fix){
            theta0 <- c(init_beta, init_scale_vertical, init_d1, init_d2)
          }else{
            stop("d1 and d2 must either both be FIXED or ESTIMATED. Try again.")
          }
        }else if(!c1_fix & !c2_fix){
          if(d1_fix & d2_fix){
            theta0 <- c(init_beta, init_scale_vertical, init_c1, init_c2)
          }else if(!d1_fix & !d2_fix){
            theta0 <- c(init_beta, init_scale_vertical, init_c1, init_d1, init_c2, init_d2)
          }else{
            stop("d1 and d2 must either both be FIXED or ESTIMATED. Try again.")
          }
        }else{
          stop("c1 and c2 must either both be FIXED or ESTIMATED. Try again.")
        }
      }else if(!a1_fix & !b1_fix & !a2_fix & !b2_fix){
        if(c1_fix & c2_fix){
          if(d1_fix & d2_fix){
            theta0 <- c(init_beta, init_scale_vertical, init_a1, init_b1, init_a2, init_b2)
          }else if(!d1_fix & !d2_fix){
            theta0 <- c(init_beta, init_scale_vertical, init_a1, init_b1, init_d1, init_a2, init_b2, init_d2)
          }else{
            stop("d1 and d2 must either both be FIXED or ESTIMATED. Try again.")
          }
        }else if(!c1_fix & !c2_fix){
          if(d1_fix & d2_fix){
            theta0 <- c(init_beta, init_scale_vertical, init_a1, init_b1, init_c1, init_a2, init_b2, init_c2)
          }else if(!d1_fix & !d2_fix){
            theta0 <- c(init_beta, init_scale_vertical, init_a1, init_b1, init_c1, init_d1, init_a2, init_b2, init_c2, init_d2)
          }else{
            stop("d1 and d2 must either both be FIXED or ESTIMATED. Try again.")
          }
        }else{
          stop("c1 and c2 must either both be FIXED or ESTIMATED. Try again.")
        }
      }else{
        stop("a1, b1, a2, and b2 must either all be FIXED or ESTIMATED. Try again.")
      }
    }else{
      stop("scale horizontal and scale vertical must either both be FIXED or ESTIMATED. Try again.")
    }
  }else if(beta_fix){
    if(scale_horizontal_fix & scale_vertical_fix){
      if(a1_fix & b1_fix & a2_fix & b2_fix){
        if(c1_fix & c2_fix){
          if(d1_fix & d2_fix){
            stop("All parameters are fixed. At least one must be estimated. Try again.")
          }else if(!d1_fix & !d2_fix){
            theta0 <- c(init_d1, init_d2)
          }else{
            stop("d1 and d2 must either both be FIXED or ESTIMATED. Try again.")
          }
        }else if(!c1_fix & !c2_fix){
          if(d1_fix & d2_fix){
            theta0 <- c(init_c1, init_c2)
          }else if(!d1_fix & !d2_fix){
            theta0 <- c(init_c1, init_d1, init_c2, init_d2)
          }else{
            stop("d1 and d2 must either both be FIXED or ESTIMATED. Try again.")
          }
        }else{
          stop("c1 and c2 must either both be FIXED or ESTIMATED. Try again.")
        }
      }else if(!a1_fix & !b1_fix & !a2_fix & !b2_fix){
        if(c1_fix & c2_fix){
          if(d1_fix & d2_fix){
            theta0 <- c(init_a1, init_b1, init_a2, init_b2)
          }else if(!d1_fix & !d2_fix){
            theta0 <- c(init_a1, init_b1, init_d1, init_a2, init_b2, init_d2)
          }else{
            stop("d1 and d2 must either both be FIXED or ESTIMATED. Try again.")
          }
        }else if(!c1_fix & !c2_fix){
          if(d1_fix & d2_fix){
            theta0 <- c(init_a1, init_b1, init_c1, init_a2, init_b2, init_c2)
          }else if(!d1_fix & !d2_fix){
            theta0 <- c(init_a1, init_b1, init_c1, init_d1, init_a2, init_b2, init_c2, init_d2)
          }else{
            stop("d1 and d2 must either both be FIXED or ESTIMATED. Try again.")
          }
        }else{
          stop("c1 and c2 must either both be FIXED or ESTIMATED. Try again.")
        }
      }else{
        stop("a1, b1, a2, and b2 must either all be FIXED or ESTIMATED. Try again.")
      }
    }else if(!scale_horizontal_fix & !scale_vertical_fix){
      if(a1_fix & b1_fix & a2_fix & b2_fix){
        if(c1_fix & c2_fix){
          if(d1_fix & d2_fix){
            theta0 <- c(init_scale_horizontal, init_scale_vertical)
          }else if(!d1_fix & !d2_fix){
            theta0 <- c(init_scale_horizontal, init_scale_vertical, init_d1, init_d2)
          }else{
            stop("d1 and d2 must either both be FIXED or ESTIMATED. Try again.")
          }
        }else if(!c1_fix & !c2_fix){
          if(d1_fix & d2_fix){
            theta0 <- c(init_scale_horizontal, init_scale_vertical, init_c1, init_c2)
          }else if(!d1_fix & !d2_fix){
            theta0 <- c(init_scale_horizontal, init_scale_vertical, init_c1, init_d1, init_c2, init_d2)
          }else{
            stop("d1 and d2 must either both be FIXED or ESTIMATED. Try again.")
          }
        }else{
          stop("c1 and c2 must either both be FIXED or ESTIMATED. Try again.")
        }
      }else if(!a1_fix & !b1_fix & !a2_fix & !b2_fix){
        if(c1_fix & c2_fix){
          if(d1_fix & d2_fix){
            theta0 <- c(init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_a2, init_b2)
          }else if(!d1_fix & !d2_fix){
            theta0 <- c(init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_d1, init_a2, init_b2, init_d2)
          }else{
            stop("d1 and d2 must either both be FIXED or ESTIMATED. Try again.")
          }
        }else if(!c1_fix & !c2_fix){
          if(d1_fix & d2_fix){
            theta0 <- c(init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_c1, init_a2, init_b2, init_c2)
          }else if(!d1_fix & !d2_fix){
            theta0 <- c(init_scale_horizontal, init_scale_vertical, init_a1, init_b1, init_c1, init_d1, init_a2, init_b2, init_c2, init_d2)
          }else{
            stop("d1 and d2 must either both be FIXED or ESTIMATED. Try again.")
          }
        }else{
          stop("c1 and c2 must either both be FIXED or ESTIMATED. Try again.")
        }
      }else{
        stop("a1, b1, a2, and b2 must either all be FIXED or ESTIMATED. Try again.")
      }
    }else if(scale_horizontal_fix & !scale_vertical_fix){
      if(a1_fix & b1_fix & a2_fix & b2_fix){
        if(c1_fix & c2_fix){
          if(d1_fix & d2_fix){
            theta0 <- c(init_scale_vertical)
          }else if(!d1_fix & !d2_fix){
            theta0 <- c(init_scale_vertical, init_d1, init_d2)
          }else{
            stop("d1 and d2 must either both be FIXED or ESTIMATED. Try again.")
          }
        }else if(!c1_fix & !c2_fix){
          if(d1_fix & d2_fix){
            theta0 <- c(init_scale_vertical, init_c1, init_c2)
          }else if(!d1_fix & !d2_fix){
            theta0 <- c(init_scale_vertical, init_c1, init_d1, init_c2, init_d2)
          }else{
            stop("d1 and d2 must either both be FIXED or ESTIMATED. Try again.")
          }
        }else{
          stop("c1 and c2 must either both be FIXED or ESTIMATED. Try again.")
        }
      }else if(!a1_fix & !b1_fix & !a2_fix & !b2_fix){
        if(c1_fix & c2_fix){
          if(d1_fix & d2_fix){
            theta0 <- c(init_scale_vertical, init_a1, init_b1, init_a2, init_b2)
          }else if(!d1_fix & !d2_fix){
            theta0 <- c(init_scale_vertical, init_a1, init_b1, init_d1, init_a2, init_b2, init_d2)
          }else{
            stop("d1 and d2 must either both be FIXED or ESTIMATED. Try again.")
          }
        }else if(!c1_fix & !c2_fix){
          if(d1_fix & d2_fix){
            theta0 <- c(init_scale_vertical, init_a1, init_b1, init_c1, init_a2, init_b2, init_c2)
          }else if(!d1_fix & !d2_fix){
            theta0 <- c(init_scale_vertical, init_a1, init_b1, init_c1, init_d1, init_a2, init_b2, init_c2, init_d2)
          }else{
            stop("d1 and d2 must either both be FIXED or ESTIMATED. Try again.")
          }
        }else{
          stop("c1 and c2 must either both be FIXED or ESTIMATED. Try again.")
        }
      }else{
        stop("a1, b1, a2, and b2 must either all be FIXED or ESTIMATED. Try again.")
      }
    }else{
      stop("scale horizontal and scale vertical must either both be FIXED or ESTIMATED. Try again.")
    }
  }
  return(theta0)
}

