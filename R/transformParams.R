transformParams <- function(theta, init_beta,
                            init_scale_horizontal, init_scale_vertical,
                            init_a1, init_b1, init_c1, init_d1,
                            init_a2, init_b2, init_c2, init_d2,
                            beta_fix, scale_horizontal_fix,
                            scale_vertical_fix,
                            a1_fix, b1_fix, c1_fix, d1_fix,
                            a2_fix, b2_fix, c2_fix, d2_fix,
                            beta_scaling, horizontal_scale_scaling,
                            vertical_scale_scaling, a1_scaling, b1_scaling,
                            c1_scaling, d1_scaling,
                            a2_scaling, b2_scaling,
                            c2_scaling, d2_scaling,
                            splines_degree, nb1, nb2){

  if(!beta_fix){

    #BETA <- 2 / (1 + exp(-theta[1] * beta_scaling)) - 1
    BETA <- 1 / (1 + exp(-theta[1] * beta_scaling))

    if(scale_horizontal_fix & scale_vertical_fix){

      SCALE_HORIZONTAL <- init_scale_horizontal
      SCALE_VERTICAL <- init_scale_vertical

      if(a1_fix & b1_fix & a2_fix & b2_fix){

        A1 <- init_a1
        B1 <- init_b1
        A2 <- init_a2
        B2 <- init_b2

        if(c1_fix & c2_fix){

          C1 <- init_c1
          C2 <- init_c2

          if(d1_fix & d2_fix){

            D1 <- init_d1
            D2 <- init_d2

          }else if(!d1_fix & !d2_fix){

            D1 <- theta[2] * d1_scaling
            D2 <- theta[3] * d2_scaling

          }

        }else if(!c1_fix & !c2_fix){

          if(d1_fix & d2_fix){

            C1 <- theta[1 + 1:nb1] * c1_scaling
            D1 <- init_d1

            C2 <- theta[1 + nb1 + 1:nb2] * c2_scaling
            D2 <- init_d2

          }else if(!d1_fix & !d2_fix){

            C1 <- theta[1 + 1:nb1] * c1_scaling
            D1 <- theta[1 + nb1 + 1] * d1_scaling

            C2 <- theta[1 + nb1 + 1 + 1:nb2] * c2_scaling
            D2 <- theta[1 + nb1 + 1 + nb2 + 1] * d2_scaling

          }

        }

      }else if(!a1_fix & !b1_fix & !a2_fix & !b2_fix){

        if(c1_fix & c2_fix){

          C1 <- init_c1
          C2 <- init_c2

          if(d1_fix & d2_fix){

            A1 <- exp(theta[2]) * a1_scaling
            B1 <- theta[3] * b1_scaling
            D1 <- init_d1

            A2 <- theta[4] * a2_scaling
            B2 <- theta[5] * b2_scaling
            D2 <- init_d2

          }else if(!d1_fix & !d2_fix){

            A1 <- exp(theta[2]) * a1_scaling
            B1 <- theta[3] * b1_scaling
            D1 <- theta[4] * d1_scaling

            A2 <- theta[5] * a2_scaling
            B2 <- theta[6] * b2_scaling
            D2 <- theta[7] * d2_scaling

          }

        }else if(!c1_fix & !c2_fix){

          if(d1_fix & d2_fix){

            A1 <- exp(theta[2]) * a1_scaling
            B1 <- theta[3] * b1_scaling
            C1 <- theta[3 + 1:nb1] * c1_scaling
            D1 <- init_d1

            A2 <- theta[3 + nb1 + 1] * a2_scaling
            B2 <- theta[3 + nb1 + 2] * b2_scaling
            C2 <- theta[3 + nb1 + 2 + 1:nb2] * c2_scaling
            D2 <- init_d2

          }else if(!d1_fix & !d2_fix){

            A1 <- exp(theta[2]) * a1_scaling
            B1 <- theta[3] * b1_scaling
            C1 <- theta[3 + 1:nb1] * c1_scaling
            D1 <- theta[3 + nb1 + 1] * d1_scaling

            A2 <- theta[3 + nb1 + 2] * a2_scaling
            B2 <- theta[3 + nb1 + 3] * b2_scaling
            C2 <- theta[3 + nb1 + 3 + 1:nb2] * c2_scaling
            D2 <- theta[3 + nb1 + 3 + nb2 + 1] * d1_scaling

          }

        }

      }

    }else if(!scale_horizontal_fix & !scale_vertical_fix){

      SCALE_HORIZONTAL <- exp(theta[2]) * horizontal_scale_scaling
      SCALE_VERTICAL <- exp(theta[3]) * vertical_scale_scaling

      if(a1_fix & b1_fix & a2_fix & b2_fix){

        A1 <- init_a1
        B1 <- init_b1
        A2 <- init_a2
        B2 <- init_b2

        if(c1_fix & c2_fix){

          C1 <- init_c1
          C2 <- init_c2

          if(d1_fix & d2_fix){

            D1 <- init_d1
            D2 <- init_d2

          }else if(!d1_fix & !d2_fix){

            D1 <- theta[4] * d1_scaling
            D2 <- theta[5] * d2_scaling

          }

        }else if(!c1_fix & !c2_fix){

          if(d1_fix & d2_fix){

            C1 <- theta[3 + 1:nb1] * c1_scaling
            D1 <- init_d1

            C2 <- theta[3 + nb1 + 1:nb2] * c2_scaling
            D2 <- init_d2

          }else if(!d1_fix & !d2_fix){

            C1 <- theta[3 + 1:nb1] * c1_scaling
            D1 <- theta[3 + nb1 + 1] * d1_scaling

            C2 <- theta[3 + nb1 + 1 + 1:nb2] * c2_scaling
            D2 <- theta[3 + nb1 + 1 + nb2 + 1] * d2_scaling

          }

        }

      }else if(!a1_fix & !b1_fix & !a2_fix & !b2_fix){

        if(c1_fix & c2_fix){

          C1 <- init_c1
          C2 <- init_c2

          if(d1_fix & d2_fix){

            A1 <- exp(theta[4]) * a1_scaling
            B1 <- theta[5] * b1_scaling
            D1 <- init_d1

            A2 <- theta[6] * a2_scaling
            B2 <- theta[7] * b2_scaling
            D2 <- init_d2

          }else if(!d1_fix & !d2_fix){

            A1 <- exp(theta[4]) * a1_scaling
            B1 <- theta[5] * b1_scaling
            D1 <- theta[6] * d1_scaling

            A2 <- theta[7] * a2_scaling
            B2 <- theta[8] * b2_scaling
            D2 <- theta[9] * d2_scaling

          }

        }else if(!c1_fix & !c2_fix){

          if(d1_fix & d2_fix){

            A1 <- exp(theta[4]) * a1_scaling
            B1 <- theta[5] * b1_scaling
            C1 <- theta[5 + 1:nb1] * c1_scaling
            D1 <- init_d1

            A2 <- theta[5 + nb1 + 1] * a2_scaling
            B2 <- theta[5 + nb1 + 2] * b2_scaling
            C2 <- theta[5 + nb1 + 2 + 1:nb2] * c2_scaling
            D2 <- init_d2

          }else if(!d1_fix & !d2_fix){

            A1 <- exp(theta[4]) * a1_scaling
            B1 <- theta[5] * b1_scaling
            C1 <- theta[5 + 1:nb1] * c1_scaling
            D1 <- theta[5 + nb1 + 1] * d1_scaling

            A2 <- theta[5 + nb1 + 2] * a2_scaling
            B2 <- theta[5 + nb1 + 3] * b2_scaling
            C2 <- theta[5 + nb1 + 3 + 1:nb2] * c2_scaling
            D2 <- theta[5 + nb1 + 3 + nb2 + 1] * d1_scaling

          }

        }

      }

    }else if(scale_horizontal_fix & !scale_vertical_fix){

      SCALE_HORIZONTAL <- init_scale_horizontal
      SCALE_VERTICAL <- exp(theta[2]) * vertical_scale_scaling

      if(a1_fix & b1_fix & a2_fix & b2_fix){

        A1 <- init_a1
        B1 <- init_b1
        A2 <- init_a2
        B2 <- init_b2

        if(c1_fix & c2_fix){

          C1 <- init_c1
          C2 <- init_c2

          if(d1_fix & d2_fix){

            D1 <- init_d1
            D2 <- init_d2

          }else if(!d1_fix & !d2_fix){

            D1 <- theta[3] * d1_scaling
            D2 <- theta[4] * d2_scaling

          }

        }else if(!c1_fix & !c2_fix){

          if(d1_fix & d2_fix){

            C1 <- theta[2 + 1:nb1] * c1_scaling
            D1 <- init_d1

            C2 <- theta[2 + nb1 + 1:nb2] * c2_scaling
            D2 <- init_d2

          }else if(!d1_fix & !d2_fix){

            C1 <- theta[2 + 1:nb1] * c1_scaling
            D1 <- theta[2 + nb1 + 1] * d1_scaling

            C2 <- theta[2 + nb1 + 1 + 1:nb2] * c2_scaling
            D2 <- theta[2 + nb1 + 1 + nb2 + 1] * d2_scaling

          }

        }

      }else if(!a1_fix & !b1_fix & !a2_fix & !b2_fix){

        if(c1_fix & c2_fix){

          C1 <- init_c1
          C2 <- init_c2

          if(d1_fix & d2_fix){

            A1 <- exp(theta[3]) * a1_scaling
            B1 <- theta[4] * b1_scaling
            D1 <- init_d1

            A2 <- theta[5] * a2_scaling
            B2 <- theta[6] * b2_scaling
            D2 <- init_d2

          }else if(!d1_fix & !d2_fix){

            A1 <- exp(theta[3]) * a1_scaling
            B1 <- theta[4] * b1_scaling
            D1 <- theta[5] * d1_scaling

            A2 <- theta[6] * a2_scaling
            B2 <- theta[7] * b2_scaling
            D2 <- theta[8] * d2_scaling

          }

        }else if(!c1_fix & !c2_fix){

          if(d1_fix & d2_fix){

            A1 <- exp(theta[3]) * a1_scaling
            B1 <- theta[4] * b1_scaling
            C1 <- theta[4 + 1:nb1] * c1_scaling
            D1 <- init_d1

            A2 <- theta[4 + nb1 + 1] * a2_scaling
            B2 <- theta[4 + nb1 + 2] * b2_scaling
            C2 <- theta[4 + nb1 + 2 + 1:nb2] * c2_scaling
            D2 <- init_d2

          }else if(!d1_fix & !d2_fix){

            A1 <- exp(theta[3]) * a1_scaling
            B1 <- theta[4] * b1_scaling
            C1 <- theta[4 + 1:nb1] * c1_scaling
            D1 <- theta[4 + nb1 + 1] * d1_scaling

            A2 <- theta[4 + nb1 + 2] * a2_scaling
            B2 <- theta[4 + nb1 + 3] * b2_scaling
            C2 <- theta[4 + nb1 + 3 + 1:nb2] * c2_scaling
            D2 <- theta[4 + nb1 + 3 + nb2 + 1] * d1_scaling
          }
        }
      }
    }
  }else if(beta_fix){

    BETA <- init_beta

    if(scale_horizontal_fix & scale_vertical_fix){

      SCALE_HORIZONTAL <- init_scale_horizontal
      SCALE_VERTICAL <- init_scale_vertical

      if(a1_fix & b1_fix & a2_fix & b2_fix){

        A1 <- init_a1
        B1 <- init_b1
        A2 <- init_a2
        B2 <- init_b2

        if(c1_fix & c2_fix){

          C1 <- init_c1
          C2 <- init_c2

          if(!d1_fix & !d2_fix){

            D1 <- theta[1] * d1_scaling
            D2 <- theta[2] * d2_scaling

          }

        }else if(!c1_fix & !c2_fix){

          if(d1_fix & d2_fix){

            C1 <- theta[1:nb1] * c1_scaling
            D1 <- init_d1

            C2 <- theta[nb1 + 1:nb2] * c2_scaling
            D2 <- init_d2

          }else if(!d1_fix & !d2_fix){

            C1 <- theta[1:nb1] * c1_scaling
            D1 <- theta[nb1 + 1] * d1_scaling

            C2 <- theta[nb1 + 1 + 1:nb2] * c2_scaling
            D2 <- theta[nb1 + 1 + nb2 + 1] * d2_scaling

          }

        }

      }else if(!a1_fix & !b1_fix & !a2_fix & !b2_fix){

        if(c1_fix & c2_fix){

          C1 <- init_c1
          C2 <- init_c2

          if(d1_fix & d2_fix){

            A1 <- exp(theta[1]) * a1_scaling
            B1 <- theta[2] * b1_scaling
            D1 <- init_d1

            A2 <- theta[3] * a2_scaling
            B2 <- theta[4] * b2_scaling
            D2 <- init_d2

          }else if(!d1_fix & !d2_fix){

            A1 <- exp(theta[1]) * a1_scaling
            B1 <- theta[2] * b1_scaling
            D1 <- theta[3] * d1_scaling

            A2 <- theta[4] * a2_scaling
            B2 <- theta[5] * b2_scaling
            D2 <- theta[6] * d2_scaling

          }

        }else if(!c1_fix & !c2_fix){

          if(d1_fix & d2_fix){

            A1 <- exp(theta[1]) * a1_scaling
            B1 <- theta[2] * b1_scaling
            C1 <- theta[2 + 1:nb1] * c1_scaling
            D1 <- init_d1

            A2 <- theta[2 + nb1 + 1] * a2_scaling
            B2 <- theta[2 + nb1 + 2] * b2_scaling
            C2 <- theta[2 + nb1 + 2 + 1:nb2] * c2_scaling
            D2 <- init_d2

          }else if(!d1_fix & !d2_fix){

            A1 <- exp(theta[1]) * a1_scaling
            B1 <- theta[2] * b1_scaling
            C1 <- theta[2 + 1:nb1] * c1_scaling
            D1 <- theta[2 + nb1 + 1] * d1_scaling

            A2 <- theta[2 + nb1 + 2] * a2_scaling
            B2 <- theta[2 + nb1 + 3] * b2_scaling
            C2 <- theta[2 + nb1 + 3 + 1:nb2] * c2_scaling
            D2 <- theta[2 + nb1 + 3 + nb2 + 1] * d1_scaling
          }
        }
      }
    }else if(!scale_horizontal_fix & !scale_vertical_fix){

      SCALE_HORIZONTAL <- exp(theta[1]) * horizontal_scale_scaling
      SCALE_VERTICAL <- exp(theta[2]) * vertical_scale_scaling

      if(a1_fix & b1_fix & a2_fix & b2_fix){

        A1 <- init_a1
        B1 <- init_b1
        A2 <- init_a2
        B2 <- init_b2

        if(c1_fix & c2_fix){

          C1 <- init_c1
          C2 <- init_c2

          if(d1_fix & d2_fix){

            D1 <- init_d1
            D2 <- init_d2

          }else if(!d1_fix & !d2_fix){

            D1 <- theta[3] * d1_scaling
            D2 <- theta[4] * d2_scaling

          }

        }else if(!c1_fix & !c2_fix){

          if(d1_fix & d2_fix){

            C1 <- theta[2 + 1:nb1] * c1_scaling
            D1 <- init_d1

            C2 <- theta[2 + nb1 + 1:nb2] * c2_scaling
            D2 <- init_d2

          }else if(!d1_fix & !d2_fix){

            C1 <- theta[2 + 1:nb1] * c1_scaling
            D1 <- theta[2 + nb1 + 1] * d1_scaling

            C2 <- theta[2 + nb1 + 1 + 1:nb2] * c2_scaling
            D2 <- theta[2 + nb1 + 1 + nb2 + 1] * d2_scaling

          }

        }

      }else if(!a1_fix & !b1_fix & !a2_fix & !b2_fix){

        if(c1_fix & c2_fix){

          C1 <- init_c1
          C2 <- init_c2

          if(d1_fix & d2_fix){

            A1 <- exp(theta[3]) * a1_scaling
            B1 <- theta[4] * b1_scaling
            D1 <- init_d1

            A2 <- theta[5] * a2_scaling
            B2 <- theta[6] * b2_scaling
            D2 <- init_d2

          }else if(!d1_fix & !d2_fix){

            A1 <- exp(theta[3]) * a1_scaling
            B1 <- theta[4] * b1_scaling
            D1 <- theta[5] * d1_scaling

            A2 <- theta[6] * a2_scaling
            B2 <- theta[7] * b2_scaling
            D2 <- theta[8] * d2_scaling

          }

        }else if(!c1_fix & !c2_fix){

          if(d1_fix & d2_fix){

            A1 <- exp(theta[3]) * a1_scaling
            B1 <- theta[4] * b1_scaling
            C1 <- theta[4 + 1:nb1] * c1_scaling
            D1 <- init_d1

            A2 <- theta[4 + nb1 + 1] * a2_scaling
            B2 <- theta[4 + nb1 + 2] * b2_scaling
            C2 <- theta[4 + nb1 + 2 + 1:nb2] * c2_scaling
            D2 <- init_d2

          }else if(!d1_fix & !d2_fix){

            A1 <- exp(theta[3]) * a1_scaling
            B1 <- theta[4] * b1_scaling
            C1 <- theta[4 + 1:nb1] * c1_scaling
            D1 <- theta[4 + nb1 + 1] * d1_scaling

            A2 <- theta[4 + nb1 + 2] * a2_scaling
            B2 <- theta[4 + nb1 + 3] * b2_scaling
            C2 <- theta[4 + nb1 + 3 + 1:nb2] * c2_scaling
            D2 <- theta[4 + nb1 + 3 + nb2 + 1] * d1_scaling
          }
        }
      }
    }else if(scale_horizontal_fix & !scale_vertical_fix){

      SCALE_HORIZONTAL <- init_scale_horizontal
      SCALE_VERTICAL <- exp(theta[1]) * vertical_scale_scaling

      if(a1_fix & b1_fix & a2_fix & b2_fix){

        A1 <- init_a1
        B1 <- init_b1
        A2 <- init_a2
        B2 <- init_b2

        if(c1_fix & c2_fix){

          C1 <- init_c1
          C2 <- init_c2

          if(d1_fix & d2_fix){

            D1 <- init_d1
            D2 <- init_d2

          }else if(!d1_fix & !d2_fix){

            D1 <- theta[2] * d1_scaling
            D2 <- theta[3] * d2_scaling

          }

        }else if(!c1_fix & !c2_fix){

          if(d1_fix & d2_fix){

            C1 <- theta[1 + 1:nb1] * c1_scaling
            D1 <- init_d1

            C2 <- theta[1 + nb1 + 1:nb2] * c2_scaling
            D2 <- init_d2

          }else if(!d1_fix & !d2_fix){

            C1 <- theta[1 + 1:nb1] * c1_scaling
            D1 <- theta[1 + nb1 + 1] * d1_scaling

            C2 <- theta[1 + nb1 + 1 + 1:nb2] * c2_scaling
            D2 <- theta[1 + nb1 + 1 + nb2 + 1] * d2_scaling
          }
        }
      }else if(!a1_fix & !b1_fix & !a2_fix & !b2_fix){

        if(c1_fix & c2_fix){

          C1 <- init_c1
          C2 <- init_c2

          if(d1_fix & d2_fix){

            A1 <- exp(theta[2]) * a1_scaling
            B1 <- theta[3] * b1_scaling
            D1 <- init_d1

            A2 <- theta[4] * a2_scaling
            B2 <- theta[5] * b2_scaling
            D2 <- init_d2

          }else if(!d1_fix & !d2_fix){

            A1 <- exp(theta[2]) * a1_scaling
            B1 <- theta[3] * b1_scaling
            D1 <- theta[4] * d1_scaling

            A2 <- theta[5] * a2_scaling
            B2 <- theta[6] * b2_scaling
            D2 <- theta[7] * d2_scaling

          }

        }else if(!c1_fix & !c2_fix){

          if(d1_fix & d2_fix){

            A1 <- exp(theta[2]) * a1_scaling
            B1 <- theta[3] * b1_scaling
            C1 <- theta[3 + 1:nb1] * c1_scaling
            D1 <- init_d1

            A2 <- theta[3 + nb1 + 1] * a2_scaling
            B2 <- theta[3 + nb1 + 2] * b2_scaling
            C2 <- theta[3 + nb1 + 2 + 1:nb2] * c2_scaling
            D2 <- init_d2

          }else if(!d1_fix & !d2_fix){

            A1 <- exp(theta[2]) * a1_scaling
            B1 <- theta[3] * b1_scaling
            C1 <- theta[3 + 1:nb1] * c1_scaling
            D1 <- theta[3 + nb1 + 1] * d1_scaling

            A2 <- theta[3 + nb1 + 2] * a2_scaling
            B2 <- theta[3 + nb1 + 3] * b2_scaling
            C2 <- theta[3 + nb1 + 3 + 1:nb2] * c2_scaling
            D2 <- theta[3 + nb1 + 3 + nb2 + 1] * d1_scaling
          }
        }
      }
    }else if(scale_horizontal_fix & !scale_vertical_fix){

      SCALE_HORIZONTAL <- init_scale_horizontal
      SCALE_VERTICAL <- exp(theta[1]) * vertical_scale_scaling

      if(a1_fix & b1_fix & a2_fix & b2_fix){

        A1 <- init_a1
        B1 <- init_b1
        A2 <- init_a2
        B2 <- init_b2

        if(c1_fix & c2_fix){

          C1 <- init_c1
          C2 <- init_c2

          if(d1_fix & d2_fix){

            D1 <- init_d1
            D2 <- init_d2

          }else if(!d1_fix & !d2_fix){

            D1 <- theta[2] * d1_scaling
            D2 <- theta[3] * d2_scaling

          }

        }else if(!c1_fix & !c2_fix){

          if(d1_fix & d2_fix){

            C1 <- theta[1 + 1:nb1] * c1_scaling
            D1 <- init_d1

            C2 <- theta[1 + nb1 + 1:nb2] * c2_scaling
            D2 <- init_d2

          }else if(!d1_fix & !d2_fix){

            C1 <- theta[1 + 1:nb1] * c1_scaling
            D1 <- theta[1 + nb1 + 1] * d1_scaling

            C2 <- theta[1 + nb1 + 1 + 1:nb2] * c2_scaling
            D2 <- theta[1 + nb1 + 1 + nb2 + 1] * d2_scaling

          }

        }

      }else if(!a1_fix & !b1_fix & !a2_fix & !b2_fix){

        if(c1_fix & c2_fix){

          C1 <- init_c1
          C2 <- init_c2

          if(d1_fix & d2_fix){

            A1 <- exp(theta[2]) * a1_scaling
            B1 <- theta[3] * b1_scaling
            D1 <- init_d1

            A2 <- theta[4] * a2_scaling
            B2 <- theta[5] * b2_scaling
            D2 <- init_d2

          }else if(!d1_fix & !d2_fix){

            A1 <- exp(theta[2]) * a1_scaling
            B1 <- theta[3] * b1_scaling
            D1 <- theta[4] * d1_scaling

            A2 <- theta[5] * a2_scaling
            B2 <- theta[6] * b2_scaling
            D2 <- theta[7] * d2_scaling
          }
        }else if(!c1_fix & !c2_fix){

          if(d1_fix & d2_fix){

            A1 <- exp(theta[2]) * a1_scaling
            B1 <- theta[3] * b1_scaling
            C1 <- theta[3 + 1:nb1] * c1_scaling
            D1 <- init_d1

            A2 <- theta[3 + nb1 + 1] * a2_scaling
            B2 <- theta[3 + nb1 + 2] * b2_scaling
            C2 <- theta[3 + nb1 + 2 + 1:nb2] * c2_scaling
            D2 <- init_d2

          }else if(!d1_fix & !d2_fix){

            A1 <- exp(theta[2]) * a1_scaling
            B1 <- theta[3] * b1_scaling
            C1 <- theta[3 + 1:nb1] * c1_scaling
            D1 <- theta[3 + nb1 + 1] * d1_scaling

            A2 <- theta[3 + nb1 + 2] * a2_scaling
            B2 <- theta[3 + nb1 + 3] * b2_scaling
            C2 <- theta[3 + nb1 + 3 + 1:nb2] * c2_scaling
            D2 <- theta[3 + nb1 + 3 + nb2 + 1] * d1_scaling
          }
        }
      }
    }
  }

  param <-  list(BETA = BETA, SCALE_HORIZONTAL = SCALE_HORIZONTAL,
                 SCALE_VERTICAL = SCALE_VERTICAL, A1 = A1, B1 = B1, C1 = C1,
                 D1 = D1, A2 = A2, B2 = B2, C2 = C2, D2 = D2)

  return(param)

}
