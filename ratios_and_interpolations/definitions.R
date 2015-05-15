# these are "definitions" of single quantities and ratios which can be used in the driver script

define.meson.quants <- function(datanames,light_masses,strange_masses,charm_masses) {
  quants <- NULL
  quants[["m_pi"]] <- list(name="m_pi", texlabel="$aM_\\pi$", m1=light_masses, m2=light_masses, datanames=datanames$ll_c, type="mps" )
  quants[["m_pi0c"]] <- list(name="m_pi0c", texlabel="$aM_{\\pi_0^{(c)}}$", m1=light_masses, m2=light_masses, datanames=datanames$ll_nud, type="mps" ) 
  quants[["m_K"]] <- list(name="m_K", texlabel="$aM_K$", m1=light_masses, m2=strange_masses, datanames=datanames$ls_c, type="mps" )
  quants[["m_D"]] <- list(name="m_D", texlabel="$aM_D$", m1=light_masses, m2=charm_masses, datanames=datanames$lc_c, type="mps" )
  quants[["m_Ds"]] <- list(name="m_Ds", texlabel="$aM_{D_s}$", m1=strange_masses, m2=charm_masses, datanames=datanames$sc_c, type="mps" )
  quants[["f_pi"]] <- list(name="f_pi", texlabel="$af_\\pi$", m1=light_masses, m2=light_masses, datanames=datanames$ll_c, type="fps" )
  quants[["f_K"]] <- list(name="f_K", texlabel="$af_K$", m1=light_masses, m2=strange_masses, datanames=datanames$ls_c, type="fps" )
  quants[["f_D"]] <- list(name="f_D", texlabel="$af_D$", m1=light_masses, m2=charm_masses, datanames=datanames$lc_c, type="fps" )
  quants[["f_Ds"]] <- list(name="f_Ds", texlabel="$af_{D_s}$", m1=strange_masses, m2=charm_masses, datanames=datanames$sc_c, type="fps" )
  quants
}

define.meson.ratios <- function(quants) {
  ratios <- NULL
  ratios[["m_pi_ov_f_pi"]] <- list( name="m_pi_ov_f_pi", texlabel="$M_\\pi/f_\\pi$", dividend=quants[["m_pi"]], divisor=quants[["f_pi"]])
  ratios[["m_K_ov_f_K"]] <- list( name="m_K_ov_f_K", texlabel="$M_K/f_K$", dividend=quants[["m_K"]], divisor=quants[["f_K"]])
  ratios[["m_D_ov_f_D"]] <- list( name="m_D_ov_f_D", texlabel="$M_D/f_D$", dividend=quants[["m_D"]], divisor=quants[["f_D"]])
  ratios[["m_Ds_ov_f_Ds"]] <- list( name="m_Ds_ov_f_Ds", texlabel="$M_{D_s}/f_{D_s}$", dividend=quants[["m_Ds"]], divisor=quants[["f_Ds"]])

  ratios[["m_K_ov_f_pi"]] <- list( name="m_K_ov_f_pi", texlabel="$M_K/f_\\pi$", dividend=quants[["m_K"]], divisor=quants[["f_pi"]])
  ratios[["m_D_ov_f_pi"]] <- list( name="m_D_ov_f_pi", texlabel="$M_D/f_\\pi$", dividend=quants[["m_D"]], divisor=quants[["f_pi"]])
  ratios[["m_Ds_ov_f_pi"]] <- list( name="m_Ds_ov_f_pi", texlabel="$M_{D_s}/f_\\pi$", dividend=quants[["m_Ds"]], divisor=quants[["f_pi"]])

  ratios[["m_Ds_ov_m_D"]] <- list( name="m_Ds_ov_m_D", texlabel="$M_{D_s}/M_D$", dividend=quants[["m_Ds"]], divisor=quants[["m_D"]])
  ratios[["m_Ds_ov_m_K"]] <- list( name="m_Ds_ov_m_K", texlabel="$M_{D_s}/M_K$", dividend=quants[["m_Ds"]], divisor=quants[["m_K"]])
  ratios[["m_Ds_ov_m_pi"]] <- list( name="m_Ds_ov_m_pi", texlabel="$M_{D_s}/M_\\pi$", dividend=quants[["m_Ds"]], divisor=quants[["m_pi"]])

  ratios[["m_D_ov_m_K"]] <- list( name="m_D_ov_m_K", texlabel="$M_{D}/M_K$", dividend=quants[["m_D"]], divisor=quants[["m_K"]])
  ratios[["m_D_ov_m_pi"]] <- list( name="m_D_ov_m_pi", texlabel="$M_{D}/M_\\pi$", dividend=quants[["m_D"]], divisor=quants[["m_pi"]])

  ratios[["m_K_ov_m_pi"]] <- list( name="m_K_ov_m_pi", texlabel="$M_K/M_\\pi$", dividend=quants[["m_K"]], divisor=quants[["m_pi"]])

  ratios[["f_K_ov_f_pi"]] <- list( name="f_K_ov_f_pi", texlabel="$f_K/f_\\pi$", dividend=quants[["f_K"]], divisor=quants[["f_pi"]])
  ratios[["f_D_ov_f_pi"]] <- list( name="f_D_ov_f_pi", texlabel="$f_D/f_\\pi$", dividend=quants[["f_D"]], divisor=quants[["f_pi"]])
  ratios[["f_Ds_ov_f_pi"]] <- list( name="f_Ds_ov_f_pi", texlabel="$f_{D_s}/f_\\pi$", dividend=quants[["f_Ds"]], divisor=quants[["f_pi"]])

  ratios[["f_D_ov_f_K"]] <- list( name="f_D_ov_f_K", texlabel="$f_D/f_K$", dividend=quants[["f_D"]], divisor=quants[["f_K"]])
  ratios[["f_Ds_ov_f_K"]] <- list( name="f_Ds_ov_f_K", texlabel="$f_{D_s}/f_K$", dividend=quants[["f_Ds"]], divisor=quants[["f_K"]])

  ratios[["f_Ds_ov_f_D"]] <- list( name="f_Ds_ov_f_D", texlabel="$f_{D_s}/f_D$", dividend=quants[["f_Ds"]], divisor=quants[["f_D"]])
  ratios
}
