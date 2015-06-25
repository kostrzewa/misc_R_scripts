source("~/code/R/misc_R_scripts/prop_error/binary_operations.R")

pi0_diffs <- function(mpi,mpi0c,mpi0f,mpi0d) {
  mpi0c_m_mpi <- compute_difference(mpi0c,mpi)
  mpi_m_mpi0f <- compute_difference(mpi,mpi0f)
  mpi_m_mpi0d <- compute_difference(mpi,mpi0d)
  mpi0c_m_mpi0f <- compute_difference(mpi0c,mpi0f)
  mpi0c_m_mpi0d <- compute_difference(mpi0c,mpi0d)

  mpi0c_m_mpi_rel <- compute_ratio(mpi0c_m_mpi,mpi)
  mpi_m_mpi0f_rel <- compute_ratio(mpi_m_mpi0f,mpi)
  mpi_m_mpi0d_rel <- compute_ratio(mpi_m_mpi0d,mpi)
  mpi0c_m_mpi0f_rel <- compute_ratio(mpi0c_m_mpi0f,mpi)
  mpi0c_m_mpi0d_rel <- compute_ratio(mpi0c_m_mpi0d,mpi)
  
  data.frame(mpi0c_m_mpi=mpi0c_m_mpi$val, dmpi0c_m_mpi=mpi0c_m_mpi$dval, 
             mpi_m_mpi0f=mpi_m_mpi0f$val, dmpi_m_mpi0f=mpi_m_mpi0f$dval,
             mpi_m_mpi0d=mpi_m_mpi0d$val, dmpi_m_mpi0d=mpi_m_mpi0d$dval,
             mpi0c_m_mpi0f=mpi0c_m_mpi0f$val, dmpi0c_m_mpi0f=mpi0c_m_mpi0f$dval,
             mpi0c_m_mpi0d=mpi0c_m_mpi0d$val, dmpi0c_m_mpi0d=mpi0c_m_mpi0d$dval,
             mpi0c_m_mpi_rel=mpi0c_m_mpi_rel$val, dmpi0c_m_mpi_rel=mpi0c_m_mpi_rel$dval,
             mpi_m_mpi0f_rel=mpi_m_mpi0f_rel$val, dmpi_m_mpi0f_rel=mpi_m_mpi0f_rel$dval,
             mpi_m_mpi0d_rel=mpi_m_mpi0d_rel$val, dmpi_m_mpi0d_rel=mpi_m_mpi0d_rel$dval,
             mpi0c_m_mpi0f_rel=mpi0c_m_mpi0f_rel$val, dmpi0c_m_mpi0f_rel=mpi0c_m_mpi0f_rel$dval,
             mpi0c_m_mpi0d_rel=mpi0c_m_mpi0d_rel$val, dmpi0c_m_mpi0d_rel=mpi0c_m_mpi0d_rel$dval)
}
