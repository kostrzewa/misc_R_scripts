source("~/code/R/misc_R_scripts/prop_error/binary_operations.R")

pi0_diffs <- function(mpi,mpi0c,mpi0f,mpi0d) {
  mpisq <- compute_product(mpi,mpi)
  mpi0csq <- compute_product(mpi0c,mpi0c)
  mpi0dsq <- compute_product(mpi0d,mpi0d)
  mpi0fsq <- compute_product(mpi0f,mpi0f)

  mpi0c_m_mpi <- compute_difference(mpi0c,mpi)
  mpi0csq_m_mpisq <- compute_difference(mpi0csq,mpisq)
  mpi_m_mpi0f <- compute_difference(mpi,mpi0f)
  mpisq_m_mpi0fsq <- compute_difference(mpisq,mpi0fsq)
  mpi_m_mpi0d <- compute_difference(mpi,mpi0d)
  mpisq_m_mpi0dsq <- compute_difference(mpisq,mpi0dsq)
  mpi0c_m_mpi0f <- compute_difference(mpi0c,mpi0f)
  mpi0csq_m_mpi0fsq <- compute_difference(mpi0csq,mpi0fsq)
  mpi0c_m_mpi0d <- compute_difference(mpi0c,mpi0d)
  mpi0csq_m_mpi0dsq <- compute_difference(mpi0csq,mpi0dsq)

  mpi0c_m_mpi_rel <- compute_ratio(mpi0c_m_mpi,mpi)
  mpi_m_mpi0f_rel <- compute_ratio(mpi_m_mpi0f,mpi)
  mpi_m_mpi0d_rel <- compute_ratio(mpi_m_mpi0d,mpi)
  mpi0c_m_mpi0f_rel <- compute_ratio(mpi0c_m_mpi0f,mpi)
  mpi0c_m_mpi0d_rel <- compute_ratio(mpi0c_m_mpi0d,mpi)
  
  data.frame(mpi=mpi$val, dmpi=mpi$dval,
             mpi0c=mpi0c$val, dmpi0c=mpi0c$dval,
             mpi0f=mpi0f$val, dmpi0f=mpi0f$dval,
             mpi0d=mpi0d$val, dmpi0d=mpi0d$dval,
             mpi0c_m_mpi=mpi0c_m_mpi$val, dmpi0c_m_mpi=mpi0c_m_mpi$dval,
             mpi0csq_m_mpisq=mpi0csq_m_mpisq$val, dmpi0csq_m_mpisq=mpi0csq_m_mpisq$dval,
             mpi_m_mpi0f=mpi_m_mpi0f$val, dmpi_m_mpi0f=mpi_m_mpi0f$dval,
             mpisq_m_mpi0fsq=mpisq_m_mpi0fsq$val, dmpisq_m_mpi0fsq=mpisq_m_mpi0fsq$dval,
             mpi_m_mpi0d=mpi_m_mpi0d$val, dmpi_m_mpi0d=mpi_m_mpi0d$dval,
             mpisq_m_mpi0dsq=mpisq_m_mpi0dsq$val, dmpisq_m_mpi0dsq=mpisq_m_mpi0dsq$dval,
             mpi0c_m_mpi0f=mpi0c_m_mpi0f$val, dmpi0c_m_mpi0f=mpi0c_m_mpi0f$dval,
             mpi0csq_m_mpi0fsq=mpi0csq_m_mpi0fsq$val, dmpi0csq_m_mpi0fsq=mpi0csq_m_mpi0fsq$dval,
             mpi0c_m_mpi0d=mpi0c_m_mpi0d$val, dmpi0c_m_mpi0d=mpi0c_m_mpi0d$dval,
             mpi0csq_m_mpi0dsq=mpi0csq_m_mpi0dsq$val, dmpi0csq_m_mpi0dsq=mpi0csq_m_mpi0dsq$dval,
             mpi0c_m_mpi_rel=mpi0c_m_mpi_rel$val, dmpi0c_m_mpi_rel=mpi0c_m_mpi_rel$dval,
             mpi_m_mpi0f_rel=mpi_m_mpi0f_rel$val, dmpi_m_mpi0f_rel=mpi_m_mpi0f_rel$dval,
             mpi_m_mpi0d_rel=mpi_m_mpi0d_rel$val, dmpi_m_mpi0d_rel=mpi_m_mpi0d_rel$dval,
             mpi0c_m_mpi0f_rel=mpi0c_m_mpi0f_rel$val, dmpi0c_m_mpi0f_rel=mpi0c_m_mpi0f_rel$dval,
             mpi0c_m_mpi0d_rel=mpi0c_m_mpi0d_rel$val, dmpi0c_m_mpi0d_rel=mpi0c_m_mpi0d_rel$dval)
}
