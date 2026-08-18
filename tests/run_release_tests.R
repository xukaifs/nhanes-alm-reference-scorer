options(stringsAsFactors = FALSE)
if (!file.exists(file.path("R", "score_conditional_alm.R"))) stop("Run from repository root.")
source(file.path("R", "score_conditional_alm.R"), encoding="UTF-8")
bundles <- load_alm_reference()
meta <- validate_reference_bundles(bundles)
stopifnot(nrow(meta)==10, setequal(meta$imputation,1:5))
stopifnot(sum(meta$sex=="Female")==5, sum(meta$sex=="Male")==5)
stopifnot(all(meta$family[meta$sex=="Female"]=="BCCG"))
stopifnot(all(meta$family[meta$sex=="Male"]=="BCT"))

fit_f <- read.csv(file.path("models","audit","corrected_H_full_fit_Female.csv"))
fit_m <- read.csv(file.path("models","audit","corrected_H_full_fit_Male.csv"))
stopifnot(nrow(fit_f)==5, nrow(fit_m)==5, all(fit_f$converged), all(fit_m$converged))
stopifnot(all(fit_f$n==7890), all(fit_m$n==8128))

example <- data.frame(
  sex=c("Female","Male","Female","Male"),
  age=c(18,69,45,60),
  height_m=c(1.60,1.75,1.60,1.75),
  alm_kg=c(15,22,15,22)
)
res <- score_conditional_alm(example,bundles=bundles,warn=FALSE)$pooled
stopifnot(nrow(res)==4, all(res$score_valid), all(is.finite(res$model_averaged_z)))

cent <- reference_centiles(example[,c("sex","age","height_m")],bundles=bundles,warn=FALSE)
stopifnot(nrow(cent)==4, all(cent$centile_valid), all(is.finite(cent$alm_p50)))

bad <- score_conditional_alm(data.frame(
  sex=c("Female","Male"), age=c(17,70), height_m=c(1.60,1.75), alm_kg=c(15,22)
),bundles=bundles,warn=FALSE)$pooled
stopifnot(all(!bad$score_valid), all(is.na(bad$model_averaged_z)))

message("v2.0.0 structural and live-scoring tests passed.")
