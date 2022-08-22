if(any(grepl("conda", .libPaths(), fixed = TRUE))){
  message("Setting libPaths")
  df = .libPaths()
  conda_i = which(grepl("conda", df, fixed = TRUE))
  .libPaths(c(df[conda_i], df[-conda_i]))
}

rmarkdown::render(
  input = snakemake@input[["markdown"]],
  clean = TRUE,
  intermediates_dir = snakemake@params[["output_dir"]],
  output_file = snakemake@output[["outfile"]],
  output_dir = snakemake@params[["output_dir"]],
  output_format = "all",
  params = list(
    rwd = snakemake@params[["rwd"]],
    out = snakemake@params[["output_name"]])
)

