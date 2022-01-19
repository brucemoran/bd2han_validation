#runder.R

run_DEP <- function(){
  source("R/DEP_BD2Han_LFQ.R")
}

render_rmd <- function(){
  lapply(dir("R", pattern = "Rmd"), function(f){
    print(paste0("Rendering: ", f))
    fname <- gsub(".Rmd", "", f)
    if(!exists(paste0("output/", fname, "/", fname, ".html"))){
      rmarkdown::render(input = paste0("R/", f), output_dir = paste0("output/",   fname))
    } else {
      print(paste0("HTML exists for: ", f))
    }
  })
}

run_DEP()
render_rmd()
