set_names <- function(variable_name){
  tools::toTitleCase(gsub("_", " ", variable_name))
}