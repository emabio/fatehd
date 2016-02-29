get.num.param <- function(pl, flag, split = " "){
  as.numeric(unlist(strsplit(grep(flag, pl, value = TRUE), split = split))[-1])
}

get.char.param <- function(pl, flag, split.flag = "^--.*--$"){
  flags.pos <- grep(split.flag, pl)
  flag.pos <- which(flags.pos == grep(flag, pl))
  as.character(pl[(flags.pos[flag.pos]+1):(flags.pos[flag.pos+1]-1)])
}