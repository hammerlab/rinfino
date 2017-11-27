
#' Helper to execute shell commands (linux)
#' @param cmdargs (str) arguments to the command
#' @param command (str) name of command to be executed (full path will be looked up using Sys.which)
#' @param internal (bool) passed to \code{system} - whether result should be returned as R object.
#' @import glue
shell <- function(cmdargs, command="bash", intern = F, ...) {
	command_path <- Sys.which(command)
	system(glue::glue('{command_path} {cmdargs}'), intern = intern, ...)
}
