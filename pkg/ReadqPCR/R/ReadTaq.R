ReadTaq <- function(filename = character()) {
# If no filename given look for .taq files in current directory
if (length(filename) == 0) {
    taqfile.path <- getwd()
    files = list.files(taqfile.path)
    filename = files[grep("\\.[tT][aA][qQ]$",files)]
    if (length(filename)>1) {
      stop("Too many Taq files in current directory","\n")
    }
}
if (length(filename) == 0) {
    stop("No Taqman filenname specified and no Taqman files in current directory:","\n")
}
if (length(filename)>1) {
    stop ("Too may Taqman files specified")
}
taq2eset(filename)
}

