### On the head node
#!/usr/bin/env Rscript
args = commandArgs(T)
rh = "/cluster/home/yjliu_jh/sbin/R/library/4.1.1"
tmp = paste(rh,"/tmp",sep="")
if(!file.exists(tmp)) {
  dir.create(tmp)
}
source("http://bioconductor.org/biocLite.R")
repos = c(biocinstallRepos(),getOption("repos"))


args <- "tanaylab/metacell"
getBiocRepos <- function() {
  
  BiocManager <- tryCatch(asNamespace("BiocManager"), error = identity)
  if (!inherits(BiocManager, "error"))
    return(BiocManager$repositories())
  
  BiocInstaller <- tryCatch(asNamespace("BiocInstaller"), error = identity)
  if (!inherits(BiocInstaller, "error"))
    return(BiocInstaller$biocinstallRepos())
  
  msg <- paste(
    "Neither BiocManager nor BiocInstaller are installed;",
    "cannot discover Bioconductor repositories"
  )
  warning(msg)
  character()
}

repos <- getBiocRepos()

pkg = utils:::getDependencies(args,,available.packages(contriburl=contrib.url(repos)))
order = download.packages(pkg, tmp, repos=repos)[,2]
write.table(order,paste(tmp,"/pkg_install_order.list",sep=""),col.name
            s=F,row.names=F,quote=F,sep="\t")

### On the compute node
#!/usr/bin/env Rscript
rh = R.home()
tmp = paste(rh,"/tmp",sep="")
if(!file.exists(tmp)) {
  quit()
}
pkg = read.table(paste(tmp,"/pkg_install_order.list",sep=""),sep="\t")
pkg = as.vector(t(pkg))
install.packages(pkg,repo=NULL)
file.remove(pkg)
file.remove(paste(tmp,"/pkg_install_order.list",sep=""))