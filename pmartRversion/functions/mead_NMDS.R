mead_NMDS <- function(XX,ZZ){
  
  library(ggplot2)
  
  # Extract component scores
  NMDS1 <- data.frame(scores(XX))$NMDS1
  NMDS2 <- data.frame(scores(XX))$NMDS2
  
  testZZ <- table(ZZ)
  if(any(testZZ < 3)){
    names <- names(which(testZZ < 3))
    ids <- which(ZZ %in% names)
    ZZ <- ZZ[-ids]
    ZZ <- droplevels(ZZ)
    NMDS1 <- NMDS1[-ids]
    NMDS2 <- NMDS2[-ids]
  }
  
  # Format treatment, "group"
  Treatment <- ZZ
  if(any(levels(Treatment) == "")){
    Treatment <- as.character(Treatment)
    Treatment <- as.factor(Treatment)
  }
  
  # Aggregate data using mean
  NMDS <- data.frame(NMDS1, NMDS2, Treatment)
  NMDS.mean = aggregate(NMDS[,1:2], list(group=Treatment), mean)
  
  # Function for drawing ellipses
  veganCovEllipse <- function(cov, center = c(0, 0), scale = 1, npoints = 100){
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
  }
  
  # Create a new dataframe with data and ellipses
  df_ell <- data.frame()
  for(g in levels(NMDS$Treatment)){
    df_ell <- rbind(df_ell, 
                    cbind(as.data.frame(with(NMDS[NMDS$Treatment==g,],
                    veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2))))),
                    group=g))
  }
  
  # Plot
  X1 <- ggplot(data = NMDS, aes(NMDS1, NMDS2)) + 
    geom_point(aes(color = Treatment), size=2.5, alpha=0.75) +
    geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2, colour=group), size=1.5, linetype=5)+
    theme_bw()+
    theme(aspect.ratio=1,
          axis.text.x=element_text(size=20),
          axis.text.y=element_text(size=20),
          axis.title.x=element_text(size=20),
          axis.title.y=element_text(size=20),
          legend.title=element_text(size=15),
          legend.text=element_text(size=15),
          panel.grid=element_blank())
  X1    
}

# mead_PCA <- function(XX, ZZ){
#   
#   plotdata <- data.frame(SampleID = rownames(scores(XX)), PC1 = data.frame(scores(XX))$NMDS1, PC2 = data.frame(scores(XX))$NMDS2, Group = ZZ)
#   
#   ggplot(plotdata, aes(x=PC1, y=PC2)) +
#     geom_point(aes(colour=Group), size=2.5, alpha=0.75)+
#     theme_bw()+
#     theme(aspect.ratio=1,
#           axis.text.x=element_text(size=20),
#           axis.text.y=element_text(size=20),
#           axis.title.x=element_text(size=20),
#           axis.title.y=element_text(size=20),
#           legend.title=element_text(size=15),
#           legend.text=element_text(size=15),
#           panel.grid=element_blank())
#   
# }

