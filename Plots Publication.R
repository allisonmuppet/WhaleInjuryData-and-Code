#### Publication plots for variable importance
# Single entanglement and vessel strike plot respectively for both species combined
# Requires combining two rfPermute objects

rm(list=ls())

output.files = dir(pattern="RData")

for (i in 1:length(output.files))  {

output <- load(output.files[i])

df1 <- as.data.frame(importance(model.entangle))
df1$Variable = row.names(importance(model.entangle))

df2 <- as.data.frame(importance(model.vessel))
df2$Variable = row.names(importance(model.vessel))

df1 <- cbind.data.frame(df1$Variable, df1$MeanDecreaseAccuracy)
df2 <- cbind.data.frame(df2$Variable, df2$MeanDecreaseAccuracy)

names(df1) <- c("Variable", "MeanDecreaseAccuracy")
names(df2) <- c("Variable", "MeanDecreaseAccuracy")

df1$Variable <- factor(df1$Variable,                                    # Factor levels in decreasing order
          levels = df1$Variable[order(df1$MeanDecreaseAccuracy, decreasing = FALSE)])

df2$Variable <- factor(df2$Variable,                                    # Factor levels in decreasing order
                       levels = df2$Variable[order(df2$MeanDecreaseAccuracy, decreasing = FALSE)])


jpeg(file=paste(species, "Variable Importance.jpg"), width=2400, height=1200, units="px", res=300)
p1 <- ggplot(data=df1, aes(x=Variable, y=MeanDecreaseAccuracy)) + geom_col() + ggtitle(paste("(A)", species, "Entanglement")) + coord_flip() +  theme(plot.title = element_text(size = 10))
p2 <- ggplot(data=df2, aes(x=Variable, y=MeanDecreaseAccuracy)) + geom_col() + ggtitle(paste("(B)", species, "Vessel")) + coord_flip() +  theme(plot.title = element_text(size = 10))

p1 <- p1 + ylab("Mean Decrease Accuracy")
p2 <- p2 + ylab("Mean Decrease Accuracy")

grid.arrange(p1, p2, nrow=1)

dev.off()

}