library(benford.analysis)
library(data.table)

property_values <- fread("property values.txt", header = T)

summary(property_values)

# Subsetting only the numeric columns since we're only interested in them
library(purrr)
property_values$BBLE <- as.numeric(property_values$BBLE)
data <- keep(property_values, is.numeric)
data <- na.omit(data)
ids <- data$BBLE
data <- within(data, rm(BBLE, BLOCK, LOT, ZIP))
rownames(data) <- ids
mahalanobis_score <- abs(data.frame(ids, mahalanobis(data, colMeans(data), cov(data), tol=1e-20)))

# z-score

scaled_data <- scale(data)
z_scores <- rowSums(scaled_data)
m_scores <- mahalanobis_score[,2]

# benford
library(benford.analysis)
benford_scores_final <- NULL
for(i in 1:8497){
        benford_scores_final[i] <- benford_scores[[i]]$MAD
}

scores <- data.frame(cbind(ids, benford_scores_final, m_scores, z_scores))
scores <- scale(scores)
scores <- data.frame(scores)
scores_ensemble <- (0.2 * scores$benford_scores_final) + (0.4 * scores$m_scores) + 
        (0.4 * scores$z_scores)
scores_ensemble <- data.frame(cbind(ids, scores_ensemble))
scores_ensemble$scores_ensemble <- scores_ensemble$scores_ensemble + 1

colnames(scores_ensemble) <- c("BBLE", "score")

property_values <- fread("property values.txt", header = T)

final_dataset <- merge(scores_ensemble, property_values, by = "BBLE")

save(final_dataset, file = "final_dataset.RData")
write.csv(scores_ensemble, file = "scores_tableau.csv")

a <- data.frame(quantile(scores_ensemble$score, seq(0.05, 1.00, 0.05)))
write.csv(a, file = "quantiles.csv")

write.csv(final_dataset[1:15,], file = "final_top_15.csv")

final_dataset <- final_dataset[order(final_dataset$score, decreasing = T),]