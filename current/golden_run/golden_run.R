library(ggplot2)

df <- read.csv("gr.csv", header = TRUE)
df["ratio"] <- df$num_of_events/df$total_q
df <- subset(df, ratio <= 5e+05)

plot <- ggplot(data = df, aes(x = run_num, y = ratio)) +
  geom_point(aes(color=file_num)) +
  xlab("Run Num") +
  ylab("Ratio") +
  ggtitle("Golden Run") +
  scale_fill_gradient(low = "#003b6f", high = "#b20000") +
  theme(axis.text.x = element_text(angle = 90, vjust = -0.1))
print(plot)


hist <- ggplot(data = df) +
  geom_density(aes(x = ratio, color = ratio), kernel = "gaussian", adjust = 0.5) +
  geom_histogram(aes(x = ratio, fill=file_num),bins = 200) +
  scale_fill_gradient(low = "#003b6f", high = "#b20000") +
  theme(axis.text.x = element_text(angle = 90, vjust = -0.1)) +
  xlab("Ratio") +
  ylab("Counts") +
  ggtitle("Golden Run")
print(hist)
