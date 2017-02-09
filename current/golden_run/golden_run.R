library(ggplot2)

df <- read.csv("physics_code/current/outputFiles/v2_all_GR.csv", header = TRUE)
df['ratio']<-df$num_of_events/df$total_q
df<-subset(df, ratio<=500000)

plot<-ggplot(data=df, aes(x=run_num, y=ratio)) +
  geom_point(aes(color=factor(file_num))) +
  xlab("Run Num") + ylab("Ratio") +
  ggtitle("Golden Run") +
  scale_fill_gradient(low="#003b6f",high="#b20000") + theme(axis.text.x = element_text(angle = 90, vjust = -0.1))
print(plot)
