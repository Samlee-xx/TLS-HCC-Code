
##### The relationship of T/B TLS with TLS maturation states ######

data = read.csv('...\\TLS APS.csv')
colnames(data) = c('TLS.subtype','Maturation.state')

t.tls = subset(data, TLS.subtype == 'B-TLS')
df_summary <- t.tls %>%
  group_by(Maturation.state) %>%
  summarise(Count = n()) %>%
  mutate(Percent = Count / sum(Count) * 100,
         Label = paste0(Maturation.state, " (", round(Percent, 1), "%)"))

state_colors <- c(
  "A" = "#F58220",  
  "S" = "#009639",  
  "P" = "#F7B589"   
)

a = ggplot(df_summary, aes(x = "", y = Percent, fill = Maturation.state)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = state_colors) +
  theme_void() +
  geom_text(aes(label = Label),
            position = position_stack(vjust = 0.5),
            size = 4.5, color = "black") +
  labs(title = "Maturation State Distribution",
       fill = "Maturation State") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

pdf(file = '...\\b-tls aps.pdf',
    height = 6,width = 6)
print(a)
dev.off()


