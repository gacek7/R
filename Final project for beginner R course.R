library("ggplot2") #rysowanie wykres�w
library("scales") #potrzebne do wyswietlenia procentow na pie chart
library("reshape2") #potrzebne do meltowania danych do heatmapy, aby zobaczy� korelacje
library("plotly") #do heatmapy
library("heatmaply") #potrzebne do heatmapy
summary(data)
#sprawdzenie i przygotowanie danych
which(data$blueAvgLevel > 18)
which(data$redAvgLevel > 18)
dim(data) #teraz sprawadzmy, ile mamy obserwacji oraz zmiennych
colSums(is.na(data)) #tutaj sprawdzamy czy s� jakie� braki w danych
#tutaj opisa� metode analizy, czemu tak a nie inaczeji
blueWin <- data[data$blueWins==1,]
redWin <- data[data$redWins==1,]
blueLose <- data[data$blueWins==0,]
redLose <- data[data$redWins==0,]
#zmniejszamy liczbe zmiennych po to, �eby w danym dataframe by�y informacje tylko o danej dru�ynie
blueWin <- blueWin[,2:20]
redWin <- redWin[,21:39]
blueLose <- blueLose[,2:20]
redLose <- redLose[,21:39]
#zmniejszamy skale dla obra�en oraz zdobytego z�ota, 1 = 1 tysi�c, b�dzie �adniej na wykresach :)
blueWin$blueChampionDamageDealt <- blueWin$blueChampionDamageDealt/1000
blueLose$blueChampionDamageDealt <- blueLose$blueChampionDamageDealt/1000
redWin$redChampionDamageDealt <- redWin$redChampionDamageDealt/1000
redLose$redChampionDamageDealt <- redLose$redChampionDamageDealt/1000
#tutaj te� zmniejszamy skale ale d�a z�ota, 1 = 1 tysi�c
blueWin$blueTotalGold <- blueWin$blueTotalGold/1000
blueLose$blueTotalGold <- blueLose$blueTotalGold/1000
redWin$redTotalGold <- redWin$redTotalGold/1000
redLose$redTotalGold <- redLose$redTotalGold/1000
#zabite miniony b�d� w dziesi�tkach
blueWin$blueTotalMinionKills <- blueWin$blueTotalMinionKills/10
redWin$redTotalMinionKills <- redWin$redTotalMinionKills/10
blueLose$blueTotalMinionKills <- blueLose$blueTotalMinionKills/10
redLose$redTotalMinionKills <- redLose$redTotalMinionKills/10
#tutaj �rednie wszystkich zmiennych dla ka�dego wariantu i dru�yny
BWmean <- apply(blueWin,2,mean)
BLmean <- apply(blueLose,2,mean)
RWmean <- apply(redWin,2,mean)
RLmean <- apply(redLose,2,mean)
#zaokr�glamy do drugiego miejsca po przecinku
BWmean <- round(BWmean,2)
BLmean <- round(BLmean,2)
RWmean <- round(RWmean,2)
RLmean <- round(RLmean,2)
#dla wygody zmienimy nazwy ze �rednimi wynikami na 0, dopiszemy je potem aby by�o czytelniej
names(BWmean) <- NULL
names(BLmean) <- NULL
names(RWmean) <- NULL
names(RLmean) <- NULL
#Raczej koniec przygotowywania danych tutaj, jak co� dopisz� to przepraszam ale z czasem coraz wi�cej pomys��w :)
#Sprawdzmy ile przeci�tnie trwa rozgrywka
mean(data$gameDuration)/60
gameduration <- data.frame(gameDuration = data$gameDuration/60)
hist(gameduration$gameDuration,main="Czas rozgrywki", freq = FALSE, breaks = 40, col="grey",xlab="Czas w minutach")
curve(dnorm(x, mean=avggameDuration, sd=sd(gameduration$gameDuration)), add=TRUE, lwd = 4, lty="solid", col = "black")
nrow(data[data$gameDuration<240,])/nrow(data) #procent sko�czonych gier przed 4 minut�
(nrow(data[data$gameDuration<1000,])-nrow(data[data$gameDuration<900,]))/nrow(data) #procent gier sko�czonych przez early surr vote
#Wygrana zale�na od strony po kt�rej sie zacze�o 
pie_chart <- data.frame(names=c("Blue Wins","Red Wins"),compare = c(sum(blueWin$blueWins),sum(redWin$redWins)))
pie <- ggplot(pie_chart,aes(x="",y=compare,fill=names)) + geom_bar(width = 1, stat="identity")
plot <- pie + coord_polar("y",start=0)
plot + scale_fill_manual(values=c("#3333FF","#FF0000")) + geom_text(aes(y = compare/2 + c(0, cumsum(compare)[-length(compare)]),
                                                                        label = percent(1-compare/sum(pie_chart[2]))), size=10)
#nasze �rednie wyniki zamieniamy w dataframe i dopisujemy nazwy zmiennych
meanOverall <- data.frame(id=1:19,team=rep(c("Blue","Red"),each=19),Win=c(BWmean,RWmean),Lose=c(BLmean,RLmean),Variable=c("Win","First Blood",
                                                                                                                          "First Tower","First Baron","First Dragon","First inhibitor","Dragon Kills","Baron Kills","Tower Kills","Inhibitor Kills",
                                                                                                                          "WardPlaced","WardKills","Kills","Deaths","Asissts","Damage Dealt","Total Gold","Total minions","Average level"))
#�rednie wyniki gdy blue wygrywa a red przegrywa
meanBluevsRed <- data.frame(id=1:19,team=rep(c("Blue","Red"),each=19),Win=c(BWmean,RLmean),Lose=c(RWmean,BLmean),Variable=c("Win","First Blood",
                                                                                                                            "First Tower","First Baron","First Dragon","First inhibitor","Dragon Kills","Baron Kills","Tower Kills","Inhibitor Kills",
                                                                                                                            "WardPlaced","WardKills","Kills","Deaths","Asissts","Damage Dealt","Total Gold","Total minions","Average level"))
#�rednie wyniki gdy red wygrywa a blue przegrywa
meanRedvsBlue <- data.frame(id=1:19,team=rep(c("Red","Blue"),each=19),Win=c(RWmean,BLmean),Lose=c(BWmean,RLmean),Variable=c("Win","First Blood",
                                                                                                                            "First Tower","First Baron","First Dragon","First inhibitor","Dragon Kills","Baron Kills","Tower Kills","Inhibitor Kills",
                                                                                                                            "WardPlaced","WardKills","Kills","Deaths","Asissts","Damage Dealt","Total Gold","Total minions","Average level"))


#obie dru�yny zwyci�skie przedstawione na wykresie, por�wnanie statystyk
plotWinWin <- ggplot(meanOverall,aes(x = id,y = Win,fill = factor(team))) +
  geom_bar(stat = "identity",colour = NA,position = "dodge") +
  scale_fill_manual(values = c("#3333FF","#FF0000")) + 
  scale_x_reverse() + 
  geom_text(aes(label = Win),position = position_dodge(1),hjust = -0.25,size = 3,fontface = "bold") +
  geom_text(aes(x = id,y = -5.5,label = Variable),vjust = .5,fontface = "bold") + 
  labs(fill = "Team") +
  guides(fill = guide_legend(reverse = TRUE)) +
  coord_flip() + 
  theme_void() +
  theme(legend.position = c(.92,.75))
plotWinWin
#Niebieska wygrywa czerwona przegrywa
plotWinLose <- ggplot(meanBluevsRed,aes(x = id,y = -Win,fill = factor(team))) +
  geom_bar(stat = "identity",colour = NA,position = "dodge") +
  scale_fill_manual(values = c("#3333FF","#FF0000")) + 
  scale_x_reverse() + 
  geom_text(aes(label = Win),position = position_dodge(1),hjust = 1.2,size = 3,fontface = "bold") +
  labs(fill = "Team") +
  guides(fill = guide_legend(reverse = TRUE)) +
  coord_flip() + 
  theme_void() +
  theme(legend.position = c(.15,.75))
plotWinLose
#Niebieska przegrywa czerwona wygrywa
plotLoseWin <- ggplot(meanRedvsBlue,aes(x = id,y = Win,fill = factor(team))) +
  geom_bar(stat = "identity",colour = NA,position = "dodge") +
  scale_fill_manual(values = c("#3333FF","#FF0000")) + 
  scale_x_reverse() + 
  geom_text(aes(label = Win),position = position_dodge(1),hjust = -0.25,size = 3,fontface = "bold") +
  geom_text(aes(x = id,y = -5.5,label = Variable),vjust = .5,fontface = "bold") + 
  labs(fill = "Team") +
  guides(fill = guide_legend(reverse = TRUE)) +
  coord_flip() + 
  theme_void() +
  theme(legend.position = c(.92,.75))
plotLoseWin
#korelacje dla wygranej niebieskiej, funkcja melt zamienia kolumne w pojedy�cz� warto��
cormat <- round(cor(redWin[2:18]),2)
melted_cormat <- melt(cormat)
#tworze macierz korelacji
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
#edytuje macierz korelacji aby dosta� �adnie zedytowane wyniki
ggheatmap + geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
#wszystkie macierze korelacji wygl�daj� bardzo podobnie, bliskie 1 wyst�puj� w tych samych miejscach