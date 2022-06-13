Waiting <- read.table("Waiting.txt", header = T)
attach(Waiting)
an <- aov(waiting ~ course+ city + course:city)
summary(an)
#verification of assumptions:
shapiro.test(an$residuals)
x11()
qqnorm(an$residuals)
qqline(an$residuals)
#we can assume gaussianity IF we can also assure homoschedasticity amongst groups
x11()
boxplot(waiting ~ course+ city)
var.test(waiting[which(city == "Iasi" & course == "Starter")], waiting[which(city == "Iasi" & course == "Main")])
var.test(waiting[which(city == "Iasi" & course == "Starter")], waiting[which(city == "Iasi" & course == "Dessert")])
var.test(waiting[which(city == "Iasi" & course == "Main")], waiting[which(city == "Iasi" & course == "Dessert")])
var.test(waiting[which(city == "Iasi" & course == "Starter")], waiting[which(city == "Bucarest" & course == "Starter")])
var.test(waiting[which(city == "Iasi" & course == "Starter")], waiting[which(city == "Bucarest" & course == "Main")])
var.test(waiting[which(city == "Iasi" & course == "Starter")], waiting[which(city == "Bucarest" & course == "Dessert")])
var.test(waiting[which(city == "Iasi" & course == "Main")], waiting[which(city == "Bucarest" & course == "Starter")])
var.test(waiting[which(city == "Iasi" & course == "Main")], waiting[which(city == "Bucarest" & course == "Main")])
var.test(waiting[which(city == "Iasi" & course == "Main")], waiting[which(city == "Bucarest" & course == "Dessert")])
var.test(waiting[which(city == "Iasi" & course == "Dessert")], waiting[which(city == "Bucarest" & course == "Starter")])
var.test(waiting[which(city == "Iasi" & course == "Dessert")], waiting[which(city == "Bucarest" & course == "Main")])
var.test(waiting[which(city == "Iasi" & course == "Dessert")], waiting[which(city == "Bucarest" & course == "Dessert")])
var.test(waiting[which(city == "Bucarest" & course == "Starter")], waiting[which(city == "Bucarest" & course == "Dessert")])
var.test(waiting[which(city == "Bucarest" & course == "Starter")], waiting[which(city == "Bucarest" & course == "Main")])
var.test(waiting[which(city == "Bucarest" & course == "Dessert")], waiting[which(city == "Bucarest" & course == "Main")])
#we have homoschedasticity

#or more simply:
fact <- with(Waiting, interaction(course, city))
bartlett.test(waiting,fact)

#So the assumptions are satisfied

summary(an)

#we eliminate the non significant factors -> one at the time!
#backward selection
an_new <- aov(waiting ~ course + course:city)
summary(an_new)

an_new <- aov(waiting ~ course)
summary(an_new)

#so we end up with a reduced model
anova(an, an_new)#still significant
shapiro.test(an_new$residuals)
bartlett.test(an_new$residuals ~ course)
#still valid

s<- sum((an_new$residuals)^2)/an_new$df

alphaB <- 0.05/6
smean <- c(courseMain = an_new$coefficients[1] + an_new$coefficients[2],
          courseStarter = an_new$coefficients[1] + an_new$coefficients[3],
          courseDessert = an_new$coefficients[1])
#CI Main - Starter
n1 <- length(which(course == "Main"))
n2 <- length(which(course == "Starter"))
n3 <- length(which(course == "Dessert"))
n1
n2
n3
n <- 60

#Since:
#smean_a - smean_b ~ N(mu_a - mu_b, simga^2*(1/n_a + 1/n_b)) = N(mu_a - mu_b, 2*simga^2/n)
#quindi, per l'intervallo di confidenza:
#(smean_a - smean_b)/sqrt(2*sigma^2/n) ~ N(0,1)
# s = SS_res/(n-g) ~ sigma^2*Chi_square(n-g)

#(smean_a - smean_b)/sqrt(2*s/n) ~ t(n-g)

band <- qt(1-alphaB/2, an_new$df)*sqrt(s*2/n)
sdeltamean <- c(smean[1] - smean[2], smean[2] - smean[3], smean[3] - smean[1])
CI1 <- cbind(inf = sdeltamean - band,
         center = sdeltamean,
         sup = sdeltamean + band)
rownames(CI1)<- c("Main - Starter", "Starter - Dessert", "Dessert - Mean")
CI1
SS <- c(sum((an_new$residuals[which(course == "Main")])^2),sum((an_new$residuals[which(course == "Starter")])^2), sum((an_new$residuals[which(course == "Dessert")])^2))
chi_up <- qchisq(1-alphaB/2, n)
chi_down <- qchisq(alphaB/2, n)
CI2 <- cbind(inf = SS/chi_up,
             center = SS/(n-1),
             sup = SS/chi_down)
rownames(CI2)<-c("Main", "Starter", "Dessert")
CI2
var(waiting[which(course == "Main")])
