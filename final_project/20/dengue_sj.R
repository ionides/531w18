library('pomp')
library('tidyverse')


sj_dengue = read_csv('http://dengueforecasting.noaa.gov/Training/San_Juan_Training_Data.csv')
sj_dengue = sj_dengue %>% mutate(total_other = other_positive_cases + additional_cases, week = 0:(nrow(sj_dengue)-1))
for(i in 7:nrow(sj_dengue)){
  unassigned = sj_dengue$total_other[i]
  sum_prev5_days = sj_dengue[c(
    'denv1_cases', 'denv2_cases','denv3_cases','denv4_cases')][(i-5):(i-1),] %>%
    apply(2, sum)
  assignment = rmultinom(1, unassigned, sum_prev5_days)
  sj_dengue[i,]$denv1_cases = sj_dengue[i,]$denv1_cases + assignment['denv1_cases', 1]
  sj_dengue[i,]$denv2_cases = sj_dengue[i,]$denv2_cases + assignment['denv2_cases', 1]
  sj_dengue[i,]$denv3_cases = sj_dengue[i,]$denv3_cases + assignment['denv3_cases', 1]
  sj_dengue[i,]$denv4_cases = sj_dengue[i,]$denv4_cases + assignment['denv4_cases', 1]
}
dd = sj_dengue %>% select(week, y1 = denv1_cases, y2 = denv2_cases, y3 = denv3_cases, y4 = denv4_cases)

# states: S, Y_1, ... Y_4, R_1, ... , R_4, Y_12, ... Y_43, R (total of 22)
# params: mu, sigma, g_ij, phi_ij (1 + 1 + 12 + 12 = total of 26, where mu and sigma known)

dengue_skel <- Csnippet('
                        int numstrains = 4;
                         double phi[4][4] = {{0, phi12, phi13, phi14}, {phi21, 0, phi23, phi24}, {phi31, phi32, 0, phi34}, {phi41, phi42, phi43, 0}};
                         double n_sec[4][4] = {{0, n12, n13, n14}, {n21, 0, n23, n24}, {n31, n32, 0, n34}, {n41, n42, n43, 0}};
                         double gammas[4][4] = {{0, gamma12, gamma13, gamma14}, {gamma21, 0, gamma23, gamma24}, {gamma31, gamma32, 0, gamma34}, {gamma41, gamma42, gamma43, 0}};
                         double lambdas[4] = {0, 0, 0, 0};
                         double n[4] = {n1, n2, n3, n4};
                         double rs[4] = {r1, r2, r3, r4};
                         double deltas[10];
                         double deltas_sec[4][4];
                         double sum_lambdas = 0;
                         double sumprod = 0;
                         double sumprodtwo = 0;
                         for(int i = 0; i < numstrains; i++){
                            sumprod = 0;
                            for(int j = 0; j < numstrains; j++){
                               if(i == j){
                                  sumprod += 0;
                               }
                               else{
                                  sumprod += phi[j][i]*n_sec[j][i];
                               }
                            }
                            lambdas[i] = bet*(n[i] + sumprod);
                            sum_lambdas += lambdas[i];
                         }
                         for(int i = 0; i < numstrains; i++){
                            sumprodtwo = 0;
                            for(int j = 0; j < numstrains; j++){
                               if(i != j){
                                  sumprodtwo += gammas[i][j]*lambdas[j];
                               }
                            }
                            deltas[i+5] = sigma*n[i] - (rs[i]*(mu + sumprodtwo));
                         }
                         deltas[0] = mu - s*sum_lambdas - mu*s;
                         for(int i = 0; i<numstrains; i++){
                            deltas[i+1] = s*lambdas[i] - (sigma + mu)*n[i];
                         }
                         double sum_ys = 0;
                         for(int i = 0; i<numstrains; i++){
                            for(int j = 0; j<numstrains; j++){
                               if(i == j) {
                                  deltas_sec[i][j] = 0;
                               }
                               else{
                                  deltas_sec[i][j] = rs[i]*gammas[i][j]*lambdas[j] - (sigma + mu)*n_sec[i][j];
                                  sum_ys += n_sec[i][j];
                               }
                            }
                         }
                         deltas[9] = sigma*sum_ys - mu*r;
                         Ds = (s + deltas[0]);
                         Dn1 = (n1 + deltas[1]);
                         Dn2 = (n2 + deltas[2]);
                         Dn3 = (n3 + deltas[3]);
                         Dn4 = (n4 + deltas[4]);
                         Dr1 = (r1 + deltas[5]);
                         Dr2 = (r2 + deltas[6]);
                         Dr3 = (r3 + deltas[7]);
                         Dr4 = (r4 + deltas[8]);
                         Dr = (r + deltas[9]);
                         Dn12 = (n12 + deltas_sec[0][1]);
                         Dn13 = (n13 + deltas_sec[0][2]);
                         Dn14 = (n14 + deltas_sec[0][3]);
                         Dn21 = (n21 + deltas_sec[1][0]);
                         Dn23 = (n23 + deltas_sec[1][2]);
                         Dn24 = (n24 + deltas_sec[1][3]);
                         Dn31 = (n31 + deltas_sec[2][0]);
                         Dn32 = (n32 + deltas_sec[2][1]);
                         Dn34 = (n34 + deltas_sec[2][3]);
                         Dn41 = (n41 + deltas_sec[3][0]);
                         Dn42 = (n42 + deltas_sec[3][1]);
                         Dn43 = (n43 + deltas_sec[3][2]);
                         ')

dengue_rprocess <- Csnippet('
                        int numstrains = 4;
                        int population = 2217968;
                        double phi[4][4] = {{0, phi12, phi13, phi14}, {phi21, 0, phi23, phi24}, {phi31, phi32, 0, phi34}, {phi41, phi42, phi43, 0}};
                        double n_sec[4][4] = {{0, n12, n13, n14}, {n21, 0, n23, n24}, {n31, n32, 0, n34}, {n41, n42, n43, 0}};
                        double gammas[4][4] = {{0, gamma12, gamma13, gamma14}, {gamma21, 0, gamma23, gamma24}, {gamma31, gamma32, 0, gamma34}, {gamma41, gamma42, gamma43, 0}};
                        double lambdas[4] = {0, 0, 0, 0};
                        double n[4] = {n1, n2, n3, n4};
                        double rs[4] = {r1, r2, r3, r4};
                        double deltas[10];
                        double deltas_sec[4][4];
                        double sum_lambdas = 0;
                        double sumprod = 0;
                        double sumprodtwo = 0;
                        double sumprodtwos[4];
                        double s_rate[5], s_trans[5], p_inf_rate[8], p_inf_trans[8], r_rate[16], r_trans[16], s_inf_rate[24], s_inf_trans[24], rec_rate[1], rec_trans[1];
                        int sum_deaths = 0;

                        // add up natural deaths from each box
                        double deaths[22];
                        deaths[0] = rbinom(s, 1-exp(-mu));
                        deaths[1] = rbinom(n1, 1-exp(-mu));
                        deaths[2] = rbinom(n2, 1-exp(-mu));
                        deaths[3] = rbinom(n3, 1-exp(-mu));
                        deaths[4] = rbinom(n4, 1-exp(-mu));
                        deaths[5] = rbinom(r1, 1-exp(-mu));
                        deaths[6] = rbinom(r2, 1-exp(-mu));
                        deaths[7] = rbinom(r3, 1-exp(-mu));
                        deaths[8] = rbinom(r4, 1-exp(-mu));
                        deaths[9] = rbinom(n12, 1-exp(-mu));
                        deaths[10] = rbinom(n13, 1-exp(-mu));
                        deaths[11] = rbinom(n14, 1-exp(-mu));
                        deaths[12] = rbinom(n21, 1-exp(-mu));
                        deaths[13] = rbinom(n23, 1-exp(-mu));
                        deaths[14] = rbinom(n24, 1-exp(-mu));
                        deaths[15] = rbinom(n31, 1-exp(-mu));
                        deaths[16] = rbinom(n32, 1-exp(-mu));
                        deaths[17] = rbinom(n34, 1-exp(-mu));
                        deaths[18] = rbinom(n41, 1-exp(-mu));
                        deaths[19] = rbinom(n42, 1-exp(-mu));
                        deaths[20] = rbinom(n43, 1-exp(-mu));
                        deaths[21] = rbinom(r, 1-exp(-mu));
                        for(int i = 0; i < 22; i++){
                           sum_deaths += deaths[i];
                        }

                        // total births
                        int births = sum_deaths;

                        // calculate foi for each strain
                        for(int i = 0; i < numstrains; i++){
                          sumprod = 0;
                          for(int j = 0; j < numstrains; j++){
                            if(i == j){
                              sumprod += 0;
                            }
                            else{
                              sumprod += phi[j][i]*n_sec[j][i];
                            }
                          }
                          lambdas[i] = bet*(n[i] + sumprod);
                          sum_lambdas += lambdas[i];
                        }
                        
                        // rates from S to primary I
                        s_rate[0] = s*lambdas[0]; // to I1
                        s_rate[1] = s*lambdas[1]; // to I2
                        s_rate[2] = s*lambdas[2]; // to I3
                        s_rate[3] = s*lambdas[3]; // to I4
                        s_rate[4] = s*mu; // natural death

                        // transitions from S to primary I
                        reulermultinom(5, s, &s_rate[0], 1, &s_trans[0]);

                        for(int i = 0; i < numstrains; i++){
                          sumprodtwo = 0;
                          for(int j = 0; j < numstrains; j++){
                            if(i != j){
                              sumprodtwo += gammas[i][j]*lambdas[j];
                            }
                          }
                          deltas[i+5] = sigma*n[i] - (rs[i]*(mu + sumprodtwo));
                          sumprodtwos[i] = sumprodtwo;
                        }
                        // rates from I_i to R_i and death
                        p_inf_rate[0] = sigma*n[0];
                        p_inf_rate[1] = mu*n[0];
                        p_inf_rate[2] = sigma*n[1];
                        p_inf_rate[3] = mu*n[1];
                        p_inf_rate[4] = sigma*n[2];
                        p_inf_rate[5] = mu*n[2];
                        p_inf_rate[6] = sigma*n[3];
                        p_inf_rate[7] = mu*n[3];

                        // transitions from I_i to R_i and death
                        reulermultinom(2, n[0], &p_inf_rate[0], 1, &p_inf_trans[0]);
                        reulermultinom(2, n[1], &p_inf_rate[2], 1, &p_inf_trans[2]);
                        reulermultinom(2, n[2], &p_inf_rate[4], 1, &p_inf_trans[4]);
                        reulermultinom(2, n[3], &p_inf_rate[6], 1, &p_inf_trans[6]);

                        deltas[0] = mu - s*sum_lambdas - mu*s;
                        for(int i = 0; i<numstrains; i++){
                          deltas[i+1] = s*lambdas[i] - (sigma + mu)*n[i];
                        }
                        double sum_ys = 0;
                        for(int i = 0; i<numstrains; i++){
                          for(int j = 0; j<numstrains; j++){
                            if(i == j) {
                              deltas_sec[i][j] = 0;
                            }
                            else{
                              deltas_sec[i][j] = rs[i]*gammas[i][j]*lambdas[j] - (sigma + mu)*n_sec[i][j];
                              sum_ys += n_sec[i][j];
                            }
                          }
                        }

                        // rates from R_i to I_ij and death
                        r_rate[0] = gamma12*lambdas[1]*r[0];
                        r_rate[1] = gamma13*lambdas[2]*r[0];
                        r_rate[2] = gamma14*lambdas[3]*r[0];
                        r_rate[3] = mu*r[0];
                        r_rate[4] = gamma21*lambdas[0]*r[1];
                        r_rate[5] = gamma23*lambdas[2]*r[1];
                        r_rate[6] = gamma24*lambdas[3]*r[1];
                        r_rate[7] = mu*r[1];
                        r_rate[8] = gamma31*lambdas[0]*r[2];
                        r_rate[9] = gamma32*lambdas[1]*r[2];
                        r_rate[10] = gamma34*lambdas[3]*r[2];
                        r_rate[11] = mu*r[2];
                        r_rate[12] = gamma41*lambdas[0]*r[3];
                        r_rate[13] = gamma42*lambdas[1]*r[3];
                        r_rate[14] = gamma43*lambdas[2]*r[3];
                        r_rate[15] = mu*r[3];
                        

                        // transitions from R_i to I_ij and death
                        reulermultinom(4, r1, &r_rate[0], 1, &r_trans[0]);
                        reulermultinom(4, r2, &r_rate[4], 1, &r_trans[4]);
                        reulermultinom(4, r3, &r_rate[8], 1, &r_trans[8]);
                        reulermultinom(4, r4, &r_rate[12], 1, &r_trans[12]);
                        
                        deltas[9] = sigma*sum_ys - mu*r;

                        // rates from I_ij to R
                        s_inf_rate[0] = sigma*n12;
                        s_inf_rate[1] = mu*n12;
                        s_inf_rate[2] = sigma*n13;
                        s_inf_rate[3] = mu*n13;
                        s_inf_rate[4] = sigma*n14;
                        s_inf_rate[5] = mu*n14;
                        s_inf_rate[6] = sigma*n21;
                        s_inf_rate[7] = mu*n21;
                        s_inf_rate[8] = sigma*n23;
                        s_inf_rate[9] = mu*n23;
                        s_inf_rate[10] = sigma*n24;
                        s_inf_rate[11] = mu*n24;
                        s_inf_rate[12] = sigma*n31;
                        s_inf_rate[13] = mu*n31;
                        s_inf_rate[14] = sigma*n32;
                        s_inf_rate[15] = mu*n32;
                        s_inf_rate[16] = sigma*n34;
                        s_inf_rate[17] = mu*n34;
                        s_inf_rate[18] = sigma*n41;
                        s_inf_rate[19] = mu*n41;
                        s_inf_rate[20] = sigma*n42;
                        s_inf_rate[21] = mu*n42;
                        s_inf_rate[22] = sigma*n43;
                        s_inf_rate[23] = mu*n43;

                        // transitions from secondary infection to permanent recovery
                        reulermultinom(2, n12, &s_inf_rate[0], 1, &s_inf_trans[0]);
                        reulermultinom(2, n13, &s_inf_rate[2], 1, &s_inf_trans[2]);                        
                        reulermultinom(2, n14, &s_inf_rate[4], 1, &s_inf_trans[4]);
                        reulermultinom(2, n21, &s_inf_rate[6], 1, &s_inf_trans[6]);
                        reulermultinom(2, n23, &s_inf_rate[8], 1, &s_inf_trans[8]);
                        reulermultinom(2, n24, &s_inf_rate[10], 1, &s_inf_trans[10]);
                        reulermultinom(2, n31 &s_inf_rate[12], 1, &s_inf_trans[12]);
                        reulermultinom(2, n32 &s_inf_rate[14], 1, &s_inf_trans[14]);
                        reulermultinom(2, n34, &s_inf_rate[16], 1, &s_inf_trans[16]);
                        reulermultinom(2, n41, &s_inf_rate[18], 1, &s_inf_trans[18]);
                        reulermultinom(2, n42, &s_inf_rate[20], 1, &s_inf_trans[20]);
                        reulermultinom(2, n43, &s_inf_rate[22], 1, &s_inf_trans[22]);

                        // rate from permanent recovery to death
                        rec_rate[0] = mu*r;

                        // transitions from permanent recovery to death
                        reulermultinom(1, r, &rec_rate[0], 1, &rec_trans[0]);
                        
                        // updates
                        s += births - s_trans[0] - s_trans[1] - s_trans[2] - s_trans[3] - s_trans[4];
                        n1 += s_trans[0] - p_inf_trans[0] - p_inf_trans[1];
                        n2 += s_trans[1] - p_inf_trans[2] - p_inf_trans[3];
                        n3 += s_trans[2] - p_inf_trans[4] - p_inf_trans[5];
                        n4 += s_trans[3] - p_inf_trans[6] - p_inf_trans[7];
                        r1 += p_inf_trans[0] - r_trans[0] - r_trans[1] - r_trans[2] - r_trans[3];
                        r2 += p_inf_trans[2] - r_trans[4] - r_trans[5] - r_trans[6] - r_trans[7];
                        r3 += p_inf_trans[4] - r_trans[8] - r_trans[9] - r_trans[10] - r_trans[11];
                        r4 += p_inf_trans[6] - r_trans[12] - r_trans[13] - r_trans[14] - r_trans[15];
                        n12 += r_trans[0] - s_inf_trans[0] - s_inf_trans[1];
                        n13 += r_trans[1] - s_inf_trans[2] - s_inf_trans[3];
                        n14 += r_trans[2] - s_inf_trans[4] - s_inf_trans[5];
                        n21 += r_trans[4] - s_inf_trans[6] - s_inf_trans[7];
                        n23 += r_trans[5] - s_inf_trans[8] - s_inf_trans[9];
                        n24 += r_trans[6] - s_inf_trans[10] - s_inf_trans[11];
                        n31 += r_trans[8] - s_inf_trans[12] - s_inf_trans[13];
                        n32 += r_trans[9] - s_inf_trans[14] - s_inf_trans[15];
                        n34 += r_trans[10] - s_inf_trans[16] - s_inf_trans[17];
                        n41 += r_trans[12] - s_inf_trans[18] - s_inf_trans[19];
                        n42 += r_trans[13] - s_inf_trans[20] - s_inf_trans[21];
                        n43 += r_trans[14] - s_inf_trans[22] - s_inf_trans[23];
                        r += s_inf_trans[0] + s_inf_trans[2] + s_inf_trans[4] + s_inf_trans[6] + s_inf_trans[8] + s_inf_trans[10] + s_inf_trans[12] + s_inf_trans[14] + s_inf_trans[16] + s_inf_trans[18] + s_inf_trans[20] + s_inf_trans[22];
                        Ds = (s + deltas[0]);
                        Dn1 = (n1 + deltas[1]);
                        Dn2 = (n2 + deltas[2]);
                        Dn3 = (n3 + deltas[3]);
                        Dn4 = (n4 + deltas[4]);
                        Dr1 = (r1 + deltas[5]);
                        Dr2 = (r2 + deltas[6]);
                        Dr3 = (r3 + deltas[7]);
                        Dr4 = (r4 + deltas[8]);
                        Dr = (r + deltas[9]);
                        Dn12 = (n12 + deltas_sec[0][1]);
                        Dn13 = (n13 + deltas_sec[0][2]);
                        Dn14 = (n14 + deltas_sec[0][3]);
                        Dn21 = (n21 + deltas_sec[1][0]);
                        Dn23 = (n23 + deltas_sec[1][2]);
                        Dn24 = (n24 + deltas_sec[1][3]);
                        Dn31 = (n31 + deltas_sec[2][0]);
                        Dn32 = (n32 + deltas_sec[2][1]);
                        Dn34 = (n34 + deltas_sec[2][3]);
                        Dn41 = (n41 + deltas_sec[3][0]);
                        Dn42 = (n42 + deltas_sec[3][1]);
                        Dn43 = (n43 + deltas_sec[3][2]);
                        ')
df = data.frame(week = c(0:1000000), y1 = c(0:1000000), y2 = c(0:1000000), y3 = c(0:1000000), y4 = c(0:1000000))
dengue_pomp <- pomp(data = df,
                times = "week",
                skeleton = map(dengue_skel, 1/7),
                paramnames = c("mu","sigma", "bet", 
                               "phi12", "phi13", "phi14",
                               "phi21", "phi23", "phi24",
                               "phi31", "phi32", "phi34",
                               "phi41", "phi42", "phi43",
                               "gamma12", "gamma13", "gamma14",
                               "gamma21", "gamma23", "gamma24",
                               "gamma31", "gamma32", "gamma34",
                               "gamma41", "gamma42", "gamma43",
                               "s.0",
                               "n1.0", "n2.0", "n3.0", "n4.0",
                               "r1.0", "r2.0", "r3.0", "r4.0",
                               "n12.0", "n13.0", "n14.0",
                               "n21.0", "n23.0", "n24.0",
                               "n31.0", "n32.0", "n34.0",
                               "n41.0", "n42.0", "n43.0",
                               "r.0"),
                statenames = c("s",
                               "n1", "n2", "n3", "n4",
                               "r1", "r2", "r3", "r4",
                               "n12", "n13", "n14",
                               "n21", "n23", "n24",
                               "n31", "n32", "n34",
                               "n41", "n42", "n43",
                               "r"),
                t0 = 0
)
  traj = trajectory(dengue_pomp, params = c(mu = 1/(70*365), sigma = 1/3.65, bet = 400/365, 
                                        phi12 = 1.9, phi13 = 1.9, phi14 = 1.9,
                                        phi21 = 1.9, phi23 = 1.9, phi24 = 1.9,
                                        phi31 = 1.9, phi32 = 1.9, phi34 = 1.9,
                                        phi41 = 1.9, phi42 = 1.9, phi43 = 1.9,
                                        gamma12 = 1.8, gamma13 = 1.8, gamma14 = 1.8,
                                        gamma21 = 1.8, gamma23 = 1.8, gamma24 = 1.8,
                                        gamma31 = 1.8, gamma32 = 1.8, gamma34 = 1.8,
                                        gamma41 = 1.8, gamma42 = 1.8, gamma43 = 1.8,
                                        s.0 = (0.1),
                                        n1.0 = (0.004), n2.0 = (0.002), n3.0 = (0.001), n4.0 = (0.003),
                                        r1.0 = 0.196, r2.0 = 0.198, r3.0 = 0.199, r4.0 = 0.197,
                                        n12.0 = 0.01, n13.0 = 0.01, n14.0 = 0.01,
                                        n21.0 = 0.01, n23.0 = 0.01, n24.0 = 0.01,
                                        n31.0 = 0.01, n32.0 = 0.01, n34.0 = 0.01,
                                        n41.0 = 0.01, n42.0 = 0.01, n43.0 = 0.01,
                                        r.0 = 0))
  
  traj2 = as.data.frame(t(as.data.frame(traj[,1,]))) %>% mutate(week = df$week)
  ggplot(data=(traj2%>% filter(week > 52000 & week < 57200))) + 
    geom_line(mapping = aes(x = week/52, y = n1 + n12 + n13 + n14, color = "Strain 1"), alpha = 0.5) + 
    geom_line(mapping = aes(x = week/52, y = n2 + n21 + n23 + n24, color = "Strain 2"), alpha = 0.5) + 
    geom_line(mapping = aes(x = week/52, y = n3 + n31 + n32 + n34, color = "Strain 3"), alpha = 0.5) + 
    geom_line(mapping = aes(x = week/52, y = n4 + n41 + n42 + n43, color = "Strain 4"), alpha = 0.5) +
    labs(x = "Year", y = "Proportion Infected") + 
    ggtitle("phi = 1.9, gamma = 1.8")
    scale_colour_manual("", 
                        breaks = c("Strain 1", "Strain 2", "Strain 3", "Strain 4"),
                        values = c("red", "blue", "green", "yellow"))
  
dengue_skel <- Csnippet('
                      int numstrains = 4;
                      double phi[4][4] = {{0, phi12, phi13, phi14}, {phi21, 0, phi23, phi24}, {phi31, phi32, 0, phi34}, {phi41, phi42, phi43, 0}};
                      double n_sec[4][4] = {{0, n12, n13, n14}, {n21, 0, n23, n24}, {n31, n32, 0, n34}, {n41, n42, n43, 0}};
                      double gammas[4][4] = {{0, gamma12, gamma13, gamma14}, {gamma21, 0, gamma23, gamma24}, {gamma31, gamma32, 0, gamma34}, {gamma41, gamma42, gamma43, 0}};
                      double lambdas[4] = {0, 0, 0, 0};
                      double n[4] = {n1, n2, n3, n4};
                      double rs[4] = {r1, r2, r3, r4};
                      double deltas[10];
                      double deltas_sec[4][4];
                      double sum_lambdas = 0;
                      double sumprod = 0;
                      double sumprodtwo = 0;
                      for(int i = 0; i < numstrains; i++){
                        sumprod = 0;
                        for(int j = 0; j < numstrains; j++){
                          if(i == j){
                            sumprod += 0;
                          }
                          else{
                            sumprod += phi[j][i]*n_sec[j][i];
                          }
                        }
                        lambdas[i] = bet*(n[i] + sumprod);
                        sum_lambdas += lambdas[i];
                      }
                      for(int i = 0; i < numstrains; i++){
                        sumprodtwo = 0;
                        for(int j = 0; j < numstrains; j++){
                          if(i != j){
                            sumprodtwo += gammas[i][j]*lambdas[j];
                          }
                        }
                        deltas[i+5] = sigma*n[i] - (rs[i]*(mu + sumprodtwo));
                      }
                      deltas[0] = mu - s*sum_lambdas - mu*s;
                      for(int i = 0; i<numstrains; i++){
                        deltas[i+1] = s*lambdas[i] - (sigma + mu)*n[i];
                      }
                      double sum_ys = 0;
                      for(int i = 0; i<numstrains; i++){
                        for(int j = 0; j<numstrains; j++){
                          if(i == j) {
                            deltas_sec[i][j] = 0;
                          }
                          else{
                            deltas_sec[i][j] = rs[i]*gammas[i][j]*lambdas[j] - (sigma + mu)*n_sec[i][j];
                            sum_ys += n_sec[i][j];
                          }
                        }
                      }
                      deltas[9] = sigma*sum_ys - mu*r;
                      Ds = (s + deltas[0]);
                      Dn1 = (n1 + deltas[1]);
                      Dn2 = (n2 + deltas[2]);
                      Dn3 = (n3 + deltas[3]);
                      Dn4 = (n4 + deltas[4]);
                      Dr1 = (r1 + deltas[5]);
                      Dr2 = (r2 + deltas[6]);
                      Dr3 = (r3 + deltas[7]);
                      Dr4 = (r4 + deltas[8]);
                      Dr = (r + deltas[9]);
                      Dn12 = (n12 + deltas_sec[0][1]);
                      Dn13 = (n13 + deltas_sec[0][2]);
                      Dn14 = (n14 + deltas_sec[0][3]);
                      Dn21 = (n21 + deltas_sec[1][0]);
                      Dn23 = (n23 + deltas_sec[1][2]);
                      Dn24 = (n24 + deltas_sec[1][3]);
                      Dn31 = (n31 + deltas_sec[2][0]);
                      Dn32 = (n32 + deltas_sec[2][1]);
                      Dn34 = (n34 + deltas_sec[2][3]);
                      Dn41 = (n41 + deltas_sec[3][0]);
                      Dn42 = (n42 + deltas_sec[3][1]);
                      Dn43 = (n43 + deltas_sec[3][2]);
                        ')
df = data.frame(week = c(0:1000000), y1 = c(0:1000000), y2 = c(0:1000000), y3 = c(0:1000000), y4 = c(0:1000000))
dengue_pomp <- pomp(data = df,
                    times = "week",
                    skeleton = map(dengue_skel, 1/7),
                    paramnames = c("mu","sigma", "bet", 
                                   "phi12", "phi13", "phi14",
                                   "phi21", "phi23", "phi24",
                                   "phi31", "phi32", "phi34",
                                   "phi41", "phi42", "phi43",
                                   "gamma12", "gamma13", "gamma14",
                                   "gamma21", "gamma23", "gamma24",
                                   "gamma31", "gamma32", "gamma34",
                                   "gamma41", "gamma42", "gamma43",
                                   "s.0",
                                   "n1.0", "n2.0", "n3.0", "n4.0",
                                   "r1.0", "r2.0", "r3.0", "r4.0",
                                   "n12.0", "n13.0", "n14.0",
                                   "n21.0", "n23.0", "n24.0",
                                   "n31.0", "n32.0", "n34.0",
                                   "n41.0", "n42.0", "n43.0",
                                   "r.0"),
                    statenames = c("s",
                                   "n1", "n2", "n3", "n4",
                                   "r1", "r2", "r3", "r4",
                                   "n12", "n13", "n14",
                                   "n21", "n23", "n24",
                                   "n31", "n32", "n34",
                                   "n41", "n42", "n43",
                                   "r"),
                    t0 = 0
)
                        
